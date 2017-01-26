!> @authors Yaman Güçlü, IPP Garching
!> @authors Edoardo Zoni, IPP Garching
!>
!> @brief
!> Program to run tests of module 'sll_m_qn_solver_2d_polar.F90'.
!>
!> @details
!> Default parameters of each test type can be overwritten (except for BCs) before
!> call to subroutine 'run_test'.
!> Tolerance is set according to the order of the radial profile of each test type
!> (expected zero or non-zero numerical error) as well as according to the mesh sizes.
!> Recall that k (Fourier mode, see default parameters) must be <= ntheta/2.
!>

program test_qn_solver_2d_polar
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_qn_solver_2d_polar, only: &
    sll_t_qn_solver_2d_polar, &
    sll_s_qn_solver_2d_polar_init, &
    sll_s_qn_solver_2d_polar_solve

  use m_test_qn_solver_2d_polar_base, only: &
    c_test_qn_solver_2d_polar_base

  use m_test_qn_solver_2d_polar_dirichlet, only: &
    t_test_qn_solver_2d_polar_dirichlet_quadratic

  use m_test_qn_solver_2d_polar_neumann_mode0, only: &
    t_test_qn_solver_2d_polar_neumann_mode0_quadratic

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type(t_test_qn_solver_2d_polar_dirichlet_quadratic)     :: test_case_dirichlet
  type(t_test_qn_solver_2d_polar_neumann_mode0_quadratic) :: test_case_neumann_mode0

  sll_int32  :: nr, nth
  sll_real64 :: error_norm, tol

  logical :: success
  success = .true.

  !-----------------------------------------------------------------------------
  ! TEST #1: Dirichlet, solver should be exact
  !-----------------------------------------------------------------------------
  nr  = 64
  nth = 32
  tol = 1.0e-11_f64

  call run_test( test_case_dirichlet, nr, nth, error_norm )

  ! Write relative error norm (global) to standard output
  write(*,"(/a)") "------------------------------------------------------------"
  write(*,"(a)")  "Homogeneous Dirichlet boundary conditions"
  write(*,"(a)")  "phi(r,theta) = (rmax-r)(r-rmin)(a + b*cos(k(theta-theta_0)))"
  write(*,"(a)")  "------------------------------------------------------------"
  write(*,"(a,e11.3)") "Relative L_inf norm of error = ", error_norm
  write(*,"(a,e11.3)") "Tolerance                    = ", tol
  if (error_norm > tol) then
    success = .false.
    write(*,"(a)") "!!! FAILED !!!"
  end if

  !-----------------------------------------------------------------------------
  ! TEST #2: Neumann mode 0, solver should be exact
  !-----------------------------------------------------------------------------
  nr  = 64
  nth = 32
  tol = 1.0e-11_f64

  call run_test( test_case_neumann_mode0, nr, nth, error_norm )

  ! Write relative error norm (global) to standard output
  write(*,"(/a)") "-----------------------------------------------------------&
       &--------------------"
  write(*,"(a)")  "Mixed Homogeneous Dirichlet / Neumann mode 0 boundary conditions"
  write(*,"(a)")  "phi(r,theta) = a(r-rmax)(r-2rmin+rmax) &
       &+ b(r-rmax)(r-rmin)cos(k(theta-theta_0))"
  write(*,"(a)")  "-----------------------------------------------------------&
       &--------------------"
  write(*,"(a,e11.3)") "Relative L_inf norm of error = ", error_norm
  write(*,"(a,e11.3)") "Tolerance                    = ", tol
  if (error_norm > tol) then
    success = .false.
    write(*,"(a/)") "!!! FAILED !!!"
  end if

  ! Check if test passed
  if(success) then
     write(*,"(/a/)") "PASSED"
  endif

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
contains
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  subroutine run_test( test_case, nr, nth, error_norm )
    class(c_test_qn_solver_2d_polar_base), intent(in   ) :: test_case
    sll_int32                            , intent(in   ) :: nr
    sll_int32                            , intent(in   ) :: nth
    sll_real64                           , intent(  out) :: error_norm

    type(sll_t_qn_solver_2d_polar) :: solver
    sll_real64                     :: max_phi
    sll_real64                     :: max_err

    sll_real64 ::  r,  th
    sll_real64 :: dr, dth
    sll_int32  :: i, j

    sll_real64 :: rlim(2)
    sll_int32  :: bcs (2)
    logical    :: adiabatic_electrons
    logical    :: use_zonal_flow
    sll_real64 :: epsilon_0

    sll_real64, allocatable :: rho   (:,:)
    sll_real64, allocatable :: phi_ex(:,:)
    sll_real64, allocatable :: phi   (:,:)

    sll_real64, allocatable :: rho_m0(:)
    sll_real64, allocatable :: b_magn(:)
    sll_real64, allocatable :: lambda(:)

    ! Extract domain limits and boundary conditions
    rlim(:) = test_case%get_rlim()
    bcs (:) = test_case%get_bcs ()

    ! Extract test-case parameters
    call test_case%get_parameters( adiabatic_electrons, use_zonal_flow, epsilon_0 )

    ! Computational grid
    dr  = (rlim(2) - rlim(1)) / nr
    dth = 2.0_f64*sll_p_pi / nth

    ! Allocate 2D distributed arrays (rho, phi, phi_ex) with layout_a
    allocate( rho   (nr+1,nth) )
    allocate( phi_ex(nr+1,nth) )
    allocate( phi   (nr+1,nth) )

    ! Load analytical solution and rho
    do j = 1, nth
      th = (j-1)*dth
      do i = 1, nr+1
        r = rlim(1) + (i-1)*dr
        phi_ex(i,j) = test_case%phi_ex( r, th )
        rho   (i,j) = test_case%rho   ( r, th )
      end do
    end do
    phi(:,:) = 0.0_f64

    ! Allocate and load 1D radial profiles (needed by solver)
    allocate( rho_m0(nr+1) )
    allocate( b_magn(nr+1) )
    allocate( lambda(nr+1) )
    !
    do i = 1, nr+1
      r = rlim(1) + (i-1)*dr
      rho_m0(i) = test_case%rho_m0( r )
      b_magn(i) = test_case%b_magn( r )
      lambda(i) = test_case%lambda( r )
    end do

    ! Initialize solver
    call sll_s_qn_solver_2d_polar_init( solver, &
      rmin           = rlim(1), &
      rmax           = rlim(2), &
      nr             = nr, &
      ntheta         = nth, &
      rho_m0         = rho_m0, &
      b_magn         = b_magn, &
      lambda         = lambda, &
      use_zonal_flow = use_zonal_flow, &
      epsilon_0      = epsilon_0     , &
      bc_rmin        = bcs(1), &
      bc_rmax        = bcs(2) )

    ! Compute numerical phi for a given rho
    call sll_s_qn_solver_2d_polar_solve( solver, rho, phi )

    ! Compute maximum norms of phi_ex and error
    max_phi = maxval(abs( phi_ex ))
    max_err = maxval(abs( phi_ex-phi ))

    ! Global relative error in maximum norm
    error_norm = max_err / max_phi

  end subroutine run_test


end program test_qn_solver_2d_polar
