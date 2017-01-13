program test_qn_solver_2d_polar
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_qn_solver_2d_polar, only: &
    sll_t_qn_solver_2d_polar, &
    sll_s_qn_solver_2d_polar_init, &
    sll_s_qn_solver_2d_polar_solve

  use m_test_case_qn_2d_base, only: &
    c_test_case_qn_solver_2d_polar

  use m_test_case_qn_2d_dirichlet, only: &
    t_test_dirichlet_zero_error

  use m_test_case_qn_2d_neumann_mode0, only: &
    t_test_neumann_mode0_zero_error

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_dirichlet, &
    sll_p_neumann, &
    sll_p_neumann_mode_0

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type(t_test_dirichlet_zero_error)     :: test_case_dirichlet
  type(t_test_neumann_mode0_zero_error) :: test_case_neumann_mode0

  sll_int32  :: nr
  sll_int32  :: nth
  sll_real64 :: error_norm
  sll_real64 :: tol

  logical :: success

  success = .true.

  !-----------------------------------------------------------------------------
  ! TEST #1: Dirichlet, solver should be exact
  !-----------------------------------------------------------------------------
  nr  = 64
  nth = 32
  tol = 1.0e-11_f64

  ! Define test case
  test_case_dirichlet%rmin                = 1.0_f64
  test_case_dirichlet%rmax                = 10.0_f64
  test_case_dirichlet%adiabatic_electrons = .true.
  test_case_dirichlet%use_zonal_flow      = .false.
  test_case_dirichlet%epsilon_0           = 1.0_f64
  test_case_dirichlet%bc_rmin             = sll_p_dirichlet
  test_case_dirichlet%bc_rmax             = sll_p_dirichlet

  call run_test( test_case_dirichlet, nr, nth, error_norm )

  ! Write relative error norm (global) to standard output
  write(*,"(/a)") "--------------------------------------------------------"
  write(*,"(a)")  "Homogeneous Dirichlet boundary conditions"
  write(*,"(a)")  "phi(r,theta) = (rmax-r)(r-rmin)(a+cos(k(theta-theta_0)))"
  write(*,"(a)")  "--------------------------------------------------------"
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

  ! Define test case
  test_case_neumann_mode0%rmin                = 1.0_f64
  test_case_neumann_mode0%rmax                = 10.0_f64
  test_case_neumann_mode0%adiabatic_electrons = .true.
  test_case_neumann_mode0%use_zonal_flow      = .false.
  test_case_neumann_mode0%epsilon_0           = 1.0_f64
  test_case_neumann_mode0%bc_rmin             = sll_p_neumann_mode_0
  test_case_neumann_mode0%bc_rmax             = sll_p_dirichlet

  call run_test( test_case_neumann_mode0, nr, nth, error_norm )

  ! Write relative error norm (global) to standard output
  write(*,"(/a)") "-----------------------------------------------------------&
       &------------------------------------------------"
  write(*,"(a)")  "Mixed Homogeneous Dirichlet / Neumann mode 0 boundary conditions"
  write(*,"(a)")  "phi(r,th) = a(rmax-r)(r+rmax-2rmin)/(rmax-rmin)^2 &
       &+ 4(rmax-r)(r-rmin)/(rmax-rmin)^2*b*cos(k(theta-theta_0))"
  write(*,"(a)")  "-----------------------------------------------------------&
       &------------------------------------------------"
  write(*,"(a,e11.3)") "Relative L_inf norm of error = ", error_norm
  write(*,"(a,e11.3)") "Tolerance                    = ", tol
  if (error_norm > tol) then
    success = .false.
    write(*,"(a/)") "!!! FAILED !!!"
  end if


  ! TODO: run convergence analysis on a more difficult test-case

  ! Check if test passed
  if(success) then
     write(*,"(/a/)") "--- PASSED ---"
  endif

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
contains
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  subroutine run_test( test_case, nr, nth, error_norm )
    class(c_test_case_qn_solver_2d_polar), intent(in   ) :: test_case
    sll_int32                            , intent(in   ) :: nr
    sll_int32                            , intent(in   ) :: nth
    sll_real64                           , intent(  out) :: error_norm

    type(sll_t_qn_solver_2d_polar) :: solver
    sll_real64                     :: max_phi
    sll_real64                     :: max_err

    sll_real64 ::  r,  th
    sll_real64 :: dr, dth
    sll_int32  :: i, j

    sll_real64, allocatable :: rhs   (:,:)
    sll_real64, allocatable :: phi_ex(:,:)
    sll_real64, allocatable :: phi   (:,:)

    sll_real64, allocatable :: rho_m0(:)
    sll_real64, allocatable :: b_magn(:)
    sll_real64, allocatable :: lambda(:)

    ! Computational grid
    dr  = (test_case%rmax - test_case%rmin) / nr
    dth = 2.0_f64*sll_p_pi / nth

    ! Allocate 2D distributed arrays (rho, phi, phi_ex) with layout_a
    allocate( rhs   (nr+1,nth) )
    allocate( phi_ex(nr+1,nth) )
    allocate( phi   (nr+1,nth) )

    ! Load analytical solution and rhs
    do j = 1, nth
      th = (j-1)*dth
      do i = 1, nr+1
        r = test_case%rmin + (i-1)*dr
        phi_ex(i,j) = test_case%phi_ex( r, th )
        rhs   (i,j) = test_case%rhs   ( r, th )
      end do
    end do
    phi(:,:) = 0.0_f64

    ! Allocate and load 1D radial profiles (needed by solver)
    allocate( rho_m0(nr+1) )
    allocate( b_magn(nr+1) )
    allocate( lambda(nr+1) )
    !
    do i = 1, nr+1
      r = test_case%rmin + (i-1)*dr
      rho_m0(i) = test_case%rho_m0( r )
      b_magn(i) = test_case%b_magn( r )
      lambda(i) = test_case%lambda( r )
    end do

    ! Initialize solver
    call sll_s_qn_solver_2d_polar_init( solver, &
      rmin     = test_case%rmin, &
      rmax     = test_case%rmax, &
      nr       = nr, &
      ntheta   = nth, &
      rho_m0   = rho_m0, &
      b_magn   = b_magn, &
      lambda   = lambda, &
      use_zonal_flow = test_case%use_zonal_flow, &
      epsilon_0      = test_case%epsilon_0     , &
      bc_rmin        = test_case%bc_rmin , &
      bc_rmax        = test_case%bc_rmax )

    ! Compute numerical phi for a given rhs
    call sll_s_qn_solver_2d_polar_solve( solver, rhs, phi )

    ! Compute maximum norms of phi_ex and error
    max_phi = maxval(abs( phi_ex ))
    max_err = maxval(abs( phi_ex-phi ))

    ! Global relative error in maximum norm
    error_norm = max_err / max_phi

  end subroutine run_test


end program test_qn_solver_2d_polar
