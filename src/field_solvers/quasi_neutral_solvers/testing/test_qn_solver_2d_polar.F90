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

program test_qn_solver_2d_polar
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"
#include "sll_errors.h"

  use sll_m_constants, only: &
    sll_p_twopi

  use sll_m_utilities, only: &
    sll_s_new_array_linspace

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

  class(c_test_qn_solver_2d_polar_base)                   , pointer :: test_case
  type (t_test_qn_solver_2d_polar_dirichlet_quadratic)    , target  :: test_case_dirichlet
  type (t_test_qn_solver_2d_polar_neumann_mode0_quadratic), target  :: test_case_neumann_mode0

  sll_int32  :: nr, nth
  sll_real64 :: error_norm, tol
  character(len=8) :: rgrid_opt

  logical :: success
  success = .true.

  !-----------------------------------------------------------------------------
  ! TEST #1: Dirichlet, solver should be exact
  !-----------------------------------------------------------------------------

  test_case => test_case_dirichlet
  nr        = 64
  nth       = 32
  tol       = 1.0e-11_f64
  rgrid_opt = "greville"

  call run_test( test_case, nr, nth, rgrid_opt, error_norm )

  ! Write relative error norm (global) to standard output
  write(*,"(/a)") "------------------------------------------------------------"
  write(*,"(a)")  "BC at r_min: homogeneous Dirichlet"
  write(*,"(a)")  "BC at r_max: homogeneous Dirichlet"
  write(*,"(a)")  "Radial grid: "// trim( rgrid_opt )
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

  test_case => test_case_neumann_mode0
  nr        = 64
  nth       = 32
  tol       = 1.0e-11_f64
  rgrid_opt = "smooth"

  call run_test( test_case, nr, nth, rgrid_opt, error_norm )

  ! Write relative error norm (global) to standard output
  write(*,"(/a)") "-----------------------------------------------------------&
       &--------------------"
  write(*,"(a)")  "BC at r_min: Neumann-mode-0"
  write(*,"(a)")  "BC at r_max: homogeneous Dirichlet"
  write(*,"(a)")  "Radial grid: "// trim( rgrid_opt )
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

  subroutine run_test( test_case, nr, nth, rgrid_opt, error_norm )
    class(c_test_qn_solver_2d_polar_base), intent(in   ) :: test_case
    sll_int32                            , intent(in   ) :: nr
    sll_int32                            , intent(in   ) :: nth
    character(len=*)                     , intent(in   ) :: rgrid_opt
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

    sll_real64, allocatable :: rgrid(:)

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

    ! Computational grid in theta
    dth = sll_p_twopi / nth

    ! Computational grid in r
    allocate( rgrid(nr+1) )
    select case (rgrid_opt)

    case ("uniform") ! uniform grid
      call sll_s_new_array_linspace( rgrid, rlim(1), rlim(2), endpoint=.true. )

    case ("smooth")  ! apply smooth coordinate transformation to uniform grid
      associate( alpha => 0.3_f64 )
        !
        ! 1. Create uniform logical grid: $\eta \in [0,1]$;
        ! 2. Apply sine transformation: 1D version of 2D transformation in
        !    P. Colella et al. JCP 230 (2011), formula (102) p. 2968;
        !    $\zeta = \eta + \alpha \sin(2\pi\eta)$, with $\zeta \in [0,1]$
        ! 3. Apply linear transformation to obtain radial grid:
        !    $r = r_{\min}(1-\zeta) + r_{\max}\zeta$, with $r \in [rmin,rmax]$.
        !
        call sll_s_new_array_linspace( rgrid, 0.0_f64, 1.0_f64, endpoint=.true. )
        rgrid = rgrid + alpha * sin( sll_p_twopi*rgrid ); rgrid(nr+1) = 1.0_f64
        rgrid = rlim(1)*(1.0_f64-rgrid) + rlim(2)*rgrid
      end associate

    case ("greville") ! similar to cubic spline with Greville's BCs
      associate( nc => nr-2 )
        dr = (rlim(2)-rlim(1))/ nc
        rgrid(1) = rlim(1)
        rgrid(2) = rlim(1) + dr/3.0_f64
        rgrid(3:nr-1) = [(rlim(1)+i*dr, i=1,nc-1)]
        rgrid(nr  ) = rlim(2) - dr/3.0_f64
        rgrid(nr+1) = rlim(2)
      end associate

    case ("default")
      SLL_ERROR("run_test","Unrecognized value for rgrid_option: "//trim(rgrid_opt))

    end select

    ! Allocate 2D distributed arrays (rho, phi, phi_ex) with layout_a
    allocate( rho   (nr+1,nth) )
    allocate( phi_ex(nr+1,nth) )
    allocate( phi   (nr+1,nth) )

    ! Load analytical solution and rho
    do j = 1, nth
      th = (j-1)*dth
      do i = 1, nr+1
        r = rgrid(i)
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
      r = rgrid(i)
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
      bc_rmax        = bcs(2), &
      rgrid          = rgrid )

    ! Compute numerical phi for a given rho
    call sll_s_qn_solver_2d_polar_solve( solver, rho, phi )

    ! Compute maximum norms of phi_ex and error
    max_phi = maxval(abs( phi_ex ))
    max_err = maxval(abs( phi_ex-phi ))

    ! Global relative error in maximum norm
    error_norm = max_err / max_phi

  end subroutine run_test

end program test_qn_solver_2d_polar
