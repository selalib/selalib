!> @authors Yaman Güçlü, IPP Garching
!> @authors Edoardo Zoni, IPP Garching
!>
!> @brief
!> Program to run tests of module 'sll_m_poisson_2d_polar.F90'.
!>
!> @details
!> Default parameters of each test type can be overwritten (except for BCs) before
!> call to subroutine 'run_test'.
!> Tolerance is set according to the order of the radial profile of each test type
!> (expected zero or non-zero numerical error) as well as according to the mesh sizes.
!> Recall that k (Fourier mode, see default parameters) must be <= ntheta/2.
!>

program test_poisson_2d_polar
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_poisson_2d_polar, only: &
    sll_t_poisson_2d_polar, &
    sll_s_poisson_2d_polar_init, &
    sll_s_poisson_2d_polar_solve

  use m_test_poisson_2d_polar_base, only: &
    c_test_poisson_2d_polar_base

  use m_test_poisson_2d_polar_dirichlet, only: &
    t_test_poisson_2d_polar_dirichlet_quadratic, &
    t_test_poisson_2d_polar_dirichlet_cubic

  use m_test_poisson_2d_polar_neumann_mode0, only: &
    t_test_poisson_2d_polar_neumann_mode0_quadratic

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type(t_test_poisson_2d_polar_dirichlet_quadratic)     :: test_case_dirichlet_zero_error
  type(t_test_poisson_2d_polar_dirichlet_cubic)         :: test_case_dirichlet
  type(t_test_poisson_2d_polar_neumann_mode0_quadratic) :: test_case_neumann_mode0_zero_error

  sll_int32  :: nr, nth
  sll_real64 :: error_norm, tol

  logical :: success
  success = .true.

  !=============================================================================
  ! TEST #1: Dirichlet, solver should be exact
  !=============================================================================
  nr  = 64
  nth = 32
  tol = 1.0e-11_f64

  call run_test( test_case_dirichlet_zero_error, nr, nth, error_norm )

  ! Write relative error norm (global) to standard output
  write(*,"(/a)") "------------------------------------------------------------"
  write(*,"(a)")  "Homogeneous Dirichlet boundary conditions"
  write(*,"(a)")  "phi(r,theta) = (r-rmax)(r-rmin)(a + b*cos(k(theta-theta_0)))"
  write(*,"(a)")  "------------------------------------------------------------"
  write(*,"(a,e11.3)") "Relative L_inf norm of error = ", error_norm
  write(*,"(a,e11.3)") "Tolerance                    = ", tol
  if (error_norm > tol) then
     success = .false.
     write(*,"(a)") "!!! FAILED !!!"
  end if

  !=============================================================================
  ! TEST #2: Dirichlet, cubic profile
  !=============================================================================
  nr  = 64
  nth = 32
  tol = 1.0e-4_f64

  call run_test( test_case_dirichlet, nr, nth, error_norm )

  ! Write relative error norm (global) to standard output
  write(*,"(/a)") "-------------------------------------------------------------"
  write(*,"(a)")  "Homogeneous Dirichlet boundary conditions               "
  write(*,"(a)")  "phi(r,theta) = r(r-rmax)(r-rmin)(a + b*cos(k(theta-theta_0)))"
  write(*,"(a)")  "-------------------------------------------------------------"
  write(*,"(a,e11.3)") "Relative L_inf norm of error = ", error_norm
  write(*,"(a,e11.3)") "Tolerance                    = ", tol
  if (error_norm > tol) then
     success = .false.
     write(*,"(a)") "!!! FAILED !!!"
  end if

  !=============================================================================
  ! TEST #3: Neumann mode 0, solver should be exact
  !=============================================================================
  nr  = 64
  nth = 32
  tol = 1.0e-11_f64

  call run_test( test_case_neumann_mode0_zero_error, nr, nth, error_norm )

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
    class(c_test_poisson_2d_polar_base), intent(in   ) :: test_case
    sll_int32                          , intent(in   ) :: nr
    sll_int32                          , intent(in   ) :: nth
    sll_real64                         , intent(  out) :: error_norm

    type(sll_t_poisson_2d_polar) :: solver
    sll_real64 :: max_phi
    sll_real64 :: max_err

    sll_real64 :: rlim(2)
    sll_int32  ::  bcs(2)
    sll_real64 ::  r,  th
    sll_real64 :: dr, dth
    sll_int32  :: i, j

    sll_real64, allocatable :: rho   (:,:)
    sll_real64, allocatable :: phi_ex(:,:)
    sll_real64, allocatable :: phi   (:,:)

    ! Extract domain limits and boundary conditions
    rlim(:) = test_case%get_rlim()
    bcs (:) = test_case%get_bcs ()

    ! Computational grid
    dr  = (rlim(2)-rlim(1))/ nr
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

    ! Initialize parallel solver
    call sll_s_poisson_2d_polar_init( solver, &
      rmin     = rlim(1), &
      rmax     = rlim(2), &
      nr       = nr, &
      ntheta   = nth, &
      bc_rmin  = bcs(1), &
      bc_rmax  = bcs(2) )

    ! Compute numerical phi for a given rho
    call sll_s_poisson_2d_polar_solve( solver, rho, phi )

    ! Compute (local) maximum norms of phi_ex and error
    max_phi = maxval(abs( phi_ex ))
    max_err = maxval(abs( phi_ex-phi ))

    ! Global relative error in maximum norm
    error_norm = max_err / max_phi

  end subroutine run_test


end program test_poisson_2d_polar
