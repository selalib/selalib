program test_poisson_2d_polar
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use iso_fortran_env, only: &
    output_unit

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_poisson_2d_polar, only: &
    sll_t_poisson_2d_polar, &
    sll_s_poisson_2d_polar_init, &
    sll_s_poisson_2d_polar_solve

  use m_test_case_poisson_2d_base, only: &
    c_test_case_poisson_2d_polar

  use m_test_case_poisson_2d_dirichlet, only: &
    t_test_dirichlet_zero_error, &
    t_test_dirichlet

  use m_test_case_poisson_2d_neumann_mode0, only: &
    t_test_neumann_mode0_zero_error

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_dirichlet, &
    sll_p_neumann, &
    sll_p_neumann_mode_0

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type(t_test_dirichlet_zero_error) :: test_case_dirichlet_zero_error
  type(t_test_dirichlet) :: test_case_dirichlet
  type(t_test_neumann_mode0_zero_error) :: test_case_neumann_mode0_zero_error

  sll_int32  :: nr
  sll_int32  :: nth
  sll_real64 :: error_norm
  sll_real64 :: tol

  logical :: success

  success = .true.

  !=============================================================================
  ! TEST #1: Dirichlet, solver should be exact
  !=============================================================================
  nr  = 64
  nth = 32
  tol = 1.0e-11_f64

  ! Define test case
  test_case_dirichlet_zero_error%rmin = 1.0_f64
  test_case_dirichlet_zero_error%rmax = 2.0_f64
  test_case_dirichlet_zero_error%bc_rmin = sll_p_dirichlet
  test_case_dirichlet_zero_error%bc_rmax = sll_p_dirichlet

  call run_test( test_case_dirichlet_zero_error, nr, nth, error_norm )

  ! Write relative error norm (global) to standard output
  write(*,"(/a)") "-------------------------------------------------------"
  write(*,"(a)")  "Homogeneous Dirichlet boundary conditions"
  write(*,"(a)")  "phi(r,theta) = (r-rmin)(r-rmax)(a+sin(k(theta-theta_0))"
  write(*,"(a)")  "-------------------------------------------------------"
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

  ! Define test case
  test_case_dirichlet%rmin = 1.0_f64
  test_case_dirichlet%rmax = 2.0_f64
  test_case_dirichlet%bc_rmin = sll_p_dirichlet
  test_case_dirichlet%bc_rmax = sll_p_dirichlet

  call run_test( test_case_dirichlet, nr, nth, error_norm )

  ! Write relative error norm (global) to standard output
  write(*,"(/a)") "--------------------------------------------------------"
  write(*,"(a)")  "Homogeneous Dirichlet boundary conditions               "
  write(*,"(a)")  "phi(r,theta) = r(r-rmin)(r-rmax)(a+sin(k(theta-theta_0))"
  write(*,"(a)")  "--------------------------------------------------------"
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

  ! Define test case
  test_case_neumann_mode0_zero_error%rmin = 1.0_f64
  test_case_neumann_mode0_zero_error%rmax = 2.0_f64
  test_case_neumann_mode0_zero_error%bc_rmin = sll_p_neumann_mode_0
  test_case_neumann_mode0_zero_error%bc_rmax = sll_p_dirichlet

  call run_test( test_case_neumann_mode0_zero_error, nr, nth, error_norm )

  ! Write relative error norm (global) to standard output
  write(*,"(/a)") "-----------------------------------------------------------&
       &--------------------"
  write(*,"(a)")  "Mixed Homogeneous Dirichlet / Neumann mode 0 boundary conditions"
  write(*,"(a)")  "phi(r,theta) = a(r-rmax)(r+rmax-2rmin) &
       &+ b(r-rmin)(r-rmax)sin(k(theta-theta_0))"
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
    class(c_test_case_poisson_2d_polar), intent(in) :: test_case
    sll_int32, intent(in) :: nr
    sll_int32, intent(in) :: nth
    sll_real64, intent(out) :: error_norm

    type(sll_t_poisson_2d_polar) :: solver
    sll_real64 :: max_phi
    sll_real64 :: max_err

    sll_real64 ::  r,  th
    sll_real64 :: dr, dth
    sll_int32 :: i, j

    sll_real64, allocatable :: rhs   (:,:)
    sll_real64, allocatable :: phi_ex(:,:)
    sll_real64, allocatable :: phi   (:,:)

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

    ! Initialize parallel solver
    call sll_s_poisson_2d_polar_init( solver, &
      rmin     = test_case%rmin, &
      rmax     = test_case%rmax, &
      nr       = nr, &
      ntheta   = nth, &
      bc_rmin  = test_case%bc_rmin , &
      bc_rmax  = test_case%bc_rmax )

    ! Compute numerical phi for a given rhs
    call sll_s_poisson_2d_polar_solve( solver, rhs, phi )

    ! Compute (local) maximum norms of phi_ex and error
    max_phi = maxval(abs( phi_ex ))
    max_err = maxval(abs( phi_ex-phi ))

    ! Global relative error in maximum norm
    error_norm = max_err / max_phi

  end subroutine run_test


end program test_poisson_2d_polar
