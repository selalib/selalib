!> @authors Yaman Güçlü, IPP Garching
!> @authors Edoardo Zoni, IPP Garching
!>
!> @brief
!> Program to run tests of module 'sll_m_qn_solver_2d_polar_par.F90'.
!>
!> @details
!> Default parameters of each test type can be overwritten (except for BCs) before
!> call to subroutine 'run_test'.
!> Tolerance is set according to the order of the radial profile of each test type
!> (expected zero or non-zero numerical error) as well as according to the mesh sizes.
!> Recall that k (Fourier mode, see default parameters) must be <= ntheta/2.

program test_qn_solver_2d_polar_par
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"
#include "sll_errors.h"

  use iso_fortran_env, only: &
    output_unit

  use sll_m_constants, only: &
    sll_p_twopi

  use sll_m_utilities, only: &
    sll_s_new_array_linspace

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_polar_origin

  use sll_m_qn_solver_2d_polar_par, only: &
    sll_t_qn_solver_2d_polar_par, &
    sll_s_qn_solver_2d_polar_par_init, &
    sll_s_qn_solver_2d_polar_par_solve

  use m_test_qn_solver_2d_polar_base, only: &
    c_test_qn_solver_2d_polar_base

  use m_test_qn_solver_2d_polar_annulus_dirichlet, only: &
    t_test_qn_solver_2d_polar_annulus_dirichlet_quadratic

  use m_test_qn_solver_2d_polar_annulus_neumann_mode0, only: &
    t_test_qn_solver_2d_polar_annulus_neumann_mode0_quadratic

  use m_test_qn_solver_2d_polar_disk_dirichlet, only: &
    t_test_qn_solver_2d_polar_disk_dirichlet_quadratic

  use sll_mpi, only: &
    mpi_max

  use sll_m_collective, only: &
    sll_t_collective_t, &
    sll_s_boot_collective, &
    sll_f_get_collective_rank, &
    sll_f_get_collective_size, &
    sll_o_collective_allreduce, &
    sll_s_halt_collective, &
    sll_v_world_collective

  use sll_m_remapper, only: &
    sll_o_initialize_layout_with_distributed_array, &
    sll_o_compute_local_sizes, &
    sll_o_local_to_global, &
    sll_f_new_layout_2d, &
    sll_t_layout_2d

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  class(c_test_qn_solver_2d_polar_base)                           , pointer :: test_case
  type (t_test_qn_solver_2d_polar_annulus_dirichlet_quadratic    ), target  :: test_case_dirichlet
  type (t_test_qn_solver_2d_polar_annulus_neumann_mode0_quadratic), target  :: test_case_neumann_mode0
  type (t_test_qn_solver_2d_polar_disk_dirichlet_quadratic       ), target  :: test_case_circle_dirichlet

  type(sll_t_collective_t), pointer :: comm
  sll_int32  :: my_rank

  sll_int32  :: nr, nth
  sll_real64 :: error_norm, tol
  character(len=8) :: rgrid_opt

  logical :: success
  success = .true.

  call sll_s_boot_collective()
  comm => sll_v_world_collective
  my_rank = sll_f_get_collective_rank( comm )

  !-----------------------------------------------------------------------------
  ! TEST #1: Dirichlet, solver should be exact
  !-----------------------------------------------------------------------------

  test_case => test_case_dirichlet
  nr        = 64
  nth       = 32
  tol       = 1.0e-11_f64
  rgrid_opt = "greville"

  call run_test( comm, test_case, nr, nth, rgrid_opt, error_norm )

  ! Write relative error norm (global) to standard output
  if (my_rank == 0) then
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
  end if

  !-----------------------------------------------------------------------------
  ! TEST #2: Neumann mode 0, solver should be exact
  !-----------------------------------------------------------------------------

  test_case => test_case_neumann_mode0
  nr        = 64
  nth       = 32
  tol       = 1.0e-11_f64
  rgrid_opt = "smooth"

  call run_test( comm, test_case, nr, nth, rgrid_opt, error_norm )

  ! Write relative error norm (global) to standard output
  if (my_rank == 0) then
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
  end if

  !-----------------------------------------------------------------------------
  ! TEST #3: Full circle (rmin=0), Dirichlet at rmax, solver should be exact
  !-----------------------------------------------------------------------------

  test_case => test_case_circle_dirichlet
  nr        = 64
  nth       = 32
  tol       = 1.0e-11_f64
  rgrid_opt = "uniform"

  call run_test( comm, test_case, nr, nth, rgrid_opt, error_norm )

  ! Write relative error norm (global) to standard output
  if (my_rank == 0) then
    write(*,"(/a)") "-----------------------------------------------------------&
         &--------------------"
    write(*,"(a)")  "BC at r_min: polar origin (i.e., full circle is simulated)"
    write(*,"(a)")  "BC at r_max: homogeneous Dirichlet"
    write(*,"(a)")  "Radial grid: "// trim( rgrid_opt )
    write(*,"(a)")  "phi(r,theta) = a (1-(r/rmax)^2) &
         &+ b 4(r/rmax)(1-r/rmax)cos(k(theta-theta_0))"
    write(*,"(a)")  "-----------------------------------------------------------&
         &--------------------"
    write(*,"(a,e11.3)") "Relative L_inf norm of error = ", error_norm
    write(*,"(a,e11.3)") "Tolerance                    = ", tol
    if (error_norm > tol) then
       success = .false.
       write(*,"(a/)") "!!! FAILED !!!"
    end if
  endif

  !-----------------------------------------------------------------------------
  ! Check if all tests have passed
  !-----------------------------------------------------------------------------
  if (my_rank == 0) then
    if(success) then
       write(*,"(/a/)") "PASSED"
    endif
  endif

  call sll_s_halt_collective()

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
contains
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  subroutine run_test( comm, test_case, nr, nth, rgrid_opt, error_norm )
    type(sll_t_collective_t)             , pointer       :: comm
    class(c_test_qn_solver_2d_polar_base), intent(in   ) :: test_case
    sll_int32                            , intent(in   ) :: nr
    sll_int32                            , intent(in   ) :: nth
    character(len=*)                     , intent(in   ) :: rgrid_opt
    sll_real64                           , intent(  out) :: error_norm

    type(sll_t_qn_solver_2d_polar_par) :: solver
    sll_real64                         :: max_phi
    sll_real64                         :: max_err

    sll_real64 ::  r,  th
    sll_real64 :: dr, dth
    sll_int32  :: i, j
    sll_int32  :: sh

    sll_real64 :: rlim(2)
    sll_int32  :: bcs (2)
    logical    :: adiabatic_electrons
    logical    :: use_zonal_flow
    sll_real64 :: epsilon_0

    sll_int32  :: num_proc
    sll_int32  :: loc_sz_r(2)
    sll_int32  :: loc_sz_a(2)
    sll_int32  :: glob_idx(2)

    sll_real64, allocatable :: rgrid(:)

    sll_real64, allocatable :: rho   (:,:)
    sll_real64, allocatable :: phi_ex(:,:)
    sll_real64, allocatable :: phi   (:,:)

    sll_real64, allocatable :: rho_m0(:)
    sll_real64, allocatable :: b_magn(:)
    sll_real64, allocatable :: lambda(:)

    type(sll_t_layout_2d), pointer :: layout_r
    type(sll_t_layout_2d), pointer :: layout_a

    ! Extract domain limits and boundary conditions
    rlim(:) = test_case%get_rlim()
    bcs (:) = test_case%get_bcs ()

    ! Extract test-case parameters
    call test_case%get_parameters( adiabatic_electrons, use_zonal_flow, epsilon_0 )

    ! Computational grid in theta
    dth = sll_p_twopi / nth

    ! Computational grid in r: handle full circle
    sh = merge( 1, 0, bcs(1) == sll_p_polar_origin )
    allocate( rgrid(nr+1-sh) )

    if (bcs(1) == sll_p_polar_origin) then
      associate( rmin => rlim(2)/real(2*nr-1,f64) )
        call sll_s_new_array_linspace( rgrid, rmin, rlim(2), endpoint=.true. )
      end associate

    else

      ! Computational grid in r: handle non-uniform spacing
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

    end if

    ! Get number of available processes
    num_proc = sll_f_get_collective_size( comm )

    ! Create 2D layout sequential in r (distributed along theta)
    layout_r => sll_f_new_layout_2d( comm )
    call sll_o_initialize_layout_with_distributed_array( &
      nr+1-sh, &
      nth , & 
      1, &
      num_proc, &
      layout_r )

    ! Create 2D layout sequential in theta (distributed along r)
    layout_a => sll_f_new_layout_2d( comm )
    call sll_o_initialize_layout_with_distributed_array( &
      nr+1-sh, &
      nth , & 
      num_proc, &
      1, &
      layout_a )

    ! Compute local size of 2D arrays in the two layouts
    call sll_o_compute_local_sizes( layout_r, loc_sz_r(1), loc_sz_r(2) )
    call sll_o_compute_local_sizes( layout_a, loc_sz_a(1), loc_sz_a(2) )

    ! Allocate 2D distributed arrays (rho, phi, phi_ex) with layout_a
    allocate( rho   (loc_sz_a(1),loc_sz_a(2)) )
    allocate( phi_ex(loc_sz_a(1),loc_sz_a(2)) )
    allocate( phi   (loc_sz_a(1),loc_sz_a(2)) )

    ! Load analytical solution and rho
    do j = 1, loc_sz_a(2)
      th = (j-1)*dth
      do i = 1, loc_sz_a(1)
        glob_idx(:) = sll_o_local_to_global( layout_a, [i,j] )
        r = rgrid(glob_idx(1))
        phi_ex(i,j) = test_case%phi_ex( r, th )
        rho   (i,j) = test_case%rho   ( r, th )
      end do
    end do
    phi(:,:) = 0.0_f64

    ! Equation parameters: allocate and load rho_m0(r) and b_magn(r)
    allocate( rho_m0(nr+1-sh) )
    allocate( b_magn(nr+1-sh) )
    do i = 1, nr+1-sh
      r = rgrid(i)
      rho_m0(i) = test_case%rho_m0( r )
      b_magn(i) = test_case%b_magn( r )
    end do

    ! Equation parameters: if required, also allocate and load b_magn(r)
    if (adiabatic_electrons) then
      allocate( lambda(nr+1-sh) )
      do i = 1, nr+1-sh
        r = rgrid(i)
        lambda(i) = test_case%lambda( r )
      end do
    end if

    ! Initialize parallel solver
    call sll_s_qn_solver_2d_polar_par_init( solver, &
      layout_r       = layout_r, &
      layout_a       = layout_a, &
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
    call sll_s_qn_solver_2d_polar_par_solve( solver, rho, phi )

    ! Compute (local) maximum norms of phi_ex and error
    max_phi = maxval(abs( phi_ex ))
    max_err = maxval(abs( phi_ex-phi ))

    ! Exchange norms across processes to obtain global norms
    call s_compute_collective_max( comm, max_phi )
    call s_compute_collective_max( comm, max_err )

    ! Global relative error in maximum norm
    error_norm = max_err / max_phi

  end subroutine run_test

  !-----------------------------------------------------------------------------
  subroutine s_compute_collective_max( comm, v )
    type(sll_t_collective_t), pointer       :: comm
    sll_real64              , intent(inout) :: v

    sll_real64 :: send_buf(1)
    sll_real64 :: recv_buf(1)

    ! Write in/out variable to sender buffer
    send_buf(1) = v

    ! Compute maximum of all v values using MPI_ALLREDUCE with MPI_MAX operation
    call sll_o_collective_allreduce( comm, send_buf, 1, mpi_max, recv_buf )

    ! Write result to in/out variable
    v = recv_buf(1)

  end subroutine s_compute_collective_max

end program test_qn_solver_2d_polar_par
