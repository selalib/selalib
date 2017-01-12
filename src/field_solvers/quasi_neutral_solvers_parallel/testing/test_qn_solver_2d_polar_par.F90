program test_qn_solver_2d_polar_par
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use iso_fortran_env, only: &
    output_unit

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_qn_solver_2d_polar_par, only: &
    sll_t_qn_solver_2d_polar_par, &
    sll_s_qn_solver_2d_polar_par_init, &
    sll_s_qn_solver_2d_polar_par_solve

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

  use sll_m_collective, only: &
    sll_t_collective_t, &
    sll_s_boot_collective, &
    sll_f_get_collective_rank, &
    sll_f_get_collective_size, &
    sll_o_collective_gather, &
    sll_o_collective_bcast, &
    sll_s_collective_barrier, &
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

  type(t_test_dirichlet_zero_error) :: test_case_dirichlet
  type(t_test_neumann_mode0_zero_error) :: test_case_neumann_mode0

  ! For MPI
  type(sll_t_collective_t), pointer :: comm
  sll_int32  :: my_rank
  sll_int32  :: nr
  sll_int32  :: nth
  sll_real64 :: error_norm
  sll_real64 :: tol

  logical :: success

  call sll_s_boot_collective()
  comm => sll_v_world_collective
  my_rank = sll_f_get_collective_rank( comm )

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

  call run_test( comm, test_case_dirichlet, nr, nth, error_norm )

  ! Write relative error norm (global) to standard output
  if (my_rank == 0) then
    write(*,"(/a)") "--- TEST 1 ---"
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

  call run_test( comm, test_case_neumann_mode0, nr, nth, error_norm )

  ! Write relative error norm (global) to standard output
  if (my_rank == 0) then
    write(*,"(/a)") "--- TEST 2 ---"
    write(*,"(a,e11.3)") "Relative L_inf norm of error = ", error_norm
    write(*,"(a,e11.3)") "Tolerance                    = ", tol
    if (error_norm > tol) then
      success = .false.
      write(*,"(a/)") "!!! FAILED !!!"
    end if
  end if

  ! TODO: run convergence analysis on a more difficult test-case

  ! Check if test passed
  if (my_rank == 0) then
    if(success) then
       write(*,"(/a/)") "--- PASSED ---"
    endif
  endif

  call sll_s_halt_collective()

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
contains
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  subroutine run_test( comm, test_case, nr, nth, error_norm )
    type(sll_t_collective_t)             , pointer       :: comm
    class(c_test_case_qn_solver_2d_polar), intent(in   ) :: test_case
    sll_int32                            , intent(in   ) :: nr
    sll_int32                            , intent(in   ) :: nth
    sll_real64                           , intent(  out) :: error_norm

    type(sll_t_qn_solver_2d_polar_par) :: solver
    sll_real64                         :: max_phi
    sll_real64                         :: max_err

    sll_real64 ::  r,  th
    sll_real64 :: dr, dth
    sll_int32  :: i, j
    sll_int32  :: num_proc
    sll_int32  :: loc_sz_r(2)
    sll_int32  :: loc_sz_a(2)
    sll_int32  :: glob_idx(2)

    sll_real64, allocatable :: rhs   (:,:)
    sll_real64, allocatable :: phi_ex(:,:)
    sll_real64, allocatable :: phi   (:,:)

    sll_real64, allocatable :: rho_m0(:)
    sll_real64, allocatable :: b_magn(:)
    sll_real64, allocatable :: lambda(:)

    type(sll_t_layout_2d), pointer :: layout_r
    type(sll_t_layout_2d), pointer :: layout_a

    ! Computational grid
    dr  = (test_case%rmax - test_case%rmin) / nr
    dth = 2.0_f64*sll_p_pi / nth

    ! Get number of available processes
    num_proc = sll_f_get_collective_size( comm )

    ! Create 2D layout sequential in r (distributed along theta)
    layout_r => sll_f_new_layout_2d( comm )
    call sll_o_initialize_layout_with_distributed_array( &
      nr+1, & 
      nth , & 
      1, &
      num_proc, &
      layout_r )

    ! Create 2D layout sequential in theta (distributed along r)
    layout_a => sll_f_new_layout_2d( comm )
    call sll_o_initialize_layout_with_distributed_array( &
      nr+1, & 
      nth , & 
      num_proc, &
      1, &
      layout_a )

    ! Compute local size of 2D arrays in the two layouts
    call sll_o_compute_local_sizes( layout_r, loc_sz_r(1), loc_sz_r(2) )
    call sll_o_compute_local_sizes( layout_a, loc_sz_a(1), loc_sz_a(2) )

    ! Allocate 2D distributed arrays (rho, phi, phi_ex) with layout_a
    allocate( rhs   (loc_sz_a(1),loc_sz_a(2)) )
    allocate( phi_ex(loc_sz_a(1),loc_sz_a(2)) )
    allocate( phi   (loc_sz_a(1),loc_sz_a(2)) )

    ! Load analytical solution and rhs
    do j = 1, loc_sz_a(2)
      th = (j-1)*dth
      do i = 1, loc_sz_a(1)
        glob_idx(:) = sll_o_local_to_global( layout_a, [i,j] )
        r = test_case%rmin + (glob_idx(1)-1)*dr
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

    ! Initialize parallel solver
    call sll_s_qn_solver_2d_polar_par_init( solver, &
      layout_r = layout_r, &
      layout_a = layout_a, &
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
    call sll_s_qn_solver_2d_polar_par_solve( solver, rhs, phi )

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

    sll_int32               :: np
    sll_int32               :: my_rank
    sll_real64              :: send_buf(1)
    sll_real64, allocatable :: recv_buf(:)

    ! Get information about parallel job
    np      = sll_f_get_collective_size( comm )
    my_rank = sll_f_get_collective_rank( comm )

    ! Write in/out variable to sender buffer
    send_buf(1) = v

    ! [ROOT only] Prepare receiver buffer
    if (my_rank == 0) then
      allocate( recv_buf(np) )
    end if

    ! Send v values to ROOT
    call sll_o_collective_gather( comm, send_buf, 1, 0, recv_buf )

    ! [ROOT only] Compute maximum of v values and prepare send buffer
    if (my_rank == 0) then
      send_buf(1) = maxval( recv_buf(:) )
      deallocate( recv_buf )
    end if

    ! Send maximum to all processes
    call sll_o_collective_bcast( comm, send_buf, 1, 0 )

    ! Write result to in/out variable
    v = send_buf(1)

  end subroutine s_compute_collective_max
  !-----------------------------------------------------------------------------


end program test_qn_solver_2d_polar_par
