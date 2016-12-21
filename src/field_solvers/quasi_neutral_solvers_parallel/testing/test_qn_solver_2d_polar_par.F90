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

  use m_test_case_2d_dirichlet_1, only: &
    t_test_dirichlet_zero_error

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_dirichlet, &
    sll_p_neumann, &
    sll_p_neumann_mode_0

  use sll_m_collective, only: &
    sll_s_boot_collective, &
    sll_f_get_collective_rank, &
    sll_f_get_collective_size, &
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

  call test_dirichlet_zero_error()

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
contains
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  subroutine test_dirichlet_zero_error()

    type(sll_t_qn_solver_2d_polar_par) :: solver
    type(t_test_dirichlet_zero_error)  :: test_case

    sll_int32  :: nr, nth
    sll_real64 ::  r,  th
    sll_real64 :: dr, dth
    sll_int32  :: i, j
    sll_int32  :: num_proc
    sll_int32  :: my_rank
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

    ! Define test case
    test_case%rmin                = 1.0_f64
    test_case%rmax                = 10.0_f64
    test_case%adiabatic_electrons = .true.
    test_case%use_zonal_flow      = .false.
    test_case%epsilon_0           = 1.0_f64
    test_case%bc_rmin             = sll_p_dirichlet
    test_case%bc_rmax             = sll_p_dirichlet

    ! Computational grid
    nr  = 256
    nth = 32
    dr  = (test_case%rmax - test_case%rmin) / nr
    dth = 2.0_f64*sll_p_pi / nth

    call sll_s_boot_collective()

    num_proc = sll_f_get_collective_size(sll_v_world_collective)  

    ! Create 2D layout sequential in r (distributed along theta)
    layout_r => sll_f_new_layout_2d( sll_v_world_collective )
    call sll_o_initialize_layout_with_distributed_array( &
      nr+1, & 
      nth , & 
      1, &
      num_proc, &
      layout_r )

    ! Create 2D layout sequential in theta (distributed along r)
    layout_a => sll_f_new_layout_2d( sll_v_world_collective )
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

    ! Allocate 1D radial profiles (needed by solver)
    allocate( rho_m0(nr+1) )
    allocate( b_magn(nr+1) )
    allocate( lambda(nr+1) )

    do i = 1, nr+1
      r = test_case%rmin + (i-1)*dr
      rho_m0(i) = test_case%rho_m0( r )
      b_magn(i) = test_case%b_magn( r )
      lambda(i) = test_case%lambda( r )
    end do

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

    call sll_s_qn_solver_2d_polar_par_solve( solver, rhs, phi )

    my_rank = sll_f_get_collective_rank( sll_v_world_collective )
    do i = 0, num_proc-1
      if (i == my_rank) then
        write(*,"(a,i0)") "processor # ", my_rank
        write(*,"(a,e15.5)") "max error = ", maxval(abs( phi-phi_ex ))
        write(*,*)
        flush( output_unit )
      end if
      call sll_s_collective_barrier( sll_v_world_collective )
    end do

    call sll_s_halt_collective()

  end subroutine test_dirichlet_zero_error

end program test_qn_solver_2d_polar_par
