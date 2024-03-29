program test_poisson_2d_fem_sps
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"

   use sll_m_working_precision, only: f64

   use sll_m_constants, only: sll_p_twopi

   use sll_m_utilities, only: sll_s_new_array_linspace

   use sll_m_boundary_condition_descriptors, only: &
      sll_p_greville, &
      sll_p_periodic

   use sll_m_bsplines, only: &
      sll_c_bsplines, &
      sll_s_bsplines_new

   use sll_m_spline_interpolator_1d, only: sll_s_spline_1d_compute_num_cells

   use sll_m_spline_2d, only: sll_t_spline_2d

   use sll_m_spline_interpolator_2d, only: sll_t_spline_interpolator_2d

   use sll_m_singular_mapping_analytic, only: sll_c_singular_mapping_analytic

   use sll_m_singular_mapping_analytic_target, only: sll_t_singular_mapping_analytic_target

   use sll_m_singular_mapping_analytic_czarny, only: sll_t_singular_mapping_analytic_czarny

   use sll_m_singular_mapping_discrete, only: sll_t_singular_mapping_discrete

   use sll_m_poisson_2d_fem_sps_dense, only: sll_t_poisson_2d_fem_sps_dense

   use sll_m_timer, only: &
      sll_t_time_mark, &
      sll_s_set_time_mark, &
      sll_f_time_elapsed_between

   use sll_m_hdf5_io_serial, only: &
      sll_t_hdf5_ser_handle, &
      sll_s_hdf5_ser_file_create, &
      sll_s_hdf5_ser_file_close, &
      sll_o_hdf5_ser_write_array, &
      sll_o_hdf5_ser_write_attribute

   implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   ! Working precision
   integer, parameter :: wp = f64

   ! To initialize B-splines (p1,p2 degrees)
   integer :: mm, n1, n2, p1, p2, ncells1, ncells2

   ! B-splines break points
   real(wp), allocatable :: breaks_eta1(:)
   real(wp), allocatable :: breaks_eta2(:)

   ! 1D B-splines
   class(sll_c_bsplines), allocatable :: bsplines_eta1
   class(sll_c_bsplines), allocatable :: bsplines_eta2

   ! Analytical and discrete mappings
   class(sll_c_singular_mapping_analytic), allocatable :: mapping_analytic
   type(sll_t_singular_mapping_discrete) :: mapping_discrete

   ! 2D splines representing right hand side and solution
   type(sll_t_spline_2d) :: spline_2d_rhs
   type(sll_t_spline_2d) :: spline_2d_phi

   ! 2D spline interpolator
   type(sll_t_spline_interpolator_2d) :: spline_interp_2d

   ! Needed for 2D interpolation of right hand side
   real(wp), allocatable :: tau_eta1(:)
   real(wp), allocatable :: tau_eta2(:)
   real(wp), allocatable :: gtau(:, :)

   ! Poisson solver
   type(sll_t_poisson_2d_fem_sps_dense) :: solver

   ! Auxiliary variables
   integer  :: i1, i2
   real(wp) :: eta(2), x(2)

   ! Timing
   type(sll_t_time_mark) :: t0, t1
   real(wp) :: dt

   ! For hdf5 I/O
   type(sll_t_hdf5_ser_handle) :: file_id
   integer                     :: h5_error

   !-----------------------------------------------------------------------------
   ! Initialize B-splines basis functions
   !-----------------------------------------------------------------------------

   call sll_s_set_time_mark(t0)

   ! Number of degrees of freedom (control points) along s and theta
   mm = 32
   n1 = mm*1
   n2 = mm*2

   ! Spline degrees along s and theta
   p1 = 3
   p2 = 5

   ! Compute number of cells from number of interpolation points along s
   call sll_s_spline_1d_compute_num_cells( &
      degree=p1, &
      bc_xmin=sll_p_greville, &
      bc_xmax=sll_p_greville, &
      nipts=n1, &
      ncells=ncells1)

   ! Construct break points along s to initialize non-uniform spline basis
   allocate (breaks_eta1(ncells1 + 1))
   call sll_s_new_array_linspace(breaks_eta1, 0.0_wp, 1.0_wp, endpoint=.true.)

   ! Create 1D spline basis along s in [0,1]
   call sll_s_bsplines_new( &
      bsplines=bsplines_eta1, &
      degree=p1, &
      periodic=.false., &
      xmin=0.0_wp, &
      xmax=1.0_wp, &
      ncells=ncells1, &
      breaks=breaks_eta1)

   ! Compute number of cells from number of interpolation points along theta
   call sll_s_spline_1d_compute_num_cells( &
      degree=p2, &
      bc_xmin=sll_p_periodic, &
      bc_xmax=sll_p_periodic, &
      nipts=n2, &
      ncells=ncells2)

   ! Construct break points along theta
   allocate (breaks_eta2(ncells2 + 1))
   call sll_s_new_array_linspace(breaks_eta2, 0.0_wp, sll_p_twopi, endpoint=.true.)

   ! Create 1D spline basis along theta in [0,2pi]
   call sll_s_bsplines_new( &
      bsplines=bsplines_eta2, &
      degree=p2, &
      periodic=.true., &
      xmin=0.0_wp, &
      xmax=sll_p_twopi, &
      ncells=ncells2)

   !-----------------------------------------------------------------------------
   ! Initialize mapping and polar B-splines
   !-----------------------------------------------------------------------------

   allocate (sll_t_singular_mapping_analytic_czarny :: mapping_analytic)

   ! Analytical mapping
   select type (mapping_analytic)
   type is (sll_t_singular_mapping_analytic_target)
      call mapping_analytic%init(x0=[0.0_wp, 0.0_wp], Delta=0.2_wp, kappa=0.3_wp)
   type is (sll_t_singular_mapping_analytic_czarny)
      call mapping_analytic%init(x0=[0.0_wp, 0.0_wp], e=1.4_wp, eps=0.3_wp)
   end select

   ! Discrete mapping
   call mapping_discrete%init(bsplines_eta1, bsplines_eta2, mapping_analytic)

   call sll_s_set_time_mark(t1)

   dt = sll_f_time_elapsed_between(t0, t1)
   write (*, '(/a,es8.1/)') " Time required for initialization of B-splines and mapping: ", dt

   !-----------------------------------------------------------------------------
   ! Interpolate right hand side in tensor-product space
   !-----------------------------------------------------------------------------

   call sll_s_set_time_mark(t0)

   ! Initialize 2D spline for right hand side
   call spline_2d_rhs%init(bsplines_eta1, bsplines_eta2)

   ! Initialize 2D spline interpolator
   call spline_interp_2d%init( &
      bsplines_eta1, &
      bsplines_eta2, &
      [sll_p_greville, sll_p_periodic], &
      [sll_p_greville, sll_p_periodic])

   ! Get interpolation points from 2D spline interpolator
   call spline_interp_2d%get_interp_points(tau_eta1, tau_eta2)

   associate (nt1 => size(tau_eta1), nt2 => size(tau_eta2))

      allocate (gtau(nt1, nt2))

      ! Evaluate right hand side on interpolation points
      do i2 = 1, nt2
         do i1 = 1, nt1
            eta = (/tau_eta1(i1), tau_eta2(i2)/)
            x = mapping_discrete%eval(eta)
            gtau(i1, i2) = rhs(x)
         end do
      end do

   end associate

   ! Compute interpolant spline
   call spline_interp_2d%compute_interpolant(spline_2d_rhs, gtau)

   call sll_s_set_time_mark(t1)

   dt = sll_f_time_elapsed_between(t0, t1)
   write (*, '(a,es8.1/)') " Time required for interpolation of right hand side: ", dt

   !-----------------------------------------------------------------------------
   ! Poisson solver
   !-----------------------------------------------------------------------------

   ! Initialize 2D spline for solution
   call spline_2d_phi%init(bsplines_eta1, bsplines_eta2)

   call sll_s_set_time_mark(t0)

   ! Initialize Poisson solver
   call solver%init(bsplines_eta1, bsplines_eta2, breaks_eta1, breaks_eta2, mapping_discrete)

   call sll_s_set_time_mark(t1)

   dt = sll_f_time_elapsed_between(t0, t1)
   write (*, '(a,es8.1/)') " Time required for initialization of Poisson solver: ", dt

   call sll_s_set_time_mark(t0)

   ! Solve
!  call solver % solve( spline_2d_rhs, spline_2d_phi ) ! Right hand side is 2D spline
   call solver%solve(rhs, spline_2d_phi) ! Right hand side is callable function

   call sll_s_set_time_mark(t1)

   dt = sll_f_time_elapsed_between(t0, t1)
   write (*, '(a,es8.1/)') " Time required for solution of Poisson equation: ", dt

   !-----------------------------------------------------------------------------
   ! HDF5 I/O
   !-----------------------------------------------------------------------------

   call sll_s_set_time_mark(t0)

   ! Create HDF5 file for output
   call sll_s_hdf5_ser_file_create('poisson_2d_fem_sps.h5', file_id, h5_error)

   ! Write useful parameters
   call sll_o_hdf5_ser_write_attribute(file_id, "/", "n1", n1, h5_error)
   call sll_o_hdf5_ser_write_attribute(file_id, "/", "n2", n2, h5_error)
   call sll_o_hdf5_ser_write_attribute(file_id, "/", "p1", p1, h5_error)
   call sll_o_hdf5_ser_write_attribute(file_id, "/", "p2", p2, h5_error)

   ! Write stiffness matrix
   call sll_o_hdf5_ser_write_array(file_id, solver%A, "/A", h5_error)

   ! Write mass matrix
   call sll_o_hdf5_ser_write_array(file_id, solver%M, "/M", h5_error)

   ! Write right hand side
   call sll_o_hdf5_ser_write_array(file_id, solver%b, "/b", h5_error)

   ! Write solution
   call sll_o_hdf5_ser_write_array(file_id, solver%x, "/x", h5_error)

   ! Write C1 projection of stiffness matrix
   call sll_o_hdf5_ser_write_array(file_id, solver%Ap, "/Ap", h5_error)

   ! Write C1 projection of mass matrix
   call sll_o_hdf5_ser_write_array(file_id, solver%Mp, "/Mp", h5_error)

!  ! Write L matrix needed for projection
!  call sll_o_hdf5_ser_write_array( file_id, solver % L, "/L", h5_error )

   ! Write reshaped solution
   call sll_o_hdf5_ser_write_array(file_id, spline_2d_phi%bcoef, "/phi", h5_error)

   ! Close HDF5 file
   call sll_s_hdf5_ser_file_close(file_id, h5_error)

   call sll_s_set_time_mark(t1)

   dt = sll_f_time_elapsed_between(t0, t1)
   write (*, '(a,es8.1/)') " Time required for writing HDF5 output: ", dt

   !-----------------------------------------------------------------------------
   ! Deallocations and free
   !-----------------------------------------------------------------------------

   deallocate (breaks_eta1)
   deallocate (breaks_eta2)

   deallocate (mapping_analytic)

   call bsplines_eta1%free()
   call bsplines_eta2%free()

   call mapping_discrete%free()

   call spline_2d_rhs%free()
   call spline_2d_phi%free()

   call solver%free()

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   SLL_PURE function rhs(x)
      real(wp), intent(in) :: x(2)
      real(wp) :: rhs

      rhs = sin(sll_p_twopi*x(1))*cos(sll_p_twopi*x(2))

   end function rhs

end program test_poisson_2d_fem_sps
