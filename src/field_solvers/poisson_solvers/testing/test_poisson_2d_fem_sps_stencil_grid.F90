program test_poisson_2d_fem_sps_stencil_grid
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

   use sll_m_spline_interpolator_2d, only: sll_t_spline_interpolator_2d

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

   ! Integer variables
   integer :: n1, n2, p1, p2, ncells1, ncells2, file_er, file_unit

   ! Character variables
   character(len=:), allocatable :: input_file

   ! Real 1D allocatables
   real(wp), allocatable :: breaks_eta1(:), breaks_eta2(:), tau_eta1(:), tau_eta2(:)

   ! Abstract polymorphic types
   class(sll_c_bsplines), allocatable :: bsplines_eta1, bsplines_eta2

   ! Concrete types
   type(sll_t_spline_interpolator_2d) :: spline_interp_2d
   type(sll_t_hdf5_ser_handle)        :: file_id

   ! Namelists in input file
   namelist /splines/ &
      n1, &
      n2, &
      p1, &
      p2

   ! Read input file
   input_file = trim("test_poisson_2d_fem_sps_stencil.nml")
   open (file=input_file, status='old', action='read', newunit=file_unit)
   read (file_unit, splines)
   close (file_unit)

   ! Create HDF5 file
   call sll_s_hdf5_ser_file_create('test_poisson_2d_fem_sps_stencil_grid.h5', file_id, file_er)

   ! HDF5 I/O
   call sll_o_hdf5_ser_write_attribute(file_id, "/", "n1", n1, file_er)
   call sll_o_hdf5_ser_write_attribute(file_id, "/", "n2", n2, file_er)
   call sll_o_hdf5_ser_write_attribute(file_id, "/", "p1", p1, file_er)
   call sll_o_hdf5_ser_write_attribute(file_id, "/", "p2", p2, file_er)

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

   ! Initialize 2D spline interpolator
   call spline_interp_2d%init( &
      bsplines_eta1, &
      bsplines_eta2, &
      [sll_p_greville, sll_p_periodic], &
      [sll_p_greville, sll_p_periodic])

   ! Get interpolation points from 2D spline interpolator
   call spline_interp_2d%get_interp_points(tau_eta1, tau_eta2)

   ! Write error
   call sll_o_hdf5_ser_write_array(file_id, tau_eta1, "/eta1_grid", file_er)
   call sll_o_hdf5_ser_write_array(file_id, tau_eta2, "/eta2_grid", file_er)

   ! Close HDF5 file
   call sll_s_hdf5_ser_file_close(file_id, file_er)

   deallocate (breaks_eta1, breaks_eta2)

   call bsplines_eta1%free()
   call bsplines_eta2%free()

end program test_poisson_2d_fem_sps_stencil_grid
