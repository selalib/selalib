program test_greville_points
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"

   use sll_m_working_precision, only: f64

   use sll_m_constants, only: &
      sll_p_pi, &
      sll_p_twopi

   use sll_m_utilities, only: sll_s_new_array_linspace

   use sll_m_boundary_condition_descriptors, only: &
      sll_p_greville, &
      sll_p_periodic

   use sll_m_bsplines, only: &
      sll_c_bsplines, &
      sll_s_bsplines_new

   use sll_m_spline_2d, only: sll_t_spline_2d

   use sll_m_spline_interpolator_1d, only: &
      sll_s_spline_1d_compute_num_cells, &
      sll_t_spline_interpolator_1d

   use sll_m_spline_interpolator_2d, only: sll_t_spline_interpolator_2d

   implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   ! Working precision
   integer, parameter :: wp = f64

   ! To initialize B-splines
   integer :: n1, n2, n4, deg1, deg2, deg4, ncells1, ncells2, ncells4

   ! To initialize non-uniform 1D B-splines along eta1
   real(wp), allocatable :: breaks_eta1(:)

   ! 1D B-splines
   class(sll_c_bsplines), allocatable :: spline_basis_eta1
   class(sll_c_bsplines), allocatable :: spline_basis_eta2
   class(sll_c_bsplines), allocatable :: spline_basis_eta4

   ! 1D/2D tensor-product spline interpolators
   type(sll_t_spline_interpolator_1d) :: spline_interpolator_1d
   type(sll_t_spline_interpolator_2d) :: spline_interpolator_2d

   ! Interpolation points and fields
   integer :: npts1, npts2, npts4
   real(wp), allocatable :: e1_node(:)
   real(wp), allocatable :: e2_node(:)
   real(wp), allocatable :: e4_node(:)

   ! Auxiliary variables
   integer  :: i1, i2, i4

   n1 = 8
   n2 = 8
   n4 = 8

   deg1 = 3
   deg2 = 3
   deg4 = 3

   ! Compute number of cells from number of interpolation points along eta1
   call sll_s_spline_1d_compute_num_cells( &
      degree=deg1, &
      bc_xmin=sll_p_greville, &
      bc_xmax=sll_p_greville, &
      nipts=n1, &
      ncells=ncells1)

   allocate (breaks_eta1(ncells1 + 1))

   call sll_s_new_array_linspace(breaks_eta1, 0.0_wp, 1.0_wp, endpoint=.true.)

   ! Create 1D spline basis along eta1 in [0,1]
   call sll_s_bsplines_new( &
      bsplines=spline_basis_eta1, &
      degree=deg1, &
      periodic=.false., &
      xmin=0.0_wp, &
      xmax=1.0_wp, &
      ncells=ncells1, &
      breaks=breaks_eta1)

   deallocate (breaks_eta1)

   ! Compute number of cells from number of interpolation points along eta2
   call sll_s_spline_1d_compute_num_cells( &
      degree=deg2, &
      bc_xmin=sll_p_periodic, &
      bc_xmax=sll_p_periodic, &
      nipts=n2, &
      ncells=ncells2)

   ! Create 1D spline basis along eta2 in [0,2pi]
   call sll_s_bsplines_new( &
      bsplines=spline_basis_eta2, &
      degree=deg2, &
      periodic=.true., &
      xmin=0.0_wp, &
      xmax=sll_p_twopi, &
      ncells=ncells2)

   ! Compute number of cells from number of interpolation points along eta4
   call sll_s_spline_1d_compute_num_cells( &
      degree=deg4, &
      bc_xmin=sll_p_greville, &
      bc_xmax=sll_p_greville, &
      nipts=n4, &
      ncells=ncells4)

   ! Create 1D spline basis along eta4 in [-1,1]
   call sll_s_bsplines_new( &
      bsplines=spline_basis_eta4, &
      degree=deg4, &
      periodic=.false., &
      xmin=-1.0_wp, &
      xmax=+1.0_wp, &
      ncells=ncells4)

   ! Initialize 1D/2D tensor-product spline interpolators
   call spline_interpolator_1d%init(spline_basis_eta4, &
                                    sll_p_greville, sll_p_greville)

   call spline_interpolator_2d%init(spline_basis_eta1, spline_basis_eta2, &
                                    [sll_p_greville, sll_p_periodic], &
                                    [sll_p_greville, sll_p_periodic])

   ! Get interpolation points and allocate 1D/2D array of values
   call spline_interpolator_1d%get_interp_points(e4_node)
   call spline_interpolator_2d%get_interp_points(e1_node, e2_node)

   npts1 = size(e1_node)
   npts2 = size(e2_node)
   npts4 = size(e4_node)

   write (*, *)
   write (*, *) "Interpolation points along eta1 in [0,1]"
   write (*, *)
   do i1 = 1, npts1
      write (*, *) e1_node(i1)
   end do

   write (*, *)
   write (*, *) "Interpolation points along eta2 in [0,2pi)"
   write (*, *)
   do i2 = 1, npts2
      write (*, *) e2_node(i2)
   end do

   write (*, *)
   write (*, *) "Interpolation points along eta4 in [-1,1]"
   write (*, *)
   do i4 = 1, npts4
      write (*, *) e4_node(i4)
   end do

   write (*, *)

   deallocate (e1_node)
   deallocate (e2_node)
   deallocate (e4_node)

   call spline_interpolator_1d%free()
   call spline_interpolator_2d%free()
   call spline_basis_eta1%free()
   call spline_basis_eta2%free()
   call spline_basis_eta4%free()

end program test_greville_points
