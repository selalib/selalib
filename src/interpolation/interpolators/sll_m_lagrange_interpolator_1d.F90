!> @ingroup interpolators
!> @brief
!> Interpolator class and methods of Lagrange 1D interpolator
!> @details
!> Implements the sll_c_interpolator_1d interface.
module sll_m_lagrange_interpolator_1d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_errors.h"
#include "sll_working_precision.h"

   use sll_m_boundary_condition_descriptors, only: &
      sll_p_one_sided, &
      sll_p_halo, &
      sll_p_periodic

   use sll_m_interpolators_1d_base, only: &
      sll_c_interpolator_1d

   use sll_m_lagrange_interpolation_1d_fast, only: &
      sll_s_lagrange_interpolation_1d_fast_disp_fixed_no_bc, &
      sll_s_lagrange_interpolation_1d_fast_disp_fixed_periodic, &
      sll_s_lagrange_interpolation_1d_fast_disp_fixed_periodicl, &
      sll_s_lagrange_interpolation_1d_fast_disp_fixed_haloc_cells

   use sll_m_lagrange_interpolation_1d, only: &
      sll_s_compute_lagrange_interpolation_1d, &
      sll_s_cubic_spline_1d_eval_array, &
      sll_f_new_lagrange_interpolation_1d, &
      sll_t_lagrange_interpolation_1d

   implicit none

   public :: &
      sll_t_lagrange_interpolator_1d, &
      sll_p_lagrange_centered, &
      sll_p_lagrange_fixed

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   !> Interpolator class of Lagrange 1D interpolator
   type, extends(sll_c_interpolator_1d) :: sll_t_lagrange_interpolator_1d
      !>
      type(sll_t_lagrange_interpolation_1d), pointer :: lagrange
      !>
      sll_int32                                    :: bc_type
      !> Number of points used for interpolation
      sll_int32                                    :: stencil_width
      !> Flag specifying how the Lagrange interpolation points should be chosen (either sll_p_lagrange_centered or sll_p_lagrange_fixed)
      sll_int32                                    :: interval_selection
   contains
      !>
      procedure, pass(interpolator) :: init => initialize_li1d_interpolator
      !>
      procedure :: compute_interpolants => compute_interpolants_li1d
      !>
      procedure :: interpolate_from_interpolant_derivatives_eta1 => interpolate_array_derivatives_li1d
      !>
      procedure :: interpolate_array => interpolate_array_li1d
      !>
      procedure :: interpolate_array_disp => interpolate_array_disp_li1d
      !>
      procedure :: interpolate_array_disp_inplace => interpolate_array_disp_inplace_li1d
      !>
      procedure :: interpolate_from_interpolant_derivative_eta1 => interpolate_derivative_eta1_li1d
      !>
      procedure :: interpolate_from_interpolant_array => interpolate_array_values_li1d
      !>
      procedure :: interpolate_from_interpolant_value => interpolate_value_li1d
      !>
      procedure, pass :: set_coefficients => set_coefficients_li1d
      !>
      procedure, pass :: get_coefficients => get_coefficients_li1d
   end type sll_t_lagrange_interpolator_1d

!PN DEFINED BUT NOT USED
! !> Deallocate the class interpolator
! interface sll_o_delete
!   module procedure delete_li1d
! end interface

   ! Flags for way how to choose the Lagrange points
   sll_int32, parameter :: sll_p_lagrange_centered = 0 !< Flag to specify Lagrange interpolation centered around the interpolation point
   sll_int32, parameter :: sll_p_lagrange_fixed = 1 !< Flag to specify Lagrange interpolation on a fixed interval centered around the point that is displaced (for interpolate_array_disp)

contains  !**********************************************************

   subroutine initialize_li1d_interpolator(interpolator, num_points, xmin, xmax, bc_type, d, periodic_last, interval_selection)
      class(sll_t_lagrange_interpolator_1d), intent(inout) :: interpolator
      sll_int32, intent(in)                        :: d, num_points, bc_type
      sll_real64, intent(in)                       :: xmin, xmax
      sll_int32, intent(in), optional              :: periodic_last
      sll_int32, intent(in), optional              :: interval_selection

      sll_int32                                    :: last

      if (present(periodic_last)) then
         last = periodic_last
      else
         last = 1
      end if

      interpolator%lagrange => sll_f_new_lagrange_interpolation_1d( &
                               num_points, &
                               xmin, &
                               xmax, &
                               bc_type, &
                               d, &
                               last)
      call sll_s_compute_lagrange_interpolation_1D(interpolator%lagrange)

      if (present(interval_selection)) then
         interpolator%interval_selection = interval_selection
      else
         interpolator%interval_selection = sll_p_lagrange_centered
      end if

      select case (interpolator%interval_selection)
      case (sll_p_lagrange_centered)
         interpolator%stencil_width = 2*d
      case (sll_p_lagrange_fixed)
         interpolator%stencil_width = 2*d + 1
      case default
         SLL_ERROR('interpolate_array_disp_li1d', 'Interval selection not implemented.')
      end select

      interpolator%bc_type = bc_type

   end subroutine

   function new_lagrange_interpolator_1d( &
      num_points, &
      xmin, &
      xmax, &
      bc_type, &
      d, &
      periodic_last) result(res)

      type(sll_t_lagrange_interpolator_1d), pointer :: res
      sll_int32, intent(in)               :: num_points
      sll_real64, intent(in)               :: xmin
      sll_real64, intent(in)               :: xmax
      sll_int32, intent(in)               :: bc_type
      sll_int32, intent(in)               :: d
      sll_int32 :: ierr
      sll_int32, intent(in), optional              :: periodic_last
      sll_int32                                    :: last

      SLL_ALLOCATE(res, ierr)

      if (present(periodic_last)) then
         last = periodic_last
      else
         last = 1
      end if

      call initialize_li1d_interpolator( &
         res, &
         num_points, &
         xmin, &
         xmax, &
         bc_type, &
         d, &
         last)

   end function new_lagrange_interpolator_1d

   subroutine interpolate_array_disp_li1d(this, num_pts, data, alpha, output_array)
      class(sll_t_lagrange_interpolator_1d), intent(inout)     :: this
      sll_real64, intent(in) :: alpha
      sll_int32, intent(in)  :: num_pts    ! size of output array
      sll_real64, dimension(:), intent(in) :: data  ! data to be interpolated points where output is desired
      sll_real64, dimension(1:num_pts), intent(out)    :: output_array

      select case (this%interval_selection)
      case (sll_p_lagrange_centered)
         call sll_s_cubic_spline_1d_eval_array(data, -alpha, this%lagrange)
         output_array = this%lagrange%data_out
      case (sll_p_lagrange_fixed)
         select case (this%bc_type)
         case (sll_p_periodic)
            if (this%lagrange%periodic_last == 0) then
     call sll_s_lagrange_interpolation_1d_fast_disp_fixed_periodic(data, output_array, alpha/this%lagrange%deta, this%stencil_width)
            else
    call sll_s_lagrange_interpolation_1d_fast_disp_fixed_periodicl(data, output_array, alpha/this%lagrange%deta, this%stencil_width)
            end if
         case (sll_p_one_sided)
        call sll_s_lagrange_interpolation_1d_fast_disp_fixed_no_bc(data, output_array, alpha/this%lagrange%deta, this%stencil_width)
         case (sll_p_halo)
  call sll_s_lagrange_interpolation_1d_fast_disp_fixed_haloc_cells(data, output_array, alpha/this%lagrange%deta, this%stencil_width)
         case default
            SLL_ERROR('interpolate_array_disp_li1d', 'Boundary type not implemented.')
         end select
      case default
         SLL_ERROR('interpolate_array_disp_li1d', 'Interval selection not implemented.')
      end select

   end subroutine interpolate_array_disp_li1d

   subroutine interpolate_array_disp_inplace_li1d(this, num_pts, data, alpha)
      class(sll_t_lagrange_interpolator_1d), intent(inout)     :: this
      sll_real64, intent(in) :: alpha
      sll_int32, intent(in)  :: num_pts    ! size of output array
      sll_real64, dimension(num_pts), intent(inout) :: data  ! data to be interpolated points where output is desired

      call sll_s_cubic_spline_1d_eval_array(data, -alpha, this%lagrange)
      data = this%lagrange%data_out

   end subroutine interpolate_array_disp_inplace_li1d

!PN DEFINED BUT NOT USED
!subroutine delete_li1d (obj)
!  class(sll_t_lagrange_interpolator_1d) :: obj
!  call delete(obj%lagrange)
!end subroutine delete_li1d

   subroutine interpolate_array_values_li1d( &
      interpolator, &
      num_pts, &
      vals_to_interpolate, &
      output_array)
      class(sll_t_lagrange_interpolator_1d), intent(inout) :: interpolator
      sll_int32, intent(in)                 :: num_pts
      sll_real64, dimension(num_pts), intent(in)   :: vals_to_interpolate
      sll_real64, dimension(num_pts), intent(out)  :: output_array
      !sll_int32 :: ierr
      output_array = 0.0_f64
      print *, 'sll_s_cubic_spline_1d_eval_array:', &
         ' not implemented for lagrange interpolation'
      print *, num_pts
      print *, maxval(vals_to_interpolate)
      output_array = 0._f64
      print *, interpolator%bc_type
      stop
   end subroutine interpolate_array_values_li1d

   subroutine interpolate_array_derivatives_li1d( &
      interpolator, &
      num_pts, &
      vals_to_interpolate, &
      output_array)
      class(sll_t_lagrange_interpolator_1d), intent(inout) :: interpolator
      sll_int32, intent(in)                 :: num_pts
      sll_real64, dimension(:), intent(in)   :: vals_to_interpolate
      sll_real64, dimension(:), intent(out)  :: output_array
      !sll_int32 :: ierr
      output_array = 0.0_f64
      print *, 'interpolate_from_interpolant_derivatives_eta1: ', &
         'not implemented for lagrange interpolation'
      print *, num_pts
      print *, maxval(vals_to_interpolate)
      print *, maxval(output_array)
      print *, interpolator%bc_type
      stop
   end subroutine interpolate_array_derivatives_li1d

   function interpolate_derivative_eta1_li1d(interpolator, eta1) result(val)
      class(sll_t_lagrange_interpolator_1d), intent(in) :: interpolator
      sll_real64             :: val
      sll_real64, intent(in) :: eta1
      print *, 'interpolate_derivative_eta1_li1d: ', &
         'not implemented for lagrange interpolation'
      print *, interpolator%bc_type
      print *, eta1
      val = 0._f64
      stop
   end function

   function interpolate_value_li1d(interpolator, eta1) result(val)
      class(sll_t_lagrange_interpolator_1d), intent(in) :: interpolator
      sll_real64 :: val
      sll_real64, intent(in) :: eta1
      print *, 'interpolate_value_li1d: ', &
         'not implemented for lagrange interpolation'
      val = 0._f64
      print *, eta1
      print *, interpolator%bc_type
      stop
   end function

   subroutine interpolate_array_li1d(this, num_pts, data, coordinates, output_array)
      class(sll_t_lagrange_interpolator_1d), intent(inout)       :: this
      !class(sll_spline_1D),  intent(in)      :: this
      sll_int32, intent(in)                 :: num_pts
      sll_real64, dimension(num_pts), intent(in)   :: coordinates
      sll_real64, dimension(:), intent(in)   :: data
      sll_real64, dimension(num_pts), intent(out)      :: output_array
      ! local variables
      !sll_int32 :: ierr
      ! lagrange interpolation only implemented for constant displacement
      print *, 'interpolate_array_li1d: ', &
         'not implemented for lagrange interpolation'
      print *, maxval(coordinates)
      print *, maxval(data)
      print *, this%bc_type
      output_array = 0._f64
      stop
   end subroutine interpolate_array_li1d

   subroutine compute_interpolants_li1d(interpolator, data_array, &
                                        eta_coords, &
                                        size_eta_coords)
      class(sll_t_lagrange_interpolator_1d), intent(inout) :: interpolator
      sll_real64, dimension(:), intent(in)           :: data_array
      sll_real64, dimension(:), intent(in), optional :: eta_coords
      sll_int32, intent(in), optional :: size_eta_coords

      print *, 'compute_interpolants_li1d:', &
         ' not implemented for lagrange interpolation'

      if (present(eta_coords) .or. present(size_eta_coords)) then
         SLL_ERROR('compute_interpolants_li1d', 'This case is not yet implemented')
      end if

      print *, maxval(data_array)
      print *, interpolator%bc_type
      stop
   end subroutine

   subroutine set_coefficients_li1d(interpolator, coeffs)
      class(sll_t_lagrange_interpolator_1d), intent(inout) :: interpolator
      sll_real64, dimension(:), intent(in), optional :: coeffs
      print *, 'set_coefficients_li1d(): ERROR: This function has not been ', &
         'implemented yet.'
      if (present(coeffs)) then
         print *, '#coeffs present but not used'
      end if
      print *, interpolator%bc_type
      stop
   end subroutine set_coefficients_li1d

   function get_coefficients_li1d(interpolator)
      class(sll_t_lagrange_interpolator_1d), intent(in) :: interpolator
      sll_real64, dimension(:), pointer            :: get_coefficients_li1d

      print *, 'get_coefficients_li1d(): ERROR: This function has not been ', &
         'implemented yet.'
      print *, interpolator%bc_type
      get_coefficients_li1d => null()
      stop
   end function get_coefficients_li1d
   !DEFINE_NULL_INTERP_1D_ARRAY_SUB(sll_t_lagrange_interpolator_1d, interpolate_array_values_li1d)
!DEFINE_NULL_INTERP_1D_ARRAY_SUB(sll_t_lagrange_interpolator_1d, interpolate_array_derivatives_li1d)
!DEFINE_NULL_INTERP_1D_POINTER_SUB(sll_t_lagrange_interpolator_1d, interpolate_pointer_derivatives_li1d)
!DEFINE_NULL_INTERP_ONE_ARG_MSG(sll_t_lagrange_interpolator_1d, interpolate_derivative_eta1_li1d)
!DEFINE_NULL_INTERP_1D_POINTER_SUB(sll_t_lagrange_interpolator_1d, interpolate_pointer_values_li1d)
!DEFINE_NULL_INTERP_ONE_ARG_MSG(sll_t_lagrange_interpolator_1d, interpolate_value_li1d)
!DEFINE_NULL_RECONSTRUCT_1D_ARRAY(sll_t_lagrange_interpolator_1d, reconstruct_array_li1d)
!DEFINE_NULL_INTERP_1D_ARRAY(sll_t_lagrange_interpolator_1d, interpolate_array_li1d)
!DEFINE_NULL_INTERP_1D_ARRAY_MSG(sll_t_lagrange_interpolator_1d, compute_interpolants_li1d)

end module sll_m_lagrange_interpolator_1d
