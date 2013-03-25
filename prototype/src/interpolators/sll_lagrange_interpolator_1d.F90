module sll_lagrange_interpolator_1d
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_interpolators_1d_base_macros.h"
#ifndef STDF95
use sll_module_interpolators_1d_base
#endif
 use sll_lagrange_interpolation
implicit none

 type,extends(sll_interpolator_1d_base) :: lagrange_1d_interpolator
   type(sll_lagrange_interpolation_1D), pointer :: lagrange
   sll_int32                                    :: bc_type 
   contains
   procedure,pass(interpolator) :: initialize => initialize_li1d_interpolator
   procedure :: compute_interpolants => compute_interpolants_li1d
   procedure :: interpolate_array_derivatives => interpolate_array_derivatives_li1d
   procedure :: interpolate_array => interpolate_array_li1d
   procedure :: interpolate_array_disp => interpolate_array_disp_li1d
   procedure :: interpolate_pointer_derivatives => interpolate_pointer_derivatives_li1d
   procedure :: interpolate_derivative_eta1 => interpolate_derivative_eta1_li1d
   procedure :: interpolate_pointer_values => interpolate_pointer_values_li1d
   procedure :: interpolate_array_values => interpolate_array_values_li1d
   procedure :: interpolate_value => interpolate_value_li1d
   procedure :: reconstruct_array => reconstruct_array_li1d
 end type lagrange_1d_interpolator

 interface delete
   module procedure delete_li1d
 end interface
contains  !**********************************************************
 
subroutine initialize_li1d_interpolator(interpolator,num_points,xmin,xmax,bc_type,d)
  class(lagrange_1d_interpolator), intent(inout) :: interpolator
    sll_int32, intent(in)                        :: d,num_points,bc_type
    sll_real64, intent(in)                       :: xmin,xmax

    interpolator%lagrange => new_lagrange_interpolation_1D( &
           num_points, &
           xmin, &
           xmax, &
           bc_type, &
           d)
end subroutine initialize_li1d_interpolator

function interpolate_array_disp_li1d(this, num_points, data, alpha) result(data_out)
  class(lagrange_1d_interpolator), intent(in)     :: this
  sll_real64, intent(in) :: alpha
  sll_int32, intent(in)  :: num_points    ! size of output array
  sll_real64, dimension(:), intent(in) :: data  ! data to be interpolated points where output is desired
  sll_real64, dimension(1:num_points)    :: data_out

call compute_lagrange_interpolation_1D(alpha,this%lagrange)
call interpolate_array_values(data,this%lagrange)
data_out=this%lagrange%data_out

end function interpolate_array_disp_li1d

subroutine delete_li1d (obj)
  class(lagrange_1d_interpolator) :: obj
  call delete(obj%lagrange)
end subroutine delete_li1d

DEFINE_NULL_INTERP_1D_ARRAY_SUB(lagrange_1d_interpolator, interpolate_array_values_li1d)
DEFINE_NULL_INTERP_1D_ARRAY_SUB(lagrange_1d_interpolator, interpolate_array_derivatives_li1d)
DEFINE_NULL_INTERP_1D_POINTER_SUB(lagrange_1d_interpolator, interpolate_pointer_derivatives_li1d)
DEFINE_NULL_INTERP_ONE_ARG_MSG(lagrange_1d_interpolator, interpolate_derivative_eta1_li1d)
DEFINE_NULL_INTERP_1D_POINTER_SUB(lagrange_1d_interpolator, interpolate_pointer_values_li1d)
DEFINE_NULL_INTERP_ONE_ARG_MSG(lagrange_1d_interpolator, interpolate_value_li1d)
DEFINE_NULL_RECONSTRUCT_1D_ARRAY(lagrange_1d_interpolator, reconstruct_array_li1d)
DEFINE_NULL_INTERP_1D_ARRAY(lagrange_1d_interpolator, interpolate_array_li1d)
DEFINE_NULL_INTERP_1D_ARRAY_MSG(lagrange_1d_interpolator, compute_interpolants_li1d)
!DEFINE_NULL_INTERPOLATE_1D_DISP(lagrange_1d_interpolator, interpolate_array_disp_li1d)

! function interpolate_array_li1d(this, num_points, data, coordinates) result(data_out)
!   class(lagrange_1d_interpolator), intent(in)     :: this
!   sll_int32, intent(in)  :: num_points    ! size of output array
!   sll_real64, dimension(:), intent(in) :: data  ! data to be interpolated points where output is desired
!   sll_real64, dimension(:), intent(in) :: coordinates  
!   sll_real64, dimension(1:num_points)    :: data_out
! 
!   call compute_lagrange_interpolation_1D(coordinates,this%lagrange)
!   call interpolate_array_values(data,this%lagrange)
!   data_out=this%lagrange%data_out
! end function interpolate_array_li1d

! subroutine compute_interpolants_li1d(interpolator, data_array)
!   class(lagrange_1d_interpolator), intent(inout) :: interpolator
!   sll_real64, dimension(:), intent(in)           :: data_array
! 
!   call compute_lagrange_interpolation_1D(data_array,interpolator%lagrange)
! end subroutine compute_interpolants_li1di


end module sll_lagrange_interpolator_1d
