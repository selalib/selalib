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
 
 subroutine initialize_li1d_interpolator(interpolator)
  class(lagrange_1d_interpolator) :: interpolator

  print*,"pas encore"
 end subroutine initialize_li1d_interpolator

!  subroutine compute_interpolants_li1d(data_array, interpolator)
!   class(lagrange_1d_interpolator)         :: interpolator
!   sll_real64, dimension(:), intent(in)    :: data_array
!   print*,"pas encore"
!  end subroutine compute_interpolants_li1d

DEFINE_NULL_INTERP_1D_ARRAY_MSG(lagrange_1d_interpolator, compute_interpolants_li1d)
DEFINE_NULL_INTERP_1D_ARRAY_SUB(lagrange_1d_interpolator, interpolate_array_derivatives_li1d)
DEFINE_NULL_INTERP_1D_ARRAY(lagrange_1d_interpolator, interpolate_array_li1d)
DEFINE_NULL_INTERPOLATE_1D_DISP(lagrange_1d_interpolator, interpolate_array_disp_li1d)
DEFINE_NULL_INTERP_1D_POINTER_SUB(lagrange_1d_interpolator, interpolate_pointer_derivatives_li1d)
DEFINE_NULL_INTERP_ONE_ARG_MSG(lagrange_1d_interpolator, interpolate_derivative_eta1_li1d)
DEFINE_NULL_INTERP_1D_POINTER_SUB(lagrange_1d_interpolator, interpolate_pointer_values_li1d)
DEFINE_NULL_INTERP_1D_ARRAY_SUB(lagrange_1d_interpolator, interpolate_array_values_li1d)
DEFINE_NULL_INTERP_ONE_ARG_MSG(lagrange_1d_interpolator, interpolate_value_li1d)
DEFINE_NULL_RECONSTRUCT_1D_ARRAY(lagrange_1d_interpolator, reconstruct_array_li1d)

 subroutine delete_li1d (obj)
  class(lagrange_1d_interpolator) :: obj
  call delete(obj%lagrange)
 end subroutine delete_li1d
end module sll_lagrange_interpolator_1d
