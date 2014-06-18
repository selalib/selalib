module sll_lagrange_interpolator_1d
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
!#include "sll_interpolators_1d_base_macros.h"
#ifndef STDF95
use sll_module_interpolators_1d_base
#endif
 use sll_lagrange_interpolation
implicit none

#ifdef STDF95
 type :: lagrange_1d_interpolator
   type(sll_lagrange_interpolation_1D), pointer :: lagrange
   sll_int32                                    :: bc_type 
#else
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
   procedure, pass :: set_coefficients => set_coefficients_li1d
   procedure, pass :: get_coefficients => get_coefficients_li1d
#endif
 end type lagrange_1d_interpolator

 interface delete
   module procedure delete_li1d
 end interface

contains  !**********************************************************
 
#ifdef STDF95
subroutine lagrange_interpolation_1d_initialize_interpolator(interpolator,num_points,xmin,xmax,bc_type,d)
  type(lagrange_1d_interpolator), intent(inout) :: interpolator
#else
subroutine initialize_li1d_interpolator(interpolator,num_points,xmin,xmax,bc_type,d)
  class(lagrange_1d_interpolator), intent(inout) :: interpolator
#endif
    sll_int32, intent(in)                        :: d,num_points,bc_type
    sll_real64, intent(in)                       :: xmin,xmax

    interpolator%lagrange => new_lagrange_interpolation_1D( &
           num_points, &
           xmin, &
           xmax, &
           bc_type, &
           d)
end subroutine 

#ifdef STDF95
function lagrange_interpolation_1d_interpolate_array_disp(this, num_points, data, alpha) result(data_out)
  type(lagrange_1d_interpolator), intent(in)     :: this
#else
function interpolate_array_disp_li1d(this, num_points, data, alpha) result(data_out)
  class(lagrange_1d_interpolator), intent(in)     :: this
#endif
  sll_real64, intent(in) :: alpha
  sll_int32, intent(in)  :: num_points    ! size of output array
  sll_real64, dimension(:), intent(in) :: data  ! data to be interpolated points where output is desired
  sll_real64, dimension(1:num_points)    :: data_out
call compute_lagrange_interpolation_1D(alpha,this%lagrange)
call interpolate_array_values(data,this%lagrange)
data_out=this%lagrange%data_out

end function 

subroutine delete_li1d (obj)
#ifdef STDF95 
  type(lagrange_1d_interpolator) :: obj
#else
  class(lagrange_1d_interpolator) :: obj
#endif
  call delete(obj%lagrange)
end subroutine delete_li1d


subroutine interpolate_array_values_li1d( &
    interpolator, &
    num_pts, &
    vals_to_interpolate, &
    output_array )
#ifdef STDF95
    type(lagrange_1d_interpolator),  intent(in) :: interpolator
#else
    class(lagrange_1d_interpolator),  intent(in) :: interpolator
#endif
    sll_int32,  intent(in)                 :: num_pts
    sll_real64, dimension(:), intent(in)   :: vals_to_interpolate
    sll_real64, dimension(:), intent(out)  :: output_array
    !sll_int32 :: ierr
    output_array = 0.0
    print*, 'interpolate_array_values:', &
         ' not implemented for lagrange interpolation'
    print *,num_pts
    print *,maxval(vals_to_interpolate)
    output_array = 0._f64
    print *,interpolator%bc_type
    stop
end subroutine interpolate_array_values_li1d


subroutine interpolate_array_derivatives_li1d( &
    interpolator, &
    num_pts, &
    vals_to_interpolate, &
    output_array )
#ifdef STDF95
    type(lagrange_1d_interpolator),  intent(in) :: interpolator
#else
    class(lagrange_1d_interpolator),  intent(in) :: interpolator
#endif
    sll_int32,  intent(in)                 :: num_pts
    sll_real64, dimension(:), intent(in)   :: vals_to_interpolate
    sll_real64, dimension(:), intent(out)  :: output_array
    !sll_int32 :: ierr
    output_array = 0.0
    print*, 'interpolate_array_derivatives: ', &
         'not implemented for lagrange interpolation'
    print *,num_pts
    print *,maxval(vals_to_interpolate)
    print *,maxval(output_array)
    print *,interpolator%bc_type
    stop
end subroutine interpolate_array_derivatives_li1d

subroutine interpolate_pointer_derivatives_li1d( &
    interpolator, &
    num_pts, &
    vals_to_interpolate, &
    output )
#ifdef STDF95
    type(lagrange_1d_interpolator),  intent(in) :: interpolator
#else
    class(lagrange_1d_interpolator),  intent(in) :: interpolator
#endif
    sll_int32,  intent(in)              :: num_pts
    sll_real64, dimension(:), pointer   :: vals_to_interpolate
    sll_real64, dimension(:), pointer   :: output
    !sll_int32 :: ierr
    print*, 'interpolate_pointer_derivatives_li1d:  ', &
         'not implemented for lagrange interpolation'
    print *,interpolator%bc_type
    print *,num_pts
    print *,maxval(vals_to_interpolate)
    print *,maxval(output)    
    stop
end subroutine interpolate_pointer_derivatives_li1d

#ifdef STDF95
  function lagrange_interpolation_1d_interpolate_derivative_eta1(interpolator, eta1) result(val)
    type(lagrange_1d_interpolator), intent(inout)  :: interpolator
#else
  function interpolate_derivative_eta1_li1d( interpolator, eta1 ) result(val)
    class(lagrange_1d_interpolator), intent(inout) :: interpolator
#endif
    sll_real64             :: val
    sll_real64, intent(in) :: eta1
     print*, 'interpolate_derivative_eta1_li1d: ', &
         'not implemented for lagrange interpolation'
    print *,interpolator%bc_type
    print *,eta1
    val = 0._f64     
    stop
  end function

subroutine interpolate_pointer_values_li1d( &
    interpolator, &
    num_pts, &
    vals_to_interpolate, &
    output )
#ifdef STDF95
    type(lagrange_1d_interpolator),  intent(in) :: interpolator
#else
    class(lagrange_1d_interpolator),  intent(in) :: interpolator
#endif
    sll_int32,  intent(in)            :: num_pts
    sll_real64, dimension(:), pointer :: vals_to_interpolate
    sll_real64, dimension(:), pointer :: output
    !sll_int32 :: ierr
    print*, 'interpolate_pointer_values_li1d: ', &
         'not implemented for lagrange interpolation'
    print *,num_pts
    print *,maxval(vals_to_interpolate)
    print *,maxval(output)
    print *,interpolator%bc_type
    stop
end subroutine interpolate_pointer_values_li1d

#ifdef STDF95
  function lagrange_interpolation_1d_interpolate_value(interpolator, eta1 ) result(val)
    type(lagrange_1d_interpolator), intent(inout) :: interpolator
#else
  function interpolate_value_li1d( interpolator, eta1 ) result(val)
    class(lagrange_1d_interpolator), intent(inout) :: interpolator
#endif
    sll_real64 :: val
    sll_real64, intent(in) :: eta1
     print*, 'interpolate_value_li1d: ', &
         'not implemented for lagrange interpolation'
    val = 0._f64
    print *,eta1
    print *,interpolator%bc_type     
    stop
  end function

  function reconstruct_array_li1d(this, num_points, data) result(res)
    ! dummy procedure
#ifdef STDF95
    type(lagrange_1d_interpolator), intent(in)      :: this
#else
    class(lagrange_1d_interpolator), intent(in)     :: this
#endif
       sll_int32, intent(in)                :: num_points! size of output array
       sll_real64, dimension(:), intent(in) :: data   ! data to be interpolated 
       sll_real64, dimension(num_points)    :: res
       print *,'#warning reconstruct_array_li1d dummy function'
       print *,num_points
       print *,maxval(data)
       print *,this%bc_type
       res(:) = 0.0_f64
  end function reconstruct_array_li1d

#ifdef STDF95
  function lagrange_interpolation_1d_interpolate_array(this, num_points, data, coordinates) &
       result(data_out)
    type(lagrange_1d_interpolator),  intent(in)       :: this
#else
  function interpolate_array_li1d(this, num_points, data, coordinates) &
       result(data_out)
    class(lagrange_1d_interpolator),  intent(in)       :: this
#endif
    !class(sll_spline_1D),  intent(in)      :: this
    sll_int32,  intent(in)                 :: num_points
    sll_real64, dimension(:), intent(in)   :: coordinates
    sll_real64, dimension(:), intent(in)   :: data
    sll_real64, dimension(num_points)      :: data_out
    ! local variables
    !sll_int32 :: ierr
    ! lagrange interpolation only implemented for constant displacement
    print*, 'interpolate_array_li1d: ', &
         'not implemented for lagrange interpolation'
    print *,maxval(coordinates)
    print *,maxval(data)
    print *,this%bc_type
    data_out = 0._f64
    stop
  end function 

#ifdef STDF95
  subroutine lagrange_interpolation_1d_compute_interpolants( interpolator, data_array,&
       eta_coords, &
       size_eta_coords)
       
    type(lagrange_1d_interpolator), intent(inout)  :: interpolator
#else
    subroutine compute_interpolants_li1d( interpolator, data_array,&
         eta_coords, &
         size_eta_coords)
      class(lagrange_1d_interpolator), intent(inout) :: interpolator
#endif
    sll_real64, dimension(:), intent(in)               :: data_array
    sll_real64, dimension(:), intent(in),optional  :: eta_coords
    sll_int32, intent(in),optional                 :: size_eta_coords
    print*, 'compute_interpolants_li1d:', &
         ' not implemented for lagrange interpolation'
    if(present(eta_coords))then
      print *,'eta_coords present but not used'
    endif     
    if(present(size_eta_coords))then
      print *,'size_eta_coords present but not used'
    endif
    print *,maxval(data_array)
    print *,interpolator%bc_type     
    stop
  end subroutine

  subroutine set_coefficients_li1d( interpolator, coeffs )
#ifdef STDF95
    type(lagrange_1d_interpolator),  intent(inout)  :: interpolator
#else
    class(lagrange_1d_interpolator),  intent(inout) :: interpolator
#endif
    sll_real64, dimension(:), intent(in), optional :: coeffs
    print *, 'set_coefficients_li1d(): ERROR: This function has not been ', &
         'implemented yet.'
    if(present(coeffs))then
      print *,'#coeffs present but not used'
    endif
    print *,interpolator%bc_type     
    stop
  end subroutine set_coefficients_li1d
  
  
  function get_coefficients_li1d(interpolator)
#ifdef STDF95
    type(lagrange_1d_interpolator),  intent(in) :: interpolator
#else
    class(lagrange_1d_interpolator),  intent(in) :: interpolator
#endif
    sll_real64, dimension(:), pointer            :: get_coefficients_li1d     
    
    print *, 'get_coefficients_li1d(): ERROR: This function has not been ', &
         'implemented yet.'
    print *,interpolator%bc_type     
    get_coefficients_li1d => null()
    stop      
  end function get_coefficients_li1d
  !DEFINE_NULL_INTERP_1D_ARRAY_SUB(lagrange_1d_interpolator, interpolate_array_values_li1d)
!DEFINE_NULL_INTERP_1D_ARRAY_SUB(lagrange_1d_interpolator, interpolate_array_derivatives_li1d)
!DEFINE_NULL_INTERP_1D_POINTER_SUB(lagrange_1d_interpolator, interpolate_pointer_derivatives_li1d)
!DEFINE_NULL_INTERP_ONE_ARG_MSG(lagrange_1d_interpolator, interpolate_derivative_eta1_li1d)
!DEFINE_NULL_INTERP_1D_POINTER_SUB(lagrange_1d_interpolator, interpolate_pointer_values_li1d)
!DEFINE_NULL_INTERP_ONE_ARG_MSG(lagrange_1d_interpolator, interpolate_value_li1d)
!DEFINE_NULL_RECONSTRUCT_1D_ARRAY(lagrange_1d_interpolator, reconstruct_array_li1d)
!DEFINE_NULL_INTERP_1D_ARRAY(lagrange_1d_interpolator, interpolate_array_li1d)
!DEFINE_NULL_INTERP_1D_ARRAY_MSG(lagrange_1d_interpolator, compute_interpolants_li1d)

end module sll_lagrange_interpolator_1d
