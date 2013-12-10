module sll_periodic_interpolator_1d
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#ifndef STDF95
use sll_module_interpolators_1d_base
#endif
use periodic_interp_module
  implicit none

#ifdef STDF95
  type                                    ::  per_1d_interpolator
#else  
  type, extends(sll_interpolator_1d_base) ::  per_1d_interpolator
#endif
    ! Be careful here. For consistency with the other interpolators
    ! num_points is the number of nodes (including both boundaries)
    ! and not the number of cells as used in the periodic interpolator module.
     sll_int32                            :: num_points ! size
     sll_real64                           :: cell_size
     sll_real64                           :: domain_size   ! length of interval
     type(periodic_interp_work), pointer  :: per_interp
#ifdef STDF95
#else
   contains
     procedure, pass(interpolator) :: initialize => initialize_per1d_interpolator
     procedure :: compute_interpolants => compute_interpolants_per1d
     procedure :: interpolate_value => interpolate_value_per1d
     procedure :: interpolate_derivative_eta1 => interpolate_deriv1_per1d
     procedure :: interpolate_array_values => interpolate_values_per1d
     procedure :: interpolate_pointer_values => interpolate_pointer_values_per1d
     procedure :: interpolate_array_derivatives => interpolate_derivatives_per1d
     procedure :: interpolate_pointer_derivatives => &
          interpolate_pointer_derivatives_per1d
     procedure, pass:: interpolate_array => per_interpolate1d
     procedure, pass:: interpolate_array_disp => per_interpolate1d_disp
     procedure, pass:: reconstruct_array
     procedure, pass :: set_coefficients => set_coefficients_per1d
     procedure, pass :: get_coefficients => get_coefficients_per1d
#endif
     
  end type per_1d_interpolator

  interface delete
     module procedure delete_per1d
  end interface delete
  
contains  ! ****************************************************************


  ! the following provides an implementation for the abstract interface 
  ! interpolate1d
  !> Define periodic interpolation of values in data define on original grid at 
  !> points coordinates
  ! Issues with the following function:
  ! - entities referenced through "this" are modified, violating the declared
  !   intent.
  ! - it is probably better to convert this into a subroutine, since data_out
  !   will be allocated on the stack (too big an array will crash the program),
  !   and some copy operation might be involved when "catching" the results.


  function new_periodic_1d_interpolator( &
    num_points, &
    xmin, &
    xmax, &
    type, &
    order)  result(res)
    
    type(per_1d_interpolator),  pointer :: res
    sll_int32,  intent(in)               :: num_points
    sll_real64, intent(in)               :: xmin
    sll_real64, intent(in)               :: xmax
    sll_int32,  intent(in)               :: type
    sll_int32,  intent(in)               :: order
    sll_int32 :: ierr
    SLL_ALLOCATE(res,ierr)
    call initialize_per1d_interpolator( &
         res, &
         num_points, &
         xmin, &
         xmax, &
         type, &
         order)
  end function new_periodic_1d_interpolator





#ifdef STDF95
  function per_interpolate_array(this, num_points, data, coordinates) &
       result(data_out)
    type(per_1d_interpolator),  intent(in)       :: this
#else
  function per_interpolate1d(this, num_points, data, coordinates) &
       result(data_out)
    class(per_1d_interpolator),  intent(in)       :: this
#endif
    !class(sll_spline_1D),  intent(in)      :: this
    sll_int32,  intent(in)                 :: num_points
    sll_real64, dimension(:), intent(in)   :: coordinates
    sll_real64, dimension(:), intent(in)   :: data
    sll_real64, dimension(num_points)      :: data_out
    ! local variables
    sll_int32 :: ierr
    ! periodic interpolation only implemented for constant displacement
    print*, 'periodic_interpolate1d: periodic interpolation not implemented', &
         ' for array of displacements'
    stop
  end function 

#ifdef STDF95
  function per_interpolate_array_at_displacement(this, num_points, &
       data, alpha) &
       result(data_out)
    type(per_1d_interpolator),  intent(in)       :: this
#else
  function per_interpolate1d_disp(this, num_points, data, alpha) &
       result(data_out)
    class(per_1d_interpolator),  intent(in)       :: this
#endif
    sll_int32,  intent(in)                 :: num_points
#ifdef STDF95
    sll_real64                :: alpha
#else
    sll_real64,  intent(in)   :: alpha
#endif
    sll_real64, dimension(:), intent(in)   :: data
    sll_real64, dimension(num_points)      :: data_out
    ! Be careful here. For consistency with the other interpolators
    ! num_points is the number of nodes (including both boundaries)
    ! and not the number of cells as used in the periodic interpolator module.
    call periodic_interp(this%per_interp, data_out, data, &
         alpha/this%cell_size)
    ! complete by periodicity
    data_out(num_points) = data_out(1)
  end function

  ! Both versions F03 and F95 of compute_interpolants_per1d should have the
  ! same name. In the F95 we should add a generic interface around this
  ! subroutine, selecting on the type of interpolator. In the F03 case the
  ! interface is the compute_interpolants routine which gets assigned to
  ! the per1d at initialization time.  
#ifdef STDF95
  subroutine periodic_compute_interpolants( interpolator, data_array,&
       eta_coords, &
       size_eta_coords)
    type(per_1d_interpolator), intent(inout)  :: interpolator
#else
  subroutine compute_interpolants_per1d( interpolator, data_array,&
       eta_coords, &
       size_eta_coords)
       
    class(per_1d_interpolator), intent(inout) :: interpolator
#endif
    sll_real64, dimension(:), intent(in)               :: data_array
    sll_real64, dimension(:), intent(in),optional  :: eta_coords
    sll_int32, intent(in),optional                 :: size_eta_coords
    print*, 'compute_interpolants_per1d:', &
         ' not implemented for periodic interpolation'
    stop
#ifdef STDF95
  end subroutine periodic_compute_interpolants
#else
  end subroutine compute_interpolants_per1d
#endif

  ! Alternative implementation for the function meant to interpolate a
  ! whole array. This implementation fixes some problems in the previous
  ! function. Furthermore, it separates the operation into the more
  ! elementary steps: one is supposed to first compute the interpolants, 
  ! then request to interpolate array values.
  subroutine interpolate_values_per1d( &
    interpolator, &
    num_pts, &
    vals_to_interpolate, &
    output_array )
#ifdef STDF95
    type(per_1d_interpolator),  intent(in) :: interpolator
#else
    class(per_1d_interpolator),  intent(in) :: interpolator
#endif
    sll_int32,  intent(in)                 :: num_pts
    sll_real64, dimension(:), intent(in)   :: vals_to_interpolate
    sll_real64, dimension(:), intent(out)  :: output_array
    sll_int32 :: ierr
    output_array = 0.0
    print*, 'interpolate_values_per1d:', &
         ' not implemented for periodic interpolation'
    stop
  end subroutine interpolate_values_per1d

  subroutine interpolate_pointer_values_per1d( &
    interpolator, &
    num_pts, &
    vals_to_interpolate, &
    output )
#ifdef STDF95
    type(per_1d_interpolator),  intent(in) :: interpolator
#else
    class(per_1d_interpolator),  intent(in) :: interpolator
#endif
    sll_int32,  intent(in)            :: num_pts
    sll_real64, dimension(:), pointer :: vals_to_interpolate
    sll_real64, dimension(:), pointer :: output
    sll_int32 :: ierr
    print*, 'interpolate_pointer_values_per1d: ', &
         'not implemented for periodic interpolation'
    stop
  end subroutine interpolate_pointer_values_per1d


  subroutine interpolate_derivatives_per1d( &
    interpolator, &
    num_pts, &
    vals_to_interpolate, &
    output_array )


#ifdef STDF95
    type(per_1d_interpolator),  intent(in) :: interpolator
#else
    class(per_1d_interpolator),  intent(in) :: interpolator
#endif
    sll_int32,  intent(in)                 :: num_pts
    sll_real64, dimension(:), intent(in)   :: vals_to_interpolate
    sll_real64, dimension(:), intent(out)  :: output_array
    sll_int32 :: ierr
    output_array = 0.0
     print*, 'interpolate_array_derivatives: ', &
         'not implemented for periodic interpolation'
    stop
  end subroutine interpolate_derivatives_per1d

  subroutine interpolate_pointer_derivatives_per1d( &
    interpolator, &
    num_pts, &
    vals_to_interpolate, &
    output )
#ifdef STDF95
    type(per_1d_interpolator),  intent(in) :: interpolator
#else
    class(per_1d_interpolator),  intent(in) :: interpolator
#endif
    sll_int32,  intent(in)              :: num_pts
    sll_real64, dimension(:), pointer   :: vals_to_interpolate
    sll_real64, dimension(:), pointer   :: output
    sll_int32 :: ierr
     print*, 'interpolate_pointer_derivatives_per1d:  ', &
         'not implemented for periodic interpolation'
    stop
  end subroutine interpolate_pointer_derivatives_per1d

#ifdef STDF95
  function periodic_interpolate_value( interpolator, eta1 ) result(val)
    type(per_1d_interpolator), intent(inout) :: interpolator
#else
  function interpolate_value_per1d( interpolator, eta1 ) result(val)
    class(per_1d_interpolator), intent(inout) :: interpolator
#endif
    sll_real64 :: val
    sll_real64, intent(in) :: eta1
     print*, 'interpolate_value_per1d: ', &
         'not implemented for periodic interpolation'
    stop
  end function
  
#ifdef STDF95
  function periodic_interpolate_derivative_eta1( interpolator, eta1 ) &
       result(val)
    type(per_1d_interpolator), intent(inout)  :: interpolator
#else
  function interpolate_deriv1_per1d( interpolator, eta1 ) result(val)
    class(per_1d_interpolator), intent(inout) :: interpolator
#endif
    sll_real64             :: val
    sll_real64, intent(in) :: eta1
     print*, 'interpolate_deriv1_per1d: ', &
         'not implemented for periodic interpolation'
    stop
#ifdef STDF95
  end function periodic_interpolate_derivative_eta1
#else
  end function interpolate_deriv1_per1d
#endif

#ifdef STDF95
  function periodic_interpolate_derivative_f95( interpolator, eta1 ) result(val)
    type(per_1d_interpolator), intent(in) :: interpolator
#else
  function interpolate_derivative_f95( interpolator, eta1 ) result(val)
    class(per_1d_interpolator), intent(in) :: interpolator
#endif
    sll_real64 :: val
    sll_real64, intent(in) :: eta1
     print*, 'interpolate_derivative_f95: ', &
         'not implemented for periodic interpolation'
    stop
#ifdef STDF95
  end function periodic_interpolate_derivative_f95
#else
  end function interpolate_derivative_f95
#endif


  ! Why is the name of this function changing depending on the standard?
  ! only one will be compiled anyway!!

  !> initialize periodic interpolator
#ifdef STDF95
  subroutine periodic_initialize( &
#else
  subroutine initialize_per1d_interpolator( &
#endif
    interpolator, &
    num_points, &
    xmin, &
    xmax, &
    type, &
    order)

#ifdef STDF95
    type(per_1d_interpolator),  intent(inout)  :: interpolator 
#else
    class(per_1d_interpolator),  intent(inout) :: interpolator 
#endif
    sll_int32,  intent(in)               :: num_points
    sll_real64, intent(in)               :: xmin
    sll_real64, intent(in)               :: xmax
    sll_int32,  intent(in)               :: type
    sll_int32,  intent(in)               :: order
    sll_int32                            :: ierr
    sll_int32  :: i  
    sll_real64 :: delta
    
    ! Be careful here. For consistency with the other interpolators
    ! num_points is the number of nodes (including both boundaries)
    ! and not the number of cells as used in the periodic interpolator module.
    interpolator%num_points = num_points - 1
    interpolator%cell_size  = (xmax-xmin) / (num_points-1)
    interpolator%domain_size = xmax-xmin

    call initialize_periodic_interp(interpolator%per_interp, num_points-1, &
         type, order)
  end subroutine

  function reconstruct_array(this, num_points, data) result(res)
    ! dummy procedure
#ifdef STDF95
    type(per_1d_interpolator), intent(in)      :: this
#else
    class(per_1d_interpolator), intent(in)     :: this
#endif
       sll_int32, intent(in)                :: num_points! size of output array
       sll_real64, dimension(:), intent(in) :: data   ! data to be interpolated 
       sll_real64, dimension(num_points)    :: res
       res(:) = 0.0_f64
  end function reconstruct_array

  subroutine delete_per1d( obj )
#ifdef STDF95
    type(per_1d_interpolator) :: obj
#else  
    class(per_1d_interpolator) :: obj
#endif
    call delete(obj%per_interp)
  end subroutine delete_per1d

  subroutine set_coefficients_per1d( interpolator, coeffs )
#ifdef STDF95
    type(per_1d_interpolator), intent(inout)  :: interpolator
#else  
    class(per_1d_interpolator), intent(inout) :: interpolator
#endif
    sll_real64, dimension(:), intent(in), optional :: coeffs
    print *, 'set_coefficients_per1d(): ERROR: This function has not been ', &
         'implemented yet.'
    stop
  end subroutine set_coefficients_per1d


  function get_coefficients_per1d(interpolator)
#ifdef STDF95
    type(per_1d_interpolator), intent(in)  :: interpolator
#else  
    class(per_1d_interpolator), intent(in) :: interpolator
#endif
    sll_real64, dimension(:), pointer            :: get_coefficients_per1d     
    
    print *, 'get_coefficients_per1d(): ERROR: This function has not been ', &
         'implemented yet.' 
  end function get_coefficients_per1d

end module sll_periodic_interpolator_1d
