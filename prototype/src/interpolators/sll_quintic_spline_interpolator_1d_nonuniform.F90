module sll_quintic_spline_interpolator_1d_nonuniform
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#ifndef STDF95
use sll_module_interpolators_1d_base
#endif
use sll_quintic_splines
  implicit none

#ifdef STDF95
  type                                    :: quintic_spline_1d_interpolator_nonuniform
#else  
  type, extends(sll_interpolator_1d_base) :: quintic_spline_1d_interpolator_nonuniform
#endif
     sll_real64, dimension(:), pointer            :: interpolation_points 
     sll_int32                                    :: num_points ! size
     sll_int32                                    :: bc_type
     type(quintic_splines_nonuniform_plan), pointer  :: spline
#ifdef STDF95
#else
   contains
     procedure, pass(interpolator) :: initialize => initialize_qs1d_interpolator_nonuniform
     procedure :: compute_interpolants => compute_interpolants_qs1d_nonuniform
     procedure :: interpolate_value => interpolate_value_qs1d_nonuniform
     procedure :: interpolate_derivative_eta1 => interpolate_value_qs1d_nonuniform ! Is not used
     procedure :: interpolate_array_values => interpolate_values_qs1d_nonuniform
     procedure :: interpolate_pointer_values => interpolate_pointer_values_qs1d_nonuniform
     procedure :: interpolate_array_derivatives => interpolate_values_qs1d_nonuniform ! Is not used
     procedure :: interpolate_pointer_derivatives => &
          interpolate_pointer_values_qs1d_nonuniform ! Is not used
     procedure, pass:: interpolate_array => spline_interpolate1d_nonuniform
     procedure, pass:: interpolate_array_disp => spline_interpolate1d_disp_nonuniform
     procedure, pass:: reconstruct_array
     procedure, pass :: set_coefficients => set_coefficients_qs1d_nonuniform
     procedure, pass :: get_coefficients => get_coefficients_qs1d_nonuniform

     !generic :: initialize => initialize_qs1d_interpolato
#endif
  end type quintic_spline_1d_interpolator_nonuniform

  interface delete_nonuniform
     module procedure delete_qs1d_nonuniform
  end interface delete_nonuniform
  
contains  ! ****************************************************************


  ! the following provides an implementation for the abstract interface 
  !interpolate1d
  !> Define spline interpolation of values in data define on original grid at 
  !> points coordinates
  ! Issues with the following function:
  ! - entities referenced through "this" are modified, violating the declared
  !   intent.
  ! - it is probably better to convert this into a subroutine, since data_out
  !   will be allocated on the stack (too big an array will crash the program),
  !   and some copy operation might be involved when "catching" the results.

#ifdef STDF95
  function quintic_spline_interpolate_array_nonuniform(this, num_points, data, coordinates) &
       result(data_out)
    type(quintic_spline_1d_interpolator_nonuniform),  intent(in)       :: this
#else
  function spline_interpolate1d_nonuniform(this, num_points, data, coordinates) &
       result(data_out)
    class(quintic_spline_1d_interpolator_nonuniform),  intent(in)       :: this
#endif

    !class(sll_spline_1D),  intent(in)      :: this
    sll_int32,  intent(in)                 :: num_points
    sll_real64, dimension(:), intent(in)   :: coordinates
    sll_real64, dimension(:), intent(in)   :: data
    sll_real64, dimension(num_points)      :: data_out
    ! compute the interpolating spline coefficients
    call compute_quintic_coeffs_nonuniform( data, this%spline )
    data_out =  quintic_splines_interpolator_nonuniform_array( &
                       coordinates, num_points, this%spline )
  end function 

#ifdef STDF95
  function quintic_spline_interpolate_array_at_displacement_nonuniform(this, num_points, &
       data, coordinates) &
       result(data_out)
    type(quintic_spline_1d_interpolator_nonuniform),  intent(in)       :: this
#else
  function spline_interpolate1d_disp_nonuniform(this, num_points, data, alpha) &
       result(data_out)
    class(quintic_spline_1d_interpolator_nonuniform),  intent(in)       :: this
#endif

    !class(sll_spline_1D),  intent(in)      :: this
    sll_int32,  intent(in)                 :: num_points
#ifdef STDF95
    sll_real64                :: alpha
#else
    sll_real64,  intent(in)   :: alpha
#endif
    sll_real64, dimension(:), intent(in)   :: data
    sll_real64, dimension(num_points)      :: data_out
    sll_real64, dimension(num_points)      :: coordinates
    sll_real64 :: length, delta
    sll_real64 :: xmin, xmax 
    sll_int32 :: i
    ! compute_quintic the interpolating spline coefficients
    call compute_quintic_coeffs_nonuniform( data, this%spline )
    ! compute array of coordinates where interpolation is performed from displacement
    length = this%interpolation_points(num_points) - &
             this%interpolation_points(1)
    delta = this%interpolation_points(2) - this%interpolation_points(1)
    xmin = this%interpolation_points(1)
    xmax = this%interpolation_points(num_points)

    if (alpha > 0 ) then 
       do i = 1, num_points
          coordinates(i) = max(this%interpolation_points(i) - alpha, xmin)
          SLL_ASSERT((xmin <=coordinates(i)).and.(coordinates(i) <= xmax))
       end do
    else
       do i = 1, num_points
          coordinates(i) = min(this%interpolation_points(i) - alpha, xmax)
          SLL_ASSERT((xmin <=coordinates(i)).and.(coordinates(i) <= xmax))
       end do
    endif

    data_out = quintic_splines_interpolator_nonuniform_array( coordinates, &
                                                num_points, this%spline )
  end function

  ! Both versions F03 and F95 of compute_interpolants_qs1d should have the
  ! same name. In the F95 we should add a generic interface around this
  ! subroutine, selecting on the type of interpolator. In the F03 case the
  ! interface is the compute_quintic_interpolants routine which gets assigned to
  ! the qs1d at initialization time.  
#ifdef STDF95
  subroutine quintic_spline_compute_interpolants_nonuniform( interpolator, data_array,&
       eta_coords, &
       size_eta_coords)

    type(quintic_spline_1d_interpolator_nonuniform), intent(inout)  :: interpolator
#else
    subroutine compute_interpolants_qs1d_nonuniform( interpolator, data_array,&
         eta_coords, &
         size_eta_coords)
    class(quintic_spline_1d_interpolator_nonuniform), intent(inout) :: interpolator
#endif
    sll_real64, dimension(:), intent(in)               :: data_array
    sll_real64, dimension(:), intent(in),optional  :: eta_coords
    sll_int32, intent(in),optional                 :: size_eta_coords
    
    if(present(eta_coords))then
      !print *,'#Warning eta_coords present but not used'
    endif
    if(present(size_eta_coords))then
      !print *,'#Warning size_eta_coords present but not used'
    endif
    
    
    call compute_quintic_coeffs_nonuniform( data_array, interpolator%spline )
#ifdef STDF95
  end subroutine quintic_spline_compute_interpolants_nonuniform
#else
  end subroutine compute_interpolants_qs1d_nonuniform
#endif

  ! Alternative implementation for the function meant to interpolate a
  ! whole array. This implementation fixes some problems in the previous
  ! function. Furthermore, it separates the operation into the more
  ! elementary steps: one is supposed to first compute the interpolants, 
  ! then request to interpolate array values.
  subroutine interpolate_values_qs1d_nonuniform( &
    interpolator, &
    num_pts, &
    vals_to_interpolate, &
    output_array )
#ifdef STDF95
    type(quintic_spline_1d_interpolator_nonuniform),  intent(in) :: interpolator
#else
    class(quintic_spline_1d_interpolator_nonuniform),  intent(in) :: interpolator
#endif
    sll_int32,  intent(in)                 :: num_pts
    sll_real64, dimension(:), intent(in)   :: vals_to_interpolate
    sll_real64, dimension(:), intent(out)  :: output_array
    output_array = quintic_splines_interpolator_nonuniform_array( &
              vals_to_interpolate, num_pts, interpolator%spline)
  end subroutine interpolate_values_qs1d_nonuniform

  subroutine interpolate_pointer_values_qs1d_nonuniform( &
    interpolator, &
    num_pts, &
    vals_to_interpolate, &
    output )
#ifdef STDF95
    type(quintic_spline_1d_interpolator_nonuniform),  intent(in) :: interpolator
#else
    class(quintic_spline_1d_interpolator_nonuniform),  intent(in) :: interpolator
#endif
    sll_int32,  intent(in)            :: num_pts
    sll_real64, dimension(:), pointer :: vals_to_interpolate
    sll_real64, dimension(:), pointer :: output
    output => quintic_splines_interpolator_nonuniform_pointer(&
          vals_to_interpolate, num_pts, interpolator%spline)
  end subroutine interpolate_pointer_values_qs1d_nonuniform

#ifdef STDF95
  function quintic_spline_non_uniform_1d_interpolate_value( interpolator, eta1 ) result(val)
    type(quintic_spline_1d_interpolator_nonuniform), intent(inout) :: interpolator
#else
  function interpolate_value_qs1d_nonuniform( interpolator, eta1 ) result(val)
    class(quintic_spline_1d_interpolator_nonuniform), intent(inout) :: interpolator
#endif

    sll_real64 :: val
    sll_real64, intent(in) :: eta1
    val = quintic_splines_interpolator_nonuniform_value( eta1, interpolator%spline )

  end function

  ! Why is the name of this function changing depending on the standard?
  ! only one will be compiled anyway!!

  !> initialize quintic spline interpolator
#ifdef STDF95
  subroutine quintic_spline_initialize_nonuniform( &
#else
  subroutine initialize_qs1d_interpolator_nonuniform( &
#endif
    knots, &
    interpolator, &
    num_points, &
    xmin, &
    xmax, &
    bc_type, &
    slope_left, &
    slope_right )

#ifdef STDF95
    type(quintic_spline_1d_interpolator_nonuniform),  intent(inout)  :: interpolator 
#else
    class(quintic_spline_1d_interpolator_nonuniform),  intent(inout) :: interpolator 
#endif
    sll_int32,  intent(in)           :: num_points
    sll_real64, intent(in)           :: xmin
    sll_real64, intent(in)           :: xmax
    sll_int32,  intent(in)           :: bc_type
    sll_real64, intent(in), optional :: slope_left
    sll_real64, intent(in), optional :: slope_right
    sll_int32                        :: ierr
    sll_int32  :: i  
    sll_real64, intent(in), dimension(:) :: knots
    print *,'#Warning bc_type provided but not used',bc_type
    print *,'#Warning xmax provided but not used',xmax
    if(present(slope_left))then
      print *,'#Warning slope_left present but not used'
    endif
    if(present(slope_right))then
      print *,'#Warning slope_right present but not used'
    endif
    interpolator%num_points = num_points
    SLL_ALLOCATE(interpolator%interpolation_points(num_points),ierr)
    interpolator%interpolation_points(1) = xmin
    do i = 2, num_points
       interpolator%interpolation_points(i) = knots(i)
    end do
    interpolator%spline => new_quintic_splines_nonuniform(knots)
  end subroutine

  function reconstruct_array(this, num_points, data) result(res)
    ! dummy procedure
#ifdef STDF95
    type(quintic_spline_1d_interpolator_nonuniform), intent(in)      :: this
#else
    class(quintic_spline_1d_interpolator_nonuniform), intent(in)     :: this
#endif
       sll_int32, intent(in)                :: num_points! size of output array
       sll_real64, dimension(:), intent(in) :: data   ! data to be interpolated 
       sll_real64, dimension(num_points)    :: res
       res(:) = 0.0_f64       
       print *,'#Warning dummy procedure reconstruct_array'
       print *,maxval(data)
       print *,this%num_points
       
  end function reconstruct_array
  
  subroutine delete_qs1d_nonuniform( obj )
#ifdef STDF95
    type(quintic_spline_1d_interpolator_nonuniform) :: obj
#else
    class(quintic_spline_1d_interpolator_nonuniform) :: obj
#endif
    call delete_quintic_splines_nonuniform(obj%spline)
  end subroutine delete_qs1d_nonuniform


  subroutine set_coefficients_qs1d_nonuniform( interpolator, coeffs )
#ifdef STDF95
    type(quintic_spline_1d_interpolator_nonuniform),intent(inout)  :: interpolator
#else
    class(quintic_spline_1d_interpolator_nonuniform),intent(inout) :: interpolator
#endif
    sll_real64, dimension(:), intent(in), optional :: coeffs
    print *, 'set_coefficients_qs1d_nonuniform(): ERROR: This function has not been ', &
         'implemented yet.'
    if(present(coeffs))then
      print *,'#coeffs present but not used'
    endif
    print *,interpolator%num_points     
    stop
  end subroutine set_coefficients_qs1d_nonuniform


  function get_coefficients_qs1d_nonuniform(interpolator)
#ifdef STDF95
    type(quintic_spline_1d_interpolator_nonuniform),intent(in)  :: interpolator
#else
    class(quintic_spline_1d_interpolator_nonuniform),intent(in) :: interpolator
#endif
    sll_real64, dimension(:), pointer   :: get_coefficients_qs1d_nonuniform    
    
    print *, 'get_coefficients_qs1d_nonuniform(): ERROR: This function has not been ', &
         'implemented yet.'
    get_coefficients_qs1d_nonuniform = 0._f64
    print *,interpolator%num_points
    stop      
  end function get_coefficients_qs1d_nonuniform


end module sll_quintic_spline_interpolator_1d_nonuniform
