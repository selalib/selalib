module sll_cubic_spline_interpolator_1d
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
use sll_interpolator_1d
use sll_splines
  implicit none
  
  type, extends(interpolator_1d_base) ::  cubic_spline_1d_interpolator
     sll_int32            :: num_points ! size
     sll_int32            :: bc_type
     type(sll_spline_1D), pointer  :: spline
   contains
     procedure, pass:: initialize_cs1d_interpolator
     procedure, pass:: interpolate_array => spline_interpolate1d
     procedure, pass:: reconstruct_array
     generic :: initialize => initialize_cs1d_interpolator
  end type cubic_spline_1d_interpolator


  
contains  ! ****************************************************************


  ! the following provides an implementation for the abstract interface interpolate1d
  !> Define spline interpolation of values in data define on original grid at points coordinates
  function spline_interpolate1d(this, num_points, data, coordinates) result(data_out)
    class(cubic_spline_1d_interpolator),  intent(in)       :: this
    !class(sll_spline_1D),  intent(in)      :: this
    sll_int32,  intent(in)                 :: num_points
    sll_real64, dimension(:), intent(in)   :: coordinates
    sll_real64, dimension(:), intent(in)   :: data
    sll_real64, dimension(num_points)      :: data_out
    ! local variables
    sll_int32 :: ierr
    ! compute the interpolating spline coefficients
    call compute_spline_1D( data, this%bc_type, this%spline )
    call interpolate_array_values( coordinates, data_out, num_points, this%spline )
    
  end function spline_interpolate1d
  
  
  !> create new spline object
  subroutine initialize_cs1d_interpolator( this, num_points, xmin, xmax, bc_type, sl, sr )
    class(cubic_spline_1d_interpolator),  intent(inout)       :: this
    sll_int32,  intent(in)               :: num_points
    sll_real64, intent(in)               :: xmin
    sll_real64, intent(in)               :: xmax
    sll_int32,  intent(in)               :: bc_type
    sll_real64, intent(in), optional     :: sl
    sll_real64, intent(in), optional     :: sr
    sll_int32                            :: ierr
    
    this%num_points = num_points
    this%bc_type = bc_type
    if (present(sl).and.present(sr)) then
       this%spline => new_spline_1D(num_points, xmin, xmax, bc_type, sl, sr )
    else
       this%spline => new_spline_1D(num_points, xmin, xmax, bc_type)
    end if
  end subroutine initialize_cs1d_interpolator

  function reconstruct_array(this, num_points, data) result(res)
    ! dummy procedure
    class(cubic_spline_1d_interpolator), intent(in)     :: this
       sll_int32, intent(in)                 :: num_points    ! size of output array
       sll_real64, dimension(:), intent(in)  :: data          ! data to be interpolated 
       sll_real64, dimension(num_points)     :: res
  end function reconstruct_array
  
end module sll_cubic_spline_interpolator_1d
