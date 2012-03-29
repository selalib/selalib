module sll_cubic_spline_interpolator_1d
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
use sll_interpolator_1d_base
use sll_splines
  implicit none
  
  type, extends(interpolator_1d) ::  cubic_spline_1d_interpolator
     sll_int32                         :: n_points ! size
     sll_real64                        :: delta    ! discretization step
     sll_real64                        :: rdelta   ! reciprocal of delta
     sll_real64                        :: xmin
     sll_real64                        :: xmax
     sll_int32                         :: bc_type  ! periodic, hermite
     sll_real64, dimension(:), pointer :: data     ! data for the spline fit
     sll_real64, dimension(:), pointer :: d        ! scratch space D (L*D = F),
                                                   ! refer to algorithm below.
                                                   ! Size depends on BC's.
     sll_real64, dimension(:), pointer :: coeffs   ! the spline coefficients
     sll_real64                        :: slope_L  ! left slope, for Hermite
     sll_real64                        :: slope_R  ! right slope, for Hermite
   contains
     procedure, pass:: interpolate_1d => spline_interpolate1d
  end type cubic_spline_1d_interpolator

 
  enum, bind(C)
     enumerator :: PERIODIC_SPLINE = 0, HERMITE_SPLINE = 1
  end enum
  
contains  ! ****************************************************************


  ! the following provides an implementation for the abstract interface interpolate1d
  !> Define spline interpolation of values in data define on original grid at points coordinates
  function spline_interpolate1d(this, num_points, data, coordinates) result(data_out)
    class(sll_spline_1D),  intent(in)       :: this
    !class(sll_spline_1D),  intent(in)      :: this
    sll_int32,  intent(in)                 :: num_points
    sll_real64, dimension(:), intent(in)   :: coordinates
    sll_real64, dimension(:), intent(in)   :: data
    sll_real64, dimension(num_points)      :: data_out
    ! local variables
    sll_int32 :: ierr
    ! compute the interpolating spline coefficients
    call compute_spline_1D( data, this%bc_type, this )
    call interpolate_array_values( coordinates, data_out, num_points, this )
    
  end function spline_interpolate1d
  
  
  !> create new spline object
  function new_spline_1D( num_points, xmin, xmax, bc_type, sl, sr )
    type(sll_spline_1D)                  :: new_spline_1D
    sll_int32,  intent(in)               :: num_points
    sll_real64, intent(in)               :: xmin
    sll_real64, intent(in)               :: xmax
    sll_int32,  intent(in)               :: bc_type
    sll_real64, intent(in), optional     :: sl
    sll_real64, intent(in), optional     :: sr
    sll_int32                            :: ierr
    !SLL_ALLOCATE( new_spline_1D, ierr )
    new_spline_1D%n_points = num_points
    new_spline_1D%xmin     = xmin
    new_spline_1D%xmax     = xmax
    new_spline_1D%delta    = (xmax - xmin)/real((num_points-1),f64)
    new_spline_1D%rdelta   = 1.0_f64/new_spline_1D%delta
    new_spline_1D%bc_type  = bc_type
    if( num_points .le. 28 ) then
       print *, 'ERROR, new_spline_1D: Because of the algorithm used, ', &
            'this function is meant to be used with arrays that are at ', &
            'least of size = 28'
       STOP 'new_spline_1D()'
    end if
    if( xmin .gt. xmax ) then
       print *, 'ERROR, new_spline_1D: xmin is greater than xmax, ', &
            'this would cause all sorts of errors.'
       STOP
    end if
    ! Some more general error checking depending on the type of boundary
    ! condition requested.
    select case (bc_type)
    case (PERIODIC_SPLINE)
       if( present(sl) .or. present(sr) ) then
          print *, 'new_spline_1D(): it is not allowed to specify the ',&
               'end ifin the case of periodic boundary conditions. ', &
               'Exiting program...'
          STOP 'new_spline_1D'
       else
          ! Assign some value, but this value should never be used in the
          ! periodic case anyway.
          new_spline_1D%slope_L = 0.0
          new_spline_1D%slope_R = 0.0
       end if
    case (HERMITE_SPLINE)
       if( present(sl) ) then
          new_spline_1D%slope_L = sl
       else
          new_spline_1D%slope_L = 0.0  ! default left slope for Hermite case
       end if
       if( present(sr) ) then
          new_spline_1D%slope_R = sr
       else
          new_spline_1D%slope_R = 0.0  ! default right slope for Hermite case
       end if
    case default
       print *, 'ERROR: compute_spline_1D(): not recognized boundary condition'
       STOP
    end select
    SLL_ALLOCATE( new_spline_1D%d(num_points),   ierr )
    ! note how the indexing of the coefficients array includes the end-
    ! points 0, num_points, num_points+1, num_points+2. These are meant to 
    ! store the boundary condition-specific data. The 'periodic' BC does
    ! not use the num_points+2 point.
    SLL_ALLOCATE( new_spline_1D%coeffs(0:num_points+2), ierr )

  end function new_spline_1D
  
end module sll_cubic_spline_interpolator_1d
