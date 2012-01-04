!> \brief  
!> The splines module provides capabilities for 1D data interpolation with cubic B-splines
!> on non uniform mesh and different boundary conditions
!> (at the time of this writing: periodic). The data to be interpolated is represented by a 
!> simple array.
!> contact: mehrenbe@math.unistra.fr for this module
!> 
module cubic_nonuniform_splines
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
  implicit none

  type cubic_nonunif_spline_1D
    sll_int32                         :: n_points ! size
    sll_real64, dimension(:), pointer :: mesh     ! the non uniform mesh
    sll_real64, dimension(:), pointer :: buf      ! memory buffer
    sll_int32                         :: sizebuf  ! size of the buffer 
    sll_real64, dimension(:), pointer :: coeffs   ! the spline coefficients
    sll_int32                         :: bc_type 
  end type cubic_nonunif_spline_1D

  enum, bind(C)
     enumerator :: PERIODIC_SPLINE = 0, HERMITE_SPLINE = 1
  end enum
  
  
contains  ! ****************************************************************
  !> create new spline object
  function new_cubic_nonunif_spline_1D( n_points, mesh, bc_type )
    type(cubic_nonunif_spline_1D), pointer        :: new_cubic_nonunif_spline_1D
    sll_int32,  intent(in)                        :: n_points
    sll_real64, dimension(:), pointer             :: mesh
    sll_int32,  intent(in)                        :: bc_type
    sll_int32                                     :: ierr
    
    SLL_ALLOCATE( new_cubic_nonunif_spline_1D, ierr )
    new_cubic_nonunif_spline_1D%n_points = n_points
    new_cubic_nonunif_spline_1D%bc_type  = bc_type
    new_cubic_nonunif_spline_1D%mesh  = mesh
    
  end function new_cubic_nonunif_spline_1D

  !> delete spline object
  subroutine delete_cubic_nonunif_spline_1D( )
  end subroutine delete_cubic_nonunif_spline_1D
  
  !> compute splines coefficients
  subroutine compute_spline_nonunif( f, bc_type, spline )
    sll_real64, dimension(:), intent(in) :: f    ! data to be fit
    sll_int32,  intent(in)               :: bc_type
    type(cubic_nonunif_spline_1D), pointer         :: spline

    select case (bc_type)
    case (PERIODIC_SPLINE)
       !call compute_spline_nonunif_1D_periodic( f, spline )
    case default
       print *, 'ERROR: compute_spline_nonunif_1D(): not recognized boundary condition'
       STOP
    end select  
  end subroutine compute_spline_nonunif
  
  !> get spline interpolate at point x
  function interpolate_value_nonunif( x, spline )
    sll_real64                                     :: interpolate_value_nonunif
    sll_real64, intent(in)                         :: x
    type(cubic_nonunif_spline_1D), pointer      :: spline
    
    interpolate_value_nonunif = 0.0_f64
  end function interpolate_value_nonunif

end module cubic_nonuniform_splines

