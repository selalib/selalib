!> @ingroup interpolators
!> @brief
!> Interpolator class and methods of hermite 1D interpolator
!> @details
!> Implements the sll_c_interpolator_1d interface.
module sll_m_hermite_interpolator_1d
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
use sll_m_interpolators_1d_base
use sll_m_hermite_interpolation_1d
implicit none
private

!> @brief
!> The hermite-based interpolator is only a wrapper around the capabilities
!! of the hermite interpolation. 
!> @details
!! All interpolators share a common interface with
!! respect to their use, as described by the interpolator_1d_base class.
!! Where the diverse interpolators diverge is in the way to initialize them.

 !> Interpolator class of Hermite 1D interpolator
 type,extends(sll_c_interpolator_1d), public :: sll_hermite_interpolator_1d
   !> PLEASE ADD DOCUMENTATION
   type(sll_hermite_interpolation_1d), pointer :: hermite
   !> PLEASE ADD DOCUMENTATION
   sll_int32                                    :: npts
   contains
   !> PLEASE ADD DOCUMENTATION
   procedure,pass(interpolator) :: initialize => initialize_hermite_interpolator_1d
   !> PLEASE ADD DOCUMENTATION
   procedure :: compute_interpolants => wrap_compute_interpolants_hermite_1d
   !> PLEASE ADD DOCUMENTATION
   procedure :: interpolate_from_interpolant_derivatives_eta1 => interpolate_array_derivatives_hi1d
   !> PLEASE ADD DOCUMENTATION
   procedure :: interpolate_array => wrap_interpolate_array_hermite_1d

   !> PLEASE ADD DOCUMENTATION
   procedure :: interpolate_array_disp => interpolate_array_disp_hi1d
   !> PLEASE ADD DOCUMENTATION
   procedure :: interpolate_array_disp_inplace => interpolate_array_disp_inplace_hi1d
   !> PLEASE ADD DOCUMENTATION
   procedure :: interpolate_from_interpolant_derivative_eta1 => interpolate_derivative_eta1_hi1d
   !> PLEASE ADD DOCUMENTATION
   procedure :: interpolate_from_interpolant_array => interpolate_array_values_hi1d
   !> PLEASE ADD DOCUMENTATION
   procedure :: interpolate_from_interpolant_value => wrap_interpolate_value_hermite_1d
   !> PLEASE ADD DOCUMENTATION
   procedure, pass :: set_coefficients => set_coefficients_hi1d
   !> PLEASE ADD DOCUMENTATION
   procedure, pass :: get_coefficients => get_coefficients_hi1d
 end type sll_hermite_interpolator_1d

!PN DEFINED BUT NOT USED
! !> Deallocate the class interpolator
! interface sll_delete
!   module procedure delete_hi1d
! end interface

 public new_hermite_interpolator_1d

contains  !**********************************************************





    !> PLEASE ADD DOCUMENTATION
  function new_hermite_interpolator_1d( &
    npts, &
    eta_min, &
    eta_max, &
    degree, &
    eta_hermite_continuity, &
    eta_bc_type, &
    const_eta_min_slope, &
    const_eta_max_slope, &
    eta_min_slopes, &
    eta_max_slopes ) &
    result(interpolator)

    type(sll_hermite_interpolator_1d), pointer :: interpolator
    sll_int32, intent(in) :: npts
    sll_real64, intent(in) :: eta_min
    sll_real64, intent(in) :: eta_max
    sll_int32, intent(in) :: degree
    sll_int32, intent(in) :: eta_hermite_continuity
    sll_int32, intent(in) :: eta_bc_type
    sll_real64, intent(in), optional :: const_eta_min_slope
    sll_real64, intent(in), optional :: const_eta_max_slope
    sll_real64, dimension(:),intent(in), optional :: eta_min_slopes
    sll_real64, dimension(:),intent(in), optional :: eta_max_slopes
    sll_int32 :: ierr
    
    SLL_ALLOCATE(interpolator,ierr)
    
    interpolator%npts = npts
    
    call interpolator%initialize( &
      npts, &
      eta_min, &
      eta_max, &
      degree, &
      eta_hermite_continuity, &
      eta_bc_type, &
      const_eta_min_slope, &
      const_eta_max_slope, &
      eta_min_slopes, &
      eta_max_slopes)    

     
  end function  new_hermite_interpolator_1d


  subroutine initialize_hermite_interpolator_1d( &
    interpolator, &    
    npts, &
    eta_min, &
    eta_max, &
    degree, &
    eta_hermite_continuity, &
    eta_bc_type, &
    const_eta_min_slope, &
    const_eta_max_slope, &
    eta_min_slopes, &
    eta_max_slopes)    

    class(sll_hermite_interpolator_1d), intent(inout) :: interpolator
    sll_int32, intent(in) :: npts
    sll_real64, intent(in) :: eta_min
    sll_real64, intent(in) :: eta_max
    sll_int32, intent(in) :: degree
    sll_int32, intent(in) :: eta_hermite_continuity
    sll_int32, intent(in) :: eta_bc_type
    sll_real64, intent(in), optional :: const_eta_min_slope
    sll_real64, intent(in), optional :: const_eta_max_slope
    sll_real64, dimension(:),intent(in), optional :: eta_min_slopes
    sll_real64, dimension(:),intent(in), optional :: eta_max_slopes
       
    interpolator%hermite  => new_hermite_interpolation_1d( &
      npts, &
      eta_min, &
      eta_max, &
      degree, &
      eta_hermite_continuity, &
      eta_bc_type, &
      const_eta_min_slope, &
      const_eta_max_slope, &
      eta_min_slopes, &
      eta_max_slopes)    
      
  end subroutine initialize_hermite_interpolator_1d


  subroutine wrap_compute_interpolants_hermite_1d( &
    interpolator, &
    data_array, &
    eta_coords, &
    size_eta_coords)
    class(sll_hermite_interpolator_1d), intent(inout) :: interpolator
    sll_real64, dimension(:), intent(in) :: data_array
    sll_real64, dimension(:), intent(in),optional   :: eta_coords
    sll_int32, intent(in), optional                 :: size_eta_coords

    if(present(eta_coords))then
      !print *,'#Warning eta_coords not used'
    endif
    if(present(size_eta_coords))then
      !print *,'#Warning size_eta_coords not used'
    endif
    call compute_interpolants_hermite_1d( interpolator%hermite, data_array )
  end subroutine wrap_compute_interpolants_hermite_1d
  
  function wrap_interpolate_value_hermite_1d( interpolator, eta1 ) result(val)
    class(sll_hermite_interpolator_1d), intent(in) :: interpolator
    sll_real64 :: val
    sll_real64, intent(in) :: eta1
    val = interpolate_value_hermite_1d( eta1, interpolator%hermite )
      
  end function wrap_interpolate_value_hermite_1d


  subroutine wrap_interpolate_array_hermite_1d( &
    this, &
    num_pts, &
    data, &
    coordinates, &
    output_array)
    class(sll_hermite_interpolator_1d),  intent(in) :: this
    sll_int32,  intent(in)                           :: num_pts
    sll_real64, dimension(num_pts), intent(in)           :: coordinates
    sll_real64, dimension(:), intent(in)           :: data
    sll_real64, dimension(num_pts), intent(out)   :: output_array
    sll_int32 :: i
    call compute_interpolants_hermite_1d( this%hermite, data )
    do i = 1, num_pts
      output_array(i) = this%interpolate_from_interpolant_value(coordinates(i))
    end do
  end subroutine wrap_interpolate_array_hermite_1d




subroutine interpolate_array_disp_hi1d(this, num_pts, data, alpha, output_array)
  class(sll_hermite_interpolator_1d), intent(in)     :: this
  sll_real64, intent(in) :: alpha
  sll_int32, intent(in)  :: num_pts    ! size of output array
  sll_real64, dimension(:), intent(in) :: data  ! data to be interpolated points where output is desired
  sll_real64, dimension(1:num_pts), intent(out) :: output_array

!call interpolate_from_interpolant_array(data,alpha,this%hermite)
!data_out=this%hermite%data_out
  print*, 'interpolate_array_disp_hi1d:', &
       ' not implemented for hermite interpolation'
  SLL_ASSERT(this%npts>0)
  output_array = 0.0_f64 * alpha + data

end subroutine interpolate_array_disp_hi1d



subroutine interpolate_array_disp_inplace_hi1d(this, num_pts, data, alpha)
  class(sll_hermite_interpolator_1d), intent(in)     :: this
  sll_real64, intent(in) :: alpha
  sll_int32, intent(in)  :: num_pts    ! size of output array
  sll_real64, dimension(num_pts), intent(inout) :: data  ! data to be interpolated points where output is desired

!call interpolate_from_interpolant_array(data,alpha,this%hermite)
!data_out=this%hermite%data_out
  print*, 'interpolate_array_disp_inplace_hi1d:', &
       ' not implemented for hermite interpolation'
  SLL_ASSERT(this%npts>0)

end subroutine interpolate_array_disp_inplace_hi1d


!PN DEFINED BUT NOT USED
!subroutine delete_hi1d (obj)
!  class(sll_hermite_interpolator_1d) :: obj
!  !should be fixed; for the moment just commented
!  !call delete(obj%hermite)
!    SLL_ASSERT(obj%npts>0)
!end subroutine delete_hi1d


subroutine interpolate_array_values_hi1d( &
    interpolator, &
    num_pts, &
    vals_to_interpolate, &
    output_array )
    class(sll_hermite_interpolator_1d),  intent(in) :: interpolator
    sll_int32,  intent(in)                 :: num_pts
    sll_real64, dimension(num_pts), intent(in)   :: vals_to_interpolate
    sll_real64, dimension(num_pts), intent(out)  :: output_array
    !sll_int32 :: ierr
    output_array = 0.0_f64
    print*, 'interpolate_from_interpolant_array:', &
         ' not implemented for hermite interpolation'
    print *,num_pts
    print *,maxval(vals_to_interpolate)
    output_array = 0._f64
    !print *,interpolator%bc_type
    stop
    SLL_ASSERT(interpolator%npts>0)
end subroutine interpolate_array_values_hi1d


subroutine interpolate_array_derivatives_hi1d( &
    interpolator, &
    num_pts, &
    vals_to_interpolate, &
    output_array )
    class(sll_hermite_interpolator_1d),  intent(in) :: interpolator
    sll_int32,  intent(in)                 :: num_pts
    sll_real64, dimension(:), intent(in)   :: vals_to_interpolate
    sll_real64, dimension(:), intent(out)  :: output_array
    !sll_int32 :: ierr
    output_array = 0.0_f64
    print*, 'interpolate_from_interpolant_derivatives_eta1: ', &
         'not implemented for hermite interpolation'
    print *,num_pts
    print *,maxval(vals_to_interpolate)
    print *,maxval(output_array)
    !print *,interpolator%bc_type
    stop
    SLL_ASSERT(interpolator%npts>0)
end subroutine interpolate_array_derivatives_hi1d


  function interpolate_derivative_eta1_hi1d( interpolator, eta1 ) result(val)
    class(sll_hermite_interpolator_1d), intent(in) :: interpolator
    sll_real64             :: val
    sll_real64, intent(in) :: eta1
     print*, 'interpolate_derivative_eta1_hi1d: ', &
         'not implemented for hermite interpolation'
    !print *,interpolator%bc_type
    print *,eta1
    val = 0._f64
    stop
    SLL_ASSERT(interpolator%npts>0)
  end function


!PN DEFINED BUT NOT USED
!  function interpolate_value_hi1d( interpolator, eta1 ) result(val)
!    class(sll_hermite_interpolator_1d), intent(in) :: interpolator
!    sll_real64 :: val
!    sll_real64, intent(in) :: eta1
!     print*, 'interpolate_value_hi1d: ', &
!         'not implemented for hermite interpolation'
!    val = 0._f64
!    print *,eta1
!    !print *,interpolator%bc_type
!    stop
!  end function


!PN DEFINED BUT NOT USED
!  function interpolate_array_hi1d(this, num_points, data, coordinates) &
!       result(data_out)
!    class(sll_hermite_interpolator_1d),  intent(in)       :: this
!    !class(sll_spline_1D),  intent(in)      :: this
!    sll_int32,  intent(in)                 :: num_points
!    sll_real64, dimension(:), intent(in)   :: coordinates
!    sll_real64, dimension(:), intent(in)   :: data
!    sll_real64, dimension(num_points)      :: data_out
!    ! local variables
!    !sll_int32 :: ierr
!    print*, 'interpolate_array_hi1d: ', &
!         'not implemented for hermite interpolation'
!    print *,maxval(coordinates)
!    print *,maxval(data)
!    !print *,this%bc_type
!    data_out = 0._f64
!    stop
!  end function

!PN DEFINED BUT NOT USED
!   subroutine compute_interpolants_hi1d( interpolator, data_array,&
!        eta_coords, &
!        size_eta_coords)
!     class(sll_hermite_interpolator_1d), intent(inout) :: interpolator
!   sll_real64, dimension(:), intent(in)               :: data_array
!   sll_real64, dimension(:), intent(in),optional  :: eta_coords
!   sll_int32, intent(in),optional                 :: size_eta_coords
!   print*, 'compute_interpolants_hi1d:', &
!        ' not implemented for hermite interpolation'
!   if(present(eta_coords))then
!     print *,'eta_coords present but not used'
!   endif
!   if(present(size_eta_coords))then
!     print *,'size_eta_coords present but not used'
!   endif
!   print *,maxval(data_array)
!   !print *,interpolator%bc_type
!   stop
! end subroutine

  subroutine set_coefficients_hi1d( interpolator, coeffs )
    class(sll_hermite_interpolator_1d),  intent(inout) :: interpolator
    sll_real64, dimension(:), intent(in), optional :: coeffs
    print *, 'set_coefficients_hi1d(): ERROR: This function has not been ', &
         'implemented yet.'
    if(present(coeffs))then
      print *,'#coeffs present but not used'
    endif
    !print *,interpolator%bc_type
    stop
    SLL_ASSERT(interpolator%npts > 0)
  end subroutine set_coefficients_hi1d


  function get_coefficients_hi1d(interpolator)
    class(sll_hermite_interpolator_1d),  intent(in) :: interpolator
    sll_real64, dimension(:), pointer            :: get_coefficients_hi1d

    print *, 'get_coefficients_hi1d(): ERROR: This function has not been ', &
         'implemented yet.'
    SLL_ASSERT(interpolator%npts>0)
    get_coefficients_hi1d => null()
    stop
  end function get_coefficients_hi1d

end module sll_m_hermite_interpolator_1d
