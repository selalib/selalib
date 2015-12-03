!> @ingroup interpolators
!> @brief
!> Interpolator class and methods of Lagrange 1D interpolator
!> @details
!> Implements the sll_c_interpolator_1d interface.
module sll_m_lagrange_interpolator_1d
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
use sll_m_interpolators_1d_base
use sll_m_lagrange_interpolation_1d

implicit none
private

 !> Interpolator class of Lagrange 1D interpolator
 type,extends(sll_c_interpolator_1d), public :: sll_lagrange_interpolator_1d
   !> PLEASE ADD DOCUMENTATION
   type(sll_lagrange_interpolation_1D), pointer :: lagrange
   !> PLEASE ADD DOCUMENTATION
   sll_int32                                    :: bc_type
   contains
   !> PLEASE ADD DOCUMENTATION
   procedure,pass(interpolator) :: initialize => initialize_li1d_interpolator
   !> PLEASE ADD DOCUMENTATION
   procedure :: compute_interpolants => compute_interpolants_li1d
   !> PLEASE ADD DOCUMENTATION
   !procedure :: interpolate_from_interpolant_derivatives_eta1 => interpolate_array_derivatives_li1d
   !> PLEASE ADD DOCUMENTATION
   procedure :: interpolate_array => interpolate_array_li1d
   !> PLEASE ADD DOCUMENTATION
   procedure :: interpolate_array_disp => interpolate_array_disp_li1d
   !> PLEASE ADD DOCUMENTATION
   procedure :: interpolate_array_disp_inplace => interpolate_array_disp_inplace_li1d
   !> PLEASE ADD DOCUMENTATION
   !procedure :: interpolate_pointer_derivatives => interpolate_pointer_derivatives_li1d
   !> PLEASE ADD DOCUMENTATION
   procedure :: interpolate_from_interpolant_derivative_eta1 => interpolate_derivative_eta1_li1d
   !> PLEASE ADD DOCUMENTATION
   !procedure :: interpolate_pointer_values => interpolate_pointer_values_li1d
   !> PLEASE ADD DOCUMENTATION
   procedure :: interpolate_from_interpolant_array => interpolate_array_values_li1d
   !> PLEASE ADD DOCUMENTATION
   procedure :: interpolate_from_interpolant_value => interpolate_value_li1d
   !> PLEASE ADD DOCUMENTATION
   !procedure :: reconstruct_array => reconstruct_array_li1d
   !> PLEASE ADD DOCUMENTATION
   procedure, pass :: set_coefficients => set_coefficients_li1d
   !> PLEASE ADD DOCUMENTATION
   procedure, pass :: get_coefficients => get_coefficients_li1d
 end type sll_lagrange_interpolator_1d

!PN DEFINED BUT NOT USED
! !> Deallocate the class interpolator
! interface sll_delete
!   module procedure delete_li1d
! end interface

 public new_lagrange_interpolator_1d

contains  !**********************************************************

subroutine initialize_li1d_interpolator(interpolator,num_points,xmin,xmax,bc_type,d, periodic_last)
  class(sll_lagrange_interpolator_1d), intent(inout) :: interpolator
    sll_int32, intent(in)                        :: d,num_points,bc_type
    sll_real64, intent(in)                       :: xmin,xmax
    sll_int32, intent(in), optional              :: periodic_last
    sll_int32                                    :: last

    if (present(periodic_last)) then
       last  = periodic_last
    else
       last = 1
    end if

    interpolator%lagrange => new_lagrange_interpolation_1D( &
           num_points, &
           xmin, &
           xmax, &
           bc_type, &
           d, &
           last)
    call compute_lagrange_interpolation_1D(interpolator%lagrange)
end subroutine

function new_lagrange_interpolator_1d( &
    num_points, &
    xmin, &
    xmax, &
    bc_type, &
    d, &
    periodic_last) result(res)

    type(sll_lagrange_interpolator_1d),  pointer :: res
    sll_int32,  intent(in)               :: num_points
    sll_real64, intent(in)               :: xmin
    sll_real64, intent(in)               :: xmax
    sll_int32,  intent(in)               :: bc_type
    sll_int32, intent(in)               :: d
    sll_int32 :: ierr
    sll_int32, intent(in), optional              :: periodic_last
    sll_int32                                    :: last

    SLL_ALLOCATE(res,ierr)

    if (present(periodic_last)) then
       last  = periodic_last
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
  class(sll_lagrange_interpolator_1d), intent(in)     :: this
  sll_real64, intent(in) :: alpha
  sll_int32, intent(in)  :: num_pts    ! size of output array
  sll_real64, dimension(:), intent(in) :: data  ! data to be interpolated points where output is desired
  sll_real64, dimension(1:num_pts), intent(out)    :: output_array

  call interpolate_from_interpolant_array(data,-alpha,this%lagrange)
  output_array=this%lagrange%data_out


  !select case (this%bc_type)
  !   case (SLL_PERIODIC)
  !      call lagrange_periodic(data, output_array, alpha, this%stencil_width)
  !   end select
  !call lagrange_halo_cells(data, output_array, alpha, this%stencil_width)

end subroutine interpolate_array_disp_li1d


subroutine interpolate_array_disp_inplace_li1d(this, num_pts, data, alpha)
  class(sll_lagrange_interpolator_1d), intent(in)     :: this
  sll_real64, intent(in) :: alpha
  sll_int32, intent(in)  :: num_pts    ! size of output array
  sll_real64, dimension(num_pts), intent(inout) :: data  ! data to be interpolated points where output is desired

  call interpolate_from_interpolant_array(data,-alpha,this%lagrange)
  data=this%lagrange%data_out

end subroutine interpolate_array_disp_inplace_li1d



!PN DEFINED BUT NOT USED
!subroutine delete_li1d (obj)
!  class(sll_lagrange_interpolator_1d) :: obj
!  call delete(obj%lagrange)
!end subroutine delete_li1d


subroutine interpolate_array_values_li1d( &
    interpolator, &
    num_pts, &
    vals_to_interpolate, &
    output_array )
    class(sll_lagrange_interpolator_1d),  intent(in) :: interpolator
    sll_int32,  intent(in)                 :: num_pts
    sll_real64, dimension(num_pts), intent(in)   :: vals_to_interpolate
    sll_real64, dimension(num_pts), intent(out)  :: output_array
    !sll_int32 :: ierr
    output_array = 0.0_f64
    print*, 'interpolate_from_interpolant_array:', &
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
    class(sll_lagrange_interpolator_1d),  intent(in) :: interpolator
    sll_int32,  intent(in)                 :: num_pts
    sll_real64, dimension(:), intent(in)   :: vals_to_interpolate
    sll_real64, dimension(:), intent(out)  :: output_array
    !sll_int32 :: ierr
    output_array = 0.0_f64
    print*, 'interpolate_from_interpolant_derivatives_eta1: ', &
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
    class(sll_lagrange_interpolator_1d),  intent(in) :: interpolator
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

  function interpolate_derivative_eta1_li1d( interpolator, eta1 ) result(val)
    class(sll_lagrange_interpolator_1d), intent(in) :: interpolator
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
    class(sll_lagrange_interpolator_1d),  intent(in) :: interpolator
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

  function interpolate_value_li1d( interpolator, eta1 ) result(val)
    class(sll_lagrange_interpolator_1d), intent(in) :: interpolator
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
    class(sll_lagrange_interpolator_1d), intent(in)     :: this
       sll_int32, intent(in)                :: num_points! size of output array
       sll_real64, dimension(:), intent(in) :: data   ! data to be interpolated
       sll_real64, dimension(num_points)    :: res
       print *,'#warning reconstruct_array_li1d dummy function'
       print *,num_points
       print *,maxval(data)
       print *,this%bc_type
       res(:) = 0.0_f64
  end function reconstruct_array_li1d

  subroutine interpolate_array_li1d(this, num_pts, data, coordinates, output_array)
    class(sll_lagrange_interpolator_1d),  intent(in)       :: this
    !class(sll_spline_1D),  intent(in)      :: this
    sll_int32,  intent(in)                 :: num_pts
    sll_real64, dimension(num_pts), intent(in)   :: coordinates
    sll_real64, dimension(:), intent(in)   :: data
    sll_real64, dimension(num_pts), intent(out)      :: output_array
    ! local variables
    !sll_int32 :: ierr
    ! lagrange interpolation only implemented for constant displacement
    print*, 'interpolate_array_li1d: ', &
         'not implemented for lagrange interpolation'
    print *,maxval(coordinates)
    print *,maxval(data)
    print *,this%bc_type
    output_array = 0._f64
    stop
  end subroutine interpolate_array_li1d

    subroutine compute_interpolants_li1d( interpolator, data_array,&
         eta_coords, &
         size_eta_coords)
      class(sll_lagrange_interpolator_1d), intent(inout) :: interpolator
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
    class(sll_lagrange_interpolator_1d),  intent(inout) :: interpolator
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
    class(sll_lagrange_interpolator_1d),  intent(in) :: interpolator
    sll_real64, dimension(:), pointer            :: get_coefficients_li1d

    print *, 'get_coefficients_li1d(): ERROR: This function has not been ', &
         'implemented yet.'
    print *,interpolator%bc_type
    get_coefficients_li1d => null()
    stop
  end function get_coefficients_li1d
  !DEFINE_NULL_INTERP_1D_ARRAY_SUB(sll_lagrange_interpolator_1d, interpolate_array_values_li1d)
!DEFINE_NULL_INTERP_1D_ARRAY_SUB(sll_lagrange_interpolator_1d, interpolate_array_derivatives_li1d)
!DEFINE_NULL_INTERP_1D_POINTER_SUB(sll_lagrange_interpolator_1d, interpolate_pointer_derivatives_li1d)
!DEFINE_NULL_INTERP_ONE_ARG_MSG(sll_lagrange_interpolator_1d, interpolate_derivative_eta1_li1d)
!DEFINE_NULL_INTERP_1D_POINTER_SUB(sll_lagrange_interpolator_1d, interpolate_pointer_values_li1d)
!DEFINE_NULL_INTERP_ONE_ARG_MSG(sll_lagrange_interpolator_1d, interpolate_value_li1d)
!DEFINE_NULL_RECONSTRUCT_1D_ARRAY(sll_lagrange_interpolator_1d, reconstruct_array_li1d)
!DEFINE_NULL_INTERP_1D_ARRAY(sll_lagrange_interpolator_1d, interpolate_array_li1d)
!DEFINE_NULL_INTERP_1D_ARRAY_MSG(sll_lagrange_interpolator_1d, compute_interpolants_li1d)

end module sll_m_lagrange_interpolator_1d
