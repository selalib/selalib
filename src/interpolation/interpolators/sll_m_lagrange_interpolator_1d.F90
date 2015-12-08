!> @ingroup interpolators
!> @brief
!> Interpolator class and methods of Lagrange 1D interpolator
!> @details
!> Implements the sll_c_interpolator_1d interface.
module sll_m_lagrange_interpolator_1d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_errors.h"
#include "sll_working_precision.h"

  use sll_m_boundary_condition_descriptors, only: &
    sll_d_one_sided, &
    sll_d_halo, &
    sll_periodic

  use sll_m_interpolators_1d_base, only: &
    sll_c_interpolator_1d

  use sll_m_lagrange_fast, only : &
    sll_s_interpolate_array_disp_lagrange_fixed_no_bc, &
    sll_s_interpolate_array_disp_lagrange_fixed_periodic, &
    sll_s_interpolate_array_disp_lagrange_fixed_halo_cells

  use sll_m_lagrange_interpolation_1d, only: &
    compute_lagrange_interpolation_1d, &
    interpolate_from_interpolant_array, &
    new_lagrange_interpolation_1d, &
    sll_lagrange_interpolation_1d

  implicit none

  public :: &
    sll_lagrange_interpolator_1d

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 !> Interpolator class of Lagrange 1D interpolator
 type,extends(sll_c_interpolator_1d) :: sll_lagrange_interpolator_1d
   !> PLEASE ADD DOCUMENTATION
   type(sll_lagrange_interpolation_1D), pointer :: lagrange
   !> PLEASE ADD DOCUMENTATION
   sll_int32                                    :: bc_type
   !> Number of points used for interpolation
   sll_int32                                    :: stencil_width 
   !> Flag specifying how the Lagrange interpolation points should be chosen (either SLL_D_INTERP_LAGRANGE_CENTERED or SLL_D_INTERP_LAGRANGE_FIXED)
   sll_int32                                    :: interval_selection
 contains
   !> PLEASE ADD DOCUMENTATION
   procedure,pass(interpolator) :: initialize => initialize_li1d_interpolator
   !> PLEASE ADD DOCUMENTATION
   procedure :: compute_interpolants => compute_interpolants_li1d
   !> PLEASE ADD DOCUMENTATION
   procedure :: interpolate_from_interpolant_derivatives_eta1 => interpolate_array_derivatives_li1d
   !> PLEASE ADD DOCUMENTATION
   procedure :: interpolate_array => interpolate_array_li1d
   !> PLEASE ADD DOCUMENTATION
   procedure :: interpolate_array_disp => interpolate_array_disp_li1d
   !> PLEASE ADD DOCUMENTATION
   procedure :: interpolate_array_disp_inplace => interpolate_array_disp_inplace_li1d
   !> PLEASE ADD DOCUMENTATION
   procedure :: interpolate_from_interpolant_derivative_eta1 => interpolate_derivative_eta1_li1d
   !> PLEASE ADD DOCUMENTATION
   procedure :: interpolate_from_interpolant_array => interpolate_array_values_li1d
   !> PLEASE ADD DOCUMENTATION
   procedure :: interpolate_from_interpolant_value => interpolate_value_li1d
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


  ! Flags for way how to choose the Lagrange points
  sll_int32, parameter :: SLL_D_INTERP_LAGRANGE_CENTERED = 0 !< Flag to specify Lagrange interpolation centered around the interpolation point
  sll_int32, parameter :: SLL_D_INTERP_LAGRANGE_FIXED    = 1 !< Flag to specify Lagrange interpolation on a fixed interval centered around the point that is displace (for interpolate_array_disp)
 

contains  !**********************************************************

subroutine initialize_li1d_interpolator(interpolator,num_points,xmin,xmax,bc_type,d, periodic_last, interval_selection)
  class(sll_lagrange_interpolator_1d), intent(inout) :: interpolator
    sll_int32, intent(in)                        :: d,num_points,bc_type
    sll_real64, intent(in)                       :: xmin,xmax
    sll_int32, intent(in), optional              :: periodic_last
    sll_int32, intent(in), optional              :: interval_selection

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
  
    if (present(interval_selection)) then
       interpolator%interval_selection = interval_selection
    else
       interpolator%interval_selection = SLL_D_INTERP_LAGRANGE_CENTERED
    end if

    select case (interpolator%interval_selection)
    case (SLL_D_INTERP_LAGRANGE_CENTERED)  
       interpolator%stencil_width = 2*d
    case (SLL_D_INTERP_LAGRANGE_FIXED)
       interpolator%stencil_width = 2*d+1
    case default
       SLL_ERROR('interpolate_array_disp_li1d', 'Interval selection not implemented.')
    end select

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

  select case (this%interval_selection)
  case (SLL_D_INTERP_LAGRANGE_CENTERED)  
     call interpolate_from_interpolant_array(data,-alpha,this%lagrange)
     output_array=this%lagrange%data_out
  case (SLL_D_INTERP_LAGRANGE_FIXED)
     select case (this%bc_type)
     case (SLL_PERIODIC)
        call sll_s_interpolate_array_disp_lagrange_fixed_periodic(data, output_array, alpha/this%lagrange%deta, this%stencil_width)
     case (SLL_D_ONE_SIDED)
        call sll_s_interpolate_array_disp_lagrange_fixed_no_bc(data, output_array, alpha/this%lagrange%deta, this%stencil_width)
     case (SLL_D_HALO)
        call sll_s_interpolate_array_disp_lagrange_fixed_halo_cells(data, output_array, alpha/this%lagrange%deta, this%stencil_width)
     case default
        SLL_ERROR('interpolate_array_disp_li1d', 'Boundary type not implemented.')
     end select
  case default
     SLL_ERROR('interpolate_array_disp_li1d', 'Interval selection not implemented.')
  end select

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
