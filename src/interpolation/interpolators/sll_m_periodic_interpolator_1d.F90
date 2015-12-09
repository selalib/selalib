!> @ingroup interpolators
!> @brief
!> Interpolator with periodic boundary conditions
!! @details
!> the following provides an implementation for the abstract interface
!! sll_c_interpolator_1d
!! Define periodic interpolation of values in data define on original grid at
!! points coordinates
module sll_m_periodic_interpolator_1d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_interpolators_1d_base, only: &
    sll_c_interpolator_1d

  use sll_m_periodic_interp, only: &
    delete, &
    initialize_periodic_interp, &
    periodic_interp, &
    periodic_interp_work

  implicit none

  public :: &
    sll_delete, &
    sll_periodic_interpolator_1d

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Periodic interpolator
  type, extends(sll_c_interpolator_1d) ::  sll_periodic_interpolator_1d
    ! Be careful here. For consistency with the other interpolators
    ! num_points is the number of nodes (including both boundaries)
    ! and not the number of cells as used in the periodic interpolator module.
     sll_int32                            :: num_points !< size
     sll_real64                           :: cell_size  !< cell size
     sll_real64                           :: domain_size!< length of interval
     type(periodic_interp_work), pointer  :: per_interp !< ???
   contains
     !>PLEASE ADD DOCUMENTATION
     procedure, pass(interpolator) :: initialize => initialize_per1d_interpolator
     !>PLEASE ADD DOCUMENTATION
     procedure :: compute_interpolants => compute_interpolants_per1d
     !>PLEASE ADD DOCUMENTATION
     procedure :: interpolate_from_interpolant_value => interpolate_value_per1d
     !>PLEASE ADD DOCUMENTATION
     procedure :: interpolate_from_interpolant_derivative_eta1 => interpolate_deriv1_per1d
     !>PLEASE ADD DOCUMENTATION
     procedure :: interpolate_from_interpolant_array => interpolate_values_per1d
     !>PLEASE ADD DOCUMENTATION
     procedure, pass:: interpolate_array => per_interpolate1d
     !>PLEASE ADD DOCUMENTATION
     procedure, pass:: interpolate_array_disp => per_interpolate1d_disp
     !>PLEASE ADD DOCUMENTATION
     procedure, pass:: interpolate_array_disp_inplace => per_interpolate1d_disp_inplace
     !>PLEASE ADD DOCUMENTATION
     procedure, pass :: set_coefficients => set_coefficients_per1d
     !>PLEASE ADD DOCUMENTATION
     procedure, pass :: get_coefficients => get_coefficients_per1d

  end type sll_periodic_interpolator_1d

  !> Deallocate the interpolator object
  interface sll_delete
     module procedure delete_per1d
  end interface sll_delete


contains  ! ****************************************************************



  !> Create a new interpolator
  function new_periodic_1d_interpolator( &
    num_points, &
    xmin, &
    xmax, &
    type, &
    order)  result(res)

    type(sll_periodic_interpolator_1d),  pointer :: res
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





  subroutine per_interpolate1d(this, num_pts, data, coordinates, output_array)
    class(sll_periodic_interpolator_1d),  intent(in)       :: this
    !class(sll_spline_1D),  intent(in)      :: this
    sll_int32,  intent(in)                 :: num_pts
    sll_real64, dimension(num_pts), intent(in)   :: coordinates
    sll_real64, dimension(:), intent(in)   :: data
    sll_real64, dimension(num_pts), intent(out)      :: output_array
    ! local variables
    !sll_int32 :: ierr
    ! periodic interpolation only implemented for constant displacement
    print*, 'periodic_interpolate1d: periodic interpolation not implemented', &
         ' for array of displacements'
    output_array = -1000000._f64
    print *,num_pts
    print *,maxval(coordinates)
    print *,maxval(data)
    print *,maxval(output_array)
    print *,this%num_points
    stop
  end subroutine per_interpolate1d

  subroutine per_interpolate1d_disp(this, num_pts, data, alpha, output_array)
    class(sll_periodic_interpolator_1d),  intent(in)       :: this
    sll_int32,  intent(in)                 :: num_pts
    sll_real64,  intent(in)   :: alpha
    sll_real64, dimension(:), intent(in)   :: data
    sll_real64, dimension(num_pts), intent(out)      :: output_array
    ! Be careful here. For consistency with the other interpolators
    ! num_points is the number of nodes (including both boundaries)
    ! and not the number of cells as used in the periodic interpolator module.
    call periodic_interp(this%per_interp, output_array, data, &
         -alpha/this%cell_size)
    ! complete by periodicity
    output_array(num_pts) = output_array(1)
  end subroutine per_interpolate1d_disp


  subroutine per_interpolate1d_disp_inplace(this, num_pts, data, alpha)
    class(sll_periodic_interpolator_1d),  intent(in)       :: this
    sll_int32,  intent(in)                 :: num_pts
    sll_real64,  intent(in)   :: alpha
    sll_real64, dimension(num_pts), intent(inout)   :: data

    ! local variable
    sll_real64 :: tmp(num_pts)
    ! Be careful here. For consistency with the other interpolators
    ! num_points is the number of nodes (including both boundaries)
    ! and not the number of cells as used in the periodic interpolator module.
    call periodic_interp(this%per_interp, tmp, data, &
         -alpha/this%cell_size)
    ! complete by periodicity
    data = tmp
    data(num_pts) = tmp(1)
  end subroutine per_interpolate1d_disp_inplace


  ! Both versions F03 and F95 of compute_interpolants_per1d should have the
  ! same name. In the F95 we should add a generic interface around this
  ! subroutine, selecting on the type of interpolator. In the F03 case the
  ! interface is the compute_interpolants routine which gets assigned to
  ! the per1d at initialization time.
  subroutine compute_interpolants_per1d( interpolator, data_array,&
       eta_coords, &
       size_eta_coords)

    class(sll_periodic_interpolator_1d), intent(inout) :: interpolator
    sll_real64, dimension(:), intent(in)               :: data_array
    sll_real64, dimension(:), intent(in),optional  :: eta_coords
    sll_int32, intent(in),optional                 :: size_eta_coords
    print*, 'compute_interpolants_per1d:', &
         ' not implemented for periodic interpolation'
    if(present(eta_coords))then
      print *,'eta_coords present but not used'
    endif
    if(present(size_eta_coords))then
      print *,'size_eta_coords present but not used'
    endif
    print *,maxval(data_array)
    print *,interpolator%num_points
    stop
  end subroutine compute_interpolants_per1d

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
    class(sll_periodic_interpolator_1d),  intent(in) :: interpolator
    sll_int32,  intent(in)                 :: num_pts
    sll_real64, dimension(num_pts), intent(in)   :: vals_to_interpolate
    sll_real64, dimension(num_pts), intent(out)  :: output_array
    !sll_int32 :: ierr
    output_array = -1000000._f64
    print*, 'interpolate_values_per1d:', &
         ' not implemented for periodic interpolation'
    print *,interpolator%num_points
    print *,num_pts
    print *,maxval(vals_to_interpolate)
    stop
  end subroutine interpolate_values_per1d


!PN DEFINED BUT NOT USED
! subroutine interpolate_derivatives_per1d( &
!   interpolator, &
!   num_pts, &
!   vals_to_interpolate, &
!   output_array )


!   class(sll_periodic_interpolator_1d),  intent(in) :: interpolator
!   sll_int32,  intent(in)                 :: num_pts
!   sll_real64, dimension(:), intent(in)   :: vals_to_interpolate
!   sll_real64, dimension(:), intent(out)  :: output_array
!   !sll_int32 :: ierr

!   print*, 'interpolate_from_interpolant_derivatives_eta1: ', &
!        'not implemented for periodic interpolation'
!   output_array = -1000000._f64
!   print *,interpolator%num_points
!   print *,num_pts
!   print *,maxval(vals_to_interpolate)
!   stop
! end subroutine interpolate_derivatives_per1d


  function interpolate_value_per1d( interpolator, eta1 ) result(val)
    class(sll_periodic_interpolator_1d), intent(in) :: interpolator
    sll_real64 :: val
    sll_real64, intent(in) :: eta1
     print*, 'interpolate_value_per1d: ', &
         'not implemented for periodic interpolation'
    val = -1000000._f64
    print *,eta1
    print *,interpolator%num_points
    stop
  end function

  function interpolate_deriv1_per1d( interpolator, eta1 ) result(val)
    class(sll_periodic_interpolator_1d), intent(in) :: interpolator
    sll_real64             :: val
    sll_real64, intent(in) :: eta1
     print*, 'interpolate_deriv1_per1d: ', &
         'not implemented for periodic interpolation'
    val = -1000000._f64
    print *,eta1
    print *,interpolator%num_points
    stop
  end function interpolate_deriv1_per1d

!PN DEFINED BUT NOT USED
! function interpolate_derivative_f95( interpolator, eta1 ) result(val)
!   class(sll_periodic_interpolator_1d), intent(in) :: interpolator
!   sll_real64 :: val
!   sll_real64, intent(in) :: eta1
!    print*, 'interpolate_derivative_f95: ', &
!        'not implemented for periodic interpolation'
!   val = -1000000._f64
!   print *,eta1
!   print *,interpolator%num_points
!   stop
! end function interpolate_derivative_f95


  ! Why is the name of this function changing depending on the standard?
  ! only one will be compiled anyway!!

  !> initialize periodic interpolator
  subroutine initialize_per1d_interpolator( &
    interpolator, &
    num_points, &
    xmin, &
    xmax, &
    type, &
    order)

    class(sll_periodic_interpolator_1d),  intent(inout) :: interpolator
    sll_int32,  intent(in)               :: num_points
    sll_real64, intent(in)               :: xmin
    sll_real64, intent(in)               :: xmax
    sll_int32,  intent(in)               :: type
    sll_int32,  intent(in)               :: order
    !sll_int32                            :: ierr
    !sll_int32  :: i
    !sll_real64 :: delta

    ! Be careful here. For consistency with the other interpolators
    ! num_points is the number of nodes (including both boundaries)
    ! and not the number of cells as used in the periodic interpolator module.
    interpolator%num_points = num_points - 1
    interpolator%cell_size  = (xmax-xmin) / (num_points-1)
    interpolator%domain_size = xmax-xmin

    call initialize_periodic_interp(interpolator%per_interp, num_points-1, &
         type, order)
  end subroutine


  subroutine delete_per1d( obj )
    class(sll_periodic_interpolator_1d) :: obj
    call delete(obj%per_interp)
  end subroutine delete_per1d

  subroutine set_coefficients_per1d( interpolator, coeffs )
    class(sll_periodic_interpolator_1d), intent(inout) :: interpolator
    sll_real64, dimension(:), intent(in), optional :: coeffs
    print *, 'set_coefficients_per1d(): ERROR: This function has not been ', &
         'implemented yet.'
    if(present(coeffs))then
      print *,'#coeffs present but not used'
    endif
    print *,interpolator%num_points
    stop
  end subroutine set_coefficients_per1d


  function get_coefficients_per1d(interpolator)
    class(sll_periodic_interpolator_1d), intent(in) :: interpolator
    sll_real64, dimension(:), pointer            :: get_coefficients_per1d

    print *, 'get_coefficients_per1d(): ERROR: This function has not been ', &
         'implemented yet.'
    print *,interpolator%num_points
    get_coefficients_per1d => null()
    stop
  end function get_coefficients_per1d

end module sll_m_periodic_interpolator_1d
