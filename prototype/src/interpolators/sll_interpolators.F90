!> \brief data type for 1D interpolation
!>
!> 


module sll_interpolator_1d
#include "sll_working_precision.h"
#include "sll_memory.h"
  implicit none

  type :: interpolator_1d
   contains
     procedure, pass :: interpolate_1d
     procedure, pass :: reconstruct_1d
  end type interpolator_1d

!!$  abstract interface
!!$     function interpolate_1d_abstract(this, num_points, data, coordinates) result(res)
!!$       use sll_working_precision
!!$       import interpolator_1d_t
!!$       class(interpolator_1d_t), intent(in)     :: this
!!$       sll_int32, intent(in)                 :: num_points    ! number of interpolation points
!!$       sll_real64, dimension(:), intent(in)  :: data          ! data to be interpolated 
!!$       sll_real64, dimension(:), intent(in)  :: coordinates   ! points where output is desired
!!$       sll_real64, dimension(num_points)     :: res
!!$     end function interpolate_1d_abstract
!!$  end interface

contains
  ! Provides a dummy interpolate_1d function defining the interface
  ! should never be used
  function interpolate_1d(this, num_points, data, coordinates) result(res)
    class(interpolator_1d), intent(in)    :: this
    sll_int32, intent(in)                 :: num_points    ! number of interpolation points
    sll_real64, dimension(:), intent(in)  :: data          ! data to be interpolated 
    sll_real64, dimension(:), intent(in)  :: coordinates   ! points where output is desired
    sll_real64, dimension(num_points)     :: res

    res = data
    print*, 'sll_interpolators_interface.F90: This is the dummy interpolator_1d function. It should never be called'
    stop
  end function interpolate_1d
  ! Reconstruction consists in computing a flux function from cell average data for a conservative scheme
  ! We provide here a reconstruct_1d function that computes the primitive and interpolates the primitive
  function reconstruct_1d(this, num_points, data, coordinates) result(res)
    class(interpolator_1d), intent(in)    :: this
    sll_int32, intent(in)                 :: num_points    ! number of interpolation points
    sll_real64, dimension(:), intent(in)  :: data          ! data to be interpolated 
    sll_real64, dimension(:), intent(in)  :: coordinates   ! points where output is desired
    sll_real64, dimension(num_points)     :: res
    ! local variables
    sll_real64, dimension(num_points) :: primitive
    sll_real64  :: delta
    sll_real64  :: eta
    sll_real64  :: avg
    sll_int32   :: i
    sll_int32   :: ierr
    sll_real64  :: coord_new, lperiod

    ! compute the primitive on the logical mesh on [0,1] assuming that data are cell values
    !--------------------------------------------------------------------------------------
    delta = 1.0_f64 / (num_points-1)
    primitive (1) = 0.0
    do i = 2, num_points
       primitive ( i ) = primitive ( i-1 ) + delta * data( i )
    end do
    ! average of dist func along the line
    avg = primitive ( num_points )
    ! modify primitive so that it becomes periodic
    eta = 0.0
    do i = 2, num_points
       eta = eta + delta
       primitive ( i ) = primitive ( i ) - avg * eta
    end do
    ! Interpolate primitive at interpolation points
    !----------------------------------------------
    res = this%interpolate_1d(num_points, primitive, coordinates)
    ! Add back value that had been removed for periodicity
    !-----------------------------------------------------
     if (coordinates(1) > 0.5_f64) then
        lperiod = -1.0_f64
     else
        lperiod = 0.0_f64
     end if
     coord_new = coordinates(1) +  lperiod
     res ( 1 ) = res ( 1 ) + avg * coord_new
     do i = 2, num_points
        ! We need here to find the points where it has been modified by periodicity
        if (coordinates(i) < coordinates(i-1)) then
           lperiod = lperiod + 1.0_f64
        end if
        coord_new = coordinates(i) +  lperiod
        res ( i ) = res ( i ) + avg * coord_new
     end do
     
   end function reconstruct_1d

 end module sll_interpolator_1d
