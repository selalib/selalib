!**************************************************************
!  Copyright INRIA
!  Authors : 
!     CALVI project team
!  
!  This code SeLaLib (for Semi-Lagrangian-Library) 
!  is a parallel library for simulating the plasma turbulence 
!  in a tokamak.
!  
!  This software is governed by the CeCILL-B license 
!  under French law and abiding by the rules of distribution 
!  of free software.  You can  use, modify and redistribute 
!  the software under the terms of the CeCILL-B license as 
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info". 
!**************************************************************

!> @file sll_timer.F90
!> @namespace sll_timer
!> @author samuel de santis samuel.desantis@math.unistra.fr
!> @brief Timing facility
!> @details <b> How to use-it </b>
!>
!>
!> First, add the line
!> @code use sll_timer @endcode
!> declare an sll_time_mark object
!> @code type(sll_time_mark) :: mark @endcode
!> 

!> We can now use the functions
!> @code
!> real :: time                            (you can use 64-bit real)
!> call sll_set_time_mark(mark1)
!> call sll_set_time_mark(mark2)
!> time = time_elapsed_between(mark1,mark2)
!> time = time_elapsed_since(mark2)
!> @endcode
!> 
module sll_timer
  implicit none

  !> @brief type use for clock reading
  type sll_time_mark
    integer(kind=selected_int_kind(10)) :: clock_ticks !< number of clocks ticks
    integer(kind=selected_int_kind(10)) :: clock_ticks_max !< the maximum value that clock_ticks can reach
    integer(kind=selected_int_kind(10)) :: clock_ticks_per_sec !< the number of cycles per second
  end type sll_time_mark
 
  contains

  !> @brief reads time parameters from system and stores in its argument.
  !> param timer_obj an object of type sll_time_mark, intent(out).
  subroutine sll_set_time_mark(timer_obj)
    type(sll_time_mark), intent(out) :: timer_obj
    call system_clock(COUNT=timer_obj%clock_ticks, &
                      COUNT_RATE=timer_obj%clock_ticks_per_sec,&
                      COUNT_MAX=timer_obj%clock_ticks_max)
  end subroutine sll_set_time_mark
 

  !> @brief Computes the time elapsed between two time marks
  !> @param t0 an object of type sll_time_mark, intent(in)
  !> @param t1 an object of type sll_time_mark, intent(in)
  !> @return the time elapsed between t0 and t1 in seconds.
  function sll_time_elapsed_between(t0,t1)
    real(kind=kind(0d0))                :: sll_time_elapsed_between
    type(sll_time_mark), intent(in)     :: t0, t1
    integer(kind=kind(t0%clock_ticks))  :: clock_ticks

    clock_ticks = t1%clock_ticks - t0%clock_ticks
!!$    if (t1%clock_ticks < t0%clock_ticks) then
!!$       clock_ticks = clock_ticks + t0%clock_ticks_max
!!$    end if
    sll_time_elapsed_between = &
         real(clock_ticks,kind=kind(0d0))/real(t0%clock_ticks_per_sec,kind(0d0))
  end function sll_time_elapsed_between

  !> @brief Computes the time elapsed since a particular time mark was set.
  !> @param t0 reference mark of type sll_time_mark. Must have been set. intent(in)
  !> @return the time elapsed since t0 was set in seconds.
  function sll_time_elapsed_since(t0)
    real(kind=kind(0d0))                :: sll_time_elapsed_since
    type(sll_time_mark), intent(in)     :: t0
    integer(kind=selected_int_kind(10)) :: t1
    integer(kind=kind(t0%clock_ticks))  :: clock_ticks

    call system_clock(count=t1)
    clock_ticks = t1 - t0%clock_ticks
!!$    if (t1%clock_ticks < t0%clock_ticks) then
!!$       clock_ticks = clock_ticks + t0%clock_ticks_max
!!$    end if
    sll_time_elapsed_since = &
         real(clock_ticks,kind=kind(0d0))/real(t0%clock_ticks_per_sec,kind(0d0))
  end function sll_time_elapsed_since


end module sll_timer
