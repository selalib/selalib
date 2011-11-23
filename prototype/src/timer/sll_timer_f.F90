!------------------------------------------------------------------------------
! Selalib
!------------------------------------------------------------------------------
!
! MODULE: sll_timer
!
!> @author
!> Module Author Name and Affiliation
!
! DESCRIPTION: 
!> timing facility.
!>
!> Selalib has its own timing facility, capable of resolving time
!! periods slightly shorter than 2 micro second.
!>
!> We define the type \a time_mark to store clock readings.
!> The available functions are:
!> - \a set_time_mark() is a function with no arguments thats returns
!! an instance of the \a time_mark type, set to a clock reading done at the
!! moment of the function call
!> - \a reset_time_mark(time_mark) is also a function, and has essentially
!! a redundant role to \a set_time_mark(). It will take a time mark and
!! return it updated with a new clock reading.
!> - \a time_elapsed_since(time_mark) is a function that takes a time mark
!! as an argument and returns the time, in seconds since the clock reading
!! contained in the passed argument. 
!>
!> <b> How to use it : </b> \n
!> ****************
!>
!> First, import the timer module.
!> \code #include "timer.h" \endcode
!> Second, declare the time marker.
!> \code time_mark :: tmark \endcode
!> Third, set the time marker and after the instruction block compute
!! the elapsed time.
!> \code
!! tmark = reset_time_mark(tmark)
!! !instruction
!! !block
!! dt = time_elapsed_since(tmark)
!! \endcode
!> The dt variable contains the execution time in second of the instruction block.
!>
! REVISION HISTORY:
! DD Mmm YYYY - Initial Version
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!------------------------------------------------------------------------------

module sll_timer
  use iso_c_binding
  implicit none

  !--------------------------------------------------------------------------- 
  ! \interface set_time_mark
  ! 
  ! @author 
  ! Routine Author Name and Affiliation.
  !
  ! DESCRIPTION:
  ! Set to a clock reading done at the moment of the function call.
  !
  ! REVISION HISTORY:
  ! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
  !
  ! @return an instance of \a time_mark type
  !---------------------------------------------------------------------------
  interface 
     type(c_ptr) function set_time_mark() bind(c, name='set_time_mark_C')
       use iso_c_binding
     end function set_time_mark
  end interface

  ! Note that the 'value' attribute is essential for this interface to work.
  ! If not present, Fortran will take the address of the pointers passed and
  ! a mess will be unleashed...
  interface
     type(c_ptr) function reset_time_mark( mark ) &
          bind(c, name='reset_time_mark_C')
       use iso_c_binding
       type(c_ptr), value :: mark
     end function reset_time_mark
  end interface

  !--------------------------------------------------------------------------- 
  ! \function time_elapsed_since
  ! \brief returns the time, in seconds, since the timer marker in argument.
  interface
     real(c_double) function time_elapsed_since( t0 ) &
          bind( c, name='time_elapsed_since_C')
       use iso_c_binding
       type(c_ptr), intent(in), value :: t0
     end function time_elapsed_since
  end interface

end module sll_timer
