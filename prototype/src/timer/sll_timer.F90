!> @file sll_timer.F90
!> @namespace sll_timer
!> @author samuel de santis samuel.desantis@math.unistra.fr
!> @brief Timing facility
!> @details <b> How to use-it </b>
!>
!> First, add the line
!> @code use sll_timer @endcode
!> declare a new time_mark pointer
!> @code type(time_mark), pointer :: mark @endcode
!> and allocate it
!> @code mark => new_time_mark() @endcode
!> @warning note that we use "=>" and not just "=".
!> \a mark is a pointer!
!>
!> We can now use the functions
!> @code
!> real :: time                            (you can use 64-bits real)
!> mark => start_time_mark(mark)
!> time = time_elapsed_between(mark1,mark2)
!> time = time_elapsed_since(mark)
!> @endcode
!> At the ending, you must free the mark with
!> @code mark => delete_mark_time(mark) @endcode
!> 
module sll_timer
  implicit none

  !> @brief type use for clock reading
  type time_mark
    integer(kind=selected_int_kind(10)) :: clocks_ticks !< number of clocks ticks
    integer(kind=selected_int_kind(10)) :: clocks_ticks_max !< the maximum value that clocks_ticks can reach
    integer(kind=selected_int_kind(10)) :: clocks_ticks_per_sec !< the number of cycles per second
  end type time_mark
 
  contains

  !> @brief allocate new time mark
  !> @return a pointer to a new time mark
  function new_time_mark()
      type(time_mark), pointer :: new_time_mark
      allocate(new_time_mark)
      CALL SYSTEM_CLOCK(COUNT_RATE=new_time_mark%clocks_ticks_per_sec,&
                         COUNT_MAX=new_time_mark%clocks_ticks_max)
  end function new_time_mark
 
  !> @brief Start the timer
  !> @param mark a pointer to the mark to initilize
  !> @return a pointer to the mark
  function start_time_mark(mark)
    type(time_mark), pointer :: start_time_mark
    type(time_mark), pointer :: mark
    if( .not. associated(mark) ) then
        print *, 'module SLL_TIMER : function START_TIME_MARK'
        print *, 'argument must be a pointer associated to time_mark object'
    endif
    call system_clock(count=mark%clocks_ticks)
    start_time_mark => mark
  end function start_time_mark

  !> @brief Reset the timer
  !> @param mark a pointer to the mark to reset
  !> @return a pointer to the mark
  function reset_time_mark(mark)
    type(time_mark), pointer :: reset_time_mark
    type(time_mark), pointer :: mark
    if( .not. associated(mark) ) then
        print *, 'module SLL_TIMER : function RESET_TIME_MARK'
        print *, 'argument must be a pointer associated to time_mark object'
    endif
    call system_clock(count=mark%clocks_ticks)
    reset_time_mark => mark
  end function reset_time_mark
  
  !> @brief Computes the time elapsed between two mark
  !> @param t0 a pointer to the first mark
  !> @param t1 a pointer to the second mark
  !> @return the time elapsed between t0 and t1 
  function time_elapsed_between(t0,t1)
    real(kind=kind(0d0)) :: time_elapsed_between
    type(time_mark), pointer, intent(in) :: t0, t1
    integer(kind=kind(t0%clocks_ticks)) :: clocks_ticks

    if( (.not. associated(t0)) .or. (.not. associated(t1)) ) then
        print *, 'module SLL_TIMER : function TIME_ELAPSED_BETWEEN'
        print *, 'argument must be a pointer associated to time_mark object'
    endif

    clocks_ticks = t1%clocks_ticks - t0%clocks_ticks
    if (t1%clocks_ticks < t0%clocks_ticks) &
         clocks_ticks = clocks_ticks + t0%clocks_ticks_max
    time_elapsed_between = real(clocks_ticks,kind=kind(time_elapsed_between)) / t0%clocks_ticks_per_sec
  end function time_elapsed_between

  !> @brief Computes the time elapsed since one mark
  !> @param t0 a pointer to the mark
  !> @return the time elapsed since t0
  function time_elapsed_since(t0)
    real(kind=kind(0d0)) :: time_elapsed_since
    type(time_mark), pointer :: t0
    type(time_mark) :: t1
    integer(kind=kind(t0%clocks_ticks)) :: clocks_ticks

    if( .not. associated(t0) ) then
        print *, 'module SLL_TIMER : function TIME_ELAPSED_SINCE'
        print *, 'argument must be a pointer associated to time_mark object'
    endif

    call system_clock(count=t1%clocks_ticks)

    clocks_ticks = t1%clocks_ticks - t0%clocks_ticks
    if (t1%clocks_ticks < t0%clocks_ticks) &
         clocks_ticks = clocks_ticks + t0%clocks_ticks_max
    time_elapsed_since = real(clocks_ticks,kind=kind(time_elapsed_since))&
                                               / t0%clocks_ticks_per_sec
  end function time_elapsed_since

  !> @brief deletes the mark in argument
  !> @param t0 a pointer to the mark to delete
  !> @return a null pointer 
  function delete_time_mark(t0)
    type(time_mark), pointer :: delete_time_mark
    type(time_mark), pointer :: t0

    deallocate(t0)
    t0 => NULL()
    delete_time_mark => NULL()
  end function delete_time_mark

end module sll_timer
