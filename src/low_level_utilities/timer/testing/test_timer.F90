program unit_test
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   use sll_m_timer, only: &
      sll_s_set_time_mark, &
      sll_f_time_elapsed_between, &
      sll_f_time_elapsed_since, &
      sll_t_time_mark

   use iso_c_binding, only: &
      c_int

   implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#define MARGIN 0.1

   interface
      integer(c_int) function c_sleep(seconds) bind(c, name='sleep')
         import
         integer(c_int), intent(in), value :: seconds
      end function c_sleep
   end interface

   type(sll_t_time_mark)  :: t0
   type(sll_t_time_mark)  :: t1
   double precision :: time
   integer :: a, b

   call system_clock(COUNT_RATE=a, COUNT_MAX=b)

   if (b .eq. 0) then
      stop 'NO CLOCK AVAILABLE'
   else
      print *, 'CLOCK AVAILABLE'
      print *, '...'
   end if

   print *, 'Start time mark t0'
   call sll_s_set_time_mark(t0)
   print *, '  ..sleep 2 seconds..  '
   call sleep(2)
   time = sll_f_time_elapsed_since(t0)
   print *, 'time elapsed since t0 : ', time

   print *, "Clock tics per second", t0%clock_ticks_per_sec

   a = int(time + MARGIN)

   print *, 'Start time mark t1'
   call sll_s_set_time_mark(t1)
   time = sll_f_time_elapsed_between(t0, t1)
   print *, 'Time between t0 and t1 (seconds)', time
   call sll_s_set_time_mark(t0)
   call sll_s_set_time_mark(t1)
   time = sll_f_time_elapsed_between(t0, t1)
   print *, 'Shortest measured time span (seconds): ', time

   if (a .eq. 2) then
      print *, 'PASSED'
   else
      print *, 'FAILED'
   end if

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine sleep(seconds)
      integer, intent(in) :: seconds
      integer :: elapsed
      elapsed = c_sleep(seconds)
   end subroutine sleep

end program unit_test
