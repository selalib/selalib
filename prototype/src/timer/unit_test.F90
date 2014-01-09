program unit_test
  use sll_timer
  implicit none

  type(sll_time_mark)  :: t0 
  type(sll_time_mark)  :: t1 
  double precision :: time
  integer :: a,b

  call system_clock(COUNT_RATE=a,&
                    COUNT_MAX=b)
  if(b .eq. 0) then
    stop 'NO CLOCK AVAILABLE'
  else
   print *, 'CLOCK AVAILABLE'
   print *,'...'
  endif

  print *, 'Start time mark t0'
  call set_time_mark(t0)
  print *, '  ..sleep 2 seconds..  '
  call sleep(2)
  time = time_elapsed_since(t0)
  print *, 'time elapsed since t0 : ',time

  a = int(time)  

  print *, 'Start time mark t1'
  call set_time_mark(t1)
  time = time_elapsed_between(t0,t1)
  print *, 'Time between t0 and t1 (seconds)',time 
  call set_time_mark(t0)
  call set_time_mark(t1)
  time = time_elapsed_between(t0,t1)
  print *, 'Shortest measured time span (seconds): ',time

  if(a .eq. 2) then
    print *, 'PASSED'
  else
    print *, 'FAILED'
  endif
end program unit_test
