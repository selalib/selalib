program unit_test
  use sll_timer
  implicit none

  type(time_mark), pointer :: t0 => NULL()
  type(time_mark), pointer :: t1 => NULL()
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

  t0 => new_time_mark()
  t1 => new_time_mark()

  print *, 'Start time mark t0'
  t0 => start_time_mark(t0)
  print *, '  ..sleep 2 seconds..  '
  call sleep(2)
  time = time_elapsed_since(t0)
  print *, 'time elapsed since t0 : ',time
  print *, 'reset time mark t0'

  a = int(time)  

  print *, 'Start time mark t1'
  t0 => reset_time_mark(t0)
  t1 => start_time_mark(t1)
  time = time_elapsed_between(t0,t1)
  print *, 'Time between t0 and t1',time 
  print *,'...'

  t0 => delete_time_mark(t0)
  t0 => delete_time_mark(t1)

  if(a .eq. 2) then
    print *, 'PASS'
  else
    print *, 'NOT PASS'
  endif
end program unit_test
