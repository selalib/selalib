program test_time
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
  use sll_timer
  use sll_fft

implicit none
sll_real64, dimension (:), allocatable :: listcoord, val, test
sll_int32 :: nbpoint
integer :: i, j, verification
type(time_mark), pointer :: t1 => NULL()
type(time_mark), pointer :: t2 => NULL()
double precision :: timer, timerforward, timeracces, timerinverse, timeraccesmacro
type(sll_fft_plan), pointer :: plan => NULL()
sll_comp64 :: mode

print*,"Nombre de point sur un coté ? :"
read*, nbpoint

!#define GET_MODE0 (mode,data) \
!      mode = cmplx(data(0),0.0_f64,kind=f64) 
#define GET_MODE0 (mode) mode=cmplx(0.0_f64,0.0_f64,kind=f64) 

#define SET_MODE0 (mode,data) \
      data(0) = real(mode,kind=f64)

!Debut du timer t1
t1 => new_time_mark()
t1 => start_time_mark(t1)

allocate (listcoord(0:nbpoint-1))
allocate (val(0:nbpoint*nbpoint-1))
allocate (test(0:nbpoint*nbpoint-1))
val=0
test=0
listcoord(0)=0
listcoord(nbpoint-1)=1
do i=1,nbpoint-2
listcoord(i)=i*1./(nbpoint-1)
end do

do i=0, nbpoint-1
 do j=0,nbpoint-1
  val(i*nbpoint+j)=func(listcoord(j),listcoord(i))
 end do
end do

mode=0

!forward = -1  inverse = 1
print*,"fft forward ..."
t2 => new_time_mark()
t2 => start_time_mark(t2)
plan => fftpack_new_plan_r2r_1d(nbpoint*nbpoint,val,test,-1)
call fftpack_apply_plan_r2r_1d(plan,val,test)

timerforward = time_elapsed_since(t2)
t2 => reset_time_mark(t2)

!Acces mod with function
do i=0,nbpoint*nbpoint-1
mode = fftpack_get_mode_real_1d(plan,test,i)
mode = mode*2
call fftpack_set_mode_real_1d(plan,test,mode,i)
end do
do i=0,nbpoint*nbpoint-1
mode = fftpack_get_mode_real_1d(plan,test,i)
mode = mode/2
call fftpack_set_mode_real_1d(plan,test,mode,i)
end do 

timeracces = time_elapsed_since(t2)
t2 => reset_time_mark(t2)

!Acces mod with macro
GET_MODE0(mode)
mode=mode*2 
SET_MODE0(mode,test)
do i=1,nbpoint*nbpoint-1
mode = fftpack_get_mode_real_1d(plan,test,i)
mode = mode*2
call fftpack_set_mode_real_1d(plan,test,mode,i)
end do
GET_MODE0(mode,test)
mode=mode/2
SET_MODE0(mode,test)
do i=1,nbpoint*nbpoint-1
mode = fftpack_get_mode_real_1d(plan,test,i)
mode = mode/2
call fftpack_set_mode_real_1d(plan,test,mode,i)
end do

timeraccesmacro = time_elapsed_since(t2)
t2 => reset_time_mark(t2)

print*,"fft inverse ..."
plan => fftpack_new_plan_r2r_1d(nbpoint*nbpoint,test,test,1,1)
call fftpack_apply_plan_r2r_1d(plan,test,test)

timerinverse = time_elapsed_since(t2)

print*,"Temps d'execution avec FFTPACK : "
print*, "fft forward 		: ", timerforward
print*, "fft inverse		: ", timerinverse
print*, "acces mod		: ", timeracces
print*, "acces mod avec macro	: ", timeraccesmacro

verification = 1
do i=0,nbpoint*nbpoint-1
  if(abs(val(i)-test(i))>0.00001) then
!    print *, "val, test, ",val(i),test(i)
    verification=0
  end if
end do

if (verification==0) then
 print*,"FAIL"

 print*,""
 print*,"val"
 print*,val(:)
 print*,"test"
 print*,test(:)

else
 print*,"TEST OK"
end if

timer = time_elapsed_since(t1)
print*, "Temps d'execution du programme : ", timer

contains

function func(x,y)
sll_real64 :: x,y,func

func=16*x*y*(x-1)*(y-1)
end function func

  function fftpack_get_mode_real_1d(plan,data,k) result(mode)
    type(sll_fft_plan), pointer :: plan
    sll_real64, dimension(0:)   :: data
    sll_int32                   :: k, n_2, n
    sll_comp64                  :: mode

    n = plan%problem_shape(1)
    n_2 = n/2 !ishft(n,-1)

    if( k .eq. 0 ) then
      mode = cmplx(data(0),0.0_f64,kind=f64)
    else if( k .eq. n_2 ) then
      mode = cmplx(data(n-1),0.0_f64,kind=f64)
    else if( k .gt. n_2 ) then
      mode = cmplx( data(2*(n-k)-1) , -data(2*(n-k)) ,kind=f64)
    else
      mode = cmplx( data(2*k-1) , data(2*k) ,kind=f64)
    endif
  end function

  subroutine fftpack_set_mode_real_1d(plan,data,new_value,k)
    type(sll_fft_plan), pointer :: plan
    sll_real64, dimension(0:)   :: data
    sll_int32                   :: k, n_2, n, index_mode
    sll_comp64                  :: new_value

    n = plan%problem_shape(1)
    n_2 = n/2 !ishft(n,-1)

    if( k .eq. 0 ) then
      data(0) = real(new_value,kind=f64)
    else if( k .eq. n_2 ) then
      data(n-1) = real(new_value,kind=f64)
    else if( k .gt. n_2 ) then
      data(2*(n-k)-1) = real(new_value,kind=f64)
      data(2*(n-k)) = -dimag(new_value)
    else
      data(2*k-1) = real(new_value,kind=f64)
      data(2*k) = dimag(new_value)
    endif
  end subroutine


 
end program
