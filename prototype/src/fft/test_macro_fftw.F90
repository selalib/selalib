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
double precision :: timer, timerforward, timeracces, timerinverse
type(sll_fft_plan), pointer :: plan => NULL()
sll_comp64 :: mode

print*,"Nombre de point sur un coté ? :"
read*, nbpoint

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

!forward = -1  inverse = 1
print*,"fft forward ..."
t2 => new_time_mark()
t2 => start_time_mark(t2)
plan => fftw_new_plan_r2r_1d(nbpoint*nbpoint,val,test,-1)
call fftw_apply_plan_r2r_1d(plan,val,test)

timerforward = time_elapsed_since(t2)
t2 => reset_time_mark(t2)

do i=0,nbpoint*nbpoint-1
mode = fftw_get_mode_real_1d(plan,test,i)
mode = mode*2
call fftw_set_mode_real_1d(plan,test,mode,i)
end do 

timeracces = time_elapsed_since(t2)
t2 => reset_time_mark(t2)

print*,"fft inverse ..."
plan => fftw_new_plan_r2r_1d(nbpoint*nbpoint,test,test,1,1)
call fftw_apply_plan_r2r_1d(plan,test,test)

timerinverse = time_elapsed_since(t2)

print*,"Temps d'execution avec SELALIB : "
print*, "fft forward : ", timerforward
print*, "acces mod   : ", timeracces
print*, "fft inverse : ", timerinverse

verification = 1
do i=0,nbpoint*nbpoint-1
  if(abs(val(i)-test(i))>0.00001) then
!    print *, "val, test, ",val(i),test(i)
    verification=0
  end if
end do

if (verification==0) then
 print*,"FAIL"

! print*,""
! print*,"val"
! print*,val(:)
! print*,"test"
! print*,test(:)

else
 print*,"TEST OK"
end if

deallocate(test)
deallocate(val)
deallocate(listcoord)

timer = time_elapsed_since(t1)
print*, "Temps d'execution du programme : ", timer

contains

function func(x,y)
sll_real64 :: x,y,func

func=16*x*y*(x-1)*(y-1)
end function func

  function fftw_get_mode_real_1d(plan,data,k) result(mode)
    type(sll_fft_plan), pointer :: plan
    sll_real64, dimension(0:)   :: data
    sll_int32                   :: k, n_2, n
    sll_comp64                  :: mode

    n = plan%problem_shape(1)
    n_2 = n/2 !ishft(n,-1)

      if( k .eq. 0 ) then
        mode = cmplx(data(0),0.0_f64,kind=f64)
      else if( k .eq. n_2 ) then
        mode = cmplx(data(n_2),0.0_f64,kind=f64)
      else if( k .gt. n_2 ) then
        !mode = complex( data(k-n_2) , -data(n-k+n_2) )
        mode = cmplx( data(n-k) , -data(k) ,kind=f64)
      else
        mode = cmplx( data(k) , data(n-k) ,kind=f64)
      endif
  end function

  subroutine fftw_set_mode_real_1d(plan,data,new_value,k)
    type(sll_fft_plan), pointer :: plan
    sll_real64, dimension(0:)   :: data
    sll_int32                   :: k, n_2, n, index_mode
    sll_comp64                  :: new_value

    n = plan%problem_shape(1)
    n_2 = n/2 !ishft(n,-1)

      if( k .eq. 0 ) then
        data(0) = real(new_value,kind=f64)
      else if( k .eq. n_2 ) then
        data(n_2) = real(new_value,kind=f64)
      else if( k .gt. n_2 ) then
        data(n-k) = real(new_value,kind=f64)
        data(k) = -dimag(new_value)
      else
        data(k) = real(new_value,kind=f64)
        data(n-k) = dimag(new_value)
      endif
  end subroutine 
end program
