program test_time
#include "sll_working_precision.h"
#include "sll_fft.h"
 use sll_timer

implicit none
sll_real64, dimension(:), allocatable ::valeur, data_in, data_out,test,test2
sll_int32 :: nbpoint,n,n_2
integer :: i,j,verification
type(sll_fft_plan), pointer :: plan => NULL()
type(time_mark), pointer :: t1 => NULL()
type(time_mark), pointer :: t2 => NULL()
type(time_mark), pointer :: t3 => NULL()
type(time_mark), pointer :: t4 => NULL()
type(time_mark), pointer :: t5 => NULL()
type(time_mark), pointer :: t6 => NULL()
double precision :: timerforward,timerinverse,timeraccess
double precision :: timeraccessmacro,timertotal,timertotalmacro
sll_comp64 :: mode


!  #define GET_MODE0(mode,data) \
!        mode = cmplx(data(0),0.0_f64,kind=f64)
!  #define SET_MODE0(new_value,data) \
!        data(0) = real(new_value,kind=f64)
!  
!  if(_DEFAULTFFTLIB==0) then
!   print*,"librairie fft : SLLFFT"
!  #define GET_MODE_N_2(mode,data) \
!          mode = cmplx(data(1),0.0_f64,kind=f64)
!  #define GET_MODE_GT_N_2(mode,data,k) \
!          mode = cmplx( data(2*(n-k)) , -data(2*(n-k)+1),kind=f64)
!  #define GET_MODE_LT_N_2(mode,data,k) \        
!          mode = cmplx( data(2*k) , data(2*k+1) ,kind=f64)
!  
!  #define SET_MODE_N_2(new_value,data) \
!          data(1) = real(new_value,kind=f64)
!  #define SET_MODE_GT_N_2(new_value,data,k) \
!          data(2*(n-k)) = real(new_value,kind=f64); \
!          data(2*(n-k)+1) = -dimag(new_value)
!  #define SET_MODE_LT_N_2(new_value,data,k) \
!          data(2*k) = real(new_value,kind=f64); \
!          data(2*k+1) = dimag(new_value)
!  
!  elseif(_DEFAULTFFTLIB==100) then
!   print*,"librairie fft : FFTPACK"
!  !#define GET_MODE_N_2(mode,data) \
!  !      mode = cmplx(data(n-1),0.0_f64,kind=f64)
!  !#define GET_MODE_GT_N_2(mode,data,k) \
!  !      mode = cmplx( data(2*(n-k)-1) , -data(2*(n-k)) ,kind=f64)
!  !#define GET_MODE_LT_N_2(mode,data,k) \   
!  !      mode = cmplx( data(2*k-1) , data(2*k) ,kind=f64)     
!  !
!  !#define SET_MODE_N_2(new_value,data) \
!  !      data(n-1) = real(new_value,kind=f64)
!  !#define SET_MODE_GT_N_2(new_value,data,k) \
!  !      data(2*(n-k)-1) = real(new_value,kind=f64); \
!  !      data(2*(n-k)) = -dimag(new_value)
!  !#define SET_MODE_LT_N_2(new_value,data,k) \
!  !      data(2*k-1) = real(new_value,kind=f64); \
!  !      data(2*k) = dimag(new_value)
!  
!  elseif(_DEFAULTFFTLIB==1000000000)then
!   print*,"librairie fft : FFTW"
!  !#define GET_MODE_N_2(mode,data) \
!  !        mode = cmplx(data(n_2),0.0_f64,kind=f64)
!  !#define GET_MODE_GT_N_2(mode,data,k) \
!  !        mode = cmplx( data(n-k) , -data(k) ,kind=f64)
!  !#define GET_MODE_LT_N_2(mode,data,k) \ 
!  !        mode = cmplx( data(k) , data(n-k) ,kind=f64)  
!  !
!  !#define SET_MODE_N_2(new_value,data) \
!  !        data(n_2) = real(new_value,kind=f64)
!  !#define SET_MODE_GT_N_2(new_value,data,k) \
!  !        data(n-k) = real(new_value,kind=f64); \
!  !        data(k) = -dimag(new_value)
!  !#define SET_MODE_LT_N_2(new_value,data,k) \
!  !        data(k) = real(new_value,kind=f64); \
!  !        data(n-k) = dimag(new_value)
!  
!  else
!   print*,"librairie fft non reconnue"
!   print*,"arret du test"
!   STOP
!  endif

TEST

if(_DEFAULTFFTLIB==0) then
 print*,"Nombre de points sur un coté ? (puissance de 2) :"
else
 print*,"Nombre de points sur un coté ? :"
endif
read*,nbpoint

allocate (valeur(1:nbpoint**2))
allocate (test(1:nbpoint**2))
allocate (test2(1:nbpoint**2))
allocate (data_in(0:nbpoint-1))
allocate (data_out(0:nbpoint-1))


do i=1,nbpoint
 do j=1,nbpoint
  valeur((i-1)*nbpoint+j)=func((i-1)*1./(nbpoint-1),(j-1)*1./(nbpoint-1))
 end do
end do

!initialisation des timers
timerforward=0
timerinverse=0
timeraccess=0
timeraccessmacro=0
timertotal=0
timertotalmacro=0
t1 => new_time_mark()
t1 => start_time_mark(t1)
t2 => new_time_mark()
t2 => start_time_mark(t2)
t3 => new_time_mark()
t3 => start_time_mark(t3)
t4 => new_time_mark()
t4 => start_time_mark(t4)
t5 => new_time_mark()
t5 => start_time_mark(t5)
t6 => new_time_mark()
t6 => start_time_mark(t6)

!forward = -1, inverse = 1
!!!Premiere partie : access mode
t5 => reset_time_mark(t5)
do j=1,nbpoint
 do i=1,nbpoint
  data_in(i-1)=valeur((j-1)*nbpoint+i)
 end do
 !fft forward
 plan => fft_new_plan(nbpoint,data_in,data_out,-1)
 t1 => reset_time_mark(t1)
 call fft_apply_plan(plan,data_in,data_out)
 timerforward = timerforward+time_elapsed_since(t1)

 !access mode SANS macro
 t3 => reset_time_mark(t3)

 do i=0,nbpoint-1
  mode = fft_get_mode(plan,data_out,i)
  mode = mode*2
  call fft_set_mode(plan,data_out,mode,i)
 end do
 do i=0,nbpoint-1
  mode = fft_get_mode(plan,data_out,i)
  mode = mode/2
  call fft_set_mode(plan,data_out,mode,i)
 end do

 timeraccess=timeraccess+time_elapsed_since(t3)

 !fft inverse
 plan => fft_new_plan(nbpoint,data_out,data_out,1,1)
 t2 => reset_time_mark(t2)
 call fft_apply_plan(plan,data_out,data_out)
 timerinverse = timerinverse + time_elapsed_since(t2)

! do i=0,nbpoint-1
!  test((j-1)*nbpoint+i+1)=data_out(i)
! end do 
end do
timertotal = time_elapsed_since(t5)

!!!!Deuxieme partie : access mode with macro
t6 => reset_time_mark(t6)
do j=1,nbpoint
 do i=1,nbpoint
  data_in(i-1)=valeur((j-1)*nbpoint+i)
 end do
 !fft forward
 plan => fft_new_plan(nbpoint,data_in,data_out,-1)
 t1 => reset_time_mark(t1)
 call fft_apply_plan(plan,data_in,data_out)
 timerforward = timerforward+time_elapsed_since(t1)

 !access mode AVEC macro
 t4 => reset_time_mark(t4)
 n = plan%problem_shape(1)
 n_2 = n/2 !ishft(n,-1)

 GET_MODE0(mode,data_out)
 mode = mode*2
 SET_MODE0(mode,data_out)
 GET_MODE_N_2(mode,data_out)
 mode=mode*2
 GET_MODE_N_2(mode,data_out)
 do i=1,n_2-1 
  GET_MODE_LT_N_2(mode,data_out,i)
  mode=mode*2
  SET_MODE_LT_N_2(mode,data_out,i)
 end do
 do i=n_2+1,nbpoint-1
  GET_MODE_GT_N_2(mode,data_out,i)
  mode=mode*2
  SET_MODE_GT_N_2(mode,data_out,i)
 end do

 GET_MODE0(mode,data_out)
 mode = mode*0.5
 SET_MODE0(mode,data_out)
 GET_MODE_N_2(mode,data_out)
 mode=mode*0.5
 GET_MODE_N_2(mode,data_out)
 do i=1,n_2-1 
  GET_MODE_LT_N_2(mode,data_out,i)
  mode=mode*0.5
  SET_MODE_LT_N_2(mode,data_out,i)
 end do
 do i=n_2+1,nbpoint-1
  GET_MODE_GT_N_2(mode,data_out,i)
  mode=mode*0.5
  SET_MODE_GT_N_2(mode,data_out,i)
 end do
 timeraccessmacro=timeraccessmacro+time_elapsed_since(t4)

 !fft inverse
 plan => fft_new_plan(nbpoint,data_out,data_out,1,1)
 t2 => reset_time_mark(t2)
 call fft_apply_plan(plan,data_out,data_out)
 timerinverse = timerinverse + time_elapsed_since(t2)

! do i=0,nbpoint-1
!  test2((j-1)*nbpoint+i+1)=data_out(i)
! end do 
end do
timertotalmacro = time_elapsed_since(t6)

!verification=1
!i=1
!do while(i<nbpoint*nbpoint+1.and.verification==1)
! if(abs(valeur(i)-test(i))>0.0001.or.abs(valeur(i)-test2(i))>0.0001) then
!  verification =0
! end if
! i=i+1
!end do
!
!print*,""
!if(verification == 0) then
! print*,"Test fail"
!else
! print*,"Test OK"
!end if

print*,""
print*,"Temps d'execution  :"
!print*," fft forward	    : ",timerforward*0.5
!print*," fft inverse	    : ",timerinverse*0.5
print*," access mode	    : ",timeraccess
print*," access mode macro  : ",timeraccessmacro
print*," temps total        : ",timertotal
print*," temps total macro  : ",timertotalmacro
print*,""
!print*," difference entre les temps : ", timertotal-timertotalmacro
!print*," difference entre les access mode : ", timeraccess-timeraccessmacro
!print*,""
print*," rapport des temps : ", timertotal/timertotalmacro
print*," rapport des temps des access : ", timeraccess/timeraccessmacro

deallocate(data_out)
deallocate(data_in)
deallocate(test)
deallocate(test2)
deallocate(valeur)
contains

function func(x,y)
sll_real64 func
real :: x,y
func = 16*x*y*(x-1)*(y-1)
end function func

end program
