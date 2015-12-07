!**************************************************************
!  Author: Jakob Ameres, jakob.ameres@tum.de
!**************************************************************

program pif_fieldsolver_test
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"


use sll_m_pif_fieldsolver
use sll_m_timer
use sll_m_utilities, only : &
     display_matrix_2d_integer
implicit none


    
!sll_real64 ::time
sll_int32 :: ierr
sll_int32 :: npart !number of particles
sll_real64, dimension(:,:), allocatable :: x!,y !particle coordinate
 type(sll_time_mark)  :: tstart, tstop
!sll_int32 :: maxmode=10
!sll_comp64, dimension(:), allocatable :: fmodes
!sll_comp64, dimension(:), allocatable :: modeone
!sll_int32, dimension(:,:), allocatable :: modes
sll_int32 :: idx, dimx
sll_real64 :: boxlen=2*sll_pi
type(pif_fieldsolver) :: SOLVER

sll_int32, dimension(:),allocatable :: stenx,stenv,stenxw
sll_int32 :: stenw


dimx=2;

!allocate stencils
SLL_ALLOCATE(stenx(dimx),ierr)
SLL_ALLOCATE(stenxw(dimx+1),ierr)
SLL_ALLOCATE(stenv(dimx),ierr)
! stenx=0; stenxw=0; stenv=0;

stenx=(/( idx,idx=1,dimx)/)
stenv=(/( idx,idx=dimx+1,2*dimx)/)
stenw=2*dimx+1
stenxw(1:dimx)=stenx
stenxw(dimx+1)=stenw


!\(idx, idx=1,dimx\);
! stenv(dimx+1:2*dimx)=1_i32;
! stenxw=stenx
! stenxw(2*dimx+1)=1_i32 !include the weights


SOLVER%dimx=dimx
!Load some particles


call SOLVER%set_box_len(boxlen)
call SOLVER%init(15)
 call display_matrix_2d_integer(transpose(SOLVER%allmodes),'i8')


npart=int(1e5,i32)

call speed_test()

do idx=1,SOLVER%problemsize()

print *,"Testing mode: ",SOLVER%allmodes(:,idx)
call mode_test(idx)
end do

contains

subroutine mode_test(modeidx)
sll_comp64, dimension(:), allocatable :: rhs, sol
sll_real64, dimension(:,:), allocatable :: xw, testx
sll_real64, dimension(:), allocatable :: y1, y2
sll_real64, dimension(:,:), allocatable :: grad
sll_int32, intent(in) :: modeidx
sll_real64, dimension(dimx) :: mode
sll_real64 :: alpha=0.5_f64
sll_int32 :: idx!, jdx

mode=(SOLVER%allmodes(:,modeidx)*SOLVER%unitmode(:))

SLL_ALLOCATE(rhs(SOLVER%problemsize()),ierr)
SLL_ALLOCATE(sol(SOLVER%problemsize()),ierr)

SLL_CLEAR_ALLOCATE(xw(1:dimx+1, 1:npart), ierr)
call random_number(xw(1:dimx,:));
xw=xw*boxlen

xw(dimx+1,:)=0.0_f64
!  do idx=1,dimx
!   xw(dimx+1,:)=xw(dimx+1,:)+ sin(xw(idx,:)*mode(idx))
!  end do

 xw(dimx+1,:)=1.0_f64+alpha*cos(matmul(mode, xw(1:dimx,:))+boxlen/3)

!scale weight
xw(dimx+1,:)=xw(dimx+1,:)*boxlen**dimx

rhs=SOLVER%get_fourier_modes(xw)/npart

!Set up test grid (randomly)
SLL_CLEAR_ALLOCATE(testx(1:dimx, 1:10000), ierr)
call random_number(testx);
testx=testx*boxlen


!test Mass
sol=SOLVER%solve_mass(rhs)

SLL_CLEAR_ALLOCATE(y1(size(testx,2)), ierr)
SLL_CLEAR_ALLOCATE(y2(size(testx,2)), ierr)

y1=1.0_f64+alpha*cos(matmul(mode, testx(1:dimx,:))+boxlen/3)
y2=SOLVER%eval_solution(testx, sol)

print *, abs(sol(modeidx))
print *,"RHO L2:", sum((y1-y2)**2)/real(size(testx),f64)


!Test components of gradient
SLL_CLEAR_ALLOCATE(grad(3,size(testx,2)), ierr)
grad=SOLVER%eval_gradient(testx,sol)
do idx=1,dimx
y1=-alpha*sin(matmul(mode, testx(1:dimx,:))+boxlen/3)*mode(idx)
print *,"NABLA RHO L2:",idx,sum((y1-grad(idx,:))**2)/real(size(testx),f64)
end do


!Test poisson solve
sol=SOLVER%solve_poisson(rhs)

y1=alpha*cos(matmul(mode, testx(1:dimx,:))+boxlen/3)/sum((mode(:)**2))
y2=SOLVER%eval_solution(testx, sol)
print *,"PHI L2:", sum((y1-y2)**2)/real(size(testx),f64)


end subroutine 


subroutine speed_test()
sll_int32 :: npart, chunksize
sll_comp64, dimension(:), allocatable :: rhs1,rhs2,rhs3,rhs4

print *,"Test for small number of particles"
npart=int(1e5,i32)
SLL_CLEAR_ALLOCATE(x(2*dimx+1, 1:npart), ierr)
call random_number(x);
x=x*boxlen

x(stenw,:)=(1.0_f64+0.5_f64*sin(1*x(1,:)-x(2,:)))*boxlen**dimx
!+sin(2*x(1,:))+sin(2*x(3,:))

SLL_ALLOCATE(rhs1(SOLVER%problemsize()),ierr)
SLL_ALLOCATE(rhs2(SOLVER%problemsize()),ierr)
SLL_ALLOCATE(rhs3(SOLVER%problemsize()),ierr)
SLL_ALLOCATE(rhs4(SOLVER%problemsize()),ierr)

print *,"---------------------------------------------------"

chunksize=256
print *, "Chunk size is set to: ", chunksize

do idx=1,4
call sll_set_time_mark(tstart)
rhs1=SOLVER%get_fourier_modes(x(stenxw,:))/npart
call sll_set_time_mark(tstop)
print *, 'Standard          ', sll_time_elapsed_between(tstart,tstop)

call sll_set_time_mark(tstart)
rhs1=SOLVER%get_fourier_modes_chunk(x(stenxw,:),chunksize)/npart
call sll_set_time_mark(tstop)
print *, 'Standard chunked: ', sll_time_elapsed_between(tstart,tstop)


call sll_set_time_mark(tstart)
rhs2=SOLVER%get_fourier_modes2(x(stenxw,:))/npart
!print *, abs(SOLVER%get_fourier_modes2(x(stenxw,:)))/npart
call sll_set_time_mark(tstop)
print *, 'Fast2:            ', sll_time_elapsed_between(tstart,tstop)

call sll_set_time_mark(tstart)
rhs2=SOLVER%get_fourier_modes2_chunk(x(stenxw,:),chunksize)/npart
call sll_set_time_mark(tstop)
print *, 'Fast2 chunked:    ', sll_time_elapsed_between(tstart,tstop)
end do

end subroutine speed_test


!  print *, SOLVER%get_fourier_modes2_chunk(x(stenxw,:),256)

! call display_matrix_2d_real(transpose(x(stenxw, :)),'F10.2')
! call display_matrix_2d_real(transpose(x),'F10.2')

!modes=generate_exponents( (/0, -2, -1/), (/ 1 , 1, 0/))
!  call display_matrix_2d_integer(transpose(SOLVER%allmodes),'i8')




end program pif_fieldsolver_test
