program particle_fourier_modes_test
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_constants, only: &
    sll_p_i1, &
    sll_p_pi

  use sll_m_timer, only: &
    sll_s_set_time_mark, &
    sll_f_time_elapsed_between, &
    sll_t_time_mark

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    
!sll_real64 ::time
sll_int32 :: ierr
sll_int32 :: npart !number of particles
sll_real64, dimension(:), allocatable :: x !particle coordinate
type(sll_t_time_mark)  :: tstart, tstop
sll_int32 :: maxmode=10
sll_comp64, dimension(:), allocatable :: fmodes
!sll_comp64, dimension(:), allocatable :: modeone
!sll_int32 :: idx

!Load some particles




print *,"Test for small number of particles"
npart=10**5
SLL_CLEAR_ALLOCATE(x(1:npart), ierr)
call random_number(x);
x=x*2*sll_p_pi

do maxmode=20,22
SLL_ALLOCATE(fmodes(1:maxmode), ierr)
fmodes = (0.0_f64,0.0_f64)


call sll_s_set_time_mark(tstart)
call calc_modes_std(maxmode,fmodes,x)
call sll_s_set_time_mark(tstop)

print *, 'Standard:    ', sll_f_time_elapsed_between(tstart,tstop)

call sll_s_set_time_mark(tstart)
call calc_modes_fast(maxmode,fmodes,x)
call sll_s_set_time_mark(tstop)

print *, 'Fast:        ', sll_f_time_elapsed_between(tstart,tstop)

call sll_s_set_time_mark(tstart)
call calc_modes_fast2(maxmode,fmodes,x)
call sll_s_set_time_mark(tstop)

print *, 'Fast2:       ', sll_f_time_elapsed_between(tstart,tstop)

call sll_s_set_time_mark(tstart)
call calc_modes_fast2_chunked(maxmode,fmodes,x,500)
call sll_s_set_time_mark(tstop)

print *, 'Fast2_chunk: ', sll_f_time_elapsed_between(tstart,tstop)


call sll_s_set_time_mark(tstart)
call calc_modes_fast_chunked(maxmode,fmodes,x,500)
call sll_s_set_time_mark(tstop)

print *, 'Fast_chunk: ', sll_f_time_elapsed_between(tstart,tstop)



SLL_DEALLOCATE_ARRAY(fmodes,ierr)
enddo
SLL_DEALLOCATE_ARRAY(x,ierr)


!-------------------------------
! 
! print *,"Test for big number of particles"
! npart=1e6
! npart=npart*128
! SLL_CLEAR_ALLOCATE(x(1:npart), ierr)
! call random_number(x);
! x=x*2*sll_p_pi


! do maxmode=1,100
! SLL_CLEAR_ALLOCATE(fmodes(1:maxmode), ierr)
! 
! call sll_s_set_time_mark(tstart)
! call calc_modes_std(maxmode,fmodes,x)
! call sll_s_set_time_mark(tstop)
! 
! print *, 'Standard:    ', sll_f_time_elapsed_between(tstart,tstop)
! call sll_s_set_time_mark(tstart)
! call calc_modes_fast2_chunked(maxmode,fmodes,x,128)
! call sll_s_set_time_mark(tstop)
! 
! print *, 'Fast2_chunk: ', sll_f_time_elapsed_between(tstart,tstop)
! 
! SLL_DEALLOCATE_ARRAY(fmodes,ierr)
! enddo

! do fmode=1,this%num_modes
!             !Be careful here, the dot_product tends to complex conjugate stuff
!             !which we don't want in this case
!             !rhs(fmode)=dot_product(exp(-sll_p_i1*fmode*ppos*2.0_f64*sll_p_pi/this%Ilength), pweight )
!             rhs(fmode)=sum(exp(-sll_p_i1*fmode*ppos*sll_p_kx/this%Ilength)*pweight)

 contains

subroutine calc_modes_std(maxmode, fmodes, x)
sll_comp64,intent(out), dimension(:):: fmodes
sll_int32 , intent(in):: maxmode
sll_real64, dimension(:), intent(in) :: x
sll_int32 :: fmode
do fmode=1,maxmode
     fmodes(fmode)=sum(exp(-sll_p_i1*fmode*x))
enddo
end subroutine calc_modes_std

subroutine calc_modes_fast(maxmode, fmodes, x)
sll_comp64,intent(out), dimension(:):: fmodes
sll_int32 , intent(in):: maxmode
sll_real64, dimension(:), intent(in) :: x
sll_comp64, dimension(size(x)) :: modeone
sll_int32 :: fmode
modeone=exp(-sll_p_i1*x)
do fmode=1,maxmode
     fmodes(fmode)=sum(modeone**fmode)
enddo
end subroutine calc_modes_fast            


subroutine calc_modes_fast_chunked(maxmode, fmodes, x,chunksize)
sll_comp64,intent(out), dimension(:):: fmodes
sll_int32 , intent(in):: maxmode
sll_real64, dimension(:), intent(in) :: x
sll_int32, intent(in):: chunksize
sll_int32 :: numx, chunk
sll_comp64, dimension(maxmode):: fmodes_chunk

numx=size(x)
do chunk=1,numx/chunksize
  call calc_modes_fast(maxmode,fmodes_chunk,x((chunk-1)*chunksize+1:chunk*chunksize))
  fmodes=fmodes+fmodes_chunk;
enddo
end subroutine calc_modes_fast_chunked



subroutine calc_modes_fast2(maxmode, fmodes, x)
sll_comp64,intent(out), dimension(:):: fmodes
sll_int32 , intent(in):: maxmode
sll_real64, dimension(:), intent(in) :: x
sll_comp64, dimension(size(x)) :: modeone
sll_comp64, dimension(size(x)) :: mode
sll_int32 :: fmode
modeone=exp(-sll_p_i1*x)
mode=modeone
do fmode=1,maxmode
      mode=mode*modeone;
     fmodes(fmode)=sum(mode)
enddo
end subroutine calc_modes_fast2            

subroutine calc_modes_fast2_chunked(maxmode, fmodes, x,chunksize)
sll_comp64,intent(out), dimension(:):: fmodes
sll_int32 , intent(in):: maxmode
sll_real64, dimension(:), intent(in) :: x
sll_int32, intent(in):: chunksize
sll_int32 :: numx, chunk
sll_comp64, dimension(maxmode):: fmodes_chunk

numx=size(x)
do chunk=1,numx/chunksize
  call calc_modes_fast2(maxmode,fmodes_chunk,x((chunk-1)*chunksize+1:chunk*chunksize))
  fmodes=fmodes+fmodes_chunk;
enddo
end subroutine calc_modes_fast2_chunked
! subroutine 
            
            
end program particle_fourier_modes_test
 
