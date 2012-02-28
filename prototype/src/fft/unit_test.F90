program unit_test
  use sll_fft
  !use sll_timer
#include "sll_working_precision.h"
  implicit none
 
  sll_real64, parameter  :: err_max = 10E-14
  sll_int32, parameter   :: n = 2**20
  sll_int32, parameter   :: hmax = 1
  sll_int32, parameter   :: imin = 10
  sll_int32, parameter   :: imax = 20
  type(sll_fft_plan), pointer :: plan
  sll_comp64, dimension(n) :: data_comp, data_copy
  sll_real64, dimension(n) :: rdata_comp, rdata_copy
  sll_real64 :: time = 0_f64
  sll_real64 :: ierr
  sll_int32 :: i,j,s,h


!  if( _DEFAULTFFTLIB == FFTW_MOD) then
!    open(1,file="time_fftw")
!  else if( _DEFAULTFFTLIB == SLLFFT_MOD) then
!    open(1,file="time_sllfft")
!  else if( _DEFAULTFFTLIB == FFTPACK_MOD) then
!    open(1,file="time_fftpack")
!  end if

  print *,'-------------------------------------------------'
  print * ,'COMPLEX TO COMPLEX'
  do i=imin,imax
   do h=1,hmax
    s = 2**i

    do j=1,s
      CALL RANDOM_COMPLEX(data_comp(j))
    enddo
    data_copy(1:s) = data_comp(1:s)
  
    plan => fft_new_plan(s,data_comp(1:s),data_comp(1:s),FFT_FORWARD)
    call fft_apply_plan(plan,data_comp(1:s),data_comp(1:s))

    time = time + fft_get_time_execution(plan)

    call fft_delete_plan(plan)

    plan => fft_new_plan(s,data_comp(1:s),data_comp(1:s),FFT_INVERSE,FFT_NORMALIZE)
    call fft_apply_plan(plan,data_comp(1:s),data_comp(1:s))
    call fft_delete_plan(plan)
    ierr = ERROR_MAX(data_comp(1:s) - data_copy(1:s))
    if( ierr > err_max ) then
      stop 'Everage error too big'
    endif
!    print * , 'Error max : ' , ierr
   enddo
!   write(1,*) time/real(hmax,kind=f64)
  enddo
  print *, 'OK'

!  close(1)


  print *,'-------------------------------------------------'
  print * ,'REAL TO REAL'
  do i=imin,imax
   do h=1,hmax
    s = 2**i

    do j=1,s
      CALL RANDOM_NUMBER(rdata_comp(j))
    enddo
    rdata_copy(1:s) = rdata_comp(1:s)
  
    plan => fft_new_plan(s,rdata_comp(1:s),rdata_comp(1:s),FFT_FORWARD)
    call fft_apply_plan(plan,rdata_comp(1:s),rdata_comp(1:s))

    time = time + fft_get_time_execution(plan)

    call fft_delete_plan(plan)

    plan => fft_new_plan(s,rdata_comp(1:s),rdata_comp(1:s),FFT_INVERSE,FFT_NORMALIZE)
    call fft_apply_plan(plan,rdata_comp(1:s),rdata_comp(1:s))
    call fft_delete_plan(plan)
    ierr = MAXVAL(rdata_comp(1:s) - rdata_copy(1:s))
    if( ierr > err_max ) then
      stop 'Everage error too big'
    endif
   enddo
  enddo
  print *, 'OK'



  print *,'-------------------------------------------------'
  print * ,'REAL TO COMPLEX and COMPLEX TO REAL '
  do i=imin,imax
   do h=1,hmax
    s = 2**i

    do j=1,s
      CALL RANDOM_NUMBER(rdata_comp(j))
    enddo
    rdata_copy(1:s) = rdata_comp(1:s)
  
    plan => fft_new_plan(s,rdata_comp(1:s),data_comp(1:s/2+1))
    call fft_apply_plan(plan,rdata_comp(1:s),data_comp(1:s/2+1))

    time = time + fft_get_time_execution(plan)

    call fft_delete_plan(plan)
    
    plan => fft_new_plan(s,data_comp(1:s/2+1),rdata_comp(1:s),FFT_NORMALIZE)
    call fft_apply_plan(plan,data_comp(1:s/2+1),rdata_comp(1:s))
    call fft_delete_plan(plan)
    ierr = MAXVAL(rdata_comp(1:s) - rdata_copy(1:s))
    if( ierr > err_max ) then
      stop 'Everage error too big'
    endif
   enddo
  enddo
  print *, 'OK'

contains

  FUNCTION ERROR_MAX(tab) RESULT(error)
    sll_comp64, DIMENSION(:) :: tab
    sll_real64 :: error

    error = MAX( MAXVAL(REAL(tab)) , MAXVAL(DIMAG(tab)) )
  END FUNCTION

  SUBROUTINE init_random_seed()
    INTEGER                            :: i, n, clock
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed

    CALL RANDOM_SEED(size = n)
    ALLOCATE(seed(n))

    CALL SYSTEM_CLOCK(COUNT=clock)

    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    CALL RANDOM_SEED(PUT = seed)

    DEALLOCATE(seed)
  END SUBROUTINE

  SUBROUTINE RANDOM_COMPLEX(c)
    COMPLEX(KIND=KIND(1.D0)) :: c
    REAL(KIND=KIND(1.D0))    :: realpart,imagpart
   
    CALL init_random_seed() 
    CALL RANDOM_NUMBER(realpart)
    CALL RANDOM_NUMBER(imagpart)
    c = COMPLEX(realpart,imagpart)
  END SUBROUTINE

end program unit_test
