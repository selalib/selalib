program unit_test
#include "sll_working_precision.h"
  use sll_fft
  implicit none

#define SLLFFT_MOD 0
#define FFTPACK_MOD 100
#define FFTW_MOD 1000000000

  type(sll_fft_plan), pointer :: p => null()

  sll_real64, parameter  :: err_max = 10E-14
  sll_int32, parameter   :: hmax = 1
  sll_int32, parameter   :: imin = 10
  sll_int32, parameter   :: imax = 15
  sll_int32, parameter   :: n = 2**imax
  sll_int32, parameter   :: m1 = 2**6
  sll_int32, parameter   :: m2 = 2**4
  sll_comp64, dimension(n) :: data_comp, data_copy
  sll_comp64, dimension(m1,m2) :: data_comp2d
  sll_real64, dimension(m1,m2) :: data_real2d
  sll_comp64, dimension(m1,m2) :: data_copy2d
  sll_real64, dimension(m1,m2) :: rdata_copy2d
  sll_real64, dimension(n) :: rdata_comp, rdata_copy, rdata
  sll_comp64 :: mode
  sll_real64 :: ierr  ! this is not used as an integer below, bad name
  sll_real64 :: err_var
  sll_int32 :: i,j,s,h,k,t, array_position, ind_mode

  call print_defaultfftlib()

! test getter and setter functions in complex case
  s = 2**imin
  do j=1,s
    CALL RANDOM_COMPLEX(data_comp(j))
  enddo

  data_copy(1:s) = data_comp(1:s)
  p => fft_new_plan(s,data_comp(1:s),data_comp(1:s),FFT_FORWARD)
  do i=0,s-1
    mode = fft_get_mode(p,data_comp(1:s),i)
    call fft_set_mode(p,data_comp(1:s),mode,i)
  enddo
  ierr = ERROR_MAX(data_comp(1:s) - data_copy(1:s))

  if( ierr .ne. 0_f64 ) then
    print *,'Average error too big',ierr
    stop
  else
    print *,'get and set mode complex ok'
  endif
  call fft_delete_plan(p)

! test getter and setter functions in real case
  s = 2**imax
  do j=1,s
    CALL RANDOM_NUMBER(rdata(j))
  enddo
  rdata_copy(1:s) = rdata(1:s)
  p => fft_new_plan(s,rdata(1:s),rdata(1:s),FFT_FORWARD)
  do i=0,s/2
    mode = fft_get_mode(p,rdata(1:s),i)
    call fft_set_mode(p,rdata(1:s),mode,i)
  enddo
  ierr = MAXVAL(ABS(rdata(1:s) - rdata_copy(1:s)))
  if( ierr .ne. 0_f64 ) then
    print *,'Average error too big',ierr
    stop
  else
    print *,'get and set mode real ok'
  endif
  call fft_delete_plan(p)


! Standard do-loop on the mode
  s = 2**imin
  do j=1,s
    CALL RANDOM_COMPLEX(data_comp(j))
  enddo
  data_copy(1:s) = data_comp(1:s)
  p => fft_new_plan(s,data_comp(1:s),data_comp(1:s),FFT_FORWARD)
  do ind_mode=0,s-1
    mode = fft_get_mode(p,data_comp(1:s),ind_mode)
  enddo
  call fft_delete_plan(p)

! optimized do-loop on the mode
  s = 2**imin
  do j=1,s
    CALL RANDOM_COMPLEX(data_comp(j))
  enddo
  data_copy(1:s) = data_comp(1:s)
  p => fft_new_plan(s,data_comp(1:s),data_comp(1:s),FFT_FORWARD)
  do array_position=0,s-1
    ind_mode = fft_ith_stored_mode(p,array_position) !The only change with the standard
                                                     !do-loop is this line.
    mode = fft_get_mode(p,data_comp(1:s),ind_mode)
  enddo
  call fft_delete_plan(p)

  print *,'-------------------------------------------------'
  print * ,'COMPLEX TO COMPLEX'
  do i=imin,imax
   do h=1,hmax
    s = 2**i

    do j=1,s
      CALL RANDOM_COMPLEX(data_comp(j))
    enddo
    data_copy(1:s) = data_comp(1:s)
  
    p => fft_new_plan(s,data_comp(1:s),data_comp(1:s),FFT_FORWARD)
    call fft_apply_plan(p,data_comp(1:s),data_comp(1:s))
    call fft_delete_plan(p)

    p => fft_new_plan(s,data_comp(1:s),data_comp(1:s),FFT_INVERSE,FFT_NORMALIZE)
    call fft_apply_plan(p,data_comp(1:s),data_comp(1:s))
    call fft_delete_plan(p)
    ierr = ERROR_MAX(data_comp(1:s) - data_copy(1:s))
    if( ierr > err_max ) then
      stop 'Average error too big'
    endif
   enddo
  enddo
  print *, 'OK'


  print *,'-------------------------------------------------'
  print * ,'REAL TO REAL'
  do i=imin,imax
   do h=1,hmax
    s = 2**i

    do j=1,s
      CALL RANDOM_NUMBER(rdata(j))
    enddo
    rdata_copy(1:s) = rdata(1:s)
  
    p => fft_new_plan(s,rdata(1:s),rdata(1:s),FFT_FORWARD)
    call fft_apply_plan(p,rdata(1:s),rdata(1:s))
    call fft_delete_plan(p)

    p => fft_new_plan(s,rdata(1:s),rdata(1:s),FFT_INVERSE,FFT_NORMALIZE)
    call fft_apply_plan(p,rdata(1:s),rdata(1:s))
    call fft_delete_plan(p)
 
    ierr = MAXVAL(ABS( rdata(1:s) - rdata_copy(1:s) ))
#ifdef FFTW_F2003
    if( ierr > err_max ) then
      print*,'Average error too big ',ierr
      stop ''
    endif
#endif
   enddo
  enddo
  print *, 'OK'


#if _DEFAULTFFTLIB!=FFTPACK_MOD
  print *,'-------------------------------------------------'
  print * ,'REAL TO COMPLEX and COMPLEX TO REAL '
  do i=imin,imax
   do h=1,hmax
    s = 2**i

    do j=1,s
      CALL RANDOM_NUMBER(rdata(j))
    enddo
    rdata_copy(1:s) = rdata(1:s)
  
    p => fft_new_plan(s,rdata(1:s),data_comp(1:s/2+1))
    call fft_apply_plan(p,rdata(1:s),data_comp(1:s/2+1))
    call fft_delete_plan(p)

    p => fft_new_plan(s,data_comp(1:s/2+1),rdata(1:s),FFT_NORMALIZE)
    call fft_apply_plan(p,data_comp(1:s/2+1),rdata(1:s))
    call fft_delete_plan(p)
    ierr = MAXVAL(ABS(rdata(1:s) - rdata_copy(1:s)))
    if( ierr > err_max ) then
      stop 'Average error too big'
    endif
   enddo
  enddo
  print *, 'OK', ierr
!----------------------------
  do i=imin,imax
   do h=1,hmax
    s = 2**i

    do j=1,s
      CALL RANDOM_NUMBER(rdata_comp(j))
    enddo
    rdata_copy(1:s) = rdata_comp(1:s)
  
    p => fft_new_plan(s,rdata_comp(1:s),data_comp(1:s/2+1),FFT_NORMALIZE)
    call fft_apply_plan(p,rdata_comp(1:s),data_comp(1:s/2+1))
    call fft_delete_plan(p)
    
    p => fft_new_plan(s,data_comp(1:s/2+1),rdata_comp(1:s))
    call fft_apply_plan(p,data_comp(1:s/2+1),rdata_comp(1:s))
    call fft_delete_plan(p)
    ierr = MAXVAL(ABS(rdata_comp(1:s) - rdata_copy(1:s)))
    if( ierr > err_max ) then
      stop 'Average error too big'
    endif
   enddo
  enddo
  print *, 'OK', ierr

  print *,'-------------------------------------------------'
  print * ,'COMPLEX TO COMPLEX 2D'
  do i=10,10
     do h=1,hmax
        s = m1!2**i
        t = m2!2**(i-1)
        
        do j=1,s
           do k=1,t
              CALL RANDOM_COMPLEX(data_comp2d(j,k))
           enddo
        enddo
        data_copy2d(1:s,1:t) = data_comp2d(1:s,1:t)
        
        p => fft_new_plan(s,t,data_comp2d(1:s,1:t),data_comp2d(1:s,1:t), &
             FFT_FORWARD)
        call fft_apply_plan(p,data_comp2d(1:s,1:t),data_comp2d(1:s,1:t))
        call fft_delete_plan(p)
        
        p => fft_new_plan(s,t,data_comp2d(1:s,1:t),data_comp2d(1:s,1:t), &
             FFT_INVERSE,FFT_NORMALIZE)
        call fft_apply_plan(p,data_comp2d(1:s,1:t),data_comp2d(1:s,1:t))
        call fft_delete_plan(p)
        ierr = 0_f64
        do j=1,t
           ierr = MAX(ERROR_MAX(data_comp2d(1:s,j) - data_copy2d(1:s,j)),ierr)
        enddo
        if( ierr > err_max ) then
           stop 'Average error too big'
        endif
     enddo
  enddo

  print * ,'COMPLEX TO COMPLEX 2D: WITH FLAGS'
  do i=10,10
     do h=1,hmax
        s = m1!2**i
        t = m2!2**(i-1)
        
        do j=1,s
           do k=1,t
              CALL RANDOM_COMPLEX(data_comp2d(j,k))
           enddo
        enddo
        data_copy2d(1:s,1:t) = data_comp2d(1:s,1:t)
        
        p => fft_new_plan(s,t,data_comp2d(1:s,1:t),data_comp2d(1:s,1:t), &
             FFT_FORWARD, FFT_ONLY_SECOND_DIRECTION)
        call fft_apply_plan(p,data_comp2d(1:s,1:t),data_comp2d(1:s,1:t))
        call fft_delete_plan(p)
        
        p => fft_new_plan(s,t,data_comp2d(1:s,1:t),data_comp2d(1:s,1:t), &
             FFT_INVERSE, FFT_NORMALIZE+FFT_ONLY_SECOND_DIRECTION)
        call fft_apply_plan(p,data_comp2d(1:s,1:t),data_comp2d(1:s,1:t))
        call fft_delete_plan(p)
        err_var = 0_f64
        do j=1,t 
           err_var = MAX(ERROR_MAX(data_comp2d(1:s,j)-data_copy2d(1:s,j)), &
                err_var)
        enddo
        print *, 'max_error = ', err_var
        if( ierr > err_max ) then
           stop 'Average error too big'
        endif
     enddo
  enddo



  print *, 'OK', ierr
#endif

#if _DEFAULTFFTLIB==SLLFFTMOD
  print *,'-------------------------------------------------'
  print * ,'COMPLEX TO COMPLEX 2D in one direction only'
  do i=10,10
     do h=1,hmax
        s = m1!2**i
        t = m2!2**(i-1)
        
        do j=1,s
           do k=1,t
              CALL RANDOM_COMPLEX(data_comp2d(j,k))
           enddo
        enddo
        data_copy2d(1:s,1:t) = data_comp2d(1:s,1:t)
        
        p => fft_new_plan(s,t,data_comp2d(1:s,1:t),data_comp2d(1:s,1:t), &
             FFT_FORWARD,FFT_ONLY_FIRST_DIRECTION)
        call fft_apply_plan(p,data_comp2d(1:s,1:t),data_comp2d(1:s,1:t))

        call fft_delete_plan(p)
        p => fft_new_plan(s,t,data_comp2d(1:s,1:t),data_comp2d(1:s,1:t), &
             FFT_FORWARD,FFT_ONLY_SECOND_DIRECTION)
        call fft_apply_plan(p,data_comp2d(1:s,1:t),data_comp2d(1:s,1:t))
        call fft_delete_plan(p)
        ! The following are 2 ways to do the same thing.
#if 1        
        p => fft_new_plan(s,t,data_comp2d(1:s,1:t),data_comp2d(1:s,1:t), &
             FFT_INVERSE,FFT_NORMALIZE)
        call fft_apply_plan(p,data_comp2d(1:s,1:t),data_comp2d(1:s,1:t))
        call fft_delete_plan(p)
#endif
#if 0
   p => fft_new_plan(s,t,data_comp2d(1:s,1:t),data_comp2d(1:s,1:t), &
             FFT_INVERSE,FFT_ONLY_SECOND_DIRECTION+FFT_NORMALIZE)
        call fft_apply_plan(p,data_comp2d(1:s,1:t),data_comp2d(1:s,1:t))
        call fft_delete_plan(p)

   p => fft_new_plan(s,t,data_comp2d(1:s,1:t),data_comp2d(1:s,1:t), &
             FFT_INVERSE,FFT_ONLY_FIRST_DIRECTION+FFT_NORMALIZE)
        call fft_apply_plan(p,data_comp2d(1:s,1:t),data_comp2d(1:s,1:t))
        call fft_delete_plan(p)
#endif
        ierr = 0_f64
        do j=1,t
           ierr = MAX(ERROR_MAX(data_comp2d(1:s,j) - data_copy2d(1:s,j)),ierr)
        enddo
        print *, 'error_max = ', ierr
        if( ierr > err_max ) then
           stop 'Average error too big'
        endif
     enddo
  enddo
  print *, 'OK', ierr
!-----------------------------
  do i=10,10
   do h=1,hmax
    s = m1!2**i
    t = m2!2**(i-1)

    do j=1,s
     do k=1,t
      CALL RANDOM_COMPLEX(data_comp2d(j,k))
     enddo
    enddo
    data_copy2d(1:s,1:t) = data_comp2d(1:s,1:t)
  
    p => fft_new_plan(s,t,data_comp2d(1:s,1:t),data_comp2d(1:s,1:t), &
                  FFT_FORWARD,FFT_ONLY_FIRST_DIRECTION + FFT_NORMALIZE)
    call fft_apply_plan(p,data_comp2d(1:s,1:t),data_comp2d(1:s,1:t))
    call fft_delete_plan(p)

    p => fft_new_plan(s,t,data_comp2d(1:s,1:t),data_comp2d(1:s,1:t), &
                  FFT_FORWARD, FFT_ONLY_SECOND_DIRECTION + FFT_NORMALIZE)
    call fft_apply_plan(p,data_comp2d(1:s,1:t),data_comp2d(1:s,1:t))
    call fft_delete_plan(p)

    p => fft_new_plan(s,t,data_comp2d(1:s,1:t),data_comp2d(1:s,1:t),FFT_INVERSE)
    call fft_apply_plan(p,data_comp2d(1:s,1:t),data_comp2d(1:s,1:t))
    call fft_delete_plan(p)
    ierr = 0_f64
    do j=1,t
      ierr = MAX(ERROR_MAX(data_comp2d(1:s,j) - data_copy2d(1:s,j)),ierr)
    enddo
    if( ierr > err_max ) then
      stop 'Average error too big'
    endif
   enddo
  enddo
  print *, 'OK', ierr
!-----------------------------
#endif

#if _DEFAULTFFTLIB!=FFTPACK_MOD
  print *,'-------------------------------------------------'
  print * ,'REAL TO COMPLEX 2D and COMPLEX TO REAL 2D'
  do i=10,10
   do h=1,hmax
    s = m1!2**i
    t = m2!2**(i-1)


    do j=1,s
     do k=1,t
      CALL RANDOM_NUMBER(data_real2d(j,k))
     enddo
    enddo
    rdata_copy2d(1:s,1:t) = data_real2d(1:s,1:t)
 
    p => fft_new_plan(s,t,data_real2d(1:s,1:t),data_comp2d(1:s/2+1,1:t))
    call fft_apply_plan(p,data_real2d(1:s,1:t),data_comp2d(1:s/2+1,1:t))
    call fft_delete_plan(p)

    p => fft_new_plan(s,t,data_comp2d(1:s/2+1,1:t),data_real2d(1:s,1:t),FFT_NORMALIZE)
    call fft_apply_plan(p,data_comp2d(1:s/2+1,1:t),data_real2d(1:s,1:t))
    call fft_delete_plan(p)
 
    ierr = 0_f64
    do j=1,t
      ierr = MAX(MAXVAL(ABS(data_real2d(1:s,j) - rdata_copy2d(1:s,j))),ierr)
    enddo
    if( ierr > err_max ) then
      print*, 'Average error too big', ierr
      stop ''
    endif
   enddo
  enddo
  print *, 'OK', ierr

!------------------------
  do i=10,10
   do h=1,hmax
    s = m1!2**i
    t = m2!2**(i-1)

    do j=1,s
     do k=1,t
      CALL RANDOM_NUMBER(data_real2d(j,k))
     enddo
    enddo
    rdata_copy2d(1:s,1:t) = data_real2d(1:s,1:t)
 
    p => fft_new_plan(s,t,data_real2d(1:s,1:t),data_comp2d(1:s/2+1,1:t),FFT_NORMALIZE)
    call fft_apply_plan(p,data_real2d(1:s,1:t),data_comp2d(1:s/2+1,1:t))
    call fft_delete_plan(p)

    p => fft_new_plan(s,t,data_comp2d(1:s/2+1,1:t),data_real2d(1:s,1:t))
    call fft_apply_plan(p,data_comp2d(1:s/2+1,1:t),data_real2d(1:s,1:t))
    call fft_delete_plan(p)

    ierr = 0_f64
    do j=1,t
      ierr = MAX(MAXVAL(ABS(data_real2d(1:s,j) - rdata_copy2d(1:s,j))),ierr)
    enddo
    if( ierr > err_max ) then
      print*, 'Average error too big', ierr
      stop ''
    endif
   enddo
  enddo
  print *, 'OK', ierr
#endif

contains

  FUNCTION ERROR_MAX(tab) RESULT(error)
    sll_comp64, DIMENSION(:) :: tab
    sll_real64 :: error

    error = MAX( MAXVAL(ABS(REAL(tab))) , MAXVAL(ABS(DIMAG(tab))) )
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
    c = CMPLX(realpart,imagpart,KIND=KIND(1.D0))
  END SUBROUTINE

end program unit_test
