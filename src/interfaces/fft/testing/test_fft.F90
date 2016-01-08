program unit_test
#include "sll_working_precision.h"
  use sll_m_fft
  implicit none

#define SLLFFT_MOD 0
#define FFTPACK_MOD 100
#define FFTW_MOD 1000000000

  type(sll_t_fft_plan), pointer :: p => null()
  type(sll_t_fft_plan), pointer :: pf => null()
  type(sll_t_fft_plan), pointer :: pb => null()

  sll_real64, parameter  :: err_max = 10E-14_f64
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
  sll_real64 :: ierr  ! this is not used as an integer below, bad name
  sll_int32 :: i,j,s,h,k,t
  sll_int32 :: rnd_seed_size
  sll_int32, allocatable :: rnd_seed(:) !< Random seed.

  ! Aligned data
  sll_comp64, dimension(:), pointer :: in
  sll_real64, dimension(:), pointer :: ar_data
 
  call sll_s_print_defaultfftlib()
  
  ! Set some random seed
  call random_seed (size=rnd_seed_size)
  allocate(rnd_seed(rnd_seed_size))
  do j=1, rnd_seed_size
     rnd_seed(j) = 100+15*j
  end do
  call random_seed (put=rnd_seed)
  

! test getter and setter functions in real case
  s = 2**imax
  do j=1,s
    CALL RANDOM_NUMBER(rdata(j))
  enddo
  rdata_copy(1:s) = rdata(1:s)
  ALLOCATE(ar_data(1:s))
  
  allocate(p)
  call sll_s_fft_init_plan_r2r_1d(p, s,rdata(1:s),rdata(1:s),sll_p_fft_forward)
  !call sll_s_fft_apply_plan_r2r_1d(p, rdata(1:s), rdata(1:s))
  do j=1,s/2+1
     data_comp(j) =  sll_f_fft_get_mode_r2c_1d(p,rdata(1:s), j-1)
     call sll_s_fft_set_mode_c2r_1d(p,ar_data, data_comp(j),j-1) 
  end do
  ierr = MAXVAL(ABS(ar_data - rdata_copy(1:s)))
  if( ierr .ne. 0_f64 ) then
    print *,'Average error too big',ierr
    stop
  else
    print *,'get and set mode real ok'
  endif
  call sll_s_fft_delete_plan(p)
  deallocate(p)
  DEALLOCATE(ar_data)


  print *,'-------------------------------------------------'
  print * ,'COMPLEX TO COMPLEX'
  do i=imin,imax
   do h=1,hmax
    s = 2**i

    do j=1,s
      CALL RANDOM_COMPLEX(data_comp(j))
    enddo
    data_copy(1:s) = data_comp(1:s)
  
    p => sll_f_fft_new_plan_c2c_1d(s,data_comp(1:s),data_comp(1:s),sll_p_fft_forward)
    call sll_s_fft_apply_plan_c2c_1d(p,data_comp(1:s),data_comp(1:s))
    call sll_s_fft_delete_plan(p)
    deallocate(p)

    p => sll_f_fft_new_plan_c2c_1d(s,data_comp(1:s),data_comp(1:s),sll_p_fft_backward,normalized=.TRUE.)
    call sll_s_fft_apply_plan_c2c_1d(p,data_comp(1:s),data_comp(1:s))
    call sll_s_fft_delete_plan(p)
    deallocate(p)
    ierr = ERROR_MAX(data_comp(1:s) - data_copy(1:s))
    if( ierr > err_max ) then
      stop 'Average error too big'
    endif
   enddo
  enddo
  print *, 'OK'


  print *,'-------------------------------------------------'
  print * ,'COMPLEX TO COMPLEX WITH ALIGNED MEMORY AND PLANNER OPTIMIZATION'

  ! Allocate aligned memory
  in => sll_f_fft_allocate_aligned_complex(m1)
  !call sll_f_fft_allocate_aligned_complex(m1, p_in, in)
 
  ! Initialize the plans 
  pf => sll_f_fft_new_plan_c2c_1d(m1,in,in,sll_p_fft_forward, aligned = .TRUE., optimization = sll_p_fft_measure)
  pb => sll_f_fft_new_plan_c2c_1d(m1,in,in,sll_p_fft_backward,normalized=.TRUE., aligned = .TRUE., optimization = sll_p_fft_patient)
 
  ! Initialize the data (note that this has to be done after initializing the plans due to the optimization level.
  do j=1,m1
     CALL RANDOM_COMPLEX(data_comp(j))
  enddo
  in = data_comp(1:m1)
    
  ! Forward transform
  call sll_s_fft_apply_plan_c2c_1d(pf,in,in)
  call sll_s_fft_delete_plan(pf)
  deallocate(pf)

  ! Backward transform
  call sll_s_fft_apply_plan_c2c_1d(pb,in,in)
  call sll_s_fft_delete_plan(pb)
  deallocate(pb)
  ierr = ERROR_MAX(data_comp(1:m1) - in)
  if( ierr > err_max ) then
     stop 'Average error too big'
  endif

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
  
    allocate(p)
    call sll_s_fft_init_plan_r2r_1d(p,s,rdata(1:s),rdata(1:s),sll_p_fft_forward)
    call sll_s_fft_apply_plan_r2r_1d(p,rdata(1:s),rdata(1:s))
    call sll_s_fft_delete_plan(p)
    deallocate(p)

    allocate(p)
    call sll_s_fft_init_plan_r2r_1d(p,s,rdata(1:s),rdata(1:s),sll_p_fft_backward,normalized = .TRUE.)
    call sll_s_fft_apply_plan_r2r_1d(p,rdata(1:s),rdata(1:s))
    call sll_s_fft_delete_plan(p)
    deallocate(p)

    ierr = MAXVAL(ABS( rdata(1:s) - rdata_copy(1:s) ))
    if( ierr > err_max ) then
      print*,'Average error too big ',ierr
      stop ''
    endif
   enddo
  enddo
  print *, 'OK'


  print *,'-------------------------------------------------'
  print * ,'REAL TO REAL WITH ALIGNED MEMORY AND PLANNER OPTIMIZATION'

  ! Allocate aligned memory
  ar_data => sll_f_fft_allocate_aligned_real(m1)
 
  ! Initialize the plans
  allocate(pf)
  call sll_s_fft_init_plan_r2r_1d(pf,m1,ar_data,ar_data,sll_p_fft_forward, aligned = .TRUE., optimization = sll_p_fft_measure)
  allocate(pb)
  call sll_s_fft_init_plan_r2r_1d(pb,m1,ar_data,ar_data,sll_p_fft_backward,normalized=.TRUE., aligned = .TRUE., optimization = sll_p_fft_patient)
 
  ! Initialize the data (note that this has to be done after initializing the plans due to the optimization level.
  do j=1,m1
     CALL RANDOM_NUMBER(rdata(j))
  enddo
  ar_data = rdata(1:m1)
    
  ! Forward transform
  call sll_s_fft_apply_plan_r2r_1d(pf,ar_data,ar_data)
  call sll_s_fft_delete_plan(pf)
  deallocate(pf)
  
  ! Backward transform
  call sll_s_fft_apply_plan_r2r_1d(pb,ar_data, ar_data)
  call sll_s_fft_delete_plan(pb)
  deallocate(pb)
  ierr = MAXVAL(ABS(rdata(1:m1) - ar_data))
  if( ierr > err_max ) then
     stop 'Average error too big'
  endif

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
  
    allocate(p)
    call sll_s_fft_init_plan_r2c_1d(p,s,rdata(1:s),data_comp(1:s/2+1))
    call sll_s_fft_apply_plan_r2c_1d(p,rdata(1:s),data_comp(1:s/2+1))
    call sll_s_fft_delete_plan(p)
    deallocate(p)

    allocate(p)
    call sll_s_fft_init_plan_c2r_1d(p,s,data_comp(1:s/2+1),rdata(1:s), normalized = .TRUE.)
    call sll_s_fft_apply_plan_c2r_1d(p,data_comp(1:s/2+1),rdata(1:s))
    call sll_s_fft_delete_plan(p)
    deallocate(p)
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
  
    allocate(p)
    call sll_s_fft_init_plan_r2c_1d(p,s,rdata_comp(1:s),data_comp(1:s/2+1), normalized = .TRUE.)
    call sll_s_fft_apply_plan_r2c_1d(p,rdata_comp(1:s),data_comp(1:s/2+1))
    call sll_s_fft_delete_plan(p)
    deallocate(p)

    allocate(p)
    call sll_s_fft_init_plan_c2r_1d(p,s,data_comp(1:s/2+1),rdata_comp(1:s))
    call sll_s_fft_apply_plan_c2r_1d(p,data_comp(1:s/2+1),rdata_comp(1:s))
    call sll_s_fft_delete_plan(p)
    deallocate(p)
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
        
        allocate(p)
        call sll_s_fft_init_plan_c2c_2d(p,s,t,data_comp2d(1:s,1:t),data_comp2d(1:s,1:t), &
             sll_p_fft_forward)
        call sll_s_fft_apply_plan_c2c_2d(p,data_comp2d(1:s,1:t),data_comp2d(1:s,1:t))
        call sll_s_fft_delete_plan(p)
        deallocate(p)

        allocate(p)
        call sll_s_fft_init_plan_c2c_2d(p,s,t,data_comp2d(1:s,1:t),data_comp2d(1:s,1:t), &
             sll_p_fft_backward,normalized = .TRUE.)
        call sll_s_fft_apply_plan_c2c_2d(p,data_comp2d(1:s,1:t),data_comp2d(1:s,1:t))
        call sll_s_fft_delete_plan(p)
        deallocate(p)
        ierr = 0._f64
        do j=1,t
           ierr = MAX(ERROR_MAX(data_comp2d(1:s,j) - data_copy2d(1:s,j)),ierr)
        enddo
        if( ierr > err_max ) then
           stop 'Average error too big'
        endif
     enddo
  enddo

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
 
    allocate(p)
    call sll_s_fft_init_plan_r2c_2d(p,s,t,data_real2d(1:s,1:t),data_comp2d(1:s/2+1,1:t))
    call sll_s_fft_apply_plan_r2c_2d(p,data_real2d(1:s,1:t),data_comp2d(1:s/2+1,1:t))
    call sll_s_fft_delete_plan(p)
    deallocate(p)

    allocate(p)
    call sll_s_fft_init_plan_c2r_2d(p,s,t,data_comp2d(1:s/2+1,1:t),data_real2d(1:s,1:t),normalized = .TRUE.)
    call sll_s_fft_apply_plan_c2r_2d(p,data_comp2d(1:s/2+1,1:t),data_real2d(1:s,1:t))
    call sll_s_fft_delete_plan(p)
    deallocate(p)
 
    ierr = 0._f64
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
 
    allocate(p)
    call sll_s_fft_init_plan_r2c_2d(p,s,t,data_real2d(1:s,1:t),data_comp2d(1:s/2+1,1:t),&
         normalized = .TRUE.)
    call sll_s_fft_apply_plan_r2c_2d(p,data_real2d(1:s,1:t),data_comp2d(1:s/2+1,1:t))
    call sll_s_fft_delete_plan(p)
    deallocate(p)

    allocate(p)
    call sll_s_fft_init_plan_c2r_2d(p,s,t,data_comp2d(1:s/2+1,1:t),data_real2d(1:s,1:t))
    call sll_s_fft_apply_plan_c2r_2d(p,data_comp2d(1:s/2+1,1:t),data_real2d(1:s,1:t))
    call sll_s_fft_delete_plan(p)
    deallocate(p)

    ierr = 0._f64
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

    error = MAX( MAXVAL(ABS(REAL(tab))) , MAXVAL(ABS(aimag(tab))) )
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
    sll_comp64, intent(out) :: c
    sll_real64 :: realpart, imagpart
   
    CALL init_random_seed() 
    CALL RANDOM_NUMBER(realpart)
    CALL RANDOM_NUMBER(imagpart)
    c = CMPLX(realpart,imagpart,kind=f64)
  END SUBROUTINE

end program unit_test
