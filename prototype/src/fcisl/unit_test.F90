program fcisl_test
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
  use sll_fcisl_module
  use sll_constants
  use sll_module_advection_1d_periodic
  
  implicit none
  
  sll_int32 :: Nc_x1
  sll_int32 :: Nc_x2
  sll_real64 :: A1
  sll_real64 :: A2
  sll_real64 :: A1_0
  sll_real64 :: A2_0
  sll_int32 :: k_mode
  sll_real64 :: x1
  sll_real64 :: x2
  sll_real64 :: dt
  sll_int32 :: nb_step
  sll_real64, dimension(:,:), allocatable :: f
  sll_real64, dimension(:,:), allocatable :: f_init
  sll_real64, dimension(:,:), allocatable :: f_exact
  sll_real64, dimension(:,:), allocatable :: f_new
  sll_real64 :: x1_min
  sll_real64 :: x1_max
  sll_real64 :: x2_min
  sll_real64 :: x2_max
  sll_real64 :: delta_x1
  sll_real64 :: delta_x2
  type(sll_oblic_derivative), pointer :: deriv
  class(sll_advection_1d_base), pointer :: adv
  sll_int32 :: degree
  sll_real64, dimension(:,:), allocatable :: phi
  sll_real64, dimension(:,:), allocatable :: phi_at_aligned
  sll_int32, dimension(:), allocatable :: spaghetti_index
  sll_real64, dimension(:,:), allocatable :: D_phi
  sll_real64, dimension(:), allocatable :: buf1d_spaghetti
  sll_real64, dimension(:), allocatable :: buf1d_spaghetto
  sll_real64 :: iota
  sll_int32 :: ierr
  sll_int32 :: shift
  sll_real64 :: iota_modif
  sll_int32 :: spaghetti_size
  sll_int32 :: num_spaghetti
  sll_int32 :: shift_guess
  logical :: use_shift_guess
  sll_int32 :: spaghetti_size_guess
  sll_int32 :: i1
  sll_int32 :: i2
  sll_real64 :: tau
  sll_int32 :: i
  sll_int32 :: s
  character(len=256) :: filename
  sll_int32 :: IO_stat
  sll_int32, parameter  :: input_file = 99


  ! namelists for data input
  namelist /params/ &
    x1_min, &
    x1_max, &
    x2_min, &
    x2_max, &
    Nc_x1, &
    Nc_x2, &
    degree, &
    iota, &
    tau, &
    dt, &
    nb_step, &
    spaghetti_size_guess, &
    shift_guess, &
    use_shift_guess

  !set default parameters
  x1_min = 0._f64
  x1_max = 1._f64
  x2_min = 0._f64
  x2_max = 1._f64
  Nc_x1 = 512
  Nc_x2 = 32
  degree = 4
  iota = 4.8_f64
  tau = 0.14_f64
  dt = 0.1_f64
  nb_step = 10  
  spaghetti_size_guess = Nc_x1
  shift_guess = 0
  use_shift_guess = .false.

  call get_command_argument(1, filename)
  if (len_trim(filename) .ne. 0)then
    print*,'#read namelist'
    open(unit = input_file, file=trim(filename)//'.nml',IOStat=IO_stat)
      if( IO_stat /= 0 ) then
        print *, '#fcisl_test failed to open file ', trim(filename)//'.nml'
        stop
      end if
      read(input_file, params) 
    close(input_file)    
  else
    print *,'#use default parameters'  
  endif
  
  print *,'#params are'
  print *,'#x1_min=',x1_min
  print *,'#x1_max=',x1_max
  print *,'#x2_min=',x2_min
  print *,'#x2_max=',x2_max
  print *,'#Nc_x1=',Nc_x1
  print *,'#Nc_x2=',Nc_x2
  print *,'#degree=',degree
  print *,'#iota=',iota
  print *,'#tau=',tau
  print *,'#dt=',dt
  print *,'#nb_step=',nb_step
  print *,'#spaghetti_size_guess=',spaghetti_size_guess
  print *,'#shift_guess=',shift_guess
  print *,'#use_shift_guess=',use_shift_guess



  
  SLL_ALLOCATE(phi(Nc_x1+1,Nc_x2+1),ierr)
  SLL_ALLOCATE(phi_at_aligned(Nc_x1+1,Nc_x2+1),ierr)
  SLL_ALLOCATE(spaghetti_index(Nc_x1+1),ierr)
  SLL_ALLOCATE(D_phi(Nc_x1+1,Nc_x2+1),ierr)
  SLL_ALLOCATE(f(Nc_x1+1,Nc_x2+1),ierr)
  SLL_ALLOCATE(f_init(Nc_x1+1,Nc_x2+1),ierr)
  SLL_ALLOCATE(f_exact(Nc_x1+1,Nc_x2+1),ierr)
  SLL_ALLOCATE(f_new(Nc_x1+1,Nc_x2+1),ierr)
  SLL_ALLOCATE(buf1d_spaghetti((Nc_x1+1)*(Nc_x2+1)),ierr)
  SLL_ALLOCATE(buf1d_spaghetto((Nc_x1+1)*(Nc_x2+1)),ierr)
  delta_x1 = (x1_max-x1_min)/real(Nc_x1,f64)
  delta_x2 = (x2_max-x2_min)/real(Nc_x2,f64)  


  do i2=1,Nc_x2+1
    x2 = x2_min+real(i2-1,f64)*delta_x2
    do i1=1,Nc_x1+1
      x1 = x1_min+real(i1-1,f64)*delta_x1
      x2 = x2_min+real(i2-1,f64)*delta_x2
      f_init(i1,i2) = sin(2._f64*sll_pi*real(k_mode,f64)*(-A2_0*x1+A1_0*x2))
      x1 = x1 - A1*real(nb_step,f64)*dt
      x2 = x2 - A2*real(nb_step,f64)*dt
      f_exact(i1,i2) = sin(2._f64*sll_pi*real(k_mode,f64)*(-A2_0*x1+A1_0*x2))
    enddo
  enddo

  
  phi = 1._f64
  call compute_oblic_shift( &
    iota, &
    Nc_x1, &
    shift, &
    iota_modif)
  !spaghetti_size_guess = Nc_x1
  print *,'#shift from iota=',shift
  print *,'#iota_modif=',iota_modif
  call compute_spaghetti_size_from_shift( &
    Nc_x1, &
    shift, &
    spaghetti_size)
  print *,'#spaghetti_size from iota modif=',spaghetti_size  

  if(use_shift_guess .eqv. .false.)then
    shift_guess = shift
  endif
  print *,'#iota_modif=',iota_modif

  print *,'#iota_guess=',iota
  print *,'#spaghetti_size_guess=',spaghetti_size_guess
  print *,'#spaghetti_size between 1 and ',Nc_x1
  call compute_spaghetti_and_shift_from_guess( &
    Nc_x1, &
    Nc_x2, &
    iota, &
    spaghetti_size_guess, &
    shift, &
    spaghetti_size)

  print *,'#shift=',shift
  print *,'#spaghetti_size=',spaghetti_size
  call compute_iota_from_shift(Nc_x1,shift,iota_modif)
  print *,'#iota_modif=',iota_modif
  !stop



!we should choose spaghetti_size_guess
!and iota_guess

!to get then spaghetti_size and shift, with iota = shift/Nc_x1
!spaghetti_size near to spaghetti_size_guess 
! and iota near to iota_guess
! this permits to avoid to declare oblic 1d advectors of several sizes
  
  adv => new_periodic_1d_advector( &
    Nc_x1, &
    x1_min, &
    x1_max, &
    SPLINE, & 
    4) 

  call compute_at_aligned( &
    phi, &
    phi_at_aligned, &
    Nc_x1, &
    Nc_x2, &
    adv, &
    x1_min, &
    x1_max, &
    iota_modif )

  print *,'#phi_bounds',minval(phi),maxval(phi)
  print *,'#phi_at_aligned_bounds',minval(phi_at_aligned),maxval(phi_at_aligned)
  
  
  
  call compute_spaghetti_size_from_shift( &
    Nc_x1, &
    shift, &
    spaghetti_size)
  print *,'#shift=',shift,spaghetti_size
  
  !do shift=-2*Nc_x1,2*Nc_x1
  call compute_spaghetti( &
    Nc_x1, &
    Nc_x2, &
    shift, &
    spaghetti_index, &
    spaghetti_size)
  !enddo
  num_spaghetti = Nc_x1/spaghetti_size
  
  print *,'#num_spaghetti=',num_spaghetti
  
  !stop  
  call load_spaghetti( &
    phi_at_aligned, &
    buf1d_spaghetti, &
    spaghetti_index, &
    Nc_x1+1, &
    Nc_x2+1) 
  s = 1  
  do i=1,num_spaghetti
    call adv%advect_1d_constant( &
        tau, &
        dt, &
        buf1d_spaghetti(s:s+spaghetti_size*Nc_x2), &
        buf1d_spaghetto(1:spaghetti_size*Nc_x2))      

    buf1d_spaghetti(s:s+spaghetti_size*Nc_x2)=buf1d_spaghetto(1:spaghetti_size*Nc_x2) 
    
    s = s+spaghetti_size*Nc_x2
  enddo   

  call unload_spaghetti( &
    buf1d_spaghetti, &
    phi_at_aligned, &
    spaghetti_index, &
    Nc_x1+1, &
    Nc_x2+1) 
  
  print *,'#maxval phi_at_aligned_shifted=',maxval(phi_at_aligned),minval(phi_at_aligned)
  
  
  stop    

  
  deriv => new_oblic_derivative( &
    degree, &
    x1_min, &
    x1_max, &
    x2_min, &
    x2_max, &
    Nc_x1, &
    Nc_x2, &
    adv)

!      call adv%advect_1d_constant( &
!        tau, &
!        -0.125_f64, &
!        phi(1:Nc_x1+1,29), &
!        deriv%buf(1,1:Nc_x1+1))      


!  call compute_oblic_derivative( &
!    deriv, &
!    tau, &
!    phi, &
!    D_phi)
  
  

  print *,'#hello'
  
  print *,'#PASSED'



end program