program aligned_derivative_2d
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
  use sll_fcisl_module
  use sll_constants
  use sll_module_advection_1d_periodic
  use sll_module_derivative_2d_oblic
  
  implicit none
  
  type(oblic_2d_derivative), pointer :: deriv
  class(sll_advection_1d_base), pointer :: adv_x1
  class(sll_advection_1d_base), pointer :: adv_x2
  sll_int32 :: i1
  sll_int32 :: i2
  sll_int32 :: Nc_x1
  sll_int32 :: Nc_x2
  sll_real64 :: A1
  sll_real64 :: A2
  sll_real64 :: A1_0
  sll_real64 :: A2_0
  sll_int32 :: k_mode
  sll_real64, dimension(:,:), allocatable :: phi
  sll_real64, dimension(:,:), allocatable :: phi_at_aligned
  sll_real64, dimension(:,:), allocatable :: D_phi_at_aligned
  sll_real64, dimension(:,:), allocatable :: Da_phi
  sll_real64, dimension(:,:), allocatable :: Dx1_phi
  sll_real64, dimension(:,:), allocatable :: phi_exact
  sll_real64, dimension(:,:), allocatable :: Dx2_phi
  sll_real64, dimension(:,:), allocatable :: Dx2_phi_exact
  !sll_real64, dimension(:) :: buf_x1
  !sll_real64, dimension(:) :: buf_x2
  sll_real64, dimension(:), allocatable :: buf_1d
  sll_real64, dimension(:), allocatable :: buf_1d_out
  sll_int32, dimension(:), allocatable :: spaghetti_index
  sll_int32 :: ierr
  sll_real64 :: x1
  sll_real64 :: x2
  sll_real64 :: x1_min
  sll_real64 :: x1_max
  sll_real64 :: x2_min
  sll_real64 :: x2_max
  sll_real64 :: delta_x1
  sll_real64 :: delta_x2
  sll_int32 :: step
  sll_real64 :: err
  sll_real64 :: alpha
  sll_int32 :: i0
  sll_real64 :: dt_loc
  sll_real64, dimension(:,:), allocatable :: buf
  sll_int32 :: d
  sll_int32 :: r
  sll_int32 :: s
  sll_real64, dimension(:), allocatable :: xx
  sll_real64, dimension(:), allocatable :: w
  sll_real64, dimension(:), allocatable :: x1_array
  sll_real64, dimension(:), allocatable :: x2_array
  sll_int32 :: ell
  sll_int32 :: i2_loc
  character(len=256) :: filename
  sll_int32 :: IO_stat
  sll_int32, parameter  :: input_file = 99
  sll_int32 :: i
  sll_int32 :: ii
  sll_int32 :: j
  sll_real64 :: iota
  sll_real64 :: iota_tau
  sll_int32 :: spaghetti_size_guess
  sll_int32 :: spaghetti_size
  sll_int32 :: shift
  sll_int32 :: num_spaghetti
  sll_real64 :: A
  character(len=256) :: advector_x1
  sll_int32 :: order_x1
  character(len=256) :: advector_x2
  sll_int32 :: order_x2
  
  
  ! namelists for data input
  namelist /params/ &
    k_mode, &
    Nc_x1, &
    Nc_x2, &
    d, &
    A1_0, &
    A2_0, &
    A1, &
    A2, &
    spaghetti_size_guess, &
    advector_x1, &
    order_x1, &
    x1_min, &
    x1_max, &
    x2_min, &
    x2_max
!    advector_x2, &
!    order_x2
  
  !initialization
  k_mode = 3
  Nc_x1 = 512
  Nc_x2 = 16
  d = 5
  spaghetti_size_guess = 512
  
  A1_0 = 3._f64
  A2_0 = 7._f64  ! we should assume A2>A1>=0
  
  A1 = 2.8357_f64
  A2 = 7.18459_f64

  !advector
  advector_x1 = "SLL_SPLINES"
  order_x1 = 4
!  advector_x2 = "SLL_SPLINES"
!  order_x2 = 4

  x1_min = 0._f64
  x1_max = 1._f64
  x2_min = 0._f64
  x2_max = 1._f64



  call get_command_argument(1, filename)
  if (len_trim(filename) .ne. 0)then
    print*,'#read namelist'
    open(unit = input_file, file=trim(filename)//'.nml',IOStat=IO_stat)
      if( IO_stat /= 0 ) then
        print *, '#aligned_translation_2d failed to open file ', trim(filename)//'.nml'
        stop
      end if
      read(input_file, params) 
    close(input_file)    
  else
    print *,'#use default parameters'  
  endif
  print *,'#k_mode=',k_mode
  print *,'#Nc_x1=',Nc_x1
  print *,'#Nc_x2=',Nc_x2
  print *,'#d=',d
  print *,'#A1_0',A1_0
  print *,'#A2_0',A2_0
  print *,'#A1=',A1
  print *,'#A2=',A2
  print *,'#advector_x1=', trim(advector_x1)
  print *,'#order_x1',order_x1
  
  delta_x1 = (x1_max-x1_min)/real(Nc_x1,f64)
  delta_x2 = (x2_max-x2_min)/real(Nc_x2,f64)  
  !we force to use even numer for having unique derivative at point
  r = -(d+1)/2
  s = (d+1)/2
  !r = -(d-1)/2
  !s = (d+1)/2
  SLL_ALLOCATE(xx(r:s),ierr)
  SLL_ALLOCATE(w(r:s),ierr)
  SLL_ALLOCATE(buf(r:s,Nc_x1+1),ierr)  
  SLL_ALLOCATE(phi(Nc_x1+1,Nc_x2+1),ierr)
  SLL_ALLOCATE(phi_at_aligned(Nc_x1+1,Nc_x2+1),ierr)
  SLL_ALLOCATE(D_phi_at_aligned(Nc_x1+1,Nc_x2+1),ierr)
  SLL_ALLOCATE(Da_phi(Nc_x1+1,Nc_x2+1),ierr)
  SLL_ALLOCATE(Dx1_phi(Nc_x1+1,Nc_x2+1),ierr)
  SLL_ALLOCATE(spaghetti_index(Nc_x1+1),ierr)
  SLL_ALLOCATE(phi_exact(Nc_x1+1,Nc_x2+1),ierr)
  SLL_ALLOCATE(Dx2_phi(Nc_x1+1,Nc_x2+1),ierr)
  SLL_ALLOCATE(Dx2_phi_exact(Nc_x1+1,Nc_x2+1),ierr)
  SLL_ALLOCATE(x1_array(Nc_x1+1),ierr)
  SLL_ALLOCATE(x2_array(Nc_x2+1),ierr)
  
  
  !print *,'#(r,s)=',r,s
  call compute_w_hermite(w,r,s)
  !print *,'#w=',w(r:s)
  
  do i=1,Nc_x1+1
    x1_array(i) = x1_min+real(i-1,f64)*delta_x1
  enddo
  do i=1,Nc_x2+1
    x2_array(i) = x2_min+real(i-1,f64)*delta_x2
  enddo

  do i1=r,s
    xx(i1) = real(i1,f64)
  enddo
  

  select case (advector_x1)
    case ("SLL_SPLINES") ! arbitrary order periodic splines
      adv_x1 => new_periodic_1d_advector( &
        Nc_x1, &
        x1_min, &
        x1_max, &
        SPLINE, & 
        order_x1) 
    case("SLL_LAGRANGE") ! arbitrary order Lagrange periodic interpolation
      adv_x1 => new_periodic_1d_advector( &
        Nc_x1, &
        x1_min, &
        x1_max, &
        LAGRANGE, & 
        order_x1)
    case default
      print*,'#advector in x1', advector_x1, ' not implemented'
      stop 
  end select    
!  select case (advector_x2)
!    case ("SLL_SPLINES") ! arbitrary order periodic splines
!      adv_x2 => new_periodic_1d_advector( &
!        Nc_x2, &
!        x2_min, &
!        x2_max, &
!        SPLINE, & 
!        order_x2) 
!    case("SLL_LAGRANGE") ! arbitrary order Lagrange periodic interpolation
!      adv_x2 => new_periodic_1d_advector( &
!        Nc_x2, &
!        x2_min, &
!        x2_max, &
!        LAGRANGE, & 
!        order_x2)
!    case default
!      print*,'#advector in x2', advector_x2, ' not implemented'
!      stop 
!  end select    



 
  do i2=1,Nc_x2+1
    x2 = x2_min+real(i2-1,f64)*delta_x2
    do i1=1,Nc_x1+1
      x1 = x1_min+real(i1-1,f64)*delta_x1
      x2 = x2_min+real(i2-1,f64)*delta_x2
      phi_exact(i1,i2) = sin(2._f64*sll_pi*real(k_mode,f64)*(-A2_0*x1/(x1_max-x1_min)+A1_0*x2/(x2_max-x2_min)))
      Dx2_phi_exact(i1,i2) = &
        A1_0*2._f64*sll_pi*real(k_mode,f64)/(x2_max-x2_min) &
        *cos(2._f64*sll_pi*real(k_mode,f64)*(-A2_0*x1/(x1_max-x1_min)+A1_0*x2/(x2_max-x2_min)))
    enddo
  enddo
  
  !classical method: derivative in z with finite difference  
  phi = phi_exact
  Dx2_phi = 0._f64
  err = 0._f64
  do i=1,Nc_x1+1
    call compute_derivative_periodic( &
    phi(i,1:Nc_x2+1), &
    Dx2_phi(i,1:Nc_x2+1), &
    Nc_x2, &
    w, &
    r, &
    s, &
    x2_max-x2_min)!call compute_derivative(phi(),Dx2_phi,Nc_x2)
  enddo   
  err = maxval(abs(Dx2_phi-Dx2_phi_exact))  
  print *,'#err for classical method=', &
    err,err/maxval(abs(Dx2_phi_exact))
  
  !new method: derivative along iota_modif
  iota = A1/A2
  call compute_spaghetti_and_shift_from_guess( &
    Nc_x1, &
    Nc_x2, &
    iota, &
    spaghetti_size_guess, &
    shift, &
    spaghetti_size)
  iota_tau = real(shift,f64)/real(Nc_x1,f64)
  
  print *,'#shift=',shift
  print *,'#spaghetti_size=',spaghetti_size
  
  
  num_spaghetti = Nc_x1/spaghetti_size
  
  print *,'#num_spaghetti=',num_spaghetti
  if(num_spaghetti*spaghetti_size .ne. Nc_x1)then
    print *,'#Problem of spaghetti size'
    stop
  endif
  
  SLL_ALLOCATE(buf_1d(spaghetti_size*Nc_x2+1),ierr)
  SLL_ALLOCATE(buf_1d_out(spaghetti_size*Nc_x2+1),ierr)
  
  ! compute phi on aligned mesh
  A = -real(shift,f64)/real(Nc_x2,f64)
  do i=1,Nc_x2+1
    call adv_x1%advect_1d_constant( &
      A, &
      real(i-1,f64)*delta_x1, &
      phi(1:Nc_x1+1,i), &
      phi_at_aligned(1:Nc_x1+1,i))      
  enddo
  
  !print *,'#diff phi-phi_at_aligned',maxval(abs(phi-phi_at_aligned))
  !print *,'#diff phi-phi_at_aligned at x2_min', &
  !  maxval(abs(phi(1:Nc_x1,1)-phi_at_aligned(1:Nc_x1,1)))
  !err = 0._f64
  !do i=1,Nc_x1+1
  !  err=max(err, &
  !    abs(phi(modulo(i+shift-1,Nc_x1)+1,Nc_x2+1)-phi_at_aligned(i,Nc_x2+1)))
  !enddo
  !print *,'#diff phi-phi_at_aligned at x2_max',err

  !  maxval(abs(phi(1:Nc_x1,Nc_x2+1)-phi_at_aligned(1:Nc_x1,Nc_x2+1)))

  !do i=1,Nc_x2+1
  !  print *,x2_min+real(i-1,f64)*delta_x2,phi(1,i),phi_at_aligned(1,i),phi_exact(1,i)
  !enddo
  !stop
  
  call compute_spaghetti( &
    Nc_x1, &
    Nc_x2, &
    shift, &
    spaghetti_index, &
    spaghetti_size)
  !do i=1,spaghetti_size
  !  print *,i,spaghetti_index(i)
  !enddo
  !stop  
  do i=1,num_spaghetti
    !load spaghetti_1d
    ell = 0
    do ii=1,spaghetti_size
      do j=1,Nc_x2
        ell = ell+1
        buf_1d(ell) = phi_at_aligned(i-1+spaghetti_index(ii),j)
      enddo
      err=0._f64
      if(ii<spaghetti_size)then
        err=max(err,abs(phi_at_aligned(i-1+spaghetti_index(ii),Nc_x2+1)- &
          phi_at_aligned(i-1+spaghetti_index(ii+1),1)))
      endif
      if(ii==spaghetti_size)then
        err=max(err,abs(phi_at_aligned(i-1+spaghetti_index(ii),Nc_x2+1)- &
          phi_at_aligned(i-1+spaghetti_index(1),1)))        
      endif      
    enddo
    !print *,'#err on spaghetti',i,err
    ell = ell+1
    buf_1d(ell) = buf_1d(1)!phi_at_aligned(i-1+spaghetti_index(spaghetti_size),Nc_x1+1)
!    do ell=1,spaghetti_size*Nc_x2+1
!      print *,x2_min+real(ell-1,f64)*delta_x2,buf_1d(ell)
!    enddo
!    stop
    !compute derivative
    call compute_derivative_periodic( &
      buf_1d(1:spaghetti_size*Nc_x2+1), &
      buf_1d_out(1:spaghetti_size*Nc_x2+1), &
      spaghetti_size*Nc_x2, &
      w, &
      r, &
      s, &
      real(spaghetti_size,f64)*(x2_max-x2_min))!call compute_derivative(phi(),Dx2_phi,Nc_x2)
    !do ii=1,spaghetti_size*Nc_x2+1
    !  print *,ii,buf_1d(ii),buf_1d_out(ii)
    !enddo
    !stop
    !unload spaghetti_1d
    ell = 0
    do ii=1,spaghetti_size
      do j=1,Nc_x2
        ell = ell+1
        D_phi_at_aligned(i-1+spaghetti_index(ii),j) = buf_1d_out(ell)
      enddo
    enddo
    do ii=1,spaghetti_size
      if(ii<spaghetti_size)then  
        D_phi_at_aligned(i-1+spaghetti_index(ii),Nc_x2+1) &
          = D_phi_at_aligned(i-1+spaghetti_index(ii+1),1)
      endif
      if(ii==spaghetti_size)then
        D_phi_at_aligned(i-1+spaghetti_index(ii),Nc_x2+1) &
          = D_phi_at_aligned(i-1+spaghetti_index(1),1)        
      endif   
    enddo
  enddo

  ! compute phi on aligned mesh
  A = real(shift,f64)/real(Nc_x2,f64)
  do i=1,Nc_x2+1
    call adv_x1%advect_1d_constant( &
      A, &
      real(i-1,f64)*delta_x1, &
      D_phi_at_aligned(1:Nc_x1+1,i), &
      Da_phi(1:Nc_x1+1,i))      
  enddo
  !Da_phi(1:Nc_x1+1,Nc_x2+1) = Da_phi(1:Nc_x1+1,1)
  
  
  do i=1,Nc_x2+1
    call compute_derivative_periodic( &
    phi(1:Nc_x1+1,i), &
    Dx1_phi(1:Nc_x1+1,i), &
    Nc_x1, &
    w, &
    r, &
    s, &
    x1_max-x1_min)!call compute_derivative(phi(),Dx2_phi,Nc_x2)
  enddo   
  
  !print *,'#x1_min,x1_max=',x1_min,x1_max
  !print *,'#x2_min,x2_max=',x2_min,x2_max
  
  !Dx2_phi = Da_phi - iota_tau*(x1_max-x1_min)/(x2_max-x2_min)*Dx1_phi
  Dx2_phi = Da_phi - iota_tau*Dx1_phi
  err = maxval(abs(Dx2_phi-Dx2_phi_exact))  
  print *,'#err for new method=', &
    err,err/maxval(abs(Dx2_phi_exact))
  !do i=1,Nc_x1+1
  !  print *,'#err',i,maxval(abs(Dx2_phi(i,:)-Dx2_phi_exact(i,:))) 
  !enddo
  !do i=1,Nc_x2+1
  !  print *,'#err2',i,maxval(abs(Dx2_phi(:,i)-Dx2_phi_exact(:,i))) 
  !enddo
  
  !do i=1,Nc_x2+1
  !  print *,x2_min+real(i-1,f64)*delta_x2,Da_phi(Nc_x1/2,i),Dx2_phi(Nc_x1/2,i), &
  !    Dx2_phi_exact(Nc_x1/2,i)
  !enddo
  
  
  !new method without spaghetti and with abstract interface 
  ! more in the spirit of Ottaviani 
  
  deriv => new_oblic_2d_derivative( &
    Nc_x1, &
    adv_x1, &
    Nc_x2, &
    x2_min, &
    x2_max, &
    r, &
    s )

  call compute_oblic_derivative_2d( &
    deriv, &
    A1, &
    A2, &
    phi, &
    Da_phi)

  Dx2_phi = (1._f64/A2)*Da_phi - (A1/A2)*Dx1_phi
  err = maxval(abs(Dx2_phi-Dx2_phi_exact))  
  print *,'#err for new method without spaghetti=', &
    err,err/maxval(abs(Dx2_phi_exact))

!do i=1,Nc_x1+1
!  print *,'#err',i,maxval(abs(Dx2_phi(i,:)-Dx2_phi_exact(i,:))) 
!enddo
!do i=1,Nc_x2+1
!  print *,'#err2',i,maxval(abs(Dx2_phi(:,i)-Dx2_phi_exact(:,i))) 
!enddo
!  
!do i=1,Nc_x2+1
!  print *,x2_min+real(i-1,f64)*delta_x2,Da_phi(Nc_x1/2,i),Dx2_phi(Nc_x1/2,i), &
!    Dx2_phi_exact(Nc_x1/2,i)
!enddo
  
  


!#ifndef NOHDF5
!      call plot_f_cartesian( &
!        0, &
!        Dx2_phi, &
!        x1_array, &
!        Nc_x1+1, &
!        x2_array, &
!        Nc_x2+1, &
!        'fold', 0._f64 )        
!!      call plot_f_cartesian( &
!!        iplot, &
!!        f_visu_light, &
!!        sim%x1_array_light, &
!!        np_x1_light, &
!!        node_positions_x2_light, &
!!        sim%num_dof_x2_light, &
!!        'light_f', time_init )        
!#endif
!
  
  
  print *,'#PASSED'
  
  


!contains
!
!#ifndef NOHDF5
!!*********************
!!*********************
!
!  !---------------------------------------------------
!  ! Save the mesh structure
!  !---------------------------------------------------
!  subroutine plot_f_cartesian( &
!    iplot, &
!    f, &
!    node_positions_x1, &
!    nnodes_x1, &
!    node_positions_x2, &
!    nnodes_x2, &
!    array_name, time)    
!    !mesh_2d)
!    use sll_xdmf
!    use sll_hdf5_io
!    sll_int32 :: file_id
!    sll_int32 :: error
!    sll_real64, dimension(:), intent(in) :: node_positions_x1
!    sll_real64, dimension(:), intent(in) :: node_positions_x2    
!     character(len=*), intent(in) :: array_name !< field name
!    sll_real64, dimension(:,:), allocatable :: x1
!    sll_real64, dimension(:,:), allocatable :: x2
!    sll_int32, intent(in) :: nnodes_x1
!    sll_int32, intent(in) :: nnodes_x2
!    sll_int32 :: i, j
!    sll_int32, intent(in) :: iplot
!    character(len=4)      :: cplot
!    sll_real64, dimension(:,:), intent(in) :: f
!    sll_real64 :: time
!    
!    if (iplot == 1) then
!
!      SLL_ALLOCATE(x1(nnodes_x1,nnodes_x2), error)
!      SLL_ALLOCATE(x2(nnodes_x1,nnodes_x2), error)
!      do j = 1,nnodes_x2
!        do i = 1,nnodes_x1
!          x1(i,j) = node_positions_x1(i) !x1_min+real(i-1,f32)*dx1
!          x2(i,j) = node_positions_x2(j) !x2_min+real(j-1,f32)*dx2
!        end do
!      end do
!      call sll_hdf5_file_create("cartesian_mesh-x1.h5",file_id,error)
!      call sll_hdf5_write_array(file_id,x1,"/x1",error)
!      call sll_hdf5_file_close(file_id, error)
!      call sll_hdf5_file_create("cartesian_mesh-x2.h5",file_id,error)
!      call sll_hdf5_write_array(file_id,x2,"/x2",error)
!      call sll_hdf5_file_close(file_id, error)
!      deallocate(x1)
!      deallocate(x2)
!
!    end if
!
!    call int2string(iplot,cplot)
!    call sll_xdmf_open(trim(array_name)//cplot//".xmf","cartesian_mesh", &
!      nnodes_x1,nnodes_x2,file_id,error)
!    write(file_id,"(a,f8.3,a)") "<Time Value='",time,"'/>"
!    call sll_xdmf_write_array(trim(array_name)//cplot,f,"values", &
!      error,file_id,"Node")
!    call sll_xdmf_close(file_id,error)
!  end subroutine plot_f_cartesian
!
!#endif
!  




end program