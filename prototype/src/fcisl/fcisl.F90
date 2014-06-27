module sll_fcisl_module
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
  use sll_logical_meshes
  use sll_module_advection_1d_base
implicit none

! a direction is chosen for interpolation
! derivative computation and advection

  type :: sll_oblic_derivative
    sll_int32 :: degree
    sll_real64, dimension(:,:), pointer :: buf
    type(sll_logical_mesh_1d), pointer :: mesh_x1
    type(sll_logical_mesh_1d), pointer :: mesh_x2
    class(sll_advection_1d_base), pointer :: adv
    sll_real64, dimension(:), pointer :: w
  end type sll_oblic_derivative

contains

  subroutine compute_oblic_shift(iota,Nc_x1,shift, iota_modif)
    sll_real64, intent(in) :: iota
    sll_int32, intent(in) :: Nc_x1
    sll_int32, intent(out) :: shift
    sll_real64, intent(out) :: iota_modif
    
    shift = floor(iota*Nc_x1+0.5)
    iota_modif = real(shift,f64)/real(Nc_x1,f64)
    
  end subroutine compute_oblic_shift

  subroutine compute_iota_from_shift(Nc_x1,shift, iota_modif)
    sll_int32, intent(in) :: Nc_x1
    sll_int32, intent(out) :: shift
    sll_real64, intent(out) :: iota_modif
    
    iota_modif = real(shift,f64)/real(Nc_x1,f64)
    
  end subroutine compute_iota_from_shift
  
  
  subroutine compute_at_aligned( &
    f_input, &
    f_output, &
    Nc_x1, &
    Nc_x2, &
    adv, &
    x1_min, &
    x1_max, &
    iota )
    sll_real64, dimension(:,:), intent(in) :: f_input
    sll_real64, dimension(:,:), intent(out) :: f_output
    sll_int32, intent(in) :: Nc_x1
    sll_int32, intent(in) :: Nc_x2
    class(sll_advection_1d_base), pointer :: adv
    sll_real64, intent(in) :: x1_min
    sll_real64, intent(in) :: x1_max
    sll_real64, intent(in) :: iota
    !local variables
    sll_int32 :: i
    sll_real64 :: A
    sll_real64 :: delta_x1
    A = iota*real(Nc_x1,f64)/real(Nc_x2,f64)
    delta_x1 = (x1_max-x1_min)/real(Nc_x1,f64)
    
    do i=1,Nc_x2+1
      call adv%advect_1d_constant( &
        A, &
        real(i-1,f64)*delta_x1, &
        f_input(1:Nc_x1+1,i), &
        f_output(1:Nc_x1+1,i))      
    enddo
  end subroutine compute_at_aligned
  
  subroutine compute_spaghetti_size_from_shift( &
    Nc_x1, &
    shift, &
    spaghetti_size)
    sll_int32, intent(in) :: Nc_x1
    sll_int32, intent(in) :: shift
    sll_int32, intent(out) :: spaghetti_size
    !local variables
    sll_int32 :: i
    sll_int32 :: i_loc
    
    spaghetti_size = 1
    do i = 1,Nc_x1
      i_loc = modulo(i*shift,Nc_x1)
      if(i_loc .ne. 0)then
        spaghetti_size = spaghetti_size+1
      else
        exit
      endif  
    enddo
  end subroutine compute_spaghetti_size_from_shift


!< compute shift and spaghetti_size
!< from iota_guess and spaghetti_size_guess
!< first search the existing spaghetti_size nearest of spaghetti_size_guess
!< for a shift between -Nc_x1 to Nc_x1
!< then looks for a shift that is near shift guess that has the same spaghetti_size   
  subroutine compute_spaghetti_and_shift_from_guess( &
    Nc_x1, &
    Nc_x2, &
    iota_guess, &
    spaghetti_size_guess, &
    shift, &
    spaghetti_size)
    sll_int32, intent(in) :: Nc_x1
    sll_int32, intent(in) :: Nc_x2
    !sll_int32, intent(in) :: shift_guess
    sll_real64, intent(in) :: iota_guess
    sll_int32, intent(in) :: spaghetti_size_guess
    sll_int32, intent(out) :: shift
    sll_int32, intent(out) :: spaghetti_size
    !local variables
    sll_int32 :: shift_guess
    sll_int32 :: shift_plus
    sll_int32 :: shift_minus
    sll_int32 :: i
    sll_int32 :: s
    sll_int32 :: i_val    
    sll_int32 :: val
    sll_int32 :: spaghetti_size_new
    s = 0
    val = 2*Nc_x1
    do i=-Nc_x1,Nc_x1
      call compute_spaghetti_size_from_shift( &
        Nc_x1, &
        i, &
        spaghetti_size)
      if (abs(spaghetti_size_guess-spaghetti_size) < abs(val))then
        val = spaghetti_size_guess-spaghetti_size
        i_val = i
        if(val == 0)then
          exit
        endif
      endif
    enddo
    
    call compute_spaghetti_size_from_shift( &
      Nc_x1, &
      i_val, &
      spaghetti_size)
    
    shift_guess = floor(iota_guess*Nc_x1)

    shift_plus = shift_guess
    shift_minus = shift_guess
    
    
    do i=0,Nc_x1
      call compute_spaghetti_size_from_shift( &
        Nc_x1, &
        shift_guess+i, &
        spaghetti_size_new)
      if(spaghetti_size_new .eq. spaghetti_size) then
        shift_plus = shift_guess+i
        exit
      endif
    enddo
    do i=0,Nc_x1              
      call compute_spaghetti_size_from_shift( &
        Nc_x1, &
        shift_guess-i, &
        spaghetti_size_new)      
      if(spaghetti_size_new .eq. spaghetti_size) then
        shift_minus = shift_guess-i
        exit
      endif        
    enddo   
    
    if(abs(shift_minus*Nc_x1-iota_guess)<abs(shift_plus*Nc_x1-iota_guess))then
      shift = shift_minus
    else
      shift = shift_plus  
    endif

    
    call compute_spaghetti_size_from_shift( &
      Nc_x1, &
      shift, &
      spaghetti_size)


!    call compute_spaghetti_size_from_shift( &
!      Nc_x1, &
!      shift_guess, &
!      spaghetti_size)
    
    
    
  end subroutine compute_spaghetti_and_shift_from_guess  
  
  subroutine compute_spaghetti( &
    Nc_x1, &
    Nc_x2, &
    shift, &
    spaghetti_index, &
    spaghetti_size)
    sll_int32, intent(in) :: Nc_x1
    sll_int32, intent(in) :: Nc_x2
    sll_int32, intent(in) :: shift
    sll_int32, dimension(:), intent(out) :: spaghetti_index
    sll_int32, intent(out) :: spaghetti_size
    !local variables
    sll_int32 :: i
    sll_int32 :: i_loc
    sll_int32 :: ii
    sll_int32 :: j
    sll_int32 :: q
    sll_int32, dimension(:), allocatable ::  check
    sll_int32 :: ierr
    
    SLL_ALLOCATE(check(Nc_x1+1),ierr)
    !0, shift,2*shift, Nc_x1*shift
    !while(modulo(k*shift,Nc_x1) .ne. 0)
    
    spaghetti_index = 0
    spaghetti_index(1) = 1
    spaghetti_size = 1
    do i = 1,Nc_x1
      i_loc = modulo(i*shift,Nc_x1)
      if(i_loc .ne. 0)then
        spaghetti_size = spaghetti_size+1
        spaghetti_index(i+1) = i_loc+1
      else
        exit
      endif  
    enddo
    !print *,'#spaghetti_size=',shift,spaghetti_size
    !print *,'#i,i_loc=',i,i_loc
    q = Nc_x1/spaghetti_size
    do j=1,q-1
      do ii=1,spaghetti_size
      spaghetti_index(ii+j*spaghetti_size) = modulo(spaghetti_index(ii)+j-1,Nc_x1)+1
      enddo
    enddo
    spaghetti_index(Nc_x1+1) = spaghetti_index(1)
    check = 0
    do i=1,Nc_x1
      if(spaghetti_index(i)<1)then
        print *,'#problem with spaghetti_index'
        stop
      endif
      if(spaghetti_index(i)>Nc_x1+1)then
        print *,'#problem with spaghetti_index'
        stop
      endif
      check(spaghetti_index(i)) = check(spaghetti_index(i))+1
      !print *,i,spaghetti_index(i)
    enddo
    if(maxval(check(1:Nc_x1)) .ne. 1)then
      print *,'#problem checking spaghetti_index'
      print *,'#maxval=',maxval(check(1:Nc_x1))
      stop
    endif
    if(minval(check(1:Nc_x1)) .ne. 1)then
      print *,'#problem checking spaghetti_index'
      print *,'#minval=',minval(check(1:Nc_x1))
      stop
    endif
    
  end subroutine compute_spaghetti

  subroutine load_spaghetti( &
    input, &
    output, &
    spaghetti_index, &
!    spaghetti_size, &
    Npts_x1, &
    Npts_x2)
    sll_real64, dimension(:,:), intent(in) :: input
    sll_real64, dimension(:), intent(out) :: output
    sll_int32, dimension(:), intent(in) :: spaghetti_index
!    sll_int32, intent(in) :: spaghetti_size
    sll_int32, intent(in) :: Npts_x1
    sll_int32, intent(in) :: Npts_x2
    !local variables
    sll_int32 :: i
    sll_int32 :: j
    sll_int32 :: s
    s = 0
    do i=1,Npts_x1
      do j=1,Npts_x2
        s = s+1
        output(s) = input(spaghetti_index(i),j)
      enddo  
    enddo
  end subroutine load_spaghetti

  subroutine unload_spaghetti( &
    input, &
    output, &
    spaghetti_index, &
!    spaghetti_size, &
    Npts_x1, &
    Npts_x2)
    sll_real64, dimension(:), intent(in) :: input
    sll_real64, dimension(:,:), intent(out) :: output
    sll_int32, dimension(:), intent(in) :: spaghetti_index
!    sll_int32, intent(in) :: spaghetti_size
    sll_int32, intent(in) :: Npts_x1
    sll_int32, intent(in) :: Npts_x2
    !local variables
    sll_int32 :: i
    sll_int32 :: j
    sll_int32 :: s
    s = 0
    do i=1,Npts_x1
      do j=1,Npts_x2
        s = s+1
        output(spaghetti_index(i),j) = input(s) 
      enddo  
    enddo
  end subroutine unload_spaghetti
  
  

  
  function new_oblic_derivative( &
    degree, &
    x1_min, &
    x1_max, &
    x2_min, &
    x2_max, &
    Nc_x1, &
    Nc_x2, &
    adv &
    ) result(res)
    type(sll_oblic_derivative), pointer :: res
    sll_int32, intent(in) :: degree
    sll_real64, intent(in) :: x1_min
    sll_real64, intent(in) :: x1_max
    sll_real64, intent(in) :: x2_min
    sll_real64, intent(in) :: x2_max
    sll_int32, intent(in) :: Nc_x1
    sll_int32, intent(in) :: Nc_x2
    class(sll_advection_1d_base), pointer :: adv
    !local variables
    sll_int32 :: ierr
    SLL_ALLOCATE(res,ierr)    
    call initialize_oblic_derivative( &
      res, &
      degree, &
      x1_min, &
      x1_max, &
      x2_min, &
      x2_max, &
      Nc_x1, &
      Nc_x2, &
      adv)
    
  end function new_oblic_derivative
  
  subroutine initialize_oblic_derivative( &
    deriv, &
    degree, &
    x1_min, &
    x1_max, &
    x2_min, &
    x2_max, &
    Nc_x1, &
    Nc_x2, &
    adv)
    type(sll_oblic_derivative) :: deriv
    sll_int32, intent(in) :: degree
    sll_real64, intent(in) :: x1_min
    sll_real64, intent(in) :: x1_max
    sll_real64, intent(in) :: x2_min
    sll_real64, intent(in) :: x2_max
    sll_int32, intent(in) :: Nc_x1
    sll_int32, intent(in) :: Nc_x2
    class(sll_advection_1d_base), pointer :: adv
    sll_int32 :: ierr
    
    deriv%degree = degree    
    deriv%mesh_x1 => new_logical_mesh_1d(Nc_x1,eta_min=x1_min,eta_max=x1_max)
    deriv%mesh_x2 => new_logical_mesh_1d(Nc_x2,eta_min=x2_min,eta_max=x2_max)
    deriv%adv => adv
    
    SLL_ALLOCATE(deriv%buf(1:Nc_x2+2*degree+1,1:Nc_x1+1),ierr) 
    SLL_ALLOCATE(deriv%w(-degree:Nc_x2+degree),ierr) 
    
    call compute_finite_difference_init(deriv%w,2*degree)
      
  end subroutine initialize_oblic_derivative


!> compute b_tau \cdot \nabla phi
!> with b_tau = (tau/sqrt(1+tau^2))*hat_theta+ (1/sqrt(1+tau^2))*hat_z
  subroutine compute_oblic_derivative( &
    deriv, &
    tau, &
    phi, &
    D_phi)
    type(sll_oblic_derivative), pointer :: deriv
    sll_real64, intent(in) :: tau
    sll_real64, dimension(:,:), intent(in) :: phi
    sll_real64, dimension(:,:), intent(out) :: D_phi
    class(sll_advection_1d_base), pointer :: adv
    !local variables
    sll_real64, dimension(:,:), pointer :: buf
    sll_int32 :: step
    sll_int32 :: Nc_x1
    sll_int32 :: Nc_x2
    sll_real64 :: delta_x2
    sll_real64, dimension(:), allocatable :: w 
    sll_int32 :: ierr
    sll_int32 :: i2_loc
    sll_int32 :: d
    
    Nc_x1 = deriv%mesh_x1%num_cells
    Nc_x2 = deriv%mesh_x2%num_cells
    d = deriv%degree
    
    
    
    buf => deriv%buf
    
    !step 1: compute phi on a field aligned mesh from the initial field aligned mesh
    ! we store on phi_store(-d:Nc_x2+d,1:Nc_x1+1)
   
    delta_x2 = deriv%mesh_x2%delta_eta
    
    if(size(phi,1)<Nc_x1+1)then
      print *,'#bad size for phi'
    endif
    if(size(phi,2)<Nc_x2+1)then
      print *,'#bad size for phi'
    endif
    
    print *,size(buf,1),size(buf,2)
       
   
    do step = -d,Nc_x2+d
      i2_loc = modulo(step+Nc_x2,Nc_x2)+1
      !deriv%buf(step,1:Nc_x1+1) = phi(1:Nc_x1+1,i2_loc)
      print *,i2_loc
      print *,real(step,f64)*delta_x2
      print *,maxval(phi(1:Nc_x1+1,i2_loc)),minval(phi(1:Nc_x1+1,i2_loc))
      print *,maxval(deriv%buf(step+d+1,1:Nc_x1+1)),minval(deriv%buf(step+d+1,1:Nc_x1+1))
      print *,tau
      print *,real(step,f64)*delta_x2
      call adv%advect_1d_constant( &
        tau, &
        real(step,f64)*delta_x2, &
        phi(1:Nc_x1+1,i2_loc), &
        deriv%buf(step+d+1,1:Nc_x1+1))      
    enddo
   
   
   
   
    !step 2: compute derivative on the field aligned mesh
   
    !step 3: compute derivative on initial cartesian mesh from field on field aligned mesh
   
   
   
  end subroutine compute_oblic_derivative


  subroutine compute_finite_difference_init(w,p)
    integer,intent(in)::p    
    real(f64),dimension(-p/2:(p+1)/2),intent(out)::w
    integer::r,s,i,j
    real(f64)::tmp

    r=-p/2
    s=(p+1)/2

!    if(modulo(p,2)==0)then
!      r=-p/2
!      s=p/2
!    else
!      r=-(p-1)/2
!      s=(p+1)/2
!    endif
        
    
    !maple code for generation of w
    !for k from r to -1 do
    !  C[k]:=product((k-j),j=r..k-1)*product((k-j),j=k+1..s):
    !  C[k]:=1/C[k]*product((-j),j=r..k-1)*product((-j),j=k+1..-1)*product((-j),j=1..s):
    !od:
    !for k from 1 to s do
    !  C[k]:=product((k-j),j=r..k-1)*product((k-j),j=k+1..s):
    !  C[k]:=1/C[k]*product((-j),j=r..-1)*product((-j),j=1..k-1)*product((-j),j=k+1..s):
    !od:
    !C[0]:=-add(C[k],k=r..-1)-add(C[k],k=1..s):
    
    do i=r,-1
      tmp=1._f64
      do j=r,i-1
        tmp=tmp*real(i-j,f64)
      enddo
      do j=i+1,s
        tmp=tmp*real(i-j,f64)
      enddo
      tmp=1._f64/tmp
      do j=r,i-1
        tmp=tmp*real(-j,f64)
      enddo
      do j=i+1,-1
        tmp=tmp*real(-j,f64)
      enddo
      do j=1,s
        tmp=tmp*real(-j,f64)
      enddo
      w(i)=tmp      
    enddo

    do i=1,s
      tmp=1._f64
      do j=r,i-1
        tmp=tmp*real(i-j,f64)
      enddo
      do j=i+1,s
        tmp=tmp*real(i-j,f64)
      enddo
      tmp=1._f64/tmp
      do j=r,-1
        tmp=tmp*real(-j,f64)
      enddo
      do j=1,i-1
        tmp=tmp*real(-j,f64)
      enddo
      do j=i+1,s
        tmp=tmp*real(-j,f64)
      enddo
      w(i)=tmp      
    enddo
    tmp=0._f64
    do i=r,-1
      tmp=tmp+w(i)
    enddo
    do i=1,s
      tmp=tmp+w(i)
    enddo
    w(0)=-tmp
    
  end subroutine compute_finite_difference_init



!  subroutine compute_field_from_phi_cartesian_1d(phi,mesh,A,interp)
!    sll_real64, dimension(:), intent(in) :: phi
!    sll_real64, dimension(:), intent(out) :: A
!    type(sll_logical_mesh_1d), pointer :: mesh
!    class(sll_interpolator_1d_base), pointer   :: interp
!    sll_int32 :: Nc_x1
!    sll_real64 :: x1_min
!    sll_real64 :: delta_x1
!    sll_real64 :: x1
!    sll_int32 :: i1
!    
!    Nc_x1 = mesh%num_cells
!    x1_min = mesh%eta_min
!    delta_x1 = mesh%delta_eta
!
!    call interp%compute_interpolants(phi)
!
!    do i1=1,Nc_x1+1
!      x1=x1_min+real(i1-1,f64)*delta_x1
!      A(i1)=interp%interpolate_derivative_eta1(x1)
!    end do
!  end subroutine compute_field_from_phi_cartesian_1d





end module sll_fcisl_module
