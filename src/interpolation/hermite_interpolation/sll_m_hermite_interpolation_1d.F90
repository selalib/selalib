module sll_m_hermite_interpolation_1d
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
  use sll_m_boundary_condition_descriptors
implicit none

!Hermite interpolation in 1d
!derivatives are given with finite stencil formulae of order p
!which can be arbitrary in each direction
!If p is odd, the reconstruction has discontinuous derivatives
!If p is even, the reconstruction has continuous derivatives
!p=6 should be quite similar to cubic splines
!do not hesitate to take large odd p, like the favourite p=17
!see also sll_m_hermite_interpolation_2d
  type :: sll_hermite_interpolation_1d
    sll_real64 :: eta_min !< eta1 min
    sll_real64 :: eta_max !< eta1 max
    sll_int32 :: Nc !< number of cells in eta1     
    sll_int32 :: degree !< interpolation degrees in eta1, eta2
    sll_int32 :: bc !<boundary conditions in eta1
                       !< can be SLL_HERMITE_PERIODIC or SLL_HERMITE_DIRICHLET
    sll_real64, dimension(:,:), pointer :: deriv
    sll_int32 :: continuity !< can be SLL_HERMITE_C0 or SLL_HERMITE_C1
    sll_int32 :: deriv_size !< 3 for SLL_HERMITE_C0 and 2 for SLL_HERMITE_C1
                            !< deriv_size = size(deriv,1)                              
  end type sll_hermite_interpolation_1d 
  
  !integer, parameter :: SLL_HERMITE_PERIODIC = 0, SLL_HERMITE_DIRICHLET = 1
  integer, parameter :: SLL_HERMITE_1d_C0 = 0, SLL_HERMITE_1d_C1 = 1

!  interface delete
!    module procedure delete_hermite_interpolation_1d
!  end interface

contains  !*****************************************************************************
  function new_hermite_interpolation_1d( &
    npts, &
    eta_min, &
    eta_max, &
    degree, &
    eta_hermite_continuity, &
    eta_bc_type, &
    const_eta_min_slope, &
    const_eta_max_slope, &
    eta_min_slopes, &
    eta_max_slopes) &
    result(interp)
    
    type(sll_hermite_interpolation_1d), pointer :: interp
    sll_int32, intent(in) :: npts
    sll_real64, intent(in) :: eta_min
    sll_real64, intent(in) :: eta_max
    sll_int32, intent(in) :: degree
    sll_int32, intent(in) :: eta_hermite_continuity
    sll_int32, intent(in) :: eta_bc_type
    sll_real64, intent(in), optional :: const_eta_min_slope
    sll_real64, intent(in), optional :: const_eta_max_slope
    sll_real64, dimension(:),intent(in), optional :: eta_min_slopes
    sll_real64, dimension(:),intent(in), optional :: eta_max_slopes
    sll_int32 :: ierr
    
    SLL_ALLOCATE(interp,ierr)
    
    call initialize_hermite_interpolation_1d( &
      interp, &
      npts, &
      eta_min, &
      eta_max, &
      degree, &
      eta_hermite_continuity, &
      eta_bc_type,   &
      const_eta_min_slope, &
      const_eta_max_slope, &
      eta_min_slopes, &
      eta_max_slopes)

  end function new_hermite_interpolation_1d

  subroutine initialize_hermite_interpolation_1d( &
    interp, &
    npts, &
    eta_min, &
    eta_max, &
    degree, &
    eta_hermite_continuity, &
    eta_bc_type,   &
    const_eta_min_slope, &
    const_eta_max_slope, &
    eta_min_slopes, &
    eta_max_slopes)

    type(sll_hermite_interpolation_1d) :: interp

    sll_int32, intent(in) :: npts
    sll_real64, intent(in) :: eta_min
    sll_real64, intent(in) :: eta_max
    sll_int32, intent(in) :: degree
    sll_int32, intent(in) :: eta_hermite_continuity
    sll_int32, intent(in) :: eta_bc_type
    sll_real64, intent(in), optional :: const_eta_min_slope
    sll_real64, intent(in), optional :: const_eta_max_slope
    sll_real64, dimension(:),intent(in), optional :: eta_min_slopes
    sll_real64, dimension(:),intent(in), optional :: eta_max_slopes
    sll_int32 :: ierr
    sll_int32 :: deriv_size
    !sll_int32 :: i
    
    interp%Nc = npts-1
    interp%degree = degree    
    interp%continuity = eta_hermite_continuity
    interp%bc = eta_bc_type
    interp%eta_min = eta_min
    interp%eta_max = eta_max
    
      select case (interp%continuity)
        case (SLL_HERMITE_1d_C0)
          interp%deriv_size = 3
        case (SLL_HERMITE_1d_C1)
          interp%deriv_size = 2
        case default
          print *,'#bad value for hermite_continuity',interp%continuity
          print *,'#in initialize_hermite_interpolation_1d'
          stop       
          SLL_ASSERT(present(const_eta_min_slope))
          SLL_ASSERT(present(const_eta_max_slope))
          SLL_ASSERT(present(eta_min_slopes))
          SLL_ASSERT(present(eta_max_slopes))
      end select
      deriv_size = interp%deriv_size
 

    SLL_ALLOCATE(interp%deriv(deriv_size,npts),ierr)
        
    
  end subroutine initialize_hermite_interpolation_1d

  subroutine compute_interpolants_hermite_1d( &
    interp, &
    f)
    type(sll_hermite_interpolation_1d) :: interp
    sll_real64, dimension(:), intent(in) :: f
    
    
    if((interp%bc==SLL_HERMITE)) then
      if(interp%continuity==SLL_HERMITE_1d_C0) then         
        call hermite_coef_nat_1d(f,interp%deriv,interp%Nc,interp%degree)
      else
        print *,'#interp%continuity=', interp%continuity
        print *,'#possible_value=', SLL_HERMITE_1d_C0
        print *,'#not implemented for the moment'
        print *,'#in compute_interpolants_hermite_1d'
        stop          
      endif
    else if((interp%bc==SLL_PERIODIC)) then
      if(interp%continuity==SLL_HERMITE_1d_C0) then         
        call hermite_coef_per_1d(f,interp%deriv,interp%Nc,interp%degree)
      else
        print *,'#interp%continuity=', interp%continuity
        print *,'#possible_value=', SLL_HERMITE_1d_C0
        print *,'#not implemented for the moment'
        print *,'#in compute_interpolants_hermite_1d'
        stop          
      endif
    
    else  
       
      print *,'#interp%bc=', interp%bc
      print *,'#possible_value=', SLL_HERMITE, SLL_PERIODIC     
      print *,'#not implemented for the moment'
      print *,'#in compute_interpolants_hermite_1d'
      stop          
    endif
    
    
  end subroutine compute_interpolants_hermite_1d



subroutine compute_w_hermite_1d(w,r,s)
    sll_int32,intent(in)::r,s
    sll_real64,dimension(r:s),intent(out)::w
    sll_int32 ::i,j
    sll_real64::tmp

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

    !print *,'w',w
    !do ii=r,s
    !  print *,ii,w(r+s-ii)
    !enddo
    !

  
  end subroutine compute_w_hermite_1d




subroutine hermite_coef_per_1d(f,buf2d,N,d)
    sll_int32,intent(in)::N,d
    sll_real64,dimension(N+1),intent(in)::f
    sll_real64,dimension(3,N+1),intent(out)::buf2d
    sll_real64 ::w_left(-d/2:(d+1)/2),w_right((-d+1)/2:d/2+1)
    sll_real64 ::tmp
    sll_int32  ::i,ii,r_left,r_right,s_left,s_right,ind 
    r_left=-d/2
    s_left=(d+1)/2
    r_right=(-d+1)/2
    s_right=d/2+1
    
    
    
    call compute_w_hermite_1d(w_left,r_left,s_left)
    if((2*(d/2)-d)==0)then
      w_right(r_right:s_right) = w_left(r_left:s_left)
    else
      w_right(r_right:s_right) = -w_left(s_left:r_left:-1)
    endif    

    
      do i=1,N+1
        buf2d(1,i)=f(modulo(i-1+N,N)+1) !f(0)
        tmp=0._f64
        do ii=r_left,s_left
          ind=modulo(i+ii-1+N,N)+1
          tmp=tmp+w_left(ii)*f(ind)
        enddo
        buf2d(2,i)=tmp !fx(0)
        tmp=0._f64
        do ii=r_right,s_right
          ind=modulo(i+ii-1+N,N)+1          
          tmp=tmp+w_right(ii)*f(ind)
        enddo
        buf2d(3,i)=tmp !fx(1)       
      enddo
      
      !print *,buf2d(1,N+1)-buf2d(1,1)
      !print *,buf2d(2,N+1)-buf2d(2,1)
      !print *,buf2d(3,N+1)-buf2d(3,1)
      
  end subroutine hermite_coef_per_1d


subroutine hermite_coef_nat_1d(f,buf2d,N,d)
    sll_int32,intent(in)::N,d
    sll_real64,dimension(N+1),intent(in)::f
    sll_real64,dimension(3,N+1),intent(out)::buf2d
    sll_real64 ::w_left(-d/2:(d+1)/2),w_right((-d+1)/2:d/2+1)
    sll_real64 ::tmp
    sll_int32  ::i,ii,r_left,r_right,s_left,s_right,ind 
    r_left=-d/2
    s_left=(d+1)/2
    r_right=(-d+1)/2
    s_right=d/2+1
    
    
    call compute_w_hermite_1d(w_left,r_left,s_left)
    if((2*(d/2)-d)==0)then
      w_right(r_right:s_right) = w_left(r_left:s_left)
    else
      w_right(r_right:s_right) = -w_left(s_left:r_left:-1)
    endif    

    
      do i=1,N+1
        buf2d(1,i)=f(i) !f(0)
        tmp=0._f64
        do ii=r_left,s_left
          !ind=modulo(i+ii-1+N,N)+1
          ind=i+ii;if(ind>N+1) ind=2*(N+1)-ind;if(ind<1) ind=2-ind
          tmp=tmp+w_left(ii)*f(ind)
        enddo
        buf2d(2,i)=tmp !fx(0)
        tmp=0._f64
        do ii=r_right,s_right
          !ind=modulo(i+ii-1+N,N)+1          
          ind=i+ii;if(ind>N+1) ind=2*(N+1)-ind;if(ind<1) ind=2-ind
          tmp=tmp+w_right(ii)*f(ind)
        enddo
        buf2d(3,i)=tmp !fx(1)       
      enddo
    
  end subroutine hermite_coef_nat_1d




  function interpolate_value_hermite_per_1d( eta, interp ) result(res)
    sll_real64,intent(in) :: eta
    type(sll_hermite_interpolation_1d), pointer :: interp
    sll_real64 :: res
    sll_real64 :: eta_tmp
    sll_int32 :: ii
    
    eta_tmp = eta
    
    call localize_per_1d(ii,eta_tmp,interp%eta_min,interp%eta_max,interp%Nc)
    call interpolate_hermite_1d(interp%deriv,ii,eta_tmp,res,interp%Nc)
    
  end function interpolate_value_hermite_per_1d

  function interpolate_value_hermite_1d( eta, interp ) result(res)
    sll_real64,intent(in) :: eta
    type(sll_hermite_interpolation_1d), pointer :: interp
    sll_real64 :: res
    sll_real64 :: eta_tmp
    sll_int32 :: ii
    
    eta_tmp = eta
    !warning: we should implement localize_1d
    !but as eta should be inside localize_nat
    !or localize_per can be used in principle
    !call localize_nat_1d(ii,eta_tmp,interp%eta_min,interp%eta_max,interp%Nc)
    call localize_per_1d(ii,eta_tmp,interp%eta_min,interp%eta_max,interp%Nc)
    call interpolate_hermite_1d(interp%deriv,ii,eta_tmp,res,interp%Nc)
    
    
  end function interpolate_value_hermite_1d





  function interpolate_value_hermite_nat_1d( eta, interp ) result(res)
    sll_real64,intent(in) :: eta
    type(sll_hermite_interpolation_1d), pointer :: interp
    sll_real64 :: res
    sll_real64 :: eta_tmp
    sll_int32 :: ii
    
    eta_tmp = eta
        
    call localize_nat_1d(ii,eta_tmp,interp%eta_min,interp%eta_max,interp%Nc)
    call interpolate_hermite_1d(interp%deriv,ii,eta_tmp,res,interp%Nc)
    
  end function interpolate_value_hermite_nat_1d


 subroutine interpolate_hermite_1d(f,i,x,fval,N)
    sll_int32,intent(in)::i,N
    !real(f64),intent(in)::xx(2),xmin(2),xmax(2)
    sll_real64,intent(in)::x
    sll_real64,intent(out)::fval
    sll_real64,dimension(0:2,0:N)::f
    !integer::i(2),i1(2),s
    sll_int32::i1!,s
    sll_real64::w(0:3)!,tmp(0:3)
    sll_real64::g(0:3)
    
    !fval =f(0,i(1),i(2))!real(i(1),f64)
    
    !return
    

      w(0)=(2._f64*x+1)*(1._f64-x)*(1._f64-x);
      w(1)=x*x*(3._f64-2._f64*x)
      w(2)=x*(1._f64-x)*(1._f64-x)
      w(3)=x*x*(x-1._f64)
      i1=i+1


    
    g(0)=f(0,i)          !f(0)
    g(1)=f(0,i1)         !f(1)
    g(2)=f(1,i)          !fx(0)
    g(3)=f(2,i)          !fx(1)

    !print *,'#in hermite:'
    !print *,'i=',i
    !print *,'x=',x
    !print *,'w=',w
    !print *,'g=',g


    fval=w(0)*g(0)+w(1)*g(1)+w(2)*g(2)+w(3)*g(3)
  
    !print *,fval !,' t',f    
  end subroutine interpolate_hermite_1d

  subroutine localize_per_1d(i,x,xmin,xmax,N)
    sll_int32,intent(out)::i
    sll_real64,intent(inout)::x
    sll_real64,intent(in)::xmin,xmax
    sll_int32,intent(in)::N
    x=(x-xmin)/(xmax-xmin)
    x=x-real(floor(x),f64)
    x=x*real(N,f64)
    i=floor(x)
    x=x-real(i,f64)
    if(i==N)then
      i=0
      x=0._f64
    endif
  end subroutine localize_per_1d
  

  subroutine localize_nat_1d(i,x,xmin,xmax,N)
    sll_int32,intent(out)::i
    sll_real64,intent(inout)::x
    sll_real64,intent(in)::xmin,xmax
    sll_int32,intent(in)::N
    x=(x-xmin)/(xmax-xmin)
    x=x*real(N,f64)
    if(x>=real(N,f64))then
      x=real(N,f64)
    endif
    if(x<=0._f64)then
      x=0._f64
    endif    
    i=floor(x)
    x=x-real(i,f64)
    if(i==N)then
      i=N-1
      x=1._f64
    endif
  end subroutine localize_nat_1d


  
    

end module
