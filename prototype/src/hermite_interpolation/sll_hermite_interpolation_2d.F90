module sll_hermite_interpolation_2d_module
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
implicit none

!Hermite interpolation in 2d
!derivatives are given with finite stencil formulae of order p
!which can be arbitrary in each direction
!If p is odd, the reconstruction has discontinuous derivatives
!If p is even, the reconstruction has continuous derivatives
!p=6 should be quite similar to cubic splines
!do not hesitate to take large odd p, like the favourite p=17
!! WARNING
! for the moment only in implementation for the case DIRICHLET x PERIODIC


  type :: sll_hermite_interpolation_2d
    sll_real64 :: eta_min(2) !< eta1 min and eta2 min
    sll_real64 :: eta_max(2) !< eta1 max and eta2 max
    sll_int32 :: Nc(2) !< number of cells in eta1 and in eta2     
    sll_int32 :: degree(2) !< interpolation degrees in eta1, eta2
    sll_int32 :: bc(2) !<boundary conditions in eta and eta2
                       !< can be SLL_HERMITE_PERIODIC or SLL_HERMITE_DIRICHLET
    sll_real64, dimension(:,:,:), pointer :: deriv
    sll_int32 :: continuity(2) !< can be SLL_HERMITE_C0 or SLL_HERMITE_C1
    sll_int32 :: deriv_size(2) !< 3 for SLL_HERMITE_C0 and 2 for SLL_HERMITE_C1
                               !< in each direction
                               !< deriv_size(1)*deriv_size(2) = size(deriv,1)                              
  end type sll_hermite_interpolation_2d 
  
  integer, parameter :: SLL_HERMITE_PERIODIC = 0, SLL_HERMITE_DIRICHLET = 1
  integer, parameter :: SLL_HERMITE_C0 = 0, SLL_HERMITE_C1 = 1

!  interface delete
!    module procedure delete_hermite_interpolation_2d
!  end interface

contains  !*****************************************************************************
  function new_hermite_interpolation_2d( &
    npts1, &
    npts2, &
    eta1_min, &
    eta1_max, &
    eta2_min, &
    eta2_max, &
    degree1, &
    degree2, &
    eta1_hermite_continuity, &
    eta2_hermite_continuity, &
    eta1_bc_type, &
    eta2_bc_type, &
    const_eta1_min_slope, &
    const_eta1_max_slope, &
    const_eta2_min_slope, &
    const_eta2_max_slope, &
    eta1_min_slopes, &
    eta1_max_slopes, &
    eta2_min_slopes, &
    eta2_max_slopes ) &
    result(interp)
    
    type(sll_hermite_interpolation_2d), pointer :: interp
    sll_int32, intent(in) :: npts1
    sll_int32, intent(in) :: npts2
    sll_real64, intent(in) :: eta1_min
    sll_real64, intent(in) :: eta1_max
    sll_real64, intent(in) :: eta2_min
    sll_real64, intent(in) :: eta2_max
    sll_int32, intent(in) :: degree1
    sll_int32, intent(in) :: degree2    
    sll_int32, intent(in) :: eta1_hermite_continuity
    sll_int32, intent(in) :: eta2_hermite_continuity
    sll_int32, intent(in) :: eta1_bc_type
    sll_int32, intent(in) :: eta2_bc_type
    sll_real64, intent(in), optional :: const_eta1_min_slope
    sll_real64, intent(in), optional :: const_eta1_max_slope
    sll_real64, intent(in), optional :: const_eta2_min_slope
    sll_real64, intent(in), optional :: const_eta2_max_slope
    sll_real64, dimension(:),intent(in), optional :: eta1_min_slopes
    sll_real64, dimension(:),intent(in), optional :: eta1_max_slopes
    sll_real64, dimension(:),intent(in), optional :: eta2_min_slopes
    sll_real64, dimension(:),intent(in), optional :: eta2_max_slopes
    sll_int32 :: ierr
    
    SLL_ALLOCATE(interp,ierr)
    
    call initialize_hermite_interpolation_2d( &
      interp, &
      npts1, &
      npts2, &
      eta1_min, &
      eta1_max, &
      eta2_min, &
      eta2_max, &
      degree1, &
      degree2, &
      eta1_hermite_continuity, &
      eta2_hermite_continuity, &
      eta1_bc_type,   &
      eta2_bc_type,   &
      const_eta1_min_slope, &
      const_eta1_max_slope, &
      const_eta2_min_slope, &
      const_eta2_max_slope, &
      eta1_min_slopes, &
      eta1_max_slopes, &
      eta2_min_slopes, &
      eta2_max_slopes )

  end function new_hermite_interpolation_2d

  subroutine initialize_hermite_interpolation_2d( &
    interp, &
    npts1, &
    npts2, &
    eta1_min, &
    eta1_max, &
    eta2_min, &
    eta2_max, &
    degree1, &
    degree2, &
    eta1_hermite_continuity, &
    eta2_hermite_continuity, &
    eta1_bc_type,   &
    eta2_bc_type,   &
    const_eta1_min_slope, &
    const_eta1_max_slope, &
    const_eta2_min_slope, &
    const_eta2_max_slope, &
    eta1_min_slopes, &
    eta1_max_slopes, &
    eta2_min_slopes, &
    eta2_max_slopes )
    type(sll_hermite_interpolation_2d) :: interp
    sll_int32, intent(in) :: npts1
    sll_int32, intent(in) :: npts2
    sll_real64, intent(in) :: eta1_min
    sll_real64, intent(in) :: eta1_max
    sll_real64, intent(in) :: eta2_min
    sll_real64, intent(in) :: eta2_max
    sll_int32, intent(in) :: degree1
    sll_int32, intent(in) :: degree2    
    sll_int32, intent(in) :: eta1_hermite_continuity
    sll_int32, intent(in) :: eta2_hermite_continuity
    sll_int32, intent(in) :: eta1_bc_type
    sll_int32, intent(in) :: eta2_bc_type
    sll_real64, intent(in), optional :: const_eta1_min_slope
    sll_real64, intent(in), optional :: const_eta1_max_slope
    sll_real64, intent(in), optional :: const_eta2_min_slope
    sll_real64, intent(in), optional :: const_eta2_max_slope
    sll_real64, dimension(:),intent(in), optional :: eta1_min_slopes
    sll_real64, dimension(:),intent(in), optional :: eta1_max_slopes
    sll_real64, dimension(:),intent(in), optional :: eta2_min_slopes
    sll_real64, dimension(:),intent(in), optional :: eta2_max_slopes
    sll_int32 :: ierr
    sll_int32 :: deriv_size
    sll_int32 :: i
    
    interp%Nc(1:2) = (/npts1-1,npts2-1/)
    interp%degree(1:2) = (/degree1,degree2/)    
    interp%continuity = (/eta1_hermite_continuity,eta2_hermite_continuity/)
    interp%bc(1:2) = (/eta1_bc_type,eta2_bc_type/)
    interp%eta_min = (/eta1_min,eta2_min/)
    interp%eta_max = (/eta1_max,eta2_max/)
    
    deriv_size = 1
    do i=1,2
      select case (interp%continuity(i))
        case (SLL_HERMITE_C0)
          interp%deriv_size(i) = 3
        case (SLL_HERMITE_C1)
          interp%deriv_size(i) = 2
        case default
          print *,'#bad value for hermite_continuity',interp%continuity
          print *,'#in initialize_hermite_interpolation_2d'
          stop       
      end select
      deriv_size = deriv_size*interp%deriv_size(i)
    enddo  

    SLL_ALLOCATE(interp%deriv(deriv_size,npts1,npts2),ierr)
        
    
  end subroutine initialize_hermite_interpolation_2d

  subroutine compute_interpolants_hermite_2d( &
    interp, &
    f)
    type(sll_hermite_interpolation_2d) :: interp
    sll_real64, dimension(:,:), intent(in) :: f
    if((interp%bc(1)==SLL_HERMITE_DIRICHLET).and.&
      (interp%bc(2)==SLL_HERMITE_PERIODIC)) then
      if((interp%continuity(1)==SLL_HERMITE_C0).and.&
        (interp%continuity(2)==SLL_HERMITE_C0)) then         
        call hermite_coef_nat_per(f,interp%deriv,interp%Nc,interp%degree)
      else
        print *,'#interp%continuity=', interp%continuity
        print *,'#possible_value=', SLL_HERMITE_C0
        print *,'#not implemented for the moment'
        print *,'#in compute_interpolants_hermite_2d'
        stop          
      endif
    else   
      print *,'#interp%bc=', interp%bc
      print *,'#possible_value=', SLL_HERMITE_DIRICHLET, SLL_HERMITE_PERIODIC     
      print *,'#not implemented for the moment'
      print *,'#in compute_interpolants_hermite_2d'
      stop          
    endif
    
    
  end subroutine compute_interpolants_hermite_2d



subroutine compute_w_hermite(w,r,s)
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

  
  end subroutine compute_w_hermite


subroutine hermite_coef_nat_per(f,buf3d,N,d)
    sll_int32,intent(in)::N(2),d(2)
    sll_real64,dimension(N(1)+1,N(2)),intent(in)::f
    sll_real64,dimension(9,N(1)+1,N(2)+1),intent(out)::buf3d
    sll_real64 ::w_left_1(-d(1)/2:(d(1)+1)/2),w_right_1((-d(1)+1)/2:d(1)/2+1)
    sll_real64 ::w_left_2(-d(2)/2:(d(2)+1)/2),w_right_2((-d(2)+1)/2:d(2)/2+1)
    sll_real64 ::tmp
    sll_int32  ::i,j,ii,r_left(2),r_right(2),s_left(2),s_right(2),ind 
    r_left=-d/2
    s_left=(d+1)/2
    r_right=(-d+1)/2
    s_right=d/2+1
    
    
    call compute_w_hermite(w_left_1,r_left(1),s_left(1))
    call compute_w_hermite(w_left_2,r_left(2),s_left(2))    
    if((2*(d(1)/2)-d(1))==0)then
      w_right_1(r_right(1):s_right(1)) = w_left_1(r_left(1):s_left(1))
    else
      w_right_1(r_right(1):s_right(1)) = -w_left_1(s_left(1):r_left(1):-1)
    endif    

    if((2*(d(2)/2)-d(2))==0)then
      w_right_2(r_right(2):s_right(2)) = w_left_2(r_left(2):s_left(2))
    else
      w_right_2(r_right(2):s_right(2)) = -w_left_2(s_left(2):r_left(2):-1)
    endif    
    
    !print *,'w(',r_left(1),':',s_left(1),')=',w_left_1(r_left(1):s_left(1))
    !print *,'w(',r_right(1),':',s_right(1),')=',w_right_1(r_right(1):s_right(1))

    
    do j=1,N(2)
      do i=1,N(1)+1
        buf3d(1,i,j)=f(i,j) !f(0,0)
        tmp=0._f64
        do ii=r_left(1),s_left(1)
          ind=i+ii;if(ind>N(1)+1) ind=2*(N(1)+1)-ind;if(ind<1) ind=2-ind
          tmp=tmp+w_left_1(ii)*f(ind,j)
        enddo
        buf3d(2,i,j)=tmp !fx(0,0)
        tmp=0._f64
        do ii=r_right(1),s_right(1)
          ind=i+ii;if(ind>N(1)+1) ind=2*(N(1)+1)-ind;if(ind<1) ind=2-ind
          tmp=tmp+w_right_1(ii)*f(ind,j)
        enddo
        buf3d(3,i,j)=tmp !fx(1,0)       
      enddo
    enddo
    do i=1,N(1)+1
      do j=1,N(2)
        tmp=0._f64
        do ii=r_left(2),s_left(2)
          ind=modulo(j+ii-1+N(2),N(2))+1
          tmp=tmp+w_left_2(ii)*f(i,ind)
        enddo
        buf3d(4,i,j)=tmp !fy(0,0)
        tmp=0._f64
        do ii=r_right(2),s_right(2)
          ind=modulo(j+ii-1+N(2),N(2))+1
          tmp=tmp+w_right_2(ii)*f(i,ind)
        enddo
        buf3d(5,i,j)=tmp !fy(0,1)               
      enddo
    enddo

    do j=1,N(2)
      do i=1,N(1)+1
        tmp=0._f64
        do ii=r_left(1),s_left(1)
          ind=i+ii;if(ind>N(1)+1) ind=2*(N(1)+1)-ind;if(ind<1) ind=2-ind
          tmp=tmp+w_left_1(ii)*buf3d(4,ind,j)
        enddo
        buf3d(6,i,j)=tmp !fxy(0,0)
        tmp=0._f64
        do ii=r_right(1),s_right(1)
          ind=i+ii;if(ind>N(1)+1) ind=2*(N(1)+1)-ind;if(ind<1) ind=2-ind
          tmp=tmp+w_right_1(ii)*buf3d(4,ind,j)
        enddo
        buf3d(7,i,j)=tmp !fxy(1,0)       
        tmp=0._f64
        do ii=r_left(1),s_left(1)
          ind=i+ii;if(ind>N(1)+1) ind=2*(N(1)+1)-ind;if(ind<1) ind=2-ind
          tmp=tmp+w_left_1(ii)*buf3d(5,ind,j)
        enddo
        buf3d(8,i,j)=tmp  !fxy(0,1)
        tmp=0._f64
        do ii=r_right(1),s_right(1)
          ind=i+ii;if(ind>N(1)+1) ind=2*(N(1)+1)-ind;if(ind<1) ind=2-ind
          tmp=tmp+w_right_1(ii)*buf3d(5,ind,j)
        enddo
        buf3d(9,i,j)=tmp !fxy(1,1)        
      enddo
    enddo

    buf3d(:,:,N(2)+1)=buf3d(:,:,1)
    
    
    !print *,'d=',d,maxval(abs(buf3d))
    
  end subroutine hermite_coef_nat_per


  function interpolate_value_hermite_2d( eta1, eta2, interp ) result(res)
    sll_real64,intent(in) :: eta1
    sll_real64,intent(in) :: eta2
    type(sll_hermite_interpolation_2d), pointer :: interp
    sll_real64 :: res
    sll_real64 :: eta_tmp(2)
    sll_int32 :: ii(2)
    
    eta_tmp(1:2) = (/eta1,eta2/)
    
    !warning: to be changed for dealing with general boundary conditions
    
    call localize_nat(ii(1),eta_tmp(1),interp%eta_min(1),interp%eta_max(1),interp%Nc(1))
    !print *,'#eta_min=',interp%eta_min
    !print *,'#eta_max=',interp%eta_max
    !print *,'#before',ii,eta_tmp(1:2),interp%Nc
    call localize_per(ii(2),eta_tmp(2),interp%eta_min(2),interp%eta_max(2),interp%Nc(2))
    !print *,'#before',ii,eta_tmp(1:2),interp%Nc
    call interpolate_hermite(interp%deriv,ii,eta_tmp,res,interp%Nc)
    !print *,'#after'
    
  end function interpolate_value_hermite_2d

 subroutine interpolate_hermite(f,i,x,fval,N)
    sll_int32,intent(in)::i(2),N(2)
    !real(f64),intent(in)::xx(2),xmin(2),xmax(2)
    sll_real64,intent(in)::x(2)
    sll_real64,intent(out)::fval
    sll_real64,dimension(0:8,0:N(1),0:N(2))::f
    !integer::i(2),i1(2),s
    sll_int32::i1(2),s
    sll_real64::w(2,0:3),tmp(0:3)
    sll_real64::g(0:3,0:3)
    
    !fval =f(0,i(1),i(2))!real(i(1),f64)
    
    !return
    
    do s=1,2
      w(s,0)=(2._f64*x(s)+1)*(1._f64-x(s))*(1._f64-x(s));
      w(s,1)=x(s)*x(s)*(3._f64-2._f64*x(s))
      w(s,2)=x(s)*(1._f64-x(s))*(1._f64-x(s))
      w(s,3)=x(s)*x(s)*(x(s)-1._f64)
      i1(s)=i(s)+1
    enddo

    
    g(0,0)=f(0,i(1),i(2))          !f(0,0)
    g(1,0)=f(0,i1(1),i(2))         !f(1,0)
    g(2,0)=f(1,i(1),i(2))          !fx(0,0)
    g(3,0)=f(2,i(1),i(2))          !fx(1,0)
    g(0,1)=f(0,i(1),i1(2))         !f(0,1) 
    g(1,1)=f(0,i1(1),i1(2))        !f(1,1)
    g(2,1)=f(1,i(1),i1(2))         !fx(0,1)
    g(3,1)=f(2,i(1),i1(2))         !fx(1,1)
    g(0,2)=f(3,i(1),i(2))          !fy(0,0)
    g(1,2)=f(3,i1(1),i(2))         !fy(1,0)
    g(2,2)=f(5,i(1),i(2))          !fxy(0,0)
    g(3,2)=f(6,i(1),i(2))          !fxy(1,0)
    g(0,3)=f(4,i(1),i(2))          !fy(0,1) 
    g(1,3)=f(4,i1(1),i(2))         !fy(1,1)
    g(2,3)=f(7,i(1),i(2))          !fxy(0,1) 
    g(3,3)=f(8,i(1),i(2))          !fxy(1,1)



    do s=0,3
      tmp(s)=w(1,0)*g(0,s)+w(1,1)*g(1,s)+w(1,2)*g(2,s)+w(1,3)*g(3,s)
    enddo  

    fval=w(2,0)*tmp(0)+w(2,1)*tmp(1)+w(2,2)*tmp(2)+w(2,3)*tmp(3)
  
    !print *,fval,' t',f    
  end subroutine interpolate_hermite

  subroutine localize_per(i,x,xmin,xmax,N)
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
  end subroutine localize_per
  

  subroutine localize_nat(i,x,xmin,xmax,N)
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
  end subroutine localize_nat


  
    

end module
