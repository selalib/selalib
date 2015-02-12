module classical_conservative_semi_lagrangian
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"

  use sll_constants
  use cubic_non_uniform_splines
  !use utils
  implicit none

contains 
  
  subroutine csl_advection_per(f,spl_per,Xstar,node_positions,N)
    !Xstar and node_positions are normalized to [0,1]
    use sll_constants
    use cubic_non_uniform_splines
    implicit none
    
    sll_real64,dimension(:),pointer::f,Xstar,node_positions
    type(cubic_nonunif_spline_1D), pointer :: spl_per
    sll_int32,intent(in):: N
    sll_real64 :: dx
    sll_int32  :: i
    sll_real64 :: M,tmp,tmp2
    dx = 1._f64/real(N,f64)
    
    
    !do i=1,N+1
    !  print *,i,node_positions(i)
    !enddo

    !do i=1,N+1
    !  print *,i,Xstar(i),f(i)
    !enddo
    
    !print *,dx
    
    do i=1,N+1
      do while (Xstar(i).gt.1._f64)
        Xstar(i) = Xstar(i)-1._f64
      end do
      do while (Xstar(i).lt.0._f64)
        Xstar(i) = Xstar(i)+1._f64
      end do    
    enddo



    !from f compute the mean
    do i=0,N-1
      f(i+1)=f(i+1)*(node_positions(i+2)-node_positions(i+1))/dx
    enddo
    
    
    !we compute the splines coefficients by solving the LU decomposition
    M=0._f64
    do i=1,N
      M=M+f(i)
    enddo
    !M=M/real(N,rk)
    do i=1,N
      f(i)=f(i)-M*(node_positions(i+1)-node_positions(i))!/dx
    enddo    
    !f_per(1)=0._f64
    !do i=2,N
    !  f_per(i)=f_per(i-1)+f(i-1)
    !enddo
    !f=f_per
    tmp=f(1)
    f(1)=0._f64
    do i=2,N
      tmp2=f(i)
      f(i)=f(i-1)+tmp
      tmp=tmp2
    enddo
    
    
    !call of compute_spline and interpolations

    call compute_spline_nonunif( f, spl_per,node_positions)
    !print *,spl_per%xmin,spl_per%xmax,node_positions(1),node_positions(N+1)
    
    
    !print *,spl_per%buf(2),spl_per%buf(3),spl_per%buf(1)
    !print *,spl_per%buf(4),spl_per%buf(5),spl_per%buf(6)
    !print *,spl_per%buf(9),spl_per%buf(7),spl_per%buf(8)
    !stop
    
    call interpolate_array_value_nonunif( Xstar, f, N, spl_per)
    
    
    tmp=f(1)!;for(i=0;i<Nx-1;i++)p[i]=p[i+1]-p[i];p[Nx-1]=tmp+M-p[Nx-1];
    do i=1,N-1
      f(i)=f(i+1)-f(i)+M*(node_positions(i+1)-node_positions(i))
    enddo
    f(N)=tmp-f(N)+M*(node_positions(1)+1._f64-node_positions(N))


    !from mean compute f
    do i=1,N
      f(i)=f(i)*dx/(node_positions(i+1)-node_positions(i))
    enddo

    f(N+1) = f(1)
    
    
    
  end subroutine csl_advection_per

subroutine poisson1dpertrap(E,L,N)
  use sll_constants
  implicit none
  sll_int,intent(in)::N
  sll_real64,dimension(N+1),intent(inout)::E
  sll_real64,intent(in)::L
  sll_int::i
  sll_real64::eold,enew,dx2,tmp
  !ensures at first that Ein is of mean zero
  
  tmp=0._f64
  do i=1,N
    tmp=tmp+E(i)
  enddo
  tmp=-tmp/real(N,f64)
  do i=1,N
    E(i)=E(i)+tmp
  enddo
  
  dx2=0.5_f64*L/real(N,f64)
  eold=E(1)
  E(1)=0._f64
  tmp=0._f64
  do i=1,N-1
    enew=E(i+1)
    E(i+1)=E(i)+dx2*(eold+enew)
    tmp=tmp+E(i+1)
    eold=enew
  enddo
  tmp=-tmp/real(N,f64)
  do i=1,N
    E(i)=E(i)+tmp
  enddo
  E(N+1)=E(1)
end subroutine poisson1dpertrap


subroutine compute_psi(a1,a2,rho,nc_eta1,nc_eta2,psi,phi_poisson,E,&
geom_x,x1n_array,x2n_array,x1c_array,x2c_array,jac_array,delta_eta1,delta_eta2,div_case)
  use sll_constants
  implicit none

  sll_int,intent(in) :: nc_eta1,nc_eta2
  sll_real64,dimension(1:nc_eta1+1) :: rho
  sll_real64,dimension(1:nc_eta1+1) :: phi_poisson
  sll_real64,dimension(1:nc_eta1+1) :: E
  sll_real64,dimension(1:nc_eta1+1,1:nc_eta2+1) :: x1n_array
  sll_real64,dimension(1:nc_eta1+1,1:nc_eta2+1) :: x2n_array
  sll_real64,dimension(1:nc_eta1+1,1:nc_eta2+1) :: x1c_array
  sll_real64,dimension(1:nc_eta1+1,1:nc_eta2+1) :: x2c_array
  sll_real64,dimension(1:nc_eta1+1,1:nc_eta2+1) :: jac_array
  sll_real64,dimension(1:nc_eta1+1,1:nc_eta2+1) :: a1
  sll_real64,dimension(1:nc_eta1+1,1:nc_eta2+1) :: a2
  sll_real64,dimension(1:nc_eta1+1,1:nc_eta2+1) :: psi
  sll_real64 :: x1_min,x1_max,x2_min,x2_max
  sll_int :: i1,i2,ii,i2p1,i2m1,i1p1,i1m1,ii1(-10:10),ii2(-10:10)
  sll_real64 :: tmp,x1,x2,phi_val,xx,a(-10:10)
  sll_real64,intent(in) :: delta_eta1,delta_eta2
  sll_real64,intent(in) :: geom_x(2,2)
  sll_int,intent(in) :: div_case
  
  !a[-2] = 0, a[-1] = -1/6, a[0] = 1, a[1] = -1/2, a[2] = -1/3
  x1_min = geom_x(1,1)
  x1_max = geom_x(2,1)
  x2_min = geom_x(1,2)
  x2_max = geom_x(2,2)



    E=rho-1._f64
    call poisson1dpertrap(E,x1_max-x1_min,nc_eta1)
    phi_poisson = E
    call poisson1dpertrap(phi_poisson,x1_max-x1_min,nc_eta1)
    tmp = phi_poisson(1)
    do i1=1,nc_eta1
      phi_poisson(i1) = -phi_poisson(i1) + tmp
    enddo
    phi_poisson(nc_eta1+1) = phi_poisson(1)
    
    do i1=1,nc_eta1!+1
      do i2=1,nc_eta2!+1
        !x1 = x1n_array(i1,i2)
        !x2 = x2n_array(i1,i2)
        x1 = x1c_array(i1,i2)
        x2 = x2c_array(i1,i2)
        phi_val = 0._f64
        xx = (x1-x1_min)/(x1_max-x1_min)-0.5_f64/real(nc_eta1,f64)
        if(xx<=0._f64)then
          xx = 0._f64
        endif
        if(xx>=1._f64)then
          xx = xx-1._f64!1._f64-1e-15_f64
        endif
        if(xx<=0._f64)then
          xx = 0._f64
        endif      
        xx = xx*real(nc_eta1,f64)
        ii = floor(xx)
        xx = xx-real(ii,f64)      
        phi_val = (1._f64-xx)*phi_poisson(ii+1)+xx*phi_poisson(ii+2)     
       psi( i1, i2 ) = ( 0.5_f64*x2**2+phi_val)!& utilisation de tableau abusive 
      enddo  
    enddo
   
     if(div_case==0)then
       !we suppose for the moment cell centered psi an derivatives
       !order 2 centered
       do i1=1,nc_eta1
         do i2=1,nc_eta2
           i1p1=modulo(i1+1-1+nc_eta1,nc_eta1)+1
           i1m1=modulo(i1-1-1+nc_eta1,nc_eta1)+1
           i2p1=modulo(i2+1-1+nc_eta2,nc_eta2)+1
           i2m1=modulo(i2-1-1+nc_eta2,nc_eta2)+1
           !a1(i1,i2)=((psi(i1,i2+1)-psi(i1,modulo(i2-1-1+nc_eta2,nc_eta2)+1))/(2._f64*delta_eta2))*(x1_max-x1_min)/jac_array(i1,i2)
           !a2(i1,i2)=-((psi(i1+1,i2)-psi(modulo(i1-1-1+nc_eta1,nc_eta1)+1,i2))/(2._f64*delta_eta1))*(x2_max-x2_min)/jac_array(i1,i2)
           !a1(i1,i2)=((psi(i1,i2p1)-psi(i1,i2m1))/(2._f64*delta_eta2))*(x1_max-x1_min)/jac_array(i1,i2)
           !a2(i1,i2)=-((psi(i1p1,i2)-psi(i1m1,i2))/(2._f64*delta_eta1))*(x2_max-x2_min)/jac_array(i1,i2)
           a1(i1,i2)=((psi(i1,i2p1)-psi(i1,i2m1))/(2._f64*delta_eta2))/jac_array(i1,i2)
           a2(i1,i2)=-((psi(i1p1,i2)-psi(i1m1,i2))/(2._f64*delta_eta1))/jac_array(i1,i2)


         enddo
       enddo
    
        do i1=1,nc_eta1+1
           a1(i1,nc_eta2+1)=a1(i1,1)
           a2(i1,nc_eta2+1)=a2(i1,1)
        enddo
        do i2=1,nc_eta2+1
           a1(nc_eta1+1,i2)=a1(1,i2)
           a2(nc_eta1+1,i2)=a2(1,i2)    
        enddo
     endif
     
     if(div_case==1)then
       !order 4 centered
       a(-2) = 1._f64/12._f64
       a(-1) =-2._f64/3._f64
       a(0) =0._f64
       a(1) =2._f64/3._f64
       a(2)=-1._f64/12._f64

       !a(-2) = 0._f64
       !a(-1) =-1._f64/2._f64
       !a(0) =0._f64
       !a(1) =1._f64/2._f64
       !a(2)=0._f64

 
     
       do i1=1,nc_eta1
         do i2=1,Nc_eta2
            do ii=-2,2
              ii1(ii) = modulo(i1+ii-1+nc_eta1,nc_eta1)+1
              ii2(ii) = modulo(i2+ii-1+nc_eta2,nc_eta2)+1
            enddo
            a1(i1,i2)=0._f64
            do ii=-2,2
              a1(i1,i2)=a1(i1,i2)+a(ii)*psi(i1,ii2(ii))
            enddo
            a2(i1,i2)=0._f64
            do ii=-2,2
              a2(i1,i2)=a2(i1,i2)-a(ii)*psi(ii1(ii),i2)
            enddo
         enddo
       enddo
        do i1=1,nc_eta1
          do i2=1,Nc_eta2
            a1(i1,i2)=(a1(i1,i2)/(delta_eta2))/jac_array(i1,i2)
            a2(i1,i2)=(a2(i1,i2)/(delta_eta1))/jac_array(i1,i2)
          enddo
        enddo
         
        do i1=1,nc_eta1+1
           a1(i1,nc_eta2+1)=a1(i1,1)
           a2(i1,nc_eta2+1)=a2(i1,1)
        enddo
        do i2=1,nc_eta2+1
           a1(nc_eta1+1,i2)=a1(1,i2)
           a2(nc_eta1+1,i2)=a2(1,i2)    
        enddo
     endif
     if(div_case==2)then
       a(-2)= 1._f64/6._f64
       a(-1) =-1._f64
       a(0) =1._f64/2._f64
       a(1) =1._f64/3._f64
       a(2) = 0._f64
     
       do i1=1,nc_eta1
         do i2=1,Nc_eta2
             a1(i1,i2)=a(-2)*psi(i1, modulo(i2-2 -1+nc_eta2,nc_eta2)+1) + &
             a(2)*psi(i1, modulo(i2+2 -1+nc_eta2,nc_eta2)+1) + a(1)*psi(i1,i2+1)+&
             a(0)*psi(i1, modulo(i2 -1+nc_eta2,nc_eta2)+1) + &
             a(-1)*psi(i1, modulo(i2-1 -1+nc_eta2,nc_eta2)+1)
             a2(i1,i2)=-(a(-2)*psi(modulo(i1-2 -1+nc_eta1,nc_eta1)+1, i2 ) +&
             a(2)*psi(modulo(i1+2-1+nc_eta1,nc_eta1)+1, i2 ) + a(1)*psi(i1+1,i2)+&
             a(0)*psi(modulo(i1 -1+nc_eta1,nc_eta1)+1, i2 ) + &
             a(-1)*psi(modulo(i1-1 -1+nc_eta1,nc_eta1)+1, i2 ))
             
         enddo
       enddo
        do i1=1,nc_eta1
          do i2=1,Nc_eta2
            a1(i1,i2)=(a1(i1,i2)/(delta_eta2))*(x1_max-x1_min)/jac_array(i1,i2)
            a2(i1,i2)=(a2(i1,i2)/(delta_eta1))*(x2_max-x2_min)/jac_array(i1,i2)
          enddo
        enddo
         
        do i1=1,nc_eta1+1
           a1(i1,nc_eta2+1)=a1(i1,1)
           a2(i1,nc_eta2+1)=a2(i1,1)
        enddo
        do i2=1,nc_eta2+1
           a1(nc_eta1+1,i2)=a1(1,i2)
           a2(nc_eta1+1,i2)=a2(1,i2)    
        enddo
     endif


     if(div_case==3)then
       do i1=1,nc_eta1
         do i2=1,Nc_eta2
           a1(i1,i2)=((psi(i1,i2+1)-psi(i1,modulo(i2-1-1+nc_eta2,nc_eta2)+1))/(2._f64*delta_eta2))*(x1_max-x1_min)/jac_array(i1,i2)
           a2(i1,i2)=-((psi(i1+1,i2)-psi(modulo(i1-1-1+nc_eta1,nc_eta1)+1,i2))/(2._f64*delta_eta1))*(x2_max-x2_min)/jac_array(i1,i2)
         enddo
       enddo
    
        do i1=1,nc_eta1+1
           a1(i1,nc_eta2+1)=a1(i1,1)
           a2(i1,nc_eta2+1)=a2(i1,1)
        enddo
        do i2=1,nc_eta2+1
           a1(nc_eta1+1,i2)=a1(1,i2)
           a2(nc_eta1+1,i2)=a2(1,i2)    
        enddo



       a(-2)= 1._f64/6._f64
       a(-1) =-1._f64
       a(0) =1._f64/2._f64
       a(1) =1._f64/3._f64
       a(2) = 0._f64
     
       do i1=1,nc_eta1
         do i2=1,Nc_eta2
            if(1+0*a1(i1,i2)>0._f64)then
            a1(i1,i2)=a(-2)*psi(i1, modulo(i2-2 -1+nc_eta2,nc_eta2)+1) + &
            a(2)*psi(i1, modulo(i2+2 -1+nc_eta2,nc_eta2)+1) + a(1)*psi(i1,i2+1)+&
             a(0)*psi(i1, modulo(i2 -1+nc_eta2,nc_eta2)+1) + &
              a(-1)*psi(i1, modulo(i2-1 -1+nc_eta2,nc_eta2)+1)
            a1(i1,i2)=(a1(i1,i2)/(delta_eta2))*(x1_max-x1_min)/jac_array(i1,i2)
            endif
            if(1+0*a2(i1,i2)>0._f64)then  
            a2(i1,i2)=-(a(-2)*psi(modulo(i1-2 -1+nc_eta1,nc_eta1)+1, i2 ) +&
             a(2)*psi(modulo(i1+2-1+nc_eta1,nc_eta1)+1, i2 ) + a(1)*psi(i1+1,i2)+&
             a(0)*psi(modulo(i1 -1+nc_eta1,nc_eta1)+1, i2 ) + &
             a(-1)*psi(modulo(i1-1 -1+nc_eta1,nc_eta1)+1, i2 ))
            a2(i1,i2)=(a2(i1,i2)/(delta_eta1))*(x2_max-x2_min)/jac_array(i1,i2)
            endif 
         enddo
       enddo


       a(-1) =-1._f64/3._f64
       a(0) =-1._f64/2._f64
       a(1) =1._f64
       a(2)=-1._f64/6._f64
       a(3)=0._f64

     
       do i1=1,nc_eta1
         do i2=1,Nc_eta2
          if(a1(i1,i2)<0._f64)then
            a1(i1,i2)=a(3)*psi(i1, modulo(i2+3 -1+nc_eta2,nc_eta2)+1) + &
            a(2)*psi(i1, modulo(i2+2 -1+nc_eta2,nc_eta2)+1) + a(1)*psi(i1,i2+1)+&
             a(0)*psi(i1, modulo(i2 -1+nc_eta2,nc_eta2)+1) + &
              a(-1)*psi(i1, modulo(i2-1 -1+nc_eta2,nc_eta2)+1)
             a1(i1,i2)=(a1(i1,i2)/(delta_eta2))*(x1_max-x1_min)/jac_array(i1,i2)
           endif     
          if(a2(i1,i2)<0._f64)then
            a2(i1,i2)=-(a(3)*psi(modulo(i1+3 -1+nc_eta1,nc_eta1)+1, i2 ) +&
             a(2)*psi(modulo(i1+2-1+nc_eta1,nc_eta1)+1, i2 ) + a(1)*psi(i1+1,i2)+&
             a(0)*psi(modulo(i1 -1+nc_eta1,nc_eta1)+1, i2 ) + &
             a(-1)*psi(modulo(i1-1 -1+nc_eta1,nc_eta1)+1, i2 ))
             a2(i1,i2)=(a2(i1,i2)/(delta_eta1))*(x2_max-x2_min)/jac_array(i1,i2)
          endif   
         enddo
       enddo



         
        do i1=1,nc_eta1+1
           a1(i1,nc_eta2+1)=a1(i1,1)
           a2(i1,nc_eta2+1)=a2(i1,1)
        enddo
        do i2=1,nc_eta2+1
           a1(nc_eta1+1,i2)=a1(1,i2)
           a2(nc_eta1+1,i2)=a2(1,i2)    
        enddo
     endif



     if(div_case==30)then



       a(-2)= 1._f64/6._f64
       a(-1) =-1._f64
       a(0) =1._f64/2._f64
       a(1) =1._f64/3._f64
       a(2) = 0._f64
     
       do i1=1,nc_eta1
         do i2=1,Nc_eta2
            a1(i1,i2)=a(-2)*psi(i1, modulo(i2-2 -1+nc_eta2,nc_eta2)+1) + &
            a(2)*psi(i1, modulo(i2+2 -1+nc_eta2,nc_eta2)+1) + a(1)*psi(i1,i2+1)+&
             a(0)*psi(i1, modulo(i2 -1+nc_eta2,nc_eta2)+1) + &
              a(-1)*psi(i1, modulo(i2-1 -1+nc_eta2,nc_eta2)+1)
            a1(i1,i2)=(a1(i1,i2)/(delta_eta2))*(x1_max-x1_min)/jac_array(i1,i2)
            a2(i1,i2)=-(a(-2)*psi(modulo(i1-2 -1+nc_eta1,nc_eta1)+1, i2 ) +&
             a(2)*psi(modulo(i1+2-1+nc_eta1,nc_eta1)+1, i2 ) + a(1)*psi(i1+1,i2)+&
             a(0)*psi(modulo(i1 -1+nc_eta1,nc_eta1)+1, i2 ) + &
             a(-1)*psi(modulo(i1-1 -1+nc_eta1,nc_eta1)+1, i2 ))
            a2(i1,i2)=(a2(i1,i2)/(delta_eta1))*(x2_max-x2_min)/jac_array(i1,i2)
          enddo
       enddo


       a(-1) =-1._f64/3._f64
       a(0) =-1._f64/2._f64
       a(1) =1._f64
       a(2)=-1._f64/6._f64
       a(3)=0._f64

     
       do i1=1,nc_eta1
         do i2=1,Nc_eta2
            tmp=a(3)*psi(i1, modulo(i2+3 -1+nc_eta2,nc_eta2)+1) + &
            a(2)*psi(i1, modulo(i2+2 -1+nc_eta2,nc_eta2)+1) + a(1)*psi(i1,i2+1)+&
             a(0)*psi(i1, modulo(i2 -1+nc_eta2,nc_eta2)+1) + &
              a(-1)*psi(i1, modulo(i2-1 -1+nc_eta2,nc_eta2)+1)
             tmp=(tmp/(delta_eta2))*(x1_max-x1_min)/jac_array(i1,i2)
            a1(i1,i2)=0.5_f64*(a1(i1,i2)+tmp)
            tmp=-(a(3)*psi(modulo(i1+3 -1+nc_eta1,nc_eta1)+1, i2 ) +&
             a(2)*psi(modulo(i1+2-1+nc_eta1,nc_eta1)+1, i2 ) + a(1)*psi(i1+1,i2)+&
             a(0)*psi(modulo(i1 -1+nc_eta1,nc_eta1)+1, i2 ) + &
             a(-1)*psi(modulo(i1-1 -1+nc_eta1,nc_eta1)+1, i2 ))
             tmp=(tmp/(delta_eta1))*(x2_max-x2_min)/jac_array(i1,i2)
             a2(i1,i2)=0.5_f64*(a2(i1,i2)+tmp)
         enddo
       enddo



         
        do i1=1,nc_eta1+1
           a1(i1,nc_eta2+1)=a1(i1,1)
           a2(i1,nc_eta2+1)=a2(i1,1)
        enddo
        do i2=1,nc_eta2+1
           a1(nc_eta1+1,i2)=a1(1,i2)
           a2(nc_eta1+1,i2)=a2(1,i2)    
        enddo

    
     endif



    if(div_case==4)then
       a(-3) = -1._f64/30._f64
       a(-2) = 1._f64/4._f64
       a(-1) = -1._f64
       a(0)  = 1._f64/3._f64
       a(1)  = 1._f64/2._f64
       a(2)  = -1._f64/20._f64
       a(3)  = 0._f64
       !a(-1) =-1._f64/3._f64
       !a(0) =-1._f64/2._f64
       !a(1) =1._f64
       !a(2)=-1._f64/6._f64
       !a(3)=0._f64

     
       do i1=1,nc_eta1
         do i2=1,Nc_eta2
            a1(i1,i2)=a(3)*psi(i1, modulo(i2+3 -1+nc_eta2,nc_eta2)+1) + &
            a(2)*psi(i1, modulo(i2+2 -1+nc_eta2,nc_eta2)+1) + a(1)*psi(i1,i2+1)+&
             a(0)*psi(i1, modulo(i2 -1+nc_eta2,nc_eta2)+1) + &
              a(-1)*psi(i1, modulo(i2-1 -1+nc_eta2,nc_eta2)+1) + &
              a(-2)*psi(i1, modulo(i2-2 -1+nc_eta2,nc_eta2)+1)+ &
              a(-3)*psi(i1, modulo(i2-3 -1+nc_eta2,nc_eta2)+1)

            a2(i1,i2)=-(a(3)*psi(modulo(i1+3 -1+nc_eta1,nc_eta1)+1, i2 ) +&
             a(2)*psi(modulo(i1+2-1+nc_eta1,nc_eta1)+1, i2 ) + a(1)*psi(i1+1,i2)+&
             a(0)*psi(modulo(i1 -1+nc_eta1,nc_eta1)+1, i2 ) + &
             a(-1)*psi(modulo(i1-1 -1+nc_eta1,nc_eta1)+1, i2 )+ &
             a(-2)*psi(modulo(i1-2 -1+nc_eta1,nc_eta1)+1, i2 )+ &
             a(-3)*psi(modulo(i1-3 -1+nc_eta1,nc_eta1)+1, i2 ))
             
         enddo
       enddo
        do i1=1,nc_eta1
          do i2=1,Nc_eta2
            a1(i1,i2)=(a1(i1,i2)/(delta_eta2))*(x1_max-x1_min)/jac_array(i1,i2)
            a2(i1,i2)=(a2(i1,i2)/(delta_eta1))*(x2_max-x2_min)/jac_array(i1,i2)
          enddo
        enddo
         
        do i1=1,nc_eta1+1
           a1(i1,nc_eta2+1)=a1(i1,1)
           a2(i1,nc_eta2+1)=a2(i1,1)
        enddo
        do i2=1,nc_eta2+1
           a1(nc_eta1+1,i2)=a1(1,i2)
           a2(nc_eta1+1,i2)=a2(1,i2)    
        enddo
     endif
      if(div_case==5)then
       a(-3) = 0._f64
       a(-2) = 1._f64/20._f64
       a(-1) = -1._f64/2._f64
       a(0)  = -1._f64/3._f64
       a(1)  = 1._f64
       a(2)  = -1._f64/4._f64
       a(3)  = 1._f64/30._f64
       !a(-1) =-1._f64/3._f64
       !a(0) =-1._f64/2._f64
       !a(1) =1._f64
       !a(2)=-1._f64/6._f64
       !a(3)=0._f64

     
       do i1=1,nc_eta1
         do i2=1,Nc_eta2
            a1(i1,i2)=a(3)*psi(i1, modulo(i2+3 -1+nc_eta2,nc_eta2)+1) + &
            a(2)*psi(i1, modulo(i2+2 -1+nc_eta2,nc_eta2)+1) + a(1)*psi(i1,i2+1)+&
             a(0)*psi(i1, modulo(i2 -1+nc_eta2,nc_eta2)+1) + &
              a(-1)*psi(i1, modulo(i2-1 -1+nc_eta2,nc_eta2)+1) + &
              a(-2)*psi(i1, modulo(i2-2 -1+nc_eta2,nc_eta2)+1)+ &
              a(-3)*psi(i1, modulo(i2-3 -1+nc_eta2,nc_eta2)+1)

            a2(i1,i2)=-(a(3)*psi(modulo(i1+3 -1+nc_eta1,nc_eta1)+1, i2 ) +&
             a(2)*psi(modulo(i1+2-1+nc_eta1,nc_eta1)+1, i2 ) + a(1)*psi(i1+1,i2)+&
             a(0)*psi(modulo(i1 -1+nc_eta1,nc_eta1)+1, i2 ) + &
             a(-1)*psi(modulo(i1-1 -1+nc_eta1,nc_eta1)+1, i2 )+ &
             a(-2)*psi(modulo(i1-2 -1+nc_eta1,nc_eta1)+1, i2 )+ &
             a(-3)*psi(modulo(i1-3 -1+nc_eta1,nc_eta1)+1, i2 ))
             
         enddo
       enddo
        do i1=1,nc_eta1
          do i2=1,Nc_eta2
            a1(i1,i2)=(a1(i1,i2)/(delta_eta2))*(x1_max-x1_min)/jac_array(i1,i2)
            a2(i1,i2)=(a2(i1,i2)/(delta_eta1))*(x2_max-x2_min)/jac_array(i1,i2)
          enddo
        enddo
         
        do i1=1,nc_eta1+1
           a1(i1,nc_eta2+1)=a1(i1,1)
           a2(i1,nc_eta2+1)=a2(i1,1)
        enddo
        do i2=1,nc_eta2+1
           a1(nc_eta1+1,i2)=a1(1,i2)
           a2(nc_eta1+1,i2)=a2(1,i2)    
        enddo
     endif



      if(div_case==50)then
       a(-3) = 0._f64
       a(-2) = 1._f64/20._f64
       a(-1) = -1._f64/2._f64
       a(0)  = -1._f64/3._f64
       a(1)  = 1._f64
       a(2)  = -1._f64/4._f64
       a(3)  = 1._f64/30._f64
       !a(-1) =-1._f64/3._f64
       !a(0) =-1._f64/2._f64
       !a(1) =1._f64
       !a(2)=-1._f64/6._f64
       !a(3)=0._f64

     
       do i1=1,nc_eta1
         do i2=1,Nc_eta2
            a1(i1,i2)=a(3)*psi(i1, modulo(i2+3 -1+nc_eta2,nc_eta2)+1) + &
            a(2)*psi(i1, modulo(i2+2 -1+nc_eta2,nc_eta2)+1) + a(1)*psi(i1,i2+1)+&
             a(0)*psi(i1, modulo(i2 -1+nc_eta2,nc_eta2)+1) + &
              a(-1)*psi(i1, modulo(i2-1 -1+nc_eta2,nc_eta2)+1) + &
              a(-2)*psi(i1, modulo(i2-2 -1+nc_eta2,nc_eta2)+1)+ &
              a(-3)*psi(i1, modulo(i2-3 -1+nc_eta2,nc_eta2)+1)

            a2(i1,i2)=-(a(3)*psi(modulo(i1+3 -1+nc_eta1,nc_eta1)+1, i2 ) +&
             a(2)*psi(modulo(i1+2-1+nc_eta1,nc_eta1)+1, i2 ) + a(1)*psi(i1+1,i2)+&
             a(0)*psi(modulo(i1 -1+nc_eta1,nc_eta1)+1, i2 ) + &
             a(-1)*psi(modulo(i1-1 -1+nc_eta1,nc_eta1)+1, i2 )+ &
             a(-2)*psi(modulo(i1-2 -1+nc_eta1,nc_eta1)+1, i2 )+ &
             a(-3)*psi(modulo(i1-3 -1+nc_eta1,nc_eta1)+1, i2 ))
             
         enddo
       enddo

       a(-3) = -1._f64/30._f64
       a(-2) = 1._f64/4._f64
       a(-1) = -1._f64
       a(0)  = 1._f64/3._f64
       a(1)  = 1._f64/2._f64
       a(2)  = -1._f64/20._f64
       a(3)  = 0._f64
       !a(-1) =-1._f64/3._f64
       !a(0) =-1._f64/2._f64
       !a(1) =1._f64
       !a(2)=-1._f64/6._f64
       !a(3)=0._f64

     
       do i1=1,nc_eta1
         do i2=1,Nc_eta2
            tmp=a(3)*psi(i1, modulo(i2+3 -1+nc_eta2,nc_eta2)+1) + &
            a(2)*psi(i1, modulo(i2+2 -1+nc_eta2,nc_eta2)+1) + a(1)*psi(i1,i2+1)+&
             a(0)*psi(i1, modulo(i2 -1+nc_eta2,nc_eta2)+1) + &
              a(-1)*psi(i1, modulo(i2-1 -1+nc_eta2,nc_eta2)+1) + &
              a(-2)*psi(i1, modulo(i2-2 -1+nc_eta2,nc_eta2)+1)+ &
              a(-3)*psi(i1, modulo(i2-3 -1+nc_eta2,nc_eta2)+1)
            a1(i1,i2)=0.5_f64*(a1(i1,i2)+tmp)
            
            tmp=-(a(3)*psi(modulo(i1+3 -1+nc_eta1,nc_eta1)+1, i2 ) +&
             a(2)*psi(modulo(i1+2-1+nc_eta1,nc_eta1)+1, i2 ) + a(1)*psi(i1+1,i2)+&
             a(0)*psi(modulo(i1 -1+nc_eta1,nc_eta1)+1, i2 ) + &
             a(-1)*psi(modulo(i1-1 -1+nc_eta1,nc_eta1)+1, i2 )+ &
             a(-2)*psi(modulo(i1-2 -1+nc_eta1,nc_eta1)+1, i2 )+ &
             a(-3)*psi(modulo(i1-3 -1+nc_eta1,nc_eta1)+1, i2 ))
             a2(i1,i2)=0.5_f64*(a2(i1,i2)+tmp)
         enddo
       enddo




        do i1=1,nc_eta1
          do i2=1,Nc_eta2
            a1(i1,i2)=(a1(i1,i2)/(delta_eta2))*(x1_max-x1_min)/jac_array(i1,i2)
            a2(i1,i2)=(a2(i1,i2)/(delta_eta1))*(x2_max-x2_min)/jac_array(i1,i2)
          enddo
        enddo
         
        do i1=1,nc_eta1+1
           a1(i1,nc_eta2+1)=a1(i1,1)
           a2(i1,nc_eta2+1)=a2(i1,1)
        enddo
        do i2=1,nc_eta2+1
           a1(nc_eta1+1,i2)=a1(1,i2)
           a2(nc_eta1+1,i2)=a2(1,i2)    
        enddo
     endif

            
         



  
end subroutine compute_psi

subroutine advect_classical_csl_1(dt,a1,f,geom_x,N_x1,N_x2,buf,node_positions_x1,&
Xstar,spl_per_x1)
!solve \partial_t f(t,x,y) + \partial_x(a1(x,y)f(t,x,y))=0 over dt
  sll_int32,intent(in) :: N_x1,N_x2
  sll_real64,intent(in) :: dt,geom_x(2,2)
  sll_real64,dimension(:),pointer :: buf,node_positions_x1,Xstar
  sll_real64,dimension(:,:),pointer:: a1,f
  type(cubic_nonunif_spline_1D), pointer :: spl_per_x1
  sll_int32 :: i1,i2,i1m1
  do i2=1,N_x2
    buf(1:N_x1) = f(1:N_x1,i2)
    do i1=1,N_x1
      i1m1=modulo(i1-1-1+N_x1,N_x1)+1
      Xstar(i1) = node_positions_x1(i1)-0.5_f64*dt*(a1(i1,i2)+a1(i1m1,i2))
    enddo    
    call csl_advection_per(buf,spl_per_x1,Xstar,node_positions_x1,N_x1)
    f(1:N_x1+1,i2) = buf(1:N_x1+1)
  enddo
  

end subroutine advect_classical_csl_1


subroutine advect_classical_csl_2(dt,a2,f,geom_x,N_x1,N_x2,buf,node_positions_x2,&
Xstar,spl_per_x2)
!solve \partial_t f(t,x,y) + \partial_x(a1(x,y)f(t,x,y))=0 over dt
  sll_int32,intent(in) :: N_x1,N_x2
  sll_real64,intent(in) :: dt,geom_x(2,2)
  sll_real64,dimension(:),pointer :: buf,node_positions_x2,Xstar
  sll_real64,dimension(:,:),pointer:: a2,f
  type(cubic_nonunif_spline_1D), pointer :: spl_per_x2
  sll_int32 :: i1,i2,i2m1
  do i1=1,N_x1
    buf(1:N_x2) = f(i1,1:N_x2)
    Xstar(1:N_x2) = node_positions_x2(1:N_x2)!-dt*a2(i1,1:N_x2)
    do i2=1,N_x2
      i2m1=modulo(i2-1-1+N_x2,N_x2)+1
      Xstar(i2) = node_positions_x2(i2)-0.5_f64*dt*(a2(i1,i2)+a2(i1,i2m1))
    enddo    

    call csl_advection_per(buf,spl_per_x2,Xstar,node_positions_x2,N_x2)
    f(i1,1:N_x2+1) = buf(1:N_x2+1)
  enddo
  

end subroutine advect_classical_csl_2

subroutine advect_classical_csl(dt,a1,a2,f,geom_x,N_x1,N_x2,buf,&
node_positions_x1,node_positions_x2,Xstar,spl_per_x1,spl_per_x2)
!solve \partial_t f(t,x,y) + \partial_x(a1(x,y)f(t,x,y))=0 over dt
  sll_int32,intent(in) :: N_x1,N_x2
  sll_real64,intent(in) :: dt,geom_x(2,2)
  sll_real64,dimension(:),pointer :: buf,node_positions_x1,node_positions_x2,Xstar
  sll_real64,dimension(:,:),pointer:: a1,a2,f
  type(cubic_nonunif_spline_1D), pointer :: spl_per_x1,spl_per_x2
  
  call advect_classical_csl_1(0.5_f64*dt,a1,f,geom_x,N_x1,N_x2,buf,node_positions_x1,&
    Xstar,spl_per_x1)
  call advect_classical_csl_2(dt,a2,f,geom_x,N_x1,N_x2,buf,node_positions_x2,&
    Xstar,spl_per_x2)
  call advect_classical_csl_1(0.5_f64*dt,a1,f,geom_x,N_x1,N_x2,buf,node_positions_x1,&
    Xstar,spl_per_x1)

end subroutine advect_classical_csl



end module classical_conservative_semi_lagrangian