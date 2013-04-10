module cg_csl_uniform_module
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"

  use sll_constants
  use cubic_non_uniform_splines
  !use utils
  implicit none

contains 
  
  subroutine csl_advection_per(f,spl_per,Xstar,node_positions,N,interp_case)
    !Xstar and node_positions are normalized to [0,1]
    use sll_constants
    use cubic_nonuniform_splines
    implicit none
    
    sll_real64,dimension(:),pointer::f,Xstar,node_positions
    type(cubic_nonunif_spline_1D), pointer :: spl_per
    sll_int32,intent(in):: N,interp_case
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
      if(abs(Xstar(i))>2._f64)then
        print *,'displacement too big in csl_advection_per',Xstar(i)
        stop
      endif
      do while (Xstar(i).ge.1._f64)
        Xstar(i) = Xstar(i)-1._f64
      end do
      do while (Xstar(i).lt.0._f64)
        Xstar(i) = Xstar(i)+1._f64
      end do
      if(Xstar(i)>=1._f64)then
        Xstar(i) = Xstar(i)-1._f64
      endif
      if((Xstar(i)>=1).or.(Xstar(i)<0))then
        print *,'#problem for Xstar in csl_advection_per',i,Xstar(i)
        stop
      endif    
    enddo

    !tmp=0._f64
    !do i=4,N-4
    !  tmp2=abs(node_positions(i)-Xstar(i)-(node_positions(i-1)-Xstar(i-1)))
    !  if(tmp2>tmp)then
    !    tmp = tmp2
    !  endif
    !enddo
    !if(tmp>0.01/real(N,f64))then
    !  do i=1,N
    !    print *,real(i-1,f64)/real(N,f64),node_positions(i),f(i),Xstar(i)
    !  enddo
    !  stop
    !endif



    if(interp_case==1)then
      call compute_spline_nonunif( f, spl_per,node_positions)
      call interpolate_array_value_nonunif( Xstar(1:N), f, N, spl_per)
      return
    endif
    
    !from f compute the mean
    do i=0,N-1
      f(i+1)=f(i+1)*(node_positions(i+2)-node_positions(i+1))/dx
    enddo
    
    
    !we compute the splines coefficients by solving the LU decomposition
    M=0._f64
    do i=1,N
      M=M+f(i)
    enddo
    !M=M/real(N,f64)
    do i=1,N
      f(i)=f(i)-M*(node_positions(i+1)-node_positions(i))!/dx
    enddo
    tmp=0._f64
    do i=1,N
      tmp=tmp+f(i)
    enddo
    if(abs(tmp)>1.e-12)then
      print *,tmp
      stop    
    endif
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
    
    !node_positions = f
        
    call interpolate_array_value_nonunif( Xstar(1:N), f(1:N), N, spl_per)
    
    
    tmp=f(1)!;for(i=0;i<Nx-1;i++)p[i]=p[i+1]-p[i];p[Nx-1]=tmp+M-p[Nx-1];
    do i=1,N-1
      f(i)=f(i+1)-f(i)+M*(node_positions(i+1)-node_positions(i))
    enddo
    f(N)=tmp-f(N)+M*(node_positions(1)+1._f64-node_positions(N))


    !from mean compute f
    do i=1,N-1
      f(i)=f(i)*dx/(node_positions(i+1)-node_positions(i))
    enddo
    f(N) = f(N)*dx/(node_positions(1)+1._f64-node_positions(N))
    f(N+1) = f(1)

    !tmp=0._f64
    !do i=1,N
    !  tmp2=min(abs(node_positions(i)-Xstar(i)),abs(node_positions(i)+1._f64-Xstar(i)))
    !  tmp2=min(tmp2,abs(node_positions(i)-1._f64-Xstar(i)))
    !  if(tmp2>tmp)then
    !    tmp = tmp2
    !  endif
    !enddo
    !tmp=0._f64
    !do i=4,N-4
    !  tmp2=abs(node_positions(i)-Xstar(i)-(node_positions(i-1)-Xstar(i-1)))
    !  if(tmp2>tmp)then
    !    tmp = tmp2
    !  endif
    !enddo
    !if(tmp>0.01/real(N,f64))then
    !  do i=1,N
    !    print *,real(i-1,f64)/real(N,f64),node_positions(i),f(i),Xstar(i)
    !  enddo
    !  stop
    !endif
    
    
    
    
  end subroutine csl_advection_per




subroutine advect_classical_csl_1(dt,a1,f,geom_x,N_x1,N_x2,buf,node_positions_x1,&
Xstar,spl_per_x1,time_case,carac_position_case,interp_case)
!solve \partial_t f(t,x,y) + \partial_x(a1(x,y)f(t,x,y))=0 over dt
  sll_int32,intent(in) :: N_x1,N_x2,time_case(2),carac_position_case,interp_case
  sll_real64,intent(in) :: dt,geom_x(2,2)
  sll_real64,dimension(:),pointer :: buf,node_positions_x1,Xstar
  sll_real64,dimension(:,:),pointer:: a1,f
  type(cubic_nonunif_spline_1D), pointer :: spl_per_x1
  sll_int32 :: i1,i2,i1m1,ix,ix1,s,i
  sll_real64 :: x,result,tmp,x1_min,x1_max,x2_min,x2_max
  
  x1_min = geom_x(1,1)
  x1_max = geom_x(2,1)
  x2_min = geom_x(1,2)
  x2_max = geom_x(2,2)
  
  do i2=1,N_x2
    buf(1:N_x1) = f(1:N_x1,i2)
    if(time_case(1)==1)then
      do i1=1,N_x1
        Xstar(i1) = node_positions_x1(i1)-dt*a1(i1,i2)/(x1_max-x1_min)
        if(carac_position_case==-1)then
          i1m1=modulo(i1-1-1+N_x1,N_x1)+1    
          Xstar(i1) = node_positions_x1(i1)-0.5_f64*dt*(a1(i1,i2)+a1(i1m1,i2))/(x1_max-x1_min)
        endif
        if(carac_position_case==1)then
          i1m1=modulo(i1+1-1+N_x1,N_x1)+1    
          Xstar(i1) = node_positions_x1(i1)-0.5_f64*dt*(a1(i1,i2)+a1(i1m1,i2))/(x1_max-x1_min)
        endif
          
      enddo      
      Xstar(N_x1+1) = Xstar(1)+1._f64
    endif        
    if(time_case(1)==2)then
      do i1=1,N_x1
        !x=i*dx-p[i];if(x>=1.)x-=1.;if(x<0.)x+=1.;if(x>=1.)x-=1.;ix=(int)(x*Nx);
        !if(x>=1. || x<0.){fprintf(stderr,"x too big/small %1.1000lg i=%d j=%d\n",x,i,j);exit(1);}
        !assert(x>=0.&&x<1.);assert(ix>=0 &&ix<Nx);x=x*Nx-ix;assert(x>=0 &&x<1.);
        !ix1=ix+1;if(ix1==Nx)ix1=0;
        !result=(1.-x)*Ey[ix+Nx*j]+x*Ey[ix1+Nx*j];p[i]=result*dtmp2*dt;
        Xstar(i1) = 0._f64
        do s=1,time_case(2)
          x = node_positions_x1(i1)-Xstar(i1)
          if(abs(x-0.5_f64)>1.5)then
            print *,'#displacement is too big in classical_csl_1',x,s
            do i=1,N_x1
              print *,node_positions_x1(i),a1(i,i2)
            enddo
            print *,node_positions_x1(N_x1+1),a1(1,i2)
            stop
          endif
          do while (x>=1._f64)
            x=x-1._f64
          enddo
          do while(x<0._f64)
            x=x+1._f64
          enddo
          if(x>=1.)then
            x=x-1._f64
          endif
          if(x<0._f64)then
            print *,'problem x should be >=0 in advec_classical_csl_1',x
            stop
          endif
          x = x*N_x1
          ix=floor(x)
          x=x-real(ix,f64)
          ix1=ix+1
          if(ix1==N_x1)then
            ix1=0
          endif
          if(ix1+1>N_x1)then
            print *,ix1+1,i1,node_positions_x1(i1)-Xstar(i1),Xstar(i1)
            print *,'Problem in advec_classical_csl_1'
            stop
          endif
          result = (1.-x)*a1(ix+1,i2)+x*a1(ix1+1,i2)
          Xstar(i1)=0.5_f64*result*dt/(x1_max-x1_min)
        enddo  
        !i1m1=modulo(i1-1-1+N_x1,N_x1)+1
        !Xstar(i1) = node_positions_x1(i1)-Xstar(i1)
        !Xstar(i1) = node_positions_x1(i1)-0.5_f64*dt*(a1(i1,i2)+a1(i1m1,i2))
      enddo
      !tmp=1._f64/real(N_x1)Xstar(2)-Xstar(1)
      !do i=1,N_x1-1
      !  if(Xstar(i+1)-Xstar(i)<tmp)then
      !    tmp=Xstar(i+1)-Xstar(i)
      !  endif
      !enddo
      if(carac_position_case==-1)then
        tmp=Xstar(N_x1)
        do i1=N_x1,2,-1
          Xstar(i1)=0.5_f64*(Xstar(i1)+Xstar(i1-1))
        enddo
        Xstar(1)=0.5_f64*(Xstar(1)+tmp)
      endif  
      if(carac_position_case==1)then
        tmp=Xstar(1)
        do i1=1,N_x1-1
          Xstar(i1)=0.5_f64*(Xstar(i1)+Xstar(i1+1))
        enddo
        Xstar(N_x1)=0.5_f64*(Xstar(N_x1-1)+tmp)
      endif
      !in Xstar are stored displacements of X_{-1/2},...X_{N-1/2}: X_star=d_{-1/2},...,d_{N-1/2}
      ! the displacement of X_{N+1/2} is the same than the displacement of X_{-1/2}
      !the primitive will be evaluated at 0dx-d_{-1/2}, dx-d_{1/2},...(N-1)dx-d_{N-1/2}
      do i1=1,N_x1
        Xstar(i1)=node_positions_x1(i1)-2._f64*Xstar(i1)
      enddo
      Xstar(N_x1+1) = Xstar(1)+ 1._f64   
    endif
    if(time_case(1)==3)then
      do i1=1,N_x1
        Xstar(i1) = node_positions_x1(i1)-dt*a1(i1,i2)/(x1_max-x1_min)
        if(carac_position_case==-1)then
          i1m1=modulo(i1-1-1+N_x1,N_x1)+1    
          Xstar(i1) = node_positions_x1(i1)-0.5_f64/real(N_x1,f64)-0.5_f64*dt*(a1(i1,i2)+a1(i1m1,i2))/(x1_max-x1_min)
        endif
        if(carac_position_case==1)then
          i1m1=modulo(i1+1-1+N_x1,N_x1)+1    
          Xstar(i1) = node_positions_x1(i1)-0.5_f64/real(N_x1,f64)-0.5_f64*dt*(a1(i1,i2)+a1(i1m1,i2))/(x1_max-x1_min)
        endif
          
      enddo      
      Xstar(N_x1+1) = Xstar(1)+1._f64
    endif        
    if(time_case(1)==4)then
      do i1=1,N_x1
        !x=i*dx-p[i];if(x>=1.)x-=1.;if(x<0.)x+=1.;if(x>=1.)x-=1.;ix=(int)(x*Nx);
        !if(x>=1. || x<0.){fprintf(stderr,"x too big/small %1.1000lg i=%d j=%d\n",x,i,j);exit(1);}
        !assert(x>=0.&&x<1.);assert(ix>=0 &&ix<Nx);x=x*Nx-ix;assert(x>=0 &&x<1.);
        !ix1=ix+1;if(ix1==Nx)ix1=0;
        !result=(1.-x)*Ey[ix+Nx*j]+x*Ey[ix1+Nx*j];p[i]=result*dtmp2*dt;
        Xstar(i1) = 0._f64
        do s=1,time_case(2)
          x = node_positions_x1(i1)+real(carac_position_case,f64)*0.5_f64/real(N_x1,f64)-Xstar(i1)
          if(abs(x-0.5_f64)>1.5)then
            print *,'#displacement is too big in classical_csl_1',x,s
            do i=1,N_x1
              print *,node_positions_x1(i),a1(i,i2)
            enddo
            print *,node_positions_x1(N_x1+1),a1(1,i2)
            stop
          endif
          !do while (x>=1._f64)
          !  x=x-1._f64
          !enddo
          !do while(x<0._f64)
          !  x=x+1._f64
          !enddo
          !if(x>=1.)then
          !  x=x-1._f64
          !endif
          !if(x<0._f64)then
          !  print *,'problem x should be >=0 in advec_classical_csl_1',x
          !  stop
          !endif
          x = x*N_x1
          ix=floor(x)
          x=x-real(ix,f64)
          ix = modulo(ix+N_x1,N_x1)
          ix1=ix+1
          if(ix1==N_x1)then
            ix1=0
          endif
          if(ix1+1>N_x1)then
            print *,ix,ix1+1,i1,node_positions_x1(i1)-Xstar(i1),Xstar(i1)
            print *,'Problem in advec_classical_csl_1'
            stop
          endif
          result = (1.-x)*a1(ix+1,i2)+x*a1(ix1+1,i2)
          Xstar(i1)=0.5_f64*result*dt/(x1_max-x1_min)
        enddo  
        !i1m1=modulo(i1-1-1+N_x1,N_x1)+1
        !Xstar(i1) = node_positions_x1(i1)-Xstar(i1)
        !Xstar(i1) = node_positions_x1(i1)-0.5_f64*dt*(a1(i1,i2)+a1(i1m1,i2))
      enddo
      do i1=1,N_x1
        !Xstar(i1)=node_positions_x1(i1)+real(carac_position_case,f64)*0.5_f64/real(N_x1,f64)&
        !  -2._f64*Xstar(i1)
        Xstar(i1)=node_positions_x1(i1)-2._f64*Xstar(i1)
      enddo
      Xstar(N_x1+1) = Xstar(1)+ 1._f64   
    endif
    if(time_case(1)==10)then
      do i1=1,N_x1
        Xstar(i1) = node_positions_x1(i1)+dt*a1(i1,i2)/(x1_max-x1_min)
        if(carac_position_case==-1)then
          i1m1=modulo(i1-1-1+N_x1,N_x1)+1    
          Xstar(i1) = node_positions_x1(i1)+0.5_f64*dt*(a1(i1,i2)+a1(i1m1,i2))/(x1_max-x1_min)
        endif
        if(carac_position_case==1)then
          i1m1=modulo(i1+1-1+N_x1,N_x1)+1    
          Xstar(i1) = node_positions_x1(i1)+0.5_f64*dt*(a1(i1,i2)+a1(i1m1,i2))/(x1_max-x1_min)
        endif
          
      enddo      
      Xstar(N_x1+1) = Xstar(1)+1._f64
    endif        



    if(time_case(1)<10)then    
      call csl_advection_per(buf,spl_per_x1,Xstar,node_positions_x1,N_x1,interp_case)
    endif
    if(time_case(1)>=10)then
      call interp1dcascons(buf,Xstar,100,N_x1)
      buf(N_x1+1)=buf(1)
    endif  
    f(1:N_x1+1,i2) = buf(1:N_x1+1)
  enddo
  f(1:N_x1+1,N_x2+1)=f(1:N_x1+1,1)

end subroutine advect_classical_csl_1


subroutine advect_classical_csl_2(dt,a2,f,geom_x,N_x1,N_x2,buf,node_positions_x2,&
Xstar,spl_per_x2,time_case,carac_position_case,interp_case)
!solve \partial_t f(t,x,y) + \partial_x(a1(x,y)f(t,x,y))=0 over dt
  sll_int32,intent(in) :: N_x1,N_x2,time_case(2),carac_position_case,interp_case
  sll_real64,intent(in) :: dt,geom_x(2,2)
  sll_real64,dimension(:),pointer :: buf,node_positions_x2,Xstar
  sll_real64,dimension(:,:),pointer:: a2,f
  type(cubic_nonunif_spline_1D), pointer :: spl_per_x2
  sll_int32 :: i1,i2,i2m1,ix,ix1,s,i
  sll_real64 :: x,result,tmp,x1_min,x1_max,x2_min,x2_max

  x1_min = geom_x(1,1)
  x1_max = geom_x(2,1)
  x2_min = geom_x(1,2)
  x2_max = geom_x(2,2)


  do i1=1,N_x1
    buf(1:N_x2) = f(i1,1:N_x2)
    if(time_case(1)==1)then
      do i2=1,N_x2
        Xstar(i2) = node_positions_x2(i2)-dt*a2(i1,i2)/(x2_max-x2_min)
        if(carac_position_case==-1)then
          i2m1=modulo(i2-1-1+N_x2,N_x2)+1    
          Xstar(i2) = node_positions_x2(i2)-0.5_f64*dt*(a2(i1,i2)+a2(i1,i2m1))/(x2_max-x2_min)
        endif  
        if(carac_position_case==1)then
          i2m1=modulo(i2+1-1+N_x2,N_x2)+1    
          Xstar(i2) = node_positions_x2(i2)-0.5_f64*dt*(a2(i1,i2)+a2(i1,i2m1))/(x2_max-x2_min)
        endif  
      enddo
      Xstar(N_x2+1)=Xstar(1)+1._f64
    endif
    if(time_case(1)==2)then  
      do i2=1,N_x2
        i2m1=modulo(i2-1-1+N_x2,N_x2)+1
        !Xstar(i2) = node_positions_x2(i2)-0.5_f64*dt*(a2(i1,i2)+a2(i1,i2m1))
        Xstar(i2) = 0._f64
        do s=1,time_case(2)
          x = node_positions_x2(i2)-Xstar(i2)
          if(abs(x-0.5_f64)>1.5)then
            print *,'displacement is too big in classical_csl_2',x
            do i=1,N_x2
              print *,node_positions_x2(i),a2(i1,i)
            enddo
            print *,node_positions_x2(N_x1+1),a2(i1,i)
            stop
          endif
          do while(x>=1._f64)
            x=x-1._f64
          enddo
          do while(x<0._f64)
            x=x+1._f64
          enddo
          if(x>=1.)then
            x=x-1._f64
          endif
          if(x<0._f64)then
            print *,'problem x should be >=0 in advec_classical_csl_2',x
            stop
          endif
          x = x*N_x2
          ix=floor(x)
          x=x-real(ix,f64)
          ix1=ix+1
          if(ix1==N_x2)then
            ix1=0
          endif
          if(ix1+1>N_x2)then
            print *,ix1+1,i2,node_positions_x2(i2)-Xstar(i2),Xstar(i2)
            print *,'Problem in advec_classical_csl_2'
            stop
          endif
          result = (1.-x)*a2(i1,ix+1)+x*a2(i1,ix1+1)
          Xstar(i2)=0.5_f64*result*dt/(x2_max-x2_min)
        enddo  
      enddo
      if(carac_position_case==-1)then
        tmp=Xstar(N_x2)
        do i2=N_x2,2,-1
          Xstar(i2)=0.5_f64*(Xstar(i2)+Xstar(i2-1))
        enddo
        Xstar(1)=0.5_f64*(Xstar(1)+tmp)
      endif
      if(carac_position_case==1)then      
        tmp=Xstar(1)
        do i2=1,N_x2-1
          Xstar(i2)=0.5_f64*(Xstar(i2)+Xstar(i2+1))
        enddo
        Xstar(N_x2)=0.5_f64*(Xstar(N_x2-1)+tmp)
      endif
      !in Xstar are stored displacements of X_{-1/2},...X_{N-1/2}: X_star=d_{-1/2},...,d_{N-1/2}
      ! the displacement of X_{N+1/2} is the same than the displacement of X_{-1/2}
      !the primitive will be evaluated at 0dx-d_{-1/2}, dx-d_{1/2},...(N-1)dx-d_{N-1/2}
      do i2=1,N_x2
        Xstar(i2)=node_positions_x2(i2)-2._f64*Xstar(i2)
      enddo    
      Xstar(N_x2+1) = Xstar(1)+ 1._f64
      
      !do i=1,N_x2+1
      !  print *,node_positions_x2(i),Xstar(i),a2(i1,i)
      !enddo
      !stop
    endif
    if(time_case(1)==3)then
      do i2=1,N_x2
        Xstar(i2) = node_positions_x2(i2)-dt*a2(i1,i2)/(x2_max-x2_min)
        if(carac_position_case==-1)then
          i2m1=modulo(i2-1-1+N_x2,N_x2)+1    
          Xstar(i2) = node_positions_x2(i2)-0.5_f64/real(N_x2,f64)-0.5_f64*dt*(a2(i1,i2)+a2(i1,i2m1))/(x2_max-x2_min)
        endif  
        if(carac_position_case==1)then
          i2m1=modulo(i2+1-1+N_x2,N_x2)+1    
          Xstar(i2) = node_positions_x2(i2)+0.5_f64/real(N_x2,f64)-0.5_f64*dt*(a2(i1,i2)+a2(i1,i2m1))/(x2_max-x2_min)
        endif  
      enddo
      Xstar(N_x2+1)=Xstar(1)+1._f64
    endif

    if(time_case(1)==4)then  
      do i2=1,N_x2
        i2m1=modulo(i2-1-1+N_x2,N_x2)+1
        !Xstar(i2) = node_positions_x2(i2)-0.5_f64*dt*(a2(i1,i2)+a2(i1,i2m1))
        Xstar(i2) = 0._f64
        do s=1,time_case(2)
          x = node_positions_x2(i2)+real(carac_position_case,f64)*0.5_f64/real(N_x2,f64)-Xstar(i2)
          if(abs(x-0.5_f64)>1.5)then
            print *,'displacement is too big in classical_csl_2',x
            do i=1,N_x2
              print *,node_positions_x2(i),a2(i1,i)
            enddo
            print *,node_positions_x2(N_x1+1),a2(i1,i)
            stop
          endif
          !do while(x>=1._f64)
          !  x=x-1._f64
          !enddo
          !do while(x<0._f64)
          !  x=x+1._f64
          !enddo
          !if(x>=1.)then
          !  x=x-1._f64
          !endif
          !if(x<0._f64)then
          !  print *,'problem x should be >=0 in advec_classical_csl_2',x
          !  stop
          !endif
          x = x*N_x2
          ix=floor(x)
          x=x-real(ix,f64)
          ix = modulo(ix+N_x2,N_x2)
          ix1=ix+1
          if(ix1==N_x2)then
            ix1=0
          endif
          if(ix1+1>N_x2)then
            print *,ix,ix1+1,i2,node_positions_x2(i2)-Xstar(i2),Xstar(i2)
            print *,'Problem in advec_classical_csl_2'
            stop
          endif
          result = (1.-x)*a2(i1,ix+1)+x*a2(i1,ix1+1)
          Xstar(i2)=0.5_f64*result*dt/(x2_max-x2_min)
        enddo  
      enddo
      do i2=1,N_x2
        !Xstar(i2)=node_positions_x2(i2)+real(carac_position_case,f64)*0.5_f64/real(N_x2,f64)&
        !  -2._f64*Xstar(i2)
        Xstar(i2)=node_positions_x2(i2)-2._f64*Xstar(i2)
      enddo    
      Xstar(N_x2+1) = Xstar(1)+ 1._f64
      
      !do i=1,N_x2+1
      !  print *,node_positions_x2(i),Xstar(i),a2(i1,i)
      !enddo
      !stop
    endif
    if(time_case(1)==10)then
      do i2=1,N_x2
        Xstar(i2) = node_positions_x2(i2)+dt*a2(i1,i2)/(x2_max-x2_min)
        if(carac_position_case==-1)then
          i2m1=modulo(i2-1-1+N_x2,N_x2)+1    
          Xstar(i2) = node_positions_x2(i2)+0.5_f64*dt*(a2(i1,i2)+a2(i1,i2m1))/(x2_max-x2_min)
        endif  
        if(carac_position_case==1)then
          i2m1=modulo(i2+1-1+N_x2,N_x2)+1    
          Xstar(i2) = node_positions_x2(i2)+0.5_f64*dt*(a2(i1,i2)+a2(i1,i2m1))/(x2_max-x2_min)
        endif  
      enddo
      Xstar(N_x2+1)=Xstar(1)+1._f64
    endif
    
    
    
    if(time_case(1)<10)then
      call csl_advection_per(buf,spl_per_x2,Xstar,node_positions_x2,N_x2,interp_case)
    endif
    if(time_case(1)>=10)then
      call interp1dcascons(buf,Xstar,100,N_x2)
      buf(N_x2+1)=buf(1)
    endif
    f(i1,1:N_x2+1) = buf(1:N_x2+1)
  enddo
  f(N_x1+1,1:N_x2+1)=f(1,1:N_x2+1)

end subroutine advect_classical_csl_2

subroutine advect_classical_csl(dt,a1,a2,f,geom_x,N_x1,N_x2,buf,&
node_positions_x1,node_positions_x2,Xstar,spl_per_x1,spl_per_x2)
!solve \partial_t f(t,x,y) + \partial_x(a1(x,y)f(t,x,y))=0 over dt
  sll_int32,intent(in) :: N_x1,N_x2
  sll_real64,intent(in) :: dt,geom_x(2,2)
  sll_real64,dimension(:),pointer :: buf,node_positions_x1,node_positions_x2,Xstar
  sll_real64,dimension(:,:),pointer:: a1,a2,f
  type(cubic_nonunif_spline_1D), pointer :: spl_per_x1,spl_per_x2
  sll_int32  :: time_case(2),carac_position_case,interp_case
  time_case(1) = 10
  time_case(2) = 20
  carac_position_case = -1
  interp_case = 0
  
  !call advect_classical_csl_1(0.5_f64*dt,a1,f,geom_x,N_x1,N_x2,buf,node_positions_x1,&
  !  Xstar,spl_per_x1,time_case)
  !call advect_classical_csl_2(dt,a2,f,geom_x,N_x1,N_x2,buf,node_positions_x2,&
  !  Xstar,spl_per_x2,time_case)
  !call advect_classical_csl_1(0.5_f64*dt,a1,f,geom_x,N_x1,N_x2,buf,node_positions_x1,&
  !  Xstar,spl_per_x1,time_case)

  call advect_classical_csl_2(0.5_f64*dt,a2,f,geom_x,N_x1,N_x2,buf,node_positions_x2,&
    Xstar,spl_per_x2,time_case,carac_position_case,interp_case)
  call advect_classical_csl_1(dt,a1,f,geom_x,N_x1,N_x2,buf,node_positions_x1,&
    Xstar,spl_per_x1,time_case,carac_position_case,interp_case)
  call advect_classical_csl_2(0.5_f64*dt,a2,f,geom_x,N_x1,N_x2,buf,node_positions_x2,&
    Xstar,spl_per_x2,time_case,carac_position_case,interp_case)


end subroutine advect_classical_csl



subroutine interp1dcascons(pp,alphax,Nbdr,Nx)
  sll_int32,intent(in) :: Nbdr,Nx
  sll_real64,dimension(:),pointer :: pp,alphax
  sll_int32 :: i,j,ix,mem=0,err
  sll_real64 :: x,dx,tmp,M,w(0:3)
  sll_real64,dimension(:),pointer :: coef,Xnode,A,B,C,ltab2,vtab2,dtab2,mtab2,p
  
  dx = 1._f64/real(Nx,f64) 
  
  SLL_ALLOCATE(coef(-Nbdr:Nx+Nbdr-1),err)
  SLL_ALLOCATE(Xnode(-Nbdr:Nx+Nbdr-1),err)
  SLL_ALLOCATE(A(0:Nx-1),err)
  SLL_ALLOCATE(B(0:Nx-1),err)
  SLL_ALLOCATE(C(0:Nx-1),err)
  SLL_ALLOCATE(ltab2(0:Nx-2),err)
  SLL_ALLOCATE(vtab2(0:Nx-2),err)
  SLL_ALLOCATE(dtab2(0:Nx-1),err)
  SLL_ALLOCATE(mtab2(0:Nx-3),err)
  SLL_ALLOCATE(p(0:Nx-1),err)
  

  !do i=0,Nx-1
  !  print *,i,pp(i+1),alphax(i+1)!,p(i)
  !enddo
  !stop


  
  !coef(-Nbdr:Nx+Nbdr-1) X(-Nbdr:Nx+Nbdr-1)
  !A(0:Nx-1) B(0:Nx-1) C(0:Nx-1) ltab2(0:Nx-2) vtab2(0:Nx-2) dtab2(0:Nx-1) mtab2(0:Nx-3)
  
  !we compute the new 1D unstructured mesh
  !//for(i=0;i<Nx;i++)X[i]=((double)i-0.5)*dx+alphax[i];
  !for(i=0;i<Nx;i++)X[i]=alphax[i];
  do i=0,Nx-1
    Xnode(i) = alphax(i+1)
    p(i) = pp(i+1)
  enddo
  
  
  !//for(i=0;i<Nx;i++)fprintf(stderr,"X[%d] %lg\n",i,X[i]);
  !for(i=0;i<Nbdr;i++)X[Nx+i]=X[i]+1.;for(i=0;i<Nbdr;i++)X[-i-1]=X[Nx-i-1]-1.;
  do i=0,Nbdr-1
    Xnode(Nx+i)=Xnode(i)+1._f64
  enddo
  do i=0,Nbdr-1
    Xnode(-i-1)=Xnode(Nx-i-1)-1._f64
  enddo

  !do i=-Nbdr,Nx+Nbdr-1
  !  print *,i,Xnode(i),i*dx
  !enddo
  !stop
  !//compute the minsize
  !tmp=X[1]-X[0];for(i=1;i<Nx;i++)if(X[i+1]-X[i]<tmp)tmp=X[i+1]-X[i];
  tmp=Xnode(1)-Xnode(0)
  do i=1,Nx-1
    if(Xnode(i+1)-Xnode(i)<tmp)then
      tmp=Xnode(i+1)-Xnode(i)
    endif  
  enddo
  if(tmp<1.e-10_f64)then
    do i=0,Nx-1
      print *,i,Xnode(i)      
    enddo
    print *,"min size of forward mesh is too small:",tmp
  endif
  !print *,tmp,1._f64/real(Nx,f64),tmp-1./real(Nx,f64)
  !stop
  !if(tmp<1.e-10){
  !  for(i=0;i<Nx;i++)fprintf(stderr,"X[%d] %lg\n",i,X[i]);
  !  fprintf(stderr,"min size of forward mesh is too small:%lg\n",tmp);exit(1);
  !}
  !//we compute the almost tridiag matrix and LU decomposition	    
  do i=0,Nx-1
    B(i)=((Xnode(i+2)-Xnode(i+1))*(Xnode(i+2)-Xnode(i+1)))/((Xnode(i+2)-Xnode(i-1))*(Xnode(i+2)-Xnode(i)))
    C(i)=((Xnode(i-1)-Xnode(i-2))*(Xnode(i-1)-Xnode(i-2)))/((Xnode(i+1)-Xnode(i-2))*(Xnode(i)-Xnode(i-2)))
    A(i)=((Xnode(i)-Xnode(i-2))*(Xnode(i+1)-Xnode(i)))/((Xnode(i+1)-Xnode(i-2))*(Xnode(i+1)-Xnode(i-1)))
    A(i)=A(i)+((Xnode(i+2)-Xnode(i))*(Xnode(i)-Xnode(i-1)))/((Xnode(i+2)-Xnode(i-1))*(Xnode(i+1)-Xnode(i-1)))
  enddo
  
  !do i=0,Nx-1
  !  print *,i,A(i),B(i),C(i)
  !enddo
  !stop
  
  dtab2(0)=A(0);vtab2(0)=B(Nx-1);ltab2(0)=B(0)/dtab2(0);mtab2(0)=C(0)/dtab2(0);
  do i=0,Nx-3
    ltab2(i)=B(i)/dtab2(i);dtab2(i+1)=A(i+1)-ltab2(i)*C(i+1)
  enddo
  do i=0,Nx-4
    vtab2(i+1)=-ltab2(i)*vtab2(i);mtab2(i+1)=-mtab2(i)*C(i+1)/dtab2(i+1);
  enddo  
  vtab2(Nx-2)=C(Nx-1)-ltab2(Nx-3)*vtab2(Nx-3);ltab2(Nx-2)=(B(Nx-2)-mtab2(Nx-3)*C(Nx-2))/dtab2(Nx-2);
  tmp=0.;
  do i=0,Nx-3
    tmp=tmp+mtab2(i)*vtab2(i)
  enddo  
  dtab2(Nx-1)=A(Nx-1)-tmp-ltab2(Nx-2)*vtab2(Nx-2);
  !//we compute the coefficients in x by solving the LU decomposition 
  M=0._f64
  do i=0,Nx-1
    M=M+p(i)
  enddo
  
  
  coef(0)=0._f64
  do i=1,Nx-1
    coef(i)=coef(i-1)+p(i-1);
  enddo
  coef(0)=coef(0)+M*B(Nx-1)
  coef(Nx-1)=coef(Nx-1)-M*C(0)
  !print *,B(Nx-1),coef(0),M
  
  !stop

  do i=1,Nx-2
    coef(i)=coef(i)-ltab2(i-1)*coef(i-1);
  enddo
  tmp=0._f64;
  do i=0,Nx-3
    tmp=tmp+coef(i)*mtab2(i)
  enddo
  coef(Nx-1)=coef(Nx-1)-(tmp+ltab2(Nx-2)*coef(Nx-2))
  coef(Nx-1)=coef(Nx-1)/dtab2(Nx-1);coef(Nx-2)=(coef(Nx-2)-coef(Nx-1)*vtab2(Nx-2))/dtab2(Nx-2);



  do i=Nx-3,0,-1
    coef(i)=(coef(i)-C(i+1)*coef(i+1)-vtab2(i)*coef(Nx-1))/dtab2(i);
  enddo
  do i=0,Nbdr-1
    coef(Nx+i)=coef(i)+M;
  enddo
  do i=0,Nbdr-1
    coef(-i-1)=coef(Nx-i-1)-M;
  enddo
  
  
  !do i=-Nbdr,Nx+Nbdr-1
  !  print *,i,coef(i)
  !enddo
  !stop  
  !//we interpolate the cumulative function on uniform mesh and get     
  do i=0,Nx-1
    tmp=(i-0.0_f64)*dx;j=i;
    if(Xnode(j)<tmp)then
      do while(Xnode(j)<tmp)
        j=j+1
      enddo
      j=j-1
    else
      do while (Xnode(j)>tmp)
        j=j-1
      enddo
    endif    
    if(.not.((Xnode(j)<=tmp) .and. (Xnode(j+1)>tmp)))then
      print *,"j= tmp= Xnode(j)= Xnode(j+1)= ",j,tmp,Xnode(j),Xnode(j+1),Xnode(j+1)-tmp;
      stop
    endif  
    !assert(Xnode(j)<=tmp && Xnode(j+1)>tmp);
    if(.not.((j>=-Nbdr).and.(j<Nx+Nbdr)))then
      print *,"j= Nx= Nbdr= ",j,Nx,Nbdr;
      stop
    endif
    !assert(j>=-Nbdr && j<Nx+Nbdr);
    w(0)=(Xnode(j+1)-tmp)*(Xnode(j+1)-tmp)*(Xnode(j+1)-tmp)/((Xnode(j+1)-Xnode(j))&
      *(Xnode(j+1)-Xnode(j-1))*(Xnode(j+1)-Xnode(j-2)));    
    w(1)=(Xnode(j+1)-tmp)*(Xnode(j+1)-tmp)*(tmp-Xnode(j-2))/((Xnode(j+1)-Xnode(j))&
      *(Xnode(j+1)-Xnode(j-1))*(Xnode(j+1)-Xnode(j-2)));
    w(1)=w(1)+(Xnode(j+2)-tmp)*(Xnode(j+1)-tmp)*(tmp-Xnode(j-1))/((Xnode(j+1)-Xnode(j))&
      *(Xnode(j+1)-Xnode(j-1))*(Xnode(j+2)-Xnode(j-1)));
    w(1)=w(1)+(Xnode(j+2)-tmp)*(Xnode(j+2)-tmp)*(tmp-Xnode(j))/((Xnode(j+1)-Xnode(j))&
      *(Xnode(j+2)-Xnode(j))*(Xnode(j+2)-Xnode(j-1)));	 
    w(2)=(Xnode(j+1)-tmp)*(tmp-Xnode(j-1))*(tmp-Xnode(j-1))/((Xnode(j+1)-Xnode(j))&
      *(Xnode(j+1)-Xnode(j-1))*(Xnode(j+2)-Xnode(j-1)));
    w(2)=w(2)+(Xnode(j+2)-tmp)*(tmp-Xnode(j-1))*(tmp-Xnode(j))/((Xnode(j+1)-Xnode(j))&
      *(Xnode(j+2)-Xnode(j))*(Xnode(j+2)-Xnode(j-1)));
    w(2)=w(2)+(Xnode(j+3)-tmp)*(tmp-Xnode(j))*(tmp-Xnode(j))/((Xnode(j+1)-Xnode(j))&
      *(Xnode(j+2)-Xnode(j))*(Xnode(j+3)-Xnode(j)));    
    w(3)=(tmp-Xnode(j))*(tmp-Xnode(j))*(tmp-Xnode(j))/((Xnode(j+1)-Xnode(j))&
      *(Xnode(j+2)-Xnode(j))*(Xnode(j+3)-Xnode(j)));
    p(i)=w(0)*coef(j-1)+w(1)*coef(j)+w(2)*coef(j+1)+w(3)*coef(j+2);
  enddo
  
  !do i=0,Nx-1
  !  print *,i,p(i)
  !enddo
  
  !stop
  tmp=p(0);
  do i=0,Nx-2
    p(i)=p(i+1)-p(i)
  enddo
  p(Nx-1)=tmp+M-p(Nx-1);
  
  
  !do i=0,Nx-1
  !  print *,i,pp(i+1),alphax(i+1),p(i)
  !enddo
  !stop
  
  do i=0,Nx-1
    pp(i+1)=p(i)
  enddo  

  SLL_DEALLOCATE(coef,err)
  SLL_DEALLOCATE(Xnode,err)
  SLL_DEALLOCATE(A,err)
  SLL_DEALLOCATE(B,err)
  SLL_DEALLOCATE(C,err)
  SLL_DEALLOCATE(ltab2,err)
  SLL_DEALLOCATE(vtab2,err)
  SLL_DEALLOCATE(dtab2,err)
  SLL_DEALLOCATE(mtab2,err)
  SLL_DEALLOCATE(p,err)
  
  
end subroutine interp1dcascons
!! begin poisson module


  subroutine poisson2dpersize(Nbuf,Nbufc,Nx,Ny)  
    sll_int32,intent(in)::Nx,Ny
    sll_int32,intent(out)::Nbuf,Nbufc(2)
    Nbuf=2*Nx+15+4*Ny+15
    Nbufc(1)=Ny;Nbufc(2)=Nx+2
  end subroutine poisson2dpersize

  subroutine poisson2dperalloc(buf,Nbuf,bufc,Nbufc,Nx,Ny)  
    sll_int32,intent(in)::Nx,Ny
    sll_int32,intent(out)::Nbuf,Nbufc(2)
    sll_real64,dimension(:),allocatable::buf    
    complex(f64),dimension(:,:,:),allocatable::bufc
    call poisson2dpersize(Nbuf,Nbufc,Nx,Ny)
    allocate(buf(Nbuf),bufc(2,Nbufc(1),Nbufc(2)))       
  end subroutine poisson2dperalloc

  subroutine poisson2dperinit(buf,bufc,Nx,Ny)  
    sll_int32,intent(in)::Nx,Ny!,Nbuf,Nbufc0,Nbufc1
    sll_real64,dimension(2*Nx+15+4*Ny+15),intent(out)::buf
    complex(f64),dimension(2,Ny,Nx+2),intent(out)::bufc
    call  dffti(Nx,buf(1:2*Nx+15))
    call  zffti(Ny,buf(2*Nx+16:2*Nx+15+4*Ny+15))
    bufc=cmplx(0._f64,0._f64)
  end subroutine poisson2dperinit


  subroutine poisson2dper(buf,Nbuf,bufc,Nbufc1,Nbufc2,Ex,Ey,Nx,Ny,Lx,Ly)
    sll_int32,intent(in)::Nbuf,Nbufc1,Nbufc2,Nx,Ny
    sll_real64,dimension(Nbuf),intent(in)::buf
    complex(f64),dimension(2,Nbufc1,Nbufc2),intent(inout)::bufc
    !real(f64),dimension(NEx0,NEx1),intent(inout)::Ex
    sll_real64,dimension(Nx,Ny)::Ex,Ey
    !f2py intent(in)::Ex,Ey
    sll_real64,intent(in)::Lx,Ly
    sll_int32::i,j
    sll_real64::re,im,tmp,kx0,ky0,kx,ky,kx2,k2

    do i=1,Ny
      call dfftf(Nx,Ex(:,i),buf(1:2*Nx+15))
    end do
    do j=1,Ny
      bufc(1,j,1)=cmplx(Ex(1,j),0._f64)
    enddo  
    do i=2,Nx/2
       do j=1,Ny
          bufc(1,j,i)=cmplx(Ex(2*i-2,j),Ex(2*i-1,j))
       end do
    end do
    do j=1,Ny
      bufc(1,j,Nx/2+1)=cmplx(Ex(2*(Nx/2),j),0._f64)
    enddo
 
    do i=1,Nx/2+1
      call zfftf(Ny,bufc(1,1:Ny,i),buf(2*Nx+16:Nbuf))
    enddo  

    
    ! calcul de la transformee de Fourier de E a partir de celle de rho
    kx0=2._f64*sll_pi/Lx
    ky0=2._f64*sll_pi/Ly
    bufc(1,1,1)=cmplx(0._f64,0._f64)
    bufc(2,1,1)=cmplx(0._f64,0._f64)
    do i=2,Nx/2+1
      !if(i/=1)then
      kx=real(i-1,f64)*kx0
      kx2=kx*kx
      k2=kx2*(real(Nx,f64)*real(Ny,f64))
      bufc(1,1,i)=-cmplx(0._f64,kx/k2)*bufc(1,1,i)
      bufc(2,1,i)=cmplx(0._f64,0._f64)
    enddo
    do i=1,Nx/2+1
      !if(i/=1)then
      kx=real(i-1,f64)*kx0
      kx2=kx*kx
      !k2=kx2*(real(Nx,f64)*real(Ny,f64))
      !bufc(1,1,i)=-cmplx(0._f64,kx/k2)*bufc(1,1,i)
      !bufc(2,1,i)=cmplx(0._f64,0._f64)
      do j=2,Ny/2+1
        ky=real(j-1,f64)*ky0
        k2=kx2+ky*ky
        k2=k2*real(Nx,f64)*real(Ny,f64)
        bufc(2,j,i)=-cmplx(0._f64,ky/k2)*bufc(1,j,i)
        bufc(1,j,i)=-cmplx(0._f64,kx/k2)*bufc(1,j,i)
      enddo
      do j=Ny/2+2,Ny
        ky=real(j-1-Ny,f64)*ky0
        k2=kx2+ky*ky
        k2=k2*real(Nx,f64)*real(Ny,f64)
        bufc(2,j,i)=-cmplx(0._f64,ky/k2)*bufc(1,j,i)
        bufc(1,j,i)=-cmplx(0._f64,kx/k2)*bufc(1,j,i)
      enddo      
    enddo

    do i=1,Nx/2+1
      call zfftb(Ny,bufc(1,1:Ny,i),buf(2*Nx+16:Nbuf))
      call zfftb(Ny,bufc(2,1:Ny,i),buf(2*Nx+16:Nbuf))
    enddo  

    do j=1,Ny
      Ex(1,j)=real(bufc(1,j,1))
      Ey(1,j)=real(bufc(2,j,1))
      do i=2,Nx/2
        Ex(2*i-2,j)=real(bufc(1,j,i))
        Ex(2*i-1,j)=aimag(bufc(1,j,i))
        Ey(2*i-2,j)=real(bufc(2,j,i))
        Ey(2*i-1,j)=aimag(bufc(2,j,i))
      enddo
      Ex(Nx,j)=real(bufc(1,j,Nx/2+1))
      Ey(Nx,j)=real(bufc(2,j,Nx/2+1))
    end do
    do i=1,Ny
      call dfftb(Nx,Ex(:,i),buf(1:2*Nx+15))
      call dfftb(Nx,Ey(:,i),buf(1:2*Nx+15))
    end do
           
  end subroutine poisson2dper

  subroutine computephipersize(Nbuf,Nbufc,Nx,Ny)  
    sll_int32,intent(in)::Nx,Ny
    sll_int32,intent(out)::Nbuf,Nbufc(2)!,Nbuf,Nbufc0,Nbufc1
    !real(f64),dimension(Nbuf),intent(out)::buf
    !complex(f64),dimension(2,Nbufc0,Nbufc1),intent(out)::bufc
    Nbuf=2*Nx+15+4*Ny+15
    Nbufc(1)=Ny;Nbufc(2)=Nx+2
  end subroutine computephipersize


  subroutine computephiperalloc(buf,bufc,Nx,Ny)  
    sll_int32,intent(in)::Nx,Ny!,Nbuf,Nbufc0,Nbufc1
    !real(f64),dimension(Nbuf),intent(out)::buf
    !complex(f64),dimension(2,Nbufc0,Nbufc1),intent(out)::bufc
    sll_real64,dimension(:),allocatable,intent(out)::buf
    complex(f64),dimension(:,:),allocatable,intent(out)::bufc
    allocate(buf(2*Nx+15+4*Ny+15),bufc(Ny,Nx+2))
  end subroutine computephiperalloc

  subroutine computephiperinit(buf,bufc,Nx,Ny)  
    sll_int32,intent(in)::Nx,Ny!,Nbuf,Nbufc0,Nbufc1
    !real(f64),dimension(Nbuf),intent(out)::buf
    !complex(f64),dimension(2,Nbufc0,Nbufc1),intent(out)::bufc
    sll_real64,dimension(2*Nx+15+4*Ny+15),intent(out)::buf
    complex(f64),dimension(Ny,Nx+2),intent(out)::bufc
    call  dffti(Nx,buf(1:2*Nx+15))
    call  zffti(Ny,buf(2*Nx+16:2*Nx+15+4*Ny+15))
    bufc=cmplx(0._f64,0._f64)
  end subroutine computephiperinit



  subroutine computephiper(buf,Nbuf,bufc,Nbufc1,Nbufc2,Phi,Nx,Ny,Lx,Ly)
    sll_int32,intent(in)::Nbuf,Nbufc1,Nbufc2,Nx,Ny
    sll_real64,dimension(Nbuf),intent(in)::buf
    complex(f64),dimension(Nbufc1,Nbufc2),intent(inout)::bufc
    !real(f64),dimension(NEx0,NEx1),intent(inout)::Ex
    sll_real64,dimension(Nx,Ny)::Phi
    !f2py intent(in)::Phi
    sll_real64,intent(in)::Lx,Ly
    sll_int32::i,j
    sll_real64::re,im,tmp,kx0,ky0,kx,ky,kx2,k2
    
    !print *,Nx
    do i=1,Ny
      call dfftf(Nx,Phi(:,i),buf(1:2*Nx+15))
    end do
    do j=1,Ny
      bufc(j,1)=cmplx(Phi(1,j),0._f64)
    enddo  
    do i=2,Nx/2
       do j=1,Ny
          bufc(j,i)=cmplx(Phi(2*i-2,j),Phi(2*i-1,j))
       end do
    end do
    do j=1,Ny
      bufc(j,Nx/2+1)=cmplx(Phi(2*(Nx/2),j),0._f64)
    enddo
 
    do i=1,Nx/2+1
      call zfftf(Ny,bufc(1:Ny,i),buf(2*Nx+16:Nbuf))
    enddo  

    
    ! calcul de la transformee de Fourier de E a partir de celle de rho
    kx0=2._f64*sll_pi/Lx
    ky0=2._f64*sll_pi/Ly
    bufc(1,1)=cmplx(0._f64,0._f64)
    do i=2,Nx/2+1
      !if(i/=1)then
      kx=real(i-1,f64)*kx0
      kx2=kx*kx
      k2=kx2*(real(Nx,f64)*real(Ny,f64))
      bufc(1,i)=-1._f64/k2*bufc(1,i)
    enddo
    do i=1,Nx/2+1
      !if(i/=1)then
      kx=real(i-1,f64)*kx0
      kx2=kx*kx
      !k2=kx2*(real(Nx,f64)*real(Ny,f64))
      !bufc(1,1,i)=-cmplx(0._f64,kx/k2)*bufc(1,1,i)
      !bufc(2,1,i)=cmplx(0._f64,0._f64)
      do j=2,Ny/2+1
        ky=real(j-1,f64)*ky0
        k2=kx2+ky*ky
        k2=k2*real(Nx,f64)*real(Ny,f64)
        bufc(j,i)=-1._f64/k2*bufc(j,i)
      enddo
      do j=Ny/2+2,Ny
        ky=real(j-1-Ny,f64)*ky0
        k2=kx2+ky*ky
        k2=k2*real(Nx,f64)*real(Ny,f64)
        bufc(j,i)=-1._f64/k2*bufc(j,i)
      enddo      
    enddo

    do i=1,Nx/2+1
      call zfftb(Ny,bufc(1:Ny,i),buf(2*Nx+16:Nbuf))
    enddo  

    do j=1,Ny
      Phi(1,j)=real(bufc(j,1))
      do i=2,Nx/2
        Phi(2*i-2,j)=real(bufc(j,i))
        Phi(2*i-1,j)=aimag(bufc(j,i))
      enddo
      Phi(Nx,j)=real(bufc(j,Nx/2+1))
    end do
    do i=1,Ny
      call dfftb(Nx,Phi(:,i),buf(1:2*Nx+15))
    end do
           
  end subroutine computephiper



  subroutine splpoissonperper2dsize(Nbufsize,Nx,Ny)
    sll_int32,intent(in)::Nx,Ny
    sll_int32,intent(out)::Nbufsize
    Nbufsize=max(Nx,Ny)+3*Nx+3*Ny
  end subroutine splpoissonperper2dsize
  
  subroutine splpoissonperper2dalloc(buf,Nx,Ny)
    sll_int32,intent(in)::Nx,Ny
    real(f64),dimension(:),allocatable,intent(out)::buf
    allocate(buf(max(Nx,Ny)+3*Nx+3*Ny))
  end subroutine splpoissonperper2dalloc

  subroutine splpoissonperper2dinit(buf,Nx,Ny)
    sll_int32,intent(in)::Nx,Ny
    sll_real64,dimension(max(Nx,Ny)+3*Nx+3*Ny),intent(out)::buf
    buf=0._f64
    call splcoefper1d0new(buf(max(Nx,Ny)+1:max(Nx,Ny)+3*Nx),Nx)
    call splcoefper1d0new(buf(max(Nx,Ny)+3*Nx+1:max(Nx,Ny)+3*Nx+3*Ny),Ny)    
  end subroutine splpoissonperper2dinit

  subroutine splpoissonperper2d(Ex,Ey,buf,sizebuf,Nx,Ny,Lx,Ly)
    sll_int32,intent(in)::Nx,Ny,sizebuf
    sll_real64,intent(in)::Lx,Ly
    sll_real64,dimension(0:Nx-1,0:Ny-1)::Ex,Ey
    sll_real64,dimension(0:sizebuf-1)::buf
    !f2py intent(in)::buf,Ex,Ey
    sll_int32::i,j
    sll_real64::dx,dy
    dx=Lx/real(Nx,f64);dy=Ly/real(Ny,f64)
    Ey=Ex    
    !periodic spline coefficients in x for Ex   
    do j=0,Ny-1
      do i=0,Nx-1;buf(i)=-(Ex(mod(i+Nx-1,Nx),j)-Ex(mod(i+1,Nx),j))*0.5_f64/dx;enddo
      call splcoefper1dnew(buf(0:Nx-1),buf(max(Nx,Ny):max(Nx,Ny)+3*Nx-1),Nx)
      do i=0,Nx-1;Ex(i,j)=buf(i);enddo      
    enddo
    !periodic spline coefficients in y for Ey   
    do i=0,Nx-1
      do j=0,Ny-1;buf(j)=-(Ey(i,mod(j+Ny-1,Ny))-Ey(i,mod(j+1,Ny)))*0.5_f64/dy;enddo
      call splcoefper1dnew(buf(0:Ny-1),buf(max(Nx,Ny)+3*Nx:sizebuf-1),Ny)
      do j=0,Ny-1;Ey(i,j)=buf(j);enddo      
    enddo
  end subroutine splpoissonperper2d


  subroutine init(E,N0,N1,k0,k1,testcase)
    sll_int32,intent(in)::N0,N1,k0,k1,testcase
    sll_real64,dimension(N0,N1),intent(out)::E
    sll_int32::i,j
    sll_real64::x,y    
    E=0._f64
    do j=1,N1
      do i=1,N0
        x=real(i-1,f64)*2._f64*sll_pi/real(N0,f64)
        y=real(j-1,f64)*2._f64*sll_pi/real(N1,f64)
        select case(testcase)
      case(1)
        E(i,j)=sin(k0*x)*sin(k1*y)
          case(2)
        E(i,j)=cos(k0*x)*sin(k1*y)
          case(3)
        E(i,j)=sin(k0*x)*cos(k1*y)
          case(4)
        E(i,j)=cos(k0*x)*cos(k1*y)
        end select    
      enddo
    enddo
    
  end subroutine init

  subroutine init2(E0,E1,N0,N1,k0,k1,testcase,L0,L1)
    sll_int32,intent(in)::N0,N1,k0,k1,testcase
    sll_real64,dimension(N0,N1),intent(out)::E0,E1
    sll_real64,intent(in)::L0,L1
    sll_int32::i,j
    sll_real64::x,y,kL0,kL1
    kL0=2._f64*sll_pi/L0*real(k0,f64)    
    kL1=2._f64*sll_pi/L1*real(k1,f64)    
    do j=1,N1
      do i=1,N0
        x=real(i-1,f64)*2._f64*sll_pi/real(N0,f64)
        y=real(j-1,f64)*2._f64*sll_pi/real(N1,f64)
        select case(testcase)
      case(1)
        E0(i,j)=-kL0/(kL0*kL0+kL1*kL1)*cos(k0*x)*sin(k1*y)
        E1(i,j)=-kL1/(kL0*kL0+kL1*kL1)*sin(k0*x)*cos(k1*y)
          case(2)
        E0(i,j)=kL0/(kL0*kL0+kL1*kL1)*sin(k0*x)*sin(k1*y)
        E1(i,j)=-kL1/(kL0*kL0+kL1*kL1)*cos(k0*x)*cos(k1*y)
          case(3)
        E0(i,j)=-kL0/(kL0*kL0+kL1*kL1)*cos(k0*x)*cos(k1*y)
        E1(i,j)=kL1/(kL0*kL0+kL1*kL1)*sin(k0*x)*sin(k1*y)
          case(4)
        E0(i,j)=kL0/(kL0*kL0+kL1*kL1)*sin(k0*x)*cos(k1*y)
        E1(i,j)=kL1/(kL0*kL0+kL1*kL1)*cos(k0*x)*sin(k1*y)
        end select    
      enddo
    enddo
    
  end subroutine init2

  

  subroutine splcoefper1d0new(luper,N)
    sll_int32,intent(in)::N
    sll_real64,dimension(0:3*N-1),intent(out)::luper
    sll_int32::i
    
    luper(0+3*0)=4._f64
    luper(2+3*0)=0.25_f64
    do i=0,N-2
      luper(1+3*i)=1._f64/luper(0+3*i)
      luper(0+3*(i+1))=4._f64-luper(1+3*i)
      luper(2+3*(i+1))=-luper(2+3*i)/luper(0+3*(i+1))
    enddo
    luper(0+3*(N-1))=luper(0+3*(N-1))-(luper(1+3*(N-2))+2._f64*luper(2+3*(N-2)))  
    do i=0,N-1
      luper(0+3*i)=1._f64/luper(0+3*i)
    enddo
  end subroutine splcoefper1d0new



  subroutine splcoefper1dnew(f,luper,N)
    sll_int32,intent(in)::N
    sll_real64,dimension(0:3*N-1),intent(in)::luper
    sll_real64,dimension(0:N-1),intent(inout)::f
    sll_int32::i
    do i=0,N-1;f(i)=6._f64*f(i);enddo;
    do i=1,N-1
      f(i)=f(i)-f(i-1)*luper(1+3*(i-1))
    enddo
    do i=0,N-2
      f(N-1)=f(N-1)-luper(2+3*i)*f(i)
    enddo
    f(N-1)=f(N-1)*luper(0+3*(N-1));f(N-2)=luper(0+3*(N-2))*(f(N-2)-(1._f64-luper(2+3*(N-3)))*f(N-1))
    do i=N-3,1,-1
      f(i)=luper(0+3*i)*(f(i)-f(i+1)+luper(2+3*(i-1))*f(N-1))
    enddo
    f(0)=luper(0+3*0)*(f(0)-f(1)-f(N-1));
  end subroutine splcoefper1dnew


  subroutine splcoefperper2dsize(sizebuf,Nx,Ny)
    sll_int32,intent(in)::Nx,Ny
    sll_int32,intent(out)::sizebuf
    sll_int32::N
    N=max(Nx,Ny)
    sizebuf=N+3*Nx+3*Ny    
  end subroutine splcoefperper2dsize

  subroutine splcoefperper2dinit(buf,Nx,Ny)
    sll_int32,intent(in)::Nx,Ny
    sll_real64,dimension(max(Nx,Ny)+3*Nx+3*Ny),intent(out)::buf
    buf=0._f64
    call splcoefper1d0new(buf(max(Nx,Ny)+1:max(Nx,Ny)+3*Nx),Nx)
    call splcoefper1d0new(buf(max(Nx,Ny)+3*Nx+1:max(Nx,Ny)+3*Nx+3*Ny),Ny)    
  end subroutine splcoefperper2dinit

  
  subroutine splcoefperper2d(f,buf,sizebuf,Nx,Ny)
    sll_int32,intent(in)::Nx,Ny,sizebuf
    sll_real64,dimension(Nx,Ny)::f
    sll_real64,dimension(sizebuf)::buf
    !f2py intent(in)::buf,f
    sll_int32::i,j    
    !periodic spline coefficients in x    
    do j=1,Ny
      do i=1,Nx;buf(i)=f(i,j);enddo
      call splcoefper1dnew(buf(1:Nx),buf(max(Nx,Ny)+1:max(Nx,Ny)+3*Nx),Nx)
      do i=1,Nx;f(i,j)=buf(i);enddo      
    enddo
    !periodic spline coefficients in y    
    do i=1,Nx
      do j=1,Ny;buf(j)=f(i,j);enddo
      call splcoefper1dnew(buf(1:Ny),buf(max(Nx,Ny)+3*Nx+1:sizebuf),Ny)
      do j=1,Ny;f(i,j)=buf(j);enddo      
    enddo
  end subroutine splcoefperper2d


  subroutine advect2dalloc(buf2d,Nbuf2d,buf,Nbuf,Nx,Ny)
    sll_int32,intent(in)::Nx,Ny
    sll_int32,intent(out)::Nbuf2d(2),Nbuf
    sll_real64,dimension(:,:),allocatable,intent(out)::buf2d
    sll_real64,dimension(:),allocatable,intent(out)::buf
    Nbuf2d(1)=Nx;Nbuf2d(2)=Ny
    Nbuf=max(Nx,Ny)+3*Nx+3*Ny
    allocate(buf2d(Nx,Ny),buf(max(Nx,Ny)+3*Nx+3*Ny))
  end subroutine advect2dalloc

  subroutine advect2dinit(buf2d,buf,Nx,Ny)
    sll_int32,intent(in)::Nx,Ny
    sll_real64,dimension(Nx,Ny),intent(out)::buf2d
    sll_real64,dimension(max(Nx,Ny)+3*Nx+3*Ny),intent(out)::buf
    buf=0._f64
    buf2d=0._f64
    call splcoefper1d0new(buf(max(Nx,Ny)+1:max(Nx,Ny)+3*Nx),Nx)
    call splcoefper1d0new(buf(max(Nx,Ny)+3*Nx+1:max(Nx,Ny)+3*Nx+3*Ny),Ny)    
    
  end subroutine advect2dinit

  subroutine advect2d(dom,f,E0,E1,N0,N1,buf,sizebuf,buf2d,dt,timecase,eps)
    sll_int32,intent(in)::N0,N1,sizebuf,timecase(0:1)
    sll_real64,dimension(0:1,0:1),intent(in)::dom
    sll_real64,intent(in)::dt,eps
    sll_real64,dimension(0:N0-1,0:N1-1)::buf2d,f,E0,E1
    !f2py intent(in)::buf2d,buf,f,E0,E1
    sll_real64,dimension(sizebuf)::buf
    sll_int32::i,j,ii,nstep=10,ix,iy
    sll_real64::fval,xx,yy,xx0,yy0,xxn,yyn,x,y,ax,ay
    
    buf2d=f
    call splcoefperper2d(buf2d,buf,sizebuf,N0,N1)
    !print*,buf2d
    !stop
    select case(timecase(0))
      case(1) !Euler
        do j=0,N1-1
          do i=0,N0-1
            xx=dom(0,0)+real(i,f64)*dom(1,0)/real(N0,f64)-E0(i,j)*dt
            yy=dom(0,1)+real(j,f64)*dom(1,1)/real(N1,f64)-E1(i,j)*dt
	    !print *,i,j,E0(i,j)!xx/dom(1,0)        
            call splperper2d(buf2d,xx,dom(0,0),dom(0,0)+dom(1,0),yy,dom(0,1),dom(0,1)+dom(1,1),f(i,j),N0,N1)
          enddo
        enddo
      case(2) !symplectic Euler
        do j=0,N1-1
          do i=0,N0-1
            xx0=dom(0,0)+real(i,f64)*dom(1,0)/real(N0,f64)!-E0(i,j)*dt
            yy0=dom(0,1)+real(j,f64)*dom(1,1)/real(N1,f64)!-E1(i,j)*dt        
            x=0._f64
	    ix=i
	    xx=xx0-E0(i,j)*dt
	    do ii=1,timecase(1)
	      x=(xx-dom(0,0))/dom(1,0)
              x=x-real(floor(x),f64)
              x=x*real(N0,f64)
              ix=floor(x)
              x=x-real(ix,f64)
              if(ix==N0)then
                ix=0
                x=0._f64
              endif
	      xxn=xx
              xx=xx0-dt*((1._f64-x)*E0(ix,j)+x*E0(mod(ix+1,N0),j))
            enddo  
	    if(abs(xxn-xx)>eps)then
	      print *,ii,xx,xxn,xx-xxn,ix,mod(ix+1,N0),x,E0(ix,j),E0(mod(ix+1,N0),j)
	      stop
	    endif
            yy=yy0-dt*((1._f64-x)*E1(ix,j)+x*E1(mod(ix+1,N0),j))  
	    call splperper2d(buf2d,xx,dom(0,0),dom(0,0)+dom(1,0),yy,dom(0,1),dom(0,1)+dom(1,1),f(i,j),N0,N1)
          enddo
        enddo	
      case(3) !symplectic Verlet
	do j=0,N1-1
          do i=0,N0-1
            xx0=dom(0,0)+real(i,f64)*dom(1,0)/real(N0,f64)!-E0(i,j)*dt
            yy0=dom(0,1)+real(j,f64)*dom(1,1)/real(N1,f64)!-E1(i,j)*dt        
            x=0._f64
	    ix=i
	    xx=xx0-E0(i,j)*0.5*dt
	    do ii=1,timecase(1)
	      x=(xx-dom(0,0))/dom(1,0)
              x=x-real(floor(x),f64)
              x=x*real(N0,f64)
              ix=floor(x)
              x=x-real(ix,f64)
              if(ix==N0)then
                ix=0
                x=0._f64
              endif
	      xxn=xx
              xx=xx0-0.5_f64*dt*((1._f64-x)*E0(ix,j)+x*E0(mod(ix+1,N0),j))
            enddo  
	    if(abs(xxn-xx)>eps)then
	      print *,'no convergence x',ii,xx,xxn,xx-xxn,ix,mod(ix+1,N0),x,E0(ix,j),E0(mod(ix+1,N0),j)
	      stop
	    endif

	    do ii=1,timecase(1)
	      y=(yy-dom(0,1))/dom(1,1)
              y=y-real(floor(y),f64)
              y=y*real(N1,f64)
              iy=floor(y)
              y=y-real(iy,f64)
              if(iy==N1)then
                iy=0
                y=0._f64
              endif
	      yyn=yy
              yy=yy0-0.5_f64*dt*((1._f64-y)*((1._f64-x)*E1(ix,iy)+x*E1(mod(ix+1,N0),iy))&
	      +y*((1._f64-x)*E1(ix,mod(iy+1,N1))+x*E1(mod(ix+1,N0),mod(iy+1,N1))))
	      yy=yy-0.5_f64*dt*((1._f64-x)*E1(ix,j)+x*E1(mod(ix+1,N0),j))
            enddo  
	    if(abs(yyn-yy)>eps)then
	      print *,'no convergence y',ii,yy,yyn,yy-yyn,iy,mod(iy+1,N1),y
	      stop
	    endif
            xx=xx-0.5_f64*dt*((1._f64-y)*((1._f64-x)*E0(ix,iy)+x*E0(mod(ix+1,N0),iy))&
	      +y*((1._f64-x)*E0(ix,mod(iy+1,N1))+x*E0(mod(ix+1,N0),mod(iy+1,N1))))  
	    call splperper2d(buf2d,xx,dom(0,0),dom(0,0)+dom(1,0),yy,dom(0,1),dom(0,1)+dom(1,1),f(i,j),N0,N1)
          enddo
        enddo
      case(4)
        do j=0,N1-1
          do i=0,N0-1
            xx0=dom(0,0)+real(i,f64)*dom(1,0)/real(N0,f64)!-E0(i,j)*dt
            yy0=dom(0,1)+real(j,f64)*dom(1,1)/real(N1,f64)!-E1(i,j)*dt        
            x=0._f64
	    ix=i
	    ax=0._f64
	    ay=0._f64
	    do ii=1,timecase(1)
	      xx=xx0-ax
	      yy=yy0-ay
	      x=(xx-dom(0,0))/dom(1,0)
              x=x-real(floor(x),f64)
              x=x*real(N0,f64)
              ix=floor(x)
              x=x-real(ix,f64)
              if(ix==N0)then
                ix=0
                x=0._f64
              endif
	      y=(yy-dom(0,1))/dom(1,1)
              y=y-real(floor(y),f64)
              y=y*real(N1,f64)
              iy=floor(y)
              y=y-real(iy,f64)
              if(iy==N1)then
                iy=0
                y=0._f64
              endif
	      ax=0.5_f64*dt*((1._f64-y)*((1._f64-x)*E0(ix,iy)+x*E0(mod(ix+1,N0),iy))&
	      +y*((1._f64-x)*E0(ix,mod(iy+1,N1))+x*E0(mod(ix+1,N0),mod(iy+1,N1))))
	      ay=0.5_f64*dt*((1._f64-y)*((1._f64-x)*E1(ix,iy)+x*E1(mod(ix+1,N0),iy))&
	      +y*((1._f64-x)*E1(ix,mod(iy+1,N1))+x*E1(mod(ix+1,N0),mod(iy+1,N1))))
	      xxn=xx
	      yyn=yy
	      xx=xx0-ax
	      yy=yy0-ay
            enddo  
	    if(abs(xxn-xx)+abs(yyn-yy)>eps)then
	      print *,ii,xx,xxn,xx-xxn,yy,yyn,yy-yyn
	      stop
	    endif
	    xx=xx0-2._f64*ax
	    yy=yy0-2._f64*ay	    
	    call splperper2d(buf2d,xx,dom(0,0),dom(0,0)+dom(1,0),yy,dom(0,1),dom(0,1)+dom(1,1),f(i,j),N0,N1)
          enddo
        enddo
        		
    end select  	   
  end subroutine advect2d



  subroutine splperper2d(f,xx,xmin,xmax,yy,ymin,ymax,fval,Nx,Ny)
    sll_int32,intent(in)::Nx,Ny
    sll_real64,intent(in)::xx,xmin,xmax,yy,ymin,ymax
    sll_real64,intent(out)::fval
    sll_real64,dimension(0:Nx-1,0:Ny-1)::f
    !f2py intent(in)::f
    sll_int32::i,j
    sll_int32::ix(0:3),iy(0:3)
    sll_real64::x,y
    sll_real64::wx(0:3),wy(0:3),tmp(0:3) 
    x=(xx-xmin)/(xmax-xmin)
    x=x-real(floor(x),f64)
    x=x*real(Nx,f64)
    i=floor(x)
    x=x-real(i,f64)
    if(i==Nx)then
      i=0
      x=0._f64
    endif


    y=(yy-ymin)/(ymax-ymin)
    y=y-real(floor(y),f64)
    y=y*real(Ny,f64)
    j=floor(y)
    y=y-real(j,f64)
   
    if(j==Ny)then
      j=0
      y=0._f64
    endif

   
    wx(0)=(1._f64/6._f64)*(1._f64-x)*(1._f64-x)*(1._f64-x);
    wx(1)=1._f64/6._f64+0.5_f64*(1._f64-x)*(-(1._f64-x)*&
    	 (1._f64-x)+(1._f64-x)+1._f64);
    wx(2)=1._f64/6._f64+0.5_f64*x*(-x*x+x+1._f64);
    wx(3)=(1._f64/6._f64)*x*x*x;

    wy(0)=(1._f64/6._f64)*(1._f64-y)*(1._f64-y)*(1._f64-y);
    wy(1)=1._f64/6._f64+0.5_f64*(1._f64-y)*(-(1._f64-y)*&
    	 (1._f64-y)+(1._f64-y)+1._f64);
    wy(2)=1._f64/6._f64+0.5_f64*y*(-y*y+y+1._f64);
    wy(3)=(1._f64/6._f64)*y*y*y;

    ix(0)=mod(i+Nx-1,Nx)
    ix(1)=i
    ix(2)=mod(i+1,Nx)
    ix(3)=mod(i+2,Nx)

    iy(0)=mod(j+Ny-1,Ny)
    iy(1)=j
    iy(2)=mod(j+1,Ny)
    iy(3)=mod(j+2,Ny)

    
    tmp(0)=wx(0)*f(ix(0),iy(0))+wx(1)*f(ix(1),iy(0))&
    +wx(2)*f(ix(2),iy(0))+wx(3)*f(ix(3),iy(0))
    tmp(1)=wx(0)*f(ix(0),iy(1))+wx(1)*f(ix(1),iy(1))&
    +wx(2)*f(ix(2),iy(1))+wx(3)*f(ix(3),iy(1))
    tmp(2)=wx(0)*f(ix(0),iy(2))+wx(1)*f(ix(1),iy(2))&
    +wx(2)*f(ix(2),iy(2))+wx(3)*f(ix(3),iy(2))
    tmp(3)=wx(0)*f(ix(0),iy(3))+wx(1)*f(ix(1),iy(3))&
    +wx(2)*f(ix(2),iy(3))+wx(3)*f(ix(3),iy(3))
    
    fval=wy(0)*tmp(0)+wy(1)*tmp(1)+wy(2)*tmp(2)+wy(3)*tmp(3)
    !print *,fval,' t',f    
  end subroutine splperper2d


  subroutine thdiagcgper(ftab,Ex,Ey,dom,Nx,Ny,nbdiag,thf,buf,Nbuf,bufc,Nbufc1,Nbufc2,modx,mody)
    sll_int32,intent(in)::Nbuf,Nbufc1,Nbufc2,Nx,Ny,nbdiag,modx,mody
    sll_real64,dimension(Nx,Ny),intent(in)::ftab,Ex,Ey
    sll_real64,dimension(Nbuf)::buf
    complex(f64),dimension(2,Nbufc1,Nbufc2)::bufc
    !f2py intent(in)::bufc,buf
    sll_real64,dimension(0:1,0:1),intent(in)::dom
    sll_real64,dimension(nbdiag),intent(out)::thf
    sll_real64::dx,dy,tmp,l1,l2,enstrophy,mass,maxi
    sll_int32::i,j,i1,j1,im1,jm1
    dx=(dom(1,0)-dom(0,0))/real(Nx,f64);dy=(dom(1,1)-dom(0,1))/real(Ny,f64)
    
    l1=0._f64;l2=0._f64;enstrophy=0._f64;mass=0._f64
    do j=1,Ny
      do i=1,Nx
	l1=l1+abs(ftab(i,j))
    	l2=l2+abs(ftab(i,j))*abs(ftab(i,j))
    	mass=mass+ftab(i,j)
      enddo
    enddo  
    
    do j=1,Ny
      do i=1,Nx
        enstrophy=enstrophy+abs(Ex(i,j))*abs(Ex(i,j))+abs(Ey(i,j))*abs(Ey(i,j))
      enddo	
    enddo	
    l1=l1*dx*dy;l2=l2*dx*dy;mass=mass*dx*dy;enstrophy=enstrophy*dx*dy


    thf(1)=mass;thf(2)=l1;thf(3)=l2;thf(4)=enstrophy;
    
    


    do i=1,Ny
      call dfftf(Nx,ftab(:,i),buf(1:2*Nx+15))
    end do
    do j=1,Ny
      bufc(1,j,1)=cmplx(ftab(1,j),0._f64)
    enddo  
    do i=2,Nx/2
       do j=1,Ny
          bufc(1,j,i)=cmplx(ftab(2*i-2,j),ftab(2*i-1,j))
       end do
    end do
    do j=1,Ny
      bufc(1,j,Nx/2+1)=cmplx(ftab(2*(Nx/2),j),0._f64)
    enddo
 
    do i=1,Nx/2+1
      call zfftf(Ny,bufc(1,1:Ny,i),buf(2*Nx+16:Nbuf))
    enddo  

    maxi=0._f64
    do j=1,Ny
      do i=1,Nx
        tmp=0._f64
	i1=i+1;im1=i-1;
        if(i==Nx)i1=1;if(i==1)im1=Nx
        j1=j+1;jm1=j-1;
        if(j==Ny)j1=1;if(j==1)jm1=Ny
        tmp=tmp+2._f64/3._f64*(Ey(i1,j)-Ey(im1,j))/dom(1,0)
        tmp=tmp+1._f64/6._f64*(Ey(i1,j1)-Ey(im1,j1))/dom(1,0)
	tmp=tmp+1._f64/6._f64*(Ey(i1,jm1)-Ey(im1,jm1))/dom(1,0)
        tmp=tmp-2._f64/3._f64*(Ex(i,j1)-Ex(i,jm1))/dom(1,1)
        tmp=tmp-1._f64/6._f64*(Ex(i1,j1)-Ex(i1,jm1))/dom(1,1)
	tmp=abs(tmp-1._f64/6._f64*(Ex(im1,j1)-Ex(im1,jm1))/dom(1,1))
	if(tmp>maxi)maxi=tmp
      enddo
    enddo

    
    
    thf(5)=sqrt(real(bufc(1,mody+1,modx+1))**2+aimag(bufc(1,mody+1,modx+1))**2)/real(Nx*Ny,f64) !mod=(1,Ny/4-1)
    thf(6)=maxi

    
    
  end subroutine thdiagcgper

  subroutine print2dper(dom,ftab,Nx,Ny,visucase,step,filename)
    sll_int32,intent(in)::Nx,Ny,visucase,step
    sll_real64,dimension(0:1,0:1),intent(in)::dom
    sll_real64,dimension(0:Nx-1,0:Ny-1),intent(in)::ftab
    character(len=*),intent(in)::filename
    character*80,str,str2
    
        
    if(visucase==0)then
      !gnuplot
      call printgp2dper(dom,ftab,Nx,Ny,step,filename)
      !write(str2,*) step
      !write(str,*) 'mv f.dat f'//trim(adjustl((str2)))//'.dat';call system(str)
      !write(str,*) 'mv f.dat '//filename//trim(adjustl((str2)))//'.dat';!call system(str)
    endif
    if(visucase==1)then
      !vtk
      call printvtk2dper(dom,ftab,Nx,Ny,step,filename)
      !write(str2,*) step
      !write(str,*) 'mv f.vtk '//filename//trim(adjustl((str2)))//'.vtk';!call system(str)
    endif
  end subroutine print2dper

  subroutine printgp2dper(dom,ftab,Nx,Ny,step,filename)
    sll_int32,intent(in)::Nx,Ny,step
    sll_real64,dimension(0:1,0:1),intent(in)::dom
    sll_real64,dimension(0:Nx-1,0:Ny-1),intent(in)::ftab
    sll_int32::i,j
    sll_real64::z(0:1),dz(0:1)
    character(len=*),intent(in)::filename
    character*80,str,str2
    write(str2,*) step
    !write(str,*) 'mv f.dat f'//trim(adjustl((str2)))//'.dat';call system(str)
    !write(str,*) 'f'//trim(adjustl((filename)))//trim(adjustl((str2)))//'.dat';!call system(str)
    str=trim(adjustl((filename)))//trim(adjustl((str2)))//'.dat';!call system(str)

    dz(0)=(dom(1,0)-dom(0,0))/real(Nx,f64);dz(1)=(dom(1,1)-dom(0,1))/real(Ny,f64)
    open(unit=900,file=str)
    do j=0,Ny-1
      do i=0,Nx-1
        z(0)=dom(0,0)+real(i,f64)*dz(0)
        z(1)=dom(0,1)+real(j,f64)*dz(1)
        write(900,*) z(0),z(1),ftab(i,j)
      enddo
      i=Nx
      z(0)=dom(0,0)+real(i,f64)*dz(0)
      z(1)=dom(0,1)+real(j,f64)*dz(1)
      write(900,*) z(0),z(1),ftab(0,j)      
      write(900,*) ''      
    enddo
    j=Ny
    do i=0,Nx-1
      z(0)=dom(0,0)+real(i,f64)*dz(0)
      z(1)=dom(0,1)+real(j,f64)*dz(1)
      write(900,*) z(0),z(1),ftab(i,0)
    enddo
    i=Nx
    z(0)=dom(0,0)+real(i,f64)*dz(0)
    z(1)=dom(0,1)+real(j,f64)*dz(1)
    write(900,*) z(0),z(1),ftab(0,0)	  
    write(900,*) ''	       
    close(900)  
  end subroutine printgp2dper

  subroutine printvtk2dper(dom,ftab,Nx,Ny,step,filename)
    sll_int32,intent(in)::Nx,Ny
    sll_real64,dimension(0:1,0:1),intent(in)::dom
    sll_real64,dimension(0:Nx-1,0:Ny-1),intent(in)::ftab
    sll_int32::i,j
    sll_int32,intent(in):: step
    sll_real64::z(0:1),dz(0:1)
    character(len=*),intent(in)::filename
    character*80,str,str2
    write(str2,*) step
    !write(str,*) 'mv f.dat f'//trim(adjustl((str2)))//'.dat';call system(str)
    write(str,*) 'f'//trim(adjustl((filename)))//trim(adjustl((str2)))//'.vtk';!call system(str)
    str=trim(adjustl((filename)))//trim(adjustl((str2)))//'.vtk';!call system(str)
    dz(0)=(dom(1,0)-dom(0,0))/real(Nx,f64);dz(1)=(dom(1,1)-dom(0,1))/real(Ny,f64)
    !open(unit=900,file='f.vtk')
    open(unit=900,file=str,form='formatted')
    write(900,'(A)')                  '# vtk DataFile Version 2.0'
    write(900,'(A)')                  'Exemple'
    write(900,'(A)')                  'ASCII'
    write(900,'(A)')                  'DATASET STRUCTURED_POINTS'
    write(900,'(A,I0,A,I0,A,I0)') 'DIMENSIONS ', Nx+1,' ', Ny+1,' ', 1
    write(900,'(A,I0,A,I0,A,I0)') 'ORIGIN ', 0,' ' , 0,' ' , 0
    !write(900,'(A,F10.4,A,F10.4,A,F10.4)') 'SPACING ', dz(0),' ', dz(1),' ', 1. 
    write(900,*) 'SPACING ', dz(0),' ', dz(1),' ', 1. 
    write(900,*)
    write(900,'(A,I0)')           'POINT_DATA ',(Nx+1)*(Ny+1)
    write(900,'(A,I0)')           'SCALARS f float ',1
    write(900,'(A)')                  'LOOKUP_TABLE default'
    
    do j=0,Ny-1
      do i=0,Nx-1
        z(0)=dom(0,0)+real(i,f64)*dz(0)
        z(1)=dom(0,1)+real(j,f64)*dz(1)
        !write(900,'(F0.8)') ftab(i,j)
        write(900,*) ftab(i,j)
      enddo
      i=Nx
      z(0)=dom(0,0)+real(i,f64)*dz(0)
      z(1)=dom(0,1)+real(j,f64)*dz(1)
      !write(900,'(F0.8)') ftab(0,j)            
      write(900,*) ftab(0,j)            
    enddo
    j=Ny
    do i=0,Nx-1
      z(0)=dom(0,0)+real(i,f64)*dz(0)
      z(1)=dom(0,1)+real(j,f64)*dz(1)
      !write(900,'(F0.8)') ftab(i,0)
      write(900,*) ftab(i,0)
    enddo
    i=Nx
    z(0)=dom(0,0)+real(i,f64)*dz(0)
    z(1)=dom(0,1)+real(j,f64)*dz(1)
    !write(900,'(F0.8)') ftab(0,0)	  	       
    write(900,*) ftab(0,0)	  	       
    close(900)  
  end subroutine printvtk2dper






end module cg_csl_uniform_module




