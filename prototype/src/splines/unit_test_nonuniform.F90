program nonuniform_spline_tester
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
!contact: mehrenbe@math.unistra.fr for this unit_test program

  use cubic_nonuniform_splines
  use numeric_constants
  implicit none



  sll_int32 :: err
  sll_int32 :: N,i,n_array,j1
  sll_real64,dimension(:), pointer :: node_positions,f_per,f_hrmt,a_in,a_out_per,a_out_hrmt,sol_per,sol_hrmt
  sll_real64,dimension(:), pointer :: node_positions_tmp,fe,Xstar,fold
  type(cubic_nonunif_spline_1D), pointer :: spl_per, spl_hrmt  
  sll_real64 :: x,val,sl,sr,xmin,xmax,dx,shift,dt,velocity,M,tmp,linf_err,nb_period
  sll_int32 :: testcase,meshcase,interpcase,nbstep,step
  
  N=16
  n_array = 100
  
  
  !parameters
  N=3200
  testcase=1
  meshcase=1
  interpcase=2
  dt=0.01_f64
  
  velocity=1e0_f64
  nb_period=1.0_f64

  
  
  !definition of the mesh
  SLL_ALLOCATE(node_positions(N+1), err)
  do i=1,N+1
    !node_positions(i)=real(i-1)/real(N)+1.0_f64/real(N)*sin(2.0_f64*sll_pi*real(i-1)/real(N))
    node_positions(i)=2.0_f64*real(i-1,f64)/real(N,f64)!sin(2.0_f64*3.1415_f64*real(i-1)/real(N))
  enddo
  node_positions(1)=5.1_f64
  do i=2,N+1
    node_positions(i)=node_positions(i-1)+1.0_f64/real(N,f64)*(1.0_f64+1.0_f64/real(N*N,f64)+sin(2.0_f64*sll_pi*real(i-1,f64)/real(N,f64)))
  enddo
  
  !print *,node_positions
  !stop
  !definition of the function
  SLL_ALLOCATE(f_hrmt(N+1), err)  
  f_hrmt = 0.0_f64
  !f_hrmt(1) = 0.1_f64
  !f_hrmt(N+1) = 0.2_f64
  do i=1,N+1
    x=node_positions(i)
    !f_hrmt(i)=sin(2.0_f64*sll_pi*x/(node_positions(N+1)-node_positions(1)))
    f_hrmt(i)=exp(x)
  enddo
  sl=2.0_f64*sll_pi
  sr=2.0_f64*sll_pi
  sl=exp(node_positions(1))
  sr=exp(node_positions(N+1))
  !sl=0.0_f64
  !sr=0.0_f64
  
  SLL_ALLOCATE(f_per(N+1), err)  
  f_per = 1.0_f64
  f_per(1) = 1.0_f64
  do i=1,N+1
    x=node_positions(i)
    !f_per(i)=sin(2.0_f64*sll_pi*x)
    f_per(i)=sin(2.0_f64*sll_pi*x/(node_positions(N+1)-node_positions(1)))
  enddo

  !f_per = f_hrmt

  print *, '#proceed to allocate the spline...'  
  spl_per =>  new_cubic_nonunif_spline_1D( N, PERIODIC_SPLINE)
  spl_hrmt =>  new_cubic_nonunif_spline_1D( N, HERMITE_SPLINE)
  call compute_spline_nonunif( f_per, spl_per, node_positions)
  call compute_spline_nonunif( f_hrmt, spl_hrmt, node_positions,sl,sr)

  
  !interpolating points
  !allocation for a_in and a_out
  SLL_ALLOCATE(a_in(n_array), err)  
  SLL_ALLOCATE(a_out_per(n_array), err)  
  SLL_ALLOCATE(a_out_hrmt(n_array), err)  
  SLL_ALLOCATE(sol_per(n_array), err)  
  SLL_ALLOCATE(sol_hrmt(n_array), err)  


  do i=1, n_array
    a_in(i) = spl_per%xmin+real(i-1,f64)/real(n_array-1,f64)*(spl_per%xmax-spl_per%xmin) 
  enddo
  a_in(1)=0.001_f64
  a_in(2)=1.0_f64
  a_in(3)=0.9_f64
  a_in(4)=0.01_f64
  a_in(5)=0.0_f64
  a_in(6)=0.1_f64
  a_in(7)=0.05_f64
  a_in(8)=0.025_f64
  a_in(9)=0.3_f64
  a_in(10)=0.2_f64
  do i=1,10
    a_in(i) = spl_per%xmin+a_in(i)*(spl_per%xmax-spl_per%xmin)
  enddo

  do i=1,n_array
    x=a_in(i)
    sol_per(i)=sin(2.0_f64*sll_pi*x/(node_positions(N+1)-node_positions(1)))
  enddo

  do i=1,n_array
    x=a_in(i)
    sol_hrmt(i)=exp(x)
  enddo

  !sol_per=sol_hrmt
  !call compute_spline_nonunif( f, spl_per)
  
  x = spl_per%xmin!1.0_f64!+1.0e-12_f64
  x = x+0.1234589_f64*(spl_per%xmax-spl_per%xmin)
  !i=N+1
  !x = node_positions(i)
  val = interpolate_value_nonunif( x, spl_per )
  print *,'#',val,x,sin(2.0_f64*sll_pi*x)-val
  val = interpolate_value_nonunif( x, spl_hrmt )
  print *,'#',val,x,sin(2.0_f64*sll_pi*x)-val
  

  call interpolate_array_value_nonunif( a_in, a_out_per,n_array, spl_per)
  call interpolate_array_value_nonunif( a_in, a_out_hrmt,n_array, spl_hrmt)
  do i=1, n_array
    !print *, a_in(i),a_out_per(i),a_out_hrmt(i),sol_per(i),sol_hrmt(i)
  enddo
  
  
  
  
  print *, '#proceed to delete the spline...'  
  call delete_cubic_nonunif_spline_1D( spl_per, err)
  call delete_cubic_nonunif_spline_1D( spl_hrmt, err)

  ! computation of conservative constant advection on non uniform mesh
  
  
  !

  SLL_ALLOCATE(node_positions_tmp(N+1), err)
  SLL_ALLOCATE(fe(N+1), err)
  SLL_ALLOCATE(Xstar(N+1), err)
  SLL_ALLOCATE(fold(N+1), err)


  xmin=0._f64
  xmax=1._f64
  dx=(xmax-xmin)/real(N,f64)
  
  !initialization of X
  node_positions(1)=0._f64
  do i=2,N+1
    node_positions(i)=xmin+real(i-1,f64)*dx
    !Xmesh(i)=xmin+(real(i,rk)+0.5*sin(2._rk*M_PI*real(i,rk)*dx))*dx
    !Xmesh(i)=Xmesh(i-1)+(2._rk+sin(2._rk*M_PI*real(i,rk)/real(N+1,rk)))*3._rk
    node_positions(i)=node_positions(i-1)+(real(i-1,f64)*real(N+1-(i-1),f64))
    !Xmesh(i)=(real(mod(i+N/2,N),rk)/real(N,rk))**2
    !if(i>0.and.Xmesh(i)<Xmesh(i-1))then
    !  Xmesh(i)=Xmesh(i)+1._rk
    !endif
  enddo
  node_positions=node_positions/node_positions(N+1)
  j1=N/3;shift=0._f64
  do i=0,N    
    node_positions_tmp(i+1)=node_positions(j1+1)+shift
    j1=j1+1
    if(j1>=N)then
      j1=j1-N;shift=shift+1._f64
    endif
  enddo
  node_positions=node_positions_tmp-node_positions_tmp(1)

  xmin=node_positions(1)
  xmax=node_positions(N+1)
  dx=(xmax-xmin)/real(N,f64)
  nbstep=int(nb_period*(xmax-xmin)/(dabs(velocity)*dt))
  nbstep=1
  !initialization of f
  do i=0,N-1
    !f(i)=0.7_rk+sin(2._rk*M_PI*(Xmesh(i)+0.3))
    f_per(i+1)=(node_positions(i+1)+node_positions(i+2))*0.5_f64!
    f_per(i+1)=0.7_f64+sin(2._f64*sll_pi*((node_positions(i+1)+node_positions(i+2))*0.5_f64+0.3_f64))
    !f(i)=sin(M_PI*(Xmesh(i)+Xmesh(i+1)))
    fe(i+1)=0.7_f64+sin(2._f64*sll_pi*((node_positions(i+1)+node_positions(i+2))*0.5_f64+0.3_f64-real(nbstep,f64)*velocity*dt))
  enddo
  
  fe(N+1)=fe(1)
  !from f compute the mean
  do i=0,N-1
    f_per(i+1)=f_per(i+1)*(node_positions(i+2)-node_positions(i+1))/dx
  enddo
  
  !initializations
  
  spl_per =>  new_cubic_nonunif_spline_1D( N, PERIODIC_SPLINE)
  call compute_spline_nonunif( f_per, spl_per, node_positions)

  
  !algorithm
  
  !definition of Xstar
  node_positions=(node_positions-xmin)/(xmax-xmin)
  Xstar=node_positions-velocity*dt/(xmax-xmin)

  !change Xstar so that it remains in (0,1)
  
  do i=1,N+1
    do while (Xstar(i).gt.1._f64)
      Xstar(i) = Xstar(i)-1._f64
    end do
    do while (Xstar(i).lt.0._f64)
      Xstar(i) = Xstar(i)+1._f64
    end do    
  enddo

  do step=1,nbstep
    fold=f_per
    !we compute the splines coefficients by solving the LU decomposition
    M=0._f64
    do i=1,N
      M=M+fold(i)
    enddo
    !M=M/real(N,rk)
    do i=1,N
      fold(i)=fold(i)-M*(node_positions(i+1)-node_positions(i))!/dx
    enddo    
    f_per(1)=0._f64
    do i=2,N
      f_per(i)=f_per(i-1)+fold(i-1)
    enddo
    fold=f_per
    
    !call of compute_spline and interpolations
    call compute_spline_nonunif( fold, spl_per)
    
    
    call interpolate_array_value_nonunif( Xstar, f_per, N, spl_per)
    
    
    tmp=f_per(1)!;for(i=0;i<Nx-1;i++)p[i]=p[i+1]-p[i];p[Nx-1]=tmp+M-p[Nx-1];
    do i=1,N-1
      f_per(i)=f_per(i+1)-f_per(i)+M*(node_positions(i+1)-node_positions(i))
    enddo
    f_per(N)=tmp-f_per(N)+M*(node_positions(1)+1._f64-node_positions(N))
  enddo

  !from mean compute f
  do i=1,N
    f_per(i)=f_per(i)*dx/(node_positions(i+1)-node_positions(i))
  enddo
  
  linf_err=0.0_f64
  f_per(N+1) = f_per(1)
  do i=1,N+1
    if(linf_err.lt.(dabs(f_per(i)-fe(i)))) then
      linf_err = dabs(f_per(i)-fe(i))
    endif
    !print *,node_positions(i), f_per(i), fe(i), f_per(i)-fe(i)
  enddo
  print *,"#err=",linf_err,nbstep
  print *,"#M=",M
end program nonuniform_spline_tester
