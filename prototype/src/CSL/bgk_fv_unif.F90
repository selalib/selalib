program bgk_fv_unif
#include "sll_working_precision.h"
#include "sll_mesh_types.h"
#include "sll_memory.h"
  use sll_constants
  use distribution_function
  use sll_diagnostics
  use sll_csl
  use sll_splines
  implicit none
  external compute_translate_nodes_periodic,compute_non_unif_integral
  external Compute_flux
  !external poisson1dpertrap
  sll_real64 :: x2_min,x2_max,x1_min,x1_max,x1,x2,delta_x1,delta_x2
  sll_real64 :: mu,xi,L,H
  !sll_real64 :: Flux_x1_p, Flux_x1_m,Flux_x2_p,Flux_x2_m
  sll_real64 :: Flux_p05,Flux_m05
  sll_real64 :: b1,b2,b3,b4, a21, a32,a42,a43,a41,a31
  sll_int  ::  rk
  sll_int    :: i,j,N_phi,err,N_x1,N_x2,i1,i2,N,nb_step
  LOGICAL :: ex
    !type(sll_distribution_function_2D_t), pointer :: dist_func
  sll_real64,dimension(:), pointer :: phi,node_positions_x1,node_positions_x2
  sll_real64,dimension(:), pointer :: abar_x1,abar_x2
  sll_real64,dimension(:), pointer :: new_node_positions,buf_1d,rho,E,E_store, rho_exact
  sll_real64,dimension(:,:), pointer :: f,f_store, chi,sigma, Flux,a1,a2,f_tmp,K1,K2,K3,K4
  sll_real64 :: phi_val,delta_x1_phi,xx,dt,alpha,val,beta
  sll_int :: i1m2,i1m1,i2m2,i2m1,i1p1
  sll_int :: ii,step,test_case,order, visu_step
  type(cubic_nonunif_spline_1D), pointer :: spl_per_x1, spl_per_x2
  
  

  N_x1 =128!,0!256! 128!128
  N_x2 =128!128!256! 128!256
  rk=2
  dt = 0.01_f64
  nb_step = 6000
  order=3
  test_case = 4
  visu_step = 10
  
  N = max(N_x1,N_x2)
  
  print *,'#max(N_x1,N_x2)=',N
  
  
  x2_max = 10._f64
  x2_min = -x2_max
  
  SLL_ALLOCATE(f(N_x1+1,N_x2+1),err)
  SLL_ALLOCATE(f_tmp(N_x1+1,N_x2+1),err)
  SLL_ALLOCATE(K1(N_x1+1,N_x2+1),err)
  SLL_ALLOCATE(K2(N_x1+1,N_x2+1),err)
  SLL_ALLOCATE(K3(N_x1+1,N_x2+1),err)
  SLL_ALLOCATE(K4(N_x1+1,N_x2+1),err)
  SLL_ALLOCATE(f_store(N_x1+1,N_x2+1),err)
  SLL_ALLOCATE(node_positions_x1(N_x1+1),err)
  SLL_ALLOCATE(node_positions_x2(N_x2+1),err)
  SLL_ALLOCATE(new_node_positions(N+1),err)
  SLL_ALLOCATE(buf_1d(N+1),err)
  SLL_ALLOCATE(rho(N_x1+1),err)
  SLL_ALLOCATE(rho_exact(N_x1+1),err)
  SLL_ALLOCATE(E(N_x1+1),err)
  SLL_ALLOCATE(E_store(N_x1+1),err)
  !SLL_ALLOCATE(chi(N_x1+1+2,N_x2+1+2),err)
  !SLL_ALLOCATE(sigma(N_x1+1+2,N_x2+1+2),err)
  SLL_ALLOCATE(a1(N_x1+1,N_x2+1),err)
  SLL_ALLOCATE(a2(N_x1+1,N_x2+1),err)
  SLL_ALLOCATE(Flux(N_x1+1,N_x2+1),err)
  !SLL_ALLOCATE(Flux_x2(N_x1+1,N_x2+1),err)
  !SLL_ALLOCATE(abar_x1(N_x1+1),err)
  !SLL_ALLOCATE(abar_x2(N_x2+1),err)
  spl_per_x1 =>  new_cubic_nonunif_spline_1D( N_x1, PERIODIC_SPLINE)
  spl_per_x2 =>  new_cubic_nonunif_spline_1D( N_x2, PERIODIC_SPLINE)
  
  
  !physical parametersreturn
  mu=0.92_f64
  xi=0.90_f64
  L=14.71_f64
  
  inquire(file='half_phi.dat', exist=ex) 
  if(.not.(ex))then
    print *,'file half_phi.dat does not exist'
    stop
  endif  
  open(unit=900,file='half_phi.dat')  
    read(900,*) N_phi,L
    N_phi = 2*N_phi
    SLL_ALLOCATE(phi(N_phi+1),err)
    do j=1,N_phi/2+1
      read(900,*) i,x1,x2
      phi(i)=x1
    enddo
    do j=N_phi/2+2,N_phi+1
      phi(j)=phi(N_phi+2-j)
    enddo
  close(900)
  
  !L = 4._f64*sll_pi
  
  x1_min = 0._f64
  x1_max = L

  if(test_case==4.or.test_case==5)then
    L = 4._f64*sll_pi
    x1_min = 0._f64
    x1_max = L
    x2_min = -6._f64
    x2_max = -x2_min
  endif



  delta_x1_phi = (x1_max-x1_min)/real(N_phi,f64)
  
  

  open(unit=900,file='phi0.dat')  
    do i1=1,N_phi+1
      x1 = x1_min+real(i1-1,f64)*delta_x1_phi
      write(900,*) x1,phi(i1)
    enddo
  close(900)
  
  

  
  print *,'#N_phi=',N_phi
  
  delta_x1 = (x1_max-x1_min)/real(N_x1,f64)
  delta_x2 = (x2_max-x2_min)/real(N_x2,f64)
   alpha=(1._f64/delta_x1)
   beta =(1._f64/delta_x2)

  do i1=1,N_x1+1
    x1 = x1_min+real(i1-1,f64)*delta_x1
    node_positions_x1(i1) = x1
  enddo
  do i2=1,N_x2+1
    x2 = x2_min+real(i2-1,f64)*delta_x2
    node_positions_x2(i2) = x2
  enddo
  
  
  
  do i1=1,N_x1+1
    do i2=1,N_x2+1
      x1 = x1_min+real(i1-1,f64)*delta_x1
      x2 = x2_min+real(i2-1,f64)*delta_x2
      phi_val = 0._f64
      xx = (x1-x1_min)/(x1_max-x1_min)
      if(xx<=0._f64)then
        xx = 0._f64
      endif
      if(xx>=1._f64)then
        xx = 1._f64-1e-15_f64
      endif
      xx = xx*real(N_phi,f64)
      ii = floor(xx)
      xx = xx-real(ii,f64)      
      phi_val = (1._f64-xx)*phi(ii+1)+xx*phi(ii+2)
      !phi_val = 0._f64
      H = 0.5_f64*x2*x2 + phi_val

      if(test_case==1)then
        val = mu/(sqrt(2._f64*sll_pi))*(2._f64-2._f64*xi)/(3._f64-2._f64*xi)*(1._f64+H/(1._f64-xi))*exp(-H)
      endif
      if(test_case==2)then
        val = mu/(sqrt(2._f64*sll_pi))*(2._f64-2._f64*xi)/(3._f64-2._f64*xi)*(1._f64+H/(1._f64-xi))*exp(-H)
        val = val*(1._f64+0.1_f64*cos(2._f64*sll_pi/L*x1))
      endif
      if(test_case==3)then
        val = exp(-0.5_f64*40._f64*((x1-.5_f64)**2+(x2-.5_f64)**2))
      endif
      if(test_case==4)then
        !linear landau damping
        val = 1._f64/(sqrt(2._f64*sll_pi))*exp(-0.5_f64*x2*x2)
        val = val*(1._f64+0.001_f64*cos(2._f64*sll_pi/L*x1))
        
      endif
      if(test_case==5)then
        !non linear landau damping
        val = 1._f64/(sqrt(2._f64*sll_pi))*exp(-0.5_f64*x2*x2)
        val = val*(1._f64+0.5_f64*cos(2._f64*sll_pi/L*x1))
      endif
      f(i1,i2) = val
    enddo
  enddo
      
   do i1=1,N_x1+1
      xx = (real(i1,f64)-0.5_f64)/real(N_x1,f64)
      if(xx<=0._f64)then
        xx = 0._f64
      endif
      if(xx>=1._f64)then
        xx = xx-1._f64!1._f64-1e-15_f64
      endif
      if(xx<=0._f64)then
        xx = 0._f64
      endif      
      xx = xx*real(N_phi,f64)
      ii = floor(xx)
      xx = xx-real(ii,f64)      
      phi_val = (1._f64-xx)*phi(ii+1)+xx*phi(ii+2)
      !rho(i1) =phi_val!****
      rho_exact(i1) = mu*(3._f64-2._f64*xi+2*phi_val)/(3._f64-2._f64*xi)*exp(-phi_val)
      
  enddo
   !Compute E0 and rho0 

    f_store = f
    !compute E
    do i1=1,N_x1+1
      buf_1d(1:N_x2+1) = f_store(i1,1:N_x2+1)
      call compute_non_unif_integral(node_positions_x2(1:N_x2+1),buf_1d(1:N_x2+1),N_x2+1,val)
      rho(i1)=val-1._f64
    enddo
    E=rho
    call poisson1dpertrap(E,x1_max-x1_min,N_x1)
    E_store = E
    
     !diagnostic
      do i1=1,N_x1+1
        x1 = x1_min+real(i1-1._f64)*delta_x1
        !write(*,*) x1,E(i1),rho(i1)!,rho_exact(i1)
      enddo
     ! stop

 
  do step=1,nb_step
    
    
    

    !Inisialize a1 and a2
      do i1=1,N_x1+1
        x1 = x1_min+real(i1-1._f64)*delta_x1
        do i2=1,N_x2+1
           x2 = x2_min+real(i2-1._f64)*delta_x2
           a1(i1,i2)=x2
           a2(i1,i2)=E_store(i1)    
          ! chi(i1,i2)   =a1(i1,i2)*f_store(i1,i2)
           !sigma(i1,i2) =a2(i1,i2)*f_store(i1,i2)
       enddo
   enddo
    

  if(rk==1) then
    f_tmp=f
    call Compute_flux(a1,a2,f_tmp,f_store,Flux,N_x1,N_x2,x1_min,x2_min,delta_x1,delta_x2,order)
   do i1=1,N_x1+1
     do i2=1,N_x2+1
        f(i1,i2)=f_store(i1,i2)+dt*Flux(i1,i2)
     enddo
   enddo
  endif
 
    if(rk==2) then
   
        !Flux=(-alpha*Flux_x1(i1,i2)-beta*Flux_x2(i1,i2))
        b1=0._f64
        b2=1._f64
        a21=1._f64/2._f64
        
        f_tmp=f
        call Compute_flux(a1,a2,f_tmp,f_store,Flux,N_x1,N_x2,x1_min,x2_min,delta_x1,delta_x2,order)
         K1=Flux
        do i1=1,N_x1+1
          do i2=1,N_x2+1
            f_tmp(i1,i2)=f_store(i1,i2)+K1(i1,i2)*dt*a21
          enddo
         enddo
        call Compute_flux(a1,a2,f_tmp,f_store,Flux,N_x1,N_x2,x1_min,x2_min,delta_x1,delta_x2,order)
         K2=Flux
   do i1=1,N_x1+1
     do i2=1,N_x2+1
        
        !K2(i1,i2)=0._f64/3._f64
        f(i1,i2)=f_store(i1,i2)+dt*(b1*K1(i1,i2)+b2*K2(i1,i2))
      !-alpha*Flux_x1(i1,i2)-beta*Flux_x2(i1,i2)
     enddo
    enddo
  endif
    if(rk==4) then
   
        !Flux=(-alpha*Flux_x1(i1,i2)-beta*Flux_x2(i1,i2))
        b1=1._f64/6._f64
        b2=1._f64/3._f64
        b3=1._f64/3._f64
        b4=1._f64/6._f64
       
        a21=1._f64/2._f64
        a32=1._f64/2._f64
        a31=0._f64
        a41=0._f64
        a42=0._f64
        a43=1._f64
        f_tmp=f_store
        call Compute_flux(a1,a2,f_tmp,f_store,Flux,N_x1,N_x2,x1_min,x2_min,delta_x1,delta_x2,order)
         K1=Flux
        do i1=1,N_x1+1
          do i2=1,N_x2+1
            f_tmp(i1,i2)=f_store(i1,i2)+K1(i1,i2)*dt*a21
          enddo
         enddo
        call Compute_flux(a1,a2,f_tmp,f_store,Flux,N_x1,N_x2,x1_min,x2_min,delta_x1,delta_x2,order)
         K2=Flux

         do i1=1,N_x1+1
          do i2=1,N_x2+1
            f_tmp(i1,i2)=f_store(i1,i2)+K1(i1,i2)*dt*a31+K2(i1,i2)*dt*a32
          enddo
         enddo
        call Compute_flux(a1,a2,f_tmp,f_store,Flux,N_x1,N_x2,x1_min,x2_min,delta_x1,delta_x2,order)
         K3=Flux
         do i1=1,N_x1+1
          do i2=1,N_x2+1
            f_tmp(i1,i2)=f_store(i1,i2)+K1(i1,i2)*dt*a41+K2(i1,i2)*dt*a42+K3(i1,i2)*dt*a43
          enddo
         enddo
        call Compute_flux(a1,a2,f_tmp,f_store,Flux,N_x1,N_x2,x1_min,x2_min,delta_x1,delta_x2,order)
         K4=Flux
   do i1=1,N_x1+1
     do i2=1,N_x2+1
        f(i1,i2)=f_store(i1,i2)+dt*(b1*K1(i1,i2)+b2*K2(i1,i2)+b3*K3(i1,i2)+b4*K4(i1,i2))
      !-alpha*Flux_x1(i1,i2)-beta*Flux_x2(i1,i2)
     enddo
    enddo
  endif
    
    f_store = f
    !compute E
    do i1=1,N_x1+1
      buf_1d(1:N_x2+1) = f_store(i1,1:N_x2+1)
      call compute_non_unif_integral(node_positions_x2(1:N_x2+1),buf_1d(1:N_x2+1),N_x2+1,val)
      rho(i1)=val-1._f64
  
    enddo
    E=rho
    call poisson1dpertrap(E,x1_max-x1_min,N_x1)
    E(N_x1+1)=E(1)
    E_store = E
    
     !diagnostic
    val=0._f64
    do i1=1,N_x1
      val = val+E_store(i1)*E_store(i1)
    enddo
    val = val/real(N_x1,f64)
    print *,step, (real(step,f64)-1._f64)*dt,val
    
     
   if(step==100000) then
   do i1=1,N_x1+1
      x1 = x1_min+real(i1-1._f64)*delta_x1
      write(*,*) x1,E(i1),rho(i1)!,rho_exact(i1)
    enddo
    stop
  endif
    !diagnostic
    ! do i1=1,N_x1
      !print *,x1_min+real(i1-1)*delta_x1,E(i1),rho(i1)
    !enddo
   
 enddo !end time step
    
    
  !enddo !end time step

  open(unit=900,file='field_final_vf_unif.dat')  
    do i1=1,N_x1+1
      x1 = x1_min+real(i1-1._f64)*delta_x1
      write(900,*) x1,E(i1),rho(i1)
    enddo
  close(900)


  
end program


 subroutine Compute_flux(a1,a2,f,f_store,Flux,N_x1,N_x2,x1_min,x2_min,delta_x1,delta_x2,order)
  use sll_constants
  implicit none
  
  
  sll_int,intent(in)::N_x1,N_x2,order
  !sll_real64,dimension(:,:), pointer :: chi,sigma
   sll_real64,dimension(N_x1+1,N_x2+1):: a1,a2,f,f_store, Flux_x1,Flux_x2,flux
  sll_real64,dimension(N_x1+1+2,N_x2+1+2) ::chi,sigma
  sll_real64,dimension(N_x1+1)::abar_x1
  sll_real   ,dimension(N_x1+1)::abar_x2!E
  sll_real64,intent(in)::x1_min,x2_min,delta_x1,delta_x2!L
  sll_int::i1,i2,i1m2,i1m1, i2m2,i2m1,i1p1,i1p2,i2p1,i2p2,i1m3
  sll_real64::x1,x2,Flux_p05,Flux_m05,alpha,beta!,eold,enew,dx2,tmp
    !compute chi(uij) and sigma(uij)
     alpha=(1._f64/delta_x1)
     beta =(1._f64/delta_x2)
   do i1=1,N_x1+1
        do i2=1,N_x2+1   
           chi(i1,i2)   =a1(i1,i2)*f(i1,i2)
           sigma(i1,i2) =a2(i1,i2)*f(i1,i2)
       enddo
   enddo

if(order==3) then
   
    chi(N_x1+2,1:N_x2+1+2)=chi(2,1:N_x2+1+2) 
    chi(N_x1+3,1:N_x2+1+2)=chi(3,1:N_x2+1+2)
    sigma(1:N_x1+1+2,N_x2+2)=sigma(1:N_x1+1+2,2)
    sigma(1:N_x1+1+2,N_x2+3)=sigma(1:N_x1+1+2,3)

  ! For fixed  i2 
 do i2=1,N_x2+1
      x2 = x2_min+real(i2-1,f64)*delta_x2
     do i1=1,N_x1
     abar_x1(i1)=(chi(i1+1,i2)-chi(i1,i2))/(f_store(i1+1,i2)-f_store(i1,i2))!a1(i1,i2)!special uniform
     enddo
     abar_x1(N_x1+1)=abar_x1(1)

  do i1=1,N_x1+1
    
       i1p1=i1+1
       i1p2=i1+2
       i1m2=i1-2
       i1m1=i1-1
     
       if (i1m1 <=0) then 
         i1m1=i1m1+N_x1
       endif
       if (i1m2 <=0) then 
         i1m2=i1m2+N_x1     
       endif
        
       !if (i1p1 >=N_x1+1) then 
        ! i1p1=i1p1+N_x1+1
       !endif
       !if (i1p2 >=N_x1+1) then 
        ! i1p2=i1p2+N_x1+1
       !endif
    
       if(abar_x1(i1)>=0) then   
        Flux_p05=-(1._f64/6._f64)*chi(i1m1,i2)+(5._f64/6._f64)*chi(i1,i2)+(1._f64/3._f64)*chi(i1p1,i2)! stencil i1-1,i1,i1+1  
        Flux_m05=-(1._f64/6._f64)*chi(i1m2,i2)+(5._f64/6._f64)*chi(i1m1,i2)+(1._f64/3._f64)*chi(i1,i2)! stencil i1-1,i1,i1+1  
        Flux_x1(i1,i2)=Flux_p05-Flux_m05
       else 
        Flux_p05=(1._f64/3._f64)*chi(i1,i2)+(5._f64/6._f64)*chi(i1p1,i2)-(1._f64/6._f64)*chi(i1p2,i2)! stencil i1,i1+1,i1+2
        Flux_m05=(1._f64/3._f64)*chi(i1m1,i2)+(5._f64/6._f64)*chi(i1,i2)-(1._f64/6._f64)*chi(i1p1,i2)! stencil i1,i1+1,i1+2
        Flux_x1(i1,i2)=Flux_p05-Flux_m05
      endif  

     
   

  enddo
 enddo 
  !For fixed i1
 do i1=1,N_x1+1
     do i2=1,N_x2
     abar_x2(i2)=(sigma(i1,i2+1)-sigma(i1,i2))/(f_store(i1,i2+1)-f_store(i1,i2))!a1(i1,i2)!(sigma(i1,i2+1)-sigma(i1,i2))/(f_store(i1,i2+1)-f_store(i1,i2))
     enddo
     !abar_x2(N_x2+1)=abar_x2(1)!special uniform
  do i2=1,N_x2+1
     i2p1=i2+1
     i2p2=i2+2
     i2m2=i2-2
     i2m1=i2-1
     
     
     if (i2m1 <=0) then 
         i2m1=i2m1+N_x2
     endif
     if (i2m2 <=0) then 
         i2m2=i2m2+N_x2      
     endif
        
      ! if (i2p1 >=N_x2+1) then 
        ! i2p1=i2p1+N_x2+1
       !endif
       !if (i2p2 >=N_x2+1) then 
       !  i2p2=i2p2+N_x2+1
    !   endif
   
    
     if(abar_x2(i2)>=0) then   
        Flux_p05=-(1._f64/6._f64)*sigma(i1,i2m1)+(5._f64/6._f64)*sigma(i1,i2)  +(1._f64/3._f64)*sigma(i1,i2p1)! stencil i2-1,i2,i2+1  
        Flux_m05=-(1._f64/6._f64)*sigma(i1,i2m2)+(5._f64/6._f64)*sigma(i1,i2m1)+(1._f64/3._f64)*sigma(i1,i2  )! stencil i2-1,i2,i2+1  
        Flux_x2(i1,i2)=Flux_p05-Flux_m05
     else 
        Flux_p05=(1._f64/3._f64)*sigma(i1,i2  )+(5._f64/6._f64)*sigma(i1,i2p1)-(1._f64/6._f64)*sigma(i1,i2p2)! stencil i2,i2+1,i2+2
        Flux_m05=(1._f64/3._f64)*sigma(i1,i2m1)+(5._f64/6._f64)*sigma(i1,i2)  -(1._f64/6._f64)*sigma(i1,i2p1)! stencil i2,i2+1,i2+2
        Flux_x2(i1,i2)=Flux_p05-Flux_m05
     endif
    
  

  enddo
 enddo
  do i1=1,N_x1+1
     do i2=1,N_x2+1
       Flux(i1,i2)=-(alpha*Flux_x1(i1,i2)+beta*Flux_x2(i1,i2))
     enddo
   enddo
endif!order


  if(order==5) then
    !stop
    !sigma(N_x1+2,1:N_x2+1+2)=sigma(2,1:N_x2+1+2) 
    !sigma(N_x1+3,1:N_x2+1+2)=sigma(3,1:N_x2+1+2)
    !sigma(1:N_x1+1+2,N_x2+2)=sigma(1:N_x1+1+2,2)
    !sigma(1:N_x1+1+2,N_x2+3)=sigma(1:N_x1+1+2,3)
    !chi(N_x1+2,1:N_x2+1+2)=chi(2,1:N_x2+1+2)
    !chi(N_x1+3,1:N_x2+1+2)=chi(3,1:N_x2+1+2)
    !chi(1:N_x1+1+2,N_x2+2)=chi(1:N_x1+1+2,2)
    !chi(1:N_x1+1+2,N_x2+3)=chi(1:N_x1+1+2,3)
  
    chi(N_x1+2,1:N_x2+1+4)=chi(2,1:N_x2+1+4) 
    chi(N_x1+3,1:N_x2+1+4)=chi(3,1:N_x2+1+4)
    chi(N_x1+4,1:N_x2+1+4)=chi(4,1:N_x2+1+4)
    !chi(1:N_x1+1+2,N_x2+2)=chi(1:N_x1+1+2,2)
    !chi(1:N_x1+1+2,N_x2+3)=chi(1:N_x1+1+2,3)
    !sigma(N_x1+2,1:N_x2+1+2)=sigma(2,1:N_x2+1+2)
    !sigma(N_x1+3,1:N_x2+1+2)=sigma(3,1:N_x2+1+2)
    sigma(1:N_x1+1+4,N_x2+2)=sigma(1:N_x1+1+4,2)
    sigma(1:N_x1+1+4,N_x2+3)=sigma(1:N_x1+1+4,3)
    sigma(1:N_x1+1+4,N_x2+4)=sigma(1:N_x1+1+4,4)

  ! For fixed  i2 
 do i2=1,N_x2+1
      x2 = x2_min+real(i2-1,f64)*delta_x2
     do i1=1,N_x1
     abar_x1(i1)=a1(i1,i2)!special uniform
     enddo
     !abar_x1(N_x1+1)=abar_x1(1)

  do i1=1,N_x1+1
    
       i1p1=i1+1
       i1p2=i1+2
       i1m3=i1-3
       i1m2=i1-2
       i1m1=i1-1
     
       if (i1m1 <=0) then 
         i1m1=i1m1+N_x1+1
       endif
       if (i1m2 <=0) then 
         i1m2=i1m2+N_x1 +2     
       endif
        if (i1m3 <=0) then 
         i1m2=i1m2+N_x1 +3     
       endif
       !if (i1p1 >=N_x1+1) then 
        ! i1p1=i1p1+N_x1+1
       !endif
       !if (i1p2 >=N_x1+1) then 
        ! i1p2=i1p2+N_x1+1
       !endif
    
       if(abar_x1(i1)>=0) then   
        Flux_p05=-(1._f64/6._f64)*chi(i1m1,i2)+(5._f64/6._f64)*chi(i1,i2)+(1._f64/3._f64)*chi(i1+1,i2)! stencil i1-1,i1,i1+1  
        Flux_m05=-(1._f64/6._f64)*chi(i1m2,i2)+(5._f64/6._f64)*chi(i1m1,i2)+(1._f64/3._f64)*chi(i1,i2)! stencil i1-1,i1,i1+1  
        Flux_x1(i1,i2)=Flux_p05-Flux_m05
       else 
        Flux_p05=(1._f64/3._f64)*chi(i1,i2)+(5._f64/6._f64)*chi(i1p1,i2)-(1._f64/6._f64)*chi(i1p2,i2)! stencil i1,i1+1,i1+2
        Flux_m05=(1._f64/3._f64)*chi(i1m1,i2)+(5._f64/6._f64)*chi(i1,i2)-(1._f64/6._f64)*chi(i1p1,i2)! stencil i1,i1+1,i1+2
        Flux_x1(i1,i2)=Flux_p05-Flux_m05
      endif  

     
   

  enddo
 enddo 
  !For fixed i1
 do i1=1,N_x1+1
     do i2=1,N_x2+1
     abar_x2(i2)=a1(i1,i2)!(sigma(i1,i2+1)-sigma(i1,i2))/(f_store(i1,i2+1)-f_store(i1,i2))
     enddo
     abar_x2(N_x2+1)=abar_x2(1)!special uniform
  do i2=1,N_x2+1
     i2p1=i2+1
     i2p2=i2+2
     i2m2=i2-2
     i2m1=i2-1
     
     
     if (i2m1 <=0) then 
         i2m1=i2m1+N_x2+1
     endif
     if (i2m2 <=0) then 
         i2m2=i2m2+N_x2 +2     
     endif
        
      ! if (i2p1 >=N_x2+1) then 
        ! i2p1=i2p1+N_x2+1
       !endif
       !if (i2p2 >=N_x2+1) then 
       !  i2p2=i2p2+N_x2+1
    !   endif
   
    
     if(abar_x2(i2)>=0) then   
        Flux_p05=-(1._f64/6._f64)*sigma(i1,i2m1)+(5._f64/6._f64)*sigma(i1,i2)  +(1._f64/3._f64)*sigma(i1,i2p1)! stencil i2-1,i2,i2+1  
        Flux_m05=-(1._f64/6._f64)*sigma(i1,i2m2)+(5._f64/6._f64)*sigma(i1,i2m1)+(1._f64/3._f64)*sigma(i1,i2  )! stencil i2-1,i2,i2+1  
        Flux_x2(i1,i2)=Flux_p05-Flux_m05
     else 
        Flux_p05=(1._f64/3._f64)*sigma(i1,i2  )+(5._f64/6._f64)*sigma(i1,i2p1)-(1._f64/6._f64)*sigma(i1,i2p2)! stencil i2,i2+1,i2+2
        Flux_m05=(1._f64/3._f64)*sigma(i1,i2m1)+(5._f64/6._f64)*sigma(i1,i2)  -(1._f64/6._f64)*sigma(i1,i2p1)! stencil i2,i2+1,i2+2
        Flux_x2(i1,i2)=Flux_p05-Flux_m05
     endif
    
  

  enddo
 enddo
  do i1=1,N_x1+1
     do i2=1,N_x2+1
       Flux(i1,i2)=-(alpha*Flux_x1(i1,i2)+beta*Flux_x2(i1,i2))
     enddo
   enddo
endif!order
end subroutine compute_flux



subroutine compute_translate_nodes_periodic(alpha,N_cells,old_node_positions,new_node_positions)
  ! compute displaced nodes in the case of a translation
  ! the nodes are put in [x_min,x_max] by periodicity
  use sll_constants
  implicit none

  sll_int,intent(in) :: N_cells
  sll_real64,dimension(1:N_cells+1) :: old_node_positions
  sll_real64,dimension(1:N_cells+1) :: new_node_positions
  sll_int :: i
  sll_real64 :: alpha,x_min,x_max,xx  
  x_min = old_node_positions(1)
  x_max = old_node_positions(N_cells+1)
  
  do i=1,N_cells+1
    xx = (old_node_positions(i)-alpha-x_min)/(x_max-x_min)
    do while(xx>=1._f64)
      xx = xx-1._f64
    enddo
    do while(xx<0._f64)
        xx = xx+1._f64
    enddo
    if(xx>1._f64)then
      print *,'Problem of localization',i,xx
      stop
    endif
    new_node_positions(i) = x_min+xx*(x_max-x_min)
  enddo
  
  
end subroutine compute_translate_nodes_periodic


subroutine compute_non_unif_integral(x_points,f_points,N_points,val)
  use sll_constants
  implicit none
  sll_real64,intent(out) :: val
  sll_int,intent(in) :: N_points
  sll_real64,dimension(1:N_points) :: x_points,f_points
  sll_int :: i
  sll_real64 :: tmp,x1,x2,fval1,fval2
  val = 0._f64
  if(N_points<=1)then
    print *,'bad value of N_points=',N_points
    stop
  endif
  do i=1,N_points-1
    x1 = x_points(i)
    x2 = x_points(i+1)
    if(x2<x1)then
      print *,i,'bad integration points x1=',x1,'x2=',x2
      stop
    endif
    fval1 = f_points(i)
    fval2 = f_points(i+1)
    tmp = 0.5_f64*(fval1+fval2)*(x2-x1)
    val=val+tmp
  enddo
  
  
end  subroutine compute_non_unif_integral


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



