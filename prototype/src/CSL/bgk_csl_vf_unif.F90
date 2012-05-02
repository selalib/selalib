program bgk_csl
#include "sll_working_precision.h"
#include "sll_mesh_types.h"
#include "sll_memory.h"
  use numeric_constants
  use distribution_function
  use sll_diagnostics
  use sll_csl
  use sll_splines
  implicit none
  external compute_translate_nodes_periodic,compute_non_unif_integral
  !external compute_flux
  !external poisson1dpertrap
  sll_real64 :: x2_min,x2_max,x1_min,x1_max,x1,x2,delta_x1,delta_x2
  sll_real64 :: mu,xi,L,H
  sll_real64 :: Flux_x1_p, Flux_x1_m,Flux_x2_p,Flux_x2_m
  !sll_int :: i1m1, i1p1, i2m2, i2p2, i1m2,i1p2,&
          !   i2m3,i2p3,i1m3,i1p3
  sll_real ::  fterm1,fterm2,fterm3,fterm4
  sll_int    :: i,j,N_phi,err,N_x1,N_x2,i1,i2,N,nb_step
  LOGICAL :: ex
    !type(sll_distribution_function_2D_t), pointer :: dist_func
  sll_real64,dimension(:), pointer :: phi,node_positions_x1,node_positions_x2
  sll_real64,dimension(:), pointer :: abar_x1,abar_x2
  sll_real64,dimension(:), pointer :: new_node_positions,buf_1d,rho,E,E_store, rho_exact
  sll_real64,dimension(:,:), pointer :: f,f_store, chi,sigma, Flux_x1,Flux_x2
  sll_real64 :: phi_val,delta_x1_phi,xx,dt,alpha,val,beta
  sll_int :: i1m3,i1m2,i1m1,i1p3,i1p2,i1p1,i2m3,i2m2,i2m1,i2p3,i2p2,i2p1
  sll_int :: ii,step,test_case,order, visu_step
  type(cubic_nonunif_spline_1D), pointer :: spl_per_x1, spl_per_x2
  
  

  N_x1 =128!256! 128!128
  N_x2 =128!256! 128!256
  dt = 0.01_f64
  nb_step = 600
  order=3
  test_case = 4
  visu_step = 10
  
  N = max(N_x1,N_x2)
  
  print *,'#max(N_x1,N_x2)=',N
  
  
  x2_max = 10._f64
  x2_min = -x2_max
  
  SLL_ALLOCATE(f(N_x1+1,N_x2+1),err)
  SLL_ALLOCATE(f_store(N_x1+1,N_x2+1),err)
  SLL_ALLOCATE(node_positions_x1(N_x1+1),err)
  SLL_ALLOCATE(node_positions_x2(N_x2+1),err)
  SLL_ALLOCATE(new_node_positions(N+1),err)
  SLL_ALLOCATE(buf_1d(N+1),err)
  SLL_ALLOCATE(rho(N_x1+1),err)
  SLL_ALLOCATE(rho_exact(N_x1+1),err)
  SLL_ALLOCATE(E(N_x1+1),err)
  SLL_ALLOCATE(E_store(N_x1+1),err)
  SLL_ALLOCATE(chi(N_x1+1+2,N_x2+1+2),err)
  SLL_ALLOCATE(sigma(N_x1+1+2,N_x2+1+2),err)
  SLL_ALLOCATE(Flux_x1(N_x1+1,N_x2+1),err)
  SLL_ALLOCATE(Flux_x2(N_x1+1,N_x2+1),err)
  SLL_ALLOCATE(abar_x1(N_x1+1),err)
  SLL_ALLOCATE(abar_x2(N_x1+1),err)
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
          !N_x1=6  
        !do i1=1,N_x1+1
           !print *,"mod2",i1-1, i1,i1+1
           !print *,"mod1", modulo(i1-1,N_x1)+1, modulo(i1,N_x1)+1, modulo(i1+1,N_x1)+1
           
        !enddo
          !print *,N_x1, modulo(1,N_x1),modulo(N_x1,N_x1),modulo(0,N_x1),modulo(-1,N_x1)
  !time iteration
    !N_x1=4
    !print *,"mod1", modulo(-2,N_x1), modulo(-1,N_x1), modulo(0,N_x1),modulo(1,N_x1),modulo(4,N_x1)
    !stop
   !stop
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
  !E(0)=0!*****
  !do i1=2,N_x1+1
  !E(i1)=-(rho(i1)-rho(i1-1))/delta_x1!*******
  !enddo
  do step=1,nb_step
    
    
    f_store = f
    !compute E
    do i1=1,N_x1+1
      buf_1d(1:N_x2+1) = f_store(i1,1:N_x2+1)
      call compute_non_unif_integral(node_positions_x2(1:N_x2+1),buf_1d(1:N_x2+1),N_x2+1,val)
      rho(i1)=val-1._f64
  
    enddo
    !rho=rho_exact
    E=rho
    call poisson1dpertrap(E,x1_max-x1_min,N_x1)
    
    E_store = E
    
     !diagnostic
    val=0._f64
    do i1=1,N_x1
      val = val+E(i1)*E(i1)
    enddo
    val = val/real(N_x1,f64)
    !print *,step, (real(step,f64)-1._f64)*dt,val
    
     
   if(step==1) then

   do i1=1,N_x1+1
      x1 = x1_min+real(i1-1,f64)*delta_x1
      write(*,*) x1,E(i1),rho(i1)!,rho_exact(i1)
    enddo
    stop
  endif
    !diagnostic
    ! do i1=1,N_x1
      !print *,x1_min+real(i1-1)*delta_x1,E(i1),rho(i1)
    !enddo

    !compute chi(uij) and sigma(uij)
      do i1=1,N_x1+1
        x1 = x1_min+real(i1-1,f64)*delta_x1
        do i2=1,N_x2+1
           x2 = x2_min+real(i2-1,f64)*delta_x2
           chi(i1,i2)   =x2*f_store(i1,i2)
           sigma(i1,i2) =E_store(i1)*f_store(i1,i2)
       enddo
   enddo
    chi(N_x1+2,1:N_x2+1+2)=chi(2,1:N_x2+1+2)
    chi(N_x1+3,1:N_x2+1+2)=chi(3,1:N_x2+1+2)
    chi(1:N_x1+1+2,N_x2+2)=chi(1:N_x1+1+2,2)
    chi(1:N_x1+1+2,N_x2+3)=chi(1:N_x1+1+2,3)
    sigma(N_x1+2,1:N_x2+1+2)=sigma(2,1:N_x2+1+2)
    sigma(N_x1+3,1:N_x2+1+2)=sigma(3,1:N_x2+1+2)
    sigma(1:N_x1+1+2,N_x2+2)=sigma(1:N_x1+1+2,2)
    sigma(1:N_x1+1+2,N_x2+3)=sigma(1:N_x1+1+2,3)
  ! For fixed  i2 
 do i2=1,N_x2+1
     do i1=1,N_x1
     abar_x1(i1)=0!(chi(i1+1,i2)-chi(i1,i2))/(f_store(i1+1,i2)-f_store(i1,i2))
     
     enddo
  
     abar_x1(N_x1+1)=abar_x1(1)
  do i1=1,N_x1+1
     !abar_x1(i1)=(chi(i1+1,i2)-chi(i1,i2))/(f_store(i1+1,i2)-f_store(i1,i2))
    if(order==3) then
     i1m2=i1-2
     i1m1=i1-1
     i1p1=i1+1
     i1p2=i1+2
     if (i1m1 <=0) then 
         i1m1=i1m1+N_x1+1
     endif
     if (i1p1 >=N_x1+1) then 
         i1p1=i1p1-N_x1-1
     endif
     if (i1m2 <=0) then 
         i1m2=i1m2+N_x1+1        
     endif
     if (i1p2 >=N_x1+1) then 
         i1p2=i1p2-N_x1-1
     endif
   
    if(abar_x1(i1m2)>=0) then     
     fterm1=chi(i1m2,i2)
    else 
     fterm1=chi(i1m2+1,i2)
    endif

     if(abar_x1(i1m1)>=0) then     
     fterm2=chi(i1m1,i2)
    else 
     fterm2=chi(i1m1+1,i2)
    endif
     
    if(abar_x1(i1)>=0) then     
     fterm3=chi(i1,i2)
    else 
     fterm3=chi(i1+1,i2)
    endif

     if(abar_x1(i1p1)>=0) then     
     fterm4=chi(i1p1,i2)
    else 
     fterm4=chi(i1p1+1,i2)
    endif
     !fterm1=chi(i1left2,i2)
     !fterm2=chi(i1left,i2)
     !fterm3=chi(i1,i2)
     !fterm4=chi(i1right,i2)
    
    
      Flux_x1(i1,i2)=(1._f64/6._f64)*fterm1-fterm2+(1._f64/2._f64)*fterm3&
       +(1._f64/3._f64)*fterm4
     !Flux_x1_p=(1._f64/6._f64)*fterm1_p-fterm2_p+(1._f64/2._f64)*fterm3_p&
      ! +(1._f64/3._f64)*fterm4_p

     !Flux_x1_m=(1._f64/6._f64)*fterm1_m-fterm2_m+(1._f64/2._f64)*fterm3_m&
     !  +(1._f64/3._f64)*fterm4_m

    ! Flux_x1(i1,i2)=(1._f64/6._f64)*chi(i1left2,i2)-chi(i1left,i2)+(1._f64/2._f64)*chi(i1,i2)&
      ! +(1._f64/3._f64)*chi(i1right,i2)
  endif   
    !Flux_x1(i1,1:N_x2+1)=-(1._f64/6._f64)*chi(i1left,1:N_x2+1)+(5._f64/6._f64)*chi(i1,1:N_x2+1)&
      ! +(1._f64/3._f64)*chi(i1right,1:N_x2+1)
    !Flux_x1(i1,1:N_x2+1)=-(1._f64/6._f64)*chi(i1-1,1:N_x2+1)+(5._f64/6._f64)*chi(i1,1:N_x2+1)+(1._f64/3._f64)*chi(i1+1,1:N_x2+1)
  if(order==5) then
     i1m3=modulo(i1-3,N_x1)
     i1m2=modulo(i1-2,N_x1)
     i1m1=modulo(i1-1,N_x1)
     i1p1=modulo(i1+1,N_x1)
     i1p2=modulo(i1+2,N_x1)
     i1p3=modulo(i1+3,N_x1)
     if(i1m3==0)  i1m3=N_x1
     if(i1m2==0)  i1m2=N_x1
     if(i1m1==0)  i1m1=N_x1
     if(i1p1==0)  i1p1=N_x1+1
     if(i1p2==0)  i1p2=N_x1+1
     if(i1p3==0)  i1p3=N_x1+1
     
     Flux_x1_p=(1._f64/30._f64)*chi(i1-2,i2)-(13._f64/60._f64)*chi(i1-1,i2)+(47._f64/60._f64)*chi(i1,i2)&
                +(9._f64/20._f64)*chi(i1+1,i2)-(1._f64/20._f64)*chi(i1+2,i2)
     Flux_x1_m=(1._f64/30._f64)*chi(i1-3,i2)-(13._f64/60._f64)*chi(i1-2,i2)+(47._f64/60._f64)*chi(i1-1,i2)&
                +(9._f64/20._f64)*chi(i1,i2)-(1._f64/20._f64)*chi(i1+1,i2)
    Flux_x1(i1,i2)=Flux_x1_p-Flux_x1_m
  endif
  enddo
 enddo 
  
 !For fixed i1
 do i1=1,N_x1+1
     do i2=1,N_x2
     abar_x2(i2)=E(i1)!(sigma(i1,i2+1)-sigma(i1,i2))/(f_store(i1,i2+1)-f_store(i1,i2))
     !print*,"aa", x2_min+(i2-1)*delta_x2,abar_x2(i2)
     enddo
    !stop
     abar_x2(N_x2+1)=abar_x2(1)
  do i2=1,N_x2+1
   if(order==3) then
     i2m2=i2-2
     i2m1=i2-1
     i2p1=i2+1
     i2p2=i2+1
     if (i2m1 <=0) then 
         i2m1=i2m1+N_x2+1        
     endif
     if (i2p1 >=N_x2+1) then 
         i2p1=i2p1-N_x2-1
     endif

     if (i2m2 <=0) then 
         i2m2=i2m2+N_x2+1        
     endif
     if (i2p2 >=N_x2+1) then 
         i2p2=i2p2-N_x2-1
     endif


     if(abar_x2(i2m2)>=0) then     
     fterm1=sigma(i1,i2m2)
    else 
     fterm1=sigma(i1,i2m2+1)
    endif

     if(abar_x2(i1m2)>=0) then     
     fterm2=sigma(i1,i2m1)
    else 
     fterm2=sigma(i1,i2m1+1)
    endif
     
    if(abar_x2(i2)>=0) then     
     fterm3=sigma(i1,i2)
    else 
     fterm3=sigma(i1,i2+1)
    endif

     if(abar_x2(i2p1)>=0) then     
     fterm4=sigma(i1,i2p1)
    else 
     fterm4=sigma(i1,i2p1+1)
    endif
     !fterm1=chi(i1left2,i2)
     !fterm2=chi(i1left,i2)
     !fterm3=chi(i1,i2)
     !fterm4=chi(i1right,i2)
    
    
      Flux_x2(i1,i2)=(1._f64/6._f64)*fterm1-fterm2+(1._f64/2._f64)*fterm3&
       +(1._f64/3._f64)*fterm4

      !Flux_x2(i1,i2)=(1._f64/6._f64)*sigma(i1,i2left2)-sigma(i1,i2left)+(1._f64/2._f64)*sigma(i1,i2)&
      !+(1._f64/3._f64)*sigma(i1,i2right)
  end if
    !Flux_x2(1:N_x1+1,i2)=-(1._f64/6._f64)*sigma(1:N_x1+1,i2left)+(5._f64/6._f64)*sigma(1:N_x1+1,i2)&
      !+(1._f64/3._f64)*sigma(1:N_x1+1,i2right)
   ! Flux_x2(1:N_x1+1,i2)=-(1._f64/6._f64)*sigma(1:N_x1+1,i2-1)+(5._f64/6._f64)*sigma(1:N_x1+1,i2)+(1._f64/3._f64)*sigma(1:N_x1+1,i2+1)
   if(order==5) then
     i2m3=modulo(i2-3,N_x2)
     i2m2=modulo(i2-2,N_x2)
     i2m1=modulo(i2-1,N_x2)
     i2p1=modulo(i2+1,N_x2)
     i2p2=modulo(i2+2,N_x2)
     i2p3=modulo(i2+3,N_x2)
     if(i2m3==0)  i2m3=N_x1
     if(i2m2==0)  i2m2=N_x1
     if(i2m1==0)  i2m1=N_x1
     if(i2p1==0)  i2p1=N_x1+1
     if(i2p2==0)  i2p2=N_x1+1
     if(i2p3==0)  i2p3=N_x1+1
     Flux_x2_p=(1._f64/30._f64)*sigma(i1,i2-2)-(13._f64/60._f64)*sigma(i1,i2-1)+(47._f64/60._f64)*sigma(i1,i2)&
                +(9._f64/20._f64)*sigma(i1,i2+1)-(1._f64/20._f64)*sigma(i1,i2+2)
     Flux_x2_m=(1._f64/30._f64)*sigma(i1,i2-3)-(13._f64/60._f64)*sigma(i1,i2-2)+(47._f64/60._f64)*sigma(i1,i2-1)&
                +(9._f64/20._f64)*sigma(i1,i2)-(1._f64/20._f64)*sigma(i1,i2+1)
    Flux_x2(i1,i2)=Flux_x2_p-Flux_x2_m
   endif
  enddo
 enddo
 alpha=(dt/delta_x1)
 beta =(dt/delta_x2)
   do i1=1,N_x1+1
     do i2=1,N_x2+1
         !Flux_x1=-(1._f64/6._f64)*chi(i1-1,i2)+(5._f64/6._f64)*chi(i1,i2)+(1._f64/3._f64)*chi(i1+1,i2)
         !Flux_x2=-(1._f64/6._f64)*sigma(i1-1,i2)+(5._f64/6._f64)*sigma(i1,i2)+(1._f64/3._f64)*sigma(i1+1,i2)
        f(i1,i2)=f_store(i1,i2)-alpha*Flux_x1(i1,i2)-beta*Flux_x2(i1,i2)
     enddo
   enddo
  !do i1 = 1,N_x1 
     ! do i2 = 1, N_x2    
       ! call sll_set_df_val(dist_func, i1, i2,f(i1,i2))
      !enddo  
    !enddo
    !if(mod(step,visu_step)==0)then
    !  call write_distribution_function ( dist_func )
    !endif
   
 enddo !end time step
    
    
  !enddo !end time step

  open(unit=900,file='field_final_vf.dat')  
    do i1=1,N_x1+1
      x1 = x1_min+real(i1-1,f64)*delta_x1
      write(900,*) x1,E(i1),rho(i1)
    enddo
  close(900)


  
end program


 !subroutine Compute_flux(E,L,N1,N2)
 ! use numeric_constants
  !implicit none
 ! sll_int,intent(in)::N1,N2
  !sll_real64,dimension(N1+1)::E
  !sll_real64,intent(in)::L
  !sll_int::i
  !sll_real64::eold,enew,dx2,tmp

!end subroutine compute_flux



subroutine compute_translate_nodes_periodic(alpha,N_cells,old_node_positions,new_node_positions)
  ! compute displaced nodes in the case of a translation
  ! the nodes are put in [x_min,x_max] by periodicity
  use numeric_constants
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
  use numeric_constants
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
  use numeric_constants
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



