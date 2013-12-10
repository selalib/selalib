program vp_non_unif_mesh
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
!contact: mehrenbe@math.unistra.fr for this  program

  use sll_constants
  use cubic_non_uniform_splines
  use bgk_mesh_construction
  use contrib_rho_module
  use classical_conservative_semi_lagrangian
  implicit none
  
  !parameters
  sll_int32 :: test_case,mesh_case,rho_case,div_case
  sll_int32 :: nb_step,N_x1,N_x2,nc_eta1,nc_eta2
  sll_real64 :: alpha_mesh,dt,velocity
  !other variables
  sll_int32 :: N,err  
  sll_real64 :: x1_min,x1_max,x2_min,x2_max
  type(cubic_nonunif_spline_1D), pointer :: spl_per_x1, spl_per_x2
   type(cubic_nonunif_spline_1D), pointer :: spl_per_x1_rho,spl_per_x1_E
  sll_real64 ::geom_x(2,2),geom_eta(2,2)
  sll_real64,dimension(:), pointer :: rho,E,E_fine,phi_poisson,buf1d
  sll_real64,dimension(:), pointer :: Xstar,node_positions_x1,node_positions_x2
  sll_real64,dimension(:), pointer :: node_positions_x1_dual,node_positions_x2_dual
  sll_real64,dimension(:), pointer :: node_positions_x1_poisson
  sll_real64,dimension(:,:),pointer::f,f_init,f_store
  !sll_real64,dimension(:,:),pointer::x2n_array,x2c_array
  !sll_real64,dimension(:,:),pointer::jac_array
  sll_real64, dimension(:,:), pointer :: a1,a2,psi
  sll_real64,dimension(:,:),pointer::integration_points
  sll_int32  :: i1,i2,i,step,N_alpha_x2,j1,N_x1_poisson,dual_case_x1,dual_case_x2
  sll_real64 :: delta_x1,delta_x2,x1,x2
  sll_real64 :: eta1_min,eta1_max,eta2_min,eta2_max,delta_eta1,delta_eta2
  sll_real64 :: val,tmp,tmp0
  sll_real64 :: alpha_x2,shift
  sll_real64,dimension(:), pointer :: node_positions_x1_tmp,node_positions_x2_tmp  
  !for mesh_case 2: [zone_1,zone_2,zone_3]=
  !                 [(N_x2-N_alpha_x2)/2,N_alpha_x2,(N_x2-N_alpha_x2)/2]=[(1-alpha)/2,alpha,(1-alpha)/2]
  

  
  
  
  test_case = 4
  mesh_case = 5
  rho_case = 5
  div_case = 1
  
  dual_case_x1 = 1
  dual_case_x2 = 1

  nb_step = 600
  N_x1 = 64
  N_x2 = 250
  
  N_x1_poisson = 2048!N_x1!
  
  alpha_x2 = 0.25_f64
  N_alpha_x2 = 100 ! for finer mesh inside should be greater than N_x2*alpha_x2
  

  alpha_mesh = 0.1e-3_f64 !0.1_f64
  dt = 0.1_f64
  velocity = 1._f64
  

  N = max(N_x1,N_x2)

  nc_eta1 = N_x1
  nc_eta2 = N_x2


  print *,'#max(N_x1,N_x2)=',N
  
  SLL_ALLOCATE(node_positions_x1_tmp(N+1), err)
  SLL_ALLOCATE(node_positions_x2_tmp(N+1), err)
  
  SLL_ALLOCATE(f(N_x1+1,N_x2+1),err)
  SLL_ALLOCATE(f_init(N_x1+1,N_x2+1),err)
  SLL_ALLOCATE(f_store(N_x1+1,N_x2+1),err)
  
  SLL_ALLOCATE(a1(N_x1+1,N_x2+1),err)
  SLL_ALLOCATE(a2(N_x1+1,N_x2+1),err)
  SLL_ALLOCATE(psi(N_x1+1,N_x2+1),err)
  
  
  SLL_ALLOCATE(integration_points(N_x1+1,N_x2+1),err)
  
  SLL_ALLOCATE(rho(N_x1+1),err)
  SLL_ALLOCATE(E_fine(N_x1_poisson+1),err)
  SLL_ALLOCATE(E(N_x1+1),err)
  SLL_ALLOCATE(node_positions_x1_poisson(N_x1_poisson+1),err)
  
  SLL_ALLOCATE(phi_poisson(N_x1+1),err)
  
  SLL_ALLOCATE(buf1d(N+1),err)
  SLL_ALLOCATE(Xstar(N+1),err)
  
  SLL_ALLOCATE(node_positions_x1(N_x1+1),err)
  SLL_ALLOCATE(node_positions_x2(N_x2+1),err)
  
  SLL_ALLOCATE(node_positions_x1_dual(N_x1+1),err)
  SLL_ALLOCATE(node_positions_x2_dual(N_x2+1),err)
  
  
  if(test_case>=4)then
    x1_min = 0._f64
    x1_max = 4._f64*sll_pi
    x2_min = -10._f64!-6._f64
    x2_max = -x2_min
  endif

  eta1_min =  0.0_f64
  eta1_max = 1.0_f64 ! 0.15_f64*x1_max! 1.0_f64
  eta2_min =  0.0_f64
  eta2_max =  1.0_f64
  
  geom_x(1,1)=x1_min
  geom_x(2,1)=x1_max
  geom_x(1,2)=x2_min
  geom_x(2,2)=x2_max
    

  geom_eta(1,1)=eta1_min
  geom_eta(2,1)=eta1_max
  geom_eta(1,2)=eta2_min
  geom_eta(2,2)=eta2_max


  delta_x1 = (x1_max-x1_min)/real(N_x1,f64)
  delta_x2 = (x2_max-x2_min)/real(N_x2,f64)

  delta_eta1 = (eta1_max-eta1_min)/real(nc_eta1,f64)
  delta_eta2 = (eta2_max-eta2_min)/real(nc_eta2,f64)

  do i=1,N_x1_poisson+1
    node_positions_x1_poisson(i) = (real(i,f64)-1._f64)/real(N_x1_poisson,f64)
  enddo

  
  !we suppose uniform mesh for the interpolation
  if(mesh_case==1)then
    do i=1,N_x1+1
      node_positions_x1(i) = x1_min+(real(i,f64)-1._f64)*delta_x1
    enddo
    do i=1,N_x2+1
      node_positions_x2(i) = x2_min+(real(i,f64)-1._f64)*delta_x2
    enddo
  endif

  if(mesh_case==2)then
    do i=1,N_x1+1
      node_positions_x1(i) = x1_min+(real(i,f64)-1._f64)*delta_x1
    enddo
    !k=floor(0.5_f64*(1._f64-alpha_x2)*N)
    !if((k<0).or.(k>=N))then
    !  print *, "bad value of alpha_x2",alpha_x2,k
    !  stop
    !endif
    ! [0,k,N-k,N]

    do i=0, (N_x2-N_alpha_x2)/2
      node_positions_x2(i+1) = real(i,f64)/real((N_x2-N_alpha_x2)/2,f64)*0.5_f64*(1._f64-alpha_x2)
    enddo
    do i=1,N_alpha_x2
      node_positions_x2(i+(N_x2-N_alpha_x2)/2+1) = &
        0.5_f64*(1._f64-alpha_x2)+alpha_x2*real(i,f64)/real(N_alpha_x2,f64)
    enddo
    do i=1,(N_x2-N_alpha_x2)/2
      node_positions_x2(i+(N_x2+N_alpha_x2)/2+1) = &
        0.5_f64*(1._f64+alpha_x2)+real(i,f64)/real((N_x2-N_alpha_x2)/2,f64)*0.5_f64*(1._f64-alpha_x2)
    enddo
    
    !do i=1,N_x2+1
    !  print *,i,node_positions_x2(i)
    !enddo
    !stop       
    do i=1,N_x2+1
      node_positions_x2(i) = x2_min+node_positions_x2(i)*(x2_max-x2_min)
    enddo
  endif


  if(mesh_case==3)then
    do i=1,N_x2+1
      node_positions_x2(i) = x2_min+(real(i,f64)-1._f64)*delta_x2
    enddo
    !k=floor(0.5_f64*(1._f64-alpha_x2)*N)
    !if((k<0).or.(k>=N))then
    !  print *, "bad value of alpha_x2",alpha_x2,k
    !  stop
    !endif
    ! [0,k,N-k,N]

    do i=0, (N_x1-N_alpha_x2)/2
      node_positions_x1(i+1) = real(i,f64)/real((N_x1-N_alpha_x2)/2,f64)&
        *0.5_f64*(1._f64-alpha_x2)
    enddo
    do i=1,N_alpha_x2
      node_positions_x1(i+(N_x1-N_alpha_x2)/2+1) = &
        0.5_f64*(1._f64-alpha_x2)+alpha_x2*real(i,f64)/real(N_alpha_x2,f64)
    enddo
    do i=1,(N_x1-N_alpha_x2)/2
      node_positions_x1(i+(N_x1+N_alpha_x2)/2+1) = &
        0.5_f64*(1._f64+alpha_x2)+real(i,f64)/real((N_x1-N_alpha_x2)/2,f64)&
        *0.5_f64*(1._f64-alpha_x2)
    enddo
    
    !do i=1,N_x2+1
    !  print *,i,node_positions_x2(i)
    !enddo
    !stop       
    do i=1,N_x1+1
      node_positions_x1(i) = x1_min+node_positions_x1(i)*(x1_max-x1_min)
    enddo
  endif


  if(mesh_case==4)then

  !initialization of X
  node_positions_x1(1)=0._f64
  do i=2,N_x1+1
    node_positions_x1(i)=real(i-1,f64)*delta_x1
    !Xmesh(i)=xmin+(real(i,rk)+0.5*sin(2._rk*M_PI*real(i,rk)*dx))*dx
    !Xmesh(i)=Xmesh(i-1)+(2._rk+sin(2._rk*M_PI*real(i,rk)/real(N+1,rk)))*3._rk
    node_positions_x1(i)=node_positions_x1(i-1)+(real(i-1,f64)*real(N_x1+1-(i-1),f64))
    !Xmesh(i)=(real(mod(i+N/2,N),rk)/real(N,rk))**2
    !if(i>0.and.Xmesh(i)<Xmesh(i-1))then
    !  Xmesh(i)=Xmesh(i)+1._rk
    !endif
  enddo
  node_positions_x1=node_positions_x1/node_positions_x1(N_x1+1)
  j1=N_x1/3;shift=0._f64
  do i=0,N_x1    
    node_positions_x1_tmp(i+1)=node_positions_x1(j1+1)+shift
    j1=j1+1
    if(j1>=N_x1)then
      j1=j1-N_x1;shift=shift+1._f64
    endif
  enddo
  node_positions_x1=node_positions_x1_tmp-node_positions_x1_tmp(1)
  
  !node_positions_x1(1)=0._f64
  !node_positions_x1(2)=1._f64
  !node_positions_x1(3)=2._f64
  !node_positions_x1(4)=node_positions_x1(3)+9.67062069240196_f64
  node_positions_x1=node_positions_x1/node_positions_x1(N_x1+1)
  
  !do i=1,N_x1+1
  !  print *,i,node_positions_x1(i)
  !enddo
  
  !stop


    do i=1,N_x1+1
      node_positions_x1(i) = x1_min+node_positions_x1(i)*(x1_max-x1_min)
    enddo


    do i=1,N_x2+1
      node_positions_x2(i) = x2_min+(real(i,f64)-1._f64)*delta_x2
    enddo
  endif
  
  
  if(mesh_case==5)then

  !initialization of X
  node_positions_x2(1)=0._f64
  do i=2,N_x2+1
    node_positions_x2(i)=real(i-1,f64)*delta_x2
    !Xmesh(i)=xmin+(real(i,rk)+0.5*sin(2._rk*M_PI*real(i,rk)*dx))*dx
    !Xmesh(i)=Xmesh(i-1)+(2._rk+sin(2._rk*M_PI*real(i,rk)/real(N+1,rk)))*3._rk
    node_positions_x2(i)=node_positions_x2(i-1)+(real(i-1,f64)*real(N_x2+1-(i-1),f64))
    !Xmesh(i)=(real(mod(i+N/2,N),rk)/real(N,rk))**2
    !if(i>0.and.Xmesh(i)<Xmesh(i-1))then
    !  Xmesh(i)=Xmesh(i)+1._rk
    !endif
  enddo
  node_positions_x2=node_positions_x2/node_positions_x2(N_x2+1)
  j1=N_x2/3;shift=0._f64
  do i=0,N_x2    
    node_positions_x2_tmp(i+1)=node_positions_x2(j1+1)+shift
    j1=j1+1
    if(j1>=N_x2)then
      j1=j1-N_x2;shift=shift+1._f64
    endif
  enddo
  node_positions_x2=node_positions_x2_tmp-node_positions_x2_tmp(1)
  
  !node_positions_x1(1)=0._f64
  !node_positions_x1(2)=1._f64
  !node_positions_x1(3)=2._f64
  !node_positions_x1(4)=node_positions_x1(3)+9.67062069240196_f64
  node_positions_x2=node_positions_x2/node_positions_x2(N_x2+1)
  
  !do i=1,N_x1+1
  !  print *,i,node_positions_x1(i)
  !enddo
  
  !stop


    do i=1,N_x2+1
      node_positions_x2(i) = x2_min+node_positions_x2(i)*(x2_max-x2_min)
    enddo


    do i=1,N_x1+1
      node_positions_x1(i) = x1_min+(real(i,f64)-1._f64)*delta_x1
    enddo
  endif


  
    
  node_positions_x1 = (node_positions_x1-x1_min)/(x1_max-x1_min)
  node_positions_x2 = (node_positions_x2-x2_min)/(x2_max-x2_min)
  
  do i1=2,N_x1
    node_positions_x1_dual(i1) =  0.5_f64*(node_positions_x1(i1)+node_positions_x1(i1-1))
  enddo
  node_positions_x1_dual(1) =  0.5_f64*(node_positions_x1(1)+node_positions_x1(N_x1)-1._f64)
  node_positions_x1_dual(N_x1+1) = node_positions_x1_dual(1) + 1._f64
  
  tmp = node_positions_x1_dual(1)
  node_positions_x1_dual = node_positions_x1_dual-tmp

  do i2=2,N_x2
    node_positions_x2_dual(i2) =  0.5_f64*(node_positions_x2(i2)+node_positions_x2(i2-1))
  enddo
  node_positions_x2_dual(1) =  0.5_f64*(node_positions_x2(1)+node_positions_x2(N_x2)-1._f64)
  node_positions_x2_dual(N_x2+1) = node_positions_x2_dual(1) + 1._f64
  
  tmp = node_positions_x2_dual(1)
  node_positions_x2_dual = node_positions_x2_dual-tmp
  
  
  
  
  !print *,node_positions_x1
  !print *,node_positions_x2
  
  
  

  !construction of mesh
  
   !call construct_bgk_mesh(nc_eta1,nc_eta2,mesh_case,&
   !  &x1n_array,x2n_array,x1c_array,x2c_array,jac_array,integration_points,&
   !  &geom_x,geom_eta,alpha_mesh,N_x1,N_x2)

  

  spl_per_x1 =>  new_cubic_nonunif_spline_1D( N_x1, SLL_PERIODIC)
  spl_per_x2 =>  new_cubic_nonunif_spline_1D( N_x2, SLL_PERIODIC)

  spl_per_x1_rho =>  new_cubic_nonunif_spline_1D( N_x1, SLL_PERIODIC)
  spl_per_x1_E =>  new_cubic_nonunif_spline_1D( N_x1_poisson, SLL_PERIODIC)
  do i1=1,nc_eta1+1
    do i2=1,nc_eta2+1
      x1 = x1_min+(x1_max-x1_min)*node_positions_x1(i1)
      x2 = x2_min+(x2_max-x2_min)*node_positions_x2(i2)
      if(test_case==3)then
        val = exp(-0.5_f64*40._f64*((x1-.5_f64)**2+(x2-.5_f64)**2))
      endif
      if(test_case==4)then
        !linear landau damping
        val = 1._f64/(sqrt(2._f64*sll_pi))*exp(-0.5_f64*x2*x2)
        !f_equil(i1,i2) = val*jac_array(i1,i2)
        val = val*(1._f64+0.001_f64*cos(2._f64*sll_pi/(x1_max-x1_min)*x1))
      endif
      if(test_case==5)then
        !non linear landau damping
        val = 1._f64/(sqrt(2._f64*sll_pi))*exp(-0.5_f64*x2*x2)
        val = val*(1._f64+0.5_f64*cos(2._f64*sll_pi/(x1_max-x1_min)*x1))
      endif
      if(test_case==6.or.test_case==7)then
        !gaussian equilibrium
        val = 1._f64/(sqrt(2._f64*sll_pi))*exp(-0.5_f64*x2*x2)
        !f_equil(i1,i2) = val*jac_array(i1,i2)
      endif
      f_init(i1,i2) = val!*jac_array(i1,i2)      
    enddo
  enddo
  
  
  
  f = f_init
  
  !compute rho init  
  integration_points(1,1:N_x2+1) = x2_min+(x2_max-x2_min)*node_positions_x2(1:N_x2+1)   
  

  do i1 = 1, N_x1+1 
    integration_points(2,1:N_x2+1) = f(i1,1:N_x2+1) 
    !rho(i1) = compute_non_unif_integral(integration_points,N_x2+1,rho_case)
  enddo
 
  
  
  
  !compute E_init   
  call compute_spline_nonunif( rho, spl_per_x1_rho,node_positions_x1)  
  

  call interpolate_array_value_nonunif( node_positions_x1_poisson, &
    E_fine, N_x1_poisson, spl_per_x1_rho)  
  
  stop
  call poisson1dpertrap(E_fine,x1_max-x1_min,N_x1_poisson)
  call compute_spline_nonunif( E_fine, spl_per_x1_E,node_positions_x1_poisson)
  call interpolate_array_value_nonunif( node_positions_x1, &
    E, N_x1, spl_per_x1_E)
  
  
    
  !do i1=1,N_x1+1
  !  print *,i1,rho(i1)
  !enddo
  !stop
    
  !call compute_rho_mapped_mesh2(rho,f_init,integration_points,rho_case,&
  !  nc_eta1,nc_eta2,geom_eta,jac_array,spl_per_x1_rho)
  !compute E init and psi
  !call compute_psi(a1,a2,rho,nc_eta1,nc_eta2,psi,phi_poisson,E,&
  !geom_x,x1n_array,x2n_array,x1c_array,x2c_array,jac_array,delta_eta1,delta_eta2,div_case)

  !call advect
  
  
  !tmp=sum(rho(1:N_x1))*delta_x1

  !do i=1,N_x1+1
  !  print *,x1_min+real(i-1,f64)*delta_x1,rho(i),E(i)
  !enddo

  !stop

  f = f_init


  !test for constant advection
  !Xstar(1:N_x1+1) = node_positions_x1(1:N_x1+1)-velocity*dt/(x1_max-x1_min)
  !f_init(1:N_x1+1,1)=1._f64
  !do i=1,N_x1+1
  !  f_init(i,1)=exp(-10._f64*(real(i-1,f64)/real(N_x1,f64)-0.5_f64)**2)
  !enddo
  
  
  do step=1,nb_step
    !print *,"#step=",i
    !f_store = f
    
    !compute field at time t_n
    !call compute_rho_mapped_mesh2(rho,f,integration_points,rho_case,&
    !  nc_eta1,nc_eta2,geom_eta,jac_array,spl_per_x1_rho)
    !call compute_psi(a1,a2,rho,nc_eta1,nc_eta2,psi,phi_poisson,E,&
    !  geom_x,x1n_array,x2n_array,x1c_array,x2c_array,jac_array,delta_eta1,delta_eta2,div_case)

    !compute rho  
    integration_points(1,1:N_x2+1) = x2_min+(x2_max-x2_min)*node_positions_x2(1:N_x2+1)   
    do i1 = 1, N_x1+1 
      integration_points(2,1:N_x2+1) = f(i1,1:N_x2+1) 
      rho(i1) = compute_non_unif_integral(integration_points,N_x2+1,rho_case)
    enddo
    !compute E
    !E=rho
    !call poisson1dpertrap(E,x1_max-x1_min,N_x1)

    tmp = 0._f64
    do i1=1,N_x1
      tmp = tmp+rho(i1)
    enddo
    tmp = -tmp/real(N_x1,f64)
    do i1=1,N_x1+1
     rho(i1) = rho(i1)+tmp
    enddo
    


    tmp = 0._f64
    do i1=1,N_x1
      tmp = tmp+rho(i1)
    enddo
    tmp = tmp*delta_x1

    
    call compute_spline_nonunif( rho, spl_per_x1_rho,node_positions_x1)  
    call interpolate_array_value_nonunif( node_positions_x1_poisson, &
      E_fine, N_x1_poisson, spl_per_x1_rho)  
    call poisson1dpertrap(E_fine,x1_max-x1_min,N_x1_poisson+1)
    call compute_spline_nonunif( E_fine, spl_per_x1_E,node_positions_x1_poisson)
    call interpolate_array_value_nonunif( node_positions_x1, &
      E, N_x1+1, spl_per_x1_E)

    !val=0._f64
    !do i1=1,N_x1
    !  val = val+E(i1)*E(i1)
    !enddo
    !val = val/real(N_x1,f64)

    val=0._f64
    do i1=1,N_x1_poisson
      val = val+E_fine(i1)*E_fine(i1)
    enddo
    val = val/real(N_x1_poisson,f64)

    !compute int fv dv  
    tmp = 0._f64
    integration_points(1,1:N_x2+1) = x2_min+(x2_max-x2_min)*node_positions_x2(1:N_x2+1)   
    do i1 = 1, N_x1+1
      do i2=1,N_x2+1 
        integration_points(2,i2) = integration_points(1,i2)*f(i1,i2)         
      enddo
      tmp = tmp+compute_non_unif_integral(integration_points,N_x2+1,rho_case)   
    enddo
    if(step==1)then
      tmp0=tmp
    endif



    print *,(real(step,f64)-1._f64)*dt,val,(tmp-tmp0)/tmp0,tmp,tmp0!-(x1_max-x1_min)

    
    ! advect in x dt/2 
    do i2=1,N_x2+1
      buf1d(1:N_x1+1) = f(1:N_x1+1,i2)
      if(dual_case_x1==0)then
        Xstar(1:N_x1+1) = node_positions_x1(1:N_x1+1) &
          -0.5_f64*(x2_min+(x2_max-x2_min)*node_positions_x2(i2))*dt/(x1_max-x1_min)
        call csl_advection_per(buf1d,spl_per_x1,Xstar,node_positions_x1,N_x1)
      endif
      if(dual_case_x1==1)then
        Xstar(1:N_x1+1) = node_positions_x1_dual(1:N_x1+1) &
          -0.5_f64*(x2_min+(x2_max-x2_min)*node_positions_x2(i2))*dt/(x1_max-x1_min)
        call csl_advection_per(buf1d,spl_per_x1,Xstar,node_positions_x1_dual,N_x1)
      endif
      f(1:N_x1+1,i2) = buf1d(1:N_x1+1)
      !print *,buf1d
      !stop

    enddo
    
    
    !compute rho  
    integration_points(1,1:N_x2+1) = x2_min+(x2_max-x2_min)*node_positions_x2(1:N_x2+1)   
    do i1 = 1, N_x1+1 
      integration_points(2,1:N_x2+1) = f(i1,1:N_x2+1) 
      rho(i1) = compute_non_unif_integral(integration_points,N_x2+1,rho_case)
    enddo
    tmp = 0._f64
    do i1=1,N_x1
      tmp = tmp+rho(i1)
    enddo
    tmp = -tmp/real(N_x1,f64)
    do i1=1,N_x1+1
     rho(i1) = rho(i1)+tmp
    enddo

    !compute E
    !E=rho
    !call poisson1dpertrap(E,x1_max-x1_min,N_x1)
    call compute_spline_nonunif( rho, spl_per_x1_rho,node_positions_x1)  
    call interpolate_array_value_nonunif( node_positions_x1_poisson, &
      E_fine, N_x1_poisson+1, spl_per_x1_rho)  
    call poisson1dpertrap(E_fine,x1_max-x1_min,N_x1_poisson)
    call compute_spline_nonunif( E_fine, spl_per_x1_E,node_positions_x1_poisson)
    call interpolate_array_value_nonunif( node_positions_x1, &
      E, N_x1+1, spl_per_x1_E)

    ! advect in v dt
    do i1=1,N_x1+1
      buf1d(1:N_x2+1) = f(i1,1:N_x2+1)
      if(dual_case_x2==0)then
        Xstar(1:N_x2+1) = node_positions_x2(1:N_x2+1) &
         -E(i1)*dt/(x2_max-x2_min)
        call csl_advection_per(buf1d,spl_per_x2,Xstar,node_positions_x2,N_x2)
      endif
      if(dual_case_x2==1)then
        Xstar(1:N_x2+1) = node_positions_x2_dual(1:N_x2+1) &
          -E(i1)*dt/(x2_max-x2_min)
        call csl_advection_per(buf1d,spl_per_x2,Xstar,node_positions_x2_dual,N_x2)
      endif
      f(i1,1:N_x2+1) = buf1d(1:N_x2+1)
    enddo
    
    ! advect in x dt/2 
    do i2=1,N_x2+1
      buf1d(1:N_x1+1) = f(1:N_x1+1,i2)
      if(dual_case_x1==0)then
        Xstar(1:N_x1+1) = node_positions_x1(1:N_x1+1) &
          -0.5_f64*(x2_min+(x2_max-x2_min)*node_positions_x2(i2))*dt/(x1_max-x1_min)
        call csl_advection_per(buf1d,spl_per_x1,Xstar,node_positions_x1,N_x1)
      endif
      if(dual_case_x1==1)then
        Xstar(1:N_x1+1) = node_positions_x1_dual(1:N_x1+1) &
          -0.5_f64*(x2_min+(x2_max-x2_min)*node_positions_x2(i2))*dt/(x1_max-x1_min)
        call csl_advection_per(buf1d,spl_per_x1,Xstar,node_positions_x1_dual,N_x1)
      endif
      f(1:N_x1+1,i2) = buf1d(1:N_x1+1)
    enddo


    
    !buf1d(1:N_x1) = f(1:N_x1,1)
    !call csl_advection_per(buf1d,spl_per_x1,Xstar,node_positions_x1,N_x1)
    !f(1:N_x1+1,1) = buf1d(1:N_x1+1)
  enddo
  !buf1d(1:N_x1+1) = f(1:N_x1+1,1)
  
  
  !do i=1,N_x1+1
  !  print *,i,f_init(i,1),buf1d(i) 
  !enddo
  

  print *,"#End of program"
  

end program


