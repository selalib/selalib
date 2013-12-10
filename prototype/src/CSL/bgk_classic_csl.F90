program bgk_classic_csl
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
  type(cubic_nonunif_spline_1D), pointer :: spl_per_x1, spl_per_x2, spl_per_x1_rho
  sll_real64 ::geom_x(2,2),geom_eta(2,2)
  sll_real64,dimension(:), pointer :: rho,E,phi_poisson,buf1d
  sll_real64,dimension(:), pointer :: Xstar,node_positions_x1,node_positions_x2
  sll_real64,dimension(:,:),pointer::f,f_init,f_store
  sll_real64,dimension(:,:),pointer::x1n_array,x2n_array,x1c_array,x2c_array
  sll_real64,dimension(:,:),pointer::jac_array
  sll_real64, dimension(:,:), pointer :: a1,a2,psi
  sll_real64,dimension(:,:,:),pointer::integration_points
  sll_int32  :: i1,i2,i,step
  sll_real64 :: delta_x1,delta_x2,x1,x2
  sll_real64 :: eta1_min,eta1_max,eta2_min,eta2_max,delta_eta1,delta_eta2,eta1c
  sll_real64 :: val

  
  
  
  test_case = 4
  mesh_case = 10
  rho_case = 2
  div_case = 1

  nb_step = 600
  N_x1 = 64
  N_x2 = 64

  alpha_mesh = 0.1e-3_f64 !0.1_f64
  dt = 0.1_f64
  velocity = 1._f64
  

  N = max(N_x1,N_x2)

  nc_eta1 = N_x1
  nc_eta2 = N_x2


  print *,'#max(N_x1,N_x2)=',N
  
  

  
  SLL_ALLOCATE(f(N_x1+1,N_x2+1),err)
  SLL_ALLOCATE(f_init(N_x1+1,N_x2+1),err)
  SLL_ALLOCATE(f_store(N_x1+1,N_x2+1),err)
  
  SLL_ALLOCATE(a1(N_x1+1,N_x2+1),err)
  SLL_ALLOCATE(a2(N_x1+1,N_x2+1),err)
  SLL_ALLOCATE(psi(N_x1+1,N_x2+1),err)
  
  
  SLL_ALLOCATE(rho(N_x1+1),err)
  SLL_ALLOCATE(E(N_x1+1),err)
  SLL_ALLOCATE(phi_poisson(N_x1+1),err)
  
  SLL_ALLOCATE(buf1d(N+1),err)
  SLL_ALLOCATE(Xstar(N+1),err)
  
  SLL_ALLOCATE(node_positions_x1(N_x1+1),err)
  SLL_ALLOCATE(node_positions_x2(N_x2+1),err)
  
  
  
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

  
  !we suppose uniform mesh for the interpolation
  
  do i=1,N_x1+1
    node_positions_x1(i) = x1_min+(real(i,f64)-1._f64)*delta_x1
  enddo
  
  node_positions_x1 = (node_positions_x1-x1_min)/(x1_max-x1_min)
  
  do i=1,N_x2+1
    node_positions_x2(i) = x2_min+(real(i,f64)-1._f64)*delta_x2
  enddo

  node_positions_x2 = (node_positions_x2-x2_min)/(x2_max-x2_min)
  
  
  !print *,node_positions_x1
  !print *,node_positions_x2
  
  !stop
  

  !construction of mesh
  
   call construct_bgk_mesh(nc_eta1,nc_eta2,mesh_case,&
     &x1n_array,x2n_array,x1c_array,x2c_array,jac_array,integration_points,&
     &geom_x,geom_eta,alpha_mesh,N_x1,N_x2)

  

  spl_per_x1 =>  new_cubic_nonunif_spline_1D( N_x1, SLL_PERIODIC)
  spl_per_x2 =>  new_cubic_nonunif_spline_1D( N_x2, SLL_PERIODIC)

  spl_per_x1_rho =>  new_cubic_nonunif_spline_1D( N_x1, SLL_PERIODIC)

  do i1=1,nc_eta1+1
    do i2=1,nc_eta2+1
      x1 = x1c_array(i1,i2)
      x2 = x2c_array(i1,i2)
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
      f_init(i1,i2) = val*jac_array(i1,i2)      
    enddo
  enddo
  
  
  !compute rho init
  call compute_rho_mapped_mesh2(rho,f_init,integration_points,rho_case,&
    nc_eta1,nc_eta2,geom_eta,jac_array,spl_per_x1_rho)
  !compute E init and psi
  call compute_psi(a1,a2,rho,nc_eta1,nc_eta2,psi,phi_poisson,E,&
  geom_x,x1n_array,x2n_array,x1c_array,x2c_array,jac_array,delta_eta1,delta_eta2,div_case)

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
  
  f=f_init
  
  do step=1,nb_step
    !print *,"#step=",i
    f_store = f
    
    !compute field at time t_n
    call compute_rho_mapped_mesh2(rho,f,integration_points,rho_case,&
      nc_eta1,nc_eta2,geom_eta,jac_array,spl_per_x1_rho)
    call compute_psi(a1,a2,rho,nc_eta1,nc_eta2,psi,phi_poisson,E,&
      geom_x,x1n_array,x2n_array,x1c_array,x2c_array,jac_array,delta_eta1,delta_eta2,div_case)

    val=0._f64
    do i1=1,N_x1
      val = val+E(i1)*E(i1)
    enddo
    val = val/real(N_x1,f64)
    print *,(real(step,f64)-1._f64)*dt,val

    
    ! advect dt/2 with field(t_n)
    call advect_classical_csl(0.5_f64*dt,a1,a2,f,geom_x,N_x1,N_x2,buf1d,&
    node_positions_x1,node_positions_x2,Xstar,spl_per_x1,spl_per_x2)

    !compute field at time t_{n+1/2}
    call compute_rho_mapped_mesh2(rho,f,integration_points,rho_case,&
      nc_eta1,nc_eta2,geom_eta,jac_array,spl_per_x1_rho)
    !compute E and psi
    call compute_psi(a1,a2,rho,nc_eta1,nc_eta2,psi,phi_poisson,E,&
      geom_x,x1n_array,x2n_array,x1c_array,x2c_array,jac_array,delta_eta1,delta_eta2,div_case)

    ! advect dt with field(t_{n+1/2})
    f = f_store
    call advect_classical_csl(dt,a1,a2,f,geom_x,N_x1,N_x2,buf1d,&
    node_positions_x1,node_positions_x2,Xstar,spl_per_x1,spl_per_x2)

    
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


