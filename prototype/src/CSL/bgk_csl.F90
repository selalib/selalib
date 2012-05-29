program bgk_csl
#include "sll_working_precision.h"
#include "sll_mesh_types.h"
#include "sll_memory.h"

  use numeric_constants
  use distribution_function
  use sll_diagnostics
  use sll_csl
  use sll_splines
  use contrib_rho_module
  implicit none
  external compute_translate_nodes_periodic,compute_non_unif_integral2
  !external poisson1dpertrap
  sll_real64 :: x2_min,x2_max,x1_min,x1_max,x1,x2,delta_x1,delta_x2,tmp
  sll_real64 :: mu,xi,L,H
  sll_int    :: i,j,N_phi,err,N_x1,N_x2,i1,i2,N,nb_step
  LOGICAL :: ex
  sll_real64,dimension(:), pointer :: phi,node_positions_x1,node_positions_x2,phi_poisson
  sll_real64,dimension(:), pointer :: new_node_positions,buf_1d,rho,E,rho_exact
  sll_real64,dimension(:,:), pointer :: f,f_init,f_equil
  sll_real64 :: phi_val,delta_x1_phi,xx,dt,alpha,val
  sll_int :: ii,step
  type(cubic_nonunif_spline_1D), pointer :: spl_per_x1, spl_per_x2

  sll_int32 :: nc_eta1, nc_eta2
  sll_real64, dimension(:,:), pointer :: x1c_array, x2c_array, jac_array
  sll_real64, dimension(:,:), pointer :: x1n_array, x2n_array
  sll_real64 :: eta1_min, eta1_max,  eta2_min, eta2_max, eta1, eta2, eta1c, eta2c
  sll_real64 :: delta_eta1, delta_eta2,alpha_mesh
  sll_int  :: mesh_case,ierr,visu_step,test_case,rho_case,phi_case
  type(geometry_2D), pointer :: geom
  type(mesh_descriptor_2D), pointer :: mesh
  type(sll_distribution_function_2D_t), pointer :: dist_func
  character(32), parameter  :: name = 'distribution_function'
  type(field_2D_vec1), pointer :: uniform_field
  type(field_2D_vec1), pointer :: uniform_field_new
  type(field_2D_vec1), pointer :: uniform_field_velocity
  type(csl_workspace), pointer :: csl_work
  sll_real64,dimension(:,:,:), pointer :: integration_points
  sll_real64,dimension(:,:), pointer :: integration_points_val
  character*80,str,str2

  mesh_case = 3
  visu_step = 10
  test_case = 4
  rho_case = 2
  phi_case = 3
     
  alpha_mesh = 1.e-2_f64 !0.1_f64
  
  N_x1 = 64
  N_x2 = 64
  dt = 0.1_f64
  nb_step = 600
  
  N = max(N_x1,N_x2)
  
  print *,'#max(N_x1,N_x2)=',N
  
  
  x2_max = 10._f64
  x2_min = -x2_max
  
  SLL_ALLOCATE(f(N_x1+1,N_x2+1),err)
  SLL_ALLOCATE(f_init(N_x1+1,N_x2+1),err)
  SLL_ALLOCATE(f_equil(N_x1+1,N_x2+1),err)
  SLL_ALLOCATE(node_positions_x1(N_x1+1),err)
  SLL_ALLOCATE(node_positions_x2(N_x2+1),err)
  SLL_ALLOCATE(new_node_positions(N+1),err)
  SLL_ALLOCATE(buf_1d(N+1),err)
  SLL_ALLOCATE(rho(N_x1+1),err)
  SLL_ALLOCATE(rho_exact(N_x1+1),err)
  SLL_ALLOCATE(E(N_x1+1),err)
  SLL_ALLOCATE(phi_poisson(N_x1+1),err)
  SLL_ALLOCATE(integration_points(3,N_x1+1,N_x2+1),err)  
  SLL_ALLOCATE(integration_points_val(2,N_x2),err)  
  
  spl_per_x1 =>  new_cubic_nonunif_spline_1D( N_x1, PERIODIC_SPLINE)
  spl_per_x2 =>  new_cubic_nonunif_spline_1D( N_x2, PERIODIC_SPLINE)
  
  f_equil=0._f64
  
  !physical parameters
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
  
  
  if(test_case>=4)then
    L = 4._f64*sll_pi
    x1_min = 0._f64
    x1_max = L
    x2_min = -6._f64
    x2_max = -x2_min
  endif

!  x1_min = 0._f64
!  x1_max = 1._f64
!  x2_min = 0._f64
!  x2_max = 1._f64
  
  
  
  
  
  
  delta_x1_phi = (x1_max-x1_min)/real(N_phi,f64)
  
  !nb_step = 10
  !dt = L/real(nb_step,f64)
  

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
  
  
  !definition of mesh

  nc_eta1 = N_x1
  nc_eta2 = N_x2

  eta1_min =  0.0_f64
  eta1_max = 1.0_f64 ! 0.15_f64*x1_max! 1.0_f64
  eta2_min =  0.0_f64
  eta2_max =  1.0_f64

  !eta1_min =  -0.5_f64
  !eta1_max = 0.5_f64 ! 0.15_f64*x1_max! 1.0_f64
  !eta2_min =  -0.5_f64
  !eta2_max =  0.5_f64


  !eta1_min =  x1_min
  !eta1_max = x1_max
  !eta2_min = x2_min
  !eta2_max = x2_max



  delta_eta1 = (eta1_max-eta1_min)/real(nc_eta1,f64)
  delta_eta2 = (eta1_max-eta1_min)/real(nc_eta2,f64)
  SLL_ALLOCATE(x1n_array(nc_eta1+1, nc_eta2+1), ierr)
  SLL_ALLOCATE(x2n_array(nc_eta1+1, nc_eta2+1), ierr)
  SLL_ALLOCATE(x1c_array(nc_eta1+1, nc_eta2+1), ierr)
  SLL_ALLOCATE(x2c_array(nc_eta1+1, nc_eta2+1), ierr)
  SLL_ALLOCATE(jac_array(nc_eta1+1, nc_eta2+1), ierr)
  
  
  if(mesh_case==1)then
    do i2=1,nc_eta2+1
      do i1=1,nc_eta1+1
        x1n_array(i1,i2) = x1_min+real(i1-1,f64)*delta_x1
        x2n_array(i1,i2) = x2_min+real(i2-1,f64)*delta_x2
        x1c_array(i1,i2) = x1_min+(real(i1,f64)-0.5_f64)*delta_x1
        x2c_array(i1,i2) = x2_min+(real(i2,f64)-0.5_f64)*delta_x2
        jac_array(i1,i2) = (x1_max-x1_min)*(x2_max-x2_min)
        !jac_array(i1,i2) = 1._f64!(x1_max-x1_min)*(x2_max-x2_min)
      enddo
    enddo
    geom => new_geometry_2D ('from_array',nc_eta1+1,nc_eta2+1, &
       x1n_array, x2n_array, x1c_array, x2c_array, jac_array,PERIODIC,PERIODIC)
    mesh => new_mesh_descriptor_2D(eta1_min, eta1_max, nc_eta1, &
       PERIODIC, eta2_min, eta2_max, nc_eta2, PERIODIC, geom)
       
    dist_func => sll_new_distribution_function_2D(mesh,CELL_CENTERED_DF, name)

    do i2=1,nc_eta2+1
      do i1=1,nc_eta1+1
        !eta1 value of intersecting point (eta2,x1)=constant
        integration_points(1,i1,i2) = (real(i1,f64)-0.5_f64)*(eta1_max-eta1_min)/real(nc_eta1,f64)
        !x2 value of intersecting point (eta2,x1)=constant
        integration_points(2,i1,i2) = x2_min+(real(i2,f64)-0.5_f64)*delta_x2
      enddo
    enddo
    
    
    
  endif
  


  

  if(mesh_case==2)then
     eta2 = 0.0_f64 
     eta2c = eta2 + 0.5_f64*delta_eta2
     do i2= 1, nc_eta2 + 1
        eta1 = 0.0_f64
        eta1c = 0.5_f64*delta_eta1
        do i1 = 1, nc_eta1 + 1
           x1n_array(i1,i2) = eta1 + alpha_mesh * sin(2*sll_pi*eta1) * sin(2*sll_pi*eta2)
           x2n_array(i1,i2) = eta2 + alpha_mesh * sin(2*sll_pi*eta1) * sin(2*sll_pi*eta2)
           x1c_array(i1,i2) = eta1c + alpha_mesh * sin(2*sll_pi*eta1c) * sin(2*sll_pi*eta2c)
           x2c_array(i1,i2) = eta2c + alpha_mesh * sin(2*sll_pi*eta1c) * sin(2*sll_pi*eta2c)
           !x1n_array(i1,i2) = (x1n_array(i1,i2) + alpha_mesh)/(1._f64+2._f64*alpha_mesh)
           !x2n_array(i1,i2) = (x2n_array(i1,i2) + alpha_mesh)/(1._f64+2._f64*alpha_mesh)
           !x1c_array(i1,i2) = (x1c_array(i1,i2) + alpha_mesh)/(1._f64+2._f64*alpha_mesh)
           !x2c_array(i1,i2) = (x2c_array(i1,i2) + alpha_mesh)/(1._f64+2._f64*alpha_mesh)
           jac_array(i1,i2) = (1.0_f64 + alpha_mesh *2._f64 *sll_pi * cos (2*sll_pi*eta1c) * sin (2*sll_pi*eta2c)) * &
             (1.0_f64 + alpha_mesh *2._f64 * sll_pi * sin (2*sll_pi*eta1c) * cos (2*sll_pi*eta2c)) - &
             alpha_mesh *2._f64 *sll_pi * sin (2*sll_pi*eta1c) * cos (2*sll_pi*eta2c) * &
             alpha_mesh *2._f64 * sll_pi * cos (2*sll_pi*eta1c) * sin (2*sll_pi*eta2c)
           eta1 = eta1 + delta_eta1
           eta1c = eta1c + delta_eta1
           x1n_array(i1,i2) = x1_min+x1n_array(i1,i2)*(x1_max-x1_min)
           x2n_array(i1,i2) = x2_min+x2n_array(i1,i2)*(x2_max-x2_min)
           x1c_array(i1,i2) = x1_min+x1c_array(i1,i2)*(x1_max-x1_min)
           x2c_array(i1,i2) = x2_min+x2c_array(i1,i2)*(x2_max-x2_min)
           jac_array(i1,i2) = jac_array(i1,i2)*(x1_max-x1_min)*(x2_max-x2_min)
        end do
        eta2 = eta2 + delta_eta2
        eta2c = eta2c + delta_eta2
     end do

    geom => new_geometry_2D ('from_array',nc_eta1+1,nc_eta2+1, &
       x1n_array, x2n_array, x1c_array, x2c_array, jac_array,PERIODIC,PERIODIC)
    mesh => new_mesh_descriptor_2D(eta1_min, eta1_max, nc_eta1, &
       PERIODIC, eta2_min, eta2_max, nc_eta2, PERIODIC, geom)
    dist_func => sll_new_distribution_function_2D(mesh,CELL_CENTERED_DF, name)

    val = 0._f64
    do i2=1,nc_eta2
      do i1=1,nc_eta1
        x1 = (real(i1,f64)-0.5_f64)/real(nc_eta1,f64)
        eta2 = (real(i2,f64)-0.5_f64)/real(nc_eta2,f64)
        tmp = alpha_mesh*sin(2._f64*sll_pi*eta2)
        do i=1,100
          val = val-(val+tmp*sin(2._f64*sll_pi*val)-x1)/&
          (1._f64+2._f64*sll_pi*tmp*cos(2._f64*sll_pi*val))
        enddo
        if(abs(val+tmp*sin(2._f64*sll_pi*val)-x1)>1.e-14)then
          print *,i1,i2,val+tmp*sin(2._f64*sll_pi*val)-x1,val
          print *,'Problem of convergence of Newton'
          stop
        endif
        !eta1 value of intersecting point (eta2,x1)=constant
        integration_points(1,i1,i2) = val
        !x2 value of intersecting point (eta2,x1)=constant
        integration_points(2,i1,i2) = x2_min+(x1-val+eta2)*(x2_max-x2_min)        
      enddo
    enddo
    !do i1=1,nc_eta1
    !  print *,i1,integration_points(i1,1), integration_points(i1,nc_eta2)
    !enddo
    !stop


  endif



  if(mesh_case==3)then
     eta2 = eta2_min 
     eta2c = eta2_min + 0.5_f64*delta_eta2
     do i2= 1, nc_eta2 + 1
        eta1 = eta1_min
        eta1c = 0.5_f64*delta_eta1
        do i1 = 1, nc_eta1 + 1
           x1n_array(i1,i2) = eta1 + alpha_mesh * sin(2*sll_pi*eta1) * sin(2*sll_pi*eta2)**2
           x2n_array(i1,i2) = eta2 + alpha_mesh * sin(2*sll_pi*eta1) * sin(2*sll_pi*eta2)
           x1c_array(i1,i2) = eta1c + alpha_mesh * sin(2*sll_pi*eta1c) * sin(2*sll_pi*eta2c)**2
           x2c_array(i1,i2) = eta2c + alpha_mesh * sin(2*sll_pi*eta1c)* sin(2*sll_pi*eta2c)
           !x1n_array(i1,i2) = (x1n_array(i1,i2) + alpha_mesh)/(1._f64+2._f64*alpha_mesh)
           !x2n_array(i1,i2) = (x2n_array(i1,i2) + alpha_mesh)/(1._f64+2._f64*alpha_mesh)
           !x1c_array(i1,i2) = (x1c_array(i1,i2) + alpha_mesh)/(1._f64+2._f64*alpha_mesh)
           !x2c_array(i1,i2) = (x2c_array(i1,i2) + alpha_mesh)/(1._f64+2._f64*alpha_mesh)
           !jac_array(i1,i2) = (1.0_f64 + alpha_mesh *2._f64 *sll_pi * cos (2*sll_pi*eta1c) * sin (2*sll_pi*eta2c)) * &
           !  (1.0_f64 + alpha_mesh *2._f64 * sll_pi * sin (2*sll_pi*eta1c) * cos (2*sll_pi*eta2c)) - &
           !  alpha_mesh *2._f64 *sll_pi * sin (2*sll_pi*eta1c) * cos (2*sll_pi*eta2c) * &
           !  alpha_mesh *2._f64 * sll_pi * cos (2*sll_pi*eta1c) * sin (2*sll_pi*eta2c)
           !val =   1.0_f64 + 2._f64*alpha_mesh *sll_pi*sin(2*sll_pi*(eta1c+eta2c))
           !if(abs(jac_array(i1,i2)-val)>1e-13)then
           !  print *,jac_array(i1,i2),val
           !  stop
           !endif
           jac_array(i1,i2) = 1._f64+2._f64*sll_pi*alpha_mesh*sin(2._f64*sll_pi*eta1c)*cos(2._f64*sll_pi*eta2c)&
           +2._f64*sll_pi*alpha_mesh*cos(2._f64*sll_pi*eta1c)&
           -2._f64*sll_pi*alpha_mesh*cos(2._f64*sll_pi*eta1c)*cos(2._f64*sll_pi*eta2c)**2&
           -4._f64*sll_pi**2*alpha_mesh**2*sin(2._f64*sll_pi*eta1c)*cos(2._f64*sll_pi*eta2c)*cos(2._f64*sll_pi*eta1c)&
           +4._f64*sll_pi**2*alpha_mesh**2*sin(2._f64*sll_pi*eta1c)*cos(2._f64*sll_pi*eta2c)**3*cos(2._f64*sll_pi*eta1c)
           eta1 = eta1 + delta_eta1
           eta1c = eta1c + delta_eta1
           x1n_array(i1,i2) = x1_min+x1n_array(i1,i2)*(x1_max-x1_min)
           x2n_array(i1,i2) = x2_min+x2n_array(i1,i2)*(x2_max-x2_min)
           x1c_array(i1,i2) = x1_min+x1c_array(i1,i2)*(x1_max-x1_min)
           x2c_array(i1,i2) = x2_min+x2c_array(i1,i2)*(x2_max-x2_min)
           jac_array(i1,i2) = jac_array(i1,i2)*(x1_max-x1_min)*(x2_max-x2_min)
        end do
        eta2 = eta2 + delta_eta2
        eta2c = eta2c + delta_eta2
     end do

    geom => new_geometry_2D ('from_array',nc_eta1+1,nc_eta2+1, &
       x1n_array, x2n_array, x1c_array, x2c_array, jac_array,PERIODIC,PERIODIC)
    mesh => new_mesh_descriptor_2D(eta1_min, eta1_max, nc_eta1, &
       PERIODIC, eta2_min, eta2_max, nc_eta2, PERIODIC, geom)

    !geom => new_geometry_2D ('from_array',nc_eta1+1,nc_eta2+1, &
    !   x1n_array, x2n_array, x1c_array, x2c_array, jac_array,PERIODIC,COMPACT)
    !mesh => new_mesh_descriptor_2D(eta1_min, eta1_max, nc_eta1, &
    !   PERIODIC, eta2_min, eta2_max, nc_eta2, COMPACT, geom)



    dist_func => sll_new_distribution_function_2D(mesh,CELL_CENTERED_DF, name)

    val = 0._f64
    do i2=1,nc_eta2
      do i1=1,nc_eta1
        x1 = (real(i1,f64)-0.5_f64)/real(nc_eta1,f64)
        eta2 = (real(i2,f64)-0.5_f64)/real(nc_eta2,f64)
        tmp = alpha_mesh*sin(2._f64*sll_pi*eta2)**2
        do i=1,100
          val = val-(val+tmp*sin(2._f64*sll_pi*val)-x1)/&
          (1._f64+2._f64*sll_pi*tmp*cos(2._f64*sll_pi*val))
        enddo
        if(abs(val+tmp*sin(2._f64*sll_pi*val)-x1)>1.e-14)then
          print *,i1,i2,val+tmp*sin(2._f64*sll_pi*val)-x1,val
          print *,'Problem of convergence of Newton'
          stop
        endif
        !eta1 value of intersecting point (eta2,x1)=constant
        integration_points(1,i1,i2) = val
        !x2 value of intersecting point (eta2,x1)=constant
        integration_points(2,i1,i2) = x2_min+((x1-val)/sin(2._f64*sll_pi*eta2)+eta2)*(x2_max-x2_min)        
      enddo
    enddo
    !do i1=1,nc_eta1
    !  print *,i1,integration_points(i1,1), integration_points(i1,nc_eta2)
    !enddo
    !stop


  endif


  open(unit=900,file='intersect_points.dat')  
    do i1=1,N_x1
      x1 = x1_min+(real(i1,f64)-0.5_f64)*delta_x1
      do i2=1,N_x2
        write(900,*) x1,integration_points(2,i1,i2),x1c_array(i1,i2),x2c_array(i1,i2)
      enddo  
    enddo
  close(900)
  
  !stop

  
  !call sll_init_distribution_function_2D( dist_func, GAUSSIAN )
  do i1=1,nc_eta1+1
    do i2=1,nc_eta2+1
      x1 = x1c_array(i1,i2)
      x2 = x2c_array(i1,i2)
      phi_val = 0._f64
      xx = (x1-x1_min)/(x1_max-x1_min)
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
        f_equil(i1,i2) = val*jac_array(i1,i2)
        val = val*(1._f64+0.001_f64*cos(2._f64*sll_pi/L*x1))
      endif
      if(test_case==5)then
        !non linear landau damping
        val = 1._f64/(sqrt(2._f64*sll_pi))*exp(-0.5_f64*x2*x2)
        val = val*(1._f64+0.5_f64*cos(2._f64*sll_pi/L*x1))
      endif
      if(test_case==6.or.test_case==7)then
        !gaussian equilibrium
        val = 1._f64/(sqrt(2._f64*sll_pi))*exp(-0.5_f64*x2*x2)
        !f_equil(i1,i2) = val*jac_array(i1,i2)
      endif

      f_init(i1,i2) = val*jac_array(i1,i2)      
      call sll_set_df_val(dist_func, i1, i2, f_init(i1,i2))      
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
      rho_exact(i1) = mu*(3._f64-2._f64*xi+2*phi_val)/(3._f64-2._f64*xi)*exp(-phi_val)
  enddo
  
  
  call write_mesh_2D(mesh)
  

  call write_distribution_function ( dist_func )


  

  ! initialize CSL  
  csl_work => new_csl_workspace( dist_func )
  uniform_field => new_field_2D_vec1(mesh)
  uniform_field_new => new_field_2D_vec1(mesh)
  uniform_field_velocity => new_field_2D_vec1(mesh)

  do i1=1,nc_eta1+1
    node_positions_x1(i1) = eta1_min+(real(i1,f64)-0.5_f64)*(eta1_max-eta1_min)/real(nc_eta1,f64)
  enddo


  do step = 1, nb_step



    do i2=1,nc_eta2
      do i1 = 1,nc_eta1
        new_node_positions(i1) = integration_points(1,i1,i2)
        if((new_node_positions(i1)>eta1_max).or.(new_node_positions(i1)<eta1_min) )then
          print *,'problem of new_node_position:',new_node_positions(i1),eta1_min,eta1_max
        endif
        if(new_node_positions(i1)<node_positions_x1(1))then
          new_node_positions(i1)=new_node_positions(i1)+eta1_max-eta1_min
        endif      
        buf_1d(i1) = sll_get_df_val(dist_func, i1, i2)/jac_array(i1,i2)!-f_equil(i1,i2)
      enddo
      buf_1d(nc_eta1+1) = buf_1d(1)

      !write(str2,*) i2
      !write(str,*) 'test'//trim(adjustl((str2)))//'.dat'
      
      !open(unit=900,file=trim(adjustl((str))))  
      !do ii=1,nc_eta1
      !  write(900,*) node_positions_x1(ii),buf_1d(ii)
      !enddo  
      !close(900)
      
      
      call compute_spline_nonunif( buf_1d, spl_per_x1, node_positions_x1)
      call interpolate_array_value_nonunif( new_node_positions, buf_1d(1:nc_eta1), nc_eta1, spl_per_x1)
      do i1 = 1,nc_eta1
        integration_points(3,i1,i2) =  buf_1d(i1)
      enddo

      !write(str2,*) i2
      !write(str,*) 'ntest'//trim(adjustl((str2)))//'.dat'
      
      !open(unit=900,file=trim(adjustl((str))))  
      !do ii=1,nc_eta1
      !  write(900,*) new_node_positions(ii),buf_1d(ii)
      !enddo  
      !close(900)
      
      
      
    enddo
    
    !if(step==11)then
    !  stop
    !endif
    !stop
  
    do i1 = 1, nc_eta1
      do i2=1,nc_eta2
        integration_points_val(1,i2) = integration_points(2,i1,i2)
        integration_points_val(2,i2) = integration_points(3,i1,i2)
      enddo
      if(rho_case==1)then
        rho(i1)= compute_non_unif_integral(integration_points_val,nc_eta2)
      endif
      if(rho_case==2)then
        rho(i1)=compute_non_unif_integral_spline(integration_points_val,nc_eta2)
      endif
      if(rho_case==3)then
        rho(i1)=compute_non_unif_integral_gaussian(integration_points_val,nc_eta2)
      endif      
      !if(test_case==4)then
       ! rho(i1) = rho(i1)+1._f64
      !endif  
    enddo  
    E=rho-1._f64
    call poisson1dpertrap(E,x1_max-x1_min,N_x1)
    
    tmp=sum(rho(1:N_x1))*delta_x1
  
    val=0._f64
    do i1=1,N_x1
      val = val+E(i1)*E(i1)
    enddo
    val = val/real(N_x1,f64)
    print *,(real(step,f64)-1._f64)*dt,val,tmp/(x1_max-x1_min)-1._f64

  
    phi_poisson = E
    call poisson1dpertrap(phi_poisson,x1_max-x1_min,N_x1)
    tmp = phi_poisson(1)
    do i1=1,N_x1
      phi_poisson(i1) = -phi_poisson(i1) + tmp
    enddo
    phi_poisson(N_x1+1) = phi_poisson(1) 
    
    if(phi_case==1)then
      do i1=2,N_x1+1
        buf_1d(i1) = 0.5_f64*(phi_poisson(i1)+phi_poisson(i1-1))
      enddo
      buf_1d(1)=buf_1d(N_x1+1)
      phi_poisson(1:N_x1+1) = buf_1d(1:N_x1+1)
    endif
    if(mod(step,visu_step)==0.or. step==1)then
      write(str2,*) step
      write(str,*) 'rho'//trim(adjustl((str2)))//'.dat'
      open(unit=900,file=trim(adjustl((str))))  
        do i1=1,N_x1
          x1 = x1_min+(real(i1,f64)-0.5_f64)*delta_x1
          write(900,*) x1,rho(i1),rho_exact(i1),phi_poisson(i1),E(i1)
        enddo
      close(900)
    endif
   if(test_case==7)then
     phi_poisson = 0._f64
   endif
    
    
!    do i1 = 1, nc_eta1+1 
!      do i2 = 1, nc_eta2+1
!        x1 = x1n_array(i1,i2)
!        x2 = x2n_array(i1,i2)
!        phi_val = 0._f64
!        xx = (x1-x1_min)/(x1_max-x1_min)
!        if(xx<=0._f64)then
!          xx = 0._f64
!        endif
!        if(xx>=1._f64)then
!          xx = xx-1._f64!1._f64-1e-15_f64
!        endif
!        if(xx<=0._f64)then
!          xx = 0._f64
!        endif      
!        xx = xx*real(N_phi,f64)
!        ii = floor(xx)
!        xx = xx-real(ii,f64)      
!        phi_val = (1._f64-xx)*phi(ii+1)+xx*phi(ii+2)     
!        FIELD_2D_AT_I( uniform_field, i1, i2 ) = ( 0.5_f64*x2**2+phi_val)!&
!        !+(x1_max-x1_min)/(real(nb_step,f64)*dt)*x2
!        !-(x2_max-x2_min)/(real(nb_step,f64)*dt)*x1
!      end do
!    end do

    !do i1=1,nc_eta1+1
    !  print *,i1,phi_poisson(i1)
    !enddo
    !stop
    call compute_spline_nonunif( phi_poisson, spl_per_x1, node_positions_x1)
    
    do i1 = 1, nc_eta1+1 
      do i2 = 1, nc_eta2+1
        x1 = x1n_array(i1,i2)
        x2 = x2n_array(i1,i2)
        phi_val = 0._f64
        if(phi_case==1)then
          xx = (x1-x1_min)/(x1_max-x1_min)
          if(xx<=0._f64)then
            xx = xx+1._f64
          endif
          if(xx>=1._f64)then
            xx = xx-1._f64!1._f64-1e-15_f64
          endif
          if(xx<=0._f64)then
            xx = 0._f64
          endif      
          xx = xx*real(N_x1,f64)
          ii = floor(xx)
          xx = xx-real(ii,f64)      
          phi_val = (1._f64-xx)*phi_poisson(ii+1)+xx*phi_poisson(ii+2)     
        endif
        if(phi_case==2)then
          xx = (x1-x1_min)/(x1_max-x1_min)-0.5_f64/real(N_x1,f64)
          if(xx<=0._f64)then
            xx = xx+1._f64
          endif
          if(xx>=1._f64)then
            xx = xx-1._f64!1._f64-1e-15_f64
          endif
          if(xx<=0._f64)then
            xx = 0._f64
          endif      
          xx = xx*real(N_x1,f64)
          ii = floor(xx)
          xx = xx-real(ii,f64)      
          phi_val = (1._f64-xx)*phi_poisson(ii+1)+xx*phi_poisson(ii+2)
          tmp = phi_val     
        endif
        if(phi_case==3)then
          xx = (x1-x1_min)/(x1_max-x1_min)!+0.5_f64/real(N_x1,f64)
          if(xx<0.5_f64/real(N_x1,f64))then
            xx = xx+1._f64
          endif
          new_node_positions(1)=xx
          call interpolate_array_value_nonunif( new_node_positions(1:1), buf_1d(1:1), 1, spl_per_x1)
          phi_val = buf_1d(1)
          !phi_val = interpolate_value_nonunif( xx, spl_per_x1 )     
        endif
        !print *,x1,xx,phi_val
        FIELD_2D_AT_I( uniform_field, i1, i2 ) = ( 0.5_f64*x2**2+phi_val)!&
        !FIELD_2D_AT_I( uniform_field_velocity, i1, i2 ) = 0.5_f64*x2**2!&
        !+(x1_max-x1_min)/(real(nb_step,f64)*dt)*x2
        !-(x2_max-x2_min)/(real(nb_step,f64)*dt)*x1
      end do
      !if(step<=2)then
      !  print *,step,1,i1,tmp,phi_val,tmp-phi_val
      !endif  
      
    end do
    
    do i1 = 1, nc_eta1 
      do i2 = 1, nc_eta2    
        f(i1,i2) = sll_get_df_val(dist_func, i1, i2)
      enddo  
    enddo
    
    
    call csl_second_order(csl_work, dist_func, uniform_field, uniform_field, dt)
    do i1=1,nc_eta1+1
      node_positions_x1(i1) = eta1_min+(real(i1,f64)-0.5_f64)*(eta1_max-eta1_min)/real(nc_eta1,f64)
    enddo
    do i2=1,nc_eta2
      do i1 = 1,nc_eta1
        new_node_positions(i1) = integration_points(1,i1,i2)
        if((new_node_positions(i1)>eta1_max).or.(new_node_positions(i1)<eta1_min) )then
          print *,'problem of new_node_position:',new_node_positions(i1),eta1_min,eta1_max
          stop
        endif
        if(new_node_positions(i1)<node_positions_x1(1))then
          new_node_positions(i1)=new_node_positions(i1)+eta1_max-eta1_min
        endif      
        buf_1d(i1) = sll_get_df_val(dist_func, i1, i2)/jac_array(i1,i2)!-f_equil(i1,i2)
      enddo
      buf_1d(nc_eta1+1) = buf_1d(1)
      call compute_spline_nonunif( buf_1d, spl_per_x1, node_positions_x1)
      call interpolate_array_value_nonunif( new_node_positions, buf_1d(1:nc_eta1), nc_eta1, spl_per_x1)
      do i1 = 1,nc_eta1
        integration_points(3,i1,i2) =  buf_1d(i1)
      enddo
    enddo
  
    do i1 = 1, nc_eta1
      do i2=1,nc_eta2
        integration_points_val(1,i2) = integration_points(2,i1,i2)
        integration_points_val(2,i2) = integration_points(3,i1,i2)
      enddo
      !rho(i1)= compute_non_unif_integral(integration_points_val,nc_eta2)
      if(rho_case==1)then
        rho(i1)= compute_non_unif_integral(integration_points_val,nc_eta2)
      endif
      if(rho_case==2)then
        rho(i1)=compute_non_unif_integral_spline(integration_points_val,nc_eta2)
      endif
      if(rho_case==3)then
        rho(i1)=compute_non_unif_integral_gaussian(integration_points_val,nc_eta2)
      endif      
      if(rho_case==4)then
        rho(i1)=compute_non_unif_integral_gaussian_sym(integration_points_val,nc_eta2)
      endif      
      !if(test_case==4)then      
      !  rho(i1) = rho(i1)+1._f64
      !endif  
    enddo
    E=rho-1._f64
    call poisson1dpertrap(E,x1_max-x1_min,N_x1)  
    phi_poisson = E
    call poisson1dpertrap(phi_poisson,x1_max-x1_min,N_x1)
    tmp = phi_poisson(1)
    do i1=1,N_x1
      phi_poisson(i1) = -phi_poisson(i1) + tmp
    enddo
    phi_poisson(N_x1+1) = phi_poisson(1) 
    if(phi_case==1)then
      do i1=2,N_x1+1
        buf_1d(i1) = 0.5_f64*(phi_poisson(i1)+phi_poisson(i1-1))
      enddo
      buf_1d(1)=buf_1d(N_x1+1)
      phi_poisson(1:N_x1+1) = buf_1d(1:N_x1+1)
    endif
  
  
   if(test_case==7)then
     phi_poisson = 0._f64
   endif

    call compute_spline_nonunif( phi_poisson, spl_per_x1, node_positions_x1)
    

    do i1 = 1, nc_eta1+1 
      do i2 = 1, nc_eta2+1
        x1 = x1n_array(i1,i2)
        x2 = x2n_array(i1,i2)
        phi_val = 0._f64
        if(phi_case==1)then
          xx = (x1-x1_min)/(x1_max-x1_min)
          if(xx<=0._f64)then
            xx = xx +1._f64
          endif
          if(xx>=1._f64)then
            xx = xx-1._f64!1._f64-1e-15_f64
          endif
          if(xx<=0._f64)then
            xx = 0._f64
          endif      
          xx = xx*real(N_x1,f64)
          ii = floor(xx)
          xx = xx-real(ii,f64)      
          phi_val = (1._f64-xx)*phi_poisson(ii+1)+xx*phi_poisson(ii+2)     
          !&
          !+(x1_max-x1_min)/(real(nb_step,f64)*dt)*x2
          !-(x2_max-x2_min)/(real(nb_step,f64)*dt)*x1
        endif  
        if(phi_case==2)then
          xx = (x1-x1_min)/(x1_max-x1_min)-0.5_f64/real(N_x1,f64)
          if(xx<=0._f64)then
            xx = xx+1._f64
          endif
          if(xx>=1._f64)then
            xx = xx-1._f64!1._f64-1e-15_f64
          endif
          if(xx<=0._f64)then
            xx = 0._f64
          endif      
          xx = xx*real(N_x1,f64)
          ii = floor(xx)
          xx = xx-real(ii,f64)      
          phi_val = (1._f64-xx)*phi_poisson(ii+1)+xx*phi_poisson(ii+2)     
          !&
          !+(x1_max-x1_min)/(real(nb_step,f64)*dt)*x2
          !-(x2_max-x2_min)/(real(nb_step,f64)*dt)*x1
          tmp = phi_val
        endif  
        if(phi_case==3)then
          xx = (x1-x1_min)/(x1_max-x1_min)!+0.5_f64/real(N_x1,f64)
          if(xx<0.5_f64/real(N_x1,f64))then
            xx = xx+1._f64
          endif
          new_node_positions(1)=xx
          call interpolate_array_value_nonunif( new_node_positions(1:1), buf_1d(1:1), 1, spl_per_x1)
          phi_val = buf_1d(1)
         
          !phi_val = interpolate_value_nonunif( xx, spl_per_x1 )     
        endif
        FIELD_2D_AT_I( uniform_field_new, i1, i2 ) = ( 0.5_f64*x2**2+phi_val)

      end do
      !if(step<=2)then
      !  print *,step,2,i1,tmp,phi_val,tmp-phi_val
      !endif  
    end do
    !if(step==2)then
    !stop
    !endif    
    
    do i1 = 1, nc_eta1 
      do i2 = 1, nc_eta2    
        call sll_set_df_val(dist_func, i1, i2,f(i1,i2))
      enddo  
    enddo

    !do i2 = 1, nc_eta2
    !  val = 0._f64
    !  do i1 = 1, nc_eta1    
    !     val=val+sll_get_df_val(dist_func, i1, i2)
    !  enddo
    !  val = val/real(nc_eta1,f64)
    !  buf_1d(i2)=0._f64!val
    !enddo


    !do i1 = 1, nc_eta1 
    !  do i2 = 1, nc_eta2
    !    val = sll_get_df_val(dist_func, i1, i2)
    !    call sll_set_df_val(dist_func, i1, i2,val-buf_1d(i2))
    !  enddo  
    !enddo


    !do i1 = 1, nc_eta1 
    !  do i2 = 1, nc_eta2    
    !    call sll_set_df_val(dist_func, i1, i2,f(i1,i2))
    !  enddo  
    !enddo


    call csl_second_order(csl_work, dist_func, uniform_field, uniform_field_new, dt)
    !call csl_second_order(csl_work, dist_func, uniform_field_velocity, uniform_field_velocity, dt)


    !do i1 = 1, nc_eta1 
    !  do i2 = 1, nc_eta2
    !    val = sll_get_df_val(dist_func, i1, i2)
    !    !call sll_set_df_val(dist_func, i1, i2,val+f_equil(i1,i2))
    !    call sll_set_df_val(dist_func, i1, i2,val+buf_1d(i2))
    !  enddo  
    !enddo

    !do i1 = 1, nc_eta1+1 
    !  do i2 = 1, nc_eta2+1
    !    val = sll_get_df_val(dist_func, i1, i2)
    !    call sll_set_df_val(dist_func, i1, i2,val-buf_1d(i2))
    !    FIELD_2D_AT_I( uniform_field, i1, i2 ) = FIELD_2D_AT_I( uniform_field, i1, i2 )&
    !    -FIELD_2D_AT_I( uniform_field_velocity, i1, i2 )
    !    FIELD_2D_AT_I( uniform_field_new, i1, i2 ) = FIELD_2D_AT_I( uniform_field_new, i1, i2 )&
    !   -FIELD_2D_AT_I( uniform_field_velocity, i1, i2 )        
    !  enddo  
    !enddo

    !do i1 = 1, nc_eta1+1 
    !  do i2 = 1, nc_eta2+1
    !    val = sll_get_df_val(dist_func, i1, i2)        
    !    call sll_set_df_val(dist_func, i1, i2,f(i1,i2))
    !    f(i1,i2)=val
    !  enddo  
    !enddo

    
    !call csl_second_order(csl_work, dist_func, uniform_field, uniform_field_new, dt)
    
    !do i1 = 1, nc_eta1+1 
    !  do i2 = 1, nc_eta2+1
    !    val = sll_get_df_val(dist_func, i1, i2)        
    !    call sll_set_df_val(dist_func, i1, i2,f(i1,i2)+val)
    !  enddo  
    !enddo

    
    
    if(mod(step,visu_step)==0)then
      call write_distribution_function ( dist_func )
    endif
  end do

  val = 0._f64
  do i1=1,nc_eta1+1
    do i2=1,nc_eta2+1
      val = max(abs(sll_get_df_val(dist_func, i1, i2)/jac_array(i1,i2)&
      -f_init(i1,i2)/jac_array(i1,i2)),val)
    enddo
  enddo  
  
  print *,'#',nc_eta1,nc_eta2,dt,nb_step,val
  
  open(unit=900,file='field_final.dat')  
    do i1=1,N_x1+1
      x1 = x1_min+real(i1-1,f64)*delta_x1
      write(900,*) x1,E(i1),rho(i1)
    enddo
  close(900)


  
end program


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


subroutine compute_non_unif_integral2(x_points,f_points,N_points,val)
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
  
  
end  subroutine compute_non_unif_integral2


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

