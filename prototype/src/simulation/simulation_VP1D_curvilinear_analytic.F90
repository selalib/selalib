module simulation_VP1D_curvilinear_analytic

!!!!!!!!!!!!!!!!!!!!!!!
!  Vlasov-Poisson 1D simulation
!  the mesh is curvilinear with analytic mapping (Colella mesh)


#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"

#ifndef STDF95
  use sll_simulation_base
#endif
  use cubic_non_uniform_splines
  use sll_constants
  implicit none

#ifdef STDF95
  type :: sll_simulation_VP1D_curvilinear_analytic
#else
  type, extends(sll_simulation_base_class) :: &
    sll_simulation_VP1D_curvilinear_analytic
#endif
    ! Numerical parameters
    sll_real64 :: dt
#ifndef STDF95
  contains
    procedure, pass(sim) :: run => run_VP1D_curvilinear_analytic
    procedure, pass(sim) :: init_from_file => VP1D_curvilinear_analytic_init
#endif
  end type sll_simulation_VP1D_curvilinear_analytic

contains

  subroutine VP1D_curvilinear_analytic_init(sim, filename)
#ifdef STDF95
    type(sll_simulation_VP1D_curvilinear_analytic), intent(inout)  :: sim
#else
    class(sll_simulation_VP1D_curvilinear_analytic), intent(inout) :: sim
#endif
    character(len=*), intent(in)                                 :: filename
    ! Declare here the variables to be read in through a namelist and that
    ! are to be kept inside the sim object. Look at the parallel vp4d simulation
    ! for an example.
    print *, 'This is a dummy function. Needs implementation.'
    
    print *,filename
    print *,sim%dt
    stop
    
  end subroutine VP1D_curvilinear_analytic_init

  ! Note that the following function has no local variables, which is silly...
  ! This just happened since the guts of the unit test were transplanted here
  ! directly, but this should be cleaned up.
  subroutine run_VP1D_curvilinear_analytic(sim)
    implicit none
#ifdef STDF95
    type(sll_simulation_VP1D_curvilinear_analytic), intent(inout) :: sim
#else
    class(sll_simulation_VP1D_curvilinear_analytic), intent(inout) :: sim
#endif
    type(cubic_nonunif_spline_1D), pointer :: spl_per_x1,spl_per_x2
    sll_int32 :: nc_eta1,nc_eta2,N_x1_poisson,N,nb_step
    sll_int32:: mesh_case,test_case,rho_case,div_Case
    sll_int32:: i,i1,i2,err,step
    sll_real64 :: x1_min,x1_max,x2_min,x2_max
    sll_real64 :: eta1_min,eta1_max,eta2_min,eta2_max
    sll_real64 :: delta_x1,delta_x2,dt,delta_eta1,delta_eta2
    sll_real64,dimension(:), pointer :: node_positions_x1,node_positions_x2
    !sll_real64,dimension(:), pointer :: node_positions_x1_unit,node_positions_x2_unit
    !sll_real64,dimension(:), pointer :: node_positions_x1_poisson
    sll_real64,dimension(:), pointer :: rho,E,phi_poisson
    sll_real64,dimension(:), pointer :: buf1d,Xstar
    !sll_real64,dimension(:), pointer :: random_vector_x1,random_vector_x2
    sll_real64,dimension(:,:), pointer :: f,f_store
    !sll_real64,dimension(:,:), pointer :: integration_points
    sll_real64,dimension(:,:),pointer::x1n_array,x2n_array,x1c_array,x2c_array
    sll_real64,dimension(:,:),pointer::jac_array,a1,a2,psi
    !sll_real64, dimension(:,:), pointer :: a1,a2,psi
    sll_real64,dimension(:,:,:),pointer::integration_points
    sll_real64 :: landau_alpha,alpha_mesh,geom_x(2,2),geom_eta(2,2)
    sll_real64 :: val,x1,x2
    
    !initialization of parameters
    nc_eta1=64
    nc_eta2=nc_eta1
    N_x1_poisson=1024
    mesh_case=2 !
    test_case=1 !1: landau damping
    rho_case=1  !1: trap 2: splines 3: gaussian 4: gaussian sym 
                !5: spline_per
    div_case=1 !            
                
    landau_alpha=1e-3_f64
    alpha_mesh=0._f64
    dt=0.1_f64
    nb_step=601
    
    sim%dt = dt
    
    N=max(nc_eta1,nc_eta2)
    
    print *,'#nc_eta1=',nc_eta1,'#nc_eta2=',nc_eta2
    print *,'#mesh_case=',mesh_case
    print *,'#test_case=',test_case
    print *,'#rho_case=',rho_case
    print *,'#landau_alpha=',landau_alpha
    print *,'#alpha_mesh=',alpha_mesh
    print *,'#N_x1_poisson=',N_x1_poisson
    print *,'#nb_step=',nb_step
    
    x1_min = 0._f64
    x1_max = 4._f64*sll_pi
    x2_min = -6._f64
    x2_max = -x2_min
    

    eta1_min =  0.0_f64
    eta1_max =  1.0_f64
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


    delta_x1 = (x1_max-x1_min)/real(nc_eta1,f64)
    delta_x2 = (x2_max-x2_min)/real(nc_eta2,f64)

    delta_eta1 = (eta1_max-eta1_min)/real(nc_eta1,f64)
    delta_eta2 = (eta2_max-eta2_min)/real(nc_eta2,f64)
    
    spl_per_x1 =>  new_cubic_nonunif_spline_1D( nc_eta1, SLL_PERIODIC)
    spl_per_x2 =>  new_cubic_nonunif_spline_1D( nc_eta2, SLL_PERIODIC)
    
    !initialization of the curvilinear mesh
    
    call construct_bgk_mesh(nc_eta1,nc_eta2,mesh_case,&
     &x1n_array,x2n_array,x1c_array,x2c_array,jac_array,integration_points,&
     &geom_x,geom_eta,alpha_mesh)

    SLL_ALLOCATE(f(nc_eta1+1,nc_eta2+1),err)
    SLL_ALLOCATE(f_store(nc_eta1+1,nc_eta2+1),err)
    SLL_ALLOCATE(a1(nc_eta1+1,nc_eta2+1),err)
    SLL_ALLOCATE(a2(nc_eta1+1,nc_eta2+1),err)
    SLL_ALLOCATE(psi(nc_eta1+1,nc_eta2+1),err)
  
  
    SLL_ALLOCATE(rho(nc_eta1+1),err)
    SLL_ALLOCATE(E(nc_eta1+1),err)
    SLL_ALLOCATE(phi_poisson(nc_eta1+1),err)
    SLL_ALLOCATE(Xstar(N+1),err)
    SLL_ALLOCATE(buf1d(N+1),err)

  
    SLL_ALLOCATE(node_positions_x1(nc_eta1+1),err)
    SLL_ALLOCATE(node_positions_x2(nc_eta2+1),err)

     
    do i=1,nc_eta1+1
      node_positions_x1(i) = (real(i,f64)-1._f64)/real(nc_eta1,f64)
    enddo

    do i=1,nc_eta2+1
      node_positions_x2(i) = (real(i,f64)-1._f64)/real(nc_eta2,f64)
    enddo

    
    !initialization of distribution function
    do i1=1,nc_eta1
      do i2=1,nc_eta2
        x1 = x1c_array(i1,i2)
        x2 = x2c_array(i1,i2)
        !Landau damping
        if(test_case==1)then
          val = 1._f64/(sqrt(2._f64*sll_pi))*exp(-0.5_f64*x2*x2)
          val = val*(1._f64+landau_alpha*cos(2._f64*sll_pi/(x1_max-x1_min)*x1))
        endif
        f(i1,i2) = val*jac_array(i1,i2)      
     enddo
   enddo    


  !compute rho
  call compute_rho_mapped_mesh2(rho,f,integration_points,rho_case,&
    nc_eta1,nc_eta2,geom_eta,jac_array,spl_per_x1)
  !compute E and psi
  call compute_psi(a1,a2,rho,nc_eta1,nc_eta2,psi,phi_poisson,E,&
  geom_x,x1n_array,x2n_array,x1c_array,x2c_array,jac_array,delta_eta1,delta_eta2,div_case)
     

  do step=1,nb_step
    !print *,"#step=",i
    f_store = f
    
    !compute field at time t_n
    call compute_rho_mapped_mesh2(rho,f,integration_points,rho_case,&
      nc_eta1,nc_eta2,geom_eta,jac_array,spl_per_x1)
    call compute_psi(a1,a2,rho,nc_eta1,nc_eta2,psi,phi_poisson,E,&
      geom_x,x1n_array,x2n_array,x1c_array,x2c_array,jac_array,delta_eta1,delta_eta2,div_case)

    val=0._f64
    do i1=1,nc_eta1
      val = val+E(i1)*E(i1)
    enddo
    val = val/real(nc_eta1,f64)
    print *,(real(step,f64)-1._f64)*dt,val

    
    ! advect dt/2 with field(t_n)
    call advect_classical_csl(0.5_f64*dt,a1,a2,f,geom_x,nc_eta1,nc_eta2,buf1d,&
    node_positions_x1,node_positions_x2,Xstar,spl_per_x1,spl_per_x2)

    !compute field at time t_{n+1/2}
    call compute_rho_mapped_mesh2(rho,f,integration_points,rho_case,&
      nc_eta1,nc_eta2,geom_eta,jac_array,spl_per_x1)
    !compute E and psi
    call compute_psi(a1,a2,rho,nc_eta1,nc_eta2,psi,phi_poisson,E,&
      geom_x,x1n_array,x2n_array,x1c_array,x2c_array,jac_array,delta_eta1,delta_eta2,div_case)

    ! advect dt with field(t_{n+1/2})
    f = f_store
    call advect_classical_csl(dt,a1,a2,f,geom_x,nc_eta1,nc_eta2,buf1d,&
    node_positions_x1,node_positions_x2,Xstar,spl_per_x1,spl_per_x2)
   enddo
    
    
    
    
  end subroutine run_VP1D_curvilinear_analytic

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



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




!!!!!!!!!!!!!!!!!!!!!!

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  functions for computing rho
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine compute_rho_mapped_mesh2&
  (rho,f,integration_points,rho_case,nc_eta1,nc_eta2,geom,jac_array,spl_per_x1)
  sll_int :: rho_case
  sll_real64,dimension(:,:,:),pointer :: integration_points
  sll_real64,dimension(:,:),pointer :: integration_points_val
  sll_real64,dimension(:), pointer :: rho
  sll_real64,dimension(:), pointer :: node_positions_x1
  sll_real64,dimension(:), pointer :: new_node_positions
  sll_real64,dimension(:), pointer :: buf_1d
  sll_real64,dimension(:,:), pointer :: jac_array
  sll_real64,dimension(:,:), pointer :: f
  sll_real64,intent(in) :: geom(2,2)
  sll_int,intent(in) :: nc_eta1,nc_eta2
  sll_real64 :: eta1_min,eta1_max,eta2_min,eta2_max
  sll_int :: i1,i2
  type(cubic_nonunif_spline_1D), pointer :: spl_per_x1
  sll_int :: err

  eta1_min = geom(1,1)
  eta1_max = geom(2,1)
  eta2_min = geom(1,2)
  eta2_max = geom(2,2)
  
  SLL_ALLOCATE(node_positions_x1(nc_eta1+1),err)
  SLL_ALLOCATE(new_node_positions(nc_eta1),err)
  SLL_ALLOCATE(buf_1d(nc_eta1+1),err)
  SLL_ALLOCATE(integration_points_val(2,nc_eta1+1),err)
  do i1=1,nc_eta1+1
    node_positions_x1(i1) = eta1_min+(real(i1,f64)-0.5_f64)*(eta1_max-eta1_min)/real(nc_eta1,f64)
    !node_positions_x1(i1) = eta1_min+(real(i1,f64)-1._f64)*(eta1_max-eta1_min)/real(nc_eta1,f64)
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
        !buf_1d(i1) = sll_get_df_val(dist_func, i1, i2)/jac_array(i1,i2)!-f_equil(i1,i2)
        buf_1d(i1) = f(i1,i2)/jac_array(i1,i2)
        
      enddo
      buf_1d(nc_eta1+1) = buf_1d(1)
      call compute_spline_nonunif( buf_1d, spl_per_x1, node_positions_x1)
      call interpolate_array_value_nonunif( new_node_positions, buf_1d(1:nc_eta1), nc_eta1, spl_per_x1)
      do i1 = 1,nc_eta1
        integration_points(3,i1,i2) =  buf_1d(i1)
      enddo
    enddo
   
   

  !Compute rho0 and E0
    do i1 = 1, nc_eta1
      do i2=1,nc_eta2
        integration_points_val(1,i2) = integration_points(2,i1,i2)
        integration_points_val(2,i2) = integration_points(3,i1,i2)
      enddo
      rho(i1)= compute_non_unif_integral(integration_points_val,nc_eta2,rho_case)
!      if(rho_case==1)then
!        rho(i1)= compute_non_unif_integral(integration_points_val,nc_eta2)
!      endif
!      if(rho_case==2)then
!        rho(i1)=compute_non_unif_integral_spline(integration_points_val,nc_eta2)
!      endif
!      if(rho_case==3)then
!        rho(i1)=compute_non_unif_integral_gaussian(integration_points_val,nc_eta2)
!      endif      
!      if(rho_case==4)then
!        rho(i1)=compute_non_unif_integral_gaussian_sym(integration_points_val,nc_eta2)
!      endif      
!      !if(test_case==4)then      
!      !  rho(i1) = rho(i1)+1._f64
!      !endif  
   enddo
   
   rho(nc_eta1+1)=rho(1)
 
     
  SLL_DEALLOCATE(node_positions_x1,err)
  SLL_DEALLOCATE(new_node_positions,err)
  SLL_DEALLOCATE(buf_1d,err)
  SLL_DEALLOCATE(integration_points_val,err)



end subroutine compute_rho_mapped_mesh2



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  functions for computing psi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


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
  
  
  if(x1n_array(1,1)==1e-23)then
    print *,'#is that possible?'
  endif
  if(x2n_array(1,1)==1e-23)then
    print *,'#is that possible?'
  endif


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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! functions for making the advection
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine advect_classical_csl_1(dt,a1,f,geom_x,N_x1,N_x2,buf,node_positions_x1,&
Xstar,spl_per_x1)
!solve \partial_t f(t,x,y) + \partial_x(a1(x,y)f(t,x,y))=0 over dt
  sll_int32,intent(in) :: N_x1,N_x2
  sll_real64,intent(in) :: dt,geom_x(2,2)
  sll_real64,dimension(:),pointer :: buf,node_positions_x1,Xstar
  sll_real64,dimension(:,:),pointer:: a1,f
  type(cubic_nonunif_spline_1D), pointer :: spl_per_x1
  sll_int32 :: i1,i2,i1m1

  if(geom_x(1,1)==1e-23)then
    print *,'#is that possible?'
  endif


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
  
  if(geom_x(1,1)==1e-23)then
    print *,'#is that possible?'
  endif
  
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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! function for the construction of the mesh
!!!!!!!!!!!!!!!!!!!!!!!
 
 
   subroutine construct_bgk_mesh(nc_eta1,nc_eta2,mesh_case,&
   &x1n_array,x2n_array,x1c_array,x2c_array,jac_array,integration_points,&
   &geom_x,geom_eta,alpha_mesh)
    use sll_constants
    implicit none
    sll_int32,intent(in)::nc_eta1,nc_eta2,mesh_case
    sll_real64,intent(in)::geom_x(2,2),geom_eta(2,2),alpha_mesh
    sll_real64,dimension(:,:),pointer::x1n_array,x2n_array,x1c_array,x2c_array
    sll_real64,dimension(:,:),pointer::jac_array
    sll_real64,dimension(:,:,:),pointer::integration_points
    sll_int32  :: i1,i2,err,i
    sll_real64 :: x1_min,x1_max,x2_min,x2_max,delta_x1,delta_x2,x1
    sll_real64 :: eta1_min,eta1_max,eta2_min,eta2_max,delta_eta1,delta_eta2,eta1,eta1c,eta2,eta2c
    sll_real64 :: val,tmp
    sll_real64 :: slope_mesh1,slope_mesh2,wk,ll,dxx,slope_mesh3
    sll_int    ::Nzon,Nzon2!,Nzon3
    sll_real64 , dimension(4)         :: ws
    sll_int , dimension(4)         :: Ns

    
    

    SLL_ALLOCATE(x1n_array(nc_eta1+1, nc_eta2+1), err)
    SLL_ALLOCATE(x2n_array(nc_eta1+1, nc_eta2+1), err)
    SLL_ALLOCATE(x1c_array(nc_eta1+1, nc_eta2+1), err)
    SLL_ALLOCATE(x2c_array(nc_eta1+1, nc_eta2+1), err)
    SLL_ALLOCATE(jac_array(nc_eta1+1, nc_eta2+1), err)
    SLL_ALLOCATE(integration_points(3,nc_eta1+1,nc_eta2+1),err)

    
    x1_min=geom_x(1,1)
    x1_max=geom_x(2,1)
    x2_min=geom_x(1,2)
    x2_max=geom_x(2,2)

    eta1_min=geom_eta(1,1)
    eta1_max=geom_eta(2,1)
    eta2_min=geom_eta(1,2)
    eta2_max=geom_eta(2,2)

        


    delta_x1 = (x1_max-x1_min)/real(nc_eta1,f64)
    delta_x2 = (x2_max-x2_min)/real(nc_eta2,f64)
    
    delta_eta1 = (eta1_max-eta1_min)/real(nc_eta1,f64)
    delta_eta2 = (eta1_max-eta1_min)/real(nc_eta2,f64)
    
    
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
    !geom => new_geometry_2D ('from_array',nc_eta1+1,nc_eta2+1, &
    !   x1n_array, x2n_array, x1c_array, x2c_array, jac_array,PERIODIC,PERIODIC)
    !mesh => new_mesh_descriptor_2D(eta1_min, eta1_max, nc_eta1, &
    !   PERIODIC, eta2_min, eta2_max, nc_eta2, PERIODIC, geom)
       
    !dist_func => sll_new_distribution_function_2D(mesh,CELL_CENTERED_DF, name)

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

    !geom => new_geometry_2D ('from_array',nc_eta1+1,nc_eta2+1, &
    !  x1n_array, x2n_array, x1c_array, x2c_array, jac_array,PERIODIC,PERIODIC)
    !mesh => new_mesh_descriptor_2D(eta1_min, eta1_max, nc_eta1, &
    !  PERIODIC, eta2_min, eta2_max, nc_eta2, PERIODIC, geom)
    !dist_func => sll_new_distribution_function_2D(mesh,CELL_CENTERED_DF, name)

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
  
        
  open(unit=900,file='xn_array.dat')  
    do i1=1,nc_eta1
      !x1 = x1_min+(real(i1,f64)-0.5_f64)*delta_x1
      x1 = x1_min+(real(i1,f64)-1._f64)*delta_x1
      do i2=1,nc_eta2
        !write(900,*) x1,integration_points(2,i1,i2),x1c_array(i1,i2),x2c_array(i1,i2)
        write(900,*) x1n_array(i1,i2),x2n_array(i1,i2)
      enddo  
    enddo
  close(900)
       
      

      !geom => new_geometry_2D ('from_array',nc_eta1+1,nc_eta2+1, &
      !   x1n_array, x2n_array, x1c_array, x2c_array, jac_array,PERIODIC,PERIODIC)
      !mesh => new_mesh_descriptor_2D(eta1_min, eta1_max, nc_eta1, &
      !   PERIODIC, eta2_min, eta2_max, nc_eta2, PERIODIC, geom)

      !geom => new_geometry_2D ('from_array',nc_eta1+1,nc_eta2+1, &
      !   x1n_array, x2n_array, x1c_array, x2c_array, jac_array,PERIODIC,COMPACT)
      !mesh => new_mesh_descriptor_2D(eta1_min, eta1_max, nc_eta1, &
      !   PERIODIC, eta2_min, eta2_max, nc_eta2, COMPACT, geom)



      !dist_func => sll_new_distribution_function_2D(mesh,CELL_CENTERED_DF, name)

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
  

   if(mesh_case==10)then
   ! 3 parts 
      slope_mesh1=1._f64/2._f64
      wk=2.82_f64
      ll=0.5_f64
      !stop
      Ns(1)=int((-wk-ll-x2_min)/delta_x2)+1
      ws(1)=x2_min+(Ns(1)-1)*delta_x2
      !Ns(2)=int((-wk+ll-x2_min)/delta_x2)+1
      !ws(2)=x2_min+(Ns(2)-1)*delta_x2
      
     ! Ns(3)=int((wk-ll-x2_min)/delta_x2)
      !ws(3)=x2_min+(Ns(3)-1)*delta_x2
      Ns(4)=int((wk+ll-x2_min)/delta_x2)+1
      ws(4)=x2_max-(ws(1)-x2_min)!x2_min+(Ns(4)-1)*delta_x2
     !iiw2=int((-wk-1-x2_min)/delta_x2)
    ! print*,ws(1)-ws(2),ws(3)-ws(4),delta_x2
     !stop
     !stop
     Nzon= int((ws(4)-ws(1))/(slope_mesh1*delta_x2))
     if(mod(nc_eta2-Nzon,2)==0) then
      Nzon=Nzon
     else
      Nzon=Nzon+1
     endif 
     ! print*,"slop",slope_mesh1,Nzon,N_x2-2*Nzon
      !stop
      !if(
      !Nzon= int(1._f64/(slope_mesh1*delta_x2))
      !slope_mesh2=(x2_max-x2_min-(ws(4)-ws(1)))/((N_x2-Nzon)*delta_x2)
        slope_mesh2=(ws(1)-x2_min)/((nc_eta2-Nzon)*delta_x2/2._f64)
        slope_mesh1=(ws(4)-ws(1))/(Nzon*delta_x2)
       !slope_mesh2=(x2_max-x2_min-2)/((N_x2-2*Nzon)*delta_x2)
       !print*,(ws(4)-ws(1))/slope_mesh1+2*(ws(1)-x2_min)/slope_mesh2,(x2_max-x2_min)!, (ws(1)-x2_min)/((N_x2-Nzon)*delta_x2/2._f64),Nzon,N_x2-Nzon
       !stop
       !(ws(2)-ws(1))-(ws(4)-ws(3)),
       !stop
      
      !do i1=1,nc_eta1+1  
       !do i2=1,nc_eta2+1  
         !x1n_array(i1,i2) = x1_min+real(i1-1,f64)*delta_x1
       !enddo
      !enddo
     
 
     do i1=1,nc_eta1+1   

       dxx=(ws(1)-x2_min)/((nc_eta2-Nzon)/2._f64)
      
       do i2=1, (nc_eta2-Nzon)/2+1

       x2n_array(i1,i2) =x2_min+real(i2-1,f64)*dxx
       x1n_array(i1,i2) = x1_min+real(i1-1,f64)*delta_x1
       x1c_array(i1,i2) = x1_min+(real(i1,f64)-0.5_f64)*delta_x1
       x2c_array(i1,i2) = x2_min+(real(i2,f64)-0.5_f64)*dxx
       !jac_array(i1,i2) = (x1_max-x1_min)*(x2_max-x2_min)*(ws(1)-x2_min)/((N_x2-Nzon)*delta_x2/2._f64)
       jac_array(i1,i2) = (x1_max-x1_min)*(ws(1)-x2_min)/(slope_mesh2)
       !jac_array(i1,i2) = (x1_max-x1_min)*(x2_max-x2_min)/slope_mesh2
       integration_points(1,i1,i2) = (real(i1,f64)-0.5_f64)*(eta1_max-eta1_min)/real(nc_eta1,f64)
       integration_points(2,i1,i2) =x2_min+(real(i2,f64)-0.5_f64)*dxx
       
      enddo
     enddo

     do i1=1,nc_eta1+1   

       dxx=(ws(4)-ws(1))/(Nzon)
      
       do i2=(nc_eta2-Nzon)/2+1, (nc_eta2-Nzon)/2+Nzon +1

       x2n_array(i1,i2) =ws(1)+real(i2-1-(nc_eta2-Nzon)/2,f64)*dxx
       x1n_array(i1,i2) = x1_min+real(i1-1,f64)*delta_x1
       x1c_array(i1,i2) = x1_min+(real(i1,f64)-0.5_f64)*delta_x1
       x2c_array(i1,i2) = ws(1)+real(i2-(nc_eta2-Nzon)/2-0.5_f64,f64)*dxx
       jac_array(i1,i2) = (x1_max-x1_min)*(ws(4)-ws(1))/(slope_mesh1)
        !jac_array(i1,i2) = (x1_max-x1_min)*(x2_max-x2_min)/slope_mesh1
                       
       integration_points(1,i1,i2) = (real(i1,f64)-0.5_f64)*(eta1_max-eta1_min)/real(nc_eta1,f64)
       integration_points(2,i1,i2) = ws(1)+real(i2-(nc_eta2-Nzon)/2-0.5_f64,f64)*dxx
      enddo
     enddo

     do i1=1,nc_eta1+1   

       dxx=(x2_max-ws(4))/((nc_eta2-Nzon)/2._f64)
      
       do i2= (nc_eta2-Nzon)/2+Nzon +1, nc_eta2+1

       x2n_array(i1,i2) =ws(4)+real(i2-1-((nc_eta2-Nzon)/2+Nzon),f64)*dxx
       x1n_array(i1,i2) = x1_min+real(i1-1,f64)*delta_x1
       x1c_array(i1,i2) = x1_min+(real(i1,f64)-0.5_f64)*delta_x1
       x2c_array(i1,i2)  =ws(4)+real(i2-((nc_eta2-Nzon)/2+Nzon)-0.5_f64,f64)*dxx
       jac_array(i1,i2) = (x1_max-x1_min)*(ws(1)-x2_min)/(slope_mesh2)
       !jac_array(i1,i2) = (x1_max-x1_min)*(x2_max-x2_min)/slope_mesh2
        
       integration_points(1,i1,i2) = (real(i1,f64)-0.5_f64)*(eta1_max-eta1_min)/real(nc_eta1,f64)
       integration_points(2,i1,i2) =ws(4)+real(i2-((nc_eta2-Nzon)/2+Nzon)-0.5_f64,f64)*dxx
      enddo
     enddo

        !eta1 value of intersecting point (eta2,x1)=constant
        !integration_points(1,i1,i2) = (real(i1,f64)-0.5_f64)*(eta1_max-eta1_min)/real(nc_eta1,f64)
        !x2 value of intersecting point (eta2,x1)=constant
        !integration_points(2,i1,i2) = x2_min+(real(i2,f64)-0.5_f64)*delta_x2
     


    !geom => new_geometry_2D ('from_array',nc_eta1+1,nc_eta2+1, &
    !   x1n_array, x2n_array, x1c_array, x2c_array, jac_array,PERIODIC,PERIODIC)
    !mesh => new_mesh_descriptor_2D(eta1_min, eta1_max, nc_eta1, &
    !   PERIODIC, eta2_min, eta2_max, nc_eta2, PERIODIC, geom)
       
    !dist_func => sll_new_distribution_function_2D(mesh,CELL_CENTERED_DF, name)

    !do i2=1,nc_eta2+1
      !do i1=1,nc_eta1+1
        !eta1 value of intersecting point (eta2,x1)=constant
        !integration_points(1,i1,i2) = (real(i1,f64)-0.5_f64)*(eta1_max-eta1_min)/real(nc_eta1,f64)
        !x2 value of intersecting point (eta2,x1)=constant
        !integration_points(2,i1,i2) = x2_min+(real(i2,f64)-0.5_f64)*delta_x2
      !enddo
    !enddo
    
   
     
    
  endif
! do i1=1,nc_eta1 +1  
!       do i2= 1, N_x2+1
!         print*, x1c_array(i1,i2),x2c_array(i1,i2),ws(1),ws(2),ws(3),ws(4)!x1c_array(i1,i2),x2c_array(i1,i2),integration_points(1,i1,i2),integration_points(2,i1,i2)  
!      enddo
!     enddo
!stop
  if(mesh_case==11)then

      slope_mesh1=1._f64/2._f64
      wk=2.82_f64
      ll=1._f64
       
      Ns(1)=int((-wk-ll-x2_min)/delta_x2)+1
      ws(1)=x2_min+(Ns(1)-1)*delta_x2
      
      ll=(int(ll/delta_x2)+1)*delta_x2
     ! print*,ll
     !stop
      !Ns(2)=int((-wk+ll-x2_min)/delta_x2)+1
      ws(2)=ws(1)+2*ll
      
     ! Ns(3)=int((wk-ll-x2_min)/delta_x2)
      !ws(3)=x2_min+(Ns(3)-1)*delta_x2
      Ns(4)=int((wk+ll-x2_min)/delta_x2)+1
      ws(4)=x2_max-(ws(1)-x2_min)!x2_min+(Ns(4)-1)*delta_x2
      ws(3)=ws(4)-2*ll
     !iiw2=int((-wk-1-x2_min)/delta_x2)
    ! print*,ws(1)-ws(2),ws(3)-ws(4),delta_x2
     !stop
     !stop
     Nzon= int((ws(2)-ws(1))/(slope_mesh1*delta_x2))

     !print*,Nzon,int((ws(4)-ws(3))/(slope_mesh1*delta_x2))
     !stop
     Nzon2=3!int((ws(3)-ws(2))/(delta_x2))

     if(mod(nc_eta2-2*Nzon-Nzon2,2)==0) then
      Nzon2=Nzon2
     else
      Nzon2=Nzon2+1
     endif 
     ! print*,"slop",slope_mesh1,Nzon,N_x2-2*Nzon
      !stop
      !if(
      !Nzon= int(1._f64/(slope_mesh1*delta_x2))
      !slope_mesh2=(x2_max-x2_min-(ws(4)-ws(1)))/((N_x2-Nzon)*delta_x2)
        slope_mesh2=(ws(1)-x2_min)/((nc_eta2-2*Nzon-Nzon2)*delta_x2/2._f64)

        slope_mesh1=(ws(2)-ws(1))/(Nzon*delta_x2)
        slope_mesh3=(ws(3)-ws(2))/(Nzon2*delta_x2)
       
       !slope_mesh2=(x2_max-x2_min-2)/((N_x2-2*Nzon)*delta_x2)
       !print*,slope_mesh1,slope_mesh2!, (ws(1)-x2_min)/((N_x2-Nzon)*delta_x2/2._f64),Nzon,N_x2-Nzon
       !(ws(2)-ws(1))-(ws(4)-ws(3)),
       !stop
      
      !do i1=1,nc_eta1+1  
       !do i2=1,nc_eta2+1  
         !x1n_array(i1,i2) = x1_min+real(i1-1,f64)*delta_x1
       !enddo
      !enddo
     
 
     do i1=1,nc_eta1+1   

       dxx=(ws(1)-x2_min)/((nc_eta2-2*Nzon-Nzon2)/2._f64)
      
       do i2=1, (nc_eta2-2*Nzon-Nzon2)/2+1

       x2n_array(i1,i2) =x2_min+real(i2-1,f64)*dxx
       x1n_array(i1,i2) = x1_min+real(i1-1,f64)*delta_x1
       x1c_array(i1,i2) = x1_min+(real(i1,f64)-0.5_f64)*delta_x1
       x2c_array(i1,i2) = x2_min+(real(i2,f64)-0.5_f64)*dxx
       
       jac_array(i1,i2) = (x1_max-x1_min)*(ws(1)-x2_min)/slope_mesh2
       
       integration_points(1,i1,i2) = (real(i1,f64)-0.5_f64)*(eta1_max-eta1_min)/real(nc_eta1,f64)
       integration_points(2,i1,i2) =x2_min+(real(i2,f64)-0.5_f64)*dxx
       !print*,x1n_array(i1,i2),x2n_array(i1,i2)
      enddo
     enddo
     
     

     do i1=1,nc_eta1+1   

       dxx=(ws(2)-ws(1))/(Nzon)
      
       do i2=(nc_eta2-2*Nzon-Nzon2)/2+1, (nc_eta2-2*Nzon-Nzon2)/2+Nzon +1

       x2n_array(i1,i2) =ws(1)+real(i2-1-(nc_eta2-2*Nzon-Nzon2)/2,f64)*dxx
       x1n_array(i1,i2) = x1_min+real(i1-1,f64)*delta_x1
       x1c_array(i1,i2) = x1_min+(real(i1,f64)-0.5_f64)*delta_x1
       x2c_array(i1,i2) = ws(1)+real(i2-(nc_eta2-2*Nzon-Nzon2)/2-0.5_f64,f64)*dxx
       jac_array(i1,i2) = (x1_max-x1_min)*(ws(2)-ws(1))/slope_mesh1
       integration_points(1,i1,i2) = (real(i1,f64)-0.5_f64)*(eta1_max-eta1_min)/real(nc_eta1,f64)
       integration_points(2,i1,i2) = ws(1)+real(i2-(nc_eta2-2*Nzon-Nzon2)/2-0.5_f64,f64)*dxx
       !print*,x1n_array(i1,i2),x2n_array(i1,i2)
      enddo
     enddo

     do i1=1,nc_eta1+1   

       dxx=(ws(3)-ws(2))/(Nzon2)
      
       do i2= (nc_eta2-2*Nzon-Nzon2)/2+Nzon +1, (nc_eta2-2*Nzon-Nzon2)/2+Nzon +Nzon2 +1

       x2n_array(i1,i2) =ws(2)+real(i2-1-((nc_eta2-2*Nzon-Nzon2)/2+Nzon),f64)*dxx
       x1n_array(i1,i2) = x1_min+real(i1-1,f64)*delta_x1
       x1c_array(i1,i2) = x1_min+(real(i1,f64)-0.5_f64)*delta_x1
       x2c_array(i1,i2) = ws(2)+real(i2-((nc_eta2-2*Nzon-Nzon2)/2+Nzon)-0.5_f64,f64)*dxx
       jac_array(i1,i2) = (x1_max-x1_min)*(ws(3)-ws(2))/slope_mesh3
       integration_points(1,i1,i2) = (real(i1,f64)-0.5_f64)*(eta1_max-eta1_min)/real(nc_eta1,f64)
       integration_points(2,i1,i2) = ws(2)+real(i2-((nc_eta2-2*Nzon-Nzon2)/2+Nzon)-0.5_f64,f64)*dxx
       !print*,x1n_array(i1,i2),x2n_array(i1,i2)
      enddo
     enddo 

       do i1=1,nc_eta1+1   

       !dxx=(ws(4)-ws(3))/(Nzon)
        dxx=(ws(2)-ws(1))/(Nzon)
       do i2=  (nc_eta2-2*Nzon-Nzon2)/2+Nzon +Nzon2 +1, (nc_eta2-2*Nzon-Nzon2)/2+2*Nzon +Nzon2 +1

       x2n_array(i1,i2) =ws(3)+real(i2-1-((nc_eta2-2*Nzon-Nzon2)/2+Nzon+Nzon2),f64)*dxx
       x1n_array(i1,i2) = x1_min+real(i1-1,f64)*delta_x1
       x1c_array(i1,i2) = x1_min+(real(i1,f64)-0.5_f64)*delta_x1
       x2c_array(i1,i2) = ws(3)+real(i2-((nc_eta2-2*Nzon-Nzon2)/2+Nzon+Nzon2)-0.5_f64,f64)*dxx
       jac_array(i1,i2) = (x1_max-x1_min)*(ws(2)-ws(1))/slope_mesh1
       integration_points(1,i1,i2) = (real(i1,f64)-0.5_f64)*(eta1_max-eta1_min)/real(nc_eta1,f64)
       integration_points(2,i1,i2) = ws(3)+real(i2-((nc_eta2-2*Nzon-Nzon2)/2+Nzon+Nzon2)-0.5_f64,f64)*dxx
      ! print*,x1n_array(i1,i2),x2n_array(i1,i2)
      enddo
     enddo 

     
     do i1=1,nc_eta1+1   

       !dxx=(x2_max-ws(4))/((N_x2-Nzon)/2._f64)
        dxx=(x2_max-ws(4))/((nc_eta2-2*Nzon-Nzon2)/2._f64)
       
       do i2= (nc_eta2-2*Nzon-Nzon2)/2+2*Nzon +Nzon2 +1, nc_eta2+1

       x2n_array(i1,i2) =ws(4)+real(i2-1-((nc_eta2-2*Nzon-Nzon2)/2+2*Nzon+Nzon2),f64)*dxx
       x1n_array(i1,i2) = x1_min+real(i1-1,f64)*delta_x1
       x1c_array(i1,i2) = x1_min+(real(i1,f64)-0.5_f64)*delta_x1
       x2c_array(i1,i2)  =ws(4)+real(i2-((nc_eta2-2*Nzon-Nzon2)/2+2*Nzon+Nzon2)-0.5_f64,f64)*dxx
       !jac_array(i1,i2) = (x1_max-x1_min)*(x2_max-x2_min)*slope_mesh2
       jac_array(i1,i2) = (x1_max-x1_min)*(ws(1)-x2_min)/slope_mesh2
       integration_points(1,i1,i2) = (real(i1,f64)-0.5_f64)*(eta1_max-eta1_min)/real(nc_eta1,f64)
       integration_points(2,i1,i2) =ws(4)+real(i2-((nc_eta2-2*Nzon-Nzon2)/2+2*Nzon+Nzon2)-0.5_f64,f64)*dxx
        !print*,x1n_array(i1,i2),x2n_array(i1,i2)
      enddo
     enddo
     !stop
        !eta1 value of intersecting point (eta2,x1)=constant
        !integration_points(1,i1,i2) = (real(i1,f64)-0.5_f64)*(eta1_max-eta1_min)/real(nc_eta1,f64)
        !x2 value of intersecting point (eta2,x1)=constant
        !integration_points(2,i1,i2) = x2_min+(real(i2,f64)-0.5_f64)*delta_x2
     


    !geom => new_geometry_2D ('from_array',nc_eta1+1,nc_eta2+1, &
    !   x1n_array, x2n_array, x1c_array, x2c_array, jac_array,PERIODIC,PERIODIC)
    !mesh => new_mesh_descriptor_2D(eta1_min, eta1_max, nc_eta1, &
    !   PERIODIC, eta2_min, eta2_max, nc_eta2, PERIODIC, geom)
       
    !dist_func => sll_new_distribution_function_2D(mesh,CELL_CENTERED_DF, name)

    !do i2=1,nc_eta2+1
      !do i1=1,nc_eta1+1
        !eta1 value of intersecting point (eta2,x1)=constant
        !integration_points(1,i1,i2) = (real(i1,f64)-0.5_f64)*(eta1_max-eta1_min)/real(nc_eta1,f64)
        !x2 value of intersecting point (eta2,x1)=constant
        !integration_points(2,i1,i2) = x2_min+(real(i2,f64)-0.5_f64)*delta_x2
      !enddo
    !enddo
    
   
     !print*,slope_mesh1,slope_mesh2,slope_mesh3,
   !print*, (2*(ws(1)-x2_min)/slope_mesh2+2*(ws(2)-ws(1))/slope_mesh1+(ws(3)-ws(2))/slope_mesh3)*(x1_max-x1_min), (x1_max-x1_min)*(x2_max-x2_min)!,2*sll_Pi
    
  endif
   

  
  
  
  
  
  
  end subroutine construct_bgk_mesh  

 
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! functions for computing integral
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  
  function compute_non_unif_integral(integration_points,N_points,rho_case)
  sll_real64 :: compute_non_unif_integral
  sll_real64,dimension(:,:),pointer :: integration_points
  sll_int,intent(in) :: N_points,rho_case
  if(rho_case==1)then
    compute_non_unif_integral= compute_non_unif_integral_trapezoid(integration_points,N_points)
  endif
  if(rho_case==2)then
    compute_non_unif_integral=compute_non_unif_integral_spline(integration_points,N_points)
  endif
  if(rho_case==3)then
    compute_non_unif_integral=compute_non_unif_integral_gaussian(integration_points,N_points)
  endif      
  if(rho_case==4)then
    compute_non_unif_integral=compute_non_unif_integral_gaussian_sym(integration_points,N_points)
  endif        
  if(rho_case==5)then
    compute_non_unif_integral=compute_non_unif_integral_spline_per(integration_points,N_points)
  endif
  
end  function compute_non_unif_integral



function compute_non_unif_integral_trapezoid(integration_points,N_points)
  sll_real64 :: compute_non_unif_integral_trapezoid
  sll_real64,dimension(:,:),pointer :: integration_points
  sll_int,intent(in) :: N_points
  sll_int :: i
  sll_real64 :: tmp,x1,x2,fval1,fval2
  compute_non_unif_integral_trapezoid = 0._f64
  if(N_points<=1)then
    print *,'bad value of N_points=',N_points
    stop
  endif
  do i=1,N_points-1
    x1 = integration_points(1,i)
    x2 = integration_points(1,i+1)
    if(x2<x1)then
      print *,i,'bad integration points x1=',x1,'x2=',x2
      stop
    endif
    fval1 = integration_points(2,i)
    fval2 = integration_points(2,i+1)
    tmp = 0.5_f64*(fval1+fval2)*(x2-x1)
    compute_non_unif_integral_trapezoid=compute_non_unif_integral_trapezoid+tmp
  enddo
  
  
end  function compute_non_unif_integral_trapezoid


function compute_non_unif_integral_spline_old(integration_points,N_points,Nb)
  sll_real64 :: compute_non_unif_integral_spline_old
  sll_real64,dimension(:,:),pointer :: integration_points
  sll_real64,dimension(:,:),pointer :: integration_points_fine
  sll_int,intent(in) :: N_points,Nb
  sll_int :: i,N_points_fine,ierr,j
  sll_real64 :: x1,x2
  type(cubic_nonunif_spline_1D), pointer :: spl
  compute_non_unif_integral_spline_old = 0._f64
  if(N_points<=1)then
    print *,'bad value of N_points=',N_points
    stop
  endif
  N_points_fine = (N_points-1)*Nb+1
  spl =>  new_cubic_nonunif_spline_1D( N_points-1, SLL_HERMITE)
  SLL_ALLOCATE(integration_points_fine(2,N_points_fine),ierr)
  do i=1,N_points-1
    x1 = integration_points(1,i)
    x2 = integration_points(1,i+1)
    if(x2<x1)then
      print *,i,'(spl) bad integration points x1=',x1,'x2=',x2
      stop
    endif
    do j=1,Nb
      integration_points_fine(1,(i-1)*Nb+j)=x1+real(j-1,f64)/real(Nb,f64)*(x2-x1)
    enddo
  enddo  
  integration_points_fine(1,N_points_fine)=integration_points(1,N_points)
  !stop  
  call compute_spline_nonunif( integration_points(2,1:N_points), spl, integration_points(1,1:N_points),0._f64,0._f64)
  call interpolate_array_value_nonunif( integration_points_fine(1,1:N_points_fine), &
  &integration_points_fine(2,1:N_points_fine),N_points_fine-1, spl)
  call delete_cubic_nonunif_spline_1D( spl, ierr)

  do i=1,N_points_fine-1
     x1 = integration_points_fine(1,i)
     x2 = integration_points_fine(1,i+1)
     if(x2<x1)then
       print *,i,'(spl2) bad integration points x1=',x1,'x2=',x2
       stop
     endif
     !print *,i,integration_points_fine(1,i)
  enddo
  
  
  compute_non_unif_integral_spline_old = compute_non_unif_integral_trapezoid(integration_points_fine,N_points_fine)
  SLL_DEALLOCATE_ARRAY(integration_points_fine,ierr)
  
end  function compute_non_unif_integral_spline_old

function compute_non_unif_integral_spline(integration_points,N_points)
  sll_real64 :: compute_non_unif_integral_spline
  sll_real64,dimension(:,:),pointer :: integration_points
  sll_real64,dimension(:,:),pointer :: integration_points_middle
  sll_int,intent(in) :: N_points
  sll_int :: i,ierr
  sll_real64 :: tmp,x1,x2,fval1,fval2,fvalm
  type(cubic_nonunif_spline_1D), pointer :: spl
  compute_non_unif_integral_spline = 0._f64

  if(N_points<=1)then
    print *,'bad value of N_points=',N_points
    stop
  endif
  spl =>  new_cubic_nonunif_spline_1D( N_points-1, SLL_HERMITE)
  SLL_ALLOCATE(integration_points_middle(2,N_points-1),ierr)
  do i=1,N_points-1
    x1 = integration_points(1,i)
    x2 = integration_points(1,i+1)
    integration_points_middle(1,i)=0.5_f64*(x1+x2)
  enddo  
  call compute_spline_nonunif( integration_points(2,1:N_points), spl, integration_points(1,1:N_points),0._f64,0._f64)
  call interpolate_array_value_nonunif( integration_points_middle(1,1:N_points-1), &
  &integration_points_middle(2,1:N_points-1),N_points-1, spl)
  call delete_cubic_nonunif_spline_1D( spl, ierr)
  
  do i=1,N_points-1
    x1 = integration_points(1,i)
    x2 = integration_points(1,i+1)
    fval1 = integration_points(2,i)
    fval2 = integration_points(2,i+1)
    fvalm = integration_points_middle(2,i)
    tmp = (fval1+4._f64*fvalm+fval2)*(x2-x1)/6._f64
    compute_non_unif_integral_spline=compute_non_unif_integral_spline+tmp
  enddo
  
  SLL_DEALLOCATE_ARRAY(integration_points_middle,ierr)
  
end  function compute_non_unif_integral_spline


function compute_non_unif_integral_spline_per(integration_points,N_points)
  sll_real64 :: compute_non_unif_integral_spline_per
  sll_real64,dimension(:,:),pointer :: integration_points
  sll_real64,dimension(:,:),pointer :: integration_points_middle
  sll_int,intent(in) :: N_points
  sll_int :: i,ierr
  sll_real64 :: tmp,x1,x2,fval1,fval2,fvalm
  type(cubic_nonunif_spline_1D), pointer :: spl
  compute_non_unif_integral_spline_per = 0._f64

  if(N_points<=1)then
    print *,'bad value of N_points=',N_points
    stop
  endif
  spl =>  new_cubic_nonunif_spline_1D( N_points-1, SLL_PERIODIC)
  SLL_ALLOCATE(integration_points_middle(2,N_points-1),ierr)
  do i=1,N_points-1
    x1 = integration_points(1,i)
    x2 = integration_points(1,i+1)
    integration_points_middle(1,i)=0.5_f64*(x1+x2)
  enddo  
  call compute_spline_nonunif( integration_points(2,1:N_points), spl, integration_points(1,1:N_points))
  call interpolate_array_value_nonunif( integration_points_middle(1,1:N_points-1), &
  &integration_points_middle(2,1:N_points-1),N_points-1, spl)
  call delete_cubic_nonunif_spline_1D( spl, ierr)
  
  do i=1,N_points-1
    x1 = integration_points(1,i)
    x2 = integration_points(1,i+1)
    fval1 = integration_points(2,i)
    fval2 = integration_points(2,i+1)
    fvalm = integration_points_middle(2,i)
    tmp = (fval1+4._f64*fvalm+fval2)*(x2-x1)/6._f64
    compute_non_unif_integral_spline_per=compute_non_unif_integral_spline_per+tmp
  enddo
  
  SLL_DEALLOCATE_ARRAY(integration_points_middle,ierr)
  
end  function compute_non_unif_integral_spline_per



function compute_non_unif_integral_gaussian(integration_points,N_points)
  sll_real64 :: compute_non_unif_integral_gaussian
  sll_real64,dimension(:,:),pointer :: integration_points
  sll_real64,dimension(:,:),pointer :: integration_points_new
  sll_int,intent(in) :: N_points
  sll_int::i,ierr
  compute_non_unif_integral_gaussian = 0.5_f64*compute_non_unif_integral_gaussian_sym(integration_points,N_points)
  
  SLL_ALLOCATE(integration_points_new(2,N_points),ierr)
  
  
  integration_points_new(1,1:N_points)=integration_points(1,1:N_points)
  do i=1,N_points
    integration_points_new(2,i)=integration_points(2,N_points-i+1)
  enddo
  
  compute_non_unif_integral_gaussian = compute_non_unif_integral_gaussian+&
  0.5_f64*compute_non_unif_integral_gaussian_sym(integration_points_new,N_points)
  
  SLL_DEALLOCATE_ARRAY(integration_points_new,ierr)
  
end  function compute_non_unif_integral_gaussian


function compute_non_unif_integral_gaussian_sym(integration_points,N_points)
  sll_real64 :: compute_non_unif_integral_gaussian_sym
  sll_real64,dimension(:,:),pointer :: integration_points
  !sll_real64,dimension(:,:),pointer :: integration_points_new
  sll_real64,dimension(:,:),allocatable :: integration_points_new
  sll_int,intent(in) :: N_points
  sll_int :: i,ierr,j,is_center_point,N_points_new
  sll_real64 :: tmp,x1,x2,x3,fval1,fval2,fval3,x4,fval4,dx_int
  sll_int :: N_int,d_gauss,j_gauss
  !sll_real64,dimension(:,:),pointer :: gauss_points
  sll_real64,dimension(:,:),allocatable :: gauss_points
  compute_non_unif_integral_gaussian_sym = 0._f64
  N_int = 2
  d_gauss = 10
  
  SLL_ALLOCATE(gauss_points(2,d_gauss+1),ierr)
  
  if(d_gauss==0)then
    gauss_points(1,1) = 0.5_f64
    gauss_points(2,1) = 1._f64
  endif

  if(d_gauss==2)then
    gauss_points(1,1) = 0.11270166537925831148_f64 
    gauss_points(1,2) = 0.50000000000000000000_f64
    gauss_points(1,3) = 0.88729833462074168852_f64
    gauss_points(2,1) = 0.27777777777777777775_f64
    gauss_points(2,2) = 0.44444444444444444445_f64
    gauss_points(2,3) = 0.27777777777777777778_f64
  endif

  if(d_gauss==4)then
    gauss_points(1,1) = 0.046910077030668003601_f64 
    gauss_points(1,2) = 0.23076534494715845448_f64
    gauss_points(1,3) = 0.50000000000000000000_f64
    gauss_points(1,4) = 0.76923465505284154552_f64
    gauss_points(1,5) = 0.95308992296933199640_f64
    gauss_points(2,1) = 0.11846344252809454382_f64
    gauss_points(2,2) = 0.23931433524968323302_f64
    gauss_points(2,3) = 0.28444444444444444382_f64
    gauss_points(2,4) = 0.23931433524968323408_f64
    gauss_points(2,5) = 0.11846344252809454332_f64
  endif

  if(d_gauss==6)then
    gauss_points(1,1) = 0.025446043828620737737_f64 
    gauss_points(1,2) = 0.12923440720030278007_f64
    gauss_points(1,3) = 0.29707742431130141655_f64
    gauss_points(1,4) = 0.50000000000000000000_f64
    gauss_points(1,5) = 0.70292257568869858345_f64
    gauss_points(1,6) = 0.87076559279969721993_f64
    gauss_points(1,7) = 0.97455395617137926226_f64
    gauss_points(2,1) = 0.064742483084434846538_f64
    gauss_points(2,2) = 0.13985269574463833578_f64
    gauss_points(2,3) = 0.19091502525255949935_f64
    gauss_points(2,4) = 0.20897959183673468797_f64
    gauss_points(2,5) = 0.19091502525255948899_f64
    gauss_points(2,6) = 0.13985269574463833022_f64
    gauss_points(2,7) = 0.064742483084434844832_f64
  endif

  if(d_gauss==8)then
    gauss_points(1,1) = 0.015919880246186955082_f64 
    gauss_points(1,2) = 0.081984446336682102850_f64
    gauss_points(1,3) = 0.19331428364970480135_f64
    gauss_points(1,4) = 0.33787328829809553548_f64
    gauss_points(1,5) = 0.50000000000000000000_f64
    gauss_points(1,6) = 0.66212671170190446452_f64
    gauss_points(1,7) = 0.80668571635029519865_f64
    gauss_points(1,8) = 0.91801555366331789715_f64
    gauss_points(1,9) = 0.98408011975381304492_f64
    gauss_points(2,1) = 0.040637194180787175822_f64
    gauss_points(2,2) = 0.090324080347428468792_f64
    gauss_points(2,3) = 0.13030534820146771147_f64
    gauss_points(2,4) = 0.15617353852000146170_f64
    gauss_points(2,5) = 0.16511967750063012022_f64
    gauss_points(2,6) = 0.15617353852000133681_f64
    gauss_points(2,7) = 0.13030534820146773892_f64
    gauss_points(2,8) = 0.090324080347428458790_f64
    gauss_points(2,9) = 0.040637194180787192773_f64
  endif

  if(d_gauss==10)then
    gauss_points(1,1) = 0.010885670926971503598_f64 
    gauss_points(1,2) = 0.056468700115952350462_f64
    gauss_points(1,3) = 0.13492399721297533795_f64
    gauss_points(1,4) = 0.24045193539659409204_f64
    gauss_points(1,5) = 0.36522842202382751383_f64
    gauss_points(1,6) = 0.50000000000000000000_f64
    gauss_points(1,7) = 0.63477157797617248617_f64
    gauss_points(1,8) = 0.75954806460340590796_f64
    gauss_points(1,9) = 0.86507600278702466205_f64
    gauss_points(1,10) = 0.94353129988404764954_f64
    gauss_points(1,11) = 0.98911432907302849640_f64
    gauss_points(2,1) = 0.027834283558086630522_f64
    gauss_points(2,2) = 0.062790184732454913642_f64
    gauss_points(2,3) = 0.093145105463867582027_f64
    gauss_points(2,4) = 0.11659688229597882542_f64
    gauss_points(2,5) = 0.13140227225513898990_f64
    gauss_points(2,6) = 0.13646254338895855961_f64
    gauss_points(2,7) = 0.13140227225511460324_f64
    gauss_points(2,8) = 0.11659688229599488815_f64
    gauss_points(2,9) = 0.093145105463868674217_f64
    gauss_points(2,10) = 0.062790184732451937905_f64
    gauss_points(2,11) = 0.027834283558088336992_f64
  endif


  
  if(N_points<=1)then
    print *,'bad value of N_points=',N_points
    stop
  endif
  i=N_points/2
  is_center_point=0
  if(2*i/=N_points)then
    is_center_point=1
  endif
  !check that the integration points are increasing
  do i=1,N_points-1
    if(integration_points(1,i+1)<=integration_points(1,i))then
      print *,'order problem for integration_points',i,integration_points(1,i),integration_points(1,i+1)
      do j=1,N_points
        print *,j,integration_points(1,j)
      enddo
      stop
    endif
  enddo
  !check for the symmetry  
  tmp=0._f64
  if(is_center_point==0)then
    do i=1,N_points/2
      if(abs(integration_points(1,N_points/2+i)+integration_points(1,N_points/2-i+1))>tmp)then
        tmp=abs(integration_points(1,N_points/2+i)+integration_points(1,N_points/2-i+1))
      endif
    enddo
    if(tmp>1.e-13)then
      print *,'integration_points are not symmetric',tmp
      do j=1,N_points
        print *,j,integration_points(1,j)
      enddo
      stop
    endif
  endif
  
  if(is_center_point==1)then
    tmp=0._f64
    tmp=abs(integration_points(1,(N_points+1)/2))
    do i=1,(N_points-1)/2
      if(abs(integration_points(1,(N_points+1)/2+i)+integration_points(1,(N_points+1)/2-i))>tmp)then
        tmp=abs(integration_points(1,(N_points+1)/2+i)+integration_points(1,(N_points+1)/2-i))
      endif
    enddo
    if(tmp>1.e-14)then
      print *,'integration_points are not symmetric',tmp
      do j=1,N_points
        print *,j,integration_points(1,j)
      enddo
      stop
    endif
  endif
  
  !we will store in a new tab so that we are in the even case
  N_points_new = N_points
  if(is_center_point==0)then
    N_points_new= N_points+1
  endif
  SLL_ALLOCATE(integration_points_new(2,N_points_new),ierr)
  if(is_center_point==1)then
    integration_points_new = integration_points
  endif
  if(is_center_point==0)then
    if(N_points/2+3>N_points)then
      print *,'N_points is too small',N_points
    endif
    integration_points_new(1:2,1:N_points/2) = integration_points(1:2,1:N_points/2)
    integration_points_new(1:2,N_points/2+2:N_points+1) = integration_points(1:2,N_points/2+1:N_points)
    !we have to predict the value of the center
    integration_points_new(1,N_points/2+1)=0._f64
    x1=integration_points(1,N_points/2+1)
    x2=integration_points(1,N_points/2+2)
    x3=integration_points(1,N_points/2+3)
    fval1=integration_points(2,N_points/2+1)*exp(0.5_f64*x1*x1)
    fval2=integration_points(2,N_points/2+2)*exp(0.5_f64*x2*x2)
    fval3=integration_points(2,N_points/2+3)*exp(0.5_f64*x3*x3)
    !print *,x1,x2,x3,fval1,fval2,fval3
    x4=0._f64
    fval4=(fval1*(x4-x2)*(x4-x3)/((x1-x2)*(x1-x3))+fval2*(x4-x1)*(x4-x3)/((x2-x1)*(x2-x3))+fval3*(x4-x1)*(x4-x2)/((x3-x1)*(x3-x2)))
    fval4=fval4*exp(-0.5_f64*x4*x4)
    integration_points_new(2,N_points/2+1)=fval4
  endif
  
  do i=1,(N_points_new-1)/4
    x1 = integration_points_new(1,(N_points_new-1)/2+2*i-1)
    x2 = integration_points_new(1,(N_points_new-1)/2+2*i)
    x3 = integration_points_new(1,(N_points_new-1)/2+2*i+1)
    fval1 = integration_points_new(2,(N_points_new-1)/2+2*i-1)*exp(0.5_f64*x1*x1)
    fval2 = integration_points_new(2,(N_points_new-1)/2+2*i)*exp(0.5_f64*x2*x2)
    fval3 = integration_points_new(2,(N_points_new-1)/2+2*i+1)*exp(0.5_f64*x3*x3)
    tmp=0._f64
    dx_int=(x3-x1)/real(N_int,f64)
    do j=1,N_int
      do j_gauss=1,d_gauss+1
        !x4 =x1+(real(j,f64)-0.5_f64)*dx_int
        x4 =x1+(real(j-1,f64)+gauss_points(1,j_gauss))*dx_int
        fval4=(fval1*(x4-x2)*(x4-x3)/((x1-x2)*(x1-x3))+fval2*(x4-x1)*(x4-x3)/((x2-x1)*(x2-x3))+fval3*(x4-x1)*(x4-x2)/((x3-x1)*(x3-x2)))
        fval4=fval4*exp(-0.5_f64*x4*x4)
        tmp=tmp+fval4*gauss_points(2,j_gauss)
      enddo
    enddo  
    tmp=tmp*dx_int
    !print *,i,x1,x2,x3
    compute_non_unif_integral_gaussian_sym =compute_non_unif_integral_gaussian_sym+tmp
    !integration_points_middle(1,i)=0.5_f64*(x1+x2)
  enddo
  j=(N_points_new-1)/4
  if(2*j/=(N_points_new-1)/2)then
    !print *,2*j,(N_points_new-1)/2
    if(2*j+1/=(N_points_new-1)/2)then
      print *,'Problem concerning N_points',2*j+1,(N_points_new-1)/2
    endif
    x1 = integration_points_new(1,N_points_new-2)
    x2 = integration_points_new(1,N_points_new-1)
    x3 = integration_points_new(1,N_points_new)
    fval1 = integration_points_new(2,N_points_new-2)*exp(0.5_f64*x1*x1)
    fval2 = integration_points_new(2,N_points_new-1)*exp(0.5_f64*x2*x2)
    fval3 = integration_points_new(2,N_points_new)*exp(0.5_f64*x3*x3)
    x4 = 2._f64*x3-x2
    fval4=(fval1*(x4-x2)*(x4-x3)/((x1-x2)*(x1-x3))+fval2*(x4-x1)*(x4-x3)/((x2-x1)*(x2-x3))+fval3*(x4-x1)*(x4-x2)/((x3-x1)*(x3-x2)))
    fval4=fval4*exp(-0.5_f64*x4*x4)


    
    x1 = integration_points_new(1,N_points_new-1)
    x2 = integration_points_new(1,N_points_new)
    x3 = 2._f64*x2-x1
    fval1 = integration_points_new(2,N_points_new-1)*exp(0.5_f64*x1*x1)
    fval2 = integration_points_new(2,N_points_new)*exp(0.5_f64*x2*x2)
    fval3 = fval4*exp(0.5_f64*x3*x3)
    
    !print *,x1,x2,x3
    !print *,fval1,fval2,fval3
    !print *,x1,x2,x3
    !print *,fval1,fval2,fval3
    
    !stop
    
    tmp=0._f64
    N_int =N_int*2
    dx_int=(x3-x1)/real(N_int,f64)
    do j=1,N_int
      do j_gauss=1,d_gauss+1
        !x4 =x1+(real(j,f64)-0.5_f64)*dx_int
        x4 =x1+(real(j-1,f64)+gauss_points(1,j_gauss))*dx_int
        fval4=(fval1*(x4-x2)*(x4-x3)/((x1-x2)*(x1-x3))+fval2*(x4-x1)*(x4-x3)/((x2-x1)*(x2-x3))&
        +fval3*(x4-x1)*(x4-x2)/((x3-x1)*(x3-x2)))
        fval4=fval4*exp(-0.5_f64*x4*x4)
        tmp=tmp+fval4*gauss_points(2,j_gauss)
      enddo
    enddo
    tmp=tmp*dx_int
    !print *,i,x1,x2,x3
    compute_non_unif_integral_gaussian_sym =compute_non_unif_integral_gaussian_sym+tmp
  endif
  
  
  
  
  compute_non_unif_integral_gaussian_sym = 2._f64*compute_non_unif_integral_gaussian_sym
  
  SLL_DEALLOCATE_ARRAY(integration_points_new,ierr)
  SLL_DEALLOCATE_ARRAY(gauss_points,ierr)

  !SLL_DEALLOCATE(integration_points_new,ierr)
  !SLL_DEALLOCATE(gauss_points,ierr)
  
end  function compute_non_unif_integral_gaussian_sym

  


end module simulation_VP1D_curvilinear_analytic
