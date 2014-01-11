module simulation_VP1D_cartesian_non_unif


!> Vlasov-Poisson 1D on uniform or nonuniform cartesian grid
!> using the Backward Semi-Lagrangian (BSL) method.
!> in development

#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_field_2d.h"


#ifndef STDF95
  use sll_simulation_base
#endif
  use cubic_non_uniform_splines
  use sll_constants
  implicit none

#ifdef STDF95
  type :: sll_simulation_VP1D_cartesian_non_unif
#else
  type, extends(sll_simulation_base_class) :: &
    sll_simulation_VP1D_cartesian_non_unif
#endif
    ! Numerical parameters
    sll_real64 :: dt
#ifndef STDF95
  contains
    procedure, pass(sim) :: run => run_VP1D_cartesian_non_unif
    procedure, pass(sim) :: init_from_file => VP1D_cart_non_unif_init
#endif
  end type sll_simulation_VP1D_cartesian_non_unif

contains

  subroutine VP1D_cart_non_unif_init(sim, filename)
#ifdef STDF95
    type(sll_simulation_VP1D_cartesian_non_unif), intent(inout)  :: sim
#else
    class(sll_simulation_VP1D_cartesian_non_unif), intent(inout) :: sim
#endif
    character(len=*), intent(in)                                 :: filename
    ! Declare here the variables to be read in through a namelist and that
    ! are to be kept inside the sim object. Look at the parallel vp4d simulation
    ! for an example.
    print *, 'This is a dummy function. Needs implementation.'
    
    print *,filename
    print *,sim%dt
    stop
    
  end subroutine VP1D_cart_non_unif_init

  ! Note that the following function has no local variables, which is silly...
  ! This just happened since the guts of the unit test were transplanted here
  ! directly, but this should be cleaned up.
  subroutine run_VP1D_cartesian_non_unif(sim)
    implicit none
#ifdef STDF95
    type(sll_simulation_VP1D_cartesian_non_unif), intent(inout) :: sim
#else
    class(sll_simulation_VP1D_cartesian_non_unif), intent(inout) :: sim
#endif



    
    type(cubic_nonunif_spline_1D), pointer :: spl_per_x1,spl_per_x2,spl_per_x1_poisson
    sll_int32 :: N_x1,N_x2,N_x1_poisson,N,nb_step
    sll_int32:: mesh_case,test_case,rho_case
    sll_int32:: i,i1,i2,err,step
    sll_real64 :: dt,x1_min,x1_max,x2_min,x2_max
    sll_real64 :: delta_x1,delta_x2,delta_x1_poisson
    sll_real64,dimension(:), pointer :: node_positions_x1,node_positions_x2
    sll_real64,dimension(:), pointer :: node_positions_x1_unit,node_positions_x2_unit
    sll_real64,dimension(:), pointer :: node_positions_x1_poisson
    sll_real64,dimension(:), pointer :: E,E_poisson,rho
    sll_real64,dimension(:), pointer :: buf1d,Xstar
    sll_real64,dimension(:), pointer :: random_vector_x1,random_vector_x2
    sll_real64,dimension(:,:), pointer :: f
    sll_real64,dimension(:,:), pointer :: integration_points
    sll_real64 :: landau_alpha,non_unif_scale_x1,non_unif_scale_x2
    sll_real64 :: val,x1,x2,tmp
    
    
    
    !initialization of parameters
    N_x1=256
    N_x2=N_x1
    N_x1_poisson=1024
    mesh_case=2 !1: uniform case 2: non_unif_scale perturbation
    test_case=1 !1: landau damping
    rho_case=1  !1: trap 2: splines 3: gaussian 4: gaussian sym 
                !5: spline_per
    landau_alpha=1e-3_f64
    dt=0.1_f64
    nb_step=601
    
    non_unif_scale_x1 = 0.1_f64
    non_unif_scale_x2 = 0.1_f64 
    
    N=max(N_x1,N_x2)
    
    print *,'#N_x1=',N_x1,'#N_x2=',N_x2
    print *,'#mesh_case=',mesh_case
    print *,'#test_case=',test_case
    print *,'#rho_case=',rho_case
    print *,'#landau_alpha=',landau_alpha
    print *,'#N_x1_poisson=',N_x1_poisson
    print *,'#nb_step=',nb_step
    x1_min = 0._f64
    x1_max = 4._f64*sll_pi
    x2_min = -6._f64
    x2_max = -x2_min
  
    delta_x1= (x1_max-x1_min)/real(N_x1,f64)
    delta_x2= (x2_max-x2_min)/real(N_x2,f64)
    
    delta_x1_poisson = (x1_max-x1_min)/real(N_x1_poisson,f64)
    
    spl_per_x1 =>  new_cubic_nonunif_spline_1D( N_x1, SLL_PERIODIC)
    spl_per_x2 =>  new_cubic_nonunif_spline_1D( N_x2, SLL_PERIODIC)

    spl_per_x1_poisson =>  new_cubic_nonunif_spline_1D( N_x1_poisson, SLL_PERIODIC)
    
    sim%dt = dt 

    SLL_ALLOCATE(node_positions_x1(N_x1+1),err)
    SLL_ALLOCATE(node_positions_x2(N_x2+1),err)
    SLL_ALLOCATE(node_positions_x1_unit(N_x1+1),err)
    SLL_ALLOCATE(node_positions_x2_unit(N_x2+1),err)
    SLL_ALLOCATE(node_positions_x1_poisson(N_x1_poisson+1),err)
    SLL_ALLOCATE(E_poisson(N_x1_poisson+1),err)
    SLL_ALLOCATE(rho(N_x1+1),err)
    SLL_ALLOCATE(E(N_x1+1),err)
    SLL_ALLOCATE(Xstar(N+1),err)
    SLL_ALLOCATE(buf1d(N+1),err)
    SLL_ALLOCATE(random_vector_x1(N_x1+1),err)
    SLL_ALLOCATE(random_vector_x2(N_x2+1),err)
    SLL_ALLOCATE(f(N_x1+1,N_x2+1),err)
    SLL_ALLOCATE(integration_points(2,N_x2+1),err)

    !begin with uniform mesh
    if(mesh_case==1)then
      do i=1,N_x1+1
        node_positions_x1(i) = x1_min+(real(i,f64)-1._f64)*delta_x1
      enddo
      do i=1,N_x2+1
        node_positions_x2(i) = x2_min+(real(i,f64)-1._f64)*delta_x2
      enddo
    endif
    

    !non uniform mesh
    call random_number(random_vector_x1)
    call random_number(random_vector_x2)
    
    random_vector_x1(1) = 0._f64
    random_vector_x1(N_x1+1) = 0._f64
    random_vector_x2(1) = 0._f64
    random_vector_x2(N_x2+1) = 0._f64


    
    if(mesh_case==2)then
    
      do i=1,N_x1+1
        node_positions_x1(i) = x1_min+&
        (real(i,f64)-1._f64+non_unif_scale_x1* random_vector_x1(i))*delta_x1
      enddo
      do i=1,N_x2+1
        node_positions_x2(i) = x2_min+&
        (real(i,f64)-1._f64+non_unif_scale_x2* random_vector_x2(i))*delta_x2
      enddo
    endif




    do i=1,N_x1_poisson+1
      node_positions_x1_poisson(i) = x1_min+(real(i,f64)-1._f64)*delta_x1_poisson
    enddo
    
    node_positions_x1_unit = (node_positions_x1-x1_min)/(x1_max-x1_min)
    node_positions_x2_unit = (node_positions_x2-x2_min)/(x2_max-x2_min)
    
    !initialization of distribution function
    do i1=1,N_x1+1
      do i2=1,N_x2+1
        x1 = node_positions_x1(i1)
        x2 = node_positions_x2(i2)
        !Landau damping
        if(test_case==1)then
          val = 1._f64/(sqrt(2._f64*sll_pi))*exp(-0.5_f64*x2*x2)
          val = val*(1._f64+landau_alpha*cos(2._f64*sll_pi/(x1_max-x1_min)*x1))
        endif
        f(i1,i2) = val      
    enddo
  enddo    
    
   !computation of rho
   integration_points(1,1:N_x2+1) = node_positions_x2(1:N_x2+1)   
   do i1 = 1, N_x1+1 
     integration_points(2,1:N_x2+1) = f(i1,1:N_x2+1) 
     rho(i1) = compute_non_unif_integral(integration_points,N_x2+1,rho_case)
   enddo
   !computation of E
   call compute_spline_nonunif( rho, spl_per_x1,node_positions_x1)  
   call interpolate_array_value_nonunif( node_positions_x1_poisson, &
    E_poisson, N_x1_poisson, spl_per_x1)  
   call poisson1dpertrap(E_poisson,x1_max-x1_min,N_x1_poisson)
   call compute_spline_nonunif( E_poisson, spl_per_x1_poisson,node_positions_x1_poisson)
   call interpolate_array_value_nonunif( node_positions_x1, &
    E, N_x1, spl_per_x1_poisson)

  do step=1,nb_step
  
    val=0._f64
    do i1=1,N_x1_poisson
      val = val+E_poisson(i1)*E_poisson(i1)
    enddo
    val = val/real(N_x1_poisson,f64)

    print *,(real(step,f64)-1._f64)*dt,val

    ! advect in x over dt/2 
    do i2=1,N_x2+1
      buf1d(1:N_x1+1) = f(1:N_x1+1,i2)
      Xstar(1:N_x1+1) = node_positions_x1(1:N_x1+1)&
      -0.5_f64*dt*node_positions_x2(i2)
      Xstar(1:N_x1+1) = (Xstar(1:N_x1+1)-x1_min)/(x1_max-x1_min)
      call csl_advection_per(buf1d,spl_per_x1,Xstar,&
        node_positions_x1_unit,N_x1)
      f(1:N_x1+1,i2) = buf1d(1:N_x1+1)
    enddo

   !computation of rho
   integration_points(1,1:N_x2+1) = node_positions_x2(1:N_x2+1)   
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

   !computation of E
   call compute_spline_nonunif( rho, spl_per_x1,node_positions_x1)  
   call interpolate_array_value_nonunif( node_positions_x1_poisson, &
    E_poisson, N_x1_poisson, spl_per_x1)  
   call poisson1dpertrap(E_poisson,x1_max-x1_min,N_x1_poisson)
   call compute_spline_nonunif( E_poisson, spl_per_x1_poisson,node_positions_x1_poisson)
   call interpolate_array_value_nonunif( node_positions_x1, &
    E, N_x1, spl_per_x1_poisson)


    ! advect in v over dt
    do i1=1,N_x1+1
      buf1d(1:N_x2+1) = f(i1,1:N_x2+1)
      Xstar(1:N_x2+1) = node_positions_x2(1:N_x2+1)-E(i1)*dt
      Xstar(1:N_x2+1)=(Xstar(1:N_x2+1)-x2_min)/(x2_max-x2_min)
      call csl_advection_per(buf1d,spl_per_x2,Xstar,&
        node_positions_x2_unit,N_x2)
      f(i1,1:N_x2+1) = buf1d(1:N_x2+1)
    enddo

    ! advect in x over dt/2 
    do i2=1,N_x2+1
      buf1d(1:N_x1+1) = f(1:N_x1+1,i2)
      Xstar(1:N_x1+1) = node_positions_x1(1:N_x1+1)&
      -0.5_f64*dt*node_positions_x2(i2)
      Xstar(1:N_x1+1) = (Xstar(1:N_x1+1)-x1_min)/(x1_max-x1_min)
      call csl_advection_per(buf1d,spl_per_x1,Xstar,&
        node_positions_x1_unit,N_x1)
      f(1:N_x1+1,i2) = buf1d(1:N_x1+1)
    enddo



  
  enddo

    
    
    
  end subroutine run_VP1D_cartesian_non_unif

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

  


end module simulation_VP1D_cartesian_non_unif
