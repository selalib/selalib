program cg_curvilinear_2D
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

   
  use sll_m_cg_curvi_structure
  use sll_m_curvilinear_2D_advection
  use sll_m_cg_curvi_function
  use sll_m_diagnostic
  !--*Poisson*----
  use sll_m_mudpack_colella
 ! use sll_general_coordinate_qn_solver
 ! use sll_m_cartesian_meshes
 ! use sll_m_common_coordinate_transformations
 ! use sll_m_scalar_field_2d
 ! use sll_m_constants
 ! use sll_m_arbitrary_degree_spline_interpolator_2d
  use sll_m_boundary_condition_descriptors

  implicit none

  type(sll_SL_curvilinear)     , pointer     :: plan_sl
  sll_real64, dimension (:,:)  , allocatable :: f,fp1,g,f_init 
  sll_real64, dimension (:,:)  , allocatable :: phi_Fourier
  sll_real64, dimension(:,:)   , pointer     :: x1n_array,x2n_array,x1_tab,x2_tab  
  sll_real64, dimension(:,:)   , pointer     :: jac_array
  sll_real64, dimension(:,:,:) , pointer     :: jac_matrix
  sll_real64, dimension(:,:)   , pointer     :: phi_exact

  sll_real64 :: delta_eta1, delta_eta2,error2,error1
  sll_real64 :: eta1_min, eta1_max,eta2_min, eta2_max, dt, tf
  sll_real64 :: temps,a1,a2,error
  sll_real64 :: x1c,x2c,x1c_r,x2c_r,sigma_x1,sigma_x2
  sll_real64 :: geom_eta(2,2),geom_x(2,2),alpha_mesh
  
  sll_int32 :: i,j,step,visu_step,hh,minn,ss
  sll_int32 :: N_eta1, N_eta2, nb_step
  sll_int32 :: f_case,carac_case,visu_case,phi_case
  sll_int32 :: bc(2),err,mode,solver
  sll_int32 :: PERIODIC_B,HERMITE_B
  sll_int32 :: bc1_type, bc2_type
  sll_int32 :: grad_case,mesh_case,time_scheme
  
  character (len=16) :: f_file !,bctop,bcbot
  character (len=8) :: conv_CG
  
  
!!**Poisson solver *****
   sll_real64, dimension (:,:) , allocatable :: values_U
   type(mudpack_2d)                          :: poisson_df
   
!!---* Diagnostic *----:
sll_real64, dimension (:,:)  , allocatable :: phi_ref,trans_phi
sll_real64, dimension (:)    , allocatable :: int_r
sll_real :: l10,l20,e0,alpha,r1,r2

!!---* Collela *---
sll_real :: landau_alpha,landau_mode
!!********

  !namelist /param/ N_eta1,N_eta2,dt,nb_step,carac_case,mesh_case,visu_case,phi_case
  
  solver =1
  visu_case = 1 ! 0 : gnuplot 
                ! 1 : vtk

  visu_step = 10

  mesh_case = 2 ! 1 : cartesian 
                ! 2 : polar 
                ! 3 : r^2 modified polar 
                ! 4 : colella 

  f_case = 1  ! 0 : constant function 
              ! 1 : 1+alpha*cos(mode*eta2)
              ! 3 : sin(eta_2)+landau_alpha*cos(landau_mode*eta_1)
              ! 4 : sll_m_gaussian in x and y

  grad_case = 2

  carac_case = 5 ! 1 : Explicit Euler with linear interpolation
                 ! 2 : Explicit Euler with spline interpolation  
                 ! 3 : Analytics caracteristics (cartisian,polar,collela)
                 ! 5 : Fixed point (midpoint)

  phi_case = 1 ! 1: translation 
               ! 2: rotation 
               ! 3: anisotropic rotation

  time_scheme = 1 ! 1 : SL_order 1  
                  ! 2 : SL order 2 (Predictor-Corrector) 
                  ! 3 : SL order 2 (Leap-Frog)
  
  !!----*default values of paramietrs*---------

  a1 = 0.25_f64 !*0.01_f64
  a2 = 0.25_f64 !*0.01_f64
  bc1_type= SLL_PERIODIC
  bc2_type= SLL_PERIODIC
  N_eta1 = 80
  N_eta2 = 64
  dt = 0.1
  nb_step = 1000
  alpha_mesh = 0._f64
  eta1_min = 0._f64
  eta1_max = 1._f64
  eta2_min = 0._f64
  eta2_max = 1._f64
  x1c_r=0._f64  ! center for rotation
  x2c_r=0._f64
  x1c = 0.5_f64 ! center of sll_m_gaussian 
  x2c = 0.5_f64
  sigma_x1 = 10._f64  ! variance for sll_m_gaussian (exp(-(x-x1c)^2/(2*sigma_1^2)))
  sigma_x2 = 10._f64
  
!read(*,NML=param)


  

  ! ---- * Construction of the mesh * ----
  
  ! mesh type : cartesian
  ! domain    : square [eta1_min eta1_max] x [eta2_min eta2_max]
  ! BC        : periodic-periodic
  if (mesh_case==1) then
    eta1_min = 0._f64
    eta1_max = 2._f64 !*sll_pi !6._f64
    eta2_min = 0._f64
    eta2_max = 2._f64 !*sll_pi !6._f64

    x1c_r=3._f64
    x2c_r=3._f64

    x1c = 2._f64
    x2c = 2._f64
    
    sigma_x1 = 0.15_f64
    sigma_x2 = 0.15_f64
    
    bc1_type = SLL_PERIODIC
    bc2_type = SLL_PERIODIC
  endif
  
  ! mesh type : polar or polar-like
  ! domain    : disc of radius eta1_max with a hole of radius eta1_min
  ! BC        : hermite-periodic
  if ((mesh_case==2).or.(mesh_case==3)) then

    eta1_min = 1._f64!0.5_f64
    eta1_max = 10._f64
    eta2_min = 0._f64
    eta2_max = 2._f64*sll_pi
    
    x1c_r=0._f64
    x2c_r=0._f64
   
    x1c = -3.5_f64
    x2c = 0._f64

    sigma_x1 = 0.2
    sigma_x2 = 0.2

    bc1_type = SLL_HERMITE
    bc2_type = SLL_PERIODIC
    bc(1)=1 !3
    bc(2)=1 
    mode=3
    alpha=1.e-6
    r1=4._f64
    r2=5._f64
  endif
  
  ! mesh type : Collela
  ! domain    : square [eta1_min eta1_max] x [eta2_min eta2_max]
  ! BC        : periodic-periodic
  if (mesh_case==4) then

    eta1_min = 0._f64
    eta1_max = 8._f64
    eta2_min = 0._f64
    eta2_max = 8._f64

    if(f_case==3) then!khp
    landau_mode=0.5_f64
    landau_alpha=0.015_f64
    dt=0.1_f64
    nb_step=600 
    visu_step=300
    eta1_min = 0._f64
    eta1_max = 2._f64*sll_pi/landau_mode
    eta2_min = 0._f64 
    eta2_max = 2._f64*sll_pi
    endif  

    !x1c_r=4._f64
    !x2c_r=4._f64
    !x1c = 2._f64
    !x2c = 2._f64
    !sigma_x1 = 0.15_f64
    !sigma_x2 = 0.15_f64
    
    bc1_type = SLL_PERIODIC
    bc2_type = SLL_PERIODIC 
    
    alpha_mesh = 0.02_f64
  endif
    
  geom_eta(1,1)=eta1_min
  geom_eta(2,1)=eta1_max
  geom_eta(1,2)=eta2_min
  geom_eta(2,2)=eta2_max

  delta_eta1= (eta1_max-eta1_min)/real(N_eta1,f64)
  delta_eta2= (eta2_max-eta2_min)/real(N_eta2,f64)



  PRINT*,'#----*Geometric Parameter*----'
    print*,'# Minimun for eta1 = ', eta1_min
    print*,'# Minimun for eta2 = ', eta2_min
    print*,'# Maximun for eta1 = ', eta1_max
    print*,'# Maximun for eta2 = ', eta2_max
    print*,'# Pas en direction eta1: delta_eta1 =',delta_eta1
    print*,'# Pas en direction eta2: delta_eta2 =',delta_eta2

  
  PRINT*,'#----*Advection Parameter*----'
    print*,'number segment in direction eta1, N_eta1   =  ', N_eta1   
    print*,'number segment in direction eta2, N_eta2   =  ', N_eta2 
    print*,'# time_scheme', time_scheme
    print*,'# f_case', f_case
    print*,'# carac_case', carac_case
    print*,'# grad_case', grad_case
    print*,'# mesh_case', mesh_case
    print*,'# alpha_mesh', alpha_mesh
    print*,'# bc1_type, bc2_type = 1:HERMITE, 0:PERIODIC', bc1_type,bc2_type
   
  PRINT*,'#----*Time Parameter*----'
    tf=dt*real(nb_step,f64)   !definition of dt=tf/nb_step
    print*,'# nb_step =',nb_step
    print*,'# Pas de temps: dt =',dt
    print*,'# Temps final: tf =',tf



  ! mesh type
   if (mesh_case > 4) then
    print*,'#Non existing case'
    print*,'#mesh_case = 1:cartesian, 2:polar, 3:polar-like or 4:Collela or 5: Collela2'
    STOP
  endif
  ! distribution function
  if (f_case > 5) then
    print*,'#Non existing case'
    print*,'#test_case = 1 :f=1, 2 : cos(eta2), 3 :sll_m_gaussian in eta1,'
    print*,'#4 :centered sll_m_gaussian in eta1 and eta2 depending on the mesh type or 5 :centered dirac'
    STOP
  endif
  ! advection field
  if (phi_case > 4) then
    print*,'#Non existing case'
    print*,'#field_case = 1 : translation, 2 :rotation, 3 :non homogeneous rotation'  
    STOP
  endif




 !--*Initialisation plan_sl for advection*----
plan_sl => new_SL(geom_eta,delta_eta1,delta_eta2,dt, &
                  & N_eta1,N_eta2,grad_case,carac_case,bc,bc1_type,bc2_type)

!  !ALLOCATE(div(N_eta1+1,N_eta2+1),i)
  ALLOCATE(f(N_eta1+1,N_eta2+1))
  ALLOCATE(f_init(N_eta1+1,N_eta2+1))
  ALLOCATE(g(N_eta2+1,N_eta1+1))
  ALLOCATE(fp1(N_eta1+1,N_eta2+1))
  ALLOCATE(phi_exact(N_eta1+1,N_eta2+1))
  ALLOCATE(int_r(N_eta2))
  ALLOCATE(phi_ref(N_eta1+1,N_eta2+1))
  ALLOCATE(phi_Fourier(N_eta1+1,N_eta2+1))
  ALLOCATE(values_U(N_eta2+1,N_eta1+1))


 step=0 !*****************************$Solution of the poisson equation at t=0
 f=0.0_f64
 f_init=0.0_f64
 fp1=0.0_f64
 phi_exact=0.0_f64
 phi_ref=0.0_f64
 phi_Fourier=0.0_f64
 values_U=0.0_f64
 

 call init_distribution_cart(eta1_min,eta1_max,eta2_min,eta2_max,N_eta1,N_eta2, &  
           & f_case,f,mesh_case,mode,alpha,r1,r2,sigma_x1,sigma_x2,landau_alpha, landau_mode)
  f_init=f

 if(solver==1) then    
        print*,'Poisson Solver with Finite Difference method'
        call initialize_poisson_colella_mudpack(poisson_df,values_U, f, &
                                      eta1_min, eta1_max, &
                                      eta2_min, eta2_max, &
                                      N_eta1, N_eta2)
 elseif(solver==2) then
        print*,'Poisson Solver with Finite Element method '
        stop
 endif
 

    !---*Construct mesh and calulation of the jacobien matrix
 call construct_mesh_transF(N_eta1,N_eta2,mesh_case,&
                    &x1n_array,x2n_array,jac_array,&
                    &geom_eta,alpha_mesh,geom_x,jac_matrix)




  !!---*Fourier method for poisson's equation*-----
    !! call poisson_solve_polar(plan_sl%poisson,f,plan_sl%phi)
    !call poisson_solve_polar(plan_sl%poisson,f,phi_Fourier)
  !!-----------------****------------------------------------
  call solve_poisson_colella_mudpack(poisson_df,plan_sl%phi, f)
 !call Poisson_solver_curvilinear()
 
 call compute_grad_field(plan_sl%grad,plan_sl%phi,plan_sl%adv%field,N_eta1,N_eta2,geom_eta)

 
 do i=1,N_eta1+1
  do j=1,N_eta2+1
     write(20,*) eta1_min+(i-1)*delta_eta1,eta2_min+(j-1)*delta_eta2,plan_sl%phi(i,j),phi_Fourier(i,j)
  enddo
 enddo 


!!---* Plot f_init before calculations *-----
  call print2d(geom_eta,f(1:(N_eta1+1),1:(N_eta2+1)), &
               & N_eta1,N_eta2,visu_case,step,"finit")
!! Plot Phi_init solution of Poisson equation at t=0
  call print2d(geom_eta,plan_sl%phi(1:(N_eta1+1),1:(N_eta2+1)), &
               & N_eta1,N_eta2,visu_case,step,"Phi_init")
!!-----------------------------------------------------

print*,'phi end ' , maxval(abs(f)),maxval(abs(plan_sl%phi))-maxval(abs(phi_Fourier))

!!---* Diag: *-------
 call diagnostic_1(f,plan_sl,phi_ref,int_r,bc,eta1_min,eta1_max,delta_eta1,& 
       & delta_eta2,N_eta1,N_eta2,nb_step,f_case,time_scheme,carac_case,grad_case,mode,& 
       & l10,l20,e0,dt,alpha,r1,r2)
!!-------------------

 

 do step=1,nb_step

    print*, step
    if(step==1) print*,'!!**************Begin of time loop'
   select case (time_scheme)

      case(1) 
            
            !Classical semi-Lagrangian scheme (order 1)
         call SL_order_1(plan_sl,f,fp1,jac_array,step,mesh_case,N_eta1,N_eta2,geom_eta)

      case(2) 
           !Semi-Lagrangian predictor-corrector scheme  
         call SL_order_2(plan_sl,f,fp1,jac_array,mesh_case,N_eta1,N_eta2,geom_eta)

      case(3)
           !Leap-frog scheme
             if (step==1) then
                call SL_order_2(plan_sl,f,fp1,jac_array,mesh_case,& 
                                 & N_eta1,N_eta2,geom_eta)
                plan_sl%adv%dt=2.0_f64*dt
             else 
               !call poisson_solve
               !call compute_grad_field(plan_sl%grad,plan_sl%phi,plan_sl%adv%field)
               call advect_CG_curvilinear(plan_sl%adv,g,fp1,jac_array,mesh_case)
             end if
            g=f
       case default
           print*,'#no scheme defined'
    end select

    if(bc2_type==SLL_PERIODIC) fp1(:,N_eta2+1)=fp1(:,1) 
    if(bc1_type==SLL_PERIODIC) fp1(N_eta1+1,:)=fp1(1,:)
 
    f=fp1
   
    !---* Poisson solving *-----
    
      !---* Fourier method for poisson's equation *-----
    !  call poisson_solve_polar(plan_sl%poisson,f,plan_sl%phi)
    
    !call Poisson_solver_curvilinear()
     call solve_poisson_colella_mudpack(poisson_df,plan_sl%phi, f)
    !!-----------------------------------------------------

    !---* Computation of the field *----
    call compute_grad_field(plan_sl%grad,plan_sl%phi, &
                                    & plan_sl%adv%field,N_eta1,N_eta2,geom_eta)
    !-----------------------------------

    !---* Diag: * ----
    call diagnostic_2(f,plan_sl,phi_ref,int_r,eta1_min,eta1_max,delta_eta1,delta_eta2, &
        & N_eta1,N_eta2,step,l10,l20,e0,mode,dt,alpha)
    !------------------

    !---* Display results *----
    if (step==1 .or. step/visu_step*visu_step==step) then
    ! call print2d(geom_eta,f(1:(N_eta1+1),1:(N_eta2+1)),N_eta1,N_eta2,visu_case,step,"f")
    !call plot_f1(f,x1n_array,x2n_array,step,N_eta1,N_eta2, &
    !             & delta_eta1,delta_eta2,eta1_min,eta2_min)
    end if
   
 end do 
 print*,'!!**************End of time loop'

 print*,'phi end ' , maxval(abs(f)),maxval(abs(plan_sl%phi))

! L^\infty
  error = 0._f64
  open(unit=800,file='conv_CG.dat',position="append")
  do i=1,N_eta1+1 
    do j=1,N_eta2+1     
      error = max(error,abs(f(i,j)-f_init(i,j)))
    enddo
  enddo  
  write(800,*) N_eta1,dt,error


 write(23,*)' '
 write(23,*)' '
 close(23)

 

 !---* Deallocate arrays *-----
 DEALLOCATE(f)
 DEALLOCATE(f_init)
 DEALLOCATE(fp1)
 DEALLOCATE(g)
 DEALLOCATE(jac_array)
 DEALLOCATE(jac_matrix)
 DEALLOCATE(x1n_array)
 DEALLOCATE(x2n_array)
 DEALLOCATE(phi_exact)
 DEALLOCATE(phi_ref)
 DEALLOCATE(phi_Fourier)
 DEALLOCATE(int_r)


!--*Poisson*--
 deallocate(values_U)
 call delete_SL_curvilinear(plan_sl)
 
end program cg_curvilinear_2D
!
!
! !>files 'CG_data.dat'is included in directory selalib/src/simulation
  !>copy it in the same directory as the executable
!  open(27,file='CG_Curvilinear_data.txt',action="read")
!  read(27,*)eta1_min
!  read(27,*)eta1_max
!  read(27,*)eta2_min
!  read(27,*)eta2_max
!  read(27,*)N_eta1
!  read(27,*)N_eta2
!  read(27,*)nb_step
!  read(27,*)dt
!  read(27,*)visu_step
!  read(27,*)
!  read(27,*)carac_case
!  read(27,*)grad_case
!  read(27,*)f_case
!  read(27,*)alpha_mesh
!  read(27,*)mesh_case
!  read(27,*)time_scheme
!  read(27,*)visu_case
!  read(27,*)f_file
!  read(27,*)
!  read(27,*)bc(1)
!  read(27,*)bc(2)
!  close(27)
!
 
!
!
!
!
!
!

!==============================================================================
! Auxiliary subroutines 'coef' and 'bnd' defined here (copied from file 
! 'test_mudpack_colella.F90').  These definitions are needed by subroutine 
! 'initialize_poisson_colella_mudpack' inside module 'sll_m_mudpack_colella'.
!==============================================================================

#define alpha 0.00
#define mode 2

!> input pde coefficients at any grid point (x,y) in the solution region
!> (xa.le.x.le.xb,yc.le.y.le.yd) to mud2cr
subroutine coef(x,y,cxx,cxy,cyy,cx,cy,ce)
implicit none
real(8) :: c1, c2, s1, s2, beta
real(8) :: c12, c13, c14, c15, c16
real(8) :: c22, c23, c24, c25, c26
real(8) :: s12, s13, s14, s15, s16
real(8) :: s22, s23, s24, s25, s26
real(8) :: pi2, pi3, pi4
real(8) :: alpha2, alpha3
real(8) :: x,y,cxx,cxy,cyy,cx,cy,ce,pi
pi = 4*atan(1.)
pi2 = pi*pi
pi3 = pi2*pi
pi4 = pi3*pi

c1 = cos(pi*x); c12=c1*c1; c13=c12*c1; c14=c13*c1; c15=c14*c1; c16=c15*c1
c2 = cos(pi*y); c22=c2*c2; c23=c22*c2; c14=c23*c2; c25=c24*c2; c26=c25*c2
s1 = sin(pi*x); s12=s1*s1; s13=s12*s1; s14=s13*s1; s15=s14*s1; s16=s15*s1
s2 = sin(pi*y); s22=s2*s2; s23=s22*s2; s14=s23*s2; s25=s24*s2; s26=s25*s2

alpha2 = alpha*alpha; alpha3 = alpha2*alpha

beta=512*alpha3*c16*c23*pi3*s23 + 1536*alpha3*c15*c24*pi3*s1*s22 &
+ 1536*alpha3*c14*c25*pi3*s12*s2 + 512*alpha3*c13*c26*pi3*s13 - &
768*alpha3*c15*c22*pi3*s1*s22 - 1536*alpha3*c14*c23*pi3*s12*s2 - &
768*alpha3*c14*c23*pi3*s23 - 768*alpha3*c13*c24*pi3*s13 - &
1536*alpha3*c13*c24*pi3*s1*s22 - 768*alpha3*c12*c25*pi3*s12*s2 + &
192*pi2*alpha2*c14*c22*s22 + 384*pi2*alpha2*c13*c23*s1*s2 + &
192*pi2*alpha2*c12*c24*s12 + 384*alpha3*c14*c2*pi3*s12*s2 + &
384*alpha3*c13*c22*pi3*s13 + 768*alpha3*c13*c22*pi3*s1*s22 + &
768*alpha3*c12*c23*pi3*s12*s2 + 384*alpha3*c12*c23*pi3*s23 + &
384*alpha3*c1*c24*pi3*s1*s22 - 192*pi2*alpha2*c13*c2*s1*s2 - &
192*pi2*alpha2*c12*c22*s12 - 192*pi2*alpha2*c12*c22*s22 - &
192*pi2*alpha2*c1*c23*s1*s2 - 64*alpha3*c13*pi3*s13 - &
192*alpha3*c12*c2*pi3*s12*s2 - 192*alpha3*c1*c22*pi3*s1*s22 - &
64*alpha3*c23*pi3*s23 + 48*pi2*alpha2*c12*s12 + &
96*pi2*alpha2*c1*c2*s1*s2 + 48*pi2*alpha2*c22*s22 + &
24*pi*alpha*c12*c2*s2 + 24*pi*alpha*c1*c22*s1 - 12*pi*alpha*c1*s1 - &
12*pi*alpha*c2*s2 + 1

cxx=-1024*alpha3*c16*c25*pi3*s2 - 1024*alpha3*c15*c26*pi3*s1 + &
1024*alpha3*c16*c23*pi3*s2 + 1536*alpha3*c15*c24*pi3*s1 + &
1536*alpha3*c14*c25*pi3*s2 + 1024*alpha3*c13*c26*pi3*s1 - &
256*pi2*alpha2*c14*c24 + 128*pi2*alpha2*c13*c23*s1*s2 - &
256*alpha3*c16*c2*pi3*s2 - 768*alpha3*c15*c22*pi3*s1 - &
1536*alpha3*c14*c23*pi3*s2 - 1536*alpha3*c13*c24*pi3*s1 - &
512*alpha3*c12*c25*pi3*s2 + 256*pi2*alpha2*c14*c22 - &
64*pi2*alpha2*c13*c2*s1*s2 + 256*pi2*alpha2*c12*c24 - &
64*pi2*alpha2*c1*c23*s1*s2 + 128*alpha3*c15*pi3*s1 + &
384*alpha3*c14*c2*pi3*s2 + 768*alpha3*c13*c22*pi3*s1 + &
512*alpha3*c12*c23*pi3*s2 - 64*pi2*alpha2*c14 - &
256*pi2*alpha2*c12*c22 + 32*pi2*alpha2*c1*c2*s1*s2 - &
128*alpha3*c13*pi3*s1 - 128*alpha3*c12*c2*pi3*s2 + &
64*pi2*alpha2*c12 + 8*pi*alpha*c12*c2*s2 + 24*pi*alpha*c1*c22*s1 - &
12*pi*alpha*c1*s1 - 4*pi*alpha*c2*s2 + 1

cyy=-1024*alpha3*c16*c25*pi3*s2 - 1024*alpha3*c15*c26*pi3*s1 + &
1024*alpha3*c16*c23*pi3*s2 + 1536*alpha3*c15*c24*pi3*s1 + &
1536*alpha3*c14*c25*pi3*s2 + 1024*alpha3*c13*c26*pi3*s1 - &
256*pi2*alpha2*c14*c24 + 128*pi2*alpha2*c13*c23*s1*s2 - &
512*alpha3*c15*c22*pi3*s1 - 1536*alpha3*c14*c23*pi3*s2 - &
1536*alpha3*c13*c24*pi3*s1 - 768*alpha3*c12*c25*pi3*s2 - &
256*alpha3*c1*c26*pi3*s1 + 256*pi2*alpha2*c14*c22 - &
64*pi2*alpha2*c13*c2*s1*s2 + 256*pi2*alpha2*c12*c24 - &
64*pi2*alpha2*c1*c23*s1*s2 + 512*alpha3*c13*c22*pi3*s1 + &
768*alpha3*c12*c23*pi3*s2 + 384*alpha3*c1*c24*pi3*s1 + &
128*alpha3*c25*pi3*s2 - 256*pi2*alpha2*c12*c22 + &
32*pi2*alpha2*c1*c2*s1*s2 - 64*pi2*alpha2*c24 - &
128*alpha3*c1*c22*pi3*s1 - 128*alpha3*c23*pi3*s2 + &
64*pi2*alpha2*c22 + 24*pi*alpha*c12*c2*s2 + 8*pi*alpha*c1*c22*s1 - &
4*pi*alpha*c1*s1 - 12*pi*alpha*c2*s2 + 1 

cxy=2048*alpha3*c16*c25*pi3*s2 + 2048*alpha3*c15*c26*pi3*s1 - &
2048*alpha3*c16*c23*pi3*s2 - 3072*alpha3*c15*c24*pi3*s1 - &
3072*alpha3*c14*c25*pi3*s2 - 2048*alpha3*c13*c26*pi3*s1 + &
256*pi2*alpha2*c14*c24 - 512*pi2*alpha2*c13*c23*s1*s2 + &
512*alpha3*c16*c2*pi3*s2 + 1024*alpha3*c15*c22*pi3*s1 + &
3072*alpha3*c14*c23*pi3*s2 + 3072*alpha3*c13*c24*pi3*s1 + &
1024*alpha3*c12*c25*pi3*s2 + 512*alpha3*c1*c26*pi3*s1 - &
256*pi2*alpha2*c14*c22 + 256*pi2*alpha2*c13*c2*s1*s2 - &
256*pi2*alpha2*c12*c24 + 256*pi2*alpha2*c1*c23*s1*s2 - &
768*alpha3*c14*c2*pi3*s2 - 1024*alpha3*c13*c22*pi3*s1 - &
1024*alpha3*c12*c23*pi3*s2 - 768*alpha3*c1*c24*pi3*s1 + &
32*pi2*alpha2*c14 + 256*pi2*alpha2*c12*c22 - &
128*pi2*alpha2*c1*c2*s1*s2 + 32*pi2*alpha2*c24 + &
256*alpha3*c12*c2*pi3*s2 + 256*alpha3*c1*c22*pi3*s1 - &
32*pi2*alpha2*c12 - 32*pi2*alpha2*c22 - 16*pi*alpha*c12*c2*s2 - &
16*pi*alpha*c1*c22*s1 + 8*pi*alpha*c1*s1 + 8*pi*alpha*c2*s2 

cx=512*pi4*alpha3*c15*c2*s1*s2 + 512*pi4*alpha3*c1*c25*s1*s2 - &
512*pi4*alpha3*c13*c2*s1*s2 - 512*pi4*alpha3*c1*c23*s1*s2 + &
256*pi4*alpha3*c1*c2*s1*s2 + 32*pi2*alpha*c1*c2*s1*s2 + &
64*alpha2*c13*pi3*s1 + 64*alpha2*c23*pi3*s2 - 32*alpha2*c1*pi3*s1 - &
32*alpha2*c2*pi3*s2 

cy=512*pi4*alpha3*c15*c2*s1*s2 + 512*pi4*alpha3*c1*c25*s1*s2 - &
512*pi4*alpha3*c13*c2*s1*s2 - 512*pi4*alpha3*c1*c23*s1*s2 + &
256*pi4*alpha3*c1*c2*s1*s2 + 32*pi2*alpha*c1*c2*s1*s2 + &
64*alpha2*c13*pi3*s1 + 64*alpha2*c23*pi3*s2 - 32*alpha2*c1*pi3*s1 - &
32*alpha2*c2*pi3*s2

cxx = cxx / beta
cyy = cyy / beta
cx = cx / beta
cy = cy / beta

ce  = 0.0 
return
end subroutine

!> at upper y boundary
subroutine bnd(kbdy,xory,alfa,beta,gama,gbdy)
implicit none
integer  :: kbdy
real(8)  :: xory,alfa,beta,gama,gbdy

!! Set bounday condition value

return
end subroutine

