program cg_curvilinear_2D
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "selalib.h"


   
  use module_cg_curvi_structure
  use curvilinear_2D_advection
  use module_cg_curvi_function
  use module_diagnostic
  !--*Poisson*----
  use sll_mudpack_colella
 ! use sll_general_coordinate_qn_solver
 ! use sll_logical_meshes
 ! use sll_common_coordinate_transformations
 ! use sll_module_scalar_field_2d_alternative
 ! use sll_constants
 ! use sll_arbitrary_degree_spline_interpolator_2d_module

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
  
  PERIODIC_B=0
  HERMITE_B=1

 
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
              ! 4 : gaussian in x and y

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
  bc1_type=PERIODIC_B
  bc2_type=PERIODIC_B
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
  x1c = 0.5_f64 ! center of gaussian 
  x2c = 0.5_f64
  sigma_x1 = 10._f64  ! variance for gaussian (exp(-(x-x1c)^2/(2*sigma_1^2)))
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
    
    bc1_type =PERIODIC_B
    bc2_type =PERIODIC_B
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

    bc1_type = HERMITE_B
    bc2_type = PERIODIC_B
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
    
    bc1_type = PERIODIC_B
    bc2_type = PERIODIC_B
    
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
    print*,'#test_case = 1 :f=1, 2 : cos(eta2), 3 :gaussian in eta1,'
    print*,'#4 :centered gaussian in eta1 and eta2 depending on the mesh type or 5 :centered dirac'
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

    if(bc2_type==PERIODIC_B) fp1(:,N_eta2+1)=fp1(:,1) 
    if(bc1_type==PERIODIC_B) fp1(N_eta1+1,:)=fp1(1,:)
 
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
! !>files 'CG_data.dat'is included in directory selalib/prototype/src/simulation
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
