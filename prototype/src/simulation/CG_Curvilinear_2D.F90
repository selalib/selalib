program cg_curvilinear_2D
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use sll_timer
  use module_cg_curvi_structure
 !use curvilinear_2D_operators
  use curvilinear_2D_advection
  use module_cg_curvi_function
  !use poisson_curvilinear
  use numeric_constants
  implicit none

  type(sll_SL_curvilinear), pointer :: plan_sl
  type(time_mark), pointer :: t1,t2
  sll_real64, dimension (:,:), allocatable :: f,fp1,g,f_init !,phi_ref,div
  sll_int32 :: i,j,step,visu_step,hh,minn,ss
  sll_int32 :: N_eta1, N_eta2, nb_step
  sll_int32 :: f_case,carac_case,visu_case,phi_case
  sll_int32 :: bc(2),err
  sll_int32 :: PERIODIC_B,HERMITE_B
  sll_int32 :: bc1_type, bc2_type
  sll_int32 :: grad_case,mesh_case,time_scheme
  sll_real64 :: delta_eta1, delta_eta2,error2,error3,error4,error5
  sll_real64 :: eta1_min, eta1_max,eta2_min, eta2_max, dt, tf
  sll_real64 :: temps,a1,a2,error
  sll_real64 :: x1c,x2c,x1c_r,x2c_r,sigma_x1,sigma_x2
  sll_real64 :: geom_eta(2,2),geom_x(2,2),alpha_mesh
  sll_real64, dimension(:,:), pointer :: x1n_array,x2n_array,x1_tab,x2_tab  !,x1c_array,x2c_array
  sll_real64, dimension(:,:), pointer :: jac_array
  sll_real64, dimension(:,:,:), pointer :: jac_matrix
  sll_real64, dimension(:,:), pointer :: phi_exact
  character (len=16) :: f_file !,bctop,bcbot
  character (len=8) :: conv_CG

  namelist /param/ N_eta1,N_eta2,dt,nb_step,carac_case,mesh_case,visu_case,phi_case
  
  PERIODIC_B=0
  HERMITE_B=1

  t1 => new_time_mark()
  t2 => new_time_mark()

 
  visu_case = 1 ! 0 : gnuplot 
                ! 1 : vtk

  visu_step = 50

  mesh_case = 2 ! 1 : cartesian 
                ! 2 : polar 
                ! 3 : r^2 modified polar 
                ! 4 : colella 

  f_case = 4  ! 1 : constant function 
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
  
  a1 = 0.25_f64 !*0.01_f64
  a2 = 0.25_f64 !*0.01_f64
  bc1_type=PERIODIC_B
  bc2_type=PERIODIC_B
  N_eta1 = 100
  N_eta2 = 100
  dt = 0.01
  nb_step = 10
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

 
  read(*,NML=param)


  ! ---- * Construction of the mesh * ----
  
  ! mesh type : cartesian
  ! domain    : square [eta1_min eta1_max] x [eta2_min eta2_max]
  ! BC        : periodic-periodic
  if (mesh_case==1) then
    eta1_min = 0._f64
    eta1_max = 6._f64
    eta2_min = 0._f64
    eta2_max = 6._f64

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

    eta1_min = 0.5_f64
    eta1_max = 8._f64
    eta2_min = 0._f64
    eta2_max = 2._f64*sll_pi
   
    x1c_r=0._f64
    x2c_r=0._f64
   
    x1c = -3.5_f64
    x2c = 0._f64

    sigma_x1 = 0.2
    sigma_x2 = 0.2


    bc1_type = PERIODIC_B
    !bc1_type = HERMITE_B
    bc2_type = PERIODIC_B
  endif
  
  ! mesh type : Collela
  ! domain    : square [eta1_min eta1_max] x [eta2_min eta2_max]
  ! BC        : periodic-periodic
  if (mesh_case==4) then

    eta1_min = 0._f64
    eta1_max = 8._f64
    eta2_min = 0._f64
    eta2_max = 8._f64

    x1c_r=4._f64
    x2c_r=4._f64


    x1c = 2._f64
    x2c = 2._f64
    
    sigma_x1 = 0.15_f64
    sigma_x2 = 0.15_f64
    
    bc1_type = PERIODIC_B
    bc2_type = PERIODIC_B
    
    alpha_mesh = 0.1_f64
  endif

  
  
  

  geom_eta(1,1)=eta1_min
  geom_eta(2,1)=eta1_max
  geom_eta(1,2)=eta2_min
  geom_eta(2,2)=eta2_max

  delta_eta1= (eta1_max-eta1_min)/real(N_eta1,f64)
  delta_eta2= (eta2_max-eta2_min)/real(N_eta2,f64)
  
  print*,'# Pas en direction eta1: delta_eta1=',delta_eta1
  print*,'# Pas en direction eta2: delta_eta2=',delta_eta2
  print*,'# carac_case',carac_case
  print*,'# grad_case', grad_case
  print*,'# mesh_case', mesh_case
  ! mesh type
   if (mesh_case > 4) then
    print*,'Non existing case'
    print*,'mesh_case = 1:cartesian, 2:polar, 3:polar-like or 4:Collela or 5: Collela2'
    STOP
  endif
  ! distribution function
  if (f_case > 5) then
    print*,'Non existing case'
    print*,'test_case = 1 :f=1, 2 : cos(eta2), 3 :gaussian in eta1,'
    print*,'4 :centered gaussian in eta1 and eta2 depending on the mesh type or 5 :centered dirac'
    STOP
  endif
  ! advection field
  if (phi_case > 4) then
    print*,'Non existing case'
    print*,'field_case = 1 : translation, 2 :rotation, 3 :non homogeneous rotation'  
    STOP
  endif

  print*,'# alpha_mesh', alpha_mesh
  print*,'# bc1_type, bc2_type = 1:HERMITE, 0:PERIODIC', bc1_type,bc2_type
 


!!$  !definition of dt=tf/nb_step
  tf=dt*real(nb_step,f64)

  print*,'# nb_step =',nb_step
  print*,'# Pas de temps: dt =',dt
  print*,'# Temps final: tf =',tf


 call construct_mesh_transF(N_eta1,N_eta2,mesh_case,&
   &x1n_array,x2n_array,jac_array,&
   &geom_eta,alpha_mesh,geom_x,jac_matrix)
 

!********************************* 
plan_sl => new_SL(geom_eta,delta_eta1,delta_eta2,dt, &
                  & N_eta1,N_eta2,grad_case,carac_case,bc,bc1_type,bc2_type)


!  !SLL_ALLOCATE(div(N_eta1+1,N_eta2+1),i)
  SLL_ALLOCATE(f(N_eta1+1,N_eta2+1),err)
  SLL_ALLOCATE(f_init(N_eta1+1,N_eta2+1),err)
  SLL_ALLOCATE(g(N_eta2+1,N_eta1+1),err)
  SLL_ALLOCATE(fp1(N_eta1+1,N_eta2+1),err)
  SLL_ALLOCATE(phi_exact(N_eta1+1,N_eta2+1),err)
  SLL_ALLOCATE(x1_tab(N_eta1+1,N_eta2+1),err)
  SLL_ALLOCATE(x2_tab(N_eta1+1,N_eta2+1),err)
!
  step=0 !*****************************$Solution of the poisson equation at t=0
  f=0.0_f64
  f_init=0.0_f64
  fp1=0.0_f64
  phi_exact=0.0_f64

  
  call init_distribution_curvilinear(N_eta1,N_eta2,f_case,f,mesh_case,&
                                     & x1n_array,x2n_array,x1c,x2c,sigma_x1,sigma_x2)
  f_init=f

  call phi_analytique(phi_exact,plan_sl%adv,phi_case,x1n_array,x2n_array,a1,a2,x1c_r,x2c_r,jac_matrix)
!  call poisson_solve_curvilinear(plan_sl%poisson,f,plan_sl%phi)
  plan_sl%phi=phi_exact
  plan_sl%adv%field=0
  call compute_grad_field(geom_eta,plan_sl%grad,plan_sl%phi,plan_sl%adv%field,bc1_type,bc2_type,N_eta1,N_eta2)

!  !write f in a file before calculations
  call print2d(geom_eta,f(1:(N_eta1+1),1:(N_eta2+1)),N_eta1,N_eta2,visu_case,step,"finit")
! !***************************************


  t1 => start_time_mark(t1)
!
 do step=1,nb_step
    if(step==1) print*,'!!**************Begin of time loop'
    select case (time_scheme)

      case(1) 
            
            !classical semi-Lagrangian scheme (order 1)
            call carac_analytique(phi_case,N_eta1,N_eta2,x1n_array,x2n_array,a1,a2,x1c_r,x2c_r,&
                                  & x1_tab,x2_tab,dt)
            call SL_order_1(plan_sl,f,fp1,jac_array,step,x1_tab,x2_tab,mesh_case)
           
      case(2) 
           !semi-Lagrangian predictor-corrector scheme  
            call carac_analytique(phi_case,N_eta1,N_eta2,x1n_array,x2n_array,a1,a2,x1c_r,x2c_r,&
                                  & x1_tab,x2_tab,dt)
            call SL_order_2(plan_sl,f,fp1,jac_array,x1_tab,x2_tab,mesh_case)

      case(3)
           !leap-frog scheme
             if (step==1) then
                call SL_order_2(plan_sl,f,fp1,jac_array,x1_tab,x2_tab,mesh_case)
                plan_sl%adv%dt=2.0_f64*dt
             else 
               !call poisson_solve_curvilinear(plan_sl%poisson,f,plan_sl%phi)
               !call compute_grad_field(plan_sl%grad,plan_sl%phi,plan_sl%adv%field)
               call advect_CG_curvilinear(plan_sl%adv,g,fp1,jac_array,x1_tab,x2_tab,mesh_case)
             end if
            g=f
      case default
           print*,'#no scheme defined'
    end select

    if(bc2_type==PERIODIC_B) fp1(:,N_eta2+1)=fp1(:,1) 
    if(bc1_type==PERIODIC_B) fp1(N_eta1+1,:)=fp1(1,:)
    
    
    f=fp1
     

    
   
    if (step/visu_step*visu_step==step) then
       call print2d(geom_eta,f(1:(N_eta1+1),1:(N_eta2+1)),N_eta1,N_eta2,visu_case,step,"f")
    !call plot_f(f,x1n_array,x2n_array,step,N_eta1,N_eta2,delta_eta1,delta_eta2,eta1_min,eta2_min)
    end if
   
 end do 
 print*,'!!**************End of time loop'


 call carac_analytique(phi_case,N_eta1,N_eta2,x1n_array,x2n_array,a1,a2,x1c_r,x2c_r,&
                       &x1_tab,x2_tab,real(nb_step,f64)*dt)
!**************
! if(mesh_case==1) then 
!   do i=1,N_eta1+1
!    do j=1,N_eta1+1
!      call correction_BC(bc1_type,bc2_type,eta1_min,eta1_max,eta2_min,eta2_max,x1_tab(i,j),x2_tab(i,j))
!    end do
!   end do
! end if
!************

 call init_distribution_curvilinear(N_eta1,N_eta2,f_case,f_init,mesh_case,&
                                    &x1_tab,x2_tab,x1c,x2c,sigma_x1,sigma_x2)
    if(bc2_type==PERIODIC_B) f(:,N_eta2+1)=f(:,1) 
    if(bc1_type==PERIODIC_B) f(N_eta1+1,:)=f(1,:)
 !write f_init in a file after calculations
  call print2d(geom_eta,f_init(1:(N_eta1+1),1:(N_eta2+1)),N_eta1,N_eta2,visu_case,nb_step,"fexact")
  call print2d(geom_eta,&
   &f(2:(N_eta1),2:(N_eta2))-f_init(2:(N_eta1),2:(N_eta2)),N_eta1-2,N_eta2-2,visu_case,nb_step,"ferr")

! L^\infty
  error = 0._f64
  error2=0._f64
  error3=0._f64
  error4=0._f64
  error5=0._f64
  open(unit=800,file='conv_CG.dat',position="append")
  do i=1,N_eta1+1 !(N_eta1/4),(3*N_eta1/4)
    do j=1,N_eta2+1 !(N_eta1/4),(3*N_eta1/4)     
      error = max(error,abs(f(i,j)-f_init(i,j)))
    enddo
  enddo  
  write(800,*) N_eta1,dt,error


    do j=1,N_eta1+1     !5,N_eta2-3
      error2 = max(error2,abs(f(1,j)))
      error3 = max(error3,abs(f(N_eta1+1,j)))
    enddo
   do i=1,N_eta1+1     !5,N_eta2-3
      error4 = max(error4,abs(f(i,1)))
      error5 = max(error5,abs(f(i,N_eta1+1)))
    enddo

  write(160,*) error2,error3,error4,error5


 write(23,*)' '
 write(23,*)' '
 close(23)

 t2 => start_time_mark(t2)
 temps=time_elapsed_between(t1,t2)
 hh=floor(temps/3600.0d0)
 minn=floor((temps-3600.0d0*real(hh))/60.0d0)
 ss=floor(temps-3600.0d0*real(hh)-60.0d0*real(minn))
 print*,'# temps pour faire la boucle en temps : ',hh,'h',minn,'min',ss,'s'


 !write the final f in a file
 !call print2d(geom_eta,f(1:(N_eta1+1),1:(N_eta2+1)),N_eta1,N_eta2,visu_case,step,"f")
 !open (unit=21,file='CGCrestart.dat')
 !write(21,*)f
 !close(21)

 !SLL_DEALLOCATE_ARRAY(div,err)
 SLL_DEALLOCATE_ARRAY(f,err)
 SLL_DEALLOCATE_ARRAY(f_init,err)
 SLL_DEALLOCATE_ARRAY(fp1,err)
 SLL_DEALLOCATE_ARRAY(g,err)
 SLL_DEALLOCATE_ARRAY(jac_array,err)
 SLL_DEALLOCATE_ARRAY(jac_matrix,err)
 SLL_DEALLOCATE_ARRAY(x1n_array,err)
 SLL_DEALLOCATE_ARRAY(x2n_array,err)
 SLL_DEALLOCATE_ARRAY(x1_tab,err)
 SLL_DEALLOCATE_ARRAY(x2_tab,err)
 SLL_DEALLOCATE_ARRAY(phi_exact,err)
 
 t1 => delete_time_mark(t1)
 t2 => delete_time_mark(t2)
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
