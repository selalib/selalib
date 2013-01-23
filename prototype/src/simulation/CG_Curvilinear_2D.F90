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
  sll_int32 :: bc1_type, bc2_type
  sll_int32 :: grad_case,mesh_case,time_scheme
  sll_real64 :: delta_eta1, delta_eta2
  sll_real64 :: eta1_min, eta1_max,eta2_min, eta2_max, dt, tf
  sll_real64 :: temps,a1,a2,error
  sll_real64 :: geom_eta(2,2),alpha_mesh
  sll_real64, dimension(:,:), pointer :: x1n_array,x2n_array,eta1_tab,eta2_tab  !,x1c_array,x2c_array
  sll_real64, dimension(:,:), pointer :: jac_array
  sll_real64, dimension(:,:), pointer :: phi_exact
  character (len=16) :: f_file !,bctop,bcbot
  character (len=8) :: conv_CG


  t1 => new_time_mark()
  t2 => new_time_mark()

  !>files 'CG_data.dat'is included in directory selalib/prototype/src/simulation
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
  visu_case = 0
  visu_step = 5
  mesh_case = 1
  f_case = 4  
  grad_case = 1
  carac_case = 1
  phi_case = 1
  time_scheme = 1
  
  a1 = 1._f64*0.01_f64
  a2 = 0._f64*0.01_f64
  bc1_type=PERIODIC_SPLINE
  bc2_type=PERIODIC_SPLINE
  N_eta1 = 30
  N_eta2 = 30
  dt = 0.01
  nb_step = 200
  alpha_mesh = 0._f64
  eta1_min = 0._f64
  eta1_max = 1._f64
  eta2_min = 0._f64
  eta2_max = 1._f64



  ! ---- * Construction of the mesh * ----
  
  ! mesh type : cartesian
  ! domain    : square [eta1_min eta1_max] x [eta2_min eta2_max]
  ! BC        : periodic-periodic
  if (mesh_case==1) then
    eta1_min = 0._f64
    eta1_max = 2._f64*sll_pi
    eta2_min = 0._f64
    eta2_max = 2._f64*sll_pi

    bc1_type = PERIODIC_SPLINE
    bc2_type = PERIODIC_SPLINE
  endif
  
  ! mesh type : polar or polar-like
  ! domain    : disc of radius eta1_max with a hole of radius eta1_min
  ! BC        : hermite-periodic
  if ((mesh_case==2).or.(mesh_case==3)) then
    eta1_min = 0.2_f64
    eta1_max = 0.8_f64
    eta2_min = 0._f64
    eta2_max = 2._f64*sll_pi

    bc1_type = HERMITE_SPLINE
    bc2_type = PERIODIC_SPLINE
  endif
  
  ! mesh type : Collela
  ! domain    : square [eta1_min eta1_max] x [eta2_min eta2_max]
  ! BC        : periodic-periodic
  if (mesh_case==4) then
    eta1_min = 0._f64
    eta1_max = 2._f64*sll_pi
    eta2_min = 0._f64
    eta2_max = 2._f64*sll_pi
    
    bc1_type = PERIODIC_SPLINE
    bc2_type = PERIODIC_SPLINE
    
    alpha_mesh = 0.01_f64
  endif

  
  
  

  geom_eta(1,1)=eta1_min
  geom_eta(2,1)=eta1_max
  geom_eta(1,2)=eta2_min
  geom_eta(2,2)=eta2_max

  delta_eta1=real(eta1_max-eta1_min,f64)/real(N_eta1,f64)
  delta_eta2=real(eta2_max-eta2_min,f64)/real(N_eta2,f64)
  
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
   &geom_eta,alpha_mesh)
 

!********************************* 
plan_sl => new_SL(eta1_min,eta1_max,eta2_min,eta2_max,delta_eta1,delta_eta2,dt, &
  & N_eta1,N_eta2,grad_case,carac_case,bc,bc1_type,bc2_type)


!  !SLL_ALLOCATE(div(N_eta1+1,N_eta2+1),i)
  SLL_ALLOCATE(f(N_eta1+1,N_eta2+1),err)
  SLL_ALLOCATE(f_init(N_eta1+1,N_eta2+1),err)
  SLL_ALLOCATE(g(N_eta2+1,N_eta1+1),err)
  SLL_ALLOCATE(fp1(N_eta1+1,N_eta2+1),err)
  SLL_ALLOCATE(phi_exact(N_eta1+1,N_eta2+1),err)
  SLL_ALLOCATE(eta1_tab(N_eta1+1,N_eta2+1),err)
  SLL_ALLOCATE(eta2_tab(N_eta1+1,N_eta2+1),err)
!
  step=0 !*****************************$Solution of the poisson equation at t=0
  f=0.0_f64
  f_init=0.0_f64
  fp1=0.0_f64
  phi_exact=0.0_f64

  call init_distribution_cartesian(eta1_min,eta1_max,eta2_min,eta2_max,N_eta1,N_eta2, & 
       & delta_eta1,delta_eta2,f_case,f,mesh_case)
  f_init=f
  call phi_analytique(phi_exact,plan_sl%adv,phi_case,x1n_array,x2n_array,a1,a2)
!  call poisson_solve_curvilivear(plan_sl%poisson,f,plan_sl%phi)
!  call compute_grad_field(plan_sl%grad,plan_sl%phi,plan_sl%adv%field)
 

!  !write f in a file before calculations
  call print2d(geom_eta,f(1:(N_eta1+1),1:(N_eta2+1)),N_eta1,N_eta2,visu_case,step,"CG")
! !***************************************


  t1 => start_time_mark(t1)
!
 do step=1,nb_step

    !print*,'#Step:',step

    select case (time_scheme)

      case(1) 
            
            !classical semi-Lagrangian scheme (order 1)
            call SL_classic(plan_sl,f,fp1,jac_array,eta1_tab,eta2_tab)
           
      case(2) 
           !semi-Lagrangian predictive-corrective scheme
            call SL_ordre_2(plan_sl,f,fp1,jac_array,eta1_tab,eta2_tab)

      case(3)
           !leap-frog scheme
             if (step==1) then
                call SL_ordre_2(plan_sl,f,fp1,jac_array,eta1_tab,eta2_tab)
                plan_sl%adv%dt=2.0_f64*dt
             else 
               !call poisson_solve_curvilinear(plan_sl%poisson,f,plan_sl%phi)
               !call compute_grad_field(plan_sl%grad,plan_sl%phi,plan_sl%adv%field)
               call advect_CG_curvilinear(plan_sl%adv,g,fp1,jac_array,eta1_tab,eta2_tab)
             end if
            g=f
      case default
           print*,'#no scheme defined'
    end select

    if(bc2_type==PERIODIC_SPLINE) fp1(:,N_eta2+1)=fp1(:,1) !!!periodic
    if(bc1_type==PERIODIC_SPLINE) fp1(N_eta1+1,:)=fp1(1,:)
    
    f=fp1

 call init_distribution_curvilinear(eta1_min,eta1_max,eta2_min,eta2_max,N_eta1,N_eta2, & 
             & delta_eta1,delta_eta2,f_case,f,mesh_case,eta1_tab,eta2_tab)

    
    if (step==1 .or. step/visu_step*visu_step==step) then
       !!call print2d(dom,f(1:(nr+1),1:(ntheta+1)),Nr,Ntheta,visu_case,step,"CGC")
    ! call plot_f(step,N_eta1,N_eta2,delta_eta1,delta_eta2,eta1_min,eta2_min)
    end if
   
 end do        !!**************End of time loop

! L^\infty
  error = 0._f64
  open(unit=800,file='conv_CG',position="append")
  do i=1,N_eta1+1
    do j=1,N_eta2+1
      error = max(error,abs(f(i,j)-f_init(i,j)))
    enddo
  enddo  
  write(800,*) N_eta1,error

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
 call print2d(geom_eta,f(1:(N_eta1+1),1:(N_eta2+1)),N_eta1,N_eta2,visu_case,step,"CGC")
 open (unit=21,file='CGCrestart.dat')
 write(21,*)f
 close(21)

 !SLL_DEALLOCATE_ARRAY(div,err)
 SLL_DEALLOCATE_ARRAY(f,err)
 SLL_DEALLOCATE_ARRAY(f_init,err)
 SLL_DEALLOCATE_ARRAY(fp1,err)
 SLL_DEALLOCATE_ARRAY(g,err)
 SLL_DEALLOCATE_ARRAY(jac_array,err)
 SLL_DEALLOCATE_ARRAY(x1n_array,err)
 SLL_DEALLOCATE_ARRAY(x2n_array,err)
 SLL_DEALLOCATE_ARRAY(eta1_tab,err)
 SLL_DEALLOCATE_ARRAY(eta2_tab,err)
 t1 => delete_time_mark(t1)
 t2 => delete_time_mark(t2)
 call delete_SL_curvilinear(plan_sl)


end program cg_curvilinear_2D
!
!
!!!$  !scheme to compute caracteristics
!!!$  ! 1 : using explicit Euler method
!!!$  ! 2 : rotation, rotation speed = -1
!!!$  ! 3 : using symplectic Euler with linear interpolation
!!!$  ! 4 : using symplectic Verlet with linear interpolation
!!!$  ! 5 : using fixed point method
!!!$  ! 6 : using modified symplectic Euler
!!!$  ! 7 : using modified symplectic Verlet
!!!$  ! 8 : using modified fixed point
!!!$  carac_case=4
!!!$
!!!$  !scheme to compute gradient
!!!$  ! 1 : finit differences in r and theta
!!!$  ! 2 : fft in theta, finit differences in r
!!!$  ! 3 : splines in r and theta
!!!$  grad=3
!!!$
!!!$  !distribution function
!!!$  ! 1 : gaussian in r, constant in theta
!!!$  ! 2 : f(r,theta)=1_[r1,r2](r)*(1+cos(theta))
!!!$  ! 3 : test distribution for poisson solver
!!!$  ! 4 : (gaussian in r)*(1+cos(theta))
!!!$  ! 5 : read f in a file with syntax : r theta x y f(i,j)
!!!$  f_case=2
!!!$  !f_file='CGfinal04.dat' !not working
!!!$
!!!$  !choose the way to compute
!!!$  ! 1 : Semi-Lagrangian scheme order 1
!!!$  ! 2 : Semi-Lagrangian scheme order 2
!!!$  ! 3 : leap-frog scheme
!!!$  scheme=2
!!!$
!!!$  !choose the visualization
!!!$  ! 0 : gnuplot
!!!$  ! 1 : vtk
!!!$  visu=0
!
!
!
!
!
!
!
