program cg_csl_uniform
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
!contact: mehrenbe@math.unistra.fr for this  program

  use sll_constants
  use cubic_non_uniform_splines
  use cg_csl_uniform_module
  implicit none
  
  !parameters
  sll_int32 :: nb_diag,visu_step,scheme,time_case(2),poisson_case,test_case
  sll_int32 :: visu_case,modx,mody
  sll_int32 :: N_x1,N_x2,nb_step
  !sll_int32 :: nb_step,N_x1,N_x2,nc_eta1,nc_eta2
  sll_real64 :: eps,dt
  sll_real64 :: landau_alpha,landau_k
  !other variables
  sll_int32 :: N,err  
  sll_real64 :: x1_min,x1_max,x2_min,x2_max,x1,x2,delta_x1,delta_x2,L_x1,L_x2
  sll_int32 :: i,i1,i2
  type(cubic_nonunif_spline_1D), pointer :: spl_per_x1, spl_per_x2
  sll_real64 ::geom_x(2,2),geom_eta(2,2)
  !sll_real64,dimension(:), pointer :: rho,E,phi_poisson,buf1d
  sll_real64,dimension(:,:), pointer :: f,E_x1,E_x2,fold
  sll_real64,dimension(:),allocatable::bufpoisson,bufsplpoisson,bufphi,bufspl,thf0,thf,thff
  complex(f64),dimension(:,:,:),allocatable::bufcpoisson
  complex(f64),dimension(:,:),allocatable::bufcphi
  sll_real64,dimension(:,:),allocatable::buf2dspl,buff
  integer::Nbuf,Nbufpoisson,Nbufcpoisson(2),Nbufphi,Nbufcphi(2),Nbufsplpoisson
  integer::Nbuf2dspl(2),Nbufspl,step
  sll_real64,dimension(:), pointer :: Xstar,node_positions_x1,node_positions_x2,buf1d


  !sll_real64,dimension(:,:),pointer::f,f_init,f_store
  !sll_real64,dimension(:,:),pointer::x1n_array,x2n_array,x1c_array,x2c_array
  !sll_real64,dimension(:,:),pointer::jac_array
  !sll_real64, dimension(:,:), pointer :: a1,a2,psi
  !sll_real64,dimension(:,:,:),pointer::integration_points
  !sll_int32  :: i1,i2,ierr,i,step
  !sll_real64 :: delta_x1,delta_x2,x1,x2,x1c,x2c
  !sll_real64 :: eta1_min,eta1_max,eta2_min,eta2_max,delta_eta1,delta_eta2,eta1,eta1c,eta2,eta2c
  !sll_real64 :: val,tmp

  nb_diag=8
  test_case=1
  visu_step=10
  scheme=20
  !pause=0.00
  !visu=0
  !execu=1
  time_case=[3,20]
  eps=1.e10_f64
  poisson_case=2
  modx=1
  mody=1
  !scratch_file="~/scratch/visu/"
  visu_case = 0
  
  SLL_ALLOCATE(thf(nb_diag),err)
  SLL_ALLOCATE(thf0(nb_diag),err)
  SLL_ALLOCATE(thff(nb_diag),err)
  
  !initialization of test_case
  if((test_case==1).or.(test_case==2)) then!khp
    landau_k=0.5_f64
    landau_alpha=0.015_f64
    N_x1=128
    N_x2=128
    dt=0.1
    nb_step=500
    x1_min = 0._f64
    x1_max = 2._f64*sll_pi/landau_k
    x2_min = 0._f64 
    x2_max = 2._f64*sll_pi
  endif  

  if(.not.((test_case==1).or.(test_case==2)))then
    print *,'test_case=',test_case,' not implemented'
    stop
  endif
  
  N = max(N_x1,N_x2)  
    
  !SLL_ALLOCATE(f(N_x1+1,N_x2+1),err)  
  !SLL_ALLOCATE(E_x1(N_x1+1,N_x2+1),err)  
  !SLL_ALLOCATE(E_x2(N_x1+1,N_x2+1),err)
  !SLL_ALLOCATE(buff(N_x1+1,N_x2+1),err)
  
  SLL_ALLOCATE(f(N_x1+1,N_x2+1),err)  
  SLL_ALLOCATE(E_x1(N_x1+1,N_x2+1),err)  
  SLL_ALLOCATE(E_x2(N_x1+1,N_x2+1),err)
  SLL_ALLOCATE(buff(N_x1+1,N_x2+1),err)
  SLL_ALLOCATE(fold(N_x1+1,N_x2+1),err)  

  SLL_ALLOCATE(buf1d(N+1),err)
  SLL_ALLOCATE(Xstar(N+1),err)
  
  SLL_ALLOCATE(node_positions_x1(N_x1+1),err)
  SLL_ALLOCATE(node_positions_x2(N_x2+1),err)



  delta_x1 = (x1_max-x1_min)/real(N_x1,f64)
  delta_x2 = (x2_max-x2_min)/real(N_x2,f64)
  L_x1 = x1_max-x1_min
  L_x2 = x2_max-x2_min
  
  geom_x(1,1) = x1_min
  geom_x(2,1) = x1_max
  geom_x(1,2) = x2_min
  geom_x(2,2) = x2_max
  
  
  
 do i=1,N_x1+1
    node_positions_x1(i) = x1_min+(real(i,f64)-1._f64)*delta_x1
  enddo
  
  node_positions_x1 = (node_positions_x1-x1_min)/(x1_max-x1_min)
  
  do i=1,N_x2+1
    node_positions_x2(i) = x2_min+(real(i,f64)-1._f64)*delta_x2
  enddo

  node_positions_x2 = (node_positions_x2-x2_min)/(x2_max-x2_min)
 
  spl_per_x1 =>  new_cubic_nonunif_spline_1D( N_x1, PERIODIC_SPLINE)
  spl_per_x2 =>  new_cubic_nonunif_spline_1D( N_x2, PERIODIC_SPLINE)
  
  
  if(test_case==1)then !khp testcase 
    do i2=1,N_x2+1
      do i1=1,N_x1+1
        x1 = x1_min+real(i1-1,f64)*delta_x1
        x2 = x2_min+real(i2-1,f64)*delta_x2         
        f(i1,i2) = sin(x2)+landau_alpha*cos(landau_k*x1)
      enddo
    enddo  
  endif

  if(test_case==2)then  !khp2 testcase
    do i2=1,N_x2+1
      do i1=1,N_x1+1
        x1 = x1_min+real(i1-1,f64)*delta_x1
        x2 = x2_min+real(i2-1,f64)*delta_x2         
        f(i1,i2) = sin(x1)*sin(x2)+&
        landau_alpha*sin(2._f64*x2)*sin(2._f64*landau_k*x1)
      enddo
    enddo  
  endif

  step=0  
  call print2dper(geom_x,f(1:N_x1,1:N_x2),N_x1,N_x2,visu_case,step,'f')


  
  call poisson2dperalloc(bufpoisson,Nbufpoisson,bufcpoisson,Nbufcpoisson,N_x1,N_x2)
  call poisson2dperinit(bufpoisson,bufcpoisson,N_x1,N_x2)
  
  
  call advect2dalloc(buf2dspl,Nbuf2dspl,bufspl,Nbufspl,N_x1,N_x2)
  call advect2dinit(buf2dspl,bufspl,N_x1,N_x2)

  call computephipersize(Nbufphi,Nbufcphi,N_x1,N_x2)
  call computephiperalloc(bufphi,bufcphi,N_x1,N_x2)
  call computephiperinit(bufphi,bufcphi,N_x1,N_x2)

  call splpoissonperper2dsize(Nbufsplpoisson,N_x1,N_x2)
  call splpoissonperper2dalloc(bufsplpoisson,N_x1,N_x2)
  call splpoissonperper2dinit(bufsplpoisson,N_x1,N_x2)
  
  E_x1=f
  if(poisson_case==1)then
    call poisson2dper(bufpoisson,Nbufpoisson,bufcpoisson,Nbufcpoisson(1),&
      Nbufcpoisson(2),E_x1(1:N_x1,1:N_x2),E_x2(1:N_x1,1:N_x2),N_x1,N_x2,L_x1,L_x2)
  endif
  
  if(poisson_case==2)then
    call computephiper(bufphi,Nbufphi,bufcphi,Nbufcphi(1),Nbufcphi(2),&
      E_x1(1:N_x1,1:N_x2),N_x1,N_x2,L_x1,L_x2)
    call splpoissonperper2d(E_x1(1:N_x1,1:N_x2),E_x2(1:N_x1,1:N_x2),&
      bufsplpoisson,Nbufsplpoisson,N_x1,N_x2,L_x1,L_x2)
  endif

  call print2dper(geom_x,E_x1(1:N_x1,1:N_x2),N_x1,N_x2,visu_case,step,'E_x')
  call print2dper(geom_x,E_x2(1:N_x1,1:N_x2),N_x1,N_x2,visu_case,step,'E_y')
  
  
  buff=f
  call thdiagcgper(buff(1:N_x1,1:N_x2),E_x1(1:N_x1,1:N_x2),E_x2(1:N_x1,1:N_x2),&
    geom_x,N_x1,N_x2,nb_diag,thf0,bufpoisson,Nbufpoisson,bufcpoisson,&
    Nbufcpoisson(1),Nbufcpoisson(2),modx,mody)

  print *,'#',thf0
  
  !stop

  if(scheme==1)then
    SLL_ALLOCATE(fold(N_x1+1,N_x2+1),err)    !time loop
    do step=1,nb_step      
      E_x1=f
      if(poisson_case==1)then
        call poisson2dper(bufpoisson,Nbufpoisson,bufcpoisson,Nbufcpoisson(1),&
          Nbufcpoisson(2),E_x1(1:N_x1,1:N_x2),E_x2(1:N_x1,1:N_x2),N_x1,N_x2,L_x1,L_x2)
      endif  
      if(poisson_case==2)then
        call computephiper(bufphi,Nbufphi,bufcphi,Nbufcphi(1),Nbufcphi(2),&
          E_x1(1:N_x1,1:N_x2),N_x1,N_x2,L_x1,L_x2)
        call splpoissonperper2d(E_x1(1:N_x1,1:N_x2),E_x2(1:N_x1,1:N_x2),&
          bufsplpoisson,Nbufsplpoisson,N_x1,N_x2,L_x1,L_x2)
      endif
      fold=f
      call thdiagcgper(fold(1:N_x1,1:N_x2),E_x1(1:N_x1,1:N_x2),E_x2(1:N_x1,1:N_x2),&
        geom_x,N_x1,N_x2,nb_diag,thf,bufpoisson,Nbufpoisson,bufcpoisson,&
        Nbufcpoisson(1),Nbufcpoisson(2),modx,mody)
      thff=(thf(2:4)-thf0(2:4))/thf0(2:4)
      ! tn,mass,(l1-l10)/l10,(l2-l20)/l20,(enstrophy-enstrophy0)/enstrophy0,
      print *,(step-1)*dt,thf(1),thff(1),thff(2),thff(3),thf(5),thf(6)
      !call advect2d(dom,f,Ey,-Ex,bufspl,buf2dspl,dt,timecase,eps)
      call advect2d(geom_x,f(1:N_x1,1:N_x2),E_x2(1:N_x1,1:N_x2),-E_x1(1:N_x1,1:N_x2),&
        N_x1,N_x2,bufspl,Nbufspl,buf2dspl,dt,time_case,eps)
      
      if(modulo(step,visu_step)==0)then
        call print2dper(geom_x,f(1:N_x1,1:N_x2),N_x1,N_x2,visu_case,step,'f')
      endif
    enddo
  endif
  
  if(scheme==2)then
    SLL_ALLOCATE(fold(N_x1+1,N_x2+1),err)    !time loop
    fold=f
    do step=1,nb_step      
      E_x1=f
      if(poisson_case==1)then
        call poisson2dper(bufpoisson,Nbufpoisson,bufcpoisson,Nbufcpoisson(1),&
          Nbufcpoisson(2),E_x1(1:N_x1,1:N_x2),E_x2(1:N_x1,1:N_x2),N_x1,N_x2,L_x1,L_x2)
      endif  
      if(poisson_case==2)then
        call computephiper(bufphi,Nbufphi,bufcphi,Nbufcphi(1),Nbufcphi(2),&
          E_x1(1:N_x1,1:N_x2),N_x1,N_x2,L_x1,L_x2)
        call splpoissonperper2d(E_x1(1:N_x1,1:N_x2),E_x2(1:N_x1,1:N_x2),&
          bufsplpoisson,Nbufsplpoisson,N_x1,N_x2,L_x1,L_x2)
      endif
      fold=f
      call thdiagcgper(fold(1:N_x1,1:N_x2),E_x1(1:N_x1,1:N_x2),E_x2(1:N_x1,1:N_x2),&
        geom_x,N_x1,N_x2,nb_diag,thf,bufpoisson,Nbufpoisson,bufcpoisson,&
        Nbufcpoisson(1),Nbufcpoisson(2),modx,mody)
      thff=(thf(2:4)-thf0(2:4))/thf0(2:4)
      ! tn,mass,(l1-l10)/l10,(l2-l20)/l20,(enstrophy-enstrophy0)/enstrophy0,
      print *,(step-1)*dt,thf(1),thff(1),thff(2),thff(3),thf(5),thf(6)
        
      fold=f
      call advect2d(geom_x,f(1:N_x1,1:N_x2),E_x2(1:N_x1,1:N_x2),-E_x1(1:N_x1,1:N_x2)&
      ,N_x1,N_x2,bufspl,Nbufspl,buf2dspl,0.5_f64*dt,&
      time_case,eps)

      E_x1=f
      if(poisson_case==1)then
        call poisson2dper(bufpoisson,Nbufpoisson,bufcpoisson,Nbufcpoisson(1),&
          Nbufcpoisson(2),E_x1(1:N_x1,1:N_x2),E_x2(1:N_x1,1:N_x2),N_x1,N_x2,L_x1,L_x2)
      endif  
      if(poisson_case==2)then
        call computephiper(bufphi,Nbufphi,bufcphi,Nbufcphi(1),Nbufcphi(2),&
          E_x1(1:N_x1,1:N_x2),N_x1,N_x2,L_x1,L_x2)
        call splpoissonperper2d(E_x1(1:N_x1,1:N_x2),E_x2(1:N_x1,1:N_x2)&
          ,bufsplpoisson,Nbufsplpoisson,N_x1,N_x2,L_x1,L_x2)
      endif
      call advect2d(geom_x,fold(1:N_x1,1:N_x2),E_x2(1:N_x1,1:N_x2),-E_x1(1:N_x1,1:N_x2),&
      N_x1,N_x2,bufspl,Nbufspl,buf2dspl,dt,time_case,eps)
      f=fold
      if(modulo(step,visu_step)==0)then
        call print2dper(geom_x,f(1:N_x1,1:N_x2),N_x1,N_x2,visu_case,step,'f')
      endif
    enddo
  endif
  
  if(scheme==10)then !Euler scheme
    SLL_ALLOCATE(fold(N_x1+1,N_x2+1),err)    !time loop
    do step=1,nb_step            
      E_x1=f
      if(poisson_case==1)then
        call poisson2dper(bufpoisson,Nbufpoisson,bufcpoisson,Nbufcpoisson(1),&
          Nbufcpoisson(2),E_x1(1:N_x1,1:N_x2),E_x2(1:N_x1,1:N_x2),N_x1,N_x2,L_x1,L_x2)
      endif  
      if(poisson_case==2)then
        call computephiper(bufphi,Nbufphi,bufcphi,Nbufcphi(1),Nbufcphi(2),&
          E_x1(1:N_x1,1:N_x2),N_x1,N_x2,L_x1,L_x2)
        call splpoissonperper2d(E_x1(1:N_x1,1:N_x2),E_x2(1:N_x1,1:N_x2),&
          bufsplpoisson,Nbufsplpoisson,N_x1,N_x2,L_x1,L_x2)
      endif
      fold=f
      call thdiagcgper(fold(1:N_x1,1:N_x2),E_x1(1:N_x1,1:N_x2),E_x2(1:N_x1,1:N_x2),&
        geom_x,N_x1,N_x2,nb_diag,thf,bufpoisson,Nbufpoisson,bufcpoisson,&
        Nbufcpoisson(1),Nbufcpoisson(2),modx,mody)
      thff=(thf(2:4)-thf0(2:4))/thf0(2:4)
      ! tn,mass,(l1-l10)/l10,(l2-l20)/l20,(enstrophy-enstrophy0)/enstrophy0,
      print *,(step-1)*dt,thf(1),thff(1),thff(2),thff(3),thf(5),thf(6)
      !call advect2d(geom_x,f,E_x2,-E_x1,N_x1,N_x2,bufspl,Nbufspl,buf2dspl,dt,&
      !time_case,eps)
      E_x1=-E_x1 !warning -E_x(x_i,y_j) stored in E_x1(i,j)
      call advect_classical_csl(dt,E_x2,E_x1,f,geom_x,N_x1,N_x2,buf1d,&
        node_positions_x1,node_positions_x2,Xstar,spl_per_x1,spl_per_x2)
      E_x1=-E_x1 !now E_x(x_i,y_j) stored in E_x1(i,j)

      
            
      if(modulo(step,visu_step)==0)then
        call print2dper(geom_x,f(1:N_x1,1:N_x2),N_x1,N_x2,visu_case,step,'f')
      endif
    enddo


  
  
  
  endif

  if(scheme==20)then !predictor corrector scheme
    SLL_ALLOCATE(fold(N_x1+1,N_x2+1),err)    !time loop
    do step=1,nb_step      
      !compute poisson
      E_x1=f
      if(poisson_case==1)then
        call poisson2dper(bufpoisson,Nbufpoisson,bufcpoisson,Nbufcpoisson(1),&
          Nbufcpoisson(2),E_x1(1:N_x1,1:N_x2),E_x2(1:N_x1,1:N_x2),N_x1,N_x2,L_x1,L_x2)
      endif  
      if(poisson_case==2)then
        call computephiper(bufphi,Nbufphi,bufcphi,Nbufcphi(1),Nbufcphi(2),&
          E_x1(1:N_x1,1:N_x2),N_x1,N_x2,L_x1,L_x2)
        call splpoissonperper2d(E_x1(1:N_x1,1:N_x2),E_x2(1:N_x1,1:N_x2),&
          bufsplpoisson,Nbufsplpoisson,N_x1,N_x2,L_x1,L_x2)
      endif
      fold=f
      !diagnostics
      call thdiagcgper(fold(1:N_x1,1:N_x2),E_x1(1:N_x1,1:N_x2),E_x2(1:N_x1,1:N_x2),&
        geom_x,N_x1,N_x2,nb_diag,thf,bufpoisson,Nbufpoisson,bufcpoisson,&
        Nbufcpoisson(1),Nbufcpoisson(2),modx,mody)
      thff=(thf(2:4)-thf0(2:4))/thf0(2:4)
      ! tn,mass,(l1-l10)/l10,(l2-l20)/l20,(enstrophy-enstrophy0)/enstrophy0,
      print *,(step-1)*dt,thf(1),thff(1),thff(2),thff(3),thf(5),thf(6)
      !call advect2d(geom_x,f,E_x2,-E_x1,N_x1,N_x2,bufspl,Nbufspl,buf2dspl,dt,&
      !time_case,eps)
      fold=f
      !compute advection of dt/2 with fold
      E_x1=-E_x1 !warning -E_x(x_i,y_j) stored in E_x1(i,j)
      call advect_classical_csl(0.5_f64*dt,E_x2,E_x1,f,geom_x,N_x1,N_x2,buf1d,&
        node_positions_x1,node_positions_x2,Xstar,spl_per_x1,spl_per_x2)
      E_x1=-E_x1 !now E_x(x_i,y_j) stored in E_x1(i,j)

      !compute field at time t_{n+1/2}
      E_x1=fold
      if(poisson_case==1)then
        call poisson2dper(bufpoisson,Nbufpoisson,bufcpoisson,Nbufcpoisson(1),&
          Nbufcpoisson(2),E_x1(1:N_x1,1:N_x2),E_x2(1:N_x1,1:N_x2),N_x1,N_x2,L_x1,L_x2)
      endif  
      if(poisson_case==2)then
        call computephiper(bufphi,Nbufphi,bufcphi,Nbufcphi(1),Nbufcphi(2),&
          E_x1(1:N_x1,1:N_x2),N_x1,N_x2,L_x1,L_x2)
        call splpoissonperper2d(E_x1(1:N_x1,1:N_x2),E_x2(1:N_x1,1:N_x2),&
          bufsplpoisson,Nbufsplpoisson,N_x1,N_x2,L_x1,L_x2)
      endif

      
      !compute advection of dt with f and E_{n+1/2}
      E_x1=-E_x1 !warning -E_x(x_i,y_j) stored in E_x1(i,j)
      call advect_classical_csl(dt,E_x2,E_x1,fold,geom_x,N_x1,N_x2,buf1d,&
        node_positions_x1,node_positions_x2,Xstar,spl_per_x1,spl_per_x2)
      E_x1=-E_x1 !now E_x(x_i,y_j) stored in E_x1(i,j)
      
      f=fold      
      if(modulo(step,visu_step)==0)then
        call print2dper(geom_x,f(1:N_x1,1:N_x2),N_x1,N_x2,visu_case,step,'f')
      endif
    enddo


  
  
  
  endif


  
  
  

  print *,"#End of program"
  

end program


