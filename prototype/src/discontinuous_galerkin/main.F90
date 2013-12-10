program VP_DG
#include "sll_working_precision.h"
#include "sll_boundary_condition_descriptors.h"

!from Madaule Eric :
!Warning : this code and all the modules it calls and I wrote are not optimized in
!term of computation time. \n
!We kept the readability of the code over computation performances

  use sll_nu_cart_mesh
  use sll_gausslobatto
  use sll_poisson4dg_1d_periodic_dg
  use sll_constants
  use sll_timestep4dg
  !this is part of FEMilaro
  use mod_sparse

  !use mod_octave_io
  !use mod_octave_io_sparse

  implicit none

  !plan for time step with DG
  type(dg_time_steping) :: dg_plan
  !distribution function at time n and n+1
  sll_real64,dimension(:,:),allocatable :: dist,distp1
  !number of time steb, number of cells un direction x and v,time step,
  !time history, snapshot of distribution and field, information to programmer
  sll_int32 :: nb_step,nx,nv,step,th,th_large,th_out
  !boudary in direction x and v
  sll_real64 :: x_min,x_max,v_min,v_max
  !mesh
  type(non_unif_cart_mesh) :: mesh
  !coefficients for fluxes
  sll_real64 :: c11,c12
  !number of Gauss-Lobatto points
  sll_int32 :: ng
  !time step and finalt time
  sll_real64 :: dt,tf
  !Gauss-Lobatto
  type(gausslobatto1D) :: gausslob!,gll
  !elements of Phi equal to 0 so Phi(x=0)=0
  sll_int32 :: xbound
  !space variables
  sll_real64 :: x,v,k,vv
  !energy and momentum to check conservation
  sll_real64 :: momentum,l1_f,l2_f,int_f
  sll_real64 :: k_en,em_en,energy,phi_jump
  !indexes for loops
  sll_int32 :: x1,x2,v1,v2
  !display parameter
  sll_int32 :: len,i,j
  sll_int32,dimension(:),allocatable :: lfor
  character(len=25) :: form,fdist,ffield
  logical :: sym_mom_comp

  ! for the python script *.py
!!$  namelist /param/ nx,nv,ng
!!$  read(*,NML=param)

  sym_mom_comp=.true. ! way to compute momentum
                       ! true for symetric computation (0->nv/2 and nv->nv/2)
                       ! false for linear computation (0->nv)

  !definition of geometry and data
  !k=0.5d0
  !k=2.0d0/13.0d0
  k=0.3d0
  !k=1.0d0/k

!!$  x_min=0.0d0
!!$  x_max=2.0d0*sll_pi ! this is to test the Poisson sover
!!$  v_min=0
!!$  v_max=sll_pi

!!$  x_min=0.0d0
!!$  x_max=4.0d0*sll_pi
!!$  v_min=-2.0d0
!!$  v_max=2.0d0

  !!!>>>Blanca's case
!!$  x_min=-sll_pi
!!$  x_max=sll_pi
!!$  v_min=-4
!!$  v_max=4
  !!!<<<

!!$  x_min=0.0d0
!!$  x_max=1.0d0
!!$  v_min=-1.0d0
!!$  v_max=1.0d0

!!$  x_min=-1.0d0
!!$  x_max=1.0d0
!!$  v_min=-1.0d0
!!$  v_max=1.0d0

  x_min=0.0d0
  !x_max=2.0d0*sll_pi/k
  x_max=20.0d0*sll_pi
  v_min=-8.0d0
  v_max=8.0d0

  nx=15
  nv=15
  ng=8

  print*,'discretization caracteristics :'
  print"(3(a5,i3))",'nx=',nx,', nv=',nv,', ng=',ng

  !xbound=ng*nx/2
  xbound=1

  allocate(dist(nx*ng,nv*ng),distp1(nx*ng,nv*ng))

  !definition or time step, delta_t and final time
  dt=0.001d0
  tf=1000.0d0
  nb_step=ceiling(tf/dt)
  th=min(20,int(0.1d0/dt))!20
  th_out=int(0.5d0/dt)
  th_large=int(10.0d0/dt) ! must be a multiple of th
  !th_large=huge(1)
  if (modulo(th_large,th)/=0) then
     print*,"WARNING : snapshot of the electric field may be miscalculated"
  end if
  tf=dt*nb_step
  print*,'size of time step      :',dt
  print*,'number of time step    :',nb_step
  print*,'final time             :',tf
  print*,'number of time history :',nb_step/th+1
  print*,'number of f/E snapshot :',nb_step/th_large+1

  len=1
  do while (nb_step>=10**len)
     len=len+1
  end do
  allocate(lfor(len))
  lfor=0
  form="(1x,i"//char(len+48)//",a1,i"//char(len+48)//",1x,a11,f7.2)"

  call init_gausslobatto_1d(ng,gausslob)
!!$print*,ng
!!$print*,gausslob%weigh
!!$print*,gausslob%node
!!$stop
  !call init_gausslobatto_1d(ng*3,gll)
  call init_nu_cart_mesh(nx,nv,mesh)
  mesh%d_etat1=(x_max-x_min)/real(nx,8)
  mesh%d_etat2=(v_max-v_min)/real(nv,8)
  !mesh%d_etat2=1.0d0
  !mesh%d_etat2(9:32)=0.25d0
  call fill_node_nuc_mesh(x_min,v_min,mesh)

!!$  open(12)
!!$  do v1=1,nv
!!$     do x1=1,nx
!!$        write(12,*)mesh%etat2(v1),mesh%etat1(x1),1.0
!!$     end do
!!$     write(12,*)""
!!$  end do
!!$  close(12)
!!$  stop

  !flux coefficients
  c12=0.5d0
  c11=real((ng-1)**2,8)/maxval(mesh%d_etat1)

  call init_timesteping_4dg(dg_plan,sll_rk4,gausslob,mesh,dt,xbound,c11,c12,alpha=0.0d0)

  !construction of distribution function
  !x_i is indexed on both mesh nodes and GLL nodes, so to have the postion in x
  !one shoul take the lower node of the element and add the part due to GLL
  !the 1.0d0 is to compensate the fact that GLL is done on [-1;1]
  !same is done for v

  !test distribution for the Poisson's problem : 
  !f(x,v)=sin(x)sin(v), (x,v)\in [0;pi]^2 (to check rhs, independants on t)
  !exact rhs = -v*cos(x)sin(x)+sin(2x)cos(v)
  dist=0.0d0
  write(fdist,*)'dist_'
  write(ffield,*)'field_'
  do i=1,len
     write(fdist,*)trim(adjustl(fdist))//char(48)
     write(ffield,*)trim(adjustl(ffield))//char(48)
     fdist=trim(adjustl(fdist))
     ffield=trim(adjustl(ffield))
  end do
  open(14,file=fdist)
  do v1=1,nv
     do v2=1,ng
        v=mesh%etat2(v1)+(1.0d0+gausslob%node(v2))/mesh%jac(nx+1,v1)
        do x1=1,nx
           do x2=1,ng
              x=mesh%etat1(x1)+(1.0d0+gausslob%node(x2))/mesh%jac(x1,nv+1)
              !test case for RHS
!!$              dist((x1-1)*ng+x2,(v1-1)*ng+v2)=sin(x)*sin(v)
!!$              dist((x1-1)*ng+x2,(v1-1)*ng+v2)=0.5d0*(+1.0d0+sin(x))*sin(v)

!!$              dist((x1-1)*ng+x2,(v1-1)*ng+v2)=(exp(-200.0d0*(v-0.8d0)**2)+ &
!!$                   & exp(-200.0d0*(v+0.8d0)**2))!*(cos(3.0d0*x)+cos(6.0d0*x)+cos(18.0d0*x))

              !!!>>>Blanca's test case
!!$              dist((x1-1)*ng+x2,(v1-1)*ng+v2)=(2.0d0-cos(2.0d0*x))* &
!!$                   & exp(-0.25d0*(4.0d0*v-1.0d0)**2)
              !!!<<<

!!$              if (abs(v)<=0.5d0) then
!!$                 dist((x1-1)*ng+x2,(v1-1)*ng+v2)=1.0d0* &
!!$                      & sin(x*sll_pi)
!!$              end if

!!$              dist((x1-1)*ng+x2,(v1-1)*ng+v2)=(x**2-1.0d0)*(v**2-1.0d0)

              !Landau
!!$              dist((x1-1)*ng+x2,(v1-1)*ng+v2)=(1.0d0-0.5d0*cos(k*x))* &
!!$                   & exp(-v**2/2.0d0)/sqrt(2.0d0*sll_pi)

              !strong oscillations two streams
!!$              dist((x1-1)*ng+x2,(v1-1)*ng+v2)=v**2/sqrt(8.0d0*sll_pi)* &
!!$                   & (2.0d0-cos(k*(x-2.0d0*sll_pi)))* &
!!$                   & exp(-v**2/2.0d0)/sqrt(2.0d0*sll_pi)

              !classical two streams instability
!!$              dist((x1-1)*ng+x2,(v1-1)*ng+v2)=(1.0d0+0.05d0*cos(k*x))/ &
!!$                   & (2.0d0*0.3d0*sqrt(2.0d0*sll_pi))*( &
!!$                   & exp(-(v-0.99d0)**2/(2.0d0*0.3d0**2))+ &
!!$                   & exp(-(v+0.99d0)**2/(2.0d0*0.3d0**2)))

              !Bump on tail
              !from michel, eric and nicolas, inria
              dist((x1-1)*ng+x2,(v1-1)*ng+v2)=(1.0d0+0.04d0*cos(k*x))/ &
                   & (10.0d0*sqrt(2.0d0*sll_pi))* &
                   & (9.0d0*exp(-v**2/2.0d0)+2.0d0*exp(-(v-4.5d0)**2/(2.0d0*0.5d0**2)))

!!$              dist((x1-1)*ng+x2,(v1-1)*ng+v2)=&!(1.0d0-0.05d0*cos(k*x))
!!$                   & 1.0d0/sqrt(2.0d0*sll_pi)* &
!!$                   & (exp(-v**2/2.0d0)+0.2d0*exp(-(v-2.0d0)**2*100))

!!$              dist((x1-1)*ng+x2,(v1-1)*ng+v2)=exp(-v**2)/sqrt(2.0d0*sll_pi)

              write(14,*)x,v,dist((x1-1)*ng+x2,(v1-1)*ng+v2)
           end do
        end do
        write(14,*)''
     end do
  end do
  close(14)

  if (dg_plan%method==sll_SSP_RK) then
     dg_plan%max0=maxval(dist)
     !prevent negative value in case the initial distribution would be negative at some point
     dg_plan%min0=max(min(minval(dist),0.0d0),0.0d0)
  end if
  if (minval(dist)<0.0d0) then
     print*,"WARNING : initial distribution is not positive"
  end if

  call param_out(x_min,x_max,v_min,v_max,nx,nv,ng,.true.,c11,c12,tf,dt,nb_step,th,th_large)

  dg_plan%t=0.0d0
  open(15,file='energy_momentum',action='write')
  write(15,*)"# t ; momentum ; total energy ; kinetic energy ; electromagnetic energy ;", &
       & " jump of phi ; ||f||_L1 ; ||f||_L2 ; min(f) ; max(f) ; integral(f)"
  momentum=0.0d0
  energy=0.0d0
  k_en=0.0d0
  em_en=0.0d0
  phi_jump=0.0d0
  l1_f=0.0d0
  l2_f=0.0d0
  int_f=0.0d0

  dg_plan%rho=0.0d0
!!$  dg_plan%norma=0.0d0
  dg_plan%bound=0
  do x1=1,nx
     do x2=1,ng
        do v1=1,nv
           do v2=1,ng
              dg_plan%rho((x1-1)*ng+x2)=dg_plan%rho((x1-1)*ng+x2)+ &
                   & dist((x1-1)*ng+x2,(v1-1)*ng+v2)*gausslob%weigh(v2)/mesh%jac(nx+1,v1)
           end do
        end do
!!$        dg_plan%norma=dg_plan%norma+dg_plan%rho((x1-1)*ng+x2)* &
!!$             & gausslob%weigh(x2)/mesh%jac(x1,nv+1)
     end do
  end do
!!$  dg_plan%norma=dg_plan%norma/(x_max-x_min)
  dg_plan%norma=1.0d0
  do x1=1,nx
     do x2=1,ng
        dg_plan%rho((x1-1)*ng+x2)=(dg_plan%rho((x1-1)*ng+x2)-dg_plan%norma)* &
             & gausslob%weigh(x2)/mesh%jac(x1,nv+1)
     end do
  end do
  dg_plan%rho(dg_plan%x0)=dg_plan%bound ! boudary condition = \Phi(0,t) = alpha

  dg_plan%phi=0.0d0
  dg_plan%field=0.0d0
  call solve(dg_plan%poisson_vp,dg_plan%rho,dg_plan%phi)
  dg_plan%field=matmul(dg_plan%poisson_vp%mat_field,dg_plan%phi)

!!$!!!to check Poisson solver
!!$l2_f=0.0d0 ! error on phi
!!$l1_f=0.0d0 ! error on E
!!$do x1=1,nx
!!$   do x2=1,ng
!!$      x=mesh%etat1(x1)+(1.0d0+gausslob%node(x2))/mesh%jac(x1,nv+1)
!!$      l2_f=l2_f+(sin(x)-  dg_plan%phi((x1-1)*ng+x2))**2*gausslob%weigh(x2)/mesh%jac(x1,nv+1)
!!$      l1_f=l1_f+(cos(x)-dg_plan%field((x1-1)*ng+x2))**2*gausslob%weigh(x2)/mesh%jac(x1,nv+1)
!!$   end do
!!$end do
!!$l2_f=sqrt(l2_f)
!!$l1_f=sqrt(l1_f)
!!$open(42,file="check_poisson",position="append")
!!$write(42,*) mesh%d_etat1(1),l2_f,l1_f
!!$close(42)
!!$!stop

  if (sym_mom_comp) then
     if (modulo(nv,2)==0) then !nv even
        do v1=1,nv/2
           do v2=1,ng
              v=mesh%etat2(v1)+(1.0d0+gausslob%node(v2))/mesh%jac(nx+1,v1)
              vv=mesh%etat2(nv+1-v1)+(1.0d0+gausslob%node(ng+1-v2))/mesh%jac(nx+1,nv+1-v1)
              do x1=1,nx
                 do x2=1,ng
                    momentum=momentum+v*dist((x1-1)*ng+x2,(v1-1)*ng+v2)* &
                         & gausslob%weigh(x2)*gausslob%weigh(v2)/mesh%jac(x1,v1)+ &
                         & vv*dist((x1-1)*ng+x2,(nv-v1)*ng+ng+1-v2)* &
                         & gausslob%weigh(x2)*gausslob%weigh(ng+1-v2)/mesh%jac(x1,nv+1-v1)
                    !kinetik energy
                    k_en=k_en+v**2*dist((x1-1)*ng+x2,(v1-1)*ng+v2)* &
                         & gausslob%weigh(x2)*gausslob%weigh(v2)/mesh%jac(x1,v1)+ &
                         & vv**2*dist((x1-1)*ng+x2,(nv-v1)*ng+ng+1-v2)* &
                         & gausslob%weigh(x2)*gausslob%weigh(ng+1-v2)/mesh%jac(x1,nv+1-v1)

                    !norms 1 and 2 of distribution
                    l1_f=l1_f+abs(dist((x1-1)*ng+x2,(v1-1)*ng+v2))*gausslob%weigh(x2)* &
                         & gausslob%weigh(v2)/mesh%jac(x1,v1)+ &
                         & abs(dist((x1-1)*ng+x2,(nv-v1)*ng+ng+1-v2))* &
                         & gausslob%weigh(x2)*gausslob%weigh(ng+1-v2)/mesh%jac(x1,nv+1-v1)
                    l2_f=l2_f+dist((x1-1)*ng+x2,(v1-1)*ng+v2)**2*gausslob%weigh(x2)* &
                         & gausslob%weigh(v2)/mesh%jac(x1,v1)+ &
                         & dist((x1-1)*ng+x2,(nv-v1)*ng+ng+1-v2)**2* &
                         & gausslob%weigh(x2)*gausslob%weigh(ng+1-v2)/mesh%jac(x1,nv+1-v1)
                    int_f=int_f+dist((x1-1)*ng+x2,(v1-1)*ng+v2)*gausslob%weigh(x2)* &
                         & gausslob%weigh(v2)/mesh%jac(x1,v1)+ &
                         & dist((x1-1)*ng+x2,(nv-v1)*ng+ng+1-v2)* &
                         & gausslob%weigh(x2)*gausslob%weigh(ng+1-v2)/mesh%jac(x1,nv+1-v1)
                 end do
              end do
           end do
        end do
     else !nv ode
        do v1=1,floor(nv/2.)-1
           do v2=1,ng
              v=mesh%etat2(v1)+(1.0d0+gausslob%node(v2))/mesh%jac(nx+1,v1)
              vv=mesh%etat2(nv+1-v1)+(1.0d0+gausslob%node(ng+1-v2))/mesh%jac(nx+1,nv+1-v1)
              do x1=1,nx
                 do x2=1,ng
                    momentum=momentum+v*dist((x1-1)*ng+x2,(v1-1)*ng+v2)* &
                         & gausslob%weigh(x2)*gausslob%weigh(v2)/mesh%jac(x1,v1)+ &
                         & vv*dist((x1-1)*ng+x2,(nv-v1)*ng+ng+1-v2)* &
                         & gausslob%weigh(x2)*gausslob%weigh(ng+1-v2)/mesh%jac(x1,nv+1-v1)
                    !kinetik energy
                    k_en=k_en+v**2*dist((x1-1)*ng+x2,(v1-1)*ng+v2)* &
                         & gausslob%weigh(x2)*gausslob%weigh(v2)/mesh%jac(x1,v1)+ &
                         & vv**2*dist((x1-1)*ng+x2,(nv-v1)*ng+ng+1-v2)* &
                         & gausslob%weigh(x2)*gausslob%weigh(ng+1-v2)/mesh%jac(x1,nv+1-v1)

                    !norms 1 and 2 of distribution
                    l1_f=l1_f+abs(dist((x1-1)*ng+x2,(v1-1)*ng+v2))*gausslob%weigh(x2)* &
                         & gausslob%weigh(v2)/mesh%jac(x1,v1)+ &
                         & abs(dist((x1-1)*ng+x2,(nv-v1)*ng+ng+1-v2))* &
                         & gausslob%weigh(x2)*gausslob%weigh(ng+1-v2)/mesh%jac(x1,nv+1-v1)
                    l2_f=l2_f+dist((x1-1)*ng+x2,(v1-1)*ng+v2)**2*gausslob%weigh(x2)* &
                         & gausslob%weigh(v2)/mesh%jac(x1,v1)+ &
                         & dist((x1-1)*ng+x2,(nv-v1)*ng+ng+1-v2)**2* &
                         & gausslob%weigh(x2)*gausslob%weigh(ng+1-v2)/mesh%jac(x1,nv+1-v1)
                    int_f=int_f+dist((x1-1)*ng+x2,(v1-1)*ng+v2)*gausslob%weigh(x2)* &
                         & gausslob%weigh(v2)/mesh%jac(x1,v1)+ &
                         & dist((x1-1)*ng+x2,(nv-v1)*ng+ng+1-v2)* &
                         & gausslob%weigh(x2)*gausslob%weigh(ng+1-v2)/mesh%jac(x1,nv+1-v1)
                 end do
              end do
           end do
        end do
        v1=floor(nv/2.)
        do v2=1,ng
           v=mesh%etat2(v1)+(1.0d0+gausslob%node(v2))/mesh%jac(nx+1,v1)
           do x1=1,nx
              do x2=1,ng
                 momentum=momentum+v*dist((x1-1)*ng+x2,(v1-1)*ng+v2)* &
                      & gausslob%weigh(x2)*gausslob%weigh(v2)/mesh%jac(x1,v1)
                 !kinetik energy
                 k_en=k_en+v**2*dist((x1-1)*ng+x2,(v1-1)*ng+v2)* &
                      & gausslob%weigh(x2)*gausslob%weigh(v2)/mesh%jac(x1,v1)

                 !norms 1 and 2 of distribution
                 l1_f=l1_f+abs(dist((x1-1)*ng+x2,(v1-1)*ng+v2))*gausslob%weigh(x2)* &
                      & gausslob%weigh(v2)/mesh%jac(x1,v1)
                 l2_f=l2_f+dist((x1-1)*ng+x2,(v1-1)*ng+v2)**2*gausslob%weigh(x2)* &
                      & gausslob%weigh(v2)/mesh%jac(x1,v1)
                 int_f=int_f+dist((x1-1)*ng+x2,(v1-1)*ng+v2)*gausslob%weigh(x2)* &
                      & gausslob%weigh(v2)/mesh%jac(x1,v1)
              end do
           end do
        end do
     end if
     l2_f=sqrt(l2_f)
  else

     do v1=1,nv
        do v2=1,ng
           v=mesh%etat2(v1)+(1.0d0+gausslob%node(v2))/mesh%jac(nx+1,v1)
           do x1=1,nx
              do x2=1,ng
                 momentum=momentum+v*dist((x1-1)*ng+x2,(v1-1)*ng+v2)* &
                      & gausslob%weigh(x2)*gausslob%weigh(v2)/mesh%jac(x1,v1)
                 !kinetik energy
                 k_en=k_en+v**2*dist((x1-1)*ng+x2,(v1-1)*ng+v2)* &
                      & gausslob%weigh(x2)*gausslob%weigh(v2)/mesh%jac(x1,v1)

                 !norms 1 and 2 of distribution
                 l1_f=l1_f+abs(dist((x1-1)*ng+x2,(v1-1)*ng+v2))*gausslob%weigh(x2)* &
                      & gausslob%weigh(v2)/mesh%jac(x1,v1)
                 l2_f=l2_f+dist((x1-1)*ng+x2,(v1-1)*ng+v2)**2*gausslob%weigh(x2)* &
                      & gausslob%weigh(v2)/mesh%jac(x1,v1)
                 int_f=int_f+dist((x1-1)*ng+x2,(v1-1)*ng+v2)*gausslob%weigh(x2)* &
                      & gausslob%weigh(v2)/mesh%jac(x1,v1)
              end do
           end do
        end do
     end do
     l2_f=sqrt(l2_f)
  end if

  open(17,file=ffield)
  do x1=1,nx
     do x2=1,ng
        !electric/potential energy
        x=mesh%etat1(x1)+(1.0d0+gausslob%node(x2))/mesh%jac(x1,nv+1)
        em_en=em_en+dg_plan%field((x1-1)*ng+x2)**2*gausslob%weigh(x2)/mesh%jac(x1,nv+1)
        write(17,*)x,dg_plan%field((x1-1)*ng+x2),dg_plan%phi((x1-1)*ng+x2), &
             & dg_plan%rho((x1-1)*ng+x2)*mesh%jac(x1,nv+1)/gausslob%weigh(x2)+dg_plan%norma
     end do
     !stabilisation term (see Blanca's DG method for the 1D VP system)
     phi_jump=phi_jump+(dg_plan%phi(x1*ng)-dg_plan%phi(modulo(x1*ng+1-1,nx*ng)+1))**2
  end do
  energy=(k_en+em_en+phi_jump*c11)/2.0d0
  write(15,*)dg_plan%t,momentum,energy,k_en,sqrt(em_en),phi_jump,l1_f,l2_f, &
       & minval(dist),maxval(dist),int_f
  close(17)
!!$!!!>>>>>>>
!!$stop!!!!!!
!!$!!!<<<<<<<
  ! time loop begin
  do step=1,nb_step
     !!!>>>Blanca"s case
     !dg_plan%bound=-sqrt(sll_pi)/8.0d0*cos(2.0d0/sll_pi*dg_plan%t)
     !!!<<<
     call time_integration_4dg(dg_plan,gausslob,mesh,dist,distp1)
     dg_plan%t=real(step,8)*dt
     dist=distp1

     if (modulo(step,th)==0 .or. step==nb_step) then
        momentum=0.0d0
        energy=0.0d0
        k_en=0.0d0
        em_en=0.0d0
        phi_jump=0.0d0
        l1_f=0.0d0
        l2_f=0.0d0
        int_f=0.0d0
        
        dg_plan%rho=0.0d0
!!$        dg_plan%norma=0.0d0
        do x1=1,nx
           do x2=1,ng
              do v1=1,nv
                 do v2=1,ng
                    dg_plan%rho((x1-1)*ng+x2)=dg_plan%rho((x1-1)*ng+x2)+ &
                      & dist((x1-1)*ng+x2,(v1-1)*ng+v2)*gausslob%weigh(v2)/mesh%jac(nx+1,v1)
                 end do
              end do
!!$              dg_plan%norma=dg_plan%norma+dg_plan%rho((x1-1)*ng+x2)* &
!!$                   & gausslob%weigh(x2)/mesh%jac(x1,nv+1)
           end do
        end do
!!$        dg_plan%norma=dg_plan%norma/(x_max-x_min)
        do x1=1,nx
           do x2=1,ng
              dg_plan%rho((x1-1)*ng+x2)=(dg_plan%rho((x1-1)*ng+x2)-dg_plan%norma)* &
                   & gausslob%weigh(x2)/mesh%jac(x1,nv+1)
           end do
        end do
        dg_plan%rho(dg_plan%x0)=dg_plan%bound ! boudary condition = \Phi(0,t) = alpha

        dg_plan%phi=0.0d0
        dg_plan%field=0.0d0
        call solve(dg_plan%poisson_vp,dg_plan%rho,dg_plan%phi)
        dg_plan%field=matmul(dg_plan%poisson_vp%mat_field,dg_plan%phi)

        if (sym_mom_comp) then
           if (modulo(nv,2)==0) then !nv even
              do v1=1,nv/2
                 do v2=1,ng
                    v=mesh%etat2(v1)+(1.0d0+gausslob%node(v2))/mesh%jac(nx+1,v1)
                    vv=mesh%etat2(nv+1-v1)+(1.0d0+gausslob%node(ng+1-v2))/mesh%jac(nx+1,nv+1-v1)
                    do x1=1,nx
                       do x2=1,ng
                          momentum=momentum+v*dist((x1-1)*ng+x2,(v1-1)*ng+v2)* &
                               & gausslob%weigh(x2)*gausslob%weigh(v2)/mesh%jac(x1,v1)+ &
                               & vv*dist((x1-1)*ng+x2,(nv-v1)*ng+ng+1-v2)* &
                               & gausslob%weigh(x2)*gausslob%weigh(ng+1-v2)/mesh%jac(x1,nv+1-v1)
                          !kinetik energy
                          k_en=k_en+v**2*dist((x1-1)*ng+x2,(v1-1)*ng+v2)* &
                               & gausslob%weigh(x2)*gausslob%weigh(v2)/mesh%jac(x1,v1)+ &
                               & vv**2*dist((x1-1)*ng+x2,(nv-v1)*ng+ng+1-v2)* &
                               & gausslob%weigh(x2)*gausslob%weigh(ng+1-v2)/mesh%jac(x1,nv+1-v1)

                          !norms 1 and 2 of distribution
                          l1_f=l1_f+abs(dist((x1-1)*ng+x2,(v1-1)*ng+v2))*gausslob%weigh(x2)* &
                               & gausslob%weigh(v2)/mesh%jac(x1,v1)+ &
                               & abs(dist((x1-1)*ng+x2,(nv-v1)*ng+ng+1-v2))* &
                               & gausslob%weigh(x2)*gausslob%weigh(ng+1-v2)/mesh%jac(x1,nv+1-v1)
                          l2_f=l2_f+dist((x1-1)*ng+x2,(v1-1)*ng+v2)**2*gausslob%weigh(x2)* &
                               & gausslob%weigh(v2)/mesh%jac(x1,v1)+ &
                               & dist((x1-1)*ng+x2,(nv-v1)*ng+ng+1-v2)**2* &
                               & gausslob%weigh(x2)*gausslob%weigh(ng+1-v2)/mesh%jac(x1,nv+1-v1)
                          int_f=int_f+dist((x1-1)*ng+x2,(v1-1)*ng+v2)*gausslob%weigh(x2)* &
                               & gausslob%weigh(v2)/mesh%jac(x1,v1)+ &
                               & dist((x1-1)*ng+x2,(nv-v1)*ng+ng+1-v2)* &
                               & gausslob%weigh(x2)*gausslob%weigh(ng+1-v2)/mesh%jac(x1,nv+1-v1)
                       end do
                    end do
                 end do
              end do
           else !nv ode
              do v1=1,floor(nv/2.)-1
                 do v2=1,ng
                    v=mesh%etat2(v1)+(1.0d0+gausslob%node(v2))/mesh%jac(nx+1,v1)
                    vv=mesh%etat2(nv+1-v1)+(1.0d0+gausslob%node(ng+1-v2))/mesh%jac(nx+1,nv+1-v1)
                    do x1=1,nx
                       do x2=1,ng
                          momentum=momentum+v*dist((x1-1)*ng+x2,(v1-1)*ng+v2)* &
                               & gausslob%weigh(x2)*gausslob%weigh(v2)/mesh%jac(x1,v1)+ &
                               & vv*dist((x1-1)*ng+x2,(nv-v1)*ng+ng+1-v2)* &
                               & gausslob%weigh(x2)*gausslob%weigh(ng+1-v2)/mesh%jac(x1,nv+1-v1)
                          !kinetik energy
                          k_en=k_en+v**2*dist((x1-1)*ng+x2,(v1-1)*ng+v2)* &
                               & gausslob%weigh(x2)*gausslob%weigh(v2)/mesh%jac(x1,v1)+ &
                               & vv**2*dist((x1-1)*ng+x2,(nv-v1)*ng+ng+1-v2)* &
                               & gausslob%weigh(x2)*gausslob%weigh(ng+1-v2)/mesh%jac(x1,nv+1-v1)

                          !norms 1 and 2 of distribution
                          l1_f=l1_f+abs(dist((x1-1)*ng+x2,(v1-1)*ng+v2))*gausslob%weigh(x2)* &
                               & gausslob%weigh(v2)/mesh%jac(x1,v1)+ &
                               & abs(dist((x1-1)*ng+x2,(nv-v1)*ng+ng+1-v2))* &
                               & gausslob%weigh(x2)*gausslob%weigh(ng+1-v2)/mesh%jac(x1,nv+1-v1)
                          l2_f=l2_f+dist((x1-1)*ng+x2,(v1-1)*ng+v2)**2*gausslob%weigh(x2)* &
                               & gausslob%weigh(v2)/mesh%jac(x1,v1)+ &
                               & dist((x1-1)*ng+x2,(nv-v1)*ng+ng+1-v2)**2* &
                               & gausslob%weigh(x2)*gausslob%weigh(ng+1-v2)/mesh%jac(x1,nv+1-v1)
                          int_f=int_f+dist((x1-1)*ng+x2,(v1-1)*ng+v2)*gausslob%weigh(x2)* &
                               & gausslob%weigh(v2)/mesh%jac(x1,v1)+ &
                               & dist((x1-1)*ng+x2,(nv-v1)*ng+ng+1-v2)* &
                               & gausslob%weigh(x2)*gausslob%weigh(ng+1-v2)/mesh%jac(x1,nv+1-v1)
                       end do
                    end do
                 end do
              end do
              v1=floor(nv/2.)
              do v2=1,ng
                 v=mesh%etat2(v1)+(1.0d0+gausslob%node(v2))/mesh%jac(nx+1,v1)
                 do x1=1,nx
                    do x2=1,ng
                       momentum=momentum+v*dist((x1-1)*ng+x2,(v1-1)*ng+v2)* &
                            & gausslob%weigh(x2)*gausslob%weigh(v2)/mesh%jac(x1,v1)
                       !kinetik energy
                       k_en=k_en+v**2*dist((x1-1)*ng+x2,(v1-1)*ng+v2)* &
                            & gausslob%weigh(x2)*gausslob%weigh(v2)/mesh%jac(x1,v1)

                       !norms 1 and 2 of distribution
                       l1_f=l1_f+abs(dist((x1-1)*ng+x2,(v1-1)*ng+v2))*gausslob%weigh(x2)* &
                            & gausslob%weigh(v2)/mesh%jac(x1,v1)
                       l2_f=l2_f+dist((x1-1)*ng+x2,(v1-1)*ng+v2)**2*gausslob%weigh(x2)* &
                            & gausslob%weigh(v2)/mesh%jac(x1,v1)
                       int_f=int_f+dist((x1-1)*ng+x2,(v1-1)*ng+v2)*gausslob%weigh(x2)* &
                            & gausslob%weigh(v2)/mesh%jac(x1,v1)
                    end do
                 end do
              end do
           end if
           l2_f=sqrt(l2_f)
        else

           do v1=1,nv
              do v2=1,ng
                 v=mesh%etat2(v1)+(1.0d0+gausslob%node(v2))/mesh%jac(nx+1,v1)
                 do x1=1,nx
                    do x2=1,ng
                       momentum=momentum+v*dist((x1-1)*ng+x2,(v1-1)*ng+v2)* &
                            & gausslob%weigh(x2)*gausslob%weigh(v2)/mesh%jac(x1,v1)
                       !kinetik energy
                       k_en=k_en+v**2*dist((x1-1)*ng+x2,(v1-1)*ng+v2)* &
                            & gausslob%weigh(x2)*gausslob%weigh(v2)/mesh%jac(x1,v1)

                       !norms 1 and 2 of distribution
                       l1_f=l1_f+abs(dist((x1-1)*ng+x2,(v1-1)*ng+v2))*gausslob%weigh(x2)* &
                            & gausslob%weigh(v2)/mesh%jac(x1,v1)
                       l2_f=l2_f+dist((x1-1)*ng+x2,(v1-1)*ng+v2)**2*gausslob%weigh(x2)* &
                            & gausslob%weigh(v2)/mesh%jac(x1,v1)
                       int_f=int_f+dist((x1-1)*ng+x2,(v1-1)*ng+v2)*gausslob%weigh(x2)* &
                            & gausslob%weigh(v2)/mesh%jac(x1,v1)
                    end do
                 end do
              end do
           end do
           l2_f=sqrt(l2_f)

        end if

        do x1=1,nx
           do x2=1,ng
              x=mesh%etat1(x1)+(1.0d0+gausslob%node(x2))/mesh%jac(x1,nv+1)
              em_en=em_en+dg_plan%field((x1-1)*ng+x2)**2*gausslob%weigh(x2)/mesh%jac(x1,nv+1)
           end do
           !stabilisation term (see Blanca's DG method for the 1D VP system)
           phi_jump=phi_jump+(dg_plan%phi(x1*ng)-dg_plan%phi(modulo(x1*ng+1-1,nx*ng)+1))**2
        end do
        energy=(k_en+em_en+phi_jump*c11)/2.0d0
        write(15,*)dg_plan%t,momentum,energy,k_en,sqrt(em_en),phi_jump,l1_f,l2_f, &
             & minval(dist),maxval(dist),int_f
     end if

     if (modulo(step,th_large)==0 .or. step==nb_step) then
        write(fdist,*)'dist_'
        write(ffield,*)'field_'
        fdist=trim(adjustl(fdist))
        ffield=trim(adjustl(ffield))
        lfor=step
        lfor(1)=step/10**(len-1)
        write(fdist,*)trim(adjustl(fdist))//char(lfor(1)+48)
        write(ffield,*)trim(adjustl(ffield))//char(lfor(1)+48)
        do i=2,len
           do j=1,i-1
              lfor(i)=lfor(i)-lfor(j)*10**(len-j)
           end do
           lfor(i)=lfor(i)/10**(len-i)
           write(fdist,*)trim(adjustl(fdist))//char(lfor(i)+48)
           write(ffield,*)trim(adjustl(ffield))//char(lfor(i)+48)
        end do

        fdist=trim(adjustl(fdist))
        ffield=trim(adjustl(ffield))

        open(16,file=fdist)
        do v1=1,nv
           do v2=1,ng
              v=mesh%etat2(v1)+(1.0d0+gausslob%node(v2))/mesh%jac(nx+1,v1)
              do x1=1,nx
                 do x2=1,ng
                    x=mesh%etat1(x1)+(1.0d0+gausslob%node(x2))/mesh%jac(x1,nv+1)
                    write(16,*)x,v,dist((x1-1)*ng+x2,(v1-1)*ng+v2)
                 end do
              end do
              write(16,*)''
           end do
        end do
        close(16)

        open(17,file=ffield)
        do x1=1,nx
           do x2=1,ng
              x=mesh%etat1(x1)+(1.0d0+gausslob%node(x2))/mesh%jac(x1,nv+1)
              write(17,*)x,dg_plan%field((x1-1)*ng+x2),dg_plan%phi((x1-1)*ng+x2), &
                   & dg_plan%rho((x1-1)*ng+x2)*mesh%jac(x1,nv+1)/gausslob%weigh(x2)+dg_plan%norma
           end do
        end do
        close(17)
     end if

     if (modulo(step,th_out)==0 .or. step==nb_step) then
        write(*,form)step,'/',nb_step,'done, time=',real(dg_plan%t,4)
     end if
  end do
  close(15)

  deallocate(dist,distp1,lfor)
  call delete(gausslob)
  call delete(dg_plan)
  call delete(mesh)

contains

  subroutine param_out(x_min,x_max,v_min,v_max,nx,nv,ng,unif,c11,c12,tf,dt,nt,th,thl)
    
    sll_real64,intent(in) :: x_min,x_max,v_min,v_max,c11,c12,tf,dt
    sll_int32,intent(in)  :: nx,nv,ng,nt,th,thl
    logical,intent(in)    :: unif

    open(20,file="parameters",action="write")!,status="new")
    write(20,*)"x bounadries :",x_min,x_max
    write(20,*)"v bounadries :",v_min,v_max
    write(20,*)""
    if (unif) then
       write(20,*)"uniforme mesh"
    else
       write(20,*)"non-uniforme mesh"
    end if
    if (bc_v==sll_periodic) then
       write(20,*)"periodicity in direction v"
    else
       write(20,*)"Dirichlet condition in v"
    end if
    write(20,*)""
    write(20,*)"number of step in direction x  :",nx
    write(20,*)"number of step in direction v  :",nv
    write(20,*)"number of Gauss-Lobatto points :",ng
    write(20,*)""
    write(20,*)"flux coefficient c11 :",c11
    write(20,*)"flux coefficient c12 :",c12
    write(20,*)""
    write(20,*)"final time :",real(tf,4)
    write(20,*)"number of time steps :",nt
    write(20,*)"size of time steps   :",real(dt,4)
    write(20,*)"frequency of time historic :",real(th*dt,4),th
    write(20,*)"number of time historic    :",nt/th+1
    write(20,*)"frequency of snapshot :",real(thl*dt,4),thl
    write(20,*)"number of snapshot    :",nt/thl+1
    write(20,*)""
    if (x_upwind) then
       write(20,*)"upwind flux in direction x"
    else
       write(20,*)"centered flux in direction x"
    end if
    if (v_upwind) then
       write(20,*)"upwind flux in direction v"
    else
       write(20,*)"centered flux in direction v"
    end if
    close(20)

  end subroutine param_out

end program VP_DG

!only to check the convergence of the Poisson solver
!and of the full solver
!those function greatly increase the execution time, so use them only if necessary
!see comments under the functions

!!$  function interp_poly_1d(x,gll,f) result(res)
!!$
!!$    sll_real64,intent(in) :: x
!!$    type(gausslobatto1d),intent(in) :: gll
!!$    sll_real,dimension(:),intent(in) :: f
!!$    sll_real64 :: res
!!$
!!$    sll_real64 :: temp
!!$    sll_int32 :: ng,i,j
!!$
!!$    ng=gll%degree+1
!!$    res=0.0d0
!!$
!!$    do i=1,ng
!!$       temp=1.0d0
!!$       do j=1,ng
!!$          if (j/=i) then
!!$             temp=temp*(x-gll%node(j))/(gll%node(i)-gll%node(j))
!!$          end if
!!$       end do
!!$       res=res+temp*f(i)
!!$    end do
!!$
!!$  end function interp_poly_1d
!!$
!!$  function interp_poly_2d(x,y,gll,f,mesh) result(res)
!!$
!!$    sll_real64,intent(in) :: x,y
!!$    type(gausslobatto1d),intent(in) :: gll
!!$    type(non_unif_cart_mesh),intent(in) :: mesh
!!$    sll_real,dimension(:,:),intent(in) :: f
!!$    sll_real64 :: res
!!$
!!$    sll_real64 :: temp
!!$    sll_int32 :: ng,i,j,k,x1,v1,nx,nv
!!$    
!!$    nx=mesh%n_etat1
!!$    nv=mesh%n_etat2
!!$
!!$    ng=gll%degree+1
!!$    res=0.0d0
!!$
!!$    x1=0
!!$    v1=0
!!$
!!$    do while(x<mesh%etat1(x1))
!!$       x1=x1+1
!!$    end do
!!$    do while(y<mesh%etat2(v1))
!!$       v1=v1+1
!!$    end do
!!$    x1=x1-1
!!$    v1=v1-1
!!$
!!$    do j=1,ng
!!$       do i=1,ng
!!$          temp=1.0d0
!!$          do k=1,ng
!!$             if (k/=j) then
!!$                temp=temp*(y-(mesh%etat2(v1)+(1.0d0+gll%node(k))/mesh%jac(nx+1,v1)))/ &
!!$                     & ((gll%node(j)-gll%node(k))/mesh%jac(nx+1,v1))
!!$             end if
!!$             if (k/=i)then
!!$                temp=temp*(x-(mesh%etat2(x1)+(1.0d0+gll%node(k))/mesh%jac(x1,nv+1)))/ &
!!$                     & ((gll%node(i)-gll%node(k))/mesh%jac(x1,nv+1))
!!$             end if
!!$          end do
!!$          res=res+temp*f((x1-1)*ng+i,(v1-1)*ng+j)
!!$       end do
!!$    end do
!!$    
!!$  end function interp_poly_2d

!!$end program VP_DG

!!$ This part was originally after the time loop.
!!$ Its purpose was for self-convergence test but it was not used.
!!$ Since it was useles in the code i wanted to remove it, but it 
!!$ still might be usefull so I did not want to loose it. Finally I
!!$ placed it there as comment.
!!$ Madaule Eric

!!$  !This part is for self-convergence tests
!!$  if (var==0) then
!!$     open(13,file='ref_data')
!!$     write(13,*)ng
!!$     write(13,*)gausslob%weigh
!!$     write(13,*)nx,nv
!!$     write(13,*)mesh%jac
!!$     close(13)
!!$     open(11,file='ref_dist')
!!$     do v1=1,nv
!!$        do v2=1,ng
!!$           do x1=1,nx
!!$              do x2=1,ng
!!$                 x=mesh%etat1(x1)+(1.0d0+gausslob%node(x2))/mesh%jac(x1,nv+1)
!!$                 v=mesh%etat2(v1)+(1.0d0+gausslob%node(v2))/mesh%jac(nx+1,v1)
!!$
!!$                 write(11,*)x,v,dist((x1-1)*ng+x2,(v1-1)*ng+v2), &
!!$                      & (2.0d0-cos(2.0d0*x-2.0d0*sll_pi*dg_plan%t))* &
!!$                      & exp(-0.25d0*(4.0d0*v-1.0d0)**2)
!!$              end do
!!$           end do
!!$           write(11,*)''
!!$        end do
!!$     end do
!!$     close(11)
!!$  else if (var==1) then
!!$     linf=0.0d0
!!$     l1=0.0d0
!!$     l2=0.0d0
!!$
!!$     open(13,file='ref_data',action='read',status='old')
!!$     read(13,*)nng
!!$     allocate(weigh(nng))
!!$     read(13,*)nnx,nnv
!!$     allocate(jac(nnx+1,nnv+1))
!!$     read(13,*)jac
!!$     close(13)
!!$
!!$     open(11,file='ref_dist',action='read',status='old')
!!$     do v1=1,nnv
!!$        do v2=1,nng
!!$           do x1=1,nnx
!!$              do x2=1,nng
!!$                 read(11,*)xx,vv,ref
!!$                 current=interp_poly_2d(xx,vv,gausslob,dist,mesh)
!!$
!!$                 linf=max(linf,abs(ref-current))
!!$                 l1=l1+abs(ref-current)/jac(x1,v1)*weigh(x2)*weigh(v2)
!!$                 l2=l2+abs(ref-current)**2/jac(x1,v1)*weigh(x2)*weigh(v2)
!!$              end do
!!$           end do
!!$        end do
!!$     end do
!!$     close(11)
!!$     l2=sqrt(l2)
!!$     open(12,file='norms',position='append')
!!$     write(12,*)'# nx,nv,ng,linf,l1,l2'
!!$     write(12,*)nx,nv,ng,linf,l1,l2
!!$     close(12)
!!$
!!$     deallocate(weigh,jac)
!!$  end if
!!$  !end of self-convergence tests
