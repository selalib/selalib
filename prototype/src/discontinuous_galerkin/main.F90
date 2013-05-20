program VP_DG
#include "sll_working_precision.h"

!from Madaule Éric :
!Warning : this code and all the modules it calls and I wrote are not optimized in
!term of computation time

  use sll_nu_cart_mesh
  use gausslobatto
  use poisson4dg
  use mod_sparse
  use sll_constants
  use timestep4dg

  !use mod_octave_io_sparse
  !use mod_umfpack

  implicit none

  !plan for time step with DG
  type(dg_time_steping) :: dg_plan
  !distribution function and right hand side
  sll_real64,dimension(:,:),allocatable :: dist,distp1
  !number of time steb, number of cells un direction x and v,time step
  sll_int32 :: nb_step,nx,nv,step
  !boudary in direction x and v
  sll_real64 :: x_min,x_max,v_min,v_max
  !mesh
  type(non_unif_cart_mesh) :: mesh
  !coefficients for fluxes
  sll_real64 :: c11,c12,c22
  !number of Gauss-Lobatto points
  sll_int32 :: ng
  !time step and finalt time
  sll_real64 :: dt,tf
  !Gauss-Lobatto
  type(gausslobatto1D) :: gausslob
  !elements of Phi equal to 0 so Phi(x=0)=0
  sll_int32 :: xbound
  !space variables
  sll_real64 :: x,v

  !indexes for loops
  sll_int32 :: x1,x2,v1,v2

  !error on Poisson
  sll_real64 :: linf,l1,l2

  !to check self-convergence
  sll_int32 :: var,nnx,nnv,nng
  sll_real64 :: ref,xx,vv,current
  sll_real64,dimension(:),allocatable :: weigh
  sll_real64,dimension(:,:),allocatable :: jac

  ! for the python script polar-exe.py
  !namelist /param/ nx,nv,ng,dt,var
  !read(*,NML=param)

  nx=20
  nv=20
  ng=5
  dt=0.1d0
  var=0

  !definition of geometry and data
!!$  x_min=0.0d0
!!$  x_max=sll_pi
!!$  v_min=0
!!$  v_max=sll_pi

!!$  x_min=0.0d0
!!$  x_max=4.0d0*sll_pi
!!$  v_min=-2.0d0
!!$  v_max=2.0d0

  !!!>>>Blanca's case
  x_min=-sll_pi
  x_max=sll_pi
  v_min=-sll_pi
  v_max=sll_pi
  !!!<<<

!!$  x_min=0.0d0
!!$  x_max=1.0d0
!!$  v_min=-1.0d0
!!$  v_max=1.0d0

!!$  x_min=-1.0d0
!!$  x_max=1.0d0
!!$  v_min=-1.0d0
!!$  v_max=1.0d0

!!$  nx=52
!!$  nv=52
!!$  ng=5

  !xbound=ng*nx/2
  xbound=1

  allocate(dist(nx*ng,nv*ng),distp1(nx*ng,nv*ng))

  !definition or time step, delta_t and final time
!!$  dt=0.0001d0
  tf=1.0d0
  nb_step=ceiling(tf/dt)
  tf=dt*nb_step
  print*,'size of time step   :',dt
  print*,'number of time step :',nb_step
  print*,'final time          :',tf

  !flux coefficients
  c22=0.0d0
  c12=0.5d0
  c11=0.1d0

  call init_gausslobatto_1d(ng,gausslob)
  !call init_gausslobatto_1d(3*ng,gll)
  call init_nu_cart_mesh(nx,nv,mesh)
  mesh%d_etat1=(x_max-x_min)/real(nx,8)
  mesh%d_etat2=(v_max-v_min)/real(nv,8)
  call fill_node_nuc_mesh(x_min,v_min,mesh)
  call init_timesteping_4dg(dg_plan,SLL_RK4,gausslob,mesh%jac(1:nx,nv+1),dt,xbound, &
       & nx,nv,ng,c11,c12)

  !dist() !construction of distribution function
  !x_i is indexed on both mesh nodes and GLL nodes, so to have the postion in x
  !one shoul take the lower node of the element and add the part due to GLL
  !the 1.0d0 is to compensate the fact that GLL is done on [-1;1]
  !same is done for v

  !test distribution : f(x,v)=sin(x)sin(v), (x,v)\in [0;pi]² (to check rhs, independants on t)
  !exact rhs = -v*cos(x)sin(x)+sin(2x)cos(v)
  dist=0.0d0
  open(14,file='dist_init')
  do v1=1,nv
     do v2=1,ng
        do x1=1,nx
           do x2=1,ng
              x=mesh%etat1(x1)+(1.0d0+gausslob%node(x2))/mesh%jac(x1,nv+1)
              v=mesh%etat2(v1)+(1.0d0+gausslob%node(v2))/mesh%jac(nx+1,v1)
!!$              dist((x1-1)*ng+x2,(v1-1)*ng+v2)=sin(x)*sin(v)

!!$              dist((x1-1)*ng+x2,(v1-1)*ng+v2)=(exp(-200.0d0*(v-0.8d0)**2)+ &
!!$                   & exp(-200.0d0*(v+0.8d0)**2))!*(cos(3.0d0*x)+cos(6.0d0*x)+cos(18.0d0*x))

              !!!>>>Blanca's test case
              dist((x1-1)*ng+x2,(v1-1)*ng+v2)=(2.0d0-cos(2.0d0*x))*exp(-0.25d0*(4.0d0*v-1.0d0)**2)
              !!!<<<

!!$              if (abs(v)<=0.5d0) then
!!$                 dist((x1-1)*ng+x2,(v1-1)*ng+v2)=1.0d0* &
!!$                      & sin(x*sll_pi)
!!$              end if

!!$              dist((x1-1)*ng+x2,(v1-1)*ng+v2)=(x**2-1.0d0)*(v**2-1.0d0)
              write(14,*)x,v,dist((x1-1)*ng+x2,(v1-1)*ng+v2)
           end do
        end do
        write(14,*)''
     end do
  end do
  close(14)

!!$  dg_plan%k(2,:,:)=dist
!!$  call rk4dg_step(dg_plan,gausslob,mesh)
!!$  
!!$  open(16,file="mneme")
!!$  write(16,*)' '
!!$  write(16,*)"#rhs dist"
!!$  do v1=1,nv
!!$     do v2=1,ng
!!$        do x1=1,nx
!!$           do x2=1,ng
!!$              x=mesh%etat1(x1)+(1.0d0+gausslob%node(x2))/mesh%jac(x1,nv+1)
!!$              v=mesh%etat2(v1)+(1.0d0+gausslob%node(v2))/mesh%jac(nx+1,v1)
!!$              write(16,*)x,v,dg_plan%rhs((x1-1)*ng+x2,(v1-1)*ng+v2), &
!!$                   & dist((x1-1)*ng+x2,(v1-1)*ng+v2)
!!$           end do
!!$        end do
!!$        write(16,*)' '
!!$     end do
!!$  end do
!!$  close(16)
  
  ! time loop
  do step=1,nb_step
     dg_plan%t=real(step-1,8)*dt
     call rk4_4dg_1d(dg_plan,gausslob,mesh,dist,distp1)
     print*,maxval(distp1)
     dist=distp1
     !print*,"step",step,'done'
  end do

  if (var==0) then
     open(13,file='ref_data')
     write(13,*)ng
     write(13,*)gausslob%weigh
     write(13,*)nx,nv
     write(13,*)mesh%jac
     close(13)
     open(11,file='ref_dist')
     do v1=1,nv
        do v2=1,ng
           do x1=1,nx
              do x2=1,ng
                 x=mesh%etat1(x1)+(1.0d0+gausslob%node(x2))/mesh%jac(x1,nv+1)
                 v=mesh%etat2(v1)+(1.0d0+gausslob%node(v2))/mesh%jac(nx+1,v1)

                 write(11,*)x,v,dist((x1-1)*ng+x2,(v1-1)*ng+v2)
              end do
           end do
           write(11,*)''
        end do
     end do
     close(11)
  else if (var==1) then
     linf=0.0d0
     l1=0.0d0
     l2=0.0d0

     open(13,file='ref_data',action='read',status='old')
     read(13,*)nng
     allocate(weigh(nng))
     read(13,*)nnx,nnv
     allocate(jac(nnx+1,nnv+1))
     read(13,*)jac
     close(13)

     open(11,file='ref_dist',action='read',status='old')
     do v1=1,nnv
        do v2=1,nng
           do x1=1,nnx
              do x2=1,nng
                 read(11,*)xx,vv,ref
                 current=interp_poly_2d(xx,vv,gausslob,dist,mesh)

                 linf=max(linf,abs(ref-current))
                 l1=l1+abs(ref-current)/jac(x1,v1)*weigh(x2)*weigh(v2)
                 l2=l2+abs(ref-current)**2/jac(x1,v1)*weigh(x2)*weigh(v2)
              end do
           end do
        end do
     end do
     close(11)
     l2=sqrt(l2)
     open(12,file='clio',position='append')
     write(12,*)'# nx,nv,ng,linf,l1,l2'
     write(12,*)nx,nv,ng,linf,l1,l2
     close(12)

     deallocate(weigh,jac)
  end if

  deallocate(dist,distp1)
  call delete_gausslobatto_1D(gausslob)
  !call delete_gausslobatto_1D(gll)
  call clear(dg_plan)
  call delete_nu_cart_mesh(mesh)

contains
!only to check the convergence of the Poisson solver
!and of the full solver

  function interp_poly_1d(x,gll,f) result(res)

    sll_real64,intent(in) :: x
    type(gausslobatto1d),intent(in) :: gll
    sll_real,dimension(:),intent(in) :: f
    sll_real64 :: res

    sll_real64 :: temp
    sll_int32 :: ng,i,j

    ng=gll%degree+1
    res=0.0d0

    do i=1,ng
       temp=1.0d0
       do j=1,ng
          if (j/=i) then
             temp=temp*(x-gll%node(j))/(gll%node(i)-gll%node(j))
          end if
       end do
       res=res+temp*f(i)
    end do

  end function interp_poly_1d

  function interp_poly_2d(x,y,gll,f,mesh) result(res)

    sll_real64,intent(in) :: x,y
    type(gausslobatto1d),intent(in) :: gll
    type(non_unif_cart_mesh),intent(in) :: mesh
    sll_real,dimension(:,:),intent(in) :: f
    sll_real64 :: res

    sll_real64 :: temp
    sll_int32 :: ng,i,j,k,x1,v1,nx,nv
    
    nx=mesh%n_etat1
    nv=mesh%n_etat2

    ng=gll%degree+1
    res=0.0d0

    x1=0
    v1=0

    do while(x<mesh%etat1(x1))
       x1=x1+1
    end do
    do while(y<mesh%etat2(v1))
       v1=v1+1
    end do
    x1=x1-1
    v1=v1-1

    do j=1,ng
       do i=1,ng
          temp=1.0d0
          do k=1,ng
             if (k/=j) then
                temp=temp*(y-(mesh%etat2(v1)+(1.0d0+gll%node(k))/mesh%jac(nx+1,v1)))/ &
                     & ((gll%node(j)-gll%node(k))/mesh%jac(nx+1,v1))
             end if
             if (k/=i)then
                temp=temp*(x-(mesh%etat2(x1)+(1.0d0+gll%node(k))/mesh%jac(x1,nv+1)))/ &
                     & ((gll%node(i)-gll%node(k))/mesh%jac(x1,nv+1))
             end if
          end do
          res=res+temp*f((x1-1)*ng+i,(v1-1)*ng+j)
       end do
    end do
    
  end function interp_poly_2d

!!$    do j=1,ng
!!$       do i=1,ng
!!$          temp=1.0d0
!!$          do k=1,ng
!!$             if (k/=j) then
!!$                temp=temp*(y-gll%node(k))/(gll%node(j)-gll%node(k))
!!$             end if
!!$             if (k/=i)then
!!$                temp=temp*(x-gll%node(k))/(gll%node(i)-gll%node(k))
!!$             end if
!!$          end do
!!$          res=res+temp*f(i,j)
!!$       end do
!!$    end do
!!$
!!$  end function interp_poly_2d

end program VP_DG
