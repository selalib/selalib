program VP_DG
#include "sll_working_precision.h"

!from Madaule Éric :
!Warning : this code and all the modules it calls and I wrote are not optimized in
!term of computation time. \n
!We kept the readability of the code over computation performances

  use sll_nu_cart_mesh
  use gausslobatto
  use poisson4dg
  use mod_sparse
  use sll_constants
  use timestep4dg

  !use mod_octave_io
  use mod_octave_io_sparse

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
  sll_real64 :: x,v,k
  !energy and momentum to check conservation
  sll_real64 :: momentum,l1_f,l2_f
  sll_real64 :: k_en,em_en,energy,phi_jump

  !indexes for loops
  sll_int32 :: x1,x2,v1,v2
  !display parameter
  sll_int32 :: len,i1,i2,i3,i4,i5,i6
  character(len=25) :: form,fdist,ffield
  sll_int32 :: irec ! for binary writting

  ! for the python script *.py
  !namelist /param/ nx,nv,ng
  !read(*,NML=param)

  !definition of geometry and data
  k=0.5d0

  x_min=0.0d0
  x_max=sll_pi
  v_min=0
  v_max=sll_pi

!!$  x_min=0.0d0
!!$  x_max=4.0d0*sll_pi
!!$  v_min=-2.0d0
!!$  v_max=2.0d0

!!$  !!!>>>Blanca's case
!!$  x_min=-sll_pi
!!$  x_max=sll_pi
!!$  v_min=-4
!!$  v_max=4
!!$  !!!<<<

!!$  x_min=0.0d0
!!$  x_max=1.0d0
!!$  v_min=-1.0d0
!!$  v_max=1.0d0

!!$  x_min=-1.0d0
!!$  x_max=1.0d0
!!$  v_min=-1.0d0
!!$  v_max=1.0d0

!!$  x_min=0.0d0
!!$  x_max=2.0d0*sll_pi/k
!!$  v_min=-10.0d0
!!$  v_max=10.0d0

  nx=30
  nv=50
  ng=5

  print*,'discretization caracteristics :'
  print"(3(a5,i3))",'nx=',nx,', nv=',nv,', ng=',ng

  !xbound=ng*nx/2
  xbound=1

  allocate(dist(nx*ng,nv*ng),distp1(nx*ng,nv*ng))

  !definition or time step, delta_t and final time
  dt=0.001d0
  tf=15.0d0
  nb_step=1!ceiling(tf/dt)
  th=20
  th_out=200
  th_large=1000
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
  form="(1x,i"//char(len+48)//",a1,i"//char(len+48)//",1x,a11,f6.2)"

  call init_gausslobatto_1d(ng,gausslob)
  call init_nu_cart_mesh(nx,nv,mesh)
  mesh%d_etat1=(x_max-x_min)/real(nx,8)
  mesh%d_etat2=(v_max-v_min)/real(nv,8)
  call fill_node_nuc_mesh(x_min,v_min,mesh)

  !flux coefficients
  c22=0.0d0
  c12=0.5d0
  c11=real(ng**2,8)/maxval(mesh%d_etat1)

  call init_timesteping_4dg(dg_plan,SLL_ee,gausslob,mesh,dt,xbound,c11,c12, &
       & alpha=0.0d0)
!       & norma=erf(5.0d0*sqrt(2.0d0)))

  !construction of distribution function
  !x_i is indexed on both mesh nodes and GLL nodes, so to have the postion in x
  !one shoul take the lower node of the element and add the part due to GLL
  !the 1.0d0 is to compensate the fact that GLL is done on [-1;1]
  !same is done for v

  !test distribution for the Poisson's problem : 
  !f(x,v)=sin(x)sin(v), (x,v)\in [0;pi]² (to check rhs, independants on t)
  !exact rhs = -v*cos(x)sin(x)+sin(2x)cos(v)
  dist=0.0d0
  inquire(iolength=irec)dist
  open(14,file='dist_000000')!,form="unformatted",access="direct",recl=irec)
  do v1=1,nv
     do v2=1,ng
        v=mesh%etat2(v1)+(1.0d0+gausslob%node(v2))/mesh%jac(nx+1,v1)
        do x1=1,nx
           do x2=1,ng
              x=mesh%etat1(x1)+(1.0d0+gausslob%node(x2))/mesh%jac(x1,nv+1)
              dist((x1-1)*ng+x2,(v1-1)*ng+v2)=sin(x)*sin(v)+1.0d0/sll_pi

!!$              dist((x1-1)*ng+x2,(v1-1)*ng+v2)=(exp(-200.0d0*(v-0.8d0)**2)+ &
!!$                   & exp(-200.0d0*(v+0.8d0)**2))!*(cos(3.0d0*x)+cos(6.0d0*x)+cos(18.0d0*x))

!!$              !!!>>>Blanca's test case
!!$              dist((x1-1)*ng+x2,(v1-1)*ng+v2)=(2.0d0-cos(2.0d0*x))* &
!!$                   & exp(-0.25d0*(4.0d0*v-1.0d0)**2)
!!$              !!!<<<

!!$              if (abs(v)<=0.5d0) then
!!$                 dist((x1-1)*ng+x2,(v1-1)*ng+v2)=1.0d0* &
!!$                      & sin(x*sll_pi)
!!$              end if

!!$              dist((x1-1)*ng+x2,(v1-1)*ng+v2)=(x**2-1.0d0)*(v**2-1.0d0)

!!$              !dist((x1-1)*ng+x2,(v1-1)*ng+v2)=(1.0d0+0.5d0*sin(k*x))* &
!!$              dist((x1-1)*ng+x2,(v1-1)*ng+v2)=(1.0d0-0.5d0*cos(k*x))* &
!!$                   !& (exp(-(v-1.0d0)**2*40.0d0)+exp(-(v+1.0d0)**2*40.0d0))
!!$                   & exp(-v**2/2.0d0)/sqrt(2.0d0*sll_pi)

!!$              dist((x1-1)*ng+x2,(v1-1)*ng+v2)=exp(-v**2)/sqrt(2.0d0*sll_pi)

              write(14,*)x,v,dist((x1-1)*ng+x2,(v1-1)*ng+v2)
              !write(14)x,v,dist((x1-1)*ng+x2,(v1-1)*ng+v2)
           end do
        end do
        write(14,*)''
     end do
  end do
  !write(14,rec=1)dist
  close(14)!,status='keep')

  dg_plan%t=0.0d0
  open(15,file='energy_momentum',action='write')
  write(15,*)"# t ; momentum ; total energy ; kinetic energy ; electromagnetic energy ;", &
       & " jump of phi ; ||f||_L1 ; ||f||_L2"
  momentum=0.0d0
  energy=0.0d0
  k_en=0.0d0
  em_en=0.0d0
  phi_jump=0.0d0
  l1_f=0.0d0
  l2_f=0.0d0

  dg_plan%rho=0.0d0
  dg_plan%norma=0.0d0
!!$  do x1=1,nx
!!$     do x2=1,ng
!!$        do v1=1,nv
!!$           do v2=1,ng
!!$              dg_plan%rho((x1-1)*ng+x2)=dg_plan%rho((x1-1)*ng+x2)+ &
!!$                   & dist((x1-1)*ng+x2,(v1-1)*ng+v2)*gausslob%weigh(v2)/mesh%jac(nx+1,v1)
!!$           end do
!!$        end do
!!$        dg_plan%rho((x1-1)*ng+x2)=(dg_plan%rho((x1-1)*ng+x2)-dg_plan%norma)* &
!!$             & gausslob%weigh(x2)/mesh%jac(x1,nv+1)
!!$     end do
!!$  end do
  do x1=1,nx*ng
     do v1=1,nv
        do v2=1,ng
           dg_plan%rho(x1)=dg_plan%rho(x1)+ &
                & dist(x1,(v1-1)*ng+v2)*gausslob%weigh(v2)/mesh%jac(nx+1,v1)
        end do
     end do
  end do
  do x1=1,nx
     do x2=1,ng
        dg_plan%norma=dg_plan%norma+dg_plan%rho((x1-1)*ng+x2)*gausslob%weigh(x2)/mesh%jac(x1,nv+1)
     end do
  end do
  dg_plan%norma=dg_plan%norma/(x_max-x_min)
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

  do v1=1,nv
     do v2=1,ng
        v=mesh%etat2(v1)+(1.0d0+gausslob%node(v2))/mesh%jac(nx+1,v1)
        do x1=1,nx
           do x2=1,ng
!!$              x=mesh%etat1(x1)+(1.0d0+gausslob%node(x2))/mesh%jac(x1,nv+1)

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
           end do
        end do
     end do
  end do
  l2_f=sqrt(l2_f)

  open(17,file='field_000000')
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
  energy=k_en+em_en+phi_jump*c11
  write(15,*)dg_plan%t,momentum,energy,k_en,em_en,phi_jump,l1_f,l2_f
  close(17)

  dg_plan%bound=0
  ! time loop begin
  do step=1,nb_step
     !dg_plan%bound=-sqrt(sll_pi)/8.0d0*cos(2.0d0/sll_pi*dg_plan%t)
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
        
        dg_plan%rho=0.0d0
        do x1=1,nx
           do x2=1,ng
              do v1=1,nv
                 do v2=1,ng
                    dg_plan%rho((x1-1)*ng+x2)=dg_plan%rho((x1-1)*ng+x2)+ &
                         & dist((x1-1)*ng+x2,(v1-1)*ng+v2)*gausslob%weigh(v2)/mesh%jac(nx+1,v1)
                 end do
              end do
              dg_plan%rho((x1-1)*ng+x2)=(dg_plan%rho((x1-1)*ng+x2)-dg_plan%norma)* &
                   & gausslob%weigh(x2)/mesh%jac(x1,nv+1)
           end do
        end do
        dg_plan%rho(dg_plan%x0)=dg_plan%bound ! boudary condition = \Phi(0,t) = alpha

        dg_plan%phi=0.0d0
        dg_plan%field=0.0d0
        call solve(dg_plan%poisson_vp,dg_plan%rho,dg_plan%phi)
        dg_plan%field=matmul(dg_plan%poisson_vp%mat_field,dg_plan%phi)

        do v1=1,nv
           do v2=1,ng
              v=mesh%etat2(v1)+(1.0d0+gausslob%node(v2))/mesh%jac(nx+1,v1)
              do x1=1,nx
                 do x2=1,ng
                    x=mesh%etat1(x1)+(1.0d0+gausslob%node(x2))/mesh%jac(x1,nv+1)

                    momentum=momentum+v*dist((x1-1)*ng+x2,(v1-1)*ng+v2)* &
                         & gausslob%weigh(x2)*gausslob%weigh(v2)/mesh%jac(x1,v1)
                    !kinetik energy
                    k_en=k_en+v**2*abs(dist((x1-1)*ng+x2,(v1-1)*ng+v2))* &
                         & gausslob%weigh(x2)*gausslob%weigh(v2)/mesh%jac(x1,v1)

                    !norms 1 and 2 of distribution
                    l1_f=l1_f+abs(dist((x1-1)*ng+x2,(v1-1)*ng+v2))*gausslob%weigh(x2)* &
                         & gausslob%weigh(v2)/mesh%jac(x1,v1)
                    l2_f=l2_f+dist((x1-1)*ng+x2,(v1-1)*ng+v2)**2*gausslob%weigh(x2)* &
                         & gausslob%weigh(v2)/mesh%jac(x1,v1)
                 end do
              end do
           end do
        end do

        do x1=1,nx
           do x2=1,ng
              x=mesh%etat1(x1)+(1.0d0+gausslob%node(x2))/mesh%jac(x1,nv+1)
              em_en=em_en+dg_plan%field((x1-1)*ng+x2)**2*gausslob%weigh(x2)/mesh%jac(x1,nv+1)
           end do
           !stabilisation term (see Blanca's DG method for the 1D VP system)
           phi_jump=phi_jump+(dg_plan%phi(x1*ng)-dg_plan%phi(modulo(x1*ng+1-1,nx*ng)+1))**2
        end do
        energy=k_en+em_en+phi_jump*c11
        write(15,*)dg_plan%t,momentum,energy,k_en,em_en,phi_jump,l1_f,l2_f
     end if

     if (modulo(step,th_large)==0 .or. step==nb_step) then
        i1=step/100000
        i2=(step-100000*i1)/10000
        i3=(step-100000*i1-10000*i2)/1000
        i4=(step-100000*i1-10000*i2-i3*1000)/100
        i5=(step-100000*i1-10000*i2-i3*1000-i4*100)/10
        i6= step-100000*i1-10000*i2-i3*1000-i4*100-i5*10
        write(ffield,*)'field_'//char(i1+48)//char(i2+48)//char(i3+48)//char(i4+48)// &
             & char(i5+48)//char(i6+48)
        write(fdist,*)'dist_'//char(i1+48)//char(i2+48)//char(i3+48)//char(i4+48)// &
             & char(i5+48)//char(i6+48)
        fdist=trim(adjustl(fdist))
        ffield=trim(adjustl(ffield))

        open(16,file=fdist)
        do v1=1,nv
           do v2=1,ng
              v=mesh%etat2(v1)+(1.0d0+gausslob%node(v2))/mesh%jac(nx+1,v1)
              do x1=1,nx
                 do x2=1,ng
                    x=mesh%etat1(x1)+(1.0d0+gausslob%node(x2))/mesh%jac(x1,nv+1)
                    write(16,*)x,v,dist((x1-1)*ng+x2,(v1-1)*ng+v2), &
                         & dg_plan%rhs((x1-1)*ng+x2,(v1-1)*ng+v2), &
                         & -v*cos(x)*sin(v)+sin(2.0d0*x)*cos(v)
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
              write(17,*)x,dg_plan%field((x1-1)*ng+x2)
           end do
        end do
        close(17)
     end if

     if (modulo(step,th_out)==0 .or. step==nb_step) then
        write(*,form)step,'/',nb_step,'done, time=',real(dg_plan%t,4)
     end if
  end do
  close(15)

  deallocate(dist,distp1)
  call delete(gausslob)
  call delete(dg_plan)
  call delete(mesh)

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
!!$ still might be usefull so I did not want to loose it, so I
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
