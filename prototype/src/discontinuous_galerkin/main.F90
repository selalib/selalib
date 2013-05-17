program VP_DG
#include "sll_working_precision.h"

  use sll_nu_cart_mesh
  use gausslobatto
  use poisson4dg
  use mod_sparse
  use sll_constants

  !use mod_octave_io_sparse
  use mod_umfpack

  implicit none

  !rho,E,phi
  sll_real64,dimension(:),allocatable :: rho,field_e,phi
  !distribution function
  sll_real64,dimension(:,:),allocatable :: dist,rhs
  !number of time steb, number of cells un direction x and v
  sll_int32 :: nx,nv,nb_step
  !boudary in direction x and v
  sll_real64 :: x_min,x_max,v_min,v_max
  !mesh
  type(non_unif_cart_mesh) :: mesh
  !coefficients for fluxes
  sll_real64 :: c11,c12,c22
  !number of Gauss-Lobatto points
  sll_int32 :: ng
  !CSC matrix array
  type(t_col) :: matvp,matvm
  type(t_col) :: fieldvp,fieldvm
  !time step and finalt time
  sll_real64 :: dt,tf
  !Gauss-Lobatto
  type(gausslobatto1D) :: gausslob

  !umfpack variables
  integer :: sys
  sll_int64 :: numeric
  !sll_int64, parameter :: umfpack_control=20,umfpack_info=90
  sll_real64 :: control(umfpack_control), info(umfpack_info)
  integer(umf_void) :: symbolic

  !indexes for loops
  sll_int32 :: i,j,k,x1,x2,v1,v2

  !error on Poisson
  !sll_real64 :: linf,l1,l2
  !variables for computation
  sll_real64 :: som1,som2

!!$  ! for the python script polar-exe.py
!!$  namelist /param/ nx,ng
!!$  read(*,NML=param)

  !definition of geometry and data
  x_min=0.0d0
  x_max=sll_pi
  v_min=0.0d0
  v_max=sll_pi
  nx=15
  nv=15
  ng=5
  call init_gausslobatto_1d(ng,gausslob)
  call init_nu_cart_mesh(nx,nv,mesh)
  mesh%d_etat1=(x_max-x_min)/real(nx,8)
  mesh%d_etat2=(v_max-v_min)/real(nv,8)
  call fill_node_nuc_mesh(x_min,v_min,mesh)

  allocate(rho(nx*ng),phi(nx*ng),field_e(nx*ng),dist(nx*ng,nv*ng),rhs(nx*ng,nv*ng))

  !definition or time step, delta_t and final time
  dt=0.01d0
  tf=1.0d0
  nb_step=ceiling(tf/dt)

  !flux coefficients
  c22=0.0d0
  c12=0.5d0
  c11=0.1d0
  !build the matrixes for the Poisson-problem
  !((D+F_E).M-¹.(D-F_\Phi)^T - C).\Phi = M.(\rho-1)
  !this is for v>0
  call poisson1d_matrix(gausslob,nx,mesh%jac(1:nx,nv+1),c11,c12,1,matvp,fieldvp)
  !this is for v<0
  c12=-0.5d0
  call poisson1d_matrix(gausslob,nx,mesh%jac(1:nx,nv+1),c11,c12,1,matvm,fieldvm)

!!$  open(17,file='melete')
!!$  call write_octave(fieldvp,'fvp',17)
!!$  close(17)

  !dist() !construction of distribution function
  !x_i is indexed on both mesh nodes and GLL nodes, so to have the postion in x
  !one shoul take the lower node of the element and add the part due to GLL
  !the 1.0d0 is to compensate the fact that GLL is done on [-1;1]
  !same is done for v and in rhs
  dist=0.0d0
  do x1=1,nx
     do v1=1,nv
        do x2=1,ng
           do v2=1,ng
              dist((i-1)*ng+j,(v1-1)*ng+v2)= &
!!$                   & (1.0d0-floor(abs(mesh%etat2(v1)+(1.0d0+gausslob%node(v2))/ &
!!$                   & mesh%jac(nx+1,v1))*4.0d0)/4.0d0)*exp((mesh%etat2(v1)+(1.0d0+ &
!!$                   & gausslob%node(v2))/mesh%jac(nx+1,v1))**2/0.1d0)* &
!!$                   & 
                   & sin(mesh%etat1(x1)+(1.0d0+gausslob%node(x2))/mesh%jac(x1,nv+1))* &
                   & sin(mesh%etat2(v1)+(1.0d0+gausslob%node(v2))/mesh%jac(nx+1,v1))
           end do
        end do
     end do
  end do

  rho=0.0d0
  do i=1,nx*ng
     do j=1,nv
        do k=1,ng
           rho(i)=rho(i)+dist(i,(j-1)*ng+k)*gausslob%weigh(k)/mesh%jac(nx+1,j)
        end do
     end do
  end do
  rho(1)=0.0d0 ! boudary = \Phi(0,t) = 0

  open(16,file='mneme')
  write(16,*)"#rho"
  do i=1,nx
     do j=1,ng
        write(16,*)mesh%etat1(i)+(1.0d0+gausslob%node(j))/mesh%jac(i,nv+1), &
             & rho((i-1)*ng+j)
     end do
  end do
  write(16,*)" "

  !rho=M*rho
  do i=1,nx
     do j=1,ng
        rho((i-1)*ng+j)=rho((i-1)*ng+j)*gausslob%weigh(j)*mesh%jac(i,nv+1)
     end do
  end do

  sys=UMFPACK_A

  call umf4def(control)
  control(1)=2
  call umf4sym (int(nx*ng,umf_int),int(nx*ng,umf_int),matvp%ap,matvp%ai,matvp%ax,symbolic,control,info)
  if (info(1) .lt. 0) then
     print *, 'Error occurred in umf4sym: ', info (1)
     stop
  end if
  !call umf4pinf(control,info);stop
  call umf4num (matvp%ap,matvp%ai,matvp%ax,symbolic,numeric,control,info)
  call umf4solr(sys,matvp%ap,matvp%ai,matvp%ax,phi,rho,numeric,control,info)

  write(16,*)' '
  write(16,*)"#phi"
  do i=1,nx
     do j=1,ng
        write(16,*)mesh%etat1(i)+(1.0d0+gausslob%node(j))/mesh%jac(i,nv+1), &
             & phi((i-1)*ng+j)
     end do
  end do
  write(16,*)" "

  field_e=matmul(fieldvp,phi)

  write(16,*)' '
  write(16,*)"#E"
  do i=1,nx
     do j=1,ng
        write(16,*)mesh%etat1(i)+(1.0d0+gausslob%node(j))/mesh%jac(i,nv+1), &
             & field_e((i-1)*ng+j)
     end do
  end do
  write(16,*)" "

  !construction or rhs
  !rhs=-v.d_x(f)+d_x(Phi).d_v(f), with d=\partial
  rhs=0.0d0
  do x1=1,nx !loop on elements in direction x
     do v1=1,nv !loop on elements in direction v
        do x2=1,ng !loop on GLL points in direction x
           do v2=1,ng !loop on GLL points in direction v
              !interior part
              som1=0.0d0
              som2=0.0d0
              do i=1,ng
                 !possible mistake on indexes on der
                 !if wrong, transpose indexes on der
                 som1=som1+dist((x1-1)*ng+i,(v1-1)*ng+v2)*gausslob%der(x2,i)
                 som2=som2+dist((x1-1)*ng+x2,(v1-1)*ng+i)*gausslob%der(v2,i)
              end do
              som1=som1*mesh%jac(x1,nv+1)
              som2=som2*mesh%jac(nx+1,v1)
              rhs((x1-1)*ng+x2,(v1-1)*ng+v2)=-(som1*gausslob%weigh(v2)* &
                   & (mesh%etat2(v1)+(1.0d0+gausslob%node(v2))/mesh%jac(nx+1,v1))- &
                   & som2*gausslob%weigh(x2)*field_e((x1-1)*ng+x2))/mesh%jac(x1,v1)

!!$              !boudaries part
!!$              if (x2==ng) then 
!!$                 if (mesh%etat2(v1)+(1.0d0+gausslob%node(v2))/mesh%jac(nx+1,v1) >= 0.0d0) then
!!$                    rhs((x1-1)*ng+x2,(v1-1)*ng+v2)=rhs((x1-1)*ng+x2,(v1-1)*ng+v2)+ &
!!$                         & (mesh%etat2(v1)+(1.0d0+gausslob%node(v2))/mesh%jac(nx+1,v1))* &
!!$                         & gausslob%weigh(v2)*dist(x1*ng,(v1-1)*ng+v2)/mesh%jac(nx+1,v1) &
!!$                         & *gausslob%weigh(x2)/mesh%jac(x1,nv+1)/4.0d0 ! added for test but wrong
!!$                 else
!!$                    rhs((x1-1)*ng+x2,(v1-1)*ng+v2)=rhs((x1-1)*ng+x2,(v1-1)*ng+v2)+ &
!!$                         & (mesh%etat2(v1)+(1.0d0+gausslob%node(v2))/mesh%jac(nx+1,v1))* &
!!$                         & gausslob%weigh(v2)*dist(modulo(x1*ng+1-1,nx*ng)+1,(v1-1)*ng+v2)/ &
!!$                         & mesh%jac(nx+1,v1) &
!!$                         & *gausslob%weigh(x2)/mesh%jac(x1,nv+1)/4.0d0
!!$                 end if
!!$              else if (x2==1) then
!!$                 if (mesh%etat2(v1)+(1.0d0+gausslob%node(v2))/mesh%jac(nx+1,v1) >= 0.0d0) then
!!$                    rhs((x1-1)*ng+x2,(v1-1)*ng+v2)=rhs((x1-1)*ng+x2,(v1-1)*ng+v2)- &
!!$                         & (mesh%etat2(v1)+(1.0d0+gausslob%node(v2))/mesh%jac(nx+1,v1))* &
!!$                         & gausslob%weigh(v2)*dist(modulo((x1-1)*ng-1,nx*ng)+1,(v1-1)*ng+v2)/ &
!!$                         & mesh%jac(nx+1,v1) &
!!$                         & *gausslob%weigh(x2)/mesh%jac(x1,nv+1)/4.0d0
!!$                 else
!!$                    rhs((x1-1)*ng+x2,(v1-1)*ng+v2)=rhs((x1-1)*ng+x2,(v1-1)*ng+v2)- &
!!$                         & (mesh%etat2(v1)+(1.0d0+gausslob%node(v2))/mesh%jac(nx+1,v1))* &
!!$                         & gausslob%weigh(v2)*dist((x1-1)*ng+1,(v1-1)*ng+v2)/mesh%jac(nx+1,v1) &
!!$                         & *gausslob%weigh(x2)/mesh%jac(x1,nv+1)/4.0d0
!!$                 end if
!!$              end if
!!$              if (v2==ng) then
!!$                 if (field_e((x1-1)*ng+x2) >= 0.0d0) then
!!$                    if (v1<nv) then
!!$                       rhs((x1-1)*ng+x2,(v1-1)*ng+v2)=rhs((x1-1)*ng+x2,(v1-1)*ng+v2)- &
!!$                            & gausslob%weigh(x2)*field_e((x1-1)*ng+x2)* &
!!$                            & dist((x1-1)*ng+x2,v1*ng+1)/mesh%jac(x1,nv+1) &
!!$                            & *gausslob%weigh(v2)/mesh%jac(nx+1,v1)/4.0d0
!!$                    end if
!!$                 else
!!$                    rhs((x1-1)*ng+x2,(v1-1)*ng+v2)=rhs((x1-1)*ng+x2,(v1-1)*ng+v2)- &
!!$                         & gausslob%weigh(x2)*field_e((x1-1)*ng+x2)* &
!!$                         & dist((x1-1)*ng+x2,v1*ng)/mesh%jac(x1,nv+1) &
!!$                         & *gausslob%weigh(v2)/mesh%jac(nx+1,v1)/4.0d0
!!$                 end if
!!$              else if (v2==1) then
!!$                 if (field_e((x1-1)*ng+x2) >= 0.0d0) then
!!$                    rhs((x1-1)*ng+x2,(v1-1)*ng+v2)=rhs((x1-1)*ng+x2,(v1-1)*ng+v2)+ &
!!$                         & gausslob%weigh(x2)*field_e((x1-1)*ng+x2)* &
!!$                         & dist((x1-1)*ng+x2,(v1-1)*ng+1)/mesh%jac(x1,nv+1) &
!!$                         & *gausslob%weigh(v2)/mesh%jac(nx+1,v1)/4.0d0
!!$                 else 
!!$                    if (v1>=2) then
!!$                       rhs((x1-1)*ng+x2,(v1-1)*ng+v2)=rhs((x1-1)*ng+x2,(v1-1)*ng+v2)+ &
!!$                            & gausslob%weigh(x2)*field_e((x1-1)*ng+x2)* &
!!$                            & dist((x1-1)*ng+x2,(v1-1)*ng)/mesh%jac(x1,nv+1) &
!!$                            & *gausslob%weigh(v2)/mesh%jac(nx+1,v1)/4.0d0
!!$                    end if
!!$                 end if
!!$              end if

              rhs((x1-1)*ng+x2,(v1-1)*ng+v2)=rhs((x1-1)*ng+x2,(v1-1)*ng+v2)*mesh%jac(x1,v1) &
                   & /(gausslob%weigh(x2)*gausslob%weigh(v2))

           end do
        end do
     end do
  end do

  write(16,*)' '
  write(16,*)"#rhs dist"
  do v1=1,nv
     do v2=1,ng
        do x1=1,nx
           do x2=1,ng
              write(16,*)mesh%etat1(x1)+(1.0d0+gausslob%node(x2))/mesh%jac(x1,nv+1), &
                   & mesh%etat2(v1)+(1.0d0+gausslob%node(v2))/mesh%jac(nx+1,v1), &
                   & rhs((x1-1)*ng+x2,(v1-1)*ng+v2), dist((x1-1)*ng+x2,(v1-1)*ng+v2)
           end do
        end do
        write(16,*)' '
     end do
  end do
  close(16)

!!$  open(12,file='clio')
!!$  call write_octave(rhs,'rhs',12)
!!$  close(12)
!!$
!!$  linf=0.0d0
!!$  l1=0.0d0
!!$  l2=0.0d0
!!$
!!$  do i=1,nx
!!$     do j=1,ng
!!$        linf=max(linf,abs(phi((i-1)*ng+j)-rho((i-1)*ng+j)))
!!$        l1=l1+abs(phi((i-1)*ng+j)-rho((i-1)*ng+j))*mesh%jac(i,nv+1)*gausslob%weigh(j)
!!$        l2=l2+abs(phi((i-1)*ng+j)-rho((i-1)*ng+j))**2*mesh%jac(i,nv+1)*gausslob%weigh(j)
!!$     end do
!!$  end do
!!$  l2=sqrt(l2)
!!$
!!$  open(14,file='aoede',position='append')
!!$  write(14,'(1A1,9X,1A2,10X,1A2,2X,1A4,22X,1A2,23X,1A2)') &
!!$       & "#","ne","ng","linf","l1","l2"
!!$  write(14,*)nx,ng,linf,l1,l2
!!$  close(14)
!!$
!!$  open(13,file='thalie')
!!$  call write_octave(pts,'r','X',13)
!!$  call write_octave(phi,'r','PHI',13)
!!$  call write_octave(rho,'r','RHO',13)
!!$  close(13)
!!$
!!$  open(11,file='olympe')
!!$  call write_octave(matvp,'vp',11)
!!$  call write_octave(matvm,'vm',11)
!!$  call write_octave(rho,'r','RHO',11)
!!$  !call write_octave(new_tri(nx*ng,1,(/(i-1,i=1,nx*ng)/),(/(0,i=1,nx*ng)/),pts),'x',11)
!!$  close(11)

  deallocate(rho,phi,field_e,dist,rhs)
  call delete_gausslobatto_1D(gausslob)
  call clear(matvm)
  call clear(matvp)
  call clear(fieldvm)
  call clear(fieldvp)
  call delete_nu_cart_mesh(mesh)

end program VP_DG
