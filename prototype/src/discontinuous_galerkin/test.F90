program test_DG
#include "sll_working_precision.h"

  use sll_nu_cart_mesh
  use gausslobatto
  use poisson4dg
  use mod_sparse
  use sll_constants
  use timestep4dg

  use mod_octave_io_sparse

  implicit none

  !searching why not periodic

  type(non_unif_cart_mesh) :: mesh
  type(gausslobatto1D) :: gll
  sll_int32 :: nx,nv,ng,x1,x2,v1,v2
  sll_real64,dimension(:),allocatable :: field
  sll_real64,dimension(:,:),allocatable :: dist,rhs
  sll_real64 :: x,v
  
  nx=10
  nv=10
  ng=5

  call init_gausslobatto_1d(ng,gll)
  call init_nu_cart_mesh(nx,nv,mesh)
  mesh%d_etat1=40.0d0/real(nx,8)
  mesh%d_etat2=40.0d0/real(nv,8)
  call fill_node_nuc_mesh(-20.0d0,-20.0d0,mesh)
  allocate(field(nx*ng),rhs(nx*ng,nv*nv),dist(nx*ng,nv*nv))

  do x1=1,nx
     do x2=1,ng
        x=mesh%etat1(x1)+(1.0d0+gll%node(x2))/mesh%jac(x1,nv+1)
        do v1=1,nv
           do v2=1,ng
              v=mesh%etat2(v1)+(1.0d0+gll%node(v2))/mesh%jac(nx+1,v1)
              dist((x1-1)*ng+x2,(v1-1)*ng+v2)=sin(x)*sin(v)
           end do
        end do
        field((x1-1)*ng+x2)=2.0d0*cos(x)
     end do
  end do

  open(12,file='thalie')
  call rhs4dg_1d(mesh,gll,field,dist,rhs)
    do v1=1,nv !loop on elements in direction v
       do v2=1,ng !loop on GLL points in direction v
          v=mesh%etat2(v1)+(1.0d0+gll%node(v2))/mesh%jac(nx+1,v1)
          do x1=1,nx !loop on elements in direction x
             do x2=1,ng !loop on GLL points in direction x
                x=mesh%etat1(x1)+(1.0d0+gll%node(x2))/mesh%jac(x1,nv+1)
                write(12,*)x,v,rhs((x1-1)*ng+x2,(v1-1)*ng+v2)
             end do
          end do
          write(12,*)''
       end do
    end do
  close(12)
  
  call clear(gll)
  call clear(mesh)
  deallocate(field,rhs,dist)

!test for Blanca : computation of integral

!!$  type(non_unif_cart_mesh) :: mesh
!!$  type(gausslobatto1D) :: gausslob
!!$  integer :: nx,nv,ng,x1,x2
!!$  sll_real64 :: x,res,f
!!$
!!$  ng=30
!!$  nx=500
!!$  nv=1
!!$
!!$  call init_gausslobatto_1d(ng,gausslob)
!!$  call init_nu_cart_mesh(nx,nv,mesh)
!!$  mesh%d_etat1=8.0d0/real(nx,8)
!!$  mesh%d_etat2=1.0d0/real(nv,8)
!!$  call fill_node_nuc_mesh(-4.0d0,0.0d0,mesh)
!!$
!!$  res=0.0d0
!!$
!!$  do x1=1,nx
!!$     do x2=1,ng
!!$        x=mesh%etat1(x1)+(1.0d0+gausslob%node(x2))/mesh%jac(x1,2)
!!$        f=exp(-0.25d0*(4.0d0*x-1.0d0)**2)
!!$        res=res+f*gausslob%weigh(x2)/mesh%jac(x1,2)
!!$     end do
!!$  end do
!!$
!!$  print*,res
!!$  print*,sqrt(sll_pi)/2.0d0

!test of 1./x*x for timestep4dg.F90
  
!!$  type(non_unif_cart_mesh) :: mesh
!!$  type(gausslobatto1D) :: gausslob
!!$  sll_int32 :: x1,x2,v1,v2
!!$  sll_real64,dimension(30*6,30*6) :: rhs
!!$  sll_int32 :: nx,nv,ng
!!$  sll_real64 :: x_min,x_max,v_min,v_max
!!$
!!$  nv=30
!!$  nx=30
!!$  ng=6
!!$
!!$  x_min=-sll_pi
!!$  x_max=sll_pi
!!$  v_min=-sll_pi
!!$  v_max=sll_pi
!!$
!!$
!!$  call init_gausslobatto_1d(ng,gausslob)
!!$  !call init_gausslobatto_1d(3*ng,gll)
!!$  call init_nu_cart_mesh(nx,nv,mesh)
!!$  mesh%d_etat1=(x_max-x_min)/real(nx,8)
!!$  mesh%d_etat2=(v_max-v_min)/real(nv,8)
!!$  call fill_node_nuc_mesh(x_min,v_min,mesh)
!!$  
!!$  rhs=1.0d0
!!$  do v1=1,nv !loop on elements in direction v
!!$     do v2=1,ng !loop on GLL points in direction v
!!$        do x1=1,nx !loop on elements in direction x
!!$           do x2=1,ng !loop on GLL points in direction x
!!$              rhs((x1-1)*ng+x2,(v1-1)*ng+v2)=rhs((x1-1)*ng+x2,(v1-1)*ng+v2)* &
!!$                   & gausslob%weigh(x2)*gausslob%weigh(v2)/mesh%jac(x1,v1)
!!$              rhs((x1-1)*ng+x2,(v1-1)*ng+v2)=rhs((x1-1)*ng+x2,(v1-1)*ng+v2)* &
!!$                   & mesh%jac(x1,v1)/(gausslob%weigh(x2)*gausslob%weigh(v2))
!!$           end do
!!$        end do
!!$     end do
!!$  end do
!!$
!!$  print*,maxval(abs(rhs))
!!$  print*,minval(abs(rhs))

!local max of der

!!$  type(gausslobatto1D) :: gausslob
!!$
!!$  call init_gausslobatto_1d(10,gausslob)
!!$  open(13,file='erato')
!!$  call write_octave(gausslob%der,'der',13)
!!$  close(13)
!!$  call delete_gausslobatto_1d(gausslob)

!check GLL nodes and weigh

!!$  integer :: i,j
!!$  type(gausslobatto1D) :: gausslob
!!$
!!$  print *, 'Exact value: '
!!$  write (*,'(e22.15)') 0.4674011002723395
!!$
!!$  print*,'Test Gauss-Lobatto'
!!$  do i=2,10
!!$     call init_gausslobatto_1d(i,gausslob)
!!$     !don't to it in real program, use transformation between real mesh and reference element
!!$     !here it is done for simplicity
!!$     gausslob%node(:)=(gausslob%node(:)+1.0d0)*(sll_pi/2.0d0)/2.0d0
!!$     write (*,'(a, i8, a, e20.12)') 'case n = ', i, ': ', &
!!$          & sum((/ (gausslob%weigh(j)*test_func(gausslob%node(j))*sll_pi/4.0d0,j=1,i) /))
!!$     call delete_gausslobatto_1d(gausslob)
!!$  end do
!!$
!!$  print*,'Test Gauss-lobatto points and weight (5 points)'
!!$  call init_gausslobatto_1d(5,gausslob)
!!$  print*,gausslob%node
!!$  print*,gausslob%weigh
!!$  call delete_gausslobatto_1d(gausslob)
!!$
!!$contains
!!$
!!$  function test_func(x)
!!$    intrinsic :: dcos
!!$    sll_real64 :: test_func
!!$    sll_real64, intent(in) :: x
!!$    test_func = x*x*dcos(x)
!!$  end function test_func

!max(w)=O(ng^a)

!!$  type(gausslobatto1D) :: gausslob
!!$  sll_real64 :: maxw,maxd
!!$  sll_int32 :: i
!!$
!!$  open(13,file='erato')
!!$
!!$  do i=2,10
!!$     call init_gausslobatto_1d(i,gausslob)
!!$     !maxw=maxval(gausslob%weigh)
!!$     gausslob%der=abs(gausslob%der)
!!$     !gausslob%der(1,1)=0.0d0
!!$     !gausslob%der(i,i)=0.0d0
!!$     !maxd=maxval(gausslob%der)
!!$     !write(13,*)i,maxw,maxd
!!$     write(13,*)i
!!$     write(13,*)gausslob%der
!!$     call delete_gausslobatto_1d(gausslob)
!!$  end do
!!$
!!$  close(13)

!transposition and derivative matrix

!!$  type(t_tri) :: d1,d2,d3
!!$  sll_int32 :: i,j,ne,ng,k
!!$  type(gausslobatto1D) :: gausslob
!!$  sll_real64 :: c12!,som
!!$
!!$  ne=5
!!$  ng=6
!!$
!!$  !c12=0.5d0
!!$  c12=0.0d0
!!$
!!$  call init_gausslobatto_1d(ng,gausslob)
!!$  d1=new_tri(ne*ng,ne*ng,ne*(ng**2+2))
!!$  d2=new_tri(ne*ng,ne*ng,ne*(ng**2+2))
!!$
!!$  print*,'der'
!!$  do i=1,ng
!!$     som=0.0d0
!!$     do j=1,ng
!!$        som=som+gausslob%der(i,j)
!!$     end do
!!$     print*,j,som
!!$  end do
!!$  print*,'end der'
!!$
!!$  d1%ti=2
!!$  d1%tj=2
!!$  d2%ti=2
!!$  d2%tj=2
!!$  d1%tx=0.0d0
!!$  d2%tx=0.0d0
!!$
!!$  do i=1,ne
!!$     !line
!!$     do j=1,ng
!!$        !column
!!$        do k=1,ng
!!$           d1%ti((i-1)*ng**2+(j-1)*ng+k)=(i-1)*ng+j
!!$           d1%tj((i-1)*ng**2+(j-1)*ng+k)=(i-1)*ng+k
!!$           d1%tx((i-1)*ng**2+(j-1)*ng+k)=gausslob%der(j,k)
!!$           d1%tx((i-1)*ng**2+(j-1)*ng+k)=gausslob%der(k,j)
!!$
!!$           d2%ti((i-1)*ng**2+(j-1)*ng+k)=(i-1)*ng+j
!!$           d2%tj((i-1)*ng**2+(j-1)*ng+k)=(i-1)*ng+k
!!$           d2%tx((i-1)*ng**2+(j-1)*ng+k)=gausslob%der(j,k)
!!$           d2%tx((i-1)*ng**2+(j-1)*ng+k)=-gausslob%der(j,k)
!!$           d2%tx((i-1)*ng**2+(j-1)*ng+k)=-gausslob%der(k,j)
!!$           
!!$           if (j==1 .and. k==1) then
!!$              d1%tx((i-1)*ng**2+(j-1)*ng+k)=d1%tx((i-1)*ng**2+(j-1)*ng+k)-(-0.5d0-c12)
!!$              d2%tx((i-1)*ng**2+(j-1)*ng+k)=d2%tx((i-1)*ng**2+(j-1)*ng+k)-0.5d0+c12
!!$              d1%tx((i-1)*ng**2+(j-1)*ng+k)=d1%tx((i-1)*ng**2+(j-1)*ng+k)+1.0d0
!!$              d2%tx((i-1)*ng**2+(j-1)*ng+k)=d2%tx((i-1)*ng**2+(j-1)*ng+k)-1.0d
!!$              d1%tx((i-1)*ng**2+(j-1)*ng+k)=d1%tx((i-1)*ng**2+(j-1)*ng+k)-(-0.5d0+c12)
!!$              d2%tx((i-1)*ng**2+(j-1)*ng+k)=d2%tx((i-1)*ng**2+(j-1)*ng+k)-0.5d0-c12
!!$           else if (j==ng .and. k==ng) then
!!$              d1%tx((i-1)*ng**2+(j-1)*ng+k)=d1%tx((i-1)*ng**2+(j-1)*ng+k)-(0.5d0+c12)
!!$              d2%tx((i-1)*ng**2+(j-1)*ng+k)=d2%tx((i-1)*ng**2+(j-1)*ng+k)+0.5d0-c12
!!$              d1%tx((i-1)*ng**2+(j-1)*ng+k)=d1%tx((i-1)*ng**2+(j-1)*ng+k)-1.0d0
!!$              d2%tx((i-1)*ng**2+(j-1)*ng+k)=d2%tx((i-1)*ng**2+(j-1)*ng+k)+1.0d0
!!$           end if
!!$        end do
!!$     end do
!!$  end do
!if (.false.) then
!!$  d1%ti(ne*ng**2+1)=ng
!!$  d1%tj(ne*ng**2+1)=ng+1
!!$  d1%tx(ne*ng**2+1)=-(0.5d0-c12)
!!$  d1%tx(ne*ng**2+1)=1.0d0
!!$
!!$  d1%ti(ne*(ng**2+2))=1
!!$  d1%tj(ne*(ng**2+2))=ne*ng
!!$  d1%tx(ne*(ng**2+2))=-(-0.5d0+c12)
!!$  d1%tx(ne*(ng**2+2))=-1.0d0
!!$  d1%tx(ne*(ng**2+2))=-(-0.5d0-c12)
!!$
!!$  d2%ti(ne*ng**2+1)=ng
!!$  d2%tj(ne*ng**2+1)=ng+1
!!$  d2%tx(ne*ng**2+1)=0.5d0+c12
!!$  d2%tx(ne*ng**2+1)=-1.0d0
!!$
!!$  d2%ti(ne*(ng**2+2))=1
!!$  d2%tj(ne*(ng**2+2))=ne*ng
!!$  d2%tx(ne*(ng**2+2))=-0.5d0-c12
!!$  d2%tx(ne*(ng**2+2))=1.0d0
!!$  d2%tx(ne*(ng**2+2))=-0.5d0+c12
!!$
!!$  do i=2,ne-1
!!$     d1%ti(ne*ng**2+i)=i*ng
!!$     d1%tj(ne*ng**2+i)=i*ng+1
!!$     d1%tx(ne*ng**2+i)=-(0.5d0-c12)
!!$     d1%tx(ne*ng**2+i)=1.0d0
!!$
!!$     d1%ti(ne*(ng**2+2)-i+1)=(i-1)*ng+1
!!$     d1%tj(ne*(ng**2+2)-i+1)=(i-1)*ng
!!$     d1%tx(ne*(ng**2+2)-i+1)=-(-0.5d0+c12)
!!$     d1%tx(ne*(ng**2+2)-i+1)=-1.0d
!!$     d1%tx(ne*(ng**2+2)-i+1)=-(-0.5d0-c12)
!!$
!!$     d2%ti(ne*ng**2+i)=i*ng
!!$     d2%tj(ne*ng**2+i)=i*ng+1
!!$     d2%tx(ne*ng**2+i)=0.5d0+c12
!!$     d2%tx(ne*ng**2+i)=-1.0d0
!!$
!!$     d2%ti(ne*(ng**2+2)-i+1)=(i-1)*ng+1
!!$     d2%tj(ne*(ng**2+2)-i+1)=(i-1)*ng
!!$     d2%tx(ne*(ng**2+2)-i+1)=-0.5d0-c12
!!$     d2%tx(ne*(ng**2+2)-i+1)=1.0d0
!!$     d2%tx(ne*(ng**2+2)-i+1)=-0.5d0+c12
!!$
!!$  end do
!!$
!!$  d1%ti(ne*ng**2+ne)=ne*ng
!!$  d1%tj(ne*ng**2+ne)=1
!!$  d1%tx(ne*ng**2+ne)=-(0.5d0-c12)
!!$  d1%tx(ne*ng**2+ne)=1.0d0
!!$
!!$  d1%ti(ne*ng**2+ne+1)=(ne-1)*ng+1
!!$  d1%tj(ne*ng**2+ne+1)=(ne-1)*ng
!!$  d1%tx(ne*ng**2+ne+1)=-(-0.5d0+c12)
!!$  d1%tx(ne*ng**2+ne+1)=-1.0d0
!!$  d1%tx(ne*ng**2+ne+1)=-(-0.5d0-c12)
!!$
!!$  d2%ti(ne*ng**2+ne)=ne*ng
!!$  d2%tj(ne*ng**2+ne)=1
!!$  d2%tx(ne*ng**2+ne)=0.5d0+c12
!!$  d2%tx(ne*ng**2+ne)=-1.0d0
!!$
!!$  d2%ti(ne*ng**2+ne+1)=(ne-1)*ng+1
!!$  d2%tj(ne*ng**2+ne+1)=(ne-1)*ng
!!$  d2%tx(ne*ng**2+ne+1)=-0.5d0-c12
!!$  d2%tx(ne*ng**2+ne+1)=1.0d0
!!$  d2%tx(ne*ng**2+ne+1)=-0.5d0+c12
!!$!end if
!!$  d1%ti=d1%ti-1
!!$  d1%tj=d1%tj-1
!!$  d2%ti=d2%ti-1
!!$  d2%tj=d2%tj-1
!!$
!!$  d2%tx=-d2%tx
!!$  d3=col2tri(tri2col(d2+transpose(d1)))
!!$  d3=col2tri(tri2col(d2+d1))
!!$
!!$  do i=1,d3%nz
!!$     if (d3%tx(i)<=0.001d0) then
!!$        d3%tx(i)=0.0d0
!!$     end if
!!$  end do
!!$
!!$  do i=1,d1%nz
!!$     if (d1%ti(i)==0) then
!!$        print*,d1%tj(i),d1%tx(i)
!!$     end if
!!$  end do
!!$
!!$  do i=1,d3%nz
!!$     print*,d3%ti(i)+1,d3%tj(i)+1,d3%tx(i)
!!$  end do
!!$
!!$  print*,d1%n,d1%m
!!$  print*,matmul((/(1.0d0,i=1,ne*ng)/),transpose(tri2col(d1)))
!!$  print*,matmul(tri2col(d1),(/(1.0d0,i=1,ne*ng)/))
!!$
!!$  open(12,file='calliope')
!!$  call write_octave(tri2col(d1),'d1',12)
!!$  call write_octave(tri2col(d2),'d2',12)
!!$  call write_octave(tri2col(d3),'d3',12)
!!$  close(12)

!mesh

!!$  !rho,E,phi
!!$  sll_real64,dimension(:),allocatable :: rho,field,phi
!!$  !distribution function
!!$  sll_real64,dimension(:,:),allocatable :: dist
!!$  !number of time steb, number of cells un direction x and v
!!$  sll_int32 :: nx,nv,nb_step
!!$  !mesh
!!$  type(non_unif_cart_mesh) :: mesh
!!$  !coefficients for fluxes
!!$  sll_real64 :: c11,c12,c21,c22
!!$
!!$  sll_int32 :: i,j
!!$
!!$  nx=30
!!$  nv=40
!!$
!!$  call init_nu_cart_mesh(nx,nv,mesh,d_etat1=0.1d0)
!!$
!!$  do i=1,nx
!!$     mesh%d_etat1(i)=exp(-0.1d0*real(nx/2-i,kind(1.0d0))**2)
!!$  end do
!!$  do i=1,nv
!!$     mesh%d_etat2(i)=exp(-0.1d0*real(nv/2-i,kind(1.0d0))**2)
!!$  end do
!!$  call fill_node_nuc_mesh(-1.0d0,-1.0d0,mesh)
!!$
!!$  open(12,file='mesh.dat')
!!$  do i=1,nx+1
!!$     do j=1,nv+1
!!$        write(12,*)mesh%etat1(i),mesh%etat2(j)
!!$     end do
!!$  end do
!!$  close(12)
!!$
!!$  call delete_nu_cart_mesh(mesh)

!derivative matrix

!!$  type(gausslobatto1d) :: gl
!!$  sll_int32 :: i
!!$  sll_int32,parameter :: n=6
!!$  sll_real64,dimension(n,n) :: d
!!$
!!$  call init_gausslobatto_1d(n,gl)
!!$  print*,gl%node
!!$
!!$  call derivative_matrix_1d(gl,d)
!!$  print*,'d'
!!$  do i=1,n
!!$     print*,d(:,i)
!!$  end do
  
end program test_DG
