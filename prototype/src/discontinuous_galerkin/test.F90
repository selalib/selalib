program test_DG
#include "sll_working_precision.h"

  use sll_nu_cart_mesh
  use gausslobatto
  use mod_sparse

  use mod_octave_io_sparse

  implicit none

!transposition and derivative matrix

  type(t_tri) :: d1,d2,d3
  sll_int32 :: i,j,ne,ng,k
  type(gausslobatto1D) :: gausslob
  sll_real64 :: c12

  ne=5
  ng=6

  !c12=-0.5d0
  c12=0.0d0

  call init_gausslobatto_1d(ng,gausslob)
  d1=new_tri(ne*ng,ne*ng,ne*(ng**2+2))
  d2=new_tri(ne*ng,ne*ng,ne*(ng**2+2))

  d1%ti=2
  d1%tj=2
  d2%ti=2
  d2%tj=2
  d1%tx=0.0d0
  d2%tx=0.0d0
 
  do i=1,ne
     !line
     do j=1,ng
        !column
        do k=1,ng
           d1%ti((i-1)*ng**2+(j-1)*ng+k)=(i-1)*ng+j
           d1%tj((i-1)*ng**2+(j-1)*ng+k)=(i-1)*ng+k
!!$           d1%tx((i-1)*ng**2+(j-1)*ng+k)=gausslob%der(j,k)
           d1%tx((i-1)*ng**2+(j-1)*ng+k)=gausslob%der(k,j)

           d2%ti((i-1)*ng**2+(j-1)*ng+k)=(i-1)*ng+j
           d2%tj((i-1)*ng**2+(j-1)*ng+k)=(i-1)*ng+k
!!$           d2%tx((i-1)*ng**2+(j-1)*ng+k)=gausslob%der(j,k)
           d2%tx((i-1)*ng**2+(j-1)*ng+k)=-gausslob%der(j,k)
           
           if (j==1 .and. k==1) then
!!$              d1%tx((i-1)*ng**2+(j-1)*ng+k)=d1%tx((i-1)*ng**2+(j-1)*ng+k)-(-0.5d0-c12)
!!$              d2%tx((i-1)*ng**2+(j-1)*ng+k)=d2%tx((i-1)*ng**2+(j-1)*ng+k)-0.5d0-c12
              d1%tx((i-1)*ng**2+(j-1)*ng+k)=d1%tx((i-1)*ng**2+(j-1)*ng+k)+1.0d0
              d2%tx((i-1)*ng**2+(j-1)*ng+k)=d2%tx((i-1)*ng**2+(j-1)*ng+k)-1.0d0
           else if (j==ng .and. k==ng) then
!!$              d1%tx((i-1)*ng**2+(j-1)*ng+k)=d1%tx((i-1)*ng**2+(j-1)*ng+k)-(0.5d0+c12)
!!$              d2%tx((i-1)*ng**2+(j-1)*ng+k)=d2%tx((i-1)*ng**2+(j-1)*ng+k)+0.5d0-c12
              d1%tx((i-1)*ng**2+(j-1)*ng+k)=d1%tx((i-1)*ng**2+(j-1)*ng+k)-1.0d0
              d2%tx((i-1)*ng**2+(j-1)*ng+k)=d2%tx((i-1)*ng**2+(j-1)*ng+k)+1.0d0
           end if
        end do
     end do
  end do

  d1%ti(ne*ng**2+1)=ng
  d1%tj(ne*ng**2+1)=ng+1
!!$  d1%tx(ne*ng**2+1)=-(0.5d0-c12)
  d1%tx(ne*ng**2+1)=1.0d0
  d1%ti(ne*(ng**2+2))=1
  d1%tj(ne*(ng**2+2))=ne*ng
!!$  d1%tx(ne*(ng**2+2))=-(-0.5d0+c12)
  d1%tx(ne*(ng**2+2))=-1.0d0

  d2%ti(ne*ng**2+1)=ng
  d2%tj(ne*ng**2+1)=ng+1
!!$  d2%tx(ne*ng**2+1)=0.5d0+c12
  d2%tx(ne*ng**2+1)=-1.0d0
  d2%ti(ne*(ng**2+2))=1
  d2%tj(ne*(ng**2+2))=ne*ng
!!$  d2%tx(ne*(ng**2+2))=-0.5d0-c12
  d2%tx(ne*(ng**2+2))=1.0d0
  do i=2,ne-1
     d1%ti(ne*ng**2+i)=i*ng
     d1%tj(ne*ng**2+i)=i*ng+1
!!$     d1%tx(ne*ng**2+i)=-(0.5d0-c12)
     d1%tx(ne*ng**2+i)=1.0d0
     d1%ti(ne*(ng**2+2)-i+1)=(i-1)*ng+1
     d1%tj(ne*(ng**2+2)-i+1)=(i-1)*ng
!!$     d1%tx(ne*(ng**2+2)-i+1)=-(-0.5d0+c12)
     d1%tx(ne*(ng**2+2)-i+1)=-1.0d0

     d2%ti(ne*ng**2+i)=i*ng
     d2%tj(ne*ng**2+i)=i*ng+1
!!$     d2%tx(ne*ng**2+i)=0.5d0+c12
     d2%tx(ne*ng**2+i)=-1.0d0
     d2%ti(ne*(ng**2+2)-i+1)=(i-1)*ng+1
     d2%tj(ne*(ng**2+2)-i+1)=(i-1)*ng
!!$     d2%tx(ne*(ng**2+2)-i+1)=-0.5d0-c12
     d2%tx(ne*(ng**2+2)-i+1)=1.0d0
  end do

  d1%ti(ne*ng**2+ne)=ne*ng
  d1%tj(ne*ng**2+ne)=1
!!$  d1%tx(ne*ng**2+ne)=-(0.5d0-c12)
  d1%tx(ne*ng**2+ne)=1.0d0
  d1%ti(ne*ng**2+ne+1)=(ne-1)*ng+1
  d1%tj(ne*ng**2+ne+1)=(ne-1)*ng
!!$  d1%tx(ne*ng**2+ne+1)=-(-0.5d0+c12)
  d1%tx(ne*ng**2+ne+1)=-1.0d0

  d2%ti(ne*ng**2+ne)=ne*ng
  d2%tj(ne*ng**2+ne)=1
!!$  d2%tx(ne*ng**2+ne)=0.5d0+c12
  d2%tx(ne*ng**2+ne)=-1.0d0
  d2%ti(ne*ng**2+ne+1)=(ne-1)*ng+1
  d2%tj(ne*ng**2+ne+1)=(ne-1)*ng
!!$  d2%tx(ne*ng**2+ne+1)=-0.5d0-c12
  d2%tx(ne*ng**2+ne+1)=1.0d0

  d1%ti=d1%ti-1
  d1%tj=d1%tj-1
  d2%ti=d2%ti-1
  d2%tj=d2%tj-1

  !d3=col2tri(tri2col(d2+transpose(d1)))
  d3=col2tri(tri2col(d2+d1))

!!$  do i=1,d3%nz
!!$     if (d3%tx(i) /= 0.0d0 )then
!!$        print*,d3%ti(i),d3%tj(i),d3%tx(i)
!!$     end if
!!$  end do

  !print*,d1%n,d1%m
  !print*,matmul((/(1.0d0,i=1,ne*ng)/),transpose(tri2col(d1)))

  open(12,file='file')
  call write_octave(tri2col(d1),'d1',12)
  call write_octave(tri2col(d2),'d2',12)
  call write_octave(tri2col(d4),'d3',12)
  close(12)

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
