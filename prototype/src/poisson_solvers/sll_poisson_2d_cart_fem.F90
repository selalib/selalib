program sll_poisson_2d_cart_fem
#include "sll_working_precision.h"
#include "sll_memory.h"
use numeric_constants
implicit none

sll_int32 :: i, j, ii, jj
sll_int32, parameter :: nx = 64, ny = 64
sll_real64 :: dimx, dimy

type :: poisson_fem
sll_real64, dimension((nx-1)*(ny-1),(nx-1)*(ny-1)) :: A, M
sll_real64, dimension(nx+1,(nx-1)*(ny-1)) :: mat
sll_real64, dimension(:), allocatable :: hx, hy    ! les h_i+1/2
end type poisson_fem

call test()

contains

subroutine test()
type( poisson_fem ) :: poisson
sll_real64, dimension(0:nx,0:ny) :: ex
sll_real64, dimension(0:nx,0:ny) :: ey
sll_real64, dimension(0:nx,0:ny) :: phi
sll_real64, dimension(0:nx,0:ny) :: rho
sll_real64, dimension(:), pointer :: x, y
sll_real64 :: dx, dy
sll_int32 :: error, mode

SLL_ALLOCATE(x(-1:nx+1),error)  
SLL_ALLOCATE(y(-1:ny+1),error) 

dimx = 2 * sll_pi
dimy = 2 * sll_pi

dx = dimx / nx
dy = dimy / ny

x(0) = 0.
y(0) = 0.

do i=1,nx
   x(i) = (i*dx) *(i*dx+1)/(1+dimx)
enddo
do j=1,ny
   y(j) = (j*dy) *(j*dy+1)/(1+dimy)
enddo


mode = 2
do i = 1, nx
   do j = 1, ny
      phi(i,j) = mode * sin(mode*x(i)) * sin(mode*y(j))
      rho(i,j) = 2_f64 * mode**3 * sin(mode*x(i))*sin(mode*y(j))
      write(10,*) x(i), y(j), rho(i,j)
   end do
   write(10,*)
end do
close(10)

call initialize(poisson, x, y)
call solve(poisson, ex, ey, rho)

do i = 1, nx
   do j = 1, ny
      write(11,*) x(i), y(j), ex(i,j),  mode**2*cos(mode*x(i))*sin(mode*y(j))
      write(12,*) x(i), y(j), ey(i,j), -mode**2*sin(mode*x(i))*cos(mode*y(j))
   end do
   write(11,*)
   write(12,*)
end do
close(11)
close(12)

end subroutine test


subroutine initialize( this, x, y )
type( poisson_fem ) :: this
sll_real64, dimension(-1:) :: x, y
sll_int32 :: i, j
sll_int32 :: error

sll_real64, dimension(4,4) :: Axelem, Ayelem, Melem
sll_real64 :: dum
sll_int32 :: Iel, info
sll_int32, dimension(4) :: Isom


SLL_ALLOCATE(this%hx(-1:nx),error)
SLL_ALLOCATE(this%hy(-1:ny),error)

do i=0,nx-1
   this%hx(i) = x(i+1)-x(i)
end do

do j=0,ny-1
   this%hy(j) = y(j+1)-y(j)
end do

this%hx(nx) = this%hx(0)  ! CL periodiques
this%hx(-1) = this%hx(nx-1)
this%hy(ny) = this%hy(0)
this%hy(-1) = this%hy(ny-1)

x(-1)   = x(0) - this%hx(nx-1)  !points utiles pour le cas period
x(nx+1) = x(nx) + this%hx(0)
y(-1)   = y(0) - this%hy(ny-1)
y(ny+1) = y(ny) + this%hy(0)

!** Construction des matrices elementaires
dum = 1.d0/6.d0
Axelem(1,1)= 2*dum ; Axelem(1,2)=-2*dum ; Axelem(1,3)= - dum ; Axelem(1,4)=   dum ;
Axelem(2,1)=-2*dum ; Axelem(2,2)= 2*dum ; Axelem(2,3)=   dum ; Axelem(2,4)= - dum ;
Axelem(3,1)= - dum ; Axelem(3,2)=   dum ; Axelem(3,3)= 2*dum ; Axelem(3,4)=-2*dum ;
Axelem(4,1)=   dum ; Axelem(4,2)=-  dum ; Axelem(4,3)=-2*dum ; Axelem(4,4)= 2*dum ;

Ayelem(1,1)= 2*dum ; Ayelem(1,2)=   dum ; Ayelem(1,3)= - dum ; Ayelem(1,4)=-2*dum ;
Ayelem(2,1)=   dum ; Ayelem(2,2)= 2*dum ; Ayelem(2,3)=-2*dum ; Ayelem(2,4)= - dum ;
Ayelem(3,1)= - dum ; Ayelem(3,2)=-2*dum ; Ayelem(3,3)= 2*dum ; Ayelem(3,4)=   dum ;
Ayelem(4,1)=-2*dum ; Ayelem(4,2)= - dum ; Ayelem(4,3)=   dum ; Ayelem(4,4)= 2*dum ;

dum = 1.d0/36.d0
Melem(1,1)=4*dum ; Melem(1,2)=2*dum ; Melem(1,3)=  dum ; Melem(1,4)=2*dum ;
Melem(2,1)=2*dum ; Melem(2,2)=4*dum ; Melem(2,3)=2*dum ; Melem(2,4)=  dum ;
Melem(3,1)=  dum ; Melem(3,2)=2*dum ; Melem(3,3)=4*dum ; Melem(3,4)=2*dum ;
Melem(4,1)=2*dum ; Melem(4,2)=  dum ; Melem(4,3)=2*dum ; Melem(4,4)=4*dum ;


this%A = 0.d0

!** Contribution des mailles interieures
do i=1,nx-2
   do j=1,ny-2
      Iel = i+(j-1)*(nx-1)
      Isom(1)=Iel; Isom(2)=Iel+1; Isom(3)=Iel+nx; Isom(4)=Iel+nx-1;
      do ii=1,4
         do jj=1,4
            this%A(Isom(ii),Isom(jj)) = this%A(Isom(ii),Isom(jj)) &
                 & + Axelem(ii,jj) * this%hy(j) / this%hx(i) &
                 & + Ayelem(ii,jj) * this%hx(i) / this%hy(j)
            this%M(Isom(ii),Isom(jj)) = this%M(Isom(ii),Isom(jj)) &
                 & + Melem(ii,jj) * this%hx(i) * this%hy(j)
         end do
      end do
   end do
end do

!** Contribution des mailles au sud et au nord
do i=1,nx-2
   Isom(3)=i+1; Isom(4)=i  !Sud
   do ii=3,4
      do jj=3,4
         this%A(Isom(ii),Isom(jj)) = this%A(Isom(ii),Isom(jj)) &
                 & + Axelem(ii,jj) * this%hy(0) / this%hx(i) &
                 & + Ayelem(ii,jj) * this%hx(i) / this%hy(0)
         this%M(Isom(ii),Isom(jj)) = this%M(Isom(ii),Isom(jj)) &
              & + Melem(ii,jj) * this%hx(i) * this%hy(0)
      end do
   end do
   Iel = (ny-2)*(nx-1)+i   !Nord
   Isom(1)=Iel; Isom(2)=Iel+1
   do ii=1,2
      do jj=1,2
         this%A(Isom(ii),Isom(jj)) = this%A(Isom(ii),Isom(jj)) &
                 & + Axelem(ii,jj) * this%hy(ny-1) / this%hx(i) &
                 & + Ayelem(ii,jj) * this%hx(i) / this%hy(ny-1)
         this%M(Isom(ii),Isom(jj)) = this%M(Isom(ii),Isom(jj)) &
              & + Melem(ii,jj) * this%hx(i) * this%hy(ny-1)
      end do
   end do
end do

!** Contribution des mailles a l'ouest et a l'est
do j=1,ny-2
   Isom(2)=1+(j-1)*(nx-1); Isom(3)=1+j*(nx-1) !Ouest
   do ii=2,3
      do jj=2,3
         this%A(Isom(ii),Isom(jj)) = this%A(Isom(ii),Isom(jj)) &
                 & + Axelem(ii,jj) * this%hy(j) / this%hx(0) &
                 & + Ayelem(ii,jj) * this%hx(0) / this%hy(j)
         this%M(Isom(ii),Isom(jj)) = this%M(Isom(ii),Isom(jj)) &
              & + Melem(ii,jj) * this%hx(0) * this%hy(j)
      end do
   end do
   Iel = j*(nx-1)                              !Est
   Isom(1)=Iel; Isom(4)=Iel+nx-1
   do ii=1,4,3
      do jj=1,4,3
         this%A(Isom(ii),Isom(jj)) = this%A(Isom(ii),Isom(jj)) &
                 & + Axelem(ii,jj) * this%hy(j) / this%hx(nx-1) &
                 & + Ayelem(ii,jj) * this%hx(nx-1) / this%hy(j)
         this%M(Isom(ii),Isom(jj)) = this%M(Isom(ii),Isom(jj)) &
              & + Melem(ii,jj) * this%hx(nx-1) * this%hy(j)
      end do
   end do
end do

!** Contribution des coins
Isom(3) = 1    !SW
this%A(1,1) = this%A(1,1) + Axelem(3,3) * this%hy(0) / this%hx(0) + Ayelem(3,3) * this%hx(0) / this%hy(0)
this%M(1,1) = this%M(1,1) + Melem(3,3) * this%hx(0) * this%hy(0)

Isom(4) = nx-1 !SE
this%A(nx-1,nx-1) = this%A(nx-1,nx-1) + Axelem(4,4) * this%hy(0) / this%hx(nx-1) &
     & + Ayelem(4,4) * this%hx(nx-1) / this%hy(0)
this%M(nx-1,nx-1) = this%M(nx-1,nx-1) + Melem(4,4) * this%hx(nx-1) * this%hy(0)

Isom(1) = (nx-1)*(ny-1)   !NE
this%A(Isom(1),Isom(1)) = this%A(Isom(1),Isom(1)) + Axelem(1,1) * this%hy(ny-1) / this%hx(nx-1) &
     & + Ayelem(1,1) * this%hx(nx-1) / this%hy(ny-1)
this%M(Isom(1),Isom(1)) = this%M(Isom(1),Isom(1)) + Melem(1,1) * this%hx(nx-1) * this%hy(ny-1)

Isom(2) = (nx-1)*(ny-2)+1 !NW
this%A(Isom(2),Isom(2)) = this%A(Isom(2),Isom(2)) + Axelem(2,2) * this%hy(ny-1) / this%hx(0) &
     & + Ayelem(2,2) * this%hx(0) / this%hy(ny-1)
this%M(Isom(2),Isom(2)) = this%M(Isom(2),Isom(2)) + Melem(2,2) * this%hx(0) * this%hy(ny-1)

this%mat = 0.d0
this%mat(nx+1,1) = this%A(1,1)
do j=2,(nx-1)*(ny-1)
   this%mat(nx+1,j) = this%A(j,j)
   this%mat(nx,j)   = this%A(j-1,j)
end do
this%mat(3,nx) = this%A(2,nx)
this%mat(2,nx) = this%A(1,nx)
do j=nx+1,(nx-1)*(ny-1)
   this%mat(3,j) = this%A(j-nx+2,j)
   this%mat(2,j) = this%A(j-nx+1,j)
   this%mat(1,j) = this%A(j-nx,j)
end do
CALL DPBTRF('U',(nx-1)*(ny-1),nx,this%mat,nx+1,info)
print*,'factorisation pour Cholesky',info

end subroutine initialize

subroutine solve( this, ex, ey, rho )
type( poisson_fem ) :: this
sll_real64, dimension(0:nx,0:ny) :: ex
sll_real64, dimension(0:nx,0:ny) :: ey
sll_real64, dimension(0:nx,0:ny) :: phi
sll_real64, dimension(0:nx,0:ny) :: rho
sll_real64, dimension((nx-1)*(ny-1)) :: b
sll_int32 :: info

!** Construction du second membre (rho a support compact --> projete)
do i=1,nx-1
   do j=1,ny-1
      b(i+(j-1)*(nx-1)) = rho(i,j)
   end do
end do

b = matmul(this%M,b)

CALL DPBTRS('U',(nx-1)*(ny-1),nx,1,this%mat,nx+1,b,(nx-1)*(ny-1),info) 
print*,'resolution par Cholesky',info

phi = 0.d0
do i=1,nx-1
   do j=1,ny-1
      phi(i,j) = b(i+(j-1)*(nx-1)) 
   end do
end do

!** Reconstruction du champ E

do i=0,nx-1
   do j=0,ny
      ex(i,j) = - (phi(i+1,j)-phi(i,j)) / this%hx(i)
   end do
end do

do i=0,nx
   do j=0,ny-1
      ey(i,j) = - (phi(i,j+1)-phi(i,j)) / this%hy(j)
   end do
end do

end subroutine solve

end program sll_poisson_2d_cart_fem
