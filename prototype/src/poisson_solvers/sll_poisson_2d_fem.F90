!>Solve Poisson equation on cartesian domain with finit elements.
!> * Compact boundary conditions.
!> * Linear system solve with lapack (Choleski)
module sll_poisson_2d_fem
#include "sll_working_precision.h"
#include "sll_memory.h"
use sll_constants
implicit none

type :: poisson_fem
   sll_real64, dimension(:,:), pointer :: A
   sll_real64, dimension(:,:), pointer :: M
   sll_real64, dimension(:,:), pointer :: mat
   sll_real64, dimension(:)  , pointer :: hx    !< step size x
   sll_real64, dimension(:)  , pointer :: hy    !< step size y
end type poisson_fem

interface initialize
   module procedure initialize_poisson_2d_fem
end interface initialize
interface solve
   module procedure solve_poisson_2d_fem
end interface solve

contains

subroutine initialize_poisson_2d_fem( this, x, y ,nx, ny)
type( poisson_fem ) :: this
sll_real64, dimension(-1:) :: x, y
sll_int32 :: i, j, ii, jj
sll_int32 :: error

sll_real64, dimension(4,4) :: Axelem
sll_real64, dimension(4,4) :: Ayelem
sll_real64, dimension(4,4) :: Melem
sll_real64 :: dum
sll_int32 :: Iel
sll_int32, dimension(4) :: Isom
sll_int32 :: nx
sll_int32 :: ny

SLL_ALLOCATE(this%hx(-1:nx),error)
SLL_ALLOCATE(this%hy(-1:ny),error)
SLL_ALLOCATE(this%A((nx-1)*(ny-1),(nx-1)*(ny-1)), error)
SLL_ALLOCATE(this%M((nx-1)*(ny-1),(nx-1)*(ny-1)), error)
SLL_ALLOCATE(this%mat(nx+1,(nx-1)*(ny-1)), error)

do i=0,nx-1
   this%hx(i) = x(i+1)-x(i)
end do

do j=0,ny-1
   this%hy(j) = y(j+1)-y(j)
end do

!Utilisé pour des conditions limites periodiques
!this%hx(nx) = this%hx(0)  
!this%hx(-1) = this%hx(nx-1)
!this%hy(ny) = this%hy(0)
!this%hy(-1) = this%hy(ny-1)
!
!x(-1)   = x(0)  - this%hx(nx-1) 
!x(nx+1) = x(nx) + this%hx(0)
!y(-1)   = y(0)  - this%hy(ny-1)
!y(ny+1) = y(ny) + this%hy(0)

!** Construction des matrices elementaires
dum = 1.d0/6.d0
Axelem(1,1)= 2*dum; Axelem(1,2)=-2*dum; Axelem(1,3)= - dum; Axelem(1,4)=   dum ;
Axelem(2,1)=-2*dum; Axelem(2,2)= 2*dum; Axelem(2,3)=   dum; Axelem(2,4)= - dum ;
Axelem(3,1)= - dum; Axelem(3,2)=   dum; Axelem(3,3)= 2*dum; Axelem(3,4)=-2*dum ;
Axelem(4,1)=   dum; Axelem(4,2)=-  dum; Axelem(4,3)=-2*dum; Axelem(4,4)= 2*dum ;

Ayelem(1,1)= 2*dum; Ayelem(1,2)=   dum; Ayelem(1,3)= - dum; Ayelem(1,4)=-2*dum ;
Ayelem(2,1)=   dum; Ayelem(2,2)= 2*dum; Ayelem(2,3)=-2*dum; Ayelem(2,4)= - dum ;
Ayelem(3,1)= - dum; Ayelem(3,2)=-2*dum; Ayelem(3,3)= 2*dum; Ayelem(3,4)=   dum ;
Ayelem(4,1)=-2*dum; Ayelem(4,2)= - dum; Ayelem(4,3)=   dum; Ayelem(4,4)= 2*dum ;

dum = 1.d0/36.d0
Melem(1,1)=4*dum; Melem(1,2)=2*dum; Melem(1,3)=  dum; Melem(1,4)=2*dum;
Melem(2,1)=2*dum; Melem(2,2)=4*dum; Melem(2,3)=2*dum; Melem(2,4)=  dum;
Melem(3,1)=  dum; Melem(3,2)=2*dum; Melem(3,3)=4*dum; Melem(3,4)=2*dum;
Melem(4,1)=2*dum; Melem(4,2)=  dum; Melem(4,3)=2*dum; Melem(4,4)=4*dum;

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
this%A(1,1) =   this%A(1,1) &
              + Axelem(3,3) * this%hy(0) / this%hx(0) &
              + Ayelem(3,3) * this%hx(0) / this%hy(0)

this%M(1,1) = this%M(1,1) + Melem(3,3) * this%hx(0) * this%hy(0)

Isom(4) = nx-1 !SE
this%A(nx-1,nx-1) =   this%A(nx-1,nx-1) &
                    + Axelem(4,4) * this%hy(0) / this%hx(nx-1) &
                    + Ayelem(4,4) * this%hx(nx-1) / this%hy(0)

this%M(nx-1,nx-1) = this%M(nx-1,nx-1) + Melem(4,4) * this%hx(nx-1) * this%hy(0)

Isom(1) = (nx-1)*(ny-1)   !NE

this%A(Isom(1),Isom(1)) =   this%A(Isom(1),Isom(1)) &
                          + Axelem(1,1) * this%hy(ny-1) / this%hx(nx-1) &
                          + Ayelem(1,1) * this%hx(nx-1) / this%hy(ny-1)

this%M(Isom(1),Isom(1)) =   this%M(Isom(1),Isom(1))  &
                          + Melem(1,1) * this%hx(nx-1) * this%hy(ny-1)

Isom(2) = (nx-1)*(ny-2)+1 !NW

this%A(Isom(2),Isom(2)) =   this%A(Isom(2),Isom(2))  &
                          + Axelem(2,2) * this%hy(ny-1) / this%hx(0) &
                          + Ayelem(2,2) * this%hx(0) / this%hy(ny-1)

this%M(Isom(2),Isom(2)) =   this%M(Isom(2),Isom(2))  &
                          + Melem(2,2) * this%hx(0) * this%hy(ny-1)

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
call dpbtrf('U',(nx-1)*(ny-1),nx,this%mat,nx+1,error)

end subroutine initialize_poisson_2d_fem

subroutine solve_poisson_2d_fem( this, ex, ey, rho, nx, ny )
type( poisson_fem ) :: this
sll_real64, dimension(:,:) :: ex
sll_real64, dimension(:,:) :: ey
sll_real64, dimension(:,:) :: rho
sll_int32 :: i, j, nx, ny
sll_real64, dimension((nx-1)*(ny-1)) :: b
sll_int32 :: error

!** Construction du second membre (rho a support compact --> projete)
do i=1,nx-1
   do j=1,ny-1
      b(i+(j-1)*(nx-1)) = rho(i,j)
   end do
end do

b = matmul(this%M,b)

call dpbtrs('U',(nx-1)*(ny-1),nx,1,this%mat,nx+1,b,(nx-1)*(ny-1),error) 

do i=1,nx-1
   do j=1,ny-1
      rho(i,j) = b(i+(j-1)*(nx-1)) 
   end do
end do

do j=1,ny-1
do i=1,nx-2
   ex(i,j) = - (rho(i+1,j)-rho(i,j)) / this%hx(i)
end do
end do

do j=1,ny-2
do i=1,nx-1
   ey(i,j) = - (rho(i,j+1)-rho(i,j)) / this%hy(j)
end do
end do

end subroutine solve_poisson_2d_fem

end module sll_poisson_2d_fem
