!>Solve Poisson equation on cartesian domain with finit elements.
!> * Compact boundary conditions.
!> * Linear system solve with lapack (Choleski)
module sll_poisson_2d_fem
#include "sll_poisson_solvers_macros.h"
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

sll_int32, private :: nx, ny
sll_int32, private :: i, j

contains

!> Initialize Poisson solver object using finite elements method.
!> Indices are shifted from [1:n+1] to [0:n] only inside this 
!> subroutine
subroutine initialize_poisson_2d_fem( this, x, y ,nn_x, nn_y, bc_x, bc_y)
type( poisson_fem ) :: this
sll_int32,  intent(in)      :: nn_x !< number of cells along x
sll_int32,  intent(in)      :: nn_y !< number of cells along y
sll_real64, dimension(nn_x) :: x    !< x nodes coordinates
sll_real64, dimension(nn_y) :: y    !< y nodes coordinates
sll_int32, optional         :: bc_x !< boundary condition type in x direction
sll_int32, optional         :: bc_y !< boundary condition type in x direction
sll_int32 :: ii, jj
sll_int32 :: error

sll_real64, dimension(4,4) :: Axelem
sll_real64, dimension(4,4) :: Ayelem
sll_real64, dimension(4,4) :: Melem
sll_real64 :: dum
sll_int32 :: iel, info
sll_int32, dimension(4) :: isom

nx = nn_x-1
ny = nn_y-1

SLL_ALLOCATE(this%hx(1:nx),error)
SLL_ALLOCATE(this%hy(1:ny),error)
SLL_ALLOCATE(this%A((nx-1)*(ny-1),(nx-1)*(ny-1)), error)
SLL_ALLOCATE(this%M((nx-1)*(ny-1),(nx-1)*(ny-1)), error)
SLL_ALLOCATE(this%mat(nx+1,(nx-1)*(ny-1)), error)

do i=1,nx
   this%hx(i) = x(i+1)-x(i)
end do

do j=1,ny
   this%hy(j) = y(j+1)-y(j)
end do

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
      iel = i+(j-1)*(nx-1)
      isom(1)=iel; isom(2)=iel+1; isom(3)=iel+nx; isom(4)=iel+nx-1;
      do ii=1,4
         do jj=1,4
            this%A(isom(ii),isom(jj)) = this%A(isom(ii),isom(jj)) &
                 & + Axelem(ii,jj) * this%hy(j) / this%hx(i) &
                 & + Ayelem(ii,jj) * this%hx(i) / this%hy(j)
            this%M(isom(ii),isom(jj)) = this%M(isom(ii),isom(jj)) &
                 & + Melem(ii,jj) * this%hx(i) * this%hy(j)
         end do
      end do
   end do
end do

if (present(bc_y) .and. bc_y == COMPACT) then
   do i=1,nx-2
      isom(3)=i+1; isom(4)=i  !Sud
      do ii=3,4
         do jj=3,4
            this%A(isom(ii),isom(jj)) = this%A(isom(ii),isom(jj)) &
                    & + Axelem(ii,jj) * this%hy(1) / this%hx(i) &
                    & + Ayelem(ii,jj) * this%hx(i) / this%hy(1)
            this%M(isom(ii),isom(jj)) = this%M(isom(ii),isom(jj)) &
                 & + Melem(ii,jj) * this%hx(i) * this%hy(1)
         end do
      end do
      iel = (ny-2)*(nx-1)+i   !Nord
      isom(1)=iel; isom(2)=iel+1
      do ii=1,2
         do jj=1,2
            this%A(isom(ii),isom(jj)) = this%A(isom(ii),isom(jj)) &
                    & + Axelem(ii,jj) * this%hy(ny-1) / this%hx(i) &
                    & + Ayelem(ii,jj) * this%hx(i) / this%hy(ny-1)
            this%M(isom(ii),isom(jj)) = this%M(isom(ii),isom(jj)) &
                 & + Melem(ii,jj) * this%hx(i) * this%hy(ny-1)
         end do
      end do
   end do
end if

if (present(bc_x) .and. bc_x == COMPACT) then
   do j=1,ny-2
      isom(2)=1+(j-1)*(nx-1); isom(3)=1+j*(nx-1) !Ouest
      do ii=2,3
         do jj=2,3
            this%A(isom(ii),isom(jj)) = this%A(isom(ii),isom(jj)) &
                    & + Axelem(ii,jj) * this%hy(j) / this%hx(1) &
                    & + Ayelem(ii,jj) * this%hx(1) / this%hy(j)
            this%M(isom(ii),isom(jj)) = this%M(isom(ii),isom(jj)) &
                 & + Melem(ii,jj) * this%hx(1) * this%hy(j)
         end do
      end do
      iel = j*(nx-1)                              !Est
      isom(1)=iel; isom(4)=iel+nx-1
      do ii=1,4,3
         do jj=1,4,3
            this%A(isom(ii),isom(jj)) = this%A(isom(ii),isom(jj)) &
                    & + Axelem(ii,jj) * this%hy(j) / this%hx(nx-1) &
                    & + Ayelem(ii,jj) * this%hx(nx-1) / this%hy(j)
            this%M(isom(ii),isom(jj)) = this%M(isom(ii),isom(jj)) &
                 & + Melem(ii,jj) * this%hx(nx-1) * this%hy(j)
         end do
      end do
   end do
end if
   
!** Contribution des coins
isom(3) = 1    !SW
this%A(1,1) =   this%A(1,1) &
              + Axelem(3,3) * this%hy(1) / this%hx(1) &
              + Ayelem(3,3) * this%hx(1) / this%hy(1)

this%M(1,1) = this%M(1,1) + Melem(3,3) * this%hx(1) * this%hy(1)

isom(4) = nx-1 !SE
this%A(nx-1,nx-1) =   this%A(nx-1,nx-1) &
                    + Axelem(4,4) * this%hy(1) / this%hx(nx-1) &
                    + Ayelem(4,4) * this%hx(nx-1) / this%hy(1)

this%M(nx-1,nx-1) = this%M(nx-1,nx-1) + Melem(4,4) * this%hx(nx-1) * this%hy(1)

isom(1) = (nx-1)*(ny-1)   !NE

this%A(isom(1),isom(1)) =   this%A(isom(1),isom(1)) &
                          + Axelem(1,1) * this%hy(ny-1) / this%hx(nx-1) &
                          + Ayelem(1,1) * this%hx(nx-1) / this%hy(ny-1)

this%M(isom(1),isom(1)) =   this%M(isom(1),isom(1))  &
                          + Melem(1,1) * this%hx(nx-1) * this%hy(ny-1)

isom(2) = (nx-1)*(ny-2)+1 !NW

this%A(isom(2),isom(2)) =   this%A(isom(2),isom(2))  &
                          + Axelem(2,2) * this%hy(ny-1) / this%hx(1) &
                          + Ayelem(2,2) * this%hx(1) / this%hy(ny-1)

this%M(isom(2),isom(2)) =   this%M(isom(2),isom(2))  &
                          + Melem(2,2) * this%hx(1) * this%hy(ny-1)

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
call dpbtrf('U',(nx-1)*(ny-1),nx,this%mat,nx+1,info)

call write_mtv_file( x, y)

end subroutine initialize_poisson_2d_fem

!> Solve the poisson equation
subroutine solve_poisson_2d_fem( this, ex, ey, rho )
type( poisson_fem )        :: this !< Poisson solver object
sll_real64, dimension(:,:) :: ex   !< x electric field
sll_real64, dimension(:,:) :: ey   !< y electric field
sll_real64, dimension(:,:) :: rho  !< charge density
sll_real64, dimension((nx-1)*(ny-1)) :: b
sll_int32 :: info

!** Construction du second membre (rho a support compact --> projete)
do i=2,nx
   do j=2,ny
      b((i-1)+(j-2)*(nx-1)) = rho(i,j)
   end do
end do

b = matmul(this%M,b)

call dpbtrs('U',(nx-1)*(ny-1),nx,1,this%mat,nx+1,b,(nx-1)*(ny-1),info) 

do i=2,nx
   do j=2,ny
      rho(i,j) = b((i-1)+(j-2)*(nx-1)) 
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

subroutine write_mtv_file( x, y )
sll_real64, dimension(:) :: x
sll_real64, dimension(:) :: y
sll_int32 :: iel
sll_real64 :: x1, y1

open(10, file="gd.mtv")
write(10,*)"$DATA=CURVE3D"
write(10,*)"%equalscale=T"
write(10,*)"%toplabel='Numero des faces ' "
   
do i=1,nx-2
   do j=1,ny-2
      write(10,*) x(i  ), y(j  ), 0.
      write(10,*) x(i+1), y(j  ), 0.
      write(10,*) x(i+1), y(j+1), 0.
      write(10,*) x(i  ), y(j+1), 0.
      write(10,*) x(i  ), y(j  ), 0.
      write(10,*)
   end do
end do

!Numeros des elements
do i=1,nx-2
   do j=1,ny-2
      iel = i+(j-1)*(nx-1)
      x1 = 0.5*(x(i)+x(i+1))
      y1 = 0.5*(y(j)+y(j+1))
      write(10,"(a)"   ,  advance="no")"@text x1="
      write(10,"(g15.3)", advance="no") x1
      write(10,"(a)"   ,  advance="no")" y1="
      write(10,"(g15.3)", advance="no") y1
      write(10,"(a)"   ,  advance="no")" z1=0. lc=4 ll='"
      write(10,"(i4)"  ,  advance="no") iel
      write(10,"(a)")"'"
   end do
end do
   
write(10,*)"$END"
close(10)

end subroutine write_mtv_file

end module sll_poisson_2d_fem
