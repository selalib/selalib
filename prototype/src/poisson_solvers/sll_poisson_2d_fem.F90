!>Solve Poisson equation on cartesian domain with finite elements.
!> * Compact boundary conditions.
!> * Linear system solve with lapack (Choleski)
module sll_poisson_2d_fem
#include "sll_poisson_solvers_macros.h"
#include "sll_working_precision.h"
#include "sll_memory.h"
use sll_constants
implicit none

!> Structure to solve Poisson equation on 2d domain. Mesh is cartesian and
!> could be irregular. Numerical method is using finite elements.
type :: poisson_fem
   sll_real64, dimension(:,:), pointer :: A   !< Mass matrix
   sll_real64, dimension(:,:), pointer :: M   !< Stiffness matrix
   sll_real64, dimension(:,:), pointer :: mat !< Matrix solve by Lapack
   sll_real64, dimension(:)  , pointer :: hx  !< step size x
   sll_real64, dimension(:)  , pointer :: hy  !< step size y
end type poisson_fem

!> Initialize the solver 
interface initialize
   module procedure initialize_poisson_2d_fem
end interface initialize
!> Compute the potential
interface solve
   module procedure solve_poisson_2d_fem
end interface solve

sll_int32, private :: nx, ny
sll_int32, private :: i, j, k
sll_int32, private :: error

contains

!> Initialize Poisson solver object using finite elements method.
!> Indices are shifted from [1:n+1] to [0:n] only inside this 
!> subroutine
subroutine initialize_poisson_2d_fem( this, x, y ,nn_x, nn_y)
type( poisson_fem ) :: this         !< solver data structure
sll_int32,  intent(in)      :: nn_x !< number of cells along x
sll_int32,  intent(in)      :: nn_y !< number of cells along y
sll_real64, dimension(nn_x) :: x    !< x nodes coordinates
sll_real64, dimension(nn_y) :: y    !< y nodes coordinates
sll_int32 :: ii, jj

sll_real64, dimension(4,4) :: Axelem
sll_real64, dimension(4,4) :: Ayelem
sll_real64, dimension(4,4) :: Melem
sll_real64 :: dum
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

!***  Interior mesh ***
do i=2,nx-1
   do j=2,ny-1
      do ii=1,4
         do jj=1,4
            this%A(som(i,j,ii),som(i,j,jj)) = this%A(som(i,j,ii),som(i,j,jj)) &
                    & + Axelem(ii,jj) * this%hy(j) / this%hx(i)   &
                    & + Ayelem(ii,jj) * this%hx(i) / this%hy(j)
            this%M(som(i,j,ii),som(i,j,jj)) = this%M(som(i,j,ii),som(i,j,jj)) &
                    & + Melem(ii,jj)  * this%hx(i) * this%hy(j)
         end do
      end do
   end do
end do

call write_mtv_file( x, y)

do i=2,nx-1
   j = 1
   isom(3)=som(i,j+1,1)
   isom(4)=som(i,j+1,2)
   do ii=3,4
      do jj=3,4
         this%A(isom(ii),isom(jj)) = this%A(isom(ii),isom(jj)) &
                 & + Axelem(ii,jj) * this%hy(j) / this%hx(i)   &
                 & + Ayelem(ii,jj) * this%hx(i) / this%hy(j)
         this%M(isom(ii),isom(jj)) = this%M(isom(ii),isom(jj)) &
                 & + Melem(ii,jj)  * this%hx(i) * this%hy(j)
      end do
   end do
   j = ny
   isom(1)=som(i,j-1,3)
   isom(2)=som(i,j-1,4)
   do ii=1,2
      do jj=1,2
         this%A(isom(ii),isom(jj)) = this%A(isom(ii),isom(jj)) &
                 & + Axelem(ii,jj) * this%hy(j) / this%hx(i)   &
                 & + Ayelem(ii,jj) * this%hx(i) / this%hy(j)
         this%M(isom(ii),isom(jj)) = this%M(isom(ii),isom(jj)) &
                 & + Melem(ii,jj)  * this%hx(i) * this%hy(j)
      end do
   end do
end do

do j=2,ny-1
   i = 1
   isom(2)=som(i+1,j,1)
   isom(3)=som(i+1,j,4)
   do ii=2,3
      do jj=2,3
         this%A(isom(ii),isom(jj)) = this%A(isom(ii),isom(jj)) &
                 & + Axelem(ii,jj) * this%hy(j) / this%hx(i) &
                 & + Ayelem(ii,jj) * this%hx(i) / this%hy(j)
         this%M(isom(ii),isom(jj)) = this%M(isom(ii),isom(jj)) &
                 & + Melem(ii,jj)  * this%hx(i) * this%hy(j)
      end do
   end do
   i = nx
   isom(1)=som(i-1,j,2)
   isom(4)=som(i-1,j,3)
   do ii=1,4,3
      do jj=1,4,3
         this%A(isom(ii),isom(jj)) = this%A(isom(ii),isom(jj)) &
                 & + Axelem(ii,jj) * this%hy(j) / this%hx(i) &
                 & + Ayelem(ii,jj) * this%hx(i) / this%hy(j)
         this%M(isom(ii),isom(jj)) = this%M(isom(ii),isom(jj)) &
                 & + Melem(ii,jj)  * this%hx(i) * this%hy(j)
      end do
   end do
end do

!Corners
i=1; j=1  !SW
isom(3) = som(i+1,j+1,1)
this%A(isom(3),isom(3)) = this%A(i,j)                           &
                        + Axelem(3,3) * this%hy(j) / this%hx(i) &
                        + Ayelem(3,3) * this%hx(i) / this%hy(j)
this%M(isom(3),isom(3)) = this%M(isom(3),isom(3))               &
                        + Melem(3,3) * this%hx(i) * this%hy(j)

i=nx; j=1 !SE
isom(4) = som(i-1,j+1,2)
this%A(isom(4),isom(4)) = this%A(isom(4),isom(4))               &
                        + Axelem(4,4) * this%hy(i) / this%hx(i) &
                        + Ayelem(4,4) * this%hx(j) / this%hy(j)
this%M(isom(4),isom(4)) = this%M(isom(4),isom(4))               &
                        + Melem(4,4) * this%hx(i) * this%hy(j)

i=nx; j=ny !NE
isom(1) = som(i-1,j-1,3)
this%A(isom(1),isom(1)) = this%A(isom(1),isom(1))               &
                        + Axelem(1,1) * this%hy(j) / this%hx(i) &
                        + Ayelem(1,1) * this%hx(i) / this%hy(j)
this%M(isom(1),isom(1)) =   this%M(isom(1),isom(1))             &
                        + Melem(1,1) * this%hx(i) * this%hy(j)

i=1; j=ny !NW
isom(2) = som(i+1,j-1,4) 
this%A(isom(2),isom(2)) = this%A(isom(2),isom(2))               &
                        + Axelem(2,2) * this%hy(j) / this%hx(i) &
                        + Ayelem(2,2) * this%hx(i) / this%hy(j)
this%M(isom(2),isom(2)) =   this%M(isom(2),isom(2))             &
                        + Melem(2,2) * this%hx(i) * this%hy(j)

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

!> Get the node index
integer function som(i, j, k)

   integer :: i, j, k

   if (k == 1) then
      som = i-1+(j-2)*(nx-1)
   else if (k == 2) then
      som = i-1+(j-2)*(nx-1)+1
   else if (k == 3) then
      som = i-1+(j-2)*(nx-1)+nx
   else if (k == 4) then
      som = i-1+(j-2)*(nx-1)+nx-1
   end if 

end function som

!> Solve the poisson equation
subroutine solve_poisson_2d_fem( this, ex, ey, rho )
type( poisson_fem )        :: this !< Poisson solver object
sll_real64, dimension(:,:) :: ex   !< x electric field
sll_real64, dimension(:,:) :: ey   !< y electric field
sll_real64, dimension(:,:) :: rho  !< charge density
sll_real64, dimension((nx-1)*(ny-1)) :: b

!** Construction du second membre (rho a support compact --> projete)
k = 0
do i=2,nx
   do j=2,ny
      k = k+1
      b(k) = rho(i,j)
   end do
end do

b = matmul(this%M,b)

call dpbtrs('U',(nx-1)*(ny-1),nx,1,this%mat,nx+1,b,(nx-1)*(ny-1),error) 

rho = 0.0
k = 0
do i=2,nx
   do j=2,ny
      k = k+1
      rho(i,j) = b(k) 
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

!> Write the Plotmtv file to plot the irregular mesh with node number and cell number
subroutine write_mtv_file( x, y )
real(8), dimension(:) :: x
real(8), dimension(:) :: y
integer :: iel, isom, nx, ny
real(8) :: x1, y1

nx = size(x)-1
ny = size(y)-1
open(10, file="mesh.mtv")
write(10,*)"$DATA=CURVE3D"
write(10,*)"%equalscale=T"
write(10,*)"%toplabel='Elements number ' "
   
do i=1,nx
   do j=1,ny
      write(10,*) x(i  ), y(j  ), 0.
      write(10,*) x(i+1), y(j  ), 0.
      write(10,*) x(i+1), y(j+1), 0.
      write(10,*) x(i  ), y(j+1), 0.
      write(10,*) x(i  ), y(j  ), 0.
      write(10,*)
   end do
end do

!Numeros des elements
iel = 0
do i=2,nx-1
   do j=2,ny-1
      iel = iel+1
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

!Numeros des noeud
do i=2,nx
   do j=2,ny
      isom = isom+1
      write(10,"(a)"   ,  advance="no")"@text x1="
      write(10,"(g15.3)", advance="no") x(i)
      write(10,"(a)"   ,  advance="no")" y1="
      write(10,"(g15.3)", advance="no") y(j)
      write(10,"(a)"   ,  advance="no")" z1=0. lc=5 ll='"
      write(10,"(i4)"  ,  advance="no") isom
      write(10,"(a)")"'"
   end do
end do
   
write(10,*)"$END"
close(10)

end subroutine write_mtv_file


end module sll_poisson_2d_fem
