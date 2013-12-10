!>Solve Poisson equation on cartesian domain with finit elements.
!> * Compact boundary conditions.
!> * Linear system solve with lapack (Choleski)
module sll_poisson_2d_periodic_fem
#include "sll_poisson_solvers_macros.h"
#include "sll_working_precision.h"
#include "sll_memory.h"
use sll_constants
implicit none

!> Structure to solve Poisson equation on 2d irregular cartesian mesh
!> with finite element numerical method
type :: poisson_2d_periodic_fem
   sll_real64, dimension(:,:), pointer :: A     !< Mass matrix
   sll_real64, dimension(:,:), pointer :: M     !< Stiffness matrix
   sll_real64, dimension(:)  , pointer :: hx    !< step size x
   sll_real64, dimension(:)  , pointer :: hy    !< step size y
   sll_int32,  dimension(:)  , pointer :: ipiv  !< Lapack array for pivoting
end type poisson_2d_periodic_fem

!> Initialize the solver
interface initialize
   module procedure initialize_poisson_2d_periodic_fem
end interface initialize
!> Compute the electric potential
interface solve
   module procedure solve_poisson_2d_periodic_fem
end interface solve

sll_int32, private :: nx, ny, nxy
sll_int32, private :: i, j, k, ii, jj
sll_int32, private :: error

private :: som, build_matrices

contains

!> Initialize Poisson solver object using finite elements method.
!> Indices are shifted from \f$ [1:n+1] \f$ to \f$ [0:n] \f$ only 
!> inside this subroutine.
subroutine initialize_poisson_2d_periodic_fem( this, x, y ,nn_x, nn_y)
type( poisson_2d_periodic_fem ) :: this !< Solver data structure
sll_int32,  intent(in)          :: nn_x !< number of cells along x
sll_int32,  intent(in)          :: nn_y !< number of cells along y
sll_real64, dimension(nn_x)     :: x    !< x nodes coordinates
sll_real64, dimension(nn_y)     :: y    !< y nodes coordinates

sll_real64, dimension(4,4) :: Axelem
sll_real64, dimension(4,4) :: Ayelem
sll_real64, dimension(4,4) :: Melem
sll_int32, dimension(4) :: isom

nx = nn_x-1
ny = nn_y-1
nxy = nx * ny

call write_mtv_periodic( x, y )

SLL_ALLOCATE(this%hx(1:nx),error)
SLL_ALLOCATE(this%hy(1:ny),error)
SLL_ALLOCATE(this%A(nxy,nxy), error)
SLL_ALLOCATE(this%M(nxy,nxy), error)

do i=1,nx
   this%hx(i) = x(i+1)-x(i)
end do

do j=1,ny
   this%hy(j) = y(j+1)-y(j)
end do

!** Construction des matrices elementaires
Axelem(1,:) = (/  2, -2, -1,  1 /)
Axelem(2,:) = (/ -2,  2,  1, -1 /)
Axelem(3,:) = (/ -1,  1,  2, -2 /)
Axelem(4,:) = (/  1, -1, -2,  2 /)

Axelem = 1.d0/6.d0 * Axelem

Ayelem(1,:) = (/  2,  1, -1, -2 /)
Ayelem(2,:) = (/  1,  2, -2, -1 /)
Ayelem(3,:) = (/ -1, -2,  2,  1 /)
Ayelem(4,:) = (/ -2, -1,  1,  2 /)

Ayelem = 1.d0/6.d0 * Ayelem

Melem(1,:) = (/ 4, 2, 1, 2/)
Melem(2,:) = (/ 2, 4, 2, 1/)
Melem(3,:) = (/ 1, 2, 4, 2/)
Melem(4,:) = (/ 2, 1, 2, 4/)

Melem = 1.d0/36.d0 * Melem

this%A = 0.d0

!***  Interior mesh ***
do j=1,ny-1
   do i=1,nx-1
      isom(1) =  som(i,j,1)
      isom(2) =  som(i,j,2)
      isom(3) =  som(i,j,3)
      isom(4) =  som(i,j,4)
      call build_matrices( this, Axelem, Ayelem, Melem, isom, i, j )
   end do
end do

do i=2,nx-1

   j = 1
   isom(1)=som(i,ny,1)
   isom(2)=som(i,ny,2)
   isom(3)=som(i,j,2)
   isom(4)=som(i,j,1)
   call build_matrices( this, Axelem, Ayelem, Melem, isom, i, j )

end do

do j=2,ny-1

   i = 1
   isom(1)=som(nx,j,1)
   isom(2)=som(i,j,1)
   isom(3)=som(i,j,4)
   isom(4)=som(nx,j,4)
   call build_matrices( this, Axelem, Ayelem, Melem, isom, i, j )

end do
   
!Corners
i=1; j=1  !SW
isom(1) = som(nx-1,ny-1,3)
isom(2) = som(i   ,ny-1,4)
isom(3) = som(i   ,j   ,1)
isom(4) = som(nx-1,j   ,2)
call build_matrices( this, Axelem, Ayelem, Melem, isom, i, j )

SLL_ALLOCATE(this%ipiv(nxy),error)

this%A(1,:) = 0.
this%A(1,1) = 1.
call DGETRF(nxy,nxy,this%A,nxy,this%ipiv,error)

end subroutine initialize_poisson_2d_periodic_fem

!> Get the node number
integer function som(i, j, k)
integer :: i, j, k

if (k == 1) then
   som = i+(j-1)*nx
else if (k == 2) then
   som = i+(j-1)*nx+1
else if (k == 3) then
   som = i+j*nx+1
else if (k == 4) then
   som = i+j*nx
end if 

end function som

!> Build matrices and factorize
subroutine build_matrices( this, Axelem, Ayelem, Melem, isom, i, j )
type( poisson_2d_periodic_fem ) :: this     !< Poisson solver object
sll_real64, dimension(:,:)      :: Axelem   !< x electric field
sll_real64, dimension(:,:)      :: Ayelem   !< y electric field
sll_real64, dimension(:,:)      :: Melem    !< charge density
sll_int32,  dimension(:)        :: isom     !< node indices
sll_int32, intent(in)           :: i        !< int(x) position on mesh
sll_int32, intent(in)           :: j        !< int(y) position on mesh

   do ii=1,4
      do jj=1,4
         this%A(isom(ii),isom(jj)) = this%A(isom(ii),isom(jj)) &
                 & + Axelem(ii,jj) * this%hy(j) / this%hx(i) &
                 & + Ayelem(ii,jj) * this%hx(i) / this%hy(j)
         this%M(isom(ii),isom(jj)) = this%M(isom(ii),isom(jj)) &
                 & + Melem(ii,jj)  * this%hx(i) * this%hy(j)
      end do
   end do

end subroutine build_matrices

!> Solve the poisson equation
subroutine solve_poisson_2d_periodic_fem( this, ex, ey, rho )
type( poisson_2d_periodic_fem ) :: this !< Poisson solver object
sll_real64, dimension(:,:)      :: ex   !< x electric field
sll_real64, dimension(:,:)      :: ey   !< y electric field
sll_real64, dimension(:,:)      :: rho  !< charge density
sll_real64, dimension(nxy)      :: b
sll_real64                      :: bmoy

!** Construction du second membre (rho a support compact --> projete)
k=0
do i=1,nx
   do j=1,ny
      k=k+1
      b(k) = rho(i,j)
   end do
end do

b = matmul(this%M,b)
b(1) = 1

call DGETRS('N',nxy,1,this%A,nxy,this%ipiv,b,nxy,error)

bmoy = sum(b) / nxy

k=0
do i=1,nx
   do j=1,ny
      k=k+1
      rho(i,j) = b(k) - bmoy
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

end subroutine solve_poisson_2d_periodic_fem

!> Write the Plotmtv file to plot mesh indices and cell numbers.
subroutine write_mtv_periodic( x, y )
integer               :: iel
real(8), dimension(:) :: x !< x node position
real(8), dimension(:) :: y !< y node position
real(8)               :: x1
real(8)               :: y1

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
do j=1,ny-1
   do i=1,nx-1
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
do i=1,nx
   do j=1,ny
      write(10,"(a)"   ,  advance="no")"@text x1="
      write(10,"(g15.3)", advance="no") x(i)
      write(10,"(a)"   ,  advance="no")" y1="
      write(10,"(g15.3)", advance="no") y(j)
      write(10,"(a)"   ,  advance="no")" z1=0. lc=5 ll='"
      write(10,"(i4)"  ,  advance="no") som(i,j,1)
      write(10,"(a)")"'"
   end do
end do
   
write(10,*)"$END"
close(10)

end subroutine write_mtv_periodic


end module sll_poisson_2d_periodic_fem
