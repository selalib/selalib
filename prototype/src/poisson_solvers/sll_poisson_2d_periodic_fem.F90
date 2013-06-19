!>Solve Poisson equation on cartesian domain with finit elements.
!> * Compact boundary conditions.
!> * Linear system solve with lapack (Choleski)
module sll_poisson_2d_periodic_fem
#include "sll_poisson_solvers_macros.h"
#include "sll_working_precision.h"
#include "sll_memory.h"
use sll_constants
implicit none

type :: poisson_2d_periodic_fem
   sll_real64, dimension(:,:), pointer :: A
   sll_real64, dimension(:,:), pointer :: M
   sll_real64, dimension(:,:), pointer :: mat
   sll_real64, dimension(:)  , pointer :: hx    !< step size x
   sll_real64, dimension(:)  , pointer :: hy    !< step size y
end type poisson_2d_periodic_fem

interface initialize
   module procedure initialize_poisson_2d_periodic_fem
end interface initialize
interface solve
   module procedure solve_poisson_2d_periodic_fem
end interface solve

sll_int32, private :: nx, ny
sll_int32, private :: i, j, k, ii, jj
sll_int32, private :: error

private :: som, build_matrices

contains

!> Initialize Poisson solver object using finite elements method.
!> Indices are shifted from [1:n+1] to [0:n] only inside this 
!> subroutine
subroutine initialize_poisson_2d_periodic_fem( this, x, y ,nn_x, nn_y)
type( poisson_2d_periodic_fem ) :: this
sll_int32,  intent(in)      :: nn_x !< number of cells along x
sll_int32,  intent(in)      :: nn_y !< number of cells along y
sll_real64, dimension(nn_x) :: x    !< x nodes coordinates
sll_real64, dimension(nn_y) :: y    !< y nodes coordinates

sll_real64, dimension(4,4) :: Axelem
sll_real64, dimension(4,4) :: Ayelem
sll_real64, dimension(4,4) :: Melem
sll_real64 :: dum
sll_int32, dimension(4) :: isom

nx = nn_x-1
ny = nn_y-1

call write_mtv_periodic( x, y )

stop

SLL_ALLOCATE(this%hx(1:nx),error)
SLL_ALLOCATE(this%hy(1:ny),error)
SLL_ALLOCATE(this%A(nx*ny,nx*ny), error)
SLL_ALLOCATE(this%M(nx*ny,nx*ny), error)
SLL_ALLOCATE(this%mat(nx+1,nx*ny), error)

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
      isom(1) =  som(i,j,1)
      isom(2) =  som(i,j,2)
      isom(3) =  som(i,j,3)
      isom(4) =  som(i,j,4)
      call build_matrices( this, Axelem, Ayelem, Melem, isom, i, j )
      print*, isom
   end do
end do

do i=2,nx-1

   j = 1
   isom(1)=som(i,ny,1)
   isom(2)=som(i,ny,2)
   isom(3)=som(i,j+1,2)
   isom(4)=som(i,j+1,1)
   call build_matrices( this, Axelem, Ayelem, Melem, isom, i, j )
   print*, isom

   !j = ny
   !isom(1)=som(i,j-1,4)
   !isom(2)=som(i,j-1,3)
   !isom(3)=som(i,2,2)
   !isom(4)=som(i,2,1)
   !print*, isom
   !call build_matrices( this, Axelem, Ayelem, Melem, isom, i, j )

end do

do j=2,ny-1

   i = 1
   isom(1)=som(nx,j,1)
   isom(2)=som(i+1,j,1)
   isom(3)=som(i+1,j,4)
   isom(4)=som(nx,j,4)
   call build_matrices( this, Axelem, Ayelem, Melem, isom, i, j )
   print*, isom

   !i = nx
   !isom(1)=som(i-1,j,2)
   !isom(2)=som(2,j,1)
   !isom(3)=som(2,j,4)
   !isom(4)=som(i-1,j,3)
   !print*, isom
   !call build_matrices( this, Axelem, Ayelem, Melem, isom, i, j )
end do
   
!Corners
i=1; j=1  !SW
isom(1) = som(nx-1,ny-1,3)
isom(2) = som(2   ,ny-1,4)
isom(3) = som(2   ,2   ,1)
isom(4) = som(nx-1,2   ,2)
call build_matrices( this, Axelem, Ayelem, Melem, isom, i, j )
print*, isom

!i=nx; j=1 !SE
!isom(1) = som(nx-1,ny-1,3)
!isom(2) = som(2   ,ny-1,4)
!isom(3) = som(2   ,2   ,1)
!isom(4) = som(nx-1,2   ,2)
!call build_matrices( this, Axelem, Ayelem, Melem, isom, i, j )

!i=nx; j=ny !NE
!isom(1) = som(nx-1,ny-1,3)
!isom(2) = som(2   ,ny-1,4)
!isom(3) = som(2   ,2   ,1)
!isom(4) = som(nx-1,2   ,2)
!call build_matrices( this, Axelem, Ayelem, Melem, isom, i, j )
!
!i=1; j=ny !NW
!isom(1) = som(nx-1,ny-1,3)
!isom(2) = som(2   ,ny-1,4)
!isom(3) = som(2   ,2   ,1)
!isom(4) = som(nx-1,2   ,2)
!call build_matrices( this, Axelem, Ayelem, Melem, isom, i, j )

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

end subroutine initialize_poisson_2d_periodic_fem

integer function som(i, j, k)

   integer :: i, j, k

   if (k == 1) then
      som = i+(j-1)*nx
   else if (k == 2) then
      som = i+(j-1)*nx+1
   else if (k == 3) then
      som = i+(j-1)*nx+nx
   else if (k == 4) then
      som = i+(j-1)*nx+nx-1
   end if 

end function som

subroutine build_matrices( this, Axelem, Ayelem, Melem, isom, i, j )
type( poisson_2d_periodic_fem )     :: this     !< Poisson solver object
sll_real64, dimension(:,:) :: Axelem   !< x electric field
sll_real64, dimension(:,:) :: Ayelem   !< y electric field
sll_real64, dimension(:,:) :: Melem    !< charge density
sll_int32, dimension(:)   :: isom
sll_int32                  :: i
sll_int32                  :: j

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
type( poisson_2d_periodic_fem )        :: this !< Poisson solver object
sll_real64, dimension(:,:)   :: ex   !< x electric field
sll_real64, dimension(:,:)   :: ey   !< y electric field
sll_real64, dimension(:,:)   :: rho  !< charge density
sll_real64, dimension(nx*ny) :: b

!** Construction du second membre (rho a support compact --> projete)
k=0
do i=1,nx
   do j=1,ny
      k=k+1
      b(k) = rho(i,j)
   end do
end do

b = matmul(this%M,b)

call dpbtrs('U',nx*ny,nx,1,this%mat,nx+1,b,nx*ny,error) 

k=0
do i=1,nx
   do j=1,ny
      k=k+1
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

end subroutine solve_poisson_2d_periodic_fem

subroutine write_mtv_periodic( x, y )
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
