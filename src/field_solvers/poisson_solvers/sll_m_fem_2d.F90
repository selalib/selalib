#ifndef DOXYGEN_SHOULD_SKIP_THIS
!> @ingroup poisson_solvers
!> @brief
!> Poisson solver using finite element
!> @details
!> Solve Poisson equation on cartesian domain with finite elements.
!> * Compact boundary conditions.
!> * Linear system solve with lapack (Choleski)
module sll_m_fem_2d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_poisson_solvers_macros.h"

! use F77_lapack, only: &
!   dpbtrf, &
!   dpbtrs

  use sll_m_utilities, only: &
    sll_o_display

  implicit none

  public :: &
    sll_o_create, &
    sll_t_fem_poisson_2d, &
    sll_o_solve

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!> @brief
!> Structure to solve Poisson equation on 2d domain. Mesh is cartesian and
!> could be irregular. 
!> @details
!> finite element numerical method with
!> Compact boundary conditions
type :: sll_t_fem_poisson_2d
   sll_real64, dimension(:,:), pointer :: A   !< Mass matrix
   sll_real64, dimension(:,:), pointer :: M   !< Stiffness matrix
   sll_real64, dimension(:,:), pointer :: mat !< Matrix solve by Lapack
   sll_real64, dimension(:)  , pointer :: hx  !< step size x
   sll_real64, dimension(:)  , pointer :: hy  !< step size y
   sll_int32                           :: nx  !< cells number along x minus 1
   sll_int32                           :: ny  !< cells number along y minus 1
   sll_int32                           :: nd
   sll_int32,  dimension(:)  , pointer :: id
   sll_real64, dimension(:)  , pointer :: xd
   sll_real64, dimension(:)  , pointer :: yd
end type sll_t_fem_poisson_2d

!> Initialize the solver
interface sll_o_create
   module procedure initialize_poisson_2d_fem
end interface sll_o_create

!> Compute the electric potential
interface sll_o_solve
   module procedure solve_poisson_2d_fem
end interface sll_o_solve


contains

!> Initialize Poisson solver object using finite elements method.
subroutine initialize_poisson_2d_fem( this, x, y ,nn_x, nn_y)
type( sll_t_fem_poisson_2d )  :: this !< solver data structure
sll_int32,  intent(in)      :: nn_x !< number of cells along x
sll_int32,  intent(in)      :: nn_y !< number of cells along y
sll_real64, dimension(nn_x) :: x    !< x nodes coordinates
sll_real64, dimension(nn_y) :: y    !< y nodes coordinates
sll_int32                   :: ii
sll_int32                   :: jj

sll_real64, dimension(4,4)  :: Axelem
sll_real64, dimension(4,4)  :: Ayelem
sll_real64, dimension(4,4)  :: Mmelem
sll_real64, dimension(2,2)  :: Mfelem
sll_int32, dimension(4)     :: isom
sll_int32                   :: i
sll_int32                   :: j
sll_int32                   :: k
sll_int32                   :: nx
sll_int32                   :: ny
sll_int32                   :: error

sll_int32                   :: kd
sll_int32                   :: nd
sll_int32                   :: n

this%nx = nn_x-1
this%ny = nn_y-1

nx = this%nx
ny = this%ny

SLL_ALLOCATE(this%hx(1:nx),error)
SLL_ALLOCATE(this%hy(1:ny),error)
SLL_ALLOCATE(this%A((nx+1)*(ny+1),(nx+1)*(ny+1)), error)
SLL_ALLOCATE(this%M((nx+1)*(ny+1),(nx+1)*(ny+1)), error)

do i=1,nx
  this%hx(i) = x(i+1)-x(i)
end do

do j=1,ny
  this%hy(j) = y(j+1)-y(j)
end do

!** Construction des matrices elementaires
Axelem(1,:) = [ 2.0_f64, -2.0_f64, -1.0_f64,  1.0_f64 ]
Axelem(2,:) = [-2.0_f64,  2.0_f64,  1.0_f64, -1.0_f64 ]
Axelem(3,:) = [-1.0_f64,  1.0_f64,  2.0_f64, -2.0_f64 ]
Axelem(4,:) = [ 1.0_f64, -1.0_f64, -2.0_f64,  2.0_f64 ]

call sll_o_display(Axelem, 'f7.3')

Axelem = Axelem/6.0_f64

Ayelem(1,:) = [ 2.0_f64,  1.0_f64, -1.0_f64, -2.0_f64 ]
Ayelem(2,:) = [ 1.0_f64,  2.0_f64, -2.0_f64, -1.0_f64 ]
Ayelem(3,:) = [-1.0_f64, -2.0_f64,  2.0_f64,  1.0_f64 ]
Ayelem(4,:) = [-2.0_f64, -1.0_f64,  1.0_f64,  2.0_f64 ]

call sll_o_display(Ayelem, 'f7.3')

Ayelem = Ayelem/6.0_f64

Mmelem(1,:) = [ 4.0_f64, 2.0_f64, 1.0_f64, 2.0_f64 ]
Mmelem(2,:) = [ 2.0_f64, 4.0_f64, 2.0_f64, 1.0_f64 ]
Mmelem(3,:) = [ 1.0_f64, 2.0_f64, 4.0_f64, 2.0_f64 ]
Mmelem(4,:) = [ 2.0_f64, 1.0_f64, 2.0_f64, 4.0_f64 ]

call sll_o_display(Mmelem, 'f7.3')

Mmelem = Mmelem/36.0_f64

Mfelem(1,:) = [ 2.0_f64, 1.0_f64]
Mfelem(2,:) = [ 1.0_f64, 2.0_f64]

Mfelem = Mfelem / 6.0_f64

this%A = 0.0_f64
this%M = 0.0_f64

!***  Interior mesh ***
do i=1,nx
  do j=1,ny
    isom(1) = som(this%nx,i,j,1)
    isom(2) = som(this%nx,i,j,2)
    isom(3) = som(this%nx,i,j,3)
    isom(4) = som(this%nx,i,j,4)
    do ii=1,4
      do jj=1,4
        this%A(isom(ii),isom(jj)) = this%A(isom(ii),isom(jj)) &
        & + Axelem(ii,jj) * this%hy(j) / this%hx(i)           &
        & + Ayelem(ii,jj) * this%hx(i) / this%hy(j)
        this%M(isom(ii),isom(jj)) = this%M(isom(ii),isom(jj)) &
        & + Mmelem(ii,jj) * this%hy(j) * this%hx(i)
      end do
    end do
  end do
end do

call write_mtv_file( x, y)
!call sll_o_display(this%A, 'f5.2')

!Set nodes dirichlet boundary conditions
this%nd = 2*(nx+ny)
allocate(this%id(this%nd))
nd = 0
do i = 1, nx
  k = i
  nd = nd+1
  this%id(nd) = k
  k = (nx+1)*(ny+1)-i+1
  nd = nd+1
  this%id(nd) = k
end do
do i = nx+1, nx*(ny+1), nx+1
  k = i
  nd = nd+1
  this%id(nd) = k
  k = i+1
  nd = nd+1
  this%id(nd) = k
end do

call sll_o_display(this%id, 'i5')

if (nx<5) call sll_o_display(this%A, 'f5.2')

do i = 1, this%nd
  k = this%id(i)
  this%A(k,k) = 1d20
end do

if (nx < 5) call sll_o_display(this%A, 'f5.2')

kd = nx+2
SLL_ALLOCATE(this%mat(kd+1,(nx+1)*(ny+1)), error)
this%mat = 0.0_f64
!do k = 0, kd   ! Loop over diagonals
!  i = kd+1-k   ! Row number in this%mat
!  do j = (nx+1)*(ny+1)-k, 1+k,-1
!    this%mat(i,j+k) = this%A(j+k,j) 
!  end do
!end do
n = (nx+1)*(ny+1)
do i = 1, n
  do j=i,min(n,i+kd)
   this%mat(kd+1+i-j,j) = this%a(i,j)
  end do
end do

if (nx < 5) call sll_o_display(this%mat, 'f5.2')

call dpbtrf('U',(nx+1)*(ny+1),kd,this%mat,kd+1,error)
print*, 'info=', error

allocate(this%xd((nx+1)*(ny+1)))
allocate(this%yd((nx+1)*(ny+1)))

k = 0
do j = 1, ny+1
do i = 1, nx+1
  k = k+1
  this%xd(k) = x(i)
  this%yd(k) = y(j)
end do
end do

 
end subroutine initialize_poisson_2d_fem

!> Get the node index
integer function som(nx, i, j, k)

integer :: i, j, k, nx

if (k == 1) then
  som = i+(j-1)*(nx+1) 
else if (k == 2) then
  som = i+(j-1)*(nx+1)+1
else if (k == 3) then
  som = i+(j-1)*(nx+1)+nx+2
else if (k == 4) then
  som = i+(j-1)*(nx+1)+nx+1
end if 

end function som

!> Solve the poisson equation
subroutine solve_poisson_2d_fem( this, ex, ey, rho )
type( sll_t_fem_poisson_2d ) :: this !< Poisson solver object
sll_real64, dimension(:,:) :: ex   !< x electric field
sll_real64, dimension(:,:) :: ey   !< y electric field
sll_real64, dimension(:,:) :: rho  !< charge density
sll_real64, dimension((this%nx+1)*(this%ny+1)) :: b

sll_int32 :: nx
sll_int32 :: ny
sll_int32 :: i, j, k
sll_int32 :: error

nx = this%nx
ny = this%ny

!** Construction du second membre 
k = 0
do j=1,ny+1
do i=1,nx+1
  k = k+1
  b(k) = rho(i,j)
end do
end do

b = matmul(this%M,b)

do i = 1, this%nd
  k = this%id(i)
  b(k) = 1.0e20 * (this%xd(k)*this%xd(k)+this%yd(k)*this%yd(k))
end do

if (nx < 5) call sll_o_display(b, 'f5.2')

call dpbtrs('U',(nx+1)*(ny+1),nx+2,1,this%mat,nx+3,b,(nx+1)*(ny+1),error) 

rho = 0.0_8
k = 0
do j=1,ny+1
do i=1,nx+1
  k = k+1
  rho(i,j) = b(k) 
end do
end do

do j=1,ny+1
do i=1,nx
  ex(i,j) = - (rho(i+1,j)-rho(i,j)) / this%hx(i)
end do
end do

do j=1,ny
do i=1,nx+1
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
sll_int32 :: i, j

nx = size(x)-1
ny = size(y)-1
open(10, file="mesh.mtv")
write(10,*)"$DATA=CURVE3D"
write(10,*)"%xmin=", 1.1*minval(x), " xmax = ", 1.1*maxval(x)
write(10,*)"%ymin=", 1.1*minval(y), " ymax = ", 1.1*maxval(y)
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
do i=1,nx
  do j=1,ny
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
isom = 0
do j=1,ny+1
  do i=1,nx+1
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


end module sll_m_fem_2d
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
