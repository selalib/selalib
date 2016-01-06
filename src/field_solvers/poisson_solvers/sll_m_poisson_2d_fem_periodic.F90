!> @ingroup poisson_solvers
!> @brief 
!> Poisson solver using finite elements
!> @details
!> Solve Poisson equation on cartesian domain with finite elements.
!> * Peridoci boundary conditions.
!> * Linear system solve with lapack (Choleski)
!> This solver is not fully tested, please use it carefully.
module sll_m_poisson_2d_fem_periodic
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_poisson_solvers_macros.h"

implicit none

public :: &
  sll_o_create, &
  sll_t_poisson_2d_fem_periodic, &
  sll_o_solve

private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!> Structure to solve Poisson equation on 2d irregular cartesian mesh
!> with finite element numerical method
type :: sll_t_poisson_2d_fem_periodic
   sll_int32                           :: nx    !< cells number along x
   sll_int32                           :: ny    !< cells number along y
   sll_real64, dimension(:,:), pointer :: A     !< Mass matrix
   sll_real64, dimension(:,:), pointer :: M     !< Stiffness matrix
   sll_real64, dimension(:)  , pointer :: hx    !< step size x
   sll_real64, dimension(:)  , pointer :: hy    !< step size y
   sll_int32,  dimension(:)  , pointer :: ipiv  !< Lapack array for pivoting
end type sll_t_poisson_2d_fem_periodic

!> Initialize the solver
interface sll_o_create
   module procedure initialize_poisson_2d_periodic_fem
end interface sll_o_create

!> Compute the electric potential
interface sll_o_solve
   module procedure solve_poisson_2d_periodic_fem
end interface sll_o_solve

contains

!> Initialize Poisson solver object using finite elements method.
!> Indices are shifted from \f$ [1:n+1] \f$ to \f$ [0:n] \f$ only 
!> inside this subroutine.
subroutine initialize_poisson_2d_periodic_fem( this, x, y ,nn_x, nn_y)
type( sll_t_poisson_2d_fem_periodic ) :: this !< Solver data structure
sll_int32,  intent(in)                :: nn_x !< number of cells along x
sll_int32,  intent(in)                :: nn_y !< number of cells along y
sll_real64, dimension(nn_x)           :: x    !< x nodes coordinates
sll_real64, dimension(nn_y)           :: y    !< y nodes coordinates

sll_real64, dimension(4,4) :: Axelem
sll_real64, dimension(4,4) :: Ayelem
sll_real64, dimension(4,4) :: Melem
sll_int32,  dimension(4)   :: isom
sll_int32 :: nxy
sll_int32 :: i
sll_int32 :: j
sll_int32 :: error

this%nx = nn_x-1
this%ny = nn_y-1
nxy = this%nx*this%ny

SLL_ALLOCATE(this%hx(1:this%nx),error)
SLL_ALLOCATE(this%hy(1:this%ny),error)
SLL_ALLOCATE(this%A(nxy,nxy), error)
SLL_ALLOCATE(this%M(nxy,nxy), error)

do i=1,this%nx
   this%hx(i) = x(i+1)-x(i)
end do

do j=1,this%ny
   this%hy(j) = y(j+1)-y(j)
end do

!** Construction des matrices elementaires
Axelem(1,:) = (/  2.0_f64, -2.0_f64, -1.0_f64,  1.0_f64 /)
Axelem(2,:) = (/ -2.0_f64,  2.0_f64,  1.0_f64, -1.0_f64 /)
Axelem(3,:) = (/ -1.0_f64,  1.0_f64,  2.0_f64, -2.0_f64 /)
Axelem(4,:) = (/  1.0_f64, -1.0_f64, -2.0_f64,  2.0_f64 /)

Axelem = 1._f64/6._f64 * Axelem

Ayelem(1,:) = (/  2.0_f64,  1.0_f64, -1.0_f64, -2.0_f64 /)
Ayelem(2,:) = (/  1.0_f64,  2.0_f64, -2.0_f64, -1.0_f64 /)
Ayelem(3,:) = (/ -1.0_f64, -2.0_f64,  2.0_f64,  1.0_f64 /)
Ayelem(4,:) = (/ -2.0_f64, -1.0_f64,  1.0_f64,  2.0_f64 /)

Ayelem = 1._f64/6._f64 * Ayelem

Melem(1,:) = (/ 4.0_f64, 2.0_f64, 1.0_f64, 2.0_f64/)
Melem(2,:) = (/ 2.0_f64, 4.0_f64, 2.0_f64, 1.0_f64/)
Melem(3,:) = (/ 1.0_f64, 2.0_f64, 4.0_f64, 2.0_f64/)
Melem(4,:) = (/ 2.0_f64, 1.0_f64, 2.0_f64, 4.0_f64/)

Melem = 1._f64/36._f64 * Melem

this%A = 0._f64
this%M = 0._f64

!***  Interior mesh ***
do j=1,this%ny-1
  do i=1,this%nx-1
    isom(1) =  som(this%nx,i,j,1)
    isom(2) =  som(this%nx,i,j,2)
    isom(3) =  som(this%nx,i,j,3)
    isom(4) =  som(this%nx,i,j,4)
    call build_matrices( this, Axelem, Ayelem, Melem, isom, i, j )
  end do
end do

do i=2,this%nx-1

  j = 1
  isom(1)=som(this%nx,i,this%ny,1)
  isom(2)=som(this%nx,i,this%ny,2)
  isom(3)=som(this%nx,i,j,2)
  isom(4)=som(this%nx,i,j,1)
  call build_matrices( this, Axelem, Ayelem, Melem, isom, i, j )

end do

do j=2,this%ny-1

  i = 1
  isom(1)=som(this%nx,this%nx,j,1)
  isom(2)=som(this%nx,i,j,1)
  isom(3)=som(this%nx,i,j,4)
  isom(4)=som(this%nx,this%nx,j,4)
  call build_matrices( this, Axelem, Ayelem, Melem, isom, i, j )

end do
   
!Corners
i=1; j=1  !SW
isom(1) = som(this%nx,this%nx-1,this%ny-1,3)
isom(2) = som(this%nx,i   ,this%ny-1,4)
isom(3) = som(this%nx,i   ,j   ,1)
isom(4) = som(this%nx,this%nx-1,j   ,2)
call build_matrices( this, Axelem, Ayelem, Melem, isom, i, j )

SLL_ALLOCATE(this%ipiv(nxy),error)

this%A(1,1) = 1.0d7
call dgetrf(nxy,nxy,this%A,nxy,this%ipiv,error)

end subroutine initialize_poisson_2d_periodic_fem

!> Get the node number
integer function som(nx, i, j, k)
integer :: nx, i, j, k, l

l = i+(j-1)*nx
if (k == 1) then
  som = l
else if (k == 2) then
  som = l+1
else if (k == 3) then
  som = l+nx+1
else if (k == 4) then
  som = l+nx
end if 

end function som

!> Build matrices and factorize
subroutine build_matrices( this, Axelem, Ayelem, Melem, isom, i, j )
type( sll_t_poisson_2d_fem_periodic ) :: this     !< Poisson solver object
sll_real64, dimension(:,:)            :: Axelem   !< x electric field
sll_real64, dimension(:,:)            :: Ayelem   !< y electric field
sll_real64, dimension(:,:)            :: Melem    !< charge density
sll_int32,  dimension(:)              :: isom     !< node indices
sll_int32, intent(in)                 :: i        !< int(x) position on mesh
sll_int32, intent(in)                 :: j        !< int(y) position on mesh

sll_int32 :: ii
sll_int32 :: jj

do ii=1,4
  do jj=1,4
    this%A(isom(ii),isom(jj)) = this%A(isom(ii),isom(jj)) &
            & + Axelem(ii,jj) * this%hy(j) / this%hx(i)   &
            & + Ayelem(ii,jj) * this%hx(i) / this%hy(j)
    this%M(isom(ii),isom(jj)) = this%M(isom(ii),isom(jj)) &
            & + Melem(ii,jj)  * this%hx(i) * this%hy(j)
  end do
end do

end subroutine build_matrices

!> Solve the poisson equation
subroutine solve_poisson_2d_periodic_fem( this, ex, ey, rho )
type( sll_t_poisson_2d_fem_periodic )  :: this !< Poisson solver object
sll_real64, dimension(:,:)             :: ex   !< x electric field
sll_real64, dimension(:,:)             :: ey   !< y electric field
sll_real64, dimension(:,:)             :: rho  !< charge density
sll_real64, dimension(this%nx*this%ny) :: b
sll_real64                             :: bmoy

sll_int32 :: nxy
sll_int32 :: i, j, k
sll_int32 :: error

nxy = this%nx * this%ny
!** Construction du second membre (rho a support compact --> projete)
k=0
do i=1,this%nx
   do j=1,this%ny
      k=k+1
      b(k) = rho(i,j)
   end do
end do

b = matmul(this%M,b)
b(1) = 1.0_f64

call dgetrs('N',nxy,1,this%A,nxy,this%ipiv,b,nxy,error)

bmoy = sum(b) / real(nxy,f64)

k=0
do i=1,this%nx
   do j=1,this%ny
      k=k+1
      rho(i,j) = b(k) - bmoy
   end do
end do

do j=1,this%ny-1
do i=1,this%nx-2
   ex(i,j) = - (rho(i+1,j)-rho(i,j)) / this%hx(i)
end do
end do

do j=1,this%ny-2
do i=1,this%nx-1
   ey(i,j) = - (rho(i,j+1)-rho(i,j)) / this%hy(j)
end do
end do

end subroutine solve_poisson_2d_periodic_fem

end module sll_m_poisson_2d_fem_periodic
