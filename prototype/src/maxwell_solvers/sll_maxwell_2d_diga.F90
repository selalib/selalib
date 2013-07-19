!>Solve Poisson equation on cartesian domain with finit elements.
!> * Compact boundary conditions.
!> * Linear system solve with lapack (Choleski)
module sll_maxwell_2d_diga
#include "sll_maxwell_solvers_macros.h"
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_constants.h"

use sll_logical_meshes
use sll_module_coordinate_transformations_2d
use sll_common_coordinate_transformations

implicit none

type :: maxwell_2d_diga
   type(sll_logical_mesh_2d), pointer :: mesh
   class(sll_coordinate_transformation_2d_analytic), pointer :: tau
   sll_int32 :: polarization
end type maxwell_2d_diga

interface initialize
   module procedure initialize_maxwell_2d_diga
end interface initialize
interface solve
   module procedure solve_maxwell_2d_diga
end interface solve

sll_int32, private :: nx, ny
sll_int32, private :: i, j, k
sll_int32, private :: error

private :: write_mtv_file

contains

!> Initialize Poisson solver object using finite elements method.
!> Indices are shifted from [1:n+1] to [0:n] only inside this 
!> subroutine
subroutine initialize_maxwell_2d_diga( this, mesh, tau, polarization)
type( maxwell_2d_diga ) :: this !< solver data object
type(sll_logical_mesh_2d), pointer :: mesh
class(sll_coordinate_transformation_2d_analytic), pointer :: tau
sll_int32 :: polarization

this%mesh => mesh
this%tau  => tau
this%polarization = polarization

call write_mtv_file( mesh, tau)

end subroutine initialize_maxwell_2d_diga

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

!> Solve the maxwell equation
subroutine solve_maxwell_2d_diga( this, ex, ey, bz, dt, jx, jy, rho )

type( maxwell_2d_diga )    :: this !< Maxwell solver object

sll_real64, dimension(:,:) :: ex   !< x electric field
sll_real64, dimension(:,:) :: ey   !< y electric field
sll_real64, dimension(:,:) :: bz   !< z magnetic field
sll_real64, intent(in)     :: dt

sll_real64, dimension(:,:), optional :: jx   !< x current field
sll_real64, dimension(:,:), optional :: jy   !< y current field
sll_real64, dimension(:,:), optional :: rho  !< charge density

end subroutine solve_maxwell_2d_diga

subroutine write_mtv_file( mesh, tau )
sll_int32 :: nc_x
sll_int32 :: nc_y
type(sll_logical_mesh_2d) :: mesh
type(sll_coordinate_transformation_2d_analytic) :: tau
integer :: iel, isom
real(8) :: x1, y1

nc_x = mesh%num_cells1
nc_y = mesh%num_cells2

open(10, file="mesh.mtv")
write(10,*)"$DATA=CURVE3D"
write(10,*)"%equalscale=T"
write(10,*)"%toplabel='Elements number ' "
   
do i=1,nc_x-1
   do j=1,nc_y-1
      write(10,*) tau%x1_at_node(i  ,j), tau%x2_at_node(i,j  ), 0.
      write(10,*) tau%x1_at_node(i+1,j), tau%x2_at_node(i,j  ), 0.
      write(10,*) tau%x1_at_node(i+1,j), tau%x2_at_node(i,j+1), 0.
      write(10,*) tau%x1_at_node(i  ,j), tau%x2_at_node(i,j+1), 0.
      write(10,*) tau%x1_at_node(i  ,j), tau%x2_at_node(i,j  ), 0.
      write(10,*)
   end do
end do

!Numeros des elements
iel = 0
do i=1,nc_x-1
   do j=1,nc_y-1
      iel = iel+1
      x1 = 0.5*(tau%x1_at_node(i,j)+tau%x1_at_node(i+1,j))
      y1 = 0.5*(tau%x2_at_node(i,j)+tau%x2_at_node(i,j+1))
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
do i=1,nc_x-1
   do j=1,nc_y-1
      isom = isom+1
      write(10,"(a)"   ,  advance="no")"@text x1="
      write(10,"(g15.3)", advance="no") tau%x1_at_node(i,j)
      write(10,"(a)"   ,  advance="no")" y1="
      write(10,"(g15.3)", advance="no") tau%x2_at_node(i,j)
      write(10,"(a)"   ,  advance="no")" z1=0. lc=5 ll='"
      write(10,"(i4)"  ,  advance="no") isom
      write(10,"(a)")"'"
   end do
end do
   
write(10,*)"$END"
close(10)

end subroutine write_mtv_file

end module sll_maxwell_2d_diga
