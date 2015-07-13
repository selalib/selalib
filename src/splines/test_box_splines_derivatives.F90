!*****************************************************************************
!> @brief
!> sll_box_splines unit test
!> @author
!> Laura S. Mendoza
!*****************************************************************************

program box_spline_tester

#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_assert.h"

use sll_utilities
use sll_constants
use sll_hex_meshes
use sll_box_splines
implicit none

type(sll_hex_mesh_2d),   pointer  :: mesh
type(sll_box_spline_2d), pointer  :: spline

sll_int32    :: ierr
sll_int32    :: num_cells
sll_int32    :: deg
sll_int32    :: rule
sll_int32    :: i
sll_real64   :: x1
sll_real64   :: x2
sll_real64, dimension(:), allocatable :: f
sll_real64, dimension(:), allocatable :: dxf
sll_real64, dimension(:), allocatable :: dyf
sll_real64, dimension(:), allocatable :: splines_on_support


! Mesh initialization
num_cells = 20
mesh => new_hex_mesh_2d(num_cells, 0._f64, 0._f64, radius = 2._f64)
call sll_display(mesh)

spline => new_box_spline_2d(mesh, SLL_DIRICHLET)

! Allocations for boxsplines and derivatives :
SLL_ALLOCATE(f(mesh%num_pts_tot),ierr)
SLL_ALLOCATE(dxf(mesh%num_pts_tot),ierr)
SLL_ALLOCATE(dyf(mesh%num_pts_tot),ierr)


! do i=1, mesh%num_pts_tot

!    x1 = mesh%global_to_x1(i)
!    x2 = mesh%global_to_x2(i)

!    ! Computing boxsplines of degree 2:
!    f(i) = compute_box_spline(spline, x1, x2, 2)
!    ! And derivatives :
!    dxf(i) = boxspline_val_der(x1, x2, 1, 1, 0)!boxspline_x1_derivative(x1, x2, 1)
!    dyf(i) = boxspline_val_der(x1, x2, 1, 0, 1)!boxspline_x2_derivative(x1, x2, 1)

! end do

!Wrtting on docs:
call write_field_hex_mesh_xmf(mesh,   f, "boxspline2")
call write_field_hex_mesh_xmf(mesh, dxf, "der1_boxspline2")
call write_field_hex_mesh_xmf(mesh, dyf, "der2_boxspline2")


SLL_DEALLOCATE_ARRAY(f, ierr)
SLL_DEALLOCATE_ARRAY(dxf, ierr)
SLL_DEALLOCATE_ARRAY(dyf, ierr)


!Writing file for CAID/DJANGO:
deg = 2
rule = 1
call write_basis_values(deg, rule)
print *, ""
print *, "Done writing CAID file : basis_value.txt"


! Computing non null splines on one cell:
SLL_ALLOCATE(splines_on_support(deg*deg*3),ierr)
splines_on_support = non_zeros_splines(mesh, 2, deg)

print *, "Non null splines of degree 1 at cell#2 (expected 1,2,3) =", &
     splines_on_support(1), &
     splines_on_support(2), &
     splines_on_support(3)

SLL_DEALLOCATE_ARRAY(splines_on_support, ierr)

end program box_spline_tester
