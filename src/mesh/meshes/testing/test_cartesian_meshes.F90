!> @internal [example]
program test_cartesian_meshes
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use sll_m_cartesian_meshes, only: &
    new_cartesian_mesh_1d, &
    new_cartesian_mesh_2d, &
    new_cartesian_mesh_3d, &
    new_cartesian_mesh_4d, &
    sll_cartesian_mesh_1d, &
    sll_cartesian_mesh_2d, &
    sll_cartesian_mesh_3d, &
    sll_cartesian_mesh_4d, &
    sll_delete, &
    sll_display, &
    operator(*)

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


type(sll_cartesian_mesh_2d), pointer :: m2d
type(sll_cartesian_mesh_3d), pointer :: m3d
type(sll_cartesian_mesh_4d), pointer :: m4d

type(sll_cartesian_mesh_1d), pointer :: m1d_x
type(sll_cartesian_mesh_1d), pointer :: m1d_y
type(sll_cartesian_mesh_2d), pointer :: m2d_xy

type(sll_cartesian_mesh_2d), pointer :: mx
type(sll_cartesian_mesh_2d), pointer :: mv
type(sll_cartesian_mesh_4d), pointer :: mxv

m1d_x => new_cartesian_mesh_1d(64)
m1d_y => new_cartesian_mesh_1d(32,eta_min=-1.0_f64,eta_max=+1.0_f64)

call sll_display(m1d_x)
call sll_display(m1d_y)

m2d_xy => m1d_x * m1d_y
call sll_display(m2d_xy)
call sll_delete(m2d_xy)

call sll_delete(m1d_x)
call sll_delete(m1d_y)

m2d => new_cartesian_mesh_2d(100,100)
m3d => new_cartesian_mesh_3d(100,100,100)
m4d => new_cartesian_mesh_4d(32,32,32,32, eta1_min=-1.0_f64, eta1_max = 2.0_f64)

call sll_display(m2d)
call sll_display(m3d)
call sll_display(m4d)

call sll_delete(m2d)
call sll_delete(m3d)
call sll_delete(m4d)

mx => new_cartesian_mesh_2d(100,100, 0.0_f64, 12.56_f64, 0.0_f64, 12.56_f64)
mv => new_cartesian_mesh_2d(64,64,-6.0_f64,6.0_f64,-6.0_f64,6.0_f64)

mxv => mx * mv
call sll_display(mxv)
call sll_delete(mxv)

call sll_delete(mx)
call sll_delete(mv)

print *, 'PASSED'

end program test_cartesian_meshes
!> @internal [example]
