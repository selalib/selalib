!> @internal [example]
program test_cartesian_meshes
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use sll_m_cartesian_meshes, only: &
    sll_f_new_cartesian_mesh_1d, &
    sll_f_new_cartesian_mesh_2d, &
    sll_f_new_cartesian_mesh_3d, &
    sll_f_new_cartesian_mesh_4d, &
    sll_t_cartesian_mesh_1d, &
    sll_t_cartesian_mesh_2d, &
    sll_t_cartesian_mesh_3d, &
    sll_t_cartesian_mesh_4d, &
    sll_o_delete, &
    sll_o_display, &
    operator(*)

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


type(sll_t_cartesian_mesh_2d), pointer :: m2d
type(sll_t_cartesian_mesh_3d), pointer :: m3d
type(sll_t_cartesian_mesh_4d), pointer :: m4d

type(sll_t_cartesian_mesh_1d), pointer :: m1d_x
type(sll_t_cartesian_mesh_1d), pointer :: m1d_y
type(sll_t_cartesian_mesh_2d), pointer :: m2d_xy

type(sll_t_cartesian_mesh_2d), pointer :: mx
type(sll_t_cartesian_mesh_2d), pointer :: mv
type(sll_t_cartesian_mesh_4d), pointer :: mxv

m1d_x => sll_f_new_cartesian_mesh_1d(64)
m1d_y => sll_f_new_cartesian_mesh_1d(32,eta_min=-1.0_f64,eta_max=+1.0_f64)

call sll_o_display(m1d_x)
call sll_o_display(m1d_y)

m2d_xy => m1d_x * m1d_y
call sll_o_display(m2d_xy)
call sll_o_delete(m2d_xy)

call sll_o_delete(m1d_x)
call sll_o_delete(m1d_y)

m2d => sll_f_new_cartesian_mesh_2d(100,100)
m3d => sll_f_new_cartesian_mesh_3d(100,100,100)
m4d => sll_f_new_cartesian_mesh_4d(32,32,32,32, eta1_min=-1.0_f64, eta1_max = 2.0_f64)

call sll_o_display(m2d)
call sll_o_display(m3d)
call sll_o_display(m4d)

call sll_o_delete(m2d)
call sll_o_delete(m3d)
call sll_o_delete(m4d)

mx => sll_f_new_cartesian_mesh_2d(100,100, 0.0_f64, 12.56_f64, 0.0_f64, 12.56_f64)
mv => sll_f_new_cartesian_mesh_2d(64,64,-6.0_f64,6.0_f64,-6.0_f64,6.0_f64)

mxv => mx * mv
call sll_o_display(mxv)
call sll_o_delete(mxv)

call sll_o_delete(mx)
call sll_o_delete(mv)

print *, 'PASSED'

end program test_cartesian_meshes
!> @internal [example]
