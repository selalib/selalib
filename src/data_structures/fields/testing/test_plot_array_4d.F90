!> @internal [example]
program unit_test_meshes
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_array_plotting, only: &
    sll_p_x1x2, &
    sll_p_x1x3, &
    sll_p_x1x4, &
    sll_p_x2x3, &
    sll_p_x2x4, &
    sll_p_x3x4, &
    sll_s_write_projection_2d

  use sll_m_cartesian_meshes, only: &
    sll_f_new_cartesian_mesh_2d, &
    sll_t_cartesian_mesh_2d, &
    sll_t_cartesian_mesh_4d, &
    sll_o_delete, &
    sll_o_display, &
    operator(*)

  use sll_m_constants, only: &
    sll_p_pi

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

type(sll_t_cartesian_mesh_2d), pointer :: mx
type(sll_t_cartesian_mesh_2d), pointer :: mv
type(sll_t_cartesian_mesh_4d), pointer :: mxv

sll_real64, allocatable :: f(:,:,:,:)

sll_int32  :: i, j, k, l
sll_int32  :: n1, n2, n3, n4
sll_int32  :: error
sll_real64 :: x, y, vx, vy, v2, kx, ky, eps
sll_int32  :: iplot = 1

mx => sll_f_new_cartesian_mesh_2d(100,100, 0.0_f64, 12.56_f64, 0.0_f64, 12.56_f64)
mv => sll_f_new_cartesian_mesh_2d(64,64,-6.0_f64,6.0_f64,-6.0_f64,6.0_f64)

mxv => mx * mv
call sll_o_display(mxv)

n1 = mx%num_cells1+1
n2 = mx%num_cells2+1
n3 = mv%num_cells1+1
n4 = mv%num_cells2+1

SLL_CLEAR_ALLOCATE(f(1:n1,1:n2,1:n3,1:n4),error)

eps = 1.0_f64

kx = 2.0_f64 * sll_p_pi / (mx%num_cells1*mx%delta_eta1)
ky = 2.0_f64 * sll_p_pi / (mx%num_cells2*mx%delta_eta2)

do l=1,n4
  do k=1,n3
    do j=1,n2
      do i=1,n1
      
        x  = mxv%eta1_min+(i-1)*mxv%delta_eta1
        y  = mxv%eta2_min+(j-1)*mxv%delta_eta2
        vx = mxv%eta3_min+(k-1)*mxv%delta_eta3
        vy = mxv%eta4_min+(l-1)*mxv%delta_eta4
      
        v2 = vx*vx+vy*vy
      
        f(i,j,k,l) = landau_cos_prod(eps, kx, ky, x, y, v2)
      
      end do
    end do
  end do
end do

call sll_s_write_projection_2d(mxv, f, 'f_x1x2', sll_p_x1x2, [32,32], iplot)
call sll_s_write_projection_2d(mxv, f, 'f_x1x3', sll_p_x1x3, [32,32], iplot)
call sll_s_write_projection_2d(mxv, f, 'f_x1x4', sll_p_x1x4, [32,32], iplot)
call sll_s_write_projection_2d(mxv, f, 'f_x2x3', sll_p_x2x3, [32,32], iplot)
call sll_s_write_projection_2d(mxv, f, 'f_x2x4', sll_p_x2x4, [32,32], iplot)
call sll_s_write_projection_2d(mxv, f, 'f_x3x4', sll_p_x3x4, [32,32], iplot)

call sll_o_delete(mxv)
call sll_o_delete(mx)
call sll_o_delete(mv)

print *, 'PASSED'


contains

  function landau_cos_prod(eps, kx, ky, x, y, v2)
    sll_real64 :: landau_cos_prod
    sll_real64, intent(in) :: x, kx
    sll_real64, intent(in) :: y, ky
    sll_real64, intent(in) :: eps, v2

    landau_cos_prod = (1._f64+eps*cos(kx*x)*cos(ky*y))/(2*sll_p_pi)*exp(-0.5_f64*v2)

  end function landau_cos_prod

end program unit_test_meshes
