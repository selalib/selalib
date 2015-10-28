!> @internal [example]
program unit_test_meshes
#include "sll_working_precision.h"
#include "sll_memory.h"

use sll_m_constants     , only: sll_pi
use sll_m_array_plotting, only: &
  write_projection_2d, &
  SLL_X1X2, &
  SLL_X1X3, &
  SLL_X1X4, &
  SLL_X2X3, &
  SLL_X2X4, &
  SLL_X3X4

use sll_m_cartesian_meshes

implicit none

type(sll_cartesian_mesh_2d), pointer :: mx
type(sll_cartesian_mesh_2d), pointer :: mv
type(sll_cartesian_mesh_4d), pointer :: mxv

sll_real64, allocatable :: f(:,:,:,:)

sll_int32  :: i, j, k, l
sll_int32  :: n1, n2, n3, n4
sll_int32  :: error
sll_real64 :: x, y, vx, vy, v2, kx, ky, eps
sll_int32  :: iplot = 1

mx => new_cartesian_mesh_2d(100,100, 0.0_f64, 12.56_f64, 0.0_f64, 12.56_f64)
mv => new_cartesian_mesh_2d(64,64,-6.0_f64,6.0_f64,-6.0_f64,6.0_f64)

mxv => mx * mv
call sll_display(mxv)

n1 = mx%num_cells1+1
n2 = mx%num_cells2+1
n3 = mv%num_cells1+1
n4 = mv%num_cells2+1

SLL_CLEAR_ALLOCATE(f(1:n1,1:n2,1:n3,1:n4),error)

eps = 1.0_f64

kx = 2.0_f64 * sll_pi / (mx%num_cells1*mx%delta_eta1)
ky = 2.0_f64 * sll_pi / (mx%num_cells2*mx%delta_eta2)

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

call write_projection_2d(mxv, f, 'f_x1x2', SLL_X1X2, [32,32], iplot)
call write_projection_2d(mxv, f, 'f_x1x3', SLL_X1X3, [32,32], iplot)
call write_projection_2d(mxv, f, 'f_x1x4', SLL_X1X4, [32,32], iplot)
call write_projection_2d(mxv, f, 'f_x2x3', SLL_X2X3, [32,32], iplot)
call write_projection_2d(mxv, f, 'f_x2x4', SLL_X2X4, [32,32], iplot)
call write_projection_2d(mxv, f, 'f_x3x4', SLL_X3X4, [32,32], iplot)

call sll_delete(mxv)
call sll_delete(mx)
call sll_delete(mv)

print *, 'PASSED'


contains

  function landau_cos_prod(eps, kx, ky, x, y, v2)
    sll_real64 :: landau_cos_prod
    sll_real64, intent(in) :: x, kx
    sll_real64, intent(in) :: y, ky
    sll_real64, intent(in) :: eps, v2

    landau_cos_prod = (1._f64+eps*cos(kx*x)*cos(ky*y))/(2*sll_pi)*exp(-0.5_f64*v2)

  end function landau_cos_prod

end program unit_test_meshes
