program test_poisson_2d_mudpack_curvilinear

#include "sll_memory.h"
#include "sll_working_precision.h"

use sll_m_boundary_condition_descriptors, only: &
  sll_p_periodic

use sll_m_cartesian_meshes, only: &
  sll_f_new_cartesian_mesh_2d, &
  sll_t_cartesian_mesh_2d

use sll_m_common_coordinate_transformations, only: &
  sll_f_identity_jac11, &
  sll_f_identity_jac12, &
  sll_f_identity_jac21, &
  sll_f_identity_jac22, &
  sll_f_identity_x1, &
  sll_f_identity_x2, &
  sll_f_sinprod_jac11, &
  sll_f_sinprod_jac12, &
  sll_f_sinprod_jac21, &
  sll_f_sinprod_jac22, &
  sll_f_sinprod_x1, &
  sll_f_sinprod_x2

use sll_m_constants, only: &
  sll_p_pi

use sll_m_mudpack_curvilinear, only: &
  sll_p_non_separable_with_cross_terms, &
  sll_p_non_separable_without_cross_terms, &
  sll_p_separable

use sll_m_poisson_2d_base, only: &
  sll_c_poisson_2d_base

use sll_m_poisson_2d_mudpack_curvilinear, only: &
  sll_f_new_poisson_2d_mudpack_curvilinear

implicit none

type(sll_t_cartesian_mesh_2d), pointer  :: mesh
class(sll_c_poisson_2d_base),  pointer  :: poisson
sll_real64                              :: x1_min
sll_real64                              :: x2_min

sll_real64                              :: kmode_x1
sll_real64                              :: kmode_x2

sll_int32                               :: nc_x1
sll_int32                               :: nc_x2
sll_real64                              :: x1_max
sll_real64                              :: x2_max     

sll_real64, dimension(:,:), allocatable :: cxx
sll_real64, dimension(:,:), allocatable :: cxy
sll_real64, dimension(:,:), allocatable :: cyy
sll_real64, dimension(:,:), allocatable :: cx
sll_real64, dimension(:,:), allocatable :: cy
sll_real64, dimension(:,:), allocatable :: ce
sll_int32                               :: ierr

sll_real64                              :: delta_x1
sll_real64                              :: delta_x2
sll_real64                              :: x1
sll_real64                              :: x2
sll_real64                              :: e
sll_int32                               :: i1 
sll_int32                               :: i2
sll_real64, dimension(:,:), pointer     :: rho
sll_real64, dimension(:,:), pointer     :: phi


nc_x1  = 32
x1_min = 0.0_f64
nc_x2  = 32
x2_min = 0.0_f64

kmode_x1 = 1.0_f64
kmode_x2 = 1.0_f64

x1_max = 2._f64*sll_p_pi/kmode_x1
x2_max = 2._f64*sll_p_pi/kmode_x2
 
mesh => sll_f_new_cartesian_mesh_2d( &
  nc_x1,                                &
  nc_x2,                                &
  eta1_min = x1_min,                    &
  eta1_max = x1_max,                    &
  eta2_min = x2_min,                    &
  eta2_max = x2_max)      

delta_x1 = mesh%delta_eta1
delta_x2 = mesh%delta_eta2
      
SLL_ALLOCATE(cxx(nc_x1+1,nc_x2+1),ierr)
SLL_ALLOCATE(cxy(nc_x1+1,nc_x2+1),ierr)
SLL_ALLOCATE(cyy(nc_x1+1,nc_x2+1),ierr)
SLL_ALLOCATE(cx( nc_x1+1,nc_x2+1),ierr)
SLL_ALLOCATE(cy( nc_x1+1,nc_x2+1),ierr)
SLL_ALLOCATE(ce( nc_x1+1,nc_x2+1),ierr)

cxx = 1._f64
cxy = 0._f64
cyy = 1._f64
cx  = 0._f64
cy  = 0._f64
ce  = 0._f64
 
poisson => sll_f_new_poisson_2d_mudpack_curvilinear( &
  x1_min,                                            &
  x1_max,                                            &
  nc_x1,                                             &
  x2_min,                                            &
  x2_max,                                            &
  nc_x2,                                             &
  sll_p_periodic,                                    & 
  sll_p_periodic,                                    & 
  sll_p_periodic,                                    & 
  sll_p_periodic,                                    &
  sll_p_non_separable_without_cross_terms,           &
  cxx_2d = cxx,                                      &
  cxy_2d = cxy,                                      &
  cyy_2d = cyy,                                      &
  cx_2d  = cx,                                       &
  cy_2d  = cy,                                       &
  ce_2d  = ce)
              
SLL_CLEAR_ALLOCATE(rho(1:nc_x1+1,1:nc_x2+1),ierr)
SLL_CLEAR_ALLOCATE(phi(1:nc_x1+1,1:nc_x2+1),ierr)

do i2=1,nc_x2+1
  x2=x2_min+real(i2-1,f64)*delta_x2
  do i1=1,nc_x1+1
    x1=x1_min+real(i1-1,f64)*delta_x1
    rho(i1,i2) = 2.0_f64*sin(x1)*sin(x2)
  end do
end do

call poisson%compute_phi_from_rho( phi, rho )            

e = 0.0_f64
do i2=1,nc_x2+1
  x2=x2_min+real(i2-1,f64)*delta_x2
  do i1=1,nc_x1+1
    x1=x1_min+real(i1-1,f64)*delta_x1
    e=e+abs(sin(x1)*sin(x2)-phi(i1,i2))
    write(10,*) x1, x2, phi(i1,i2), rho(i1,i2)
  end do
  write(10,*)
end do

write(*,"(' error on phi : ', g15.5)") e / (nc_x1*nc_x2)

if ( e / (nc_x1*nc_X2) < 0.0001) then
  print*, 'PASSED'
else
  print*, 'FAILED'
end if


end program test_poisson_2d_mudpack_curvilinear
