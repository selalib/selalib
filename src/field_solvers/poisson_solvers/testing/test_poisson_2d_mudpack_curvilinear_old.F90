program test_poisson_2d_mudpack_curvilinear_old

#include "sll_memory.h"
#include "sll_working_precision.h"

use sll_m_boundary_condition_descriptors, only: &
  sll_p_dirichlet,                              &
  sll_p_hermite,                                &
  sll_p_periodic

use sll_m_cartesian_meshes, only: &
  sll_f_new_cartesian_mesh_2d, &
  sll_t_cartesian_mesh_2d

use sll_m_constants, only: &
  sll_p_pi

use sll_m_coordinate_transformation_2d_base, only: &
  sll_c_coordinate_transformation_2d_base

use sll_m_common_coordinate_transformations, only: &
    sll_f_polar_jac11, &
    sll_f_polar_jac12, &
    sll_f_polar_jac21, &
    sll_f_polar_jac22, &
    sll_f_polar_x1, &
    sll_f_polar_x2

use sll_m_coordinate_transformations_2d, only: &
  sll_f_new_coordinate_transformation_2d_analytic

use sll_m_poisson_2d_base, only: &
  sll_c_poisson_2d_base

use sll_m_poisson_2d_mudpack_curvilinear_old, only: &
    sll_f_new_poisson_2d_mudpack_curvilinear_old

use sll_m_cubic_spline_interpolator_2d, only: &
  sll_f_new_cubic_spline_interpolator_2d

use sll_m_interpolators_2d_base, only: &
  sll_c_interpolator_2d

implicit none

type(sll_t_cartesian_mesh_2d),                  pointer :: mesh
class(sll_c_coordinate_transformation_2d_base), pointer :: tau
class(sll_c_poisson_2d_base),                   pointer :: poisson
class(sll_c_interpolator_2d),                   pointer :: phi_interp2d

sll_real64 :: r_min
sll_real64 :: r_max
    
sll_int32  :: nc_x1
sll_int32  :: nc_x2
sll_real64 :: x1_min
sll_real64 :: x1_max     
sll_real64 :: x2_min
sll_real64 :: x2_max     

sll_real64, dimension(:,:), pointer :: b11
sll_real64, dimension(:,:), pointer :: b12
sll_real64, dimension(:,:), pointer :: b21
sll_real64, dimension(:,:), pointer :: b22
sll_real64, dimension(:,:), pointer :: c

sll_real64, allocatable :: rho(:,:)
sll_real64, allocatable :: phi(:,:)

sll_int32  :: ierr
sll_int32  :: i1, i2
sll_real64 :: x1, x2
sll_real64 :: delta_x1, delta_x2
sll_real64 :: e, r

nc_x1 = 32
r_min = 0.01_f64
r_max = 2.00_f64
nc_x2 = 32

x1_min = r_min
x1_max = r_max
x2_min = 0._f64
x2_max = 2._f64*sll_p_pi

mesh => sll_f_new_cartesian_mesh_2d( &
  nc_x1,                             &
  nc_x2,                             &
  eta1_min = x1_min,                 &
  eta1_max = x1_max,                 &
  eta2_min = x2_min,                 &
  eta2_max = x2_max)      

delta_x1 = mesh%delta_eta1
delta_x2 = mesh%delta_eta2

tau => sll_f_new_coordinate_transformation_2d_analytic( &
  "analytic_polar_transformation", &
  mesh, &
  sll_f_polar_x1, &
  sll_f_polar_x2, &
  sll_f_polar_jac11, &
  sll_f_polar_jac12, &
  sll_f_polar_jac21, &
  sll_f_polar_jac22, &
  params=(/0._f64,0._f64,0._f64,0._f64/))  

SLL_CLEAR_ALLOCATE(b11(1:nc_x1+1,1:nc_x2+1),ierr)
SLL_CLEAR_ALLOCATE(b12(1:nc_x1+1,1:nc_x2+1),ierr)
SLL_CLEAR_ALLOCATE(b21(1:nc_x1+1,1:nc_x2+1),ierr)
SLL_CLEAR_ALLOCATE(b22(1:nc_x1+1,1:nc_x2+1),ierr)
SLL_CLEAR_ALLOCATE(c(  1:nc_x1+1,1:nc_x2+1),ierr)

b11 = 1._f64
b22 = 1._f64
b12 = 0._f64
b21 = 0._f64
c   = 0._f64

phi_interp2d => sll_f_new_cubic_spline_interpolator_2d( &
  nc_x1+1,                                              &
  nc_x2+1,                                              &
  x1_min,                                               &
  x1_max,                                               &
  x2_min,                                               &
  x2_max,                                               &
  sll_p_periodic,                                       &
  sll_p_periodic)         

poisson => sll_f_new_poisson_2d_mudpack_curvilinear_old( &
 tau,                                                    &
 x1_min,                                                 &
 x1_max,                                                 &
 nc_x1,                                                  &
 x2_min,                                                 &
 x2_max,                                                 &
 nc_x2,                                                  &
 sll_p_dirichlet,                                        &
 sll_p_dirichlet,                                        &
 sll_p_periodic,                                         &
 sll_p_periodic,                                         &
 sll_p_hermite,                                          &
 sll_p_periodic,                                         &
 b11,                                                    &
 b12,                                                    &
 b21,                                                    &
 b22,                                                    &
 c)

SLL_CLEAR_ALLOCATE(rho(1:nc_x1+1,1:nc_x2+1),ierr)
SLL_CLEAR_ALLOCATE(phi(1:nc_x1+1,1:nc_x2+1),ierr)

do i2=1,nc_x2+1
  x2=x2_min+real(i2-1,f64)*delta_x2
  do i1=1,nc_x1+1
    x1=x1_min+real(i1-1,f64)*delta_x1
    r = sqrt(tau%x1(x1,x2)**2+tau%x2(x1,x2)**2)
    rho(i1,i2) = 4._f64 * sll_p_pi * f(r)
  end do
end do

call poisson%compute_phi_from_rho( phi, rho )

e = 0.0_f64
do i2=1,nc_x2+1
  x2=x2_min+real(i2-1,f64)*delta_x2
  do i1=1,nc_x1+1
    x1=x1_min+real(i1-1,f64)*delta_x1
    r = sqrt(tau%x1(x1,x2)**2+tau%x2(x1,x2)**2)
    e=e+abs(u(r)-phi(i1,i2))
    !write(10,*) tau%x1(x1,x2), tau%x2(x1,x2), phi(i1,i2), u(r)
  end do
  !write(10,*)
end do

write(*,"(' error on phi : ', g15.5)") e / (nc_x1*nc_x2)

if ( e / (nc_x1*nc_x2) < 0.1 ) then
  print*, 'PASSED'
else
  print*, 'FAILED'
end if


contains

!Charge density is a solid cylinder of radius 0.2
function f(r)

  sll_real64 :: f
  sll_real64 :: r

  if ( 0._f64 <= r .and. r <= 0.2_f64 ) then
    f = 1.0_f64
  else 
    f = 0.0_f64
  end if

end function f

!We have the equation :
! -4 pi f(r) = Laplacian(u(r))
!If
!  f(r) = rho                    for 0   <= r <= 0.2
!  f(r) = 0.0                    for 0.2 <  r <= 1
!Then
!  u(r) = -pi * f(r) * r^2 + a_1  for 0   <= r <= 0.2
!  u(r) = a_2 * ln(r)            for 0.2 <  r <= 1
!
!  a_1 =  0.04 * pi * f(r) * (-2*ln(0.2)+1)
!  a_2 = -0.08 * pi * f(r)

function u(r)

  sll_real64 :: u
  sll_real64 :: r
  sll_real64 :: pi
  sll_real64 :: a_0, a_1, a_2, a_3

  pi  =  4.0_f64 * atan(1._f64)
  a_0 =  0.0_f64
  a_1 =  0.04_f64 * pi * f(r) * (-2.0_f64*log(0.2_f64)+1.0_f64)
  a_2 = -0.08_f64 * pi !* f(r)
  a_3 =  0.0_f64

  if (0._f64 < r .and. r <= 0.2_f64) then
    u = -pi * f(r) * r*r + a_0*log(r) + a_1 
  else if ( 0.2_f64 < r .and. r <= 1._f64) then
    u = a_2 * log(r) + a_3
  else 
    u = 0._f64
  end if

end function u

end program test_poisson_2d_mudpack_curvilinear_old
