program test_maxwell_2d_diga_wave
!--------------------------------------------------------------------------
!  test 2D Maxwell solver based on discontinuous galerkine on a mapped mesh
!--------------------------------------------------------------------------
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_constants.h"
#include "sll_maxwell_solvers_macros.h"
#include "sll_file_io.h"

use sll_logical_meshes
use sll_module_coordinate_transformations_2d
use sll_common_coordinate_transformations
use sll_dg_fields
use sll_maxwell_2d_diga
use sll_maxwell_solvers_base

implicit none

!=====================================!
! Simulation parameters               !
!=====================================!
sll_int32, parameter :: nc_eta1 = 10  !
sll_int32, parameter :: nc_eta2 = 10  !
sll_int32, parameter :: degree  = 4   !
!=====================================!

sll_int32  :: nstep
sll_real64 :: eta1_max, eta1_min
sll_real64 :: eta2_max, eta2_min
sll_real64 :: delta_eta1, delta_eta2

type(sll_logical_mesh_2d), pointer :: mesh
class(sll_coordinate_transformation_2d_analytic), pointer :: tau

type(maxwell_2d_diga)   :: maxwell_TE

type(dg_field), pointer :: ex, ex0, dx, sx
type(dg_field), pointer :: ey, ey0, dy, sy
type(dg_field), pointer :: bz, bz0, dz, sz

sll_real64  :: time
sll_int32   :: istep
sll_real64  :: dt
sll_real64  :: cfl = 0.5
character(len=4) :: cstep
!init functions
sll_real64, external :: gaussian

mesh => new_logical_mesh_2d(nc_eta1, nc_eta2, &
                            eta1_min=-5._f64, eta1_max=5._f64, &
                            eta2_min=-5._f64, eta2_max=5._f64)

write(*,"(3f8.3,i4)") mesh%eta1_min,mesh%eta1_max,mesh%delta_eta1,mesh%num_cells1
write(*,"(3f8.3,i4)") mesh%eta2_min,mesh%eta2_max,mesh%delta_eta2,mesh%num_cells2

eta1_min = mesh%eta1_min
eta1_max = mesh%eta1_max
eta2_min = mesh%eta2_min
eta2_max = mesh%eta2_max

delta_eta1 = mesh%delta_eta1
delta_eta2 = mesh%delta_eta2

! "Identity transformation";
tau => new_coordinate_transformation_2d_analytic( &
       "identity_transformation",                 &
       mesh,                                      &
       identity_x1,                               &
       identity_x2,                               &
       identity_jac11,                            &
       identity_jac12,                            &
       identity_jac21,                            &
       identity_jac22,                            &
       SLL_NULL_REAL64 )

! "Affine transformation";
!
! x1 = (B1-A1)*(cos(alpha)*eta1-sin(alpha)*eta2) + A1
! x2 = (B2-A2)*(sin(alpha)*eta1+cos(alpha)*eta2) + A2
!
!tau => new_coordinate_transformation_2d_analytic( &
!       "affine_transformation",                   &
!       mesh,                                      &
!       affine_x1,                                 &
!       affine_x2,                                 &
!       affine_jac11,                              &
!       affine_jac12,                              &
!       affine_jac21,                              &
!       affine_jac22,                              &
!       (/0.0_f64,2.0_f64,0.0_f64,1.0_f64,0.0_f64/) )

! "Homography transformation"
!
!        x1 = (a*eta1*b*eta2+c)/(g*eta1+h*eta2+1) 
!        x2 = (d*eta1*e*eta2+f)/(g*eta1+h*eta2+1) 
!
!  a = fixed scale factor in x1 direction with scale x2 unchanged.
!  b = scale factor in x1 direction proportional to x2 distance from origin.
!  c = origin translation in x1 direction.
!  d = scale factor in x2 direction proportional to x1 distance from origin.
!  e = fixed scale factor in x2 direction with scale x1 unchanged.
!  f = origin translation in x2 direction.
!  g = proportional scale factors x1 and x2 in function of x1.
!  h = proportional scale factors x1 and x2 in function of x2.
!   
!tau => new_coordinate_transformation_2d_analytic( &
!       "homography_transformation",               &
!       mesh,                                      &
!       homography_x1,                             &
!       homography_x2,                             &
!       homography_jac11,                          &
!       homography_jac12,                          &
!       homography_jac21,                          &
!       homography_jac22,                          &
!       [1.0_f64,1.0_f64,1.0_f64, &
!        0.0_f64,1.0_f64,0.0_f64, &
!        0.0_f64,0.0_f64] )
!
! "Colella transformation";
! sinusoidal product (see P. Colella et al. JCP 230 (2011) formula 
! (102) p 2968):
!
! x1 = eta1 + 0.1 * sin(2*pi*eta1) * sin(2*pi*eta2)
! x2 = eta2 + 0.1 * sin(2*pi*eta1) * sin(2*pi*eta2)

!tau => new_coordinate_transformation_2d_analytic( &
!       "collela_transformation",                  &
!       mesh,                                      &
!       sinprod_x1,                                &
!       sinprod_x2,                                &
!       sinprod_jac11,                             &
!       sinprod_jac12,                             &
!       sinprod_jac21,                             &
!       sinprod_jac22,                             &  
!       (/0.1_f64,0.1_f64,1.0_f64,1.0_f64/) )

!tau => new_coordinate_transformation_2d_analytic( &
!       "rubber_sheeting_transformation",          &
!       mesh,                                      &
!       rubber_sheeting_x1,                        &
!       rubber_sheeting_x2,                        &
!       rubber_sheeting_jac11,                     &
!       rubber_sheeting_jac12,                     &
!       rubber_sheeting_jac21,                     &
!       rubber_sheeting_jac22,                     &
!       [-1.0_f64,2.0_f64,0.0_f64,0.0_f64, &
!         0.0_f64,0.0_f64,1.0_f64,0.0_f64] )

call tau%write_to_file(SLL_IO_MTV)

time = 0.0_f64

ex  => new_dg_field(degree,tau) 
ey  => new_dg_field(degree,tau) 
bz  => new_dg_field(degree,tau,gaussian) 

ex0 => new_dg_field(degree,tau) 
ey0 => new_dg_field(degree,tau) 
bz0 => new_dg_field(degree,tau) 

dx  => new_dg_field(degree,tau) 
dy  => new_dg_field(degree,tau) 
dz  => new_dg_field(degree,tau) 

sx  => new_dg_field(degree,tau) 
sy  => new_dg_field(degree,tau) 
sz  => new_dg_field(degree,tau) 

dt = cfl/sqrt(1./(delta_eta1/(degree+1))**2+1./(delta_eta2/(degree+1))**2)
nstep = 200

call initialize(maxwell_TE,        &
                tau,               &
                degree,            &
                TE_POLARIZATION,   &
                SLL_SILVER_MULLER, &  !South
                SLL_SILVER_MULLER, &  !East
                SLL_SILVER_MULLER, &  !North
                SLL_SILVER_MULLER, &  !West
                SLL_UNCENTERED )      !Flux


do istep = 1, nstep !*** Loop over time

   write(*,"(' istep = ',I6,' time = ',g12.3,' sec')") istep, time

   call rksetup()

   call solve(maxwell_TE, ex, ey, bz, dx, dy, dz)

   call accumulate(1._f64/6._f64)
   call rkstage(0.5_f64)

   call solve(maxwell_TE, ex, ey, bz, dx, dy, dz)
   call accumulate(1._f64/3._f64)
   call rkstage(0.5_f64)

   call solve(maxwell_TE, ex, ey, bz, dx, dy, dz)
   call accumulate(1._f64/3._f64)
   call rkstage(1.0_f64)

   call solve(maxwell_TE, ex, ey, bz, dx, dy, dz)
   call accumulate(1._f64/6._f64)

   call rkstep()

   time = time + dt
   
   call int2string(istep, cstep)
   call bz%write_to_file('bz')
   !call bz%write_to_file('bz'//cstep, SLL_IO_XDMF, time)

end do ! next time step

print"(a)", 'PASSED'

contains

subroutine rksetup()

   sx%array = 0.0_f64
   sy%array = 0.0_f64
   sz%array = 0.0_f64

   ex0%array = ex%array 
   ey0%array = ey%array 
   bz0%array = bz%array

end subroutine rksetup

subroutine rkstage(coef)

   sll_real64, intent(in) :: coef

   ex%array = ex0%array + coef * dt * dx%array
   ey%array = ey0%array + coef * dt * dy%array
   bz%array = bz0%array + coef * dt * dz%array

end subroutine rkstage

subroutine accumulate(coef)

   sll_real64, intent(in) :: coef

   sx%array = sx%array + coef * dx%array
   sy%array = sy%array + coef * dy%array
   sz%array = sz%array + coef * dz%array

end subroutine accumulate

subroutine rkstep()

   ex%array = ex0%array + dt * sx%array
   ey%array = ey0%array + dt * sy%array
   bz%array = bz0%array + dt * sz%array

end subroutine rkstep

end program test_maxwell_2d_diga_wave
