!>$L_x,L_y$ domain dimensions and M,N are integers.
!>$
!>\omega = \sqrt{(\frac{M\pi}{L_x})^2+(\frac{N\pi}{L_y})^2}
!>$
!>$
!>B_z(x,y,t) =   - \cos(M \pi \frac{x}{L_x})  \cos(N \pi \frac{y}{L_y}) \cos(\omega t) 
!>$
!>$
!>E_x(x,y,t) = \frac{c^2 N \pi }{\omega Ly} cos(M \pi \frac{x}{L_x}) \sin(N \pi  \frac{y}{L_y}) \sin(\omega t) 
!>$
!>$
!>E_y(x,y,t) = - \frac{c^2 M \pi }{\omega Lx} \sin (M \pi \frac{x}{L_x}) \cos (N \pi  \frac{y}{L_y}) \sin(\omega t) 
!>$
!>

program test_maxwell_2d_diga
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

implicit none

!=====================================!
! Simulation parameters               !
!=====================================!
sll_int32, parameter :: nc_eta1 = 2   !
sll_int32, parameter :: nc_eta2 = 2   !
sll_int32, parameter :: degree  = 1   !
!=====================================!

sll_int32  :: nstep
sll_real64 :: eta1_max, eta1_min
sll_real64 :: eta2_max, eta2_min
sll_real64 :: delta_eta1, delta_eta2

type(sll_logical_mesh_2d), pointer :: mesh
class(sll_coordinate_transformation_2d_analytic), pointer :: tau
class(sll_coordinate_transformation_2d_analytic), pointer :: colella

type(maxwell_2d_diga)   :: maxwell_TE

type(dg_field), pointer :: ex, ex0, dx, sx
type(dg_field), pointer :: ey, ey0, dy, sy
type(dg_field), pointer :: bz, bz0, dz, sz
type(dg_field), pointer :: exact

sll_real64  :: time
sll_int32   :: istep
sll_real64  :: dt
sll_real64  :: cfl = 0.5
sll_real64  :: error

!init functions
sll_real64, external :: sol_bz, sol_ex, sol_ey, uniform_x, uniform_y

mesh => new_logical_mesh_2d(nc_eta1, nc_eta2, &
                            eta1_min=-0._f64, eta1_max=2._f64, &
                            eta2_min=-0._f64, eta2_max=2._f64)

write(*,"(3f8.3,i4)") mesh%eta1_min,mesh%eta1_max,mesh%delta_eta1,mesh%num_cells1
write(*,"(3f8.3,i4)") mesh%eta2_min,mesh%eta2_max,mesh%delta_eta2,mesh%num_cells2

eta1_min = mesh%eta1_min
eta1_max = mesh%eta1_max
eta2_min = mesh%eta2_min
eta2_max = mesh%eta2_max

delta_eta1 = mesh%delta_eta1
delta_eta2 = mesh%delta_eta2

! "Colella transformation";
! sinusoidal product (see P. Colella et al. JCP 230 (2011) formula 
! (102) p 2968):
!
! x1 = eta1 + 0.1 * sin(2*pi*eta1) * sin(2*pi*eta2)
! x2 = eta2 + 0.1 * sin(2*pi*eta1) * sin(2*pi*eta2)

colella => new_coordinate_transformation_2d_analytic( &
       "collela_transformation",                      &
       mesh,                                          &
       sinprod_x1,                                    &
       sinprod_x2,                                    &
       sinprod_jac11,                                 &
       sinprod_jac12,                                 &
       sinprod_jac21,                                 &
       sinprod_jac22,                                 &  
       (/0.1_f64,0.1_f64,1.0_f64,1.0_f64/) )

call colella%write_to_file(SLL_IO_GNUPLOT)

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

call tau%write_to_file(SLL_IO_MTV)

time  = 0.0_f64

ex => new_dg_field(degree,tau,uniform_x) 
ey => new_dg_field(degree,tau,uniform_y) 
bz => new_dg_field(degree,tau,uniform_x) 

ex0 => new_dg_field(degree,tau) 
ey0 => new_dg_field(degree,tau) 
bz0 => new_dg_field(degree,tau) 

dx => new_dg_field(degree,tau) 
dy => new_dg_field(degree,tau) 
dz => new_dg_field(degree,tau) 

sx => new_dg_field(degree,tau) 
sy => new_dg_field(degree,tau) 
sz => new_dg_field(degree,tau) 

exact => new_dg_field( degree, tau)

dt = cfl/sqrt(1./(delta_eta1*delta_eta1)+1./(delta_eta2*delta_eta2))

dt = 0.1
nstep = ceiling(sll_pi/dt)
nstep = 20
print*, 'dt = ', dt

call initialize(maxwell_TE, tau, degree, TE_POLARIZATION, &
                SLL_DIRICHLET, SLL_DIRICHLET )

call ex%write_to_file('ex')
call ey%write_to_file('ey')
call bz%write_to_file('bz')

do istep = 1, nstep !*** Loop over time

   write(*,"(10x,' istep = ',I6,' time = ',g12.3,' sec')") istep, time

   call rksetup()

   call solve(maxwell_TE, ex, ey, bz, dx, dy, dz)
   call dx%write_to_file('dx')
   call dy%write_to_file('dy')
   call dz%write_to_file('dz')
   stop
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

   !call check_error_ex(time)
   !call check_error_ey(time)
   !call check_error_bz(time)


end do ! next time step

print*, maxval(bz%array)-exp(1.0_f64)

contains

subroutine check_error_ex(time)

   sll_real64, intent(in) :: time

   call exact%set_value(sol_ex, time)

   error = maxval(abs(ex%array-exact%array))
   write(*,"(' EX erreur = ',g15.5)",advance="no") error

end subroutine check_error_ex

subroutine check_error_ey(time)

   sll_real64, intent(in) :: time

   call exact%set_value(sol_ey, time)

   error = maxval(abs(ey%array-exact%array))
   write(*,"(' EY erreur = ',g15.5)",advance="no") error

end subroutine check_error_ey

subroutine check_error_bz(time)

   sll_real64, intent(in) :: time

   call exact%set_value(sol_bz, time)

   error = maxval(abs(bz%array-exact%array))
   write(*,"(' BZ erreur = ',g15.5)",advance="no") error

end subroutine check_error_bz


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

end program test_maxwell_2d_diga

