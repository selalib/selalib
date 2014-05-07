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
sll_int32, parameter :: nstep   = 1   !
sll_int32, parameter :: nc_eta1 = 20  !
sll_int32, parameter :: nc_eta2 = 20  !
sll_int32, parameter :: mode    = 2   !
sll_int32, parameter :: degree  = 2   !
!=====================================!

sll_real64 :: eta1_max, eta1_min
sll_real64 :: eta2_max, eta2_min
sll_real64 :: delta_eta1, delta_eta2

type(sll_logical_mesh_2d), pointer :: mesh
class(sll_coordinate_transformation_2d_analytic), pointer :: tau
class(sll_coordinate_transformation_2d_analytic), pointer :: colella

type(maxwell_2d_diga)   :: maxwell_TE

type(dg_field), pointer :: ex
type(dg_field), pointer :: ey
type(dg_field), pointer :: bz
type(dg_field), pointer :: bz_exact

sll_real64  :: time
sll_int32   :: istep
sll_real64  :: dt
sll_real64  :: cfl = 0.5
sll_real64  :: error

!init functions
sll_real64, external :: sol_bz, sol_ex, sol_ey

mesh => new_logical_mesh_2d(nc_eta1, nc_eta2, &
                            eta1_min=0._f64, eta1_max=2.*sll_pi, &
                            eta2_min=0._f64, eta2_max=2.*sll_pi)

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
       "collela_transformation", &
       mesh, &
       sinprod_x1, &
       sinprod_x2, &
       sinprod_jac11, &
       sinprod_jac12, &
       sinprod_jac21, &
       sinprod_jac22, & 
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

ex => new_dg_field( degree, tau, sol_ex) 
ey => new_dg_field( degree, tau, sol_ey) 
bz => new_dg_field( degree, tau, sol_bz) 
bz_exact => new_dg_field( degree, tau, sol_bz)

call bz%write_to_file('bz')

dt = cfl  / sqrt (1./(delta_eta1*delta_eta1)+1./(delta_eta2*delta_eta2))

time  = 0.0_f64

call initialize(maxwell_TE, tau, degree, TE_POLARIZATION)

do istep = 1, nstep !*** Loop over time


   call solve(maxwell_TE, ex, ey, bz, dt)

   time = time + dt
  
   call bz_exact%set_value(sol_bz, time)
   error = maxval(bz%array-bz_exact%array)

   write(*,"(10x,' istep = ',I6)",advance="no") istep
   write(*,"(' time = ',g12.3,' sec')",advance="no") time
   write(*,"(' erreur = ',g15.5)") error
   call ex%write_to_file('ex')
   call ey%write_to_file('ey')

   call ex%set_value(sol_ex, time)
   call ey%set_value(sol_ey, time)

   call ex%write_to_file('sol_ex')
   call ey%write_to_file('sol_ey')

end do ! next time step

print*,'PASSED'

end program test_maxwell_2d_diga

