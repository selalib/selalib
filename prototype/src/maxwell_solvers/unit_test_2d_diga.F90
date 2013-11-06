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

use sll_logical_meshes
use sll_module_coordinate_transformations_2d
use sll_common_coordinate_transformations
use sll_maxwell_2d_diga

implicit none

sll_real64 :: eta1_max, eta1_min
sll_real64 :: eta2_max, eta2_min
sll_real64 :: delta_eta1, delta_eta2

sll_int32  :: nc_eta1, nc_eta2
sll_int32  :: error

type(sll_logical_mesh_2d), pointer :: mesh
class(sll_coordinate_transformation_2d_analytic), pointer :: tau

type(maxwell_2d_diga)                   :: maxwell_TE
type(maxwell_2d_diga)                   :: maxwell_TM

sll_real64, dimension(:,:), allocatable :: ex
sll_real64, dimension(:,:), allocatable :: ey
sll_real64, dimension(:,:), allocatable :: bz
sll_real64, dimension(:,:), allocatable :: bz_exact

sll_real64, dimension(:,:), allocatable :: bx
sll_real64, dimension(:,:), allocatable :: by
sll_real64, dimension(:,:), allocatable :: ez
sll_real64, dimension(:,:), allocatable :: ez_exact

sll_int32                               :: i, j
sll_real64                              :: omega
sll_real64                              :: time
sll_int32                               :: istep, nstep
sll_real64                              :: err_te
sll_real64                              :: err_tm
sll_real64                              :: dt
sll_real64                              :: cfl = 0.5
sll_int32                               :: degree = 2
sll_int32,  parameter                   :: mode = 2

nc_eta1 = 2; nc_eta2 = 1

!mesh => new_logical_mesh_2d(nc_eta1, nc_eta2)
mesh => new_logical_mesh_2d(nc_eta1, nc_eta2, &
                            eta1_min=-1._f64, eta2_min=-1._f64)

write(*,*) mesh%eta1_min,mesh%eta1_max,mesh%delta_eta1,mesh%num_cells1
write(*,*) mesh%eta2_min,mesh%eta2_max,mesh%delta_eta2,mesh%num_cells2

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

!tau => new_coordinate_transformation_2d_analytic( &
!       "collela_transformation", &
!       mesh, &
!       sinprod_x1, &
!       sinprod_x2, &
!       sinprod_jac11, &
!       sinprod_jac12, &
!       sinprod_jac21, &
!       sinprod_jac22 )

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

call initialize(maxwell_TE, tau, degree, TE_POLARIZATION)

stop
call initialize(maxwell_TM, tau, degree, TM_POLARIZATION)

dt = cfl  / sqrt (1./(delta_eta1*delta_eta1)+1./(delta_eta2*delta_eta2))
nstep = 100

time  = 0.

omega = sqrt( (mode*sll_pi/(nc_eta1*delta_eta1))**2   &
        &    +(mode*sll_pi/(nc_eta2*delta_eta2))**2)

SLL_CLEAR_ALLOCATE(ex(1:nc_eta1+1,1:nc_eta2+1),error)
SLL_CLEAR_ALLOCATE(ey(1:nc_eta1+1,1:nc_eta2+1),error)
SLL_CLEAR_ALLOCATE(bz(1:nc_eta1+1,1:nc_eta2+1),error)

SLL_CLEAR_ALLOCATE(bx(1:nc_eta1+1,1:nc_eta2+1),error)
SLL_CLEAR_ALLOCATE(by(1:nc_eta1+1,1:nc_eta2+1),error)
SLL_CLEAR_ALLOCATE(ez(1:nc_eta1+1,1:nc_eta2+1),error)

SLL_CLEAR_ALLOCATE(bz_exact(1:nc_eta1+1,1:nc_eta2+1),error)
SLL_CLEAR_ALLOCATE(ez_exact(1:nc_eta1+1,1:nc_eta2+1),error)

do istep = 1, nstep !*** Loop over time

   time = time + 0.5_f64*dt

   do j = 1, nc_eta2+1
   do i = 1, nc_eta1+1
      bz_exact(i,j) =   - cos(mode*sll_pi*(i-0.5_f64)*delta_eta1)    &
                        * cos(mode*sll_pi*(j-0.5_f64)*delta_eta2)    &
                        * cos(omega*time)
   end do
   end do

   ez_exact = - bz_exact

   if (istep == 1) then

      bz = bz_exact
      ez = ez_exact

   else

      err_te = maxval(bz-bz_exact)
      err_tm = maxval(ez-ez_exact)
   
      write(*,"(10x,' istep = ',I6)",advance="no") istep
      write(*,"(' time = ',g12.3,' sec')",advance="no") time
      write(*,"(' erreur L2 = ',2g15.5)") err_te, err_tm

   end if

   !call plot_fields('ez',ez, ez_exact, istep, time)

   call solve(maxwell_TE, ex, ey, bz, dt)
   call solve(maxwell_TM, bx, by, ez, dt)

   time = time + 0.5_f64*dt


end do ! next time step

print*,'PASSED'

DEALLOCATE(ex)
DEALLOCATE(ey)
DEALLOCATE(bz)
DEALLOCATE(bz_exact)
call delete(mesh)

end program test_maxwell_2d_diga
