program test_dg_fields
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_constants.h"
#include "sll_file_io.h"

use sll_logical_meshes
use sll_module_coordinate_transformations_2d
use sll_common_coordinate_transformations
use sll_dg_fields

implicit none

!=====================================!
! Simulation parameters               !
!=====================================!
sll_int32, parameter :: nc_eta1 = 10  !
sll_int32, parameter :: nc_eta2 = 20  !
sll_int32, parameter :: degree  = 3   !
!=====================================!

sll_real64 :: eta1_max, eta1_min
sll_real64 :: eta2_max, eta2_min
sll_real64 :: delta_eta1, delta_eta2

type(sll_logical_mesh_2d), pointer :: mesh
class(sll_coordinate_transformation_2d_analytic), pointer :: tau
class(sll_coordinate_transformation_2d_analytic), pointer :: collela

type(dg_field), pointer :: ex
type(dg_field), pointer :: bz

sll_real64, external :: gaussian, add

mesh => new_logical_mesh_2d(nc_eta1, nc_eta2, &
                            eta1_min=-1._f64, eta1_max=1._f64, &
                            eta2_min=-1._f64, eta2_max=1._f64)

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

collela => new_coordinate_transformation_2d_analytic( &
       "collela_transformation", &
       mesh, &
       sinprod_x1, &
       sinprod_x2, &
       sinprod_jac11, &
       sinprod_jac12, &
       sinprod_jac21, &
       sinprod_jac22, & 
       (/0.1_f64,0.1_f64,2.0_f64,2.0_f64/) )

call collela%write_to_file(SLL_IO_GNUPLOT)

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

ex => new_dg_field( degree, tau, add) 
call ex%write_to_file('ex', SLL_IO_GMSH)
call ex%write_to_file('ex', SLL_IO_MTV)
call ex%write_to_file('ex', SLL_IO_XDMF)

bz => new_dg_field( degree, collela, gaussian) 
call bz%write_to_file('bz', SLL_IO_GMSH)
call bz%write_to_file('bz', SLL_IO_MTV)
call bz%write_to_file('bz', SLL_IO_XDMF)

print*,'PASSED'

end program test_dg_fields
