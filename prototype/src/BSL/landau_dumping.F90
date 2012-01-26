program landau_dumping
#include "sll_assert.h"
#include "sll_working_precision.h"
#include "sll_mesh_types.h"
#include "sll_memory.h"
#include "sll_poisson_solvers.h"

use numeric_constants
use distribution_function
use sll_diagnostics
use sll_bsl
use sll_poisson_2d_periodic

implicit none
  
sll_int32  :: nc_eta1, nc_eta2, nc_eta3, nc_eta4
sll_int32  :: i_step, n_step
sll_real64 :: eta1_min, eta1_max, eta2_min, eta2_max
sll_real64 :: eta3_min, eta3_max, eta4_min, eta4_max
sll_real64 :: eta1, eta2, eta3, eta4
sll_int32  :: ix, jx, iv, jv

type(geometry_2D), pointer :: geom_x, geom_v
type(mesh_descriptor_2D), pointer :: mesh_x, mesh_v
type(sll_distribution_function_4D_t), pointer :: dist_func
type(bsl_workspace_4d), pointer :: bsl_work

procedure(scalar_function_2D), pointer :: x1, x2, v1, v2
sll_real64 :: delta_x1, delta_x2, delta_v1, delta_v2
sll_real64 :: delta_t

type(poisson_2d_periodic), pointer :: poisson
type(field_2D_vec2),       pointer :: exy, exy_exact
type(field_2D_vec1),       pointer :: rho
sll_int32                          :: error
character(len=4) :: counter

eta1_min =  0.0_f64; eta1_max =  2.0_f64 * sll_pi
eta2_min =  0.0_f64; eta2_max =  2.0_f64 * sll_pi

geom_x => new_geometry_2D('cartesian')

nc_eta1 = 64; nc_eta2 = 64

mesh_x => new_mesh_descriptor_2D(eta1_min, eta1_max, nc_eta1, &
          PERIODIC, eta2_min, eta2_max, nc_eta2, PERIODIC, geom_x)

call write_mesh_2D(mesh_x,"mesh_x")

rho       => new_field_2D_vec1(mesh_x)
exy       => new_field_2D_vec2(mesh_x)
exy_exact => new_field_2D_vec2(mesh_x)

poisson   => new_poisson_2d_periodic(exy)

eta3_min = -6.0_f64; eta3_max =  6.0_f64 
eta4_min = -6.0_f64; eta4_max =  6.0_f64 

geom_v => new_geometry_2D('cartesian')
nc_eta3 = 31; nc_eta4 = 31

mesh_v => new_mesh_descriptor_2D(eta3_min, eta3_max, nc_eta3, &
          PERIODIC, eta4_min, eta4_max, nc_eta4, PERIODIC, geom_v)

dist_func => sll_new_distribution_function_4D(mesh_x, mesh_v, NODE_CENTERED_DF, 'f')

call sll_init_distribution_function_4D( dist_func, LANDAU)

x1 => dist_func%field%descriptor_x%geom%x1
x2 => dist_func%field%descriptor_x%geom%x2
v1 => dist_func%field%descriptor_v%geom%x1
v2 => dist_func%field%descriptor_v%geom%x2

delta_x1 = GET_MESH_DELTA_ETA1(mesh_x)
delta_x2 = GET_MESH_DELTA_ETA2(mesh_x)

delta_v1 = GET_MESH_DELTA_ETA1(mesh_v)
delta_v2 = GET_MESH_DELTA_ETA2(mesh_v)

call write_distribution_function(dist_func)

! initialize BSL  
bsl_work => new_bsl_workspace( dist_func)

! run BSL method using 10 time steps and second order splitting
n_step = 10
delta_t = 1.0_f64/n_step
do i_step = 1, n_step
   call bsl_second_order(bsl_work, dist_func, exy, exy, delta_t)
   call write_distribution_function(dist_func)
   call compute_rho(dist_func, rho)
   call write_vec1d(rho%data,mesh_x%nc_eta1+1,mesh_x%nc_eta2+1,"rho","mesh_x",0)
   call int2string(dist_func%plot_counter,counter)
   call solve_poisson_2d_periodic(poisson,exy,rho,error)
   call write_vec2d(exy%data%v1,exy%data%v2,mesh_x%nc_eta1+1, &
                    mesh_x%nc_eta2+1,"exy"//counter,"mesh_x",0)
end do

call delete_poisson_2d_periodic(poisson)
call delete_field_2D_vec1( rho )
call delete_field_2D_vec2( exy )
call delete_field_2D_vec2( exy_exact )
print *, 'Successful, exiting program.'
  
end program landau_dumping
