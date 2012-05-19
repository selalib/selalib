program test_poisson_1d
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_field_1d.h"
#include "sll_poisson_solvers.h"
use numeric_constants
use sll_module_mapped_meshes_1d
use geometry_functions

use sll_poisson_1D_periodic

implicit none

type (sll_mapped_mesh_1d_analytic), target :: mesh_1d
class(sll_mapped_mesh_1d_base), pointer    :: m

type (scalar_field_1d)     :: ex
type (scalar_field_1d)     :: ex_exact
type (scalar_field_1d)     :: rho

type (poisson_1d_periodic)         :: poisson

sll_int32   :: nc_eta1
sll_real64  :: eta1_min, eta1_max
sll_real64  :: delta_eta1
sll_int32   :: error
sll_int32   :: mode
sll_int32   :: i

eta1_min = 0.0; eta1_max = 2*sll_pi;
nc_eta1 = 128

call initialize_mesh_1d_analytic( &
    mesh_1d,           &
    "mesh_1d",          &
    nc_eta1 + 1,    &
    linear_map_poisson_f,        &
    linear_map_poisson_jac_f )

call mesh_1d%write_to_file()
m => mesh_1d

call initialize_scalar_field_1d( &
    rho, &
    "rho", &
    m, &
    NODE_CENTERED_FIELD, &
    linear_map_poisson_f)

call initialize_scalar_field_1d( &
    ex, &
    "ex", &
    m, &
    NODE_CENTERED_FIELD, &
    linear_map_poisson_f)

call initialize_scalar_field_1d( &
    ex_exact, &
    "ex_exact", &
    m, &
    NODE_CENTERED_FIELD, &
    linear_map_poisson_f)

call new(poisson, nc_eta1, error) 

mode = 7
delta_eta1 = mesh_1d%delta_eta1
do i=1,nc_eta1+1
   rho%data(i)      =  mode**2*sin(mode*(i-1)*delta_eta1)
   ex_exact%data(i) = -mode*cos(mode*(i-1)*delta_eta1)
end do

call solve(poisson, ex, rho)
    
print*,'mode=',mode,'   error=',maxval(abs(FIELD_DATA(ex)-FIELD_DATA(ex_exact)))

end program test_poisson_1d
