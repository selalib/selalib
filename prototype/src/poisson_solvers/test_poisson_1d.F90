program test_poisson_1d
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
!#include "sll_field_1d.h"
#include "sll_poisson_solvers.h"
use sll_constants
use sll_module_mapped_meshes_1d
use geometry_functions
use sll_splines
use sll_poisson_1D_periodic

implicit none

type (sll_mapped_mesh_1d_analytic), target :: mesh_1d
#ifdef STDF95
type(sll_mapped_mesh_1d_analytic), pointer    :: m
#else
class(sll_mapped_mesh_1d_base), pointer    :: m
#endif

sll_real64, dimension(:), allocatable :: ex
sll_real64, dimension(:), allocatable :: ex_exact
sll_real64, dimension(:), allocatable :: rho

type (poisson_1d_periodic)         :: poisson

sll_int32   :: nc_eta1
sll_real64  :: eta1_min, eta1_max
sll_int32   :: error
sll_int32   :: mode
sll_int32   :: i

nc_eta1 = 128

call initialize_mesh_1d_analytic( mesh_1d, "mesh_1d",        &
                                  nc_eta1 + 1,               &
                                  linear_map_poisson_f,      &
                                  linear_map_poisson_jac_f )

#ifdef STDF95
call mma_write_to_file(mesh_1d)
#else
call mesh_1d%write_to_file()
#endif
m => mesh_1d

SLL_ALLOCATE(rho(nc_eta1+1),error)
SLL_ALLOCATE(ex(nc_eta1+1),error)
SLL_ALLOCATE(ex_exact(nc_eta1+1),error)

mode = 4
do i=1,nc_eta1+1
#ifdef STDF95
   rho(i)      =  mode**2*sin(mode*x1_node_analytic_1d(mesh_1d, i))
   ex_exact(i) = -mode*cos(mode*x1_node_analytic_1d(mesh_1d, i))
#else
   rho(i)      =  mode**2*sin(mode*mesh_1d%x1_at_node(i))
   ex_exact(i) = -mode*cos(mode*mesh_1d%x1_at_node(i))
#endif
end do

#ifdef STDF95
eta1_min = x1_node_analytic_1d(mesh_1d, 1)
eta1_max = x1_node_analytic_1d(mesh_1d, nc_eta1+1)
#else
eta1_min = mesh_1d%x1_at_node(1)
eta1_max = mesh_1d%x1_at_node(nc_eta1+1)
#endif
call new(poisson, eta1_min, eta1_max, nc_eta1, error) 

call solve(poisson, ex, rho)
    
print*,'mode=',mode,'   error=',maxval(abs(ex-ex_exact))

end program test_poisson_1d
