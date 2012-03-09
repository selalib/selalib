program test_poisson_1d
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_mesh_types.h"
#include "sll_poisson_solvers.h"
use numeric_constants

use sll_poisson_1D_periodic

implicit none

type (mesh_descriptor_1D), pointer :: geometry

type (field_1D_vec1), pointer      :: ex
type (field_1D_vec1), pointer      :: ex_exact
type (field_1D_vec1), pointer      :: rho

type (poisson_1d_periodic)         :: poisson

sll_int32   :: nc_eta1
sll_real64  :: eta1_min, eta1_max
sll_real64  :: delta_eta1
sll_int32   :: error
sll_int32   :: mode
sll_int32   :: i

eta1_min = 0.0; eta1_max = 2*sll_pi;
nc_eta1 = 128

geometry => new_mesh_descriptor_1D( eta1_min, eta1_max, nc_eta1, PERIODIC )
rho      => new_field_1D_vec1( geometry )
ex       => new_field_1D_vec1( geometry )
ex_exact => new_field_1D_vec1( geometry )

call new(poisson, nc_eta1, error) 

mode = 7
delta_eta1 = geometry%delta_eta1
do i=1,nc_eta1+1
   rho%data(i)      =  mode**2*sin(mode*(i-1)*delta_eta1)
   ex_exact%data(i) = -mode*cos(mode*(i-1)*delta_eta1)
end do

call solve(poisson, ex, rho)
    
print*,'mode=',mode,'   error=',maxval(abs(FIELD_DATA(ex)-FIELD_DATA(ex_exact)))

end program test_poisson_1d
