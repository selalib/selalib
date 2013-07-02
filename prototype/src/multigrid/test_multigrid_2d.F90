program test_multigrid_cartesian
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_constants.h"
#include "sll_multigrid.h"

use sll_remapper
use sll_collective
use sll_gnuplot_parallel
use sll_boundary_condition_descriptors

implicit none

sll_int32                               :: ncx, ncy
sll_int32                               :: nx_loc, ny_loc
sll_int32                               :: error
sll_real64                              :: xmin, xmax
sll_real64                              :: ymin, ymax
sll_real64                              :: dx, dy
sll_real64, dimension(:,:), allocatable :: x, y
sll_real64, dimension(:,:), allocatable :: rho
sll_real64, dimension(:,:), allocatable :: phi
sll_real64, dimension(:,:), allocatable :: r
sll_int32                               :: i, j
sll_real64                              :: err
sll_int32, dimension(1:2)               :: global
sll_int32                               :: gi, gj
sll_int32                               :: myrank
type(layout_2D), pointer                :: layout
sll_real32                              :: ok 
sll_real32, dimension(1)                :: prod4test
type(sll_multigrid_solver_2d)           :: poisson_periodic
type(sll_multigrid_solver_2d)           :: poisson_dirichlet

ok = 1.0

!Boot parallel environment
call sll_boot_collective()

! Number of cells is equal to number of points in this case
ncx  = 64
ncy  = 64
xmin = 0.0_f64; xmax = 2*sll_pi
ymin = 0.0_f64; ymax = 2*sll_pi

myrank = sll_get_collective_rank(sll_world_collective)

dx = (xmax-xmin)/ncx
dy = (ymax-ymin)/ncy

! Layout and local sizes for FFTs in x-direction
layout => new_layout_2D( sll_world_collective )
call initialize_layout_with_distributed_2D_array( ncx, ncy, 2, 2, layout)
call sll_view_lims_2D( layout )

call compute_local_sizes_2d( layout, nx_loc, ny_loc )

SLL_CLEAR_ALLOCATE(rho(1:nx_loc+2,1:ny_loc+2), error)
SLL_CLEAR_ALLOCATE(phi(1:nx_loc+2,1:ny_loc+2), error)
SLL_CLEAR_ALLOCATE(r(1:nx_loc+2,1:ny_loc+2), error)
SLL_CLEAR_ALLOCATE(x(1:nx_loc,1:ny_loc), error)
SLL_CLEAR_ALLOCATE(y(1:nx_loc,1:ny_loc), error)

! initialize reference array
do j=1,ny_loc
   do i=1,nx_loc
      global = local_to_global_2D( layout, (/i, j/))
      gi = global(1)
      gj = global(2)
      x(i,j)  = (gi-1)*dx
      y(i,j)  = (gj-1)*dy
   end do
end do

call initialize(poisson_periodic, layout, xmin, xmax, ncx, SLL_PERIODIC, &
                                 ymin, ymax, ncy, SLL_PERIODIC )

rho(2:nx_loc+1,2:ny_loc+1)  = -2.0_f64*cos(x)*cos(y)
r = 1.
call solve(poisson_periodic, rho, phi, r)

call sll_gnuplot_curv_2d_parallel(x,y, &
     phi(2:nx_loc+1,2:ny_loc+1), "cos", 1, error)  
call sll_gnuplot_curv_2d_parallel(x,y, &
     cos(x)*cos(y)             , "cos", 2, error)  

err = sum(abs(cos(x)*cos(y)-phi(2:nx_loc+1,2:ny_loc+1)))/(ncx*ncy)

call flush(6); print*, ' ------------------'
call flush(6); print*, ' myrank ', myrank
call flush(6); print*, ' local average error:', err
call flush(6); print*, ' dx*dy =', dx*dy
call flush(6); print*, ' ------------------'

call delete(poisson_periodic)


call initialize(poisson_dirichlet, layout, &
                xmin, xmax, ncx, SLL_DIRICHLET, &
                ymin, ymax, ncy, SLL_DIRICHLET )

rho(2:nx_loc+1,2:ny_loc+1)  = -2.0_f64*sin(x)*sin(y)
r = 1.
call solve(poisson_dirichlet, rho, phi, r)

call sll_gnuplot_curv_2d_parallel(x,y, &
     phi(2:nx_loc+1,2:ny_loc+1), "resultat_sin", 1, error)  
call sll_gnuplot_curv_2d_parallel(x,y, &
     sin(x)*sin(y),              "solution_sin", 1, error)  

err = sum(abs(sin(x)*sin(y)-phi(2:nx_loc+1,2:ny_loc+1)))/(ncx*ncy)

call flush(6); print*, ' ------------------'
call flush(6); print*, ' myrank ', myrank
call flush(6); print*, ' local average error:', err
call flush(6); print*, ' dx*dy =', dx*dy
call flush(6); print*, ' ------------------'

call delete(poisson_dirichlet)

 
SLL_DEALLOCATE_ARRAY(phi, error)
SLL_DEALLOCATE_ARRAY(rho, error)

call sll_collective_reduce(sll_world_collective, (/ ok /), &
                           1, MPI_PROD, 0, prod4test )

if (myrank==0 .and. prod4test(1)==1.) then
   call flush(6)
   print*, 'PASSED'
else
   print*, 'Test stopped'
endif

call sll_halt_collective()

end program test_multigrid_cartesian
