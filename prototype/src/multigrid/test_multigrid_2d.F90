program test_multigrid_cartesian
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_constants.h"
#include "sll_multigrid.h"

use sll_remapper
use sll_collective
use sll_gnuplot_parallel

implicit none

sll_int32                                    :: ncx, ncy
sll_int32                                    :: nx_loc, ny_loc
sll_int32                                    :: error
sll_real64                                   :: Lx, Ly
sll_real64                                   :: dx, dy
sll_real64                                   :: x, y
sll_real64, dimension(:,:), allocatable      :: rho
sll_real64, dimension(:,:), allocatable      :: phi_an
sll_real64, dimension(:,:), allocatable      :: phi
sll_real64, dimension(:,:), allocatable      :: r
sll_int32                                    :: i, j
sll_real64                                   :: average_err
sll_int32, dimension(1:2)                    :: global
sll_int32                                    :: gi, gj
sll_int32                                    :: myrank
type(layout_2D), pointer                     :: layout
sll_int64                                    :: colsz ! collective size
sll_int32                                    :: e
sll_real32                                   :: ok 
sll_real32, dimension(1)                     :: prod4test
type(sll_multigrid_solver_2d)                :: poisson

ok = 1.0

!Boot parallel environment
call sll_boot_collective()

! Number of cells is equal to number of points in this case
ncx = 512
ncy = 512
Lx  = 2.0*sll_pi
Ly  = 2.0*sll_pi

colsz  = sll_get_collective_size(sll_world_collective)
myrank = sll_get_collective_rank(sll_world_collective)

dx = Lx/ncx
dy = Ly/ncy

colsz  = sll_get_collective_size(sll_world_collective)
e      = int(log(real(colsz))/log(2.))-1
print *, 'running on ', 2**e, 'processes'

! Layout and local sizes for FFTs in x-direction
layout => new_layout_2D( sll_world_collective )
call initialize_layout_with_distributed_2D_array( ncx, ncy, 2, 2, layout)
call sll_view_lims_2D( layout )

call initialize_2d(poisson, layout, 0._f64, Lx, ncx, 0._f64, Ly, ncy)

call compute_local_sizes_2d( layout, nx_loc, ny_loc )

SLL_ALLOCATE(rho(nx_loc+2,ny_loc+2), error)
SLL_ALLOCATE(phi_an(nx_loc+2,ny_loc+2), error)
SLL_ALLOCATE(phi(nx_loc+2,ny_loc+2), error)
SLL_ALLOCATE(r(nx_loc+2,ny_loc+2), error)

! initialize reference array
do j=1,ny_loc
   do i=1,nx_loc
      global = local_to_global_2D( layout, (/i, j/))
      gi = global(1)
      gj = global(2)
      x  = (gi-1)*dx
      y  = (gj-1)*dy
      phi_an(i,j) =  cos(x)*sin(y) 
      rho(i,j)    = -2.0_f64*phi_an(i,j)
      r(i,j)      = 1.
   end do
end do

call solve_2d(poisson, rho, phi, r)

global = local_to_global_2D( layout, (/1, 1/))

call sll_gnuplot_rect_2d_parallel(dble(global(1)-1.0)*dx,dx, &
                                  dble(global(2)-1.0)*dy,dy, &
                                  phi, "potential", 1, error)  

average_err  = 0.d0
do j=1,ny_loc
   do i=1,nx_loc
      average_err  = average_err + abs(phi_an(i,j) - phi(i,j))
   enddo
enddo
     
average_err  = average_err/(ncx*ncy)

call flush(6); print*, ' ------------------'
call flush(6); print*, ' myrank ', myrank
call flush(6); print*, 'local average error:', average_err
call flush(6); print*, 'dx*dy =', dx*dy
call flush(6); print*, ' ------------------'

 
SLL_DEALLOCATE_ARRAY(phi, error)
SLL_DEALLOCATE_ARRAY(rho, error)
SLL_DEALLOCATE_ARRAY(phi_an, error)

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
