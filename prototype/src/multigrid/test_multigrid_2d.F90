#define MPI_MASTER 0

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
sll_real64, dimension(:,:), allocatable :: cof
sll_int32                               :: i, j
sll_real64                              :: err
sll_int32, dimension(1:2)               :: global
sll_int32                               :: gi, gj
sll_int32                               :: myrank
type(layout_2D), pointer                :: layout
type(sll_multigrid_solver_2d)           :: poisson_periodic
type(sll_multigrid_solver_2d)           :: poisson_dirichlet
type(sll_multigrid_solver_2d)           :: poisson_polar
sll_real64, dimension(:), allocatable   :: r
sll_real64, dimension(:), allocatable   :: theta
sll_int32                               :: nr, nr_loc
sll_int32                               :: ntheta, ntheta_loc
sll_real64                              :: r_min, r_max, delta_r
sll_real64                              :: theta_min, theta_max, delta_theta
sll_int32, parameter                    :: n=2

!Boot parallel environment
call sll_boot_collective()

! Number of cells is equal to number of points in this case
ncx  = 512
ncy  = 512
xmin = 0.0_f64; xmax = 2*sll_pi
ymin = 0.0_f64; ymax = 2*sll_pi

myrank = sll_get_collective_rank(sll_world_collective)

dx = (xmax-xmin)/ncx
dy = (ymax-ymin)/ncy

layout => new_layout_2D( sll_world_collective )
call initialize_layout_with_distributed_2D_array( ncx, ncy, 2, 2, layout)

if (myrank == MPI_MASTER) call sll_view_lims_2D( layout )

call compute_local_sizes_2d( layout, nx_loc, ny_loc )

SLL_CLEAR_ALLOCATE(rho(1:nx_loc+2,1:ny_loc+2), error)
SLL_CLEAR_ALLOCATE(phi(1:nx_loc+2,1:ny_loc+2), error)
SLL_CLEAR_ALLOCATE(cof(1:nx_loc+2,1:ny_loc+2), error)
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

call test_periodic()
!call test_dirichlet()

call sll_delete(layout)

!SLL_DEALLOCATE_ARRAY(phi, error)
!SLL_DEALLOCATE_ARRAY(rho, error)
!SLL_DEALLOCATE_ARRAY(cof, error)
!
!r_min       = 1.0_f64
!r_max       = 2.0_f64
!theta_min   = 0.0_f64
!theta_max   = 2.0_f64 * sll_pi
!nr          = 32
!ntheta      = 128
!delta_r     = (r_max-r_min)/nr
!delta_theta = 2*sll_pi/ntheta
!
!layout => new_layout_2D( sll_world_collective )
!call initialize_layout_with_distributed_2D_array( nr, ntheta, 2, 2, layout)
!
!if (myrank == MPI_MASTER) then
!   print*, ' polar '
!   call sll_view_lims_2D( layout )
!end if
!call flush(6)
!call compute_local_sizes_2d( layout, nr_loc, ntheta_loc )
!
!SLL_CLEAR_ALLOCATE(rho(1:nr_loc+2,1:ntheta_loc+2),error)
!SLL_CLEAR_ALLOCATE(phi(1:nr_loc+2,1:ntheta_loc+2),error)
!SLL_CLEAR_ALLOCATE(cof(1:nr_loc+2,1:ntheta_loc+2),error)
!SLL_CLEAR_ALLOCATE(r(1:nr_loc+1),error)
!SLL_CLEAR_ALLOCATE(theta(1:ntheta_loc+1),error)
!
!do i = 1, nr_loc
!   global = local_to_global_2D( layout, (/i, 1/))
!   gi = global(1)
!   r(i)=r_min+(gi-1)*delta_r
!end do
!
!do j = 1, ntheta_loc
!   global = local_to_global_2D( layout, (/1, j/))
!   gj = global(2)
!   theta(j)=(gj-1)*delta_theta
!end do
!
!call test_polar()
!
!SLL_DEALLOCATE_ARRAY(phi, error)
!SLL_DEALLOCATE_ARRAY(rho, error)
!SLL_DEALLOCATE_ARRAY(cof, error)

call flush(6)
if (myrank==0) then
   print*, 'PASSED'
endif

call sll_halt_collective()

contains

subroutine test_periodic()

call initialize(poisson_periodic, layout, xmin, xmax, ncx, SLL_PERIODIC, &
                                 ymin, ymax, ncy, SLL_PERIODIC )

rho(2:nx_loc+1,2:ny_loc+1)  = -2.0_f64*cos(x)*cos(y)
cof = 1.
call solve(poisson_periodic, rho, phi, cof)

call sll_gnuplot_curv_2d_parallel(x,y, &
     phi(2:nx_loc+1,2:ny_loc+1), "phi", 1, error)  
call sll_gnuplot_curv_2d_parallel(x,y, &
     cos(x)*cos(y)             , "cos", 1, error)  

err = sum(abs(cos(x)*cos(y)-phi(2:nx_loc+1,2:ny_loc+1)))/(ncx*ncy)

call flush(6); print*, ' - periodic -------'
call flush(6); print*, ' myrank ', myrank
call flush(6); print*, ' local average error:', err
call flush(6); print*, ' dx*dy =', dx*dy
call flush(6); print*, ' ------------------'

call delete(poisson_periodic)

end subroutine test_periodic

subroutine test_dirichlet()

call initialize(poisson_dirichlet, layout, &
                xmin, xmax, ncx, SLL_DIRICHLET, &
                ymin, ymax, ncy, SLL_DIRICHLET )

rho(2:nx_loc+1,2:ny_loc+1)  = -2.0_f64*sin(x)*sin(y)
cof = 1.0_f64
call solve(poisson_dirichlet, rho, phi, cof)

call sll_gnuplot_curv_2d_parallel(x,y, &
     phi(2:nx_loc+1,2:ny_loc+1), "resultat_sin", 1, error)  
call sll_gnuplot_curv_2d_parallel(x,y, &
     sin(x)*sin(y),              "solution_sin", 1, error)  

err = sum(abs(sin(x)*sin(y)-phi(2:nx_loc+1,2:ny_loc+1)))/(ncx*ncy)

call flush(6); print*, ' - dirichlet ------'
call flush(6); print*, ' myrank ', myrank
call flush(6); print*, ' local average error:', err
call flush(6); print*, ' dx*dy =', dx*dy
call flush(6); print*, ' ------------------'

call delete(poisson_dirichlet)

end subroutine test_dirichlet

subroutine test_polar()

call initialize(poisson_polar, layout,                       &
                r_min,     r_max,     nr,     SLL_DIRICHLET, &
                theta_min, theta_max, ntheta, SLL_PERIODIC )

do i =1,nr_loc
   do j=1,ntheta_loc
      rho(i,j) = f_sin(r(i), theta(j))
      cof(i,j) = r(i)
   end do
end do

call solve(poisson_polar, rho, phi, cof)

call sll_gnuplot_curv_2d_parallel(x,y, &
     phi(2:nx_loc+1,2:ny_loc+1), "cos", 1, error)  
call sll_gnuplot_curv_2d_parallel(x,y, &
     cos(x)*cos(y)             , "cos", 2, error)  

err = 0.0_f64
do j=1,ntheta
   do i=1,nr
      err = err + abs(phi(i,j)-(r(i)-r_min)*(r(i)-r_max)*sin(n*theta(j))*r(i))
   end do
end do

call flush(6); print*, ' -- polar ---------'
call flush(6); print*, ' myrank ', myrank
call flush(6); print*, ' local average error:', err/(nr*ntheta)
call flush(6); print*, ' dx*dy =', dx*dy
call flush(6); print*, ' ------------------'

call delete(poisson_polar)

end subroutine test_polar

sll_real64 function f_cos( r, theta )

   !sage: assume(r>=1)
   !sage: assume(r<=2)
   !sage: phi = (r-r_min)*(r-r_max)*r*cos(n*theta)
   !sage: diff(r*diff(phi,r),r)/r + diff(phi,theta,theta)/(r*r)

   sll_real64 :: r
   sll_real64 :: theta

   f_cos = -(r-r_max)*(r-r_min)*n*n*cos(n*theta)/r &
           + ((r-r_max)*(r-r_min)*cos(n*theta)  &
           + (r-r_max)*r*cos(n*theta) + (r-r_min)*r*cos(n*theta) &
           + 2*((r-r_max)*cos(n*theta) + (r-r_min)*cos(n*theta) &
           + r*cos(n*theta))*r)/r


end function f_cos

sll_real64 function f_sin( r, theta)

   !sage: assume(r>=1)
   !sage: assume(r<=2)
   !sage: phi = (r-r_min)*(r-r_max)*r*sin(n*theta)
   !sage: diff(r*diff(phi,r),r)/r + diff(phi,theta,theta)/(r*r)

   sll_real64 :: r
   sll_real64 :: theta
   
   f_sin = -(r-r_max)*(r-r_min)*n*n*sin(n*theta)/r &
         + ((r-r_max)*(r-r_min)*sin(n*theta) &
         + (r-r_max)*r*sin(n*theta) + (r-r_min)*r*sin(n*theta) &
         + 2*((r-r_max)*sin(n*theta) + (r-r_min)*sin(n*theta)  &
         + r*sin(n*theta))*r)/r

end function f_sin

end program test_multigrid_cartesian
