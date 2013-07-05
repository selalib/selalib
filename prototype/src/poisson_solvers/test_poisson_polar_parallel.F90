program test_poisson_polar_parallel
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_constants.h"

use sll_poisson_polar_parallel
use sll_collective
use sll_remapper
use sll_gnuplot_parallel

implicit none

type(sll_poisson_polar) :: poisson
sll_real64, dimension(:,:), allocatable :: rhs
sll_real64, dimension(:,:), allocatable :: phi
sll_real64, dimension(:,:), allocatable :: phi_cos
sll_real64, dimension(:,:), allocatable :: phi_sin
sll_real64, dimension(:),   allocatable :: r
sll_real64, dimension(:),   allocatable :: theta
sll_real64, dimension(:,:), allocatable :: x
sll_real64, dimension(:,:), allocatable :: y

sll_int32  :: i, j
sll_int32  :: nr, nc_r
sll_int32  :: ntheta, nc_theta
sll_real64 :: r_min, r_max, delta_r
sll_real64 :: theta_min, theta_max, delta_theta
sll_int32  :: error

sll_int32, parameter  :: n = 4

type(layout_2D), pointer  :: layout_r
type(layout_2D), pointer  :: layout_theta
sll_int32, dimension(2)   :: global
sll_int32                 :: gi, gj
sll_int32                 :: psize
sll_int32                 :: prank
sll_int32                 :: e
sll_int32                 :: nproc_r
sll_int32                 :: nproc_theta
sll_int32                 :: nr_loc
sll_int32                 :: ntheta_loc

!Boot parallel environment
call sll_boot_collective()

psize = sll_get_collective_size(sll_world_collective)
prank = sll_get_collective_rank(sll_world_collective)

e      = int(log(real(psize))/log(2.))
print *, 'running on ', 2**e, 'processes'

! Layout and local sizes for FFTs in x-direction
layout_r     => new_layout_2D( sll_world_collective )
layout_theta => new_layout_2D( sll_world_collective )
nproc_r     = 2**e
nproc_theta = 2**e

nc_r = 64
nc_theta = 128

call initialize_layout_with_distributed_2D_array( nc_r,     &
                                                  nc_theta, &
                                                    1,      &
                                               nproc_r,     &
                                              layout_r )

call initialize_layout_with_distributed_2D_array( nc_r,     &
                                                  nc_theta, &
                                               nproc_theta, &
                                                    1,      &
                                              layout_theta )
call flush(6)
call sll_view_lims_2D(layout_r)
call flush(6)
call sll_view_lims_2D(layout_theta)
call flush(6)

r_min   = 1.0_f64
r_max   = 2.0_f64

theta_min  = 0.0_f64
theta_max  = 2.0_f64 * sll_pi

nr     = nc_r+1
ntheta = nc_theta+1
delta_r     = (r_max-r_min)/real(nr-1,f64)
delta_theta = 2.0_f64*sll_pi/real(ntheta-1,f64)

call compute_local_sizes(layout_r, nr_loc, ntheta_loc )
SLL_CLEAR_ALLOCATE(rhs(1:nr_loc,1:ntheta_loc),error)
SLL_CLEAR_ALLOCATE(phi(1:nr_loc,1:ntheta_loc),error)
SLL_CLEAR_ALLOCATE(phi_cos(1:nr_loc,1:ntheta_loc),error)
SLL_CLEAR_ALLOCATE(phi_sin(1:nr_loc,1:ntheta_loc),error)
SLL_CLEAR_ALLOCATE(r(1:nr_loc),error)
SLL_CLEAR_ALLOCATE(theta(1:ntheta_loc),error)

do i = 1, nr_loc
   global = local_to_global_2D( layout_r, (/i, 1/))
   gi = global(1)
   r(i)=r_min+(gi-1)*delta_r
end do
do j = 1, ntheta_loc
   global = local_to_global_2D( layout_r, (/1, j/))
   gj = global(2)
   theta(j)=(gj-1)*delta_theta
end do

SLL_CLEAR_ALLOCATE(x(1:nr,1:ntheta),error)
SLL_CLEAR_ALLOCATE(y(1:nr,1:ntheta),error)

do j=1,ntheta_loc
   do i=1,nr_loc
      phi_cos(i,j) = (r(i)-r_min)*(r(i)-r_max)*cos(n*theta(j))*r(i)
      phi_sin(i,j) = (r(i)-r_min)*(r(i)-r_max)*sin(n*theta(j))*r(i)
      x(i,j)   = r(i)*cos(theta(j))
      y(i,j)   = r(i)*sin(theta(j))
   end do
end do
call sll_gnuplot_2d_parallel(x, y, phi_sin, 'phi_sin', 1, error)
call sll_gnuplot_2d_parallel(x, y, phi_cos, 'phi_cos', 1, error)

call initialize( poisson,layout_r, layout_theta, &
                 r_min,r_max,nr-1,ntheta-1, DIRICHLET, DIRICHLET)

do i =1,nr_loc
   do j=1,ntheta_loc
      rhs(i,j) = - f_sin(r(i), theta(j))
   end do
end do

call solve(poisson, rhs, phi)

call error_max(phi_sin,phi,1e-4_f64)

contains

subroutine error_max(phi, phi_exact, tolmax)

sll_real64, intent(in), dimension(:,:) :: phi
sll_real64, intent(in), dimension(:,:) :: phi_exact
sll_real64 :: errmax, tolmax 

errmax = maxval(abs(phi_exact-phi))
write(*,201) errmax
if ( errmax > tolmax ) then
   print*,'FAILED'
   stop
else
   print*,'PASSED'
end if

201 format(' maximum error  =  ',e10.3)

end subroutine error_max


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

end program test_poisson_polar_parallel
