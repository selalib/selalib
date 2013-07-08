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
sll_real64, dimension(:),   allocatable :: t
sll_real64, dimension(:,:), allocatable :: x
sll_real64, dimension(:,:), allocatable :: y

sll_int32  :: i, j
sll_int32  :: nc_r
sll_int32  :: nc_t
sll_real64 :: r_min, r_max, delta_r
sll_real64 :: t_min, t_max, delta_t
sll_int32  :: error

sll_int32, parameter  :: n = 4

type(layout_2D), pointer  :: layout_t
sll_int32, dimension(2)   :: global
sll_int32                 :: gi, gj
sll_int32                 :: psize
sll_int32                 :: prank
sll_int32                 :: e
sll_int32                 :: nproc_t
sll_int32                 :: nr_loc
sll_int32                 :: nt_loc

!Boot parallel environment
call sll_boot_collective()

psize    = sll_get_collective_size(sll_world_collective)
prank    = sll_get_collective_rank(sll_world_collective)

e        = int(log(real(psize))/log(2.))
print *, 'running on ', 2**e, 'processes'

! Layout and local sizes for FFTs in x-direction
layout_t => new_layout_2D( sll_world_collective )
nproc_t  = 2**e
nc_r     = 64
nc_t     = 128

r_min   = 1.0_f64
r_max   = 2.0_f64

t_min  = 0.0_f64
t_max  = 2.0_f64 * sll_pi

delta_r = (r_max-r_min)/real(nc_r-1,f64)
delta_t = 2.0_f64*sll_pi/real(nc_t-1,f64)

call initialize_layout_with_distributed_2D_array( nc_r+1,   &
                                                  nc_t+1,   &
                                                  nproc_t,  &
                                                    1,      &
                                                  layout_t )
call flush(6)
call sll_view_lims_2D(layout_t)
call flush(6)

call compute_local_sizes(layout_t, nr_loc, nt_loc )

SLL_CLEAR_ALLOCATE(rhs(1:nr_loc,1:nt_loc),error)
SLL_CLEAR_ALLOCATE(phi(1:nr_loc,1:nt_loc),error)
SLL_CLEAR_ALLOCATE(phi_cos(1:nr_loc,1:nt_loc),error)
SLL_CLEAR_ALLOCATE(phi_sin(1:nr_loc,1:nt_loc),error)
SLL_CLEAR_ALLOCATE(r(1:nr_loc),error)
SLL_CLEAR_ALLOCATE(t(1:nt_loc),error)
SLL_CLEAR_ALLOCATE(x(1:nr_loc,1:nt_loc),error)
SLL_CLEAR_ALLOCATE(y(1:nr_loc,1:nt_loc),error)

do i = 1, nr_loc
   global = local_to_global_2D( layout_t, (/i, 1/))
   gi = global(1)
   r(i)=r_min+(gi-1)*delta_r
end do

do j = 1, nt_loc
   global = local_to_global_2D( layout_t, (/1, j/))
   gj = global(2)
   t(j)=(gj-1)*delta_t
end do

do j=1,nt_loc
   do i=1,nr_loc
      phi_cos(i,j) = (r(i)-r_min)*(r(i)-r_max)*cos(n*t(j))*r(i)
      phi_sin(i,j) = (r(i)-r_min)*(r(i)-r_max)*sin(n*t(j))*r(i)
      x(i,j)   = r(i)*cos(t(j))
      y(i,j)   = r(i)*sin(t(j))
   end do
end do

call initialize( poisson, layout_t, &
                 r_min, r_max, nc_r, nc_t,    &
                 SLL_DIRICHLET, SLL_DIRICHLET)


do i =1,nr_loc
   do j=1,nt_loc
      rhs(i,j) = - f_sin(r(i), t(j))
   end do
end do


call solve(poisson, rhs, phi)

call sll_gnuplot_2d_parallel(x, y, phi_sin, 'phi_sin',  1, error)
call sll_gnuplot_2d_parallel(x, y, phi,     'solution', 1, error)

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


sll_real64 function f_cos( r, t )

   !sage: assume(r>=1)
   !sage: assume(r<=2)
   !sage: phi = (r-r_min)*(r-r_max)*r*cos(n*t)
   !sage: diff(r*diff(phi,r),r)/r + diff(phi,t,t)/(r*r)

   sll_real64 :: r
   sll_real64 :: t

   f_cos = -(r-r_max)*(r-r_min)*n*n*cos(n*t)/r &
           + ((r-r_max)*(r-r_min)*cos(n*t)  &
           + (r-r_max)*r*cos(n*t) + (r-r_min)*r*cos(n*t) &
           + 2*((r-r_max)*cos(n*t) + (r-r_min)*cos(n*t) &
           + r*cos(n*t))*r)/r


end function f_cos

sll_real64 function f_sin( r, t)

   !sage: assume(r>=1)
   !sage: assume(r<=2)
   !sage: phi = (r-r_min)*(r-r_max)*r*sin(n*t)
   !sage: diff(r*diff(phi,r),r)/r + diff(phi,t,t)/(r*r)

   sll_real64 :: r
   sll_real64 :: t
   
   f_sin = -(r-r_max)*(r-r_min)*n*n*sin(n*t)/r &
         + ((r-r_max)*(r-r_min)*sin(n*t) &
         + (r-r_max)*r*sin(n*t) + (r-r_min)*r*sin(n*t) &
         + 2*((r-r_max)*sin(n*t) + (r-r_min)*sin(n*t)  &
         + r*sin(n*t))*r)/r

end function f_sin

end program test_poisson_polar_parallel
