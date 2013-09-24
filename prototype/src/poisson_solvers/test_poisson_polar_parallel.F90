program test_poisson_polar_parallel
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_constants.h"

use sll_collective
use sll_gnuplot_parallel
use sll_poisson_polar_parallel

implicit none

type(sll_poisson_polar) :: poisson
sll_real64, dimension(:,:), allocatable :: rhs
sll_real64, dimension(:,:), allocatable :: phi
sll_real64, dimension(:,:), allocatable :: phi_cos
sll_real64, dimension(:,:), allocatable :: phi_sin
sll_real64, dimension(:),   allocatable :: r
sll_real64, dimension(:),   allocatable :: a
sll_real64, dimension(:,:), allocatable :: x
sll_real64, dimension(:,:), allocatable :: y


type(layout_2D), pointer :: layout_r ! sequential in r direction
type(layout_2D), pointer :: layout_a ! sequential in theta direction

sll_int32, dimension(2)   :: global
sll_int32                 :: gi, gj
sll_int32  :: i, j
sll_int32  :: nc_r, nr
sll_int32  :: nc_a, na
sll_real64 :: r_min, r_max, delta_r
sll_real64 :: a_min, a_max, delta_a
sll_int32  :: error
sll_int32  :: psize
sll_int32  :: prank
sll_int32  :: nr_loc
sll_int32  :: na_loc

sll_int32, parameter  :: n = 4

print*,'Testing the Poisson solver in 2D, polar coordinate'

r_min   = 1.0_f64
r_max   = 2.0_f64

a_min   = 0.0_f64
a_max   = 2.0_f64 * sll_pi

nc_r    = 256
nc_a    = 1024
nr      = nc_r+1
na      = nc_a+1
delta_r = (r_max-r_min)/real(nr-1,f64)
delta_a = 2.0_f64*sll_pi/real(na-1,f64)

!Boot parallel environment
call sll_boot_collective()

psize  = sll_get_collective_size(sll_world_collective)
prank  = sll_get_collective_rank(sll_world_collective)

layout_r => new_layout_2D( sll_world_collective )
layout_a => new_layout_2D( sll_world_collective )

call initialize_layout_with_distributed_2D_array( nr,     &
                                                  na,     &
                                                  1,      &
                                                  psize,  &
                                                  layout_r )

call initialize_layout_with_distributed_2D_array( nr,     &
                                                  na,     &
                                                  psize,  &
                                                  1,      &
                                                  layout_a )
call flush(6)
if (prank == 0) then
   call sll_view_lims_2D(layout_a)
   call sll_view_lims_2D(layout_a)
end if
call flush(6)

call compute_local_sizes(layout_a, nr_loc, na_loc )
SLL_CLEAR_ALLOCATE(rhs(1:nr_loc,1:na_loc),error)
SLL_CLEAR_ALLOCATE(phi(1:nr_loc,1:na_loc),error)
SLL_CLEAR_ALLOCATE(phi_cos(1:nr_loc,1:na_loc),error)
SLL_CLEAR_ALLOCATE(phi_sin(1:nr_loc,1:na_loc),error)
SLL_CLEAR_ALLOCATE(r(1:nr_loc),error)
SLL_CLEAR_ALLOCATE(a(1:na_loc),error)
SLL_CLEAR_ALLOCATE(x(1:nr_loc,1:na_loc),error)
SLL_CLEAR_ALLOCATE(y(1:nr_loc,1:na_loc),error)

do i = 1, nr_loc
   global = local_to_global_2D( layout_a, (/i, 1/))
   gi = global(1)
   r(i)=r_min+(gi-1)*delta_r
end do

do j = 1, na_loc
   global = local_to_global_2D( layout_a, (/1, j/))
   gj = global(2)
   a(j)=(gj-1)*delta_a
end do

do j=1,na_loc
   do i=1,nr_loc
      phi_cos(i,j) = (r(i)-r_min)*(r(i)-r_max)*cos(n*a(j))*r(i)
      phi_sin(i,j) = (r(i)-r_min)*(r(i)-r_max)*sin(n*a(j))*r(i)
      x(i,j)   = r(i)*cos(a(j))
      y(i,j)   = r(i)*sin(a(j))
   end do
end do

call initialize( poisson,   &
                 layout_r,      &
                 layout_a,      &
                 r_min,         &
                 r_max,         &
                 nc_r,          &
                 nc_a,          &
                 SLL_DIRICHLET, &
                 SLL_DIRICHLET)

do i =1,nr_loc
   do j=1,na_loc
      rhs(i,j) = - f_sin(r(i), a(j))
   end do
end do

call solve_poisson_polar(poisson, rhs, phi)

call sll_gnuplot_2d_parallel(x, y, phi_sin, 'phi_sin',  1, error)
call sll_gnuplot_2d_parallel(x, y, phi,     'solution', 1, error)

call error_max(phi_sin,phi,1e-4_f64)
call sll_halt_collective()

contains

subroutine error_max(phi, phi_exact, tolmax)

sll_real64, intent(in), dimension(:,:) :: phi
sll_real64, intent(in), dimension(:,:) :: phi_exact
sll_real64 :: errmax, tolmax 

errmax = maxval(abs(phi_exact-phi))
write(*,201) errmax
if ( errmax > tolmax ) then
   print*,'FAILED'
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
