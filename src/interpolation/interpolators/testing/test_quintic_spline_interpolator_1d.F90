program quintic_spline_interpolator_1d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_dirichlet

  use sll_m_quintic_spline_interpolator_1d, only: &
    sll_s_set_values_at_boundary, &
    sll_t_quintic_spline_interpolator_1d

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

type(sll_t_quintic_spline_interpolator_1d) :: spline

sll_real64                            :: error
sll_real64, allocatable, dimension(:) :: point
sll_real64, allocatable, dimension(:) :: pdata  
sll_real64, allocatable, dimension(:) :: fdata
sll_real64, allocatable, dimension(:) :: coord
sll_real64, allocatable, dimension(:) :: gdata

sll_int32 :: ierr, i
sll_int32, parameter :: n = 32
sll_int32, parameter :: m = 512
sll_real64  :: x_min, x_max, delta

SLL_ALLOCATE(coord(n), ierr)
SLL_ALLOCATE(pdata(n), ierr)
SLL_ALLOCATE(fdata(m), ierr)
SLL_ALLOCATE(gdata(m), ierr)
SLL_ALLOCATE(point(m), ierr)

print *, 'initialize data and interpolation_points array'
x_min = 0.0_f64
x_max = 1.0_f64
delta = (x_max - x_min ) / real(n-1,f64) 
do i=1,n
  coord(i) = x_min + (i-1)*delta
  pdata(i) = f(coord(i))
end do

print*, 'Quintic spline interpolation'
call spline%initialize(n, x_min, x_max, sll_p_dirichlet, sll_p_dirichlet )

call sll_s_set_values_at_boundary( spline,     &
                             f( x_min),  &
                             f( x_max),  &
                             df(x_min),  &
                             df(x_max))

call spline%compute_interpolants(pdata, coord, n )
  

delta = (x_max - x_min ) / real(m-1,f64) 
do i=1,m
  point(i) = x_min + (i-1)*delta
  gdata(i) = f(point(i))
  fdata(i) = spline%interpolate_from_interpolant_value(point(i))
  write(47,*) point(i), fdata(i), gdata(i)
end do

error = maxval(abs(fdata-gdata))
print*, 'Error=',error
print*, 'Successful, exiting program.'
print*, 'PASSED'

contains

function f(x)

  real(8) :: x
  real(8) :: f
  
  f = exp(x)
  
end function

function df(x)

  real(8) :: x
  real(8) :: df
  
  df = exp(x)
  
end function

end program quintic_spline_interpolator_1d
