program diagnostics_tester
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"

  use sll_splines
  use numeric_constants
  use sll_diagnostics
  implicit none
  sll_int32 :: nx, nv, i, j
  sll_real64, allocatable, dimension(:,:) :: x
  sll_real64, allocatable, dimension(:,:) :: v
  sll_real64, allocatable, dimension(:,:) :: f
  sll_real64 :: angle, xt, vt, pi, R
  
  pi = 4.*atan(1.)

  nx = 128
  nv = 64

  ! Create the coordinate data.
  allocate(x(nx,nv))
  allocate(v(nx,nv))

  do j = 1, nv
     vt = real(j-1)/(nv-1)
     angle = vt * 2. * pi
     do i = 1, nx
        xt = real(i-1) / float(nx-1)
        R = (1.-xt)*2. + xt*5.
        x(i,j) = R * cos(angle)
        v(i,j) = R * sin(angle)
     end do
  end do 

  ! Create the scalar data.
  allocate(f(nx-1,nv-1))
  do j = 1, nv-1
     do i = 1, nx-1
        f(i,j) = i*sin(float(j-1))
     end do
  end do
 
  call write_fxv( f, x, v)

  call write_hdf5_data()

end program diagnostics_tester
