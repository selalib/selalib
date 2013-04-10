program test_sll_odd_degree_spline_interpolator_1d
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
  use sll_constants
  use util_constants

#ifndef STDF95
  use sll_module_interpolators_1d_base
#endif
  use sll_odd_degree_spline_interpolator_1d
    implicit none

#ifdef STDF95
  type(odd_degree_spline_1d_interpolator), pointer  :: interp
#else
  class(sll_interpolator_1d_base), pointer     :: interp
#endif

  type(odd_degree_spline_1d_interpolator), target :: spline

  sll_real64                            :: error, error_disp
  sll_real64                            :: phase
  sll_real64, allocatable, dimension(:) :: interpolation_points
  sll_real64, allocatable, dimension(:) :: data  
  sll_real64, allocatable, dimension(:) :: out, out_disp
  sll_real64, allocatable, dimension(:) :: coordinates_d, coordinates_disp
  sll_real64, allocatable, dimension(:) :: data_interp, data_interp_disp
  sll_int32                             :: ierr, i, degree, degree_max
  sll_int32, parameter                  :: n = 512
  sll_real64                            :: x_min, x_max, delta
  sll_real64                            :: mu, alpha

  SLL_ALLOCATE(data(n), ierr)
  SLL_ALLOCATE(out(n), ierr)
  SLL_ALLOCATE(out_disp(n), ierr)
  SLL_ALLOCATE(interpolation_points(n), ierr)
  SLL_ALLOCATE(coordinates_d(n), ierr)
  SLL_ALLOCATE(coordinates_disp(n), ierr)
  SLL_ALLOCATE(data_interp(n), ierr)
  SLL_ALLOCATE(data_interp_disp(n), ierr)



  print *, 'initialize data and interpolation_points array'

  degree_max = 11
  x_min = -10.0_f64
  x_max = 10.0_f64
  mu = (x_min+x_max)/2
  delta = (x_max - x_min ) / real(n-1,f64) 
  phase = 0.0_f64

  do degree=1,degree_max,2

  call random_number(alpha)
  alpha = alpha * delta

  do i=1,n
     coordinates_d(i) = x_min + (i-1)*delta 
     interpolation_points(i) = coordinates_d(i) 
     data(i)        = exp( - ( coordinates_d(i) - mu )**2  )
     data_interp(i) = exp( - ( interpolation_points(i) - mu )**2  )
  end do

  print*, 'odd_degree spline interpolation'
#ifdef STDF95
!  call odd_degree_spline_initialize(spline, n, x_min, x_max, 2, degree )
  call odd_degree_spline_initialize(spline, n, x_min, x_max, degree )
#else
!  call spline%initialize(n, x_min, x_max, 2, degree )
  call spline%initialize(n, x_min, x_max, degree )
#endif

  interp =>  spline
#ifdef STDF95
  out = odd_degree_spline_interpolate_array(interp, n, data, interpolation_points)
  out_disp = odd_degree_spline_interpolate_array_at_displacement(interp, n, &
                                                   data, coordinates_disp)
#else
  out = interp%interpolate_array(n, data, interpolation_points)
  out_disp = interp%interpolate_array_disp(n, data, alpha)
#endif

    if (alpha > 0 ) then 
       do i = 1, n
          coordinates_disp(i) = max(interpolation_points(i) - alpha, x_min)
          SLL_ASSERT((x_min <=coordinates_disp(i)).and.(coordinates_disp(i) <= x_max))
       end do
    else
       do i = 1, n
          coordinates_disp(i) = min(interpolation_points(i) - alpha, x_max)
          SLL_ASSERT((x_min <=coordinates_disp(i)).and.(coordinates_disp(i) <= x_max))
       end do
    endif

  do i=1,n
     data_interp_disp(i) = exp( - ( coordinates_disp(i) - mu )**2  )
  end do

  error = 0.0_f64
  error_disp = 0.0_f64

  do i=1,n   
     error = error + ( data_interp(i) - out(i) )**2
     error_disp = error_disp + ( data_interp_disp(i) - out_disp(i) )**2
     !print*, coordinates_disp(i), data_interp_disp(i) , out_disp(i)
  end do

  error = sqrt( error / sum( data_interp**2 ) )
  error_disp = sqrt( error_disp/sum( data_interp_disp**2 ) )

  print*, ''
  print*, '    splines degree =', degree
  print*, '    error =', error
  print*, '    error_disp =', error_disp
  print*, ''

  enddo

  SLL_DEALLOCATE_ARRAY(data, ierr)
  SLL_DEALLOCATE_ARRAY(out, ierr)
  SLL_DEALLOCATE_ARRAY(interpolation_points, ierr)
  SLL_DEALLOCATE_ARRAY(coordinates_d, ierr)
  SLL_DEALLOCATE_ARRAY(coordinates_disp, ierr)
  SLL_DEALLOCATE_ARRAY(data_interp, ierr)
  SLL_DEALLOCATE_ARRAY(data_interp_disp, ierr)

  print *, 'PASSED'

end program test_sll_odd_degree_spline_interpolator_1d
