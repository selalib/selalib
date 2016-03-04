program cubic_spline_interpolator_1d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_periodic

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_cubic_spline_interpolator_1d, only: &
    sll_t_cubic_spline_interpolator_1d

  use sll_m_interpolators_1d_base, only: &
    sll_c_interpolator_1d

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    sll_int32, parameter  :: n1 = 20
    sll_int32, parameter  :: n2 = 64
    sll_int32, parameter  :: m = 512
    sll_real64, parameter :: tol1 = 1.d-4 ! tolerance for error (n1)
    sll_real64, parameter :: tol2 = 1.d-6 ! tolerance for error (n2)
    sll_real64            :: error1
    sll_real64            :: error2

    call test_for_given_n( n1, error1 )
    call test_for_given_n( n2, error2 )

    if (error1 > tol1) then
       print*, 'Failed for n=', n1
       print*, 'FAILED.'
    elseif (error2 > tol2) then
       print*, 'Failed for n=', n2
       print*, 'FAILED.'
    else
       print*, 'Successful, exiting program.'
       print*, 'PASSED'
    end if

contains

  subroutine test_for_given_n (n, error)
    sll_int32,  intent(in)  :: n
    sll_real64, intent(out) :: error

    class(sll_c_interpolator_1d), pointer       :: interp
    type(sll_t_cubic_spline_interpolator_1d), target :: spline

    sll_real64, allocatable, dimension(:) :: point
    sll_real64, allocatable, dimension(:) :: pdata  
    sll_real64, allocatable, dimension(:) :: fdata
    sll_real64, allocatable, dimension(:) :: coord
    sll_real64, allocatable, dimension(:) :: gdata
    
    sll_int32 :: ierr, i
    
    sll_real64  :: x_min, x_max, delta
    
    SLL_ALLOCATE(coord(n), ierr)
    SLL_ALLOCATE(pdata(n), ierr)
    SLL_ALLOCATE(point(m), ierr)
    SLL_ALLOCATE(fdata(m), ierr)
    SLL_ALLOCATE(gdata(m), ierr)
    
    print*, 'Initialize data and point array'
    x_min = 0.0_f64
    x_max = 2.0_f64 * sll_p_pi
    delta = (x_max - x_min ) / real(n-1,f64) 
    do i=1,n
       coord(i) = (i-1)*delta
       pdata(i) = f(coord(i))
       !print*, i,coord(i), pdata(i)
    end do
    
    delta = (x_max - x_min ) / real(m-1,f64) 
    do i=1,m
       point(i) = (i - 1) * delta
       gdata(i) = f(point(i))
    end do
    
    print*, 'Cubic spline interpolation'
    call spline%initialize(n, x_min, x_max, sll_p_periodic )
    
    interp => spline
    call interp%compute_interpolants(pdata)
    
    do i = 1, m
       fdata(i) = interp%interpolate_from_interpolant_value(point(i))
    end do
    
    error = maxval(abs(gdata-fdata))
    print*, 'error=', error
    

  end subroutine test_for_given_n

  function f(x)

    sll_real64 :: x
    sll_real64 :: f

    f = 2.0_f64*(sin(x) + 2.5_f64 + cos(x))

  end function f

end program cubic_spline_interpolator_1d
