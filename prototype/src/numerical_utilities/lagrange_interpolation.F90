module lagrange_interpolation
#include "sll_working_precision.h"
#include "sll_assert.h"
  implicit none

contains

  ! lagrange_interpolate returns the value of y(x), using a Lagrange 
  ! interpolation based on the given array values of x(i) and y(x(i)).
  ! The interpolating polynomial is of degree 'degree'.
  function lagrange_interpolate( x, degree, xi, yi )
    sll_real64, intent(in)               :: x
    sll_int32, intent(in)                :: degree
    sll_real64, intent(in), dimension(:) :: xi
    sll_real64, intent(in), dimension(:) :: yi
    sll_real64                           :: lagrange_interpolate
    sll_real64, dimension(1:degree+1)    :: p
    sll_int32                            :: i  
    sll_int32                            :: deg
    sll_int32                            :: m  ! step size

    SLL_ASSERT( size(xi) >= degree+1 )
    SLL_ASSERT( size(yi) >= degree+1 )
    SLL_ASSERT( degree   >= 0 )

    if( (x < xi(1)) .or. (x > xi(degree+1)) ) then
       print *, 'lagrange_interpolate() warning: x is outside of the range ', &
            'xi given as argument.'
       print *, 'x = ', x, 'xmin = ', xi(1), 'xmax = ',xi(degree+1)
    end if

    ! Load the local array with the values of y(i), which are also the values
    ! of the Lagrange polynomials of order 0.
    do i=1,degree+1
       p(i) = yi(i)
    end do

    ! Build the Lagrange polynomials up to the desired degree. Here we use
    ! Neville's recursive relation, starting from the zeroth-order 
    ! polynomials and working towards the desired degree. 
    m = 1
    do deg=degree, 1, -1
       do i=1,deg
          p(i) = (1.0_f64/(xi(i)-xi(i+m)))*((x-xi(i+m))*p(i)+(xi(i)-x)*p(i+1))
       end do
       m = m+1
    end do
    lagrange_interpolate = p(1)
  end function lagrange_interpolate



end module lagrange_interpolation
