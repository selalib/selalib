!> @ingroup integration
!> @brief 
!> Rectangle integration
!> @details
!> Low-level mathematical utility 
!> that applies the 
!> Rectangle method to compute numeric integrals.
module sll_m_rectangle_integration
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  implicit none

  public :: &
    rectangle_integrate_1d

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifndef DOXYGEN_SHOULD_SKIP_THIS
abstract interface
   !> 1d real function
   function function_1d(x)
      use sll_m_working_precision ! can't pass a header file because the
                                ! preprocessor prevents double inclusion.
                                ! This is very rare.
      sll_real64             :: function_1d
      sll_real64, intent(in) :: x
   end function function_1d
end interface
#endif

!> Integrate numerically with Gauss-Lobatto formula
interface rectangle_integrate_1d
  module procedure rectangle_integral_1d 
end interface

contains


  !> @brief Integrate with rectangle formula
  !> @details To integrate the function \f$ f(x) \f$
  !> , we use the rectangle formula 
  !> \f[ \int_{-1}^1 f(x)dx \approx \sum_{k=1}^{n} w_k f(x_k) \f]
  !> where n and x represents the desired number of points and their 
  !> positions.
  !>
  !> @param f     the function to be integrated
  !> @param[in] x positions of points
  !> @param[in] n the number of points
  !> @return The value of the integral
  function rectangle_integral_1d( f, x, n )
    sll_real64                :: rectangle_integral_1d
    procedure(function_1d)    :: f
    sll_int32,  intent(in)    :: n 
    sll_real64, dimension(n)  :: x
    sll_int32                 :: k
    sll_real64                :: ans

    ans = 0.0_f64
    do k=1,n-1
       ans = ans + f(x(k))*(x(k+1)-x(k))
    end do
    rectangle_integral_1d = ans

  end function rectangle_integral_1d

  !> Returns a 1d array of size (n) containing rectangle 
  !> integration weights in the interval [x(1),x(n)].
  !> @param[in] n Number of gauss points.
  !> @param[in] x Point poisitions in interval.
  !> @return    w Array containing weights.
  function rectangle_weights( n, x ) result(w)
    sll_int32,  intent(in)           :: n 
    sll_real64, dimension(n)         :: x
    sll_real64, dimension(n)         :: w
    sll_int32                        :: k
    
    do k = 1, n-1
      w(k) = x(k+1) - x(k)
    end do
    w(n) = 0.0_f64

  end function rectangle_weights
  
end module sll_m_rectangle_integration
