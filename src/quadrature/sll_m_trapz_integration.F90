!> @ingroup integration
!> @brief 
!> Trapezoid formula for numerical integration
!> @details
!> Low-level mathematical utility 
!> that applies the 
!> Trapezoid formula to compute numeric integrals.
module sll_m_trapz_integration
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  implicit none

  public :: &
    trapz_integrate_1d

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
interface trapz_integrate_1d
  module procedure trapz_integral_1d 
end interface

contains


  !> @brief Integrate with trapz formula
  !> @details To integrate the function \f$ f(x) \f$
  !> , we use the trapz formula 
  !> \f[ \int_{-1}^1 f(x)dx \approx \sum_{k=1}^{n} w_k f(x_k) \f]
  !> where n and x represents the desired number of points and their 
  !> positions.
  !>
  !> @param f     the function to be integrated
  !> @param[in] x positions of points
  !> @param[in] n the number of points
  !> @return The value of the integral
  function trapz_integral_1d( f, x, n )
    sll_real64                :: trapz_integral_1d
    procedure(function_1d)    :: f
    sll_int32,  intent(in)    :: n 
    sll_real64, dimension(n)  :: x
    sll_int32                 :: k
    sll_real64                :: ans

    ans = 0.0_f64
    do k=1,n-1
       ans = ans + 0.5*(f(x(k))+f(x(k+1)))*(x(k+1)-x(k))
    end do
    trapz_integral_1d = ans

  end function trapz_integral_1d

  !> Returns a 1d array of size (n) containing trapz 
  !> integration weights in the interval [x(1),x(n)].
  !> @param[in] n Number of gauss points.
  !> @param[in] x Point poisitions in interval.
  !> @return    w Array containing weights.
  function trapz_weights( n, x ) result(w)
    sll_int32,  intent(in)           :: n 
    sll_real64, dimension(n)         :: x
    sll_real64, dimension(n)         :: w
    sll_int32                        :: k
    
    w(1) = 0.5_f64*(x(2)-x(1))
    do k = 2, n-1
      w(k) = 0.5_f64*(x(k+1)-x(k-1))
    enddo  
    w(n) = 0.5_f64*(x(n)-x(n-1))

  end function trapz_weights
  
end module sll_m_trapz_integration
