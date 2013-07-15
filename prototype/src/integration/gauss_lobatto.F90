module gauss_lobatto_integration
#include "sll_working_precision.h"
#include "sll_assert.h"
  
implicit none

abstract interface
   function function_1D(x)
      use sll_working_precision ! can't pass a header file because the
                                ! preprocessor prevents double inclusion.
                                ! This is very rare.
      sll_real64             :: function_1D
      sll_real64, intent(in) :: x
   end function function_1D
end interface

interface gauss_lobatto_integrate_1d
  module procedure gauss_lobatto_integral_1d 
end interface

contains


  !> @brief Gauss-Lobatto Quadrature.
  !> @details To integrate the function \f$ f(x) \f$
  !> (real-valued and of a single, real-valued argument x)
  !> over the interval \f$ [a,b] \f$, we use the Gauss-Lobatto formula 
  !> \f[ \int_{-1}^1 f(x)dx \approx \sum_{k=1}^{n} w_k f(x_k) \f]
  !> where n represents the desired number of Gauss points.
  !>
  !> the function maps the interval \f$ [-1,1] \f$ into the
  !> arbitrary interval \f$ [a,b] \f$.
  !>
  !> To be considered is to split this function into degree-specific
  !> functions to avoid the select statement.
  !> @param f the function to be integrated
  !> @param[in] a left-bound of the definition interval of f  
  !> @param[in] b right-bound of the definition interval of f 
  !> @param[in] n the desired number of Gauss points
  !> @return The value of the integral
  function gauss_lobatto_integral_1D( f, a, b, n )
    sll_real64                :: gauss_lobatto_integral_1D
    procedure(function_1D)    :: f
    sll_real64, intent(in)    :: a
    sll_real64, intent(in)    :: b
    sll_int32,  intent(in)    :: n 
    sll_real64, dimension(n)  :: xk
    sll_real64, dimension(n)  :: wk
    sll_int32                 :: k
    sll_int32                 :: err
    sll_real64                :: alpha(0:n-1), beta(0:n-1)
    sll_real64                :: de(n), da(n), db(n)
    sll_real64                :: ans
    sll_real64                :: x
    sll_real64                :: c1
    sll_real64                :: c2

    xk(:) = 0.0_f64
    wk(:) = 0.0_f64

    alpha = 0.0_f64
    do k = 0, n-1
       beta(k)=real(k,kind(n-1))**2/((2.0d0*k+1)*(2.0d0*k-1))
    end do

    !for Gauss-Legendre and Gauss-Lobatto, beta(0)=int(dlambda)
    !see Algorithm xxx - ORTHPOL: A package of routines for  generating orthogonal
    !polynomials and Gauss-type quadrature rules by _Walter Gautschi_
    beta(0)=2.0d0

    call dlob(n-2,alpha,beta,-1.0_f64,1._f64,xk,wk,err,de,da,db)

    ans = 0.0
    ! need to map the interval [-1,1] into the interval [a,b]
    c1 = 0.5_f64*(b-a)
    c2 = 0.5_f64*(b+a)
    do k=1,n
       x = c1*xk(k) + c2
       ans = ans + f(x)*wk(k)
    end do
    gauss_lobatto_integral_1D = c1*ans

  end function gauss_lobatto_integral_1D



  subroutine test_gauss_lobatto( f, a, b, n )
    procedure(function_1D)    :: f
    sll_real64, intent(in)    :: a
    sll_real64, intent(in)    :: b
    sll_int32,  intent(in)    :: n 
    sll_real64, dimension(n)  :: xk
    sll_real64, dimension(n)  :: wk
    sll_int32                 :: k
    sll_int32                 :: err
    sll_real64                :: alpha(0:n-1), beta(0:n-1)
    sll_real64                :: de(n), da(n), db(n)

    xk(:) = 0.0_f64
    wk(:) = 0.0_f64

    alpha = 0.0_f64
    do k = 0, n-1
       beta(k)=real(k,kind(n-1))**2/((2.0d0*k+1)*(2.0d0*k-1))
    end do

    !for Gauss-Legendre and Gauss-Lobatto, beta(0)=int(dlambda)
    !see Algorithm xxx - ORTHPOL: A package of routines for  generating orthogonal
    !polynomials and Gauss-type quadrature rules by _Walter Gautschi_
    beta(0)=2.0d0

    call dlob(n-2,alpha,beta,a,b,xk,wk,err,de,da,db)


    write(*,*) "xk=",xk(:)
    write(*,*) "wk=",wk(:)


  end subroutine test_gauss_lobatto


end module gauss_lobatto_integration
