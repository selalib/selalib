
!***********************************************************************************
!
! Selalib      
! Module: newton_raphson.F90  
!   
!> @author                    
!> Aliou DIOUF (aliou.l.diouf@inria.fr)
!                                  
!> @brief 
!> Solves a non linear system of equations by using Newton-Raphson algorithm.
!
!> @details
!> Suppose you have to solve the following system: \n
!> \f$ \f$\n
!> \f$ \left \{ \begin{array}{ccccccccccccccccc}
!> f_1(x_1,x_2,\cdots,x_n) = y_1 \\\
!> f_2(x_1,x_2,\cdots,x_n) = y_2 \\\
!> \vdots \\\
!> f_n(x_1,x_2,\cdots,x_n) = y_n \\\
!> \end{array}
!> \right.
!> \f$\n
!> \f$ \f$\n
!> To solve it, the principle of Newton-Raphson method is:\n
!> - Pose:\n 
!> \f$ x=(x_1,x_2,\cdots,x_n) \f$ \n 
!> \f$ f_i-y_i=F_i \f$ \n 
!> \f$ F=(F_1,F_2,\cdots,F_n)\f$ \n
!> - We have:\n
!> \f$F_i(x+\delta x) = F_i (x) + \displaystyle\sum_{j=1}^n\frac{\partial F_i}
!> {\partial x_j}\delta x_j + O(\delta x^2) \f$\n
!> The matrix of partial derivatives appearing in this equation is the Jacobian
!> matrix J:\n
!> \f$ J_{ij}=\frac{\partial F_i}{\partial x_j}\f$\n
!> - By neglecting terms of order \f$\delta x^2 \f$ and higher and by setting 
!> \f$ F(x + \delta x) = 0\f$, we obtain a set of linear equations for the corrections
!> \f$\delta x\f$ that move each function closer to zero simultaneously, namely\n
!> \f$ J\times\delta x=-F \f$ \n
!> This linear system yields the corrections \f$\delta x\f$ which is then added to the 
!> solution vector,\n
!> \f$ x_{new} = x_{old} + \delta x\f$\n
!> The process is iterated to convergence to the final solution x of the non linear
!> system.
!
!**************************************************************************************

module newton_raphson
#include "sll_working_precision.h"
implicit none

#ifdef STDF95
#else
  abstract interface
     function function_1D(x)
       use sll_working_precision ! can't pass a header file because the
                                 ! preprocessor prevents double inclusion.
                                 ! This is very rare.
       sll_real64             :: function_1D
       sll_real64, intent(in) :: x
     end function function_1D
  end interface
#endif

  interface newton_raphson_1D
     module procedure newton_raphson_1D_real
  end interface

contains

  !> @brief
  !> Solves a non linear equation by using Newton-Raphson algorithm.\n

  !> @details
  !> Here we have to solve: \n
  !> \f$ f(x) = y \f$
  
  !> @param[in] x0         : initial (guessed) solution \n
  !> @param[in] f, y       : f - y is the zeroed function. f is an abstract function \n
  !> @param[in] jac        : the jacobian of f (f'(x) in 1D). jac is an abstract function \n
  !> @param[in] tolx, tolf : the error telerance for x and f respectively \n
  !> @param[out] x         : such f(x) - y = 0 \n
  !> @param[out] error     : |f(x) - y|

  subroutine newton_raphson_1D_real(y, f, jac, x, error, x0, tolx, tolf)
  
    sll_real64, intent(in)            :: y
#ifdef STDF95
    sll_real64                        :: f, jac
#else
    procedure(function_1D)            :: f
    procedure(function_1D)            :: jac
#endif
    sll_real64, intent(out)           :: x
    sll_real64, intent(out)           :: error
    sll_real64, optional              :: x0
    sll_real64, optional              :: tolx
    sll_real64, optional              :: tolf
    sll_real64                        :: x0_local
    sll_real64                        :: tolx_local
    sll_real64                        :: tolf_local
    sll_real64                        :: dx


    ! Default values of optional variables x0, tolx, tolf
    if (present(x0)) then
       x0_local = x0
    else
       x0_local= 1.d0
    endif
    if (present(tolx)) then
       tolx_local = tolx
    else
       tolx_local = 1.e-16
    endif
    if (present(tolx)) then
       tolf_local = tolf
    else
       tolf_local = epsilon(1.0_f64)
    endif
  
    ! Initialization
    x = x0_local
    if (jac(x)==0.d0) then
       print*, '(Bad initial (guessed) solution. Please choose another:'
       stop
    endif
    dx = 1.d0

    convergence: do while ( abs(f(x))>=tolf_local .and. abs(dx)>=tolx_local )
       dx = -(f(x) - y) / jac(x)
       x = x + dx     
    enddo convergence

    error = abs(f(x)-y)
  
  end subroutine newton_raphson_1D_real


end module newton_raphson
