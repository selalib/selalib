
!****************************************************
!
! Selalib      
! Module: test_func_module.F90
!
!> @brief 
!> Module of test functions for sll_splines unit test
!   
!> @authors                    
!> Aliou DIOUF (aliou.l.diouf@inria.fr)
!                                  
!****************************************************

module test_func_module
#include "sll_working_precision.h"
  use util_constants

! Test functions modules

  implicit none
  contains


    function f(x, function_num )
     
      sll_real64, intent(in) :: x
      sll_real64             :: f
      sll_int32              :: function_num

      function_test: select case(function_num)
        case(1)
           f = cos(x)
        case(2)
           f = sin(x)
        case(3)
           f = cos(x) + sin(x)
        case(4)
           f = cos(x) * sin(x)
        case(5)
           f = exp( cos(x) )
        case(6)
           f = exp( sin(x) )
        case(7)
           f = exp(cos(x)) + exp( sin(x) )
        case(8)
           f = exp(cos(x)) * exp(sin(x))
        case(9)
           f = cos(sin(x))
        case(10)
           f = sin(cos(x))
        case(11)
           f = cos(x)**2
        case(12)
           f = sin(x)**2
        case default
           if (x==real(function_num-13,f64)*(XMAX-XMIN)/real(NP-1,f64)+XMIN .or. &
               ! real(i-1-12,f64)*(XMAX-XMIN)/real(NP-1,f64)+XMIN) then
               x==real(function_num-13+NP,f64)*(XMAX-XMIN)/real(NP-1,f64)+XMIN) then ! for periodic BC
              f = 1.d0
           else
              f = 0.d0
           endif
      end select function_test

    end function f


    function fprime(x, function_num ) 
     
      sll_real64, intent(in) :: x
      sll_real64             :: fprime
      sll_int32              :: function_num

      function_test: select case(function_num)
        case(1)
           fprime = -sin(x)
        case(2)
           fprime = cos(x)
        case(3)
           fprime = cos(x) - sin(x)
        case(4)
           fprime = cos(2*x)
        case(5)
           fprime = -sin(x) * exp(cos(x))
        case(6)
           fprime = cos(x) * exp(sin(x))
        case(7)
           fprime = cos(x)*exp(sin(x)) - sin(x)*exp(cos(x))
        case(8)
           fprime = (cos(x)-sin(x))*exp(cos(x)+sin(x))
        case(9)
           fprime = -cos(x)*sin(sin(x))
        case(10)
           fprime = -sin(x)*cos(cos(x))
        case(11)
           fprime = -sin(2*x)
        case(12)
           fprime = sin(2*x)
        case default
           fprime = 0.d0
      end select function_test

    end function fprime


end module test_func_module
