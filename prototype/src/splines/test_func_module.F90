module test_func_module
#include "sll_working_precision.h"

! Test functions modules

  implicit none
  contains


    function f(x)
     
      sll_real64, intent(in) :: x
      sll_real64 ::  f
      sll_int32 :: function_num

      open (unit=10, file='function_num.dat')
         read (10,*) function_num
      close(10)

      function_test: select case(function_num)
        case(1)
           f = x
        case(2)
           f = x**2
        case(3)
           f = exp(x)
        case(4)
           f = cos(x)
        case(5)
           f = sin(x)
        case(6)
           f = exp(cos(x))
        case(7)
           f = exp(sin(x))
        case(8)
           f = cos(exp(x))
        case(9)
           f = sin(exp(x))
        case default
           print*, 'Function not available'
           stop
      end select function_test

    end function f


    function fprime(x) 
     
      sll_real64, intent(in) :: x
      sll_real64 :: fprime
      sll_int32 :: function_num

      open (unit=10, file='function_num.dat')
         read (10,*) function_num
      close(10)

      function_test: select case(function_num)
        case(1)
           fprime = 1.0_f64
        case(2)
           fprime = 2*x
        case(3)
           fprime = exp(x)
        case(4)
           fprime = -sin(x)
        case(5)
           fprime = cos(x)
        case(6)
           fprime = -sin(x)*exp(cos(x))
        case(7)
           fprime = cos(x)*exp(sin(x))
        case(8)
           fprime = -sin(exp(x))*exp(x)
        case(9)
           fprime = cos(exp(x))*exp(x)
        case default
           print*, 'Function not available'
           stop
      end select function_test

    end function fprime


end module test_func_module
