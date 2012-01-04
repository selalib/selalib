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
           print*, 'Function not available'
           stop
      end select function_test

    end function fprime


end module test_func_module
