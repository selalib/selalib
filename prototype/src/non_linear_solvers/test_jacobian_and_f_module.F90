module test_jacobian_and_f_module
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
           f = x**2
        case(2)
           f = exp(x) - exp(2.d0)
        case(3)
           f = x - 1
        case(4)
           f = log(x**2 + 1)
        case(5)
           f = 2*x - 1
        case(6)
           f = log(exp(x)+.5)
        case(7)
           f = x**16
        case default
           print*, 'Function not available'
           stop
      end select function_test

    end function f


    function jac(x) 
     
      sll_real64, intent(in) :: x
      sll_real64 :: jac
      sll_int32 :: function_num

      open (unit=10, file='function_num.dat')
         read (10,*) function_num
      close(10)

      function_test: select case(function_num)
        case(1)
           jac = 2*x
        case(2)
           jac = exp(x)
        case(3)
           jac = 1.d0
        case(4)
           jac = 2*x/(x**2 + 1)
        case(5)
           jac = 2.d0
        case(6)
           jac = exp(x)/(exp(x)+.5)
        case(7)
           jac = 16*x**15
        case default
           print*, 'Function not available'
           stop
      end select function_test

    end function jac

    function exact_sol()   
   
      sll_real64 :: exact_sol
      sll_int32 :: function_num

      open (unit=10, file='function_num.dat')
         read (10,*) function_num
      close(10)

      function_test: select case(function_num)
        case(1)
           exact_sol = 0.d0
        case(2)
           exact_sol = 2.d0
        case(3)
           exact_sol = 1.d0
        case(4)
           exact_sol = 0.d0
        case(5)
           exact_sol = .5d0
        case(6)
           exact_sol = -log(2.d0)
        case(7)
           exact_sol = 0.d0
        case default
           print*, 'Function not available'
           stop
      end select function_test

    end function exact_sol

end module test_jacobian_and_f_module
