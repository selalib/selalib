
!***********************************************************************************
!
! Selalib      
! Module: unit_test.F90
!
!> @brief 
!> newton_raphson unit test
!   
!> @author                    
!> Aliou DIOUF (aliou.l.diouf@inria.fr)
!                                  
!************************************************************************************

program newton_raphson_tester 
#include "sll_working_precision.h"
  use OMP_LIB
  implicit none
  
  sll_real64 :: y
  sll_real64 :: tolx
  sll_real64 :: tolf
  sll_int32  :: i
  sll_int32  :: n_test
  sll_int32  :: nb_functions
  sll_int32  :: ok
  sll_int32  :: nb_cores 
  
  nb_cores = 16!OMP_GET_NUM_PROCS ()
  n_test = 1
  nb_functions = 7
  ok = 1 
  y = 0.d0
  tolx = 1.e-16
  tolf = epsilon(1.0_f64)

  if (n_test <= nb_cores) then
     !$OMP PARALLEL num_threads(n_test)	
     !$OMP DO SCHEDULE(STATIC,1)
     do i = 1, n_test
        call test_process(y, nb_functions, tolx, tolf, ok)
     enddo
     !$OMP END DO
     !$OMP END PARALLEL
  else
     !$OMP PARALLEL	
     !$OMP DO SCHEDULE(STATIC,(n_test)/nb_cores)
     do i = 1, n_test
        call test_process(y, nb_functions, tolx, tolf, ok)
     enddo
     !$OMP END DO
     !$OMP END PARALLEL
  endif
  
  if (ok == 1) then
     print*, 'PASS: Newton-Raphson unit test'
  endif

end program newton_raphson_tester

subroutine test_process(y, nb_functions, tolx, tolf, ok)
   use test_jacobian_and_f_module
   use newton_raphson
   
   sll_real64 :: y
   sll_int32  :: nb_functions
   sll_real64 :: tolx
   sll_real64 :: tolf
   sll_int32  :: ok
   sll_int32  :: i
   sll_real64 :: x
   sll_real64 :: error
   sll_real64 :: x0
   
   call random_number(x0)
   
   do i = 1, nb_functions        
      open (unit=10, file = 'function_num.dat')
      write (10,*) i
      close(10)        
      call newton_raphson_1D(y, f, jac, x, error, x0) 
      if (abs(x-exact_sol())>=tolx .and. error>=tolf) then
         ok = 0
         print*, 'Test stopped by failure'
         print*, 'x =', x, 'f =', f(x)
         stop
      endif
print*, 'x=',x
   enddo
   
 end subroutine test_process
