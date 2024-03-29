!**************************************************************
!
! Selalib 2012
! Module: test_sll_penta_diagonal.F90
!
!> @brief
!> Selalib sll_m_penta_diagonal.F90 tester
!
!> Last modification: September 20, 2012
!
!> @authors
!> Aliou DIOUF (aliou.l.diouf@inria.fr)
!
!**************************************************************

program test_sll_penta_diagonal

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

   use sll_m_penta_diagonal, only: &
      sll_o_create, &
      sll_o_delete, &
      sll_t_penta_diagonal_solver, &
      sll_o_solve

   implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   sll_int32                               :: n_max, nb_test
   sll_int32                               :: n, i_test, i, ierr
   sll_real64                              :: a, b, c
   sll_real64, dimension(:), allocatable   :: f
   sll_real64, dimension(:), allocatable   :: x, x_exact
   sll_real64                              :: error, norm
   type(sll_t_penta_diagonal_solver), pointer :: plan

   print *, ' '
   print *, 'Testing Toeplitz penta-diagonal solver...'
   print *, ' '

   n_max = 1000
   nb_test = 10

   do i_test = 1, nb_test

      do n = 3, n_max

         a = real(n*i_test, f64)
         b = a/4.
         c = a/4.-1

         SLL_ALLOCATE(f(n), ierr)
         SLL_ALLOCATE(x(n), ierr)
         SLL_ALLOCATE(x_exact(n), ierr)

         call random_number(x_exact)

         f(1) = a*x_exact(1) + b*x_exact(2) + c*x_exact(3)
         f(n) = c*x_exact(n - 2) + b*x_exact(n - 1) + a*x_exact(n)

         if (n > 3) then
            f(2) = b*x_exact(1) + a*x_exact(2) + b*x_exact(3) + c*x_exact(4)
            f(n - 1) = c*x_exact(n - 3) + b*x_exact(n - 2) + &
                       a*x_exact(n - 1) + b*x_exact(n)
         else
            f(2) = b*x_exact(1) + a*x_exact(2) + b*x_exact(3)
         end if

         do i = 3, n - 2
            f(i) = c*x_exact(i - 2) + b*x_exact(i - 1) + a*x_exact(i) + &
                   b*x_exact(i + 1) + c*x_exact(i + 2)
         end do

         plan => sll_o_create(n)
         call sll_o_solve(a, b, c, f, plan)
         x = plan%solution

         error = 0.0_f64
         norm = 0.0_f64
         do i = 1, n
            error = error + (x_exact(i) - x(i))**2
            norm = norm + (x_exact(i))**2
         end do

         error = error/norm; 
         print *, 'Test', i_test, ', nb_points =', n, ', error = ', error

         if (error > 1.0e-15_f64) then
            print *, 'Toeplitz penta-diagonal solver: FAILED'
            stop
         end if

         SLL_DEALLOCATE_ARRAY(f, ierr)
         SLL_DEALLOCATE_ARRAY(x, ierr)
         SLL_DEALLOCATE_ARRAY(x_exact, ierr)
         call sll_o_delete(plan); 
      end do
   end do

   print *, ' '
   print *, 'Toeplitz penta-diagonal solver: PASSED'
   print *, ' '

end program test_sll_penta_diagonal
