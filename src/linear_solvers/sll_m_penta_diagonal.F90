!**************************************************************
!  Copyright INRIA
!  Authors :
!     CALVI project team
!
!  This code SeLaLib (for Semi-Lagrangian-Library)
!  is a parallel library for simulating the plasma turbulence
!  in a tokamak.
!
!  This software is governed by the CeCILL-B license
!  under French law and abiding by the rules of distribution
!  of free software.  You can  use, modify and redistribute
!  the software under the terms of the CeCILL-B license as
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info".
!**************************************************************

!> @ingroup utilities
!> @brief Toeplitz penta-diagonal system solver
!> @details
!> This module are not fully tested, please use it carefully.
module sll_m_penta_diagonal
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

   implicit none

   public :: &
      sll_o_create, &
      sll_o_delete, &
      sll_t_penta_diagonal_solver, &
      sll_o_solve

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   !> Initialize the penta diagonal solver
   type sll_t_penta_diagonal_solver
      sll_int32                         :: n
      sll_real64, dimension(:), pointer :: e1
      sll_real64, dimension(:), pointer :: e2
      sll_real64, dimension(:), pointer :: y
      sll_real64, dimension(:), pointer :: z
      sll_real64, dimension(:), pointer :: w
      sll_real64, dimension(:), pointer :: solution
   end type sll_t_penta_diagonal_solver

   interface sll_o_create
      module procedure new_penta_diagonal
   end interface sll_o_create

   interface sll_o_solve
      module procedure solve_penta_diagonal
   end interface sll_o_solve

   interface sll_o_delete
      module procedure delete_penta_diagonal
   end interface sll_o_delete

contains

   !> allocate the solveR
   function new_penta_diagonal(n) result(plan)

      sll_int32                               :: n, ierr
      type(sll_t_penta_diagonal_solver), pointer :: plan

      if (n < 3) then
         print *, 'Matrix size must be at least 3x3'
         print *, 'Exiting...'
         stop
      end if

      ! Plan allocation
      SLL_ALLOCATE(plan, ierr)

      ! Plan components allocation
      plan%n = n; 
      SLL_ALLOCATE(plan%e1(n), ierr)
      SLL_ALLOCATE(plan%e2(n), ierr)
      SLL_ALLOCATE(plan%y(n), ierr)
      SLL_ALLOCATE(plan%z(n), ierr)
      SLL_ALLOCATE(plan%w(n), ierr)
      SLL_ALLOCATE(plan%solution(n), ierr)

   end function new_penta_diagonal

   !> The solution will be set in plan%solution
   subroutine solve_penta_diagonal(a, b, c, f, plan)

      sll_real64                              :: a, b, c
      sll_real64, dimension(:)                :: f
      type(sll_t_penta_diagonal_solver), pointer :: plan
      sll_real64                              :: s, t, p, l1, l2
      sll_real64                              :: d, d1, d2
      sll_int32                               :: n, i
      sll_real64                              :: sign_of_a = 0.0_f64

      if (abs(a) <= 2.0_f64*(abs(b) + abs(c))) then
         print *, 'a, b, and c must be such that: |a| > 2(|b|+|c|)'
         print *, a, b, c
         print *, 'Exiting...'
         stop
      end if

      n = plan%n

      s = (a/2 + c)*(a/2 + c) - b*b
      t = a*a/2 - b*b - 2*c*c

      if (a /= 0.0_f64) then
         sign_of_a = a/abs(a)
      end if

      p = (a - 2*c)/4 + sign_of_a*sqrt(sign_of_a*(a - 2*c)* &
                                       sqrt(s) + t)/2.0_f64 + sign_of_a*sqrt(s)/2.0_f64
      l1 = b/(p + c)
      l2 = c/p

      do i = 1, n
         plan%y(i) = f(i)/p
         plan%e1(i) = 0.0_f64
         plan%e2(i) = 0.0_f64
      end do
      plan%e1(1) = 1.0_f64
      plan%e2(2) = 1.0_f64

      plan%y = solve_subsystem(l1, l2, plan%y, n)
      plan%z = solve_subsystem(l1, l2, plan%e1, n)
      plan%w = solve_subsystem(l1, l2, plan%e2, n)

      d = 1.0_f64 + l1*l2*(plan%z(2) + plan%w(1)) + l2*l2*(plan%z(1) + plan%w(2)) + &
          l1*l1*plan%z(1) + l2**4*(plan%z(1)*plan%w(2) - plan%z(2)*plan%w(1))

      d1 = l2**4*(plan%y(1)*plan%w(2) - plan%y(2)*plan%w(1)) + &
           (l1*l1 + l2*l2)*plan%y(1) + l1*l2*plan%y(2)

      d2 = l2**4*(plan%y(2)*plan%z(1) - plan%y(1)*plan%z(2)) + &
           l1*l2*plan%y(1) + l2*l2*plan%y(2)

      do i = 1, n
         plan%solution(i) = plan%y(i) - (d1*plan%z(i) + d2*plan%w(i))/d
      end do

   end subroutine solve_penta_diagonal

   function solve_subsystem(l1, l2, b, n) result(x)

      sll_real64               :: l1, l2
      sll_real64, dimension(:) :: b
      sll_int32                :: n, i
      sll_real64, dimension(n) :: x, y
      sll_real64 :: tmp

      y(1) = b(1)
      y(2) = b(2) - l1*y(1)

      do i = 3, n
         tmp = b(i) - (l2*y(i - 2) + l1*y(i - 1))
         if (abs(tmp) < 1e-30_f64) then
            tmp = 0.0_f64
         end if
         y(i) = tmp
      end do

      x(n) = y(n)
      x(n - 1) = y(n - 1) - l1*y(n)
      do i = n - 2, 1, -1
         tmp = y(i) - (l1*x(i + 1) + l2*x(i + 2))
         if (abs(tmp) < 1e-30_f64) then
            tmp = 0.0_f64
         end if
         x(i) = tmp
      end do

   end function solve_subsystem

   !> deallocat the solver
   subroutine delete_penta_diagonal(plan)

      type(sll_t_penta_diagonal_solver), pointer :: plan
      sll_int32                               :: ierr

      ! Plan components deallocation
      SLL_DEALLOCATE_ARRAY(plan%e1, ierr)
      SLL_DEALLOCATE_ARRAY(plan%e2, ierr)
      SLL_DEALLOCATE_ARRAY(plan%y, ierr)
      SLL_DEALLOCATE_ARRAY(plan%z, ierr)
      SLL_DEALLOCATE_ARRAY(plan%w, ierr)
      SLL_DEALLOCATE_ARRAY(plan%solution, ierr)

      ! Plan deallocation
      SLL_DEALLOCATE_ARRAY(plan, ierr)

   end subroutine delete_penta_diagonal

end module sll_m_penta_diagonal

