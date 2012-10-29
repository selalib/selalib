!**************************************************************
!
! Selalib 2012     
! Module: sll_toep_penta_diagonal.F90
!
!> @brief 
!> Selalib Toeplitz penta-diagonal system solver
!
!> Last modification: September 20, 2012
!   
!> @authors                    
!> Aliou DIOUF (aliou.l.diouf@inria.fr)
!                                  
!**************************************************************

module sll_toep_penta_diagonal
#include "sll_working_precision.h"
#include "sll_memory.h"
implicit none

  type toep_penta_diagonal_plan  
    sll_int32                             :: n
    sll_real64, dimension(:), pointer :: e1
    sll_real64, dimension(:), pointer :: e2
    sll_real64, dimension(:), pointer :: y
    sll_real64, dimension(:), pointer :: z
    sll_real64, dimension(:), pointer :: w
    sll_real64, dimension(:), pointer :: x
    sll_real64, dimension(:), pointer :: for_subsys
  end type toep_penta_diagonal_plan

contains 

  function new_toep_penta_diagonal(n) result(plan)

    sll_int32                               :: n, ierr
    type(toep_penta_diagonal_plan), pointer :: plan

    if (n<3) then
      print*, 'Matrix size must be at least 3x3'
      print*, 'Exiting...'
      stop
    endif

    ! Plan allocation
    SLL_ALLOCATE(plan, ierr)

    ! Plan components allocation
    plan%n  = n; 
    SLL_ALLOCATE(plan%e1(n), ierr)
    SLL_ALLOCATE(plan%e2(n), ierr)
    SLL_ALLOCATE(plan%y(n), ierr)
    SLL_ALLOCATE(plan%x(n), ierr)

  end function new_toep_penta_diagonal


  function solve_toep_penta_diagonal(a, b, c, f, plan) result(x)

    sll_real64                              :: a, b, c
    sll_real64, dimension(:)                :: f
    type(toep_penta_diagonal_plan), pointer :: plan
    sll_real64                              :: s, t, p, l1, l2
    sll_real64                              :: d, d1, d2
    sll_real64, dimension(plan%n)           :: e1, e2, x, y, z, w
    sll_int32                               :: n, i, sign_of_a=0

    if ( abs(a) <= 2*(abs(b)+abs(c)) ) then   
      print*, 'a, b, and c must be such that: |a| > 2(|b|+|c|)'
      print*, a, b, c
      print*, 'Exiting...'
      stop
    endif

    e1  = plan%e1
    e2  = plan%e2
    y   = plan%y
    x   = plan%x
    n   = plan%n

    s = (a/2+c)*(a/2+c) - b*b
    t = a*a/2 - b*b - 2*c*c

    if (a/=0.d0) then
       sign_of_a = int(a/abs(a))
    endif

    p = (a-2*c)/4 + sign_of_a*sqrt(sign_of_a*(a-2*c)* &
                     sqrt(s)+t)/2 + sign_of_a*sqrt(s)/2
    l1 = b/(p+c)
    l2 = c/p

    do i=1,n
      y(i) = f(i)/p
      e1(i) = 0.d0
      e2(i) = 0.d0
    enddo
    e1(1) = 1.d0
    e2(2) = 1.d0

    y = solve_subsystem(l1, l2, y, n)
    z = solve_subsystem(l1, l2, e1,  n)
    w = solve_subsystem(l1, l2, e2,  n)

    d = 1.d0 + l1*l2*(z(2)+w(1)) + l2*l2*(z(1)+w(2)) + &
                l1*l1*z(1) + l2**4*(z(1)*w(2)-z(2)*w(1))

    d1 = l2**4*(y(1)*w(2)-y(2)*w(1)) + &
         (l1*l1+l2*l2)*y(1) + l1*l2*y(2)

    d2 = l2**4*(y(2)*z(1)-y(1)*z(2)) + &
                 l1*l2*y(1) + l2*l2*y(2)

    do i=1,n
      x(i) = y(i) - (d1*z(i)+d2*w(i))/d   
    enddo

  end function solve_toep_penta_diagonal


  function solve_subsystem(l1, l2, b, n) result (x)

    sll_real64               :: l1, l2
    sll_real64, dimension(:) :: b
    sll_int32                :: n, i
    sll_real64, dimension(n) :: x, y
     
    y(1) = b(1)
    y(2) = b(2) - l1*y(1)
    do i=3,n
      y(i) = b(i) - ( l2*y(i-2) + l1*y(i-1) )
    enddo

    x(n) = y(n)
    x(n-1) = y(n-1) - l1*y(n)
    do i=n-2,1,-1
      x(i) = y(i) - ( l1*x(i+1) + l2*x(i+2) )
    enddo

  end function solve_subsystem


  subroutine delete_toep_penta_diagonal(plan)

    type(toep_penta_diagonal_plan), pointer :: plan
    sll_int32                               :: ierr

    ! Plan components deallocation 
    SLL_DEALLOCATE_ARRAY(plan%e1, ierr)
    SLL_DEALLOCATE_ARRAY(plan%e2, ierr)
    SLL_DEALLOCATE_ARRAY(plan%y, ierr)

    ! Plan deallocation
    SLL_DEALLOCATE_ARRAY(plan, ierr)
 
  end subroutine delete_toep_penta_diagonal

end module sll_toep_penta_diagonal

