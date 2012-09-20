!**************************************************************
!
! Selalib 2012     
! Module: sll_quintic_splines.F90
!
!> @brief 
!> Selalib quintic splines interpolator
!
!> Last modification: September 20, 2012
!   
!> @authors                    
!> Aliou DIOUF (aliou.l.diouf@inria.fr)
!                                  
!**************************************************************

module sll_quintic_splines
#include "sll_working_precision.h"
#include "sll_memory.h"
use sll_toep_penta_diagonal
use arbitrary_degree_splines
implicit none

  type quintic_splines_plan
    sll_int32                             :: n
    sll_real64                            :: xmin, xmax
    sll_real64, dimension(:), allocatable :: b_at_node
    sll_real64, dimension(:), allocatable :: coeffs
  end type quintic_splines_plan

contains

  function new_quintic_splines(n, xmin, xmax, f) result(plan)

    sll_int32                         :: n, ierr
    sll_real64                        :: xmin, xmax
    sll_real64, dimension(:)          :: f
    type(quintic_splines_plan),pointer :: plan

    ! Plan allocation
    SLL_ALLOCATE(plan, ierr)
    ! plan component allocation
    SLL_ALLOCATE(plan%coeffs(n+6), ierr)

    plan%n = n
    plan%xmin = xmin
    plan%xmax = xmax
    plan%b_at_node = uniform_b_splines_at_x( 5, 0.d0 )
    call compute_coeffs(f, plan)

  end function new_quintic_splines


  subroutine compute_coeffs(f, plan_splines)

  ! f is the vector of the values of the function 
  !  in the nodes of the mesh*/

    sll_real64, dimension(:)                :: f
    type(quintic_splines_plan), pointer     :: plan_splines
    sll_real64                              :: a, b, c
    sll_real64                              :: xmin, xmax, h
    sll_int32                               :: n
    type(toep_penta_diagonal_plan), pointer :: plan_pent

    n = plan_splines%n
    xmin = plan_splines%xmin
    xmax = plan_splines%xmax
    h = (xmax-xmin)/n

    a = plan_splines%b_at_node(3)
    b = plan_splines%b_at_node(2)
    c = plan_splines%b_at_node(1)

    plan_pent => new_toep_penta_diagonal(n+6)
    plan_splines%coeffs = solve_toep_penta_diagonal(a, b, c, f, plan_pent)
    call delete_toep_penta_diagonal(plan_pent)

  end subroutine compute_coeffs

  function quintic_splines(x, plan_splines) result(s)
  ! The interpolator spline function

    type(quintic_splines_plan), pointer :: plan_splines
    sll_int32                          :: n, left, j
    sll_real64                         :: x, xmin, xmax
    sll_real64                         :: h, s, t0
    sll_real64, dimension(6)           :: b

    xmin = plan_splines%xmin
    xmax = plan_splines%xmax
    n = plan_splines%n
    h = (xmax-xmin)/n

    t0 = (x-xmin)/h
    left = int(t0) ! Determine the leftmost support index 'i' of x
    t0 = t0 - left ! compute normalized_offset

    b = uniform_b_splines_at_x( 5, t0 )
    s = 0

    do j=left-5,left
      if( (j>=-5) .and. (j<=n) ) then
        s = s + plan_splines%coeffs(j+6) * b(j-left+6)
      endif
    enddo

  end function quintic_splines


  subroutine delete_quintic_splines(plan)

    type(quintic_splines_plan), pointer :: plan
    sll_int32                          :: ierr

    SLL_DEALLOCATE_ARRAY(plan%coeffs, ierr)
    SLL_DEALLOCATE_ARRAY(plan, ierr)
 
  end subroutine delete_quintic_splines

end module sll_quintic_splines
