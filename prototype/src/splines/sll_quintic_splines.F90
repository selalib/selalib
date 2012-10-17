!**************************************************************
!
! Selalib 2012     
! Module: sll_quintic_splines.F90
!
!> @brief 
!> Selalib quintic splines interpolator
!
!> Last modification: October 03, 2012
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

  type quintic_splines_plan_uniform
    sll_int32                             :: nb_cells
    sll_real64                            :: xmin
    sll_real64                            :: xmax
    sll_real64, dimension(:), allocatable :: b_at_node
    sll_real64, dimension(:), allocatable :: coeffs
  end type quintic_splines_plan_uniform

  type quintic_splines_plan_non_uni
    sll_real64, dimension(:), allocatable     :: coeffs
    type(arbitrary_degree_spline_1d), pointer :: spline_obj
  end type quintic_splines_plan_non_uni  

  interface compute_coeffs
     module procedure compute_coeffs_uniform, compute_coeffs_non_uni
  end interface compute_coeffs

  interface quintic_splines
     module procedure quintic_splines_uniform, quintic_splines_non_uni
  end interface quintic_splines

  interface delete_quintic_splines
     module procedure delete_quintic_splines_uniform, delete_quintic_splines_non_uni
  end interface delete_quintic_splines

contains

  ! *************************************************************************
  !
  !                    UNIFORM QUINTIC SPLINES STUFFS
  !
  ! *************************************************************************

  function new_quintic_splines_uniform(nb_cells, xmin, xmax, f) result(plan)

    sll_int32                                   :: nb_cells, ierr
    sll_real64                                  :: xmin, xmax
    sll_real64, dimension(:)                    :: f
    type(quintic_splines_plan_uniform), pointer :: plan

    ! Plan allocation
    SLL_ALLOCATE(plan, ierr)
    ! plan component allocation
    SLL_ALLOCATE(plan%coeffs(nb_cells+6), ierr)

    plan%nb_cells = nb_cells
    plan%xmin = xmin
    plan%xmax = xmax
    plan%b_at_node = uniform_b_splines_at_x( 5, 0.d0 )
    call compute_coeffs_uniform(f, plan)

  end function new_quintic_splines_uniform


  subroutine compute_coeffs_uniform(f, plan_splines)

  ! f is the vector of the values of the function 
  !  in the nodes of the mesh*/

    sll_real64, dimension(:)                        :: f
    type(quintic_splines_plan_uniform), pointer     :: plan_splines
    sll_real64                                      :: a, b, c
    sll_int32                                       :: nb_cells
    type(toep_penta_diagonal_plan), pointer         :: plan_pent

    nb_cells = plan_splines%nb_cells
    a = plan_splines%b_at_node(3)
    b = plan_splines%b_at_node(2)
    c = plan_splines%b_at_node(1)

    plan_pent => new_toep_penta_diagonal(nb_cells+6)
    plan_splines%coeffs = solve_toep_penta_diagonal(a, b, c, f, plan_pent)
    call delete_toep_penta_diagonal(plan_pent)

  end subroutine compute_coeffs_uniform

  function quintic_splines_uniform(x, plan_splines) result(s)
  ! The interpolator spline function

    type(quintic_splines_plan_uniform), pointer :: plan_splines
    sll_int32                                   :: nb_cells, left, j
    sll_real64                                  :: x, xmin, xmax
    sll_real64                                  :: h, s, t0
    sll_real64, dimension(6)                    :: b

    xmin = plan_splines%xmin
    xmax = plan_splines%xmax
    nb_cells = plan_splines%nb_cells
    h = (xmax-xmin)/nb_cells

    t0 = (x-xmin)/h
    left = int(t0) ! Determine the leftmost support index 'i' of x
    t0 = t0 - left ! compute normalized_offset

    b = uniform_b_splines_at_x( 5, t0 )
    s = 0

    do j=left-5,left
      if( (j>=-5) .and. (j<=nb_cells) ) then
        s = s + plan_splines%coeffs(j+6) * b(j-left+6)
      endif
    enddo

  end function quintic_splines_uniform


  subroutine delete_quintic_splines_uniform(plan)

    type(quintic_splines_plan_uniform), pointer :: plan
    sll_int32                                   :: ierr

    SLL_DEALLOCATE_ARRAY(plan%coeffs, ierr)
    SLL_DEALLOCATE_ARRAY(plan, ierr)
 
  end subroutine delete_quintic_splines_uniform

  ! *************************************************************************
  !
  !                  NON UNIFORM QUINTIC SPLINES STUFFS
  !
  ! *************************************************************************

  ! nb_pts = n + 1
  function new_quintic_splines_non_uni(knots, f) result(plan)

    sll_int32                                   :: nb_pts, ierr
    sll_real64, dimension(:), intent(in)        :: knots
    sll_real64, dimension(:)                    :: f
    type(quintic_splines_plan_non_uni), pointer :: plan

    ! Plan allocation
    SLL_ALLOCATE(plan, ierr)
    ! plan component allocation
    nb_pts = size(knots)
    SLL_ALLOCATE(plan%coeffs(nb_pts+5), ierr)
    plan%spline_obj=>new_arbitrary_degree_spline_1d(5, knots, nb_pts, 1)
    call compute_coeffs_non_uni(f, plan)

  end function new_quintic_splines_non_uni


  subroutine compute_coeffs_non_uni(f, plan_splines)

  ! f is the vector of the values of the function 
  !  in the nodes of the mesh*/

    sll_real64, dimension(:)                        :: f
    type(quintic_splines_plan_non_uni), pointer     :: plan_splines
    sll_real64                                      :: a, b, c
    sll_int32                                       :: nb_pts
    type(toep_penta_diagonal_plan), pointer         :: plan_pent
    sll_real64, dimension(6)                        :: b_at_x

    nb_pts = plan_splines%spline_obj%num_pts
    b_at_x = b_splines_at_x( plan_splines%spline_obj, nb_pts-1, &
                            plan_splines%spline_obj%xmax )
    a = b_at_x(3)
    b = b_at_x(2)
    c = b_at_x(1)

    plan_pent => new_toep_penta_diagonal(nb_pts+5)
    plan_splines%coeffs = solve_toep_penta_diagonal(a, b, c, f, plan_pent)
    call delete_toep_penta_diagonal(plan_pent)

  end subroutine compute_coeffs_non_uni

  function quintic_splines_non_uni(x, plan_splines) result(s)
  ! The interpolator spline function

    type(quintic_splines_plan_non_uni), pointer            :: plan_splines
    sll_int32                                              :: nb_pts, cell, left, j
    sll_real64                                             :: x, s
    sll_real64, dimension(6)                               :: b
    sll_real64, dimension(plan_splines%spline_obj%num_pts) :: knots

    knots = plan_splines%spline_obj%k
    nb_pts = plan_splines%spline_obj%num_pts

    cell = 1
    do while(  ( (x<knots(cell)) .or. (x>knots(cell+1)) ) .and. (cell<nb_pts-1)  )
         cell = cell + 1
    enddo

    left = cell - 1
    s = 0
    do j=left-5,left
      if( (j>=-5) .and. (j<=nb_pts-1) ) then
        s = s + plan_splines%coeffs(j+6) * b(j-left+6)
      endif
    enddo

  end function quintic_splines_non_uni


  subroutine delete_quintic_splines_non_uni(plan)

    type(quintic_splines_plan_non_uni), pointer :: plan
    sll_int32                                   :: ierr

    SLL_DEALLOCATE_ARRAY(plan%coeffs, ierr)
    call delete_arbitrary_order_spline_1d( plan%spline_obj )
    SLL_DEALLOCATE_ARRAY(plan, ierr)
 
  end subroutine delete_quintic_splines_non_uni

end module sll_quintic_splines
