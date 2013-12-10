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

!**************************************************************
!
! Selalib 2012     
! Module: sll_quintic_splines.F90
!
!> @brief 
!> Selalib quintic splines interpolator
!
!> Last modification: Nov. 29, 2012
!   
!> @authors                    
!> Aliou DIOUF (aliou.l.diouf@inria.fr)
!                                  
!**************************************************************

module sll_quintic_splines
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
use sll_toep_penta_diagonal
use sll_arbitrary_degree_splines
implicit none

  type quintic_splines_uniform_plan
    sll_int32                               :: num_pts
    sll_real64                              :: xmin
    sll_real64                              :: xmax
    sll_real64, dimension(6)                :: b_at_node
    sll_real64, dimension(:), pointer       :: coeffs
    type(toep_penta_diagonal_plan), pointer :: plan_pentadiagonal
  end type quintic_splines_uniform_plan

  type quintic_splines_nonuniform_plan
    sll_int32                                 :: num_pts
    sll_real64                                :: xmin
    sll_real64                                :: xmax
    sll_real64, dimension(:), pointer         :: coeffs
    sll_real64, dimension(:,:), pointer       :: matrix
    sll_real64, dimension(:), pointer         :: ipiv ! for matrix LU solving
    type(arbitrary_degree_spline_1d), pointer :: spline_obj
  end type quintic_splines_nonuniform_plan  

  interface compute_quintic_coeffs
     module procedure compute_quintic_coeffs_uniform, &
          compute_quintic_coeffs_nonuniform
  end interface compute_quintic_coeffs

  interface quintic_splines
     module procedure quintic_splines_interpolator_uniform_value, &
                         quintic_splines_interpolator_nonuniform_value
  end interface quintic_splines

  interface delete_quintic_splines
     module procedure delete_quintic_splines_uniform, &
          delete_quintic_splines_nonuniform
  end interface delete_quintic_splines

contains

  ! *************************************************************************
  !
  !                       UNIFORM QUINTIC SPLINES
  !
  ! *************************************************************************

  function new_quintic_splines_uniform(num_pts, xmin, xmax) result(plan)
    type(quintic_splines_uniform_plan), pointer :: plan
    sll_int32, intent(in)                       :: num_pts
    sll_real64                                  :: xmin
    sll_real64                                  :: xmax
    sll_int32                                   :: ierr
    !sll_int32,  intent(in), optional            :: bc_type

    if( num_pts < 6 ) then
       print *, 'ERROR, new_quintic_splines_uniform: Because of the algorithm used, ', &
                 'this function is meant to be used with arrays that are at ', &
                 'least of size = 6'
       STOP 'new_quintic_splines_uniform'
    endif

    SLL_ALLOCATE(plan, ierr)
    SLL_ALLOCATE(plan%coeffs(-5:num_pts-1), ierr)

    plan%num_pts = num_pts
    plan%xmin = xmin
    plan%xmax = xmax
    plan%b_at_node = uniform_b_splines_at_x( 5, 0.d0 )
    plan%plan_pentadiagonal => new_toep_penta_diagonal(num_pts+5)

  end function new_quintic_splines_uniform


  subroutine compute_quintic_coeffs_uniform(f, plan)

  ! f is the vector of the values of the function 
  !  in the nodes of the mesh*/

    sll_real64, dimension(:)                    :: f
    type(quintic_splines_uniform_plan), pointer :: plan
    sll_real64                                  :: a, b, c
    sll_int32                                   :: num_pts

    num_pts = plan%num_pts 
    ! a is the value to be duplicated in the principal diagonal of the matrix
    ! b for the diagonals -1 and 1
    ! c for the diagonals -2 and 2
    a = plan%b_at_node(3)
    b = plan%b_at_node(2)
    c = plan%b_at_node(1)

    ! To solve the linear system, we need to include the values of f outside 
    ! the domain which are 0 because f is compact. We storage all values of 
    ! f (inside the domain and outside) in the right side vector of the linear system.
    ! This vector is inout for the linear system solver
    plan%coeffs = 0.d0
    plan%coeffs (-3:num_pts-4) = f

    call solve_toep_penta_diagonal(a, b, c, plan%coeffs, plan%plan_pentadiagonal)
    plan%coeffs = plan%plan_pentadiagonal%solution

  end subroutine compute_quintic_coeffs_uniform

  function quintic_splines_interpolator_uniform_value(x, plan) result(s)

    type(quintic_splines_uniform_plan), pointer :: plan
    sll_int32                                   :: n, left
    sll_real64                                  :: x, xmin, xmax
    sll_real64                                  :: h, s, t0
    sll_real64, dimension(6)                    :: b

    xmin = plan%xmin
    xmax = plan%xmax
    n = plan%num_pts - 1
    h = (xmax-xmin)/n

    ! Run some checks on the arguments.
    SLL_ASSERT(associated(plan))
    SLL_ASSERT(x >= xmin)
    SLL_ASSERT(x <= xmax)

    t0 = (x-xmin)/h
    left = int(t0) ! Determine the leftmost support index 'i' of x
    t0 = t0 - left ! compute normalized_offset
    b = uniform_b_splines_at_x( 5, t0 )

    s = dot_product( plan%coeffs(left-5:left), b )

  end function quintic_splines_interpolator_uniform_value

  function quintic_splines_interpolator_uniform_array(array, num_pts, plan) result(res)
  
    sll_real64, dimension(:)                    :: array
    type(quintic_splines_uniform_plan), pointer :: plan
    sll_int32                                   :: i, num_pts
    sll_real64, dimension(num_pts)              :: res

    do i=1,num_pts
       res(i) = quintic_splines_interpolator_uniform_value(array(i), plan)
    enddo

  end function quintic_splines_interpolator_uniform_array

  function quintic_splines_interpolator_uniform_pointer(ptr, num_pts, plan) result(res)
  
    sll_real64, dimension(:), pointer           :: ptr
    type(quintic_splines_uniform_plan), pointer :: plan
    sll_int32                                   :: i, num_pts
    sll_real64, dimension(:), pointer           :: res

    res => ptr

    do i=1,num_pts
       res(i) = quintic_splines_interpolator_uniform_value(ptr(i), plan)
    enddo

  end function quintic_splines_interpolator_uniform_pointer

  subroutine delete_quintic_splines_uniform(plan)

    type(quintic_splines_uniform_plan), pointer :: plan
    sll_int32                                   :: ierr

    SLL_ASSERT(associated(plan))
    call delete_toep_penta_diagonal(plan%plan_pentadiagonal)
    SLL_DEALLOCATE_ARRAY(plan%coeffs, ierr)
    SLL_DEALLOCATE_ARRAY(plan, ierr)
 
  end subroutine delete_quintic_splines_uniform

  ! *************************************************************************
  !
  !                     NON UNIFORM QUINTIC SPLINES
  !
  ! *************************************************************************

  ! num_pts = nb_cells + 1
  function new_quintic_splines_nonuniform(knots) result(plan)

    sll_int32                                      :: num_pts, ierr
    sll_real64, dimension(:), intent(in)           :: knots
    sll_real64, dimension(:), allocatable          :: knots_fictive
    type(quintic_splines_nonuniform_plan), pointer :: plan
    sll_real64, dimension(:,:), allocatable        :: A
    sll_int32                                      :: i, j, m, cell
    sll_int32                                      :: KL, KU, LDAB
    sll_real64, dimension(6)                       :: b_at_x

    num_pts = size(knots)
    m = num_pts + 5

    if( num_pts < 6 ) then
       print *, 'ERROR, new_quintic_splines_nonuniform: Because of the algorithm used, ', &
                 'this function is meant to be used with arrays that are at ', &
                 'least of size = 6'
       STOP 'new_quintic_splines_nonuniform'
    endif

    ! Plan allocation
    SLL_ALLOCATE(plan, ierr)

    plan%num_pts = num_pts
    plan%xmin = knots(1)
    plan%xmax = knots(num_pts)

    ! plan component allocation
    SLL_ALLOCATE(plan%coeffs(-5:num_pts-1), ierr)
    SLL_ALLOCATE(plan%ipiv(m), ierr)
    SLL_ALLOCATE(knots_fictive(-5:num_pts+4), ierr)
    SLL_ALLOCATE(A(m,m), ierr)
    SLL_ALLOCATE(plan%matrix(m,m), ierr)

    do i=-5,-1
       knots_fictive(i) = 2*knots(1) - knots(-i+1)
    enddo
    knots_fictive(0:num_pts-1) = knots
    do i=num_pts,num_pts+4
       knots_fictive(i) = 2*knots(num_pts) - knots(2*num_pts-i-1) 
    enddo

    plan%spline_obj=>new_arbitrary_degree_spline_1d(5, knots_fictive, &
                                                size(knots_fictive), 1)

    KL = 2 ! for LAPACK use
    KU = 2 ! for LAPACK use
    A = 0.d0

    do i=-2,num_pts+2
       cell = i+6
       if ( i+6 == size (knots_fictive) ) then
          cell = cell - 1
       endif
       b_at_x = b_splines_at_x( plan%spline_obj, cell, knots_fictive(i) )
       do j= -KL, KU
          if ( (i+3+j>0) .and. (i+3+j<=m) ) then
             A(i+3,i+3+j) = b_at_x(j+KL+1) 
          endif
       enddo
    enddo    

    ! For linear system solving with LAPACK
    do j=1,m
       do i=max(1,j-ku), min(m,j+kl)
          plan%matrix(kl+ku+1+i-j,j) = A(i,j)
       enddo
    enddo
    LDAB = size(plan%matrix,1)
    ! LAPACK's LU factorization
    call DGBTRF( m, m, KL, KU, plan%matrix, LDAB, plan%IPIV, ierr )

    SLL_DEALLOCATE_ARRAY(A, ierr)
    SLL_DEALLOCATE_ARRAY(knots_fictive, ierr)

  end function new_quintic_splines_nonuniform


  subroutine compute_quintic_coeffs_nonuniform(f, plan)

  ! f is the vector of the values of the function 
  !  in the nodes of the mesh*/

    sll_real64, dimension(:)                       :: f
    type(quintic_splines_nonuniform_plan), pointer :: plan
    sll_int32                                      :: n, ierr
    sll_int32                                      :: KL, KU, LDAB, m

    n = plan%num_pts - 1
    m = plan%num_pts + 5
    KL = 2 ! for LAPACK use
    KU = 2 ! for LAPACK use
    LDAB = size(plan%matrix,1) ! for LAPACK use

    ! Solve the linear system with LAPACK's LU
    plan%coeffs = 0.d0
    plan%coeffs(-3:n-3) = f

    call DGBTRS( 'N', m, KL, KU, 1, plan%matrix, LDAB, plan%IPIV, plan%coeffs, m, ierr)

  end subroutine compute_quintic_coeffs_nonuniform

  function quintic_splines_interpolator_nonuniform_value(x, plan) &
    result(s)

    type(quintic_splines_nonuniform_plan), pointer :: plan
    sll_int32                                      :: cell, left
    sll_real64                                     :: x, s
    sll_real64, dimension(6)                       :: b

    ! Run some checks on the arguments.
    SLL_ASSERT(associated(plan))
    SLL_ASSERT(x >= plan%xmin)
    SLL_ASSERT(x <= plan%xmax)

    ! This is for finding the cell containing x. It has to be improved regarding the computing time:
    !--------------------------------------------------------------------------
    cell = 1
    do while( x < plan%spline_obj%k(cell) .or. x >= plan%spline_obj%k(cell+1) )
         cell = cell + 1
    enddo
    !--------------------------------------------------------------------------

    b = b_splines_at_x( plan%spline_obj, cell, x )
    left = cell - 6 ! left in the real mesh and not in the fictive one.
    ! In fact, cell is the cell number of x in the fictive mesh and left
    ! have to be for the left node in the real mesh.

    s = dot_product( plan%coeffs(left-5:left), b )

  end function quintic_splines_interpolator_nonuniform_value

  function quintic_splines_interpolator_nonuniform_array(array, &
                            num_pts, plan) result(res)
  
    sll_real64, dimension(:)                    :: array
    type(quintic_splines_nonuniform_plan), pointer :: plan
    sll_int32                                   :: i, num_pts
    sll_real64, dimension(num_pts)              :: res

    do i=1,num_pts
       res(i) = quintic_splines_interpolator_nonuniform_value( &
                                      array(i), plan)
    enddo

  end function quintic_splines_interpolator_nonuniform_array

  function quintic_splines_interpolator_nonuniform_pointer(ptr, &
                            num_pts, plan) result(res)
  
    sll_real64, dimension(:), pointer           :: ptr
    type(quintic_splines_nonuniform_plan), pointer :: plan
    sll_int32                                   :: i, num_pts
    sll_real64, dimension(:), pointer           :: res

    res => ptr

    do i=1,num_pts
       res(i) = quintic_splines_interpolator_nonuniform_value( &
                                        ptr(i), plan)
    enddo

  end function quintic_splines_interpolator_nonuniform_pointer

  subroutine delete_quintic_splines_nonuniform(plan)

    type(quintic_splines_nonuniform_plan), pointer :: plan
    sll_int32                                   :: ierr

    SLL_ASSERT(associated(plan))
    SLL_DEALLOCATE_ARRAY(plan%coeffs, ierr)
    call delete( plan%spline_obj )
    SLL_DEALLOCATE_ARRAY(plan, ierr)
 
  end subroutine delete_quintic_splines_nonuniform

end module sll_quintic_splines
