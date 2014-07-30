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

!***************************************************************************
!
! Selalib 2012     
! Module: sll_odd_degree_splines_uniform.F90
!
!> @brief 
!> Selalib odd degree splines interpolator
!
!> Start date: July 26, 2012
!> Last modification: Nov 30, 2012
!   
!> @authors                    
!> Aliou DIOUF (aliou.l.diouf@inria.fr)
!                                  
!***************************************************************************

module sll_odd_degree_splines

#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_assert.h"
use sll_arbitrary_degree_splines
  implicit none

  type odd_degree_splines_uniform_plan
    sll_int32                             :: num_pts
    sll_int32                             :: degree
    sll_real64                            :: xmin
    sll_real64                            :: xmax
    sll_real64, dimension(:), pointer     :: coeffs
    sll_real64, dimension(:,:), pointer   :: matrix
   ! matrix will be the result of Choleski factorization
   end type odd_degree_splines_uniform_plan

  type odd_degree_splines_nonuniform_plan
    sll_int32                             :: num_pts
    sll_int32                             :: degree
    sll_real64                            :: xmin
    sll_real64                            :: xmax
    sll_real64, dimension(:), pointer       :: coeffs
    sll_real64, dimension(:), pointer       :: b_spline
    sll_real64, dimension(:,:), pointer     :: matrix
    sll_real64, dimension(:), pointer       :: ipiv ! for matrix LU solving
    type(arbitrary_degree_spline_1d), pointer :: spline_obj
  end type odd_degree_splines_nonuniform_plan 

  interface compute_odd_degree_coeffs
     module procedure compute_odd_degree_coeffs_uniform, &
          compute_odd_degree_coeffs_nonuniform
  end interface compute_odd_degree_coeffs

  interface odd_degree_splines
     module procedure odd_degree_splines_interpolator_uniform_value, &
                         odd_degree_splines_interpolator_nonuniform_value
  end interface odd_degree_splines

  interface delete_odd_degree_splines
     module procedure delete_odd_degree_splines_uniform, &
          delete_odd_degree_splines_nonuniform
  end interface delete_odd_degree_splines

contains 

  function new_odd_degree_splines_uniform(num_pts,degree,xmin,xmax)result(plan)

    type(odd_degree_splines_uniform_plan), pointer :: plan
    sll_int32                                      :: ierr, num_pts, degree
    sll_real64                                     :: xmin, xmax
    sll_real64, dimension(:), allocatable          :: matrix_elements
!    sll_real64, dimension(:,:), allocatable        :: A
    sll_int32                                      :: KD, i, j, m, LDAB

    if (mod(degree,2) == 0) then
       print*, 'This needs to run with odd spline degree'
       print*, 'Exiting...'
       stop
    endif

    m = num_pts+degree
    KD = degree/2 ! called KD for lapack use

    SLL_ALLOCATE(plan, ierr)
    SLL_ALLOCATE(matrix_elements(degree+1), ierr)
!    SLL_ALLOCATE(A(m,m), ierr)
    !SLL_ALLOCATE(plan%matrix(m,m), ierr)
!  print *,'#before allocation',KD+1,m
    SLL_ALLOCATE(plan%matrix(KD+1,m), ierr) 
    !please allocate no full matrix, when not necessary...
!  print *,'#after allocation',KD+1,m
    SLL_ALLOCATE(plan%coeffs(-degree:num_pts-1), ierr)

    plan%num_pts = num_pts
    plan%degree = degree
    plan%xmin = xmin
    plan%xmax = xmax

    matrix_elements = uniform_b_splines_at_x( degree, 0.d0 ) 

!    A = 0.d0
!    do i=1,m
!       do j= -KD, KD
!          if ( (i+j>0) .and. (i+j<=m) ) then
!             A(i,i+j) = matrix_elements(j+KD+1) !i+j=k -> j=k-i if(k>0 and k<=m)
!          endif
!       enddo
!    enddo

    ! For Lapack use

    do j=1,m
       do i=j,min(m,j+KD)
          plan%matrix(1+i-j,j) = matrix_elements(j-i+KD+1)!A(i,j)
       enddo
    enddo

    ! Cholesky factorization with
    LDAB = size(plan%matrix,1)
    call DPBTRF( 'L', m, KD, plan%matrix, LDAB, ierr )

    SLL_DEALLOCATE_ARRAY(matrix_elements, ierr)
!    SLL_DEALLOCATE_ARRAY(A, ierr)

  end function new_odd_degree_splines_uniform


  subroutine compute_odd_degree_coeffs_uniform(f, plan)

    ! f is the vector of the values of the function 
    !  in the nodes of the mesh

    sll_real64, dimension(:)                        :: f
    type(odd_degree_splines_uniform_plan), pointer  :: plan
    sll_int32                                       :: deg, n, m
    sll_int32                                       :: KD, LDAB, ierr
    
    deg = plan%degree
    n = plan%num_pts - 1

    ! plan%coeffs is the solution vector of the linear system
    ! It is inout for lapack routine
    plan%coeffs = 0.d0
    plan%coeffs(deg/2-deg:deg/2-deg+n) = f

    ! Solve the linear system with LAPACK

    m = size(plan%coeffs)
    KD = deg/2
    LDAB = size(plan%matrix,1)

    ! Solve the linear system with Cholesky factorization
    call DPBTRS( 'L', m, KD, 1, plan%matrix, LDAB, plan%coeffs, m, ierr )

  end subroutine compute_odd_degree_coeffs_uniform



  !> \f[ 
  !> s(x_i) = f(x_i) 
  !       <=>
  !> c_{i-degree}*b_{i-degree} + ... + c_i*b_i = f(x_i), i=-degree..n
  !> c_j=0 if j<i-degree or j>n  
  !> \f]
  function odd_degree_splines_interpolator_uniform_value(x, plan) result(s)

    sll_real64                                     :: x, xmin, xmax
    sll_real64                                     :: h, s
    type(odd_degree_splines_uniform_plan), pointer :: plan
    sll_int32                                      :: n, left, deg
    sll_real64, dimension(plan%degree+1)           :: b
    sll_real64                                     :: t0

    xmin = plan%xmin
    xmax = plan%xmax
    n = plan%num_pts - 1
    deg = plan%degree
    h = (xmax-xmin)/n

    ! Run some checks on the arguments.
    SLL_ASSERT(associated(plan))
    SLL_ASSERT(x >= xmin)
    SLL_ASSERT(x <= xmax)

    t0 = (x-xmin)/h
    left = int(t0) ! Determine the leftmost support index 'i' of x
    t0 = t0 - left ! compute normalized_offset
    b = uniform_b_splines_at_x( deg, t0 )

    s = dot_product( plan%coeffs(left-deg:left), b )

  end function odd_degree_splines_interpolator_uniform_value

  function odd_degree_splines_interpolator_uniform_array(array, &
                               num_pts, plan_splines) result(res)
  
    sll_real64, dimension(:)                       :: array
    type(odd_degree_splines_uniform_plan), pointer :: plan_splines
    sll_int32                                      :: i, num_pts
    sll_real64, dimension(num_pts)                 :: res

    do i=1,num_pts
       res(i) = odd_degree_splines_interpolator_uniform_value( &
                                         array(i), plan_splines)
    enddo

  end function odd_degree_splines_interpolator_uniform_array

  function odd_degree_splines_interpolator_uniform_pointer(ptr, &
                               num_pts, plan_splines) result(res)
  
    sll_real64, dimension(:), pointer              :: ptr
    type(odd_degree_splines_uniform_plan), pointer :: plan_splines
    sll_int32                                      :: i, num_pts
    sll_real64, dimension(:), pointer              :: res

    res => ptr

    do i=1,num_pts
       res(i) = odd_degree_splines_interpolator_uniform_value( &
                                           ptr(i), plan_splines)
    enddo

  end function odd_degree_splines_interpolator_uniform_pointer


  subroutine delete_odd_degree_splines_uniform(plan)

    type(odd_degree_splines_uniform_plan), pointer :: plan
    sll_int32                                      :: ierr

    SLL_ASSERT(associated(plan))
    SLL_DEALLOCATE_ARRAY(plan%coeffs, ierr)
    SLL_DEALLOCATE_ARRAY(plan%matrix, ierr)
    SLL_DEALLOCATE_ARRAY(plan, ierr)

  end subroutine delete_odd_degree_splines_uniform

  ! *************************************************************************
  !
  !                  NON UNIFORM ODD DEGREE SPLINES STUFFS
  !
  ! *************************************************************************

  ! num_pts = nb_cells + 1
  function new_odd_degree_splines_nonuniform(degree, knots) result(plan)

    sll_int32                                         :: degree , num_pts, ierr
    sll_real64, dimension(:), intent(in)              :: knots
    sll_real64, dimension(:), allocatable             :: knots_fictive
    type(odd_degree_splines_nonuniform_plan), pointer :: plan
    sll_real64, dimension(:,:), allocatable           :: A
    sll_int32                                         :: i, ii, j, m, cell
    sll_int32                                         :: KL, KU, LDAB

    num_pts = size(knots)
    m = num_pts + degree

    if( num_pts < degree+1 ) then
       print *, 'ERROR, new_odd_degree_splines_nonuniform: Because of the algorithm used, ', &
                 'this function is meant to be used with arrays that are at ', &
                 'least of size = ', degree+1
       STOP 'new_odd_degree_splines_nonuniform'
    endif

    ! Plan allocation
    SLL_ALLOCATE(plan, ierr)

    plan%num_pts = num_pts
    plan%xmin = knots(1)
    plan%xmax = knots(num_pts)
    plan%degree = degree

    ! plan component allocation
    SLL_ALLOCATE(plan%coeffs(-degree:num_pts-1), ierr)
    SLL_ALLOCATE(plan%ipiv(m), ierr)
    SLL_ALLOCATE(knots_fictive(-degree:num_pts-1+degree+1/degree), ierr) ! When degree = 1,
    ! adding 1 point beyond xmax would be not enough
    ! So instead of 1 point beyond xmax, we add 2 points, for degree 1
    SLL_ALLOCATE(A(m,m), ierr)
    SLL_ALLOCATE(plan%matrix(m,m), ierr) !A and plan%matrix(m,m) 
    ! should not be full matrix...
    SLL_ALLOCATE(plan%b_spline(degree+1), ierr)

    do i=-degree,-1
       knots_fictive(i) = 2*knots(1) - knots(-i+1)
    enddo
    knots_fictive(0:num_pts-1) = knots
    do i=num_pts,num_pts-1+degree+1/degree
       knots_fictive(i) = 2*knots(num_pts) - knots(2*num_pts-i-1+1/degree) + 1/degree
    enddo

    plan%spline_obj=>new_arbitrary_degree_spline_1d(degree, knots_fictive, &
                                                     size(knots_fictive), 1)

    KL = degree/2 ! for LAPACK use
    KU = degree/2 ! for LAPACK use
    A = 0.d0

    do i=degree/2-degree+1,num_pts+degree/2
       cell = i+degree+1
       if ( i+degree+1 == size (knots_fictive) ) then
          cell = cell - 1
       endif
       plan%b_spline = b_splines_at_x( plan%spline_obj, cell, knots_fictive(i) )
       ii = i + degree - degree/2
       do j= -KL, KU
          if ( (ii+j>0) .and. (ii+j<=m) ) then
             A(ii,ii+j) = plan%b_spline(j+KL+1) !
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

  end function new_odd_degree_splines_nonuniform


  subroutine compute_odd_degree_coeffs_nonuniform(f, plan)

  ! f is the vector of the values of the function 
  !  in the nodes of the mesh*/

    sll_real64, dimension(:)                          :: f
    type(odd_degree_splines_nonuniform_plan), pointer :: plan
    sll_int32                                         :: n, ierr, deg
    sll_int32                                         :: KL, KU, LDAB, m

    deg = plan%degree
    n = plan%num_pts - 1
    m = size(plan%coeffs)
    KL = deg/2 ! for LAPACK use
    KU = deg/2 ! for LAPACK use
    LDAB = size(plan%matrix,1) ! for LAPACK use

    ! Solve the linear system with LAPACK's LU
    plan%coeffs = 0.d0
    plan%coeffs(deg/2-deg:deg/2-deg+n) = f

    call DGBTRS( 'N', m, KL, KU, 1, plan%matrix, LDAB, plan%IPIV, plan%coeffs, m, ierr)

  end subroutine compute_odd_degree_coeffs_nonuniform

  function odd_degree_splines_interpolator_nonuniform_value(x, plan) result(s)

    type(odd_degree_splines_nonuniform_plan), pointer :: plan
    sll_int32                                         :: n, cell, left
    sll_real64                                        :: x, s
    sll_int32                                         :: deg

    ! Run some checks on the arguments.
    SLL_ASSERT(associated(plan))
    SLL_ASSERT(x >= plan%spline_obj%xmin)
    SLL_ASSERT(x <= plan%spline_obj%xmax)

    n = plan%spline_obj%num_pts - 1
    deg = plan%degree

    ! This is for finding the cell containing x. It has to be improved regarding the computing time:
    !--------------------------------------------------------------------------
    cell = 1
    do while( x < plan%spline_obj%k(cell) .or. x >= plan%spline_obj%k(cell+1) )
         cell = cell + 1
    enddo
    !--------------------------------------------------------------------------

    plan%b_spline = b_splines_at_x( plan%spline_obj, cell, x )
    left = cell - (deg+1) ! left in the real mesh and not in the fictive one.
    ! In fact, cell is the cell number of x in the fictive mesh and left
    ! have to be for the left node in the real mesh.

    s = dot_product( plan%coeffs(left-deg:left), plan%b_spline )

  end function odd_degree_splines_interpolator_nonuniform_value

  function odd_degree_splines_interpolator_nonuniform_array(array, &
                            num_pts, plan_splines) result(res)
  
    sll_real64, dimension(:)                    :: array
    type(odd_degree_splines_nonuniform_plan), pointer :: plan_splines
    sll_int32                                   :: i, num_pts
    sll_real64, dimension(num_pts)              :: res

    do i=1,num_pts
       res(i) = odd_degree_splines_interpolator_nonuniform_value( &
                                         array(i), plan_splines)
    enddo

  end function odd_degree_splines_interpolator_nonuniform_array

  function odd_degree_splines_interpolator_nonuniform_pointer(ptr, &
                            num_pts, plan_splines) result(res)
  
    sll_real64, dimension(:), pointer           :: ptr
    type(odd_degree_splines_nonuniform_plan), pointer :: plan_splines
    sll_int32                                   :: i, num_pts
    sll_real64, dimension(:), pointer           :: res

    res => ptr

    do i=1,num_pts
       res(i) = odd_degree_splines_interpolator_nonuniform_value( &
                                           ptr(i), plan_splines)
    enddo

  end function odd_degree_splines_interpolator_nonuniform_pointer

  subroutine delete_odd_degree_splines_nonuniform(plan)

    type(odd_degree_splines_nonuniform_plan), pointer :: plan
    sll_int32                                   :: ierr

    SLL_ASSERT(associated(plan))
    SLL_DEALLOCATE_ARRAY(plan%coeffs, ierr)
    call delete( plan%spline_obj )
    SLL_DEALLOCATE_ARRAY(plan, ierr)
 
  end subroutine delete_odd_degree_splines_nonuniform

end module sll_odd_degree_splines
