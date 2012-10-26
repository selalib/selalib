!***************************************************************************
!
! Selalib 2012     
! Module: sll_odd_degree_splines_uniform.F90
!
!> @brief 
!> Selalib odd degree splines interpolator
!
!> Start date: July 26, 2012
!> Last modification: October 25, 2012
!   
!> @authors                    
!> Aliou DIOUF (aliou.l.diouf@inria.fr)
!                                  
!***************************************************************************

module sll_odd_degree_splines

#include "sll_memory.h"
#include "sll_working_precision.h"
use arbitrary_degree_splines
  implicit none

  type odd_degree_splines_uniform_plan
    sll_int32                             :: num_pts
    sll_int32                             :: degree
    sll_real64                            :: xmin
    sll_real64                            :: xmax
<<<<<<< HEAD
#ifdef STDF95
    sll_real64, dimension(:), pointer :: b_at_node
    sll_real64, dimension(:), pointer :: coeffs
#else
    sll_real64, dimension(:), allocatable :: b_at_node
    sll_real64, dimension(:), allocatable :: coeffs
#endif
  end type odd_degree_splines_plan
=======
    sll_real64, dimension(:), pointer     :: b_at_node
    sll_real64, dimension(:), allocatable :: coeffs
  end type odd_degree_splines_uniform_plan
>>>>>>> origin/prototype-devel

contains 

  function new_odd_degree_splines_uniform(num_pts,degree,xmin,xmax)result(plan)

    type(odd_degree_splines_uniform_plan), pointer :: plan
    sll_int32                                      :: ierr, num_pts, degree
    sll_real64                                     :: xmin, xmax

    if (mod(degree,2) == 0) then
       print*, 'This needs to run with odd spline degree'
       print*, 'Exiting...'
       stop
    endif

    SLL_ALLOCATE(plan, ierr)
    SLL_ALLOCATE(plan%b_at_node(degree+1), ierr)
    SLL_ALLOCATE(plan%coeffs(num_pts+degree), ierr)

    plan%num_pts = num_pts
    plan%degree = degree
    plan%xmin = xmin
    plan%xmax = xmax
    plan%b_at_node = uniform_b_splines_at_x( degree, 0.d0 )

  end function new_odd_degree_splines_uniform


  subroutine compute_coeffs_uniform(f, plan)

    ! f is the vector of the values of the function 
    !  in the nodes of the mesh

    sll_real64, dimension(:)                        :: f
    type(odd_degree_splines_uniform_plan), pointer  :: plan
    sll_real64, dimension(plan%num_pts+plan%degree) :: g
    sll_real64                                      :: xmin, xmax, h
    sll_int32                                       :: degree, n, m, i, j
    sll_real64, dimension(size(g),size(g))          :: A, AB
    sll_int32                                       :: KD, LDAB, ierr
    
    degree = plan%degree
    xmin = plan%xmin
    xmax = plan%xmax
    n = plan%num_pts - 1
    h = (xmax-xmin)/n

    g = 0.d0
    g(degree/2+1:n+degree/2+1) = f

    ! Solve the linear system with LAPACK

    m = size(g)
    KD = degree/2 ! called KD for lapack use

    A = 0.d0
    do i=1,m
       do j= -KD, KD
          if ( (i+j>0) .and. (i+j<=m) ) then
             A(i,i+j) = plan%b_at_node(j+KD+1) 
          endif
       enddo
    enddo

    do j=1,m
       do i=j,min(m,j+KD)
          AB(1+i-j,j) = A(i,j)
       enddo
    enddo

    LDAB = size(AB,1)
    ! Cholesky factorization
    call DPBTRF( 'L', m, KD, AB, LDAB, ierr )
    ! Solve the linear system with Cholesky factorization
    plan%coeffs = g
    call DPBTRS( 'L', m, KD, 1, AB, LDAB, plan%coeffs, m, ierr )
!print*, sum(f-matmul(A,plan%coeffs))
  end subroutine compute_coeffs_uniform



  !> s(x_i) = f(x_i) 
  !       <=>
  !> c_{i-degree}*b_{i-degree} + ... + c_i*b_i = f(x_i), i=-degree..n
  !> c_j=0 if j<i-degree or j>n  
  function odd_degree_splines_interpolator_uniform_value(x, plan) result(s)

    sll_real64                                     :: x, xmin, xmax
    sll_real64                                     :: h, s
    type(odd_degree_splines_uniform_plan), pointer :: plan
    sll_int32                                      :: n, j, left, degree
    sll_real64, dimension(plan%degree+1)           :: b
    sll_real64                                     :: t0

    xmin = plan%xmin
    xmax = plan%xmax
    n = plan%num_pts - 1
    degree = plan%degree
    h = (xmax-xmin)/n

    t0 = (x-xmin)/h
    left = int(t0) ! Determine the leftmost support index 'i' of x
    t0 = t0 - left ! compute normalized_offset

    b = uniform_b_splines_at_x( degree, t0 )
    s = 0.d0

    do j=left-degree,left
       if ( (j>=-degree) .and. (j<=n) ) then
          s = s + plan%coeffs(j+degree+1) * b(j-left+degree+1)
       endif
    enddo

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

    SLL_DEALLOCATE_ARRAY(plan%coeffs, ierr)
    SLL_DEALLOCATE_ARRAY(plan%b_at_node, ierr)
    SLL_DEALLOCATE_ARRAY(plan, ierr)

  end subroutine delete_odd_degree_splines_uniform

end module sll_odd_degree_splines
