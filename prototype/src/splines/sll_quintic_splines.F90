!**************************************************************
!
! Selalib 2012     
! Module: sll_quintic_splines.F90
!
!> @brief 
!> Selalib quintic splines interpolator
!
!> Last modification: Nov. 20, 2012
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
use arbitrary_degree_splines
implicit none

  type quintic_splines_uniform_plan
    sll_int32                         :: num_pts
    sll_real64                        :: xmin
    sll_real64                        :: xmax
    sll_real64, dimension(6)          :: b_at_node
    sll_real64, dimension(:), pointer :: coeffs
    type(toep_penta_diagonal_plan), pointer :: plan_pentadiagonal
  end type quintic_splines_uniform_plan

  type quintic_splines_nonuniform_plan
    sll_real64, dimension(:), pointer :: coeffs
    type(arbitrary_degree_spline_1d), pointer :: spline_obj
  end type quintic_splines_nonuniform_plan  

  interface compute_quintic_coeffs
     module procedure compute_quintic_coeffs_uniform, &
          compute_quintic_coeffs_nonuniform
  end interface compute_quintic_coeffs

  ! What is this?
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

    SLL_ALLOCATE(plan, ierr)
    SLL_ALLOCATE(plan%coeffs(num_pts+5), ierr)

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
    plan%coeffs(3:num_pts+2) = f

    plan%coeffs = solve_toep_penta_diagonal(a, b, c, plan%coeffs, plan%plan_pentadiagonal)

  end subroutine compute_quintic_coeffs_uniform

  function quintic_splines_interpolator_uniform_value(x, plan) result(s)

    type(quintic_splines_uniform_plan), pointer :: plan
    sll_int32                                   :: n, left, j
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
    s = 0.d0

    do j=left-5,left
      if( (j>=-5) .and. (j<=n) ) then
        s = s + plan%coeffs(j+6) * b(j-left+6)
      endif
    enddo

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

    sll_int32                                   :: num_pts, ierr
    sll_real64, dimension(:), intent(in)        :: knots
    type(quintic_splines_nonuniform_plan), pointer :: plan

    ! Plan allocation
    SLL_ALLOCATE(plan, ierr)
    ! plan component allocation
    num_pts = size(knots)
    SLL_ALLOCATE(plan%coeffs(num_pts+5), ierr)

    plan%spline_obj=>new_arbitrary_degree_spline_1d(5, knots, num_pts, 1)

  end function new_quintic_splines_nonuniform


  subroutine compute_quintic_coeffs_nonuniform(f, plan)

  ! f is the vector of the values of the function 
  !  in the nodes of the mesh*/

    sll_real64, dimension(:)                       :: f
    type(quintic_splines_nonuniform_plan), pointer :: plan
    sll_real64, dimension(size(f)+3)               :: g, ipiv
    sll_real64, dimension(size(f)+3,size(f)+3)     :: A, AB 
    sll_int32                                      :: num_pts, i, ierr
    sll_real64, dimension(6)                       :: b_at_x
    sll_int32                                      :: KL, KU, m, j, LDAB

    num_pts = plan%spline_obj%num_pts
    m = size(g)

    ! To solve the linear system, we need to include the values of f outside 
    ! the domain which are 0 because f is compact. We storage all values of 
    ! f (inside the domain and outside) in g.
    g = 0.d0
    g(2:num_pts+1) = f
    m = size(g)

    ! called for lapack use
    KL = 2
    KU = 2
    A = 0.d0

    do i=-1, num_pts+3

       if ( i <= 1 ) then
          b_at_x = b_splines_at_x( plan%spline_obj, 1, plan%spline_obj%k(1) )
       elseif ( i >= num_pts ) then
          b_at_x = b_splines_at_x( plan%spline_obj, num_pts-1, plan%spline_obj%k(num_pts) )
print*, b_at_x; stop
       else
          b_at_x = b_splines_at_x( plan%spline_obj, i, plan%spline_obj%k(i) )

       endif

       do j= -KL, KL
          if ( (i+KL+j>=1) .and. (i+KL+j<=m) ) then
             A(i+KL,i+KL+j) = b_at_x(j+KL+1) 
          endif
       enddo
    enddo

    ! Solve the linear system with LAPACK
    do j=1,m
       do i=max(1,j-ku), min(m,j+kl)
          AB(kl+ku+1+i-j,j) = A(i,j)
       enddo
    enddo

    LDAB = size(AB,1)
    ! LU factorization
    call DGBTRF( m, m, KL, KU, AB, LDAB, IPIV, ierr )
if (ierr/=0) print*, a, ierr; stop
    ! Solve the linear system with Cholesky factorization
    plan%coeffs = g
    call DGBTRS( 'N', m, KL, KU, 1, AB, LDAB, IPIV, plan%coeffs, m, ierr)

  end subroutine compute_quintic_coeffs_nonuniform

  function quintic_splines_interpolator_nonuniform_value(x, plan_splines) &
    result(s)

    type(quintic_splines_nonuniform_plan), pointer            :: plan_splines
    sll_int32                                              :: n, cell, left, j
    sll_real64                                             :: x, s
    sll_real64, dimension(6)                               :: b
    sll_real64, dimension(plan_splines%spline_obj%num_pts) :: knots
    sll_int32                                              :: ierr

    ! Run some checks on the arguments.
    SLL_ASSERT(associated(plan_splines))
    SLL_ASSERT(x >= plan_splines%spline_obj%xmin)
    SLL_ASSERT(x <= plan_splines%spline_obj%xmax)

    n = plan_splines%spline_obj%num_pts - 1
    knots = plan_splines%spline_obj%k(1:n+1)

    !cell = 1
    !do while( ( (x<knots(cell)) .or. (x>knots(cell+1)) ) .and. (cell<n) )
    !     cell = cell + 1
    !enddo
    call find_cell(x, knots, (/(j, j=0,n)/), cell, ierr)

    left = cell - 1
    s = 0
    do j=left-5,left
      if( (j>=-5) .and. (j<=n) ) then
        s = s + plan_splines%coeffs(j+6) * b(j-left+6)
      endif
    enddo

  end function quintic_splines_interpolator_nonuniform_value

  !> indices is array containing the indices of the mesh: 0, 1,..., num_pts-1
  recursive subroutine find_cell(x, knots, indices, cell, ierr)

    double precision               :: x
    double precision, dimension(:) :: knots
    integer, dimension(:)          :: indices
    integer                        :: num_pts, n, cell, ierr

    num_pts = size(knots)
    n = num_pts / 2
    ierr = 0

    if ( num_pts > 2 ) then
       call find_cell( x, knots(1:n), indices(1:n), cell, ierr )
       if (ierr==0) then
          call find_cell( x, knots(n+1:num_pts), &
                indices(n+1:num_pts), cell, ierr )
       endif
    elseif (num_pts==2) then
       if ( (knots(1)<=x) .and. (x<=knots(2)) ) then
          cell = indices(2)
          ierr = 1
       endif
    endif

  end subroutine find_cell

  function quintic_splines_interpolator_nonuniform_array(array, &
                            num_pts, plan_splines) result(res)
  
    sll_real64, dimension(:)                    :: array
    type(quintic_splines_nonuniform_plan), pointer :: plan_splines
    sll_int32                                   :: i, num_pts
    sll_real64, dimension(num_pts)              :: res

    do i=1,num_pts
       res(i) = quintic_splines_interpolator_nonuniform_value( &
                                      array(i), plan_splines)
    enddo

  end function quintic_splines_interpolator_nonuniform_array

  function quintic_splines_interpolator_nonuniform_pointer(ptr, &
                            num_pts, plan_splines) result(res)
  
    sll_real64, dimension(:), pointer           :: ptr
    type(quintic_splines_nonuniform_plan), pointer :: plan_splines
    sll_int32                                   :: i, num_pts
    sll_real64, dimension(:), pointer           :: res

    res => ptr

    do i=1,num_pts
       res(i) = quintic_splines_interpolator_nonuniform_value( &
                                        ptr(i), plan_splines)
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
