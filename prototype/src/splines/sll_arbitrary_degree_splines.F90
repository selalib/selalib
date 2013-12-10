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

! This module attempts to replace the functionalities of the arbitrary
! degree spline functions by deBoor. 

module sll_arbitrary_degree_splines
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
  implicit none

  ! With these read-only parameters, we mimic the behavior of an enumerator
  ! or object macros. These are the alternatives with a different 
  ! implementation.
  sll_int32, parameter :: PERIODIC_ARBITRARY_DEG_SPLINE = 0
  sll_int32, parameter :: OPEN_ARBITRARY_DEG_SPLINE     = 1 

  type :: arbitrary_degree_spline_1d
     sll_int32  :: num_pts
     sll_int32  :: bc_type
     sll_int32  :: degree
     sll_real64 :: xmin
     sll_real64 :: xmax
     sll_real64, dimension(:), pointer :: k   ! local copy of knots array
  end type arbitrary_degree_spline_1d

  interface delete
     module procedure delete_arbitrary_order_spline_1d
  end interface delete

contains

  function new_arbitrary_degree_spline_1d( degree, knots, num_pts, bc_type )
    type(arbitrary_degree_spline_1d), pointer :: new_arbitrary_degree_spline_1d
    sll_int32, intent(in)                    :: degree
    sll_real64, dimension(:), intent(in)     :: knots
    sll_int32,  intent(in)                   :: num_pts
    sll_int32,  intent(in)                   :: bc_type
    sll_int32                                :: i
    sll_int32                                :: ierr
    sll_real64                               :: period

    if( size(knots) < num_pts ) then
       print *, 'ERROR. new_arbitrary_degree_spline_1d(): ', &
            'size of given knots array is smaller than the stated ', &
            'number of points.'
       print *, 'size(knots) = ', size(knots), 'num_pts = ', num_pts
       STOP
    end if
    if( degree < 0 ) then
       print *, 'ERROR. new_arbitrary_degree_spline_1d(): ', &
            'only positive integer values for degree are allowed, ', &
            'given: ', degree
    end if
    select case( bc_type )
       case(:-1)
          print *, 'ERROR. new_arbitrary_degree_spline_1d(): ', &
               'invalid boundary condition type supplied.'
          STOP
       case(2:)
          print *, 'ERROR. new_arbitrary_degree_spline_1d(): ', &
               'invalid boundary condition type supplied.'
          STOP
    end select

    SLL_ALLOCATE( new_arbitrary_degree_spline_1d, ierr )
    new_arbitrary_degree_spline_1d%num_pts  = num_pts
    new_arbitrary_degree_spline_1d%bc_type  = bc_type
    new_arbitrary_degree_spline_1d%degree   = degree
    new_arbitrary_degree_spline_1d%xmin     = knots(1)
    new_arbitrary_degree_spline_1d%xmax     = knots(num_pts)

    ! 'period' is useful in the case of periodic boundary conditions.
    period = knots(num_pts) - knots(1)

    ! Build the local copy of the knots. Here we simply pad the local array
    ! with an amount of nodes that depends on the degree of the requested 
    ! spline. We aim at setting up the indexing in such a way that the 
    ! original indexing of 'knots' is preserved, i.e.: knots(i) = k(i), at
    ! least whenever the scope of the indices defined here is active.
    SLL_ALLOCATE(new_arbitrary_degree_spline_1d%k(1-degree:num_pts+degree),ierr)
    do i=1,num_pts
       new_arbitrary_degree_spline_1d%k(i) = knots(i)
    end do


    ! Fill out the extra points at both ends of the local knots array with
    ! values proper to the boundary condition requested.
    select case( bc_type )
       case(PERIODIC_ARBITRARY_DEG_SPLINE)
          ! The logic behind the periodic boundary condition is the following.
          ! The given knots array has minimum (knots(1)) and maximum (knots(n))
          ! values at either end. This defines a length 'L'. If interpreted 
          ! as a periodic space, this is also the period. Thus, as we extend
          ! the number of knots at both ends of the given array, we use the
          ! periodicity condition to fill out the new values:
          !
          !                    .
          !                    .
          !                    .
          !           knots(-1) = knots(n-2) - L
          !           knots( 0) = knots(n-1) - L
          !                    .
          !                    .
          !                    .
          !           knots(n+1) = knots(1) + L
          !           knots(n+2) = knots(2) + L
          !                    .
          !                    .
          !                    .
          !

          ! Note that the following loop is not supposed to get executed
          ! if degree == 0.
          do i=1,degree
             ! Fill out the extra nodes on the left
             new_arbitrary_degree_spline_1d%k(1-i) = knots(num_pts-i) - period
             ! Fill out the extra nodes on the right
             new_arbitrary_degree_spline_1d%k(num_pts+i) = knots(i+1) + period
          end do

       case(OPEN_ARBITRARY_DEG_SPLINE)
          ! The 'open' boundary condition simply extends the new values
          ! of the local array at both ends with repeated endpoint values.
          ! That is
          !
          !     ... = knots(-2) = knots(-1) = knots(0) = knots(1)
          !
          ! and
          !
          !    knots(n+1) = knots(n+2) = knots(n+3) = ... =  knots(n)
          do i=1-degree,0
             new_arbitrary_degree_spline_1d%k(i) = knots(1)
          end do
          do i=num_pts+1,num_pts+degree
             new_arbitrary_degree_spline_1d%k(i) = knots(num_pts)
          end do
    end select
  end function new_arbitrary_degree_spline_1d

  ! b_splines_at_x() returns the values of all the B-splines of a given 
  ! degree that have support in cell 'cell' and evaluated at the point 'x'. 
  ! In other words, if B[j,i](x) is the spline of degree 'j' whose leftmost 
  ! support is at cell 'i' and evaluated at 'x', then b_splines_at_x returns 
  ! the sequence (in the form of an array):
  ! 
  ! B[j,i-degree](x), B[j,i-degree+1](x), B[j,i-degree+2](x), ..., B[j,i](x)
  !
  ! Implementation notes: 
  !
  ! It would have been very simple and convenient to implement this with
  ! a recursion, since:
  !
  !               x-t(i)                       t(i+j+1)-x
  ! B[j,i](x) = ----------- * B[j-1,i](x) + ----------------* B[j-1,i+1](x)
  !             t(i+j)-t(i)                  t(i+j+1)-t(i+1)
  !
  ! and
  !
  ! B[0,i] = 1 if t(i) <= x < t(i+1), and 0 otherwise.
  !
  ! More generally:
  !
  ! if t(i) <= x < t(i+j+1) then the formula above applies but B[j,i](x) = 0
  ! otherwise.
  !
  ! The problem with the above recursion is that it will end up computing the
  ! splines of lower orders very redundantly, much like the problem of 
  ! computing a Fibonacci sequence with a recursion. For such a critical
  ! function (this will be present in inner loops), this is not acceptable.
  !
  ! Here we try a different approach but still use the idea of the 
  ! recursion formula above. We use the fact that for the desired sequence
  ! of splines of degree 'J', we need information within 2*J+1 cells:
  ! (i-J):(i+J). We populate these with the values of the B[0,i](x) splines
  ! and iteratively build the higher order splines as needed.

  !> b_splines_at_x( spline_obj, cell, x ) computes the values of all the
  !> splines which have support in 'cell' and evaluates them at point 'x',
  !> which is supposed to be in cell. The spline object should have already
  !> been initialized and will contain information on the spline degree
  !> to use and the type of boundary condition desired.
  function b_splines_at_x( spline_obj, cell, x )
    type(arbitrary_degree_spline_1d), pointer      :: spline_obj
    sll_int32, intent(in)                          :: cell
    sll_real64, intent(in)                         :: x
    sll_int32                                      :: deg
    sll_real64, dimension(1:spline_obj%degree+1)   :: b_splines_at_x
    sll_real64, dimension(1:2*spline_obj%degree+1) :: splines
    sll_int32                                      :: i
    sll_int32                                      :: j
    sll_int32                                      :: last
    sll_real64                                     :: ti     ! t(i)
    sll_real64                                     :: tip1   ! t(i+1)
    sll_real64                                     :: tipj   ! t(i+j)
    sll_real64                                     :: tipjp1 ! t(i+j+1)
    sll_real64                                     :: fac1
    sll_real64                                     :: fac2
    sll_real64                                     :: term1
    sll_real64                                     :: term2
    sll_int32                                      :: current

    ! Run some checks on the arguments.
    SLL_ASSERT(associated(spline_obj))
    SLL_ASSERT(x >= spline_obj%xmin)
    SLL_ASSERT(x <= spline_obj%xmax)
    SLL_ASSERT(cell >= 1)
    SLL_ASSERT(cell <= spline_obj%num_pts - 1)
    ! This is checked always.
    if( .not. ((x >= spline_obj%k(cell)).and.(x <= spline_obj%k(cell+1)))) then
       print *, 'ERROR. b_splines_at_x(): the given value of x is not ', &
            'inside the specified cell.'
       STOP
    end if
    deg = spline_obj%degree

    ! Build the zeroth-order splines. The middle cell of the splines array
    ! corresponds to the 'cell' given as argument. So, 
    ! splines(spline_degree+1) = 1.0.
    splines(:)     = 0.0
    splines(deg+1) = 1.0

    ! Build the higher order splines. 
    last = 2*deg
    do j=1,deg
       do i=1,last
          current    = cell - deg + i - 1
          ti         = spline_obj%k(current)
          tip1       = spline_obj%k(current+1)
          tipj       = spline_obj%k(current+j)
          tipjp1     = spline_obj%k(current+j+1)
          ! This is a dangerous situation for which we need some sort of
          ! protection: What guarantees are there that these denominators
          ! will not be zero?? This should probably be error-checked, else
          ! one can just end up with an array of NaN's.
          ! 10-29-12: well, apparently this can be protected against if we
          ! force the terms to zero whenever this happens because the 
          ! corresponding spline equals zero, hence the indetermination is 
          ! always forced to zero.
#if 0
          if( tipj == ti ) then
             print *, 'calculation is shot... NaNs will be found'
             print *, 'x = ', x
             print *, 'current = ', current
             print *, 'deg = ', deg
             print *, 'cell = ', cell
             print *, 'knot(current) = ', spline_obj%k(current)
             print *, 'knot(current+1) = ', spline_obj%k(current+1)
             print *, 'knot(current+j) = ', spline_obj%k(current+j),j
          end if

          if( tipjp1 == tip1 ) then
             print *, 'calculation is shot... NaNs will be found'
             print *, 'x = ', x
             print *, 'current = ', current
             print *, 'deg = ', deg
             print *, 'cell = ', cell
             print *, 'knot(current) = ', spline_obj%k(current)
             print *, 'knot(current+1) = ', spline_obj%k(current+1)
             print *, 'knot(current+j) = ', spline_obj%k(current+j),j
             print *, 'knot(current+j+1) = ', spline_obj%k(current+j+1)
          end if



#endif
  !
  !               x-t(i)                       t(i+j+1)-x
  ! B[j,i](x) = ----------- * B[j-1,i](x) + ----------------* B[j-1,i+1](x)
  !             t(i+j)-t(i)                  t(i+j+1)-t(i+1)

          if(tipj==ti)then
            fac1=0._f64
          elseif (abs(tipj-ti)<1.e-15) then
            print *,tipj-ti
          else     
            fac1       = (x - ti)/(tipj - ti)            
          endif  
          if(tipjp1==tip1)then
            fac2=0._f64
          elseif (abs(tipjp1-tip1)<1.e-15) then
            print *,tipjp1-tip1
          else
            fac2       = (tipjp1 - x)/(tipjp1 - tip1)
          endif
          
          
          ! Super-ugly step to eliminate those cases where the denominator is
          ! zero but the spline value is also zero, thus avoiding the 
          ! indeterminate result ( NaN's ). 
          if( splines(i) .ne. 0.0 ) then
             term1 = fac1*splines(i)
          else
             term1 = 0.0_f64
          end if
          if( splines(i+1) .ne. 0.0 ) then
             term2 = fac2*splines(i+1)
          else
             term2 = 0.0_f64
          end if
          splines(i) = term1 + term2
!          splines(i) = fac1*splines(i) + fac2*splines(i+1)
       end do
       last = last - 1
    end do
    b_splines_at_x(1:deg+1) = splines(1:deg+1)
  end function b_splines_at_x

  ! b_spline_derivatives_at_x() returns an array with the derivative values of 
  ! the B-splines of a requested order that are supported in 'cell' and 
  ! evaluated at 'x'. The return value has the format:
  !
  ! B'[j,i-degree](x), B'[j,i-degree+1](x), B'[j,i-degree+2](x),..., B'[j,i](x)
  !
  ! where 'j' is the requested degree of the spline.
  function b_spline_derivatives_at_x( spline_obj, cell, x )
    type(arbitrary_degree_spline_1d), pointer      :: spline_obj
    sll_int32, intent(in)                          :: cell
    sll_real64, intent(in)                         :: x
    sll_int32                                      :: deg
    sll_real64, dimension(1:spline_obj%degree+1)   :: b_spline_derivatives_at_x
    sll_real64, dimension(1:2*spline_obj%degree+1) :: derivs
    sll_real64, dimension(1:2*spline_obj%degree+1) :: splines 
    sll_int32                                      :: i
    sll_int32                                      :: j
    sll_int32                                      :: last
    sll_real64                                     :: ti     ! t(i)
    sll_real64                                     :: tip1   ! t(i+1)
    sll_real64                                     :: tipj   ! t(i+j)
    sll_real64                                     :: tipjp1 ! t(i+j+1)
    sll_real64                                     :: fac1
    sll_real64                                     :: fac2
    sll_real64                                     :: term1
    sll_real64                                     :: term2
    sll_real64                                     :: delta_left
    sll_real64                                     :: delta1, delta2
    sll_real64                                     :: delta_right
    sll_int32                                      :: current
    sll_int32                                      :: num_pts

    ! Run some checks on the arguments.
    SLL_ASSERT(associated(spline_obj))
    SLL_ASSERT(x >= spline_obj%xmin)
    SLL_ASSERT(x <= spline_obj%xmax)
    SLL_ASSERT(cell >= 1)
    SLL_ASSERT(cell <= spline_obj%num_pts - 1)
    ! This is checked always.
    if( .not. (x >= spline_obj%k(cell)) .and. (x <= spline_obj%k(cell+1))) then
       print *, 'ERROR. b_splines_at_x(): the given value of x is not ', &
            'inside the specified cell.'
       STOP
    end if
    deg = spline_obj%degree
    num_pts = spline_obj%num_pts

    ! for OPEN BC's, in a similar way that the end nodes are replicated beyond
    ! the domain to carry out the calculations, for the derivatives we also 
    ! extend the values of the cell spacings at both ends
    delta_left  = spline_obj%k(2)       - spline_obj%k(1)
    delta_right = spline_obj%k(num_pts) - spline_obj%k(num_pts-1)

    ! FIXME, MORTAL SIN: HERE WE HAVE DUPLICATED CODE WITH THE PREVIOUS
    ! FUNCTION AND EXPECT TO DUPLICATE THIS EVEN MORE WITH A FUNCTION
    ! THAT FURTHER COMBINES VALUES AND DERIVATIVES. THIS IS NOT ACCEPTABLE.
    ! EVENTUALLY THIS CODE SEGMENT SHOULD BE MACROIZED. LEAVE AS IS FOR THE
    ! MOMENT SINCE WE NEED TO FIX THE ISSUE OF POSSIBLY ZERO VALUES IN
    ! DENOMINATORS AT LEAST. PRODUCTION VERSION SHOULD NOT HAVE DUPLICATED
    ! CODE IN THIS CASE.

    ! Build the zeroth-order splines. The middle cell of the splines array
    ! corresponds to the 'cell' given as argument. So, 
    ! splines(deg+1) = 1.0.
    splines(:)     = 0.0
    splines(deg+1) = 1.0

    ! Build the higher order splines.
    last = 2*deg  
    do j=1,deg-1 ! we stop earlier to compute derivatives
       do i=1,last
          current    = cell - deg + i - 1
          ti         = spline_obj%k(current)
          tip1       = spline_obj%k(current+1)
          tipj       = spline_obj%k(current+j)
          tipjp1     = spline_obj%k(current+j+1)
          ! This is a dangerous situation for which we need some sort of
          ! protection: What guarantees are there that these denominators
          ! will not be zero?? This should probably be error-checked, else
          ! one can just end up with an array of NaN's.
          if(tipj==ti)then
            !print *,(x-ti)*splines(i)
            fac1  = 0._f64
          else
            !print *,(tipjp1 - x)*splines(i+1)
            fac1       = (x - ti)/(tipj - ti)
          endif
          !fac1       = (x - ti)/(tipj - ti)
          if(tipjp1==tip1)then
            !print *,(tipjp1 - x)*splines(i+1)
            fac2  = 0._f64
          else          
            fac2       = (tipjp1 - x)/(tipjp1 - tip1)
          endif  
          ! AGAIN: Super-ugly step to eliminate those cases where the 
          ! denominator is zero but the spline value is also zero, thus 
          ! avoiding the indeterminate result ( NaN's ). 
          if( splines(i) .ne. 0.0 ) then
             term1 = fac1*splines(i)
          else
             term1 = 0.0_f64
          end if
          if( splines(i+1) .ne. 0.0 ) then
             term2 = fac2*splines(i+1)
          else
             term2 = 0.0_f64
          end if
          splines(i) = term1 + term2
!          splines(i) = fac1*splines(i) + fac2*splines(i+1)
       end do
       last = last - 1
    end do
    ! At this moment we have an array with values of the splines up to the
    ! order spline_degree - 1. Proceed to compute the derivatives of order
    ! spline_degree.
    if(deg>0)then
      do i=1,last
       current = cell - deg + i - 1
       delta1 = spline_obj%k(current+deg) - spline_obj%k(current)
       delta2 = spline_obj%k(current+deg+1) - spline_obj%k(current+1)
       call check_if_delta_is_equal_to_zero(delta1, current, spline_obj)
       call check_if_delta_is_equal_to_zero(delta2, current, spline_obj)
       derivs(i) = ( deg*splines(i)/delta1 - deg*splines(i+1)/delta2 )
       !delta_x = spline_obj%k(current+1) - spline_obj%k(current)
       !if( delta_x == 0.0 ) then
       !   if( current .le. 1 ) then
       !      delta_x = delta_left
       !   else if( current .ge. spline_obj%num_pts ) then
       !      delta_x = delta_right
       !   end if
       !end if
       !derivs(i) = (splines(i) - splines(i+1))/delta_x
      end do
    endif
    b_spline_derivatives_at_x(1:deg+1) = derivs(1:deg+1)

  end function b_spline_derivatives_at_x

  subroutine check_if_delta_is_equal_to_zero(delta, current, spline_obj)

    sll_real64, intent(inout)                 :: delta
    sll_int32, intent(in)                     :: current
    sll_int32                                 :: num_pts
    type(arbitrary_degree_spline_1d), pointer :: spline_obj
    sll_real64                                :: delta_left, delta_right

    num_pts = spline_obj%num_pts
    delta_left  = spline_obj%k(2)       - spline_obj%k(1)
    delta_right = spline_obj%k(num_pts) - spline_obj%k(num_pts-1)

    if ( delta == 0.0 ) then
       if ( current .le. 1 ) then
          delta = delta_left
       else if( current .ge. num_pts ) then
          delta = delta_right
       end if
    end if

  end subroutine check_if_delta_is_equal_to_zero

  function b_splines_and_derivs_at_x( spline_obj, cell, x )
    type(arbitrary_degree_spline_1d), pointer      :: spline_obj
    sll_int32, intent(in)                          :: cell
    sll_real64, intent(in)                         :: x
    sll_int32                                      :: deg
    sll_real64, dimension(2,1:spline_obj%degree+1) :: b_splines_and_derivs_at_x
    sll_real64, dimension(1:2*spline_obj%degree+1) :: derivs
    sll_real64, dimension(1:2*spline_obj%degree+1) :: splines 
    sll_int32                                      :: i
    sll_int32                                      :: j
    sll_int32                                      :: last
    sll_real64                                     :: ti     ! t(i)
    sll_real64                                     :: tip1   ! t(i+1)
    sll_real64                                     :: tipj   ! t(i+j)
    sll_real64                                     :: tipjp1 ! t(i+j+1)
    sll_real64                                     :: fac1
    sll_real64                                     :: fac2
    sll_real64                                     :: term1
    sll_real64                                     :: term2
    sll_real64                                     :: delta1, delta2
    sll_real64                                     :: delta_left
    sll_real64                                     :: delta_right
    sll_int32                                      :: current
    sll_int32                                      :: num_pts

    ! Run some checks on the arguments.
    SLL_ASSERT(associated(spline_obj))
    SLL_ASSERT(x >= spline_obj%xmin)
    SLL_ASSERT(x <= spline_obj%xmax)
    SLL_ASSERT(cell >= 1)
    SLL_ASSERT(cell <= spline_obj%num_pts - 1)
    ! This is checked always.
    if( .not. (x >= spline_obj%k(cell)) .and. (x <= spline_obj%k(cell+1))) then
       print *, 'ERROR. b_splines_at_x(): the given value of x is not ', &
            'inside the specified cell.'
       STOP
    end if

    deg = spline_obj%degree
    num_pts = spline_obj%num_pts

    ! for OPEN BC's, in a similar way that the end nodes are replicated beyond
    ! the domain to carry out the calculations, for the derivatives we also 
    ! extend the values of the cell spacings at both ends
    delta_left  = spline_obj%k(2)       - spline_obj%k(1)
    delta_right = spline_obj%k(num_pts) - spline_obj%k(num_pts-1)

    ! FIXME, MORTAL SIN: HERE WE HAVE DUPLICATED CODE WITH THE PREVIOUS
    ! FUNCTION AND EXPECT TO DUPLICATE THIS EVEN MORE WITH A FUNCTION
    ! THAT FURTHER COMBINES VALUES AND DERIVATIVES. THIS IS NOT ACCEPTABLE.
    ! EVENTUALLY THIS CODE SEGMENT SHOULD BE MACROIZED. LEAVE AS IS FOR THE
    ! MOMENT SINCE WE NEED TO FIX THE ISSUE OF POSSIBLY ZERO VALUES IN
    ! DENOMINATORS AT LEAST. PRODUCTION VERSION SHOULD NOT HAVE DUPLICATED
    ! CODE IN THIS CASE.

    ! Build the zeroth-order splines. The middle cell of the splines array
    ! corresponds to the 'cell' given as argument. So, 
    ! splines(deg+1) = 1.0.
    splines(:)     = 0.0
    splines(deg+1) = 1.0

    ! Build the higher order splines.
    last = 2*deg  
    do j=1,deg-1 ! we stop earlier to compute derivatives
       do i=1,last
          current    = cell - deg + i - 1
          ti         = spline_obj%k(current)
          tip1       = spline_obj%k(current+1)
          tipj       = spline_obj%k(current+j)
          tipjp1     = spline_obj%k(current+j+1)
          ! This is a dangerous situation for which we need some sort of
          ! protection: What guarantees are there that these denominators
          ! will not be zero?? This should probably be error-checked, else
          ! one can just end up with an array of NaN's.
          if(tipj==ti)then
          else
            fac1       = (x - ti)/(tipj - ti)
          endif
          if(tipjp1==tip1)then
            fac2 = 0._f64
          else
            fac2       = (tipjp1 - x)/(tipjp1 - tip1)
          endif  
          ! AGAIN: Super-ugly step to eliminate those cases where the 
          ! denominator is zero but the spline value is also zero, thus 
          ! avoiding the indeterminate result ( NaN's ). 
          if( splines(i) .ne. 0.0 ) then
             term1 = fac1*splines(i)
          else
             term1 = 0.0_f64
          end if
          if( splines(i+1) .ne. 0.0 ) then
             term2 = fac2*splines(i+1)
          else
             term2 = 0.0_f64
          end if
          splines(i) = term1 + term2
!          splines(i) = fac1*splines(i) + fac2*splines(i+1)
       end do
       last = last - 1
    end do
    print *, 'array to compute derivatives: ',splines(1:last)
    ! At this moment we have an array with values of the splines up to the
    ! order spline_degree - 1. Proceed to compute the derivatives of order
    ! spline_degree.
    if(deg>0)then
      do i=1,last
       current = cell - deg + i - 1
       delta1 = spline_obj%k(current+deg) - spline_obj%k(current)
       delta2 = spline_obj%k(current+deg+1) - spline_obj%k(current+1) 
       derivs(i) = ( deg*splines(i)/delta1 - deg*splines(i+1)/delta2 )
       call check_if_delta_is_equal_to_zero(delta1, current, spline_obj)
       call check_if_delta_is_equal_to_zero(delta2, current, spline_obj)
       !delta_x = spline_obj%k(current+1) - spline_obj%k(current)
       !if( delta_x == 0.0 ) then
       !   if( current .le. 1 ) then
       !      delta_x = delta_left
       !   else if( current .ge. spline_obj%num_pts ) then
       !      delta_x = delta_right
       !   end if
       !end if
       !derivs(i) = (splines(i) - splines(i+1))/delta_x
      end do
    endif
    b_splines_and_derivs_at_x(2,1:deg+1) = derivs(1:deg+1)
    ! Finish computing the splines for the desired degree
    if(deg>0)then
      j = deg
      do i=1,last
       current    = cell - deg + i - 1
       ti         = spline_obj%k(current)
       tip1       = spline_obj%k(current+1)
       tipj       = spline_obj%k(current+j)
       tipjp1     = spline_obj%k(current+j+1)
       ! This is a dangerous situation for which we need some sort of
       ! protection: What guarantees are there that these denominators
       ! will not be zero?? This should probably be error-checked, else
       ! one can just end up with an array of NaN's.
       if(tipj==ti)then
         fac1 = 0._f64
       else
         fac1       = (x - ti)/(tipj - ti)
       endif
       !fac1       = (x - ti)/(tipj - ti)
       if(tipjp1==tip1)then
         fac2 = 0._f64
       else
         fac2       = (tipjp1 - x)/(tipjp1 - tip1)
       endif  
          ! AGAIN: Super-ugly step to eliminate those cases where the 
          ! denominator is zero but the spline value is also zero, thus 
          ! avoiding the indeterminate result ( NaN's ). 
          if( splines(i) .ne. 0.0 ) then
             term1 = fac1*splines(i)
          else
             term1 = 0.0_f64
          end if
          if( splines(i+1) .ne. 0.0 ) then
             term2 = fac2*splines(i+1)
          else
             term2 = 0.0_f64
          end if
          splines(i) = term1 + term2
!       splines(i) = fac1*splines(i) + fac2*splines(i+1)
      end do
    endif  
    b_splines_and_derivs_at_x(1,1:deg+1) = splines(1:deg+1)
  end function b_splines_and_derivs_at_x

  subroutine delete_arbitrary_order_spline_1d( spline )
    type(arbitrary_degree_spline_1d), pointer :: spline
    sll_int32                    :: ierr
    if( .not. associated(spline) ) then
       print *, 'ERROR. delete_arbitrary_order_spline_1d(): given spline ', &
            'pointer is not associated.'
       STOP
    end if
    SLL_DEALLOCATE( spline%k, ierr )
    SLL_DEALLOCATE( spline, ierr )
  end subroutine delete_arbitrary_order_spline_1d





  ! *************************************************************************
  !
  !                    UNIFORM B-SPLINE FUNCTIONS
  !
  ! *************************************************************************

  !> returns an array with the values of the b-splines of the 
  !> requested degree, evaluated at a given cell offset. The cell size is
  !> normalized between 0 and 1, thus the offset given must be a number
  !> between 0 and 1.

  function uniform_b_splines_at_x( spline_degree, normalized_offset ) result(out)
    !new uniform_b_splines_at_x
    !out(1:d+1)= B_d(-(d+1)/2+d+x),...,B_d(-(d+1)/2+x) with d=spline_degree and x=normalized_offset
    !where B_d=B_{d-1}*B_0 and B_0=1_[-1/2,1/2] and * is convolution
    !complexity is lower than deboor splines in this uniform setting (translation of a unique B-spline)
    !the following code can be used for comparison with deboor
    !do i=-d,d+1
    !t(i+d+1)=real(i,8)
    !enddo
    !call bsplvb(t,d+1,1,normalized_offset,d+1,out)
    !we also have the property (from the symmetry of the B-spline)
    !out(1:d+1)= B_d(-(d+1)/2+xx),...,B_d(-(d+1)/2+d+xx),..., where xx=1-normalized_offset
    
    implicit none
    sll_int32, intent(in)                      :: spline_degree
    sll_real64, intent(in)                     :: normalized_offset
    sll_real64, dimension(1:spline_degree+1)   :: out
    sll_int32                                  :: degree,i
    sll_real64                                 :: degree_real,i_real
    sll_real64                                 :: r_degree
    sll_real64                                 :: fac1,fac2
    sll_real64                                 :: x,xx,tmp1,tmp2
    !SLL_ASSERT( spline_degree >= 0 )
    !SLL_ASSERT( normalized_offset >= 0.0_f64 )
    !SLL_ASSERT( normalized_offset <= 1.0_f64 )
    out(1) = 1._f64
    xx=normalized_offset
    x=1._f64-normalized_offset
    !exchange x and xx for reverse
    do degree=1,spline_degree
      degree_real = real(degree,f64)
      r_degree    = 1._f64/degree_real
      fac1        = x*r_degree
      tmp1        = out(1)
      out(1)      = fac1*out(1)
      do i=1,degree-1
        i_real = real(i,f64)
        fac1   = (x+i_real)*r_degree
        fac2   = (xx+degree_real-i_real)*r_degree        
        tmp2   = fac1*out(i+1)+fac2*tmp1
        tmp1   = out(i+1)
        out(i+1) = tmp2
      enddo
      fac2            = xx*r_degree 
      out(degree+1)     = fac2*tmp1
    enddo
  end function uniform_b_splines_at_x





  function uniform_b_splines_at_x_old( spline_degree, normalized_offset )
    sll_int32, intent(in)                      :: spline_degree
    sll_real64, dimension(1:spline_degree+1)   :: uniform_b_splines_at_x_old
    sll_real64, intent(in)                     :: normalized_offset
    sll_real64, dimension(1:2*spline_degree+1) :: splines
    sll_int32                                  :: i
    sll_int32                                  :: j
    sll_int32                                  :: last
    sll_real64                                 :: jreal
    sll_real64                                 :: r_jreal
    sll_real64                                 :: fac1
    sll_real64                                 :: fac2
    sll_real64                                 :: temp

    SLL_ASSERT( spline_degree >= 0 )
    SLL_ASSERT( normalized_offset >= 0.0_f64 )
    SLL_ASSERT( normalized_offset <= 1.0_f64 )

    ! Build the zeroth-order splines. The middle cell of the splines array
    ! corresponds to the 'cell' that contains 'x'. So for example, if a cubic
    ! spline is requested, the zeroth-order splines will be:
    !
    !        0, 0, 0, 1, 0, 0, 0
    !
    ! And from this we recursively build the higher degree splines.
    splines(:)               = 0.0_f64
    splines(spline_degree+1) = 1.0_f64
    
    ! Build the higher order splines. 
    last = 2*spline_degree  
    do j=1,spline_degree
       jreal      = real(j,f64)
       r_jreal    = 1.0_f64/jreal
       do i=1,last
          temp       = real(spline_degree - i, f64)
          fac1       = (temp + normalized_offset + 1.0_f64)*r_jreal
          fac2       = (-temp + jreal - normalized_offset)*r_jreal
          splines(i) = fac1*splines(i) + fac2*splines(i+1)
       end do
       last = last - 1
    end do
    uniform_b_splines_at_x_old(1:spline_degree+1) = splines(1:spline_degree+1)
  end function uniform_b_splines_at_x_old





  !> returns an array with the values of the b-spline derivatives of the 
  !> requested degree, evaluated at a given cell offset. The cell size is
  !> normalized between 0 and 1, hence the results must be divided by the
  !> real cell size to scale back the results.
  function uniform_b_spline_derivatives_at_x( spline_degree, normalized_offset )
    sll_int32, intent(in)                      :: spline_degree
    sll_real64, dimension(1:spline_degree+1):: uniform_b_spline_derivatives_at_x
    sll_real64, intent(in)                     :: normalized_offset
    sll_real64, dimension(1:2*spline_degree+1) :: splines
    sll_real64, dimension(spline_degree+1)     :: derivs
    sll_int32                                  :: i
    sll_int32                                  :: j
    sll_int32                                  :: last
    sll_real64                                 :: jreal
    sll_real64                                 :: r_jreal
    sll_real64                                 :: fac1
    sll_real64                                 :: fac2
    sll_real64                                 :: temp

    SLL_ASSERT( spline_degree >= 0 )
    SLL_ASSERT( normalized_offset >= 0.0_f64 )
    SLL_ASSERT( normalized_offset <= 1.0_f64 )

    ! Build the zeroth-order splines. The middle cell of the splines array
    ! corresponds to the 'cell' that contains 'x'. So for example, if a cubic
    ! spline is requested, the zeroth-order splines will be:
    !
    !        0, 0, 0, 1, 0, 0, 0
    !
    ! And from this we recursively build the higher degree splines.
    splines(:)               = 0.0_f64
    splines(spline_degree+1) = 1.0_f64
    
    ! Build the higher order splines. 
    last = 2*spline_degree  
    do j=1,spline_degree - 1 ! we stop earlier to compute the derivatives
       jreal      = real(j,f64)
       r_jreal    = 1.0_f64/jreal
       do i=1,last
          temp       = real(spline_degree - i, f64)
          fac1       = (temp + normalized_offset + 1.0_f64)*r_jreal
          fac2       = (-temp + jreal - normalized_offset)*r_jreal
          splines(i) = fac1*splines(i) + fac2*splines(i+1)
       end do
       last = last - 1
    end do
    ! check: splines must add to 1.0
    !     print *,  'sum = ', sum(splines(1:spline_degree+1))
    ! At this moment, we have an array with values of the splines up to
    ! the degree 'spline_degree -1'. We proceed to compute the derivatives
    ! of degree 'spline_degree'.
    ! print *, 'splines for derivatives:', splines(:)
    derivs(:) = 0._f64
    do i=1,last
       derivs(i) = splines(i) - splines(i+1)
    end do
    uniform_b_spline_derivatives_at_x(1:spline_degree+1) = &
         derivs(1:spline_degree+1)
  end function uniform_b_spline_derivatives_at_x

  !> returns an array with the values of the b-spline derivatives of the 
  !> requested degree, evaluated at a given cell offset. The cell size is
  !> normalized between 0 and 1, hence the results must be divided by the
  !> real cell size to scale back the results.
  function uniform_b_splines_and_derivs_at_x( degree, normalized_offset )
    sll_int32, intent(in)                 :: degree
    sll_real64, dimension(2,1:degree+1)   :: uniform_b_splines_and_derivs_at_x
    sll_real64, intent(in)                :: normalized_offset
    sll_real64, dimension(1:2*degree+1)   :: splines
    sll_real64, dimension(degree+1)       :: derivs
    sll_int32                             :: i
    sll_int32                             :: j
    sll_int32                             :: last
    sll_real64                            :: jreal
    sll_real64                            :: r_jreal
    sll_real64                            :: fac1
    sll_real64                            :: fac2
    sll_real64                            :: temp

    SLL_ASSERT( degree >= 0 )
    SLL_ASSERT( normalized_offset >= 0.0_f64 )
    SLL_ASSERT( normalized_offset <= 1.0_f64 )

    ! Build the zeroth-order splines. The middle cell of the splines array
    ! corresponds to the 'cell' that contains 'x'. So for example, if a cubic
    ! spline is requested, the zeroth-order splines will be:
    !
    !        0, 0, 0, 1, 0, 0, 0
    !
    ! And from this we recursively build the higher degree splines.
    splines(:)               = 0.0_f64
    splines(degree+1) = 1.0_f64
    
    ! Build the higher order splines. 
    last = 2*degree  
    do j=1,degree - 1 ! we stop earlier to compute the derivatives
       jreal      = real(j,f64)
       r_jreal    = 1.0_f64/jreal
       do i=1,last
          temp       = real(degree - i, f64)
          fac1       = (temp + normalized_offset + 1.0_f64)*r_jreal
          fac2       = (-temp + jreal - normalized_offset)*r_jreal
          splines(i) = fac1*splines(i) + fac2*splines(i+1)
       end do
       last = last - 1
    end do
    ! At this moment, we have an array with values of the splines up to
    ! the degree 'degree -1'. We proceed to compute the derivatives
    ! of degree 'degree'.
    ! print *, 'splines for derivatives:', splines(:)
    do i=1,last
       derivs(i) = splines(i) - splines(i+1)
    end do
    ! Fill out the derivatives in the output array
    uniform_b_splines_and_derivs_at_x(2,1:degree+1) = derivs(1:degree+1)

    ! Do the last pass to complete the splines of the desired degree
    if(degree>0)then
      j = degree
      jreal      = real(j,f64)
      r_jreal    = 1.0_f64/jreal
      do i=1,last
       temp       = real(degree - i, f64)
       fac1       = (temp + normalized_offset + 1.0_f64)*r_jreal
       fac2       = (-temp + jreal - normalized_offset)*r_jreal
       splines(i) = fac1*splines(i) + fac2*splines(i+1)
      end do
    endif
    
    uniform_b_splines_and_derivs_at_x(1,1:degree+1) = splines(1:degree+1)
  end function uniform_b_splines_and_derivs_at_x






end module sll_arbitrary_degree_splines
