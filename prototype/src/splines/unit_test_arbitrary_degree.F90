program arbitrary_degree_spline_tester

#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"

  use sll_arbitrary_degree_splines

  implicit none

  logical                                :: passed_test
  passed_test = .true.

  print *, '*****************************************************************'
  print *, ' Testing arbitrary degree splines module: '
  print *, '*****************************************************************'

  print *, '*****************************************************************'
  print *, ' non uniform periodic '
  print *, '*****************************************************************'
  call test_nonuniform_arb_deg_splines_periodic( passed_test )
  print *, '*****************************************************************'
  print *, ' uniform B-splines randomly '
  print *, '*****************************************************************'
  call test_uniform_b_splines_randomly( passed_test )
  print *, '*****************************************************************'
  print *, ' non uniform open '
  print *, '*****************************************************************'
  call test_nonuniform_arb_deg_splines_open( passed_test )

  if (passed_test) then
     print *, 'PASSED'
  else
     print *, 'FAILED'
  end if

contains

  subroutine test_nonuniform_arb_deg_splines_periodic( passed_test )

    logical, intent(inout) :: passed_test
    sll_real64, dimension(:), allocatable :: knots
    sll_int32  :: i,j
    sll_int32  :: num_pts
    sll_int32  :: degree
    sll_real64 :: min_val
    sll_int32  :: ierr
    sll_real64 :: rnd
    sll_real64 :: step
    sll_int32  :: cell
    sll_real64 :: x
    sll_real64 :: acc, acc2
    sll_real64 :: criterion
    sll_int32  :: num_tests
    sll_real64, dimension(:), allocatable :: answer1
    sll_real64, dimension(:), allocatable :: answer2
    sll_real64, dimension(:,:), allocatable :: answer3
    type(arbitrary_degree_spline_1d), pointer :: spline

    num_tests = 10
    criterion = 1.0e-15
    degree    = 3
    min_val   = 0.0
    num_pts   = 10
    step      = 1.0
    SLL_ALLOCATE(knots(num_pts),ierr)
    SLL_ALLOCATE(answer1(degree+1),ierr)
    SLL_ALLOCATE(answer2(degree+1),ierr)
    SLL_ALLOCATE(answer3(2,degree+1),ierr)

    ! fill knots array. Try first a uniform set of knots to compare with the
    ! uniform spline functions.
    knots(1) = min_val
    do i=2,num_pts
       call random_number(rnd)
       knots(i) = knots(i-1) + rnd !step
    end do
    print *, ' '
    print *, 'knots array = ', knots(:)
    print *, ' '
    ! fill spline object
    spline => new_arbitrary_degree_spline_1d( &
         degree, &
         knots, &
         num_pts, &
         PERIODIC_ARBITRARY_DEG_SPLINE )

    do j=1,num_tests
       call random_number(rnd)
       x = min_val + rnd*(knots(num_pts)-min_val)
       cell = find_index( x, knots, num_pts )
       acc = 0.0_f64
       
       ! test spline values
       answer1(:) = b_splines_at_x(spline, cell, x)
       acc = abs(1.0_f64 - sum(answer1(1:degree+1)))
       passed_test = passed_test .and. (acc < criterion)
       if( passed_test .eqv. .false. ) then
          print *, 'nonuniform splines test failure, spline values case:'
          print *, 'cell = ', cell, 'x = ', x
          print *, 'nonuniform: ', answer1(:)
          print *, 'accumulator = ', acc
          print*, 'Exiting...'
          stop
       end if
       
       ! test spline derivatives
       acc = 0.0_f64
       answer1(:) = b_spline_derivatives_at_x(spline, cell, x)
       acc = abs(sum(answer1(1:degree+1)))
       passed_test = passed_test .and. (acc < criterion)
       if( passed_test .eqv. .false. ) then
          print *, 'periodic nonuniform splines test failure, derivatives case:'
          print *, 'cell = ', cell, 'x = ', x
          print *, 'nonuniform: ', answer1(:)
          print *, 'accumulator = ', acc
          print*, 'Exiting...'
          stop
       end if
       
       ! test values and derivatives
       acc  = 0.0_f64
       acc2 = 0.0_f64
       answer3(:,:) = b_splines_and_derivs_at_x(spline, cell, x)
       acc  = acc +  abs(sum(answer3(1,1:degree+1)))
       acc2 = acc2 + abs(sum(answer3(2,1:degree+1)))
       
       passed_test = passed_test .and. ((1.0_f64 - acc) < criterion)
       if( passed_test .eqv. .false. ) then
          print *, 'periodic nonuniform splines test failure, ',&
               'values and derivatives case:'
          print *, 'cell = ', cell, 'x = ', x
          print *, 'nonuniform: ', answer3(:,:)
          print *, 'uniform: ', answer1(:)
          print *, 'accumulator = ', acc
          print*, 'Exiting...'
          stop
       end if
       passed_test = passed_test .and. (acc2 < criterion)
       if( passed_test .eqv. .false. ) then
          print *, 'periodic nonuniform splines test failure, ',&
               'values and derivatives case:'
          print *, 'cell = ', cell, 'x = ', x
          print *, 'nonuniform: ', answer3(:,:)
          print *, 'accumulator = ', acc
          print*, 'Exiting...'
          stop
       end if

    end do

    SLL_DEALLOCATE_ARRAY(answer1, ierr)
    SLL_DEALLOCATE_ARRAY(answer2, ierr)
    SLL_DEALLOCATE_ARRAY(answer3, ierr)

    call delete(spline)
  end subroutine test_nonuniform_arb_deg_splines_periodic

  subroutine test_uniform_b_splines_randomly( passed_flag )
    logical, intent(inout)      :: passed_flag
    sll_real64                  :: criterion
    sll_real64                  :: argument
    sll_real64                  :: argument_copy
    sll_real64, dimension(:), allocatable   :: results
    sll_real64, dimension(:), allocatable   :: derivatives
    sll_real64, dimension(:,:), allocatable :: sp_and_derivs
    sll_int32                   :: num_tests
    sll_int32                   :: i
    sll_int32                   :: j
    sll_int32                   :: max_degree
    sll_int32                   :: ierr

    criterion          = 1.0e-14
    argument           = 0.0_f64
    num_tests          = 10000
    argument_copy      = argument
    max_degree         = 12
 
    do j=0, max_degree
       do i=1,num_tests
          SLL_CLEAR_ALLOCATE(results(1:j+1),ierr)
          SLL_CLEAR_ALLOCATE(derivatives(1:j+1),ierr)
          SLL_CLEAR_ALLOCATE(sp_and_derivs(1:2,1:j+1),ierr)
          call random_number(argument)

          results(:) = uniform_b_splines_at_x(j, argument)
          derivatives(:) = uniform_b_spline_derivatives_at_x(j, argument)
          sp_and_derivs(:,:) = uniform_b_splines_and_derivs_at_x(j, argument)

          passed_flag = passed_flag .and. (abs(sum(results)-1.0_f64).le.criterion)

          if( passed_flag .eqv. .false. ) then
             print *, 'Test_uniform_b_splines, values: wrong result for x = ', &
                  argument
             print *, 'Degree = ', j, 'Reduction = ', sum(results)
             print*, 'Exiting...'
             stop
          end if

          passed_flag = passed_flag .and. &
               (abs(sum(derivatives)-0.0_f64).le.criterion)
          if( passed_flag .eqv. .false. ) then
             print *, 'Test_uniform_b_splines, derivs: wrong result for x = ', &
                  argument
             print *, 'Degree = ', j, 'Reduction = ', sum(derivatives)
             print*, 'Exiting...'
             stop
          end if

          passed_flag = passed_flag .and. &
               (abs(sum( sp_and_derivs(1,:) ) - 1.0_f64).le.criterion) .and. &
               (abs(sum( sp_and_derivs(2,:) ) - 0.0_f64).le.criterion)
          if( passed_flag .eqv. .false. ) then
             print *, 'Test_uniform_b_splines, vals and derivs: ', &
                  'wrong result for x = ', argument
             print *, 'Degree = ', j, 'Reduction = ', sum(sp_and_derivs(1,:)) &
                                                    , sum(sp_and_derivs(2,:))
             print*, 'Exiting...'
             stop
          end if
          results(:)         = 0.0_f64
          derivatives(:)     = 0.0_f64
          sp_and_derivs(:,:) = 0.0_f64

          SLL_DEALLOCATE_ARRAY(results, ierr)
          SLL_DEALLOCATE_ARRAY(derivatives, ierr)
          SLL_DEALLOCATE_ARRAY(sp_and_derivs, ierr)
       end do
    end do
  end subroutine test_uniform_b_splines_randomly


  ! The case of 'open' boundary condition yields spline values different
  ! than in the 'periodic' case. Since we can not compare with the uniform
  ! splines anymore to get the right answer, we need a different criterion.
  ! For lack of something better, at this moment we only check that the 
  ! different spline values, when added, will equal 1.0.
  subroutine test_nonuniform_arb_deg_splines_open( passed_test )
    logical, intent(inout) :: passed_test
    sll_real64, dimension(:), allocatable :: knots
    sll_int32  :: i,j
    sll_int32  :: num_pts
    sll_int32  :: degree
    sll_real64 :: min_val
    sll_int32  :: ierr
    sll_real64 :: rnd
    sll_real64 :: step
    sll_int32  :: cell
    sll_real64 :: x
    sll_real64 :: acc, acc2
    sll_real64 :: criterion
    sll_int32  :: num_tests
    sll_real64, dimension(:), allocatable     :: answer
    sll_real64, dimension(:,:), allocatable   :: answer2
    type(arbitrary_degree_spline_1d), pointer :: spline

    num_tests = 10 !100000
    criterion = 1.0e-15
    degree  = 3
    min_val = 0.0
    num_pts = 10
    step    = 1.0
    SLL_ALLOCATE(knots(num_pts),ierr)
    SLL_ALLOCATE(answer(degree+1),ierr)
    SLL_ALLOCATE(answer2(2,degree+1),ierr)

    ! fill knots array. Try first a uniform set of knots to compare with the
    ! uniform spline functions.
    knots(1) = min_val
    do i=2,num_pts
       call random_number(rnd)
       knots(i) = knots(i-1) + step + rnd
    end do
    print *, 'knots array = ', knots(:)

    ! fill spline object
    spline => new_arbitrary_degree_spline_1d( &
         degree, &
         knots, &
         num_pts, &
         OPEN_ARBITRARY_DEG_SPLINE )

    do j=1,num_tests
       call random_number(rnd)
       x = min_val + rnd*(knots(num_pts)-min_val)
       cell = find_index( x, knots, num_pts )
       acc = 0.0_f64
       
       ! test spline values
       answer(:) = b_splines_at_x(spline, cell, x)
       acc = sum(answer(1:degree+1))
       passed_test = passed_test .and. (abs(1.0_f64 - acc) < criterion)
       if( passed_test .eqv. .false. ) then
          print *, 'nonuniform splines test failure:'
          print *, 'cell = ', cell, 'x = ', x
          print *, 'nonuniform: ', answer(:)
          print *, 'accumulator = ', acc
          print *, 'Exiting...'
          stop
       end if
       
       ! test spline derivatives
       acc = 0.0_f64
       answer(:) = b_spline_derivatives_at_x(spline, cell, x)
       acc = sum(answer(1:degree+1))
       passed_test = passed_test .and. (abs(acc) < criterion)
       if( passed_test .eqv. .false. ) then
          print *, 'nonuniform splines test failure, derivatives case:'
          print *, 'cell = ', cell, 'x = ', x
          print *, 'derivatives: ', answer(:)
          print *, 'accumulator = ', acc
          print*, 'Exiting...'
          stop
       end if
       
       ! test values and derivatives
       acc  = 0.0_f64
       acc2 = 0.0_f64
       answer2(:,:) = b_splines_and_derivs_at_x(spline, cell, x)
       acc  = sum(answer2(1,1:degree+1))
       acc2 = sum(answer2(2,1:degree+1))
       passed_test = passed_test .and. (abs(1.0_f64 - acc) < criterion)
       passed_test = passed_test .and. (abs(acc2) < criterion)
       if( passed_test .eqv. .false. ) then
          print *, 'nonuniform splines test failure, ',&
               'values and derivatives case:'
          print *, 'cell = ', cell, 'x = ', x
          print *, 'splines: ', answer2(1,:)
          print *, 'derivatives: ', answer2(2,:)
          print *, 'sum of splines = ', acc
          print *, 'sum of derivatives = ', acc2
          print*, 'Exiting...'
          stop
       end if
    end do
    SLL_DEALLOCATE_ARRAY(answer, ierr)
    SLL_DEALLOCATE_ARRAY(answer2, ierr)
    call delete(spline)
  end subroutine test_nonuniform_arb_deg_splines_open

  ! given an array of values ordered in ascending order and a value 'x', 
  ! find index returns the index 'i' of the array such that:
  !                  array(i) <= x <= array(i+1)
  ! fortran-indexing is assumed. If x is not between array(1) and array(n),
  ! where 'n' is the array length, then the value -1 is returned.
  function find_index( x, array, array_length  )
    sll_int32 :: find_index
    sll_real64, intent(in) :: x
    sll_real64, intent(in), dimension(:) :: array
    sll_int32, intent(in)  :: array_length
    sll_int32 :: i

    if( x < array(1) ) then
       find_index = -1
       return
    end if

    do i=1,array_length-1
       if( (array(i) <= x) .and. (x <= array(i+1)) ) then
          find_index = i
          return
       end if
    end do

    if( x > array(array_length) ) then
       find_index = -1
       return
    end if
  end function find_index

end program arbitrary_degree_spline_tester
