program test_arbitrary_degree_splines
!-----------------
! Contact: Eric Sonnendrucker

#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"

  use sll_m_arbitrary_degree_splines
  use sll_m_timer

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
  print *, '*****************************************************************'
  print *, ' test CPU time '
  print *, '*****************************************************************'
  call test_cpu_time

  if (passed_test) then
     print *, 'PASSED'
  else
     print *, 'FAILED'
  end if

contains

  subroutine test_nonuniform_arb_deg_splines_periodic( passed_test )

    logical, intent(inout) :: passed_test
    sll_real64, dimension(:), allocatable :: grid
    sll_int32  :: i,j
    sll_int32  :: num_pts
    sll_int32  :: degree
    sll_real64 :: min_val
    sll_int32  :: ierr
    sll_real64 :: rnd
    sll_int32  :: cell
    sll_real64 :: x
    sll_real64 :: acc, acc2
    sll_real64 :: criterion
    sll_int32  :: num_tests
    sll_real64, dimension(:), allocatable :: x_test
    sll_real64, dimension(:,:), allocatable :: expected1
    sll_real64, dimension(:,:), allocatable :: expected2
    sll_real64, dimension(:), allocatable :: answer1
    sll_real64, dimension(:), allocatable :: answer2
    sll_real64, dimension(:,:), allocatable :: answer3
    type(arbitrary_degree_spline_1d), pointer :: spline

    ! Test on random grid for random degree
    num_tests = 10
    criterion = 1.0d-15
    min_val   = 0.0_f64
    num_pts   = 10

    degree    = 6
    print *, " ------------- Degree = ", degree, " -----------------"

    SLL_ALLOCATE(grid(num_pts),ierr)
    SLL_ALLOCATE(answer1(degree+1),ierr)
    SLL_ALLOCATE(answer2(degree+1),ierr)
    SLL_ALLOCATE(answer3(2,degree+1),ierr)

    ! --------- 1D SPLINE INITIALIZATION ON NON UNIFORM MESH ----
    ! Creating non uniform mesh....
    grid(1) = min_val
    do i=2,num_pts
       call random_number(rnd)
       grid(i) = grid(i-1) + rnd !step
    end do
    ! ..... non uniform mesh done
    ! creating non uniform 1d spline 
    spline => new_arbitrary_degree_spline_1d( &
         degree, &
         grid, &
         num_pts, &
         PERIODIC_ARBITRARY_DEG_SPLINE )
    ! --------- INITIALIZATION DONE ------------

    ! To compensate random factor, we do the test more than once:
    do j=1,num_tests
       ! We compute a point randomly on the mesh:
       call random_number(rnd)
       x = min_val + rnd*(grid(num_pts)-min_val)
       ! We find its cell:
       cell = find_cell( spline, x )
       ! initialization accumulator:
       acc = 0.0_f64
       ! computing all non zero splines at point x:
       answer1(:) = b_splines_at_x(spline, cell, x)
       ! computing all non zero spline derivatives at point x:
       answer2(:) = b_spline_derivatives_at_x(spline, cell, x)
       ! computing both all non zero splines and derivatives at point x:
       answer3(:,:) = b_splines_and_derivs_at_x(spline, cell, x)
       ! testing partition of unity property of b-splines:
       acc = abs(1.0_f64 - sum(answer1(1:degree+1)))
       passed_test = passed_test .and. (acc < criterion)
       ! corresponding prints:
       if( passed_test .eqv. .false. ) then
          print *, 'nonuniform splines test failure, spline values case:'
          print *, 'cell = ', cell, 'x = ', x
          print *, 'nonuniform: ', answer1(:)
          print *, 'accumulator = ', acc
          print*, 'Exiting...'
          stop
       else
          print *, "test =", j, 'cell = ', cell, 'x=', x
          print*, "   sum of non null splines at x (expected 1) =", &
               sum(answer1(1:degree+1)), "... PASSED"
       end if

       ! test spline derivatives
       acc = 0.0_f64
       ! test of sum all derivatives (should be 0):
       acc = abs(sum(answer2(1:degree+1)))
       passed_test = passed_test .and. (acc < criterion)
       if( passed_test .eqv. .false. ) then
          print *, 'periodic nonuniform splines test failure, derivatives case:'
          print *, 'cell = ', cell, 'x = ', x
          print *, 'nonuniform: ', answer1(:)
          print *, 'accumulator = ', acc
          print*, 'Exiting...'
          stop
       else 
          print *,  "   sum of non null splines derivatives at x", &
               " (expected 0) =", sum(answer1(1:degree+1)), "... PASSED"
       end if

       ! test function b_splines_and_derivs_at_x
       ! check that spline values are the same as those from b_splines_at_x
       ! and derivatives are the same as those from b_spline_derivatives_at_x
       acc  = 0.0_f64
       acc2 = 0.0_f64
       do i=1, degree+1
          acc  = max(acc,  abs(answer1(i)-answer3(1,i)))
          acc2 = max(acc2, abs(answer2(i)-answer3(2,i)))
       end do

       passed_test = passed_test .and. (acc < criterion)
       if( passed_test .eqv. .false. ) then
          print *, 'periodic nonuniform splines test failure:'          
          print*,  'values of splines computed in b_splines_at_x and ',&
               ' b_spline_derivatives_at_x are not the same'
          print *, 'cell = ', cell, 'x = ', x
          print *, 'b_splines_at_x: ', answer1(:)
          print *, 'b_spline_and_derivs_at_x: ', answer3(1,:)
          print *, 'accumulator = ', acc
          print*, 'Exiting...'
          stop
       end if
       passed_test = passed_test .and. (acc2 < criterion)
       if( passed_test .eqv. .false. ) then
          print *, 'periodic nonuniform splines test failure:'
          print*,  'values of derivatives computed in ', &
               'b_spline_derivatives_at_x and ',&
               ' b_spline_and_derivs_at_x are not the same'
          print *, 'cell = ', cell, 'x = ', x
          print *, 'b_spline_derivatives_at_x: ', answer2(:)
          print *, 'b_spline_and_derivs_at_x: ', answer3(2,:)
          print *, 'accumulator = ', acc2
          print*, 'Exiting...'
          stop
       end if

    end do

    SLL_DEALLOCATE_ARRAY(grid,ierr)
    SLL_DEALLOCATE_ARRAY(answer1, ierr)
    SLL_DEALLOCATE_ARRAY(answer2, ierr)
    SLL_DEALLOCATE_ARRAY(answer3, ierr)

    call sll_delete(spline)

    ! Test on given grid with known answer
    criterion = 1.0d-15
    min_val   = 0.0_f64
    num_pts   = 5
    num_tests = 7
    degree    = 3
    print *, " ------------- Degree = ", degree, " -----------------"

    SLL_ALLOCATE(grid(num_pts),ierr)
    SLL_ALLOCATE(x_test(num_tests),ierr)
    SLL_ALLOCATE(expected1(degree+1,num_tests),ierr)
    SLL_ALLOCATE(expected2(degree+1,num_tests),ierr)
    SLL_ALLOCATE(answer1(degree+1),ierr)
    SLL_ALLOCATE(answer2(degree+1),ierr)
    SLL_ALLOCATE(answer3(2,degree+1),ierr)

    ! --------- 1D SPLINE INITIALIZATION ON NON UNIFORM MESH ----
    ! Define non uniform grid
    grid = (/0.0_f64, 2.0_f64, 3.0_f64, 4.5_f64, 5.0_f64 /)
    ! creating non uniform 1d spline 
    spline => new_arbitrary_degree_spline_1d( &
         degree, &
         grid, &
         num_pts, &
         PERIODIC_ARBITRARY_DEG_SPLINE )
    ! --------- INITIALIZATION DONE ------------
    ! array of points where splines will be evaluated
    x_test = (/0._f64, 0.5_f64, 1.2_f64, 2.001_f64, 3._f64, 4.0_f64, 5.0_f64 /)
    !x_test = (/0._f64, 0.5_f64, 1._f64, 2._f64, 3._f64, 4.5_f64, 5._f64 /)
    ! Expected values of bsplines and derivatives at these points
    ! at point 0.0
    expected1(1,1) = 0.4_f64
    expected1(2,1) = 0.5714285714285714_f64
    expected1(3,1) = 0.028571428571428574_f64
    expected1(4,1) = 0.0_f64
    ! at point 0.5
    expected1(1,2) = 0.16874999999999999_f64
    expected1(2,2) = 0.644345238095238_f64
    expected1(3,2) = 0.18227513227513226_f64
    expected1(4,2) = 0.004629629629629629_f64  
    ! at point 1.2
    expected1(1,3) = 0.0256_f64
    expected1(2,3) = 0.42742857142857144_f64
    expected1(3,3) = 0.4829714285714285_f64 
    expected1(4,3) = 0.06399999999999999_f64
    ! at point 2.001
    expected1(1,4) = 0.0949526665714286_f64
    expected1(2,4) = 0.6083063706285714_f64
    expected1(3,4) = 0.2967409626666666_f64 
    expected1(4,4) = 1.3333333333328926d-10
    ! at point 3.0
    expected1(1,5) = 0.2_f64
    expected1(2,5) = 0.66666666666666666_f64
    expected1(3,5) = 0.13333333333333333_f64 
    expected1(4,5) = 0.0_f64
    ! at point 4.0
    expected1(1,6) = 0.007407407407407407_f64
    expected1(2,6) = 0.25925925925925924_f64
    expected1(3,6) = 0.65_f64 
    expected1(4,6) = 0.08333333333333333_f64
    ! at point 5.0
    expected1(1,7) = 0.0_f64
    expected1(2,7) = 0.4_f64
    expected1(3,7) = 0.5714285714285714_f64
    expected1(4,7) = 0.028571428571428574_f64   
    ! derivatives
    ! at point 0.0
    expected2(1,1) = -0.6_f64
    expected2(2,1) = 0.42857142857142866_f64
    expected2(3,1) = 0.17142857142857143_f64
    expected2(4,1) = 0.0_f64
    ! at point 0.5
    expected2(1,2) = -0.3375_f64
    expected2(2,2) = -0.0982142857142857_f64
    expected2(3,2) = 0.4079365079365079_f64
    expected2(4,2) = 0.027777777777777777_f64  
    ! point 1.2
    expected2(1,3) = -0.096_f64
    expected2(2,3) = -0.44571428571428573_f64
    expected2(3,3) = 0.38171428571428573_f64 
    expected2(4,3) = 0.16_f64
    ! at point 2.001
    expected2(1,4) = -0.2851431428571429_f64
    expected2(2,4) = -0.15974525714285698_f64
    expected2(3,4) = 0.44488799999999984_f64 
    expected2(4,4) = 3.999999999999119d-07
    ! at point 3.0
    expected2(1,5) = -0.4_f64
    expected2(2,5) = 0.0_f64
    expected2(3,5) = 0.4_f64 
    expected2(4,5) = 0.0_f64
    ! at point 4.0
    expected2(1,6) = -0.04444444444444444_f64
    expected2(2,6) = -0.55555555555555555_f64
    expected2(3,6) = 0.35_f64 
    expected2(4,6) = 0.25_f64
    ! at point 5.0
    expected2(1,7) = 0.0_f64
    expected2(2,7) = -0.6_f64
    expected2(3,7) = 0.42857142857142866_f64
    expected2(4,7) = 0.17142857142857143_f64   

    do j=1,num_tests
       x= x_test(j)

       ! We find its cell:
       cell = find_cell( spline, x )
       ! initialization accumulator:
       acc = 0.0_f64
       ! computing all non zero splines at point x:
       answer1(:) = b_splines_at_x(spline, cell, x)
       ! check difference with expected values
       do i = 1, degree + 1
          acc = max(acc,abs(answer1(i)-expected1(i,j)))
       end do
       passed_test = passed_test .and. (acc < criterion)
       ! corresponding prints:
       if( passed_test .eqv. .false. ) then
          print *, 'nonuniform splines test failure, spline values case:'
          print *, 'cell = ', cell, ',   x = ', x
          do i = 1, degree + 1
             print*, 'computed: ', answer1(i), 'expected = ', expected1(i,j)
          end do
          print*, 'Exiting...'
          stop
       else
          print *, 'test:', j, 'Evaluation of splines',  &
              ' at point ', x_test(j), "... PASSED"
       end if

       ! test spline derivatives
       ! computing derivatives of all non zero splines:
       answer2(:) = b_spline_derivatives_at_x(spline, cell, x)
       ! check difference with expected values
       acc = 0.0_f64
       do i = 1, degree + 1
          acc = max(acc,abs(answer2(i)-expected2(i,j)))
       end do
       passed_test = passed_test .and. (acc < criterion)
       if( passed_test .eqv. .false. ) then
          print *, 'periodic nonuniform splines test failure:'
          print *, 'computation of derivatives:'
          print *, 'cell = ', cell, 'x = ', x
          do i = 1, degree + 1
             print*, 'computed: ', answer2(i), 'expected = ', expected2(i,j)
          end do
          print *, 'error = ', acc
          print*, 'Exiting...'
          stop
       else
          print *, '                  Evaluation of splines derivatives', &
              ' at point ', x_test(j), "... PASSED"
       end if
    end do

    SLL_DEALLOCATE_ARRAY(answer1, ierr)
    SLL_DEALLOCATE_ARRAY(answer2, ierr)
    SLL_DEALLOCATE_ARRAY(answer3, ierr)

    call sll_delete(spline)

  end subroutine test_nonuniform_arb_deg_splines_periodic

  subroutine test_uniform_b_splines_randomly( passed_flag )
    logical, intent(inout)      :: passed_flag
    ! local variables
    sll_real64                  :: criterion
    sll_real64                  :: argument
    sll_real64                  :: argument_copy
    sll_real64, dimension(:), allocatable   :: results
    sll_real64, dimension(:), allocatable   :: derivatives
    sll_real64, dimension(:,:), allocatable :: sp_and_derivs
    sll_real64, dimension(:), allocatable   :: results_n
    sll_real64, dimension(:), allocatable   :: derivatives_n
    sll_real64, dimension(:,:), allocatable :: sp_and_derivs_n
    sll_int32                   :: num_tests
    sll_int32                   :: i
    sll_int32                   :: j
    sll_int32                   :: max_degree
    sll_int32                   :: ierr
    sll_int32, parameter        :: num_pts = 12
    sll_real64, dimension(num_pts)  :: grid 
    type(arbitrary_degree_spline_1d), pointer :: spline

    criterion          = 1.0d-14
    argument           = 0.0_f64
    num_tests          = 100
    argument_copy      = argument
    max_degree         = 10

    ! Define uniform grid for nonuniform spline
    grid = (/0.0_f64, 1.0_f64, 2.0_f64, 3.0_f64, 4.0_f64, 5.0_f64, 6.0_f64, &
         7.0_f64, 8.0_f64, 9.0_f64, 10.0_f64, 11.0_f64 /)

    do j=1, max_degree
       print*, '--------- degree=', j 
       SLL_CLEAR_ALLOCATE(results(1:j+1),ierr)
       SLL_CLEAR_ALLOCATE(derivatives(1:j+1),ierr)
       SLL_CLEAR_ALLOCATE(sp_and_derivs(1:2,1:j+1),ierr)
       SLL_CLEAR_ALLOCATE(results_n(1:j+1),ierr)
       SLL_CLEAR_ALLOCATE(derivatives_n(1:j+1),ierr)
       SLL_CLEAR_ALLOCATE(sp_and_derivs_n(1:2,1:j+1),ierr)

       ! creating non uniform 1d spline on uniform grid for comparison
       spline => new_arbitrary_degree_spline_1d( &
            j, &
            grid, &
            num_pts, &
            PERIODIC_ARBITRARY_DEG_SPLINE )

       do i=1,num_tests
          ! draw random number between 0 and 1
          call random_number(argument)
          ! compute values with uniform splines
          results(:) = uniform_b_splines_at_x(j, argument)
          derivatives(:) = uniform_b_spline_derivatives_at_x(j, argument)
          sp_and_derivs(:,:) = uniform_b_splines_and_derivs_at_x(j, argument)
          ! compute values with non uniform splines (cell is always 1)
          results_n(:) = b_splines_at_x(spline, 1, argument)
          derivatives_n(:) = b_spline_derivatives_at_x(spline, 1, argument)
          sp_and_derivs_n(:,:) = b_splines_and_derivs_at_x(spline, 1, argument)

          passed_flag = passed_flag .and. &
               (maxval(abs(results-results_n)).le.criterion)

          if( passed_flag .eqv. .false. ) then
             print *, 'Test_uniform_b_splines, values: wrong result for x = ', &
                  argument
             print *, 'Degree = ', j
             print*, 'uniform splines', results
             print*, 'nonuniform splines', results_n
             print*, 'Exiting...'
             stop
          end if

          passed_flag = passed_flag .and. &
               (maxval(abs(derivatives-derivatives_n)).le.criterion)
          if( passed_flag .eqv. .false. ) then
             print *, 'Test_uniform_b_splines, derivativess: ', &
                  'wrong result for x = ', &
                  argument
             print *, 'Degree = ', j
             print*, 'uniform splines', derivatives
             print*, 'nonuniform splines', derivatives_n
             print*, 'Exiting...'
             stop
          end if

          passed_flag = passed_flag .and. &
               (maxval(abs(sp_and_derivs(1,:)-sp_and_derivs_n(1,:) )) &
               .le.criterion) .and. &
               (maxval(abs(sp_and_derivs(2,:)-sp_and_derivs_n(2,:) )) &
               .le.criterion)
          if( passed_flag .eqv. .false. ) then
             print *, 'Test_uniform_b_splines, vals and derivs: ', &
                  'wrong result for x = ', argument
             print*, 'uniform splines values',  sp_and_derivs(1,:)
             print*, 'nonuniform splines values', sp_and_derivs_n(1,:)
             print*, 'uniform splines derivatives',  sp_and_derivs(2,:)
             print*, 'nonuniform spline derivs', sp_and_derivs_n(2,:)
             print*, 'Exiting...'
             stop
          end if
          results(:)         = 0.0_f64
          derivatives(:)     = 0.0_f64
          sp_and_derivs(:,:) = 0.0_f64

       end do
       SLL_DEALLOCATE_ARRAY(results, ierr)
       SLL_DEALLOCATE_ARRAY(derivatives, ierr)
       SLL_DEALLOCATE_ARRAY(sp_and_derivs, ierr)
       SLL_DEALLOCATE_ARRAY(results_n, ierr)
       SLL_DEALLOCATE_ARRAY(derivatives_n, ierr)
       SLL_DEALLOCATE_ARRAY(sp_and_derivs_n, ierr)
       call sll_delete(spline)
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
    criterion = 1.0d-15
    degree  = 3
    min_val = 0.0_f64
    num_pts = 10
    step    = 1.0_f64
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
    !print *, 'knots array = ', knots(:)

    ! fill spline object
    spline => new_arbitrary_degree_spline_1d( &
         degree, &
         knots, &
         num_pts, &
         OPEN_ARBITRARY_DEG_SPLINE )

    do j=1,num_tests
       call random_number(rnd)
       x = min_val + rnd*(knots(num_pts)-min_val)
       cell = find_cell( spline, x )
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
    call sll_delete(spline)
  end subroutine test_nonuniform_arb_deg_splines_open

  subroutine test_cpu_time
    ! local variables
    sll_real64, dimension(:), allocatable :: grid
    sll_int32  :: i,j
    sll_int32  :: num_pts
    sll_int32  :: degree
    sll_real64 :: min_val
    sll_int32  :: num_tests
    sll_real64 :: rnd
    sll_int32  :: ierr
    sll_int32, dimension(:), allocatable  :: cells
    sll_real64, dimension(:), allocatable :: x
    sll_real64, dimension(:), allocatable :: xx
    sll_real64, dimension(:), allocatable :: answer1
    sll_real64, dimension(:), allocatable :: answer2
    sll_real64, dimension(:,:), allocatable :: answer3
    type(arbitrary_degree_spline_1d), pointer :: spline
    type(sll_time_mark)  :: t0 
    sll_real64 :: time 

    ! Test performance of nonuniform arbitrary degree spline evaluation
    num_pts   = 10
    min_val = 0.0_f64
    degree    = 7
    num_tests = 100000
    print *, "Test performance of spline evaluation for"
    print *, " Spline degree = ", degree
    print *, " -----------------------------------------------------" 
    SLL_ALLOCATE(grid(num_pts),ierr)
    SLL_ALLOCATE(x(num_tests),ierr)
    SLL_ALLOCATE(xx(num_tests),ierr)
    SLL_ALLOCATE(cells(num_tests),ierr)
    SLL_ALLOCATE(answer1(degree+1),ierr)
    SLL_ALLOCATE(answer2(degree+1),ierr)
    SLL_ALLOCATE(answer3(2,degree+1),ierr)


    ! --------- 1D SPLINE INITIALIZATION ON NON UNIFORM MESH ----
    ! Creating non uniform mesh....
    grid(1) = min_val
    do i=2,num_pts
       call random_number(rnd)
       grid(i) = grid(i-1) + rnd !step
    end do
    ! ..... non uniform mesh done
    ! creating non uniform 1d spline 
    spline => new_arbitrary_degree_spline_1d( &
         degree, &
         grid, &
         num_pts, &
         PERIODIC_ARBITRARY_DEG_SPLINE )
    ! --------- INITIALIZATION DONE ------------

    do j=1,num_tests
       ! We compute a point randomly on the mesh:
       call random_number(rnd)
       x(j) = min_val + rnd*(grid(num_pts)-min_val)
       xx(j) = rnd
       ! We find its cell:
       cells(j) = find_cell( spline, x(j) )
    end do

    ! computing all non zero splines at all points in x:
    call sll_set_time_mark(t0)
    do j=1,num_tests
       answer1(:) = b_splines_at_x(spline, cells(j), x(j))
    end do
    time = sll_time_elapsed_since(t0)
    print *, 'Computing time for  b_splines_at_x: ', time

    call sll_set_time_mark(t0)
    do j=1,num_tests
       call compute_b_spline_at_x_mm(spline%knots, cells(j), x(j), degree, &
            answer1)
    end do
    time = sll_time_elapsed_since(t0)
    print *, 'Computing time for  compute_b_splines_at_x_mm: ', time
    ! computing all non zero spline derivatives at point x:
    call sll_set_time_mark(t0)
    do j=1,num_tests
       answer2(:) = b_spline_derivatives_at_x(spline, cells(j), x(j))
    end do
    time = sll_time_elapsed_since(t0)
    print *, 'Computing time for  b_spline_derivatives_at_x: ', time

    call sll_set_time_mark(t0)
    do j=1,num_tests
       call compute_b_spline_and_deriv_at_x_mm(spline%knots, cells(j), x(j), &
            degree, answer3)
    end do
    time = sll_time_elapsed_since(t0)
    print *, 'Computing time for  compute_b_splines_at_x_mm: ', time

    ! computing both all non zero splines and derivatives at point x:
    call sll_set_time_mark(t0)
    do j=1,num_tests
       answer3(:,:) = b_splines_and_derivs_at_x(spline, cells(j), x(j))
    end do
    time = sll_time_elapsed_since(t0)
    print *, 'Computing time for  b_splines_and_derivs_at_x: ', time
    ! computing all non zero uniform splines  at point x:
    call sll_set_time_mark(t0)
    do j=1,num_tests
       answer1(:) = uniform_b_splines_at_x(degree, xx(j))
    end do
    time = sll_time_elapsed_since(t0)
    print *, 'Computing time for  uniform_b_splines_at_x: ', time
    ! computing all non zero uniform splines derivatives at point x:
    call sll_set_time_mark(t0)
    do j=1,num_tests
       answer2(:) = uniform_b_spline_derivatives_at_x(degree, xx(j))
    end do
    time = sll_time_elapsed_since(t0)
    print *, 'Computing time for  uniform_b_spline_derivatives_at_x: ', time
    ! computing all non zero uniform splines and derivatives at point x:
    call sll_set_time_mark(t0)
    do j=1,num_tests
       answer3(:,:) = uniform_b_splines_and_derivs_at_x(degree, xx(j))
    end do
    time = sll_time_elapsed_since(t0)
    print *, 'Computing time for  uniform_b_splines_and_derivs_at_x: ', time

  end subroutine test_cpu_time

end program test_arbitrary_degree_splines
