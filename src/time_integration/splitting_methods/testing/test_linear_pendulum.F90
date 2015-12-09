!> @ingroup operator_splitting
!> @brief Unit test for operator splitting. Order check for linear pendulum 
!> 

program test_linear_pendulum
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_linear_pendulum_operators, only: &
    check_order

  use sll_m_operator_splitting, only: &
    sll_lie_tv, &
    sll_lie_vt, &
    sll_order6_tvt, &
    sll_order6_vtv, &
    sll_strang_tvt, &
    sll_strang_vtv, &
    sll_triple_jump_tvt, &
    sll_triple_jump_vtv

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
  ! variables
  logical :: test_passed
  logical :: all_tests_passed = .true.
  sll_real64 :: steps_fine
  sll_int32  :: expected_order
  
  ! test SLL_LIE_TV
  steps_fine = real(200,f64)
  expected_order = 1
  test_passed = .true.
  call check_order(SLL_LIE_TV,steps_fine, expected_order, test_passed)
  if (.not.(test_passed)) then
     print*, 'Problem with order of SLL_LIE_TV'
     all_tests_passed = .false.
  end if

  ! test SLL_LIE_VT
  steps_fine = real(200,f64)
  expected_order = 1
  test_passed = .true.
  call check_order(SLL_LIE_VT,steps_fine, expected_order, test_passed)
   if (.not.(test_passed)) then
     print*, 'Problem with order of SLL_LIE_VT'
     all_tests_passed = .false.
  end if

  ! test SLL_STRANG_TVT
  steps_fine = real(100,f64)
  expected_order = 2
  test_passed = .true.
  call check_order(SLL_STRANG_TVT,steps_fine, expected_order, test_passed)
  if (.not.(test_passed)) then
     print*, 'Problem with order of SLL_STRANG_TVT'
     all_tests_passed = .false.
  end if

  ! test SLL_STRANG_VTV
  steps_fine = real(100,f64)
  expected_order = 2
  test_passed = .true.
  call check_order(SLL_STRANG_VTV,steps_fine, expected_order, test_passed)
  if (.not.(test_passed)) then
     print*, 'Problem with order of SLL_STRANG_VTV'
     all_tests_passed = .false.
  end if

  ! test SLL_TRIPLE_JUMP_TVT
  steps_fine = real(64,f64)
  expected_order = 4
  test_passed = .true.
  call check_order(SLL_TRIPLE_JUMP_TVT,steps_fine, expected_order, test_passed)
  if (.not.(test_passed)) then
     print*, 'Problem with order of SLL_TRIPLE_JUMP_TVT'
     all_tests_passed = .false.
  end if
  
  ! test SLL_TRIPLE_JUMP_VTV
  steps_fine = real(64,f64)
  expected_order = 4
  test_passed = .true.
  call check_order(SLL_TRIPLE_JUMP_VTV,steps_fine, expected_order, test_passed)
  if (.not.(test_passed)) then
     print*, 'Problem with order of SLL_TRIPLE_JUMP_VTV'
     all_tests_passed = .false.
  end if
  
  ! test SLL_ORDER6_TVT 
  steps_fine = real(20,f64)
  expected_order = 6
  test_passed = .true.
  call check_order(SLL_ORDER6_TVT,steps_fine, expected_order, test_passed)
  if (.not.(test_passed)) then
     print*, 'Problem with order of SLL_ORDER6_TVT'
     all_tests_passed = .false.
  end if

  ! test SLL_ORDER6_VTV 
  steps_fine = real(20,f64)
  expected_order = 6
  test_passed = .true.
  call check_order(SLL_ORDER6_VTV,steps_fine, expected_order, test_passed)
  if (.not.(test_passed)) then
     print*, 'Problem with order of SLL_ORDER6_VTV'
     all_tests_passed = .false.
  end if

  if (test_passed) print*, 'PASSED'

end program test_linear_pendulum

