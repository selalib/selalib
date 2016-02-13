!> @ingroup sll_t_operator_splitting
!> @brief Unit test for operator splitting. Order check for linear pendulum 
!> 

program test_linear_pendulum
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_linear_pendulum_operators, only: &
    sll_s_check_order

  use sll_m_operator_splitting, only: &
    sll_p_lie_tv, &
    sll_p_lie_vt, &
    sll_p_order6_tvt, &
    sll_p_order6_vtv, &
    sll_p_strang_tvt, &
    sll_p_strang_vtv, &
    sll_p_triple_jump_tvt, &
    sll_p_triple_jump_vtv

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
  ! variables
  logical :: test_passed
  logical :: all_tests_passed = .true.
  sll_real64 :: steps_fine
  sll_int32  :: expected_order
  
  ! test sll_p_lie_tv
  steps_fine = real(200,f64)
  expected_order = 1
  test_passed = .true.
  call sll_s_check_order(sll_p_lie_tv,steps_fine, expected_order, test_passed)
  if (.not.(test_passed)) then
     print*, 'Problem with order of sll_p_lie_tv'
     all_tests_passed = .false.
  end if

  ! test sll_p_lie_vt
  steps_fine = real(200,f64)
  expected_order = 1
  test_passed = .true.
  call sll_s_check_order(sll_p_lie_vt,steps_fine, expected_order, test_passed)
   if (.not.(test_passed)) then
     print*, 'Problem with order of sll_p_lie_vt'
     all_tests_passed = .false.
  end if

  ! test sll_p_strang_tvt
  steps_fine = real(100,f64)
  expected_order = 2
  test_passed = .true.
  call sll_s_check_order(sll_p_strang_tvt,steps_fine, expected_order, test_passed)
  if (.not.(test_passed)) then
     print*, 'Problem with order of sll_p_strang_tvt'
     all_tests_passed = .false.
  end if

  ! test sll_p_strang_vtv
  steps_fine = real(100,f64)
  expected_order = 2
  test_passed = .true.
  call sll_s_check_order(sll_p_strang_vtv,steps_fine, expected_order, test_passed)
  if (.not.(test_passed)) then
     print*, 'Problem with order of sll_p_strang_vtv'
     all_tests_passed = .false.
  end if

  ! test sll_p_triple_jump_tvt
  steps_fine = real(64,f64)
  expected_order = 4
  test_passed = .true.
  call sll_s_check_order(sll_p_triple_jump_tvt,steps_fine, expected_order, test_passed)
  if (.not.(test_passed)) then
     print*, 'Problem with order of sll_p_triple_jump_tvt'
     all_tests_passed = .false.
  end if
  
  ! test sll_p_triple_jump_vtv
  steps_fine = real(64,f64)
  expected_order = 4
  test_passed = .true.
  call sll_s_check_order(sll_p_triple_jump_vtv,steps_fine, expected_order, test_passed)
  if (.not.(test_passed)) then
     print*, 'Problem with order of sll_p_triple_jump_vtv'
     all_tests_passed = .false.
  end if
  
  ! test sll_p_order6_tvt 
  steps_fine = real(20,f64)
  expected_order = 6
  test_passed = .true.
  call sll_s_check_order(sll_p_order6_tvt,steps_fine, expected_order, test_passed)
  if (.not.(test_passed)) then
     print*, 'Problem with order of sll_p_order6_tvt'
     all_tests_passed = .false.
  end if

  ! test sll_p_order6_vtv 
  steps_fine = real(20,f64)
  expected_order = 6
  test_passed = .true.
  call sll_s_check_order(sll_p_order6_vtv,steps_fine, expected_order, test_passed)
  if (.not.(test_passed)) then
     print*, 'Problem with order of sll_p_order6_vtv'
     all_tests_passed = .false.
  end if

  if (test_passed) print*, 'PASSED'

end program test_linear_pendulum

