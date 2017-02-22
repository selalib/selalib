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

! PROGRAM: test_characteristics_1d_trapezoid
!
! DESCRIPTION:
!> @ingroup characteristics
!> @brief   unit test for 'characteristics_1D_trapezoid' integrator
!> @author  Yaman Güçlü

program test_characteristics_1d_trapezoid
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"
#include "sll_errors.h"

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_periodic

  use sll_m_characteristics_1d_base, only: &
    sll_c_characteristics_1d_base

  use sll_m_characteristics_1d_trapezoid, only: &
    sll_f_new_trapezoid_1d_charac

  use sll_m_cubic_spline_interpolator_1d, only: &
    sll_f_new_cubic_spline_interpolator_1d

  use sll_m_interpolators_1d_base, only: &
    sll_c_interpolator_1d

  use m_characteristics_1d_test_helper, only: &
    sll_c_ch1d_test_case, &
    sll_t_ch1d_constant_flow,  &
    sll_t_ch1d_stationary_compression_wave, &
    sll_s_test_characteristics_1d

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  class(sll_c_interpolator_1d),         pointer :: A_interp => null()
  class(sll_c_characteristics_1d_base), pointer :: trap     => null()
  class(sll_c_ch1d_test_case),      allocatable :: test_case

  sll_int32  :: num_cells
  sll_real64 :: dt
  sll_real64 :: tol
  sll_real64 :: error
  logical    :: passed
  sll_int32  :: i

  !-----------------------------------------------------------------------------
  ! Discretization parameter
  num_cells = 16

  !-----------------------------------------------------------------------------
  passed = .true.

  ! Constant advection: numerical results should be "exact"
  dt  = 0.7_f64
  tol = 1e-14_f64
  allocate( sll_t_ch1d_constant_flow :: test_case )
  call s_init_characteristics_1d_trapezoid( trap, A_interp, test_case, num_cells )
  error = sll_s_test_characteristics_1d( trap, test_case, num_cells, dt )
  print *, "Constant advection test-case: Tolerance =", tol
  print *, "max(abs(error)) = ", error
  if( error > tol ) passed = .false.
  deallocate( trap      )
  deallocate( A_interp  )
  deallocate( test_case )

  ! Stationary compression wave: numerical result has error O(dt^3)
  dt  = 0.2_f64
  tol = 2e-3_f64
  allocate( sll_t_ch1d_stationary_compression_wave :: test_case )
  call s_init_characteristics_1d_trapezoid( trap, A_interp, test_case, num_cells )
  do i=1,5
    error = sll_s_test_characteristics_1d( trap, test_case, num_cells, dt )
    print *, "Stationary compression wave test-case: Tolerance =", tol
    print *, "max(abs(error)) = ", error
    if( error > tol ) passed = .false.
    dt  = dt  /  2.0_f64
    tol = tol / (2.0_f64**3)
  end do
  deallocate( trap      )
  deallocate( A_interp  )
  deallocate( test_case )

  ! Unit test output
  if( passed ) then
    print *, '#PASSED'
  endif

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
contains
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine s_init_characteristics_1d_trapezoid( &
      trap, A_interp, test_case, num_cells )
    class(sll_c_characteristics_1d_base), pointer    :: trap
    class(sll_c_interpolator_1d),         pointer    :: A_interp
    class(sll_c_ch1d_test_case),          intent(in) :: test_case
    sll_int32,                            intent(in) :: num_cells

    character(len=*), parameter :: this_sub_name = &
      "s_init_characteristic_1d_trapezoid"

    if( associated( trap ) ) then
      SLL_ERROR( this_sub_name, "pointer argument 'trap' already associated" )
    else if( associated( A_interp ) ) then
      SLL_ERROR( this_sub_name, "pointer argument 'A_interp' already associated" )
    end if

    ! Periodic interpolator for flow field
    A_interp => sll_f_new_cubic_spline_interpolator_1d( &
      num_points = num_cells+1,         &
      xmin       = test_case%eta_min(), &
      xmax       = test_case%eta_max(), &
      bc_type    = sll_p_periodic )

    ! OBJECT TO BE TESTED: trapezoidal integrator of characteristic trajectories
    trap => sll_f_new_trapezoid_1d_charac( &
      Npts     = num_cells+1,         &
      A_interp = A_interp,            &
      bc_type  = sll_p_periodic,      &
      eta_min  = test_case%eta_min(), &
      eta_max  = test_case%eta_max() )

  end subroutine s_init_characteristics_1d_trapezoid


end program test_characteristics_1d_trapezoid
