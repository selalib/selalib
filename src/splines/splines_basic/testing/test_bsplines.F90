program test_arbitrary_degree_splines
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_working_precision, only: f64

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_periodic, &
    sll_p_open

  use sll_m_low_level_bsplines, only: &
    sll_t_bsplines                        , &
    sll_s_bsplines_init_from_grid         , &
    sll_s_bsplines_free                   , &
    sll_s_bsplines_eval_basis             , &
    sll_s_bsplines_eval_deriv             , &
    sll_s_bsplines_eval_basis_and_deriv   , &
    sll_s_bsplines_eval_basis_and_n_derivs, &
    sll_f_find_cell, &
    sll_s_bsplines_eval_basis_mm          , &
    sll_s_bsplines_eval_basis_and_deriv_mm, &
    sll_s_uniform_bsplines_eval_basis          , &
    sll_s_uniform_bsplines_eval_deriv          , &
    sll_s_uniform_bsplines_eval_basis_and_deriv, &
    sll_s_uniform_bsplines_eval_basis_and_n_derivs

  use sll_m_timer, only: &
    sll_s_set_time_mark, &
    sll_f_time_elapsed_since, &
    sll_t_time_mark

  ! new interface
  use sll_m_bsplines_base, only: &
    sll_c_bsplines

  use sll_m_bsplines_uniform, only: &
    sll_t_bsplines_uniform

  use sll_m_bsplines_non_uniform, only: &
    sll_t_bsplines_non_uniform

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Working precision
  integer, parameter :: wp = f64

  ! Output format for timing
  character(len=*), parameter :: timing_fmt = "('Computing time for  ',a,':',es12.2)"

  character(len=54) :: subr
  logical           :: passed_test

  passed_test = .true.

  write(*,*)
  write(*,'(a)') '*****************************************************************'
  write(*,'(a)') ' Compare uniform and non-uniform B-splines on uniform grid'
  write(*,'(a)') '*****************************************************************'
  write(*,*)
  call compare_uniform_and_non_uniform_bsplines( passed_test )
  write(*,*)
  write(*,'(a)') '*****************************************************************'
  write(*,'(a)') ' Non-uniform B-splines, periodic '
  write(*,'(a)') '*****************************************************************'
  write(*,*)
  call test_non_uniform_bsplines_periodic( passed_test )
  write(*,*)
  write(*,'(a)') '*****************************************************************'
  write(*,'(a)') ' Non-uniform B-splines, open '
  write(*,'(a)') '*****************************************************************'
  write(*,*)
  call test_non_uniform_bsplines_open( passed_test )
  write(*,*)
  write(*,'(a)') '*****************************************************************'
  write(*,'(a)') ' Test performance (CPU timing) '
  write(*,'(a)') '*****************************************************************'
  write(*,*)
  call test_performance

  if (passed_test) then
     write(*,*) 'PASSED'
  else
     write(*,*) 'FAILED'
  end if

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !----------------------------------------------------------------------------
  subroutine compare_uniform_and_non_uniform_bsplines( passed_test )
    logical, intent(inout) :: passed_test

    integer  :: i, p, jmin
    integer  :: ncells, num_tests
    integer  :: degree, max_degree
    real(wp) :: x
    real(wp) :: tolerance, tolerance_derivs

    real(wp), allocatable :: grid(:)
    real(wp), allocatable :: values_old(:)
    real(wp), allocatable :: values_new(:)
    real(wp), allocatable :: values_old_nu(:)
    real(wp), allocatable :: values_new_nu(:)
    real(wp), allocatable :: derivs_old(:)
    real(wp), allocatable :: derivs_new(:)
    real(wp), allocatable :: derivs_old_nu(:)
    real(wp), allocatable :: derivs_new_nu(:)
    real(wp), allocatable :: values_and_deriv(:,:)
    real(wp), allocatable :: values_and_deriv_nu(:,:)
    real(wp), allocatable :: values_and_n_derivs_old(:,:)
    real(wp), allocatable :: values_and_n_derivs_new(:,:)
    real(wp), allocatable :: values_and_n_derivs_old_nu(:,:)
    real(wp), allocatable :: values_and_n_derivs_new_nu(:,:)

    ! B-splines, old type
    type(sll_t_bsplines) :: bsplines
    ! Uniform and non-uniform B-splines, new interface
    type(sll_t_bsplines_uniform)     :: bsplines_uniform
    type(sll_t_bsplines_non_uniform) :: bsplines_non_uniform

    ncells     = 11
    num_tests  = 100
    max_degree = 9
    tolerance  = 1.0e-15_wp

    ! define uniform grid
    allocate( grid(ncells+1) )
    grid(:) = [(0.0_wp+real(i,wp), i=0,ncells)]

    do degree = 1, max_degree

      write(*,'(a,i0)',advance='no') 'degree = ', degree
      allocate(values_old   (degree+1)); values_old    = 0.0_wp
      allocate(values_new   (degree+1)); values_new    = 0.0_wp
      allocate(values_old_nu(degree+1)); values_old_nu = 0.0_wp
      allocate(values_new_nu(degree+1)); values_new_nu = 0.0_wp

      allocate(derivs_old   (degree+1)); derivs_old    = 0.0_wp
      allocate(derivs_new   (degree+1)); derivs_new    = 0.0_wp
      allocate(derivs_old_nu(degree+1)); derivs_old_nu = 0.0_wp
      allocate(derivs_new_nu(degree+1)); derivs_new_nu = 0.0_wp

      allocate(values_and_deriv   (0:1,degree+1)); values_and_deriv    = 0.0_wp
      allocate(values_and_deriv_nu(0:1,degree+1)); values_and_deriv_nu = 0.0_wp

      allocate(values_and_n_derivs_old   (0:degree,degree+1)); values_and_n_derivs_old    = 0.0_wp
      allocate(values_and_n_derivs_new   (0:degree,degree+1)); values_and_n_derivs_new    = 0.0_wp
      allocate(values_and_n_derivs_old_nu(0:degree,degree+1)); values_and_n_derivs_old_nu = 0.0_wp
      allocate(values_and_n_derivs_new_nu(0:degree,degree+1)); values_and_n_derivs_new_nu = 0.0_wp

      ! initialize B-splines, old type
      call sll_s_bsplines_init_from_grid( bsplines, degree, grid, sll_p_periodic, sll_p_periodic )
      ! initialize uniform B-splines, new interface
      call bsplines_uniform % init( degree, periodic=.true., xmin=grid(1), xmax=grid(ncells+1), ncells=ncells )
      ! initialize non-uniform B-splines, new interface
      call bsplines_non_uniform % init( degree, periodic=.true., breaks=grid )

      do i=1,num_tests

        ! draw random number between 0 and 1
        call random_number(x)

        ! compute with uniform B-splines
        call sll_s_uniform_bsplines_eval_basis( degree, x, values_old )
        call sll_s_uniform_bsplines_eval_deriv( degree, x, derivs_old )
        call sll_s_uniform_bsplines_eval_basis_and_deriv( degree, x, values_and_deriv )
        call sll_s_uniform_bsplines_eval_basis_and_n_derivs( degree, x, degree, values_and_n_derivs_old )

        ! compute with uniform B-splines, new interface
        call bsplines_uniform % eval_basis( x, values_new, jmin )
        call bsplines_uniform % eval_deriv( x, derivs_new, jmin )
        call bsplines_uniform % eval_basis_and_n_derivs( x, degree, values_and_n_derivs_new, jmin )

        ! compute with non uniform B-splines (cell is always 1)
        call sll_s_bsplines_eval_basis( bsplines, 1, x, values_old_nu )
        call sll_s_bsplines_eval_deriv( bsplines, 1, x, derivs_old_nu )
        call sll_s_bsplines_eval_basis_and_deriv( bsplines, 1, x, values_and_deriv_nu )
        call sll_s_bsplines_eval_basis_and_n_derivs( bsplines, 1, x, degree, values_and_n_derivs_old_nu )

        ! compute with non uniform B-splines, new interface
        call bsplines_non_uniform % eval_basis( x, values_new_nu, jmin )
        call bsplines_non_uniform % eval_deriv( x, derivs_new_nu, jmin )
        call bsplines_non_uniform % eval_basis_and_n_derivs( x, degree, values_and_n_derivs_new_nu, jmin )

        passed_test = passed_test .and. &
          maxval( abs( values_old-values_new    ) ) < tolerance .and. &
          maxval( abs( values_old-values_old_nu ) ) < tolerance .and. &
          maxval( abs( values_old-values_new_nu ) ) < tolerance

        !TODO: improve output
        if( .not. passed_test ) then
          write(*,*) 'compare_uniform_and_non_uniform_bsplines, values: wrong result for x = ', x
          write(*,*) 'degree = ', degree
          write(*,*) 'uniform B-splines,     old type:', values_old
          write(*,*) 'uniform B-splines,     new type:', values_new
          write(*,*) 'non-uniform B-splines, old type:', values_old_nu
          write(*,*) 'non-uniform B-splines, new type:', values_new_nu
          write(*,*) 'exiting...'
          stop
        end if

        passed_test = passed_test .and. &
          maxval( abs( derivs_old-derivs_new    ) ) < tolerance .and. &
          maxval( abs( derivs_old-derivs_old_nu ) ) < tolerance .and. &
          maxval( abs( derivs_old-derivs_new_nu ) ) < tolerance

        if( .not. passed_test ) then
          write(*,*) 'compare_uniform_and_non_uniform_bsplines, derivs: ', 'wrong result for x = ', x
          write(*,*) 'degree = ', degree
          write(*,*) 'uniform B-splines,     old type:', derivs_old
          write(*,*) 'uniform B-splines,     new type:', derivs_new
          write(*,*) 'non-uniform B-splines, old type:', derivs_old_nu
          write(*,*) 'non-uniform B-splines, new type:', derivs_new_nu
          write(*,*) 'exiting...'
          stop
        end if

        passed_test = passed_test .and. &
          maxval( abs( values_and_deriv(0,:)-values_and_deriv_nu(0,:) ) ) < tolerance .and. &
          maxval( abs( values_and_deriv(1,:)-values_and_deriv_nu(1,:) ) ) < tolerance

        if( .not. passed_test ) then
          write(*,*) 'compare_uniform_and_non_uniform_bsplines, values and derivs: ', 'wrong result for x = ', x
          write(*,*) 'uniform B-splines,     values:', values_and_deriv   (0,:)
          write(*,*) 'non-uniform B-splines, values:', values_and_deriv_nu(0,:)
          write(*,*) 'uniform B-splines,     derivs:', values_and_deriv   (1,:)
          write(*,*) 'non-uniform B-splines, derivs:', values_and_deriv_nu(1,:)
          write(*,*) 'exiting...'
          stop
        end if

        do p = 0, degree
          tolerance_derivs = tolerance * 10.0_wp**p
          passed_test = passed_test .and. &
            maxval( abs( values_and_n_derivs_old(p,:)-values_and_n_derivs_new(p,:)    ) ) < tolerance_derivs &
            .and. &
            maxval( abs( values_and_n_derivs_old(p,:)-values_and_n_derivs_old_nu(p,:) ) ) < tolerance_derivs &
            .and. &
            maxval( abs( values_and_n_derivs_old(p,:)-values_and_n_derivs_new_nu(p,:) ) ) < tolerance_derivs

          if( .not. passed_test ) then
             write(*,*) 'compare_uniform_and_non_uniform_bsplines, values and n derivs: ', 'wrong result for x = ', x
             write(*,*) 'uniform B-splines,     old type:', p, values_and_n_derivs_old(p,:)
             write(*,*) 'uniform B-splines,     new type:', p, values_and_n_derivs_new(p,:)
             write(*,*) 'non-uniform B-splines, old type:', p, values_and_n_derivs_old_nu(p,:)
             write(*,*) 'non-uniform B-splines, new type:', p, values_and_n_derivs_new_nu(p,:)
             write(*,*) 'exiting...'
             stop
          end if
        end do

      end do

      if( passed_test ) write(*,*) '  OK'

      deallocate(values_old)
      deallocate(values_new)
      deallocate(values_old_nu)
      deallocate(values_new_nu)
      deallocate(derivs_old)
      deallocate(derivs_new)
      deallocate(derivs_old_nu)
      deallocate(derivs_new_nu)
      deallocate(values_and_deriv   )
      deallocate(values_and_deriv_nu)
      deallocate(values_and_n_derivs_old)
      deallocate(values_and_n_derivs_new)
      deallocate(values_and_n_derivs_old_nu)
      deallocate(values_and_n_derivs_new_nu)

      ! Free B-splines, old type
      call sll_s_bsplines_free( bsplines )
      ! Free uniform B-splines, new interface
      call bsplines_uniform % free()
      ! Free non-uniform B-splines, new interface
      call bsplines_non_uniform % free()

    end do

  end subroutine compare_uniform_and_non_uniform_bsplines

  !----------------------------------------------------------------------------
  subroutine test_non_uniform_bsplines_periodic( passed_test )
    logical, intent(inout) :: passed_test

    integer  :: i, j, jmin, cell
    integer  :: num_pts, num_tests
    integer  :: degree
    real(wp) :: x, xmin, rnd
    real(wp) :: tolerance
    real(wp), allocatable :: grid(:)
    real(wp), allocatable :: x_test(:)
    real(wp), allocatable :: expected1(:,:)
    real(wp), allocatable :: expected2(:,:)
    real(wp), allocatable :: values_old(:)
    real(wp), allocatable :: values_new(:)
    real(wp), allocatable :: derivs_old(:)
    real(wp), allocatable :: derivs_new(:)
    real(wp), allocatable :: values_and_deriv(:,:)
    real(wp), allocatable :: values_and_n_derivs_old(:,:)
    real(wp), allocatable :: values_and_n_derivs_new(:,:)

    ! B-splines, old type
    type(sll_t_bsplines) :: bsplines
    ! Non-uniform B-splines, new interface
    type(sll_t_bsplines_non_uniform) :: bsplines_non_uniform

    !--------------------------------------------------------------------------
    ! Test on random grid for given degree
    !--------------------------------------------------------------------------

    num_tests = 10
    num_pts   = 10
    xmin      = 0.0_wp
    tolerance = 1.0e-15_wp

    degree = 6
    write(*,'(a)') "Test on random grid for given degree"
    write(*,'(a)') "------------------------------------"
    write(*,*)
    write(*,'(a,i0,a)') "degree = ", degree, ":"
    write(*,*)

    allocate(grid(num_pts))
    allocate(values_old(degree+1))
    allocate(values_new(degree+1))
    allocate(derivs_old(degree+1))
    allocate(derivs_new(degree+1))
    allocate(values_and_deriv(2,degree+1))
    allocate(values_and_n_derivs_old(2,degree+1))
    allocate(values_and_n_derivs_new(2,degree+1))

    ! create non-uniform mesh
    grid(1) = xmin
    do i = 2, num_pts
      call random_number(rnd)
      grid(i) = grid(i-1) + rnd !step
    end do

    ! initialize B-splines, old type
    call sll_s_bsplines_init_from_grid( bsplines, degree, grid, sll_p_periodic, sll_p_periodic )
    ! initialize non-uniform B-splines, new interface
    call bsplines_non_uniform % init( degree, periodic=.true., breaks=grid )

    do j = 1, num_tests

      ! compute a point randomly on the mesh
      call random_number(rnd)
      x = xmin + rnd*(grid(num_pts)-xmin)

      cell = sll_f_find_cell( bsplines, x )

      ! compute with non-uniform B-splines, old type
      call sll_s_bsplines_eval_basis( bsplines, cell, x, values_old )
      call sll_s_bsplines_eval_deriv( bsplines, cell, x, derivs_old )
      call sll_s_bsplines_eval_basis_and_deriv( bsplines, cell, x, values_and_deriv )
      call sll_s_bsplines_eval_basis_and_n_derivs( bsplines, cell, x, 1, values_and_n_derivs_old )

      ! compute with non uniform B-splines, new interface
      call bsplines_non_uniform % eval_basis( x, values_new, jmin )
      call bsplines_non_uniform % eval_deriv( x, derivs_new, jmin )
      call bsplines_non_uniform % eval_basis_and_n_derivs( x, 1, values_and_n_derivs_new, jmin )

      ! test partition of unity
      passed_test = passed_test .and. &
        abs( 1.0_wp - sum( values_old(1:degree+1) ) ) < tolerance .and. &
        abs( 1.0_wp - sum( values_new(1:degree+1) ) ) < tolerance

      if( .not. passed_test ) then
        write(*,*) 'test_non_uniform_bsplines_periodic, values: wrong result for x = ', x
        write(*,*) 'non-uniform B-splines, old type, sum of values: ', sum( values_old(1:degree+1) )
        write(*,*) 'non-uniform B-splines, new type, sum of values: ', sum( values_new(1:degree+1) )
        write(*,*) 'exiting...'
        stop
      else
        write(*,'(a,i2,a,i2,a)') "test ", j, ", cell ", cell, ": testing partition of unity   OK"
      end if

      ! test sum of all derivatives (should be 0)
      passed_test = passed_test .and. &
        abs( sum( derivs_old(1:degree+1) ) ) < tolerance .and. &
        abs( sum( derivs_new(1:degree+1) ) ) < tolerance

      if( .not. passed_test ) then
        write(*,*) 'test_non_uniform_bsplines_periodic, derivs: wrong result for x = ', x
        write(*,*) 'non-uniform B-splines, old type, sum of derivs: ', sum( derivs_old(1:degree+1) )
        write(*,*) 'non-uniform B-splines, new type, sum of derivs: ', sum( derivs_new(1:degree+1) )
        write(*,*) 'exiting...'
        stop
      else
        write(*,'(a)') "                  testing sum of derivatives   OK"
      end if

      ! test methods to obtain both values and derivs
      passed_test = passed_test .and. &
        maxval( abs( values_old(:)-values_and_deriv(1,:) ) ) < tolerance .and. &
        maxval( abs( values_new(:)-values_and_deriv(1,:) ) ) < tolerance

      if( .not. passed_test ) then
        write(*,*) 'test_non_uniform_bsplines_periodic, values and derivs: wrong result for x = ', x
        write(*,*) 'non-uniform B-splines, old type, values: ', values_old
        write(*,*) 'non-uniform B-splines, new type, values: ', values_new
        write(*,*) 'non-uniform B-splines, old type, values_and_deriv: ', values_and_deriv(1,:)
        write(*,*) 'exiting...'
        stop
      else
        write(*,'(a)') "                  consistency checks           OK"
      end if

      passed_test = passed_test .and. &
        maxval( abs( derivs_old(:)-values_and_deriv(2,:) ) ) < tolerance .and. &
        maxval( abs( derivs_new(:)-values_and_deriv(2,:) ) ) < tolerance

      if( .not. passed_test ) then
        write(*,*) 'test_non_uniform_bsplines_periodic, values and derivs: wrong result for x = ', x
        write(*,*) 'non-uniform B-splines, old type, values: ', derivs_old
        write(*,*) 'non-uniform B-splines, new type, values: ', derivs_new
        write(*,*) 'non-uniform B-splines, old type, values_and_deriv: ', values_and_deriv(2,:)
        write(*,*) 'exiting...'
        stop
      else
        write(*,'(a)') "                  consistency checks           OK"
      end if

      passed_test = passed_test .and. &
        maxval( abs( values_and_deriv(2,:)-values_and_n_derivs_old(2,:) ) ) < tolerance .and. &
        maxval( abs( values_and_deriv(2,:)-values_and_n_derivs_new(2,:) ) ) < tolerance

      if( .not. passed_test ) then
        write(*,*) 'test_non_uniform_bsplines_periodic, values and derivs: wrong result for x = ', x
        write(*,*) 'non-uniform B-splines, old type, values_and_deriv: ', values_and_deriv(2,:)
        write(*,*) 'non-uniform B-splines, old type, values_and_n_derivs: ', values_and_n_derivs_old(2,:)
        write(*,*) 'non-uniform B-splines, new type, values_and_n_derivs: ', values_and_n_derivs_new(2,:)
        write(*,*) 'exiting...'
        stop
      else
        write(*,'(a)') "                  consistency checks           OK"
      end if

    end do

    deallocate(grid)
    deallocate(values_old)
    deallocate(values_new)
    deallocate(derivs_old)
    deallocate(derivs_new)
    deallocate(values_and_deriv)
    deallocate(values_and_n_derivs_old)
    deallocate(values_and_n_derivs_new)

    call sll_s_bsplines_free(bsplines)

    !--------------------------------------------------------------------------
    ! Test on given grid with known answer
    !--------------------------------------------------------------------------

    num_tests = 7
    num_pts   = 5
    xmin      = 0.0_wp
    tolerance = 1.0e-15_wp

    degree = 3
    write(*,*)
    write(*,'(a)') "Test on given grid with known answer"
    write(*,'(a)') "------------------------------------"
    write(*,*)
    write(*,'(a,i0,a)') "degree = ", degree, ":"
    write(*,*)

    allocate(grid(num_pts))
    allocate(values_old(degree+1))
    allocate(values_new(degree+1))
    allocate(derivs_old(degree+1))
    allocate(derivs_new(degree+1))
    allocate(values_and_deriv(2,degree+1))

    allocate(x_test(num_tests))
    allocate(expected1(degree+1,num_tests))
    allocate(expected2(degree+1,num_tests))

    ! define non-uniform grid
    grid = (/0.0_wp, 2.0_wp, 3.0_wp, 4.5_wp, 5.0_wp /)

    ! initialize B-splines, old type
    call sll_s_bsplines_init_from_grid( bsplines, degree, grid, sll_p_periodic, sll_p_periodic )
    ! initialize non-uniform B-splines, new interface
    call bsplines_non_uniform % init( degree, periodic=.true., breaks=grid )

    ! array of points where splines will be evaluated
    x_test = (/0._wp, 0.5_wp, 1.2_wp, 2.001_wp, 3._wp, 4.0_wp, 5.0_wp /)

    ! Expected values of bsplines and derivatives at these points
    ! expected(:,i): values expected at point x_test(i)
    ! values:
    expected1(1,1) = 0.4_wp
    expected1(2,1) = 0.5714285714285714_wp
    expected1(3,1) = 0.028571428571428574_wp
    expected1(4,1) = 0.0_wp
    !
    expected1(1,2) = 0.16874999999999999_wp
    expected1(2,2) = 0.644345238095238_wp
    expected1(3,2) = 0.18227513227513226_wp
    expected1(4,2) = 0.004629629629629629_wp
    !
    expected1(1,3) = 0.0256_wp
    expected1(2,3) = 0.42742857142857144_wp
    expected1(3,3) = 0.4829714285714285_wp
    expected1(4,3) = 0.06399999999999999_wp
    !
    expected1(1,4) = 0.0949526665714286_wp
    expected1(2,4) = 0.6083063706285714_wp
    expected1(3,4) = 0.2967409626666666_wp
    expected1(4,4) = 1.3333333333328926e-10_wp
    !
    expected1(1,5) = 0.2_wp
    expected1(2,5) = 0.6666666666666666_wp
    expected1(3,5) = 0.13333333333333333_wp
    expected1(4,5) = 0.0_wp
    !
    expected1(1,6) = 0.007407407407407407_wp
    expected1(2,6) = 0.25925925925925924_wp
    expected1(3,6) = 0.65_wp
    expected1(4,6) = 0.08333333333333333_wp
    !
    expected1(1,7) = 0.0_wp
    expected1(2,7) = 0.4_wp
    expected1(3,7) = 0.5714285714285714_wp
    expected1(4,7) = 0.028571428571428574_wp
    ! derivatives:
    expected2(1,1) = -0.6_wp
    expected2(2,1) = 0.42857142857142866_wp
    expected2(3,1) = 0.17142857142857143_wp
    expected2(4,1) = 0.0_wp
    !
    expected2(1,2) = -0.3375_wp
    expected2(2,2) = -0.0982142857142857_wp
    expected2(3,2) = 0.4079365079365079_wp
    expected2(4,2) = 0.027777777777777777_wp
    !
    expected2(1,3) = -0.096_wp
    expected2(2,3) = -0.44571428571428573_wp
    expected2(3,3) = 0.38171428571428573_wp
    expected2(4,3) = 0.16_wp
    !
    expected2(1,4) = -0.2851431428571429_wp
    expected2(2,4) = -0.15974525714285698_wp
    expected2(3,4) = 0.44488799999999984_wp
    expected2(4,4) = 3.999999999999119e-07_wp
    !
    expected2(1,5) = -0.4_wp
    expected2(2,5) = 0.0_wp
    expected2(3,5) = 0.4_wp
    expected2(4,5) = 0.0_wp
    !
    expected2(1,6) = -0.04444444444444444_wp
    expected2(2,6) = -0.55555555555555555_wp
    expected2(3,6) = 0.35_wp
    expected2(4,6) = 0.25_wp
    !
    expected2(1,7) = 0.0_wp
    expected2(2,7) = -0.6_wp
    expected2(3,7) = 0.42857142857142866_wp
    expected2(4,7) = 0.17142857142857143_wp

    do j = 1, num_tests

      x = x_test(j)

      cell = sll_f_find_cell( bsplines, x )

      ! compute with non-uniform B-splines, old type
      call sll_s_bsplines_eval_basis( bsplines, cell, x, values_old )
      call sll_s_bsplines_eval_deriv( bsplines, cell, x, derivs_old )

      ! compute with non uniform B-splines, new interface
      call bsplines_non_uniform % eval_basis( x, values_new, jmin )
      call bsplines_non_uniform % eval_deriv( x, derivs_new, jmin )

      ! values: check difference with expected values
      passed_test = passed_test .and. &
        maxval( abs( values_old(:)-expected1(:,j) ) ) < tolerance .and. &
        maxval( abs( values_new(:)-expected1(:,j) ) ) < tolerance

      if( .not. passed_test ) then
        write(*,*) 'test_non_uniform_bsplines_periodic, values: wrong result for x = ', x
        write(*,*) 'non-uniform B-splines, old type, values: ', values_old
        write(*,*) 'non-uniform B-splines, new type, values: ', values_new
        write(*,*) 'non-uniform B-splines, expected results: ', expected1(:,j)
        write(*,*) 'exiting...'
        stop
      else
        write(*,'(a,i2,a,i2,a)') "test ", j, ", cell ", cell, ": comparing values      of B-splines   OK"
      end if

      ! derivs: check difference with expected values
      passed_test = passed_test .and. &
        maxval( abs( derivs_old(:)-expected2(:,j) ) ) < tolerance .and. &
        maxval( abs( derivs_new(:)-expected2(:,j) ) ) < tolerance

      if( .not. passed_test ) then
        write(*,*) 'test_non_uniform_bsplines_periodic, derivs: wrong result for x = ', x
        write(*,*) 'non-uniform B-splines, old type, values: ', derivs_old
        write(*,*) 'non-uniform B-splines, new type, values: ', derivs_new
        write(*,*) 'non-uniform B-splines, expected results: ', expected2(:,j)
        write(*,*) 'exiting...'
        stop
      else
        write(*,'(a)') "                  comparing derivatives of B-splines   OK"
      end if

    end do

    deallocate(values_old)
    deallocate(values_new)
    deallocate(derivs_old)
    deallocate(derivs_new)
    deallocate(values_and_deriv)

    call sll_s_bsplines_free( bsplines )
    call bsplines_non_uniform % free()

  end subroutine test_non_uniform_bsplines_periodic

  !----------------------------------------------------------------------------
  ! The case of 'open' boundary condition yields spline values different
  ! than in the 'periodic' case. Since we can not compare with the uniform
  ! splines anymore to get the right answer, we need a different criterion.
  ! For lack of something better, at the moment we only check the property of
  ! the partition of unity.
  subroutine test_non_uniform_bsplines_open( passed_test )
    logical, intent(inout) :: passed_test

    integer  :: i, j, jmin, cell
    integer  :: num_pts, num_tests
    integer  :: degree
    real(wp) :: x, xmin, rnd
    real(wp) :: step
    real(wp) :: tolerance
    real(wp), allocatable :: grid(:)
    real(wp), allocatable :: values_old(:)
    real(wp), allocatable :: values_new(:)
    real(wp), allocatable :: derivs_old(:)
    real(wp), allocatable :: derivs_new(:)
    real(wp), allocatable :: values_and_deriv(:,:)

    ! B-splines, old type
    type(sll_t_bsplines) :: bsplines
    ! Non-uniform B-splines, new interface
    type(sll_t_bsplines_non_uniform) :: bsplines_non_uniform

    num_tests = 10
    num_pts   = 10
    xmin      = 0.0_wp
    step      = 1.0_wp
    tolerance = 1.0d-15

    degree  = 3
    write(*,'(a)') "Test on random grid for given degree"
    write(*,'(a)') "------------------------------------"
    write(*,*)
    write(*,'(a,i0,a)') "degree = ", degree, ":"
    write(*,*)

    allocate(grid(num_pts))
    allocate(values_old(degree+1))
    allocate(values_new(degree+1))
    allocate(derivs_old(degree+1))
    allocate(derivs_new(degree+1))
    allocate(values_and_deriv(2,degree+1))

    ! create non-uniform mesh
    grid(1) = xmin
    do i = 2, num_pts
      call random_number(rnd)
      grid(i) = grid(i-1) + step + rnd
    end do

    ! initialize B-splines, old type
    call sll_s_bsplines_init_from_grid( bsplines, degree, grid, sll_p_open, sll_p_open )
    ! initialize non-uniform B-splines, new interface
    call bsplines_non_uniform % init( degree, periodic=.false., breaks=grid )

    do j = 1, num_tests

      call random_number(rnd)
      x = xmin + rnd*(grid(num_pts)-xmin)
      cell = sll_f_find_cell( bsplines, x )

      ! compute with non-uniform B-splines, old type
      call sll_s_bsplines_eval_basis( bsplines, cell, x, values_old )
      call sll_s_bsplines_eval_deriv( bsplines, cell, x, derivs_old )
      call sll_s_bsplines_eval_basis_and_deriv( bsplines, cell, x, values_and_deriv )

      ! compute with non uniform B-splines, new interface
      call bsplines_non_uniform % eval_basis( x, values_new, jmin )
      call bsplines_non_uniform % eval_deriv( x, derivs_new, jmin )

      ! test partition of unity
      passed_test = passed_test .and. &
        abs( 1.0_wp - sum( values_old(1:degree+1) ) ) < tolerance .and. &
        abs( 1.0_wp - sum( values_new(1:degree+1) ) ) < tolerance

      if( .not. passed_test ) then
        write(*,*) 'test_non_uniform_bsplines_open, values: wrong result for x = ', x
        write(*,*) 'non-uniform B-splines, old type, sum of values: ', sum( values_old(1:degree+1) )
        write(*,*) 'non-uniform B-splines, new type, sum of values: ', sum( values_new(1:degree+1) )
        write(*,*) 'exiting...'
        stop
      else
        write(*,'(a,i2,a,i2,a)') "test ", j, ", cell ", cell, ": testing partition of unity   OK"
      end if

      ! test sum of all derivatives (should be 0)
      passed_test = passed_test .and. &
        abs( sum( derivs_old(1:degree+1) ) ) < tolerance .and. &
        abs( sum( derivs_new(1:degree+1) ) ) < tolerance

      if( .not. passed_test ) then
        write(*,*) 'test_non_uniform_bsplines_open, derivs: wrong result for x = ', x
        write(*,*) 'non-uniform B-splines, old type, sum of derivs: ', sum( derivs_old(1:degree+1) )
        write(*,*) 'non-uniform B-splines, new type, sum of derivs: ', sum( derivs_new(1:degree+1) )
        write(*,*) 'exiting...'
        stop
      else
        write(*,'(a)') "                  testing sum of derivatives   OK"
      end if

      ! test values and derivatives
      passed_test = passed_test .and. &
        abs( 1.0_wp - sum( values_and_deriv(1,1:degree+1) ) ) < tolerance .and. &
        abs(          sum( values_and_deriv(2,1:degree+1) ) ) < tolerance

      if( .not. passed_test ) then
        write(*,*) 'test_non_uniform_bsplines_periodic, values and derivs: wrong result for x = ', x
        write(*,*) 'non-uniform B-splines, old type, values_and_deriv(1,:): ', values_and_deriv(1,:)
        write(*,*) 'non-uniform B-splines, old type, values_and_deriv(2,:): ', values_and_deriv(2,:)
        write(*,*) 'exiting...'
        stop
      else
        write(*,'(a)') "                  consistency checks           OK"
      end if

    end do

    deallocate(grid)
    deallocate(values_old)
    deallocate(values_new)
    deallocate(derivs_old)
    deallocate(derivs_new)
    deallocate(values_and_deriv)

    call sll_s_bsplines_free( bsplines )
    call bsplines_non_uniform % free()

  end subroutine test_non_uniform_bsplines_open

  !----------------------------------------------------------------------------
  subroutine test_performance

    real(wp), allocatable :: grid(:)
    integer  :: i, j, jmin
    integer  :: ncells, num_tests, num_derivs
    integer  :: degree
    real(wp) :: xmin, rnd

    integer , allocatable :: cells(:)
    real(wp), allocatable :: x(:)
    real(wp), allocatable :: xx(:)

    real(wp), allocatable :: values(:)
    real(wp), allocatable :: derivs(:)
    real(wp), allocatable :: values_and_deriv(:,:)
    real(wp), allocatable :: values_and_n_derivs(:,:)

    ! B-splines, old type
    type(sll_t_bsplines) :: bsplines
    ! Uniform and non-uniform B-splines, new interface
    type(sll_t_bsplines_uniform)     :: bsplines_uniform
    type(sll_t_bsplines_non_uniform) :: bsplines_non_uniform

    type(sll_t_time_mark) :: t0
    real(wp) :: time

    ! Test performance of nonuniform arbitrary degree spline evaluation
    num_tests  = 1000000
    ncells     = 11
    num_derivs = 2
    xmin       = 0.0_wp

    degree = 6
    write(*,'(a,i0,a)') "degree = ", degree, ":"
    write(*,*)

    allocate(x (num_tests))
    allocate(xx(num_tests))
    allocate(cells(num_tests))

    allocate(grid(ncells+1))

    allocate(values(degree+1))
    allocate(derivs(degree+1))
    allocate(values_and_deriv(2,degree+1))
    allocate(values_and_n_derivs(0:num_derivs,degree+1))

    !--------------------------------------------------------------------------
    ! Compare performance of non-uniform B-splines
    !--------------------------------------------------------------------------

    write(*,'(a)') "Compare performance of non-uniform B-splines"
    write(*,'(a)') "--------------------------------------------"
    write(*,*)

    ! create non-uniform mesh
    grid(1) = xmin
    do i = 2, ncells+1
       call random_number(rnd)
       grid(i) = grid(i-1) + rnd !step
    end do

    ! initialize B-splines, old type
    call sll_s_bsplines_init_from_grid( bsplines, degree, grid, sll_p_periodic, sll_p_periodic )
    ! initialize uniform B-splines, new interface
    call bsplines_uniform % init( degree, periodic=.true., xmin=grid(1), xmax=grid(ncells+1), ncells=ncells )
    ! initialize non-uniform B-splines, new interface
    call bsplines_non_uniform % init( degree, periodic=.true., breaks=grid )

    do j = 1, num_tests
      call random_number(rnd)
      x(j)  = xmin + rnd*(grid(ncells+1)-xmin)
      xx(j) = rnd
      cells(j) = sll_f_find_cell( bsplines, x(j) )
    end do

    !--------------------------------------------------------------------------
    ! computing all non zero splines at all points in x:
    call sll_s_set_time_mark(t0)
    do j = 1, num_tests
      call sll_s_bsplines_eval_basis( bsplines, cells(j), x(j), values )
    end do
    time = sll_f_time_elapsed_since(t0) / real(num_tests,wp)
    subr = "sll_s_bsplines_eval_basis"
    write(*,timing_fmt) adjustl( subr ), time

    call sll_s_set_time_mark(t0)
    do j = 1, num_tests
      call sll_s_bsplines_eval_basis_mm( bsplines%knots, cells(j), x(j), degree, values )
    end do
    time = sll_f_time_elapsed_since(t0) / real(num_tests,wp)
    subr = "sll_s_bsplines_eval_basis_mm"
    write(*,timing_fmt) adjustl( subr ), time

    ! computing with non-uniform B-splines, new interface
    call sll_s_set_time_mark( t0 )
    do j = 1, num_tests
      call bsplines_non_uniform % eval_basis( xx(j), values, jmin )
    end do
    time = sll_f_time_elapsed_since(t0) / real(num_tests,wp)
    subr = "sll_t_bsplines_non_uniform % eval_basis()"
    write(*,timing_fmt) adjustl( subr ), time
    write(*,*)

    !--------------------------------------------------------------------------
    ! computing both all non zero splines and derivatives at point x:
    call sll_s_set_time_mark(t0)
    do j = 1, num_tests
      call sll_s_bsplines_eval_basis_and_deriv( bsplines, cells(j), x(j), values_and_deriv )
    end do
    time = sll_f_time_elapsed_since(t0) / real(num_tests,wp)
    subr = "sll_s_bsplines_eval_basis_and_deriv"
    write(*,timing_fmt) adjustl( subr ), time

    call sll_s_set_time_mark(t0)
    do j = 1, num_tests
      call sll_s_bsplines_eval_basis_and_deriv_mm( bsplines%knots, cells(j), x(j), degree, values_and_deriv )
    end do
    time = sll_f_time_elapsed_since(t0) / real(num_tests,wp)
    subr = "sll_s_bsplines_eval_basis_and_deriv_mm"
    write(*,timing_fmt) adjustl( subr ), time
    write(*,*)

    !--------------------------------------------------------------------------
    ! computing all non zero spline derivatives at point x:
    call sll_s_set_time_mark(t0)
    do j = 1, num_tests
      call sll_s_bsplines_eval_deriv( bsplines, cells(j), x(j), derivs )
    end do
    time = sll_f_time_elapsed_since(t0) / real(num_tests,wp)
    subr = "sll_s_bsplines_eval_deriv"
    write(*,timing_fmt) adjustl( subr ), time

    ! computing with non-uniform B-splines, new interface
    call sll_s_set_time_mark( t0 )
    do j = 1, num_tests
      call bsplines_non_uniform % eval_deriv( xx(j), derivs, jmin )
    end do
    time = sll_f_time_elapsed_since(t0) / real(num_tests,wp)
    subr = "sll_t_bsplines_non_uniform % eval_deriv()"
    write(*,timing_fmt) adjustl( subr ), time
    write(*,*)

    !--------------------------------------------------------------------------
    ! computing both all non zero splines and all derivatives at point x:
    call sll_s_set_time_mark(t0)
    do j = 1, num_tests
      call sll_s_bsplines_eval_basis_and_n_derivs( bsplines, cells(j), x(j), 1, values_and_deriv )
    end do
    time = sll_f_time_elapsed_since(t0) / real(num_tests,wp)
    subr = "sll_s_bsplines_eval_basis_and_n_derivs"
    write(*,timing_fmt) adjustl( subr ), time

    ! computing with non-uniform B-splines, new interface
    call sll_s_set_time_mark( t0 )
    do j = 1, num_tests
      call bsplines_non_uniform % eval_basis_and_n_derivs( xx(j), num_derivs, values_and_n_derivs, jmin )
    end do
    time = sll_f_time_elapsed_since(t0) / real(num_tests,wp)
    subr = "sll_t_bsplines_non_uniform % eval_basis_and_n_derivs()"
    write(*,timing_fmt) adjustl( subr ), time
    write(*,*)

    !--------------------------------------------------------------------------
    ! Compare performance of uniform B-splines
    !--------------------------------------------------------------------------

    write(*,'(a)') "Compare performance of uniform B-splines"
    write(*,'(a)') "----------------------------------------"
    write(*,*)

    !--------------------------------------------------------------------------
    ! computing all non zero uniform splines at point x:
    call sll_s_set_time_mark(t0)
    do j = 1, num_tests
      call sll_s_uniform_bsplines_eval_basis( degree, xx(j), values )
    end do
    time = sll_f_time_elapsed_since(t0) / real(num_tests,wp)
    subr = "sll_s_uniform_bsplines_eval_basis"
    write(*,timing_fmt) adjustl( subr ), time

    ! computing with uniform B-splines, new interface
    call sll_s_set_time_mark( t0 )
    do j = 1, num_tests
      call bsplines_uniform % eval_basis( xx(j), values, jmin )
    end do
    time = sll_f_time_elapsed_since(t0) / real(num_tests,wp)
    subr = "sll_t_bsplines_uniform % eval_basis()"
    write(*,timing_fmt) adjustl( subr ), time
    write(*,*)

    !--------------------------------------------------------------------------
    ! computing all non zero uniform splines derivatives at point x:
    call sll_s_set_time_mark(t0)
    do j = 1, num_tests
      call sll_s_uniform_bsplines_eval_deriv( degree, xx(j), derivs )
    end do
    time = sll_f_time_elapsed_since(t0) / real(num_tests,wp)
    subr = "sll_s_uniform_bsplines_eval_deriv"
    write(*,timing_fmt) adjustl( subr ), time

    ! computing with uniform B-splines, new interface
    call sll_s_set_time_mark( t0 )
    do j = 1, num_tests
      call bsplines_uniform % eval_deriv( xx(j), derivs, jmin )
    end do
    time = sll_f_time_elapsed_since(t0) / real(num_tests,wp)
    subr = "sll_t_bsplines_uniform % eval_deriv()"
    write(*,timing_fmt) adjustl( subr ), time
    write(*,*)

    !--------------------------------------------------------------------------
    ! computing all non zero uniform splines and derivatives at point x:
    call sll_s_set_time_mark(t0)
    do j = 1, num_tests
      call sll_s_uniform_bsplines_eval_basis_and_deriv( degree, xx(j), values_and_deriv )
    end do
    time = sll_f_time_elapsed_since(t0) / real(num_tests,wp)
    subr = "sll_s_uniform_bsplines_eval_basis_and_deriv"
    write(*,timing_fmt) adjustl( subr ), time
    write(*,*)

    !--------------------------------------------------------------------------
    ! computing all non zero uniform splines and all derivatives at point x:
    call sll_s_set_time_mark(t0)
    do j = 1, num_tests
      call sll_s_uniform_bsplines_eval_basis_and_n_derivs( degree, xx(j), num_derivs, values_and_n_derivs )
    end do
    time = sll_f_time_elapsed_since(t0) / real(num_tests,wp)
    subr = "sll_s_uniform_bsplines_eval_basis_and_n_derivs"
    write(*,timing_fmt) adjustl( subr ), time

    ! computing with uniform B-splines, new interface
    call sll_s_set_time_mark( t0 )
    do j = 1, num_tests
      call bsplines_uniform % eval_basis_and_n_derivs( xx(j), num_derivs, values_and_n_derivs, jmin )
    end do
    time = sll_f_time_elapsed_since(t0) / real(num_tests,wp)
    subr = "sll_t_bsplines_uniform % eval_basis_and_n_derivs()"
    write(*,timing_fmt) adjustl( subr ), time
    write(*,*)

    deallocate(values)
    deallocate(derivs)
    deallocate(values_and_deriv   )
    deallocate(values_and_n_derivs)

    ! Free B-splines, old type
    call sll_s_bsplines_free( bsplines )
    ! Free uniform B-splines, new interface
    call bsplines_uniform % free()
    ! Free non-uniform B-splines, new interface
    call bsplines_non_uniform % free()

  end subroutine test_performance

end program test_arbitrary_degree_splines
