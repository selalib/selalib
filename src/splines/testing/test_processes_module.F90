module test_processes_module
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_hermite, &
    sll_p_periodic

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_cubic_splines, only: &
    sll_s_compute_cubic_spline_1d, &
    sll_s_compute_cubic_spline_2d, &
    sll_f_interpolate_derivative, &
    sll_s_interpolate_from_interpolant_array, &
    sll_f_interpolate_from_interpolant_value, &
    sll_f_interpolate_value_2d, &
    sll_f_new_cubic_spline_1d, &
    sll_f_new_cubic_spline_2d, &
    sll_t_cubic_spline_1d, &
    sll_t_cubic_spline_2d, &
    sll_o_delete

  use test_func_module, only: &
    f, &
    fprime

  use util_constants, only: &
    np, &
    npx1, &
    npx2, &
    r1, &
    r2, &
    x1max, &
    x1min, &
    x2max, &
    x2min, &
    xmax, &
    xmin

  implicit none

  public :: &
    coscos, &
    deriv1_polar_x, &
    deriv1_polar_y, &
    deriv2_polar_x, &
    deriv2_polar_y, &
    dmycos, &
    fxy, &
    interpolator_tester_1d_hrmt, &
    interpolator_tester_1d_prdc, &
    interpolator_tester_2d_hrmt_prdc, &
    interpolator_tester_2d_prdc_prdc, &
    line, &
    mcossin, &
    msincos, &
    mycos, &
    plane, &
    plane2, &
    plane2_deriv, &
    plane3, &
    plane3_deriv_x, &
    plane3_deriv_y, &
    plane_deriv, &
    polar_x, &
    polar_y, &
    sincos_prod, &
    sinsin, &
    test_2d_spline_hrmt_hrmt, &
    test_2d_spline_hrmt_hrmt_no_slopes, &
    test_2d_spline_hrmt_prdc, &
    test_2d_spline_hrmt_prdc_no_slopes, &
    test_2d_spline_prdc_hrmt, &
    test_2d_spline_prdc_hrmt_no_slopes, &
    test_spline_1d_hrmt

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  abstract interface
     function fx(x)
       use sll_m_working_precision
       sll_real64 :: fx
       sll_real64, intent(in) :: x
     end function fx
  end interface

  abstract interface
     function fxy (x,y)
       use sll_m_working_precision
       sll_real64 :: fxy
       sll_real64, intent(in) :: x
       sll_real64, intent(in) :: y
     end function fxy
  end interface

  abstract interface
     function spline_interpolator_1d(  x, spline )
       use sll_m_working_precision
       use sll_m_cubic_splines
       sll_real64 :: spline_interpolator_1d
       sll_real64, intent(in) :: x
       type(sll_t_cubic_spline_1d), pointer :: spline
     end function spline_interpolator_1d
  end interface


  abstract interface
     function spline_interpolator_2d(  x, y, spline )
       use sll_m_working_precision
       use sll_m_cubic_splines
       sll_real64 :: spline_interpolator_2d
       sll_real64, intent(in) :: x
       sll_real64, intent(in) :: y
       type(sll_t_cubic_spline_2d), pointer :: spline
     END function spline_interpolator_2d
  end interface

contains


  subroutine test_process_1d(i_test, ok)
    implicit none

    sll_int32                             :: i_test, ok, i, err 
    type(sll_t_cubic_spline_1d), pointer         :: sp1
    type(sll_t_cubic_spline_1d), pointer         :: sp2
    sll_real64                            :: x
    sll_real64                            :: phase
    sll_real64, allocatable, dimension(:) :: coordinates
    sll_real64, allocatable, dimension(:) :: data  
    sll_real64, allocatable, dimension(:) :: deriv
    sll_real64, allocatable, dimension(:) :: out
    sll_real64                            :: accumulator1, accumulator2  
    sll_real64                            :: accumulator3, accumulator4
    sll_real64                            :: accumulator5, accumulator6

    accumulator1 = 0.0_f64
    accumulator2 = 0.0_f64
    accumulator3 = 0.0_f64
    accumulator4 = 0.0_f64
    accumulator5 = 0.0_f64
    accumulator6 = 0.0_f64

    print *, 'Spline module unit tester'
    print *, 'allocate data array'
    
    SLL_ALLOCATE(data(np), err)
    SLL_ALLOCATE(deriv(np), err)
    SLL_ALLOCATE(out(np), err)
    SLL_ALLOCATE(coordinates(np), err)
    
    print *, 'initialize data and coordinates array'
    do i=1,np
       phase          = real(i-1,f64)*(xmax-xmin)/real(np-1,f64) + xmin
       data(i)        = f(phase, i_test)
       deriv(i)       = fprime(phase, i_test)
       coordinates(i) = phase
    end do
     
    print *, 'proceed to allocate the spline...'
    sp1 =>  sll_f_new_cubic_spline_1d( np, xmin, xmax, sll_p_periodic )
    call sll_s_compute_cubic_spline_1d( data, sp1 )
    sp2 =>  sll_f_new_cubic_spline_1d( &
         np, &
         xmin, &
         xmax, &
         sll_p_hermite, &
         fprime(xmin, i_test), &
         fprime(xmax, i_test) )
    call sll_s_compute_cubic_spline_1d( data, sp2 )
    
    print *, 'cumulative errors: '
    print *, 'periodic case, np-1 points: '
    print *, 'interpolating individual values from 1 to np-1:'
    do i=1, np-1
       x = real(i-1,f64)*(xmax-xmin)/real(np-1,f64)+xmin
       accumulator1 = accumulator1 + abs(data(i) - sll_f_interpolate_from_interpolant_value(x, sp1))
    end do
    print *, 'checking periodicity:'
    print *, 'difference between values at points 1 and np: ', &
         abs(data(1) - sll_f_interpolate_from_interpolant_value(xmax,sp1))
    print *, 'interpolating the whole array:'
    call sll_s_interpolate_from_interpolant_array(coordinates, out, np-1, sp1)
    do i=1, np-1
       accumulator3 = accumulator3 + abs(data(i) - out(i))
    end do
    
    print *, 'hermite case, np points: '
    do i=1, np
       x = real(i-1,f64)*(xmax-xmin)/real(np-1,f64) + xmin
       accumulator2 = accumulator2 + abs(data(i) - sll_f_interpolate_from_interpolant_value(x, sp2))
    end do
    call sll_s_interpolate_from_interpolant_array(coordinates, out, np, sp2)
    do i=1, np
       accumulator4 = accumulator4 + abs(data(i) - out(i))
    end do
     
    print *, '----------------------------------------------------'
    print *, 'RESULTS: '
    print *, 'Periodic case: '
    print *, 'average error at the nodes (single values) = '
    print *, accumulator1/real(np,f64)
    print *, 'average error at the nodes (whole array) = '
    print *, accumulator3/real(np,f64)
    write (*,'(a,f8.5)')   'original data(0)    = ', data(1)
    write (*,'(a,f20.15)') &
         'interpolated        = ', sll_f_interpolate_from_interpolant_value( xmin,sp1)
    
    write (*,'(a,f20.15)')   'original data((np-1)/4) = ', data((np-1)/4)
    write (*,'(a,f20.15)') &
         'interpolated        = ', sll_f_interpolate_from_interpolant_value( (xmax-xmin)/4.0+xmin,sp1)
     
    if ( (accumulator1/real(np,f64) >= 1.0d-14) .or. &
         (accumulator3/real(np,f64) >= 1.0d-14) ) then 
       ok = 0
       print*, 'i_test =', i_test
       print *, 'splines unit test stopped by periodic spline1d test failure'
       stop
    end if
    
    print *, '**************************** '
    print *, 'Hermite case: '
    print *, 'average error at the nodes = '
    print *, accumulator2/real(np,f64)
    write (*,'(a,f8.5)')   'original data(0)    = ', data(1)
    write (*,'(a,f20.15)') &
         'interpolated        = ', sll_f_interpolate_from_interpolant_value( 0.0_f64,sp2)
    write (*,'(a,f20.15)')   'original data((np-1)/4) = ', data((np-1)/4+1)
    write (*,'(a,f20.15)') &
         'interpolated        = ', sll_f_interpolate_from_interpolant_value( (xmax-xmin)/4.0,sp2)
     print *, 'spline coefficients: '
     
    if ( (accumulator2/real(np,f64) >= 1.0d-14) .or. &
         (accumulator4/real(np,f64) >= 1.0d-14) ) then
       ok = 0
       print*, 'i_test =', i_test
       print *, 'splines unit test stopped by Hermite slpine1d test failure'
       stop
    end if
    if (i_test<=12) then 
       print *, '---------------------------------------------'
       print *, 'DERIVATIVES TEST'
       do i=1, np
          phase = real(i-1,f64)*(xmax-xmin)/real(np-1,f64)+xmin
          accumulator5 = accumulator5 + abs(deriv(i) - &
                         sll_f_interpolate_derivative(phase, sp1))
          accumulator6 = accumulator6 + abs(deriv(i) - &
                         sll_f_interpolate_derivative(phase, sp2))
       end do
    
       if ( &
            (accumulator5/real(np,f64) > 1.0d-6) .or. &
            (accumulator6/real(np,f64) > 1.0d-6) ) then
          ok = 0
          print *, 'average error at the nodes (single values, periodic) = '
          print *, accumulator5/real(np,f64)
          print *, 'average error at the nodes (single values, hermite) = '
          print *, accumulator6/real(np,f64)
          print *, 'i_test =', i_test
          print *, 'splines unit test stopped by DERIVATIVES TEST failure'
          stop
       end if
    endif
     
    call sll_o_delete(sp1)
    call sll_o_delete(sp2)
    
    SLL_DEALLOCATE_ARRAY(data, err)
    SLL_DEALLOCATE_ARRAY(deriv, err)
    SLL_DEALLOCATE_ARRAY(out, err)
    SLL_DEALLOCATE_ARRAY(coordinates, err)
  
  end subroutine test_process_1d
    
  subroutine test_process_2d(i_test, j_test, ok)
    use test_func_module
    implicit none

    sll_int32                               :: i_test, j_test, ok
    sll_int32                               :: i, j, err
    type(sll_t_cubic_spline_2d), pointer            :: sp2d_1
    type(sll_t_cubic_spline_2d), pointer            :: sp2d_2
    type(sll_t_cubic_spline_2d), pointer            :: sp2d_3
    type(sll_t_cubic_spline_2d), pointer            :: sp2d_4
    sll_real64                              :: phase_x1, phase_x2
    sll_real64                              :: acc_2D, x1, x2
    sll_real64, allocatable, dimension(:,:) :: data_2d
    sll_real64, allocatable, dimension(:)   :: coordinates_i
    sll_real64, allocatable, dimension(:)   :: coordinates_j

    SLL_ALLOCATE(data_2d(npx1, npx2), err)        
    SLL_ALLOCATE(coordinates_i(npx1),err)
    SLL_ALLOCATE(coordinates_j(npx2),err)
        
    do i=1, npx1
       coordinates_i(i) = real(i-1,f64)*(x1max-x1min)/real(npx1-1,f64)+x1min
    end do
        
    do j=1, npx2
       coordinates_j(j) = real(j-1,f64)*(x2max-x2min)/real(npx2-1,f64)+x2min
    end do
        
    do j=1,npx2
       phase_x2 = coordinates_j(j)
       do i=1,npx1
          phase_x1 = coordinates_i(i)
          data_2d(i,j) = f(phase_x1, i_test)*f(phase_x2, j_test)
       end do
    end do
        
    ! Test the periodic-periodic spline2d 
    sp2d_1 => sll_f_new_cubic_spline_2d( npx1, npx2, &
         x1min, x1max, &
         x2min, x2max, &
         sll_p_periodic, sll_p_periodic )
    call sll_s_compute_cubic_spline_2d( data_2d, sp2d_1 )
    acc_2D = 0.0_f64
    do j=1, npx2
       do i=1, npx1
          x1 = coordinates_i(i)
          x2 = coordinates_j(j)
          acc_2D = acc_2D + &
               abs(data_2d(i,j) - sll_f_interpolate_value_2d(x1,x2,sp2d_1))
       end do
    end do
    if (acc_2D/(npx1*npx2)>=1.d-14) then
       ok = 0
       print*, '(i_test, j_test) = (', i_test, ',', j_test,')' 
       print *, 'Average cumulative error, spline2d, periodic-periodic: ', &
         acc_2D/(npx1*npx2)
       print *, 'splines unit test stopped by periodic-periodic spline2d test failure'
       stop
    endif
    ! Test the Hermite-periodic spline2d
    !print *, 'Testing hermite-periodic spline. Test 1'
    sp2d_2 => sll_f_new_cubic_spline_2d( npx1, npx2, &
         x1min, x1max, &
         x2min, x2max, &
         sll_p_hermite, sll_p_periodic, &
         const_slope_x1_min = fprime(x1min, i_test), &
         const_slope_x1_max = fprime(x1max, i_test) )
    call sll_s_compute_cubic_spline_2d( data_2d, sp2d_2 )
    acc_2D = 0.0_f64
    do j=1, npx2
       do i=1, npx1
          x1 = coordinates_i(i)
          x2 = coordinates_j(j)
          acc_2D = acc_2D + &
               abs(data_2d(i,j) - sll_f_interpolate_value_2d(x1,x2,sp2d_2))
       end do
    end do
    if (acc_2D/(npx1*npx2)>=1.d-14) then
       ok = 0
       print*, '(i_test, j_test) = (', i_test, ',', j_test,')' 
       print *, 'Average cumulative error, spline2d, Hermite-periodic: ', &
         acc_2D/(npx1*npx2)
       print *, 'splines unit test stopped by Hermite-periodic spline2d test failure'
       stop
    endif
    ! Test the periodic-Hermite spline2d
    sp2d_3 => sll_f_new_cubic_spline_2d( npx1, npx2, &
         x1min, x1max, &
         x2min, x2max, &
         sll_p_periodic, sll_p_hermite, &
         const_slope_x2_min = fprime(x2min, j_test), &
         const_slope_x2_max = fprime(x2max, j_test) )
    call sll_s_compute_cubic_spline_2d( data_2d, sp2d_3 )
    acc_2D = 0.0_f64
    do j=1, npx2
       do i=1, npx1
          x1 = coordinates_i(i)
          x2 = coordinates_j(j)
          acc_2D = acc_2D + &
               abs(data_2d(i,j) - sll_f_interpolate_value_2d(x1,x2,sp2d_3))
       end do
    end do
    if (acc_2D/(npx1*npx2)>=1.d-14) then
       ok = 0
       print*, '(i_test, j_test) = (', i_test, ',', j_test,')' 
       print *, 'Average cumulative error, spline2d, periodic-Hermite: ', &
         acc_2D/(npx1*npx2)
       print *, 'splines unit test stopped by periodic-Hermite spline2d test failure'
       stop
    endif
    ! Test the Hermite-Hermite spline2d
    sp2d_4 => sll_f_new_cubic_spline_2d( npx1, npx2, &
         x1min, x1max, &
         x2min, x2max, &
         sll_p_hermite, sll_p_hermite,&
         fprime(x1min, i_test), fprime(x1max, i_test), &
         fprime(x2min, j_test), fprime(x2max, j_test) )   
    call sll_s_compute_cubic_spline_2d( data_2d, sp2d_4 )
    acc_2D = 0.0_f64
    do j=1, npx2
       do i=1, npx1
          x1 = coordinates_i(i)
          x2 = coordinates_j(j)
          acc_2D = acc_2D + &
               abs(data_2d(i,j) - sll_f_interpolate_value_2d(x1,x2,sp2d_4))
       end do
    end do
    if (acc_2D/(npx1*npx2)>=1.d-14) then
       ok = 0
       print*, '(i_test, j_test) = (', i_test, ',', j_test,')' 
       print *, 'Average cumulative error, spline2d, Hermite-Hermite: ', &
         acc_2D/(npx1*npx2)
       print *, 'splines unit test stopped by Hermite-Hermite spline2d test failure'
       stop
    endif
     
    call sll_o_delete(sp2d_1)
    call sll_o_delete(sp2d_2)
    call sll_o_delete(sp2d_3)
    call sll_o_delete(sp2d_4)
    
    SLL_DEALLOCATE_ARRAY(data_2d, err)
    SLL_DEALLOCATE_ARRAY(coordinates_i,err)
    SLL_DEALLOCATE_ARRAY(coordinates_j,err)

  end subroutine test_process_2d

  function sincos_prod( x, y ) result(val)
    sll_real64 :: val
    sll_real64, intent(in) :: x, y
    val = sin(x)*cos(y)
  end function sincos_prod

  !*********************************************************************
  ! The above routines should be converted into individualized routines
  ! like the ones below, designed to test individual utilities.
  !
  !*********************************************************************

  subroutine interpolator_tester_1d_prdc( &
    func, &           ! function to produce the data
    result_f, &       ! function to produce the right result
    interpolator_f, & ! function used to interpolate the data
    xmin, &           ! extent of the domain over where the spline is tested
    xmax, &
    npts, &
    test_passed, &
    criterion )
    procedure(fx)                     :: func
    procedure(fx)                     :: result_f
    procedure(spline_interpolator_1d) :: interpolator_f
    sll_real64, intent(in)            :: xmin
    sll_real64, intent(in)            :: xmax
    sll_int32, intent(in)             :: npts
    logical, intent(out)              :: test_passed
    sll_real64, intent(in), optional  :: criterion
    sll_real64                        :: max_tolerable_err
    sll_real64, allocatable, dimension(:) :: data_in
    sll_real64, allocatable, dimension(:) :: correct_data_out
    sll_int32  :: ierr, i
    sll_real64 :: h1
    sll_real64 :: x1
    sll_real64 :: acc, val
    type(sll_t_cubic_spline_1d), pointer :: spline
    sll_real64 :: average_error
 
    if(present(criterion)) then
       max_tolerable_err = criterion
    else
       max_tolerable_err = 1.0d-14
    end if
    h1 = (xmax - xmin)/real(npts-1,f64)
    acc = 0.0_f64

    ! allocate arrays and initialize them
    SLL_ALLOCATE(data_in(npts),ierr)
    SLL_ALLOCATE(correct_data_out(npts), ierr)

    do i=0,npts-1
       x1 = xmin + real(i,f64)*h1
       data_in(i+1) = func(x1)
       correct_data_out(i+1) = result_f(x1)
    end do

    spline => sll_f_new_cubic_spline_1d( &
      npts, &
      xmin, &
      xmax, &
      sll_p_periodic )

    call sll_s_compute_cubic_spline_1d(data_in, spline)
    do i=0,npts-1
       x1 = xmin + real(i,f64)*h1 
       val = interpolator_f(x1,spline)
       acc = acc + abs(val-correct_data_out(i+1))
       !          print *, '(i) = ',i+1, 'correct value = ', &
       !               correct_data_out(i+1), '. Calculated = ', val   
    end do
    average_error = acc/(real(npts,f64))
    print *, 'Average error = ', average_error
    if( average_error .le. max_tolerable_err ) then
       test_passed = .true.
    else
       test_passed = .false.
    end if
    ! deallocate memory
    call sll_o_delete(spline)
    SLL_DEALLOCATE_ARRAY(data_in, ierr)
    SLL_DEALLOCATE_ARRAY(correct_data_out, ierr)
  end subroutine interpolator_tester_1d_prdc


  subroutine interpolator_tester_1d_hrmt( &
    func, &           ! function to produce the data
    result_f, &       ! function to produce the right result
    interpolator_f, & ! function used to interpolate the data
    xmin, &           ! extent of the domain over where the spline is tested
    xmax, &
    test_passed, &
    criterion )
    procedure(fx)                     :: func
    procedure(fx)                     :: result_f
    procedure(spline_interpolator_1d) :: interpolator_f
    sll_real64, intent(in)            :: xmin
    sll_real64, intent(in)            :: xmax
    logical, intent(out)              :: test_passed
    sll_int32, parameter              :: npts_min = 9
    sll_int32, parameter              :: npts_max = 129
    sll_int32                         :: npts
    sll_real64, intent(in), optional  :: criterion
    sll_real64                        :: max_tolerable_err
    sll_real64, allocatable, dimension(:) :: data_in
    sll_real64, allocatable, dimension(:) :: correct_data_out
    sll_int32  :: ierr, i
    sll_real64 :: h1
    sll_real64 :: x1
    sll_real64 :: acc, val
    type(sll_t_cubic_spline_1d), pointer :: spline
    sll_real64 :: average_error
 
    if(present(criterion)) then
       max_tolerable_err = criterion
    else
       max_tolerable_err = 1.0d-14
    end if

    test_passed = .true.

    do npts=npts_min, npts_max
       h1 = (xmax - xmin)/real(npts-1,f64)
       acc = 0.0_f64

       ! allocate arrays and initialize them
       SLL_CLEAR_ALLOCATE(data_in(1:npts),ierr)
       SLL_CLEAR_ALLOCATE(correct_data_out(1:npts), ierr)

       do i=0,npts-1
          x1 = xmin + real(i,f64)*h1
          data_in(i+1) = func(x1)
          correct_data_out(i+1) = result_f(x1)
       end do

       spline => sll_f_new_cubic_spline_1d( &
            npts, &
            xmin, &
            xmax, &
            sll_p_hermite )

       call sll_s_compute_cubic_spline_1d(data_in, spline)
       do i=0,npts-2
          x1 = xmin + real(i,f64)*h1 
          val = interpolator_f(x1,spline)
          acc = acc + abs(val-correct_data_out(i+1))
          !print *, '(i) = ',i+1, 'correct value = ', &
          !correct_data_out(i+1), '. Calculated = ', val, 'delta = ', h1, 'delta^4 = ', h1**4
       end do
       ! Do the last point separately because due to roundoff, the last value
       ! is outside of the specified domain.
       val = interpolator_f(xmax, spline)
       acc = acc + abs(val-correct_data_out(i+1))
       !print *, '(i) = ',npts, 'correct value = ', correct_data_out(npts), &
       !     '. Calculated = ', val, 'delta = ', h1, 'delta^4 = ', h1**4
       average_error = acc/(real(npts,f64))
       print *, 'Test for num. points: ', npts,'Average error = ', average_error
       if( average_error .le. h1**4 ) then
          test_passed = test_passed .and. .true.
       else
          test_passed = .false.
          print *, 'Failure in case: num. points = ', npts, 'Average error = ', &
               average_error, 'delta^4 = ', h1**4
       end if
       ! deallocate memory
       call sll_o_delete(spline)
       SLL_DEALLOCATE_ARRAY(data_in, ierr)
       SLL_DEALLOCATE_ARRAY(correct_data_out, ierr)
    end do
  end subroutine interpolator_tester_1d_hrmt



  ! Magic numbers inside this function should be converted to routine parameters
  subroutine test_spline_1d_hrmt ( &
    func_1d, &
    test_passed )
    ! This function does not explicitly specifies the slopes since these are
    ! computed by default inside the splines.
    procedure(fx)          :: func_1d
    logical, intent(out)   :: test_passed
    logical    :: local_test_passed
    sll_real64, allocatable, dimension(:) :: data_in
    sll_int32  :: ierr, i
    sll_real64 :: h1
    sll_real64 :: x1
    sll_real64 :: acc, val
    type(sll_t_cubic_spline_1d), pointer :: spline
    sll_real64 :: average_error
    sll_int32, parameter :: np_min = 9
    sll_int32, parameter :: np_max = 129
    sll_int32            :: npts

    local_test_passed = .true.
    print *, 'x1min, x1max = ', x1min, x1max
    ! Run the test over a range of data sizes.
    print*, np_min , np_max
    do npts=np_min, np_max
       h1 = (x1max - x1min)/real(npts-1,f64)
       acc = 0.0_f64
       print*, npts
       SLL_ALLOCATE(data_in(npts),ierr)
       do i=0,npts-1
          x1 = x1min + real(i,f64)*h1 
          data_in(i+1) = func_1d(x1)
       end do
       spline => sll_f_new_cubic_spline_1d( &
            npts, &
            x1min, &
            x1max, &
            sll_p_hermite )
       
       call sll_s_compute_cubic_spline_1d( data_in, spline )
       acc = 0.0_f64
       do i=0,npts-2 ! last point excluded and done separately...
          x1 = x1min + real(i,f64)*h1 
          val = sll_f_interpolate_from_interpolant_value(x1,spline)
          !print *,'x = ', x1, 'true data: ',data_in(i+1), 'interpolated: ', val
          acc = acc + abs(val-data_in(i+1))  
       end do
       ! Do the last point separately because due to roundoff error, it ends
       ! up out of the range inside the loop.
       val = sll_f_interpolate_from_interpolant_value(x1max, spline)
       acc = acc + abs(val-data_in(npts))    

       average_error = acc/(real(npts,f64))
       print *, 'test_spline_1d_hrmt(): average error = ', average_error, &
            'problem size = ', npts, 'points.'
       if( average_error .le. 5.0d-14 ) then
          local_test_passed = local_test_passed .and. .true.
       else
          local_test_passed = .false.
          print *, 'test_spline_1d_hrmt(): TEST FAILED'
          print *, 'problem size: ', npts, 'average error = ', average_error
       end if
       SLL_DEALLOCATE_ARRAY(data_in, ierr)
    end do
    test_passed = local_test_passed
  end subroutine test_spline_1d_hrmt


  subroutine interpolator_tester_2d_prdc_prdc( &
    func_2d, &            ! function to generate input data array
    partial_x1, &         ! function to generate 'correct answer' array
    interpolation_func, & ! function to test
    test_passed )

    procedure(fxy)                          :: func_2d   
    procedure(fxy)                          :: partial_x1 
    procedure(spline_interpolator_2d)       :: interpolation_func
    logical, intent(out)                    :: test_passed
    sll_real64, allocatable, dimension(:,:) :: data_in
    sll_real64, allocatable, dimension(:,:) :: correct_data_out
    sll_int32  :: ierr, i, j
    sll_real64 :: h1, h2 ! cell spacings
    sll_real64 :: x1, x2
    sll_real64 :: acc, val
    type(sll_t_cubic_spline_2d), pointer :: spline
    sll_real64 :: average_error
    h1 = (x1max - x1min)/real(npx1-1,f64)
    h2 = (x2max - x2min)/real(npx2-1,f64)
    acc = 0.0_f64

    ! allocate arrays and initialize them
    SLL_ALLOCATE(data_in(npx1,npx2),ierr)
    SLL_ALLOCATE(correct_data_out(npx1,npx2), ierr)
    do j=0,npx2-1
       do i=0,npx1-1
          x1 = x1min + real(i,f64)*h1 
          x2 = x2min + real(j,f64)*h2
          data_in(i+1,j+1) = func_2d(x1,x2)
          correct_data_out(i+1,j+1) = partial_x1(x1,x2)
       end do
    end do

    spline => sll_f_new_cubic_spline_2d( &
      npx1, &
      npx2, &
      x1min, &
      x1max, &
      x2min, &
      x2max, &
      sll_p_periodic, &
      sll_p_periodic )

    call sll_s_compute_cubic_spline_2d(data_in,spline)

    do j=0,npx2-1
       do i=0,npx1-1
          x1 = x1min + real(i,f64)*h1 
          x2 = x2min + real(j,f64)*h2
          val = interpolation_func(x1,x2,spline)
          acc = acc + abs(val-correct_data_out(i+1,j+1))  
 !         print *, '(i,j) = ',i+1,j+1, 'correct value = ', &
 !              correct_data_out(i+1,j+1), '. Calculated = ', val     
       end do
    end do
    average_error = acc/(real(npx1*npx2,f64))
    print *, 'Average error = ', average_error
    if( average_error .le. 1.0d-5 ) then
       test_passed = .true.
    else
       test_passed = .false.
    end if
  end subroutine interpolator_tester_2d_prdc_prdc

  subroutine interpolator_tester_2d_hrmt_prdc( &
    func_2d, &            ! function to generate input data array
    transformed_func, &   ! function to generate 'correct answer' array
    interpolation_func, & ! interpolation function to test
    slope_min_func, &
    slope_max_func, &
    test_passed )

    procedure(fxy)                          :: func_2d   
    procedure(fxy)                          :: transformed_func
    procedure(spline_interpolator_2d)          :: interpolation_func
    procedure(fxy)                          :: slope_min_func
    procedure(fxy)                          :: slope_max_func
    logical, intent(out)                    :: test_passed
    sll_real64, allocatable, dimension(:,:) :: data_in
    sll_real64, allocatable, dimension(:,:) :: correct_data_out
    sll_real64, allocatable, dimension(:)   :: slopes_min, slopes_max
    sll_int32  :: ierr, i, j, imax, jmax
    sll_real64 :: h1, h2 ! cell spacings
    sll_real64 :: x1, x2
    sll_real64 :: acc, val, err, min_err, max_err
    type(sll_t_cubic_spline_2d), pointer :: spline
    sll_real64 :: average_error
    h1  = 1.0_f64/real(npx1-1,f64)
    h2  = 1.0_f64/real(npx2-1,f64)
    acc = 0.0_f64

    max_err = 0.0_f64
    min_err = 1.0_f64
    ! allocate arrays and initialize them
    SLL_ALLOCATE(data_in(npx1,npx2),ierr)
    SLL_ALLOCATE(correct_data_out(npx1,npx2), ierr)
    SLL_ALLOCATE(slopes_min(npx2), ierr)
    SLL_ALLOCATE(slopes_max(npx2), ierr)

    do j=0,npx2-1
       do i=0,npx1-1
          x1 = real(i,f64)*h1 
          x2 = real(j,f64)*h2
          data_in(i+1,j+1) = func_2d(x1,x2)
          correct_data_out(i+1,j+1) = transformed_func(x1,x2)
       end do
    end do

    do j=0,npx2-1
       x2 = real(j,f64)*h2
       slopes_min(j+1) = slope_min_func(0.0_f64,x2)
       slopes_max(j+1) = slope_max_func(1.0_f64,x2)
    end do

    spline => sll_f_new_cubic_spline_2d( &
      npx1, &
      npx2, &
      0.0_f64, &
      1.0_f64, &
      0.0_f64, &
      1.0_f64, &
      sll_p_hermite, &
      sll_p_periodic, &
      x1_min_slopes=slopes_min, &
      x1_max_slopes=slopes_max )

    call sll_s_compute_cubic_spline_2d(data_in,spline)

    do j=0,npx2-1
       do i=0,npx1-1
          x1 = real(i,f64)*h1 
          x2 = real(j,f64)*h2
          val = interpolation_func(x1,x2,spline)
          err = abs(val-correct_data_out(i+1,j+1))  
          acc = acc + err
          if( err .gt. max_err ) then
             max_err = err
             imax = i+1
             jmax = j+1
          end if
          if( err .lt. min_err ) then
             min_err = err
          end if
!          print *, '(i,j) = ',i+1,j+1, 'correct value = ', &
!               correct_data_out(i+1,j+1), '. Calculated = ', val     
       end do
    end do
    average_error = acc/(real(npx1*npx2,f64))
    print *, 'interpolator_tester_2d_hrmt_prdc results. Average error = ', &
         average_error, 'Max err = ', max_err, 'Min err = ', min_err
    print *, 'max error at (i,j) = ', imax, jmax
    if( average_error .le. 1.0d-5 ) then
       test_passed = .true.
    else
       test_passed = .false.
    end if
  end subroutine interpolator_tester_2d_hrmt_prdc


  subroutine test_2d_spline_hrmt_prdc( &
    transform_func, &
    eta1_min_slope_func, &
    eta1_max_slope_func, &
    test_passed )

    procedure(fxy) :: transform_func
    procedure(fxy) :: eta1_min_slope_func
    procedure(fxy) :: eta1_max_slope_func
    logical, intent(out) :: test_passed

    sll_int32 :: i, j, ierr
    sll_real64, dimension(:,:), allocatable :: data
    sll_real64, dimension(:), allocatable   :: eta1_min_slopes
    sll_real64, dimension(:), allocatable   :: eta1_max_slopes
    sll_real64 :: h1, h2, eta1, eta2, acc, val, true_val, ave_err
    type(sll_t_cubic_spline_2d), pointer :: spline

    h1 = 1.0_f64/(npx1-1)
    h2 = 1.0_f64/(npx2-1)
    ! allocate and fill out the data
    SLL_ALLOCATE(data(npx1,npx2),ierr)
    do j=0,npx2-1
       eta2 = real(j,f64)*h2
       do i=0,npx1-1
          eta1 = real(i,f64)*h1
          data(i+1,j+1) = transform_func(eta1,eta2)
       end do
    end do

    ! allocate and fill out the bc data
    SLL_ALLOCATE(eta1_min_slopes(npx2),ierr)
    SLL_ALLOCATE(eta1_max_slopes(npx2),ierr)
    do j=0,npx2-1
       eta2 = real(j,f64)*h2
       eta1_min_slopes(j+1) = eta1_min_slope_func(0.0_f64,eta2)
       eta1_max_slopes(j+1) = eta1_max_slope_func(1.0_f64,eta2)
    end do

    spline =>sll_f_new_cubic_spline_2d( &
         npx1, &
         npx2, &
         0.0_f64, &
         1.0_f64, &
         0.0_f64, &
         1.0_f64, &
         sll_p_hermite, &
         sll_p_periodic, &
         x1_min_slopes=eta1_min_slopes, &
         x1_max_slopes=eta1_max_slopes )

    call sll_s_compute_cubic_spline_2d( data, spline )

    ! compare results
    acc = 0.0_f64
    do j=0,npx2-1
       eta2 = real(j,f64)*h2
       do i=0,npx1-1
          eta1     = real(i,f64)*h1
          val      = sll_f_interpolate_value_2d( eta1, eta2, spline )
          true_val = transform_func( eta1, eta2 )
          acc      = acc + abs(true_val - val)
          if(abs(true_val-val).gt.1.0d-14) then
             print *, 'i,j = ',i,j, 'eta1,eta2 = ', eta1, eta2, 'true, val = ',&
                  true_val, val
             !print *, spline%coeffs(npx1,j-1:j+1)
          end if
       end do
    end do
    ave_err = acc/real(npx1*npx2,f64) 
    print *, 'test_2d_spline_hrmt_prdc results. Average error = ', ave_err
    if( ave_err .le. 1.0d-14 ) then
       test_passed = .true.
    else
       test_passed = .false.
    end if
    SLL_DEALLOCATE_ARRAY(data,ierr)
    SLL_DEALLOCATE_ARRAY(eta1_min_slopes,ierr)
    SLL_DEALLOCATE_ARRAY(eta1_max_slopes,ierr)
  end subroutine test_2d_spline_hrmt_prdc

  ! same as test_2d_spline_hrmt_prdc but internally it does not pass the
  ! slope values to the spline, thus testing the default values.
  subroutine test_2d_spline_hrmt_prdc_no_slopes( &
    transform_func, &
    eta1_min_slope_func, &
    eta1_max_slope_func, &
    test_passed )

    procedure(fxy) :: transform_func
    procedure(fxy) :: eta1_min_slope_func
    procedure(fxy) :: eta1_max_slope_func
    logical, intent(out) :: test_passed

    sll_int32 :: i, j, ierr
    sll_real64, dimension(:,:), allocatable :: data
    sll_real64, dimension(:), allocatable   :: eta1_min_slopes
    sll_real64, dimension(:), allocatable   :: eta1_max_slopes
    sll_real64 :: h1, h2, eta1, eta2, acc, val, true_val, ave_err
    type(sll_t_cubic_spline_2d), pointer :: spline

    h1 = 1.0_f64/(npx1-1)
    h2 = 1.0_f64/(npx2-1)
    ! allocate and fill out the data
    SLL_ALLOCATE(data(npx1,npx2),ierr)
    do j=0,npx2-1
       eta2 = real(j,f64)*h2
       do i=0,npx1-1
          eta1 = real(i,f64)*h1
          data(i+1,j+1) = transform_func(eta1,eta2)
       end do
    end do

    ! allocate and fill out the bc data
    SLL_ALLOCATE(eta1_min_slopes(npx2),ierr)
    SLL_ALLOCATE(eta1_max_slopes(npx2),ierr)
    do j=0,npx2-1
       eta2 = real(j,f64)*h2
       eta1_min_slopes(j+1) = eta1_min_slope_func(0.0_f64,eta2)
       eta1_max_slopes(j+1) = eta1_max_slope_func(1.0_f64,eta2)
    end do

    spline =>sll_f_new_cubic_spline_2d( &
         npx1, &
         npx2, &
         0.0_f64, &
         1.0_f64, &
         0.0_f64, &
         1.0_f64, &
         sll_p_hermite, &
         sll_p_periodic ) !, &
!         x1_min_slopes=eta1_min_slopes, &
!         x1_max_slopes=eta1_max_slopes )

    call sll_s_compute_cubic_spline_2d( data, spline )

    ! compare results
    acc = 0.0_f64
    do j=0,npx2-1
       eta2 = real(j,f64)*h2
       do i=0,npx1-1
          eta1     = real(i,f64)*h1
          val      = sll_f_interpolate_value_2d( eta1, eta2, spline )
          true_val = transform_func( eta1, eta2 )
          acc      = acc + abs(true_val - val)
          if(abs(true_val-val).gt.1.0d-14) then
             print *, 'i,j = ',i,j, 'eta1,eta2 = ', eta1, eta2, 'true, val = ',&
                  true_val, val
             !print *, spline%coeffs(npx1,j-1:j+1)
          end if
       end do
    end do
    ave_err = acc/real(npx1*npx2,f64) 
    print *, 'test_2d_spline_hrmt_prdc results. Average error = ', ave_err
    if( ave_err .le. 1.0d-14 ) then
       test_passed = .true.
    else
       test_passed = .false.
    end if
    SLL_DEALLOCATE_ARRAY(data,ierr)
    SLL_DEALLOCATE_ARRAY(eta1_min_slopes,ierr)
    SLL_DEALLOCATE_ARRAY(eta1_max_slopes,ierr)
  end subroutine test_2d_spline_hrmt_prdc_no_slopes



  subroutine test_2d_spline_prdc_hrmt( &
    transform_func, &
    eta2_min_slope_func, &
    eta2_max_slope_func, &
    test_passed )

    procedure(fxy) :: transform_func
    procedure(fxy) :: eta2_min_slope_func
    procedure(fxy) :: eta2_max_slope_func
    logical, intent(out) :: test_passed

    sll_int32 :: i, j, ierr
    sll_real64, dimension(:,:), allocatable :: data
    sll_real64, dimension(:), allocatable   :: eta2_min_slopes
    sll_real64, dimension(:), allocatable   :: eta2_max_slopes
    sll_real64 :: h1, h2, eta1, eta2, acc, val, true_val, ave_err
    type(sll_t_cubic_spline_2d), pointer :: spline

    h1 = 1.0_f64/(npx1-1)
    h2 = 1.0_f64/(npx2-1)
    ! allocate and fill out the data
    SLL_ALLOCATE(data(npx1,npx2),ierr)
    do j=0,npx2-1
       eta2 = real(j,f64)*h2
       do i=0,npx1-1
          eta1 = real(i,f64)*h1
          data(i+1,j+1) = transform_func(eta1,eta2)
       end do
    end do

    ! allocate and fill out the bc data
    SLL_ALLOCATE(eta2_min_slopes(npx1),ierr)
    SLL_ALLOCATE(eta2_max_slopes(npx1),ierr)
    do i=0,npx1-1
       eta1 = real(i,f64)*h1
       eta2 = 0.0_f64
       eta2_min_slopes(i+1) = eta2_min_slope_func(eta1,eta2)
       eta2_max_slopes(i+1) = eta2_max_slope_func(eta1,eta2)
    end do

    spline =>sll_f_new_cubic_spline_2d( &
         npx1, &
         npx2, &
         0.0_f64, &
         1.0_f64, &
         0.0_f64, &
         1.0_f64, &
         sll_p_periodic, &
         sll_p_hermite, & 
         x2_min_slopes=eta2_min_slopes, &
         x2_max_slopes=eta2_max_slopes )

    call sll_s_compute_cubic_spline_2d( data, spline )

    ! compare results
    acc = 0.0_f64
    do j=0,npx2-1
       eta2 = real(j,f64)*h2
       do i=0,npx1-1
          eta1     = real(i,f64)*h1
          val      = sll_f_interpolate_value_2d( eta1, eta2, spline )
          true_val = transform_func( eta1, eta2 )
          acc      = acc + abs(true_val - val)
          if(abs(true_val-val).gt.1.0d-14) then
             print *, 'i,j = ',i,j, 'eta1,eta2 = ', eta1, eta2, 'true, val = ',&
                  true_val, val
             !print *, spline%coeffs(npx1,j-1:j+1)
          end if
       end do
    end do
    ave_err = acc/real(npx1*npx2,f64) 
    print *, 'test_2d_spline_prdc_hrmt results. Average error = ', ave_err
    if( ave_err .le. 1.0d-14 ) then
       test_passed = .true.
    else
       test_passed = .false.
    end if
    SLL_DEALLOCATE_ARRAY(data,ierr)
    SLL_DEALLOCATE_ARRAY(eta2_min_slopes,ierr)
    SLL_DEALLOCATE_ARRAY(eta2_max_slopes,ierr)
  end subroutine test_2d_spline_prdc_hrmt

  ! like test_2d_spline_prdc_hrmt but tests the default slope values
  subroutine test_2d_spline_prdc_hrmt_no_slopes( &
    transform_func, &
    eta2_min_slope_func, &
    eta2_max_slope_func, &
    test_passed )

    procedure(fxy) :: transform_func
    procedure(fxy) :: eta2_min_slope_func
    procedure(fxy) :: eta2_max_slope_func
    logical, intent(out) :: test_passed

    sll_int32 :: i, j, ierr
    sll_real64, dimension(:,:), allocatable :: data
    sll_real64, dimension(:), allocatable   :: eta2_min_slopes
    sll_real64, dimension(:), allocatable   :: eta2_max_slopes
    sll_real64 :: h1, h2, eta1, eta2, acc, val, true_val, ave_err
    type(sll_t_cubic_spline_2d), pointer :: spline

    h1 = 1.0_f64/(npx1-1)
    h2 = 1.0_f64/(npx2-1)
    ! allocate and fill out the data
    SLL_ALLOCATE(data(npx1,npx2),ierr)
    do j=0,npx2-1
       eta2 = real(j,f64)*h2
       do i=0,npx1-1
          eta1 = real(i,f64)*h1
          data(i+1,j+1) = transform_func(eta1,eta2)
       end do
    end do

    ! allocate and fill out the bc data
    SLL_ALLOCATE(eta2_min_slopes(npx1),ierr)
    SLL_ALLOCATE(eta2_max_slopes(npx1),ierr)
    do i=0,npx1-1
       eta1 = real(i,f64)*h1
       eta2 = 0.0_f64
       eta2_min_slopes(i+1) = eta2_min_slope_func(eta1,eta2)
       eta2_max_slopes(i+1) = eta2_max_slope_func(eta1,eta2)
    end do

    spline =>sll_f_new_cubic_spline_2d( &
         npx1, &
         npx2, &
         0.0_f64, &
         1.0_f64, &
         0.0_f64, &
         1.0_f64, &
         sll_p_periodic, &
         sll_p_hermite ) !, & 
!         x2_min_slopes=eta2_min_slopes, &
!         x2_max_slopes=eta2_max_slopes )

    call sll_s_compute_cubic_spline_2d( data, spline )

    ! compare results
    acc = 0.0_f64
    do j=0,npx2-1
       eta2 = real(j,f64)*h2
       do i=0,npx1-1
          eta1     = real(i,f64)*h1
          val      = sll_f_interpolate_value_2d( eta1, eta2, spline )
          true_val = transform_func( eta1, eta2 )
          acc      = acc + abs(true_val - val)
          if(abs(true_val-val).gt.1.0d-14) then
             print *, 'i,j = ',i,j, 'eta1,eta2 = ', eta1, eta2, 'true, val = ',&
                  true_val, val
             !print *, spline%coeffs(npx1,j-1:j+1)
          end if
       end do
    end do
    ave_err = acc/real(npx1*npx2,f64) 
    print *, 'test_2d_spline_prdc_hrmt results. Average error = ', ave_err
    if( ave_err .le. 1.0d-14 ) then
       test_passed = .true.
    else
       test_passed = .false.
    end if
    SLL_DEALLOCATE_ARRAY(data,ierr)
    SLL_DEALLOCATE_ARRAY(eta2_min_slopes,ierr)
    SLL_DEALLOCATE_ARRAY(eta2_max_slopes,ierr)
  end subroutine test_2d_spline_prdc_hrmt_no_slopes


  subroutine test_2d_spline_hrmt_hrmt( &
    transform_func, &
    eta1_min_slope_func, &
    eta1_max_slope_func, &
    eta2_min_slope_func, &
    eta2_max_slope_func, &
    test_passed )
 
    procedure(fxy) :: transform_func
    procedure(fxy) :: eta1_min_slope_func
    procedure(fxy) :: eta1_max_slope_func
    procedure(fxy) :: eta2_min_slope_func
    procedure(fxy) :: eta2_max_slope_func
    logical, intent(out) :: test_passed

    sll_int32 :: i, j, ierr
    sll_real64, dimension(:,:), allocatable :: data
    sll_real64, dimension(:), allocatable   :: eta1_min_slopes
    sll_real64, dimension(:), allocatable   :: eta1_max_slopes
    sll_real64, dimension(:), allocatable   :: eta2_min_slopes
    sll_real64, dimension(:), allocatable   :: eta2_max_slopes
    sll_real64 :: h1, h2, eta1, eta2, acc, val, true_val, ave_err
    type(sll_t_cubic_spline_2d), pointer :: spline

    h1 = 1.0_f64/(npx1-1)
    h2 = 1.0_f64/(npx2-1)
    ! allocate and fill out the data
    SLL_ALLOCATE(data(npx1,npx2),ierr)
    do j=0,npx2-1
       eta2 = real(j,f64)*h2
       do i=0,npx1-1
          eta1 = real(i,f64)*h1
          data(i+1,j+1) = transform_func(eta1,eta2)
       end do
    end do

    ! allocate and fill out the bc data
    SLL_ALLOCATE(eta1_min_slopes(npx1),ierr)
    SLL_ALLOCATE(eta1_max_slopes(npx1),ierr)
    SLL_ALLOCATE(eta2_min_slopes(npx1),ierr)
    SLL_ALLOCATE(eta2_max_slopes(npx1),ierr)

    do i=0,npx1-1
       eta1 = real(i,f64)*h1
       eta2 = 0.0_f64
       eta2_min_slopes(i+1) = eta2_min_slope_func(eta1,eta2)
       eta2_max_slopes(i+1) = eta2_max_slope_func(eta1,eta2)
    end do

    do j=0,npx2-1
       eta1 = 0.0_f64
       eta2 = real(j,f64)*h2
       eta1_min_slopes(j+1) = eta1_min_slope_func(eta1,eta2)
       eta1_max_slopes(j+1) = eta1_max_slope_func(eta1,eta2)
    end do

    spline =>sll_f_new_cubic_spline_2d( &
         npx1, &
         npx2, &
         0.0_f64, &
         1.0_f64, &
         0.0_f64, &
         1.0_f64, &
         sll_p_hermite, &
         sll_p_hermite, & 
         x1_min_slopes=eta1_min_slopes, &
         x1_max_slopes=eta1_max_slopes, &
         x2_min_slopes=eta2_min_slopes, &
         x2_max_slopes=eta2_max_slopes )

    call sll_s_compute_cubic_spline_2d( data, spline )

    ! compare results
    acc = 0.0_f64
    do j=0,npx2-1
       eta2 = real(j,f64)*h2
       do i=0,npx1-1
          eta1     = real(i,f64)*h1
          val      = sll_f_interpolate_value_2d( eta1, eta2, spline )
          true_val = transform_func( eta1, eta2 )
          acc      = acc + abs(true_val - val)
          if(abs(true_val-val).gt.1.0d-14) then
             print *, 'i,j = ',i,j, 'eta1,eta2 = ', eta1, eta2, 'true, val = ',&
                  true_val, val
             !print *, spline%coeffs(npx1,j-1:j+1)
          end if
       end do
    end do
    ave_err = acc/real(npx1*npx2,f64) 
    print *, 'test_2d_spline_hrmt_hrmt results. Average error = ', ave_err
    if( ave_err .le. 1.0d-14 ) then
       test_passed = .true.
    else
       test_passed = .false.
    end if
    SLL_DEALLOCATE_ARRAY(data,ierr)
    SLL_DEALLOCATE_ARRAY(eta2_min_slopes,ierr)
    SLL_DEALLOCATE_ARRAY(eta2_max_slopes,ierr)
  end subroutine test_2d_spline_hrmt_hrmt

  subroutine test_2d_spline_hrmt_hrmt_no_slopes( transform_func, test_passed )
    procedure(fxy) :: transform_func
    logical, intent(out) :: test_passed

    sll_int32 :: i, j, im, jm, ierr
    sll_real64, dimension(:,:), allocatable :: data
    sll_real64 :: h1, h2, eta1, eta2, acc, val, true_val, ave_err, max_err
    type(sll_t_cubic_spline_2d), pointer :: spline

    h1 = 1.0_f64/(npx1-1)
    h2 = 1.0_f64/(npx2-1)
    ! allocate and fill out the data
    SLL_ALLOCATE(data(npx1,npx2),ierr)
    do j=0,npx2-1
       eta2 = real(j,f64)*h2
       do i=0,npx1-1
          eta1 = real(i,f64)*h1
          data(i+1,j+1) = transform_func(eta1,eta2)
       end do
    end do

    spline =>sll_f_new_cubic_spline_2d( &
         npx1, &
         npx2, &
         0.0_f64, &
         1.0_f64, &
         0.0_f64, &
         1.0_f64, &
         sll_p_hermite, &
         sll_p_hermite )

    call sll_s_compute_cubic_spline_2d( data, spline )

    ! compare results
    max_err = 0.0_f64
    acc     = 0.0_f64
    do j=0,npx2-1
       eta2 = real(j,f64)*h2
       do i=0,npx1-1
          eta1     = real(i,f64)*h1
          val      = sll_f_interpolate_value_2d( eta1, eta2, spline )
          true_val = transform_func( eta1, eta2 )
          if( abs(val - true_val) > max_err ) then
             max_err = abs(val-true_val)
             im = i
             jm = j
          end if
          acc      = acc + abs(true_val - val)
          if(abs(true_val-val).gt.1.0d-14) then
             print *, 'i,j = ',i,j, 'eta1,eta2 = ', eta1, eta2, 'true, val = ',&
                  true_val, val
             !print *, spline%coeffs(npx1,j-1:j+1)
          end if
       end do
    end do
    ave_err = acc/real(npx1*npx2,f64) 
    print *, 'test_2d_spline_hrmt_hrmt results. Average error = ', ave_err
    if( ave_err .le. 1.0d-14 ) then
       test_passed = .true.
    else
       test_passed = .false.
    end if
    print *, 'hrmt_hrmt maximum error found = ', max_err, 'at i,j = ', &
         im+1, jm+1
    SLL_DEALLOCATE_ARRAY(data,ierr)
  end subroutine test_2d_spline_hrmt_hrmt_no_slopes



  function line(x)
    sll_real64 :: line
    sll_real64, intent(in) :: x
    line = x
  end function line

  function plane( x, y )
    sll_real64 :: plane
    sll_real64, intent(in) :: x, y
    plane = 2.0*x + 1.0+0._f64*y
  end function plane

  function plane_deriv( x, y )
    sll_real64 :: plane_deriv
    sll_real64, intent(in) :: x, y
    plane_deriv = 2.0+ 0._f64*x+0._f64*y
  end function plane_deriv

  function plane2( x, y)
    sll_real64 :: plane2
    sll_real64, intent(in) :: x, y
    plane2 = 2.0*y + 1.0+ 0._f64*x
  end function plane2

  function plane2_deriv( x, y )
    sll_real64 :: plane2_deriv
    sll_real64, intent(in) :: x, y
    plane2_deriv = 2.0+ 0._f64*x+0._f64*y
  end function plane2_deriv

  function plane3( x, y)
    sll_real64 :: plane3
    sll_real64, intent(in) :: x, y
    plane3 = 1.0 + 2.0*x + 2.0*y 
  end function plane3

  function plane3_deriv_x( x, y )
    sll_real64 :: plane3_deriv_x
    sll_real64, intent(in) :: x, y
    plane3_deriv_x = 2.0+ 0._f64*x+0._f64*y
  end function plane3_deriv_x

  function plane3_deriv_y( x, y )
    sll_real64 :: plane3_deriv_y
    sll_real64, intent(in) :: x, y
    plane3_deriv_y = 2.0+ 0._f64*x+0._f64*y
  end function plane3_deriv_y

  function polar_x( eta1, eta2 )
    sll_real64 :: polar_x
    sll_real64, intent(in) :: eta1, eta2
    polar_x = (r1 + eta1*(r2-r1))*cos(2*sll_p_pi*eta2)
  end function polar_x

  function polar_y( eta1, eta2 )
    sll_real64 :: polar_y
    sll_real64, intent(in) :: eta1, eta2
    polar_y = (r1 + eta1*(r2-r1))*sin(2*sll_p_pi*eta2)
  end function polar_y

  function deriv1_polar_x( eta1, eta2 )
    sll_real64 :: deriv1_polar_x
    sll_real64, intent(in) :: eta1, eta2
    deriv1_polar_x = (r2-r1)*cos(2.0_f64*sll_p_pi*eta2)+0._f64*eta1
  end function deriv1_polar_x

  function deriv2_polar_x( eta1, eta2 )
    sll_real64 :: deriv2_polar_x
    sll_real64, intent(in) :: eta1, eta2
    deriv2_polar_x = -(r1+eta1*(r2-r1))*sin(2.0_f64*sll_p_pi*eta2)*2.0_f64*sll_p_pi
  end function deriv2_polar_x

  function deriv1_polar_y( eta1, eta2 )
    sll_real64 :: deriv1_polar_y
    sll_real64, intent(in) :: eta1, eta2
    deriv1_polar_y = (r2-r1)*sin(2.0_f64*sll_p_pi*eta2)+0._f64*eta1
  end function deriv1_polar_y

  function deriv2_polar_y( eta1, eta2 )
    sll_real64 :: deriv2_polar_y
    sll_real64, intent(in) :: eta1, eta2
    deriv2_polar_y = (r1+eta1*(r2-r1))*cos(2.0_f64*sll_p_pi*eta2)*2.0_f64*sll_p_pi
  end function deriv2_polar_y

  function mycos(x)
    intrinsic :: cos
    sll_real64 :: mycos
    sll_real64, intent(in) :: x
    mycos = cos(x)
  end function mycos

  function dmycos(x)
    intrinsic :: sin
    sll_real64 :: dmycos
    sll_real64, intent(in) :: x
    dmycos = -sin(x)
  end function dmycos

  function coscos(x,y)
    sll_real64 :: coscos
    sll_real64, intent(in) :: x, y
    coscos = cos(x)*cos(y)
  end function coscos

  function msincos(x,y)
    sll_real64 :: msincos
    sll_real64, intent(in) :: x, y
    msincos = -cos(y)*sin(x)
  end function msincos

  function mcossin(x,y)
    sll_real64 :: mcossin
    sll_real64, intent(in) :: x, y
    mcossin = -cos(x)*sin(y)
  end function mcossin

  function sinsin(x,y)
    sll_real64 :: sinsin
    sll_real64, intent(in) :: x, y
    sinsin = sin(2.0_f64*sll_p_pi*x)*sin(2.0_f64*sll_p_pi*y)
  end function sinsin

end module test_processes_module
