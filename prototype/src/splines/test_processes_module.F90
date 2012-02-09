
!****************************************************
!
! Selalib      
! Module: test_processes_module.F90
!
!> @brief 
!> Module of test processes for sll_splines unit test
!   
!> @authors                    
!> Aliou DIOUF (aliou.l.diouf@inria.fr)
!                                  
!****************************************************

module test_processes_module
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"
  use sll_splines
  use numeric_constants
  use util_constants
    use test_func_module
  implicit none
contains


  subroutine test_process_1d(i_test, ok)
    implicit none

    sll_int32                             :: i_test, ok, i, err 
    type(sll_spline_1d), pointer          :: sp1
    type(sll_spline_1d), pointer          :: sp2
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
    
    SLL_ALLOCATE(data(NP), err)
    SLL_ALLOCATE(deriv(NP), err)
    SLL_ALLOCATE(out(NP), err)
    SLL_ALLOCATE(coordinates(NP), err)
    
    print *, 'initialize data and coordinates array'
    do i=1,NP
       phase          = real(i-1,f64)*(XMAX-XMIN)/real(NP-1,f64) + XMIN
       data(i)        = f(phase, i_test)
       deriv(i)       = fprime(phase, i_test)
       coordinates(i) = phase
    end do
     
    print *, 'proceed to allocate the spline...'
    sp1 =>  new_spline_1D( NP, XMIN, XMAX, PERIODIC_SPLINE )
    call compute_spline_1D( data, PERIODIC_SPLINE, sp1 )
    sp2 =>  new_spline_1D( NP, XMIN, XMAX, HERMITE_SPLINE, fprime(XMIN, i_test), fprime(XMAX, i_test) )
    call compute_spline_1D( data, HERMITE_SPLINE, sp2 )
    
    print *, 'Contents of the spline 1:'
    print *, sp1%xmin
    print *, sp1%xmax
    print *, sp1%delta
    print *, sp1%rdelta
    print *, sp1%bc_type
#if PRINT_SPLINE_COEFFS
    print *, sp1%coeffs(:)
#endif
    print *, 'Contents of the spline 2:'
    print *, sp2%xmin
    print *, sp2%xmax
    print *, sp2%delta
    print *, sp2%rdelta
    print *, sp2%bc_type
#if PRINT_SPLINE_COEFFS
    print *, sp2%coeffs(:)
#endif
     
    print *, 'cumulative errors: '
    print *, 'periodic case, NP-1 points: '
    print *, 'interpolating individual values from 1 to NP-1:'
    do i=1, NP-1
       x = real(i-1,f64)*(XMAX-XMIN)/real(NP-1,f64)+XMIN
       accumulator1 = accumulator1 + abs(data(i) - interpolate_value(x, sp1))
    end do
    print *, 'checking periodicity:'
    print *, 'difference between values at points 1 and NP: ', &
         abs(data(1) - interpolate_value(XMAX,sp1))
    print *, 'interpolating the whole array:'
    call interpolate_array_values(coordinates, out, NP-1, sp1)
    do i=1, NP-1
       accumulator3 = accumulator3 + abs(data(i) - out(i))
    end do
    
    print *, 'hermite case, NP points: '
    do i=1, NP
       x = real(i-1,f64)*(XMAX-XMIN)/real(NP-1,f64) + XMIN
       accumulator2 = accumulator2 + abs(data(i) - interpolate_value(x, sp2))
    end do
    call interpolate_array_values(coordinates, out, NP, sp2)
    do i=1, NP
       accumulator4 = accumulator4 + abs(data(i) - out(i))
    end do
     
    print *, '----------------------------------------------------'
    print *, 'RESULTS: '
    print *, 'Periodic case: '
    print *, 'average error at the nodes (single values) = '
    print *, accumulator1/real(NP,f64)
    print *, 'average error at the nodes (whole array) = '
    print *, accumulator3/real(NP,f64)
    write (*,'(a,f8.5)')   'original data(0)    = ', data(1)
    write (*,'(a,f20.15)') &
         'interpolated        = ', interpolate_value( XMIN,sp1)
    
    write (*,'(a,f20.15)')   'original data((NP-1)/4) = ', data((NP-1)/4)
    write (*,'(a,f20.15)') &
         'interpolated        = ', interpolate_value( (XMAX-XMIN)/4.0+XMIN,sp1)
     
     
#if TEST_INTEGRATION
    print *, 'integrating the periodic spline...'
    print *, gauss_legendre_integrate_1D( interpolate_value, sp1, XMIN, XMAX,4)
#endif
     
#if PRINT_SPLINE_COEFFS
    print *, 'spline coefficients: '
    print *, sp1%coeffs(:)
#endif
     
    if ( (accumulator1/real(NP,f64) >= 1.0e-15) .or. &
          (accumulator3/real(NP,f64) >= 1.0e-15) ) then 
       ok = 0
       print*, 'i_test =', i_test
       print *, 'splines unit test stopped by periodic spline1d test failure'
       stop
    end if
    
    print *, '**************************** '
    print *, 'Hermite case: '
    print *, 'average error at the nodes = '
    print *, accumulator2/real(NP,f64)
    write (*,'(a,f8.5)')   'original data(0)    = ', data(1)
    write (*,'(a,f20.15)') &
         'interpolated        = ', interpolate_value( 0.0_f64,sp2)
    write (*,'(a,f20.15)')   'original data((NP-1)/4) = ', data((NP-1)/4+1)
    write (*,'(a,f20.15)') &
         'interpolated        = ', interpolate_value( (XMAX-XMIN)/4.0,sp2)
     print *, 'spline coefficients: '
#if PRINT_SPLINE_COEFFS
    print *, sp2%coeffs(:)
#endif
     
    if ( (accumulator2/real(NP,f64) >= 1.0e-15) .or. &
         (accumulator4/real(NP,f64) >= 1.0e-15) ) then
       ok = 0
       print*, 'i_test =', i_test
       print *, 'splines unit test stopped by Hermite slpine1d test failure'
       stop
    end if
    if (i_test<=12) then 
       print *, '---------------------------------------------'
       print *, 'DERIVATIVES TEST'
       do i=1, NP
          phase = real(i-1,f64)*(XMAX-XMIN)/real(NP-1,f64)+XMIN
          accumulator5 = accumulator5 + abs(deriv(i) - &
                         interpolate_derivative(phase, sp1))
          accumulator6 = accumulator6 + abs(deriv(i) - &
                         interpolate_derivative(phase, sp2))
       end do
    
       if ( &
            (accumulator5/real(NP,f64) > 1.0e-6) .or. &
            (accumulator6/real(NP,f64) > 1.0e-6) ) then
          ok = 0
          print *, 'average error at the nodes (single values, periodic) = '
          print *, accumulator5/real(NP,f64)
          print *, 'average error at the nodes (single values, hermite) = '
          print *, accumulator6/real(NP,f64)
          print *, 'i_test =', i_test
          print *, 'splines unit test stopped by DERIVATIVES TEST failure'
          stop
       end if
    endif
     
    call delete_spline_1D(sp1)
    call delete_spline_1D(sp2)
    
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
    type(sll_spline_2d), pointer            :: sp2d_1
    type(sll_spline_2d), pointer            :: sp2d_2
    type(sll_spline_2d), pointer            :: sp2d_3
    type(sll_spline_2d), pointer            :: sp2d_4
    sll_real64                              :: phase_x1, phase_x2
    sll_real64                              :: acc_2D, x1, x2
    sll_real64, allocatable, dimension(:,:) :: data_2d
    sll_real64, allocatable, dimension(:)   :: coordinates_i
    sll_real64, allocatable, dimension(:)   :: coordinates_j

    SLL_ALLOCATE(data_2d(NPX1, NPX2), err)        
    SLL_ALLOCATE(coordinates_i(NPX1),err)
    SLL_ALLOCATE(coordinates_j(NPX2),err)
        
    do i=1, NPX1
       coordinates_i(i) = real(i-1,f64)*(X1MAX-X1MIN)/real(NPX1-1,f64)+X1MIN
    end do
        
    do j=1, NPX2
       coordinates_j(j) = real(j-1,f64)*(X2MAX-X2MIN)/real(NPX2-1,f64)+X2MIN
    end do
        
    do j=1,NPX2
       phase_x2 = coordinates_j(j)
       do i=1,NPX1
          phase_x1 = coordinates_i(i)
          data_2d(i,j) = f(phase_x1, i_test)*f(phase_x2, j_test)
       end do
    end do
        
    ! Test the periodic-periodic spline2d 
    sp2d_1 => new_spline_2D( NPX1, NPX2, &
         X1MIN, X1MAX, &
         X2MIN, X2MAX, &
         PERIODIC_SPLINE, PERIODIC_SPLINE )
    call compute_spline_2D_prdc_prdc( data_2d, sp2d_1 )
    acc_2D = 0.0
    do j=1, NPX2
       do i=1, NPX1
          x1 = coordinates_i(i)
          x2 = coordinates_j(j)
          acc_2D = acc_2D + &
               abs(data_2d(i,j) - interpolate_value_2D(x1,x2,sp2d_1))
       end do
    end do
    if (acc_2D/(NPX1*NPX2)>=1.e-14) then
       ok = 0
       print*, '(i_test, j_test) = (', i_test, ',', j_test,')' 
       print *, 'Average cumulative error, spline2d, periodic-periodic: ', &
         acc_2D/(NPX1*NPX2)
       print *, 'splines unit test stopped by periodic-periodic spline2d test failure'
       stop
    endif
    ! Test the Hermite-periodic spline2d
    sp2d_2 => new_spline_2D( NPX1, NPX2, &
         X1MIN, X1MAX, &
         X2MIN, X2MAX, &
         HERMITE_SPLINE, PERIODIC_SPLINE, &
         x1_min_slope = fprime(X1MIN, i_test), &
         x1_max_slope = fprime(X1MAX, i_test) )
    call compute_spline_2D_hrmt_prdc( data_2d, sp2d_2 )
    acc_2D = 0.0
    do j=1, NPX2
       do i=1, NPX1
          x1 = coordinates_i(i)
          x2 = coordinates_j(j)
          acc_2D = acc_2D + &
               abs(data_2d(i,j) - interpolate_value_2D(x1,x2,sp2d_2))
       end do
    end do
    if (acc_2D/(NPX1*NPX2)>=1.e-14) then
       ok = 0
       print*, '(i_test, j_test) = (', i_test, ',', j_test,')' 
       print *, 'Average cumulative error, spline2d, Hermite-periodic: ', &
         acc_2D/(NPX1*NPX2)
       print *, 'splines unit test stopped by Hermite-periodic spline2d test failure'
       stop
    endif
    ! Test the periodic-Hermite spline2d
    sp2d_3 => new_spline_2D( NPX1, NPX2, &
         X1MIN, X1MAX, &
         X2MIN, X2MAX, &
         PERIODIC_SPLINE, HERMITE_SPLINE, &
         x2_min_slope = fprime(X2MIN, j_test), &
         x2_max_slope = fprime(X2MAX, j_test) )
    call compute_spline_2D_prdc_hrmt( data_2d, sp2d_3 )
    acc_2D = 0.0
    do j=1, NPX2
       do i=1, NPX1
          x1 = coordinates_i(i)
          x2 = coordinates_j(j)
          acc_2D = acc_2D + &
               abs(data_2d(i,j) - interpolate_value_2D(x1,x2,sp2d_3))
       end do
    end do
    if (acc_2D/(NPX1*NPX2)>=1.e-14) then
       ok = 0
       print*, '(i_test, j_test) = (', i_test, ',', j_test,')' 
       print *, 'Average cumulative error, spline2d, periodic-Hermite: ', &
         acc_2D/(NPX1*NPX2)
       print *, 'splines unit test stopped by periodic-Hermite spline2d test failure'
       stop
    endif
    ! Test the Hermite-Hermite spline2d
    sp2d_4 => new_spline_2D( NPX1, NPX2, &
         X1MIN, X1MAX, &
         X2MIN, X2MAX, &
         HERMITE_SPLINE, HERMITE_SPLINE,&
         fprime(X1MIN, i_test), fprime(X1MAX, i_test), &
         fprime(X2MIN, j_test), fprime(X2MAX, j_test) )   
    call compute_spline_2D_hrmt_hrmt( data_2d, sp2d_4 )
    acc_2D = 0.0
    do j=1, NPX2
       do i=1, NPX1
          x1 = coordinates_i(i)
          x2 = coordinates_j(j)
          acc_2D = acc_2D + &
               abs(data_2d(i,j) - interpolate_value_2D(x1,x2,sp2d_4))
       end do
    end do
    if (acc_2D/(NPX1*NPX2)>=1.e-14) then
       ok = 0
       print*, '(i_test, j_test) = (', i_test, ',', j_test,')' 
       print *, 'Average cumulative error, spline2d, Hermite-Hermite: ', &
         acc_2D/(NPX1*NPX2)
       print *, 'splines unit test stopped by Hermite-Hermite spline2d test failure'
       stop
    endif
     
    call delete(sp2d_1)
    call delete(sp2d_2)
    call delete(sp2d_3)
    call delete(sp2d_4)
    
    SLL_DEALLOCATE_ARRAY(data_2d, err)
    SLL_DEALLOCATE_ARRAY(coordinates_i,err)
    SLL_DEALLOCATE_ARRAY(coordinates_j,err)

  end subroutine test_process_2d

end module test_processes_module
