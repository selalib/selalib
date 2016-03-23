program test_cubic_splines
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"
  use sll_m_boundary_condition_descriptors, only: &
    sll_p_hermite, &
    sll_p_periodic

  use sll_m_constants, only : &
       sll_p_twopi, &
       sll_p_pi
 use sll_m_cubic_splines, only: &
    sll_s_compute_cubic_spline_1d, &
    sll_s_compute_cubic_spline_2d, &
    sll_f_interpolate_derivative, &
    sll_s_interpolate_from_interpolant_array, &
    sll_f_interpolate_from_interpolant_value, &
    sll_f_interpolate_value_2d, &
    sll_f_interpolate_x1_derivative_2d, &
    sll_f_interpolate_x2_derivative_2d, &
    sll_f_new_cubic_spline_1d, &
    sll_f_new_cubic_spline_2d, &
    sll_t_cubic_spline_1d, &
    sll_t_cubic_spline_2d, &
    sll_o_delete

  use test_func_module, only: &
    f, &
    fprime

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  logical :: passed

  passed = .true.

  ! Test 1d cubic splines with periodic boundary conditions
  call test_cubic_spline_1d( 0, passed )
  ! Test 1d cubic splines with Hermite boundary condtions
  call test_cubic_spline_1d( 1, passed )

  ! Test 2d cubic splines with periodic-periodic boundary conditions
  call test_cubic_spline_2d( 0, passed )
  ! Test 2d cubic splines with periodic-Hermite boundary condtions
  call test_cubic_spline_2d( 1, passed )
  ! Test 2d cubic splines with Hermite-periodic boundary conditions
  call test_cubic_spline_2d( 2, passed )
  ! Test 2d cubic splines with Hermite-Hermite boundary condtions
  call test_cubic_spline_2d( 3, passed )

  if (passed .eqv. .true. ) then
     print*, 'PASSED.'
  else
     print*, 'FAILED.'
  end if

contains

  subroutine test_cubic_spline_1d( bc, passed )
    integer(kind=4), intent( in    ) :: bc ! specify boundary condtions: 0 for periodic, 1 for Hermite
    logical,         intent( inout ) :: passed

    sll_int32, parameter :: np = 32

    sll_int32 :: i
    sll_real64 :: xmin, xmax
    sll_real64 :: delta_x
    
    sll_real64 :: x(np+1)
    sll_real64 :: data(np+1)
    sll_real64 :: deriv(np+1)
    sll_real64 :: data_interp(np)
    sll_real64 :: deriv_interp(np)
    sll_real64 :: interp_ngrid(2)
    sll_real64 :: ref_ngrid(2)
    sll_real64 :: x_ngrid
    type(sll_t_cubic_spline_1d), pointer :: sp1
    
    xmin = 0.0_f64
    xmax = sll_p_twopi
    
    delta_x = (xmax-xmin)/real(np,f64)
    
    do i=1,np+1
       x(i) = real(i-1, f64)* delta_x
       data(i) = exp(sin(x(i)))
       deriv(i) = cos(x(i)) * exp(sin(x(i)))
    end do
    
    if ( bc == 0) then! periodic boundary conditions       
       print*, 'Cubic spline 1d, periodic boundary conditions:'
       sp1 =>  sll_f_new_cubic_spline_1d( np+1, xmin, xmax, sll_p_periodic )
    elseif( bc == 1) then ! Hermite boundary conditions
       print*, 'Cubic spline 1d, Hermite boundary conditions:'
       sp1 =>  sll_f_new_cubic_spline_1d( &
            np+1, &
            xmin, &
            xmax, &
            sll_p_hermite, &
            deriv(1), &
            deriv(np+1) )
    end if
    call sll_s_compute_cubic_spline_1d( data, sp1 )
    
    do i=1,np
       data_interp(i) = sll_f_interpolate_from_interpolant_value(x(i), sp1)
       deriv_interp(i) = sll_f_interpolate_derivative(x(i), sp1)
    end do
    
    x_ngrid = (real(np/2,f64)-0.5_f64)*delta_x
    ref_ngrid(1) = exp(sin(x_ngrid))
    interp_ngrid(1) = sll_f_interpolate_from_interpolant_value(x_ngrid, sp1)
    ref_ngrid(2) = cos(x_ngrid)*exp(sin(x_ngrid))
    interp_ngrid(2) = sll_f_interpolate_derivative(x_ngrid, sp1)
    
    print*, 'Interpolation error for value at', x_ngrid, ':', ref_ngrid(1)-interp_ngrid(1)
    print*, 'Interpolation error for derivative at', x_ngrid, ':',ref_ngrid(2)-interp_ngrid(2)
    print*, 'Maximum error of values at grid points:', maxval(abs(data(1:np)-data_interp))
    print*, 'Maximum error of derivative at grid points:', maxval(abs(deriv(1:np)-deriv_interp))
    
    
    if ( maxval(abs(data(1:np)-data_interp))> 1d-14 ) then
       print*, 'Maximum error at grid points too large.'
       passed = .false.
    end if
    if ( maxval(abs(deriv(1:np)-deriv_interp))> 3d-4 ) then
       print*, 'Maximum error of derivative at grid points too large.'
       passed = .false.
    end if
    if (abs(ref_ngrid(1)-interp_ngrid(1)) > 2d-5 ) then
       print*, 'Error at non-grid point too large.'
       passed = .false.
    end if
    if (abs(ref_ngrid(2)-interp_ngrid(2)) > 3d-5 ) then
       print*, 'Error of derivative at non-grid point too large.'
       passed = .false.
    end if
    
    call sll_o_delete(sp1)

  end subroutine test_cubic_spline_1d


  subroutine test_cubic_spline_2d( bc, passed )
    integer(kind=4), intent( in    ) :: bc ! specify boundary condtions: 0 for periodic, 1 for Hermite
    logical,         intent( inout ) :: passed

    sll_int32, parameter :: np1 = 28
    sll_int32, parameter :: np2 = 32

    sll_int32 :: i, j
    sll_real64 :: xmin, xmax
    sll_real64 :: delta_x
    sll_real64 :: ymin, ymax
    sll_real64 :: delta_y
    
    sll_real64 :: x(np1+1), y(np2+1)
    sll_real64 :: data(np1+1, np2+1)
    sll_real64 :: deriv_x(np1+1, np2+1), deriv_y(np1+1, np2+1)
    sll_real64 :: data_interp(np1, np2)
    sll_real64 :: deriv_x_interp(np1, np2), deriv_y_interp(np1, np2)
    sll_real64 :: interp_ngrid(3)
    sll_real64 :: ref_ngrid(3)
    sll_real64 :: x_ngrid, y_ngrid
    type(sll_t_cubic_spline_2d), pointer :: sp2
    
    print*, 'Cubic splines 2d:'
    xmin = 0.0_f64
    xmax = sll_p_twopi
    
    delta_x = (xmax-xmin)/real(np1,f64)

    do i=1,np1+1
       x(i) = real(i-1, f64)* delta_x
    end do
    ymin = 0.0_f64
    ymax = sll_p_twopi
    
    delta_y = (ymax-ymin)/real(np2,f64)

    do j=1,np2+1
       y(j) = real(j-1, f64)* delta_y
    end do

    do i=1, np1+1
       do j=1, np2+1
          data(i,j) = exp( cos(x(i))* sin(y(j)) )
          deriv_x(i,j) = -data(i,j) * sin( x(i) ) * sin( y(j) )
          deriv_y(i,j) = data(i,j) * cos( x(i) ) * cos( y(j) )
       end do
    end do
    
    if (bc == 0) then ! periodic-periodic
       print*, 'Cubic spline 2d, periodic-periodic:'
       sp2 => sll_f_new_cubic_spline_2d( np1+1, np2+1, &
            xmin, xmax, &
            ymin, ymax, &
            sll_p_periodic, sll_p_periodic)
    elseif (bc == 1) then ! periodic-Hermite
       print*, 'Cubic spline 2d, periodic-periodic:'
       sp2 => sll_f_new_cubic_spline_2d( np1+1, np2+1, &
            xmin, xmax, &
            ymin, ymax, &
            sll_p_periodic, sll_p_hermite, &
            x2_min_slopes = deriv_y(:,1), &
            x2_max_slopes = deriv_y(:,1) )
    elseif(bc == 2) then ! Hermite-periodic
       print*, 'Cubic spline 2d, periodic-periodic:'
       sp2 => sll_f_new_cubic_spline_2d( np1+1, np2+1, &
            xmin, xmax, &
            ymin, ymax, &
            sll_p_hermite, sll_p_periodic, &
            x1_min_slopes = deriv_x(1,:), &
            x1_max_slopes = deriv_x(1,:) )
    elseif( bc == 3 ) then ! Hermite-Hermite
       sp2 => sll_f_new_cubic_spline_2d( np1+1, np2+1, &
            xmin, xmax, &
            ymin, ymax, &
            sll_p_hermite, sll_p_hermite, &
            x1_min_slopes = deriv_x(1,:), &
            x1_max_slopes = deriv_x(1,:), &
            x2_min_slopes = deriv_y(:,1), &
            x2_max_slopes = deriv_y(:,1) )
    else
       print*, 'Boundary type not implemented.'
       passed = .false.
    end if

    call sll_s_compute_cubic_spline_2d( data, sp2 )

    do i=1,np1
       do j=1, np2
          data_interp(i,j) = sll_f_interpolate_value_2d(x(i), y(j), sp2)
          deriv_x_interp(i,j) = sll_f_interpolate_x1_derivative_2d(x(i), y(j), sp2)
          deriv_y_interp(i,j) = sll_f_interpolate_x2_derivative_2d(x(i), y(j), sp2)
       end do
    end do
    
    x_ngrid = (real(np1/2,f64)-0.5_f64)*delta_x
    y_ngrid = (real(np2/2,f64)+0.5_f64)*delta_y
    ref_ngrid(1) =  exp( cos(x_ngrid)* sin(y_ngrid) )
    interp_ngrid(1) = sll_f_interpolate_value_2d(x_ngrid, y_ngrid, sp2)
    ref_ngrid(2) = -ref_ngrid(1) * sin( x_ngrid ) * sin( y_ngrid )
    interp_ngrid(2) = sll_f_interpolate_x1_derivative_2d(x_ngrid, y_ngrid, sp2)
    ref_ngrid(3) = ref_ngrid(1) * cos( x_ngrid ) * cos( y_ngrid )
    interp_ngrid(3) = sll_f_interpolate_x2_derivative_2d(x_ngrid, y_ngrid, sp2)
    
    print*, 'Interpolation error for value at (', x_ngrid, ',', y_ngrid, '):', ref_ngrid(1)-interp_ngrid(1)
    print*, 'Interpolation error for x1 derivative at (', x_ngrid, ',', y_ngrid, '):' ,ref_ngrid(2)-interp_ngrid(2)
    print*, 'Interpolation error for x2 derivative at (', x_ngrid, ',', y_ngrid, '):' ,ref_ngrid(3)-interp_ngrid(3)
    print*, 'Maximum error of values at grid points:', maxval(abs(data(1:np1, 1:np2)-data_interp))
    print*, 'Maximum error of x1 derivative at grid points:', maxval(abs(deriv_x(1:np1, 1:np2)-deriv_x_interp))
    print*, 'Maximum error of x2 derivative at grid points:', maxval(abs(deriv_y(1:np1, 1:np2)-deriv_y_interp))

    if ( maxval(abs(data(1:np1, 1:np2)-data_interp))> 1d-14 ) then
       print*, 'Maximum error at grid points too large.'
       passed = .false.
    end if
    if ( maxval(abs(deriv_x(1:np1, 1:np2)-deriv_x_interp))> 4d-4 ) then
       print*, 'Maximum error of x1 derivative at grid points too large.'
       passed = .false.
    end if
    if ( maxval(abs(deriv_y(1:np1, 1:np2)-deriv_y_interp))> 3d-4 ) then
       print*, 'Maximum error of x1 derivative at grid points too large.'
       passed = .false.
    end if
    if (abs(ref_ngrid(1)-interp_ngrid(1)) > 2d-5 ) then
       print*, 'Error at non-grid point too large.'
       passed = .false.
    end if
    if (abs(ref_ngrid(2)-interp_ngrid(2)) > 3d-6 ) then
       print*, 'Error of x1 derivative at non-grid point too large.'
       passed = .false.
    end if    
    if (abs(ref_ngrid(3)-interp_ngrid(3)) > 4d-5 ) then
       print*, 'Error of x2 derivative at non-grid point too large.'
       passed = .false.
    end if
    

  end subroutine test_cubic_spline_2d


end program test_cubic_splines
