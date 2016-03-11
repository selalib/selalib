!PN This is a general test for bsplines 1d and 2d
!PN This test is used to test performance and precision
program test_bsplines_1d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_boundary_condition_descriptors, only: &
       sll_p_hermite, &
       sll_p_greville, &
       sll_p_periodic

  use sll_m_bspline_interpolation, only: &
    sll_s_compute_bspline_1d, &
 !   sll_o_compute_bspline_2d, &
    sll_s_interpolate_array_derivatives_1d, &
    sll_s_interpolate_array_values_1d, &
 !   sll_s_interpolate_array_values_2d, &
    sll_f_interpolate_derivative_1d, &
    sll_f_interpolate_value_1d, &
    !   sll_f_interpolate_value_2d, &
    sll_s_bspline_interpolation_1d_init, &
    sll_s_bspline_interpolation_1d_free, &
    sll_s_bspline_interpolation_2d_init, &
    sll_t_bspline_interpolation_1d, &
    sll_t_bspline_interpolation_2d 


  use sll_m_constants, only: &
    sll_p_pi

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sll_real64                    :: err1
sll_real64                    :: err2
sll_real64                    :: err3
sll_real64                    :: err4
sll_real64                    :: err5
sll_real64                    :: err6
sll_int32                     :: i
sll_int32                     :: j
sll_int32,  parameter         :: nstep = 1
sll_real64,  parameter        :: tol = 1.0d-3    ! tolerance for tests
!sll_int32,  parameter         :: m = 2
sll_real64                    :: t0, t1, t2, t3, t4, t5
sll_int32                     :: ierr
sll_int32                     :: deg
logical                                :: passed_test
passed_test = .true.
  
print*,'***************************************************************'
print*,'*** 1D PERIODIC ***'
print*,'***************************************************************'
do deg=3,6
!   call test_process_1d(sll_p_periodic,deg, passed_test)
end do
print*,'***************************************************************'
print*,'*** 1D GREVILLE ***'
print*,'***************************************************************'
do deg=3,3
!   call test_process_1d(sll_p_greville,deg, passed_test)
end do
print*,'***************************************************************'
print*,'*** 1D HERMITE ***'
print*,'***************************************************************'
do deg=3,3
   call test_process_1d(sll_p_hermite,deg,passed_test)
end do
print*,'***************************************************************'
print*,'*** 2D PERIODIC ***'
print*,'***************************************************************'
!call test_process_2d(sll_p_periodic,sll_p_periodic, passed_test)
print*,'***************************************************************'
print*,'*** 2D OPEN  ***'
print*,'***************************************************************'
!call test_process_2d(sll_p_open,sll_p_open, passed_test)
if (passed_test) then
   print *, 'PASSED'
else
   print *, 'FAILED'
end if

contains

  subroutine test_process_1d(bc_type,deg,passed_test)

    type(sll_t_bspline_interpolation_1d) :: bspline_1d
    sll_int32, intent(in) :: bc_type
    logical :: passed_test
    sll_int32  :: deg
    sll_real64 :: x_min   = 0.0_f64
    sll_real64 :: x_max   = 1.0_f64


    sll_real64, dimension(:), allocatable :: x
    sll_real64, dimension(:), allocatable :: xx
    sll_real64, dimension(:), allocatable :: y
    sll_real64, dimension(:), allocatable :: gtau
    sll_real64, dimension(:), allocatable :: htau

    sll_int32,  parameter                 :: npts = 5 ! defines bsplines
    sll_int32,  parameter                 :: n = 27    ! number of evaluation points
    sll_real64                            :: h

    SLL_ALLOCATE(x(n),ierr)
    SLL_ALLOCATE(y(n),ierr)
    SLL_ALLOCATE(xx(n),ierr)

    ! define uniform evaluation grid
    h = (x_max-x_min)/(n-1)
    do i = 1, n
       x(i) = x_min + (i-1)*h
    end do
    ! define set of random numbers for evaluation
    call random_number(xx)
    xx = xx * (x_max-x_min)

    print*,'+++++++++++++++++++++++++++++++++'
    print*,'*** Spline degree = ', deg
    print*,'+++++++++++++++++++++++++++++++++'
    call sll_s_bspline_interpolation_1d_init( bspline_1d, npts, deg, x_min, x_max, bc_type)
    print*, 'bspline_interpolation_1d_init constructed'

    SLL_ALLOCATE(gtau(bspline_1d%n),ierr)
    print*, '------------------------------------------'
    print*, 'Test on cosinus at uniformly spaced points'
    print*, '------------------------------------------'
    gtau = cos(2*sll_p_pi*bspline_1d%tau)
    !print*, 'tau:  ', bspline_1d%tau
    !print*, 'gtau: ', gtau

    call cpu_time(t0)
    call sll_s_compute_bspline_1d(bspline_1d, gtau)
    call cpu_time(t1)
    do j = 1,n
       !print*, x_max, x(j), bspline_1d%bsp%knots(npts)
       y(j) = sll_f_interpolate_value_1d( bspline_1d, x(j))
       !print*, j, x(j), cos(2*sll_p_pi*x(j)), y(j), y(j) - cos(2*sll_p_pi*x(j))
    end do
    err1 = maxval(abs(y-cos(2*sll_p_pi*x)))
    if (err1 > tol) then
       print*,'-------> Test failed in sll_f_interpolate_value_1d'
       passed_test = .false.
    end if
    print*, " sll_f_interpolate_value_1d: average error =             ", &
         sum(abs(y-cos(2*sll_p_pi*x)))/real(n,f64)
    print*, " sll_f_interpolate_value_1d: maximum error =             ", &
         maxval(abs(y-cos(2*sll_p_pi*x)))
    call cpu_time(t2)
    !.......
    do j = 1,nstep
       call sll_s_interpolate_array_values_1d( bspline_1d, n, x, y)
    end do
    err1 = maxval(abs(y-cos(2*sll_p_pi*x)))
    if (err1 > tol) then
       print*,'-------> Test failed in sll_s_interpolate_array_values_1d'
       passed_test = .false.
    end if
    print*, " sll_s_interpolate_array_values_1d: average error =      ", &
         sum(abs(y-cos(2*sll_p_pi*x)))/real(n,f64)
    print*, " sll_s_interpolate_array_values_1d: maximum error =      ", &
         maxval(abs(y-cos(2*sll_p_pi*x)))
    call cpu_time(t3)
    !......
    do j = 1,n
       y(j) = sll_f_interpolate_derivative_1d( bspline_1d, x(j))
    end do
    err1 = maxval(abs(y+2*sll_p_pi*sin(2*sll_p_pi*x)))
    if (err1 > tol) then
       print*,'-------> Test failed in sll_f_interpolate_derivative_1d'
       passed_test = .false.
    end if
    print*, " sll_f_interpolate_derivative_1d: average error =        ", &
         sum(abs(y+2*sll_p_pi*sin(2*sll_p_pi*x)))/real(n,f64)
    print*, " sll_f_interpolate_derivative_1d: maximum error =        ", &
         maxval(abs(y+2*sll_p_pi*sin(2*sll_p_pi*x)))
    call cpu_time(t4)
    !......
    do j = 1,nstep
       call sll_s_interpolate_array_derivatives_1d( bspline_1d, n, x, y)
    end do
    err1 = maxval(abs(y+2*sll_p_pi*sin(2*sll_p_pi*x)))
    if (err1 > tol) then
       print*,'-------> Test failed in sll_s_interpolate_array_derivatives_1d'
       passed_test = .false.
    end if
    print*, " sll_s_interpolate_array_derivatives_1d: average error = ", &
         sum(abs(y+2*sll_p_pi*sin(2*sll_p_pi*x)))/real(n,f64)
    print*, " sll_s_interpolate_array_derivatives_1d: maximum error = ", &
         maxval(abs(y+2*sll_p_pi*sin(2*sll_p_pi*x)))
    call cpu_time(t5)

    print*, ' ------------------------------------------------------- '
    print*, ' CPU time '
    print*, ' time spent to compute interpolants           : ', t1-t0
    print*, ' time spent to interpolate single values      : ', t2-t1
    print*, ' time spent to interpolate array values       : ', t3-t2
    print*, ' time spent to interpolate single derivatives : ', t4-t3
    print*, ' time spent to interpolate array derivatives  : ', t5-t4
    print*, ' ------------------------------------------------------- '

    print*, 'Test on sinus at random points'
    print*, ' -----------------------------'
    SLL_ALLOCATE(htau(bspline_1d%n),ierr)
    htau = sin(2*sll_p_pi*bspline_1d%tau)
    call sll_s_compute_bspline_1d(bspline_1d, htau)

    do j = 1,n
       y(j) = sll_f_interpolate_value_1d( bspline_1d, xx(j))
    end do
    err1 = maxval(abs(y-sin(2*sll_p_pi*xx)))
    if (err1 > tol) then
       print*,'-------> Test failed in sll_f_interpolate_value_1d'
       passed_test = .false.
    end if
    print*, " sll_f_interpolate_value_1d: average error =             ", &
         sum(abs(y-sin(2*sll_p_pi*xx)))/real(n,f64)
    print*, " sll_f_interpolate_value_1d: maximum error =             ", &
         maxval(abs(y-sin(2*sll_p_pi*xx)))
    call cpu_time(t2)
    !.......
    do j = 1,nstep
       call sll_s_interpolate_array_values_1d( bspline_1d, n, xx, y)
    end do
    err1 = maxval(abs(y-sin(2*sll_p_pi*xx)))
    if (err1 > tol) then
       print*,'-------> Test failed in sll_s_interpolate_array_values_1d'
       passed_test = .false.
    end if
    print*, " sll_s_interpolate_array_values_1d: average error =      ", &
         sum(abs(y-sin(2*sll_p_pi*xx)))/real(n,f64)
    print*, " sll_s_interpolate_array_values_1d: maximum error =      ", &
         maxval(abs(y-sin(2*sll_p_pi*xx)))
    call cpu_time(t3)
    !......
    do j = 1,n
       y(j) = sll_f_interpolate_derivative_1d( bspline_1d, xx(j))
    end do
    err1 = maxval(abs(y-2*sll_p_pi*cos(2*sll_p_pi*xx)))
    if (err1 > tol) then
       print*,'-------> Test failed in sll_f_interpolate_derivative_1d'
       passed_test = .false.
    end if
    print*, " sll_f_interpolate_derivative_1d: average error =        ", &
         sum(abs(y-2*sll_p_pi*cos(2*sll_p_pi*xx)))/real(n,f64)
    print*, " sll_f_interpolate_derivative_1d: maximum error =        ", &
         maxval(abs(y-2*sll_p_pi*cos(2*sll_p_pi*xx)))
    call cpu_time(t4)
    !......
    do j = 1,nstep
       call sll_s_interpolate_array_derivatives_1d( bspline_1d, n, xx, y)
    end do
    err1 = maxval(abs(y-2*sll_p_pi*cos(2*sll_p_pi*xx)))
    if (err1 > tol) then
       print*,'-------> Test failed in sll_s_interpolate_array_derivatives_1d'
       passed_test = .false.
    end if
    print*, " sll_s_interpolate_array_derivatives_1d: average error = ", &
         sum(abs(y-2*sll_p_pi*cos(2*sll_p_pi*xx)))/real(n,f64)
    print*, " sll_s_interpolate_array_derivatives_1d: maximum error = ", &
         maxval(abs(y-2*sll_p_pi*cos(2*sll_p_pi*xx)))
    call cpu_time(t5)

    !......
    SLL_DEALLOCATE_ARRAY(gtau,ierr)
    SLL_DEALLOCATE_ARRAY(htau,ierr)
    call sll_s_bspline_interpolation_1d_free(bspline_1d)

  end subroutine test_process_1d

!!$subroutine test_process_2d(bc1_type, bc2_type)
!!$
!!$ type(sll_t_bspline_interpolation_2d) :: bspline_2d  
!!$  sll_int32, intent(in)   :: bc1_type
!!$  sll_int32, intent(in)   :: bc2_type
!!$
!!$  sll_real64, allocatable :: ftau(:,:)
!!$  sll_real64, allocatable :: gtau(:,:)
!!$  sll_real64, allocatable :: tau1(:,:)
!!$  sll_real64, allocatable :: tau2(:,:)
!!$
!!$  sll_int32,  parameter   :: n1 = 64
!!$  sll_int32,  parameter   :: n2 = 64
!!$  sll_real64              :: sl1(n2)      ! slopes at boundaries
!!$  sll_real64              :: sr1(n2)      ! slopes at boundaries
!!$  sll_real64              :: sl2(n1)      ! slopes at boundaries
!!$  sll_real64              :: sr2(n1)      ! slopes at boundaries
!!$
!!$  sll_real64              :: f
!!$
!!$  sll_real64, parameter   :: dpi = 2*sll_p_pi
!!$
!!$  sll_real64, parameter   :: x1_min = 0.0_f64
!!$  sll_real64, parameter   :: x1_max = 1.0_f64
!!$  sll_real64, parameter   :: x2_min = 0.0_f64
!!$  sll_real64, parameter   :: x2_max = 1.0_f64
!!$
!!$  bspline_2d => sll_f_new_bspline_2d( n1, d, x1_min, x1_max, bc1_type, &
!!$       n2, d, x2_min, x2_max, bc2_type  )
!!$  print*, 'bspline_2d allocated'
!!$
!!$  allocate(ftau(n1,n2), gtau(n1,n2))
!!$  allocate(tau1(n1,n2), tau2(n1,n2))
!!$
!!$  do j = 1, n2
!!$     do i = 1, n1
!!$        tau1(i,j) = bspline_2d%bs1%tau(i)
!!$        tau2(i,j) = bspline_2d%bs2%tau(j)
!!$     end do
!!$  end do
!!$
!!$  ftau = cos(2*sll_p_pi*tau1) * cos(2*sll_p_pi*tau2)
!!$
!!$  sl1 = - sin(2*sll_p_pi*x1_min) * 2.0_f64*sll_p_pi * cos(2*sll_p_pi*tau2( 1, :))
!!$  sr1 = - sin(2*sll_p_pi*x1_max) * 2.0_f64*sll_p_pi * cos(2*sll_p_pi*tau2(n1, :))
!!$  sl2 = - sin(2*sll_p_pi*x2_min) * 2.0_f64*sll_p_pi * cos(2*sll_p_pi*tau1( :, 1))
!!$  sr2 = - sin(2*sll_p_pi*x2_max) * 2.0_f64*sll_p_pi * cos(2*sll_p_pi*tau1( :,n2))
!!$
!!$  call cpu_time(t0)
!!$  call sll_o_compute_bspline_2d(bspline_2d, ftau, sl1, sr1, sl2, sr2)
!!$  call cpu_time(t1)
!!$  do j = 1,nstep
!!$     call sll_s_interpolate_array_values_2d( bspline_2d, n1, n2, ftau, gtau, 0, 0)
!!$  end do
!!$  err1 = sum(abs(gtau-cos(dpi*tau1)*cos(dpi*tau2)))/real(n1*n2,f64)
!!$  err2 = maxval(abs(gtau-cos(dpi*tau1)*cos(dpi*tau2)))
!!$  call cpu_time(t2)
!!$  do j = 1,nstep
!!$     call sll_s_interpolate_array_values_2d( bspline_2d, n1, n2, ftau, gtau, 1, 0)
!!$  end do
!!$  err3 = sum(abs(gtau+dpi*sin(dpi*tau1)*cos(dpi*tau2)))/real(n1*n2,f64)
!!$  err4 = maxval(abs(gtau+dpi*sin(dpi*tau1)*cos(dpi*tau2)))
!!$  call cpu_time(t3)
!!$  do j = 1,nstep
!!$     call sll_s_interpolate_array_values_2d( bspline_2d, n1, n2, ftau, gtau, 0, 1)
!!$  end do
!!$  err5 = sum(abs(gtau+dpi*sin(dpi*tau2)*cos(dpi*tau1)))/real(n1*n2,f64)
!!$  err6 = maxval(abs(gtau+dpi*sin(dpi*tau2)*cos(dpi*tau1)))
!!$  call cpu_time(t4)
!!$
!!$  print*, "-----------------------------------------------------------"
!!$  print*, "ERRORS USING ARRAYS FUNCTIONS -----------------------------"
!!$  print*, "-----------------------------------------------------------"
!!$  print*, " average                                   : ", err1
!!$  print*, " maximum                                   : ", err2
!!$  print*, " average x1 derivatives                    : ", err3
!!$  print*, " maximum x1 derivatives                    : ", err4
!!$  print*, " average x2 derivatives                    : ", err5
!!$  print*, " maximum x2 derivatives                    : ", err6
!!$  print*, "-----------------------------------------------------------"
!!$  print*, "CPU TIME USING ARRAYS FUNCTIONS ---------------------------"
!!$  print*, "-----------------------------------------------------------"
!!$  print*, ' compute interpolants             : ', t1-t0
!!$  print*, ' interpolate array values         : ', t2-t1
!!$  print*, ' interpolate array x1 derivatives : ', t3-t2
!!$  print*, ' interpolate array x2 derivatives : ', t4-t3
!!$  print*, "-----------------------------------------------------------"
!!$
!!$  call cpu_time(t0)
!!$  call sll_o_compute_bspline_2d(bspline_2d, ftau, sl1, sr1, sl2, sr2)
!!$  call cpu_time(t1)
!!$  err1 = 0.0_f64
!!$  do j = 1, n2
!!$     do i = 1, n1
!!$        f = sll_f_interpolate_value_2d(bspline_2d,tau1(i,j),tau2(i,j),0,0)
!!$        err1 = err1 + abs(f-ftau(i,j))
!!$        write(10,*) tau1(i,j), tau2(i,j), f, ftau(i,j)
!!$     end do
!!$     write(10,*) 
!!$  end do
!!$  call cpu_time(t2)
!!$  err2 = 0.0_f64
!!$  do j = 1, n2
!!$     do i = 1, n1
!!$        f = sll_f_interpolate_value_2d(bspline_2d,tau1(i,j),tau2(i,j),1,0)
!!$        err2 = err2 + abs(f+dpi*sin(dpi*tau1(i,j))*cos(dpi*tau2(i,j)))
!!$     end do
!!$  end do
!!$  call cpu_time(t3)
!!$  err3 = 0.0_f64
!!$  do j = 1, n2
!!$     do i = 1, n1
!!$        f = sll_f_interpolate_value_2d(bspline_2d,tau1(i,j),tau2(i,j),0,1)
!!$        err3 = err3 + abs(f+dpi*sin(dpi*tau2(i,j))*cos(dpi*tau1(i,j)))
!!$     end do
!!$  end do
!!$  call cpu_time(t4)
!!$
!!$  print*, "----------------------------------------------------------"
!!$  print*, "ERRORS USING POINT VALUE FUNCTION ------------------------"
!!$  print*, "----------------------------------------------------------"
!!$  print*, " values                             : ", err1/(n1*n2)
!!$  print*, " x1 derivatives                     : ", err2/(n1*n2)
!!$  print*, " x2 derivatives                     : ", err3/(n1*n2)
!!$  print*, "----------------------------------------------------------"
!!$  print*, "CPU TIME USING POINT VALUE FUCNTION ----------------------"
!!$  print*, "----------------------------------------------------------"
!!$  print*, ' compute interpolants      : ', t1-t0
!!$  print*, ' interpolate_from_interpolant_value         : ', t2-t1
!!$  print*, ' interpolate_x1_derivative : ', t3-t2
!!$  print*, ' interpolate_x2_derivative : ', t4-t3
!!$  print*, "----------------------------------------------------------"
!!$
!!$  deallocate(tau1, tau2)
!!$  deallocate(ftau, gtau)
!!$
!!$end subroutine test_process_2d

end program test_bsplines_1d
