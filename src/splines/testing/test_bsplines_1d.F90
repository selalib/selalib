!PN This is a general test for bsplines 1d and 2d
!PN This test is used to test performance and precision
program test_bsplines_1d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_hermite, &
    sll_p_periodic

  use sll_m_bsplines, only: &
    sll_s_compute_bspline_1d, &
    sll_o_compute_bspline_2d, &
    sll_s_interpolate_array_derivatives_1d, &
    sll_s_interpolate_array_values_1d, &
    sll_s_interpolate_array_values_2d, &
    sll_f_interpolate_derivative_1d, &
    sll_f_interpolate_value_1d, &
    sll_f_interpolate_value_2d, &
    sll_f_new_bspline_1d, &
    sll_f_new_bspline_2d, &
    sll_t_bspline_1d, &
    sll_t_bspline_2d, &
    sll_s_update_bspline_1d

  use sll_m_constants, only: &
    sll_p_pi

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

type(sll_t_bspline_1d), pointer :: bspline_1d
type(sll_t_bspline_2d), pointer :: bspline_2d

sll_real64                    :: err1
sll_real64                    :: err2
sll_real64                    :: err3
sll_real64                    :: err4
sll_real64                    :: err5
sll_real64                    :: err6
sll_int32                     :: i
sll_int32                     :: j
sll_int32,  parameter         :: d = 3
sll_int32,  parameter         :: nstep = 1

sll_int32,  parameter         :: m = 2
sll_real64                    :: t0, t1, t2, t3, t4
sll_int32                     :: ierr

print*,'***************************************************************'
print*,'*** 1D PERIODIC ***'
print*,'***************************************************************'
call test_process_1d(sll_p_periodic)
print*,'***************************************************************'
print*,'*** 1D HERMITE ***'
print*,'***************************************************************'
call test_process_1d(sll_p_hermite)
print*,'***************************************************************'
print*,'*** 2D PERIODIC ***'
print*,'***************************************************************'
call test_process_2d(sll_p_periodic,sll_p_periodic)
print*,'***************************************************************'
print*,'*** 2D HERMITE  ***'
print*,'***************************************************************'
call test_process_2d(sll_p_hermite,sll_p_hermite)
print*,'PASSED'

contains

subroutine test_process_1d(bc_type)

  sll_int32, intent(in) :: bc_type

  sll_real64 :: tau_min   = 0.0_f64
  sll_real64 :: tau_max   = 1.0_f64
  sll_real64 :: slope_min = 0.0_f64
  sll_real64 :: slope_max = 0.0_f64

  sll_real64, dimension(:), allocatable :: x
  sll_real64, dimension(:), allocatable :: y
  sll_real64, dimension(:), allocatable :: gtau
  sll_real64, dimension(:), allocatable :: htau

  sll_int32,  parameter                 :: n = 1024
  sll_real64                            :: h
  
  SLL_ALLOCATE(x(n),ierr)
  SLL_ALLOCATE(y(n),ierr)
  SLL_ALLOCATE(gtau(n),ierr)
  SLL_ALLOCATE(htau(n),ierr)
  
  h = 1.0_f64/(n-1)
  do i = 1, n
    x(i) = (i-1)*h
  end do

  if (bc_type == sll_p_periodic) print*, "Periodic Bspline"
  bspline_1d => sll_f_new_bspline_1d( n, d, tau_min, tau_max, bc_type)
  print*, 'bspline_1d allocated'
  
  gtau = cos(2*sll_p_pi*bspline_1d%tau)
  slope_min = -sin(2.0_f64*sll_p_pi*tau_min)*2.0_f64*sll_p_pi
  slope_max = -sin(2.0_f64*sll_p_pi*tau_max)*2.0_f64*sll_p_pi
  call cpu_time(t0)
  call sll_s_compute_bspline_1d(bspline_1d, gtau, slope_min, slope_max)
  call cpu_time(t1)
  do j = 1,nstep
    call sll_s_interpolate_array_values_1d( bspline_1d, n, x, y)
  end do
  print*, " average values error = ", sum(abs(y-cos(2*sll_p_pi*x)))/real(n,f64)
  print*, " maximum values error = ", maxval(abs(y-cos(2*sll_p_pi*x)))
  call cpu_time(t2)
  do j = 1,nstep
    call sll_s_interpolate_array_derivatives_1d( bspline_1d, n, x, y)
  end do
  print*, " average derivatives error = ", sum(abs(y+2*sll_p_pi*sin(2*sll_p_pi*x)))/real(n,f64)
  print*, " maximum derivatives error = ", maxval(abs(y+2*sll_p_pi*sin(2*sll_p_pi*x)))
  call cpu_time(t3)

  print*, ' ------------------------------------------------------- '
  print*, ' CPU time '
  print*, ' time spent to compute interpolants          : ', t1-t0
  print*, ' time spent to interpolate array values      : ', t2-t1
  print*, ' time spent to interpolate array derivatives : ', t3-t2
  print*, ' ------------------------------------------------------- '
  
  htau = sin(2*sll_p_pi*bspline_1d%tau)
  call sll_s_update_bspline_1d(bspline_1d, htau)
  call random_number(x)
  x = x * (tau_max-tau_min)
  call sll_s_interpolate_array_values_1d( bspline_1d, n, x, y)
  print*, "L2 norm error = ", sqrt(sum((y-sin(2*sll_p_pi*x))**2*h))
  call sll_s_interpolate_array_derivatives_1d( bspline_1d, n, x, y)
  print*, "H1 norm error = ", sqrt(sum((y-2*sll_p_pi*cos(2*sll_p_pi*x))**2*h))
  
  call cpu_time(t0)
  do j = 1,nstep
    err1 = 0.0_f64
    do i = 1, n
      err1 = err1 + &
        abs(sll_f_interpolate_value_1d(bspline_1d,x(i))-sin(2*sll_p_pi*x(i))) 
    end do
  end do
  call cpu_time(t1)
  do j = 1,nstep
    err2 = 0.0_f64
    do i = 1, n
      err2 = err2 + &
        abs(sll_f_interpolate_derivative_1d(bspline_1d,x(i))-2*sll_p_pi*cos(2*sll_p_pi*x(i))) 
    end do
  end do
  call cpu_time(t2)
  
  print*, "-------------------------------------------------"
  print*, " values error = ", err1 / n
  print*, " derivatives error = ", err2 / n
  print*, ' time spent in interpolate_from_interpolant_value      : ', t1-t0
  print*, ' time spent in interpolate_derivative : ', t2-t1
  print*, "-------------------------------------------------"

end subroutine test_process_1d

subroutine test_process_2d(bc1_type, bc2_type)

  sll_int32, intent(in)   :: bc1_type
  sll_int32, intent(in)   :: bc2_type

  sll_real64, allocatable :: ftau(:,:)
  sll_real64, allocatable :: gtau(:,:)
  sll_real64, allocatable :: tau1(:,:)
  sll_real64, allocatable :: tau2(:,:)

  sll_int32,  parameter   :: n1 = 64
  sll_int32,  parameter   :: n2 = 64
  sll_real64              :: sl1(n2)      ! slopes at boundaries
  sll_real64              :: sr1(n2)      ! slopes at boundaries
  sll_real64              :: sl2(n1)      ! slopes at boundaries
  sll_real64              :: sr2(n1)      ! slopes at boundaries

  sll_real64              :: f

  sll_real64, parameter   :: dpi = 2*sll_p_pi

  sll_real64, parameter   :: x1_min = 0.0_f64
  sll_real64, parameter   :: x1_max = 1.0_f64
  sll_real64, parameter   :: x2_min = 0.0_f64
  sll_real64, parameter   :: x2_max = 1.0_f64
  
  bspline_2d => sll_f_new_bspline_2d( n1, d, x1_min, x1_max, bc1_type, &
                                n2, d, x2_min, x2_max, bc2_type  )
  print*, 'bspline_2d allocated'

  allocate(ftau(n1,n2), gtau(n1,n2))
  allocate(tau1(n1,n2), tau2(n1,n2))
  
  do j = 1, n2
    do i = 1, n1
      tau1(i,j) = bspline_2d%bs1%tau(i)
      tau2(i,j) = bspline_2d%bs2%tau(j)
    end do
  end do
  
  ftau = cos(2*sll_p_pi*tau1) * cos(2*sll_p_pi*tau2)

  sl1 = - sin(2*sll_p_pi*x1_min) * 2.0_f64*sll_p_pi * cos(2*sll_p_pi*tau2( 1, :))
  sr1 = - sin(2*sll_p_pi*x1_max) * 2.0_f64*sll_p_pi * cos(2*sll_p_pi*tau2(n1, :))
  sl2 = - sin(2*sll_p_pi*x2_min) * 2.0_f64*sll_p_pi * cos(2*sll_p_pi*tau1( :, 1))
  sr2 = - sin(2*sll_p_pi*x2_max) * 2.0_f64*sll_p_pi * cos(2*sll_p_pi*tau1( :,n2))

  call cpu_time(t0)
  call sll_o_compute_bspline_2d(bspline_2d, ftau, sl1, sr1, sl2, sr2)
  call cpu_time(t1)
  do j = 1,nstep
    call sll_s_interpolate_array_values_2d( bspline_2d, n1, n2, ftau, gtau, 0, 0)
  end do
  err1 = sum(abs(gtau-cos(dpi*tau1)*cos(dpi*tau2)))/real(n1*n2,f64)
  err2 = maxval(abs(gtau-cos(dpi*tau1)*cos(dpi*tau2)))
  call cpu_time(t2)
  do j = 1,nstep
    call sll_s_interpolate_array_values_2d( bspline_2d, n1, n2, ftau, gtau, 1, 0)
  end do
  err3 = sum(abs(gtau+dpi*sin(dpi*tau1)*cos(dpi*tau2)))/real(n1*n2,f64)
  err4 = maxval(abs(gtau+dpi*sin(dpi*tau1)*cos(dpi*tau2)))
  call cpu_time(t3)
  do j = 1,nstep
    call sll_s_interpolate_array_values_2d( bspline_2d, n1, n2, ftau, gtau, 0, 1)
  end do
  err5 = sum(abs(gtau+dpi*sin(dpi*tau2)*cos(dpi*tau1)))/real(n1*n2,f64)
  err6 = maxval(abs(gtau+dpi*sin(dpi*tau2)*cos(dpi*tau1)))
  call cpu_time(t4)

  print*, "-----------------------------------------------------------"
  print*, "ERRORS USING ARRAYS FUNCTIONS -----------------------------"
  print*, "-----------------------------------------------------------"
  print*, " average                                   : ", err1
  print*, " maximum                                   : ", err2
  print*, " average x1 derivatives                    : ", err3
  print*, " maximum x1 derivatives                    : ", err4
  print*, " average x2 derivatives                    : ", err5
  print*, " maximum x2 derivatives                    : ", err6
  print*, "-----------------------------------------------------------"
  print*, "CPU TIME USING ARRAYS FUNCTIONS ---------------------------"
  print*, "-----------------------------------------------------------"
  print*, ' compute interpolants             : ', t1-t0
  print*, ' interpolate array values         : ', t2-t1
  print*, ' interpolate array x1 derivatives : ', t3-t2
  print*, ' interpolate array x2 derivatives : ', t4-t3
  print*, "-----------------------------------------------------------"
  
  call cpu_time(t0)
  call sll_o_compute_bspline_2d(bspline_2d, ftau, sl1, sr1, sl2, sr2)
  call cpu_time(t1)
  err1 = 0.0_f64
  do j = 1, n2
    do i = 1, n1
      f = sll_f_interpolate_value_2d(bspline_2d,tau1(i,j),tau2(i,j),0,0)
      err1 = err1 + abs(f-ftau(i,j))
      write(10,*) tau1(i,j), tau2(i,j), f, ftau(i,j)
    end do
    write(10,*) 
  end do
  call cpu_time(t2)
  err2 = 0.0_f64
  do j = 1, n2
    do i = 1, n1
      f = sll_f_interpolate_value_2d(bspline_2d,tau1(i,j),tau2(i,j),1,0)
      err2 = err2 + abs(f+dpi*sin(dpi*tau1(i,j))*cos(dpi*tau2(i,j)))
    end do
  end do
  call cpu_time(t3)
  err3 = 0.0_f64
  do j = 1, n2
    do i = 1, n1
      f = sll_f_interpolate_value_2d(bspline_2d,tau1(i,j),tau2(i,j),0,1)
      err3 = err3 + abs(f+dpi*sin(dpi*tau2(i,j))*cos(dpi*tau1(i,j)))
    end do
  end do
  call cpu_time(t4)
  
  print*, "----------------------------------------------------------"
  print*, "ERRORS USING POINT VALUE FUNCTION ------------------------"
  print*, "----------------------------------------------------------"
  print*, " values                             : ", err1/(n1*n2)
  print*, " x1 derivatives                     : ", err2/(n1*n2)
  print*, " x2 derivatives                     : ", err3/(n1*n2)
  print*, "----------------------------------------------------------"
  print*, "CPU TIME USING POINT VALUE FUCNTION ----------------------"
  print*, "----------------------------------------------------------"
  print*, ' compute interpolants      : ', t1-t0
  print*, ' interpolate_from_interpolant_value         : ', t2-t1
  print*, ' interpolate_x1_derivative : ', t3-t2
  print*, ' interpolate_x2_derivative : ', t4-t3
  print*, "----------------------------------------------------------"

  deallocate(tau1, tau2)
  deallocate(ftau, gtau)

end subroutine test_process_2d

end program test_bsplines_1d
