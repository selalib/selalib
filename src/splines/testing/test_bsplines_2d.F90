!PN This is a general test for bsplines 1d and 2d
!PN This test is used to test performance and precision
program test_bsplines_2d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

use sll_m_boundary_condition_descriptors, only: &
     sll_p_hermite, &
     sll_p_greville, &
     sll_p_periodic, &
     sll_p_mirror

use sll_m_bspline_interpolation

use sll_m_constants, only: pi => sll_p_pi

implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sll_real64                    :: err1
sll_int32                     :: i
sll_int32                     :: j
sll_int32,  parameter         :: nstep = 1
sll_real64,  parameter        :: tol = 1.0d-3    ! tolerance for tests
sll_real64                    :: t0, t1, t2, t3, t4, t5, t6, t7
sll_int32                     :: ierr
sll_int32                     :: deg
logical                       :: passed_test

passed_test = .true.
  
print*,'***************************************************************'
print*,'*** 2D PERIODIC-PERIODIC ***'
print*,'***************************************************************'
do deg=3,9
   call test_process_2d(sll_p_periodic,deg, passed_test)
end do
print*,'***************************************************************'
print*,'*** 2D GREVILLE-GREVILLE ***'
print*,'***************************************************************'
do deg=3,9
   call test_process_2d(sll_p_greville,deg, passed_test)
end do
print*,'***************************************************************'
print*,'*** 2D HERMITE-HERMITE ***'
print*,'***************************************************************'
do deg=3,9
   call test_process_2d(sll_p_hermite,deg, passed_test)
end do
print*,'***************************************************************'
print*,'*** 1D HERMITE WITH MIRROR KNOT POINTS***'
print*,'***************************************************************'
do deg=3,9
   call test_process_2d(sll_p_hermite,deg, passed_test, sll_p_mirror)
end do

if (passed_test) then
  print *, 'PASSED'
else
  print *, 'FAILED'
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine test_process_2d(bc_type,deg,passed_test,spline_bc_type)

type(sll_t_bspline_interpolation_2d) :: bspline_2d
sll_int32, intent(in)                :: bc_type
logical                              :: passed_test
sll_int32, optional                  :: spline_bc_type

! local variables
sll_int32  :: deg
sll_real64 :: x1_min   = 0.0_f64
sll_real64 :: x1_max   = 1.0_f64
sll_real64 :: x2_min   = 0.0_f64
sll_real64 :: x2_max   = 1.0_f64


sll_real64, dimension(:,:), allocatable :: x1
sll_real64, dimension(:,:), allocatable :: x2
sll_real64, dimension(:,:), allocatable :: y
sll_real64, dimension(:,:), allocatable :: xx
sll_real64, dimension(:,:), allocatable :: yy
sll_real64, dimension(:,:), allocatable :: gtau
sll_real64, dimension(:,:), allocatable :: htau
sll_real64, dimension(:), pointer       :: taux
sll_real64, dimension(:), pointer       :: tauy

sll_real64, dimension(:),   allocatable :: bc

sll_int32,  parameter                   :: npts1 = 51 ! defines bsplines
sll_int32,  parameter                   :: npts2 = 51 ! defines bsplines
sll_int32,  parameter                   :: n1 = 27    ! nb of evaluation pts
sll_int32,  parameter                   :: n2 = 27    ! nb of evaluation pts
sll_real64                              :: h1, h2

SLL_ALLOCATE(x1(n1,n2),ierr)
SLL_ALLOCATE(x2(n1,n2),ierr)
SLL_ALLOCATE(y(n1,n2),ierr)

! define uniform evaluation grid
h1 = (x1_max-x1_min)/(n1-1)
h2 = (x2_max-x2_min)/(n2-1)
do j = 1, n2
  do i = 1, n1
    x1(i,j) = x1_min + (i-1)*h1
    x2(i,j) = x2_min + (j-1)*h2
  end do
end do


print*,'+++++++++++++++++++++++++++++++++'
print*,'*** Spline degree = ', deg
print*,'+++++++++++++++++++++++++++++++++'

if (present(spline_bc_type)) then
  call sll_s_bspline_interpolation_2d_init( bspline_2d, &
     npts1, npts2, deg, deg, x1_min, x2_min, x1_max, x2_max, &
     bc_type, bc_type, spline_bc_type, spline_bc_type )
else 
  call sll_s_bspline_interpolation_2d_init( bspline_2d, &
     npts1, npts2, deg, deg, x1_min, x2_min, x1_max, x2_max, &
     bc_type, bc_type)
end if

taux => bspline_2d%bs1%tau
tauy => bspline_2d%bs2%tau

print*, 'bspline_interpolation_2d_init constructed'

SLL_ALLOCATE(gtau(bspline_2d%bs1%n, bspline_2d%bs2%n),ierr)
SLL_ALLOCATE(bc(2*(deg/2)),ierr)

print*, '------------------------------------------'
print*, 'Test on cosinus at uniformly spaced points'
print*, '------------------------------------------'
do j = 1, size(tauy)
  do i = 1, size(taux)
    gtau(i,j) = cos(2*pi*taux(i)) * cos(2*pi*tauy(j)) 
  end do
end do

call cpu_time(t0)

if (bc_type == sll_p_hermite) then
  ! Compute boundary conditions for Hermite case
  bc=0.0_f64
  if (modulo(deg,2) == 0) then
    bc(1) = 1.0_f64 
    do i=3,deg/2,2
      bc(i)=-4*pi**2*bc(i-2)
    end do
  else
    bc(2) = -4*pi**2
    do i=4,deg/2,2
      bc(i)=-4*pi**2*bc(i-2)
    end do
  end if
  call sll_s_compute_bspline_2d(bspline_2d, gtau, bc, bc, bc, bc)
else    
  call sll_s_compute_bspline_2d(bspline_2d, gtau)
end if

call cpu_time(t1)

do j = 1,n2
  do i = 1,n1
    y(i,j) = sll_f_interpolate_value_2d( bspline_2d, x1(i,j),x2(i,j))
  end do
end do

err1 = maxval(abs(y-cos(2*pi*x1)*cos(2*pi*x2)))
if (err1 > tol) then
  print*,'-------> Test failed in sll_f_interpolate_value_2d'
  passed_test = .false.
end if
print*, " sll_f_interpolate_value_2d: average error =             ", &
     sum(abs(y-cos(2*pi*x1)*cos(2*pi*x2)))/real(n1*n2,f64)
print*, " sll_f_interpolate_value_2d: maximum error =             ", &
     maxval(abs(y-cos(2*pi*x1)*cos(2*pi*x2)))

if (.not. passed_test) stop

call cpu_time(t2)

!.......
do j = 1,nstep
  call sll_s_interpolate_array_values_2d( bspline_2d, n1, n2, x1, x2, y)
end do
err1 = maxval(abs(y-cos(2*pi*x1)*cos(2*pi*x2)))
if (err1 > tol) then
  print*,'-------> Test failed in sll_s_interpolate_array_values_2d'
  passed_test = .false.
end if
print*, " sll_s_interpolate_array_values_2d: average error =      ", &
    sum(abs(y-cos(2*pi*x1)*cos(2*pi*x2)))/real(n1*n2,f64)
print*, " sll_s_interpolate_array_values_2d: maximum error =      ", &
    maxval(abs(y-cos(2*pi*x1)*cos(2*pi*x2)))

if (.not. passed_test) stop
call cpu_time(t3)

!......
do j = 1,n2
  do i = 1,n1
    y(i,j) = sll_f_interpolate_derivative_x1_2d( bspline_2d, x1(i,j), x2(i,j))
  end do
end do
err1 = maxval(abs(y+2*pi*sin(2*pi*x1)*cos(2*pi*x2)))
if (err1 > tol) then
  print*,'-------> Test failed in sll_f_interpolate_derivative_x1_2d'
  passed_test = .false.
end if
print*, " sll_f_interpolate_derivative_x1_2d: average error =        ", &
     sum(abs(y+2*pi*sin(2*pi*x1)*cos(2*pi*x2)))/real(n1*n2,f64)
print*, " sll_f_interpolate_derivative_x1_2d: maximum error =        ", &
     maxval(abs(y+2*pi*sin(2*pi*x1)*cos(2*pi*x2)))

if (.not. passed_test) stop
call cpu_time(t4)

!......
do j = 1,n2
  do i = 1,n1
    y(i,j) = sll_f_interpolate_derivative_x2_2d( bspline_2d, x1(i,j), x2(i,j))
  end do
end do
err1 = maxval(abs(y+2*pi*sin(2*pi*x2)*cos(2*pi*x1)))
if (err1 > tol) then
  print*,'-------> Test failed in sll_f_interpolate_derivative_x2_2d'
  passed_test = .false.
end if
print*, " sll_f_interpolate_derivative_x2_2d: average error =        ", &
     sum(abs(y+2*pi*sin(2*pi*x2)*cos(2*pi*x1)))/real(n1*n2,f64)
print*, " sll_f_interpolate_derivative_x2_2d: maximum error =        ", &
     maxval(abs(y+2*pi*sin(2*pi*x2)*cos(2*pi*x1)))

if (.not. passed_test) stop
call cpu_time(t5)

!......
do j = 1,nstep
  call sll_s_interpolate_array_derivatives_x1_2d( bspline_2d, n1, n2, x1, x2, y)
end do
err1 = maxval(abs(y+2*pi*sin(2*pi*x1)*cos(2*pi*x2)))
if (err1 > tol) then
  print*,'-------> Test failed in sll_s_interpolate_array_derivatives_x1_2d'
  passed_test = .false.
end if
print*, " sll_s_interpolate_array_derivatives_x1_2d: average error = ", &
    sum(abs(y+2*pi*sin(2*pi*x1)*cos(2*pi*x2)))/real(n1*n2,f64)
print*, " sll_s_interpolate_array_derivatives_x1_2d: maximum error = ", &
    maxval(abs(y+2*pi*sin(2*pi*x1)*cos(2*pi*x2)))

if (.not. passed_test) stop
call cpu_time(t6)

!......
do j = 1,nstep
  call sll_s_interpolate_array_derivatives_x2_2d( bspline_2d, n1, n2, x1, x2, y)
end do
err1 = maxval(abs(y+2*pi*sin(2*pi*x2)*cos(2*pi*x1)))
if (err1 > tol) then
  print*,'-------> Test failed in sll_s_interpolate_array_derivatives_x2_2d'
  passed_test = .false.
end if
print*, " sll_s_interpolate_array_derivatives_x2_2d: average error = ", &
    sum(abs(y+2*pi*sin(2*pi*x2)*cos(2*pi*x1)))/real(n1*n2,f64)
print*, " sll_s_interpolate_array_derivatives_x2_2d: maximum error = ", &
    maxval(abs(y+2*pi*sin(2*pi*x2)*cos(2*pi*x1)))

if (.not. passed_test) stop
call cpu_time(t7)

print*, ' ------------------------------------------------------- '
print*, ' CPU time '
print*, ' time spent to compute interpolants              : ', t1-t0
print*, ' time spent to interpolate single values         : ', t2-t1
print*, ' time spent to interpolate array values          : ', t3-t2
print*, ' time spent to interpolate single x1 derivatives : ', t4-t3
print*, ' time spent to interpolate single x2 derivatives : ', t5-t4
print*, ' time spent to interpolate array  x1 derivatives : ', t6-t5
print*, ' time spent to interpolate array  x2 derivatives : ', t7-t6
print*, ' ------------------------------------------------------- '

print*, 'Test on sinus at random points'
print*, ' -----------------------------'

SLL_ALLOCATE(xx(n1,n2),ierr)
SLL_ALLOCATE(yy(n1,n2),ierr)
call random_number(xx)
call random_number(yy)
xx = xx * (x1_max-x1_min)
yy = yy * (x2_max-x2_min)

SLL_ALLOCATE(htau(bspline_2d%bs1%n, bspline_2d%bs2%n),ierr)

do j = 1, bspline_2d%bs2%n
  do i = 1, bspline_2d%bs1%n
    htau(i,j) = sin(2*pi*taux(i))*sin(2*pi*tauy(j))
    write(12,*) taux(i), tauy(j), htau(i,j) 
  end do
  write(12,*)
end do
close(12)

if (bc_type == sll_p_hermite) then
  ! Compute boundary conditions for Hermite case
  bc=0.0_f64
  if (modulo(deg,2) == 0) then
    bc(2) = 2*pi
    do i=4,deg/2,2
      bc(i)=-4*pi**2*bc(i-2)
    end do
  else
    bc(1) = 2*pi
    do i=3,deg/2,2
      bc(i)=-4*pi**2*bc(i-2)
    end do
  end if
  call sll_s_compute_bspline_2d(bspline_2d, htau, bc, bc, bc, bc)
else    
  call sll_s_compute_bspline_2d(bspline_2d, htau)
end if

do j = 1,n2
  do i = 1,n1
    y(i,j) = sll_f_interpolate_value_2d( bspline_2d, xx(i,j), yy(i,j))
  end do
end do
err1 = maxval(abs(y-sin(2*pi*xx)*sin(2*pi*yy)))
if (err1 > tol) then
  print*,'-------> Test failed in sll_f_interpolate_value_2d'
  passed_test = .false.
end if
print*, " sll_f_interpolate_value_2d: average error =             ", &
     sum(abs(y-sin(2*pi*xx)*sin(2*pi*yy)))/real(n1*n2,f64)
print*, " sll_f_interpolate_value_2d: maximum error =             ", &
     maxval(abs(y-sin(2*pi*xx)*sin(2*pi*yy)))

if (.not. passed_test) stop

call cpu_time(t2)

!.......
do j = 1,nstep
  call sll_s_interpolate_array_values_2d( bspline_2d, n1, n2, xx, yy, y)
end do
err1 = maxval(abs(y-sin(2*pi*xx)*sin(2*pi*yy)))
if (err1 > tol) then
  print*,'-------> Test failed in sll_s_interpolate_array_values_2d'
  passed_test = .false.
end if
print*, " sll_s_interpolate_array_values_2d: average error =      ", &
    sum(abs(y-sin(2*pi*xx)*sin(2*pi*yy)))/real(n1*n2,f64)
print*, " sll_s_interpolate_array_values_2d: maximum error =      ", &
     maxval(abs(y-sin(2*pi*xx)*sin(2*pi*yy)))
if (.not. passed_test) stop

call cpu_time(t3)

!......
do j = 1,n2
  do i = 1,n1
    y(i,j) = sll_f_interpolate_derivative_x1_2d( bspline_2d, xx(i,j), yy(i,j))
  end do
end do
err1 = maxval(abs(y-2*pi*cos(2*pi*xx)*sin(2*pi*yy)))
if (err1 > tol) then
  print*,'-------> Test failed in sll_f_interpolate_derivative_x1_2d'
  passed_test = .false.
end if
print*, " sll_f_interpolate_derivative_x1_2d: average error =        ", &
     sum(abs(y-2*pi*cos(2*pi*xx)*sin(2*pi*yy)))/real(n1*n2,f64)
print*, " sll_f_interpolate_derivative_x1_2d: maximum error =        ", &
     maxval(abs(y-2*pi*cos(2*pi*xx)*sin(2*pi*yy)))

if (.not. passed_test) stop
call cpu_time(t4)

!......
do j = 1,n2
  do i = 1,n1
    y(i,j) = sll_f_interpolate_derivative_x2_2d( bspline_2d, xx(i,j), yy(i,j))
  end do
end do
err1 = maxval(abs(y-2*pi*cos(2*pi*yy)*sin(2*pi*xx)))
if (err1 > tol) then
  print*,'-------> Test failed in sll_f_interpolate_derivative_x2_2d'
  passed_test = .false.
end if
print*, " sll_f_interpolate_derivative_x2_2d: average error =        ", &
     sum(abs(y-2*pi*cos(2*pi*yy)*sin(2*pi*xx)))/real(n1*n2,f64)
print*, " sll_f_interpolate_derivative_x2_2d: maximum error =        ", &
     maxval(abs(y-2*pi*cos(2*pi*yy)*sin(2*pi*xx)))
if (.not. passed_test) stop

call cpu_time(t5)

!......
do j = 1,nstep
  call sll_s_interpolate_array_derivatives_x1_2d( bspline_2d, n1, n2, xx, yy, y)
end do
err1 = maxval(abs(y-2*pi*cos(2*pi*xx)*sin(2*pi*yy)))
if (err1 > tol) then
  print*,'-------> Test failed in sll_s_interpolate_array_derivatives_x1_2d'
  passed_test = .false.
end if
print*, " sll_s_interpolate_array_derivatives_x1_2d: average error = ", &
     sum(abs(y-2*pi*cos(2*pi*xx)*sin(2*pi*yy)))/real(n1*n2,f64)
print*, " sll_s_interpolate_array_derivatives_x1_2d: maximum error = ", &
     maxval(abs(y-2*pi*cos(2*pi*xx)*sin(2*pi*yy)))

if (.not. passed_test) stop
call cpu_time(t6)

!......
do j = 1,nstep
  call sll_s_interpolate_array_derivatives_x2_2d( bspline_2d, n1, n2, xx, yy, y)
end do
err1 = maxval(abs(y-2*pi*cos(2*pi*yy)*sin(2*pi*xx)))
if (err1 > tol) then
  print*,'-------> Test failed in sll_s_interpolate_array_derivatives_x2_2d'
  passed_test = .false.
end if
print*, " sll_s_interpolate_array_derivatives_x2_2d: average error = ", &
     sum(abs(y-2*pi*cos(2*pi*yy)*sin(2*pi*xx)))/real(n1*n2,f64)
print*, " sll_s_interpolate_array_derivatives_x2_2d: maximum error = ", &
     maxval(abs(y-2*pi*cos(2*pi*yy)*sin(2*pi*xx)))

if (.not. passed_test) stop
call cpu_time(t7)
!......
print*, ' ------------------------------------------------------- '
print*, ' CPU time '
print*, ' time spent to compute interpolants              : ', t1-t0
print*, ' time spent to interpolate single values         : ', t2-t1
print*, ' time spent to interpolate array values          : ', t3-t2
print*, ' time spent to interpolate single x1 derivatives : ', t4-t3
print*, ' time spent to interpolate single x2 derivatives : ', t5-t4
print*, ' time spent to interpolate array  x1 derivatives : ', t6-t5
print*, ' time spent to interpolate array  x2 derivatives : ', t7-t6
print*, ' ------------------------------------------------------- '

SLL_DEALLOCATE_ARRAY(gtau,ierr)
SLL_DEALLOCATE_ARRAY(htau,ierr)
call sll_s_bspline_interpolation_2d_free(bspline_2d)

end subroutine test_process_2d

end program test_bsplines_2d
