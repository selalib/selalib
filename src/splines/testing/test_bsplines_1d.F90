!PN This is a general test for bsplines 1d and 2d
!PN This test is used to test performance and precision
program test_bspline_1d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

use sll_m_boundary_condition_descriptors, only: &
     sll_p_hermite, &
     sll_p_greville, &
     sll_p_periodic, &
     sll_p_mirror

use sll_m_spline_1d_non_uniform, only: &
  sll_s_bspline_1d_compute_interpolant, &
  sll_s_bspline_1d_eval_array_deriv, &
  sll_s_bspline_1d_eval_array, &
  sll_f_bspline_1d_eval_deriv, &
  sll_f_bspline_1d_eval, &
  sll_s_bspline_1d_init, &
  sll_s_bspline_1d_free, &
  sll_t_bspline_1d

use sll_m_constants, only: &
  sll_p_pi

implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sll_real64                    :: err1
sll_int32                     :: i
sll_int32                     :: j
sll_int32,  parameter         :: nstep = 1
sll_real64                    :: t0, t1, t2, t3, t4, t5
sll_int32                     :: ierr
sll_int32                     :: deg
logical                       :: passed_test
sll_real64                    :: tol(9)

! Set tolerances for various spline degrees
tol(1) = 4e-01_f64
tol(2) = 1e-02_f64
tol(3) = 1e-03_f64
tol(4) = 1e-05_f64
tol(5) = 1e-07_f64
tol(6) = 1e-08_f64
tol(7) = 1e-10_f64
tol(8) = 1e-11_f64
tol(9) = 1e-12_f64

! Initialize PASSED/FAILED condition
passed_test = .true.

print*,'***************************************************************'
print*,'*** 1D PERIODIC ***'
print*,'***************************************************************'
do deg=1,9
   call test_process_1d( sll_p_periodic, deg, tol(deg), passed_test)
end do
print*,
print*,'***************************************************************'
print*,'*** 1D GREVILLE ***'
print*,'***************************************************************'
do deg=1,9
   call test_process_1d( sll_p_greville, deg, tol(deg), passed_test )
end do
print*,
print*,'***************************************************************'
print*,'*** 1D HERMITE ***'
print*,'***************************************************************'
do deg=1,9
   call test_process_1d( sll_p_hermite, deg, tol(deg), passed_test )
end do
print*,
!print*,'***************************************************************'
!print*,'*** 1D HERMITE  WITH MIRROR KNOT POINTS***'
!print*,'***************************************************************'
!do deg=1,9
!   call test_process_1d( sll_p_hermite, deg, tol(deg), passed_test, sll_p_mirror )
!end do

if (passed_test) then
   print *, 'PASSED'
else
   print *, 'FAILED'
end if

contains

  subroutine test_process_1d( bc_type, deg, tol, passed_test, spline_bc_type )

    sll_int32 , intent(in   )           :: bc_type
    sll_int32 , intent(in   )           :: deg
    sll_real64, intent(in   )           :: tol
    logical   , intent(inout)           :: passed_test
    sll_int32 , intent(in   ), optional :: spline_bc_type

    ! local variables
    type(sll_t_bspline_1d) :: bspline_1d
    sll_real64 :: x_min   = 0.0_f64
    sll_real64 :: x_max   = 1.0_f64


    sll_real64, dimension(:), allocatable :: x
    sll_real64, dimension(:), allocatable :: xx
    sll_real64, dimension(:), allocatable :: y
    sll_real64, dimension(:), allocatable :: gtau
    sll_real64, dimension(:), allocatable :: htau
    sll_real64, dimension(:), allocatable :: bc

    sll_int32,  parameter                 :: npts = 51 ! defines bsplines
    sll_int32,  parameter                 :: n = 27    ! number of evaluation points
    sll_real64                            :: h

    SLL_ALLOCATE(x(n),ierr)
    SLL_ALLOCATE(y(n),ierr)
    SLL_ALLOCATE(xx(n),ierr)

    ! define uniform evaluation grid
    h = (x_max-x_min)/real(n-1,f64)
    do i = 1, n
       x(i) = x_min + real(i-1,f64)*h
    end do
    ! define set of random numbers for evaluation
    call random_number(xx)
    xx = x_min + xx * (x_max-x_min)

    print*,'+++++++++++++++++++++++++++++++++'
    print*,'*** Spline degree = ', deg
    print*,'+++++++++++++++++++++++++++++++++'
    if (present(spline_bc_type)) then
       call sll_s_bspline_1d_init( bspline_1d, npts, deg, x_min, x_max, &
            bc_type, bc_type, spline_bc_type )
    else
       call sll_s_bspline_1d_init( bspline_1d, npts, deg, x_min, x_max, &
            bc_type, bc_type )
    end if
    print*, 'bspline_1d_init constructed'

    SLL_ALLOCATE(gtau(size(bspline_1d%tau)),ierr)
    SLL_ALLOCATE(bc(deg/2),ierr)
    print*, '------------------------------------------'
    print*, 'Test on cosinus at uniformly spaced points'
    print*, '------------------------------------------'
    gtau = cos(2*sll_p_pi*bspline_1d%tau)
    !print*, 'tau:  ', bspline_1d%tau
    !print*, 'gtau: ', gtau

    call cpu_time(t0)
    if (bc_type == sll_p_hermite) then
       ! Compute boundary conditions for Hermite case
       bc=0.0_f64
       if (modulo(deg,2) == 0) then !------> even degree
          bc(1) = 1.0_f64
          do i=3,deg/2,2
             bc(i)=-4*sll_p_pi**2*bc(i-2)
          end do
       else !------------------------------>  odd degree
         if (deg>3) then
           bc(2) = -4*sll_p_pi**2
           do i=4,deg/2,2
             bc(i)=-4*sll_p_pi**2*bc(i-2)
           end do
         end if
       end if
       call sll_s_bspline_1d_compute_interpolant(bspline_1d, gtau, bc, bc)
    else
       call sll_s_bspline_1d_compute_interpolant(bspline_1d, gtau)
    end if
    call cpu_time(t1)
    do j = 1,n
       y(j) = sll_f_bspline_1d_eval( bspline_1d, x(j))
    end do
    err1 = maxval(abs(y-cos(2*sll_p_pi*x)))
    if (err1 > tol) then
       print*,'-------> Test failed in sll_f_bspline_1d_eval'
       passed_test = .false.
    end if
    print*, " sll_f_bspline_1d_eval              average error = ", &
         sum(abs(y-cos(2*sll_p_pi*x)))/real(n,f64)
    print*, " sll_f_bspline_1d_eval              maximum error = ", &
         maxval(abs(y-cos(2*sll_p_pi*x)))
    call cpu_time(t2)
    !.......
    do j = 1,nstep
       call sll_s_bspline_1d_eval_array( bspline_1d, x, y)
    end do
    err1 = maxval(abs(y-cos(2*sll_p_pi*x)))
    if (err1 > tol) then
       print*,'-------> Test failed in sll_s_bspline_1d_eval_array'
       passed_test = .false.
    end if
    print*, " sll_s_bspline_1d_eval_array       average error = ", &
         sum(abs(y-cos(2*sll_p_pi*x)))/real(n,f64)
    print*, " sll_s_bspline_1d_eval_array       maximum error = ", &
         maxval(abs(y-cos(2*sll_p_pi*x)))
    call cpu_time(t3)
    !......
    do j = 1,n
       y(j) = sll_f_bspline_1d_eval_deriv( bspline_1d, x(j))
    end do
    err1 = maxval(abs(y+2*sll_p_pi*sin(2*sll_p_pi*x)))
    if (err1 > tol) then
       print*,'-------> Test failed in sll_f_bspline_1d_eval_deriv'
       passed_test = .false.
    end if
    print*, " sll_f_bspline_1d_eval_deriv         average error = ", &
         sum(abs(y+2*sll_p_pi*sin(2*sll_p_pi*x)))/real(n,f64)
    print*, " sll_f_bspline_1d_eval_deriv         maximum error = ", &
         maxval(abs(y+2*sll_p_pi*sin(2*sll_p_pi*x)))
    call cpu_time(t4)
    !......
    do j = 1,nstep
       call sll_s_bspline_1d_eval_array_deriv( bspline_1d, x, y)
    end do
    err1 = maxval(abs(y+2*sll_p_pi*sin(2*sll_p_pi*x)))
    if (err1 > tol) then
       print*,'-------> Test failed in sll_s_bspline_1d_eval_array_deriv'
       passed_test = .false.
    end if
    print*, " sll_s_bspline_1d_eval_array_deriv  average error = ", &
         sum(abs(y+2*sll_p_pi*sin(2*sll_p_pi*x)))/real(n,f64)
    print*, " sll_s_bspline_1d_eval_array_deriv  maximum error = ", &
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
    SLL_ALLOCATE(htau(size(bspline_1d%tau)),ierr)
    htau = sin(2*sll_p_pi*bspline_1d%tau)
    if (bc_type == sll_p_hermite) then
       ! Compute boundary conditions for Hermite case
       bc=0.0_f64
       if (modulo(deg,2) == 0) then !------> even degree
         if (deg>2) then
           bc(2) = 2*sll_p_pi
           do i=4,deg/2,2
             bc(i)=-4*sll_p_pi**2*bc(i-2)
           end do
         end if
       else !------------------------------>  odd degree
         if (deg>1) then
           bc(1) = 2*sll_p_pi
           do i=3,deg/2,2
              bc(i)=-4*sll_p_pi**2*bc(i-2)
           end do
         end if
       end if
       call sll_s_bspline_1d_compute_interpolant(bspline_1d, htau, bc, bc)
    else
       call sll_s_bspline_1d_compute_interpolant(bspline_1d, htau)
    end if

    do j = 1,n
       y(j) = sll_f_bspline_1d_eval( bspline_1d, xx(j))
    end do
    err1 = maxval(abs(y-sin(2*sll_p_pi*xx)))
    if (err1 > tol) then
       print*,'-------> Test failed in sll_f_bspline_1d_eval'
       passed_test = .false.
    end if
    print*, " sll_f_bspline_1d_eval              average error = ", &
         sum(abs(y-sin(2*sll_p_pi*xx)))/real(n,f64)
    print*, " sll_f_bspline_1d_eval              maximum error = ", &
         maxval(abs(y-sin(2*sll_p_pi*xx)))
    call cpu_time(t2)
    !.......
    do j = 1,nstep
       call sll_s_bspline_1d_eval_array( bspline_1d, xx, y)
    end do
    err1 = maxval(abs(y-sin(2*sll_p_pi*xx)))
    if (err1 > tol) then
       print*,'-------> Test failed in sll_s_bspline_1d_eval_array'
       passed_test = .false.
    end if
    print*, " sll_s_bspline_1d_eval_array       average error = ", &
         sum(abs(y-sin(2*sll_p_pi*xx)))/real(n,f64)
    print*, " sll_s_bspline_1d_eval_array       maximum error = ", &
         maxval(abs(y-sin(2*sll_p_pi*xx)))
    call cpu_time(t3)
    !......
    do j = 1,n
       y(j) = sll_f_bspline_1d_eval_deriv( bspline_1d, xx(j))
    end do
    err1 = maxval(abs(y-2*sll_p_pi*cos(2*sll_p_pi*xx)))
    if (err1 > tol) then
       print*,'-------> Test failed in sll_f_bspline_1d_eval_deriv'
       passed_test = .false.
    end if
    print*, " sll_f_bspline_1d_eval_deriv         average error = ", &
         sum(abs(y-2*sll_p_pi*cos(2*sll_p_pi*xx)))/real(n,f64)
    print*, " sll_f_bspline_1d_eval_deriv         maximum error = ", &
         maxval(abs(y-2*sll_p_pi*cos(2*sll_p_pi*xx)))
    call cpu_time(t4)
    !......
    do j = 1,nstep
       call sll_s_bspline_1d_eval_array_deriv( bspline_1d, xx, y)
    end do
    err1 = maxval(abs(y-2*sll_p_pi*cos(2*sll_p_pi*xx)))
    if (err1 > tol) then
       print*,'-------> Test failed in sll_s_bspline_1d_eval_array_deriv'
       passed_test = .false.
    end if
    print*, " sll_s_bspline_1d_eval_array_deriv  average error = ", &
         sum(abs(y-2*sll_p_pi*cos(2*sll_p_pi*xx)))/real(n,f64)
    print*, " sll_s_bspline_1d_eval_array_deriv  maximum error = ", &
         maxval(abs(y-2*sll_p_pi*cos(2*sll_p_pi*xx)))
    call cpu_time(t5)

    !......
    SLL_DEALLOCATE_ARRAY(gtau,ierr)
    SLL_DEALLOCATE_ARRAY(htau,ierr)
    call sll_s_bspline_1d_free(bspline_1d)

  end subroutine test_process_1d

end program test_bspline_1d
