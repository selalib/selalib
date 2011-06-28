module sll_splines
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
  use sll_tridiagonal
  implicit none

  ! The sll_spline_1D provides the services:
  ! - sll_initialize_spline_1d( array, 
  !                             n, 
  !                             xmin, 
  !                             xmax, 
  !                             boundary_cond, 
  !                             spline_obj )
  ! - sll_delete_spline_1d(spline_ptr)
  ! - 
  ! The type should be an 'enumeration' especially if we do things with 
  ! Fortran 2003
  type sll_spline_1D
     sll_int32                         :: n_points ! size
     sll_real64                        :: delta    ! discretization step
     sll_real64                        :: rdelta   ! reciprocal of delta
     sll_real64                        :: xmin
     sll_real64                        :: xmax
     sll_int32                         :: bc_type  ! natural, periodic
     sll_real64, dimension(:), pointer :: data     ! data for the spline fit
     sll_real64, dimension(:), pointer :: c        ! the spline coefficients
  end type sll_spline_1D
  
  
contains  ! ****************************************************************

  ! Tentative, this should be done with some kind of enumeration, especially
  ! if we are considering using F2003, the problem with the object macro is 
  ! that to expose this, we need to redefine these macros in the header file.
#define PERIODIC_SPLINE 0
#define HERMITE_SPLINE  1
#define NATURAL_SPLINE  2


  function new_spline_1D( data, num_pts, xmin, xmax, bc_type )
    type(sll_spline_1D), pointer         :: new_spline_1D
    sll_real64, dimension(:), intent(in), target :: data
    sll_int32,  intent(in)               :: num_pts
    sll_real64, intent(in)               :: xmin
    sll_real64, intent(in)               :: xmax
    sll_int32,  intent(in)               :: bc_type
    sll_int32                            :: ierr
    SLL_ALLOCATE( new_spline_1D, ierr )
    new_spline_1D%n_points = num_pts
    new_spline_1D%xmin     = xmin
    new_spline_1D%xmax     = xmax
    new_spline_1D%delta    = (xmax - xmin)/real(num_pts-1,f64)
    new_spline_1D%rdelta   = 1.0_f64/new_spline_1D%delta
    new_spline_1D%data     => data
    SLL_ALLOCATE( new_spline_1D%c(num_pts), ierr )
    call compute_spline( data, num_pts, bc_type, new_spline_1D )
  end function new_spline_1D

  ! - data: the array whose data must be fit with the cubic spline.
  ! - n: length of the data array that must be fit with the spline.
  ! - bc_type: an integer flag describing the type of boundary conditions 
  !   desired.
  ! - spline_obj: the spline object to be initialized.
  ! This version assumes that the data are uniformly spaced.
  !
  ! The idea behind this spline implementation is the factorization of the 
  ! array:
  !
  !          + 4  1              1 +
  !          | 1  4  1             |
  !          |    .  .  .          |
  !          |       .  .  .       |
  !          |          .  .  .    |
  !          |             1  4  1 |
  !          + 1              1  4 +
  !
  ! in the form:
  !
  !            A = L*L^t
  !
  ! where:
  !
  !          + a                 b +
  !          | b  a                |
  !          |    .  .             |
  !   L =    |       .  .          |    (zeros not shown)
  !          |          .  .       |
  !          |             b  a    |
  !          +                b  a +
  !
  ! This factorization is achieved for the values:
  !
  ! a² = (2+sqrt(3))/6, b² = (2-sqrt(3))/6.
  !
  ! The spline coefficients C are thus the solution of:
  !
  !                  L*L^t*C = F.
  !
  ! Hence, we solve first for D in:
  !
  !                  L*D = F
  !
  ! This solution is achieved in a two-step process. First we compute the
  ! first term d1, which is given by:
  !
  ! d1 =
  !
  ! 1           b           b                       i  b  i
  !---(f(x1) - ---f(xN) + (---)²*f(xN-1) + ... + (-1)(---) (f(xN-i+1)-bd(N-i)))
  ! a           a           a                          a
  !
  ! The series converges (since a > b) and we can approximate the series
  ! by a partial sum:
  !
  !       1                            b
  ! d1 = ---(f(x1) + SUM(i=1,M)(-1)^i(---)^i*f(xN-i+1))
  !       a                            a
  !
  ! The rest of the terms can be found by:
  !
  ! d(i) = 1/a*(f(x i) - b*d(i-1))
  !
  ! Once D is known, the same procedure can be used to compute C in
  !
  !                   L^t*C = D
  !
  ! The same process is carried out. First, the last coefficient is 
  ! calculated by:
  !
  !      c(N) = 1/a*(d(N) + SUM(i=1,M) (-1)^i*(b/a)^i*d(i))
  !
  ! And the rest of the terms, starting with C(N-1) and working backwards:
  !
  !  c(i) = 1/a*(d(i) - b*c(i+1))

  subroutine compute_spline( f, np, bc_type, spline )
    sll_real64, dimension(:), intent(in) :: f    ! data to be fit
    sll_int32,  intent(in)            :: np
    sll_int32,  intent(in)            :: bc_type
    type(sll_spline_1D), pointer      :: spline
    sll_real64, dimension(:), pointer :: coeffs
    sll_int32                         :: i
    sll_int32                         :: err
    sll_real64, parameter             :: a=sqrt((2.0_f64+sqrt(3.0_f64))/6.0_f64)
    sll_real64, parameter             :: r_a = 1.0_f64/a
    sll_real64, parameter             :: b=sqrt((2.0_f64-sqrt(3.0_f64))/6.0_f64)
    sll_real64, parameter             :: b_a = b/a
    sll_real64                        :: coeff_tmp   = 1.0_f64
    sll_real64                        :: d1
    ! consider passing this as an argument or as a scratch space contained
    ! in the spline object.
    sll_real64, dimension(:), allocatable :: d
    SLL_ALLOCATE(d(np), err)
    coeffs     => spline%c

    ! Compute d(1):
#define NUM_TERMS 30
    ! we can check later how many points one will not even affect
    ! the result in the accumulator
    d1 =  f(1)
    do i = 0, NUM_TERMS-1  ! if NUM_TERMS == 0, then only f(np) is considered.
       coeff_tmp = coeff_tmp*(-b_a)
       d1 = d1 + coeff_tmp*f(np-i)
    end do
    ! Fill the d array with the intermediate result
    d(1) = d1*r_a
    do i = 2,np
       d(i) = r_a*(f(i) - b*d(i-1))
    end do
    ! Compute the coefficients. Start with first term
    d1        = d(np)
    coeff_tmp = 1.0_f64
    do i = 1, NUM_TERMS
       coeff_tmp = coeff_tmp*(-b_a)
       d1 = d1 + coeff_tmp*d(i)
    end do
    coeffs(np) = d1*r_a
    ! rest of the coefficients:
    do i = np-1, 1, -1
       coeffs(i) = r_a*(d(i) - b*coeffs(i+1))
    end do
    SLL_DEALLOCATE_ARRAY( d, err )
  end subroutine compute_spline

  ! Here we use the cubic B-spline centered at node 'i', supported on
  ! four cells, two on each side of node 'i':
  !
  !                S_(i-1)    S_i    S_(i+1)  S_(i+2)
  !              +        +        +        +        +
  !              |        |        |        |        |
  !
  !                                *
  !                             *     *
  !                           *         *
  !                          *           *
  !                         *             *
  !                       *                 *
  !                   *                         *  
  !      --------*--------+--------+--------+--------*-----
  !             i-2      i-1       i       i+1      i+2
  !
  ! From left to right, the curves that define the above shape piecewise are:
  ! 
  ! S_(i-1)(x)  defined in [x_(i-2), x_(i-1)]:
  !
  ! = - 1/(4h^3)*(x_(i-2)-x)^3
  !
  ! S_i(x)  defined in [x_(i-1), x_i]:
  !
  ! = 1/(4h^3)*[3*(x_(i-1)-x)^3 + 3*h*(x_(i-1)-x)^2 - 3*h^2*(x_(i-1)-x) + h^3]
  !
  ! S_(i+1)(x)  defined in [x_i, x_(i+1)]:
  !
  ! = 1/(4h^3)*[3*(x-x_(i+1))^3 + 3*h*(x-x_(i+1))^2 - 3*h^2*(x-x_(i+1)) + h^3]
  !
  ! S_(i+2)(x)  defined in [x_(i+1), x_(i+2)]:
  !
  ! = - 1/(4h^3)*(x - x_(i+2))^3
  !
  ! Thus any given point will have a value defined by the 4 contributions of
  ! the splines that are supported in the same interval as the point. 
  ! In this implementation, we choose a normalized representation, in which
  ! h = 1 and a given point is represented by the number of the cell it lives
  ! in and the local offset within this cell. The offset is also normalized to
  ! unity.
  ! 
  ! The image of a given point 'x', under a function modeled by the splines
  ! would be:
  !
  ! S(x) = C_(i-1)*S_(i+2) + C_i*S_(i+1) + C_(i+1)*S_i + C_(i+2)*S_(i-1)
  !
  ! In normalized form, this can be written:
  !
  ! S(x) = 
  !
  ! 1/6*C_(i-1)*(1-dx)^3 + 1/6*C_i*[-3*(1-dx)^3 + 3*(1-dx)^2 +3*(1-dx) + 1] +
  !
  ! 1/6*C_(i+1)*[-3*dx^3 + 3*dx^2 + 3*dx + 1] + 1/6*C_(i+2)*dx^3
  !
  ! The above is the actual interpolation. Here we attempt to minimize the 
  ! number of operations by carrying out a few factorizations. Hence, what is
  ! implemented may be something more akin to:
  !
  ! 1/6*[(1-dx)*[(1-dx)*[(1-dx)*[C_(i-1) - 3*C_i] + 3*C_i] + 3*C_i] + C_i +
  ! 
  !     dx*[dx*[dx*[C_(i+2) - 3*C_(i+1)] + 3*C_(i+1)] + 3*C_(i+1)] + C_(i+1)]
  ! 
  ! As a check, when the point has a zero offset (dx = 0), the above formula
  ! reduces correctly to:
  !
  ! 1/6*[C_(i-1) + 4*C_i + C_(i+1)]
  !
  ! Similarly, if the offset is unity (dx = 1), the above formula reduces to
  !
  ! 1/6*[C_i + 4*C_(i+1) + C_(i+2)]
  ! 
  ! as needed.
  function interpolate_value( x, spline )
    sll_real64                        :: interpolate_value
    intrinsic                         :: int, real
    sll_real64, intent(in)            :: x
    type(sll_spline_1D), pointer      :: spline
    sll_real64, dimension(:), pointer :: coeffs
    sll_real64                        :: rh   ! reciprocal of cell spacing
    sll_int32                         :: cell
    sll_real64                        :: dx
    sll_real64                        :: cdx  ! 1-dx
    sll_real64                        :: t0
    sll_real64                        :: t1
    sll_real64                        :: t2
    sll_real64                        :: t3
    sll_real64                        :: t4
    sll_real64                        :: cim1 ! C_(i-1)
    sll_real64                        :: ci   ! C_i
    sll_real64                        :: cip1 ! C_(i+1)
    sll_real64                        :: cip2 ! C_(i+2)
    sll_int32                         :: num_points
    ! FIXME: arg checks here
    num_points = spline%n_points
    rh     = spline%rdelta
    coeffs => spline%c
    ! find the cell and offset for x
    t0   = x*rh
    cell = int(t0) + 1
    dx   = t0 - real(cell-1)
    cdx  = 1.0_f64 - dx
    !  write (*,'(a,i8, a, f20.12)') 'cell = ', cell, ',   dx = ', dx
    if (cell .eq. 1) then
          cim1 = coeffs(num_points)
          ci   = coeffs(1)
          cip1 = coeffs(2)
          cip2 = coeffs(3)
       else if (cell .eq. (num_points-1)) then   ! last cell
          cim1 = coeffs(cell-1)
          ci   = coeffs(cell)
          cip1 = coeffs(1)  ! last and first points are equal for perodic
          cip2 = coeffs(2)       ! wraps around the end of the array
       else
          cim1 = coeffs(cell-1)
          ci   = coeffs(cell)
          cip1 = coeffs(cell+1)
          cip2 = coeffs(cell+2)
       end if
    t1   = 3.0_f64*ci
    t3   = 3.0_f64*cip1
    t2   = cdx*(cdx*(cdx*(cim1 - t1) + t1) + t1) + ci
    t4   =  dx*( dx*( dx*(cip2 - t3) + t3) + t3) + cip1
    interpolate_value = (1.0_f64/6.0_f64)*(t2 + t4)
  end function interpolate_value
end module sll_splines
