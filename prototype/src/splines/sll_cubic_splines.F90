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

!> \file sll_cubic_splines.F90
!> \namespace sll_cubic_splines
!> \brief  
!> provides capabilities for data interpolation with 
!> cubic B-splines and different boundary conditions
!> \details
!> (at the time of this writing: periodic, hermite). The data to be 
!> interpolated is represented by a simple array.  The spline coefficients 
!> and other information are stored in a spline object, which is also used 
!> to interpolate the fitted data.
!> 
module sll_cubic_splines
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
  use sll_tridiagonal  ! Used for 'slow' algorithm implementation
  use sll_boundary_condition_descriptors
  implicit none
  
  type  ::  sll_cubic_spline_1D
     sll_int32                         :: n_points ! size
     sll_real64                        :: delta    ! discretization step
     sll_real64                        :: rdelta   ! reciprocal of delta
     sll_real64                        :: xmin
     sll_real64                        :: xmax
     sll_int32                         :: bc_type  ! periodic, hermite
     ! check: the following pointer is not used!
     sll_real64, dimension(:), pointer :: data     ! data for the spline fit
     sll_real64, dimension(:), pointer :: d        ! scratch space D (L*D = F),
                                                   ! refer to algorithm below.
                                                   ! Size depends on BC's.
     sll_real64, dimension(:), pointer :: coeffs   ! the spline coefficients
     sll_real64                        :: slope_L  ! left slope, for Hermite
     sll_real64                        :: slope_R  ! right slope, for Hermite
     logical                           :: compute_slope_L
     logical                           :: compute_slope_R
     ! Data required for the 'slow' algorithm based on a standard
     ! tridiagonal system solution. Note that we use the same nomenclature
     ! as in the sll_tridiagonal module.
     logical                           :: use_fast_algorithm
     sll_real64, dimension(:), pointer :: a
     sll_real64, dimension(:), pointer :: cts
     sll_int32, dimension(:), pointer  :: ipiv 
     sll_real64, dimension(:), pointer :: f_aux ! Hermite needs extended f array
  end type sll_cubic_spline_1D

  ! Are x1 and x2 the coordinates that we should use? Or are eta1 and eta2
  ! the more logical ones given that the splines live in an uniform mesh??
  ! Need to define this and then be consistent throughout.
  type sll_cubic_spline_2D
     sll_int32                           :: num_pts_x1
     sll_int32                           :: num_pts_x2
     sll_real64                          :: x1_delta
     sll_real64                          :: x1_rdelta
     sll_real64                          :: x2_delta
     sll_real64                          :: x2_rdelta
     sll_real64                          :: x1_min
     sll_real64                          :: x1_max
     sll_real64                          :: x2_min
     sll_real64                          :: x2_max
     sll_int32                           :: x1_bc_type
     sll_int32                           :: x2_bc_type
     ! if data is not used, it should be deleted make a decision...
     sll_real64, dimension(:,:), pointer :: data   ! data for the spline fit
     sll_real64, dimension(:), pointer   :: d1     ! scratch space D (L*D = F),
                                                   ! refer to algorithm below.
                                                   ! Size depends on BC's.
     sll_real64, dimension(:), pointer   :: d2     ! Second scratch space
     sll_real64, dimension(:,:), pointer :: coeffs   ! the spline coefficients
     sll_real64, dimension(:), pointer   :: x1_min_slopes
     sll_real64, dimension(:), pointer   :: x1_max_slopes
     sll_real64, dimension(:), pointer   :: x2_min_slopes
     sll_real64, dimension(:), pointer   :: x2_max_slopes
     sll_real64, dimension(:), pointer   :: x1_min_slopes_coeffs
     sll_real64, dimension(:), pointer   :: x1_max_slopes_coeffs
     sll_real64, dimension(:), pointer   :: x2_min_slopes_coeffs
     sll_real64, dimension(:), pointer   :: x2_max_slopes_coeffs
     logical                             :: compute_slopes_x1_min
     logical                             :: compute_slopes_x1_max
     logical                             :: compute_slopes_x2_min
     logical                             :: compute_slopes_x2_max
  end type sll_cubic_spline_2D


  interface delete
     module procedure delete_spline_1D, delete_spline_2D
  end interface

  ! Some useful macros that should probably be put in a different file to 
  ! make them more widely accessible. Here we use them to compute default
  ! values for the slopes. 
  ! Actually the proper way to factor this is with a function that takes 
  ! an array segment of a certain length and returns the slope. Such
  ! function(s) could be made to work in 1D or 2D by passing the right
  ! array segment. Keep this in mind...

    ! (-25/12, 4, -3, 4/3, -1/4) stencil
#define FORWARD_FD_5PT( f, r_delta ) \
    r_delta*(-(25.0_f64/12.0_f64)*f(1) + 4.0_f64*f(2) -3.0_f64*f(3) + (4.0_f64/3.0_f64)*f(4) - 0.25_f64*f(5))

    ! (1/4, -4/3, 3, -4, 25/12) stencil
#define BACKWARD_FD_5PT( f, r_delta, np ) \
    r_delta*(0.25_f64*f(np-4) - (4.0_f64/3.0_f64)*f(np-3) + 3.0_f64*f(np-2) - 4.0_f64*f(np-1) + (25.0_f64/12.0_f64)*f(np) )

    ! (-3/2, 2, -1/2 ) stencil
#define FORWARD_FD_3PT( f, r_delta ) \
    r_delta*(-1.5_f64*f(1) + 2.0_f64*f(2) - 0.5_f64*f(3))

    ! ( 1/2, -2, 3/2 ) stencil
#define BACKWARD_FD_3PT( f, r_delta, np ) \
    r_delta*(0.5_f64*f(np-2)-2.0_f64*f(np-1) +1.5_f64*f(np))


contains  ! ****************************************************************

  ! The following implementation embodies the algorithm described in
  ! Eric Sonnendrucker's "A possibly faster algorithm for cubic splines on
  ! a uniform grid" (unpublished).
  
  ! The array of spline coefficients has NP+3 elements. The extra elements
  ! at the ends (i.e.: 0, NP+1, NP+2) store coefficients whose values are
  ! determined by the type of boundary condition used. This is invisible
  ! to the user, who should not be concerned with this implementation detail.

  ! The following parameter determines the problem size at which the alternative
  ! spline algorithm is used.
#define NUM_TERMS 27  
  !> create new spline object
  function new_spline_1D( num_points, xmin, xmax, bc_type, sl, sr )
    type(sll_cubic_spline_1D), pointer         :: new_spline_1D
    sll_int32,  intent(in)               :: num_points
    sll_real64, intent(in)               :: xmin
    sll_real64, intent(in)               :: xmax
    sll_int32,  intent(in)               :: bc_type
    sll_real64, intent(in), optional     :: sl
    sll_real64, intent(in), optional     :: sr
    sll_int32                            :: ierr
    sll_int32                            :: i

    SLL_ALLOCATE( new_spline_1D, ierr )
    new_spline_1D%n_points = num_points
    new_spline_1D%xmin     = xmin
    new_spline_1D%xmax     = xmax
    new_spline_1D%delta    = (xmax - xmin)/real((num_points-1),f64)
    new_spline_1D%rdelta   = 1.0_f64/new_spline_1D%delta
    new_spline_1D%bc_type  = bc_type
    if( num_points .lt. NUM_TERMS ) then
       new_spline_1D%use_fast_algorithm = .false.
    else
       new_spline_1D%use_fast_algorithm = .true.
    end if
    if( xmin .gt. xmax ) then
       print *, 'ERROR, new_spline_1D: xmin is greater than xmax, ', &
            'this would cause all sorts of errors.'
       STOP
    end if
    ! Some more general error checking depending on the type of boundary
    ! condition requested.
    select case (bc_type)
    case (SLL_PERIODIC)
       if( present(sl) .or. present(sr) ) then
          print *, 'new_spline_1D(): it is not allowed to specify the ',&
               'end ifin the case of periodic boundary conditions. ', &
               'Exiting program...'
          STOP 'new_spline_1D'
       else
          ! Assign some value, but this value should never be used in the
          ! periodic case anyway.
          new_spline_1D%slope_L = 0.0
          new_spline_1D%slope_R = 0.0
       end if
       if( new_spline_1D%use_fast_algorithm .eqv. .false. ) then
          SLL_ALLOCATE(new_spline_1D%a(3*num_points),ierr)
          SLL_ALLOCATE(new_spline_1D%cts(7*num_points),ierr)
          SLL_ALLOCATE(new_spline_1D%ipiv(num_points),ierr)
          new_spline_1D%f_aux => null() ! not needed in periodic case
          ! Initialize and factorize the tridiagonal system. See detailed
          ! comment below regarding the structure of this matrix.
          new_spline_1D%a(1) = 1.0_f64/6.0_f64
          new_spline_1D%a(2) = 4.0_f64/6.0_f64
          new_spline_1D%a(3) = 1.0_f64/6.0_f64
          new_spline_1D%a(3*num_points-2) = 1.0_f64/6.0_f64
          new_spline_1D%a(3*num_points-1) = 4.0_f64/6.0_f64
          new_spline_1D%a(3*num_points  ) = 1.0_f64/6.0_f64
          do i=1,num_points-2
             new_spline_1D%a(3*i+1) = 1.0_f64/6.0_f64
             new_spline_1D%a(3*i+2) = 4.0_f64/6.0_f64
             new_spline_1D%a(3*i+3) = 1.0_f64/6.0_f64
          end do
          call setup_cyclic_tridiag( &
               new_spline_1D%a, &
               num_points, &
               new_spline_1D%cts, &
               new_spline_1D%ipiv )
       else
          new_spline_1D%a => null()
          new_spline_1D%cts => null()
          new_spline_1D%ipiv => null()
          new_spline_1D%f_aux => null()
       end if
    case (SLL_HERMITE)
       if( present(sl) ) then
          new_spline_1D%slope_L = sl
          new_spline_1D%compute_slope_L = .false.
       else
          new_spline_1D%slope_L = 0.0           ! just a filler value
          new_spline_1D%compute_slope_L = .true.
       end if
       if( present(sr) ) then
          new_spline_1D%slope_R = sr
          new_spline_1D%compute_slope_R = .false.
       else
          new_spline_1D%slope_R = 0.0           ! just a filler value
          new_spline_1D%compute_slope_R = .true.
       end if
       if( new_spline_1D%use_fast_algorithm .eqv. .false. ) then
          SLL_ALLOCATE(new_spline_1D%a(3*(num_points+2)),ierr)
          SLL_ALLOCATE(new_spline_1D%cts(7*(num_points+2)),ierr)
          SLL_ALLOCATE(new_spline_1D%ipiv(num_points+2),ierr)
          SLL_ALLOCATE(new_spline_1D%f_aux(num_points+2),ierr)
          ! Initialize and factorize the tridiagonal system. See detailed
          ! comment below regarding the structure of this matrix. Matrix 'a'
          ! is (np+2)X(np+2)
          new_spline_1D%a(1) = 0.0_f64
          new_spline_1D%a(2) = 4.0_f64/6.0_f64
          new_spline_1D%a(3) = 2.0_f64/6.0_f64
          new_spline_1D%a(3*(num_points+2)-2) = 2.0_f64/6.0_f64
          new_spline_1D%a(3*(num_points+2)-1) = 4.0_f64/6.0_f64
          new_spline_1D%a(3*(num_points+2)  ) = 0.0_f64
          do i=1,num_points
             new_spline_1D%a(3*i+1) = 1.0_f64/6.0_f64
             new_spline_1D%a(3*i+2) = 4.0_f64/6.0_f64
             new_spline_1D%a(3*i+3) = 1.0_f64/6.0_f64
          end do
          call setup_cyclic_tridiag( &
               new_spline_1D%a, &
               num_points+2, &
               new_spline_1D%cts, &
               new_spline_1D%ipiv )
       else
          new_spline_1D%a => null()
          new_spline_1D%cts => null()
          new_spline_1D%ipiv => null()
          new_spline_1D%f_aux => null()
       end if
    case default
       print *, 'ERROR: new_spline_1D(): not recognized boundary condition'
       STOP
    end select
    if( new_spline_1D%use_fast_algorithm .eqv. .true. ) then
       SLL_ALLOCATE( new_spline_1D%d(num_points),   ierr )
    else
       new_spline_1D%d => null()
    end if
    ! note how the indexing of the coefficients array includes the end-
    ! points 0, num_points, num_points+1, num_points+2. These are meant to 
    ! store the boundary condition-specific data. The 'periodic' BC does
    ! not use the num_points+2 point.
    SLL_ALLOCATE( new_spline_1D%coeffs(0:num_points+2), ierr )
  end function new_spline_1D
  
  ! - data: the array whose data must be fit with the cubic spline.
  ! - np: (number of points; length of the data array that must be fit with 
  !   the spline.
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
  ! 1           b           b  2                   i  b  i
  !---(f(x1) - ---f(xN) + (---)*f(xN-1) + ... + (-1)(---) (f(xN-i+1)-bd(N-i)))
  ! a           a           a                         a
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
  !
  ! The algorithm above is not implemented whenever the number of points is
  ! smaller than 28. In such cases we fall back to a more straightforward but
  ! also more costly implementation using a standard tridiagonal system
  ! solver.
  !
  ! In the periodic case, we are solving the problem A*C=F, where
  !
  !                + 4  1              1 +         1
  !                | 1  4  1             |
  !             1  |    .  .  .          |         .
  !       A =  --- |       .  .  .       |         .
  !             6  |          .  .  .    |         .
  !                |             1  4  1 |
  !                + 1              1  4 +      num_points
  !
  ! In the Hermite BC case, the problem is modified as follows:
  !
  !    + 4  2                + + c_0    +   + 6*f_0+2*delta*f_0´          +  0
  !    | 1  4  1             | |        |   | 6*f_1                       |  1
  !    |    .  .  .          | |    .   |   |    .                        |  .
  ! A =|       .  .  .       |*|    .   | = |    .                        |  .
  !    |          .  .  .    | |    .   |   |    .                        |  .
  !    |             1  4  1 | |        |   | 6*f_np                      |  np
  !    +                2  4 + +c_(np+1)+   + 6*f_(np+1)+2*delta*f_(np+1)'+ np+1



  !> compute spline coefficients

  !> compute_spline_1D() computes the spline coefficients using the parameters
  !> initially given to 'spline' and using 'f' as the data.
  subroutine compute_spline_1D( f, spline )
    sll_real64, dimension(:), intent(in) :: f    ! data to be fit
    sll_int32                            :: bc_type
    type(sll_cubic_spline_1D), pointer         :: spline
    ! Note that this function does no error checking and basically
    ! outsources this task to the functions it is wrapping around.
    ! This is so because those functions can be used independently
    ! (if the user wants to avoid the overhead of calling this
    ! wrapper function), so in any case, the error checking of
    ! the arguments will be carried out at least once.
    bc_type = spline%bc_type;
    select case (bc_type)
    case (SLL_PERIODIC)
       call compute_spline_1D_periodic( f, spline )
    case (SLL_HERMITE)
       call compute_spline_1D_hermite( f, spline )
    case default
       print *, 'ERROR: compute_spline_1D(): not recognized boundary condition'
       STOP
    end select
  end subroutine compute_spline_1D


  ! The following auxiliary functions:
  ! compute_spline_1D_periodic_aux() 
  ! compute_spline_1D_hermite_aux()
  ! are the fundamental building blocks. These are meant to do the work
  ! needed to compute the splines. Other functions are essentially 
  ! wrappers around these.  Clients of these routines are responsible for 
  ! all error-checking.
  ! f: data for which the cubic spline fit is desired.
  ! num_pts: number of points that compose the data
  ! d: scratch array, size num_pts, needed to compute the coefficients.
  ! coeffs: output of the computation, size 0:num_pts+2
  subroutine compute_spline_1D_periodic_aux( f, num_pts, d, coeffs )
    sll_real64, dimension(:), pointer :: f
    sll_int32, intent(in)             :: num_pts
    sll_real64, dimension(:), pointer :: d
    sll_real64, dimension(:), pointer :: coeffs
#ifdef STDF95
    sll_real64, parameter             :: a = 0.78867513459481287_f64
    sll_real64, parameter             :: r_a = 1.2679491924311228_f64
    sll_real64, parameter             :: b = 0.21132486540518716_f64
    sll_real64, parameter             :: b_a = 0.26794919243112275_f64
#else
    sll_real64, parameter             :: a=sqrt((2.0_f64+sqrt(3.0_f64))/6.0_f64)
    sll_real64, parameter             :: r_a = 1.0_f64/a
    sll_real64, parameter             :: b=sqrt((2.0_f64-sqrt(3.0_f64))/6.0_f64)
    sll_real64, parameter             :: b_a = b/a
#endif
    sll_real64                        :: coeff_tmp
    sll_real64                        :: d1
    sll_int32                         :: i
    sll_int32                         :: np
    SLL_ASSERT( size(f) .ge. num_pts )
    SLL_ASSERT( size(d) .ge. num_pts )
    SLL_ASSERT( size(coeffs) .ge. num_pts )
    SLL_ASSERT( num_pts .gt. NUM_TERMS)

    np     =  num_pts
    ! Compute d(1):
    d1 =  f(1)
    coeff_tmp = 1.0_f64
    do i = 0, NUM_TERMS-1  ! if NUM_TERMS == 0, only f(nc) is considered.
!    do i = 1, NUM_TERMS 
       coeff_tmp = coeff_tmp*(-b_a)
       d1 = d1 + coeff_tmp*f(np-1-i)
    end do
    ! Fill the d array with the intermediate result
    d(1) = d1*r_a
    do i = 2,np-1
       d(i) = r_a*(f(i) - b*d(i-1))
    end do
    ! Compute the coefficients. Start with first term
    d1        = d(np-1)
    coeff_tmp = 1.0_f64
    do i = 1, NUM_TERMS
       coeff_tmp = coeff_tmp*(-b_a)
       d1 = d1 + coeff_tmp*d(i)
    end do
!    coeffs(np-1) = d1*r_a
    coeffs(np) = d1*r_a
    ! rest of the coefficients:
    do i = np-2, 1, -1
!       coeffs(i) = r_a*(d(i) - b*coeffs(i+1))
       coeffs(i+1) = r_a*(d(i) - b*coeffs(i+2))
    end do
!!$    coeffs(0)    = coeffs(np-1)
!!$    coeffs(np)   = coeffs(1)
!!$    coeffs(np+1) = coeffs(2)
!!$    coeffs(np+2) = coeffs(3)
    coeffs(1)    = coeffs(np)
    coeffs(np+1)   = coeffs(2)
    coeffs(np+2) = coeffs(3)
    coeffs(np+3) = coeffs(4)

  end subroutine compute_spline_1D_periodic_aux

  subroutine compute_spline_1D_hermite_aux( &
    f,       &
    num_pts, &
    d,       &
    slope_l, &
    slope_r, &
    delta,   &
    coeffs )

    sll_real64, dimension(:), pointer :: f
    sll_int32, intent(in)             :: num_pts
    sll_real64, dimension(:), pointer :: d
    sll_real64, intent(in)            :: slope_l
    sll_real64, intent(in)            :: slope_r
    sll_real64, intent(in)            :: delta
    sll_real64, dimension(:), pointer :: coeffs
    sll_int32                         :: i
    sll_int32                         :: np
#ifdef STDF95
    sll_real64, parameter             :: a = 0.78867513459481287_f64
    sll_real64, parameter             :: r_a = 1.2679491924311228_f64
    sll_real64, parameter             :: b = 0.21132486540518716_f64
    sll_real64, parameter             :: b_a = 0.26794919243112275_f64
    sll_real64, parameter             :: ralpha = 1.8612097182041993_f64 
#else
    sll_real64, parameter             :: a=sqrt((2.0_f64+sqrt(3.0_f64))/6.0_f64)
    sll_real64, parameter             :: r_a = 1.0_f64/a
    sll_real64, parameter             :: b=sqrt((2.0_f64-sqrt(3.0_f64))/6.0_f64)
    sll_real64, parameter             :: b_a = b/a
    sll_real64, parameter             :: ralpha = sqrt(6.0_f64/sqrt(3.0_f64))
#endif
    sll_real64                        :: coeff_tmp
    sll_real64                        :: d1
    sll_real64                        :: f1   ! to store modified value of f(1)
    sll_real64                        :: fnp  ! for modified value of f(np)
    np      =  num_pts
    ! For Hermitian boundary conditions with non-zero slope, we can use the
    ! same algorithm than for the zero-slope case, with the difference that
    ! we tweak the values of the source term, f, right at the endpoints.
    !
    ! This can be derived from the formula for S'(x_N), which yields 
    !
    !             S'(x_N) = 1/(2*h)*[-C_(N-1) + C_(N+1)]
    !
    ! For a slope = 0, this reduces to C_(N+1) = C_(N-1). For slope != 0, we
    ! have an extra constant term (2*h*S') that we can absorb with the 
    ! source term.
    f1  = f(1) 
    fnp = f(np) - delta * slope_r / 3.0_f64

    ! Compute d(1):
    d1 =  f1
    coeff_tmp = 1.0_f64

    ! Since we want to consider the case in which there is a given slope
    ! at point 1, we assume that all points to the left of point 1 are
    ! displaced linearly with an offset given by a line with the given slope.
    ! This requires the subtraction of the  2*slope*delta*(i-1) term...
    do i = 2, NUM_TERMS  ! if NUM_TERMS == 2, only f(2) is considered.
       coeff_tmp = coeff_tmp*(-b_a)
       d1 = d1 + coeff_tmp*(f(i)-2.0*slope_l*delta*(i-1)) 
    end do
    ! Fill the d array with the intermediate result
    d(1) = d1*r_a
    do i = 2,np-1
       d(i) = r_a*(f(i) - b*d(i-1))
    end do
    d(np) = ralpha*(0.5_f64*fnp - b*d(np-1))
    ! Compute the coefficients. Start with first term
!    coeffs(np) = ralpha*d(np)
    coeffs(np+1) = ralpha*d(np)
    do i = np-1, 1, -1
!       coeffs(i) = r_a*(d(i) - b*coeffs(i+1))
       coeffs(i+1) = r_a*(d(i) - b*coeffs(i+2))
    end do
!!$    coeffs(0)    = coeffs(2)    - 2.0 * delta * slope_l
!!$    coeffs(np+1) = coeffs(np-1) + 2.0 * delta * slope_r
!!$    coeffs(np+2) = 0.0 !coeffs(np-2)  ! not used
! debugging...
    coeffs(1)    = coeffs(3)  - 2.0 * delta * slope_l
    coeffs(np+2) = coeffs(np) + 2.0 * delta * slope_r
!    coeffs(1)    = coeffs(3)  + 2.0 * delta * slope_l
!    coeffs(np+2) = coeffs(np) - 2.0 * delta * slope_r
    coeffs(np+3) = 0.0 !coeffs(np-2)  ! not used
  end subroutine compute_spline_1D_hermite_aux


  subroutine compute_spline_1D_periodic( f, spline )
    sll_real64, dimension(:), intent(in), target :: f    ! data to be fit
    type(sll_cubic_spline_1D), pointer         :: spline
    sll_real64, dimension(:), pointer    :: coeffs
    sll_int32                         :: np
    sll_real64, dimension(:), pointer :: d
    sll_real64, dimension(:), pointer :: fp

    if( .not. associated(spline) ) then
       ! FIXME: THROW ERROR
       print *, 'ERROR: compute_spline_1D_periodic(): ', &
            'uninitialized spline object passed as argument. Exiting... '
       STOP
    end if
    if( .not. (size(f) .ge. spline%n_points ) ) then
       ! FIXME: THROW ERROR
       print *, 'ERROR: compute_spline_1D_periodic(): '
       write (*,'(a, i8, a, i8)') 'spline object needs data of size >= ', &
            spline%n_points, ' . Passed size: ', size(f)
       STOP
    end if
    np = spline%n_points
    if( spline%use_fast_algorithm .eqv. .false. ) then
       call solve_cyclic_tridiag_double( &
            spline%cts, &
            spline%ipiv, &
            f, &
            spline%n_points, &
            spline%coeffs(1:np) )
       ! and set the periodic BC by setting the coefficient values
       spline%coeffs(0)    = spline%coeffs(np-1)
       spline%coeffs(np+1) = spline%coeffs(2)
       spline%coeffs(np+2) = spline%coeffs(3)
    else
       fp     => f
       d      => spline%d
       coeffs => spline%coeffs(0:np+2)
       ! Remember that now coeffs(1) refers to spline%coeffs(0)
       call compute_spline_1D_periodic_aux( fp, np, d, coeffs )
    end if
  end subroutine compute_spline_1D_periodic

  subroutine compute_spline_1D_hermite( f, spline )
    sll_real64, dimension(:), intent(in), target :: f    ! data to be fit
    type(sll_cubic_spline_1D), pointer         :: spline
    sll_real64, dimension(:), pointer :: coeffs
    sll_int32                         :: np
    sll_real64, dimension(:), pointer :: fp
    sll_real64, dimension(:), pointer :: d
    sll_real64                        :: slope_l
    sll_real64                        :: slope_r
    sll_real64                        :: delta
    sll_real64                        :: r_delta ! reciprocal

    if( .not. associated(spline) ) then
       ! FIXME: THROW ERROR
       print *, 'ERROR: compute_spline_1D_hermite(): ', &
            'uninitialized spline object passed as argument. Exiting... '
       STOP
    end if
    if( .not. (size(f) .ge. spline%n_points ) ) then
       ! FIXME: THROW ERROR
       print *, 'ERROR: compute_spline_1D_hermite(): '
       write (*,'(a, i8, a, i8)') 'spline object needs data of size >= ', &
            spline%n_points, ' . Passed size: ', size(f)
       STOP
    end if
    fp      => f
    np      =  spline%n_points
    d       => spline%d
    coeffs  => spline%coeffs(0:)
    delta   = spline%delta
    r_delta = 1.0_f64/delta

    ! Estimate numerically the values of the slopes based on the given
    ! values of 'f' if the user did not provide values initially.
    if( spline%compute_slope_L .eqv. .true. ) then
       slope_l = FORWARD_FD_5PT( f, r_delta ) 
    else
       slope_l = spline%slope_L
    end if
    if( spline%compute_slope_R .eqv. .true. ) then
       slope_r = BACKWARD_FD_5PT( f, r_delta, np )
    else
       slope_r = spline%slope_R
    end if

    if( spline%use_fast_algorithm .eqv. .false. ) then
       ! load the source term.
       spline%f_aux(2:np+1) = f(1:np)
       ! set the end-points to reflect the slope information
       spline%f_aux(1)    = f(1)  + (1.0_f64/3.0_f64)*delta*slope_l
       spline%f_aux(np+2) = f(np) - (1.0_f64/3.0_f64)*delta*slope_r
       call solve_cyclic_tridiag_double( &
            spline%cts, &
            spline%ipiv, &
            spline%f_aux, &
            np+2, &
            spline%coeffs )
    else
       call compute_spline_1D_hermite_aux( &
            fp,np,d, slope_l,slope_r, delta, coeffs)
    end if
  end subroutine compute_spline_1D_hermite


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
  ! = - 1/(6h^3)*(x_(i-2)-x)^3
  !
  ! S_i(x)  defined in [x_(i-1), x_i]:
  !
  ! = 1/(6h^3)*[3*(x_(i-1)-x)^3 + 3*h*(x_(i-1)-x)^2 - 3*h^2*(x_(i-1)-x) + h^3]
  !
  ! S_(i+1)(x)  defined in [x_i, x_(i+1)]:
  !
  ! = 1/(6h^3)*[3*(x-x_(i+1))^3 + 3*h*(x-x_(i+1))^2 - 3*h^2*(x-x_(i+1)) + h^3]
  !
  ! S_(i+2)(x)  defined in [x_(i+1), x_(i+2)]:
  !
  ! = - 1/(6h^3)*(x - x_(i+2))^3
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
  ! In normalized form (h=1), this can be written:
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

  ! interpolate_value_aux() is meant to be used as a private function. It
  ! aims at providing a reusable code to carry out interpolations in the
  ! other module functions (like higher dimensions). This is why this 
  ! function does no argument checking. Arguments should be checked by 
  ! the caller.
  function interpolate_value_aux( x, xmin, rh, coeffs )
    sll_real64                        :: interpolate_value_aux
    sll_real64, intent(in)            :: x
    sll_real64, intent(in)            :: xmin
    sll_real64, dimension(:), pointer :: coeffs
    sll_real64, intent(in)            :: rh   ! reciprocal of cell spacing
    sll_int32                         :: cell
    sll_real64                        :: dx
    sll_real64                        :: cdx  ! 1-dx
    sll_real64                        :: t0   ! temp/scratch variables ...
    sll_real64                        :: t1
    sll_real64                        :: t2
    sll_real64                        :: t3
    sll_real64                        :: t4
    sll_real64                        :: cim1 ! C_(i-1)
    sll_real64                        :: ci   ! C_i
    sll_real64                        :: cip1 ! C_(i+1)
    sll_real64                        :: cip2 ! C_(i+2)
    ! find the cell and offset for x
    t0        = (x-xmin)*rh
    cell      = int(t0) + 1
    dx        = t0 - real(cell-1)
    cdx       = 1.0_f64 - dx
!!$    cim1      = coeffs(cell-1)
!!$    ci        = coeffs(cell)
!!$    cip1      = coeffs(cell+1)
!!$    cip2      = coeffs(cell+2)
    cim1      = coeffs(cell)
    ci        = coeffs(cell+1)
    cip1      = coeffs(cell+2)
    cip2      = coeffs(cell+3)
    !    print *, 'intepolate_value_aux(): coefficients:'
    !    print *, cim1, ci, cip1, cip2, cdx, dx
    t1        = 3.0_f64*ci
    t3        = 3.0_f64*cip1
    t2        = cdx*(cdx*(cdx*(cim1 - t1) + t1) + t1) + ci
    t4        =  dx*( dx*( dx*(cip2 - t3) + t3) + t3) + cip1
    ! print *, 't2 and t4: ', t2, t4
    interpolate_value_aux = (1.0_f64/6.0_f64)*(t2 + t4)
    !print *, interpolate_value_aux
  end function interpolate_value_aux

  
  !> get spline interpolate at point x
  function interpolate_value( x, spline )
    sll_real64                        :: interpolate_value
    intrinsic                         :: associated, int, real
#ifdef STDF95
    sll_real64               :: x
    type(sll_cubic_spline_1D)      :: spline
#else
    sll_real64, intent(in)            :: x
    type(sll_cubic_spline_1D), pointer      :: spline
#endif
    sll_real64, dimension(:), pointer :: coeffs
    sll_real64                        :: xmin
    sll_real64                        :: rh   ! reciprocal of cell spacing
    
    ! We set these as assertions since we want the flexibility of turning
    ! them off.
    SLL_ASSERT( (x .ge. spline%xmin) .and. (x .le. spline%xmax) )
#ifdef STDF95
#else
    SLL_ASSERT( associated(spline) )
#endif
    xmin = spline%xmin
    rh   = spline%rdelta
    coeffs => spline%coeffs(0:spline%n_points+2)
    interpolate_value = interpolate_value_aux( x, xmin, rh, coeffs )
  end function interpolate_value

!  !> Just a copy of the function interpolate_value but with a different name.
!  !> Need because, i want a function in my interpolator interface called interpolate value
!  !> and i have a polymorphic error.
!  function interpolate_value_1D( x, spline )
!    sll_real64                        :: interpolate_value_1D
!    intrinsic                         :: associated, int, real
!    sll_real64, intent(in)            :: x
!    type(sll_cubic_spline_1D), pointer      :: spline
!    sll_real64, dimension(:), pointer :: coeffs
!    sll_real64                        :: xmin
!    sll_real64                        :: rh   ! reciprocal of cell spacing
!    ! We set these as assertions since we want the flexibility of turning
!    ! them off.
!    SLL_ASSERT( (x .ge. spline%xmin) .and. (x .le. spline%xmax) )
!    SLL_ASSERT( associated(spline) )
!    xmin = spline%xmin
!    rh        = spline%rdelta
!    coeffs => spline%coeffs
!    interpolate_value_1D = interpolate_value_aux( x, xmin, rh, coeffs )
!  end function interpolate_value_1D
  
  !> interpolates the values given as an array of points as input.
  subroutine interpolate_array_values( a_in, a_out, n, spline )
    intrinsic                               :: associated, int, real
    sll_int32, intent(in)                   :: n
    sll_real64, dimension(1:n), intent(in)  :: a_in
    sll_real64, dimension(1:n), intent(out) :: a_out
    type(sll_cubic_spline_1D), pointer            :: spline
    sll_real64, dimension(:), pointer       :: coeffs
    sll_real64                              :: rh   ! reciprocal of cell spacing
    sll_int32                               :: cell
    sll_real64                              :: dx
    sll_real64                              :: cdx  ! 1-dx
    sll_real64                              :: t0   ! temp/scratch variables ...
    sll_real64                              :: t1
    sll_real64                              :: t2
    sll_real64                              :: t3
    sll_real64                              :: t4
    sll_real64                              :: cim1 ! C_(i-1)
    sll_real64                              :: ci   ! C_i
    sll_real64                              :: cip1 ! C_(i+1)
    sll_real64                              :: cip2 ! C_(i+2)
    sll_int32                               :: num_cells
    sll_real64                              :: x
    sll_int32                               :: i
    SLL_ASSERT( associated(spline) )
    ! FIXME: arg checks here
    num_cells = spline%n_points-1
    rh        = spline%rdelta
    coeffs    => spline%coeffs
    ! find the cell and offset for x
    do i=1,n
       x        = a_in(i)
       if(.not.( (x .ge. spline%xmin) .and. (x .le. spline%xmax) ))then
         print*, 'splines', x,  spline%xmin, spline%xmax
      endif
       !print*, 'splines', x,  spline%xmin, spline%xmax
       SLL_ASSERT( (x .ge. spline%xmin) .and. (x .le. spline%xmax) )
       t0       = (x-spline%xmin)*rh
       cell     = int(t0) + 1
       dx       = t0 - real(cell-1)
       cdx      = 1.0_f64 - dx
       !  write (*,'(a,i8, a, f20.12)') 'cell = ', cell, ',   dx = ', dx
       cim1     = coeffs(cell-1)
       ci       = coeffs(cell)
       cip1     = coeffs(cell+1)
       cip2     = coeffs(cell+2)
       t1       = 3.0_f64*ci
       t3       = 3.0_f64*cip1
       t2       = cdx*(cdx*(cdx*(cim1 - t1) + t1) + t1) + ci
       t4       =  dx*( dx*( dx*(cip2 - t3) + t3) + t3) + cip1
       a_out(i) = (1.0_f64/6.0_f64)*(t2 + t4)
       !print*,'interpolate_array_values', i, a_out(i)
    end do
  end subroutine interpolate_array_values

  !> interpolates the values given as a pointer to an array of points.
  ! FIXME: The following function is not in the unit test.
  subroutine interpolate_pointer_values( ptr_in, ptr_out, n, spline )
    intrinsic                               :: associated, int, real
    sll_int32, intent(in)                   :: n
    sll_real64, dimension(:), pointer       :: ptr_in
    sll_real64, dimension(:), pointer       :: ptr_out
    type(sll_cubic_spline_1D), pointer            :: spline
    sll_real64, dimension(:), pointer       :: coeffs
    sll_real64                              :: rh   ! reciprocal of cell spacing
    sll_int32                               :: cell
    sll_real64                              :: dx
    sll_real64                              :: cdx  ! 1-dx
    sll_real64                              :: t0   ! temp/scratch variables ...
    sll_real64                              :: t1
    sll_real64                              :: t2
    sll_real64                              :: t3
    sll_real64                              :: t4
    sll_real64                              :: cim1 ! C_(i-1)
    sll_real64                              :: ci   ! C_i
    sll_real64                              :: cip1 ! C_(i+1)
    sll_real64                              :: cip2 ! C_(i+2)
    sll_int32                               :: num_cells
    sll_real64                              :: x
    sll_int32                               :: i
    SLL_ASSERT( associated(spline) )
    SLL_ASSERT( associated(ptr_in) )
    SLL_ASSERT( associated(ptr_out) )
    ! FIXME: arg checks here
    num_cells = spline%n_points-1
    rh        = spline%rdelta
    coeffs    => spline%coeffs
    ! find the cell and offset for x
    do i=1,n
       x        = ptr_in(i)
       !print*, 'splines', x,  spline%xmin, spline%xmax
       SLL_ASSERT( (x .ge. spline%xmin) .and. (x .le. spline%xmax) )
       t0       = (x-spline%xmin)*rh
       cell     = int(t0) + 1
       dx       = t0 - real(cell-1)
       cdx      = 1.0_f64 - dx
       !  write (*,'(a,i8, a, f20.12)') 'cell = ', cell, ',   dx = ', dx
       cim1     = coeffs(cell-1)
       ci       = coeffs(cell)
       cip1     = coeffs(cell+1)
       cip2     = coeffs(cell+2)
       t1       = 3.0_f64*ci
       t3       = 3.0_f64*cip1
       t2       = cdx*(cdx*(cdx*(cim1 - t1) + t1) + t1) + ci
       t4       =  dx*( dx*( dx*(cip2 - t3) + t3) + t3) + cip1
       ptr_out(i) = (1.0_f64/6.0_f64)*(t2 + t4)
       !print*,'interpolate_array_values', i, a_out(i)
    end do
  end subroutine interpolate_pointer_values

  ! interpolate_derivative_aux() is a private function aimed at abstracting
  ! away the capability of computing the derivative at a point, given the
  ! array of cubic spline coefficients.

  ! seems that I need to change this interface to include the number of points
  function interpolate_derivative_aux( x, xmin, rh, coeffs )
    sll_real64                        :: interpolate_derivative_aux
    intrinsic                         :: int, real
    sll_real64, intent(in)            :: x
    sll_real64, intent(in)            :: xmin
    sll_real64, intent(in)            :: rh   ! reciprocal of cell spacing
    sll_real64, dimension(:), pointer :: coeffs
    sll_int32                         :: cell
    sll_real64                        :: dx
    sll_real64                        :: t0   ! temp/scratch variables ...
    sll_real64                        :: t1
    sll_real64                        :: t2
    sll_real64                        :: t3
    sll_real64                        :: cim1 ! C_(i-1)
    sll_real64                        :: ci   ! C_i
    sll_real64                        :: cip1 ! C_(i+1)
    sll_real64                        :: cip2 ! C_(i+2)

    ! find the cell and offset for x
    t0        = (x-xmin)*rh
    cell      = int(t0) + 1
    dx        = t0 - real(cell-1)
    ! write (*,'(a,i8, a, f20.12)') 'cell = ', cell, ',   dx = ', dx
!!$    cim1      = coeffs(cell-1)
!!$    ci        = coeffs(cell)
!!$    cip1      = coeffs(cell+1)
!!$    cip2      = coeffs(cell+2)
    cim1      = coeffs(cell)
    ci        = coeffs(cell+1)
    cip1      = coeffs(cell+2)
    cip2      = coeffs(cell+3)
    t1 = 2.0_f64*(cim1 - 2.0_f64*ci + cip1)
    t2 = -cim1 + 3.0_f64*(ci - cip1) + cip2
    t3 =  cip1 - cim1
    interpolate_derivative_aux = 0.5_f64*rh*(dx*(t1 + dx*t2) + t3)
  end function interpolate_derivative_aux


  function interpolate_derivative( x, spline )
    sll_real64                        :: interpolate_derivative
    intrinsic                         :: associated
    sll_real64, intent(in)            :: x
    sll_real64, dimension(:), pointer :: coeffs
#ifdef STDF95
    type(sll_cubic_spline_1D)               :: spline
#else
    type(sll_cubic_spline_1D), pointer      :: spline
#endif

    ! We set these as assertions since we want the flexibility of turning
    ! them off.
    SLL_ASSERT( (x .ge. spline%xmin) .and. (x .le. spline%xmax) )
#ifdef STDF95
#else
    SLL_ASSERT( associated(spline) )
#endif
    coeffs => spline%coeffs(0:spline%n_points+2)
    interpolate_derivative = interpolate_derivative_aux( &
         x, &
         spline%xmin, &
         spline%rdelta, &
         coeffs)
  end function interpolate_derivative

  subroutine interpolate_array_derivatives( &
    array_in, &
    num_pts, &
    array_out, &
    spline )

    intrinsic :: associated
    sll_real64, dimension(:), intent(in)  :: array_in
    sll_int32, intent(in)                 :: num_pts
    sll_real64, dimension(:), intent(out) :: array_out
    type(sll_cubic_spline_1d), pointer          :: spline
    sll_real64, dimension(:), pointer     :: coeffs
    sll_int32 :: i

    SLL_ASSERT( num_pts .le. size(array_in) )
    SLL_ASSERT( associated(spline) )

    coeffs => spline%coeffs(0:spline%n_points+2)
    do i=1,num_pts
       SLL_ASSERT((array_in(i).ge.spline%xmin).and.(array_in(i).le.spline%xmax))
       array_out(i) = interpolate_derivative_aux( &
            array_in(i), spline%xmin, spline%rdelta, coeffs )
    end do
  end subroutine interpolate_array_derivatives

  ! FIXME: The following subroutine is not in the unit test
  subroutine interpolate_pointer_derivatives( &
    ptr_in, &
    num_pts, &
    ptr_out, &
    spline )

    intrinsic :: associated
    sll_real64, dimension(:), pointer  :: ptr_in
    sll_int32, intent(in)              :: num_pts
    sll_real64, dimension(:), pointer  :: ptr_out
    type(sll_cubic_spline_1d), pointer       :: spline
    sll_real64, dimension(:), pointer  :: coeffs
    sll_int32 :: i

    SLL_ASSERT( num_pts .le. size(ptr_in) )
    SLL_ASSERT( associated(spline) )
    SLL_ASSERT( associated(ptr_in) )
    SLL_ASSERT( associated(ptr_out))
    coeffs => spline%coeffs(0:spline%n_points+2)

    do i=1,num_pts
       SLL_ASSERT((ptr_in(i).ge.spline%xmin) .and. (ptr_in(i).le.spline%xmax))
       ptr_out(i) = interpolate_derivative_aux( &
            ptr_in(i), spline%xmin, spline%rdelta, coeffs )
    end do
  end subroutine interpolate_pointer_derivatives


  subroutine delete_spline_1D( spline )
    type(sll_cubic_spline_1D), pointer :: spline
    sll_int32                    :: ierr
    ! Fixme: some error checking, whether the spline pointer is associated
    ! for instance
    SLL_ASSERT( associated(spline) )
    if( spline%use_fast_algorithm .eqv. .true. ) then
       SLL_DEALLOCATE( spline%d, ierr )
    end if
    SLL_DEALLOCATE( spline%coeffs, ierr )
    spline%data => null()
    if( spline%use_fast_algorithm .eqv. .false. ) then
       SLL_DEALLOCATE( spline%a, ierr )
       SLL_DEALLOCATE( spline%cts, ierr )
       SLL_DEALLOCATE( spline%ipiv, ierr )
       if( spline%bc_type == SLL_HERMITE ) then
          SLL_DEALLOCATE( spline%f_aux, ierr )
       end if
    end if
    SLL_DEALLOCATE( spline, ierr )
  end subroutine delete_spline_1D

  !-----------------------------------------------------------------------
  !
  ! Functions and subroutines for the 2D spline.
  !
  !----------------------------------------------------------------------

  function new_spline_2D( &
    num_pts_x1,   &
    num_pts_x2,   &
    x1_min,       &
    x1_max,       &
    x2_min,       &
    x2_max,       &
    x1_bc_type,   &
    x2_bc_type,   &
    const_slope_x1_min, &
    const_slope_x1_max, &
    const_slope_x2_min, &
    const_slope_x2_max, &
    x1_min_slopes, &
    x1_max_slopes, &
    x2_min_slopes, &
    x2_max_slopes )

    type(sll_cubic_spline_2D), pointer         :: new_spline_2D
    sll_int32,  intent(in)               :: num_pts_x1
    sll_int32,  intent(in)               :: num_pts_x2
    sll_real64, intent(in)               :: x1_min
    sll_real64, intent(in)               :: x1_max
    sll_real64, intent(in)               :: x2_min
    sll_real64, intent(in)               :: x2_max
    sll_int32,  intent(in)               :: x1_bc_type
    sll_int32,  intent(in)               :: x2_bc_type
    sll_real64, intent(in), optional     :: const_slope_x1_min
    sll_real64, intent(in), optional     :: const_slope_x1_max
    sll_real64, intent(in), optional     :: const_slope_x2_min
    sll_real64, intent(in), optional     :: const_slope_x2_max
    sll_real64, intent(in), dimension(:), optional :: x1_min_slopes
    sll_real64, intent(in), dimension(:), optional :: x1_max_slopes
    sll_real64, intent(in), dimension(:), optional :: x2_min_slopes
    sll_real64, intent(in), dimension(:), optional :: x2_max_slopes
    sll_int32                            :: bc_selector
    sll_int32                            :: ierr
    sll_int32                            :: i
    SLL_ALLOCATE( new_spline_2D, ierr )
    new_spline_2D%num_pts_x1 = num_pts_x1
    new_spline_2D%num_pts_x2 = num_pts_x2
    new_spline_2D%x1_min     = x1_min
    new_spline_2D%x1_max     = x1_max
    new_spline_2D%x2_min     = x2_min
    new_spline_2D%x2_max     = x2_max
    new_spline_2D%x1_delta   = (x1_max - x1_min)/real((num_pts_x1-1),f64)
    new_spline_2D%x2_delta   = (x2_max - x2_min)/real((num_pts_x2-1),f64)
    new_spline_2D%x1_rdelta  = 1.0_f64/new_spline_2D%x1_delta
    new_spline_2D%x2_rdelta  = 1.0_f64/new_spline_2D%x2_delta
    new_spline_2D%x1_bc_type = x1_bc_type
    new_spline_2D%x2_bc_type = x2_bc_type
    if( (num_pts_x1 .le. NUM_TERMS) .or. (num_pts_x2 .le. NUM_TERMS) ) then
       print *, 'ERROR, new_spline_2D: Because of the algorithm used, this ', &
       'function is meant to be used with arrays that are at least of size = 28'
       STOP 'new_spline_2D()'
    end if
    if( (x1_min .gt. x1_max) .or. (x2_min .gt. x2_max) ) then
       print *, 'ERROR, new_spline_1D: one of the xmin is greater than the ', &
       'corresponding xmax, this would cause all sorts of errors.'
       STOP
    end if

    ! Check that slope arrays are of the right size. Consider making this 
    ! something more permanent than an assertion. Does this compile if the
    ! assertions are turned off???
    if( present(x1_min_slopes) ) then
       SLL_ASSERT(size(x1_min_slopes) .ge. num_pts_x2 )
    end if
    if( present(x1_max_slopes) ) then
       SLL_ASSERT(size(x1_max_slopes) .ge. num_pts_x2 )
    end if
    if( present(x2_min_slopes) ) then
       SLL_ASSERT(size(x2_min_slopes) .ge. num_pts_x1 )
    end if
    if( present(x2_max_slopes) ) then
       SLL_ASSERT(size(x2_max_slopes) .ge. num_pts_x1 )
    end if

    SLL_ALLOCATE( new_spline_2D%d1(num_pts_x1),   ierr )
    SLL_ALLOCATE( new_spline_2D%d2(num_pts_x2),   ierr )

    ! Treat the bc_selector variable essentially like a bit field, to 
    ! accumulate the information on the different boundary conditions
    ! given. This scheme allows to add more types of boundary conditions
    ! if necessary.
    bc_selector = 0
    if( x1_bc_type .eq. SLL_PERIODIC ) then 
       bc_selector = bc_selector + 1
    end if
    if( x1_bc_type .eq. SLL_HERMITE ) then
       bc_selector = bc_selector + 2
    end if
    if( x2_bc_type .eq. SLL_PERIODIC ) then 
       bc_selector = bc_selector + 4
    end if
    if( x2_bc_type .eq. SLL_HERMITE ) then
       bc_selector = bc_selector + 8
    end if
    select case (bc_selector)
    case ( 5 ) 
       ! both boundary condition types are periodic
       if( &
          present(x1_min_slopes) .or. present(x1_max_slopes) .or. &
          present(x2_min_slopes) .or. present(x2_max_slopes) .or. &
          present(const_slope_x1_min) .or. present(const_slope_x1_max) .or. &
          present(const_slope_x2_min) .or. present(const_slope_x2_max) )then

          print *, 'new_spline_2D(): it is not allowed to specify the end', &
               'slopes in the case of doubly periodic boundary conditions.', &
               'Exiting program...'
          STOP 'new_spline_2D'
       else
          new_spline_2d%x1_min_slopes => null()
          new_spline_2d%x1_max_slopes => null()
          new_spline_2d%x2_min_slopes => null()
          new_spline_2d%x2_max_slopes => null()
          new_spline_2d%x1_min_slopes_coeffs => null()
          new_spline_2d%x1_max_slopes_coeffs => null()
          new_spline_2d%x2_min_slopes_coeffs => null()
          new_spline_2d%x2_max_slopes_coeffs => null()
       end if
    case ( 6 ) 
       ! Hermite condition in X1 and periodic in X2 
       if( &
          present(x2_min_slopes) .or. present(x2_max_slopes) .or. &
          present(const_slope_x2_min) .or. present(const_slope_x2_max) ) then
          print *, 'new_spline_2D(): hermite-periodic case, it is not ', &
               'allowed to specify the end slopes in the case of periodic ', &
               'boundary conditions.', &
               'Exiting program...'
          STOP 'new_spline_2D'
       end if
       if( present(const_slope_x1_min) .and. present(x1_min_slopes) ) then
          print *, 'new_spline_2D(): hermite-periodic-case, it is not ', &
               'allowed to specify simultaneously a constant value for ', &
               'the slopes at x1_min and an array-specified set of slopes.'
          STOP 'new_spline_2D'
       end if
       if( present(const_slope_x1_max) .and. present(x1_max_slopes) ) then
          print *, 'new_spline_2D(): hermite-periodic-case, it is not ', &
               'allowed to specify simultaneously a constant value for ', &
               'the slopes at x1_max and an array-specified set of slopes.'
          STOP 'new_spline_2D'
       end if

       ! X2 slope arrays are not needed
       new_spline_2d%x2_min_slopes => null()
       new_spline_2d%x2_max_slopes => null()
       new_spline_2d%x2_min_slopes_coeffs => null()
       new_spline_2d%x2_max_slopes_coeffs => null()

       ! But X1 slopes are.
       SLL_ALLOCATE(new_spline_2d%x1_min_slopes(num_pts_x2),ierr)
       SLL_ALLOCATE(new_spline_2d%x1_max_slopes(num_pts_x2),ierr)

       ! And since the X1 direction splines are computed second, then we
       ! need to convert the slope information into spline coefficient
       ! information.
       SLL_ALLOCATE(new_spline_2d%x1_min_slopes_coeffs(0:num_pts_x2+2),ierr)
       SLL_ALLOCATE(new_spline_2d%x1_max_slopes_coeffs(0:num_pts_x2+2),ierr)

       ! NOTE: because we are using the spline coefficients directly and not
       ! the slope values, the slope values are redundant and at this point
       ! we could deallocate those arrays...

       ! The following macro is obviously intended only for use within
       ! new_spline_2D(). But this should be replaced with a subroutine
       ! whenever possible.
#define FILL_SLOPES(const_opt, input_opt, numpts, output, slopes) \
       if( present(input_opt) ) then;                           \
          do i=1,numpts;                                        \
             new_spline_2D%output(i) = input_opt(i);            \
          end do;                                               \
          new_spline_2D%slopes = .false.;                       \
       else if( present(const_opt) ) then;                      \
          do i=1,numpts;                                        \
             new_spline_2D%output(i) = const_opt;               \
          end do;                                               \
          new_spline_2D%slopes = .false.;                       \
       else;                                                    \
          new_spline_2D%slopes = .true.;                        \
       end if

       ! Set the values of the slopes at x1_min     
       FILL_SLOPES(const_slope_x1_min,x1_min_slopes,num_pts_x2,x1_min_slopes,compute_slopes_x1_min)

       ! Set the values of the slopes at x1_max
       FILL_SLOPES(const_slope_x1_max,x1_max_slopes,num_pts_x2,x1_max_slopes,compute_slopes_x1_max)

    case( 9 )
       ! Periodic in X1 and Hermite in X2
       if( &
          present(x1_min_slopes) .or. present(x1_max_slopes) .or. &
          present(const_slope_x1_min) .or. present(const_slope_x1_max) ) then
          print *, 'new_spline_2D(): periodic-hermite case, it is not ', &
               'allowed to specify the end slopes in the case of periodic ', &
               'boundary conditions.', &
               'Exiting program...'
          STOP 'new_spline_2D'
       end if
       if( present(const_slope_x2_min) .and. present(x2_min_slopes) ) then
          print *, 'new_spline_2D(): periodic-hermite case, it is not ', &
               'allowed to specify simultaneously a constant value for ', &
               'the slopes at x2_min and an array-specified set of slopes.'
          STOP 'new_spline_2D'
       end if
       if( present(const_slope_x2_max) .and. present(x2_max_slopes) ) then
          print *, 'new_spline_2D(): periodic-hermite case, it is not ', &
               'allowed to specify simultaneously a constant value for ', &
               'the slopes at x2_max and an array-specified set of slopes.'
          STOP 'new_spline_2D'
       end if

       ! X1 slope arrays are not needed
       new_spline_2d%x1_min_slopes => null()
       new_spline_2d%x1_max_slopes => null()
       new_spline_2d%x1_min_slopes_coeffs => null()
       new_spline_2d%x1_max_slopes_coeffs => null()

       ! But X2 slopes are.
       SLL_ALLOCATE(new_spline_2d%x2_min_slopes(num_pts_x1),ierr)
       SLL_ALLOCATE(new_spline_2d%x2_max_slopes(num_pts_x1),ierr)
       ! Except that the slope information is used directly, and not as
       ! spline coefficients data.
       new_spline_2d%x2_min_slopes_coeffs => null()
       new_spline_2d%x2_max_slopes_coeffs => null()

       ! Set the values of the slopes at x2_min     
       FILL_SLOPES(const_slope_x2_min,x2_min_slopes,num_pts_x1,x2_min_slopes,compute_slopes_x2_min)

       ! Set the values of the slopes at x2_max
       FILL_SLOPES(const_slope_x2_max,x2_max_slopes,num_pts_x1,x2_max_slopes,compute_slopes_x2_max)

    case( 10 )
       ! Hermite conditions in both, X1 and X2
       if( present(const_slope_x1_min) .and. present(x1_min_slopes) ) then
          print *, 'new_spline_2D(): hermite-hermite case, it is not ', &
               'allowed to specify simultaneously a constant value for ', &
               'the slopes at x1_min and an array-specified set of slopes.'
          STOP 'new_spline_2D'
       end if
       if( present(const_slope_x1_max) .and. present(x1_max_slopes) ) then
          print *, 'new_spline_2D(): hermite-hermite case, it is not ', &
               'allowed to specify simultaneously a constant value for ', &
               'the slopes at x2_max and an array-specified set of slopes.'
          STOP 'new_spline_2D'
       end if
       if( present(const_slope_x2_min) .and. present(x2_min_slopes) ) then
          print *, 'new_spline_2D(): hermite-hermite case, it is not ', &
               'allowed to specify simultaneously a constant value for ', &
               'the slopes at x2_min and an array-specified set of slopes.'
          STOP 'new_spline_2D'
       end if
       if( present(const_slope_x2_max) .and. present(x2_max_slopes) ) then
          print *, 'new_spline_2D(): hermite-hermite case, it is not ', &
               'allowed to specify simultaneously a constant value for ', &
               'the slopes at x2_max and an array-specified set of slopes.'
          STOP 'new_spline_2D'
       end if

       ! Both, X1 and X2 slope arrays are needed.
       SLL_ALLOCATE(new_spline_2d%x1_min_slopes(num_pts_x2),ierr)
       SLL_ALLOCATE(new_spline_2d%x1_max_slopes(num_pts_x2),ierr)
       SLL_ALLOCATE(new_spline_2d%x2_min_slopes(num_pts_x1),ierr)
       SLL_ALLOCATE(new_spline_2d%x2_max_slopes(num_pts_x1),ierr)
       ! Given the order in which we compute the splines, first along X2 and 
       ! then along X1, the coefficients-as-slopes data would never be used
       ! in the X2 direction. Consider eliminating these pointers from the 
       ! object.
       new_spline_2d%x2_min_slopes_coeffs => null()
       new_spline_2d%x2_max_slopes_coeffs => null()
       SLL_ALLOCATE(new_spline_2d%x1_min_slopes_coeffs(0:num_pts_x2+2),ierr)
       SLL_ALLOCATE(new_spline_2d%x1_max_slopes_coeffs(0:num_pts_x2+2),ierr)

       ! Set the values of the slopes at x1_min     
       FILL_SLOPES(const_slope_x1_min,x1_min_slopes,num_pts_x2,x1_min_slopes,compute_slopes_x1_min)

       ! Set the values of the slopes at x1_max
       FILL_SLOPES(const_slope_x1_max,x1_max_slopes,num_pts_x2,x1_max_slopes,compute_slopes_x1_max)

       ! Set the values of the slopes at x2_min     
       FILL_SLOPES(const_slope_x2_min,x2_min_slopes,num_pts_x1,x2_min_slopes,compute_slopes_x2_min)

       ! Set the values of the slopes at x2_max
       FILL_SLOPES(const_slope_x2_max,x2_max_slopes,num_pts_x1,x2_max_slopes,compute_slopes_x2_max)

    case default
       print *, 'ERROR: new_spline_2D(): ', &
            'did not recognize given boundary conditions.'
       STOP
    end select
    ! Reminder: Fortran arrays are column-major ordered...
    ! Note: The indexing of the coefficients array includes the end-
    ! points 0, num_points, num_points+1, num_points+2. These are meant to 
    ! store the boundary condition-specific data. The 'periodic' BC does
    ! not use the num_points+2 point.
    SLL_ALLOCATE( new_spline_2D%coeffs(0:num_pts_x1+2,0:num_pts_x2+2), ierr )
  end function new_spline_2D

  subroutine compute_spline_2D_prdc_prdc( data, spline )
    sll_real64, dimension(:,:), intent(in), target :: data  ! data to be fit
    type(sll_cubic_spline_2D), pointer         :: spline
    sll_real64, dimension(:), pointer    :: coeffs
    sll_int32                            :: npx1
    sll_int32                            :: npx2
    sll_real64, dimension(:), pointer    :: d1
    sll_real64, dimension(:), pointer    :: d2
    sll_real64, dimension(:), pointer    :: datap ! 1D data slice pointer
    sll_int32                            :: i
    sll_int32                            :: j
    if( .not. associated(spline) ) then
       ! FIXME: THROW ERROR
       print *, 'ERROR: compute_spline_2D_prdc_prdc(): ', &
            'uninitialized spline object passed as argument. Exiting... '
       STOP
    end if
    if( (size(data,1) .lt. spline%num_pts_x1 ) ) then
       ! FIXME: THROW ERROR
       print *, 'ERROR: compute_spline_2D_prdc_prdc(): '
       write (*,'(a, i8, a, i8)') 'spline object needs data of size >= ', &
            spline%num_pts_x1, ' . Passed size: ', size(data,1)
       STOP
    end if
    if( (size(data,2) .lt. spline%num_pts_x2 ) ) then
       ! FIXME: THROW ERROR
       print *, 'ERROR: compute_spline_2D_prdc_prdc(): '
       write (*,'(a, i8, a, i8)') 'spline object needs data of size >= ', &
            spline%num_pts_x2, ' . Passed size: ', size(data,2)
       STOP
    end if
    npx1   =  spline%num_pts_x1
    npx2   =  spline%num_pts_x2
    d1     => spline%d1
    d2     => spline%d2
    ! build splines along the x2 direction. Note: due to Fortran's 
    ! column-major ordering, this uses long strides in memory.
    do i=1,npx1 
       datap  => data(i,1:npx2)
       coeffs => spline%coeffs(i,0:npx2+2)
       call compute_spline_1D_periodic_aux( datap, npx2, spline%d2, coeffs )
    end do
    ! build splines along the x1 direction. Note: due to Fortran's 
    ! column-major ordering, this involves short strides in memory.
    ! Note that the data are the spline coefficients computed in the
    ! previous step, so the array dimensions are slightly bigger than in
    ! the original data.
    do j=0,npx2+2  
       datap  => spline%coeffs(1:npx1,j)
       ! same trick regarding the starting point of this pointer. This is
       ! not good.
       coeffs => spline%coeffs(0:npx1+2,j)
       call compute_spline_1D_periodic_aux( datap, npx1, d1, coeffs )
    end do
  end subroutine compute_spline_2D_prdc_prdc

  subroutine compute_spline_2D_hrmt_prdc( data, spline )
    sll_real64, dimension(:,:), intent(in), target :: data  ! data to be fit
    type(sll_cubic_spline_2D), pointer         :: spline
    sll_real64, dimension(:), pointer    :: coeffs
    sll_int32                            :: npx1
    sll_int32                            :: npx2
    sll_real64, dimension(:), pointer    :: d1
    sll_real64, dimension(:), pointer    :: d2
    sll_real64, dimension(:), pointer    :: datap ! 1D data slice pointer
    sll_real64, dimension(:), pointer    :: coeffs_ptr1
    sll_real64, dimension(:), pointer    :: coeffs_ptr2
    sll_real64                           :: min_slope ! slopes at endpoints
    sll_real64                           :: max_slope
    sll_int32                            :: i
    sll_int32                            :: j
    sll_real64                           :: r_x1_delta

    if( .not. associated(spline) ) then
       ! FIXME: THROW ERROR
       print *, 'ERROR: compute_spline_2D_prdc_prdc(): ', &
            'uninitialized spline object passed as argument. Exiting... '
       STOP
    end if
    if( (size(data,1) .lt. spline%num_pts_x1 ) ) then
       ! FIXME: THROW ERROR
       print *, 'ERROR: compute_spline_2D_prdc_prdc(): '
       write (*,'(a, i8, a, i8)') 'spline object needs data of size >= ', &
            spline%num_pts_x1, ' . Passed size: ', size(data,1)
       STOP
    end if
    if( (size(data,2) .lt. spline%num_pts_x2 ) ) then
       ! FIXME: THROW ERROR
       print *, 'ERROR: compute_spline_2D_prdc_prdc(): '
       write (*,'(a, i8, a, i8)') 'spline object needs data of size >= ', &
            spline%num_pts_x2, ' . Passed size: ', size(data,2)
       STOP
    end if
    npx1   =  spline%num_pts_x1
    npx2   =  spline%num_pts_x2
    d1     => spline%d1
    d2     => spline%d2
    ! build splines along the x2 direction (the periodic direction). 
    ! Note: due to Fortran's column-major ordering, this uses long strides 
    ! in memory.
    do i=1,npx1 
       datap  => data(i,1:npx2)
       ! Intentionally, we make coeffs point to the index 1 of the coefficients
       ! array. coeffs(0) is still a valid call, which happens inside
       ! compute_spline_1D_periodic_aux(). This is an ugly trick and
       ! this demonstrates that the _aux() function is broken. (We need 
       ! knowledge of its internals to use it properly). This should be 
       ! fixed.
       coeffs    => spline%coeffs(i,0:npx2+2)
       call compute_spline_1D_periodic_aux( datap, npx2, d2, coeffs )
    end do
    ! build splines along the x1 direction (the Hermite direction). 
    ! Note: due to Fortran's column-major ordering, this involves short 
    ! strides in memory.
    !
    ! Note also that the data are the spline coefficients computed in the
    ! previous step, so the array dimensions are slightly bigger than in
    ! the original data.

    ! First, compute the spline coefficients along the lines (:,0) and
    ! (:,npx2+1)
#if 0
    datap => spline%coeffs(1:npx1,0)
    coeffs => spline%coeffs(0:npx1+2,0)
    min_slope = spline%x1_min_slopes(npx2-1)
    max_slope = spline%x1_max_slopes(npx2-1)
    call compute_spline_1D_hermite_aux( &
         datap, &
         npx1, &
         spline%d1, &
         min_slope, &  
         max_slope,       &
         spline%x1_delta, &
         coeffs )

    datap => spline%coeffs(1:npx1,npx2+1)
    coeffs => spline%coeffs(0:npx1+2,npx2+1)
    min_slope = spline%x1_min_slopes(2)
    max_slope = spline%x1_max_slopes(2)
    call compute_spline_1D_hermite_aux( &
         datap, &
         npx1, &
         spline%d1, &
         min_slope, &  
         max_slope,       &
         spline%x1_delta, &
         coeffs )
#endif
    r_x1_delta = 1.0_f64/spline%x1_delta
    if( spline%compute_slopes_x1_min .eqv. .true. ) then
       ! compute default value for the derivative at the first point based 
       ! on the given data. Can't use the macros here... go figure...
       do j=1,npx2
          spline%x1_min_slopes(j) = &
               r_x1_delta*( -(25.0_f64/12.0_f64)*data(1,j) + &
                                         4.0_f64*data(2,j) - &
                                         3.0_f64*data(3,j) + &
                               (4.0_f64/3.0_f64)*data(4,j) - &
                                        0.25_f64*data(5,j) )
! Earlier version (delete when happy with substitute):
!!$          r_x1_delta*(-1.5_f64*data(1,j) + &
!!$               2.0_f64*data(2,j) - &
!!$               0.5_f64*data(3,j) )
       end do
    end if
    if( spline%compute_slopes_x1_max .eqv. .true. ) then
       ! estimate the derivative at the last point.
       do j=1,npx2
          spline%x1_max_slopes(j) = &
               r_x1_delta*(-0.25_f64*data(npx1-4,j) + &
                   (4.0_f64/3.0_f64)*data(npx1-3,j) + &
                             3.0_f64*data(npx1-2,j) - &
                             4.0_f64*data(npx1-1,j) + &
                 (25.0_f64/12.0_f64)*data(npx1,j) )
!!$! Previous estimate:
!!$          r_x1_delta*( 0.5_f64*data(npx1-2,j) - &
!!$               2.0_f64*data(npx1-1,j) + &
!!$               1.5_f64*data(npx1,j) )
       end do
    end if
    ! At this point, the values of the slopes are available because either
    ! the caller has supplied the values or they have been computed 
    ! numerically. 
    ! Compute the spline coefficients for the available slope values
    coeffs_ptr1 => spline%x1_min_slopes_coeffs(0:npx2+2) 
    call compute_spline_1D_periodic_aux( &
         spline%x1_min_slopes, &
         npx2, &
         spline%d2, &
         coeffs_ptr1 )
    coeffs_ptr2 => spline%x1_max_slopes_coeffs(0:npx2+2)
    call compute_spline_1D_periodic_aux( &
         spline%x1_max_slopes, &
         npx2, &
         spline%d2, &
         coeffs_ptr2 )

    do j=0,npx2+2  
       datap  => spline%coeffs(1:npx1,j)
       coeffs => spline%coeffs(0:npx1+2,j)
       min_slope = spline%x1_min_slopes_coeffs(j)
       max_slope = spline%x1_max_slopes_coeffs(j)
       call compute_spline_1D_hermite_aux( &
            datap, &
            npx1, &
            spline%d1, &
            min_slope, &  
            max_slope,       &
            spline%x1_delta, &
            coeffs )
!       call compute_spline_1D_periodic_aux( datap, npx1, d1, coeffs )
    end do
  end subroutine compute_spline_2D_hrmt_prdc

  subroutine compute_spline_2D_prdc_hrmt( data, spline )
    sll_real64, dimension(:,:), intent(in), target :: data  ! data to be fit
    type(sll_cubic_spline_2D), pointer         :: spline
    sll_real64, dimension(:), pointer    :: coeffs
    sll_int32                            :: npx1
    sll_int32                            :: npx2
    sll_real64, dimension(:), pointer    :: d1
    sll_real64, dimension(:), pointer    :: d2
    sll_real64, dimension(:), pointer    :: datap ! 1D data slice pointer
    sll_real64                           :: min_slope ! slopes at endpoints
    sll_real64                           :: max_slope
    sll_int32                            :: i
    sll_int32                            :: j
    sll_real64                           :: r_x2_delta

    if( .not. associated(spline) ) then
       ! FIXME: THROW ERROR
       print *, 'ERROR: compute_spline_2D_prdc_prdc(): ', &
            'uninitialized spline object passed as argument. Exiting... '
       STOP
    end if
    if( (size(data,1) .lt. spline%num_pts_x1 ) ) then
       ! FIXME: THROW ERROR
       print *, 'ERROR: compute_spline_2D_prdc_prdc(): '
       write (*,'(a, i8, a, i8)') 'spline object needs data of size >= ', &
            spline%num_pts_x1, ' . Passed size: ', size(data,1)
       STOP
    end if
    if( (size(data,2) .lt. spline%num_pts_x2 ) ) then
       ! FIXME: THROW ERROR
       print *, 'ERROR: compute_spline_2D_prdc_prdc(): '
       write (*,'(a, i8, a, i8)') 'spline object needs data of size >= ', &
            spline%num_pts_x2, ' . Passed size: ', size(data,2)
       STOP
    end if
    npx1   =  spline%num_pts_x1
    npx2   =  spline%num_pts_x2
    d1     => spline%d1
    d2     => spline%d2
    r_x2_delta = 1.0_f64/spline%x2_delta

    ! Numerically compute the values of the slopes if the user has not given
    ! them.
    if( spline%compute_slopes_x2_min .eqv. .true. ) then
       ! forward difference scheme to estimate the derivative
       ! at the first point based on the given data.
       do i=1,npx1
          spline%x2_min_slopes(i) = &
               r_x2_delta*(-(25.0_f64/12.0_f64)*data(i,1) + &
               4.0_f64*data(i,2) - & 
               3.0_f64*data(i,3) + &
               (4.0_f64/3.0_f64)*data(i,4) - &
               0.25_f64*data(i,5))

!!$          ! previous estimate:
!!$          r_x2_delta*(-1.5_f64*data(i,1) + &
!!$               2.0_f64*data(i,2) - &
!!$               0.5_f64*data(i,3) )
       end do
    end if
    if( spline%compute_slopes_x2_max .eqv. .true. ) then
       ! backward difference scheme to estimate the derivative
       ! at the last point.
       do i=1,npx1
          spline%x2_max_slopes(i) = &
               r_x2_delta*(-0.25_f64*data(i,npx2-4) + &
                   (4.0_f64/3.0_f64)*data(i,npx2-3) + &
                             3.0_f64*data(i,npx2-2) - &
                             4.0_f64*data(i,npx2-1) + &
                 (25.0_f64/12.0_f64)*data(i,npx2) )

!!$          ! previous estimate (delete eventually)
!!$          r_x2_delta*( 0.5_f64*data(i,npx2-2) - &
!!$               2.0_f64*data(i,npx2-1) + &
!!$               1.5_f64*data(i,npx2) )
       end do
    end if

    ! build splines along the x2 direction (hermite direction). Note: due 
    ! to Fortran's column-major ordering, this uses long strides in memory.
    do i=1,npx1 
       datap  => data(i,1:npx2)
       coeffs => spline%coeffs(i,0:npx2+2)
       min_slope = spline%x2_min_slopes(i)
       max_slope = spline%x2_max_slopes(i)
       call compute_spline_1D_hermite_aux( &
            datap, &
            npx2, &
            spline%d2, &
            min_slope, &
            max_slope, &
            spline%x2_delta, &
            coeffs )
    end do
    ! build splines along the x1 direction. Note: due to Fortran's 
    ! column-major ordering, this involves short strides in memory.
    ! Note that the data are the spline coefficients computed in the
    ! previous step, so the array dimensions are slightly bigger than in
    ! the original data.
    do j=0,npx2+1  
       datap  => spline%coeffs(1:npx1,j)
       coeffs    => spline%coeffs(0:npx1+2,j)
       call compute_spline_1D_periodic_aux( datap, npx1, spline%d1, coeffs )
    end do
  end subroutine compute_spline_2D_prdc_hrmt

  subroutine compute_spline_2D_hrmt_hrmt( data, spline )
    sll_real64, dimension(:,:), intent(in), target :: data  ! data to be fit
    type(sll_cubic_spline_2D), pointer         :: spline
    sll_real64, dimension(:), pointer    :: coeffs
    sll_int32                            :: npx1
    sll_int32                            :: npx2
    sll_real64, dimension(:), pointer    :: d1
    sll_real64, dimension(:), pointer    :: d2
    sll_real64, dimension(:), pointer    :: datap ! 1D data slice pointer
    sll_real64                           :: min_slope ! slopes at endpoints
    sll_real64                           :: max_slope
    sll_int32                            :: i
    sll_int32                            :: j
    sll_real64                           :: r_x1_delta ! reciprocal of x1_delta
    sll_real64                           :: r_x2_delta ! reciprocal of x2_delta

    if( .not. associated(spline) ) then
       ! FIXME: THROW ERROR
       print *, 'ERROR: compute_spline_2D_prdc_prdc(): ', &
            'uninitialized spline object passed as argument. Exiting... '
       STOP
    end if
    if( (size(data,1) .lt. spline%num_pts_x1 ) ) then
       ! FIXME: THROW ERROR
       print *, 'ERROR: compute_spline_2D_prdc_prdc(): '
       write (*,'(a, i8, a, i8)') 'spline object needs data of size >= ', &
            spline%num_pts_x1, ' . Passed size: ', size(data,1)
       STOP
    end if
    if( (size(data,2) .lt. spline%num_pts_x2 ) ) then
       ! FIXME: THROW ERROR
       print *, 'ERROR: compute_spline_2D_prdc_prdc(): '
       write (*,'(a, i8, a, i8)') 'spline object needs data of size >= ', &
            spline%num_pts_x2, ' . Passed size: ', size(data,2)
       STOP
    end if
    npx1   =  spline%num_pts_x1
    npx2   =  spline%num_pts_x2
    d1     => spline%d1
    d2     => spline%d2
    r_x1_delta = 1.0_f64/spline%x1_delta
    r_x2_delta = 1.0_f64/spline%x2_delta

    ! Compute the values of the slopes in case that they were not given.
    if( spline%compute_slopes_x1_min .eqv. .true. ) then
       ! forward difference scheme to estimate the derivative
       ! at the first point based on the given data.
       do j=1,npx2
          spline%x1_min_slopes(j) = &
               r_x1_delta*(-(25.0_f64/12.0_f64)*data(1,j) + &
                                        4.0_f64*data(2,j) - &
                                        3.0_f64*data(3,j) + &
                              (4.0_f64/3.0_f64)*data(4,j) - &
                                       0.25_f64*data(5,j))

!!$          ! previous estimate:
!!$          r_x1_delta*(-1.5_f64*data(1,j) + &
!!$               2.0_f64*data(2,j) - &
!!$               0.5_f64*data(3,j) )
       end do
    end if
    if( spline%compute_slopes_x1_max .eqv. .true. ) then
       ! backward difference scheme to estimate the derivative
       ! at the last point.
       do j=1,npx2
          spline%x1_max_slopes(j) = &
               r_x1_delta*( -0.25_f64*data(npx1-4,j) + &
                    (4.0_f64/3.0_f64)*data(npx1-3,j) + &
                              3.0_f64*data(npx1-2,j) - &
                              4.0_f64*data(npx1-1,j) + &
                  (25.0_f64/12.0_f64)*data(npx1,j) )

!!$          ! previous estimate:
!!$          r_x1_delta*( 0.5_f64*data(npx1-2,j) - &
!!$               2.0_f64*data(npx1-1,j) + &
!!$               1.5_f64*data(npx1,j) )
       end do
    end if

    if( spline%compute_slopes_x2_min .eqv. .true. ) then
       ! forward difference scheme to estimate the derivative
       ! at the first point based on the given data.
       do i=1,npx1
          spline%x2_min_slopes(i) = &
               r_x2_delta*(-(25.0_f64/12.0_f64)*data(i,1) + &
               4.0_f64*data(i,2) - &
               3.0_f64*data(i,3) + &
               (4.0_f64/3.0_f64)*data(i,4) - &
               0.25_f64*data(i,5))

!!$          !previous estimate:
!!$          r_x2_delta*(-1.5_f64*data(i,1) + &
!!$               2.0_f64*data(i,2) - &
!!$               0.5_f64*data(i,3) )
       end do
    end if
    if( spline%compute_slopes_x2_max .eqv. .true. ) then
       ! backward difference scheme to estimate the derivative
       ! at the last point.
       do i=1,npx1
          spline%x2_max_slopes(i) = &
               r_x2_delta*(-0.25_f64*data(i,npx2-4) + &
                   (4.0_f64/3.0_f64)*data(i,npx2-3) + &
                             3.0_f64*data(i,npx2-2) - &
                             4.0_f64*data(i,npx2-1) + &
                 (25.0_f64/12.0_f64)*data(i,npx2) )

!!$          !previous estimate
!!$          r_x2_delta*( 0.5_f64*data(i,npx2-2) - &
!!$               2.0_f64*data(i,npx2-1) + &
!!$               1.5_f64*data(i,npx2) )
       end do
    end if


    ! build splines along the x2 direction. Note: due to Fortran's 
    ! column-major ordering, this uses long strides in memory.
    do i=1,npx1 
       datap  => data(i,1:npx2)
       ! Intentionally, we make coeffs point to the index 1 of the coefficients
       ! array. coeffs(0) is still a valid call, which is what 
       ! compute_spline_1D_periodic_aux() does. This is an ugly trick and
       ! this demonstrates that the _aux() function is broken. (We need 
       ! knowledge of its internals to use it properly). This should be 
       ! fixed.
       coeffs => spline%coeffs(i,0:npx2+2) 
       min_slope = spline%x2_min_slopes(i)
       max_slope = spline%x2_max_slopes(i)
       call compute_spline_1D_hermite_aux( &
            datap, &
            npx2, &
            spline%d2, &
            min_slope, &  
            max_slope, &
            spline%x2_delta, &
            coeffs )
    end do
    ! build splines along the x1 direction. Note: due to Fortran's 
    ! column-major ordering, this involves short strides in memory.
    ! Note that the data are the spline coefficients computed in the
    ! previous step, so the array dimensions are slightly bigger than in
    ! the original data.

    ! First we compute out of the loop the splines for the coefficients 
    ! in the (:,0) and (:,npx2+1) rows. TO DO THIS PROPERLY WE SHOULD
    ! INTRODUCE AN ESTIMATE OF THE SLOPES, WHICH FOR THESE VALUES FALL 
    ! OUT OF RANGE. Here we just "reflect" the values around the edge point.
    ! ... it's better than nothing.
    datap => spline%coeffs(1:npx1,0)
    coeffs => spline%coeffs(0:npx1+2,0)
    min_slope = spline%x1_min_slopes(2)
    max_slope = spline%x1_max_slopes(2)
    call compute_spline_1D_hermite_aux( &
         datap, &
         npx1, &
         spline%d1, &
         min_slope, &  
         max_slope,       &
         spline%x1_delta, &
         coeffs )

    datap => spline%coeffs(1:npx1,npx2+1)
    coeffs => spline%coeffs(0:npx1+2,npx2+1)
    min_slope = spline%x1_min_slopes(npx2-1)
    max_slope = spline%x1_max_slopes(npx2-1)
    call compute_spline_1D_hermite_aux( &
         datap, &
         npx1, &
         spline%d1, &
         min_slope, &  
         max_slope,       &
         spline%x1_delta, &
         coeffs )

    do j=1,npx2  
       datap  => spline%coeffs(1:npx1,j)
       coeffs => spline%coeffs(0:npx1+2,j)
       min_slope = spline%x1_min_slopes(j)
       max_slope = spline%x1_max_slopes(j)
       call compute_spline_1D_hermite_aux( &
            datap, &
            npx1, &
            spline%d1, &
            min_slope, &  
            max_slope, &
            spline%x1_delta, &
            coeffs )
    end do
  end subroutine compute_spline_2D_hrmt_hrmt

  subroutine compute_spline_2D( data, spline )
    sll_real64, dimension(:,:), intent(in), target :: data  ! data to be fit
    type(sll_cubic_spline_2D), pointer         :: spline
    sll_int32 :: bc1
    sll_int32 :: bc2
    sll_int32 :: bc_selector
    if( .not. associated(spline) ) then
       ! FIXME: THROW ERROR
       print *, 'ERROR: compute_spline_2D(): ', &
            'uninitialized spline object passed as argument. Exiting... '
       STOP
    end if
 
    bc1 = spline%x1_bc_type
    bc2 = spline%x2_bc_type

    ! Treat the bc_selector variable essentially like a bit field, to 
    ! accumulate the information on the different boundary conditions
    ! given. This scheme allows to add more types of boundary conditions
    ! if necessary.
    bc_selector = 0

    ! We make every case explicit to facilitate adding more BC types in
    ! the future.
    if( bc1 .eq. SLL_PERIODIC ) then 
       bc_selector = bc_selector + 1
    end if

    if( bc1 .eq. SLL_HERMITE ) then 
       bc_selector = bc_selector + 2
    end if

    if( bc2 .eq. SLL_PERIODIC ) then 
       bc_selector = bc_selector + 4
    end if

    if( bc2 .eq. SLL_HERMITE ) then
       bc_selector = bc_selector + 8
    end if

    select case (bc_selector)
       case ( 5 ) 
          ! both boundary condition types are periodic
          call compute_spline_2D_prdc_prdc( data, spline )
       case ( 6 )
          ! hermite in X1 and periodic in X2
          call compute_spline_2D_hrmt_prdc( data, spline )
       case ( 9 ) 
          ! periodic condition in X1 and hermite in X2 
          call compute_spline_2D_prdc_hrmt( data, spline )
       case( 10 )
          ! Hermite conditions in both, X1 and X2
          call compute_spline_2D_hrmt_hrmt( data, spline )
       case default
          print *, 'ERROR: compute_spline_2D(): ', &
            'did not recognize given boundary condition combination.'
       STOP
    end select
  end subroutine compute_spline_2D

  ! deposit_value_2D(): given a spline that describes the decomposition 
  ! of the distribution function at time t^n, and two 2D arrays x1 and x2 
  ! where the foot of the forward characteristics are stored, returns
  ! a 2D array a_out which is the updated distribution function at time t^{n+1}
  !
  ! the boundary conditions are taken into account and any type of BC are 
  ! allowed
  subroutine deposit_value_2D(x1, x2, spline, a_out)
    intrinsic :: real, int
    sll_real64, dimension(1:,1:), intent(in)      :: x1
    sll_real64, dimension(1:,1:), intent(in)      :: x2
    type(sll_cubic_spline_2D), pointer                  :: spline
    sll_real64, dimension(:,:),intent(out)        :: a_out

    sll_real64  :: cij   ! C_ij
    sll_real64  :: x1_min
    sll_real64  :: x2_min
    sll_real64  :: dx1
    sll_real64  :: dx2
    sll_real64  :: cdx1  ! 1-dx1
    sll_real64  :: cdx2  ! 1-dx2
    sll_int32   :: cell1
    sll_int32   :: cell2
    sll_real64  :: rh1
    sll_real64  :: rh2
    sll_int32   :: n1
    sll_int32   :: n2
    
		! local variables
    sll_int32   :: i1
    sll_int32   :: i2
    
    sll_int32   :: nt1
    sll_int32   :: nt2
		
    sll_real64  :: svalx1, svalx2, svalx3, svalx4
    sll_real64  :: svaly1, svaly2, svaly3, svaly4
    sll_real64  :: ax1, ax2, ax3, ay1, ay2, ay3
		
    sll_real64 :: t1,t2
    sll_int32  :: ipm1,ip,ipp1,ipp2
    sll_int32  :: jpm1,jp,jpp1,jpp2
    
    sll_int32 :: bc1
    sll_int32 :: bc2
    
    if( .not. associated(spline) ) then
       ! FIXME: THROW ERROR
       print *, 'ERROR: deposit_value_2D(): ', &
            'uninitialized spline object passed as argument. Exiting... '
       STOP
    end if
    
    if ((size(x1,1).ne.spline%num_pts_x1).or.(size(x1,2).ne.spline%num_pts_x2)) then
       ! FIXME: THROW ERROR
       print *, 'ERROR: deposit_value_2D(): '
       write (*,'(a, i8, i8, a, i8, i8, a)') 'array of feets of characteristics needs data of size = (', &
            spline%num_pts_x1, spline%num_pts_x2,') . Passed size: (', size(x1,1), size(x1,2),')'
       STOP
    end if
    if ((size(x2,1).ne.spline%num_pts_x1).or.(size(x2,2).ne.spline%num_pts_x2)) then
       ! FIXME: THROW ERROR
       print *, 'ERROR: deposit_value_2D(): '
       write (*,'(a, i8, i8, a, i8, i8, a)') 'array of feets of characteristics needs data of size = (', &
            spline%num_pts_x1, spline%num_pts_x2,') . Passed size: (', size(x2,1), size(x2,2),')'
       STOP
    end if
    if ((size(a_out,1).ne.spline%num_pts_x1).or.(size(a_out,2).ne.spline%num_pts_x2)) then
       ! FIXME: THROW ERROR
       print *, 'ERROR: deposit_value_2D(): '
       write (*,'(a, i8, i8, a, i8, i8, a)') 'array of feets of characteristics needs data of size = (', &
            spline%num_pts_x1, spline%num_pts_x2,') . Passed size: (', size(a_out,1), size(a_out,2),')'
       STOP
    end if

    bc1 = spline%x1_bc_type
    bc2 = spline%x2_bc_type

    x1_min     = spline%x1_min
    x2_min     = spline%x2_min
    rh1        = spline%x1_rdelta
    rh2        = spline%x2_rdelta
				
    n1         = spline%num_pts_x1
    n2         = spline%num_pts_x2
    
    if( bc1 .eq. SLL_PERIODIC ) then 
      nt1 = n1-1
    end if
    if( bc1 .eq. SLL_HERMITE ) then 
      nt1 = n1
    end if
    if( bc2 .eq. SLL_PERIODIC ) then 
      nt2 = n2-1
    end if
    if( bc2 .eq. SLL_HERMITE ) then 
      nt2 = n2
    end if

    a_out = 0._f64
		
    do i1 = 1,nt1                       
      do i2 = 1,nt2
		
        ! find the cell and offset for x1
        t1          = (x1(i1,i2)-x1_min)*rh1
        cell1       = floor(t1) + 1
        dx1         = t1- real(cell1-1)
        cdx1        = 1.0_f64 - dx1
		
        ! find the cell and offset for x2
        t2         = (x2(i1,i2)-x2_min)*rh2
        cell2      = floor(t2) + 1
        dx2        = t2 - real(cell2-1)
        cdx2       = 1.0_f64 - dx2
      
        ! checking if the particle is well localized
        ! the particle is in the cell (cell1,cell2)
        ! placed in ((cell1-1)/rh1+x1_min,(cell2-1)/rh2+x2_min)
        if (((cell1-1)/rh1+x1_min>x1(i1,i2)).or.(x1(i1,i2)>cell1/rh1+x1_min)) then
          print*,'problem with the localization of r', (cell1-1)/rh1+x1_min,x1(i1,i2),cell1/rh1+x1_min,dx1
        end if
        if (((cell2-1)/rh2+x2_min>x2(i1,i2)).or.(x2(i1,i2)>cell2/rh2+x2_min)) then
          print*,'problem with the localization of theta', (cell2-1)/rh2+x2_min,x2(i1,i2),cell2/rh2+x2_min,dx2
        end if

        cij = spline%coeffs(i1,i2)
      
        ! index depending on the BC 
        if( bc1 .eq. SLL_PERIODIC ) then 
          ipm1 = mod(cell1+n1-3,n1-1)+1
          ip   = mod(cell1+n1-2,n1-1)+1
          ipp1 = mod(cell1+n1-1,n1-1)+1
          ipp2 = mod(cell1+n1  ,n1-1)+1
        end if
        if( bc1 .eq. SLL_HERMITE ) then 
          ipm1=cell1-1
          ip  =cell1
          ipp1=cell1+1
          ipp2=cell1+2
        end if
        if( bc2 .eq. SLL_PERIODIC ) then 
          jpm1 = mod(cell2+n2-3,n2-1)+1
          jp   = mod(cell2+n2-2,n2-1)+1
          jpp1 = mod(cell2+n2-1,n2-1)+1
          jpp2 = mod(cell2+n2  ,n2-1)+1
        end if
        if( bc2 .eq. SLL_HERMITE ) then 
          jpm1=cell2-1
          jp  =cell2
          jpp1=cell2+1
          jpp2=cell2+2
        end if
        			
        ax1 = cdx1 
        ax2 = ax1*ax1
        ax3 = ax2*ax1
        
        ay1 = cdx2
        ay2 = ay1*ay1
        ay3 = ay2*ay1
									
        svalx1 = ax3/6._f64
        svalx2 = (0.5_f64*(-ax3+ax2+ax1) + 1._f64/6._f64)
        svalx3 = (2._f64/3._f64 - ax2 + 0.5_f64*ax3)
        svalx4 = ((1._f64-ax3)/6._f64 + 0.5_f64*(ax2-ax1))
		
        svaly1 = ay3/6._f64
        svaly2 = 0.5_f64*(-ay3+ay2+ay1) + 1._f64/6._f64
        svaly3 = 2._f64/3.0_f64 - ay2 + 0.5_f64*ay3
        svaly4 = (1._f64-ay3)/6._f64 + 0.5_f64*(ay2-ay1)
			
        if (ipm1.ge.1) then
          if (jpm1.ge.1) then
            a_out(ipm1,jpm1) = a_out(ipm1,jpm1) + cij*svalx1*svaly1
          end if
          a_out(ipm1,jp)     = a_out(ipm1,jp)   + cij*svalx1*svaly2
          if (jpp1.le.n2) then
            a_out(ipm1,jpp1) = a_out(ipm1,jpp1) + cij*svalx1*svaly3
          end if
          if (jpp2.le.n2) then
            a_out(ipm1,jpp2) = a_out(ipm1,jpp2) + cij*svalx1*svaly4
          end if
        end if
        
        if (jpm1.ge.1) then
          a_out(ip,jpm1) = a_out(ip,jpm1) + cij*svalx2*svaly1
        end if
        a_out(ip,jp)     = a_out(ip,jp)   + cij*svalx2*svaly2
        if (jpp1.le.n2) then
          a_out(ip,jpp1) = a_out(ip,jpp1) + cij*svalx2*svaly3
        end if
        if (jpp2.le.n2) then
          a_out(ip,jpp2) = a_out(ip,jpp2) + cij*svalx2*svaly4
        end if
            
        if (ipp1.le.n1) then
          if (jpm1.ge.1) then
            a_out(ipp1,jpm1) = a_out(ipp1,jpm1) + cij*svalx3*svaly1
          end if
          a_out(ipp1,jp)     = a_out(ipp1,jp)   + cij*svalx3*svaly2
          if (jpp1.le.n2) then
            a_out(ipp1,jpp1) = a_out(ipp1,jpp1) + cij*svalx3*svaly3
          end if
          if (jpp2.le.n2) then
            a_out(ipp1,jpp2) = a_out(ipp1,jpp2) + cij*svalx3*svaly4
          end if
        end if
        
        if (ipp2.le.n1) then
          if (jpm1.ge.1) then
            a_out(ipp2,jpm1) = a_out(ipp2,jpm1) + cij*svalx4*svaly1
          end if
          a_out(ipp2,jp)     = a_out(ipp2,jp)   + cij*svalx4*svaly2
          if (jpp1.le.n2) then
            a_out(ipp2,jpp1) = a_out(ipp2,jpp1) + cij*svalx4*svaly3
          end if
          if (jpp2.le.n2) then
            a_out(ipp2,jpp2) = a_out(ipp2,jpp2) + cij*svalx4*svaly4
          end if
        end if
              
        if (bc1.eq.SLL_HERMITE) then 
           if (i1.eq.1) then
              if (jpm1.ge.1) then
                 a_out(1,jpm1) = a_out(1,jpm1) + spline%coeffs(0,i2)/6._f64*svaly1
              end if
              a_out(1,jp) = a_out(1,jp) + spline%coeffs(0,i2)/6._f64*svaly2
            if (jpp1.le.n2) then
              a_out(1,jpp1) = a_out(1,jpp1) + spline%coeffs(0,i2)/6._f64*svaly3
            end if
            if (jpp2.le.n2) then
              a_out(1,jpp2) = a_out(1,jpp2) + spline%coeffs(0,i2)/6._f64*svaly4
            end if
          end if
      
          if (i1.eq.n1) then
            if (jpm1.ge.1) then
              a_out(n1,jpm1) = a_out(n1,jpm1) + spline%coeffs(n1+1,i2)/6._f64*svaly1
            end if
            a_out(n1,jp)     = a_out(n1,jp)   + spline%coeffs(n1+1,i2)/6._f64*svaly2
            if (jpp1.le.n2) then
              a_out(n1,jpp1) = a_out(n1,jpp1) + spline%coeffs(n1+1,i2)/6._f64*svaly3
            end if
            if (jpp2.le.n2) then
              a_out(n1,jpp2) = a_out(n1,jpp2) + spline%coeffs(n1+1,i2)/6._f64*svaly4
            end if
          end if
        end if
      
        if (bc2.eq.SLL_HERMITE) then 
          if (i2.eq.1) then
            if (ipm1.ge.1) then
              a_out(ipm1,1) = a_out(ipm1,1) + spline%coeffs(i1,0)/6._f64*svalx1
            end if
            a_out(ip,1) = a_out(ip,1) + spline%coeffs(i1,0)/6._f64*svalx2
            if (ipp1.ne.n1) then
              a_out(ipp1,1) = a_out(ipp1,1) + spline%coeffs(i1,0)/6._f64*svalx3
            end if
            if (ipp2.ne.n1) then
              a_out(ipp2,1) = a_out(ipp2,1) + spline%coeffs(i1,0)/6._f64*svalx4
            end if
          end if
      
          if (i2.eq.n2) then
            if (ipm1.ge.1) then
              a_out(ipm1,n2) = a_out(ipm1,n2) + spline%coeffs(i1,n2+1)/6._f64*svalx1
            end if
            a_out(ip,n1)     = a_out(ip,n1)   + spline%coeffs(i1,n2+1)/6._f64*svalx2
            if (ipp1.ne.n1) then
              a_out(ipp1,n2) = a_out(ipp1,n2) + spline%coeffs(i1,n2+1)/6._f64*svalx3
            end if
            if (ipp2.ne.n1) then
              a_out(ipp2,n2) = a_out(ipp2,n2) + spline%coeffs(i1,n2+1)/6._f64*svalx4
            end if
          end if
        end if
        
      end do
    end do
    
    if( bc1 .eq. SLL_PERIODIC ) then 
      a_out(n1,:) = a_out(1,:)
    end if
    if( bc2 .eq. SLL_PERIODIC ) then 
      a_out(:,n2) = a_out(:,1)
    end if
    			 
  end subroutine deposit_value_2D


  function interpolate_value_2D( x1, x2, spline )
    sll_real64                          :: interpolate_value_2D
    intrinsic                           :: associated, int, real
    sll_real64, intent(in)              :: x1
    sll_real64, intent(in)              :: x2
    type(sll_cubic_spline_2D), pointer        :: spline
    sll_real64                          :: rh1   ! reciprocal of cell spacing
    sll_real64                          :: rh2   ! reciprocal of cell spacing
    sll_int32                           :: cell
    sll_real64                          :: dx
    sll_real64                          :: cdx  ! 1-dx
    sll_real64                          :: t0   ! temp/scratch variables ...
    sll_real64                          :: t1
    sll_real64                          :: t2
    sll_real64                          :: t3
    sll_real64                          :: t4
    sll_real64                          :: cim1 ! C_(i-1)
    sll_real64                          :: ci   ! C_i
    sll_real64                          :: cip1 ! C_(i+1)
    sll_real64                          :: cip2 ! C_(i+2)
    sll_real64                          :: x1_min
    sll_real64                          :: x2_min
    sll_int32                           :: num_pts_x1
    sll_int32                           :: num_pts_x2
    sll_real64, dimension(:), pointer   :: coeffs_line_jm1
    sll_real64, dimension(:), pointer   :: coeffs_line_j
    sll_real64, dimension(:), pointer   :: coeffs_line_jp1
    sll_real64, dimension(:), pointer   :: coeffs_line_jp2
    ! We set these as assertions since we want the flexibility of turning
    ! them off.
    SLL_ASSERT( (x1 .ge. spline%x1_min) .and. (x1 .le. spline%x1_max) )
    SLL_ASSERT( (x2 .ge. spline%x2_min) .and. (x2 .le. spline%x2_max) )
    SLL_ASSERT( associated(spline) )
    x1_min     = spline%x1_min
    x2_min     = spline%x2_min
    num_pts_x1 = spline%num_pts_x1
    num_pts_x2 = spline%num_pts_x2
    rh1        = spline%x1_rdelta
    rh2        = spline%x2_rdelta
    ! find the cell and offset for x2
    t0         = (x2-x2_min)*rh2
    cell       = int(t0) + 1
    dx         = t0 - real(cell-1)
    cdx        = 1.0_f64 - dx
    !  write (*,'(a,i8, a, f20.12)') 'cell = ', cell, ',   dx = ', dx
    ! interpolate the coefficients along the line of constant x1. These 
    ! computations are independent from one another. A little problem is
    ! the redundancy in the computation of the cell and offset along each
    ! of the constant x2 lines, as this will be done by each call of 
    ! interpolate_value_aux(). This suggests that the proper refactoring
    ! of this function would have the cell and offset as arguments.
    coeffs_line_jm1 => spline%coeffs(0:num_pts_x1+2, cell-1)
    coeffs_line_j   => spline%coeffs(0:num_pts_x1+2, cell)
    coeffs_line_jp1 => spline%coeffs(0:num_pts_x1+2, cell+1)
    coeffs_line_jp2 => spline%coeffs(0:num_pts_x1+2, cell+2)
    cim1      = interpolate_value_aux(x1, x1_min, rh1, coeffs_line_jm1)
    ci        = interpolate_value_aux(x1, x1_min, rh1, coeffs_line_j  )
    cip1      = interpolate_value_aux(x1, x1_min, rh1, coeffs_line_jp1)
    cip2      = interpolate_value_aux(x1, x1_min, rh1, coeffs_line_jp2)
    t1        = 3.0_f64*ci
    t3        = 3.0_f64*cip1
    t2        = cdx*(cdx*(cdx*(cim1 - t1) + t1) + t1) + ci
    t4        =  dx*( dx*( dx*(cip2 - t3) + t3) + t3) + cip1
    interpolate_value_2D = (1.0_f64/6.0_f64)*(t2 + t4)
  end function interpolate_value_2D

  ! interpolate_x1_derivative_2D(): given discrete data f(i,j) that are
  ! described by a 2-dimensional cubic spline fit s(x1,x2), where the
  ! continuous variables x1 and x2 are within the original limits of i and j
  ! respectively, interpolate_x1_derivative() returns the value of
  !
  !         partial s
  !       -------------
  !         partial x1
  !
  ! evaluated at the point (x1,x2). (Sorry for the ambiguous use of x1)
  function interpolate_x1_derivative_2D( x1, x2, spline )
    sll_real64                          :: interpolate_x1_derivative_2D
    intrinsic                           :: associated, int, real
    sll_real64, intent(in)              :: x1
    sll_real64, intent(in)              :: x2
#ifdef STDF95
    type(sll_cubic_spline_2D)                 :: spline
#else
    type(sll_cubic_spline_2D), pointer        :: spline
#endif
    sll_real64                          :: rh1   ! reciprocal of cell spacing
    sll_real64                          :: rh2   ! reciprocal of cell spacing
    sll_int32                           :: cell
    sll_real64                          :: dx
    sll_real64                          :: cdx  ! 1-dx
    sll_real64                          :: t0   ! temp/scratch variables ...
    sll_real64                          :: t1
    sll_real64                          :: t2
    sll_real64                          :: t3
    sll_real64                          :: t4
    sll_real64                          :: cim1 ! C_(i-1)
    sll_real64                          :: ci   ! C_i
    sll_real64                          :: cip1 ! C_(i+1)
    sll_real64                          :: cip2 ! C_(i+2)
    sll_real64                          :: x1_min
    sll_real64                          :: x2_min
    sll_int32                           :: num_pts_x1
    sll_int32                           :: num_pts_x2
    sll_real64, dimension(:), pointer   :: coeffs_line_jm1
    sll_real64, dimension(:), pointer   :: coeffs_line_j
    sll_real64, dimension(:), pointer   :: coeffs_line_jp1
    sll_real64, dimension(:), pointer   :: coeffs_line_jp2
    ! We set these as assertions since we want the flexibility of turning
    ! them off.
    SLL_ASSERT( (x1 .ge. spline%x1_min) .and. (x1 .le. spline%x1_max) )
    SLL_ASSERT( (x2 .ge. spline%x2_min) .and. (x2 .le. spline%x2_max) )
#ifdef STDF95
#else
    SLL_ASSERT( associated(spline) )
#endif
    x1_min     = spline%x1_min
    x2_min     = spline%x2_min
    num_pts_x1 = spline%num_pts_x1
    num_pts_x2 = spline%num_pts_x2
    rh1        = spline%x1_rdelta
    rh2        = spline%x2_rdelta
    ! find the cell and offset for x2
    t0         = (x2-x2_min)*rh2
    cell       = int(t0) + 1
    dx         = t0 - real(cell-1)
    cdx        = 1.0_f64 - dx
    !  write (*,'(a,i8, a, f20.12)') 'cell = ', cell, ',   dx = ', dx
    ! interpolate the coefficients along the line of constant x1. These 
    ! computations are independent from one another. A little problem is
    ! the redundancy in the computation of the cell and offset along each
    ! of the constant x2 lines, as this will be done by each call of 
    ! interpolate_value_aux(). This suggests that the proper refactoring
    ! of this function would have the cell and offset as arguments.
    coeffs_line_jm1 => spline%coeffs(0:num_pts_x1+2, cell-1)
    coeffs_line_j   => spline%coeffs(0:num_pts_x1+2, cell)
    coeffs_line_jp1 => spline%coeffs(0:num_pts_x1+2, cell+1)
    coeffs_line_jp2 => spline%coeffs(0:num_pts_x1+2, cell+2)
    cim1      = interpolate_derivative_aux(x1, x1_min, rh1, coeffs_line_jm1)
    ci        = interpolate_derivative_aux(x1, x1_min, rh1, coeffs_line_j  )
    cip1      = interpolate_derivative_aux(x1, x1_min, rh1, coeffs_line_jp1)
    cip2      = interpolate_derivative_aux(x1, x1_min, rh1, coeffs_line_jp2)
    t1        = 3.0_f64*ci
    t3        = 3.0_f64*cip1
    t2        = cdx*(cdx*(cdx*(cim1 - t1) + t1) + t1) + ci
    t4        =  dx*( dx*( dx*(cip2 - t3) + t3) + t3) + cip1
    interpolate_x1_derivative_2D = (1.0_f64/6.0_f64)*(t2 + t4)
  end function interpolate_x1_derivative_2D

  ! interpolate_x2_derivative_2D(): given discrete data f(i,j) that are
  ! described by a 2-dimensional cubic spline fit s(x1,x2), where the
  ! continuous variables x1 and x2 are within the original limits of i and j
  ! respectively, interpolate_x1_derivative() returns the value of
  !
  !         partial s
  !       -------------
  !         partial x2
  !
  ! evaluated at the point (x1,x2). (Sorry for the ambiguous use of x1)
  function interpolate_x2_derivative_2D( x1, x2, spline )
    sll_real64                          :: interpolate_x2_derivative_2D
    intrinsic                           :: associated, int, real
    sll_real64, intent(in)              :: x1
    sll_real64, intent(in)              :: x2
#ifdef STDF95
    type(sll_cubic_spline_2D)                 :: spline
#else
    type(sll_cubic_spline_2D), pointer        :: spline
#endif
    sll_real64                          :: rh1   ! reciprocal of cell spacing
    sll_real64                          :: rh2   ! reciprocal of cell spacing
    sll_int32                           :: cell
    sll_real64                          :: dx
    sll_real64                          :: cdx  ! 1-dx
    sll_real64                          :: t0   ! temp/scratch variables ...
    sll_real64                          :: t1
    sll_real64                          :: t2
    sll_real64                          :: t3
    sll_real64                          :: cim1 ! C_(i-1)
    sll_real64                          :: ci   ! C_i
    sll_real64                          :: cip1 ! C_(i+1)
    sll_real64                          :: cip2 ! C_(i+2)
    sll_real64                          :: x1_min
    sll_real64                          :: x2_min
    sll_int32                           :: num_pts_x1
    sll_int32                           :: num_pts_x2
    sll_real64, dimension(:), pointer   :: coeffs_line_jm1
    sll_real64, dimension(:), pointer   :: coeffs_line_j
    sll_real64, dimension(:), pointer   :: coeffs_line_jp1
    sll_real64, dimension(:), pointer   :: coeffs_line_jp2
    ! We set these as assertions since we want the flexibility of turning
    ! them off.
    SLL_ASSERT( (x1 .ge. spline%x1_min) .and. (x1 .le. spline%x1_max) )
    SLL_ASSERT( (x2 .ge. spline%x2_min) .and. (x2 .le. spline%x2_max) )
#ifdef STDF95
#else
    SLL_ASSERT( associated(spline) )
#endif
    x1_min     = spline%x1_min
    x2_min     = spline%x2_min
    num_pts_x1 = spline%num_pts_x1
    num_pts_x2 = spline%num_pts_x2
    rh1        = spline%x1_rdelta
    rh2        = spline%x2_rdelta
    ! find the cell and offset for x2
    t0         = (x2-x2_min)*rh2
    cell       = int(t0) + 1
    dx         = t0 - real(cell-1)
    cdx        = 1.0_f64 - dx
    !  write (*,'(a,i8, a, f20.12)') 'cell = ', cell, ',   dx = ', dx
    ! interpolate the coefficients along the line of constant x1. These 
    ! computations are independent from one another. A little problem is
    ! the redundancy in the computation of the cell and offset along each
    ! of the constant x2 lines, as this will be done by each call of 
    ! interpolate_value_aux(). This suggests that the proper refactoring
    ! of this function would have the cell and offset as arguments.
    coeffs_line_jm1 => spline%coeffs(0:num_pts_x1+2, cell-1)
    coeffs_line_j   => spline%coeffs(0:num_pts_x1+2, cell)
    coeffs_line_jp1 => spline%coeffs(0:num_pts_x1+2, cell+1)
    coeffs_line_jp2 => spline%coeffs(0:num_pts_x1+2, cell+2)
    cim1      = interpolate_value_aux(x1, x1_min, rh1, coeffs_line_jm1)
    ci        = interpolate_value_aux(x1, x1_min, rh1, coeffs_line_j  )
    cip1      = interpolate_value_aux(x1, x1_min, rh1, coeffs_line_jp1)
    cip2      = interpolate_value_aux(x1, x1_min, rh1, coeffs_line_jp2)
    t1 = 2.0_f64*(cim1 - 2.0_f64*ci + cip1)
    t2 = -cim1 + 3.0_f64*(ci - cip1) + cip2
    t3 =  cip1 - cim1
    interpolate_x2_derivative_2D = 0.5_f64*rh2*(dx*(t1 + dx*t2) + t3)
  end function interpolate_x2_derivative_2D


  subroutine delete_spline_2D( spline )
    type(sll_cubic_spline_2D), pointer :: spline
    sll_int32                    :: ierr
    ! Fixme: some error checking, whether the spline pointer is associated
    ! for instance
    if( .not. associated(spline) ) then
       print *, 'delete_spline_2D(): passed spline is not associated'
       STOP
    end if
    SLL_DEALLOCATE( spline%d1, ierr )
    SLL_DEALLOCATE( spline%d2, ierr )
    SLL_DEALLOCATE( spline%coeffs, ierr )
    spline%data => null()
    SLL_DEALLOCATE( spline, ierr )
  end subroutine delete_spline_2D


#undef NUM_TERMS
end module sll_cubic_splines

