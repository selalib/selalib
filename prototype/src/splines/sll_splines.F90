!> \brief  
!> The splines module provides capabilities for 1D data interpolation with cubic B-splines
!> and different boundary conditions
!>
!> (at the time of this writing: periodic, hermite). The data to be interpolated is represented by a 
!> simple array.  The spline coefficients and other information are stored in a spline object, 
!> which is also used to interpolate the fitted data.
!> 
module sll_splines
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
  implicit none
  
  type sll_spline_1D
     sll_int32                         :: n_points ! size
     sll_real64                        :: delta    ! discretization step
     sll_real64                        :: rdelta   ! reciprocal of delta
     sll_real64                        :: xmin
     sll_real64                        :: xmax
     sll_int32                         :: bc_type  ! periodic, hermite
     sll_real64, dimension(:), pointer :: data     ! data for the spline fit
     sll_real64, dimension(:), pointer :: d        ! scratch space D (L*D = F),
                                                   ! refer to algorithm below.
                                                   ! Size depends on BC's.
     sll_real64, dimension(:), pointer :: coeffs   ! the spline coefficients
     sll_real64                        :: slope_L  ! left slope, for Hermite
     sll_real64                        :: slope_R  ! right slope, for Hermite
  end type sll_spline_1D

  ! Are x1 and x2 the coordinates that we should use? Or are eta1 and eta2
  ! the more logical ones given that the splines live in an uniform mesh??
  ! Need to define this and then be consistent throughout.
  type sll_spline_2D
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
     sll_int32                           :: bc_mix_identifier
     ! if data is not used, it should be deleted make a decision...
     sll_real64, dimension(:,:), pointer :: data   ! data for the spline fit
     sll_real64, dimension(:), pointer   :: d1     ! scratch space D (L*D = F),
                                                   ! refer to algorithm below.
                                                   ! Size depends on BC's.
     sll_real64, dimension(:), pointer   :: d2     ! Second scratch space
     sll_real64, dimension(:,:), pointer :: coeffs   ! the spline coefficients
     sll_real64                          :: x1_min_slope  ! for Hermite BCs
     sll_real64                          :: x1_max_slope  ! for Hermite BCs
     sll_real64                          :: x2_min_slope  ! for Hermite BCs
     sll_real64                          :: x2_max_slope  ! for Hermite BCs
  end type sll_spline_2D


  ! At all cost we want to avoid the use of obscure numeric flags to
  ! induce some function behavior. Here we use a Fortran2003 feature, the
  ! enumeration. It is thus possible to give descriptive names to flags,
  ! instead of using some numeric code that one needs to look up somewhere.

  enum, bind(C)
     enumerator :: PERIODIC_SPLINE = 0, HERMITE_SPLINE = 1
  end enum
  
  interface delete
     module procedure delete_spline_1D, delete_spline_2D
  end interface

contains  ! ****************************************************************

  ! The following could be changed to a direct access eventually, since it
  ! is hard to conceive that the slopes would ever be anything different than
  ! a double precision value at the top of the sll_spline_1D object. For now,
  ! these are the only slots inside the spline that are meant to be modified
  ! outside of the initialization or spline computation functions.
  
  !> set_slope_left
  !> \param[in] spline object
  !> \param[in] value  set derivative on left hand side to value
  subroutine set_slope_left(spline, value)
    type(sll_spline_1D), pointer :: spline
    sll_real64, intent(in)       :: value
    if( .not. associated(spline) ) then
       print *, 'set_slope_left(): not associated spline objet passed.'
       STOP
    end if
    spline%slope_l = value
  end subroutine set_slope_left
  
  !> set slope on right hand side for hermite boundary conditions
  !> \param[in] spline object
  !> \param[in] value to which set the derivative on the right hand side
  subroutine set_slope_right(spline, value)
    type(sll_spline_1D), pointer :: spline
    sll_real64, intent(in)       :: value
    if( .not. associated(spline) ) then
       print *, 'set_slope_left(): not associated spline objet passed.'
       STOP
    end if
    spline%slope_r = value
  end subroutine set_slope_right
  
  
  ! The following implementation embodies the algorithm described in
  ! Eric Sonnendrucker's "A possibly faster algorithm for cubic splines on
  ! a uniform grid" (unpublished).
  
  ! The array of spline coefficients has NP+3 elements. The extra elements
  ! at the ends (i.e.: 0, NP+1, NP+2) store coefficients whose values are
  ! determined by the type of boundary condition used. This is invisible
  ! to the user, who should not be concerned with this implementation detail.
  
  !> create new spline object
  function new_spline_1D( num_points, xmin, xmax, bc_type, sl, sr )
    type(sll_spline_1D), pointer         :: new_spline_1D
    sll_int32,  intent(in)               :: num_points
    sll_real64, intent(in)               :: xmin
    sll_real64, intent(in)               :: xmax
    sll_int32,  intent(in)               :: bc_type
    sll_real64, intent(in), optional     :: sl
    sll_real64, intent(in), optional     :: sr
    sll_int32                            :: ierr
    SLL_ALLOCATE( new_spline_1D, ierr )
    new_spline_1D%n_points = num_points
    new_spline_1D%xmin     = xmin
    new_spline_1D%xmax     = xmax
    new_spline_1D%delta    = (xmax - xmin)/real((num_points-1),f64)
    new_spline_1D%rdelta   = 1.0_f64/new_spline_1D%delta
    new_spline_1D%bc_type  = bc_type
    if( num_points .le. 28 ) then
       print *, 'ERROR, new_spline_1D: Because of the algorithm used, ', &
            'this function is meant to be used with arrays that are at ', &
            'least of size = 28'
       STOP 'new_spline_1D()'
    end if
    if( xmin .gt. xmax ) then
       print *, 'ERROR, new_spline_1D: xmin is greater than xmax, ', &
            'this would cause all sorts of errors.'
       STOP
    end if
    ! Some more general error checking depending on the type of boundary
    ! condition requested.
    select case (bc_type)
    case (PERIODIC_SPLINE)
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
    case (HERMITE_SPLINE)
       if( present(sl) ) then
          new_spline_1D%slope_L = sl
       else
          new_spline_1D%slope_L = 0.0  ! default left slope for Hermite case
       end if
       if( present(sr) ) then
          new_spline_1D%slope_R = sr
       else
          new_spline_1D%slope_R = 0.0  ! default right slope for Hermite case
       end if
    case default
       print *, 'ERROR: compute_spline_1D(): not recognized boundary condition'
       STOP
    end select
    SLL_ALLOCATE( new_spline_1D%d(num_points),   ierr )
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
  !> compute spline coefficients
  subroutine compute_spline_1D( f, bc_type, spline )
    sll_real64, dimension(:), intent(in) :: f    ! data to be fit
    sll_int32,  intent(in)               :: bc_type
    type(sll_spline_1D), pointer         :: spline
    ! Note that this function does no error checking and basically
    ! outsources this task to the functions it is wrapping around.
    ! This is so because those functions can be used independently
    ! (if the user wants to avoid the overhead of calling this
    ! wrapper function), so in any case, the error checking of
    ! the arguments will be carried out at least once.
    select case (bc_type)
    case (PERIODIC_SPLINE)
       call compute_spline_1D_periodic( f, spline )
    case (HERMITE_SPLINE)
       call compute_spline_1D_hermite( f, spline )
    case default
       print *, 'ERROR: compute_spline_1D(): not recognized boundary condition'
       STOP
    end select
  end subroutine compute_spline_1D

#define NUM_TERMS 27
  ! The following auxiliary functions:
  ! compute_spline_1D_periodic_aux() 
  ! compute_spline_1D_hermite_aux()
  ! are the fundamental building blocks. These are meant to do the work
  ! needed to compute the splines. Other functions are essentially 
  ! wrappers around these.  Clients of these routines are responsible for 
  ! all error-checking.
  subroutine compute_spline_1D_periodic_aux( f, num_pts, d, coeffs )
    sll_real64, dimension(:), pointer :: f
    sll_int32, intent(in)             :: num_pts
    sll_real64, dimension(:), pointer :: d
    sll_real64, dimension(:), pointer :: coeffs
    sll_real64, parameter             :: a=sqrt((2.0_f64+sqrt(3.0_f64))/6.0_f64)
    sll_real64, parameter             :: r_a = 1.0_f64/a
    sll_real64, parameter             :: b=sqrt((2.0_f64-sqrt(3.0_f64))/6.0_f64)
    sll_real64, parameter             :: b_a = b/a
    sll_real64                        :: coeff_tmp
    sll_real64                        :: d1
    sll_int32                         :: i
    sll_int32                         :: np
    SLL_ASSERT( size(f) .ge. num_pts )
    SLL_ASSERT( size(d) .ge. num_pts )
    SLL_ASSERT( size(coeffs) .ge. num_pts )
    SLL_ASSERT( (num_pts .ge. 0) .and. (num_pts .lt. NUM_TERMS))
    np     =  num_pts
    ! Compute d(1):
    d1 =  f(1)
    coeff_tmp = 1.0_f64
    do i = 0, NUM_TERMS-1  ! if NUM_TERMS == 0, only f(nc) is considered.
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
    coeffs(np-1) = d1*r_a
    ! rest of the coefficients:
    do i = np-2, 1, -1
       coeffs(i) = r_a*(d(i) - b*coeffs(i+1))
    end do
    coeffs(0)    = coeffs(np-1)
    coeffs(np)   = coeffs(1)
    coeffs(np+1) = coeffs(2)
    coeffs(np+2) = coeffs(3)
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
    sll_real64, parameter             :: a=sqrt((2.0_f64+sqrt(3.0_f64))/6.0_f64)
    sll_real64, parameter             :: r_a = 1.0_f64/a
    sll_real64, parameter             :: b=sqrt((2.0_f64-sqrt(3.0_f64))/6.0_f64)
    sll_real64, parameter             :: b_a = b/a
    sll_real64, parameter             :: ralpha = sqrt(6.0_f64/sqrt(3.0_f64))
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
    coeffs(np) = ralpha*d(np)
    do i = np-1, 1, -1
       coeffs(i) = r_a*(d(i) - b*coeffs(i+1))
    end do
    coeffs(0)    = coeffs(2)    - 2.0 * delta * slope_l
    coeffs(np+1) = coeffs(np-1) + 2.0 * delta * slope_r
    coeffs(np+2) = 0.0 !coeffs(np-2)  ! not used
  end subroutine compute_spline_1D_hermite_aux


  subroutine compute_spline_1D_periodic( f, spline )
    sll_real64, dimension(:), intent(in), target :: f    ! data to be fit
    type(sll_spline_1D), pointer         :: spline
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
    fp     => f
    np     =  spline%n_points
    d      => spline%d
    coeffs => spline%coeffs
    call compute_spline_1D_periodic_aux( fp, np, d, coeffs )
  end subroutine compute_spline_1D_periodic



  subroutine compute_spline_1D_hermite( f, spline )
    sll_real64, dimension(:), intent(in), target :: f    ! data to be fit
    type(sll_spline_1D), pointer      :: spline
    sll_real64, dimension(:), pointer :: coeffs
    sll_int32                         :: np
    sll_real64, dimension(:), pointer :: fp
    sll_real64, dimension(:), pointer :: d
    sll_real64                        :: slope_l
    sll_real64                        :: slope_r
    sll_real64                        :: delta

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
    coeffs  => spline%coeffs
    slope_l = spline%slope_L
    slope_r = spline%slope_R
    delta   = spline%delta
    call compute_spline_1D_hermite_aux(fp,np,d, slope_l,slope_r, delta, coeffs)
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
    cim1      = coeffs(cell-1)
    ci        = coeffs(cell)
    cip1      = coeffs(cell+1)
    cip2      = coeffs(cell+2)
    t1        = 3.0_f64*ci
    t3        = 3.0_f64*cip1
    t2        = cdx*(cdx*(cdx*(cim1 - t1) + t1) + t1) + ci
    t4        =  dx*( dx*( dx*(cip2 - t3) + t3) + t3) + cip1
    interpolate_value_aux = (1.0_f64/6.0_f64)*(t2 + t4)
  end function interpolate_value_aux
  
  !> get spline interpolate at point x
  function interpolate_value( x, spline )
    sll_real64                        :: interpolate_value
    intrinsic                         :: associated, int, real
    sll_real64, intent(in)            :: x
    type(sll_spline_1D), pointer      :: spline
    sll_real64, dimension(:), pointer :: coeffs
    sll_real64                        :: xmin
    sll_real64                        :: rh   ! reciprocal of cell spacing
    ! We set these as assertions since we want the flexibility of turning
    ! them off.
    SLL_ASSERT( (x .ge. spline%xmin) .and. (x .le. spline%xmax) )
    SLL_ASSERT( associated(spline) )
    xmin = spline%xmin
    rh        = spline%rdelta
    coeffs => spline%coeffs
    interpolate_value = interpolate_value_aux( x, xmin, rh, coeffs )
  end function interpolate_value
  
  !> get spline interpolate at array of points
  subroutine interpolate_array_values( a_in, a_out, n, spline )
    intrinsic                               :: associated, int, real
    sll_int32, intent(in)                   :: n
    sll_real64, dimension(1:n), intent(in)  :: a_in
    sll_real64, dimension(1:n), intent(out) :: a_out
    type(sll_spline_1D), pointer            :: spline
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

  function interpolate_derivative( x, spline )
    sll_real64                        :: interpolate_derivative
    intrinsic                         :: associated, int, real
    sll_real64, intent(in)            :: x
    type(sll_spline_1D), pointer      :: spline
    sll_real64, dimension(:), pointer :: coeffs
    sll_real64                        :: rh   ! reciprocal of cell spacing
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
    sll_int32                         :: num_cells
    ! We set these as assertions since we want the flexibility of turning
    ! them off.
    SLL_ASSERT( (x .ge. spline%xmin) .and. (x .le. spline%xmax) )
    SLL_ASSERT( associated(spline) )
    num_cells = spline%n_points-1
    rh        = spline%rdelta
    coeffs    => spline%coeffs
    ! find the cell and offset for x
    t0        = (x-spline%xmin)*rh
    cell      = int(t0) + 1
    dx        = t0 - real(cell-1)
    ! write (*,'(a,i8, a, f20.12)') 'cell = ', cell, ',   dx = ', dx
    cim1      = coeffs(cell-1)
    ci        = coeffs(cell)
    cip1      = coeffs(cell+1)
    cip2      = coeffs(cell+2)
    t1 = 2.0_f64*(cim1 - 2.0_f64*ci + cip1)
    t2 = -cim1 + 3.0_f64*(ci - cip1) + cip2
    t3 =  cip1 - cim1
    interpolate_derivative = 0.5_f64*rh*(dx*(t1 + dx*t2) + t3)
  end function interpolate_derivative

  subroutine delete_spline_1D( spline )
    type(sll_spline_1D), pointer :: spline
    sll_int32                    :: ierr
    ! Fixme: some error checking, whether the spline pointer is associated
    ! for instance
    SLL_ASSERT( associated(spline) )
    SLL_DEALLOCATE( spline%d, ierr )
    SLL_DEALLOCATE( spline%coeffs, ierr )
    spline%data => null()
    SLL_DEALLOCATE( spline, ierr )
  end subroutine delete_spline_1D

  !-----------------------------------------------------------------------
  !
  ! Functions and subroutines for the 2D spline.
  !
  !----------------------------------------------------------------------

  ! Provide an modification function for the values of the slopes at the
  ! endpoints of the 2D spline. Note that there is full redundancy in these
  ! routines. A single change in any of these should justify the creation of
  ! a macro.
  subroutine set_x1_min_slope(spline, value)
    type(sll_spline_2D), pointer :: spline
    sll_real64, intent(in)       :: value
    if( .not. associated(spline) ) then
       print *, 'set_x1_min_slope(): not associated spline objet passed.'
       STOP
    end if
    spline%x1_min_slope = value
  end subroutine set_x1_min_slope

  subroutine set_x1_max_slope(spline, value)
    type(sll_spline_2D), pointer :: spline
    sll_real64, intent(in)       :: value
    if( .not. associated(spline) ) then
       print *, 'set_x1_max_slope(): not associated spline objet passed.'
       STOP
    end if
    spline%x1_max_slope = value
  end subroutine set_x1_max_slope

  subroutine set_x2_min_slope(spline, value)
    type(sll_spline_2D), pointer :: spline
    sll_real64, intent(in)       :: value
    if( .not. associated(spline) ) then
       print *, 'set_x2_min_slope(): not associated spline objet passed.'
       STOP
    end if
    spline%x2_min_slope = value
  end subroutine set_x2_min_slope

  subroutine set_x2_max_slope(spline, value)
    type(sll_spline_2D), pointer :: spline
    sll_real64, intent(in)       :: value
    if( .not. associated(spline) ) then
       print *, 'set_x2_max_slope(): not associated spline objet passed.'
       STOP
    end if
    spline%x2_max_slope = value
  end subroutine set_x2_max_slope

  function new_spline_2D( &
    num_pts_x1,   &
    num_pts_x2,   &
    x1_min,       &
    x1_max,       &
    x2_min,       &
    x2_max,       &
    x1_bc_type,   &
    x2_bc_type,   &
    x1_min_slope, &
    x1_max_slope, &
    x2_min_slope, &
    x2_max_slope &
    )

    type(sll_spline_2D), pointer         :: new_spline_2D
    sll_int32,  intent(in)               :: num_pts_x1
    sll_int32,  intent(in)               :: num_pts_x2
    sll_real64, intent(in)               :: x1_min
    sll_real64, intent(in)               :: x1_max
    sll_real64, intent(in)               :: x2_min
    sll_real64, intent(in)               :: x2_max
    sll_int32,  intent(in)               :: x1_bc_type
    sll_int32,  intent(in)               :: x2_bc_type
    sll_real64, intent(in), optional     :: x1_min_slope
    sll_real64, intent(in), optional     :: x1_max_slope
    sll_real64, intent(in), optional     :: x2_min_slope
    sll_real64, intent(in), optional     :: x2_max_slope
    sll_int32                            :: bc_selector
    sll_int32                            :: ierr
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

    ! Treat the bc_selector variable essentially like a bit field, to 
    ! accumulate the information on the different boundary conditions
    ! given. This scheme allows to add more types of boundary conditions
    ! if necessary.
    bc_selector = 0
    if( x1_bc_type .eq. HERMITE_SPLINE ) then 
       bc_selector = bc_selector + 1
    else if ( x2_bc_type .eq. HERMITE_SPLINE ) then
       bc_selector = bc_selector + 2
    end if
    new_spline_2D%bc_mix_identifier = bc_selector

    select case (bc_selector)
    case ( 0 ) 
       ! both boundary condition types are periodic
       if( &
          present(x1_min_slope) .or. present(x1_max_slope) .or. &
          present(x2_min_slope) .or. present(x2_max_slope) ) then

          print *, 'new_spline_2D(): it is not allowed to specify the end', &
               'slopes in the case of periodic boundary conditions.', &
               'Exiting program...'
          STOP 'new_spline_2D'
       else
          ! Assign some value, but this value should never be used in the
          ! full periodic case anyway.
          new_spline_2D%x1_min_slope = 0.0
          new_spline_2D%x1_max_slope = 0.0
          new_spline_2D%x2_min_slope = 0.0
          new_spline_2D%x2_max_slope = 0.0
       end if
    case ( 1 ) 
       ! Hermite condition in X1 and periodic in X2 
       if( present(x2_min_slope) .or. present(x2_max_slope) ) then
          print *, 'new_spline_2D(): it is not allowed to specify the end', &
               'slopes in the case of periodic boundary conditions.', &
               'Exiting program...'
          STOP 'new_spline_2D'
       end if
       if ( present(x1_min_slope) ) then
          ! need to check if the slopes are present, and to apply default 
          ! values if not.
          new_spline_2D%x1_min_slope = x1_min_slope
       else
          ! apply default value for the slope
          new_spline_2D%x1_min_slope = 0.0
       end if
       if ( present(x1_max_slope) ) then
          new_spline_2D%x1_max_slope = x1_max_slope
       else
          ! apply default value
          new_spline_2D%x1_max_slope = 0.0
       end if
       new_spline_2D%x2_min_slope = 0.0
       new_spline_2D%x2_max_slope = 0.0
    case( 2 )
       ! Periodic in X1 and Hermite in X2
       if( present(x1_min_slope) .or. present(x1_max_slope) ) then
          print *, 'new_spline_2D(): it is not allowed to specify the end', &
               'slopes in the case of periodic boundary conditions.', &
               'Exiting program...'
          STOP 'new_spline_2D'
       end if
       if ( present(x2_min_slope) ) then
          ! need to check if the slopes are present, and to apply default 
          ! values if not.
          new_spline_2D%x2_min_slope = x2_min_slope
       else
          ! apply default value for the slope
          new_spline_2D%x2_min_slope = 0.0
       end if
       if ( present(x2_max_slope) ) then
          new_spline_2D%x2_max_slope = x2_max_slope
       else
          ! apply default value
          new_spline_2D%x2_max_slope = 0.0
       end if
       ! set the periodic conditions that will not be used
       new_spline_2D%x1_min_slope = 0.0
       new_spline_2D%x1_max_slope = 0.0
    case( 3 )
       ! Hermite conditions in both, X1 and X2
       if ( present(x1_min_slope) ) then
          ! need to check if the slopes are present, and to apply default 
          ! values if not.
          new_spline_2D%x1_min_slope = x1_min_slope
       else
          ! apply default value for the slope
          new_spline_2D%x1_min_slope = 0.0
       end if
       if ( present(x1_max_slope) ) then
          ! need to check if the slopes are present, and to apply default 
          ! values if not.
          new_spline_2D%x1_max_slope = x1_max_slope
       else
          ! apply default value for the slope
          new_spline_2D%x1_max_slope = 0.0
       end if
       if ( present(x2_min_slope) ) then
          ! need to check if the slopes are present, and to apply default 
          ! values if not.
          new_spline_2D%x2_min_slope = x2_min_slope
       else
          ! apply default value for the slope
          new_spline_2D%x2_min_slope = 0.0
       end if
       if ( present(x2_max_slope) ) then
          new_spline_2D%x2_max_slope = x2_max_slope
       else
          ! apply default value
          new_spline_2D%x2_max_slope = 0.0
       end if
    case default
       print *, 'ERROR: new_spline_2D(): ', &
            'did not recognize given boundary conditions.'
       STOP
    end select
    SLL_ALLOCATE( new_spline_2D%d1(num_pts_x1),   ierr )
    SLL_ALLOCATE( new_spline_2D%d2(num_pts_x2),   ierr )
    ! Reminder: Fortran arrays are column-major ordered...
    ! Note: The indexing of the coefficients array includes the end-
    ! points 0, num_points, num_points+1, num_points+2. These are meant to 
    ! store the boundary condition-specific data. The 'periodic' BC does
    ! not use the num_points+2 point.
    SLL_ALLOCATE( new_spline_2D%coeffs(0:num_pts_x1+2,0:num_pts_x2+2), ierr )
  end function new_spline_2D

  subroutine compute_spline_2D_prdc_prdc( data, spline )
    sll_real64, dimension(:,:), intent(in), target :: data  ! data to be fit
    type(sll_spline_2D), pointer         :: spline
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
    if( .not. (size(data,1) .ge. spline%num_pts_x1 ) ) then
       ! FIXME: THROW ERROR
       print *, 'ERROR: compute_spline_2D_prdc_prdc(): '
       write (*,'(a, i8, a, i8)') 'spline object needs data of size >= ', &
            spline%num_pts_x1, ' . Passed size: ', size(data,1)
       STOP
    end if
    if( .not. (size(data,2) .ge. spline%num_pts_x2 ) ) then
       ! FIXME: THROW ERROR
       print *, 'ERROR: compute_spline_2D_prdc_prdc(): '
       write (*,'(a, i8, a, i8)') 'spline object needs data of size >= ', &
            spline%num_pts_x2, ' . Passed size: ', size(data,2)
       STOP
    end if
print *, 'assigning local variables'
    npx1   =  spline%num_pts_x1
    npx2   =  spline%num_pts_x2
    d1     => spline%d1
    d2     => spline%d2
    ! build splines along the x2 direction. Note: due to Fortran's 
    ! column-major ordering, this uses long strides in memory.
    do j=1,npx1
       datap  => data(1:npx2,j)
       coeffs => spline%coeffs(1:npx2+2,j)
       call compute_spline_1D_periodic_aux( datap, npx1, d1, coeffs )
    end do
    ! build splines along the x1 direction. Note: due to Fortran's 
    ! column-major ordering, this involves short strides in memory.
    do i=1,npx2
       datap  => data(i,1:npx1)
       coeffs => spline%coeffs(i,1:npx1+2)
       call compute_spline_1D_periodic_aux( datap, npx1, d1, coeffs )
    end do

  end subroutine compute_spline_2D_prdc_prdc


  function interpolate_value_2D( x1, x2, spline )
    sll_real64                          :: interpolate_value_2D
    intrinsic                           :: associated, int, real
    sll_real64, intent(in)              :: x1
    sll_real64, intent(in)              :: x2
    type(sll_spline_2D), pointer        :: spline
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
    num_pts_x1 = spline%num_pts_x1
    num_pts_x2 = spline%num_pts_x2
    rh1        = spline%x1_rdelta
    rh2        = spline%x2_rdelta
    ! find the cell and offset for x2
    t0         = x2*rh2
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



  subroutine delete_spline_2D( spline )
    type(sll_spline_2D), pointer :: spline
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

  ! b_splines_at_x() returns the values of all the B-splines of a given 
  ! degree that have support in cell 'cell' and evaluated at the point 'x'. 
  ! In other words, if B[j,i](x) is the spline of degree 'j' whose leftmost 
  ! support is at cell 'i' and evaluated at 'x', then b_splines_at_x returns 
  ! the sequence (in the form of an array):
  ! 
  ! B[j,i-degree](x), B[j,i-degree+1](x), B[j,i-degree+2](x), ..., B[j,i](x)
  !
  ! Implementation notes: 
  !
  ! It would have been very simple and convenient to implement this with
  ! a recursion, since:
  !
  !               x-t(i)                       t(i+j+1)-x
  ! B[j,i](x) = ----------- * B[j-1,i](x) + ----------------* B[j-1,i+1](x)
  !             t(i+j)-t(i)                  t(i+j+1)-t(i+1)
  !
  ! and
  !
  ! B[0,i] = 1 if t(i) <= x < t(i+1), and 0 otherwise.
  !
  ! More generally:
  !
  ! if t(i) <= x < t(i+j+1) then the formula above applies but B[j,i](x) = 0
  ! otherwise.
  !
  ! The problem with the above recursion is that it will end up computing the
  ! splines of lower orders very redundantly, much like the problem of 
  ! computing a Fibonacci sequence with a recursion. For such a critical
  ! function (this will be present in inner loops), this is not acceptable.
  !
  ! Here we try a different approach but still use the idea of the 
  ! recursion formula above. We use the fact that for the desired sequence
  ! of splines of degree 'J', we need information within 2*J+1 cells:
  ! (i-J):(i+J). We populate these with the values of the B[0,i](x) splines
  ! and iteratively build the higher order splines as needed.

  function b_splines_at_x( knot_positions, spline_degree, cell, x )
    sll_int32, intent(in)                      :: spline_degree
    sll_real64, dimension(1:spline_degree+1)   :: b_splines_at_x
    sll_real64, dimension(:), intent(in)       :: knot_positions
    sll_int32, intent(in)                      :: cell
    sll_real64, intent(in)                     :: x
    sll_real64, dimension(1:2*spline_degree+1) :: splines
    sll_int32                                  :: i
    sll_int32                                  :: j
    sll_int32                                  :: last
    sll_real64                                 :: ti     ! t(i)
    sll_real64                                 :: tip1   ! t(i+1)
    sll_real64                                 :: tipj   ! t(i+j)
    sll_real64                                 :: tipjp1 ! t(i+j+1)
    sll_real64                                 :: fac1
    sll_real64                                 :: fac2
    sll_int32                                  :: current
    ! what argument checking to do here? Assertions only...

    ! Build the zeroth-order splines. The middle cell of the splines array
    ! corresponds to the 'cell' given as argument. So for example,
    ! splines(spline_degree+1) will have the value 1.0 if 'x' is inside 'cell'.
    do i=1,2*spline_degree+1
       if((x .ge. knot_positions(cell - spline_degree + i - 1)) .and. &
          (x .lt. knot_positions(cell - spline_degree + i))) then
          splines(i) = 1.0
       else
          splines(i) = 0.0
       end if
    end do
    
    ! Build the higher order splines. All of this work is redundant in 
    ! case that 'x' is not in the support of any of the splines. We might
    ! want to check for this condition and return accordingly.
    last = 2*spline_degree  
    do j=1,spline_degree
       do i=1,last
          current    = cell - spline_degree + i - 1
          ti         = knot_positions(current)
          tip1       = knot_positions(current+1)
          tipj       = knot_positions(current+j)
          tipjp1     = knot_positions(current+j+1)
          ! This is a dangerous situation for which we need some sort of
          ! protection: What guarantees are there that these denominators
          ! will not be zero?? This should probably be error-checked, else
          ! one can just end up with an array of NaN's.
          fac1       = (x - ti)/(tipj - ti)
          fac2       = (tipjp1 - x)/(tipjp1 - tip1)
          splines(i) = fac1*splines(i) + fac2*splines(i+1)
       end do
       last = last - 1
    end do
    b_splines_at_x(1:spline_degree+1) = splines(1:spline_degree+1)
  end function b_splines_at_x

  ! b_spline_derivatives_at_x() returns an array with the values of the
  ! B-splines of a requested order that are supported in 'cell' and evaluated
  ! at 'x'. The return value has the format:
  !
  ! B'[j,i-degree](x), B'[j,i-degree+1](x), B'[j,i-degree+2](x),..., B'[j,i](x)
  !
  ! where 'j' is the requested degree of the spline.
  function b_spline_derivatives_at_x( knots, spline_degree, cell, x )
    sll_int32, intent(in)                    :: spline_degree
    sll_real64, dimension(spline_degree+1)   :: b_spline_derivatives_at_x
    sll_real64, dimension(:), intent(in)     :: knots
    sll_int32, intent(in)                    :: cell
    sll_real64, intent(in)                   :: x
    sll_real64, dimension(spline_degree+1)   :: derivs
    sll_real64, dimension(2*spline_degree+1) :: splines 
    sll_int32                                :: current
    sll_real64                               :: delta_x
    sll_real64                               :: fac1
    sll_real64                               :: fac2
    sll_int32                                :: i
    sll_int32                                :: j
    sll_int32                                :: last
    sll_real64                               :: ti
    sll_real64                               :: tip1
    sll_real64                               :: tipj
    sll_real64                               :: tipjp1
    ! FIXME: ARGUMENT CHECKS
    ! Compute derivatives of the splines of order spline_degree.
    ! what argument checking to do here? Use assertions only...

    ! FIXME, MORTAL SIN: HERE WE HAVE DUPLICATED CODE WITH THE PREVIOUS
    ! FUNCTION AND EXPECT TO DUPLICATE THIS EVEN MORE WITH A FUNCTION
    ! THAT FURTHER COMBINES VALUES AND DERIVATIVES. THIS IS NOT ACCEPTABLE.
    ! EVENTUALLY THIS CODE SEGMENT SHOULD BE MACROIZED. LEAVE AS IS FOR THE
    ! MOMENT SINCE WE NEED TO FIX THE ISSUE OF POSSIBLY ZERO VALUES IN
    ! DENOMINATORS AT LEAST. PRODUCTION VERSION SHOULD NOT HAVE DUPLICATED
    ! CODE IN THIS CASE.

    ! Build the zeroth-order splines. The middle cell of the splines array
    ! corresponds to the 'cell' given as argument. So for example,
    ! splines(spline_degree+1) will have the value 1.0 if 'x' is inside 'cell'.
    do i=1,2*spline_degree+1
       if((x .ge. knots(cell - spline_degree + i - 1)) .and. &
          (x .lt. knots(cell - spline_degree + i))) then
          splines(i) = 1.0
       else
          splines(i) = 0.0
       end if
    end do

    ! Build the higher order splines. All of this work is redundant in 
    ! case that 'x' is not in the support of any of the splines. We might
    ! want to check for this condition and return accordingly.
    last = 2*spline_degree  
    do j=1,spline_degree-1 ! we stop earlier to compute derivatives
       do i=1,last
          current    = cell - spline_degree + i - 1
          ti         = knots(current)
          tip1       = knots(current+1)
          tipj       = knots(current+j)
          tipjp1     = knots(current+j+1)
          ! This is a dangerous situation for which we need some sort of
          ! protection: What guarantees are there that these denominators
          ! will not be zero?? This should probably be error-checked, else
          ! one can just end up with an array of NaN's.
          fac1       = (x - ti)/(tipj - ti)
          fac2       = (tipjp1 - x)/(tipjp1 - tip1)
          splines(i) = fac1*splines(i) + fac2*splines(i+1)
       end do
       last = last - 1
    end do
    ! At this moment we have an array with values of the splines up to the
    ! order spline_degree - 1. Proceed to compute the derivatives of order
    ! spline_degree.
    do i=1,last
       current = cell - spline_degree + i - 1
       delta_x = knots(current+1) - knots(current)
       derivs(i) = (splines(i) - splines(i+1))/delta_x
    end do
    b_spline_derivatives_at_x(1:spline_degree+1) = derivs(1:spline_degree+1)
  end function b_spline_derivatives_at_x

#undef NUM_TERMS
end module sll_splines

