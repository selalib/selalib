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

  ! At all cost we want to avoid the use of obscure numeric flags to
  ! induce some function behavior. Here we use a Fortran2003 feature, the
  ! enumeration. It is thus possible to give descriptive names to flags,
  ! instead of using some numeric code that one needs to look up somewhere.

  enum, bind(C)
     enumerator :: PERIODIC_SPLINE = 0, HERMITE_SPLINE = 1
  end enum
  
  interface delete
     module procedure delete_spline_1D
  end interface

contains  ! ****************************************************************

  ! The following could be changed to a direct access eventually, since it
  ! is hard to conceive that the slopes would ever be anything different than
  ! a double precision value at the top of the sll_spline_1D object. For now,
  ! these are the only slots inside the spline that are meant to be modified
  ! outside of the initialization or spline computation functions.

  subroutine set_slope_left(spline, value)
    type(sll_spline_1D), pointer :: spline
    sll_real64, intent(in)       :: value
    if( .not. associated(spline) ) then
       print *, 'set_slope_left(): not associated spline objet passed.'
       STOP
    end if
    spline%slope_l = value
  end subroutine set_slope_left

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
       print *, 'ERROR, new_spline_1D: Because of the algorithm used, this function is meant to be used with arrays that are at least of size = 28'
       STOP 'new_spline_1D()'
    end if
    if( xmin .gt. xmax ) then
       print *, 'ERROR, new_spline_1D: xmin is greater than xmax, this would cause all sorts of errors.'
       STOP
    end if
    ! Some more general error checking depending on the type of boundary
    ! condition requested.
    select case (bc_type)
    case (PERIODIC_SPLINE)
       if( present(sl) .or. present(sr) ) then
          print *, 'new_spline_1D(): it is not allowed to specify the end slopes in the case of periodic boundary conditions. Exiting program...'
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
  subroutine compute_spline_1D_periodic( f, spline )
    sll_real64, dimension(:), intent(in) :: f    ! data to be fit
    type(sll_spline_1D), pointer         :: spline
    sll_real64, dimension(:), pointer    :: coeffs
    sll_int32                         :: i
    sll_int32                         :: np
    sll_real64, parameter             :: a=sqrt((2.0_f64+sqrt(3.0_f64))/6.0_f64)
    sll_real64, parameter             :: r_a = 1.0_f64/a
    sll_real64, parameter             :: b=sqrt((2.0_f64-sqrt(3.0_f64))/6.0_f64)
    sll_real64, parameter             :: b_a = b/a
    sll_real64                        :: coeff_tmp
    sll_real64                        :: d1
    sll_real64, dimension(:), pointer :: d

    if( .not. associated(spline) ) then
       ! FIXME: THROW ERROR
       print *, 'ERROR: compute_spline_1D_periodic(): uninitialized spline object passed as argument. Exiting... '
       STOP
    end if
    if( .not. (size(f) .ge. spline%n_points ) ) then
       ! FIXME: THROW ERROR
       print *, 'ERROR: compute_spline_1D_periodic(): '
       write (*,'(a, i8, a, i8)') 'spline object needs data of size >= ', &
            spline%n_points, ' . Passed size: ', size(f)
       STOP
    end if
    np     =  spline%n_points
    d      => spline%d
    coeffs => spline%coeffs
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
  end subroutine compute_spline_1D_periodic

  subroutine compute_spline_1D_hermite( f, spline )
    sll_real64, dimension(:), intent(in) :: f    ! data to be fit
    type(sll_spline_1D), pointer      :: spline
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
    sll_real64, dimension(:), pointer :: d
    sll_real64                        :: slope_l
    sll_real64                        :: slope_r
    sll_real64                        :: delta
    sll_real64                        :: f1   ! to store modified value of f(1)
    sll_real64                        :: fnp  ! for modified value of f(np)

    if( .not. associated(spline) ) then
       ! FIXME: THROW ERROR
       print *, 'ERROR: compute_spline_1D_hermite(): uninitialized spline object passed as argument. Exiting... '
       STOP
    end if
    if( .not. (size(f) .ge. spline%n_points ) ) then
       ! FIXME: THROW ERROR
       print *, 'ERROR: compute_spline_1D_hermite(): '
       write (*,'(a, i8, a, i8)') 'spline object needs data of size >= ', &
            spline%n_points, ' . Passed size: ', size(f)
       STOP
    end if
    np      =  spline%n_points
    d       => spline%d
    coeffs  => spline%coeffs
    slope_l = spline%slope_L
    slope_r = spline%slope_R
    delta   = spline%delta
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
  end subroutine compute_spline_1D_hermite

#undef NUM_TERMS


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
  function interpolate_value( x, spline )
    sll_real64                        :: interpolate_value
    intrinsic                         :: associated, int, real
    sll_real64, intent(in)            :: x
    type(sll_spline_1D), pointer      :: spline
    sll_real64, dimension(:), pointer :: coeffs
    sll_real64                        :: rh   ! reciprocal of cell spacing
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
    sll_int32                         :: num_cells
    ! We set these as assertions since we want the flexibility of turning
    ! them off.
    SLL_ASSERT( (x .ge. spline%xmin) .and. (x .le. spline%xmax) )
    SLL_ASSERT( associated(spline) )
    ! FIXME: arg checks here
    num_cells = spline%n_points-1
    rh        = spline%rdelta
    coeffs    => spline%coeffs
    ! find the cell and offset for x
    t0        = x*rh
    cell      = int(t0) + 1
    dx        = t0 - real(cell-1)
    cdx       = 1.0_f64 - dx
    !  write (*,'(a,i8, a, f20.12)') 'cell = ', cell, ',   dx = ', dx
    cim1      = coeffs(cell-1)
    ci        = coeffs(cell)
    cip1      = coeffs(cell+1)
    cip2      = coeffs(cell+2)
    t1        = 3.0_f64*ci
    t3        = 3.0_f64*cip1
    t2        = cdx*(cdx*(cdx*(cim1 - t1) + t1) + t1) + ci
    t4        =  dx*( dx*( dx*(cip2 - t3) + t3) + t3) + cip1
    interpolate_value = (1.0_f64/6.0_f64)*(t2 + t4)
  end function interpolate_value

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


  ! July 25, 2011: Tried to use the type-bound procedure approach
  ! to write a one-parameter function like:
  !
  ! function interpolate_val( this, x)
  !
  ! which, being type-bound, could be called like:
  !
  ! spline_object%interpolate(x)
  !
  ! In principle, if the above would be recognized as the function signature
  ! by the compiler, this could mean that the simple gauss_legendre-type 
  ! integrator could be used to integrate the spline function.
  ! This approach failed on various grounds:
  ! gfortran still does not implement polymorphic types (class construct).
  ! This is not a showstopper. However, upon implementing the type-bound
  ! function and passing it to gauss_legendre_integral_1D(), I received an
  ! error saying that the type-bound function was expected with arguments.
  ! It could be that this can be fixed by setting the appropriate function
  ! pointers, but at this moment I decided not to go through this and 
  ! abandon this approach until it is truly important to have something like
  ! this.
  ! I leave this here for now but eventually this should be deleted.
#if 0
  function interpolate_val( this, x )
    sll_real64                        :: interpolate_val
    intrinsic                         :: associated, int, real
    sll_real64, intent(in)            :: x
    type(sll_spline_1D)               :: this ! class is not recognized...
    sll_real64, dimension(:), pointer :: coeffs
    sll_real64                        :: rh   ! reciprocal of cell spacing
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
    sll_int32                         :: num_cells

    SLL_ASSERT( (x .ge. this%xmin) .and. (x .le. this%xmax) )
!    SLL_ASSERT( associated(this) )
    ! FIXME: arg checks here
    num_cells = this%n_cells
    rh        = this%rdelta
    coeffs    => this%c
    ! find the cell and offset for x
    t0        = x*rh
    cell      = int(t0) + 1
    dx        = t0 - real(cell-1)
    cdx       = 1.0_f64 - dx
    !  write (*,'(a,i8, a, f20.12)') 'cell = ', cell, ',   dx = ', dx
    cim1      = coeffs(cell-1)
    ci        = coeffs(cell)
    cip1      = coeffs(cell+1)
    cip2      = coeffs(cell+2)
    t1        = 3.0_f64*ci
    t3        = 3.0_f64*cip1
    t2        = cdx*(cdx*(cdx*(cim1 - t1) + t1) + t1) + ci
    t4        =  dx*( dx*( dx*(cip2 - t3) + t3) + t3) + cip1
    interpolate_val = (1.0_f64/6.0_f64)*(t2 + t4)
  end function interpolate_val
#endif

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
end module sll_splines
