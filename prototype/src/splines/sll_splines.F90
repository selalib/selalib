module sll_splines
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
  implicit none

  type sll_spline_1D
     sll_int32                         :: n_cells  ! size
     sll_real64                        :: delta    ! discretization step
     sll_real64                        :: rdelta   ! reciprocal of delta
     sll_real64                        :: xmin
     sll_real64                        :: xmax
     sll_int32                         :: bc_type  ! periodic, hermite
     sll_real64, dimension(:), pointer :: data     ! data for the spline fit
     sll_real64, dimension(:), pointer :: d        ! scratch space D (L*D = F),
                                                   ! refer to algorithm below.
                                                   ! Size depends on BC's.
     sll_real64, dimension(:), pointer :: c        ! the spline coefficients
  end type sll_spline_1D

  ! At all cost we want to avoid the use of obscure numeric flags to
  ! induce some function behavior. Here we use a Fortran2003 feature, the
  ! enumeration. It is thus possible to give descriptive names to flags,
  ! instead of using some numeric code that one needs to look up somewhere.

  enum, bind(C)
     enumerator :: PERIODIC_SPLINE = 0, HERMITE_SPLINE = 1
  end enum
  
contains  ! ****************************************************************

  ! The following implementation embodies the algorithm described in
  ! Eric Sonnendrucker's "A possibly faster algorithm for cubic splines on
  ! a uniform grid" (unpublished).

  ! The array of spline coefficients has NC+4 elements. The extra elements
  ! at the ends (i.e.: 0, NC+2, NC+3) store coefficients whose values are
  ! determined by the type of boundary condition used. This is invisible
  ! to the user, who should not be concerned with this implementation detail.

  function new_spline_1D( num_cells, xmin, xmax, bc_type )
    type(sll_spline_1D), pointer         :: new_spline_1D
    !sll_real64, dimension(:), intent(in), target :: data
    sll_int32,  intent(in)               :: num_cells
    sll_real64, intent(in)               :: xmin
    sll_real64, intent(in)               :: xmax
    sll_int32,  intent(in)               :: bc_type
    sll_int32                            :: ierr
    sll_int32                            :: num_points
    SLL_ALLOCATE( new_spline_1D, ierr )
    new_spline_1D%n_cells = num_cells
    new_spline_1D%xmin    = xmin
    new_spline_1D%xmax    = xmax
    new_spline_1D%delta   = (xmax - xmin)/real(num_cells,f64)
    new_spline_1D%rdelta  = 1.0_f64/new_spline_1D%delta
    new_spline_1D%bc_type = bc_type
    !new_spline_1D%data    => data
    select case (bc_type)
    case (PERIODIC_SPLINE)
       num_points = num_cells
    case (HERMITE_SPLINE)
       num_points = num_cells + 1
    case default
       ! FIXME, throw error
       print *, 'new_spline_1D error: bc_type not recognized'
    end select
    SLL_ALLOCATE( new_spline_1D%d(num_points),   ierr )
    ! note how the indexing of the coefficients array includes the end-
    ! points 0, num_cells+1, num_cells+2, num_cells+3. These are meant to 
    ! store the boundary condition-specific data. The 'periodic' BC does
    ! not use the num_cells+3 point.
    SLL_ALLOCATE( new_spline_1D%c(0:num_cells+3), ierr )
    !call compute_spline_1D( data, num_cells, bc_type, new_spline_1D )
  end function new_spline_1D

  ! - data: the array whose data must be fit with the cubic spline.
  ! - nc: (number of cells; length of the data array that must be fit with 
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

  subroutine compute_spline_1D( f, nc, bc_type, spline )
    sll_real64, dimension(:), intent(in) :: f    ! data to be fit
    sll_int32,  intent(in)            :: nc
    sll_int32,  intent(in)            :: bc_type
    type(sll_spline_1D), pointer      :: spline

    select case (bc_type)
    case (PERIODIC_SPLINE)
       call compute_spline_1D_periodic( f, nc, spline )
    case (HERMITE_SPLINE)
       call compute_spline_1D_hermite( f, nc, spline )
    end select
  end subroutine compute_spline_1D

#define NUM_TERMS 27
  subroutine compute_spline_1D_periodic( f, nc, spline )
    sll_real64, dimension(:), intent(in) :: f    ! data to be fit
    sll_int32,  intent(in)            :: nc
    type(sll_spline_1D), pointer      :: spline
    sll_real64, dimension(:), pointer :: coeffs
    sll_int32                         :: i
!    sll_int32                         :: err
    sll_real64, parameter             :: a=sqrt((2.0_f64+sqrt(3.0_f64))/6.0_f64)
    sll_real64, parameter             :: r_a = 1.0_f64/a
    sll_real64, parameter             :: b=sqrt((2.0_f64-sqrt(3.0_f64))/6.0_f64)
    sll_real64, parameter             :: b_a = b/a
    sll_real64                        :: coeff_tmp
    sll_real64                        :: d1
    sll_real64, dimension(:), pointer :: d
 ! CHECK INPUT, nc is redundant, change interface
    SLL_ASSERT(nc == spline%n_cells)
    SLL_ASSERT(associated(spline))
    SLL_ASSERT(size(f) >= spline%n_cells)
    
    d      => spline%d
    coeffs => spline%c
    ! Compute d(1):
    d1 =  f(1)
    coeff_tmp = 1.0_f64
    do i = 0, NUM_TERMS-1  ! if NUM_TERMS == 0, only f(nc) is considered.
       coeff_tmp = coeff_tmp*(-b_a)
       d1 = d1 + coeff_tmp*f(nc-i)
    end do
    ! Fill the d array with the intermediate result
    d(1) = d1*r_a
    do i = 2,nc
       d(i) = r_a*(f(i) - b*d(i-1))
    end do
    ! Compute the coefficients. Start with first term
    d1        = d(nc)
    coeff_tmp = 1.0_f64
    do i = 1, NUM_TERMS
       coeff_tmp = coeff_tmp*(-b_a)
       d1 = d1 + coeff_tmp*d(i)
    end do
    coeffs(nc) = d1*r_a
    ! rest of the coefficients:
    do i = nc-1, 1, -1
       coeffs(i) = r_a*(d(i) - b*coeffs(i+1))
    end do
    coeffs(0)    = coeffs(nc)
    coeffs(nc+1) = coeffs(1)
    coeffs(nc+2) = coeffs(2)
    coeffs(nc+3) = coeffs(3)
  end subroutine compute_spline_1D_periodic

  subroutine compute_spline_1D_hermite( f, nc, spline )
    sll_real64, dimension(:), intent(in) :: f    ! data to be fit
    sll_int32,  intent(in)            :: nc
    type(sll_spline_1D), pointer      :: spline
    sll_real64, dimension(:), pointer :: coeffs
    sll_int32                         :: i
    sll_real64, parameter             :: a=sqrt((2.0_f64+sqrt(3.0_f64))/6.0_f64)
    sll_real64, parameter             :: r_a = 1.0_f64/a
    sll_real64, parameter             :: b=sqrt((2.0_f64-sqrt(3.0_f64))/6.0_f64)
    sll_real64, parameter             :: b_a = b/a
    sll_real64, parameter             :: ralpha = sqrt(6.0_f64/sqrt(3.0_f64))
    sll_real64                        :: coeff_tmp
    sll_real64                        :: d1
    sll_real64, dimension(:), pointer :: d

    ! CHECK INPUT, nc is redundant, change interface, see periodic
    d      => spline%d
    coeffs => spline%c
    ! Compute d(1):
    d1 =  f(1)
    coeff_tmp = 1.0_f64
    do i = 2, NUM_TERMS  ! if NUM_TERMS == 2, only f(2) is considered.
       coeff_tmp = coeff_tmp*(-b_a)
       d1 = d1 + coeff_tmp*f(i)
    end do
    ! Fill the d array with the intermediate result
    d(1) = d1*r_a
    do i = 2,nc
       d(i) = r_a*(f(i) - b*d(i-1))
    end do
    d(nc+1) = ralpha*(0.5_f64*f(nc+1) - b*d(nc))
    ! Compute the coefficients. Start with first term
    coeffs(nc+1) = ralpha*d(nc+1)
    do i = nc, 1, -1
       coeffs(i) = r_a*(d(i) - b*coeffs(i+1))
    end do
    coeffs(0)    = coeffs(2)
    coeffs(nc+2) = coeffs(nc)
    coeffs(nc+3) = coeffs(nc-1)
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
    SLL_ASSERT( (x .ge. spline%xmin) .and. (x .le. spline%xmax) )
    SLL_ASSERT( associated(spline) )
    ! FIXME: arg checks here
    num_cells = spline%n_cells
    rh        = spline%rdelta
    coeffs    => spline%c
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
    num_cells = spline%n_cells
    rh        = spline%rdelta
    coeffs    => spline%c
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
    SLL_DEALLOCATE( spline%c, ierr )
    spline%data => null()
    SLL_DEALLOCATE( spline, ierr )
  end subroutine delete_spline_1D
end module sll_splines
