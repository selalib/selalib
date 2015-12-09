!> @ingroup splines
!> Contains bsplines implementation
module sll_m_bsplines

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

! use F77_deboor, only: &
!   banfac, &
!   banslv

  use sll_m_boundary_condition_descriptors, only: &
    sll_periodic

  use sll_m_deboor_splines_1d, only: &
    bsplvb, &
    bsplvd, &
    deboor_type, &
    interv

  use sll_m_fornberg, only: &
    apply_fd

  implicit none

  public :: &
    compute_bspline_1d, &
    compute_bspline_2d, &
    delete_bspline_1d, &
    interpolate_array_derivatives_1d, &
    interpolate_array_values_1d, &
    interpolate_array_values_2d, &
    interpolate_derivative_1d, &
    interpolate_value_1d, &
    interpolate_value_2d, &
    new_bspline_1d, &
    new_bspline_2d, &
    sll_bspline_1d, &
    sll_bspline_2d, &
    update_bspline_1d

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
!> @brief 
!> basic type for one-dimensional B-spline data. 
!> @details This should be
!> treated as an opaque type. It is better to use it through interpolator
type :: sll_bspline_1d

  sll_int32                 :: n
  sll_int32                 :: k
  sll_real64, pointer       :: tau(:)
  sll_real64, pointer       :: t(:)
  sll_real64, pointer       :: q(:)
  sll_real64, pointer       :: bcoef(:)
  sll_int32                 :: bc_type
  sll_real64, pointer       :: a(:,:)
  sll_real64, pointer       :: dbiatx(:,:)
  sll_real64                :: length
  sll_real64                :: vl
  sll_real64                :: vr
  sll_real64                :: sl
  sll_real64                :: sr
  logical                   :: compute_vl
  logical                   :: compute_vr
  logical                   :: compute_sl
  logical                   :: compute_sr
  sll_real64, allocatable   :: aj(:)
  sll_real64, allocatable   :: dl(:)
  sll_real64, allocatable   :: dr(:)

end type sll_bspline_1d

!> @brief 
!> basic type for two-dimensional B-spline data. 
!> @details 
!> treated as an opaque type. No access to its internals is directly allowed.
type :: sll_bspline_2d

  type(sll_bspline_1d), pointer :: bs1
  type(sll_bspline_1d), pointer :: bs2
  sll_real64,           pointer :: bcoef(:,:)
  sll_real64,           pointer :: x1_min_slopes(:) => null() 
  sll_real64,           pointer :: x1_max_slopes(:) => null()
  sll_real64,           pointer :: x2_min_slopes(:) => null()
  sll_real64,           pointer :: x2_max_slopes(:) => null()

end type sll_bspline_2d

interface compute_bspline_2d
  module procedure compute_bspline_2d_with_constant_slopes
  module procedure compute_bspline_2d_with_variable_slopes
end interface compute_bspline_2d


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!> @brief Returns a pointer to a heap-allocated bspline object.
!> @param[in] num_points Number of points where the data to be 
!> interpolated are represented.
!> @param[in] degree Spline degree 
!> @param[in] xmin Minimum value of the abscissae where the data are meant 
!> to be interpolated.
!> @param[in] xmax Maximum value of the abscissae where the data are meant 
!> to be interpolated.
!> @param[in] bc_type A boundary condition specifier. Must be one of the
!> symbols defined in the SLL_BOUNDARY_CONDITION_DESCRIPTORS module.
!> @param[in] sl OPTIONAL: The value of the slope at xmin, for use in the case
!> of hermite boundary conditions.
!> @param[in] sr OPTIONAL: The value of the slope at xmin, for use in the case
!> of hermite boundary conditions.
!> @return a pointer to a heap-allocated cubic spline object.
function new_bspline_1d( num_points, degree, xmin, xmax, bc_type, sl, sr )

  type(sll_bspline_1D), pointer    :: new_bspline_1d
  sll_int32,  intent(in)           :: num_points
  sll_int32,  intent(in)           :: degree
  sll_real64, intent(in)           :: xmin
  sll_real64, intent(in)           :: xmax
  sll_int32,  intent(in)           :: bc_type
  sll_real64, intent(in), optional :: sl
  sll_real64, intent(in), optional :: sr

  sll_int32                        :: ierr
  sll_int32                        :: i
  sll_int32                        :: n
  sll_int32                        :: k
  sll_int32, parameter             :: m=2
  sll_real64                       :: delta

  SLL_ALLOCATE( new_bspline_1d, ierr )

  k = degree+1
  n = num_points

  SLL_ALLOCATE(new_bspline_1d%a(k,k), ierr)

  new_bspline_1d%bc_type = bc_type
  new_bspline_1d%n       = num_points
  new_bspline_1d%k       = degree+1
  new_bspline_1d%length  = xmax - xmin

  delta = new_bspline_1d%length / (n-1)

  SLL_ALLOCATE(new_bspline_1d%tau(n), ierr)
  do i = 1, n
    new_bspline_1d%tau(i) = xmin + (i-1) * delta
  end do

  if ( bc_type == SLL_PERIODIC) then

    SLL_ALLOCATE(new_bspline_1d%t(n+k), ierr)
    SLL_ALLOCATE(new_bspline_1d%bcoef(n),  ierr)
    SLL_ALLOCATE(new_bspline_1d%q(1:(2*k-1)*n), ierr)

    new_bspline_1d%t(1:k) = xmin
    if ( mod(k,2) == 0 ) then
      do i = k+1,n
        new_bspline_1d%t(i) = new_bspline_1d%tau(i-k/2) 
      end do
    else
      do i = k+1, n
        new_bspline_1d%t(i) = 0.5*(new_bspline_1d%tau(i  -(k-1)/2) &
                                  +new_bspline_1d%tau(i-1-(k-1)/2))
      end do
    end if
    new_bspline_1d%t(n+1:n+k) = xmax
       
  else

    if (present(sl)) then
      new_bspline_1d%sl = sl
      new_bspline_1d%compute_sl = .false.
    else
      new_bspline_1d%compute_sl = .true.
    end if
    if (present(sr)) then
      new_bspline_1d%sr = sr
      new_bspline_1d%compute_sr = .false.
    else
      new_bspline_1d%compute_sr = .true.
    end if

    SLL_ALLOCATE(new_bspline_1d%t(n+k+m),           ierr)
    SLL_ALLOCATE(new_bspline_1d%bcoef(n+m),         ierr)
    SLL_ALLOCATE(new_bspline_1d%q(1:(2*k-1)*(n+m)), ierr)

    new_bspline_1d%t(1:k)   = xmin
    if ( mod(k,2) == 0) then
      new_bspline_1d%t(k+1:n+m) = new_bspline_1d%tau(2:n-1)
    else
      do i = k+1, n+m
        new_bspline_1d%t(i) = 0.5*(new_bspline_1d%tau(i-k)+ &
                                   new_bspline_1d%tau(i-k+1))
      end do
    end if
    new_bspline_1d%t(n+m+1:n+m+k) = xmax

  end if

  allocate(new_bspline_1d%dbiatx(k,m))
  allocate(new_bspline_1d%aj(k))
  allocate(new_bspline_1d%dl(k))
  allocate(new_bspline_1d%dr(k))

end function new_bspline_1d

!> @brief Returns a pointer to a heap-allocated bspline object.
!> @param[in] nx1 Number of points where the data to be interpolated are 
!>            represented.
!> @param[in] degree1 Spline degree 
!> @param[in] x1min Minimum value of the abscissae where the data are meant 
!> to be interpolated.
!> @param[in] x1max Maximum value of the abscissae where the data are meant 
!> to be interpolated.
!> @param[in] bc1 A boundary condition specifier. Must be one of the
!> symbols defined in the SLL_BOUNDARY_CONDITION_DESCRIPTORS module.
!> @param[in] sl1 OPTIONAL: The value of the slope at xmin, for use in the case
!> of hermite boundary conditions.
!> @param[in] sr1 OPTIONAL: The value of the slope at xmin, for use in the case
!> of hermite boundary conditions.
!> @param[in] nx2 Number of points where the data to be interpolated are 
!>            represented.
!> @param[in] degree2 Spline degree 
!> @param[in] x2min Minimum value of the abscissae where the data are meant 
!> to be interpolated.
!> @param[in] x2max Maximum value of the abscissae where the data are meant 
!> to be interpolated.
!> @param[in] bc2 A boundary condition specifier. Must be one of the
!> symbols defined in the sll_m_boundary_condition_descriptors module.
!> @param[in] sl2 OPTIONAL: The value of the slope at xmin, for use in the case
!> of hermite boundary conditions.
!> @param[in] sr2 OPTIONAL: The value of the slope at xmin, for use in the case
!> of hermite boundary conditions.
!> @return a pointer to a heap-allocated cubic spline object.
function new_bspline_2d( nx1, degree1, x1_min, x1_max, bc1, &
                         nx2, degree2, x2_min, x2_max, bc2, &
                         sl1, sr1, sl2, sr2  )

  type(sll_bspline_2d), pointer    :: new_bspline_2d

  sll_int32,  intent(in)           :: nx1
  sll_int32,  intent(in)           :: degree1
  sll_real64, intent(in)           :: x1_min
  sll_real64, intent(in)           :: x1_max
  sll_int32,  intent(in)           :: bc1
  sll_int32,  intent(in)           :: nx2
  sll_int32,  intent(in)           :: degree2
  sll_real64, intent(in)           :: x2_min
  sll_real64, intent(in)           :: x2_max
  sll_int32,  intent(in)           :: bc2
  sll_real64, intent(in), optional :: sl1
  sll_real64, intent(in), optional :: sr1
  sll_real64, intent(in), optional :: sl2
  sll_real64, intent(in), optional :: sr2

  sll_int32                        :: n1
  sll_int32                        :: n2
  sll_int32                        :: ierr

  SLL_ALLOCATE(new_bspline_2d, ierr )
  SLL_ALLOCATE(new_bspline_2d%x1_min_slopes(1:nx2), ierr)
  SLL_ALLOCATE(new_bspline_2d%x1_max_slopes(1:nx2), ierr)
  SLL_ALLOCATE(new_bspline_2d%x2_min_slopes(1:nx1), ierr)
  SLL_ALLOCATE(new_bspline_2d%x2_max_slopes(1:nx1), ierr)

  if (present(sl1) .and. present(sr1)) then
    new_bspline_2d%bs1 => new_bspline_1d(nx1,degree1,x1_min,x1_max,bc1,sl1,sr1)
  else
    new_bspline_2d%bs1 => new_bspline_1d(nx1,degree1,x1_min,x1_max,bc1)
  end if

  if (present(sl1) .and. present(sr1)) then
    new_bspline_2d%bs2 => new_bspline_1d(nx2,degree2,x2_min,x2_max,bc2,sl2,sr2)
  else
    new_bspline_2d%bs2 => new_bspline_1d(nx2,degree2,x2_min,x2_max,bc2)
  end if

  n1 = size(new_bspline_2d%bs1%bcoef)
  n2 = size(new_bspline_2d%bs2%bcoef)
  SLL_CLEAR_ALLOCATE(new_bspline_2d%bcoef(1:n1,1:n2), ierr)

  if (present(sl1)) new_bspline_2d%x1_min_slopes(:) = sl1
  if (present(sr1)) new_bspline_2d%x1_max_slopes(:) = sr1
  if (present(sl2)) new_bspline_2d%x2_min_slopes(:) = sl2
  if (present(sr2)) new_bspline_2d%x2_max_slopes(:) = sr2

end function new_bspline_2d


subroutine build_system_with_derivative(this)

  type(sll_bspline_1d)    :: this 

  sll_real64              :: taui
  sll_int32               :: kpkm2
  sll_int32               :: left
  sll_int32               :: n
  sll_int32               :: k
  sll_int32               :: iflag
  sll_int32               :: mflag
  sll_int32               :: i
  sll_int32               :: j
  sll_int32               :: jj
  sll_int32               :: l
  sll_int32, parameter    :: m=2
  type(deboor_type)       :: db
  
  n = this%n
  k = this%k

  kpkm2       = 2*(k-1)
  left        = k
  this%q      = 0.0_f64
  this%dbiatx = 0.0_f64
  
  SLL_ASSERT(m < n) 

  l = 0 ! line number of the built system matrix

  db%ilo = k

  do i = 1, n
      
    taui = this%tau(i)
    call interv( db, this%t, n+m+k, taui, left, mflag )

    if (i < n) then

      call bsplvb ( db, this%t, k, 1, taui, left, this%bcoef )
      jj = i-left+1+(left-k)*(k+k-1)+l
      do j = 1, k
        jj = jj + kpkm2
        this%q(jj) = this%bcoef(j)
      end do
   
      if ( i == 1 ) then   
        call bsplvd( db, this%t, k, taui, left, this%a, this%dbiatx, 2)
        l = l + 1
        jj = i-left+1+(left-k)*(k+k-1)+l
        do j = 1, k
          jj = jj + kpkm2
          this%q(jj) = this%dbiatx(j,2)
        end do
      end if

    else

      call bsplvd( db, this%t, k, taui, left, this%a, this%dbiatx, 2)
      jj = i-left+1+(left-k)*(k+k-1)+l
      do j = 1, k
        jj = jj + kpkm2
        this%q(jj) = this%dbiatx(j,2)
      end do
      l = l + 1
      
      call bsplvb ( db, this%t, k, 1, taui, left, this%bcoef )
      jj = i-left+1+(left-k)*(k+k-1)+l
      do j = 1, k
        jj = jj + kpkm2 
        this%q(jj) = this%bcoef(j)
      end do

    end if
 
  end do
  
  !Obtain factorization of A, stored again in Q.

  call banfac ( this%q, k+k-1, n+m, k-1, k-1, iflag )

end subroutine build_system_with_derivative

subroutine build_system_periodic(this)

  type(sll_bspline_1d)    :: this 

  sll_real64              :: taui
  sll_int32               :: kpkm2
  sll_int32               :: left
  sll_int32               :: n
  sll_int32               :: k
  sll_int32               :: iflag
  sll_int32               :: i
  sll_int32               :: j
  sll_int32               :: jj
  sll_int32               :: ilp1mx
  type(deboor_type)       :: db
  
  !PN Warning:
  !PN The system built for periodic boundary conditions is wrong
  !PN This needs improvements, it simplifies the splines
  !PN computed in deboor_splines directory.

  n = this%n
  k = this%k

  kpkm2       = 2*(k-1)
  left        = k
  this%q      = 0.0_f64
  this%dbiatx = 0.0_f64
  
  do i = 1, n
    taui   = this%tau(i)
    ilp1mx = min ( i + k, n + 1 )
    left   = max ( left, i )
    if ( taui < this%t(left) ) stop '  The linear system is not invertible!'
    do while ( this%t(left+1) <= taui )
      left = left + 1
      if ( left < ilp1mx ) cycle
      left = left - 1
      if ( this%t(left+1) < taui ) stop ' The linear system is not invertible!'
      exit
    end do
    call bsplvb ( db, this%t, k, 1, taui, left, this%bcoef )
    jj = i-left+1+(left-k)*(k+k-1)
    do j = 1, k
        jj = jj + kpkm2
        this%q(jj) = this%bcoef(j)
    end do
  end do
  
  call banfac ( this%q, k+k-1, n, k-1, k-1, iflag )

end subroutine build_system_periodic

!> @brief
!>  produces the B-spline coefficients of an interpolating spline.
!> @details
!>  The spline is of order K with knots T(1:N+K), and takes on the 
!>  value GTAU(I) at TAU(I), for I = 1 to N.
!>  The I-th equation of the linear system A * BCOEF = B 
!>  for the B-spline coefficients of the interpolant enforces interpolation
!>  at TAU(1:N).
!>  The matrix A is generated row by row and stored, diagonal by diagonal,
!>  in the rows of the array Q, with the main diagonal going
!>  into row K.  
!>  The banded system is then solved by a call to BANFAC, which 
!>  constructs the triangular factorization for A and stores it again in
!>  Q, followed by a call to BANSLV, which then obtains the solution
!>  BCOEF by substitution.
!>  The B-coefficients for the interpolant 
!>  of an additional data set can be obtained without going through all 
!>  the calculations in this routine, simply by calling the subroutine
!>  update_bspline_1d
!>
!>  @param[inout], the bspline object with knots and data point abscissas.
!>  The triangular factorization of the coefficient matrix of the linear 
!>  system for the B-coefficients is computed.
!>  @param[in]  gtau(N), the data ordinates.
!>  @param[in]  slope_min, the derivative at left boundary.
!>  @param[in]  slope_max, the derivative at right boundary.
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
subroutine compute_bspline_1d(this, gtau, slope_min, slope_max)

  type(sll_bspline_1d)   :: this 
  sll_real64, intent(in) :: gtau(:)
  sll_real64, optional   :: slope_min
  sll_real64, optional   :: slope_max

  if ( this%bc_type == SLL_PERIODIC) then
    call build_system_periodic(this)
  else
    call build_system_with_derivative(this)
  end if

  if (present(slope_min) .and. present(slope_max)) then
    call update_bspline_1d( this, gtau, slope_min, slope_max)
  else
    call update_bspline_1d( this, gtau)
  end if

end subroutine compute_bspline_1d

subroutine compute_bspline_2d_with_constant_slopes(this, gtau, &
  sl1_l, sl1_r, sl2_l, sl2_r)

  type(sll_bspline_2d)    :: this 
  sll_real64, intent(in)  :: gtau(:,:)
  sll_real64, optional    :: sl1_l
  sll_real64, optional    :: sl1_r
  sll_real64, optional    :: sl2_l
  sll_real64, optional    :: sl2_r

  if ( this%bs1%bc_type == SLL_PERIODIC) then
    call build_system_periodic(this%bs1)
  else
    call build_system_with_derivative(this%bs1)
  end if

  if ( this%bs2%bc_type == SLL_PERIODIC) then
    call build_system_periodic(this%bs2)
  else
    call build_system_with_derivative(this%bs2)
  end if

  if (present(sl1_l)) this%x1_min_slopes = sl1_l
  if (present(sl1_r)) this%x1_max_slopes = sl1_r
  if (present(sl2_l)) this%x2_min_slopes = sl2_l
  if (present(sl2_r)) this%x2_max_slopes = sl2_r

  call update_bspline_2d(this, gtau)

end subroutine compute_bspline_2d_with_constant_slopes

subroutine compute_bspline_2d_with_variable_slopes(this, gtau, &
  sl1_l, sl1_r, sl2_l, sl2_r)

  type(sll_bspline_2d)    :: this 
  sll_real64, intent(in)  :: gtau(:,:)
  sll_real64, intent(in)  :: sl1_l(:)
  sll_real64, intent(in)  :: sl1_r(:)
  sll_real64, intent(in)  :: sl2_l(:)
  sll_real64, intent(in)  :: sl2_r(:)

  if ( this%bs1%bc_type == SLL_PERIODIC) then
    call build_system_periodic(this%bs1)
  else
    call build_system_with_derivative(this%bs1)
  end if

  if ( this%bs2%bc_type == SLL_PERIODIC) then
    call build_system_periodic(this%bs2)
  else
    call build_system_with_derivative(this%bs2)
  end if

  SLL_ASSERT(size(sl1_l) == this%bs2%n)
  this%x1_min_slopes = sl1_l
  SLL_ASSERT(size(sl1_r) == this%bs2%n)
  this%x1_max_slopes = sl1_r
  SLL_ASSERT(size(sl2_l) == this%bs1%n)
  this%x2_min_slopes = sl2_l
  SLL_ASSERT(size(sl2_r) == this%bs1%n)
  this%x2_max_slopes = sl2_r

  call update_bspline_2d(this, gtau)

end subroutine compute_bspline_2d_with_variable_slopes


!> @brief 
!> update 2 values before computing bsplines coefficients
!> @details
!> If the points positions did not change use this function instead
!> of compute_bspline_2d. You still need to call compute_bspline_2d at
!> the beginning to build the linear system.
subroutine update_bspline_2d(this, gtau, sl1_l, sl1_r, sl2_l, sl2_r)

  type(sll_bspline_2d)    :: this 
  sll_real64, intent(in)  :: gtau(:,:)
  sll_real64, optional    :: sl1_l
  sll_real64, optional    :: sl1_r
  sll_real64, optional    :: sl2_l
  sll_real64, optional    :: sl2_r
  sll_real64, allocatable :: bwork(:,:)
  sll_real64, allocatable :: coeff1(:)
  sll_real64, allocatable :: coeff2(:)

  sll_int32               :: i
  sll_int32               :: j
  sll_int32               :: n1
  sll_int32               :: n2
  sll_int32               :: bc1
  sll_int32               :: bc2
  sll_int32               :: ierr

  sll_int32               :: m1
  
  n1 = this%bs1%n
  n2 = this%bs2%n

  bc1 = this%bs1%bc_type
  bc2 = this%bs2%bc_type

  if ( bc1 == SLL_PERIODIC ) then
    m1 = 0
  else
    m1 = 2
  end if

  SLL_CLEAR_ALLOCATE(bwork(1:n2,1:n1+m1),ierr)

  if (present(sl1_l)) this%x1_min_slopes = sl1_l
  if (present(sl1_r)) this%x1_max_slopes = sl1_r
  if (present(sl2_l)) this%x2_min_slopes = sl2_l
  if (present(sl2_r)) this%x2_max_slopes = sl2_r

  do j = 1, n2
    call update_bspline_1d( this%bs1, gtau(:,j), &
      this%x1_min_slopes(j), this%x1_max_slopes(j))
    bwork(j,:) = this%bs1%bcoef
  end do
  
  SLL_CLEAR_ALLOCATE(coeff1(1:n1+m1), ierr)
  SLL_CLEAR_ALLOCATE(coeff2(1:n1+m1), ierr)

  coeff1(1:n1) = this%x2_min_slopes(:)
  coeff2(1:n1) = this%x2_max_slopes(:)

  do i = 1, n1+m1
    call update_bspline_1d( this%bs2, bwork(:,i), coeff1(i), coeff2(i))
    this%bcoef(i,:) = this%bs2%bcoef(:)
  end do

  deallocate(bwork)

end subroutine update_bspline_2d

!>@brief
!> produces the B-spline coefficients of an interpolating spline.
!> @details
!> If interpolants are already computed and if you change only 
!> the data set and knots and points positions did not change
!> use this routine.
subroutine update_bspline_1d(this, gtau, slope_min, slope_max)

  type(sll_bspline_1d)    :: this 
  sll_real64, intent(in)  :: gtau(:)
  sll_real64, optional    :: slope_min
  sll_real64, optional    :: slope_max


  sll_int32               :: n
  sll_int32               :: k
  sll_int32, parameter    :: m = 2
  sll_real64              :: slope(0:1)

  n = this%n
  k = this%k

  SLL_ASSERT(size(gtau) == n)

  if (this%bc_type == SLL_PERIODIC) then

    this%bcoef = gtau
    call banslv ( this%q, k+k-1, n, k-1, k-1, this%bcoef )

  else

    this%bcoef(1)   = gtau(1)

    if(this%compute_sl) then
      call apply_fd(k+1,1,this%tau(1:k+1),gtau(1:k+1), &
        this%tau(1),slope(0:1))
      this%bcoef(2) = slope(1)
    else if (present(slope_min)) then
      this%bcoef(2) = slope_min
    else
      this%bcoef(2) = this%sl
    end if
    this%bcoef(3:n) = gtau(2:n-1)
    if(this%compute_sr) then
      call apply_fd(k+1,1,this%tau(n-k-1:n),gtau(n-k-1:n), &
        this%tau(n),slope(0:1))
      this%bcoef(n+1) = slope(1)
    else if (present(slope_max)) then
      this%bcoef(n+1) = slope_max
    else
      this%bcoef(n+1) = this%sr
    end if
    this%bcoef(n+2) = gtau(n)

    call banslv ( this%q, k+k-1, n+m, k-1, k-1, this%bcoef )

  end if
  
end subroutine update_bspline_1d


!> @brief returns the value of the image of an abscissae,
!> The spline coefficients
!> used are stored in the spline object pointer.
!> @param[in] x input double-precison element containing the 
!> abscissae to be interpolated.
!> @param[out] y output double-precision element containing the 
!> results of the interpolation.
!> @param[inout] spline the spline object pointer, duly initialized and 
!> already operated on by the compute_bspline_1d() subroutine.
function interpolate_value_1d( this, x) result(y)

  type(sll_bspline_1d)    :: this 
  sll_real64, intent(in)  :: x
  sll_real64              :: y
  
  sll_int32               :: i
  sll_int32               :: j
  sll_int32               :: ilo
  sll_int32               :: jc
  sll_int32               :: jcmax
  sll_int32               :: jcmin
  sll_int32               :: jj
  sll_int32               :: mflag
  sll_int32               :: k
  sll_int32               :: n
  sll_int32               :: nmk
  sll_int32,  parameter   :: m = 2
  type(deboor_type)       :: db
  
  k = this%k
  n = this%n  

  if (this%bc_type == SLL_PERIODIC) then
    nmk = n+k
  else
    nmk = n+m+k
  end if
  
  call interv( db, this%t, nmk, x, i, mflag )
  
  y = this%bcoef(i)
  
  if ( mflag /= 0 ) return
  
  if ( k <= i ) then
    do j = 1, k-1
      this%dl(j) = x - this%t(i+1-j)
    end do
    jcmin = 1
  else
    jcmin = 1-(i-k)
    do j = 1, i
      this%dl(j) = x - this%t(i+1-j)
    end do
    do j = i, k-1
      this%aj(k-j) = 0.0_f64
      this%dl(j) = this%dl(i)
    end do
  end if
  
  if ( nmk-k < i ) then
    jcmax = nmk-i
    do j = 1, jcmax
      this%dr(j) = this%t(i+j) - x
    end do
    do j = jcmax, k-1
      this%aj(j+1) = 0.0_f64
      this%dr(j) = this%dr(jcmax)
    end do
  else
    jcmax = k
    do j = 1, k-1
      this%dr(j) = this%t(i+j) - x
    end do
  end if
  
  do jc = jcmin, jcmax
    this%aj(jc) = this%bcoef(i-k+jc)
  end do
  
  do j = 1, k-1
    ilo = k-j
    do jj = 1, k-j
      this%aj(jj) = (this%aj(jj+1)*this%dl(ilo)+this%aj(jj)*this%dr(jj)) &
                    /(this%dl(ilo)+this%dr(jj))
      ilo = ilo - 1
    end do
  end do
  
  y = this%aj(1)
  
end function interpolate_value_1d

!> @brief returns the values of the images of a collection of 
!> abscissae, represented by a 1D array in another output array. 
!> The spline coefficients used are stored in the spline object pointer.
!> @param[in] x input double-precison element array containing the 
!> abscissae to be interpolated.
!> @param[out] y output double-precision element array containing the 
!> results of the interpolation.
!> @param[in] n the number of elements of the input array which are to be
!> interpolated.
!> @param[inout] spline the spline object pointer, duly initialized and 
!> already operated on by the compute_bspline_1d() subroutine.
subroutine interpolate_array_values_1d( this, n, x, y)

  type(sll_bspline_1d), pointer     :: this 
  sll_int32,            intent(in)  :: n
  sll_real64,           intent(in)  :: x(n)
  sll_real64,           intent(out) :: y(n)
  
  call interpolate_array_derivatives_1d_aux( this, n, x, y, jderiv=0)

end subroutine interpolate_array_values_1d

!> @brief returns the values of the derivatives evaluated at a 
!> collection of abscissae stored by a 1D array in another output 
!> array. The spline coefficients used are stored in the spline 
!> object pointer.
!> @param[in] x input double-precison element array containing the 
!> abscissae to be interpolated.
!> @param[out] y output double-precision element array containing the 
!> results of the interpolation.
!> @param[in] n the number of elements of the input array which are to be
!> interpolated.
!> @param[inout] spline the spline object pointer, duly initialized and 
!> already operated on by the compute_bspline_1d() subroutine.
subroutine interpolate_array_derivatives_1d( this, n, x, y)
  type(sll_bspline_1d)    :: this 
  sll_int32,  intent(in)  :: n
  sll_real64, intent(in)  :: x(n)
  sll_real64, intent(out) :: y(n)

  call interpolate_array_derivatives_1d_aux( this, n, x, y, jderiv=1)

end subroutine interpolate_array_derivatives_1d

subroutine interpolate_array_derivatives_1d_aux( this, n, x, y, jderiv)

  type(sll_bspline_1d)    :: this 
  sll_int32,  intent(in)  :: n
  sll_real64, intent(in)  :: x(n)
  sll_real64, intent(out) :: y(n)
  
  sll_int32               :: i
  sll_int32               :: j
  sll_int32               :: jlo
  sll_int32               :: jc
  sll_int32               :: jcmax
  sll_int32               :: jcmin
  sll_int32               :: jj
  sll_int32               :: mflag
  sll_int32               :: ierr
  sll_int32               :: k
  sll_int32               :: l
  sll_int32               :: nmk
  sll_int32,  parameter   :: m = 2
  sll_real64, allocatable :: aj(:)
  sll_real64, allocatable :: dl(:)
  sll_real64, allocatable :: dr(:)
  sll_real64              :: xi
  sll_int32               :: jderiv
  type(deboor_type)       :: db
  
  k = this%k

  if (this%bc_type == SLL_PERIODIC) then
    nmk = n+k
  else
    nmk = n+m+k
  end if
  
  SLL_ALLOCATE(aj(k), ierr)
  SLL_ALLOCATE(dl(k), ierr)
  SLL_ALLOCATE(dr(k), ierr)
  
  db%ilo = k
  
  do l = 1, n
  
    xi = x(l)
    if ( k <= jderiv ) cycle
    !
    !  Find I so that 1 <= I < N+K and T(I) < T(I+1) and T(I) <= X < T(I+1).
    !
    !  If no such I can be found, X lies outside the support of the
    !  spline F and  BVALUE = 0.  The asymmetry in this choice of I makes F
    !  right continuous.
    !
    call interv ( db, this%t, nmk, xi, i, mflag )
    
    if ( mflag /= 0 ) return
    !
    !  If K = 1 (and JDERIV = 0), BVALUE = BCOEF(I).
    !
    if ( k <= 1 ) then
      y(l) = this%bcoef(i)
      cycle
    end if
    !
    !  Store the K B-spline coefficients relevant for the knot interval
    !  ( T(I),T(I+1) ) in AJ(1),...,AJ(K) and compute DL(J) = X - T(I+1-J),
    !  DR(J) = T(I+J)-X, J=1,...,K-1.  Set any of the AJ not obtainable
    !  from input to zero.
    !
    !  Set any T's not obtainable equal to T(1) or to T(N+K) appropriately.
    !
    jcmin = 1
    
    if ( k <= i ) then
      do j = 1, k-1
        dl(j) = x(l) - this%t(i+1-j)
      end do
    else
      jcmin = 1 - ( i - k )
      do j = 1, i
        dl(j) = xi - this%t(i+1-j)
      end do
      do j = i, k-1
        aj(k-j) = 0.0_f64
        dl(j) = dl(i)
      end do
    end if
    
    jcmax = k
    if ( nmk-k < i ) then
      jcmax = k + n + m - i
      do j = 1, jcmax
        dr(j) = this%t(i+j) - xi
      end do
      do j = jcmax, k-1
        aj(j+1) = 0.0_8
        dr(j) = dr(jcmax)
      end do
    else
      do j = 1, k-1
        dr(j) = this%t(i+j) - xi
      end do
    end if
    
    do jc = jcmin, jcmax
      aj(jc) = this%bcoef(i-k+jc)
    end do
    !  Difference the coefficients JDERIV times.
    do j = 1, jderiv
      jlo = k - j
      do jj = 1, k - j
        aj(jj) = ((aj(jj+1)-aj(jj))/(dl(jlo)+dr(jj)))*(k-j)
        jlo = jlo - 1
      end do
    end do
    !
    !  Compute value at X in (T(I),T(I+1)) of JDERIV-th derivative,
    !  given its relevant B-spline coefficients in AJ(1),...,AJ(K-JDERIV).
    !
    do j = jderiv+1, k-1
      jlo = k-j
      do jj = 1, k-j
        aj(jj) = (aj(jj+1)*dl(jlo)+aj(jj)*dr(jj))/(dl(jlo)+dr(jj))
        jlo = jlo - 1
      end do
    end do
    
    y(l) = aj(1)

  end do

  deallocate(aj,dl,dr)

end subroutine interpolate_array_derivatives_1d_aux

function interpolate_derivative_1d( this, x) result(y)

  type(sll_bspline_1d)    :: this 
  sll_real64, intent(in)  :: x
  sll_real64              :: y
  
  sll_int32               :: i
  sll_int32               :: j
  sll_int32               :: ilo
  sll_int32               :: jc
  sll_int32               :: jcmax
  sll_int32               :: jcmin
  sll_int32               :: jj
  sll_int32               :: mflag
  sll_int32               :: k
  sll_int32               :: n
  sll_int32               :: nmk
  sll_int32,  parameter   :: m = 2
  sll_int32,  parameter   :: jderiv = 1
  type(deboor_type)       :: db
  
  k = this%k
  n = this%n

  if (this%bc_type == SLL_PERIODIC) then
    nmk = n+k
  else
    nmk = n+m+k
  end if
  
  !  Find I so that 1 <= I < N+K and T(I) < T(I+1) and T(I) <= X < T(I+1).
  !
  !  If no such I can be found, X lies outside the support of the
  !  spline F and  BVALUE = 0.  The asymmetry in this choice of I makes F
  !  right continuous.
  !
  call interv ( db, this%t, nmk, x, i, mflag )
  
  y = this%bcoef(i)
  
  if ( k <= jderiv ) return
  
  if ( mflag /= 0 ) return
  !
  !  If K = 1 (and JDERIV = 0), BVALUE = BCOEF(I).
  !
  if ( k <= 1 ) return
  !
  !  Store the K B-spline coefficients relevant for the knot interval
  !  ( T(I),T(I+1) ) in AJ(1),...,AJ(K) and compute DL(J) = X - T(I+1-J),
  !  DR(J) = T(I+J)-X, J=1,...,K-1.  Set any of the AJ not obtainable
  !  from input to zero.
  !
  !  Set any T's not obtainable equal to T(1) or to T(N+K) appropriately.
  !
  jcmin = 1
  
  if ( k <= i ) then
    do j = 1, k-1
      this%dl(j) = x - this%t(i+1-j)
    end do
  else
    jcmin = 1-(i-k)
    do j = 1, i
      this%dl(j) = x - this%t(i+1-j)
    end do
    do j = i, k-1
      this%aj(k-j) = 0.0_f64
      this%dl(j) = this%dl(i)
    end do
  end if
  
  jcmax = k
  if ( nmk-k < i ) then
    jcmax = nmk-i
    do j = 1, jcmax
      this%dr(j) = this%t(i+j) - x
    end do
    do j = jcmax, k-1
      this%aj(j+1) = 0.0_8
      this%dr(j) = this%dr(jcmax)
    end do
  else
    do j = 1, k-1
      this%dr(j) = this%t(i+j) - x
    end do
  end if
  
  do jc = jcmin, jcmax
    this%aj(jc) = this%bcoef(i-k+jc)
  end do
  !  Difference the coefficients JDERIV times.
  do j = 1, jderiv
    ilo = k - j
    do jj = 1, k - j
      this%aj(jj) = ((this%aj(jj+1)-this%aj(jj)) &
                    /(this%dl(ilo)+this%dr(jj)))*(k-j)
      ilo = ilo-1
    end do
  end do
  !
  !  Compute value at X in (T(I),T(I+1)) of JDERIV-th derivative,
  !  given its relevant B-spline coefficients in AJ(1),...,AJ(K-JDERIV).
  !
  do j = jderiv+1, k-1
    ilo = k-j
    do jj = 1, k-j
      this%aj(jj) = (this%aj(jj+1)*this%dl(ilo)+this%aj(jj)*this%dr(jj)) &
        /(this%dl(ilo)+this%dr(jj))
      ilo = ilo-1
    end do
  end do
  
  y = this%aj(1)
  
end function interpolate_derivative_1d

subroutine delete_bspline_1d( spline )
  type(sll_bspline_1d), pointer :: spline
  print*, associated(spline)
end subroutine delete_bspline_1d 

subroutine delete_bspline_2D( spline )
  type(sll_bspline_2D), pointer :: spline
  print*, associated(spline)
end subroutine delete_bspline_2D 

subroutine interpolate_array_values_2d(this, n1, n2, x, y, ideriv, jderiv)

type(sll_bspline_2d)    :: this
sll_int32               :: n1
sll_int32               :: n2
sll_real64, intent(in)  :: x(:,:)
sll_real64, intent(out) :: y(:,:)
sll_int32               :: ideriv
sll_int32               :: jderiv

sll_int32               :: i
sll_int32               :: j, jj
sll_int32               :: jc, jcmin, jcmax

sll_real64, allocatable :: ajx(:), ajy(:)
sll_real64, allocatable :: dlx(:), dly(:)
sll_real64, allocatable :: drx(:), dry(:)
sll_real64, allocatable :: wrk(:)

sll_int32               :: nx, kx, ny, ky
sll_int32               :: left, leftx, lefty
sll_int32               :: jlo
sll_int32               :: klo
sll_int32               :: llo
sll_int32               :: mflag
sll_int32               :: ierr
sll_int32               :: jjj
sll_int32               :: kkk
sll_int32               :: nmkx
sll_int32               :: nmky

sll_real64              :: xi
sll_real64              :: xj
sll_real64, pointer     :: tx(:)
sll_real64, pointer     :: ty(:)
type(deboor_type)       :: db

nx   =  this%bs1%n
ny   =  this%bs2%n
SLL_ASSERT(n1 <= nx)
SLL_ASSERT(n2 <= ny)
SLL_ASSERT(n1 == size(x,1))
SLL_ASSERT(n2 == size(x,2))
kx   =  this%bs1%k
ky   =  this%bs2%k
tx   => this%bs1%t
ty   => this%bs2%t

if (this%bs1%bc_type == SLL_PERIODIC) then
  nmkx = nx+kx
else
  nmkx = nx+kx+2
end if
if (this%bs2%bc_type == SLL_PERIODIC) then
  nmky = ny+ky
else
  nmky = ny+ky+2
end if

SLL_CLEAR_ALLOCATE(ajx(1:kx),ierr)
SLL_CLEAR_ALLOCATE(dlx(1:kx),ierr)
SLL_CLEAR_ALLOCATE(drx(1:kx),ierr)
SLL_CLEAR_ALLOCATE(ajy(1:ky),ierr)
SLL_CLEAR_ALLOCATE(dly(1:ky),ierr)
SLL_CLEAR_ALLOCATE(dry(1:ky),ierr)
SLL_CLEAR_ALLOCATE(wrk(1:nmkx),ierr)

jlo = ky
do j=1,n2
  db%ilo = kx
  klo = jlo
  xj  = this%bs2%tau(j)
  call interv(db,ty,nmky,xj,lefty,mflag)
  do i=1,nx
    xi = this%bs1%tau(i)
    call interv(db,tx,nmkx,xi,leftx,mflag)
    do jj=1,ky
      jcmin = 1
      if ( kx <= leftx ) then
        do jjj = 1, kx-1
          dlx(jjj) = xi - tx(leftx+1-jjj)
        end do
      else
        jcmin = 1-(leftx-kx)
        do jjj = 1, leftx
          dlx(jjj) = xi - tx(leftx+1-jjj)
        end do
        do jjj = leftx, kx-1
          ajx(kx-jjj) = 0.0_f64
          dlx(jjj) = dlx(leftx)
        end do
      end if
      jcmax = kx
      if ( nmkx-kx < leftx ) then
        jcmax = nmkx-leftx
        do jjj = 1, nmkx-leftx
          drx(jjj) = tx(leftx+jjj) - xi
        end do
        do jjj = nmkx-leftx, kx-1
          ajx(jjj+1) = 0.0_f64
          drx(jjj) = drx(nmkx-leftx)
        end do
      else
        do jjj = 1, kx-1
          drx(jjj) = tx(leftx+jjj) - xi
        end do
      end if
      do jc = jcmin, jcmax
        ajx(jc) = this%bcoef(leftx-kx+jc,lefty-ky+jj)
      end do
      do jjj = 1, ideriv
        llo = kx - jjj
        do kkk = 1, kx - jjj
          ajx(kkk) = ((ajx(kkk+1)-ajx(kkk))/(dlx(llo)+drx(kkk)))*(kx-jjj)
          llo = llo-1
        end do
      end do
      do jjj = ideriv+1, kx-1
        llo = kx-jjj
        do kkk = 1, kx-jjj
          ajx(kkk) = (ajx(kkk+1)*dlx(llo)+ajx(kkk)*drx(kkk)) &
                     /(dlx(llo)+drx(kkk))
          llo = llo - 1
        end do
      end do
      wrk(jj) = ajx(1)
    end do
    call interv(db,ty(lefty-ky+1:nmky),ky+ky,xj,left,mflag)
    jcmin = 1
    if ( ky <= left ) then
      do jjj = 1, ky-1
        dly(jjj) = xj - ty(lefty-ky+left+1-jjj)
      end do
    else
      jcmin = 1-(left-ky)
      do jjj = 1, left
        dly(jjj) = xj - ty(lefty-ky+left+1-jjj)
      end do
      do jjj = left, ky-1
        ajy(ky-jjj) = 0.0_f64
        dly(jjj) = dly(left)
      end do
    end if
    jcmax = ky
    if ( ky < left ) then
      jcmax = ky+ky-left
      do jjj = 1, ky+ky-left
        dry(jjj) = ty(lefty-ky+left+jjj) - xj
      end do
      do jjj = ky+ky-left, ky-1
        ajy(jjj+1) = 0.0_f64
        dry(jjj) = dry(ky+ky-left)
      end do
    else
      do jjj = 1, ky-1
        dry(jjj) = ty(lefty-ky+left+jjj) - xj
      end do
    end if
    do jc = jcmin, jcmax
      ajy(jc) = wrk(left-ky+jc)
    end do
    do jjj = 1, jderiv
      llo = ky - jjj
      do kkk = 1, ky - jjj
        ajy(kkk) = ((ajy(kkk+1)-ajy(kkk))/(dly(llo)+dry(kkk)))*(ky-jjj)
        llo = llo-1
      end do
    end do
    do jjj = jderiv+1, ky-1
      llo = ky-jjj
      do kkk = 1, ky-jjj
        ajy(kkk) = (ajy(kkk+1)*dly(llo)+ajy(kkk)*dry(kkk))/(dly(llo)+dry(kkk))
        llo = llo - 1
      end do
    end do
    y(i,j) = ajy(1)
  end do
end do

deallocate(ajx)
deallocate(dlx)
deallocate(drx)
deallocate(ajy)
deallocate(dly)
deallocate(dry)
deallocate(wrk)

end subroutine interpolate_array_values_2d

function interpolate_value_2d(this, xi, xj, ideriv, jderiv ) result (y)

type(sll_bspline_2d)    :: this
sll_real64, intent(in)  :: xi
sll_real64, intent(in)  :: xj
sll_int32,  intent(in)  :: ideriv
sll_int32,  intent(in)  :: jderiv
sll_real64              :: y

sll_int32               :: jj
sll_int32               :: jc, jcmin, jcmax
sll_int32               :: nx, kx, ny, ky
sll_int32               :: left, leftx, lefty
sll_int32               :: llo
sll_int32               :: mflag
sll_int32               :: jjj
sll_int32               :: kkk
sll_int32               :: nmkx
sll_int32               :: nmky

sll_real64, pointer     :: tx(:)
sll_real64, pointer     :: ty(:)

sll_real64, allocatable :: work(:)
type(deboor_type)       :: db

nx   =  this%bs1%n
ny   =  this%bs2%n
kx   =  this%bs1%k
ky   =  this%bs2%k
tx   => this%bs1%t
ty   => this%bs2%t

allocate(work(size(this%bs1%bcoef)))
work = 0.0_f64

if (this%bs1%bc_type == SLL_PERIODIC) then
  nmkx = nx+kx
else
  nmkx = nx+kx+2
end if
if (this%bs1%bc_type == SLL_PERIODIC) then
  nmky = ny+ky
else
  nmky = ny+ky+2
end if

call interv(db,tx,nmkx,xi,leftx,mflag)
call interv(db,ty,nmky,xj,lefty,mflag)

do jj=1,ky
  jcmin = 1
  if ( kx <= leftx ) then
    do jjj = 1, kx-1
      this%bs1%dl(jjj) = xi - tx(leftx+1-jjj)
    end do
  else
    jcmin = 1-(leftx-kx)
    do jjj = 1, leftx
      this%bs1%dl(jjj) = xi - tx(leftx+1-jjj)
    end do
    do jjj = leftx, kx-1
      this%bs1%aj(kx-jjj) = 0.0_f64
      this%bs1%dl(jjj) = this%bs1%dl(leftx)
    end do
  end if
  jcmax = kx
  if ( nmkx-kx < leftx ) then
    jcmax = nmkx-leftx
    do jjj = 1, nmkx-leftx
      this%bs1%dr(jjj) = tx(leftx+jjj) - xi
    end do
    do jjj = nmkx-leftx, kx-1
      this%bs1%aj(jjj+1) = 0.0_f64
      this%bs1%dr(jjj) = this%bs1%dr(nmkx-leftx)
    end do
  else
    do jjj = 1, kx-1
      this%bs1%dr(jjj) = tx(leftx+jjj) - xi
    end do
  end if
  do jc = jcmin, jcmax
    this%bs1%aj(jc) = this%bcoef(leftx-kx+jc,lefty-ky+jj)
  end do
  do jjj = 1, ideriv
    llo = kx - jjj
    do kkk = 1, kx - jjj
      this%bs1%aj(kkk) = ((this%bs1%aj(kkk+1)-this%bs1%aj(kkk)) &
        /(this%bs1%dl(llo)+this%bs1%dr(kkk)))*(kx-jjj)
      llo = llo-1
    end do
  end do
  do jjj = ideriv+1, kx-1
    llo = kx-jjj
    do kkk = 1, kx-jjj
      this%bs1%aj(kkk) = (this%bs1%aj(kkk+1)*this%bs1%dl(llo)+ &
                          this%bs1%aj(kkk  )*this%bs1%dr(kkk))/  &
                         (this%bs1%dl(llo  )+this%bs1%dr(kkk))
      llo = llo - 1
    end do
  end do
  work(jj) = this%bs1%aj(1)
end do

!klo = this%bs2%ilo
call interv(db,ty(lefty-ky+1:nmky),ky+ky,xj,left,mflag)

jcmin = 1
if ( ky <= left ) then
  do jjj = 1, ky-1
    this%bs2%dl(jjj) = xj - ty(lefty-ky+left+1-jjj)
  end do
else
  jcmin = 1-(left-ky)
  do jjj = 1, left
    this%bs2%dl(jjj) = xj - ty(lefty-ky+left+1-jjj)
  end do
  do jjj = left, ky-1
    this%bs2%aj(ky-jjj) = 0.0_f64
    this%bs2%dl(jjj) = this%bs2%dl(left)
  end do
end if
jcmax = ky
if ( ky < left ) then
  jcmax = ky+ky-left
  do jjj = 1, ky+ky-left
    this%bs2%dr(jjj) = ty(lefty-ky+left+jjj) - xj
  end do
  do jjj = ky+ky-left, ky-1
    this%bs2%aj(jjj+1) = 0.0_f64
    this%bs2%dr(jjj) = this%bs2%dr(ky+ky-left)
  end do
else
  do jjj = 1, ky-1
    this%bs2%dr(jjj) = ty(lefty-ky+left+jjj) - xj
  end do
end if
do jc = jcmin, jcmax
  this%bs2%aj(jc) = work(left-ky+jc)
end do
do jjj = 1, jderiv
  llo = ky - jjj
  do kkk = 1, ky - jjj
    this%bs2%aj(kkk) = ((this%bs2%aj(kkk+1)-this%bs2%aj(kkk)) &
      /(this%bs2%dl(llo)+this%bs2%dr(kkk)))*(ky-jjj)
    llo = llo-1
  end do
end do
do jjj = jderiv+1, ky-1
  llo = ky-jjj
  do kkk = 1, ky-jjj
    this%bs2%aj(kkk) = (this%bs2%aj(kkk+1)*this%bs2%dl(llo)+ &
                        this%bs2%aj(kkk)*this%bs2%dr(kkk))/  &
                       (this%bs2%dl(llo)+this%bs2%dr(kkk))
    llo = llo - 1
  end do
end do
y = this%bs2%aj(1)

end function interpolate_value_2d

subroutine interpolate_array_x1_derivatives_2d(this, n1, n2, x, y)

type(sll_bspline_2d)    :: this
sll_int32,  intent(in)  :: n1
sll_int32,  intent(in)  :: n2
sll_real64, intent(in)  :: x(:,:)
sll_real64, intent(out) :: y(:,:)

SLL_ASSERT(this%bs1%n>0)
SLL_ASSERT(n1 == size(x,1))
SLL_ASSERT(n2 == size(y,2))
y = 0.0_f64
call interpolate_array_values_2d(this, n1, n2, x, y, 1, 0)

end subroutine interpolate_array_x1_derivatives_2d

subroutine interpolate_array_x2_derivatives_2d(this, n1, n2, x, y)

type(sll_bspline_2d)    :: this
sll_int32,  intent(in)  :: n1
sll_int32,  intent(in)  :: n2
sll_real64, intent(in)  :: x(:,:)
sll_real64, intent(out) :: y(:,:)

SLL_ASSERT(this%bs1%n>0)
SLL_ASSERT(n1 == size(x,1))
SLL_ASSERT(n2 == size(y,2))
y = 0.0_f64
call interpolate_array_values_2d(this, n1, n2, x, y, 0, 1)

end subroutine interpolate_array_x2_derivatives_2d

end module sll_m_bsplines
