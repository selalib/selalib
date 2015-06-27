!> @ingroup splines
!> Contains bsplines documentation
module sll_bsplines

#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_boundary_condition_descriptors.h"
#include "sll_utilities.h"

implicit none 

private
  

!> @brief 
!> basic type for one-dimensional B-spline data. 
!> @details This should be
!> treated as an opaque type. No access to its internals is directly allowed.
type, public :: sll_bspline_1d

  sll_int32                 :: n
  sll_int32                 :: k
  sll_real64, pointer       :: tau(:)
  sll_real64, pointer       :: t(:)
  sll_real64, pointer       :: q(:)
  sll_real64, pointer       :: bcoef(:)
  sll_int32                 :: bc_type
  sll_real64, pointer       :: dbiatx(:,:)
  sll_real64, dimension(20) :: deltal
  sll_real64, dimension(20) :: deltar
  sll_int32                 :: j = 1
  sll_real64                :: length
  sll_int32                 :: ilo
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

type, public :: sll_bspline_2d

  sll_int32                 :: n
  sll_int32                 :: k
  sll_real64, pointer       :: tau(:)
  sll_real64, pointer       :: t(:)
  sll_real64, pointer       :: q(:)
  sll_real64, pointer       :: bcoef(:)
  sll_int32                 :: bc_type
  sll_real64, pointer       :: dbiatx(:,:)
  sll_real64, pointer       :: a(:,:)
  sll_real64, dimension(20) :: deltal
  sll_real64, dimension(20) :: deltar
  sll_int32                 :: j = 1
  sll_real64                :: length
  sll_int32                 :: ilo

end type sll_bspline_2d

public :: new_bspline_1d
public :: delete_bspline_1d
public :: compute_bspline_1d
public :: update_bspline_1d
public :: interpolate_value
public :: interpolate_derivative
public :: interpolate_array_values
public :: interpolate_array_derivatives

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!> @brief Returns a pointer to a heap-allocated bspline object.
!> @param[in] num_points Number of points where the data to be 
!> interpolated are represented.
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

  SLL_ALLOCATE( new_bspline_1d, ierr )

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

  call initialize_bspline_1d( new_bspline_1d,  &
                              num_points,      &
                              degree+1,        &
                              xmin,            &
                              xmax,            &
                              bc_type          ) 

end function new_bspline_1d

subroutine initialize_bspline_1d(this, n, k, tau_min, tau_max, bc_type)

  type(sll_bspline_1d)              :: this
  sll_int32         , intent(in)    :: n
  sll_int32         , intent(in)    :: k
  sll_real64        , intent(in)    :: tau_min
  sll_real64        , intent(in)    :: tau_max
  sll_int32         , intent(in)    :: bc_type

  sll_int32                         :: i
  sll_int32                         :: ierr
  sll_int32, parameter              :: m=2
  sll_real64                        :: delta_tau

  this%bc_type = bc_type
  this%n       = n
  this%k       = k
  this%length  = tau_max - tau_min
  delta_tau    = this%length / (n-1)

  SLL_ALLOCATE(this%tau(n), ierr)
  do i = 1, n
    this%tau(i) = tau_min + (i-1) * delta_tau
  end do

  if ( bc_type == SLL_PERIODIC) then

    SLL_ALLOCATE(this%t(n+k), ierr)
    SLL_ALLOCATE(this%bcoef(n),  ierr)
    SLL_ALLOCATE(this%q(1:(2*k-1)*n), ierr)

    this%t(1:k)     = tau_min
    if ( mod(k,2) == 0 ) then
      do i = k+1,n
        this%t(i) = this%tau(i-k/2) 
      end do
    else
      do i = k+1, n
        this%t(i) = 0.5*(this%tau(i-(k-1)/2)+this%tau(i-1-(k-1)/2))
      end do
    end if
    this%t(n+1:n+k) = tau_max

  else

    SLL_ALLOCATE(this%t(n+k+m), ierr)
    SLL_ALLOCATE(this%bcoef(n+m),  ierr)
    SLL_ALLOCATE(this%q(1:(2*k-1)*(n+m)), ierr)

    this%t(1:k)         = tau_min
    this%t(k+1:n+m)     = this%tau(2:n-1)
    this%t(n+m+1:n+m+k) = tau_max

  end if

  allocate(this%dbiatx(k,m))
  allocate(this%aj(k))
  allocate(this%dl(k))
  allocate(this%dr(k))

end subroutine initialize_bspline_1d

subroutine build_system(this)

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
  sll_int32               :: ilo
  
  n = this%n
  k = this%k

  kpkm2       = 2*(k-1)
  left        = k
  this%q      = 0.0_f64
  this%dbiatx = 0.0_f64
  
  SLL_ASSERT(m < n) 

  l = 0 ! index for the derivative

  ilo = k

  do i = 1, n
      
    taui = this%tau(i)
    call interv( this%t, n+m+k, taui, left, ilo, mflag )

    if (i < n) then

      call bsplvb ( this, k, 1, taui, left, this%bcoef )
      jj = i-left+1+(left-k)*(k+k-1)+l
      do j = 1, k
        jj = jj + kpkm2
        this%q(jj) = this%bcoef(j)
      end do
   
      if ( i == 1 ) then   
        call bsplvd( this, k, taui, left, 2)
        l = l + 1
        jj = i-left+1+(left-k)*(k+k-1)+l
        do j = 1, k
          jj = jj + kpkm2
          this%q(jj) = this%dbiatx(j,2)
        end do
      end if

    else

      call bsplvd( this, k, taui, left, 2)
      jj = i-left+1+(left-k)*(k+k-1)+l
      do j = 1, k
        jj = jj + kpkm2
        this%q(jj) = this%dbiatx(j,2)
      end do
      l = l + 1
      
      call bsplvb ( this, k, 1, taui, left, this%bcoef )
      jj = i-left+1+(left-k)*(k+k-1)+l
      do j = 1, k
        jj = jj + kpkm2 
        this%q(jj) = this%bcoef(j)
      end do

    end if
 
  end do
  
  !Obtain factorization of A, stored again in Q.

  call banfac ( this%q, k+k-1, n+m, k-1, k-1, iflag )

end subroutine build_system

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
      
    taui = this%tau(i)
    ilp1mx = min ( i + k, n + 1 )
    left = max ( left, i )
    if ( taui < this%t(left) ) stop '  The linear system is not invertible!'

    do while ( this%t(left+1) <= taui )
      
      left = left + 1
      
      if ( left < ilp1mx ) cycle
      
      left = left - 1
      
      if ( this%t(left+1) < taui ) stop '  The linear system is not invertible!'
      
      exit
      
    end do

    call bsplvb ( this, k, 1, taui, left, this%bcoef )
    jj = i-left+1+(left-k)*(k+k-1)
    do j = 1, k
        jj = jj + kpkm2
        this%q(jj) = this%bcoef(j)
    end do
   
  end do
  
  !Obtain factorization of A, stored again in Q.

  call banfac ( this%q, k+k-1, n, k-1, k-1, iflag )

end subroutine build_system_periodic

!> @brief
!>  produces the B-spline coefficients of an interpolating spline.
!> @details
!>   The spline is of order K with knots T(1:N+K), and takes on the 
!>   value GTAU(I) at TAU(I), for I = 1 to N.
!>   The I-th equation of the linear system A * BCOEF = B 
!>   for the B-spline coefficients of the interpolant enforces interpolation
!>   at TAU(1:N).
!>   The matrix A is generated row by row and stored, diagonal by diagonal,
!>   in the rows of the array Q, with the main diagonal going
!>   into row K.  
!>   The banded system is then solved by a call to BANFAC, which 
!>   constructs the triangular factorization for A and stores it again in
!>   Q, followed by a call to BANSLV, which then obtains the solution
!>   BCOEF by substitution.
!>   The B-coefficients for the interpolant 
!>   of an additional data set can be obtained without going through all 
!>   the calculations in this routine, simply by calling the subroutine
!>   update_bspline_1d
!>
!>   @param[inout], the bspline object with knots and data point abscissas.
!>   The triangular factorization of the coefficient matrix of the linear 
!>   system for the B-coefficients is computed.
!>   @param[in]  gtau(N), the data ordinates.
!>   @param[in]  slope_min, the derivative at left boundary.
!>   @param[in]  slope_max, the derivative at right boundary.
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
    call build_system(this)
  end if

  if (present(slope_min) .and. present(slope_max)) then
    call update_bspline_1d( this, gtau, slope_min, slope_max)
  else
    call update_bspline_1d( this, gtau)
  end if

  this%ilo = this%k

end subroutine compute_bspline_1d

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

    if (present(slope_min)) then
      this%bcoef(2)   = slope_min
    else
      if(this%compute_sl) then
        call apply_fd(k+1,1,this%tau(1:k+1),gtau(1:k+1),this%tau(1),slope(0:1))
        this%bcoef(2) = slope(1)
      else
        this%bcoef(2) = this%sl
      end if
    end if
    this%bcoef(3:n) = gtau(2:n-1)
    if (present(slope_max)) then
      this%bcoef(n+1) = slope_max
    else
      if(this%compute_sr) then
        call apply_fd(k+1,1,this%tau(n-k-1:n),gtau(n-k-1:n),this%tau(n),slope(0:1))
        this%bcoef(n+1) = slope(1)
      else
        this%bcoef(n+1) = this%sr
      end if
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
function interpolate_value( this, x) result(y)

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
  
  k = this%k
  n = this%n  

  if (this%bc_type == SLL_PERIODIC) then
    nmk = n+k
  else
    nmk = n+m+k
  end if
  
  call interv( this%t, nmk, x, i, this%ilo, mflag )
  
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
      this%aj(jj) = (this%aj(jj+1)*this%dl(ilo)+this%aj(jj)*this%dr(jj))/(this%dl(ilo)+this%dr(jj))
      ilo = ilo - 1
    end do
  end do
  
  y = this%aj(1)
  
end function interpolate_value

!> @brief returns the values of the derivatives evaluated at a 
!> collection of abscissae stored by a 1D array in another output 
!> array. The spline coefficients
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
subroutine interpolate_array_values( this, n, x, y)

  type(sll_bspline_1d), pointer     :: this 
  sll_int32,            intent(in)  :: n
  sll_real64,           intent(in)  :: x(n)
  sll_real64,           intent(out) :: y(n)
  
  sll_int32               :: i
  sll_int32               :: j
  sll_int32               :: ilo
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
  
  k = this%k
  if (this%bc_type == SLL_PERIODIC) then
    nmk = n+k
  else
    nmk = n+m+k
  end if
  
  SLL_CLEAR_ALLOCATE(aj(1:k), ierr)
  SLL_CLEAR_ALLOCATE(dl(1:k), ierr)
  SLL_CLEAR_ALLOCATE(dr(1:k), ierr)

  ilo = k
  i   = k
  
  do l = 1, n
  
    xi = x(l)
  
    call interv ( this%t, nmk, xi, i, ilo, mflag )
    
    if ( mflag /= 0 ) return
  
    if ( k <= 1 ) then
      y(l) = this%bcoef(i)
      cycle
    end if
    
    if ( k <= i ) then
      do j = 1, k-1
        dl(j) = xi - this%t(i+1-j)
      end do
      jcmin = 1
    else
      jcmin = 1-(i-k)
      do j = 1, i
        dl(j) = xi - this%t(i+1-j)
      end do
      do j = i, k-1
        aj(k-j) = 0.0_f64
        dl(j) = dl(i)
      end do
    end if
    
    if ( nmk-k < i ) then
      jcmax = nmk-i
      do j = 1, jcmax
        dr(j) = this%t(i+j) - xi
      end do
      do j = jcmax, k-1
        aj(j+1) = 0.0_f64
        dr(j) = dr(jcmax)
      end do
    else
      jcmax = k
      do j = 1, k-1
        dr(j) = this%t(i+j) - xi
      end do
    end if
    
    do jc = jcmin, jcmax
      aj(jc) = this%bcoef(i-k+jc)
    end do
    
    do j = 1, k-1
      jlo = k-j
      do jj = 1, k-j
        aj(jj) = (aj(jj+1)*dl(jlo)+aj(jj)*dr(jj))/(dl(jlo)+dr(jj))
        jlo = jlo - 1
      end do
    end do
    
    y(l) = aj(1)
  
  end do
  
  deallocate(aj,dl,dr)

end subroutine interpolate_array_values

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
subroutine interpolate_array_derivatives( this, n, x, y)

  type(sll_bspline_1d)    :: this 
  sll_int32,  intent(in)  :: n
  sll_real64, intent(in)  :: x(n)
  sll_real64, intent(out) :: y(n)
  
  sll_int32               :: i
  sll_int32               :: j
  sll_int32               :: ilo
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
  sll_int32               :: jderiv = 1
  
  k = this%k

  if (this%bc_type == SLL_PERIODIC) then
    nmk = n+k
  else
    nmk = n+m+k
  end if
  
  SLL_ALLOCATE(aj(k), ierr)
  SLL_ALLOCATE(dl(k), ierr)
  SLL_ALLOCATE(dr(k), ierr)
  
  ilo = k
  
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
    call interv ( this%t, nmk, xi, i, ilo, mflag )
    
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

end subroutine interpolate_array_derivatives

function interpolate_derivative( this, x) result(y)

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
  call interv ( this%t, nmk, x, i, this%ilo, mflag )
  
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
      this%aj(jj) = ((this%aj(jj+1)-this%aj(jj))/(this%dl(ilo)+this%dr(jj)))*(k-j)
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
      this%aj(jj) = (this%aj(jj+1)*this%dl(ilo)+this%aj(jj)*this%dr(jj))/(this%dl(ilo)+this%dr(jj))
      ilo = ilo-1
    end do
  end do
  
  y = this%aj(1)
  
end function interpolate_derivative

!>@brief
!> Calculates the nonvanishing B-splines and derivatives at X.
!> @details
!>   Values at X of all the relevant B-splines of order K:K+1-NDERIV
!>   are generated via BSPLVB and stored temporarily in DBIATX.
!>   Then the B-spline coefficients of the required derivatives
!>   of the B-splines of interest are generated by differencing,
!>   each from the preceding one of lower order, and combined with
!>   the values of B-splines of corresponding order in DBIATX
!>   to produce the desired values.
!
!
!> Parameters:
!> @param[in] T(LEFT+K), the knot sequence.  It is assumed that
!>   T(LEFT) < T(LEFT+1).  Also, the output is correct only if
!>   T(LEFT) <= X <= T(LEFT+1).
!> @param[in] K, the order of the B-splines to be evaluated.
!> @param[in] X, the point at which these values are sought.
!> @param[in] LEFT, indicates the left endpoint of the interval of
!>   interest.  The K B-splines whose support contains the interval
!>   ( T(LEFT), T(LEFT+1) ) are to be considered.
!
!>   DBIATX(K,NDERIV).  DBIATX(I,M) contains
!>   the value of the (M-1)st derivative of the (LEFT-K+I)-th B-spline
!>   of order K for knot sequence T, I=M,...,K, M=1,...,NDERIV.
!>   @parama[in] NDERIV, indicates that values of B-splines and their
!>   derivatives up to but not including the NDERIV-th are asked for.
!
!   Reference: 
!   Carl DeBoor,
!   A Practical Guide to Splines,
!   Springer, 2001,
!   ISBN: 0387953663.
subroutine bsplvd ( this, k, x, left, nderiv )
    
  type(sll_bspline_1d)      :: this
  sll_int32,  intent(in)    :: k
  sll_real64, intent(in)    :: x
  sll_int32,  intent(in)    :: left
  sll_int32,  intent(in)    :: nderiv

  sll_real64, allocatable   :: a(:,:)
  sll_real64 :: factor
  sll_real64 :: fkp1mm
  sll_int32  :: i
  sll_int32  :: ideriv
  sll_int32  :: il
  sll_int32  :: j
  sll_int32  :: jlow
  sll_int32  :: jp1mid
  sll_int32  :: ldummy
  sll_int32  :: m
  sll_int32  :: mhigh
  
  allocate(a(k,k))

  mhigh = max ( min ( nderiv, k ), 1 )
  !
  !  MHIGH is usually equal to NDERIV.
  !
  call bsplvb ( this, k+1-mhigh, 1, x, left, this%dbiatx(:,1) )
  
  if ( mhigh == 1 ) return
  !
  !  The first column of DBIATX always contains the B-spline values
  !  for the current order.  These are stored in column K+1-current
  !  order before BSPLVB is called to put values for the next
  !  higher order on top of it.
  !
  ideriv = mhigh
  do m = 2, mhigh
    jp1mid = 1
    do j = ideriv, k
       this%dbiatx(j,ideriv) = this%dbiatx(jp1mid,1)
       jp1mid = jp1mid + 1
    end do
    ideriv = ideriv - 1
    call bsplvb ( this, k+1-ideriv, 2, x, left, this%dbiatx(:,1) )
  end do
  !
  !  At this point, B(LEFT-K+I, K+1-J)(X) is in DBIATX(I,J) for
  !  I=J,...,K and J=1,...,MHIGH ('=' NDERIV).
  !
  !  In particular, the first column of DBIATX is already in final form.
  !
  !  To obtain corresponding derivatives of B-splines in subsequent columns,
  !  generate their B-representation by differencing, then evaluate at X.
  !
  jlow = 1
  do i = 1, k
    a(jlow:k,i) = 0.0D+00
    jlow = i
    a(i,i) = 1.0D+00
  end do
  !
  !  At this point, A(.,J) contains the B-coefficients for the J-th of the
  !  K B-splines of interest here.
  !
  do m = 2, mhigh
    fkp1mm = real(k+1-m, kind = f64 )
    il = left
    i = k
    !
    !  For J = 1,...,K, construct B-coefficients of (M-1)st derivative of
    !  B-splines from those for preceding derivative by differencing
    !  and store again in  A(.,J).  The fact that  A(I,J) = 0 for
    !  I < J is used.
    !
    do ldummy = 1, k+1-m
      factor = fkp1mm / ( this%t(il+k+1-m) - this%t(il) )
      !  The assumption that T(LEFT) < T(LEFT+1) makes denominator
      !  in FACTOR nonzero.
      a(i,1:i) = ( a(i,1:i) - a(i-1,1:i) ) * factor
      il = il - 1
      i = i - 1
    end do
    !  For I = 1,...,K, combine B-coefficients A(.,I) with B-spline values
    !  stored in DBIATX(.,M) to get value of (M-1)st derivative of
    !  I-th B-spline (of interest here) at X, and store in DBIATX(I,M).
    !
    !  Storage of this value over the value of a B-spline
    !  of order M there is safe since the remaining B-spline derivatives
    !  of the same order do not use this value due to the fact
    !  that  A(J,I) = 0  for J < I.
    do i = 1, k
      jlow = max ( i, m )
      this%dbiatx(i,m) = dot_product ( a(jlow:k,i), this%dbiatx(jlow:k,m) )
    end do
  end do

  deallocate(a)

end subroutine bsplvd
  
!***********************************************************************
!> @brief
!> Evaluates B-splines at a point X with a given knot sequence.
!> @details
!>   Evaluates all possibly nonzero B-splines at X of order
!>     JOUT = MAX ( JHIGH, (J+1)*(INDEX-1) )
!>   with knot sequence T. The recurrence relation
!>   \f[
!> B_{i,j+1}(x) = \frac{x-t_i}{t_{i+j}-t_i} B_{i,j}(x) + 
!> \frac{t_{i+j+1}-x}{t_{i+j+1}-t_{i+1}} B_{i+1,j}(x)
!> \f]
!>   is used to generate B(LEFT-J:LEFT,J+1)(X) from B(LEFT-J+1:LEFT,J)(X)
!>   storing the new values in BIATX over the old.
!
!>   The facts that
!>     B(I,1)(X) = 1  if  T(I) <= X < T(I+1)
!>   and that
!>     B(I,J)(X) = 0  unless  T(I) <= X < T(I+J)
!>   are used.
!>    @param[in] T(LEFT+JOUT), the knot sequence.  T is assumed to
!>    be nondecreasing, and also, T(LEFT) must be strictly less than
!>    T(LEFT+1).
!>
!>    @param[in] JHIGH, INDEX, determine the order
!>    JOUT = max ( JHIGH, (J+1)*(INDEX-1) )
!>    of the B-splines whose values at X are to be returned.
!>    INDEX is used to avoid recalculations when several
!>    columns of the triangular array of B-spline values are
!>    needed, for example, in BVALUE or in BSPLVD.
!>    If INDEX = 1, the calculation starts from scratch and the entire
!>    triangular array of B-spline values of orders
!>    1, 2, ...,JHIGH is generated order by order, that is,
!>    column by column.
!>    If INDEX = 2, only the B-spline values of order J+1, J+2, ..., JOUT
!>    are generated, the assumption being that BIATX, J,
!>    DELTAL, DELTAR are, on entry, as they were on exit
!>    at the previous call.  In particular, if JHIGH = 0,
!>    then JOUT = J+1, that is, just the next column of B-spline
!>    values is generated.
!>    Warning: the restriction  JOUT <= JMAX (= 20) is
!>    imposed arbitrarily by the dimension statement for DELTAL
!>    and DELTAR, but is nowhere checked for.
!>
!>    @param[in] x, the point at which the B-splines are to be evaluated.
!>    @param[in] left, an integer chosen so that T(LEFT) <= X <= T(LEFT+1).
!>    @param[out] biatx(jout), with biatx(i) containing the
!>    value at X of the polynomial of order JOUT which agrees
!>    with the B-spline B(LEFT-JOUT+I,JOUT,T) on the interval
!>    (T(LEFT),T(LEFT+1)).
!>

!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663.
!
subroutine bsplvb ( this, jhigh, index, x, left, biatx )
  
  type(sll_bspline_1d)    :: this
  sll_int32,  intent(in)  :: jhigh
  sll_int32,  intent(in)  :: index
  sll_real64, intent(in)  :: x
  sll_int32,  intent(in)  :: left
  sll_real64, intent(out) :: biatx(:) ! (jhigh)
  
  sll_int32               :: i
  sll_real64              :: saved
  sll_real64              :: term
  
  if ( index == 1 ) then
    this%j = 1
    biatx(1) = 1.0_f64
    if ( jhigh <= this%j ) return
  end if
  
  SLL_ASSERT ( this%t(left+1) > this%t(left) )
  
  !if ( this%t(left+1) <= this%t(left) ) then
  !  print*,'x=',x
  !  write ( *, '(a)' ) ' '
  !  write ( *, '(a)' ) 'BSPLVB - Fatal error!'
  !  write ( *, '(a)' ) '  It is required that T(LEFT) < T(LEFT+1).'
  !  write ( *, '(a,i8)' ) '  But LEFT = ', left
  !  write ( *, '(a,g14.6)' ) '  T(LEFT) =   ', this%t(left)
  !  write ( *, '(a,g14.6)' ) '  T(LEFT+1) = ', this%t(left+1)
  !  stop
  !end if
  
  do
     
    this%deltar(this%j) = this%t(left+this%j) - x
    this%deltal(this%j) = x - this%t(left+1-this%j)
    
    saved = 0.0_f64
    do i = 1, this%j
       term = biatx(i) / ( this%deltar(i) + this%deltal(this%j+1-i) )
       biatx(i) = saved + this%deltar(i) * term
       saved = this%deltal(this%j+1-i) * term
    end do
  
    biatx(this%j+1) = saved
    this%j = this%j + 1
    
    if ( jhigh <= this%j ) exit
  
  end do

end subroutine bsplvb
  
!> Brackets a real value in an ascending vector of values.
!> @details
!!
!!  The XT array is a set of increasing values.  The goal of the routine
!!  is to determine the largest index I so that XT(I) <= X.
!!
!!  The routine is designed to be efficient in the common situation
!!  that it is called repeatedly, with X taken from an increasing
!!  or decreasing sequence.
!!
!!  This will happen when a piecewise polynomial is to be graphed.
!!  The first guess for LEFT is therefore taken to be the value
!!  returned at the previous call and stored in the local variable ILO.
!!
!!  A first check ascertains that ILO < LXT.  This is necessary
!!  since the present call may have nothing to do with the previous
!!  call.  Then, if
!!
!!    XT(ILO) <= X < XT(ILO+1),
!!
!!  we set LEFT = ILO and are done after just three comparisons.
!!
!!  Otherwise, we repeatedly double the difference ISTEP = IHI - ILO
!!  while also moving ILO and IHI in the direction of X, until
!!
!!    XT(ILO) <= X < XT(IHI)
!!
!!  after which we use bisection to get, in addition, ILO + 1 = IHI.
!!  The value LEFT = ILO is then returned.
!!
!!  Reference:
!!
!!    Carl DeBoor,
!!    A Practical Guide to Splines,
!!    Springer, 2001,
!!    ISBN: 0387953663.
!!
!! @param[in] xt(lxt), a nondecreasing sequence of values.
!! @param[in] lxt, the dimension of xt.
!! @param[in] x, the point whose location with
!!    respect to the sequence XT is to be determined.
!! @param[out] left, the index of the bracketing value:
!!      1     if             X  <  XT(1)
!!      I     if   XT(I)  <= X  < XT(I+1)
!!      LXT   if  XT(LXT) <= X
!! @param[out] mflag, indicates whether X lies within the
!!    range of the data.
!!    -1:            X  <  XT(1)
!!     0: XT(I)   <= X  < XT(I+1)
!!    +1: XT(LXT) <= X
!<
subroutine interv( xt, lxt, x, left, ilo, mflag )
 
  sll_real64, intent(in)    :: xt(:)
  sll_int32,  intent(in)    :: lxt
  sll_real64, intent(in)    :: x
  sll_int32,  intent(out)   :: left
  sll_int32,  intent(inout) :: ilo
  sll_int32,  intent(out)   :: mflag
  
  sll_int32                 :: ihi
  sll_int32                 :: istep
  sll_int32                 :: middle
  sll_real64                :: xtmax

  xtmax = xt(lxt)
  
  ihi = ilo + 1
  if ( ihi >= lxt ) then
    if ( x >= xtmax ) goto 110
    ilo = lxt - 1
    ihi = lxt
  end if
  if ( xt(ihi) <= x ) goto 20
  if ( xt(ilo) <= x ) then
    mflag = 0
    left = ilo
    return
  end if
  !
  !  Now X < XT(ILO).  Decrease ILO to capture X.
  !
  istep = 1
  10  continue
  ihi = ilo
  ilo = ihi - istep
  if ( 1 < ilo ) then
    if ( xt(ilo) <= x ) goto 50
    istep = istep * 2
    goto 10
  end if
  ilo = 1
  if ( x < xt(1) ) then
    mflag = -1
    left = 1
    return
  end if
  goto 50
  !
  !  Now XT(IHI) <= X.  Increase IHI to capture X.
  !
  20 continue
  istep = 1
  30 continue
  ilo = ihi
  ihi = ilo + istep
  if ( ihi < lxt ) then
    if ( x < xt(ihi) ) goto 50
    istep = istep * 2
    goto 30
  end if
  if ( xtmax <= x ) goto 110
  !
  !  Now XT(ILO) < = X < XT(IHI).  Narrow the interval.
  !
  ihi = lxt
  50  continue
  do
    middle = (ilo+ihi)/2
    if ( middle == ilo ) then
      mflag = 0
      left = ilo
      return
    end if
    !
    !  It is assumed that MIDDLE = ILO in case IHI = ILO+1.
    !
    if ( xt(middle) <= x ) then
      ilo = middle
    else
      ihi = middle
    end if
     
  end do
  !
  !  Set output and return.
  !
  110 continue
  mflag = 1
  if ( x == xtmax ) mflag = 0

  do left = lxt, 1, -1
    if ( xt(left) < xtmax ) return
  end do

end subroutine interv

!> @brief Returns the interpolated value of the derivative in the x1 
!> direction at the point
!> (x1,x2) using the spline decomposition stored in the spline object.
!> @param[in] x1 first coordinate.
!> @param[in] x2 second coordinate.
!> @param[in] spline pointer to spline object.
!> @returns the interpolated value of the derivative in the x1 
function interpolate_x1_derivative_2D( x1, x2, spline )
  sll_real64                          :: interpolate_x1_derivative_2D
  sll_real64, intent(in)              :: x1
  sll_real64, intent(in)              :: x2
  type(sll_bspline_2D), pointer       :: spline

end function interpolate_x1_derivative_2D

  ! interpolate_x2_derivative_2D(): given discrete data f(i,j) that are
  ! described by a 2-dimensional bspline fit s(x1,x2), where the
  ! continuous variables x1 and x2 are within the original limits of i and j
  ! respectively, interpolate_x1_derivative() returns the value of
  !
  !         partial s
  !       -------------
  !         partial x2
  !
  ! evaluated at the point (x1,x2). (Sorry for the ambiguous use of x1)

  !> @brief 
  !> Returns the interpolated value of the derivative 
  !> @details
  !> in the x2 direction at the point(x1,x2) using the spline 
  !> decomposition stored in the spline object.
  !> @param[in] x1 first coordinate.
  !> @param[in] x2 second coordinate.
  !> @param[in] spline pointer to spline object.
  !> @return interpolate_x2_derivative_2D  interpolated value of the derivative 
  function interpolate_x2_derivative_2D( x1, x2, spline )
    sll_real64                          :: interpolate_x2_derivative_2D
    intrinsic                           :: associated, int, real
    sll_real64, intent(in)              :: x1
    sll_real64, intent(in)              :: x2
    type(sll_bspline_2D), pointer  :: spline
  end function interpolate_x2_derivative_2D


  subroutine delete_bspline_1d( spline )
    type(sll_bspline_1d), pointer :: spline
  end subroutine delete_bspline_1d 
  subroutine delete_bspline_2D( spline )
    type(sll_bspline_2D), pointer :: spline
  end subroutine delete_bspline_2D 

end module sll_bsplines
