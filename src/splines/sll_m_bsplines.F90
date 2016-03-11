!> @ingroup splines
!> Implements arbitrary degree bspline implementation at knot averages (Greville points)
!> given a B-Spline object from sll_m_arbitrary_degree_splines
module sll_m_bspline_interpolation

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_boundary_condition_descriptors, only: &
       sll_p_periodic, &
       sll_p_hermite, &
       sll_p_greville

  use sll_m_arbitrary_degree_splines, only: &
       sll_t_arbitrary_degree_spline_1d, &
       sll_s_arbitrary_degree_spline_1d_init, &
       sll_s_arbitrary_degree_spline_1d_free, &
       sll_f_find_cell, &
       sll_s_splines_at_x, &
       sll_s_spline_derivatives_at_x, &
       sll_s_splines_and_n_derivs_at_x, &
       sll_s_uniform_b_splines_at_x

  use schur_complement

!ES  use sll_m_fornberg, only: &
!ES    sll_s_apply_fd

  implicit none

  public :: &
       sll_t_bspline_interpolation_1d, &
       sll_s_bspline_interpolation_1d_init, &
       sll_s_compute_bspline_1d, &
       sll_f_interpolate_value_1d, &
       sll_s_interpolate_array_values_1d, &
       sll_f_interpolate_derivative_1d, &
       sll_s_interpolate_array_derivatives_1d, &
       sll_s_bspline_interpolation_1d_free, &
       sll_t_bspline_interpolation_2d, &
       sll_s_bspline_interpolation_2d_init !, &
!    sll_o_compute_bspline_2d, &
!    sll_f_interpolate_value_2d, &
!    sll_s_interpolate_array_values_2d

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
!> @brief 
!> basic type for one-dimensional B-spline data. 
!> @details This should be
!> treated as an opaque type. It is better to use it through interpolator
type :: sll_t_bspline_interpolation_1d
  type (sll_t_arbitrary_degree_spline_1d) :: bsp
  sll_int32                 :: n        ! dimension of spline space
  sll_int32                 :: deg      ! degree of spline
  sll_real64, pointer       :: tau(:)   ! n interpolation points
  sll_real64, pointer       :: bcoef(:) ! bspline coefficients
  sll_int32                 :: bc_type
  sll_real64, pointer       :: q(:,:)   ! triangular factorization of coefficient matrix
  sll_real64, pointer       :: values(:) 
  sll_real64, pointer       :: bsdx(:,:)
  sll_real64                :: length
  sll_int32                 :: offset
  type(schur_complement_solver)  :: schur
!!$  sll_real64                :: vl
!!$  sll_real64                :: vr
!!$  sll_real64                :: sl
!!$  sll_real64                :: sr
!!$  logical                   :: compute_vl
!!$  logical                   :: compute_vr
!!$  logical                   :: compute_sl
!!$  logical                   :: compute_sr
!!$  sll_real64, allocatable   :: aj(:)
!!$  sll_real64, allocatable   :: dl(:)
!!$  sll_real64, allocatable   :: dr(:)

end type sll_t_bspline_interpolation_1d

!> @brief 
!> basic type for two-dimensional B-spline data. 
!> @details 
!> treated as an opaque type. No access to its internals is directly allowed.
type :: sll_t_bspline_interpolation_2d

  type(sll_t_bspline_interpolation_1d) :: bs1
  type(sll_t_bspline_interpolation_1d) :: bs2
  sll_real64,           pointer :: bcoef(:,:)

end type sll_t_bspline_interpolation_2d

!!$interface sll_o_compute_bspline_2d
!!$  module procedure compute_bspline_2d_with_constant_slopes
!!$  module procedure compute_bspline_2d_with_variable_slopes
!!$end interface sll_o_compute_bspline_2d


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!> @brief Constructor for sll_t_bspline_interpolation_1d object
!> @param[in] num_points Number of points where the data to be 
!> interpolated are represented.
!> @param[in] degree Spline degree 
!> @param[in] xmin Minimum value of the abscissae where the data are meant 
!> to be interpolated.
!> @param[in] xmax Maximum value of the abscissae where the data are meant 
!> to be interpolated.
!> @param[in] bc_type A boundary condition specifier. Can be either
!> sll_p_periodic for periodic splines or sll_p_open for open knots  
  subroutine sll_s_bspline_interpolation_1d_init( self, num_points, degree, xmin, xmax, bc_type)
    type(sll_t_bspline_interpolation_1d) :: self
    sll_int32,  intent(in)               :: num_points
    sll_int32,  intent(in)               :: degree
    sll_real64, intent(in)               :: xmin
    sll_real64, intent(in)               :: xmax
    sll_int32,  intent(in)               :: bc_type

    ! local variables
    sll_real64, dimension(num_points)    :: grid
    sll_int32                            :: ierr
    sll_int32                            :: i
    sll_int32                            :: j
    sll_int32                            :: iflag
    !  sll_int32                            :: n
    !  sll_int32                            :: k
    !  sll_int32, parameter                 :: m=0 !ES m=2
    sll_real64                           :: delta

    ! set first attributes
    self%bc_type = bc_type
    if (bc_type == sll_p_periodic) then
       self%n = num_points - 1 ! dimension of periodic spline space
       self%offset = degree/2  ! offset needed for periodic spline evaluation
    else
       self%n = num_points + degree - 1 ! dimension of non periodic spline space
       self%offset = 0  
    end if
    self%deg       = degree
    self%length  = xmax - xmin
    ! define grid for knots
    delta = self%length / (num_points-1)
    do i=1,num_points-1
       grid(i) = xmin + (i-1) * delta
    end do
    grid(num_points) = xmax
    ! construct a sll_t_arbitrary_degree_spline_1d object
    call sll_s_arbitrary_degree_spline_1d_init(self%bsp, degree, grid, num_points, bc_type )

    SLL_ALLOCATE(self%bcoef(self%n),  ierr)
    SLL_ALLOCATE(self%values(self%deg+1), ierr)
    SLL_ALLOCATE(self%bsdx(self%deg+1,self%deg+1), ierr)
    
    ! Allocate tau which contains interpolation points
    SLL_ALLOCATE(self%tau(self%n), ierr)
    select case (bc_type)
    case (sll_p_periodic)
       if (modulo(degree,2) == 0) then
          ! for even degree interpolation points are cell midpoints
          do i = 1, self%n
             self%tau(i) = xmin + (i-0.5_f64) * delta
          end do
       else
          ! for odd degree interpolation points are grid points excluding last point
          do i = 1, self%n
             self%tau(i) = xmin + (i-1) * delta
          end do
       end if
    case (sll_p_hermite)
       if (modulo(degree,2) == 0) then
          ! for even degree interpolation points are cell midpoints
          do i = 1, num_points - 1
             self%tau(i) = xmin + (i-0.5_f64) * delta
          end do
       else
          ! for odd degree interpolation points are grid points including last point
          do i = 1, num_points 
             self%tau(i) = xmin + (i-1) * delta
          end do
       end if
    case (sll_p_greville)
       do i = 1, self%n
          ! Set interpolation points to be Greville points 
          self%tau(i) = 0.0_f64
          do j=1-degree,0
             self%tau(i) = self%tau(i) + self%bsp%knots(i+j)
          end do
          self%tau(i) = self%tau(i) / degree
       end do
    end select
    !print*, 'tau ', self%tau
    !print*, 'knots', self%bsp%knots
    ! Assemble banded matrix (B_j(tau(i))) for spline interpolation
    select case (bc_type)
    case (sll_p_periodic) 
       call build_system_periodic(self) 
    case (sll_p_greville)
       call build_system_greville(self)
    case (sll_p_hermite)
       call build_system_with_derivative(self)
    end select
    

    !  allocate(self%dbiatx(k,m))
    !  allocate(self%aj(k))
    !  allocate(self%dl(k))
    !  allocate(self%dr(k))

  end subroutine sll_s_bspline_interpolation_1d_init

  !> @brief private subroutine for assembling and factorizing 
  !> linear system needed for periodic spline interpolation
  !> @param[inout] self bspline interpolation object
  subroutine build_system_periodic(self)
    type(sll_t_bspline_interpolation_1d)   :: self
    ! local variables
    sll_int32                              :: i
    sll_int32                              :: j
    sll_int32                              :: ierr
    sll_int32                              :: sizeq

    ! evaluate bsplines at interpolation points
    if (modulo(self%deg,2) == 0) then
       ! for even degree interpolation points are cell midpoints
       call sll_s_uniform_b_splines_at_x(self%deg, 0.5_f64, self%values)
    else
       ! for odd degree interpolation points are grid points
       call sll_s_uniform_b_splines_at_x(self%deg, 0.0_f64, self%values)
    end if
    sizeq = self%deg/2
    ! assemble matrix q for linear system
    SLL_ALLOCATE(self%q(2*sizeq+1,self%n), ierr)
    do j=1,self%n
       do i=1, 2*sizeq+1
          self%q(i,j) = self%values(i)
       end do
    end do

    ! Perform LU decomposition of matrix q
    call schur_complement_fac(self%schur, self%n, sizeq, self%q)

  end subroutine build_system_periodic

  !> @brief private subroutine for assembling and factorizing 
  !> linear system needed for spline interpolation at Greville points
  !> @param[inout] self bspline interpolation object
  subroutine build_system_greville(self)
    type(sll_t_bspline_interpolation_1d)   :: self
    ! local variables
    sll_int32                            :: iflag
    sll_int32                            :: i
    sll_int32                            :: j
    sll_int32                            :: k
    sll_int32                            :: ib
    sll_int32                            :: icell
    sll_real64                           :: x
    
    k = self%deg
    SLL_ALLOCATE(self%q(2*k+1,self%n), iflag)

    do i=1,self%n
       x = self%tau(i)
       icell = sll_f_find_cell(self%bsp, x )
       call sll_s_splines_at_x(self%bsp, icell, x, self%values)
       do j=1,k+1
!          print*, i, icell, j,  icell+j-1, icell+j-i+k , x, self%values(j)
          ib = icell+j-1
          self%q( i-icell-j+k+2,ib) = self%values(j)
       end do
    end do
    ! Perform LU decomposition of matrix q
    call banfac ( self%q, 2*k+1, self%n, k, k, iflag )

  end subroutine build_system_greville
  
  !> @brief private subroutine for assembling and factorizing 
  !> linear system needed for spline interpolation with Hermite boundary conditions
  !> @param[inout] self bspline interpolation object
  !> @param[in] nder number of derivatives to be computed
  subroutine build_system_with_derivative(self)
    type(sll_t_bspline_interpolation_1d)   :: self

    ! local variables
    sll_int32                              :: nder
    sll_int32                              :: iflag
    sll_int32                              :: i
    sll_int32                              :: j
    sll_int32                              :: ib
    sll_int32                              :: offset
    sll_int32                              :: mid
    sll_int32                              :: k
    sll_int32                              :: icell
    sll_real64                             :: x

    k = self%deg
    
    ! number of derivative values needed on boundary depending on spline degree
    nder = k/2

    SLL_ALLOCATE(self%q(2*(k+nder)+1,self%n), iflag)
    
    ! For even degree splines interpolation points are at cell midpoints
    ! value of function on the boundary needed as additional boundary conditions
    ! for odd degree splines  baundary is in the interpolation points
    ! only derivative values are needed as boundary conditions
    if (modulo(k,2)==0) then ! spline degree even
       offset = 0
       mid = k
    else
       offset = 1
       mid = k + 1
    end if
   
    ! boundary conditions at xmin
    x = self%bsp%xmin
    icell = sll_f_find_cell(self%bsp, x )
    call sll_s_splines_and_n_derivs_at_x( self%bsp, icell, x , nder , self%bsdx)
    do i=1, nder
       do j=1,k+1
          ib = icell+j-1
          !  i-icell-j+k+2
          self%q(i-ib+mid,ib) = self%bsdx(i+offset,j)
       end do
    end do
    ! interpolation points
    do i=nder+1,self%n - nder
       x = self%tau(i-nder)
       icell = sll_f_find_cell(self%bsp, x )
       call sll_s_splines_at_x(self%bsp, icell, x, self%values)
       do j=1,k+1
          ib = icell+j-1 
          print*, i,x,icell, i-nder -ib+mid, ib
          self%q(i - nder -ib + mid,ib) = self%values(j)
       end do
    end do
    ! boundary conditions at xmax
    x = self%bsp%xmax
    icell = sll_f_find_cell(self%bsp, x )
    call sll_s_splines_and_n_derivs_at_x( self%bsp, icell, x , nder , self%bsdx)
    do i= 1, nder
       do j=1,k+1
          ib = icell+j-1
          self%q( self%n-nder+i-ib+mid,ib) = self%bsdx(i+offset,j)
       end do
    end do
    do i=1,2*k+1
       print*,self%q(i,:)
    end do
    ! Perform LU decomposition of matrix q
    call banfac ( self%q, 2*k+1, self%n, k, k, iflag )

  end subroutine build_system_with_derivative
  
  subroutine sll_s_compute_bspline_1d(self, gtau, val_min, val_max)

    type(sll_t_bspline_interpolation_1d)   :: self
    sll_real64, intent(in) :: gtau(:)
    sll_real64, optional   :: val_min(:)
    sll_real64, optional   :: val_max(:)
    ! local variables
    sll_int32   :: ncond
    
    select case (self%bc_type)
    case(sll_p_periodic)
       self%bcoef = gtau
       call schur_complement_slv ( self%schur, self%n, self%deg/2, self%q,  self%bcoef )
    case (sll_p_greville)
       self%bcoef = gtau
       call banslv ( self%q, 2*self%deg+1, self%n, self%deg, self%deg, self%bcoef )
    case (sll_p_hermite)
       ! number of needed conditions at boundary
       ncond = self%deg/2
       if (present(val_min)) then
          self%bcoef(1:ncond) = val_min(1:ncond)
       else  ! set needed boundary values to 0
          self%bcoef(1:ncond) = 0.0_f64
       end if
       self%bcoef(ncond+1:self%n-ncond) = gtau
       if (present(val_max)) then
          self%bcoef(self%n-ncond+1:self%n) = val_max(1:ncond)
       else ! set needed boundary values to 0
          self%bcoef(self%n-ncond+1:self%n) = 0.0_f64
       end if
       call banslv ( self%q, 2*self%deg+1, self%n, self%deg, self%deg, self%bcoef )
    end select
    !  call sll_s_update_bspline_1d( self, gtau)

  end subroutine sll_s_compute_bspline_1d

!> @brief returns the value of the image of an abscissae,
!> The spline coefficients
!> used are stored in the spline object pointer.
!> @param[in] x input double-precison element containing the 
!> abscissae to be interpolated.
!> @param[out] y output double-precision element containing the 
!> results of the interpolation.
!> @param[inout] spline the spline object pointer, duly initialized and 
!> already operated on by the sll_s_compute_bspline_1d() subroutine.
  function sll_f_interpolate_value_1d( self, x) result(y)

    type(sll_t_bspline_interpolation_1d)    :: self 
    sll_real64, intent(in)  :: x
    sll_real64              :: y
    !local variables
    sll_int32               :: icell
    sll_int32               :: ib
    sll_int32               :: j
    sll_real64              :: val

    ! get bspline values at x
    icell =  sll_f_find_cell( self%bsp, x )
    call sll_s_splines_at_x(self%bsp, icell, x, self%values)
    y = 0.0_f64
    do j=1,self%deg+1
       !       ib = mod(icell+j-self%deg/2-2+self%n,self%n) + 1
       ib = mod(icell+j-2-self%offset+self%n,self%n) + 1
       y = y + self%values(j)*self%bcoef(ib)
    !        print*, icell,j, ib, self%values(j), y
    enddo

  end function sll_f_interpolate_value_1d

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
!> already operated on by the sll_s_compute_bspline_1d() subroutine.
subroutine sll_s_interpolate_array_values_1d( self, n, x, y)

  type(sll_t_bspline_interpolation_1d)     :: self
  sll_int32,            intent(in)  :: n
  sll_real64,           intent(in)  :: x(n)
  sll_real64,           intent(out) :: y(n)
  ! local variables    
  sll_real64 :: val
  sll_real64 :: xx
  sll_int32  :: i
  sll_int32  :: j
  sll_int32  :: ib
  sll_int32  :: icell

  do i= 1, n
     xx = x(i)
     val = 0.0_f64
     ! get bspline values at x
     icell =  sll_f_find_cell( self%bsp, xx )
     call sll_s_splines_at_x(self%bsp, icell, xx, self%values)
     do j=1, self%deg+1
        !ib = mod(icell+j-self%deg/2-2+self%n,self%n) + 1
        ib = mod(icell+j-2-self%offset+self%n,self%n) + 1
        val = val + self%values(j)*self%bcoef(ib)
     enddo
     y(i) = val 
  end do

end subroutine sll_s_interpolate_array_values_1d


!> @brief returns the values of the derivatives evaluated at a point.
!> The spline coefficients used are stored in the spline object.
!> @param[inout] self the spline object pointer, duly initialized and 
!> already operated on by the sll_s_compute_bspline_1d() subroutine.
!> @param[in] x  abscissa to be interpolated.
!> @param[out] y  result of the interpolation.
function sll_f_interpolate_derivative_1d( self, x) result(y)

  type(sll_t_bspline_interpolation_1d)    :: self 
  sll_real64, intent(in)  :: x
  sll_real64              :: y

    !local variables
  sll_int32               :: icell
  sll_int32               :: ib
  sll_int32               :: j
  sll_real64              :: val
  
  ! get bspline derivatives at x
  icell =  sll_f_find_cell( self%bsp, x )
  call sll_s_spline_derivatives_at_x(self%bsp, icell, x, self%values)
  y = 0.0_f64
  do j=1,self%deg+1
     !ib = mod(icell+j-self%deg/2-2+self%n,self%n) + 1
     ib = mod(icell+j-2-self%offset+self%n,self%n) + 1
     y = y + self%values(j)*self%bcoef(ib)
!     print*, icell,j, ib, self%bcoef(ib), self%values(j), y
  enddo
  ! scale
  !y = y / self%length
end function sll_f_interpolate_derivative_1d


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
!> already operated on by the sll_s_compute_bspline_1d() subroutine.
subroutine sll_s_interpolate_array_derivatives_1d( self, n, x, y)
  type(sll_t_bspline_interpolation_1d)    :: self 
  sll_int32,  intent(in)  :: n
  sll_real64, intent(in)  :: x(n)
  sll_real64, intent(out) :: y(n)

  ! local variables    
  sll_real64 :: val
  sll_real64 :: xx
  sll_int32  :: i
  sll_int32  :: j
  sll_int32  :: ib
  sll_int32  :: icell

  do i= 1, n
     xx = x(i)
     val = 0.0_f64
     ! get bspline derivatives at xx
     icell =  sll_f_find_cell( self%bsp, xx )
     call sll_s_spline_derivatives_at_x(self%bsp, icell, xx, self%values)
     do j=1, self%deg+1
        !ib = mod(icell+j-self%deg/2-2+self%n,self%n) + 1
        ib = mod(icell+j-2-self%offset+self%n,self%n) + 1
        val = val + self%values(j)*self%bcoef(ib)
     enddo
     y(i) = val 
  end do

end subroutine sll_s_interpolate_array_derivatives_1d

!> @brief Destructor. Frees all pointers of object
!> @param[inout] self object to be freed
subroutine sll_s_bspline_interpolation_1d_free( self )
  type(sll_t_bspline_interpolation_1d) :: self
  sll_int32 :: ierr

  ! deallocate all pointers
  SLL_DEALLOCATE_ARRAY(self%bcoef,ierr)
  SLL_DEALLOCATE_ARRAY(self%tau,ierr)
  SLL_DEALLOCATE_ARRAY(self%q,ierr)
  SLL_DEALLOCATE_ARRAY(self%values,ierr)
  SLL_DEALLOCATE_ARRAY(self%bsdx,ierr)

  ! free attribute objects
  call sll_s_arbitrary_degree_spline_1d_free(self%bsp)
  call schur_complement_free(self%schur)
  
end subroutine sll_s_bspline_interpolation_1d_free

  
!>@brief
!> produces the B-spline coefficients of an interpolating spline.
!> @details
!> If interpolants are already computed and if you change only 
!> the data set and knots and points positions did not change
!> use self routine.
!>
!!$subroutine sll_s_update_bspline_1d(self, gtau) !, slope_min, slope_max)
!!$
!!$  type(sll_t_bspline_interpolation_1d)    :: self 
!!$  sll_real64, intent(in)  :: gtau(:)
!!$!  sll_real64, optional    :: slope_min
!!$!  sll_real64, optional    :: slope_max
!!$
!!$
!!$  sll_int32               :: n
!!$  sll_int32               :: k
!!$  sll_int32, parameter    :: m = 2
!!$  sll_real64              :: slope(0:1)
!!$
!!$  n = self%n
!!$  k = self%k
!!$
!!$  SLL_ASSERT(size(gtau) == n)
!!$
!!$  if (self%bc_type == sll_p_periodic) then
!!$
!!$    self%bcoef = gtau
!!$    call banslv ( self%q, k+k-1, n, k-1, k-1, self%bcoef )
!!$
!!$  else
!!$
!!$    self%bcoef(1)   = gtau(1)
!!$
!!$    if(self%compute_sl) then
!!$      call sll_s_apply_fd(k+1,1,self%tau(1:k+1),gtau(1:k+1), &
!!$        self%tau(1),slope(0:1))
!!$      self%bcoef(2) = slope(1)
!!$    else if (present(slope_min)) then
!!$      self%bcoef(2) = slope_min
!!$    else
!!$      self%bcoef(2) = self%sl
!!$    end if
!!$    self%bcoef(3:n) = gtau(2:n-1)
!!$    if(self%compute_sr) then
!!$      call sll_s_apply_fd(k+1,1,self%tau(n-k-1:n),gtau(n-k-1:n), &
!!$        self%tau(n),slope(0:1))
!!$      self%bcoef(n+1) = slope(1)
!!$    else if (present(slope_max)) then
!!$      self%bcoef(n+1) = slope_max
!!$    else
!!$      self%bcoef(n+1) = self%sr
!!$    end if
!!$    self%bcoef(n+2) = gtau(n)
!!$     
!!$     self%bcoef = gtau
!!$!ES ?     call banslv ( self%q, k+k-1, n+m, k-1, k-1, self%bcoef )
!!$
!!$  end if
!!$  
!!$end subroutine sll_s_update_bspline_1d

  
  ! *************************************************************************
  !
  !                   2D SPLINE INTERPOLATION
  !
  ! *************************************************************************

!> @brief Initialises a 2D spline interpolation object.
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
!> @return a spline interpolation object.
subroutine sll_s_bspline_interpolation_2d_init( self, nx1, degree1, x1_min, x1_max, bc1, &
                         nx2, degree2, x2_min, x2_max, bc2 )

  type(sll_t_bspline_interpolation_2d)    :: self

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

  sll_int32                        :: n1
  sll_int32                        :: n2
  sll_int32                        :: ierr

  call sll_s_bspline_interpolation_1d_init(self%bs1,nx1,degree1,x1_min,x1_max,bc1)
  call sll_s_bspline_interpolation_1d_init(self%bs2,nx2,degree2,x2_min,x2_max,bc2)

  n1 = size(self%bs1%bcoef)
  n2 = size(self%bs2%bcoef)
  SLL_CLEAR_ALLOCATE(self%bcoef(1:n1,1:n2), ierr)

end subroutine sll_s_bspline_interpolation_2d_init
  
!!$subroutine compute_bspline_2d_with_constant_slopes(self, gtau, &
!!$  sl1_l, sl1_r, sl2_l, sl2_r)
!!$
!!$  type(sll_t_bspline_interpolation_2d)    :: self 
!!$  sll_real64, intent(in)  :: gtau(:,:)
!!$  sll_real64, optional    :: sl1_l
!!$  sll_real64, optional    :: sl1_r
!!$  sll_real64, optional    :: sl2_l
!!$  sll_real64, optional    :: sl2_r
!!$
!!$  if ( self%bs1%bc_type == sll_p_periodic) then
!!$    call build_system_periodic(self%bs1)
!!$  else
!!$    call build_system_with_derivative(self%bs1)
!!$  end if
!!$
!!$  if ( self%bs2%bc_type == sll_p_periodic) then
!!$    call build_system_periodic(self%bs2)
!!$  else
!!$    call build_system_with_derivative(self%bs2)
!!$  end if
!!$
!!$  if (present(sl1_l)) self%x1_min_slopes = sl1_l
!!$  if (present(sl1_r)) self%x1_max_slopes = sl1_r
!!$  if (present(sl2_l)) self%x2_min_slopes = sl2_l
!!$  if (present(sl2_r)) self%x2_max_slopes = sl2_r
!!$
!!$  call update_bspline_2d(self, gtau)
!!$
!!$end subroutine compute_bspline_2d_with_constant_slopes
!!$
!!$subroutine compute_bspline_2d_with_variable_slopes(self, gtau, &
!!$  sl1_l, sl1_r, sl2_l, sl2_r)
!!$
!!$  type(sll_t_bspline_interpolation_2d)    :: self 
!!$  sll_real64, intent(in)  :: gtau(:,:)
!!$  sll_real64, intent(in)  :: sl1_l(:)
!!$  sll_real64, intent(in)  :: sl1_r(:)
!!$  sll_real64, intent(in)  :: sl2_l(:)
!!$  sll_real64, intent(in)  :: sl2_r(:)
!!$
!!$  if ( self%bs1%bc_type == sll_p_periodic) then
!!$    call build_system_periodic(self%bs1)
!!$  else
!!$    call build_system_with_derivative(self%bs1)
!!$  end if
!!$
!!$  if ( self%bs2%bc_type == sll_p_periodic) then
!!$    call build_system_periodic(self%bs2)
!!$  else
!!$    call build_system_with_derivative(self%bs2)
!!$  end if
!!$
!!$  SLL_ASSERT(size(sl1_l) == self%bs2%n)
!!$  self%x1_min_slopes = sl1_l
!!$  SLL_ASSERT(size(sl1_r) == self%bs2%n)
!!$  self%x1_max_slopes = sl1_r
!!$  SLL_ASSERT(size(sl2_l) == self%bs1%n)
!!$  self%x2_min_slopes = sl2_l
!!$  SLL_ASSERT(size(sl2_r) == self%bs1%n)
!!$  self%x2_max_slopes = sl2_r
!!$
!!$  call update_bspline_2d(self, gtau)
!!$
!!$end subroutine compute_bspline_2d_with_variable_slopes


!!$!> @brief 
!!$!> update 2 values before computing bsplines coefficients
!!$!> @details
!!$!> If the points positions did not change use self function instead
!!$!> of sll_o_compute_bspline_2d. You still need to call sll_o_compute_bspline_2d at
!!$!> the beginning to build the linear system.
!!$subroutine update_bspline_2d(self, gtau) !, sl1_l, sl1_r, sl2_l, sl2_r)
!!$
!!$  type(sll_t_bspline_interpolation_2d)    :: self 
!!$  sll_real64, intent(in)  :: gtau(:,:)
!  sll_real64, optional    :: sl1_l
!  sll_real64, optional    :: sl1_r
!  sll_real64, optional    :: sl2_l
!  sll_real64, optional    :: sl2_r
!!$  sll_real64, allocatable :: bwork(:,:)
!!$  sll_real64, allocatable :: coeff1(:)
!!$  sll_real64, allocatable :: coeff2(:)
!!$
!!$  sll_int32               :: i
!!$  sll_int32               :: j
!!$  sll_int32               :: n1
!!$  sll_int32               :: n2
!!$  sll_int32               :: bc1
!!$  sll_int32               :: bc2
!!$  sll_int32               :: ierr
!!$
!!$  sll_int32               :: m1
!!$  
!!$  n1 = self%bs1%n
!!$  n2 = self%bs2%n
!!$
!!$  bc1 = self%bs1%bc_type
!!$  bc2 = self%bs2%bc_type
!!$
!!$  if ( bc1 == sll_p_periodic ) then
!!$    m1 = 0
!!$  else
!!$    m1 = 2
!!$  end if
!!$
!!$  SLL_CLEAR_ALLOCATE(bwork(1:n2,1:n1+m1),ierr)
!!$
!  if (present(sl1_l)) self%x1_min_slopes = sl1_l
!  if (present(sl1_r)) self%x1_max_slopes = sl1_r
!  if (present(sl2_l)) self%x2_min_slopes = sl2_l
!  if (present(sl2_r)) self%x2_max_slopes = sl2_r
!!$
!!$  do j = 1, n2
!!$    call sll_s_update_bspline_1d( self%bs1, gtau(:,j)) !, &
!!$!ES      self%x1_min_slopes(j), self%x1_max_slopes(j))
!!$    bwork(j,:) = self%bs1%bcoef
!!$  end do
!!$  
!  SLL_CLEAR_ALLOCATE(coeff1(1:n1+m1), ierr)
!  SLL_CLEAR_ALLOCATE(coeff2(1:n1+m1), ierr)

 ! coeff1(1:n1) = self%x2_min_slopes(:)
!  coeff2(1:n1) = self%x2_max_slopes(:)

!  do i = 1, n1+m1
!    call sll_s_update_bspline_1d( self%bs2, bwork(:,i), coeff1(i), coeff2(i))
!    self%bcoef(i,:) = self%bs2%bcoef(:)
!  end do
!!$
!!$  deallocate(bwork)
!!$
!!$end subroutine update_bspline_2d






subroutine delete_bspline_2D( spline )
  type(sll_t_bspline_interpolation_2d) :: spline
  !print*, associated(spline)
end subroutine delete_bspline_2D 

!!$subroutine sll_s_interpolate_array_values_2d(self, n1, n2, x, y, ideriv, jderiv)
!!$
!!$type(sll_t_bspline_interpolation_2d)    :: self
!!$sll_int32               :: n1
!!$sll_int32               :: n2
!!$sll_real64, intent(in)  :: x(:,:)
!!$sll_real64, intent(out) :: y(:,:)
!!$sll_int32               :: ideriv
!!$sll_int32               :: jderiv
!!$
!!$sll_int32               :: i
!!$sll_int32               :: j, jj
!!$sll_int32               :: jc, jcmin, jcmax
!!$
!!$sll_real64, allocatable :: ajx(:), ajy(:)
!!$sll_real64, allocatable :: dlx(:), dly(:)
!!$sll_real64, allocatable :: drx(:), dry(:)
!!$sll_real64, allocatable :: wrk(:)
!!$
!!$sll_int32               :: nx, kx, ny, ky
!!$sll_int32               :: left, leftx, lefty
!!$sll_int32               :: jlo
!!$sll_int32               :: klo
!!$sll_int32               :: llo
!!$sll_int32               :: mflag
!!$sll_int32               :: ierr
!!$sll_int32               :: jjj
!!$sll_int32               :: kkk
!!$sll_int32               :: nmkx
!!$sll_int32               :: nmky
!!$
!!$sll_real64              :: xi
!!$sll_real64              :: xj
!!$sll_real64, pointer     :: tx(:)
!!$sll_real64, pointer     :: ty(:)
!!$type(sll_t_deboor_type)       :: db
!!$
!!$nx   =  self%bs1%n
!!$ny   =  self%bs2%n
!!$SLL_ASSERT(n1 <= nx)
!!$SLL_ASSERT(n2 <= ny)
!!$SLL_ASSERT(n1 == size(x,1))
!!$SLL_ASSERT(n2 == size(x,2))
!!$kx   =  self%bs1%deg + 1
!!$ky   =  self%bs2%deg + 1
!!$tx   => self%bs1%t
!!$ty   => self%bs2%t
!!$
!!$if (self%bs1%bc_type == sll_p_periodic) then
!!$  nmkx = nx+kx
!!$else
!!$  nmkx = nx+kx+2
!!$end if
!!$if (self%bs2%bc_type == sll_p_periodic) then
!!$  nmky = ny+ky
!!$else
!!$  nmky = ny+ky+2
!!$end if
!!$
!!$SLL_CLEAR_ALLOCATE(ajx(1:kx),ierr)
!!$SLL_CLEAR_ALLOCATE(dlx(1:kx),ierr)
!!$SLL_CLEAR_ALLOCATE(drx(1:kx),ierr)
!!$SLL_CLEAR_ALLOCATE(ajy(1:ky),ierr)
!!$SLL_CLEAR_ALLOCATE(dly(1:ky),ierr)
!!$SLL_CLEAR_ALLOCATE(dry(1:ky),ierr)
!!$SLL_CLEAR_ALLOCATE(wrk(1:nmkx),ierr)
!!$
!!$jlo = ky
!!$do j=1,n2
!!$  db%ilo = kx
!!$  klo = jlo
!!$  xj  = self%bs2%tau(j)
!!$  call sll_s_interv(db,ty,nmky,xj,lefty,mflag)
!!$  do i=1,nx
!!$    xi = self%bs1%tau(i)
!!$    call sll_s_interv(db,tx,nmkx,xi,leftx,mflag)
!!$    do jj=1,ky
!!$      jcmin = 1
!!$      if ( kx <= leftx ) then
!!$        do jjj = 1, kx-1
!!$          dlx(jjj) = xi - tx(leftx+1-jjj)
!!$        end do
!!$      else
!!$        jcmin = 1-(leftx-kx)
!!$        do jjj = 1, leftx
!!$          dlx(jjj) = xi - tx(leftx+1-jjj)
!!$        end do
!!$        do jjj = leftx, kx-1
!!$          ajx(kx-jjj) = 0.0_f64
!!$          dlx(jjj) = dlx(leftx)
!!$        end do
!!$      end if
!!$      jcmax = kx
!!$      if ( nmkx-kx < leftx ) then
!!$        jcmax = nmkx-leftx
!!$        do jjj = 1, nmkx-leftx
!!$          drx(jjj) = tx(leftx+jjj) - xi
!!$        end do
!!$        do jjj = nmkx-leftx, kx-1
!!$          ajx(jjj+1) = 0.0_f64
!!$          drx(jjj) = drx(nmkx-leftx)
!!$        end do
!!$      else
!!$        do jjj = 1, kx-1
!!$          drx(jjj) = tx(leftx+jjj) - xi
!!$        end do
!!$      end if
!!$      do jc = jcmin, jcmax
!!$        ajx(jc) = self%bcoef(leftx-kx+jc,lefty-ky+jj)
!!$      end do
!!$      do jjj = 1, ideriv
!!$        llo = kx - jjj
!!$        do kkk = 1, kx - jjj
!!$          ajx(kkk) = ((ajx(kkk+1)-ajx(kkk))/(dlx(llo)+drx(kkk)))*(kx-jjj)
!!$          llo = llo-1
!!$        end do
!!$      end do
!!$      do jjj = ideriv+1, kx-1
!!$        llo = kx-jjj
!!$        do kkk = 1, kx-jjj
!!$          ajx(kkk) = (ajx(kkk+1)*dlx(llo)+ajx(kkk)*drx(kkk)) &
!!$                     /(dlx(llo)+drx(kkk))
!!$          llo = llo - 1
!!$        end do
!!$      end do
!!$      wrk(jj) = ajx(1)
!!$    end do
!!$    call sll_s_interv(db,ty(lefty-ky+1:nmky),ky+ky,xj,left,mflag)
!!$    jcmin = 1
!!$    if ( ky <= left ) then
!!$      do jjj = 1, ky-1
!!$        dly(jjj) = xj - ty(lefty-ky+left+1-jjj)
!!$      end do
!!$    else
!!$      jcmin = 1-(left-ky)
!!$      do jjj = 1, left
!!$        dly(jjj) = xj - ty(lefty-ky+left+1-jjj)
!!$      end do
!!$      do jjj = left, ky-1
!!$        ajy(ky-jjj) = 0.0_f64
!!$        dly(jjj) = dly(left)
!!$      end do
!!$    end if
!!$    jcmax = ky
!!$    if ( ky < left ) then
!!$      jcmax = ky+ky-left
!!$      do jjj = 1, ky+ky-left
!!$        dry(jjj) = ty(lefty-ky+left+jjj) - xj
!!$      end do
!!$      do jjj = ky+ky-left, ky-1
!!$        ajy(jjj+1) = 0.0_f64
!!$        dry(jjj) = dry(ky+ky-left)
!!$      end do
!!$    else
!!$      do jjj = 1, ky-1
!!$        dry(jjj) = ty(lefty-ky+left+jjj) - xj
!!$      end do
!!$    end if
!!$    do jc = jcmin, jcmax
!!$      ajy(jc) = wrk(left-ky+jc)
!!$    end do
!!$    do jjj = 1, jderiv
!!$      llo = ky - jjj
!!$      do kkk = 1, ky - jjj
!!$        ajy(kkk) = ((ajy(kkk+1)-ajy(kkk))/(dly(llo)+dry(kkk)))*(ky-jjj)
!!$        llo = llo-1
!!$      end do
!!$    end do
!!$    do jjj = jderiv+1, ky-1
!!$      llo = ky-jjj
!!$      do kkk = 1, ky-jjj
!!$        ajy(kkk) = (ajy(kkk+1)*dly(llo)+ajy(kkk)*dry(kkk))/(dly(llo)+dry(kkk))
!!$        llo = llo - 1
!!$      end do
!!$    end do
!!$    y(i,j) = ajy(1)
!!$  end do
!!$end do
!!$
!!$deallocate(ajx)
!!$deallocate(dlx)
!!$deallocate(drx)
!!$deallocate(ajy)
!!$deallocate(dly)
!!$deallocate(dry)
!!$deallocate(wrk)
!!$
!!$end subroutine sll_s_interpolate_array_values_2d
!!$
!!$function sll_f_interpolate_value_2d(self, xi, xj, ideriv, jderiv ) result (y)
!!$
!!$type(sll_t_bspline_interpolation_2d)    :: self
!!$sll_real64, intent(in)  :: xi
!!$sll_real64, intent(in)  :: xj
!!$sll_int32,  intent(in)  :: ideriv
!!$sll_int32,  intent(in)  :: jderiv
!!$sll_real64              :: y
!!$
!!$sll_int32               :: jj
!!$sll_int32               :: jc, jcmin, jcmax
!!$sll_int32               :: nx, kx, ny, ky
!!$sll_int32               :: left, leftx, lefty
!!$sll_int32               :: llo
!!$sll_int32               :: mflag
!!$sll_int32               :: jjj
!!$sll_int32               :: kkk
!!$sll_int32               :: nmkx
!!$sll_int32               :: nmky
!!$
!!$sll_real64, pointer     :: tx(:)
!!$sll_real64, pointer     :: ty(:)
!!$
!!$sll_real64, allocatable :: work(:)
!!$type(sll_t_deboor_type)       :: db
!!$
!!$nx   =  self%bs1%n
!!$ny   =  self%bs2%n
!!$kx   =  self%bs1%deg + 1
!!$ky   =  self%bs2%deg + 1
!!$tx   => self%bs1%t
!!$ty   => self%bs2%t
!!$
!!$allocate(work(size(self%bs1%bcoef)))
!!$work = 0.0_f64
!!$
!!$if (self%bs1%bc_type == sll_p_periodic) then
!!$  nmkx = nx+kx
!!$else
!!$  nmkx = nx+kx+2
!!$end if
!!$if (self%bs1%bc_type == sll_p_periodic) then
!!$  nmky = ny+ky
!!$else
!!$  nmky = ny+ky+2
!!$end if
!!$
!!$call sll_s_interv(db,tx,nmkx,xi,leftx,mflag)
!!$call sll_s_interv(db,ty,nmky,xj,lefty,mflag)
!!$
!!$do jj=1,ky
!!$  jcmin = 1
!!$  if ( kx <= leftx ) then
!!$    do jjj = 1, kx-1
!!$      self%bs1%dl(jjj) = xi - tx(leftx+1-jjj)
!!$    end do
!!$  else
!!$    jcmin = 1-(leftx-kx)
!!$    do jjj = 1, leftx
!!$      self%bs1%dl(jjj) = xi - tx(leftx+1-jjj)
!!$    end do
!!$    do jjj = leftx, kx-1
!!$      self%bs1%aj(kx-jjj) = 0.0_f64
!!$      self%bs1%dl(jjj) = self%bs1%dl(leftx)
!!$    end do
!!$  end if
!!$  jcmax = kx
!!$  if ( nmkx-kx < leftx ) then
!!$    jcmax = nmkx-leftx
!!$    do jjj = 1, nmkx-leftx
!!$      self%bs1%dr(jjj) = tx(leftx+jjj) - xi
!!$    end do
!!$    do jjj = nmkx-leftx, kx-1
!!$      self%bs1%aj(jjj+1) = 0.0_f64
!!$      self%bs1%dr(jjj) = self%bs1%dr(nmkx-leftx)
!!$    end do
!!$  else
!!$    do jjj = 1, kx-1
!!$      self%bs1%dr(jjj) = tx(leftx+jjj) - xi
!!$    end do
!!$  end if
!!$  do jc = jcmin, jcmax
!!$    self%bs1%aj(jc) = self%bcoef(leftx-kx+jc,lefty-ky+jj)
!!$  end do
!!$  do jjj = 1, ideriv
!!$    llo = kx - jjj
!!$    do kkk = 1, kx - jjj
!!$      self%bs1%aj(kkk) = ((self%bs1%aj(kkk+1)-self%bs1%aj(kkk)) &
!!$        /(self%bs1%dl(llo)+self%bs1%dr(kkk)))*(kx-jjj)
!!$      llo = llo-1
!!$    end do
!!$  end do
!!$  do jjj = ideriv+1, kx-1
!!$    llo = kx-jjj
!!$    do kkk = 1, kx-jjj
!!$      self%bs1%aj(kkk) = (self%bs1%aj(kkk+1)*self%bs1%dl(llo)+ &
!!$                          self%bs1%aj(kkk  )*self%bs1%dr(kkk))/  &
!!$                         (self%bs1%dl(llo  )+self%bs1%dr(kkk))
!!$      llo = llo - 1
!!$    end do
!!$  end do
!!$  work(jj) = self%bs1%aj(1)
!!$end do
!!$
!!$!klo = self%bs2%ilo
!!$call sll_s_interv(db,ty(lefty-ky+1:nmky),ky+ky,xj,left,mflag)
!!$
!!$jcmin = 1
!!$if ( ky <= left ) then
!!$  do jjj = 1, ky-1
!!$    self%bs2%dl(jjj) = xj - ty(lefty-ky+left+1-jjj)
!!$  end do
!!$else
!!$  jcmin = 1-(left-ky)
!!$  do jjj = 1, left
!!$    self%bs2%dl(jjj) = xj - ty(lefty-ky+left+1-jjj)
!!$  end do
!!$  do jjj = left, ky-1
!!$    self%bs2%aj(ky-jjj) = 0.0_f64
!!$    self%bs2%dl(jjj) = self%bs2%dl(left)
!!$  end do
!!$end if
!!$jcmax = ky
!!$if ( ky < left ) then
!!$  jcmax = ky+ky-left
!!$  do jjj = 1, ky+ky-left
!!$    self%bs2%dr(jjj) = ty(lefty-ky+left+jjj) - xj
!!$  end do
!!$  do jjj = ky+ky-left, ky-1
!!$    self%bs2%aj(jjj+1) = 0.0_f64
!!$    self%bs2%dr(jjj) = self%bs2%dr(ky+ky-left)
!!$  end do
!!$else
!!$  do jjj = 1, ky-1
!!$    self%bs2%dr(jjj) = ty(lefty-ky+left+jjj) - xj
!!$  end do
!!$end if
!!$do jc = jcmin, jcmax
!!$  self%bs2%aj(jc) = work(left-ky+jc)
!!$end do
!!$do jjj = 1, jderiv
!!$  llo = ky - jjj
!!$  do kkk = 1, ky - jjj
!!$    self%bs2%aj(kkk) = ((self%bs2%aj(kkk+1)-self%bs2%aj(kkk)) &
!!$      /(self%bs2%dl(llo)+self%bs2%dr(kkk)))*(ky-jjj)
!!$    llo = llo-1
!!$  end do
!!$end do
!!$do jjj = jderiv+1, ky-1
!!$  llo = ky-jjj
!!$  do kkk = 1, ky-jjj
!!$    self%bs2%aj(kkk) = (self%bs2%aj(kkk+1)*self%bs2%dl(llo)+ &
!!$                        self%bs2%aj(kkk)*self%bs2%dr(kkk))/  &
!!$                       (self%bs2%dl(llo)+self%bs2%dr(kkk))
!!$    llo = llo - 1
!!$  end do
!!$end do
!!$y = self%bs2%aj(1)
!!$
!!$end function sll_f_interpolate_value_2d
!!$
!!$subroutine interpolate_array_x1_derivatives_2d(self, n1, n2, x, y)
!!$
!!$type(sll_t_bspline_interpolation_2d)    :: self
!!$sll_int32,  intent(in)  :: n1
!!$sll_int32,  intent(in)  :: n2
!!$sll_real64, intent(in)  :: x(:,:)
!!$sll_real64, intent(out) :: y(:,:)
!!$
!!$SLL_ASSERT(self%bs1%n>0)
!!$SLL_ASSERT(n1 == size(x,1))
!!$SLL_ASSERT(n2 == size(y,2))
!!$y = 0.0_f64
!!$call sll_s_interpolate_array_values_2d(self, n1, n2, x, y, 1, 0)
!!$
!!$end subroutine interpolate_array_x1_derivatives_2d
!!$
!!$subroutine interpolate_array_x2_derivatives_2d(self, n1, n2, x, y)
!!$
!!$type(sll_t_bspline_interpolation_2d)    :: self
!!$sll_int32,  intent(in)  :: n1
!!$sll_int32,  intent(in)  :: n2
!!$sll_real64, intent(in)  :: x(:,:)
!!$sll_real64, intent(out) :: y(:,:)
!!$
!!$SLL_ASSERT(self%bs1%n>0)
!!$SLL_ASSERT(n1 == size(x,1))
!!$SLL_ASSERT(n2 == size(y,2))
!!$y = 0.0_f64
!!$call sll_s_interpolate_array_values_2d(self, n1, n2, x, y, 0, 1)
!!$
!!$end subroutine interpolate_array_x2_derivatives_2d

end module sll_m_bspline_interpolation
