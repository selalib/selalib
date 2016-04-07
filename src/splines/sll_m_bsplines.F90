!> @ingroup splines
!> Implements arbitrary degree bspline interpolation on a uniform grid
!> given a B-Spline object from sll_m_arbitrary_degree_splines
module sll_m_bspline_interpolation

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_boundary_condition_descriptors, only: &
       sll_p_periodic, &
       sll_p_hermite, &
       sll_p_greville, &
       sll_p_mirror

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

  implicit none

  public ::                                       &
       sll_t_bspline_interpolation_1d,            &
       sll_t_bspline_interpolation_2d,            &
       sll_s_bspline_interpolation_1d_init,       &
       sll_s_bspline_interpolation_2d_init,       &
       sll_s_bspline_interpolation_1d_free,       &
       sll_s_bspline_interpolation_2d_free,       &
       sll_s_compute_bspline_1d,                  &
       sll_s_compute_bspline_2d,                  &
       sll_f_interpolate_value_1d,                &
       sll_f_interpolate_value_2d,                &
       sll_s_interpolate_array_values_1d,         &
       sll_s_interpolate_array_values_2d,         &
       sll_f_interpolate_derivative_1d,           &
       sll_f_interpolate_derivative_x1_2d,        &
       sll_f_interpolate_derivative_x2_2d,        &
       sll_s_interpolate_array_derivatives_1d,    &
       sll_s_interpolate_array_derivatives_x1_2d, &
       sll_s_interpolate_array_derivatives_x2_2d

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
!> @brief 
!> basic type for one-dimensional B-spline interpolation. 
type :: sll_t_bspline_interpolation_1d

  type (sll_t_arbitrary_degree_spline_1d) :: bsp
  sll_int32                               :: n        ! dimension of spline space
  sll_int32                               :: deg      ! degree of spline
  sll_real64, pointer                     :: tau(:)   ! n interpolation points
  sll_real64, pointer                     :: bcoef(:) ! bspline coefficients
  sll_int32                               :: bc_type  ! boundary condition
  sll_real64, pointer                     :: q(:,:)   ! triangular factorization
  sll_real64, pointer                     :: values(:) 
  sll_real64, pointer                     :: bsdx(:,:)
  sll_int32,  pointer                     :: ipiv(:)
  sll_real64                              :: length
  sll_int32                               :: offset
  type(schur_complement_solver)           :: schur
  sll_real64, pointer                     :: bc_left(:)
  sll_real64, pointer                     :: bc_right(:)

end type sll_t_bspline_interpolation_1d

!> @brief 
!> basic type for two-dimensional B-spline data. 
!> @details 
!> treated as an opaque type. No access to its internals is directly allowed.
type :: sll_t_bspline_interpolation_2d

  type(sll_t_bspline_interpolation_1d) :: bs1
  type(sll_t_bspline_interpolation_1d) :: bs2
  sll_real64, pointer                  :: bcoef(:,:)
  sll_real64, pointer                  :: bwork(:,:)

end type sll_t_bspline_interpolation_2d

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
  subroutine sll_s_bspline_interpolation_1d_init( self, &
    num_points,     &
    degree,         &
    xmin,           &
    xmax,           &
    bc_type,        &
    spline_bc_type, &
    bc_left,        &
    bc_right)

    type(sll_t_bspline_interpolation_1d) :: self
    sll_int32,  intent(in)               :: num_points
    sll_int32,  intent(in)               :: degree
    sll_real64, intent(in)               :: xmin
    sll_real64, intent(in)               :: xmax
    sll_int32,  intent(in)               :: bc_type
    sll_int32,  optional                 :: spline_bc_type
    sll_real64, optional                 :: bc_left(:)
    sll_real64, optional                 :: bc_right(:)
    
    ! local variables
    sll_real64, dimension(num_points)    :: grid
    sll_int32                            :: ierr
    sll_int32                            :: i
    sll_int32                            :: j
    sll_real64                           :: delta

    ! knot point mirroring works only for Hermite boundary conditions. Check this
    if (present(spline_bc_type).and.(spline_bc_type == sll_p_mirror)) then
       if (.not.(bc_type == sll_p_hermite)) then
          print*, 'Mirror knot sequence at boundary only works with Hermite BC'
          stop
       end if
    end if
    
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
    if (present(spline_bc_type)) then
       call sll_s_arbitrary_degree_spline_1d_init(self%bsp, degree, grid, num_points, &
            bc_type, spline_bc_type )
    else
       call sll_s_arbitrary_degree_spline_1d_init(self%bsp, degree, grid, num_points, &
            bc_type)
    end if

    SLL_ALLOCATE(self%bcoef(self%n),  ierr)
    SLL_ALLOCATE(self%values(self%deg+1), ierr)
    SLL_ALLOCATE(self%bsdx(self%deg/2+1,self%deg+1), ierr)
    SLL_ALLOCATE(self%ipiv(self%n), ierr)
    
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
       ! deal with possible boundary conditions
       if (present(bc_left).and. present(bc_right)) then
          SLL_ASSERT((size(bc_left)==self%deg/2))
          SLL_ASSERT((size(bc_right)==self%deg/2))
          SLL_ALLOCATE(self%bc_left(self%deg/2),ierr)
          self%bc_left = bc_left
          SLL_ALLOCATE(self%bc_right(self%deg/2),ierr)
          self%bc_right = bc_right
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

    ! Assemble banded matrix (B_j(tau(i))) for spline interpolation
    select case (bc_type)
    case (sll_p_periodic) 
       call build_system_periodic(self) 
    case (sll_p_greville)
       call build_system_greville(self)
    case (sll_p_hermite)
       call build_system_with_derivative(self)
    end select
    
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
    type(sll_t_bspline_interpolation_1d) :: self

    ! local variables
    sll_int32                            :: iflag
    sll_int32                            :: j
    sll_int32                            :: k
    sll_int32                            :: ii
    sll_int32                            :: jj
    sll_int32                            :: icell
    sll_real64                           :: x
    
    ! The matrix q is a banded matrix using the storage required by banfac (De Boor)
    ! It has k bands above diagonal, k bands below and the diagonal itself
    ! The term A(ii,jj) of the full matrix is stored in q(ii-jj+k+1,jj)
    ! The Bspline interpolation matrix at Greville points x_ii is
    ! A(ii,jj) = B_jj(x_ii)
    k = self%deg - 1 
    SLL_ALLOCATE(self%q(2*k+1,self%n), iflag)
    self%q = 0.0_f64
    
    ! Treat i=1 separately
    x =  self%bsp%xmin
    icell = 1
    ii = 1   ! first Greville point
    call sll_s_splines_at_x(self%bsp, icell, x, self%values)
    ! iterate only to k+1 as last bspline is 0
    do j=1,k+1
       jj = icell+j-1
       self%q(ii-jj+k+1,jj) = self%values(j)
    end do

    do ii=2,self%n - 1
       x = self%tau(ii)
       icell = sll_f_find_cell(self%bsp, x )
       call sll_s_splines_at_x(self%bsp, icell, x, self%values)
       do j=1,k+2
          jj = icell+j-1
          self%q(ii-jj+k+1,jj) = self%values(j)
       end do
    end do
    ! Treat i=self%n separately
    x =  self%bsp%xmax
    icell = self%n - self%deg
    ii = self%n  ! last Greville point
    call sll_s_splines_at_x(self%bsp, icell, x, self%values)
    ! iterate only to k as first bspline is 0
    do j=2,k+2
       jj = icell+j-1
       self%q(ii-jj+k+1,jj) = self%values(j)
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
    sll_int32                              :: nbc
    sll_int32                              :: iflag
    sll_int32                              :: i
    sll_int32                              :: j
    sll_int32                              :: ii
    sll_int32                              :: jj
    sll_int32                              :: offset
    sll_int32                              :: k
    sll_int32                              :: icell
    sll_real64                             :: x
    
    ! number of boundary conditions neededdepending on spline degree
    k = self%deg-1
    nbc = self%deg/2

    ! The matrix q is a banded matrix using the storage required by DGBTRF (LAPACK)
    ! It has 2*k bands above diagonal, k bands below and the diagonal itself
    ! k additional bands above diagonal needed for pivoting
    ! The term A(ii,jj) of the full matrix is stored in q(ii-jj+2*k+1,jj)
    ! The Bspline interpolation matrix at Greville points x_ii is
    ! A(ii,jj) = B_jj(x_ii)
    SLL_ALLOCATE(self%q(3*k+1,self%n), iflag)
    self%q = 0.0_f64
     
    ! For even degree splines interpolation points are at cell midpoints
    ! value of function on the boundary needed as additional boundary conditions
    ! for odd degree splines  baundary is in the interpolation points
    ! only derivative values are needed as boundary conditions
    if (modulo(k+1,2)==0) then ! spline degree even
       offset = 0
    else
       offset = 1
    end if
   
    ! boundary conditions at xmin
    ii=0
    x = self%bsp%xmin
    icell = 1
    call sll_s_splines_and_n_derivs_at_x( self%bsp, icell, x , nbc , self%bsdx)
    do i=1, nbc
       ! iterate only to k+1 as last bspline is 0
       ii=ii+1
       do j=1,k+1
          jj = icell+j-1
          self%q(ii-jj+2*k+1,jj) = self%bsdx(i+offset,j)
       end do
    end do
    ! interpolation points
    do i=1,self%n - 2*nbc
       ii = ii + 1
       x = self%tau(i)
       icell = sll_f_find_cell(self%bsp, x )
       call sll_s_splines_at_x(self%bsp, icell, x, self%values)
       do j=1,k+2
          jj = icell+j-1
          self%q(ii-jj+2*k+1,jj) = self%values(j)
       end do
    end do
    ! boundary conditions at xmax
    x = self%bsp%xmax
    icell = self%n - self%deg
    call sll_s_splines_and_n_derivs_at_x( self%bsp, icell, x , nbc , self%bsdx)
    do i= 1, nbc
       ii = ii + 1
       do j=2,k+2
          jj = icell+j-1
          self%q( ii-jj+2*k+1,jj) = self%bsdx(i+offset,j)
       end do
    end do

    ! Perform LU decomposition of matrix q with Lapack
    call dgbtrf(self%n,self%n,k,k,self%q,3*k+1,self%ipiv,iflag)

  end subroutine build_system_with_derivative
  
  subroutine sll_s_compute_bspline_1d(self, gtau, val_min, val_max)

    type(sll_t_bspline_interpolation_1d)   :: self
    sll_real64, intent(in) :: gtau(:)
    sll_real64, optional   :: val_min(:)
    sll_real64, optional   :: val_max(:)
    ! local variables
    sll_int32   :: ncond
    sll_int32   :: k
    sll_int32   :: iflag
    
    select case (self%bc_type)
    case(sll_p_periodic)
       self%bcoef = gtau
       call schur_complement_slv ( self%schur, self%n, self%deg/2, self%q,  self%bcoef )
    case (sll_p_greville)
       self%bcoef = gtau
       call banslv ( self%q, 2*self%deg-1, self%n, self%deg-1, self%deg-1, self%bcoef )
    case (sll_p_hermite)
       ! number of needed conditions at boundary
       ncond = self%deg/2
       if (present(val_min)) then
          self%bcoef(1:ncond) = val_min(1:ncond)
       else  ! set needed boundary values to 0
          self%bcoef(1:ncond) = 0.0_f64
       end if
       self%bcoef(ncond+1:self%n-ncond) = gtau(1:self%n-2*ncond)
       if (present(val_max)) then
          self%bcoef(self%n-ncond+1:self%n) = val_max(1:ncond)
       else ! set needed boundary values to 0
          self%bcoef(self%n-ncond+1:self%n) = 0.0_f64
       end if
       ! Use Lapack to solve banded system
       k = self%deg-1
       call dgbtrs('N',self%n,k,k,1,self%q,3*k+1,self%ipiv,self%bcoef, self%n, iflag)
    end select

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

    ! get bspline values at x
    icell =  sll_f_find_cell( self%bsp, x )
    call sll_s_splines_at_x(self%bsp, icell, x, self%values)
    y = 0.0_f64
    do j=1,self%deg+1
       ib = mod(icell+j-2-self%offset+self%n,self%n) + 1
       y = y + self%values(j)*self%bcoef(ib)
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

  type(sll_t_bspline_interpolation_1d) :: self
  sll_int32,            intent(in)     :: n
  sll_real64,           intent(in)     :: x(n)
  sll_real64,           intent(out)    :: y(n)

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

  type(sll_t_bspline_interpolation_1d) :: self 
  sll_real64, intent(in)               :: x
  sll_real64                           :: y

  !local variables
  sll_int32                            :: icell
  sll_int32                            :: ib
  sll_int32                            :: j
 
  ! get bspline derivatives at x
  icell =  sll_f_find_cell( self%bsp, x )
  call sll_s_spline_derivatives_at_x(self%bsp, icell, x, self%values)
  y = 0.0_f64
  do j=1,self%deg+1
     ib = mod(icell+j-2-self%offset+self%n,self%n) + 1
     y = y + self%values(j)*self%bcoef(ib)
  enddo

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
  
  ! *************************************************************************
  !
  !                   2D SPLINE INTERPOLATION
  !
  ! *************************************************************************

!> @brief Initialises a 2D spline interpolation object.
!> @param[in] nx1 Number of points where the data to be interpolated are 
!>            represented.
!> @param[in] degree1 Spline degree 
!> @param[in] x1_min Minimum value of the abscissae where the data are meant 
!> to be interpolated.
!> @param[in] x1_max Maximum value of the abscissae where the data are meant 
!> to be interpolated.
!> @param[in] bc1 A boundary condition specifier. Must be one of the
!> symbols defined in the SLL_BOUNDARY_CONDITION_DESCRIPTORS module.
!> @param[in] nx2 Number of points where the data to be interpolated are 
!>            represented.
!> @param[in] degree2 Spline degree 
!> @param[in] x2_min Minimum value of the abscissae where the data are meant 
!> to be interpolated.
!> @param[in] x2_max Maximum value of the abscissae where the data are meant 
!> to be interpolated.
!> @param[in] bc2 A boundary condition specifier. Must be one of the
!> symbols defined in the sll_m_boundary_condition_descriptors module.
!> @param[in] spline_bc_type1 A boundary condition specifier (see sll_s_bspline_interpolation_1d_init). 
!> @param[in] spline_bc_type2 A boundary condition specifier (see sll_s_bspline_interpolation_1d_init). 
!> @param[in] bc_left1 value of function on the west boundary
!> @param[in] bc_left2 value of function on the south boundary
!> @param[in] bc_right1 value of function on the east boundary
!> @param[in] bc_right2 value of the function on the north boundary
!> @return a spline interpolation object.
subroutine sll_s_bspline_interpolation_2d_init( &
  self,            &
  nx1,             &
  nx2,             &
  degree1,         &
  degree2,         &
  x1_min,          &
  x2_min,          &
  x1_max,          &
  x2_max,          &
  bc1,             &
  bc2,             &
  spline_bc_type1, &
  spline_bc_type2, &
  bc_left1,        &
  bc_left2,        &
  bc_right1,       &
  bc_right2        )

  type(sll_t_bspline_interpolation_2d) :: self

  sll_int32,  intent(in)               :: nx1
  sll_int32,  intent(in)               :: degree1
  sll_real64, intent(in)               :: x1_min
  sll_real64, intent(in)               :: x1_max
  sll_int32,  intent(in)               :: bc1
  sll_int32,  intent(in)               :: nx2
  sll_int32,  intent(in)               :: degree2
  sll_real64, intent(in)               :: x2_min
  sll_real64, intent(in)               :: x2_max
  sll_int32,  intent(in)               :: bc2

  sll_int32,  optional                 :: spline_bc_type1
  sll_int32,  optional                 :: spline_bc_type2
  sll_real64, optional, pointer        :: bc_left1(:)
  sll_real64, optional, pointer        :: bc_left2(:)
  sll_real64, optional, pointer        :: bc_right1(:)
  sll_real64, optional, pointer        :: bc_right2(:)

  sll_int32                            :: n1
  sll_int32                            :: n2
  sll_int32                            :: ierr

  if (present(spline_bc_type1)) then
    call sll_s_bspline_interpolation_1d_init(self%bs1,nx1,degree1, &
      x1_min,x1_max,bc1, spline_bc_type1)
  else
    call sll_s_bspline_interpolation_1d_init(self%bs1,nx1,degree1, &
      x1_min,x1_max,bc1)
  end if
  if (present(spline_bc_type2)) then
    call sll_s_bspline_interpolation_1d_init(self%bs2,nx2,degree2, &
      x2_min,x2_max,bc2, spline_bc_type2)
  else
    call sll_s_bspline_interpolation_1d_init(self%bs2,nx2,degree2, &
      x2_min,x2_max,bc2)
  end if

  n1 = self%bs1%n
  n2 = self%bs2%n
  SLL_CLEAR_ALLOCATE(self%bwork(1:n2,1:n1), ierr)
  SLL_CLEAR_ALLOCATE(self%bcoef(1:n1,1:n2), ierr)

end subroutine sll_s_bspline_interpolation_2d_init
  
subroutine sll_s_compute_bspline_2d(self, gtau, &
  val1_min, val1_max, val2_min, val2_max)

  type(sll_t_bspline_interpolation_2d) :: self 
  sll_real64, intent(in)               :: gtau(:,:)
  sll_real64, intent(in), optional     :: val1_min(:,:)
  sll_real64, intent(in), optional     :: val1_max(:,:)
  sll_real64, intent(in), optional     :: val2_min(:,:)
  sll_real64, intent(in), optional     :: val2_max(:,:)

  sll_int32                            :: i
  sll_int32                            :: j

  if( present(val1_min) .and. present(val1_max)) then
    do j = 1, size(gtau,2)
      call sll_s_compute_bspline_1d( self%bs1, gtau(:,j), &
        val1_min(:,j), val1_max(:,j))
      self%bwork(j,:) = self%bs1%bcoef(:)
    end do

  else
    do j = 1, size(gtau,2)
      call sll_s_compute_bspline_1d( self%bs1, gtau(:,j))
      self%bwork(j,:) = self%bs1%bcoef(:)
    end do
  end if

  if( present(val2_min) .and. present(val2_max)) then
    do i = 1, size(self%bs1%bcoef)
      call sll_s_compute_bspline_1d( self%bs2, self%bwork(:,i), &
        val2_min(:,i), val2_max(:,i))
      self%bcoef(i,:) = self%bs2%bcoef(:)
    end do
  else
    do i = 1, size(self%bs1%bcoef)
      call sll_s_compute_bspline_1d( self%bs2, self%bwork(:,i))
      self%bcoef(i,:) = self%bs2%bcoef(:)
    end do
  end if

end subroutine sll_s_compute_bspline_2d

subroutine free_bspline_2D( spline )
  type(sll_t_bspline_interpolation_2d) :: spline
end subroutine free_bspline_2D 

subroutine sll_s_interpolate_array_values_2d(self, n1, n2, x1, x2, y )

type(sll_t_bspline_interpolation_2d) :: self
sll_int32                            :: n1
sll_int32                            :: n2
sll_real64, intent(in)               :: x1(:,:)
sll_real64, intent(in)               :: x2(:,:)
sll_real64, intent(out)              :: y(:,:)

sll_real64, allocatable :: work(:)
sll_int32               :: icell, jcell
sll_int32               :: i1, i2, k1, k2
sll_int32               :: ib, jb
sll_int32               :: j, k

SLL_ASSERT(n1 == size(x1,1) .and. n1 == size(x2,1) )
SLL_ASSERT(n2 == size(x1,2) .and. n2 == size(x2,2) )

k1 = self%bs1%deg+1
k2 = self%bs2%deg+1

allocate(work(k2))

do i2 = 1, n2
  do i1 = 1, n1
    icell = sll_f_find_cell( self%bs1%bsp, x1(i1,i2) )
    jcell = sll_f_find_cell( self%bs2%bsp, x2(i1,i2) )
    work = 0.0_f64
    do j = 1, k2
      jb = mod(jcell+j-2-self%bs2%offset+self%bs2%n,self%bs2%n) + 1
      call sll_s_splines_at_x(self%bs1%bsp, icell, x1(i1,i2), self%bs1%values)
      do k=1, k1
        ib = mod(icell+k-2-self%bs1%offset+self%bs1%n,self%bs1%n) + 1
        work(j) = work(j) + self%bs1%values(k)*self%bcoef(ib,jb)
      end do
    end do
    call sll_s_splines_at_x(self%bs2%bsp, jcell, x2(i1,i2), self%bs2%values)
    y(i1,i2) = 0.0_f64
    do k=1,k2
      jb = mod(jcell+k-2-self%bs2%offset+self%bs2%n,self%bs2%n) + 1
      y(i1,i2) = y(i1,i2) + self%bs2%values(k)*work(k)
    enddo
  end do
end do

deallocate(work)

end subroutine sll_s_interpolate_array_values_2d

function sll_f_interpolate_value_2d(self, xi, xj ) result (y)

type(sll_t_bspline_interpolation_2d) :: self
sll_real64, intent(in)               :: xi
sll_real64, intent(in)               :: xj
sll_real64                           :: y

sll_int32                            :: j
sll_int32                            :: k
sll_int32                            :: ib
sll_int32                            :: jb
sll_int32                            :: k1
sll_int32                            :: k2
sll_int32                            :: icell
sll_int32                            :: jcell
sll_real64, allocatable              :: work(:)

k1 = self%bs1%deg+1
k2 = self%bs2%deg+1
allocate(work(k2))
work = 0.0_f64

icell = sll_f_find_cell( self%bs1%bsp, xi )
jcell = sll_f_find_cell( self%bs2%bsp, xj )

do j = 1, k2
  jb = mod(jcell+j-2-self%bs2%offset+self%bs2%n,self%bs2%n) + 1
  call sll_s_splines_at_x(self%bs1%bsp, icell, xi, self%bs1%values)
  do k=1, k1
    ib = mod(icell+k-2-self%bs1%offset+self%bs1%n,self%bs1%n) + 1
    work(j) = work(j) + self%bs1%values(k)*self%bcoef(ib,jb)
  end do
end do

call sll_s_splines_at_x(self%bs2%bsp, jcell, xj, self%bs2%values)
y = 0.0_f64
do k=1,k2
  jb = mod(jcell+k-2-self%bs2%offset+self%bs2%n,self%bs2%n) + 1
  y = y + self%bs2%values(k)*work(k)
enddo

deallocate(work)

end function sll_f_interpolate_value_2d

function sll_f_interpolate_derivative_x1_2d(self, x1, x2 ) result(y)

type(sll_t_bspline_interpolation_2d) :: self
sll_real64, intent(in)               :: x1
sll_real64, intent(in)               :: x2
sll_real64                           :: y

sll_real64, allocatable :: work(:)
sll_int32               :: i, j
sll_int32               :: k1, k2, icell, jcell
sll_int32               :: ib, jb

k1 = self%bs1%deg+1
k2 = self%bs2%deg+1

icell =  sll_f_find_cell( self%bs1%bsp, x1 )
jcell =  sll_f_find_cell( self%bs2%bsp, x2 )

allocate(work(k2))
work = 0.0_f64
do j = 1, k2
  jb = mod(jcell+j-2-self%bs2%offset+self%bs2%n,self%bs2%n) + 1
  call sll_s_spline_derivatives_at_x(self%bs1%bsp, icell, x1, self%bs1%values)
  do i = 1, k1
    ib = mod(icell+i-2-self%bs1%offset+self%bs1%n,self%bs1%n) + 1
    work(j) = work(j) + self%bs1%values(i)*self%bcoef(ib,jb)
  enddo
end do

call sll_s_splines_at_x(self%bs2%bsp, jcell, x2, self%bs2%values)
y = 0.0_f64
do j = 1, k2
  jb = mod(jcell+j-2-self%bs2%offset+self%bs2%n,self%bs2%n) + 1
  y = y + self%bs2%values(j)*work(j)
enddo

deallocate(work)

end function sll_f_interpolate_derivative_x1_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function sll_f_interpolate_derivative_x2_2d(self, x1, x2 ) result(y)

type(sll_t_bspline_interpolation_2d) :: self
sll_real64, intent(in)               :: x1
sll_real64, intent(in)               :: x2
sll_real64                           :: y

sll_real64, allocatable :: work(:)
sll_int32               :: i, j
sll_int32               :: k1, k2, icell, jcell
sll_int32               :: ib, jb

k1 = self%bs1%deg+1
k2 = self%bs2%deg+1

icell =  sll_f_find_cell( self%bs1%bsp, x1 )
jcell =  sll_f_find_cell( self%bs2%bsp, x2 )

allocate(work(k2))
work = 0.0_f64
do j = 1, k2
  jb = mod(jcell+j-2-self%bs2%offset+self%bs2%n,self%bs2%n) + 1
  call sll_s_splines_at_x(self%bs1%bsp, icell, x1, self%bs1%values)
  do i = 1, k1
    ib = mod(icell+i-2-self%bs1%offset+self%bs1%n,self%bs1%n) + 1
    work(j) = work(j) + self%bs1%values(i)*self%bcoef(ib,jb)
  enddo
end do

call sll_s_spline_derivatives_at_x(self%bs2%bsp, jcell, x2, self%bs2%values)
y = 0.0_f64
do j = 1, k2
  jb = mod(jcell+j-2-self%bs2%offset+self%bs2%n,self%bs2%n) + 1
  y = y + self%bs2%values(j)*work(j)
enddo

deallocate(work)

end function sll_f_interpolate_derivative_x2_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine sll_s_bspline_interpolation_2d_free(self ) 

type(sll_t_bspline_interpolation_2d)    :: self


end subroutine sll_s_bspline_interpolation_2d_free

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine sll_s_interpolate_array_derivatives_x1_2d( self, n1, n2, x1, x2, y)

type(sll_t_bspline_interpolation_2d) :: self
sll_int32                            :: n1
sll_int32                            :: n2
sll_real64                           :: x1(:,:)
sll_real64                           :: x2(:,:)
sll_real64                           :: y(:,:)

sll_int32 :: i1, i2
sll_real64, allocatable :: work(:)
sll_int32               :: i, j
sll_int32               :: k1, k2, icell, jcell
sll_int32               :: ib, jb
sll_real64              :: xi, xj, yy

k1 = self%bs1%deg+1
k2 = self%bs2%deg+1
allocate(work(k2))

do i2 = 1, n2
  do i1 = 1, n1

    xi = x1(i1,i2)
    xj = x2(i1,i2)
    icell =  sll_f_find_cell( self%bs1%bsp, xi )
    jcell =  sll_f_find_cell( self%bs2%bsp, xj )
    
    work = 0.0_f64
    do j = 1, k2
      jb = mod(jcell+j-2-self%bs2%offset+self%bs2%n,self%bs2%n) + 1
      call sll_s_spline_derivatives_at_x(self%bs1%bsp, icell, xi, self%bs1%values)
      do i = 1, k1
        ib = mod(icell+i-2-self%bs1%offset+self%bs1%n,self%bs1%n) + 1
        work(j) = work(j) + self%bs1%values(i)*self%bcoef(ib,jb)
      enddo
    end do
    call sll_s_splines_at_x(self%bs2%bsp, jcell, xj, self%bs2%values)
    yy = 0.0_f64
    do j = 1, k2
      jb = mod(jcell+j-2-self%bs2%offset+self%bs2%n,self%bs2%n) + 1
      yy = yy + self%bs2%values(j)*work(j)
    enddo
    
    y(i1,i2) = yy
    
  end do
end do

deallocate(work)

end subroutine sll_s_interpolate_array_derivatives_x1_2d

subroutine sll_s_interpolate_array_derivatives_x2_2d( self, n1, n2, x1, x2, y)

type(sll_t_bspline_interpolation_2d) :: self
sll_int32                            :: n1
sll_int32                            :: n2
sll_real64                           :: x1(:,:)
sll_real64                           :: x2(:,:)
sll_real64                           :: y(:,:)

sll_real64, allocatable :: work(:)
sll_int32               :: i, j, i1, i2
sll_int32               :: k1, k2, icell, jcell
sll_int32               :: ib, jb
sll_real64              :: xi, xj, yy

k1 = self%bs1%deg+1
k2 = self%bs2%deg+1
allocate(work(k2))

do i2 = 1, n2
  do i1 = 1, n1
    xi = x1(i1,i2)
    xj = x2(i1,i2)
    icell =  sll_f_find_cell( self%bs1%bsp, xi )
    jcell =  sll_f_find_cell( self%bs2%bsp, xj )
    
    work = 0.0_f64
    do j = 1, k2
      jb = mod(jcell+j-2-self%bs2%offset+self%bs2%n,self%bs2%n) + 1
      call sll_s_splines_at_x(self%bs1%bsp, icell, xi, self%bs1%values)
      do i = 1, k1
        ib = mod(icell+i-2-self%bs1%offset+self%bs1%n,self%bs1%n) + 1
        work(j) = work(j) + self%bs1%values(i)*self%bcoef(ib,jb)
      enddo
    end do
    
    call sll_s_spline_derivatives_at_x(self%bs2%bsp, jcell, xj, self%bs2%values)
    yy = 0.0_f64
    do j = 1, k2
      jb = mod(jcell+j-2-self%bs2%offset+self%bs2%n,self%bs2%n) + 1
      yy = yy + self%bs2%values(j)*work(j)
    enddo

    y(i1,i2) = yy

  end do
end do

deallocate(work)

end subroutine sll_s_interpolate_array_derivatives_x2_2d

end module sll_m_bspline_interpolation
