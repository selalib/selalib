!> @ingroup splines
!> @brief   Interpolator for 2D polar splines of arbitrary degree,
!>          on uniform and non-uniform grids (directions are independent)
!> @author  Yaman Güçlü  - IPP Garching

module sll_m_polar_spline_interpolator_2d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"

  use sll_m_working_precision, only: f64

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_periodic, &
    sll_p_hermite, &
    sll_p_greville

  use sll_m_bsplines, only: &
    sll_c_bsplines, &
    sll_s_bsplines_new_mirror_copy

  use sll_m_spline_1d, only: sll_t_spline_1d
  use sll_m_spline_2d, only: sll_t_spline_2d

  use sll_m_spline_interpolator_1d, only: &
    sll_t_spline_interpolator_1d, &
    sll_s_spline_1d_compute_num_cells

  use sll_m_spline_interpolator_2d, only: &
    sll_t_spline_interpolator_2d, &
    sll_t_spline_2d_boundary_data, &
    sll_s_spline_2d_compute_num_cells

  implicit none

  public :: &
    sll_t_polar_spline_interpolator_2d, &
    sll_s_polar_spline_2d_compute_num_cells

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Working precision
  integer, parameter :: wp = f64

  !-----------------------------------------------------------------------------
  !> 2D tensor-product spline interpolator
  !-----------------------------------------------------------------------------
  type :: sll_t_polar_spline_interpolator_2d

    private

    class(sll_c_bsplines), pointer :: bspl1 => null()
    class(sll_c_bsplines), pointer :: bspl2 => null()

    integer                            :: nbc_rmax
    integer                            :: nipts(2)
    integer                            :: i1_start
    integer                            :: even_deg1
    real(wp), allocatable              :: ext_gtau(:,:)
    type(sll_t_spline_2d)              :: ext_spline
    type(sll_t_spline_interpolator_2d) :: ext_interp
    class(sll_c_bsplines), allocatable :: ext_bspl1

  contains

    procedure :: init                => s_polar_spline_interpolator_2d__init
    procedure :: free                => s_polar_spline_interpolator_2d__free
    procedure :: get_interp_points   => s_polar_spline_interpolator_2d__get_interp_points
    procedure :: compute_interpolant => s_polar_spline_interpolator_2d__compute_interpolant

  end type sll_t_polar_spline_interpolator_2d

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !-----------------------------------------------------------------------------
  !> @brief      Calculate number of cells from number of interpolation points
  !> @details    Important for parallelization: for given spline degree and BCs,
  !>             calculate the numbers of grid cells along x1 and x2 that yield
  !>             the desired number of interpolation points along x1 and x2
  !>
  !> @param[in]  degree   spline degrees along x1 and x2
  !> @param[in]  bc_xmin  boundary conditions at x1_min and x2_min
  !> @param[in]  bc_xmax  boundary conditions at x1_max and x2_max
  !> @param[in]  nipts    desired number of interpolation points along x1 and x2
  !> @param[out] ncells   calculated number of grid cells along x1 and x2
  !-----------------------------------------------------------------------------
  subroutine sll_s_polar_spline_2d_compute_num_cells( &
      degree , &
      bc_rmax, &
      nipts  , &
      ncells )

    integer, intent(in   ) :: degree (2)
    integer, intent(in   ) :: bc_rmax
    integer, intent(in   ) :: nipts  (2)
    integer, intent(  out) :: ncells (2)

    ! Along radial direction r, on extended domain [-R,R]
    !
    ! Option 1: breakpoint at r=0
!    call sll_s_spline_1d_compute_num_cells( &
!      degree(1), bc_rmax, bc_rmax, 2*nipts(1)-1, ncells(1) )
    !
    ! Option 2: no breakpoint at r=0
!    call sll_s_spline_1d_compute_num_cells( &
!      degree(1), bc_rmax, bc_rmax, 2*nipts(1), ncells(1) )
    !
    ! Option 3: no breakpoint at r=0, and use radial basis
    call sll_s_spline_1d_compute_num_cells( &
      degree(1), bc_rmax, bc_rmax, 2*nipts(1), ncells(1) )
    ncells(1) = (ncells(1)+1)/2

    ! Along angular direction theta
    call sll_s_spline_1d_compute_num_cells( &
      degree(2), sll_p_periodic, sll_p_periodic, nipts(2), ncells(2) )

  end subroutine sll_s_polar_spline_2d_compute_num_cells

  !-----------------------------------------------------------------------------
  !> @brief      Initialize a 2D tensor-product spline interpolator object
  !> @param[out] self     2D tensor-product spline interpolator
  !> @param[in]  bspl1    B-splines (basis) along x1 direction
  !> @param[in]  bspl2    B-splines (basis) along x2 direction
  !> @param[in]  bc_xmin  boundary conditions at x1_min and x2_min
  !> @param[in]  bc_xmax  boundary conditions at x1_max and x2_max
  !-----------------------------------------------------------------------------
  subroutine s_polar_spline_interpolator_2d__init( &
    self   , &
    bspl1  , &
    bspl2  , &
    bc_rmax )

    class(sll_t_polar_spline_interpolator_2d), intent(  out) :: self
    class(sll_c_bsplines),       target, intent(in   ) :: bspl1
    class(sll_c_bsplines),       target, intent(in   ) :: bspl2
    integer                            , intent(in   ) :: bc_rmax

    ! Check requirements on radial basis
    SLL_ASSERT( bspl1 % radial )
    SLL_ASSERT( .not. bspl1 % periodic )

    ! Check requirements on angular basis
    SLL_ASSERT( bspl2 % periodic )
    SLL_ASSERT( modulo( bspl2 % ncells, 2 ) == 0 )

    ! Check boundary condition
    SLL_ASSERT( any( bc_rmax == [sll_p_greville, sll_p_hermite] ) )

    ! Store pointers to B-splines
    self % bspl1 => bspl1
    self % bspl2 => bspl2

    ! Create B-spline on extended radial domain [-R,R]
    call sll_s_bsplines_new_mirror_copy( self % bspl1, self % ext_bspl1 )

    ! Initialize 2D spline with polar base
    call self % ext_spline % init( self % ext_bspl1, self % bspl2 )

    ! Initialize 2D interpolator on extended domain [-R,R] X [0,2\pi]
    call self % ext_interp % init( &
      self % ext_bspl1, &
      self %     bspl2, &
      bc_xmin = [bc_rmax, sll_p_periodic], &
      bc_xmax = [bc_rmax, sll_p_periodic] )

    ! Determine number of additional boundary conditions at r=r_max (Hermite)
    self % nbc_rmax = merge( self%bspl1%degree/2, 0, bc_rmax == sll_p_hermite )

    ! Allocate extended radial array of function values
    associate( n1 => self % ext_bspl1 % nbasis, &
               n2 => self %     bspl2 % nbasis, &
               a1 => self % nbc_rmax )

      allocate( self % ext_gtau (n1-2*a1, n2) )

    end associate

    ! Determine size of 2D grid of interpolation points in physical domain
    self % nipts(1) = (size( self%ext_gtau, 1 ) + 1) / 2 
    self % nipts(2) =  size( self%ext_gtau, 2 )

    ! Determine first physical point in extended radial domain
    self % i1_start = size( self%ext_gtau, 1 ) - self%nipts(1) + 1

    ! Store info about whether spline degree along r is even or odd
    self % even_deg1 = 1 - modulo( self%bspl1%degree, 2 )

  end subroutine s_polar_spline_interpolator_2d__init

  !-----------------------------------------------------------------------------
  !> @brief        Destroy local objects and free allocated memory
  !> @param[inout] self  2D tensor-product spline interpolator
  !-----------------------------------------------------------------------------
  subroutine s_polar_spline_interpolator_2d__free( self )

    class(sll_t_polar_spline_interpolator_2d), intent(inout) :: self

    ! Destroy local objects: 1D splines and interpolators along x1 and x2
    call self % ext_interp % free()
    call self % ext_bspl1  % free()

    ! Deallocate 2D array of interpolation points
    deallocate( self % ext_gtau )

  end subroutine s_polar_spline_interpolator_2d__free

  !-----------------------------------------------------------------------------
  !> @brief      Get coordinates of interpolation points (2D tensor grid)
  !> @param[in]  self  2D tensor-product spline interpolator
  !> @param[out] tau1  x1 coordinates of interpolation points
  !> @param[out] tau2  x2 coordinates of interpolation points
  !-----------------------------------------------------------------------------
  subroutine s_polar_spline_interpolator_2d__get_interp_points( self, tau1, tau2 )

    class(sll_t_polar_spline_interpolator_2d), intent(in   ) :: self
    real(wp),                     allocatable, intent(  out) :: tau1(:)
    real(wp),                     allocatable, intent(  out) :: tau2(:)

    real(wp), allocatable :: ext_tau1(:)
    real(wp), allocatable :: ext_tau2(:)

    ! Obtain interpolation points on extended domain
    call self % ext_interp % get_interp_points( ext_tau1, ext_tau2 )

    ! Determine interpolation points on physical domain
    associate( n1    => self%nipts(1), &
               n2    => self%nipts(2), &
               first => self%i1_start  )

      allocate( tau1(n1), source=ext_tau1(first:) )
      allocate( tau2(n2), source=ext_tau2(:)      )

    end associate

  end subroutine s_polar_spline_interpolator_2d__get_interp_points

  !-----------------------------------------------------------------------------
  !> @brief        Compute interpolating 2D spline
  !> @details      Compute coefficients of 2D tensor-product spline that
  !>               interpolates function values on grid. If Hermite BCs are used,
  !>               function derivatives at appropriate boundaries are also needed.
  !>
  !> @param[inout] self           2D tensor-product spline interpolator
  !> @param[inout] spline         2D tensor-product spline
  !> @param[in]    gtau           function values of interpolation points
  !> @param[in]    boundary_data  (optional) structure with boundary conditions
  !-----------------------------------------------------------------------------
  subroutine s_polar_spline_interpolator_2d__compute_interpolant( self, &
      spline, gtau, derivs_rmax )

    class(sll_t_polar_spline_interpolator_2d), intent(inout)           :: self
    type (sll_t_spline_2d)                   , intent(inout)           :: spline
    real(wp)                                 , intent(in   )           :: gtau(:,:)
    real(wp)                                 , intent(in   ), optional :: derivs_rmax(:,:)

    character(len=*), parameter :: &
      this_sub_name = "sll_t_polar_spline_interpolator_2d % compute_interpolant"

    class(sll_t_spline_2d_boundary_data), allocatable :: boundary_data
    integer :: i1, i2, j1, j2, k, shift

    SLL_ASSERT( spline % belongs_to_space( self % bspl1, self % bspl2 ) )
    SLL_ASSERT( size( gtau, 1 ) == self % nipts(1) )
    SLL_ASSERT( size( gtau, 2 ) == self % nipts(2) )

    ! Copy physical domain
    self % ext_gtau (self%i1_start:,:) = gtau(:,:)

    shift = self % nipts(2) / 2

    ! Fill in r<0 region with mirrored data
    do i2 = 1, self % nipts(2)
      j2 = 1 + modulo( i2-1+shift, self%nipts(2) )
      do j1 = 1, self%i1_start-1
        i1 = self%i1_start-j1
        self % ext_gtau(j1,j2) = gtau(i1,i2)
      end do
    end do

    ! Boundary data
    ! TODO: allocate in init()
    if (present( derivs_rmax )) then
      SLL_ASSERT( size( derivs_rmax, 1 ) == self%nbc_rmax )
      SLL_ASSERT( size( derivs_rmax, 2 ) == self%nipts(2) )

      allocate( boundary_data )
      allocate( boundary_data % derivs_x1_min (self%nbc_rmax, self%nipts(2)) )
      allocate( boundary_data % derivs_x1_max (self%nbc_rmax, self%nipts(2)) )

      ! Derivatives at r=r_max
      boundary_data % derivs_x1_max(:,:) = derivs_rmax(:,:)

      ! Derivatives at r=r_min: must apply shift in theta and change sign if odd
      do i2 = 1, self % nipts(2)
        j2 = 1 + modulo( i2-1+shift, self%nipts(2) )
        do k = 1, self%nbc_rmax
          boundary_data % derivs_x1_min(k,j2) = (-1)**(k-self%even_deg1) * derivs_rmax(k,i2)
        end do
      end do
    end if

    ! Compute interpolant on extended domain
    call self % ext_interp % compute_interpolant( self%ext_spline, self%ext_gtau, boundary_data )

    ! Copy data back onto 2D polar spline
    associate( start => 1 + size( self%ext_spline%bcoef, 1 ) - size( spline%bcoef, 1 ) )

      spline % bcoef(:,:) = self % ext_spline % bcoef(start:,:)

    end associate

  end subroutine s_polar_spline_interpolator_2d__compute_interpolant

end module sll_m_polar_spline_interpolator_2d
