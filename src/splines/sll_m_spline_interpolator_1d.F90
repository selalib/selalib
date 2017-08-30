module sll_m_spline_interpolator_1d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"

  use sll_m_working_precision, only: f64

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_periodic, &
    sll_p_hermite , &
    sll_p_greville

  use sll_m_bsplines_base, only: &
    sll_c_bsplines

  use sll_m_spline_1d_new, only: &
    sll_t_spline_1d

  use sll_m_spline_matrix, only: &
    sll_c_spline_matrix, &
    sll_s_spline_matrix_new

  implicit none

  public :: &
    sll_t_spline_interpolator_1d

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Working precision
  integer, parameter :: wp = f64

  !> Allowed boundary conditions
  integer, parameter :: &
    allowed_bcs(*) = [sll_p_periodic, sll_p_hermite, sll_p_greville]

  !> 1D spline interpolator
  type :: sll_t_spline_interpolator_1d

    ! Private attributes
    class(sll_c_bsplines),          pointer, private :: bspl => null()
    integer                                , private ::  bc_xmin
    integer                                , private ::  bc_xmax
    integer                                , private :: nbc_xmin
    integer                                , private :: nbc_xmax
    integer                                , private :: odd
    real(wp)                               , private :: dx
    real(wp),                   allocatable, private :: tau(:)
    class(sll_c_spline_matrix), allocatable, private :: matrix

  contains

    procedure :: init                => s_spline_interpolator_1d__init
    procedure :: free                => s_spline_interpolator_1d__free
    procedure :: get_interp_points   => s_spline_interpolator_1d__get_interp_points
    procedure :: compute_interpolant => s_spline_interpolator_1d__compute_interpolant

  end type sll_t_spline_interpolator_1d

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine s_spline_interpolator_1d__init( self, bspl, bc_xmin, bc_xmax )

    class(sll_t_spline_interpolator_1d), intent(  out) :: self
    class(sll_c_bsplines),       target, intent(in   ) :: bspl
    integer                            , intent(in   ) :: bc_xmin
    integer                            , intent(in   ) :: bc_xmax

    integer :: kl, ku
    character(len=32) :: matrix_type

    ! Sanity checks
    SLL_ASSERT( any( bc_xmin == allowed_bcs ) )
    SLL_ASSERT( any( bc_xmax == allowed_bcs ) )

    ! Save pointer to B-splines
    ! (later needed to verify 1D spline input to 'compute_interpolant')
    self % bspl => bspl

    ! Save other data
    self%bc_xmin  = bc_xmin
    self%bc_xmax  = bc_xmax
    self%nbc_xmin = merge( bspl%degree/2, 0, bc_xmin == sll_p_hermite )
    self%nbc_xmax = merge( bspl%degree/2, 0, bc_xmax == sll_p_hermite )
    self%odd      = modulo( bspl%degree, 2 )

    ! Save average cell size for normalization of derivatives
    self%dx = (bspl%xmax - bspl%xmin) / bspl%ncells

    ! Compute interpolation points and number of diagonals in linear system
    ! TODO: implement subroutines below
    if (bspl % uniform) then
      call s_compute_interpolation_points_uniform( self, self%tau )
      call s_compute_num_diags_uniform( self, kl, ku )
    else
      call s_compute_interpolation_points_non_uniform( self, self%tau )
      call s_compute_num_diags_non_uniform( self, kl, ku )
    end if

    ! Determine matrix storage type and allocate matrix
    associate( n => self % bspl % nbasis )

      if (bc_xmin == sll_p_periodic) then
        if (kl+1+ku >= n) then
          matrix_type = "dense"
        else
          matrix_type = "periodic_banded"
        end if
      else
        matrix_type = "banded"
      end if

      call sll_s_spline_matrix_new( self%matrix, matrix_type, n, kl, ku )

    end associate ! n

    ! Fill in entries A_ij of linear system A*x = b for interpolation
    call s_build_system( self, self%matrix )

    ! Factorize matrix A to speed-up solution for any given b
    call self % matrix % factorize()

  end subroutine s_spline_interpolator_1d__init

  !-----------------------------------------------------------------------------
  subroutine s_spline_interpolator_1d__free( self )

    class(sll_t_spline_interpolator_1d), intent(inout) :: self

    nullify   ( self % bspl )
    deallocate( self % tau  )

  end subroutine s_spline_interpolator_1d__free

  !-----------------------------------------------------------------------------
  subroutine s_spline_interpolator_1d__get_interp_points( self, tau )

    class(sll_t_spline_interpolator_1d), intent(in   ) :: self
    real(wp),               allocatable, intent(  out) :: tau(:)

    SLL_ASSERT( allocated( self%tau ) )
    allocate( tau(size(self%tau)), source=self%tau )

  end subroutine s_spline_interpolator_1d__get_interp_points

  !-----------------------------------------------------------------------------
  !> @brief        Compute interpolating 1D spline
  !> Computes coefficients of 1D spline that interpolates function on grid.
  !> If Hermite BCs are used, derivatives at boundary are also needed.
  !>
  !> @param[inout] self       1D spline interpolator object
  !> @param[inout] spline     1D spline
  !> @param[in]    gtau       function values of interpolation points
  !> @param[in]    derivs_xmin (optional) array with boundary conditions at xmin
  !> @param[in]    derivs_xmax (optional) array with boundary conditions at xmax
  !-----------------------------------------------------------------------------
  subroutine s_spline_interpolator_1d__compute_interpolant( self, &
      spline, gtau, derivs_xmin, derivs_xmax )

    class(sll_t_spline_interpolator_1d), intent(in   ) :: self
    type (sll_t_spline_1d)             , intent(inout) :: spline
    real(wp)                           , intent(in   ) :: gtau(:)
    real(wp),                  optional, intent(in   ) :: derivs_xmin(:)
    real(wp),                  optional, intent(in   ) :: derivs_xmax(:)

    character(len=*), parameter :: &
      this_sub_name = "sll_t_spline_interpolator_1d % compute_interpolant"

    integer :: i

    SLL_ASSERT( size(gtau) == self%bspl%nbasis - self%nbc_xmin - self%nbc_xmax )
    SLL_ASSERT( spline % belongs_to_space( self % bspl ) )

    associate ( n        => self%bspl%nbasis, &
                deg      => self%bspl%degree, &
                nbc_xmin => self%nbc_xmin   , &
                nbc_xmax => self%nbc_xmax   , &
                bcoef    => spline%bcoef    )

      ! Special case: linear spline
      if (deg == 1) then
        bcoef(1:n) = gtau(1:n)
        ! Periodic only: "wrap around" coefficients onto extended array
        if (self%bspl%periodic) then
          bcoef(0)   = bcoef(n)
          bcoef(n+1) = bcoef(1)
        end if
        return
      end if

      ! Hermite boundary conditions at xmin, if any
      ! NOTE: For consistency with the linear system, the i-th derivative
      !       provided by the user must be multiplied by dx^i
      if ( self%bc_xmin == sll_p_hermite ) then
        if ( present( derivs_xmin ) ) then
          bcoef(1:nbc_xmin) = &
                     [(derivs_xmin(i)*self%dx**(i+self%odd-1), i=nbc_xmin,1,-1)]
        else
          SLL_ERROR(this_sub_name,"Hermite BC at xmin requires derivatives")
        end if
      end if

      ! Interpolation points
      bcoef(nbc_xmin+1:n-nbc_xmax) = gtau(:)

      ! Hermite boundary conditions at xmax, if any
      ! NOTE: For consistency with the linear system, the i-th derivative
      !       provided by the user must be multiplied by dx^i
      if ( self%bc_xmax == sll_p_hermite ) then
        if ( present( derivs_xmax ) ) then
          bcoef(n-nbc_xmax+1:n) = &
                        [(derivs_xmax(i)*self%dx**(i+self%odd-1), i=1,nbc_xmax)]
        else
          SLL_ERROR(this_sub_name,"Hermite BC at xmax requires derivatives")
        end if
      end if

      ! Solve linear system and compute coefficients
      call self % matrix % solve_inplace( bcoef(1:n) )

      ! Periodic only: "wrap around" coefficients onto extended array
      if (self%bc_xmin == sll_p_periodic) then
        associate( g => 1+deg/2 )
          bcoef(1-g:0)   = bcoef(n-g+1:n)
          bcoef(n+1:n+g) = bcoef(1:g)
        end associate
      end if

    end associate ! n, deg, bcoef

  end subroutine s_spline_interpolator_1d__compute_interpolant

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !!                            PRIVATE SUBROUTINES
  !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! TODO: move subroutine to init()
  !-----------------------------------------------------------------------------
  !> @brief        Private subroutine for assembling and factorizing linear
  !>               system for any kind of boundary conditions at xmin and xmax
  !> @param[in]    self   1D spline interpolator object
  !> @param[inout] matrix generic 'spline' matrix (dense/banded/periodic-banded)
  !-----------------------------------------------------------------------------
  subroutine s_build_system( self, matrix )

    class(sll_t_spline_interpolator_1d), intent(in   ) :: self
    class(sll_c_spline_matrix)         , intent(inout) :: matrix

    integer  :: i,j,d,s
    integer  :: j0,d0
    integer  :: order
    real(wp) :: x
    real(wp) :: values(self%bspl%degree+1)
    real(wp), allocatable :: derivs(:,:)

    ! NEW
    integer :: jmin

    associate( n        => self%bspl%nbasis, &
               deg      => self%bspl%degree, &
               nbc_xmin => self%nbc_xmin   , &
               nbc_xmax => self%nbc_xmax   )

      if ( any( [self%bc_xmin,self%bc_xmax] == sll_p_hermite ) ) then
        allocate ( derivs (0:deg/2, 1:deg+1) )
      end if

      ! Hermite boundary conditions at xmin, if any
      if ( self%bc_xmin == sll_p_hermite ) then
        x = self%bspl%xmin
        call self % bspl % eval_basis_and_n_derivs( x, nbc_xmin, derivs, jmin )

        ! In order to improve the condition number of the matrix, we normalize all
        ! derivatives by multiplying the i-th derivative by dx^i
        associate( h => [(self%dx**i, i=1, ubound(derivs,1))] )
          do j = lbound(derivs,2), ubound(derivs,2)
            derivs(1:,j) = derivs(1:,j) * h(1:)
          end do
        end associate

        do i = 1, nbc_xmin
          ! iterate only to deg as last bspline is 0
          order = nbc_xmin-i+self%odd
          do j = 1, deg
            call matrix % set_element( i, j, derivs(order,j) )
          end do
        end do

      end if

      ! Interpolation points
      do i = nbc_xmin+1, n-nbc_xmax
        x = self%tau(i-nbc_xmin)
        call self % bspl % eval_basis( x, values, jmin )
        do s = 1, deg+1
          j = modulo( jmin+s-1, n ) + 1 
          call matrix % set_element( i, j, values(s) )
        end do
      end do

      ! Hermite boundary conditions at xmax, if any
      if ( self%bc_xmax == sll_p_hermite ) then
        x = self%bspl%xmax
        call self % bspl % eval_basis_and_n_derivs( x, nbc_xmin, derivs, jmin )

        ! In order to improve the condition number of the matrix, we normalize all
        ! derivatives by multiplying the i-th derivative by dx^i
        associate( h => [(self%dx**i, i=1, ubound(derivs,1))] )
          do j = lbound(derivs,2), ubound(derivs,2)
            derivs(1:,j) = derivs(1:,j) * h(1:)
          end do
        end associate

        do i = n-nbc_xmax+1, n
          order = i-(n-nbc_xmax+1)+self%odd
          j0 = n-deg
          d0 = 1
          do s = 1, deg
            j = j0 + s
            d = d0 + s
            call matrix % set_element( i, j, derivs(order,d) )
          end do
        end do

      end if

      if ( allocated( derivs ) ) deallocate( derivs )

    end associate ! n, deg

  end subroutine s_build_system

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !!                       PRIVATE SUBROUTINES, UNIFORM
  !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !-----------------------------------------------------------------------------
  ! TODO: implement subroutine
  subroutine s_compute_interpolation_points_uniform( self, tau )
    class(sll_t_spline_interpolator_1d), intent(in   ) :: self
    real(wp),               allocatable, intent(  out) :: tau(:)
  end subroutine s_compute_interpolation_points_uniform

  !-----------------------------------------------------------------------------
  ! TODO: implement subroutine
  subroutine s_compute_num_diags_uniform( self, kl, ku )
    class(sll_t_spline_interpolator_1d), intent(in   ) :: self
    integer                            , intent(  out) :: kl
    integer                            , intent(  out) :: ku
  end subroutine s_compute_num_diags_uniform

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !!                       PRIVATE SUBROUTINES, NON-UNIFORM
  !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !-----------------------------------------------------------------------------
  ! TODO: implement subroutine
  subroutine s_compute_interpolation_points_non_uniform( self, tau )
    class(sll_t_spline_interpolator_1d), intent(in   ) :: self
    real(wp),               allocatable, intent(  out) :: tau(:)
  end subroutine s_compute_interpolation_points_non_uniform

  !-----------------------------------------------------------------------------
  ! TODO: implement subroutine
  subroutine s_compute_num_diags_non_uniform( self, kl, ku )
    class(sll_t_spline_interpolator_1d), intent(in   ) :: self
    integer                            , intent(  out) :: kl
    integer                            , intent(  out) :: ku
  end subroutine s_compute_num_diags_non_uniform

end module sll_m_spline_interpolator_1d
