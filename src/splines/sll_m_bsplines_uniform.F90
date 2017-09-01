module sll_m_bsplines_uniform
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"

  use sll_m_working_precision, only: f64
  use sll_m_bsplines_base    , only: sll_c_bsplines

  implicit none

  public :: &
    sll_t_bsplines_uniform

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Working precision
  integer, parameter :: wp = f64

  !> Abstract type, B-splines
  type, extends(sll_c_bsplines) :: sll_t_bsplines_uniform

    ! Private components (public ones defined in base class)
    real(wp), private :: inv_dx

  contains
    procedure :: init                    => s_bsplines_uniform__init
    procedure :: free                    => s_bsplines_uniform__free
    procedure :: find_cell               => f_bsplines_uniform__find_cell
    procedure :: eval_basis              => s_bsplines_uniform__eval_basis
    procedure :: eval_deriv              => s_bsplines_uniform__eval_deriv
    procedure :: eval_basis_and_n_derivs => s_bsplines_uniform__eval_basis_and_n_derivs

  end type sll_t_bsplines_uniform

  ! Inverse of integers for later use (max spline degree = 32)
  integer             :: index
  real(wp), parameter :: inv(*) = [(1.0_wp/real(index,wp), index=1,32)]

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Module subroutine, private
  SLL_PURE subroutine s_bsplines_uniform__get_icell_and_offset( self, x, icell, x_offset )
    class(sll_t_bsplines_uniform), intent(in   ) :: self
    real(wp)                     , intent(in   ) :: x
    integer                      , intent(  out) :: icell
    real(wp)                     , intent(  out) :: x_offset

    SLL_ASSERT( x >= self % xmin )
    SLL_ASSERT( x <= self % xmax )

    if      (x == self%xmin) then; icell = 1          ; x_offset = 0.0_wp
    else if (x == self%xmax) then; icell = self%ncells; x_offset = 1.0_wp
    else
      x_offset = (x-self%xmin) * self%inv_dx  ! 0 <= x_offset <= num_cells
      icell    = int( x_offset )
      x_offset = x_offset - real( icell, wp ) ! 0 <= x_offset < 1
      icell    = icell + 1
    end if

  end subroutine s_bsplines_uniform__get_icell_and_offset

  !-----------------------------------------------------------------------------
  subroutine s_bsplines_uniform__init( self, degree, periodic, xmin, xmax, ncells )
    class(sll_t_bsplines_uniform), intent(  out) :: self
    integer                      , intent(in   ) :: degree
    logical                      , intent(in   ) :: periodic
    real(wp)                     , intent(in   ) :: xmin
    real(wp)                     , intent(in   ) :: xmax
    integer                      , intent(in   ) :: ncells

    ! Public attributes
    self % degree   = degree
    self % periodic = periodic
    self % uniform  = .true.
    self % ncells   = ncells
    self % nbasis   = merge( ncells  , ncells+degree, periodic )
    self % offset   = merge( degree/2, 0            , periodic )
    self % xmin     = xmin
    self % xmax     = xmax

    ! Private attributes
    self % inv_dx   = real(ncells,wp) / (xmax-xmin)

  end subroutine s_bsplines_uniform__init

  !-----------------------------------------------------------------------------
  !> Free storage
  subroutine s_bsplines_uniform__free( self )
    class(sll_t_bsplines_uniform), intent(inout) :: self
  end subroutine s_bsplines_uniform__free

  !----------------------------------------------------------------------------
  function f_bsplines_uniform__find_cell( self, x ) result( icell )
    class(sll_t_bsplines_uniform), intent(in) :: self
    real(wp)                     , intent(in) :: x
    integer :: icell

    icell = -1 ! to avoid warning at compile time
    SLL_ERROR("sll_t_bsplines_uniform % find_cell","procedure not implemented")

  end function f_bsplines_uniform__find_cell

  !-----------------------------------------------------------------------------
  !> Evaluate value at x of all basis functions with support in local cell
  !> (jrange identifies index of basis functions)
  SLL_PURE subroutine s_bsplines_uniform__eval_basis( self, x, values, jmin )
    class(sll_t_bsplines_uniform), intent(in   ) :: self
    real(wp)                     , intent(in   ) :: x
    real(wp)                     , intent(  out) :: values(0:)
    integer                      , intent(  out) :: jmin

    integer  :: icell
    real(wp) :: x_offset

    integer  :: j, r
    real(wp) :: j_real
    real(wp) :: xx
    real(wp) :: temp
    real(wp) :: saved

    ! Check on inputs
    SLL_ASSERT( size(values) == 1+self%degree )

    ! 1. Compute cell index 'icell' and x_offset
    call s_bsplines_uniform__get_icell_and_offset( self, x, icell, x_offset )

    ! 2. Compute index range of B-splines with support over cell 'icell'
    jmin = icell - self%offset

    ! 3. Compute values of aforementioned B-splines
    associate( bspl => values, spline_degree => self%degree )

      bspl(0) = 1.0_wp
      do j = 1, spline_degree
         j_real = real(j,wp)
         xx     = -x_offset
         saved  = 0.0_wp
         do r = 0, j-1
            xx      = xx + 1.0_wp
            temp    = bspl(r) * inv(j)
            bspl(r) = saved + xx * temp
            saved   = (j_real - xx) * temp
         end do
         bspl(j) = saved
      end do

    end associate  ! bspl, spline_degree

  end subroutine s_bsplines_uniform__eval_basis

  !-----------------------------------------------------------------------------
  !> Evaluate derivative at x of all basis functions with support in local cell
  SLL_PURE subroutine s_bsplines_uniform__eval_deriv( self, x, derivs, jmin )
    class(sll_t_bsplines_uniform), intent(in   ) :: self
    real(wp)                     , intent(in   ) :: x
    real(wp)                     , intent(  out) :: derivs(0:)
    integer                      , intent(  out) :: jmin

    integer  :: icell
    real(wp) :: x_offset

    integer  :: j, r
    real(wp) :: j_real
    real(wp) :: xx
    real(wp) :: temp
    real(wp) :: saved

    real(wp) :: bj, bjm1 

    ! Check on inputs
    SLL_ASSERT( size(derivs) == 1+self%degree )

    ! 1. Compute cell index 'icell' and x_offset
    call s_bsplines_uniform__get_icell_and_offset( self, x, icell, x_offset )

    ! 2. Compute index range of B-splines with support over cell 'icell'
    jmin = icell - self%offset

    ! 3. Compute derivatives of aforementioned B-splines
    !    Derivatives are normalized, hence they should be divided by dx
    associate( bspl => derivs, spline_degree => self%degree, start => self%inv_dx )

      ! only need splines of lower degree to compute derivatives
      bspl(0) = start
      do j = 1, spline_degree - 1
         j_real = real(j,wp)
         xx     = -x_offset
         saved  = 0.0_wp
         do r = 0, j-1
            xx      = xx + 1.0_wp
            temp    = bspl(r) * inv(j)
            bspl(r) = saved + xx * temp
            saved   = (j_real - xx) * temp
         end do
         bspl(j) = saved
      end do

      ! compute derivatives
      bjm1 = bspl(0)
      bj = bjm1
      bspl(0) = -bjm1
      do j=1, spline_degree - 1
         bj = bspl(j)
         bspl(j) = bjm1 - bj
         bjm1 = bj
      end do
      bspl(spline_degree) = bj

    end associate  ! bspl, spline_degree

  end subroutine s_bsplines_uniform__eval_deriv

  !-----------------------------------------------------------------------------
  !> Evaluate value and n derivatives at x of all basis functions with support in local cell
  SLL_PURE subroutine s_bsplines_uniform__eval_basis_and_n_derivs( self, x, n, derivs, jmin )
    class(sll_t_bsplines_uniform), intent(in   ) :: self
    real(wp)                     , intent(in   ) :: x
    integer                      , intent(in   ) :: n
    real(wp)                     , intent(  out) :: derivs(0:,0:)
    integer                      , intent(  out) :: jmin

    integer  :: icell
    real(wp) :: x_offset

    integer  :: j, r
    real(wp) :: j_real
    real(wp) :: xx
    real(wp) :: temp
    real(wp) :: saved

    integer  :: k, s1, s2, rk, pk, j1, j2
    real(wp) :: d

    ! GFortran: to allocate on stack use -fstack-arrays
    real(wp) :: ndu (0:self%degree, 0:self%degree)
    real(wp) :: a   (0:1          , 0:self%degree)

    ! Check on inputs
    SLL_ASSERT( n >= 0 )
    SLL_ASSERT( size(derivs,1) == 1+n           )
    SLL_ASSERT( size(derivs,2) == 1+self%degree )

    ! 1. Compute cell index 'icell' and x_offset
    call s_bsplines_uniform__get_icell_and_offset( self, x, icell, x_offset )

    ! 2. Compute index range of B-splines with support over cell 'icell'
    jmin = icell - self%offset

    ! 3. Recursively evaluate B-splines (see "sll_s_uniform_bsplines_eval_basis")
    !    up to self%degree, and store them all in the upper-right triangle of ndu
    associate( spline_degree => self%degree, bspl => derivs )

      ndu(0,0) = 1.0_wp
      do j = 1, spline_degree
         j_real = real(j,wp)
         xx     = -x_offset
         saved  = 0.0_wp
         do r = 0, j-1
            xx       = xx + 1.0_wp
            temp     = ndu(r,j-1) * inv(j)
            ndu(r,j) = saved + xx * temp
            saved    = (j_real - xx) * temp
         end do
         ndu(j,j) = saved
      end do
      bspl(0,:) = ndu(:,spline_degree)

    end associate  ! spline_degree, bspl

    ! 4. Use equation 2.10 in "The NURBS Book" to compute n derivatives
    associate( deg => self%degree, bsdx => derivs )

      do r = 0, deg
         s1 = 0
         s2 = 1
         a(0,0) = 1.0_wp
         do k = 1, n
            d  = 0.0_wp
            rk = r-k
            pk = deg-k
            if (r >= k) then
               a(s2,0) = a(s1,0) * inv(pk+1)
               d = a(s2,0) * ndu(rk,pk)
            end if
            if (rk > -1) then
               j1 = 1
            else
               j1 = -rk
            end if
            if (r-1 <= pk) then
               j2 = k-1
            else
               j2 = deg-r
            end if
            do j = j1, j2
               a(s2,j) = (a(s1,j) - a(s1,j-1)) * inv(pk+1)
               d = d + a(s2,j) * ndu(rk+j,pk)
            end do
            if (r <= pk) then
               a(s2,k) = - a(s1,k-1) * inv(pk+1)
               d = d + a(s2,k) * ndu(r,pk)
            end if
            bsdx(k,r) = d
            j  = s1
            s1 = s2
            s2 = j
         end do
      end do

      ! Multiply result by correct factors:
      ! deg!/(deg-n)! = deg*(deg-1)*...*(deg-n+1)
      ! k-th derivatives are normalized, hence they should be divided by dx^k
      d = real(deg,wp) * self%inv_dx
      do k = 1, n
         bsdx(k,:) = bsdx(k,:) * d
         d = d * real(deg-k,wp) * self%inv_dx
      end do

    end associate  ! deg, bsdx

  end subroutine s_bsplines_uniform__eval_basis_and_n_derivs

end module sll_m_bsplines_uniform
