!> @ingroup splines
!> @brief   Access point to B-splines (base class and factory)
!> @author  Yaman Güçlü - IPP Garching

module sll_m_bsplines
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"

  use sll_m_working_precision   , only: f64
  use sll_m_bsplines_base       , only: sll_c_bsplines
  use sll_m_bsplines_uniform    , only: sll_t_bsplines_uniform
  use sll_m_bsplines_non_uniform, only: sll_t_bsplines_non_uniform

  implicit none

  public :: &
    sll_c_bsplines, &
    sll_s_bsplines_new

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Working precision
  integer, parameter :: wp = f64

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !-----------------------------------------------------------------------------
  !> @brief      Allocate and initialize uniform or non-uniform B-splines
  !> @param[out] bsplines  allocatable polymorphic object
  !> @param[in]  degree    spline degree
  !> @param[in]  periodic  .true. if domain is periodic
  !> @param[in]  xmin      x coordinate of left  boundary of domain
  !> @param[in]  xmax      x coordinate of right boundary of domain
  !> @param[in]  ncells    number of cells in domain (one polynomial per cell)
  !> @param[in]  breaks    list of breakpoints (only for non-uniform B-splines)
  !-----------------------------------------------------------------------------
  subroutine sll_s_bsplines_new( &
      bsplines, &
      degree  , &
      periodic, &
      xmin    , &
      xmax    , &
      ncells  , &
      breaks  )

    class(sll_c_bsplines), allocatable, intent(  out) :: bsplines
    integer                           , intent(in   ) :: degree
    logical                           , intent(in   ) :: periodic
    real(wp)                          , intent(in   ) :: xmin
    real(wp)                          , intent(in   ) :: xmax
    integer                           , intent(in   ) :: ncells
    real(wp),                 optional, intent(in   ) :: breaks(:)

    logical :: uniform

    ! Sanity checks
    SLL_ASSERT( degree > 0  )
    SLL_ASSERT( ncells > 0  )
    SLL_ASSERT( xmin < xmax )

    ! Determine if B-splines are uniform based on 'breaks' optional argument
    if (present( breaks )) then
      uniform = (size( breaks ) == 0) ! Fall back to uniform if array is empty
    else
      uniform = .true.
    end if

    ! Non-uniform case: perform sanity checks on breakpoints
    if (.not. uniform) then
      SLL_ASSERT( size( breaks ) == ncells+1     )
      SLL_ASSERT( xmin == breaks(1)              )
      SLL_ASSERT( xmax == breaks(size( breaks )) )
    end if

    ! Allocate B-splines object to correct type
    if (uniform) then
      allocate( sll_t_bsplines_uniform     :: bsplines )
    else
      allocate( sll_t_bsplines_non_uniform :: bsplines )
    end if

    ! Initialize B-splines object with type-specific constructor
    select type( bsplines )
    type is( sll_t_bsplines_uniform )
      call bsplines % init( degree, periodic, xmin, xmax, ncells )
    type is( sll_t_bsplines_non_uniform )
      call bsplines % init( degree, periodic, breaks )
    end select

  end subroutine sll_s_bsplines_new

end module sll_m_bsplines
