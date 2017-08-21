!> @ingroup splines
!> @brief   Access point to 1D splines (base class and factory)
!> @author  Yaman Güçlü - IPP Garching

module sll_m_spline_1d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"

  use sll_m_working_precision    , only: f64
  use sll_m_spline_1d_base       , only: sll_c_spline_1d
  use sll_m_spline_1d_uniform    , only: sll_t_spline_1d_uniform
  use sll_m_spline_1d_non_uniform, only: sll_t_spline_1d_non_uniform

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_periodic, &
    sll_p_hermite , &
    sll_p_greville

  implicit none

  public :: &
    sll_c_spline_1d, &
    sll_s_spline_1d_compute_num_cells, &
    sll_s_spline_1d_new

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Working precision
  integer, parameter :: wp = f64

  !> Allowed boundary conditions
  integer, parameter :: &
    allowed_bcs(*) = [sll_p_periodic, sll_p_hermite, sll_p_greville]

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !-----------------------------------------------------------------------------
  !> @brief     Calculate number of cells from number of interpolation points
  !> @details   Important for parallelization: for given spline degree and BCs,
  !>            calculate the number of grid cells that yields the desired
  !>            number of interpolation points
  !>
  !> @param[in]  degree   spline degree
  !> @param[in]  bc_xmin  boundary condition type at left  boundary (x=xmin)
  !> @param[in]  bc_xmax  boundary condition type at right boundary (x=xmax)
  !> @param[in]  nipts    desired number of interpolation points
  !> @param[out] ncells   number of cells in domain
  !-----------------------------------------------------------------------------
  subroutine sll_s_spline_1d_compute_num_cells( &
      degree , &
      bc_xmin, &
      bc_xmax, &
      nipts  , &
      ncells )

    integer, intent(in   ) :: degree
    integer, intent(in   ) :: bc_xmin
    integer, intent(in   ) :: bc_xmax
    integer, intent(in   ) :: nipts
    integer, intent(  out) :: ncells

    ! Sanity checks
    SLL_ASSERT( degree > 0  )
    SLL_ASSERT( any( bc_xmin == allowed_bcs ) )
    SLL_ASSERT( any( bc_xmax == allowed_bcs ) )

    if (any( [bc_xmin,bc_xmax]==sll_p_periodic ) .and. bc_xmin /= bc_xmax) then
      SLL_ERROR( "sll_s_spline_1d_compute_num_cells", "Incompatible BCs" )
    end if

    if (bc_xmin == sll_p_periodic) then
      ncells = nipts
    else
      associate( nbc_xmin => merge( degree/2, 0, bc_xmin == sll_p_hermite ), &
                 nbc_xmax => merge( degree/2, 0, bc_xmax == sll_p_hermite ) )
        ncells = nipts + nbc_xmin + nbc_xmax - degree
      end associate
    end if

  end subroutine sll_s_spline_1d_compute_num_cells

  !-----------------------------------------------------------------------------
  !> @brief      Allocate and initialize uniform or non-uniform 1D spline
  !> @param[out] spline   allocatable polymorphic object
  !> @param[in]  degree   spline degree
  !> @param[in]  ncells   number of cells in domain (one polynomial per cell)
  !> @param[in]  xmin     x coordinate of left  boundary of domain
  !> @param[in]  xmax     x coordinate of right boundary of domain
  !> @param[in]  bc_xmin  boundary condition type at x=xmin
  !> @param[in]  bc_xmax  boundary condition type at x=xmax
  !> @param[in]  breaks   list of break points (only for non-uniform spline)
  !-----------------------------------------------------------------------------
  subroutine sll_s_spline_1d_new( &
      spline , &
      degree , &
      ncells , &
      xmin   , &
      xmax   , &
      bc_xmin, &
      bc_xmax, &
      breaks )

    class(sll_c_spline_1d), allocatable, intent(  out) :: spline
    integer                            , intent(in   ) :: degree
    integer                            , intent(in   ) :: ncells
    real(wp)                           , intent(in   ) :: xmin
    real(wp)                           , intent(in   ) :: xmax
    integer                            , intent(in   ) :: bc_xmin
    integer                            , intent(in   ) :: bc_xmax
    real(wp),                  optional, intent(in   ) :: breaks(:)

    ! Sanity checks
    SLL_ASSERT( degree > 0  )
    SLL_ASSERT( ncells > 0  )
    SLL_ASSERT( xmin < xmax )
    SLL_ASSERT( any( bc_xmin == allowed_bcs ) )
    SLL_ASSERT( any( bc_xmax == allowed_bcs ) )

    ! Allocate 1D spline to correct type
    if (present( breaks ) .and. size( breaks ) > 0) then

      SLL_ASSERT( size( breaks ) == ncells+1     )
      SLL_ASSERT( xmin == breaks(1)              )
      SLL_ASSERT( xmax == breaks(size( breaks )) )

      allocate( sll_t_spline_1d_non_uniform :: spline )

    else

      allocate( sll_t_spline_1d_uniform     :: spline )

    end if

    ! Properly initialize 1D spline
    select type( spline )

    type is( sll_t_spline_1d_uniform )

      call spline % init( degree, ncells, xmin, xmax, bc_xmin, bc_xmax )

    type is( sll_t_spline_1d_non_uniform )

      call spline % init( degree, breaks, bc_xmin, bc_xmax )

    end select

  end subroutine sll_s_spline_1d_new

end module sll_m_spline_1d
