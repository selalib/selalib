!> @ingroup splines
!> @brief   Access point to 1D splines (base class and factory)
!> @author  Yaman Güçlü  - IPP Garching
!> @author  Edoardo Zoni - IPP Garching

module sll_m_spline_2d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"

  use sll_m_working_precision, only: &
    f64

  use sll_m_spline_1d, only: &
    sll_s_spline_1d_compute_num_cells

  use sll_m_spline_2d_base, only: &
    sll_c_spline_2d, &
    sll_t_spline_2d_boundary_data

  use sll_m_spline_2d_uniform, only: &
    sll_t_spline_2d_uniform

  use sll_m_spline_2d_non_uniform, only: &
    sll_t_spline_2d_non_uniform

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_periodic, &
    sll_p_hermite , &
    sll_p_greville

  implicit none

  public :: &
    sll_c_spline_2d, &
    sll_s_spline_2d_compute_num_cells, &
    sll_s_spline_2d_new, &
    sll_t_spline_2d_boundary_data

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
  subroutine sll_s_spline_2d_compute_num_cells( &
      degree , &
      bc_xmin, &
      bc_xmax, &
      nipts  , &
      ncells )

    integer, intent(in   ) :: degree (2)
    integer, intent(in   ) :: bc_xmin(2)
    integer, intent(in   ) :: bc_xmax(2)
    integer, intent(in   ) :: nipts  (2)
    integer, intent(  out) :: ncells (2)

    integer :: i

    do i = 1, 2
      call sll_s_spline_1d_compute_num_cells( &
        degree(i), bc_xmin(i), bc_xmax(i), nipts(i), ncells(i) )
    end do

  end subroutine sll_s_spline_2d_compute_num_cells

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
  subroutine sll_s_spline_2d_new( &
      spline , &
      degree , &
      ncells , &
      xmin   , &
      xmax   , &
      bc_xmin, &
      bc_xmax, &
      breaks1, &
      breaks2 )

    class(sll_c_spline_2d), allocatable, intent(  out) :: spline
    integer                            , intent(in   ) :: degree (2)
    integer                            , intent(in   ) :: ncells (2)
    real(wp)                           , intent(in   ) :: xmin   (2)
    real(wp)                           , intent(in   ) :: xmax   (2)
    integer                            , intent(in   ) :: bc_xmin(2)
    integer                            , intent(in   ) :: bc_xmax(2)
    real(wp),                  optional, intent(in   ) :: breaks1(:)
    real(wp),                  optional, intent(in   ) :: breaks2(:)

    ! Sanity checks
    SLL_ASSERT( all( degree > 0 )  )
    SLL_ASSERT( all( ncells > 0 )  )
    SLL_ASSERT( all( xmin < xmax ) )
    SLL_ASSERT( any( bc_xmin(1) == allowed_bcs ) )
    SLL_ASSERT( any( bc_xmax(1) == allowed_bcs ) )
    SLL_ASSERT( any( bc_xmin(2) == allowed_bcs ) )
    SLL_ASSERT( any( bc_xmax(2) == allowed_bcs ) )

    ! Allocate 2D spline to correct type
    if (present( breaks1 ) .and. size( breaks1 ) > 0 .and. &
        present( breaks2 ) .and. size( breaks2 ) > 0) then

      SLL_ASSERT( size( breaks1 ) == ncells(1)+1      )
      SLL_ASSERT( xmin(1) == breaks1(1)               )
      SLL_ASSERT( xmax(1) == breaks1(size( breaks1 )) )

      SLL_ASSERT( size( breaks2 ) == ncells(2)+1      )
      SLL_ASSERT( xmin(2) == breaks2(1)               )
      SLL_ASSERT( xmax(2) == breaks2(size( breaks2 )) )

      allocate( sll_t_spline_2d_non_uniform :: spline )

    else

      allocate( sll_t_spline_2d_uniform     :: spline )

    end if

    ! Properly initialize 1D spline
    select type( spline )
    type is( sll_t_spline_2d_uniform )
      call spline % init( degree, ncells, xmin, xmax, bc_xmin, bc_xmax )
    type is( sll_t_spline_2d_non_uniform )
      call spline % init( degree(:), breaks1, breaks2, bc_xmin(:), bc_xmax(:) )
    end select

  end subroutine sll_s_spline_2d_new

end module sll_m_spline_2d
