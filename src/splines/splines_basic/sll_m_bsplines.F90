!> @ingroup splines
!> @brief   Access point to B-splines of arbitrary degree providing factory function
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
    sll_s_bsplines_new, &
    sll_s_bsplines_new_2d_polar, &
    sll_s_bsplines_new_mirror_copy

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

    class(sll_c_bsplines), allocatable, intent(inout) :: bsplines
    integer                           , intent(in   ) :: degree
    logical                           , intent(in   ) :: periodic
    real(wp)                          , intent(in   ) :: xmin
    real(wp)                          , intent(in   ) :: xmax
    integer                           , intent(in   ) :: ncells
    real(wp),                 optional, intent(in   ) :: breaks(:)

    logical :: uniform

    ! Sanity checks
    SLL_ASSERT( .not. allocated( bsplines ) )
    SLL_ASSERT( degree > 0  )
    SLL_ASSERT( ncells > 0  )
    SLL_ASSERT( xmin < xmax )

    ! Determine if B-splines are uniform based on 'breaks' optional argument
    uniform = .not. present( breaks )

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

    ! Set 'radial' flag to false
    bsplines % radial = .false.

  end subroutine sll_s_bsplines_new

  !-----------------------------------------------------------------------------
  !> @brief      Allocate and initialize two B-splines for use on polar grid
  !-----------------------------------------------------------------------------
  subroutine sll_s_bsplines_new_2d_polar( &
      bsplines_radial , &
      bsplines_angular, &
      degree          , &
      ncells          , &
      max_radius      , &
      theta_lims      , &
      breaks_radial   , &
      breaks_angular  )

    class(sll_c_bsplines), allocatable, intent(inout) :: bsplines_radial
    class(sll_c_bsplines), allocatable, intent(inout) :: bsplines_angular
    integer                           , intent(in   ) :: degree(2)
    integer                           , intent(in   ) :: ncells(2)
    real(wp)                          , intent(in   ) :: max_radius
    real(wp)                          , intent(in   ) :: theta_lims(2)
    real(wp),                 optional, intent(in   ) :: breaks_radial (:)
    real(wp),                 optional, intent(in   ) :: breaks_angular(:)

    character(len=*), parameter :: this_sub_name = "sll_s_bsplines_new_2d_polar"

    logical  :: uniform(2)
    real(wp) :: dr

    ! Sanity checks
    SLL_ASSERT( .not. allocated( bsplines_radial  ) )
    SLL_ASSERT( .not. allocated( bsplines_angular ) )
    SLL_ASSERT( degree(1) > 0 )
    SLL_ASSERT( degree(2) > 0 )
    SLL_ASSERT( ncells(1) > 0 )
    SLL_ASSERT( ncells(2) > 0 )
    SLL_ASSERT( theta_lims(1) < theta_lims(2) )

    ! Specific checks for polar geometry
    SLL_ASSERT( max_radius > 0.0_wp )
    SLL_ASSERT( modulo( ncells(2), 2 ) == 0 )

    ! We do not accept non-uniform splines yet
    if (present( breaks_radial )) then
      uniform(1) = .false.
      SLL_ERROR( this_sub_name, "non-uniform option not yet implemented." )
    else
      uniform(1) = .true.
    end if

    if (present( breaks_angular )) then
      uniform(2) = .false.
      SLL_ERROR( this_sub_name, "non-uniform option not yet implemented." )
    else
      uniform(2) = .true.
    end if

    dr = max_radius / (real(ncells(1),wp)-0.5_wp)

    call sll_s_bsplines_new( &
      bsplines = bsplines_radial, &
      degree   = degree(1) , &
      periodic = .false.   , &
      xmin     = -0.5_wp*dr, &
      xmax     = max_radius, &
      ncells   = ncells(1) , &
      breaks   = breaks_radial )

    call sll_s_bsplines_new( &
      bsplines = bsplines_angular, &
      degree   = degree(2)    , &
      periodic = .true.       , &
      xmin     = theta_lims(1), &
      xmax     = theta_lims(2), &
      ncells   = ncells(2)    , &
      breaks   = breaks_angular )

    ! Set 'radial' flag to true
    bsplines_radial % radial = .true.

  end subroutine sll_s_bsplines_new_2d_polar

  !-----------------------------------------------------------------------------
  !> @brief  Create B-splines by mirroring radial basis onto extended domain [-R,R]
  !-----------------------------------------------------------------------------
  subroutine sll_s_bsplines_new_mirror_copy( radial_basis, extended_basis )
    class(sll_c_bsplines),              intent(in   ) :: radial_basis
    class(sll_c_bsplines), allocatable, intent(inout) :: extended_basis

    character(len=*), parameter :: this_sub_name = "sll_s_bsplines_new_mirror_copy"

    SLL_ASSERT( radial_basis % radial )
    SLL_ASSERT( .not. allocated( extended_basis ) )

    select type( radial_basis )

    type is( sll_t_bsplines_uniform )
      call sll_s_bsplines_new( &
        bsplines =  extended_basis, &
        degree   =  radial_basis % degree, &
        periodic = .false., &
        xmin     = -radial_basis % xmax, &
        xmax     =  radial_basis % xmax, &
        ncells   =  radial_basis % ncells*2 - 1 )

    class default
      SLL_ERROR( this_sub_name, "non-uniform option not yet implemented." )

    end select

  end subroutine sll_s_bsplines_new_mirror_copy

end module sll_m_bsplines
