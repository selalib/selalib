module sll_m_simulation_state
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"

  use sll_m_working_precision, only: f64

  use sll_m_bsplines_base, only: sll_c_bsplines

  use sll_m_spline_2d, only: sll_t_spline_2d

  use sll_m_singular_mapping_discrete, only: sll_t_singular_mapping_discrete

  use sll_m_point_charge, only: sll_t_point_charge

  use sll_m_electric_field, only: sll_t_electric_field

  implicit none

  public :: sll_t_simulation_state

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Working precision
  integer, parameter :: wp = f64

  type :: sll_t_simulation_state

    real(wp), allocatable :: rho(:,:)

    class(sll_c_bsplines), pointer :: bsplines_eta1
    class(sll_c_bsplines), pointer :: bsplines_eta2

    type(sll_t_singular_mapping_discrete), pointer :: mapping_discrete

    type(sll_t_spline_2d)      :: spline_2d_rho
    type(sll_t_spline_2d)      :: spline_2d_phi
    type(sll_t_electric_field) :: electric_field

    integer :: nc
    logical :: point_charges_present
    type(sll_t_point_charge), allocatable :: point_charges(:)

    logical :: evolve_background

  contains

    procedure :: init => s_simulation_state__init
    procedure :: copy => s_simulation_state__copy
    procedure :: free => s_simulation_state__free

  end type sll_t_simulation_state

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine s_simulation_state__init( &
    self            , &
    bsplines_eta1   , &
    bsplines_eta2   , &
    mapping_discrete, &
    ntau1           , &
    ntau2           , &
    nc              , &
    intensity       , &
    location        , &
    evolve_background )
    class(sll_t_simulation_state)                 , intent(inout) :: self
    class(sll_c_bsplines)                , pointer, intent(in   ) :: bsplines_eta1
    class(sll_c_bsplines)                , pointer, intent(in   ) :: bsplines_eta2
    type(sll_t_singular_mapping_discrete), pointer, intent(in   ) :: mapping_discrete
    integer                                       , intent(in   ) :: ntau1
    integer                                       , intent(in   ) :: ntau2
    integer                                       , intent(in   ) :: nc
    real(wp)                                      , intent(in   ) :: intensity(:)
    real(wp)                                      , intent(in   ) :: location(:,:)
    logical                                       , intent(in   ) :: evolve_background

    integer :: ic

    self % bsplines_eta1 => bsplines_eta1
    self % bsplines_eta2 => bsplines_eta2

    self % mapping_discrete => mapping_discrete

    allocate( self % rho( ntau1, ntau2+1 ) )

    call self % spline_2d_rho  % init( self % bsplines_eta1, self % bsplines_eta2 )
    call self % spline_2d_phi  % init( self % bsplines_eta1, self % bsplines_eta2 )
    call self % electric_field % init( self % mapping_discrete, self % spline_2d_phi )

    SLL_ASSERT( nc == 0 .or. nc > 0 )

    self % nc = nc

    if ( nc == 0 ) then
      self % point_charges_present = .false.
    else
      self % point_charges_present = .true.
    end if

    if ( self % point_charges_present ) then

      SLL_ASSERT( nc == size( intensity ) )
      SLL_ASSERT( nc == size( location, 2 ) )
      SLL_ASSERT( 2  == size( location, 1 ) )

      allocate( self % point_charges( nc ) )
      do ic = 1, nc
        call self % point_charges(ic) % init( intensity(ic), location(:,ic) )
      end do

    end if

    self % evolve_background = evolve_background

  end subroutine s_simulation_state__init

  !-----------------------------------------------------------------------------
  subroutine s_simulation_state__copy( self, sim_state_copy )
    class(sll_t_simulation_state), intent(in   ) :: self
    class(sll_t_simulation_state), intent(  out) :: sim_state_copy

    integer :: ic

    sim_state_copy % bsplines_eta1 => self % bsplines_eta1
    sim_state_copy % bsplines_eta2 => self % bsplines_eta2

    sim_state_copy % mapping_discrete => self % mapping_discrete

    allocate( sim_state_copy % rho( size( self % rho, 1 ), size( self % rho, 2 ) ) )

    call sim_state_copy % spline_2d_rho  % init( sim_state_copy % bsplines_eta1, sim_state_copy % bsplines_eta2 )
    call sim_state_copy % spline_2d_phi  % init( sim_state_copy % bsplines_eta1, sim_state_copy % bsplines_eta2 )
    call sim_state_copy % electric_field % init( sim_state_copy % mapping_discrete, sim_state_copy % spline_2d_phi )

    sim_state_copy % nc = self % nc

    if ( sim_state_copy % nc == 0 ) then
      sim_state_copy % point_charges_present = .false.
    else
      sim_state_copy % point_charges_present = .true.
    end if

    if ( sim_state_copy % point_charges_present ) then

      allocate( sim_state_copy % point_charges( sim_state_copy % nc ) )
      do ic = 1, sim_state_copy % nc
        sim_state_copy % point_charges(ic) % intensity = self % point_charges(ic) % intensity
        sim_state_copy % point_charges(ic) % location  = self % point_charges(ic) % location
      end do

    end if

    sim_state_copy % evolve_background = self % evolve_background

  end subroutine s_simulation_state__copy

  !-----------------------------------------------------------------------------
  subroutine s_simulation_state__free( self )
    class(sll_t_simulation_state), intent(inout) :: self

    deallocate( self % rho )

    if ( self % point_charges_present ) deallocate( self % point_charges )

  end subroutine s_simulation_state__free

end module sll_m_simulation_state
