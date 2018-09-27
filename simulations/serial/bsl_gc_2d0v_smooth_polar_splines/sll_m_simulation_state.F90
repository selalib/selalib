module sll_m_simulation_state
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"

  use sll_m_working_precision, only: f64

  use sll_m_spline_2d, only: sll_t_spline_2d

  use sll_m_point_charge, only: sll_t_point_charge

  implicit none

  public :: sll_t_simulation_state

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Working precision
  integer, parameter :: wp = f64

  type :: sll_t_simulation_state

    real(wp), allocatable :: rho(:,:)

    type(sll_t_spline_2d), pointer :: spline_2d_rho => null()
    type(sll_t_spline_2d), pointer :: spline_2d_phi => null()

    type(sll_t_point_charge), allocatable :: point_charges(:)

  contains

    procedure :: init => s_simulation_state__init
    procedure :: copy => s_simulation_state__copy
    procedure :: free => s_simulation_state__free

  end type sll_t_simulation_state

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine s_simulation_state__init( &
    self         , &
    ntau1        , &
    ntau2        , &
    nc           , &
    intensity    , &
    location     , &
    spline_2d_rho, &
    spline_2d_phi )
    class(sll_t_simulation_state)   , intent(inout) :: self
    integer                         , intent(in   ) :: ntau1
    integer                         , intent(in   ) :: ntau2
    integer                         , intent(in   ) :: nc
    real(wp)                        , intent(in   ) :: intensity(:)
    real(wp)                        , intent(in   ) :: location(:,:)
    type(sll_t_spline_2d), target   , intent(in   ) :: spline_2d_rho
    type(sll_t_spline_2d), target   , intent(in   ) :: spline_2d_phi

    integer :: ic

    allocate( self % rho( ntau1, ntau2+1 ) )

    self % spline_2d_rho     => spline_2d_rho
    self % spline_2d_phi     => spline_2d_phi

    SLL_ASSERT( nc == 0 .or. nc == size( intensity ) )
    SLL_ASSERT( nc == 0 .or. nc == size( location, 2 ) )
    SLL_ASSERT( 2  == size( location, 1 ) )

    if ( nc /= 0 ) then
      allocate( self % point_charges( nc ) )
      do ic = 1, nc
        call self % point_charges(ic) % init( intensity(ic), location(:,ic) )
      end do
    end if

  end subroutine s_simulation_state__init

  !-----------------------------------------------------------------------------
  subroutine s_simulation_state__copy( self, sim_state_copy )
    class(sll_t_simulation_state), intent(in   ) :: self
    class(sll_t_simulation_state), intent(inout) :: sim_state_copy

    associate( rho => self % rho, point_charges => self % point_charges )
      allocate( sim_state_copy % rho( size( rho, 1 ), size( rho, 2 ) ) , source = rho )
      allocate( sim_state_copy % point_charges( size( point_charges ) ), source = point_charges )
    end associate

  end subroutine s_simulation_state__copy

  !-----------------------------------------------------------------------------
  subroutine s_simulation_state__free( self )
    class(sll_t_simulation_state), intent(inout) :: self

    deallocate( self % rho )

    nullify( self % spline_2d_rho )
    nullify( self % spline_2d_phi )

    deallocate( self % point_charges )

  end subroutine s_simulation_state__free

end module sll_m_simulation_state
