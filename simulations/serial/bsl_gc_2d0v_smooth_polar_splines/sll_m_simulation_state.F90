module sll_m_simulation_state
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"

  use sll_m_working_precision, only: f64

  use sll_m_spline_2d, only: sll_t_spline_2d

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

    type(sll_t_spline_2d)      :: spline_2d_rho
    type(sll_t_spline_2d)      :: spline_2d_phi
    type(sll_t_electric_field) :: electric_field

    integer :: nc
    logical :: point_charges_present
    type(sll_t_point_charge), allocatable :: point_charges(:)

  contains

    procedure :: init => s_simulation_state__init
    procedure :: free => s_simulation_state__free

  end type sll_t_simulation_state

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine s_simulation_state__init( self, ntau1, ntau2, nc, intensity, location )
    class(sll_t_simulation_state)     , intent(inout) :: self
    integer                           , intent(in   ) :: ntau1
    integer                           , intent(in   ) :: ntau2
    integer                           , intent(in   ) :: nc
    real(wp)                          , intent(in   ) :: intensity(:)
    real(wp)                          , intent(in   ) :: location(:,:)

    integer :: ic

    allocate( self % rho( ntau1, ntau2+1 ) )

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

  end subroutine s_simulation_state__init

  !-----------------------------------------------------------------------------
  subroutine s_simulation_state__free( self )
    class(sll_t_simulation_state), intent(inout) :: self

    deallocate( self % rho )

    if ( self % point_charges_present ) deallocate( self % point_charges )

  end subroutine s_simulation_state__free

end module sll_m_simulation_state
