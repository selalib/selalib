module sll_m_pic_base

#include "sll_working_precision.h"

  implicit none
  private

  !============================================================================
  ! Particle species
  !============================================================================
  type, public :: sll_species

    character(len=64) :: name !< species name
    sll_real64        :: q    !< charge of a single particle 
    sll_real64        :: m    !< mass   of a single particle

  contains
    procedure         :: q_over_m  !< charge over mass ratio

  end type sll_species


  !============================================================================
  ! Particle group
  !============================================================================
  type, public, abstract :: sll_particle_group_base

    class( sll_species ), pointer :: species
    sll_int32                     :: id
    
  contains
    ! Getters
    procedure( i_get_coords ), deferred :: get_x
    procedure( i_get_coords ), deferred :: get_v
    procedure( i_get_scalar ), deferred :: get_charge
    procedure( i_get_scalar ), deferred :: get_mass

    ! Setters
    procedure( i_set_coords ), deferred :: set_x
    procedure( i_set_coords ), deferred :: set_v

  end type sll_particle_group_base

  !----------------------------------------------------------------------------
  abstract interface
   pure function i_get_scalar( self, i ) result( r )
    use sll_working_precision
    import sll_particle_group_base
    class( sll_particle_group_base ), intent( in ) :: self
    sll_int32                       , intent( in ) :: i
    sll_real64 :: r
   end function i_get_scalar
  end interface

  !----------------------------------------------------------------------------
  abstract interface
   pure function i_get_coords( self, i ) result( r )
    use sll_working_precision
    import sll_particle_group_base
    class( sll_particle_group_base ), intent( in ) :: self
    sll_int32                       , intent( in ) :: i
    sll_real64 :: r(3)
   end function i_get_coords
  end interface

  !----------------------------------------------------------------------------
  abstract interface
   subroutine i_set_coords( self, i, x )
    use sll_working_precision
    import sll_particle_group_base
    class( sll_particle_group_base ), intent( inout ) :: self
    sll_int32                       , intent( in    ) :: i
    sll_real64                      , intent( in    ) :: x(3)
   end subroutine i_set_coords
  end interface

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
contains
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !----------------------------------------------------------------------------
  pure function q_over_m( self ) result( r )
    class( sll_species ), intent( in ) :: self
    sll_real64 :: r

    r = self%q / self%m
  end function q_over_m

end module sll_m_pic_base
