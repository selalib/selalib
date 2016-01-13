module sll_m_particle_group_base

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  implicit none

  public :: &
    sll_c_particle_group_base, &
    sll_t_species

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !============================================================================
  ! Particle species
  !============================================================================
  type :: sll_t_species

    character(len=64) :: name !< species name
    sll_real64        :: q    !< charge of a single particle
    sll_real64        :: m    !< mass   of a single particle

  contains
    procedure         :: q_over_m  !< charge over mass ratio
    procedure         :: initialize => initialize_species !< Constructor

  end type sll_t_species


  !============================================================================
  ! Particle group
  !============================================================================
  type, abstract :: sll_c_particle_group_base

    class( sll_t_species ), allocatable :: species
    sll_int32                     :: id
    sll_int32                     :: n_particles !< number of particles local to the processor
    sll_int32                     :: n_total_particles !< number of particles in total simulation    
    sll_int32                     :: n_weights !< number of weights per particle

  contains
    ! Getters
    procedure( i_get_coords ), deferred :: get_x
    procedure( i_get_coords ), deferred :: get_v
    procedure( i_get_scalar ), deferred :: get_charge
    procedure( i_get_scalar ), deferred :: get_mass
    procedure( i_get_array  ), deferred :: get_weights

    ! Setters
    procedure( i_set_coords ), deferred :: set_x
    procedure( i_set_coords ), deferred :: set_v
    procedure( i_set_array  ), deferred :: set_weights
    procedure( set_scalar ),   deferred :: set_common_weight

    procedure( empty ),        deferred :: delete

!    ! Getters for the whole group
!    procedure( get_all_coords), deferred :: get_all_x
!    procedure( get_all_coords), deferred :: get_all_v

  end type sll_c_particle_group_base

  !----------------------------------------------------------------------------
  abstract interface
   pure function i_get_scalar( self, i , i_weight) result( r )
    use sll_m_working_precision
    import sll_c_particle_group_base
    class( sll_c_particle_group_base ), intent( in ) :: self
    sll_int32                       , intent( in ) :: i
    sll_int32, optional             , intent( in ) :: i_weight
    sll_real64 :: r
   end function i_get_scalar
  end interface

  !----------------------------------------------------------------------------
  abstract interface
   pure function i_get_coords( self, i ) result( r )
    use sll_m_working_precision
    import sll_c_particle_group_base
    class( sll_c_particle_group_base ), intent( in ) :: self
    sll_int32                       , intent( in ) :: i
    sll_real64 :: r(3)
   end function i_get_coords
  end interface

  !----------------------------------------------------------------------------
  abstract interface
   pure function i_get_array( self, i ) result( r )
    use sll_m_working_precision
    import sll_c_particle_group_base
    class( sll_c_particle_group_base ), intent( in ) :: self
    sll_int32                       , intent( in ) :: i
    sll_real64 :: r(self%n_weights)
  end function i_get_array
  end interface

  !----------------------------------------------------------------------------
  abstract interface
   subroutine i_set_coords( self, i, x )
    use sll_m_working_precision
    import sll_c_particle_group_base
    class( sll_c_particle_group_base ), intent( inout ) :: self
    sll_int32                       , intent( in    ) :: i
    sll_real64                      , intent( in    ) :: x(3)
   end subroutine i_set_coords
  end interface

!----------------------------------------------------------------------------
  abstract interface
   subroutine i_set_scalar( self, i, x )
    use sll_m_working_precision
    import sll_c_particle_group_base
    class( sll_c_particle_group_base ), intent( inout ) :: self
    sll_int32                       , intent( in    ) :: i
    sll_real64                      , intent( in    ) :: x
  end subroutine i_set_scalar
  end interface

!----------------------------------------------------------------------------
  abstract interface
   subroutine set_scalar( self, x )
    use sll_m_working_precision
    import sll_c_particle_group_base
    class( sll_c_particle_group_base ), intent( inout ) :: self
    sll_real64                      , intent( in    ) :: x
  end subroutine set_scalar
  end interface

  !----------------------------------------------------------------------------
  abstract interface
   subroutine i_set_array( self, i, x )
    use sll_m_working_precision
    import sll_c_particle_group_base
    class( sll_c_particle_group_base ), intent( inout ) :: self
    sll_int32                       , intent( in    ) :: i
    sll_real64                      , intent( in    ) :: x(self%n_weights)
  end subroutine i_set_array
  end interface

 !---------------------------------------------------------------------------!  
  abstract interface
     subroutine empty(self)
       import sll_c_particle_group_base
       class( sll_c_particle_group_base ), intent( inout ) :: self

     end subroutine empty
  end interface

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
contains
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  !----------------------------------------------------------------------------
  pure function q_over_m( self ) result( r )
    class( sll_t_species ), intent( in ) :: self
    sll_real64 :: r

    r = self%q / self%m
  end function q_over_m

  
  !----------------------------------------------------------------------------
  subroutine initialize_species(&
       self, &
       species_charge,     &
       species_mass        &
       )
    class(sll_t_species), intent ( out ) :: self
    sll_real64, intent ( in )   :: species_charge
    sll_real64, intent ( in )   :: species_mass

    sll_int32  :: ierr

    self%q = species_charge
    self%m = species_mass

  end subroutine initialize_species


end module sll_m_particle_group_base
