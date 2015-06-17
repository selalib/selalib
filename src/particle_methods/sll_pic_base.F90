module sll_module_pic_base

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
    procedure           :: q_over_m  !< charge over mass ratio
    procedure           :: species_new

  end type sll_species


  !============================================================================
  ! Particle group
  !============================================================================
  type, public, abstract :: sll_particle_group_base

    class( sll_species ), pointer :: species
    sll_int32                     :: id

    sll_int32                     :: dimension_x
    sll_int32                     :: dimension_v
    
  contains
    ! Getters
    procedure( i_get_coords ), deferred :: get_x
    procedure( i_get_coords ), deferred :: get_v
    procedure( i_get_scalar ), deferred :: get_charge
    procedure( i_get_scalar ), deferred :: get_mass

    procedure( i_get_integer), deferred :: get_cell_index

    ! Setters
    procedure( i_set_coords ), deferred :: set_x
    procedure( i_set_coords ), deferred :: set_v
    procedure( i_set_scalar ), deferred :: set_particle_weight
    procedure( set_scalar   ), deferred :: set_common_weight

    ! Charge deposition
    procedure( dep_charge_2d), deferred :: deposit_charge_2d

    ! Initializers
    procedure( set_landau_params    ), deferred :: set_landau_parameters
    procedure( rand_init            ),  deferred :: random_initializer
    procedure( cart_init            ),  deferred :: cartesian_initializer

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
   pure function :: i_get_integer(self, i ) result( i_out )
    use sll_working_precision
    import sll_particle_group_base
    class( sll_particle_group_base ),   intent( in )    ::  self
    sll_int32,                          intent( in )    ::  i
    sll_int32  ::  i_out

   end function i_get_integer
  end interface

  !----------------------------------------------------------------------------
  abstract interface
   subroutine i_set_scalar( self, i, s )
    use sll_working_precision
    import sll_particle_group_base
    class( sll_particle_group_base ), intent( inout ) :: self
    sll_int32                       , intent( in    ) :: i
    sll_real64                      , intent( in    ) :: s
   end subroutine i_set_scalar
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

  !----------------------------------------------------------------------------
  abstract interface
   subroutine rand_init( self, nb_particles, initial_density_identifier, rand_seed, rank, world_size)
    use sll_working_precision
    import sll_particle_group_base
    class( sll_particle_group_base ), intent( inout ) :: self
    sll_int32                       , intent( in    ) :: nb_particles
    sll_int32                       , intent( in    ) :: initial_density_identifier
    sll_int32, dimension(:)         , intent( in ), optional :: rand_seed
    sll_int32                       , intent( in ), optional :: rank, world_size
   end subroutine rand_init
  end interface

  !----------------------------------------------------------------------------
  abstract interface
   subroutine cart_init( self, nb_cells_x, nb_cells_v, x_min, x_max, v_min, v_max )
    use sll_working_precision
    import sll_particle_group_base
    class( sll_particle_group_base ), intent( inout ) :: self
    sll_int32                       , intent( in    ) :: nb_cells_x(3)
    sll_int32                       , intent( in    ) :: nb_cells_v(3)
    sll_real64                      , intent( in    ) :: x_min(3)
    sll_real64                      , intent( in    ) :: x_max(3)
    sll_real64                      , intent( in    ) :: v_min(3)
    sll_real64                      , intent( in    ) :: v_max(3)
   end subroutine cart_init
  end interface

  !----------------------------------------------------------------------------
  abstract interface
   subroutine set_landau_params( self, thermal_speed, alpha, k_landau )
    use sll_working_precision
    import sll_particle_group_base
    class( sll_particle_group_base ), intent( inout ) :: self
    sll_real64                      , intent( in    ) :: thermal_speed
    sll_real64                      , intent( in    ) :: alpha
    sll_real64                      , intent( in    ) :: k_landau

   end subroutine set_landau_params
  end interface


  !----------------------------------------------------------------------------
  abstract interface
   subroutine dep_charge_2d( self, charge_accumulator )
    use sll_working_precision
    import sll_particle_group_base
    class( sll_particle_group_base ),           intent( inout ) :: self
    type( sll_charge_accumulator_2d ), pointer, intent( inout ) :: charge_accumulator
   end subroutine dep_charge_2d
  end interface


!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
contains
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  pure function species_new( &
      species_charge,     &
      species_mass,       &
  ) result(res)

    type(sll_species), pointer :: res
    res%q = species_charge
    res%m = species_mass

  end function species_new


  !----------------------------------------------------------------------------
  pure function q_over_m( self ) result( r )
    class( sll_species ), intent( in ) :: self
    sll_real64 :: r

    r = self%q / self%m
  end function q_over_m

end module sll_module_pic_base
