!> @ingroup particle_methods
!> @brief
!> Common interface for all groups of particles
module sll_module_pic_base

#include "sll_working_precision.h"
#include "sll_memory.h"

  implicit none
  private

  public :: species_new, apply_periodic_bc_on_cartesian_mesh_2d, x_is_in_domain_2d

  !============================================================================
  !> @brief Particle species
  !============================================================================
  type, public :: sll_species
    !      character(len=64) :: name !< species name
      sll_real64        :: q    !< charge of a single particle
      sll_real64        :: m    !< mass   of a single particle

  contains
    procedure, pass(self) :: q_over_m => get_q_over_m
!    procedure  :: species_new ! => create_species_new      !! is that a good syntax for the constructor interface?

  end type sll_species

!  interface sll_species
!    procedure, public species_new
!  end interface

  !============================================================================
  !> @brief Particle group
  !============================================================================
  type, public, abstract :: sll_particle_group_base

    class( sll_species ), pointer :: species
    sll_int32                     :: id

    sll_int32                     :: dimension_x
    sll_int32                     :: dimension_v

    logical                       :: domain_is_periodic(3)

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
    procedure( set_landau_params ),  deferred :: set_landau_parameters
    procedure( init              ),  deferred :: initializer

    ! Visualize
    procedure( vis ),  deferred :: visualize_f_slice_x_vx

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
   pure function i_get_integer(self, i ) result( i_out )
    use sll_working_precision
    import sll_particle_group_base
    class( sll_particle_group_base ),   intent( in )    ::  self
    sll_int32,                          intent( in )    ::  i
    sll_int32  ::  i_out

   end function i_get_integer
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
   subroutine set_scalar( self, s )
    use sll_working_precision
    import sll_particle_group_base
    class( sll_particle_group_base ), intent( inout ) :: self
    sll_real64                      , intent( in    ) :: s
   end subroutine set_scalar
  end interface

  !----------------------------------------------------------------------------
  abstract interface
   subroutine init( self, initial_density_identifier, rand_seed, rank, world_size)
    use sll_working_precision
    import sll_particle_group_base
    class( sll_particle_group_base ), intent( inout ) :: self
    sll_int32                       , intent( in    ) :: initial_density_identifier
    sll_int32, dimension(:)         , intent( in ), optional :: rand_seed
    sll_int32                       , intent( in ), optional :: rank, world_size
   end subroutine init
  end interface

  !----------------------------------------------------------------------------
  abstract interface
   subroutine vis( self, array_name, iplot )
    use sll_working_precision
    import sll_particle_group_base
    class( sll_particle_group_base ), intent( inout ) :: self
    character(len=*),                   intent(in)      :: array_name   !< field name
    sll_int32,                          intent(in)      :: iplot        !< plot counter

   end subroutine vis
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
   subroutine dep_charge_2d( self, charge_accumulator, target_total_charge )
    use sll_working_precision
#include "sll_accumulators.h"
    import sll_particle_group_base
    class( sll_particle_group_base ),           intent( inout ) :: self
    type( sll_charge_accumulator_2d ), pointer, intent( inout ) :: charge_accumulator
    sll_real64,                                 intent(in), optional :: target_total_charge
   end subroutine dep_charge_2d
  end interface


!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
contains
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  function species_new( &
      species_charge,     &
      species_mass        &
  ) result(res)
    sll_real64, intent ( in )   :: species_charge
    sll_real64, intent ( in )   :: species_mass
    type(sll_species), pointer  :: res
    sll_int32  :: ierr

    SLL_ALLOCATE( res, ierr )
    !    res%name =
    res%q = species_charge
    res%m = species_mass

  end function species_new


  !----------------------------------------------------------------------------
  function get_q_over_m( self ) result( r )
    class( sll_species ), intent( in ) :: self
    sll_real64 :: r

    r = self%q / self%m

  end function get_q_over_m



  ! todo: put this in the right module (with the meshes...)

  ! tells whether the given point is in the given domain, with boolean arguments for the domain periodicity
  ! (taken from previous function in_bounds_periodic)
  function x_is_in_domain_2d( x, y, mesh, x_periodic, y_periodic ) result(res)

    use sll_cartesian_meshes
    sll_real64,                     intent( in )            :: x, y
    type(sll_cartesian_mesh_2d),    intent( in ), pointer   :: mesh
    logical,                        intent( in )            :: x_periodic
    logical,                        intent( in )            :: y_periodic
    logical     :: res

    res = ( x >= mesh%eta1_min )                                                                                    &
          .and.                                                                                                     &
          ( ( x < mesh%eta1_max .and. x_periodic ) .or. ( x <= mesh%eta1_max .and. .not. x_periodic ) )             &
          .and.                                                                                                     &
          ( y >= mesh%eta2_min )                                                                                    &
          .and.                                                                                                     &
          ( ( y < mesh%eta2_max .and. y_periodic ) .or. ( y <= mesh%eta2_max .and. .not. y_periodic) )

  end function x_is_in_domain_2d

  ! <<apply_periodic_bc_on_cartesian_mesh_2d>>
  
  ! todo: put this in the right module (with the meshes...)

  subroutine apply_periodic_bc_on_cartesian_mesh_2d( mesh, x, y )

    use sll_cartesian_meshes
    ! [[file:../working_precision/sll_working_precision.h]]
    use sll_working_precision

    type(sll_cartesian_mesh_2d), pointer :: mesh
    sll_real64, intent(inout) :: x
    sll_real64, intent(inout) :: y

    x = mesh%eta1_min + modulo(x - mesh%eta1_min, mesh%eta1_max - mesh%eta1_min)
    y = mesh%eta2_min + modulo(y - mesh%eta2_min, mesh%eta2_max - mesh%eta2_min)
  end subroutine apply_periodic_bc_on_cartesian_mesh_2d

end module sll_module_pic_base
