!> @ingroup particle_methods
!> @brief
!> Common interface for groups of particles with remapping capabilities (like FSL or LTPIC)

module sll_m_remapped_pic_base

#include "sll_working_precision.h"
#include "sll_memory.h"

  implicit none
  private

  public :: temp_species_new

  !============================================================================
  !> @brief Particle species (temporary: should use the common type for that)
  !============================================================================
  type, public :: sll_temp_species
    !      character(len=64) :: name !< species name
      sll_real64        :: q    !< charge of a single particle
      sll_real64        :: m    !< mass   of a single particle

  contains
    procedure, pass(self) :: q_over_m => get_q_over_m
!    procedure  :: temp_species_new ! => create_temp_species_new      !! is that a good syntax for the constructor interface?

  end type sll_temp_species

!  interface sll_temp_species
!    procedure, public temp_species_new
!  end interface

  !============================================================================
  !> @brief Particle group
  !============================================================================
  type, public, abstract :: sll_c_remapped_particle_group       ! previous name sll_particle_group_base

    class( sll_temp_species ), pointer :: species
    sll_int32                     :: id

    sll_int32                     :: dimension_x
    sll_int32                     :: dimension_v
    sll_int32                     :: number_particles

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

    ! Remapping
    procedure( rem ), deferred :: remap

    ! Visualize
    procedure( vis ),  deferred :: visualize_f_slice_x_vx

  end type sll_c_remapped_particle_group

  !----------------------------------------------------------------------------
  abstract interface
   pure function i_get_scalar( self, i ) result( r )
    use sll_m_working_precision
    import sll_c_remapped_particle_group
    class( sll_c_remapped_particle_group ), intent( in ) :: self
    sll_int32                       , intent( in ) :: i
    sll_real64 :: r
   end function i_get_scalar
  end interface

  !----------------------------------------------------------------------------
  abstract interface
   pure function i_get_coords( self, i ) result( r )
    use sll_m_working_precision
    import sll_c_remapped_particle_group
    class( sll_c_remapped_particle_group ), intent( in ) :: self
    sll_int32                       , intent( in ) :: i
    sll_real64 :: r(3)
   end function i_get_coords
  end interface

  !----------------------------------------------------------------------------
  abstract interface
   pure function i_get_integer(self, i ) result( i_out )
    use sll_m_working_precision
    import sll_c_remapped_particle_group
    class( sll_c_remapped_particle_group ),   intent( in )    ::  self
    sll_int32,                          intent( in )    ::  i
    sll_int32  ::  i_out

   end function i_get_integer
  end interface

  !----------------------------------------------------------------------------
  abstract interface
   subroutine i_set_coords( self, i, x )
    use sll_m_working_precision
    import sll_c_remapped_particle_group
    class( sll_c_remapped_particle_group ), intent( inout ) :: self
    sll_int32                       , intent( in    ) :: i
    sll_real64                      , intent( in    ) :: x(3)
   end subroutine i_set_coords
  end interface

  !----------------------------------------------------------------------------
  abstract interface
   subroutine i_set_scalar( self, i, s )
    use sll_m_working_precision
    import sll_c_remapped_particle_group
    class( sll_c_remapped_particle_group ), intent( inout ) :: self
    sll_int32                       , intent( in    ) :: i
    sll_real64                      , intent( in    ) :: s
   end subroutine i_set_scalar
  end interface

  !----------------------------------------------------------------------------
  abstract interface
   subroutine set_scalar( self, s )
    use sll_m_working_precision
    import sll_c_remapped_particle_group
    class( sll_c_remapped_particle_group ), intent( inout ) :: self
    sll_real64                      , intent( in    ) :: s
   end subroutine set_scalar
  end interface

  !----------------------------------------------------------------------------
  abstract interface
   subroutine init( self, initial_density_identifier, target_total_charge, enforce_total_charge, rand_seed, rank, world_size)
    use sll_m_working_precision
    import sll_c_remapped_particle_group
    class( sll_c_remapped_particle_group ), intent( inout ) :: self
    sll_int32                       , intent( in    ) :: initial_density_identifier
    sll_real64,                       intent( in )    :: target_total_charge
    logical,                          intent( in )    :: enforce_total_charge
    sll_int32, dimension(:)         , intent( in ), optional :: rand_seed
    sll_int32                       , intent( in ), optional :: rank, world_size

   end subroutine init
  end interface

  !----------------------------------------------------------------------------
  abstract interface
   subroutine rem( self, target_total_charge, enforce_total_charge )
    use sll_m_working_precision
    import sll_c_remapped_particle_group
    class( sll_c_remapped_particle_group ), intent( inout ) :: self
    sll_real64,                       intent( in )    :: target_total_charge
    logical,                          intent( in )    :: enforce_total_charge

   end subroutine rem
  end interface

  !----------------------------------------------------------------------------
  abstract interface
   subroutine no_arg( self )
    use sll_m_working_precision
    import sll_c_remapped_particle_group
    class( sll_c_remapped_particle_group ), intent( inout ) :: self

   end subroutine no_arg
  end interface

  !----------------------------------------------------------------------------
  abstract interface
   !> plots a 2d slice, using a 4d computational grid
   subroutine vis(self, array_name, plot_np_x, plot_np_y, plot_np_vx, plot_np_vy, iplot)
    use sll_m_working_precision
    import sll_c_remapped_particle_group
    class( sll_c_remapped_particle_group ), intent( inout ) :: self
    character(len=*),                   intent(in)      :: array_name   !< field name
    sll_int32,                          intent(in)      :: plot_np_x    !< nb of points in the x  plotting grid (see comment above)
    sll_int32,                          intent(in)      :: plot_np_y    !< nb of points in the y  plotting grid (see comment above)
    sll_int32,                          intent(in)      :: plot_np_vx   !< nb of points in the vx plotting grid (see comment above)
    sll_int32,                          intent(in)      :: plot_np_vy   !< nb of points in the vy plotting grid (see comment above)
    sll_int32,                          intent(in)      :: iplot        !< plot counter

   end subroutine vis
  end interface

  !----------------------------------------------------------------------------
  abstract interface
   subroutine set_landau_params( self, thermal_speed, alpha, k_landau )
    use sll_m_working_precision
    import sll_c_remapped_particle_group
    class( sll_c_remapped_particle_group ), intent( inout ) :: self
    sll_real64                      , intent( in    ) :: thermal_speed
    sll_real64                      , intent( in    ) :: alpha
    sll_real64                      , intent( in    ) :: k_landau

   end subroutine set_landau_params
  end interface

  !----------------------------------------------------------------------------
  abstract interface
   subroutine set_hat_f0_params( self, x0, y0, vx0, vy0, r_x, r_y, r_vx, r_vy, basis_height, shift )
    use sll_m_working_precision
    import sll_c_remapped_particle_group
    class( sll_c_remapped_particle_group ), intent( inout ) :: self
    sll_real64,                     intent(in)      :: x0
    sll_real64,                     intent(in)      :: y0
    sll_real64,                     intent(in)      :: vx0
    sll_real64,                     intent(in)      :: vy0
    sll_real64,                     intent(in)      :: r_x
    sll_real64,                     intent(in)      :: r_y
    sll_real64,                     intent(in)      :: r_vx
    sll_real64,                     intent(in)      :: r_vy
    sll_real64,                     intent(in)      :: basis_height
    sll_real64,                     intent(in)      :: shift

   end subroutine set_hat_f0_params
  end interface


  !----------------------------------------------------------------------------
  abstract interface
   subroutine dep_charge_2d( self, charge_accumulator, target_total_charge, enforce_total_charge )
    use sll_m_working_precision
#include "sll_accumulators.h"
    import sll_c_remapped_particle_group
    class( sll_c_remapped_particle_group ),     intent( inout ) :: self
    type( sll_charge_accumulator_2d ), pointer, intent( inout ) :: charge_accumulator
    sll_real64,                                 intent(in)      :: target_total_charge
    logical,                                    intent(in)      :: enforce_total_charge
   end subroutine dep_charge_2d
  end interface

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
contains
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  function temp_species_new( &
      species_charge,     &
      species_mass        &
  ) result(res)
    sll_real64, intent ( in )   :: species_charge
    sll_real64, intent ( in )   :: species_mass
    type(sll_temp_species), pointer  :: res
    sll_int32  :: ierr

    SLL_ALLOCATE( res, ierr )
    !    res%name =
    res%q = species_charge
    res%m = species_mass

  end function temp_species_new


  !----------------------------------------------------------------------------
  function get_q_over_m( self ) result( r )
    class( sll_temp_species ), intent( in ) :: self
    sll_real64 :: r

    r = self%q / self%m

  end function get_q_over_m

end module sll_m_remapped_pic_base
