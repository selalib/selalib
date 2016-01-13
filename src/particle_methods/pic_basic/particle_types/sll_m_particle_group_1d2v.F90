!> @ingroup particle_types
!> @author Katharina Kormann, IPP
!> @brief Simple particle group type for 1d2v.
!> @details ...
module sll_m_particle_group_1d2v

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_particle_group_base, only: &
    sll_c_particle_group_base, &
    sll_f_species_new

  implicit none

  public :: &
    sll_f_new_particle_group_1d2v, &
    sll_t_particle_group_1d2v

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!> Simple version of a PIC particle group in 1d2v
type, extends(sll_c_particle_group_base) :: sll_t_particle_group_1d2v
   !sll_int32               :: n_particles !< number of particle
   sll_real64, pointer :: particle_array(:,:) !< array of particles
   sll_real64 :: common_weight = 1.0_f64

contains
    ! Getters
    procedure :: get_x  => get_x_1d2v !> Get the values of the coordinate of a particle
    procedure :: get_v  => get_v_1d2v !> Get the values of the velocity of a particle
    procedure :: get_charge => get_charge_1d2v !> Get the charge(s)
    procedure :: get_mass => get_mass_1d2v !> Get the mass(es)
    procedure :: get_weights => get_weights_1d2v !> Get the common weight (not used for this particle group)

    ! Setters
    procedure :: set_x => set_x_1d2v !> Set the values of the coordinate of a particle
    procedure :: set_v => set_v_1d2v !> Set the values of the velocity of a particle
    procedure :: set_weights => set_weight_1d2v !> Set the weight(s) of a particle
    procedure :: set_common_weight => set_common_weight_1d2v !> Set the common weight for the particle

    ! Initializer
    procedure :: initialize => initialize_particle_group_1d2v !> Initialization function
    procedure :: delete => delete_particle_group_1d2v !> Destructor

   
 end type sll_t_particle_group_1d2v

contains

  !----------------------------------------------------------------------!
  !> Destructor
  subroutine delete_particle_group_1d2v(self)
    class( sll_t_particle_group_1d2v ), intent(inout) :: self

    deallocate(self%particle_array)

  end subroutine delete_particle_group_1d2v

    !----------------------------------------------------------------------!
  subroutine initialize_particle_group_1d2v (self, n_particles, n_total_particles, charge, mass, n_weights)
    class( sll_t_particle_group_1d2v ), intent( inout) :: self !< particle group 
    sll_int32                       , intent( in )    :: n_particles !< number of particles local to the processor
    sll_int32                       , intent( in )    :: n_total_particles !< number of particles in total simulation
    sll_real64                      , intent( in )    :: charge !< charge of the particle species
    sll_real64                      , intent( in )    :: mass   !< mass of the particle species(self, n_particles)
    sll_int32                       , intent( in )    :: n_weights !< number of weights

    sll_int32                                         :: ierr 

    self%n_particles = n_particles
    self%n_total_particles = n_total_particles

    SLL_ALLOCATE(self%particle_array(3+n_weights, n_particles), ierr) 

    self%species => sll_f_species_new( charge, mass )

    self%n_weights = n_weights

  end subroutine initialize_particle_group_1d2v

  !----------------------------------------------------------------------!
  !> Constructor
  function sll_f_new_particle_group_1d2v(n_particles, n_total_particles, charge, mass, n_weights) result(self)
    class( sll_t_particle_group_1d2v ),  pointer :: self
    sll_int32                       , intent( in )    :: n_particles !< number of particles local to the processor
    sll_int32                       , intent( in )    :: n_total_particles !< number of particles in total simulation
    sll_real64                      , intent( in )    :: charge !< charge of the particle species
    sll_real64                      , intent( in )    :: mass   !< mass of the particle species
    sll_int32                       , intent( in )    :: n_weights !< number of weights
    
    sll_int32                                         :: ierr

    SLL_ALLOCATE(self, ierr)

    call self%initialize( n_particles, n_total_particles, charge, mass, n_weights )
    
  end function sll_f_new_particle_group_1d2v
  
  !----------------------------------------------------------------------!
  pure function get_x_1d2v( self, i ) result( r )
    class( sll_t_particle_group_1d2v ), intent( in ) :: self !< particle group
    sll_int32                       , intent( in ) :: i !< no. of the particle
    sll_real64 :: r(3) !< position of particle i

    r = 1.0_f64
    r(1) = self%particle_array(1, i)
    
  end function get_x_1d2v

  !----------------------------------------------------------------------!
  pure function get_v_1d2v( self, i ) result( r )
    class( sll_t_particle_group_1d2v ), intent( in ) :: self !< particle group
    sll_int32                       , intent( in ) :: i !< no. of the particle
    sll_real64 :: r(3)

    r = 1.0_f64
    r(1:2) = self%particle_array( 2:3, i )
    
  end function get_v_1d2v

! Old version without optional argument
!!$  !----------------------------------------------------------------------!
!!$  pure function get_charge_1d2v( self, i ) result (r)
!!$        class( sll_t_particle_group_1d2v ), intent( in ) :: self !< particle group
!!$    sll_int32                       , intent( in ) :: i !< no. of the particle
!!$    sll_real64 :: r(self%n_weights) !< charges(s) of particle i
!!$
!!$    r = self%species%q  * self%particle_array(i, 4:3+self%n_weights)
!!$
!!$  end function get_charge_1d2v
!!$
!!$
!!$  !----------------------------------------------------------------------!
!!$  pure function get_mass_1d2v( self, i) result (r)
!!$        class( sll_t_particle_group_1d2v ), intent( in ) :: self !< particle group
!!$    sll_int32                       , intent( in ) :: i !< no. of the particle
!!$    sll_real64 :: r(self%n_weights) !< masses(s) of particle i
!!$
!!$    r = self%species%m * self%particle_array(i, 4:3+self%n_weights)
!!$
!!$  end function get_mass_1d2v

  !----------------------------------------------------------------------!
  pure function get_charge_1d2v( self, i , i_weight) result (r)
        class( sll_t_particle_group_1d2v ), intent( in ) :: self !< particle group
    sll_int32                           , intent( in ) :: i !< no. of the particle
    sll_int32, optional                 , intent( in ) :: i_weight
    sll_real64 :: r !< charges(s) of particle i

    sll_int32 :: i_wi

    i_wi = 1
    if(present(i_weight)) i_wi = i_weight
    r = self%species%q  * self%particle_array(3+i_wi, i) * self%common_weight

  end function get_charge_1d2v


  !----------------------------------------------------------------------!
  pure function get_mass_1d2v( self, i, i_weight) result (r)
        class( sll_t_particle_group_1d2v ), intent( in ) :: self !< particle group
    sll_int32                           , intent( in ) :: i !< no. of the particle
    sll_int32, optional                 , intent( in ) :: i_weight 
    sll_real64 :: r !< masses(s) of particle i

    sll_int32 :: i_wi

    i_wi = 1
    if(present(i_weight)) i_wi = i_weight
    r = self%species%m * self%particle_array( 3+i_wi, i) * self%common_weight

  end function get_mass_1d2v

  !----------------------------------------------------------------------!
  pure function get_weights_1d2v( self, i) result (r)
        class( sll_t_particle_group_1d2v ), intent( in ) :: self !< particle group
    sll_int32                       , intent( in ) :: i !< no. of the particle
    sll_real64 :: r(self%n_weights) !< weight(s) of particle i

    r = self%particle_array(4:3+self%n_weights, i)

  end function get_weights_1d2v

  !----------------------------------------------------------------------!
  subroutine set_x_1d2v( self, i, x )
    class( sll_t_particle_group_1d2v ), intent( inout ) :: self !< particle group
    sll_int32                       , intent( in ) :: i !< no. of the particle
    sll_real64                      , intent( in):: x(3) !< first component holds the value of the position to be set

    self%particle_array(1, i) = x(1)
    
  end subroutine set_x_1d2v

  subroutine set_v_1d2v( self, i, x )
    class( sll_t_particle_group_1d2v ), intent( inout ) :: self !< particle group
    sll_int32                       , intent( in ) :: i !< no. of the particle
    sll_real64                      , intent( in):: x(3) !< first two components hold the values of the velocity to be set

    self%particle_array(2:3, i) = x(1:2)
    
  end subroutine set_v_1d2v
  
  !----------------------------------------------------------------------!
  subroutine set_weight_1d2v( self, i, x )
    class( sll_t_particle_group_1d2v ), intent( inout ) :: self !< particle group
    sll_int32                       , intent( in ) :: i !< no. of the particle
    sll_real64                      , intent( in):: x(self%n_weights) !< particle weight(s)

    self%particle_array(4:3+self%n_weights, i) = x
    
  end subroutine set_weight_1d2v

  
  !----------------------------------------------------------------------!
  subroutine set_common_weight_1d2v( self, x )
    class( sll_t_particle_group_1d2v ), intent( inout ) :: self !< particle group
    sll_real64                      , intent( in):: x !< common weight

    self%common_weight = x
    
  end subroutine set_common_weight_1d2v




end module sll_m_particle_group_1d2v
