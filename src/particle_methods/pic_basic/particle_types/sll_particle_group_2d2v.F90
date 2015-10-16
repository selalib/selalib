!> @ingroup particle_types
!> @author Katharina Kormann, IPP
!> @brief Simple particle group type for 2d2v.
!> @details ...
module sll_m_particle_group_2d2v

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use sll_working_precision
  use sll_m_pic_base

  implicit none
  private
  
  public :: sll_new_particle_group_2d2v

!> Simple version of a PIC particle group in 2d2v
type, public, extends(sll_particle_group_base) :: sll_particle_group_2d2v
   !sll_int32               :: n_particles !< number of particle
   sll_real64, pointer :: particle_array(:,:) !< array of particles

contains
    ! Getters
    procedure :: get_x  => get_x_2d2v !> Get the values of the coordinate of a particle
    procedure :: get_v  => get_v_2d2v !> Get the values of the velocity of a particle
    procedure :: get_charge => get_charge_2d2v !> Get the charge
    procedure :: get_mass => get_mass_2d2v !> Get the mass
    procedure :: get_weights => get_weights_2d2v !> Get the particle weight(s)
    procedure :: get_common_weight => get_common_weight_2d2v !> Get the common weight (not used for this particle group)

    ! Setters
    procedure :: set_x => set_x_2d2v !> Set the values of the coordinate of a particle
    procedure :: set_v => set_v_2d2v !> Set the values of the velocity of a particle
    procedure :: set_weights => set_weights_2d2v !> Set the weight(s) of a particle
    procedure :: set_common_weight => set_common_weight_2d2v !> Set the common weight for the particle

    ! Initializer
    procedure :: initialize => initialize_particle_group_2d2v  !> Initialization function
   
end type sll_particle_group_2d2v

contains

  !----------------------------------------------------------------------!
  !> Constructor
  function sll_new_particle_group_2d2v(n_particles, n_total_particles, charge, mass) result(self)
    class( sll_particle_group_2d2v ),  pointer :: self
    sll_int32                       , intent( in )    :: n_particles !< number of particles local to the processor
    sll_int32                       , intent( in )    :: n_total_particles !< number of particles in total simulation
    sll_real64                      , intent( in )    :: charge !< charge of the particle species
    sll_real64                      , intent( in )    :: mass   !< mass of the particle species
    
    sll_int32                                         :: ierr

    SLL_ALLOCATE(self, ierr)
    self%n_particles = n_particles
    self%n_total_particles = n_total_particles
    SLL_ALLOCATE(self%particle_array(self%n_particles,5), ierr) 
    self%species => species_new( charge, mass)

    self%n_weights = 1

  end function sll_new_particle_group_2d2v
  
  !----------------------------------------------------------------------!
  pure function get_x_2d2v( self, i ) result( r )
    class( sll_particle_group_2d2v ), intent( in ) :: self  !< particle group
    sll_int32                       , intent( in ) :: i !< no. of the particle
    sll_real64 :: r(3) !< first two components hold the value of the particle position

    r = 1.0_f64
    r(1:2) = self%particle_array(i, 1:2)
    
  end function get_x_2d2v

  !----------------------------------------------------------------------!
  pure function get_v_2d2v( self, i ) result( r )
    class( sll_particle_group_2d2v ), intent( in ) :: self  !< particle group
    sll_int32                       , intent( in ) :: i !< no. of the particle
    sll_real64 :: r(3) !< first two components hold the value of the particle velocity

    r = 1.0_f64
    r(1:2) = self%particle_array(i, 3:4)
    
  end function get_v_2d2v

  !----------------------------------------------------------------------!
  pure function get_charge_2d2v( self, i ) result (r)
        class( sll_particle_group_2d2v ), intent( in ) :: self  !< particle group
    sll_int32                       , intent( in ) :: i !< no. of the particle
    sll_real64 :: r !< charge

    r = self%species%q

  end function get_charge_2d2v

  !----------------------------------------------------------------------!
  pure function get_mass_2d2v( self, i) result (r)
        class( sll_particle_group_2d2v ), intent( in ) :: self  !< particle group
    sll_int32                       , intent( in ) :: i !< no. of the particle
    sll_real64 :: r !< mass

    r = self%species%m

  end function get_mass_2d2v
  !----------------------------------------------------------------------!

  pure function get_common_weight_2d2v( self, i) result (r)
        class( sll_particle_group_2d2v ), intent( in ) :: self  !< particle group
    sll_int32                       , intent( in ) :: i !< no. of the particle
    sll_real64 :: r !< common weight

    r = 1.0_f64

  end function get_common_weight_2d2v

  !----------------------------------------------------------------------!
  pure function get_weights_2d2v( self, i) result (r)
        class( sll_particle_group_2d2v ), intent( in ) :: self  !< particle group
    sll_int32                       , intent( in ) :: i !< no. of the particle
    sll_real64 :: r(self%n_weights) !< particle weight(s)

    r = self%species%q_over_m() * self%particle_array(i,5)

  end function get_weights_2d2v


  !----------------------------------------------------------------------!
  subroutine set_x_2d2v( self, i, x )
    class( sll_particle_group_2d2v ), intent( inout ) :: self  !< particle group
    sll_int32                       , intent( in ) :: i !< no. of the particle
    sll_real64                      , intent( in):: x(3) !< components 1 and 2 hold the particle position to be set

    self%particle_array(i, 1:2) = x(1:2)
    
  end subroutine set_x_2d2v

  !----------------------------------------------------------------------!
  subroutine set_v_2d2v( self, i, x )
    class( sll_particle_group_2d2v ), intent( inout ) :: self  !< particle group
    sll_int32                       , intent( in ) :: i !< no. of the particle
    sll_real64                      , intent( in):: x(3) !< component 1 and 2 hold the particle velocity to be set

    self%particle_array(i, 3:4) = x(1:2)
    
  end subroutine set_v_2d2v
  
  !----------------------------------------------------------------------!
  subroutine set_weights_2d2v( self, i, x )
    class( sll_particle_group_2d2v ), intent( inout ) :: self  !< particle group
    sll_int32                       , intent( in ) :: i !< no. of the particle
    sll_real64                      , intent( in):: x(self%n_weights) !< particle weight(s) to be set

    self%particle_array(i, 5) = x(1)
    
  end subroutine set_weights_2d2v

  !----------------------------------------------------------------------!
  subroutine set_common_weight_2d2v( self, i, x )
    class( sll_particle_group_2d2v ), intent( inout ) :: self  !< particle group
    sll_int32                       , intent( in ) :: i !< no. of the particle
    sll_real64                      , intent( in):: x !< common weight

  end subroutine set_common_weight_2d2v

  !----------------------------------------------------------------------!
  subroutine initialize_particle_group_2d2v (self, n_particles, n_total_particles, charge, mass)
    class( sll_particle_group_2d2v ),  intent( inout ) :: self  !< particle group
    sll_int32                       , intent( in )    :: n_particles !< number of particles local to the processor
    sll_int32                       , intent( in )    :: n_total_particles !< number of particles in total simulation
    sll_real64                      , intent( in )    :: charge !< charge of the particle species
    sll_real64                      , intent( in )    :: mass   !< mass of the particle species
     
    sll_int32                                         :: ierr

    self%n_particles = n_particles
    self%n_total_particles = n_total_particles
    SLL_ALLOCATE(self%particle_array(self%n_particles,5), ierr) 
    self%species => species_new( charge, mass)

  end subroutine initialize_particle_group_2d2v

  
end module sll_m_particle_group_2d2v
