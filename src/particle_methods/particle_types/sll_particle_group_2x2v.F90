!> @ingroup particle_methods
!> @author Katharina Kormann
!> @brief Simple particle group type for 2x2v.
!> @details ...
module sll_m_particle_group_2x2v

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use sll_working_precision
  use sll_module_pic_base

  implicit none

!> Simple version of a PIC particle group in 2x2v
type, extends(sll_particle_group_base) :: sll_particle_group_2x2v
   !sll_int32               :: n_particles !< number of particle
   sll_real64, pointer :: particle_array(:,:) !< array of particles

contains
    ! Getters
    procedure :: get_x  => get_x_2x2v
    procedure :: get_v  => get_v_2x2v
    procedure :: get_charge => get_charge_2x2v
    procedure :: get_mass => get_mass_2x2v

    ! Setters
    procedure :: set_x => set_x_2x2v
    procedure :: set_v => set_v_2x2v
    procedure :: set_weight => set_weight_2x2v

    ! Initializer
    procedure :: initialize => initialize_particle_group_2x2v 


   
end type sll_particle_group_2x2v

contains

  !----------------------------------------------------------------------!
  !> Constructor
  function sll_new_particle_group_2x2v(n_particles, n_total_particles, charge, mass) result(self)
    class( sll_particle_group_2x2v ),  pointer :: self
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

  end function sll_new_particle_group_2x2v
  
  !----------------------------------------------------------------------!
  pure function get_x_2x2v( self, i ) result( r )
    class( sll_particle_group_2x2v ), intent( in ) :: self 
    sll_int32                       , intent( in ) :: i
    sll_real64 :: r(3)

    r(1:2) = self%particle_array(i, 1:2)
    
  end function get_x_2x2v

  !----------------------------------------------------------------------!
  pure function get_v_2x2v( self, i ) result( r )
    class( sll_particle_group_2x2v ), intent( in ) :: self
    sll_int32                       , intent( in ) :: i
    sll_real64 :: r(3)

    r(1:2) = self%particle_array(i, 3:4)
    
  end function get_v_2x2v

  !----------------------------------------------------------------------!
  pure function get_charge_2x2v( self, i ) result (r)
        class( sll_particle_group_2x2v ), intent( in ) :: self
    sll_int32                       , intent( in ) :: i
    sll_real64 :: r

    r = self%species%q * self%particle_array(i,5)

  end function get_charge_2x2v

  !----------------------------------------------------------------------!
  pure function get_mass_2x2v( self, i) result (r)
        class( sll_particle_group_2x2v ), intent( in ) :: self
    sll_int32                       , intent( in ) :: i
    sll_real64 :: r

    r = self%species%m * self%particle_array(i,5)

  end function get_mass_2x2v

  !----------------------------------------------------------------------!
  subroutine set_x_2x2v( self, i, x )
    class( sll_particle_group_2x2v ), intent( inout ) :: self
    sll_int32                       , intent( in ) :: i
    sll_real64                      , intent( in):: x(3)

    self%particle_array(i, 1:2) = x(1:2)
    
  end subroutine set_x_2x2v

  subroutine set_v_2x2v( self, i, x )
    class( sll_particle_group_2x2v ), intent( inout ) :: self
    sll_int32                       , intent( in ) :: i
    sll_real64                      , intent( in):: x(3)

    self%particle_array(i, 3:4) = x(1:2)
    
  end subroutine set_v_2x2v
  
  !----------------------------------------------------------------------!
  subroutine set_weight_2x2v( self, i, x )
    class( sll_particle_group_2x2v ), intent( inout ) :: self
    sll_int32                       , intent( in ) :: i
    sll_real64                      , intent( in):: x

    self%particle_array(i, 5) = x
    
  end subroutine set_weight_2x2v

  !----------------------------------------------------------------------!
  subroutine initialize_particle_group_2x2v (self, n_particles)
     class( sll_particle_group_2x2v ),  intent( inout ) :: self
     sll_int32                       , intent( in )    :: n_particles

     sll_int32                                         :: ierr

    self%n_particles = n_particles
    SLL_ALLOCATE(self%particle_array(n_particles,5), ierr) 
    self%species%q = 1.0
    self%species%m = 1.0

  end subroutine initialize_particle_group_2x2v


end module sll_m_particle_group_2x2v
