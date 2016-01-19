!> @ingroup particle_types
!> @author Katharina Kormann, IPP
!> @brief Simple particle group type for 2d2v.
!> @details ...
module sll_m_particle_group_2d2v

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_particle_group_base, only: &
    sll_c_particle_group_base, &
    sll_t_species

  implicit none

  public :: &
    sll_f_new_particle_group_2d2v, &
    sll_s_new_particle_group_2d2v, &
    sll_s_new_particle_group_2d2v_ptr, &
    sll_t_particle_group_2d2v

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!> Simple version of a PIC particle group in 2d2v
type, extends(sll_c_particle_group_base) :: sll_t_particle_group_2d2v
   !sll_int32               :: n_particles !< number of particle
   sll_real64, allocatable :: particle_array(:,:) !< array of particles
   sll_real64 :: common_weight = 1.0_f64

contains
    ! Getters
    procedure :: get_x  => get_x_2d2v !> Get the values of the coordinate of a particle
    procedure :: get_v  => get_v_2d2v !> Get the values of the velocity of a particle
    procedure :: get_charge => get_charge_2d2v !> Get the charge(s)
    procedure :: get_mass => get_mass_2d2v !> Get the mass(es)
    procedure :: get_weights => get_weights_2d2v !> Get weight(s) of the particle

    ! Setters
    procedure :: set_x => set_x_2d2v !> Set the values of the coordinate of a particle
    procedure :: set_v => set_v_2d2v !> Set the values of the velocity of a particle
    procedure :: set_weights => set_weights_2d2v !> Set the weight(s) of a particle
    procedure :: set_common_weight => set_common_weight_2d2v !> Set the common weight for the particle

    ! Initializer
    procedure :: initialize => initialize_particle_group_2d2v  !> Initialization function
    procedure :: delete => delete_particle_group_2d2v !> Destructor


end type sll_t_particle_group_2d2v

contains

  !----------------------------------------------------------------------!
  !> Destructor
  subroutine delete_particle_group_2d2v(self)
    class( sll_t_particle_group_2d2v ), intent( inout ) :: self  !< particle group

    deallocate(self%particle_array)

  end subroutine delete_particle_group_2d2v

  !----------------------------------------------------------------------!
  !> Initialize particle group
  subroutine initialize_particle_group_2d2v (self, n_particles, n_total_particles, charge, mass, n_weights)
    class( sll_t_particle_group_2d2v ), intent( inout ) :: self  !< particle group
    sll_int32                       , intent( in )    :: n_particles !< number of particles local to the processor
    sll_int32                       , intent( in )    :: n_total_particles !< number of particles in total simulation
    sll_real64                      , intent( in )    :: charge !< charge of the particle species
    sll_real64                      , intent( in )    :: mass   !< mass of the particle species
    sll_int32                       , intent( in )    :: n_weights !< number of weights
     
    sll_int32                                         :: ierr

    self%n_particles = n_particles
    self%n_total_particles = n_total_particles

    SLL_ALLOCATE(self%particle_array(4+n_weights, self%n_particles), ierr) 

    allocate(self%species, stat=ierr)
    SLL_ASSERT( ierr == 0)
    call self%species%initialize( charge, mass)

    self%n_weights = n_weights

  end subroutine initialize_particle_group_2d2v


  !----------------------------------------------------------------------!
  !> Constructor (legacy version)
  function sll_f_new_particle_group_2d2v(n_particles, n_total_particles, charge, mass, n_weights) result(self)
    class( sll_t_particle_group_2d2v ),  pointer        :: self
    sll_int32                       , intent( in )    :: n_particles !< number of particles local to the processor
    sll_int32                       , intent( in )    :: n_total_particles !< number of particles in total simulation
    sll_real64                      , intent( in )    :: charge !< charge of the particle species
    sll_real64                      , intent( in )    :: mass   !< mass of the particle species
    sll_int32                       , intent(in)      :: n_weights !< number of weights
    
    sll_int32                                         :: ierr

    SLL_ALLOCATE(self, ierr)
    call self%initialize(n_particles, n_total_particles, charge, mass, n_weights)

  end function sll_f_new_particle_group_2d2v

  !----------------------------------------------------------------------!
  !> Constructor for abstract type
  subroutine sll_s_new_particle_group_2d2v_ptr(particle_group, n_particles, n_total_particles, charge, mass, n_weights)
    class( sll_c_particle_group_base ),  pointer, intent( out )  :: particle_group
    sll_int32                       , intent( in )    :: n_particles !< number of particles local to the processor
    sll_int32                       , intent( in )    :: n_total_particles !< number of particles in total simulation
    sll_real64                      , intent( in )    :: charge !< charge of the particle species
    sll_real64                      , intent( in )    :: mass   !< mass of the particle species
    sll_int32                       , intent(in)      :: n_weights !< number of weights
    
    sll_int32                                         :: ierr

    SLL_ALLOCATE( sll_t_particle_group_2d2v :: particle_group, ierr)
    
    select type( particle_group )
    type is ( sll_t_particle_group_2d2v )
       call particle_group%initialize(n_particles, n_total_particles, charge, mass, n_weights)
    end select

  end subroutine sll_s_new_particle_group_2d2v_ptr


  !----------------------------------------------------------------------!
  !> Constructor for abstract type
  subroutine sll_s_new_particle_group_2d2v(particle_group, n_particles, n_total_particles, charge, mass, n_weights)
    class( sll_c_particle_group_base ),  allocatable, intent( out )  :: particle_group
    sll_int32                       , intent( in )    :: n_particles !< number of particles local to the processor
    sll_int32                       , intent( in )    :: n_total_particles !< number of particles in total simulation
    sll_real64                      , intent( in )    :: charge !< charge of the particle species
    sll_real64                      , intent( in )    :: mass   !< mass of the particle species
    sll_int32                       , intent(in)      :: n_weights !< number of weights
    
    sll_int32                                         :: ierr

    SLL_ALLOCATE( sll_t_particle_group_2d2v :: particle_group, ierr)
    
    select type( particle_group )
    type is ( sll_t_particle_group_2d2v )
       call particle_group%initialize(n_particles, n_total_particles, charge, mass, n_weights)
    end select

  end subroutine sll_s_new_particle_group_2d2v
  
  !----------------------------------------------------------------------!
  !> Get positions of particle \a i
  pure function get_x_2d2v( self, i ) result( r )
    class( sll_t_particle_group_2d2v ), intent( in ) :: self  !< particle group
    sll_int32                       , intent( in ) :: i !< no. of the particle
    sll_real64 :: r(3) !< first two components hold the value of the particle position

    r = 1.0_f64
    r(1:2) = self%particle_array(1:2, i)
    
  end function get_x_2d2v

  !----------------------------------------------------------------------!
  !> Get velocities of particle \a i
  pure function get_v_2d2v( self, i ) result( r )
    class( sll_t_particle_group_2d2v ), intent( in ) :: self  !< particle group
    sll_int32                       , intent( in ) :: i !< no. of the particle
    sll_real64 :: r(3) !< first two components hold the value of the particle velocity

    r = 1.0_f64
    r(1:2) = self%particle_array(3:4, i)
    
  end function get_v_2d2v

! Old version without optional argument
!!$  !----------------------------------------------------------------------!
!!$  pure function get_charge_2d2v( self, i ) result (r)
!!$        class( sll_t_particle_group_2d2v ), intent( in ) :: self  !< particle group
!!$    sll_int32                       , intent( in ) :: i !< no. of the particle
!!$    sll_real64 :: r(self%n_weights) !< particle charge(s)
!!$
!!$    r = self%species%q * self%particle_array(i,5:4+self%n_weights)
!!$
!!$  end function get_charge_2d2v
!!$
!!$  !----------------------------------------------------------------------!
!!$  pure function get_mass_2d2v( self, i) result (r)
!!$        class( sll_t_particle_group_2d2v ), intent( in ) :: self  !< particle group
!!$    sll_int32                       , intent( in ) :: i !< no. of the particle
!!$    sll_real64 :: r(self%n_weights) !< particle mass(es)
!!$
!!$    r = self%species%m * self%particle_array(i,5:4+self%n_weights)
!!$
!!$  end function get_mass_2d2v

  !----------------------------------------------------------------------!
  !> Get charge of particle \a i ( q * particle_weight)
  pure function get_charge_2d2v( self, i , i_weight) result (r)
        class( sll_t_particle_group_2d2v ), intent( in ) :: self !< particle group
    sll_int32                           , intent( in ) :: i !< no. of the particle
    sll_int32, optional                 , intent( in ) :: i_weight !< index of weight to be used (default: 1)
    sll_real64 :: r !< charges(s) of particle i

    sll_int32 :: i_wi

    i_wi = 1
    if(present(i_weight)) i_wi = i_weight
    r = self%species%q  * self%particle_array(4+i_wi, i) * self%common_weight

  end function get_charge_2d2v


  !----------------------------------------------------------------------!
  !> Get mass of particle \a i ( m * particle_weight)
  pure function get_mass_2d2v( self, i, i_weight) result (r)
        class( sll_t_particle_group_2d2v ), intent( in ) :: self !< particle group
    sll_int32                           , intent( in ) :: i !< no. of the particle
    sll_int32, optional                 , intent( in ) :: i_weight !< index of weight to be used (default: 1)
    sll_real64 :: r !< masses(s) of particle i

    sll_int32 :: i_wi

    i_wi = 1
    if(present(i_weight)) i_wi = i_weight
    r = self%species%m * self%particle_array(4+i_wi, i) * self%common_weight

  end function get_mass_2d2v

 
  !----------------------------------------------------------------------!
  !> Get weights of particle \a i 
  pure function get_weights_2d2v( self, i) result (r)
        class( sll_t_particle_group_2d2v ), intent( in ) :: self  !< particle group
    sll_int32                       , intent( in ) :: i !< no. of the particle
    sll_real64 :: r(self%n_weights) !< particle mass(es)

    r = self%species%m * self%particle_array(5:4+self%n_weights, i)

  end function get_weights_2d2v

  !----------------------------------------------------------------------!
  !> Set positions of particle \a i 
  subroutine set_x_2d2v( self, i, x )
    class( sll_t_particle_group_2d2v ), intent( inout ) :: self  !< particle group
    sll_int32                       , intent( in ) :: i !< no. of the particle
    sll_real64                      , intent( in):: x(3) !< components 1 and 2 hold the particle position to be set

    self%particle_array(1:2, i) = x(1:2)
    
  end subroutine set_x_2d2v

  !----------------------------------------------------------------------!
  !> Set velocities of particle \a i 
  subroutine set_v_2d2v( self, i, x )
    class( sll_t_particle_group_2d2v ), intent( inout ) :: self  !< particle group
    sll_int32                       , intent( in ) :: i !< no. of the particle
    sll_real64                      , intent( in):: x(3) !< component 1 and 2 hold the particle velocity to be set

    self%particle_array(3:4, i) = x(1:2)
    
  end subroutine set_v_2d2v
  
  !----------------------------------------------------------------------!
  !> Set weights of particle \a i 
  subroutine set_weights_2d2v( self, i, x )
    class( sll_t_particle_group_2d2v ), intent( inout ) :: self  !< particle group
    sll_int32                       , intent( in ) :: i !< no. of the particle
    sll_real64                      , intent( in):: x(self%n_weights) !< particle weight(s) to be set

    self%particle_array(5:4+self%n_weights, i) = x
    
  end subroutine set_weights_2d2v

  !----------------------------------------------------------------------!
  !> Set the common weight
  subroutine set_common_weight_2d2v( self, x )
    class( sll_t_particle_group_2d2v ), intent( inout ) :: self  !< particle group
    sll_real64                      , intent( in):: x !< common weight

    self%common_weight = x

  end subroutine set_common_weight_2d2v


  
end module sll_m_particle_group_2d2v
