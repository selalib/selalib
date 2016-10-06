!> @ingroup particle_groups
!> @author Benedikt Perse, Ipp
!> @brief Simple particle group type for 1d1v.
!> @details ...
module sll_m_particle_group_1d1v

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_particle_group_base, only: &
    sll_c_particle_group_base, &
    sll_t_species

  implicit none

public :: &
    sll_s_new_particle_group_1d1v, &
    sll_s_new_particle_group_1d1v_ptr, &
    sll_t_particle_group_1d1v


  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Simple version of a PIC particle group in 1d1v
type, extends(sll_c_particle_group_base) :: sll_t_particle_group_1d1v
   !sll_int32               :: n_particles !< number of particle
   sll_real64, allocatable :: particle_array(:,:) !< array of particles
   sll_real64 :: common_weight = 1.0_f64

contains
    ! Getters
    procedure :: get_x  => get_x_1d1v !> Get the values of the coordinate of a particle
    procedure :: get_v  => get_v_1d1v !> Get the values of the velocity of a particle
    procedure :: get_charge => get_charge_1d1v !> Get the charge(s)
    procedure :: get_mass => get_mass_1d1v !> Get the mass(es)
    procedure :: get_weights => get_weights_1d1v !> Get the common weight (not used for this particle group)

    ! Setters
    procedure :: set_x => set_x_1d1v !> Set the values of the coordinate of a particle
    procedure :: set_v => set_v_1d1v !> Set the values of the velocity of a particle
    procedure :: set_weights => set_weight_1d1v !> Set the weight(s) of a particle
    procedure :: set_common_weight => set_common_weight_1d1v !> Set the common weight for the particle

    ! Initializer
    procedure :: init => initialize_particle_group_1d1v !> Initialization function
    procedure :: free => delete_particle_group_1d1v !> Destructor

    procedure :: print => print_particle_group_1d1v !> print particle array to file
   
 end type sll_t_particle_group_1d1v

contains

!----------------------------------------------------------------------!
  !> Destructor
  subroutine delete_particle_group_1d1v(self)
    class( sll_t_particle_group_1d1v ), intent(inout) :: self !< particle group

    deallocate(self%particle_array)

  end subroutine delete_particle_group_1d1v

    !----------------------------------------------------------------------!
  !> Initialization of the particle group
  subroutine initialize_particle_group_1d1v ( &
       self, &
       n_particles, &
       n_total_particles, &
       charge, &
       mass, &
       n_weights)
    class( sll_t_particle_group_1d1v ), intent( inout)  :: self !< particle group 
    sll_int32                         , intent( in )    :: n_particles !< number of particles local to the processor
    sll_int32                         , intent( in )    :: n_total_particles !< number of particles in total simulation
    sll_real64                        , intent( in )    :: charge !< charge of the particle species
    sll_real64                        , intent( in )    :: mass   !< mass of the particle species(self, n_particles)
    sll_int32                         , intent( in )    :: n_weights !< number of weights

    sll_int32                                           :: ierr 

    self%n_particles = n_particles
    self%n_total_particles = n_total_particles

    SLL_ALLOCATE(self%particle_array(2+n_weights, n_particles), ierr) 

    allocate(self%species, stat=ierr)
    SLL_ASSERT( ierr == 0)
    call self%species%init( charge, mass )

    self%n_weights = n_weights

  end subroutine initialize_particle_group_1d1v



  !----------------------------------------------------------------------!
  !> Constructor for pointer
  subroutine sll_s_new_particle_group_1d1v_ptr(&
       particle_group, &
       n_particles, &
       n_total_particles, &
       charge, &
       mass, &
       n_weights)
    class( sll_c_particle_group_base ),  pointer, intent( out ) :: particle_group !< abstract particle group
    sll_int32                                   , intent( in )  :: n_particles !< number of particles local to the processor
    sll_int32                                   , intent( in )  :: n_total_particles !< number of particles in total simulation
    sll_real64                                  , intent( in )  :: charge !< charge of the particle species
    sll_real64                                  , intent( in )  :: mass   !< mass of the particle species
    sll_int32                                   , intent( in )  :: n_weights !< number of weights
    
    sll_int32                                                   :: ierr

    SLL_ALLOCATE( sll_t_particle_group_1d1v :: particle_group, ierr)

    select type (particle_group)
    type is (sll_t_particle_group_1d1v)
       call particle_group%init( n_particles, n_total_particles, charge, mass, n_weights )
    end select

  end subroutine sll_s_new_particle_group_1d1v_ptr

 !----------------------------------------------------------------------!
  !> Constructor for allocatable
  subroutine sll_s_new_particle_group_1d1v(&
       particle_group, &
       n_particles, &
       n_total_particles, &
       charge, &
       mass, &
       n_weights)
    class( sll_c_particle_group_base ), allocatable,  intent( out ) :: particle_group !< abstract particle group
    sll_int32                                      ,  intent( in )  :: n_particles !< number of particles local to the processor
    sll_int32                                      , intent( in )   :: n_total_particles !< number of particles in total simulation
    sll_real64                                     , intent( in )   :: charge !< charge of the particle species
    sll_real64                                     , intent( in )   :: mass   !< mass of the particle species
    sll_int32                                      , intent( in )   :: n_weights !< number of weights
    

    allocate( sll_t_particle_group_1d1v :: particle_group)

    select type (particle_group)
    type is (sll_t_particle_group_1d1v)
       call particle_group%init( n_particles, n_total_particles, charge, mass, n_weights )
    end select
    
  end subroutine sll_s_new_particle_group_1d1v


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Get position
pure function get_x_1d1v(self,i) result(r)
 class( sll_t_particle_group_1d1v ), intent( in ) :: self !< particle group
    sll_int32                         , intent( in ) :: i !< no. of the particle
    sll_real64 :: r(3) !< position of particle i

    r = 1.0_f64
    r(1) = self%particle_array(1, i)

 end function get_x_1d1v

  !----------------------------------------------------------------------!
  !> Get velocity
  pure function get_v_1d1v( self, i ) result( r )
    class( sll_t_particle_group_1d1v ), intent( in ) :: self !< particle group
    sll_int32                         , intent( in ) :: i !< no. of the particle
    sll_real64 :: r(3)

    r = 1.0_f64
    r(1) = self%particle_array( 2, i )
    
  end function get_v_1d1v

  !----------------------------------------------------------------------!
  !> Get charge of particle (q * particle_weight)
  pure function get_charge_1d1v( self, i , i_weight) result (r)
        class( sll_t_particle_group_1d1v ), intent( in ) :: self !< particle group
    sll_int32                             , intent( in ) :: i !< no. of the particle
    sll_int32, optional                   , intent( in ) :: i_weight !< index of weight to be used (default: 1)
    sll_real64 :: r !< charges(s) of particle i

    sll_int32 :: i_wi

    i_wi = 1
    if(present(i_weight)) i_wi = i_weight
    r = self%species%q  * self%particle_array(2+i_wi, i) * self%common_weight

  end function get_charge_1d1v

!----------------------------------------------------------------------!
  !> Get mass of particle (m * particle_weight)
  pure function get_mass_1d1v( self, i, i_weight) result (r)
        class( sll_t_particle_group_1d1v ), intent( in ) :: self !< particle group
    sll_int32                             , intent( in ) :: i !< no. of the particle
    sll_int32, optional                   , intent( in ) :: i_weight !< index of weight to be used (default: 1)
    sll_real64 :: r !< masses(s) of particle i

    sll_int32 :: i_wi

    i_wi = 1
    if(present(i_weight)) i_wi = i_weight
    r = self%species%m * self%particle_array( 2+i_wi, i) * self%common_weight

  end function get_mass_1d1v

  !----------------------------------------------------------------------!
  !> Get particle weights
  pure function get_weights_1d1v( self, i) result (r)
        class( sll_t_particle_group_1d1v ), intent( in ) :: self !< particle group
    sll_int32                             , intent( in ) :: i !< no. of the particle
    sll_real64                                           :: r(self%n_weights) !< weight(s) of particle i

    r = self%particle_array(3:2+self%n_weights, i)

  end function get_weights_1d1v

 !----------------------------------------------------------------------!
  !> Set position of particle \a i
  subroutine set_x_1d1v( self, i, x )
    class( sll_t_particle_group_1d1v ), intent( inout ) :: self !< particle group
    sll_int32                         , intent( in )    :: i !< no. of the particle
    sll_real64                        , intent( in)     :: x(3) !< first component holds the value of the position to be set

    self%particle_array(1, i) = x(1)
    
  end subroutine set_x_1d1v

  !----------------------------------------------------------------------!
  !> Set velocity of particle \a i
  subroutine set_v_1d1v( self, i, x )
    class( sll_t_particle_group_1d1v ), intent( inout ) :: self !< particle group
    sll_int32                         , intent( in )    :: i !< no. of the particle
    sll_real64                        , intent( in)     :: x(3) !< first component holds the values of the velocity to be set

    self%particle_array(2, i) = x(1)
    
  end subroutine set_v_1d1v
  
  !----------------------------------------------------------------------!
  !> Set weights of particle \a i
  subroutine set_weight_1d1v( self, i, x )
    class( sll_t_particle_group_1d1v ), intent( inout ) :: self !< particle group
    sll_int32                         , intent( in )    :: i !< no. of the particle
    sll_real64                        , intent( in)     :: x(self%n_weights) !< particle weight(s)

    self%particle_array(3:2+self%n_weights, i) = x
    
  end subroutine set_weight_1d1v

  
  !----------------------------------------------------------------------!
  !> Set the common weight
  subroutine set_common_weight_1d1v( self, x )
    class( sll_t_particle_group_1d1v ), intent( inout ) :: self !< particle group
    sll_real64                      , intent( in):: x !< common weight

    self%common_weight = x
    
  end subroutine set_common_weight_1d1v


  !----------------------------------------------------------------------!
  !> Print particle array to file
  subroutine print_particle_group_1d1v(self, filename)
    class( sll_t_particle_group_1d1v ), intent(in) :: self !< particle group
    character(len=*), intent(in) :: filename !< name of output file 
    sll_int32 :: file_id

    open(newunit=file_id,file=filename)
    write(file_id,*) self%particle_array
    close(file_id)
    
  end subroutine print_particle_group_1d1v



end module sll_m_particle_group_1d1v
