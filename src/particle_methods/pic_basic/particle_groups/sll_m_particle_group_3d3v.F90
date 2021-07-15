!> @ingroup particle_groups
!> @author Katharina Kormann, IPP
!> @brief Simple particle group group for 3d3v.
!> @details ...
module sll_m_particle_group_3d3v

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_particle_group_base, only: &
       sll_c_particle_group_base, &
       sll_t_species

  implicit none

  public :: &
       sll_s_new_particle_group_3d3v, &
       sll_s_new_particle_group_3d3v_ptr, &
       sll_t_particle_group_3d3v

  private
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Simple version of a PIC particle group in 2d2v
  type, extends(sll_c_particle_group_base) :: sll_t_particle_group_3d3v
     !sll_int32               :: n_particles !< number of particle
     sll_real64, allocatable :: particle_array(:,:) !< array of particles
     sll_real64 :: common_weight = 1.0_f64

   contains
     ! Getters
     procedure :: get_x  => get_x_3d3v !> Get the values of the coordinate of a particle
     procedure :: get_v  => get_v_3d3v !> Get the values of the velocity of a particle
     procedure :: get_charge => get_charge_3d3v !> Get the charge(s)
     procedure :: get_mass => get_mass_3d3v !> Get the mass(es)
     procedure :: get_weights => get_weights_3d3v !> Get weight(s) of the particle
     procedure :: get_common_weight => get_common_weight_3d3v

     ! Setters
     procedure :: set_x => set_x_3d3v !> Set the values of the coordinate of a particle
     procedure :: set_v => set_v_3d3v !> Set the values of the velocity of a particle
     procedure :: set_weights => set_weights_3d3v !> Set the weight(s) of a particle
     procedure :: set_common_weight => set_common_weight_3d3v !> Set the common weight for the particle

     ! Initializer
     procedure :: init => initialize_particle_group_3d3v  !> Initialization function
     procedure :: free => delete_particle_group_3d3v !> Destructor

     procedure :: print => print_particle_group_3d3v
     procedure :: read  => read_particle_group_3d3v

  end type sll_t_particle_group_3d3v

contains

  !----------------------------------------------------------------------!
  !> Destructor
  subroutine delete_particle_group_3d3v(self)
    class( sll_t_particle_group_3d3v ), intent( inout ) :: self  !< particle group

    deallocate(self%particle_array)

  end subroutine delete_particle_group_3d3v

  !----------------------------------------------------------------------!
  !> Initialize particle group
  subroutine initialize_particle_group_3d3v (&
       self, &
       n_particles, &
       n_total_particles, &
       charge, &
       mass, &
       n_weights)
    class( sll_t_particle_group_3d3v ), intent( inout ) :: self  !< particle group
    sll_int32                         , intent( in )    :: n_particles !< number of particles local to the processor
    sll_int32                         , intent( in )    :: n_total_particles !< number of particles in total simulation
    sll_real64                        , intent( in )    :: charge !< charge of the particle species
    sll_real64                        , intent( in )    :: mass   !< mass of the particle species
    sll_int32                         , intent( in )    :: n_weights !< number of weights

    sll_int32                                           :: ierr

    self%n_particles = n_particles
    self%n_total_particles = n_total_particles

    SLL_ALLOCATE(self%particle_array(6+n_weights, self%n_particles), ierr) 

    allocate(self%species, stat=ierr)
    SLL_ASSERT( ierr == 0)
    call self%species%init( charge, mass)

    self%n_weights = n_weights

  end subroutine initialize_particle_group_3d3v


  !----------------------------------------------------------------------!
  !> Constructor for abstract type
  subroutine sll_s_new_particle_group_3d3v_ptr(&
       particle_group, &
       n_particles, &
       n_total_particles, &
       charge, &
       mass, &
       n_weights)
    class( sll_c_particle_group_base ),  pointer, intent( out )  :: particle_group
    sll_int32                                   , intent( in )   :: n_particles !< number of particles local to the processor
    sll_int32                                   , intent( in )   :: n_total_particles !< number of particles in total simulation
    sll_real64                                  , intent( in )   :: charge !< charge of the particle species
    sll_real64                                  , intent( in )   :: mass   !< mass of the particle species
    sll_int32                                   , intent(in)     :: n_weights !< number of weights

    sll_int32                                                    :: ierr

    SLL_ALLOCATE( sll_t_particle_group_3d3v :: particle_group, ierr)

    select type( particle_group )
    type is ( sll_t_particle_group_3d3v )
       call particle_group%init(n_particles, n_total_particles, charge, mass, n_weights)
    end select

  end subroutine sll_s_new_particle_group_3d3v_ptr


  !----------------------------------------------------------------------!
  !> Constructor for abstract type
  subroutine sll_s_new_particle_group_3d3v(&
       particle_group, &
       n_particles, &
       n_total_particles, &
       charge, &
       mass, &
       n_weights)
    class( sll_c_particle_group_base ),  allocatable, intent( out )   :: particle_group
    sll_int32                                       , intent( in )    :: n_particles !< number of particles local to the processor
    sll_int32                                       , intent( in )    :: n_total_particles !< number of particles in total simulation
    sll_real64                                      , intent( in )    :: charge !< charge of the particle species
    sll_real64                                      , intent( in )    :: mass   !< mass of the particle species
    sll_int32                                       , intent(in)      :: n_weights !< number of weights

    sll_int32                                                         :: ierr

    SLL_ALLOCATE( sll_t_particle_group_3d3v :: particle_group, ierr)

    select type( particle_group )
    type is ( sll_t_particle_group_3d3v )
       call particle_group%init(n_particles, n_total_particles, charge, mass, n_weights)
    end select

  end subroutine sll_s_new_particle_group_3d3v

  !----------------------------------------------------------------------!
  !> Get positions of particle \a i
  pure function get_x_3d3v( self, i ) result( r )
    class( sll_t_particle_group_3d3v ), intent( in ) :: self  !< particle group
    sll_int32                         , intent( in ) :: i !< no. of the particle
    sll_real64 :: r(3) !< first two components hold the value of the particle position

    r = self%particle_array(1:3, i)

  end function get_x_3d3v

  !----------------------------------------------------------------------!
  !> Get velocities of particle \a i
  pure function get_v_3d3v( self, i ) result( r )
    class( sll_t_particle_group_3d3v ), intent( in ) :: self  !< particle group
    sll_int32                         , intent( in ) :: i !< no. of the particle
    sll_real64 :: r(3) !< first two components hold the value of the particle velocity

    r = self%particle_array(4:6, i)

  end function get_v_3d3v


  !----------------------------------------------------------------------!
  !> Get charge of particle \a i ( q * particle_weight)
  pure function get_charge_3d3v( self, i , i_weight) result (r)
    class( sll_t_particle_group_3d3v ), intent( in ) :: self !< particle group
    sll_int32                             , intent( in ) :: i !< no. of the particle
    sll_int32, optional                   , intent( in ) :: i_weight !< index of weight to be used (default: 1)
    sll_real64 :: r !< charges(s) of particle i

    sll_int32 :: i_wi

    i_wi = 1
    if(present(i_weight)) i_wi = i_weight
    r = self%species%q  * self%particle_array(6+i_wi, i) * self%common_weight

  end function get_charge_3d3v


  !----------------------------------------------------------------------!
  !> Get mass of particle \a i ( m * particle_weight)
  pure function get_mass_3d3v( self, i, i_weight) result (r)
    class( sll_t_particle_group_3d3v ), intent( in ) :: self !< particle group
    sll_int32                             , intent( in ) :: i !< no. of the particle
    sll_int32, optional                   , intent( in ) :: i_weight !< index of weight to be used (default: 1)
    sll_real64 :: r !< masses(s) of particle i

    sll_int32 :: i_wi

    i_wi = 1
    if(present(i_weight)) i_wi = i_weight
    r = self%species%m * self%particle_array(6+i_wi, i) * self%common_weight

  end function get_mass_3d3v


  !----------------------------------------------------------------------!
  !> Get weights of particle \a i 
  pure function get_weights_3d3v( self, i) result (r)
    class( sll_t_particle_group_3d3v ), intent( in ) :: self  !< particle group
    sll_int32                             , intent( in ) :: i !< no. of the particle
    sll_real64 :: r(self%n_weights) !< particle mass(es)

    r = self%particle_array(7:6+self%n_weights, i)

  end function get_weights_3d3v

  !----------------------------------------------------------------------!
  !> Set the common weight
  pure function get_common_weight_3d3v( self )result(r)
    class( sll_t_particle_group_3d3v ), intent( in ) :: self !< particle group
    sll_real64                                       :: r !< common weight

    r = self%common_weight

  end function get_common_weight_3d3v

  !----------------------------------------------------------------------!
  !> Set positions of particle \a i 
  subroutine set_x_3d3v( self, i, x )
    class( sll_t_particle_group_3d3v ), intent( inout ) :: self  !< particle group
    sll_int32                       , intent( in )      :: i !< no. of the particle
    sll_real64                      , intent( in)       :: x(3) !< components 1 and 2 hold the particle position to be set

    self%particle_array(1:3, i) = x

  end subroutine set_x_3d3v

  !----------------------------------------------------------------------!
  !> Set velocities of particle \a i 
  subroutine set_v_3d3v( self, i, x )
    class( sll_t_particle_group_3d3v ), intent( inout ) :: self  !< particle group
    sll_int32                         , intent( in )    :: i !< no. of the particle
    sll_real64                        , intent( in)     :: x(3) !< component 1 and 2 hold the particle velocity to be set

    self%particle_array(4:6, i) = x

  end subroutine set_v_3d3v

  !----------------------------------------------------------------------!
  !> Set weights of particle \a i 
  subroutine set_weights_3d3v( self, i, x )
    class( sll_t_particle_group_3d3v ), intent( inout ) :: self  !< particle group
    sll_int32                       , intent( in )      :: i !< no. of the particle
    sll_real64                      , intent( in)       :: x(self%n_weights) !< particle weight(s) to be set

    self%particle_array(7:6+self%n_weights, i) = x

  end subroutine set_weights_3d3v

  !----------------------------------------------------------------------!
  !> Set the common weight
  subroutine set_common_weight_3d3v( self, x )
    class( sll_t_particle_group_3d3v ), intent( inout ) :: self  !< particle group
    sll_real64                      , intent( in)       :: x !< common weight

    self%common_weight = x

  end subroutine set_common_weight_3d3v

  !----------------------------------------------------------------------!
  !> Print particle array
  subroutine print_particle_group_3d3v(self, filename)
    class( sll_t_particle_group_3d3v ), intent(in) :: self
    character(len=*), intent(in) :: filename
    sll_int32 :: file_id
    sll_int32 :: i

    open(newunit=file_id,file=filename//'_weight')
    write(file_id,*) self%common_weight
    close(file_id)

    open(newunit=file_id,file=filename)
    write(file_id,*) self%common_weight
    do i=1, self%n_particles
       write(file_id,*) self%particle_array(:,i)
    end do
    close(file_id)

  end subroutine print_particle_group_3d3v

  !----------------------------------------------------------------------!
  !> Read particle array from file
  subroutine read_particle_group_3d3v(self, filename)
    class( sll_t_particle_group_3d3v ), intent(inout) :: self !< particle group
    character(len=*), intent(in) :: filename !< name of output file
    sll_int32 :: file_id
    sll_int32 :: i, ierr

    open(newunit=file_id,file=filename//'_weight', status='old', iostat=ierr)
    if( ierr .ne. 0 ) then
       SLL_ERROR('read_particle_group_3d3v', 'File not found: '//filename//'_weight')
    end if
    read(file_id,*) self%common_weight
    close(file_id)

    open(newunit=file_id,file=filename, status='old', iostat=ierr)
    if( ierr .ne. 0 ) then
       SLL_ERROR('read_particle_group_3d3v', 'File not found: '//filename)
    end if
    read(file_id,*) self%common_weight
    do i = 1, self%n_particles
       read(file_id,*) self%particle_array(:,i)
       !print*, 'a', i, self%n_particles
    end do
    close(file_id)

  end subroutine read_particle_group_3d3v

end module sll_m_particle_group_3d3v
