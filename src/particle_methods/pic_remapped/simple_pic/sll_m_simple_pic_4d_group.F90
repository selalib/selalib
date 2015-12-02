!**************************************************************
!  Copyright INRIA
!  Authors : 
!     ???
!
!  This code SeLaLib (for Semi-Lagrangian-Library) 
!  is a parallel library for simulating the plasma turbulence 
!  in a tokamak.
!  
!  This software is governed by the CeCILL-B license 
!  under French law and abiding by the rules of distribution 
!  of free software.  You can  use, modify and redistribute 
!  the software under the terms of the CeCILL-B license as 
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info". 
!**************************************************************

!> @ingroup particle_methods

!> @author MCP ALH

!> @brief Module for groups of particles of type sll_simple_pic_4d_particle <!--
!> [[file:simple_pic_4d_particle.F90::sll_simple_pic_4d_particle]] -->

module sll_m_simple_pic_4d_group

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

! #include "particle_representation.h"   NEEDED?

use sll_m_accumulators
  use sll_m_working_precision
  use sll_m_simple_pic_4d_particle
  use sll_m_cartesian_meshes
  use sll_m_remapped_pic_base
  use sll_m_pic_random_initializers, only: sll_pic_4d_random_unweighted_initializer_landau_f0
  use sll_m_remapped_pic_utilities, only:x_is_in_domain_2d, apply_periodic_bc_on_cartesian_mesh_2d
  implicit none

  !> Group of sll_m_simple_pic_4d_particle::sll_simple_pic_4d_particle

  ! [[file:simple_pic_4d_particle.F90::sll_simple_pic_4d_particle]]
  
  type, extends(sll_c_remapped_particle_group) :: sll_simple_pic_4d_group

    !> @name The particles
    !> @{
    ! sll_int32                                                   :: number_particles
    type(sll_simple_pic_4d_particle),   dimension(:), pointer   :: particle_list
    sll_real64                                                  :: common_weight
    !> @}

    !> @name The physical mesh used eg in the Poisson solver
    !> @{
    type(sll_cartesian_mesh_2d), pointer    :: space_mesh_2d
    !> @}

    !> @name The initial density (put this in a separate object ?)
    !> @{
    sll_real64      :: thermal_speed
    sll_real64      :: alpha
    sll_real64      :: k_landau
    !> @}

  contains

    !> @name Getters
    !> @{
    procedure :: get_x          => simple_pic_4d_get_x
    procedure :: get_v          => simple_pic_4d_get_v
    procedure :: get_charge     => simple_pic_4d_get_charge
    procedure :: get_mass       => simple_pic_4d_get_mass
    procedure :: get_cell_index => simple_pic_4d_get_cell_index
    !> @}
    
    !> @name Setters
    !> @{
    procedure :: set_x                  => simple_pic_4d_set_x
    procedure :: set_v                  => simple_pic_4d_set_v
    procedure :: set_common_weight      => simple_pic_4d_set_common_weight
    procedure :: set_particle_weight    => simple_pic_4d_set_particle_weight      ! not to be called for this class
    !> @}
    
    !> @name Initializers
    !> @{
    procedure, pass(self) :: set_landau_parameters  => simple_pic_4d_set_landau_parameters
    procedure :: initializer                        => simple_pic_4d_initializer
    !> @}

    procedure :: deposit_charge_2d          => simple_pic_4d_deposit_charge_2d
    procedure :: remap                      => simple_pic_4d_remap
    procedure :: visualize_f_slice_x_vx     => simple_pic_4d_visualize_f_slice_x_vx


  end type sll_simple_pic_4d_group

  interface sll_delete
     module procedure sll_simple_pic_4d_group_delete
  end interface sll_delete

contains

  !----------------------------------------------------------------------------
  pure function simple_pic_4d_get_common_charge( self ) result( r )
    class( sll_simple_pic_4d_group ), intent( in ) :: self
    sll_int32  :: dummy_i
    sll_real64 :: r

    dummy_i = 0

    ! same charge for all particles
    r = simple_pic_4d_get_charge(self, dummy_i)

  end function simple_pic_4d_get_common_charge

  !----------------------------------------------------------------------------
  pure function simple_pic_4d_get_charge( self, i ) result( r )
    class( sll_simple_pic_4d_group ), intent( in ) :: self
    sll_int32                       , intent( in ) :: i
    sll_real64 :: r
#ifdef DEBUG
    sll_int32 :: dummy
    dummy = i
#endif

    ! same charge for all particles
    r = self%species%q * self%common_weight

  end function simple_pic_4d_get_charge


  !----------------------------------------------------------------------------
  pure function simple_pic_4d_get_mass( self, i ) result( r )
    class( sll_simple_pic_4d_group ), intent( in ) :: self
    sll_int32                       , intent( in ) :: i
    sll_real64 :: r
#ifdef DEBUG
    sll_int32 :: dummy
    dummy = i
#endif

    ! same mass for all particles
    r = self%species%m * self%common_weight

  end function simple_pic_4d_get_mass


  !----------------------------------------------------------------------------
  pure function simple_pic_4d_get_x( self, i ) result( r )
    class( sll_simple_pic_4d_group ), intent( in ) :: self
    sll_int32                       , intent( in ) :: i
    sll_real64 :: r(3)

    ! get x
    r(1) = self%space_mesh_2d%eta1_min + self%space_mesh_2d%delta_eta1*(                            &
                real(self%particle_list(i)%offset_x + self%particle_list(i)%i_cell_x - 1, f64)      &
            )
    ! get y
    r(2) = self%space_mesh_2d%eta2_min + self%space_mesh_2d%delta_eta2*(                            &
                real(self%particle_list(i)%offset_y + self%particle_list(i)%i_cell_y - 1, f64)      &
            )

  end function simple_pic_4d_get_x


  !----------------------------------------------------------------------------
  pure function simple_pic_4d_get_v( self, i ) result( r )
    class( sll_simple_pic_4d_group ), intent( in ) :: self
    sll_int32                       , intent( in ) :: i
    sll_real64 :: r(3)

    ! get vx
    r(1) = self%particle_list(i)%vx
    ! get vy
    r(2) = self%particle_list(i)%vy

  end function simple_pic_4d_get_v


  !----------------------------------------------------------------------------
  ! get the cartesian cell index (here i_out to match the abstract interface), 1 <= i_out <= num_cells_x * num_cells_y
  pure function simple_pic_4d_get_cell_index(self, i) result(i_out)
    class(sll_simple_pic_4d_group),  intent( in )   ::  self
    sll_int32,                      intent( in )    ::  i       !> particle index
    sll_int32                                       ::  i_out   !> cell index
    sll_int32 ::  i_cell_x, i_cell_y
    sll_int32 ::  num_cells_x, num_cells_y

    i_cell_x    = self%particle_list(i)%i_cell_x
    i_cell_y    = self%particle_list(i)%i_cell_y
    num_cells_x = self%space_mesh_2d%num_cells1
    num_cells_y = self%space_mesh_2d%num_cells2

    i_out = 1 + modulo(i_cell_x - 1,  num_cells_x) + modulo(i_cell_y - 1,  num_cells_y) * num_cells_x

  end function simple_pic_4d_get_cell_index



  !----------------------------------------------------------------------------
  ! transforms a standard particle position (x,y) in (i_cell_x, i_cell_y, offset_x, offset_y) and
  ! sets the particle field accordingly.
  ! -> here the indices i_cell_x and i_cell_y do not need to be within [1, mesh%num_cells1] or [1, mesh%num_cells2]
  !    so that: - in periodic domains, the flows are better represented (no information is lost using modulo)
  !             - in non-periodic domains we can track outside particles (markers)
  !
  ! note: the integer index of the physical cell (used eg for the Poisson solver) is then obtained with get_poisson_cell_index
  subroutine simple_pic_4d_set_x( self, i, x )
    class( sll_simple_pic_4d_group ), intent( inout ) :: self
    sll_int32                       , intent( in    ) :: i
    sll_real64                      , intent( in    ) :: x(3)

    type(sll_cartesian_mesh_2d),      pointer :: space_mesh_2d
    type(sll_simple_pic_4d_particle), pointer :: particle
    sll_int32               :: i_cell_x, i_cell_y
    sll_real32              :: offset_x, offset_y
    sll_real64              :: temp

    space_mesh_2d => self%space_mesh_2d
    particle => self%particle_list(i)

    temp = (x(1) - space_mesh_2d%eta1_min) / space_mesh_2d%delta_eta1
    i_cell_x  = 1 + int(floor(temp))
    offset_x = real(temp - real(i_cell_x - 1,f64), f32)

    temp = (x(2) - space_mesh_2d%eta2_min) / space_mesh_2d%delta_eta2
    i_cell_y  = 1 + int(floor(temp))
    offset_y = real(temp - real(i_cell_y - 1,f64), f32)

    SLL_ASSERT(offset_x >= 0)
    SLL_ASSERT(offset_x <= 1 )
    SLL_ASSERT(offset_y >= 0)
    SLL_ASSERT(offset_y <= 1 )

    particle%i_cell_x = i_cell_x
    particle%i_cell_y = i_cell_y
    particle%offset_x = offset_x
    particle%offset_y = offset_y

  end subroutine simple_pic_4d_set_x


  !----------------------------------------------------------------------------
  subroutine simple_pic_4d_set_v( self, i, x )
    class( sll_simple_pic_4d_group ), intent( inout ) :: self
    sll_int32                       , intent( in    ) :: i
    sll_real64                      , intent( in    ) :: x(3)  !> this is the velocity, but argument name in abstract interface is x

    type(sll_simple_pic_4d_particle), pointer :: particle

    particle => self%particle_list(i)
    particle%vx = x(1)
    particle%vy = x(2)

  end subroutine simple_pic_4d_set_v


  !----------------------------------------------------------------------------
  subroutine simple_pic_4d_set_common_weight( self, s )
    class( sll_simple_pic_4d_group ), intent( inout ) :: self
    sll_real64                      , intent( in    ) :: s

    self%common_weight = s

  end subroutine simple_pic_4d_set_common_weight


  !----------------------------------------------------------------------------
  subroutine simple_pic_4d_set_particle_weight( self, i, s )
    class( sll_simple_pic_4d_group ), intent( inout ) :: self
    sll_int32                       , intent( in    ) :: i
    sll_real64                      , intent( in    ) :: s

    print*, "Error (8654354237645) -- this subroutine is not implemented for simple_pic_4d_group objects"
    print*, i, s, storage_size(self)
    stop

  end subroutine simple_pic_4d_set_particle_weight


  !----------------------------------------------------------------------------
  subroutine simple_pic_4d_set_landau_parameters( self, thermal_speed, alpha, k_landau )
    class( sll_simple_pic_4d_group ), intent( inout ) :: self
    sll_real64                      , intent( in    ) :: thermal_speed
    sll_real64                      , intent( in    ) :: alpha
    sll_real64                      , intent( in    ) :: k_landau

    self%thermal_speed = thermal_speed
    self%alpha = alpha
    self%k_landau = k_landau

  end subroutine simple_pic_4d_set_landau_parameters


  !----------------------------------------------------------------------------
   subroutine simple_pic_4d_initializer( self, initial_density_identifier, rand_seed, rank, world_size )
    class( sll_simple_pic_4d_group ), intent( inout ) :: self
    sll_int32                       , intent( in    ) :: initial_density_identifier
    sll_int32, dimension(:)         , intent( in ), optional :: rand_seed
    sll_int32                       , intent( in ), optional :: rank, world_size
    !sll_int32                       :: ierr

    ! for the moment we only use the landau damping initial density
    ! so we don't use the initial_density_identifier, but eventually it should say which initial density is used

    call sll_pic_4d_random_unweighted_initializer_landau_f0 (   &
      self%thermal_speed, self%alpha, self%k_landau,            & ! -> these parameters should be members of the initializer object
      self,                                                     &
      self%space_mesh_2d,                                       &
      self%number_particles,                                    &
      rand_seed, rank, world_size                               &
    )
   return
   !PN ADD TO PREVENT WARNING
   print*, initial_density_identifier

   end subroutine simple_pic_4d_initializer


  subroutine simple_pic_4d_deposit_charge_2d( self, charge_accumulator, target_total_charge )
    class( sll_simple_pic_4d_group ),           intent( inout )  :: self
    type( sll_charge_accumulator_2d ), pointer, intent( inout ) :: charge_accumulator
    sll_real64,                                 intent(in), optional :: target_total_charge             ! for a check (if present)

    type( charge_accumulator_cell_2d ), pointer             :: charge_accumulator_cell
    type(sll_simple_pic_4d_particle), pointer :: particle
    sll_int32       :: i_part, i_cell
    sll_real64      :: xy_part(3)
    sll_real64      :: particle_charge
    sll_real64      :: deposited_charge
    sll_real64      :: dx, dy

    call reset_charge_accumulator_to_zero ( charge_accumulator )

    deposited_charge = 0.0_f64

    particle_charge = simple_pic_4d_get_common_charge(self)
    do i_part = 1, self%number_particles

      particle => self%particle_list( i_part )
      dx = real(particle%offset_x,f64)
      dy = real(particle%offset_y,f64)
      i_cell = self%get_cell_index(i_part)
      charge_accumulator_cell => charge_accumulator%q_acc(i_cell)

      xy_part = self%get_x( i_part )  ! x and y
      if( x_is_in_domain_2d(    xy_part(1), xy_part(2),         &
                                self%space_mesh_2d,             &
                                self%domain_is_periodic(1),   &
                                self%domain_is_periodic(2) ))then

        charge_accumulator_cell%q_sw = charge_accumulator_cell%q_sw + particle_charge * (1.0_f64 - dx) * (1.0_f64 - dy)
        charge_accumulator_cell%q_se = charge_accumulator_cell%q_se + particle_charge *            dx  * (1.0_f64 - dy)
        charge_accumulator_cell%q_nw = charge_accumulator_cell%q_nw + particle_charge * (1.0_f64 - dx) *            dy
        charge_accumulator_cell%q_ne = charge_accumulator_cell%q_ne + particle_charge *            dx  *            dy
        deposited_charge = deposited_charge + particle_charge

      else ! particle not in domain (should store the reference for later processing)
        print*, "Error (097647687): for the moment every particle should be in the (periodic) 2d domain..."
        stop
      end if

      if( present(target_total_charge) )then
        if( abs( deposited_charge - target_total_charge ) > 0.0000001 * abs(target_total_charge) )then
           print*, "Warning (8756537654) in [simple_pic_4d_deposit_charge_2d]: deposited_charge and target_total_charge differ"
           print*, "Warning (8756537654) deposited_charge    = ", deposited_charge
           print*, "Warning (8756537654) target_total_charge = ", target_total_charge
        end if
      end if

    end do

  end subroutine simple_pic_4d_deposit_charge_2d

  !----------------------------------------------------------------------------
  subroutine simple_pic_4d_remap(self)

    class( sll_simple_pic_4d_group ),   intent( inout ) :: self

    print*, " ------------ ------------ ------------ ------------ ------------ ------------ ------------ ------------ ------------"
    print*, " WARNING (765764768675) -- remap routine called for a group of simple_pic_4d particles has no effect...              "
    print*, " ------------ ------------ ------------ ------------ ------------ ------------ ------------ ------------ ------------"
    print*, storage_size(self)

  end subroutine simple_pic_4d_remap

  !----------------------------------------------------------------------------
  subroutine simple_pic_4d_visualize_f_slice_x_vx(self, array_name, iplot)

    class( sll_simple_pic_4d_group ),   intent( inout ) :: self
    character(len=*),                   intent(in)      :: array_name !< field name
    sll_int32,                          intent(in)      :: iplot      !< plot counter
    !character(len=4)                                    :: cplot

    print*, " ------------ ------------ ------------ ------------ ------------ ------------ ------------ ------------ ------------"
    print*, " WARNING (956542375763) -- this function (simple_pic_4d_visualize_f_slice_x_vx) does nothing, need to be implemented!"
    print*, " ------------ ------------ ------------ ------------ ------------ ------------ ------------ ------------ ------------"
    print*, len(array_name), iplot, storage_size(self)

  end subroutine simple_pic_4d_visualize_f_slice_x_vx

  !----------------------------------------------------------------------------
  ! Constructor
  function sll_simple_pic_4d_group_new( &
        species_charge,         &
        species_mass,           &
        particle_group_id,      &
        domain_is_x_periodic,   &
        domain_is_y_periodic,   &
        number_particles,       &
        space_mesh_2d ) result(res)

    type( sll_simple_pic_4d_group ), pointer :: res

    sll_int32, intent(in)   :: particle_group_id
    sll_real64, intent(in)  :: species_charge
    sll_real64, intent(in)  :: species_mass
    logical, intent(in)     :: domain_is_x_periodic
    logical, intent(in)     :: domain_is_y_periodic
    sll_int32, intent( in ) :: number_particles

    type(sll_cartesian_mesh_2d), pointer, intent(in) :: space_mesh_2d
    sll_int32               :: ierr

    SLL_ALLOCATE( res, ierr )

    ! the group
    res%species => temp_species_new( species_charge, species_mass )

    res%id = particle_group_id
    res%dimension_x = 2
    res%dimension_v = 2

    res%number_particles = number_particles
    SLL_ALLOCATE( res%particle_list(number_particles), ierr )

    ! physical 2d mesh, used eg in the Poisson solver
    if (.not.associated(space_mesh_2d) ) then
       print*, 'ERROR (876876456), given space_mesh_2d is not associated'
    end if
    res%space_mesh_2d => space_mesh_2d
    res%domain_is_periodic(1) = domain_is_x_periodic
    res%domain_is_periodic(2) = domain_is_y_periodic

    end function sll_simple_pic_4d_group_new


  !----------------------------------------------------------------------------
  ! Destructor
  subroutine sll_simple_pic_4d_group_delete(particle_group)
    class(sll_simple_pic_4d_group), pointer :: particle_group
    sll_int32 :: ierr

    if(.not. associated(particle_group) ) then
       print *, 'sll_simple_pic_4d_group_delete(): ERROR (986876), passed group was not associated.'
    end if
    SLL_DEALLOCATE(particle_group%particle_list, ierr)
    SLL_DEALLOCATE(particle_group%space_mesh_2d, ierr)
    SLL_DEALLOCATE(particle_group, ierr)

  end subroutine sll_simple_pic_4d_group_delete


end module sll_m_simple_pic_4d_group
