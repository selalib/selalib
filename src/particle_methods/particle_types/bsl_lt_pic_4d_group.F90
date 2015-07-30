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

!> @brief Module for groups of particles of type sll_bsl_lt_pic_4d_particle <!--
!> [[file:bsl_lt_pic_4d_particle.F90::sll_bsl_lt_pic_4d_particle]] -->

module sll_bsl_lt_pic_4d_group_module

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_accumulators.h"
#include "sll_errors.h"

! #include "particle_representation.h"   NEEDED?

  use sll_constants, only: sll_pi
  use sll_working_precision
  use sll_cartesian_meshes
  use sll_module_pic_base
  use sll_bsl_lt_pic_4d_particle_module
  use sll_bsl_lt_pic_4d_utilities_module

  implicit none

  !> Group of @ref sll_bsl_lt_pic_4d_particle
  type, extends(sll_particle_group_base) :: sll_bsl_lt_pic_4d_group

    !> @name The particles
    !> @{
    sll_int32                                                   :: spline_degree
    sll_int32                                                   :: number_parts_x
    sll_int32                                                   :: number_parts_y
    sll_int32                                                   :: number_parts_vx
    sll_int32                                                   :: number_parts_vy
    sll_int32                                                   :: number_particles
    type(sll_bsl_lt_pic_4d_particle),   dimension(:), pointer   :: particle_list
    !> @}

    !> @name The physical mesh used eg in the Poisson solver
    !> @{
    type(sll_cartesian_mesh_2d), pointer    :: space_mesh_2d
    !> @}

    !> @name The remapping grid in phase space and quasi-interpolation coefficients (for cubic spline particle shapes)
    !> @{
    type(sll_cartesian_mesh_4d),                pointer         :: remapping_grid
    sll_real64, dimension(:,:,:,:),             pointer         :: target_values
    sll_real64, dimension(:),                   pointer         :: lt_pic_interpolation_coefs
    !> @}

    !> @name The initial density (at some point this should be put in a separate initializer object)
    !> @{
    sll_real64      :: thermal_speed
    sll_real64      :: alpha
    sll_real64      :: k_landau
    !> @}

  contains

    !> @name Getters
    !> @{
    procedure :: get_x          => bsl_lt_pic_4d_get_x
    procedure :: get_v          => bsl_lt_pic_4d_get_v
    procedure :: get_charge     => bsl_lt_pic_4d_get_charge
    procedure :: get_mass       => bsl_lt_pic_4d_get_mass
    procedure :: get_cell_index => bsl_lt_pic_4d_get_cell_index
    !> @}
    
    !> @name Setters
    !> @{
    procedure :: set_x                  => bsl_lt_pic_4d_set_x
    procedure :: set_v                  => bsl_lt_pic_4d_set_v

    ! todo: use only one function with particle index as optional parameter
    procedure :: set_common_weight      => bsl_lt_pic_4d_set_common_weight     ! not to be called for this class
    procedure :: set_particle_weight    => bsl_lt_pic_4d_set_particle_weight
    !> @}
    
    !> @name Initializers
    !> @{
    procedure, pass(self) :: set_landau_parameters  => bsl_lt_pic_4d_set_landau_parameters
    procedure             :: initializer     => bsl_lt_pic_4d_initializer
    !> @}
    
    procedure :: deposit_charge_2d          => bsl_lt_pic_4d_deposit_charge_2d

    procedure :: bsl_lt_pic_4d_initializer_landau_f0
    procedure :: bsl_lt_pic_4d_write_landau_density_on_remap_grid
    procedure :: bsl_lt_pic_4d_compute_new_particles

  end type sll_bsl_lt_pic_4d_group

  interface sll_delete
     module procedure sll_bsl_lt_pic_4d_group_delete
  end interface sll_delete

  !! MCP (July 16) -- this is to make the subroutine external, in a separate file, but does not work yet --

!  interface
!    subroutine bsl_lt_pic_4d_initializer( self, initial_density_identifier, rand_seed, rank, world_size )
!    use sll_working_precision
!    import sll_bsl_lt_pic_4d_group
!    class( sll_bsl_lt_pic_4d_group ), intent( inout ) :: self
!    sll_int32                       , intent( in    ) :: initial_density_identifier
!    sll_int32, dimension(:)         , intent( in ), optional :: rand_seed
!    sll_int32                       , intent( in ), optional :: rank, world_size
!    !    sll_int32                       :: ierr
!    !
!    !    call self%initializer_landau_f0 (              &
!    !      self%thermal_speed, self%alpha, self%k_landau )        ! -> these parameters should be members of the initializer object
!
!    end subroutine bsl_lt_pic_4d_initializer
!  end interface

contains

  !----------------------------------------------------------------------------
  pure function bsl_lt_pic_4d_get_charge( self, i ) result( r )
    class( sll_bsl_lt_pic_4d_group ), intent( in ) :: self
    sll_int32                       , intent( in ) :: i
    sll_real64 :: r

    r = self%species%q * self%particle_list(i)%weight

  end function bsl_lt_pic_4d_get_charge


  !----------------------------------------------------------------------------
  pure function bsl_lt_pic_4d_get_mass( self, i ) result( r )
    class( sll_bsl_lt_pic_4d_group ), intent( in ) :: self
    sll_int32                       , intent( in ) :: i
    sll_real64 :: r

    r = self%species%m * self%particle_list(i)%weight

  end function bsl_lt_pic_4d_get_mass


  !----------------------------------------------------------------------------
  pure function bsl_lt_pic_4d_get_x( self, i ) result( r )
    class( sll_bsl_lt_pic_4d_group ), intent( in ) :: self
    sll_int32                       , intent( in ) :: i
    sll_real64 :: r(3)

    ! get x
    r(1) = self%space_mesh_2d%eta1_min + self%space_mesh_2d%delta_eta1*(                            &
                self%particle_list(i)%offset_x + real(self%particle_list(i)%i_cell_x - 1, f64)      &
            )
    ! get y
    r(2) = self%space_mesh_2d%eta2_min + self%space_mesh_2d%delta_eta2*(                            &
                self%particle_list(i)%offset_y + real(self%particle_list(i)%i_cell_y - 1, f64)      &
            )

  end function bsl_lt_pic_4d_get_x


  !----------------------------------------------------------------------------
  pure function bsl_lt_pic_4d_get_v( self, i ) result( r )
    class( sll_bsl_lt_pic_4d_group ), intent( in ) :: self
    sll_int32                       , intent( in ) :: i
    sll_real64 :: r(3)

    ! get vx
    r(1) = self%particle_list(i)%vx
    ! get vy
    r(2) = self%particle_list(i)%vy

  end function bsl_lt_pic_4d_get_v


  !----------------------------------------------------------------------------
  ! get the cartesian cell index (here i_out to match the abstract interface), 1 <= i_out <= num_cells_x * num_cells_y
  !
  ! same function in the simple_pic_4d group -- todo: possible to use the same function?
  pure function bsl_lt_pic_4d_get_cell_index(self, i) result(i_out)
    class(sll_bsl_lt_pic_4d_group),  intent( in )   ::  self
    sll_int32,                      intent( in )    ::  i       !> particle index
    sll_int32                                       ::  i_out   !> cell index
    sll_int32 ::  i_cell_x, i_cell_y
    sll_int32 ::  num_cells_x, num_cells_y

    i_cell_x    = self%particle_list(i)%i_cell_x
    i_cell_y    = self%particle_list(i)%i_cell_y
    num_cells_x = self%space_mesh_2d%num_cells1
    num_cells_y = self%space_mesh_2d%num_cells2

    i_out = 1 + modulo(i_cell_x - 1,  num_cells_x) + modulo(i_cell_y - 1,  num_cells_y) * num_cells_x

  end function bsl_lt_pic_4d_get_cell_index



  !----------------------------------------------------------------------------
  ! transforms a standard particle position (x,y) in (i_cell_x, i_cell_y, offset_x, offset_y) and
  ! sets the particle field accordingly.
  ! -> here the indices i_cell_x and i_cell_y do not need to be within [1, mesh%num_cells1] or [1, mesh%num_cells2]
  !    so that: - in periodic domains, the flows are better represented (no information is lost using modulo)
  !             - in non-periodic domains we can track outside particles (markers)
  !
  ! note: the integer index of the physical cell (used eg for the Poisson solver) is then obtained with get_poisson_cell_index
  !
  ! same function in the simple_pic_4d group -- todo: possible to use the same function?
  subroutine bsl_lt_pic_4d_set_x( self, i, x )
    class( sll_bsl_lt_pic_4d_group ), intent( inout ) :: self
    sll_int32                       , intent( in    ) :: i
    sll_real64                      , intent( in    ) :: x(3)

    type(sll_cartesian_mesh_2d),      pointer :: space_mesh_2d
    type(sll_bsl_lt_pic_4d_particle), pointer :: particle
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

  end subroutine bsl_lt_pic_4d_set_x


  !----------------------------------------------------------------------------
  subroutine bsl_lt_pic_4d_set_v( self, i, x )
    class( sll_bsl_lt_pic_4d_group ), intent( inout ) :: self
    sll_int32                       , intent( in    ) :: i
    sll_real64                      , intent( in    ) :: x(3)  !> this is the velocity, but argument name in abstract interface is x

    type(sll_bsl_lt_pic_4d_particle), pointer :: particle

    particle => self%particle_list(i)
    particle%vx = x(1)
    particle%vy = x(2)

  end subroutine bsl_lt_pic_4d_set_v


  !----------------------------------------------------------------------------
  subroutine bsl_lt_pic_4d_set_common_weight( self, s )
    class( sll_bsl_lt_pic_4d_group ), intent( inout ) :: self
    sll_real64                      , intent( in    ) :: s

    print*, "Error (9O8657864) -- this subroutine is not implemented for bsl_lt_pic_4d_group objects"
    stop

  end subroutine bsl_lt_pic_4d_set_common_weight


  !----------------------------------------------------------------------------
  subroutine bsl_lt_pic_4d_set_particle_weight( self, i, s )
    class( sll_bsl_lt_pic_4d_group ), intent( inout ) :: self
    sll_int32                       , intent( in    ) :: i
    sll_real64                      , intent( in    ) :: s

    self%particle_list(i)%weight = s

  end subroutine bsl_lt_pic_4d_set_particle_weight


  !----------------------------------------------------------------------------
  subroutine bsl_lt_pic_4d_set_landau_parameters( self, thermal_speed, alpha, k_landau )
    class( sll_bsl_lt_pic_4d_group ), intent( inout ) :: self
    sll_real64                      , intent( in    ) :: thermal_speed
    sll_real64                      , intent( in    ) :: alpha
    sll_real64                      , intent( in    ) :: k_landau

    self%thermal_speed = thermal_speed
    self%alpha = alpha
    self%k_landau = k_landau

  end subroutine bsl_lt_pic_4d_set_landau_parameters


  !----------------------------------------------------------------------------
    !   subroutine bsl_lt_pic_4d_initializer( self, initial_density_identifier, rand_seed, rank, world_size )
    !    class( sll_bsl_lt_pic_4d_group ), intent( inout ) :: self
    !    sll_int32                       , intent( in    ) :: initial_density_identifier
    !    sll_int32, dimension(:)         , intent( in ), optional :: rand_seed
    !    sll_int32                       , intent( in ), optional :: rank, world_size
    !    sll_int32                       :: ierr
    !
    !    call self%initializer_landau_f0 (              &
    !      self%thermal_speed, self%alpha, self%k_landau )        ! -> these parameters should be members of the initializer object
    !
    !   end subroutine bsl_lt_pic_4d_initializer

  subroutine bsl_lt_pic_4d_deposit_charge_2d( self, charge_accumulator )
    class( sll_bsl_lt_pic_4d_group ),           intent( inout )  :: self
    type( sll_charge_accumulator_2d ), pointer, intent( inout ) :: charge_accumulator

    type( charge_accumulator_cell_2d ), pointer             :: charge_accumulator_cell
    type(sll_bsl_lt_pic_4d_particle), pointer :: particle
    sll_int32       :: i_part, i_cell
    sll_real64      :: xy_part(3)
    sll_real64      :: particle_charge
    sll_real64      :: dx, dy

    call reset_charge_accumulator_to_zero ( charge_accumulator )


    print *, "error (8768764876598759764) please implement this"
    stop

    do i_part = 1, self%number_particles

      particle => self%particle_list( i_part )
      dx = particle%offset_x
      dy = particle%offset_y
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

      else ! particle not in domain (should store the reference for later processing)
        print*, "Error (097647687): for the moment every particle should be in the (periodic) 2d domain..."
        stop
      end if

    end do

  end subroutine bsl_lt_pic_4d_deposit_charge_2d

  !----------------------------------------------------------------------------
  ! Constructor
  !> @brief Constructor for a group of bsl_lt_pic_4d particles
  function sll_bsl_lt_pic_4d_group_new( &
        species_charge,         &
        species_mass,           &
        particle_group_id,      &
        domain_is_x_periodic,   &
        domain_is_y_periodic,   &
        spline_degree,          &
        number_parts_x,         &
        number_parts_y,         &
        number_parts_vx,        &
        number_parts_vy,        &
        remap_grid_vx_min,      &
        remap_grid_vx_max,      &
        remap_grid_vy_min,      &
        remap_grid_vy_max,      &
        space_mesh_2d ) result(res)

    type( sll_bsl_lt_pic_4d_group ), pointer :: res

    sll_real64, intent(in)  :: species_charge
    sll_real64, intent(in)  :: species_mass
    sll_int32, intent(in)   :: particle_group_id
    logical, intent(in)     :: domain_is_x_periodic
    logical, intent(in)     :: domain_is_y_periodic
    sll_int32, intent(in)   :: spline_degree
    sll_int32, intent(in)   :: number_parts_x
    sll_int32, intent(in)   :: number_parts_y
    sll_int32, intent(in)   :: number_parts_vx
    sll_int32, intent(in)   :: number_parts_vy
    sll_real64, intent(in)  :: remap_grid_vx_min
    sll_real64, intent(in)  :: remap_grid_vx_max
    sll_real64, intent(in)  :: remap_grid_vy_min
    sll_real64, intent(in)  :: remap_grid_vy_max
    type(sll_cartesian_mesh_2d), pointer, intent(in) :: space_mesh_2d

    sll_int32               :: remap_grid_number_cells_x
    sll_int32               :: remap_grid_number_cells_y
    sll_int32               :: remap_grid_number_cells_vx
    sll_int32               :: remap_grid_number_cells_vy
    sll_real64              :: remap_grid_x_min
    sll_real64              :: remap_grid_x_max
    sll_real64              :: remap_grid_y_min
    sll_real64              :: remap_grid_y_max

    sll_int32               :: ierr
    character(len=*), parameter :: this_fun_name = "sll_bsl_lt_pic_4d_group_new"
    character(len=128)       :: err_msg

    SLL_ALLOCATE( res, ierr )

    !> create the species object for this particle group
    res%species => species_new( species_charge, species_mass )

    res%id = particle_group_id
    res%dimension_x = 2
    res%dimension_v = 2

    !> create the particle list
    res%spline_degree = spline_degree
    res%number_parts_x   = number_parts_x
    res%number_parts_y   = number_parts_y
    res%number_parts_vx  = number_parts_vx
    res%number_parts_vy  = number_parts_vy
    res%number_particles = number_parts_x * number_parts_y * number_parts_vx * number_parts_vy

    SLL_ALLOCATE( res%particle_list(res%number_particles), ierr )

    !> assign the physical 2d mesh (used eg in the Poisson solver)
    if (.not.associated(space_mesh_2d) ) then
       err_msg = 'Error: given space_mesh_2d is not associated'
       SLL_ERROR( this_fun_name, err_msg )
    end if
    res%space_mesh_2d => space_mesh_2d
    res%domain_is_periodic(1) = domain_is_x_periodic
    res%domain_is_periodic(2) = domain_is_y_periodic

    !> create the particle list and the array of target values
    SLL_ALLOCATE( res%target_values(number_parts_x, number_parts_y, number_parts_vx, number_parts_vy), ierr )

    !> create the ''remapping grid'', a cartesian phase space grid used to initialize and remap the particles
    remap_grid_x_min = space_mesh_2d%eta1_min
    remap_grid_x_max = space_mesh_2d%eta1_max
    remap_grid_y_min = space_mesh_2d%eta2_min
    remap_grid_y_max = space_mesh_2d%eta2_max

    if( domain_is_x_periodic )then
        remap_grid_number_cells_x = number_parts_x
    else
        remap_grid_number_cells_x = number_parts_x - 1
    end if
    if( domain_is_y_periodic )then
        remap_grid_number_cells_y = number_parts_y
    else
        remap_grid_number_cells_y = number_parts_y - 1
    end if
    remap_grid_number_cells_vx = number_parts_vx - 1
    remap_grid_number_cells_vy = number_parts_vy - 1
    res%remapping_grid => new_cartesian_mesh_4d(remap_grid_number_cells_x,        &
                                                remap_grid_number_cells_y,        &
                                                remap_grid_number_cells_vx,       &
                                                remap_grid_number_cells_vy,       &
                                                remap_grid_x_min,   &
                                                remap_grid_x_max,   &
                                                remap_grid_y_min,   &
                                                remap_grid_y_max,   &
                                                remap_grid_vx_min,  &
                                                remap_grid_vx_max,  &
                                                remap_grid_vy_min,  &
                                                remap_grid_vy_max   &
                                              )

    !> store the quasi-interpolation coefficients used in the remappings and initialisation, in the case of cubic spline particles
    if( spline_degree == 3) then
        SLL_ALLOCATE( res%lt_pic_interpolation_coefs(-1:1), ierr )
           res%lt_pic_interpolation_coefs(-1) = -1./6.
           res%lt_pic_interpolation_coefs(0)  =  8./6.
           res%lt_pic_interpolation_coefs(1)  =  -1./6.
    end if

    end function sll_bsl_lt_pic_4d_group_new


  subroutine bsl_lt_pic_4d_initializer( self, initial_density_identifier, rand_seed, rank, world_size )

    class( sll_bsl_lt_pic_4d_group ), intent( inout ) :: self
    sll_int32                       , intent( in    ) :: initial_density_identifier
    sll_int32, dimension(:)         , intent( in ), optional :: rand_seed
    sll_int32                       , intent( in ), optional :: rank, world_size
    sll_int32                       :: ierr

    call self%bsl_lt_pic_4d_initializer_landau_f0 (     &
      self%thermal_speed, self%alpha, self%k_landau     &
    )        ! -> these parameters should be members of the initializer object

  end subroutine bsl_lt_pic_4d_initializer


  ! initialize the bsl_lt_pic group with the landau f0 distribution
  subroutine bsl_lt_pic_4d_initializer_landau_f0 (          &
              p_group, thermal_speed, alpha, k_landau       &
          )
    class(sll_bsl_lt_pic_4d_group), intent(inout)       :: p_group
    sll_real64, intent(in)                              :: thermal_speed, alpha, k_landau


    call p_group%bsl_lt_pic_4d_write_landau_density_on_remap_grid( thermal_speed, alpha, k_landau )
    call p_group%bsl_lt_pic_4d_compute_new_particles()

  end subroutine bsl_lt_pic_4d_initializer_landau_f0

  subroutine bsl_lt_pic_4d_write_landau_density_on_remap_grid(    &
              p_group,                              &
              thermal_speed, alpha, k_landau        &
              )

    class(sll_bsl_lt_pic_4d_group), intent(inout)    :: p_group
    sll_real64, intent(in)                          :: thermal_speed, alpha, k_landau

    sll_int32 :: j_x
    sll_int32 :: j_y
    sll_int32 :: j_vx
    sll_int32 :: j_vy
    sll_int32 :: number_particles
    sll_int32 :: number_parts_x
    sll_int32 :: number_parts_y
    sll_int32 :: number_parts_vx
    sll_int32 :: number_parts_vy
    sll_real64 :: h_parts_x
    sll_real64 :: h_parts_y
    sll_real64 :: h_parts_vx
    sll_real64 :: h_parts_vy
    sll_real64 :: parts_x_min
    sll_real64 :: parts_y_min
    sll_real64 :: parts_vx_min
    sll_real64 :: parts_vy_min
    sll_real64 :: one_over_thermal_velocity
    sll_real64 :: one_over_two_pi
    sll_real64 :: x_j
    sll_real64 :: y_j
    sll_real64 :: vx_j
    sll_real64 :: vy_j
    type(sll_cartesian_mesh_2d),      pointer  :: m2d
    sll_real64 :: f_x, f_y, f_vx, f_vy

    number_particles = p_group%number_particles
    one_over_thermal_velocity = 1./thermal_speed
    one_over_two_pi = 1./(2*sll_pi)

    number_parts_x  = p_group%number_parts_x
    number_parts_y  = p_group%number_parts_y
    number_parts_vx = p_group%number_parts_vx
    number_parts_vy = p_group%number_parts_vy

    h_parts_x    = p_group%remapping_grid%delta_eta1
    h_parts_y    = p_group%remapping_grid%delta_eta2
    h_parts_vx   = p_group%remapping_grid%delta_eta3
    h_parts_vy   = p_group%remapping_grid%delta_eta4

    parts_x_min    = p_group%remapping_grid%eta1_min
    parts_y_min    = p_group%remapping_grid%eta2_min
    parts_vx_min   = p_group%remapping_grid%eta3_min
    parts_vy_min   = p_group%remapping_grid%eta4_min

    ! Poisson mesh associated to the particles
    m2d => p_group%space_mesh_2d

    ! compute the values of f0 on the (cartesian, phase-space) remapping grid
    x_j = parts_x_min
    do j_x = 1, number_parts_x
      f_x = eval_landau_fx(alpha, k_landau, x_j)
      y_j = parts_y_min
      do j_y = 1, number_parts_y
        vx_j = parts_vx_min
        do j_vx = 1, number_parts_vx
          f_vx = one_over_thermal_velocity * exp(-0.5*(vx_j*one_over_thermal_velocity)**2)
          vy_j = parts_vy_min
          do j_vy = 1, number_parts_vy
            f_vy = one_over_thermal_velocity * exp(-0.5*(vy_j*one_over_thermal_velocity)**2)
            p_group%target_values(j_x,j_y,j_vx,j_vy) = one_over_two_pi * f_x * f_vx * f_vy
            vy_j = vy_j + h_parts_vy
          end do
          vx_j = vx_j + h_parts_vx
        end do
        y_j = y_j + h_parts_y
      end do
      x_j = x_j + h_parts_x
    end do

  end subroutine bsl_lt_pic_4d_write_landau_density_on_remap_grid


  function eval_landau_fx(alpha, kx, x)
    sll_real64 :: alpha, kx, x
    sll_real64 :: eval_landau_fx
    eval_landau_fx = 1._f64 + alpha * cos(kx * x)
  end function eval_landau_fx


  ! position the particle on the cartesian remapping grid
  ! and compute new weights in order to approximate the point values stored in the remapping grid
  !  subroutine sll_compute_new_lt_particles_4d( &        ! old name
  subroutine bsl_lt_pic_4d_compute_new_particles( &
              p_group )

    class(sll_bsl_lt_pic_4d_group), intent(inout) :: p_group
    sll_int32 :: k, k_ngb
    sll_int32 :: k_temp_debug
    sll_int32 :: j_x
    sll_int32 :: j_y
    sll_int32 :: j_vx
    sll_int32 :: j_vy
    sll_int32 :: j_aux_x
    sll_int32 :: j_aux_y
    sll_int32 :: j_aux_vx
    sll_int32 :: j_aux_vy
    sll_int32 :: l_x
    sll_int32 :: l_y
    sll_int32 :: l_vx
    sll_int32 :: l_vy
    sll_int32 :: ierr
    sll_int32 :: number_parts_x
    sll_int32 :: number_parts_y
    sll_int32 :: number_parts_vx
    sll_int32 :: number_parts_vy
    sll_real64 :: h_parts_x
    sll_real64 :: h_parts_y
    sll_real64 :: h_parts_vx
    sll_real64 :: h_parts_vy
    sll_real64 :: parts_x_min
    sll_real64 :: parts_y_min
    sll_real64 :: parts_vx_min
    sll_real64 :: parts_vy_min
    sll_real64 :: w_k
    sll_real64 :: x_j
    sll_real64 :: y_j
    sll_real64 :: vx_j
    sll_real64 :: vy_j
    sll_real32 :: d_vol
    type(sll_cartesian_mesh_2d),    pointer :: m2d
    sll_int32,  dimension(:,:,:,:), pointer :: particle_indices    !  why pointer ?
    sll_real64, dimension(3)                :: coords

    number_parts_x  = p_group%number_parts_x
    number_parts_y  = p_group%number_parts_y
    number_parts_vx = p_group%number_parts_vx
    number_parts_vy = p_group%number_parts_vy

    h_parts_x    = p_group%remapping_grid%delta_eta1
    h_parts_y    = p_group%remapping_grid%delta_eta2
    h_parts_vx   = p_group%remapping_grid%delta_eta3
    h_parts_vy   = p_group%remapping_grid%delta_eta4

    d_vol = real( h_parts_x*h_parts_y*h_parts_vx*h_parts_vy,f32)

    parts_x_min    = p_group%remapping_grid%eta1_min
    parts_y_min    = p_group%remapping_grid%eta2_min
    parts_vx_min   = p_group%remapping_grid%eta3_min
    parts_vy_min   = p_group%remapping_grid%eta4_min

    ! Poisson mesh associated to the particles
    m2d => p_group%space_mesh_2d

    SLL_ALLOCATE( particle_indices(number_parts_x, number_parts_y, number_parts_vx, number_parts_vy), ierr )
    particle_indices(:,:,:,:) = 0
    coords(:) = 0

    ! compute the particle weights from the values of f0 on the (cartesian, phase-space) remapping grid
    k_temp_debug = 0
    x_j = parts_x_min
    do j_x = 1, number_parts_x
      y_j = parts_y_min
      do j_y = 1, number_parts_y
        vx_j = parts_vx_min
        do j_vx = 1, number_parts_vx
          vy_j = parts_vy_min
          do j_vy = 1, number_parts_vy

            k_temp_debug = k_temp_debug + 1
            call get_particle_index_from_initial_position_on_cartesian_grid(            &
                j_x, j_y, j_vx, j_vy,                                                   &
                number_parts_x, number_parts_y, number_parts_vx, number_parts_vy,       &
                k                                                                       &
            )

            SLL_ASSERT(k == k_temp_debug)

            if( p_group%spline_degree == 1 )then
                w_k = d_vol * real( p_group%target_values(j_x,j_y,j_vx,j_vy) ,f32)
            else if( p_group%spline_degree == 3 )then
                w_k = 0
                do l_x = -1, 1
                    j_aux_x = j_x + l_x
                    if( p_group%domain_is_periodic(1) )then
                        if( j_aux_x < 1 ) j_aux_x = j_aux_x + number_parts_x
                        if( j_aux_x > number_parts_x ) j_aux_x = j_aux_x - number_parts_x
                    end if
                    if( j_aux_x >= 1 .and. j_aux_x <= number_parts_x )then
                        do l_y = -1, 1
                            j_aux_y = j_y + l_y
                            if( p_group%domain_is_periodic(2) )then
                                if( j_aux_y < 1 ) j_aux_y = j_aux_y + number_parts_y
                                if( j_aux_y > number_parts_y ) j_aux_y = j_aux_y - number_parts_y
                            end if
                            if( j_aux_y >= 1 .and. j_aux_y <= number_parts_y )then
                                do l_vx = -1, 1
                                    j_aux_vx = j_vx + l_vx
                                    if( j_aux_vx >= 1 .and. j_aux_vx <= number_parts_vx )then
                                        do l_vy = -1, 1
                                            j_aux_vy = j_vy + l_vy
                                            if( j_aux_vy >= 1 .and. j_aux_vy <= number_parts_vy )then
                                                ! MCP: here by discarding outside loop instances we assume
                                                ! that non-periodic bc = zero bc. We should instead
                                                ! keep those loop instances and use the specified bounddary condition when
                                                ! the node is in some "fat boundary zone" )
                                                w_k = w_k +         &
                                                    real(           &
                                                    p_group%lt_pic_interpolation_coefs(l_x )  *  &
                                                    p_group%lt_pic_interpolation_coefs(l_y )  *  &
                                                    p_group%lt_pic_interpolation_coefs(l_vx)  *  &
                                                    p_group%lt_pic_interpolation_coefs(l_vy)  *  &
                                                    p_group%target_values(j_aux_x,j_aux_y,j_aux_vx,j_aux_vy) ,f32)
                                            end if
                                        end do
                                    end if
                                end do
                            end if
                        end do
                    end if
                end do
                w_k = d_vol * w_k
            else
               print *, 'sll_bsl_lt_pic_initialize_some4Dfunction(): ERROR, value of p_group%spline_degree ', &
                        ' is invalid: ', p_group%spline_degree
               STOP
            end if

            !> set the position, velocity and weight for the k-th particle
            call p_group%set_particle_weight( k, w_k )

            coords(0) = x_j
            coords(1) = y_j
            call p_group%set_x( k, coords )

            coords(0) = vx_j
            coords(1) = vy_j
            call p_group%set_v( k, coords )

            !            call global_to_cell_offset_extended(    &
            !                    x_j, y_j, &
            !                    m2d,      &
            !                    p_group%particle_list(k)%i_cell_x, &
            !                    p_group%particle_list(k)%i_cell_y, &
            !                    p_group%particle_list(k)%dx, &
            !                    p_group%particle_list(k)%dy )

            !            p_group%particle_list(k)%vx = vx_j
            !            p_group%particle_list(k)%vy = vy_j


            !> set the particle connectivity
            particle_indices( j_x, j_y, j_vx, j_vy ) = k

            if(j_x == 1)then
                ! [neighbor index = own index] means: no neighbor
                ! in the x-periodic case this will be changed when dealing the last particle in the x dimension
                p_group%particle_list(k)%ngb_xleft_index = k
            else
                ! set the connectivity (in both directions) with left neighbor
                k_ngb = particle_indices(j_x-1,j_y,j_vx,j_vy)
                p_group%particle_list(k)%ngb_xleft_index = k_ngb
                p_group%particle_list(k_ngb)%ngb_xright_index = k
                if(j_x == number_parts_x)then
                    if( p_group%domain_is_periodic(1) )then
                        ! set the connectivity (in both directions) with right neighbor
                        k_ngb = particle_indices(1,j_y,j_vx,j_vy)
                        p_group%particle_list(k)%ngb_xright_index = k_ngb
                        p_group%particle_list(k_ngb)%ngb_xleft_index = k
                    else
                        ! [neighbor index = own index] means: no neighbor
                        p_group%particle_list(k)%ngb_xright_index = k
                    end if
                end if
            end if
            if(j_y == 1)then
                ! [neighbor index = own index] means: no neighbor
                ! in the y-periodic case this will be changed when dealing the last particle in the y dimension
                p_group%particle_list(k)%ngb_yleft_index = k
            else
                ! set the connectivity (in both directions) with left neighbor
                k_ngb = particle_indices(j_x,j_y-1,j_vx,j_vy)
                p_group%particle_list(k)%ngb_yleft_index = k_ngb
                p_group%particle_list(k_ngb)%ngb_yright_index = k
                if(j_y == number_parts_y)then
                    if( p_group%domain_is_periodic(2) )then
                        ! set the connectivity (in both directions) with right neighbor
                        k_ngb = particle_indices(j_x,1,j_vx,j_vy)
                        p_group%particle_list(k)%ngb_yright_index = k_ngb
                        p_group%particle_list(k_ngb)%ngb_yleft_index = k
                    else
                        ! [neighbor index = own index] means: no neighbor
                        p_group%particle_list(k)%ngb_yright_index = k
                    end if
                end if
            end if
            if(j_vx == 1)then
                ! [neighbor index = own index] means: no neighbor
                p_group%particle_list(k)%ngb_vxleft_index = k
            else
                ! set the connectivity (in both directions) with left neighbor
                k_ngb = particle_indices(j_x,j_y,j_vx-1,j_vy)
                p_group%particle_list(k)%ngb_vxleft_index = k_ngb
                p_group%particle_list(k_ngb)%ngb_vxright_index = k
                if(j_vx == number_parts_vx)then
                    ! [neighbor index = own index] means: no neighbor
                    p_group%particle_list(k)%ngb_vxright_index = k
                end if
            end if
            if(j_vy == 1)then
                ! [neighbor index = own index] means: no neighbor
                p_group%particle_list(k)%ngb_vyleft_index = k
            else
                ! set the connectivity (in both directions) with left neighbor
                k_ngb = particle_indices(j_x,j_y,j_vx,j_vy-1)
                p_group%particle_list(k)%ngb_vyleft_index = k_ngb
                p_group%particle_list(k_ngb)%ngb_vyright_index = k
                if(j_vy == number_parts_vy)then
                    ! [neighbor index = own index] means: no neighbor
                    p_group%particle_list(k)%ngb_vyright_index = k
                end if
            end if

            vy_j = vy_j + h_parts_vy
          end do
          vx_j = vx_j + h_parts_vx
        end do
        y_j = y_j + h_parts_y
      end do
      x_j = x_j + h_parts_x
    end do

  end subroutine bsl_lt_pic_4d_compute_new_particles


!!!!!!!!!!!! ---------------------------------------------!!!!!!!!!!!! ---------------------------------------------!!!!!!!!!!!! ---------------------------------------------
!!
!!          OLD version
!!
!!!!!!!!!!!! ---------------------------------------------!!!!!!!!!!!! ---------------------------------------------!!!!!!!!!!!! ---------------------------------------------

    !  function sll_lt_pic_4d_group_new( &
    !!  function new_lt_particle_4d_group( &   ! old name
    !        spline_degree,     &
    !        number_parts_x,    &
    !        number_parts_y,    &
    !        number_parts_vx,   &
    !        number_parts_vy,   &
    !        remap_grid_vx_min, &
    !        remap_grid_vx_max, &
    !        remap_grid_vy_min, &
    !        remap_grid_vy_max, &
    !        !        particle_array_size, &
    !        !        guard_list_size,     &
    !        qoverm,              &
    !        domain_is_x_periodic,&
    !        domain_is_y_periodic,&
    !        mesh ) result(res)
    !
    !    type(sll_lt_pic_4d_group), pointer :: res
    !
    !    sll_int32, intent(in)   :: spline_degree
    !    sll_int32, intent(in)   :: number_parts_x
    !    sll_int32, intent(in)   :: number_parts_y
    !    sll_int32, intent(in)   :: number_parts_vx
    !    sll_int32, intent(in)   :: number_parts_vy
    !    sll_real64, intent(in)  :: remap_grid_vx_min
    !    sll_real64, intent(in)  :: remap_grid_vx_max
    !    sll_real64, intent(in)  :: remap_grid_vy_min
    !    sll_real64, intent(in)  :: remap_grid_vy_max
    !    sll_real64, intent(in)  :: qoverm
    !    logical, intent(in)      :: domain_is_x_periodic
    !    logical, intent(in)      :: domain_is_y_periodic
    !
    !    !        sll_int32, intent(in)   :: particle_array_size
    !    !        sll_int32, intent(in)   :: guard_list_size
    !    sll_int32   :: particle_array_size
    !    sll_int32   :: guard_list_size
    !
    !    type(sll_cartesian_mesh_2d), pointer :: mesh
    !    sll_int32   :: ierr
    !    sll_int32   :: number_particles
    !    sll_int32   :: remap_grid_number_cells_x
    !    sll_int32   :: remap_grid_number_cells_y
    !    sll_int32   :: remap_grid_number_cells_vx
    !    sll_int32   :: remap_grid_number_cells_vy
    !    sll_real64  :: remap_grid_x_min
    !    sll_real64  :: remap_grid_x_max
    !    sll_real64  :: remap_grid_y_min
    !    sll_real64  :: remap_grid_y_max
    !
    !    number_particles = number_parts_x * number_parts_y * number_parts_vx * number_parts_vy
    !    particle_array_size = number_particles
    !    guard_list_size = number_particles
    !    !    if( number_particles > particle_array_size ) then
    !    !       print *, 'sll_lt_pic_4d_group_new(): ERROR (code=6454357),  number_particles should not ', &
    !    !            'be greater than the requested memory size, particle_array_size.'
    !    !       print *, 'note: number_particles, particle_array_size = ', number_particles, particle_array_size
    !    !       STOP
    !    !    end if
    !
    !    SLL_ALLOCATE( res, ierr )
    !    res%spline_degree = spline_degree
    !    res%number_parts_x   = number_parts_x
    !    res%number_parts_y   = number_parts_y
    !    res%number_parts_vx  = number_parts_vx
    !    res%number_parts_vy  = number_parts_vy
    !    res%number_particles = number_particles
    !    res%active_particles = number_particles
    !    res%guard_list_size  = guard_list_size
    !    res%qoverm           = qoverm
    !    res%domain_is_x_periodic = domain_is_x_periodic
    !    res%domain_is_y_periodic = domain_is_y_periodic
    !    res%track_markers_outside_domain = .false.
    !    res%use_exact_f0 = .false.
    !    SLL_ALLOCATE( res%p_list(particle_array_size), ierr )
    !    SLL_ALLOCATE( res%p_guard(guard_list_size), ierr )
    !    SLL_ALLOCATE( res%target_values(number_parts_x, number_parts_y, number_parts_vx, number_parts_vy), ierr )
    !
    !    ! MCP [DEBUG]
    !    SLL_ALLOCATE( res%debug_bsl_remap(number_parts_x, number_parts_y, number_parts_vx, number_parts_vy,4,3), ierr )
    !
    !    ! physical mesh, used for Poisson solver
    !    if (.not.associated(mesh) ) then
    !       print*, 'sll_lt_pic_4d_group_new(): ERROR, passed mesh not associated'
    !    endif
    !    res%mesh => mesh
    !
    !    ! phase space grid to initialize and remap the particles
    !    remap_grid_x_min = mesh%eta1_min
    !    remap_grid_x_max = mesh%eta1_max
    !    remap_grid_y_min = mesh%eta2_min
    !    remap_grid_y_max = mesh%eta2_max
    !
    !    if( domain_is_x_periodic )then
    !        remap_grid_number_cells_x = number_parts_x
    !    else
    !        remap_grid_number_cells_x = number_parts_x - 1
    !    end if
    !    if( domain_is_y_periodic )then
    !        remap_grid_number_cells_y = number_parts_y
    !    else
    !        remap_grid_number_cells_y = number_parts_y - 1
    !    end if
    !    remap_grid_number_cells_vx = number_parts_vx - 1
    !    remap_grid_number_cells_vy = number_parts_vy - 1
    !    res%remapping_grid => new_cartesian_mesh_4d(  remap_grid_number_cells_x,        &
    !                                                remap_grid_number_cells_y,        &
    !                                                remap_grid_number_cells_vx,       &
    !                                                remap_grid_number_cells_vy,       &
    !                                                remap_grid_x_min,   &
    !                                                remap_grid_x_max,   &
    !                                                remap_grid_y_min,   &
    !                                                remap_grid_y_max,   &
    !                                                remap_grid_vx_min,  &
    !                                                remap_grid_vx_max,  &
    !                                                remap_grid_vy_min,  &
    !                                                remap_grid_vy_max   &
    !                                              )
    !
    !
    !
    !    if( spline_degree == 3) then  ! ( otherwise, no need... )
    !        SLL_ALLOCATE( res%ltpic_interpolation_coefs(-1:1), ierr )
    !           res%ltpic_interpolation_coefs(-1) = -1./6.
    !           res%ltpic_interpolation_coefs(0)  =  8./6.
    !           res%ltpic_interpolation_coefs(1)  =  -1./6.
    !    end if
    !
    !
    !
    !  end function sll_lt_pic_4d_group_new



  !----------------------------------------------------------------------------
  ! Destructor
  subroutine sll_bsl_lt_pic_4d_group_delete(particle_group)
    class(sll_bsl_lt_pic_4d_group), pointer :: particle_group
    sll_int32 :: ierr

    if(.not. associated(particle_group) ) then
       print *, 'sll_bsl_lt_pic_4d_group_delete(): ERROR (9087987578996), passed group was not associated.'
    end if
    SLL_DEALLOCATE(particle_group%particle_list, ierr)
    SLL_DEALLOCATE(particle_group%space_mesh_2d, ierr)
    SLL_DEALLOCATE(particle_group, ierr)

  end subroutine sll_bsl_lt_pic_4d_group_delete


end module sll_bsl_lt_pic_4d_group_module
