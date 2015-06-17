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


module sll_simple_pic_4d_group_module

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
! #include "particle_representation.h"   NEEDED?

  use sll_working_precision
  use sll_simple_pic_4d_particle_module
  use sll_cartesian_meshes
  use sll_module_pic_base

  implicit none


  type, extends(sll_particle_group_base) :: sll_simple_pic_4d_group

    !    ! the group
    !    class( sll_species ), pointer   :: species
        !    sll_int32                       :: id

    ! the particles
    sll_int32                                                   :: number_particles
    type(sll_simple_pic_4d_particle),   dimension(:), pointer   :: particle_list
    sll_real64                                                  :: common_weight

    ! the physical mesh used eg in the Poisson solver
    type(sll_cartesian_mesh_2d), pointer    :: space_mesh_2d
    logical                                 :: domain_is_x_periodic
    logical                                 :: domain_is_y_periodic

    ! the initial density (put this in a separate object ?)
    sll_real64      :: thermal_speed
    sll_real64      :: alpha
    sll_real64      :: k_landau

  contains

    ! Getters
    procedure :: get_x          => simple_pic_4d_get_x
    procedure :: get_v          => simple_pic_4d_get_v
    procedure :: get_charge     => simple_pic_4d_get_charge
    procedure :: get_mass       => simple_pic_4d_get_mass
    procedure :: get_cell_index => simple_pic_4d_get_cell_index

    ! Setters
    procedure :: set_x                  => simple_pic_4d_set_x
    procedure :: set_v                  => simple_pic_4d_set_v
    procedure :: set_common_weight      => simple_pic_4d_set_common_weight
    procedure :: set_particle_weight    => simple_pic_4d_set_particle_weight

    ! Initializers
    procedure :: set_landau_params      => simple_pic_4d_set_landau_params
    procedure :: random_initializer     => simple_pic_4d_random_initializer
    procedure :: cartesian_initializer  => simple_pic_4d_cartesian_initializer


    procedure :: deposit_charge_2d          => simple_pic_4d_deposit_charge_2d

  end type sll_simple_pic_4d_group

  interface sll_delete
     module procedure sll_simple_pic_4d_group_delete
  end interface sll_delete

contains


  !----------------------------------------------------------------------------
  pure function simple_pic_4d_get_charge( self, i ) result( r )
    class( sll_simple_pic_4d_group ), intent( in ) :: self
    sll_int32                       , intent( in ) :: i
    sll_real64 :: r

    r = self%species%q * self%common_weight

  end function simple_pic_4d_get_charge


  !----------------------------------------------------------------------------
  pure function simple_pic_4d_get_mass( self, i ) result( r )
    class( sll_simple_pic_4d_group ), intent( in ) :: self
    sll_int32                       , intent( in ) :: i
    sll_real64 :: r

    r = self%species%m * self%common_weight

  end function simple_pic_4d_get_mass


  !----------------------------------------------------------------------------
  pure function simple_pic_4d_get_x( self, i ) result( r )
    class( sll_simple_pic_4d_group ), intent( in ) :: self
    sll_int32                       , intent( in ) :: i
    sll_real64 :: r(3)

    type(sll_cartesian_mesh_2d),      pointer :: space_mesh_2d
    type(sll_simple_pic_4d_particle), pointer :: particle

    space_mesh_2d => self%space_mesh_2d
    particle => self%particle_list(i)

    ! get x
    r(1) = space_mesh_2d%eta1_min + space_mesh_2d%delta_eta1*( particle%offset_x + real(particle%i_cell_x - 1, f64) )
    ! get y
    r(2) = space_mesh_2d%eta2_min + space_mesh_2d%delta_eta2*( particle%offset_y + real(particle%i_cell_y - 1, f64) )

  end function simple_pic_4d_get_x


  !----------------------------------------------------------------------------
  pure function simple_pic_4d_get_v( self, i ) result( r )
    class( sll_simple_pic_4d_group ), intent( in ) :: self
    sll_int32                       , intent( in ) :: i
    sll_real64 :: r(3)

    type(sll_simple_pic_4d_particle), pointer :: particle

    particle => self%particle_list(i)

    ! get vx
    r(1) = particle%v_x
    ! get vy
    r(2) = particle%v_y

  end function simple_pic_4d_get_v


  !----------------------------------------------------------------------------
  pure function ::  simple_pic_4d_get_cell_index(particle_group, i_part) res(i_cell)
    type(sll_simple_pic_4d_group),  intent( in )    ::  particle_group
    sll_int32,                      intent( in )    ::  i_part
    sll_int32,                      intent( out )   ::  i_cell
    sll_int32 ::  i_cell_x, i_cell_y
    sll_int32 ::  num_cells_x, num_cells_y

    i_cell_x = particle_group%particle_list(i_part)%i_cell_x
    i_cell_y = particle_group%particle_list(i_part)%i_cell_y
    num_cells_x = particle_group%space_mesh_2d%num_cells1
    num_cells_y = particle_group%space_mesh_2d%num_cells2

    i_cell = 1 + modulo(i_cell_x - 1,  num_cells_x) + modulo(i_cell_y - 1,  num_cells_y) * num_cells_x

    SLL_ASSERT( i_cell >= 1)
    SLL_ASSERT( i_cell <= num_cells_x * num_cells_y )

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
  subroutine simple_pic_4d_set_v( self, i, v )
    class( sll_simple_pic_4d_group ), intent( inout ) :: self
    sll_int32                       , intent( in    ) :: i
    sll_real64                      , intent( in    ) :: v(3)

    type(sll_simple_pic_4d_particle), pointer :: particle

    particle => self%particle_list(i)
    particle%v_x = v(1)
    particle%v_y = v(2)

  end subroutine simple_pic_4d_set_v


  !----------------------------------------------------------------------------
  subroutine simple_pic_4d_set_common_weight( self, w )
    class( sll_simple_pic_4d_group ), intent( inout ) :: self
    sll_real64                      , intent( in    ) :: w

    self%common_weight = w

  end subroutine simple_pic_4d_set_common_weight


  !----------------------------------------------------------------------------
  subroutine simple_pic_4d_set_particle_weight( self, i, w )
    class( sll_simple_pic_4d_group ), intent( inout ) :: self
    sll_int32                       , intent( in    ) :: i
    sll_real64                      , intent( in    ) :: w

    print*, "Error (8654354237645) -- this subroutine is not implemented for simple_pic_4d_group objects"
    stop

  end subroutine simple_pic_4d_set_particle_weight


  !----------------------------------------------------------------------------
  subroutine simple_pic_4d_set_landau_params( self, thermal_speed, alpha, k_landau )
    class( sll_simple_pic_4d_group ), intent( inout ) :: self
    sll_real64                      , intent( in    ) :: thermal_speed
    sll_real64                      , intent( in    ) :: alpha
    sll_real64                      , intent( in    ) :: k_landau

    self%thermal_speed = thermal_speed
    self%alpha = alpha
    self%k_landau = k_landau

  end subroutine simple_pic_4d_set_landau_params


  !----------------------------------------------------------------------------
   subroutine simple_pic_4d_random_initializer( self, initial_density_identifier, rand_seed, rank, world_size )
    class( sll_simple_pic_4d_group ), intent( inout ) :: self
    sll_int32                       , intent( in    ) :: initial_density_identifier
    sll_int32, dimension(:)         , intent( in ), optional :: rand_seed
    sll_int32                       , intent( in ), optional :: rank, world_size

    ! for the moment we only use the landau damping initial density
    ! so we don't use the initial_density_identifier, but eventually it should say which initial density is used

    call sll_pic_4d_random_unweighted_initializer_landau_f0 (   &
      self%thermal_speed, self%alpha, self%k_landau,            & ! -> these parameters should be members of the initializer object
      self,                                                     &
      self%space_mesh_2d,                                       &
      self%number_particles                                     &
      rand_seed, rank, world_size                               &
    )

   end subroutine simple_pic_4d_random_initializer



  !----------------------------------------------------------------------------
   subroutine simple_pic_4d_cartesian_initializer( self, nb_cells_x, nb_cells_v, x_min, x_max, v_min, v_max )
    class( sll_simple_pic_4d_group ), intent( inout ) :: self
    sll_int32                       , intent( in    ) :: nb_cells_x(3)
    sll_int32                       , intent( in    ) :: nb_cells_v(3)
    sll_real64                      , intent( in    ) :: x_min(3)
    sll_real64                      , intent( in    ) :: x_max(3)
    sll_real64                      , intent( in    ) :: v_min(3)
    sll_real64                      , intent( in    ) :: v_max(3)

    print*, "Error (987856) -- this subroutine is not implemented for simple_pic_4d_group objects"
    stop

   end subroutine simple_pic_4d_cartesian_initializer


  subroutine simple_pic_4d_deposit_charge_2d( self, charge_accumulator )
    class( sll_simple_pic_4d_group ),           intent( inout ) :: self
    type( sll_charge_accumulator_2d ), pointer, intent( inout ) :: charge_accumulator

    sll_int32,      :: i_part
    sll_real64,     :: x_part(3)
    sll_real64,     :: particle_charge
    sll_real64,     :: tmp1, tmp2

    call reset_charge_accumulator_to_zero ( charge_accumulator )

    particle_charge = particle_group%get_species_charge
    do i_part = 1, self%number_particles

      particle => self%particle_list( i_part )
      dx = particle%offset_x
      dy = particle%offset_y
      i_cell = get_cell_index(particle_group, i_part)
      charge_accumulator_cell => charge_accumulator%q_acc(i_cell)

      xy_part = self%get_x( i_part )  ! x and y
      if( x_is_in_domain_2d( x_part, self%space_mesh_2d, self%domain_is_x_periodic, self%domain_is_y_periodic ) )then

        tmp1 = (1.0_f64 - dx)
        tmp2 = (1.0_f64 - dy)
        charge_accumulator_cell%q_sw = charge_accumulator_cell%q_sw + q*tmp1*tmp2
        charge_accumulator_cell%q_se = charge_accumulator_cell%q_se + q*  dx*tmp2
        charge_accumulator_cell%q_nw = charge_accumulator_cell%q_nw + q*tmp1*  dy
        charge_accumulator_cell%q_ne = charge_accumulator_cell%q_ne + q*  dx*  dy

      else ! particle not in domain (should store the reference for later processing)
        print*, "Error (097647687): for the moment every particle should be in the (periodic) 2d domain..."
        stop
      end if

    end do

  end subroutine simple_pic_4d_deposit_charge_2d

  !----------------------------------------------------------------------------
  ! Constructor
  function sll_simple_pic_4d_group_new( &
        number_particles,       &
        species_charge,         &
        species_mass,           &
        particle_group_id,      &
        domain_is_x_periodic,   &
        domain_is_y_periodic,   &
        space_mesh_2d ) result(res)

    type(sll_lt_pic_4d_group), pointer :: res

    sll_real64, intent(in)  :: species_charge
    sll_real64, intent(in)  :: species_mass
    type(sll_cartesian_mesh_2d), pointer :: space_mesh_2d

    logical, intent(in)      :: domain_is_x_periodic
    logical, intent(in)      :: domain_is_y_periodic

    SLL_ALLOCATE( res, ierr )

    ! the group
    res%species => species_new( species_charge, species_mass )
    res%id = particle_group_id
    res%dimension_x = 2
    res%dimension_v = 2

    ! the particles
    res%number_particles = number_particles
    SLL_ALLOCATE( res%particle_list(number_particles), ierr )

    ! physical 2d mesh, used eg in the Poisson solver
    if (.not.associated(space_mesh_2d) ) then
       print*, 'ERROR (876876456), given space_mesh_2d is not associated'
    end if
    res%space_mesh_2d => space_mesh_2d
    res%domain_is_x_periodic = domain_is_x_periodic
    res%domain_is_y_periodic = domain_is_y_periodic

    end function sll_simple_pic_4d_group_new


  !----------------------------------------------------------------------------
  ! Destructor
  subroutine sll_simple_pic_4d_group_delete(particle_group)
    type(sll_simple_pic_4d_group), pointer :: particle_group
    sll_int32 :: ierr

    if(.not. associated(particle_group) ) then
       print *, 'sll_simple_pic_4d_group_delete(): ERROR (986876), passed group was not associated.'
    end if
    SLL_DEALLOCATE(particle_group%particle_list, ierr)
    SLL_DEALLOCATE(particle_group%mesh, ierr)
    SLL_DEALLOCATE(p_group, ierr)

  end subroutine sll_simple_pic_4d_group_delete



  ! tells whether the given point is in the given domain, with boolean arguments for the domain periodicity
  ! (taken from previous function in_bounds_periodic)
  function x_is_in_domain_2d( x, mesh, x_periodic, y_periodic ) result(res)

    sll_real64,                     intent( in )    :: x(3)             ! dimension 3 to fit the general interface
    logical,                        intent( in )    :: x_periodic
    logical,                        intent( in )    :: y_periodic
    type(sll_cartesian_mesh_2d),    intent( in ), pointer :: mesh
    logical     :: res

    res = ( x(1) >= mesh%eta1_min )                                                                               &
          .and.                                                                                                   &
          ( ( x(1) < mesh%eta1_max .and. x_periodic ) .or. ( x(1) <= mesh%eta1_max .and. .not. x_periodic ) )     &
          .and.                                                                                                   &
          ( x(2) >= mesh%eta2_min )                                                                               &
          .and.                                                                                                   &
          ( ( x(2) < mesh%eta2_max .and. y_periodic ) .or. ( x(2) <= mesh%eta2_max .and. .not. y_periodic) )

  end function x_is_in_domain_2d


end module sll_simple_pic_4d_group_module


































! --------------   ancien module ----------------


     sll_int32 :: spline_degree
     sll_int32 :: number_parts_x
     sll_int32 :: number_parts_y
     sll_int32 :: number_parts_vx
     sll_int32 :: number_parts_vy
     sll_int32 :: number_particles
     sll_int32 :: active_particles
     sll_int32 :: guard_list_size
     sll_real64 :: qoverm 
     logical    :: domain_is_x_periodic
     logical    :: domain_is_y_periodic
     logical    :: track_markers_outside_domain   !! default value is true
     logical    :: use_exact_f0                   ! if false, interpolate f0 from its values on the initial/remapped particle grid
     type(sll_cartesian_mesh_2d), pointer                     :: mesh

     ! <<sll_lt_pic_4d_group-p_list>> uses [[file:lt_pic_4d_particle.F90::sll_lt_pic_4d_particle]]
     type(sll_lt_pic_4d_particle), dimension(:), pointer        :: p_list
     type(sll_lt_pic_4d_particle_guard_ptr), dimension(:), pointer :: p_guard

     ! num_postprocess_particles: an array indexed by the thread number (if any),
     ! of the number of particles to post-process after the main loop
     sll_int32, dimension(:), pointer :: num_postprocess_particles

     ! <<sll_lt_pic_4d_group-remapping_grid>> uses [[file:../meshes/sll_cartesian_meshes.F90::sll_cartesian_mesh_4d]]
     type(sll_cartesian_mesh_4d),               pointer       :: remapping_grid
!     sll_real64, dimension(:,:,:,:),          pointer       :: remapped_values     ! old name
     sll_real64, dimension(:,:,:,:),          pointer       :: target_values

     ! MCP [DEBUG]
     sll_real64, dimension(:,:,:,:,:,:),          pointer       :: debug_bsl_remap

     sll_real64, dimension(:),                pointer       :: ltpic_interpolation_coefs

    ! temp fields (landau damping)
    sll_real64 :: thermal_speed
    sll_real64 :: alpha_landau
    sll_real64 :: k_landau


  end type sll_lt_pic_4d_group


  interface sll_delete
     module procedure sll_lt_pic_4d_group_delete
!     module procedure delete_lt_particle_4d_group  ! old name
  end interface sll_delete

contains


  function sll_lt_pic_4d_group_new( &
!  function new_lt_particle_4d_group( &   ! old name
        spline_degree,     &
        number_parts_x,    &
        number_parts_y,    &
        number_parts_vx,   &
        number_parts_vy,   &
        remap_grid_vx_min, &
        remap_grid_vx_max, &
        remap_grid_vy_min, &
        remap_grid_vy_max, &       
        !        particle_array_size, &
        !        guard_list_size,     &
        qoverm,              &
        domain_is_x_periodic,&
        domain_is_y_periodic,&
        mesh ) result(res)

    type(sll_lt_pic_4d_group), pointer :: res

    sll_int32, intent(in)   :: spline_degree
    sll_int32, intent(in)   :: number_parts_x
    sll_int32, intent(in)   :: number_parts_y
    sll_int32, intent(in)   :: number_parts_vx
    sll_int32, intent(in)   :: number_parts_vy
    sll_real64, intent(in)  :: remap_grid_vx_min   
    sll_real64, intent(in)  :: remap_grid_vx_max   
    sll_real64, intent(in)  :: remap_grid_vy_min   
    sll_real64, intent(in)  :: remap_grid_vy_max   
    sll_real64, intent(in)  :: qoverm
    logical, intent(in)      :: domain_is_x_periodic
    logical, intent(in)      :: domain_is_y_periodic
    
    !        sll_int32, intent(in)   :: particle_array_size
    !        sll_int32, intent(in)   :: guard_list_size
    sll_int32   :: particle_array_size
    sll_int32   :: guard_list_size

    type(sll_cartesian_mesh_2d), pointer :: mesh
    sll_int32   :: ierr
    sll_int32   :: number_particles
    sll_int32   :: remap_grid_number_cells_x
    sll_int32   :: remap_grid_number_cells_y
    sll_int32   :: remap_grid_number_cells_vx
    sll_int32   :: remap_grid_number_cells_vy
    sll_real64  :: remap_grid_x_min   
    sll_real64  :: remap_grid_x_max   
    sll_real64  :: remap_grid_y_min   
    sll_real64  :: remap_grid_y_max   
        
    number_particles = number_parts_x * number_parts_y * number_parts_vx * number_parts_vy
    particle_array_size = number_particles
    guard_list_size = number_particles
    !    if( number_particles > particle_array_size ) then
    !       print *, 'sll_lt_pic_4d_group_new(): ERROR (code=6454357),  number_particles should not ', &
    !            'be greater than the requested memory size, particle_array_size.'
    !       print *, 'note: number_particles, particle_array_size = ', number_particles, particle_array_size
    !       STOP
    !    end if

    SLL_ALLOCATE( res, ierr )
    res%spline_degree = spline_degree
    res%number_parts_x   = number_parts_x
    res%number_parts_y   = number_parts_y
    res%number_parts_vx  = number_parts_vx
    res%number_parts_vy  = number_parts_vy
    res%number_particles = number_particles
    res%active_particles = number_particles
    res%guard_list_size  = guard_list_size
    res%qoverm           = qoverm
    res%domain_is_x_periodic = domain_is_x_periodic
    res%domain_is_y_periodic = domain_is_y_periodic
    res%track_markers_outside_domain = .false.
    res%use_exact_f0 = .false.
    SLL_ALLOCATE( res%p_list(particle_array_size), ierr )
    SLL_ALLOCATE( res%p_guard(guard_list_size), ierr )
    SLL_ALLOCATE( res%target_values(number_parts_x, number_parts_y, number_parts_vx, number_parts_vy), ierr )

    ! MCP [DEBUG]
    SLL_ALLOCATE( res%debug_bsl_remap(number_parts_x, number_parts_y, number_parts_vx, number_parts_vy,4,3), ierr )

    ! physical mesh, used for Poisson solver
    if (.not.associated(mesh) ) then
       print*, 'sll_lt_pic_4d_group_new(): ERROR, passed mesh not associated'
    endif
    res%mesh => mesh

    ! phase space grid to initialize and remap the particles    
    remap_grid_x_min = mesh%eta1_min
    remap_grid_x_max = mesh%eta1_max
    remap_grid_y_min = mesh%eta2_min
    remap_grid_y_max = mesh%eta2_max
    
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
    res%remapping_grid => new_cartesian_mesh_4d(  remap_grid_number_cells_x,        &
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
  
 

    if( spline_degree == 3) then  ! ( otherwise, no need... )   
        SLL_ALLOCATE( res%ltpic_interpolation_coefs(-1:1), ierr )
           res%ltpic_interpolation_coefs(-1) = -1./6.
           res%ltpic_interpolation_coefs(0)  =  8./6.
           res%ltpic_interpolation_coefs(1)  =  -1./6.
    end if

  
  
  end function sll_lt_pic_4d_group_new


  subroutine sll_lt_pic_4d_group_delete(p_group)
    type(sll_lt_pic_4d_group), pointer :: p_group
    sll_int32 :: ierr

    if(.not. associated(p_group) ) then
       print *, 'delete_lt_particle_group_4d(): ERROR, passed group was not ', &
            'associated.'
    end if
    SLL_DEALLOCATE(p_group%p_list, ierr)
    SLL_DEALLOCATE(p_group%p_guard, ierr)
    SLL_DEALLOCATE(p_group%target_values, ierr)

    ! MCP [DEBUG]
    SLL_DEALLOCATE( p_group%debug_bsl_remap, ierr )

    SLL_DEALLOCATE(p_group, ierr)
    
  end subroutine sll_lt_pic_4d_group_delete


end module sll_lt_pic_4d_grou_modulee
