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

!> @brief Module for groups of particles of type sll_bsl_lt_pic_4d_particle

module sll_m_bsl_lt_pic_4d_group

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_accumulators.h"
#include "sll_errors.h"

  use sll_m_constants, only: sll_pi
  use sll_m_working_precision
  use sll_m_cartesian_meshes
  use sll_m_sparse_grid_4d, only: sparse_grid_interpolator_4d
  use sll_m_remapped_pic_base
  use sll_m_bsl_lt_pic_4d_particle
  use sll_m_bsl_lt_pic_4d_utilities !, only: int_list_element, int_list_element_ptr, add_element_in_list
  use sll_m_remapped_pic_utilities, only: x_is_in_domain_2d, apply_periodic_bc_on_cartesian_mesh_2d
  use sll_m_gnuplot
  use sll_m_sobol, only: i8_sobol

  implicit none

  ! types of interpolation for the remapped f
  sll_int32, parameter :: SLL_BSL_LT_PIC_REMAP_WITH_SPLINES = 0
  sll_int32, parameter :: SLL_BSL_LT_PIC_REMAP_WITH_SPARSE_GRIDS = 1

  ! possible scenarios for the reconstruction routine
  sll_int32, parameter :: SLL_BSL_LT_PIC_DEPOSIT_F = 0
  sll_int32, parameter :: SLL_BSL_LT_PIC_REMAP_F = 1
  sll_int32, parameter :: SLL_BSL_LT_PIC_WRITE_F_ON_GIVEN_GRID = 2
  sll_int32, parameter :: SLL_BSL_LT_PIC_SET_WEIGHTS_ON_DEPOSITION_PARTICLES = 3

  ! possible densities for f0
  sll_int32, parameter :: SLL_BSL_LT_PIC_LANDAU_F0 = 0
  sll_int32, parameter :: SLL_BSL_LT_PIC_HAT_F0 = 1

  ! types of deposition particles
  sll_int32, parameter :: SLL_BSL_LT_PIC_FIXED_GRID = 0
  sll_int32, parameter :: SLL_BSL_LT_PIC_TRANSPORTED_RANDOM = 1


  !> Group of @ref sll_bsl_lt_pic_4d_particle
  type, extends(sll_c_remapped_particle_group) :: sll_bsl_lt_pic_4d_group

    !> @name The markers (particles pushed forward, carry no weights)
    !> @{
    sll_int32                                                   :: number_flow_markers_x
    sll_int32                                                   :: number_flow_markers_y
    sll_int32                                                   :: number_flow_markers_vx
    sll_int32                                                   :: number_flow_markers_vy
    sll_int32                                                   :: number_flow_markers
    type(sll_cartesian_mesh_4d), pointer                        :: initial_markers_grid
    type(sll_bsl_lt_pic_4d_particle),   dimension(:), pointer   :: markers_list
    !> @}

    !> @name The flow grid (4d cartesian cells where the flow is linearized)
    !> @{
    type(sll_cartesian_mesh_4d), pointer    :: flow_grid
    !> @}

    !> @name The physical mesh used eg in the Poisson solver
    !> @{
    type(sll_cartesian_mesh_2d), pointer    :: space_mesh_2d
    !> @}

    !> @name The deposition particles (will be created on the fly in each cell of the flow_grid, when depositing the charge)
    !> @{
    sll_int32                                 :: deposition_particles_type          ! fixed_grid or transported
    type(sll_cartesian_mesh_4d), pointer      :: deposition_grid                    !< used if type = fixed_grid
    sll_int32                                 :: number_moving_deposition_particles !<  = 0  if  type == fixed_grid
    sll_real64, dimension(:,:), allocatable   :: deposition_particles_eta           !< used if type = transported
    sll_real64, dimension(:), allocatable     :: deposition_particles_weight        !< used if type = transported (weight = charge)
    !> @}

    !> @name General parameters for the interpolation of the remapped density f
    !> @{
    sll_int32                                                   :: remapped_f_interpolation_type  !> 0 = splines, 1 = sparse grid
    sll_int32                                                   :: remapped_f_interpolation_degree
    sll_real64, dimension(4)                                    :: remapping_grid_eta_min   ! x-y domain: same as for Poisson grid
    sll_real64, dimension(4)                                    :: remapping_grid_eta_max   ! vx-vy domain specified in addition
    !> @}

    !> @name The sparse grid object used for the interpolation of the remapped density f, if remapping with sparse grid
    !> @{
    type(sparse_grid_interpolator_4d)                           :: sparse_grid_interpolator
    sll_int32,  dimension(4)                                    :: sparse_grid_max_levels
    sll_real64, dimension(:), allocatable                       :: remapped_f_sparse_grid_coefficients
    ! maybe more stuff is needed here
    !> @}

    !> @name The cartesian grid used for the interpolation of the remapped density f, if remapping with splines
    !> @{
    type(sll_cartesian_mesh_4d),                pointer         :: remapping_cart_grid
    sll_real64, dimension(:,:,:,:), allocatable                 :: remapped_f_splines_coefficients
    sll_int32                                                   :: remapping_cart_grid_number_cells_x
    sll_int32                                                   :: remapping_cart_grid_number_cells_y
    sll_int32                                                   :: remapping_cart_grid_number_cells_vx
    sll_int32                                                   :: remapping_cart_grid_number_cells_vy
    !> quasi-interpolation coefs are needed for cubic (or higher degree) spline interpolation
    sll_real64, dimension(:),                   pointer         :: lt_pic_interpolation_coefs

    !> @name The initial density (at some point this should be put in a separate initializer object)
    !> @{
    sll_real64                      :: thermal_speed
    sll_real64                      :: alpha
    sll_real64                      :: k_landau

    sll_real64                      :: hat_f0_x0
    sll_real64                      :: hat_f0_y0
    sll_real64                      :: hat_f0_vx0
    sll_real64                      :: hat_f0_vy0
    sll_real64                      :: hat_f0_r_x
    sll_real64                      :: hat_f0_r_y
    sll_real64                      :: hat_f0_r_vx
    sll_real64                      :: hat_f0_r_vy
    sll_real64                      :: hat_f0_basis_height
    sll_real64                      :: hat_f0_shift
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
    procedure :: set_x                      => bsl_lt_pic_4d_set_x
    procedure :: set_v                      => bsl_lt_pic_4d_set_v

    ! todo: use only one function with particle index as optional parameter?
    procedure :: set_common_weight          => bsl_lt_pic_4d_set_common_weight     ! not to be called for this class
    procedure :: set_particle_weight        => bsl_lt_pic_4d_set_particle_weight
    !> @}
    
    !> @name Initializers
    !> @{
    procedure :: set_landau_parameters      =>  bsl_lt_pic_4d_set_landau_parameters
    procedure :: set_hat_f0_parameters      =>  bsl_lt_pic_4d_set_hat_f0_parameters
    procedure :: initializer                =>  bsl_lt_pic_4d_initializer
    !> @}
    
    procedure :: deposit_charge_2d          => bsl_lt_pic_4d_deposit_charge_2d
    procedure :: remap                      => bsl_lt_pic_4d_remap
    procedure :: visualize_f_slice_x_vx     => bsl_lt_pic_4d_visualize_f_slice_x_vx

    procedure :: remapping_cart_grid_number_nodes_x
    procedure :: remapping_cart_grid_number_nodes_y
    procedure :: remapping_cart_grid_number_nodes_vx
    procedure :: remapping_cart_grid_number_nodes_vy

    procedure :: bsl_lt_pic_4d_write_hat_density_on_remapping_grid        !> this evaluates an analytic f0
    procedure :: bsl_lt_pic_4d_write_landau_density_on_remapping_grid     !> this evaluates an analytic f0
    procedure :: bsl_lt_pic_4d_remap_f       !> this evaluates f with the bs_lt_pic method and compute the new interpolation coefs

!    procedure :: bsl_lt_pic_4d_compute_new_remapped_f_coefficients
    procedure :: bsl_lt_pic_4d_compute_new_spline_coefs

    procedure :: bsl_lt_pic_4d_reset_markers_position
    procedure :: bsl_lt_pic_4d_set_markers_connectivity

    procedure :: bsl_lt_pic_set_deposition_particles_coordinates
    procedure :: bsl_lt_pic_set_deposition_particles_weights

    procedure :: bsl_lt_pic_4d_write_f_on_grid_or_deposit
    procedure :: bsl_lt_pic_4d_interpolate_value_of_remapped_f

    procedure :: get_ltp_deformation_matrix
    procedure :: periodic_correction

  end type sll_bsl_lt_pic_4d_group

  interface sll_delete
     module procedure sll_bsl_lt_pic_4d_group_delete
  end interface sll_delete


  !! MCP (July 16) -- this is to make the subroutine external, in a separate file, but does not work yet --

!  interface
!    subroutine bsl_lt_pic_4d_initializer( self, initial_density_identifier, rand_seed, rank, world_size )
!    use sll_m_working_precision
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

    r = self%species%q * self%markers_list(i)%weight

  end function bsl_lt_pic_4d_get_charge


  !----------------------------------------------------------------------------
  pure function bsl_lt_pic_4d_get_mass( self, i ) result( r )
    class( sll_bsl_lt_pic_4d_group ), intent( in ) :: self
    sll_int32                       , intent( in ) :: i
    sll_real64 :: r

    r = self%species%m * self%markers_list(i)%weight

  end function bsl_lt_pic_4d_get_mass


  !----------------------------------------------------------------------------------------------------------------------------
  ! gets the physical coordinates of a 'particle', which can be of two types:
  !   1) either a flow marker (for i = 1, ... self%number_flow_markers)
  !   2) or a deposition particle (for i = self%number_flow_markers+1, ... self%number_moving_deposition_particles)

  pure function bsl_lt_pic_4d_get_x( self, i ) result( r )
    class( sll_bsl_lt_pic_4d_group ), intent( in ) :: self
    sll_int32                       , intent( in ) :: i
    sll_real64 :: r(3)

    if( i < 1 .or. i > self%number_flow_markers + self%number_moving_deposition_particles )then
      ! returning an error would be nice...
      r(1000) = 0   ! will this be caught as an error ?
      return
    end if

    if( i >= 1 .and. i <= self%number_flow_markers )then
      ! then the particle is a flow marker

      ! get x
      r(1) = self%space_mesh_2d%eta1_min + &
             self%space_mesh_2d%delta_eta1*(                            &
             real(self%markers_list(i)%offset_x + self%markers_list(i)%i_cell_x - 1, f64)      )
      ! get y
      r(2) = self%space_mesh_2d%eta2_min + self%space_mesh_2d%delta_eta2*( &
             real(self%markers_list(i)%offset_y + self%markers_list(i)%i_cell_y - 1, f64)      )

    else if( i >= self%number_flow_markers + 1 .and. i <= self%number_flow_markers + self%number_moving_deposition_particles )then
      ! then the particle is a deposition particle

      r(1) = self%deposition_particles_eta(i - self%number_flow_markers, 1)
      r(2) = self%deposition_particles_eta(i - self%number_flow_markers, 2)

    end if

    r(3) = 0.0_f64

  end function bsl_lt_pic_4d_get_x


  !----------------------------------------------------------------------------------------------------------------------------
  ! gets the velocity coordinates of a 'particle', which can be of two types:
  !   1) either a flow marker (for i = 1, ... self%number_flow_markers)
  !   2) or a deposition particle (for i = self%number_flow_markers+1, ... self%number_moving_deposition_particles)

  pure function bsl_lt_pic_4d_get_v( self, i ) result( r )
    class( sll_bsl_lt_pic_4d_group ), intent( in ) :: self
    sll_int32                       , intent( in ) :: i
    sll_real64 :: r(3)

    if( i < 1 .or. i > self%number_flow_markers + self%number_moving_deposition_particles )then
      ! returning an error would be nice...
      r(1000) = 0   ! will this be caught as an error ?
      return
    end if

    if( i >= 1 .and. i <= self%number_flow_markers )then
      ! then the particle is a flow marker

      ! get vx
      r(1) = self%markers_list(i)%vx
      ! get vy
      r(2) = self%markers_list(i)%vy

    else if( i >= self%number_flow_markers + 1 .and. i <= self%number_flow_markers + self%number_moving_deposition_particles )then
      ! then the particle is a deposition particle

      r(1) = self%deposition_particles_eta(i - self%number_flow_markers, 3)
      r(2) = self%deposition_particles_eta(i - self%number_flow_markers, 4)

    end if

    r(3) = 0.0_f64

  end function bsl_lt_pic_4d_get_v


  !----------------------------------------------------------------------------------------------------------------------------
  ! get the cartesian cell index of a 'particle', which can be of two types:
  !   1) either a flow marker (for i = 1, ... self%number_flow_markers)
  !   2) or a deposition particle (for i = self%number_flow_markers+1, ... self%number_moving_deposition_particles)
  !
  ! note: (here i_out to match the abstract interface), 1 <= i_out <= num_cells_x * num_cells_y
  !
  ! note: almost the same function in the simple_pic_4d group -- maybe use the same function?

  pure function bsl_lt_pic_4d_get_cell_index(self, i) result(i_out)
    class(sll_bsl_lt_pic_4d_group),  intent( in )   ::  self
    sll_int32,                      intent( in )    ::  i       !> particle index
    sll_int32                                       ::  i_out   !> cell index
    sll_int32  ::  i_cell_x, i_cell_y
    sll_int32  ::  num_cells_x, num_cells_y
    sll_real64 ::  x_part
    sll_real64 ::  y_part
    sll_real64 ::  tmp !, dx, dy

    if( i < 1 .or. i > self%number_flow_markers + self%number_moving_deposition_particles )then
      ! returning an error would be nice...
      i_out = -1000   ! will this be caught as an error ?
      return
    end if

    if( i >= 1 .and. i <= self%number_flow_markers )then
      ! then the particle is a flow marker

      i_cell_x    = self%markers_list(i)%i_cell_x
      i_cell_y    = self%markers_list(i)%i_cell_y

    else if( i >= self%number_flow_markers + 1 .and. i <= self%number_flow_markers + self%number_moving_deposition_particles )then
      ! then the particle is a deposition particle

      ! find poisson (x-)cell containing this deposition particle, and relative position in the cell
      x_part = self%deposition_particles_eta(i - self%number_flow_markers, 1)
      tmp = ( x_part - self%space_mesh_2d%eta1_min ) / self%space_mesh_2d%delta_eta1
      i_cell_x = int( tmp ) + 1
      !! NOTE: if we need the relative position within the cell (between 0 and 1), it is: dx = tmp - (i_cell_x - 1)

      ! find poisson (y-)cell containing this node, seen as a deposition particle, and relative position in the cell
      y_part = self%deposition_particles_eta(i - self%number_flow_markers, 2)
      tmp = ( y_part - self%space_mesh_2d%eta2_min ) / self%space_mesh_2d%delta_eta2
      i_cell_y = int( tmp ) + 1
      !! NOTE: if we need the relative position within the cell (between 0 and 1), it is: dy = tmp - (i_cell_y - 1)

    end if

    ! then get the integer index for this cell
    num_cells_x = self%space_mesh_2d%num_cells1
    num_cells_y = self%space_mesh_2d%num_cells2
    i_out = 1 + modulo(i_cell_x - 1,  num_cells_x) + modulo(i_cell_y - 1,  num_cells_y) * num_cells_x   ! often denoted i_cell

  end function bsl_lt_pic_4d_get_cell_index



  !----------------------------------------------------------------------------------------------------------------------------
  ! sets the physical coordinates of a 'particle', which can be of two types:
  !   1) either a flow marker (for i = 1, ... self%number_flow_markers)
  !   2) or a deposition particle (for i = self%number_flow_markers+1, ... self%number_moving_deposition_particles)
  !
  ! In the first case, because of the data structure we need to transform a standard particle position (x,y) in
  ! (i_cell_x, i_cell_y, offset_x, offset_y).
  ! -> here the indices i_cell_x and i_cell_y do not need to be within [1, mesh%num_cells1] or [1, mesh%num_cells2]
  !    so that: - in periodic domains, the flows are better represented (no information is lost using modulo)
  !             - in non-periodic domains we can track outside particles (markers)
  !
  ! note: the integer index of the physical cell (used eg for the Poisson solver) is then obtained with get_poisson_cell_index
  !
  ! same function in the simple_pic_4d group -- (possible to use the same function?)

  subroutine bsl_lt_pic_4d_set_x( self, i, x )
    class( sll_bsl_lt_pic_4d_group ), intent( inout ) :: self
    sll_int32                       , intent( in    ) :: i
    sll_real64                      , intent( in    ) :: x(3)

    type(sll_cartesian_mesh_2d),      pointer :: space_mesh_2d
    type(sll_bsl_lt_pic_4d_particle), pointer :: particle
    sll_int32               :: i_cell_x, i_cell_y
    sll_real32              :: offset_x, offset_y
    sll_real64              :: temp

    if( i < 1 .or. i > self%number_flow_markers + self%number_moving_deposition_particles )then
      ! returning an error would be nice...
      temp = x(1000)   ! will this be caught as an error ?
      return
    end if

    if( i >= 1 .and. i <= self%number_flow_markers )then
      ! then we set the physical coordinates of a flow marker
      ! (maybe change the name of the structure sll_bsl_lt_pic_4d_particle -> sll_bsl_lt_pic_4d_flow_marker ?)

      space_mesh_2d => self%space_mesh_2d
      particle => self%markers_list(i)

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

    else if( i >= self%number_flow_markers + 1 .and. i <= self%number_flow_markers + self%number_moving_deposition_particles )then
      ! then we set the physical coordinates of a deposition particle
      ! (maybe use a structure closer to the simple_pic particles, cell-based ?)

      self%deposition_particles_eta(i - self%number_flow_markers, 1) = x(1)
      self%deposition_particles_eta(i - self%number_flow_markers, 2) = x(2)

    end if


  end subroutine bsl_lt_pic_4d_set_x


  !----------------------------------------------------------------------------
  ! sets the velocity coordinates of a 'particle', which can be of two types:
  !   1) either a flow marker (for i = 1, ... self%number_flow_markers)
  !   2) or a deposition particle (for i = self%number_flow_markers+1, ... self%number_moving_deposition_particles)

 subroutine bsl_lt_pic_4d_set_v( self, i, x )
    class( sll_bsl_lt_pic_4d_group ), intent( inout ) :: self
    sll_int32                       , intent( in    ) :: i
    sll_real64                      , intent( in    ) :: x(3)  !> this is the velocity, but argument name in abstract interface is x

    type(sll_bsl_lt_pic_4d_particle), pointer :: particle

    if( i < 1 .or. i > self%number_flow_markers + self%number_moving_deposition_particles )then
      ! returning an error would be nice...
      particle => self%markers_list(i)   ! will this be caught as an error ?
      return
    end if

    if( i >= 1 .and. i <= self%number_flow_markers )then
      ! then we set the velocity coordinates of a flow marker
      ! (maybe change the name of the structure sll_bsl_lt_pic_4d_particle -> sll_bsl_lt_pic_4d_flow_marker ?)

      particle => self%markers_list(i)
      particle%vx = x(1)
      particle%vy = x(2)

    else if( i >= self%number_flow_markers + 1 .and. i <= self%number_flow_markers + self%number_moving_deposition_particles )then
      ! then we set the physical coordinates of a deposition particle

      self%deposition_particles_eta(i - self%number_flow_markers, 3) = x(1)
      self%deposition_particles_eta(i - self%number_flow_markers, 4) = x(2)

    end if

  end subroutine bsl_lt_pic_4d_set_v


  !----------------------------------------------------------------------------
  subroutine bsl_lt_pic_4d_set_common_weight( self, s )
    class( sll_bsl_lt_pic_4d_group ), intent( inout ) :: self
    sll_real64                      , intent( in    ) :: s

    print*, "Error (9O8657864) -- this subroutine is not implemented for sll_bsl_lt_pic_4d_group objects", s, storage_size(self)
    stop

  end subroutine bsl_lt_pic_4d_set_common_weight


  !----------------------------------------------------------------------------
  subroutine bsl_lt_pic_4d_set_particle_weight( self, i, s )
    class( sll_bsl_lt_pic_4d_group ), intent( inout ) :: self
    sll_int32                       , intent( in    ) :: i
    sll_real64                      , intent( in    ) :: s

    print*, "Error (97658758) -- this subroutine is not implemented for sll_bsl_lt_pic_4d_group objects", s, storage_size(self)
    stop

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
  subroutine bsl_lt_pic_4d_set_hat_f0_parameters( self, x0, y0, vx0, vy0, r_x, r_y, r_vx, r_vy, basis_height, shift )
    class(sll_bsl_lt_pic_4d_group), intent(inout)   :: self
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

    self%hat_f0_x0 = x0
    self%hat_f0_y0 = y0
    self%hat_f0_vx0 = vx0
    self%hat_f0_vy0 = vy0
    self%hat_f0_r_x = r_x
    self%hat_f0_r_y = r_y
    self%hat_f0_r_vx = r_vx
    self%hat_f0_r_vy = r_vy
    self%hat_f0_basis_height = basis_height
    self%hat_f0_shift = shift

  end subroutine bsl_lt_pic_4d_set_hat_f0_parameters


  !--------------------------------------------------------------------------------------------------------------------------
  !> This subroutine places the deposition particles with a quasi-random (Sobol) sequence. It does not compute the weights.
  !> The weights are computed in an external call to the 'write_f_on_grid_or_deposit' subroutine
  subroutine bsl_lt_pic_set_deposition_particles_coordinates(self, rank)
    class(sll_bsl_lt_pic_4d_group), intent(inout)   :: self
    sll_int32, intent(in), optional                 :: rank
    sll_int32                                       :: this_rank
    sll_int32                                       :: i_part
    sll_int32                                       :: i_dim
    sll_int32                                       :: ierr
    sll_int64                                       :: sobol_seed
    sll_real64                                      :: rdn(4)

    SLL_ASSERT( self%deposition_particles_type == SLL_BSL_LT_PIC_TRANSPORTED_RANDOM )

    ! initial value of the seed (it is incremented by one in the call to i8_sobol)
    if(present(rank))then
      this_rank = rank
    else
      this_rank = 0
    end if
    sobol_seed = 10 + this_rank * self%number_moving_deposition_particles

    do i_part = 1, self%number_moving_deposition_particles
      ! Generate Sobol numbers on [0,1]
      call i8_sobol(int(4,8), sobol_seed, rdn)

      ! Transform rdn to the proper intervals
      do i_dim = 1, 4
        self%deposition_particles_eta(i_part, i_dim) =   self%remapping_grid_eta_min(i_dim) * (1 - rdn(i_dim)) &
                                                       + self%remapping_grid_eta_max(i_dim) * rdn(i_dim)
      end do
    end do

  end subroutine


  !----------------------------------------------------------------------------------
  ! do not change the position of the deposition particles, but compute (and set)
  ! their weights using the BSL_LT_PIC reconstruction

  subroutine bsl_lt_pic_set_deposition_particles_weights(self)
    class( sll_bsl_lt_pic_4d_group ),           intent( inout ) :: self
    type(sll_charge_accumulator_2d),  pointer :: void_charge_accumulator
    type(sll_cartesian_mesh_4d),      pointer :: void_grid_4d
    sll_real64, dimension(:,:),       pointer :: void_array_2d

    sll_real64    :: dummy_total_charge
    logical       :: enforce_total_charge
    sll_int32     :: scenario

    nullify(void_charge_accumulator)
    nullify(void_grid_4d)
    nullify(void_array_2d)

    scenario = SLL_BSL_LT_PIC_SET_WEIGHTS_ON_DEPOSITION_PARTICLES
    dummy_total_charge = 0.0_f64
    enforce_total_charge = .false.

    SLL_ASSERT( self%deposition_particles_type == SLL_BSL_LT_PIC_TRANSPORTED_RANDOM )

    call self%bsl_lt_pic_4d_write_f_on_grid_or_deposit(void_charge_accumulator,             &
                                                          scenario,                         &
                                                          void_grid_4d,                     &
                                                          void_array_2d,                    &
                                                          dummy_total_charge,               &
                                                          enforce_total_charge              &
                                                          )

  end subroutine


  !----------------------------------------------------------------------------
  ! deposit charge carried by the bsl_lt_pic_4d particles on a 2d mesh

  subroutine bsl_lt_pic_4d_deposit_charge_2d( self, charge_accumulator, target_total_charge, enforce_total_charge)
    class( sll_bsl_lt_pic_4d_group ),           intent( inout ) :: self
    type( sll_charge_accumulator_2d ), pointer, intent( inout ) :: charge_accumulator
    sll_real64,                                 intent(in)      :: target_total_charge
    logical,                                    intent(in)      :: enforce_total_charge

    type(charge_accumulator_cell_2d), pointer :: charge_accumulator_cell
    type(sll_cartesian_mesh_4d),     pointer  :: void_grid_4d        ! make this argument optional ?
    sll_real64, dimension(:,:),      pointer  :: void_array_2d       ! make this argument optional ?
    sll_int32     :: scenario

    sll_int32  :: i_cell_x
    sll_int32  :: i_cell_y
    sll_int32  :: i_cell
    sll_int32  :: i_part
    sll_real64 :: particle_charge
    sll_real64 :: dx
    sll_real64 :: dy
    sll_real64 :: tmp
    sll_real64 :: x_part
    sll_real64 :: y_part


    if( self%deposition_particles_type == SLL_BSL_LT_PIC_FIXED_GRID )then
      ! then we compute new weights on the deposition grid, as follows
      nullify(void_grid_4d)
      nullify(void_array_2d)
      scenario = SLL_BSL_LT_PIC_DEPOSIT_F
      call reset_charge_accumulator_to_zero ( charge_accumulator )

      call self%bsl_lt_pic_4d_write_f_on_grid_or_deposit(charge_accumulator,                                &
                                                         scenario,                                          &
                                                         void_grid_4d,                                      &
                                                         void_array_2d,                                     &
                                                         target_total_charge,                               &
                                                         enforce_total_charge                               &
                                                         )

    else
      ! then we simply deposit the charge carried by the deposition particles

      SLL_ASSERT( self%deposition_particles_type == SLL_BSL_LT_PIC_TRANSPORTED_RANDOM )
      do i_part = 1, self%number_moving_deposition_particles

        particle_charge = self%deposition_particles_weight(i_part)
        x_part = self%deposition_particles_eta(i_part, 1)
        y_part = self%deposition_particles_eta(i_part, 2)

        ! todo: use the interface function for the computation of the Poisson cell index? (but we need to return the cell_offset)

        ! find poisson (x-)cell containing this deposition particle, and relative position in the cell
        tmp = ( x_part - self%space_mesh_2d%eta1_min ) / self%space_mesh_2d%delta_eta1
        i_cell_x = int( tmp ) + 1
        dx = tmp - (i_cell_x - 1)  ! x-offset in cell (between 0 and 1)

        ! find poisson (y-)cell containing this node, seen as a deposition particle, and relative position in the cell
        tmp = ( y_part - self%space_mesh_2d%eta2_min ) / self%space_mesh_2d%delta_eta2
        i_cell_y = int( tmp ) + 1
        dy = tmp - (i_cell_y - 1)  ! y-offset in cell (between 0 and 1)

        ! set the proper accumulator cell for the deposition
        i_cell = i_cell_x + (i_cell_y - 1) * self%space_mesh_2d%num_cells1   !  (see global_to_cell_offset)

        charge_accumulator_cell => charge_accumulator%q_acc(i_cell)

        if( x_is_in_domain_2d(    x_part, y_part,                 &
                                  self%space_mesh_2d,             &
                                  self%domain_is_periodic(1),     &
                                  self%domain_is_periodic(2) ))then

          charge_accumulator_cell%q_sw = charge_accumulator_cell%q_sw + particle_charge * (1.0_f64 - dx) * (1.0_f64 - dy)
          charge_accumulator_cell%q_se = charge_accumulator_cell%q_se + particle_charge *            dx  * (1.0_f64 - dy)
          charge_accumulator_cell%q_nw = charge_accumulator_cell%q_nw + particle_charge * (1.0_f64 - dx) *            dy
          charge_accumulator_cell%q_ne = charge_accumulator_cell%q_ne + particle_charge *            dx  *            dy
          ! counter ???    deposited_charge = deposited_charge + particle_charge

        else
          ! particle not in domain (should store the reference for later processing)
          print*, "Error (09864543254786875): for the moment every particle should be in the (periodic) 2d domain..."
          stop
        end if

      end do

    end if
  end subroutine bsl_lt_pic_4d_deposit_charge_2d


  !> bsl_lt_pic_4d_visualize_f_slice_x_vx  plots an approximation of  f_x_vx = \int \int f(x,y,v_x,v_y) d y d v_y
  !>   - the plot is done on a 2d grid, but uses a 4d grid to evaluate f
  !>   - grid dimensions: we give the number of points, and the boundaries are given by the remapping domain
  !>   - calls sll_gnuplot_2d to write the data file

  subroutine bsl_lt_pic_4d_visualize_f_slice_x_vx(self, array_name, plot_np_x, plot_np_y, plot_np_vx, plot_np_vy, iplot)

    class( sll_bsl_lt_pic_4d_group ),   intent( inout ) :: self
    character(len=*),                   intent(in)      :: array_name   !< field name
    sll_int32,                          intent(in)      :: plot_np_x    !< nb of points in the x  plotting grid (see comment above)
    sll_int32,                          intent(in)      :: plot_np_y    !< nb of points in the y  plotting grid (see comment above)
    sll_int32,                          intent(in)      :: plot_np_vx   !< nb of points in the vx plotting grid (see comment above)
    sll_int32,                          intent(in)      :: plot_np_vy   !< nb of points in the vy plotting grid (see comment above)
    sll_int32,                          intent(in)      :: iplot        !< plot counter

    sll_int32 :: ierr

    sll_real64, dimension(:,:),       pointer :: x_vx_grid_values
    type(sll_cartesian_mesh_4d),      pointer :: plotting_grid_4d
    type(sll_charge_accumulator_2d),  pointer :: void_charge_accumulator
    sll_int32     :: scenario
    sll_real64    :: dummy_total_charge
    logical       :: enforce_total_charge

    nullify(void_charge_accumulator)

    !    print *, " plot A"

    SLL_ALLOCATE( x_vx_grid_values(plot_np_x, plot_np_vx), ierr)

    plotting_grid_4d => new_cartesian_mesh_4d(plot_np_x  - 1,                 &
                                              plot_np_y  - 1,                 &
                                              plot_np_vx - 1,                 &
                                              plot_np_vy - 1,                 &
                                              self%remapping_grid_eta_min(1), &
                                              self%remapping_grid_eta_max(1), &
                                              self%remapping_grid_eta_min(2), &
                                              self%remapping_grid_eta_max(2), &
                                              self%remapping_grid_eta_min(3), &
                                              self%remapping_grid_eta_max(3), &
                                              self%remapping_grid_eta_min(4), &
                                              self%remapping_grid_eta_max(4)  &
                                            )

    scenario = SLL_BSL_LT_PIC_WRITE_F_ON_GIVEN_GRID
    dummy_total_charge = 0.0_f64
    enforce_total_charge = .false.

    call self%bsl_lt_pic_4d_write_f_on_grid_or_deposit(void_charge_accumulator,   &
                                                       scenario,                  &
                                                       plotting_grid_4d,          &
                                                       x_vx_grid_values,          &
                                                       dummy_total_charge,        &
                                                       enforce_total_charge       &
                                                       )

    ! print *, "plot T"

    call sll_gnuplot_2d(self%remapping_grid_eta_min(1), &
                        self%remapping_grid_eta_max(1), &
                        plot_np_x,                      &     ! (note: this is indeed the nb of plotted points, not 'cells')
                        self%remapping_grid_eta_min(3), &
                        self%remapping_grid_eta_max(3), &
                        plot_np_vx,                     &     ! (same comment)
                        x_vx_grid_values,               &
                        array_name,                     &
                        iplot,                          &
                        ierr )

  end subroutine bsl_lt_pic_4d_visualize_f_slice_x_vx

  function remapping_cart_grid_number_nodes_x(p_group) result(val)
    class(sll_bsl_lt_pic_4d_group), intent(in)  :: p_group
    sll_int32                                  :: val

    if( p_group%domain_is_periodic(1) )then
        val = p_group%remapping_cart_grid_number_cells_x
    else
        val = p_group%remapping_cart_grid_number_cells_x + 1
    end if
  end function

  function remapping_cart_grid_number_nodes_y(p_group) result(val)
    class(sll_bsl_lt_pic_4d_group), intent(in)  :: p_group
    sll_int32                                  :: val

    if( p_group%domain_is_periodic(2) )then
        val = p_group%remapping_cart_grid_number_cells_y
    else
        val = p_group%remapping_cart_grid_number_cells_y + 1
    end if
  end function

  function remapping_cart_grid_number_nodes_vx(p_group) result(val)
    class(sll_bsl_lt_pic_4d_group), intent(in)  :: p_group
    sll_int32                                  :: val

    val = p_group%remapping_cart_grid_number_cells_vx + 1
  end function

  function remapping_cart_grid_number_nodes_vy(p_group) result(val)
    class(sll_bsl_lt_pic_4d_group), intent(in)  :: p_group
    sll_int32                                  :: val

    val = p_group%remapping_cart_grid_number_cells_vy + 1
  end function

  !----------------------------------------------------------------------------
  ! Constructor
  !> @brief Constructor for a group of bsl_lt_pic_4d particles
  function sll_bsl_lt_pic_4d_group_new(             &
        species_charge,                             &
        species_mass,                               &
        particle_group_id,                          &
        domain_is_x_periodic,                       &
        domain_is_y_periodic,                       &
        remap_f_type,                               &
        remap_degree,                               &
        remapping_grid_vx_min,                      &
        remapping_grid_vx_max,                      &
        remapping_grid_vy_min,                      &
        remapping_grid_vy_max,                      &
        remapping_cart_grid_number_cells_x,         &   ! for splines
        remapping_cart_grid_number_cells_y,         &   ! for splines
        remapping_cart_grid_number_cells_vx,        &   ! for splines
        remapping_cart_grid_number_cells_vy,        &   ! for splines
        remapping_sparse_grid_max_levels,           &   ! for the sparse grid: for now, same level in each dimension
        deposition_particles_type,                  &
        minimum_number_of_deposition_particles,     &   ! lower bound if deposition particles = fixed_grid, actual nb if transported
        number_flow_markers_x,                      &
        number_flow_markers_y,                      &
        number_flow_markers_vx,                     &
        number_flow_markers_vy,                     &
        flow_grid_number_cells_x,                   &
        flow_grid_number_cells_y,                   &
        flow_grid_number_cells_vx,                  &
        flow_grid_number_cells_vy,                  &
        space_mesh_2d ) result(res)

    type( sll_bsl_lt_pic_4d_group ), pointer :: res

    sll_real64,               intent(in)  :: species_charge
    sll_real64,               intent(in)  :: species_mass
    sll_int32,                intent(in)  :: particle_group_id
    logical,                  intent(in)  :: domain_is_x_periodic
    logical,                  intent(in)  :: domain_is_y_periodic

    sll_int32,                intent(in)  :: remap_f_type
    sll_int32,                intent(in)  :: remap_degree
    sll_real64,               intent(in)  :: remapping_grid_vx_min
    sll_real64,               intent(in)  :: remapping_grid_vx_max
    sll_real64,               intent(in)  :: remapping_grid_vy_min
    sll_real64,               intent(in)  :: remapping_grid_vy_max
    sll_int32,                intent(in)  :: remapping_cart_grid_number_cells_x
    sll_int32,                intent(in)  :: remapping_cart_grid_number_cells_y
    sll_int32,                intent(in)  :: remapping_cart_grid_number_cells_vx
    sll_int32,                intent(in)  :: remapping_cart_grid_number_cells_vy
    sll_int32,  dimension(4), intent(in)  :: remapping_sparse_grid_max_levels
    sll_int32,                intent(in)  :: deposition_particles_type
    sll_int32,                intent(in)  :: minimum_number_of_deposition_particles
    sll_int32,                intent(in)  :: number_flow_markers_x
    sll_int32,                intent(in)  :: number_flow_markers_y
    sll_int32,                intent(in)  :: number_flow_markers_vx
    sll_int32,                intent(in)  :: number_flow_markers_vy
    sll_int32,                intent(in)  :: flow_grid_number_cells_x
    sll_int32,                intent(in)  :: flow_grid_number_cells_y
    sll_int32,                intent(in)  :: flow_grid_number_cells_vx
    sll_int32,                intent(in)  :: flow_grid_number_cells_vy

    type(sll_cartesian_mesh_2d), pointer, intent(in) :: space_mesh_2d

    sll_int32               :: number_cells_initial_markers_grid_x
    sll_int32               :: number_cells_initial_markers_grid_y
    sll_int32               :: number_cells_initial_markers_grid_vx
    sll_int32               :: number_cells_initial_markers_grid_vy

    sll_int32               :: ierr
    character(len=*), parameter :: this_fun_name = "sll_bsl_lt_pic_4d_group_new"
    character(len=128)      :: err_msg

    sll_int32  :: nb_deposition_particles_x
    sll_int32  :: nb_deposition_particles_y
    sll_int32  :: nb_deposition_particles_vx
    sll_int32  :: nb_deposition_particles_vy

    sll_int32  :: deposition_grid_num_cells_x
    sll_int32  :: deposition_grid_num_cells_y
    sll_int32  :: deposition_grid_num_cells_vx
    sll_int32  :: deposition_grid_num_cells_vy

    sll_real64 :: deposition_grid_x_min
    sll_real64 :: deposition_grid_x_max
    sll_real64 :: deposition_grid_y_min
    sll_real64 :: deposition_grid_y_max
    sll_real64 :: deposition_grid_vx_min
    sll_real64 :: deposition_grid_vx_max
    sll_real64 :: deposition_grid_vy_min
    sll_real64 :: deposition_grid_vy_max

    sll_real64 :: h_deposition_grid_x
    sll_real64 :: h_deposition_grid_y

    sll_int32  :: cst_int
    sll_real64 :: cst_real, ratio_vx, ratio_vy
    sll_real64 :: tmp, tmp1, tmp2

    if (.not.associated(space_mesh_2d) ) then
       err_msg = 'Error: given space_mesh_2d is not associated'
       SLL_ERROR( this_fun_name, err_msg )
    end if

    SLL_ALLOCATE( res, ierr )

    !> create the species object for this particle group
    res%species => temp_species_new( species_charge, species_mass )

    res%id = particle_group_id
    res%dimension_x = 2
    res%dimension_v = 2

    !> A. discretization of the flow:
    !>    - A.1 list of marker coordinates (pushed forward)
    !>    - A.2 cartesian grid of initial markers
    !>    - A.3 flow grid: 4d cartesian cells where the flow is linearized

    !> A.1 list of marker coordinates
    res%number_flow_markers_x  = number_flow_markers_x
    res%number_flow_markers_y  = number_flow_markers_y
    res%number_flow_markers_vx = number_flow_markers_vx
    res%number_flow_markers_vy = number_flow_markers_vy
    res%number_flow_markers    = number_flow_markers_x * number_flow_markers_y * number_flow_markers_vx * number_flow_markers_vy

    SLL_ALLOCATE( res%markers_list(res%number_flow_markers), ierr )

    !> A.2 cartesian grid of initial markers
    if( domain_is_x_periodic )then
        number_cells_initial_markers_grid_x = number_flow_markers_x
    else
        number_cells_initial_markers_grid_x = number_flow_markers_x - 1
    end if
    if( domain_is_y_periodic )then
        number_cells_initial_markers_grid_y = number_flow_markers_y
    else
        number_cells_initial_markers_grid_y = number_flow_markers_y - 1
    end if
    number_cells_initial_markers_grid_vx = number_flow_markers_vx - 1
    number_cells_initial_markers_grid_vy = number_flow_markers_vy - 1

    res%initial_markers_grid => new_cartesian_mesh_4d( &
      number_cells_initial_markers_grid_x, &
      number_cells_initial_markers_grid_y, &
      number_cells_initial_markers_grid_vx, &
      number_cells_initial_markers_grid_vy, &
      space_mesh_2d%eta1_min, &
      space_mesh_2d%eta1_max, &
      space_mesh_2d%eta2_min, &
      space_mesh_2d%eta2_max, &
      remapping_grid_vx_min, &
      remapping_grid_vx_max, &
      remapping_grid_vy_min, &
      remapping_grid_vy_max )


    !> A.2 flow grid
    res%flow_grid => new_cartesian_mesh_4d( &
      flow_grid_number_cells_x, &
      flow_grid_number_cells_y, &
      flow_grid_number_cells_vx, &
      flow_grid_number_cells_vy, &
      space_mesh_2d%eta1_min, &
      space_mesh_2d%eta1_max, &
      space_mesh_2d%eta2_min, &
      space_mesh_2d%eta2_max, &
      remapping_grid_vx_min, &
      remapping_grid_vx_max, &
      remapping_grid_vy_min, &
      remapping_grid_vy_max )


    !> B. discretization of the field (Poisson grid)
    !> physical 2d mesh (used eg in the Poisson solver)

    res%space_mesh_2d => space_mesh_2d
    res%domain_is_periodic(1) = domain_is_x_periodic
    res%domain_is_periodic(2) = domain_is_y_periodic


    !> C. discretization of the remapped f:
    !>    - C.0 size of the remapping grid used in the interpolator for the remapped_f (for splines or sparse grid)
    !>    - C.1 interpolator for splines
    !>      C.1.a  cartesian remapping grid
    !>      C.1.b  array of interpolation coefficients for remapped_f
    !>    - C.2 interpolator for sparse grids
    !>      C.2.a  sparse remapping grid
    !>      C.2.b  array of interpolation coefficients for remapped_f

    !> C.0 size of the remapping grid
    res%remapping_grid_eta_min(1) = space_mesh_2d%eta1_min
    res%remapping_grid_eta_max(1) = space_mesh_2d%eta1_max
    res%remapping_grid_eta_min(2) = space_mesh_2d%eta2_min
    res%remapping_grid_eta_max(2) = space_mesh_2d%eta2_max
    res%remapping_grid_eta_min(3) = remapping_grid_vx_min
    res%remapping_grid_eta_max(3) = remapping_grid_vx_max
    res%remapping_grid_eta_min(4) = remapping_grid_vy_min
    res%remapping_grid_eta_max(4) = remapping_grid_vy_max

    res%remapped_f_interpolation_type = remap_f_type
    res%remapped_f_interpolation_degree = remap_degree

    if( res%remapped_f_interpolation_type == SLL_BSL_LT_PIC_REMAP_WITH_SPLINES )then

      ! C.1 interpolator for splines

      SLL_ASSERT( remap_degree==1 )

      ! C.1.a  cartesian remapping grid
      res%remapping_cart_grid_number_cells_x  = remapping_cart_grid_number_cells_x
      res%remapping_cart_grid_number_cells_y  = remapping_cart_grid_number_cells_y
      res%remapping_cart_grid_number_cells_vx = remapping_cart_grid_number_cells_vx
      res%remapping_cart_grid_number_cells_vy = remapping_cart_grid_number_cells_vy

      res%remapping_cart_grid => new_cartesian_mesh_4d(       &
                                                  remapping_cart_grid_number_cells_x,        &
                                                  remapping_cart_grid_number_cells_y,        &
                                                  remapping_cart_grid_number_cells_vx,       &
                                                  remapping_cart_grid_number_cells_vy,       &
                                                  res%remapping_grid_eta_min(1),   &
                                                  res%remapping_grid_eta_max(1),   &
                                                  res%remapping_grid_eta_min(2),   &
                                                  res%remapping_grid_eta_max(2),   &
                                                  res%remapping_grid_eta_min(3),  &
                                                  res%remapping_grid_eta_max(3),  &
                                                  res%remapping_grid_eta_min(4),  &
                                                  res%remapping_grid_eta_max(4)   &
                                                )

      !> store the quasi-interpolation coefficients used in the remappings and initialisation, in the case of cubic spline particles
      if( res%remapped_f_interpolation_degree == 3) then
          SLL_ALLOCATE( res%lt_pic_interpolation_coefs(-1:1), ierr )
             res%lt_pic_interpolation_coefs(-1) = -1.0_f64/6.0_f64
             res%lt_pic_interpolation_coefs(0)  =  8.0_f64/6.0_f64
             res%lt_pic_interpolation_coefs(1)  = -1.0_f64/6.0_f64
      end if


      ! C.1.b  array of spline coefficients for remapped_f
      allocate(res%remapped_f_splines_coefficients( res%remapping_cart_grid_number_nodes_x(),     &
                                                    res%remapping_cart_grid_number_nodes_y(),     &
                                                    res%remapping_cart_grid_number_nodes_vx(),    &
                                                    res%remapping_cart_grid_number_nodes_vy() )   &
               , stat=ierr)
      call test_error_code(ierr, 'Memory allocation Failure.', __FILE__, __LINE__)


    else if( res%remapped_f_interpolation_type == SLL_BSL_LT_PIC_REMAP_WITH_SPARSE_GRIDS )then
      ! C.2 interpolator for sparse grids

      ! C.2.a  sparse remapping grid
      print*, "sll_bsl_lt_pic_4d_group_new - C.2", remapping_sparse_grid_max_levels
      res%sparse_grid_max_levels(1) = remapping_sparse_grid_max_levels(1)
      res%sparse_grid_max_levels(2) = remapping_sparse_grid_max_levels(2)
      res%sparse_grid_max_levels(3) = remapping_sparse_grid_max_levels(3)
      res%sparse_grid_max_levels(4) = remapping_sparse_grid_max_levels(4)
      call res%sparse_grid_interpolator%initialize(         &
                  res%sparse_grid_max_levels,               &
                  res%remapped_f_interpolation_degree,      &
                  res%remapped_f_interpolation_degree+1,    &
                  0,                                        &  ! interpolation_type for the sparse grid (splines or Lagrange)
                  res%remapping_grid_eta_min,               &
                  res%remapping_grid_eta_max                &
                  )

      ! C.2.b  array of sparse grid coefficients for remapped_f
      SLL_ALLOCATE( res%remapped_f_sparse_grid_coefficients(res%sparse_grid_interpolator%size_basis), ierr )
      res%remapped_f_sparse_grid_coefficients = 0.0_f64




    !> D. discretization of the deposited f -- uses deposition particles which can either be:
    !     - D.1 on a fixed grid (with new charges computed at every time step)
    !     - D.2 or transported with the flow (like sdt particles), and re-initialized from time to time, using BSL_LT_PIC techniques

    res%deposition_particles_type = deposition_particles_type

    if( res%deposition_particles_type == SLL_BSL_LT_PIC_FIXED_GRID )then

      ! D.1 deposition particles on a fixed (cartesian) grid denoted 'deposition_grid', with the following properties:
      !
      !   - it matches the poisson grid in the sense that its nb of cells satisfies (with dg_np_d = deposition_grid_num_points_d)
      !     dg_np_x ~ cst_int * p_group%space_mesh_2d%num_cells1  and
      !     dg_np_y ~ cst_int * p_group%space_mesh_2d%num_cells2  for some  cst_int
      !
      !   - the griding in vx, vy is obtained from that of the initial markers, using the linear scaling
      !     dg_np_vx / (dg_np_x * dg_np_y) = (approx) number_flow_markers_vx / (number_flow_markers_x * number_flow_markers_y)
      !     dg_np_vy / (dg_np_x * dg_np_y) = (approx) number_flow_markers_vy / (number_flow_markers_x * number_flow_markers_y)
      !
      !   - its number of nodes satisfies
      !     dg_np = dg_np_x * dg_np_y * dg_nc_vx * dg_nc_vy >= minimum_number_of_deposition_particles

      ! we will have  dg_np_vx ~ cst_int * cst_int * ratio_vx
      ratio_vx = res%number_flow_markers_vx * 1./ (res%number_flow_markers_x * res%number_flow_markers_y)       &
                                           * res%space_mesh_2d%num_cells1 * res%space_mesh_2d%num_cells2
      ! and           dg_np_vy ~ cst_int * cst_int * ratio_vy
      ratio_vy = res%number_flow_markers_vy * 1./ (res%number_flow_markers_x * res%number_flow_markers_y)       &
                                           * res%space_mesh_2d%num_cells1 * res%space_mesh_2d%num_cells2

      ! and           dg_np ~ cst_int * cst_int * cst_int * cst_int * num_cells1 * num_cells2 * ratio_vx * ratio_vy
      !                     >=  minimum_number_of_deposition_particles

      ! cst_real is the float approx of cst_int above
      cst_real = (real(minimum_number_of_deposition_particles, f64)   &
             / real( ratio_vx * ratio_vy * res%space_mesh_2d%num_cells1 * res%space_mesh_2d%num_cells2, f64)) ** (1./6)

      cst_int = int(ceiling( cst_real ))


      nb_deposition_particles_x  = max( cst_int * res%space_mesh_2d%num_cells1, 2 )
      nb_deposition_particles_y  = max( cst_int * res%space_mesh_2d%num_cells2, 2 )
      nb_deposition_particles_vx = max( int(ceiling( cst_real * cst_real * ratio_vx )), 2 )
      nb_deposition_particles_vy = max( int(ceiling( cst_real * cst_real * ratio_vy )), 2 )

      deposition_grid_num_cells_x  = nb_deposition_particles_x  - 1
      deposition_grid_num_cells_y  = nb_deposition_particles_y  - 1
      deposition_grid_num_cells_vx = nb_deposition_particles_vx - 1
      deposition_grid_num_cells_vy = nb_deposition_particles_vy - 1

      tmp=real(nb_deposition_particles_x * nb_deposition_particles_y * nb_deposition_particles_vx * nb_deposition_particles_vy, f64)
      print*, "[", this_fun_name, "] will use ", tmp, "deposition particles"
      print*, "[bsl_lt_pic_4d_write_f_on_grid_or_deposit -- DEPOSIT_F] should be at least ", minimum_number_of_deposition_particles
      SLL_ASSERT( tmp >= minimum_number_of_deposition_particles )

      ! then we position the grid of deposition cells so that every deposition particle is _inside_ a poisson cell
      ! (so that we do not have to do something special for periodic boundary conditions)

      h_deposition_grid_x = res%space_mesh_2d%delta_eta1 / cst_int    ! distance between two depos. particles in x dimension
      h_deposition_grid_y = res%space_mesh_2d%delta_eta2 / cst_int    ! same in y
      deposition_grid_x_min  = res%space_mesh_2d%eta1_min + 0.5 * h_deposition_grid_x
      deposition_grid_x_max  = res%space_mesh_2d%eta1_max - 0.5 * h_deposition_grid_x
      deposition_grid_y_min  = res%space_mesh_2d%eta2_min + 0.5 * h_deposition_grid_y
      deposition_grid_y_max  = res%space_mesh_2d%eta2_max - 0.5 * h_deposition_grid_y

      ! in velocity the bounds are those of the remapping grid
      deposition_grid_vx_min = res%remapping_grid_eta_min(3)
      deposition_grid_vx_max = res%remapping_grid_eta_max(3)
      deposition_grid_vy_min = res%remapping_grid_eta_min(4)
      deposition_grid_vy_max = res%remapping_grid_eta_max(4)

      res%deposition_grid => new_cartesian_mesh_4d( deposition_grid_num_cells_x,        &
                                                    deposition_grid_num_cells_y,        &
                                                    deposition_grid_num_cells_vx,       &
                                                    deposition_grid_num_cells_vy,       &
                                                    deposition_grid_x_min,   &
                                                    deposition_grid_x_max,   &
                                                    deposition_grid_y_min,   &
                                                    deposition_grid_y_max,   &
                                                    deposition_grid_vx_min,  &
                                                    deposition_grid_vx_max,  &
                                                    deposition_grid_vy_min,  &
                                                    deposition_grid_vy_max   &
                                                   )

    else if( res%deposition_particles_type == SLL_BSL_LT_PIC_TRANSPORTED_RANDOM )then

      ! D.2 deposition particles will be transported with the flow (like sdt particles)
      !     and re-initialized from time to time, using BSL_LT_PIC techniques

      res%number_moving_deposition_particles = minimum_number_of_deposition_particles
      SLL_ALLOCATE( res%deposition_particles_eta(res%number_moving_deposition_particles, 4), ierr )
      SLL_ALLOCATE( res%deposition_particles_weight(res%number_moving_deposition_particles), ierr )

    else
       err_msg = 'Error: incorrect value for deposition_particles_type'
       SLL_ERROR( this_fun_name, err_msg )

    end if

    res%number_particles = res%number_flow_markers + res%number_moving_deposition_particles

    else
      SLL_ERROR("sll_bsl_lt_pic_4d_group_new", "ahem, a test must be broken -- you should not be reading this :)")
    end if

  end function sll_bsl_lt_pic_4d_group_new


  !> initializes the interpolation coefficients for f0 on the remapping grid, and the flow markers
  !> Note: since no interpolations are needed in the evaluation of f0, the arrays of interpolation coefficients can be used
  !> to store the nodal values of f0 (in the remapping routines this is not possible)
  subroutine bsl_lt_pic_4d_initializer( self, initial_density_identifier, rand_seed, rank, world_size )

    class( sll_bsl_lt_pic_4d_group ), intent( inout ) :: self
    sll_int32                       , intent( in    ) :: initial_density_identifier
    sll_int32, dimension(:)         , intent( in ), optional :: rand_seed
    sll_int32                       , intent( in ), optional :: rank, world_size

    !> A. initialize the remapping tool:
    !>    - A.1  write the nodal values of f0 on the arrays of interpolation coefs
    print *, "bsl_lt_pic_4d_initializer -- step A"
    if( initial_density_identifier == SLL_BSL_LT_PIC_LANDAU_F0 )then
      call self%bsl_lt_pic_4d_write_landau_density_on_remapping_grid( self%thermal_speed, self%alpha, self%k_landau )
    else if( initial_density_identifier == SLL_BSL_LT_PIC_HAT_F0 )then
      call self%bsl_lt_pic_4d_write_hat_density_on_remapping_grid(                              &
                      self%hat_f0_x0, self%hat_f0_y0, self%hat_f0_vx0, self%hat_f0_vy0,         &
                      self%hat_f0_r_x, self%hat_f0_r_y, self%hat_f0_r_vx, self%hat_f0_r_vy,     &
                                                                   self%hat_f0_basis_height, self%hat_f0_shift)
    else
      SLL_ERROR( "bsl_lt_pic_4d_initializer", "wrong value for initial_density_identifier" )
    end if

    !>    - A.2  compute the interpolation coefs for remapped_f (using the nodal values stored in the arrays of interpolation coefs)

    if( self%remapped_f_interpolation_type == SLL_BSL_LT_PIC_REMAP_WITH_SPLINES )then
      call self%bsl_lt_pic_4d_compute_new_spline_coefs()
    else if( self%remapped_f_interpolation_type == SLL_BSL_LT_PIC_REMAP_WITH_SPARSE_GRIDS )then
      call self%sparse_grid_interpolator%compute_hierarchical_surplus(      &
                self%remapped_f_sparse_grid_coefficients                    &
           )
    else
      SLL_ERROR( "bsl_lt_pic_4d_initializer", "wrong value for remapped_f_interpolation_type" )
    end if

    !> B. initialize the (flow) markers
    print *, "bsl_lt_pic_4d_initializer -- step B"

    call self%bsl_lt_pic_4d_reset_markers_position()
    call self%bsl_lt_pic_4d_set_markers_connectivity()

    !> D. if deposition particles are transported, initialize them -- this must be done after the remapping tool is operational!

    print *, "bsl_lt_pic_4d_initializer -- step C"

    if( self%deposition_particles_type == SLL_BSL_LT_PIC_TRANSPORTED_RANDOM )then
      if(present(rank))then
        call self%bsl_lt_pic_set_deposition_particles_coordinates(rank)
      else
        call self%bsl_lt_pic_set_deposition_particles_coordinates()
      end if
      call self%bsl_lt_pic_set_deposition_particles_weights()
    end if

    return

    !PN ADD TO PREVENT WARNING
    SLL_ASSERT(present(world_size))
    SLL_ASSERT(present(rand_seed))

  end subroutine bsl_lt_pic_4d_initializer


  subroutine bsl_lt_pic_4d_write_landau_density_on_remapping_grid(    &
              p_group,                              &
              thermal_speed, alpha, k_landau        &
              )

    class(sll_bsl_lt_pic_4d_group), intent(inout)    :: p_group
    sll_real64, intent(in)                          :: thermal_speed, alpha, k_landau

    sll_int32 :: j
    sll_int32 :: j_x
    sll_int32 :: j_y
    sll_int32 :: j_vx
    sll_int32 :: j_vy
    sll_int32 :: number_nodes_x
    sll_int32 :: number_nodes_y
    sll_int32 :: number_nodes_vx
    sll_int32 :: number_nodes_vy
    sll_real64 :: h_x
    sll_real64 :: h_y
    sll_real64 :: h_vx
    sll_real64 :: h_vy
    sll_real64 :: x_min
    sll_real64 :: y_min
    sll_real64 :: vx_min
    sll_real64 :: vy_min
    sll_real64 :: x_j
    sll_real64 :: vx_j
    sll_real64 :: vy_j

    sll_real64 :: f_x
    sll_real64 :: f_vx
    sll_real64 :: f_vy
    sll_real64 :: one_over_thermal_velocity
    sll_real64 :: one_over_two_pi

    sll_real64 :: total_density ! DEBUG

    one_over_thermal_velocity = 1./thermal_speed
    one_over_two_pi = 1./(2*sll_pi)


    ! print*, "bsl_lt_pic_4d_write_landau_density_on_remapping_grid -- a "

    if( p_group%remapped_f_interpolation_type == SLL_BSL_LT_PIC_REMAP_WITH_SPLINES )then

      number_nodes_x  = p_group%remapping_cart_grid_number_nodes_x()
      number_nodes_y  = p_group%remapping_cart_grid_number_nodes_y()
      number_nodes_vx = p_group%remapping_cart_grid_number_nodes_vx()
      number_nodes_vy = p_group%remapping_cart_grid_number_nodes_vy()

      h_x    = p_group%remapping_cart_grid%delta_eta1
      h_y    = p_group%remapping_cart_grid%delta_eta2
      h_vx   = p_group%remapping_cart_grid%delta_eta3
      h_vy   = p_group%remapping_cart_grid%delta_eta4

      x_min    = p_group%remapping_cart_grid%eta1_min
      y_min    = p_group%remapping_cart_grid%eta2_min
      vx_min   = p_group%remapping_cart_grid%eta3_min
      vy_min   = p_group%remapping_cart_grid%eta4_min

      ! compute the values of f0 on the (cartesian, phase-space) remapping grid

      SLL_ASSERT( size(p_group%remapped_f_splines_coefficients,1) == number_nodes_x )
      SLL_ASSERT( size(p_group%remapped_f_splines_coefficients,2) == number_nodes_y )
      SLL_ASSERT( size(p_group%remapped_f_splines_coefficients,3) == number_nodes_vx )
      SLL_ASSERT( size(p_group%remapped_f_splines_coefficients,4) == number_nodes_vy )

      ! print*, "bsl_lt_pic_4d_write_landau_density_on_remapping_grid -- assert ok "

      do j_x = 1, number_nodes_x
        x_j = x_min + (j_x-1) * h_x
        f_x = eval_landau_fx(alpha, k_landau, x_j)

        do j_y = 1, number_nodes_y
          ! (density does not depend on y)

          do j_vx = 1, number_nodes_vx
            vx_j = vx_min + (j_vx-1) * h_vx
            f_vx = one_over_thermal_velocity * exp(-0.5*(vx_j*one_over_thermal_velocity)**2)

            do j_vy = 1, number_nodes_vy
              vy_j = vy_min + (j_vy-1) * h_vy
              f_vy = one_over_thermal_velocity * exp(-0.5*(vy_j*one_over_thermal_velocity)**2)

              ! here we store a nodal value but later this array will indeed store spline coefficients
              p_group%remapped_f_splines_coefficients(j_x,j_y,j_vx,j_vy) = one_over_two_pi * f_x * f_vx * f_vy
              ! todo: use temp array (for name clarity?)

            end do
          end do
        end do
      end do

    else

      !      print*, "bsl_lt_pic_4d_write_landau_density_on_remapping_grid -- b"
      !      print*, p_group%remapped_f_interpolation_type
      !      print*, "bsl_lt_pic_4d_write_landau_density_on_remapping_grid -- c"
      !      print*, p_group%sparse_grid_interpolator%size_basis
      !      print*, size(p_group%remapped_f_sparse_grid_coefficients,1)
      !      print*, "bsl_lt_pic_4d_write_landau_density_on_remapping_grid -- d"
      SLL_ASSERT( p_group%remapped_f_interpolation_type == SLL_BSL_LT_PIC_REMAP_WITH_SPARSE_GRIDS )
      SLL_ASSERT( size(p_group%remapped_f_sparse_grid_coefficients,1) == p_group%sparse_grid_interpolator%size_basis )

      !      print*, "bsl_lt_pic_4d_write_landau_density_on_remapping_grid -- WARNING DEBUG 1908987987876768687687568765876986987 "
      !      print*, "bsl_lt_pic_4d_write_landau_density_on_remapping_grid -- WARNING --- INITIALIZING WITH 1 ------------------- "
      !      print*, "bsl_lt_pic_4d_write_landau_density_on_remapping_grid -- WARNING DEBUG 1908987987876768687687568765876986987 "
      !      print*, "bsl_lt_pic_4d_write_landau_density_on_remapping_grid -- WARNING --- INITIALIZING WITH 1 ------------------- "
      !      print*, "bsl_lt_pic_4d_write_landau_density_on_remapping_grid -- WARNING DEBUG 1908987987876768687687568765876986987 "
      !      print*, "bsl_lt_pic_4d_write_landau_density_on_remapping_grid -- WARNING --- INITIALIZING WITH 1 ------------------- "
      !      print*, "bsl_lt_pic_4d_write_landau_density_on_remapping_grid -- WARNING DEBUG 1908987987876768687687568765876986987 "
      !      print*, "bsl_lt_pic_4d_write_landau_density_on_remapping_grid -- WARNING --- INITIALIZING WITH 1 ------------------- "
      !      print*, "bsl_lt_pic_4d_write_landau_density_on_remapping_grid -- WARNING DEBUG 1908987987876768687687568765876986987 "
      !      print*, "bsl_lt_pic_4d_write_landau_density_on_remapping_grid -- WARNING --- INITIALIZING WITH 1 ------------------- "
      !      print*, "bsl_lt_pic_4d_write_landau_density_on_remapping_grid -- WARNING DEBUG 1908987987876768687687568765876986987 "
      !      print*, "bsl_lt_pic_4d_write_landau_density_on_remapping_grid -- WARNING --- INITIALIZING WITH 1 ------------------- "
      !      print*, "bsl_lt_pic_4d_write_landau_density_on_remapping_grid -- WARNING DEBUG 1908987987876768687687568765876986987 "
      !      print*, "bsl_lt_pic_4d_write_landau_density_on_remapping_grid -- WARNING --- INITIALIZING WITH 1 ------------------- "
      !      print*, "bsl_lt_pic_4d_write_landau_density_on_remapping_grid -- WARNING DEBUG 1908987987876768687687568765876986987 "
      !      print*, "bsl_lt_pic_4d_write_landau_density_on_remapping_grid -- WARNING --- INITIALIZING WITH 1 ------------------- "
      !      print*, "bsl_lt_pic_4d_write_landau_density_on_remapping_grid -- WARNING DEBUG 1908987987876768687687568765876986987 "
      !      print*, "bsl_lt_pic_4d_write_landau_density_on_remapping_grid -- WARNING --- INITIALIZING WITH 1 ------------------- "

      do j=1, p_group%sparse_grid_interpolator%size_basis
          x_j  = p_group%sparse_grid_interpolator%hierarchy(j)%coordinate(1)
          vx_j = p_group%sparse_grid_interpolator%hierarchy(j)%coordinate(3)
          vy_j = p_group%sparse_grid_interpolator%hierarchy(j)%coordinate(4)

          f_x = 1._f64 + alpha * cos(k_landau * x_j)  ! eval_landau_fx(alpha, k_landau, x_j)
          f_vx = one_over_thermal_velocity * exp(-0.5*(vx_j * one_over_thermal_velocity)**2)
          f_vy = one_over_thermal_velocity * exp(-0.5*(vy_j * one_over_thermal_velocity)**2)

          ! here we store a nodal value but later this array will indeed store sparse grid coefficients
          p_group%remapped_f_sparse_grid_coefficients(j) = one_over_two_pi * f_x * f_vx * f_vy

          !                 p_group%remapped_f_sparse_grid_coefficients(j) = 1.

      end do


      !      print*, "bsl_lt_pic_4d_write_landau_density_on_remapping_grid -- DEBUG 09890908987876988097986 "
      !      print*, "bsl_lt_pic_4d_write_landau_density_on_remapping_grid -- DEBUG 09890908987876988097986 "
      !      print*, "bsl_lt_pic_4d_write_landau_density_on_remapping_grid -- DEBUG 09890908987876988097986 "
      !      print*, "bsl_lt_pic_4d_write_landau_density_on_remapping_grid -- DEBUG 09890908987876988097986 "
      !      print*, "bsl_lt_pic_4d_write_landau_density_on_remapping_grid -- DEBUG 09890908987876988097986 "
      !      print*, "bsl_lt_pic_4d_write_landau_density_on_remapping_grid -- DEBUG 09890908987876988097986 "
      !      print*, "bsl_lt_pic_4d_write_landau_density_on_remapping_grid -- DEBUG 09890908987876988097986 "
      !      print*, "bsl_lt_pic_4d_write_landau_density_on_remapping_grid -- DEBUG 09890908987876988097986 "
      !      print*, "bsl_lt_pic_4d_write_landau_density_on_remapping_grid -- DEBUG 09890908987876988097986 "
      !      print*, "bsl_lt_pic_4d_write_landau_density_on_remapping_grid -- DEBUG 09890908987876988097986 "
      !      print*, "bsl_lt_pic_4d_write_landau_density_on_remapping_grid -- DEBUG 09890908987876988097986 "
      !      print*, "bsl_lt_pic_4d_write_landau_density_on_remapping_grid -- DEBUG 09890908987876988097986 "
      !
      !      number_nodes_x  = 50
      !      number_nodes_y  = 1
      !      number_nodes_vx = 51
      !      number_nodes_vy = 51
      !
      !      h_x    = 4*sll_pi / number_nodes_x
      !      h_y    = 1./number_nodes_y
      !      h_vx   = 10./(number_nodes_vx-1)
      !      h_vy   = 10./(number_nodes_vy-1)
      !
      !
      !      x_min    = 0
      !      y_min    = 0
      !      vx_min   = -5
      !      vy_min   = -5
      !
      !      total_density = 0
      !
      !      print *, "alpha, k_landau, one_over_thermal_velocity, one_over_two_pi = ", &
      !                    alpha, k_landau, one_over_thermal_velocity, one_over_two_pi
      !
      !      do j_x = 1, number_nodes_x
      !        x_j = x_min + (j_x-1) * h_x
      !        f_x = 1._f64 + alpha * cos(k_landau * x_j) ! eval_landau_fx(alpha, k_landau, x_j)
      !
      !        print *, " DIFF -- ", 1._f64 + alpha * cos(k_landau * x_j),  eval_landau_fx(alpha, k_landau, x_j),    &
      !                  1._f64 + alpha * cos(k_landau * x_j) - eval_landau_fx(alpha, k_landau, x_j)
      !        do j_y = 1, number_nodes_y
      !          ! (density does not depend on y)
      !          do j_vx = 1, number_nodes_vx
      !            vx_j = vx_min + (j_vx-1) * h_vx
      !            f_vx = one_over_thermal_velocity * exp(-0.5*(vx_j*one_over_thermal_velocity)**2)
      !
      !            do j_vy = 1, number_nodes_vy
      !              vy_j = vy_min + (j_vy-1) * h_vy
      !              f_vy = one_over_thermal_velocity * exp(-0.5*(vy_j*one_over_thermal_velocity)**2)
      !
      !              total_density = total_density + ( one_over_two_pi * f_x * f_vx * f_vy ) * h_x * h_y * h_vx * h_vy
      !
      !            end do
      !          end do
      !        end do
      !      end do
      !
      !      print*, "bsl_lt_pic_4d_write_landau_density_on_remapping_grid -- CHECKING TOTAL DENSITY:  ", total_density

    end if

  end subroutine bsl_lt_pic_4d_write_landau_density_on_remapping_grid


  !  subroutine bsl_lt_pic_4d_write_landau_density_on_sparse_grid(    &
  !              p_group,                              &
  !              thermal_speed, alpha, k_landau        &
  !              )
  !
  !    class(sll_bsl_lt_pic_4d_group), intent(inout)   :: p_group
  !    sll_real64, intent(in)                          :: thermal_speed, alpha, k_landau
  !
  !    sll_int32 :: j
  !    sll_real64 :: one_over_thermal_velocity
  !    sll_real64 :: one_over_two_pi
  !    sll_real64 :: x_j
  !    sll_real64 :: vx_j
  !    sll_real64 :: vy_j
  !    sll_real64 :: f_x
  !    sll_real64 :: f_vx
  !    sll_real64 :: f_vy
  !
  !    one_over_thermal_velocity = 1./thermal_speed
  !    one_over_two_pi = 1./(2*sll_pi)
  !
  !    do j=1, p_group%sparse_grid_interpolator%size_basis
  !        x_j = p_group%sparse_grid_interpolator%hierarchy(j)%coordinate(1)
  !        vx_j = p_group%sparse_grid_interpolator%hierarchy(j)%coordinate(3)
  !        vy_j = p_group%sparse_grid_interpolator%hierarchy(j)%coordinate(4)
  !
  !        f_x = eval_landau_fx(alpha, k_landau, x_j)
  !        f_vx = one_over_thermal_velocity * exp(-0.5*(vx_j * one_over_thermal_velocity)**2)
  !        f_vy = one_over_thermal_velocity * exp(-0.5*(vy_j * one_over_thermal_velocity)**2)
  !
  !        p_group%remapped_f_sparse_grid_coefficients(j) = one_over_two_pi * f_x * f_vx * f_vy
  !    end do
  !
  !  end subroutine bsl_lt_pic_4d_write_landau_density_on_sparse_grid


  ! <<bsl_lt_pic_4d_write_hat_density_on_remapping_grid>>
  subroutine bsl_lt_pic_4d_write_hat_density_on_remapping_grid ( &
        p_group,                &
        x0, y0, vx0, vy0,       &
        r_x, r_y, r_vx, r_vy,   &
        basis_height, shift &
      )

    class(sll_bsl_lt_pic_4d_group), intent(inout)     :: p_group
    sll_real64, intent(in)                            :: x0, y0, vx0, vy0
    sll_real64, intent(in)                            :: r_x, r_y, r_vx, r_vy
    sll_real64, intent(in)                            :: basis_height, shift

    sll_int32 :: j
    sll_int32 :: j_x
    sll_int32 :: j_y
    sll_int32 :: j_vx
    sll_int32 :: j_vy
    sll_int32 :: number_nodes_x
    sll_int32 :: number_nodes_y
    sll_int32 :: number_nodes_vx
    sll_int32 :: number_nodes_vy
    sll_real64 :: h_x
    sll_real64 :: h_y
    sll_real64 :: h_vx
    sll_real64 :: h_vy
    sll_real64 :: x_min
    sll_real64 :: y_min
    sll_real64 :: vx_min
    sll_real64 :: vy_min
    sll_real64 :: x_j
    sll_real64 :: y_j
    sll_real64 :: vx_j
    sll_real64 :: vy_j

    if( p_group%remapped_f_interpolation_type == SLL_BSL_LT_PIC_REMAP_WITH_SPLINES )then

      number_nodes_x  = p_group%remapping_cart_grid_number_nodes_x()
      number_nodes_y  = p_group%remapping_cart_grid_number_nodes_y()
      number_nodes_vx = p_group%remapping_cart_grid_number_nodes_vx()
      number_nodes_vy = p_group%remapping_cart_grid_number_nodes_vy()

      h_x    = p_group%remapping_cart_grid%delta_eta1
      h_y    = p_group%remapping_cart_grid%delta_eta2
      h_vx   = p_group%remapping_cart_grid%delta_eta3
      h_vy   = p_group%remapping_cart_grid%delta_eta4

      x_min    = p_group%remapping_cart_grid%eta1_min
      y_min    = p_group%remapping_cart_grid%eta2_min
      vx_min   = p_group%remapping_cart_grid%eta3_min
      vy_min   = p_group%remapping_cart_grid%eta4_min

      ! compute the values of f0 on the (cartesian, phase-space) remapping grid

      do j_x = 1, number_nodes_x
        x_j = x_min + (j_x-1) * h_x

        do j_y = 1, number_nodes_y
          y_j = y_min + (j_y-1) * h_y

          do j_vx = 1, number_nodes_vx
            vx_j = vx_min + (j_vx-1) * h_vx

            do j_vy = 1, number_nodes_vy
              vy_j = vy_min + (j_vy-1) * h_vy

              ! here we store a nodal value but later this array will indeed store spline coefficients
              p_group%remapped_f_splines_coefficients(j_x,j_y,j_vx,j_vy) = eval_hat_function(x0,y0,vx0,vy0,r_x,r_y,r_vx,r_vy,   &  !todo: use temp array
                                                                         basis_height, shift,                                   &
                                                                         x_j, y_j, vx_j, vy_j)
            end do
          end do
        end do
      end do

    else

      SLL_ASSERT( p_group%remapped_f_interpolation_type == SLL_BSL_LT_PIC_REMAP_WITH_SPARSE_GRIDS )

      do j=1, p_group%sparse_grid_interpolator%size_basis
        x_j  = p_group%sparse_grid_interpolator%hierarchy(j)%coordinate(1)
        y_j  = p_group%sparse_grid_interpolator%hierarchy(j)%coordinate(2)
        vx_j = p_group%sparse_grid_interpolator%hierarchy(j)%coordinate(3)
        vy_j = p_group%sparse_grid_interpolator%hierarchy(j)%coordinate(4)

        ! here we store a nodal value but later this array will indeed store sparse grid coefficients
        p_group%remapped_f_sparse_grid_coefficients(j) = eval_hat_function(x0,y0,vx0,vy0,r_x,r_y,r_vx,r_vy,     &
                                                             basis_height, shift,                               &
                                                             x_j, y_j, vx_j, vy_j)
      end do

    end if

  end subroutine bsl_lt_pic_4d_write_hat_density_on_remapping_grid


  !  subroutine bsl_lt_pic_4d_write_hat_density_on_remap_grid ( &
  !        p_group,                &
  !        x0, y0, vx0, vy0,       &
  !        r_x, r_y, r_vx, r_vy,   &
  !        basis_height, hat_shift &
  !      )
  !
  !    class(sll_bsl_lt_pic_4d_group), intent(inout)       :: p_group
  !    sll_real64, intent(in)                          :: x0, y0, vx0, vy0
  !    sll_real64, intent(in)                          :: r_x, r_y, r_vx, r_vy
  !    sll_real64, intent(in)                          :: basis_height, hat_shift
  !
  !    sll_int32 :: j_x
  !    sll_int32 :: j_y
  !    sll_int32 :: j_vx
  !    sll_int32 :: j_vy
  !    sll_int32 :: number_particles
  !    sll_int32 :: number_parts_x
  !    sll_int32 :: number_parts_y
  !    sll_int32 :: number_parts_vx
  !    sll_int32 :: number_parts_vy
  !    sll_real64 :: h_parts_x
  !    sll_real64 :: h_parts_y
  !    sll_real64 :: h_parts_vx
  !    sll_real64 :: h_parts_vy
  !    sll_real64 :: parts_x_min
  !    sll_real64 :: parts_y_min
  !    sll_real64 :: parts_vx_min
  !    sll_real64 :: parts_vy_min
  !    sll_real64 :: x_j
  !    sll_real64 :: y_j
  !    sll_real64 :: vx_j
  !    sll_real64 :: vy_j
  !    !sll_real64 :: f_x, f_y, f_vx, f_vy
  !
  !    number_particles = p_group%number_particles
  !
  !    number_parts_x  = p_group%number_parts_x
  !    number_parts_y  = p_group%number_parts_y
  !    number_parts_vx = p_group%number_parts_vx
  !    number_parts_vy = p_group%number_parts_vy
  !
  !    h_parts_x    = p_group%remapping_cart_grid%delta_eta1
  !    h_parts_y    = p_group%remapping_cart_grid%delta_eta2
  !    h_parts_vx   = p_group%remapping_cart_grid%delta_eta3
  !    h_parts_vy   = p_group%remapping_cart_grid%delta_eta4
  !
  !    parts_x_min    = p_group%remapping_cart_grid%eta1_min
  !    parts_y_min    = p_group%remapping_cart_grid%eta2_min
  !    parts_vx_min   = p_group%remapping_cart_grid%eta3_min
  !    parts_vy_min   = p_group%remapping_cart_grid%eta4_min
  !
  !    ! compute the values of f0 on the (cartesian, phase-space) remapping grid
  !    x_j = parts_x_min
  !    do j_x = 1, number_parts_x
  !      y_j = parts_y_min
  !      do j_y = 1, number_parts_y
  !        vx_j = parts_vx_min
  !        do j_vx = 1, number_parts_vx
  !          vy_j = parts_vy_min
  !          do j_vy = 1, number_parts_vy
  !            p_group%remapped_f_cart_grid_values(j_x,j_y,j_vx,j_vy) = eval_hat_function(x0,y0,vx0,vy0,r_x,r_y,r_vx,r_vy, &
  !                                                                         basis_height, hat_shift,         &
  !                                                                         x_j, y_j, vx_j, vy_j)
  !            vy_j = vy_j + h_parts_vy
  !          end do
  !          vx_j = vx_j + h_parts_vx
  !        end do
  !        y_j = y_j + h_parts_y
  !      end do
  !      x_j = x_j + h_parts_x
  !    end do
  !
  !  end subroutine bsl_lt_pic_4d_write_hat_density_on_remap_grid



  !> compute the spline coefficients for the function values stored in the cartesian remapping grid
  !> nodal values are either given in the array nodal_values_on_remapping_cart_grid (if present),
  !> or stored in the array remapped_f_splines_coefficients
  subroutine bsl_lt_pic_4d_compute_new_spline_coefs( &
              p_group, nodal_values_on_remapping_cart_grid )

    class(sll_bsl_lt_pic_4d_group),       intent(inout) :: p_group
    sll_real64, dimension(:,:,:,:),   target,    intent( in ),  optional  :: nodal_values_on_remapping_cart_grid
    sll_real64, dimension(:,:,:,:),   allocatable, target   :: tmp_nodal_values
    sll_real64, dimension(:,:,:,:),   pointer       :: tmp_nodal_values_ptr
    sll_int32 :: ierr

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
    sll_int32 :: number_nodes_x
    sll_int32 :: number_nodes_y
    sll_int32 :: number_nodes_vx
    sll_int32 :: number_nodes_vy
    sll_real64 :: h_x
    sll_real64 :: h_y
    sll_real64 :: h_vx
    sll_real64 :: h_vy
    sll_real64 :: x_min
    sll_real64 :: y_min
    sll_real64 :: vx_min
    sll_real64 :: vy_min
    sll_real64 :: x_j
    sll_real64 :: y_j
    sll_real64 :: vx_j
    sll_real64 :: vy_j

    sll_real64 :: d_vol
    sll_real64 :: coef

    number_nodes_x  = p_group%remapping_cart_grid_number_nodes_x()
    number_nodes_y  = p_group%remapping_cart_grid_number_nodes_y()
    number_nodes_vx = p_group%remapping_cart_grid_number_nodes_vx()
    number_nodes_vy = p_group%remapping_cart_grid_number_nodes_vy()

    h_x    = p_group%remapping_cart_grid%delta_eta1
    h_y    = p_group%remapping_cart_grid%delta_eta2
    h_vx   = p_group%remapping_cart_grid%delta_eta3
    h_vy   = p_group%remapping_cart_grid%delta_eta4

    x_min    = p_group%remapping_cart_grid%eta1_min
    y_min    = p_group%remapping_cart_grid%eta2_min
    vx_min   = p_group%remapping_cart_grid%eta3_min
    vy_min   = p_group%remapping_cart_grid%eta4_min

    d_vol = h_x * h_y * h_vx * h_vy

    SLL_ASSERT( p_group%remapped_f_interpolation_type == SLL_BSL_LT_PIC_REMAP_WITH_SPLINES )

    if( present(nodal_values_on_remapping_cart_grid) )then
      SLL_ASSERT( size(nodal_values_on_remapping_cart_grid,1) == number_nodes_x  )
      SLL_ASSERT( size(nodal_values_on_remapping_cart_grid,2) == number_nodes_y  )
      SLL_ASSERT( size(nodal_values_on_remapping_cart_grid,3) == number_nodes_vx )
      SLL_ASSERT( size(nodal_values_on_remapping_cart_grid,4) == number_nodes_vy )
      tmp_nodal_values_ptr => nodal_values_on_remapping_cart_grid
    else
      SLL_ALLOCATE( tmp_nodal_values(number_nodes_x, number_nodes_y, number_nodes_vx, number_nodes_vy), ierr )
      tmp_nodal_values = p_group%remapped_f_splines_coefficients
      tmp_nodal_values_ptr => tmp_nodal_values
    end if

    ! compute the spline coefs from the values of the remapped f on the (cartesian, phase-space) remapping grid

    if( p_group%remapped_f_interpolation_degree == 1 )then

      p_group%remapped_f_splines_coefficients = d_vol * tmp_nodal_values_ptr

    else

      ! for higher order splines first store the nodal values

      do j_x = 1, number_nodes_x
        x_j = x_min + (j_x-1) * h_x

        do j_y = 1, number_nodes_y
          y_j = y_min + (j_y-1) * h_y

          do j_vx = 1, number_nodes_vx
            vx_j = vx_min + (j_vx-1) * h_vx

            do j_vy = 1, number_nodes_vy
              vy_j = vy_min + (j_vy-1) * h_vy

              SLL_ASSERT( p_group%remapped_f_interpolation_degree == 3 )

              coef = 0.0_f64
              do l_x = -1, 1
                j_aux_x = j_x + l_x
                if( p_group%domain_is_periodic(1) )then
                  if( j_aux_x < 1 ) j_aux_x = j_aux_x + number_nodes_x
                  if( j_aux_x > number_nodes_x ) j_aux_x = j_aux_x - number_nodes_x
                end if
                if( j_aux_x >= 1 .and. j_aux_x <= number_nodes_x )then
                  do l_y = -1, 1
                    j_aux_y = j_y + l_y
                    if( p_group%domain_is_periodic(2) )then
                      if( j_aux_y < 1 ) j_aux_y = j_aux_y + number_nodes_y
                      if( j_aux_y > number_nodes_y ) j_aux_y = j_aux_y - number_nodes_y
                    end if
                    if( j_aux_y >= 1 .and. j_aux_y <= number_nodes_y )then
                      do l_vx = -1, 1
                        j_aux_vx = j_vx + l_vx
                        if( j_aux_vx >= 1 .and. j_aux_vx <= number_nodes_vx )then
                          do l_vy = -1, 1
                            j_aux_vy = j_vy + l_vy
                            if( j_aux_vy >= 1 .and. j_aux_vy <= number_nodes_vy )then
                              ! MCP: here by discarding outside loop instances we assume
                              ! that non-periodic bc = zero bc. We should instead
                              ! keep those loop instances and use the specified boundary condition when
                              ! the node is in some "fat boundary zone" )
                              coef = coef +                                               &
                                        real(                                             &
                                            p_group%lt_pic_interpolation_coefs(l_x )  *   &
                                            p_group%lt_pic_interpolation_coefs(l_y )  *   &
                                            p_group%lt_pic_interpolation_coefs(l_vx)  *   &
                                            p_group%lt_pic_interpolation_coefs(l_vy)  *   &
                                            tmp_nodal_values_ptr(j_aux_x,j_aux_y,j_aux_vx,j_aux_vy) &
                                            ,f64)
                            end if
                          end do
                        end if
                      end do
                    end if
                  end do
                end if
              end do
              coef = d_vol * coef
              p_group%remapped_f_splines_coefficients(j_x, j_y, j_vx, j_vy) = coef
            end do
          end do
        end do
      end do
    end if

  end subroutine bsl_lt_pic_4d_compute_new_spline_coefs


  !> set the connectivity arrays so that every marker knows its neighbors on the initial (cartesian) grid
  subroutine bsl_lt_pic_4d_set_markers_connectivity( p_group )
    class(sll_bsl_lt_pic_4d_group),intent(inout) :: p_group

    sll_int32  :: k
    sll_int32  :: k_ngb
    sll_int32  :: k_check

    sll_int32  :: j_x
    sll_int32  :: j_y
    sll_int32  :: j_vx
    sll_int32  :: j_vy
    sll_int32  :: number_flow_markers_x
    sll_int32  :: number_flow_markers_y
    sll_int32  :: number_flow_markers_vx
    sll_int32  :: number_flow_markers_vy
    sll_real64 :: h_x
    sll_real64 :: h_y
    sll_real64 :: h_vx
    sll_real64 :: h_vy
    sll_real64 :: x_min
    sll_real64 :: y_min
    sll_real64 :: vx_min
    sll_real64 :: vy_min
    sll_real64 :: x_j
    sll_real64 :: y_j
    sll_real64 :: vx_j
    sll_real64 :: vy_j
    sll_int32,  dimension(:,:,:,:), allocatable    :: markers_indices
    sll_int32 :: ierr

    number_flow_markers_x  = p_group%number_flow_markers_x
    number_flow_markers_y  = p_group%number_flow_markers_y
    number_flow_markers_vx = p_group%number_flow_markers_vx
    number_flow_markers_vy = p_group%number_flow_markers_vy

    h_x    = p_group%initial_markers_grid%delta_eta1
    h_y    = p_group%initial_markers_grid%delta_eta2
    h_vx   = p_group%initial_markers_grid%delta_eta3
    h_vy   = p_group%initial_markers_grid%delta_eta4

    x_min    = p_group%initial_markers_grid%eta1_min
    y_min    = p_group%initial_markers_grid%eta2_min
    vx_min   = p_group%initial_markers_grid%eta3_min
    vy_min   = p_group%initial_markers_grid%eta4_min

    SLL_ALLOCATE(markers_indices(number_flow_markers_x, number_flow_markers_y, number_flow_markers_vx, number_flow_markers_vy),ierr)
    markers_indices(:,:,:,:) = 0

    k_check = 0
    do j_x = 1, number_flow_markers_x
      x_j = x_min + (j_x-1) * h_x

      do j_y = 1, number_flow_markers_y
        y_j = y_min + (j_y-1) * h_y

        do j_vx = 1, number_flow_markers_vx
          vx_j = vx_min + (j_vx-1) * h_vx

          do j_vy = 1, number_flow_markers_vy
            vy_j = vy_min + (j_vy-1) * h_vy

            k_check = k_check + 1
            k = marker_index_from_initial_position_on_cartesian_grid(                         &
                    j_x, j_y, j_vx, j_vy,                                                     &
                    number_flow_markers_x, number_flow_markers_y, number_flow_markers_vx, number_flow_markers_vy  &
                )
            SLL_ASSERT(k == k_check)

            markers_indices( j_x, j_y, j_vx, j_vy ) = k

            if(j_x == 1)then
                ! [neighbor index = own index] means: no neighbor
                ! in the x-periodic case this will be changed when dealing with the last marker in the x dimension
                p_group%markers_list(k)%ngb_xleft_index = k
            else
                ! set the connectivity (in both directions) with left neighbor
                k_ngb = markers_indices(j_x-1,j_y,j_vx,j_vy)
                p_group%markers_list(k)%ngb_xleft_index = k_ngb
                p_group%markers_list(k_ngb)%ngb_xright_index = k
                if(j_x == number_flow_markers_x)then
                    if( p_group%domain_is_periodic(1) )then
                        ! set the connectivity (in both directions) with right neighbor
                        k_ngb = markers_indices(1,j_y,j_vx,j_vy)
                        p_group%markers_list(k)%ngb_xright_index = k_ngb
                        p_group%markers_list(k_ngb)%ngb_xleft_index = k
                    else
                        ! [neighbor index = own index] means: no neighbor
                        p_group%markers_list(k)%ngb_xright_index = k
                    end if
                end if
            end if
            if(j_y == 1)then
                ! [neighbor index = own index] means: no neighbor
                ! in the y-periodic case this will be changed when dealing with the last marker in the y dimension
                p_group%markers_list(k)%ngb_yleft_index = k
            else
                ! set the connectivity (in both directions) with left neighbor
                k_ngb = markers_indices(j_x,j_y-1,j_vx,j_vy)
                p_group%markers_list(k)%ngb_yleft_index = k_ngb
                p_group%markers_list(k_ngb)%ngb_yright_index = k
                if(j_y == number_flow_markers_y)then
                    if( p_group%domain_is_periodic(2) )then
                        ! set the connectivity (in both directions) with right neighbor
                        k_ngb = markers_indices(j_x,1,j_vx,j_vy)
                        p_group%markers_list(k)%ngb_yright_index = k_ngb
                        p_group%markers_list(k_ngb)%ngb_yleft_index = k
                    else
                        ! [neighbor index = own index] means: no neighbor
                        p_group%markers_list(k)%ngb_yright_index = k
                    end if
                end if
            end if
            if(j_vx == 1)then
                ! [neighbor index = own index] means: no neighbor
                p_group%markers_list(k)%ngb_vxleft_index = k
            else
                ! set the connectivity (in both directions) with left neighbor
                k_ngb = markers_indices(j_x,j_y,j_vx-1,j_vy)
                p_group%markers_list(k)%ngb_vxleft_index = k_ngb
                p_group%markers_list(k_ngb)%ngb_vxright_index = k
                if(j_vx == number_flow_markers_vx)then
                    ! [neighbor index = own index] means: no neighbor
                    p_group%markers_list(k)%ngb_vxright_index = k
                end if
            end if
            if(j_vy == 1)then
                ! [neighbor index = own index] means: no neighbor
                p_group%markers_list(k)%ngb_vyleft_index = k
            else
                ! set the connectivity (in both directions) with left neighbor
                k_ngb = markers_indices(j_x,j_y,j_vx,j_vy-1)
                p_group%markers_list(k)%ngb_vyleft_index = k_ngb
                p_group%markers_list(k_ngb)%ngb_vyright_index = k
                if(j_vy == number_flow_markers_vy)then
                    ! [neighbor index = own index] means: no neighbor
                    p_group%markers_list(k)%ngb_vyright_index = k
                end if
            end if

          end do
        end do
      end do
    end do

  end subroutine bsl_lt_pic_4d_set_markers_connectivity


  !> reset the markers on the initial (markers) grid
  subroutine bsl_lt_pic_4d_reset_markers_position( p_group )
    class(sll_bsl_lt_pic_4d_group),intent(inout) :: p_group

    sll_int32 :: k
    sll_int32 :: k_check
    sll_real64, dimension(3)      :: coords

    sll_int32 :: j_x
    sll_int32 :: j_y
    sll_int32 :: j_vx
    sll_int32 :: j_vy
    sll_int32 :: number_flow_markers_x
    sll_int32 :: number_flow_markers_y
    sll_int32 :: number_flow_markers_vx
    sll_int32 :: number_flow_markers_vy
    sll_real64 :: h_x
    sll_real64 :: h_y
    sll_real64 :: h_vx
    sll_real64 :: h_vy
    sll_real64 :: x_min
    sll_real64 :: y_min
    sll_real64 :: vx_min
    sll_real64 :: vy_min
    sll_real64 :: x_j
    sll_real64 :: y_j
    sll_real64 :: vx_j
    sll_real64 :: vy_j

    number_flow_markers_x  = p_group%number_flow_markers_x
    number_flow_markers_y  = p_group%number_flow_markers_y
    number_flow_markers_vx = p_group%number_flow_markers_vx
    number_flow_markers_vy = p_group%number_flow_markers_vy

    h_x    = p_group%initial_markers_grid%delta_eta1
    h_y    = p_group%initial_markers_grid%delta_eta2
    h_vx   = p_group%initial_markers_grid%delta_eta3
    h_vy   = p_group%initial_markers_grid%delta_eta4

    x_min    = p_group%initial_markers_grid%eta1_min
    y_min    = p_group%initial_markers_grid%eta2_min
    vx_min   = p_group%initial_markers_grid%eta3_min
    vy_min   = p_group%initial_markers_grid%eta4_min

    k_check = 0
    do j_x = 1, number_flow_markers_x
      x_j = x_min + (j_x-1) * h_x

      do j_y = 1, number_flow_markers_y
        y_j = y_min + (j_y-1) * h_y

        do j_vx = 1, number_flow_markers_vx
          vx_j = vx_min + (j_vx-1) * h_vx

          do j_vy = 1, number_flow_markers_vy
            vy_j = vy_min + (j_vy-1) * h_vy

            k_check = k_check + 1
            k = marker_index_from_initial_position_on_cartesian_grid(                         &
                    j_x, j_y, j_vx, j_vy,                                                     &
                    number_flow_markers_x, number_flow_markers_y, number_flow_markers_vx, number_flow_markers_vy  &
                )
            SLL_ASSERT(k == k_check)

            coords(1) = x_j
            coords(2) = y_j
            call p_group%set_x( k, coords )

            coords(1) = vx_j
            coords(2) = vy_j
            call p_group%set_v( k, coords )

          end do
        end do
      end do
    end do

  end subroutine bsl_lt_pic_4d_reset_markers_position


  !> compute new interpolation coefficients so that the remapped_f is an approximation of the current f
  !> and
  !> reset the markers on the initial grid
  subroutine bsl_lt_pic_4d_remap ( self )
    class(sll_bsl_lt_pic_4d_group),intent(inout) :: self

    !> A. write the nodal values of the remapped_f (evaluated with the bsl_lt_pic method) and compute the new interpolation coefs
    print *, "bsl_lt_pic_4d_remap -- step A"
    call self%bsl_lt_pic_4d_remap_f()

    !> B. and reset the (flow) markers on the initial (cartesian) grid -- no need to reset their connectivity
    print *, "bsl_lt_pic_4d_remap -- step B"
    call self%bsl_lt_pic_4d_reset_markers_position()

    !> D. if deposition particles are transported, initialize them -- this must be done after the remapping tool is operational!

    print *, "bsl_lt_pic_4d_remap -- step C"
    if( self%deposition_particles_type == SLL_BSL_LT_PIC_TRANSPORTED_RANDOM )then
      call self%bsl_lt_pic_set_deposition_particles_coordinates()    ! todo: try without resetting the particles coordinates?
      call self%bsl_lt_pic_set_deposition_particles_weights()
    end if

  end subroutine bsl_lt_pic_4d_remap



  ! <<bsl_lt_pic_4d_remap_f>> <<ALH>> reconstructs f (with the bsl_lt_pic approach) on the remapping grid,
  ! and compute the new interpolation coefficients

  subroutine bsl_lt_pic_4d_remap_f( p_group )

    class(sll_bsl_lt_pic_4d_group),intent(inout) :: p_group

    type(sll_charge_accumulator_2d),pointer :: void_charge_accumulator
    type(sll_cartesian_mesh_4d),    pointer :: void_grid_4d
    sll_real64, dimension(:,:),     pointer :: void_array_2d

    sll_real64    :: dummy_total_charge
    logical       :: enforce_total_charge
    sll_int32     :: scenario

    nullify(void_charge_accumulator)
    nullify(void_grid_4d)
    nullify(void_array_2d)

    scenario = SLL_BSL_LT_PIC_REMAP_F
    dummy_total_charge = 0.0_f64
    enforce_total_charge = .false.

    ! (no need to distinguish between sparse grid or spline remapping here, this is done inside the function below)
    call p_group%bsl_lt_pic_4d_write_f_on_grid_or_deposit(void_charge_accumulator,          &
                                                          scenario,                         &
                                                          void_grid_4d,                     &
                                                          void_array_2d,                    &
                                                          dummy_total_charge,               &
                                                          enforce_total_charge              &
                                                          )

  end subroutine bsl_lt_pic_4d_remap_f


  !> separate interpolation routine for the remapped f
  function bsl_lt_pic_4d_interpolate_value_of_remapped_f ( p_group, eta ) result(val)
    class(sll_bsl_lt_pic_4d_group), intent(inout)  :: p_group
    sll_real64, dimension(4),       intent(in)  :: eta           !< Position where to interpolate
    sll_real64                                  :: val

    if( p_group%remapped_f_interpolation_type == SLL_BSL_LT_PIC_REMAP_WITH_SPLINES )then

      SLL_ERROR("bsl_lt_pic_4d_interpolate_value_of_remapped_f", "this part not implemented yet")

      ! write it here the spline interpolation ? or use existing modules in Selalib ?

    else if( p_group%remapped_f_interpolation_type == SLL_BSL_LT_PIC_REMAP_WITH_SPARSE_GRIDS )then

      val = p_group%sparse_grid_interpolator%interpolate_value(p_group%remapped_f_sparse_grid_coefficients, eta)

    else

      SLL_ERROR("bsl_lt_pic_4d_interpolate_value_of_remapped_f", "wrong value for p_group%remapped_f_interpolation_type")

    end if
  end function


#define UPDATE_CLOSEST_MARKER_ARRAYS_USING_NEIGHBOR_CELLS(djx,djy,djvx,djvy)                                            \
    do;                                                                                                                 \
        k_neighbor = closest_marker(j_x+(djx), j_y+(djy), j_vx+(djvx), j_vy+(djvy));                                    \
;                                                                                                                       \
        if(k_neighbor /= 0) then;  do          ;                                                                        \
            coords = p_group%get_x(k_neighbor) ;                                                                        \
            x = coords(1) ;                                                                                             \
            y = coords(2) ;                                                                                             \
            coords = p_group%get_v(k_neighbor) ;                                                                        \
            vx = coords(1) ;                                                                                            \
            vy = coords(2) ;                                                                                            \
            call periodic_correction(p_group,x,y) ;                                                                     \
            x_aux = x - p_group%flow_grid%eta1_min;                                                                     \
            y_aux = y - p_group%flow_grid%eta2_min;                                                                     \
            vx_aux = vx - p_group%flow_grid%eta3_min;                                                                   \
            vy_aux = vy - p_group%flow_grid%eta4_min;                                                                   \
            call update_closest_marker_arrays(k_neighbor,                                                               \
                                                x_aux, y_aux, vx_aux, vy_aux,                                           \
                                                j_x, j_y, j_vx, j_vy,                                                   \
                                                h_flow_grid_x,                                                          \
                                                h_flow_grid_y,                                                          \
                                                h_flow_grid_vx,                                                         \
                                                h_flow_grid_vy,                                                         \
                                                closest_marker,                                                         \
                                                closest_marker_distance) ;                                              \
        exit;                                                                                                           \
        end do;                                                                                                         \
        end if;                                                                                                         \
    exit;                                                                                                               \
    end do



  !> new version of the write_f_on_grid_or_deposit routine, with call to separate interpolation routine for the remapped f

  ! todo: update the description below
  ! <<bsl_lt_pic_4d_write_f_on_grid_or_deposit>> <<ALH>> has two scenarios:
  !  - 1.  the "write f" scenario:
  !        write the density on the (phase-space) remapping grid, using the method described
  !        in the "BSL-remapping" notes (version of december 2, 2014)
  !
  !        -- this function should be a faster alternative to [[bsl_lt_pic_4d_write_f_on_remapping_grid]] --
  !
  !        Note: the (x,y)-projection of the remapping grid may be larger than the "Poisson" 2d mesh associated with the
  !        particle group (in particular if the (x,y) domain is not periodic)
  !
  !  - 2.  the "deposition" scenario:
  !        deposit the charge on the Poisson cells given in the given charge_accumulator_2d object
  !
  !  In every case, this routine computes (approximated) values of the density f(t_n) on 'virtual' nodes (sometimes called
  !  'virtual particles') located on a cartesian grid of the 4d phase-space.
  !  For different reasons, these virtual nodes are gathered in 'virtual' cells, and the given arguments
  !  n_virtual_x, n_virtual_y, n_virtual_vx, n_virtual_vy is the number of virtual nodes per virtual cell, in every dimension
  !
  !  -> in the "write f" scenario, these virtual nodes are in fact the nodes of a given grid -- either the remapping one,
  !     or some given grid -- hence they are not really virtual.
  !     On the contrary the virtual cells are really virtual, in the sense that they are just a way to gather the computations
  !
  !  -> in the "deposition" scenario, the virtual nodes are really virtual: they correspond to temporary particles which
  !     are deposited with a standard PIC procedure. And the virtual cells are "half-virtual" in the sense that their (x,y)
  !     projection coincides with the cells of the Poisson mesh, whereas in the velocity dimensions they are created to gather
  !     the computations, just as in the "remapping" scenario
  !
  !     In particular, taking a larger value for n_virtual has the following effect:
  !     -> in the  "write f" scenario, larger virtual cells will be used to compute the approximated values of f(t_n) on the
  !        nodes of the (remapping or given) grid. This will speed-up the code and is morally ok if the characteristic flow
  !        is smooth
  !     -> in the "deposition" scenario, finer grids of virtual point particles will be (temporarily) created and
  !        deposited. This will slow down the code and is morally required if the density f(t_n) is not locally smooth
  !
  !  - Note: in the "write f" scenario, the grid is known through:
  !        the number of grid points (not cells) in every dimension
  !        and the max and min coordinates of these points, in every dimension
  !
  !  - target_total_charge is an optional argument that is given to make the deposition method conservative
  !    (note: it may be used also in the 'write_f' scenario, but one has to define what conservative means in this case)
  !
  !  Note: This routine is an evolution from sll_lt_pic_4d_write_bsl_f_on_remap_grid (which will be eventually discarded)


  subroutine bsl_lt_pic_4d_write_f_on_grid_or_deposit(p_group,                    &
                                                      charge_accumulator,         &
                                                      scenario,                   &
                                                      given_grid_4d,              &
                                                      given_array_2d,             &
                                                      target_total_charge,        &
                                                      enforce_total_charge       &
                                                      )

    class(sll_bsl_lt_pic_4d_group),           intent(inout) :: p_group          !> particle group (with markers and remapping grid)
    type(sll_charge_accumulator_2d), pointer, intent(inout) :: charge_accumulator
    sll_int32,                                intent(in)    :: scenario
    type(sll_cartesian_mesh_4d),     pointer, intent(in)    :: given_grid_4d
    sll_real64, dimension(:,:),      pointer, intent(inout) :: given_array_2d   ! assumed in x, vx for now
    sll_real64,                               intent(in)    :: target_total_charge
    logical,                                  intent(in)    :: enforce_total_charge

    type(charge_accumulator_cell_2d),  pointer :: charge_accumulator_cell

    sll_real64 :: deposited_charge
    sll_real64 :: charge_correction_factor

    logical :: create_deposition_particles_on_a_grid
    logical :: reconstruct_f_on_g_grid

    ! index of marker closest to the center of each flow cell

    sll_int32,          dimension(:,:,:,:), allocatable     :: closest_marker
    sll_real64,         dimension(:,:,:,:), allocatable     :: closest_marker_distance

    ! temporary arrays to store the reconstructed data, in the remapping case
    sll_real64, dimension(:,:,:,:), allocatable, target     :: tmp_f_values_on_remapping_cart_grid
    sll_real64, dimension(:),       allocatable             :: tmp_f_values_on_remapping_sparse_grid

    ! array of integer linked lists (declared below) useful when remapping on sparse grids
    ! (can also simplify the code when remapping on cartesian grids which do not match with the flow cells? but maybe too costly)
    type(int_list_element_ptr), dimension(:,:,:,:), allocatable     :: nodes_in_flow_cell
    type(int_list_element),     pointer                             :: new_int_list_element, head

    sll_int32  :: k ! marker index
    sll_int32  :: k_neighbor
    sll_int32  :: i_x,i_y,i_vx,i_vy  ! grid node indices
    sll_int32  :: j_x,j_y,j_vx,j_vy  ! flow cell indices
    sll_int32  :: m_x,m_y,m_vx,m_vy

    sll_int32  :: i_min_x
    sll_int32  :: i_min_y
    sll_int32  :: i_min_vx
    sll_int32  :: i_min_vy

    sll_int32  :: i_max_x
    sll_int32  :: i_max_y
    sll_int32  :: i_max_vx
    sll_int32  :: i_max_vy

    sll_int32  :: node_counter
    sll_int32  :: node_index

    sll_int32  :: k_marker_closest_to_first_corner

    sll_real64 :: deposition_dvol

    ! <<g>> cartesian grid pointer to the remapping grid
    type(sll_cartesian_mesh_4d),pointer :: g

    sll_int32  :: g_num_points_x
    sll_int32  :: g_num_points_y
    sll_int32  :: g_num_points_vx
    sll_int32  :: g_num_points_vy

    sll_real64 :: g_grid_x_min
    sll_real64 :: g_grid_y_min
    sll_real64 :: g_grid_vx_min
    sll_real64 :: g_grid_vy_min

    sll_real64 :: h_g_grid_x
    sll_real64 :: h_g_grid_y
    sll_real64 :: h_g_grid_vx
    sll_real64 :: h_g_grid_vy


    ! the flow cells are the cells where the flow will be linearized    ---   note: previous name "virtual cells" was too abstract
    ! cf [[file:~/mcp/maltpic/ltpic-bsl.tex::N*]]
    ! same as \delta{x,y,vx,vy} in [[file:~/mcp/maltpic/ltpic-bsl.tex::h_parts_x]]

    sll_real64 :: h_flow_grid_x
    sll_real64 :: h_flow_grid_y
    sll_real64 :: h_flow_grid_vx
    sll_real64 :: h_flow_grid_vy

    sll_int32 :: flow_grid_num_cells_x
    sll_int32 :: flow_grid_num_cells_y
    sll_int32 :: flow_grid_num_cells_vx
    sll_int32 :: flow_grid_num_cells_vy

    sll_real64 :: flow_grid_x_min
    sll_real64 :: flow_grid_y_min
    sll_real64 :: flow_grid_vx_min
    sll_real64 :: flow_grid_vy_min

    ! the markers are initially distributed on a cartesian grid, then pushed forward to represent (and approximate) the flow
    ! cf [[file:~/mcp/maltpic/ltpic-bsl.tex::h_parts_x]]

    sll_real64 :: h_markers_x
    sll_real64 :: h_markers_y
    sll_real64 :: h_markers_vx
    sll_real64 :: h_markers_vy

    sll_real64 :: markers_x_min
    sll_real64 :: markers_y_min
    sll_real64 :: markers_vx_min
    sll_real64 :: markers_vy_min

    !    sll_real64 :: h_virtual_parts_x
    !    sll_real64 :: h_virtual_parts_y
    !    sll_real64 :: h_virtual_parts_vx
    !    sll_real64 :: h_virtual_parts_vy
    !
    !    sll_real64 :: inv_h_virtual_parts_x
    !    sll_real64 :: inv_h_virtual_parts_y
    !    sll_real64 :: inv_h_virtual_parts_vx
    !    sll_real64 :: inv_h_virtual_parts_vy

!    sll_int32 :: number_virtual_particles_x
!    sll_int32 :: number_virtual_particles_y
!    sll_int32 :: number_virtual_particles_vx
!    sll_int32 :: number_virtual_particles_vy

!    sll_real64 :: phase_space_virtual_dvol
!
!    sll_real64 :: virtual_parts_x_min
!    sll_real64 :: virtual_parts_y_min
!    sll_real64 :: virtual_parts_vx_min
!    sll_real64 :: virtual_parts_vy_min

!    sll_real64 :: virtual_grid_x_min
!    sll_real64 :: virtual_grid_x_max
!    sll_real64 :: virtual_grid_y_min
!    sll_real64 :: virtual_grid_y_max
!    sll_real64 :: virtual_grid_vx_min
!    sll_real64 :: virtual_grid_vx_max
!    sll_real64 :: virtual_grid_vy_min
!    sll_real64 :: virtual_grid_vy_max
!
    sll_int32 :: number_of_deposition_particles_per_flow_cell

    sll_real64 :: x
    sll_real64 :: y
    sll_real64 :: vx
    sll_real64 :: vy

    sll_real64 :: closest_marker_distance_to_first_corner
    sll_real64 :: marker_distance_to_first_corner

    ! working space
    sll_real64 :: tmp, tmp1, tmp2

    sll_real64 :: mesh_period_x
    sll_real64 :: mesh_period_y
    ! sll_real64 :: inv_period_x
    ! sll_real64 :: inv_period_y

    ! results from [[get_ltp_deformation_matrix]]

    sll_real64 :: d11,d12,d13,d14 ! coefs of matrix D (backward Jacobian)
    sll_real64 :: d21,d22,d23,d24
    sll_real64 :: d31,d32,d33,d34
    sll_real64 :: d41,d42,d43,d44

    sll_real64, dimension(3)  :: coords
    sll_real64, dimension(4)  :: eta        !> coordinates in the logical (cartesian) 4d space

    ! coordinates of particle k at time n and time 0
    sll_real64 :: x_k,y_k,vx_k,vy_k
    sll_real64 :: x_k_t0,y_k_t0,vx_k_t0,vy_k_t0

    sll_real64 :: x_to_xk, y_to_yk, vx_to_vxk, vy_to_vyk

    sll_real64 :: d1_x, d1_y, d1_vx, d1_vy
    sll_real64 :: d2_x, d2_y, d2_vx, d2_vy
    sll_real64 :: d3_x, d3_y, d3_vx, d3_vy
    sll_real64 :: d4_x, d4_y, d4_vx, d4_vy

    sll_real64 :: cell_offset_x
    sll_real64 :: cell_offset_y

    sll_real64 :: deposition_particle_charge_factor
    sll_real64 :: phase_space_volume

    sll_int32  :: nodes_number
    sll_real64 :: reconstructed_f_value
    sll_real64 :: reconstructed_charge

    ! coordinates of a reconstruction point at time 0, absolute
    sll_real64 :: x_t0
    sll_real64 :: y_t0
    sll_real64 :: vx_t0
    sll_real64 :: vy_t0

    ! coordinates of a reconstruction point at time 0, relative to the nearby reference marker (= closest to cell center)
    sll_real64 :: x_t0_to_xk_t0
    sll_real64 :: y_t0_to_yk_t0
    sll_real64 :: vx_t0_to_vxk_t0
    sll_real64 :: vy_t0_to_vyk_t0

    sll_int32  :: ierr
    sll_int32  :: i_cell_x
    sll_int32  :: i_cell_y
    sll_int32  :: i_cell

    ! temporary workspace
    sll_real64 :: x_aux
    sll_real64 :: y_aux
    sll_real64 :: vx_aux
    sll_real64 :: vy_aux

    sll_real64 :: debug_charge
    sll_int32  :: debug_count


    ! --- end of declarations!

    debug_charge = 0.0_f64
    debug_count = 0

    ! getting the parameters of the flow grid
    flow_grid_x_min    = p_group%flow_grid%eta1_min
    flow_grid_y_min    = p_group%flow_grid%eta2_min
    flow_grid_vx_min   = p_group%flow_grid%eta3_min
    flow_grid_vy_min   = p_group%flow_grid%eta4_min

    h_flow_grid_x  = p_group%flow_grid%delta_eta1
    h_flow_grid_y  = p_group%flow_grid%delta_eta2
    h_flow_grid_vx = p_group%flow_grid%delta_eta3
    h_flow_grid_vy = p_group%flow_grid%delta_eta4

    flow_grid_num_cells_x  = p_group%flow_grid%num_cells1
    flow_grid_num_cells_y  = p_group%flow_grid%num_cells2
    flow_grid_num_cells_vx = p_group%flow_grid%num_cells3
    flow_grid_num_cells_vy = p_group%flow_grid%num_cells4

    ! initialize the multi-purpose g grid
    nullify(g)
    reconstruct_f_on_g_grid = .false.
    g_num_points_x  = 0
    g_num_points_y  = 0
    g_num_points_vx = 0
    g_num_points_vy = 0

    deposited_charge = 0.0_f64

    !> A.  preparation of the point sets where f will be reconstructed, depending on the different scenarios
    if( scenario == SLL_BSL_LT_PIC_DEPOSIT_F )then

      ! choose between cartesian grid or random bunches of deposition (=virtual) particles
      create_deposition_particles_on_a_grid = .true.

      if( create_deposition_particles_on_a_grid )then

        reconstruct_f_on_g_grid = .true.
        g => p_group%deposition_grid

        ! the boundary nodes of the deposition grid are inside the domain, even with periodic boundary conditions
        g_num_points_x  = g%num_cells1 + 1
        g_num_points_y  = g%num_cells2 + 1
        g_num_points_vx = g%num_cells3 + 1
        g_num_points_vy = g%num_cells4 + 1


      else

        ! the deposition particles will be created to deposit their charge but not stored in memory
        !        number_of_deposition_particles_per_flow_cell = number_of_deposition_particles / (  flow_grid_num_cells_x    &
  !                                                                                               * flow_grid_num_cells_y    &
    !                                                                                             * flow_grid_num_cells_vx   &
    !                                                                                             * flow_grid_num_cells_vy )

        SLL_ASSERT( .not. reconstruct_f_on_g_grid )
        SLL_ERROR("bsl_lt_pic_4d_write_f_on_grid_or_deposit", " this part not implemented yet...")

      end if

    else if(  scenario == SLL_BSL_LT_PIC_REMAP_F                                                  &
              .and. p_group%remapped_f_interpolation_type == SLL_BSL_LT_PIC_REMAP_WITH_SPLINES    &
      )then

      ! with splines, the remapping grid is cartesian

      reconstruct_f_on_g_grid = .true.
      g => p_group%remapping_cart_grid

      g_num_points_x  = p_group%remapping_cart_grid_number_nodes_x()
      g_num_points_y  = p_group%remapping_cart_grid_number_nodes_y()
      g_num_points_vx = p_group%remapping_cart_grid_number_nodes_vx()
      g_num_points_vy = p_group%remapping_cart_grid_number_nodes_vy()

      ! allocate temp array to store the nodal values of the remapped f
      SLL_ALLOCATE(tmp_f_values_on_remapping_cart_grid(g_num_points_x, g_num_points_y, g_num_points_vx, g_num_points_vy), ierr)
      tmp_f_values_on_remapping_cart_grid = 0.0_f64

    else if( (scenario == SLL_BSL_LT_PIC_SET_WEIGHTS_ON_DEPOSITION_PARTICLES)   &
        .or.                                                                    &
        (scenario == SLL_BSL_LT_PIC_REMAP_F .and. p_group%remapped_f_interpolation_type == SLL_BSL_LT_PIC_REMAP_WITH_SPARSE_GRIDS) &
      )then

      ! common preparation step for sparse grid remapping and weights computing for unstructured cloud of deposition particles
      ! => prepare the array of linked lists that will store the node indices contained in the flow cells (one list per cell)
      allocate(nodes_in_flow_cell(flow_grid_num_cells_x,   &
                                  flow_grid_num_cells_y,   &
                                  flow_grid_num_cells_vx,  &
                                  flow_grid_num_cells_vy)  &
             , stat=ierr)
      call test_error_code(ierr, 'Memory allocation Failure.', __FILE__, __LINE__)

      do j_x = 1, flow_grid_num_cells_x
        do j_y = 1, flow_grid_num_cells_y
          do j_vx = 1, flow_grid_num_cells_vx
            do j_vy = 1, flow_grid_num_cells_vy
              nullify(nodes_in_flow_cell(j_x,j_y,j_vx,j_vy)%pointed_element)
            end do
          end do
        end do
      end do


      if( scenario == SLL_BSL_LT_PIC_REMAP_F )then
        SLL_ASSERT( p_group%remapped_f_interpolation_type == SLL_BSL_LT_PIC_REMAP_WITH_SPARSE_GRIDS )
        !        nodes_coordinate_list => p_group%sparse_grid_interpolator%hierarchy(node_index)%coordinate
        nodes_number = p_group%sparse_grid_interpolator%size_basis

        ! allocate temp array to store the nodal values of the remapped f
        !        print*, "7646547 --- will allocate with ...%size_basis = ", p_group%sparse_grid_interpolator%size_basis
        SLL_ALLOCATE(tmp_f_values_on_remapping_sparse_grid(p_group%sparse_grid_interpolator%size_basis), ierr)
        tmp_f_values_on_remapping_sparse_grid = 0.0_f64

      else if( scenario == SLL_BSL_LT_PIC_SET_WEIGHTS_ON_DEPOSITION_PARTICLES )then
        SLL_ASSERT( p_group%deposition_particles_type == SLL_BSL_LT_PIC_TRANSPORTED_RANDOM )
        !        nodes_coordinate_list => p_group%sparse_grid_interpolator%hierarchy(node_index)%coordinate
        nodes_number = p_group%number_moving_deposition_particles

        phase_space_volume =    (p_group%remapping_grid_eta_max(4) - p_group%remapping_grid_eta_min(4))    &
                              * (p_group%remapping_grid_eta_max(3) - p_group%remapping_grid_eta_min(3))    &
                              * (p_group%remapping_grid_eta_max(2) - p_group%remapping_grid_eta_min(2))    &
                              * (p_group%remapping_grid_eta_max(1) - p_group%remapping_grid_eta_min(1))

        deposition_particle_charge_factor = phase_space_volume * p_group%species%q / p_group%number_moving_deposition_particles

        ! reset the weights of the deposition particles, because maybe not every deposition particle weight will be set
        p_group%deposition_particles_weight = 0

      end if

      ! then loop to store the sparse grid node indices in linked lists corresponding to the flow cells that contain them
      do node_index = 1, nodes_number

        ! get node coordinates:
        if( scenario == SLL_BSL_LT_PIC_REMAP_F )then
          x = p_group%sparse_grid_interpolator%hierarchy(node_index)%coordinate(1)
          y = p_group%sparse_grid_interpolator%hierarchy(node_index)%coordinate(2)
          vx = p_group%sparse_grid_interpolator%hierarchy(node_index)%coordinate(3)
          vy = p_group%sparse_grid_interpolator%hierarchy(node_index)%coordinate(4)

        else if( scenario == SLL_BSL_LT_PIC_SET_WEIGHTS_ON_DEPOSITION_PARTICLES )then
          x = p_group%deposition_particles_eta(node_index, 1)
          y = p_group%deposition_particles_eta(node_index, 2)
          vx = p_group%deposition_particles_eta(node_index, 3)
          vy = p_group%deposition_particles_eta(node_index, 4)

        end if

        ! find the index (j_x,j_y,j_vx,j_vy) of the flow cell containing this node (same piece of code as below)
        x_aux = x - flow_grid_x_min
        j_x = int( x_aux / h_flow_grid_x ) + 1

        y_aux = y - flow_grid_y_min
        j_y = int( y_aux / h_flow_grid_y ) + 1

        vx_aux = vx - flow_grid_vx_min
        j_vx = int( vx_aux / h_flow_grid_vx ) + 1

        vy_aux = vy - flow_grid_vy_min
        j_vy = int( vy_aux / h_flow_grid_vy ) + 1

        ! discard if flow cell is off-bounds
        if(  j_x >= 1 .and. j_x <= flow_grid_num_cells_x .and. &
             j_y >= 1 .and. j_y <= flow_grid_num_cells_y .and. &
             j_vx >= 1 .and. j_vx <= flow_grid_num_cells_vx .and. &
             j_vy >= 1 .and. j_vy <= flow_grid_num_cells_vy  )then

          ! increment the proper linked list
          SLL_ALLOCATE( new_int_list_element, ierr )
          new_int_list_element%value = node_index
          head => nodes_in_flow_cell(j_x,j_y,j_vx,j_vy)%pointed_element
          nodes_in_flow_cell(j_x,j_y,j_vx,j_vy)%pointed_element => add_element_in_list(head, new_int_list_element)

        end if

      end do

      SLL_ASSERT( .not. reconstruct_f_on_g_grid )

    else if( scenario == SLL_BSL_LT_PIC_WRITE_F_ON_GIVEN_GRID )then

      ! then use the given 4d grid and write values in given (x, vx for now) array given_array_2d
      reconstruct_f_on_g_grid = .true.
      g => given_grid_4d

      g_num_points_x = given_grid_4d%num_cells1 + 1
      g_num_points_y = given_grid_4d%num_cells2 + 1
      g_num_points_vx = given_grid_4d%num_cells3 + 1
      g_num_points_vy = given_grid_4d%num_cells4 + 1

      ! for now we assume that given_array_2d is in (x, vx) space -- and matches the given grid in 4d where f should be evaluated
      SLL_ASSERT( size(given_array_2d,1) == g_num_points_x )
      SLL_ASSERT( size(given_array_2d,2) == g_num_points_vx )
      given_array_2d(:,:) = 0.0_f64

    else

      SLL_ERROR("bsl_lt_pic_4d_write_f_on_grid_or_deposit", "unknown value for parameter scenario")

    end if

    if( reconstruct_f_on_g_grid )then

      SLL_ASSERT( associated(g) )

      SLL_ASSERT( g_num_points_x > 1 )
      SLL_ASSERT( g_num_points_y > 1 )
      SLL_ASSERT( g_num_points_vx > 1 )
      SLL_ASSERT( g_num_points_vy > 1 )

      h_g_grid_x  = g%delta_eta1
      h_g_grid_y  = g%delta_eta2
      h_g_grid_vx = g%delta_eta3
      h_g_grid_vy = g%delta_eta4

      g_grid_x_min  = g%eta1_min
      g_grid_y_min  = g%eta2_min
      g_grid_vx_min = g%eta3_min
      g_grid_vy_min = g%eta4_min

      if( scenario == SLL_BSL_LT_PIC_DEPOSIT_F )then
        deposition_dvol = h_g_grid_x * h_g_grid_y * h_g_grid_vx * h_g_grid_vy
      end if
    end if

    !> B. Preparatory work for the linearization of the flow on the flow cells:
    !>    - find out the closest marker to each cell center,
    !>      by looping over all markers and noting which flow cell contains it.
    !>      (The leftmost flow cell in each dimension may not be complete.)

    SLL_ALLOCATE(closest_marker(flow_grid_num_cells_x,flow_grid_num_cells_y,flow_grid_num_cells_vx,flow_grid_num_cells_vy),ierr)
    closest_marker(:,:,:,:) = 0

    allocate(closest_marker_distance(flow_grid_num_cells_x,   &
                                     flow_grid_num_cells_y,   &
                                     flow_grid_num_cells_vx,  &
                                     flow_grid_num_cells_vy)  &
               , stat=ierr)
    call test_error_code(ierr, 'Memory allocation Failure.', __FILE__, __LINE__)
    closest_marker_distance(:,:,:,:) = 0.0_f64

    closest_marker_distance_to_first_corner = 1d30
    k_marker_closest_to_first_corner = 0

    do k=1, p_group%number_flow_markers    ! [[file:../pic_particle_types/lt_pic_4d_group.F90::number_particles]]

      ! find absolute (x_k,y_k,vx_k,vy_k) coordinates for k-th marker.
      coords = p_group%get_x(k)
      x_k = coords(1)
      y_k = coords(2)
      coords = p_group%get_v(k)
      vx_k = coords(1)
      vy_k = coords(2)

      ! which _virtual_ cell is this particle in?

      x_aux = x_k - flow_grid_x_min
      j_x = int( x_aux / h_flow_grid_x ) + 1

      y_aux = y_k - flow_grid_y_min
      j_y = int( y_aux / h_flow_grid_y ) + 1

      vx_aux = vx_k - flow_grid_vx_min
      j_vx = int( vx_aux / h_flow_grid_vx ) + 1

      vy_aux = vy_k - flow_grid_vy_min
      j_vy = int( vy_aux / h_flow_grid_vy ) + 1

      ! discard this marker if the flow cell is off-bounds
      if(  j_x >= 1 .and. j_x <= flow_grid_num_cells_x .and. &
           j_y >= 1 .and. j_y <= flow_grid_num_cells_y .and. &
           j_vx >= 1 .and. j_vx <= flow_grid_num_cells_vx .and. &
           j_vy >= 1 .and. j_vy <= flow_grid_num_cells_vy  )then

        call update_closest_marker_arrays(k,                              &
                                          x_aux, y_aux, vx_aux, vy_aux,   &
                                          j_x, j_y, j_vx, j_vy,           &
                                          h_flow_grid_x,                  &
                                          h_flow_grid_y,                  &
                                          h_flow_grid_vx,                 &
                                          h_flow_grid_vy,                 &
                                          closest_marker,                 &
                                          closest_marker_distance)
      end if

      marker_distance_to_first_corner = abs(x_aux) + abs(y_aux) + abs(vx_aux) + abs(vy_aux)
      if( marker_distance_to_first_corner < closest_marker_distance_to_first_corner )then
        closest_marker_distance_to_first_corner = marker_distance_to_first_corner
        k_marker_closest_to_first_corner = k
      end if
    end do

    closest_marker(1,1,1,1) = k_marker_closest_to_first_corner

    if( .not. ( p_group%domain_is_periodic(1) .and. p_group%domain_is_periodic(2) ) )then
      print*, "WARNING -- STOP -- verify that the non-periodic case is well implemented"
      stop
    end if

    if(p_group%domain_is_periodic(1)) then
      ! here the domain corresponds to the Poisson mesh
      mesh_period_x = p_group%space_mesh_2d%eta1_max - p_group%space_mesh_2d%eta1_min
    else
      mesh_period_x = 0.0_f64
    end if

    if(p_group%domain_is_periodic(2)) then
      ! here the domain corresponds to the Poisson mesh
      mesh_period_y = p_group%space_mesh_2d%eta2_max - p_group%space_mesh_2d%eta2_min
    else
      mesh_period_y = 0.0_f64
    end if


    !> C. loop on the flow cells (main loop) -- on each flow cell, we
    !>   - C.1 linearize the flow using the position of the markers
    !>   - C.2 find the relevant points where f should be reconstructed
    !>   - C.3 reconstruct f on these points (using the affine backward flow and the interpolation tool for the remapped_f)
    !>   - C.4 write the resulting f value or deposit the deposition particle just created (depending on the scenario)

    ! <<loop_on_flow_cells>> [[file:~/mcp/maltpic/ltpic-bsl.tex::algo:pic-vr:loop_over_all_cells]]

    ! cell size of the initial_markers_grid, for finite differencing of the flow  - same as in [[write_f_on_remap_grid-h_parts_x]]
    h_markers_x    = p_group%initial_markers_grid%delta_eta1
    h_markers_y    = p_group%initial_markers_grid%delta_eta2
    h_markers_vx   = p_group%initial_markers_grid%delta_eta3
    h_markers_vy   = p_group%initial_markers_grid%delta_eta4

    markers_x_min    = p_group%initial_markers_grid%eta1_min
    markers_y_min    = p_group%initial_markers_grid%eta2_min
    markers_vx_min   = p_group%initial_markers_grid%eta3_min
    markers_vy_min   = p_group%initial_markers_grid%eta4_min

    if( (scenario == SLL_BSL_LT_PIC_REMAP_F)                              &
                .and.                                                     &
                (p_group%remapped_f_interpolation_type == SLL_BSL_LT_PIC_REMAP_WITH_SPARSE_GRIDS) )then
      node_counter = 0
    end if

    do j_x = 1, flow_grid_num_cells_x
      do j_y = 1, flow_grid_num_cells_y
        do j_vx = 1,flow_grid_num_cells_vx
          do j_vy = 1,flow_grid_num_cells_vy

            ! [[file:~/mcp/maltpic/ltpic-bsl.tex::algo:pic-vr:find_closest_real_particle]] Find the marker
            ! which is closest to the flow cell center.  Note: speed-wise, it may be necessary to find a way not to scan
            ! all the markers for every cell.  We avoid scanning all the markers for each cell by using the
            ! precomputed array [[closest_marker]]. Flow cells which do not contain any marker are skipped.

            k = closest_marker(j_x,j_y,j_vx,j_vy)

            if(k == 0) then
              if( j_x > 1 )then
                UPDATE_CLOSEST_MARKER_ARRAYS_USING_NEIGHBOR_CELLS(-1,0,0,0)
              end if
              if( j_x < flow_grid_num_cells_x )then
                UPDATE_CLOSEST_MARKER_ARRAYS_USING_NEIGHBOR_CELLS( 1,0,0,0)
              end if

              if( j_y > 1 )then
                UPDATE_CLOSEST_MARKER_ARRAYS_USING_NEIGHBOR_CELLS(0,-1,0,0)
              end if
              if( j_y < flow_grid_num_cells_y )then
                UPDATE_CLOSEST_MARKER_ARRAYS_USING_NEIGHBOR_CELLS(0, 1,0,0)
              end if

              if( j_vx > 1 )then
                UPDATE_CLOSEST_MARKER_ARRAYS_USING_NEIGHBOR_CELLS(0,0,-1,0)
              end if
              if( j_vx < flow_grid_num_cells_vx )then
                UPDATE_CLOSEST_MARKER_ARRAYS_USING_NEIGHBOR_CELLS(0,0, 1,0)
              end if

              if( j_vy > 1 )then
                UPDATE_CLOSEST_MARKER_ARRAYS_USING_NEIGHBOR_CELLS(0,0,0,-1)
              end if
              if( j_vy < flow_grid_num_cells_vy )then
                UPDATE_CLOSEST_MARKER_ARRAYS_USING_NEIGHBOR_CELLS(0,0,0, 1)
              end if
            end if

            k = closest_marker(j_x,j_y,j_vx,j_vy)
            SLL_ASSERT(k /= 0)

            !>   - C.1 linearize the flow using the position of the markers

            ! [[file:~/mcp/maltpic/ltpic-bsl.tex::hat-bz*]] In this flow cell we will use the
            ! k-th backward flow, obtained with the deformation matrix for the k-th marker see [[get_ltp_deformation_matrix]].
            ! Call below uses parameters inspired from [[sll_lt_pic_4d_write_f_on_remap_grid-get_ltp_deformation_matrix]]

            call p_group%get_ltp_deformation_matrix (       &
                 k,                                         &
                 mesh_period_x,                             &
                 mesh_period_y,                             &
                 h_markers_x,                               &
                 h_markers_y,                               &
                 h_markers_vx,                              &
                 h_markers_vy,                              &
                 x_k,y_k,vx_k,vy_k,                         &
                 d11,d12,d13,d14,                           &
                 d21,d22,d23,d24,                           &
                 d31,d32,d33,d34,                           &
                 d41,d42,d43,d44                            &
                 )

            ! Find position of marker k at time 0
            ! [[get_initial_position_on_cartesian_grid_from_marker_index]]

            call get_initial_position_on_cartesian_grid_from_marker_index(k,                    &
                 p_group%number_flow_markers_x, p_group%number_flow_markers_y,                  &
                 p_group%number_flow_markers_vx, p_group%number_flow_markers_vy,                &
                 m_x,m_y,m_vx,m_vy)

            x_k_t0  = markers_x_min  + (m_x-1)  * h_markers_x
            y_k_t0  = markers_y_min  + (m_y-1)  * h_markers_y
            vx_k_t0 = markers_vx_min + (m_vx-1) * h_markers_vx
            vy_k_t0 = markers_vy_min + (m_vy-1) * h_markers_vy


            ! comment below (and pointers) should be udated
            ! <<loop_on_virtual_particles_in_one_flow_cell>>
            ! [[file:~/mcp/maltpic/ltpic-bsl.tex::algo:pic-vr:find_f0_for_each_virtual_particle]] Loop over all
            ! flow particles in the cell to compute the value of f0 at that point (Following
            ! [[file:~/mcp/maltpic/ltpic-bsl.tex::BSL_remapping_algo]])
            ! [[file:~/mcp/maltpic/ltpic-bsl.tex::algo:pic-vr:create_virtual_particles]] Create a temporary set of
            ! virtual particles inside the cell.

            !>  - C.2 find the relevant points (x, y, vx, vy) where f should be reconstructed

            !>    - C.2.a first we treat the case of [remapping with a sparse grid]
            !>            or [computing the weights of a cloud of deposition particles]: nodes are stored in linked lists

            if( ( (scenario == SLL_BSL_LT_PIC_REMAP_F)                                                    &
                  .and.                                                                                   &
                  (p_group%remapped_f_interpolation_type == SLL_BSL_LT_PIC_REMAP_WITH_SPARSE_GRIDS) )     &
                .or.                                                                                      &
                  (scenario == SLL_BSL_LT_PIC_SET_WEIGHTS_ON_DEPOSITION_PARTICLES)                        &
              )then

              SLL_ASSERT( .not. reconstruct_f_on_g_grid )

              new_int_list_element => nodes_in_flow_cell(j_x,j_y,j_vx,j_vy)%pointed_element

              do while( associated(new_int_list_element) )
                node_index = new_int_list_element%value

                ! here we reconstruct f on the sparse grid nodes
                if( scenario == SLL_BSL_LT_PIC_REMAP_F )then
                  x  = p_group%sparse_grid_interpolator%hierarchy(node_index)%coordinate(1)
                  y  = p_group%sparse_grid_interpolator%hierarchy(node_index)%coordinate(2)
                  vx = p_group%sparse_grid_interpolator%hierarchy(node_index)%coordinate(3)
                  vy = p_group%sparse_grid_interpolator%hierarchy(node_index)%coordinate(4)
                else
                  SLL_ASSERT( scenario == SLL_BSL_LT_PIC_SET_WEIGHTS_ON_DEPOSITION_PARTICLES )
                  x  = p_group%deposition_particles_eta(node_index, 1)
                  y  = p_group%deposition_particles_eta(node_index, 2)
                  vx = p_group%deposition_particles_eta(node_index, 3)
                  vy = p_group%deposition_particles_eta(node_index, 4)
                end if

                x_to_xk   = x - x_k
                y_to_yk   = y - y_k
                vx_to_vxk = vx - vx_k
                vy_to_vyk = vy - vy_k

                ! find z_t0 = (x_t0, y_t0, vx_t0, vy_t0), the position of the node at time = 0 using the affine bwd flow
                x_t0_to_xk_t0   = d11 * x_to_xk + d12 * y_to_yk + d13 * vx_to_vxk + d14 * vy_to_vyk
                y_t0_to_yk_t0   = d21 * x_to_xk + d22 * y_to_yk + d23 * vx_to_vxk + d24 * vy_to_vyk
                vx_t0_to_vxk_t0 = d31 * x_to_xk + d32 * y_to_yk + d33 * vx_to_vxk + d34 * vy_to_vyk
                vy_t0_to_vyk_t0 = d41 * x_to_xk + d42 * y_to_yk + d43 * vx_to_vxk + d44 * vy_to_vyk

                x_t0  = x_t0_to_xk_t0   + x_k_t0
                y_t0  = y_t0_to_yk_t0   + y_k_t0
                vx_t0 = vx_t0_to_vxk_t0 + vx_k_t0
                vy_t0 = vy_t0_to_vyk_t0 + vy_k_t0

                ! put back z_t0 inside domain (if needed) to enforce periodic boundary conditions
                call periodic_correction(p_group,x_t0,y_t0)

                ! now (x_t0, y_t0, vx_t0, vy_t0) is the (approx) position of the node at time t=0
                eta(1) = x_t0
                eta(2) = y_t0
                eta(3) = vx_t0
                eta(4) = vy_t0

                ! interpolation here may use sparse grid or splines, depending on the method chosen for p_group
                reconstructed_f_value = p_group%bsl_lt_pic_4d_interpolate_value_of_remapped_f(eta)

                ! here we store a nodal value but later this array will indeed store sparse grid coefficients
                if( scenario == SLL_BSL_LT_PIC_REMAP_F )then
                  tmp_f_values_on_remapping_sparse_grid(node_index) = reconstructed_f_value
                else
                  SLL_ASSERT( scenario == SLL_BSL_LT_PIC_SET_WEIGHTS_ON_DEPOSITION_PARTICLES )
                  SLL_ASSERT( p_group%deposition_particles_type == SLL_BSL_LT_PIC_TRANSPORTED_RANDOM )
                  p_group%deposition_particles_weight(node_index) = reconstructed_f_value * deposition_particle_charge_factor
                end if

                ! [DEBUG]
                if( .false. )then
                  print*, "[REMAP] reconstructing for node_index = ", node_index
                  print*, "from:                  x, y, vx, vy = ", x, y, vx, vy
                  print*, "using             x_k,y_k,vx_k,vy_k = ", x_k,y_k,vx_k,vy_k
                  print*, "x_to_xk , y_to_yk , vx_to_vxk , vy_to_vyk = ", x_to_xk , y_to_yk , vx_to_vxk , vy_to_vyk
                  print*, " d11,d12,d13,d14  = ", d11,d12,d13,d14
                  print*, " d21,d22,d23,d24  = ", d21,d22,d23,d24
                  print*, " d31,d32,d33,d34  = ", d31,d32,d33,d34
                  print*, " d41,d42,d43,d44  = ", d41,d42,d43,d44

                  print*, "and: eta_t0_to_etak_t0 = ",  x_t0_to_xk_t0, y_t0_to_yk_t0, vx_t0_to_vxk_t0, vy_t0_to_vyk_t0
                  print*, "on:                    x_t0, y_t0, vx_t0, vy_t0 = ", x_t0, y_t0, vx_t0, vy_t0
                  print*, "on:                                         eta = ", eta

                  print*, "value:  reconstructed_f_value = ", reconstructed_f_value

                  stop
                end if

                node_counter = node_counter + 1

                new_int_list_element => new_int_list_element%next

              end do

            !>    - C.2.b next we treat the case of a deposition scenario

            else if( (scenario == SLL_BSL_LT_PIC_DEPOSIT_F) .and. (.not. create_deposition_particles_on_a_grid) )then

              SLL_ASSERT( .not. reconstruct_f_on_g_grid )

              SLL_ERROR("bsl_lt_pic_4d_write_f_on_grid_or_deposit", "this part not implemented yet")

              ! here we should
              !         - create a bunch of random deposition particles within this flow cell,
              !         - then for each deposition particle:
              !         - reconstruct the value of f there to get their charge
              !         - and deposit these charges on the accumulator cells

            !>    - C.2.c finally we treat the case of a remapping with splines or writing on a given grid
            !>      (in both cases the nodes are on the g grid constructed above)

            else if( (scenario == SLL_BSL_LT_PIC_WRITE_F_ON_GIVEN_GRID)                                       &
                     .or.                                                                                     &
                     ( (scenario == SLL_BSL_LT_PIC_REMAP_F)                                                   &
                       .and. (p_group%remapped_f_interpolation_type == SLL_BSL_LT_PIC_REMAP_WITH_SPLINES ) )  &
                     .or.                                                                                     &
                     ( (scenario == SLL_BSL_LT_PIC_DEPOSIT_F)                                 &
                       .and. create_deposition_particles_on_a_grid ) &
                       )then

              SLL_ASSERT( reconstruct_f_on_g_grid )

              !              if( (scenario == SLL_BSL_LT_PIC_DEPOSIT_F)                                 &
              !                       .and. create_deposition_particles_on_a_grid )then
              !                print *, "-- DEBUG CHECK 876475645"
              !              end if

              ! Now we loop over grid points inside this flow cell:
              !
              ! points in the grid are of the form  d(i) = g_grid_d_min + (i-1) * h_g_grid_d,   i = 1, .. g_num_points_d
              ! (with  d = x, y, vx or vy  and  g_num_points_d = g_num_cells_d  or  g_num_cells_d+1  , depending on the periodicity)
              ! and this flow cell has the form     [flow_grid_d_min + (j-1) * h_flow_grid_d, flow_grid_d_min + j * h_flow_grid_d[
              ! so eta_d(i) is in this flow cell if  i_min <= i <= i_max
              ! where i_min = ceiling( (flow_grid_min - g_grid_d_min + (j-1)*h_flow_grid_d)/h_g_grid_d + 1 )
              ! and   i_max = ceiling( (flow_grid_min - g_grid_d_min + (j)  *h_flow_grid_d)/h_g_grid_d + 1 ) - 1

              i_min_x  = int(ceiling( (flow_grid_x_min  - g_grid_x_min  + (j_x-1)  * h_flow_grid_x) / h_g_grid_x + 1  ))
              i_min_y  = int(ceiling( (flow_grid_y_min  - g_grid_y_min  + (j_y-1)  * h_flow_grid_y) / h_g_grid_y + 1  ))
              i_min_vx = int(ceiling( (flow_grid_vx_min - g_grid_vx_min + (j_vx-1) * h_flow_grid_vx)/ h_g_grid_vx + 1 ))
              i_min_vy = int(ceiling( (flow_grid_vy_min - g_grid_vy_min + (j_vy-1) * h_flow_grid_vy)/ h_g_grid_vy + 1 ))

              i_max_x  = int(ceiling( (flow_grid_x_min  - g_grid_x_min  + (j_x)   * h_flow_grid_x) / h_g_grid_x + 1  )) - 1
              i_max_y  = int(ceiling( (flow_grid_y_min  - g_grid_y_min  + (j_y)   * h_flow_grid_y) / h_g_grid_y + 1  )) - 1
              i_max_vx = int(ceiling( (flow_grid_vx_min - g_grid_vx_min + (j_vx)  * h_flow_grid_vx)/ h_g_grid_vx + 1 )) - 1
              i_max_vy = int(ceiling( (flow_grid_vy_min - g_grid_vy_min + (j_vy)  * h_flow_grid_vy)/ h_g_grid_vy + 1 )) - 1

              i_min_x  = max(i_min_x,  1)
              i_min_y  = max(i_min_y,  1)
              i_min_vx = max(i_min_vx, 1)
              i_min_vy = max(i_min_vy, 1)

              i_max_x  = min(i_max_x,  g_num_points_x)
              i_max_y  = min(i_max_y,  g_num_points_y)
              i_max_vx = min(i_max_vx, g_num_points_vx)
              i_max_vy = min(i_max_vy, g_num_points_vy)

              do i_x = i_min_x, i_max_x
                x = g_grid_x_min + (i_x-1)*h_g_grid_x
                x_to_xk = x - x_k
                d1_x = d11 * x_to_xk
                d2_x = d21 * x_to_xk
                d3_x = d31 * x_to_xk
                d4_x = d41 * x_to_xk

                if( scenario == SLL_BSL_LT_PIC_DEPOSIT_F )then
                  ! find poisson (x-)cell containing this node, seen as a deposition particle, and relative position in the cell
                  tmp = ( x - p_group%space_mesh_2d%eta1_min ) / p_group%space_mesh_2d%delta_eta1
                  i_cell_x = int( tmp ) + 1
                  cell_offset_x = tmp - (i_cell_x-1)  ! between 0 and 1

                  ! SLL_ASSERT( abs(h_deposition_grid_x - h_g_grid_x) < 0.00001 * h_g_grid_x )
                end if

                do i_y = i_min_y, i_max_y
                  y = g_grid_y_min + (i_y-1)*h_g_grid_y
                  y_to_yk = y - y_k
                  d1_y = d12 * y_to_yk
                  d2_y = d22 * y_to_yk
                  d3_y = d32 * y_to_yk
                  d4_y = d42 * y_to_yk

                  if( scenario == SLL_BSL_LT_PIC_DEPOSIT_F )then
                    ! find poisson (y-)cell containing this node, seen as a deposition particle, and relative position in the cell
                    tmp = ( y - p_group%space_mesh_2d%eta2_min ) / p_group%space_mesh_2d%delta_eta2
                    i_cell_y = int( tmp ) + 1
                    cell_offset_y = tmp - (i_cell_y-1)  ! between 0 and 1

                    ! set the proper accumulator cell for the deposition
                    i_cell = i_cell_x + (i_cell_y-1) * p_group%space_mesh_2d%num_cells1   !  (see global_to_cell_offset)
                    charge_accumulator_cell => charge_accumulator%q_acc(i_cell)

                  end if

                  do i_vx = i_min_vx, i_max_vx
                    vx = g_grid_vx_min + (i_vx-1)*h_g_grid_vx
                    vx_to_vxk = vx - vx_k
                    d1_vx = d13 * vx_to_vxk
                    d2_vx = d23 * vx_to_vxk
                    d3_vx = d33 * vx_to_vxk
                    d4_vx = d43 * vx_to_vxk


                    do i_vy = i_min_vy, i_max_vy
                      vy = g_grid_vy_min + (i_vy-1)*h_g_grid_vy
                      vy_to_vyk = vy - vy_k
                      d1_vy = d14 * vy_to_vyk
                      d2_vy = d24 * vy_to_vyk
                      d3_vy = d34 * vy_to_vyk
                      d4_vy = d44 * vy_to_vyk

                      !> C.3 reconstruct the value of f at z_i = (x,y,vx,vy)
                      !>    for this we need z_t0, the position of the z_i at time = 0, using the affine backward flow

                      x_t0_to_xk_t0   = d1_x + d1_y + d1_vx + d1_vy
                      y_t0_to_yk_t0   = d2_x + d2_y + d2_vx + d2_vy
                      vx_t0_to_vxk_t0 = d3_x + d3_y + d3_vx + d3_vy
                      vy_t0_to_vyk_t0 = d4_x + d4_y + d4_vx + d4_vy

                      x_t0  = x_t0_to_xk_t0   + x_k_t0
                      y_t0  = y_t0_to_yk_t0   + y_k_t0
                      vx_t0 = vx_t0_to_vxk_t0 + vx_k_t0
                      vy_t0 = vy_t0_to_vyk_t0 + vy_k_t0

                      ! put back z_t0 inside domain (if needed) to enforce periodic boundary conditions
                      call periodic_correction(p_group,x_t0,y_t0)

                      ! now (x_t0, y_t0, vx_t0, vy_t0) is the (approx) position of the node z_i at time t=0
                      eta(1) = x_t0
                      eta(2) = y_t0
                      eta(3) = vx_t0
                      eta(4) = vy_t0

                      ! interpolation here may use sparse grid or splines, depending on the method chosen for p_group
                      reconstructed_f_value = p_group%bsl_lt_pic_4d_interpolate_value_of_remapped_f(eta)

                      ! [DEBUG]
                      if( .false. )then
                        if( scenario == SLL_BSL_LT_PIC_DEPOSIT_F )then
                          SLL_ASSERT( create_deposition_particles_on_a_grid )
                          print *, " -- DEBUG CHECK  --  8755648765  --  reconstructed_f_value = ", reconstructed_f_value
                        end if
                      end if

                      ! [DEBUG]
                      if( .false. )then
                        if( scenario == SLL_BSL_LT_PIC_WRITE_F_ON_GIVEN_GRID )then
                          print*, "[WRITE ON GIVEN GRID]    reconstructing "
                          print*, "on:                       eta = ", eta
                          print*, "value:  reconstructed_f_value = ", reconstructed_f_value
                        end if
                      end if

                      if( reconstructed_f_value /= 0 )then

                        if( scenario == SLL_BSL_LT_PIC_DEPOSIT_F )then

                          SLL_ASSERT( create_deposition_particles_on_a_grid )

                          reconstructed_charge = reconstructed_f_value * deposition_dvol * p_group%species%q

                          debug_count = debug_count + 1
                          debug_charge = debug_charge + deposition_dvol * p_group%species%q

                          if( .false. )then
                            print *, "[DEBUG] -- [deposit with] ", reconstructed_f_value * deposition_dvol * p_group%species%q
                            print *, "[DEBUG] -- reconstructed_charge = ", reconstructed_charge
                          end if
                          tmp1 = (1.0_f64 - cell_offset_x)
                          tmp2 = (1.0_f64 - cell_offset_y)

                          charge_accumulator_cell%q_sw = charge_accumulator_cell%q_sw             &
                                  + reconstructed_charge * tmp1 * tmp2

                          charge_accumulator_cell%q_se = charge_accumulator_cell%q_se             &
                                  + reconstructed_charge *  cell_offset_x * tmp2

                          charge_accumulator_cell%q_nw = charge_accumulator_cell%q_nw             &
                                  + reconstructed_charge * tmp1 *  cell_offset_y

                          charge_accumulator_cell%q_ne = charge_accumulator_cell%q_ne             &
                                  + reconstructed_charge *  cell_offset_x *  cell_offset_y

                          ! count the total charge
                          deposited_charge = deposited_charge + reconstructed_charge

                        else

                          if( scenario == SLL_BSL_LT_PIC_REMAP_F )then

                            SLL_ASSERT( p_group%remapped_f_interpolation_type == SLL_BSL_LT_PIC_REMAP_WITH_SPLINES )

                            SLL_ASSERT(i_x  <= p_group%remapping_cart_grid_number_nodes_x())
                            SLL_ASSERT(i_y  <= p_group%remapping_cart_grid_number_nodes_y())
                            SLL_ASSERT(i_vx <= p_group%remapping_cart_grid_number_nodes_vx())
                            SLL_ASSERT(i_vy <= p_group%remapping_cart_grid_number_nodes_vy())

                            ! store the nodal value on the temporary array
                            tmp_f_values_on_remapping_cart_grid(i_x,i_y,i_vx,i_vy) = reconstructed_f_value

                          else if( scenario == SLL_BSL_LT_PIC_WRITE_F_ON_GIVEN_GRID )then

                            SLL_ASSERT( i_x  <= size(given_array_2d,1) )
                            SLL_ASSERT( i_vx <= size(given_array_2d,2) )

                            given_array_2d(i_x,i_vx) = given_array_2d(i_x,i_vx)         &
                                    + reconstructed_f_value * h_g_grid_y * h_g_grid_vy

                          else

                            SLL_ERROR("bsl_lt_pic_4d_write_f_on_grid_or_deposit", "ahem... you should not be reading this :)")

                          end if
                        end if
                      else
                        print *, "654654535466545434564 -- ZERO VALUE !"
                      end if
                    ! this is the end of the (fourfold) loop on the grid nodes
                    end do
                  end do
                end do
              end do
            else
              SLL_ERROR("bsl_lt_pic_4d_write_f_on_grid_or_deposit", "broken test on scenarios -- you should not be reading this :)")
            end if
          ! and this is the end of (fourfold) loop on the flow cells
          end do
        end do
      end do
    end do

    !> in the remapping case, last step is to compute the new remapping coefs from the nodal values on the remapping grid
    if( scenario == SLL_BSL_LT_PIC_REMAP_F )then

      if( p_group%remapped_f_interpolation_type == SLL_BSL_LT_PIC_REMAP_WITH_SPLINES )then

        ! computing the new spline coefs
        call p_group%bsl_lt_pic_4d_compute_new_spline_coefs( tmp_f_values_on_remapping_cart_grid )

      else if( p_group%remapped_f_interpolation_type == SLL_BSL_LT_PIC_REMAP_WITH_SPARSE_GRIDS )then

        ! computing the new sparse grid (hierarchical) coefs
        call p_group%sparse_grid_interpolator%compute_hierarchical_surplus(   &
                tmp_f_values_on_remapping_sparse_grid                         &
           )
        p_group%remapped_f_sparse_grid_coefficients = tmp_f_values_on_remapping_sparse_grid

        ! [DEBUG]
        if( .false. )then
          print*, "CHECK - remap on sparse grid"
          print*, node_counter, p_group%sparse_grid_interpolator%size_basis
          do node_index = 1, p_group%sparse_grid_interpolator%size_basis
            print*, "reading: ", node_index, p_group%remapped_f_sparse_grid_coefficients(node_index)
          end do
        end if

      else
        SLL_ERROR("writing f on the remap grid", "broken test")
      end if
    end if


    !print *,  " [DEPOSIT_CHARGE] DEBUG  -- 9878768786877 --- debug_count  = ", debug_count
    !print *,  " [DEPOSIT_CHARGE] DEBUG  -- 9878768786877 --- g_num_points:  ", &
    !        g_num_points_x * g_num_points_y * g_num_points_vx * g_num_points_vy
    !print *,  " [DEPOSIT_CHARGE] DEBUG  -- 9878768786877 --- g_num_points modified :  ", &
    !        g_num_points_x* g_num_points_y * (g_num_points_vx -1) * (g_num_points_vy - 1)
    !
    !print *,  " [DEPOSIT_CHARGE] DEBUG  -- 9878768786877 --- debug_charge = ", debug_charge
    !print *,  " [DEPOSIT_CHARGE] DEBUG  -- --- (scenario == SLL_BSL_LT_PIC_DEPOSIT_F) = ", (scenario == SLL_BSL_LT_PIC_DEPOSIT_F)
    !print *,  " [DEPOSIT_CHARGE] DEBUG  -- 9878768786877 --- deposited_charge = ", deposited_charge

    if( (scenario == SLL_BSL_LT_PIC_DEPOSIT_F) .and. enforce_total_charge )then


      if( deposited_charge == 0 )then
        print *, "WARNING (76576537475) -- total deposited charge is zero, which is strange..."
        print *, "                      -- (no charge correction in this case) "
      else
        charge_correction_factor = target_total_charge / deposited_charge

        print *, "[Enforcing charge]: target_total_charge, deposited_charge, charge_correction_factor = ", &
          target_total_charge, deposited_charge, charge_correction_factor

        do i_cell_x = 1, p_group%space_mesh_2d%num_cells1
          do i_cell_y = 1, p_group%space_mesh_2d%num_cells2

            ! index of the Poisson cell
            i_cell = i_cell_x + (i_cell_y-1) * p_group%space_mesh_2d%num_cells1
            charge_accumulator_cell => charge_accumulator%q_acc(i_cell)
            charge_accumulator_cell%q_sw = charge_accumulator_cell%q_sw * charge_correction_factor
            charge_accumulator_cell%q_se = charge_accumulator_cell%q_se * charge_correction_factor
            charge_accumulator_cell%q_nw = charge_accumulator_cell%q_nw * charge_correction_factor
            charge_accumulator_cell%q_ne = charge_accumulator_cell%q_ne * charge_correction_factor

          end do
        end do
      end if
    end if

  end subroutine bsl_lt_pic_4d_write_f_on_grid_or_deposit



  ! <<onestep>> <<ALH>> utility function for finding the neighbours of a particle, used by [[ONESTEPMACRO]]. "dim"
  ! corresponds to one of x,y,vx,vy.

  subroutine onestep(dim,dim_t0,kprime,p_list,h_markers_dim)

    sll_int :: dim
    sll_real64 :: dim_t0
    sll_int32 :: neighbour

    type(sll_bsl_lt_pic_4d_particle), dimension(:), pointer,intent(in) :: p_list

    !sll_int32 :: ngb_dim_right_index
    !sll_int32 :: ngb_dim_left_index
    sll_real64 :: h_markers_dim
    sll_int32 :: kprime
    sll_int32 :: j,jumps

    ! <<up>> means that kprime needs to go up ie increase in coordinate

    logical :: up

    ! Move to a closer neighbour only if dim_t0 is not located in a cell of size h_markers_dim and with a left bound of
    ! dim_t0

    if(kprime == 0)return

    ! How many jumps do we need to do in that direction to reduce the distance 'dim_t0' to a minimum?

    jumps = int(abs(dim_t0/h_markers_dim))

    ! dim_t0 < 0 means that the virtual particle is at the left of kprime (dim_t0 is a relative coordinate). If dim_t0
    ! is negative, add one step to move to a positive relative signed distance between dim_t0 and kprime.

    up = .true.
    if(dim_t0 < 0) then
       jumps = jumps + 1
       up = .false.
    end if

    ! resulting signed distance between marker at t0 and kprime

    if(up) then
       dim_t0 = dim_t0 - jumps * h_markers_dim
    else
       dim_t0 = dim_t0 + jumps * h_markers_dim
    endif

    ! do as many jumps as required through the neighbour pointers in the given dimension (1:x, 2:y, 3:vx, 4:vy). kprime
    ! can become zero (ie fall out of the domain) in non-periodic dimensions.

    j = 1
    do while(j<=jumps .and. kprime/=0)

       ! going through neighbours

       select case (dim)
#define ALONG_X 1
       case(ALONG_X)
          if(up) then
             neighbour = p_list(kprime)%ngb_xright_index
          else
             neighbour = p_list(kprime)%ngb_xleft_index
          endif
#define ALONG_Y 2
       case(ALONG_Y)
          if(up) then
             neighbour = p_list(kprime)%ngb_yright_index
          else
             neighbour = p_list(kprime)%ngb_yleft_index
          endif
#define ALONG_VX 3
       case(ALONG_VX)
          if(up) then
             neighbour = p_list(kprime)%ngb_vxright_index
          else
             neighbour = p_list(kprime)%ngb_vxleft_index
          endif
#define ALONG_VY 4
       case(ALONG_VY)
          if(up) then
             neighbour = p_list(kprime)%ngb_vyright_index
          else
             neighbour = p_list(kprime)%ngb_vyleft_index
          endif
       case default
          SLL_ASSERT(.false.)
       end select

       ! The convention in [[bsl_lt_pic_4d_compute_new_particles]] is that if there is no neighbour then a neighbour
       ! index is equal to the particle index

       if(neighbour/=kprime) then
          kprime = neighbour
       else
          kprime = 0
       endif
       j = j + 1
    end do
  end subroutine onestep


  ! <<get_ltp_deformation_matrix>>
  !> Compute the coefficients of the particle 'deformation' matrix
  !! which approximates the Jacobian matrix of the exact backward flow
  !! (defined between the current time and that of the particle initialization -- or the last particle remapping)
  !! at the current particle position.
  !! It is computed in two steps:
  !!    * first we compute a finite-difference approximation of the forward Jacobian matrix by using the fact that the
  !!      initial (or remapped) particles were located on a cartesian grid in phase space. Specifically, the entries of
  !!      the forward Jacobian matrix (say, (J^n_k)_xy = d(F^n_x)/dy (z^0_k) -- with z is the phase space coordinate)
  !!      are approximated by finite differences: here, by  DF / 2*h_y
  !!      with DF = F^n_x(z^0_k + (0,h_y,0,0)) - F^n_x(z^0_k - (0,h_y,0,0))
  !!    * then we invert that approximated forward Jacobian matrix
  !!
  !! Note: when computing finite differences in a periodic dimension (say x), one must be careful since two values of F_x
  !!    can be close in the periodic interval but distant by almost L_x in the (stored) [0,L_x[ representation.
  !!    To account for this (frequent) phenomenon we do the following:
  !!    when the difference DF (see example above) is larger than L_x/2, we make it smaller by adding +/- L_x to it.
  !!    Note that here we could very well be making the slope smaller than what it should actually be: indeed if the function
  !!    F^n_x is having a steep slope at z^0_k which adds one (or several) periods L_x to DF then our modification will
  !!    artificially lower the slope. But this is impossible to prevent in full generality (indeed: a steep slope crossing the
  !!    0 or L_x value will be lowered anyway in the [0,L_x[ representation) and should not be frequent (indeed it only happens
  !!    when F^n_x has high derivatives and the particle grid is coarse, which will lead to bad approximations anyhow).

  subroutine get_ltp_deformation_matrix (       &
                        p_group,                &
                        k,                      &
                        mesh_period_x,          &
                        mesh_period_y,          &
                        h_markers_x,              &
                        h_markers_y,              &
                        h_markers_vx,             &
                        h_markers_vy,             &
                        x_k, y_k,               &
                        vx_k, vy_k,             &
                        d11,d12,d13,d14,        &
                        d21,d22,d23,d24,        &
                        d31,d32,d33,d34,        &
                        d41,d42,d43,d44         &
                        )

        class(sll_bsl_lt_pic_4d_group),intent(inout) :: p_group
        sll_int32, intent(in) :: k

        sll_real64, intent(in)  :: mesh_period_x
        sll_real64, intent(in)  :: mesh_period_y

        sll_real64, intent(in)  :: h_markers_x
        sll_real64, intent(in)  :: h_markers_y
        sll_real64, intent(in)  :: h_markers_vx
        sll_real64, intent(in)  :: h_markers_vy

        sll_real64, intent(out) :: x_k, y_k         ! particle center in physical space
        sll_real64, intent(out) :: vx_k, vy_k       ! particle center in velocity space
        sll_real64, intent(out) :: d11,d12,d13,d14  ! coefs of matrix D (backward Jacobian)
        sll_real64, intent(out) :: d21,d22,d23,d24
        sll_real64, intent(out) :: d31,d32,d33,d34
        sll_real64, intent(out) :: d41,d42,d43,d44

        sll_int32   :: k_ngb
        sll_real64  :: x_k_left,  x_k_right
        sll_real64  :: y_k_left,  y_k_right
        sll_real64  :: vx_k_left, vx_k_right
        sll_real64  :: vy_k_left, vy_k_right
        sll_real64, dimension(3)  :: coords

        sll_real64  :: j11,j12,j13,j14   ! coefs of matrix J = D^-1 (forward Jacobian)
        sll_real64  :: j21,j22,j23,j24
        sll_real64  :: j31,j32,j33,j34
        sll_real64  :: j41,j42,j43,j44
        sll_real64  :: factor, det_J, inv_det_J

        logical domain_is_x_periodic
        logical domain_is_y_periodic

        domain_is_x_periodic = p_group%domain_is_periodic(1)
        domain_is_y_periodic = p_group%domain_is_periodic(2)

        coords(:) = 0.0_f64

        coords = p_group%get_x(k)
        x_k = coords(1)
        y_k = coords(2)

        coords = p_group%get_v(k)
        vx_k = coords(1)
        vy_k = coords(2)

        ! Compute the forward Jacobian matrix J_k

        ! ------   d/d_x terms
        factor = 1./(2*h_markers_x)

        k_ngb  = p_group%markers_list(k)%ngb_xright_index
        if( k_ngb == k )then
           ! no right neighbor is available, use a non-centered finite difference
           factor = 2*factor
           x_k_right   = x_k
           y_k_right   = y_k
           vx_k_right  = vx_k
           vy_k_right  = vy_k
        else
            coords = p_group%get_x(k_ngb)
            x_k_right = coords(1)
            y_k_right = coords(2)
            coords = p_group%get_v(k_ngb)
            vx_k_right = coords(1)
            vy_k_right = coords(2)

            if( domain_is_x_periodic .and. x_k_right < x_k - 0.5*mesh_period_x ) x_k_right = x_k_right + mesh_period_x
            if( domain_is_x_periodic .and. x_k_right > x_k + 0.5*mesh_period_x ) x_k_right = x_k_right - mesh_period_x
            if( domain_is_y_periodic .and. y_k_right < y_k - 0.5*mesh_period_y ) y_k_right = y_k_right + mesh_period_y
            if( domain_is_y_periodic .and. y_k_right > y_k + 0.5*mesh_period_y ) y_k_right = y_k_right - mesh_period_y
        end if


        k_ngb  = p_group%markers_list(k)%ngb_xleft_index
        if( k_ngb == k )then
           ! no left neighbor is available, use a non-centered finite difference
           factor = 2*factor
           x_k_left   = x_k
           y_k_left   = y_k
           vx_k_left  = vx_k
           vy_k_left  = vy_k
        else
            coords = p_group%get_x(k_ngb)
            x_k_left = coords(1)
            y_k_left = coords(2)
            coords = p_group%get_v(k_ngb)
            vx_k_left = coords(1)
            vy_k_left = coords(2)

            if( domain_is_x_periodic .and. x_k_left < x_k - 0.5*mesh_period_x ) x_k_left = x_k_left + mesh_period_x
            if( domain_is_x_periodic .and. x_k_left > x_k + 0.5*mesh_period_x ) x_k_left = x_k_left - mesh_period_x
            if( domain_is_y_periodic .and. y_k_left < y_k - 0.5*mesh_period_y ) y_k_left = y_k_left + mesh_period_y
            if( domain_is_y_periodic .and. y_k_left > y_k + 0.5*mesh_period_y ) y_k_left = y_k_left - mesh_period_y
        end if

           j11 = factor * ( x_k_right  - x_k_left  )
           j21 = factor * ( y_k_right  - y_k_left  )
           j31 = factor * ( vx_k_right - vx_k_left )
           j41 = factor * ( vy_k_right - vy_k_left )

        ! ------   d/d_y terms
        factor = 1./(2*h_markers_y)

        k_ngb  = p_group%markers_list(k)%ngb_yright_index
        if( k_ngb == k )then
           ! no right neighbor is available, use a non-centered finite difference
           factor = 2*factor
           x_k_right   = x_k
           y_k_right   = y_k
           vx_k_right  = vx_k
           vy_k_right  = vy_k
        else
            coords = p_group%get_x(k_ngb)
            x_k_right = coords(1)
            y_k_right = coords(2)
            coords = p_group%get_v(k_ngb)
            vx_k_right = coords(1)
            vy_k_right = coords(2)

            if( domain_is_x_periodic .and. x_k_right < x_k - 0.5*mesh_period_x ) x_k_right = x_k_right + mesh_period_x
            if( domain_is_x_periodic .and. x_k_right > x_k + 0.5*mesh_period_x ) x_k_right = x_k_right - mesh_period_x
            if( domain_is_y_periodic .and. y_k_right < y_k - 0.5*mesh_period_y ) y_k_right = y_k_right + mesh_period_y
            if( domain_is_y_periodic .and. y_k_right > y_k + 0.5*mesh_period_y ) y_k_right = y_k_right - mesh_period_y
        end if

        k_ngb  = p_group%markers_list(k)%ngb_yleft_index
        if( k_ngb == k )then
           ! no left neighbor is available, use a non-centered finite difference
           factor = 2*factor
           x_k_left   = x_k
           y_k_left   = y_k
           vx_k_left  = vx_k
           vy_k_left  = vy_k
        else
            coords = p_group%get_x(k_ngb)
            x_k_left = coords(1)
            y_k_left = coords(2)
            coords = p_group%get_v(k_ngb)
            vx_k_left = coords(1)
            vy_k_left = coords(2)

            if( domain_is_x_periodic .and. x_k_left < x_k - 0.5*mesh_period_x ) x_k_left = x_k_left + mesh_period_x
            if( domain_is_x_periodic .and. x_k_left > x_k + 0.5*mesh_period_x ) x_k_left = x_k_left - mesh_period_x
            if( domain_is_y_periodic .and. y_k_left < y_k - 0.5*mesh_period_y ) y_k_left = y_k_left + mesh_period_y
            if( domain_is_y_periodic .and. y_k_left > y_k + 0.5*mesh_period_y ) y_k_left = y_k_left - mesh_period_y
        end if

           j12 = factor * ( x_k_right  - x_k_left  )
           j22 = factor * ( y_k_right  - y_k_left  )
           j32 = factor * ( vx_k_right - vx_k_left )
           j42 = factor * ( vy_k_right - vy_k_left )


        ! ------   d/d_vx terms
        factor = 1./(2*h_markers_vx)

        k_ngb  = p_group%markers_list(k)%ngb_vxright_index
        if( k_ngb == k )then
           ! no right neighbor is available, use a non-centered finite difference
           factor = 2*factor
           x_k_right   = x_k
           y_k_right   = y_k
           vx_k_right  = vx_k
           vy_k_right  = vy_k
        else
            coords = p_group%get_x(k_ngb)
            x_k_right = coords(1)
            y_k_right = coords(2)
            coords = p_group%get_v(k_ngb)
            vx_k_right = coords(1)
            vy_k_right = coords(2)

            if( domain_is_x_periodic .and. x_k_right < x_k - 0.5*mesh_period_x ) x_k_right = x_k_right + mesh_period_x
            if( domain_is_x_periodic .and. x_k_right > x_k + 0.5*mesh_period_x ) x_k_right = x_k_right - mesh_period_x
            if( domain_is_y_periodic .and. y_k_right < y_k - 0.5*mesh_period_y ) y_k_right = y_k_right + mesh_period_y
            if( domain_is_y_periodic .and. y_k_right > y_k + 0.5*mesh_period_y ) y_k_right = y_k_right - mesh_period_y
        end if

        k_ngb  = p_group%markers_list(k)%ngb_vxleft_index
        if( k_ngb == k )then
            ! no left neighbor is available, use a non-centered finite difference
            factor = 2*factor
            x_k_left   = x_k
            y_k_left   = y_k
            vx_k_left  = vx_k
            vy_k_left  = vy_k
        else
            coords = p_group%get_x(k_ngb)
            x_k_left = coords(1)
            y_k_left = coords(2)
            coords = p_group%get_v(k_ngb)
            vx_k_left = coords(1)
            vy_k_left = coords(2)

            if( domain_is_x_periodic .and. x_k_left < x_k - 0.5*mesh_period_x ) x_k_left = x_k_left + mesh_period_x
            if( domain_is_x_periodic .and. x_k_left > x_k + 0.5*mesh_period_x ) x_k_left = x_k_left - mesh_period_x
            if( domain_is_y_periodic .and. y_k_left < y_k - 0.5*mesh_period_y ) y_k_left = y_k_left + mesh_period_y
            if( domain_is_y_periodic .and. y_k_left > y_k + 0.5*mesh_period_y ) y_k_left = y_k_left - mesh_period_y
        end if

           j13 = factor * ( x_k_right  - x_k_left  )
           j23 = factor * ( y_k_right  - y_k_left  )
           j33 = factor * ( vx_k_right - vx_k_left )
           j43 = factor * ( vy_k_right - vy_k_left )


        ! ------   d/d_vy terms
        factor = 1./(2*h_markers_vy)

        k_ngb  = p_group%markers_list(k)%ngb_vyright_index
        if( k_ngb == k )then
           ! no right neighbor is available, use a non-centered finite difference
           factor = 2*factor
           x_k_right   = x_k
           y_k_right   = y_k
           vx_k_right  = vx_k
           vy_k_right  = vy_k
        else
            coords = p_group%get_x(k_ngb)
            x_k_right = coords(1)
            y_k_right = coords(2)
            coords = p_group%get_v(k_ngb)
            vx_k_right = coords(1)
            vy_k_right = coords(2)

            if( domain_is_x_periodic .and. x_k_right < x_k - 0.5*mesh_period_x ) x_k_right = x_k_right + mesh_period_x
            if( domain_is_x_periodic .and. x_k_right > x_k + 0.5*mesh_period_x ) x_k_right = x_k_right - mesh_period_x
            if( domain_is_y_periodic .and. y_k_right < y_k - 0.5*mesh_period_y ) y_k_right = y_k_right + mesh_period_y
            if( domain_is_y_periodic .and. y_k_right > y_k + 0.5*mesh_period_y ) y_k_right = y_k_right - mesh_period_y
        end if

        k_ngb  = p_group%markers_list(k)%ngb_vyleft_index
        if( k_ngb == k )then
            ! no left neighbor is available, use a non-centered finite difference
            factor = 2*factor
            x_k_left   = x_k
            y_k_left   = y_k
            vx_k_left  = vx_k
            vy_k_left  = vy_k
        else
            coords = p_group%get_x(k_ngb)
            x_k_left = coords(1)
            y_k_left = coords(2)
            coords = p_group%get_v(k_ngb)
            vx_k_left = coords(1)
            vy_k_left = coords(2)

            if( domain_is_x_periodic .and. x_k_left < x_k - 0.5*mesh_period_x ) x_k_left = x_k_left + mesh_period_x
            if( domain_is_x_periodic .and. x_k_left > x_k + 0.5*mesh_period_x ) x_k_left = x_k_left - mesh_period_x
            if( domain_is_y_periodic .and. y_k_left < y_k - 0.5*mesh_period_y ) y_k_left = y_k_left + mesh_period_y
            if( domain_is_y_periodic .and. y_k_left > y_k + 0.5*mesh_period_y ) y_k_left = y_k_left - mesh_period_y
        end if

           j14 = factor * ( x_k_right  - x_k_left  )
           j24 = factor * ( y_k_right  - y_k_left  )
           j34 = factor * ( vx_k_right - vx_k_left )
           j44 = factor * ( vy_k_right - vy_k_left )


           ! Compute D_k the inverse of J_k

        det_J = j11 * j22 * j33 * j44 &
              + j11 * j23 * j34 * j42 &
               + j11 * j24 * j32 * j43 &
              + j12 * j21 * j34 * j43 &
              + j12 * j23 * j31 * j44 &
              + j12 * j24 * j33 * j41 &
              + j13 * j21 * j32 * j44 &
              + j13 * j22 * j34 * j41 &
              + j13 * j24 * j31 * j42 &
              + j14 * j21 * j33 * j42 &
              + j14 * j22 * j31 * j43 &
              + j14 * j23 * j32 * j41 &
              - j11 * j22 * j34 * j43 &
              - j11 * j23 * j32 * j44 &
              - j11 * j24 * j33 * j42 &
              - j12 * j21 * j33 * j44 &
              - j12 * j23 * j34 * j41 &
              - j12 * j24 * j31 * j43 &
              - j13 * j21 * j34 * j42 &
              - j13 * j22 * j31 * j44 &
              - j13 * j24 * j32 * j41 &
              - j14 * j21 * j32 * j43 &
              - j14 * j22 * j33 * j41 &
              - j14 * j23 * j31 * j42
        !print*,'det  =',det
        inv_det_J = 1./det_J

        d11 = j22 * j33 * j44 &
            + j23 * j34 * j42 &
            + j24 * j32 * j43 &
            - j22 * j34 * j43 &
            - j23 * j32 * j44 &
            - j24 * j33 * j42
        d11 = inv_det_J * d11

        d12 = j12 * j34 * j43 &
            + j13 * j32 * j44 &
            + j14 * j33 * j42 &
            - j12 * j33 * j44 &
            - j13 * j34 * j42 &
            - j14 * j32 * j43
        d12 = inv_det_J * d12

        d13 = j12 * j23 * j44 &
            + j13 * j24 * j42 &
            + j14 * j22 * j43 &
            - j12 * j24 * j43 &
            - j13 * j22 * j44 &
            - j14 * j23 * j42
        d13 = inv_det_J * d13

        d14 = j12 * j24 * j33 &
            + j13 * j22 * j34 &
            + j14 * j23 * j32 &
            - j12 * j23 * j34 &
            - j13 * j24 * j32 &
            - j14 * j22 * j33
        d14 = inv_det_J * d14

        d21 = j21 * j34 * j43 &
            + j23 * j31 * j44 &
            + j24 * j33 * j41 &
            - j21 * j33 * j44 &
            - j23 * j34 * j41 &
            - j24 * j31 * j43
        d21 = inv_det_J * d21

        d22 = j11 * j33 * j44 &
            + j13 * j34 * j41 &
            + j14 * j31 * j43 &
            - j11 * j34 * j43 &
            - j13 * j31 * j44 &
            - j14 * j33 * j41
        d22 = inv_det_J * d22

        d23 = j11 * j24 * j43 &
            + j13 * j21 * j44 &
            + j14 * j23 * j41 &
            - j11 * j23 * j44 &
            - j13 * j24 * j41 &
            - j14 * j21 * j43
        d23 = inv_det_J * d23

        d24 = j11 * j23 * j34 &
            + j13 * j24 * j31 &
            + j14 * j21 * j33 &
            - j11 * j24 * j33 &
            - j13 * j21 * j34 &
            - j14 * j23 * j31
        d24 = inv_det_J * d24

        d31 = j21 * j32 * j44 &
            + j22 * j34 * j41 &
            + j24 * j31 * j42 &
            - j21 * j34 * j42 &
            - j22 * j31 * j44 &
            - j24 * j32 * j41
        d31 = inv_det_J * d31

        d32 = j11 * j34 * j42 &
            + j12 * j31 * j44 &
            + j14 * j32 * j41 &
            - j11 * j32 * j44 &
            - j12 * j34 * j41 &
            - j14 * j31 * j42
        d32 = inv_det_J * d32

        d33 = j11 * j22 * j44 &
            + j12 * j24 * j41 &
            + j14 * j21 * j42 &
            - j11 * j24 * j42 &
            - j12 * j21 * j44 &
            - j14 * j22 * j41
        d33 = inv_det_J * d33

        d34 = j11 * j24 * j32 &
            + j12 * j21 * j34 &
            + j14 * j22 * j31 &
            - j11 * j22 * j34 &
            - j12 * j24 * j31 &
            - j14 * j21 * j32
        d34 = inv_det_J * d34

        d41 = j21 * j33 * j42 &
            + j22 * j31 * j43 &
            + j23 * j32 * j41 &
            - j21 * j32 * j43 &
            - j22 * j33 * j41 &
            - j23 * j31 * j42
        d41 = inv_det_J * d41

        d42 = j11 * j32 * j43 &
            + j12 * j33 * j41 &
            + j13 * j31 * j42 &
            - j11 * j33 * j42 &
            - j12 * j31 * j43 &
            - j13 * j32 * j41
        d42 = inv_det_J * d42

        d43 = j11 * j23 * j42 &
            + j12 * j21 * j43 &
            + j13 * j22 * j41 &
            - j11 * j22 * j43 &
            - j12 * j23 * j41 &
            - j13 * j21 * j42
        d43 = inv_det_J * d43

        d44 = j11 * j22 * j33 &
            + j12 * j23 * j31 &
            + j13 * j21 * j32 &
            - j11 * j23 * j32 &
            - j12 * j21 * j33 &
            - j13 * j22 * j31
        d44 = inv_det_J * d44

end subroutine get_ltp_deformation_matrix


! puts the point (x,y) back into the computational domain if periodic in x or y (or both)
! otherwise, does nothing
subroutine periodic_correction(p_group, x, y)
    class(sll_bsl_lt_pic_4d_group),  intent(in)    :: p_group
    sll_real64, intent(inout) :: x
    sll_real64, intent(inout) :: y
    type(sll_cartesian_mesh_2d), pointer :: mesh

    mesh => p_group%space_mesh_2d

    if( p_group%domain_is_periodic(1)                        &
        .and.                                               &
        ( (x < mesh%eta1_min) .or. (x >= mesh%eta1_max) )   &
      ) then
          call apply_periodic_bc_x( mesh, x)
    end if
    if( p_group%domain_is_periodic(2)                        &
        .and.                                               &
        ( (y < mesh%eta2_min) .or. (y >= mesh%eta2_max) )   &
      ) then
          call apply_periodic_bc_y( mesh, y)
    end if
  end subroutine periodic_correction


  !----------------------------------------------------------------------------
  ! Destructor
  subroutine sll_bsl_lt_pic_4d_group_delete(particle_group)
    class(sll_bsl_lt_pic_4d_group), pointer :: particle_group
    sll_int32 :: ierr

    if(.not. associated(particle_group) ) then
       print *, 'sll_bsl_lt_pic_4d_group_delete(): ERROR (9087987578996), passed group was not associated.'
    end if
    SLL_DEALLOCATE(particle_group%markers_list, ierr)
    SLL_DEALLOCATE(particle_group%space_mesh_2d, ierr)
    SLL_DEALLOCATE(particle_group, ierr)

  end subroutine sll_bsl_lt_pic_4d_group_delete


end module sll_m_bsl_lt_pic_4d_group
