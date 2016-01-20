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
  use sll_m_bsl_lt_pic_4d_utilities !, only: int_list_element, int_list_element_ptr, add_element_in_int_list
  use sll_m_remapped_pic_utilities, only:           &
            x_is_in_domain_2d,                      &
            apply_periodic_bc_on_cartesian_mesh_2d, &
            get_inverse_matrix_with_given_size,     &
            get_4d_cell_containing_point
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
  sll_int32, parameter :: SLL_BSL_LT_PIC_BASIC = 0
  sll_int32, parameter :: SLL_BSL_LT_PIC_FLEXIBLE = 1

  ! types of deposition particles positions  //  and of flow markers
  sll_int32, parameter :: SLL_BSL_LT_PIC_STRUCTURED = 0
  sll_int32, parameter :: SLL_BSL_LT_PIC_UNSTRUCTURED  = 1

  ! types of deposition particles movement
  sll_int32, parameter :: SLL_BSL_LT_PIC_FIXED = 0
  sll_int32, parameter :: SLL_BSL_LT_PIC_PUSHED = 1

  ! old values:
  !  sll_int32, parameter :: SLL_BSL_LT_PIC_FIXED_GRID = 0
  !  sll_int32, parameter :: SLL_BSL_LT_PIC_TRANSPORTED_RANDOM = 1

  !> Group of @ref sll_bsl_lt_pic_4d_particle
  type, extends(sll_c_remapped_particle_group) :: sll_bsl_lt_pic_4d_group

    !> @name The markers (particles pushed forward, carry no weights) -- structured case
    !> @{
    sll_int32                                                   :: flow_markers_type        ! structured or unstructured
    ! structured flow markers always start from a cartesian grid in phase space (at initialization and remapping steps)
    sll_int32                                                   :: number_flow_markers_x
    sll_int32                                                   :: number_flow_markers_y
    sll_int32                                                   :: number_flow_markers_vx
    sll_int32                                                   :: number_flow_markers_vy
    sll_int32                                                   :: number_struct_flow_markers
    type(sll_cartesian_mesh_4d), pointer                        :: initial_markers_grid
    type(sll_bsl_lt_pic_4d_particle),   dimension(:), pointer   :: struct_markers_list
    ! When using unstructured flow markers, we store their indices in chained lists attached to the flow cells
    sll_int32                                                   :: nb_unstruct_markers_per_cell
    sll_int32                                                   :: max_nb_unstruct_markers
    sll_int32                                                   :: number_flow_markers      !< used for struct or unstruct markers

    type(int_list_element_ptr), dimension(:,:,:,:), allocatable :: unstruct_markers_in_flow_cell   !< to track markers in each cell
    type(int_list_element),     pointer                         :: unstruct_markers_outside_flow_grid   !< to track markers outside
    sll_real64, dimension(:,:), allocatable                     :: unstruct_markers_eta
    sll_real64, dimension(:,:), allocatable                     :: unstruct_markers_eta_at_remapping_time
    sll_int32,  dimension(:),   allocatable                     :: unstruct_markers_relevant_neighbor
    !> @}

    !> @name The flow grid (4d cartesian cells where the flow is linearized)
    !> @{
    type(sll_cartesian_mesh_4d), pointer    :: flow_grid
    sll_real64                              :: flow_grid_h    !< average step in grid
    sll_real64                              :: flow_grid_a1   !< anisotropic parameters: h_dim = flow_grid_h * flow_grid_a_dim
    sll_real64                              :: flow_grid_a2
    sll_real64                              :: flow_grid_a3
    sll_real64                              :: flow_grid_a4
    !> @}

    !> @name The physical mesh used eg in the Poisson solver
    !> @{
    type(sll_cartesian_mesh_2d), pointer    :: space_mesh_2d
    !> @}

    !> @name The deposition particles (will be created on the fly in each cell of the flow_grid, when depositing the charge)
    !> @{
    sll_int32                                 :: deposition_particles_type       !< basic (=first implementation) or flexible (new)
    sll_int32                                 :: deposition_particles_pos_type   !< structured  or  unstructured (random)
    sll_int32                                 :: deposition_particles_move_type  !< fixed (new at each step) or pushed (until remap)
    type(sll_cartesian_mesh_4d), pointer      :: deposition_particles_grid                 !< used if type = struct_grid
    sll_int32                                 :: number_deposition_particles        !< number of deposition particles
    sll_int32                                 :: number_moving_deposition_particles !< number of pushed deposition particles
    sll_real64, dimension(:,:), allocatable   :: deposition_particles_eta           !< used for if type = flexible
    sll_real64, dimension(:), allocatable     :: deposition_particles_weight        !< used for if type = flexible (weight = charge)
    !> @}

    !> This is how the deposition particles are used, depending on their type:
    !> if type = flexible, we allow deposition particles to be :
    !>  - move_type = fixed (always on the same grid, structured or not) or pushed (until remapping step)
    !>  - pos_type = structured or unstructured (random), when reset
    !>
    !> if move_type = fixed,
    !>    - positions and weights of deposition particles are reset in the deposition routine (called at each time step)
    !> if move_type = pushed,
    !>    - positions and weights of deposition particles are reset in the remapping routine (called at remapping step)
    !>      in particular, computing the weights can be done with direct interpolation (no bsl_lt_pic reconstruction)
    !>
    !> The "basic" type corresponds to fixed and structured deposition particles, first implemented.
    !> In this case, the deposition particles are not stored but only computed inside the "write_f_or_deposit" routine


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

    procedure :: bsl_lt_pic_4d_initialize_unstruct_markers         !> creates quasi-random distribution of unstructured markers
    procedure :: bsl_lt_pic_4d_prepare_unstruct_markers_for_flow_jacobians  !> build simplexes of relevant markers in each flow cell

    procedure :: reset_deposition_particles_coordinates
    procedure :: reset_deposition_particles_weights_with_bsl_reconstruction
    procedure :: reset_deposition_particles_weights_with_direct_interpolation
    procedure :: get_deposition_particle_charge_factor

    procedure :: bsl_lt_pic_4d_write_f_on_grid_or_deposit
    procedure :: bsl_lt_pic_4d_interpolate_value_of_remapped_f

    procedure :: update_flow_cell_lists_with_new_marker_position

    procedure :: get_ltp_deformation_matrix                           !> the local bwd flow using structured flow markers
    procedure :: get_deformation_matrix_from_unstruct_markers_in_cell !> the local bwd flow using unstructured flow markers
                                                                      ! name was get_affine_back_flow_in_cell_from_unstruct_markers
    procedure :: anisotropic_flow_grid_scalar_product
    procedure :: anisotropic_flow_grid_distance
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

    r = self%species%q * self%struct_markers_list(i)%weight

  end function bsl_lt_pic_4d_get_charge


  !----------------------------------------------------------------------------
  pure function bsl_lt_pic_4d_get_mass( self, i ) result( r )
    class( sll_bsl_lt_pic_4d_group ), intent( in ) :: self
    sll_int32                       , intent( in ) :: i
    sll_real64 :: r

    r = self%species%m * self%struct_markers_list(i)%weight

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
      ! returning an error would be nice but this is not possible in a pure function
      r = 1d30
      return
    end if

    if( i >= 1 .and. i <= self%number_flow_markers )then
      ! then the particle is a flow marker

      if( self%flow_markers_type == SLL_BSL_LT_PIC_STRUCTURED )then
        ! get x
        r(1) = self%space_mesh_2d%eta1_min + &
               self%space_mesh_2d%delta_eta1*(                            &
               real(self%struct_markers_list(i)%offset_x + self%struct_markers_list(i)%i_cell_x - 1, f64)      )
        ! get y
        r(2) = self%space_mesh_2d%eta2_min + self%space_mesh_2d%delta_eta2*( &
               real(self%struct_markers_list(i)%offset_y + self%struct_markers_list(i)%i_cell_y - 1, f64)      )
      else
        ! then self%flow_markers_type == SLL_BSL_LT_PIC_UNSTRUCTURED
        r(1) = self%unstruct_markers_eta(i, 1)
        r(2) = self%unstruct_markers_eta(i, 2)
      end if

    else if( i >= self%number_flow_markers + 1 .and. i <= self%number_flow_markers + self%number_moving_deposition_particles )then

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
      ! returning an error would be nice but this is not possible in a pure function
      r = 1d30
      return
    end if

    if( i >= 1 .and. i <= self%number_flow_markers )then
      ! then the particle is a flow marker

      if( self%flow_markers_type == SLL_BSL_LT_PIC_STRUCTURED )then
        ! get vx
        r(1) = self%struct_markers_list(i)%vx
        ! get vy
        r(2) = self%struct_markers_list(i)%vy
      else
        ! then self%flow_markers_type == SLL_BSL_LT_PIC_UNSTRUCTURED
        ! get vx
        r(1) = self%unstruct_markers_eta(i, 3)
        ! get vy
        r(2) = self%unstruct_markers_eta(i, 4)
      end if

    else if( i >= self%number_flow_markers + 1 .and. i <= self%number_flow_markers + self%number_moving_deposition_particles )then
      ! then the particle is a (flexible, pushed) deposition particle

      ! get vx
      r(1) = self%deposition_particles_eta(i - self%number_flow_markers, 3)
      ! get vy
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
    logical    ::  use_x_and_y_part

    if( i < 1 .or. i > self%number_flow_markers + self%number_moving_deposition_particles )then
      ! returning an error would be nice...
      i_out = 1000000000   ! will this be caught as an error ?
      return
    end if

    if( i >= 1 .and. i <= self%number_flow_markers )then
      ! then the particle is a flow marker

      use_x_and_y_part = .true.
      if( self%flow_markers_type == SLL_BSL_LT_PIC_STRUCTURED )then
        i_cell_x    = self%struct_markers_list(i)%i_cell_x
        i_cell_y    = self%struct_markers_list(i)%i_cell_y
        use_x_and_y_part = .false.
      else
        ! then self%flow_markers_type == SLL_BSL_LT_PIC_UNSTRUCTURED
        x_part = self%unstruct_markers_eta(i, 1)
        y_part = self%unstruct_markers_eta(i, 2)
      end if

    else if( i >= self%number_flow_markers + 1 .and. i <= self%number_flow_markers + self%number_moving_deposition_particles )then
      ! then the particle is a deposition particle

      x_part = self%deposition_particles_eta(i - self%number_flow_markers, 1)
      y_part = self%deposition_particles_eta(i - self%number_flow_markers, 2)

    end if

    if( use_x_and_y_part )then
      ! compute the x and y indices of the cell
      tmp = ( x_part - self%space_mesh_2d%eta1_min ) / self%space_mesh_2d%delta_eta1
      i_cell_x = int( tmp ) + 1
      ! NOTE: if we need the relative position within the cell (between 0 and 1), it is: dx = tmp - (i_cell_x - 1)
      tmp = ( y_part - self%space_mesh_2d%eta2_min ) / self%space_mesh_2d%delta_eta2
      i_cell_y = int( tmp ) + 1
      ! NOTE: if we need the relative position within the cell (between 0 and 1), it is: dy = tmp - (i_cell_y - 1)
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
    sll_int32                :: i_cell_x, i_cell_y
    sll_real32               :: offset_x, offset_y
    sll_real64               :: temp
    sll_real64, dimension(4) :: eta_marker
    sll_int32 :: old_j_x, old_j_y, old_j_vx, old_j_vy
    logical   :: marker_is_outside

    !print *, " *********** *********** *********** *********** *********** *********** *********** ***********         i = ", i
    !print *, " number_flow_markers = ", self%number_flow_markers
    !print *, " number_moving_deposition_particles = ", self%number_moving_deposition_particles
    !print *, " sum: ", self%number_flow_markers + self%number_moving_deposition_particles

    SLL_ASSERT( i >= 1 .and. i <= self%number_flow_markers + self%number_moving_deposition_particles )

    if( i >= 1 .and. i <= self%number_flow_markers )then
      ! then we set the physical coordinates of a flow marker
      ! (maybe change the name of the structure sll_bsl_lt_pic_4d_particle -> sll_bsl_lt_pic_4d_flow_marker ?)

      if( self%flow_markers_type == SLL_BSL_LT_PIC_STRUCTURED )then
        space_mesh_2d => self%space_mesh_2d
        particle => self%struct_markers_list(i)

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
      else
        ! then self%flow_markers_type == SLL_BSL_LT_PIC_UNSTRUCTURED

        ! we also need to update the flow cell's lists:
        ! store previous flow cell index
        eta_marker = self%unstruct_markers_eta(i, :)
        call get_4d_cell_containing_point(eta_marker, self%flow_grid, old_j_x, old_j_y, old_j_vx, old_j_vy, marker_is_outside)

        ! set x and y
        self%unstruct_markers_eta(i, 1) = x(1)
        self%unstruct_markers_eta(i, 2) = x(2)

        ! update flow cell lists
        call self%update_flow_cell_lists_with_new_marker_position(i, old_j_x, old_j_y, old_j_vx, old_j_vy)

      end if

    else if( i >= self%number_flow_markers + 1 .and. i <= self%number_flow_markers + self%number_moving_deposition_particles )then

      ! then the particle is a (pushed, flexible) deposition particle
      SLL_ASSERT( self%deposition_particles_move_type == SLL_BSL_LT_PIC_PUSHED )
      SLL_ASSERT( self%deposition_particles_type == SLL_BSL_LT_PIC_FLEXIBLE )
      SLL_ASSERT( self%number_moving_deposition_particles == self%number_deposition_particles )

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
    sll_real64, dimension(4) :: eta_marker
    sll_int32  :: old_j_x, old_j_y, old_j_vx, old_j_vy
    logical    :: marker_is_outside

    SLL_ASSERT( i >= 1 .and. i <= self%number_flow_markers + self%number_moving_deposition_particles )

    if( i >= 1 .and. i <= self%number_flow_markers )then
      ! then we set the velocity coordinates of a flow marker
      ! (maybe change the name of the structure sll_bsl_lt_pic_4d_particle -> sll_bsl_lt_pic_4d_flow_marker ?)


      if( self%flow_markers_type == SLL_BSL_LT_PIC_STRUCTURED )then
        particle => self%struct_markers_list(i)
        particle%vx = x(1)
        particle%vy = x(2)
      else
        ! then self%flow_markers_type == SLL_BSL_LT_PIC_UNSTRUCTURED

        ! we also need to update the flow cell's lists:
        ! store previous flow cell index
        eta_marker = self%unstruct_markers_eta(i, :)
        call get_4d_cell_containing_point(eta_marker, self%flow_grid, old_j_x, old_j_y, old_j_vx, old_j_vy, marker_is_outside)

        ! set vx and vy
        self%unstruct_markers_eta(i, 3) = x(1)
        self%unstruct_markers_eta(i, 4) = x(2)

        ! update flow cell lists
        call self%update_flow_cell_lists_with_new_marker_position(i, old_j_x, old_j_y, old_j_vx, old_j_vy)

      end if

    else if( i >= self%number_flow_markers + 1 .and. i <= self%number_flow_markers + self%number_moving_deposition_particles )then
      ! then we set the physical coordinates of a deposition particle

      ! then the particle is a (pushed, flexible) deposition particle
      SLL_ASSERT( self%deposition_particles_move_type == SLL_BSL_LT_PIC_PUSHED )
      SLL_ASSERT( self%deposition_particles_type == SLL_BSL_LT_PIC_FLEXIBLE )
      SLL_ASSERT( self%number_moving_deposition_particles == self%number_deposition_particles )

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

    print*, "Error (97658758) -- this subroutine is not implemented for sll_bsl_lt_pic_4d_group objects", i, s, storage_size(self)
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
  !> This subroutine places the deposition particles with the sepcified method. It does not compute the weights.
  !> The weights are computed in an external call to the 'write_f_on_grid_or_deposit' subroutine
  subroutine reset_deposition_particles_coordinates(self, rank)
    class(sll_bsl_lt_pic_4d_group), intent(inout)   :: self
    sll_int32, intent(in), optional                 :: rank
    sll_int32                                       :: this_rank
    sll_int32                                       :: i_x, i_y, i_vx, i_vy
    sll_int32                                       :: i_part
    sll_int32                                       :: i_dim
    sll_int64                                       :: sobol_seed
    sll_real64                                      :: rdn(4)
    sll_real64                                      :: eta_part(4)

    SLL_ASSERT( self%deposition_particles_type == SLL_BSL_LT_PIC_FLEXIBLE )

    if( self%deposition_particles_pos_type == SLL_BSL_LT_PIC_UNSTRUCTURED )then

      ! uses a quasi-random (Sobol) sequence
      ! initial value of the seed (it is incremented by one in the call to i8_sobol)
      if(present(rank))then
        this_rank = rank
      else
        this_rank = 0
      end if
      sobol_seed = int(10 + this_rank * self%number_deposition_particles, 8)

      do i_part = 1, self%number_deposition_particles
        ! Generate Sobol numbers on [0,1]
        call i8_sobol(int(4,8), sobol_seed, rdn)

        ! Transform rdn to the proper intervals
        do i_dim = 1, 4
          self%deposition_particles_eta(i_part, i_dim) =   self%remapping_grid_eta_min(i_dim) * (1 - rdn(i_dim)) &
                                                         + self%remapping_grid_eta_max(i_dim) * rdn(i_dim)
        end do
      end do

    else

      SLL_ASSERT( self%deposition_particles_pos_type == SLL_BSL_LT_PIC_STRUCTURED )

      ! fill in the grid of deposition particles (in every dimension we use the end nodes, even with periodic boundaries)
      i_part = 0
      do i_x = 1, self%deposition_particles_grid%num_cells1 + 1
        eta_part(1) = self%deposition_particles_grid%eta1_min + (i_x-1) * self%deposition_particles_grid%delta_eta1
        do i_y = 1, self%deposition_particles_grid%num_cells2 + 1
          eta_part(2) = self%deposition_particles_grid%eta2_min + (i_y-1) * self%deposition_particles_grid%delta_eta2
          do i_vx = 1, self%deposition_particles_grid%num_cells3 + 1
            eta_part(3) = self%deposition_particles_grid%eta3_min + (i_vx-1) * self%deposition_particles_grid%delta_eta3
            do i_vy = 1, self%deposition_particles_grid%num_cells4 + 1
              eta_part(4) = self%deposition_particles_grid%eta4_min + (i_vy-1) * self%deposition_particles_grid%delta_eta4
              i_part = i_part + 1
              self%deposition_particles_eta(i_part, :) = eta_part
            end do
          end do
        end do
      end do

      SLL_ASSERT( i_part == self%number_deposition_particles )

    end if

  end subroutine


  !----------------------------------------------------------------------------------
  ! do not change the position of the deposition particles, but compute (and set)
  ! their weights using the BSL_LT_PIC reconstruction

  subroutine reset_deposition_particles_weights_with_bsl_reconstruction(self)
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

    ! for basic deposition particles, the reconstruction is always done inside the write_f_on_grid routine
    SLL_ASSERT( self%deposition_particles_type == SLL_BSL_LT_PIC_FLEXIBLE )

    ! for pushed deposition particles, the reconstruction is always done at the remapping step, using a direct interpolation
    SLL_ASSERT( self%deposition_particles_type == SLL_BSL_LT_PIC_FIXED )

    ! reset the weights of the deposition particles, because maybe not every deposition particle weight will be set
    self%deposition_particles_weight = 0.0d0

    call self%bsl_lt_pic_4d_write_f_on_grid_or_deposit(void_charge_accumulator,             &
                                                          scenario,                         &
                                                          void_grid_4d,                     &
                                                          void_array_2d,                    &
                                                          dummy_total_charge,               &
                                                          enforce_total_charge              &
                                                          )

  end subroutine reset_deposition_particles_weights_with_bsl_reconstruction

  !----------------------------------------------------------------------------------
  ! do not change the position of the deposition particles, but compute (and set)
  ! their weights using a direct interpolation with the remapping tool (ie with flow = Id)

  subroutine reset_deposition_particles_weights_with_direct_interpolation(self)
    class( sll_bsl_lt_pic_4d_group ),           intent( inout ) :: self

    sll_real64    :: eta(4)
    sll_int32     :: i_part
    sll_real64    :: deposition_particle_charge_factor

    ! for basic deposition particles, the reconstruction is always done inside the write_f_on_grid routine
    SLL_ASSERT( self%deposition_particles_type == SLL_BSL_LT_PIC_FLEXIBLE )

    ! for fixed deposition particles, the reconstruction is done at each time step using a bsl_lt_pic reconstruction
    SLL_ASSERT( self%deposition_particles_type == SLL_BSL_LT_PIC_PUSHED )

    ! reset the weights of the deposition particles, because maybe not every deposition particle weight will be set
    self%deposition_particles_weight = 0.0d0

    deposition_particle_charge_factor = self%get_deposition_particle_charge_factor()

    do i_part = 1, self%number_deposition_particles
      eta = self%deposition_particles_eta(i_part, :)
      self%deposition_particles_weight(i_part) = deposition_particle_charge_factor                                    &
                                                  * self%bsl_lt_pic_4d_interpolate_value_of_remapped_f(eta)
    end do

  end subroutine reset_deposition_particles_weights_with_direct_interpolation


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


    if( self%deposition_particles_type == SLL_BSL_LT_PIC_BASIC )then
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
      SLL_ASSERT( self%deposition_particles_type == SLL_BSL_LT_PIC_FLEXIBLE )

      ! then we deposit the charge carried by the deposition particles with a quite standard technique

      ! but first: for flexible deposition particles of "fixed" type, we re-initialize their positions and weights
      if( self%deposition_particles_move_type == SLL_BSL_LT_PIC_FIXED )then
        call self%reset_deposition_particles_coordinates()
        call self%reset_deposition_particles_weights_with_bsl_reconstruction()
      end if

      do i_part = 1, self%number_deposition_particles

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

  function get_deposition_particle_charge_factor(p_group) result(val)
    class(sll_bsl_lt_pic_4d_group), intent(in)  :: p_group
    sll_real64 :: val

    sll_int32  :: nodes_number
    sll_real64 :: phase_space_volume

    nodes_number = p_group%number_deposition_particles

    phase_space_volume =    (p_group%remapping_grid_eta_max(4) - p_group%remapping_grid_eta_min(4))    &
                          * (p_group%remapping_grid_eta_max(3) - p_group%remapping_grid_eta_min(3))    &
                          * (p_group%remapping_grid_eta_max(2) - p_group%remapping_grid_eta_min(2))    &
                          * (p_group%remapping_grid_eta_max(1) - p_group%remapping_grid_eta_min(1))

    val = phase_space_volume * p_group%species%q / p_group%number_deposition_particles

  end function get_deposition_particle_charge_factor

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
        deposition_particles_pos_type,              &
        deposition_particles_move_type,             &
        number_deposition_particles,                &   ! (in a previous implementation this was only a lower bound)
        nb_deposition_particles_per_cell_x,         &
        nb_deposition_particles_per_cell_y,         &
        nb_deposition_particles_vx,                 &
        nb_deposition_particles_vy,                 &
        flow_markers_type,                          &
        number_flow_markers_x,                      &
        number_flow_markers_y,                      &
        number_flow_markers_vx,                     &
        number_flow_markers_vy,                     &
        nb_unstruct_markers_per_cell,               &
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
    sll_int32,                intent(in)  :: deposition_particles_pos_type
    sll_int32,                intent(in)  :: deposition_particles_move_type
    sll_int32,                intent(in)  :: number_deposition_particles
    sll_int32,                intent(in)  :: nb_deposition_particles_per_cell_x   !< used with deposition particles on a struct grid
    sll_int32,                intent(in)  :: nb_deposition_particles_per_cell_y   !< used with deposition particles on a struct grid
    sll_int32,                intent(in)  :: nb_deposition_particles_vx           !< used with deposition particles on a struct grid
    sll_int32,                intent(in)  :: nb_deposition_particles_vy           !< used with deposition particles on a struct grid
    sll_int32,                intent(in)  :: flow_markers_type
    sll_int32,                intent(in)  :: number_flow_markers_x
    sll_int32,                intent(in)  :: number_flow_markers_y
    sll_int32,                intent(in)  :: number_flow_markers_vx
    sll_int32,                intent(in)  :: number_flow_markers_vy
    sll_int32,                intent(in)  :: nb_unstruct_markers_per_cell
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

    sll_int32  :: effective_nb_deposition_particles_x
    sll_int32  :: effective_nb_deposition_particles_y
    sll_int32  :: effective_nb_deposition_particles_vx
    sll_int32  :: effective_nb_deposition_particles_vy

    sll_int32  :: deposition_particles_grid_num_cells_x
    sll_int32  :: deposition_particles_grid_num_cells_y
    sll_int32  :: deposition_particles_grid_num_cells_vx
    sll_int32  :: deposition_particles_grid_num_cells_vy

    sll_real64 :: deposition_particles_grid_x_min
    sll_real64 :: deposition_particles_grid_x_max
    sll_real64 :: deposition_particles_grid_y_min
    sll_real64 :: deposition_particles_grid_y_max
    sll_real64 :: deposition_particles_grid_vx_min
    sll_real64 :: deposition_particles_grid_vx_max
    sll_real64 :: deposition_particles_grid_vy_min
    sll_real64 :: deposition_particles_grid_vy_max

    sll_real64 :: h_deposition_particles_grid_x
    sll_real64 :: h_deposition_particles_grid_y

    logical    :: use_deposition_particles_grid
    sll_int32  :: effective_nb_deposition_particles_on_grid

    sll_int32  :: cst_int
    sll_real64 :: cst_real, ratio_vx, ratio_vy
    sll_real64 :: tmp !, tmp1, tmp2
    logical    :: derive_deposition_particles_grid_from_other_parameters

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
    !>    - A.1 flow grid: 4d cartesian cells where the flow is linearized
    !>    - A.2 if structured flow marker:
    !>      A.2.a list of marker coordinates (pushed forward)
    !>      A.2.b cartesian grid of initial markers
    !>    - A.3 if unstructured flow marker:
    !>      A.3.a allocate arrays for unstructured flow markers


    !> A.1 flow grid
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

    ! set parameters of the anisotropic distance in the flow_grid:
    !   we define  flow_grid_h  and  flow_grid_a1, .. , flow_grid_a4
    !   such that  flow_grid%delta_eta1 = flow_grid_a1 * flow_grid_h  and similarly in the other directions
    res%flow_grid_h = (1./4)*real(  res%flow_grid%delta_eta1  &
                                  + res%flow_grid%delta_eta2  &
                                  + res%flow_grid%delta_eta3  &
                                  + res%flow_grid%delta_eta4, f64)
    res%flow_grid_a1 = res%flow_grid%delta_eta1 / res%flow_grid_h
    res%flow_grid_a2 = res%flow_grid%delta_eta2 / res%flow_grid_h
    res%flow_grid_a3 = res%flow_grid%delta_eta3 / res%flow_grid_h
    res%flow_grid_a4 = res%flow_grid%delta_eta4 / res%flow_grid_h

    !> A.1-2 flow markers:

    res%flow_markers_type = flow_markers_type

    if( res%flow_markers_type == SLL_BSL_LT_PIC_STRUCTURED )then
      !> A.2.a list of marker coordinates (pushed forward)
      res%number_flow_markers_x  = number_flow_markers_x
      res%number_flow_markers_y  = number_flow_markers_y
      res%number_flow_markers_vx = number_flow_markers_vx
      res%number_flow_markers_vy = number_flow_markers_vy
      res%number_struct_flow_markers =   number_flow_markers_x  &
                                       * number_flow_markers_y  &
                                       * number_flow_markers_vx &
                                       * number_flow_markers_vy

      SLL_ALLOCATE( res%struct_markers_list(res%number_struct_flow_markers), ierr )

      !> A.2.b cartesian grid of initial markers
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

        res%number_flow_markers = res%number_struct_flow_markers

    else

      SLL_ASSERT( res% flow_markers_type == SLL_BSL_LT_PIC_UNSTRUCTURED )

      !>      A.3.a initialize parameters of unstructured flow markers

      if( nb_unstruct_markers_per_cell < 5 )then
        err_msg = "Error (875765786): we need at least 5 markers per cell (to create local simplexes defining the linearized flow)"
        SLL_ERROR(this_fun_name, err_msg)
      end if

      res%nb_unstruct_markers_per_cell = nb_unstruct_markers_per_cell
      res%max_nb_unstruct_markers = res%nb_unstruct_markers_per_cell * flow_grid_number_cells_x   &
                                                                     * flow_grid_number_cells_y   &
                                                                     * flow_grid_number_cells_vx  &
                                                                     * flow_grid_number_cells_vy

      !>      A.3.b allocate arrays for unstructured flow markers

      ! (we do not use the macro SLL_ALLOCATE here because the line is too long)
      allocate(res%unstruct_markers_in_flow_cell( flow_grid_number_cells_x,     &
                                                  flow_grid_number_cells_y,     &
                                                  flow_grid_number_cells_vx,    &
                                                  flow_grid_number_cells_vy )   &
               , stat=ierr)
      call test_error_code(ierr, 'Memory allocation Failure.', __FILE__, __LINE__)
      ! no need to do something for the list unstruct_markers_outside_flow_grid, it will just be nullified at the initialization

      SLL_ALLOCATE( res%unstruct_markers_eta(res%max_nb_unstruct_markers, 4) , ierr )
      SLL_ALLOCATE( res%unstruct_markers_eta_at_remapping_time(res%max_nb_unstruct_markers, 4) , ierr )
      SLL_ALLOCATE( res%unstruct_markers_relevant_neighbor(res%max_nb_unstruct_markers) , ierr )

      ! for the moment we are conservative and set this to the max number, to avoid the risk of forgetting some markers
      ! in the external push loops that would use a static (say, initial) number of particles (or markers)
      res%number_flow_markers = res%max_nb_unstruct_markers

    end if


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
      print*, "[", this_fun_name, "] - sparse grid levels for the remapping tool:", remapping_sparse_grid_max_levels
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

    end if


    !> D. discretization of the deposited f -- uses deposition particles which can be of several types (see comments on top of file)

    res%deposition_particles_type = deposition_particles_type

    if( res%deposition_particles_type == SLL_BSL_LT_PIC_BASIC )then
      res%deposition_particles_pos_type = SLL_BSL_LT_PIC_STRUCTURED
      res%deposition_particles_move_type = SLL_BSL_LT_PIC_FIXED
      ! Note: results should be the same as with FLEXIBLE particles of type STRUCTURED + FIXED,
      ! however BASIC particles are created on the fly (during deposition) and not stored
    else
      res%deposition_particles_pos_type  = deposition_particles_pos_type
      res%deposition_particles_move_type = deposition_particles_move_type
    end if

    use_deposition_particles_grid = ( res%deposition_particles_type == SLL_BSL_LT_PIC_BASIC   &
                                      .or. res%deposition_particles_pos_type == SLL_BSL_LT_PIC_STRUCTURED )

    derive_deposition_particles_grid_from_other_parameters = .false.

    if( use_deposition_particles_grid )then

      ! deposition particles will be reset on a cartesian grid denoted 'deposition_particles_grid'.

      if( derive_deposition_particles_grid_from_other_parameters )then

        ! this is how the deposition grid was created before January 17, 2016, just kept for memory (discard at some point)
        ! here we need parameters of a struct grid for the flow markers
        SLL_ASSERT( res%flow_markers_type == SLL_BSL_LT_PIC_STRUCTURED )

        ! In this implementation, we create the grid of deposition particles with the following properties:
        !
        !   - it matches the poisson grid in the sense that its nb of cells satisfies (with dg_np_d = deposition_particles_grid_num_points_d)
        !     dg_np_x ~ cst_int * p_group%space_mesh_2d%num_cells1  and
        !     dg_np_y ~ cst_int * p_group%space_mesh_2d%num_cells2  for some  cst_int
        !
        !   - the griding in vx, vy is obtained from that of the initial markers, using the linear scaling
        !     dg_np_vx / (dg_np_x * dg_np_y) = (approx) number_flow_markers_vx / (number_flow_markers_x * number_flow_markers_y)
        !     dg_np_vy / (dg_np_x * dg_np_y) = (approx) number_flow_markers_vy / (number_flow_markers_x * number_flow_markers_y)
        !
        !   - its number of nodes satisfies
        !     dg_np = dg_np_x * dg_np_y * dg_nc_vx * dg_nc_vy >= number_deposition_particles

        ! we will have  dg_np_vx ~ cst_int * cst_int * ratio_vx
        ratio_vx = real(res%number_flow_markers_vx * 1./ (res%number_flow_markers_x * res%number_flow_markers_y)          &
                                              * res%space_mesh_2d%num_cells1 * res%space_mesh_2d%num_cells2   ,f64)
        ! and           dg_np_vy ~ cst_int * cst_int * ratio_vy
        ratio_vy = real(res%number_flow_markers_vy * 1./ (res%number_flow_markers_x * res%number_flow_markers_y)          &
                                              * res%space_mesh_2d%num_cells1 * res%space_mesh_2d%num_cells2   ,f64)

        ! and           dg_np ~ cst_int * cst_int * cst_int * cst_int * num_cells1 * num_cells2 * ratio_vx * ratio_vy
        !                     >=  number_deposition_particles

        ! cst_real is the float approx of cst_int above
        cst_real = (real(number_deposition_particles, f64)   &
               / real( ratio_vx * ratio_vy * res%space_mesh_2d%num_cells1 * res%space_mesh_2d%num_cells2, f64)) ** (1./6)

        cst_int = int(ceiling( cst_real ))

        effective_nb_deposition_particles_x  = max( cst_int * res%space_mesh_2d%num_cells1, 2 )
        effective_nb_deposition_particles_y  = max( cst_int * res%space_mesh_2d%num_cells2, 2 )
        effective_nb_deposition_particles_vx = max( int(ceiling( cst_real * cst_real * ratio_vx )), 2 )
        effective_nb_deposition_particles_vy = max( int(ceiling( cst_real * cst_real * ratio_vy )), 2 )

        h_deposition_particles_grid_x = res%space_mesh_2d%delta_eta1 / cst_int    ! distance between two deposition particles in x dimension
        h_deposition_particles_grid_y = res%space_mesh_2d%delta_eta2 / cst_int    ! same in y

        effective_nb_deposition_particles_on_grid =   effective_nb_deposition_particles_x  * effective_nb_deposition_particles_y   &
                                                    * effective_nb_deposition_particles_vx * effective_nb_deposition_particles_vy
        print*, "[", this_fun_name, "] will use ", effective_nb_deposition_particles_on_grid, "deposition particles"
        print*, "[", this_fun_name, "] should be at least ", number_deposition_particles
        SLL_ASSERT( effective_nb_deposition_particles_on_grid >= number_deposition_particles )

      else

        ! this is now the default construction for the grid of deposition particles

        effective_nb_deposition_particles_x  = max( res%space_mesh_2d%num_cells1 * nb_deposition_particles_per_cell_x, 2 )
        effective_nb_deposition_particles_y  = max( res%space_mesh_2d%num_cells2 * nb_deposition_particles_per_cell_y, 2 )
        effective_nb_deposition_particles_vx = max( nb_deposition_particles_vx, 2 )
        effective_nb_deposition_particles_vy = max( nb_deposition_particles_vy, 2 )

        ! distance between two deposition particles will be:
        h_deposition_particles_grid_x = res%space_mesh_2d%delta_eta1 / nb_deposition_particles_per_cell_x
        h_deposition_particles_grid_y = res%space_mesh_2d%delta_eta2 / nb_deposition_particles_per_cell_y

        effective_nb_deposition_particles_on_grid =   effective_nb_deposition_particles_x  * effective_nb_deposition_particles_y   &
                                                    * effective_nb_deposition_particles_vx * effective_nb_deposition_particles_vy

        print*, "[", this_fun_name, "] (default build) will use ", effective_nb_deposition_particles_on_grid, "deposition particles"
      end if

      ! we create the deposition grid so that every deposition particle is _inside_ a poisson cell
      ! (so that we do not have to do something special for periodic boundary conditions)

      deposition_particles_grid_num_cells_x  = effective_nb_deposition_particles_x  - 1
      deposition_particles_grid_num_cells_y  = effective_nb_deposition_particles_y  - 1
      deposition_particles_grid_num_cells_vx = effective_nb_deposition_particles_vx - 1
      deposition_particles_grid_num_cells_vy = effective_nb_deposition_particles_vy - 1

      deposition_particles_grid_x_min  = res%space_mesh_2d%eta1_min + 0.5 * h_deposition_particles_grid_x
      deposition_particles_grid_x_max  = res%space_mesh_2d%eta1_max - 0.5 * h_deposition_particles_grid_x
      deposition_particles_grid_y_min  = res%space_mesh_2d%eta2_min + 0.5 * h_deposition_particles_grid_y
      deposition_particles_grid_y_max  = res%space_mesh_2d%eta2_max - 0.5 * h_deposition_particles_grid_y

      ! in velocity the bounds are those of the remapping grid
      deposition_particles_grid_vx_min = res%remapping_grid_eta_min(3)
      deposition_particles_grid_vx_max = res%remapping_grid_eta_max(3)
      deposition_particles_grid_vy_min = res%remapping_grid_eta_min(4)
      deposition_particles_grid_vy_max = res%remapping_grid_eta_max(4)

      res%deposition_particles_grid => new_cartesian_mesh_4d( deposition_particles_grid_num_cells_x,        &
                                                    deposition_particles_grid_num_cells_y,        &
                                                    deposition_particles_grid_num_cells_vx,       &
                                                    deposition_particles_grid_num_cells_vy,       &
                                                    deposition_particles_grid_x_min,   &
                                                    deposition_particles_grid_x_max,   &
                                                    deposition_particles_grid_y_min,   &
                                                    deposition_particles_grid_y_max,   &
                                                    deposition_particles_grid_vx_min,  &
                                                    deposition_particles_grid_vx_max,  &
                                                    deposition_particles_grid_vy_min,  &
                                                    deposition_particles_grid_vy_max   &
                                                   )

      tmp = abs(h_deposition_particles_grid_x - res%deposition_particles_grid%delta_eta1)
      SLL_ASSERT( tmp < 0.00001 * h_deposition_particles_grid_x )
      tmp = abs(h_deposition_particles_grid_y - res%deposition_particles_grid%delta_eta2)
      SLL_ASSERT( tmp < 0.00001 * h_deposition_particles_grid_y )

    end if

    if( res%deposition_particles_type == SLL_BSL_LT_PIC_FLEXIBLE )then

      ! Flexible deposition particles will have their weights and coordinates stored in the arrays initialized below
      !     they may be transported with the flow (like sdt particles) and re-initialized on remapping steps
      !     or stay on the grid and have new weights computed at each time step

      if( res%deposition_particles_pos_type == SLL_BSL_LT_PIC_STRUCTURED )then
        res%number_deposition_particles        = effective_nb_deposition_particles_on_grid
      else
        SLL_ASSERT( res%deposition_particles_pos_type == SLL_BSL_LT_PIC_UNSTRUCTURED )
        res%number_deposition_particles        = number_deposition_particles
      end if

      if( res%deposition_particles_move_type == SLL_BSL_LT_PIC_FIXED )then
        res%number_moving_deposition_particles = 0
      else
        SLL_ASSERT( res%deposition_particles_move_type == SLL_BSL_LT_PIC_PUSHED )
        res%number_moving_deposition_particles = res%number_deposition_particles
      end if

      SLL_ALLOCATE( res%deposition_particles_eta(res%number_deposition_particles, 4), ierr )
      SLL_ALLOCATE( res%deposition_particles_weight(res%number_deposition_particles), ierr )

    else
       SLL_ASSERT( res%deposition_particles_type == SLL_BSL_LT_PIC_BASIC )

       res%number_deposition_particles        = number_deposition_particles
       res%number_moving_deposition_particles = 0
    end if

    SLL_ASSERT( res%number_deposition_particles >= 0 )
    SLL_ASSERT( res%number_moving_deposition_particles*(res%number_moving_deposition_particles-res%number_deposition_particles)==0 )

    ! the variable "number_particles" is used in the interface to push particles
    res%number_particles = res%number_flow_markers + res%number_moving_deposition_particles

    !    else
    !      err_msg = "ahem, a test must be broken -- you should not be reading this :)"
    !      SLL_ERROR( this_fun_name, err_msg )
    !    end if

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
    print *, "bsl_lt_pic_4d_initializer -- step A: initialize the remapping tool"
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
    print *, "bsl_lt_pic_4d_initializer -- step B: initialize the flow markers"
    if( self%flow_markers_type == SLL_BSL_LT_PIC_STRUCTURED )then
      call self%bsl_lt_pic_4d_reset_markers_position()
      call self%bsl_lt_pic_4d_set_markers_connectivity()
    else
      SLL_ASSERT( self%flow_markers_type == SLL_BSL_LT_PIC_UNSTRUCTURED )
      call self%bsl_lt_pic_4d_initialize_unstruct_markers()
      call self%bsl_lt_pic_4d_prepare_unstruct_markers_for_flow_jacobians()
    end if

    !> C. if deposition particles are pushed, we initialize them now -- this requires that the remapping tool is initialized
    if( self%deposition_particles_type == SLL_BSL_LT_PIC_FLEXIBLE     &
        .and. self%deposition_particles_move_type == SLL_BSL_LT_PIC_PUSHED )then
      print *, "bsl_lt_pic_4d_initializer -- step C: initialize the deposition cells: "
      print *, "bsl_lt_pic_4d_initializer -- (C) will create ", self%number_deposition_particles, "deposition_particles..."
      ! if deposition particles are fixed, then they are initialized at each time step, in the deposition routine
      if(present(rank))then
        call self%reset_deposition_particles_coordinates(rank)
      else
        call self%reset_deposition_particles_coordinates()
      end if
      ! since the remapping tool has been set, computing the weights can be done with straightforward interpolation (flow = Id)
      call self%reset_deposition_particles_weights_with_direct_interpolation()
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

    !    sll_real64 :: total_density ! DEBUG

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
                p_group%struct_markers_list(k)%ngb_xleft_index = k
            else
                ! set the connectivity (in both directions) with left neighbor
                k_ngb = markers_indices(j_x-1,j_y,j_vx,j_vy)
                p_group%struct_markers_list(k)%ngb_xleft_index = k_ngb
                p_group%struct_markers_list(k_ngb)%ngb_xright_index = k
                if(j_x == number_flow_markers_x)then
                    if( p_group%domain_is_periodic(1) )then
                        ! set the connectivity (in both directions) with right neighbor
                        k_ngb = markers_indices(1,j_y,j_vx,j_vy)
                        p_group%struct_markers_list(k)%ngb_xright_index = k_ngb
                        p_group%struct_markers_list(k_ngb)%ngb_xleft_index = k
                    else
                        ! [neighbor index = own index] means: no neighbor
                        p_group%struct_markers_list(k)%ngb_xright_index = k
                    end if
                end if
            end if
            if(j_y == 1)then
                ! [neighbor index = own index] means: no neighbor
                ! in the y-periodic case this will be changed when dealing with the last marker in the y dimension
                p_group%struct_markers_list(k)%ngb_yleft_index = k
            else
                ! set the connectivity (in both directions) with left neighbor
                k_ngb = markers_indices(j_x,j_y-1,j_vx,j_vy)
                p_group%struct_markers_list(k)%ngb_yleft_index = k_ngb
                p_group%struct_markers_list(k_ngb)%ngb_yright_index = k
                if(j_y == number_flow_markers_y)then
                    if( p_group%domain_is_periodic(2) )then
                        ! set the connectivity (in both directions) with right neighbor
                        k_ngb = markers_indices(j_x,1,j_vx,j_vy)
                        p_group%struct_markers_list(k)%ngb_yright_index = k_ngb
                        p_group%struct_markers_list(k_ngb)%ngb_yleft_index = k
                    else
                        ! [neighbor index = own index] means: no neighbor
                        p_group%struct_markers_list(k)%ngb_yright_index = k
                    end if
                end if
            end if
            if(j_vx == 1)then
                ! [neighbor index = own index] means: no neighbor
                p_group%struct_markers_list(k)%ngb_vxleft_index = k
            else
                ! set the connectivity (in both directions) with left neighbor
                k_ngb = markers_indices(j_x,j_y,j_vx-1,j_vy)
                p_group%struct_markers_list(k)%ngb_vxleft_index = k_ngb
                p_group%struct_markers_list(k_ngb)%ngb_vxright_index = k
                if(j_vx == number_flow_markers_vx)then
                    ! [neighbor index = own index] means: no neighbor
                    p_group%struct_markers_list(k)%ngb_vxright_index = k
                end if
            end if
            if(j_vy == 1)then
                ! [neighbor index = own index] means: no neighbor
                p_group%struct_markers_list(k)%ngb_vyleft_index = k
            else
                ! set the connectivity (in both directions) with left neighbor
                k_ngb = markers_indices(j_x,j_y,j_vx,j_vy-1)
                p_group%struct_markers_list(k)%ngb_vyleft_index = k_ngb
                p_group%struct_markers_list(k_ngb)%ngb_vyright_index = k
                if(j_vy == number_flow_markers_vy)then
                    ! [neighbor index = own index] means: no neighbor
                    p_group%struct_markers_list(k)%ngb_vyright_index = k
                end if
            end if

          end do
        end do
      end do
    end do

  end subroutine bsl_lt_pic_4d_set_markers_connectivity


  !> initialize the markers on the initial (markers) grid
  subroutine bsl_lt_pic_4d_initialize_unstruct_markers( p_group )
    class(sll_bsl_lt_pic_4d_group),intent(inout) :: p_group

    type(int_list_element),  pointer   :: new_int_list_element, head
    sll_real64, dimension(4) :: flow_cell_eta_min
    sll_real64, dimension(4) :: flow_cell_eta_max

    logical    :: local_verbose
    sll_int64  :: sobol_seed
    sll_real64 :: rdn(4)
    sll_int32  :: ierr
    sll_int32  :: j_x, j_y, j_vx, j_vy
    sll_int32  :: i_dim
    sll_int32  :: i_marker
    sll_int32  :: i_local_marker
    sll_int32  :: nb_unstruct_markers_per_cell

    ! place a few markers in each flow cell with quasi-random sequences
    nb_unstruct_markers_per_cell = p_group%nb_unstruct_markers_per_cell

    i_marker = 1      ! index of marker in global list
    sobol_seed = int(675,8)  ! initial value of the seed (it is incremented by one in the call to i8_sobol)

    local_verbose = .true.

    do j_x = 1, p_group%flow_grid%num_cells1
      if( local_verbose )then
        print *, "loading unstructured flow markers in flow cells: j_x = ", j_x, "/", p_group%flow_grid%num_cells1, "..."
      end if
      flow_cell_eta_min(1) = p_group%flow_grid%eta1_min + (j_x-1) * p_group%flow_grid%delta_eta1
      flow_cell_eta_max(1) = p_group%flow_grid%eta1_min + (j_x)   * p_group%flow_grid%delta_eta1

      do j_y = 1, p_group%flow_grid%num_cells2
        flow_cell_eta_min(2) = p_group%flow_grid%eta2_min + (j_y-1) * p_group%flow_grid%delta_eta2
        flow_cell_eta_max(2) = p_group%flow_grid%eta2_min + (j_y)   * p_group%flow_grid%delta_eta2

        do j_vx = 1, p_group%flow_grid%num_cells3
          flow_cell_eta_min(3) = p_group%flow_grid%eta3_min + (j_vx-1) * p_group%flow_grid%delta_eta3
          flow_cell_eta_max(3) = p_group%flow_grid%eta3_min + (j_vx)   * p_group%flow_grid%delta_eta3

          do j_vy = 1, p_group%flow_grid%num_cells4
            flow_cell_eta_min(4) = p_group%flow_grid%eta4_min + (j_vy-1) * p_group%flow_grid%delta_eta4
            flow_cell_eta_max(4) = p_group%flow_grid%eta4_min + (j_vy)   * p_group%flow_grid%delta_eta4

            do i_local_marker = 1, nb_unstruct_markers_per_cell

              ! Generate 4 Sobol numbers on [0,1]
              call i8_sobol(int(4,8), sobol_seed, rdn)

              ! Transform rdn to the proper intervals
              do i_dim = 1, 4
                p_group%unstruct_markers_eta(i_marker, i_dim) =   flow_cell_eta_min(i_dim) * (1 - rdn(i_dim)) &
                                                                + flow_cell_eta_max(i_dim) * rdn(i_dim)
              end do

              ! increment the proper linked list (see other uses in this module)
              SLL_ALLOCATE( new_int_list_element, ierr )
              new_int_list_element%value = i_marker
              head => p_group%unstruct_markers_in_flow_cell(j_x,j_y,j_vx,j_vy)%pointed_element
              p_group%unstruct_markers_in_flow_cell(j_x,j_y,j_vx,j_vy)%pointed_element &
                    => add_element_in_int_list(head, new_int_list_element)

              ! increment the global marker index
              i_marker = i_marker + 1

            end do
          end do
        end do
      end do
    end do

    ! initialize the linked list of outside markers
    nullify(p_group%unstruct_markers_outside_flow_grid)

  end subroutine bsl_lt_pic_4d_initialize_unstruct_markers

  ! ------------------------------------------------------------------------------------------------------------------------------
  !> prepare a relevant set of 5 (D+1) markers to be able to later compute the flow Jacobian matrix in each cell
  !> after this routine, each cell will have a list of 5 markers defining a simplex of sufficient volume
  !>
  !> the rule for selecting 5 relevant markers (x^1, .. x^5) in a cell Omega_j is the following:
  !>    - the marker #1 is relevant if its distance to the cell's center C^j is > diam(Omega_j)/4.
  !>      Note: considering that   diam_d(Omega)_j = h * a_d   for  d = 1, .., D (here D=4)  we will use the anisotropic l2 distance
  !>        (|| x ||_a)^2 := <x, x>_a       where      <x, y >_a  := sum_{d = 1, .. D}  (x_d / a_d) * (y_d / a_d)
  !>        and test if
  !>        || ( (x^1_d - C^j_d) / diam_d(Omega_j) )_{d = 1, .. D} ||_a > h/4
  !>    - for i = 2, .. 5 the marker #i is relevant if its a-distance to the space < x^2-x^1, .. x^{i-1}-x^1 >  is > diam(Omega_j)/4
  !>      again in the anisotropic l2 distance.
  !>      To test this we loop over the markers x in the cell and compute their anisotropic l2 projection
  !>      in the affine space,   Px = x^1 + sum_l=2^{i-1} c_{l-1} (x^l - x^1)   defined by
  !>        < Px - x^1, x^m - x^1 >_a = < x - x^1, x^m - x^1 >_a      for    m = 2, .. i-1
  !>      (for i=2 we have Px = x^1) and test whether || x - Px ||_a > h/4
  !>      Note: if one marker was too close from a low-dimensional space we can discard it from further searches
  !>    - at each stage, if no marker is relevant we add a new marker with quasi-random position that makes it relevant in the cell
  !> Note: Maybe 1/4 is not the best choice.
  !> Note: We will also remove some markers in the cells if they are non-relevant and too many
  !>
  !> Note: in constructor we should have defined the parameters  flow_grid_a1, .. flow_grid_a4 and flow_grid_h  such that
  !> p_group%flow_grid%delta_eta1 = p_group%flow_grid_a1 * p_group%flow_grid_h
  !> p_group%flow_grid%delta_eta2 = p_group%flow_grid_a2 * p_group%flow_grid_h
  !> p_group%flow_grid%delta_eta3 = p_group%flow_grid_a3 * p_group%flow_grid_h
  !> p_group%flow_grid%delta_eta4 = p_group%flow_grid_a4 * p_group%flow_grid_h

  subroutine bsl_lt_pic_4d_prepare_unstruct_markers_for_flow_jacobians( p_group )
    class(sll_bsl_lt_pic_4d_group),intent(inout) :: p_group


    type(int_list_element),       pointer           :: list_of_marker_indices_to_be_discarded
    type(int_list_element),       pointer           :: new_int_list_element, head_int_list
    type(int_list_element),       pointer           :: new_aux_int_list_element
    type(marker_list_element),    pointer           :: list_of_markers_to_be_added
    type(marker_list_element),    pointer           :: new_marker_list_element, head_marker_list

    sll_real64, dimension(4)    :: flow_cell_eta_min
    sll_real64, dimension(4)    :: flow_cell_eta_mid
    sll_real64, dimension(4)    :: flow_cell_eta_max
    sll_real64, dimension(4)    :: eta_marker
    sll_real64, dimension(4)    :: eta_projected
    sll_real64, dimension(5,4)  :: eta_relevant_marker
    sll_real64, dimension(4)    :: eta_marker_to_remove
    sll_real64, dimension(4)    :: projection_coef
    sll_real64, dimension(4)    :: aux_b
    sll_real64, dimension(4,4)  :: aux_matrix
    sll_real64, dimension(4,4)  :: projection_matrix

    logical    :: relevant_marker_found
    !    logical    :: marker_found_and_removed
    logical    :: ok_flag
    logical    :: local_verbose
    logical    :: marker_is_outside
    sll_int32  :: ierr
    sll_int64  :: sobol_seed
    sll_real64 :: rdn(4)
    sll_real64 :: a_distance
    sll_int32  :: j_x, j_y, j_vx, j_vy
    sll_int32  :: old_j_x, old_j_y, old_j_vx, old_j_vy
    sll_int32  :: l, m
    sll_int32  :: i_dim
    sll_int32  :: i_relevant
    sll_int32  :: i_marker
    sll_int32  :: matrix_size
    sll_int32  :: nb_unrelevant_markers_in_this_cell
    sll_int32  :: nb_relevant_markers_in_this_cell
    sll_int32  :: i_first_relevant_marker
    sll_int32  :: i_last_relevant_marker
    sll_int32  :: nb_markers_created

    local_verbose = .false.

    SLL_ASSERT( p_group%flow_markers_type == SLL_BSL_LT_PIC_UNSTRUCTURED )

    ! list of markers that will be temporary stored for removal
    nullify( list_of_marker_indices_to_be_discarded )

    ! list of markers that will be temporary stored for insertion
    nullify( list_of_markers_to_be_added )

    ! reset the array: 0 means no relevant neighbor
    p_group%unstruct_markers_relevant_neighbor(:) = 0

    sobol_seed = int(876,8) ! seed for newly added markers (we could reuse the one from initial loading...)
    nb_markers_created = 0    ! to count

    do j_x = 1, p_group%flow_grid%num_cells1

      if( local_verbose )then
        print *, "[preparing the unstruct markers], A. searching/creating relevant markers in flow cells: j_x = ", &
                      j_x, "/", p_group%flow_grid%num_cells1, "..."
      end if

      flow_cell_eta_min(1) = p_group%flow_grid%eta1_min + (j_x-1) * p_group%flow_grid%delta_eta1
      flow_cell_eta_max(1) = p_group%flow_grid%eta1_min + (j_x)   * p_group%flow_grid%delta_eta1
      flow_cell_eta_mid(1) = flow_cell_eta_min(1) + 0.5*p_group%flow_grid%delta_eta1

      do j_y = 1, p_group%flow_grid%num_cells2
        flow_cell_eta_min(2) = p_group%flow_grid%eta2_min + (j_y-1) * p_group%flow_grid%delta_eta2
        flow_cell_eta_max(2) = p_group%flow_grid%eta2_min + (j_y)   * p_group%flow_grid%delta_eta2
        flow_cell_eta_mid(2) = flow_cell_eta_min(2) + 0.5*p_group%flow_grid%delta_eta2

        do j_vx = 1, p_group%flow_grid%num_cells3
          flow_cell_eta_min(3) = p_group%flow_grid%eta3_min + (j_vx-1) * p_group%flow_grid%delta_eta3
          flow_cell_eta_max(3) = p_group%flow_grid%eta3_min + (j_vx)   * p_group%flow_grid%delta_eta3
          flow_cell_eta_mid(3) = flow_cell_eta_min(3) + 0.5*p_group%flow_grid%delta_eta3

          do j_vy = 1, p_group%flow_grid%num_cells4
            flow_cell_eta_min(4) = p_group%flow_grid%eta4_min + (j_vy-1) * p_group%flow_grid%delta_eta4
            flow_cell_eta_max(4) = p_group%flow_grid%eta4_min + (j_vy)   * p_group%flow_grid%delta_eta4
            flow_cell_eta_mid(4) = flow_cell_eta_min(4) + 0.5*p_group%flow_grid%delta_eta4

            if( local_verbose )then
              print *, "[pum], A. -- Going in flow cell ", j_x, j_y, j_vx, j_vy, " ++++++++++++++++++++++++++++++++++++++++++++++++"
            end if

            do i_relevant = 1, 5

              if( local_verbose )then
                print *, "[pum], A. -- LOOKING FOR NEW RELEVANT MARKER: i_relevant = ", i_relevant
                print *, "           relevant markers already found: "
                do m = 1, i_relevant - 1
                  print *, " m = ", m
                  print *, " eta_relevant_marker(m) = ", eta_relevant_marker(m,:)
                  print *, " eta_relevant_marker(m)-eta_relevant_marker(1) = ", eta_relevant_marker(m,:)-eta_relevant_marker(1,:)
                end do
              end if

              ! start search for the relevant pointer #i_relevant, using the criteria specified above

              if( i_relevant <= 2 )then
                ! for the 1st relevant marker we look at every marker in cell, and restart once afterwards
                new_int_list_element => p_group%unstruct_markers_in_flow_cell(j_x,j_y,j_vx,j_vy)%pointed_element
              else
                ! we can discard markers that are too close from the first relevant one.
                ! Moreover we will need the projection matrix for i_relevant > 2
                aux_matrix = 0.0d0
                do l = 2, i_relevant-1
                  do m = 2, i_relevant-1
                    ! aux_matrix will be seen as a matrix of size (i_relevant-2) x (i_relevant-2)
                    aux_matrix(m-1,l-1) = p_group%anisotropic_flow_grid_scalar_product( &
                      eta_relevant_marker(l, :) - eta_relevant_marker(1, :),        &
                      eta_relevant_marker(m, :) - eta_relevant_marker(1, :)         &
                    )
                  end do
                end do
                ! compute the inverse matrix
                matrix_size = i_relevant-2
                call get_inverse_matrix_with_given_size(matrix_size, aux_matrix, projection_matrix, ok_flag)
                SLL_ASSERT( ok_flag )
              end if

              ! loop over the markers in cell to search for the relevant marker #i_relevant
              ! (do not start over the loop, as markers already seen can be discarded: either relevant already, or too close)
              relevant_marker_found = .false.
              do while( .not. relevant_marker_found )

                if( associated(new_int_list_element) )then
                  ! there is a marker in cell: take it
                  i_marker = new_int_list_element%value
                  eta_marker = p_group%unstruct_markers_eta(i_marker, :)
                  new_int_list_element => new_int_list_element%next

                  if( local_verbose )then
                    print *, "[pum], A. -- trying marker in cell,   i_marker = ", i_marker, "-------------"
                    print *, " eta               = ", eta_marker
                    print *, " flow_cell_eta_min = ", flow_cell_eta_min
                    print *, " flow_cell_eta_mid = ", flow_cell_eta_mid
                    print *, " flow_cell_eta_max = ", flow_cell_eta_max
                  end if


                else
                  ! no more markers in cell: create a new one with quasi-random coordinates
                  call i8_sobol(int(4,8), sobol_seed, rdn) ! 4 Sobol numbers in [0,1]

                  ! Transform rdn to the proper intervals
                  do i_dim = 1, 4
                    eta_marker(i_dim) =   flow_cell_eta_min(i_dim) * (1 - rdn(i_dim)) &
                                        + flow_cell_eta_max(i_dim) * rdn(i_dim)
                  end do
                  i_marker = 0  ! local trick to mean that this marker has no global index, hence is not stored in the global array
                  nb_markers_created = nb_markers_created + 1
                  if( local_verbose )then
                    print *, "[pum], A. -- created new marker ---------------------------------------------------------------------"
                    print *, " eta               = ", eta_marker
                    print *, " flow_cell_eta_min = ", flow_cell_eta_min
                    print *, " flow_cell_eta_mid = ", flow_cell_eta_mid
                    print *, " flow_cell_eta_max = ", flow_cell_eta_max
                  end if
                end if

                ! compute the projection projected_marker on the affine space of previous relevant markers
                if( i_relevant == 1 )then
                  ! for the first relevant marker we measure its distance to the center (midpoint) of the cell
                  eta_projected = flow_cell_eta_mid
                else if( i_relevant == 2 )then
                  ! trivial projection for the second relevant marker
                  eta_projected = eta_relevant_marker(1, :)
                else
                  ! getting the projected marker with the projection matrix
                  do l = 2, i_relevant-1
                    aux_b(l-1) = p_group%anisotropic_flow_grid_scalar_product(    &
                      eta_marker                - eta_relevant_marker(1, :),    &
                      eta_relevant_marker(l, :) - eta_relevant_marker(1, :)     &
                    )
                  end do
                  ! compute projection coefs, ie the c_{l-1} such that Px = x^1 + sum_{l=2}^{i-1} c_{l-1} (x^l - x^1)
                  ! (see subroutine description above)
                  eta_projected = eta_relevant_marker(1, :)
                  do l = 2, i_relevant-1
                    projection_coef(l-1) = 0.0d0
                    do m = 2, i_relevant-1
                      projection_coef(l-1) = projection_coef(l-1) + projection_matrix(l-1,m-1) * aux_b(m-1)
                    end do
                    eta_projected = eta_projected + projection_coef(l-1) * (eta_relevant_marker(l, :) - eta_relevant_marker(1, :))
                  end do
                end if

                ! measure the distance between new marker and its projection

                a_distance = p_group%anisotropic_flow_grid_distance(eta_marker, eta_projected)
                if( local_verbose )then
                  print *, "[pum], A. -- i_relevant, p_group%flow_grid_h = ", i_relevant, p_group%flow_grid_h
                  print *, "[pum], A. -- p_group%flow_grid%delta_eta1 = ", p_group%flow_grid%delta_eta1
                  print *, "[pum], A. -- p_group%flow_grid%delta_eta2 = ", p_group%flow_grid%delta_eta2
                  print *, "[pum], A. -- p_group%flow_grid%delta_eta3 = ", p_group%flow_grid%delta_eta3
                  print *, "[pum], A. -- p_group%flow_grid%delta_eta4 = ", p_group%flow_grid%delta_eta4
                  print *, "           eta_projected     = ", eta_projected
                  print *, "[pum], A. ++ distance(eta, eta_projected)    = ", a_distance
                end if

                if( a_distance > p_group%flow_grid_h / 4 )then
                  ! then this marker is relevant, store its coordinates to continue the search
                  eta_relevant_marker(i_relevant,:) = eta_marker
                  relevant_marker_found = .true.
                  if( local_verbose ) print *, "[pum], A. -- found relevant marker !"
                  if( i_marker > 0 )then
                    ! this marker is already in the cell, and we know its index
                    p_group%unstruct_markers_relevant_neighbor(i_marker) = i_marker  ! means: relevant but neighbor unknown
                    if( local_verbose ) print *, "[pum], A. -- marker already present: i_marker = ", i_marker

                  else
                    if( local_verbose ) print *, "[pum], A. -- (newly created marker)"
                    ! it is a new marker: we must add its coordinates in the markers list
                    ! -> store cell index and marker coordinates in a new element of the temporary insertion list
                    SLL_ALLOCATE( new_marker_list_element, ierr )
                    new_marker_list_element%cell_j_x = j_x
                    new_marker_list_element%cell_j_y = j_y
                    new_marker_list_element%cell_j_vx = j_vx
                    new_marker_list_element%cell_j_vy = j_vy
                    new_marker_list_element%eta = eta_marker
                    new_marker_list_element%flag = 1      ! means: relevant
                    ! increment the proper linked list
                    head_marker_list => list_of_markers_to_be_added
                    list_of_markers_to_be_added => add_element_in_marker_list(head_marker_list, new_marker_list_element)
                  end if
                else
                  ! if the test is failed then this marker is not relevant, just ignore it. Its index may be stored later for removal
                  if( local_verbose ) print *, "[pum], A. --  relevant marker not found. try again..."
                end if

              end do   ! looping over the (remaining) markers in cell to find the relevant marker #i_relevant

            end do   ! loop over the index i_relevant

            ! now that we have found 5 relevant markers we should remove (ie, mark for removal) or insert (ie, mark for insertion)
            ! some unrelevant markers in the cell.
            ! Specifically, we should have M-5 unrelevant markers per cell, where M is the prescribed (max) nb of markers per cell
            !
            ! Note: for this we can run through the list unstruct_markers_in_flow_cell(j_x,j_y,j_vx,j_vy) which contains
            ! all the unrelevant markers in the cell (some new relevant markers may be not therein yet, but we know their number)
            nb_unrelevant_markers_in_this_cell = 0
            new_int_list_element => p_group%unstruct_markers_in_flow_cell(j_x,j_y,j_vx,j_vy)%pointed_element
            do while( associated(new_int_list_element) )
              ! there is a marker in cell: check status
              i_marker = new_int_list_element%value
              if( p_group%unstruct_markers_relevant_neighbor(i_marker) == 0 )then
                ! then this marker is not relevant
                if( nb_unrelevant_markers_in_this_cell == p_group%nb_unstruct_markers_per_cell - 5 )then
                  ! we have already reached the nb of unrelevant markers in cell: add its index in the removal list
                  SLL_ALLOCATE( new_aux_int_list_element, ierr )
                  new_aux_int_list_element%value = i_marker
                  ! increment the proper linked list
                  head_int_list => list_of_marker_indices_to_be_discarded
                  list_of_marker_indices_to_be_discarded => add_element_in_int_list(      &
                                                              head_int_list,              &
                                                              new_aux_int_list_element)
                else
                  nb_unrelevant_markers_in_this_cell = nb_unrelevant_markers_in_this_cell + 1
                end if
              else
                SLL_ASSERT( p_group%unstruct_markers_relevant_neighbor(i_marker) == i_marker )
              end if
              new_int_list_element => new_int_list_element%next
            end do
            do while( nb_unrelevant_markers_in_this_cell < p_group%nb_unstruct_markers_per_cell - 5 )

              ! not enough markers in cell: create a new one with quasi-random coordinates
              call i8_sobol(int(4,8), sobol_seed, rdn) ! 4 Sobol numbers in [0,1]

              ! Transform rdn to the proper intervals
              do i_dim = 1, 4
                eta_marker(i_dim) =   flow_cell_eta_min(i_dim) * (1 - rdn(i_dim)) &
                                    + flow_cell_eta_max(i_dim) * rdn(i_dim)
              end do

              ! it is a new marker: we must add its coordinates in the markers list
              ! -> store cell index and marker coordinates in a new element of the temporary insertion list
              SLL_ALLOCATE( new_marker_list_element, ierr )
              new_marker_list_element%cell_j_x = j_x
              new_marker_list_element%cell_j_y = j_y
              new_marker_list_element%cell_j_vx = j_vx
              new_marker_list_element%cell_j_vy = j_vy
              new_marker_list_element%eta = eta_marker
              new_marker_list_element%flag = 0      ! means: not relevant
              ! increment the proper linked list
              head_marker_list => list_of_markers_to_be_added
              list_of_markers_to_be_added => add_element_in_marker_list(head_marker_list, new_marker_list_element)

              nb_markers_created = nb_markers_created + 1
              nb_unrelevant_markers_in_this_cell = nb_unrelevant_markers_in_this_cell + 1
            end do

            SLL_ASSERT( nb_unrelevant_markers_in_this_cell == p_group%nb_unstruct_markers_per_cell - 5 )
          end do  ! 4-fold loop over the flow cells
        end do
      end do
    end do

    ! next step, we add in the removal list every marker outside the flow grid (do that now so that they are removed first)
    new_int_list_element => p_group%unstruct_markers_outside_flow_grid
    do while( associated(new_int_list_element))
      i_marker = new_int_list_element%value

      ! the  marker # i_marker is outside the cell domain: add its index in the removal list (as done above)
      SLL_ALLOCATE( new_aux_int_list_element, ierr )
      new_aux_int_list_element%value = i_marker
      ! increment the proper linked list
      head_int_list => list_of_marker_indices_to_be_discarded
      list_of_marker_indices_to_be_discarded => add_element_in_int_list(      &
                                                  head_int_list,              &
                                                  new_aux_int_list_element)

      new_int_list_element => new_int_list_element%next
    end do

    if( local_verbose )then
      print *, "preparing the unstruct markers  A-end -------------- -------------- -------------- -------------- -------------- "
      print *, "preparing the unstruct markers,                    nb_markers_created = ", nb_markers_created
      print *, "preparing the unstruct markers,   -> nb_markers_created per flow cell = ", nb_markers_created    &
          * 1./(p_group%flow_grid%num_cells1*p_group%flow_grid%num_cells2*p_group%flow_grid%num_cells3*p_group%flow_grid%num_cells4)
      print *, "preparing the unstruct markers  A-end -------------- -------------- -------------- -------------- -------------- "
    end if

    ! now add in the main array the markers that have been stored in the temporary list (using those that were stored for removal)
    if( local_verbose )then
      print *, "preparing the unstruct markers, B. adding in main array the markers that have been stored in the temporary list... "
    end if

    new_marker_list_element => list_of_markers_to_be_added
    new_int_list_element => list_of_marker_indices_to_be_discarded
    do while( associated(new_marker_list_element) )
      ! there should be one marker index stored for removal
      SLL_ASSERT( associated(new_int_list_element) )
      i_marker = new_int_list_element%value
      ! it should not be flagged as relevant
      SLL_ASSERT( p_group%unstruct_markers_relevant_neighbor(i_marker) == 0 )
      ! store index of the flow cell containing the obsolete marker, to remove its index from cell list
      eta_marker_to_remove = p_group%unstruct_markers_eta(i_marker,:)
      call get_4d_cell_containing_point(    &
              eta_marker_to_remove,                     &
              p_group%flow_grid,                        &
              old_j_x, old_j_y, old_j_vx, old_j_vy,     &
              marker_is_outside                         &
           )
      ! store this marker coordinates in main array
      p_group%unstruct_markers_eta(i_marker,:) = new_marker_list_element%eta
      ! test whether the new (relevant) marker is in a different cell than the old (unrelevant) one
      j_x  = new_marker_list_element%cell_j_x
      j_y  = new_marker_list_element%cell_j_y
      j_vx = new_marker_list_element%cell_j_vx
      j_vy = new_marker_list_element%cell_j_vy
      call p_group%update_flow_cell_lists_with_new_marker_position(   &
              i_marker,                             &
              old_j_x, old_j_y, old_j_vx, old_j_vy, &
              j_x, j_y, j_vx, j_vy                  &
           )
      ! see whether this marker is relevant, and flag it as such using a temporary value in unstruct_markers_relevant_neighbor
      if( new_marker_list_element%flag == 1 )then
        ! the marker is relevant
        p_group%unstruct_markers_relevant_neighbor(i_marker) = i_marker
      else
        ! the marker is not relevant
        p_group%unstruct_markers_relevant_neighbor(i_marker) = 0
      end if
      ! move forward in the two temporary lists
      new_marker_list_element => new_marker_list_element%next
      new_int_list_element => new_int_list_element%next
    end do

    ! there should be no markers left in the removal list, since every cell now contains the maximum number of markers per cell
    SLL_ASSERT( .not. associated(new_int_list_element) )

    ! Final step: now that we have 5 relevant markers in each cell, we can build the local sets of relevant markers which need
    ! to be accessible from every marker.
    ! This will be done by storing the indices of the relevant markers inside the 'unstruct_markers_relevant_neighbor' array
    ! in such a way that given any marker with index i_marker, the indices of the 5 relevant markers in the flow cell are:
    ! i_relevant_marker_1 = p_group%unstruct_markers_relevant_neighbor(i_marker)
    ! i_relevant_marker_2 = p_group%unstruct_markers_relevant_neighbor(i_relevant_marker_1)
    ! i_relevant_marker_3 = p_group%unstruct_markers_relevant_neighbor(i_relevant_marker_2)
    ! i_relevant_marker_4 = p_group%unstruct_markers_relevant_neighbor(i_relevant_marker_3)
    ! i_relevant_marker_5 = p_group%unstruct_markers_relevant_neighbor(i_relevant_marker_4)
    !
    ! Note that here, i_marker may or may not be relevant itself, but it will always point to a relevant marker
    ! In particular we have
    !   i_marker == i_relevant_marker_5                                                             if i_marker was relevant,
    ! and
    ! i_relevant_marker_1 == p_group%unstruct_markers_relevant_neighbor(i_relevant_marker_5)         otherwise
    do j_x = 1, p_group%flow_grid%num_cells1
      if( local_verbose )then
        print *, "preparing the unstruct markers, C. building the loops of relevant markers in the flow cells: j_x = ", &
                      j_x, "/", p_group%flow_grid%num_cells1, "..."
      end if

      do j_y = 1, p_group%flow_grid%num_cells2
        do j_vx = 1, p_group%flow_grid%num_cells3
          do j_vy = 1, p_group%flow_grid%num_cells4

            ! first loop through the markers in the cell, search for one that is relevant:
            ! called 'first' here even if it was not first found above (order does not matter now, this is just to create a loop)
            i_first_relevant_marker = 0
            new_int_list_element => p_group%unstruct_markers_in_flow_cell(j_x,j_y,j_vx,j_vy)%pointed_element
            do while( associated(new_int_list_element) .and. i_first_relevant_marker == 0 )
              i_marker = new_int_list_element%value
              if( p_group%unstruct_markers_relevant_neighbor(i_marker) > 0 )then
                i_first_relevant_marker = i_marker
              end if
              new_int_list_element => new_int_list_element%next
            end do
            SLL_ASSERT( i_first_relevant_marker > 0 )

            ! second loop through the markers in the cell to create the loop
            i_last_relevant_marker = i_first_relevant_marker
            nb_relevant_markers_in_this_cell = 0      ! this counter just to check
            nb_unrelevant_markers_in_this_cell = 0    ! this counter just to check

            new_int_list_element => p_group%unstruct_markers_in_flow_cell(j_x,j_y,j_vx,j_vy)%pointed_element
            do while( associated(new_int_list_element))
              i_marker = new_int_list_element%value
              if( p_group%unstruct_markers_relevant_neighbor(i_marker) > 0 )then
                ! then this marker is relevant, will point to the last relevant marker seen (unless it is the first, for the loop)
                SLL_ASSERT( p_group%unstruct_markers_relevant_neighbor(i_marker) == i_marker )
                if( .not. i_marker == i_first_relevant_marker )then
                  ! will set value of p_group%unstruct_markers_relevant_neighbor(i_marker) at the end, for the loop
                  p_group%unstruct_markers_relevant_neighbor(i_marker) = i_last_relevant_marker
                  i_last_relevant_marker = i_marker
                end if
                nb_relevant_markers_in_this_cell = nb_relevant_markers_in_this_cell + 1
              else
                ! then this marker is not relevant, may point to any relevant marker
                p_group%unstruct_markers_relevant_neighbor(i_marker) = i_first_relevant_marker
                nb_unrelevant_markers_in_this_cell = nb_unrelevant_markers_in_this_cell + 1
              end if
              ! also set the remapped coordinates field (not affected by push) to store the position at remapping time
              p_group%unstruct_markers_eta_at_remapping_time(i_marker,:) = p_group%unstruct_markers_eta(i_marker,:)
              new_int_list_element => new_int_list_element%next
            end do
            ! close the loop by making the first relevant marker point to the last one just seen
            p_group%unstruct_markers_relevant_neighbor(i_first_relevant_marker) = i_last_relevant_marker

            if( local_verbose )then
              print *, "cell = ", j_x, j_y, j_vx, j_vy
              print *, "nb_relevant_markers_in_this_cell = ", nb_relevant_markers_in_this_cell
              print *, "nb_unrelevant_markers_in_this_cell = ", nb_unrelevant_markers_in_this_cell
            end if
            SLL_ASSERT( nb_relevant_markers_in_this_cell == 5 )
            SLL_ASSERT( nb_unrelevant_markers_in_this_cell + 5 == p_group%nb_unstruct_markers_per_cell )

            ! todo: add a test when computing the jacobian matrices, to check that the determinant is large enough

          end do  ! 4-fold loop over the flow cells
        end do
      end do
    end do

  end subroutine bsl_lt_pic_4d_prepare_unstruct_markers_for_flow_jacobians

  ! ------------------------------------------------------------------------------------------------------------------------
  ! This subroutine updates the list of markers in the flow cells ('unstruct_markers_in_flow_cell') given a marker
  ! (of global index 'i_marker') that was in some cell old_j and is now in a (possibly different) cell j.
  ! Note: markers outside the flow grid are listed in a special list: unstruct_markers_outside_flow_grid
  !
  subroutine update_flow_cell_lists_with_new_marker_position(p_group,   &
                                                             i_marker,  &
                                                             old_j_x, old_j_y, old_j_vx, old_j_vy, &
                                                             new_j_x, new_j_y, new_j_vx, new_j_vy)
    class(sll_bsl_lt_pic_4d_group), intent(inout) :: p_group
    sll_int32,                      intent(in)  :: i_marker
    sll_int32,                      intent(in)  :: old_j_x, old_j_y, old_j_vx, old_j_vy
    sll_int32,                      intent(in), optional  :: new_j_x, new_j_y, new_j_vx, new_j_vy
    type(int_list_element),       pointer           :: new_int_list_element, head_int_list
    sll_real64, dimension(4)    :: eta_marker
    sll_int32  :: j_x, j_y, j_vx, j_vy
    sll_int32  :: ierr
    logical    :: eta_is_outside_grid
    logical    :: old_eta_is_outside_grid
    logical    :: old_and_new_etas_are_in_same_cell
    logical    :: marker_found_and_removed

    if( .not. present(new_j_x) )then
      ! this is usually the case when the marker has just been moved and we do not know in which flow cell it is
      SLL_ASSERT( (.not. present(new_j_y)) .and. (.not. present(new_j_vx)) .and. (.not. present(new_j_vy)) )
      eta_marker = p_group%unstruct_markers_eta(i_marker, :)
      call get_4d_cell_containing_point(eta_marker, p_group%flow_grid, j_x, j_y, j_vx, j_vy, eta_is_outside_grid)
    else
      SLL_ASSERT( (present(new_j_y)) .and. (present(new_j_vx)) .and. (present(new_j_vy)) )
      j_x = new_j_x
      j_y = new_j_y
      j_vx = new_j_vx
      j_vy = new_j_vy
      eta_is_outside_grid = (    &
             j_x  < 1 .or. j_x >  p_group%flow_grid%num_cells1 &
        .or. j_y  < 1 .or. j_y >  p_group%flow_grid%num_cells2 &
        .or. j_vx < 1 .or. j_vx > p_group%flow_grid%num_cells3 &
        .or. j_vy < 1 .or. j_vy > p_group%flow_grid%num_cells4 &
      )
    end if

    old_eta_is_outside_grid = (  &
           old_j_x  < 1 .or. old_j_x >  p_group%flow_grid%num_cells1 &
      .or. old_j_y  < 1 .or. old_j_y >  p_group%flow_grid%num_cells2 &
      .or. old_j_vx < 1 .or. old_j_vx > p_group%flow_grid%num_cells3 &
      .or. old_j_vy < 1 .or. old_j_vy > p_group%flow_grid%num_cells4 &
    )

    ! see whether old and current marker positions are in different cells (the exterior domain is treated as one cell)
    old_and_new_etas_are_in_same_cell =                                                       &
        (j_x == old_j_x .and. j_y == old_j_y .and. j_vx == old_j_vx .and. j_vy == old_j_vy)   &
        .or.                                                                                  &
        (eta_is_outside_grid .and. old_eta_is_outside_grid)

    if( .not. old_and_new_etas_are_in_same_cell )then
      ! Then we must change the corresponding lists:
      ! 1. add the element i_marker in the list of the new cell
      SLL_ALLOCATE( new_int_list_element, ierr )
      new_int_list_element%value = i_marker
      if( eta_is_outside_grid )then
        head_int_list => p_group%unstruct_markers_outside_flow_grid
        p_group%unstruct_markers_outside_flow_grid &
                        => add_element_in_int_list(head_int_list, new_int_list_element)
      else
        head_int_list => p_group%unstruct_markers_in_flow_cell(j_x, j_y, j_vx, j_vy)%pointed_element
        p_group%unstruct_markers_in_flow_cell(j_x,j_y,j_vx,j_vy)%pointed_element &
                        => add_element_in_int_list(head_int_list, new_int_list_element)
      end if
      ! 2. remove the element i_marker from the list of the old cell
      marker_found_and_removed = .false.
      if( old_eta_is_outside_grid )then
        head_int_list => p_group%unstruct_markers_outside_flow_grid
      else
        head_int_list => p_group%unstruct_markers_in_flow_cell(old_j_x, old_j_y, old_j_vx, old_j_vy)%pointed_element
      end if
      ! test the first element in the list (requires separate treatment)
      if( head_int_list%value == i_marker )then
        ! first element is i_marker, then the cell list should point to the second one
        if( old_eta_is_outside_grid )then
          p_group%unstruct_markers_outside_flow_grid &
                => head_int_list%next
        else
          p_group%unstruct_markers_in_flow_cell(old_j_x, old_j_y, old_j_vx, old_j_vy)%pointed_element &
                => head_int_list%next
        end if
        marker_found_and_removed = .true.
      end if
      ! test the other elements in the list if needed
!      print*, "DEBUG 654654654 --                 looking for i_marker = ", i_marker
      do while( associated(head_int_list%next) .and. (.not. marker_found_and_removed) )
!        print*, "associated(head_int_list%next) = ", associated(head_int_list%next)
!        print*, "head_int_list%next%value = ", head_int_list%next%value
        if( head_int_list%next%value == i_marker )then
          ! remove the next element in list
          head_int_list%next => head_int_list%next%next
!          print*, "element FOUND and REMOVED. Now,"
!          print*, " associated(head_int_list%next) = ", associated(head_int_list%next)
!          if( associated(head_int_list%next) ) print*, "head_int_list%next%value = ", head_int_list%next%value
          marker_found_and_removed = .true.
        else
!          print*, "moving forward..."
!          if( associated(head_int_list%next) ) print*, " [BEFORE MOVING FWD] head_int_list%next%value = ", head_int_list%next%value
          head_int_list => head_int_list%next
        end if
      end do
      SLL_ASSERT( marker_found_and_removed )
      ! 3. (optional) check it is well removed:
      if( old_eta_is_outside_grid )then
        head_int_list => p_group%unstruct_markers_outside_flow_grid
      else
        head_int_list => p_group%unstruct_markers_in_flow_cell(old_j_x, old_j_y, old_j_vx, old_j_vy)%pointed_element
      end if
      do while( associated(head_int_list) )
        SLL_ASSERT( head_int_list%value .ne. i_marker )
        head_int_list => head_int_list%next
      end do
    end if

  end subroutine update_flow_cell_lists_with_new_marker_position

  ! ------------------------------------------------------------------------------------------------------------------------
  ! This subroutine computes
  !   - the current coordinates of the marker closest from cell center: x_k, y_k, vx_k, vy_k
  !   - the past (last remapping time) coordinates of this marker : x_k_t0, y_k_t0, vx_k_t0, vy_k_t0
  !   - the coefficients of the deformation matrix (backward Jacobian) : d11, d12, .. d44
  !
  ! To do so we proceed as follows:
  !   1.  we find the marker that is closest from the cell center, and its group of relevant neighbors with positions x^1, .. x^5
  !       (they have been prepared in the subroutine bsl_lt_pic_4d_prepare_unstruct_markers_for_flow_jacobians)
  !   2.  then we define the forward Jacobian matrix: J = (J_{l,m})_{l,m = 1, .. 4} by the relations
  !           J (hat x^{r+1} - hat x^1) = x^{r+1} - x^1     for   r = 1, .. 4
  !       where x^r is the position of the r-th relevant marker at the current time
  !       and hat x^r is its position at the last remapping time
  !       In particular, if we set
  !           hat_delta_eta(m,r) = (hat x^{r+1} - hat x^1)_m
  !       and
  !           delta_eta(m,r)     = (x^{r+1} - x^1)_m
  !       for  m,r = 1, .. 4, then we have
  !           J * hat_delta_eta = delta_eta
  !       which yields D = J^{-1} = hat_delta_eta * (delta_eta)^{-1}

  subroutine get_deformation_matrix_from_unstruct_markers_in_cell(  &
                      p_group,                                      &
                      j_x, j_y, j_vx, j_vy,                         &
                      mesh_period_x, mesh_period_y,                 &
                      x_k, y_k, vx_k, vy_k,                         &
                      x_k_t0, y_k_t0, vx_k_t0, vy_k_t0,             &
                      d11, d12, d13, d14,                           &
                      d21, d22, d23, d24,                           &
                      d31, d32, d33, d34,                           &
                      d41, d42, d43, d44                            &
                   )

    class(sll_bsl_lt_pic_4d_group),intent(inout) :: p_group
    sll_int32,  intent(in)  :: j_x, j_y, j_vx, j_vy                 !< indices of the flow cell
    sll_real64, intent(in)  :: mesh_period_x
    sll_real64, intent(in)  :: mesh_period_y

    sll_real64, intent(out) :: x_k, y_k, vx_k, vy_k                 !< marker coordinates at current time
    sll_real64, intent(out) :: x_k_t0, y_k_t0, vx_k_t0, vy_k_t0     !< marker coordinates in the past (last remapping time)
    sll_real64, intent(out) :: d11,d12,d13,d14  ! coefs of deformation matrix D (backward Jacobian)
    sll_real64, intent(out) :: d21,d22,d23,d24
    sll_real64, intent(out) :: d31,d32,d33,d34
    sll_real64, intent(out) :: d41,d42,d43,d44

    type(int_list_element),  pointer  :: head_int_list

    sll_real64, dimension(4)    :: flow_cell_eta_min
    sll_real64, dimension(4)    :: flow_cell_eta_mid
    sll_real64, dimension(4)    :: flow_cell_eta_max
    sll_real64, dimension(4)    :: eta_marker
    !    sll_real64, dimension(4)    :: eta_projected
    sll_real64, dimension(4,4)  :: hat_delta_eta
    sll_real64, dimension(4,4)  :: delta_eta
    sll_real64, dimension(4,4)  :: inverse_delta_eta
    sll_real64, dimension(4,4)  :: bwd_jacobian

    logical    :: ok_flag
    logical    :: domain_is_x_periodic
    logical    :: domain_is_y_periodic

    sll_int32  :: cell_step
    sll_int32  :: nb_markers_in_cell
    sll_int32  :: j_aux_x, j_aux_y, j_aux_vx, j_aux_vy
    sll_int32  :: i_closest_marker
    sll_int32  :: i_marker
    sll_int32  :: i_first_relevant_marker
    sll_int32  :: i_next_relevant_marker
    sll_int32  :: matrix_size
    sll_int32  :: l, m, r
    sll_real64 :: closest_marker_distance
    sll_real64 :: this_marker_distance

    domain_is_x_periodic = p_group%domain_is_periodic(1)
    domain_is_y_periodic = p_group%domain_is_periodic(2)

    SLL_ASSERT( p_group%flow_markers_type == SLL_BSL_LT_PIC_UNSTRUCTURED )

    flow_cell_eta_min(1) = p_group%flow_grid%eta1_min + (j_x-1) * p_group%flow_grid%delta_eta1
    flow_cell_eta_max(1) = p_group%flow_grid%eta1_min + (j_x)   * p_group%flow_grid%delta_eta1
    flow_cell_eta_mid(1) = flow_cell_eta_min(1) + 0.5*p_group%flow_grid%delta_eta1

    flow_cell_eta_min(2) = p_group%flow_grid%eta2_min + (j_y-1) * p_group%flow_grid%delta_eta2
    flow_cell_eta_max(2) = p_group%flow_grid%eta2_min + (j_y)   * p_group%flow_grid%delta_eta2
    flow_cell_eta_mid(2) = flow_cell_eta_min(2) + 0.5*p_group%flow_grid%delta_eta2

    flow_cell_eta_min(3) = p_group%flow_grid%eta3_min + (j_vx-1) * p_group%flow_grid%delta_eta3
    flow_cell_eta_max(3) = p_group%flow_grid%eta3_min + (j_vx)   * p_group%flow_grid%delta_eta3
    flow_cell_eta_mid(3) = flow_cell_eta_min(3) + 0.5*p_group%flow_grid%delta_eta3

    flow_cell_eta_min(4) = p_group%flow_grid%eta4_min + (j_vy-1) * p_group%flow_grid%delta_eta4
    flow_cell_eta_max(4) = p_group%flow_grid%eta4_min + (j_vy)   * p_group%flow_grid%delta_eta4
    flow_cell_eta_mid(4) = flow_cell_eta_min(4) + 0.5*p_group%flow_grid%delta_eta4

    ! 1-a. loop through the markers in cell to find the one closest to the center
    i_closest_marker = 0
    closest_marker_distance = 1d30
    SLL_ASSERT( j_x  >= 1 .and. j_x  <= p_group%flow_grid%num_cells1 )
    SLL_ASSERT( j_y  >= 1 .and. j_y  <= p_group%flow_grid%num_cells2 )
    SLL_ASSERT( j_vx >= 1 .and. j_vx <= p_group%flow_grid%num_cells3 )
    SLL_ASSERT( j_vy >= 1 .and. j_vy <= p_group%flow_grid%num_cells4 )
    head_int_list => p_group%unstruct_markers_in_flow_cell(j_x,j_y,j_vx,j_vy)%pointed_element
    nb_markers_in_cell = 0
    do while( associated(head_int_list) )
      i_marker = head_int_list%value
      eta_marker = p_group%unstruct_markers_eta(i_marker, :)
      this_marker_distance = p_group%anisotropic_flow_grid_distance(eta_marker, flow_cell_eta_mid)
      if( this_marker_distance < closest_marker_distance )then
        i_closest_marker = i_marker
        closest_marker_distance = this_marker_distance
      end if
      nb_markers_in_cell = nb_markers_in_cell + 1
      head_int_list => head_int_list%next
    end do

    ! if the given cell had no marker (this is always possible), we take any marker in a neighboring cell
    if( i_closest_marker == 0 )then
      SLL_ASSERT( nb_markers_in_cell == 0 )
      cell_step = 0
      do while( i_closest_marker == 0 )
        cell_step = cell_step + 1
        do j_aux_x = j_x - cell_step, j_x + cell_step
          if( i_closest_marker > 0 ) EXIT
          if( j_aux_x < 1 .or. j_aux_x > p_group%flow_grid%num_cells1 ) CYCLE

          do j_aux_y = j_y - cell_step, j_y + cell_step
            if( i_closest_marker > 0 ) EXIT
            if( j_aux_y < 1 .or. j_aux_y > p_group%flow_grid%num_cells2 ) CYCLE

            do j_aux_vx = j_vx - cell_step, j_vx + cell_step
              if( i_closest_marker > 0 ) EXIT
              if( j_aux_vx < 1 .or. j_aux_vx > p_group%flow_grid%num_cells3 ) CYCLE

              do j_aux_vy = j_vy - cell_step, j_vy + cell_step
                if( i_closest_marker > 0 ) EXIT
                if( j_aux_vy < 1 .or. j_aux_vy > p_group%flow_grid%num_cells4 ) CYCLE

                ! take the first marker in the first cell at distance cell_step from given cell (here l_1 but could be l_inf)
                if( abs(j_aux_x - j_x) + abs(j_aux_y - j_y) + abs(j_aux_vx - j_vx) + abs(j_aux_vy - j_vy) == cell_step )then
                  head_int_list => p_group%unstruct_markers_in_flow_cell(j_aux_x,j_aux_y,j_aux_vx,j_aux_vy)%pointed_element
                  if( associated(head_int_list) )then
                    ! then there is a marker in that cell, just take it
                    i_closest_marker = head_int_list%value
                  end if
                end if
              end do
            end do
          end do
        end do
      end do
    end if

    SLL_ASSERT( i_closest_marker > 0 )

    ! 1-b. get the 5 relevant markers close to this marker, and store their positions at current and remapping time
    ! reminder: hat_delta_eta(m,r) = (hat x^{r+1} - hat x^1)_m

    ! (i_closest_marker may not be relevant, but its relevant_neighbor is always one)
    i_first_relevant_marker = p_group%unstruct_markers_relevant_neighbor(i_closest_marker)
    i_next_relevant_marker = i_first_relevant_marker
    do r = 1, 4
      i_next_relevant_marker = p_group%unstruct_markers_relevant_neighbor(i_next_relevant_marker)
      hat_delta_eta(:,r) =  p_group%unstruct_markers_eta_at_remapping_time(i_next_relevant_marker, :)  &
                          - p_group%unstruct_markers_eta_at_remapping_time(i_first_relevant_marker,:)
      delta_eta(:,r)     =  p_group%unstruct_markers_eta(i_next_relevant_marker, :)  &
                          - p_group%unstruct_markers_eta(i_first_relevant_marker,:)
      ! rectify the values in periodic case: difference should not be too large
      if( domain_is_x_periodic )then
        if( delta_eta(1,r) > +0.5*mesh_period_x ) delta_eta(1,r) = delta_eta(1,r) - mesh_period_x
        if( delta_eta(1,r) < -0.5*mesh_period_x ) delta_eta(1,r) = delta_eta(1,r) + mesh_period_x
      end if
      if( domain_is_y_periodic )then
        if( delta_eta(2,r) > +0.5*mesh_period_y ) delta_eta(2,r) = delta_eta(2,r) - mesh_period_y
        if( delta_eta(2,r) < -0.5*mesh_period_y ) delta_eta(2,r) = delta_eta(2,r) + mesh_period_y
      end if
      ! note: no need to rectify hat_delta_eta which corresponds to relative positions at remapping time, within flow cells
    end do

    ! check that the relevant markers form a loop of length 5
    SLL_ASSERT( p_group%unstruct_markers_relevant_neighbor(i_next_relevant_marker) == i_first_relevant_marker )

    ! 2. compute D = J^{-1} = hat_delta_eta * (delta_eta)^{-1}
    matrix_size = 4
    call get_inverse_matrix_with_given_size( matrix_size, delta_eta, inverse_delta_eta, ok_flag )
    SLL_ASSERT( ok_flag )
    bwd_jacobian = 0.0d0
    do l = 1, 4
      do m = 1, 4
        do r = 1, 4
          bwd_jacobian(l,m) =  hat_delta_eta(l,r) * inverse_delta_eta(r,m)
        end do
      end do
    end do

    ! then store the results in output variables
    x_k  = p_group%unstruct_markers_eta(i_closest_marker,1)
    y_k  = p_group%unstruct_markers_eta(i_closest_marker,2)
    vx_k = p_group%unstruct_markers_eta(i_closest_marker,3)
    vy_k = p_group%unstruct_markers_eta(i_closest_marker,4)

    x_k_t0  = p_group%unstruct_markers_eta_at_remapping_time(i_closest_marker,1)
    y_k_t0  = p_group%unstruct_markers_eta_at_remapping_time(i_closest_marker,2)
    vx_k_t0 = p_group%unstruct_markers_eta_at_remapping_time(i_closest_marker,3)
    vy_k_t0 = p_group%unstruct_markers_eta_at_remapping_time(i_closest_marker,4)

    d11 = bwd_jacobian(1,1)
    d12 = bwd_jacobian(1,2)
    d13 = bwd_jacobian(1,3)
    d14 = bwd_jacobian(1,4)

    d21 = bwd_jacobian(2,1)
    d22 = bwd_jacobian(2,2)
    d23 = bwd_jacobian(2,3)
    d24 = bwd_jacobian(2,4)

    d31 = bwd_jacobian(3,1)
    d32 = bwd_jacobian(3,2)
    d33 = bwd_jacobian(3,3)
    d34 = bwd_jacobian(3,4)

    d41 = bwd_jacobian(4,1)
    d42 = bwd_jacobian(4,2)
    d43 = bwd_jacobian(4,3)
    d44 = bwd_jacobian(4,4)

  end subroutine get_deformation_matrix_from_unstruct_markers_in_cell

  function anisotropic_flow_grid_scalar_product(p_group, eta_a, eta_b) result(val)
    class(sll_bsl_lt_pic_4d_group), intent(inout) :: p_group
    sll_real64, dimension(4),       intent(in)    :: eta_a
    sll_real64, dimension(4),       intent(in)    :: eta_b
    sll_real64                                    :: val

    val =   ( eta_a(1) / p_group%flow_grid_a1 * eta_b(1) / p_group%flow_grid_a1 )   &
          + ( eta_a(2) / p_group%flow_grid_a2 * eta_b(2) / p_group%flow_grid_a2 )   &
          + ( eta_a(3) / p_group%flow_grid_a3 * eta_b(3) / p_group%flow_grid_a3 )   &
          + ( eta_a(4) / p_group%flow_grid_a4 * eta_b(4) / p_group%flow_grid_a4 )

  end function

  function anisotropic_flow_grid_distance(p_group, eta_a, eta_b) result(val)
    class(sll_bsl_lt_pic_4d_group), intent(inout) :: p_group
    sll_real64, dimension(4),       intent(in)    :: eta_a
    sll_real64, dimension(4),       intent(in)    :: eta_b
    sll_real64                                    :: val

    val = sqrt(p_group%anisotropic_flow_grid_scalar_product(eta_a - eta_b, eta_a - eta_b))

  end function


  !> reset the structured markers on the initial (markers) grid
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

    SLL_ASSERT( p_group%flow_markers_type == SLL_BSL_LT_PIC_STRUCTURED )

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
  subroutine bsl_lt_pic_4d_remap( self )
    class(sll_bsl_lt_pic_4d_group),intent(inout) :: self

    !> A. write the nodal values of the remapped_f (evaluated with the bsl_lt_pic method) and compute the new interpolation coefs
    print *, "bsl_lt_pic_4d_remap -- step A"
    call self%bsl_lt_pic_4d_remap_f()

    !> B. and reset the (flow) markers on the initial (cartesian) grid -- no need to reset their connectivity
    print *, "bsl_lt_pic_4d_remap -- step B"
    if( self%flow_markers_type == SLL_BSL_LT_PIC_STRUCTURED )then
      call self%bsl_lt_pic_4d_reset_markers_position()
    else
      SLL_ASSERT( self%flow_markers_type == SLL_BSL_LT_PIC_UNSTRUCTURED )
      call self%bsl_lt_pic_4d_prepare_unstruct_markers_for_flow_jacobians()
    end if

    !> C. if deposition particles are pushed, we remap them now
    if( self%deposition_particles_type == SLL_BSL_LT_PIC_FLEXIBLE     &
        .and. self%deposition_particles_move_type == SLL_BSL_LT_PIC_PUSHED )then
      print *, "bsl_lt_pic_4d_remap -- step C"
      print *, "bsl_lt_pic_4d_remap -- (C) will reset ", self%number_deposition_particles, "deposition_particles..."
      ! if deposition particles are fixed, then they are initialized at each time step, in the deposition routine
      call self%reset_deposition_particles_coordinates()
      ! since the remapping tool has been reset, computing the weights can be done with straightforward interpolation (flow = Id)
      call self%reset_deposition_particles_weights_with_direct_interpolation()
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

    if( scenario == SLL_BSL_LT_PIC_SET_WEIGHTS_ON_DEPOSITION_PARTICLES )then
        ! then the deposition particles should be of flexible and fixed type
        ! (indeed in the pushed type the weights are set after remapping with a simple interpolation, no BSL_LT_PIC reconstruction)
        SLL_ASSERT( p_group%deposition_particles_type == SLL_BSL_LT_PIC_FLEXIBLE )
        SLL_ASSERT( p_group%deposition_particles_move_type == SLL_BSL_LT_PIC_FIXED )
    end if

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
        g => p_group%deposition_particles_grid

        ! the boundary nodes of the deposition grid are inside the domain, even with periodic boundary conditions
        g_num_points_x  = g%num_cells1 + 1
        g_num_points_y  = g%num_cells2 + 1
        g_num_points_vx = g%num_cells3 + 1
        g_num_points_vy = g%num_cells4 + 1


      else

        ! the deposition particles will be created to deposit their charge but not stored in memory
        !        number_of_deposition_particles_per_flow_cell = number_deposition_particles / (  flow_grid_num_cells_x    &
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

        !  nodes_coordinate_list => p_group%sparse_grid_interpolator%hierarchy(node_index)%coordinate
        deposition_particle_charge_factor = p_group%get_deposition_particle_charge_factor()

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
          nodes_in_flow_cell(j_x,j_y,j_vx,j_vy)%pointed_element => add_element_in_int_list(head, new_int_list_element)

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

    !> B. Preparatory work for the linearization of the flow on the flow cells -- in the case of structured markers
    !>    - find out the closest marker to each cell center,
    !>      by looping over all markers and noting which flow cell contains it.
    !>      (The leftmost flow cell in each dimension may not be complete.)
    if( p_group%flow_markers_type == SLL_BSL_LT_PIC_STRUCTURED )then
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
    end if

    if( .not. ( p_group%domain_is_periodic(1) .and. p_group%domain_is_periodic(2) ) )then
      print*, "WARNING -- STOP -- verify that the non-periodic case is well implemented"
      stop
    end if

    !> C. loop on the flow cells (main loop) -- on each flow cell, we
    !>   - C.1 linearize the flow using the position of the markers
    !>   - C.2 find the relevant points where f should be reconstructed
    !>   - C.3 reconstruct f on these points (using the affine backward flow and the interpolation tool for the remapped_f)
    !>   - C.4 write the resulting f value or deposit the deposition particle just created (depending on the scenario)

    ! <<loop_on_flow_cells>> [[file:~/mcp/maltpic/ltpic-bsl.tex::algo:pic-vr:loop_over_all_cells]]

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

            if( p_group%flow_markers_type == SLL_BSL_LT_PIC_STRUCTURED )then
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

            else
              ! unstructured flow markers
              SLL_ASSERT( p_group%flow_markers_type == SLL_BSL_LT_PIC_UNSTRUCTURED )

              ! using the lists of markers and relevant markers in cells, compute:
              !   - the current coordinates of the marker closest from cell center: x_k, y_k, vx_k, vy_k
              !   - the past (last remapping time) coordinates of this marker : x_k_t0, y_k_t0, vx_k_t0, vy_k_t0
              !   - the coefficients of the deformation matrix (backward Jacobian) : d11, d12, .. d44
              call p_group%get_deformation_matrix_from_unstruct_markers_in_cell(    &
                      j_x, j_y, j_vx, j_vy,                                         &
                      mesh_period_x, mesh_period_y,                                 &
                      x_k, y_k, vx_k, vy_k,                                         &
                      x_k_t0, y_k_t0, vx_k_t0, vy_k_t0,                             &
                      d11, d12, d13, d14,                                           &
                      d21, d22, d23, d24,                                           &
                      d31, d32, d33, d34,                                           &
                      d41, d42, d43, d44                                            &
                   )

            end if

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

        k_ngb  = p_group%struct_markers_list(k)%ngb_xright_index
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


        k_ngb  = p_group%struct_markers_list(k)%ngb_xleft_index
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

        k_ngb  = p_group%struct_markers_list(k)%ngb_yright_index
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

        k_ngb  = p_group%struct_markers_list(k)%ngb_yleft_index
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

        k_ngb  = p_group%struct_markers_list(k)%ngb_vxright_index
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

        k_ngb  = p_group%struct_markers_list(k)%ngb_vxleft_index
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

        k_ngb  = p_group%struct_markers_list(k)%ngb_vyright_index
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

        k_ngb  = p_group%struct_markers_list(k)%ngb_vyleft_index
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
    SLL_DEALLOCATE(particle_group%struct_markers_list, ierr)
    SLL_DEALLOCATE(particle_group%space_mesh_2d, ierr)
    SLL_DEALLOCATE(particle_group, ierr)

  end subroutine sll_bsl_lt_pic_4d_group_delete


end module sll_m_bsl_lt_pic_4d_group
