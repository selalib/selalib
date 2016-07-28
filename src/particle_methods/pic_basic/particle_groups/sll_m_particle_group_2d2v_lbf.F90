!**************************************************************
!  Copyright INRIA
!  Authors : 
!     MCP
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

!> @author MCP

!> @brief Module for a particle group with linearized-backward-flow (lbf) resamplings

module sll_m_particle_group_2d2v_lbf

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_cartesian_meshes, only: &
    sll_f_new_cartesian_mesh_4d, &
    sll_t_cartesian_mesh_4d

  use sll_m_particle_group_base, only: &
    sll_c_particle_group_base, &
    sll_t_species

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_initial_distribution, only: &
     sll_c_distribution_params

  use sll_m_sparse_grid_4d, only: &
    sll_t_sparse_grid_interpolator_4d

  use sll_m_particle_lbf_utilities, only: &
    sll_t_int_list_element, &
    sll_t_int_list_element_ptr, &
    sll_f_add_element_in_int_list, &
    sll_s_convert_4d_index_to_1d, &
    sll_s_convert_1d_index_to_4d, &
    sll_s_update_closest_particle_arrays, &
    sll_s_get_1d_cell_containing_point

    ! sll_f_eval_hat_function, &

  implicit none

  public :: &
    sll_t_particle_group_2d2v_lbf, &
    sll_s_new_particle_group_2d2v_lbf_ptr

  private

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> possible values for the parameter reconstruction_set_type in the reconstruction routine
  sll_int32, parameter :: sll_p_lbf_remapping_grid = 1
  sll_int32, parameter :: sll_p_lbf_given_grid = 2

  !> Group of @ref sll_t_particle_group_2d2v_lbf
  type, extends(sll_c_particle_group_base) :: sll_t_particle_group_2d2v_lbf

   !> @name The particles
   !> @{
    sll_int32                                                   :: n_particles_x
    sll_int32                                                   :: n_particles_y
    sll_int32                                                   :: n_particles_vx
    sll_int32                                                   :: n_particles_vy
    sll_real64, allocatable                                     :: particle_array(:,:) !< array of particles
    sll_real64                                                  :: common_weight      !< not needed for now -> put in base class?
    !> @}

    !> @name The lbf grid
    !> @{
    !>   for simplicity, this grid is used to
    !>   - define the backward linearized flow
    !>   - and sample and resample the particles (with one particle at center of each cell)
    !>   - hence its domain is also that of the sparse grid
    type(sll_t_cartesian_mesh_4d), pointer           :: lbf_grid
    logical                                          :: domain_is_periodic(2)
    !> @}


    !> @name Sparse grid interpolation of the remapped density f, dimensions of remapping grid = particle grid
    !> @{
    sll_int32                                        :: remapped_f_interpolation_degree
    type(sll_t_sparse_grid_interpolator_4d)          :: sparse_grid_interpolator
    sll_int32,  dimension(4)                         :: sparse_grid_max_levels
    sll_real64, dimension(:), allocatable            :: remapped_f_sparse_grid_coefficients
    sll_real64, dimension(:), allocatable            :: tmp_f_values_on_remapping_sparse_grid   !< used only during remapping
    !> @}

  contains

    !> @name Getters
    !> @{
    procedure :: get_x          => get_x_2d2v_lbf         !> Get the position of a particle
    procedure :: get_v          => get_v_2d2v_lbf         !> Get the velocity of a particle
    procedure :: get_charge     => get_charge_2d2v_lbf    !> Get the charge(s) of a particle
    procedure :: get_mass       => get_mass_2d2v_lbf      !> Get the mass(es) of a particle
    procedure :: get_weights    => get_weights_2d2v_lbf   !> Get weight(s) of a particle
    !> @}
    
    !> @name Setters
    !> @{
    procedure :: set_x              => set_x_2d2v_lbf  !> Set the position of a particle
    procedure :: set_v              => set_v_2d2v_lbf  !> Set the velocity of a particle
    procedure :: set_weights        => set_weights_2d2v_lbf         !> Set the weight(s) of a particle
    procedure :: set_common_weight  => set_common_weight_2d2v_lbf   !> Set the common weight for the particles
    !> @}

    !> Initializer and destructor
    procedure :: init => initialize_particle_group_2d2v_lbf   !> Initialization function
    procedure :: free => delete_particle_group_2d2v_lbf       !> Destructor

    !> @name Sampling and resampling
    !> @{
    !    procedure :: write_hat_density_on_remapping_grid        !> this evaluates an analytic f0
    !    procedure :: write_landau_density_on_remapping_grid     !> this evaluates an analytic f0
    !
    procedure, private :: write_known_density_on_remapping_grid     !> this evaluates some known distribution
    procedure, private :: reset_particles_positions
    procedure, private :: reset_particles_weights_with_direct_interpolation
    procedure, private :: reconstruct_f_lbf
    procedure, private :: reconstruct_f_lbf_on_remapping_grid
    procedure :: reconstruct_f_lbf_on_given_grid
    procedure :: sample
    procedure :: resample
    !> @}

    procedure, private :: interpolate_value_of_remapped_f
    procedure, private :: get_ltp_deformation_matrix              !> FD approx of local bwd flow Jacobian matrix
    procedure, private :: get_neighbor_index

  end type sll_t_particle_group_2d2v_lbf


contains

  !> Getters -------------------------------------------------------------------------------------------------------------------

  !----------------------------------------------------------------------------------------------------------------------------
  !> Get the physical coordinates of a particle
  pure function get_x_2d2v_lbf( self, i ) result( r )
    class( sll_t_particle_group_2d2v_lbf ), intent( in ) :: self  !< particle group
    sll_int32                             , intent( in ) :: i     !< particle index
    sll_real64 :: r(3)  !< first two components hold the value of the particle position

    r = 1.0_f64
    r(1:2) = self%particle_array(1:2, i)

  end function get_x_2d2v_lbf

  !----------------------------------------------------------------------------------------------------------------------------
  !> Get the velocity of a particle
  pure function get_v_2d2v_lbf( self, i ) result( r )
    class( sll_t_particle_group_2d2v_lbf ), intent( in ) :: self  !< particle group
    sll_int32                             , intent( in ) :: i     !< particle index
    sll_real64 :: r(3)  !< first two components hold the value of the particle velocity

    r = 1.0_f64
    r(1:2) = self%particle_array(3:4, i)

  end function get_v_2d2v_lbf



  !----------------------------------------------------------------------!
  !> Get charge of a particle (q * particle_weight)
  pure function get_charge_2d2v_lbf( self, i , i_weight) result (r)
        class( sll_t_particle_group_2d2v_lbf ), intent( in ) :: self      !< particle group
    sll_int32                                 , intent( in ) :: i         !< particle index
    sll_int32, optional                       , intent( in ) :: i_weight  !< weight index to be used (default: 1) - not needed now
    sll_real64 :: r     !< charge of particle i

    sll_int32 :: i_wi

    i_wi = 1
    if(present(i_weight)) i_wi = i_weight   !< particles in this class have only one weight so i_weight is not needed
    r = self%species%q  * self%particle_array(4+i_wi, i) * self%common_weight

  end function get_charge_2d2v_lbf


  !----------------------------------------------------------------------!
  !> Get mass of a particle ( m * particle_weight)
  pure function get_mass_2d2v_lbf( self, i, i_weight) result (r)
        class( sll_t_particle_group_2d2v_lbf ), intent( in ) :: self      !< particle group
    sll_int32                                 , intent( in ) :: i         !< particle index
    sll_int32, optional                       , intent( in ) :: i_weight  !< weight index to be used (default: 1) - not needed now
    sll_real64 :: r !< mass of particle i

    sll_int32 :: i_wi

    i_wi = 1
    if(present(i_weight)) i_wi = i_weight
    r = self%species%m  * self%particle_array(4+i_wi, i) * self%common_weight

  end function get_mass_2d2v_lbf


  !----------------------------------------------------------------------!
  !> Get weights of a particle
  pure function get_weights_2d2v_lbf( self, i) result (r)
    class( sll_t_particle_group_2d2v_lbf ), intent( in ) :: self  !< particle group
    sll_int32                             , intent( in ) :: i !< no. of the particle
    sll_real64 :: r(self%n_weights) !< particle weight(s)

    r = self%particle_array(5:4+self%n_weights, i)

  end function get_weights_2d2v_lbf

  !> Setters -------------------------------------------------------------------------------------------------------------------

  !----------------------------------------------------------------------!
  !> Set the position of a particle
  subroutine set_x_2d2v_lbf( self, i, x )
    class( sll_t_particle_group_2d2v_lbf ), intent( inout ) :: self   !< particle group
    sll_int32                         , intent( in )      :: i      !< particle index
    sll_real64                        , intent( in )      :: x(3)   !< components 1 and 2 hold the particle position to be set

    self%particle_array(1:2, i) = x(1:2)

  end subroutine set_x_2d2v_lbf

  !----------------------------------------------------------------------!
  !> Set the velocity of a particle
  subroutine set_v_2d2v_lbf( self, i, x )
    class( sll_t_particle_group_2d2v_lbf ), intent( inout ) :: self   !< particle group
    sll_int32                             , intent( in )    :: i      !< particle index
    sll_real64                            , intent( in )    :: x(3)   !< component 1 and 2 hold the particle velocity to be set

    self%particle_array(3:4, i) = x(1:2)

  end subroutine set_v_2d2v_lbf

  !----------------------------------------------------------------------!
  !> Set the weights of a particle
  subroutine set_weights_2d2v_lbf( self, i, x )
    class( sll_t_particle_group_2d2v_lbf ), intent( inout ) :: self   !< particle group
    sll_int32                             , intent( in )      :: i      !< particle index
    sll_real64                            , intent( in )      :: x(self%n_weights) !< particle weight(s) to be set

    self%particle_array(5:4+self%n_weights, i) = x

  end subroutine set_weights_2d2v_lbf

  !----------------------------------------------------------------------!
  !> Set the common weight
  subroutine set_common_weight_2d2v_lbf( self, x )
    class( sll_t_particle_group_2d2v_lbf ), intent( inout ) :: self  !< particle group
    sll_real64                            ,  intent( in )   :: x !< common weight

    self%common_weight = x

  end subroutine set_common_weight_2d2v_lbf


  !> Initializer ---------------------------------------------------------------------------------------------------------------
  !> Initialize particle group
  subroutine initialize_particle_group_2d2v_lbf  (  &
    self,            &
    species_charge,    &
    species_mass,      &
    domain_is_x_periodic,    &
    domain_is_y_periodic,    &
    remap_degree,    &
    remapping_grid_eta_min, &
    remapping_grid_eta_max, &
    remapping_sparse_grid_max_levels, &
    n_particles_x,  &
    n_particles_y,  &
    n_particles_vx, &
    n_particles_vy, &
    n_weights   &
  )

    class( sll_t_particle_group_2d2v_lbf ), intent( inout ) :: self  !< particle group

    sll_real64,               intent(in)  :: species_charge
    sll_real64,               intent(in)  :: species_mass
    logical,                  intent(in)  :: domain_is_x_periodic
    logical,                  intent(in)  :: domain_is_y_periodic
    sll_int32,                intent(in)  :: remap_degree
    sll_real64, dimension(4), intent(in)  :: remapping_grid_eta_min
    sll_real64, dimension(4), intent(in)  :: remapping_grid_eta_max
    sll_int32,  dimension(4), intent(in)  :: remapping_sparse_grid_max_levels
    sll_int32,                intent(in)  :: n_particles_x
    sll_int32,                intent(in)  :: n_particles_y
    sll_int32,                intent(in)  :: n_particles_vx
    sll_int32,                intent(in)  :: n_particles_vy
    sll_int32,                intent(in)  :: n_weights      !< number of weights per particle (only 1 for now)

    character(len=*), parameter :: this_fun_name = "initialize_particle_group_2d2v_lbf"
    sll_int32               :: ierr

    print*, "[", this_fun_name, "] - initializing the particle group... "

    self%domain_is_periodic(1) = domain_is_x_periodic
    self%domain_is_periodic(2) = domain_is_y_periodic

    self%n_particles_x = n_particles_x
    self%n_particles_y = n_particles_y
    self%n_particles_vx = n_particles_vx
    self%n_particles_vy = n_particles_vy
    self%n_particles = n_particles_x * n_particles_y * n_particles_vx * n_particles_vy
    self%n_total_particles = self%n_particles   !< sequential runs for now
    self%n_weights = n_weights    !< one weight per particle for now

    !> create the species object for this particle group
    allocate(self%species, stat=ierr)
    SLL_ASSERT( ierr == 0)
    call self%species%init( species_charge, species_mass )

    SLL_ALLOCATE( self%particle_array(4+n_weights, self%n_particles), ierr)
    SLL_ASSERT( ierr == 0)

    !> A. discretization of the flow:
    !>    the flow grid is the 4d cartesian grid where the backward flow is linearized
    !>    in this deterministic case, the particles are initially located at the center of each flow cell
    self%lbf_grid => sll_f_new_cartesian_mesh_4d( &
      n_particles_x, &
      n_particles_y, &
      n_particles_vx, &
      n_particles_vy, &
      remapping_grid_eta_min(1), &
      remapping_grid_eta_max(1), &
      remapping_grid_eta_min(2), &
      remapping_grid_eta_max(2), &
      remapping_grid_eta_min(3), &
      remapping_grid_eta_max(3), &
      remapping_grid_eta_min(4), &
      remapping_grid_eta_max(4) )

    !> B. discretization of the remapped f, with sparse grids
    self%remapped_f_interpolation_degree = remap_degree

    print*, "[", this_fun_name, "] - sparse grid levels for the remapping tool:", remapping_sparse_grid_max_levels
    self%sparse_grid_max_levels = remapping_sparse_grid_max_levels
    call self%sparse_grid_interpolator%initialize( &
                  self%sparse_grid_max_levels,   &
                  self%remapped_f_interpolation_degree,   &
                  self%remapped_f_interpolation_degree+1,    &
                  0,    &                                     !< interpolation_type for the sparse grid (splines or Lagrange)
                  remapping_grid_eta_min,    &
                  remapping_grid_eta_max    &
                  )

    !> C.  array of sparse grid coefficients for remapped_f
    SLL_ALLOCATE( self%remapped_f_sparse_grid_coefficients(self%sparse_grid_interpolator%size_basis), ierr )
    self%remapped_f_sparse_grid_coefficients = 0.0_f64

    !> this array is to store temporary values of remapped f, while those of previously remapped f still needed
    SLL_ALLOCATE( self%tmp_f_values_on_remapping_sparse_grid(self%sparse_grid_interpolator%size_basis), ierr)

  end subroutine initialize_particle_group_2d2v_lbf

  !-----------------------------------------------------------------------------------------------------------------------------
  !> Constructor for abstract type
  subroutine sll_s_new_particle_group_2d2v_lbf_ptr(&
      particle_group, &
      species_charge,    &
      species_mass,      &
      domain_is_x_periodic,    &
      domain_is_y_periodic,    &
      remap_degree,    &
      remapping_grid_eta_min, &
      remapping_grid_eta_max, &
      remapping_sparse_grid_max_levels, &
      n_particles_x,  &
      n_particles_y,  &
      n_particles_vx, &
      n_particles_vy &
    )
    class( sll_c_particle_group_base ),  pointer, intent( out )  :: particle_group
    sll_real64,               intent(in)  :: species_charge
    sll_real64,               intent(in)  :: species_mass
    logical,                  intent(in)  :: domain_is_x_periodic
    logical,                  intent(in)  :: domain_is_y_periodic
    sll_int32,                intent(in)  :: remap_degree
    sll_real64, dimension(4), intent(in)  :: remapping_grid_eta_min
    sll_real64, dimension(4), intent(in)  :: remapping_grid_eta_max
    sll_int32,  dimension(4), intent(in)  :: remapping_sparse_grid_max_levels
    sll_int32,                intent(in)  :: n_particles_x
    sll_int32,                intent(in)  :: n_particles_y
    sll_int32,                intent(in)  :: n_particles_vx
    sll_int32,                intent(in)  :: n_particles_vy

    sll_int32 :: ierr
    sll_int32 :: n_weights

    SLL_ALLOCATE( sll_t_particle_group_2d2v_lbf :: particle_group, ierr)
    n_weights = 1

    select type( particle_group )
    type is ( sll_t_particle_group_2d2v_lbf )
      call particle_group%init( &
        species_charge,    &
        species_mass,      &
        domain_is_x_periodic,    &
        domain_is_y_periodic,    &
        remap_degree,    &
        remapping_grid_eta_min, &
        remapping_grid_eta_max, &
        remapping_sparse_grid_max_levels, &
        n_particles_x,  &
        n_particles_y,  &
        n_particles_vx, &
        n_particles_vy, &
        n_weights &
      )
    end select

  end subroutine sll_s_new_particle_group_2d2v_lbf_ptr

  !> Destructor ----------------------------------------------------------------------------------------------------------------
  subroutine delete_particle_group_2d2v_lbf( self )
    class( sll_t_particle_group_2d2v_lbf ), intent( inout ) :: self  !< particle group

    deallocate(self%particle_array)
    deallocate(self%remapped_f_sparse_grid_coefficients)
    deallocate(self%tmp_f_values_on_remapping_sparse_grid)

  end subroutine

  !> ---------------------------------------------------------------------------------------------------------------------------
  !> Below this line are functions specific to the lbf particles with resampling -----------------------------------------------



  !> sample (layer for resample procedure) -------------------------------------------------------------------------------------
  !> this routine essentially does 2 things:
  !>  * locates the particles on an initial cartesian grid,
  !>  * computes interpolation (sparse grid) coefficients to approximate the initial distribution
  subroutine sample( self, target_total_charge, enforce_total_charge, init_f_params, rand_seed, rank, world_size )
    class(sll_t_particle_group_2d2v_lbf),   intent( inout ) :: self
    sll_real64,                       intent( in )    :: target_total_charge
    logical,                          intent( in )    :: enforce_total_charge
    class(sll_c_distribution_params), intent( in )    :: init_f_params               !< parameters of the initial distribution
    sll_int32, dimension(:)         , intent( in ), optional :: rand_seed
    sll_int32                       , intent( in ), optional :: rank, world_size

    if( present(rand_seed) )then
      SLL_ASSERT( present( rank ) )
      SLL_ASSERT( present( world_size ) )
      call self%resample( target_total_charge, enforce_total_charge, init_f_params, rand_seed, rank, world_size )
    else
      call self%resample( target_total_charge, enforce_total_charge, init_f_params )
    end if
  end subroutine

  !> resample (and sample) -----------------------------------------------------------------------------------------------------
  !> this routine essentially does 2 things:
  !>  * (re)set the particles on the initial grid
  !>  * computes new interpolation (sparse grid) coefficients so that the remapped_f is an approximation of the target density f
  !>    -- note: if not initial step (true remap), the evaluation of the transported f relies on the lbf method
  subroutine resample( self, target_total_charge, enforce_total_charge, init_f_params, rand_seed, rank, world_size )
    class(sll_t_particle_group_2d2v_lbf),   intent( inout ) :: self
    sll_real64,                       intent( in )    :: target_total_charge
    logical,                          intent( in )    :: enforce_total_charge
    class(sll_c_distribution_params), intent( in ), optional  :: init_f_params   ! < use for initial sampling
    sll_int32, dimension(:)         , intent( in ), optional :: rand_seed
    sll_int32                       , intent( in ), optional :: rank, world_size

    logical :: initial_step

    initial_step = present(init_f_params)

    !> A. write the nodal values of the target density f, then compute the new interpolation coefs
    print *, "particle_group_2d2v_lbf%resample -- step A (LBF approximation of transported density on remapping grid)"

    !>    - A.1  write the nodal values of f on the arrays of interpolation coefs
    if( initial_step )then
      !> write initial density f0 on remapping grid
      call self%write_known_density_on_remapping_grid(init_f_params)

    else
      !> reconstruct transported density fn on remapping grid with the lbf method
      call self%reconstruct_f_lbf_on_remapping_grid()
    end if

    !>  compute the sparse grid interpolation coefs for remapped_f, using the nodal values stored in the tmp array
    call self%sparse_grid_interpolator%compute_hierarchical_surplus( self%tmp_f_values_on_remapping_sparse_grid )
    self%remapped_f_sparse_grid_coefficients = self%tmp_f_values_on_remapping_sparse_grid

    !> B. reset the particles on the cartesian grid
    print *, "particle_group_2d2v_lbf%resample -- step B (deterministic resampling of particles using the remapped density)"
    call self%reset_particles_positions
    !> since the remapping tool has been reset, computing the weights can be done with straightforward interpolation (flow = Id)
    call self%reset_particles_weights_with_direct_interpolation( target_total_charge, enforce_total_charge )

    !> not needed for now, but prevents warnings
    if(present(rank))then
      SLL_ASSERT(present(world_size))
      SLL_ASSERT(present(rand_seed))
    end if

  end subroutine resample

  !> write_known_density_on_remapping_grid -------------------------------------------------------------------------------------
  !> this routine writes a given distribution on the interpolation (sparse grid) nodes to further approximate it
  subroutine write_known_density_on_remapping_grid(    &
      self,        &
      distribution_params       &
    )

    class(sll_t_particle_group_2d2v_lbf), intent(inout) :: self
    class(sll_c_distribution_params),     intent(in)    :: distribution_params

    sll_int32 :: j
    sll_real64 :: x_j(1:2)
    sll_real64 :: v_j(1:2)

    SLL_ASSERT( size(self%tmp_f_values_on_remapping_sparse_grid,1) == self%sparse_grid_interpolator%size_basis )
    self%tmp_f_values_on_remapping_sparse_grid = 0.0_f64

    do j=1, self%sparse_grid_interpolator%size_basis
        x_j = self%sparse_grid_interpolator%hierarchy(j)%coordinate(1:2)  !< here the 2d position
        v_j = self%sparse_grid_interpolator%hierarchy(j)%coordinate(3:4)
        self%tmp_f_values_on_remapping_sparse_grid(j) = distribution_params%eval_xv_density( x_j, v_j )
    end do

  end subroutine write_known_density_on_remapping_grid


  !> reset_particles_positions -------------------------------------------------------------------------------------------------
  !> reset the particles on the initial (cartesian) grid
  !> -- we use the (multi-purpose) lbf grid
  !> -- NOTE: particles will be located at the center of each cell (to avoid a particular treatment of periodic domains)
  subroutine reset_particles_positions( self )
    class(sll_t_particle_group_2d2v_lbf),intent(inout) :: self

    sll_int32 :: k
    sll_int32 :: k_check
    sll_real64, dimension(3)      :: coords

    sll_int32 :: j_x
    sll_int32 :: j_y
    sll_int32 :: j_vx
    sll_int32 :: j_vy
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

    SLL_ASSERT( self%n_particles_x  == self%lbf_grid%num_cells1 )
    SLL_ASSERT( self%n_particles_y  == self%lbf_grid%num_cells2 )
    SLL_ASSERT( self%n_particles_vx == self%lbf_grid%num_cells3 )
    SLL_ASSERT( self%n_particles_vy == self%lbf_grid%num_cells4 )

    h_x    = self%lbf_grid%delta_eta1
    h_y    = self%lbf_grid%delta_eta2
    h_vx   = self%lbf_grid%delta_eta3
    h_vy   = self%lbf_grid%delta_eta4

    x_min    = self%lbf_grid%eta1_min
    y_min    = self%lbf_grid%eta2_min
    vx_min   = self%lbf_grid%eta3_min
    vy_min   = self%lbf_grid%eta4_min

    k_check = 0
    do j_x = 1, self%n_particles_x
      x_j = x_min + (j_x-0.5) * h_x

      do j_y = 1, self%n_particles_y
        y_j = y_min + (j_y-0.5) * h_y

        do j_vx = 1, self%n_particles_vx
          vx_j = vx_min + (j_vx-0.5) * h_vx

          do j_vy = 1, self%n_particles_vy
            vy_j = vy_min + (j_vy-0.5) * h_vy

            k_check = k_check + 1
            call sll_s_convert_4d_index_to_1d( &
              k,  &
              j_x, j_y, j_vx, j_vy,  &
              self%n_particles_x, self%n_particles_y, self%n_particles_vx, self%n_particles_vy  &
            )
            SLL_ASSERT(k == k_check)

            coords(1) = x_j
            coords(2) = y_j
            call self%set_x( k, coords )

            coords(1) = vx_j
            coords(2) = vy_j
            call self%set_v( k, coords )

          end do
        end do
      end do
    end do

  end subroutine reset_particles_positions

  !> get_neighbor_index --------------------------------------------------------------------------------------------------------
  !> returns the index of the specified particle neighbor (defines the particle connectivity)
  function get_neighbor_index( self, k, dim, dir ) result(k_ngb)
    class(sll_t_particle_group_2d2v_lbf),intent(inout) :: self

    sll_int32, intent(in)  :: k     !< particle index
    sll_int32, intent(in)  :: dim   !< neighbor dimension
    sll_int32, intent(in)  :: dir   !< neighbor direction (-1: left, 1:right)
    sll_int32  :: k_ngb   !< neighbor particle index

    sll_int32  :: j_x
    sll_int32  :: j_y
    sll_int32  :: j_vx
    sll_int32  :: j_vy
    sll_int32  :: j_ngb

    call sll_s_convert_1d_index_to_4d( &
      j_x, j_y, j_vx, j_vy,        &
      k, &
      self%n_particles_x, self%n_particles_y,  &
      self%n_particles_vx, self%n_particles_vy  &
    )

    select case( dim )
      case( 1 )
        j_ngb = j_x + dir
        ! correct neighbor index if out of bounds
        if(j_ngb < 1 .or. j_ngb > self%n_particles_x )then
          if( self%domain_is_periodic(1) )then
            j_ngb = modulo(j_ngb-1, self%n_particles_x)+1
          else
            SLL_ASSERT( (j_x == 1 .and. dir == -1) .or. (j_x == self%n_particles_x .and. dir == 1) )
            k_ngb = k ! no neighbor, return the index itself
            return
            SLL_ERROR( "get_neighbor_index", "CHECK 1 -- should not be here")
          end if
        end if
        call sll_s_convert_4d_index_to_1d(  &
            k_ngb, &
            j_ngb, j_y, j_vx, j_vy, &
            self%n_particles_x, self%n_particles_y,  &
            self%n_particles_vx, self%n_particles_vy  &
        )
        return

      case( 2 )
        j_ngb = j_y + dir
        ! correct neighbor index if out of bounds
        if(j_ngb < 1 .or. j_ngb > self%n_particles_y )then
          if( self%domain_is_periodic(2) )then
            j_ngb = modulo(j_ngb-1, self%n_particles_y)+1
          else
            SLL_ASSERT( (j_y == 1 .and. dir == -1) .or. (j_y == self%n_particles_y .and. dir == 1) )
            k_ngb = k ! no neighbor, return the index itself
            return
            SLL_ERROR( "get_neighbor_index", "CHECK 2 -- should not be here")
          end if
        end if
        call sll_s_convert_4d_index_to_1d(  &
            k_ngb, &
            j_x, j_ngb, j_vx, j_vy, &
            self%n_particles_x, self%n_particles_y,  &
            self%n_particles_vx, self%n_particles_vy  &
        )
        return

      case( 3 )
        j_ngb = j_vx + dir
        ! correct neighbor index if out of bounds
        if(j_ngb < 1 .or. j_ngb > self%n_particles_vx )then
          ! no periodicity along vx
          SLL_ASSERT( (j_vx == 1 .and. dir == -1) .or. (j_vx == self%n_particles_vx .and. dir == 1) )
          k_ngb = k ! no neighbor, return the index itself
          return
          SLL_ERROR( "get_neighbor_index", "CHECK 3 -- should not be here")
        end if
        call sll_s_convert_4d_index_to_1d(  &
            k_ngb, &
            j_x, j_y, j_ngb, j_vy, &
            self%n_particles_x, self%n_particles_y,  &
            self%n_particles_vx, self%n_particles_vy  &
        )
        return

      case( 4 )
        j_ngb = j_vy + dir
        ! correct neighbor index if out of bounds
        if(j_ngb < 1 .or. j_ngb > self%n_particles_vy )then
          ! no periodicity along vy
          SLL_ASSERT( (j_vy == 1 .and. dir == -1) .or. (j_vy == self%n_particles_vy .and. dir == 1) )
          k_ngb = k ! no neighbor, return the index itself
          return
          SLL_ERROR( "get_neighbor_index", "CHECK 4 -- should not be here")
        end if
        call sll_s_convert_4d_index_to_1d(  &
            k_ngb, &
            j_x, j_y, j_vx, j_ngb, &
            self%n_particles_x, self%n_particles_y,  &
            self%n_particles_vx, self%n_particles_vy  &
        )
        return

      case default
          SLL_ERROR("get_neighbor_index", "wrong value for dim")
    end select

  end function

  !> reset_particles_weights_with_direct_interpolation -------------------------------------------------------------------------
  !> compute (and set) the particle weights using a direct interpolation (here, sparse grid)
  !> -- does not change the position of the deposition particles
  subroutine reset_particles_weights_with_direct_interpolation(  &
      self,                             &
      target_total_charge,              &
      enforce_total_charge              &
  )
    class( sll_t_particle_group_2d2v_lbf ),           intent( inout ) :: self
    sll_real64,                                 intent( in )    :: target_total_charge
    logical,                                    intent( in )    :: enforce_total_charge
    sll_real64    :: eta(4)
    sll_int32     :: i_part
    sll_real64    :: phase_space_volume
    sll_real64    :: particle_density_factor
    sll_real64    :: point_density, total_computed_density
    sll_real64    :: total_computed_charge
    sll_real64    :: charge_correction_factor

    !> set the particle group common_weight to 1 (will be a factor of each particle weight)
    call self%set_common_weight(1.0_f64)

    phase_space_volume =    (self%lbf_grid%eta4_max - self%lbf_grid%eta4_min)    &
                          * (self%lbf_grid%eta3_max - self%lbf_grid%eta3_min)    &
                          * (self%lbf_grid%eta2_max - self%lbf_grid%eta2_min)    &
                          * (self%lbf_grid%eta1_max - self%lbf_grid%eta1_min)
    particle_density_factor = phase_space_volume / self%n_particles
    total_computed_density = 0.0d0
    do i_part = 1, self%n_particles
      eta = self%particle_array(1:4, i_part)
      point_density = particle_density_factor * self%interpolate_value_of_remapped_f(eta)
      self%particle_array(5, i_part) = point_density   !> 5 is the first weight index
      total_computed_density = total_computed_density + point_density
    end do

    if( enforce_total_charge )then
      total_computed_charge = self%species%q * total_computed_density
      if( total_computed_density == 0 )then
        print *, "WARNING (876786587689) -- total computed_density is zero, which is strange..."
        print *, "                       -- (no charge correction in this case) "
      else
        charge_correction_factor = target_total_charge / total_computed_charge

        print *, "[Enforcing charge in reset_particles_weights_with_direct_interpolation] ... "
        print *, "   ...   target_total_charge, total_computed_charge, charge_correction_factor = ", &
            target_total_charge, total_computed_charge, charge_correction_factor

        ! todo: try with    self%particle_array(5, :) = self%particle_array(5, :) * charge_correction_factor
        do i_part = 1, self%n_particles
          self%particle_array(5, i_part) = self%particle_array(5, i_part) * charge_correction_factor !> 5 is the first weight index
        end do
      end if
    end if

  end subroutine reset_particles_weights_with_direct_interpolation

  !> interpolate_value_of_remapped_f -------------------------------------------------------------------------------------------
  !> computes the value of the (last) remapped density f using the interpolation tool
  !>  -- for now, assuming only sparse grid interpolation
  function interpolate_value_of_remapped_f ( self, eta ) result(val)
    class(sll_t_particle_group_2d2v_lbf), intent(inout)  :: self
    sll_real64, dimension(4),       intent(in)  :: eta           !< Position where to interpolate
    sll_real64                                  :: val

    val = self%sparse_grid_interpolator%interpolate_from_interpolant_value(self%remapped_f_sparse_grid_coefficients, eta)

  end function

  !> macro for the lbf approximation of a transported density (reconstruct_f_lbf) below ----------------------------------------
#define UPDATE_CLOSEST_PARTICLE_ARRAYS_USING_NEIGHBOR_CELLS(djx,djy,djvx,djvy)                                            \
    do;                                                                                                                 \
        k_neighbor = closest_particle(j_x+(djx), j_y+(djy), j_vx+(djvx), j_vy+(djvy));                                    \
;                                                                                                                       \
        if(k_neighbor /= 0) then;  do          ;                                                                        \
            coords = self%get_x(k_neighbor) ;                                                                        \
            x = coords(1) ;                                                                                             \
            y = coords(2) ;                                                                                             \
            coords = self%get_v(k_neighbor) ;                                                                        \
            vx = coords(1) ;                                                                                            \
            vy = coords(2) ;                                                                                            \
            call periodic_correction(self,x,y) ;                                                                     \
            call sll_s_update_closest_particle_arrays(k_neighbor,                                                       \
                x - self%lbf_grid%eta1_min, y - self%lbf_grid%eta2_min, vx - self%lbf_grid%eta3_min, vy - self%lbf_grid%eta4_min, \
                j_x, j_y, j_vx, j_vy,                                                   \
                h_flow_grid_x,                                                          \
                h_flow_grid_y,                                                          \
                h_flow_grid_vx,                                                         \
                h_flow_grid_vy,                                                         \
                closest_particle,                                                       \
                closest_particle_distance) ;                                            \
        exit;                                                                           \
        end do;                                                                         \
        end if;                                                                         \
    exit;                                                                               \
    end do



  !> reconstruct_f_lbf  ----------------------------------------------------------------------------------------------------------
  !> reconstruct point values of the transported f using a bsl approximation based on linearized backward flows (lbf)
  !> values are reconstructed on different grids, depending on reconstruction_set_type:
  !>  - sll_p_lbf_remapping_grid: the remapping grid (for now, a sparse grid)
  !>  - sll_p_lbf_given_grid: a 4d cartesian grid

  subroutine reconstruct_f_lbf( &
    self,    &
    reconstruction_set_type,    &
    given_grid_4d,   &
    given_array_4d, &
    reconstruct_f_on_last_node,  &
    target_total_charge,   &
    enforce_total_charge   &
  )
    class(sll_t_particle_group_2d2v_lbf),     intent(inout) :: self          !> particle group
    sll_int32,                                intent(in)    :: reconstruction_set_type
    type(sll_t_cartesian_mesh_4d),   pointer, intent(in)    :: given_grid_4d
    sll_real64, dimension(:,:,:,:),  pointer, intent(inout) :: given_array_4d
    logical,                                  intent(in)    :: reconstruct_f_on_last_node(4)
    sll_real64,                               intent(in)    :: target_total_charge
    logical,                                  intent(in)    :: enforce_total_charge

    sll_real64 :: deposited_charge
    sll_real64 :: charge_correction_factor

    !> index of particle closest to the center of each flow cell
    sll_int32,          dimension(:,:,:,:), allocatable     :: closest_particle
    sll_real64,         dimension(:,:,:,:), allocatable     :: closest_particle_distance

    !> array of integer linked lists (declared below) useful when remapping on sparse grids
    !> (can also simplify the code when remapping on cartesian grids which do not match with the flow cells? but maybe too costly)
    type(sll_t_int_list_element_ptr), dimension(:,:,:,:), allocatable     :: nodes_in_flow_cell
    type(sll_t_int_list_element),     pointer                             :: new_int_list_element, head

    sll_int32  :: k ! particle index
    sll_int32  :: k_neighbor
    sll_int32  :: i_x, i_y, i_vx, i_vy  !< grid node indices
    sll_int32  :: j_x, j_y, j_vx, j_vy  !< flow cell indices
    sll_int32  :: m_x, m_y, m_vx, m_vy
    sll_int32  :: j_tmp

    sll_int32  :: i_min_x
    sll_int32  :: i_min_y
    sll_int32  :: i_min_vx
    sll_int32  :: i_min_vy

    sll_int32  :: i_max_x
    sll_int32  :: i_max_y
    sll_int32  :: i_max_vx
    sll_int32  :: i_max_vy

    sll_int32  :: node_index
    sll_int32  :: k_particle_closest_to_first_corner

    !> <<g>> cartesian grid pointer to the remapping grid
    type(sll_t_cartesian_mesh_4d),pointer :: g

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

    !> the flow cells are the cells where the flow will be linearized
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

    sll_int32  :: flow_grid_j_x_min
    sll_int32  :: flow_grid_j_x_max
    sll_int32  :: flow_grid_j_y_min
    sll_int32  :: flow_grid_j_y_max
    sll_int32  :: flow_grid_j_vx_min
    sll_int32  :: flow_grid_j_vx_max
    sll_int32  :: flow_grid_j_vy_min
    sll_int32  :: flow_grid_j_vy_max

    sll_real64 :: x
    sll_real64 :: y
    sll_real64 :: vx
    sll_real64 :: vy

    sll_real64 :: closest_particle_distance_to_first_corner
    sll_real64 :: particle_distance_to_first_corner

    sll_real64 :: mesh_period_x
    sll_real64 :: mesh_period_y

    !> results from [[get_ltp_deformation_matrix]]
    sll_real64 :: d11,d12,d13,d14 !< coefs of matrix D (backward Jacobian)
    sll_real64 :: d21,d22,d23,d24
    sll_real64 :: d31,d32,d33,d34
    sll_real64 :: d41,d42,d43,d44

    sll_real64, dimension(3)  :: coords
    sll_real64, dimension(4)  :: eta        !< coordinates in the logical (cartesian) 4d space

    !> coordinates of particle k at time n and time 0
    sll_real64 :: x_k,y_k,vx_k,vy_k
    sll_real64 :: x_k_t0,y_k_t0,vx_k_t0,vy_k_t0

    sll_real64 :: x_to_xk, y_to_yk, vx_to_vxk, vy_to_vyk

    sll_real64 :: d1_x, d1_y, d1_vx, d1_vy
    sll_real64 :: d2_x, d2_y, d2_vx, d2_vy
    sll_real64 :: d3_x, d3_y, d3_vx, d3_vy
    sll_real64 :: d4_x, d4_y, d4_vx, d4_vy

    sll_int32  :: nodes_number
    sll_real64 :: reconstructed_f_value

    !> coordinates of a reconstruction point at time 0, absolute
    sll_real64 :: x_t0
    sll_real64 :: y_t0
    sll_real64 :: vx_t0
    sll_real64 :: vy_t0

    !> coordinates of a reconstruction point at time 0, relative to the nearby reference particle (= closest to cell center)
    sll_real64 :: x_t0_to_xk_t0
    sll_real64 :: y_t0_to_yk_t0
    sll_real64 :: vx_t0_to_vxk_t0
    sll_real64 :: vy_t0_to_vyk_t0

    sll_int32  :: ierr

    !> getting the parameters of the flow grid (where the bwd flow is linearized)
    flow_grid_x_min    = self%lbf_grid%eta1_min
    flow_grid_y_min    = self%lbf_grid%eta2_min
    flow_grid_vx_min   = self%lbf_grid%eta3_min
    flow_grid_vy_min   = self%lbf_grid%eta4_min

    h_flow_grid_x  = self%lbf_grid%delta_eta1
    h_flow_grid_y  = self%lbf_grid%delta_eta2
    h_flow_grid_vx = self%lbf_grid%delta_eta3
    h_flow_grid_vy = self%lbf_grid%delta_eta4

    flow_grid_num_cells_x  = self%lbf_grid%num_cells1
    flow_grid_num_cells_y  = self%lbf_grid%num_cells2
    flow_grid_num_cells_vx = self%lbf_grid%num_cells3
    flow_grid_num_cells_vy = self%lbf_grid%num_cells4

    !> reconstruction will be done through local linearization of the bwd flow in flow cells.
    !> these flow cells are defined by the lbf_grid but we may not need all of them (if reconstructing f on a given grid, see below)
    flow_grid_j_x_min = 1
    flow_grid_j_x_max = flow_grid_num_cells_x
    flow_grid_j_y_min = 1
    flow_grid_j_y_max = flow_grid_num_cells_y
    flow_grid_j_vx_min = 1
    flow_grid_j_vx_max = flow_grid_num_cells_vx
    flow_grid_j_vy_min = 1
    flow_grid_j_vy_max = flow_grid_num_cells_vy

    deposited_charge = 0.0_f64

    !> A.  preparation of the point sets where f will be reconstructed, depending on the cases
    if( reconstruction_set_type == sll_p_lbf_remapping_grid )then

      !> => prepare the array of linked lists that will store the node indices contained in the flow cells (one list per cell)
      allocate(nodes_in_flow_cell(flow_grid_num_cells_x,   &
                                  flow_grid_num_cells_y,   &
                                  flow_grid_num_cells_vx,  &
                                  flow_grid_num_cells_vy)  &
             , stat=ierr)
      call sll_s_test_error_code(ierr, 'Memory allocation Failure.', __FILE__, __LINE__)

      do j_x = 1, flow_grid_num_cells_x
        do j_y = 1, flow_grid_num_cells_y
          do j_vx = 1, flow_grid_num_cells_vx
            do j_vy = 1, flow_grid_num_cells_vy
              nullify(nodes_in_flow_cell(j_x,j_y,j_vx,j_vy)%pointed_element)
            end do
          end do
        end do
      end do

      nodes_number = self%sparse_grid_interpolator%size_basis
      SLL_ASSERT( size(self%tmp_f_values_on_remapping_sparse_grid,1) == nodes_number )
      self%tmp_f_values_on_remapping_sparse_grid = 0.0_f64

      !> then loop to store the sparse grid node indices in linked lists corresponding to the flow cells that contain them
      do node_index = 1, nodes_number

        !> get node coordinates:
        x = self%sparse_grid_interpolator%hierarchy(node_index)%coordinate(1)
        y = self%sparse_grid_interpolator%hierarchy(node_index)%coordinate(2)
        vx = self%sparse_grid_interpolator%hierarchy(node_index)%coordinate(3)
        vy = self%sparse_grid_interpolator%hierarchy(node_index)%coordinate(4)

        !> find the index (j_x,j_y,j_vx,j_vy) of the flow cell containing this node (same piece of code as below)
        call sll_s_get_1d_cell_containing_point(j_x,  x,  flow_grid_x_min,  h_flow_grid_x )
        call sll_s_get_1d_cell_containing_point(j_y,  y,  flow_grid_y_min,  h_flow_grid_y )
        call sll_s_get_1d_cell_containing_point(j_vx, vx, flow_grid_vx_min, h_flow_grid_vx)
        call sll_s_get_1d_cell_containing_point(j_vy, vy, flow_grid_vy_min, h_flow_grid_vy)

        !> discard if flow cell is off-bounds
        if(  j_x  >= flow_grid_j_x_min  .and. j_x  <= flow_grid_j_x_max  .and. &
             j_y  >= flow_grid_j_y_min  .and. j_y  <= flow_grid_j_y_max  .and. &
             j_vx >= flow_grid_j_vx_min .and. j_vx <= flow_grid_j_vx_max .and. &
             j_vy >= flow_grid_j_vy_min .and. j_vy <= flow_grid_j_vy_max  )then

          !> increment the proper linked list
          SLL_ALLOCATE( new_int_list_element, ierr )
          new_int_list_element%value = node_index
          head => nodes_in_flow_cell(j_x,j_y,j_vx,j_vy)%pointed_element
          nodes_in_flow_cell(j_x,j_y,j_vx,j_vy)%pointed_element => sll_f_add_element_in_int_list(head, new_int_list_element)
        end if

      end do

      nullify(g)
      g_num_points_x  = 0
      g_num_points_y  = 0
      g_num_points_vx = 0
      g_num_points_vy = 0

    else
      SLL_ASSERT( reconstruction_set_type == sll_p_lbf_given_grid )

      !> then use the given 4d grid and write values in given array
      g => given_grid_4d

      !> if the given grid is strictly contained in the lbf grid, we can reduce the number of flow cells where f is reconstructed
      call sll_s_get_1d_cell_containing_point( j_tmp, g%eta1_min, flow_grid_x_min, h_flow_grid_x )
      flow_grid_j_x_min = max(flow_grid_j_x_min, j_tmp)
      call sll_s_get_1d_cell_containing_point( j_tmp, g%eta1_max, flow_grid_x_min, h_flow_grid_x )
      flow_grid_j_x_max = min(flow_grid_j_x_max, j_tmp)

      call sll_s_get_1d_cell_containing_point( j_tmp, g%eta2_min, flow_grid_y_min, h_flow_grid_y )
      flow_grid_j_y_min = max(flow_grid_j_y_min, j_tmp)
      call sll_s_get_1d_cell_containing_point( j_tmp, g%eta2_max, flow_grid_y_min, h_flow_grid_y )
      flow_grid_j_y_max = min(flow_grid_j_y_max, j_tmp)

      call sll_s_get_1d_cell_containing_point( j_tmp, g%eta3_min, flow_grid_vx_min, h_flow_grid_vx )
      flow_grid_j_vx_min = max(flow_grid_j_vx_min, j_tmp)
      call sll_s_get_1d_cell_containing_point( j_tmp, g%eta3_max, flow_grid_vx_min, h_flow_grid_vx )
      flow_grid_j_vx_max = min(flow_grid_j_vx_max, j_tmp)

      call sll_s_get_1d_cell_containing_point( j_tmp, g%eta4_min, flow_grid_vy_min, h_flow_grid_vy )
      flow_grid_j_vy_min = max(flow_grid_j_vy_min, j_tmp)
      call sll_s_get_1d_cell_containing_point( j_tmp, g%eta4_max, flow_grid_vy_min, h_flow_grid_vy )
      flow_grid_j_vy_max = min(flow_grid_j_vy_max, j_tmp)

      g_num_points_x = g%num_cells1
      if( reconstruct_f_on_last_node(1) ) g_num_points_x = g_num_points_x + 1
      SLL_ASSERT( size(given_array_4d, 1) == g_num_points_x  )

      g_num_points_y = g%num_cells2
      if( reconstruct_f_on_last_node(2) ) g_num_points_y = g_num_points_y + 1
      SLL_ASSERT( size(given_array_4d, 2) == g_num_points_y  )

      g_num_points_vx = g%num_cells3
      if( reconstruct_f_on_last_node(3) ) g_num_points_vx = g_num_points_vx + 1
      SLL_ASSERT( size(given_array_4d, 3) == g_num_points_vx )

      g_num_points_vy = g%num_cells4
      if( reconstruct_f_on_last_node(4) ) g_num_points_vy = g_num_points_vy + 1
      SLL_ASSERT( size(given_array_4d, 4) == g_num_points_vy )
      given_array_4d(:,:,:,:) = 0.0_f64

      h_g_grid_x  = g%delta_eta1
      h_g_grid_y  = g%delta_eta2
      h_g_grid_vx = g%delta_eta3
      h_g_grid_vy = g%delta_eta4

      g_grid_x_min  = g%eta1_min
      g_grid_y_min  = g%eta2_min
      g_grid_vx_min = g%eta3_min
      g_grid_vy_min = g%eta4_min

    end if

    !> B. Preparatory work for the linearization of the flow on the flow cells:
    !>    find out the closest particle to each cell center,
    !>    by looping over all particles and noting which flow cell contains it.
    !>    (The leftmost flow cell in each dimension may not be complete.)
    SLL_ALLOCATE(closest_particle(flow_grid_num_cells_x,flow_grid_num_cells_y,flow_grid_num_cells_vx,flow_grid_num_cells_vy),ierr)
    closest_particle(:,:,:,:) = 0

    allocate(closest_particle_distance(flow_grid_num_cells_x,   &
                                       flow_grid_num_cells_y,   &
                                       flow_grid_num_cells_vx,  &
                                       flow_grid_num_cells_vy)  , stat=ierr)
    call sll_s_test_error_code(ierr, 'Memory allocation Failure.', __FILE__, __LINE__)
    closest_particle_distance(:,:,:,:) = 0.0_f64

    closest_particle_distance_to_first_corner = 1d30
    k_particle_closest_to_first_corner = 0

    do k=1, self%n_particles

      !> find absolute (x_k,y_k,vx_k,vy_k) coordinates for k-th particle.
      coords = self%get_x(k)
      x_k = coords(1)
      y_k = coords(2)
      coords = self%get_v(k)
      vx_k = coords(1)
      vy_k = coords(2)

      !> which flow cell is this particle in?
      call sll_s_get_1d_cell_containing_point(j_x, x_k,  flow_grid_x_min,  h_flow_grid_x )
      call sll_s_get_1d_cell_containing_point(j_y, y_k,  flow_grid_y_min,  h_flow_grid_y )
      call sll_s_get_1d_cell_containing_point(j_vx, vx_k, flow_grid_vx_min, h_flow_grid_vx)
      call sll_s_get_1d_cell_containing_point(j_vy, vy_k, flow_grid_vy_min, h_flow_grid_vy)

      !> discard this particle if off-bounds
      if(  j_x  >= flow_grid_j_x_min  .and. j_x  <= flow_grid_j_x_max  .and. &
           j_y  >= flow_grid_j_y_min  .and. j_y  <= flow_grid_j_y_max  .and. &
           j_vx >= flow_grid_j_vx_min .and. j_vx <= flow_grid_j_vx_max .and. &
           j_vy >= flow_grid_j_vy_min .and. j_vy <= flow_grid_j_vy_max )then

        call sll_s_update_closest_particle_arrays(  &
          k,  &
          x_k  - flow_grid_x_min, &
          y_k  - flow_grid_y_min, &
          vx_k - flow_grid_vx_min, &
          vy_k - flow_grid_vy_min, &
          j_x, j_y, j_vx, j_vy,  &
          h_flow_grid_x,  &
          h_flow_grid_y,  &
          h_flow_grid_vx,  &
          h_flow_grid_vy,  &
          closest_particle,  &
          closest_particle_distance &
        )
      end if

      particle_distance_to_first_corner =   abs(x_k  - flow_grid_x_min )  &
                                          + abs(y_k  - flow_grid_y_min )  &
                                          + abs(vx_k - flow_grid_vx_min)  &
                                          + abs(vy_k - flow_grid_vy_min)
      if( particle_distance_to_first_corner < closest_particle_distance_to_first_corner )then
        closest_particle_distance_to_first_corner = particle_distance_to_first_corner
        k_particle_closest_to_first_corner = k
      end if
    end do

    closest_particle(1,1,1,1) = k_particle_closest_to_first_corner

    if( .not. ( self%domain_is_periodic(1) .and. self%domain_is_periodic(2) ) )then
      print*, "WARNING -- STOP -- verify that the non-periodic case is well implemented"
      stop
    end if

    !> C. loop on the flow cells (main loop):
    !>   on each flow cell, we
    !>   - C.1 linearize the flow using the position of the particles
    !>   - C.2 find the relevant points where f should be reconstructed
    !>   - C.3 reconstruct f on these points (using the affine backward flow and the interpolation tool for the remapped_f)
    !>   - C.4 write the resulting f value on the proper place

    if(self%domain_is_periodic(1)) then
      !> here the domain corresponds to the Poisson mesh
      mesh_period_x = self%lbf_grid%eta1_max - self%lbf_grid%eta1_min
    else
      mesh_period_x = 0.0_f64
    end if

    if(self%domain_is_periodic(2)) then
      !> here the domain corresponds to the Poisson mesh
      mesh_period_y = self%lbf_grid%eta2_max - self%lbf_grid%eta2_min
    else
      mesh_period_y = 0.0_f64
    end if

    do j_x = flow_grid_j_x_min, flow_grid_j_x_max
      do j_y = flow_grid_j_y_min, flow_grid_j_y_max
        do j_vx = flow_grid_j_vx_min, flow_grid_j_vx_max
          do j_vy = flow_grid_j_vy_min, flow_grid_j_vy_max

            !> Find the particle which is closest to the cell center -- and if we do not have it yet, look on the neighboring cells
            k = closest_particle(j_x,j_y,j_vx,j_vy)

            if(k == 0) then
              if( j_x > 1 )then
                UPDATE_CLOSEST_PARTICLE_ARRAYS_USING_NEIGHBOR_CELLS(-1,0,0,0)
              end if
              if( j_x < flow_grid_num_cells_x )then
                UPDATE_CLOSEST_PARTICLE_ARRAYS_USING_NEIGHBOR_CELLS( 1,0,0,0)
              end if

              if( j_y > 1 )then
                UPDATE_CLOSEST_PARTICLE_ARRAYS_USING_NEIGHBOR_CELLS(0,-1,0,0)
              end if
              if( j_y < flow_grid_num_cells_y )then
                UPDATE_CLOSEST_PARTICLE_ARRAYS_USING_NEIGHBOR_CELLS(0, 1,0,0)
              end if

              if( j_vx > 1 )then
                UPDATE_CLOSEST_PARTICLE_ARRAYS_USING_NEIGHBOR_CELLS(0,0,-1,0)
              end if
              if( j_vx < flow_grid_num_cells_vx )then
                UPDATE_CLOSEST_PARTICLE_ARRAYS_USING_NEIGHBOR_CELLS(0,0, 1,0)
              end if

              if( j_vy > 1 )then
                UPDATE_CLOSEST_PARTICLE_ARRAYS_USING_NEIGHBOR_CELLS(0,0,0,-1)
              end if
              if( j_vy < flow_grid_num_cells_vy )then
                UPDATE_CLOSEST_PARTICLE_ARRAYS_USING_NEIGHBOR_CELLS(0,0,0, 1)
              end if
            end if

            k = closest_particle(j_x,j_y,j_vx,j_vy)
            SLL_ASSERT(k /= 0)

            !>   - C.1 linearize the flow using the position of the particles

            !> In this flow cell we will use the k-th backward flow, requires the deformation matrix of the k-th particle

            call self%get_ltp_deformation_matrix (       &
                 k,                                         &
                 mesh_period_x,                             &
                 mesh_period_y,                             &
                 h_flow_grid_x,                               &
                 h_flow_grid_y,                               &
                 h_flow_grid_vx,                              &
                 h_flow_grid_vy,                              &
                 x_k,y_k,vx_k,vy_k,                         &
                 d11,d12,d13,d14,                           &
                 d21,d22,d23,d24,                           &
                 d31,d32,d33,d34,                           &
                 d41,d42,d43,d44                            &
                 )

            !> Find position of particle k at time 0
            call sll_s_convert_1d_index_to_4d(  &
              m_x,m_y,m_vx,m_vy, &
              k, &
              self%n_particles_x, self%n_particles_y,  &
              self%n_particles_vx, self%n_particles_vy  &
            )
            !> initially the particles are on the center of the lbf cells (also the flow grid)
            x_k_t0  = flow_grid_x_min  + (m_x-0.5)  * h_flow_grid_x
            y_k_t0  = flow_grid_y_min  + (m_y-0.5)  * h_flow_grid_y
            vx_k_t0 = flow_grid_vx_min + (m_vx-0.5) * h_flow_grid_vx
            vy_k_t0 = flow_grid_vy_min + (m_vy-0.5) * h_flow_grid_vy


            !>  - C.2 find the relevant points (x, y, vx, vy) where f should be reconstructed
            !>    - C.2.a first we treat the case of [remapping with a sparse grid]
            !>            or [computing the weights of a cloud of deposition particles]: nodes are stored in linked lists

            if( reconstruction_set_type == sll_p_lbf_remapping_grid )then

              new_int_list_element => nodes_in_flow_cell(j_x,j_y,j_vx,j_vy)%pointed_element

              do while( associated(new_int_list_element) )
                node_index = new_int_list_element%value

                x  = self%sparse_grid_interpolator%hierarchy(node_index)%coordinate(1)
                y  = self%sparse_grid_interpolator%hierarchy(node_index)%coordinate(2)
                vx = self%sparse_grid_interpolator%hierarchy(node_index)%coordinate(3)
                vy = self%sparse_grid_interpolator%hierarchy(node_index)%coordinate(4)

                !> node coordinates relative to current particle position
                x_to_xk   = x - x_k
                y_to_yk   = y - y_k
                vx_to_vxk = vx - vx_k
                vy_to_vyk = vy - vy_k

                !> find z_t0 = (x_t0, y_t0, vx_t0, vy_t0), the position of the node at time = 0 using the affine bwd flow
                x_t0_to_xk_t0   = d11 * x_to_xk + d12 * y_to_yk + d13 * vx_to_vxk + d14 * vy_to_vyk
                y_t0_to_yk_t0   = d21 * x_to_xk + d22 * y_to_yk + d23 * vx_to_vxk + d24 * vy_to_vyk
                vx_t0_to_vxk_t0 = d31 * x_to_xk + d32 * y_to_yk + d33 * vx_to_vxk + d34 * vy_to_vyk
                vy_t0_to_vyk_t0 = d41 * x_to_xk + d42 * y_to_yk + d43 * vx_to_vxk + d44 * vy_to_vyk

                x_t0  = x_t0_to_xk_t0   + x_k_t0
                y_t0  = y_t0_to_yk_t0   + y_k_t0
                vx_t0 = vx_t0_to_vxk_t0 + vx_k_t0
                vy_t0 = vy_t0_to_vyk_t0 + vy_k_t0

                !> put back z_t0 inside domain (if needed) to enforce periodic boundary conditions
                call periodic_correction(self,x_t0,y_t0)

                !> now (x_t0, y_t0, vx_t0, vy_t0) is the (approx) position of the node at time t=0
                eta(1) = x_t0
                eta(2) = y_t0
                eta(3) = vx_t0
                eta(4) = vy_t0

                !> interpolate and store nodal value on ad-hoc array
                reconstructed_f_value = self%interpolate_value_of_remapped_f(eta)
                self%tmp_f_values_on_remapping_sparse_grid(node_index) = reconstructed_f_value

                new_int_list_element => new_int_list_element%next

              end do

            !>    - C.2.b next we treat the case of a remapping with splines or writing on a given grid
            !>      (in both cases the nodes are on the g grid constructed above)

            else
              SLL_ASSERT( reconstruction_set_type == sll_p_lbf_given_grid )

              !> loop over the grid points inside this flow cell:
              !>
              !> points in the grid are of the form  d(i) = g_grid_d_min + (i-1) * h_g_grid_d,   i = 1, .. g_num_points_d
              !> (with  d = x, y, vx or vy  and  g_num_points_d = g_num_cells_d  or  g_num_cells_d+1 , depending on the periodicity)
              !> and this flow cell has the form     [flow_grid_d_min + (j-1) * h_flow_grid_d, flow_grid_d_min + j * h_flow_grid_d[
              !> so eta_d(i) is in this flow cell if  i_min <= i <= i_max
              !> where i_min = ceiling( (flow_grid_min - g_grid_d_min + (j-1)*h_flow_grid_d)/h_g_grid_d + 1 )
              !> and   i_max = ceiling( (flow_grid_min - g_grid_d_min + (j)  *h_flow_grid_d)/h_g_grid_d + 1 ) - 1

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

                do i_y = i_min_y, i_max_y
                  y = g_grid_y_min + (i_y-1)*h_g_grid_y
                  y_to_yk = y - y_k
                  d1_y = d12 * y_to_yk
                  d2_y = d22 * y_to_yk
                  d3_y = d32 * y_to_yk
                  d4_y = d42 * y_to_yk

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

                      !> put back z_t0 inside domain (if needed) to enforce periodic boundary conditions
                      call periodic_correction(self,x_t0,y_t0)

                      !> now (x_t0, y_t0, vx_t0, vy_t0) is the (approx) position of the node z_i at time t=0
                      eta(1) = x_t0
                      eta(2) = y_t0
                      eta(3) = vx_t0
                      eta(4) = vy_t0

                      !> interpolation here may use sparse grid or splines, depending on the method chosen for self
                      reconstructed_f_value = self%interpolate_value_of_remapped_f(eta)

                      ! [DEBUG]
                      if( .false. )then
                        if( reconstruction_set_type == sll_p_lbf_given_grid )then
                          print*, "[WRITE ON GIVEN GRID]    reconstructing "
                          print*, "on:                       eta = ", eta
                          print*, "value:  reconstructed_f_value = ", reconstructed_f_value
                        end if
                      end if

                      !> C.4  write the reconstructed value at proper place

                      if( reconstructed_f_value /= 0 )then

                        SLL_ASSERT( 1 <= i_x  .and. i_x  <= g_num_points_x  )
                        SLL_ASSERT( 1 <= i_y  .and. i_y  <= g_num_points_y  )
                        SLL_ASSERT( 1 <= i_vx .and. i_vx <= g_num_points_vx )
                        SLL_ASSERT( 1 <= i_vy .and. i_vy <= g_num_points_vy )

                        if( abs(given_array_4d(i_x, i_y, i_vx, i_vy)) > 0.00001 )then
                          print *, "Warning -- 987698666999979979 -- stored value should be zero"
                          print *, "given_array_4d(i_x, i_y, i_vx, i_vy) = ", given_array_4d(i_x, i_y, i_vx, i_vy)
                          print *, "with (i_x, i_y, i_vx, i_vy) = ", i_x, i_y, i_vx, i_vy
                        end if
                        given_array_4d(i_x, i_y, i_vx, i_vy) = reconstructed_f_value

                      else
                        print *, "Warning -- 654654535466545434564 -- just reconstructed a Zero value"
                      end if
                    !> this is the end of the (fourfold) loop on the grid nodes
                    end do
                  end do
                end do
              end do
            end if
          !> and this is the end of (fourfold) loop on the flow cells
          end do
        end do
      end do
    end do

    if( enforce_total_charge )then
      SLL_ASSERT( reconstruction_set_type == sll_p_lbf_given_grid )
      if( deposited_charge == 0 )then
        print *, "WARNING (76576537475) -- total deposited charge is zero, which is strange..."
        print *, "                      -- (no charge correction in this case) "
      else
        charge_correction_factor = target_total_charge / deposited_charge
        print *, "[Enforcing charge]: target_total_charge, deposited_charge, charge_correction_factor = ", &
                    target_total_charge, deposited_charge, charge_correction_factor
        do i_x = 1, g_num_points_x
          do i_y = 1, g_num_points_y
            do i_vx = 1, g_num_points_vx
              do i_vy = 1, g_num_points_vy
                given_array_4d(i_x, i_y, i_vx, i_vy) = given_array_4d(i_x, i_y, i_vx, i_vy) &
                                                       * charge_correction_factor
              end do
            end do
          end do
        end do
      end if
    end if

  end subroutine reconstruct_f_lbf

  !> reconstruct_f_lbf_on_remapping_grid ---------------------------------------------------------------------------------------
  !> layer for reconstruct_f_lbf, when reconstructing f on the remapping grid
  subroutine reconstruct_f_lbf_on_remapping_grid( self )
    class(sll_t_particle_group_2d2v_lbf),     intent(inout) :: self          !> particle group
    type(sll_t_cartesian_mesh_4d),   pointer  :: void_grid_4d
    sll_real64, dimension(:,:,:,:),  pointer  :: void_array_4d
    logical    :: dummy_reconstruct_f_on_last_node(4)
    sll_real64 :: dummy_target_total_charge
    logical    :: dummy_enforce_total_charge

    nullify(void_grid_4d)
    nullify(void_array_4d)
    dummy_reconstruct_f_on_last_node = .false.
    dummy_target_total_charge = 0.0_f64
    dummy_enforce_total_charge = .false.

    call self%reconstruct_f_lbf( &
      sll_p_lbf_remapping_grid, &
      void_grid_4d, &
      void_array_4d, &
      dummy_reconstruct_f_on_last_node, &
      dummy_target_total_charge,  &
      dummy_enforce_total_charge )

  end subroutine reconstruct_f_lbf_on_remapping_grid


  !> reconstruct_f_lbf_on_remapping_grid ---------------------------------------------------------------------------------------
  !> layer for reconstruct_f_lbf, when reconstructing f on a given grid
  subroutine reconstruct_f_lbf_on_given_grid( &
    self,    &
    given_grid_4d,   &
    given_array_4d, &
    reconstruct_f_on_last_node,  &
    target_total_charge,   &
    enforce_total_charge   &
  )
    class(sll_t_particle_group_2d2v_lbf),     intent(inout) :: self          !> particle group
    type(sll_t_cartesian_mesh_4d),   pointer, intent(in)    :: given_grid_4d
    sll_real64, dimension(:,:,:,:),  pointer, intent(inout) :: given_array_4d
    logical,                                  intent(in)    :: reconstruct_f_on_last_node(4)
    sll_real64,                               intent(in)    :: target_total_charge
    logical,                                  intent(in)    :: enforce_total_charge

    call self%reconstruct_f_lbf( &
      sll_p_lbf_given_grid, &
      given_grid_4d, &
      given_array_4d, &
      reconstruct_f_on_last_node, &
      target_total_charge,  &
      enforce_total_charge )

  end subroutine reconstruct_f_lbf_on_given_grid


  !> get_ltp_deformation_matrix ------------------------------------------------------------------------------------------------
  !> Compute the coefficients of the particle 'deformation' matrix
  !> which approximates the Jacobian matrix of the exact backward flow
  !> (defined between the current time and that of the particle initialization -- or the last particle remapping)
  !> at the current particle position.
  !> It is computed in two steps:
  !>    * first we compute a finite-difference approximation of the forward Jacobian matrix by using the fact that the
  !>      initial (or remapped) particles were located on a cartesian grid in phase space. Specifically, the entries of
  !>      the forward Jacobian matrix (say, (J^n_k)_xy = d(F^n_x)/dy (z^0_k) -- with z is the phase space coordinate)
  !>      are approximated by finite differences: here, by  DF / 2*h_y
  !>      with DF = F^n_x(z^0_k + (0,h_y,0,0)) - F^n_x(z^0_k - (0,h_y,0,0))
  !>    * then we invert that approximated forward Jacobian matrix
  !>
  !> Note: when computing finite differences in a periodic dimension (say x), one must be careful since two values of F_x
  !>    can be close in the periodic interval but distant by almost L_x in the (stored) [0,L_x[ representation.
  !>    To account for this (frequent) phenomenon we do the following:
  !>    when the difference DF (see example above) is larger than L_x/2, we make it smaller by adding +/- L_x to it.
  !>    Note that here we could very well be making the slope smaller than what it should actually be: indeed if the function
  !>    F^n_x is having a steep slope at z^0_k which adds one (or several) periods L_x to DF then our modification will
  !>    artificially lower the slope. But this is impossible to prevent in full generality (indeed: a steep slope crossing the
  !>    0 or L_x value will be lowered anyway in the [0,L_x[ representation) and should not be frequent (indeed it only happens
  !>    when F^n_x has high derivatives and the particle grid is coarse, which will lead to bad approximations anyhow).
  !  <<get_ltp_deformation_matrix>>
  subroutine get_ltp_deformation_matrix (       &
                        self,                &
                        k,                      &
                        mesh_period_x,          &
                        mesh_period_y,          &
                        h_particles_x,              &
                        h_particles_y,              &
                        h_particles_vx,             &
                        h_particles_vy,             &
                        x_k, y_k,               &
                        vx_k, vy_k,             &
                        d11,d12,d13,d14,        &
                        d21,d22,d23,d24,        &
                        d31,d32,d33,d34,        &
                        d41,d42,d43,d44         &
                        )

        class(sll_t_particle_group_2d2v_lbf),intent(inout) :: self
        sll_int32, intent(in) :: k

        sll_real64, intent(in)  :: mesh_period_x
        sll_real64, intent(in)  :: mesh_period_y

        sll_real64, intent(in)  :: h_particles_x
        sll_real64, intent(in)  :: h_particles_y
        sll_real64, intent(in)  :: h_particles_vx
        sll_real64, intent(in)  :: h_particles_vy

        sll_real64, intent(out) :: x_k, y_k         !< particle center in physical space
        sll_real64, intent(out) :: vx_k, vy_k       !< particle center in velocity space
        sll_real64, intent(out) :: d11,d12,d13,d14  !< coefs of matrix D (backward Jacobian)
        sll_real64, intent(out) :: d21,d22,d23,d24
        sll_real64, intent(out) :: d31,d32,d33,d34
        sll_real64, intent(out) :: d41,d42,d43,d44

        sll_int32   :: k_ngb
        sll_real64  :: x_k_left,  x_k_right
        sll_real64  :: y_k_left,  y_k_right
        sll_real64  :: vx_k_left, vx_k_right
        sll_real64  :: vy_k_left, vy_k_right
        sll_real64, dimension(3)  :: coords

        sll_real64  :: j11,j12,j13,j14   !< coefs of matrix J = D^-1 (forward Jacobian)
        sll_real64  :: j21,j22,j23,j24
        sll_real64  :: j31,j32,j33,j34
        sll_real64  :: j41,j42,j43,j44
        sll_real64  :: factor, det_J, inv_det_J

        ! for clarity:
        sll_int32, parameter :: dim_x  = 1
        sll_int32, parameter :: dim_y  = 2
        sll_int32, parameter :: dim_vx = 3
        sll_int32, parameter :: dim_vy = 4
        sll_int32, parameter :: left_ngb  = -1
        sll_int32, parameter :: right_ngb =  1

        logical domain_is_x_periodic
        logical domain_is_y_periodic

        domain_is_x_periodic = self%domain_is_periodic(1)
        domain_is_y_periodic = self%domain_is_periodic(2)

        coords(:) = 0.0_f64

        coords = self%get_x(k)
        x_k = coords(1)
        y_k = coords(2)

        coords = self%get_v(k)
        vx_k = coords(1)
        vy_k = coords(2)

        !> Compute the forward Jacobian matrix J_k

        !> ------   d/d_x terms
        factor = 1./(2*h_particles_x)

        k_ngb = self%get_neighbor_index(k, dim_x, right_ngb)
        if( k_ngb == k )then
           !> no right neighbor is available, use a non-centered finite difference
           factor = 2*factor
           x_k_right   = x_k
           y_k_right   = y_k
           vx_k_right  = vx_k
           vy_k_right  = vy_k
        else
            coords = self%get_x(k_ngb)
            x_k_right = coords(1)
            y_k_right = coords(2)
            coords = self%get_v(k_ngb)
            vx_k_right = coords(1)
            vy_k_right = coords(2)

            if( domain_is_x_periodic .and. x_k_right < x_k - 0.5*mesh_period_x ) x_k_right = x_k_right + mesh_period_x
            if( domain_is_x_periodic .and. x_k_right > x_k + 0.5*mesh_period_x ) x_k_right = x_k_right - mesh_period_x
            if( domain_is_y_periodic .and. y_k_right < y_k - 0.5*mesh_period_y ) y_k_right = y_k_right + mesh_period_y
            if( domain_is_y_periodic .and. y_k_right > y_k + 0.5*mesh_period_y ) y_k_right = y_k_right - mesh_period_y
        end if


        k_ngb = self%get_neighbor_index(k, dim_x, left_ngb)
        if( k_ngb == k )then
           !> no left neighbor is available, use a non-centered finite difference
           factor = 2*factor
           x_k_left   = x_k
           y_k_left   = y_k
           vx_k_left  = vx_k
           vy_k_left  = vy_k
        else
            coords = self%get_x(k_ngb)
            x_k_left = coords(1)
            y_k_left = coords(2)
            coords = self%get_v(k_ngb)
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

        !> ------   d/d_y terms
        factor = 1./(2*h_particles_y)

        k_ngb = self%get_neighbor_index(k, dim_y, right_ngb)
        if( k_ngb == k )then
           !> no right neighbor is available, use a non-centered finite difference
           factor = 2*factor
           x_k_right   = x_k
           y_k_right   = y_k
           vx_k_right  = vx_k
           vy_k_right  = vy_k
        else
            coords = self%get_x(k_ngb)
            x_k_right = coords(1)
            y_k_right = coords(2)
            coords = self%get_v(k_ngb)
            vx_k_right = coords(1)
            vy_k_right = coords(2)

            if( domain_is_x_periodic .and. x_k_right < x_k - 0.5*mesh_period_x ) x_k_right = x_k_right + mesh_period_x
            if( domain_is_x_periodic .and. x_k_right > x_k + 0.5*mesh_period_x ) x_k_right = x_k_right - mesh_period_x
            if( domain_is_y_periodic .and. y_k_right < y_k - 0.5*mesh_period_y ) y_k_right = y_k_right + mesh_period_y
            if( domain_is_y_periodic .and. y_k_right > y_k + 0.5*mesh_period_y ) y_k_right = y_k_right - mesh_period_y
        end if

        k_ngb = self%get_neighbor_index(k, dim_y, left_ngb)
        if( k_ngb == k )then
           !> no left neighbor is available, use a non-centered finite difference
           factor = 2*factor
           x_k_left   = x_k
           y_k_left   = y_k
           vx_k_left  = vx_k
           vy_k_left  = vy_k
        else
            coords = self%get_x(k_ngb)
            x_k_left = coords(1)
            y_k_left = coords(2)
            coords = self%get_v(k_ngb)
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


        !> ------   d/d_vx terms
        factor = 1./(2*h_particles_vx)

        k_ngb = self%get_neighbor_index(k, dim_vx, right_ngb)
        if( k_ngb == k )then
           !> no right neighbor is available, use a non-centered finite difference
           factor = 2*factor
           x_k_right   = x_k
           y_k_right   = y_k
           vx_k_right  = vx_k
           vy_k_right  = vy_k
        else
            coords = self%get_x(k_ngb)
            x_k_right = coords(1)
            y_k_right = coords(2)
            coords = self%get_v(k_ngb)
            vx_k_right = coords(1)
            vy_k_right = coords(2)

            if( domain_is_x_periodic .and. x_k_right < x_k - 0.5*mesh_period_x ) x_k_right = x_k_right + mesh_period_x
            if( domain_is_x_periodic .and. x_k_right > x_k + 0.5*mesh_period_x ) x_k_right = x_k_right - mesh_period_x
            if( domain_is_y_periodic .and. y_k_right < y_k - 0.5*mesh_period_y ) y_k_right = y_k_right + mesh_period_y
            if( domain_is_y_periodic .and. y_k_right > y_k + 0.5*mesh_period_y ) y_k_right = y_k_right - mesh_period_y
        end if

        k_ngb = self%get_neighbor_index(k, dim_vx, left_ngb)
        if( k_ngb == k )then
            !> no left neighbor is available, use a non-centered finite difference
            factor = 2*factor
            x_k_left   = x_k
            y_k_left   = y_k
            vx_k_left  = vx_k
            vy_k_left  = vy_k
        else
            coords = self%get_x(k_ngb)
            x_k_left = coords(1)
            y_k_left = coords(2)
            coords = self%get_v(k_ngb)
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


        !> ------   d/d_vy terms
        factor = 1./(2*h_particles_vy)

        k_ngb = self%get_neighbor_index(k, dim_vy, right_ngb)
        if( k_ngb == k )then
           !> no right neighbor is available, use a non-centered finite difference
           factor = 2*factor
           x_k_right   = x_k
           y_k_right   = y_k
           vx_k_right  = vx_k
           vy_k_right  = vy_k
        else
            coords = self%get_x(k_ngb)
            x_k_right = coords(1)
            y_k_right = coords(2)
            coords = self%get_v(k_ngb)
            vx_k_right = coords(1)
            vy_k_right = coords(2)

            if( domain_is_x_periodic .and. x_k_right < x_k - 0.5*mesh_period_x ) x_k_right = x_k_right + mesh_period_x
            if( domain_is_x_periodic .and. x_k_right > x_k + 0.5*mesh_period_x ) x_k_right = x_k_right - mesh_period_x
            if( domain_is_y_periodic .and. y_k_right < y_k - 0.5*mesh_period_y ) y_k_right = y_k_right + mesh_period_y
            if( domain_is_y_periodic .and. y_k_right > y_k + 0.5*mesh_period_y ) y_k_right = y_k_right - mesh_period_y
        end if

        k_ngb = self%get_neighbor_index(k, dim_vy, left_ngb)
        if( k_ngb == k )then
            !> no left neighbor is available, use a non-centered finite difference
            factor = 2*factor
            x_k_left   = x_k
            y_k_left   = y_k
            vx_k_left  = vx_k
            vy_k_left  = vy_k
        else
            coords = self%get_x(k_ngb)
            x_k_left = coords(1)
            y_k_left = coords(2)
            coords = self%get_v(k_ngb)
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


           !> Compute D_k the inverse of J_k

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

  !> periodic_correction -------------------------------------------------------------------------------------------------------
  !> puts the point (x,y) back into the computational domain if periodic in x or y (or both) otherwise, does nothing
  subroutine periodic_correction(p_group, x, y)
    class(sll_t_particle_group_2d2v_lbf),  intent(in)    :: p_group
    sll_real64, intent(inout) :: x
    sll_real64, intent(inout) :: y
    sll_real64 :: eta_max, eta_min

    if( p_group%domain_is_periodic(1) )then
      eta_min = p_group%lbf_grid%eta1_min
      eta_max = p_group%lbf_grid%eta1_max
      if( (x < eta_min) .or. (x >= eta_max) ) x = eta_min + modulo(x - eta_min, eta_max - eta_min)
    end if

    if( p_group%domain_is_periodic(2) )then
      eta_min = p_group%lbf_grid%eta2_min
      eta_max = p_group%lbf_grid%eta2_max
      if( (y < eta_min) .or. (y >= eta_max) ) y = eta_min + modulo(y - eta_min, eta_max - eta_min)
    end if
  end subroutine periodic_correction


end module sll_m_particle_group_2d2v_lbf
