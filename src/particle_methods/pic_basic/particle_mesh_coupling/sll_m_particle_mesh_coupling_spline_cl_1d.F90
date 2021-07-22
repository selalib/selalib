!> @ingroup particle_mesh_coupling
!> @author Benedikt Perse, IPP
!> @brief Kernel smoother for 2d with splines of arbitrary degree placed on a uniform mesh.
!> @details Spline with index i starts at point i
module sll_m_particle_mesh_coupling_spline_cl_1d

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_errors.h"

  use sll_m_low_level_bsplines, only: &
       sll_s_uniform_bsplines_eval_basis

  use sll_m_gauss_legendre_integration, only : &
       sll_f_gauss_legendre_points_and_weights

  use sll_m_particle_mesh_coupling_base_1d, only: &
       sll_p_collocation, &
       sll_c_particle_mesh_coupling_1d, &
       sll_p_galerkin

  use sll_m_particle_group_base, only: &
       sll_c_particle_group_base

  use sll_m_splines_pp, only: &
       sll_t_spline_pp_1d, &
       sll_s_spline_pp_horner_m_1d, &
       sll_s_spline_pp_horner_primitive_1d, &
       sll_f_spline_pp_horner_1d, &
       sll_s_spline_pp_free_1d, &
       sll_s_spline_pp_init_1d
  implicit none

  public :: &
       sll_t_particle_mesh_coupling_spline_cl_1d, &
       sll_s_new_particle_mesh_coupling_spline_cl_1d_ptr, &
       sll_s_new_particle_mesh_coupling_spline_cl_1d

  private
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !>  Spline kernel smoother in1d.
  type, extends(sll_c_particle_mesh_coupling_1d) :: sll_t_particle_mesh_coupling_spline_cl_1d
     type(sll_t_spline_pp_1d) :: spline_pp !< 1d pp-spline
     
     ! Information about the particles
     sll_int32  :: no_particles !< Number of particles of underlying PIC method (processor local)
     sll_int32  :: n_span !< Number of intervals where spline non zero (spline_degree + 1)
     sll_real64 :: scaling !< Scaling factor depending on whether Galerkin or collocation
     sll_int32  :: n_quad_points !< Number of quadrature points
     sll_real64, allocatable :: spline_val(:) !< scratch data for spline evaluation
     sll_real64, allocatable :: spline_val_more(:) !< more scratch data for spline evaluation
     ! For 1d quadrature 
     sll_int32 :: n_quad_points_line !< no. of quadrature points for the line integral
     sll_real64, allocatable :: quad_xw_line(:,:) !< quadrature weights and points for the line integral
     sll_real64, allocatable :: quad_xw(:,:) !< quadrature weights and points

   contains
     procedure :: add_charge => add_charge_single_spline_cl_1d !> Add charge of one particle
     procedure :: add_particle_mass => add_particle_mass_spline_cl_1d !> Add the contribution of one particle to the approximate mass matrix defined by the sum over all particles using the symmetry
     procedure :: add_particle_mass_full => add_particle_mass_spline_cl_1d_full !> Add the contribution of one particle to the approximate mass matrix defined by the sum over all particles
     procedure :: evaluate => evaluate_field_single_spline_cl_1d !> Evaluate spline function with given coefficients
     procedure :: evaluate_multiple => evaluate_multiple_spline_cl_1d !> Evaluate multiple spline functions with given coefficients
     procedure :: add_current_update_v => add_current_update_v_spline_cl_1d !> Add contribution of pne particle to the current density and update velocity
     procedure :: add_current => add_current_spline_cl_1d !> Add contribution of one particle to the current density (integrated over x)

     procedure :: init => init_spline_cl_1d !> Constructor
     procedure :: free => free_spline_cl_1d !> Destructor

  end type sll_t_particle_mesh_coupling_spline_cl_1d

contains


  !---------------------------------------------------------------------------!
  !> Add charge of one particle
  subroutine add_charge_single_spline_cl_1d(self, position, marker_charge, rho_dofs)
    class( sll_t_particle_mesh_coupling_spline_cl_1d ), intent(inout)   :: self !< kernel smoother object
    sll_real64,                               intent( in )    :: position(self%dim) !< Position of the particle
    sll_real64,                               intent( in )    :: marker_charge !< Particle weights time charge
    sll_real64,                               intent( inout ) :: rho_dofs(self%n_dofs) !< Coefficient vector of the charge distribution
    !local variables
    sll_int32 :: i
    sll_int32 :: box
    sll_real64 :: xi

    call convert_x_to_xbox( self, position, xi, box )
    call sll_s_uniform_bsplines_eval_basis_clamped( self%spline_pp, self%n_cells, self%spline_degree, xi, box, self%spline_val )

    do i = 1, self%n_span
       rho_dofs(box+i-1) = rho_dofs(box+i-1) +&
            (marker_charge * self%spline_val(i)* self%scaling)
    end do


  end subroutine add_charge_single_spline_cl_1d


  !---------------------------------------------------------------------------!
  !> Add the contribution of one particle to the approximate mass matrix
  subroutine add_particle_mass_spline_cl_1d(self, position, marker_charge, particle_mass)
    class( sll_t_particle_mesh_coupling_spline_cl_1d ), intent(inout)   :: self !< kernel smoother object
    sll_real64,                               intent( in )    :: position(self%dim) !< Position of the particle
    sll_real64,                               intent( in )    :: marker_charge !< Particle weights time charge
    sll_real64,                               intent( inout ) :: particle_mass(:,:) !< Coefficient vector of the charge distribution

    !local variables
    sll_int32 :: i, column
    sll_int32 :: box
    sll_real64 :: xi

    SLL_ASSERT( size(particle_mass,1) == self%n_span )
    SLL_ASSERT( size(particle_mass,2) == self%n_dofs )

    call convert_x_to_xbox( self, position, xi, box )

    call sll_s_uniform_bsplines_eval_basis_clamped( self%spline_pp, self%n_cells, self%spline_degree, xi, box, self%spline_val )

    do i = 1, self%n_span
       do column = i, self%n_span
          particle_mass(column-i+1, box+i-1) = particle_mass( column-i+1, box+i-1) + &
               marker_charge * self%spline_val(i) * self%spline_val(column)
          ! Note: No scaling since Galerkin function
       end do
    end do

  end subroutine add_particle_mass_spline_cl_1d


  !----------------------------------------------------------------------!
  !> Add the contribution of one particle to the approximate mass matrix
  subroutine add_particle_mass_spline_cl_1d_full(self, position, marker_charge, particle_mass)
    class( sll_t_particle_mesh_coupling_spline_cl_1d ), intent(inout)   :: self !< kernel smoother object
    sll_real64,                               intent( in )    :: position(self%dim) !< Position of the particle
    sll_real64,                               intent( in )    :: marker_charge !< Particle weights time charge
    sll_real64,                               intent( inout ) :: particle_mass(:,:) !< Coefficient vector of the charge distribution
    !local variables
    sll_int32 :: i, column, ind
    sll_int32 :: box
    sll_real64 :: xi

    SLL_ASSERT( size(particle_mass,1) == 2*self%spline_degree+1 )
    SLL_ASSERT( size(particle_mass,2) == self%n_dofs )

    call convert_x_to_xbox( self, position, xi, box )

    call sll_s_uniform_bsplines_eval_basis_clamped( self%spline_pp, self%n_cells, self%spline_degree, xi, box, self%spline_val )

    do i = 1, self%n_span
       ind=1+(self%n_span-i)
       do column = 1, self%n_span
          particle_mass(ind, box+i-1) = particle_mass( ind, box+i-1) + &
               marker_charge * self%spline_val(i) * self%spline_val(column)
          ind = ind+1
       end do
    end do

  end subroutine add_particle_mass_spline_cl_1d_full


  !---------------------------------------------------------------------------!
  !> Evaluate field at at position \a position
  subroutine evaluate_field_single_spline_cl_1d(self, position, field_dofs, field_value)
    class (sll_t_particle_mesh_coupling_spline_cl_1d), intent( inout ) :: self !< Kernel smoother object 
    sll_real64,                              intent( in )    :: position(self%dim) !< Position of the particle
    sll_real64,                              intent( in )    :: field_dofs(self%n_dofs) !< Coefficient vector for the field DoFs
    sll_real64,                              intent( out )   :: field_value !< Value(s) of the electric fields at given position
    !local variables
    sll_int32 :: i
    sll_int32 :: box
    sll_real64 :: xi

    call convert_x_to_xbox( self, position, xi, box )

    call sll_s_uniform_bsplines_eval_basis_clamped( self%spline_pp, self%n_cells, self%spline_degree, xi, box, self%spline_val )


    field_value = 0.0_f64
    do i = 1, self%n_span
       field_value = field_value + &
            field_dofs(box+i-1) *  &
            self%spline_val(i)
    end do

  end subroutine evaluate_field_single_spline_cl_1d


  !---------------------------------------------------------------------------!
  !> Evaluate several fields at position \a position
  subroutine evaluate_multiple_spline_cl_1d(self, position, components, field_dofs, field_value)
    class (sll_t_particle_mesh_coupling_spline_cl_1d), intent( inout ) :: self !< Kernel smoother object 
    sll_real64,                              intent( in )    :: position(self%dim) !< Position of the particle
    sll_int32,                               intent(in)      :: components(:) !< Components of field_dofs that shall be updated
    sll_real64,                              intent( in )    :: field_dofs(:,:) !< Coefficient vector for the field DoFs
    sll_real64,                              intent(out)     :: field_value(:) !< Value(s) of the electric fields at given position
    !local variables
    sll_int32 :: i
    sll_int32 :: box
    sll_real64 :: xi

    SLL_ASSERT( size(field_dofs,1) == self%n_dofs )
    SLL_ASSERT( size(field_dofs,2) == size(field_value) )

    call convert_x_to_xbox( self, position, xi, box )

    call sll_s_uniform_bsplines_eval_basis_clamped( self%spline_pp, self%n_cells, self%spline_degree, xi, box, self%spline_val )


    field_value = 0.0_f64
    do i = 1, self%n_span

       field_value = field_value + &
            field_dofs(box+i-1,components) *  &
            self%spline_val(i)
    end do

  end subroutine evaluate_multiple_spline_cl_1d


  !---------------------------------------------------------------------------!
  !> Add current with integration over x
  subroutine add_current_spline_cl_1d( self, position_old, position_new, marker_charge, j_dofs )
    class( sll_t_particle_mesh_coupling_spline_cl_1d ), intent(inout)   :: self !< kernel smoother object
    sll_real64,                               intent( in )    :: position_old(self%dim) !< Position of the particle
    sll_real64,                               intent( in )    :: position_new(self%dim) !< Position of the particle
    sll_real64,                               intent( in )    :: marker_charge !< Particle weights time charge
    sll_real64,                               intent( inout ) :: j_dofs(self%n_dofs) !< Coefficient vector of the current density
    !local variables
    sll_real64 :: xold, xnew, xnewtilde, xbox
    sll_int32 :: boxold, boxnew, sigma_counter, boxdiff, increment, box
    sll_real64 :: sigma_r, sigma_l, sigma, sigma_next, weight
    sll_int32 :: q, bl
    sll_real64, allocatable :: vals(:) 

    call convert_x_to_xbox( self, position_old, xold, boxold )
    call convert_x_to_xbox( self, position_new, xnew, boxnew )

    ! Now we need to compute the normalized 1d line along which to integrate:
    boxdiff = boxnew-boxold
    xnewtilde = xnew + real(boxdiff,f64)
    sigma_l = 0.0_f64
    box = boxold ! We start in the old box
    sigma_counter = 0

    if (boxdiff > 0 ) then
       allocate( vals(boxdiff+1) )
       do bl = 1, boxdiff
          vals(bl) = (real(bl,f64)  - xold)/(xnewtilde-xold)
       end do
       vals(boxdiff+1) = 1.0_f64
       sigma_next = vals(1)
       increment = 1
    elseif (boxdiff < 0 ) then
       allocate ( vals(-boxdiff+1) )
       do bl = boxdiff+1, 0
          vals(-bl+1) = (real(bl,f64)  - xold)/(xnewtilde-xold)
       end do
       vals(-boxdiff+1) = 1.0_f64
       sigma_next = vals(1)
       increment = -1
    else
       sigma_next = 1.0_f64
       increment = 0
    end if


    sigma_r = 0.0_f64
    do while ( sigma_r < 1.0_f64 )
       sigma_r = sigma_next

       do q = 1, self%n_quad_points_line
          sigma = sigma_l + (sigma_r-sigma_l) * self%quad_xw_line(1,q)
          xbox = xold* (1.0_f64 - sigma) + xnewtilde * sigma  - real(sigma_counter*increment, f64)
          if (xbox > 1.0_f64 ) then
             xbox = min(xbox, 1._f64)
          elseif (xbox < 0.0_f64 ) then
             xbox = max(xbox, 0._f64)
          end if

          weight = self%quad_xw_line(2,q)*(sigma_r-sigma_l)
          call integrate_spline_cl_1d(self, box, xbox, marker_charge*weight, j_dofs )
       end do
       if (sigma_r < 1.0_f64 ) then
          ! Update the
          sigma_counter = sigma_counter+1
          sigma_next = vals(sigma_counter+1)
          box = box + increment
          sigma_l = sigma_r
       end if
    end do

  end subroutine add_current_spline_cl_1d


  !> Helper function for \a add_current
  subroutine integrate_spline_cl_1d( self, box, xbox, weight, j_dofs )
    class( sll_t_particle_mesh_coupling_spline_cl_1d ), intent(inout)   :: self !< kernel smoother object
    sll_int32,  intent( in    ) :: box
    sll_real64, intent( in    ) :: xbox
    sll_real64, intent( in    ) :: weight
    sll_real64, intent( inout ) :: j_dofs(:)
    !local variables
    sll_int32 :: i

    call sll_s_uniform_bsplines_eval_basis_clamped( self%spline_pp, self%n_cells, self%spline_degree, xbox, box, self%spline_val  )

    do i = 1, self%spline_degree+1
       j_dofs(box+i-1) = j_dofs(box+i-1) + &
            weight * self%spline_val(i) 
    end do

  end subroutine integrate_spline_cl_1d


  !---------------------------------------------------------------------------!
  !> Add current for one particle and update v (according to H_p1 part in Hamiltonian splitting)
  subroutine add_current_update_v_spline_cl_1d (self, position_old, position_new, marker_charge, qoverm, bfield_dofs, vi, j_dofs)
    class(sll_t_particle_mesh_coupling_spline_cl_1d), intent(inout) :: self !< kernel smoother object
    sll_real64, intent(in)    :: position_old(self%dim) !< Position at time t
    sll_real64, intent(in)    :: position_new(self%dim) !< Position at time t + \Delta t
    sll_real64, intent(in)    :: marker_charge !< Particle weight time charge
    sll_real64, intent(in)    :: qoverm !< charge to mass ration
    sll_real64, intent(in)    :: bfield_dofs(:) !< Coefficient of B-field expansion
    sll_real64, intent(inout) :: vi(:) !< Velocity of the particles
    sll_real64, intent(inout) :: j_dofs(:) !< Coefficients of current expansion
    ! local variables
    sll_int32  :: index_old, index_new, ind
    sll_real64 :: r_old, r_new

    call convert_x_to_xbox( self, position_old, r_old, index_old )
    call convert_x_to_xbox( self, position_new, r_new, index_new )

    if (index_old == index_new) then
       if (r_old < r_new) then
          call update_jv( self, r_old, r_new, index_old, marker_charge, &
               qoverm, 1.0_f64, vi(2), j_dofs, bfield_dofs)
       else
          call update_jv( self, r_new, r_old, index_old, marker_charge, qoverm, &
               -1.0_f64, vi(2), j_dofs, bfield_dofs)
       end if
    elseif (index_old < index_new) then
       call update_jv ( self, r_old, 1.0_f64, index_old, marker_charge, &
            qoverm, 1.0_f64, vi(2), j_dofs, bfield_dofs)
       call update_jv ( self, 0.0_f64, r_new, index_new, marker_charge, &
            qoverm, 1.0_f64, vi(2), j_dofs, bfield_dofs)
       do ind = index_old+1, index_new-1
          call update_jv ( self, 0.0_f64, 1.0_f64, ind, marker_charge, &
               qoverm, 1.0_f64, vi(2), j_dofs, bfield_dofs)
       end do
    else
       call update_jv ( self, r_new, 1.0_f64, index_new, marker_charge, qoverm, &
            -1.0_f64, vi(2), j_dofs, bfield_dofs)
       call update_jv ( self, 0.0_f64, r_old, index_old, marker_charge, qoverm, &
            -1.0_f64, vi(2), j_dofs, bfield_dofs)
       do ind = index_new+1, index_old-1
          call update_jv ( self, 0.0_f64, 1.0_f64, ind, marker_charge, qoverm, &
               -1.0_f64, vi(2), j_dofs, bfield_dofs)
       end do
    end if

  end subroutine add_current_update_v_spline_cl_1d

  !> Helper function for \a add_current_update_v.
  subroutine update_jv(self, lower, upper, index, marker_charge, qoverm, sign, vi, j_dofs, bfield_dofs)
    class(sll_t_particle_mesh_coupling_spline_cl_1d), intent(inout) :: self !< time splitting object 
    sll_real64,                             intent(in)    :: lower
    sll_real64,                             intent(in)    :: upper
    sll_int32,                              intent(in)    :: index
    sll_real64,                             intent(in)    :: marker_charge
    sll_real64,                             intent(in)    :: qoverm
    sll_real64,                             intent(in)    :: sign
    sll_real64,                             intent(inout) :: vi
    sll_real64,                             intent(in)    :: bfield_dofs(:)
    sll_real64,                             intent(inout) :: j_dofs(:)
    !Local variables
    sll_int32  :: ind, i_grid, j
    sll_real64 :: c1, c2

    c1 =  0.5_f64*(upper-lower)
    c2 =  0.5_f64*(upper+lower)

    call sll_s_uniform_bsplines_eval_basis_clamped( self%spline_pp, self%n_cells, self%spline_degree, c1*self%quad_xw(1,1)+c2, index, self%spline_val )
    self%spline_val = self%spline_val * (self%quad_xw(2,1)*c1)
    do j=2,self%n_quad_points
       call sll_s_uniform_bsplines_eval_basis_clamped( self%spline_pp, self%n_cells, self%spline_degree, c1*self%quad_xw(1,j)+c2, index, self%spline_val_more )
       self%spline_val = self%spline_val + self%spline_val_more * (self%quad_xw(2,j)*c1)
    end do
    self%spline_val = self%spline_val * (sign*self%delta_x)

    ind = 1
    do i_grid = index, index + self%spline_degree
       j_dofs(i_grid) = j_dofs(i_grid) + &
            (marker_charge*self%spline_val(ind)* self%scaling)
       vi = vi - qoverm* self%spline_val(ind)*bfield_dofs(i_grid)
       ind = ind + 1
    end do

  end subroutine update_jv


  !< Helper function to convert x to xbox
  subroutine convert_x_to_xbox( self, position, xi, box )
    class( sll_t_particle_mesh_coupling_spline_cl_1d ), intent(inout)   :: self !< kernel smoother object
    sll_real64,                               intent( in )    :: position(self%dim) !< Position of the particle
    sll_real64,                               intent( out )    :: xi !< Position of the particle
    sll_int32,                                intent( out )    :: box !< Position of the particle

    xi = (position(1) - self%domain(1)) /self%delta_x
    box = floor( xi ) + 1
    xi = xi - real(box-1, f64)

    if( box == self%n_cells + 1 ) then
       if( xi == 0._f64)then 
          xi = 1._f64
          box = self%n_cells
       else
          print*, 'box,x, xbox', box, position(1), xi
          SLL_ERROR('convert_x_to_xbox', 'box1 value to high' ) 
       end if
    else if( box > self%n_cells + 1 ) then
       print*, 'box,x, xbox', box, position(1), xi
       SLL_ERROR('convert_x_to_xbox', 'box1 value to high' ) 
    else if( box <= 0 ) then
       print*, 'box,x, xbox', box, position(1), xi
       SLL_ERROR('convert_x_to_xbox', 'box1 value to low' ) 
    end if

  end subroutine convert_x_to_xbox

  !> Helper function to evaluate uniform clamped basis splines
  subroutine sll_s_uniform_bsplines_eval_basis_clamped( spline, n_cells, degree, xi, box, spline_val )
    type(sll_t_spline_pp_1d), intent( in    ) :: spline
    sll_int32,                intent( in    ) :: n_cells
    sll_int32,                intent( in    ) :: degree
    sll_real64,               intent( in    ) :: xi
    sll_int32,                intent( in    ) :: box
    sll_real64,               intent(   out ) :: spline_val(:)
    !local variables
    sll_int32 :: i

    if(box >= 1 .and. box <= degree-1)then
       do i=1, degree+1
          spline_val(i) = sll_f_spline_pp_horner_1d( degree, spline%poly_coeffs_boundary_left(:,:,box), xi, i)
       end do
    else if (box >= degree .and. box <= n_cells-degree+1)then
       call sll_s_uniform_bsplines_eval_basis( degree, xi, spline_val )
    else if(box >= n_cells-degree+2 .and. box <= n_cells)then
       do i=1, degree+1
          spline_val(i) = sll_f_spline_pp_horner_1d( degree, spline%poly_coeffs_boundary_right(:,:,box-1-n_cells+degree), xi, i)
       end do
    else
       print*, 'box, xbox', box, xi
       SLL_ERROR( "uniform_bsplines_eval_basis_clamped", "box out of range" )
    end if

  end subroutine sll_s_uniform_bsplines_eval_basis_clamped


  !> Initializer
  subroutine init_spline_cl_1d( self, domain, n_cells, no_particles, spline_degree, smoothing_type, boundary )
    class( sll_t_particle_mesh_coupling_spline_cl_1d), intent(out)  :: self !< kernel smoother object
    sll_real64,                              intent(in) :: domain(2) !< x_min and x_max of the domain
    sll_int32,                               intent(in) :: n_cells !< number of DoFs (spline coefficients)
    sll_int32,                               intent(in) :: no_particles !< number of particles
    sll_int32,                               intent(in) :: spline_degree !< Degree of smoothing kernel spline
    sll_int32,                               intent(in) :: smoothing_type !< Define if Galerkin or collocation smoothing for right scaling in accumulation routines
    sll_int32,                               intent(in) :: boundary !< spline boundary conditions
    !local variables
    sll_int32 :: ierr


    self%dim = 1

    ! Store grid information
    self%domain = domain
    self%n_cells = n_cells
    self%n_dofs = n_cells+spline_degree
    self%delta_x = (self%domain(2)-self%domain(1))/real(n_cells, f64)

    ! Store basis function information
    self%no_particles = no_particles

    ! Initialize information on the spline
    self%spline_degree = spline_degree
    self%n_span = spline_degree + 1

    ! Initialize information on smoothing type
    if (smoothing_type == sll_p_collocation) then
       self%scaling = 1.0_f64/self%delta_x
    elseif (smoothing_type == sll_p_galerkin) then
       self%scaling = 1.0_f64
    else
       print*, 'Smoothing Type ', smoothing_type, ' not implemented for kernel_smoother_spline_1d.'
    end if

    self%n_quad_points = (self%spline_degree+2)/2

    ALLOCATE( self%spline_val(self%n_span), stat = ierr)
    SLL_ASSERT( ierr == 0 )
    ALLOCATE( self%spline_val_more(self%n_span), stat = ierr )
    SLL_ASSERT( ierr == 0 )
    ALLOCATE( self%quad_xw(2,self%n_quad_points), stat = ierr )
    SLL_ASSERT( ierr == 0 )

    call sll_s_spline_pp_init_1d( self%spline_pp, spline_degree, n_cells, boundary)

    self%n_quad_points_line = (self%spline_degree+2)/2

    allocate( self%quad_xw_line(2,self%n_quad_points_line) )
    ! normalized Gauss Legendre points and weights
    self%quad_xw_line = sll_f_gauss_legendre_points_and_weights(self%n_quad_points_line)

    self%quad_xw_line(1,:) = 0.5_f64*(self%quad_xw_line(1,:)+1.0_f64)
    self%quad_xw_line(2,:) = 0.5_f64*self%quad_xw_line(2,:)


  end subroutine init_spline_cl_1d


  !> Finalizer
  subroutine free_spline_cl_1d(self)
    class (sll_t_particle_mesh_coupling_spline_cl_1d), intent( inout ) :: self !< Kernel smoother object 

    deallocate(self%spline_val)
    call sll_s_spline_pp_free_1d(self%spline_pp)

  end subroutine free_spline_cl_1d


  !-------------------------------------------------------------------------------------------
  !< Constructor for abstract type
  subroutine sll_s_new_particle_mesh_coupling_spline_cl_1d_ptr(smoother, domain, n_cells, no_particles, spline_degree, smoothing_type, boundary)
    class( sll_c_particle_mesh_coupling_1d), pointer, intent(out) :: smoother !< kernel smoother object
    sll_int32,                              intent(in)  :: n_cells !< number of DoFs (spline coefficients)
    sll_real64,                             intent(in)  :: domain(2) !< x_min and x_max of the domain
    sll_int32,                              intent(in)  :: no_particles !< number of particles
    sll_int32,                              intent(in)  :: spline_degree !< Degree of smoothing kernel spline
    sll_int32,                              intent(in)  :: smoothing_type !< Define if Galerkin or collocation smoothing for right scaling in accumulation routines
    sll_int32,                               intent(in) :: boundary !< spline boundary conditions
    !local variables
    sll_int32 :: ierr


    SLL_ALLOCATE( sll_t_particle_mesh_coupling_spline_cl_1d :: smoother , ierr)
    SLL_ASSERT( ierr == 0)

    select type( smoother )
    type is ( sll_t_particle_mesh_coupling_spline_cl_1d )
       call smoother%init( domain, n_cells, no_particles, spline_degree, smoothing_type, boundary )
    end select

  end subroutine sll_s_new_particle_mesh_coupling_spline_cl_1d_ptr


  !----------------------------------------------------------------------------------
  !< Constructor for abstract type
  subroutine sll_s_new_particle_mesh_coupling_spline_cl_1d(smoother, domain, n_cells, no_particles, spline_degree, smoothing_type, boundary)
    class( sll_c_particle_mesh_coupling_1d), allocatable, intent(out):: smoother !< kernel smoother object
    sll_int32,                                  intent(in) :: n_cells !< number of DoFs (spline coefficients)
    sll_real64,                                 intent(in) :: domain(2) !< x_min and x_max of the domain
    sll_int32,                                  intent(in) :: no_particles !< number of particles
    sll_int32,                                  intent(in) :: spline_degree !< Degree of smoothing kernel spline
    sll_int32,                                  intent(in) :: smoothing_type !< Define if Galerkin or collocation smoothing for right scaling in accumulation routines
    sll_int32,                                  intent(in) :: boundary !< spline boundary conditions
    !local variables
    sll_int32 :: ierr


    allocate( sll_t_particle_mesh_coupling_spline_cl_1d :: smoother , stat=ierr)
    SLL_ASSERT( ierr == 0)

    select type( smoother )
    type is ( sll_t_particle_mesh_coupling_spline_cl_1d )
       call smoother%init( domain, n_cells, no_particles, spline_degree, smoothing_type, boundary )
    end select

  end subroutine sll_s_new_particle_mesh_coupling_spline_cl_1d


end module sll_m_particle_mesh_coupling_spline_cl_1d
