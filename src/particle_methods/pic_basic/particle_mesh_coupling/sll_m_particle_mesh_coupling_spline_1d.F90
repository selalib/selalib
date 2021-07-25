!> @ingroup particle_mesh_coupling
!> @author Katharina Kormann, IPP
!> @brief Kernel smoother for 1d with splines of arbitrary degree placed on a uniform mesh.
!> @details Spline with index i starts at point i
module sll_m_particle_mesh_coupling_spline_1d

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_errors.h"

  use sll_m_gauss_legendre_integration, only : &
       sll_f_gauss_legendre_points_and_weights

  use sll_m_low_level_bsplines, only: &
       sll_s_uniform_bsplines_eval_basis, &
       sll_s_uniform_bsplines_eval_deriv

  use sll_m_particle_mesh_coupling_base_1d, only: &
       sll_p_collocation, &
       sll_c_particle_mesh_coupling_1d, &
       sll_p_galerkin

  use sll_m_particle_group_base, only: &
       sll_c_particle_group_base

  use sll_m_splines_pp, only: &
       sll_t_spline_pp_1d, &
       sll_s_spline_pp_init_1d, &
       sll_s_spline_pp_free_1d, &
       sll_f_spline_pp_horner_1d, &
       sll_s_spline_pp_horner_primitive_1d

  implicit none

  public :: &
       sll_t_particle_mesh_coupling_spline_1d, &
       sll_s_new_particle_mesh_coupling_spline_1d_ptr, &
       sll_s_new_particle_mesh_coupling_spline_1d

  private
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !>  Spline kernel smoother in1d.
  type, extends(sll_c_particle_mesh_coupling_1d) :: sll_t_particle_mesh_coupling_spline_1d
     type(sll_t_spline_pp_1d) :: spline_pp !< 1d pp-spline
     
     ! Information about the particles
     sll_int32  :: no_particles !< Number of particles of underlying PIC method (processor local)
     sll_int32  :: n_span !< Number of intervals where spline non zero (spline_degree + 1)
     sll_real64 :: scaling !< Scaling factor depending on whether Galerkin or collocation
     sll_int32  :: n_quad_points !< Number of quadrature points
     sll_int32  :: n_quadp1_points !< Number of quadrature points for non-linear disgradE and smoothing (add_charge and evaluate)
     sll_int32  :: n_quadp2_points !< Number of quadrature points for degree increased by 2
     sll_int32  :: n_quads2_points !< Number of quadrature points for  smoothing (add_current)

     sll_real64, allocatable :: spline_val(:) !< scratch data for spline evaluation
     sll_real64, allocatable :: spline_pol_val(:) !< scratch data for spline evaluation
     sll_real64, allocatable :: spline_val_more(:) !< more scratch data for spline evaluation
     sll_real64, allocatable :: quad_xw(:,:) !< quadrature weights and points
     sll_real64, allocatable :: quadp1_xw(:,:) !< quadrature weights and points
     sll_real64, allocatable :: quadp2_xw(:,:) !< quadrature weights and points
     sll_real64, allocatable :: quads1_xw(:,:) !< quadrature weights and points for smoothing (evaluate and add_charge)
     sll_real64, allocatable :: quads2_xw(:,:) !< quadrature weights and points for smoothing (add_current)

     sll_int32 :: n_quad_points_line !< no. of quadrature points for the line integral
     sll_real64, allocatable :: quad_xw_line(:,:) !< quadrature weights and points for the line integral

   contains
     procedure :: add_charge => add_charge_single_spline_1d !> Add charge of one particle
     procedure :: add_particle_mass => add_particle_mass_spline_1d !> Add the contribution of one particle to the approximate mass matrix defined by the sum over all particles using the symmetry
     procedure :: add_particle_mass_full => add_particle_mass_spline_1d_full !> Add the contribution of one particle to the approximate mass matrix defined by the sum over all particles
     procedure :: evaluate => evaluate_field_single_spline_1d !> Evaluate spline function with given coefficients
     procedure :: evaluate_multiple => evaluate_multiple_spline_1d !> Evaluate multiple spline functions with given coefficients
     procedure :: add_current_update_v => add_current_update_v_spline_pp_1d !> Add contribution of pne particle to the current density and update velocity
     procedure :: add_current => add_current_spline_1d !> Add contribution of one particle to the current density (integrated over x)

     procedure :: init => init_spline_1d !> Constructor
     procedure :: free => free_spline_1d !> Destructor

     procedure :: add_current_evaluate => add_current_evaluate_spline_1d !> Combines add_current and evaluate_int
     procedure :: evaluate_int => evaluate_int_spline_1d !> Evaluates the integral int_{poisition_old}^{position_new} field(x) d x
     procedure :: evaluate_int_quad => evaluate_int_quad_spline_1d !> For explicit subcycling propagator, evaluate_int based on quadrature formula for integration
     procedure :: evaluate_int_linear_quad => evaluate_int_linear_quad_spline_1d !> Evaluates the integrals with tau function for implicit subcycling

     procedure :: evaluate_simple => evaluate_field_single_spline_1d !> The pointer points to the evaluation functions "evaluate". This pointer is provided, since the derive type (smooth version) needs to access this evaluation but not as its own evaluation function.

     procedure :: evaluate_int_subnewton => evaluate_int_spline_1d_subnewton !> For explicit subcycling propagator, evaluates the integral needed for Newtons's method.
     procedure :: evaluate_int_linear_quad_subnewton => evaluate_int_linear_quad_subnewton_spline_1d !> For implicit subcycling propagator
     procedure :: add_current_split => add_current_split_spline_1d !> Add current version with tau and (1-tau) which is needed for the implicit subcycling scheme
     procedure :: add_charge_int => add_charge_int_spline_1d !> For explicit version of the subcycling propagator (to evaluate the charge as it comes out integrated from the formulation)

     procedure :: update_jv !> helper function to compute the integral of j using Gauss quadrature (needs to be provided for the derived type with smoothing)

  end type sll_t_particle_mesh_coupling_spline_1d

contains


  !---------------------------------------------------------------------------!
  !> Add charge of one particle
  subroutine add_charge_single_spline_1d(self, position, marker_charge, rho_dofs)
    class( sll_t_particle_mesh_coupling_spline_1d ), intent(inout)   :: self !< kernel smoother object
    sll_real64,                               intent( in )    :: position(self%dim) !< Position of the particle
    sll_real64,                               intent( in )    :: marker_charge !< Particle weights time charge
    sll_real64,                               intent( inout ) :: rho_dofs(self%n_cells) !< Coefficient vector of the charge distribution
    !local variables
    sll_int32 :: i1
    sll_int32 :: index1d, box
    sll_real64 :: xi 

    call convert_x_to_xbox( self, position, xi, box )
    box = box - self%spline_degree
    call sll_s_uniform_bsplines_eval_basis(self%spline_degree, xi, self%spline_val)
    !alternative version with horner and pp-form
    !call sll_s_spline_pp_horner_m_1d(self%spline_pp, self%spline_val, self%spline_degree, xi)

    do i1 = 1, self%n_span
       index1d = modulo(box+i1-2,self%n_cells)+1
       rho_dofs(index1d) = rho_dofs(index1d) +&
            (marker_charge * self%spline_val(i1)* self%scaling)
    end do

  end subroutine add_charge_single_spline_1d


  !---------------------------------------------------------------------------!
  !> Add the contribution of one particle to the approximate mass matrix
  subroutine add_particle_mass_spline_1d(self, position, marker_charge, particle_mass)
    class( sll_t_particle_mesh_coupling_spline_1d ), intent(inout)   :: self !< kernel smoother object
    sll_real64,                               intent( in )    :: position(self%dim) !< Position of the particle
    sll_real64,                               intent( in )    :: marker_charge !< Particle weights time charge
    sll_real64,                               intent( inout ) :: particle_mass(:,:) !< Coefficient vector of the charge distribution

    !local variables
    sll_int32 :: i1, column
    sll_int32 :: index1d, box
    sll_real64 :: xi

    SLL_ASSERT( size(particle_mass,1) == self%n_span )
    SLL_ASSERT( size(particle_mass,2) == self%n_cells )

    call convert_x_to_xbox( self, position, xi, box )
    box = box - self%spline_degree
    call sll_s_uniform_bsplines_eval_basis(self%spline_degree, xi, self%spline_val)
    !alternative version with horner and pp-form
    !call sll_s_spline_pp_horner_m_1d(self%spline_pp, self%spline_val, self%spline_degree, xi)

    do i1 = 1, self%n_span
       index1d = modulo(box+i1-2,self%n_cells)+1
       do column = i1, self%n_span
          particle_mass(column-i1+1, index1d) = particle_mass( column-i1+1, index1d) + &
               marker_charge * self%spline_val(i1) * self%spline_val(column)
          ! Note: No scaling since Galerkin function
       end do
    end do

  end subroutine add_particle_mass_spline_1d

  !----------------------------------------------------------------------!
  !> Add the contribution of one particle to the approximate mass matrix
  subroutine add_particle_mass_spline_1d_full(self, position, marker_charge, particle_mass)
    class( sll_t_particle_mesh_coupling_spline_1d ), intent(inout)   :: self !< kernel smoother object
    sll_real64,                               intent( in )    :: position(self%dim) !< Position of the particle
    sll_real64,                               intent( in )    :: marker_charge !< Particle weights time charge
    sll_real64,                               intent( inout ) :: particle_mass(:,:) !< Coefficient vector of the charge distribution

    !local variables
    sll_int32 :: i1, column, ind
    sll_int32 :: index1d, box
    sll_real64 :: xi

    SLL_ASSERT( size(particle_mass,1) == 2*self%spline_degree+1 )
    SLL_ASSERT( size(particle_mass,2) == self%n_cells )

    call convert_x_to_xbox( self, position, xi, box )
    box = box - self%spline_degree
    call sll_s_uniform_bsplines_eval_basis(self%spline_degree, xi, self%spline_val)
    !alternative version with horner and pp-form
    !call sll_s_spline_pp_horner_m_1d(self%spline_pp, self%spline_val, self%spline_degree, xi)
    do i1 = 1, self%n_span
       index1d = modulo(box+i1-2,self%n_cells)+1
       ind=1+(self%n_span-i1)
       do column = 1, self%n_span
          particle_mass(ind, index1d) = particle_mass( ind, index1d) + &
               marker_charge * self%spline_val(i1) * self%spline_val(column)
          ind = ind+1
       end do
    end do

  end subroutine add_particle_mass_spline_1d_full


  !---------------------------------------------------------------------------!
  !> Evaluate field at at position \a position
  subroutine evaluate_field_single_spline_1d(self, position, field_dofs, field_value)
    class (sll_t_particle_mesh_coupling_spline_1d), intent( inout ) :: self !< Kernel smoother object 
    sll_real64,                              intent( in )    :: position(self%dim) !< Position of the particle
    sll_real64,                              intent( in )    :: field_dofs(self%n_cells) !< Coefficient vector for the field DoFs
    sll_real64,                              intent( out )   :: field_value !< Value(s) of the electric fields at given position

    !local variables
    sll_int32 :: i1
    sll_int32 :: index1d, box
    sll_real64 :: xi

    call convert_x_to_xbox( self, position, xi, box )
    box = box - self%spline_degree
    ! Version based on the uniform bspline functions
    call sll_s_uniform_bsplines_eval_basis(self%spline_degree, xi, self%spline_val)
    !alternative version with horner and pp-form
    !call sll_s_spline_pp_horner_m_1d(self%spline_pp, self%spline_val, self%spline_degree, xi)

    field_value = 0.0_f64
    do i1 = 1, self%n_span
       index1d = modulo(box+i1-2, self%n_cells)+1
       field_value = field_value + &
            field_dofs(index1d) *  &
            self%spline_val(i1)
    end do

  end subroutine evaluate_field_single_spline_1d


  !---------------------------------------------------------------------------!
  !> Evaluate several fields at position \a position
  subroutine evaluate_multiple_spline_1d(self, position, components, field_dofs, field_value)
    class (sll_t_particle_mesh_coupling_spline_1d), intent( inout ) :: self !< Kernel smoother object 
    sll_real64,                              intent( in )    :: position(self%dim) !< Position of the particle
    sll_int32,                               intent(in)      :: components(:) !< Components of field_dofs that shall be updated
    sll_real64,                              intent( in )    :: field_dofs(:,:) !< Coefficient vector for the field DoFs
    sll_real64,                              intent(out)     :: field_value(:) !< Value(s) of the electric fields at given position

    !local variables
    sll_int32 :: i1
    sll_int32 :: index1d, box
    sll_real64 :: xi

    SLL_ASSERT( size(field_dofs,1) == self%n_cells )
    SLL_ASSERT( size(field_dofs,2) == size(field_value) )

    call convert_x_to_xbox( self, position, xi, box )

    box = box - self%spline_degree    
    call sll_s_uniform_bsplines_eval_basis(self%spline_degree, xi, self%spline_val)
    !alternative version with horner and pp-form
    !call sll_s_spline_pp_horner_m_1d(self%spline_pp, self%spline_val, self%spline_degree, xi)

    field_value = 0.0_f64
    do i1 = 1, self%n_span
       index1d = modulo(box+i1-2, self%n_cells)+1
       field_value = field_value + &
            field_dofs(index1d,components) *  &
            self%spline_val(i1)
    end do

  end subroutine evaluate_multiple_spline_1d


  !> Add current with integration over x
  subroutine add_current_spline_1d( self, position_old, position_new, marker_charge, j_dofs )
    class( sll_t_particle_mesh_coupling_spline_1d ), intent(inout)   :: self !< kernel smoother object
    sll_real64,                               intent( in )    :: position_old(self%dim) !< Position of the particle
    sll_real64,                               intent( in )    :: position_new(self%dim) !< Position of the particle
    sll_real64,                               intent( in )    :: marker_charge !< Particle weights time charge
    sll_real64,                               intent( inout ) :: j_dofs(self%n_cells) !< Coefficient vector of the current density

    sll_real64 :: dx_new, dx_old
    sll_int32 :: box_new, box_old, i, index, j, ind

    call convert_x_to_xbox( self, position_new, dx_new, box_new )
    call convert_x_to_xbox( self, position_old, dx_old, box_old )

    !print*, dx_old, box_old, dx_new, box_new

    do i=1, self%n_span
       self%spline_val(i) =  sll_f_spline_pp_horner_1d(self%spline_degree+1, self%spline_pp%poly_coeffs_fpa, dx_old, i) &
            * self%delta_x
       self%spline_val_more(i) =  sll_f_spline_pp_horner_1d(self%spline_degree+1, self%spline_pp%poly_coeffs_fpa, dx_new, i) &
            * self%delta_x
    end do

    ! New version to avoid the cancellation
    if (box_old<box_new) then
       i=0
       do ind = box_old-self%spline_degree, min(box_new-self%spline_degree-1, box_old )
          index = modulo(ind-1,self%n_cells)+1
          i = i+1
          j_dofs(index) = j_dofs(index) + (marker_charge * &
               (self%delta_x-self%spline_val(i)))
       end do
       j=0
       do ind = box_new-self%spline_degree, box_old
          j = j+1
          i = i+1
          index = modulo(ind-1,self%n_cells)+1
          j_dofs(index) = j_dofs(index) + marker_charge * &
               ((self%spline_val_more(j) - self%spline_val(i)))
       end do
       do ind =  max(box_old+1, box_new-self%spline_degree), box_new
          j = j+1
          index = modulo(ind-1,self%n_cells)+1
          j_dofs(index) = j_dofs(index) + (marker_charge * &
               (self%spline_val_more(j) ))
       end do

       do ind = box_old+1, box_new-self%spline_degree-1
          index = modulo(ind-1, self%n_cells)+1
          j_dofs(index)  = j_dofs(index) + marker_charge * self%delta_x
       end do
    else

       j=0
       do ind = box_new-self%spline_degree, min(box_old-self%spline_degree-1, box_new)
          index = modulo(ind-1,self%n_cells)+1
          j = j+1
          j_dofs(index) = j_dofs(index) + (marker_charge * &
               (-self%delta_x+self%spline_val_more(j)))
       end do
       i=0
       do ind = box_old-self%spline_degree, box_new
          j = j+1
          i = i+1
          index = modulo(ind-1,self%n_cells)+1
          j_dofs(index) = j_dofs(index) + (marker_charge * &
               (self%spline_val_more(j) - self%spline_val(i)))
       end do
       do ind =  max(box_new+1, box_old-self%spline_degree), box_old
          i = i+1
          index = modulo(ind-1,self%n_cells)+1
          j_dofs(index) = j_dofs(index) + (marker_charge * &
               (-self%spline_val(i) ))
       end do

       do ind = box_new+1, box_old-self%spline_degree-1
          index = modulo(ind-1, self%n_cells)+1
          j_dofs(index)  = j_dofs(index) - marker_charge * self%delta_x
       end do
    end if
    
  end subroutine add_current_spline_1d

  !> Combines add_current and evaluate_int
  subroutine add_current_evaluate_spline_1d(self, position_old, position_new, marker_charge, vbar, field_dofs, j_dofs, field)
    class( sll_t_particle_mesh_coupling_spline_1d ), intent(inout)   :: self !< kernel smoother object
    sll_real64,                               intent( in )    :: position_old(self%dim) !< Position of the particle
    sll_real64,                               intent( in )    :: position_new(self%dim) !< Position of the particle
    sll_real64,                               intent( in )    :: marker_charge !< Particle weights time charge
    sll_real64,                               intent( in )    :: vbar !< Particle weights time charge
    sll_real64,                               intent( in    ) :: field_dofs(self%n_cells) !< Coefficient vector of the current density
    sll_real64,                               intent( inout ) :: j_dofs(self%n_cells) !< Coefficient vector of the current density
    sll_real64,                               intent( out   ) :: field !< Efield

    sll_real64 :: dx_new, dx_old
    sll_int32 :: box_new, box_old, i, index, ind, j

    call convert_x_to_xbox( self, position_new, dx_new, box_new )
    call convert_x_to_xbox( self, position_old, dx_old, box_old )


    do i=1, self%n_span
       self%spline_val(i) =  sll_f_spline_pp_horner_1d(self%spline_degree+1, self%spline_pp%poly_coeffs_fpa, dx_old, i) &
            * self%delta_x
       self%spline_val_more(i) =  sll_f_spline_pp_horner_1d(self%spline_degree+1, self%spline_pp%poly_coeffs_fpa, dx_new, i) &
            * self%delta_x
    end do

    field = 0.0_f64


    ! New version to avoid the cancellation
    if (box_old<box_new) then
       i=0
       do ind = box_old-self%spline_degree, min(box_new-self%spline_degree-1, box_old)
          index = modulo(ind-1,self%n_cells)+1
          i = i+1
          j_dofs(index) = j_dofs(index) + (marker_charge * &
               (self%delta_x-self%spline_val(i)))
          field = field + ((self%delta_x - self%spline_val(i))) * field_dofs(index)
       end do
       j=0
       do ind = box_new-self%spline_degree, box_old
          j = j+1
          i = i+1
          index = modulo(ind-1,self%n_cells)+1
          j_dofs(index) = j_dofs(index) + marker_charge * &
               ((self%spline_val_more(j) - self%spline_val(i)))
          field = field + ((self%spline_val_more(j) - self%spline_val(i))) * field_dofs(index)
       end do

       do ind =  max(box_old+1, box_new-self%spline_degree), box_new
          j = j+1
          index = modulo(ind-1,self%n_cells)+1
          j_dofs(index) = j_dofs(index) + (marker_charge * &
               (self%spline_val_more(j) ))
          field = field + (self%spline_val_more(j) ) * field_dofs(index)
       end do

       do ind = box_old+1, box_new-self%spline_degree-1
          index = modulo(ind-1, self%n_cells)+1
          j_dofs(index)  = j_dofs(index) + marker_charge * self%delta_x
          field = field + self%delta_x * field_dofs(index)
       end do
    else

       j=0
       do ind = box_new-self%spline_degree, min(box_old-self%spline_degree-1, box_new)
          index = modulo(ind-1,self%n_cells)+1
          j = j+1
          j_dofs(index) = j_dofs(index) + (marker_charge * &
               (-self%delta_x+self%spline_val_more(j)))
          field = field + (-self%delta_x+self%spline_val_more(j)) * field_dofs(index)
       end do
       i=0
       do ind = box_old-self%spline_degree, box_new
          j = j+1
          i = i+1
          index = modulo(ind-1,self%n_cells)+1
          j_dofs(index) = j_dofs(index) + (marker_charge * &
               (self%spline_val_more(j) - self%spline_val(i)))
          field = field + (self%spline_val_more(j) - self%spline_val(i)) * field_dofs(index)
       end do
       do ind =  max(box_new+1, box_old-self%spline_degree), box_old
          i = i+1
          index = modulo(ind-1,self%n_cells)+1
          j_dofs(index) = j_dofs(index) + (marker_charge * &
               (-self%spline_val(i) ))
          field = field + (-self%spline_val(i) ) * field_dofs(index)
       end do

       do ind = box_new+1, box_old-self%spline_degree-1
          index = modulo(ind-1, self%n_cells)+1
          j_dofs(index)  = j_dofs(index) - marker_charge * self%delta_x
          field = field - self%delta_x * field_dofs(index)
       end do

    end if

    field = field/vbar

  end subroutine add_current_evaluate_spline_1d

  !> Evaluates the integral int_{poisition_old}^{position_new} field(x) d x
  subroutine evaluate_int_spline_1d(self, position_old, position_new, vbar, field_dofs, field)
    class( sll_t_particle_mesh_coupling_spline_1d ), intent(inout)   :: self !< kernel smoother object
    sll_real64,                               intent( in )    :: position_old(self%dim) !< Position of the particle
    sll_real64,                               intent( in )    :: position_new(self%dim) !< Position of the particle
    sll_real64,                               intent( in )    :: vbar !< Particle weights time charge
    sll_real64,                               intent( in    ) :: field_dofs(self%n_cells) !< Coefficient vector of the current density
    sll_real64,                               intent( out   ) :: field !< Efield

    sll_real64 :: dx_new, dx_old
    sll_int32 :: box_new, box_old, i, index, ind, j

    call convert_x_to_xbox( self, position_new, dx_new, box_new )
    call convert_x_to_xbox( self, position_old, dx_old, box_old )


    do i=1, self%n_span
       self%spline_val(i) =  sll_f_spline_pp_horner_1d(self%spline_degree+1, self%spline_pp%poly_coeffs_fpa, dx_old, i) &
            * self%delta_x
       self%spline_val_more(i) =  sll_f_spline_pp_horner_1d(self%spline_degree+1, self%spline_pp%poly_coeffs_fpa, dx_new, i) &
            * self%delta_x
    end do

    field = 0.0_f64


    ! New version to avoid the cancellation
    if (box_old<box_new) then
       i=0
       do ind = box_old-self%spline_degree, min(box_new-self%spline_degree-1, box_old)
          index = modulo(ind-1,self%n_cells)+1
          i = i+1
          field = field + ((self%delta_x - self%spline_val(i))) * field_dofs(index)
       end do
       j=0
       do ind = box_new-self%spline_degree, box_old
          j = j+1
          i = i+1
          index = modulo(ind-1,self%n_cells)+1
          field = field + ((self%spline_val_more(j) - self%spline_val(i))) * field_dofs(index)
       end do
       do ind =  max(box_old+1, box_new-self%spline_degree), box_new
          j = j+1
          index = modulo(ind-1,self%n_cells)+1
          field = field + (self%spline_val_more(j) ) * field_dofs(index)
       end do
       do ind = box_old+1, box_new-self%spline_degree-1
          index = modulo(ind-1, self%n_cells)+1
          field = field + self%delta_x * field_dofs(index)
       end do
    else

       j=0
       do ind = box_new-self%spline_degree, min(box_old-self%spline_degree-1, box_new)
          index = modulo(ind-1,self%n_cells)+1
          j = j+1
          field = field + (-self%delta_x+self%spline_val_more(j)) * field_dofs(index)
       end do
       i=0
       do ind = box_old-self%spline_degree, box_new
          j = j+1
          i = i+1
          index = modulo(ind-1,self%n_cells)+1
          field = field + (self%spline_val_more(j) - self%spline_val(i)) * field_dofs(index)
       end do
       do ind =  max(box_new+1,box_old-self%spline_degree), box_old
          i = i+1
          index = modulo(ind-1,self%n_cells)+1
          field = field + (-self%spline_val(i) ) * field_dofs(index)
       end do

       do ind = box_new+1, box_old-self%spline_degree-1
          index = modulo(ind-1, self%n_cells)+1
          field = field - self%delta_x * field_dofs(index)
       end do

    end if

    field = field/vbar

  end subroutine evaluate_int_spline_1d

  !> Evaluates an integrated field (between position_old and position_new) but with a quadrature formula for the integration
  subroutine evaluate_int_quad_spline_1d(self, position_old, position_new, field_dofs, field)
    class( sll_t_particle_mesh_coupling_spline_1d ), intent(inout)   :: self !< kernel smoother object
    sll_real64,                               intent( in )    :: position_old(self%dim) !< Position of the particle
    sll_real64,                               intent( in )    :: position_new(self%dim) !< Position of the particle
    sll_real64,                               intent( in    ) :: field_dofs(self%n_cells) !< Coefficient vector of the current density
    sll_real64,                               intent( out   ) :: field !< Efield
    !local variables
    sll_real64 :: dx_new, dx_old
    sll_int32 :: box_new, box_old, ind
    sll_real64 :: rinterval, lowert, uppert

    call convert_x_to_xbox( self, position_new, dx_new, box_new )
    call convert_x_to_xbox( self, position_old, dx_old, box_old )

    rinterval = self%delta_x/(position_new(1)-position_old(1))

    field = 0.0_f64

    if ( box_old == box_new ) then
       call  update_v_partial(self, dx_new, dx_old, 1.0_f64, 0.0_f64, box_old, field_dofs, field )
    elseif ( box_old < box_new ) then
       lowert = 0.0_f64
       uppert = (1.0_f64-dx_old)*rinterval
       call update_v_partial( self, 1.0_f64, dx_old, uppert, lowert, box_old, field_dofs, field )
       do ind = box_old+1, box_new-1
          lowert = uppert
          uppert = uppert + rinterval
          call update_v_partial( self, 1.0_f64, 0.0_f64, uppert, lowert, ind, field_dofs, field )
       end do
       call update_v_partial( self, dx_new, 0.0_f64, 1.0_f64, uppert, box_new, field_dofs, field )
    else
       lowert = 1.0_f64
       uppert = 1.0_f64 + (rinterval-dx_new*rinterval)
       call update_v_partial( self, 1.0_f64, dx_new, uppert, lowert, box_new, field_dofs, field )
       do ind = box_new+1, box_old-1
          lowert = uppert
          uppert = lowert + rinterval
          call update_v_partial( self, 1.0_f64, 0.0_f64, uppert, lowert, ind, field_dofs, field )
       end do
       call update_v_partial( self, dx_old, 0.0_f64, 0.0_f64, uppert, box_old, field_dofs, field )
       field = - field ! Account for the fact that we have integrated the wrong direction
    end if

  end subroutine evaluate_int_quad_spline_1d

  ! Evaluation of the integral needed for Newton's method in subcycling
  subroutine evaluate_int_spline_1d_subnewton(self, position_old, position_new, field_dofs, field)
    class( sll_t_particle_mesh_coupling_spline_1d ), intent(inout)   :: self !< kernel smoother object
    sll_real64,                               intent( in )    :: position_old(self%dim) !< Position of the particle
    sll_real64,                               intent( in )    :: position_new(self%dim) !< Position of the particle
    sll_real64,                               intent( in    ) :: field_dofs(self%n_cells) !< Coefficient vector of the current density
    sll_real64,                               intent( out   ) :: field
    !local variables
    sll_real64 :: dx_new, dx_old
    sll_int32 :: box_new, box_old, ind
    sll_real64 :: rinterval, lowert, uppert

    call convert_x_to_xbox( self, position_new, dx_new, box_new )
    call convert_x_to_xbox( self, position_old, dx_old, box_old )

    rinterval = self%delta_x/(position_new(1)-position_old(1))

    field = 0.0_f64

    if ( box_old == box_new ) then
       call  update_v_partial_newton(self, dx_new, dx_old, 1.0_f64, 0.0_f64, box_old, field_dofs, field )
    elseif ( box_old < box_new ) then
       lowert = 0.0_f64
       uppert = (1.0_f64-dx_old)*rinterval
       call update_v_partial_newton( self, 1.0_f64, dx_old, uppert, lowert, box_old, field_dofs, field )
       do ind = box_old+1, box_new-1
          lowert = uppert
          uppert = uppert + rinterval
          call update_v_partial_newton( self, 1.0_f64, 0.0_f64, uppert, lowert, ind, field_dofs, field )
       end do
       call update_v_partial_newton( self, dx_new, 0.0_f64, 1.0_f64, uppert, box_new, field_dofs, field )
    else
       lowert = 1.0_f64
       uppert = 1.0_f64 + (rinterval-dx_new*rinterval)
       call update_v_partial_newton( self, 1.0_f64, dx_new, uppert, lowert, box_new, field_dofs, field )
       do ind = box_new+1, box_old-1
          lowert = uppert
          uppert = lowert + rinterval
          call update_v_partial_newton( self, 1.0_f64, 0.0_f64, uppert, lowert, ind, field_dofs, field )
       end do
       call update_v_partial_newton( self, dx_old, 0.0_f64, 0.0_f64, uppert, box_old, field_dofs, field )
       field = - field ! Account for the fact that we have integrated the wrong direction
    end if

  end subroutine evaluate_int_spline_1d_subnewton


  !> Add current for one particle and update v (according to H_p1 part in Hamiltonian splitting)
  subroutine add_current_update_v_spline_pp_1d (self, position_old, position_new, marker_charge, qoverm, bfield_dofs, vi, j_dofs)
    class(sll_t_particle_mesh_coupling_spline_1d), intent(inout) :: self !< kernel smoother object
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

    index_old = index_old - 1
    index_new = index_new - 1

    if (index_old == index_new) then
       call update_jv_pp( self, r_old, r_new, index_old, marker_charge, &
            qoverm,  vi(2), j_dofs, bfield_dofs)
    elseif (index_old < index_new) then
       call update_jv_pp ( self, r_old, 1.0_f64, index_old, marker_charge, &
            qoverm, vi(2), j_dofs, bfield_dofs)
       call update_jv_pp ( self, 0.0_f64, r_new, index_new, marker_charge, &
            qoverm, vi(2), j_dofs, bfield_dofs)
       do ind = index_old+1, index_new-1
          call update_jv_pp ( self, 0.0_f64, 1.0_f64, ind, marker_charge, &
               qoverm, vi(2), j_dofs, bfield_dofs)
       end do
    else
       call update_jv_pp ( self, 1.0_f64,r_new,  index_new, marker_charge, qoverm, &
            vi(2), j_dofs, bfield_dofs)
       call update_jv_pp ( self, r_old, 0.0_f64, index_old, marker_charge, qoverm, &
            vi(2), j_dofs, bfield_dofs)
       do ind = index_new+1, index_old-1
          call update_jv_pp ( self, 1.0_f64,0.0_f64,  ind, marker_charge, qoverm, &
               vi(2), j_dofs, bfield_dofs)
       end do
    end if

  end subroutine add_current_update_v_spline_pp_1d

  !> Helper function for \a add_current_update_v.
  subroutine update_jv_pp(self, lower, upper, index, marker_charge, qoverm, vi, j_dofs, bfield_dofs)
    class(sll_t_particle_mesh_coupling_spline_1d), intent(inout) :: self !< time splitting object 
    sll_real64,                             intent(in)    :: lower
    sll_real64,                             intent(in)    :: upper
    sll_int32,                              intent(in)    :: index
    sll_real64,                             intent(in)    :: marker_charge
    sll_real64,                             intent(in)    :: qoverm
    sll_real64,                             intent(inout) :: vi
    sll_real64,                             intent(in)    :: bfield_dofs(self%n_cells)
    sll_real64,                             intent(inout) :: j_dofs(self%n_cells)
    !Local variables
    sll_int32  :: ind, i_grid, i_mod, n_cells

    n_cells = self%n_cells
    !Evaluation of the primitive integral at the lower and upper bound of the gridcell
    call sll_s_spline_pp_horner_primitive_1d(self%spline_val, self%spline_degree, self%spline_pp%poly_coeffs_fp, lower) 
    call sll_s_spline_pp_horner_primitive_1d(self%spline_val_more, self%spline_degree, self%spline_pp%poly_coeffs_fp, upper) 

    self%spline_val = (self%spline_val_more - self%spline_val) *(self%delta_x) 


    ind = 1
    do i_grid = index - self%spline_degree , index
       i_mod = modulo(i_grid, n_cells ) + 1
       j_dofs(i_mod) = j_dofs(i_mod) + &
            (marker_charge*self%spline_val(ind)* self%scaling)
       vi = vi - qoverm* self%spline_val(ind)*bfield_dofs(i_mod)
       ind = ind + 1
    end do

  end subroutine update_jv_pp

  !---------------------------------------------------------------------------!
  !> Add current for one particle and update v (according to H_p1 part in Hamiltonian splitting)
  subroutine add_current_update_v_spline_1d (self, position_old, position_new, marker_charge, qoverm, bfield_dofs, vi, j_dofs)
    class(sll_t_particle_mesh_coupling_spline_1d), intent(inout) :: self !< kernel smoother object
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

    index_old = index_old - 1
    index_new = index_new - 1

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


  end subroutine add_current_update_v_spline_1d

  !> Helper function for \a add_current_update_v.
  subroutine update_jv(self, lower, upper, index, marker_charge, qoverm, sign, vi, j_dofs, bfield_dofs)
    class(sll_t_particle_mesh_coupling_spline_1d), intent(inout) :: self !< time splitting object 
    sll_real64,                             intent(in)    :: lower
    sll_real64,                             intent(in)    :: upper
    sll_int32,                              intent(in)    :: index
    sll_real64,                             intent(in)    :: marker_charge
    sll_real64,                             intent(in)    :: qoverm
    sll_real64,                             intent(in)    :: sign
    sll_real64,                             intent(inout) :: vi
    sll_real64,                             intent(in)    :: bfield_dofs(self%n_cells)
    sll_real64,                             intent(inout) :: j_dofs(self%n_cells)
    !Local variables
    sll_int32  :: ind, i_grid, i_mod, n_cells, j
    sll_real64 :: c1, c2

    n_cells = self%n_cells

    c1 =  0.5_f64*(upper-lower)
    c2 =  0.5_f64*(upper+lower)


    call sll_s_uniform_bsplines_eval_basis(self%spline_degree, c1*self%quad_xw(1,1)+c2, &
         self%spline_val)
    self%spline_val = self%spline_val * (self%quad_xw(2,1)*c1)
    do j=2,self%n_quad_points
       call sll_s_uniform_bsplines_eval_basis(self%spline_degree, c1*self%quad_xw(1,j)+c2, &
            self%spline_val_more)
       self%spline_val = self%spline_val + self%spline_val_more * (self%quad_xw(2,j)*c1)
    end do
    self%spline_val = self%spline_val * (sign*self%delta_x)

    ind = 1
    do i_grid = index - self%spline_degree , index
       i_mod = modulo(i_grid, n_cells ) + 1 
       j_dofs(i_mod) = j_dofs(i_mod) + &
            (marker_charge*self%spline_val(ind)* self%scaling)
       vi = vi - qoverm* self%spline_val(ind)*bfield_dofs(i_mod)
       ind = ind + 1
    end do

  end subroutine update_jv

  !> Helper function for evaluate_int_quad (despite the name it evaluates a field and does not update v)
  subroutine update_v_partial(self, upper, lower, uppert, lowert, index, field_dofs, bint )
    class( sll_t_particle_mesh_coupling_spline_1d ), intent(inout)   :: self !< kernel smoother object
    sll_real64, intent( in    ) :: upper
    sll_real64, intent( in    ) :: lower
    sll_real64, intent( in    ) :: uppert
    sll_real64, intent( in    ) :: lowert
    sll_int32,  intent( in    ) :: index
    sll_real64, intent( in    ) :: field_dofs(self%n_cells) !< Coefficient vector of the current density
    sll_real64, intent( inout ) :: bint
    !local variables    
    sll_int32 :: j, ind, i_mod, i_grid, n_cells
    sll_real64 :: c1, c2, c1t, c2t


    n_cells = self%n_cells

    c1 = 0.5_f64*(upper-lower)
    c2 = 0.5_f64*(upper+lower)

    c1t = 0.5_f64*(uppert-lowert)
    c2t = 0.5_f64*(uppert+lowert)

    call sll_s_uniform_bsplines_eval_basis ( self%spline_degree, c1*self%quadp1_xw(1,1) + c2, &
         self%spline_val)
    self%spline_val = self%spline_val * (self%quadp1_xw(2,1)*c1t) * ( c1t* self%quadp1_xw(1,1) + c2t )

    do j=2, self%n_quadp1_points
       call sll_s_uniform_bsplines_eval_basis( self%spline_degree, c1*self%quadp1_xw(1,j)+c2, &
            self%spline_val_more )
       self%spline_val = self%spline_val + self%spline_val_more * (self%quadp1_xw(2,j)*c1t) *  ( c1t* self%quadp1_xw(1,j) + c2t )
    end do
    self%spline_val = self%spline_val !* self%delta_x(1)

    ind = 1
    do i_grid = index - self%spline_degree, index
       i_mod = modulo(i_grid-1, n_cells ) + 1
       bint = bint + self%spline_val(ind) * field_dofs( i_mod )
       ind = ind+1
    end do

  end subroutine update_v_partial

  

  !> Helper function for evaluate_int_quad_subnewton (despite the name it evaluates a field and does not update v)
  subroutine update_v_partial_newton(self, upper, lower, uppert, lowert, index, field_dofs, bint )
    class( sll_t_particle_mesh_coupling_spline_1d ), intent(inout)   :: self !< kernel smoother object
    sll_real64, intent( in    ) :: upper
    sll_real64, intent( in    ) :: lower
    sll_real64, intent( in    ) :: uppert
    sll_real64, intent( in    ) :: lowert
    sll_int32,  intent( in    ) :: index
    sll_real64, intent( in    ) :: field_dofs(self%n_cells) !< Coefficient vector of the current density
    sll_real64, intent( inout ) :: bint

    sll_int32 :: j, ind, i_mod, i_grid, n_cells
    sll_real64 :: c1, c2, c1t, c2t, tau

    n_cells = self%n_cells

    c1 = 0.5_f64*(upper-lower)
    c2 = 0.5_f64*(upper+lower)

    c1t = 0.5_f64*(uppert-lowert)
    c2t = 0.5_f64*(uppert+lowert)


    call sll_s_uniform_bsplines_eval_deriv ( self%spline_degree, c1*self%quadp1_xw(1,1) + c2, &
         self%spline_val)
    tau = ( c1t* self%quadp1_xw(1,1) + c2t )
    self%spline_val = self%spline_val * (self%quadp1_xw(2,1)*c1t) * ( 1.0_f64 - tau ) * tau

    do j=2, self%n_quadp1_points
       call sll_s_uniform_bsplines_eval_deriv( self%spline_degree, c1*self%quadp1_xw(1,j)+c2, &
            self%spline_val_more )
       tau = ( c1t* self%quadp1_xw(1,j) + c2t )
       self%spline_val = self%spline_val + self%spline_val_more * (self%quadp1_xw(2,j)*c1t) *  ( 1.0_f64 - tau ) * tau
    end do
    self%spline_val = self%spline_val !* self%delta_x(1)

    ind = 1
    do i_grid = index - self%spline_degree, index
       i_mod = modulo(i_grid-1, n_cells ) + 1
       bint = bint + self%spline_val(ind) * field_dofs( i_mod )
       ind = ind+1
    end do

  end subroutine update_v_partial_newton


  !> Helper function that identifies in which box the particle is found and its normalized poistion in the box
  subroutine convert_x_to_xbox( self, position, xi, box )
    class( sll_t_particle_mesh_coupling_spline_1d ), intent(inout)   :: self !< kernel smoother object
    sll_real64,                               intent( in )    :: position(self%dim) !< Position of the particle
    sll_real64,                               intent( out )    :: xi !< Position of the particle
    sll_int32,                                intent( out )    :: box !< Position of the particle

    xi = (position(1) - self%domain(1)) /self%delta_x
    box = floor( xi ) + 1
    xi = xi - real(box-1, f64)

  end subroutine convert_x_to_xbox


  !> Initializer
  subroutine init_spline_1d( self, domain, n_cells, no_particles, spline_degree, smoothing_type )
    class( sll_t_particle_mesh_coupling_spline_1d), intent(out)  :: self !< kernel smoother object
    sll_int32,                               intent(in) :: n_cells !< number of DoFs (spline coefficients)
    sll_real64,                              intent(in) :: domain(2) !< x_min and x_max of the domain
    sll_int32,                               intent(in) :: no_particles !< number of particles
    sll_int32,                               intent(in) :: spline_degree !< Degree of smoothing kernel spline
    sll_int32,                               intent(in) :: smoothing_type !< Define if Galerkin or collocation smoothing for right scaling in accumulation routines 
    !local variables
    sll_int32 :: ierr


    self%dim = 1

    ! Store grid information
    self%domain = domain
    self%n_cells = n_cells
    self%n_dofs = n_cells
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

    ! normalized Gauss Legendre points and weights
    self%quad_xw = sll_f_gauss_legendre_points_and_weights(self%n_quad_points)
    call sll_s_spline_pp_init_1d( self%spline_pp, spline_degree, n_cells)

    ! For nonlinear disgradE
    self%n_quadp1_points = (self%spline_degree+3)/2
    allocate( self%quadp1_xw(2, self%n_quadp1_points) )
    allocate( self%spline_pol_val(self%n_span) )
    ! normalized Gauss Legendre points and weights
    self%quadp1_xw = sll_f_gauss_legendre_points_and_weights(self%n_quadp1_points)

    ! For implicit subcycling
    self%n_quadp2_points = (self%spline_degree+4)/2 
    allocate( self%quadp2_xw(2, self%n_quadp2_points) )
    ! normalized Gauss Legendre points and weights
    self%quadp2_xw = sll_f_gauss_legendre_points_and_weights(self%n_quadp2_points)


    ! For smoothed evaluate and add_charge
    allocate( self%quads1_xw(2, self%n_quadp1_points) )
    ! normalized Gauss Legendre points and weights to [0,1]
    self%quads1_xw = self%quadp1_xw
    self%quads1_xw(1,:) = (self%quads1_xw(1,:) + 1.0_f64)*0.5_f64
    self%quads1_xw(2,:) = self%quads1_xw(2,:)*0.5_f64



    ! For smoothed add_current
    self%n_quads2_points = (self%spline_degree+4)/2
    allocate( self%quads2_xw(2, self%n_quads2_points) )
    ! normalized Gauss Legendre points and weights to [0,1]
    self%quads2_xw = sll_f_gauss_legendre_points_and_weights(self%n_quads2_points)
    self%quads2_xw(1,:) = (self%quads2_xw(1,:) + 1.0_f64)*0.5_f64
    self%quads2_xw(2,:) = self%quads2_xw(2,:)*0.5_f64

    self%n_quad_points_line = (self%spline_degree+2)/2

    allocate( self%quad_xw_line(2,self%n_quad_points_line) )
    ! normalized Gauss Legendre points and weights
    self%quad_xw_line = sll_f_gauss_legendre_points_and_weights(self%n_quad_points_line)

    self%quad_xw_line(1,:) = 0.5_f64*(self%quad_xw_line(1,:)+1.0_f64)
    self%quad_xw_line(2,:) = 0.5_f64*self%quad_xw_line(2,:)

  end subroutine init_spline_1d

  !> Destructor
  subroutine free_spline_1d(self)
    class (sll_t_particle_mesh_coupling_spline_1d), intent( inout ) :: self !< Kernel smoother object 

    deallocate(self%spline_val)
    deallocate(self%spline_val_more)
    deallocate(self%quad_xw)
    call sll_s_spline_pp_free_1d(self%spline_pp)

  end subroutine free_spline_1d


  !-------------------------------------------------------------------------------------------
  !< Constructor for abstract type
  subroutine sll_s_new_particle_mesh_coupling_spline_1d_ptr(smoother, domain, n_cells, no_particles, spline_degree, smoothing_type)
    class( sll_c_particle_mesh_coupling_1d), pointer, intent(out) :: smoother !< kernel smoother object
    sll_int32,                              intent(in)  :: n_cells !< number of DoFs (spline coefficients)
    sll_real64,                             intent(in)  :: domain(2) !< x_min and x_max of the domain
    sll_int32,                              intent(in)  :: no_particles !< number of particles
    sll_int32,                              intent(in)  :: spline_degree !< Degree of smoothing kernel spline
    sll_int32,                              intent(in)  :: smoothing_type !< Define if Galerkin or collocation smoothing for right scaling in accumulation routines 

    !local variables
    sll_int32 :: ierr


    SLL_ALLOCATE( sll_t_particle_mesh_coupling_spline_1d :: smoother , ierr)
    SLL_ASSERT( ierr == 0)

    select type( smoother )
    type is ( sll_t_particle_mesh_coupling_spline_1d )
       call smoother%init( domain, n_cells, no_particles, spline_degree, smoothing_type )
    end select

  end subroutine sll_s_new_particle_mesh_coupling_spline_1d_ptr

  !----------------------------------------------------------------------------------
  !< Constructor for abstract type
  subroutine sll_s_new_particle_mesh_coupling_spline_1d(smoother, domain, n_cells, no_particles, spline_degree, smoothing_type)
    class( sll_c_particle_mesh_coupling_1d), allocatable, intent(out):: smoother !< kernel smoother object
    sll_int32,                                  intent(in) :: n_cells !< number of DoFs (spline coefficients)
    sll_real64,                                 intent(in) :: domain(2) !< x_min and x_max of the domain
    sll_int32,                                  intent(in) :: no_particles !< number of particles
    sll_int32,                                  intent(in) :: spline_degree !< Degree of smoothing kernel spline
    sll_int32,                                  intent(in) :: smoothing_type !< Define if Galerkin or collocation smoothing for right scaling in accumulation routines 

    !local variables
    sll_int32 :: ierr


    allocate( sll_t_particle_mesh_coupling_spline_1d :: smoother , stat=ierr)
    SLL_ASSERT( ierr == 0)

    select type( smoother )
    type is ( sll_t_particle_mesh_coupling_spline_1d )
       call smoother%init( domain, n_cells, no_particles, spline_degree, smoothing_type )
    end select

  end subroutine sll_s_new_particle_mesh_coupling_spline_1d




  !> Add current version with tau and (1-tau) which is needed for the implicit subcycling scheme
  subroutine add_current_split_spline_1d(self, position_old, position_new, iter, total_iter, marker_charge, j_dofs )
    class( sll_t_particle_mesh_coupling_spline_1d ), intent(inout)   :: self !< kernel smoother object
    sll_real64,                               intent( in )    :: position_old(self%dim) !< Position of the particle
    sll_real64,                               intent( in )    :: position_new(self%dim) !< Position of the particle
    sll_int32,                                intent( in    ) :: iter
    sll_int32,                                intent( in    ) :: total_iter
    sll_real64,                               intent( in )    :: marker_charge !< Particle weights time charge
    sll_real64,                               intent( inout ) :: j_dofs(self%n_cells,2) !< Coefficient vector of the current density
    !local variables
    sll_real64 :: dx_new, dx_old
    sll_int32 :: box_new, box_old, ind
    sll_real64 :: rinterval, lowert, uppert

    dx_new = (position_new(1)-self%domain(1))/self%delta_x
    box_new = floor(dx_new)+1
    dx_new = dx_new - real(box_new-1,f64)
    dx_old = (position_old(1)-self%domain(1))/self%delta_x
    box_old = floor(dx_old)+1
    dx_old = dx_old - real(box_old-1,f64)

    if ( box_old == box_new ) then
       call  update_current_partial(self, dx_new, dx_old, 1.0_f64, 0.0_f64, box_old, &
            iter, total_iter, marker_charge, j_dofs )
    elseif ( box_old < box_new ) then
       rinterval = self%delta_x/(position_new(1)-position_old(1))
       lowert = 0.0_f64
       uppert = (1.0_f64-dx_old)*rinterval
       call update_current_partial( self, 1.0_f64, dx_old, uppert, lowert, box_old, &
            iter, total_iter, marker_charge, j_dofs )
       do ind = box_old+1, box_new-1
          lowert = uppert
          uppert = uppert + rinterval
          call update_current_partial( self, 1.0_f64, 0.0_f64, uppert, lowert, ind, &
               iter, total_iter, marker_charge, j_dofs )
       end do
       call update_current_partial( self, dx_new, 0.0_f64, 1.0_f64, uppert, box_new, &
            iter, total_iter, marker_charge, j_dofs )
    else
       rinterval = self%delta_x/(position_new(1)-position_old(1))
       lowert = 1.0_f64
       uppert = 1.0_f64 + (rinterval - dx_new*rinterval)
       ! Note the negative sign in the marker_charge accounts for the fact that we are integrating in the wrong direction
       call update_current_partial( self, 1.0_f64, dx_new, uppert, lowert, box_new, &
            iter, total_iter, -marker_charge, j_dofs )
       do ind = box_new+1, box_old-1
          lowert = uppert
          uppert = lowert + rinterval
          call update_current_partial( self, 1.0_f64, 0.0_f64, uppert, lowert, ind, &
               iter, total_iter, -marker_charge, j_dofs )
       end do
       call update_current_partial( self, dx_old, 0.0_f64, 0.0_f64, uppert, box_old, &
            iter, total_iter, -marker_charge, j_dofs )
    end if

  end subroutine add_current_split_spline_1d


  !> Helper function for add_current_split
  subroutine update_current_partial(self, upper, lower, uppert, lowert, index, iter, total_iter,  marker_charge, j_dofs )
    class( sll_t_particle_mesh_coupling_spline_1d ), intent(inout)   :: self !< kernel smoother object
    sll_real64, intent( in    ) :: upper
    sll_real64, intent( in    ) :: lower
    sll_real64, intent( in    ) :: uppert
    sll_real64, intent( in    ) :: lowert
    sll_int32,  intent( in    ) :: index
    sll_int32,  intent( in    ) :: iter
    sll_int32,  intent( in    ) :: total_iter
    sll_real64, intent( in    ) :: marker_charge !< Particle weights time charge
    sll_real64, intent( inout ) :: j_dofs(self%n_cells,2) !< Coefficient vector of the current density


    sll_int32 :: j, ind, i_mod, i_grid, n_cells
    sll_real64 :: c1, c2, c1t, c2t
    sll_real64 :: time

    n_cells = self%n_cells

    c1 = 0.5_f64*(upper-lower)
    c2 = 0.5_f64*(upper+lower)

    c1t = 0.5_f64*(uppert-lowert)
    c2t = 0.5_f64*(uppert+lowert)

    call sll_s_uniform_bsplines_eval_basis ( self%spline_degree, &
         c1*self%quadp1_xw(1,1) + c2, &
         self%spline_val)
    time = (c1t* self%quadp1_xw(1,1) + c2t + real(iter,f64))/real(total_iter, f64)
    self%spline_pol_val = self%spline_val * (self%quadp1_xw(2,1)*c1t) * &
         ( 1.0_f64 - time )
    self%spline_val = self%spline_val * (self%quadp1_xw(2,1)*c1t) * time

    do j=2, self%n_quadp1_points
       call sll_s_uniform_bsplines_eval_basis( self%spline_degree, &
            c1*self%quadp1_xw(1,j)+c2, &
            self%spline_val_more )
       time = (c1t* self%quadp1_xw(1,j) + c2t + real(iter,f64))/real(total_iter, f64)
       self%spline_val = self%spline_val + &
            self%spline_val_more * (self%quadp1_xw(2,j)*c1t) *  time
       self%spline_pol_val = self%spline_pol_val + &
            self%spline_val_more * (self%quadp1_xw(2,j)*c1t) *  ( 1.0_f64 - time )
    end do

    ind = 1
    do i_grid = index - self%spline_degree, index
       i_mod = modulo(i_grid-1, n_cells ) + 1
       j_dofs(i_mod,1) = j_dofs(i_mod,1) + marker_charge * self%spline_val(ind)
       j_dofs(i_mod,2) = j_dofs(i_mod,2) + marker_charge * self%spline_pol_val(ind)
       ind = ind+1
    end do

  end subroutine update_current_partial

  !> Evaluates the integrals with tau function for implicit subcycling
  subroutine evaluate_int_linear_quad_spline_1d(self, position_old, position_new, iter, total_iter, field_dofs_1, field_dofs_2, field)
    class( sll_t_particle_mesh_coupling_spline_1d ), intent(inout)   :: self !< kernel smoother object
    sll_real64,                               intent( in )    :: position_old(self%dim) !< Position of the particle
    sll_real64,                               intent( in )    :: position_new(self%dim) !< Position of the particle
    sll_real64,                               intent( in    ) :: field_dofs_1(self%n_cells) !< Coefficient vector of the current density
    sll_real64,                               intent( in    ) :: field_dofs_2(self%n_cells) !< Coefficient vector of the current density
    sll_int32,  intent( in    ) :: iter
    sll_int32,  intent( in    ) :: total_iter
    sll_real64,                               intent( out   ) :: field(2)
    !local variables
    sll_real64 :: dx_new, dx_old
    sll_int32 :: box_new, box_old, ind
    sll_real64 :: rinterval, lowert, uppert

    print*, position_new(1)
    print*, self%domain(1)
    print*, self%delta_x
    
    dx_new = (position_new(1)-self%domain(1))/self%delta_x
    box_new = floor(dx_new)+1
    dx_new = dx_new - real(box_new-1,f64)
    dx_old = (position_old(1)-self%domain(1))/self%delta_x
    box_old = floor(dx_old)+1
    dx_old = dx_old - real(box_old-1,f64)

    if ( box_old ==  box_new ) then
       rinterval = 0.0_f64
    else
       rinterval = self%delta_x/(position_new(1)-position_old(1))
    end if
    field = 0.0_f64

    if ( box_old == box_new ) then
       call  update_field_linear_partial(self, dx_new, dx_old, 1.0_f64, 0.0_f64, box_old, &
            field_dofs_1, field_dofs_2, iter, total_iter, field )
    elseif ( box_old < box_new ) then
       lowert = 0.0_f64
       uppert = (1.0_f64-dx_old)*rinterval
       call update_field_linear_partial( self, 1.0_f64, dx_old, uppert, lowert, box_old, &
            field_dofs_1, field_dofs_2, iter, total_iter, field )
       do ind = box_old+1, box_new-1
          lowert = uppert
          uppert = uppert + rinterval
          call update_field_linear_partial( self, 1.0_f64, 0.0_f64, uppert, lowert, ind, &
               field_dofs_1, field_dofs_2, iter, total_iter, field )
       end do
       call update_field_linear_partial( self, dx_new, 0.0_f64, 1.0_f64, uppert, box_new, &
            field_dofs_1, field_dofs_2, iter, total_iter, field )
    else
       lowert = 1.0_f64
       uppert = 1.0_f64 + (rinterval-dx_new*rinterval)
       call update_field_linear_partial( self, 1.0_f64, dx_new, uppert, lowert, box_new, &
            field_dofs_1, field_dofs_2, iter, total_iter, field )
       do ind = box_new+1, box_old-1
          lowert = uppert
          uppert = lowert + rinterval
          call update_field_linear_partial( self, 1.0_f64, 0.0_f64, uppert, lowert, ind, &
               field_dofs_1, field_dofs_2, iter, total_iter, field )
       end do
       call update_field_linear_partial( self, dx_old, 0.0_f64, 0.0_f64, uppert, box_old, &
            field_dofs_1, field_dofs_2, iter, total_iter, field )
       field = - field ! Account for the fact that we have integrated the wrong direction
    end if

  end subroutine evaluate_int_linear_quad_spline_1d


  !> Helper function for evaluate_int_linear_quad (implements part within one cell)
  subroutine update_field_linear_partial(self, upper, lower, uppert, lowert, index, field_dofs_1, field_dofs_2, iter, total_iter, bint )
    class( sll_t_particle_mesh_coupling_spline_1d ), intent(inout)   :: self !< kernel smoother object
    sll_real64, intent( in    ) :: upper
    sll_real64, intent( in    ) :: lower
    sll_real64, intent( in    ) :: uppert
    sll_real64, intent( in    ) :: lowert
    sll_int32,  intent( in    ) :: index
    sll_real64, intent( in    ) :: field_dofs_1(self%n_cells) !< Coefficient vector of the current density
    sll_real64, intent( in    ) :: field_dofs_2(self%n_cells) !< Coefficient vector of the current density
    sll_int32,  intent( in    ) :: iter
    sll_int32,  intent( in    ) :: total_iter
    sll_real64, intent( inout ) :: bint(2)

    sll_int32 :: j, ind, i_mod, i_grid, n_cells
    sll_real64 :: c1, c2, c1t, c2t
    sll_real64 :: tau, time, weight, field_dof_t

    n_cells = self%n_cells

    c1 = 0.5_f64*(upper-lower)
    c2 = 0.5_f64*(upper+lower)

    c1t = 0.5_f64*(uppert-lowert)
    c2t = 0.5_f64*(uppert+lowert)

    do j=1, self%n_quadp2_points
       tau =  c1t* self%quadp2_xw(1,j) + c2t
       time = ( tau +real(iter,f64))/real(total_iter, f64)
       weight =  (self%quadp2_xw(2,j)*c1t)
       call sll_s_uniform_bsplines_eval_basis( self%spline_degree, c1*self%quadp2_xw(1,j)+c2, &
            self%spline_val )
       ind = 1
       do i_grid = index - self%spline_degree, index
          i_mod = modulo(i_grid-1, n_cells ) + 1
          field_dof_t = field_dofs_1( i_mod ) + time * ( field_dofs_2( i_mod ) - field_dofs_1( i_mod ) )
          bint(1) = bint(1) + self%spline_val(ind) * field_dof_t * weight * tau
          bint(2) = bint(2) + self%spline_val(ind) * field_dof_t * weight * (1.0_f64 - tau )
          ind = ind+1
       end do
    end do

  end subroutine update_field_linear_partial

  !> For implicit subcycling propagator, evaluates the integral needed for Newtons's method in the evaluate_int_linear_quad part.
  subroutine evaluate_int_linear_quad_subnewton_spline_1d(self, position_old, position_new, iter, total_iter, field_dofs_1, field_dofs_2, field)
    class( sll_t_particle_mesh_coupling_spline_1d ), intent(inout)   :: self !< kernel smoother object
    sll_real64,                               intent( in )    :: position_old(self%dim) !< Position of the particle
    sll_real64,                               intent( in )    :: position_new(self%dim) !< Position of the particle
    sll_real64,                               intent( in    ) :: field_dofs_1(self%n_cells) !< Coefficient vector of the current density
    sll_real64,                               intent( in    ) :: field_dofs_2(self%n_cells) !< Coefficient vector of the current density
    sll_int32,  intent( in    ) :: iter
    sll_int32,  intent( in    ) :: total_iter
    sll_real64,                               intent( out   ) :: field
    !local variables
    sll_real64 :: dx_new, dx_old
    sll_int32 :: box_new, box_old, ind
    sll_real64 :: rinterval, lowert, uppert

    dx_new = (position_new(1)-self%domain(1))/self%delta_x
    box_new = floor(dx_new)+1
    dx_new = dx_new - real(box_new-1,f64)
    dx_old = (position_old(1)-self%domain(1))/self%delta_x
    box_old = floor(dx_old)+1
    dx_old = dx_old - real(box_old-1,f64)

    rinterval = self%delta_x/(position_new(1)-position_old(1))

    field = 0.0_f64

    if ( box_old == box_new ) then
       call  update_field_linear_partial_newton(self, dx_new, dx_old, 1.0_f64, 0.0_f64, box_old, &
            field_dofs_1, field_dofs_2, iter, total_iter, field )
    elseif ( box_old < box_new ) then
       lowert = 0.0_f64
       uppert = (1.0_f64-dx_old)*rinterval
       call update_field_linear_partial_newton( self, 1.0_f64, dx_old, uppert, lowert, box_old, &
            field_dofs_1, field_dofs_2, iter, total_iter, field )
       do ind = box_old+1, box_new-1
          lowert = uppert
          uppert = uppert + rinterval
          call update_field_linear_partial_newton( self, 1.0_f64, 0.0_f64, uppert, lowert, ind, &
               field_dofs_1, field_dofs_2, iter, total_iter, field )
       end do
       call update_field_linear_partial_newton( self, dx_new, 0.0_f64, 1.0_f64, uppert, box_new, &
            field_dofs_1, field_dofs_2, iter, total_iter, field )
    else
       lowert = 1.0_f64
       uppert = 1.0_f64 + (rinterval-dx_new*rinterval)
       call update_field_linear_partial_newton( self, 1.0_f64, dx_new, uppert, lowert, box_new, &
            field_dofs_1, field_dofs_2, iter, total_iter, field )
       do ind = box_new+1, box_old-1
          lowert = uppert
          uppert = lowert + rinterval
          call update_field_linear_partial_newton( self, 1.0_f64, 0.0_f64, uppert, lowert, ind, &
               field_dofs_1, field_dofs_2, iter, total_iter, field )
       end do
       call update_field_linear_partial_newton( self, dx_old, 0.0_f64, 0.0_f64, uppert, box_old, &
            field_dofs_1, field_dofs_2, iter, total_iter, field )
       field = - field ! Account for the fact that we have integrated the wrong direction
    end if

  end subroutine evaluate_int_linear_quad_subnewton_spline_1d


  !> Helper function for evaluate_int_linear_quad_subnewton (takes care of the per cell computations)
  subroutine update_field_linear_partial_newton(self, upper, lower, uppert, lowert, index, field_dofs_1, field_dofs_2, iter, total_iter, bint )
    class( sll_t_particle_mesh_coupling_spline_1d ), intent(inout)   :: self !< kernel smoother object
    sll_real64, intent( in    ) :: upper
    sll_real64, intent( in    ) :: lower
    sll_real64, intent( in    ) :: uppert
    sll_real64, intent( in    ) :: lowert
    sll_int32,  intent( in    ) :: index
    sll_real64, intent( in    ) :: field_dofs_1(self%n_cells) !< Coefficient vector of the current density
    sll_real64, intent( in    ) :: field_dofs_2(self%n_cells) !< Coefficient vector of the current density
    sll_int32,  intent( in    ) :: iter
    sll_int32,  intent( in    ) :: total_iter
    sll_real64, intent( inout ) :: bint

    sll_int32 :: j, ind, i_mod, i_grid, n_cells
    sll_real64 :: c1, c2, c1t, c2t
    sll_real64 :: tau, time, weight, field_dof_t

    n_cells = self%n_cells

    c1 = 0.5_f64*(upper-lower)
    c2 = 0.5_f64*(upper+lower)

    c1t = 0.5_f64*(uppert-lowert)
    c2t = 0.5_f64*(uppert+lowert)

    do j=1, self%n_quadp2_points
       tau =  c1t* self%quadp2_xw(1,j) + c2t
       time = ( tau +real(iter,f64))/real(total_iter, f64)
       weight =  (self%quadp2_xw(2,j)*c1t)
       call sll_s_uniform_bsplines_eval_deriv ( self%spline_degree, c1*self%quadp2_xw(1,1) + c2, &
            self%spline_val)
       ind = 1
       do i_grid = index - self%spline_degree, index
          i_mod = modulo(i_grid-1, n_cells ) + 1
          field_dof_t = field_dofs_1( i_mod ) + time * ( field_dofs_2( i_mod ) - field_dofs_1( i_mod ) )
          bint = bint + self%spline_val(ind) * field_dof_t * weight * (1.0_f64 - tau ) * tau
          ind = ind+1
       end do
    end do

  end subroutine update_field_linear_partial_newton

  !> For explicit version of the subcycling propagator (to evaluate the charge as it comes out integrated from the formulation)
  subroutine add_charge_int_spline_1d(self, position_old, position_new, marker_charge, j_dofs )
    class( sll_t_particle_mesh_coupling_spline_1d ), intent(inout)   :: self !< kernel smoother object
    sll_real64,                               intent( in )    :: position_old(self%dim) !< Position of the particle
    sll_real64,                               intent( in )    :: position_new(self%dim) !< Position of the particle
    sll_real64,                               intent( in )    :: marker_charge !< Particle weights time charge
    sll_real64,                               intent( inout ) :: j_dofs(self%n_cells) !< Coefficient vector of the current density
    !local variables
    sll_real64 :: dx_new, dx_old
    sll_int32 :: box_new, box_old, ind
    sll_real64 :: rinterval, lowert, uppert
    
    dx_new = (position_new(1)-self%domain(1))/self%delta_x
    box_new = floor(dx_new)+1
    dx_new = dx_new - real(box_new-1,f64)
    dx_old = (position_old(1)-self%domain(1))/self%delta_x
    box_old = floor(dx_old)+1
    dx_old = dx_old - real(box_old-1,f64)

    if ( box_old == box_new ) then
       call  update_charge_int_partial(self, dx_new, dx_old, 1.0_f64, 0.0_f64, box_old, &
            marker_charge, j_dofs )
    elseif ( box_old < box_new ) then
       rinterval = self%delta_x/(position_new(1)-position_old(1))
       lowert = 0.0_f64
       uppert = (1.0_f64-dx_old)*rinterval
       call update_charge_int_partial( self, 1.0_f64, dx_old, uppert, lowert, box_old, &
            marker_charge, j_dofs )
       do ind = box_old+1, box_new-1
          lowert = uppert
          uppert = uppert + rinterval
          call update_charge_int_partial( self, 1.0_f64, 0.0_f64, uppert, lowert, ind, &
               marker_charge, j_dofs )
       end do
       call update_charge_int_partial( self, dx_new, 0.0_f64, 1.0_f64, uppert, box_new, &
            marker_charge, j_dofs )
    else
       rinterval = self%delta_x/(position_new(1)-position_old(1))
       lowert = 1.0_f64
       uppert = 1.0_f64 + (rinterval - dx_new*rinterval)
       ! Note the negative sign in the marker_charge accounts for the fact that we are integrating in the wrong direction
       call update_charge_int_partial( self, 1.0_f64, dx_new, uppert, lowert, box_new, &
            -marker_charge, j_dofs )
       do ind = box_new+1, box_old-1
          lowert = uppert
          uppert = lowert + rinterval
          call update_charge_int_partial( self, 1.0_f64, 0.0_f64, uppert, lowert, ind, &
               -marker_charge, j_dofs )
       end do
       call update_charge_int_partial( self, dx_old, 0.0_f64, 0.0_f64, uppert, box_old, &
            -marker_charge, j_dofs )
    end if

  end subroutine add_charge_int_spline_1d


  !> Helper function for add_charge_int (takes care of the per cell computations)
  subroutine update_charge_int_partial(self, upper, lower, uppert, lowert, index, marker_charge, j_dofs )
    class( sll_t_particle_mesh_coupling_spline_1d ), intent(inout)   :: self !< kernel smoother object
    sll_real64, intent( in    ) :: upper
    sll_real64, intent( in    ) :: lower
    sll_real64, intent( in    ) :: uppert
    sll_real64, intent( in    ) :: lowert
    sll_int32,  intent( in    ) :: index
    sll_real64, intent( in    ) :: marker_charge !< Particle weights time charge
    sll_real64, intent( inout ) :: j_dofs(self%n_cells) !< Coefficient vector of the current density
    !local variables
    sll_int32 :: j, ind, i_mod, i_grid, n_cells
    sll_real64 :: c1, c2, c1t, c2t
   
    n_cells = self%n_cells

    c1 = 0.5_f64*(upper-lower)
    c2 = 0.5_f64*(upper+lower)

    c1t = 0.5_f64*(uppert-lowert)
    c2t = 0.5_f64*(uppert+lowert)

    call sll_s_uniform_bsplines_eval_basis ( self%spline_degree, &
         c1*self%quadp1_xw(1,1) + c2, &
         self%spline_val)
    self%spline_val = self%spline_val * (self%quadp1_xw(2,1)*c1t) 

    do j=2, self%n_quadp1_points
       call sll_s_uniform_bsplines_eval_basis( self%spline_degree, &
            c1*self%quadp1_xw(1,j)+c2, &
            self%spline_val_more )
       self%spline_val = self%spline_val + &
            self%spline_val_more * (self%quadp1_xw(2,j)*c1t)
    end do

    ind = 1
    do i_grid = index - self%spline_degree, index
       i_mod = modulo(i_grid-1, n_cells ) + 1
       j_dofs(i_mod) = j_dofs(i_mod) + marker_charge * self%spline_val(ind)
       ind = ind+1
    end do

  end subroutine update_charge_int_partial


end module sll_m_particle_mesh_coupling_spline_1d
