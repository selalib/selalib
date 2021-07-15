!> @ingroup particle_mesh_coupling
!> @author Katharina Kormann, IPP
!> @brief Base class for kernel smoothers for accumulation and field evaluation in PIC.
!> @details This base class gives an abstract interface to the basic functions for accumulation of charge and current densities as well as the evaluation of a function at particle positions.
module sll_m_particle_mesh_coupling_base_1d

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  implicit none

  public :: &
       sll_p_collocation, &
       sll_p_galerkin, &
       sll_c_particle_mesh_coupling_1d

  private
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Define parameters to set if Galerkin or collocation scaling should be used in accumulation routines
  sll_int32, parameter :: sll_p_galerkin = 0
  sll_int32, parameter :: sll_p_collocation = 1

  !> Basic type of a kernel smoother used for PIC simulations
  type, abstract :: sll_c_particle_mesh_coupling_1d
     sll_int32              :: dim 
     sll_int32              :: n_cells !< Number of grid cells
     sll_int32              :: n_dofs  !< Number of degrees of freedom of the smoothing kernels.
     sll_int32              :: spline_degree !< Degree of smoothing kernel spline

   contains

     procedure(add_single), deferred        :: add_charge !> Add the contribution of one particle to the charge density
     procedure(add_array), deferred         :: add_particle_mass !> Add the contribution of one particle to the approximate mass matrix defined by the sum over all particles
     procedure(add_array), deferred         :: add_particle_mass_full !> Add the contribution of one particle to the approximate mass matrix defined by the sum over all particles
     procedure(eval_single), deferred       :: evaluate  !> Evaluate spline function with given coefficients
     procedure(eval_multiple), deferred     :: evaluate_multiple  !> Evaluate multiple spline function with given coefficients
     procedure(add_update) , deferred       :: add_current_update_v !> Add contribution of one particle to the current density and update velocity
     procedure(add_current), deferred       :: add_current !> Add contribution of one particle to the current density (integrated over x)
     procedure(empty), deferred             :: free !< Destructor


  end type sll_c_particle_mesh_coupling_1d

  !---------------------------------------------------------------------------!
  abstract interface
     subroutine add_single(self, position, marker_charge, rho_dofs) 
       use sll_m_working_precision
       import sll_c_particle_mesh_coupling_1d
       class (sll_c_particle_mesh_coupling_1d), intent( inout ) :: self !< Kernel smoother object
       sll_real64,                    intent( in )    :: position(self%dim) !< Position of the particle
       sll_real64,                    intent( in )    :: marker_charge !< Particle weight times charge
       sll_real64,                    intent( inout ) :: rho_dofs(self%n_dofs) !< Coefficient vector of the charge distribution

     end subroutine add_single
  end interface

  !---------------------------------------------------------------------------!
  abstract interface
     subroutine add_array(self, position, marker_charge, particle_mass) 
       use sll_m_working_precision
       import sll_c_particle_mesh_coupling_1d
       class (sll_c_particle_mesh_coupling_1d), intent( inout ) :: self !< Kernel smoother object
       sll_real64,                    intent( in )    :: position(self%dim) !< Position of the particle
       sll_real64,                    intent( in )    :: marker_charge !< Particle weight times charge
       sll_real64,                    intent( inout ) :: particle_mass(:, :) !< Coefficient vector of the charge distribution

     end subroutine add_array
  end interface

  !---------------------------------------------------------------------------!
  abstract interface
     subroutine eval_single(self, position, field_dofs, field_value)
       use sll_m_working_precision
       import sll_c_particle_mesh_coupling_1d
       class (sll_c_particle_mesh_coupling_1d), intent( inout ) :: self !< Kernel smoother object 
       sll_real64,                    intent( in )    :: position(self%dim) !< Position of the particle
       sll_real64,                    intent( in )    :: field_dofs(self%n_dofs) !< Coefficient vector for the field DoFs
       sll_real64, intent(out) :: field_value !< Value(s) of the electric fields at given position
     end subroutine eval_single
  end interface

  !---------------------------------------------------------------------------!
  abstract interface
     subroutine eval_multiple(self, position, components, field_dofs, field_value)
       use sll_m_working_precision
       import sll_c_particle_mesh_coupling_1d
       class (sll_c_particle_mesh_coupling_1d), intent( inout ) :: self !< Kernel smoother object 
       sll_real64,                    intent( in )    :: position(self%dim) !< Position of the particle
       sll_int32,                     intent(in)      :: components(:)
       sll_real64,                    intent( in )    :: field_dofs(:,:) !< Coefficient vector for the field DoFs
       sll_real64,                    intent(out)     :: field_value(:) !< Value(s) of the electric fields at given position
     end subroutine eval_multiple
  end interface

  !---------------------------------------------------------------------------!
  abstract interface
     subroutine add_current(self, position_old, position_new, marker_charge, j_dofs) 
       use sll_m_working_precision
       import sll_c_particle_mesh_coupling_1d
       class (sll_c_particle_mesh_coupling_1d), intent( inout ) :: self !< Kernel smoother object
       sll_real64,                    intent( in )    :: position_old(self%dim) !< Position of the particle
       sll_real64,                    intent( in )    :: position_new(self%dim) !< Position of the particle
       sll_real64,                    intent( in )    :: marker_charge !< Particle weight times charge
       sll_real64,                    intent( inout ) :: j_dofs(self%n_dofs) !< Coefficient vector of the charge distribution

     end subroutine add_current
  end interface

  !---------------------------------------------------------------------------!
  abstract interface
     subroutine add_update (self, &
          position_old, &
          position_new, &
          marker_charge, &
          qoverm, &
          bfield_dofs, &
          vi, &
          j_dofs)
       use sll_m_working_precision
       import sll_c_particle_mesh_coupling_1d
       class(sll_c_particle_mesh_coupling_1d), intent(inout) :: self !< kernel smoother object
       sll_real64, intent(in)    :: position_old(self%dim) !< Position at time t
       sll_real64, intent(in)    :: position_new(self%dim) !< Position at time t+\Delta t
       sll_real64, intent(in)    :: marker_charge          !< Particle weight times charge
       sll_real64, intent(in)    :: qoverm   !< charge to mass ratio
       sll_real64, intent(in)    :: bfield_dofs(:) !< values of the B-field at the dofs
       sll_real64, intent(inout) :: vi(:) !< Velocity of the particle
       sll_real64, intent(inout) :: j_dofs(:) !< Current at the DoFs

     end subroutine add_update
  end interface

  !---------------------------------------------------------------------------!  
  abstract interface
     subroutine empty(self)
       import sll_c_particle_mesh_coupling_1d
       class (sll_c_particle_mesh_coupling_1d), intent( inout ) :: self !< Kernel smoother object 

     end subroutine empty
  end interface

end module sll_m_particle_mesh_coupling_base_1d
