!> @ingroup particle_mesh_coupling
!> @author Katharina Kormann, IPP
!> @brief Base class for kernel smoothers for accumulation and field evaluation in PIC.
!> @details This base class gives an abstract interface to the basic functions for accumulation of charge and current densities as well as the evaluation of a function at particle positions.
module sll_m_particle_mesh_coupling_base_3d

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use sll_m_splines_pp, only : &
       sll_t_spline_pp_3d

  implicit none

  public :: &
       sll_c_particle_mesh_coupling_3d

  private
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Define parameters to set if Galerkin or collocation scaling should be used in accumulation routines
  sll_int32, parameter :: sll_p_galerkin = 0
  sll_int32, parameter :: sll_p_collocation = 1

  !> Basic type of a kernel smoother used for PIC simulations
  type, abstract :: sll_c_particle_mesh_coupling_3d
     sll_int32                :: dim = 3
     sll_real64               :: delta_x(3)  !< Value of grid spacing along all directions.
     sll_real64               :: domain(3,2) !< Definition of the domain: domain(1,1) = x1_min  domain(1,2) = x1_max
     sll_int32                :: n_total  !< product of the cell points in all three directions
     sll_int32                :: spline_degree(3) !< Degree of smoothing kernel spline
     sll_int32                :: n_cells(3) !< Number of cells points per dimension for use on tensor product cells based smoothing kernels.
     sll_int32                :: n_dofs(3) !< Number of degrees of freedom

     type(sll_t_spline_pp_3d) :: spline_pp_0 !< 3d pp-spline of given spline degree
     type(sll_t_spline_pp_3d) :: spline_pp_1 !< 3d pp-spline of given spline degree-1

   contains

     procedure(add_single), deferred           :: add_charge !> Add the contribution of one particle to the charge density
     procedure(add_array), deferred            :: add_particle_mass !> Collect diagonal entry of particle mass matrix
     procedure(add_array_mixed), deferred      :: add_particle_mass_od !> Collect off-diagonal entry of particle mass matrix
     procedure(eval_single), deferred          :: evaluate !> Evaluate spline function with given coefficients
     procedure(eval_multiple), deferred        :: evaluate_multiple !> Evaluate multiple spline function with given coefficients
     procedure(add_current), deferred          :: add_current !> Add contribution of one particle to the current density (integrated over x)
     procedure(add_current_evaluate), deferred :: add_current_evaluate !> Add contribution of one particle to the current density (integrated over x)
     procedure(add_update), deferred           :: add_current_update_v_component1 !> Add contribution of one particle to the current density in x1 direction and update velocity
     procedure(add_update), deferred           :: add_current_update_v_component2 !> Add contribution of one particle to the current density in x2 direction and update velocity
     procedure(add_update), deferred           :: add_current_update_v_component3 !> Add contribution of one particle to the current density in x3 direction and update velocity
     procedure(empty), deferred                :: free !< Destructor


  end type sll_c_particle_mesh_coupling_3d

  !---------------------------------------------------------------------------!
  abstract interface
     subroutine add_single(self, position, marker_charge, degree, rho_dofs)
       use sll_m_working_precision
       import sll_c_particle_mesh_coupling_3d
       class (sll_c_particle_mesh_coupling_3d), intent( inout ) :: self !< Kernel smoother object
       sll_real64,                               intent( in )    :: position(3) !< Position of the particle
       sll_real64,                               intent( in )    :: marker_charge !< Particle weights time charge
       sll_int32,                                intent( in )    :: degree(3) !< Spline degree along each dimension
       sll_real64,                               intent( inout ) :: rho_dofs(:)

     end subroutine add_single
  end interface

  !---------------------------------------------------------------------------!
  abstract interface
     subroutine add_array(self, position, marker_charge, degree, particle_mass )
       use sll_m_working_precision
       import sll_c_particle_mesh_coupling_3d
       class (sll_c_particle_mesh_coupling_3d), intent( inout ) :: self !< Kernel smoother object
       sll_real64,                               intent( in )    :: position(3) !< Position of the particle
       sll_real64,                               intent( in )    :: marker_charge !< Particle weights time charge
       sll_int32,                                intent( in )    :: degree(3) !< Spline degree along each dimension
       sll_real64,                    intent( inout ) :: particle_mass(:, :) !< Coefficient vector of the charge distribution

     end subroutine add_array
  end interface

  !---------------------------------------------------------------------------!
  abstract interface
     subroutine add_array_mixed(self, position, marker_charge, degree1, degree2, particle_mass )
       use sll_m_working_precision
       import sll_c_particle_mesh_coupling_3d
       class (sll_c_particle_mesh_coupling_3d), intent( inout ) :: self !< Kernel smoother object
       sll_real64,                               intent( in )    :: position(3) !< Position of the particle
       sll_real64,                               intent( in )    :: marker_charge !< Particle weights time charge
       sll_int32,                                intent( in )    :: degree1(3), degree2(3) !< Spline degree along each dimension
       sll_real64,                    intent( inout ) :: particle_mass(:, :) !< Coefficient vector of the charge distribution

     end subroutine add_array_mixed
  end interface

  !---------------------------------------------------------------------------!
  abstract interface
     subroutine eval_single(self, position, degree, field_dofs, field_value)
       use sll_m_working_precision
       import sll_c_particle_mesh_coupling_3d
       class (sll_c_particle_mesh_coupling_3d), intent( inout ) :: self !< Kernel smoother object 
       sll_real64,                              intent( in )    :: position(3) !< Position of the particle
       sll_int32 ,                              intent( in )    :: degree(3) !< Spline degree of the various components
       sll_real64,                              intent( in )    :: field_dofs(:) !< Coefficient vector for the field DoFs
       sll_real64,                              intent( out )   :: field_value !< Value(s) of the fields at given position
     end subroutine eval_single
  end interface

  !---------------------------------------------------------------------------!
  abstract interface
     subroutine eval_multiple(self, position, components, field_dofs, field_value)
       use sll_m_working_precision
       import sll_c_particle_mesh_coupling_3d
       class (sll_c_particle_mesh_coupling_3d), intent( inout ) :: self !< Kernel smoother object 
       sll_real64,                    intent( in )    :: position(3) !< Position of the particle
       sll_int32,                     intent(in)      :: components(:)
       sll_real64,                    intent( in )    :: field_dofs(:,:) !< Coefficient vector for the field DoFs
       sll_real64,                    intent(out)     :: field_value(:) !< Value(s) of the electric fields at given position
     end subroutine eval_multiple
  end interface

  !---------------------------------------------------------------------------!
  abstract interface
     subroutine add_current(self, position_old, position_new, xdot, j_dofs) 
       use sll_m_working_precision
       import sll_c_particle_mesh_coupling_3d
       class (sll_c_particle_mesh_coupling_3d), intent( inout ) :: self !< Kernel smoother object
       sll_real64,                    intent( in )    :: position_old(3) !< Position of the particle
       sll_real64,                    intent( in )    :: position_new(3) !< Position of the particle
       sll_real64,                    intent( in )    :: xdot(3) !< Particle weight times charge times velocity
       sll_real64,                    intent( inout ) :: j_dofs(:) !< Coefficient vector of the charge distribution

     end subroutine add_current
  end interface

  !---------------------------------------------------------------------------!
  abstract interface
     subroutine add_current_evaluate(self, position_old, position_new, xdot, efield_dofs, j_dofs, efield_val) 
       use sll_m_working_precision
       import sll_c_particle_mesh_coupling_3d
       class (sll_c_particle_mesh_coupling_3d), intent( inout ) :: self !< Kernel smoother object
       sll_real64,                               intent( in )    :: position_old(3) !< Position of the particle
       sll_real64,                               intent( in )    :: position_new(3) !< Position of the particle
       sll_real64,                               intent( in )    :: xdot(3) !< velocity
       sll_real64,                               intent( in )    :: efield_dofs(:) !< electric field dofs
       sll_real64,                               intent( inout ) :: j_dofs(:) !< current dofs
       sll_real64,                               intent( out   ) :: efield_val(3) !< electric field value

     end subroutine add_current_evaluate
  end interface

  !---------------------------------------------------------------------------!
  abstract interface
     subroutine add_update (self, position_old, position_new, marker_charge, &
          qoverm, bfield_dofs, vi, j_dofs)
       use sll_m_working_precision
       import sll_c_particle_mesh_coupling_3d
       class(sll_c_particle_mesh_coupling_3d), intent(inout) :: self !< kernel smoother object
       sll_real64, intent(in)    :: position_old(3) !< Position at time t
       sll_real64, intent(in)    :: position_new !< Position at time t+\Delta t
       sll_real64, intent(in)    :: marker_charge          !< Particle weight times charge
       sll_real64, intent(in)    :: qoverm   !< charge to mass ratio
       sll_real64, intent(in)    :: bfield_dofs(:) !< values of the B-field at the dofs
       sll_real64, intent(inout) :: vi(3) !< Velocity of the particle
       sll_real64, intent(inout) :: j_dofs(:) !< Current at the DoFs

     end subroutine add_update
  end interface

  !---------------------------------------------------------------------------!  
  abstract interface
     subroutine empty(self)
       import sll_c_particle_mesh_coupling_3d
       class (sll_c_particle_mesh_coupling_3d), intent( inout ) :: self !< Kernel smoother object 

     end subroutine empty
  end interface

end module sll_m_particle_mesh_coupling_base_3d
