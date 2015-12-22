!> @ingroup kernel_smoothers
!> @author Katharina Kormann, IPP
!> @brief Base class for kernel smoothers for accumulation and field evaluation in PIC.
!> @details This base class gives an abstract interface to the basic functions for accumulation of charge and current densities as well as the evaluation of a function at particle positions.
module sll_m_kernel_smoother_base

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use sll_m_particle_group_base, only: &
    sll_c_particle_group_base

  implicit none

  public :: &
    sll_p_collocation, &
    sll_p_galerkin, &
    sll_c_kernel_smoother_base

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  ! Define parameters to set if Galerkin or collocation scaling should be used in accumulation routines
  sll_int32, parameter :: sll_p_galerkin = 0
  sll_int32, parameter :: sll_p_collocation = 1

  !> Basic type of a kernel smoother used for PIC simulations
  type, abstract :: sll_c_kernel_smoother_base
     sll_int32              :: n_dofs  !< Number of degrees of freedom of the smoothing kernels.
     sll_int32, allocatable :: n_grid(:) !< Number of grid points per dimension for use on tensor product grid based smoothing kernels.
     
   contains
     procedure(update_this), deferred           :: compute_shape_factors !< Prepare for the accumulation by computing the shape factors
     procedure(update_dofs), deferred           :: accumulate_rho_from_klimontovich !< Accumulate the charge density
     procedure(update_dofs_component), deferred :: accumulate_j_from_klimontovich !< Accumulate the current density
     procedure(evaluate_particle), deferred     :: evaluate_kernel_function_particle !< Evaluate function for a certain particle based on the precomputed shape factors
     procedure    :: evaluate_kernel_function_particles !< Evaluate function for all particle based on the precomputed shape factors
  end type sll_c_kernel_smoother_base

!---------------------------------------------------------------------------!
  abstract interface
     subroutine update_this(this, particle_group)
       use sll_m_working_precision
       import sll_c_particle_group_base
       import sll_c_kernel_smoother_base
       class( sll_c_kernel_smoother_base), intent(inout) :: this !< Kernel smoother object.
       class( sll_c_particle_group_base), intent(in)     :: particle_group !< Particle group object.
     end subroutine update_this
  end interface
  
!---------------------------------------------------------------------------!
  abstract interface
     subroutine update_dofs(this, particle_group, rho_dofs)       
       use sll_m_working_precision
       import sll_c_particle_group_base
       import sll_c_kernel_smoother_base
       class( sll_c_kernel_smoother_base), intent(in)    :: this !< Kernel smoother object.
       class( sll_c_particle_group_base), intent(in)     :: particle_group !< Particle group object.
       sll_real64, intent(inout)                       :: rho_dofs(:) !< Degrees of freedom in kernel representation (can be point values or weights in a basis function representation).
     end subroutine update_dofs
  end interface

!---------------------------------------------------------------------------!
  abstract interface
     subroutine update_dofs_component(this, &
          particle_group, &
          j_dofs, &
          component)       
       use sll_m_working_precision
       import sll_c_particle_group_base
       import sll_c_kernel_smoother_base
       class( sll_c_kernel_smoother_base), intent(in)    :: this !< Kernel smoother object.
       class( sll_c_particle_group_base), intent(in)     :: particle_group !< Particle group object.
       sll_real64, intent(inout)                       :: j_dofs(:)!< Degrees of freedom in kernel representation (can be point values or weights in a basis function representation).
       sll_int32, intent (in)                          :: component !< Component of the current density that should be evaluated.
     end subroutine update_dofs_component
  end interface

!---------------------------------------------------------------------------!
  abstract interface
     subroutine evaluate_particle(this, rho_dofs, i_part, particle_value)       
       use sll_m_working_precision
       import sll_c_kernel_smoother_base
       class( sll_c_kernel_smoother_base), intent(in) :: this !< Kernel smoother object.
       sll_real64, intent(in)                       :: rho_dofs(:) !< Degrees of freedom in kernel representation.
       sll_int32, intent(in)                        :: i_part !< particle number
       sll_real64, intent(out)                      :: particle_value !< Value of the function at the position of particle \a i_part
     end subroutine evaluate_particle
  end interface

contains
  
  !---------------------------------------------------------------------------!
  subroutine evaluate_kernel_function_particles(this, particle_group, rho_dofs, particle_values)       
    class( sll_c_kernel_smoother_base), intent(in)    :: this !< Kernel smoother object.
    class( sll_c_particle_group_base), intent(in)     :: particle_group !< Particle group object.
    sll_real64, intent(in)                       :: rho_dofs(:) !< Degrees of freedom in kernel representation.
    sll_real64, intent(out)                      :: particle_values(:) !< Values of the function represented by \a rho_dofs at particle positions.
    
    !local variables
    sll_int32 :: i_part
    
    do i_part = 1, particle_group%n_particles
       call this%evaluate_kernel_function_particle(rho_dofs, i_part, particle_values(i_part))
    end do
    
  end subroutine evaluate_kernel_function_particles



end module sll_m_kernel_smoother_base
