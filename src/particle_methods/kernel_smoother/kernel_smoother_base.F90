!> @ingroup particle_methods
!> @author Katharina Kormann
!> @brief Base class for kernel smoothers for accumulation and field evaluation in PIC.
!> @details This base class gives an abstract interface to the basic functions for accumulation of charge and current densities as well as the evaluation of a function at particle positions.
module sll_m_kernel_smoother_base

#include "sll_working_precision.h"

  use sll_module_pic_base

  implicit none
  private

  !> Basic type of a kernel smoother used for PIC simulations
  type, public, abstract :: sll_kernel_smoother_base
     sll_int32              :: n_dofs  !< Number of degrees of freedom of the smoothing kernels.
     sll_int32, allocatable :: n_grid(:) !< Number of grid points per dimension for use on tensor product grid based smoothing kernels.
     
   contains
     procedure(update_this), deferred           :: compute_shape_factors !< Prepare for the accumulation by computing the shape factors
     procedure(update_dofs), deferred           :: accumulate_rho_from_klimontovich !< Accumulate the charge density
     procedure(update_dofs_component), deferred :: accumulate_j_from_klimontovich !< Accumulate the current density
     procedure(evaluate), deferred              :: evaluate_kernel_function !< Evaluate a function based on the given degrees of freedom
  end type sll_kernel_smoother_base

!---------------------------------------------------------------------------!
  abstract interface
     subroutine update_this(this, particle_group)
       use sll_working_precision
       import sll_particle_group_base
       import sll_kernel_smoother_base
       class( sll_kernel_smoother_base), intent(inout) :: this !< Kernel smoother object.
       class( sll_particle_group_base), intent(in)     :: particle_group !< Particle group object.
     end subroutine update_this
  end interface
  
!---------------------------------------------------------------------------!
  abstract interface
     subroutine update_dofs(this, particle_group, rho_dofs)       
       use sll_working_precision
       import sll_particle_group_base
       import sll_kernel_smoother_base
       class( sll_kernel_smoother_base), intent(in)    :: this !< Kernel smoother object.
       class( sll_particle_group_base), intent(in)     :: particle_group !< Particle group object.
       sll_real64, intent(inout)                       :: rho_dofs(:) !< Degrees of freedom in kernel representation (can be point values or weights in a basis function representation).
     end subroutine update_dofs
  end interface

!---------------------------------------------------------------------------!
  abstract interface
     subroutine update_dofs_component(this, &
          particle_group, &
          j_dofs, &
          component)       
       use sll_working_precision
       import sll_particle_group_base
       import sll_kernel_smoother_base
       class( sll_kernel_smoother_base), intent(in)    :: this !< Kernel smoother object.
       class( sll_particle_group_base), intent(in)     :: particle_group !< Particle group object.
       sll_real64, intent(inout)                       :: j_dofs(:)!< Degrees of freedom in kernel representation (can be point values or weights in a basis function representation).
       sll_int32, intent (in)                          :: component !< Component of the current density that should be evaluated.
     end subroutine update_dofs_component
  end interface

!---------------------------------------------------------------------------!
  abstract interface
     subroutine evaluate(this, particle_group, rho_dofs, particle_values)       
       use sll_working_precision
       import sll_particle_group_base
       import sll_kernel_smoother_base
       class( sll_kernel_smoother_base), intent(in)    :: this !< Kernel smoother object.
       class( sll_particle_group_base), intent(in)     :: particle_group !< Particle group object.
       sll_real64, intent(in)                       :: rho_dofs(:) !< Degrees of freedom in kernel representation.
       sll_real64, intent(out)                      :: particle_values(:) !< Values of the function represented by \a rho_dofs at particle positions.
     end subroutine evaluate
  end interface


end module sll_m_kernel_smoother_base
