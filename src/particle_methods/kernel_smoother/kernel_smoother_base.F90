module sll_m_kernel_smoother_base

#include "sll_working_precision.h"
! TODO: Klaere passendes Format fuer dofs und particle_values. Might be abstract to be dimension independent or array but than tensor product structure and dimension dependent.

  use sll_module_pic_base

  implicit none
  private

  type, public, abstract :: sll_kernel_smoother_base
     sll_int32 :: n_dofs
     sll_int32, allocatable :: n_grid(:)
     
   contains
     procedure(update_this), deferred           :: compute_shape_factors
     procedure(update_dofs), deferred           :: accumulate_rho_from_klimontovich
     procedure(update_dofs_component), deferred :: accumulate_j_from_klimontovich
     procedure(evaluate), deferred              :: evaluate_kernel_function 
  end type sll_kernel_smoother_base

!---------------------------------------------------------------------------!
  abstract interface
     subroutine update_this(this, particle_group)
       use sll_working_precision
       import sll_particle_group_base
       import sll_kernel_smoother_base
       class( sll_kernel_smoother_base), intent(inout) :: this
       class( sll_particle_group_base), intent(in)     :: particle_group
     end subroutine update_this
  end interface
  
!---------------------------------------------------------------------------!
  abstract interface
     subroutine update_dofs(this, particle_group, rho_dofs)       
       use sll_working_precision
       import sll_particle_group_base
       import sll_kernel_smoother_base
       class( sll_kernel_smoother_base), intent(in)    :: this
       class( sll_particle_group_base), intent(in)     :: particle_group
       sll_real64, intent(inout)                       :: rho_dofs(:)
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
       class( sll_kernel_smoother_base), intent(in)    :: this
       class( sll_particle_group_base), intent(in)     :: particle_group
       sll_real64, intent(inout)                       :: j_dofs(:)
       sll_int32, intent (in)                          :: component
     end subroutine update_dofs_component
  end interface

!---------------------------------------------------------------------------!
  abstract interface
     subroutine evaluate(this, particle_group, rho_dofs, particle_values)       
       use sll_working_precision
       import sll_particle_group_base
       import sll_kernel_smoother_base
       class( sll_kernel_smoother_base), intent(in)    :: this
       class( sll_particle_group_base), intent(in)     :: particle_group
       sll_real64, intent(in)                       :: rho_dofs(:)
       sll_real64, intent(out)                      :: particle_values(:)
     end subroutine evaluate
  end interface


end module sll_m_kernel_smoother_base
