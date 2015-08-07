!> @ingroup particle_methods
!> @brief
!> Interface for combined Maxwell solver and kernel smoother for PIC methods
!> @details
!> Derived types can be implemented directly or initialized using a Maxwell solver and a kernel smoother.

module sll_m_pic_maxwell_base
#include "sll_working_precision.h"

  use sll_module_pic_base

  implicit none
  private
  
  type, public, abstract :: sll_pic_maxwell_base
     sll_int32              :: n_dofs  !< Number of degrees of freedom of the smoothing kernels.
     sll_int32, allocatable :: n_grid(:) !< Number of grid points per dimension for use on tensor product grid based smoothing kernels.
     

   contains
     ! Maxwell solver functions
     procedure(compute_field1_from_field2_pic), deferred :: &
          compute_E_from_B !< Solve E and B part of Amperes law with B constant in time
     procedure(compute_field1_from_field2_pic), deferred :: &
          compute_B_from_E !< Solve Faraday equation with E constant in time
     procedure(signature_compute_E_from_rho_pic), deferred :: &
          compute_E_from_rho !< Solve E from rho using Poisson
     !procedure(signature_solve), deferred :: &
     !     solve !< Solve Amperes law and Faraday equation
     ! Kernel smoother functions
     procedure(update_this_pic), deferred           :: compute_shape_factors !< Prepare for the accumulation by computing the shape factors
     procedure(update_dofs_pic), deferred           :: accumulate_rho_from_klimontovich !< Accumulate the charge density
     procedure(update_dofs_component_pic), deferred :: accumulate_j_from_klimontovich !< Accumulate the current density
     procedure(evaluate_pic), deferred              :: evaluate_kernel_function !< Evaluate a function based on the given degrees of freedom
     

  end type sll_pic_maxwell_base

  ! Interfaces for Maxwell
  abstract interface 
     subroutine compute_field1_from_field2_pic(this, delta_t, field_in, field_out)
     use sll_working_precision
     import sll_pic_maxwell_base
     
     class(sll_pic_maxwell_base) :: this
     sll_real64, intent(in)     :: delta_t
     sll_real64, intent(in)     :: field_in(:)
     sll_real64, intent(inout)  :: field_out(:)
   end subroutine compute_field1_from_field2_pic
  end interface

  abstract interface    
    subroutine signature_compute_E_from_rho_pic(this, E, rho )
      use sll_working_precision
      import sll_pic_maxwell_base       
      class(sll_pic_maxwell_base) :: this
      sll_real64,dimension(:),intent(in) :: rho
      sll_real64,dimension(:),intent(out) :: E
    end subroutine signature_compute_E_from_rho_pic
  end interface

  ! Interfaces for kernel smoother
  !---------------------------------------------------------------------------!
  abstract interface
     subroutine update_this_pic(this, particle_group)
       use sll_working_precision
       import sll_particle_group_base
       import sll_pic_maxwell_base
       class( sll_pic_maxwell_base), intent(inout) :: this !< Kernel smoother object.
       class( sll_particle_group_base), intent(in)     :: particle_group !< Particle group object.
     end subroutine update_this_pic
  end interface
  
!---------------------------------------------------------------------------!
  abstract interface
     subroutine update_dofs_pic(this, particle_group, rho_dofs)       
       use sll_working_precision
       import sll_particle_group_base
       import sll_pic_maxwell_base
       class( sll_pic_maxwell_base), intent(in)    :: this !< Kernel smoother object.
       class( sll_particle_group_base), intent(in)     :: particle_group !< Particle group object.
       sll_real64, intent(inout)                       :: rho_dofs(:) !< Degrees of freedom in kernel representation (can be point values or weights in a basis function representation).
     end subroutine update_dofs_pic
  end interface
!!$
!---------------------------------------------------------------------------!
  abstract interface
     subroutine update_dofs_component_pic(this, &
          particle_group, &
          j_dofs, &
          component)       
       use sll_working_precision
       import sll_particle_group_base
       import sll_pic_maxwell_base
       class( sll_pic_maxwell_base), intent(in)    :: this !< Kernel smoother object.
       class( sll_particle_group_base), intent(in)     :: particle_group !< Particle group object.
       sll_real64, intent(inout)                       :: j_dofs(:)!< Degrees of freedom in kernel representation (can be point values or weights in a basis function representation).
       sll_int32, intent (in)                          :: component !< Component of the current density that should be evaluated.
     end subroutine update_dofs_component_pic
  end interface

!---------------------------------------------------------------------------!
  abstract interface
     subroutine evaluate_pic(this, particle_group, rho_dofs, particle_values)       
       use sll_working_precision
       import sll_particle_group_base
       import sll_pic_maxwell_base
       class( sll_pic_maxwell_base), intent(in)    :: this !< Kernel smoother object.
       class( sll_particle_group_base), intent(in)     :: particle_group !< Particle group object.
       sll_real64, intent(in)                       :: rho_dofs(:) !< Degrees of freedom in kernel representation.
       sll_real64, intent(out)                      :: particle_values(:) !< Values of the function represented by \a rho_dofs at particle positions.
     end subroutine evaluate_pic
  end interface


end module sll_m_pic_maxwell_base
