module sll_m_linear_operator_GTM
#include "sll_working_precision.h"
  
  use sll_m_linear_operator_abstract, only : &
       sll_t_linear_operator_abstract

  use sll_m_linear_operator_block, only : &
       sll_t_linear_operator_block

  use sll_m_spline_fem_utilities, only : &
       sll_s_multiply_gt

  implicit none

  public :: sll_t_linear_operator_GTM

  private

  type, extends(sll_t_linear_operator_abstract) :: sll_t_linear_operator_GTM
     type(sll_t_linear_operator_block), pointer :: mass   !< block mass matrix
     sll_int32  :: n_dofs(3)  !< number of degrees of freedom
     sll_real64 :: delta_x(3) !< cell size
     sll_int32  :: n_total !< product of number of degrees of freedom


   contains
     procedure :: create => create_linear_operator_GTM
     procedure :: free => free_GTM
     procedure :: dot => dot_GTM
     procedure :: print_info => print_info_GTM

  end type sll_t_linear_operator_GTM

contains


  subroutine create_linear_operator_GTM( self, mass, n_dofs, delta_x)
    class(sll_t_linear_operator_GTM), intent( inout ) :: self !< Linear operator object
    type(sll_t_linear_operator_block), target :: mass   !< block mass matrix
    sll_int32                                 :: n_dofs(3) !< number of degrees of freedom
    sll_real64                                :: delta_x(3) !< cell size

    self%mass    => mass
    self%n_dofs  = n_dofs
    self%n_total = product(n_dofs)
    self%delta_x = delta_x

    self%n_rows  = self%n_total
    self%n_cols  = 3*self%n_total

    self%n_global_rows = self%n_rows
    self%n_global_cols = self%n_cols
  end subroutine create_linear_operator_GTM


  subroutine free_GTM( self )
    class(sll_t_linear_operator_GTM), intent( inout ) :: self !< Linear operator object
    self%mass => null()

  end subroutine free_GTM


  subroutine dot_GTM( self, x, y )
    class(sll_t_linear_operator_GTM), intent( in ) :: self !< Linear operator object
    sll_real64, intent( in    ) :: x(:) !< Inputvariable
    sll_real64, intent(   out ) :: y(:) !< Outputvariable
    !local variable
    sll_real64 :: scratch(3*self%n_total)
  
    call self%mass%dot(x, scratch)
    call sll_s_multiply_gt(self%n_dofs, self%delta_x, scratch, y)

  end subroutine dot_GTM

  subroutine print_info_GTM( self )
    class(sll_t_linear_operator_GTM), intent(in) :: self !< Linear operator object
  end subroutine print_info_GTM


end module sll_m_linear_operator_GTM
