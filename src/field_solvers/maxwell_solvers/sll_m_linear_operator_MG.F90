module sll_m_linear_operator_MG
#include "sll_working_precision.h"
  
  use sll_m_linear_operator_abstract, only : &
       sll_t_linear_operator_abstract

  use sll_m_linear_operator_block, only : &
       sll_t_linear_operator_block

  use sll_m_spline_fem_utilities, only : &
       sll_s_multiply_g

  implicit none

  public :: sll_t_linear_operator_MG

  private

  type, extends(sll_t_linear_operator_abstract) :: sll_t_linear_operator_MG
     type(sll_t_linear_operator_block), pointer :: mass   !< block mass matrix
     sll_int32  :: n_dofs(3)  !< number of degrees of freedom
     sll_real64 :: delta_x(3) !< cell size
     sll_int32  :: n_total !< product of number of degrees of freedom


   contains
     procedure :: create => create_linear_operator_MG
     procedure :: free => free_MG
     procedure :: dot => dot_MG
     procedure :: print_info => print_info_MG

  end type sll_t_linear_operator_MG

contains


  subroutine create_linear_operator_MG( self, mass, n_dofs, delta_x)
    class(sll_t_linear_operator_MG), intent( inout ) :: self !< Linear operator object
    type(sll_t_linear_operator_block), target :: mass   !< block mass matrix
    sll_int32                                 :: n_dofs(3) !< number of degrees of freedom
    sll_real64                                :: delta_x(3) !< cell size

    self%mass    => mass
    self%n_dofs  = n_dofs
    self%n_total = product(n_dofs)
    self%delta_x = delta_x

    self%n_rows  = 3*self%n_total
    self%n_cols  = self%n_total

    self%n_global_rows = self%n_rows
    self%n_global_cols = self%n_cols
  end subroutine create_linear_operator_MG


  subroutine free_MG( self )
    class(sll_t_linear_operator_MG), intent( inout ) :: self !< Linear operator object
    self%mass => null()

  end subroutine free_MG


  subroutine dot_MG( self, x, y )
    class(sll_t_linear_operator_MG), intent( in ) :: self !< Linear operator object
    sll_real64, intent( in    ) :: x(:) !< Inputvariable
    sll_real64, intent(   out ) :: y(:) !< Outputvariable
    !local variable
    sll_real64 :: scratch(3*self%n_total)

    call sll_s_multiply_g(self%n_dofs, self%delta_x, x, scratch)
    call self%mass%dot(scratch, y)

  end subroutine dot_MG

  subroutine print_info_MG( self )
    class(sll_t_linear_operator_MG), intent(in) :: self !< Linear operator object
  end subroutine print_info_MG


end module sll_m_linear_operator_MG
