module sll_m_linear_operator_poisson_3d
#include "sll_working_precision.h"
  
  use sll_m_linear_operator_abstract, only : &
       sll_t_linear_operator_abstract

  use sll_m_linear_operator_block, only : &
       sll_t_linear_operator_block

  use sll_m_spline_fem_utilities, only : &
       sll_s_multiply_g, &
       sll_s_multiply_gt

  implicit none

  public :: sll_t_linear_operator_poisson_3d

  private

  type, extends(sll_t_linear_operator_abstract) :: sll_t_linear_operator_poisson_3d
     type(sll_t_linear_operator_block), pointer :: mass   !< block mass matrix
     sll_int32  :: n_dofs(3)  !< number of degrees of freedom
     sll_real64 :: delta_x(3) !< cell size
     sll_int32  :: n_total !< product of number of degrees of freedom


   contains
     procedure :: create => create_linear_operator_poisson_3d
     procedure :: free => free_poisson_3d
     procedure :: dot => dot_poisson_3d
     procedure :: print_info => print_info_poisson_3d

  end type sll_t_linear_operator_poisson_3d

contains


  subroutine create_linear_operator_poisson_3d( self, mass, n_dofs, delta_x)
    class(sll_t_linear_operator_poisson_3d), intent( inout ) :: self !< Linear operator object
    type(sll_t_linear_operator_block), target :: mass   !< block mass matrix
    sll_int32                                 :: n_dofs(3) !< number of degrees of freedom
    sll_real64                                :: delta_x(3) !< cell size

    self%mass    => mass
    self%n_dofs  = n_dofs
    self%n_total = product(n_dofs)
    self%delta_x = delta_x

    self%n_rows  = self%n_total
    self%n_cols  = self%n_total

    self%n_global_rows = self%n_rows
    self%n_global_cols = self%n_cols
  end subroutine create_linear_operator_poisson_3d


  subroutine free_poisson_3d( self )
    class(sll_t_linear_operator_poisson_3d), intent( inout ) :: self !< Linear operator object
    self%mass => null()

  end subroutine free_poisson_3d


  subroutine dot_poisson_3d( self, x, y )
    class(sll_t_linear_operator_poisson_3d), intent( in ) :: self !< Linear operator object
    sll_real64, intent( in    ) :: x(:) !< Inputvariable
    sll_real64, intent(   out ) :: y(:) !< Outputvariable
    !local variable
    sll_real64 :: scratch(3*self%n_total),  scratch1(3*self%n_total)

    call sll_s_multiply_g(self%n_dofs, self%delta_x, x, scratch)
    call self%mass%dot(scratch, scratch1)
    call sll_s_multiply_gt(self%n_dofs, self%delta_x, scratch1, y)

  end subroutine dot_poisson_3d

  subroutine print_info_poisson_3d( self )
    class(sll_t_linear_operator_poisson_3d), intent(in) :: self !< Linear operator object
  end subroutine print_info_poisson_3d


end module sll_m_linear_operator_poisson_3d
