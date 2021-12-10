module sll_m_linear_operator_GTM_cl
#include "sll_working_precision.h"

  use sll_m_linear_operator_abstract, only : &
       sll_t_linear_operator_abstract

  use sll_m_linear_operator_block, only : &
       sll_t_linear_operator_block

  use sll_m_spline_fem_utilities_3d_clamped, only : &
       sll_s_multiply_gt_clamped

  implicit none

  public :: sll_t_linear_operator_GTM_cl

  private

  type, extends(sll_t_linear_operator_abstract) :: sll_t_linear_operator_GTM_cl
     type(sll_t_linear_operator_block), pointer :: mass   !< block mass matrix
     sll_int32  :: n_dofs(3)  !< number of degrees of freedom
     sll_real64 :: delta_x(3) !< cell size
     sll_int32 :: s_deg_0(3) !< spline degree 0-forms 
     sll_int32 :: n_total0, n_total1  !< product of number of degrees of freedom


   contains
     procedure :: create => create_linear_operator_GTM_cl
     procedure :: free => free_GTM_cl
     procedure :: dot => dot_GTM_cl
     procedure :: print_info => print_info_GTM_cl

  end type sll_t_linear_operator_GTM_cl

contains


  subroutine create_linear_operator_GTM_cl( self, mass, n_dofs, delta_x, s_deg_0 )
    class(sll_t_linear_operator_GTM_cl), intent( inout ) :: self !< Linear operator object
    type(sll_t_linear_operator_block), target :: mass   !< block mass matrix
    sll_int32                                 :: n_dofs(3) !< number of degrees of freedom
    sll_real64                                :: delta_x(3) !< cell size
    sll_int32                                 :: s_deg_0(3) !< spline degree 0-forms 

    self%mass    => mass
    self%n_dofs  = n_dofs
    self%s_deg_0 = s_deg_0
    self%n_total0 = product(self%n_dofs)
    self%n_total1 = (self%n_dofs(1)-1)*self%n_dofs(2)*self%n_dofs(3)
    self%delta_x = delta_x

    self%n_rows  = self%n_total0
    self%n_cols  = self%n_total1+2*self%n_total0

    self%n_global_rows = self%n_rows
    self%n_global_cols = self%n_cols
  end subroutine create_linear_operator_GTM_cl


  subroutine free_GTM_cl( self )
    class(sll_t_linear_operator_GTM_cl), intent( inout ) :: self !< Linear operator object
    self%mass => null()

  end subroutine free_GTM_cl


  subroutine dot_GTM_cl( self, x, y )
    class(sll_t_linear_operator_GTM_cl), intent( in ) :: self !< Linear operator object
    sll_real64, intent( in    ) :: x(:) !< Inputvariable
    sll_real64, intent(   out ) :: y(:) !< Outputvariable
    !local variable
    sll_real64 :: scratch(self%n_total1+2*self%n_total0)

    call self%mass%dot(x, scratch)
    call sll_s_multiply_gt_clamped(self%n_dofs, self%delta_x, self%s_deg_0, scratch, y)

  end subroutine dot_GTM_cl

  subroutine print_info_GTM_cl( self )
    class(sll_t_linear_operator_GTM_cl), intent(in) :: self !< Linear operator object
  end subroutine print_info_GTM_cl


end module sll_m_linear_operator_GTM_cl
