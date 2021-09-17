module sll_m_linear_operator_poisson_clamped_1d
#include "sll_working_precision.h"
  
  use sll_m_linear_operator_abstract, only : &
       sll_t_linear_operator_abstract

  use sll_m_matrix_csr, only: &
       sll_t_matrix_csr

  use sll_m_spline_fem_utilities_3d_clamped, only: &
       sll_s_multiply_g_clamped_1d, &
       sll_s_multiply_gt_clamped_1d
  
  implicit none

  public :: sll_t_linear_operator_poisson_clamped_1d

  private

  type, extends(sll_t_linear_operator_abstract) :: sll_t_linear_operator_poisson_clamped_1d
     type(sll_t_matrix_csr), pointer :: mass   !< block mass matrix
     sll_int32  :: s_deg_0 !< spline degree 0-forms
     sll_int32  :: n_dofs  !< number of degrees of freedom
     sll_real64 :: delta_x !< cell size
     
   contains
     procedure :: create     => create_linear_operator_poisson_clamped_1d
     procedure :: free       => free_poisson_clamped_1d
     procedure :: dot        => dot_poisson_clamped_1d
     procedure :: print_info => print_info_poisson_clamped_1d

  end type sll_t_linear_operator_poisson_clamped_1d


contains


  subroutine create_linear_operator_poisson_clamped_1d( self, mass, s_deg_0, n_dofs, delta_x)
    class(sll_t_linear_operator_poisson_clamped_1d), intent( inout ) :: self !< Linear operator object
    type(sll_t_matrix_csr), target   :: mass   !< mass matrix
    sll_int32, intent(in)  :: s_deg_0 !< spline degree 0-forms
    sll_int32, intent(in)  :: n_dofs  !< number of degrees of freedom
    sll_real64, intent(in) :: delta_x !< cell size

    self%mass    => mass
    self%s_deg_0 = s_deg_0
    self%n_dofs  = n_dofs
    self%delta_x = delta_x

    self%n_rows  = self%n_dofs
    self%n_cols  = self%n_dofs

    self%n_global_rows = self%n_rows
    self%n_global_cols = self%n_cols
  end subroutine create_linear_operator_poisson_clamped_1d


  subroutine free_poisson_clamped_1d( self )
    class(sll_t_linear_operator_poisson_clamped_1d), intent( inout ) :: self !< Linear operator object
    self%mass => null()

  end subroutine free_poisson_clamped_1d


  subroutine dot_poisson_clamped_1d( self, x, y )
    class(sll_t_linear_operator_poisson_clamped_1d), intent( in ) :: self !< Linear operator object
    sll_real64, intent( in    ) :: x(:) !< Inputvariable
    sll_real64, intent(   out ) :: y(:) !< Outputvariable
    !local variable
    sll_real64 :: scratch1(self%n_dofs-1), scratch2(self%n_dofs-1)

    call sll_s_multiply_g_clamped_1d(self%n_dofs, self%delta_x, self%s_deg_0, x, scratch1)
    call self%mass%dot(scratch1, scratch2)
    call sll_s_multiply_gt_clamped_1d(self%n_dofs, self%delta_x, self%s_deg_0, scratch2, y)
  end subroutine dot_poisson_clamped_1d


  subroutine print_info_poisson_clamped_1d( self )
    class(sll_t_linear_operator_poisson_clamped_1d), intent(in) :: self !< Linear operator object
  end subroutine print_info_poisson_clamped_1d

end module sll_m_linear_operator_poisson_clamped_1d
