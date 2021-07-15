module sll_m_linear_operator_schur_eb_cl_1d
#include "sll_working_precision.h"
  
  use sll_m_linear_operator_abstract

  use sll_m_matrix_csr, only: &
       sll_t_matrix_csr

  use sll_m_spline_fem_utilities_3d_clamped, only: &
       sll_s_multiply_g_clamped_1d, &
       sll_s_multiply_gt_clamped_1d

 
  implicit none

  public :: sll_t_linear_operator_schur_eb_cl_1d

  private
  type, extends(sll_t_linear_operator_abstract) :: sll_t_linear_operator_schur_eb_cl_1d
     
     type(sll_t_matrix_csr), pointer :: mass0 !< spline mass matrix for 0-form
     type(sll_t_matrix_csr), pointer :: mass1 !< spline mass matrix for 1-form
     sll_int32 :: n_dofs !< number of degrees of freedom
     sll_real64 :: delta_x !< cell size
     sll_real64 :: sign = 1.0_f64 !< sign
     sll_int32  :: s_deg_0 !< spline degree 0-forms 
     
   contains
     procedure :: create => create_linear_operator_schur_eb_cl_1d
     procedure :: free => free_schur_eb_cl_1d
     procedure :: dot => dot_schur_eb_cl_1d
     procedure :: print_info => print_info_schur_eb_cl_1d

  end type sll_t_linear_operator_schur_eb_cl_1d


contains

  subroutine create_linear_operator_schur_eb_cl_1d( self, mass0, mass1, n_dofs, delta_x, s_deg_0 )
    class(sll_t_linear_operator_schur_eb_cl_1d), intent( inout ) :: self !< Linear operator object
    type(sll_t_matrix_csr), target :: mass0 !< spline mass matrix for 0-form
    type(sll_t_matrix_csr), target :: mass1 !< spline mass matrix for 1-form
    sll_int32  :: n_dofs !< number of degrees of freedom
    sll_real64 :: delta_x !< cell size
    sll_int32  :: s_deg_0 !< spline degree 0-forms 

    self%mass0  => mass0
    self%mass1  => mass1
    self%n_dofs = n_dofs
    self%delta_x = delta_x
    self%s_deg_0 = s_deg_0

    self%n_rows = self%n_dofs
    self%n_cols = self%n_dofs

    self%n_global_rows = self%n_rows
    self%n_global_cols = self%n_cols

  end subroutine create_linear_operator_schur_eb_cl_1d

  subroutine free_schur_eb_cl_1d( self )
    class(sll_t_linear_operator_schur_eb_cl_1d), intent( inout ) :: self !< Linear operator object

    self%mass0=> null()
    self%mass1=> null()

  end subroutine free_schur_eb_cl_1d


  subroutine dot_schur_eb_cl_1d ( self, x, y )
    class(sll_t_linear_operator_schur_eb_cl_1d), intent( in ) :: self !< Linear operator object
    sll_real64, intent( in    ) :: x(:) !< Inputvariable
    sll_real64, intent(   out ) :: y(:) !< Outputvariable

    !local variables
    sll_real64 ::  scratch0(self%n_dofs), scratch1(self%n_dofs-1), scratch2(self%n_dofs-1)

    ! Compute D x
    call sll_s_multiply_g_clamped_1d( self%n_dofs, self%delta_x, self%s_deg_0, x, scratch1)
    
    ! Compute M1 D x
    call self%mass1%dot( scratch1, scratch2 )

    ! Compute D^T M1 D x + M_b D x
    call sll_s_multiply_gt_clamped_1d( self%n_dofs, self%delta_x, self%s_deg_0, scratch2, scratch0 )

    !scratch0(1)= scratch0(1)+scratch1(1)
    !scratch0(self%n_dofs)= scratch0(self%n_dofs)-scratch1(self%n_dofs-1)

    ! Compute M0 x
    call self%mass0%dot(  x, y )

    ! Sum up the two parts
    y = y + self%sign * scratch0

  end subroutine dot_schur_eb_cl_1d

  subroutine print_info_schur_eb_cl_1d( self )
    class(sll_t_linear_operator_schur_eb_cl_1d), intent(in) :: self !< Linear operator object

  end subroutine print_info_schur_eb_cl_1d


end module sll_m_linear_operator_schur_eb_cl_1d
