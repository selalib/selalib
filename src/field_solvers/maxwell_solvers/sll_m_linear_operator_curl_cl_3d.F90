module sll_m_linear_operator_curl_cl_3d
#include "sll_working_precision.h"

  use sll_m_linear_operator_abstract

  use sll_m_spline_fem_utilities_3d_clamped, only : &
       sll_s_multiply_g_clamped, &
       sll_s_multiply_gt_clamped, &
       sll_s_multiply_c_clamped, &
       sll_s_multiply_ct_clamped

  use sll_m_linear_operator_block, only : &
       sll_t_linear_operator_block

  implicit none

  public :: sll_t_linear_operator_curl_cl_3d

  private
  type, extends(sll_t_linear_operator_abstract) :: sll_t_linear_operator_curl_cl_3d
     type(sll_t_linear_operator_block), pointer :: mass1   !< block mass matrix
     type(sll_t_linear_operator_block), pointer :: mass2   !< block mass matrix
     sll_int32 :: n_total0, n_total1  !< product of number of degrees of freedom
     sll_int32 :: n_dofs(3) !< number of degrees of freedom
     sll_int32 :: s_deg_0(3) !< spline degree 0-forms 
     sll_real64 :: delta_x(3) !< cell size
     sll_real64 :: epsilon = 1._f64
    

   contains
     procedure :: create => create_linear_operator_curl_cl_3d
     procedure :: free => free_curl_cl_3d
     procedure :: dot => dot_curl_cl_3d
     procedure :: print_info => print_info_curl_cl_3d

  end type sll_t_linear_operator_curl_cl_3d


contains

  subroutine create_linear_operator_curl_cl_3d( self, mass1, mass2, n_dofs, delta_x, s_deg_0 )
    class(sll_t_linear_operator_curl_cl_3d), intent( inout ) :: self !< Linear operator object
    type(sll_t_linear_operator_block), target :: mass1   !< block mass matrix
    type(sll_t_linear_operator_block), target :: mass2   !< block mass matrix
    sll_int32  :: n_dofs(3) !< number of degrees of freedom
    sll_real64 :: delta_x(3) !< cell size
    sll_int32  ::  s_deg_0(3) !< spline degree 0-forms 

    self%mass1 => mass1
    self%mass2 => mass2

    self%n_dofs = n_dofs
    self%delta_x = delta_x
    self%s_deg_0 = s_deg_0
    self%n_total0 = product(self%n_dofs)
    self%n_total1 = (self%n_dofs(1)-1)*self%n_dofs(2)*self%n_dofs(3)

    self%n_rows = self%n_total1+2*self%n_total0
    self%n_cols = self%n_total1+2*self%n_total0

    self%n_global_rows = self%n_rows
    self%n_global_cols = self%n_cols

  end subroutine create_linear_operator_curl_cl_3d

  subroutine free_curl_cl_3d( self )
    class(sll_t_linear_operator_curl_cl_3d), intent( inout ) :: self !< Linear operator object

    self%mass1=> null()
    self%mass2=> null()

  end subroutine free_curl_cl_3d


  subroutine dot_curl_cl_3d ( self, x, y )
    class(sll_t_linear_operator_curl_cl_3d), intent( in ) :: self !< Linear operator object
    sll_real64, intent( in    ) :: x(:) !< Inputvariable
    sll_real64, intent(   out ) :: y(:) !< Outputvariable
    !local variables
     sll_real64 :: scratch0(self%n_total1+2*self%n_total0), scratch1(self%n_total1+2*self%n_total0), scratch2(self%n_total0+2*self%n_total1), scratch3(self%n_total0+2*self%n_total1)

    ! Compute C x
    call sll_s_multiply_c_clamped( self%n_dofs, self%delta_x, self%s_deg_0, x, scratch2 )

    ! Compute M2 C x
    call self%mass2%dot( scratch2, scratch3  )

    ! Compute C^T M2 C x
    call sll_s_multiply_ct_clamped( self%n_dofs, self%delta_x, self%s_deg_0, scratch3, y )
!!$
    ! Compute M1 x
    call self%mass1%dot( x, scratch1 )

    ! Compute G^T M1 x
    call sll_s_multiply_gt_clamped( self%n_dofs, self%delta_x, self%s_deg_0, scratch1, scratch2(1:self%n_total0) )

     ! Compute G G^T M1 x
    call sll_s_multiply_g_clamped( self%n_dofs, self%delta_x, self%s_deg_0, scratch2(1:self%n_total0), scratch1 )

    ! Compute M1 G G^T M_1 x
    call self%mass1%dot( scratch1, scratch0 )

    ! y = (C^T M2 C + epsilon M1 G G^T M_1) x
    y = y + self%epsilon*scratch0
    
  end subroutine dot_curl_cl_3d


  subroutine print_info_curl_cl_3d( self )
    class(sll_t_linear_operator_curl_cl_3d), intent(in) :: self !< Linear operator object

  end subroutine print_info_curl_cl_3d


end module sll_m_linear_operator_curl_cl_3d
