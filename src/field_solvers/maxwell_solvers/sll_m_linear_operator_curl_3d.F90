module sll_m_linear_operator_curl_3d
#include "sll_working_precision.h"

  use sll_m_linear_operator_abstract

  use sll_m_spline_fem_utilities, only : &
       sll_s_multiply_g, &
       sll_s_multiply_gt, &
       sll_s_multiply_c, &
       sll_s_multiply_ct

  use sll_m_linear_operator_block, only : &
       sll_t_linear_operator_block

  implicit none

  public :: sll_t_linear_operator_curl_3d

  private
  type, extends(sll_t_linear_operator_abstract) :: sll_t_linear_operator_curl_3d
     type(sll_t_linear_operator_block), pointer :: mass1   !< block mass matrix
     type(sll_t_linear_operator_block), pointer :: mass2   !< block mass matrix
     sll_int32 :: n_dofs(3) !< number of degrees of freedom
     sll_real64 :: delta_x(3) !< cell size
     sll_int32 :: n_total !< product of number of degrees of freedom

   contains
     procedure :: create => create_linear_operator_curl_3d
     procedure :: free => free_curl_3d
     procedure :: dot => dot_curl_3d
     procedure :: print_info => print_info_curl_3d

  end type sll_t_linear_operator_curl_3d


contains

  subroutine create_linear_operator_curl_3d( self, mass1, mass2, n_dofs, delta_x )
    class(sll_t_linear_operator_curl_3d), intent( inout ) :: self !< Linear operator object
    type(sll_t_linear_operator_block), target :: mass1   !< block mass matrix
    type(sll_t_linear_operator_block), target :: mass2   !< block mass matrix
    sll_int32  :: n_dofs(3) !< number of degrees of freedom
    sll_real64 :: delta_x(3) !< cell size

    self%mass1 => mass1
    self%mass2 => mass2

    self%n_total = product(n_dofs)
    self%n_dofs = n_dofs
    self%delta_x = delta_x

    self%n_rows = 3*self%n_total
    self%n_cols = 3*self%n_total

    self%n_global_rows = self%n_rows
    self%n_global_cols = self%n_cols

  end subroutine create_linear_operator_curl_3d

  subroutine free_curl_3d( self )
    class(sll_t_linear_operator_curl_3d), intent( inout ) :: self !< Linear operator object

    self%mass1=> null()
    self%mass2=> null()

  end subroutine free_curl_3d


  subroutine dot_curl_3d ( self, x, y )
    class(sll_t_linear_operator_curl_3d), intent( in ) :: self !< Linear operator object
    sll_real64, intent( in    ) :: x(:) !< Inputvariable
    sll_real64, intent(   out ) :: y(:) !< Outputvariable
    !local variables
    sll_real64 :: scratch1(3*self%n_total), scratch2(3*self%n_total)

    ! Compute C x
    call sll_s_multiply_c( self%n_dofs, self%delta_x, x, scratch1 )

    ! Compute M2 C x
    call self%mass2%dot( scratch1, scratch2  )

    ! Compute C^T M2 C x
    call sll_s_multiply_ct( self%n_dofs, self%delta_x, scratch2, y )
!!$
    ! Compute M1 x
    call self%mass1%dot( x, scratch1 )

    ! Compute G^T M1 x
    call sll_s_multiply_gt( self%n_dofs, self%delta_x, scratch1, scratch2(1:self%n_total) )

     ! Compute G G^T M1 x
    call sll_s_multiply_g( self%n_dofs, self%delta_x, scratch2(1:self%n_total), scratch1 )

    ! Compute M1 G G^T M_1 x
    call self%mass1%dot( scratch1, scratch2 )

    ! Sum up the two parts
    y = y + scratch2
    
  end subroutine dot_curl_3d


  subroutine print_info_curl_3d( self )
    class(sll_t_linear_operator_curl_3d), intent(in) :: self !< Linear operator object

  end subroutine print_info_curl_3d


end module sll_m_linear_operator_curl_3d
