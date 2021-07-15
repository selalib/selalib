module sll_m_linear_operator_schur_eb_cl_3d
#include "sll_working_precision.h"

  use sll_m_linear_operator_abstract

  use sll_m_spline_fem_utilities_3d_clamped, only : &
       sll_s_multiply_c_clamped, &
       sll_s_multiply_ct_clamped

  use sll_m_linear_operator_block, only : &
       sll_t_linear_operator_block

!!$  use sll_m_matrix_csr, only: &
!!$       sll_t_matrix_csr
  
  implicit none

  public :: sll_t_linear_operator_schur_eb_cl_3d
  
  private
  type, extends(sll_t_linear_operator_abstract) :: sll_t_linear_operator_schur_eb_cl_3d
  
     type(sll_t_linear_operator_block), pointer :: mass1   !< block mass matrix
     type(sll_t_linear_operator_block), pointer :: mass2   !< block mass matrix
     !type(sll_t_matrix_csr), pointer            :: mass1d(:,:)
     sll_int32 :: n_total0, n_total1  !< product of number of degrees of freedom
     sll_int32 :: n_dofs(3) !< number of degrees of freedom
     sll_int32 :: s_deg_0(3) !< spline degree 0-forms 
     sll_real64 :: sign = 1.0_f64 !< sign
     sll_real64 :: delta_x(3) !< cell size
     
   contains
     procedure :: create => create_linear_operator_schur_eb_cl_3d
     procedure :: free => free_schur_eb_cl_3d
     procedure :: dot => dot_schur_eb_cl_3d
     procedure :: print_info => print_info_schur_eb_cl_3d

  end type sll_t_linear_operator_schur_eb_cl_3d


contains

  subroutine create_linear_operator_schur_eb_cl_3d( self, mass1, mass2, n_dofs, delta_x, s_deg_0 )
    class(sll_t_linear_operator_schur_eb_cl_3d), intent( inout ) :: self !< Linear operator object
    type(sll_t_linear_operator_block), target :: mass1   !< block mass matrix
    type(sll_t_linear_operator_block), target :: mass2   !< block mass matrix
    !type(sll_t_matrix_csr), target            :: mass1d(3,3)
    sll_int32  :: n_dofs(3) !< number of degrees of freedom
    sll_real64 :: delta_x(3) !< cell size
    sll_int32  ::  s_deg_0(3) !< spline degree 0-forms 
    
    self%mass1  => mass1
    self%mass2  => mass2
    !self%mass1d => mass1d

    self%n_dofs = n_dofs
    self%delta_x = delta_x
    self%s_deg_0 = s_deg_0
    self%n_total0 = product(self%n_dofs)
    self%n_total1 = (self%n_dofs(1)-1)*self%n_dofs(2)*self%n_dofs(3)
    
    
    
    self%n_rows = self%n_total1+2*self%n_total0
    self%n_cols = self%n_total1+2*self%n_total0
    
    self%n_global_rows = self%n_rows
    self%n_global_cols = self%n_cols
    
  end subroutine create_linear_operator_schur_eb_cl_3d

  subroutine free_schur_eb_cl_3d( self )
    class(sll_t_linear_operator_schur_eb_cl_3d), intent( inout ) :: self !< Linear operator object

    self%mass1=> null()
    self%mass2=> null()
    !self%mass1d=> null()
    
  end subroutine free_schur_eb_cl_3d
  
  
  subroutine dot_schur_eb_cl_3d ( self, x, y )
    class(sll_t_linear_operator_schur_eb_cl_3d), intent( in ) :: self !< Linear operator object
    sll_real64, intent( in    ) :: x(:) !< Inputvariable
    sll_real64, intent(   out ) :: y(:) !< Outputvariable
    !local variables
    sll_real64 :: scratch0(1:self%n_total0), scratch1(self%n_total1+2*self%n_total0), scratch2(self%n_total0+2*self%n_total1), scratch3(self%n_total0+2*self%n_total1)
    
    ! Compute C x
    call sll_s_multiply_c_clamped( self%n_dofs, self%delta_x, self%s_deg_0, x, scratch2 )
    
    ! Compute M2 C x
    call self%mass2%dot( scratch2, scratch3 )
    
    ! Compute C^T M2 C x
    call sll_s_multiply_ct_clamped( self%n_dofs, self%delta_x, self%s_deg_0, scratch3, scratch1 )
!!$
!!$    ! Compute C^T M2 C x + M_b C x
!!$    call multiply_mass_2dkron( self, [2,2,1], scratch2(1+self%n_total0+self%n_total1:self%n_total0+2*self%n_total1), scratch0 )
!!$    scratch1(1+self%n_total1:self%n_total1+self%n_total0)=scratch1(1+self%n_total1:self%n_total1+self%n_total0) - scratch0 
!!$    call multiply_mass_2dkron( self, [2,1,2], scratch2(1+self%n_total0:self%n_total0+self%n_total1), scratch0 )
!!$    scratch1(1+self%n_total1+self%n_total0:self%n_total1+2*self%n_total0)=scratch1(1+self%n_total1+self%n_total0:self%n_total1+2*self%n_total0) + scratch0
   
    
    ! Compute M1 x
    call self%mass1%dot( x, y )
    
    ! Sum up the two parts
    y = y + self%sign * scratch1
    
  end subroutine dot_schur_eb_cl_3d

  
  subroutine print_info_schur_eb_cl_3d( self )
    class(sll_t_linear_operator_schur_eb_cl_3d), intent(in) :: self !< Linear operator object
    
  end subroutine print_info_schur_eb_cl_3d

!!$
!!$  !> Multiply by the mass matrix 
!!$  subroutine multiply_mass_2dkron(  self, deg, coefs_in, coefs_out )
!!$    class(sll_t_linear_operator_schur_eb_cl) :: self
!!$    sll_int32,  intent( in   )  :: deg(:) !< \a deg(i) specifies the degree of the 1d mass matrix in dimension \a i (Note: 1 for 0-form, 2 for 1-form, 3 for 0-1-form mix)
!!$    sll_real64, intent( in    )  :: coefs_in(:)
!!$    sll_real64, intent(   out )  :: coefs_out(:)  
!!$    ! Local variables
!!$    sll_int32 :: j,k
!!$    sll_real64 :: work3d(self%n_dofs(1), self%n_dofs(2), self%n_dofs(3))
!!$    sll_real64 :: work_d2_in(self%n_dofs(2)), work_d2_out(self%n_dofs(2))
!!$    sll_real64 :: work_d3_in(self%n_dofs(3)), work_d3_out(self%n_dofs(3))
!!$    
!!$    coefs_out=0._f64
!!$    if (deg(1)==2)then
!!$       work3d=0._f64
!!$       do k=1,self%n_dofs(3)
!!$          do j=1,self%n_dofs(2)
!!$             work3d(1,j,k) = -coefs_in(1+(j-1)*(self%n_dofs(1)-1)+(k-1)*(self%n_dofs(1)-1)*self%n_dofs(2))
!!$             work3d(self%n_dofs(1),j,k) = coefs_in(self%n_dofs(1)-1+(j-1)*(self%n_dofs(1)-1)+(k-1)*(self%n_dofs(1)-1)*self%n_dofs(2))
!!$          end do
!!$       end do
!!$
!!$       do k=1,self%n_dofs(3)
!!$          work_d2_in = work3d(1,:,k)
!!$          call self%mass1d(deg(2),2)%dot( work_d2_in, work_d2_out )
!!$          work3d(1,:,k) = work_d2_out
!!$
!!$          work_d2_in = work3d(self%n_dofs(1),:,k)
!!$          call self%mass1d(deg(2),2)%dot( work_d2_in, work_d2_out )
!!$          work3d(self%n_dofs(1),:,k) = work_d2_out
!!$       end do
!!$
!!$       do j=1,self%n_dofs(2)
!!$          work_d3_in = work3d(1,j,:)
!!$          call self%mass1d(deg(3),3)%dot( work_d3_in, work_d3_out )
!!$          do k=1,self%n_dofs(3)
!!$             coefs_out(1+(j-1)*self%n_dofs(1)+(k-1)*self%n_dofs(1)*self%n_dofs(2)) = work_d3_out(k)
!!$          end do
!!$          work_d3_in = work3d(self%n_dofs(1),j,:)
!!$          call self%mass1d(deg(3),3)%dot( work_d3_in, work_d3_out )
!!$          do k=1,self%n_dofs(3)
!!$             coefs_out(self%n_dofs(1)+(j-1)*self%n_dofs(1)+(k-1)*self%n_dofs(1)*self%n_dofs(2)) = work_d3_out(k)
!!$          end do
!!$       end do
!!$    else
!!$       print*, 'multiply_mass_2dkron not implemented for this degree'
!!$       stop
!!$    end if
!!$   
!!$  end subroutine multiply_mass_2dkron

  

  
end module sll_m_linear_operator_schur_eb_cl_3d
