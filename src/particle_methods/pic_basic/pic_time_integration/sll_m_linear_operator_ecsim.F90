module sll_m_linear_operator_ecsim
#include "sll_working_precision.h"
  
  use sll_m_linear_operator_abstract
  
  use sll_m_spline_fem_utilities, only : &
       sll_s_spline_fem_mass_line

  implicit none
  
  public :: sll_t_linear_operator_ecsim
  
  private

  type, extends(sll_t_linear_operator_abstract) :: sll_t_linear_operator_ecsim
     sll_int32 :: degree !< Spline degree
     sll_real64 :: dt !< Timestep
     sll_real64 :: dx !< Grid spacing
     sll_real64, allocatable :: mass_line_0(:) !< Entries of mass matrix M0
     sll_real64, allocatable :: mass_line_1(:) !< Entries of mass matrix M1
     sll_real64, pointer :: m1(:,:), m2(:,:), m4(:,:) !< Blockmatrices of particle-mass matrix 
     
     
   contains
     procedure :: create => create_linear_operator_ecsim
     procedure :: free => free_ecsim
     procedure :: dot  => dot_mono_r2r_ecsim
     procedure :: print_info => print_info_ecsim
     
  end type sll_t_linear_operator_ecsim
  
contains
  
  subroutine create_linear_operator_ecsim( self, n_dof,degree, mass_line_0,mass_line_1, m1,m2, m4)
    class(sll_t_linear_operator_ecsim), intent( inout ) :: self !< Linear operator object
    sll_int32 :: n_dof !< Number of degrees of freedom
    sll_int32 :: degree !< Splinedegree
    sll_real64 :: mass_line_0(:) !< Entries of mass matrix M0
    sll_real64 :: mass_line_1(:) !< Entries of mass matrix M1
    sll_real64,target :: m1(:,:) !< Blockmatrix of particle-mass matrix
    sll_real64,target :: m2(:,:) !< Blockmatrix of particle-mass matrix
    sll_real64,target :: m4(:,:) !< Blockmatrix of particle-mass matrix
    
    
    self%n_dof=n_dof
    self%degree=degree
    self%dt=0.0_f64
    self%m1=>m1
    self%m2=>m2
    self%m4=>m4
    allocate(self%mass_line_0(degree+1))
    allocate(self%mass_line_1(degree))
    self%mass_line_0=mass_line_0
    self%mass_line_1=mass_line_1

    
    self%n_rows = 3*n_dof
    self%n_cols = 3*n_dof
    
    self%n_global_rows = self%n_rows
    self%n_global_cols = self%n_cols
  end subroutine create_linear_operator_ecsim
  

  subroutine free_ecsim( self )
    class(sll_t_linear_operator_ecsim), intent( inout ) :: self !< Linear operator object
    deallocate(self%mass_line_0)
    deallocate(self%mass_line_1)
    self%m1=>null()
    self%m2=>null()
    self%m4=>null()
    
  end subroutine free_ecsim

  subroutine dot_mono_r2r_ecsim( self, x, y )
    class(sll_t_linear_operator_ecsim), intent( in ) :: self !< Linear operator object
    sll_real64, intent( in    ) :: x(:) !< Inputvariable
    sll_real64, intent(   out ) :: y(:) !< Outputvariable
    ! local variables
    sll_int32 :: row, column
    sll_real64 :: scratch(self%n_dof+1) 
    scratch=0.0_f64
    y=0.0_f64
    do row=1, self%n_dof
       ! Monoplication of the inputvector with entries of the mass and particle-mass matrices
       y(row) = (self%mass_line_1(1)+0.5_f64*self%dt*self%m1(1,row))*x(row)
       y(self%n_dof+row)= (self%mass_line_0(1)+0.5_f64*self%dt*self%m4(1,row))*x(self%n_dof+row)

       scratch(row)=self%mass_line_1(1)*x(2*self%n_dof+row)
       do column = 2, self%degree
          y(row) = y(row) +&
               self%mass_line_1(column) * &
               (x(modulo(row+column-2,self%n_dof)+1) +&
               x(modulo(row-column,self%n_dof)+1))+&
               0.5_f64*self%dt*self%m1(column,row) * &
               x(modulo(row+column-2,self%n_dof)+1)+& 
               0.5_f64*self%dt*self%m1(column,modulo(row-column,self%n_dof)+1) * &
               x(modulo(row-column,self%n_dof)+1) 

          scratch(row)=scratch(row)+&
               self%mass_line_1(column) * &  
               (x(2*self%n_dof+ modulo(row+column-2,self%n_dof)+1) +&
               x(2*self%n_dof+ modulo(row-column,self%n_dof)+1))
       end do
       do column = 2, self%degree+1
          y(self%n_dof+row) = y(self%n_dof+row) +&
               self%mass_line_0(column) * &
               (x(self%n_dof+modulo(row+column-2,self%n_dof)+1) +&
               x(self%n_dof+ modulo(row-column,self%n_dof)+1))+&
            0.5_f64*self%dt*self%m4(column,row)*&
            x(self%n_dof+modulo(row+column-2,self%n_dof)+1)+&
            0.5_f64*self%dt*self%m4(column,modulo(row-column,self%n_dof)+1)*&
            x(self%n_dof+modulo(row-column,self%n_dof)+1)
       end do

       do column = 1, 2*self%degree
          y(row) = y(row) - &             
               0.5_f64*self%dt*self%m2(7-column,modulo(column-4+row-1,self%n_dof)+1)*&
               x(self%n_dof+modulo(row+column-self%degree-2,self%n_dof)+1) 

          y(row+self%n_dof) = y(row+self%n_dof) + &
               0.5_f64*self%dt*self%m2(column,row) * &
               x(modulo(row+column-self%degree-1,self%n_dof)+1)
       end do
    end do
    ! Update of the outputvector with derivatives of the inputvector 
    scratch(self%n_dof+1)=scratch(1)
    do row=1, self%n_dof
       y(row+self%n_dof)=  y(row+self%n_dof)-&
            0.5_f64*self%dt/self%dx*(scratch(row)-scratch(row+1))

       y(2*self%n_dof+row)=x(2*self%n_dof+row)+ &
            0.5_f64*self%dt/self%dx*(x(self%n_dof+row)-x(self%n_dof+modulo(row-2,self%n_dof)+1))

    end do

  end subroutine dot_mono_r2r_ecsim

  
  subroutine print_info_ecsim( self )
    class(sll_t_linear_operator_ecsim), intent(in) :: self !< Linear operator object
  end subroutine print_info_ecsim

end module sll_m_linear_operator_ecsim
