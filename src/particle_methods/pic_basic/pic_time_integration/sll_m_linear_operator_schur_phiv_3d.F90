module sll_m_linear_operator_schur_phiv_3d
#include "sll_working_precision.h"
  use sll_m_linear_operator_abstract

  use sll_m_maxwell_3d_base, only: &
       sll_c_maxwell_3d_base

  use sll_m_linear_operator_block, only : &
       sll_t_linear_operator_block

  implicit none

  public :: sll_t_linear_operator_schur_phiv_3d
  
  private
  type, extends(sll_t_linear_operator_abstract) :: sll_t_linear_operator_schur_phiv_3d
  
     class(sll_c_maxwell_3d_base), pointer :: maxwell_solver      !< Maxwell solver
     type(sll_t_linear_operator_block), pointer ::  particle_mass !< Particle mass
     sll_int32 :: n_total0, n_total1 !< total number of Dofs for 0- and 1-form
     sll_real64 :: sign = 1.0_f64 !> sign of the schur operator

   contains
     procedure :: create => create_linear_operator_schur_phiv_3d
     procedure :: free => free_schur_phiv_3d
     procedure :: dot => dot_schur_phiv_3d
     procedure :: print_info => print_info_schur_phiv_3d

  end type sll_t_linear_operator_schur_phiv_3d


contains

  subroutine create_linear_operator_schur_phiv_3d( self, maxwell_solver, particle_mass)
    class(sll_t_linear_operator_schur_phiv_3d), intent( inout ) :: self !< Schur operator
    class(sll_c_maxwell_3d_base), target :: maxwell_solver      !< Maxwell solver
    type(sll_t_linear_operator_block), target ::  particle_mass !< Particle mass


    self%particle_mass => particle_mass
    self%maxwell_solver => maxwell_solver

   self%n_total0 = maxwell_solver%n_total0
   self%n_total1 = maxwell_solver%n_total1
    
    self%n_rows = self%n_total0
    self%n_cols = self%n_total0
    
    self%n_global_rows = self%n_rows
    self%n_global_cols = self%n_cols
    
  end subroutine create_linear_operator_schur_phiv_3d

  subroutine free_schur_phiv_3d( self )
    class(sll_t_linear_operator_schur_phiv_3d), intent( inout ) :: self !< Schur operator

    self%maxwell_solver => null()
    self%particle_mass => null()

  end subroutine free_schur_phiv_3d
  
  
  subroutine dot_schur_phiv_3d ( self, x, y )
    class(sll_t_linear_operator_schur_phiv_3d), intent( in ) :: self !< Schur operator
    sll_real64, intent( in    ) :: x(:) !< Inputvariable
    sll_real64, intent(   out ) :: y(:) !< Outpoutvariable

    !local variables
    sll_real64  :: f(self%n_total1+2*self%n_total0), g(self%n_total1+2*self%n_total0), z(self%n_total0)
    !M_0 phi+ dt^2/4 q^2/m G^T M^\star G phi
    call self%maxwell_solver%multiply_g( x, f )
    call self%particle_mass%dot( f, g )
    call self%maxwell_solver%multiply_gt( g, z )
    call self%maxwell_solver%multiply_mass([0], x, y)
    y= y + self%sign*z
    
  end subroutine dot_schur_phiv_3d

  subroutine print_info_schur_phiv_3d( self )
    class(sll_t_linear_operator_schur_phiv_3d), intent(in) :: self !< Schur operator
    
  end subroutine print_info_schur_phiv_3d

  
end module sll_m_linear_operator_schur_phiv_3d
