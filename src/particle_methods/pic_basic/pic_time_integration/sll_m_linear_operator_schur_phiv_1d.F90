module sll_m_linear_operator_schur_phiv_1d
#include "sll_working_precision.h"

  use sll_m_linear_operator_abstract, only: &
       sll_t_linear_operator_abstract

  use sll_m_maxwell_1d_base, only: &
       sll_c_maxwell_1d_base

  use sll_m_particle_mass_1d_base, only : &
       sll_c_particle_mass_1d_base

  
  implicit none

  public :: sll_t_linear_operator_schur_phiv_1d
  
  private
  type, extends(sll_t_linear_operator_abstract) :: sll_t_linear_operator_schur_phiv_1d
  
     class(sll_c_maxwell_1d_base), pointer :: maxwell_solver      !< Maxwell solver
     class(sll_c_particle_mass_1d_base), pointer ::  particle_mass !< Particle mass
     sll_int32  :: n_dofs0, n_dofs1 !< number of Dofs for 0- and 1-form
     sll_real64 :: sign = 1.0_f64 !< sign of the schur operator
     sll_int32  :: degree !< spline degree

   contains
     procedure :: create => create_linear_operator_schur_phiv_1d
     procedure :: free => free_schur_phiv_1d
     procedure :: dot => dot_schur_phiv_1d
     procedure :: print_info => print_info_schur_phiv_1d

  end type sll_t_linear_operator_schur_phiv_1d


contains

  subroutine create_linear_operator_schur_phiv_1d( self, maxwell_solver, particle_mass )
    class(sll_t_linear_operator_schur_phiv_1d), intent( inout ) :: self !< Schur operator
    class(sll_c_maxwell_1d_base), target :: maxwell_solver      !< Maxwell solver
    class( sll_c_particle_mass_1d_base), target ::  particle_mass !< Particle mass
         
    self%particle_mass => particle_mass
    self%maxwell_solver => maxwell_solver

    seLf%degree = self%maxwell_solver%s_deg_1
    self%n_dofs0 = self%maxwell_solver%n_dofs0
    self%n_dofs1 = self%maxwell_solver%n_dofs1
    
    self%n_rows = self%n_dofs0
    self%n_cols = self%n_dofs0
    
    self%n_global_rows = self%n_rows
    self%n_global_cols = self%n_cols
    
  end subroutine create_linear_operator_schur_phiv_1d

  subroutine free_schur_phiv_1d( self )
    class(sll_t_linear_operator_schur_phiv_1d), intent( inout ) :: self !< Schur operator

    self%maxwell_solver => null()
    self%particle_mass => null()

  end subroutine free_schur_phiv_1d
  
  
  subroutine dot_schur_phiv_1d ( self, x, y )
    class(sll_t_linear_operator_schur_phiv_1d), intent( in ) :: self !< Schur operator
    sll_real64, intent( in    ) :: x(:) !< Inputvariable
    sll_real64, intent(   out ) :: y(:) !< Outpoutvariable
    !local variables
    sll_real64  :: g(self%n_dofs1), h(self%n_dofs1), z(self%n_dofs0)
    
    !M_0 phi+ dt^2/4 q^2/m G^T M^\star G phi
    call self%maxwell_solver%multiply_g( x, g )
    call self%particle_mass%dot( g, h )
    call self%maxwell_solver%multiply_gt( h, z )
    call self%maxwell_solver%multiply_mass( x, y, self%maxwell_solver%s_deg_0 )
    y = y + self%sign*z

  end subroutine dot_schur_phiv_1d
  
  subroutine print_info_schur_phiv_1d( self )
    class(sll_t_linear_operator_schur_phiv_1d), intent(in) :: self !< Schur operator
    
  end subroutine print_info_schur_phiv_1d
 
end module sll_m_linear_operator_schur_phiv_1d
