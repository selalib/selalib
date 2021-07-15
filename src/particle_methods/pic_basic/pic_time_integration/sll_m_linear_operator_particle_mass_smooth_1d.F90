module sll_m_linear_operator_particle_mass_smooth_1d
#include "sll_working_precision.h"

  use sll_m_particle_mass_1d_base, only: &
       sll_c_particle_mass_1d_base


  implicit none

  public :: sll_t_linear_operator_particle_mass_smooth_1d
  
  private

  type, extends(sll_c_particle_mass_1d_base) :: sll_t_linear_operator_particle_mass_smooth_1d
    

   contains
     procedure :: create => create_linear_operator_particle_mass_smooth_1d
     procedure :: free => free_particle_mass_smooth_1d
     procedure :: dot => dot_particle_mass_smooth_1d
     procedure :: print_info => print_info_particle_mass_smooth_1d

  end type sll_t_linear_operator_particle_mass_smooth_1d


contains

  subroutine create_linear_operator_particle_mass_smooth_1d( self, degree, n_dofs )
    class(sll_t_linear_operator_particle_mass_smooth_1d), intent( inout ) :: self !< Particle mass
    sll_int32,  intent( in ) :: degree !< spline degree
    sll_int32,  intent( in ) :: n_dofs !< number of degrees of freedom
       
    self%degree = degree
    self%n_dofs = n_dofs

    allocate( self%particle_mass( 2*self%degree+5, self%n_dofs) )
    self%particle_mass = 0._f64
    self%n_rows = self%n_dofs
    self%n_cols = self%n_dofs
    
    self%n_global_rows = self%n_rows
    self%n_global_cols = self%n_cols
    
  end subroutine create_linear_operator_particle_mass_smooth_1d

  subroutine free_particle_mass_smooth_1d( self )
    class(sll_t_linear_operator_particle_mass_smooth_1d), intent( inout ) :: self !< Particle mass

    deallocate( self%particle_mass )

  end subroutine free_particle_mass_smooth_1d
  
  
  subroutine dot_particle_mass_smooth_1d ( self, x, y )
    class(sll_t_linear_operator_particle_mass_smooth_1d), intent( in ) :: self !< Particle mass
    sll_real64, intent( in    ) :: x(:) !< Inputvariable
    sll_real64, intent(   out ) :: y(:) !< Outpoutvariable
    !local variables
    sll_int32 :: ind
    sll_int32 :: i, j
    
    do i = 1, self%n_dofs
      y(i) = 0.0_f64
      ind = 1
      do j = i-self%degree-2, i+self%degree+2
        y(i) = y(i) + self%sign * self%particle_mass( ind, i) * x(1 + modulo(j-1,self%n_dofs))
        ind = ind+1
      end do
    end do
                    
  
    
  end subroutine dot_particle_mass_smooth_1d


  subroutine print_info_particle_mass_smooth_1d( self )
    class(sll_t_linear_operator_particle_mass_smooth_1d), intent(in) :: self !< Particle mass
  end subroutine print_info_particle_mass_smooth_1d
  
end module sll_m_linear_operator_particle_mass_smooth_1d
