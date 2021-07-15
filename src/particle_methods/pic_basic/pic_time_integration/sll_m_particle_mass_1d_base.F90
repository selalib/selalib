module sll_m_particle_mass_1d_base
#include "sll_working_precision.h"

  use sll_m_linear_operator_abstract, only: &
       sll_t_linear_operator_abstract  

  implicit none

  public :: sll_c_particle_mass_1d_base

  private

  type, abstract, extends(sll_t_linear_operator_abstract) :: sll_c_particle_mass_1d_base

     sll_int32 :: degree !< spline degree
     sll_int32 :: n_dofs !< degrees of freedom
     sll_real64, allocatable :: particle_mass(:,:) !< array containing the particle mass
     sll_real64 :: sign = 1.0_f64 !< force sign of the particle mass



   contains
     procedure(create_particle_mass_1d), deferred :: create 

  end type sll_c_particle_mass_1d_base

  abstract interface
     subroutine create_particle_mass_1d( self, degree, n_dofs )
       use sll_m_working_precision
       import sll_c_particle_mass_1d_base
       class(sll_c_particle_mass_1d_base), intent( inout ) :: self
       sll_int32,  intent( in ) :: degree
       sll_int32,  intent( in ) :: n_dofs

     end subroutine create_particle_mass_1d
  end interface

end module sll_m_particle_mass_1d_base
