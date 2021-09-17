module sll_m_particle_mass_3d_base
#include "sll_working_precision.h"

  use sll_m_linear_operator_abstract, only: &
       sll_t_linear_operator_abstract  

  implicit none

  public :: sll_c_particle_mass_3d_base

  private

  type, abstract, extends(sll_t_linear_operator_abstract) :: sll_c_particle_mass_3d_base
     sll_int32 :: degree(3) !< spline degree
     sll_int32 :: n_dofs(3) !< degrees of freedom
     sll_int32 :: n_total   !< product of degrees of freedom
     sll_real64, allocatable :: particle_mass(:,:) !< array containing the particle mass
     sll_real64 :: sign = 1.0_f64 !< force sign of the particle mass
     

   contains
     procedure(create_particle_mass_3d), deferred :: create 

  end type sll_c_particle_mass_3d_base

  abstract interface
     subroutine create_particle_mass_3d( self, n_cells, degree, degree2 )
       use sll_m_working_precision
       import sll_c_particle_mass_3d_base
       class(sll_c_particle_mass_3d_base), intent( inout ) :: self
       sll_int32,  intent( in ) :: n_cells(3)
       sll_int32,  intent( in ) :: degree(3)
       sll_int32, optional,  intent( in ) :: degree2(3)
       
     end subroutine create_particle_mass_3d
  end interface

end module sll_m_particle_mass_3d_base
