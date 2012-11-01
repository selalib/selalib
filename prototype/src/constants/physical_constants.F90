
module physical_constants
#include "sll_working_precision.h"
  implicit none

  sll_real64, parameter :: sll_c = 2.99792458D8            ! speed of light in vacuum (def) m/s  
  sll_real64, parameter :: sll_epsilon_0 = 8.854187817D-12 ! permittivity of free space F/m      
  sll_real64, parameter :: sll_mu_0 = 12.566370614D-7      ! permeability of free space N/A^2      
  sll_real64, parameter :: sll_e_charge = 1.60217733D-19   ! electron charge magnitude (49) C      
  sll_real64, parameter :: sll_e_mass = 9.1093897D-31      ! electron mass (54) kg                  
  sll_real64, parameter :: sll_proton_mass = 1.6726231D-27 ! proton mass (10) kg                    
  sll_real64, parameter :: sll_g = 9.80665D0               ! standard grav. accel., sea level m/s^2 

end module physical_constants
