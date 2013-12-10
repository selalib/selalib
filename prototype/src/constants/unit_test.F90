program unit_test
use sll_constants

implicit none

print *,"pi=",sll_pi

print*, "speed of light in vacuum (def) m/s    :", sll_c 
print*, "permittivity of free space F/m        :", sll_epsilon_0 
print*, "permeability of free space N/A^2      :", sll_mu_0 
print*, "electron charge magnitude (49) C      :", sll_e_charge 
print*, "electron mass (54) kg                 :", sll_e_mass 
print*, "proton mass (10) kg                   :", sll_proton_mass 
print*, "standard grav. accel., sea level m/s^2:", sll_g 

end program
