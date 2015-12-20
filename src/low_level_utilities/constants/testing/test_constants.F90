program unit_test
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  use sll_m_constants, only: &
    sll_p_c, &
    sll_p_charge, &
    sll_p_mass, &
    sll_p_epsilon_0, &
    sll_p_g, &
    sll_p_mu_0, &
    sll_p_pi, &
    sll_p_proton_mass

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

print *,"pi=",sll_p_pi

print*, "speed of light in vacuum (def) m/s    :", sll_p_c 
print*, "permittivity of free space F/m        :", sll_p_epsilon_0 
print*, "permeability of free space N/A^2      :", sll_p_mu_0 
print*, "electron charge magnitude (49) C      :", sll_p_charge 
print*, "electron mass (54) kg                 :", sll_p_mass 
print*, "proton mass (10) kg                   :", sll_p_proton_mass 
print*, "standard grav. accel., sea level m/s^2:", sll_p_g 

end program
