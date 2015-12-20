program test_descriptors
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  use sll_m_descriptors, only: &
    sll_p_landau_diag, &
    sll_t_vlasovpoisson_sim

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

type(sll_t_vlasovpoisson_sim) :: testcase


call testcase%parse('  sll_p_landau_diag       ')
 print *, testcase%name()
 print *, sll_p_landau_diag%name(), sll_p_landau_diag%id

print *, "PASSED"

print *, trim('  sll_p_landau_diag       ')


end program
