  !---------------------------------------------------------------------------
  !  Module for the Vlasov-Poisson system in Fourier coefficient
  !  Laurent Navoret 2012-10
  !---------------------------------------------------------------------------
! Sample computation with the following characteristics:
! - vlasov-poisson fourier in velocity space
! - 1Dx1D cartesian


program sim_eul_vp_1d1v_cart_fe_ft
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"

  use sll_m_sim_eul_vp_1d1v_cart_fe_ft, only: &
    sll_t_simulation_vp1d_fourier_fem

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    type(sll_t_simulation_vp1d_fourier_fem) :: simulation

    print *, 'Begin of vp1d_fourier_fem test'
    call simulation%run( )

    print *, 'reached end of sll_m_sim_vp1d_fourier_fem test'
    print *, 'PASSED'
    !call delete_vp4d_par_cart(simulation)

end program sim_eul_vp_1d1v_cart_fe_ft


