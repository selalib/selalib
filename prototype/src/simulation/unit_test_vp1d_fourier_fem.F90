  !---------------------------------------------------------------------------
  !  Module for the Vlasov-Poisson system in Fourier coefficient
  !  Laurent Navoret 2012-10
  !---------------------------------------------------------------------------
! Sample computation with the following characteristics:
! - vlasov-poisson fourier in velocity space
! - 1Dx1D cartesian


program unit_test_vp1d_fourier_fem
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use simulation_vp1d_fourier_fem
  implicit none

    type(sll_simulation_vp1d_fourier_fem) :: simulation

    print *, 'Begin of vp1d_fourier_fem test'
    call simulation%run( )

    print *, 'reached end of simulation_vp1d_fourier_fem test'
    print *, 'PASSED'
    !call delete_vp4d_par_cart(simulation)

end program unit_test_vp1d_fourier_fem


