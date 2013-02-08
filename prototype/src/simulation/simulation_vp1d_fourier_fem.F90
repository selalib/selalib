  !---------------------------------------------------------------------------
  !  Module for the Vlasov-Poisson system in Fourier coefficient
  !  Laurent Navoret 2012-10
  !---------------------------------------------------------------------------

module simulation_vp1d_fourier_fem

#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"

#ifndef STDF95
  use sll_simulation_base
#endif
  implicit none

#ifdef STDF95
  type :: sll_simulation_vp1d_fourier_fem
#else
  type, extends(sll_simulation_base_class) :: &
    sll_simulation_vp1d_fourier_fem
#endif
    ! Numerical parameters
    sll_real64 :: dt
#ifndef STDF95
  contains
    procedure, pass(sim) :: run => run_vp1d_fourier_fem
#endif
  end type sll_simulation_vp1d_fourier_fem

contains


  ! Note that the following function has no local variables, which is silly...
  ! This just happened since the guts of the unit test were transplanted here
  ! directly, but this should be cleaned up.
  subroutine run_vp1d_fourier_fem(sim)
#ifdef STDF95
    type(sll_simulation_vp1d_fourier_fem), intent(inout) :: sim
#else
    class(sll_simulation_vp1d_fourier_fem), intent(inout) :: sim
#endif
  end subroutine run_vp1d_fourier_fem


end module simulation_vp1d_fourier_fem


