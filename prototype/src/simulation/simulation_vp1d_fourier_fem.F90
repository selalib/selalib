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
    procedure, pass(sim) :: init_from_file => VP1D_fourier_fem_init
#endif
  end type sll_simulation_vp1d_fourier_fem

contains

  subroutine VP1D_fourier_fem_init(sim, filename)
#ifdef STDF95
    type(sll_simulation_vp1d_fourier_fem), intent(inout)  :: sim
#else
    class(sll_simulation_vp1d_fourier_fem), intent(inout) :: sim
#endif
    character(len=*), intent(in)                                 :: filename
    ! Declare here the variables to be read in through a namelist and that
    ! are to be kept inside the sim object. Look at the parallel vp4d simulation
    ! for an example.
    print *, 'This is a dummy function. Needs implementation.'
  end subroutine VP1D_fourier_fem_init


  subroutine run_vp1d_fourier_fem(sim)
#ifdef STDF95
    type(sll_simulation_vp1d_fourier_fem), intent(inout) :: sim
#else
    class(sll_simulation_vp1d_fourier_fem), intent(inout) :: sim
#endif
  end subroutine run_vp1d_fourier_fem


end module simulation_vp1d_fourier_fem


