!---------------------------------------------------------------------------
!  Module for the Vlasov-Poisson system in Fourier coefficient
!  Laurent Navoret 2012-10
!---------------------------------------------------------------------------

module sll_m_sim_eul_vp_1d1v_cart_fe_ft

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

   use sll_m_sim_base, only: &
      sll_c_simulation_base_class

   implicit none

   public :: &
      sll_t_simulation_vp1d_fourier_fem

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   type, extends(sll_c_simulation_base_class) :: &
      sll_t_simulation_vp1d_fourier_fem
      ! Numerical parameters
      sll_real64 :: dt
   contains
      procedure, pass(sim) :: run => run_vp1d_fourier_fem
      procedure, pass(sim) :: init_from_file => VP1D_fourier_fem_init
   end type sll_t_simulation_vp1d_fourier_fem

contains

   subroutine VP1D_fourier_fem_init(sim, filename)
      class(sll_t_simulation_vp1d_fourier_fem), intent(inout) :: sim
      character(len=*), intent(in)                                 :: filename
      ! Declare here the variables to be read in through a namelist and that
      ! are to be kept inside the sim object. Look at the parallel vp4d simulation
      ! for an example.
      print *, 'This is a dummy function. Needs implementation.'
   end subroutine VP1D_fourier_fem_init

   subroutine run_vp1d_fourier_fem(sim)
      class(sll_t_simulation_vp1d_fourier_fem), intent(inout) :: sim
   end subroutine run_vp1d_fourier_fem

end module sll_m_sim_eul_vp_1d1v_cart_fe_ft
