!> @ingroup distribution function
!> @author Benedikt Perse
!> @brief functions for initial profile of the particle distribution function
module sll_m_profile_functions
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"
#include "sll_errors.h"

  use sll_m_constants, only : &
       sll_p_pi, sll_p_twopi

  implicit none

  public :: &
       sll_t_profile_functions

  private
  
  type :: sll_t_profile_functions
     sll_real64 :: a = 0._f64
     sll_real64 :: rbar = 0.5_f64
     sll_real64 :: rhobar = 1._f64
     sll_real64 :: kr = 0._f64
     sll_real64 :: wr = 1._f64
     sll_real64 :: tibar = 1._f64
     sll_real64 :: kti = 0._f64
     sll_real64 :: wti = 1._f64
     sll_real64 :: tebar = 1._f64
     sll_real64 :: kte = 0._f64
     sll_real64 :: wte = 1._f64

   contains
     procedure :: init => init_profile_functions
     procedure :: rho_0
     procedure :: drho_0
     procedure :: T_i
     procedure :: dT_i
     procedure :: T_e
     procedure :: radial_distrib
     !procedure :: dradial_distrib
  end type sll_t_profile_functions

contains


  subroutine init_profile_functions(self, input_file)
    class(sll_t_profile_functions),   intent( inout ) :: self
    sll_int32, intent(in) :: input_file
    !local variables
    sll_int32 :: io_stat
    sll_real64 :: a, rbar, rhobar, kr, omegar, tibar, kti, omegati, tebar, kte, omegate

    namelist /profile_parameters/ a, rbar, rhobar, kr, omegar, tibar, kti, omegati, tebar, kte, omegate

    read(input_file, profile_parameters,IOStat=io_stat)
!!$    open(newunit = input_file, file=trim(filename), status='old',IOStat=io_stat)
    if (io_stat /= 0) then
       self%a = 0._f64
       self%rbar = 0.5_f64
       self%rhobar = 1._f64
       self%kr = 0._f64
       self%wr = 1._f64
       self%tibar = 1._f64
       self%kti = 0._f64
       self%wti = 1._f64
       self%tebar = 1._f64
       self%kte = 0._f64
       self%wte = 1._f64
       
    else
       self%a = a
       self%rbar = rbar
       self%rhobar = rhobar
       self%kr = kr
       self%wr = omegar
       self%tibar = tibar
       self%kti = kti
       self%wti = omegati
       self%tebar = tebar
       self%kte = kte
       self%wte = omegate
    end if

  end subroutine init_profile_functions

  function rho_0(self, r)
    class(sll_t_profile_functions),   intent( inout ) :: self
    sll_real64 :: rho_0
    sll_real64, intent(in) :: r

    rho_0 = self%rhobar * exp( -self%kr*self%wr*tanh(self%a*(r-self%rbar)/self%wr) )
  end function rho_0


  function drho_0(self, r)
    class(sll_t_profile_functions),   intent( inout ) :: self
    sll_real64 :: drho_0
    sll_real64, intent(in) :: r
    
    drho_0 = -self%kr/(cosh(self%a*(r-self%rbar)/self%wr)**2) *self%rho_0(r)
  end function drho_0
  
  function T_i(self, r)
    class(sll_t_profile_functions),   intent( inout ) :: self
    sll_real64 :: T_i
    sll_real64, intent(in) :: r

    T_i = self%tibar*(1._f64-self%kti*self%wti*tanh(self%a*(r-self%rbar)/self%wti) )
    !T_i = self%tibar * exp( -self%kti*self%wti*tanh(self%a*(r-self%rbar)/self%wti) )
    if( T_i< 0._f64) then
       SLL_ERROR('T_i','error in temperature in T_i')
    end if
  end function T_i

  function dT_i(self, r)
    class(sll_t_profile_functions),   intent( inout ) :: self
    sll_real64 :: dT_i
    sll_real64, intent(in) :: r

    dT_i = -self%tibar*self%kti/cosh(self%a*(r-self%rbar)/self%wti)**2
    !dT_i = -self%kti/(cosh(self%a*(r-self%rbar)/self%wti)**2)*self%T_i(r)
  
  end function dT_i

  function T_e(self, r)
    class(sll_t_profile_functions),   intent( inout ) :: self
    sll_real64 :: T_e
    sll_real64, intent(in) :: r

    T_e = self%tebar * exp( -self%kte*self%wte*tanh(self%a*(r-self%rbar)/self%wte) )
    if( T_e< 0._f64) then
       SLL_ERROR('T_e','error in temperature in T_e')
    end if

  end function T_e

  function radial_distrib(self, r)
    class(sll_t_profile_functions),   intent( inout ) :: self
    sll_real64 :: radial_distrib
    sll_real64, intent(in) :: r

    radial_distrib = sin(sll_p_pi*r)
    !radial_distrib = exp(-self%wti*(self%a*(r-self%rbar))**2/(4._f64*self%wr) )
  end function radial_distrib

!!$  function dradial_distrib(self, r)
!!$    class(sll_t_profile_functions),   intent( inout ) :: self
!!$    sll_real64 :: dradial_distrib
!!$    sll_real64, intent(in) :: r
!!$
!!$    dradial_distrib = 0._f64!-self%a*2._f64*(self%a*(r-rbar))/(self%wr/self%wti)*self%radial_distrib(r)
!!$  end function dradial_distrib

end module sll_m_profile_functions
