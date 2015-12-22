!**************************************************************
!  Copyright INRIA
!  Authors : 
!     CALVI project team
!  
!  This code SeLaLib (for Semi-Lagrangian-Library) 
!  is a parallel library for simulating the plasma turbulence 
!  in a tokamak.
!  
!  This software is governed by the CeCILL-B license 
!  under French law and abiding by the rules of distribution 
!  of free software.  You can  use, modify and redistribute 
!  the software under the terms of the CeCILL-B license as 
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info". 
!**************************************************************

module sll_m_characteristics_1d_base
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_working_precision.h"

  implicit none

  public :: &
    sll_f_process_outside_point_periodic, &
    sll_f_process_outside_point_set_to_limit, &
    sll_i_signature_process_outside_point_1d, &
    sll_c_characteristics_1d_base

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  ! For computing the characteristics in 1d
  type, abstract :: sll_c_characteristics_1d_base 
  contains
    procedure(signature_compute_characteristics_1d), deferred, pass(charac) :: &
      compute_characteristics
  
  end type sll_c_characteristics_1d_base
  
  abstract interface
    !solves eta'(t) = A(eta1(t))
    !A <=> A
    !dt <=> dt  
    !eta(dt) <=> input
    !eta(0) <=> output

    subroutine signature_compute_characteristics_1d(&
        charac, &
        A, &
        dt, &
        input, &
        output)       
      use sll_m_working_precision
      import sll_c_characteristics_1d_base       
      class(sll_c_characteristics_1d_base) :: charac
      sll_real64,dimension(:),intent(in) :: A
      sll_real64,intent(in) :: dt
      sll_real64,dimension(:),intent(in) :: input
      sll_real64,dimension(:),intent(out) :: output

    end subroutine signature_compute_characteristics_1d
  end interface

  abstract interface
    ! change the value of eta when eta<=eta_min or eta>=eta_max
    ! depending on boundary conditions
    function sll_i_signature_process_outside_point_1d( eta, eta_min, eta_max ) result(eta_out)
      use sll_m_working_precision
      sll_real64, intent(in)  :: eta
      sll_real64, intent(in) :: eta_min
      sll_real64, intent(in) :: eta_max
      sll_real64 :: eta_out      
    end function sll_i_signature_process_outside_point_1d
  end interface
  
contains

  ! periodic case
  ! called when bc_type = SLL_PERIODIC
  function sll_f_process_outside_point_periodic( eta, eta_min, eta_max ) result(eta_out)
      use sll_m_working_precision
      sll_real64, intent(in)  :: eta
      sll_real64, intent(in) :: eta_min
      sll_real64, intent(in) :: eta_max
      sll_real64 :: eta_out

      eta_out = (eta-eta_min)/(eta_max-eta_min)      
      eta_out = eta_out-real(floor(eta_out),kind=f64)
      if(eta_out==1._f64)then
        eta_out = 0._f64
      endif      
      if(.not.((eta_out>=0).and.(eta_out<1)))then
        print *,'#eta=',eta
        print *,'#eta_min=',eta_min
        print *,'#eta_max=',eta_max
        print *,'#(eta-eta_min)/(eta_max-eta_min)=',(eta-eta_min)/(eta_max-eta_min)
        print *,'#floor(-1e-19)',floor(-1e-19)
        print *,'#eta_out=',eta_out
      endif      
      SLL_ASSERT((eta_out>=0).and.(eta_out<1))
      eta_out = eta_min+eta_out*(eta_max-eta_min) 
      SLL_ASSERT((eta_out>=eta_min).and.(eta_out<eta_max))      
      
  end function sll_f_process_outside_point_periodic

  ! set to limit case
  ! called when bc_type = SLL_SET_TO_LIMIT
  
  function sll_f_process_outside_point_set_to_limit( eta, eta_min, eta_max ) result(eta_out)
      use sll_m_working_precision
      sll_real64, intent(in)  :: eta
      sll_real64, intent(in) :: eta_min
      sll_real64, intent(in) :: eta_max
      sll_real64 :: eta_out
      
      eta_out = (eta-eta_min)/(eta_max-eta_min)      
      if(eta_out>1)then
        eta_out = 1._f64
      endif
      if(eta_out<0)then
        eta_out = 0._f64
      endif
      SLL_ASSERT((eta_out>=0).and.(eta_out<=1))      
      eta_out = eta_min+eta_out*(eta_max-eta_min) 
      SLL_ASSERT((eta_out>=eta_min).and.(eta_out<=eta_max))      
      
  end function sll_f_process_outside_point_set_to_limit
  
  
  
  
end module sll_m_characteristics_1d_base
