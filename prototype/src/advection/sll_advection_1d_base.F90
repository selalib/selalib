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

module sll_module_advection_1d_base
#include "sll_working_precision.h"
#include "sll_assert.h"
  implicit none
   !solves \partial_t f +A\partial_x f = 0
   ! A <=> A
   ! dt <=> dt  
   ! f(dt) <=> input
   ! f(0) <=> output
  type, abstract :: sll_advection_1d_base 
  contains
    procedure(signature_advect_1d_constant), deferred, pass(adv) :: &
      advect_1d_constant
    procedure(signature_advect_1d), deferred, pass(adv) :: &
      advect_1d
  
  end type sll_advection_1d_base

  type :: sll_advection_1d_base_ptr 
    class(sll_advection_1d_base), pointer :: ptr
  end type sll_advection_1d_base_ptr




! For A constant in 1d

 abstract interface
    subroutine signature_advect_1d_constant(&
      adv, &
      A, &
      dt, &
      input, &
      output)
      use sll_working_precision
      import sll_advection_1d_base       
      class(sll_advection_1d_base) :: adv
      sll_real64, intent(in) :: A
      sll_real64, intent(in) :: dt 
      sll_real64, dimension(:), intent(in) :: input
      sll_real64, dimension(:), intent(out) :: output
    end subroutine signature_advect_1d_constant
  end interface

! for A non constant in 1d

 abstract interface
    subroutine signature_advect_1d(&
       adv, &
       A, &
       dt, &
       input, &
       output)
      use sll_working_precision
      import sll_advection_1d_base       
      class(sll_advection_1d_base) :: adv
      sll_real64, dimension(:), intent(in) :: A
      sll_real64, intent(in) :: dt 
      sll_real64, dimension(:), intent(in) :: input
      sll_real64, dimension(:), intent(out) :: output
   end subroutine signature_advect_1d
 end interface






end module sll_module_advection_1d_base
