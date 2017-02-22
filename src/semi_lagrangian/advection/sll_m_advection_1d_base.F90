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

!> @ingroup advection
!> @brief
!> Abstract class for advection
!> @details
!> Solves \f$ \partial_t f + A \partial_x f = 0 \f$
module sll_m_advection_1d_base
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  implicit none

  public :: &
    sll_c_advection_1d_base, &
    sll_t_advection_1d_base_ptr

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

type, abstract :: sll_c_advection_1d_base 

contains

  procedure(signature_advect_1d_constant), deferred, pass(adv) :: advect_1d_constant
  procedure(signature_advect_1d),          deferred, pass(adv) :: advect_1d
  procedure(signature_advect_1d_delete),   deferred, pass(adv) :: delete

end type sll_c_advection_1d_base

type :: sll_t_advection_1d_base_ptr 
  class(sll_c_advection_1d_base), pointer :: ptr
end type sll_t_advection_1d_base_ptr

#ifndef DOXYGEN_SHOULD_SKIP_THIS
abstract interface

  subroutine signature_advect_1d_constant( adv,   &
                                           a,     &
                                           dt,    &
                                           input, &
                                           output)

    use sll_m_working_precision
    import sll_c_advection_1d_base       

    class(sll_c_advection_1d_base)          :: adv
    sll_real64,               intent(in)  :: a
    sll_real64,               intent(in)  :: dt 
    sll_real64, dimension(:), intent(in)  :: input
    sll_real64, dimension(:), intent(out) :: output

  end subroutine signature_advect_1d_constant

end interface

abstract interface

  subroutine signature_advect_1d( adv,       &
                                    A,       &
                                    dt,      &
                                    input,   &
                                    output)

    use sll_m_working_precision
    import sll_c_advection_1d_base       

    class(sll_c_advection_1d_base)          :: adv
    sll_real64, dimension(:), intent(in)  :: a
    sll_real64,               intent(in)  :: dt 
    sll_real64, dimension(:), intent(in)  :: input
    sll_real64, dimension(:), intent(out) :: output

  end subroutine signature_advect_1d

end interface

abstract interface

  subroutine signature_advect_1d_delete( adv )

    import sll_c_advection_1d_base       
    class(sll_c_advection_1d_base), intent(inout) :: adv

  end subroutine signature_advect_1d_delete

end interface


#endif /* DOXYGEN_SHOULD_SKIP_THIS */

end module sll_m_advection_1d_base
