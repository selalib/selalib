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

module sll_m_fftw3
#ifdef FFTW_F2003
   use, intrinsic :: iso_c_binding
   include 'fftw3.f03'
#else
#ifdef __INTEL_COMPILER
   include "fftw/fftw3.f"
#else
   include "fftw3.f"
#endif
#endif
end module sll_m_fftw3
