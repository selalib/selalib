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

module sll_strings
implicit none

interface operator (.conc.)
   module procedure sll_concatenate_strings
end interface

contains

function sll_concatenate_strings(s1, s2)
  character(len=*), intent(in) :: s1
  character(len=*), intent(in) :: s2
  character(len=(len(s1)+len(s2))) :: sll_concatenate_strings
  sll_concatenate_strings = s1//s2
end function sll_concatenate_strings



end module sll_strings

