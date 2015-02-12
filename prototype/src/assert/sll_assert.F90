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


module sll_assert
  implicit none

contains

  ! This routine just takes the informations about the assertion error and 
  ! writes them on the sreen. 
  ! This function is only meant to be used by the assert macro. No Doxygen 
  ! documentation needed.
subroutine sll_assertion(msg, file, line)
  intrinsic :: trim
  character(len=*), intent(in) :: msg
  character(len=*), intent(in) :: file
  integer, intent(in)          :: line
  character(len=8)            :: local_line
  write(local_line, '(i8)') line ! hoping that I could trim this later, but no..
  write (*,'(a, a, a, a, a)') msg, ': Assertion error triggered in file ', file, ' in line ', trim(local_line)
  stop ':  ASSERTION FAILURE'
end subroutine sll_assertion

end module sll_assert
