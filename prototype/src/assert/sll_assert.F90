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

!------------------------------------------------------------------------------
! SELALIB
!------------------------------------------------------------------------------
!
! MODULE: sll_assert
!
!> @author
!> Module Author Name and Affiliation
!
! DESCRIPTION: 
!> This is a very small but very useful capability that permits
!! the developer to devise many sorts of defensive programming techniques.
!!
!! The simple idea is that of declaring a condition that is expected to be 
!! true, thus triggering a descriptive error if the condition is violated.
!
! REVISION HISTORY:
! DD Mmm YYYY - Initial Version
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!------------------------------------------------------------------------------

module sll_assert
  implicit none

contains

   !---------------------------------------------------------------------------  
   !> @author 
   !> Routine Author Name.
   !
   ! DESCRIPTION: 
   !> Output pattern in the case of assertion error.
   !!
   !! This routine just takes the informations about the assertion error and writes them on the sreen. 
   !
   ! REVISION HISTORY:
   ! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
   !
   !> @param[in] msg short description of the error
   !> @param[in] file file name in which the error occurred
   !> @param[in] line line number of the call
   !---------------------------------------------------------------------------

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
