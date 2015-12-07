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
! MODULE: sll_m_memory
!
!> @ingroup memory
!> @author
!> Edwin Chacon-Golcher INRIA
!> @brief
!> Implements the error testing function and other related functionalities.
module sll_m_memory
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  implicit none

  public :: &
    test_error_code

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Useful for 'alternative implementation'. Ignore for now.
#if 0
 interface sll_allocate
     module procedure sll_allocate_int_1d, sll_allocate_real_1d, &
          sll_allocate_complex_1d
  end interface
#endif

contains 

   !---------------------------------------------------------------------------  
   !> @brief
   !> Test if there is an error and write it on the screen.
   !> @details
   !! We assume that an error code is set to 0
   !! if there is no problem at the function output.
   !! Thus, we check the value of the error code and if it's different to 0,
   !! the routine writes an error message on the screen. 
   !> @param[in] err_code the error code to be tested
   !> @param[in] descriptor short description of the error
   !> @param[in] file_name file name in which the error occurred
   !> @param[in] line_number line number of the call
   !---------------------------------------------------------------------------

  subroutine test_error_code(err_code, descriptor, file_name, line_number)
    integer                      :: err_code
    character(len=*), intent(in) :: descriptor
    character(len=*), intent(in) :: file_name
    integer, intent(in)          :: line_number
    if (err_code .ne. 0) then
       write(*, '(a, a, i8)') descriptor, ' Triggered in FILE '//file_name// &
            ', in LINE: ', line_number
       stop 'ERROR: test_error_code(): exiting program'
    end if
  end subroutine test_error_code


end module sll_m_memory
