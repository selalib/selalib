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

module sll_fft_utils
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_utilities.h"

  implicit none
  
contains
  !****************************************************************************
  ! How-to manipulate flags ?
  ! -------------------------
  !
  ! The differents options are stored in the integer plan%style.
  ! Each bit of the representation of this integer correspond to an option.
  !  bit          5   4   3   2   1   0 
  !  ...|---|---|---|---|---|---|---|---|
  !     |   |   | 0 | 0 | 0 | 0 | 0 | 1 |
  !  ...|---|---|---|---|---|---|---|---|
  !
  ! If the bit is equal to 1 the flag is enabled otherwise if 0 the flag is 
  ! disabled.
  ! In the example above the bit 0 is set to 1, so the flag FFT_NORMALIZE is 
  ! enabled.

  ! Return NO if the flag is disabled in the integer or YES if enabled.
  ! This function should be renamed: fft_flag_is_present or fft_is_flag_present
  function fft_is_present_flag(intege,bit) result(bool)
    sll_int32, intent(in)                   :: bit, intege
    logical                                 :: bool
    sll_int32                               :: m
   
    SLL_ASSERT( is_power_of_two( int(bit,kind=i64) ) )

    m = iand(intege,bit)
    if( m .eq. bit ) then
      bool = .true.
    else
      bool = .false.
    endif 
  end function
  ! END "How-to manipulate flags ?" section
  ! ***************************************************************************

end module sll_fft_utils
