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

module sll_utilities
#include "sll_working_precision.h"
  implicit none
!  intrinsic :: selected_int_kind ! this line gives an error, why?

  ! Tentative implementation of a standard-compliant way to get the
  ! memory footprint of a variable. This is our yardstick...
  !
  ! selected_int_kind(r) returns the default integer scalar that is the kind 
  ! type parameter value for an integer data type able to represent all 
  ! integer values n in the range -10^(-r) < n < 10^r, where r is a scalar 
  ! integer. If more than one is available, a kind with least decimal exponent 
  ! range is chosen (and least kind value if several have least decimal 
  ! exponent range). If no corresponding kind is availalble, the result is -1. 
  ! (Metcalf & Reid. F90/95 2nd ed. p. 176).
  !
  ! We are, maybe dangerously, relying on the common practice of many compilers
  ! of using the kind values to indicate the number of bytes of storage
  ! occupied by a value. But this is not mandated by the standard. For 
  ! instance a compiler that only has available 4-byte integers may still
  ! support kind values of 1, 2 and 4 to 'ease portability' from other 
  ! platforms. The size of k1 will become our basic 'yardstick' to measure
  ! the size of a memory footprint of a variable. When we ask for the size
  ! of 'var', the answer will be given in terms of how many 'yardsticks'
  ! are needed to represent 'var'. The implications of this assumption
  ! need to be checked further.

  integer, parameter :: byte_size = selected_int_kind(0)

contains

  function is_power_of_two( n )
    intrinsic             :: not, iand
    sll_int64, intent(in) :: n
    logical               :: is_power_of_two
    if( (n>0) .and. (0 .eq. (iand(n,(n-1)))) ) then
       is_power_of_two = .true.
    else
       is_power_of_two = .false.
    end if
  end function is_power_of_two

  function is_even( n )
    intrinsic             :: modulo    
    sll_int32, intent(in) :: n
    logical               :: is_even
    if( modulo(n,2) .eq. 0 ) then
       is_even = .true.
    else
       is_even = .false.
    end if
  end function is_even

  subroutine int2string( istep, cstep )
    integer, intent(in) :: istep
    character(len=4), intent(out) :: cstep
    character(len=1) :: aa,bb,cc,dd
    integer :: kk1, kk2, kk3, kk4

    kk1 = istep/1000
    aa  = char(kk1 + 48)
    kk2 = (istep - kk1*1000)/100
    bb  = char(kk2 + 48)
    kk3 = (istep - (kk1*1000) - (kk2*100))/10
    cc  = char(kk3 + 48)
    kk4 = (istep - (kk1*1000) - (kk2*100) - (kk3*10))/1
    dd  = char(kk4 + 48)
    cstep = aa//bb//cc//dd

  end subroutine int2string

!> Get a file unit number free before creating a file
  subroutine sll_new_file_id(file_id, error)
   
    sll_int32, intent(out) :: error     !< error code
    sll_int32, intent(out) :: file_id   !< file unit number
    logical                :: lopen    
      
    error=1

    do 100 file_id=20,99
  
       inquire(unit=file_id,opened=lopen)
       if(lopen) then
          cycle
       else
          open(file_id,status='SCRATCH',err=100)
          close(file_id,status='DELETE',err=100)
          error=0
          exit
       end if
 
    100 continue

    !SLL_ASSERT(error == 0)
   
  end subroutine sll_new_file_id

end module sll_utilities
