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

!> @namespace sll_utilities
!> @author Selalib team
!> @brief 
!> Library with some useful numerical utlities
!>
!> - Add  :
!> \code
!> #include "sll_utilities.h"
!> \endcode

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

  interface sll_display
     module procedure display_matrix_2d_integer
     module procedure display_matrix_2d_real
     module procedure display_vector_integer
     module procedure display_vector_real
  end interface sll_display

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

!> Convert an integer < 9999 to a 4 characters string
  subroutine int2string( istep, cstep )
    integer, intent(in) :: istep             !< input integer
    character(len=4), intent(out) :: cstep   !< output string
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

  subroutine display_vector_real(array, real_format)

   sll_real64, dimension(:) :: array
   character(len=*)         :: real_format
   character(len=20)        :: display_format
   sll_int32                :: n
   sll_int32                :: i

   n = size(array,1)

   write(display_format, "('(''|''',i4,a,''' |'')')") n, real_format

   write(*,*)
   write(*,display_format) (array(i), i = 1, n)
   write(*,*)

  end subroutine display_vector_real

  subroutine display_vector_integer(array, integer_format)

   sll_int32, dimension(:) :: array
   character(len=*)        :: integer_format
   character(len=20)       :: display_format
   sll_int32               :: n
   sll_int32               :: i

   n = size(array)

   write(display_format, "('(''|''',i4,a,''' |'')')") n, integer_format

   write(*,*)
   write(*,display_format) (array(i), i = 1, n)
   write(*,*)

  end subroutine display_vector_integer


!> Outputs an error message:
!>   - PRTFIL : unit number for print-out
!>   - SEVRTY : 'W' - Warning 'F' - Fatal
  subroutine display_matrix_2d_real(array, real_format)

   sll_real64, dimension(:,:) :: array
   character(len=*)           :: real_format
   character(len=20)          :: display_format
   sll_int32                  :: n1
   sll_int32                  :: n2
   sll_int32                  :: i
   sll_int32                  :: j

   n1 = size(array,1)
   n2 = size(array,2)

   write(display_format, "('(''|''',i4,a,''' |'')')") n2, real_format

   write(*,*)
   do i = 1, n1
      write(*,display_format) (array(i,j), j = 1, n2)
   end do
   write(*,*)

  end subroutine display_matrix_2d_real

  subroutine display_matrix_2d_integer(array, integer_format)

   sll_int32, dimension(:,:) :: array
   character(len=*)          :: integer_format
   character(len=20)         :: display_format
   sll_int32                 :: n1
   sll_int32                 :: n2
   sll_int32                 :: i
   sll_int32                 :: j

   n1 = size(array,1)
   n2 = size(array,2)

   write(display_format, "('(''|''',i4,a,''' |'')')") n2, integer_format

   write(*,*)
   do i = 1, n1
      write(*,display_format) (array(i,j), j = 1, n2)
   end do
   write(*,*)

  end subroutine display_matrix_2d_integer


!> Outputs an error message:
!>   - PRTFIL : unit number for print-out
!>   - SEVRTY : 'W' - Warning 'F' - Fatal
!>   - WHERE  : in which program or subroutine
!>   - ErrMsg : error message
subroutine errout( prtfil, sevrty, lwhere, iline, ErrMsg )

sll_int32, intent(in) ::  prtfil      !< output file unit number
character(len=1),intent(in) :: sevrty !< "W" or "F" : Warning or Fatal
character(len=*),intent(in) :: lwhere !< subroutine where the error appends
sll_int32       ,intent(in) :: iline  !< line number 
character(len=*),intent(in) :: ErrMsg !< error message
character(len=4)            :: cline  
    
write( prtfil, * )
select case ( sevrty )  !     *** Severity ***
case ( 'W' )
   write(prtfil,"(/10x,a)") '*** WARNING ***'
case ( 'F' )
   write(prtfil,"(/10x,a)") '*** FATAL ERROR ***'
case default
   write(prtfil,"(/10x,a)") '*** FATAL ERROR ***'
   write(prtfil,"(/10x,a)") &
        'Error handler (ERROUT) called with unknown severity level: ', &
        SEVRTY
end select
call int2string(iline,cline)
write( prtfil,"(/10x,a)") &
     'Generated by program or subroutine: ', trim(lwhere)//":"//cline
write( prtfil,"(/10x,a)") trim(ErrMsg)
write( prtfil,"(/10x,a)")
    
! return or stop depending on severity
if ( sevrty == 'W' ) then
   return
else
   stop 'Fatal Error: See print file for details'
end if
    
end subroutine errout

!> Subroutine to open data file for slv2d and
!> create the thf.dat file to write results
subroutine initialize_file(data_file_id, thf_file_id)
  sll_int32 :: data_file_id !< namelist file for slv2d
  sll_int32 :: thf_file_id  !< thf file for energy plot
  sll_int32 :: error
  character(len=72) :: filename
  integer :: IO_stat 
    
  call getarg( 1, filename)

  call sll_new_file_id(data_file_id, error)
  open(data_file_id,file=trim(filename),IOStat=IO_stat)
  if (IO_stat/=0) STOP "Miss argument file.nml"

  call sll_new_file_id(thf_file_id, error)
  open(thf_file_id,file="thf.dat",IOStat=IO_stat)
  if (IO_stat/=0) STOP "erreur d'ouverture du fichier thf.dat"

end subroutine initialize_file
  
!> Routine from slv2d to write diagnostics
subroutine time_history(file_id, desc, fformat, array)
   sll_int32 :: file_id !< file unit number
   character(3) :: desc !< name of the diagnostics
   character(14) :: fformat !< fortran output format
   sll_real64, dimension(:) :: array !< data array
    
   if (desc(1:3)=="thf") then
      !print *,'array', array
      write(file_id,fformat) array
   else
      write(*,*) desc," not recognized"
   endif
    
end subroutine time_history

end module sll_utilities
