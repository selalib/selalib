!------------------------------------------------------------------------------
! SELALIB
!------------------------------------------------------------------------------
!
! MODULE: sll_io
!
!> @author
!> Pierre Navaro
!> Edwin
!
! DESCRIPTION: 
!
!> @brief
!> Implements the functions to display selalib objects
!>
!>@details
!> 
!
! REVISION HISTORY:
! 30 03 2012 - Initial Version
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!------------------------------------------------------------------------------
module sll_io
  
#include "sll_working_precision.h"
#include "sll_assert.h"
use sll_xdmf
use sll_gnuplot

#ifdef STDF95
    integer, parameter :: SLL_IO_XDMF = 0, &
                  SLL_IO_VTK  = 1, &
                  SLL_IO_GNUPLOT = 2
#else
  enum, bind(C)
    enumerator :: SLL_IO_XDMF = 0, &
                  SLL_IO_VTK  = 1, &
                  SLL_IO_GNUPLOT = 2
  end enum
#endif

contains

subroutine sll_new_file_id(file_id, error)

   sll_int32, intent(out)               :: error
   sll_int32, intent(out)               :: file_id
   logical                              :: lopen
   
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

   SLL_ASSERT(error == 0)
   

end subroutine sll_new_file_id

  
end module sll_io
