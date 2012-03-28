!------------------------------------------------------------------------------
! SELALIB
!------------------------------------------------------------------------------
!
! MODULE: sll_gnuplot
!
!> @author
!> Pierre Navaro
!>
!
! DESCRIPTION: 
!
!> @brief
!> Implements the functions to write data file plotable by GNUplot
!>
!>@details
!>
!> External links:
!> - http://www.gnuplot.info
!
! REVISION HISTORY:
! DD Mmm YYYY - Initial Version
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!------------------------------------------------------------------------------
module sll_gnuplot
#include "sll_working_precision.h"
#include "sll_assert.h"

implicit none

contains  

subroutine sll_gnuplot_write(array,array_name,error)

   sll_real64, dimension(:), intent(in) :: array
   character(len=*), intent(in)         :: array_name
   sll_int32, intent(out)               :: error
   sll_int32                            :: file_id
   sll_int32                            :: ipoints
   sll_int32                            :: npoints
   logical                              :: lopen
   
   npoints = size(array)

   error=0

   do 100 file_id=20,99

      inquire(unit=file_id,opened=lopen)
      if(lopen) then
         cycle
      else
         open(file_id,status='SCRATCH',err=100)
         close(file_id,status='DELETE',err=100)
         goto 200
      end if
 
   100 continue
   error=1
   200 continue
   error=0

   open(file_id,FILE=array_name//".gnu",FORM='FORMATTED',IOSTAT=error)
   rewind(file_id)

   write(file_id,"(a)")"plot '-' using 1:2 t '"//array_name//"' with linesp"
   do ipoints = 1, npoints
      write(file_id,"(g15.3)") array(ipoints)
   end do
   write(file_id,"(a)")"e"
   write(file_id,"(a)")"pause -1 'Hit return to quit'"
   close(file_id)

end subroutine sll_gnuplot_write

end module sll_gnuplot
