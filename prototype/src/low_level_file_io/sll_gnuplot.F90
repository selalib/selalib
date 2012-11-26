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
use sll_io
use sll_ascii_io

implicit none

interface sll_gnuplot_field_2d
module procedure sll_gnuplot_rect_2d
end interface

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

   open(file_id,FILE=trim(array_name)//".gnu",FORM='FORMATTED',IOSTAT=error)
   rewind(file_id)

   write(file_id,"(a)")"plot '-' t '"//trim(array_name)//"' with linesp"
   do ipoints = 1, npoints
      write(file_id,"(g15.3)") array(ipoints)
   end do
   write(file_id,"(a)")"e"
   write(file_id,"(a)")"pause -1 'Hit return to quit'"
   close(file_id)

end subroutine sll_gnuplot_write

subroutine sll_gnuplot_rect_2d(xmin, xmax, nx, ymin, ymax, ny, array, array_name, iplot, error)  
sll_real64 :: xmin, xmax, ymin, ymax
sll_int32  :: nx, ny
sll_real64, dimension(nx,ny) :: array
character(len=*) :: array_name
character(len=4) :: fin
sll_int32 :: iplot
sll_int32, save :: gnu_id
sll_int32 :: file_id, error
sll_int32 :: i, j
sll_real64 :: dx, dy, x, y

call int2string(iplot, fin)

if ( iplot == 1 ) then
   call sll_new_file_id(gnu_id, error)
end if

open(gnu_id,file=array_name//".gnu", position="append")
write(gnu_id,*)"splot  '"//array_name//"_"//fin//".dat' w l"
close(gnu_id)

call sll_ascii_file_create(array_name//'_'//fin//'.dat', file_id, error )
dx = (xmax-xmin)/(nx-1)
dy = (ymax-ymin)/(ny-1)
x = xmin
do i=1,nx
   y = ymin
   do j=1,ny
      write(file_id,*) x,y,array(i,j)  
      y = y + dy
   end do
   x = x + dx
   write(file_id,*)
enddo
close(file_id)

end subroutine sll_gnuplot_rect_2d

end module sll_gnuplot
