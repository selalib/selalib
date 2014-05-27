subroutine gp_plot2d(field, x, y, sx, ex, sy, ey, fieldname)
#include "sll_working_precision.h"
#include "sll_utilities.h"
   use sll_collective

   implicit none

   sll_int32        :: i, j
   sll_int32        :: iproc, prank, psize
   character(len=*) :: fieldname
   character(len=4) :: crank
   sll_int32,  intent(in) :: sx, ex, sy, ey
   sll_real64, intent(in) :: x(sx-1:ex+1,sy-1:ey+1)
   sll_real64, intent(in) :: y(sx-1:ex+1,sy-1:ey+1)
   sll_real64, intent(in) :: field(sx-1:ex+1,sy-1:ey+1)

   !write domains
   psize = sll_get_collective_size(sll_world_collective)
   prank = sll_get_collective_rank(sll_world_collective)
   call int2string(prank, crank)
   open( 80, file = fieldname//crank//".dat" )
      do i=sx,ex
         do j=sy,ey
            write(80,"(4e15.5)") x(i,j), y(i,j), field(i,j)
         end do
         write(80,*) 
      end do
   close(80)
   
   !write master file
   if (prank == 0) then
      open( 90, file = fieldname//'.gnu')
      write(90,"(a)",advance='no')"splot '"//fieldname//"0000.dat' w lines"
      do iproc = 1, psize - 1
         call int2string(iproc, crank)
         write(90,"(a)",advance='no') ",'"//fieldname//crank//".dat' w lines" 
      end do
      write(90,*)
      close(90)
   end if

end subroutine gp_plot2d
