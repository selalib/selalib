module diagnostics
use mpi
#include "sll_working_precision.h"
#include "sll_constants.h"
#include "sll_utilities.h"
implicit none


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine gp_plot2d(p, hxi, hyi, sx, ex, sy, ey, wk )

   sll_int32,  intent(in) :: sx, ex, sy, ey
   sll_real64, intent(in), dimension(sx-1:ex+1,sy-1:ey+1) :: p
   sll_real64, intent(in) :: hxi, hyi
   sll_real64, intent(in) :: wk
   sll_int32  :: i, j
   sll_int32  :: ierr
   sll_real64 :: xi, yj, errloc, cx, cy, exact
   sll_int32  :: prank, psize, iproc
   character(len=4) :: crank

   call MPI_COMM_RANK(MPI_COMM_WORLD,prank,ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD,psize,ierr)
   
   cx=2.0d0*sll_pi*wk
   cy=2.0d0*sll_pi*wk

   !write domains
   call int2string(prank, crank)
   open( 80, file = "p"//crank//".dat" )
      do i=sx,ex
         xi=(float(i)-1.5d0)/hxi
         do j=sy,ey
            yj=(float(j)-1.5d0)/hyi
            exact=sin(cx*xi)*sin(cy*yj)
            errloc=abs(p(i,j)-exact)
            write(80,"(4e15.5)") xi, yj, errloc, p(i,j)
         end do
         write(80,*) 
      end do
   close(80)
   
   !write master file
   if (prank == 0) then
      open( 90, file = 'p.gnu')
      write(90,"(a)",advance='no')"splot 'p0000.dat' w lines"
      do iproc = 1, psize - 1
         call int2string(iproc, crank)
         write(90,"(a)",advance='no') ",'p"//crank//".dat' w lines" 
      end do
      write(90,*)
      close(90)
   end if

end subroutine gp_plot2d

end module diagnostics
