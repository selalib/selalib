subroutine gxch1pla(sx,ex,sy,ey,sz,ez,a,comm3d,neighbor, &
                    bd,planetype,req,ireq,IOUT)

use mpi
implicit none 
integer :: sx,ex,sy,ey,sz,ez,comm3d
integer :: neighbor(26),bd(26),planetype(3),req(52),ireq,IOUT
real(8) :: a(sx-1:ex+1,sy-1:ey+1,sz-1:ez+1)
!------------------------------------------------------------------------
! Subroutine to exchange planes of thickness 1 of boundary data between
! 'myid' and the 6 processes in the directions i+, i-, j+, j-, k+, k-. 
! Can be used to communicate scalar (p, r,...) as well as vector 
! quantities ((u,v,w), (ut,vt,wt)). The appropriate datatypes have 
! been defined in the 'type_mpi' subroutine.
!
! plane i
!       |     7     |
!      ---------------                k
!       |           |               ^ 
!       |           |               |
!      5|    myid   |1              |
!       |           |              ------>
!       |           |               |     j
!      ---------------
!       |     3     |
!
! plane i-1                  plane i+1
!       |           |              |           | 
!      ---------------            --------------- 
!       |           |              |           |
!       |           |              |           |
!       |     9     |              |    18     |
!       |           |              |           |
!       |           |              |           |
!      ---------------            ---------------
!       |           |              |           |
!
! Code      : tmgd3, test program for 3D parallel multigrid solver
! Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
! Called in : setdens,gbackin, gadvect, gdrbsor, gprbsor, gveloc, 
!             mgdrelax, mgdrestr, mgdrtrsf, 
! Calls     : MPI_ISEND, MPI_IRECV (non-blocking version)
!             MPI_SENDRECV (blocking version)
!------------------------------------------------------------------------
integer ierr
!--------------------------non-blocking----------------------------------
! planes i=constant
!
! send to 18
!
if (bd(18).eq.0) then
  ireq=ireq+1
  call MPI_ISEND(a(ex,sy,sz),1,planetype(1),neighbor(18), &
                 0,comm3d,req(ireq),ierr)
end if
!
! receive from 9
!
if (bd(9).eq.0) then
  ireq=ireq+1
  call MPI_IRECV(a(sx-1,sy,sz),1,planetype(1),neighbor(9), &
                 0,comm3d,req(ireq),ierr)
end if
!
! send to 9
!
if (bd(9).eq.0) then
  ireq=ireq+1
  call MPI_ISEND(a(sx,sy,sz),1,planetype(1),neighbor(9), &
                 1,comm3d,req(ireq),ierr)
end if
!
! receive from 18
!
if (bd(18).eq.0) then
  ireq=ireq+1
  call MPI_IRECV(a(ex+1,sy,sz),1,planetype(1),neighbor(18), &
                 1,comm3d,req(ireq),ierr)
end if
!------------------------------------------------------------------------
! planes j=constant
!
! send to 1
!
if (bd(1).eq.0) then
  ireq=ireq+1
  call MPI_ISEND(a(sx,ey,sz),1,planetype(2),neighbor(1), &
                 2,comm3d,req(ireq),ierr)
end if
!
! receive from 5
!
if (bd(5).eq.0) then
  ireq=ireq+1
  call MPI_IRECV(a(sx,sy-1,sz),1,planetype(2),neighbor(5), &
                 2,comm3d,req(ireq),ierr)
end if
!
! send to 5
!
if (bd(5).eq.0) then
  ireq=ireq+1
  call MPI_ISEND(a(sx,sy,sz),1,planetype(2),neighbor(5), &
                 3,comm3d,req(ireq),ierr)
end if
!
! receive from 1
!
if (bd(1).eq.0) then
  ireq=ireq+1
  call MPI_IRECV(a(sx,ey+1,sz),1,planetype(2),neighbor(1), &
                 3,comm3d,req(ireq),ierr)
end if
!------------------------------------------------------------------------
! planes k=constant
!
! send to 7
!
if (bd(7).eq.0) then
  ireq=ireq+1
  call MPI_ISEND(a(sx,sy,ez),1,planetype(3),neighbor(7), &
                 4,comm3d,req(ireq),ierr)
end if
!
! receive from 3
!
if (bd(3).eq.0) then
  ireq=ireq+1
  call MPI_IRECV(a(sx,sy,sz-1),1,planetype(3),neighbor(3), &
                 4,comm3d,req(ireq),ierr)
end if
!
! send to 3
!
if (bd(3).eq.0) then
  ireq=ireq+1
  call MPI_ISEND(a(sx,sy,sz),1,planetype(3),neighbor(3), &
                 5,comm3d,req(ireq),ierr)
end if
!
! receive from 7
!
if (bd(7).eq.0) then
  ireq=ireq+1
  call MPI_IRECV(a(sx,sy,ez+1),1,planetype(3),neighbor(7), &
                 5,comm3d,req(ireq),ierr)
end if
return
end
