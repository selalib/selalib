subroutine grid1_type(gtype,realtype,sx,ex,sy,ey,sz,ez,IOUT)

implicit none 
include "mpif.h"
integer :: gtype(7),realtype,sx,ex,sy,ey,sz,ez,IOUT
!------------------------------------------------------------------------
! Define the 7 derived datatypes needed to communicate the boundary 
! data of (sx-1:ex+1,sy-1:ey+1,sz-1:ez+1) arrays between 'myid' and
! its 26 neighbors.
!
! gtype(l): l=1 -> i=const planes
!           l=2 -> j=const planes
!           l=3 -> k=const planes
!           l=4 -> (i=const,j=const) lines
!           l=5 -> (i=const,k=const) lines
!           l=6 -> (j=const,k=const) lines
!           l=7 -> (i=const,j=const,k=const) corners
!
! Code      : tmgd3, test program for 3D parallel multigrid solver
! Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
! Called in : mgdinit
! Calls     : MPI_TYPE_CONTIGUOUS, MPI_TYPE_COMMIT, MPI_TYPE_VECTOR,
!             MPI_TYPE_EXTENT, MPI_TYPE_HVECTOR
!------------------------------------------------------------------------
integer :: i,ierr

!------------------------------------------------------------------------
! datatype for one 1*1 corner (i=const,j=const,k=const)
!
call MPI_TYPE_CONTIGUOUS(1,realtype,gtype(7),ierr)
call MPI_TYPE_COMMIT(gtype(7),ierr)
!------------------------------------------------------------------------
! datatype for one (i=const,j=const) line 
!
call MPI_TYPE_VECTOR(ez-sz+1,1,(ex-sx+3)*(ey-sy+3),realtype,gtype(4),ierr)
call MPI_TYPE_COMMIT(gtype(4),ierr)
!
! datatype for one (i=const,k=const) line
!
call MPI_TYPE_VECTOR(ey-sy+1,1,ex-sx+3,realtype,gtype(5),ierr)
call MPI_TYPE_COMMIT(gtype(5),ierr)
!
! datatype for one (j=const,k=const) line
!
call MPI_TYPE_CONTIGUOUS(ex-sx+1,realtype,gtype(6),ierr)
call MPI_TYPE_COMMIT(gtype(6),ierr)
!------------------------------------------------------------------------
! datatype for one i=const plane
!
call MPI_TYPE_EXTENT(realtype,i,ierr)
call MPI_TYPE_HVECTOR(ez-sz+1,1,(ex-sx+3)*(ey-sy+3)*i,gtype(5),gtype(1),ierr)
call MPI_TYPE_COMMIT(gtype(1),ierr)
!
! datatype for one j=const plane
!
call MPI_TYPE_VECTOR(ez-sz+1,ex-sx+1,(ex-sx+3)*(ey-sy+3),realtype,gtype(2),ierr)
call MPI_TYPE_COMMIT(gtype(2),ierr)
!
! datatype for one k=const plane
!
call MPI_TYPE_VECTOR(ey-sy+1,ex-sx+1,ex-sx+3,realtype,gtype(3),ierr)
call MPI_TYPE_COMMIT(gtype(3),ierr)

return
end
