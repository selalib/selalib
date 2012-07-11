# 1 "grid1_type.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "grid1_type.F"
      subroutine grid1_type(gtype,realtype,sx,ex,sy,ey,sz,ez,IOUT)

# 1 "compdir.inc" 1


      implicit none 
# 3 "grid1_type.F" 2
      include "mpif.h"
      integer gtype(7),realtype,sx,ex,sy,ey,sz,ez,IOUT
c------------------------------------------------------------------------
c Define the 7 derived datatypes needed to communicate the boundary 
c data of (sx-1:ex+1,sy-1:ey+1,sz-1:ez+1) arrays between 'myid' and
c its 26 neighbors.
c
c gtype(l): l=1 -> i=const planes
c           l=2 -> j=const planes
c           l=3 -> k=const planes
c           l=4 -> (i=const,j=const) lines
c           l=5 -> (i=const,k=const) lines
c           l=6 -> (j=const,k=const) lines
c           l=7 -> (i=const,j=const,k=const) corners
c
c Code      : tmgd3, test program for 3D parallel multigrid solver
c Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
c Called in : mgdinit
c Calls     : MPI_TYPE_CONTIGUOUS, MPI_TYPE_COMMIT, MPI_TYPE_VECTOR,
c             MPI_TYPE_EXTENT, MPI_TYPE_HVECTOR
c------------------------------------------------------------------------
      integer i,ierr




c------------------------------------------------------------------------
c datatype for one 1*1 corner (i=const,j=const,k=const)
c
      call MPI_TYPE_CONTIGUOUS(1,realtype,gtype(7),ierr)
      call MPI_TYPE_COMMIT(gtype(7),ierr)
c------------------------------------------------------------------------
c datatype for one (i=const,j=const) line 
c
      call MPI_TYPE_VECTOR(ez-sz+1,1,(ex-sx+3)*(ey-sy+3),realtype,
     1                     gtype(4),ierr)
      call MPI_TYPE_COMMIT(gtype(4),ierr)
c
c datatype for one (i=const,k=const) line
c
      call MPI_TYPE_VECTOR(ey-sy+1,1,ex-sx+3,realtype,gtype(5),ierr)
      call MPI_TYPE_COMMIT(gtype(5),ierr)
c
c datatype for one (j=const,k=const) line
c
      call MPI_TYPE_CONTIGUOUS(ex-sx+1,realtype,gtype(6),ierr)
      call MPI_TYPE_COMMIT(gtype(6),ierr)
c------------------------------------------------------------------------
c datatype for one i=const plane
c
      call MPI_TYPE_EXTENT(realtype,i,ierr)
      call MPI_TYPE_HVECTOR(ez-sz+1,1,(ex-sx+3)*(ey-sy+3)*i,
     1                      gtype(5),gtype(1),ierr)
      call MPI_TYPE_COMMIT(gtype(1),ierr)
c
c datatype for one j=const plane
c
      call MPI_TYPE_VECTOR(ez-sz+1,ex-sx+1,(ex-sx+3)*(ey-sy+3),
     1                     realtype,gtype(2),ierr)
      call MPI_TYPE_COMMIT(gtype(2),ierr)
c
c datatype for one k=const plane
c
      call MPI_TYPE_VECTOR(ey-sy+1,ex-sx+1,ex-sx+3,
     1                     realtype,gtype(3),ierr)
      call MPI_TYPE_COMMIT(gtype(3),ierr)
c



      return
      end
