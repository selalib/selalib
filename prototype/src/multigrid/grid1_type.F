      subroutine grid1_type(itype,jtype,ijtype,realtype,sx,ex,sy,ey,
     1                      IOUT)
# include "compdir.inc"
      include "mpif.h"
      integer itype,jtype,ijtype,realtype,sx,ex,sy,ey,IOUT
c------------------------------------------------------------------------
c Define the 3 derived datatypes needed to communicate the boundary
c data of (sx-1:ex+1,sy-1:ey+1) arrays between 'myid' and its 8
c neighbors
c
c Code      : tmgd2
c Called in : mgdinit
c Calls     : MPI_TYPE_CONTIGUOUS, MPI_TYPE_COMMIT, MPI_TYPE_VECTOR
c------------------------------------------------------------------------
      integer ier
# if cdebug
      double precision tinitial
      tinitial=MPI_WTIME()
# endif
c
c datatype for one row
c
      call MPI_TYPE_CONTIGUOUS(ex-sx+1,realtype,itype,ierr)
      call MPI_TYPE_COMMIT(itype,ierr)
c
c datatype for one column
c
      call MPI_TYPE_VECTOR(ey-sy+1,1,ex-sx+3,realtype,jtype,ierr)
      call MPI_TYPE_COMMIT(jtype,ierr)
c
c datatype for one 1*1 corner
c
      call MPI_TYPE_CONTIGUOUS(1,realtype,ijtype,ierr)
      call MPI_TYPE_COMMIT(ijtype,ierr)
c
# if cdebug
      timing(7)=timing(7)+MPI_WTIME()-tinitial
# endif
      return
      end
