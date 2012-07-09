      subroutine mgdrsetf(sxf,exf,syf,eyf,szf,ezf,rf,r,IOUT)
# include "compdir.inc"
      include "mpif.h"
      integer sxf,exf,syf,eyf,szf,ezf,IOUT
      REALN rf(sxf-1:exf+1,syf-1:eyf+1,szf-1:ezf+1),
     1      r(sxf-1:exf+1,syf-1:eyf+1,szf-1:ezf+1)
c------------------------------------------------------------------------
c Set the fine grid values of the density in the work vector
c
c Code      : mgd3, 3-D parallel multigrid solver
c Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
c Called in : mgdsolver
c Calls     : --
c------------------------------------------------------------------------
      integer i,j,k
# if cdebug
      double precision tinitial
      tinitial=MPI_WTIME()
# endif
c
      do k=szf-1,ezf+1
        do j=syf-1,eyf+1
          do i=sxf-1,exf+1
            rf(i,j,k)=r(i,j,k)
          end do
        end do
      end do
c
# if cdebug
      timing(85)=timing(85)+MPI_WTIME()-tinitial
# endif
      return
      end
