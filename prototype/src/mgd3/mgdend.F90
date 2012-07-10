      subroutine mgdend(ngrid)
# include "compdir.inc"
      include "mpif.h"
      integer ngrid
c
      integer nxk,nyk,nzk,sxk,exk,syk,eyk,szk,ezk
      integer kpbgn,kcbgn,kdatatype
      integer sxi,exi,syi,eyi,szi,ezi
      integer nxr,nyr,nzr,sxr,exr,syr,eyr,szr,ezr
      integer rdatatype
      common/mgd/nxk(20),nyk(20),nzk(20),
     1           sxk(20),exk(20),syk(20),eyk(20),szk(20),ezk(20),
     2           kpbgn(20),kcbgn(20),kdatatype(7,20),
     3           sxi(20),exi(20),syi(20),eyi(20),szi(20),ezi(20),
     4           nxr(20),nyr(20),nzr(20),sxr(20),exr(20),syr(20),
     5           eyr(20),szr(20),ezr(20),rdatatype(7,20)
c------------------------------------------------------------------------
c Free the MPI datatypes associated witht the multigrid code
c
c Code      : mgd3, 3-D parallel multigrid solver
c Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
c Called in : main
c Calls     : MPI_TYPE_FREE
c------------------------------------------------------------------------
      integer j,k,ierr
# if cdebug
      double precision tinitial
      tinitial=MPI_WTIME()
# endif
c
      do k=1,ngrid-1
        do j=1,7
          call MPI_TYPE_FREE(kdatatype(j,k),ierr)
# if !WMGD
          call MPI_TYPE_FREE(rdatatype(j,k),ierr)
# endif
        end do
      end do
c
# if cdebug
      timing(97)=timing(97)+MPI_WTIME()-tinitial
# endif
      return
      end
