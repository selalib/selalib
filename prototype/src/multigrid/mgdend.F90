subroutine mgdend(ngrid)
use mgd3
#ifdef INTEL
implicit none
include "mpif.h"
#else
use mpi
implicit none
#endif

#include "mgd3.h"

integer :: ngrid

!------------------------------------------------------------------------
! Free the MPI datatypes associated witht the multigrid code
!
! Code      : mgd3, 3-D parallel multigrid solver
! Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
! Called in : main
! Calls     : MPI_TYPE_FREE
!------------------------------------------------------------------------

integer j,k,ierr
# if cdebug
double precision tinitial
tinitial=MPI_WTIME()
# endif

do k=1,ngrid-1
  do j=1,7
    call MPI_TYPE_FREE(kdatatype(j,k),ierr)
# if !WMGD
    call MPI_TYPE_FREE(rdatatype(j,k),ierr)
# endif
  end do
end do

# if cdebug
timing(97)=timing(97)+MPI_WTIME()-tinitial
# endif
return
end
      subroutine mgdend(ngrid)
# include "compdir.inc"
      include "mpif.h"
      integer ngrid
c
c common for multigrid indices and datatypes
c
      integer nxk,nyk,sxk,exk,syk,eyk,kpbgn,kcbgn
      integer ikdatatype,jkdatatype,ijkdatatype
      integer sxi,exi,syi,eyi
      integer nxr,nyr,sxr,exr,syr,eyr
      integer irdatatype,jrdatatype,ijrdatatype
      common/mgd/nxk(20),nyk(20),sxk(20),exk(20),syk(20),eyk(20),
     1           kpbgn(20),kcbgn(20),ikdatatype(20),jkdatatype(20),
     2           ijkdatatype(20),sxi(20),exi(20),syi(20),eyi(20),
     3           nxr(20),nyr(20),sxr(20),exr(20),syr(20),eyr(20),
     4           irdatatype(20),jrdatatype(20),ijrdatatype(20)

c------------------------------------------------------------------------
c Free the MPI datatypes associated witht the multigrid code
c
c Code      : mgd2, 2-D parallel multigrid solver
c Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
c Called in : main
c Calls     : MPI_TYPE_FREE
c------------------------------------------------------------------------
      integer k,ierr
# if cdebug
      double precision tinitial
      tinitial=MPI_WTIME()
# endif
c
      do k=1,ngrid-1
        call MPI_TYPE_FREE(ikdatatype(k),ierr)
        call MPI_TYPE_FREE(jkdatatype(k),ierr)
        call MPI_TYPE_FREE(ijkdatatype(k),ierr)
# if !WMGD
        call MPI_TYPE_FREE(irdatatype(k),ierr)
        call MPI_TYPE_FREE(jrdatatype(k),ierr)
        call MPI_TYPE_FREE(ijrdatatype(k),ierr)
# endif
      end do
c
# if cdebug
      timing(97)=timing(97)+MPI_WTIME()-tinitial
# endif
      return
      end
