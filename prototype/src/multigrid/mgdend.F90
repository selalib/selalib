module mgdend
#include "sll_working_precision.h"
#ifdef INTEL
implicit none
include "mpif.h"
#else
use mpi
implicit none
#endif

contains

subroutine mgdend_3d(ngrid)
use mgd3

#include "mgd3.h"

sll_int32 :: ngrid

!------------------------------------------------------------------------
! Free the MPI datatypes associated witht the multigrid code
!
! Code      : mgd3, 3-D parallel multigrid solver
! Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
! Called in : main
! Calls     : MPI_TYPE_FREE
!------------------------------------------------------------------------

sll_int32 :: j,k,ierr
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

end subroutine


subroutine mgdend_2d(ngrid)
#include "mgd2.h"

sll_int32 :: ngrid
!
! common for multigrid indices and datatypes
!
sll_int32 :: nxk,nyk,sxk,exk,syk,eyk,kpbgn,kcbgn
sll_int32 :: ikdatatype,jkdatatype,ijkdatatype
sll_int32 :: sxi,exi,syi,eyi
sll_int32 :: nxr,nyr,sxr,exr,syr,eyr
sll_int32 :: irdatatype,jrdatatype,ijrdatatype
common/mgd/nxk(20),nyk(20),sxk(20),exk(20),syk(20),eyk(20),   &
           kpbgn(20),kcbgn(20),ikdatatype(20),jkdatatype(20), &
           ijkdatatype(20),sxi(20),exi(20),syi(20),eyi(20),   &
           nxr(20),nyr(20),sxr(20),exr(20),syr(20),eyr(20),   &
           irdatatype(20),jrdatatype(20),ijrdatatype(20)

!------------------------------------------------------------------------
! Free the MPI datatypes associated witht the multigrid code
!
! Code      : mgd2, 2-D parallel multigrid solver
! Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
! Called in : main
! Calls     : MPI_TYPE_FREE
!------------------------------------------------------------------------
sll_int32 :: k,ierr
# if cdebug
double precision tinitial
tinitial=MPI_WTIME()
# endif

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

# if cdebug
timing(97)=timing(97)+MPI_WTIME()-tinitial
# endif

end subroutine

end module mgdend
