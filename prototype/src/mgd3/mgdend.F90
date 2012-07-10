subroutine mgdend(ngrid)
use mgd3
# include "compdir.inc"
include "mpif.h"

integer :: ngrid
!integer :: nxk,nyk,nzk,sxk,exk,syk,eyk,szk,ezk
!integer :: kpbgn,kcbgn,kdatatype
!integer :: sxi,exi,syi,eyi,szi,ezi
!integer :: nxr,nyr,nzr,sxr,exr,syr,eyr,szr,ezr
!integer :: rdatatype
!common/mgd/nxk(20),nyk(20),nzk(20),                          &
!           sxk(20),exk(20),syk(20),eyk(20),szk(20),ezk(20),  &
!           kpbgn(20),kcbgn(20),kdatatype(7,20),              &
!           sxi(20),exi(20),syi(20),eyi(20),szi(20),ezi(20),  &
!           nxr(20),nyr(20),nzr(20),sxr(20),exr(20),syr(20),  &
!           eyr(20),szr(20),ezr(20),rdatatype(7,20)
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
