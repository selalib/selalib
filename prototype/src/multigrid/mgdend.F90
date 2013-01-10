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
