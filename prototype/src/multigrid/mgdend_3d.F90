!> Free the MPI datatypes associated witht the multigrid code
!> Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
subroutine mgdend(ngrid)
use mgd3
use mpi
#include "sll_working_precision.h"
implicit none
#include "mgd3.h"
integer :: ngrid
sll_int32 :: j,k,ierr

do k=1,ngrid-1
  do j=1,7
    call MPI_TYPE_FREE(kdatatype(j,k),ierr)
# if !WMGD
    call MPI_TYPE_FREE(rdatatype(j,k),ierr)
# endif
  end do
end do

return
end
