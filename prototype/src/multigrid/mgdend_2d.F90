!> Free the MPI datatypes associated witht the multigrid code
subroutine mgdend(ngrid)
#include "sll_working_precision.h"
use sll_mgd2, only:ikdatatype,jkdatatype,ijkdatatype, &
                   irdatatype,jrdatatype,ijrdatatype

#include "mgd2.h"


sll_int32 :: k,ierr
sll_int32 :: ngrid

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

call MPI_FINALIZE(ierr)

end subroutine
