!> For the old version of the multigrid code, set the fine grid values 
!> of the density in the work vector
!> Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
subroutine mgdrsetf(sxf,exf,syf,eyf,rf,r)
use mpi
#include "sll_working_precision.h"
#include "mgd2.h"
sll_int32  :: sxf,exf,syf,eyf
sll_real64 :: rf(sxf-1:exf+1,syf-1:eyf+1),r(sxf-1:exf+1,syf-1:eyf+1)
sll_int32  :: i,j
# if DEBUG
sll_real64 :: tinitial
tinitial=MPI_WTIME()
# endif

do j=syf-1,eyf+1
  do i=sxf-1,exf+1
    rf(i,j)=r(i,j)
  end do
end do

# if DEBUG
timing(85)=timing(85)+MPI_WTIME()-tinitial
# endif

end subroutine
