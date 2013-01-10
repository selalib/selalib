subroutine mgdsetf(sxf,exf,syf,eyf,phi,rhs,phif,rhsf)
#include "sll_working_precision.h"
#include "mgd2.h"
sll_int32  :: sxf,exf,syf,eyf
sll_real64 :: phi(sxf-1:exf+1,syf-1:eyf+1),rhs(sxf-1:exf+1,syf-1:eyf+1)
sll_real64 :: phif(sxf-1:exf+1,syf-1:eyf+1),rhsf(sxf-1:exf+1,syf-1:eyf+1)
!------------------------------------------------------------------------
! Set the fine grid values in the work vector
!
! Code      : mgd2, 2-D parallel multigrid solver
! Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
! Called in : mgdsolver
! Calls     : --
!------------------------------------------------------------------------
sll_int32  :: i,j
# if cdebug
double precision tinitial
tinitial=MPI_WTIME()
# endif

do j=syf-1,eyf+1
  do i=sxf-1,exf+1
    phi(i,j)=phif(i,j)
    rhs(i,j)=rhsf(i,j)
  end do
end do

# if cdebug
timing(88)=timing(88)+MPI_WTIME()-tinitial
# endif

end subroutine
