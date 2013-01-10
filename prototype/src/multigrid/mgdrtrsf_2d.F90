!------------------------------------------------------------------------
!> For the old version of the multigrid code, transfer values of the 
!> density from a finer to a coarser grid level. It is necessary to
!> exchange the boundary density data because the grid "shifts" to
!> the right as it becomes coarser. (In the new version of the
!> multigrid code, there is no such shift, hence no communication is 
!> needed).
!
!> Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
!------------------------------------------------------------------------
subroutine mgdrtrsf(sxc,exc,syc,eyc,nxc,nyc,rc,            &
                       sxf,exf,syf,eyf,nxf,nyf,rf,            &
                       comm2d,myid,neighbor,bd,itype,jtype)
#include "sll_working_precision.h"
#include "mgd2.h"
sll_int32  :: sxc,exc,syc,eyc,nxc,nyc,sxf,exf,syf,eyf,nxf,nyf
sll_int32  :: comm2d,myid,neighbor(8),bd(8),itype,jtype

sll_real64 :: rc(sxc-1:exc+1,syc-1:eyc+1)
sll_real64 :: rf(sxf-1:exf+1,syf-1:eyf+1)

sll_int32  :: i,j,ic,jc,i1,i2,j1,j2

# if cdebug
double precision tinitial
tinitial=MPI_WTIME()
# endif

if (nxc.lt.nxf) then
  i1=1
  i2=0
else 
  i1=0
  i2=1
end if
if (nyc.lt.nyf) then
  j1=1
  j2=0
else
  j1=0
  j2=1
end if
do jc=syc,eyc
  j=j1*(2*jc-1)+j2*jc
  do ic=sxc,exc
    i=i1*(2*ic-1)+i2*ic
    rc(ic,jc)=rf(i,j)
  end do
end do
!
! exchange the boundary values (need only lines, not corner)
!
call gxch1lin(rc,comm2d,sxc,exc,syc,eyc,neighbor,bd,itype,jtype)

# if cdebug
timing(86)=timing(86)+MPI_WTIME()-tinitial
# endif

end subroutine
