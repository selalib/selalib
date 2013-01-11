subroutine mgdrestr(sxc,exc,syc,eyc,nxc,nyc,phic,rhsc,   &
                    sxf,exf,syf,eyf,nxf,nyf,phif,cof,    &
                    resf,iresw,comm2d,myid,neighbor,bd,  &
                    itype,jtype,ijtype)

use mpi
#include "sll_working_precision.h"
#include "mgd2.h"

sll_int32  :: sxc,exc,syc,eyc,nxc,nyc,iresw
sll_int32  :: sxf,exf,syf,eyf,nxf,nyf
sll_int32  :: comm2d,myid,neighbor(8),bd(8),itype,jtype,ijtype
sll_real64 :: phic(sxc-1:exc+1,syc-1:eyc+1)
sll_real64 :: rhsc(sxc-1:exc+1,syc-1:eyc+1)
sll_real64 :: phif(sxf-1:exf+1,syf-1:eyf+1)
sll_real64 :: resf(sxf-1:exf+1,syf-1:eyf+1)
sll_real64 :: cof(sxf-1:exf+1,syf-1:eyf+1,6)
!------------------------------------------------------------------------
! Calculate the residual and restrict it to the coarser level. In
! the new version, the restriction involves 4 points. In the old
! version, it involves 9 points (5 if half-weighting).
!
! Code      : mgd2, 2-D parallel multigrid solver
! Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
! Called in : mgdkcyc
! Calls     : gxch1lin, gxch1cor
!------------------------------------------------------------------------
sll_int32 :: i,j,isrt,jsrt,iinc,jinc,ic,jc
# if DEBUG
sll_real64 :: tinitial
tinitial=MPI_WTIME()
# endif
!------------------------------------------------------------------------
do jc=syc-1,eyc+1
  do ic=sxc-1,exc+1
    phic(ic,jc)=0.0d0
    rhsc(ic,jc)=0.0d0
  end do
end do
!
! calculate residual
!
do j=syf,eyf
  do i=sxf,exf
    resf(i,j)=cof(i,j,6)-(cof(i,j,1)*phif(i-1,j) &
                         +cof(i,j,2)*phif(i+1,j) &
                         +cof(i,j,3)*phif(i,j-1) &
                         +cof(i,j,4)*phif(i,j+1) &
                         +cof(i,j,5)*phif(i,j))
  end do
end do
# if WMGD
!------------------------------------------------------------------------
! new version: calculate the right-hand side at the coarser grid
! level from the averages of the values at the 4 surrounding points;
! if there is no coarsifying in one direction, only 2 points are
! used; no exchange of boundary data is necessary
!
if (nxc.eq.nxf) then
  do jc=syc,eyc
    j=2*jc-2
    do ic=sxc,exc
      i=ic
      rhsc(ic,jc)=0.5d0*(resf(i,j)+resf(i,j+1))
    end do
  end do
else if (nyc.eq.nyf) then
  do jc=syc,eyc
    j=jc
    do ic=sxc,exc
      i=2*ic-2
      rhsc(ic,jc)=0.5d0*(resf(i,j)+resf(i+1,j))
    end do
  end do
else
  do jc=syc,eyc
    j=2*jc-2
    do ic=sxc,exc
      i=2*ic-2
      rhsc(ic,jc)=0.25d0*(resf(i,j)+resf(i+1,j) &
                   +resf(i,j+1)+resf(i+1,j+1))
    end do
  end do
end if
# else
!------------------------------------------------------------------------
! old version: have to exchange boundary data; if full-weighting, 
! need to exchange also corner points
!
call gxch1lin(resf,comm2d,sxf,exf,syf,eyf,neighbor,bd,itype,jtype)
if (iresw.eq.1) then
  call gxch1cor(resf,comm2d,sxf,exf,syf,eyf,neighbor,bd,ijtype)
end if
!
! restrict it to coarser level
!
if (nxc.lt.nxf) then
  isrt=2*sxc-1
  iinc=2
else
  isrt=sxc
  iinc=1
end if
if (nyc.lt.nyf) then
  jsrt=2*syc-1
  jinc=2
else
  jsrt=syc
  jinc=1
end if

if (iresw.eq.1) then
!
! use full weighting
!
  j=jsrt
  do jc=syc,eyc
    i=isrt
    do ic=sxc,exc
      rhsc(ic,jc)=0.25d0*resf(i,j)                        &
                 +0.125d0*(resf(i+1,j)+resf(i-1,j)        &
                          +resf(i,j+1)+resf(i,j-1))       &
                 +0.0625d0*(resf(i+1,j+1)+resf(i+1,j-1)   &
                           +resf(i-1,j-1)+resf(i-1,j+1))
      i=i+iinc
    end do
    j=j+jinc
  end do
else if (iresw.eq.2) then
!
! use half-weighting
!
  j=jsrt
  do jc=syc,eyc
    i=isrt
    do ic=sxc,exc
      rhsc(ic,jc)=0.5d0*resf(i,j)                     &
                 +0.125d0*(resf(i+1,j)+resf(i-1,j)    &
                          +resf(i,j+1)+resf(i,j-1))
      i=i+iinc
    end do
    j=j+jinc
  end do
end if
# endif

# if DEBUG
timing(91)=timing(91)+MPI_WTIME()-tinitial
# endif

end subroutine
