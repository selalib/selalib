!> Add correction from coarse grid level to fine grid level. Uses
!> bilinear interpolation for the old version of the multigrid code,
!> and area weighting for its new version.
!>
!> Tested for the case where coarsifying takes place in all directions
subroutine mgdcor(sxf,exf,syf,eyf,nxf,nyf,phif,  &
                  sxc,exc,syc,eyc,nxc,nyc,phic,  &
                  sx1,ex1,sy1,ey1,bd,phibc)

#include "sll_working_precision.h"
#include "mgd2.h"

sll_int32  :: sxf,exf,syf,eyf,nxf,nyf
sll_int32  :: sxc,exc,syc,eyc,nxc,nyc,sx1,ex1,sy1,ey1,bd(8)
sll_real64 :: phif(sxf-1:exf+1,syf-1:eyf+1)
sll_real64 :: phic(sxc-1:exc+1,syc-1:eyc+1),phibc(4)

sll_int32 :: i,j,ic,jc,i1,i2,j1,j2
# if cdebug
sll_real64 ::  precision tinitial
tinitial=MPI_WTIME()
# endif

# if WMGD
!------------------------------------------------------------------------
! new version: the correction is the weighted average of either two 
! or four points at the coarser grid level depending on whether 
! coarsifying takes place in all directions or not
!
if (nxf.eq.nxc) then
  do jc=syc-1,eyc
    j=2*jc-1
    do i=sxf-1,exf+1
      ic=i
      phif(i,j)=phif(i,j)+(3.0d0*phic(ic,jc)+phic(ic,jc+1))/4.0d0
      phif(i,j+1)=phif(i,j+1)+(phic(ic,jc)+3.0d0*phic(ic,jc+1))/4.0d0
    end do
  end do
else if (nyf.eq.nyc) then
  do j=syf-1,eyf+1
    jc=j
    do ic=sxc-1,exc
      i=2*ic-1
      phif(i,j)=phif(i,j)+(3.0d0*phic(ic,jc)+phic(ic+1,jc))/4.0d0
      phif(i+1,j)=phif(i+1,j)+(phic(ic,jc)+3.0d0*phic(ic+1,jc))/4.0d0
    end do
  end do
else
  do jc=syc-1,eyc
    j=2*jc-1
    do ic=sxc-1,exc
      i=2*ic-1
      phif(i,j)=phif(i,j)+  &
        (9.0d0*phic(ic,jc)+3.0d0*phic(ic+1,jc)  &
        +3.0d0*phic(ic,jc+1)+phic(ic+1,jc+1))/16.0d0
      phif(i+1,j)=phif(i+1,j)+  &
        (3.0d0*phic(ic,jc)+9.0d0*phic(ic+1,jc)  &
        +phic(ic,jc+1)+3.0d0*phic(ic+1,jc+1))/16.0d0
      phif(i,j+1)=phif(i,j+1)+  &
        (3.0d0*phic(ic,jc)+phic(ic+1,jc)  &
        +9.0d0*phic(ic,jc+1)+3.0d0*phic(ic+1,jc+1))/16.0d0
      phif(i+1,j+1)=phif(i+1,j+1)+  &
        (phic(ic,jc)+3.0d0*phic(ic+1,jc)  &
        +3.0d0*phic(ic,jc+1)+9.0d0*phic(ic+1,jc+1))/16.0d0
    end do
  end do
end if
!
! impose Neumann and Dirichlet boundary conditions
! TEMP: periodicity is not enforced to save one call to gxch1lin;
! check whether it has an impact or not...
!
      call mgdbdry(sxf,exf,syf,eyf,phif,bd,phibc,IOUT)
# else
!------------------------------------------------------------------------
! old version
!
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
!
! identity at the points of the fine grid which have odd indices 
! in i and j
!
do jc=sy1,ey1
  j=j1*(2*jc-1)+j2*jc
  do ic=sx1,ex1
    i=i1*(2*ic-1)+i2*ic
    phif(i,j)=phif(i,j)+phic(ic,jc)
  end do
end do
!
! interpolation of the two neighboring values for the points of the 
! fine grid with even index for i and odd index for j
!
if (nxc.lt.nxf) then
  do jc=sy1,ey1
    j=j1*(2*jc-1)+j2*jc
    do ic=sxc-1,exc
      i=2*ic
      phif(i,j)=phif(i,j)+0.5d0*(phic(ic,jc)+phic(ic+1,jc))
    end do
  end do
end if
!
! interpolation of the two neighboring values for the points of the 
! fine grid with odd index for i and even index for j
!
if (nyc.lt.nyf) then
  do jc=syc-1,eyc
    j=2*jc
    do ic=sx1,ex1
      i=i1*(2*ic-1)+i2*ic
      phif(i,j)=phif(i,j)+0.5d0*(phic(ic,jc)+phic(ic,jc+1))
    end do
  end do
end if
!
! interpolation of the four neighboring values for the points of the
! fine grid with even indices for i and j
!
if (nxc.lt.nxf.and.nyc.lt.nyf) then
  do jc=syc-1,eyc
    j=2*jc
    do ic=sxc-1,exc
      i=2*ic
      phif(i,j)=phif(i,j)+0.25d0*(phic(ic,jc)+phic(ic+1,jc) &
                                 +phic(ic,jc+1)+phic(ic+1,jc+1))
    end do
  end do
end if
# endif

# if cdebug
timing(92)=timing(92)+MPI_WTIME()-tinitial
# endif

end subroutine
