subroutine mgdcor(sxf,exf,syf,eyf,szf,ezf,nxf,nyf,nzf,phif, &
                  sxc,exc,syc,eyc,szc,ezc,nxc,nyc,nzc,phic, &
                  sx1,ex1,sy1,ey1,sz1,ez1,bd,phibc,IOUT)
use mpi
implicit none
#include "mgd3.h"
integer :: sxf,exf,syf,eyf,szf,ezf,nxf,nyf,nzf
integer :: sxc,exc,syc,eyc,szc,ezc,nxc,nyc,nzc
integer :: sx1,ex1,sy1,ey1,sz1,ez1,bd(26),IOUT
real(8) :: phif(sxf-1:exf+1,syf-1:eyf+1,szf-1:ezf+1)
real(8) :: phic(sxc-1:exc+1,syc-1:eyc+1,szc-1:ezc+1),phibc(6)
!------------------------------------------------------------------------
! Add correction from coarse grid level to fine grid level. Uses
! bilinear interpolation for the old version of the multigrid code,
! and volume weighting for its new version.
!
! Code      : mgd3, 3-D parallel multigrid solver
! Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
! Called in : mgdkcyc
! Calls     : mgdbdry
!------------------------------------------------------------------------
integer i,j,k,ic,jc,kc,i1,i2,j1,j2,k1,k2

# if WMGD
!------------------------------------------------------------------------
! new version: the correction is the weighted average of either two 
! or four points at the coarser grid level depending on whether 
! coarsifying takes place in all directions or not
!
! first, general case where coarsifying takes place in all directions;
! take the averages from the 8 surrounding points
!
if ((nxc.lt.nxf).and.(nyc.lt.nyf).and.(nzc.lt.nzf)) then
  do kc=szc-1,ezc
    k=2*kc-1
    do jc=syc-1,eyc
      j=2*jc-1
      do ic=sxc-1,exc
        i=2*ic-1
        phif(i,j,k)=phif(i,j,k)+                                 &
          (27.0d0*phic(ic,jc,kc)+9.0d0*phic(ic,jc+1,kc)          &
          +9.0d0*phic(ic,jc,kc+1)+3.0d0*phic(ic,jc+1,kc+1)       &
          +9.0d0*phic(ic+1,jc,kc)+3.0d0*phic(ic+1,jc+1,kc)       &
          +3.0d0*phic(ic+1,jc,kc+1)+phic(ic+1,jc+1,kc+1))        &
          /64.0d0
        phif(i+1,j,k)=phif(i+1,j,k)+                             &
          (9.0d0*phic(ic,jc,kc)+3.0d0*phic(ic,jc+1,kc)           &
          +3.0d0*phic(ic,jc,kc+1)+phic(ic,jc+1,kc+1)             &
          +27.0d0*phic(ic+1,jc,kc)+9.0d0*phic(ic+1,jc+1,kc)      &
          +9.0d0*phic(ic+1,jc,kc+1)+3.0d0*phic(ic+1,jc+1,kc+1))  &
          /64.0d0
        phif(i,j+1,k)=phif(i,j+1,k)+                             &
          (9.0d0*phic(ic,jc,kc)+27.0d0*phic(ic,jc+1,kc)          &
          +3.0d0*phic(ic,jc,kc+1)+9.0d0*phic(ic,jc+1,kc+1)       &
          +3.0d0*phic(ic+1,jc,kc)+9.0d0*phic(ic+1,jc+1,kc)       &
          +phic(ic+1,jc,kc+1)+3.0d0*phic(ic+1,jc+1,kc+1))        &
          /64.0d0
        phif(i,j,k+1)=phif(i,j,k+1)+                             &
          (9.0d0*phic(ic,jc,kc)+3.0d0*phic(ic,jc+1,kc)           &
          +27.0d0*phic(ic,jc,kc+1)+9.0d0*phic(ic,jc+1,kc+1)      &
          +3.0d0*phic(ic+1,jc,kc)+phic(ic+1,jc+1,kc)             &
          +9.0d0*phic(ic+1,jc,kc+1)+3.0d0*phic(ic+1,jc+1,kc+1))  &
          /64.0d0
        phif(i+1,j+1,k)=phif(i+1,j+1,k)+                         &
          (3.0d0*phic(ic,jc,kc)+9.0d0*phic(ic,jc+1,kc)           &
          +phic(ic,jc,kc+1)+3.0d0*phic(ic,jc+1,kc+1)             &
          +9.0d0*phic(ic+1,jc,kc)+27.0d0*phic(ic+1,jc+1,kc)      &
          +3.0d0*phic(ic+1,jc,kc+1)+9.0d0*phic(ic+1,jc+1,kc+1))  &
          /64.0d0
        phif(i+1,j,k+1)=phif(i+1,j,k+1)+                         &
          (3.0d0*phic(ic,jc,kc)+phic(ic,jc+1,kc)                 &
          +9.0d0*phic(ic,jc,kc+1)+3.0d0*phic(ic,jc+1,kc+1)       &
          +9.0d0*phic(ic+1,jc,kc)+3.0d0*phic(ic+1,jc+1,kc)       &
          +27.0d0*phic(ic+1,jc,kc+1)+9.0d0*phic(ic+1,jc+1,kc+1)) &
          /64.0d0
        phif(i,j+1,k+1)=phif(i,j+1,k+1)+                         &
          (3.0d0*phic(ic,jc,kc)+9.0d0*phic(ic,jc+1,kc)           &
          +9.0d0*phic(ic,jc,kc+1)+27.0d0*phic(ic,jc+1,kc+1)      &
          +phic(ic+1,jc,kc)+3.0d0*phic(ic+1,jc+1,kc)             &
          +3.0d0*phic(ic+1,jc,kc+1)+9.0d0*phic(ic+1,jc+1,kc+1))  &
          /64.0d0
        phif(i+1,j+1,k+1)=phif(i+1,j+1,k+1)+                     &
          (phic(ic,jc,kc)+3.0d0*phic(ic,jc+1,kc)                 &
          +3.0d0*phic(ic,jc,kc+1)+9.0d0*phic(ic,jc+1,kc+1)       &
          +3.0d0*phic(ic+1,jc,kc)+9.0d0*phic(ic+1,jc+1,kc)       &
          +9.0d0*phic(ic+1,jc,kc+1)+27.0d0*phic(ic+1,jc+1,kc+1))
          /64.0d0
      end do
    end do
  end do
else
!
! no coarsifying in two directions; take the averages from the 2 
! surrounding points
!
  if ((nxf.eq.nxc).and.(nyf.eq.nyc)) then
    do kc=szc-1,ezc
      k=2*kc-1
      do j=syf-1,eyf+1
        jc=j
        do i=sxf-1,exf+1
          ic=i
          phif(i,j,k)=phif(i,j,k)+(3.0d0*phic(ic,jc,kc)+phic(ic,jc,kc+1))/4.0d0
          phif(i,j,k+1)=phif(i,j,k+1)+(phic(ic,jc,kc)+3.0d0*phic(ic,jc,kc+1))/4.0d0
        end do
      end do
    end do
  else if ((nxf.eq.nxc).and.(nzf.eq.nzc)) then
    do k=szf-1,ezf+1
      kc=k
      do jc=syc-1,eyc
        j=2*jc-1
        do i=sxf-1,exf+1
          ic=i
          phif(i,j,k)=phif(i,j,k)+(3.0d0*phic(ic,jc,kc)+phic(ic,jc+1,kc))/4.0d0
          phif(i,j+1,k)=phif(i,j+1,k)+(phic(ic,jc,kc)+3.0d0*phic(ic,jc+1,kc))/4.0d0
        end do
      end do
    end do
  else if ((nyf.eq.nyc).and.(nzf.eq.nzc)) then
    do k=szf-1,ezf+1
      kc=k
      do j=syf-1,eyf+1
        jc=j
        do ic=sxc-1,exc
          i=2*ic-1
          phif(i,j,k)=phif(i,j,k)+(3.0d0*phic(ic,jc,kc)+phic(ic+1,jc,kc))/4.0d0
          phif(i+1,j,k)=phif(i+1,j,k)+(phic(ic,jc,kc)+3.0d0*phic(ic+1,jc,kc))/4.0d0
        end do
      end do
    end do

! no coarsifying in one direction; take the averages from the 4
! surrounding points
!
  else if (nxf.eq.nxc) then
    do kc=szc-1,ezc
      k=2*kc-1
      do jc=syc-1,eyc
        j=2*jc-1
        do i=sxf-1,exf+1
          ic=i
          phif(i,j,k)=phif(i,j,k)+                            &
            (9.0d0*phic(ic,jc,kc)+3.0d0*phic(ic,jc+1,kc)      &
            +3.0d0*phic(ic,jc,kc+1)+phic(ic,jc+1,kc+1))       &
            /16.0d0
          phif(i,j+1,k)=phif(i,j+1,k)+                        &
            (3.0d0*phic(ic,jc,kc)+9.0d0*phic(ic,jc+1,kc)      &
            +phic(ic,jc,kc+1)+3.0d0*phic(ic,jc+1,kc+1))       &
            /16.0d0
          phif(i,j,k+1)=phif(i,j,k+1)+                        &
            (3.0d0*phic(ic,jc,kc)+phic(ic,jc+1,kc)            &
            +9.0d0*phic(ic,jc,kc+1)+3.0d0*phic(ic,jc+1,kc+1)) &
            /16.0d0
          phif(i,j+1,k+1)=phif(i,j+1,k+1)+                    &
            (phic(ic,jc,kc)+3.0d0*phic(ic,jc+1,kc)            &
            +3.0d0*phic(ic,jc,kc+1)+9.0d0*phic(ic,jc+1,kc+1)) &
            /16.0d0
        end do
      end do
    end do
  else if (nyf.eq.nyc) then
    do kc=szc-1,ezc
      k=2*kc-1
      do j=syf-1,eyf+1
        jc=j
        do ic=sxc-1,exc
          i=2*ic-1
          phif(i,j,k)=phif(i,j,k)+                            &
            (9.0d0*phic(ic,jc,kc)+3.0d0*phic(ic+1,jc,kc)      &
            +3.0d0*phic(ic,jc,kc+1)+phic(ic+1,jc,kc+1))       &
            /16.0d0
          phif(i+1,j,k)=phif(i+1,j,k)+                        &
            (3.0d0*phic(ic,jc,kc)+9.0d0*phic(ic+1,jc,kc)      &
            +phic(ic,jc,kc+1)+3.0d0*phic(ic+1,jc,kc+1))       &
            /16.0d0
          phif(i,j,k+1)=phif(i,j,k+1)+                        &
            (3.0d0*phic(ic,jc,kc)+phic(ic+1,jc,kc)            &
            +9.0d0*phic(ic,jc,kc+1)+3.0d0*phic(ic+1,jc,kc+1)) &
            /16.0d0
          phif(i+1,j,k+1)=phif(i+1,j,k+1)+                    &
            (phic(ic,jc,kc)+3.0d0*phic(ic+1,jc,kc)            &
            +3.0d0*phic(ic,jc,kc+1)+9.0d0*phic(ic+1,jc,kc+1)) &
            /16.0d0
        end do
      end do
    end do
  else if (nzf.eq.nzc) then
    do k=szf-1,ezf+1
      kc=k
      do jc=syc-1,eyc
        j=2*jc-1
        do ic=sxc-1,exc
          i=2*ic-1
          phif(i,j,k)=phif(i,j,k)+                            &
            (9.0d0*phic(ic,jc,kc)+3.0d0*phic(ic+1,jc,kc)      &
            +3.0d0*phic(ic,jc+1,kc)+phic(ic+1,jc+1,kc))       &
            /16.0d0
          phif(i+1,j,k)=phif(i+1,j,k)+                        &
            (3.0d0*phic(ic,jc,kc)+9.0d0*phic(ic+1,jc,kc)      &
            +phic(ic,jc+1,kc)+3.0d0*phic(ic+1,jc+1,kc))       &
            /16.0d0
          phif(i,j+1,k)=phif(i,j+1,k)+                        &
            (3.0d0*phic(ic,jc,kc)+phic(ic+1,jc,kc)            &
            +9.0d0*phic(ic,jc+1,kc)+3.0d0*phic(ic+1,jc+1,kc)) &
            /16.0d0
          phif(i+1,j+1,k)=phif(i+1,j+1,k)+                    &
            (phic(ic,jc,kc)+3.0d0*phic(ic+1,jc,kc)            &
            +3.0d0*phic(ic,jc+1,kc)+9.0d0*phic(ic+1,jc+1,kc)) &
            /16.0d0
        end do
      end do
    end do
  end if
end if

! impose Neumann and Dirichlet boundary conditions
! TEMP: periodicity is not enforced to save one call to gxch1lin;
! check whether it has an impact or not...
!
      call mgdbdry(sxf,exf,syf,eyf,szf,ezf,phif,bd,phibc,IOUT)
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
if (nzc.lt.nzf) then
  k1=1
  k2=0
else
  k1=0
  k2=1
end if
!
! identity at the points of the fine grid which have odd indices 
! in i, j, and k
!
do kc=sz1,ez1
  k=k1*(2*kc-1)+k2*kc
  do jc=sy1,ey1
    j=j1*(2*jc-1)+j2*jc
    do ic=sx1,ex1
      i=i1*(2*ic-1)+i2*ic
      phif(i,j,k)=phif(i,j,k)+phic(ic,jc,kc)
    end do
  end do
end do

! interpolation of the two neighboring values for the points of the 
! fine grid with even index for i and odd indices for j and k
!
if (nxc.lt.nxf) then
  do kc=sz1,ez1
    k=k1*(2*kc-1)+k2*kc
    do jc=sy1,ey1
      j=j1*(2*jc-1)+j2*jc
      do ic=sxc-1,exc
        i=2*ic
        phif(i,j,k)=phif(i,j,k)+0.5d0*(phic(ic,jc,kc)+phic(ic+1,jc,kc))
      end do
    end do
  end do
end if
!
! interpolation of the two neighboring values for the points of the 
! fine grid with even index for j and odd indices for i and k
!
if (nyc.lt.nyf) then
  do kc=sz1,ez1
    k=k1*(2*kc-1)+k2*kc
    do jc=syc-1,eyc
      j=2*jc
      do ic=sx1,ex1
        i=i1*(2*ic-1)+i2*ic
        phif(i,j,k)=phif(i,j,k)+0.5d0*(phic(ic,jc,kc)+phic(ic,jc+1,kc))
      end do
    end do
  end do
end if
!
! interpolation of the two neighboring values for the points of the
! fine grid with even index for k and odd indices for i and j
!
if (nzc.lt.nzf) then
  do kc=szc-1,ezc
    k=2*kc
    do jc=sy1,ey1
      j=j1*(2*jc-1)+j2*jc
      do ic=sx1,ex1
        i=i1*(2*ic-1)+i2*ic
        phif(i,j,k)=phif(i,j,k)+0.5d0*(phic(ic,jc,kc)+phic(ic,jc,kc+1))
      end do
    end do
  end do
end if
!
! interpolation of the four neighboring values for the points of the
! fine grid with even indices for i and j and odd index for k
!
if (nxc.lt.nxf.and.nyc.lt.nyf) then
  do kc=sz1,ez1
    k=k1*(2*kc-1)+k2*kc
    do jc=syc-1,eyc
      j=2*jc
      do ic=sxc-1,exc
        i=2*ic
        phif(i,j,k)=phif(i,j,k)                                  &
                   +0.25d0*(phic(ic,jc,kc)+phic(ic+1,jc,kc)      &
                           +phic(ic,jc+1,kc)+phic(ic+1,jc+1,kc))
      end do
    end do
  end do
end if
!
! interpolation of the four neighboring values for the points of the
! fine grid with even indices for i and k and odd index for j
!
if (nxc.lt.nxf.and.nzc.lt.nzf) then
  do kc=szc-1,ezc
    k=2*kc
    do jc=sy1,ey1
      j=j1*(2*jc-1)+j2*jc
      do ic=sxc-1,exc
        i=2*ic
        phif(i,j,k)=phif(i,j,k)                                  &
                   +0.25d0*(phic(ic,jc,kc)+phic(ic+1,jc,kc)      &
                           +phic(ic,jc,kc+1)+phic(ic+1,jc,kc+1))
      end do
    end do
  end do
end if
!
! interpolation of the four neighboring values for the points of the
! fine grid with even indices for j and k and odd index for i
!
if (nyc.lt.nyf.and.nzc.lt.nzf) then
  do kc=szc-1,ezc
    k=2*kc
    do jc=syc-1,eyc
      j=2*jc
      do ic=sx1,ex1
        i=i1*(2*ic-1)+i2*ic
        phif(i,j,k)=phif(i,j,k)                                  &
                   +0.25d0*(phic(ic,jc,kc)+phic(ic,jc+1,kc)      &
                           +phic(ic,jc,kc+1)+phic(ic,jc+1,kc+1))
      end do
    end do
  end do
end if
!
! interpolation of the eight neighboring values for the points of the
! fine grid with even indices for i, j, and k
!
if (nxc.lt.nxf.and.nyc.lt.nyf.and.nzc.lt.nzf) then
  do kc=szc-1,ezc
    k=2*kc
    do jc=syc-1,eyc
      j=2*jc
      do ic=sxc-1,exc
        i=2*ic
        phif(i,j,k)=phif(i,j,k)                                  &
               +0.125d0*(phic(ic,jc,kc)+phic(ic+1,jc,kc)         &
                        +phic(ic,jc+1,kc)+phic(ic,jc,kc+1)       &
                        +phic(ic+1,jc+1,kc)+phic(ic+1,jc,kc+1)   &
                        +phic(ic,jc+1,kc+1)+phic(ic+1,jc+1,kc+1))
      end do
    end do
  end do
end if
# endif

return
end
