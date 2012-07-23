subroutine mgdrestr(sxc,exc,syc,eyc,szc,ezc,nxc,nyc,nzc,phic,    &
                    rhsc,sxf,exf,syf,eyf,szf,ezf,nxf,nyf,nzf,    &
                    phif,cof,resf,iresw,comm3dp,comm3dl,comm3dc, &
                    neighbor,bd,datatypes,IOUT)

use mpi
implicit none
# include "mgd3.h"
integer :: sxc,exc,syc,eyc,szc,ezc,nxc,nyc,nzc,IOUT
integer :: sxf,exf,syf,eyf,szf,ezf,nxf,nyf,nzf
integer :: iresw,comm3dp,comm3dl,comm3dc
integer :: neighbor(26),bd(26),datatypes(7)
real(8) :: phic(sxc-1:exc+1,syc-1:eyc+1,szc-1:ezc+1)
real(8) :: rhsc(sxc-1:exc+1,syc-1:eyc+1,szc-1:ezc+1)
real(8) :: phif(sxf-1:exf+1,syf-1:eyf+1,szf-1:ezf+1)
real(8) :: resf(sxf-1:exf+1,syf-1:eyf+1,szf-1:ezf+1)
real(8) :: cof(sxf-1:exf+1,syf-1:eyf+1,szf-1:ezf+1,8)
!------------------------------------------------------------------------
! Calculate the residual and restrict it to the coarser level. Full
! weighting (iresw=1) is the best method. A half-weighting method
! is also offered (iresw=2) but I advise against using it as it is
! less robust.
!
! Code      : mgd3, 3-D parallel multigrid solver
! Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
! Called in : mgdkcyc
! Calls     : gxch1pla, gxch1lin, gxch1cor, 
!             MPI_WAITALL (non-blocking version)
!------------------------------------------------------------------------
integer :: i,j,k,isrt,jsrt,ksrt,iinc,jinc,kinc,ic,jc,kc
# if !WMGD
integer ireq,req(52)
integer status(MPI_STATUS_SIZE,52),ierr
# endif
!------------------------------------------------------------------------
do kc=szc-1,ezc+1
  do jc=syc-1,eyc+1
    do ic=sxc-1,exc+1
      phic(ic,jc,kc)=0.0d0
      rhsc(ic,jc,kc)=0.0d0
    end do
  end do
end do
!
! calculate residual
!
do k=szf,ezf
  do j=syf,eyf
    do i=sxf,exf
      resf(i,j,k)=cof(i,j,k,8)-(cof(i,j,k,1)*phif(i-1,j,k)  &
                               +cof(i,j,k,2)*phif(i+1,j,k)  &
                               +cof(i,j,k,3)*phif(i,j-1,k)  &
                               +cof(i,j,k,4)*phif(i,j+1,k)  &
                               +cof(i,j,k,5)*phif(i,j,k-1)  &
                               +cof(i,j,k,6)*phif(i,j,k+1)  &
                               +cof(i,j,k,7)*phif(i,j,k))
    end do
  end do
end do
# if WMGD
!------------------------------------------------------------------------
! new version: calculate the right-hand side at the coarser grid
! level from the averages of the values at the 8 surrounding points;
! if there is no coarsifying in one direction, only 2 or 4 points are
! used; no exchange of boundary data is necessary
!
! first, general case where coarsifying takes place in all directions;
! calculate the value from the 8 surrounding points
!
if ((nxc.lt.nxf).and.(nyc.lt.nyf).and.(nzc.lt.nzf)) then
  do kc=szc,ezc
    k=2*kc-2
    do jc=syc,eyc
      j=2*jc-2
      do ic=sxc,exc
        i=2*ic-2
        rhsc(ic,jc,kc)=0.125d0*(resf(i,j,k)+resf(i,j+1,k)        &
                               +resf(i,j,k+1)+resf(i,j+1,k+1)    &
                               +resf(i+1,j,k)+resf(i+1,j+1,k+1)  &
                               +resf(i+1,j+1,k)+resf(i+1,j,k+1))
      end do
    end do
  end do
else if
!
! no coarsifying in two directions; calculate the value from the 2 
! surrounding points
!
  if ((nxf.eq.nxc).and.(nyf.eq.nyc)) then
    do kc=szc,ezc
      k=2*kc-2
      do jc=syc,eyc
        j=jc
        do ic=sxc,exc
          i=ic
          rhsc(ic,jc,kc)=0.5d0*(resf(i,j,k)+resf(i,j,k+1))
        end do
      end do
    end do
  else if ((nxf.eq.nxc).and.(nzf.eq.nzc)) then
    do kc=szc,ezc
      k=kc
      do jc=syc,eyc
        j=2*jc-2
        do ic=sxc,exc
          i=ic
          rhsc(ic,jc,kc)=0.5d0*(resf(i,j,k)+resf(i,j+1,k))
        end do
      end do
    end do
  else if ((nyf.eq.nyc).and.(nzf.eq.nzc)) then
    do kc=szc,ezc
      k=kc
      do jc=syc,eyc
        j=jc
        do ic=sxc,exc
          i=2*ic-2
          rhsc(ic,jc,kc)=0.5d0*(resf(i,j,k)+resf(i+1,j,k))
        end do
      end do
    end do
!
! no coarsifying in one direction; calculate the value from the 4
! surrounding points
!
  else if (nxf.eq.nxc) then
    do kc=szc,ezc
      k=2*kc-2
      do jc=syc,eyc
        j=2*jc-2
        do ic=sxc,exc
          i=ic
          rhsc(ic,jc,kc)=0.25d0*(resf(i,j,k)+resf(i,j+1,k)
                          +resf(i,j,k+1)+resf(i,j+1,k+1))
        end do
      end do
    end do
  else if (nyf.eq.nyc) then
    do kc=szc,ezc
      k=2*kc-2
      do jc=syc,eyc
        j=jc
        do ic=sxc,exc
          i=2*ic-2
          rhsc(ic,jc,kc)=0.25d0*(resf(i,j,k)+resf(i,j,k+1)
                                +resf(i+1,j,k)+resf(i+1,j,k+1))
        end do
      end do
    end do
  else if (nzf.eq.nzc) then
    do kc=szc,ezc
      k=kc
      do jc=syc,eyc
        j=2*jc-2
        do ic=sxc,exc
          i=2*ic-2
          rhsc(ic,jc,kc)=0.25d0*(resf(i,j,k)+resf(i+1,j,k)
                                +resf(i,j+1,k)+resf(i+1,j+1,k))
        end do
      end do
    end do
  end if
end if
# else
!------------------------------------------------------------------------
! old version: have to exchange boundary data; if full-weighting, 
! need to exchange also lines and corners
!

ireq=0

call gxch1pla(sxf,exf,syf,eyf,szf,ezf,resf,comm3dp,neighbor,   &
              bd,datatypes(1),req,ireq,IOUT)
if (iresw.eq.1) then
  call gxch1lin(sxf,exf,syf,eyf,szf,ezf,resf,comm3dl,neighbor, &
                bd,datatypes(4),req,ireq,IOUT)
  call gxch1cor(sxf,exf,syf,eyf,szf,ezf,resf,comm3dc,neighbor, &
                bd,datatypes(7),req,ireq,IOUT)
end if

call MPI_WAITALL(ireq,req,status,ierr)

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
if (nzc.lt.nzf) then
  ksrt=2*szc-1
  kinc=2
else
  ksrt=szc
  kinc=1
end if

if (iresw.eq.1) then
!
! use full weighting
!
  k=ksrt
  do kc=szc,ezc
    j=jsrt
    do jc=syc,eyc
      i=isrt
      do ic=sxc,exc
        rhsc(ic,jc,kc)=0.125d0*resf(i,j,k)                         &
              +0.0625d0*(resf(i+1,j,k)+resf(i-1,j,k)               &
                        +resf(i,j+1,k)+resf(i,j-1,k)               &
                        +resf(i,j,k+1)+resf(i,j,k-1))              &
              +0.03125d0*(resf(i+1,j+1,k)+resf(i+1,j-1,k)          &
                         +resf(i+1,j,k+1)+resf(i+1,j,k-1)          &
                         +resf(i,j+1,k+1)+resf(i,j+1,k-1)          &
                         +resf(i,j-1,k+1)+resf(i,j-1,k-1)          &
                         +resf(i-1,j+1,k)+resf(i-1,j-1,k)          &
                         +resf(i-1,j,k+1)+resf(i-1,j,k-1))         &
              +0.015625d0*(resf(i+1,j+1,k+1)+resf(i+1,j+1,k-1)     &
                          +resf(i+1,j-1,k+1)+resf(i+1,j-1,k-1)     &
                          +resf(i-1,j+1,k+1)+resf(i-1,j+1,k-1)     &
                          +resf(i-1,j-1,k+1)+resf(i-1,j-1,k-1))
        i=i+iinc
      end do
      j=j+jinc
    end do
    k=k+kinc
  end do
else if (iresw.eq.2) then

!
! use half-weighting
!

  k=ksrt
  do kc=szc,ezc
    j=jsrt
    do jc=syc,eyc
      i=isrt
      do ic=sxc,exc
        rhsc(ic,jc,kc)=0.5d0*resf(i,j,k)                          &
              +0.25d0*(resf(i+1,j,k)+resf(i-1,j,k)                &
                      +resf(i,j+1,k)+resf(i,j-1,k)                &
                      +resf(i,j,k+1)+resf(i,j,k-1))/3.0d0
        i=i+iinc
      end do
      j=j+jinc
    end do
    k=k+kinc
  end do
end if
# endif
return
end
