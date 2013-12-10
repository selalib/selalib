subroutine mgdkcyc(work,res,kcur,kcycle,iprer,ipost,iresw, &
                   comm3dp,comm3dl,comm3dc,neighbor,bd,phibc)

use mpi
use mgd3, only: sxk,exk,syk,eyk,ezk,szk,nxk,nyk,nzk, &
                sxi,exi,syi,eyi,ezi,szi,iout,kpbgn,  &
                kcbgn, kdatatype

implicit none
#include "mgd3.h"
integer :: kcur,kcycle,iprer,ipost,iresw
integer :: comm3dp,comm3dl,comm3dc
integer :: neighbor(26),bd(26)
real(8) :: work(*),res(*),phibc(6,*)

!------------------------------------------------------------------------
! Do one multigrid K-cycle
! K=1 -> V-cycle
! K=2 -> W-cycle
!
! Code      : mgd3, 3-D parallel multigrid solver
! Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
! Called in : mgdsolver
! Calls     : mgdrelax, mgdrestr, mgdcor
!------------------------------------------------------------------------

integer :: sxf,exf,syf,eyf,szf,ezf,nxf,nyf,nzf,ipf,icf
integer :: sxc,exc,syc,eyc,szc,ezc,nxc,nyc,nzc,ipc,irc
integer :: klevel,kount(20),l,nrel
integer :: sx1,ex1,sy1,ey1,sz1,ez1

klevel=kcur

! pre-relax at current finest grid level

sxf=sxk(klevel)
exf=exk(klevel)
syf=syk(klevel)
eyf=eyk(klevel)
szf=szk(klevel)
ezf=ezk(klevel)
nxf=nxk(klevel)
nyf=nyk(klevel)
nzf=nzk(klevel)
ipf=kpbgn(klevel)
icf=kcbgn(klevel)
call mgdrelax(sxf,exf,syf,eyf,szf,ezf,work(ipf),work(icf),iprer, &
              comm3dp,neighbor,bd,phibc(1,klevel),kdatatype(1,klevel),IOUT)

! if at coarsest level, post-relax
!
if (kcur.eq.1) goto 5
!
! calculate residual and restrict it to kcur-1
! (note: the residuals are stored in the memory space used by the
! the rhs array rhsf, which is therefore erased)
!
sxc=sxk(klevel-1)
exc=exk(klevel-1)
syc=syk(klevel-1)
eyc=eyk(klevel-1)
szc=szk(klevel-1)
ezc=ezk(klevel-1)
nxc=nxk(klevel-1)
nyc=nyk(klevel-1)
nzc=nzk(klevel-1)
ipc=kpbgn(klevel-1)
irc=kcbgn(klevel-1)+7*(exc-sxc+3)*(eyc-syc+3)*(ezc-szc+3)
call mgdrestr(sxc,exc,syc,eyc,szc,ezc,nxc,nyc,nzc,work(ipc),  &
              work(irc),sxf,exf,syf,eyf,szf,ezf,nxf,nyf,nzf,  &
              work(ipf),work(icf),res,iresw,comm3dp,comm3dl,  &
              comm3dc,neighbor,bd,kdatatype(1,klevel),IOUT)

! set counter for grid levels to zero
!
do l=1,kcur
  kount(l)=0
end do
!
! set new level and continue K-cycling
!
klevel=kcur-1
nrel=iprer
!
!------------------------------------------------------------------------
! K-cycle control point
!
10    continue
!
! post-relax when kcur revisited
!
if (klevel.eq.kcur) goto 5
!
! count "hit" at current level
!
kount(klevel)=kount(klevel)+1
!
! relax at current level
!
sxf=sxk(klevel)
exf=exk(klevel)
syf=syk(klevel)
eyf=eyk(klevel)
szf=szk(klevel)
ezf=ezk(klevel)
nxf=nxk(klevel)
nyf=nyk(klevel)
nzf=nzk(klevel)
ipf=kpbgn(klevel)
icf=kcbgn(klevel)
call mgdrelax(sxf,exf,syf,eyf,szf,ezf,work(ipf),work(icf),nrel,  &
              comm3dp,neighbor,bd,phibc(1,klevel),kdatatype(1,klevel),IOUT)

if (kount(klevel).eq.(kcycle+1)) then
!
! K-cycle(iprer,ipost) complete at klevel
! inject correction to finer grid
!      
  sxf=sxk(klevel+1)
  exf=exk(klevel+1)
  syf=syk(klevel+1)
  eyf=eyk(klevel+1)
  szf=szk(klevel+1)
  ezf=ezk(klevel+1)
  nxf=nxk(klevel+1)
  nyf=nyk(klevel+1)
  nzf=nzk(klevel+1)
  ipf=kpbgn(klevel+1)
  sxc=sxk(klevel)
  exc=exk(klevel)
  syc=syk(klevel)
  eyc=eyk(klevel)
  szc=szk(klevel)
  ezc=ezk(klevel)
  nxc=nxk(klevel)
  nyc=nyk(klevel)
  nzc=nzk(klevel)
  ipc=kpbgn(klevel)
  sx1=sxi(klevel)
  ex1=exi(klevel)
  sy1=syi(klevel)
  ey1=eyi(klevel)
  sz1=szi(klevel)
  ez1=ezi(klevel)
  call mgdcor(sxf,exf,syf,eyf,szf,ezf,nxf,nyf,nzf,work(ipf),     &
              sxc,exc,syc,eyc,szc,ezc,nxc,nyc,nzc,work(ipc),     &
              sx1,ex1,sy1,ey1,sz1,ez1,bd,phibc(1,klevel+1),IOUT)
!
! reset counter to zero at klevel
!
  kount(klevel)=0
!
! ascend to next higher level and set to post-relax there
!
  klevel=klevel+1
  nrel=ipost
  goto 10
        
else
!
! K-cycle not complete so descend unless at coarsest
!
  if (klevel.gt.1) then
    sxc=sxk(klevel-1)
    exc=exk(klevel-1)
    syc=syk(klevel-1)
    eyc=eyk(klevel-1)
    szc=szk(klevel-1)
    ezc=ezk(klevel-1)
    nxc=nxk(klevel-1)
    nyc=nyk(klevel-1)
    nzc=nzk(klevel-1)
    ipc=kpbgn(klevel-1)
    irc=kcbgn(klevel-1)+7*(exc-sxc+3)*(eyc-syc+3)*(ezc-szc+3)
    call mgdrestr(sxc,exc,syc,eyc,szc,ezc,nxc,nyc,nzc,work(ipc), &
                  work(irc),sxf,exf,syf,eyf,szf,ezf,nxf,nyf,nzf, &
                  work(ipf),work(icf),res,iresw,comm3dp,comm3dl, &
                  comm3dc,neighbor,bd,kdatatype(1,klevel),IOUT)
!
! pre-relax at next coarser level
!
    klevel=klevel-1
    nrel=iprer
    goto 10

  else
!
! post-relax at coarsest level
!
          
    sxf=sxk(klevel)
    exf=exk(klevel)
    syf=syk(klevel)
    eyf=eyk(klevel)
    szf=szk(klevel)
    ezf=ezk(klevel)
    ipf=kpbgn(klevel)
    icf=kcbgn(klevel)
    call mgdrelax(sxf,exf,syf,eyf,szf,ezf,work(ipf),work(icf),  &
                  ipost,comm3dp,neighbor,bd,phibc(1,klevel),    &
                  kdatatype(1,klevel),IOUT)
!
! inject correction to grid level 2
!
    sxf=sxk(2)
    exf=exk(2)
    syf=syk(2)
    eyf=eyk(2)
    szf=szk(2)
    ezf=ezk(2)
    nxf=nxk(2)
    nyf=nyk(2)
    nzf=nzk(2)
    ipf=kpbgn(2)
    sxc=sxk(1)
    exc=exk(1)
    syc=syk(1)
    eyc=eyk(1)
    szc=szk(1)
    ezc=ezk(1)
    nxc=nxk(1)
    nyc=nyk(1)
    nzc=nzk(1)
    ipc=kpbgn(1)
    sx1=sxi(1)
    ex1=exi(1)
    sy1=syi(1)
    ey1=eyi(1)
    sz1=szi(1)
    ez1=ezi(1)
    call mgdcor(sxf,exf,syf,eyf,szf,ezf,nxf,nyf,nzf,work(ipf),  &
                sxc,exc,syc,eyc,szc,ezc,nxc,nyc,nzc,work(ipc),  &
                sx1,ex1,sy1,ey1,sz1,ez1,bd,phibc(1,2),IOUT)
!
! set to post-relax at level 2
!
    nrel=ipost
    klevel=2
    goto 10

  end if

end if

5     continue
!
!------------------------------------------------------------------------
! post-relax at kcur level
!
sxf=sxk(kcur)
exf=exk(kcur)
syf=syk(kcur)
eyf=eyk(kcur)
szf=szk(kcur)
ezf=ezk(kcur)
ipf=kpbgn(kcur)
icf=kcbgn(kcur)
call mgdrelax(sxf,exf,syf,eyf,szf,ezf,work(ipf),work(icf),ipost, &
              comm3dp,neighbor,bd,phibc(1,kcur),kdatatype(1,kcur),IOUT)

return
end
