subroutine mgdkcyc(work,res,kcur,kcycle,iprer,ipost,iresw, &
                   comm3dp,comm3dl,comm3dc,neighbor,bd,phibc)

use mpi
use mgd3

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

# if cdebug
timing(89)=timing(89)+MPI_WTIME()-tinitial
# endif
return
end
      subroutine mgdkcyc(work,res,kcur,kcycle,iprer,ipost,iresw,
     1                   comm2d,myid,neighbor,bd,phibc,IOUT)
# include "compdir.inc"
      include "mpif.h"
      integer kcur,kcycle,iprer,ipost,iresw,IOUT
      integer comm2d,myid,neighbor(8),bd(8)
      REALN work(*),res(*),phibc(4,*)
c 
      integer nxk,nyk,sxk,exk,syk,eyk,kpbgn,kcbgn
      integer ikdatatype,jkdatatype,ijkdatatype
      integer sxi,exi,syi,eyi
      integer nxr,nyr,sxr,exr,syr,eyr
      integer irdatatype,jrdatatype,ijrdatatype
      common/mgd/nxk(20),nyk(20),sxk(20),exk(20),syk(20),eyk(20),
     1           kpbgn(20),kcbgn(20),ikdatatype(20),jkdatatype(20),
     2           ijkdatatype(20),sxi(20),exi(20),syi(20),eyi(20),
     3           nxr(20),nyr(20),sxr(20),exr(20),syr(20),eyr(20),
     4           irdatatype(20),jrdatatype(20),ijrdatatype(20)
c------------------------------------------------------------------------
c Do one multigrid K-cycle
c K=1 -> V-cycle
c K=2 -> W-cycle
c
c Code      : mgd2, 2-D parallel multigrid solver
c Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
c Called in : mgdsolver
c Calls     : mgdrelax, mgdrestr, mgdcor
c------------------------------------------------------------------------
      integer sxf,exf,syf,eyf,nxf,nyf,ipf,icf
      integer sxc,exc,syc,eyc,nxc,nyc,ipc,irc
      integer klevel,itype,jtype,ijtype,kount(20),l,nrel
      integer sx1,ex1,sy1,ey1
# if cdebug
      double precision tinitial
      tinitial=MPI_WTIME()
# endif
c
      klevel=kcur
c
c pre-relax at current finest grid level
c
      sxf=sxk(klevel)
      exf=exk(klevel)
      syf=syk(klevel)
      eyf=eyk(klevel)
      nxf=nxk(klevel)
      nyf=nyk(klevel)
      ipf=kpbgn(klevel)
      icf=kcbgn(klevel)
      itype=ikdatatype(klevel)
      jtype=jkdatatype(klevel)
      call mgdrelax(sxf,exf,syf,eyf,work(ipf),work(icf),iprer,
     1              comm2d,myid,neighbor,bd,phibc(1,klevel),
     2              itype,jtype,IOUT)
c
c if at coarsest level, post-relax
c
      if (kcur.eq.1) goto 5
c
c calculate residual and restrict it to kcur-1 
c (note: the residuals are stored in the memory space used by the
c the rhs array rhsf, which is therefore erased)
c
      sxc=sxk(klevel-1)
      exc=exk(klevel-1)
      syc=syk(klevel-1)
      eyc=eyk(klevel-1)
      nxc=nxk(klevel-1)
      nyc=nyk(klevel-1)
      ipc=kpbgn(klevel-1)
      irc=kcbgn(klevel-1)+5*(exc-sxc+3)*(eyc-syc+3)
      itype=ikdatatype(klevel)
      jtype=jkdatatype(klevel)
      ijtype=ijkdatatype(klevel)
      call mgdrestr(sxc,exc,syc,eyc,nxc,nyc,work(ipc),work(irc),
     1              sxf,exf,syf,eyf,nxf,nyf,work(ipf),work(icf),
     2              res,iresw,comm2d,myid,neighbor,bd,
     3              itype,jtype,ijtype,IOUT)
c
c set counter for grid levels to zero
c
      do l=1,kcur
        kount(l)=0
      end do
c
c set new level and continue K-cycling
c
      klevel=kcur-1
      nrel=iprer
c
c------------------------------------------------------------------------
c K-cycle control point
c
10    continue
c
c post-relax when kcur revisited
c
      if (klevel.eq.kcur) goto 5
c
c count "hit" at current level
c
      kount(klevel)=kount(klevel)+1
c
c relax at current level
c
      sxf=sxk(klevel)
      exf=exk(klevel)
      syf=syk(klevel)
      eyf=eyk(klevel)
      nxf=nxk(klevel)
      nyf=nyk(klevel)
      ipf=kpbgn(klevel)
      icf=kcbgn(klevel)
      itype=ikdatatype(klevel)
      jtype=jkdatatype(klevel)
      call mgdrelax(sxf,exf,syf,eyf,work(ipf),work(icf),nrel,
     1              comm2d,myid,neighbor,bd,phibc(1,klevel),
     2              itype,jtype,IOUT)

      if (kount(klevel).eq.(kcycle+1)) then
c
c K-cycle(iprer,ipost) complete at klevel
c inject correction to finer grid
c      
        sxf=sxk(klevel+1)
        exf=exk(klevel+1)
        syf=syk(klevel+1)
        eyf=eyk(klevel+1)
        nxf=nxk(klevel+1)
        nyf=nyk(klevel+1)
        ipf=kpbgn(klevel+1)
        sxc=sxk(klevel)
        exc=exk(klevel)
        syc=syk(klevel)
        eyc=eyk(klevel)
        nxc=nxk(klevel)
        nyc=nyk(klevel)
        ipc=kpbgn(klevel)
        sx1=sxi(klevel)
        ex1=exi(klevel)
        sy1=syi(klevel)
        ey1=eyi(klevel)
        call mgdcor(sxf,exf,syf,eyf,nxf,nyf,work(ipf),
     1              sxc,exc,syc,eyc,nxc,nyc,work(ipc),
     2              sx1,ex1,sy1,ey1,bd,phibc(1,klevel+1),IOUT)
c
c reset counter to zero at klevel
c
        kount(klevel)=0
c
c ascend to next higher level and set to post-relax there
c
        klevel=klevel+1
        nrel=ipost
        goto 10
        
      else
c
c K-cycle not complete so descend unless at coarsest
c
        if (klevel.gt.1) then
          sxc=sxk(klevel-1)
          exc=exk(klevel-1)
          syc=syk(klevel-1)
          eyc=eyk(klevel-1)
          nxc=nxk(klevel-1)
          nyc=nyk(klevel-1)
          ipc=kpbgn(klevel-1)
          irc=kcbgn(klevel-1)+5*(exc-sxc+3)*(eyc-syc+3)
          itype=ikdatatype(klevel)
          jtype=jkdatatype(klevel)
          ijtype=ijkdatatype(klevel)
          call mgdrestr(sxc,exc,syc,eyc,nxc,nyc,work(ipc),work(irc),
     1                  sxf,exf,syf,eyf,nxf,nyf,work(ipf),work(icf),
     2                  res,iresw,comm2d,myid,neighbor,bd,
     3                  itype,jtype,ijtype,IOUT)
c
c pre-relax at next coarser level
c
          klevel=klevel-1
          nrel=iprer
          goto 10

        else
c
c post-relax at coarsest level
c
          sxf=sxk(klevel)
          exf=exk(klevel)
          syf=syk(klevel)
          eyf=eyk(klevel)
          ipf=kpbgn(klevel)
          icf=kcbgn(klevel)
          itype=ikdatatype(klevel)
          jtype=jkdatatype(klevel)
          call mgdrelax(sxf,exf,syf,eyf,work(ipf),work(icf),ipost,
     1                  comm2d,myid,neighbor,bd,phibc(1,klevel),
     2                  itype,jtype,IOUT)
c
c inject correction to grid level 2
c
          sxf=sxk(2)
          exf=exk(2)
          syf=syk(2)
          eyf=eyk(2)
          nxf=nxk(2)
          nyf=nyk(2)
          ipf=kpbgn(2)
          sxc=sxk(1)
          exc=exk(1)
          syc=syk(1)
          eyc=eyk(1)
          nxc=nxk(1)
          nyc=nyk(1)
          ipc=kpbgn(1)
          sx1=sxi(1)
          ex1=exi(1)
          sy1=syi(1)
          ey1=eyi(1)
          call mgdcor(sxf,exf,syf,eyf,nxf,nyf,work(ipf),
     1                sxc,exc,syc,eyc,nxc,nyc,work(ipc),
     2                sx1,ex1,sy1,ey1,bd,phibc(1,2),IOUT)
c
c set to post-relax at level 2
c
          nrel=ipost
          klevel=2
          goto 10

        end if

      end if

5     continue
c
c------------------------------------------------------------------------
c post-relax at kcur level
c
      sxf=sxk(kcur)
      exf=exk(kcur)
      syf=syk(kcur)
      eyf=eyk(kcur)
      ipf=kpbgn(kcur)
      icf=kcbgn(kcur)
      itype=ikdatatype(kcur)
      jtype=jkdatatype(kcur)
      call mgdrelax(sxf,exf,syf,eyf,work(ipf),work(icf),ipost,
     1              comm2d,myid,neighbor,bd,phibc(1,kcur),
     2              itype,jtype,IOUT)
c
# if cdebug
      timing(89)=timing(89)+MPI_WTIME()-tinitial
# endif
      return
      end
