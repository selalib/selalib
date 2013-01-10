subroutine mgdkcyc(work,res,kcur,kcycle,iprer,ipost,iresw, &
                   comm2d,myid,neighbor,bd,phibc)

#include "sll_working_precision.h"
#include "mgd2.h"

sll_int32  :: kcur,kcycle,iprer,ipost,iresw
sll_int32  :: comm2d,myid,neighbor(8),bd(8)
sll_real64 :: work(*),res(*),phibc(4,*)

sll_int32 :: nxk,nyk,sxk,exk,syk,eyk,kpbgn,kcbgn
sll_int32 :: ikdatatype,jkdatatype,ijkdatatype
sll_int32 :: sxi,exi,syi,eyi
sll_int32 :: nxr,nyr,sxr,exr,syr,eyr
sll_int32 :: irdatatype,jrdatatype,ijrdatatype
common/mgd/nxk(20),nyk(20),sxk(20),exk(20),syk(20),eyk(20),   &
           kpbgn(20),kcbgn(20),ikdatatype(20),jkdatatype(20), &
           ijkdatatype(20),sxi(20),exi(20),syi(20),eyi(20),   &
           nxr(20),nyr(20),sxr(20),exr(20),syr(20),eyr(20),   &
           irdatatype(20),jrdatatype(20),ijrdatatype(20)

!------------------------------------------------------------------------
! Do one multigrid K-cycle
! K=1 -> V-cycle
! K=2 -> W-cycle
!
! Code      : mgd2, 2-D parallel multigrid solver
! Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
! Called in : mgdsolver
! Calls     : mgdrelax, mgdrestr, mgdcor
!------------------------------------------------------------------------
sll_int32 :: sxf,exf,syf,eyf,nxf,nyf,ipf,icf
sll_int32 :: sxc,exc,syc,eyc,nxc,nyc,ipc,irc
sll_int32 :: klevel,itype,jtype,ijtype,kount(20),l,nrel
sll_int32 :: sx1,ex1,sy1,ey1
# if cdebug
double precision tinitial
tinitial=MPI_WTIME()
# endif

klevel=kcur
!
! pre-relax at current finest grid level
!
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
call mgdrelax(sxf,exf,syf,eyf,work(ipf),work(icf),iprer, &
              comm2d,myid,neighbor,bd,phibc(1,klevel),   &
              itype,jtype)
!
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
nxc=nxk(klevel-1)
nyc=nyk(klevel-1)
ipc=kpbgn(klevel-1)
irc=kcbgn(klevel-1)+5*(exc-sxc+3)*(eyc-syc+3)
itype=ikdatatype(klevel)
jtype=jkdatatype(klevel)
ijtype=ijkdatatype(klevel)
call mgdrestr(sxc,exc,syc,eyc,nxc,nyc,work(ipc),work(irc), &
              sxf,exf,syf,eyf,nxf,nyf,work(ipf),work(icf), &
              res,iresw,comm2d,myid,neighbor,bd,           &
              itype,jtype,ijtype)
!
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
nxf=nxk(klevel)
nyf=nyk(klevel)
ipf=kpbgn(klevel)
icf=kcbgn(klevel)
itype=ikdatatype(klevel)
jtype=jkdatatype(klevel)
call mgdrelax(sxf,exf,syf,eyf,work(ipf),work(icf),nrel, &
              comm2d,myid,neighbor,bd,phibc(1,klevel),itype,jtype)

if (kount(klevel).eq.(kcycle+1)) then
!
! K-cycle(iprer,ipost) complete at klevel
! inject correction to finer grid
!      
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
  call mgdcor(sxf,exf,syf,eyf,nxf,nyf,work(ipf),         &
              sxc,exc,syc,eyc,nxc,nyc,work(ipc),         &
              sx1,ex1,sy1,ey1,bd,phibc(1,klevel+1))
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

! K-cycle not complete so descend unless at coarsest

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
    call mgdrestr(sxc,exc,syc,eyc,nxc,nyc,work(ipc),work(irc),  &
                  sxf,exf,syf,eyf,nxf,nyf,work(ipf),work(icf),  &
                  res,iresw,comm2d,myid,neighbor,bd,            &
                  itype,jtype,ijtype)
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
    ipf=kpbgn(klevel)
    icf=kcbgn(klevel)
    itype=ikdatatype(klevel)
    jtype=jkdatatype(klevel)
    call mgdrelax(sxf,exf,syf,eyf,work(ipf),work(icf),ipost, &
                  comm2d,myid,neighbor,bd,phibc(1,klevel),   &
                  itype,jtype)
!
! inject correction to grid level 2
!
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
    call mgdcor(sxf,exf,syf,eyf,nxf,nyf,work(ipf),  &
                sxc,exc,syc,eyc,nxc,nyc,work(ipc),  &
                sx1,ex1,sy1,ey1,bd,phibc(1,2))
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
ipf=kpbgn(kcur)
icf=kcbgn(kcur)
itype=ikdatatype(kcur)
jtype=jkdatatype(kcur)

call mgdrelax(sxf,exf,syf,eyf,work(ipf),work(icf),ipost, &
              comm2d,myid,neighbor,bd,phibc(1,kcur),     &
              itype,jtype)

# if cdebug
timing(89)=timing(89)+MPI_WTIME()-tinitial
# endif

end subroutine
