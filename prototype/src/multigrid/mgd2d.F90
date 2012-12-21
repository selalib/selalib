subroutine mgdsolver(isol,sx,ex,sy,ey,phif,rhsf,r,ngrid,work, &
                     maxcy,tolmax,kcycle,iprer,ipost,iresw, &
                     xl,yl,rro,nx,ny,comm2d,myid,neighbor, &
                     bd,phibc,iter,nprscr,IOUT,nerror)

#include "mgd2.h"
include "mpif.h"
integer isol,sx,ex,sy,ey,ngrid,nx,ny,IOUT
integer maxcy,kcycle,iprer,ipost,iresw
REALN phif(sx-1:ex+1,sy-1:ey+1),rhsf(sx-1:ex+1,sy-1:ey+1)
REALN r(sx-1:ex+1,sy-1:ey+1)
REALN work(*),tolmax,xl,yl,rro,phibc(4,20)
integer comm2d,myid,neighbor(8),bd(8),iter,nerror
logical nprscr

integer nxk,nyk,sxk,exk,syk,eyk,kpbgn,kcbgn
integer ikdatatype,jkdatatype,ijkdatatype
integer sxi,exi,syi,eyi
integer nxr,nyr,sxr,exr,syr,eyr
integer irdatatype,jrdatatype,ijrdatatype
common/mgd/nxk(20),nyk(20),sxk(20),exk(20),syk(20),eyk(20),  &
           kpbgn(20),kcbgn(20),ikdatatype(20),jkdatatype(20),  &
           ijkdatatype(20),sxi(20),exi(20),syi(20),eyi(20),  &
           nxr(20),nyr(20),sxr(20),exr(20),syr(20),eyr(20),  &
           irdatatype(20),jrdatatype(20),ijrdatatype(20)

!------------------------------------------------------------------------
! Parallel multigrid solver in 2-D cartesian coordinates for the
! elliptic equation:      div(cof(grad(phif)))=rhsf
!
! isol=1 -> density
! isol=2 -> pressure
! 
! Written for periodic, wall (Neumann), and constant value
! (Dirichlet) BCs. Tested roughly for all these BCs. There are 
! two versions of the multigrid code, which are separated by the 
! compiler directive WMGD set in 'compdir.inc'. The old version 
! (WMGD=0) corresponds to the original, "traditional" grid setup, 
! and works well when all boundary conditions are periodic. When one 
! of the BCs is not periodic, must compile with the new version 
! (WMGD=1), which uses a different grid setup and new restriction 
! and correction operators (see 'mgdinit' for more details). It is 
! less accurate (or ,equivalently, slower to converge) for the case 
! of all-periodic BCs than the old version, but is better able to 
! handle wall BCs.
!
! Notes: - the values of the rhs contained in the array rhsf are
!          transferred to the work vector and the memory of the 
!          rhsf array is then used to store the residuals at the
!          different levels; therefore, it does not preserve its
!          initial values
!        - with the initialization and clean-up routines mgdinit
!          and mgdend, this multigrid code is self-standing
!
! Code      : mgd2, 2-D parallel multigrid solver
! Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
! Called in : main
! Calls     : mgdrpde, mgdpfpde, mgdphpde, mgdrsetf, mgdppde, mgdrtrsf,
!              -> discretize the pde
!             mgdsetf
!               -> set the initial guesses and the right-hand side
!             mgdkcyc, mgderr,
!               -> do the actual cycling
!             gscale, gxch1lin, gxch1cor,
!               -> rescale pressure and density around average values
!------------------------------------------------------------------------
REALN avo,acorr
integer sxf,exf,syf,eyf,nxf,nyf,sxc,exc,syc,eyc,nxc,nyc
integer ipf,irf,ipc,irc,sxm,exm,sym,eym,nxm,nym,ip,ic,kcur
integer itype,jtype,ijtype,lev,ir1,ir2
# if cdebug
double precision tinitial
tinitial=MPI_WTIME()
# endif
!------------------------------------------------------------------------
! discretize pde at all levels
!
if (isol.eq.1) then
!
! density: only have to set geometric factors
!
  do k=1,ngrid
    sxm=sxk(k)
    exm=exk(k)
    sym=syk(k)
    eym=eyk(k)
    nxm=nxk(k)
    nym=nyk(k)
    ic=kcbgn(k)
    call mgdrpde(sxm,exm,sym,eym,nxm,nym,work(ic),xl,yl,bd,IOUT)
  end do
else
!
! pressure: have to do a lot more work. The coefficients in the new
!           and old versions of the multigrid code are located at
!           different places on the grid.
! Note: the code requires that corasifying takes place in all 
! directions between levels ngrid and ngrid-1
!
! determine coefficients at finest grid level (mid-point values) from
! two neighboring density points
!
  sxf=sxk(ngrid)
  exf=exk(ngrid)
  syf=syk(ngrid)
  eyf=eyk(ngrid)
  nxf=nxk(ngrid)
  nyf=nyk(ngrid)
  icf=kcbgn(ngrid)
  call mgdpfpde(sxf,exf,syf,eyf,nxf,nyf,work(icf),r,xl,yl,bd,IOUT)
# if WMGD
!
! new version: determine coefficients at coarser grid levels from
! four neighboring density points; no communication of boundary data
! is involved because of the grid setup, under the condition that
! mod(nx,nxprocs)=mod(ny,nyprocs)=0
!
  do k=ngrid-1,1,-1
    sxm=sxk(k)
    exm=exk(k)
    sym=syk(k)
    eym=eyk(k)
    nxm=nxk(k)
    nym=nyk(k)
    ic=kcbgn(k)
    call mgdphpde(sxm,exm,sym,eym,nxm,nym,work(ic), &
                  sx,ex,sy,ey,nxf,nyf,r,bd,xl,yl,IOUT)
  end do
# else
  if (ngrid.gt.1) then
!
! old version: use two locations ir1 and ir2 in the work vector to
! store the density on the different grid levels. First set r in 
! work(ir1) at the finest grid level; this is used to determine the 
! coefficients at level k=ngrid-1
!
    ir1=kpbgn(ngrid)
    ir2=kcbgn(ngrid)+5*(exf-sxf+3)*(eyf-syf+3)
    sxf=sxr(ngrid-1)
    exf=exr(ngrid-1)
    syf=syr(ngrid-1)
    eyf=eyr(ngrid-1)
    nxf=nxr(ngrid-1)
    nyf=nyr(ngrid-1)
    lev=1
    call mgdrsetf(sxf,exf,syf,eyf,work(ir1),r,IOUT)
!
! for the levels k=ngrid-1,1, calculate the coefficients from the
! densities stored in the arrays work(ir1) and work(ir2)
! for the levels k=ngrid-1,2, transfer to the level below the values
! of the density needed there to determine the coefficients; exchange
! of the boundary density data is necessary
!
    do k=ngrid-1,1,-1
      sxm=sxk(k)
      exm=exk(k)
      sym=syk(k)
      eym=eyk(k)
      nxm=nxk(k)
      nym=nyk(k)
      ic=kcbgn(k)
      if (lev.eq.1) then
        call mgdppde(sxm,exm,sym,eym,nxm,nym,work(ic), &
                     sxf,exf,syf,eyf,work(ir1),xl,yl,bd,IOUT)
      else
        call mgdppde(sxm,exm,sym,eym,nxm,nym,work(ic), &
                     sxf,exf,syf,eyf,work(ir2),xl,yl,bd,IOUT)
      end if
      if (k.gt.1) then
        sxc=sxr(k-1)
        exc=exr(k-1)
        syc=syr(k-1)
        eyc=eyr(k-1)
        nxc=nxr(k-1)
        nyc=nyr(k-1)
        itype=irdatatype(k-1)
        jtype=jrdatatype(k-1)
        if (lev.eq.1) then
          call mgdrtrsf(sxc,exc,syc,eyc,nxc,nyc,work(ir2), &
                        sxf,exf,syf,eyf,nxf,nyf,work(ir1), &
                        comm2d,myid,neighbor,bd,itype,jtype,IOUT)
          lev=2
        else
          call mgdrtrsf(sxc,exc,syc,eyc,nxc,nyc,work(ir1), &
                        sxf,exf,syf,eyf,nxf,nyf,work(ir2), &
                        comm2d,myid,neighbor,bd,itype,jtype,IOUT)
          lev=1
        end if
        sxf=sxc
        exf=exc
        syf=syc
        eyf=eyc
        nxf=nxc
        nyf=nyc
      end if
    end do
  end if
# endif
end if
!------------------------------------------------------------------------
! set phi,rhsf in work at the finest grid level
!
sxf=sxk(ngrid)
exf=exk(ngrid)
syf=syk(ngrid)
eyf=eyk(ngrid)
ipf=kpbgn(ngrid)
irf=kcbgn(ngrid)+5*(exf-sxf+3)*(eyf-syf+3)
call mgdsetf(sxf,exf,syf,eyf,work(ipf),work(irf),phif,rhsf,IOUT)
!------------------------------------------------------------------------
! cycling at kcur=ngrid level
!
kcur=ngrid
do iter=1,maxcy
  call mgdkcyc(work,rhsf,kcur,kcycle,iprer,ipost,iresw, &
               comm2d,myid,neighbor,bd,phibc,IOUT)
  sxm=sxk(ngrid)
  exm=exk(ngrid)
  sym=syk(ngrid)
  eym=eyk(ngrid)
  ip=kpbgn(ngrid)
  call mgderr(relmax,sxm,exm,sym,eym,phif,work(ip),comm2d,IOUT)
  if (relmax.le.tolmax) goto 1000
end do
!------------------------------------------------------------------------
! if not converged in maxcy cycles, issue an error message and quit
!
if (myid.eq.0) write(IOUT,100) maxcy,relmax
100   format('WARNING: failed to achieve convergence in ',i5, &
       ' cycles  error=',e12.5)
nerror=1
return
!------------------------------------------------------------------------
! converged
!
1000  continue
!
! rescale phif
!
if (isol.eq.1) then
  avo=rro
else
  avo=0.0d0
end if
call gscale(sx,ex,sy,ey,phif,avo,acorr,comm2d,nx,ny,IOUT)
!
! exchange boundary data and impose periodic BCs
!
call gxch1lin(phif,comm2d,sx,ex,sy,ey,neighbor,bd, &
              ikdatatype(ngrid),jkdatatype(ngrid),IOUT)
call gxch1cor(phif,comm2d,sxm,exm,sym,eym,neighbor,bd, &
              ijkdatatype(ngrid),IOUT)
# if WMGD
!
! impose wall and Dirichlet BCs
!
call mgdbdry(sx,ex,sy,ey,phif,bd,vbc,IOUT)
# endif
!
if (isol.eq.1) then
  if (nprscr.and.myid.eq.0) write(IOUT,110) relmax,iter,acorr
# if cdebug
  timing(95)=timing(95)+MPI_WTIME()-tinitial
# endif
else
  if (nprscr.and.myid.eq.0) write(IOUT,120) relmax,iter,acorr
# if cdebug
  timing(96)=timing(96)+MPI_WTIME()-tinitial
# endif
end if
!
return
110 format('  R MGD     err=',e8.3,' iters=',i5,' rcorr=',e9.3)
120 format('  P MGD     err=',e8.3,' iters=',i5,' pcorr=',e9.3)
end
subroutine gerr(sx,ex,sy,ey,p,comm2d,wk,hxi,hyi,pi,nx,ny,IOUT)
# include "compdir.inc"
include "mpif.h"
integer sx,ex,sy,ey,comm2d,IOUT,nx,ny
REALN p(sx-1:ex+1,sy-1:ey+1),wk,hxi,hyi,pi
!-----------------------------------------------------------------------
! Calculate the error between the numerical and exact solution to
! the test problem.
!
! Code      : tmgd2
! Called in : main
! Calls     : MPI_ALLREDUCE
!-----------------------------------------------------------------------
integer i,j,ierr
REALN errloc,err,cx,cy,exact
!
! calculate local error
!
cx=2.0d0*pi*wk
cy=2.0d0*pi*wk
errloc=0.0d0
do j=sy,ey
  yj=(float(j)-1.5d0)/hyi
  do i=sx,ex
    xi=(float(i)-1.5d0)/hxi
    exact=sin(cx*xi)*sin(cy*yj)
    errloc=errloc+abs(p(i,j)-exact)
  end do
end do
!
! calculate global error
!
# if double_precision
call MPI_ALLREDUCE(errloc,err,1,MPI_DOUBLE_PRECISION,MPI_SUM,
                   comm2d,ierr)
# else
call MPI_ALLREDUCE(errloc,err,1,MPI_REAL,MPI_SUM,comm2d,ierr)
# endif
write(IOUT,100) errloc/float(nx*ny),err/float(nx*ny)
100   format(/,'Local error: ',e13.6,'  total error: ',e13.6,/)

return
end
subroutine ginit(sx,ex,sy,ey,p,r,f,wk,hxi,hyi,pi,IOUT)
# include "mgd2.h"
include "mpif.h"
integer sx,ex,sy,ey,IOUT
REALN p(sx-1:ex+1,sy-1:ey+1),r(sx-1:ex+1,sy-1:ey+1)
REALN f(sx-1:ex+1,sy-1:ey+1),hxi,hyi,wk,pi
!-----------------------------------------------------------------------
! Initialize the pressure, density, and right-hand side of the
! elliptic equation div(1/r*grad(p))=f
!
! Code      : tmgd2
! Called in : main
! Calls     : --
!-----------------------------------------------------------------------
integer i,j
REALN cnst,cx,cy,xi,yj
!
do j=sy-1,ey+1
  do i=sx-1,ex+1
    p(i,j)=0.0d0
    r(i,j)=1.0d0
    f(i,j)=0.0d0
  end do
end do
cnst=-8.0d0*(pi*wk)**2
cx=2.0d0*pi*wk
cy=2.0d0*pi*wk
do j=sy,ey
  yj=(float(j)-1.5d0)/hyi
  do i=sx,ex
    xi=(float(i)-1.5d0)/hxi
    f(i,j)=cnst*sin(cx*xi)*sin(cy*yj)
  end do
end do
!
return
end

subroutine grid1_type(itype,jtype,ijtype,realtype,sx,ex,sy,ey,IOUT)
# include "compdir.inc"
include "mpif.h"
integer itype,jtype,ijtype,realtype,sx,ex,sy,ey,IOUT
!------------------------------------------------------------------------
! Define the 3 derived datatypes needed to communicate the boundary
! data of (sx-1:ex+1,sy-1:ey+1) arrays between 'myid' and its 8
! neighbors
!
! Code      : tmgd2
! Called in : mgdinit
! Calls     : MPI_TYPE_CONTIGUOUS, MPI_TYPE_COMMIT, MPI_TYPE_VECTOR
!------------------------------------------------------------------------
integer ier
# if cdebug
double precision tinitial
tinitial=MPI_WTIME()
# endif
!
! datatype for one row
!
call MPI_TYPE_CONTIGUOUS(ex-sx+1,realtype,itype,ierr)
call MPI_TYPE_COMMIT(itype,ierr)
!
! datatype for one column
!
call MPI_TYPE_VECTOR(ey-sy+1,1,ex-sx+3,realtype,jtype,ierr)
call MPI_TYPE_COMMIT(jtype,ierr)
!
! datatype for one 1*1 corner
!
call MPI_TYPE_CONTIGUOUS(1,realtype,ijtype,ierr)
call MPI_TYPE_COMMIT(ijtype,ierr)
!
# if cdebug
timing(7)=timing(7)+MPI_WTIME()-tinitial
# endif
return
end
subroutine MPE_DECOMP1D(n,numprocs,myid,s,e)
# include "compdir.inc"
integer n, numprocs, myid, s, e
integer nlocal
integer deficit
!------------------------------------------------------------------------
!  From the MPE library
!  This file contains a routine for producing a decomposition of a 1-d 
!  array when given a number of processors.  It may be used in "direct" 
!  product decomposition.  The values returned assume a "global" domain 
!  in [1:n]
!
! Code      : tmgd2
! Called in : main
! Calls     : --
!------------------------------------------------------------------------
nlocal  = n / numprocs
s	      = myid * nlocal + 1
deficit = mod(n,numprocs)
s	      = s + min(myid,deficit)
if (myid .lt. deficit) then
    nlocal = nlocal + 1
endif
e = s + nlocal - 1
if (e .gt. n .or. myid .eq. numprocs-1) e = n
return
end

subroutine mgdinit(vbc,phibc,ixp,jyq,iex,jey,ngrid,nxp2,nyp2,  &
                   sx,ex,sy,ey,realtype,nxprocs,nyprocs,nwork, &
                   ibdry,jbdry,myid,IOUT,nerror)
# include "compdir.inc"
include "mpif.h"
integer ixp,jyq,iex,jey,ngrid,nxp2,nyp2,sx,ex,sy,ey
integer realtype,nxprocs,nyprocs,nwork,ibdry,jbdry
integer myid,IOUT,nerror
REALN vbc(4),phibc(4,20)

integer nxk,nyk,sxk,exk,syk,eyk,kpbgn,kcbgn
integer ikdatatype,jkdatatype,ijkdatatype
integer sxi,exi,syi,eyi
integer nxr,nyr,sxr,exr,syr,eyr
integer irdatatype,jrdatatype,ijrdatatype
common/mgd/nxk(20),nyk(20),sxk(20),exk(20),syk(20),eyk(20),     &
           kpbgn(20),kcbgn(20),ikdatatype(20),jkdatatype(20),   &
           ijkdatatype(20),sxi(20),exi(20),syi(20),eyi(20),     &
           nxr(20),nyr(20),sxr(20),exr(20),syr(20),eyr(20),     &
           irdatatype(20),jrdatatype(20),ijrdatatype(20)
!------------------------------------------------------------------------
! Initialize the parallel multigrid solver: subdomain indices and
! MPI datatypes.
!
! The multigrid code comes in two versions. With the WMGD compiler
! directive set to 0, the grid setup is vertex-centered:
!
! WMGD=0
! 
!  |------|-----|-----|-----|-----|            fine
!  1      2     3     4     5     6 
!
!  |------------|-----------|-----------|      coarse
!  1            2           3           4
!
! With WMGD set to 1, it is cell-centered:
!
! WMGD=1   
!           |                       |
!        |--|--|-----|-----|-----|--|--|       fine
!        1  |  2     3     4     5  |  6
!           |                       |
!     |-----|-----|-----------|-----|-----|    coarse
!     1     |     2           3     |     4
!           |                       |
!          wall                    wall
!
! For WMGD=0, the restriction and correction operators are standard
! (choice of full or half weighting for the restriction, bilinear
! interpolation for the correction). This works fine for periodic
! boundary conditions. However, when there are Neumann (wall) or
! Dirichlet BCs, this grid setup results in a loss of accuracy near
! the boundaries when the grid is staggered (the discretization of
! the relaxation operator is first-order locally there). With the
! grid setup corresponding to WMGD=1, accuracy remains second-order
! all the time. As the grid gets coarser, it remains centered on the
! domain instead of "shifting to the right". This option works for
! periodic, Neumann, and Dirichlet BCs, although only periodic and
! Neumann BCs have been tested thoroughly. There is one catch, though.
! For a problem with purely periodic BCs, WMGD=0 converges in less
! cycles than WMGD=1 and requires less CPU time (the penalty is
! apparently between 10 and 50%). This can be attributed to the loss
! of accuracy in the restriction and correction operators due to the
! fact that WMGD=0 uses a support of 3 points in each direction 
! whereas WMGD=1 uses only 2 points.
!
! Both versions offer the option to coarsify in one direction and
! not the other, except at the finest grid level, where coarsifying
! MUST take place along all axes. However, it is possible to have
! ngrid=iex=jey=1, with ixp=nx and jyq=ny. In this case, the code
! never enters 'mgdrestr' and 'mgdcor' and all it does is Gauss-Seidel
! iterate at the finest grid level. This can be useful as a preliminary
! check.
!
! Note: some memory could be saved by noting that the cof arrays
! need be dimensioned (sxm:exm,sym:eym) and not
! (sxm-1:exm+1,sym-1:eym+1)... Probably not too difficult to
! make the change
!
! Code      : mgd2, 2-D parallel multigrid solver
! Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
! Called in : main
! Calls     : grid1_type
!------------------------------------------------------------------------
integer i,j,k,nxf,nyf,nxm,nym,kps,sxm,exm,sym,eym,ierr,nxc,nyc
# if cdebug
double precision tinitial
tinitial=MPI_WTIME()
# endif
!------------------------------------------------------------------------
! set /mgd/ variables to zero
!
do k=1,20
  nxk(k)=0
  nyk(k)=0
  sxk(k)=0
  exk(k)=0
  syk(k)=0
  eyk(k)=0
  kpbgn(k)=0
  kcbgn(k)=0
  ikdatatype(k)=MPI_DATATYPE_NULL
  jkdatatype(k)=MPI_DATATYPE_NULL
  ijkdatatype(k)=MPI_DATATYPE_NULL
  irdatatype(k)=MPI_DATATYPE_NULL
  jrdatatype(k)=MPI_DATATYPE_NULL
  ijrdatatype(k)=MPI_DATATYPE_NULL
  sxi(k)=0
  exi(k)=0
  syi(k)=0
  eyi(k)=0
  nxr(k)=0
  nyr(k)=0
  sxr(k)=0
  exr(k)=0
  syr(k)=0
  eyr(k)=0
end do
!------------------------------------------------------------------------
! make a number of checks
!
# if WMGD
!
! check that, for the new version, the number of processes in each
! direction divides the number of points in that direction (the
! determination of the pressure coefficients in 'mgdphpde' and the
! restriction in 'mgdrestr' would not be complete and would require
! doing some inter-process data communication - warning: would
! be very complex because of the pressure coefficients)
!
if (mod(nxp2-2,nxprocs).ne.0) then
  write(IOUT,100) nxp2-2,nxprocs
  nerror=1
  return
end if
if (mod(nyp2-2,nyprocs).ne.0) then
  write(IOUT,110) nyp2-2,nyprocs
  nerror=1
  return
end if
100 format(/,'ERROR in mgdinit: nx=',i3,' is not a multiple of ', &
           'nxprocs=',i3,/,'cannot use the new version of the ', &
           'multigrid code',/)
110 format(/,'ERROR in mgdinit: ny=',i3,' is not a multiple of ', &
           'nyprocs=',i3,/,'cannot use the new version of the ', &
           'multigrid code',/)
# else
!
! check that the old version is not used with non-periodic BCs
!
if (ibdry.ne.0.or.jbdry.ne.0) then
  write(IOUT,120) ibdry,jbdry
  nerror=1
  return
end if
120 format(/,'ERROR in mgdinit: ibdry=',i2,' jbdry=',i2, &
           /,'cannot use the old version of the multigrid code', &
           /,'boundary conditions that are not periodic', &
           /,'-> change compiler directive to 1 in compdir.inc', &
           /,'   and recompile the multigrid code',/)
# endif
!
! check that the dimensions are correct
!
i=ixp*2**(iex-1)+1
if ((nxp2-1).ne.i) then
  write(IOUT,130) nxp2-1,i
  nerror=1
  return
end if
j=jyq*2**(jey-1)+1
if ((nyp2-1).ne.j) then
  write(IOUT,140) nyp2-1,j
  nerror=1
  return
end if
130 format(/,'ERROR in mgdinit: nxp1=',i3,' <> ixp*2**(iex-1)+1=', &
           i3,/,'-> adjust the multigrid parameters ixp and iex', &
           ' in main',/)
140 format(/,'ERROR in mgdinit: nyp1=',i3,' <> jyq*2**(jey-1)+1=', &
           i3,/,'-> adjust the multigrid parameters jyq and jey', &
           ' in main',/)
!
! check that the number of points at the coarser level is not smaller
! than the number of processes in either direction
!
if (ixp.lt.nxprocs) then
  write(IOUT,150) ixp,nxprocs
  nerror=1
  return
end if
if (jyq.lt.nyprocs) then
  write(IOUT,160) jyq,nyprocs
  nerror=1
  return
end if
150 format(/,'ERROR in mgdinit: ixp=',i3,' < nxprocs=',i3,/, &
           ' there must be at least one grid point at the ', &
           'coarsest grid level',/, &
           '-> increase ixp and decrease iex correspondingly', &
           ' in main',/)
160 format(/,'ERROR in mgdinit: jyq=',i3,' < nyprocs=',i3,/, &
           ' there must be at least one grid point at the ', &
           'coarsest grid level',/, &
           '-> increase jyq and decrease jey correspondingly', &
           ' in main',/)
!
! check that coarsifying takes place in all directions at the finest
! grid level
!
if (ngrid.gt.1) then
  if (iex.eq.1) then
    write(IOUT,170) ngrid,iex
    nerror=1
    return
  end if
  if (jey.eq.1) then
    write(IOUT,180) ngrid,jey
    nerror=1
    return
  end if
end if
170 format(/,'ERROR in mgdinit: ngrid=',i3,' iex=',i3, &
           /,'no coarsifying at the finest grid level in x-direction', &
           /,'this is not allowed by the mutligrid code',/)
180 format(/,'ERROR in mgdinit: ngrid=',i3,' jey=',i3, &
           /,'no coarsifying at the finest grid level in y-direction', &
           /,'this is not allowed by the mutligrid code',/)
!------------------------------------------------------------------------
! define all grid levels
! I have adopted the same notations as in Mudpack as far as possible.
! When a confusion was possible, I added a suffix 'm' to the name
! of the variables. For example, nxm is nxp2-1 for the multigrid
! code whereas nx means nxp2-2 in the rest of the code.
!
do k=1,ngrid
  nxk(k)=ixp*2**(max(k+iex-ngrid,1)-1)+1
  nyk(k)=jyq*2**(max(k+jey-ngrid,1)-1)+1
end do
!
! for all grid levels, set the indices of the subdomain the process
! 'myid' work on, as well as the datatypes needed for the exchange
! of boundary data
!
nxf=nxk(ngrid)
nyf=nyk(ngrid)
sxk(ngrid)=sx
exk(ngrid)=ex
syk(ngrid)=sy
eyk(ngrid)=ey
call grid1_type(ikdatatype(ngrid),jkdatatype(ngrid), &
                ijkdatatype(ngrid),realtype,sx,ex,sy,ey,IOUT)
do k=ngrid-1,1,-1
  nxm=nxk(k)
  nym=nyk(k)
  if (nxm.lt.nxf) then
    sxk(k)=sxk(k+1)/2+1
    exk(k)=(exk(k+1)-1)/2+1
  else
    sxk(k)=sxk(k+1)
    exk(k)=exk(k+1)
  end if
  if (nym.lt.nyf) then
    syk(k)=syk(k+1)/2+1
    eyk(k)=(eyk(k+1)-1)/2+1
  else
    syk(k)=syk(k+1)
    eyk(k)=eyk(k+1)
  end if
  nxf=nxm
  nyf=nym
  call grid1_type(ikdatatype(k),jkdatatype(k),ijkdatatype(k), &
                  realtype,sxk(k),exk(k),syk(k),eyk(k),IOUT)
end do
# if xdebug1
!
! print out the indices and determine the size of the MPI messages
! as a rough check
! 
write(IOUT,*) 'size of the multigrid phi-messages'
do k=ngrid,1,-1
  call MPI_TYPE_SIZE(ikdatatype(k),nsiz1,ierr)
  call MPI_TYPE_SIZE(jkdatatype(k),nsiz2,ierr)
  call MPI_TYPE_SIZE(ijkdatatype(k),nsiz3,ierr)
  write(IOUT,*) 'myid: ',myid,' k=',k,' sxk=',sxk(k),' exk=', &
                exk(k),' syk=',syk(k),' eyk=',eyk(k), &
                ' size of datatypes: ',nsiz1,nsiz2,nsiz3
end do
# endif
!
! set work space indices for phi, cof at each grid level
! check that there is sufficient work space
!
kps=1
do k=ngrid,1,-1
  sxm=sxk(k)
  exm=exk(k)
  sym=syk(k)
  eym=eyk(k)
  kpbgn(k)=kps
  kcbgn(k)=kpbgn(k)+(exm-sxm+3)*(eym-sym+3)
  kps=kcbgn(k)+6*(exm-sxm+3)*(eym-sym+3)
end do
if (kps.gt.nwork) then
  write(IOUT,200) kps,nwork,myid
  nerror=1
  return
else
   write(IOUT,210) kps,nwork
end if

200 format(/,'ERROR in mgdinit: not enough work space',/, &
        ' kps=',i10,' nwork=',i10,' myid: ',i3,/, &
         '-> put the formula for nwork in main in ', &
         'comments',/,'   and set nwork to the value of kps',/)
210  format(/,'WARNING in mgdinit: kps=',i10,' nwork=',i10, &
         /,'can optimize the amount of memory needed by ', &
           'the multigrid code',/,'by putting the formula ', &
           'for nwork into comments and setting',/,'nwork ', &
           'to the value of kps',/)
# if WMGD
!------------------------------------------------------------------------
! For the new version of the multigrid code, set the boundary values 
! to be used for the Dirichlet boundaries. It is possible to assign
! 4 different constant values to the 4 different sides. The values are
! assigned at the finest grid level, zero is assigned at all levels
! below
! 
! vbc, phibc:
!
!       -----vbc(4)------ 
!       |               |
!       |               |
!     vbc(3)----------vbc(1)
!       |               |
!       |               |
!       ------vbc(2)-----
!
do j=1,4
  phibc(j,ngrid)=vbc(j)
end do
do k=ngrid-1,1,-1
  do j=1,4
    phibc(j,k)=0.0d0
  end do
end do
# else
!------------------------------------------------------------------------
! set indices for range of ic and jc on coarser grid level which are
! supported on finer grid level, i.e. for which the points
! (x(2*ic-1,2*jc-1),y(2*ic-1,2*jc-1)) are defined in the subdomain
! of process 'myid'. This allows to avoid having any communication
! after the interpolation (or prolongation) step; it should be used
! in that operation only. 
!
! example:
!   a)  - - -|-----------|-----------|-----------|
!                                  exf=8      exf+1=9
!
!       - ---------------|-----------------------|
!                      exc=4                  exc+1=5
!
!   b)  - - -|-----------|
!          exf=5      exf+1=6
!
!       - - -|-----------------------|
!          exc=3                  exc+1=4
!
      
sxi(ngrid)=sxk(ngrid)-1
exi(ngrid)=exk(ngrid)+1
syi(ngrid)=syk(ngrid)-1
eyi(ngrid)=eyk(ngrid)+1
nxf=nxk(ngrid)
nyf=nyk(ngrid)
do k=ngrid-1,1,-1
  nxc=nxk(k)
  nyc=nyk(k)
  if (nxc.lt.nxf) then
    sxi(k)=sxk(k)-1+mod(sxk(k+1),2)
    exi(k)=exk(k)+1-mod(exk(k+1),2)
  else
    sxi(k)=sxk(k)-1
    exi(k)=exk(k)+1
  end if
  if (nyc.lt.nyf) then
    syi(k)=syk(k)-1+mod(syk(k+1),2)
    eyi(k)=eyk(k)+1-mod(eyk(k+1),2)
  else
    syi(k)=syk(k)-1
    eyi(k)=eyk(k)+1
  end if
# if xdebug1
  write(IOUT,*) 'myid: ',myid,' k=',k,' sxi=',sxi(k),' exi=', &
                exi(k),' syi=',syi(k),' eyi=',eyi(k)
# endif
  nxf=nxc
  nyf=nyc
end do
!------------------------------------------------------------------------
! set indices for determining the coefficients in the elliptic
! equation div(cof*grad(P))=rhs. Used only when solving for the
! pressure. When setting these coefficients at level k, need
! the values of the density at midpoints, i.e. at level k+1
! (if coarsifying takes place between the levels k and k+1).
! If coarsifying took place at all levels, the array cof could
! be used as temporary storage space for the densities, with
! cof(*,*,6) at level k+1 giving the values cof(*,*,1->5) at level
! k, and the already defined datatypes could be used for the 
! exchange of the boundary values. However, this does not work
! in case there is coarsifying in one direction between two levels.
!
! Example: - - -|----------|----------|----------|- -  
!               3          4          5          6    \
!                                                      | coarsifying
!                          |                     |     |
!                          |                     |    /
!                          V                     V
!          - - -|---------------------|---------------------|- -
!               2          r          3          r   \      4
!                                                     |
!                          |                     |    | no coarsifying
!                          |                     |    /
!                          V                     V
!          - - -|---------------------|---------------------|- -
!               2          r          3          r          4
!
! At the finest grid level, the coefficients are determined by a 
! special procedure directly from r, so that no value needs to be
! assigned to the indices at level ngrid.
!
do k=ngrid-1,1,-1
  nxf=nxk(k+1)
  nyf=nyk(k+1)
  nxc=nxk(k)
  nyc=nyk(k)
  if (nxc.lt.nxf) then
    sxr(k)=sxk(k+1)
    exr(k)=exk(k+1)
    nxr(k)=nxk(k+1)
    sxm=sxr(k)
    exm=exr(k)
    nxm=nxr(k)
  else
    if (k.eq.(ngrid-1)) then
      write(IOUT,300) ngrid,ngrid-1
300         format('ERROR in mgdinit: no coarsifying between level ', &
             i3,' and level ',i3,/,', the current version of ', &
             'the code cannot cope with that',/, &
             ' -> decrease the value of ixp and/or jyq and',/, &
             '    increase the value of iex and/or jey')
      nerror=1
      return
    end if
    sxr(k)=sxm
    exr(k)=exm
    nxr(k)=nxm
  end if
  if (nyc.lt.nyf) then
    syr(k)=syk(k+1)
    eyr(k)=eyk(k+1)
    nyr(k)=nyk(k+1)
    sym=syr(k)
    eym=eyr(k)
    nym=nyr(k)
  else
    if (k.eq.(ngrid-1)) then
      write(IOUT,300) ngrid,ngrid-1
      nerror=1
      return
    end if
    syr(k)=sym
    eyr(k)=eym
    nyr(k)=nym
  end if
  call grid1_type(irdatatype(k),jrdatatype(k),ijrdatatype(k), &
                  realtype,sxr(k),exr(k),syr(k),eyr(k),IOUT)
end do
# if xdebug1
!
! print out the indices and determine the size of the MPI messages
! as a rough check
! 
write(IOUT,*) 'size of the r-messages'
do k=ngrid-1,1,-1
  call MPI_TYPE_SIZE(irdatatype(k),nsiz1,ierr)
  call MPI_TYPE_SIZE(jrdatatype(k),nsiz2,ierr)
  call MPI_TYPE_SIZE(ijrdatatype(k),nsiz3,ierr)
  write(IOUT,*) 'myid: ',myid,' k=',k,' sxr=',sxr(k),' exr=', &
                exr(k),' syr=',syr(k),' eyr=',eyr(k), &
                ' size of datatypes: ',nsiz1,nsiz2,nsiz3
end do
# endif
# endif

# if cdebug
timing(81)=timing(81)+MPI_WTIME()-tinitial
# endif
return
end

