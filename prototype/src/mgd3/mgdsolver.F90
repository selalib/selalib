subroutine mgdsolver(isol,sx,ex,sy,ey,sz,ez,phif,rhsf,r,ngrid, &
                     work,maxcy,tolmax,kcycle,iprer,ipost,     &
                     iresw,xl,yl,zl,rro,nx,ny,nz,comm3d,       &
                     comm3dp,comm3dl,comm3dc,myid,neighbor,    &
                     bd,phibc,iter,nprscr,IOUT,nerror)


use mgd3
# include "compdir.inc"
include "mpif.h"
integer :: sx,ex,sy,ey,sz,ez,ngrid,IOUT
integer :: maxcy,kcycle,iprer,ipost,iresw
real(8) :: phif(sx-1:ex+1,sy-1:ey+1,sz-1:ez+1)
real(8) :: rhsf(sx-1:ex+1,sy-1:ey+1,sz-1:ez+1)
real(8) :: r(sx-1:ex+1,sy-1:ey+1,sz-1:ez+1),rro
real(8) :: work(*),tolmax,xl,yl,zl,phibc(6,20)
integer :: isol,comm3d,comm3dp,comm3dl,comm3dc,myid
integer :: neighbor(26),bd(26),iter,nerror
logical :: nprscr
integer :: k
integer :: nx, ny, nz, icf
real(8) :: relmax

!------------------------------------------------------------------------
! Parallel multigrid solver in 3-D cartesian coordinates for the 
! elliptic equation:      div(cof(grad(phif)))=rhsf
!
! isol=1 -> density
! isol=2 -> pressure
! 
! Written for periodic, wall (Neumann), and constant value (Dirichlet)
! BCs. There are two versions of the multigrid code, which are 
! separated by the compiler directive WMGD set in 'compdir.inc'. 
! The old version (WMGD=0) corresponds to the original, "traditional" 
! grid setup, and works well when all boundary conditions are 
! periodic. When one of the BCs is not periodic, must compile with 
! the new version (WMGD=1), which uses a different grid setup and 
! new restriction and correction operators (see 'mgdinit' for more 
! details). It is less accurate (or ,equivalently, slower to converge) 
! for the case of all-periodic BCs than the old version, but is better 
! able to handle wall BCs.
!
! Notes: - the values of the rhs contained in the array rhsf are
!          transferred to the work vector and the memory of the 
!          rhsf array is then used to store the residuals at the
!          different levels; therefore, it does not preserve its
!          initial values
!        - with the initialization and clean-up routines mgdinit
!          and mgdend, this multigrid code is self-standing
!
! Code      : mgd3, 3-D parallel multigrid solver
! Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
! Called in : main
! Calls     : mgdrpde, mgdpfpde, mgdphpde, mgdrsetf, mgdppde, mgdrtrsf,
!              -> discretize the pde
!             mgdsetf
!               -> set the initial guesses and the right-hand side
!             mgdkcyc, mgderr
!               -> do the actual cycling
!             gscale, gxch1pla, gxch1lin, gxch1cor,
!               -> rescale pressure and density around average values
!             gxch1pla, gxch1lin, gxch1cor, MPI_WAITALL (non-blocking
!             version)
!------------------------------------------------------------------------
real(8) :: avo,acorr
integer :: sxf,exf,syf,eyf,szf,ezf,nxf,nyf,nzf
integer :: sxc,exc,syc,eyc,szc,ezc,nxc,nyc,nzc
integer :: sxm,exm,sym,eym,szm,ezm,nxm,nym,nzm
integer :: ipf,irf,ip,ic,kcur,lev,ir1,ir2
integer :: ireq,req(52)
integer :: status(MPI_STATUS_SIZE,52),ierr

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
    szm=szk(k)
    ezm=ezk(k)
    nxm=nxk(k)
    nym=nyk(k)
    nzm=nzk(k)
    ic=kcbgn(k)
    call mgdrpde(sxm,exm,sym,eym,szm,ezm,nxm,nym,nzm,work(ic),xl,yl,zl,IOUT)
  end do
else
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
  szf=szk(ngrid)
  ezf=ezk(ngrid)
  nxf=nxk(ngrid)
  nyf=nyk(ngrid)
  nzf=nzk(ngrid)
  icf=kcbgn(ngrid)
  call mgdpfpde(sxf,exf,syf,eyf,szf,ezf,nxf,nyf,nzf,work(icf),r,xl,yl,zl,IOUT)
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
    szm=szk(k)
    ezm=ezk(k)
    nxm=nxk(k)
    nym=nyk(k)
    nzm=nzk(k)
    ic=kcbgn(k)
    call mgdphpde(sxm,exm,sym,eym,szm,ezm,nxm,nym,nzm,work(ic),  &
                  sx,ex,sy,ey,sz,ez,nxf,nyf,nzf,r,bd,xl,yl,zl,   &
                  IOUT)
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
    ir2=kcbgn(ngrid)+7*(exf-sxf+3)*(eyf-syf+3)*(ezf-szf+3)
    sxf=sxr(ngrid-1)
    exf=exr(ngrid-1)
    syf=syr(ngrid-1)
    eyf=eyr(ngrid-1)
    szf=szr(ngrid-1)
    ezf=ezr(ngrid-1)
    nxf=nxr(ngrid-1)
    nyf=nyr(ngrid-1)
    nzf=nzr(ngrid-1)
    lev=1
    call mgdrsetf(sxf,exf,syf,eyf,szf,ezf,work(ir1),r,IOUT)
!
! for the levels k=ngrid-1,1, calculate the coefficients from the
! densities stored in the arrays work(ir1) and work(ir2)
! for the levels k=ngrid-1,2, transfer to the level below the values
! of the density needed there to determine the coefficients
!
    do k=ngrid-1,1,-1
      sxm=sxk(k)
      exm=exk(k)
      sym=syk(k)
      eym=eyk(k)
      szm=szk(k)
      ezm=ezk(k)
      nxm=nxk(k)
      nym=nyk(k)
      nzm=nzk(k)
      ic=kcbgn(k)
      if (lev.eq.1) then
        call mgdppde(sxm,exm,sym,eym,szm,ezm,nxm,nym,nzm,work(ic), &
                     sxf,exf,syf,eyf,szf,ezf,work(ir1),            &
                     xl,yl,zl,IOUT)
      else
        call mgdppde(sxm,exm,sym,eym,szm,ezm,nxm,nym,nzm,work(ic), &
                     sxf,exf,syf,eyf,szf,ezf,work(ir2),xl,yl,zl,IOUT)
      end if
      if (k.gt.1) then
        sxc=sxr(k-1)
        exc=exr(k-1)
        syc=syr(k-1)
        eyc=eyr(k-1)
        szc=szr(k-1)
        ezc=ezr(k-1)
        nxc=nxr(k-1)
        nyc=nyr(k-1)
        nzc=nzr(k-1)
        if (lev.eq.1) then
          call mgdrtrsf(sxc,exc,syc,eyc,szc,ezc,nxc,nyc,nzc, &
                        work(ir2),sxf,exf,syf,eyf,szf,ezf,   &
                        nxf,nyf,nzf,work(ir1),comm3dp,myid,  &
                        neighbor,bd,rdatatype(1,k-1),IOUT)
          lev=2
        else
          call mgdrtrsf(sxc,exc,syc,eyc,szc,ezc,nxc,nyc,nzc, &
                        work(ir1),sxf,exf,syf,eyf,szf,ezf,   &
                        nxf,nyf,nzf,work(ir2),comm3dp,myid,  &
                        neighbor,bd,rdatatype(1,k-1),IOUT)
          lev=1
        end if
        sxf=sxc
        exf=exc
        syf=syc
        eyf=eyc
        szf=szc
        ezf=ezc
        nxf=nxc
        nyf=nyc
        nzf=nzc
      end if
    end do
  end if
# endif
end if

!------------------------------------------------------------------------
! set phi,rhsf in work
!
sxf=sxk(ngrid)
exf=exk(ngrid)
syf=syk(ngrid)
eyf=eyk(ngrid)
szf=szk(ngrid)
ezf=ezk(ngrid)
ipf=kpbgn(ngrid)
irf=kcbgn(ngrid)+7*(exf-sxf+3)*(eyf-syf+3)*(ezf-szf+3)
call mgdsetf(sxf,exf,syf,eyf,szf,ezf,work(ipf),work(irf),phif,rhsf,IOUT)
!------------------------------------------------------------------------
! cycling at kcur=ngrid level
!
kcur=ngrid
do iter=1,maxcy
  call mgdkcyc(work,rhsf,kcur,kcycle,iprer,ipost,iresw,         &
               comm3dp,comm3dl,comm3dc,neighbor,bd,phibc,IOUT)
  sxm=sxk(ngrid)
  exm=exk(ngrid)
  sym=syk(ngrid)
  eym=eyk(ngrid)
  szm=szk(ngrid)
  ezm=ezk(ngrid)
  ip=kpbgn(ngrid)
  call mgderr(relmax,sxm,exm,sym,eym,szm,ezm,phif,work(ip),comm3d,IOUT)
  if (relmax.le.tolmax) goto 1000
end do
!------------------------------------------------------------------------
! If not converged in maxcy, issue an error message and quit
!
write(IOUT,100) maxcy,relmax
100 format('WARNING: failed to achieve convergence in ',i5,' cycles  error=',e12.5)
nerror=1
return
!------------------------------------------------------------------------
! converged
!
1000 continue
!
! rescale phif
!
if (isol.eq.1) then
  avo=rro
else
  avo=0.0d0
end if
call gscale(sx,ex,sy,ey,sz,ez,phif,avo,acorr,comm3d,nx,ny,nz,isol,IOUT)
!
! exchange boundary data and impose periodic BCs
!

ireq=0

call gxch1pla(sxm,exm,sym,eym,szm,ezm,phif,comm3dp,neighbor, &
              bd,kdatatype(1,ngrid),req,ireq,IOUT)
call gxch1lin(sxm,exm,sym,eym,szm,ezm,phif,comm3dl,neighbor, &
              bd,kdatatype(4,ngrid),req,ireq,IOUT)
call gxch1cor(sxm,exm,sym,eym,szm,ezm,phif,comm3dc,neighbor, &
              bd,kdatatype(7,ngrid),req,ireq,IOUT)

call MPI_WAITALL(ireq,req,status,ierr)

# if WMGD
!
! impose wall and Dirichlet BCs
!
call mgdbdry(sx,ex,sy,ey,sz,ez,phif,bd,vbc,IOUT)
# endif

if (isol.eq.1) then
   if (nprscr.and.myid.eq.0) write(IOUT,110) relmax,iter,acorr
   110 format('  R MGD     err=',e8.3,' iters=',i5,' rcorr=',e9.3)
else
  if (nprscr.and.myid.eq.0) write(IOUT,120) relmax,iter,acorr
   120 format('  P MGD     err=',e8.3,' iters=',i5,' pcorr=',e9.3)
end if

return
end
