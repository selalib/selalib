      subroutine mgdsolver(isol,sx,ex,sy,ey,sz,ez,phif,rhsf,r,ngrid,
     1                     work,maxcy,tolmax,kcycle,iprer,ipost,
     2                     iresw,xl,yl,zl,rro,nx,ny,nz,comm3d,
     3                     comm3dp,comm3dl,comm3dc,myid,neighbor,
     4                     bd,phibc,iter,nprscr,IOUT,nerror)
# include "compdir.inc"
      include "mpif.h"
      integer sx,ex,sy,ey,sz,ez,ngrid,IOUT
      integer maxcy,kcycle,iprer,ipost,iresw
      REALN phif(sx-1:ex+1,sy-1:ey+1,sz-1:ez+1)
      REALN rhsf(sx-1:ex+1,sy-1:ey+1,sz-1:ez+1)
      REALN r(sx-1:ex+1,sy-1:ey+1,sz-1:ez+1),rro
      REALN work(*),tolmax,xl,yl,zl,phibc(6,20)
      integer isol,comm3d,comm3dp,comm3dl,comm3dc,myid
      integer neighbor(26),bd(26),iter,nerror
      logical nprscr
c 
      integer nxk,nyk,nzk,sxk,exk,syk,eyk,szk,ezk
      integer kpbgn,kcbgn,kdatatype
      integer sxi,exi,syi,eyi,szi,ezi
      integer nxr,nyr,nzr,sxr,exr,syr,eyr,szr,ezr
      integer rdatatype
      common/mgd/nxk(20),nyk(20),nzk(20),
     1           sxk(20),exk(20),syk(20),eyk(20),szk(20),ezk(20),
     2           kpbgn(20),kcbgn(20),kdatatype(7,20),
     3           sxi(20),exi(20),syi(20),eyi(20),szi(20),ezi(20),
     4           nxr(20),nyr(20),nzr(20),sxr(20),exr(20),syr(20),
     5           eyr(20),szr(20),ezr(20),rdatatype(7,20)
c------------------------------------------------------------------------
c Parallel multigrid solver in 3-D cartesian coordinates for the 
c elliptic equation:      div(cof(grad(phif)))=rhsf
c
c isol=1 -> density
c isol=2 -> pressure
c 
c Written for periodic, wall (Neumann), and constant value (Dirichlet)
c BCs. There are two versions of the multigrid code, which are 
c separated by the compiler directive WMGD set in 'compdir.inc'. 
c The old version (WMGD=0) corresponds to the original, "traditional" 
c grid setup, and works well when all boundary conditions are 
c periodic. When one of the BCs is not periodic, must compile with 
c the new version (WMGD=1), which uses a different grid setup and 
c new restriction and correction operators (see 'mgdinit' for more 
c details). It is less accurate (or ,equivalently, slower to converge) 
c for the case of all-periodic BCs than the old version, but is better 
c able to handle wall BCs.
c
c Notes: - the values of the rhs contained in the array rhsf are
c          transferred to the work vector and the memory of the 
c          rhsf array is then used to store the residuals at the
c          different levels; therefore, it does not preserve its
c          initial values
c        - with the initialization and clean-up routines mgdinit
c          and mgdend, this multigrid code is self-standing
c
c Code      : mgd3, 3-D parallel multigrid solver
c Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
c Called in : main
c Calls     : mgdrpde, mgdpfpde, mgdphpde, mgdrsetf, mgdppde, mgdrtrsf,
c              -> discretize the pde
c             mgdsetf
c               -> set the initial guesses and the right-hand side
c             mgdkcyc, mgderr
c               -> do the actual cycling
c             gscale, gxch1pla, gxch1lin, gxch1cor,
c               -> rescale pressure and density around average values
c             gxch1pla, gxch1lin, gxch1cor, MPI_WAITALL (non-blocking
c             version)
c------------------------------------------------------------------------
      REALN avo,acorr
      integer sxf,exf,syf,eyf,szf,ezf,nxf,nyf,nzf
      integer sxc,exc,syc,eyc,szc,ezc,nxc,nyc,nzc
      integer sxm,exm,sym,eym,szm,ezm,nxm,nym,nzm
      integer ipf,irf,ipc,irc,ip,ic,kcur,lev,ir1,ir2
      integer ireq,req(52)
# if NBLOCKGR
      integer status(MPI_STATUS_SIZE,52),ierr
# endif
# if cdebug
      double precision tinitial
# if NBLOCKGR
      double precision tmpi
# endif
      tinitial=MPI_WTIME()
# endif
c------------------------------------------------------------------------
c discretize pde at all levels
c
      if (isol.eq.1) then
c
c density: only have to set geometric factors
c
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
          call mgdrpde(sxm,exm,sym,eym,szm,ezm,nxm,nym,nzm,work(ic),
     1                 xl,yl,zl,IOUT)
        end do
      else
c pressure: have to do a lot more work. The coefficients in the new
c           and old versions of the multigrid code are located at
c           different places on the grid.
c Note: the code requires that corasifying takes place in all 
c directions between levels ngrid and ngrid-1
c
c determine coefficients at finest grid level (mid-point values) from
c two neighboring density points
c
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
        call mgdpfpde(sxf,exf,syf,eyf,szf,ezf,nxf,nyf,nzf,work(icf),r,
     1                xl,yl,zl,IOUT)
# if WMGD
c
c new version: determine coefficients at coarser grid levels from
c four neighboring density points; no communication of boundary data
c is involved because of the grid setup, under the condition that
c mod(nx,nxprocs)=mod(ny,nyprocs)=0
c
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
          call mgdphpde(sxm,exm,sym,eym,szm,ezm,nxm,nym,nzm,work(ic),
     1                  sx,ex,sy,ey,sz,ez,nxf,nyf,nzf,r,bd,xl,yl,zl,
     2                  IOUT)
        end do
# else
        if (ngrid.gt.1) then
c
c old version: use two locations ir1 and ir2 in the work vector to
c store the density on the different grid levels. First set r in 
c work(ir1) at the finest grid level; this is used to determine the 
c coefficients at level k=ngrid-1
c
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
c
c for the levels k=ngrid-1,1, calculate the coefficients from the
c densities stored in the arrays work(ir1) and work(ir2)
c for the levels k=ngrid-1,2, transfer to the level below the values
c of the density needed there to determine the coefficients
c
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
              call mgdppde(sxm,exm,sym,eym,szm,ezm,nxm,nym,nzm,work(ic),
     1                     sxf,exf,syf,eyf,szf,ezf,work(ir1),
     2                     xl,yl,zl,IOUT)
            else
              call mgdppde(sxm,exm,sym,eym,szm,ezm,nxm,nym,nzm,work(ic),
     1                     sxf,exf,syf,eyf,szf,ezf,work(ir2),
     2                     xl,yl,zl,IOUT)
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
                call mgdrtrsf(sxc,exc,syc,eyc,szc,ezc,nxc,nyc,nzc,
     1                        work(ir2),sxf,exf,syf,eyf,szf,ezf,
     2                        nxf,nyf,nzf,work(ir1),comm3dp,myid,
     3                        neighbor,bd,rdatatype(1,k-1),IOUT)
                lev=2
              else
                call mgdrtrsf(sxc,exc,syc,eyc,szc,ezc,nxc,nyc,nzc,
     1                        work(ir1),sxf,exf,syf,eyf,szf,ezf,
     2                        nxf,nyf,nzf,work(ir2),comm3dp,myid,
     3                        neighbor,bd,rdatatype(1,k-1),IOUT)
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
c------------------------------------------------------------------------
c set phi,rhsf in work
c
      sxf=sxk(ngrid)
      exf=exk(ngrid)
      syf=syk(ngrid)
      eyf=eyk(ngrid)
      szf=szk(ngrid)
      ezf=ezk(ngrid)
      ipf=kpbgn(ngrid)
      irf=kcbgn(ngrid)+7*(exf-sxf+3)*(eyf-syf+3)*(ezf-szf+3)
      call mgdsetf(sxf,exf,syf,eyf,szf,ezf,work(ipf),work(irf),
     1             phif,rhsf,IOUT)
c------------------------------------------------------------------------
c cycling at kcur=ngrid level
c
      kcur=ngrid
      do iter=1,maxcy
        call mgdkcyc(work,rhsf,kcur,kcycle,iprer,ipost,iresw,
     1               comm3dp,comm3dl,comm3dc,neighbor,bd,
     2               phibc,IOUT)
        sxm=sxk(ngrid)
        exm=exk(ngrid)
        sym=syk(ngrid)
        eym=eyk(ngrid)
        szm=szk(ngrid)
        ezm=ezk(ngrid)
        ip=kpbgn(ngrid)
        call mgderr(relmax,sxm,exm,sym,eym,szm,ezm,phif,work(ip),
     1              comm3d,IOUT)
        if (relmax.le.tolmax) goto 1000
      end do
c------------------------------------------------------------------------
c If not converged in maxcy, issue an error message and quit
c
      write(IOUT,100) maxcy,relmax
100   format('WARNING: failed to achieve convergence in ',i5,
     1       ' cycles  error=',e12.5)
      nerror=1
      return
c------------------------------------------------------------------------
c converged
c
1000  continue
c
c rescale phif
c
      if (isol.eq.1) then
        avo=rro
      else
        avo=0.0d0
      end if
      call gscale(sx,ex,sy,ey,sz,ez,phif,avo,acorr,comm3d,nx,ny,nz,
     1            isol,IOUT)
c
c exchange boundary data and impose periodic BCs
c
# if NBLOCKGR
      ireq=0
# endif
      call gxch1pla(sxm,exm,sym,eym,szm,ezm,phif,comm3dp,neighbor,
     1              bd,kdatatype(1,ngrid),req,ireq,IOUT)
      call gxch1lin(sxm,exm,sym,eym,szm,ezm,phif,comm3dl,neighbor,
     1              bd,kdatatype(4,ngrid),req,ireq,IOUT)
      call gxch1cor(sxm,exm,sym,eym,szm,ezm,phif,comm3dc,neighbor,
     1              bd,kdatatype(7,ngrid),req,ireq,IOUT)
# if NBLOCKGR
# if cdebug
      tmpi=MPI_WTIME()
# endif
      call MPI_WAITALL(ireq,req,status,ierr)
# if cdebug
      nwaitall=nwaitall+1
      twaitall=twaitall+MPI_WTIME()-tmpi
# endif
# endif
# if WMGD
c
c impose wall and Dirichlet BCs
c
      call mgdbdry(sx,ex,sy,ey,sz,ez,phif,bd,vbc,IOUT)
# endif
c
      if (isol.eq.1) then
        if (nprscr.and.myid.eq.0) write(IOUT,110) relmax,iter,acorr
110     format('  R MGD     err=',e8.3,' iters=',i5,' rcorr=',e9.3)
# if cdebug
        timing(95)=timing(95)+MPI_WTIME()-tinitial
# endif
      else
        if (nprscr.and.myid.eq.0) write(IOUT,120) relmax,iter,acorr
120     format('  P MGD     err=',e8.3,' iters=',i5,' pcorr=',e9.3)
# if cdebug
        timing(96)=timing(96)+MPI_WTIME()-tinitial
# endif
      end if
c
      return
      end
