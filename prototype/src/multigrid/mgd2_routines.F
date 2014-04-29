
      subroutine mgdbdry(sxm,exm,sym,eym,phi,bd,phibc,IOUT)
#include "compdir.inc"
      include "mpif.h"
      integer sxm,exm,sym,eym,bd(8),IOUT
      REALN phi(sxm-1:exm+1,sym-1:eym+1),phibc(4)
c------------------------------------------------------------------------
c Enforce the Neumann and Dirichlet boundary conditions
c
c Code      : mgd2, 2-D parallel multigrid solver
c Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
c Called in : mgdrelax, mgdsolver
c Calls     : --
c------------------------------------------------------------------------
      integer i,j
# if cdebug
      double precision tinitial
      tinitial=MPI_WTIME()
# endif
c
      if (bd(1).eq.1) then
        do j=sym-1,eym+1
          phi(exm+1,j)= 0.0 !phi(exm,j)
        end do
      else if (bd(1).eq.2) then
        do j=sym-1,eym+1
          phi(exm+1,j)=2.0d0*phibc(1)-phi(exm,j)
        end do
      end if
      if (bd(5).eq.1) then
        do j=sym-1,eym+1
          phi(sxm-1,j)= 0.0 !phi(sxm,j)
        end do
      else if (bd(5).eq.2) then
        do j=sym-1,eym+1
          phi(sxm-1,j)=2.0d0*phibc(3)-phi(sxm,j)
        end do
      end if
      if (bd(3).eq.1) then
        do i=sxm-1,exm+1
          phi(i,sym-1)= 0.0 !phi(i,sym)
        end do
      else if (bd(3).eq.2) then
        do i=sxm-1,exm+1
          phi(i,sym-1)=2.0d0*phibc(2)-phi(i,sym)
        end do
      end if
      if (bd(7).eq.1) then
        do i=sxm-1,exm+1
          phi(i,eym+1)= 0.0 !phi(i,eym)
        end do
      else if (bd(7).eq.2) then
        do i=sxm-1,exm+1
          phi(i,eym+1)=2.0d0*phibc(4)-phi(i,eym)
        end do
      end if
c
# if cdebug
      timing(93)=timing(93)+MPI_WTIME()-tinitial
# endif
      return
      end

      subroutine mgdcor(sxf,exf,syf,eyf,nxf,nyf,phif,
     1                  sxc,exc,syc,eyc,nxc,nyc,phic,
     2                  sx1,ex1,sy1,ey1,bd,phibc,IOUT)
#include "compdir.inc"
      include "mpif.h"
      integer sxf,exf,syf,eyf,nxf,nyf,IOUT
      integer sxc,exc,syc,eyc,nxc,nyc,sx1,ex1,sy1,ey1,bd(8)
      REALN phif(sxf-1:exf+1,syf-1:eyf+1)
      REALN phic(sxc-1:exc+1,syc-1:eyc+1),phibc(4)
c------------------------------------------------------------------------
c Add correction from coarse grid level to fine grid level. Uses
c bilinear interpolation for the old version of the multigrid code,
c and area weighting for its new version.
c
c Tested for the case where coarsifying takes place in all directions
c
c Code      : mgd2, 2-D parallel multigrid solver
c Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
c Called in : mgdkcyc
c Calls     : mgdbdry
c------------------------------------------------------------------------
      integer i,j,ic,jc,i1,i2,j1,j2
# if cdebug
      double precision tinitial
      tinitial=MPI_WTIME()
# endif
c
# if WMGD
c------------------------------------------------------------------------
c new version: the correction is the weighted average of either two 
c or four points at the coarser grid level depending on whether 
c coarsifying takes place in all directions or not
c
      if (nxf.eq.nxc) then
        do jc=syc-1,eyc
          j=2*jc-1
          do i=sxf-1,exf+1
            ic=i
            phif(i,j)=phif(i,j)+
     1        (3.0d0*phic(ic,jc)+phic(ic,jc+1))/4.0d0
            phif(i,j+1)=phif(i,j+1)+
     1        (phic(ic,jc)+3.0d0*phic(ic,jc+1))/4.0d0
          end do
        end do
      else if (nyf.eq.nyc) then
        do j=syf-1,eyf+1
          jc=j
          do ic=sxc-1,exc
            i=2*ic-1
            phif(i,j)=phif(i,j)+
     1        (3.0d0*phic(ic,jc)+phic(ic+1,jc))/4.0d0
            phif(i+1,j)=phif(i+1,j)+
     1        (phic(ic,jc)+3.0d0*phic(ic+1,jc))/4.0d0
          end do
        end do
      else
        do jc=syc-1,eyc
          j=2*jc-1
          do ic=sxc-1,exc
            i=2*ic-1
            phif(i,j)=phif(i,j)+
     1        (9.0d0*phic(ic,jc)+3.0d0*phic(ic+1,jc)
     2        +3.0d0*phic(ic,jc+1)+phic(ic+1,jc+1))/16.0d0
            phif(i+1,j)=phif(i+1,j)+
     1        (3.0d0*phic(ic,jc)+9.0d0*phic(ic+1,jc)
     2        +phic(ic,jc+1)+3.0d0*phic(ic+1,jc+1))/16.0d0
            phif(i,j+1)=phif(i,j+1)+
     1        (3.0d0*phic(ic,jc)+phic(ic+1,jc)
     2        +9.0d0*phic(ic,jc+1)+3.0d0*phic(ic+1,jc+1))/16.0d0
            phif(i+1,j+1)=phif(i+1,j+1)+
     1        (phic(ic,jc)+3.0d0*phic(ic+1,jc)
     2        +3.0d0*phic(ic,jc+1)+9.0d0*phic(ic+1,jc+1))/16.0d0
          end do
        end do
      end if
c
c impose Neumann and Dirichlet boundary conditions
C TEMP: periodicity is not enforced to save one call to gxch1lin;
c check whether it has an impact or not...
c
      call mgdbdry(sxf,exf,syf,eyf,phif,bd,phibc,IOUT)
# else
c------------------------------------------------------------------------
c old version
c
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
c
c identity at the points of the fine grid which have odd indices 
c in i and j
c
      do jc=sy1,ey1
        j=j1*(2*jc-1)+j2*jc
        do ic=sx1,ex1
          i=i1*(2*ic-1)+i2*ic
          phif(i,j)=phif(i,j)+phic(ic,jc)
        end do
      end do
c
c interpolation of the two neighboring values for the points of the 
c fine grid with even index for i and odd index for j
c
      if (nxc.lt.nxf) then
        do jc=sy1,ey1
          j=j1*(2*jc-1)+j2*jc
          do ic=sxc-1,exc
            i=2*ic
            phif(i,j)=phif(i,j)+0.5d0*(phic(ic,jc)+phic(ic+1,jc))
          end do
        end do
      end if
c
c interpolation of the two neighboring values for the points of the 
c fine grid with odd index for i and even index for j
c
      if (nyc.lt.nyf) then
        do jc=syc-1,eyc
          j=2*jc
          do ic=sx1,ex1
            i=i1*(2*ic-1)+i2*ic
            phif(i,j)=phif(i,j)+0.5d0*(phic(ic,jc)+phic(ic,jc+1))
          end do
        end do
      end if
c
c interpolation of the four neighboring values for the points of the
c fine grid with even indices for i and j
c
      if (nxc.lt.nxf.and.nyc.lt.nyf) then
        do jc=syc-1,eyc
          j=2*jc
          do ic=sxc-1,exc
            i=2*ic
            phif(i,j)=phif(i,j)+0.25d0*(phic(ic,jc)+phic(ic+1,jc)
     1                                 +phic(ic,jc+1)+phic(ic+1,jc+1))
          end do
        end do
      end if
# endif
c
# if cdebug
      timing(92)=timing(92)+MPI_WTIME()-tinitial
# endif
      return
      end

      subroutine mgdend(ngrid)
#include "compdir.inc"
      include "mpif.h"
      integer ngrid
c
c common for multigrid indices and datatypes
c
      integer nxk,nyk,sxk,exk,syk,eyk,kpbgn,kcbgn
      integer ikdatatype,jkdatatype,ijkdatatype
      integer sxi,exi,syi,eyi
      integer nxr,nyr,sxr,exr,syr,eyr
      integer irdatatype,jrdatatype,ijrdatatype
      integer ierr
      common/mgd/nxk(20),nyk(20),sxk(20),exk(20),syk(20),eyk(20),
     1           kpbgn(20),kcbgn(20),ikdatatype(20),jkdatatype(20),
     2           ijkdatatype(20),sxi(20),exi(20),syi(20),eyi(20),
     3           nxr(20),nyr(20),sxr(20),exr(20),syr(20),eyr(20),
     4           irdatatype(20),jrdatatype(20),ijrdatatype(20)

c------------------------------------------------------------------------
c Free the MPI datatypes associated witht the multigrid code
c
c Code      : mgd2, 2-D parallel multigrid solver
c Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
c Called in : main
c Calls     : MPI_TYPE_FREE
c------------------------------------------------------------------------
      integer k
# if cdebug
      double precision tinitial
      tinitial=MPI_WTIME()
# endif
c
      do k=1,ngrid-1
        call MPI_TYPE_FREE(ikdatatype(k),ierr)
        call MPI_TYPE_FREE(jkdatatype(k),ierr)
        call MPI_TYPE_FREE(ijkdatatype(k),ierr)
# if !WMGD
        call MPI_TYPE_FREE(irdatatype(k),ierr)
        call MPI_TYPE_FREE(jrdatatype(k),ierr)
        call MPI_TYPE_FREE(ijrdatatype(k),ierr)
# endif
      end do
c
# if cdebug
      timing(97)=timing(97)+MPI_WTIME()-tinitial
# endif
      return
      end

      subroutine mgderr(relmax,sxm,exm,sym,eym,phio,phin,comm2d,IOUT)
#include "compdir.inc"
      include "mpif.h"
      integer sxm,exm,sym,eym,comm2d,IOUT
      integer ierr
      REALN relmax
      REALN phio(sxm-1:exm+1,sym-1:eym+1)
      REALN phin(sxm-1:exm+1,sym-1:eym+1)
c------------------------------------------------------------------------
c Calculate the error between the new and old iterates of phi and 
c save the new iterate into the phio array.
c
c Code      : mgd2, 2-D parallel multigrid solver
c Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
c Called in : mgdsolver
c Calls     : MPI_ALLREDUCE
c------------------------------------------------------------------------
      REALN phloc,reloc
      integer i,j
# if cdebug
      double precision tinitial
      tinitial=MPI_WTIME()
# endif
c
c calculate local error
c
      phloc=0.0d0
      reloc=0.0d0
      do j=sym,eym
        do i=sxm,exm
          phloc=max(phloc,abs(phin(i,j)))
          reloc=max(reloc,abs(phin(i,j)-phio(i,j)))
        end do
      end do
      if (phloc.gt.0.0) then
        reloc=reloc/phloc
      else
        reloc=0.0d0
      end if
c
c global reduce across all processes
c
      call MPI_ALLREDUCE(reloc,relmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,
     1                   comm2d,ierr)
# if cdebug
      nallreduce=nallreduce+1
# endif
c
c save new values into ouput array
c
      do j=sym-1,eym+1
        do i=sxm-1,exm+1
          phio(i,j)=phin(i,j)
        end do
      end do
c
# if cdebug
      timing(94)=timing(94)+MPI_WTIME()-tinitial
# endif
      return
      end


      subroutine mgdinit(vbc,phibc,ixp,jyq,iex,jey,ngrid,nxp2,nyp2,
     1                   sx,ex,sy,ey,realtype,nxprocs,nyprocs,nwork,
     3                   ibdry,jbdry,myid,IOUT,nerror)
#include "compdir.inc"
      include "mpif.h"
      integer ixp,jyq,iex,jey,ngrid,nxp2,nyp2,sx,ex,sy,ey
      integer realtype,nxprocs,nyprocs,nwork,ibdry,jbdry
      integer myid,IOUT,nerror
      REALN vbc(4),phibc(4,20)
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
c Initialize the parallel multigrid solver: subdomain indices and
c MPI datatypes.
c
c The multigrid code comes in two versions. With the WMGD compiler
c directive set to 0, the grid setup is vertex-centered:
c
c WMGD=0
c 
c  |------|-----|-----|-----|-----|            fine
c  1      2     3     4     5     6 
c
c  |------------|-----------|-----------|      coarse
c  1            2           3           4
c
c With WMGD set to 1, it is cell-centered:
c
c WMGD=1   
c           |                       |
c        |--|--|-----|-----|-----|--|--|       fine
c        1  |  2     3     4     5  |  6
c           |                       |
c     |-----|-----|-----------|-----|-----|    coarse
c     1     |     2           3     |     4
c           |                       |
c          wall                    wall
c
c For WMGD=0, the restriction and correction operators are standard
c (choice of full or half weighting for the restriction, bilinear
c interpolation for the correction). This works fine for periodic
c boundary conditions. However, when there are Neumann (wall) or
c Dirichlet BCs, this grid setup results in a loss of accuracy near
c the boundaries when the grid is staggered (the discretization of
c the relaxation operator is first-order locally there). With the
c grid setup corresponding to WMGD=1, accuracy remains second-order
c all the time. As the grid gets coarser, it remains centered on the
c domain instead of "shifting to the right". This option works for
c periodic, Neumann, and Dirichlet BCs, although only periodic and
c Neumann BCs have been tested thoroughly. There is one catch, though.
c For a problem with purely periodic BCs, WMGD=0 converges in less
c cycles than WMGD=1 and requires less CPU time (the penalty is
c apparently between 10 and 50%). This can be attributed to the loss
c of accuracy in the restriction and correction operators due to the
c fact that WMGD=0 uses a support of 3 points in each direction 
c whereas WMGD=1 uses only 2 points.
c
c Both versions offer the option to coarsify in one direction and
c not the other, except at the finest grid level, where coarsifying
c MUST take place along all axes. However, it is possible to have
c ngrid=iex=jey=1, with ixp=nx and jyq=ny. In this case, the code
c never enters 'mgdrestr' and 'mgdcor' and all it does is Gauss-Seidel
c iterate at the finest grid level. This can be useful as a preliminary
c check.
c
c Note: some memory could be saved by noting that the cof arrays
c need be dimensioned (sxm:exm,sym:eym) and not
c (sxm-1:exm+1,sym-1:eym+1)... Probably not too difficult to
c make the change
c
c Code      : mgd2, 2-D parallel multigrid solver
c Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
c Called in : main
c Calls     : grid1_type
c------------------------------------------------------------------------

      integer i,j,k,nxf,nyf,nxm,nym,kps,sxm,exm,sym,eym,nxc,nyc
# if cdebug
      double precision tinitial
      tinitial=MPI_WTIME()
# endif
c------------------------------------------------------------------------
c set /mgd/ variables to zero
c
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
c------------------------------------------------------------------------
c make a number of checks
c
# if WMGD
c
c check that, for the new version, the number of processes in each
c direction divides the number of points in that direction (the
c determination of the pressure coefficients in 'mgdphpde' and the
c restriction in 'mgdrestr' would not be complete and would require
c doing some inter-process data communication - warning: would
c be very complex because of the pressure coefficients)
c
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
100   format(/,'ERROR in mgdinit: nx=',i3,' is not a multiple of ',
     1       'nxprocs=',i3,/,'cannot use the new version of the ',
     2       'multigrid code',/)
110   format(/,'ERROR in mgdinit: ny=',i3,' is not a multiple of ',
     1       'nyprocs=',i3,/,'cannot use the new version of the ',
     2       'multigrid code',/)
# else
c
c check that the old version is not used with non-periodic BCs
c
      if (ibdry.ne.0.or.jbdry.ne.0) then
        write(IOUT,120) ibdry,jbdry
        nerror=1
        return
      end if
120   format(/,'ERROR in mgdinit: ibdry=',i2,' jbdry=',i2,
     1       /,'cannot use the old version of the multigrid code',
     2       /,'boundary conditions that are not periodic',
     3       /,'-> change compiler directive to 1 in compdir.inc',
     4       /,'   and recompile the multigrid code',/)
# endif
c
c check that the dimensions are correct
c
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
130   format(/,'ERROR in mgdinit: nxp1=',i3,' <> ixp*2**(iex-1)+1=',
     1       i3,/,'-> adjust the multigrid parameters ixp and iex',
     2       ' in main',/)
140   format(/,'ERROR in mgdinit: nyp1=',i3,' <> jyq*2**(jey-1)+1=',
     1       i3,/,'-> adjust the multigrid parameters jyq and jey',
     2       ' in main',/)
c
c check that the number of points at the coarser level is not smaller
c than the number of processes in either direction
c
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
150   format(/,'ERROR in mgdinit: ixp=',i3,' < nxprocs=',i3,/,
     1       ' there must be at least one grid point at the ',
     2       'coarsest grid level',/,
     3       '-> increase ixp and decrease iex correspondingly',
     4       ' in main',/)
160   format(/,'ERROR in mgdinit: jyq=',i3,' < nyprocs=',i3,/,
     1       ' there must be at least one grid point at the ',
     2       'coarsest grid level',/,
     3       '-> increase jyq and decrease jey correspondingly',
     4       ' in main',/)
c
c check that coarsifying takes place in all directions at the finest
c grid level
c
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
170   format(/,'ERROR in mgdinit: ngrid=',i3,' iex=',i3,
     1       /,'no coarsifying at the finest grid level in x-direction',
     2       /,'this is not allowed by the mutligrid code',/)
180   format(/,'ERROR in mgdinit: ngrid=',i3,' jey=',i3,
     1       /,'no coarsifying at the finest grid level in y-direction',
     2       /,'this is not allowed by the mutligrid code',/)
c------------------------------------------------------------------------
c define all grid levels
c I have adopted the same notations as in Mudpack as far as possible.
c When a confusion was possible, I added a suffix 'm' to the name
c of the variables. For example, nxm is nxp2-1 for the multigrid
c code whereas nx means nxp2-2 in the rest of the code.
c
      do k=1,ngrid
        nxk(k)=ixp*2**(max(k+iex-ngrid,1)-1)+1
        nyk(k)=jyq*2**(max(k+jey-ngrid,1)-1)+1
      end do
c
c for all grid levels, set the indices of the subdomain the process
c 'myid' work on, as well as the datatypes needed for the exchange
c of boundary data
c
      nxf=nxk(ngrid)
      nyf=nyk(ngrid)
      sxk(ngrid)=sx
      exk(ngrid)=ex
      syk(ngrid)=sy
      eyk(ngrid)=ey
      call grid1_type(ikdatatype(ngrid),jkdatatype(ngrid),
     1                ijkdatatype(ngrid),realtype,sx,ex,sy,ey,IOUT)
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
        call grid1_type(ikdatatype(k),jkdatatype(k),ijkdatatype(k),
     1                  realtype,sxk(k),exk(k),syk(k),eyk(k),IOUT)
      end do
# if xdebug1
c
c print out the indices and determine the size of the MPI messages
c as a rough check
c 
      write(IOUT,*) 'size of the multigrid phi-messages'
      do k=ngrid,1,-1
        call MPI_TYPE_SIZE(ikdatatype(k),nsiz1,ierr)
        call MPI_TYPE_SIZE(jkdatatype(k),nsiz2,ierr)
        call MPI_TYPE_SIZE(ijkdatatype(k),nsiz3,ierr)
        write(IOUT,*) 'myid: ',myid,' k=',k,' sxk=',sxk(k),' exk=',
     1                exk(k),' syk=',syk(k),' eyk=',eyk(k),
     2                ' size of datatypes: ',nsiz1,nsiz2,nsiz3
      end do
# endif
c
c set work space indices for phi, cof at each grid level
c check that there is sufficient work space
c
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
200     format(/,'ERROR in mgdinit: not enough work space',/,
     1         ' kps=',i10,' nwork=',i10,' myid: ',i3,/,
     2         '-> put the formula for nwork in main in ',
     3         'comments',/,'   and set nwork to the value of kps',/)
      else
!        write(IOUT,210) kps,nwork
!210     format(/,'WARNING in mgdinit: kps=',i10,' nwork=',i10,
!     1         /,'can optimize the amount of memory needed by ',
!     2           'the multigrid code',/,'by putting the formula ',
!     3           'for nwork into comments and setting',/,'nwork ',
!     4           'to the value of kps',/)
         nwork = kps

      end if
# if WMGD
c------------------------------------------------------------------------
c For the new version of the multigrid code, set the boundary values 
c to be used for the Dirichlet boundaries. It is possible to assign
c 4 different constant values to the 4 different sides. The values are
c assigned at the finest grid level, zero is assigned at all levels
c below
c 
c vbc, phibc:
c
c       -----vbc(4)------ 
c       |               |
c       |               |
c     vbc(3)----------vbc(1)
c       |               |
c       |               |
c       ------vbc(2)-----
c
      do j=1,4
        phibc(j,ngrid)=vbc(j)
      end do
      do k=ngrid-1,1,-1
        do j=1,4
          phibc(j,k)=0.0d0
        end do
      end do
# else
c------------------------------------------------------------------------
c set indices for range of ic and jc on coarser grid level which are
c supported on finer grid level, i.e. for which the points
c (x(2*ic-1,2*jc-1),y(2*ic-1,2*jc-1)) are defined in the subdomain
c of process 'myid'. This allows to avoid having any communication
c after the interpolation (or prolongation) step; it should be used
c in that operation only. 
c
c example:
c   a)  - - -|-----------|-----------|-----------|
c                                  exf=8      exf+1=9
c
c       - ---------------|-----------------------|
c                      exc=4                  exc+1=5
c
c   b)  - - -|-----------|
c          exf=5      exf+1=6
c
c       - - -|-----------------------|
c          exc=3                  exc+1=4
c
      
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
        write(IOUT,*) 'myid: ',myid,' k=',k,' sxi=',sxi(k),' exi=',
     1                exi(k),' syi=',syi(k),' eyi=',eyi(k)
# endif
        nxf=nxc
        nyf=nyc
      end do
c------------------------------------------------------------------------
c set indices for determining the coefficients in the elliptic
c equation div(cof*grad(P))=rhs. Used only when solving for the
c pressure. When setting these coefficients at level k, need
c the values of the density at midpoints, i.e. at level k+1
c (if coarsifying takes place between the levels k and k+1).
c If coarsifying took place at all levels, the array cof could
c be used as temporary storage space for the densities, with
c cof(*,*,6) at level k+1 giving the values cof(*,*,1->5) at level
c k, and the already defined datatypes could be used for the 
c exchange of the boundary values. However, this does not work
c in case there is coarsifying in one direction between two levels.
c
c Example: - - -|----------|----------|----------|- -  
c               3          4          5          6    \
c                                                      | coarsifying
c                          |                     |     |
c                          |                     |    /
c                          V                     V
c          - - -|---------------------|---------------------|- -
c               2          r          3          r   \      4
c                                                     |
c                          |                     |    | no coarsifying
c                          |                     |    /
c                          V                     V
c          - - -|---------------------|---------------------|- -
c               2          r          3          r          4
c
c At the finest grid level, the coefficients are determined by a 
c special procedure directly from r, so that no value needs to be
c assigned to the indices at level ngrid.
c
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
300         format('ERROR in mgdinit: no coarsifying between level ',
     1             i3,' and level ',i3,/,', the current version of ',
     2             'the code cannot cope with that',/,
     3             ' -> decrease the value of ixp and/or jyq and',/,
     4             '    increase the value of iex and/or jey')
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
        call grid1_type(irdatatype(k),jrdatatype(k),ijrdatatype(k),
     1                  realtype,sxr(k),exr(k),syr(k),eyr(k),IOUT)
      end do
# if xdebug1
c
c print out the indices and determine the size of the MPI messages
c as a rough check
c 
      write(IOUT,*) 'size of the r-messages'
      do k=ngrid-1,1,-1
        call MPI_TYPE_SIZE(irdatatype(k),nsiz1,ierr)
        call MPI_TYPE_SIZE(jrdatatype(k),nsiz2,ierr)
        call MPI_TYPE_SIZE(ijrdatatype(k),nsiz3,ierr)
        write(IOUT,*) 'myid: ',myid,' k=',k,' sxr=',sxr(k),' exr=',
     1                exr(k),' syr=',syr(k),' eyr=',eyr(k),
     2                ' size of datatypes: ',nsiz1,nsiz2,nsiz3
      end do
# endif
# endif
c
# if cdebug
      timing(81)=timing(81)+MPI_WTIME()-tinitial
# endif
      return
      end

      subroutine mgdkcyc(work,res,kcur,kcycle,iprer,ipost,iresw,
     1                   comm2d,myid,neighbor,bd,phibc,IOUT)
#include "compdir.inc"
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

      subroutine mgdpfpde(sxf,exf,syf,eyf,nxf,nyf,cof,r,xl,yl,bd,IOUT)
#include "compdir.inc"
      include "mpif.h"
      integer sxf,exf,syf,eyf,nxf,nyf,bd(8),IOUT
      REALN cof(sxf-1:exf+1,syf-1:eyf+1,6)
      REALN r(sxf-1:exf+1,syf-1:eyf+1),xl,yl
c------------------------------------------------------------------------
c Determine coefficients for the pressure equation at the finest grid
c level. These coefficients involve densities half-way between the
c pressure and density nodes. Works for periodic, Neumann, and
c Dirichlet boundary conditions.
c
c cof array:
c
c         cof(4)
c           |
c           |
c cof(1)--cof(5)--cof(2)
c           |
c           |
c         cof(3)
c
c Code      : mgd2, 2-D parallel multigrid solver
c Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
c Called in : mgdsolver
c Calls     : --
c------------------------------------------------------------------------
      REALN dlx,todlxx,dly,todlyy,rij
      integer i,j
# if cdebug
      double precision tinitial
      tinitial=MPI_WTIME()
# endif
c
      dlx=xl/float(nxf-1)
      todlxx=2.0d0/(dlx*dlx)
      dly=yl/float(nyf-1)
      todlyy=2.0d0/(dly*dly)
      do j=syf,eyf
        do i=sxf,exf
          rij=r(i,j)
          cof(i,j,1)=todlxx/(rij+r(i-1,j))
          cof(i,j,2)=todlxx/(rij+r(i+1,j))
          cof(i,j,3)=todlyy/(rij+r(i,j-1))
          cof(i,j,4)=todlyy/(rij+r(i,j+1))
        end do
      end do
c
c enforce wall BCs 
c
      if (bd(1).eq.1) then
        do j=syf,eyf
          cof(exf,j,2)=0.0d0
        end do
      end if
      if (bd(5).eq.1) then
        do j=syf,eyf
          cof(sxf,j,1)=0.0d0
        end do
      end if
      if (bd(3).eq.1) then
        do i=sxf,exf
          cof(i,syf,3)=0.0d0
        end do
      end if
      if (bd(7).eq.1) then
        do i=sxf,exf
          cof(i,eyf,4)=0.0d0
        end do
      end if
c
c calculate diagonal term
c
      do j=syf,eyf
        do i=sxf,exf
          cof(i,j,5)=-(cof(i,j,1)+cof(i,j,2)+cof(i,j,3)+cof(i,j,4))
        end do
      end do
c
# if cdebug
      timing(84)=timing(84)+MPI_WTIME()-tinitial
# endif
      return
      end

      subroutine mgdphpde(sxm,exm,sym,eym,nxm,nym,cof,
     1                    sx,ex,sy,ey,nxf,nyf,r,bd,xl,yl,IOUT)
#include "compdir.inc"
      include "mpif.h"
      integer sxm,exm,sym,eym,nxm,nym,sx,ex,sy,ey,nxf,nyf,bd(8),IOUT
      REALN cof(sxm-1:exm+1,sym-1:eym+1,6)
      REALN r(sx-1:ex+1,sy-1:ey+1),xl,yl
c------------------------------------------------------------------------
c For the new version of the multigrid code, determine the coefficients
c for the pressure equation at all grid levels except the finest one.
c The coefficients are determined directly from the density array r
c through some manipulation of indices and are values at (i+1/2,j+1/2)
c points. Works for periodic, Neumann, and Dirichlet boundary
c conditions.
c
c cof array:
c
c         cof(4)
c           |
c           |
c cof(1)--cof(5)--cof(2)
c           |
c           |
c         cof(3)
c
c Code      : mgd2, 2-D parallel multigrid solver
c Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
c Called in : mgdsolver
c Calls     : --
c------------------------------------------------------------------------
      REALN dlx,fodlxx,dly,fodlyy
      integer i,j,im,jm,is,js,istep,jstep
# if cdebug
      double precision tinitial
      tinitial=MPI_WTIME()
# endif
c
c calculate off-diagonal terms
c
      dlx=xl/float(nxm-1)
      fodlxx=4.0d0/(dlx*dlx)
      dly=yl/float(nym-1)
      fodlyy=4.0d0/(dly*dly)
      istep=(nxf-1)/(nxm-1)
      jstep=(nyf-1)/(nym-1)
      do j=sym,eym
        jm=2*jstep*j-3*(jstep-1)
        do i=sxm,exm
          im=2*istep*i-3*(istep-1)
          is=(im-istep)/2
          js=jm/2
          cof(i,j,1)=fodlxx/(r(is,js)+r(is+1,js)
     1                      +r(is,js+1)+r(is+1,js+1))
          is=(im+istep)/2
          cof(i,j,2)=fodlxx/(r(is,js)+r(is+1,js)
     1                      +r(is,js+1)+r(is+1,js+1))
          is=im/2
          js=(jm-jstep)/2
          cof(i,j,3)=fodlyy/(r(is,js)+r(is+1,js)
     1                      +r(is,js+1)+r(is+1,js+1))
          js=(jm+jstep)/2
          cof(i,j,4)=fodlyy/(r(is,js)+r(is+1,js)
     1                      +r(is,js+1)+r(is+1,js+1))
        end do
      end do
c
c enforce wall BCs
c
      if (bd(1).eq.1) then
        do j=sym,eym
          cof(exm,j,2)=0.0d0
        end do
      end if
      if (bd(5).eq.1) then
        do j=sym,eym
          cof(sxm,j,1)=0.0d0
        end do
      end if
      if (bd(3).eq.1) then
        do i=sxm,exm
          cof(i,sym,3)=0.0d0
        end do
      end if
      if (bd(7).eq.1) then
        do i=sxm,exm
          cof(i,eym,4)=0.0d0
        end do
      end if
c
c calculate diagonal term
c
      do j=sym,eym
        do i=sxm,exm
          cof(i,j,5)=-(cof(i,j,1)+cof(i,j,2)+cof(i,j,3)+cof(i,j,4))
        end do
      end do
c
# if cdebug
      timing(83)=timing(83)+MPI_WTIME()-tinitial
# endif
      return
      end

      subroutine mgdppde(sxm,exm,sym,eym,nxm,nym,cof,
     1                   sxf,exf,syf,eyf,rf,xl,yl,bd,IOUT)
#include "compdir.inc"
      include "mpif.h"
      integer sxm,exm,sym,eym,nxm,nym,sxf,exf,syf,eyf,bd(8),IOUT
      REALN cof(sxm-1:exm+1,sym-1:eym+1,6)
      REALN rf(sxf-1:exf+1,syf-1:eyf+1),xl,yl
c------------------------------------------------------------------------
c For the old version of the multigrid code, determine coefficients 
c for the pressure equation at all grid levels but the finest one.
c The coefficients are determined from the values of the density
c at integer nodes (i,j). Works only for periodic boundary conditions.
c
c cof array:
c
c         cof(4)
c           |
c           |
c cof(1)--cof(5)--cof(2)
c           |
c           |
c         cof(3)
c
c Code      : mgd2, 2-D parallel multigrid solver
c Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
c Called in : mgdsolver
c Calls     : --
c------------------------------------------------------------------------
      REALN dlx,odlxx,dly,odlyy
      integer i,j,is,js
      double precision c1, c2, c3, c4
# if cdebug
      double precision tinitial
      tinitial=MPI_WTIME()
# endif
c
c calculate off-diagonal terms
c
      dlx=xl/float(nxm-1)
      odlxx=1.0d0/(dlx*dlx)
      dly=yl/float(nym-1)
      odlyy=1.0d0/(dly*dly)
      do j=sym,eym
        js=2*j-1
        do i=sxm,exm
          is=2*i-1
          C1=odlxx/rf(is-1,js)
          C2=odlxx/rf(is+1,js)
          C3=odlyy/rf(is,js-1)
          C4=odlyy/rf(is,js+1)
          cof(i,j,1)=C1
          cof(i,j,2)=C2
          cof(i,j,3)=C3
          cof(i,j,4)=C4
          cof(i,j,5)=-(C1+C2+C3+C4)
        end do
      end do
c
# if cdebug
      timing(87)=timing(87)+MPI_WTIME()-tinitial
# endif
      return
      end

      subroutine mgdrelax(sxm,exm,sym,eym,phi,cof,iters,comm2d,myid,
     1                    neighbor,bd,phibc,itype,jtype,IOUT)
#include "compdir.inc"
      include "mpif.h"
      integer sxm,exm,sym,eym,iters,IOUT
      integer comm2d,myid,neighbor(8),bd(8),itype,jtype
      REALN phi(sxm-1:exm+1,sym-1:eym+1)
      REALN cof(sxm-1:exm+1,sym-1:eym+1,6),phibc(4)
c------------------------------------------------------------------------
c Gauss-Seidel point relaxation with Red & Black ordering. Works for
c periodic, Neumann, and Dirichlet boundary conditions.
c
c Code      : mgd2, 2-D parallel multigrid solver
c Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
c Called in : mgdkcyc
c Calls     : mgdbdry, gxch1lin
c------------------------------------------------------------------------
      integer rb,it,ipass,i,j
# if cdebug
      double precision tinitial
      tinitial=MPI_WTIME()
# endif
c
c do iters sweeps in the subdomain
c
      do it=1,iters
        rb=mod(sxm,2)
        do ipass=1,2
          do j=sym,eym
            do i=sxm+rb,exm,2
              phi(i,j)=(cof(i,j,6)-(cof(i,j,1)*phi(i-1,j)
     1                             +cof(i,j,2)*phi(i+1,j)
     2                             +cof(i,j,3)*phi(i,j-1)
     3                             +cof(i,j,4)*phi(i,j+1)))/cof(i,j,5)
            end do
            rb=1-rb
          end do
          rb=1-mod(sxm,2)
# if WMGD
c
c new version: impose Neumann and Dirichlet boundary conditions
c
          call mgdbdry(sxm,exm,sym,eym,phi,bd,phibc,IOUT)
# endif
        end do
      end do
c
c Exchange boundary data only once at the end. Since the number
c of relaxation sweeps at each level is characteristically small
c (1 or 2 are common values), this does not damage the convergence
c rate too badly. Overall, I have found a significant reduction
c in execution time. This also imposes the periodic BCs.
c
      call gxch1lin(phi,comm2d,sxm,exm,sym,eym,neighbor,bd,
     1              itype,jtype,IOUT)
# if WMGD
c
c new version: impose Neumann and Dirichlet boundary conditions
c
      call mgdbdry(sxm,exm,sym,eym,phi,bd,phibc,IOUT)
# endif 
c
# if cdebug
      timing(90)=timing(90)+MPI_WTIME()-tinitial
# endif
      return
      end

      subroutine mgdrestr(sxc,exc,syc,eyc,nxc,nyc,phic,rhsc,
     1                    sxf,exf,syf,eyf,nxf,nyf,phif,cof,
     2                    resf,iresw,comm2d,myid,neighbor,bd,
     3                    itype,jtype,ijtype,IOUT)
#include "compdir.inc"
      include "mpif.h"
      integer sxc,exc,syc,eyc,nxc,nyc,iresw
      integer sxf,exf,syf,eyf,nxf,nyf,IOUT
      integer comm2d,myid,neighbor(8),bd(8),itype,jtype,ijtype
      REALN phic(sxc-1:exc+1,syc-1:eyc+1)
      REALN rhsc(sxc-1:exc+1,syc-1:eyc+1)
      REALN phif(sxf-1:exf+1,syf-1:eyf+1)
      REALN resf(sxf-1:exf+1,syf-1:eyf+1)
      REALN cof(sxf-1:exf+1,syf-1:eyf+1,6)
c------------------------------------------------------------------------
c Calculate the residual and restrict it to the coarser level. In
c the new version, the restriction involves 4 points. In the old
c version, it involves 9 points (5 if half-weighting).
c
c Code      : mgd2, 2-D parallel multigrid solver
c Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
c Called in : mgdkcyc
c Calls     : gxch1lin, gxch1cor
c------------------------------------------------------------------------
      integer i,j,isrt,jsrt,iinc,jinc,ic,jc
# if cdebug
      double precision tinitial
      tinitial=MPI_WTIME()
# endif
c------------------------------------------------------------------------
      do jc=syc-1,eyc+1
        do ic=sxc-1,exc+1
          phic(ic,jc)=0.0d0
          rhsc(ic,jc)=0.0d0
        end do
      end do
c
c calculate residual
c
      do j=syf,eyf
        do i=sxf,exf
          resf(i,j)=cof(i,j,6)-(cof(i,j,1)*phif(i-1,j)
     1                         +cof(i,j,2)*phif(i+1,j)
     2                         +cof(i,j,3)*phif(i,j-1)
     3                         +cof(i,j,4)*phif(i,j+1)
     4                         +cof(i,j,5)*phif(i,j))
        end do
      end do
# if WMGD
c------------------------------------------------------------------------
c new version: calculate the right-hand side at the coarser grid
c level from the averages of the values at the 4 surrounding points;
c if there is no coarsifying in one direction, only 2 points are
c used; no exchange of boundary data is necessary
c
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
            rhsc(ic,jc)=0.25d0*(resf(i,j)+resf(i+1,j)
     1                         +resf(i,j+1)+resf(i+1,j+1))
          end do
        end do
      end if
# else
c------------------------------------------------------------------------
c old version: have to exchange boundary data; if full-weighting, 
c need to exchange also corner points
c
      call gxch1lin(resf,comm2d,sxf,exf,syf,eyf,neighbor,bd,
     1              itype,jtype,IOUT)
      if (iresw.eq.1) then
        call gxch1cor(resf,comm2d,sxf,exf,syf,eyf,neighbor,bd,
     1                ijtype,IOUT)
      end if
c
c restrict it to coarser level
c
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
c
      if (iresw.eq.1) then
c
c use full weighting
c
        j=jsrt
        do jc=syc,eyc
          i=isrt
          do ic=sxc,exc
            rhsc(ic,jc)=0.25d0*resf(i,j)
     1                 +0.125d0*(resf(i+1,j)+resf(i-1,j)
     2                          +resf(i,j+1)+resf(i,j-1))
     3                 +0.0625d0*(resf(i+1,j+1)+resf(i+1,j-1)
     4                           +resf(i-1,j-1)+resf(i-1,j+1))
            i=i+iinc
          end do
          j=j+jinc
        end do
      else if (iresw.eq.2) then
c
c use half-weighting
c
        j=jsrt
        do jc=syc,eyc
          i=isrt
          do ic=sxc,exc
            rhsc(ic,jc)=0.5d0*resf(i,j)
     1                 +0.125d0*(resf(i+1,j)+resf(i-1,j)
     2                          +resf(i,j+1)+resf(i,j-1))
            i=i+iinc
          end do
          j=j+jinc
        end do
      end if
# endif
c
# if cdebug
      timing(91)=timing(91)+MPI_WTIME()-tinitial
# endif
      return
      end

      subroutine mgdrpde(sxm,exm,sym,eym,nxm,nym,cof,xl,yl,bd,IOUT)
#include "compdir.inc"
      include "mpif.h"
      integer sxm,exm,sym,eym,nxm,nym,bd(8),IOUT
      REALN cof(sxm-1:exm+1,sym-1:eym+1,6),xl,yl
c------------------------------------------------------------------------
c Discretize the pde: set the coefficients of the cof matrix. Works
c for periodic, Neumann, and Dirichlet boundary conditions.
c
c cof array:
c
c         cof(4)
c           |
c           |
c cof(1)--cof(5)--cof(2)
c           |
c           |
c         cof(3)
c
c Code      : mgd2, 2-D parallel multigrid solver
c Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
c Called in : mgdsolver
c Calls     : --
c------------------------------------------------------------------------
      REALN dlx,odlxx,dly,odlyy
      integer i,j
# if cdebug
      double precision tinitial
      tinitial=MPI_WTIME()
# endif
c
c calculate off-diagonal terms
c
      dlx=xl/float(nxm-1)
      odlxx=1.0d0/(dlx*dlx)
      dly=yl/float(nym-1)
      odlyy=1.0d0/(dly*dly)
      do j=sym,eym
        do i=sxm,exm
          cof(i,j,1)=odlxx
          cof(i,j,2)=odlxx
          cof(i,j,3)=odlyy
          cof(i,j,4)=odlyy
        end do
      end do
c
c enforce Neumann BCs
c
      if (bd(1).eq.1) then
        do j=sym,eym
          cof(exm,j,2)=0.0d0
        end do
      end if
      if (bd(5).eq.1) then
        do j=sym,eym
          cof(sxm,j,1)=0.0d0
        end do
      end if
      if (bd(3).eq.1) then
        do i=sxm,exm
          cof(i,sym,3)=0.0d0
        end do
      end if
      if (bd(7).eq.1) then
        do i=sxm,exm
          cof(i,eym,4)=0.0d0
        end do
      end if
c
c calculate diagonal term
c
      do j=sym,eym
        do i=sxm,exm
          cof(i,j,5)=-(cof(i,j,1)+cof(i,j,2)+cof(i,j,3)+cof(i,j,4))
        end do
      end do
c
# if cdebug
      timing(82)=timing(82)+MPI_WTIME()-tinitial
# endif
      return
      end

      subroutine mgdrsetf(sxf,exf,syf,eyf,rf,r,IOUT)
#include "compdir.inc"
      include "mpif.h"
      integer sxf,exf,syf,eyf,IOUT
      REALN rf(sxf-1:exf+1,syf-1:eyf+1),r(sxf-1:exf+1,syf-1:eyf+1)
c------------------------------------------------------------------------
c For the old version of the multigrid code, set the fine grid values 
c of the density in the work vector
c
c Code      : mgd2, 2-D parallel multigrid solver
c Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
c Called in : mgdsolver
c Calls     : --
c------------------------------------------------------------------------
      integer i,j
# if cdebug
      double precision tinitial
      tinitial=MPI_WTIME()
# endif
c
      do j=syf-1,eyf+1
        do i=sxf-1,exf+1
          rf(i,j)=r(i,j)
        end do
      end do
c
# if cdebug
      timing(85)=timing(85)+MPI_WTIME()-tinitial
# endif
      return
      end

      subroutine mgdrtrsf(sxc,exc,syc,eyc,nxc,nyc,rc,
     1                    sxf,exf,syf,eyf,nxf,nyf,rf,
     2                    comm2d,myid,neighbor,bd,itype,jtype,IOUT)
#include "compdir.inc"
      include "mpif.h"
      integer sxc,exc,syc,eyc,nxc,nyc,sxf,exf,syf,eyf,nxf,nyf,IOUT
      integer comm2d,myid,neighbor(8),bd(8),itype,jtype
      REALN rc(sxc-1:exc+1,syc-1:eyc+1)
      REALN rf(sxf-1:exf+1,syf-1:eyf+1)
c------------------------------------------------------------------------
c For the old version of the multigrid code, transfer values of the 
c density from a finer to a coarser grid level. It is necessary to
c exchange the boundary density data because the grid "shifts" to
c the right as it becomes coarser. (In the new version of the
c multigrid code, there is no such shift, hence no communication is 
c needed).
c
c Code      : mgd2, 2-D parallel multigrid solver
c Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
c Called in : mgdsolver
c Calls     : gxch1lin
c------------------------------------------------------------------------
      integer i,j,ic,jc,i1,i2,j1,j2
# if cdebug
      double precision tinitial
      tinitial=MPI_WTIME()
# endif
c
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
      do jc=syc,eyc
        j=j1*(2*jc-1)+j2*jc
        do ic=sxc,exc
          i=i1*(2*ic-1)+i2*ic
          rc(ic,jc)=rf(i,j)
        end do
      end do
c
c exchange the boundary values (need only lines, not corner)
c
      call gxch1lin(rc,comm2d,sxc,exc,syc,eyc,neighbor,bd,
     1              itype,jtype,IOUT)
c
# if cdebug
      timing(86)=timing(86)+MPI_WTIME()-tinitial
# endif
      return
      end

      subroutine mgdsetf(sxf,exf,syf,eyf,phi,rhs,phif,rhsf,IOUT)
#include "compdir.inc"
      include "mpif.h"
      integer sxf,exf,syf,eyf,IOUT
      REALN phi(sxf-1:exf+1,syf-1:eyf+1)
      REALN rhs(sxf-1:exf+1,syf-1:eyf+1)
      REALN phif(sxf-1:exf+1,syf-1:eyf+1)
      REALN rhsf(sxf-1:exf+1,syf-1:eyf+1)
c------------------------------------------------------------------------
c Set the fine grid values in the work vector
c
c Code      : mgd2, 2-D parallel multigrid solver
c Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
c Called in : mgdsolver
c Calls     : --
c------------------------------------------------------------------------
      integer i,j
# if cdebug
      double precision tinitial
      tinitial=MPI_WTIME()
# endif
c
      do j=syf-1,eyf+1
        do i=sxf-1,exf+1
          phi(i,j)=phif(i,j)
          rhs(i,j)=rhsf(i,j)
        end do
      end do
c
# if cdebug
      timing(88)=timing(88)+MPI_WTIME()-tinitial
# endif
      return
      end

      subroutine mgdsolver(isol,sx,ex,sy,ey,phif,rhsf,r,ngrid,work,
     1                     maxcy,tolmax,kcycle,iprer,ipost,iresw,
     2                     xl,yl,rro,nx,ny,comm2d,myid,neighbor,
     3                     bd,phibc,iter,nprscr,IOUT,nerror)
#include "compdir.inc"
      include "mpif.h"
      integer isol,sx,ex,sy,ey,ngrid,nx,ny,IOUT
      integer maxcy,kcycle,iprer,ipost,iresw
      REALN phif(sx-1:ex+1,sy-1:ey+1),rhsf(sx-1:ex+1,sy-1:ey+1)
      REALN r(sx-1:ex+1,sy-1:ey+1)
      REALN work(*),tolmax,xl,yl,rro,phibc(4,20)
      integer comm2d,myid,neighbor(8),bd(8),iter,nerror
      logical nprscr
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
c Parallel multigrid solver in 2-D cartesian coordinates for the
c elliptic equation:      div(cof(grad(phif)))=rhsf
c
c isol=1 -> density
c isol=2 -> pressure
c 
c Written for periodic, wall (Neumann), and constant value
c (Dirichlet) BCs. Tested roughly for all these BCs. There are 
c two versions of the multigrid code, which are separated by the 
c compiler directive WMGD set in 'compdir.inc'. The old version 
c (WMGD=0) corresponds to the original, "traditional" grid setup, 
c and works well when all boundary conditions are periodic. When one 
c of the BCs is not periodic, must compile with the new version 
c (WMGD=1), which uses a different grid setup and new restriction 
c and correction operators (see 'mgdinit' for more details). It is 
c less accurate (or ,equivalently, slower to converge) for the case 
c of all-periodic BCs than the old version, but is better able to 
c handle wall BCs.
c
c Notes: - the values of the rhs contained in the array rhsf are
c          transferred to the work vector and the memory of the 
c          rhsf array is then used to store the residuals at the
c          different levels; therefore, it does not preserve its
c          initial values
c        - with the initialization and clean-up routines mgdinit
c          and mgdend, this multigrid code is self-standing
c
c Code      : mgd2, 2-D parallel multigrid solver
c Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
c Called in : main
c Calls     : mgdrpde, mgdpfpde, mgdphpde, mgdrsetf, mgdppde, mgdrtrsf,
c              -> discretize the pde
c             mgdsetf
c               -> set the initial guesses and the right-hand side
c             mgdkcyc, mgderr,
c               -> do the actual cycling
c             gscale, gxch1lin, gxch1cor,
c               -> rescale pressure and density around average values
c------------------------------------------------------------------------
      REALN avo,acorr
      integer sxf,exf,syf,eyf,nxf,nyf,sxc,exc,syc,eyc,nxc,nyc
      integer ipf,irf,ipc,irc,sxm,exm,sym,eym,nxm,nym,ip,ic,kcur
      integer itype,jtype,ijtype,lev,ir1,ir2
      integer icf, k
      REALN relmax
# if cdebug
      double precision tinitial

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
          nxm=nxk(k)
          nym=nyk(k)
          ic=kcbgn(k)
          call mgdrpde(sxm,exm,sym,eym,nxm,nym,work(ic),xl,yl,bd,IOUT)
        end do
      else
c
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
        nxf=nxk(ngrid)
        nyf=nyk(ngrid)
        icf=kcbgn(ngrid)
        call mgdpfpde(sxf,exf,syf,eyf,nxf,nyf,work(icf),r,xl,yl,bd,
     1                IOUT)
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
          nxm=nxk(k)
          nym=nyk(k)
          ic=kcbgn(k)
          call mgdphpde(sxm,exm,sym,eym,nxm,nym,work(ic),
     1                  sx,ex,sy,ey,nxf,nyf,r,bd,xl,yl,IOUT)
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
          ir2=kcbgn(ngrid)+5*(exf-sxf+3)*(eyf-syf+3)
          sxf=sxr(ngrid-1)
          exf=exr(ngrid-1)
          syf=syr(ngrid-1)
          eyf=eyr(ngrid-1)
          nxf=nxr(ngrid-1)
          nyf=nyr(ngrid-1)
          lev=1
          call mgdrsetf(sxf,exf,syf,eyf,work(ir1),r,IOUT)
c
c for the levels k=ngrid-1,1, calculate the coefficients from the
c densities stored in the arrays work(ir1) and work(ir2)
c for the levels k=ngrid-1,2, transfer to the level below the values
c of the density needed there to determine the coefficients; exchange
c of the boundary density data is necessary
c
          do k=ngrid-1,1,-1
            sxm=sxk(k)
            exm=exk(k)
            sym=syk(k)
            eym=eyk(k)
            nxm=nxk(k)
            nym=nyk(k)
            ic=kcbgn(k)
            if (lev.eq.1) then
              call mgdppde(sxm,exm,sym,eym,nxm,nym,work(ic),
     1                     sxf,exf,syf,eyf,work(ir1),xl,yl,bd,IOUT)
            else
              call mgdppde(sxm,exm,sym,eym,nxm,nym,work(ic),
     1                     sxf,exf,syf,eyf,work(ir2),xl,yl,bd,IOUT)
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
                call mgdrtrsf(sxc,exc,syc,eyc,nxc,nyc,work(ir2),
     1                        sxf,exf,syf,eyf,nxf,nyf,work(ir1),
     2                        comm2d,myid,neighbor,bd,itype,jtype,IOUT)
                lev=2
              else
                call mgdrtrsf(sxc,exc,syc,eyc,nxc,nyc,work(ir1),
     1                        sxf,exf,syf,eyf,nxf,nyf,work(ir2),
     2                        comm2d,myid,neighbor,bd,itype,jtype,IOUT)
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
c------------------------------------------------------------------------
c set phi,rhsf in work at the finest grid level
c
      sxf=sxk(ngrid)
      exf=exk(ngrid)
      syf=syk(ngrid)
      eyf=eyk(ngrid)
      ipf=kpbgn(ngrid)
      irf=kcbgn(ngrid)+5*(exf-sxf+3)*(eyf-syf+3)
      call mgdsetf(sxf,exf,syf,eyf,work(ipf),work(irf),phif,rhsf,IOUT)
c------------------------------------------------------------------------
c cycling at kcur=ngrid level
c
      kcur=ngrid
      do iter=1,maxcy
        call mgdkcyc(work,rhsf,kcur,kcycle,iprer,ipost,iresw,
     1               comm2d,myid,neighbor,bd,phibc,IOUT)
        sxm=sxk(ngrid)
        exm=exk(ngrid)
        sym=syk(ngrid)
        eym=eyk(ngrid)
        ip=kpbgn(ngrid)
        call mgderr(relmax,sxm,exm,sym,eym,phif,work(ip),comm2d,IOUT)
        if (relmax.le.tolmax) goto 1000
      end do
c------------------------------------------------------------------------
c if not converged in maxcy cycles, issue an error message and quit
c
      if (myid.eq.0) write(IOUT,100) maxcy,relmax
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

      call gscale(sx,ex,sy,ey,phif,avo,acorr,comm2d,nx,ny,IOUT)
c
c exchange boundary data and impose periodic BCs
c
      call gxch1lin(phif,comm2d,sx,ex,sy,ey,neighbor,bd,
     1              ikdatatype(ngrid),jkdatatype(ngrid),IOUT)
      call gxch1cor(phif,comm2d,sxm,exm,sym,eym,neighbor,bd,
     1              ijkdatatype(ngrid),IOUT)
# if WMGD
c
c impose wall and Dirichlet BCs
c
      call mgdbdry(sx,ex,sy,ey,phif,bd,phibc(:,1),IOUT)
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

      subroutine grid1_type(itype,jtype,ijtype,realtype,sx,ex,sy,ey,
     1                      IOUT)
#include "compdir.inc"
      include "mpif.h"
      integer ierr
      integer itype,jtype,ijtype,realtype,sx,ex,sy,ey,IOUT
c------------------------------------------------------------------------
c Define the 3 derived datatypes needed to communicate the boundary
c data of (sx-1:ex+1,sy-1:ey+1) arrays between 'myid' and its 8
c neighbors
c
c Code      : tmgd2
c Called in : mgdinit
c Calls     : MPI_TYPE_CONTIGUOUS, MPI_TYPE_COMMIT, MPI_TYPE_VECTOR
c------------------------------------------------------------------------
      integer ier
# if cdebug
      double precision tinitial
      tinitial=MPI_WTIME()
# endif
c
c datatype for one row
c
      call MPI_TYPE_CONTIGUOUS(ex-sx+1,realtype,itype,ierr)
      call MPI_TYPE_COMMIT(itype,ierr)
c
c datatype for one column
c
      call MPI_TYPE_VECTOR(ey-sy+1,1,ex-sx+3,realtype,jtype,ierr)
      call MPI_TYPE_COMMIT(jtype,ierr)
c
c datatype for one 1*1 corner
c
      call MPI_TYPE_CONTIGUOUS(1,realtype,ijtype,ierr)
      call MPI_TYPE_COMMIT(ijtype,ierr)
c
# if cdebug
      timing(7)=timing(7)+MPI_WTIME()-tinitial
# endif
      return
      end

      subroutine gxch1cor(a,comm2d,sx,ex,sy,ey,neighbor,bd,
     1                    ijdatatype,IOUT)
#include "compdir.inc"
      include "mpif.h"
      integer sx,ex,sy,ey,IOUT
      REALN a(sx-1:ex+1,sy-1:ey+1)
      integer comm2d,neighbor(8),bd(8),ijdatatype
      integer ierr
c------------------------------------------------------------------------
c Subroutine to exchange one corner point of boundary data between 
c "diagonally" neighboring processes. This subroutineccan be used to 
c exchange scalar as well as vector variables; it suffices to pass the 
c right argument for the datatypes:
c ijdatatype -> r, p, tmp...
c ij11datatype -> (u,v)
c ij12datatype -> (ut,vt)
c
c 'neighbor' and 'bd' arrays:
c
c     6 |           | 8
c       |           |
c    ------------------
c       |           |
c       |           |
c       |   myid    |   
c       |           |
c       |           |
c    ------------------
c       |           |
c     4 |           | 2
c
c Code      : tmgd2
c Called in : mgdrestr, mgdsolver
c Calls     : MPI_ISEND, MPI_IRECV, MPI_WAITALL (non-blocking version)
c             MPI_SENDRECV (blocking version)
c----------------------------------------------------------------------- 
# if NBLOCKGR
      integer req(8),status(MPI_STATUS_SIZE,8),ireq
# else
      integer status(MPI_STATUS_SIZE)
# endif
# if cdebug
# if NBLOCKGR
      integer nc
# endif
      double precision tinitial
      tinitial=MPI_WTIME()
# endif
# if NBLOCKGR
c--------------------------non-blocking----------------------------------
      ireq=0
c
c send to 2
c
      if (bd(2).eq.0) then
        ireq=ireq+1
        call MPI_ISEND(a(ex,sy),1,ijdatatype,neighbor(2),
     1                 0,comm2d,req(ireq),ierr)
      end if
c
c receive from 6
c
      if (bd(6).eq.0) then
        ireq=ireq+1
        call MPI_IRECV(a(sx-1,ey+1),1,ijdatatype,neighbor(6),
     1                 0,comm2d,req(ireq),ierr)
      end if
c
c send to 4
c
      if (bd(4).eq.0) then
        ireq=ireq+1
        call MPI_ISEND(a(sx,sy),1,ijdatatype,neighbor(4),
     1                 1,comm2d,req(ireq),ierr)
      end if
c
c receive from 8
c
      if (bd(8).eq.0) then
        ireq=ireq+1
        call MPI_IRECV(a(ex+1,ey+1),1,ijdatatype,neighbor(8),
     1                 1,comm2d,req(ireq),ierr)
      end if
c
c send to 6
c
      if (bd(6).eq.0) then
        ireq=ireq+1
        call MPI_ISEND(a(sx,ey),1,ijdatatype,neighbor(6),
     1                 1,comm2d,req(ireq),ierr)
      end if
c
c receive from 2
c
      if (bd(2).eq.0) then
        ireq=ireq+1
        call MPI_IRECV(a(ex+1,sy-1),1,ijdatatype,neighbor(2),
     1                 1,comm2d,req(ireq),ierr)
      end if
c
c send to 8
c
      if (bd(8).eq.0) then
        ireq=ireq+1
        call MPI_ISEND(a(ex,ey),1,ijdatatype,neighbor(8),
     1                 0,comm2d,req(ireq),ierr)
      end if
c
c receive from 4
c
      if (bd(4).eq.0) then
        ireq=ireq+1
        call MPI_IRECV(a(sx-1,sy-1),1,ijdatatype,neighbor(4),
     1                 0,comm2d,req(ireq),ierr)
      end if
c
c wait for all the messages to be sent and received before going on.
c
      call MPI_WAITALL(ireq,req,status,ierr)
# if cdebug
      nc=4-(bd(2)+bd(4)+bd(6)+bd(8))
      nisend(2,1)=nisend(2,1)+nc
      nirecv(2,1)=nirecv(2,1)+nc
      nwaitall=nwaitall+1
# endif
# else
c----------------------------blocking------------------------------------
c send to 2 and receive from 6
c
      call MPI_SENDRECV(a(ex,sy),1,ijdatatype,neighbor(2),0,
     1                  a(sx-1,ey+1),1,ijdatatype,neighbor(6),0,
     2                  comm2d,status,ierr)
c
c send to 4 and receive from 8
c
      call MPI_SENDRECV(a(sx,sy),1,ijdatatype,neighbor(4),1,
     1                  a(ex+1,ey+1),1,ijdatatype,neighbor(8),1,
     2                  comm2d,status,ierr)
c
c send to 6 and receive from 2
c
      call MPI_SENDRECV(a(sx,ey),1,ijdatatype,neighbor(6),1,
     1                  a(ex+1,sy-1),1,ijdatatype,neighbor(2),1,
     2                  comm2d,status,ierr)
c
c send to 8 and receive from 4
c
      call MPI_SENDRECV(a(ex,ey),1,ijdatatype,neighbor(8),0,
     1                  a(sx-1,sy-1),1,ijdatatype,neighbor(4),0,
     2                  comm2d,status,ierr)
# if cdebug
      nsendrecv(2,1)=nsendrecv(2,1)+4
# endif
# endif
c
# if cdebug
      timing(60)=timing(60)+MPI_WTIME()-tinitial
# endif
      return
      end

      subroutine gxch1lin(a,comm2d,sx,ex,sy,ey,neighbor,bd,
     1                    idatatype,jdatatype,IOUT)
#include "compdir.inc"
      include "mpif.h"
      integer sx,ex,sy,ey,IOUT
      REALN a(sx-1:ex+1,sy-1:ey+1)
      integer comm2d,neighbor(8),bd(8),idatatype,jdatatype
      integer ierr
c------------------------------------------------------------------------
c Subroutine to exchange one-lines (one row and one column) of 
c boundary data between "directly" neighboring processes. This subroutine 
c can be used to exchange scalar as well as vector variables; it suffices 
c to pass the right argument for the datatypes: 
c idatatype, jdatatype -> r, p, tmp...
c i11datatype, j11datatype -> (u,v)
c i12datatype, j12datatype -> (ut,vt)
c
c If other quantities need to be passed, appropriate datatypes for
c them have to be defined first in type_mpi.
c
c 'neighbor' and 'bd' arrays:
c
c       |     7     | 
c       |           |
c    ------------------
c       |           |
c       |           |
c     5 |   myid    | 1 
c       |           |
c       |           |
c    ------------------
c       |           |
c       |     3     |  
c
c Code      : tmgd2
c Called in : mgdrelax, mgdrelax, mgdrestr, mgdrtrsf
c Calls     : MPI_ISEND, MPI_IRECV, MPI_WAITALL (non-blocking version)
c             MPI_SENDRECV (blocking version)
c------------------------------------------------------------------------
# if NBLOCKGR
      integer req(8),status(MPI_STATUS_SIZE,8),ireq
# else
      integer status(MPI_STATUS_SIZE)
# endif
# if cdebug
# if NBLOCKGR
      integer nc
# endif
      double precision tinitial
      tinitial=MPI_WTIME()
# endif
# if NBLOCKGR
c--------------------------non-blocking----------------------------------
      ireq=0
c
c send to 1
c
      if (bd(1).eq.0) then
        ireq=ireq+1
        call MPI_ISEND(a(ex,sy),1,jdatatype,neighbor(1),
     1                 0,comm2d,req(ireq),ierr)
      end if
c
c receive from 5
c
      if (bd(5).eq.0) then
        ireq=ireq+1
        call MPI_IRECV(a(sx-1,sy),1,jdatatype,neighbor(5),
     1                 0,comm2d,req(ireq),ierr)
      end if
c
c send to 3
c
      if (bd(3).eq.0) then
        ireq=ireq+1
        call MPI_ISEND(a(sx,sy),1,idatatype,neighbor(3),
     1                 1,comm2d,req(ireq),ierr)
      end if
c
c receive from 7
c 
      if (bd(7).eq.0) then
        ireq=ireq+1
        call MPI_IRECV(a(sx,ey+1),1,idatatype,neighbor(7),
     1                 1,comm2d,req(ireq),ierr)
      end if
c
c send to 5
c
      if (bd(5).eq.0) then
        ireq=ireq+1
        call MPI_ISEND(a(sx,sy),1,jdatatype,neighbor(5),
     1                 1,comm2d,req(ireq),ierr)
      end if
c
c receive from 1
c
      if (bd(1).eq.0) then
        ireq=ireq+1
        call MPI_IRECV(a(ex+1,sy),1,jdatatype,neighbor(1),
     1                 1,comm2d,req(ireq),ierr)
      end if
c
c send to 7
c
      if (bd(7).eq.0) then
        ireq=ireq+1
        call MPI_ISEND(a(sx,ey),1,idatatype,neighbor(7),
     1                 0,comm2d,req(ireq),ierr)
      end if
c
c receive from 3
c
      if (bd(3).eq.0) then
        ireq=ireq+1
        call MPI_IRECV(a(sx,sy-1),1,idatatype,neighbor(3),
     1                 0,comm2d,req(ireq),ierr)
      end if
c
c wait for all the messages to be sent and received before going on.
c
      call MPI_WAITALL(ireq,req,status,ierr)
# if cdebug
      nc=4-(bd(1)+bd(3)+bd(5)+bd(7))
      nisend(1,1)=nisend(1,1)+nc
      nirecv(1,1)=nirecv(1,1)+nc
      nwaitall=nwaitall+1
# endif
# else
c----------------------------blocking------------------------------------
c send to 1 and receive from 5
c
      call MPI_SENDRECV(a(ex,sy),1,jdatatype,neighbor(1),0,
     1                  a(sx-1,sy),1,jdatatype,neighbor(5),0,
     2                  comm2d,status,ierr)
c
c send to 3 and receive from 7
c
      call MPI_SENDRECV(a(sx,sy),1,idatatype,neighbor(3),1,
     1                  a(sx,ey+1),1,idatatype,neighbor(7),1,
     2                  comm2d,status,ierr)
c
c send to 5 and receive from 1
c
      call MPI_SENDRECV(a(sx,sy),1,jdatatype,neighbor(5),1,
     1                  a(ex+1,sy),1,jdatatype,neighbor(1),1,
     2                  comm2d,status,ierr)
c
c send to 7 and receive from 3
c
      call MPI_SENDRECV(a(sx,ey),1,idatatype,neighbor(7),0,
     1                  a(sx,sy-1),1,idatatype,neighbor(3),0,
     2                  comm2d,status,ierr)
# if cdebug
      nsendrecv(1,1)=nsendrecv(1,1)+4
# endif
# endif
c
# if cdebug
      timing(59)=timing(59)+MPI_WTIME()-tinitial
# endif
      return
      end

      subroutine gscale(sx,ex,sy,ey,a,avo,acorr,comm2d,nx,ny,IOUT)
# include "compdir.inc"
      include "mpif.h"
      integer sx,ex,sy,ey,nx,ny,IOUT
      REALN a(sx-1:ex+1,sy-1:ey+1),avo,acorr
      integer comm2d
c------------------------------------------------------------------------
c Rescale the field a so that its average inside the domain
c remains constant and equal to avo. For the density,avo should
c be rro, this ensures conservation of mass. For the pressure,
c avo should be 0 so that the average pressure does not drift
c away from 0, which is the initial value.
c
c Code      : tmgd2
c Called in : mgdsolver
c Calls     : MPI_ALLREDUCE
c------------------------------------------------------------------------
      REALN avloc,av
      integer i,j,ierr
# if cdebug
      double precision tinitial
      tinitial=MPI_WTIME()
# endif
c
c determine average value
c
      avloc=0.0d0
      do j=sy,ey
        do i=sx,ex
          avloc=avloc+a(i,j)
        end do
      end do
c
c global reduce across all process
c
# if double_precision
      call MPI_ALLREDUCE(avloc,av,1,MPI_DOUBLE_PRECISION,MPI_SUM,
     1                   comm2d,ierr)
# else
      call MPI_ALLREDUCE(avloc,av,1,MPI_REAL,MPI_SUM,comm2d,ierr)
# endif
# if cdebug
      nallreduce=nallreduce+1
# endif
   
      if ( abs(av) > epsilon(1d0)) av=av/float(nx*ny)
c
c do correction
c
      acorr=avo-av
      do j=sy,ey
        do i=sx,ex
          a(i,j)=a(i,j)+acorr
        end do
      end do
c
# if cdebug
      timing(49)=timing(49)+MPI_WTIME()-tinitial
# endif
      return
      end 

