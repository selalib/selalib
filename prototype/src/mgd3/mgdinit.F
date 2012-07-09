      subroutine mgdinit(vbc,phibc,ixp,jyq,kzr,iex,jey,kez,ngrid,nxp2,
     1                   nyp2,nzp2,sx,ex,sy,ey,sz,ez,realtype,nxprocs,
     2                   nyprocs,nzprocs,nwork,ibdry,jbdry,kbdry,myid,
     3                   IOUT,nerror)
# include "compdir.inc"
      include "mpif.h"
      integer ixp,jyq,kzr,iex,jey,kez,ngrid,nxp2,nyp2,nzp2
      integer sx,ex,sy,ey,sz,ez,realtype,nxprocs,nyprocs,nzprocs
      integer nwork,ibdry,jbdry,kbdry,myid,IOUT,nerror
      REALN vbc(6),phibc(6,20)
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
c Initialize the parallel multigrid solver: subdomain indices,
c MPI datatypes, boundary values for Dirichlet boundaries.
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
c ngrid=iex=jey=kzr=1, with ixp=nx, jyq=ny, and kez=nz. In this case, 
c the code never enters 'mgdrestr' and 'mgdcor' and all it does is 
c Gauss-Seidel iterate at the finest grid level. This can be useful 
c as a preliminary check.
c
c Note: some memory could be saved by noting that the cof arrays
c need be dimensioned (sxm:exm,sym:eym,szm:ezm) and not
c (sxm-1:exm+1,sym-1:eym+1,szm-1:ezm+1)... Probably not too difficult 
c to make the change
c
c Code      : mgd3, 3-D parallel multigrid solver
c Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
c Called in : main
c Calls     : grid1_type
c------------------------------------------------------------------------
      integer i,j,k,nxf,nyf,nzf,nxm,nym,nzm,kps,nxc,nyc,nzc
      integer sxm,exm,sym,eym,szm,ezm
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
        nzk(k)=0
        sxk(k)=0
        exk(k)=0
        syk(k)=0
        eyk(k)=0
        szk(k)=0
        ezk(k)=0
        kpbgn(k)=0
        kcbgn(k)=0
        do j=1,7
          kdatatype(j,k)=MPI_DATATYPE_NULL
          rdatatype(j,k)=MPI_DATATYPE_NULL
        end do
        sxi(k)=0
        exi(k)=0
        syi(k)=0
        eyi(k)=0
        szi(k)=0
        ezi(k)=0
        nxr(k)=0
        nyr(k)=0
        nzr(k)=0
        sxr(k)=0
        exr(k)=0
        syr(k)=0
        eyr(k)=0
        szr(k)=0
        ezr(k)=0
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
      if (mod(nzp2-2,nzprocs).ne.0) then
        write(IOUT,120) nzp2-2,nzprocs
        nerror=1
        return
      end if
100   format(/,'ERROR in mgdinit: nx=',i3,' is not a multiple of ',
     1       'nxprocs=',i3,/,'cannot use the new version of the ',
     2       'multigrid code',/)
110   format(/,'ERROR in mgdinit: ny=',i3,' is not a multiple of ',
     1       'nyprocs=',i3,/,'cannot use the new version of the ',
     2       'multigrid code',/)
120   format(/,'ERROR in mgdinit: nz=',i3,' is not a multiple of ',
     1       'nzprocs=',i3,/,'cannot use the new version of the ',
     2       'multigrid code',/)
# else
c
c check that the old version is not used with non-periodic BCs
c
      if (ibdry.ne.0.or.jbdry.ne.0.or.kbdry.ne.0) then
        write(IOUT,130) ibdry,jbdry,kbdry
        nerror=1
        return
      end if
130   format(/,'ERROR in mgdinit: ibdry=',i2,' jbdry=',i2,' kbdry=',i2,
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
        write(IOUT,140) nxp2-1,i
        nerror=1
        return
      end if
      j=jyq*2**(jey-1)+1
      if ((nyp2-1).ne.j) then
        write(IOUT,150) nyp2-1,j
        nerror=1
        return
      end if
      k=kzr*2**(kez-1)+1
      if ((nzp2-1).ne.k) then
        write(IOUT,160) nzp2-1,k
        nerror=1
        return
      end if
140   format(/,'ERROR in mgdinit: nxp1=',i3,' <> ixp*2**(iex-1)+1=',
     1       i3,/,'-> adjust the multigrid parameters ixp and iex',
     2       ' in main',/)
150   format(/,'ERROR in mgdinit: nyp1=',i3,' <> jyq*2**(jey-1)+1=',
     1       i3,/,'-> adjust the multigrid parameters jyq and jey',
     2       ' in main',/)
160   format(/,'ERROR in mgdinit: nzp1=',i3,' <> kzr*2**(kez-1)+1=',
     1       i3,/,'-> adjust the multigrid parameters kzr and kez',
     2       ' in main',/)
c
c check that the number of points at the coarser level is not smaller
c than the number of processes in either direction
c
      if (ixp.lt.nxprocs) then
        write(IOUT,170) ixp,nxprocs
        nerror=1
        return
      end if
      if (jyq.lt.nyprocs) then
        write(IOUT,180) jyq,nyprocs
        nerror=1
        return
      end if
      if (kzr.lt.nzprocs) then
        write(IOUT,190) kzr,nzprocs
        nerror=1
        return
      end if
170   format(/,'ERROR in mgdinit: ixp=',i3,' < nxprocs=',i3,/,
     1       ' there must be at least one grid point at the ',
     2       'coarsest grid level',/,
     3       '-> increase ixp and decrease iex correspondingly',
     4       ' in main',/)
180   format(/,'ERROR in mgdinit: jyq=',i3,' < nyprocs=',i3,/,
     1       ' there must be at least one grid point at the ',
     2       'coarsest grid level',/,
     3       '-> increase jyq and decrease jey correspondingly',
     4       ' in main',/)
190   format(/,'ERROR in mgdinit: kzr=',i3,' < nzprocs=',i3,/,
     1       ' there must be at least one grid point at the ',
     2       'coarsest grid level',/,
     3       '-> increase kzr and decrease kez correspondingly',
     4       ' in main',/)
c
c check that coarsifying takes place in all directions at the finest
c grid level
c
      if (ngrid.gt.1) then
        if (iex.eq.1) then
          write(IOUT,200) ngrid,iex
          nerror=1
          return
        end if
        if (jey.eq.1) then
          write(IOUT,210) ngrid,jey
          nerror=1
          return
        end if
        if (kez.eq.1) then
          write(IOUT,220) ngrid,kez
          nerror=1
          return
        end if
      end if
200   format(/,'ERROR in mgdinit: ngrid=',i3,' iex=',i3,
     1       /,'no coarsifying at the finest grid level in x-direction',
     2       /,'this is not allowed by the mutligrid code',/)
210   format(/,'ERROR in mgdinit: ngrid=',i3,' jey=',i3,
     1       /,'no coarsifying at the finest grid level in y-direction',
     2       /,'this is not allowed by the mutligrid code',/)
220   format(/,'ERROR in mgdinit: ngrid=',i3,' kez=',i3,
     1       /,'no coarsifying at the finest grid level in z-direction',
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
        nzk(k)=kzr*2**(max(k+kez-ngrid,1)-1)+1
      end do
c
c for all grid levels, set the indices of the subdomain the process
c 'myid' work on, as well as the datatypes needed for the exchange
c of boundary data:
c
      nxf=nxk(ngrid)
      nyf=nyk(ngrid)
      nzf=nzk(ngrid)
      sxk(ngrid)=sx
      exk(ngrid)=ex
      syk(ngrid)=sy
      eyk(ngrid)=ey
      szk(ngrid)=sz
      ezk(ngrid)=ez
      call grid1_type(kdatatype(1,ngrid),realtype,sxk(ngrid),exk(ngrid),
     1                syk(ngrid),eyk(ngrid),szk(ngrid),ezk(ngrid),IOUT)
      do k=ngrid-1,1,-1
        nxm=nxk(k)
        nym=nyk(k)
        nzm=nzk(k)
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
        if (nzm.lt.nzf) then
          szk(k)=szk(k+1)/2+1
          ezk(k)=(ezk(k+1)-1)/2+1
        else
          szk(k)=szk(k+1)
          ezk(k)=ezk(k+1)
        end if
        nxf=nxm
        nyf=nym
        nzf=nzm
        call grid1_type(kdatatype(1,k),realtype,sxk(k),exk(k),
     1                  syk(k),eyk(k),szk(k),ezk(k),IOUT)
      end do
c
c set work space indices for phi, cof at each grid level, and check
c that there is sufficient work space
c
      kps=1
      do k=ngrid,1,-1
        sxm=sxk(k)
        exm=exk(k)
        sym=syk(k)
        eym=eyk(k)
        szm=szk(k)
        ezm=ezk(k)
        kpbgn(k)=kps
        kcbgn(k)=kpbgn(k)+(exm-sxm+3)*(eym-sym+3)*(ezm-szm+3)
        kps=kcbgn(k)+8*(exm-sxm+3)*(eym-sym+3)*(ezm-szm+3)
      end do
      if (kps.gt.nwork) then
        write(IOUT,230) kps,nwork,myid
        nerror=1
        return
230     format(/,'ERROR in mgdinit: not enough work space',/,
     1         ' kps=',i12,' nwork=',i12,' myid: ',i3,/,
     2         ' -> put the formula for nwork in main in ',
     3         'comments',/,'    and set nwork to the value of kps',/)
      else
        write(IOUT,240) kps,nwork
240     format(/,'WARNING in mgdinit: kps=',i10,' nwork=',i10,
     1         /,'can optimize the amount of memory needed by ',
     2           'the multigrid code',/,'by putting the formula ',
     3           'for nwork into comments and setting',/,'nwork ',
     4           'to the value of kps',/)
      end if
# if WMGD
c------------------------------------------------------------------------
c For the new version of the multigrid code, set the boundary values 
c to be used for the Dirichlet boundaries. It is possible to assign
c 6 different constant values to the 6 different sides. The values are
c assigned at the finest grid level, zero is assigned at all levels
c below
c 
c vbc, phibc:
c                                      k
c                   vbc(6)             
c                   /                  ^  ^
c       -----vbc(4)/-----              | / i
c       |         /     |              |/
c       |        /      |         -----/-----> j
c     vbc(3)----+-----vbc(1)          /|
c       |      /        |            / |
c       |     /         |           /  |
c       -----/vbc(2)----|
c           /
c         vbc(5)
c
      do j=1,6
        phibc(j,ngrid)=vbc(j)
      end do
      do k=ngrid-1,1,-1
        do j=1,6
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
      szi(ngrid)=szk(ngrid)-1
      ezi(ngrid)=ezk(ngrid)+1
      nxf=nxk(ngrid)
      nyf=nyk(ngrid)
      nzf=nzk(ngrid)
      do k=ngrid-1,1,-1
        nxc=nxk(k)
        nyc=nyk(k)
        nzc=nzk(k)
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
        if (nzc.lt.nzf) then
          szi(k)=szk(k)-1+mod(szk(k+1),2)
          ezi(k)=ezk(k)+1-mod(ezk(k+1),2)
        else
          szi(k)=szk(k)-1
          ezi(k)=ezk(k)+1
        end if
        nxf=nxc
        nyf=nyc
        nzf=nzc
      end do
c------------------------------------------------------------------------
c set indices for determining the coefficients in the elliptic
c equation div(cof*grad(P))=rhs. Used only when solving for the
c pressure. When setting these coefficients at level k, need
c the values of the density at midpoints, i.e. at level k+1
c (if coarsifying takes place between the levels k and k+1).
c If coarsifying took place at all levels, the array cof could
c be used as temporary storage space for the densities, with
c cof(*,*,8) at level k+1 giving the values cof(*,*,1->7) at level
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
        nzf=nzk(k+1)
        nxc=nxk(k)
        nyc=nyk(k)
        nzc=nzk(k)
        if (nxc.lt.nxf) then
          sxr(k)=sxk(k+1)
          exr(k)=exk(k+1)
          nxr(k)=nxk(k+1)
          sxm=sxr(k)
          exm=exr(k)
          nxm=nxr(k)
        else
          if (k.eq.(ngrid-1)) then
            write(IOUT,250) ngrid,ngrid-1
250         format('ERROR in mgdinit: no coarsifying between level ',
     1             i3,' and level ',i3,/,', the current version of ',
     2             'the code cannot cope with that',/,
     3             ' -> change the parameters (ixp,jyq,kzr) and',
     4             ' (iex,jyq,kez) in main',/)
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
            write(IOUT,250) ngrid,ngrid-1
            nerror=1
            return
          end if
          syr(k)=sym
          eyr(k)=eym
          nyr(k)=nym
        end if
        if (nzc.lt.nzf) then
          szr(k)=szk(k+1)
          ezr(k)=ezk(k+1)
          nzr(k)=nzk(k+1)
          szm=szr(k)
          ezm=ezr(k)
          nzm=nzr(k)
        else
          if (k.eq.(ngrid-1)) then
            write(IOUT,250) ngrid,ngrid-1
            nerror=1
            return
          end if
          szr(k)=szm
          ezr(k)=ezm
          nzr(k)=nzm
        end if
        call grid1_type(rdatatype(1,k),realtype,sxr(k),exr(k),
     1                  syr(k),eyr(k),szr(k),ezr(k),IOUT)
      end do
# endif
c
# if cdebug
      timing(81)=timing(81)+MPI_WTIME()-tinitial
# endif
      return
      end

