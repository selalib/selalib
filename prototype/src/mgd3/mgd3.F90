subroutine mgdinit(vbc,phibc,ixp,jyq,kzr,iex,jey,kez,ngrid,nxp2,&
                   nyp2,nzp2,sx,ex,sy,ey,sz,ez,realtype,nxprocs,&
                   nyprocs,nzprocs,nwork,ibdry,jbdry,kbdry,myid,&
                   IOUT,nerror)
# include "compdir.inc"
include "mpif.h"
integer :: ixp,jyq,kzr,iex,jey,kez,ngrid,nxp2,nyp2,nzp2
integer :: sx,ex,sy,ey,sz,ez,realtype,nxprocs,nyprocs,nzprocs
integer :: nwork,ibdry,jbdry,kbdry,myid,IOUT,nerror
REALN   :: vbc(6),phibc(6,20)
 
integer :: nxk,nyk,nzk,sxk,exk,syk,eyk,szk,ezk
integer :: kpbgn,kcbgn,kdatatype
integer :: sxi,exi,syi,eyi,szi,ezi
integer :: nxr,nyr,nzr,sxr,exr,syr,eyr,szr,ezr
integer :: rdatatype
common/mgd/nxk(20),nyk(20),nzk(20),                          &
           sxk(20),exk(20),syk(20),eyk(20),szk(20),ezk(20),  &
           kpbgn(20),kcbgn(20),kdatatype(7,20),              &
           sxi(20),exi(20),syi(20),eyi(20),szi(20),ezi(20),  &
           nxr(20),nyr(20),nzr(20),sxr(20),exr(20),syr(20),  &
           eyr(20),szr(20),ezr(20),rdatatype(7,20)
!------------------------------------------------------------------------
! Initialize the parallel multigrid solver: subdomain indices,
! MPI datatypes, boundary values for Dirichlet boundaries.
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
! ngrid=iex=jey=kzr=1, with ixp=nx, jyq=ny, and kez=nz. In this case, 
! the code never enters 'mgdrestr' and 'mgdcor' and all it does is 
! Gauss-Seidel iterate at the finest grid level. This can be useful 
! as a preliminary check.
!
! Note: some memory could be saved by noting that the cof arrays
! need be dimensioned (sxm:exm,sym:eym,szm:ezm) and not
! (sxm-1:exm+1,sym-1:eym+1,szm-1:ezm+1)... Probably not too difficult 
! to make the change
!
! Code      : mgd3, 3-D parallel multigrid solver
! Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
! Called in : main
! Calls     : grid1_type
!------------------------------------------------------------------------
integer :: i,j,k,nxf,nyf,nzf,nxm,nym,nzm,kps,nxc,nyc,nzc
integer :: sxm,exm,sym,eym,szm,ezm
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
if (mod(nzp2-2,nzprocs).ne.0) then
  write(IOUT,120) nzp2-2,nzprocs
  nerror=1
  return
end if
100  format(/,'ERROR in mgdinit: nx=',i3,' is not a multiple of ', &
     &       'nxprocs=',i3,/,'cannot use the new version of the ', &
     &       'multigrid code',/)
110  format(/,'ERROR in mgdinit: ny=',i3,' is not a multiple of ', &
     &       'nyprocs=',i3,/,'cannot use the new version of the ', &
     &       'multigrid code',/)
120  format(/,'ERROR in mgdinit: nz=',i3,' is not a multiple of ', &
     &       'nzprocs=',i3,/,'cannot use the new version of the ', &
     &       'multigrid code',/)
# else

!
! check that the old version is not used with non-periodic BCs
!
if (ibdry.ne.0.or.jbdry.ne.0.or.kbdry.ne.0) then
  write(IOUT,130) ibdry,jbdry,kbdry
  nerror=1
  return
end if
130  format(/,'ERROR in mgdinit: ibdry=',i2,' jbdry=',i2,' kbdry=',i2, &
     &       /,'cannot use the old version of the multigrid code',     &
     &       /,'boundary conditions that are not periodic',            &
     &       /,'-> change compiler directive to 1 in compdir.inc',     &
     &       /,'   and recompile the multigrid code',/)
# endif

!
! check that the dimensions are correct
!
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
140  format(/,'ERROR in mgdinit: nxp1=',i3,' <> ixp*2**(iex-1)+1=', &
     &       i3,/,'-> adjust the multigrid parameters ixp and iex', &
     &       ' in main',/)
150  format(/,'ERROR in mgdinit: nyp1=',i3,' <> jyq*2**(jey-1)+1=', &
     &       i3,/,'-> adjust the multigrid parameters jyq and jey', &
     &       ' in main',/)
160  format(/,'ERROR in mgdinit: nzp1=',i3,' <> kzr*2**(kez-1)+1=', &
     &       i3,/,'-> adjust the multigrid parameters kzr and kez', &
     &       ' in main',/)
!
! check that the number of points at the coarser level is not smaller
! than the number of processes in either direction
!
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
170  format(/,'ERROR in mgdinit: ixp=',i3,' < nxprocs=',i3,/,    &
     &       ' there must be at least one grid point at the ',   &
     &       'coarsest grid level',/,                            &
     &       '-> increase ixp and decrease iex correspondingly', &
     &       ' in main',/)
180  format(/,'ERROR in mgdinit: jyq=',i3,' < nyprocs=',i3,/,    &
     &       ' there must be at least one grid point at the ',   &
     &       'coarsest grid level',/,                            &
     &       '-> increase jyq and decrease jey correspondingly', &
     &       ' in main',/)
190  format(/,'ERROR in mgdinit: kzr=',i3,' < nzprocs=',i3,/,    &
     &       ' there must be at least one grid point at the ',   &
     &       'coarsest grid level',/,                            &
     &       '-> increase kzr and decrease kez correspondingly', &
     &       ' in main',/)
!
! check that coarsifying takes place in all directions at the finest
! grid level
!
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
200  format(/,'ERROR in mgdinit: ngrid=',i3,' iex=',i3,                  &
     &       /,'no coarsifying at the finest grid level in x-direction', &
     &       /,'this is not allowed by the mutligrid code',/)
210  format(/,'ERROR in mgdinit: ngrid=',i3,' jey=',i3,                  &
     &       /,'no coarsifying at the finest grid level in y-direction', &
     &       /,'this is not allowed by the mutligrid code',/)
220  format(/,'ERROR in mgdinit: ngrid=',i3,' kez=',i3,                  &
     &       /,'no coarsifying at the finest grid level in z-direction', &
     &       /,'this is not allowed by the mutligrid code',/)
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
  nzk(k)=kzr*2**(max(k+kez-ngrid,1)-1)+1
end do
!
! for all grid levels, set the indices of the subdomain the process
! 'myid' work on, as well as the datatypes needed for the exchange
! of boundary data:
!
nxf=nxk(ngrid)
nyf=nyk(ngrid)
nzf=nzk(ngrid)
sxk(ngrid)=sx
exk(ngrid)=ex
syk(ngrid)=sy
eyk(ngrid)=ey
szk(ngrid)=sz
ezk(ngrid)=ez
call grid1_type(kdatatype(1,ngrid),realtype,sxk(ngrid),exk(ngrid), &
                syk(ngrid),eyk(ngrid),szk(ngrid),ezk(ngrid),IOUT)
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
  call grid1_type(kdatatype(1,k),realtype,sxk(k),exk(k), &
                  syk(k),eyk(k),szk(k),ezk(k),IOUT)
end do
!
! set work space indices for phi, cof at each grid level, and check
! that there is sufficient work space
!
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
230 format(/,'ERROR in mgdinit: not enough work space',/,            &
    &        ' kps=',i12,' nwork=',i12,' myid: ',i3,/,               &
    &        ' -> put the formula for nwork in main in ',            &
    &        'comments',/,'    and set nwork to the value of kps',/)
else
  write(IOUT,240) kps,nwork
240  format(/,'WARNING in mgdinit: kps=',i10,' nwork=',i10,          &
     &         /,'can optimize the amount of memory needed by ',     &
     &           'the multigrid code',/,'by putting the formula ',   &
     &           'for nwork into comments and setting',/,'nwork ',   &
     &           'to the value of kps',/)
end if
# if WMGD
!------------------------------------------------------------------------
! For the new version of the multigrid code, set the boundary values 
! to be used for the Dirichlet boundaries. It is possible to assign
! 6 different constant values to the 6 different sides. The values are
! assigned at the finest grid level, zero is assigned at all levels
! below
! 
! vbc, phibc:
!                                      k
!                   vbc(6)             
!                   /                  ^  ^
!       -----vbc(4)/-----              | / i
!       |         /     |              |/
!       |        /      |         -----/-----> j
!     vbc(3)----+-----vbc(1)          /|
!       |      /        |            / |
!       |     /         |           /  |
!       -----/vbc(2)----|
!           /
!         vbc(5)
!
do j=1,6
  phibc(j,ngrid)=vbc(j)
end do
do k=ngrid-1,1,-1
  do j=1,6
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
!------------------------------------------------------------------------
! set indices for determining the coefficients in the elliptic
! equation div(cof*grad(P))=rhs. Used only when solving for the
! pressure. When setting these coefficients at level k, need
! the values of the density at midpoints, i.e. at level k+1
! (if coarsifying takes place between the levels k and k+1).
! If coarsifying took place at all levels, the array cof could
! be used as temporary storage space for the densities, with
! cof(*,*,8) at level k+1 giving the values cof(*,*,1->7) at level
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
 250  format('ERROR in mgdinit: no coarsifying between level ', &
      &      i3,' and level ',i3,/,', the current version of ', &
      &      'the code cannot cope with that',/,                &
      &      ' -> change the parameters (ixp,jyq,kzr) and',     &
      &      ' (iex,jyq,kez) in main',/)
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
  call grid1_type(rdatatype(1,k),realtype,sxr(k),exr(k), &
                  syr(k),eyr(k),szr(k),ezr(k),IOUT)
end do
# endif

# if cdebug
timing(81)=timing(81)+MPI_WTIME()-tinitial
# endif
return
end

