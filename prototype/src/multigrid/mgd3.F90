module mgd3
#include "sll_working_precision.h"
use mpi 
implicit none

#include "mgd3.h"

   sll_int32 , dimension(20)   :: nxk,nyk,nzk,sxk,exk,syk,eyk,szk,ezk
   sll_int32 , dimension(20)   :: kpbgn,kcbgn
   sll_int32 , dimension(7,20) :: kdatatype
   sll_int32 , dimension(20)   :: sxi,exi,syi,eyi,szi,ezi
   sll_int32 , dimension(20)   :: nxr,nyr,nzr,sxr,exr,syr,eyr,szr,ezr
   sll_int32 , dimension(7,20) :: rdatatype
   sll_int32 , parameter       :: iout=6

   type, public :: block
      sll_int32  :: id
      sll_int32  :: sx, ex
      sll_int32  :: sy, ey
      sll_int32  :: sz, ez
      sll_int32  :: ixp, jyq, kzr
      sll_int32  :: iex, jey, kez
      sll_int32  :: ngrid
      sll_int32  :: neighbor(26),bd(26)
   end type block

   type, public :: mg_solver
      sll_int32  :: nx, ny, nz
      sll_int32  :: ibdry
      sll_int32  :: jbdry
      sll_int32  :: kbdry
      sll_int32  :: nxprocs
      sll_int32  :: nyprocs
      sll_int32  :: nzprocs
      sll_real64 :: vbc(6),phibc(6,20)
      sll_int32  :: comm3d,comm3dp,comm3dl,comm3dc
      sll_real64 :: tolmax
      sll_real64 :: xl,yl,zl
      sll_int32  :: maxcy, kcycle, iprer, ipost, iresw
      sll_int32  :: isol
   end type mg_solver

   sll_real64, private, allocatable :: work(:)

   private grid1_type
   sll_real64, parameter, private :: rro=1.0d0

contains

subroutine initialize_mgd3(my_block,my_mg,nerror)

implicit none

integer         :: ngrid
integer         :: ixp,jyq,kzr,iex,jey,kez,nxp2,nyp2,nzp2
integer         :: sx,ex,sy,ey,sz,ez
integer         :: nerror,m0,m1
type(block)     :: my_block
type(mg_solver) :: my_mg
 
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
character(len=1)  :: num(10)
data (num(i),i=1,10)/'0','1','2','3','4','5','6','7','8','9'/


!-----------------------------------------------------------------------
! open file for output of messages and check that the number of 
! processes is correct
!
m1=mod(my_block%id,10)+1
m0=mod(my_block%id/10,10)+1

!------------------------------------------------------------------------
! set /mgd/ variables to zero
!

nxk=0
nyk=0
nzk=0
sxk=0
exk=0
syk=0
eyk=0
szk=0
ezk=0
kpbgn=0
kcbgn=0

do k=1,20
  do j=1,7
    kdatatype(j,k)=MPI_DATATYPE_NULL
    rdatatype(j,k)=MPI_DATATYPE_NULL
  end do
end do

sxi=0
exi=0
syi=0
eyi=0
szi=0
ezi=0
nxr=0
nyr=0
nzr=0
sxr=0
exr=0
syr=0
eyr=0
szr=0
ezr=0

sx = my_block%sx ; ex = my_block%ex
sy = my_block%sy ; ey = my_block%ey
sz = my_block%sz ; ez = my_block%ez

!
! check that the dimensions are correct
!

ixp = my_block%ixp 
jyq = my_block%jyq 
kzr = my_block%kzr 
iex = my_block%iex 
jey = my_block%jey 
kez = my_block%kez 

ngrid = my_block%ngrid

nxp2 = my_mg%nx+2
nyp2 = my_mg%ny+2
nzp2 = my_mg%nz+2

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
if (mod(nxp2-2,my_mg%nxprocs).ne.0) then
  write(IOUT,100) nxp2-2,my_mg%nxprocs
  nerror=1
  return
end if
if (mod(nyp2-2,my_mg%nyprocs).ne.0) then
  write(IOUT,110) nyp2-2,my_mg%nyprocs
  nerror=1
  return
end if
if (mod(nzp2-2,my_mg%nzprocs).ne.0) then
  write(IOUT,120) nzp2-2,my_mg%nzprocs
  nerror=1
  return
end if
# else

!
! check that the old version is not used with non-periodic BCs
!
if (my_mg%ibdry.ne.0.or.my_mg%jbdry.ne.0.or.my_mg%kbdry.ne.0) then
  write(IOUT,130) my_mg%ibdry,my_mg%jbdry,my_mg%kbdry
  nerror=1
  return
end if
# endif


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

!
! check that the number of points at the coarser level is not smaller
! than the number of processes in either direction
!
if (ixp.lt.my_mg%nxprocs) then
  write(IOUT,170) ixp,my_mg%nxprocs
  nerror=1
  return
end if
if (jyq.lt.my_mg%nyprocs) then
  write(IOUT,180) jyq,my_mg%nyprocs
  nerror=1
  return
end if
if (kzr.lt.my_mg%nzprocs) then
  write(IOUT,190) kzr,my_mg%nzprocs
  nerror=1
  return
end if
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

call grid1_type(kdatatype(1,ngrid),sxk(ngrid),exk(ngrid), &
                syk(ngrid),eyk(ngrid),szk(ngrid),ezk(ngrid))

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
  call grid1_type(kdatatype(1,k),sxk(k),exk(k), &
                  syk(k),eyk(k),szk(k),ezk(k))
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

!if (kps.gt.nwork) then
!  write(IOUT,230) kps,nwork,my_block%id
!  nerror=1
!  return
!else
!  write(IOUT,240) kps,nwork
!end if

allocate(work(kps))

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
  my_mg%phibc(j,ngrid)=my_mg%vbc(j)
end do
do k=ngrid-1,1,-1
  do j=1,6
    my_mg%phibc(j,k)=0.0d0
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
  call grid1_type(rdatatype(1,k),sxr(k),exr(k), &
                  syr(k),eyr(k),szr(k),ezr(k))
end do
# endif

return

#if WMGD
100 format(/,'ERROR in mgdinit: nx=',i3,' is not a multiple of ', &
    &       'nxprocs=',i3,/,'cannot use the new version of the ', &
    &       'multigrid code',/)
110 format(/,'ERROR in mgdinit: ny=',i3,' is not a multiple of ', &
    &       'nyprocs=',i3,/,'cannot use the new version of the ', &
    &       'multigrid code',/)
120 format(/,'ERROR in mgdinit: nz=',i3,' is not a multiple of ', &
    &       'nzprocs=',i3,/,'cannot use the new version of the ', &
    &       'multigrid code',/)
#endif
130 format(/,'ERROR in mgdinit: ibdry=',i2,' jbdry=',i2,' kbdry=',i2, &
    &       /,'cannot use the old version of the multigrid code',     &
    &       /,'boundary conditions that are not periodic',            &
    &       /,'-> change compiler directive to 1 in compdir.inc',     &
    &       /,'   and recompile the multigrid code',/)
140 format(/,'ERROR in mgdinit: nxp1=',i3,' <> ixp*2**(iex-1)+1=', &
    &       i3,/,'-> adjust the multigrid parameters ixp and iex', &
    &       ' in main',/)
150 format(/,'ERROR in mgdinit: nyp1=',i3,' <> jyq*2**(jey-1)+1=', &
    &       i3,/,'-> adjust the multigrid parameters jyq and jey', &
    &       ' in main',/)
160 format(/,'ERROR in mgdinit: nzp1=',i3,' <> kzr*2**(kez-1)+1=', &
    &       i3,/,'-> adjust the multigrid parameters kzr and kez', &
    &       ' in main',/)   
170 format(/,'ERROR in mgdinit: ixp=',i3,' < nxprocs=',i3,/,    &
    &       ' there must be at least one grid point at the ',   &
    &       'coarsest grid level',/,                            &
    &       '-> increase ixp and decrease iex correspondingly', &
    &       ' in main',/)
180 format(/,'ERROR in mgdinit: jyq=',i3,' < nyprocs=',i3,/,    &
    &       ' there must be at least one grid point at the ',   &
    &       'coarsest grid level',/,                            &
    &       '-> increase jyq and decrease jey correspondingly', &
    &       ' in main',/)
190 format(/,'ERROR in mgdinit: kzr=',i3,' < nzprocs=',i3,/,    &
    &       ' there must be at least one grid point at the ',   &
    &       'coarsest grid level',/,                            &
    &       '-> increase kzr and decrease kez correspondingly', &
    &       ' in main',/)
200 format(/,'ERROR in mgdinit: ngrid=',i3,' iex=',i3,                  &
    &       /,'no coarsifying at the finest grid level in x-direction', &
    &       /,'this is not allowed by the mutligrid code',/)
210 format(/,'ERROR in mgdinit: ngrid=',i3,' jey=',i3,                  &
    &       /,'no coarsifying at the finest grid level in y-direction', &
    &       /,'this is not allowed by the mutligrid code',/)
220 format(/,'ERROR in mgdinit: ngrid=',i3,' kez=',i3,                  &
    &       /,'no coarsifying at the finest grid level in z-direction', &
    &       /,'this is not allowed by the mutligrid code',/)
!230 format(/,'ERROR in mgdinit: not enough work space',/,            &
!    &        ' kps=',i12,' nwork=',i12,' myid: ',i3,/,               &
!    &        ' -> put the formula for nwork in main in ',            &
!    &        'comments',/,'    and set nwork to the value of kps',/)
!240 format(/,'WARNING in mgdinit: kps=',i10,' nwork=',i10,          &
!    &         /,'can optimize the amount of memory needed by ',     &
!    &           'the multigrid code',/,'by putting the formula ',   &
!    &           'for nwork into comments and setting',/,'nwork ',   &
!    &           'to the value of kps',/)
250 format('ERROR in mgdinit: no coarsifying between level ', &
    &      i3,' and level ',i3,/,', the current version of ', &
    &      'the code cannot cope with that',/,                &
    &      ' -> change the parameters (ixp,jyq,kzr) and',     &
    &      ' (iex,jyq,kez) in main',/)

end subroutine initialize_mgd3


subroutine grid1_type(gtype,sx,ex,sy,ey,sz,ez)

implicit none 
integer :: gtype(7),realtype,sx,ex,sy,ey,sz,ez
!------------------------------------------------------------------------
! Define the 7 derived datatypes needed to communicate the boundary 
! data of (sx-1:ex+1,sy-1:ey+1,sz-1:ez+1) arrays between 'myid' and
! its 26 neighbors.
!
! gtype(l): l=1 -> i=const planes
!           l=2 -> j=const planes
!           l=3 -> k=const planes
!           l=4 -> (i=const,j=const) lines
!           l=5 -> (i=const,k=const) lines
!           l=6 -> (j=const,k=const) lines
!           l=7 -> (i=const,j=const,k=const) corners
!
! Code      : tmgd3, test program for 3D parallel multigrid solver
! Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
! Called in : mgdinit
! Calls     : MPI_TYPE_CONTIGUOUS, MPI_TYPE_COMMIT, MPI_TYPE_VECTOR,
!             MPI_TYPE_EXTENT, MPI_TYPE_HVECTOR
!------------------------------------------------------------------------
integer :: ierr
INTEGER(KIND=MPI_ADDRESS_KIND) :: i, LB

realtype = MPI_REAL8

!------------------------------------------------------------------------
! datatype for one 1*1 corner (i=const,j=const,k=const)
!
call MPI_TYPE_CONTIGUOUS(1,realtype,gtype(7),ierr)
call MPI_TYPE_COMMIT(gtype(7),ierr)
!------------------------------------------------------------------------
! datatype for one (i=const,j=const) line 
!
call MPI_TYPE_VECTOR(ez-sz+1,1,(ex-sx+3)*(ey-sy+3),realtype,gtype(4),ierr)
call MPI_TYPE_COMMIT(gtype(4),ierr)
!
! datatype for one (i=const,k=const) line
!
call MPI_TYPE_VECTOR(ey-sy+1,1,ex-sx+3,realtype,gtype(5),ierr)
call MPI_TYPE_COMMIT(gtype(5),ierr)
!
! datatype for one (j=const,k=const) line
!
call MPI_TYPE_CONTIGUOUS(ex-sx+1,realtype,gtype(6),ierr)
call MPI_TYPE_COMMIT(gtype(6),ierr)
!------------------------------------------------------------------------
! datatype for one i=const plane
!
!MPI_TYPE_EXTENT(DATATYPE, EXTENT, IERROR)
!    INTEGER    DATATYPE, EXTENT, IERROR
!call MPI_TYPE_EXTENT(realtype,i,ierr)

!MPI_TYPE_GET_EXTENT(DATATYPE, LB, EXTENT, IERROR)
!    INTEGER    DATATYPE, IERROR
!    INTEGER(KIND=MPI_ADDRESS_KIND) LB, EXTENT

call MPI_TYPE_GET_EXTENT(realtype, LB, i, IERR)

!MPI_TYPE_HVECTOR(COUNT, BLOCKLENGTH, STRIDE, OLDTYPE, NEWTYPE,
!        IERROR)
!    INTEGER    COUNT, BLOCKLENGTH, STRIDE, OLDTYPE
!    INTEGER    NEWTYPE, IERROR
!call MPI_TYPE_HVECTOR(ez-sz+1,1,(ex-sx+3)*(ey-sy+3)*i,gtype(5),gtype(1),ierr)

!MPI_TYPE_CREATE_HVECTOR(COUNT, BLOCKLENGTH, STRIDE, OLDTYPE,
!    NEWTYPE, IERROR)
!    INTEGER    COUNT, BLOCKLENGTH, OLDTYPE, NEWTYPE, IERROR
!    INTEGER(KIND=MPI_ADDRESS_KIND) STRIDE

call MPI_TYPE_CREATE_HVECTOR(ez-sz+1,1,(ex-sx+3)*(ey-sy+3)*i,gtype(5),gtype(1),ierr)

call MPI_TYPE_COMMIT(gtype(1),ierr)
!
! datatype for one j=const plane
!
call MPI_TYPE_VECTOR(ez-sz+1,ex-sx+1,(ex-sx+3)*(ey-sy+3),realtype,gtype(2),ierr)
call MPI_TYPE_COMMIT(gtype(2),ierr)
!
! datatype for one k=const plane
!
call MPI_TYPE_VECTOR(ey-sy+1,ex-sx+1,ex-sx+3,realtype,gtype(3),ierr)
call MPI_TYPE_COMMIT(gtype(3),ierr)

end subroutine grid1_type

subroutine solve(phif,rhsf,r,my_block,my_mg,nerror)

implicit none
integer :: sx,ex,sy,ey,sz,ez,ngrid
integer :: maxcy,kcycle,iprer,ipost,iresw
real(8) :: phif(:,:,:)
real(8) :: rhsf(:,:,:)
real(8) :: r(:,:,:)

real(8) :: xl,yl,zl,phibc(6,20)
integer :: isol,comm3d,comm3dp,comm3dl,comm3dc,myid
integer :: neighbor(26),bd(26),iter,nerror
integer :: k
integer :: nx, ny, nz, icf
real(8) :: relmax
type(block) :: my_block
type(mg_solver) :: my_mg


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

comm3d   = my_mg%comm3d
comm3dp  = my_mg%comm3dp
comm3dl  = my_mg%comm3dl
comm3dc  = my_mg%comm3dc
phibc    = my_mg%phibc
maxcy    = my_mg%maxcy
kcycle   = my_mg%kcycle
iprer    = my_mg%iprer
ipost    = my_mg%ipost
iresw    = my_mg%iresw
kcycle   = my_mg%kcycle
nx       = my_mg%nx
ny       = my_mg%ny
nz       = my_mg%nz

xl       = my_mg%xl
yl       = my_mg%yl
zl       = my_mg%zl

sx       = my_block%sx
sy       = my_block%sy
sz       = my_block%sz

ex       = my_block%ex
ey       = my_block%ey
ez       = my_block%ez

ngrid    = my_block%ngrid


myid     = my_block%id

neighbor = my_block%neighbor
bd       = my_block%bd

isol     = my_mg%isol

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
               comm3dp,comm3dl,comm3dc,neighbor,bd,phibc)
  sxm=sxk(ngrid)
  exm=exk(ngrid)
  sym=syk(ngrid)
  eym=eyk(ngrid)
  szm=szk(ngrid)
  ezm=ezk(ngrid)
  ip=kpbgn(ngrid)
  call mgderr(relmax,sxm,exm,sym,eym,szm,ezm,phif,work(ip),comm3d,IOUT)
  if (relmax.le.my_mg%tolmax) goto 1000
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
call gscale(sx,ex,sy,ey,sz,ez,phif,avo,acorr,comm3d,nx,ny,nz)
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
   if (myid.eq.0) write(IOUT,110) relmax,iter,acorr
   110 format('  R MGD     err=',e11.3,' iters=',i5,' rcorr=',e11.3)
else
  if (myid.eq.0) write(IOUT,120) relmax,iter,acorr
   120 format('  P MGD     err=',e11.3,' iters=',i5,' pcorr=',e11.3)
end if

return
end subroutine solve


!------------------------------------------------------------------------
!>  From the MPE library
!>  This file contains a routine for producing a decomposition of a 1-d 
!>  array when given a number of processors.  It may be used in "direct" 
!>  product decomposition.  The values returned assume a "global" domain 
!>  in [1:n]
subroutine MPE_DECOMP1D(n,numprocs,myid,s,e)
implicit none 

integer :: n, numprocs, myid, s, e
integer :: nlocal
integer :: deficit
!------------------------------------------------------------------------
nlocal  = n / numprocs
s = myid * nlocal + 1
deficit = mod(n,numprocs)
s = s + min(myid,deficit)
if (myid .lt. deficit) then
   nlocal = nlocal + 1
endif
e = s + nlocal - 1
if (e .gt. n .or. myid .eq. numprocs-1) e = n
return

end subroutine MPE_DECOMP1D

end module mgd3
