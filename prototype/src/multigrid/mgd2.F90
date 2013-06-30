!> Multigrid solver in 2 dimensions on a cartesian grid
!> \authors Pierre Navaro
!> \remark this a refactored version of Bernard Bunner code 
!> (http://cs-www.cs.yale.edu/homes/douglas-craig/mgnet-codes-bunner.html)
module mgd2
#include "sll_working_precision.h"
use mpi
use gxch1_2d

implicit none

sll_int32, private :: sx, ex, sy, ey

logical ::  WMGD 
sll_int32 :: nxk(20),nyk(20),sxk(20),exk(20),syk(20),eyk(20)
sll_int32 :: kpbgn(20),kcbgn(20),ikdatatype(20),jkdatatype(20)
sll_int32 :: ijkdatatype(20),sxi(20),exi(20),syi(20),eyi(20)
sll_int32 :: nxr(20),nyr(20),sxr(20),exr(20),syr(20),eyr(20)
sll_int32 :: irdatatype(20),jrdatatype(20),ijrdatatype(20)

type, public :: block
   sll_int32  :: id
   sll_int32  :: sx, ex
   sll_int32  :: sy, ey
   sll_int32  :: ixp, jyq
   sll_int32  :: iex, jey
   sll_int32  :: ngrid
   sll_int32  :: neighbor(8),bd(8)
end type block

type, public :: mg_solver
   sll_int32  :: nx, ny
   sll_int32  :: ibdry
   sll_int32  :: jbdry
   sll_int32  :: nxprocs
   sll_int32  :: nyprocs
   sll_real64 :: vbc(4),phibc(4,20)
   sll_int32  :: comm2d
   sll_real64 :: tolmax
   sll_real64 :: xl,yl
   sll_int32  :: maxcy
   sll_int32  :: kcycle
   sll_int32  :: iprer
   sll_int32  :: ipost
   sll_int32  :: iresw
   sll_int32  :: isol
end type mg_solver

contains

!>  \function
!>  From the MPE library \n
!>  This file contains a routine for producing a decomposition of a 1-d  \n
!>  array when given a number of processors.  It may be used in "direct"  \n
!>  product decomposition.  The values returned assume a "global" domain  \n
!>  in [1:n]
subroutine MPE_DECOMP1D(n,numprocs,myid,s,e)

   sll_int32  :: n
   sll_int32  :: numprocs
   sll_int32  :: myid
   sll_int32  :: s
   sll_int32  :: e
   sll_int32  :: nlocal
   sll_int32  :: deficit

   nlocal  = n / numprocs
   s = myid * nlocal + 1
   deficit = mod(n,numprocs)
   s = s + min(myid,deficit)
   if (myid < deficit) then
      nlocal = nlocal + 1
   endif
   e = s + nlocal - 1
   if (e > n .or. myid == numprocs-1) e = n
   return

end subroutine

!> Initialize the parallel multigrid solver: subdomain indices and
!> MPI datatypes.
!>
!> The multigrid code comes in two versions. With the WMGD compiler
!> directive set to 0, the grid setup is vertex-centered:
!>
!> WMGD=0
!> 
!>  |------|-----|-----|-----|-----|            fine
!>  1      2     3     4     5     6 
!>
!>  |------------|-----------|-----------|      coarse
!>  1            2           3           4
!>
!> With WMGD set to 1, it is cell-centered:
!>
!> WMGD=1   
!>           |                       |
!>        |--|--|-----|-----|-----|--|--|       fine
!>        1  |  2     3     4     5  |  6
!>           |                       |
!>     |-----|-----|-----------|-----|-----|    coarse
!>     1     |     2           3     |     4
!>           |                       |
!>          wall                    wall
!>
!> For WMGD=0, the restriction and correction operators are standard
!> (choice of full or half weighting for the restriction, bilinear
!> interpolation for the correction). This works fine for periodic
!> boundary conditions. However, when there are Neumann (wall) or
!> Dirichlet BCs, this grid setup results in a loss of accuracy near
!> the boundaries when the grid is staggered (the discretization of
!> the relaxation operator is first-order locally there). With the
!> grid setup corresponding to WMGD=1, accuracy remains second-order
!> all the time. As the grid gets coarser, it remains centered on the
!> domain instead of "shifting to the right". This option works for
!> periodic, Neumann, and Dirichlet BCs, although only periodic and
!> Neumann BCs have been tested thoroughly. There is one catch, though.
!> For a problem with purely periodic BCs, WMGD=0 converges in less
!> cycles than WMGD=1 and requires less CPU time (the penalty is
!> apparently between 10 and 50%). This can be attributed to the loss
!> of accuracy in the restriction and correction operators due to the
!> fact that WMGD=0 uses a support of 3 points in each direction 
!> whereas WMGD=1 uses only 2 points.
!>
!> Both versions offer the option to coarsify in one direction and
!> not the other, except at the finest grid level, where coarsifying
!> MUST take place along all axes. However, it is possible to have
!> ngrid=iex=jey=1, with ixp=nx and jyq=ny. In this case, the code
!> never enters 'mgdrestr' and 'mgdcor' and all it does is Gauss-Seidel
!> iterate at the finest grid level. This can be useful as a preliminary
!> check.
!>
!> Note: some memory could be saved by noting that the cof arrays
!> need be dimensioned (sxm:exm,sym:eym) and not
!> (sxm-1:exm+1,sym-1:eym+1)... Probably not too difficult to
!> make the change
!>
!> Bernard Bunner (bunner@engin.umich.edu), January 1998
subroutine initialize_mgd2(my_block, my_mg, nerror)

sll_int32       :: nerror
type(block)     :: my_block
type(mg_solver) :: my_mg

!sll_int32  :: nxk,nyk,sxk,exk,syk,eyk,kpbgn,kcbgn
!sll_int32  :: ikdatatype,jkdatatype,ijkdatatype
!sll_int32  :: sxi,exi,syi,eyi
!sll_int32  :: nxr,nyr,sxr,exr,syr,eyr
!sll_int32  :: irdatatype,jrdatatype,ijrdatatype
sll_int32  :: i,j,k,nxf,nyf,nxm,nym,kps,sxm,exm,sym,eym,ierr,nxc,nyc

! set /mgd/ variables to zero
nerror = 0
do k=1,20
  nxk(k)         = 0
  nyk(k)         = 0
  sxk(k)         = 0
  exk(k)         = 0
  syk(k)         = 0
  eyk(k)         = 0
  kpbgn(k)       = 0
  kcbgn(k)       = 0
  ikdatatype(k)  = MPI_DATATYPE_NULL
  jkdatatype(k)  = MPI_DATATYPE_NULL
  ijkdatatype(k) = MPI_DATATYPE_NULL
  irdatatype(k)  = MPI_DATATYPE_NULL
  jrdatatype(k)  = MPI_DATATYPE_NULL
  ijrdatatype(k) = MPI_DATATYPE_NULL
  sxi(k)         = 0
  exi(k)         = 0
  syi(k)         = 0
  eyi(k)         = 0
  nxr(k)         = 0
  nyr(k)         = 0
  sxr(k)         = 0
  exr(k)         = 0
  syr(k)         = 0
  eyr(k)         = 0
end do


if (my_mg%ibdry.ne.0.or.my_mg%jbdry.ne.0) then

  write(6,120) my_mg%ibdry,my_mg%jbdry
  WMGD = .true.
  ! check that the number of processes in each
  ! direction divides the number of points in that direction 
  if (mod(my_mg%nx,my_mg%nxprocs).ne.0) then
    write(6,100) my_mg%nx,my_mg%nxprocs
    nerror=1
    return
  end if
  if (mod(my_mg%ny,my_mg%nyprocs).ne.0) then
    write(6,110) my_mg%ny,my_mg%nyprocs
    nerror=1
    return
  end if

else

  WMGD = .false.

endif

100 format(/,'ERROR in mgdinit: nx=',i3,' is not a multiple of ', &
           'nxprocs=',i3,/,'cannot use the new version of the ',  &
           'multigrid code',/)
110 format(/,'ERROR in mgdinit: ny=',i3,' is not a multiple of ', &
           'nyprocs=',i3,/,'cannot use the new version of the ',  &
           'multigrid code',/)
120  format(/,'NOTE in mgdinit: ibdry=',i2,' jbdry=',i2, &
            /,'boundary conditions that are not periodic',/)
!------------------------------------------------------------------------
! define all grid levels
! I have adopted the same notations as in Mudpack as far as possible.
! When a confusion was possible, I added a suffix 'm' to the name
! of the variables. For example, nxm is nx+1 for the multigrid
! code whereas nx means nxp2-2 in the rest of the code.
!
do k=1,my_block%ngrid
  nxk(k)=my_block%ixp*2**(max(k+my_block%iex-my_block%ngrid,1)-1)+1
  nyk(k)=my_block%jyq*2**(max(k+my_block%jey-my_block%ngrid,1)-1)+1
end do

! for all grid levels, set the indices of the subdomain the process
! 'myid' work on, as well as the datatypes needed for the exchange
! of boundary data

nxf=nxk(my_block%ngrid)
nyf=nyk(my_block%ngrid)
sxk(my_block%ngrid)=my_block%sx
exk(my_block%ngrid)=my_block%ex
syk(my_block%ngrid)=my_block%sy
eyk(my_block%ngrid)=my_block%ey
call grid1_type(my_block,ikdatatype(my_block%ngrid), &
                         jkdatatype(my_block%ngrid), &
                         ijkdatatype(my_block%ngrid),&
                         my_block%sx, &
                         my_block%ex, &
                         my_block%sy, &
                         my_block%ey)
do k=my_block%ngrid-1,1,-1
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
  call grid1_type(my_block,ikdatatype(k),jkdatatype(k),ijkdatatype(k), &
                  sxk(k),exk(k),syk(k),eyk(k))
end do
!
! set work space indices for phi, cof at each grid level
! check that there is sufficient work space
!
kps=1
do k=my_block%ngrid,1,-1
  sxm=sxk(k)
  exm=exk(k)
  sym=syk(k)
  eym=eyk(k)
  kpbgn(k)=kps
  kcbgn(k)=kpbgn(k)+(exm-sxm+3)*(eym-sym+3)
  kps=kcbgn(k)+6*(exm-sxm+3)*(eym-sym+3)
end do

!if (kps.gt.nwork) then
!  write(6,200) kps,nwork,my_block%id
!  nerror=1
!  return
!else
!  write(6,210) kps,nwork
!end if

!200 format(/,'ERROR in mgdinit: not enough work space',/, &
!    ' kps=',i10,' nwork=',i10,' myid: ',i3,/, &
!    '-> put the formula for nwork in main in ', &
!    'comments',/,'   and set nwork to the value of kps',/)
!210 format(/,'WARNING in mgdinit: kps=',i10,' nwork=',i10, &
!         /,'can optimize the amount of memory needed by ', &
!           'the multigrid code',/,'by putting the formula ', &
!           'for nwork into comments and setting',/,'nwork ', &
!           'to the value of kps',/)

if (WMGD) then
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
   my_mg%phibc(:,:)=0.0d0
   my_mg%phibc(:,my_block%ngrid)=my_mg%vbc(:)
else  
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
      
sxi(my_block%ngrid)=sxk(my_block%ngrid)-1
exi(my_block%ngrid)=exk(my_block%ngrid)+1
syi(my_block%ngrid)=syk(my_block%ngrid)-1
eyi(my_block%ngrid)=eyk(my_block%ngrid)+1
nxf=nxk(my_block%ngrid)
nyf=nyk(my_block%ngrid)
do k=my_block%ngrid-1,1,-1
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
do k=my_block%ngrid-1,1,-1
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
    if (k.eq.(my_block%ngrid-1)) then
      write(6,300) my_block%ngrid,my_block%ngrid-1
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
    if (k.eq.(my_block%ngrid-1)) then
      write(6,300) my_block%ngrid,my_block%ngrid-1
      nerror=1
      return
    end if
    syr(k)=sym
    eyr(k)=eym
    nyr(k)=nym
  end if
  call grid1_type(my_block,irdatatype(k),jrdatatype(k),ijrdatatype(k), &
                  sxr(k),exr(k),syr(k),eyr(k))
end do
endif

return
300 format('ERROR in mgdinit: no coarsifying between level ', &
           i3,' and level ',i3,/,', the current version of ', &
           'the code cannot cope with that',/, &
           ' -> decrease the value of ixp and/or jyq and',/, &
           '    increase the value of iex and/or jey')
end subroutine

!> Parallel multigrid solver in 2-D cartesian coordinates for the
!> elliptic equation:      div(cof(grad(phif)))=rhsf
!>
!> Written for periodic, wall (Neumann), and constant value
!> (Dirichlet) BCs. Tested roughly for all these BCs. There are 
!> two versions of the multigrid code, which are separated by the 
!> compiler directive WMGD set in 'compdir.inc'. The old version 
!> (WMGD=0) corresponds to the original, "traditional" grid setup, 
!> and works well when all boundary conditions are periodic. When one 
!> of the BCs is not periodic, must compile with the new version 
!> (WMGD=1), which uses a different grid setup and new restriction 
!> and correction operators (see 'mgdinit' for more details). It is 
!> less accurate (or ,equivalently, slower to converge) for the case 
!> of all-periodic BCs than the old version, but is better able to 
!> handle wall BCs.
!>
!> Notes: - the values of the rhs contained in the array rhsf are
!>          transferred to the work vector and the memory of the 
!>          rhsf array is then used to store the residuals at the
!>          different levels; therefore, it does not preserve its
!>          initial values
!>        - with the initialization and clean-up routines mgdinit
!>          and mgdend, this multigrid code is self-standing
!>
subroutine mgd2_solver(my_block,my_mg,phif,rhsf,r,work, &
                     iter,nprscr,nerror)
sll_real64      :: phif(:,:),rhsf(:,:)
sll_real64      :: r(:,:)
sll_real64      :: work(*)
sll_int32       :: iter,nerror
type(block)     :: my_block
type(mg_solver) :: my_mg

logical    :: nprscr

sll_int32  :: k, icf
sll_real64 :: avo,acorr, relmax
sll_int32  :: sxf,exf,syf,eyf,nxf,nyf,sxc,exc,syc,eyc,nxc,nyc
sll_int32  :: ipf,irf,ipc,irc,sxm,exm,sym,eym,nxm,nym,ip,ic,kcur
sll_int32  :: itype,jtype,ijtype,lev,ir1,ir2


! discretize pde at all levels
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
  sxf=sxk(my_block%ngrid)
  exf=exk(my_block%ngrid)
  syf=syk(my_block%ngrid)
  eyf=eyk(my_block%ngrid)
  nxf=nxk(my_block%ngrid)
  nyf=nyk(my_block%ngrid)
  icf=kcbgn(my_block%ngrid)
  call mgdpfpde(sxf,exf,syf,eyf,nxf,nyf,work(icf),r,my_mg%xl,my_mg%yl,my_block%bd)
if (WMGD) then
!
! new version: determine coefficients at coarser grid levels from
! four neighboring density points; no communication of boundary data
! is involved because of the grid setup, under the condition that
! mod(nx,nxprocs)=mod(ny,nyprocs)=0
!
  do k=my_block%ngrid-1,1,-1
    sxm=sxk(k)
    exm=exk(k)
    sym=syk(k)
    eym=eyk(k)
    nxm=nxk(k)
    nym=nyk(k)
    ic=kcbgn(k)
    call mgdphpde(sxm,exm,sym,eym,nxm,nym,work(ic), &
                  sx,ex,sy,ey,nxf,nyf,my_block%bd,my_mg%xl,my_mg%yl)
  end do
else  
  if (my_block%ngrid.gt.1) then
!
! old version: use two locations ir1 and ir2 in the work vector to
! store the density on the different grid levels. First set r in 
! work(ir1) at the finest grid level; this is used to determine the 
! coefficients at level k=ngrid-1
!
    ir1=kpbgn(my_block%ngrid)
    ir2=kcbgn(my_block%ngrid)+5*(exf-sxf+3)*(eyf-syf+3)
    sxf=sxr(my_block%ngrid-1)
    exf=exr(my_block%ngrid-1)
    syf=syr(my_block%ngrid-1)
    eyf=eyr(my_block%ngrid-1)
    nxf=nxr(my_block%ngrid-1)
    nyf=nyr(my_block%ngrid-1)
    lev=1
    call mgdrsetf(sxf,exf,syf,eyf,work(ir1),r)
!
! for the levels k=ngrid-1,1, calculate the coefficients from the
! densities stored in the arrays work(ir1) and work(ir2)
! for the levels k=ngrid-1,2, transfer to the level below the values
! of the density needed there to determine the coefficients; exchange
! of the boundary density data is necessary
!
    do k=my_block%ngrid-1,1,-1
      sxm=sxk(k)
      exm=exk(k)
      sym=syk(k)
      eym=eyk(k)
      nxm=nxk(k)
      nym=nyk(k)
      ic=kcbgn(k)
      if (lev.eq.1) then
        call mgdppde(sxm,exm,sym,eym,nxm,nym,work(ic), &
                     sxf,exf,syf,eyf,work(ir1),        &
                     my_mg%xl,my_mg%yl,my_block%bd)
      else
        call mgdppde(sxm,exm,sym,eym,nxm,nym,work(ic), &
                     sxf,exf,syf,eyf,work(ir2),        &
                     my_mg%xl,my_mg%yl,my_block%bd)
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
                        my_mg%comm2d,my_block%id,          &
                        my_block%neighbor,my_block%bd,     &
                        itype,jtype)
          lev=2
        else
          call mgdrtrsf(sxc,exc,syc,eyc,nxc,nyc,work(ir1), &
                        sxf,exf,syf,eyf,nxf,nyf,work(ir2), &
                        my_mg%comm2d,my_block%id,          &
                        my_block%neighbor,my_block%bd,     &
                        itype,jtype)
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
endif

!------------------------------------------------------------------------
! set phi,rhsf in work at the finest grid level
!
sxf = sxk(my_block%ngrid)
exf = exk(my_block%ngrid)
syf = syk(my_block%ngrid)
eyf = eyk(my_block%ngrid)
ipf = kpbgn(my_block%ngrid)
irf = kcbgn(my_block%ngrid)+5*(exf-sxf+3)*(eyf-syf+3)
call mgdsetf(sxf,exf,syf,eyf,work(ipf),work(irf),phif,rhsf)
!------------------------------------------------------------------------
! cycling at kcur=ngrid level
!
kcur=my_block%ngrid
do iter=1,my_mg%maxcy
  call mgdkcyc(work,rhsf,kcur,my_mg%kcycle,my_mg%iprer,     &
               my_mg%ipost,my_mg%iresw, my_mg%comm2d,       &
               my_block%id,my_block%neighbor,my_block%bd,   &
               my_mg%phibc)
  sxm=sxk(my_block%ngrid)
  exm=exk(my_block%ngrid)
  sym=syk(my_block%ngrid)
  eym=eyk(my_block%ngrid)
  ip=kpbgn(my_block%ngrid)
  call mgderr(relmax,sxm,exm,sym,eym,phif,work(ip),my_mg%comm2d)
  if (relmax.le.my_mg%tolmax) goto 1000
end do
!------------------------------------------------------------------------
! if not converged in maxcy cycles, issue an error message and quit
!
if (my_block%id.eq.0) write(6,100) my_mg%maxcy,relmax
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
avo=0.0d0
call gscale(my_block%sx,my_block%ex,my_block%sy,my_block%ey, &
            phif,avo,acorr,my_mg%comm2d,my_mg%nx,my_mg%ny)
!
! exchange boundary data and impose periodic BCs
!
call gxch1lin(phif,my_mg%comm2d,my_block%sx,my_block%ex,     &
              my_block%sy,my_block%ey,my_block%neighbor,     &
              my_block%bd, ikdatatype(my_block%ngrid),       &
              jkdatatype(my_block%ngrid))

call gxch1cor(phif,my_mg%comm2d,sxm,exm,sym,eym,             &
              my_block%neighbor,my_block%bd,                 &
              ijkdatatype(my_block%ngrid))
if (WMGD) then
!
! impose wall and Dirichlet BCs
!
call mgdbdry(sx,ex,sy,ey,phif,my_block%bd,my_mg%phibc)
endif

if (nprscr.and.my_block%id.eq.0) write(6,120) relmax,iter,acorr
120 format('  P MGD     err=',e8.3,' iters=',i5,' pcorr=',e9.3)

return
end subroutine

!> Define the 3 derived datatypes needed to communicate the boundary
!> data of (sx-1:ex+1,sy-1:ey+1) arrays between 'myid' and its 8
!> neighbors
subroutine grid1_type(my_block,itype,jtype,ijtype,sx,ex,sy,ey)
sll_int32  :: itype,jtype,ijtype,realtype,sx,ex,sy,ey
sll_int32  :: ierr
type(block) :: my_block

realtype = MPI_REAL8
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

return
end subroutine


!------------------------------------------------------------------------
!> Determine coefficients for the pressure equation at the finest grid
!> level. These coefficients involve densities half-way between the
!> pressure and density nodes. Works for periodic, Neumann, and
!> Dirichlet boundary conditions.
!>
!> cof array:
!>
!>         cof(4)
!>           |
!>           |
!> cof(1)--cof(5)--cof(2)
!>           |
!>           |
!>         cof(3)
!>
! Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
!------------------------------------------------------------------------
subroutine mgdpfpde(sxf,exf,syf,eyf,nxf,nyf,cof,r,xl,yl,bd)

sll_int32  :: sxf,exf,syf,eyf,nxf,nyf,bd(8)
sll_real64 :: cof(sxf-1:exf+1,syf-1:eyf+1,6)
sll_real64 :: r(sxf-1:exf+1,syf-1:eyf+1),xl,yl
sll_real64 :: dlx,todlxx,dly,todlyy,rij
sll_int32  :: i,j

dlx=xl/float(nxf-1)
dly=yl/float(nyf-1)

   todlxx=2.0d0/(dlx*dlx)
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

! enforce wall BCs 

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

! calculate diagonal term

do j=syf,eyf
  do i=sxf,exf
    cof(i,j,5)=-(cof(i,j,1)+cof(i,j,2)+cof(i,j,3)+cof(i,j,4))
  end do
end do

end subroutine


!> For the new version of the multigrid code, determine the coefficients
!> for the pressure equation at all grid levels except the finest one.
!> The coefficients are determined directly from the density array r
!> through some manipulation of indices and are values at (i+1/2,j+1/2)
!> points. Works for periodic, Neumann, and Dirichlet boundary
!> conditions.
!>
!> cof array:
!>
!>         cof(4)
!>           |
!>           |
!> cof(1)--cof(5)--cof(2)
!>           |
!>           |
!>         cof(3)
!>
!> Bernard Bunner 
subroutine mgdphpde(sxm,exm,sym,eym,nxm,nym,cof,sx,ex,sy,ey,nxf,nyf,bd,xl,yl)


sll_int32  :: sxm,exm,sym,eym,nxm,nym,sx,ex,sy,ey,nxf,nyf,bd(8)
sll_real64 :: cof(sxm-1:exm+1,sym-1:eym+1,6)
sll_real64 :: xl,yl
sll_real64 :: dlx,fodlxx,dly,fodlyy
sll_int32  :: i,j,im,jm,is,js,istep,jstep

! calculate off-diagonal terms

dlx=xl/float(nxm-1)
dly=yl/float(nym-1)

istep=(nxf-1)/(nxm-1)
jstep=(nyf-1)/(nym-1)

fodlxx=1.0d0/(dlx*dlx)
fodlyy=1.0d0/(dly*dly)

do j=sym,eym
  jm=2*jstep*j-3*(jstep-1)
  do i=sxm,exm
    im=2*istep*i-3*(istep-1)
    is=(im-istep)/2
    js=jm/2
    cof(i,j,1)=fodlxx
    is=(im+istep)/2
    cof(i,j,2)=fodlxx
    is=im/2
    js=(jm-jstep)/2
    cof(i,j,3)=fodlyy
    js=(jm+jstep)/2
    cof(i,j,4)=fodlyy
  end do
end do

! enforce wall BCs

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

! calculate diagonal term

do j=sym,eym
  do i=sxm,exm
    cof(i,j,5)=-(cof(i,j,1)+cof(i,j,2)+cof(i,j,3)+cof(i,j,4))
  end do
end do

end subroutine

!> For the old version of the multigrid code, determine coefficients 
!> for the pressure equation at all grid levels but the finest one.
!> The coefficients are determined from the values of the density
!> at integer nodes (i,j). Works only for periodic boundary conditions.
!>
!> cof array:
!>
!>         cof(4)
!>           |
!>           |
!> cof(1)--cof(5)--cof(2)
!>           |
!>           |
!>         cof(3)
!>
!> Bernard Bunner 
subroutine mgdppde(sxm,exm,sym,eym,nxm,nym,cof,     &
                   sxf,exf,syf,eyf,rf,xl,yl,bd)

sll_int32  :: sxm,exm,sym,eym,nxm,nym,sxf,exf,syf,eyf,bd(8)
sll_real64 :: cof(sxm-1:exm+1,sym-1:eym+1,6)
sll_real64 :: rf(sxf-1:exf+1,syf-1:eyf+1),xl,yl
sll_real64 :: dlx,odlxx,dly,odlyy
sll_int32  :: i,j,is,js
sll_real64 :: c1, c2, c3, c4

! calculate off-diagonal terms
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

end subroutine

!> Discretize the pde: set the coefficients of the cof matrix. Works
!> for periodic, Neumann, and Dirichlet boundary conditions.
!>
!> cof array:
!>
!>         cof(4)
!>           |
!>           |
!> cof(1)--cof(5)--cof(2)
!>           |
!>           |
!>         cof(3)
!>
subroutine mgdrpde(sxm,exm,sym,eym,nxm,nym,cof,xl,yl,bd)
sll_int32  :: sxm,exm,sym,eym,nxm,nym,bd(8)
sll_real64 :: cof(sxm-1:exm+1,sym-1:eym+1,6),xl,yl
sll_real64 :: dlx,odlxx,dly,odlyy
sll_int32  :: i,j

! calculate off-diagonal terms
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

!
! enforce Neumann BCs
!
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
!
! calculate diagonal term
!
do j=sym,eym
  do i=sxm,exm
    cof(i,j,5)=-(cof(i,j,1)+cof(i,j,2)+cof(i,j,3)+cof(i,j,4))
  end do
end do

end subroutine

!> For the old version of the multigrid code, set the fine grid values 
!> of the density in the work vector
!> Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
subroutine mgdrsetf(sxf,exf,syf,eyf,rf,r)

sll_int32  :: sxf,exf,syf,eyf
sll_real64 :: rf(sxf-1:exf+1,syf-1:eyf+1),r(sxf-1:exf+1,syf-1:eyf+1)
sll_int32  :: i,j

do j=syf-1,eyf+1
  do i=sxf-1,exf+1
    rf(i,j)=r(i,j)
  end do
end do

end subroutine

!------------------------------------------------------------------------
!> For the old version of the multigrid code, transfer values of the 
!> density from a finer to a coarser grid level. It is necessary to
!> exchange the boundary density data because the grid "shifts" to
!> the right as it becomes coarser. (In the new version of the
!> multigrid code, there is no such shift, hence no communication is 
!> needed).
!
!> Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
!------------------------------------------------------------------------
subroutine mgdrtrsf(sxc,exc,syc,eyc,nxc,nyc,rc,            &
                       sxf,exf,syf,eyf,nxf,nyf,rf,            &
                       comm2d,myid,neighbor,bd,itype,jtype)
sll_int32  :: sxc,exc,syc,eyc,nxc,nyc,sxf,exf,syf,eyf,nxf,nyf
sll_int32  :: comm2d,myid,neighbor(8),bd(8),itype,jtype

sll_real64 :: rc(sxc-1:exc+1,syc-1:eyc+1)
sll_real64 :: rf(sxf-1:exf+1,syf-1:eyf+1)

sll_int32  :: i,j,ic,jc,i1,i2,j1,j2

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

! exchange the boundary values (need only lines, not corner)
call gxch1lin(rc,comm2d,sxc,exc,syc,eyc,neighbor,bd,itype,jtype)

end subroutine

subroutine mgdsetf(sxf,exf,syf,eyf,phi,rhs,phif,rhsf)
sll_int32  :: sxf,exf,syf,eyf
sll_real64 :: phi(sxf-1:exf+1,syf-1:eyf+1),rhs(sxf-1:exf+1,syf-1:eyf+1)
sll_real64 :: phif(sxf-1:exf+1,syf-1:eyf+1),rhsf(sxf-1:exf+1,syf-1:eyf+1)
!------------------------------------------------------------------------
! Set the fine grid values in the work vector
!
! Code      : mgd2, 2-D parallel multigrid solver
! Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
! Called in : mgdsolver
! Calls     : --
!------------------------------------------------------------------------
sll_int32  :: i,j

do j=syf-1,eyf+1
  do i=sxf-1,exf+1
    phi(i,j)=phif(i,j)
    rhs(i,j)=rhsf(i,j)
  end do
end do

end subroutine

!> Free the MPI datatypes associated witht the multigrid code
subroutine mgdend(ngrid)

sll_int32 :: k,ierr
sll_int32 :: ngrid

do k=1,ngrid-1
   call MPI_TYPE_FREE(ikdatatype(k),ierr)
   call MPI_TYPE_FREE(jkdatatype(k),ierr)
   call MPI_TYPE_FREE(ijkdatatype(k),ierr)
 if (.not. WMGD) then
   call MPI_TYPE_FREE(irdatatype(k),ierr)
   call MPI_TYPE_FREE(jrdatatype(k),ierr)
   call MPI_TYPE_FREE(ijrdatatype(k),ierr)
 endif
end do

!call MPI_FINALIZE(ierr)


end subroutine

!> Add correction from coarse grid level to fine grid level. Uses
!> bilinear interpolation for the old version of the multigrid code,
!> and area weighting for its new version.
!>
!> Tested for the case where coarsifying takes place in all directions
subroutine mgdcor(sxf,exf,syf,eyf,nxf,nyf,phif,  &
                  sxc,exc,syc,eyc,nxc,nyc,phic,  &
                  sx1,ex1,sy1,ey1,bd,phibc)

sll_int32  :: sxf,exf,syf,eyf,nxf,nyf
sll_int32  :: sxc,exc,syc,eyc,nxc,nyc,sx1,ex1,sy1,ey1,bd(8)
sll_real64 :: phif(sxf-1:exf+1,syf-1:eyf+1)
sll_real64 :: phic(sxc-1:exc+1,syc-1:eyc+1),phibc(4)

sll_int32 :: i,j,ic,jc,i1,i2,j1,j2

if (WMGD) then
!------------------------------------------------------------------------
! new version: the correction is the weighted average of either two 
! or four points at the coarser grid level depending on whether 
! coarsifying takes place in all directions or not
!
if (nxf.eq.nxc) then
  do jc=syc-1,eyc
    j=2*jc-1
    do i=sxf-1,exf+1
      ic=i
      phif(i,j)=phif(i,j)+(3.0d0*phic(ic,jc)+phic(ic,jc+1))/4.0d0
      phif(i,j+1)=phif(i,j+1)+(phic(ic,jc)+3.0d0*phic(ic,jc+1))/4.0d0
    end do
  end do
else if (nyf.eq.nyc) then
  do j=syf-1,eyf+1
    jc=j
    do ic=sxc-1,exc
      i=2*ic-1
      phif(i,j)=phif(i,j)+(3.0d0*phic(ic,jc)+phic(ic+1,jc))/4.0d0
      phif(i+1,j)=phif(i+1,j)+(phic(ic,jc)+3.0d0*phic(ic+1,jc))/4.0d0
    end do
  end do
else
  do jc=syc-1,eyc
    j=2*jc-1
    do ic=sxc-1,exc
      i=2*ic-1
      phif(i,j)=phif(i,j)+  &
        (9.0d0*phic(ic,jc)+3.0d0*phic(ic+1,jc)  &
        +3.0d0*phic(ic,jc+1)+phic(ic+1,jc+1))/16.0d0
      phif(i+1,j)=phif(i+1,j)+  &
        (3.0d0*phic(ic,jc)+9.0d0*phic(ic+1,jc)  &
        +phic(ic,jc+1)+3.0d0*phic(ic+1,jc+1))/16.0d0
      phif(i,j+1)=phif(i,j+1)+  &
        (3.0d0*phic(ic,jc)+phic(ic+1,jc)  &
        +9.0d0*phic(ic,jc+1)+3.0d0*phic(ic+1,jc+1))/16.0d0
      phif(i+1,j+1)=phif(i+1,j+1)+  &
        (phic(ic,jc)+3.0d0*phic(ic+1,jc)  &
        +3.0d0*phic(ic,jc+1)+9.0d0*phic(ic+1,jc+1))/16.0d0
    end do
  end do
end if
!
! impose Neumann and Dirichlet boundary conditions
! TEMP: periodicity is not enforced to save one call to gxch1lin;
! check whether it has an impact or not...
!
      call mgdbdry(sxf,exf,syf,eyf,phif,bd,phibc)
else  
!------------------------------------------------------------------------
! old version
!
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
!
! identity at the points of the fine grid which have odd indices 
! in i and j
!
do jc=sy1,ey1
  j=j1*(2*jc-1)+j2*jc
  do ic=sx1,ex1
    i=i1*(2*ic-1)+i2*ic
    phif(i,j)=phif(i,j)+phic(ic,jc)
  end do
end do
!
! interpolation of the two neighboring values for the points of the 
! fine grid with even index for i and odd index for j
!
if (nxc.lt.nxf) then
  do jc=sy1,ey1
    j=j1*(2*jc-1)+j2*jc
    do ic=sxc-1,exc
      i=2*ic
      phif(i,j)=phif(i,j)+0.5d0*(phic(ic,jc)+phic(ic+1,jc))
    end do
  end do
end if
!
! interpolation of the two neighboring values for the points of the 
! fine grid with odd index for i and even index for j
!
if (nyc.lt.nyf) then
  do jc=syc-1,eyc
    j=2*jc
    do ic=sx1,ex1
      i=i1*(2*ic-1)+i2*ic
      phif(i,j)=phif(i,j)+0.5d0*(phic(ic,jc)+phic(ic,jc+1))
    end do
  end do
end if
!
! interpolation of the four neighboring values for the points of the
! fine grid with even indices for i and j
!
if (nxc.lt.nxf.and.nyc.lt.nyf) then
  do jc=syc-1,eyc
    j=2*jc
    do ic=sxc-1,exc
      i=2*ic
      phif(i,j)=phif(i,j)+0.25d0*(phic(ic,jc)+phic(ic+1,jc) &
                                 +phic(ic,jc+1)+phic(ic+1,jc+1))
    end do
  end do
end if
endif

end subroutine


subroutine mgdrelax(sxm,exm,sym,eym,phi,cof,iters,comm2d,myid, &
                    neighbor,bd,phibc,itype,jtype)


sll_int32 :: sxm,exm,sym,eym,iters
sll_int32 :: comm2d,myid,neighbor(8),bd(8),itype,jtype
sll_real64 :: phi(sxm-1:exm+1,sym-1:eym+1)
sll_real64 :: cof(sxm-1:exm+1,sym-1:eym+1,6),phibc(4)
!------------------------------------------------------------------------
! Gauss-Seidel point relaxation with Red & Black ordering. Works for
! periodic, Neumann, and Dirichlet boundary conditions.
!
! Code      : mgd2, 2-D parallel multigrid solver
! Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
! Called in : mgdkcyc
! Calls     : mgdbdry, gxch1lin
!------------------------------------------------------------------------
sll_int32 :: rb,it,ipass,i,j
!
! do iters sweeps in the subdomain
!
do it=1,iters
  rb=mod(sxm,2)
  do ipass=1,2
    do j=sym,eym
      do i=sxm+rb,exm,2
        phi(i,j)=(cof(i,j,6)-(cof(i,j,1)*phi(i-1,j)               &
                             +cof(i,j,2)*phi(i+1,j)               &
                             +cof(i,j,3)*phi(i,j-1)               &
                             +cof(i,j,4)*phi(i,j+1)))/cof(i,j,5)
      end do
      rb=1-rb
    end do
    rb=1-mod(sxm,2)
if (WMGD) then
!
! new version: impose Neumann and Dirichlet boundary conditions
!
    call mgdbdry(sxm,exm,sym,eym,phi,bd,phibc)
endif
  end do
end do
!
! Exchange boundary data only once at the end. Since the number
! of relaxation sweeps at each level is characteristically small
! (1 or 2 are common values), this does not damage the convergence
! rate too badly. Overall, I have found a significant reduction
! in execution time. This also imposes the periodic BCs.
!
call gxch1lin(phi,comm2d,sxm,exm,sym,eym,neighbor,bd,itype,jtype)
if (WMGD) then
!
! new version: impose Neumann and Dirichlet boundary conditions
!
call mgdbdry(sxm,exm,sym,eym,phi,bd,phibc)
endif 

end subroutine

!> Calculate the residual and restrict it to the coarser level. In
!> the new version, the restriction involves 4 points. In the old
!> version, it involves 9 points (5 if half-weighting).
!>
!> Bernard Bunner
subroutine mgdrestr(sxc,exc,syc,eyc,nxc,nyc,phic,rhsc,   &
                    sxf,exf,syf,eyf,nxf,nyf,phif,cof,    &
                    resf,iresw,comm2d,myid,neighbor,bd,  &
                    itype,jtype,ijtype)

sll_int32  :: sxc,exc,syc,eyc,nxc,nyc,iresw
sll_int32  :: sxf,exf,syf,eyf,nxf,nyf
sll_int32  :: comm2d,myid,neighbor(8),bd(8),itype,jtype,ijtype
sll_real64 :: phic(sxc-1:exc+1,syc-1:eyc+1)
sll_real64 :: rhsc(sxc-1:exc+1,syc-1:eyc+1)
sll_real64 :: phif(sxf-1:exf+1,syf-1:eyf+1)
sll_real64 :: resf(sxf-1:exf+1,syf-1:eyf+1)
sll_real64 :: cof(sxf-1:exf+1,syf-1:eyf+1,6)
sll_int32  :: i,j,isrt,jsrt,iinc,jinc,ic,jc
!------------------------------------------------------------------------
do jc=syc-1,eyc+1
  do ic=sxc-1,exc+1
    phic(ic,jc)=0.0d0
    rhsc(ic,jc)=0.0d0
  end do
end do
!
! calculate residual
!
do j=syf,eyf
  do i=sxf,exf
    resf(i,j)=cof(i,j,6)-(cof(i,j,1)*phif(i-1,j) &
                         +cof(i,j,2)*phif(i+1,j) &
                         +cof(i,j,3)*phif(i,j-1) &
                         +cof(i,j,4)*phif(i,j+1) &
                         +cof(i,j,5)*phif(i,j))
  end do
end do
if (WMGD) then
!------------------------------------------------------------------------
! new version: calculate the right-hand side at the coarser grid
! level from the averages of the values at the 4 surrounding points;
! if there is no coarsifying in one direction, only 2 points are
! used; no exchange of boundary data is necessary
!
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
      rhsc(ic,jc)=0.25d0*(resf(i,j)+resf(i+1,j) &
                   +resf(i,j+1)+resf(i+1,j+1))
    end do
  end do
end if
else  
!------------------------------------------------------------------------
! old version: have to exchange boundary data; if full-weighting, 
! need to exchange also corner points
!
call gxch1lin(resf,comm2d,sxf,exf,syf,eyf,neighbor,bd,itype,jtype)
if (iresw.eq.1) then
  call gxch1cor(resf,comm2d,sxf,exf,syf,eyf,neighbor,bd,ijtype)
end if
!
! restrict it to coarser level
!
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

if (iresw.eq.1) then
!
! use full weighting
!
  j=jsrt
  do jc=syc,eyc
    i=isrt
    do ic=sxc,exc
      rhsc(ic,jc)=0.25d0*resf(i,j)                        &
                 +0.125d0*(resf(i+1,j)+resf(i-1,j)        &
                          +resf(i,j+1)+resf(i,j-1))       &
                 +0.0625d0*(resf(i+1,j+1)+resf(i+1,j-1)   &
                           +resf(i-1,j-1)+resf(i-1,j+1))
      i=i+iinc
    end do
    j=j+jinc
  end do
else if (iresw.eq.2) then
!
! use half-weighting
!
  j=jsrt
  do jc=syc,eyc
    i=isrt
    do ic=sxc,exc
      rhsc(ic,jc)=0.5d0*resf(i,j)                     &
                 +0.125d0*(resf(i+1,j)+resf(i-1,j)    &
                          +resf(i,j+1)+resf(i,j-1))
      i=i+iinc
    end do
    j=j+jinc
  end do
end if
endif

end subroutine


!> Calculate the error between the new and old iterates of phi and 
!> save the new iterate into the phio array.
subroutine mgderr(relmax,sxm,exm,sym,eym,phio,phin,comm2d)

sll_int32  :: sxm,exm,sym,eym,comm2d
sll_real64 :: relmax
sll_real64 :: phio(sxm-1:exm+1,sym-1:eym+1),phin(sxm-1:exm+1,sym-1:eym+1)
sll_real64 :: phloc,reloc
sll_int32  :: i,j,ierr
!
! calculate local error
!
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

! global reduce across all processes
call MPI_ALLREDUCE(reloc,relmax,1,MPI_REAL8,MPI_MAX,comm2d,ierr)

! save new values into ouput array
do j=sym-1,eym+1
  do i=sxm-1,exm+1
    phio(i,j)=phin(i,j)
  end do
end do

end subroutine

!> Enforce the Neumann and Dirichlet boundary conditions
subroutine mgdbdry(sxm,exm,sym,eym,phi,bd,phibc)


sll_int32  :: sxm,exm,sym,eym,bd(8)
sll_real64 :: phi(sxm-1:exm+1,sym-1:eym+1),phibc(4)
sll_int32 :: i,j

if (bd(1).eq.1) then
   do j=sym-1,eym+1
      phi(exm+1,j)=phi(exm,j)
   end do
else if (bd(1).eq.2) then
   do j=sym-1,eym+1
      phi(exm+1,j)=2.0d0*phibc(1)-phi(exm,j)
   end do
end if
if (bd(5).eq.1) then
   do j=sym-1,eym+1
      phi(sxm-1,j)=phi(sxm,j)
   end do
else if (bd(5).eq.2) then
   do j=sym-1,eym+1
      phi(sxm-1,j)=2.0d0*phibc(3)-phi(sxm,j)
   end do
end if
if (bd(3).eq.1) then
   do i=sxm-1,exm+1
      phi(i,sym-1)=phi(i,sym)
   end do
else if (bd(3).eq.2) then
   do i=sxm-1,exm+1
      phi(i,sym-1)=2.0d0*phibc(2)-phi(i,sym)
   end do
end if
if (bd(7).eq.1) then
   do i=sxm-1,exm+1
      phi(i,eym+1)=phi(i,eym)
   end do
else if (bd(7).eq.2) then
   do i=sxm-1,exm+1
      phi(i,eym+1)=2.0d0*phibc(4)-phi(i,eym)
   end do
end if

end subroutine

!> Do one multigrid K-cycle
!> K=1 -> V-cycle
!> K=2 -> W-cycle
subroutine mgdkcyc(work,res,kcur,kcycle,iprer,ipost,iresw, &
                   comm2d,myid,neighbor,bd,phibc)

sll_int32  :: kcur,kcycle,iprer,ipost,iresw
sll_int32  :: comm2d,myid,neighbor(8),bd(8)
sll_real64 :: work(*),res(*),phibc(4,*)

sll_int32 :: sxf,exf,syf,eyf,nxf,nyf,ipf,icf
sll_int32 :: sxc,exc,syc,eyc,nxc,nyc,ipc,irc
sll_int32 :: klevel,itype,jtype,ijtype,kount(20),l,nrel
sll_int32 :: sx1,ex1,sy1,ey1

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

end subroutine
end module mgd2
