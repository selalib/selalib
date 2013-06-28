!> Multigrid solver in 2 dimensions on a cartesian grid
!> \authors Pierre Navaro
!> \remark this a refactored version of Bernard Bunner code 
!> (http://cs-www.cs.yale.edu/homes/douglas-craig/mgnet-codes-bunner.html)
module mgd2
#include "sll_working_precision.h"
use mpi
use gxch1_2d

implicit none

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
   sll_int32  :: maxcy, kcycle, iprer, ipost, iresw
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
!> Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
subroutine initialize_mgd2(my_block, my_mg, nerror)
#include "mgd2.h"

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
100 format(/,'ERROR in mgdinit: nx=',i3,' is not a multiple of ', &
           'nxprocs=',i3,/,'cannot use the new version of the ',  &
           'multigrid code',/)
110 format(/,'ERROR in mgdinit: ny=',i3,' is not a multiple of ', &
           'nyprocs=',i3,/,'cannot use the new version of the ',  &
           'multigrid code',/)
# else
!
! check that the old version is not used with non-periodic BCs
!
if (my_mg%ibdry.ne.0.or.my_mg%jbdry.ne.0) then
  write(6,120) my_mg%ibdry,my_mg%jbdry
  nerror=1
  return
end if
120  format(/,'ERROR in mgdinit: ibdry=',i2,' jbdry=',i2, &
            /,'cannot use the old version of the multigrid code', &
            /,'boundary conditions that are not periodic', &
            /,'-> change compiler directive to 1 in compdir.inc', &
            /,'   and recompile the multigrid code',/)
# endif
!
! check that the dimensions are correct
!
i=my_block%ixp*2**(my_block%iex-1)+1
if ((my_mg%nx+1).ne.i) then
  write(6,130) my_mg%nx+1,i
  nerror=1
  return
end if
j=my_block%jyq*2**(my_block%jey-1)+1
if ((my_mg%ny+1).ne.j) then
  write(6,140) my_mg%ny+1,j
  nerror=1
  return
end if
130  format(/,'ERROR in mgdinit: nxp1=',i3,' <> ixp*2**(iex-1)+1=', &
            i3,/,'-> adjust the multigrid parameters ixp and iex', &
            ' in main',/)
140  format(/,'ERROR in mgdinit: nyp1=',i3,' <> jyq*2**(jey-1)+1=', &
            i3,/,'-> adjust the multigrid parameters jyq and jey', &
            ' in main',/)
!
! check that the number of points at the coarser level is not smaller
! than the number of processes in either direction
!
if (my_block%ixp.lt.my_mg%nxprocs) then
  write(6,150) my_block%ixp,my_mg%nxprocs
  nerror=1
  return
end if
if (my_block%jyq.lt.my_mg%nyprocs) then
  write(6,160) my_block%jyq,my_mg%nyprocs
  nerror=1
  return
end if
150  format(/,'ERROR in mgdinit: ixp=',i3,' < nxprocs=',i3,/, &
            ' there must be at least one grid point at the ', &
            'coarsest grid level',/, &
            '-> increase ixp and decrease iex correspondingly', &
            ' in main',/)
160  format(/,'ERROR in mgdinit: jyq=',i3,' < nyprocs=',i3,/, &
            ' there must be at least one grid point at the ', &
            'coarsest grid level',/, &
            '-> increase jyq and decrease jey correspondingly', &
            ' in main',/)
!
! check that coarsifying takes place in all directions at the finest
! grid level
!
if (my_block%ngrid.gt.1) then
  if (my_block%iex.eq.1) then
    write(6,170) my_block%ngrid,my_block%iex
    nerror=1
    return
  end if
  if (my_block%jey.eq.1) then
    write(6,180) my_block%ngrid,my_block%jey
    nerror=1
    return
  end if
end if
170  format(/,'ERROR in mgdinit: ngrid=',i3,' iex=',i3, &
            /,'no coarsifying at the finest grid level in x-direction', &
            /,'this is not allowed by the mutligrid code',/)
180  format(/,'ERROR in mgdinit: ngrid=',i3,' jey=',i3, &
            /,'no coarsifying at the finest grid level in y-direction', &
            /,'this is not allowed by the mutligrid code',/)
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
# if xdebug1
!
! print out the indices and determine the size of the MPI messages
! as a rough check
! 
write(6,*) 'size of the multigrid phi-messages'
do k=my_block%ngrid,1,-1
  call MPI_TYPE_SIZE(ikdatatype(k),nsiz1,ierr)
  call MPI_TYPE_SIZE(jkdatatype(k),nsiz2,ierr)
  call MPI_TYPE_SIZE(ijkdatatype(k),nsiz3,ierr)
  write(6,*) 'myid: ',my_block%id,' k=',k,' sxk=',sxk(k),' exk=', &
                exk(k),' syk=',syk(k),' eyk=',eyk(k), &
                ' size of datatypes: ',nsiz1,nsiz2,nsiz3
end do
# endif
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
  my_mg%phibc(j,my_block%ngrid)=my_mg%vbc(j)
end do
do k=my_block%ngrid-1,1,-1
  do j=1,4
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
# if xdebug1
  write(6,*) 'myid: ',my_block%id,' k=',k,' sxi=',sxi(k),' exi=', &
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
300   format('ERROR in mgdinit: no coarsifying between level ', &
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
# if xdebug1
!
! print out the indices and determine the size of the MPI messages
! as a rough check
! 
write(6,*) 'size of the r-messages'
do k=my_block%ngrid-1,1,-1
  call MPI_TYPE_SIZE(irdatatype(k),nsiz1,ierr)
  call MPI_TYPE_SIZE(jrdatatype(k),nsiz2,ierr)
  call MPI_TYPE_SIZE(ijrdatatype(k),nsiz3,ierr)
  write(6,*) 'myid: ',my_block%id,' k=',k,' sxr=',sxr(k),' exr=', &
                exr(k),' syr=',syr(k),' eyr=',eyr(k), &
                ' size of datatypes: ',nsiz1,nsiz2,nsiz3
end do
# endif
# endif

return
end subroutine

!>------------------------------------------------------------------------
!> Parallel multigrid solver in 2-D cartesian coordinates for the
!> elliptic equation:      div(cof(grad(phif)))=rhsf
!>
!> isol=1 -> density
!> isol=2 -> pressure
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
!> Code      : mgd2, 2-D parallel multigrid solver
!> Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
!> Called in : main
!> Calls     : mgdrpde, mgdpfpde, mgdphpde, mgdrsetf, mgdppde, mgdrtrsf,
!>              -> discretize the pde
!>             mgdsetf
!>               -> set the initial guesses and the right-hand side
!>             mgdkcyc, mgderr,
!>               -> do the actual cycling
!>             gscale, gxch1lin, gxch1cor,
!>               -> rescale pressure and density around average values
!>------------------------------------------------------------------------
subroutine mgd2_solver(my_block,my_mg,phif,rhsf,r,work, &
                     rro,iter,nprscr,nerror)
#include "mgd2.h"
sll_real64      :: phif(:,:),rhsf(:,:)
sll_real64      :: r(:,:)
sll_real64      :: work(*),rro
sll_int32       :: iter,nerror
type(block)     :: my_block
type(mg_solver) :: my_mg

logical    :: nprscr

sll_int32  :: k, icf
sll_real64 :: avo,acorr, relmax
sll_int32  :: sxf,exf,syf,eyf,nxf,nyf,sxc,exc,syc,eyc,nxc,nyc
sll_int32  :: ipf,irf,ipc,irc,sxm,exm,sym,eym,nxm,nym,ip,ic,kcur
sll_int32  :: itype,jtype,ijtype,lev,ir1,ir2
# if DEBUG
sll_real64 :: tinitial
tinitial=MPI_WTIME()
# endif
!------------------------------------------------------------------------
! discretize pde at all levels
!
if (my_mg%isol.eq.1) then
!
! density: only have to set geometric factors
!
  do k=1,my_block%ngrid
    sxm=sxk(k)
    exm=exk(k)
    sym=syk(k)
    eym=eyk(k)
    nxm=nxk(k)
    nym=nyk(k)
    ic=kcbgn(k)
    call mgdrpde(sxm,exm,sym,eym,nxm,nym,work(ic), &
                 my_mg%xl,my_mg%yl,my_block%bd)
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
  sxf=sxk(my_block%ngrid)
  exf=exk(my_block%ngrid)
  syf=syk(my_block%ngrid)
  eyf=eyk(my_block%ngrid)
  nxf=nxk(my_block%ngrid)
  nyf=nyk(my_block%ngrid)
  icf=kcbgn(my_block%ngrid)
  call mgdpfpde(sxf,exf,syf,eyf,nxf,nyf,work(icf),r,my_mg%xl,my_mg%yl,my_block%bd)
# if WMGD
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
                  sx,ex,sy,ey,nxf,nyf,r,my_block%bd,my_mg%xl,my_mg%yl)
  end do
# else
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
# endif
      end if
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
if (my_mg%isol.eq.1) then
  avo=rro
else
  avo=0.0d0
end if
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
# if WMGD
!
! impose wall and Dirichlet BCs
!
call mgdbdry(sx,ex,sy,ey,phif,my_block%bd,my_mg%vbc)
# endif

if (my_mg%isol.eq.1) then
  if (nprscr.and.my_block%id.eq.0) write(6,110) relmax,iter,acorr
110     format('  R MGD     err=',e8.3,' iters=',i5,' rcorr=',e9.3)
else
  if (nprscr.and.my_block%id.eq.0) write(6,120) relmax,iter,acorr
120     format('  P MGD     err=',e8.3,' iters=',i5,' pcorr=',e9.3)
end if

return
end subroutine

!> Define the 3 derived datatypes needed to communicate the boundary
!> data of (sx-1:ex+1,sy-1:ey+1) arrays between 'myid' and its 8
!> neighbors
subroutine grid1_type(my_block,itype,jtype,ijtype,sx,ex,sy,ey)
#include "mgd2.h"
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

end module mgd2
