module sll_multigrid_2d
#include "sll_working_precision.h"
#include "sll_constants.h"
#include "sll_memory.h"

use sll_boundary_condition_descriptors
use sll_multigrid_base
use sll_collective, only: sll_boot_collective,      &
                          sll_world_collective,     &
                          sll_halt_collective,      &
                          sll_get_collective_rank,  &
                          sll_get_collective_size

use sll_gnuplot_parallel
use sll_remapper
use mgd2

implicit none

type, public, extends(sll_multigrid_solver) :: sll_multigrid_solver_2d

   sll_real64 :: x_min
   sll_real64 :: x_max
   sll_real64 :: y_min
   sll_real64 :: y_max

   sll_real64 :: xl
   sll_real64 :: yl
   
   sll_int32  :: nc_x
   sll_int32  :: nc_y
   sll_int32  :: nprocs_x
   sll_int32  :: nprocs_y

   type(block)     :: block

end type sll_multigrid_solver_2d

interface initialize
   module procedure initialize_multigrid_2d
end interface initialize
interface solve
   module procedure solve_multigrid_2d
end interface solve
interface delete
   module procedure delete_multigrid_2d
end interface delete

sll_int32, private :: i, j, k
sll_int32, private :: error
sll_real64, dimension(:), allocatable :: work

contains

subroutine delete_multigrid_2d(this)
type(sll_multigrid_solver_2d) :: this

SLL_DEALLOCATE_ARRAY(work, error)

end subroutine delete_multigrid_2d

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
subroutine initialize_multigrid_2d(this, layout,             &
                                   x_min, x_max, nc_x, bc_x, &
                                   y_min, y_max, nc_y, bc_y, &
                                   vbc_x, vbc_y)

type(sll_multigrid_solver_2d) :: this
type(layout_2d), pointer      :: layout
sll_real64, intent(in)        :: x_min
sll_real64, intent(in)        :: x_max
sll_real64, intent(in)        :: y_min
sll_real64, intent(in)        :: y_max
sll_int32 , intent(in)        :: nc_x
sll_int32 , intent(in)        :: nc_y
sll_int32 , intent(in)        :: bc_x
sll_int32 , intent(in)        :: bc_y
sll_int32 , optional          :: vbc_x(2)
sll_int32 , optional          :: vbc_y(2)

sll_int32  :: nxp2,nyp2,nx,ny,nxprocs,nyprocs,ixp,jyq,iex,jey
sll_int32  :: nwork
sll_int32  :: ngrid,nxdim,nydim

logical   :: periods(2)
sll_int32 :: numprocs,myid,sx,ex,sy,ey,neighbor(8),bd(8)
sll_int32 :: realtype
sll_int32 :: ierr,m0,m1,dims(2)
sll_int32 :: nbrright,nbrbottom,nbrleft,nbrtop

integer :: statut(MPI_STATUS_SIZE)

sll_int32  :: nxf,nyf,nxm,nym,sxm,exm,sym,eym,nxc,nyc, kps

this%layout => layout
numprocs    = sll_get_collective_size(sll_world_collective)
myid        = sll_get_collective_rank(sll_world_collective)
this%comm2d = sll_world_collective%comm
realtype    = MPI_REAL8

nxprocs       = int(sqrt(real(numprocs,f64)))
nyprocs       = nxprocs
this%nprocs_x = nxprocs
this%nprocs_y = nyprocs
this%x_min    = x_min
this%x_max    = x_max
this%y_min    = y_min
this%y_max    = y_max
this%nc_x     = nc_x
this%nc_y     = nc_y

if (nxprocs*nyprocs /= numprocs) goto 1000

nx = nc_x
ny = nc_y

nxp2 = nc_x+2
nyp2 = nc_y+2

sx = get_layout_2d_i_min(this%layout,myid)+1
ex = get_layout_2d_i_max(this%layout,myid)+1
sy = get_layout_2d_j_min(this%layout,myid)+1
ey = get_layout_2d_j_max(this%layout,myid)+1

ixp = 4    ! the resolution of the global domain,
jyq = 4    ! such that :
iex = ceiling(log((nx-1.)/ixp)/log(2.))+1
jey = ceiling(log((ny-1.)/jyq)/log(2.))+1

this%maxcy  = 300     ! maximum number of multigrid cycle
this%tolmax = 1.0d-05 ! desired accuracy
this%kcycle = 1       ! 1 : V-cycles 2 : W-cycles
this%iprer  = 2       ! relaxation steps before restricting the residual
this%ipost  = 1       ! relaxation steps after correction
this%iresw  = 1       ! 1 : fully weighted residuals 2 : half-weighted residuals

nwork  = (4*nx*ny*8)/(3*nxprocs*nyprocs)+(64*(nx+ny))/3+(32*4)/3

ngrid  = max(iex,jey)
nxdim  = ceiling(float(nxp2-2)/float(nxprocs))+2
nydim  = ceiling(float(nyp2-2)/float(nyprocs))+2

if ((ex-sx+3) > nxdim) then
  write(6,110) myid,nxdim,ex-sx+3
  goto 1000
end if
if ((ey-sy+3).gt.nydim) then
  write(6,120) myid,nydim,ey-sy+3
  goto 1000
end if

if (.not. allocated(work)) then
   SLL_CLEAR_ALLOCATE(work(1:nwork), error)
end if

! open file for output of messages and check that the number of 
! processes is correct

m1 = mod(myid,10)+1
m0 = mod(myid/10,10)+1
if (numprocs.ne.(nxprocs*nyprocs)) then
  write(6,100)
  goto 1000
end if

! create Cartesian topology with periodic boundary conditions

dims(1)    = nxprocs
dims(2)    = nyprocs
periods(1) = .true. !non periodic in r
periods(2) = .true.  !periodic in theta
call MPI_CART_CREATE(MPI_COMM_WORLD,2,dims,periods,.true.,this%comm2d,ierr)
call MPI_COMM_RANK(this%comm2d,myid,ierr)

if(bc_x == SLL_PERIODIC) then
   this%ibdry  = 0    ! Periodic bc in direction x
else if(bc_x == SLL_DIRICHLET) then
   this%ibdry  = 1    ! Periodic bc in direction x
   if (present(vbc_x)) then
      this%vbc(1:2) = vbc_x ! Value for x_min, x_max
   else
      this%vbc(1:2) = 0.0_f64
   end if
else
   print*, " multigrid solver wrong boundary conditions in x"
   goto 1000
end if

if(bc_y == SLL_PERIODIC) then
   this%jbdry  = 0    ! Periodic bc in direction y
else if(bc_y == SLL_DIRICHLET) then
   this%jbdry  = 1    ! Periodic bc in direction y
   if (present(vbc_y)) then
      this%vbc(3:4) = vbc_y ! Value for y_min, y_max
   else
      this%vbc(3:4) = 0.0_f64
   end if
else
   print*, " multigrid solver wrong boundary conditions in y"
   goto 1000
end if


bd(:)  = 0

! find the neighbors conventions
!
!     6 |     7     | 8
!       |           |
!     -----------------
!       |           |
!       |           |
!     5 |   myid    | 1 
!       |           |
!       |           |
!     -----------------
!       |           |
!     4 |     3     | 2

call MPI_CART_SHIFT(this%comm2d,0,1,nbrleft,nbrright,ierr)
call MPI_CART_SHIFT(this%comm2d,1,1,nbrbottom,nbrtop,ierr)

neighbor(1)=nbrright
neighbor(3)=nbrbottom
neighbor(5)=nbrleft
neighbor(7)=nbrtop

call MPI_BARRIER(this%comm2d,ierr)
call MPI_SENDRECV(neighbor(3),1,MPI_INTEGER,nbrright,0, &
                  neighbor(4),1,MPI_INTEGER,nbrleft,0, &
                  this%comm2d,statut,ierr)
call MPI_SENDRECV(neighbor(7),1,MPI_INTEGER,nbrright,1, &
                  neighbor(6),1,MPI_INTEGER,nbrleft,1, &
                  this%comm2d,statut,ierr)
call MPI_SENDRECV(neighbor(3),1,MPI_INTEGER,nbrleft,0, &
                  neighbor(2),1,MPI_INTEGER,nbrright,0, &
                  this%comm2d,statut,ierr)
call MPI_SENDRECV(neighbor(7),1,MPI_INTEGER,nbrleft,1, &
                  neighbor(8),1,MPI_INTEGER,nbrright,1, &
                  this%comm2d,statut,ierr)
!do i=1,8
!  write(6,*) 'neighbor: ',neighbor(i),' bd: ',bd(i)
!end do

this%block%id       = myid
this%block%sx       = sx
this%block%ex       = ex
this%block%sy       = sy
this%block%ey       = ey
this%block%ngrid    = ngrid
this%block%neighbor = neighbor
this%block%bd       = bd

this%block%ixp      = ixp
this%block%jyq      = jyq
this%block%iex      = iex
this%block%jey      = jey

!
! check that the dimensions are correct
!
i=ixp*2**(iex-1)+1
if ((this%nc_x+1).ne.i) then
  write(6,130) this%nc_x+1,i
  goto 1000
end if
j=jyq*2**(jey-1)+1
if ((this%nc_y+1).ne.j) then
  write(6,140) this%nc_y+1,j
  goto 1000
end if
!
! check that the number of points at the coarser level is not smaller
! than the number of processes in either direction
!
if (ixp.lt.this%nprocs_x) then
  write(6,150) ixp,this%nprocs_x
  goto 1000
end if
if (jyq.lt.this%nprocs_y) then
  write(6,160) jyq,this%nprocs_y
  goto 1000
end if
! check that coarsifying takes place in all directions at the finest
! grid level
!
if (ngrid.gt.1) then
  if (iex.eq.1) then
    write(6,170) ngrid,iex
    goto 1000
  end if
  if (jey.eq.1) then
    write(6,180) ngrid,jey
    goto 1000
  end if
end if

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

if (this%ibdry.ne.0.or.this%jbdry.ne.0) then

  write(6,220) this%ibdry,this%jbdry
  WMGD = .true.
  ! check that the number of processes in each
  ! direction divides the number of points in that direction 
  if (mod(this%nc_x,this%nprocs_x).ne.0) then
    write(6,230) this%nc_x,this%nprocs_x
    goto 1000
  end if
  if (mod(this%nc_y,this%nprocs_y).ne.0) then
    write(6,210) this%nc_y,this%nprocs_y
    goto 1000  
  end if

else

  WMGD = .false.

endif

!------------------------------------------------------------------------
! define all grid levels
! I have adopted the same notations as in Mudpack as far as possible.
! When a confusion was possible, I added a suffix 'm' to the name
! of the variables. For example, nxm is nx+1 for the multigrid
! code whereas nx means nxp2-2 in the rest of the code.
!
do k=1,this%block%ngrid
  nxk(k)=this%block%ixp*2**(max(k+this%block%iex-this%block%ngrid,1)-1)+1
  nyk(k)=this%block%jyq*2**(max(k+this%block%jey-this%block%ngrid,1)-1)+1
end do

! for all grid levels, set the indices of the subdomain the process
! 'myid' work on, as well as the datatypes needed for the exchange
! of boundary data

nxf=nxk(this%block%ngrid)
nyf=nyk(this%block%ngrid)
sxk(this%block%ngrid)=this%block%sx
exk(this%block%ngrid)=this%block%ex
syk(this%block%ngrid)=this%block%sy
eyk(this%block%ngrid)=this%block%ey
call grid1_type(ikdatatype(this%block%ngrid), &
                jkdatatype(this%block%ngrid), &
                ijkdatatype(this%block%ngrid),&
                this%block%sx, &
                this%block%ex, &
                this%block%sy, &
                this%block%ey)
do k=this%block%ngrid-1,1,-1
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
                  sxk(k),exk(k),syk(k),eyk(k))
end do
!
! set work space indices for phi, cof at each grid level
! check that there is sufficient work space
!
kps=1
do k=this%block%ngrid,1,-1
  sxm=sxk(k)
  exm=exk(k)
  sym=syk(k)
  eym=eyk(k)
  kpbgn(k)=kps
  kcbgn(k)=kpbgn(k)+(exm-sxm+3)*(eym-sym+3)
  kps=kcbgn(k)+6*(exm-sxm+3)*(eym-sym+3)
end do

!if (kps.gt.nwork) then
!  write(6,200) kps,nwork,this%block%id
!  error=1
!  retur
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
   this%phibc(:,:)=0.0d0
   this%phibc(:,this%block%ngrid)=this%vbc(:)
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
      
sxi(this%block%ngrid)=sxk(this%block%ngrid)-1
exi(this%block%ngrid)=exk(this%block%ngrid)+1
syi(this%block%ngrid)=syk(this%block%ngrid)-1
eyi(this%block%ngrid)=eyk(this%block%ngrid)+1
nxf=nxk(this%block%ngrid)
nyf=nyk(this%block%ngrid)
do k=this%block%ngrid-1,1,-1
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
do k=this%block%ngrid-1,1,-1
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
    if (k.eq.(this%block%ngrid-1)) then
      write(6,300) this%block%ngrid,this%block%ngrid-1
      goto 1000
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
    if (k.eq.(this%block%ngrid-1)) then
      write(6,300) this%block%ngrid,this%block%ngrid-1
      goto 1000
    end if
    syr(k)=sym
    eyr(k)=eym
    nyr(k)=nym
  end if
  call grid1_type(irdatatype(k),jrdatatype(k),ijrdatatype(k), &
                  sxr(k),exr(k),syr(k),eyr(k))
end do
endif

return
1000 continue
write(6,200)
call sll_halt_collective()
close(8)
stop

100 format(/,'ERROR: numprocs <> (nxprocs*nyprocs)',/)
110 format(/,'ERROR: process:',i3,' nxdim=',i4,' < ex-sx+3=',i4,/)
120 format(/,'ERROR: process:',i3,' nydim=',i4,' < ey-sy+3=',i4,/) 
130 format(/,'nxp1=',i3,' <> ixp*2**(iex-1)+1=',                       &
          i3,/,'-> adjust the multigrid parameters ixp and iex',/)
140 format(/,'nyp1=',i3,' <> jyq*2**(jey-1)+1=',                       &
          i3,/,'-> adjust the multigrid parameters jyq and jey',/)
150 format(/,'ixp=',i3,' < nxprocs=',i3,/,                             &
            ' there must be at least one grid point at the ',          &
            'coarsest grid level',/,                                   &
            '-> increase ixp and decrease iex correspondingly',/)
160 format(/,'jyq=',i3,' < nyprocs=',i3,/,                             &
            ' there must be at least one grid point at the ',          &
            'coarsest grid level',/,                                   &
            '-> increase jyq and decrease jey correspondingly',/)
170 format(/,'ngrid=',i3,' iex=',i3,                                   &
           /,'no coarsifying at the finest grid level in x-direction', &
           /,'this is not allowed by the mutligrid code',/)
180 format(/,'ngrid=',i3,' jey=',i3,                                   &
           /,'no coarsifying at the finest grid level in y-direction', &
           /,'this is not allowed by the mutligrid code',/)
200 format(/,'ERROR in multigrid code initialization',/)
230 format(/,'ERROR in mgdinit: nx=',i3,' is not a multiple of ', &
           'nxprocs=',i3,/,'cannot use the new version of the ',  &
           'multigrid code',/)
210 format(/,'ERROR in mgdinit: ny=',i3,' is not a multiple of ', &
           'nyprocs=',i3,/,'cannot use the new version of the ',  &
           'multigrid code',/)
220  format(/,'NOTE in mgdinit: ibdry=',i2,' jbdry=',i2, &
            /,'boundary conditions that are not periodic',/)
300 format('no coarsifying between level ', &
           i3,' and level ',i3,/,', the current version of ', &
           'the code cannot cope with that',/, &
           ' -> decrease the value of ixp and/or jyq and',/, &
           '    increase the value of iex and/or jey')


end subroutine initialize_multigrid_2d

!> Parallel multigrid solver in 2-D cartesian coordinates for the
!> elliptic equation:      div(cof(grad(phi)))=rhs
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
!> Notes: - the values of the rhs contained in the array rhs are
!>          transferred to the work vector and the memory of the 
!>          rhs array is then used to store the residuals at the
!>          different levels; therefore, it does not preserve its
!>          initial values
!>        - with the initialization and clean-up routines mgdinit
!>          and mgdend, this multigrid code is self-standing
!>
subroutine solve_multigrid_2d(this, rhs, phi, r)

type(sll_multigrid_solver_2d) :: this
sll_real64, dimension(:,:)    :: rhs
sll_real64, dimension(:,:)    :: phi
sll_real64, dimension(:,:)    :: r

sll_int32  :: iter
sll_real64 :: avo, acorr, relmax
logical    :: nprscr = .true.
sll_int32  :: sx, ex, sy, ey
sll_int32  :: exc, exm, exf, sxc, sxm, sxf, nxc, nxm, nxf
sll_int32  :: eyc, eym, eyf, syc, sym, syf, nyc, nym, nyf
sll_int32  :: ic, icf, ip, ipf, ir1, ir2, irf, itype, jtype
sll_int32  :: lev, kcur
!-----------------------------------------------------------------------
! initialize problem
! xl,yl are the dimensions of the domain
! wk is the wavenumber (must be an integer value)
! rro is the average density
!
r = 1.0


this%xl     = this%x_max-this%x_min
this%yl     = this%y_max-this%y_min

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

sxf=sxk(this%block%ngrid)
exf=exk(this%block%ngrid)
syf=syk(this%block%ngrid)
eyf=eyk(this%block%ngrid)
nxf=nxk(this%block%ngrid)
nyf=nyk(this%block%ngrid)
icf=kcbgn(this%block%ngrid)

call mgdpfpde(sxf,exf,syf,eyf,nxf,nyf,work(icf),r,this%xl,this%yl,this%block%bd)

if (WMGD) then
!
! new version: determine coefficients at coarser grid levels from
! four neighboring density points; no communication of boundary data
! is involved because of the grid setup, under the condition that
! mod(nx,nxprocs)=mod(ny,nyprocs)=0
!
  do k=this%block%ngrid-1,1,-1
    sxm=sxk(k)
    exm=exk(k)
    sym=syk(k)
    eym=eyk(k)
    nxm=nxk(k)
    nym=nyk(k)
    ic=kcbgn(k)
    call mgdphpde(sxm,exm,sym,eym,nxm,nym,work(ic), &
                  sx,ex,sy,ey,nxf,nyf,this%block%bd,this%xl,this%yl)
  end do
else  
  if (this%block%ngrid.gt.1) then
!
! old version: use two locations ir1 and ir2 in the work vector to
! store the density on the different grid levels. First set r in 
! work(ir1) at the finest grid level; this is used to determine the 
! coefficients at level k=ngrid-1
!
    ir1=kpbgn(this%block%ngrid)
    ir2=kcbgn(this%block%ngrid)+5*(exf-sxf+3)*(eyf-syf+3)
    sxf=sxr(this%block%ngrid-1)
    exf=exr(this%block%ngrid-1)
    syf=syr(this%block%ngrid-1)
    eyf=eyr(this%block%ngrid-1)
    nxf=nxr(this%block%ngrid-1)
    nyf=nyr(this%block%ngrid-1)
    lev=1
    call mgdrsetf(sxf,exf,syf,eyf,work(ir1),r)
!
! for the levels k=ngrid-1,1, calculate the coefficients from the
! densities stored in the arrays work(ir1) and work(ir2)
! for the levels k=ngrid-1,2, transfer to the level below the values
! of the density needed there to determine the coefficients; exchange
! of the boundary density data is necessary
!
    do k=this%block%ngrid-1,1,-1
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
                     this%xl,this%yl,this%block%bd)
      else
        call mgdppde(sxm,exm,sym,eym,nxm,nym,work(ic), &
                     sxf,exf,syf,eyf,work(ir2),        &
                     this%xl,this%yl,this%block%bd)
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
                        this%comm2d,this%block%id,          &
                        this%block%neighbor,this%block%bd,     &
                        itype,jtype)
          lev=2
        else
          call mgdrtrsf(sxc,exc,syc,eyc,nxc,nyc,work(ir1), &
                        sxf,exf,syf,eyf,nxf,nyf,work(ir2), &
                        this%comm2d,this%block%id,          &
                        this%block%neighbor,this%block%bd,     &
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
! set phi,rhs in work at the finest grid level
!
sxf = sxk(this%block%ngrid)
exf = exk(this%block%ngrid)
syf = syk(this%block%ngrid)
eyf = eyk(this%block%ngrid)
ipf = kpbgn(this%block%ngrid)
irf = kcbgn(this%block%ngrid)+5*(exf-sxf+3)*(eyf-syf+3)
call mgdsetf(sxf,exf,syf,eyf,work(ipf),work(irf),phi,rhs)
!------------------------------------------------------------------------
! cycling at kcur=ngrid level
!
kcur=this%block%ngrid
do iter=1,this%maxcy
  call mgdkcyc(work,rhs,kcur,this%kcycle,this%iprer,     &
               this%ipost,this%iresw, this%comm2d,       &
               this%block%id,this%block%neighbor,this%block%bd,   &
               this%phibc)
  sxm=sxk(this%block%ngrid)
  exm=exk(this%block%ngrid)
  sym=syk(this%block%ngrid)
  eym=eyk(this%block%ngrid)
  ip=kpbgn(this%block%ngrid)
  call mgderr(relmax,sxm,exm,sym,eym,phi,work(ip),this%comm2d)
  if (relmax.le.this%tolmax) goto 1000
end do
!------------------------------------------------------------------------
! if not converged in maxcy cycles, issue an error message and quit
!
if (this%block%id.eq.0) write(6,100) this%maxcy,relmax
100   format('WARNING: failed to achieve convergence in ',i5, &
       ' cycles  error=',e12.5)
goto 1000
!------------------------------------------------------------------------
! converged
!
1000  continue
!
! rescale phi
!
avo=0.0d0
call gscale(this%block%sx,this%block%ex,this%block%sy,this%block%ey, &
            phi,avo,acorr,this%comm2d,this%nc_x,this%nc_y)

! exchange boundary data and impose periodic BCs

call gxch1lin(phi,this%comm2d,this%block%sx,this%block%ex,     &
              this%block%sy,this%block%ey,this%block%neighbor,     &
              this%block%bd, ikdatatype(this%block%ngrid),       &
              jkdatatype(this%block%ngrid))

call gxch1cor(phi,this%comm2d,sxm,exm,sym,eym,             &
              this%block%neighbor,this%block%bd,                 &
              ijkdatatype(this%block%ngrid))
if (WMGD) then
!
! impose wall and Dirichlet BCs
!
call mgdbdry(sx,ex,sy,ey,phi,this%block%bd,this%phibc)
endif

if (nprscr.and.this%block%id.eq.0) write(6,120) relmax,iter,acorr
120 format('  P MGD     err=',e12.3,' iters=',i5,' pcorr=',e12.3)

if (error == 1) then
   write(6,210)
   call sll_halt_collective()
end if

210 format(/,'ERROR in multigrid code solver',/)

end subroutine solve_multigrid_2d

end module sll_multigrid_2d
