module sll_multigrid_parallel

use sll_collective, only: sll_boot_collective,      &
                          sll_world_collective,     &
                          sll_halt_collective,      &
                          sll_get_collective_rank,  &
                          sll_get_collective_size

use sll_gnuplot_parallel
use sll_remapper
use mgd2

#include "sll_working_precision.h"
#include "sll_constants.h"
#include "sll_memory.h"
#include "mgd2.h"

contains

subroutine initialize_multigrid(layout, &
                                eta1_min, eta1_max, nnodes_eta1, &
                                eta2_min, eta2_max, nnodes_eta2, &
                                geometry)

type(layout_2d), pointer :: layout
sll_real64, intent(in)   :: eta1_min
sll_real64, intent(in)   :: eta1_max
sll_real64, intent(in)   :: eta2_min
sll_real64, intent(in)   :: eta2_max
sll_int32 , intent(in)   :: nnodes_eta1
sll_int32 , intent(in)   :: nnodes_eta2
sll_int32 , intent(in)   :: geometry
sll_int32 , parameter    :: CARTESIAN = 1, POLAR = 2

sll_int32  :: nxp2,nyp2,nx,ny,nxprocs,nyprocs,ixp,jyq,iex,jey
sll_int32  :: maxcy,kcycle,iprer,ipost,iresw,nwork
sll_int32  :: ngrid,nxdim,nydim
sll_real64 :: tolmax

logical           :: nprscr,periods(2)
sll_int32         :: numprocs,myid,comm2d,sx,ex,sy,ey,neighbor(8),bd(8),iter
sll_int32         :: realtype,ibdry,jbdry
sll_int32         :: ierr,nerror,m0,m1,dims(2),coords(2),i,j
sll_int32         :: nbrright,nbrbottom,nbrleft,nbrtop
sll_real64        :: xl,yl,hxi,hyi,vbc(4),phibc(4,20),wk,rro
character(len=1)  :: num(10)
data (num(i),i=1,10)/'0','1','2','3','4','5','6','7','8','9'/

sll_int32 :: error

type(block)     :: my_block
type(mg_solver) :: my_mg

integer :: statut(MPI_STATUS_SIZE)

sll_real64, allocatable  :: p(:,:)
sll_real64, allocatable  :: r(:,:)
sll_real64, allocatable  :: f(:,:)
sll_real64, allocatable  :: work(:)
sll_real64, allocatable  :: xcoord(:,:)
sll_real64, allocatable  :: ycoord(:,:)

nxp2 = nnodes_eta1 + 2
nyp2 = nnodes_eta2 + 2

nx = nxp2 - 2
ny = nyp2 - 2

nxprocs = 2
nyprocs = 2

ixp = 4    ! the resolution of the global domain,
jyq = 4    ! such that :
iex = 5    ! nx = ixp*2**(iex-1)
jey = 5    ! ny = jyq*2**(jey-1)

maxcy  = 300     ! maximum number of multigrid cycle
tolmax = 1.0d-05 ! desired accuracy
kcycle = 1       ! 1 : V-cycles 2 : W-cycles
iprer  = 2       ! relaxation steps before restricting the residual
ipost  = 1       ! relaxation steps after correction
iresw  = 1       ! 1 : fully weighted residuals 2 : half-weighted residuals

nwork  = (4*nx*ny*8)/(3*nxprocs*nyprocs)+(64*(nx+ny))/3+(32*4)/3

ngrid  = max(iex,jey)
nxdim  = int(float(nxp2-2)/float(nxprocs)+0.99)+2
nydim  = int(float(nyp2-2)/float(nyprocs)+0.99)+2

SLL_CLEAR_ALLOCATE(p(1:nxdim,1:nydim), error)
SLL_CLEAR_ALLOCATE(r(1:nxdim,1:nydim), error)
SLL_CLEAR_ALLOCATE(f(1:nxdim,1:nydim), error)
SLL_CLEAR_ALLOCATE(work(1:nwork), error)

numprocs  = sll_get_collective_size(sll_world_collective)
myid      = sll_get_collective_rank(sll_world_collective)
comm2d    = sll_world_collective%comm
realtype  = MPI_REAL8

! open file for output of messages and check that the number of 
! processes is correct

m1 = mod(myid,10)+1
m0 = mod(myid/10,10)+1
if (numprocs.ne.(nxprocs*nyprocs)) then
  write(6,100)
  stop
end if

! create Cartesian topology with periodic boundary conditions

dims(1)    = nxprocs
dims(2)    = nyprocs
periods(1) = .false. !non periodic in r
periods(2) = .true.  !periodic in theta
call MPI_CART_CREATE(MPI_COMM_WORLD,2,dims,periods,.true.,comm2d,ierr)
call MPI_COMM_RANK(comm2d,myid,ierr)

ibdry  = 1    ! Dirchlet bc in direction 1
jbdry  = 0    ! Periodic bc in direction 2

vbc(1) = 0.0d0 ! Value for eta_1 = eta1_min
vbc(2) = 0.0d0 ! Value for eta_1 = eta2_max
vbc(3) = 0.0d0 ! Value for eta_2 = eta2_min
vbc(4) = 0.0d0 ! Value for eta_2 = eta2_max

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

call MPI_CART_SHIFT(comm2d,0,1,nbrleft,nbrright,ierr)
call MPI_CART_SHIFT(comm2d,1,1,nbrbottom,nbrtop,ierr)

neighbor(1)=nbrright
neighbor(3)=nbrbottom
neighbor(5)=nbrleft
neighbor(7)=nbrtop

call MPI_BARRIER(comm2d,ierr)
call MPI_SENDRECV(neighbor(3),1,MPI_INTEGER,nbrright,0, &
                  neighbor(4),1,MPI_INTEGER,nbrleft,0, &
                  comm2d,statut,ierr)
call MPI_SENDRECV(neighbor(7),1,MPI_INTEGER,nbrright,1, &
                  neighbor(6),1,MPI_INTEGER,nbrleft,1, &
                  comm2d,statut,ierr)
call MPI_SENDRECV(neighbor(3),1,MPI_INTEGER,nbrleft,0, &
                  neighbor(2),1,MPI_INTEGER,nbrright,0, &
                  comm2d,statut,ierr)
call MPI_SENDRECV(neighbor(7),1,MPI_INTEGER,nbrleft,1, &
                  neighbor(8),1,MPI_INTEGER,nbrright,1, &
                  comm2d,statut,ierr)
!-----------------------------------------------------------------------
! find indices of subdomain and check that dimensions of arrays are
! sufficient
!
call MPI_CART_GET(comm2d,2,dims,periods,coords,ierr)
call MPE_DECOMP1D(nxp2-2,dims(1),coords(1),sx,ex)
sx=sx+1
ex=ex+1
call MPE_DECOMP1D(nyp2-2,dims(2),coords(2),sy,ey)
sy=sy+1
ey=ey+1
if ((ex-sx+3) > nxdim) then
  write(6,110) myid,nxdim,ex-sx+3
  nerror=1
  stop
end if
if ((ey-sy+3).gt.nydim) then
  write(6,120) myid,nydim,ey-sy+3
  nerror=1
  stop
end if
write(6,*) 'sx=',sx,' ex=',ex,' sy=',sy,' ey=',ey
do i=1,8
  write(6,*) 'neighbor: ',neighbor(i),' bd: ',bd(i)
end do

my_block%id       = myid
my_block%sx       = sx
my_block%ex       = ex
my_block%sy       = sy
my_block%ey       = ey
my_block%ngrid    = ngrid
my_block%neighbor = neighbor
my_block%bd       = bd

my_block%ixp      = ixp
my_block%jyq      = jyq
my_block%iex      = iex
my_block%jey      = jey

my_mg%vbc         = vbc
my_mg%phibc       = phibc
my_mg%nx          = nx
my_mg%ny          = ny
my_mg%nxprocs     = nxprocs
my_mg%nyprocs     = nyprocs
my_mg%ibdry       = ibdry
my_mg%jbdry       = jbdry
my_mg%comm2d      = comm2d

if (geometry == CARTESIAN) then
   call initialize_mgd2(my_block,my_mg,nerror)
else if (geometry == POLAR) then
   call initialize_mgd_polar(my_block,my_mg,nerror)
end if
if (nerror.eq.1) goto 1000
!-----------------------------------------------------------------------
! initialize problem
! xl,yl are the dimensions of the domain
! wk is the wavenumber (must be an integer value)
! rro is the average density
! 1/hxi,1/hyi are the spatial resolutions
!
xl=1.0d0
yl=1.0d0
wk=2.0d0
rro=1.0d0
hxi=float(nx)/xl
hyi=float(ny)/yl
write(6,*) 'hxi=',hxi,' hyi=',hyi
call ginit(sx,ex,sy,ey,p,r,f,wk,hxi,hyi)

my_mg%xl     = xl
my_mg%yl     = yl
my_mg%tolmax = tolmax
my_mg%kcycle = kcycle
my_mg%maxcy  = maxcy
my_mg%iprer  = iprer
my_mg%ipost  = ipost
my_mg%iresw  = iresw
my_mg%isol   = 2

if (geometry == CARTESIAN) then
   call mgd2_solver(my_block,my_mg,p,f,r,work, &
                  rro, iter,.true.,nerror)
   if (nerror.eq.1) goto 1000
   ! compare numerical and exact solutions
   call gerr(sx,ex,sy,ey,p,comm2d,wk,hxi,hyi,nx,ny)
   
   call sll_gnuplot_rect_2d_parallel(dble(sx-2.0)/hxi, dble(1)/hxi, &
                                     dble(sy-2.0)/hyi, dble(1)/hyi, &
                                     p, "potential", 1, error)  
  
else if (geometry == POLAR) then
   call mgd_polar_solver(my_block,my_mg,p,f,r,work, &
                  rro, iter,.true.,nerror)

   if (nerror.eq.1) goto 1000
   ! compare numerical and exact solutions
   call gerr(sx,ex,sy,ey,p,comm2d,wk,hxi,hyi,nx,ny)
   
   call sll_gnuplot_curv_2d_parallel(xcoord, ycoord, &
                                     p, "potential", 1, error)  
  
end if

return
call sll_halt_collective()
1000  write(6,200)
close(8)
stop

100 format(/,'ERROR: numprocs <> (nxprocs*nyprocs)',/)
110 format(/,'ERROR: process:',i3,' nxdim=',i4,' < ex-sx+3=',i4,/,  &
       ' -> put the parameter formula for nxdim in main.F in ',     &
            'comments and',/,'    assign to nxdim the maximum ',    &
            'value of ex-sx+3',/)
120 format(/,'ERROR: process:',i3,' nydim=',i4,' < ey-sy+3=',i4,/,  &
           ' -> put the parameter formula for nydim in main.F in ', &
           'comments and'/,'     assign to nydim the maximum ',     &
           'value of ey-sy+3',/)
200 format(/,'ERROR in multigrid code',/)

end subroutine initialize_multigrid

!> Initialize the pressure, density, and right-hand side of the
!> elliptic equation div(1/r*grad(p))=f
subroutine ginit(sx,ex,sy,ey,p,r,f,wk,hxi,hyi)

sll_int32  :: sx,ex,sy,ey
sll_real64 :: p(sx-1:ex+1,sy-1:ey+1),r(sx-1:ex+1,sy-1:ey+1)
sll_real64 :: f(sx-1:ex+1,sy-1:ey+1),hxi,hyi,wk
sll_int32  :: i,j
sll_real64 :: cnst,cx,cy,xi,yj

p = 0.0d0
r = 1.0d0
f = 0.0d0

cnst = -8.0d0*(sll_pi*wk)**2
cx   = 2.0d0*sll_pi*wk
cy   = 2.0d0*sll_pi*wk

do j=sy,ey
  do i=sx,ex
    xi=(i-1.5)/hxi
    yj=(j-1.5)/hyi
    f(i,j)=cnst*sin(cx*xi)*sin(cy*yj)
  end do
end do

end subroutine

!> Calculate the error between the numerical and exact solution to
!> the test problem.
subroutine gerr(sx,ex,sy,ey,p,comm2d,wk,hxi,hyi,nx,ny)
sll_int32  :: sx,ex,sy,ey,comm2d,nx,ny
sll_real64 :: p(sx-1:ex+1,sy-1:ey+1),wk,hxi,hyi,xi,yj
sll_int32  :: i,j,ierr
sll_real64 :: errloc,err,cx,cy,exact

! calculate local error
cx=2.0d0*sll_pi*wk
cy=2.0d0*sll_pi*wk
errloc=0.0d0
do j=sy,ey
  yj=(j-1.5)/hyi
  do i=sx,ex
    xi=(i-1.5)/hxi
    exact=sin(cx*xi)*sin(cy*yj)
    errloc=errloc+abs(p(i,j)-exact)
  end do
end do

! calculate global error
call MPI_ALLREDUCE(errloc,err,1,MPI_REAL8,MPI_SUM,comm2d,ierr)
write(6,100) errloc/float(nx*ny),err/float(nx*ny)
100   format(/,'Local error: ',e13.6,'  total error: ',e13.6,/)

end subroutine

end module sll_multigrid_parallel
