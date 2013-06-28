module sll_multigrid_2d

use sll_multigrid_base
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

implicit none

private :: gerr, ginit

type, public, extends(sll_multigrid_solver) :: sll_multigrid_solver_2d

   sll_real64 :: eta1_min
   sll_real64 :: eta1_max
   sll_real64 :: eta2_min
   sll_real64 :: eta2_max
   
   sll_int32 :: nc_eta1
   sll_int32 :: nc_eta2
   sll_int32 :: nprocs_1
   sll_int32 :: nprocs_2

end type sll_multigrid_solver_2d

interface initialize
   module procedure initialize_2d
end interface initialize

sll_int32, private :: error

contains

subroutine initialize_2d(this, eta1_min, eta1_max, nc_eta1, &
                               eta2_min, eta2_max, nc_eta2)

type(sll_multigrid_solver_2d) :: this
sll_real64, intent(in)        :: eta1_min
sll_real64, intent(in)        :: eta1_max
sll_real64, intent(in)        :: eta2_min
sll_real64, intent(in)        :: eta2_max
sll_int32 , intent(in)        :: nc_eta1
sll_int32 , intent(in)        :: nc_eta2

sll_int32  :: nxp2,nyp2,nx,ny,nxprocs,nyprocs,ixp,jyq,iex,jey
sll_int32  :: nwork
sll_int32  :: ngrid,nxdim,nydim
sll_real64 :: tolmax

logical           :: nprscr,periods(2)
sll_int32         :: numprocs,myid,sx,ex,sy,ey,neighbor(8),bd(8),iter
sll_int32         :: realtype
sll_int32         :: ierr,m0,m1,dims(2),coords(2),i,j
sll_int32         :: nbrright,nbrbottom,nbrleft,nbrtop
sll_real64        :: xl,yl,hxi,hyi,wk,rro
character(len=1)  :: num(10)

data (num(i),i=1,10)/'0','1','2','3','4','5','6','7','8','9'/


type(block)     :: my_block
type(mg_solver) :: my_mg

integer :: statut(MPI_STATUS_SIZE)

sll_real64, allocatable :: p(:,:)
sll_real64, allocatable :: r(:,:)
sll_real64, allocatable :: f(:,:)
sll_real64, allocatable :: xcoord(:,:)
sll_real64, allocatable :: ycoord(:,:)

numprocs    = sll_get_collective_size(sll_world_collective)
myid        = sll_get_collective_rank(sll_world_collective)
this%comm2d = sll_world_collective%comm
realtype    = MPI_REAL8

nxprocs = int(sqrt(real(numprocs,f64)))
nyprocs = nxprocs
this%nprocs_1 = nxprocs
this%nprocs_2 = nyprocs

if (nxprocs*nyprocs /= numprocs) goto 1000

nx = nc_eta1
ny = nc_eta2

nxp2 = nx+2
nyp2 = nx+2

! Layout and local sizes for FFTs in x-direction
this%layout => new_layout_2D( sll_world_collective )
call initialize_layout_with_distributed_2D_array( nc_eta1, nc_eta2, &
       nxprocs, nyprocs, this%layout)

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
nxdim  = int(float(nxp2-2)/float(nxprocs)+0.99)+2
nydim  = int(float(nyp2-2)/float(nyprocs)+0.99)+2

if ((ex-sx+3) > nxdim) then
  write(6,110) myid,nxdim,ex-sx+3
  goto 1000
end if
if ((ey-sy+3).gt.nydim) then
  write(6,120) myid,nydim,ey-sy+3
  goto 1000
end if
write(6,*) 'sx=',sx,' ex=',ex,' sy=',sy,' ey=',ey

SLL_CLEAR_ALLOCATE(p(1:nxdim,1:nydim), error)
SLL_CLEAR_ALLOCATE(r(1:nxdim,1:nydim), error)
SLL_CLEAR_ALLOCATE(f(1:nxdim,1:nydim), error)
SLL_CLEAR_ALLOCATE(this%work(1:nwork), error)


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
periods(1) = .true. !non periodic in r
periods(2) = .true.  !periodic in theta
call MPI_CART_CREATE(MPI_COMM_WORLD,2,dims,periods,.true.,this%comm2d,ierr)
call MPI_COMM_RANK(this%comm2d,myid,ierr)

this%ibdry  = 0    ! Periodic bc in direction 1
this%jbdry  = 0    ! Periodic bc in direction 2

this%vbc(1) = 0.0d0 ! Value for eta_1 = eta1_min
this%vbc(2) = 0.0d0 ! Value for eta_1 = eta2_max
this%vbc(3) = 0.0d0 ! Value for eta_2 = eta2_min
this%vbc(4) = 0.0d0 ! Value for eta_2 = eta2_max

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

my_mg%vbc         = this%vbc
my_mg%phibc       = this%phibc
my_mg%nx          = nx
my_mg%ny          = ny
my_mg%nxprocs     = nxprocs
my_mg%nyprocs     = nyprocs
my_mg%ibdry       = this%ibdry
my_mg%jbdry       = this%jbdry
my_mg%comm2d      = this%comm2d

call initialize_mgd2(my_block,my_mg,error)
if (error == 1) write(6,200)
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
my_mg%tolmax = this%tolmax
my_mg%kcycle = this%kcycle
my_mg%maxcy  = this%maxcy
my_mg%iprer  = this%iprer
my_mg%ipost  = this%ipost
my_mg%iresw  = this%iresw
my_mg%isol   = 2

call mgd2_solver(my_block,my_mg,p,f,r,this%work, &
                 rro, iter,.true.,error)
if (error == 1) write(6,210)
! compare numerical and exact solutions
call gerr(sx,ex,sy,ey,p,this%comm2d,wk,hxi,hyi,nx,ny)

call sll_gnuplot_rect_2d_parallel(dble(sx-2.0)/hxi, dble(1)/hxi, &
                                  dble(sy-2.0)/hyi, dble(1)/hyi, &
                                  p, "potential", 1, error)  
  
return
1000 continue
call sll_halt_collective()
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
200 format(/,'ERROR in multigrid code initialization',/)
210 format(/,'ERROR in multigrid code solver',/)

end subroutine initialize_2d

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

end module sll_multigrid_2d
