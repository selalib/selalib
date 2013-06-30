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

implicit none

type, public, extends(sll_multigrid_solver) :: sll_multigrid_solver_2d

   sll_real64 :: x_min
   sll_real64 :: x_max
   sll_real64 :: y_min
   sll_real64 :: y_max
   
   sll_int32 :: nc_x
   sll_int32 :: nc_y
   sll_int32 :: nprocs_x
   sll_int32 :: nprocs_y

   type(block)     :: block
   type(mg_solver) :: mg

end type sll_multigrid_solver_2d

interface initialize
   module procedure initialize_2d
end interface initialize

sll_int32, private :: i, j, k
sll_int32, private :: error

contains

subroutine initialize_2d(this, layout,                &
                         x_min, x_max, nc_x, &
                         y_min, y_max, nc_y)

type(sll_multigrid_solver_2d) :: this
type(layout_2d), pointer      :: layout
sll_real64, intent(in)        :: x_min
sll_real64, intent(in)        :: x_max
sll_real64, intent(in)        :: y_min
sll_real64, intent(in)        :: y_max
sll_int32 , intent(in)        :: nc_x
sll_int32 , intent(in)        :: nc_y

sll_int32  :: nxp2,nyp2,nx,ny,nxprocs,nyprocs,ixp,jyq,iex,jey
sll_int32  :: nwork
sll_int32  :: ngrid,nxdim,nydim

logical   :: periods(2)
sll_int32 :: numprocs,myid,sx,ex,sy,ey,neighbor(8),bd(8)
sll_int32 :: realtype
sll_int32 :: ierr,m0,m1,dims(2)
sll_int32 :: nbrright,nbrbottom,nbrleft,nbrtop

integer :: statut(MPI_STATUS_SIZE)

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

nxp2 = nx+2
nyp2 = nx+2

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

SLL_CLEAR_ALLOCATE(this%work(1:nwork), error)

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

this%ibdry  = 0    ! Periodic bc in direction 1
this%jbdry  = 0    ! Periodic bc in direction 2

this%vbc(1) = 0.0d0 ! Value for eta_1 = x_min
this%vbc(2) = 0.0d0 ! Value for eta_1 = y_max
this%vbc(3) = 0.0d0 ! Value for eta_2 = y_min
this%vbc(4) = 0.0d0 ! Value for eta_2 = y_max

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

this%mg%vbc         = this%vbc
this%mg%phibc       = this%phibc
this%mg%nx          = nx
this%mg%ny          = ny
this%mg%nxprocs     = nxprocs
this%mg%nyprocs     = nyprocs
this%mg%ibdry       = this%ibdry
this%mg%jbdry       = this%jbdry
this%mg%comm2d      = this%comm2d
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

call initialize_mgd2(this%block,this%mg,error)
if (error == 1) then
   goto 1000
end if

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

end subroutine initialize_2d

subroutine solve_2d(this, f, p, r)
sll_real64, dimension(:,:) :: f
sll_real64, dimension(:,:) :: p
sll_real64, dimension(:,:) :: r
type(sll_multigrid_solver_2d) :: this
sll_int32  :: iter
!-----------------------------------------------------------------------
! initialize problem
! xl,yl are the dimensions of the domain
! wk is the wavenumber (must be an integer value)
! rro is the average density
!
r = 1.0

this%mg%xl     = this%x_max-this%x_min
this%mg%yl     = this%y_max-this%y_min
this%mg%tolmax = this%tolmax
this%mg%kcycle = this%kcycle
this%mg%maxcy  = this%maxcy
this%mg%iprer  = this%iprer
this%mg%ipost  = this%ipost
this%mg%iresw  = this%iresw
this%mg%isol   = 2

call mgd2_solver(this%block,this%mg,p,f,r,this%work,iter,.true.,error)

if (error == 1) then
   write(6,210)
   call sll_halt_collective()
end if

210 format(/,'ERROR in multigrid code solver',/)

end subroutine solve_2d

end module sll_multigrid_2d
