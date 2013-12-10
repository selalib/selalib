!-----------------------------------------------------------------------
! Test problem for 3D multigrid parallel code mgd3. mgd3 solves the
! non-separable elliptic equation div(1/r*grad(p))=f on a rectangular
! domain with a staggered uniform grid and either periodic or Neumann
! (set normal derivative of p at boundary to 0) or Dirichlet (set value
! of p at boundary) boundary conditions.
!
! Here, r is constant and equal to 1 everywhere. The right-hand side
! f is chosen so that it corresponds to 
!                p=sin(2*pi*wk*x)*sin(2*pi*wk*y)*sin(2*pi*wk*z),
! where k is the wavenumber. The program computes the numerical 
! solution and compares it to the exact expression.
!
! The program also shows how to implement mgd3. The user must define
! a few parameters describing the resolution and the number of levels
! to operate on, four communicators and a datatype.
!
! An important feature is that before doing any calculation using
! mgdsolver, the user must call mgdinit to initialize the multigrid
! parameters. Once this is done, mgdsolver can be called any number
! of times, the variables are not overwritten.
!
! Input     : none
! Outputs   : messages -> out* files
! Code      : tmgd3, 3-D parallel multigrid solver
! Calls     : MPI_INIT, MPI_COMM_RANK, MPI_COMM_SIZE, MPI_CART_CREATE,
!             MPI_COMM_RANK, MPI_CART_SHIFT, MPI_BARRIER, MPI_SENDRECV,
!             MPI_CART_GET, MPE_DECOMP1D, 
!             mgdinit, ginit, mgdsolver, gerr
!-----------------------------------------------------------------------
program test_multigrid_3d
use hdf5, only:HID_T,HSIZE_T,HSSIZE_T
use mgd3
use sll_collective, only: sll_boot_collective,      &
                          sll_world_collective,     &
                          sll_halt_collective,      &
                          sll_get_collective_rank,  &
                          sll_get_collective_size,  &
                          sll_collective_barrier
use sll_hdf5_io_parallel
use sll_xml_io
use sll_constants
#include "mgd3.h"
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_utilities.h"
implicit none

logical :: periods(3)
sll_int32 :: nxdim, nydim, nzdim
sll_int32 :: numprocs, comm3d, comm3dp, comm3dl, comm3dc
sll_int32 :: sx, ex, sy, ey, sz, ez, neighbor(26), bd(26)
sll_int32 :: ngb(3), myid
sll_int32 :: ierr, nerror, dims(3), coords(3)

sll_real64, allocatable :: p(:,:,:), f(:,:,:), r(:,:,:)

sll_real64 :: hxi, hyi, hzi, wk

type(block) :: my_block
type(mg_solver) :: my_mg

sll_int32   :: error
sll_real64  :: tcpu1
sll_real64  :: tcpu2

integer(HID_T)                  :: file_id
integer(HSIZE_T), dimension(3)  :: datadims
integer(HSSIZE_T), dimension(3) :: offset 
sll_real64, dimension(:,:,:), allocatable :: x, y, z
sll_int32 :: i, j, k

! Boot parallel environment
call sll_boot_collective()

numprocs  = sll_get_collective_size(sll_world_collective)
myid      = sll_get_collective_rank(sll_world_collective)
comm3d    = sll_world_collective%comm

if( myid .eq. 0) print *, ' '

!Mutigrid solver parameters

my_mg%ibdry=0
my_mg%jbdry=0
my_mg%kbdry=0

my_mg%nxprocs=2
my_mg%nyprocs=2
my_mg%nzprocs=1

my_mg%tolmax = 1e-6

my_mg%maxcy=300
my_mg%kcycle=1
my_mg%iprer=2
my_mg%ipost=1
my_mg%iresw=1

my_mg%nx=32
my_mg%ny=32
my_mg%nz=32

my_block%ixp=4
my_block%jyq=4
my_block%kzr=4

my_block%iex=4
my_block%jey=4
my_block%kez=4

my_block%ngrid=max(my_block%iex,my_block%jey,my_block%kez)

my_block%id = myid


if (numprocs.ne.(my_mg%nxprocs*my_mg%nyprocs*my_mg%nzprocs)) then
  write(6,100)
100     format(/,'ERROR: numprocs <> (nxprocs*nyprocs*nzprocs)',/)
  stop
end if



!-----------------------------------------------------------------------
! create Cartesian topology with periodic boundary conditions
!
dims(1)=my_mg%nxprocs
dims(2)=my_mg%nyprocs
dims(3)=my_mg%nzprocs

periods(1)=.true.
periods(2)=.true.
periods(3)=.true.

call MPI_CART_CREATE(MPI_COMM_WORLD,3,dims,periods,.true.,comm3d,ierr)
call MPI_COMM_DUP(comm3d,comm3dp,ierr)
call MPI_COMM_DUP(comm3d,comm3dl,ierr)
call MPI_COMM_DUP(comm3d,comm3dc,ierr)
call MPI_COMM_RANK(comm3d,my_block%id,ierr)

my_mg%vbc(:)=0.0d0
bd(:)=0
!-----------------------------------------------------------------------
! find the neighbors
! convention
!
! plane i
!      6|     7     |8
!      ---------------                k
!       |           |               ^ 
!       |           |               |
!      5|    myid   |1              |
!       |           |              ------>
!       |           |               |     j
!      ---------------
!      4|     3     |2
!
! plane i-1                  plane i+1
!     15|    16     |17          24|     25    |26 
!      ---------------            --------------- 
!       |           |              |           |
!       |           |              |           |
!     14|     9     |10          23|     18    |19
!       |           |              |           |
!       |           |              |           |
!      ---------------            ---------------
!     13|    12     |11          22|     21    |20
!
call MPI_CART_GET(comm3d,3,dims,periods,coords,ierr)

!
! neighbors which share a plane with 'myid'
!
ngb=(/coords(1),coords(2)+1,coords(3)/)
call MPI_CART_RANK(comm3d,ngb,neighbor(1),ierr)
ngb=(/coords(1),coords(2),coords(3)-1/)
call MPI_CART_RANK(comm3d,ngb,neighbor(3),ierr)
ngb=(/coords(1),coords(2)-1,coords(3)/)
call MPI_CART_RANK(comm3d,ngb,neighbor(5),ierr)
ngb=(/coords(1),coords(2),coords(3)+1/)
call MPI_CART_RANK(comm3d,ngb,neighbor(7),ierr)
ngb=(/coords(1)-1,coords(2),coords(3)/)
call MPI_CART_RANK(comm3d,ngb,neighbor(9),ierr)
ngb=(/coords(1)+1,coords(2),coords(3)/)
call MPI_CART_RANK(comm3d,ngb,neighbor(18),ierr)
!
! neighbors which share a line with 'myid'
!
ngb=(/coords(1),coords(2)+1,coords(3)-1/)
call MPI_CART_RANK(comm3d,ngb,neighbor(2),ierr)
ngb=(/coords(1),coords(2)-1,coords(3)-1/)
call MPI_CART_RANK(comm3d,ngb,neighbor(4),ierr)
ngb=(/coords(1),coords(2)-1,coords(3)+1/)
call MPI_CART_RANK(comm3d,ngb,neighbor(6),ierr)
ngb=(/coords(1),coords(2)+1,coords(3)+1/)
call MPI_CART_RANK(comm3d,ngb,neighbor(8),ierr)
ngb=(/coords(1)-1,coords(2)+1,coords(3)/)
call MPI_CART_RANK(comm3d,ngb,neighbor(10),ierr)
ngb=(/coords(1)-1,coords(2),coords(3)-1/)
call MPI_CART_RANK(comm3d,ngb,neighbor(12),ierr)
ngb=(/coords(1)-1,coords(2)-1,coords(3)/)
call MPI_CART_RANK(comm3d,ngb,neighbor(14),ierr)
ngb=(/coords(1)-1,coords(2),coords(3)+1/)
call MPI_CART_RANK(comm3d,ngb,neighbor(16),ierr)
ngb=(/coords(1)+1,coords(2)+1,coords(3)/)
call MPI_CART_RANK(comm3d,ngb,neighbor(19),ierr)
ngb=(/coords(1)+1,coords(2),coords(3)-1/)
call MPI_CART_RANK(comm3d,ngb,neighbor(21),ierr)
ngb=(/coords(1)+1,coords(2)-1,coords(3)/)
call MPI_CART_RANK(comm3d,ngb,neighbor(23),ierr)
ngb=(/coords(1)+1,coords(2),coords(3)+1/)
call MPI_CART_RANK(comm3d,ngb,neighbor(25),ierr)
!
! neighbors which share a corner with 'myid'
!
ngb=(/coords(1)-1,coords(2)+1,coords(3)-1/)
call MPI_CART_RANK(comm3d,ngb,neighbor(11),ierr)
ngb=(/coords(1)-1,coords(2)-1,coords(3)-1/)
call MPI_CART_RANK(comm3d,ngb,neighbor(13),ierr)
ngb=(/coords(1)-1,coords(2)-1,coords(3)+1/)
call MPI_CART_RANK(comm3d,ngb,neighbor(15),ierr)
ngb=(/coords(1)-1,coords(2)+1,coords(3)+1/)
call MPI_CART_RANK(comm3d,ngb,neighbor(17),ierr)
ngb=(/coords(1)+1,coords(2)+1,coords(3)-1/)
call MPI_CART_RANK(comm3d,ngb,neighbor(20),ierr)
ngb=(/coords(1)+1,coords(2)-1,coords(3)-1/)
call MPI_CART_RANK(comm3d,ngb,neighbor(22),ierr)
ngb=(/coords(1)+1,coords(2)-1,coords(3)+1/)
call MPI_CART_RANK(comm3d,ngb,neighbor(24),ierr)
ngb=(/coords(1)+1,coords(2)+1,coords(3)+1/)
call MPI_CART_RANK(comm3d,ngb,neighbor(26),ierr)
!-----------------------------------------------------------------------
! find indices of subdomain and check that dimensions of arrays are
! sufficient
!
call MPE_DECOMP1D(my_mg%nx,dims(1),coords(1),sx,ex)
sx=sx+1
ex=ex+1
call MPE_DECOMP1D(my_mg%ny,dims(2),coords(2),sy,ey)
sy=sy+1
ey=ey+1
call MPE_DECOMP1D(my_mg%nz,dims(3),coords(3),sz,ez)
sz=sz+1
ez=ez+1

nxdim=int(float(my_mg%nx)/float(my_mg%nxprocs)+0.99)+2
nydim=int(float(my_mg%ny)/float(my_mg%nyprocs)+0.99)+2
nzdim=int(float(my_mg%nz)/float(my_mg%nzprocs)+0.99)+2

if ((ex-sx+3).gt.nxdim) then
  write(6,110) my_block%id,nxdim,ex-sx+3
  nerror=1
  stop
end if
if ((ey-sy+3).gt.nydim) then
  write(6,120) my_block%id,nydim,ey-sy+3
  nerror=1
  stop
end if
if ((ez-sz+3).gt.nzdim) then
  write(6,130) my_block%id,nzdim,ez-sz+3
  nerror=1
  stop
end if
!write(6,*) 'sx=',sx,' ex=',ex,' sy=',sy,' ey=',ey
!do i=1,26
!  write(6,*) 'neighbor: ',neighbor(i),' bd: ',bd(i)
!end do

!-----------------------------------------------------------------------
! initialize mgd3
!-----------------------------------------------------------------------

my_block%sx = sx; my_block%ex = ex
my_block%sy = sy; my_block%ey = ey
my_block%sz = sz; my_block%ez = ez

allocate(p(sx-1:ex+1,sy-1:ey+1,sz-1:ez+1))
allocate(f(sx-1:ex+1,sy-1:ey+1,sz-1:ez+1))
allocate(r(sx-1:ex+1,sy-1:ey+1,sz-1:ez+1))

call initialize_mgd3(my_block,my_mg,nerror)

if (nerror.eq.1) goto 1000
!-----------------------------------------------------------------------
! initialize problem
! xl,yl,zl are the dimensions of the domain
! wk is the wavenumber (must be an integer value)
! rro is the average density
! 1/hxi,1/hyi,1/hzi are the spatial resolutions
!-----------------------------------------------------------------------

my_mg%xl=1.0d0
my_mg%yl=1.0d0
my_mg%zl=1.0d0
wk=2.0d0
hxi=float(my_mg%nx)/my_mg%xl
hyi=float(my_mg%ny)/my_mg%yl
hzi=float(my_mg%nz)/my_mg%zl

call ginit()

call write_xdmf_3d(myid,numprocs,f,sx,ex,sy,ey,sz,ez,hxi,hyi,hzi,nerror)

SLL_ALLOCATE(x(sx:ex,sy:ey,sz:ez),error)
SLL_ALLOCATE(y(sx:ex,sy:ey,sz:ez),error)
SLL_ALLOCATE(z(sx:ex,sy:ey,sz:ez),error)

tcpu1 = MPI_WTIME()

! initialize the local data    
do k=sz,ez
   do j=sy,ey
      do i=sx,ex
         x(i,j,k)=(float(i))/sngl(hxi)
         y(i,j,k)=(float(j))/sngl(hyi)
         z(i,j,k)=(float(k))/sngl(hzi)
      end do
   end do
end do

datadims = (/my_mg%nx,my_mg%ny,my_mg%nz/)
offset   = (/sx-2,sy-2,sz-2/)

call sll_hdf5_file_create('grid3d.h5',file_id, error)
call sll_hdf5_write_array(file_id,datadims,offset,x,'x',error)
call sll_hdf5_write_array(file_id,datadims,offset,y,'y',error)
call sll_hdf5_write_array(file_id,datadims,offset,z,'z',error)
call sll_hdf5_write_array(file_id,datadims,offset,f(sx:ex,sy:ey,sz:ez),'array',error)
call sll_hdf5_file_close(file_id, error)

if (myid == 0) then
   call sll_xml_file_create("grid3d.xmf",file_id,error)
   call sll_xml_grid_geometry(file_id, 'grid3d.h5', my_mg%nx, &
                                       'grid3d.h5', my_mg%ny, &
                                       'grid3d.h5', my_mg%nz, &
                                       'x', 'y', 'z' )
   call sll_xml_field(file_id,'values', "grid3d.h5:/array", &
                      my_mg%nx,my_mg%ny,my_mg%nz,'HDF','Node')
   call sll_xml_file_close(file_id,error)
end if

call sll_collective_barrier(sll_world_collective)
tcpu2 = MPI_WTIME()
if (myid == 0) &
   write(*,"(//10x,' Temps CPU = ', G15.3, ' sec' )") (tcpu2-tcpu1)*numprocs
!-----------------------------------------------------------------------
! solve using mgd3
!-----------------------------------------------------------------------

my_mg%comm3d  = comm3d
my_mg%comm3dp = comm3dp
my_mg%comm3dl = comm3dl
my_mg%comm3dc = comm3dc

my_block%neighbor = neighbor
my_block%bd = bd

my_mg%isol = 2

call solve(p,f,r,my_block,my_mg,nerror)

call write_xdmf_3d(myid,numprocs,p,sx,ex,sy,ey,sz,ez,hxi,hyi,hzi,nerror)

if (nerror.eq.1) goto 1000

!-----------------------------------------------------------------------
! compare numerical and exact solutions
!-----------------------------------------------------------------------

call gerr()
!-----------------------------------------------------------------------
close(8)

if (myid==0) print*,"PASSED"

call sll_halt_collective()
print*,"PASSED"

stop
1000  write(6,200)
200   format(/,'ERROR in multigrid code',/)
close(8)
stop

110 format(/,'ERROR: process:',i3,' nxdim=',i4,' < ex-sx+3=',i4,/, &
    ' -> put the parameter formula for nxdim in main.F in ',       &
    'comments and',/,'    assign to nxdim the maximum ',           &
    'value of ex-sx+3',/)
120 format(/,'ERROR: process:',i3,' nydim=',i4,' < ey-sy+3=',i4,/, &
    ' -> put the parameter formula for nydim in main.F in ',       &
    'comments and'/,'     assign to nydim the maximum ',           &
    'value of ey-sy+3',/)
130 format(/,'ERROR: process:',i3,' nzdim=',i4,' < ez-sz+3=',i4,/, &
    ' -> put the parameter formula for nzdim in main.F in ',       &
    'comments and'/,'     assign to nzdim the maximum ',           &
    'value of ez-sz+3',/)

contains

subroutine ginit()
implicit none

!-----------------------------------------------------------------------
! Initialize the pressure, density, and right-hand side of the
! elliptic equation div(1/r*grad(p))=f
!
! Code      : tmgd3, test program for 3D parallel multigrid solver
! Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
! Called in : main
! Calls     : --
!-----------------------------------------------------------------------

sll_int32  :: i,j,k
sll_real64 :: cnst,cx,cy,cz,xi,yj,zk

p = 0.0d0
r = 1.0d0
f = 0.0d0

cnst=-12.0d0*(sll_pi*wk)**2
cx=2.0d0*sll_pi*wk
cy=2.0d0*sll_pi*wk
cz=2.0d0*sll_pi*wk
do k=sz,ez
  zk=(float(k)-1.5d0)/hzi
  do j=sy,ey
    yj=(float(j)-1.5d0)/hyi
    do i=sx,ex
      xi=(float(i)-1.5d0)/hxi
      f(i,j,k)=cnst*sin(cx*xi)*sin(cy*yj)*sin(cz*zk)
    end do
  end do
end do

return
end subroutine ginit

subroutine gerr()
implicit none
!-----------------------------------------------------------------------
! Calculate the error between the numerical and exact solution to
! the test problem.
!
! Code      : tmgd3, test program for 3D parallel multigrid solver
! Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
! Called in : main
! Calls     : MPI_ALLREDUCE
!-----------------------------------------------------------------------
sll_int32  :: i,j,k,ierr
sll_real64 :: errloc,err,cx,cy,cz,exact,zk,yj,xi
!
! calculate local error
!
cx=2.0d0*sll_pi*wk
cy=2.0d0*sll_pi*wk
cz=2.0d0*sll_pi*wk
errloc=0.0d0
do k=sz,ez
  zk=(float(k)-1.5d0)/hzi
  do j=sy,ey
    yj=(float(j)-1.5d0)/hyi
    do i=sx,ex
      xi=(float(i)-1.5d0)/hxi
      exact=sin(cx*xi)*sin(cy*yj)*sin(cz*zk)
      errloc=errloc+abs(p(i,j,k)-exact)
    end do
  end do
end do
!
! calculate global error
!
call MPI_ALLREDUCE(errloc,err,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm3d,ierr)
write(6,100) my_block%id,errloc/float(my_mg%nx*my_mg%ny*my_mg%nz),err/float(my_mg%nx*my_mg%ny*my_mg%nz)
100   format(/,i3,3x,'Local error: ',e13.6,'  total error: ',e13.6,/)

return
end subroutine gerr

end program test_multigrid_3d
