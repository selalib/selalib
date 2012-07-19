program test_mgd3
use mgd3
implicit none
#include "mgd3.h"
include "mpif.h"
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
!
! parameters
!
integer, parameter :: maxcy=300
integer, parameter :: kcycle=1, iprer=2, ipost=1, iresw=1
integer, parameter :: iout=60
real(8), parameter :: tolmax=1.0d-05
!
! variables
!
logical :: periods(3)
integer :: nwork,nxdim,nydim,nzdim
integer :: numprocs,comm3d,comm3dp,comm3dl,comm3dc
integer :: sx,ex,sy,ey,sz,ez,neighbor(26),bd(26),iter
integer :: ngb(3)
integer :: ierr,nerror,m0,m1,dims(3),coords(3),i

real(8), allocatable :: p(:,:,:), f(:,:,:), r(:,:,:)
real(8), allocatable :: work(:)

real(8) :: xl,yl,zl,hxi,hyi,hzi,wk,rro,pi
character(len=20) :: outfile
character(len=1)  :: num(10)

type(block) :: my_block
type(mg_solver) :: my_mg

data (num(i),i=1,10)/'0','1','2','3','4','5','6','7','8','9'/

pi=4.0d0*atan(1.0d0)

!Mutigrid solver parameters

my_mg%ibdry=0
my_mg%jbdry=0
my_mg%kbdry=0

my_mg%nxprocs=2
my_mg%nyprocs=2
my_mg%nzprocs=1

my_block%nx=32
my_block%ny=32
my_block%nz=32

my_block%ixp=4
my_block%jyq=4
my_block%kzr=4

my_block%iex=4
my_block%jey=4
my_block%kez=4

my_block%ngrid=max(my_block%iex,my_block%jey,my_block%kez)

!-----------------------------------------------------------------------
! initialize MPI and create a datatype for real numbers
!
call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,my_block%id,ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)

nwork=9*int(8.0/7.0*float(my_block%nx*my_block%ny*my_block%nz)           &
      /float(my_mg%nxprocs*my_mg%nyprocs*my_mg%nzprocs))        &
      +27*int(4.0/3.0*(                       &
      float(my_block%nx*my_block%ny)/float(my_mg%nxprocs*my_mg%nyprocs)     &
      +float(my_block%nx*my_block%nz)/float(my_mg%nxprocs*my_mg%nzprocs)    &
      +float(my_block%ny*my_block%nz)/float(my_mg%nyprocs*my_mg%nzprocs)))  &
      +81*int(float(my_block%nx)/float(my_mg%nxprocs)        &
      +float(my_block%ny)/float(my_mg%nyprocs)               &
      +float(my_block%nz)/float(my_mg%nzprocs))+243*my_block%ngrid

allocate(work(nwork))

!-----------------------------------------------------------------------
! open file for output of messages and check that the number of 
! processes is correct
!
outfile='out'
outfile(4:5)='_'
m1=mod(my_block%id,10)+1
m0=mod(my_block%id/10,10)+1
outfile(5:6)=num(m0)
outfile(6:7)=num(m1)
open(iout,file=outfile,status='unknown',form='formatted')
if (numprocs.ne.(my_mg%nxprocs*my_mg%nyprocs*my_mg%nzprocs)) then
  write(iout,100)
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
call MPE_DECOMP1D(my_block%nx,dims(1),coords(1),sx,ex)
sx=sx+1
ex=ex+1
call MPE_DECOMP1D(my_block%ny,dims(2),coords(2),sy,ey)
sy=sy+1
ey=ey+1
call MPE_DECOMP1D(my_block%nz,dims(3),coords(3),sz,ez)
sz=sz+1
ez=ez+1

nxdim=int(float(my_block%nx)/float(my_mg%nxprocs)+0.99)+2
nydim=int(float(my_block%ny)/float(my_mg%nyprocs)+0.99)+2
nzdim=int(float(my_block%nz)/float(my_mg%nzprocs)+0.99)+2

if ((ex-sx+3).gt.nxdim) then
  write(iout,110) my_block%id,nxdim,ex-sx+3
  nerror=1
  return
end if
if ((ey-sy+3).gt.nydim) then
  write(iout,120) my_block%id,nydim,ey-sy+3
  nerror=1
  return
end if
if ((ez-sz+3).gt.nzdim) then
  write(iout,130) my_block%id,nzdim,ez-sz+3
  nerror=1
  return
end if
write(iout,*) 'sx=',sx,' ex=',ex,' sy=',sy,' ey=',ey
do i=1,26
  write(iout,*) 'neighbor: ',neighbor(i),' bd: ',bd(i)
end do

!-----------------------------------------------------------------------
! initialize mgd3
!-----------------------------------------------------------------------

my_block%sx = sx; my_block%ex = ex
my_block%sy = sy; my_block%ey = ey
my_block%sz = sz; my_block%ez = ez

allocate(p(sx-1:ex+1,sy-1:ey+1,sz-1:ez+1))
allocate(f(sx-1:ex+1,sy-1:ey+1,sz-1:ez+1))
allocate(r(sx-1:ex+1,sy-1:ey+1,sz-1:ez+1))

call initialize(my_block,nwork,my_mg,iout,nerror)

if (nerror.eq.1) goto 1000
!-----------------------------------------------------------------------
! initialize problem
! xl,yl,zl are the dimensions of the domain
! wk is the wavenumber (must be an integer value)
! rro is the average density
! 1/hxi,1/hyi,1/hzi are the spatial resolutions
!-----------------------------------------------------------------------

xl=1.0d0
yl=1.0d0
zl=1.0d0
wk=5.0d0
rro=1.0d0
hxi=float(my_block%nx)/xl
hyi=float(my_block%ny)/yl
hzi=float(my_block%nz)/zl



call ginit()

!-----------------------------------------------------------------------
! solve using mgd3
!-----------------------------------------------------------------------

call solve(2,sx,ex,sy,ey,sz,ez,p,f,r,my_block%ngrid,work,&
           maxcy,tolmax,kcycle,iprer,ipost,iresw,&
           xl,yl,zl,rro,my_block%nx,my_block%ny,my_block%nz,comm3d,comm3dp,&
           comm3dl,comm3dc,my_block%id,neighbor,bd,my_mg%phibc,&
           iter,.true.,iout,nerror)

if (nerror.eq.1) goto 1000

!-----------------------------------------------------------------------
! compare numerical and exact solutions
!-----------------------------------------------------------------------

call gerr()
!-----------------------------------------------------------------------
close(8)
call MPI_FINALIZE(ierr)
stop
1000  write(iout,200)
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

integer :: i,j,k
real(8) :: cnst,cx,cy,cz,xi,yj,zk

do k=sz-1,ez+1
  do j=sy-1,ey+1
    do i=sx-1,ex+1
      p(i,j,k)=0.0d0
      r(i,j,k)=1.0d0
      f(i,j,k)=0.0d0
    end do
  end do
end do
cnst=-12.0d0*(pi*wk)**2
cx=2.0d0*pi*wk
cy=2.0d0*pi*wk
cz=2.0d0*pi*wk
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
integer :: i,j,k,ierr
real(8) :: errloc,err,cx,cy,cz,exact,zk,yj,xi
!
! calculate local error
!
cx=2.0d0*pi*wk
cy=2.0d0*pi*wk
cz=2.0d0*pi*wk
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
write(6,100) my_block%id,errloc/float(my_block%nx*my_block%ny*my_block%nz),err/float(my_block%nx*my_block%ny*my_block%nz)
100   format(/,i3,3x,'Local error: ',e13.6,'  total error: ',e13.6,/)

return
end


subroutine MPE_DECOMP1D(n,numprocs,myid,s,e)
implicit none 

integer :: n, numprocs, myid, s, e
integer :: nlocal
integer :: deficit
!------------------------------------------------------------------------
!  From the MPE library
!  This file contains a routine for producing a decomposition of a 1-d 
!  array when given a number of processors.  It may be used in "direct" 
!  product decomposition.  The values returned assume a "global" domain 
!  in [1:n]
!
! Code      : tmgd3
! Called in : main
! Calls     : --
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
end

end program test_mgd3
