program tmgd3
# include "compdir.inc"
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
integer, parameter :: nxp2=34, nyp2=34, nzp2=34
integer, parameter :: nx=nxp2-2, ny=nyp2-2, nz=nzp2-2
integer, parameter :: nxprocs=2, nyprocs=2, nzprocs=1
integer, parameter :: ixp=4, jyq=4, kzr=4, iex=4, jey=4, kez=4
integer, parameter :: ngrid=max(iex,jey,kez)
integer, parameter :: maxcy=300
integer, parameter :: kcycle=1, iprer=2, ipost=1, iresw=1
integer, parameter :: nwork=9*int(8.0/7.0*float(nx*ny*nz)     &
                      /float(nxprocs*nyprocs*nzprocs))        &
                      +27*int(4.0/3.0*(                       &
                      float(nx*ny)/float(nxprocs*nyprocs)     &
                      +float(nx*nz)/float(nxprocs*nzprocs)    &
                      +float(ny*nz)/float(nyprocs*nzprocs)))  &
                      +81*int(float(nx)/float(nxprocs)        &
                      +float(ny)/float(nyprocs)               &
                      +float(nz)/float(nzprocs))+243*ngrid
integer, parameter :: nxdim=int(float(nxp2-2)/float(nxprocs)+0.99)+2
integer, parameter :: nydim=int(float(nyp2-2)/float(nyprocs)+0.99)+2
integer, parameter :: nzdim=int(float(nzp2-2)/float(nzprocs)+0.99)+2
integer, parameter :: IOUT=6
real(8), parameter :: tolmax=1.0d-05
!
! variables
!
logical :: nprscr,periods(3)
integer :: numprocs,myid,comm3d,comm3dp,comm3dl,comm3dc
integer :: sx,ex,sy,ey,sz,ez,neighbor(26),bd(26),iter
integer :: realtype,ibdry,jbdry,kbdry,ngb(3)
integer :: ierr,nerror,m0,m1,dims(3),coords(3),i,j,k
integer :: status(MPI_STATUS_SIZE)
REALN p(nxdim,nydim,nzdim),r(nxdim,nydim,nzdim)
REALN f(nxdim,nydim,nzdim),work(nwork)

REALN xl,yl,zl,hxi,hyi,hzi,vbc(6),phibc(6,20),wk,rro,err,pi
character*20 outfile
character*1 num(10)
data (num(i),i=1,10)/'0','1','2','3','4','5','6','7','8','9'/
pi=4.0d0*atan(1.0d0)
!-----------------------------------------------------------------------
! initialize MPI and create a datatype for real numbers
!
call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)
# if double_precision
      realtype=MPI_DOUBLE_PRECISION
# else
      realtype=MPI_REAL
# endif
!-----------------------------------------------------------------------
! open file for output of messages and check that the number of 
! processes is correct
!
outfile='out'
outfile(4:5)='_'
m1=mod(myid,10)+1
m0=mod(myid/10,10)+1
outfile(5:6)=num(m0)
outfile(6:7)=num(m1)
open(IOUT,file=outfile,status='unknown',form='formatted')
if (numprocs.ne.(nxprocs*nyprocs*nzprocs)) then
  write(IOUT,100)
100     format(/,'ERROR: numprocs <> (nxprocs*nyprocs*nzprocs)',/)
  stop
end if
!-----------------------------------------------------------------------
! create Cartesian topology with periodic boundary conditions
!
dims(1)=nxprocs
dims(2)=nyprocs
dims(3)=nzprocs
periods(1)=.true.
periods(2)=.true.
periods(3)=.true.
call MPI_CART_CREATE(MPI_COMM_WORLD,3,dims,periods,.true.,comm3d,ierr)
call MPI_COMM_DUP(comm3d,comm3dp,ierr)
call MPI_COMM_DUP(comm3d,comm3dl,ierr)
call MPI_COMM_DUP(comm3d,comm3dc,ierr)
call MPI_COMM_RANK(comm3d,myid,ierr)
ibdry=0
jbdry=0
kbdry=0
vbc(1)=0.0d0
vbc(2)=0.0d0
vbc(3)=0.0d0
vbc(4)=0.0d0
vbc(5)=0.0d0
vbc(6)=0.0d0
do i=1,26
  bd(i)=0
end do
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
ngb(1)=coords(1)
ngb(2)=coords(2)+1
ngb(3)=coords(3)
call MPI_CART_RANK(comm3d,ngb,neighbor(1),ierr)
ngb(1)=coords(1)
ngb(2)=coords(2)
ngb(3)=coords(3)-1
call MPI_CART_RANK(comm3d,ngb,neighbor(3),ierr)
ngb(1)=coords(1)
ngb(2)=coords(2)-1
ngb(3)=coords(3)
call MPI_CART_RANK(comm3d,ngb,neighbor(5),ierr)
ngb(1)=coords(1)
ngb(2)=coords(2)
ngb(3)=coords(3)+1
call MPI_CART_RANK(comm3d,ngb,neighbor(7),ierr)
ngb(1)=coords(1)-1
ngb(2)=coords(2)
ngb(3)=coords(3)
call MPI_CART_RANK(comm3d,ngb,neighbor(9),ierr)
ngb(1)=coords(1)+1
ngb(2)=coords(2)
ngb(3)=coords(3)
call MPI_CART_RANK(comm3d,ngb,neighbor(18),ierr)
!
! neighbors which share a line with 'myid'
!
ngb(1)=coords(1)
ngb(2)=coords(2)+1
ngb(3)=coords(3)-1
call MPI_CART_RANK(comm3d,ngb,neighbor(2),ierr)
ngb(1)=coords(1)
ngb(2)=coords(2)-1
ngb(3)=coords(3)-1
call MPI_CART_RANK(comm3d,ngb,neighbor(4),ierr)
ngb(1)=coords(1)
ngb(2)=coords(2)-1
ngb(3)=coords(3)+1
call MPI_CART_RANK(comm3d,ngb,neighbor(6),ierr)
ngb(1)=coords(1)
ngb(2)=coords(2)+1
ngb(3)=coords(3)+1
call MPI_CART_RANK(comm3d,ngb,neighbor(8),ierr)
ngb(1)=coords(1)-1
ngb(2)=coords(2)+1
ngb(3)=coords(3)
call MPI_CART_RANK(comm3d,ngb,neighbor(10),ierr)
ngb(1)=coords(1)-1
ngb(2)=coords(2)
ngb(3)=coords(3)-1
call MPI_CART_RANK(comm3d,ngb,neighbor(12),ierr)
ngb(1)=coords(1)-1
ngb(2)=coords(2)-1
ngb(3)=coords(3)
call MPI_CART_RANK(comm3d,ngb,neighbor(14),ierr)
ngb(1)=coords(1)-1
ngb(2)=coords(2)
ngb(3)=coords(3)+1
call MPI_CART_RANK(comm3d,ngb,neighbor(16),ierr)
ngb(1)=coords(1)+1
ngb(2)=coords(2)+1
ngb(3)=coords(3)
call MPI_CART_RANK(comm3d,ngb,neighbor(19),ierr)
ngb(1)=coords(1)+1
ngb(2)=coords(2)
ngb(3)=coords(3)-1
call MPI_CART_RANK(comm3d,ngb,neighbor(21),ierr)
ngb(1)=coords(1)+1
ngb(2)=coords(2)-1
ngb(3)=coords(3)
call MPI_CART_RANK(comm3d,ngb,neighbor(23),ierr)
ngb(1)=coords(1)+1
ngb(2)=coords(2)
ngb(3)=coords(3)+1
call MPI_CART_RANK(comm3d,ngb,neighbor(25),ierr)
!
! neighbors which share a corner with 'myid'
!
ngb(1)=coords(1)-1
ngb(2)=coords(2)+1
ngb(3)=coords(3)-1
call MPI_CART_RANK(comm3d,ngb,neighbor(11),ierr)
ngb(1)=coords(1)-1
ngb(2)=coords(2)-1
ngb(3)=coords(3)-1
call MPI_CART_RANK(comm3d,ngb,neighbor(13),ierr)
ngb(1)=coords(1)-1
ngb(2)=coords(2)-1
ngb(3)=coords(3)+1
call MPI_CART_RANK(comm3d,ngb,neighbor(15),ierr)
ngb(1)=coords(1)-1
ngb(2)=coords(2)+1
ngb(3)=coords(3)+1
call MPI_CART_RANK(comm3d,ngb,neighbor(17),ierr)
ngb(1)=coords(1)+1
ngb(2)=coords(2)+1
ngb(3)=coords(3)-1
call MPI_CART_RANK(comm3d,ngb,neighbor(20),ierr)
ngb(1)=coords(1)+1
ngb(2)=coords(2)-1
ngb(3)=coords(3)-1
call MPI_CART_RANK(comm3d,ngb,neighbor(22),ierr)
ngb(1)=coords(1)+1
ngb(2)=coords(2)-1
ngb(3)=coords(3)+1
call MPI_CART_RANK(comm3d,ngb,neighbor(24),ierr)
ngb(1)=coords(1)+1
ngb(2)=coords(2)+1
ngb(3)=coords(3)+1
call MPI_CART_RANK(comm3d,ngb,neighbor(26),ierr)
!-----------------------------------------------------------------------
! find indices of subdomain and check that dimensions of arrays are
! sufficient
!
call MPE_DECOMP1D(nxp2-2,dims(1),coords(1),sx,ex)
sx=sx+1
ex=ex+1
call MPE_DECOMP1D(nyp2-2,dims(2),coords(2),sy,ey)
sy=sy+1
ey=ey+1
call MPE_DECOMP1D(nzp2-2,dims(3),coords(3),sz,ez)
sz=sz+1
ez=ez+1
if ((ex-sx+3).gt.nxdim) then
  write(IOUT,110) myid,nxdim,ex-sx+3
  nerror=1
  return
end if
if ((ey-sy+3).gt.nydim) then
  write(IOUT,120) myid,nydim,ey-sy+3
  nerror=1
  return
end if
if ((ez-sz+3).gt.nzdim) then
  write(IOUT,130) myid,nzdim,ez-sz+3
  nerror=1
  return
end if
write(IOUT,*) 'sx=',sx,' ex=',ex,' sy=',sy,' ey=',ey
do i=1,26
  write(IOUT,*) 'neighbor: ',neighbor(i),' bd: ',bd(i)
end do
110   format(/,'ERROR: process:',i3,' nxdim=',i4,' < ex-sx+3=',i4,/, &
       ' -> put the parameter formula for nxdim in main.F in ', &
       'comments and',/,'    assign to nxdim the maximum ', &
       'value of ex-sx+3',/)
120   format(/,'ERROR: process:',i3,' nydim=',i4,' < ey-sy+3=',i4,/, &
            ' -> put the parameter formula for nydim in main.F in ', &
            'comments and'/,'     assign to nydim the maximum ', &
            'value of ey-sy+3',/)
130   format(/,'ERROR: process:',i3,' nzdim=',i4,' < ez-sz+3=',i4,/, &
            ' -> put the parameter formula for nzdim in main.F in ', &
            'comments and'/,'     assign to nzdim the maximum ', &
            'value of ez-sz+3',/)
!-----------------------------------------------------------------------
! initialize mgd3
!
call mgdinit(vbc,phibc,ixp,jyq,kzr,iex,jey,kez,ngrid,nxp2,	&
             nyp2,nzp2,sx,ex,sy,ey,sz,ez,realtype,nxprocs,	&
             nyprocs,nzprocs,nwork,ibdry,jbdry,kbdry,myid,	&
             IOUT,nerror)
if (nerror.eq.1) goto 1000
!-----------------------------------------------------------------------
! initialize problem
! xl,yl,zl are the dimensions of the domain
! wk is the wavenumber (must be an integer value)
! rro is the average density
! 1/hxi,1/hyi,1/hzi are the spatial resolutions
!
xl=1.0d0
yl=1.0d0
zl=1.0d0
wk=5.0d0
rro=1.0d0
hxi=float(nxp2-2)/xl
hyi=float(nyp2-2)/yl
hzi=float(nzp2-2)/zl
write(IOUT,*) 'hxi=',hxi,' hyi=',hyi,' hzi=',hzi
call ginit(sx,ex,sy,ey,sz,ez,p,r,f,wk,hxi,hyi,hzi,pi,IOUT)
!-----------------------------------------------------------------------
! solve using mgd3
!
call mgdsolver(2,sx,ex,sy,ey,sz,ez,p,f,r,ngrid,work,	&
               maxcy,tolmax,kcycle,iprer,ipost,iresw,	&
               xl,yl,zl,rro,nx,ny,nz,comm3d,comm3dp,	&
               comm3dl,comm3dc,myid,neighbor,bd,phibc,	&
               iter,.true.,IOUT,nerror)
if (nerror.eq.1) goto 1000
!-----------------------------------------------------------------------
! compare numerical and exact solutions
!
call gerr(sx,ex,sy,ey,sz,ez,p,comm3d,wk,hxi,hyi,hzi,pi,nx,ny,nz,IOUT)
!-----------------------------------------------------------------------
close(8)
call MPI_FINALIZE(ierr)
stop
1000  write(IOUT,200)
200   format(/,'ERROR in multigrid code',/)
close(8)
stop
end program tmgd3
