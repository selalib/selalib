program test_multigrid_2d

use mpi
use diagnostics
# include "compdir.inc"

!-----------------------------------------------------------------------
! Test problem for 2D multigrid parallel code mgd2. mgd2 solves the
! non-separable elliptic equation div(1/r*grad(p))=f on a rectangular
! domain with a staggered uniform grid and either periodic or Neumann
! (set normal derivative of p at boundary to 0) or Dirichlet (set value
! of p at boundary) boundary conditions.
!
! Here, r is constant and equal to 1 everywhere. The right-hand side
! f is chosen so that it corresponds to 
!                p=sin(2*pi*wk*x)*sin(2*pi*wk*y),
! where k is the wavenumber. The program computes the numerical 
! solution and compares it to the exact expression.
!
! The program also shows how to implement mgd2. The user must define
! a few parameters describing the resolution and the number of levels
! to operate on, a communicator and a datatype.
!
! An important feature is that before doing any calculation using
! mgdsolver, the user must call mgdinit to initialize the multigrid
! parameters. Once this is done, mgdsolver can be called any number
! of times, the variables are not overwritten.
!
! Input     : none
! Outputs   : messages -> out* files
! Code      : tmgd2, 2-D parallel multigrid solver
! Calls     : MPI_INIT, MPI_COMM_RANK, MPI_COMM_SIZE, MPI_CART_CREATE,
!             MPI_COMM_RANK, MPI_CART_SHIFT, MPI_BARRIER, MPI_SENDRECV,
!             MPI_CART_GET, MPE_DECOMP1D, 
!             mgdinit, ginit, mgdsolver, gerr
!-----------------------------------------------------------------------
!
! parameters
!
      integer nxp2,nyp2,nx,ny,nxprocs,nyprocs,ixp,jyq,iex,jey
      integer maxcy,kcycle,iprer,ipost,iresw,nwork
      integer ngrid,nxdim,nydim,IOUT
      REALN tolmax
      parameter( nxp2=258, nyp2=258)
      parameter( nx=nxp2-2, ny=nyp2-2)
      parameter( nxprocs=4, nyprocs=2)
      parameter( ixp=4, jyq=4, iex=7, jey=7)
      parameter( maxcy=300, tolmax=1.0d-05)
      parameter( kcycle=1, iprer=2, ipost=1, iresw=1)
      parameter( nwork=(4*nx*ny*8)/(3*nxprocs*nyprocs)  &
     &                +(64*(nx+ny))/3+(32*4)/3)
      parameter( ngrid=max(iex,jey))
      parameter( nxdim=int(float(nxp2-2)/float(nxprocs)+0.99)+2)
      parameter( nydim=int(float(nyp2-2)/float(nyprocs)+0.99)+2)
      parameter( IOUT=6)
!
! variables
!
      logical periods(2)
      integer numprocs,myid,comm2d,sx,ex,sy,ey,neighbor(8),bd(8),iter
      integer realtype,ibdry,jbdry
      integer nerror,m0,m1,dims(2),coords(2),i
      integer status(MPI_STATUS_SIZE),nbrright,nbrbottom,nbrleft,nbrtop
      REALN p(nxdim,nydim),r(nxdim,nydim),f(nxdim,nydim),work(nwork)
      REALN xl,yl,hxi,hyi,vbc(4),phibc(4,20),wk,rro,pi
      character(len=20) outfile
      character(len=1) num(10)
      data (num(i),i=1,10)/'0','1','2','3','4','5','6','7','8','9'/
      pi=4.0d0*atan(1.0d0)

!-----------------------------------------------------------------------
! initialize MPI and create a datatype for real numbers
!
      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)
      realtype=MPI_DOUBLE_PRECISION

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
      !open(IOUT,file=outfile,status='unknown',form='formatted')
      if (numprocs.ne.(nxprocs*nyprocs)) then
        write(IOUT,100)
100     format(/,'ERROR: numprocs <> (nxprocs*nyprocs)',/)
        stop
      end if
!-----------------------------------------------------------------------
! create Cartesian topology with periodic boundary conditions
!
      dims(1)=nxprocs
      dims(2)=nyprocs
      periods(1)=.true.
      periods(2)=.true.
      call MPI_CART_CREATE(MPI_COMM_WORLD,2,dims,periods,.true., &
     &                     comm2d,ierr)
      call MPI_COMM_RANK(comm2d,myid,ierr)
      ibdry=0
      jbdry=0
      vbc(1)=0.0d0
      vbc(2)=0.0d0
      vbc(3)=0.0d0
      vbc(4)=0.0d0
      do i=1,8
        bd(i)=0
      end do
!-----------------------------------------------------------------------
! find the neighbors
! conventions
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
!
      call MPI_CART_SHIFT(comm2d,0,1,nbrleft,nbrright,ierr)
      call MPI_CART_SHIFT(comm2d,1,1,nbrbottom,nbrtop,ierr)
      neighbor(1)=nbrright
      neighbor(3)=nbrbottom
      neighbor(5)=nbrleft
      neighbor(7)=nbrtop
      call MPI_BARRIER(comm2d,ierr)
      call MPI_SENDRECV(neighbor(3),1,MPI_INTEGER,nbrright,0, &
     &                  neighbor(4),1,MPI_INTEGER,nbrleft,0, &
     &                  comm2d,status,ierr)
      call MPI_SENDRECV(neighbor(7),1,MPI_INTEGER,nbrright,1, &
     &                  neighbor(6),1,MPI_INTEGER,nbrleft,1, &
     &                  comm2d,status,ierr)
      call MPI_SENDRECV(neighbor(3),1,MPI_INTEGER,nbrleft,0, &
     &                  neighbor(2),1,MPI_INTEGER,nbrright,0, &
     &                  comm2d,status,ierr)
      call MPI_SENDRECV(neighbor(7),1,MPI_INTEGER,nbrleft,1, &
     &                  neighbor(8),1,MPI_INTEGER,nbrright,1, &
     &                  comm2d,status,ierr)
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
      if ((ex-sx+3).gt.nxdim) then
        write(IOUT,110) myid,nxdim,ex-sx+3
        nerror=1
        goto 1000
      end if
      if ((ey-sy+3).gt.nydim) then
        write(IOUT,120) myid,nydim,ey-sy+3
        nerror=1
        goto 1000
      end if
      write(IOUT,*) 'sx=',sx,' ex=',ex,' sy=',sy,' ey=',ey
      do i=1,8
        write(IOUT,*) 'neighbor: ',neighbor(i),' bd: ',bd(i)
      end do
110   format(/,'ERROR: process:',i3,' nxdim=',i4,' < ex-sx+3=',i4,/, &
     &       ' -> put the parameter formula for nxdim in main.F in ', &
     &       'comments and',/,'    assign to nxdim the maximum ', &
     &       'value of ex-sx+3',/)
120   format(/,'ERROR: process:',i3,' nydim=',i4,' < ey-sy+3=',i4,/, &
     &       ' -> put the parameter formula for nydim in main.F in ', &
     &       'comments and'/,'     assign to nydim the maximum ', &
     &       'value of ey-sy+3',/)
!-----------------------------------------------------------------------
! initialize mgd2
!
      call mgdinit(vbc,phibc,ixp,jyq,iex,jey,ngrid,nxp2,nyp2, &
     &             sx,ex,sy,ey,realtype,nxprocs,nyprocs,nwork, &
     &             ibdry,jbdry,myid,IOUT,nerror)
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
      hxi=float(nxp2-2)/xl
      hyi=float(nyp2-2)/yl
      write(IOUT,*) 'hxi=',hxi,' hyi=',hyi
      call ginit(sx,ex,sy,ey,p,r,f,wk,hxi,hyi,pi)
!-----------------------------------------------------------------------
! solve using mgd2
!
      call mgdsolver(2,sx,ex,sy,ey,p,f,r,ngrid,work, &
     &               maxcy,tolmax,kcycle,iprer,ipost,iresw, &
     &               xl,yl,rro,nx,ny,comm2d,myid,neighbor, &
     &               bd,phibc,iter,.true.,IOUT,nerror)

      if (nerror.eq.1) goto 1000
!-----------------------------------------------------------------------
! compare numerical and exact solutions
!
      call gerr(sx,ex,sy,ey,p,comm2d,wk,hxi,hyi,pi,nx,ny,IOUT)
      call gp_plot2d(myid,numprocs,p,nxdim,nydim,hxi,hyi,sx,sy)
!-----------------------------------------------------------------------
      close(8)
      call mpi_finalize(nerror)
      stop
1000  write(IOUT,200)
200   format(/,'ERROR in multigrid code',/)
      close(8)
      stop
      end

subroutine ginit(sx,ex,sy,ey,p,r,f,wk,hxi,hyi,pi)
# include "compdir.inc"
integer sx,ex,sy,ey
REALN p(sx-1:ex+1,sy-1:ey+1),r(sx-1:ex+1,sy-1:ey+1), &
     &      f(sx-1:ex+1,sy-1:ey+1),hxi,hyi,wk,pi
!-----------------------------------------------------------------------
! Initialize the pressure, density, and right-hand side of the
! elliptic equation div(1/r*grad(p))=f
!
!-----------------------------------------------------------------------
integer i,j
REALN cnst,cx,cy,xi,yj
!
do j=sy-1,ey+1
  do i=sx-1,ex+1
    p(i,j)=0.0d0
    r(i,j)=1.0d0
    f(i,j)=0.0d0
  end do
end do
cnst=-8.0d0*(pi*wk)**2
cx=2.0d0*pi*wk
cy=2.0d0*pi*wk
do j=sy,ey
  yj=(float(j)-1.5d0)/hyi
  do i=sx,ex
    xi=(float(i)-1.5d0)/hxi
    f(i,j)=cnst*sin(cx*xi)*sin(cy*yj)
  end do
end do

ierr = 0
!
return
end

subroutine gerr(sx,ex,sy,ey,p,comm2d,wk,hxi,hyi,pi,nx,ny,IOUT)
# include "compdir.inc"
include "mpif.h"
integer sx,ex,sy,ey,comm2d,IOUT,nx,ny
REALN p(sx-1:ex+1,sy-1:ey+1),wk,hxi,hyi,pi
!-----------------------------------------------------------------------
! Calculate the error between the numerical and exact solution to
! the test problem.
!
! Code      : tmgd2
! Called in : main
! Calls     : MPI_ALLREDUCE
!-----------------------------------------------------------------------
integer i,j
REALN errloc,err,cx,cy,exact, xi, yj
!
! calculate local error
!
cx=2.0d0*pi*wk
cy=2.0d0*pi*wk
errloc=0.0d0
do j=sy,ey
  yj=(float(j)-1.5d0)/hyi
  do i=sx,ex
    xi=(float(i)-1.5d0)/hxi
    exact=sin(cx*xi)*sin(cy*yj)
    errloc=errloc+abs(p(i,j)-exact)
  end do
end do

! calculate global error
!
call MPI_ALLREDUCE(errloc,err,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
     &                   comm2d,ierr)
write(IOUT,100) errloc/float(nx*ny),err/float(nx*ny)
100   format(/,'Local error: ',e13.6,'  total error: ',e13.6,/)
!
ierr = 0

return
end

subroutine gscale(sx,ex,sy,ey,a,avo,acorr,comm2d,nx,ny)
# include "compdir.inc"
include "mpif.h"
integer sx,ex,sy,ey,nx,ny
REALN a(sx-1:ex+1,sy-1:ey+1),avo,acorr
integer comm2d
!------------------------------------------------------------------------
! Rescale the field a so that its average inside the domain
! remains constant and equal to avo. For the density,avo should
! be rro, this ensures conservation of mass. For the pressure,
! avo should be 0 so that the average pressure does not drift
! away from 0, which is the initial value.
!
! Code      : tmgd2
! Called in : mgdsolver
! Calls     : MPI_ALLREDUCE
!------------------------------------------------------------------------
REALN avloc,av
integer i,j
# if cdebug
double precision tinitial
tinitial=MPI_WTIME()
# endif
!
! determine average value
!
avloc=0.0d0
do j=sy,ey
  do i=sx,ex
    avloc=avloc+a(i,j)
  end do
end do

! global reduce across all process
!
call MPI_ALLREDUCE(avloc,av,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
                   comm2d,ierr)
# if cdebug
nallreduce=nallreduce+1
# endif
av=av/float(nx*ny)
!
! do correction
!
acorr=avo-av
do j=sy,ey
  do i=sx,ex
    a(i,j)=a(i,j)+acorr
  end do
end do

# if cdebug
timing(49)=timing(49)+MPI_WTIME()-tinitial
# endif

ierr = 0

return
end 

subroutine MPE_DECOMP1D(n,numprocs,myid,s,e)
# include "compdir.inc"
integer n, numprocs, myid, s, e
integer nlocal
integer deficit
!------------------------------------------------------------------------
!  From the MPE library
!  This file contains a routine for producing a decomposition of a 1-d 
!  array when given a number of processors.  It may be used in "direct" 
!  product decomposition.  The values returned assume a "global" domain 
!  in [1:n]
!
! Code      : tmgd2
! Called in : main
! Calls     : --
!------------------------------------------------------------------------
nlocal  = n / numprocs
s      = myid * nlocal + 1
deficit = mod(n,numprocs)
s      = s + min(myid,deficit)
if (myid .lt. deficit) then
    nlocal = nlocal + 1
endif
e = s + nlocal - 1
if (e .gt. n .or. myid .eq. numprocs-1) e = n

ierr = 0 
return
end
