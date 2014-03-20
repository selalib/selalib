program test_multigrid_2d

use mpi
use diagnostics
use sll_multigrid_2d
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
      type(mgd2_solver) :: mgd

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


      pi=4.0d0*atan(1.0d0)

!-----------------------------------------------------------------------
! initialize MPI and create a datatype for real numbers
!
      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)
      realtype=MPI_DOUBLE_PRECISION

      if (numprocs.ne.(nxprocs*nyprocs)) then
        write(IOUT,100)
100     format(/,'ERROR: numprocs <> (nxprocs*nyprocs)',/)
        call mpi_finalize(ierr)
        stop
      end if


      call initialize(mgd, &
                      0._f64, 1._f64, nx, nxprocs, SLL_PERIODIC, &
                      0._f64, 1._f64, ny, nyprocs, SLL_PERIODIC)
     

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
integer i,j, ierr
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

return
end

subroutine gscale(sx,ex,sy,ey,a,avo,acorr,comm2d,nx,ny)
# include "compdir.inc"
include "mpif.h"
integer sx,ex,sy,ey,nx,ny
REALN a(sx-1:ex+1,sy-1:ey+1),avo,acorr
integer comm2d, ierr
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

return
end
