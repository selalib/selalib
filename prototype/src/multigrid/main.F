      program main
      use mpi
      use diagnostics
# include "compdir.inc"
c-----------------------------------------------------------------------
c Test problem for 2D multigrid parallel code mgd2. mgd2 solves the
c non-separable elliptic equation div(1/r*grad(p))=f on a rectangular
c domain with a staggered uniform grid and either periodic or Neumann
c (set normal derivative of p at boundary to 0) or Dirichlet (set value
c of p at boundary) boundary conditions.
c
c Here, r is constant and equal to 1 everywhere. The right-hand side
c f is chosen so that it corresponds to 
c                p=sin(2*pi*wk*x)*sin(2*pi*wk*y),
c where k is the wavenumber. The program computes the numerical 
c solution and compares it to the exact expression.
c
c The program also shows how to implement mgd2. The user must define
c a few parameters describing the resolution and the number of levels
c to operate on, a communicator and a datatype.
c
c An important feature is that before doing any calculation using
c mgdsolver, the user must call mgdinit to initialize the multigrid
c parameters. Once this is done, mgdsolver can be called any number
c of times, the variables are not overwritten.
c
c Input     : none
c Outputs   : messages -> out* files
c Code      : tmgd2, 2-D parallel multigrid solver
c Calls     : MPI_INIT, MPI_COMM_RANK, MPI_COMM_SIZE, MPI_CART_CREATE,
c             MPI_COMM_RANK, MPI_CART_SHIFT, MPI_BARRIER, MPI_SENDRECV,
c             MPI_CART_GET, MPE_DECOMP1D, 
c             mgdinit, ginit, mgdsolver, gerr
c-----------------------------------------------------------------------
c
c parameters
c
      integer nxp2,nyp2,nx,ny,nxprocs,nyprocs,ixp,jyq,iex,jey
      integer maxcy,kcycle,iprer,ipost,iresw,nwork
      integer ngrid,nxdim,nydim,IOUT
      REALN tolmax
      parameter( nxp2=258, nyp2=258)
      parameter( nx=nxp2-2, ny=nyp2-2)
      parameter( nxprocs=2, nyprocs=2)
      parameter( ixp=4, jyq=4, iex=7, jey=7)
      parameter( maxcy=300, tolmax=1.0d-05)
      parameter( kcycle=1, iprer=2, ipost=1, iresw=1)
      parameter( nwork=(4*nx*ny*8)/(3*nxprocs*nyprocs)
     1                +(64*(nx+ny))/3+(32*4)/3)
      parameter( ngrid=max(iex,jey))
      parameter( nxdim=int(float(nxp2-2)/float(nxprocs)+0.99)+2)
      parameter( nydim=int(float(nyp2-2)/float(nyprocs)+0.99)+2)
      parameter( IOUT=6)
c
c variables
c
      logical nprscr,periods(2)
      integer numprocs,myid,comm2d,sx,ex,sy,ey,neighbor(8),bd(8),iter
      integer realtype,ibdry,jbdry
      integer ierr,nerror,m0,m1,dims(2),coords(2),i,j
      integer status(MPI_STATUS_SIZE),nbrright,nbrbottom,nbrleft,nbrtop
      REALN p(nxdim,nydim),r(nxdim,nydim),f(nxdim,nydim),work(nwork)
      REALN xl,yl,hxi,hyi,vbc(4),phibc(4,20),wk,rro,err,pi
      character(len=20) outfile
      character(len=1) num(10)
      data (num(i),i=1,10)/'0','1','2','3','4','5','6','7','8','9'/
      pi=4.0d0*atan(1.0d0)
c-----------------------------------------------------------------------
c initialize MPI and create a datatype for real numbers
c
      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)
# if double_precision
      realtype=MPI_DOUBLE_PRECISION
# else
      realtype=MPI_REAL
# endif
c-----------------------------------------------------------------------
c open file for output of messages and check that the number of 
c processes is correct
c
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
c-----------------------------------------------------------------------
c create Cartesian topology with periodic boundary conditions
c
      dims(1)=nxprocs
      dims(2)=nyprocs
      periods(1)=.true.
      periods(2)=.true.
      call MPI_CART_CREATE(MPI_COMM_WORLD,2,dims,periods,.true.,
     1                     comm2d,ierr)
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
c-----------------------------------------------------------------------
c find the neighbors
c conventions
c
c     6 |     7     | 8
c       |           |
c     -----------------
c       |           |
c       |           |
c     5 |   myid    | 1 
c       |           |
c       |           |
c     -----------------
c       |           |
c     4 |     3     | 2
c
      call MPI_CART_SHIFT(comm2d,0,1,nbrleft,nbrright,ierr)
      call MPI_CART_SHIFT(comm2d,1,1,nbrbottom,nbrtop,ierr)
      neighbor(1)=nbrright
      neighbor(3)=nbrbottom
      neighbor(5)=nbrleft
      neighbor(7)=nbrtop
      call MPI_BARRIER(comm2d,ierr)
      call MPI_SENDRECV(neighbor(3),1,MPI_INTEGER,nbrright,0,
     1                  neighbor(4),1,MPI_INTEGER,nbrleft,0,
     1                  comm2d,status,ierr)
      call MPI_SENDRECV(neighbor(7),1,MPI_INTEGER,nbrright,1,
     1                  neighbor(6),1,MPI_INTEGER,nbrleft,1,
     1                  comm2d,status,ierr)
      call MPI_SENDRECV(neighbor(3),1,MPI_INTEGER,nbrleft,0,
     1                  neighbor(2),1,MPI_INTEGER,nbrright,0,
     1                  comm2d,status,ierr)
      call MPI_SENDRECV(neighbor(7),1,MPI_INTEGER,nbrleft,1,
     1                  neighbor(8),1,MPI_INTEGER,nbrright,1,
     1                  comm2d,status,ierr)
c-----------------------------------------------------------------------
c find indices of subdomain and check that dimensions of arrays are
c sufficient
c
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
110   format(/,'ERROR: process:',i3,' nxdim=',i4,' < ex-sx+3=',i4,/,
     1       ' -> put the parameter formula for nxdim in main.F in ',
     2       'comments and',/,'    assign to nxdim the maximum ',
     3       'value of ex-sx+3',/)
120   format(/,'ERROR: process:',i3,' nydim=',i4,' < ey-sy+3=',i4,/,
     1       ' -> put the parameter formula for nydim in main.F in ',
     2       'comments and'/,'     assign to nydim the maximum ',
     3       'value of ey-sy+3',/)
c-----------------------------------------------------------------------
c initialize mgd2
c
      call mgdinit(vbc,phibc,ixp,jyq,iex,jey,ngrid,nxp2,nyp2,
     1             sx,ex,sy,ey,realtype,nxprocs,nyprocs,nwork,
     2             ibdry,jbdry,myid,IOUT,nerror)
      if (nerror.eq.1) goto 1000
c-----------------------------------------------------------------------
c initialize problem
c xl,yl are the dimensions of the domain
c wk is the wavenumber (must be an integer value)
c rro is the average density
c 1/hxi,1/hyi are the spatial resolutions
c
      xl=1.0d0
      yl=1.0d0
      wk=5.0d0
      rro=1.0d0
      hxi=float(nxp2-2)/xl
      hyi=float(nyp2-2)/yl
      write(IOUT,*) 'hxi=',hxi,' hyi=',hyi
      call ginit(sx,ex,sy,ey,p,r,f,wk,hxi,hyi,pi,IOUT)
c-----------------------------------------------------------------------
c solve using mgd2
c
      call mgdsolver(2,sx,ex,sy,ey,p,f,r,ngrid,work,
     1               maxcy,tolmax,kcycle,iprer,ipost,iresw,
     2               xl,yl,rro,nx,ny,comm2d,myid,neighbor,
     3               bd,phibc,iter,.true.,IOUT,nerror)
      if (nerror.eq.1) goto 1000
c-----------------------------------------------------------------------
c compare numerical and exact solutions
c
      call gerr(sx,ex,sy,ey,p,comm2d,wk,hxi,hyi,pi,nx,ny,IOUT)
      call gp_plot2d(myid,numprocs,p,nxdim,nydim,hxi,hyi,sx,sy)
c-----------------------------------------------------------------------
      close(8)
      call mpi_finalize(nerror)
      stop
1000  write(IOUT,200)
200   format(/,'ERROR in multigrid code',/)
      close(8)
      stop
      end
