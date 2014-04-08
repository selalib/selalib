module sll_multigrid_2d
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_boundary_condition_descriptors.h"

   use sll_collective

   implicit none

   type, public :: multigrid_2d

      sll_int32  :: sx
      sll_int32  :: ex
      sll_int32  :: sy
      sll_int32  :: ey
      sll_int32  :: comm2d
      sll_int32  :: pdims(2)
      sll_int32  :: coords(2)

      sll_real64, pointer :: work(:)

   end type multigrid_2d

   interface initialize
      module procedure initialize_multigrid_2d
   end interface initialize

   interface solve
      module procedure solve_multigrid_2d
   end interface solve

   public :: initialize, solve, write_topology

   private

   sll_int32  :: bd(8)
   sll_int32  :: nx
   sll_int32  :: ny
   sll_real64 :: phibc(4,20)
   sll_real64 :: wk
   sll_int32  :: myid
   sll_int32  :: neighbor(8)
   sll_int32  :: ngrid
   sll_int32  :: nerror
   sll_real64 :: hxi
   sll_real64 :: hyi
   sll_real64 :: xl
   sll_real64 :: yl

   sll_int32,  parameter :: IOUT    = 6

   sll_int32,  parameter :: maxcy   = 2000     !Max number of mg cycles
   sll_real64, parameter :: tolmax  = 1.0d-06 !Desired accuracy
   sll_int32,  parameter :: kcycle  = 2       !1: V-cycles 2: W-cycles
   sll_int32,  parameter :: iprer   = 3       !number of relaxation sweep
   sll_int32,  parameter :: ipost   = 1
   sll_int32,  parameter :: iresw   = 1       !1: fully weighted residual
                                              !2: half-weighted residual

contains

subroutine initialize_multigrid_2d( this,            &
                       x_min, x_max, y_min, y_max,   &
                       nxprocs, nyprocs, nc_x, nc_y, &
                       bctype_x, bctype_y,           &
                       bcvalue_xmin, bcvalue_xmax,   &
                       bcvalue_ymin, bcvalue_ymax  )

   type(multigrid_2d)     :: this

   sll_real64, intent(in) :: x_min
   sll_real64, intent(in) :: x_max
   sll_real64, intent(in) :: y_min
   sll_real64, intent(in) :: y_max

   sll_int32, intent(in)  :: nxprocs
   sll_int32, intent(in)  :: nyprocs

   sll_int32, intent(in)  :: nc_x
   sll_int32, intent(in)  :: nc_y

   sll_int32              :: bctype_x
   sll_int32              :: bctype_y

   sll_real64, optional   :: bcvalue_xmin
   sll_real64, optional   :: bcvalue_xmax
   sll_real64, optional   :: bcvalue_ymin
   sll_real64, optional   :: bcvalue_ymax

   sll_int32 :: nxdim
   sll_int32 :: nydim
   sll_int32 :: ixp
   sll_int32 :: jyq 

   logical    :: periods(2)

   sll_int32  :: iex, jey

   sll_int32  :: numprocs
   sll_int32  :: ibdry
   sll_int32  :: jbdry
   sll_int32  :: i
   sll_int32  :: status(MPI_STATUS_SIZE)
   sll_int32  :: nbrright
   sll_int32  :: nbrbottom
   sll_int32  :: nbrleft
   sll_int32  :: nbrtop
   sll_int32  :: nwork
   sll_int32  :: ierr

   sll_int32  :: sx
   sll_int32  :: ex
   sll_int32  :: sy
   sll_int32  :: ey
   sll_int32  :: comm2d
   sll_int32  :: coords(2)

   sll_real64, dimension(4) :: vbc

   nx = nc_x
   ny = nc_y

   ixp = nx !2*nxprocs 
   jyq = ny !2*nyprocs

   iex   = ceiling(log((nx-1.)/ixp)/log(2.))+1
   jey   = ceiling(log((ny-1.)/jyq)/log(2.))+1
   ngrid = max(iex,jey)

   write(*,*) " ngrid = ", ngrid

   nwork = (4*nx*ny*8)/(3*nxprocs*nyprocs)+(64*(nx+ny))/3+(32*4)/3
   nxdim = int(float(nx)/float(nxprocs)+0.99)+2
   nydim = int(float(ny)/float(nyprocs)+0.99)+2

   xl = x_max - x_min
   yl = y_max - y_min
   
   hxi = real(nx,f64)
   hyi = real(ny,f64)

!-----------------------------------------------------------------------
! initialize MPI and create a datatype for real numbers
!

   call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)

   if (numprocs.ne.(nxprocs*nyprocs)) then
      write(IOUT,"(/,'ERROR: numprocs <> (nxprocs*nyprocs)',/)")
      call MPI_FINALIZE(ierr)
      stop
   end if

   vbc(1:4) = 0.0d0

   if (bctype_x == SLL_PERIODIC) then

      ibdry       = 0
      periods(1)  = .true.

   else if (bctype_x == SLL_NEUMANN) then

      ibdry       = 1
      periods(1)  = .false.

   else if (bctype_x == SLL_DIRICHLET) then

      ibdry       = 2
      periods(1)  = .false.
      if (present(bcvalue_xmax) .and. present(bcvalue_xmin)) then
         vbc(1)   = bcvalue_xmax
         vbc(3)   = bcvalue_xmin
      end if

   end if 

   if (bctype_y == SLL_PERIODIC) then

      jbdry       = 0
      periods(2)  = .true.

   else if (bctype_y == SLL_NEUMANN) then

      jbdry       = 1
      periods(2)  = .false.

   else if (bctype_y == SLL_DIRICHLET) then

      jbdry       = 2
      periods(2)  = .false.
      if (present(bcvalue_ymax) .and. present(bcvalue_ymin)) then
         vbc(2)   = bcvalue_ymin
         vbc(4)   = bcvalue_ymax
      end if

   end if 

   bd(1:8)  = 0

   ! create Cartesian topology 
   neighbor = MPI_PROC_NULL
   this%pdims(1) = nxprocs
   this%pdims(2) = nyprocs
   call MPI_CART_CREATE(MPI_COMM_WORLD,2,this%pdims,periods,.true.,comm2d,ierr)

   !  For the Dirichlet boundaries, it is possible to assign
   !  4 different constant values to the 4 different sides. The values are
   !  assigned at the finest grid level, zero is assigned at all levels below
   ! 
   !  vbc:
   !
   !       -----vbc(4)------ 
   !       |               |
   !       |               |
   !     vbc(3)----------vbc(1)
   !       |               |
   !       |               |
   !       ------vbc(2)-----
   !
   !

   !-----------------------------------------------------------------------
   ! neighbors conventions
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
   !-----------------------------------------------------------------------

   call MPI_CART_SHIFT(comm2d,0,1,nbrleft,nbrright,ierr)
   call MPI_CART_SHIFT(comm2d,1,1,nbrbottom,nbrtop,ierr)
   neighbor(1)=nbrright
   neighbor(3)=nbrbottom
   neighbor(5)=nbrleft
   neighbor(7)=nbrtop

   write(*,*) myid, neighbor(1), neighbor(3), neighbor(5), neighbor(7)
   call flush(6)

   CALL MPI_CART_COORDS(comm2d,myid,2,coords,ierr)
   this%coords = coords

   if (nbrright /= MPI_PROC_NULL .and. nbrbottom /= MPI_PROC_NULL) &
      CALL MPI_CART_RANK(comm2d,(/coords(1)+1,coords(2)-1/),neighbor(2),ierr)
   if (nbrleft /= MPI_PROC_NULL .and. nbrbottom /= MPI_PROC_NULL) &
      CALL MPI_CART_RANK(comm2d,(/coords(1)-1,coords(2)-1/),neighbor(4),ierr)
 
   if (nbrleft /= MPI_PROC_NULL .and. nbrtop /= MPI_PROC_NULL) &
      CALL MPI_CART_RANK(comm2d,(/coords(1)-1,coords(2)+1/),neighbor(6),ierr)
   if (nbrright /= MPI_PROC_NULL .and. nbrtop /= MPI_PROC_NULL) &
      CALL MPI_CART_RANK(comm2d,(/coords(1)+1,coords(2)+1/),neighbor(8),ierr)

   !-----------------------------------------------------------------------
   ! find indices of subdomain and check that dimensions of arrays are
   ! sufficient
   !

   call MPI_CART_GET(comm2d,2,this%pdims,periods,coords,ierr)

   call MPE_DECOMP1D(nx,this%pdims(1),coords(1),sx,ex)
   sx=sx+1
   ex=ex+1
   call MPE_DECOMP1D(ny,this%pdims(2),coords(2),sy,ey)
   sy=sy+1
   ey=ey+1

   if ((ex-sx+3) > nxdim) then
      write(IOUT,110) myid,nxdim,ex-sx+3
      nerror=1
      write(IOUT,"(/,'ERROR in multigrid code : 178',/)")
   end if

   if ((ey-sy+3) > nydim) then
      write(IOUT,120) myid,nydim,ey-sy+3
      nerror=1
      write(IOUT,"(/,'ERROR in multigrid code : 183',/)")
   end if

   write(IOUT,*) 'sx=',sx,' ex=',ex,' sy=',sy,' ey=',ey

   where(neighbor >= 0)
      bd = 0
   elsewhere
      bd = 1
   end where

   do i=1,8
      write(IOUT,*) myid, 'neighbor: ',neighbor(i),' bd: ',bd(i)
   end do

   call flush(6)

   call mgdinit(   vbc,         &
                   phibc,       &
                   ixp,         &
                   jyq,         &
                   iex,         &
                   jey,         &
                   ngrid,       &
                   nx+2,        &
                   ny+2,        &
                   sx,          &
                   ex,          &
                   sy,          &
                   ey,          &
                   MPI_REAL8,   &
                   nxprocs,     &
                   nyprocs,     &
                   nwork,       &
                   ibdry,       &
                   jbdry,       &
                   myid,        &
                   IOUT,        &
                   nerror)

   SLL_CLEAR_ALLOCATE(this%work(1:nwork),ierr)

   this%ex = ex
   this%sx = sx
   this%ey = ey
   this%sy = sy
   this%comm2d = comm2d

110 format(/,'ERROR: process:',i3,' nxdim=',i4,' < ex-sx+3=',i4,/, &
    &       ' -> put the parameter formula for nxdim in main.F in ', &
    &       'comments and',/,'    assign to nxdim the maximum ', &
    &       'value of ex-sx+3',/)
120 format(/,'ERROR: process:',i3,' nydim=',i4,' < ey-sy+3=',i4,/, &
    &       ' -> put the parameter formula for nydim in main.F in ', &
    &       'comments and'/,'     assign to nydim the maximum ', &
    &       'value of ey-sy+3',/)
end subroutine initialize_multigrid_2d

subroutine solve_multigrid_2d(this, p, f, r)

   type(multigrid_2d)         :: this
   sll_real64, dimension(:,:) :: p
   sll_real64, dimension(:,:) :: f
   sll_real64, dimension(:,:) :: r

   sll_int32  :: iter

   sll_real64, parameter :: rro = 0.0_f64 !Average density

!-----------------------------------------------------------------------
! solve using mgd2
!
   call mgdsolver(1,this%sx,this%ex,this%sy,this%ey,p,f,r,ngrid,this%work, &
  &               maxcy,tolmax,kcycle,iprer,ipost,iresw, &
  &               xl,yl,rro,nx,ny,this%comm2d,myid,neighbor, &
  &               bd,phibc,iter,.true.,IOUT,nerror)


end subroutine solve_multigrid_2d

subroutine MPE_DECOMP1D(n,numprocs,myid,s,e)

   sll_int32 :: n, numprocs, myid, s, e
   sll_int32 :: nlocal
   sll_int32 :: deficit

   !------------------------------------------------------------------------
   !  From the MPE library
   !  This file contains a routine for producing a decomposition of a 1-d 
   !  array when given a number of processors.  It may be used in "direct" 
   !  product decomposition.  The values returned assume a "global" domain 
   !  in [1:n]
   !------------------------------------------------------------------------

   nlocal  = n / numprocs
   s      = myid * nlocal + 1
   deficit = mod(n,numprocs)
   s      = s + min(myid,deficit)
   if (myid  < deficit) then
       nlocal = nlocal + 1
   endif
   e = s + nlocal - 1
   if (e  >  n .or. myid == numprocs-1) e = n

end subroutine mpe_decomp1d

subroutine write_topology( this )

   type(multigrid_2d) :: this

   sll_real64 :: xp, dx
   sll_real64 :: yp, dy
   sll_int32  :: prank
   sll_int32  :: psize
   sll_int32  :: code
   sll_int32  :: iproc, jproc
   sll_int32  :: tag = 1111
   sll_int32  :: statut(MPI_STATUS_SIZE)

   call MPI_COMM_RANK(MPI_COMM_WORLD,prank,code)
   call MPI_COMM_SIZE(MPI_COMM_WORLD,psize,code)

   xp = real(this%coords(1),4)/this%pdims(1) * xl 
   yp = real(this%coords(2),4)/this%pdims(2) * yl 

   dx = 1. / this%pdims(1) * xl
   dy = 1. / this%pdims(2) * yl
   
   if ( prank == 0) then
   
      open(4,file="mpi.mtv")
      write(4,*)"$DATA=CONTOUR"
      write(4,*)"%equalscale=T"
      write(4,*)"%interp     = 0"
      write(4,*)"%contstyle  = 0"
      write(4,*)"%meshplot   = on" 
      write(4,*)"%hiddenline = off" 
      write(4,*)"%xmin=",0., " xmax = ", xl
      write(4,*)"%ymin=",0., " ymax = ", yl
      write(4,*)"% nx   = ", this%pdims(1)+1
      write(4,*)"% ny   = ", this%pdims(2)+1
      do jproc = 1, this%pdims(2)+1
         write(4,*)(0., iproc=1,this%pdims(1)+1)
      end do

      do iproc = 0, psize-1
         if (iproc > 0) then
            call MPI_RECV(xp,1,MPI_REAL8,iproc,tag,MPI_COMM_WORLD,statut,code)
            call MPI_RECV(yp,1,MPI_REAL8,iproc,tag,MPI_COMM_WORLD,statut,code)
         end if
         write(4,111)xp+.5*dx,yp+.5*dy,iproc
      end do
   
      write(4,"('$END')")
      close(4)

   else

      call MPI_SEND(xp, 1, MPI_REAL8, 0, tag, MPI_COMM_WORLD, code)
      call MPI_SEND(yp, 1, MPI_REAL8, 0, tag, MPI_COMM_WORLD, code)

   end if

111 format("@text x1=",f7.3,' y1=',f7.3," z1=0.0 linelabel='",i4,"'")

end subroutine write_topology

end module sll_multigrid_2d
