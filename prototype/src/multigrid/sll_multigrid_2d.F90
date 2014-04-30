module sll_multigrid_2d
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_boundary_condition_descriptors.h"

   use sll_collective
   use sll_ascii_io

   implicit none

   type, public :: multigrid_2d

      sll_int32  :: sx
      sll_int32  :: ex
      sll_int32  :: sy
      sll_int32  :: ey
      sll_real64 :: xl
      sll_real64 :: yl
      sll_int32  :: comm2d
      sll_int32  :: pdims(2)
      sll_int32  :: coords(2)
      sll_int32  :: neighbor(8)

      sll_real64, pointer :: work(:)

   end type multigrid_2d

   interface initialize
      module procedure initialize_multigrid_2d
      module procedure initialize_multigrid_2d_periodic
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
   sll_int32  :: ngrid
   sll_int32  :: nerror

   sll_int32,  parameter :: IOUT    = 6
   sll_int32,  parameter :: maxcy   = 5000    !Max number of mg cycles
   sll_real64, parameter :: tolmax  = 1.0d-07 !Desired accuracy
   sll_int32,  parameter :: kcycle  = 2       !1: V-cycles 2: W-cycles
   sll_int32,  parameter :: iprer   = 5       !number of relaxation sweep
   sll_int32,  parameter :: ipost   = 1
   sll_int32,  parameter :: iresw   = 1       !1: fully weighted residual
                                              !2: half-weighted residual

contains

subroutine initialize_multigrid_2d_periodic( this,     &
                                             sx,       &
                                             ex,       &
                                             sy,       &
                                             ey,       &
                                             xl,       &
                                             yl,       &
                                             nxprocs,  &
                                             nyprocs,  &
                                             nc_x,     &
                                             nc_y,     &
                                             comm2d,   &
                                             neighbor  )

   type(multigrid_2d)     :: this

   sll_int32, intent(in)  :: sx
   sll_int32, intent(in)  :: ex
   sll_int32, intent(in)  :: sy
   sll_int32, intent(in)  :: ey

   sll_real64, intent(in) :: xl
   sll_real64, intent(in) :: yl

   sll_int32, intent(in)  :: nxprocs
   sll_int32, intent(in)  :: nyprocs

   sll_int32, intent(in)  :: nc_x
   sll_int32, intent(in)  :: nc_y
   sll_int32, intent(in)  :: comm2d
   sll_int32, intent(in)  :: neighbor(8)

   sll_int32 :: nxdim
   sll_int32 :: nydim
   sll_int32 :: ixp
   sll_int32 :: jyq 

   sll_int32  :: iex, jey

   sll_int32  :: ibdry
   sll_int32  :: jbdry
   sll_int32  :: nwork
   sll_int32  :: ierr

   sll_real64, dimension(4) :: vbc

   this%xl = xl
   this%yl = yl

   nx = nc_x
   ny = nc_y

   ixp = 4*nxprocs 
   jyq = 4*nyprocs

   iex   = ceiling(log((nx-1.)/ixp)/log(2.))+1
   jey   = ceiling(log((ny-1.)/jyq)/log(2.))+1
   ngrid = max(iex,jey)

   nwork = (4*nx*ny*8)/(3*nxprocs*nyprocs)+(64*(nx+ny))/3+(32*4)/3
   nxdim = int(float(nx)/float(nxprocs)+0.99)+2
   nydim = int(float(ny)/float(nyprocs)+0.99)+2

   call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)

   vbc(1:4) = 0.0d0

   ibdry    = 0
   jbdry    = 0

   bd(1:8)  = 0

   this%neighbor = neighbor

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

   this%ex     = ex
   this%sx     = sx
   this%ey     = ey
   this%sy     = sy
   this%comm2d = comm2d

end subroutine initialize_multigrid_2d_periodic

subroutine initialize_multigrid_2d( this,                         &
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

   ixp = 4*nxprocs 
   jyq = 4*nyprocs

   iex   = ceiling(log((nx-1.)/ixp)/log(2.))+1
   jey   = ceiling(log((ny-1.)/jyq)/log(2.))+1
   ngrid = max(iex,jey)

   nwork = (4*nx*ny*8)/(3*nxprocs*nyprocs)+(64*(nx+ny))/3+(32*4)/3
   nxdim = int(float(nx)/float(nxprocs)+0.99)+2
   nydim = int(float(ny)/float(nyprocs)+0.99)+2

   this%xl = x_max - x_min
   this%yl = y_max - y_min
   
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
   this%neighbor = MPI_PROC_NULL
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
   this%neighbor(1)=nbrright
   this%neighbor(3)=nbrbottom
   this%neighbor(5)=nbrleft
   this%neighbor(7)=nbrtop

   CALL MPI_CART_COORDS(comm2d,myid,2,coords,ierr)
   this%coords = coords

   if (nbrright /= MPI_PROC_NULL .and. nbrbottom /= MPI_PROC_NULL) &
      CALL MPI_CART_RANK(comm2d,(/coords(1)+1,coords(2)-1/),this%neighbor(2),ierr)
   if (nbrleft /= MPI_PROC_NULL .and. nbrbottom /= MPI_PROC_NULL) &
      CALL MPI_CART_RANK(comm2d,(/coords(1)-1,coords(2)-1/),this%neighbor(4),ierr)
 
   if (nbrleft /= MPI_PROC_NULL .and. nbrtop /= MPI_PROC_NULL) &
      CALL MPI_CART_RANK(comm2d,(/coords(1)-1,coords(2)+1/),this%neighbor(6),ierr)
   if (nbrright /= MPI_PROC_NULL .and. nbrtop /= MPI_PROC_NULL) &
      CALL MPI_CART_RANK(comm2d,(/coords(1)+1,coords(2)+1/),this%neighbor(8),ierr)

   !-----------------------------------------------------------------------
   ! find indices of subdomain and check that dimensions of arrays are
   ! sufficient
   !

   call MPI_CART_GET(comm2d,2,this%pdims,periods,coords,ierr)

   call mpe_decomp1d(nx,this%pdims(1),coords(1),sx,ex)
   sx=sx+1
   ex=ex+1
   call mpe_decomp1d(ny,this%pdims(2),coords(2),sy,ey)
   sy=sy+1
   ey=ey+1

#ifdef DEBUG
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
#endif

   bd(1) = ibdry
   bd(2) = jbdry
   bd(3) = ibdry
   bd(4) = jbdry

   where(this%neighbor >= 0)
      bd = 0
   end where

#ifdef DEBUG
   call MPI_COMM_RANK(comm2d,prank,ierr)
   call MPI_COMM_SIZE(comm2d,psize,ierr)
   do iproc=0, psize-1
      if (iproc == prank) then
         print *, ' Rank ', iproc
         print *, ' neighbor ', this%neighbor
         print *, ' bd       ', bd
         print *, ' vbc ', vbc
         write(IOUT,*) 'sx=',sx,' ex=',ex,' sy=',sy,' ey=',ey
         call MPI_Barrier(MPI_COMM_WORLD, ierr)
      end if
   enddo
#endif

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

#ifdef DEBUG
110 format(/,'ERROR: process:',i3,' nxdim=',i4,' < ex-sx+3=',i4,/, &
    &       ' -> put the parameter formula for nxdim in main.F in ', &
    &       'comments and',/,'    assign to nxdim the maximum ', &
    &       'value of ex-sx+3',/)
120 format(/,'ERROR: process:',i3,' nydim=',i4,' < ey-sy+3=',i4,/, &
    &       ' -> put the parameter formula for nydim in main.F in ', &
    &       'comments and'/,'     assign to nydim the maximum ', &
    &       'value of ey-sy+3',/)
#endif

end subroutine initialize_multigrid_2d

subroutine solve_multigrid_2d(this, p, f, r)

   type(multigrid_2d)         :: this
   sll_real64, dimension(:,:) :: p
   sll_real64, dimension(:,:) :: f
   sll_real64, dimension(:,:) :: r
   sll_int32                  :: iter
   sll_real64, parameter      :: rro = 0.0_f64 !Average density

!-----------------------------------------------------------------------
! solve using mgd2
!
   call mgdsolver(1,             &
  &               this%sx,       &
                  this%ex,       &
                  this%sy,       &
                  this%ey,       &
                  p,             &
                  f,             &   
                  r,             &
                  ngrid,         &
                  this%work,     &
                  maxcy,         &
                  tolmax,        &
                  kcycle,        &
                  iprer,         &
                  ipost,         &
                  iresw,         &
                  this%xl,       &
                  this%yl,       &
                  rro,           &
                  nx,            &
                  ny,            &
                  this%comm2d,   &
                  myid,          &
                  this%neighbor, &
                  bd,            &
                  phibc,         &
                  iter,          &
                  .true.,        &
                  IOUT,          &
                  nerror)


end subroutine solve_multigrid_2d

subroutine mpe_decomp1d(n,numprocs,myid,s,e)

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
   s       = myid * nlocal + 1
   deficit = mod(n,numprocs)
   s       = s + min(myid,deficit)
   if (myid  < deficit) then
       nlocal = nlocal + 1
   endif
   e = s + nlocal - 1
   if (e  >  n .or. myid == numprocs-1) e = n

end subroutine mpe_decomp1d

subroutine write_topology( this )

   type(multigrid_2d) :: this

   sll_int32  :: coords(2)
   sll_real64 :: dx
   sll_real64 :: dy
   sll_int32  :: prank
   sll_int32  :: psize
   sll_int32  :: code
   sll_int32  :: iproc, jproc
   sll_int32  :: file_id
   sll_int32  :: error

   call MPI_COMM_RANK(this%comm2d,prank,code)
   call MPI_COMM_SIZE(this%comm2d,psize,code)

   dx = 1.0_f64 / this%pdims(1) * this%xl
   dy = 1.0_f64 / this%pdims(2) * this%yl
   
   if ( prank == 0) then
   
      call sll_ascii_file_create("mpi_topology.mtv", file_id, error)
      write(file_id,*)"$DATA=CONTOUR"
      write(file_id,*)"%equalscale=T"
      write(file_id,*)"%interp     = 0"
      write(file_id,*)"%contstyle  = 0"
      write(file_id,*)"%meshplot   = on" 
      write(file_id,*)"%hiddenline = off" 
      write(file_id,*)"%xmin=",0., " xmax = ", this%xl
      write(file_id,*)"%ymin=",0., " ymax = ", this%yl
      write(file_id,*)"% nx   = ", this%pdims(1)+1
      write(file_id,*)"% ny   = ", this%pdims(2)+1
      do jproc = 1, this%pdims(2)+1
         write(file_id,*)(0., iproc=1,this%pdims(1)+1)
      end do

      do iproc = 0, psize-1
         CALL MPI_CART_COORDS(this%comm2d,iproc,2,coords,error)
         write(file_id,111)coords(1)+.5*dx,coords(2)+.5*dy,iproc
      end do
   
      write(file_id,"('$END')")
      close(file_id)

   end if

111 format("@text x1=",f7.3,' y1=',f7.3," z1=0.0 linelabel='",i4,"'")

end subroutine write_topology

end module sll_multigrid_2d
