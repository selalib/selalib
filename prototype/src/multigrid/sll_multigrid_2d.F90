module sll_multigrid_2d
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_boundary_condition_descriptors.h"
use mpi
use diagnostics
#include "compdir.inc"

type :: mgd2_solver

   sll_int32  :: ngrid
   sll_int32  :: sx
   sll_int32  :: ex
   sll_int32  :: sy
   sll_int32  :: ey
   sll_int32  :: nwork 
   sll_int32  :: nx
   sll_int32  :: ny
   sll_int32  :: comm2d
   sll_int32  :: neighbor(8)
   sll_int32  :: bd(8)
   sll_int32  :: phibc(4,20)

   sll_real64 :: xl
   sll_real64 :: yl

   sll_real64, pointer :: work(:)

end type mgd2_solver

integer :: ierr

contains

!> initialize mgd2 solver
subroutine initialize( this, &
                       x_min, x_max, nx, nxprocs, bc_x,  &
                       y_min, y_max, ny, nyprocs, bc_y,  &
                       vbc_x, vbc_y  )

   type(mgd2_solver) :: this

   sll_int32, intent(in)  :: nx
   sll_int32, intent(in)  :: ny
   sll_real64, intent(in) :: x_min
   sll_real64, intent(in) :: x_max
   sll_real64, intent(in) :: y_min
   sll_real64, intent(in) :: y_max
   sll_int32, intent(in)  :: bc_x
   sll_int32, intent(in)  :: bc_y
   sll_int32, intent(in)  :: nxprocs
   sll_int32, intent(in)  :: nyprocs
   sll_real64, optional   :: vbc_x(2)
   sll_real64, optional   :: vbc_y(2)

   sll_int32 :: i
   sll_int32 :: myid
   sll_int32 :: numprocs
   sll_int32 :: iex
   sll_int32 :: jey
   sll_int32 :: sx
   sll_int32 :: ex
   sll_int32 :: sy
   sll_int32 :: ey
   sll_int32 :: ibdry
   sll_int32 :: jbdry
   sll_int32 :: nerror
   sll_int32 :: ierr
   sll_int32 :: neighbor(8)
   sll_int32 :: dims(2)
   sll_int32 :: coords(2)
   sll_int32 :: status(MPI_STATUS_SIZE)
   sll_int32 :: nbrright,nbrbottom,nbrleft,nbrtop
   sll_int32 :: nxdim
   sll_int32 :: nydim

   sll_int32, parameter :: iout = 6
   sll_int32, parameter :: ixp=4, jyq=4

   logical   :: periods(2)
   sll_real64 :: vbc(4)


   call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)

   iex = ceiling(log((nx-1.)/ixp)/log(2.))+1
   jey = ceiling(log((ny-1.)/jyq)/log(2.))+1

   this%nx = nx
   this%ny = ny
   this%xl = x_max - x_min
   this%yl = y_max - y_min

   this%nwork = (4*nx*ny*8)/(3*nxprocs*nyprocs)+(64*(nx+ny))/3+(32*4)/3

   SLL_CLEAR_ALLOCATE(this%work(1:this%nwork), ierr)

   this%ngrid = max(iex,jey)

   nxdim=int(float(nx)/float(nxprocs)+0.99)+2
   nydim=int(float(ny)/float(nyprocs)+0.99)+2

   if (numprocs.ne.(nxprocs*nyprocs)) then
      write(IOUT,"(/,'ERROR: numprocs <> (nxprocs*nyprocs)',/)")
      call MPI_FINALIZE(ierr)
   end if

   ! create Cartesian topology with periodic boundary conditions

   dims(1)    = nxprocs
   dims(2)    = nyprocs
   periods(1) = .true.
   periods(2) = .true.
   call MPI_CART_CREATE(MPI_COMM_WORLD,2,dims,periods,.true.,this%comm2d,ierr)

   ibdry    = 0
   jbdry    = 0
   vbc(1:4) = 0.0d0
   this%bd(1:8)  = 0

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

   call MPI_CART_SHIFT(this%comm2d,0,1,nbrleft,nbrright,ierr)
   call MPI_CART_SHIFT(this%comm2d,1,1,nbrbottom,nbrtop,ierr)
   neighbor(1)=nbrright
   neighbor(3)=nbrbottom
   neighbor(5)=nbrleft
   neighbor(7)=nbrtop
   call MPI_BARRIER(this%comm2d,ierr)
   call MPI_SENDRECV(neighbor(3),1,MPI_INTEGER,nbrright,0, &
   &                  neighbor(4),1,MPI_INTEGER,nbrleft,0, &
   &                  this%comm2d,status,ierr)
   call MPI_SENDRECV(neighbor(7),1,MPI_INTEGER,nbrright,1, &
   &                  neighbor(6),1,MPI_INTEGER,nbrleft,1, &
   &                  this%comm2d,status,ierr)
   call MPI_SENDRECV(neighbor(3),1,MPI_INTEGER,nbrleft,0, &
   &                  neighbor(2),1,MPI_INTEGER,nbrright,0, &
   &                  this%comm2d,status,ierr)
   call MPI_SENDRECV(neighbor(7),1,MPI_INTEGER,nbrleft,1, &
   &                  neighbor(8),1,MPI_INTEGER,nbrright,1, &
   &                  this%comm2d,status,ierr)

   !-----------------------------------------------------------------------
   ! find indices of subdomain and check that dimensions of arrays are
   ! sufficient
   !

   call MPI_CART_GET(this%comm2d,2,dims,periods,coords,ierr)
   call MPE_DECOMP1D(nx,dims(1),coords(1),sx,ex)
   sx=sx+1
   ex=ex+1
   call MPE_DECOMP1D(ny,dims(2),coords(2),sy,ey)
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

   this%sx = sx
   this%ex = ex
   this%sy = sy
   this%ey = ey

   do i=1,8
      write(IOUT,*) 'neighbor: ',neighbor(i),' bd: ',this%bd(i)
   end do

   this%neighbor = neighbor

110   format(/,'ERROR: process:',i3,' nxdim=',i4,' < ex-sx+3=',i4,/, &
     &       ' -> put the parameter formula for nxdim in main.F in ', &
     &       'comments and',/,'    assign to nxdim the maximum ', &
     &       'value of ex-sx+3',/)
120   format(/,'ERROR: process:',i3,' nydim=',i4,' < ey-sy+3=',i4,/, &
     &       ' -> put the parameter formula for nydim in main.F in ', &
     &       'comments and'/,'     assign to nydim the maximum ', &
     &       'value of ey-sy+3',/)

   if(bc_x == SLL_PERIODIC) then
      ibdry  = 0    ! Periodic bc in direction x
   else if(bc_x == SLL_DIRICHLET) then
      ibdry  = 1    ! Periodic bc in direction x
      if (present(vbc_x)) then
         vbc(1:2) = vbc_x ! Value for x_min, x_max
      else
         vbc(1:2) = 0.0_f64
      end if
   else
      print*, " multigrid solver wrong boundary conditions in x"
   end if
   
   if(bc_y == SLL_PERIODIC) then
      jbdry  = 0    ! Periodic bc in direction y
   else if(bc_y == SLL_DIRICHLET) then
      jbdry  = 1    ! Periodic bc in direction y
      if (present(vbc_y)) then
         vbc(3:4) = vbc_y ! Value for y_min, y_max
      else
         vbc(3:4) = 0.0_f64
      end if
   else
      print*, " multigrid solver wrong boundary conditions in y"
   end if

   call mgdinit(   vbc,              &
                   this%phibc,       &
                   ixp,              &
                   jyq,              &
                   iex,              &
                   jey,              &
                   this%ngrid,       &
                   nx+2,             &
                   ny+2,             &
                   this%sx,          &
                   this%ex,          &
                   this%sy,          &
                   this%ey,          &
                   MPI_REAL8,        &
                   nxprocs,          &
                   nyprocs,          &
                   this%nwork,       &
                   ibdry,            &
                   jbdry,            &
                   myid,             &
                   IOUT,             &
                   nerror)

   return

   1000 call mpi_finalize(ierr)

   stop

end subroutine initialize


!> solve Poisson using mgd2
subroutine solve(this, p, f, r)

   type(mgd2_solver) :: this
   sll_real64, dimension(:,:) :: p
   sll_real64, dimension(:,:) :: f
   sll_real64, dimension(:,:) :: r
   sll_int32,  parameter      :: iout   = 6
   sll_int32,  parameter      :: maxcy  = 300
   sll_int32,  parameter      :: kcycle = 1
   sll_int32,  parameter      :: iprer  = 2
   sll_int32,  parameter      :: iresw  = 1
   sll_int32,  parameter      :: ipost  = 1
   sll_real64, parameter      :: tolmax = 1.0d-05
   sll_real64, parameter      :: rro    = 1.0_f64
   sll_int32                  :: myid
   sll_int32                  :: ierr
   sll_int32                  :: iter
   sll_int32                  :: nerror

   call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)

   call mgdsolver(2,               &
                  this%sx,         &
                  this%ex,         &
                  this%sy,         &
                  this%ey,         &
                  p,               &
                  f,               &
                  r,               &
                  this%ngrid,      &
                  this%work,       &
                  maxcy,           &
                  tolmax,          &
                  kcycle,          &
                  iprer,           &
                  ipost,           &
                  iresw,           &
                  this%xl,         &
                  this%yl,         &
                  rro,             &
                  this%nx,         &
                  this%ny,         &
                  this%comm2d,     &
                  myid,            &
                  this%neighbor,   &
                  this%bd,         &
                  this%phibc,      &
                  iter,            &
                  .true.,          &
                  IOUT,            &
                  nerror)

end subroutine solve

end module sll_multigrid_2d
