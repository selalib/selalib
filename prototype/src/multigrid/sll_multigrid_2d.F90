module sll_multigrid_2d
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_boundary_condition_descriptors.h"

   use sll_collective
   use diagnostics

   implicit none

   sll_real64, dimension(:),   allocatable :: work
   sll_int32  :: bd(8)
   sll_int32  :: nx
   sll_int32  :: ny
   sll_real64 :: phibc(4,20)
   sll_real64 :: wk
   sll_int32  :: comm2d
   sll_real64, parameter :: tolmax  = 1.0d-05
   sll_int32  :: sx
   sll_int32  :: ex
   sll_int32  :: sy
   sll_int32  :: ey
   sll_int32,  parameter :: IOUT    = 6
   sll_int32  :: myid
   sll_int32  :: neighbor(8)
   sll_int32  :: ngrid
   sll_int32  :: nerror

contains

subroutine initialize( nxprocs, nyprocs, nxdim, nydim  )

   sll_int32, intent(in)  :: nxprocs
   sll_int32, intent(in)  :: nyprocs

   sll_int32,  parameter :: ixp     = 4
   sll_int32,  parameter :: jyq     = 4

!
! variables
!
   logical    :: periods(2)

   sll_int32  :: iex, jey

   sll_int32  :: numprocs
   sll_int32  :: ibdry
   sll_int32  :: jbdry
   sll_int32  :: dims(2),coords(2),i
   sll_int32  :: status(MPI_STATUS_SIZE)
   sll_int32  :: nbrright
   sll_int32  :: nbrbottom
   sll_int32  :: nbrleft
   sll_int32  :: nbrtop
   sll_int32  :: nwork
   sll_int32  :: nxdim
   sll_int32  :: nydim
   sll_int32  :: ierr

   sll_real64 :: vbc(4)



   iex   = ceiling(log((nx-1.)/ixp)/log(2.))+1
   jey   = ceiling(log((ny-1.)/jyq)/log(2.))+1
   ngrid = max(iex,jey)

   nwork = (4*nx*ny*8)/(3*nxprocs*nyprocs)+(64*(nx+ny))/3+(32*4)/3
   nxdim = int(float(nx)/float(nxprocs)+0.99)+2
   nydim = int(float(ny)/float(nyprocs)+0.99)+2

   SLL_CLEAR_ALLOCATE(work(1:nwork),ierr)

!-----------------------------------------------------------------------
! initialize MPI and create a datatype for real numbers
!
   call MPI_INIT(ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)

   if (numprocs.ne.(nxprocs*nyprocs)) then
     write(IOUT,"(/,'ERROR: numprocs <> (nxprocs*nyprocs)',/)")
     call mpi_finalize(ierr)
     stop
   end if

   call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)

   if (numprocs.ne.(nxprocs*nyprocs)) then
      write(IOUT,"(/,'ERROR: numprocs <> (nxprocs*nyprocs)',/)")
      call MPI_FINALIZE(ierr)
   end if

   ! create Cartesian topology with periodic boundary conditions

   dims(1)    = nxprocs
   dims(2)    = nyprocs
   periods(1) = .true.
   periods(2) = .true.
   call MPI_CART_CREATE(MPI_COMM_WORLD,2,dims,periods,.true.,comm2d,ierr)

   ibdry    = 0
   jbdry    = 0
   vbc(1:4) = 0.0d0
   bd(1:8)  = 0

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
   call MPE_DECOMP1D(nx,dims(1),coords(1),sx,ex)
   sx=sx+1
   ex=ex+1
   call MPE_DECOMP1D(ny,dims(2),coords(2),sy,ey)
   sy=sy+1
   ey=ey+1
   if ((ex-sx+3).gt.nxdim) then
      write(IOUT,110) myid,nxdim,ex-sx+3
      nerror=1
      write(IOUT,"(/,'ERROR in multigrid code : 178',/)")
   end if
   if ((ey-sy+3).gt.nydim) then
      write(IOUT,120) myid,nydim,ey-sy+3
      nerror=1
      write(IOUT,"(/,'ERROR in multigrid code : 183',/)")
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

   ibdry  = 0    ! Periodic bc in direction x
   jbdry  = 0    ! Periodic bc in direction y

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

end subroutine initialize

subroutine solve(p, f, r)
   sll_real64, dimension(:,:) :: p
   sll_real64, dimension(:,:) :: f
   sll_real64, dimension(:,:) :: r

   sll_real64 :: hxi
   sll_real64 :: hyi
   sll_real64 :: xl
   sll_real64 :: yl
   sll_real64 :: rro
   sll_int32  :: iter
   sll_int32,  parameter :: kcycle  = 1
   sll_int32,  parameter :: iprer   = 2
   sll_int32,  parameter :: ipost   = 1
   sll_int32,  parameter :: iresw   = 1
   sll_int32,  parameter :: maxcy   = 300
!-----------------------------------------------------------------------
! initialize problem
! xl,yl are the dimensions of the domain
! wk is the wavenumber (must be an integer value)
! rro is the average density
! 1/hxi,1/hyi are the spatial resolutions
!
   xl  = 1.0d0
   yl  = 1.0d0
   wk  = 2.0d0
   rro = 1.0d0
   hxi = float(nx)/xl
   hyi = float(ny)/yl
   write(IOUT,*) 'hxi=',hxi,' hyi=',hyi

   call ginit(sx,ex,sy,ey,p,r,f,wk,hxi,hyi,sll_pi)

!-----------------------------------------------------------------------
! solve using mgd2
!
   call mgdsolver(2,sx,ex,sy,ey,p,f,r,ngrid,work, &
  &               maxcy,tolmax,kcycle,iprer,ipost,iresw, &
  &               xl,yl,rro,nx,ny,comm2d,myid,neighbor, &
  &               bd,phibc,iter,.true.,IOUT,nerror)

   if (nerror.eq.1)  then
      write(IOUT,"(/,'ERROR in multigrid code : 261',/)")
   end if
!-----------------------------------------------------------------------
! compare numerical and exact solutions
!
   call gerr(sx,ex,sy,ey,p,comm2d,wk,hxi,hyi,sll_pi,nx,ny,IOUT)
!-----------------------------------------------------------------------
   return



end subroutine solve

end module sll_multigrid_2d
