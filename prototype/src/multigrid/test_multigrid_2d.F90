program test_multigrid_2d

use sll_collective
use diagnostics
!use sll_multigrid_2d
use sll_boundary_condition_descriptors
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_constants.h"
use sll_multigrid_2d

implicit none

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
!             mgdinit, mgdsolver, gerr
!-----------------------------------------------------------------------
!
! parameters
!
   sll_int32,  parameter :: nxprocs = 4
   sll_int32,  parameter :: nyprocs = 2

   sll_real64, dimension(:,:), allocatable :: p
   sll_real64, dimension(:,:), allocatable :: f
   sll_real64, dimension(:,:), allocatable :: r

   sll_int32 :: nxdim
   sll_int32 :: nydim

   sll_int32 :: ierr

   nx = 128
   ny = 64

   call initialize( nxprocs, nyprocs, nxdim, nydim  )

   SLL_CLEAR_ALLOCATE(p(1:nxdim,1:nydim),ierr)
   SLL_CLEAR_ALLOCATE(f(1:nxdim,1:nydim),ierr)
   SLL_CLEAR_ALLOCATE(r(1:nxdim,1:nydim),ierr)

   call solve(p, f, r)

   call gp_plot2d(p, dble(nx), dble(ny), sx, ex, sy, ey, wk )

   call mpi_finalize(nerror)
   stop

end program test_multigrid_2d

subroutine ginit(sx,ex,sy,ey,p,r,f,wk,hxi,hyi,pi)

   use sll_working_precision

   sll_int32  :: sx,ex,sy,ey
   sll_real64 :: p(sx-1:ex+1,sy-1:ey+1)
   sll_real64 :: r(sx-1:ex+1,sy-1:ey+1)
   sll_real64 :: f(sx-1:ex+1,sy-1:ey+1)
   sll_real64 :: hxi,hyi,wk,pi
   !-----------------------------------------------------------------------
   ! Initialize the pressure, density, and right-hand side of the
   ! elliptic equation div(1/r*grad(p))=f
   !
   !-----------------------------------------------------------------------
   sll_int32  :: i,j
   sll_real64 :: cnst,cx,cy,xi,yj
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

   use mpi
   use sll_working_precision

   sll_int32 :: sx,ex,sy,ey,comm2d,IOUT,nx,ny
   sll_real64 :: p(sx-1:ex+1,sy-1:ey+1),wk,hxi,hyi,pi
   !-----------------------------------------------------------------------
   ! Calculate the error between the numerical and exact solution to
   ! the test problem.
   !
   ! Code      : tmgd2
   ! Called in : main
   ! Calls     : MPI_ALLREDUCE
   !-----------------------------------------------------------------------
   sll_int32 :: i,j, ierr
   sll_real64 :: errloc,err,cx,cy,exact, xi, yj
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

   use mpi
   use sll_working_precision

   sll_int32 :: sx,ex,sy,ey,nx,ny
   sll_real64 :: a(sx-1:ex+1,sy-1:ey+1),avo,acorr
   sll_int32 :: comm2d, ierr
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
   sll_real64 :: avloc,av
   sll_int32 :: i,j
#if cdebug
   sll_real64 :: tinitial
   tinitial=MPI_WTIME()
#endif
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
#if cdebug
   nallreduce=nallreduce+1
#endif
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
   
#if cdebug
   timing(49)=timing(49)+MPI_WTIME()-tinitial
#endif
   
   
   return

end 

subroutine MPE_DECOMP1D(n,numprocs,myid,s,e)

   use sll_working_precision
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
   if (myid .lt. deficit) then
       nlocal = nlocal + 1
   endif
   e = s + nlocal - 1
   if (e .gt. n .or. myid .eq. numprocs-1) e = n

   return

end
