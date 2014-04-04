program test_multigrid_2d

use sll_collective
use sll_multigrid_2d
use sll_boundary_condition_descriptors
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_constants.h"
#include "sll_utilities.h"

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
   sll_int32,  parameter :: nxprocs = 2
   sll_int32,  parameter :: nyprocs = 2

   sll_real64, dimension(:,:), allocatable :: p
   sll_real64, dimension(:,:), allocatable :: f
   sll_real64, dimension(:,:), allocatable :: r
   sll_real64, dimension(:,:), allocatable :: x
   sll_real64, dimension(:,:), allocatable :: y

   sll_int32 :: nxdim
   sll_int32 :: nydim

   sll_real64 :: cnst,cx,cy,wk
   sll_real64 :: err, errloc

   type(multigrid_2d) :: mgd_solver

   sll_int32, parameter :: nx = 128
   sll_int32, parameter :: ny = 64

   sll_int32 :: prank
   sll_int32 :: psize
   sll_int32 :: error
   sll_int32 :: i, j

   call sll_boot_collective() 

   prank = sll_get_collective_rank(sll_world_collective)
   psize = sll_get_collective_size(sll_world_collective)

   call initialize( mgd_solver,                         &
                    0.0_f64, 1.0_f64, 0.0_f64, 1.0_f64, &
                    nxprocs, nyprocs, nx, ny, nxdim, nydim  )

   SLL_CLEAR_ALLOCATE(p(1:nxdim,1:nydim),error)
   SLL_CLEAR_ALLOCATE(f(1:nxdim,1:nydim),error)
   SLL_CLEAR_ALLOCATE(r(1:nxdim,1:nydim),error)
   SLL_CLEAR_ALLOCATE(x(1:nxdim,1:nydim),error)
   SLL_CLEAR_ALLOCATE(y(1:nxdim,1:nydim),error)

   ! initialize problem
   ! xl,yl are the dimensions of the domain
   ! wk is the wavenumber (must be an integer value)
   ! rro is the average density
   ! 1/hxi,1/hyi are the spatial resolutions

   wk  = 2.0d0

   !-----------------------------------------------------------------------
   ! Initialize the pressure, density, and right-hand side of the
   ! elliptic equation div(1/r*grad(p))=f
   !-----------------------------------------------------------------------
   !
   r    =  1.0_f64
   cnst = -8.0_f64*(sll_pi*wk)**2

   cx   = 2.0d0*sll_pi*wk
   cy   = 2.0d0*sll_pi*wk

   do j=1, nydim
     do i=1, nxdim
       x(i,j) =(float(mgd_solver%sx+i)-1.5d0)/nx
       y(i,j) =(float(mgd_solver%sy+j)-1.5d0)/ny
     end do
   end do

   f = cnst*sin(cx*x)*sin(cy*y)

   call solve(mgd_solver, p, f, r)

   call gp_plot2d()
   ! compare numerical and exact solutions

   !-----------------------------------------------------------------------
   ! Calculate the error between the numerical and exact solution to
   ! the test problem.
   !
   ! calculate local error
   !

   errloc=sum(abs(p-sin(cx*x)*sin(cy*y)))
   
   ! calculate global error
   !
   call MPI_ALLREDUCE(errloc,err,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
        &                   mgd_solver%comm2d,error)

   write(*,100) errloc/float(nx*ny),err/float(nx*ny)
   
   call sll_halt_collective()

   stop

100 format(/,'Local error: ',e13.6,'  total error: ',e13.6,/)

contains

   subroutine gp_plot2d()

   sll_int32  :: iproc
   character(len=4) :: crank

   !write domains
   call int2string(prank, crank)
   open( 80, file = "p"//crank//".dat" )
      do i=1,nxdim
         do j=1,nydim
            write(80,"(4e15.5)") x(i,j), y(i,j), p(i,j)
         end do
         write(80,*) 
      end do
   close(80)
   
   !write master file
   if (prank == 0) then
      open( 90, file = 'p.gnu')
      write(90,"(a)",advance='no')"splot 'p0000.dat' w lines"
      do iproc = 1, psize - 1
         call int2string(iproc, crank)
         write(90,"(a)",advance='no') ",'p"//crank//".dat' w lines" 
      end do
      write(90,*)
      close(90)
   end if

end subroutine gp_plot2d

end program test_multigrid_2d


