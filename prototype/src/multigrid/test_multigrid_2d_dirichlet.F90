#define MPI_MASTER 0

program test_multigrid_2d_dirichlet

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
!-----------------------------------------------------------------------

   sll_real64, dimension(:,:), allocatable :: p
   sll_real64, dimension(:,:), allocatable :: q
   sll_real64, dimension(:,:), allocatable :: f
   sll_real64, dimension(:,:), allocatable :: r
   sll_real64, dimension(:,:), allocatable :: x
   sll_real64, dimension(:,:), allocatable :: y

   sll_real64 :: err, errloc

   type(multigrid_2d) :: solver

   sll_int32 :: prank
   sll_int32 :: narg
   sll_int32 :: psize
   sll_int32 :: error
   sll_int32 :: i, j, sx, ex, sy, ey

!
! parameters
!
   sll_int32             :: nxprocs
   sll_int32             :: nyprocs 
   sll_int32,  parameter :: nx = 64
   sll_int32,  parameter :: ny = 64
   sll_real64, parameter :: x_min = -1.0_f64
   sll_real64, parameter :: x_max = +1.0_f64
   sll_real64, parameter :: y_min = -1.0_f64
   sll_real64, parameter :: y_max = +1.0_f64
   character(len=4)      :: buffer

   sll_real64 :: delta_x
   sll_real64 :: delta_y

   sll_real64 :: start_time, end_time, total_time

   call sll_boot_collective() 

   prank = sll_get_collective_rank(sll_world_collective)
   psize = sll_get_collective_size(sll_world_collective)

   narg = command_argument_count()

   if (narg /= 2) then
      if (prank == MPI_MASTER) then
         print*, 'give nxprocs and nyprocs'
         print*, 'try'
         print*, 'mpirun -np 4 test_multigrid_2d_dirichlet 2 2 '
      end if
      call sll_halt_collective()
      stop
   end if

   call get_command_argument(1, buffer); read(buffer,"(i4)") nxprocs
   call get_command_argument(2, buffer); read(buffer,"(i4)") nyprocs

   start_time = MPI_WTIME() 

   call initialize( solver,                               &
                    x_min, x_max, y_min, y_max,           &
                    nxprocs, nyprocs, nx, ny,             &
                    SLL_DIRICHLET, SLL_DIRICHLET,         &
                    0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64)

   call write_topology( solver )

   sx = solver%sx
   ex = solver%ex
   sy = solver%sy
   ey = solver%ey

   SLL_CLEAR_ALLOCATE(p(sx-1:ex+1,sy-1:ey+1),error)
   SLL_CLEAR_ALLOCATE(q(sx-1:ex+1,sy-1:ey+1),error)
   SLL_CLEAR_ALLOCATE(f(sx-1:ex+1,sy-1:ey+1),error)
   SLL_CLEAR_ALLOCATE(r(sx-1:ex+1,sy-1:ey+1),error)
   SLL_CLEAR_ALLOCATE(x(sx-1:ex+1,sy-1:ey+1),error)
   SLL_CLEAR_ALLOCATE(y(sx-1:ex+1,sy-1:ey+1),error)

   ! initialize problem
   ! xl,yl are the dimensions of the domain
   ! wk is the wavenumber (must be an integer value)
   ! rro is the average density
   ! 1/hxi,1/hyi are the spatial resolutions

   !-----------------------------------------------------------------------
   ! Initialize the pressure, density, and right-hand side of the
   ! elliptic equation div(1/r*grad(p))=f
   !-----------------------------------------------------------------------

   r    =  1.0_f64

   delta_x = (x_max - x_min) / nx
   delta_y = (y_max - y_min) / ny

   do j=sy-1, ey+1
     do i=sx-1, ex+1
       x(i,j) = x_min + (i-1.5d0)*delta_x
       y(i,j) = y_min + (j-1.5d0)*delta_y
     end do
   end do

   q = exp(-(x*x+y*y)/0.04_f64)

   do j=sy, ey
     do i=sx, ex
        f(i,j) = (q(i-1,j)-2.*q(i,j)+q(i+1,j))/(delta_x*delta_x) + &
                 (q(i,j-1)-2.*q(i,j)+q(i,j+1))/(delta_y*delta_y)
     end do
   end do

   p = 0.0_f64
   call gp_plot2d(f,x,y,sx,ex,sy,ey,'f')

   call solve(solver, p, f, r)

   end_time = MPI_WTIME() 

   call MPI_ALLREDUCE(end_time-start_time,total_time,1,MPI_REAL8,MPI_SUM, &
                      solver%comm2d,error)

   call gp_plot2d(p,  x,y,sx,ex,sy,ey,'p')
   call gp_plot2d(q,  x,y,sx,ex,sy,ey,'q')
   call gp_plot2d(p-q,x,y,sx,ex,sy,ey,'e')
   call gp_plot2d(f,  x,y,sx,ex,sy,ey,'r')

   !-----------------------------------------------------------------------
   ! Calculate the error between the numerical and exact solution to
   ! the test problem.
   !-----------------------------------------------------------------------

   errloc=sum(abs(p(sx:ex,sy:ey)-q(sx:ex,sy:ey)))
   
   ! calculate global error
   call MPI_ALLREDUCE(errloc,err,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
                      solver%comm2d,error)

   write(*,100) errloc/float(nx*ny),err/float(nx*ny)

   call flush(6)
   if (prank == MPI_MASTER) &
   write(*,"('proc:',i3,' time: ', f12.4)") prank, total_time
   call flush(6)
   call sll_halt_collective()

   stop

100 format(/,'Local error: ',e13.6,'  total error: ',e13.6,/)

end program test_multigrid_2d_dirichlet
