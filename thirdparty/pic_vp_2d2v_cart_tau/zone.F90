module m_zone
#include "sll_working_precision.h"
use sll_m_working_precision
use sll_m_constants
use mpi

private

type, public :: mesh_fields
  sll_real64, dimension(:,:), pointer :: ex
  sll_real64, dimension(:,:), pointer :: ey
  sll_real64, dimension(:,:), pointer :: r0
end type mesh_fields

logical :: relativ 

sll_int32,  public :: nx         ! Number of cells along x
sll_int32,  public :: ny         ! Number of cells along y
sll_int32,  public :: nstep      ! Time step number
sll_int32,  public :: nbpart     ! Number of particles
sll_real64, public :: dimx       ! Domain length along x
sll_real64, public :: dimy       ! Domain length along y
sll_real64, public :: dx         ! dimx / nx
sll_real64, public :: dy         ! dimy / ny
sll_real64, public :: dt         ! Time step
sll_real64, public :: alpha      ! Perturbation amplitude
sll_real64, public :: kx         ! Perturbation wave number along x
sll_real64, public :: ky         ! Perturbation wave number along y
sll_int32,  public :: ntau       ! discretization in fast scale direction
sll_real64, public :: ep         ! epsilon (strength of the magnetic field)
logical,    public :: plot       ! rho plot (true or false)
sll_real64, public :: tfinal     ! time (max)
sll_real64, public :: reset      ! time period for particles reset

sll_int32  :: npm        ! Number of particles by cell
sll_int32  :: nstepmax   ! Time step number (max)

integer, dimension(MPI_STATUS_SIZE) :: stat

integer :: nproc=1, rang=0
integer :: tag=1111
integer, public :: code

public readin, init_mpi, mpi_global_master, finish_mpi

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine readin( filename )

implicit none

character(len=*) :: filename

namelist/donnees/   nx,  & ! number of points along x direction
                    ny,  & ! number of points along y direction
                tfinal,  & ! duree maxi
              nstepmax,  & ! nbre d'iterations maxi
                  ntau,  & ! number of tau iterations
                    ep,  & ! epsilon
                   npm,  & ! number of particles by cell
                    kx,  & ! wave number along x
                    ky,  & ! wave number along y
                    dt,  & ! time step
                  plot,  & ! true or false
                 alpha,  & ! amplitude of perturbation
                 reset     ! time period for particles reset

nstepmax   = 2000          ! nbre d'iterations maxi
nx         = 120           ! nombre de pts suivant x
ny         = 10            ! nombre de pts suivant y
ntau       = 32
npm        = 50
ep         = 0.1_f64
alpha      = 0.05_f64
kx         = 0.5_f64
ky         = 1.0_f64 
dt         = 0.1_f64
plot       = .false.
tfinal     = nstepmax * dt
reset      = tfinal

write(*,"(a)") "Input file name :"// filename
open(93,file=filename,status='old')
read(93,donnees) 
close(93)

dimx  = 2.0_f64 * sll_p_pi / kx
dimy  = 2.0_f64 * sll_p_pi / ky  

dx = dimx / nx
dy = dimy / ny

nbpart = npm * nx * ny

nstep = min(floor(tfinal/dt),nstepmax)

end subroutine readin

!*********************************************************************

subroutine init_mpi( prank , psize)
integer :: prank
integer :: psize
integer :: iproc

call MPI_INIT(code)
call MPI_COMM_RANK(MPI_COMM_WORLD,prank,code)
call MPI_COMM_SIZE(MPI_COMM_WORLD,psize,code)

end subroutine init_mpi

!*********************************************************************

subroutine finish_mpi()

call MPI_BARRIER(MPI_COMM_WORLD,code)
call MPI_FINALIZE(code)
stop

end subroutine finish_mpi

!*********************************************************************

subroutine mpi_global_master()

call MPI_BCAST(nx,    1,MPI_INTEGER,0,MPI_COMM_WORLD,code)
call MPI_BCAST(ny,    1,MPI_INTEGER,0,MPI_COMM_WORLD,code)
call MPI_BCAST(nstep, 1,MPI_INTEGER,0,MPI_COMM_WORLD,code)
call MPI_BCAST(nbpart,1,MPI_INTEGER,0,MPI_COMM_WORLD,code)
call MPI_BCAST(ntau,  1,MPI_INTEGER,0,MPI_COMM_WORLD,code)
call MPI_BCAST(dimx,  1,MPI_REAL8  ,0,MPI_COMM_WORLD,code)
call MPI_BCAST(dimy,  1,MPI_REAL8  ,0,MPI_COMM_WORLD,code)
call MPI_BCAST(dx,    1,MPI_REAL8  ,0,MPI_COMM_WORLD,code)
call MPI_BCAST(dy,    1,MPI_REAL8  ,0,MPI_COMM_WORLD,code)
call MPI_BCAST(alpha, 1,MPI_REAL8  ,0,MPI_COMM_WORLD,code)
call MPI_BCAST(kx,    1,MPI_REAL8  ,0,MPI_COMM_WORLD,code)
call MPI_BCAST(ky,    1,MPI_REAL8  ,0,MPI_COMM_WORLD,code)
call MPI_BCAST(dt,    1,MPI_REAL8  ,0,MPI_COMM_WORLD,code)
call MPI_BCAST(ep,    1,MPI_REAL8  ,0,MPI_COMM_WORLD,code)
call MPI_BCAST(reset, 1,MPI_REAL8  ,0,MPI_COMM_WORLD,code)

end subroutine mpi_global_master


end module m_zone
