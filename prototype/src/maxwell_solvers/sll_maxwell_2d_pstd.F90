!------------------------------------------------------------------------------
! SELALIB
!------------------------------------------------------------------------------
!
!> @namespace sll_maxwelle_2d_pstd
!
!> @author
!> Pierre Navaro
!>
!
! DESCRIPTION: 
!
!> @brief
!> Implements the Maxwell solver in 2D with periodic boundary conditions
!>
!>@details
!>This module depends on:
!> - memory
!> - precision
!> - assert 
!> - numerical_utilities
!> - constants
!> - mesh_types
!> - diagnostics
!> - sll_utilities
!>
! REVISION HISTORY:
! 09 01 2012 - Initial Version
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!------------------------------------------------------------------------------


#define FFTW_ALLOCATE(array,sz_array,p_array)                             \
self%sz_array = int((nx/2+1)*ny,C_SIZE_T);               \
self%p_array = fftw_alloc_complex(self%sz_array);        \
call c_f_pointer(self%p_array, self%array, [nx/2+1,ny])  \


program sll_maxwell_2d_pstd

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

use numeric_constants

use, intrinsic :: iso_c_binding
implicit none
include 'fftw3.f03'

type :: maxwell
   sll_int32                                          :: nx
   sll_int32                                          :: ny
   sll_real64                                         :: dx
   sll_real64                                         :: dy
   sll_real64, dimension(:,:), allocatable            :: kx
   sll_real64, dimension(:,:), allocatable            :: ky
   type(C_PTR)                                        :: fw
   type(C_PTR)                                        :: bw
   complex(C_DOUBLE_COMPLEX), dimension(:,:), pointer :: hxt
   complex(C_DOUBLE_COMPLEX), dimension(:,:), pointer :: hyt
   complex(C_DOUBLE_COMPLEX), dimension(:,:), pointer :: hzt
   complex(C_DOUBLE_COMPLEX), dimension(:,:), pointer :: ext
   complex(C_DOUBLE_COMPLEX), dimension(:,:), pointer :: eyt
   complex(C_DOUBLE_COMPLEX), dimension(:,:), pointer :: ezt
   integer(C_SIZE_T)                                  :: sz_hxt
   integer(C_SIZE_T)                                  :: sz_hyt
   integer(C_SIZE_T)                                  :: sz_ezt
   integer(C_SIZE_T)                                  :: sz_ext
   integer(C_SIZE_T)                                  :: sz_eyt
   integer(C_SIZE_T)                                  :: sz_hzt
   type(C_PTR)                                        :: p_hxt
   type(C_PTR)                                        :: p_hyt
   type(C_PTR)                                        :: p_hzt
   type(C_PTR)                                        :: p_ext
   type(C_PTR)                                        :: p_eyt
   type(C_PTR)                                        :: p_ezt
end type maxwell

sll_int32, parameter :: nthreads = 2

sll_real64, dimension(:,:), allocatable :: x
sll_real64, dimension(:,:), allocatable :: y

sll_real64, dimension(:,:), allocatable :: hx
sll_real64, dimension(:,:), allocatable :: hy
sll_real64, dimension(:,:), allocatable :: ez

sll_real64, dimension(:,:), allocatable :: ex
sll_real64, dimension(:,:), allocatable :: ey
sll_real64, dimension(:,:), allocatable :: hz
sll_real64 :: xmin, xmax, ymin, ymax
sll_real64 :: omega, dt, cfl, time, dx, dy, dimx, dimy, tfinal
sll_int32  :: nstep, istep, iplot, nstepmax, error
sll_int32  :: nx, ny, i, j, md=2, nd=2, ik, jk

type(maxwell) :: this

nx = 128
ny = 128

SLL_ALLOCATE(ex(nx,ny), error)
SLL_ALLOCATE(ey(nx,ny), error)
SLL_ALLOCATE(ez(nx,ny), error)
SLL_ALLOCATE(hx(nx,ny), error)
SLL_ALLOCATE(hx(nx,ny), error)
SLL_ALLOCATE(hz(nx,ny), error)

xmin = 0.
xmax = sll_pi
ymin = 0.
ymax = sll_pi

dimx = xmax - xmin
dimy = ymax - ymin

call init_maxwell_solver(this, xmin, xmax, nx, &
                               ymin, ymax, ny, &
                               hx, hy, ez, error )                                

do j = 1, ny
   do i = 1, nx
      x(i,j) = xmin + (i-0.5) / (xmax - xmin)
      y(i,j) = ymin + (j-0.5) / (ymax - ymin)
   end do
end do

dt = cfl  / sqrt (1./(dx*dx)+1./(dy*dy)) 
nstep = floor(tfinal/dt)
write(*,*)
write(*,*) " dx = ", dx
write(*,*) " dy = ", dy
write(*,*) " dt = ", dt
write(*,*)
if( nstep > nstepmax ) nstep = nstepmax
write(*,*) " Nombre d'iteration nstep = ", nstep

time  = 0.
iplot = 0

omega = sqrt(2._f64)

ex = 0.0d0
ey = 0.0d0
hz = - cos(x)*cos(y)*cos(omega*(-0.5*dt))

!Polarisation TE
!ex =  cos(x)*sin(y)*sin(omega*time)/omega
!ey = -sin(x)*cos(y)*sin(omega*time)/omega
!hz = -cos(x)*cos(y)*cos(omega*time)
!
!Polarisation TM
!ez =  cos(x)*cos(y)*cos(omega*time)
!hx =  cos(x)*sin(y)*sin(omega*time)/omega
!hy = -sin(x)*cos(y)*sin(omega*time)/omega

contains

subroutine init_maxwell_solver(self, xmin, xmax, nx, &
                                     ymin, ymax, ny, &
                                     hx, hy, ez, error )
   type(maxwell)                         :: self
   sll_real64                                :: xmin
   sll_real64                                :: xmax
   sll_real64                                :: ymin
   sll_real64                                :: ymax
   sll_int32                                 :: nx
   sll_int32                                 :: ny
   sll_real64, dimension(:,:), intent(inout) :: hx
   sll_real64, dimension(:,:), intent(inout) :: hy
   sll_real64, dimension(:,:), intent(inout) :: ez
   !sll_real64, dimension(:,:), intent(inout) :: ex
   !sll_real64, dimension(:,:), intent(inout) :: ey
   !sll_real64, dimension(:,:), intent(inout) :: hz
   sll_int32                                 :: error
   sll_real64                                :: dx
   sll_real64                                :: dy
   sll_real64                                :: kx0, kx
   sll_real64                                :: ky0, ky

   self%dx = (xmax-xmin) / (nx-1d0)
   self%dy = (ymax-ymin) / (ny-1d0)

   FFTW_ALLOCATE(hxt,sz_hxt,p_hxt)
   FFTW_ALLOCATE(hyt,sz_hyt,p_hyt)
   FFTW_ALLOCATE(ezt,sz_ezt,p_ezt)

   call dfftw_init_threads(error)
   if (error == 0) stop 'FFTW CAN''T USE THREADS'
   call dfftw_plan_with_nthreads(nthreads)
   
   self%fw = fftw_plan_dft_r2c_2d(ny, nx, ez, self%ezt, FFTW_ESTIMATE)
   self%bw = fftw_plan_dft_c2r_2d(ny, nx, self%ezt, ez, FFTW_ESTIMATE)

   SLL_ALLOCATE(self%kx(nx/2+1,ny), error)
   SLL_ALLOCATE(self%ky(nx/2+1,ny), error)
   
   dx  = self%dx
   dy  = self%dy
   kx0 = 2._f64*sll_pi/(nx*dx)
   ky0 = 2._f64*sll_pi/(ny*dy)

   do ik=1,nx/2+1
      kx = (ik-1)*kx0
      do jk = 1, ny/2
         ky = (jk-1)*ky0
         self%kx(jk,ik) = kx
         self%ky(jk,ik) = ky
      end do
      do jk = ny/2+1 , ny     
         ky = (jk-1-ny)*ky0
         self%kx(jk,ik) = kx
         self%ky(jk,ik) = ky
      end do
   end do
   self%kx(1,1) = 1.0_f64

end subroutine init_maxwell_solver

!> Solve faraday 
subroutine solve_faraday(self, hx, hy, ez)

   type(maxwell),intent(inout)  :: self
   sll_real64, dimension(:,:), intent(inout) :: hx
   sll_real64, dimension(:,:), intent(out)   :: hy
   sll_real64, dimension(:,:), intent(out)   :: ez
   sll_int32  :: nx, ny
   sll_real64 :: dx, dy

   nx = self%nx
   ny = self%ny

   call fftw_execute_dft_r2c(self%fw, ez, self%ezt)

   call fftw_execute_dft_c2r(self%bw, self%ezt, ez)

   hx(i,j) = hx(i,j) - dt 

   self%ext(1,1) = 0.0_f64
   self%eyt(1,1) = 0.0_f64
   self%ext = -cmplx(0.0_f64,self%kx,kind=f64)*self%ext
   self%eyt = -cmplx(0.0_f64,self%ky,kind=f64)*self%eyt



end subroutine solve_faraday

!> Solve ampere
subroutine solve_ampere_maxwell(self, hx, hy, ez)

   type(maxwell),intent(inout)  :: self
   sll_int32                    :: nx
   sll_int32                    :: ny
   sll_real64, dimension(:,:)   :: hx
   sll_real64, dimension(:,:)   :: hy
   sll_real64, dimension(:,:)   :: ez

   call fftw_execute_dft_r2c(self%fw, hx, self%hxt)
   call fftw_execute_dft_r2c(self%fw, hy, self%hyt)


   call fftw_execute_dft_c2r(self%bw, self%hxt, hx)
   call fftw_execute_dft_c2r(self%bw, self%hyt, hy)
   call fftw_execute_dft_c2r(self%bw, self%ezt, ez)

   nx = self%nx
   ny = self%ny

   hx = hx / (nx*ny)     ! normalize
   hy = hy / (nx*ny)     ! normalize
   hz = hz / (nx*ny)     ! normalize

end subroutine solve_ampere_maxwell


subroutine free_maxwell(self)
type(maxwell) :: self
sll_int32       :: error

if (c_associated(self%p_hxt)) call fftw_free(self%p_hxt)
if (c_associated(self%p_hyt)) call fftw_free(self%p_hyt)
if (c_associated(self%p_ezt)) call fftw_free(self%p_ezt)

call dfftw_destroy_plan(self%fw)
call dfftw_destroy_plan(self%bw)
if (nthreads > 1) then
   call dfftw_cleanup_threads(error)
end if

end subroutine free_maxwell

end program sll_maxwell_2d_pstd
