!------------------------------------------------------------------------------
! SELALIB
!------------------------------------------------------------------------------
!
!> @namespace sll_maxwell_2d_pstd
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


#define FFTW_ALLOCATE(array,array_size,sz_array,p_array) \
sz_array = int((array_size/2+1),C_SIZE_T);       \
p_array = fftw_alloc_complex(sz_array);        \
call c_f_pointer(p_array, array, [array_size/2+1])  \


module sll_maxwell_2d_pstd

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

use, intrinsic :: iso_c_binding
use numeric_constants

implicit none

interface initialize
 module procedure new_maxwell_2d_pstd
end interface
interface solve
 module procedure solve_maxwell_2d_pstd
end interface
interface free
 module procedure free_maxwell_2d_pstd
end interface

type, public :: maxwell_pstd
   sll_int32                                          :: nx
   sll_int32                                          :: ny
   sll_real64, dimension(:,:), allocatable            :: rot
   sll_real64, dimension(:), allocatable              :: kx
   sll_real64, dimension(:), allocatable              :: ky
   type(C_PTR)                                        :: fwx, fwy
   type(C_PTR)                                        :: bwx, bwy
   complex(C_DOUBLE_COMPLEX), dimension(:),   pointer :: hxt_x, hxt_y
   complex(C_DOUBLE_COMPLEX), dimension(:),   pointer :: hyt_x, hyt_y
   complex(C_DOUBLE_COMPLEX), dimension(:),   pointer :: hzt_x, hzt_y
   complex(C_DOUBLE_COMPLEX), dimension(:),   pointer :: ext_x, ext_y
   complex(C_DOUBLE_COMPLEX), dimension(:),   pointer :: eyt_x, eyt_y
   complex(C_DOUBLE_COMPLEX), dimension(:),   pointer :: ezt_x, ezt_y
   integer(C_SIZE_T)                                  :: sz_hxt_x, sz_hxt_y
   integer(C_SIZE_T)                                  :: sz_hyt_x, sz_hyt_y
   integer(C_SIZE_T)                                  :: sz_ezt_x, sz_ezt_y
   integer(C_SIZE_T)                                  :: sz_ext_x, sz_ext_y
   integer(C_SIZE_T)                                  :: sz_eyt_x, sz_eyt_y
   integer(C_SIZE_T)                                  :: sz_hzt_x, sz_hzt_y
   type(C_PTR)                                        :: p_hxt_y
   type(C_PTR)                                        :: p_hyt_x
   type(C_PTR)                                        :: p_hzt
   type(C_PTR)                                        :: p_ext
   type(C_PTR)                                        :: p_eyt
   type(C_PTR)                                        :: p_ezt_x, p_ezt_y
end type maxwell_pstd

sll_int32, private :: i, j

include 'fftw3.f03'

contains

subroutine new_maxwell_2d_pstd(self, xmin, xmax, nx, &
                                     ymin, ymax, ny, error )
   type(maxwell_pstd)                        :: self
   sll_real64                                :: xmin
   sll_real64                                :: xmax
   sll_real64                                :: ymin
   sll_real64                                :: ymax
   sll_int32                                 :: nx
   sll_int32                                 :: ny
   sll_int32                                 :: error
   sll_real64                                :: dx
   sll_real64                                :: dy
   sll_real64                                :: kx0
   sll_real64                                :: ky0

   self%nx = nx
   self%ny = ny

   FFTW_ALLOCATE(self%hxt_y,ny/2+1,self%sz_hxt_y,self%p_hxt_y)
   FFTW_ALLOCATE(self%hyt_x,nx/2+1,self%sz_hyt_x,self%p_hyt_x)
   FFTW_ALLOCATE(self%ezt_x,nx/2+1,self%sz_ezt_x,self%p_ezt_x)
   FFTW_ALLOCATE(self%ezt_y,ny/2+1,self%sz_ezt_y,self%p_ezt_y)

   SLL_ALLOCATE(self%rot(nx,ny), error)

   !call dfftw_init_threads(error)
   !if (error == 0) stop 'FFTW CAN''T USE THREADS'
   !call dfftw_plan_with_nthreads(nthreads)
   
   self%fwx = fftw_plan_dft_r2c_1d(nx, self%rot(:,1), self%ezt_x, FFTW_MEASURE)
   self%bwx = fftw_plan_dft_c2r_1d(nx, self%ezt_x, self%rot(:,1), FFTW_MEASURE)
   self%fwy = fftw_plan_dft_r2c_1d(ny, self%rot(1,:), self%ezt_y, FFTW_MEASURE)
   self%bwy = fftw_plan_dft_c2r_1d(ny, self%ezt_y, self%rot(1,:), FFTW_MEASURE)

   SLL_ALLOCATE(self%kx(nx/2+1), error)
   SLL_ALLOCATE(self%ky(ny/2+1), error)
   
   dx = (xmax-xmin) / nx
   dy = (ymax-ymin) / ny

   kx0 = 2._f64*sll_pi/(nx*dx)
   ky0 = 2._f64*sll_pi/(ny*dy)

   do i=1,nx/2+1
      self%kx(i) = (i-1)*kx0
   end do
   self%kx(1) = 1.0_f64
   do j=1,ny/2+1
      self%ky(j) = (j-1)*ky0
   end do
   self%ky(1) = 1.0_f64

end subroutine new_maxwell_2d_pstd

subroutine solve_maxwell_2d_pstd(self, hx, hy, ez, dt)

   type(maxwell_pstd)          :: self
   sll_real64 , intent(inout), dimension(:,:)   :: hx, hy, ez
   sll_real64 , intent(in)   :: dt

   !H(n-1/2)--> H(n+1/2) sur les pts interieurs   
   call faraday_pstd(self, hx, hy, ez, dt)   

   !call cl_periodiques(self, hx, hy, ez, dt)

   !E(n)-->E(n+1) sur les pts interieurs
   call ampere_maxwell_pstd(self, hx, hy, ez, dt) 

end subroutine solve_maxwell_2d_pstd

!> Solve faraday 
subroutine faraday_pstd(self, hx, hy, ez, dt)

   type(maxwell_pstd),intent(inout)  :: self
   sll_real64, dimension(:,:), intent(inout) :: hx
   sll_real64, dimension(:,:), intent(inout) :: hy
   sll_real64, dimension(:,:), intent(inout) :: ez
   sll_int32                                 :: nx
   sll_int32                                 :: ny
   sll_real64, intent(in)                    :: dt

   nx = self%nx
   ny = self%ny

   do i = 1, nx
      call fftw_execute_dft_r2c(self%fwy, ez(i,1:ny), self%ezt_y)
      self%ezt_y = -cmplx(0.0_f64,self%ky,kind=f64)*self%ezt_y
      call fftw_execute_dft_c2r(self%bwy, self%ezt_y, self%rot(i,:))
      hx(i,1:ny) = hx(i,1:ny) - dt * self%rot(i,:) / ny
   end do

   hx(nx+1,:) = hx(1,:) 
   hx(:,ny+1) = hx(:,1) 

   do j = 1, ny
      call fftw_execute_dft_r2c(self%fwx, ez(1:nx,j), self%ezt_x)
      self%ezt_x = -cmplx(0.0_f64,self%kx,kind=f64)*self%ezt_x
      call fftw_execute_dft_c2r(self%bwx, self%ezt_x, self%rot(:,j))
      hy(1:nx,j) = hy(1:nx,j) + dt * self%rot(:,j) / nx
   end do

   hy(nx+1,:) = hy(1,:)
   hy(:,ny+1) = hy(:,1)

end subroutine faraday_pstd

!> Solve ampere
subroutine ampere_maxwell_pstd(self, hx, hy, ez, dt)

   type(maxwell_pstd),intent(inout)  :: self
   sll_int32                    :: nx
   sll_int32                    :: ny
   sll_real64, dimension(:,:)   :: hx
   sll_real64, dimension(:,:)   :: hy
   sll_real64, dimension(:,:)   :: ez
   sll_real64                   :: dt

   nx = self%nx
   ny = self%ny

   do j = 1, ny
      call fftw_execute_dft_r2c(self%fwx, hy(1:nx,j), self%hyt_x)
      self%hyt_x = -cmplx(0.0_f64,self%kx,kind=f64)*self%hyt_x
      call fftw_execute_dft_c2r(self%bwx, self%hyt_x, self%rot(:,j))
      ez(1:nx,j) = ez(1:nx,j) + dt * self%rot(:,j) / nx
   end do

   ez(nx+1,:)   = ez(1,:)

   do i = 1, nx
      call fftw_execute_dft_r2c(self%fwy, hx(i,1:ny), self%hxt_y)
      self%hxt_y = -cmplx(0.0_f64,self%ky,kind=f64)*self%hxt_y
      call fftw_execute_dft_c2r(self%bwy, self%hxt_y, self%rot(i,:))
      ez(i,1:ny) = ez(i,1:ny) - dt * self%rot(i,:) / ny
   end do

   ez(:,ny+1)  = ez(:,1)

end subroutine ampere_maxwell_pstd


subroutine free_maxwell_2d_pstd(self)
type(maxwell_pstd) :: self
!sll_int32       :: error

if (c_associated(self%p_hxt_y)) call fftw_free(self%p_hxt_y)
if (c_associated(self%p_hyt_x)) call fftw_free(self%p_hyt_x)
if (c_associated(self%p_ezt_x)) call fftw_free(self%p_ezt_x)
if (c_associated(self%p_ezt_y)) call fftw_free(self%p_ezt_y)

call dfftw_destroy_plan(self%fwx)
call dfftw_destroy_plan(self%fwy)
call dfftw_destroy_plan(self%bwx)
call dfftw_destroy_plan(self%bwy)
!if (nthreads > 1) then
!   call dfftw_cleanup_threads(error)
!end if

end subroutine free_maxwell_2d_pstd

end module sll_maxwell_2d_pstd
