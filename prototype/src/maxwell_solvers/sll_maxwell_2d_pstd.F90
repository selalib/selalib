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


#define FFTW_ALLOCATE(array,array_size,sz_array,p_array)  \
sz_array = int((array_size/2+1),C_SIZE_T);                \
p_array = fftw_alloc_complex(sz_array);                   \
call c_f_pointer(p_array, array, [array_size/2+1])        \

#define D_DX(field)                                           \
call fftw_execute_dft_r2c(self%fwx, field, self%tmp_x);       \
self%tmp_x = -cmplx(0.0_f64,self%kx,kind=f64)*self%tmp_x;     \
call fftw_execute_dft_c2r(self%bwx, self%tmp_x, self%d_dx);   \
self%d_dx = self%d_dx / nx

#define D_DY(field)                                           \
call fftw_execute_dft_r2c(self%fwy, field, self%tmp_y);       \
self%tmp_y = -cmplx(0.0_f64,self%ky,kind=f64)*self%tmp_y;     \
call fftw_execute_dft_c2r(self%bwy, self%tmp_y, self%d_dy);    \
self%d_dy = self%d_dy / ny


module sll_maxwell_2d_pstd

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

use, intrinsic :: iso_c_binding
use numeric_constants

implicit none

interface initialize
 module procedure new_maxwell_2d_pstd
end interface initialize
interface solve
 module procedure solve_maxwell_2d
end interface solve
interface free
 module procedure free_maxwell_2d_pstd
end interface free

enum, bind(C)
   enumerator :: TE_POLARIZATION = 0, TM_POLARIZATION = 1
end enum

type, public :: maxwell_pstd
   sll_int32                                          :: nx
   sll_int32                                          :: ny
   sll_real64, dimension(:), allocatable              :: d_dx
   sll_real64, dimension(:), allocatable              :: d_dy
   sll_real64, dimension(:), allocatable              :: kx
   sll_real64, dimension(:), allocatable              :: ky
   type(C_PTR)                                        :: fwx, fwy
   type(C_PTR)                                        :: bwx, bwy
   complex(C_DOUBLE_COMPLEX), dimension(:),   pointer :: tmp_x, tmp_y
   integer(C_SIZE_T)                                  :: sz_tmp_x, sz_tmp_y
   type(C_PTR)                                        :: p_tmp_x, p_tmp_y
   sll_int32                                          :: polarization
end type maxwell_pstd

sll_int32, private :: i, j

include 'fftw3.f03'

contains



subroutine new_maxwell_2d_pstd(self, xmin, xmax, nx, &
                                     ymin, ymax, ny, &
                                     polarization, error )
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
   sll_int32                                 :: polarization

   self%nx = nx
   self%ny = ny
   self%polarization = polarization

   FFTW_ALLOCATE(self%tmp_x,nx/2+1,self%sz_tmp_x,self%p_tmp_x)
   FFTW_ALLOCATE(self%tmp_y,ny/2+1,self%sz_tmp_y,self%p_tmp_y)
   SLL_ALLOCATE(self%d_dx(nx), error)
   SLL_ALLOCATE(self%d_dy(ny), error)

   !call dfftw_init_threads(error)
   !if (error == 0) stop 'FFTW CAN''T USE THREADS'
   !call dfftw_plan_with_nthreads(nthreads)
   
   self%fwx = fftw_plan_dft_r2c_1d(nx, self%d_dx,  self%tmp_x, FFTW_ESTIMATE)
   self%bwx = fftw_plan_dft_c2r_1d(nx, self%tmp_x, self%d_dx,  FFTW_ESTIMATE)
   self%fwy = fftw_plan_dft_r2c_1d(ny, self%d_dy,  self%tmp_y, FFTW_ESTIMATE)
   self%bwy = fftw_plan_dft_c2r_1d(ny, self%tmp_y, self%d_dy,  FFTW_ESTIMATE)

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

subroutine solve_maxwell_2d(self, fx, fy, fz, dt)

   type(maxwell_pstd)          :: self
   sll_real64 , intent(inout), dimension(:,:)   :: fx
   sll_real64 , intent(inout), dimension(:,:)   :: fy
   sll_real64 , intent(inout), dimension(:,:)   :: fz
   sll_real64 , intent(in)   :: dt

   IF ( self%polarization == TM_POLARIZATION) then
      call faraday_tm(self, fx, fy, fz, dt)   
      call bc_periodic(self, fx, fy, fz)
      call ampere_maxwell_tm(self, fx, fy, fz, dt) 
   end if

   IF ( self%polarization == TE_POLARIZATION) then
      call faraday_te(self, fx, fy, fz, dt)   
      call bc_periodic(self, fx, fy, fz)
      call ampere_maxwell_te(self, fx, fy, fz, dt) 
   end if

end subroutine solve_maxwell_2d


subroutine bc_periodic(self, fx, fy, fz)
   type(maxwell_pstd)          :: self
   sll_real64 , intent(inout), dimension(:,:) :: fx
   sll_real64 , intent(inout), dimension(:,:) :: fy
   sll_real64 , intent(inout), dimension(:,:) :: fz
   sll_int32 :: nx, ny

   nx = self%nx
   ny = self%ny

   fx(:,ny+1) = fx(:,1) 
   fy(nx+1,:) = fy(1,:)
   fz(nx+1,:) = fz(1,:)
   fz(:,ny+1) = fz(:,1)

end subroutine bc_periodic


!> Solve faraday 
subroutine faraday_tm(self, hx, hy, ez, dt)

   type(maxwell_pstd),intent(inout)  :: self
   sll_real64, dimension(:,:), intent(inout) :: hx
   sll_real64, dimension(:,:), intent(inout) :: hy
   sll_real64, dimension(:,:), intent(inout) :: ez
   sll_int32                                 :: nx
   sll_int32                                 :: ny
   sll_real64, intent(in)                    :: dt

   nx = self%nx
   ny = self%ny

   do i = 1, nx+1
      D_DY(ez(i,1:ny))
      hx(i,1:ny) = hx(i,1:ny) - dt * self%d_dy
   end do


   do j = 1, ny+1
      D_DX(ez(1:nx,j))
      hy(1:nx,j) = hy(1:nx,j) + dt * self%d_dx
   end do


end subroutine faraday_tm

!> Solve faraday 
subroutine faraday_te(self, ex, ey, hz, dt)

   type(maxwell_pstd),intent(inout)  :: self
   sll_real64, dimension(:,:), intent(inout) :: ex
   sll_real64, dimension(:,:), intent(inout) :: ey
   sll_real64, dimension(:,:), intent(inout) :: hz
   sll_int32                                 :: nx
   sll_int32                                 :: ny
   sll_real64, intent(in)                    :: dt

   nx = self%nx
   ny = self%ny

   do i = 1, nx+1
      D_DY(ex(i,1:ny))
      hz(i,1:ny) = hz(i,1:ny) + dt * self%d_dy
   end do

   do j = 1, ny+1
      D_DX(ey(1:nx,j))
      hz(1:nx,j) = hz(1:nx,j) - dt * self%d_dx
   end do

end subroutine faraday_te

!> Solve ampere
subroutine ampere_maxwell_tm(self, hx, hy, ez, dt)

   type(maxwell_pstd),intent(inout)  :: self
   sll_int32                    :: nx
   sll_int32                    :: ny
   sll_real64, dimension(:,:)   :: hx
   sll_real64, dimension(:,:)   :: hy
   sll_real64, dimension(:,:)   :: ez
   sll_real64                   :: dt

   nx = self%nx
   ny = self%ny

   do j = 1, ny+1
      D_DX(hy(1:nx,j))
      ez(1:nx,j) = ez(1:nx,j) + dt * self%d_dx
   end do


   do i = 1, nx+1
      D_DY(hx(i,1:ny))
      ez(i,1:ny) = ez(i,1:ny) - dt * self%d_dy
   end do


end subroutine ampere_maxwell_tm

!> Solve ampere
subroutine ampere_maxwell_te(self, ex, ey, hz, dt)

   type(maxwell_pstd),intent(inout)  :: self
   sll_int32                    :: nx
   sll_int32                    :: ny
   sll_real64, dimension(:,:)   :: ex
   sll_real64, dimension(:,:)   :: ey
   sll_real64, dimension(:,:)   :: hz
   sll_real64                   :: dt

   nx = self%nx
   ny = self%ny

   do j = 1, ny+1
      D_DX(hz(1:nx,j))
      ey(1:nx,j) = ey(1:nx,j) - dt * self%d_dx
   end do


   do i = 1, nx+1
      D_DY(hz(i,1:ny))
      ex(i,1:ny) = ex(i,1:ny) + dt * self%d_dy
   end do


end subroutine ampere_maxwell_te

subroutine free_maxwell_2d_pstd(self)
type(maxwell_pstd) :: self
!sll_int32       :: error

if (c_associated(self%p_tmp_x)) call fftw_free(self%p_tmp_x)
if (c_associated(self%p_tmp_y)) call fftw_free(self%p_tmp_y)

call dfftw_destroy_plan(self%fwx)
call dfftw_destroy_plan(self%fwy)
call dfftw_destroy_plan(self%bwx)
call dfftw_destroy_plan(self%bwy)
!if (nthreads > 1) then
!   call dfftw_cleanup_threads(error)
!end if

end subroutine free_maxwell_2d_pstd

end module sll_maxwell_2d_pstd
