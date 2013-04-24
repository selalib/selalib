!**************************************************************
!  Copyright INRIA
!  Authors : 
!     Pierre Navaro 
!  
!  This code SeLaLib (for Semi-Lagrangian-Library) 
!  is a parallel library for simulating the plasma turbulence 
!  in a tokamak.
!  
!  This software is governed by the CeCILL-B license 
!  under French law and abiding by the rules of distribution 
!  of free software.  You can  use, modify and redistribute 
!  the software under the terms of the CeCILL-B license as 
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info". 
!**************************************************************


#define FFTW_ALLOCATE(array,array_size,sz_array,p_array)      \
sz_array = int((array_size/2+1),C_SIZE_T);                    \
p_array = fftw_alloc_complex(sz_array);                       \
call c_f_pointer(p_array, array, [array_size/2+1])            \

#define D_DX(field)                                           \
call fftw_execute_dft_r2c(self%fwx, field, self%tmp_x);       \
self%tmp_x = -cmplx(0.0_f64,self%kx,kind=f64)*self%tmp_x;     \
call fftw_execute_dft_c2r(self%bwx, self%tmp_x, self%d_dx);   \
self%d_dx = self%d_dx / nx

#define D_DY(field)                                           \
call fftw_execute_dft_r2c(self%fwy, field, self%tmp_y);       \
self%tmp_y = -cmplx(0.0_f64,self%ky,kind=f64)*self%tmp_y;     \
call fftw_execute_dft_c2r(self%bwy, self%tmp_y, self%d_dy);   \
self%d_dy = self%d_dy / ny


!> @author
!> Pierre Navaro
!> @brief
!> Implements the Maxwell solver in 2D with periodic boundary conditions
!> with PSTD method.
!> @details
!> Field derivative is made using Fast Fourier Transform.
!>
!>\f$ \displaystyle
!>\frac{\partial \psi}{\partial x} \Big|_i = F_x^{-1} [ -jk_xF_x(\psi)]_i,
!>\f$
!>
!> For Maxwell system the scheme is
!>
!>\f$ \displaystyle
!>H_u\Big|^{n+1/2}_{i,j,k} = H_u\Big|^{n-1/2}_{i,j,k}  - \frac{\Delta t}{\mu} \Big\{F_v^{-1}[-jk_vF_v(E_w)]|_{i,j,k}^n -F_w^{-1}[-jk_wF_w(E_v)]|_{i,j,k}^{n}\Big\},
!>\f$
!>
!>\f$ \displaystyle
!>E_u\Big|^{n+1}_{i,j,k} = E_u\Big|^{n}_{i,j,k}  + \frac{\Delta t}{\varepsilon} \Big\{F_v^{-1}[-jk_vF_v(H_w)]|_{i,j,k}^{n+1/2} -F_w^{-1}[-jk_wF_w(H_v)]|_{i,j,k}^{n+1/2}\Big\},
!>\f$
!>
!>where \f$(u,v,w) = (x,y,z),(y,z,x),(z,x,y)\f$


module sll_maxwell_2d_pstd
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_maxwell_solvers_macros.h"


use, intrinsic :: iso_c_binding
use sll_constants

implicit none
private

!> Initialize maxwell solver 2d cartesian periodic with PSTD scheme
interface new
 module procedure new_maxwell_2d_pstd
end interface new

!> Solve maxwell solver 2d cartesian periodic with PSTD scheme
interface solve
 module procedure solve_maxwell_2d
end interface solve

!> Delete maxwell solver 2d cartesian periodic with PSTD scheme
interface free
 module procedure free_maxwell_2d_pstd
end interface free

public :: new, free, solve, ampere_te, faraday_te, ampere_tm, faraday_tm

!> Maxwell solver object
type, public :: maxwell_pstd
   sll_int32                          :: nx           !< x nodes number
   sll_int32                          :: ny           !< y nodes number
   sll_real64, dimension(:), pointer  :: d_dx         !< field x derivative
   sll_real64, dimension(:), pointer  :: d_dy         !< field y derivative
   sll_real64, dimension(:), pointer  :: kx           !< x wave number
   sll_real64, dimension(:), pointer  :: ky           !< y wave number
   type(C_PTR)                        :: fwx          !< forward fft plan along x
   type(C_PTR)                        :: fwy          !< forward fft plan along y
   type(C_PTR)                        :: bwx          !< backward fft plan along x
   type(C_PTR)                        :: bwy          !< backward fft plan along y
   type(C_PTR)                        :: p_tmp_x      !< pointer for memory allocation
   type(C_PTR)                        :: p_tmp_y      !< pointer for memory allocation
   complex(C_DOUBLE_COMPLEX), pointer :: tmp_x(:)     !< x fft transform
   complex(C_DOUBLE_COMPLEX), pointer :: tmp_y(:)     !< y fft transform
   integer(C_SIZE_T)                  :: sz_tmp_x     !< size for memory allocation
   integer(C_SIZE_T)                  :: sz_tmp_y     !< size for memory allocation
   sll_int32                          :: polarization !< TE or TM
   sll_real64                         :: e_0          !< electric conductivity
   sll_real64                         :: mu_0         !< magnetic permeability
end type maxwell_pstd

sll_int32, private :: i, j

include 'fftw3.f03'

contains

!> Initialize 2d maxwell solver on cartesian mesh with PSTD scheme
subroutine new_maxwell_2d_pstd(self,xmin,xmax,nx,ymin,ymax,ny,polarization)

   type(maxwell_pstd) :: self         !< maxwell object
   sll_real64         :: xmin         !< xmin
   sll_real64         :: xmax         !< xmax
   sll_real64         :: ymin         !< ymin
   sll_real64         :: ymax         !< ymax
   sll_int32          :: nx           !< x nodes number
   sll_int32          :: ny           !< y nodes number
   sll_int32          :: polarization !< TE or TM
   sll_int32          :: error        !< error code
   sll_real64         :: dx           !< x space step
   sll_real64         :: dy           !< y space step
   sll_real64         :: kx0
   sll_real64         :: ky0
   integer(C_SIZE_T)  :: sz_tmp_x
   integer(C_SIZE_T)  :: sz_tmp_y

   self%nx = nx
   self%ny = ny
   self%polarization = polarization

   self%e_0  = 1._f64
   self%mu_0 = 1._f64

   FFTW_ALLOCATE(self%tmp_x,nx,sz_tmp_x,self%p_tmp_x)
   FFTW_ALLOCATE(self%tmp_y,ny,sz_tmp_y,self%p_tmp_y)
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

   do i=2,nx/2+1
      self%kx(i) = (i-1)*kx0
   end do
   self%kx(1) = 1.0_f64
   do j=2,ny/2+1
      self%ky(j) = (j-1)*ky0
   end do
   self%ky(1) = 1.0_f64

end subroutine new_maxwell_2d_pstd

!> this routine exists only for testing purpose. Use ampere and faraday
!> in your appication.
subroutine solve_maxwell_2d(self, fx, fy, fz, dt)

   type(maxwell_pstd), intent(inout)          :: self !< maxwell object
   sll_real64 , intent(inout), dimension(:,:) :: fx   !< Ex or Bx
   sll_real64 , intent(inout), dimension(:,:) :: fy   !< Ey or By
   sll_real64 , intent(inout), dimension(:,:) :: fz   !< Bz or Ez
   sll_real64 , intent(in)                    :: dt   !< time step

   IF ( self%polarization == TM_POLARIZATION) then
      call faraday_tm(self, fx, fy, fz, 0.5*dt)   
      call bc_periodic(self, fx, fy, fz)
      call ampere_tm(self, fx, fy, fz, dt) 
      call bc_periodic(self, fx, fy, fz)
      call faraday_tm(self, fx, fy, fz, 0.5*dt)   
      call bc_periodic(self, fx, fy, fz)
   end if

   IF ( self%polarization == TE_POLARIZATION) then
      call faraday_te(self, fx, fy, fz, 0.5*dt)   
      call bc_periodic(self, fx, fy, fz)
      call ampere_te(self, fx, fy, fz, dt) 
      call bc_periodic(self, fx, fy, fz)
      call faraday_te(self, fx, fy, fz, 0.5*dt)   
      call bc_periodic(self, fx, fy, fz)
   end if

end subroutine solve_maxwell_2d


!> Impose periodic boundary conditions
subroutine bc_periodic(self, fx, fy, fz)

   type(maxwell_pstd), intent(inout)          :: self !< maxwell object
   sll_real64 , intent(inout), dimension(:,:) :: fx   !< Ex or Bx
   sll_real64 , intent(inout), dimension(:,:) :: fy   !< Ey or By
   sll_real64 , intent(inout), dimension(:,:) :: fz   !< Bz or Ez
   sll_int32 :: nx, ny

   nx = self%nx
   ny = self%ny

   fx(:,ny+1) = fx(:,1) 
   fy(nx+1,:) = fy(1,:)
   fz(nx+1,:) = fz(1,:)
   fz(:,ny+1) = fz(:,1)

end subroutine bc_periodic


!> Solve faraday equation  (hx,hy,ez)
subroutine faraday_tm(self, hx, hy, ez, dt)

   type(maxwell_pstd),intent(inout)          :: self    !< Maxwell object
   sll_real64, dimension(:,:), intent(inout) :: hx      !< Magnetic field x
   sll_real64, dimension(:,:), intent(inout) :: hy      !< Magnetic field y
   sll_real64, dimension(:,:), intent(inout) :: ez      !< Electric field z
   sll_int32                                 :: nx      !< x nodes number
   sll_int32                                 :: ny      !< y nodes number
   sll_real64, intent(in)                    :: dt      !< time step
   sll_real64                                :: dt_mu

   nx = self%nx
   ny = self%ny

   dt_mu = dt / self%mu_0 

   do i = 1, nx
      D_DY(ez(i,1:ny))
      hx(i,1:ny) = hx(i,1:ny) - dt_mu * self%d_dy
   end do

   do j = 1, ny
      D_DX(ez(1:nx,j))
      hy(1:nx,j) = hy(1:nx,j) + dt_mu * self%d_dx
   end do

end subroutine faraday_tm

!> Solve faraday equation (ex,ey,hz)
subroutine faraday_te(self, ex, ey, hz, dt)

   type(maxwell_pstd),intent(inout)          :: self   !< maxwell object
   sll_real64, dimension(:,:), intent(inout) :: ex     !< electric field x
   sll_real64, dimension(:,:), intent(inout) :: ey     !< electric field y
   sll_real64, dimension(:,:), intent(inout) :: hz     !< magnetic field z
   sll_int32                                 :: nx     !< x nodes number
   sll_int32                                 :: ny     !< y nodes number
   sll_real64, intent(in)                    :: dt     !< time step
   sll_real64                                :: dt_mu

   nx = self%nx
   ny = self%ny

   dt_mu = dt / self%mu_0 

   do i = 1, nx
      D_DY(ex(i,1:ny))
      hz(i,1:ny) = hz(i,1:ny) + dt_mu * self%d_dy
   end do

   do j = 1, ny
      D_DX(ey(1:nx,j))
      hz(1:nx,j) = hz(1:nx,j) - dt_mu * self%d_dx
   end do

end subroutine faraday_te

!> Solve ampere maxwell equation (hx,hy,ez)
subroutine ampere_tm(self, hx, hy, ez, dt, jz)

   type(maxwell_pstd),intent(inout)      :: self   !< maxwell object
   sll_int32                             :: nx     !< x nodes number
   sll_int32                             :: ny     !< y nodes number
   sll_real64, dimension(:,:)            :: hx     !< magnetic field x
   sll_real64, dimension(:,:)            :: hy     !< magnetic field y
   sll_real64, dimension(:,:)            :: ez     !< electric field z
   sll_real64                            :: dt     !< time step
   sll_real64                            :: dt_e
   sll_real64, dimension(:,:), optional  :: jz     !< current z

   nx = self%nx
   ny = self%ny

   dt_e = dt / self%e_0

   do j = 1, ny
      D_DX(hy(1:nx,j))
      ez(1:nx,j) = ez(1:nx,j) + dt_e * self%d_dx
   end do

   do i = 1, nx
      D_DY(hx(i,1:ny))
      ez(i,1:ny) = ez(i,1:ny) - dt_e * self%d_dy
   end do

   if (present(jz)) then
      ez = ez - dt_e * jz 
   end if

end subroutine ampere_tm

!> Solve ampere maxwell equation (ex,ey,hz)
subroutine ampere_te(self, ex, ey, hz, dt, jx, jy)

   type(maxwell_pstd),intent(inout)      :: self   !< maxwell equation
   sll_real64, dimension(:,:)            :: ex     !< electric field x
   sll_real64, dimension(:,:)            :: ey     !< electric field y
   sll_real64, dimension(:,:)            :: hz     !< magnetic field z
   sll_real64                            :: dt     !< time step
   sll_real64, dimension(:,:), optional  :: jx     !< current x
   sll_real64, dimension(:,:), optional  :: jy     !< current y

   sll_real64                            :: dt_e
   sll_int32                             :: nx
   sll_int32                             :: ny

   nx = self%nx
   ny = self%ny

   dt_e = dt / self%e_0

   do j = 1, ny
      D_DX(hz(1:nx,j))
      ey(1:nx,j) = ey(1:nx,j) - dt_e * self%d_dx
   end do

   do i = 1, nx
      D_DY(hz(i,1:ny))
      ex(i,1:ny) = ex(i,1:ny) + dt_e * self%d_dy
   end do

   if (present(jx) .and. present(jy)) then
      ex = ex - dt_e * jx 
      ey = ey - dt_e * jy 
   end if

end subroutine ampere_te

!> delete maxwell solver object
subroutine free_maxwell_2d_pstd(self)
type(maxwell_pstd) :: self

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
