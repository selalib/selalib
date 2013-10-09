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

#include "sll_fftw.h"

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

#define D_DZ(field)                                           \
call fftw_execute_dft_r2c(self%fwz, field, self%tmp_z);       \
self%tmp_z = -cmplx(0.0_f64,self%kz,kind=f64)*self%tmp_z;     \
call fftw_execute_dft_c2r(self%bwz, self%tmp_z, self%d_dz);   \
self%d_dz = self%d_dz / nz

!> @author
!> Pierre Navaro
!> @brief
!> Implements the Maxwell solver in 3D with periodic boundary conditions
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
module sll_maxwell_3d_pstd
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_maxwell_solvers_macros.h"
#include "sll_constants.h"
use fftw3

implicit none
private

!> Initialize maxwell solver 2d cartesian periodic with PSTD scheme
interface new
 module procedure new_maxwell_3d_pstd
end interface new

!> Solve maxwell solver 2d cartesian periodic with PSTD scheme
interface solve
 module procedure solve_maxwell_3d
end interface solve

!> Delete maxwell solver 2d cartesian periodic with PSTD scheme
interface free
 module procedure free_maxwell_3d_pstd
end interface free

public :: new, free, solve, ampere, faraday

!> Maxwell solver object
type, public :: maxwell_pstd_3d
   sll_int32                          :: nx           !< x nodes number
   sll_int32                          :: ny           !< y nodes number
   sll_int32                          :: nz           !< y nodes number
   sll_real64, dimension(:), pointer  :: d_dx         !< field x derivative
   sll_real64, dimension(:), pointer  :: d_dy         !< field y derivative
   sll_real64, dimension(:), pointer  :: d_dz         !< field y derivative
   sll_real64, dimension(:), pointer  :: kx           !< x wave number
   sll_real64, dimension(:), pointer  :: ky           !< y wave number
   sll_real64, dimension(:), pointer  :: kz           !< y wave number
   fftw_plan                          :: fwx          !< forward fft plan along x
   fftw_plan                          :: fwy          !< forward fft plan along y
   fftw_plan                          :: fwz          !< forward fft plan along y
   fftw_plan                          :: bwx          !< backward fft plan along x
   fftw_plan                          :: bwy          !< backward fft plan along y
   fftw_plan                          :: bwz          !< backward fft plan along y
   fftw_plan                          :: p_tmp_x      !< pointer for memory allocation
   fftw_plan                          :: p_tmp_y      !< pointer for memory allocation
   fftw_plan                          :: p_tmp_z      !< pointer for memory allocation
   fftw_comp                , pointer :: tmp_x(:)     !< x fft transform
   fftw_comp                , pointer :: tmp_y(:)     !< y fft transform
   fftw_comp                , pointer :: tmp_z(:)     !< y fft transform
   fftw_int                           :: sz_tmp_x     !< size for memory allocation
   fftw_int                           :: sz_tmp_y     !< size for memory allocation
   fftw_int                           :: sz_tmp_z     !< size for memory allocation
   sll_int32                          :: polarization !< TE or TM
   sll_real64                         :: e_0          !< electric conductivity
   sll_real64                         :: mu_0         !< magnetic permeability
end type maxwell_pstd_3d

sll_int32, private :: i, j, k

contains

!> Initialize 2d maxwell solver on cartesian mesh with PSTD scheme
subroutine new_maxwell_3d_pstd(self,xmin,xmax,nx, &
                                    ymin,ymax,ny, &
                                    zmin,zmax,nz )

   type(maxwell_pstd_3d) :: self         !< maxwell object
   sll_real64            :: xmin         !< xmin
   sll_real64            :: xmax         !< xmax
   sll_real64            :: ymin         !< ymin
   sll_real64            :: ymax         !< ymax
   sll_real64            :: zmin         !< zmin
   sll_real64            :: zmax         !< zmax
   sll_int32             :: nx           !< x nodes number
   sll_int32             :: ny           !< y nodes number
   sll_int32             :: nz           !< z nodes number
   sll_int32             :: error        !< error code
   sll_real64            :: dx           !< x space step
   sll_real64            :: dy           !< y space step
   sll_real64            :: dz           !< z space step
   sll_real64            :: kx0
   sll_real64            :: ky0
   sll_real64            :: kz0

   fftw_int              :: sz_tmp_x
   fftw_int              :: sz_tmp_y
   fftw_int              :: sz_tmp_z

   self%nx = nx
   self%ny = ny
   self%nz = nz

   self%e_0  = 1._f64
   self%mu_0 = 1._f64

   FFTW_ALLOCATE(self%tmp_x,nx/2+1,sz_tmp_x,self%p_tmp_x)
   FFTW_ALLOCATE(self%tmp_y,ny/2+1,sz_tmp_y,self%p_tmp_y)
   FFTW_ALLOCATE(self%tmp_z,nz/2+1,sz_tmp_z,self%p_tmp_z)

   SLL_ALLOCATE(self%d_dx(nx), error)
   SLL_ALLOCATE(self%d_dy(ny), error)
   SLL_ALLOCATE(self%d_dz(nz), error)

   !call dfftw_init_threads(error)
   !if (error == 0) stop 'FFTW CAN''T USE THREADS'
   !call dfftw_plan_with_nthreads(nthreads)
   
   NEW_FFTW_PLAN_R2C_1D(self%fwx, nx, self%d_dx,  self%tmp_x)
   NEW_FFTW_PLAN_C2R_1D(self%bwx, nx, self%tmp_x, self%d_dx)
   NEW_FFTW_PLAN_R2C_1D(self%fwy, ny, self%d_dy,  self%tmp_y)
   NEW_FFTW_PLAN_C2R_1D(self%bwy, ny, self%tmp_y, self%d_dy)
   NEW_FFTW_PLAN_R2C_1D(self%fwz, nz, self%d_dz,  self%tmp_z)
   NEW_FFTW_PLAN_C2R_1D(self%bwz, nz, self%tmp_z, self%d_dz)

   SLL_ALLOCATE(self%kx(nx/2+1), error)
   SLL_ALLOCATE(self%ky(ny/2+1), error)
   SLL_ALLOCATE(self%kz(nz/2+1), error)
   
   dx = (xmax-xmin) / nx
   dy = (ymax-ymin) / ny
   dz = (zmax-zmin) / nz

   kx0 = 2._f64*sll_pi/(nx*dx)
   ky0 = 2._f64*sll_pi/(ny*dy)
   kz0 = 2._f64*sll_pi/(nz*dz)

   do i=2,nx/2+1
      self%kx(i) = (i-1)*kx0
   end do
   self%kx(1) = 1.0_f64
   do j=2,ny/2+1
      self%ky(j) = (j-1)*ky0
   end do
   self%ky(1) = 1.0_f64
   do k=2,nz/2+1
      self%kz(k) = (k-1)*kz0
   end do
   self%kz(1) = 1.0_f64

end subroutine new_maxwell_3d_pstd

!> this routine exists only for testing purpose. Use ampere and faraday
!> in your appication.
subroutine solve_maxwell_3d(self, ex, ey, ez, bx, by, bz, dt)

   type(maxwell_pstd_3d), intent(inout)           :: self !< maxwell object
   sll_real64, intent(inout), dimension(:,:,:) :: ex   !< Ex
   sll_real64, intent(inout), dimension(:,:,:) :: ey   !< Ey
   sll_real64, intent(inout), dimension(:,:,:) :: ez   !< Ez
   sll_real64, intent(inout), dimension(:,:,:) :: bx   !< Bx
   sll_real64, intent(inout), dimension(:,:,:) :: by   !< By
   sll_real64, intent(inout), dimension(:,:,:) :: bz   !< Bz
   sll_real64, intent(in)                      :: dt   !< time step

   call faraday(self, ex, ey, ez, bx, by, bz, 0.5*dt)   
   call ampere(self, ex, ey, ez, bx, by, bz, dt) 
   call faraday(self, ex, ey, ez, bx, by, bz, 0.5*dt)   

end subroutine solve_maxwell_3d


!!> Impose periodic boundary conditions
!subroutine bc_periodic(self, ex, ey, ez, hx, hy, hz)
!
!   type(maxwell_pstd_3d), intent(inout)          :: self !< maxwell object
!   sll_real64 , intent(inout), dimension(:,:) :: ex   !< Ex
!   sll_real64 , intent(inout), dimension(:,:) :: ey   !< Ey
!   sll_real64 , intent(inout), dimension(:,:) :: ez   !< Ez
!   sll_real64 , intent(inout), dimension(:,:) :: hx   !< Bx
!   sll_real64 , intent(inout), dimension(:,:) :: hy   !< By
!   sll_real64 , intent(inout), dimension(:,:) :: hz   !< Bz
!   sll_int32 :: nx, ny, nz
!
!   nx = self%nx
!   ny = self%ny
!   nz = self%nz
!
!   !fx(:,ny+1) = fx(:,1) 
!   !fy(nx+1,:) = fy(1,:)
!   !fz(nx+1,:) = fz(1,:)
!   !fz(:,ny+1) = fz(:,1)
!
!end subroutine bc_periodic


!> Solve faraday equation  (hx,hy,ez)
subroutine faraday(self, hx, hy, hz, ex, ey, ez, dt)

   type(maxwell_pstd_3d),intent(inout)         :: self  !< Maxwell object
   sll_real64, dimension(:,:,:), intent(inout) :: hx    !< Magnetic field x
   sll_real64, dimension(:,:,:), intent(inout) :: hy    !< Magnetic field y
   sll_real64, dimension(:,:,:), intent(inout) :: hz    !< Magnetic field z
   sll_real64, dimension(:,:,:), intent(inout) :: ex    !< Electric field x
   sll_real64, dimension(:,:,:), intent(inout) :: ey    !< Electric field y
   sll_real64, dimension(:,:,:), intent(inout) :: ez    !< Electric field z
   sll_int32                                   :: nx    !< x nodes number
   sll_int32                                   :: ny    !< y nodes number
   sll_int32                                   :: nz    !< z nodes number
   sll_real64, intent(in)                      :: dt    !< time step
   sll_real64                                  :: dt_mu

   nx = self%nx
   ny = self%ny
   nz = self%nz

   dt_mu = dt / self%mu_0 

   do k = 1, nz
      do i = 1, nx
         D_DY(ez(i,:,k))
         hx(i,:,k) = hx(i,:,k) - dt_mu * self%d_dy
         D_DY(ex(i,:,k))
         hz(i,:,k) = hz(i,:,k) + dt_mu * self%d_dy
      end do
   end do

   do j = 1, ny
      do i = 1, nx
         D_DZ(ey(i,j,:))
         hx(i,j,:) = hx(i,j,:) + dt_mu * self%d_dz
         D_DZ(ex(i,j,:))
         hy(i,j,:) = hy(i,j,:) - dt_mu * self%d_dz
      end do
   end do

   do k = 1, nz
      do j = 1, ny
         D_DX(ez(:,j,k))
         hy(:,j,k) = hy(:,j,k) + dt_mu * self%d_dx
         D_DX(ey(:,j,k))
         hz(:,j,k) = hz(:,j,k) - dt_mu * self%d_dx
      end do
   end do
   
end subroutine faraday

!> Solve ampere maxwell equation (hx,hy,ez)
subroutine ampere(self, hx, hy, hz, ex, ey, ez, dt, jx, jy, jz)

   type(maxwell_pstd_3d),intent(inout)    :: self !< maxwell object
   sll_real64, dimension(:,:,:)           :: hx   !< magnetic field x
   sll_real64, dimension(:,:,:)           :: hy   !< magnetic field y
   sll_real64, dimension(:,:,:)           :: hz   !< magnetic field z
   sll_real64, dimension(:,:,:)           :: ex   !< electric field x
   sll_real64, dimension(:,:,:)           :: ey   !< electric field y
   sll_real64, dimension(:,:,:)           :: ez   !< electric field z
   sll_real64                             :: dt   !< time step
   sll_real64, dimension(:,:,:), optional :: jx   !< current z
   sll_real64, dimension(:,:,:), optional :: jy   !< current z
   sll_real64, dimension(:,:,:), optional :: jz   !< current z

   sll_int32                              :: nx   !< x nodes number
   sll_int32                              :: ny   !< y nodes number
   sll_int32                              :: nz   !< z nodes number
   sll_real64                             :: dt_e

   nx = self%nx
   ny = self%ny
   nz = self%nz

   dt_e = dt / self%e_0

   do k = 1, nz
   do i = 1, nx
      D_DY(hz(i,:,k))
      ex(i,:,k) = ex(i,:,k) + dt_e * self%d_dy
      D_DY(hx(i,:,k))
      ez(i,:,k) = ez(i,:,k) - dt_e * self%d_dy
   end do
   end do

   do j = 1, ny
   do i = 1, nx
      D_DZ(hy(i,j,:))
      ex(i,j,:) = ex(i,j,:) - dt_e * self%d_dz
      D_DZ(hx(i,j,:))
      ey(i,j,:) = ey(i,j,:) + dt_e * self%d_dz
   end do
   end do

   do k = 1, nz
   do j = 1, ny
      D_DX(hz(:,j,k))
      ey(:,j,k) = ey(:,j,k) - dt_e * self%d_dx
      D_DX(hy(:,j,k))
      ez(:,j,k) = ez(:,j,k) + dt_e * self%d_dx
   end do
   end do

   if (present(jx) .and. present(jy) .and. present(jz)) then
      ex = ex - dt_e * jx 
      ey = ey - dt_e * jy 
      ez = ez - dt_e * jz 
   end if

end subroutine ampere

!> delete maxwell solver object
subroutine free_maxwell_3d_pstd(self)
type(maxwell_pstd_3d) :: self

#ifdef FFTW_F2003
if (c_associated(self%p_tmp_x)) call fftw_free(self%p_tmp_x)
if (c_associated(self%p_tmp_y)) call fftw_free(self%p_tmp_y)
if (c_associated(self%p_tmp_z)) call fftw_free(self%p_tmp_z)
#endif

call fftw_destroy_plan(self%fwx)
call fftw_destroy_plan(self%fwy)
call fftw_destroy_plan(self%fwz)
call fftw_destroy_plan(self%bwx)
call fftw_destroy_plan(self%bwy)
call fftw_destroy_plan(self%bwz)

!if (nthreads > 1) then
!   call dfftw_cleanup_threads(error)
!end if

end subroutine free_maxwell_3d_pstd

end module sll_maxwell_3d_pstd
