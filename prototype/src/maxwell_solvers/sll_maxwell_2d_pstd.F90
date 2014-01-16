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
self%d_dx = self%d_dx / nc_x

#define D_DY(field)                                           \
call fftw_execute_dft_r2c(self%fwy, field, self%tmp_y);       \
self%tmp_y = -cmplx(0.0_f64,self%ky,kind=f64)*self%tmp_y;     \
call fftw_execute_dft_c2r(self%bwy, self%tmp_y, self%d_dy);   \
self%d_dy = self%d_dy / nc_y

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
#include "sll_constants.h"
#include "sll_maxwell_solvers_macros.h"

use fftw3

implicit none
private

!> Initialize maxwell solver 2d cartesian periodic with PSTD scheme
interface initialize
 module procedure new_maxwell_2d_pstd
end interface initialize

!> Solve maxwell solver 2d cartesian periodic with PSTD scheme
interface solve
 module procedure solve_maxwell_2d_pstd
end interface solve

interface ampere
 module procedure ampere_2d_pstd
end interface ampere

interface faraday
 module procedure faraday_2d_pstd
end interface faraday


!> Delete maxwell solver 2d cartesian periodic with PSTD scheme
interface delete
 module procedure free_maxwell_2d_pstd
end interface delete

public :: initialize, delete, solve, ampere, faraday

!> Maxwell solver object
type, public :: maxwell_2d_pstd
   sll_int32                          :: nc_x         !< x cells number
   sll_int32                          :: nc_y         !< y cells number
   sll_real64, dimension(:), pointer  :: d_dx         !< field x derivative
   sll_real64, dimension(:), pointer  :: d_dy         !< field y derivative
   sll_real64, dimension(:), pointer  :: kx           !< x wave number
   sll_real64, dimension(:), pointer  :: ky           !< y wave number
   fftw_plan                          :: fwx          !< forward fft plan along x
   fftw_plan                          :: fwy          !< forward fft plan along y
   fftw_plan                          :: bwx          !< backward fft plan along x
   fftw_plan                          :: bwy          !< backward fft plan along y
   fftw_plan                          :: p_tmp_x      !< pointer for memory allocation
   fftw_plan                          :: p_tmp_y      !< pointer for memory allocation
   fftw_comp                , pointer :: tmp_x(:)     !< x fft transform
   fftw_comp                , pointer :: tmp_y(:)     !< y fft transform
   fftw_int                           :: sz_tmp_x     !< size for memory allocation
   fftw_int                           :: sz_tmp_y     !< size for memory allocation
   sll_int32                          :: polarization !< TE or TM
   sll_real64                         :: e_0          !< electric conductivity
   sll_real64                         :: mu_0         !< magnetic permeability

end type maxwell_2d_pstd

sll_int32, private :: i, j

contains

!> Initialize 2d maxwell solver on cartesian mesh with PSTD scheme
subroutine new_maxwell_2d_pstd(self,xmin,xmax,nc_x,ymin,ymax,nc_y,polarization)

   type(maxwell_2d_pstd) :: self         !< maxwell object
   sll_real64         :: xmin         !< xmin
   sll_real64         :: xmax         !< xmax
   sll_real64         :: ymin         !< ymin
   sll_real64         :: ymax         !< ymax
   sll_int32          :: nc_x         !< x cells number
   sll_int32          :: nc_y         !< y cells number
   sll_int32          :: polarization !< TE or TM
   sll_int32          :: error        !< error code
   sll_real64         :: dx           !< x space step
   sll_real64         :: dy           !< y space step
   sll_real64         :: kx0
   sll_real64         :: ky0

   self%nc_x = nc_x
   self%nc_y = nc_y
   self%polarization = polarization

   self%e_0  = 1._f64
   self%mu_0 = 1._f64

   SLL_ALLOCATE(self%d_dx(nc_x), error)
   SLL_ALLOCATE(self%d_dy(nc_y), error)

   FFTW_ALLOCATE(self%tmp_x,nc_x/2+1,self%sz_tmp_x,self%p_tmp_x)
   FFTW_ALLOCATE(self%tmp_y,nc_y/2+1,self%sz_tmp_y,self%p_tmp_y)

   NEW_FFTW_PLAN_R2C_1D(self%fwx, nc_x, self%d_dx,  self%tmp_x)
   NEW_FFTW_PLAN_C2R_1D(self%bwx, nc_x, self%tmp_x, self%d_dx)
   NEW_FFTW_PLAN_R2C_1D(self%fwy, nc_y, self%d_dy,  self%tmp_y)
   NEW_FFTW_PLAN_C2R_1D(self%bwy, nc_y, self%tmp_y, self%d_dy)

   !call dfftw_init_threads(error)
   !if (error == 0) stop 'FFTW CAN''T USE THREADS'
   !call dfftw_plan_with_nthreads(nthreads)
   

   SLL_ALLOCATE(self%kx(nc_x/2+1), error)
   SLL_ALLOCATE(self%ky(nc_y/2+1), error)
   
   dx = (xmax-xmin) / nc_x
   dy = (ymax-ymin) / nc_y

   kx0 = 2._f64*sll_pi/(nc_x*dx)
   ky0 = 2._f64*sll_pi/(nc_y*dy)

   do i=2,nc_x/2+1
      self%kx(i) = (i-1)*kx0
   end do
   self%kx(1) = 1.0_f64
   do j=2,nc_y/2+1
      self%ky(j) = (j-1)*ky0
   end do
   self%ky(1) = 1.0_f64

end subroutine new_maxwell_2d_pstd

!> this routine exists only for testing purpose. Use ampere and faraday
!> in your appication.
subroutine solve_maxwell_2d_pstd(self, fx, fy, fz, dt)

   type(maxwell_2d_pstd), intent(inout)          :: self !< maxwell object
   sll_real64 , intent(inout), dimension(:,:) :: fx   !< Ex or Bx
   sll_real64 , intent(inout), dimension(:,:) :: fy   !< Ey or By
   sll_real64 , intent(inout), dimension(:,:) :: fz   !< Bz or Ez
   sll_real64 , intent(in)                    :: dt   !< time step

   IF ( self%polarization == TM_POLARIZATION) then
      call faraday_tm_2d_pstd(self, fx, fy, fz, 0.5*dt)   
      call bc_periodic_2d_pstd(self, fx, fy, fz)
      call ampere_tm_2d_pstd(self, fx, fy, fz, dt) 
      call bc_periodic_2d_pstd(self, fx, fy, fz)
      call faraday_tm_2d_pstd(self, fx, fy, fz, 0.5*dt)   
      call bc_periodic_2d_pstd(self, fx, fy, fz)
   end if

   IF ( self%polarization == TE_POLARIZATION) then
      call faraday_te_2d_pstd(self, fx, fy, fz, 0.5*dt)   
      call bc_periodic_2d_pstd(self, fx, fy, fz)
      call ampere_te_2d_pstd(self, fx, fy, fz, dt) 
      call bc_periodic_2d_pstd(self, fx, fy, fz)
      call faraday_te_2d_pstd(self, fx, fy, fz, 0.5*dt)   
      call bc_periodic_2d_pstd(self, fx, fy, fz)
   end if

end subroutine solve_maxwell_2d_pstd


!> Impose periodic boundary conditions
subroutine bc_periodic_2d_pstd(self, fx, fy, fz)

   type(maxwell_2d_pstd), intent(inout)          :: self !< maxwell object
   sll_real64 , intent(inout), dimension(:,:) :: fx   !< Ex or Bx
   sll_real64 , intent(inout), dimension(:,:) :: fy   !< Ey or By
   sll_real64 , intent(inout), dimension(:,:) :: fz   !< Bz or Ez
   sll_int32 :: nc_x, nc_y

   nc_x = self%nc_x
   nc_y = self%nc_y

   fx(:,nc_y+1) = fx(:,1) 
   fy(nc_x+1,:) = fy(1,:)
   fz(nc_x+1,:) = fz(1,:)
   fz(:,nc_y+1) = fz(:,1)

end subroutine bc_periodic_2d_pstd

!> Solve Faraday equation
subroutine faraday_2d_pstd(self, fx, fy, fz, dt)
   type(maxwell_2d_pstd),intent(inout)       :: self    !< Maxwell object
   sll_real64, dimension(:,:), intent(inout) :: fx      !< field x
   sll_real64, dimension(:,:), intent(inout) :: fy      !< field y
   sll_real64, dimension(:,:), intent(inout) :: fz      !< field z
   sll_real64 , intent(in)                   :: dt      !< time step

   if ( self%polarization == TM_POLARIZATION) then
      call faraday_tm_2d_pstd(self, fx, fy, fz, dt)
   end if

   if ( self%polarization == TE_POLARIZATION) then
      call faraday_te_2d_pstd(self, fx, fy, fz, dt)
   end if

end subroutine faraday_2d_pstd

subroutine ampere_2d_pstd(self, fx, fy, fz, dt, sx, sy)
   type(maxwell_2d_pstd),intent(inout)       :: self    !< Maxwell object
   sll_real64, dimension(:,:), intent(inout) :: fx      !< field x
   sll_real64, dimension(:,:), intent(inout) :: fy      !< field y
   sll_real64, dimension(:,:), intent(inout) :: fz      !< field z
   sll_real64 , intent(in)                   :: dt      !< time step
   sll_real64, dimension(:,:), optional      :: sx      !< source x
   sll_real64, dimension(:,:), optional      :: sy      !< source y

   if ( self%polarization == TM_POLARIZATION .and. &
        present(sx) .and. .not. present(sy)) then
      call ampere_tm_2d_pstd(self, fx, fy, fz, dt, sx)
   else if ( self%polarization == TE_POLARIZATION .and. &
        present(sx) .and. present(sy)) then
      call ampere_te_2d_pstd(self, fx, fy, fz, dt, sx, sy)
   end if

end subroutine ampere_2d_pstd

!> Solve faraday equation  (hx,hy,ez)
subroutine faraday_tm_2d_pstd(self, hx, hy, ez, dt)

   type(maxwell_2d_pstd),intent(inout)       :: self    !< Maxwell object
   sll_real64, dimension(:,:), intent(inout) :: hx      !< Magnetic field x
   sll_real64, dimension(:,:), intent(inout) :: hy      !< Magnetic field y
   sll_real64, dimension(:,:), intent(inout) :: ez      !< Electric field z
   sll_int32                                 :: nc_x    !< x cells number
   sll_int32                                 :: nc_y    !< y cells number
   sll_real64, intent(in)                    :: dt      !< time step
   sll_real64                                :: dt_mu

   nc_x = self%nc_x
   nc_y = self%nc_y

   dt_mu = dt / self%mu_0 

   do i = 1, nc_x
      D_DY(ez(i,1:nc_y))
      hx(i,1:nc_y) = hx(i,1:nc_y) - dt_mu * self%d_dy
   end do

   do j = 1, nc_y
      D_DX(ez(1:nc_x,j))
      hy(1:nc_x,j) = hy(1:nc_x,j) + dt_mu * self%d_dx
   end do

end subroutine faraday_tm_2d_pstd

!> Solve faraday equation (ex,ey,hz)
subroutine faraday_te_2d_pstd(self, ex, ey, hz, dt)

   type(maxwell_2d_pstd),intent(inout)       :: self   !< maxwell object
   sll_real64, dimension(:,:), intent(inout) :: ex     !< electric field x
   sll_real64, dimension(:,:), intent(inout) :: ey     !< electric field y
   sll_real64, dimension(:,:), intent(inout) :: hz     !< magnetic field z
   sll_int32                                 :: nc_x   !< x cells number
   sll_int32                                 :: nc_y   !< y cells number
   sll_real64, intent(in)                    :: dt     !< time step
   sll_real64                                :: dt_mu

   nc_x = self%nc_x
   nc_y = self%nc_y

   dt_mu = dt / self%mu_0 

   do i = 1, nc_x
      D_DY(ex(i,1:nc_y))
      hz(i,1:nc_y) = hz(i,1:nc_y) + dt_mu * self%d_dy
   end do

   do j = 1, nc_y
      D_DX(ey(1:nc_x,j))
      hz(1:nc_x,j) = hz(1:nc_x,j) - dt_mu * self%d_dx
   end do

end subroutine faraday_te_2d_pstd

!> Solve ampere maxwell equation (hx,hy,ez)
subroutine ampere_tm_2d_pstd(self, hx, hy, ez, dt, jz)

   type(maxwell_2d_pstd),intent(inout)   :: self   !< maxwell object
   sll_int32                             :: nc_x   !< x cells number
   sll_int32                             :: nc_y   !< y cells number
   sll_real64, dimension(:,:)            :: hx     !< magnetic field x
   sll_real64, dimension(:,:)            :: hy     !< magnetic field y
   sll_real64, dimension(:,:)            :: ez     !< electric field z
   sll_real64                            :: dt     !< time step
   sll_real64                            :: dt_e
   sll_real64, dimension(:,:), optional  :: jz     !< current z

   nc_x = self%nc_x
   nc_y = self%nc_y

   dt_e = dt / self%e_0

   do j = 1, nc_y
      D_DX(hy(1:nc_x,j))
      ez(1:nc_x,j) = ez(1:nc_x,j) + dt_e * self%d_dx
   end do

   do i = 1, nc_x
      D_DY(hx(i,1:nc_y))
      ez(i,1:nc_y) = ez(i,1:nc_y) - dt_e * self%d_dy
   end do

   if (present(jz)) then
      ez = ez - dt_e * jz 
   end if

end subroutine ampere_tm_2d_pstd

!> Solve ampere maxwell equation (ex,ey,hz)
subroutine ampere_te_2d_pstd(self, ex, ey, hz, dt, jx, jy)

   type(maxwell_2d_pstd),intent(inout)      :: self   !< maxwell equation
   sll_real64, dimension(:,:)            :: ex     !< electric field x
   sll_real64, dimension(:,:)            :: ey     !< electric field y
   sll_real64, dimension(:,:)            :: hz     !< magnetic field z
   sll_real64                            :: dt     !< time step
   sll_real64, dimension(:,:), optional  :: jx     !< current x
   sll_real64, dimension(:,:), optional  :: jy     !< current y

   sll_real64                            :: dt_e
   sll_int32                             :: nc_x
   sll_int32                             :: nc_y

   nc_x = self%nc_x
   nc_y = self%nc_y

   dt_e = dt / self%e_0

   do j = 1, nc_y
      D_DX(hz(1:nc_x,j))
      ey(1:nc_x,j) = ey(1:nc_x,j) - dt_e * self%d_dx
   end do

   do i = 1, nc_x
      D_DY(hz(i,1:nc_y))
      ex(i,1:nc_y) = ex(i,1:nc_y) + dt_e * self%d_dy
   end do

   if (present(jx) .and. present(jy)) then
      ex = ex - dt_e * jx 
      ey = ey - dt_e * jy 
   end if

end subroutine ampere_te_2d_pstd

!> delete maxwell solver object
subroutine free_maxwell_2d_pstd(self)
type(maxwell_2d_pstd) :: self

#ifdef FFTW_F2003
if (c_associated(self%p_tmp_x)) call fftw_free(self%p_tmp_x)
if (c_associated(self%p_tmp_y)) call fftw_free(self%p_tmp_y)
#endif

call fftw_destroy_plan(self%fwx)
call fftw_destroy_plan(self%fwy)
call fftw_destroy_plan(self%bwx)
call fftw_destroy_plan(self%bwy)
!if (nthreads > 1) then
!   call dfftw_cleanup_threads(error)
!end if

end subroutine free_maxwell_2d_pstd

end module sll_maxwell_2d_pstd
