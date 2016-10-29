!**************************************************************
!  Copyright INRIA
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
!
!
!  Contact : Pierre Navaro http://wwww-irma.u-strasbg.fr/~navaro
!
!**************************************************************

#include "sll_fftw.h"

#define D_DX(FIELD)                                           \
self%d_dx = FIELD;                                            \
call sll_s_fft_exec_r2c_1d(self%fwx, self%d_dx, self%tmp_x);   \
self%tmp_x(2:nc_x/2+1)=-cmplx(0.0_f64,self%kx(2:nc_x/2+1),kind=f64)*self%tmp_x(2:nc_x/2+1); \
call sll_s_fft_exec_c2r_1d(self%bwx, self%tmp_x, self%d_dx);   \
self%d_dx = self%d_dx / nc_x

#define D_DY(FIELD)                                           \
self%d_dy = FIELD;                                            \
call sll_s_fft_exec_r2c_1d(self%fwy, self%d_dy, self%tmp_y);   \
self%tmp_y(2:nc_y/2+1)=-cmplx(0.0_f64,self%ky(2:nc_y/2+1),kind=f64)*self%tmp_y(2:nc_y/2+1); \
call sll_s_fft_exec_c2r_1d(self%bwy, self%tmp_y, self%d_dy);   \
self%d_dy = self%d_dy / nc_y

!> @ingroup maxwell_solvers
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
module sll_m_maxwell_2d_pstd
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_maxwell_solvers_macros.h"

use sll_m_constants, only: sll_p_pi
use sll_m_fft, only: sll_t_fft, &
  sll_s_fft_free,               &
  sll_s_fft_init_c2r_1d,        &
  sll_s_fft_init_r2c_1d,        &
  sll_s_fft_exec_c2r_1d,        &
  sll_s_fft_exec_r2c_1d

  implicit none

  public :: sll_t_maxwell_2d_pstd, &
  sll_s_init_maxwell_2d_pstd,      &
  sll_s_solve_maxwell_2d_pstd,     &
  sll_s_solve_ampere_2d_pstd,      &
  sll_s_free_maxwell_2d_pstd

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!> @brief 
!> Maxwell solver object
!> @details
!> We solve Maxwell system with PSTD numerical method. The type contains
!> information about FFT, mesh and physical properties.
type :: sll_t_maxwell_2d_pstd

   private
   sll_int32           :: nc_eta1      !< x cells number
   sll_int32           :: nc_eta2      !< y cells number
   sll_int32           :: polarization !< TE or TM
   sll_real64          :: e_0          !< electric conductivity
   sll_real64          :: mu_0         !< magnetic permeability
   sll_real64          :: c            !< speed of light
   sll_real64          :: eta1_min     !< left side 
   sll_real64          :: eta1_max     !< right side
   sll_real64          :: delta_eta1   !< step size
   sll_real64          :: eta2_min     !< bottom side
   sll_real64          :: eta2_max     !< top side
   sll_real64          :: delta_eta2   !< step size
   sll_real64, pointer :: d_dx(:)      !< field x derivative 
   sll_real64, pointer :: d_dy(:)      !< field y derivative
   sll_real64, pointer :: kx(:)        !< x wave number
   sll_real64, pointer :: ky(:)        !< y wave number
   type(sll_t_fft)     :: fwx          !< forward fft plan along x
   type(sll_t_fft)     :: fwy          !< forward fft plan along y
   type(sll_t_fft)     :: bwx          !< backward fft plan along x
   type(sll_t_fft)     :: bwy          !< backward fft plan along y
   sll_comp64, pointer :: tmp_x(:)     !< x fft transform
   sll_comp64, pointer :: tmp_y(:)     !< y fft transform

end type sll_t_maxwell_2d_pstd


!> Solve maxwell solver 2d cartesian periodic with PSTD scheme

!> Solve ampere equation using maxwell solver 2d cartesian periodic with PSTD scheme

!> Solve faraday equation using solver 2d cartesian periodic with PSTD scheme

!> Delete maxwell solver 2d cartesian periodic with PSTD scheme


contains

!> Initialize 2d maxwell solver on cartesian mesh with PSTD scheme
subroutine sll_s_init_maxwell_2d_pstd(self,xmin,xmax,nc_x,ymin,ymax,nc_y,polarization)

   type(sll_t_maxwell_2d_pstd) :: self         !< maxwell object
   sll_real64, intent(in)      :: xmin         !< x min
   sll_real64, intent(in)      :: xmax         !< x max
   sll_real64, intent(in)      :: ymin         !< y min
   sll_real64, intent(in)      :: ymax         !< y max
   sll_int32,  intent(in)      :: nc_x         !< x cells number
   sll_int32,  intent(in)      :: nc_y         !< y cells number
   sll_int32,  intent(in)      :: polarization !< TE or TM

   sll_int32                   :: error        !< error code
   sll_real64                  :: dx           !< x space step
   sll_real64                  :: dy           !< y space step
   sll_real64                  :: kx0
   sll_real64                  :: ky0
   sll_int32                   :: i, j

   self%nc_eta1 = nc_x
   self%nc_eta2 = nc_y
   self%polarization = polarization

   self%e_0  = 1._f64
   self%mu_0 = 1._f64

   SLL_ALLOCATE(self%d_dx(nc_x), error)
   SLL_ALLOCATE(self%d_dy(nc_y), error)

   SLL_ALLOCATE(self%tmp_x(nc_x/2+1), error)
   SLL_ALLOCATE(self%tmp_y(nc_y/2+1), error)

   call sll_s_fft_init_r2c_1d(self%fwx, nc_x, self%d_dx,  self%tmp_x)
   call sll_s_fft_init_c2r_1d(self%bwx, nc_x, self%tmp_x, self%d_dx)
   call sll_s_fft_init_r2c_1d(self%fwy, nc_y, self%d_dy,  self%tmp_y)
   call sll_s_fft_init_c2r_1d(self%bwy, nc_y, self%tmp_y, self%d_dy)

   SLL_ALLOCATE(self%kx(nc_x/2+1), error)
   SLL_ALLOCATE(self%ky(nc_y/2+1), error)
   
   dx = (xmax-xmin) / nc_x
   dy = (ymax-ymin) / nc_y

   kx0 = 2._f64*sll_p_pi/(nc_x*dx)
   ky0 = 2._f64*sll_p_pi/(nc_y*dy)

   do i=2,nc_x/2+1
      self%kx(i) = (i-1)*kx0
   end do
   self%kx(1) = 1.0_f64
   do j=2,nc_y/2+1
      self%ky(j) = (j-1)*ky0
   end do
   self%ky(1) = 1.0_f64

end subroutine sll_s_init_maxwell_2d_pstd

!> self routine exists only for testing purpose. Use ampere and faraday
!> in your appication.
subroutine sll_s_solve_maxwell_2d_pstd(self, fx, fy, fz, dt)

   type(sll_t_maxwell_2d_pstd), intent(inout)          :: self !< maxwell object
   sll_real64 , intent(inout), dimension(:,:) :: fx   !< Ex or Bx
   sll_real64 , intent(inout), dimension(:,:) :: fy   !< Ey or By
   sll_real64 , intent(inout), dimension(:,:) :: fz   !< Bz or Ez
   sll_real64 , intent(in)                    :: dt   !< time step

   IF ( self%polarization == TM_POLARIZATION) then
      call faraday_tm_2d_pstd(self, fx, fy, fz, 0.5*dt)   
      call ampere_tm_2d_pstd(self, fx, fy, fz, dt) 
      call faraday_tm_2d_pstd(self, fx, fy, fz, 0.5*dt)   
   end if

   IF ( self%polarization == TE_POLARIZATION) then
      call faraday_te_2d_pstd(self, fx, fy, fz, 0.5*dt)   
      call ampere_te_2d_pstd(self, fx, fy, fz, dt) 
      call faraday_te_2d_pstd(self, fx, fy, fz, 0.5*dt)   
   end if

end subroutine sll_s_solve_maxwell_2d_pstd

!> Solve Faraday equation
subroutine sll_s_solve_faraday_2d_pstd(self, fx, fy, fz, dt)
   type(sll_t_maxwell_2d_pstd),intent(inout) :: self    !< Maxwell object
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

end subroutine sll_s_solve_faraday_2d_pstd

subroutine sll_s_solve_ampere_2d_pstd(self, fx, fy, fz, dt, sx, sy)
   type(sll_t_maxwell_2d_pstd),intent(inout) :: self    !< Maxwell object
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

end subroutine sll_s_solve_ampere_2d_pstd

!> Solve faraday equation  (hx,hy,ez)
subroutine faraday_tm_2d_pstd(self, hx, hy, ez, dt)

   type(sll_t_maxwell_2d_pstd),intent(inout) :: self    !< Maxwell object
   sll_real64, dimension(:,:), intent(inout) :: hx      !< Magnetic field x
   sll_real64, dimension(:,:), intent(inout) :: hy      !< Magnetic field y
   sll_real64, dimension(:,:), intent(inout) :: ez      !< Electric field z
   sll_int32                                 :: nc_x    !< x cells number
   sll_int32                                 :: nc_y    !< y cells number
   sll_real64, intent(in)                    :: dt      !< time step

   sll_real64                                :: dt_mu
   sll_int32                                 :: i, j

   nc_x = self%nc_eta1
   nc_y = self%nc_eta2

   dt_mu = dt / self%mu_0 

   do i = 1, nc_x
      D_DY(ez(i,1:nc_y))
      hx(i,1:nc_y) = hx(i,1:nc_y) - dt_mu * self%d_dy
   end do

   do j = 1, nc_y
      D_DX(ez(1:nc_x,j))
      hy(1:nc_x,j) = hy(1:nc_x,j) + dt_mu * self%d_dx
   end do

   if (size(hx,2) == nc_y+1) hx(:,nc_y+1) = hx(:,1)
   if (size(hy,1) == nc_x+1) hy(nc_x+1,:) = hy(1,:)

end subroutine faraday_tm_2d_pstd

!> Solve faraday equation (ex,ey,hz)
subroutine faraday_te_2d_pstd(self, ex, ey, hz, dt)

   type(sll_t_maxwell_2d_pstd),intent(inout) :: self   !< maxwell object
   sll_real64, dimension(:,:), intent(inout) :: ex     !< electric field x
   sll_real64, dimension(:,:), intent(inout) :: ey     !< electric field y
   sll_real64, dimension(:,:), intent(inout) :: hz     !< magnetic field z
   sll_int32                                 :: nc_x   !< x cells number
   sll_int32                                 :: nc_y   !< y cells number
   sll_real64, intent(in)                    :: dt     !< time step

   sll_real64                                :: dt_mu
   sll_int32                                 :: i, j

   nc_x = self%nc_eta1
   nc_y = self%nc_eta2

   dt_mu = dt / self%mu_0 

   do i = 1, nc_x
      D_DY(ex(i,1:nc_y))
      hz(i,1:nc_y) = hz(i,1:nc_y) + dt_mu * self%d_dy
   end do

   do j = 1, nc_y
      D_DX(ey(1:nc_x,j))
      hz(1:nc_x,j) = hz(1:nc_x,j) - dt_mu * self%d_dx
   end do

   if (size(hz,1) == nc_x+1) hz(nc_x+1,:) = hz(1,:)
   if (size(hz,2) == nc_y+1) hz(:,nc_y+1) = hz(:,1)

end subroutine faraday_te_2d_pstd

!> Solve ampere maxwell equation (hx,hy,ez)
subroutine ampere_tm_2d_pstd(self, hx, hy, ez, dt, jz)

   type(sll_t_maxwell_2d_pstd),intent(inout) :: self   !< maxwell object
   sll_int32                                 :: nc_x   !< x cells number
   sll_int32                                 :: nc_y   !< y cells number
   sll_real64, dimension(:,:)                :: hx     !< magnetic field x
   sll_real64, dimension(:,:)                :: hy     !< magnetic field y
   sll_real64, dimension(:,:)                :: ez     !< electric field z
   sll_real64                                :: dt     !< time step
   sll_real64, dimension(:,:), optional      :: jz     !< current z

   sll_real64                                :: dt_e
   sll_int32                                 :: i, j

   nc_x = self%nc_eta1
   nc_y = self%nc_eta2

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

   if (size(ez,1) == nc_x+1) ez(nc_x+1,:) = ez(1,:)
   if (size(ez,2) == nc_y+1) ez(:,nc_y+1) = ez(:,1)

end subroutine ampere_tm_2d_pstd

!> Solve ampere maxwell equation (ex,ey,hz)
subroutine ampere_te_2d_pstd(self, ex, ey, hz, dt, jx, jy)

   type(sll_t_maxwell_2d_pstd),intent(inout) :: self   !< maxwell equation
   sll_real64, dimension(:,:)                :: ex     !< electric field x
   sll_real64, dimension(:,:)                :: ey     !< electric field y
   sll_real64, dimension(:,:)                :: hz     !< magnetic field z
   sll_real64                                :: dt     !< time step
   sll_real64, dimension(:,:), optional      :: jx     !< current x
   sll_real64, dimension(:,:), optional      :: jy     !< current y

   sll_real64                                :: dt_e
   sll_int32                                 :: nc_x
   sll_int32                                 :: nc_y
   sll_int32                                 :: i, j

   nc_x = self%nc_eta1
   nc_y = self%nc_eta2

   dt_e = dt / self%e_0

   do j = 1, nc_y
      D_DX(hz(1:nc_x,j))
      ey(1:nc_x,j) = ey(1:nc_x,j) - dt_e * self%d_dx
   end do

   do i = 1, nc_x
      D_DY(hz(i,1:nc_y))
      ex(i,1:nc_y) = ex(i,1:nc_y) + dt_e * self%d_dy
   end do

   If (present(jx) .and. present(jy)) then
      ex = ex - dt_e * jx 
      ey = ey - dt_e * jy 
   end if

   if (size(ex,1) == nc_x+1) ex(nc_x+1,:) = ex(1,:)
   if (size(ey,2) == nc_y+1) ey(:,nc_y+1) = ey(:,1)

end subroutine ampere_te_2d_pstd

!> delete maxwell solver object
subroutine sll_s_free_maxwell_2d_pstd(self)
type(sll_t_maxwell_2d_pstd) :: self

call sll_s_fft_free(self%fwx)
call sll_s_fft_free(self%fwy)
call sll_s_fft_free(self%bwx)
call sll_s_fft_free(self%bwy)

end subroutine sll_s_free_maxwell_2d_pstd

end module sll_m_maxwell_2d_pstd
