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
!**************************************************************
!
!
!  Contact : Pierre Navaro http://wwww-irma.u-strasbg.fr/~navaro
!

#include "sll_fftw.h"

#define D_DX(field)                                           \
self%d_dx = field;                                            \
call sll_s_fft_exec_r2c_1d(self%fwx, self%d_dx, self%tmp_x);   \
self%tmp_x(2:nc_x/2+1)=-cmplx(0.0_f64,self%kx(2:nc_x/2+1),kind=f64)*self%tmp_x(2:nc_x/2+1);     \
call sll_s_fft_exec_c2r_1d(self%bwx, self%tmp_x, self%d_dx);   \
self%d_dx = self%d_dx / nc_x

#define D_DY(field)                                           \
self%d_dy = field;                                            \
call sll_s_fft_exec_r2c_1d(self%fwy, self%d_dy, self%tmp_y);   \
self%tmp_y(2:nc_y/2+1)=-cmplx(0.0_f64,self%ky(2:nc_y/2+1),kind=f64)*self%tmp_y(2:nc_y/2+1);     \
call sll_s_fft_exec_c2r_1d(self%bwy, self%tmp_y, self%d_dy);   \
self%d_dy = self%d_dy / nc_y

#define D_DZ(field)                                           \
self%d_dz = field;                                            \
call sll_s_fft_exec_r2c_1d(self%fwz, self%d_dz, self%tmp_z);   \
self%tmp_z(2:nc_z/2+1)=-cmplx(0.0_f64,self%kz(2:nc_z/2+1),kind=f64)*self%tmp_z(2:nc_z/2+1);     \
call sll_s_fft_exec_c2r_1d(self%bwz, self%tmp_z, self%d_dz);   \
self%d_dz = self%d_dz / nc_z

!> @ingroup maxwell_solvers
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
!>H_u\Big|^{n+1/2}_{i,j,k} = H_u\Big|^{n-1/2}_{i,j,k}  
!> - \frac{\Delta t}{\mu} \Big\{F_v^{-1}[-jk_vF_v(E_w)]|_{i,j,k}^n 
!> -F_w^{-1}[-jk_wF_w(E_v)]|_{i,j,k}^{n}\Big\},
!>\f$
!>
!>\f$ \displaystyle
!>E_u\Big|^{n+1}_{i,j,k} = E_u\Big|^{n}_{i,j,k}  
!> + \frac{\Delta t}{\varepsilon} \Big\{F_v^{-1}[-jk_vF_v(H_w)]|_{i,j,k}^{n+1/2} 
!> -F_w^{-1}[-jk_wF_w(H_v)]|_{i,j,k}^{n+1/2}\Big\},
!>\f$
!>
!>where \f$(u,v,w) = (x,y,z),(y,z,x),(z,x,y)\f$
module sll_m_maxwell_3d_pstd
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_maxwell_solvers_macros.h"

use sll_m_constants, only: sll_p_pi

use sll_m_fft

implicit none

public :: &
    faraday, ampere, &
    sll_o_delete, &
    sll_o_create

private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!> Initialize maxwell solver 2d cartesian periodic with PSTD scheme
interface sll_o_create
 module procedure new_maxwell_3d_pstd
end interface sll_o_create

!> Solve Ampere equation 3d cartesian periodic with PSTD scheme
interface sll_solve_ampere
 module procedure ampere
end interface sll_solve_ampere

!> Solve Faraday equation 3d cartesian periodic with PSTD scheme
interface sll_solve_faraday
 module procedure faraday
end interface sll_solve_faraday

!> Delete maxwell solver 3d cartesian periodic with PSTD scheme
interface sll_o_delete
 module procedure free_maxwell_3d_pstd
end interface sll_o_delete


!> Maxwell solver object
type, public :: sll_t_maxwell_3d_pstd
  private
  sll_int32                          :: nc_x         !< x cells number
  sll_int32                          :: nc_y         !< y cells number
  sll_int32                          :: nc_z         !< z cells number
  sll_real64, dimension(:), pointer  :: d_dx         !< field x derivative
  sll_real64, dimension(:), pointer  :: d_dy         !< field y derivative
  sll_real64, dimension(:), pointer  :: d_dz         !< field y derivative
  sll_real64, dimension(:), pointer  :: kx           !< x wave number
  sll_real64, dimension(:), pointer  :: ky           !< y wave number
  sll_real64, dimension(:), pointer  :: kz           !< z wave number
  type(sll_t_fft)                    :: fwx          !< forward fft plan along x
  type(sll_t_fft)                    :: fwy          !< forward fft plan along y
  type(sll_t_fft)                    :: fwz          !< forward fft plan along y
  type(sll_t_fft)                    :: bwx          !< backward fft plan along x
  type(sll_t_fft)                    :: bwy          !< backward fft plan along y
  type(sll_t_fft)                    :: bwz          !< backward fft plan along y
  sll_comp64               , pointer :: tmp_x(:)     !< x fft transform
  sll_comp64               , pointer :: tmp_y(:)     !< y fft transform
  sll_comp64               , pointer :: tmp_z(:)     !< y fft transform
  sll_int32                          :: polarization !< TE or TM
  sll_real64                         :: e_0          !< electric conductivity
  sll_real64                         :: mu_0         !< magnetic permeability
end type sll_t_maxwell_3d_pstd

contains

!> Initialize 2d maxwell solver on cartesian mesh with PSTD scheme
subroutine new_maxwell_3d_pstd(self,xmin,xmax,nc_x, &
                                    ymin,ymax,nc_y, &
                                    zmin,zmax,nc_z )

  type(sll_t_maxwell_3d_pstd) :: self   !< maxwell object
  sll_real64, intent(in)      :: xmin   !< x min
  sll_real64, intent(in)      :: xmax   !< x max
  sll_real64, intent(in)      :: ymin   !< y min
  sll_real64, intent(in)      :: ymax   !< y max
  sll_real64, intent(in)      :: zmin   !< z min
  sll_real64, intent(in)      :: zmax   !< z max
  sll_int32 , intent(in)      :: nc_x   !< x cells number
  sll_int32 , intent(in)      :: nc_y   !< y cells number
  sll_int32 , intent(in)      :: nc_z   !< z cells number

  sll_int32                   :: error  !< error code
  sll_real64                  :: dx     !< x space step
  sll_real64                  :: dy     !< y space step
  sll_real64                  :: dz     !< z space step
  sll_real64                  :: kx0
  sll_real64                  :: ky0
  sll_real64                  :: kz0

  sll_int32                   :: i, j, k

  self%nc_x = nc_x
  self%nc_y = nc_y
  self%nc_z = nc_z

  self%e_0  = 1.0_f64
  self%mu_0 = 1.0_f64

  SLL_ALLOCATE(self%tmp_x(nc_x/2+1), error)
  SLL_ALLOCATE(self%tmp_y(nc_y/2+1), error)
  SLL_ALLOCATE(self%tmp_z(nc_z/2+1), error)

  SLL_ALLOCATE(self%d_dx(nc_x), error)
  SLL_ALLOCATE(self%d_dy(nc_y), error)
  SLL_ALLOCATE(self%d_dz(nc_z), error)

  call sll_s_fft_init_r2c_1d(self%fwx, nc_x, self%d_dx,  self%tmp_x)
  call sll_s_fft_init_c2r_1d(self%bwx, nc_x, self%tmp_x, self%d_dx)
  call sll_s_fft_init_r2c_1d(self%fwy, nc_y, self%d_dy,  self%tmp_y)
  call sll_s_fft_init_c2r_1d(self%bwy, nc_y, self%tmp_y, self%d_dy)
  call sll_s_fft_init_r2c_1d(self%fwz, nc_z, self%d_dz,  self%tmp_z)
  call sll_s_fft_init_c2r_1d(self%bwz, nc_z, self%tmp_z, self%d_dz)

  SLL_ALLOCATE(self%kx(nc_x/2+1), error)
  SLL_ALLOCATE(self%ky(nc_y/2+1), error)
  SLL_ALLOCATE(self%kz(nc_z/2+1), error)
  
  dx = (xmax-xmin) / nc_x
  dy = (ymax-ymin) / nc_y
  dz = (zmax-zmin) / nc_z

  kx0 = 2._f64*sll_p_pi/(nc_x*dx)
  ky0 = 2._f64*sll_p_pi/(nc_y*dy)
  kz0 = 2._f64*sll_p_pi/(nc_z*dz)

  do i=2,nc_x/2+1
     self%kx(i) = (i-1)*kx0
  end do
  self%kx(1) = 1.0_f64
  do j=2,nc_y/2+1
     self%ky(j) = (j-1)*ky0
  end do
  self%ky(1) = 1.0_f64
  do k=2,nc_z/2+1
     self%kz(k) = (k-1)*kz0
  end do
  self%kz(1) = 1.0_f64

end subroutine new_maxwell_3d_pstd

!> Solve faraday equation  (hx,hy,ez)
subroutine faraday(self, ex, ey, ez, hx, hy, hz, dt)

  type(sll_t_maxwell_3d_pstd),  intent(inout) :: self  !< Maxwell object
  sll_real64, dimension(:,:,:), intent(inout) :: ex    !< Electric field x
  sll_real64, dimension(:,:,:), intent(inout) :: ey    !< Electric field y
  sll_real64, dimension(:,:,:), intent(inout) :: ez    !< Electric field z
  sll_real64, dimension(:,:,:), intent(inout) :: hx    !< Magnetic field x
  sll_real64, dimension(:,:,:), intent(inout) :: hy    !< Magnetic field y
  sll_real64, dimension(:,:,:), intent(inout) :: hz    !< Magnetic field z
  sll_int32                                   :: nc_x  !< x cells number
  sll_int32                                   :: nc_y  !< y cells number
  sll_int32                                   :: nc_z  !< z cells number
  sll_real64, intent(in)                      :: dt    !< time step

  sll_real64 :: dt_mu
  sll_int32  :: i, j, k

  nc_x = self%nc_x
  nc_y = self%nc_y
  nc_z = self%nc_z

  dt_mu = dt / self%mu_0 

  do k = 1, nc_z+1
     do i = 1, nc_x+1
        D_DY(ez(i,1:nc_y,k))
        hx(i,1:nc_y,k) = hx(i,1:nc_y,k) - dt_mu * self%d_dy
        D_DY(ex(i,1:nc_y,k))
        hz(i,1:nc_y,k) = hz(i,1:nc_y,k) + dt_mu * self%d_dy
     end do
  end do

  hx(:,nc_y+1,:) = hx(:,1,:)
  hz(:,nc_y+1,:) = hz(:,1,:)

  do j = 1, nc_y+1
     do i = 1, nc_x+1
        D_DZ(ey(i,j,1:nc_z))
        hx(i,j,1:nc_z) = hx(i,j,1:nc_z) + dt_mu * self%d_dz
        D_DZ(ex(i,j,1:nc_z))
        hy(i,j,1:nc_z) = hy(i,j,1:nc_z) - dt_mu * self%d_dz
     end do
  end do

  hx(:,:,nc_z+1) = hx(:,:,1)
  hy(:,:,nc_z+1) = hy(:,:,1)

  do k = 1, nc_z+1
     do j = 1, nc_y+1
        D_DX(ez(1:nc_x,j,k))
        hy(1:nc_x,j,k) = hy(1:nc_x,j,k) + dt_mu * self%d_dx
        D_DX(ey(1:nc_x,j,k))
        hz(1:nc_x,j,k) = hz(1:nc_x,j,k) - dt_mu * self%d_dx
     end do
  end do

  hy(nc_x+1,:,:) = hx(1,:,:)
  hz(nc_x+1,:,:) = hy(1,:,:)

end subroutine faraday

!> Solve ampere maxwell equation (hx,hy,ez)
subroutine ampere(self, hx, hy, hz, ex, ey, ez, dt, jx, jy, jz)

  type(sll_t_maxwell_3d_pstd),intent(inout) :: self !< maxwell object
  sll_real64, dimension(:,:,:)              :: hx   !< magnetic field x
  sll_real64, dimension(:,:,:)              :: hy   !< magnetic field y
  sll_real64, dimension(:,:,:)              :: hz   !< magnetic field z
  sll_real64, dimension(:,:,:)              :: ex   !< electric field x
  sll_real64, dimension(:,:,:)              :: ey   !< electric field y
  sll_real64, dimension(:,:,:)              :: ez   !< electric field z
  sll_real64                                :: dt   !< time step
  sll_real64, dimension(:,:,:), optional    :: jx   !< current x
  sll_real64, dimension(:,:,:), optional    :: jy   !< current y
  sll_real64, dimension(:,:,:), optional    :: jz   !< current z

  sll_int32                                 :: nc_x !< x cells number
  sll_int32                                 :: nc_y !< y cells number
  sll_int32                                 :: nc_z !< z cells number

  sll_real64 :: dt_e
  sll_int32  :: i, j, k

  nc_x = self%nc_x
  nc_y = self%nc_y
  nc_z = self%nc_z

  dt_e = dt / self%e_0

  do k = 1, nc_z
  do i = 1, nc_x
    D_DY(hz(i,1:nc_y,k))
    ex(i,1:nc_y,k) = ex(i,1:nc_y,k) + dt_e * self%d_dy
    D_DY(hx(i,1:nc_y,k))
    ez(i,1:nc_y,k) = ez(i,1:nc_y,k) - dt_e * self%d_dy
  end do
  end do

  ex(:,nc_y+1,:) = ex(:,1,:)
  ez(:,nc_y+1,:) = ez(:,1,:)

  do j = 1, nc_y
  do i = 1, nc_x
    D_DZ(hy(i,j,1:nc_z))
    ex(i,j,1:nc_z) = ex(i,j,1:nc_z) - dt_e * self%d_dz
    D_DZ(hx(i,j,1:nc_z))
    ey(i,j,1:nc_z) = ey(i,j,1:nc_z) + dt_e * self%d_dz
  end do
  end do

  ex(:,:,nc_z+1) = ex(:,:,1)
  ey(:,:,nc_z+1) = ey(:,:,1)

  do k = 1, nc_z
  do j = 1, nc_y
    D_DX(hz(1:nc_x,j,k))
    ey(1:nc_x,j,k) = ey(1:nc_x,j,k) - dt_e * self%d_dx
    D_DX(hy(1:nc_x,j,k))
    ez(1:nc_x,j,k) = ez(1:nc_x,j,k) + dt_e * self%d_dx
  end do
  end do

  ey(nc_x+1,:,:) = ex(1,:,:)
  ez(nc_x+1,:,:) = ey(1,:,:)

  if (present(jx) .and. present(jy) .and. present(jz)) then
    ex = ex - dt_e * jx 
    ey = ey - dt_e * jy 
    ez = ez - dt_e * jz 
  end if


end subroutine ampere

!> delete maxwell solver object
subroutine free_maxwell_3d_pstd(self)

  type(sll_t_maxwell_3d_pstd) :: self
  
  call sll_s_fft_free(self%fwx)
  call sll_s_fft_free(self%fwy)
  call sll_s_fft_free(self%fwz)
  call sll_s_fft_free(self%bwx)
  call sll_s_fft_free(self%bwy)
  call sll_s_fft_free(self%bwz)

end subroutine free_maxwell_3d_pstd

end module sll_m_maxwell_3d_pstd
