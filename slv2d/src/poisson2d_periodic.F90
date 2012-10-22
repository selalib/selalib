module poisson2d_periodic

#ifdef _FFTW
#include "selalib.h"

use geometry_module
use, intrinsic :: iso_c_binding
implicit none
include 'fftw3.f03'

interface new
   module procedure new_potential
   module procedure new_e_fields
end interface
interface solve
   module procedure solve_potential
   module procedure solve_e_fields
end interface


type :: poisson2d
   type(geometry)                       :: geom
   sll_comp64, dimension(:,:), pointer  :: rhot
   sll_comp64, dimension(:,:), pointer  :: ext
   sll_comp64, dimension(:,:), pointer  :: eyt
   sll_real64, dimension(:,:), pointer  :: kx, ky, k2
   type(C_PTR)                          :: fw, bw

end type

sll_int32, parameter :: nthreads = 2

contains

subroutine new_potential(self, rho, geomx, error )
   type(poisson2d)                           :: self
   sll_real64, dimension(:,:), intent(inout) :: rho
   type(geometry),intent(in)                 :: geomx
   sll_int32                                 :: error
   sll_int32                                 :: nx, ny
   sll_real64                                :: dx, dy
   sll_real64                                :: kx0, kx
   sll_real64                                :: ky0, ky
   sll_int32                                 :: ik, jk

   self%geom = geomx
   nx = geomx%nx
   ny = geomx%ny

   SLL_ALLOCATE(self%k2(nx/2+1,ny), error)
   SLL_ALLOCATE(self%rhot(nx/2+1,ny), error)
   call dfftw_init_threads(error)
   if (error == 0) stop 'FFTW CAN''T USE THREADS'
   call dfftw_plan_with_nthreads(nthreads)
   
   self%fw = fftw_plan_dft_r2c_2d(ny, nx, rho, self%rhot, FFTW_ESTIMATE)
   self%bw = fftw_plan_dft_c2r_2d(ny, nx, self%rhot, rho, FFTW_ESTIMATE)

   dx = self%geom%dx
   dy = self%geom%dy
   kx0=2._f64*sll_pi/(nx*dx)
   ky0=2._f64*sll_pi/(ny*dy)

   do ik=1,nx/2+1
      kx  = (ik-1)*kx0
      do jk = 1, ny/2
         ky  = (jk-1)*ky0
         self%k2(ik,jk) = kx*kx+ky*ky
      end do
      do jk = ny/2+1,ny     
         ky= (jk-1-ny)*ky0
         self%k2(ik,jk) = kx*kx+ky*ky
      end do
   end do
   self%k2(1,1) = 1.0_f64

end subroutine new_potential

subroutine new_e_fields(self, ex, ey, geomx, error )
   type(poisson2d) :: self
   sll_real64, dimension(:,:), intent(inout) :: ex
   sll_real64, dimension(:,:), intent(inout) :: ey
   type(geometry),intent(in)  :: geomx
   sll_int32  :: error
   sll_int32                                 :: nx, ny
   sll_int32  :: ik, jk
   sll_real64 :: kx1, kx0, ky0
   sll_real64                                :: dx, dy

   self%geom = geomx
   nx = geomx%nx
   ny = geomx%ny
   dx = geomx%dx
   dy = geomx%dy

   SLL_ALLOCATE(self%rhot(nx/2+1,ny), error)
   SLL_ALLOCATE(self%ext(nx/2+1,ny), error)
   SLL_ALLOCATE(self%eyt(nx/2+1,ny), error)
   SLL_ALLOCATE(self%kx (nx/2+1,ny), error)
   SLL_ALLOCATE(self%ky (nx/2+1,ny), error)
   SLL_ALLOCATE(self%k2 (nx/2+1,ny), error)

   self%fw = fftw_plan_dft_r2c_2d(ny, nx, ex, self%ext, FFTW_ESTIMATE)
   self%bw = fftw_plan_dft_c2r_2d(ny, nx, self%eyt, ey, FFTW_ESTIMATE)

   kx0 = 2._f64*sll_pi/(nx*dx)
   ky0 = 2._f64*sll_pi/(ny*dy)
   
   do ik=1,nx/2+1
      kx1 = (ik-1)*kx0
      do jk = 1, ny/2
         self%kx(ik,jk) = kx1
         self%ky(ik,jk) = (jk-1)*ky0
      end do
      do jk = ny/2+1 , ny     
         self%kx(ik,jk) = kx1
         self%ky(ik,jk) = (jk-1-ny)*ky0
      end do
   end do
   self%kx(1,1) = 1.0_f64
   
   self%k2 = self%kx*self%kx+self%ky*self%ky
   self%kx = self%kx/self%k2
   self%ky = self%ky/self%k2

   SLL_DEALLOCATE(self%k2, error)

end subroutine new_e_fields

!> Solve Poisson equation on 2D mesh with periodic boundary conditions. 
!> return potential.
subroutine solve_potential(self, rho, phi)

   type(poisson2d),intent(inout)  :: self
   sll_real64, dimension(:,:), intent(inout) :: rho
   sll_real64, dimension(:,:), intent(out)   :: phi
   sll_int32                                 :: nx, ny

   call fftw_execute_dft_r2c(self%fw, rho, self%rhot)

   self%rhot = self%rhot / self%k2

   call fftw_execute_dft_c2r(self%bw, self%rhot, phi)

   nx = self%geom%nx
   ny = self%geom%ny

   phi = phi / (nx*ny)     ! normalize

end subroutine solve_potential

!> Solve Poisson equation on 2D mesh with periodic boundary conditions. 
!> return electric fields.
subroutine solve_e_fields(self,e_x,e_y,rho,nrj)

   type(poisson2d),intent(inout)  :: self
   sll_real64, dimension(:,:), intent(inout) :: rho
   sll_real64, dimension(:,:), intent(out)   :: e_x
   sll_real64, dimension(:,:), intent(out)   :: e_y
   sll_real64, optional                      :: nrj
   sll_int32  :: nx, ny
   sll_real64 :: dx, dy

   nx = self%geom%nx
   ny = self%geom%ny

   call fftw_execute_dft_r2c(self%fw, rho, self%rhot)

   self%ext(1,1) = 0.0_f64
   self%eyt(1,1) = 0.0_f64
   self%ext = -cmplx(0.0_f64,self%kx,kind=f64)*self%rhot
   self%eyt = -cmplx(0.0_f64,self%ky,kind=f64)*self%rhot

   call fftw_execute_dft_c2r(self%bw, self%ext, e_x)
   call fftw_execute_dft_c2r(self%bw, self%eyt, e_y)

   e_x = e_x / (nx*ny)
   e_y = e_y / (nx*ny)

   if (present(nrj)) then 
      dx = self%geom%dx
      dy = self%geom%dy
      nrj=sum(e_x*e_x+e_y*e_y)*dx*dy
      if (nrj>1.e-30) then 
         nrj=0.5_wp*log(nrj)
      else
         nrj=-10**9
      endif
   end if

end subroutine solve_e_fields

subroutine free_poisson(self)
type(poisson2d) :: self
sll_int32       :: error
call dfftw_destroy_plan(self%fw)
call dfftw_destroy_plan(self%bw)
if (nthreads > 1) then
   call dfftw_cleanup_threads(error)
end if

end subroutine


#endif
end module poisson2d_periodic
