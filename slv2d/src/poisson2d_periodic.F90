module poisson2dpp_seq

#include "selalib.h"

use geometry_module
use, intrinsic :: iso_c_binding
implicit none
include 'fftw3.f03'

interface new
   module procedure initialize
end interface
interface dealloc
   module procedure dealloc_poisson2dpp
end interface
interface solve
   module procedure solve_potential
   module procedure solve_e_fields
end interface

type(C_PTR), private :: fw, bw
sll_int32, private :: nx, ny
sll_comp64, dimension(:,:), allocatable :: rhot
sll_comp64, dimension(:,:), allocatable :: ext
sll_comp64, dimension(:,:), allocatable :: eyt
sll_real64, private                     :: dx,dy
sll_real64, private                     :: kx0, kx
sll_real64, private                     :: ky0, ky
sll_comp64, dimension(:,:), pointer     :: rhst, ext, eyt
sll_real64, dimension(:,:), pointer     :: kx, ky, k2

contains

subroutine initialize(geomx, rho, error )
   sll_real64, dimension(:,:), intent(in) :: rho
   type(geometry),intent(in)  :: geomx
   sll_int32,  intent(out), optional :: error

   nc_x = geomx%nx
   nc_y = geomx%ny

   dx   = geomx%nx
   dy   = geomx%ny

   SLL_ALLOCATE(phi(nc_x,nc_y), error)

   SLL_ALLOCATE(rhst(nc_y,nc_x/2+1), error)
   SLL_ALLOCATE(ext (nc_y,nc_x/2+1), error)
   SLL_ALLOCATE(eyt (nc_y,nc_x/2+1), error)
   SLL_ALLOCATE(kx  (nc_y,nc_x/2+1), error)
   SLL_ALLOCATE(ky  (nc_y,nc_x/2+1), error)
   SLL_ALLOCATE(k2  (nc_y,nc_x/2+1), error)

   fw = fftw_plan_dft_r2c_2d(ny, nx, rho(1:nx,1:ny), rhot, FFTW_ESTIMATE);
   bw = fftw_plan_dft_c2r_2d(ny, nx, rhot, phi(1:nx,1:ny), FFTW_ESTIMATE)

   call wave_number_vectors(this)

end subroutine initialize

!> Solve Poisson equation on 2D mesh with periodic boundary conditions. 
!> return potential.
subroutine solve_potential(phi,rho)

   sll_real64, dimension(:,:), intent(in)    :: rho
   sll_real64, dimension(:,:), intent(out)   :: phi

   call fftw_execute_dft_r2c(fw, rho, rhst)
   rhst = rhst / k2
   call fftw_execute_dft_c2r(bw, rhst, phi)

   phi = phi / (nc_x*nc_y)     ! normalize FFTs

end subroutine solve_potential

!> Solve Poisson equation on 2D mesh with periodic boundary conditions. 
!> return electric fields.
subroutine solve_e_fields(rho,e_x,e_y,nrj)

   sll_real64, dimension(:,:), intent(in)   :: rho
   sll_real64, dimension(:,:), intent(out)  :: e_x
   sll_real64, dimension(:,:), intent(out)  :: e_y
   sll_real64, optional                     :: nrj

   call fftw_execute_dft_r2c(fw, rho, rhst)

   ext(1,1) = 0.0_f64
   eyt(1,1) = 0.0_f64
   ext = -cmplx(0.0_f64,kx/k2,kind=f64)*rhst
   eyt = -cmplx(0.0_f64,ky/k2,kind=f64)*rhst

   call fftw_execute_dft_c2r(bw, ext, e_x)
   call fftw_execute_dft_c2r(bw, eyt, e_y)

   e_x = e_x / (nc_x*nc_y)
   e_y = e_y / (nc_x*nc_y)

   if (present(nrj)) then 
      nrj=sum(e_x*e_x+e_y*e_y)*dx*dy
      if (nrj>1.e-30) then 
         nrj=0.5_wp*log(nrj)
      else
         nrj=-10**9
      endif
   end if

end subroutine solve_e_fields

subroutine wave_number_vectors()

   sll_int32  :: ik, jk
   sll_real64 :: kx1, kx0, ky0
   
   kx0 = 2._f64*sll_pi/(nx*dx)
   ky0 = 2._f64*sll_pi/(ny*dy)
   
   do ik=1,nc_x/2+1
      kx1 = (ik-1)*kx0
      do jk = 1, nc_y/2
         kx(jk,ik) = kx1
         ky(jk,ik) = (jk-1)*ky0
      end do
      do jk = nc_y/2+1 , nc_y     
         this%kx(jk,ik) = kx1
         this%ky(jk,ik) = (jk-1-nc_y)*ky0
      end do
   end do
   kx(1,1) = 1.0_f64
   
   k2 = kx*kx+ky*ky

end subroutine wave_number_vectors

end module poisson2dpp_seq
