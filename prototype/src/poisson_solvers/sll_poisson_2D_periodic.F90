module sll_poisson_2D_periodic
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_mesh_types.h"

use numeric_constants
use fft_module

implicit none
private

public :: new, free, solve

type, public :: poisson2dp
sll_comp64, dimension(:,:), pointer :: rhot, ext, eyt
type(fftclass)                      :: fftx, ffty
logical                             :: transpose
type(mesh_descriptor_2d)            :: mesh
end type poisson2dp

sll_real64, parameter, private      :: zero = 0.0_f64

interface new
   module procedure new_poisson2dp
end interface
interface free
   module procedure free_poisson2dp
end interface
interface solve
   module procedure solve_poisson2dp
end interface

contains

subroutine new_poisson2dp(this,rho,mesh,error)
type(poisson2dp),intent(out)           :: this
sll_real64, dimension(:,:),intent(in)  :: rho
type(mesh_descriptor_2D), intent(in)   :: mesh
sll_int32, intent(out)                 :: error
sll_int32                              :: nx
sll_int32                              :: ny

this%transpose=.false.

this%mesh = mesh
nx = mesh%nc_eta1 + 1
ny = mesh%nc_eta2 + 1

! la taille totale de la zone en kx est nx/2
! ATTENTION : les tableaux concernees par ce decoupage sont rhot,
! ext et eyt. Ce sont des complexes. Pour cette raison leur taille est
! la moitie de celle des tableaux reels correspondants.

SLL_ALLOCATE(this%rhot(ny,nx/2+1), error)
SLL_ALLOCATE(this%ext(ny,nx/2+1),  error)
SLL_ALLOCATE(this%eyt(ny,nx/2+1),  error)

call initfft(this%fftx, rho      , nx)
call initfft(this%ffty, this%rhot, ny)

end subroutine new_poisson2dp

subroutine free_poisson2dp(this)
type(poisson2dp),intent(out) :: this
deallocate(this%rhot,this%ext,this%eyt)
end subroutine free_poisson2dp

subroutine solve_poisson2dp(this,ex,ey,rho,phi,error)
type(poisson2dp), intent(inout)  :: this
sll_real64, dimension(:,:), intent(in)  :: rho
sll_real64, dimension(:,:), intent(out) :: ex, ey
sll_real64, dimension(:,:), intent(out) :: phi
sll_int32 , intent(out)                 :: error
sll_int32                               :: nx,ny
sll_real64                              :: dx,dy
sll_real64                              :: kx0, kx, kx2
sll_real64                              :: ky0, ky, ky2
sll_int32                               :: ik, jk, i, j

nx = this%mesh%nc_eta1 + 1
ny = this%mesh%nc_eta2 + 1
dx = this%mesh%delta_eta1
dy = this%mesh%delta_eta2

phi = rho
call fft(this%fftx,phi)

do j=1,ny
   this%rhot(j,1) = dcmplx(phi(1,j),0.0)
   do i=2, nx/2
      this%rhot(j,i)=dcmplx(phi(2*i-2,j),phi(2*i-1,j))
   end do
   this%rhot(j,nx/2+1) = dcmplx(phi(2*nx/2,j),0.0)
end do

call fft(this%ffty,this%rhot)

kx0=2._f64*sll_pi/(nx*dx)
ky0=2._f64*sll_pi/(ny*dy)

this%ext(1,1)=0.
this%eyt(1,1)=0.

jk = 1
do ik=2,nx/2+1
   kx= (ik-1)*kx0
   kx2 = kx*kx
   this%ext(jk,ik)  = -dcmplx(zero,kx/kx2)*this%rhot(jk,ik)
   this%eyt(jk,ik)  = 0.
   this%rhot(jk,ik) = this%rhot(jk,ik)/kx2
end do

do ik=1,nx/2+1
   kx= (ik-1)*kx0
   kx2 = kx*kx
   do jk = 2, ny/2+1
      ky  = (jk-1)*ky0
      ky2 = kx2 +ky*ky
      this%ext(jk,ik)=-dcmplx(zero,kx/ky2)*this%rhot(jk,ik)
      this%eyt(jk,ik)=-dcmplx(zero,ky/ky2)*this%rhot(jk,ik)
      this%rhot(jk,ik)= this%rhot(jk,ik)/ky2
   end do
   do jk = ny/2+2 , ny     
      ky= (jk-1-ny)*ky0
      ky2= kx2 +ky*ky

      this%ext(jk,ik)=-dcmplx(zero,kx/ky2)*this%rhot(jk,ik)
      this%eyt(jk,ik)=-dcmplx(zero,ky/ky2)*this%rhot(jk,ik)
      this%rhot(jk,ik)= this%rhot(jk,ik)/ky2
   end do
end do

call fftinv(this%ffty,this%ext)
call fftinv(this%ffty,this%eyt)
call fftinv(this%ffty,this%rhot)

do j=1,ny

   ex(1,j)  = dble(this%ext(j,1))
   ey(1,j)  = dble(this%eyt(j,1))
   phi(1,j)  = dble(this%rhot(j,1))

   do i=2,nx/2

      ex(2*i-2,j) = dble(this%ext(j,i))
      ex(2*i-1,j) = dimag(this%ext(j,i))

      ey(2*i-2,j) = dble(this%eyt(j,i))
      ey(2*i-1,j) = dimag(this%eyt(j,i))

      phi(2*i-2,j) = dble(this%rhot(j,i))
      phi(2*i-1,j) = dimag(this%rhot(j,i))

   end do

   ex(nx,j) = dble(this%ext(j,nx/2+1))
   ey(nx,j) = dble(this%eyt(j,nx/2+1))
   phi(nx,j) = dble(this%rhot(j,nx/2+1))

end do

call fftinv(this%fftx,ex)
call fftinv(this%fftx,ey)
call fftinv(this%fftx,phi)

end subroutine solve_poisson2dp

end module sll_poisson_2D_periodic
