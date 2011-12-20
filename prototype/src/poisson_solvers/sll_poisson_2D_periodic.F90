module sll_poisson_2D_periodic
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_mesh_types.h"

use numeric_constants
use fft_module

implicit none
private

public :: new_poisson_2d_periodic
public :: solve
public :: delete

type, public :: poisson_2d_periodic
  type(field_2d_vec1), pointer        :: sol, rhs
  sll_comp64, dimension(:,:), pointer :: rhst, ext, eyt
  type(fftclass)                      :: fftx, ffty
  type(mesh_descriptor_2d), pointer   :: descriptor
end type poisson_2d_periodic

sll_real64, parameter, private      :: zero = 0.0_f64

interface delete
   module procedure delete_poisson_2d_periodic
end interface
interface solve
   module procedure solve_poisson_2d_periodic
end interface

contains

function new_poisson_2d_periodic(mesh)

   type(poisson_2d_periodic),pointer :: new_poisson_2d_periodic
   type(mesh_descriptor_2d), pointer :: mesh
   sll_int32                         :: error
   sll_int32                         :: nx
   sll_int32                         :: ny

   nx = GET_MESH_NC_ETA1(mesh) + 1
   ny = GET_MESH_NC_ETA2(mesh) + 1

   SLL_ALLOCATE(new_poisson_2d_periodic, error)
   SLL_ALLOCATE(new_poisson_2d_periodic%descriptor, error)
   SLL_ALLOCATE(new_poisson_2d_periodic%rhst(ny,nx/2+1), error)
   SLL_ALLOCATE(new_poisson_2d_periodic%ext(ny,nx/2+1),  error)
   SLL_ALLOCATE(new_poisson_2d_periodic%eyt(ny,nx/2+1),  error)

   new_poisson_2d_periodic%descriptor => mesh
   
   call initdfft(new_poisson_2d_periodic%fftx, nx)
   call initcfft(new_poisson_2d_periodic%ffty, ny)

end function new_poisson_2d_periodic

subroutine delete_poisson_2d_periodic(this)

   type(poisson_2d_periodic),intent(out) :: this
   deallocate(this%rhst,this%ext,this%eyt)

end subroutine delete_poisson_2d_periodic

subroutine solve_poisson_2d_periodic(this,ex,ey,rhs,sol,error)

   type(poisson_2d_periodic), intent(inout)  :: this
   sll_real64, dimension(:,:), intent(in)    :: rhs
   sll_real64, dimension(:,:), intent(out)   :: ex, ey
   sll_real64, dimension(:,:), intent(out)   :: sol
   sll_int32 , intent(out)                   :: error
   sll_int32                                 :: nx,ny
   sll_real64                                :: dx,dy
   sll_real64                                :: kx0, kx, kx2
   sll_real64                                :: ky0, ky, ky2
   sll_int32                                 :: ik, jk, i, j

   nx = GET_FIELD_NC_ETA1(this) + 1
   ny = GET_FIELD_NC_ETA2(this) + 1
   dx = GET_FIELD_DELTA_ETA1(this)
   dy = GET_FIELD_DELTA_ETA2(this)

   SLL_ASSERT(size(sol,1) >= size(rhs,1))
   SLL_ASSERT(size(sol,2) >= size(rhs,2))

   sol = rhs
   call fft(this%fftx,sol)
   
   do j=1,ny
      this%rhst(j,1) = cmplx(sol(1,j),0.0)
      do i=2, nx/2
         this%rhst(j,i)=cmplx(sol(2*i-2,j),sol(2*i-1,j))
      end do
      this%rhst(j,nx/2+1) = cmplx(sol(2*nx/2,j),0.0)
   end do

   call fft(this%ffty,this%rhst)

   kx0=2._f64*sll_pi/(nx*dx)
   ky0=2._f64*sll_pi/(ny*dy)

   this%ext(1,1)=0.
   this%eyt(1,1)=0.

   jk = 1
   do ik=2,nx/2+1
      kx= (ik-1)*kx0
      kx2 = kx*kx
      this%ext(jk,ik)  = -cmplx(zero,kx/kx2)*this%rhst(jk,ik)
      this%eyt(jk,ik)  = 0.
      this%rhst(jk,ik) = this%rhst(jk,ik)/kx2
   end do

   do ik=1,nx/2+1
      kx= (ik-1)*kx0
      kx2 = kx*kx
      do jk = 2, ny/2+1
         ky  = (jk-1)*ky0
         ky2 = kx2 +ky*ky
         this%ext(jk,ik)=-cmplx(zero,kx/ky2)*this%rhst(jk,ik)
         this%eyt(jk,ik)=-cmplx(zero,ky/ky2)*this%rhst(jk,ik)
         this%rhst(jk,ik)= this%rhst(jk,ik)/ky2
      end do
      do jk = ny/2+2 , ny     
         ky= (jk-1-ny)*ky0
         ky2= kx2 +ky*ky
   
         this%ext(jk,ik)=-cmplx(zero,kx/ky2)*this%rhst(jk,ik)
         this%eyt(jk,ik)=-cmplx(zero,ky/ky2)*this%rhst(jk,ik)
         this%rhst(jk,ik)= this%rhst(jk,ik)/ky2
      end do
   end do

   call fftinv(this%ffty,this%ext)
   call fftinv(this%ffty,this%eyt)
   call fftinv(this%ffty,this%rhst)

   do j=1,ny
   
      ex(1,j)  = dble(this%ext(j,1))
      ey(1,j)  = dble(this%eyt(j,1))
      sol(1,j)  = dble(this%rhst(j,1))
   
      do i=2,nx/2
   
         ex(2*i-2,j) = dble(this%ext(j,i))
         ex(2*i-1,j) = aimag(this%ext(j,i))
   
         ey(2*i-2,j) = dble(this%eyt(j,i))
         ey(2*i-1,j) = aimag(this%eyt(j,i))
   
         sol(2*i-2,j) = dble(this%rhst(j,i))
         sol(2*i-1,j) = aimag(this%rhst(j,i))
   
      end do

      ex(nx,j) = dble(this%ext(j,nx/2+1))
      ey(nx,j) = dble(this%eyt(j,nx/2+1))
      sol(nx,j) = dble(this%rhst(j,nx/2+1))

   end do

   call fftinv(this%fftx,ex)
   call fftinv(this%fftx,ey)
   call fftinv(this%fftx,sol)

   error = 0

end subroutine solve_poisson_2d_periodic

end module sll_poisson_2D_periodic
