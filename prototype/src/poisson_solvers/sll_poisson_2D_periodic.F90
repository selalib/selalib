!------------------------------------------------------------------------------
! SELALIB
!------------------------------------------------------------------------------
!
! MODULE: sll_poisson_2d_periodic
!
!> @author
!> Pierre Navaro
!>
!
! DESCRIPTION: 
!
!> @brief
!> Implements the Poisson solver in 2D with periodic boundary conditions
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
public :: solve_poisson_2d_periodic
public :: delete_poisson_2d_periodic

!> Object with data to solve Poisson equation on 2d domain with
!> periodic boundary conditions
type, public :: poisson_2d_periodic
  type(field_2d_vec1), pointer        :: sol, rhs
  sll_comp64, dimension(:,:), pointer :: rhst, ext, eyt
  type(fftclass)                      :: fftx, ffty
  type(mesh_descriptor_2d), pointer   :: descriptor
end type poisson_2d_periodic

!> Create an object to solve Poisson equation on 2D mesh with periodic
!> boundary conditions:
!> To solve the potential input parameter is a sll_real64 2d array.
!> To solve the potential derivatives (electric field), input parameter is a field_2d_vec2.
interface new_poisson_2d_periodic
   module procedure new_poisson_2d_periodic_E_fields
   module procedure new_poisson_2d_periodic_potential
end interface

!> Solve Poisson equation on 2D mesh with periodic 
!> boundary conditions. 
interface solve_poisson_2d_periodic
   module procedure solve_poisson_2d_periodic_E_fields
   module procedure solve_poisson_2d_periodic_potential
end interface

contains


function new_poisson_2d_periodic_potential(potential)

   type(poisson_2d_periodic),pointer :: new_poisson_2d_periodic_potential
   type(field_2D_vec1),      pointer :: potential
   type(mesh_descriptor_2d), pointer :: mesh
   sll_int32                         :: error
   sll_int32                         :: nx
   sll_int32                         :: ny

   mesh => potential%descriptor
   nx = GET_MESH_NC_ETA1(mesh) + 1
   ny = GET_MESH_NC_ETA2(mesh) + 1

   SLL_ALLOCATE(new_poisson_2d_periodic_potential, error)
   SLL_ALLOCATE(new_poisson_2d_periodic_potential%descriptor, error)
   SLL_ALLOCATE(new_poisson_2d_periodic_potential%rhst(ny,nx/2+1), error)

   new_poisson_2d_periodic_potential%descriptor => mesh
   
   call initdfft(new_poisson_2d_periodic_potential%fftx, nx)
   call initcfft(new_poisson_2d_periodic_potential%ffty, ny)

end function new_poisson_2d_periodic_potential

function new_poisson_2d_periodic_e_fields(e_fields)

   type(poisson_2d_periodic),pointer :: new_poisson_2d_periodic_e_fields
   type(field_2D_vec2), pointer :: e_fields
   type(mesh_descriptor_2d), pointer :: mesh
   sll_int32                         :: error
   sll_int32                         :: nx
   sll_int32                         :: ny

   mesh => e_fields%descriptor
   nx = GET_MESH_NC_ETA1(mesh) + 1
   ny = GET_MESH_NC_ETA2(mesh) + 1

   SLL_ALLOCATE(new_poisson_2d_periodic_e_fields, error)
   SLL_ALLOCATE(new_poisson_2d_periodic_e_fields%descriptor, error)
   SLL_ALLOCATE(new_poisson_2d_periodic_e_fields%rhst(ny,nx/2+1), error)
   SLL_ALLOCATE(new_poisson_2d_periodic_e_fields%ext(ny,nx/2+1),  error)
   SLL_ALLOCATE(new_poisson_2d_periodic_e_fields%eyt(ny,nx/2+1),  error)

   new_poisson_2d_periodic_e_fields%descriptor => mesh
   
   call initdfft(new_poisson_2d_periodic_e_fields%fftx, nx)
   call initcfft(new_poisson_2d_periodic_e_fields%ffty, ny)

end function new_poisson_2d_periodic_e_fields

!> Delete the Poisson object
subroutine delete_poisson_2d_periodic(this)

   type(poisson_2d_periodic),intent(out) :: this
   deallocate(this%rhst)

end subroutine delete_poisson_2d_periodic

subroutine solve_poisson_2d_periodic_potential(this,sol,rhs,error)

   type(poisson_2d_periodic), intent(inout) :: this
   sll_real64, dimension(:,:), intent(in)   :: rhs
   sll_real64, dimension(:,:), intent(out)  :: sol
   sll_int32 , intent(out)                  :: error
   sll_int32                                :: nx,ny
   sll_real64                               :: dx,dy
   sll_real64                               :: kx0, kx, kx2
   sll_real64                               :: ky0, ky, ky2
   sll_int32                                :: ik, jk, i, j

   nx = GET_FIELD_NC_ETA1(this) + 1
   ny = GET_FIELD_NC_ETA2(this) + 1
   dx = GET_FIELD_DELTA_ETA1(this)
   dy = GET_FIELD_DELTA_ETA2(this)

   sol = rhs
   call fft(this%fftx,sol)
   
   do j=1,ny
      this%rhst(j,1) = cmplx(sol(1,j),0._f64)
      do i=2, nx/2
         this%rhst(j,i)=cmplx(sol(2*i-2,j),sol(2*i-1,j))
      end do
      this%rhst(j,nx/2+1) = cmplx(sol(2*nx/2,j),0._f64)
   end do

   call fft(this%ffty,this%rhst)

   kx0=2._f64*sll_pi/(nx*dx)
   ky0=2._f64*sll_pi/(ny*dy)

   jk = 1
   do ik=2,nx/2+1
      kx= (ik-1)*kx0
      kx2 = kx*kx
      this%rhst(jk,ik) = this%rhst(jk,ik)/kx2
   end do

   do ik=1,nx/2+1
      kx= (ik-1)*kx0
      kx2 = kx*kx
      do jk = 2, ny/2+1
         ky  = (jk-1)*ky0
         ky2 = kx2 +ky*ky
         this%rhst(jk,ik)= this%rhst(jk,ik)/ky2
      end do
      do jk = ny/2+2 , ny     
         ky= (jk-1-ny)*ky0
         ky2= kx2 +ky*ky
         this%rhst(jk,ik)= this%rhst(jk,ik)/ky2
      end do
   end do

   call fftinv(this%ffty,this%rhst)

   do j=1,ny
      sol(1,j) = dble(this%rhst(j,1))
      do i=2,nx/2
         sol(2*i-2,j) = dble(this%rhst(j,i))
         sol(2*i-1,j) = aimag(this%rhst(j,i))
      end do
      sol(nx,j) = dble(this%rhst(j,nx/2+1))
   end do

   call fftinv(this%fftx,sol)

   error = 0

end subroutine solve_poisson_2d_periodic_potential

subroutine solve_poisson_2d_periodic_E_fields(this,e_fields,rhs,error)

   type(poisson_2d_periodic), intent(inout) :: this
   type(field_2d_vec1), intent(in)          :: rhs
   type(field_2d_vec2), intent(out)         :: e_fields
   sll_int32 , intent(out)                  :: error
   sll_int32                                :: nx,ny
   sll_real64                               :: dx,dy
   sll_real64                               :: kx0, kx, kx2
   sll_real64                               :: ky0, ky, ky2
   sll_int32                                :: ik, jk, i, j

   sll_real64, dimension(:,:), allocatable  :: tmp
   sll_real64, parameter :: zero = 0.0_f64

   nx = GET_FIELD_NC_ETA1(this) + 1
   ny = GET_FIELD_NC_ETA2(this) + 1
   dx = GET_FIELD_DELTA_ETA1(this)
   dy = GET_FIELD_DELTA_ETA2(this)

   SLL_ALLOCATE(tmp(nx,ny),error)

   tmp = rhs%data
   call fft(this%fftx,tmp)
   
   do j=1,ny
      this%rhst(j,1) = cmplx(tmp(1,j),0.0)
      do i=2, nx/2
         this%rhst(j,i) = cmplx(tmp(2*i-2,j),tmp(2*i-1,j))
      end do
      this%rhst(j,nx/2+1) = cmplx(tmp(2*nx/2,j),0.0_f64)
   end do

   call fft(this%ffty,this%rhst)

   kx0=2._f64*sll_pi/(nx*dx)
   ky0=2._f64*sll_pi/(ny*dy)

   this%ext(1,1) = 0.0_f64
   this%eyt(1,1) = 0.0_f64

   jk = 1
   do ik=2,nx/2+1
      kx  = (ik-1)*kx0
      kx2 = kx*kx
      this%ext(jk,ik)  = -cmplx(zero,kx/kx2)*this%rhst(jk,ik)
      this%eyt(jk,ik)  = 0.0_f64
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

   do j=1,ny
   
      e_fields%data(1,j)%v1  = dble(this%ext(j,1))
      e_fields%data(1,j)%v2  = dble(this%eyt(j,1))
   
      do i=2,nx/2
   
         e_fields%data(2*i-2,j)%v1 = dble(this%ext(j,i))
         e_fields%data(2*i-1,j)%v1 = aimag(this%ext(j,i))
   
         e_fields%data(2*i-2,j)%v2 = dble(this%eyt(j,i))
         e_fields%data(2*i-1,j)%v2 = aimag(this%eyt(j,i))
   
      end do

      e_fields%data(nx,j)%v1 = dble(this%ext(j,nx/2+1))
      e_fields%data(nx,j)%v2 = dble(this%eyt(j,nx/2+1))

   end do

   call fftinv(this%fftx,e_fields%data%v1)
   call fftinv(this%fftx,e_fields%data%v2)

   error = 0

end subroutine solve_poisson_2d_periodic_E_fields

end module sll_poisson_2D_periodic
