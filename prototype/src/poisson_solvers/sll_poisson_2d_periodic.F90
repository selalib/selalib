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

module sll_poisson_2d_periodic

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_mesh_2d.h"

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
  type(scalar_field_2d), pointer        :: sol, rhs
  sll_comp64, dimension(:,:), pointer :: rhst, ext, eyt
  type(fftclass)                      :: fftx, ffty
  class(sll_mapped_mesh_2d_base), pointer   :: mesh
  sll_real64, dimension(:,:), pointer :: kx, ky, k2
end type poisson_2d_periodic

!> Create an object to solve Poisson equation on 2D mesh with periodic
!> boundary conditions:
!> To solve the potential input parameter is a sll_real64 2d array.
!> To solve the potential derivatives (electric field), input parameter 
!> is a field_2d_vec2.
interface new_poisson_2d_periodic
   module procedure new_poisson_2d_periodic_E_fields
   module procedure new_poisson_2d_periodic_potential
end interface

!> Solve Poisson equation on 2D mesh with periodic boundary conditions. 
interface solve_poisson_2d_periodic
   module procedure solve_poisson_2d_periodic_E_fields
   module procedure solve_poisson_2d_periodic_potential
end interface

contains

function new_poisson_2d_periodic_potential(potential)

   type(poisson_2d_periodic),pointer :: new_poisson_2d_periodic_potential
   type(scalar_field_2d),      pointer :: potential
   class(sll_mapped_mesh_2d_base), pointer :: mesh
   sll_int32                         :: error
   sll_int32                         :: ncx
   sll_int32                         :: ncy

   mesh => potential%mesh
   ncx = GET_MESH_NC_ETA1(mesh)
   ncy = GET_MESH_NC_ETA2(mesh)

   SLL_ALLOCATE(new_poisson_2d_periodic_potential,                   error)
   !SLL_ALLOCATE(new_poisson_2d_periodic_potential%mesh,        error)
   SLL_ALLOCATE(new_poisson_2d_periodic_potential%rhst(ncy,ncx/2+1), error)
   SLL_ALLOCATE(new_poisson_2d_periodic_potential%kx(ncy,ncx/2+1),   error)
   SLL_ALLOCATE(new_poisson_2d_periodic_potential%ky(ncy,ncx/2+1),   error)
   SLL_ALLOCATE(new_poisson_2d_periodic_potential%k2(ncy,ncx/2+1),   error)

   new_poisson_2d_periodic_potential%mesh => mesh

   call initdfft(new_poisson_2d_periodic_potential%fftx, ncx)
   call initcfft(new_poisson_2d_periodic_potential%ffty, ncy)

   call wave_number_vectors(new_poisson_2d_periodic_potential)

end function new_poisson_2d_periodic_potential

function new_poisson_2d_periodic_e_fields(e_fields)

   type(poisson_2d_periodic), pointer :: new_poisson_2d_periodic_e_fields
   type(vector_field_2d),       pointer :: e_fields
   class(sll_mapped_mesh_2d_base),  pointer :: mesh
   sll_int32                          :: ncx
   sll_int32                          :: ncy
   sll_int32                          :: error

   mesh => e_fields%mesh
   ncx = GET_MESH_NC_ETA1(mesh)
   ncy = GET_MESH_NC_ETA2(mesh)

   SLL_ALLOCATE(new_poisson_2d_periodic_e_fields,                   error)
   !SLL_ALLOCATE(new_poisson_2d_periodic_e_fields%mesh,        error)
   SLL_ALLOCATE(new_poisson_2d_periodic_e_fields%rhst(ncy,ncx/2+1), error)
   SLL_ALLOCATE(new_poisson_2d_periodic_e_fields%ext(ncy,ncx/2+1),  error)
   SLL_ALLOCATE(new_poisson_2d_periodic_e_fields%eyt(ncy,ncx/2+1),  error)
   SLL_ALLOCATE(new_poisson_2d_periodic_e_fields%kx(ncy,ncx/2+1),   error)
   SLL_ALLOCATE(new_poisson_2d_periodic_e_fields%ky(ncy,ncx/2+1),   error)
   SLL_ALLOCATE(new_poisson_2d_periodic_e_fields%k2(ncy,ncx/2+1),   error)

   new_poisson_2d_periodic_e_fields%mesh => mesh
   
   call initdfft(new_poisson_2d_periodic_e_fields%fftx, ncx)
   call initcfft(new_poisson_2d_periodic_e_fields%ffty, ncy)

   call wave_number_vectors(new_poisson_2d_periodic_e_fields)

end function new_poisson_2d_periodic_e_fields

subroutine solve_poisson_2d_periodic_potential(this,sol,rhs,error)

   type(poisson_2d_periodic), intent(inout) :: this
   sll_real64, dimension(:,:), intent(in)   :: rhs
   sll_real64, dimension(:,:), intent(out)  :: sol
   sll_int32 , intent(out)                  :: error
   sll_int32                                :: ncx,ncy
   sll_int32                                :: i, j

   ncx = GET_FIELD_NC_ETA1(this)
   ncy = GET_FIELD_NC_ETA2(this)

   sol(1:ncx,1:ncy) = rhs(1:ncx,1:ncy)
   do j=1,ncy
      call dfftf(ncx,sol(1:ncx,j),this%fftx%coefd)
   end do

   call transpose_r2c(sol(1:ncx,1:ncy), this%rhst)

   do i=1,ncx/2+1
      call zfftf( ncy, this%rhst(:,i), this%ffty%coefcd)
   end do

   this%rhst = this%rhst / (this%kx*this%kx+this%ky*this%ky)

   do i=1,ncx/2+1
      call zfftb( ncy, this%rhst(:,i),  this%ffty%coefcd )
   end do

   call transpose_c2r(this%rhst, sol(1:ncx,1:ncy))

   do j=1,ncy
      call dfftb( ncx, sol(1:ncx,j),  this%fftx%coefd )
   end do

   sol(1:ncx,1:ncy) = sol(1:ncx,1:ncy) / (ncx*ncy)     ! normalize FFTs

   sol(ncx+1,:) = sol(1,:)
   sol(:,ncy+1) = sol(:,1)

   error = 0

end subroutine solve_poisson_2d_periodic_potential

subroutine solve_poisson_2d_periodic_E_fields(this,e_fields,rhs,error)

   type(poisson_2d_periodic), intent(inout) :: this
   type(scalar_field_2d), intent(in)          :: rhs
   type(vector_field_2d), intent(out)         :: e_fields
   sll_int32 , intent(out)                  :: error
   sll_int32                                :: ncx,ncy
   sll_int32                                :: i, j

   ncx = GET_FIELD_NC_ETA1(rhs)
   ncy = GET_FIELD_NC_ETA2(rhs)

   do j=1,ncy
      call dfftf(ncx, rhs%data(1:ncx,j),this%fftx%coefd)
   end do

   call transpose_r2c(rhs%data(1:ncx,1:ncy), this%rhst)

   do i=1,ncx/2+1
      call zfftf( ncy, this%rhst(:,i), this%ffty%coefcd)
   end do

   this%ext(1,1) = 0.0_f64
   this%eyt(1,1) = 0.0_f64

   this%ext = -cmplx(0.0_f64,this%kx/this%k2,kind=f64)*this%rhst
   this%eyt = -cmplx(0.0_f64,this%ky/this%k2,kind=f64)*this%rhst

   do i=1,ncx/2+1
      call zfftb( ncy, this%ext(:,i),  this%ffty%coefcd )
      call zfftb( ncy, this%eyt(:,i),  this%ffty%coefcd )
   end do

   call transpose_c2r(this%ext, e_fields%data(1:ncx,1:ncy)%v1)
   call transpose_c2r(this%eyt, e_fields%data(1:ncx,1:ncy)%v2)

   do j=1,ncy
      call dfftb( ncx, e_fields%data(1:ncx,j)%v1,  this%fftx%coefd )
      call dfftb( ncx, e_fields%data(1:ncx,j)%v2,  this%fftx%coefd )
   end do

   e_fields%data(1:ncx,1:ncy)%v1 = e_fields%data(1:ncx,1:ncy)%v1 / (ncx*ncy)
   e_fields%data(1:ncx,1:ncy)%v2 = e_fields%data(1:ncx,1:ncy)%v2 / (ncx*ncy)

   e_fields%data(ncx+1,:)%v1 = e_fields%data(1,:)%v1
   e_fields%data(:,ncy+1)%v1 = e_fields%data(:,1)%v1
   e_fields%data(ncx+1,:)%v2 = e_fields%data(1,:)%v2
   e_fields%data(:,ncy+1)%v2 = e_fields%data(:,1)%v2

   error = 0

end subroutine solve_poisson_2d_periodic_E_fields

subroutine wave_number_vectors(this)

   type(poisson_2d_periodic) :: this
   sll_int32  :: ik, jk
   sll_int32  :: ncx, ncy
   sll_real64 :: dx, dy
   sll_real64 :: kx, ky, kx0, ky0
   
   ncx = GET_FIELD_NC_ETA1(this)
   ncy = GET_FIELD_NC_ETA2(this)
   
   dx = GET_FIELD_DELTA_ETA1(this)
   dy = GET_FIELD_DELTA_ETA2(this)
   
   kx0 = 2._f64*sll_pi/(ncx*dx)
   ky0 = 2._f64*sll_pi/(ncy*dy)
   
   do ik=1,ncx/2+1
      kx  = (ik-1)*kx0
      do jk = 1, ncy/2
         ky  = (jk-1)*ky0
         this%kx(jk,ik) = kx
         this%ky(jk,ik) = ky
      end do
      do jk = ncy/2+1 , ncy     
         ky  = (jk-1-ncy)*ky0
         this%kx(jk,ik) = kx
         this%ky(jk,ik) = ky
      end do
   end do
   this%kx(1,1) = 1.0_f64
   
   this%k2 = this%kx*this%kx+this%ky*this%ky

end subroutine wave_number_vectors

subroutine transpose_r2c(real_array, comp_array)

   sll_real64, dimension(:,:), intent(in)  :: real_array
   sll_comp64, dimension(:,:), intent(out) :: comp_array
   sll_int32 :: i, j, n1, n2

   n1 = size(real_array,1)
   n2 = size(real_array,2)

   SLL_ASSERT(size(comp_array,1)==n2)
   SLL_ASSERT(size(comp_array,2)==n1/2+1)

   do j=1,n2
      comp_array(j,1) = cmplx(real_array(1,j),0._f64,kind=f64)
      do i=2, n1/2
         comp_array(j,i) = cmplx(real_array(2*i-2,j),real_array(2*i-1,j),kind=f64)
      end do
      comp_array(j,n1/2+1) = cmplx(real_array(n1,j),0._f64,kind=f64)
   end do

end subroutine transpose_r2c

!> convert complex array to real and transpose
subroutine transpose_c2r(comp_array, real_array)

   sll_comp64, dimension(:,:), intent(in) :: comp_array
   sll_real64, dimension(:,:), intent(out)  :: real_array
   sll_int32 :: i, j, n1, n2

   n1 = size(real_array,1)
   n2 = size(real_array,2)

   SLL_ASSERT((n2==size(comp_array,1)))
   SLL_ASSERT((size(comp_array,2)==n1/2+1))

   do j=1,n2
      real_array(1,j) = real(comp_array(j,1),kind=f64)
      do i=2,n1/2
         real_array(2*i-2,j) = real(comp_array(j,i),kind=f64)
         real_array(2*i-1,j) = dimag(comp_array(j,i))
      end do
      real_array(n1,j) = real(comp_array(j,n1/2+1),kind=f64)
   end do

end subroutine transpose_c2r

!> Delete the Poisson object
subroutine delete_poisson_2d_periodic(this)

   type(poisson_2d_periodic), pointer :: this
   this => null()

end subroutine delete_poisson_2d_periodic

end module sll_poisson_2D_periodic

