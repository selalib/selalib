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

use numeric_constants
use fft_module

implicit none
private



!> Object with data to solve Poisson equation on 2d domain with
!> periodic boundary conditions
type, public :: poisson_2d_periodic
  sll_int32                           :: nc_x, nc_y
  sll_real64                          :: x_min, x_max
  sll_real64                          :: y_min, y_max
  sll_comp64, dimension(:,:), pointer :: rhst, ext, eyt
  sll_real64, dimension(:,:), pointer :: kx, ky, k2
  type(fftclass)                      :: fftx, ffty
contains
   procedure :: initialize
   procedure :: solve_e_fields
   procedure :: solve_potential
   generic   :: solve => solve_e_fields, solve_potential
end type poisson_2d_periodic

!> Create an object to solve Poisson equation on 2D mesh with periodic
!> boundary conditions:
!> To solve the potential input parameter is a sll_real64 2d array.
!> To solve the potential derivatives (electric field), input parameter 
!> is a field_2d_vec2.

!> Solve Poisson equation on 2D mesh with periodic boundary conditions. 

contains

subroutine initialize(this, x_min, x_max, nc_x, &
               y_min, y_max, nc_y, error )

   class(poisson_2d_periodic)         :: this
   sll_int32, intent(in)              :: nc_x
   sll_int32, intent(in)              :: nc_y
   sll_real64, intent(in)             :: x_min, x_max
   sll_real64, intent(in)             :: y_min, y_max
   sll_int32, intent(out), optional   :: error
   
   this%nc_x = nc_x
   this%nc_y = nc_y

   this%x_min = x_min
   this%x_max = x_max
   this%y_min = y_min
   this%y_max = y_max

   SLL_ALLOCATE(this%rhst(nc_y,nc_x/2+1), error)
   SLL_ALLOCATE(this%ext(nc_y,nc_x/2+1), error)
   SLL_ALLOCATE(this%eyt(nc_y,nc_x/2+1), error)
   SLL_ALLOCATE(this%kx(  nc_y,nc_x/2+1), error)
   SLL_ALLOCATE(this%ky(  nc_y,nc_x/2+1), error)
   SLL_ALLOCATE(this%k2(  nc_y,nc_x/2+1), error)

   call initdfft(this%fftx, nc_x)
   call initcfft(this%ffty, nc_y)

   call wave_number_vectors(this)

end subroutine initialize

subroutine solve_potential(this,sol,rhs)

   class(poisson_2d_periodic)                :: this
   sll_real64, dimension(:,:), intent(in)    :: rhs
   sll_real64, dimension(:,:), intent(out)   :: sol
   sll_int32                                 :: nc_x, nc_y
   sll_int32                                 :: i, j

   nc_x = this%nc_x
   nc_y = this%nc_y


   sol(1:nc_x,1:nc_y) = rhs(1:nc_x,1:nc_y)
   do j=1,nc_y
      call dfftf(nc_x, sol(1:nc_x,j), this%fftx%coefd)
   end do

   call transpose_r2c(sol(1:nc_x,1:nc_y), this%rhst)

   do i=1,nc_x/2+1
      call zfftf( nc_y, this%rhst(:,i), this%ffty%coefcd)
   end do

   this%rhst = this%rhst / this%k2

   do i=1,nc_x/2+1
      call zfftb( nc_y, this%rhst(:,i), this%ffty%coefcd )
   end do

   call transpose_c2r(this%rhst, sol(1:nc_x,1:nc_y))

   do j=1,nc_y
      call dfftb( nc_x, sol(1:nc_x,j),  this%fftx%coefd )
   end do

   sol(1:nc_x,1:nc_y) = sol(1:nc_x,1:nc_y) / (nc_x*nc_y)     ! normalize FFTs

   sol(nc_x+1,:) = sol(1,:)
   sol(:,nc_y+1) = sol(:,1)

end subroutine solve_potential

subroutine solve_e_fields(this,field_x,field_y,rhs)

   class(poisson_2d_periodic)               :: this
   sll_real64, dimension(:,:), intent(in)   :: rhs
   sll_real64, dimension(:,:), intent(out)  :: field_x
   sll_real64, dimension(:,:), intent(out)  :: field_y
   sll_int32                                :: nc_x, nc_y
   sll_int32                                :: i, j

   nc_x = this%nc_x
   nc_y = this%nc_y

   this%rhst = 0.0_f64
   this%ext  = 0.0_f64
   this%eyt  = 0.0_f64
   field_x   = 0.0_f64
   field_y   = 0.0_f64

   do j=1,nc_y
      call dfftf(nc_x, rhs(1:nc_x,j), this%fftx%coefd)
   end do

   call transpose_r2c(rhs(1:nc_x,1:nc_y), this%rhst)

   do i=1,nc_x/2+1
      call zfftf( nc_y, this%rhst(:,i), this%ffty%coefcd)
   end do

   this%ext(1,1) = 0.0_f64
   this%eyt(1,1) = 0.0_f64
   this%ext = -cmplx(0.0_f64,this%kx/this%k2,kind=f64)*this%rhst
   this%eyt = -cmplx(0.0_f64,this%ky/this%k2,kind=f64)*this%rhst

   do i=1,nc_x/2+1
      call zfftb( nc_y, this%ext(:,i), this%ffty%coefcd )
      call zfftb( nc_y, this%eyt(:,i), this%ffty%coefcd )
   end do

   call transpose_c2r(this%ext, field_x(1:nc_x,1:nc_y))
   call transpose_c2r(this%eyt, field_y(1:nc_x,1:nc_y))

   do j=1,nc_y
      call dfftb( nc_x, field_x(1:nc_x,j), this%fftx%coefd )
      call dfftb( nc_x, field_y(1:nc_x,j), this%fftx%coefd )
   end do

   field_x(1:nc_x,1:nc_y) = field_x(1:nc_x,1:nc_y) / (nc_x*nc_y)
   field_y(1:nc_x,1:nc_y) = field_y(1:nc_x,1:nc_y) / (nc_x*nc_y)

   field_x(nc_x+1,:) = field_x(1,:)
   field_x(:,nc_y+1) = field_x(:,1)
   field_y(nc_x+1,:) = field_y(1,:)
   field_y(:,nc_y+1) = field_y(:,1)

end subroutine solve_e_fields

subroutine wave_number_vectors(this)

   type(poisson_2d_periodic) :: this
   sll_int32  :: ik, jk
   sll_int32  :: nc_x, nc_y
   sll_real64 :: kx, ky, kx0, ky0
   
   nc_x = this%nc_x
   nc_y = this%nc_y
   
   kx0 = 2._f64*sll_pi/(this%x_max-this%x_min)
   ky0 = 2._f64*sll_pi/(this%y_max-this%y_min)
   
   do ik=1,nc_x/2+1
      kx  = (ik-1)*kx0
      do jk = 1, nc_y/2
         ky  = (jk-1)*ky0
         this%kx(jk,ik) = kx
         this%ky(jk,ik) = ky
      end do
      do jk = nc_y/2+1 , nc_y     
         ky  = (jk-1-nc_y)*ky0
         this%kx(jk,ik) = kx
         this%ky(jk,ik) = ky
      end do
   end do
   this%kx(1,1) = 1.0_f64
   
   this%k2 = this%kx*this%kx+this%ky*this%ky

end subroutine wave_number_vectors

!> convert real array to complex and transpose
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

end module sll_poisson_2D_periodic

