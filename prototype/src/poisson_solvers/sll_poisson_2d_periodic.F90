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

#ifdef _FFTPACK

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

#else

use, intrinsic :: iso_c_binding
implicit none
include 'fftw3.f03'

interface solve
   module procedure solve_potential
   module procedure solve_e_fields
end interface
interface free
   module procedure free_poisson
end interface

type :: poisson_2d_periodic
   sll_real64, dimension(:,:), pointer  :: kx, ky, k2
   type(C_PTR)                          :: fw, bw
   complex(C_DOUBLE_COMPLEX), dimension(:,:), pointer :: rhot
   complex(C_DOUBLE_COMPLEX), dimension(:,:), pointer :: ext
   complex(C_DOUBLE_COMPLEX), dimension(:,:), pointer :: eyt
   integer(C_SIZE_T) :: sz_rhot, sz_ext, sz_eyt
   type(C_PTR) :: p_rhot, p_ext, p_eyt
   sll_int32 :: nc_x, nc_y
   sll_real64  :: dx, dy
end type poisson_2d_periodic


#endif

contains

#ifdef _FFTPACK
!> Create an object to solve Poisson equation on 2D mesh with periodic
!> boundary conditions:
subroutine initialize(this, x_min, x_max, nc_x, &
                      y_min, y_max, nc_y, error )

   class(poisson_2d_periodic)        :: this
   sll_int32,  intent(in)            :: nc_x
   sll_int32,  intent(in)            :: nc_y
   sll_real64, intent(in)            :: x_min, x_max
   sll_real64, intent(in)            :: y_min, y_max
   sll_int32,  intent(out), optional :: error
   
   this%nc_x = nc_x
   this%nc_y = nc_y

   this%x_min = x_min
   this%x_max = x_max
   this%y_min = y_min
   this%y_max = y_max

   SLL_ALLOCATE(this%rhst(nc_y,nc_x/2+1), error)
   SLL_ALLOCATE(this%ext (nc_y,nc_x/2+1), error)
   SLL_ALLOCATE(this%eyt (nc_y,nc_x/2+1), error)
   SLL_ALLOCATE(this%kx  (nc_y,nc_x/2+1), error)
   SLL_ALLOCATE(this%ky  (nc_y,nc_x/2+1), error)
   SLL_ALLOCATE(this%k2  (nc_y,nc_x/2+1), error)

   call initdfft(this%fftx, nc_x)
   call initcfft(this%ffty, nc_y)

   call wave_number_vectors(this)

end subroutine initialize

!> Solve Poisson equation on 2D mesh with periodic boundary conditions. 
!> return potential.
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

!> Solve Poisson equation on 2D mesh with periodic boundary conditions. 
!> return electric fields.
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

   class(poisson_2d_periodic) :: this
   sll_int32  :: ik, jk
   sll_int32  :: nc_x, nc_y
   sll_real64 :: kx, ky, kx0, ky0
   
   nc_x = this%nc_x
   nc_y = this%nc_y
   
   kx0 = 2._f64*sll_pi/(this%x_max-this%x_min)
   ky0 = 2._f64*sll_pi/(this%y_max-this%y_min)
   
   do ik=1,nc_x/2+1
      kx = (ik-1)*kx0
      do jk = 1, nc_y/2
         ky = (jk-1)*ky0
         this%kx(jk,ik) = kx
         this%ky(jk,ik) = ky
      end do
      do jk = nc_y/2+1 , nc_y     
         ky = (jk-1-nc_y)*ky0
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

   sll_comp64, dimension(:,:), intent(in)  :: comp_array
   sll_real64, dimension(:,:), intent(out) :: real_array
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

#else

subroutine initialize(self, x_min, x_max, nc_x, &
                      y_min, y_max, nc_y, rho, error )

   type(poisson_2d_periodic) :: self
   sll_real64, dimension(:,:), intent(inout) :: rho
   sll_real64, intent(in) :: x_min, x_max, y_min, y_max
   sll_int32  :: error
   sll_int32                                 :: nc_x, nc_y
   sll_int32  :: ik, jk
   sll_real64 :: kx1, kx0, ky0

   self%nc_x = nc_x
   self%nc_y = nc_y

   self%dx   = (x_max-x_min) / nc_x
   self%dy   = (y_max-y_min) / nc_y

   self%sz_rhot = int((nc_x/2+1)*nc_y,C_SIZE_T)
   self%p_rhot = fftw_alloc_complex(self%sz_rhot)
   call c_f_pointer(self%p_rhot, self%rhot, [nc_x/2+1,nc_y])

   self%sz_ext = int((nc_x/2+1)*nc_y,C_SIZE_T)
   self%p_ext = fftw_alloc_complex(self%sz_ext)
   call c_f_pointer(self%p_ext, self%ext, [nc_x/2+1,nc_y])

   self%sz_eyt = int((nc_x/2+1)*nc_y,C_SIZE_T)
   self%p_eyt = fftw_alloc_complex(self%sz_eyt)
   call c_f_pointer(self%p_eyt, self%eyt, [nc_x/2+1,nc_y])

   SLL_ALLOCATE(self%kx (nc_x/2+1,nc_y), error)
   SLL_ALLOCATE(self%ky (nc_x/2+1,nc_y), error)
   SLL_ALLOCATE(self%k2 (nc_x/2+1,nc_y), error)

   !call dfftw_init_threads(error)
   !if (error == 0) stop 'FFTW CAN''T USE THREADS'
   !call dfftw_plan_with_nthreads(nthreads)

   self%fw = fftw_plan_dft_r2c_2d(nc_y,nc_x,rho(1:nc_x,1:nc_y), &
                                  self%ext,FFTW_ESTIMATE)
   self%bw = fftw_plan_dft_c2r_2d(nc_y,nc_x,self%eyt,           &
                                  rho(1:nc_x,1:nc_y),FFTW_ESTIMATE)

   kx0 = 2._f64*sll_pi/(x_max-x_min)
   ky0 = 2._f64*sll_pi/(y_max-y_min)
   
   do ik=1,nc_x/2+1
      kx1 = (ik-1)*kx0
      do jk = 1, nc_y/2
         self%kx(ik,jk) = kx1
         self%ky(ik,jk) = (jk-1)*ky0
      end do
      do jk = nc_y/2+1 , nc_y     
         self%kx(ik,jk) = kx1
         self%ky(ik,jk) = (jk-1-nc_y)*ky0
      end do
   end do
   self%kx(1,1) = 1.0_f64
   
   self%k2 = self%kx*self%kx+self%ky*self%ky
   self%kx = self%kx/self%k2
   self%ky = self%ky/self%k2

   !SLL_DEALLOCATE(self%k2, error)

end subroutine initialize

!> Solve Poisson equation on 2D mesh with periodic boundary conditions. 
!> return potential.
subroutine solve_potential(self, rho, phi)

   type(poisson_2d_periodic),intent(inout)  :: self
   sll_real64, dimension(:,:), intent(inout) :: rho
   sll_real64, dimension(:,:), intent(out)   :: phi
   sll_int32                                 :: nc_x, nc_y

   call fftw_execute_dft_r2c(self%fw, rho, self%rhot)

   self%rhot = self%rhot / self%k2

   call fftw_execute_dft_c2r(self%bw, self%rhot, phi)

   nc_x = self%nc_x
   nc_y = self%nc_y

   phi = phi / (nc_x*nc_y)     ! normalize

end subroutine solve_potential

!> Solve Poisson equation on 2D mesh with periodic boundary conditions. 
!> return electric fields.
subroutine solve_e_fields(self,e_x,e_y,rho,nrj)

   type(poisson_2d_periodic),intent(inout)  :: self
   sll_real64, dimension(:,:), intent(inout) :: rho
   sll_real64, dimension(:,:), intent(out)   :: e_x
   sll_real64, dimension(:,:), intent(out)   :: e_y
   sll_real64, optional                      :: nrj
   sll_int32  :: nc_x, nc_y
   sll_real64 :: dx, dy

   nc_x = self%nc_x
   nc_y = self%nc_y

   call fftw_execute_dft_r2c(self%fw, rho, self%rhot)

   self%ext(1,1) = 0.0_f64
   self%eyt(1,1) = 0.0_f64
   self%ext = -cmplx(0.0_f64,self%kx,kind=f64)*self%rhot
   self%eyt = -cmplx(0.0_f64,self%ky,kind=f64)*self%rhot

   call fftw_execute_dft_c2r(self%bw, self%ext, e_x)
   call fftw_execute_dft_c2r(self%bw, self%eyt, e_y)

   e_x = e_x / (nc_x*nc_y)
   e_y = e_y / (nc_x*nc_y)

   if (present(nrj)) then 
      dx = self%dx
      dy = self%dy
      nrj=sum(e_x*e_x+e_y*e_y)*dx*dy
      if (nrj>1.e-30) then 
         nrj=0.5_f64*log(nrj)
      else
         nrj=-10**9
      endif
   end if

end subroutine solve_e_fields

subroutine free_poisson(self)
type(poisson_2d_periodic) :: self
call fftw_free(self%p_rhot)
if (c_associated(self%p_ext)) call fftw_free(self%p_ext)
if (c_associated(self%p_eyt)) call fftw_free(self%p_eyt)
call dfftw_destroy_plan(self%fw)
call dfftw_destroy_plan(self%bw)
!if (nthreads > 1) then
   !call dfftw_cleanup_threads(error)
!end if

end subroutine


#endif
end module sll_poisson_2D_periodic
