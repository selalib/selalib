!> @author
!> Pierre Navaro
!> @brief
!> Implements the Poisson solver in 2D with periodic boundary conditions
!> @details
!> This module uses FFTPACK library
module sll_poisson_2d_periodic

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_constants.h"

use sll_poisson_solvers

implicit none
private

sll_int32, private :: i, j

!> fft type use to do fft with fftpack library
type, public :: fftclass
   sll_real64, dimension(:), pointer :: coefc !< data for complex fft
   sll_real64, dimension(:), pointer :: work  !< work data for fft
   sll_real64, dimension(:), pointer :: workc !< workc complex
   sll_real64, dimension(:), pointer :: coefd !< data for double fft
   sll_real64, dimension(:), pointer :: workd !< work data
   sll_real64, dimension(:), pointer :: coefcd!< data for complex fft
   sll_int32  :: n  !< number of samples in each sequence
end type fftclass

!> Object with data to solve Poisson equation on 2d domain with
!> periodic boundary conditions
type, public :: poisson_2d_periodic
  sll_int32   :: nc_x  !< number of cells direction x
  sll_int32   :: nc_y  !< number of cells direction y
  sll_real64  :: dx    !< step size direction x
  sll_real64  :: dy    !< step size direction y
  sll_real64  :: x_min !< left corner direction x
  sll_real64  :: x_max !< right corner direction x
  sll_real64  :: y_min !< left corner direction y
  sll_real64  :: y_max !< right corner direction y
  sll_comp64, dimension(:,:), pointer :: rhst !< rhs fft
  sll_comp64, dimension(:,:), pointer :: ext  !< x electric field fft
  sll_comp64, dimension(:,:), pointer :: eyt  !< y electric field fft
  sll_real64, dimension(:,:), pointer :: kx   !< wave number x
  sll_real64, dimension(:,:), pointer :: ky   !< wave number y
  sll_real64, dimension(:,:), pointer :: k2   !< \f$ k_x^2+k_y^2 \f$
  type(fftclass)                      :: fftx !< fft plan in direction x
  type(fftclass)                      :: ffty !< fft plan in direction y
end type poisson_2d_periodic

!> Create a new poisson solver on 1d mesh
interface new
  module procedure new_poisson_2d_periodic_fftpack
end interface


!> Initialize
interface initialize
  module procedure initialize_poisson_2d_periodic_fftpack
end interface

!> Solve
interface solve
   module procedure solve_potential_poisson_2d_periodic_fftpack
   module procedure solve_e_fields_poisson_2d_periodic_fftpack
end interface

!interface delete
   !module procedure free_poisson_2d_periodic_fftpack
!end interface


interface fft
   module procedure doubfft, doubcfft
end interface
interface fftinv
   module procedure doubfftinv,  doubcfftinv
end interface

public :: initialize, new, solve

contains




  !> Create a new solver
  function new_poisson_2d_periodic_fftpack(&
    x_min, &
    x_max, &
    nc_x, &
    y_min, &
    y_max, &
    nc_y, &
    error) &
    result(this)
   type(poisson_2d_periodic),pointer :: this   !< self object
   sll_int32,  intent(in)    :: nc_x   !< number of cells direction x
   sll_int32,  intent(in)    :: nc_y   !< number of cells direction y
   sll_real64, intent(in)    :: x_min  !< left corner direction x
   sll_real64, intent(in)    :: x_max  !< right corner direction x
   sll_real64, intent(in)    :: y_min  !< left corner direction y
   sll_real64, intent(in)    :: y_max  !< right corner direction y
   sll_int32,  intent(out)   :: error  !< error code

   SLL_ALLOCATE(this, error)
   call initialize_poisson_2d_periodic_fftpack( &
           this, x_min, x_max, nc_x, y_min, y_max, nc_y, error )

  end function new_poisson_2d_periodic_fftpack 



!> Create an object to solve Poisson equation on 2D mesh with periodic
!> boundary conditions:
subroutine initialize_poisson_2d_periodic_fftpack( &
           this, x_min, x_max, nc_x, y_min, y_max, nc_y, error )

   type(poisson_2d_periodic) :: this   !< self object
   sll_int32,  intent(in)    :: nc_x   !< number of cells direction x
   sll_int32,  intent(in)    :: nc_y   !< number of cells direction y
   sll_real64, intent(in)    :: x_min  !< left corner direction x
   sll_real64, intent(in)    :: x_max  !< right corner direction x
   sll_real64, intent(in)    :: y_min  !< left corner direction y
   sll_real64, intent(in)    :: y_max  !< right corner direction y
   sll_int32,  intent(out)   :: error  !< error code
   
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

end subroutine initialize_poisson_2d_periodic_fftpack

!> Solve Poisson equation on 2D mesh with periodic boundary conditions. 
!> return potential.
subroutine solve_potential_poisson_2d_periodic_fftpack(this,sol,rhs)

   type(poisson_2d_periodic)               :: this !< self object
   sll_real64, dimension(:,:), intent(in)  :: rhs  !< charge density
   sll_real64, dimension(:,:), intent(out) :: sol  !< electric potential
   sll_int32                               :: nc_x !< number of cells direction x
   sll_int32                               :: nc_y !< number of cells direction y
   sll_int32                               :: i, j

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

end subroutine solve_potential_poisson_2d_periodic_fftpack

!> Solve Poisson equation on 2D mesh with periodic boundary conditions. 
!> return electric fields.
subroutine solve_e_fields_poisson_2d_periodic_fftpack(this,field_x,field_y,rhs,nrj)

   type(poisson_2d_periodic)               :: this    !< self object
   sll_real64, dimension(:,:), intent(in)  :: rhs     !< charge density
   sll_real64, dimension(:,:), intent(out) :: field_x !< electric field direction x
   sll_real64, dimension(:,:), intent(out) :: field_y !< electric field direction y
   sll_int32                               :: nc_x    !< number of cells direction x
   sll_int32                               :: nc_y    !< number of cells direction y
   sll_int32                               :: i, j
   sll_real64, optional                    :: nrj     !< \f$ \sqrt{e_x^2+e_y^2} \f$

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

   if (present(nrj)) then 
      nrj=sum(field_x(1:nc_x,1:nc_y)*field_x(1:nc_x,1:nc_y) &
        +field_y(1:nc_x,1:nc_y)*field_y(1:nc_x,1:nc_y))*this%dx*this%dy
!      if (nrj>1.e-30) then 
!         !nrj=0.5_f64*log(nrj)
!      else
!         nrj=-10**9
!      endif
   end if

end subroutine solve_e_fields_poisson_2d_periodic_fftpack

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

subroutine initdfft(this,l)

   type(fftclass) :: this
   sll_int32 :: l 
   this%n = l 
   allocate(this%coefd(2*this%n+15))
   call dffti(this%n,this%coefd)

end subroutine initdfft

subroutine initcfft(this,l)

   type(fftclass) :: this
   sll_int32 :: l 
   this%n = l
   allocate(this%coefcd(4*this%n+15))
   call zffti(this%n,this%coefcd)

end subroutine initcfft

subroutine doubfft(this,array)

   type(fftclass) :: this
   sll_real64, dimension(:,:) :: array

   do i=1, size(array,2)   ! number of 1d transforms
      call dfftf( this%n, array(:,i), this%coefd)
   end do

   array = array /this%n      ! normalize FFT

end subroutine doubfft

subroutine doubcfft(this,array)

   type(fftclass) :: this
   sll_comp64, dimension(:,:) :: array

   do i=1, size(array,2)   ! number of 1d transforms
      call zfftf( this%n, array(:,i), this%coefcd)
   end do

   array = array /this%n      ! normalize FFT

end subroutine doubcfft

subroutine doubfftinv(this,array)

   type(fftclass) :: this
   sll_real64, dimension(:,:) :: array

   do i=1, size(array,2)   ! number of 1d transforms
      call dfftb( this%n, array(:,i),  this%coefd )
   end do

end subroutine doubfftinv

subroutine doubcfftinv(this,array)

   type(fftclass) :: this
   sll_comp64, dimension(:,:) :: array

   do i=1, size(array,2)   ! number of 1d transforms
      call zfftb( this%n, array(:,i),  this%coefcd )
   end do

end subroutine doubcfftinv

end module sll_poisson_2D_periodic
