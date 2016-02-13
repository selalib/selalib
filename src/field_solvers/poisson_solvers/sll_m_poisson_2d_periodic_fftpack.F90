#ifndef DOXYGEN_SHOULD_SKIP_THIS

!> @ingroup poisson_solvers
!> @brief
!> Implements the Poisson solver in 2D with periodic boundary conditions
!> @details
!> This module uses FFTPACK library
module sll_m_poisson_2d_periodic_fftpack

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

! use F77_fftpack, only: &
!   dfftb, &
!   dfftf, &
!   dffti, &
!   zfftb, &
!   zfftf, &
!   zffti

  use sll_m_constants, only: &
    sll_p_pi

  implicit none

  public :: &
    sll_o_initialize, &
    sll_o_new, &
    sll_t_poisson_2d_periodic_fftpack, &
    sll_o_solve

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!> fft type use to do fft with fftpack library
type :: fftclass
   sll_real64, dimension(:), pointer :: coefc !< data for complex fft
   sll_real64, dimension(:), pointer :: work  !< work data for fft
   sll_real64, dimension(:), pointer :: workc !< workc complex
   sll_real64, dimension(:), pointer :: coefd !< data for double fft
   sll_real64, dimension(:), pointer :: workd !< work data
   sll_real64, dimension(:), pointer :: coefcd!< data for complex fft
   sll_int32  :: n  !< number of samples in each sequence
end type fftclass

!> Object with data to sll_o_solve Poisson equation on 2d domain with
!> periodic boundary conditions
type :: sll_t_poisson_2d_periodic_fftpack
  sll_int32   :: nc_x  !< number of cells direction x
  sll_int32   :: nc_y  !< number of cells direction y
  sll_real64  :: x_min !< left corner direction x
  sll_real64  :: x_max !< right corner direction x
  sll_real64  :: y_min !< left corner direction y
  sll_real64  :: y_max !< right corner direction y
  sll_real64  :: dx    !< step size direction x
  sll_real64  :: dy    !< step size direction y
  sll_comp64, dimension(:,:), pointer :: rhst !< rhs fft
  sll_comp64, dimension(:,:), pointer :: ext  !< x electric field fft
  sll_comp64, dimension(:,:), pointer :: eyt  !< y electric field fft
  sll_real64, dimension(:,:), pointer :: kx   !< wave number x
  sll_real64, dimension(:,:), pointer :: ky   !< wave number y
  sll_real64, dimension(:,:), pointer :: k2   !< \f$ k_x^2+k_y^2 \f$
  type(fftclass)                      :: fftx !< fft plan in direction x
  type(fftclass)                      :: ffty !< fft plan in direction y
end type sll_t_poisson_2d_periodic_fftpack

!> Create a sll_o_new poisson solver on 1d mesh
interface sll_o_new
  module procedure new_poisson_2d_periodic_fftpack
end interface


!> sll_o_initialize
interface sll_o_initialize
  module procedure initialize_poisson_2d_periodic_fftpack
end interface

!> sll_o_solve
interface sll_o_solve
   module procedure solve_potential_poisson_2d_periodic_fftpack
   module procedure solve_e_fields_poisson_2d_periodic_fftpack
end interface

interface delete
   module procedure free_poisson_2d_periodic_fftpack
end interface

!PN DEFINED BUT NOT USED
!PN interface fft
!PN    module procedure doubfft, doubcfft
!PN end interface
!PN interface fftinv
!PN    module procedure doubfftinv,  doubcfftinv
!PN end interface


contains

  !> Create a sll_o_new solver
  !> @return
  function new_poisson_2d_periodic_fftpack(&
    x_min, &
    x_max, &
    nc_x, &
    y_min, &
    y_max, &
    nc_y, &
    error) &
    result(self)

   type(sll_t_poisson_2d_periodic_fftpack),pointer :: self   !< self object
   sll_int32,  intent(in)    :: nc_x   !< number of cells direction x
   sll_int32,  intent(in)    :: nc_y   !< number of cells direction y
   sll_real64, intent(in)    :: x_min  !< left corner direction x
   sll_real64, intent(in)    :: x_max  !< right corner direction x
   sll_real64, intent(in)    :: y_min  !< left corner direction y
   sll_real64, intent(in)    :: y_max  !< right corner direction y
   sll_int32,  intent(out)   :: error  !< error code

   SLL_ALLOCATE(self, error)
   call initialize_poisson_2d_periodic_fftpack( &
           self, x_min, x_max, nc_x, y_min, y_max, nc_y, error )

  end function new_poisson_2d_periodic_fftpack 


!> delete sll_t_poisson_2d_periodic_fftpack
subroutine free_poisson_2d_periodic_fftpack( self, error )
   type(sll_t_poisson_2d_periodic_fftpack) :: self   !< self object
   sll_int32,  intent(out)   :: error  !< error code

   SLL_DEALLOCATE(self%rhst, error)
   SLL_DEALLOCATE(self%ext, error)
   SLL_DEALLOCATE(self%eyt, error)
   SLL_DEALLOCATE(self%kx, error)
   SLL_DEALLOCATE(self%ky, error)
   SLL_DEALLOCATE(self%k2, error)
   
   SLL_DEALLOCATE(self%fftx%coefc, error)
   SLL_DEALLOCATE(self%fftx%work, error)
   SLL_DEALLOCATE(self%fftx%workc, error)
   SLL_DEALLOCATE(self%fftx%coefd, error)
   SLL_DEALLOCATE(self%fftx%workd, error)
   SLL_DEALLOCATE(self%fftx%coefcd, error)

   SLL_DEALLOCATE(self%ffty%coefc, error)
   SLL_DEALLOCATE(self%ffty%work, error)
   SLL_DEALLOCATE(self%ffty%workc, error)
   SLL_DEALLOCATE(self%ffty%coefd, error)
   SLL_DEALLOCATE(self%ffty%workd, error)
   SLL_DEALLOCATE(self%ffty%coefcd, error)

      
end subroutine free_poisson_2d_periodic_fftpack

!> Create an object to sll_o_solve Poisson equation on 2D mesh with periodic
!> boundary conditions:
subroutine initialize_poisson_2d_periodic_fftpack( &
           self, x_min, x_max, nc_x, y_min, y_max, nc_y, error )

   type(sll_t_poisson_2d_periodic_fftpack) :: self   !< self object
   sll_int32,  intent(in)    :: nc_x   !< number of cells direction x
   sll_int32,  intent(in)    :: nc_y   !< number of cells direction y
   sll_real64, intent(in)    :: x_min  !< left corner direction x
   sll_real64, intent(in)    :: x_max  !< right corner direction x
   sll_real64, intent(in)    :: y_min  !< left corner direction y
   sll_real64, intent(in)    :: y_max  !< right corner direction y
   sll_int32,  intent(out)   :: error  !< error code
   
   self%nc_x = nc_x
   self%nc_y = nc_y

   self%x_min = x_min
   self%x_max = x_max
   self%y_min = y_min
   self%y_max = y_max
   self%dx   = (x_max-x_min) / real(nc_x, f64)
   self%dy   = (y_max-y_min) / real(nc_y, f64)

   SLL_ALLOCATE(self%rhst(nc_y,nc_x/2+1), error)
   SLL_ALLOCATE(self%ext (nc_y,nc_x/2+1), error)
   SLL_ALLOCATE(self%eyt (nc_y,nc_x/2+1), error)
   SLL_ALLOCATE(self%kx  (nc_y,nc_x/2+1), error)
   SLL_ALLOCATE(self%ky  (nc_y,nc_x/2+1), error)
   SLL_ALLOCATE(self%k2  (nc_y,nc_x/2+1), error)

#ifdef DEBUG
   print*, " FFTPACK version of poisson 2d periodic solver "
#endif

   call initdfft(self%fftx, nc_x)
   call initcfft(self%ffty, nc_y)

   call wave_number_vectors(self)

end subroutine initialize_poisson_2d_periodic_fftpack

!> sll_o_solve Poisson equation on 2D mesh with periodic boundary conditions. 
!> return potential.
subroutine solve_potential_poisson_2d_periodic_fftpack(self,sol,rhs)

   type(sll_t_poisson_2d_periodic_fftpack)       :: self !< self object
   sll_real64, dimension(:,:), intent(in)  :: rhs  !< charge density
   sll_real64, dimension(:,:), intent(out) :: sol  !< electric potential
   sll_int32                               :: nc_x !< number of cells direction x
   sll_int32                               :: nc_y !< number of cells direction y
   sll_int32                               :: i, j

   nc_x = self%nc_x
   nc_y = self%nc_y

   sol(1:nc_x,1:nc_y) = rhs(1:nc_x,1:nc_y)
   do j=1,nc_y
      call dfftf(nc_x, sol(1:nc_x,j), self%fftx%coefd)
   end do

   call transpose_r2c(sol(1:nc_x,1:nc_y), self%rhst)

   do i=1,nc_x/2+1
      call zfftf( nc_y, self%rhst(:,i), self%ffty%coefcd)
   end do

   self%rhst = self%rhst / self%k2

   do i=1,nc_x/2+1
      call zfftb( nc_y, self%rhst(:,i), self%ffty%coefcd )
   end do

   call transpose_c2r(self%rhst, sol(1:nc_x,1:nc_y))

   do j=1,nc_y
      call dfftb( nc_x, sol(1:nc_x,j),  self%fftx%coefd )
   end do

   sol(1:nc_x,1:nc_y) = sol(1:nc_x,1:nc_y) / real(nc_x*nc_y, f64)     ! normalize FFTs

   if (size(sol,1) == nc_x+1) sol(nc_x+1,:) = sol(1,:)
   if (size(sol,2) == nc_y+1) sol(:,nc_y+1) = sol(:,1)

end subroutine solve_potential_poisson_2d_periodic_fftpack

!> sll_o_solve Poisson equation on 2D mesh with periodic boundary conditions. 
!> return electric fields.
subroutine solve_e_fields_poisson_2d_periodic_fftpack(self,field_x,field_y,rhs,nrj)
! THIS routine changes the RHS despite its declaration as intent(in)
! It should be fixed !!!
   type(sll_t_poisson_2d_periodic_fftpack) :: self    !< self object
   sll_real64, dimension(:,:), intent(in)  :: rhs     !< charge density
   sll_real64, dimension(:,:), intent(out) :: field_x !< electric field direction x
   sll_real64, dimension(:,:), intent(out) :: field_y !< electric field direction y
   sll_int32                               :: nc_x    !< number of cells direction x
   sll_int32                               :: nc_y    !< number of cells direction y
   sll_int32                               :: i, j
   sll_real64, optional                    :: nrj     !< \f$ \sqrt{e_x^2+e_y^2} \f$

   nc_x = self%nc_x
   nc_y = self%nc_y

   self%rhst = cmplx(0.0_f64,0.0,kind=f64)
   self%ext  = cmplx(0.0_f64,0.0,kind=f64)
   self%eyt  = cmplx(0.0_f64,0.0,kind=f64)
   field_x   = 0.0_f64
   field_y   = 0.0_f64

   do j=1,nc_y
      call dfftf(nc_x, rhs(1:nc_x,j), self%fftx%coefd)
   end do

   call transpose_r2c(rhs(1:nc_x,1:nc_y), self%rhst)

   do i=1,nc_x/2+1
      call zfftf( nc_y, self%rhst(:,i), self%ffty%coefcd)
   end do

   self%ext(1,1) = (0.0_f64,0.0_f64)
   self%eyt(1,1) = (0.0_f64,0.0_f64)
   self%ext = -cmplx(0.0_f64,self%kx/self%k2,kind=f64)*self%rhst
   self%eyt = -cmplx(0.0_f64,self%ky/self%k2,kind=f64)*self%rhst

   do i=1,nc_x/2+1
      call zfftb( nc_y, self%ext(:,i), self%ffty%coefcd )
      call zfftb( nc_y, self%eyt(:,i), self%ffty%coefcd )
   end do

   call transpose_c2r(self%ext, field_x(1:nc_x,1:nc_y))
   call transpose_c2r(self%eyt, field_y(1:nc_x,1:nc_y))

   do j=1,nc_y
      call dfftb( nc_x, field_x(1:nc_x,j), self%fftx%coefd )
      call dfftb( nc_x, field_y(1:nc_x,j), self%fftx%coefd )
   end do

   field_x(1:nc_x,1:nc_y) = field_x(1:nc_x,1:nc_y) / (nc_x*nc_y)
   field_y(1:nc_x,1:nc_y) = field_y(1:nc_x,1:nc_y) / (nc_x*nc_y)

   if (size(field_x,1) == nc_x+1) field_x(nc_x+1,:) = field_x(1,:)
   if (size(field_x,2) == nc_y+1) field_x(:,nc_y+1) = field_x(:,1)
   if (size(field_y,1) == nc_x+1) field_y(nc_x+1,:) = field_y(1,:)
   if (size(field_y,2) == nc_y+1) field_y(:,nc_y+1) = field_y(:,1)

   if (present(nrj)) then 
      nrj=sum(field_x(1:nc_x,1:nc_y)*field_x(1:nc_x,1:nc_y) &
        +field_y(1:nc_x,1:nc_y)*field_y(1:nc_x,1:nc_y))*self%dx*self%dy
   end if

end subroutine solve_e_fields_poisson_2d_periodic_fftpack

subroutine wave_number_vectors(self)

   type(sll_t_poisson_2d_periodic_fftpack) :: self
   sll_int32  :: ik, jk
   sll_int32  :: nc_x, nc_y
   sll_real64 :: kx, ky, kx0, ky0
   
   nc_x = self%nc_x
   nc_y = self%nc_y
   
   kx0 = 2._f64*sll_p_pi/(self%x_max-self%x_min)
   ky0 = 2._f64*sll_p_pi/(self%y_max-self%y_min)
   
   do ik=1,nc_x/2+1
      kx = (ik-1)*kx0
      do jk = 1, nc_y/2
         ky = (jk-1)*ky0
         self%kx(jk,ik) = kx
         self%ky(jk,ik) = ky
      end do
      do jk = nc_y/2+1 , nc_y     
         ky = (jk-1-nc_y)*ky0
         self%kx(jk,ik) = kx
         self%ky(jk,ik) = ky
      end do
   end do
   self%kx(1,1) = 1.0_f64
   
   self%k2 = self%kx*self%kx+self%ky*self%ky

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
         real_array(2*i-1,j) = aimag(comp_array(j,i))
      end do
      real_array(n1,j) = real(comp_array(j,n1/2+1),kind=f64)
   end do

end subroutine transpose_c2r

subroutine initdfft(self,l)

   type(fftclass) :: self
   sll_int32 :: l 
   self%n = l 
   allocate(self%coefd(2*self%n+15))
   call dffti(self%n,self%coefd)

end subroutine initdfft

subroutine initcfft(self,l)

   type(fftclass) :: self
   sll_int32 :: l 
   self%n = l
   allocate(self%coefcd(4*self%n+15))
   call zffti(self%n,self%coefcd)

end subroutine initcfft

!PN DEFINED BUT NOT USED
!PN subroutine doubfft(self,array)
!PN 
!PN    type(fftclass) :: self
!PN    sll_real64, dimension(:,:) :: array
!PN    sll_int32 :: i
!PN 
!PN    do i=1, size(array,2)   ! number of 1d transforms
!PN       call dfftf( self%n, array(:,i), self%coefd)
!PN    end do
!PN 
!PN    array = array /self%n      ! normalize FFT
!PN 
!PN end subroutine doubfft
!PN 
!PN subroutine doubcfft(self,array)
!PN 
!PN    type(fftclass) :: self
!PN    sll_comp64, dimension(:,:) :: array
!PN    sll_int32 :: i
!PN 
!PN    do i=1, size(array,2)   ! number of 1d transforms
!PN       call zfftf( self%n, array(:,i), self%coefcd)
!PN    end do
!PN 
!PN    array = array /self%n      ! normalize FFT
!PN 
!PN end subroutine doubcfft
!PN 
!PN subroutine doubfftinv(self,array)
!PN 
!PN    type(fftclass) :: self
!PN    sll_real64, dimension(:,:) :: array
!PN    sll_int32 :: i
!PN 
!PN    do i=1, size(array,2)   ! number of 1d transforms
!PN       call dfftb( self%n, array(:,i),  self%coefd )
!PN    end do
!PN 
!PN end subroutine doubfftinv
!PN 
!PN subroutine doubcfftinv(self,array)
!PN 
!PN    type(fftclass) :: self
!PN    sll_comp64, dimension(:,:) :: array
!PN    sll_int32 :: i
!PN 
!PN    do i=1, size(array,2)   ! number of 1d transforms
!PN       call zfftb( self%n, array(:,i),  self%coefcd )
!PN    end do
!PN 
!PN end subroutine doubcfftinv

end module sll_m_poisson_2d_periodic_fftpack

#endif /* DOXYGEN_SHOULD_SKIP_THIS */
