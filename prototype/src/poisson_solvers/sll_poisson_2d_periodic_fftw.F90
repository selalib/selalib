!**************************************************************
!  Copyright INRIA
!  Authors : 
!     CALVI project team
!  
!  This code SeLaLib (for Semi-Lagrangian-Library) 
!  is a parallel library for simulating the plasma turbulence 
!  in a tokamak.
!  
!  This software is governed by the CeCILL-B license 
!  under French law and abiding by the rules of distribution 
!  of free software.  You can  use, modify and redistribute 
!  the software under the terms of the CeCILL-B license as 
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info". 
!**************************************************************

!> @brief  
!> Implements the Poisson solver in 2D with periodic boundary conditions
!> @details
!> This module uses fftw library
module sll_poisson_2d_periodic_fftw

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_constants.h"
#include "sll_fftw.h"

use fftw3

implicit none

!> Create a new poisson solver on 1d mesh
interface new
  module procedure new_poisson_2d_periodic_fftw
end interface

!> Initialize the Poisson solver using fftw library
interface initialize
  module procedure initialize_poisson_2d_periodic_fftw
end interface

!> get the potential or electric fields
interface solve
   module procedure solve_potential_poisson_2d_periodic_fftw
   module procedure solve_e_fields_poisson_2d_periodic_fftw
end interface

!> Delete the Poisson solver object
interface delete
   module procedure delete_poisson_2d_periodic_fftw
end interface

!> derived type to solve the Poisson equation on 2d regular cartesian mesh 
!> with periodic boundary conditions on both sides
type, public :: poisson_2d_periodic_fftw

   sll_real64, dimension(:,:), pointer :: kx       !< wave number in x
   sll_real64, dimension(:,:), pointer :: ky       !< wave number in y
   sll_real64, dimension(:,:), pointer :: k2       !< \f[ k_x^2 + k_y^2 \f]
   sll_int32                           :: nc_x     !< cells number in x
   sll_int32                           :: nc_y     !< cells number in y
   sll_real64                          :: dx       !< x step size
   sll_real64                          :: dy       !< y step size
   fftw_plan                           :: fw       !< forward fftw plan
   fftw_plan                           :: bw       !< backward fftw plan
   fftw_comp                , pointer  :: rht(:,:) !< fft(rho)
   fftw_comp                , pointer  :: ext(:,:) !< fft(ex)
   fftw_comp                , pointer  :: eyt(:,:) !< fft(ey)
   fftw_int                            :: sz_fft   !< fft size
   fftw_plan                           :: p_rht    !< C array pointer
   fftw_plan                           :: p_ext    !< C array pointer
   fftw_plan                           :: p_eyt    !< C array pointer

end type poisson_2d_periodic_fftw

public initialize, new, solve, delete

contains


  !> Create a new solver
  !> @return
  function new_poisson_2d_periodic_fftw(&
    x_min, &
    x_max, &
    nc_x, &
    y_min, &
    y_max, &
    nc_y, &
    error) &
    result(this)
   type(poisson_2d_periodic_fftw),pointer :: this   !< self object
   sll_int32,  intent(in)    :: nc_x   !< number of cells direction x
   sll_int32,  intent(in)    :: nc_y   !< number of cells direction y
   sll_real64, intent(in)    :: x_min  !< left corner direction x
   sll_real64, intent(in)    :: x_max  !< right corner direction x
   sll_real64, intent(in)    :: y_min  !< left corner direction y
   sll_real64, intent(in)    :: y_max  !< right corner direction y
   sll_int32,  intent(out)   :: error  !< error code

   SLL_ALLOCATE(this, error)
   call initialize_poisson_2d_periodic_fftw( &
           this, x_min, x_max, nc_x, y_min, y_max, nc_y, error )

  end function new_poisson_2d_periodic_fftw 





!> Initialize the Poisson solver
subroutine initialize_poisson_2d_periodic_fftw(self, &
                      x_min, x_max, nc_x, &
                      y_min, y_max, nc_y, error )

   type(poisson_2d_periodic_fftw) :: self   !< Self data object
   sll_real64, intent(in)    :: x_min  !< left corner direction x
   sll_real64, intent(in)    :: x_max  !< right corner direction x
   sll_real64, intent(in)    :: y_min  !< left corner direction y
   sll_real64, intent(in)    :: y_max  !< right corner direction y
   sll_int32                 :: error  !< error code
   sll_int32                 :: nc_x   !< number of cells direction x
   sll_int32                 :: nc_y   !< number of cells direction y
   sll_int32                 :: ik, jk
   sll_real64                :: kx1, kx0, ky0

   sll_real64, dimension(:,:), allocatable :: tmp

   self%nc_x = nc_x
   self%nc_y = nc_y

   self%dx   = (x_max-x_min) / nc_x
   self%dy   = (y_max-y_min) / nc_y


   SLL_ALLOCATE(tmp(1:nc_x,1:nc_y),error)

#ifdef FFTW_F2003

   self%sz_fft = int((nc_x/2+1)*nc_y,C_SIZE_T)

   self%p_rht = fftw_alloc_complex(self%sz_fft)
   call c_f_pointer(self%p_rht, self%rht, [nc_x/2+1,nc_y])

   self%p_ext = fftw_alloc_complex(self%sz_fft)
   call c_f_pointer(self%p_ext, self%ext, [nc_x/2+1,nc_y])

   self%p_eyt = fftw_alloc_complex(self%sz_fft)
   call c_f_pointer(self%p_eyt, self%eyt, [nc_x/2+1,nc_y])

   self%fw = fftw_plan_dft_r2c_2d(nc_y,nc_x,tmp,self%ext,FFTW_PATIENT)
   self%bw = fftw_plan_dft_c2r_2d(nc_y,nc_x,self%eyt,tmp,FFTW_PATIENT)

#else

   !call dfftw_init_threads(error)
   !if (error == 0) stop 'FFTW CAN''T USE THREADS'
   !call dfftw_plan_with_nthreads(nthreads)

   SLL_ALLOCATE(self%rht(1:nc_x/2+1,1:nc_y),error)
   SLL_ALLOCATE(self%ext(1:nc_x/2+1,1:nc_y),error)
   SLL_ALLOCATE(self%eyt(1:nc_x/2+1,1:nc_y),error)
   call dfftw_plan_dft_r2c_2d(self%fw,nc_x,nc_y,tmp,self%rht,FFTW_PATIENT)
   call dfftw_plan_dft_c2r_2d(self%bw,nc_x,nc_y,self%rht,tmp,FFTW_PATIENT)
   
#endif

   deallocate(tmp)

   SLL_ALLOCATE(self%kx(nc_x/2+1,nc_y), error)
   SLL_ALLOCATE(self%ky(nc_x/2+1,nc_y), error)
   SLL_ALLOCATE(self%k2(nc_x/2+1,nc_y), error)

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

end subroutine initialize_poisson_2d_periodic_fftw

!> Solve Poisson equation on 2D mesh with periodic boundary conditions. 
!> return potential.
subroutine solve_potential_poisson_2d_periodic_fftw(self, phi, rho)

   type(poisson_2d_periodic_fftw)   :: self !< self data object
   sll_real64, dimension(:,:), intent(in) :: rho  !< charge density
   sll_real64, dimension(:,:), intent(out)   :: phi  !< electric potential
   sll_int32                                 :: nc_x !< number of cells direction x
   sll_int32                                 :: nc_y !< number of cells direction y

   nc_x = self%nc_x
   nc_y = self%nc_y

   phi(1:nc_x,1:nc_y) = rho(1:nc_x,1:nc_y)
   call fftw_execute_dft_r2c(self%fw, phi(1:nc_x,1:nc_y), self%rht)

   self%rht = self%rht / self%k2

   call fftw_execute_dft_c2r(self%bw, self%rht, phi(1:nc_x,1:nc_y))

   nc_x = self%nc_x
   nc_y = self%nc_y

   phi = phi / (nc_x*nc_y)     ! normalize

   !Node centered case
   if(size(phi,1) == nc_x+1) phi(nc_x+1,:) = phi(1,:)
   if(size(phi,2) == nc_y+1) phi(:,nc_y+1) = phi(:,1)

end subroutine solve_potential_poisson_2d_periodic_fftw

!> Solve Poisson equation on 2D mesh with periodic boundary conditions. 
!> return electric fields.
subroutine solve_e_fields_poisson_2d_periodic_fftw(self,e_x,e_y,rho,nrj)

   type(poisson_2d_periodic_fftw),intent(inout)   :: self !< Self data object
   sll_real64, dimension(:,:), intent(in) :: rho  !< Charge density
   sll_real64, dimension(:,:), intent(out)   :: e_x  !< Electric field x
   sll_real64, dimension(:,:), intent(out)   :: e_y  !< Electric field y
   sll_real64, optional                      :: nrj  !< Energy 
   sll_int32  :: nc_x, nc_y
   sll_real64 :: dx, dy

   nc_x = self%nc_x
   nc_y = self%nc_y

   e_x(1:nc_x,1:nc_y) = rho(1:nc_x,1:nc_y)
   call fftw_execute_dft_r2c(self%fw, e_x(1:nc_x,1:nc_y), self%rht)

   self%ext(1,1) = 0.0_f64
   self%eyt(1,1) = 0.0_f64
   self%ext = -cmplx(0.0_f64,self%kx,kind=f64)*self%rht
   self%eyt = -cmplx(0.0_f64,self%ky,kind=f64)*self%rht

   call fftw_execute_dft_c2r(self%bw, self%ext, e_x(1:nc_x,1:nc_y))
   call fftw_execute_dft_c2r(self%bw, self%eyt, e_y(1:nc_x,1:nc_y))

   e_x = e_x / (nc_x*nc_y)
   e_y = e_y / (nc_x*nc_y)

   !Node centered case
   if (size(e_x,1) == nc_x+1) e_x(nc_x+1,:) = e_x(1,:)
   if (size(e_x,2) == nc_y+1) e_x(:,nc_y+1) = e_x(:,1)
   if (size(e_y,1) == nc_x+1) e_y(nc_x+1,:) = e_y(1,:)
   if (size(e_y,2) == nc_y+1) e_y(:,nc_y+1) = e_y(:,1)

   if (present(nrj)) then 
      dx = self%dx
      dy = self%dy
      nrj=sum(e_x(1:nc_x,1:nc_y)*e_x(1:nc_x,1:nc_y) &
        +e_y(1:nc_x,1:nc_y)*e_y(1:nc_x,1:nc_y))*dx*dy
      !if (nrj>1.e-30) then 
      !   nrj=0.5_f64*log(nrj)
      !else
      !   nrj=-10**9
      !endif
   end if

end subroutine solve_e_fields_poisson_2d_periodic_fftw

!> Delete the Poisson object
subroutine delete_poisson_2d_periodic_fftw(self)

   type(poisson_2d_periodic_fftw) :: self

#ifdef FFTW_F2003
   call fftw_free(self%p_rht)
   if (c_associated(self%p_ext)) call fftw_free(self%p_ext)
   if (c_associated(self%p_eyt)) call fftw_free(self%p_eyt)
#endif

   call fftw_destroy_plan(self%fw)
   call fftw_destroy_plan(self%bw)

   !if (nthreads > 1) then
      !call dfftw_cleanup_threads(error)
   !end if

end subroutine delete_poisson_2d_periodic_fftw

end module sll_poisson_2D_periodic_fftw
