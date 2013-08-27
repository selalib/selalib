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

!> @author
!> Pierre Navaro
!> @brief
!> Implements the Poisson solver in 2D with periodic boundary conditions
!> @details
!> This module depends on:
!> - sll_memory
!> - sll_precision
!> - sll_assert 
!> - sll_constants
!> - sll_utilities
!>

module sll_poisson_2d_periodic

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_constants.h"

use, intrinsic :: iso_c_binding
implicit none
include 'fftw3.f03'

interface initialize
  module procedure initialize_poisson_2d_periodic_fftw
end interface

interface solve
   module procedure solve_potential_poisson_2d_periodic_fftw
   module procedure solve_e_fields_poisson_2d_periodic_fftw
end interface

interface delete
   module procedure free_poisson_2d_periodic_fftw
end interface

type, public :: poisson_2d_periodic
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

public initialize, solve, delete

contains

subroutine initialize_poisson_2d_periodic_fftw(self, &
                      x_min, x_max, nc_x, &
                      y_min, y_max, nc_y, error )

   type(poisson_2d_periodic) :: self
   sll_real64, intent(in) :: x_min, x_max, y_min, y_max
   sll_int32  :: error
   sll_int32  :: nc_x, nc_y
   sll_int32  :: ik, jk
   sll_real64 :: kx1, kx0, ky0
   sll_real64, dimension(:,:), allocatable :: tmp

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

   SLL_ALLOCATE(tmp(1:nc_x,1:nc_y),error)
   self%fw = fftw_plan_dft_r2c_2d(nc_y,nc_x,tmp,self%ext,FFTW_MEASURE)
   self%bw = fftw_plan_dft_c2r_2d(nc_y,nc_x,self%eyt,tmp,FFTW_MEASURE)
   deallocate(tmp)

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

   type(poisson_2d_periodic),intent(inout)  :: self
   sll_real64, dimension(:,:), intent(inout) :: rho
   sll_real64, dimension(:,:), intent(out)   :: phi
   sll_int32                                 :: nc_x, nc_y

   nc_x = self%nc_x
   nc_y = self%nc_y

   call fftw_execute_dft_r2c(self%fw, rho(1:nc_x,1:nc_y), self%rhot)

   self%rhot = self%rhot / self%k2

   call fftw_execute_dft_c2r(self%bw, self%rhot, phi(1:nc_x,1:nc_y))

   nc_x = self%nc_x
   nc_y = self%nc_y

   phi = phi / (nc_x*nc_y)     ! normalize

   phi(nc_x+1,:) = phi(1,:)
   phi(:,nc_y+1) = phi(:,1)

end subroutine solve_potential_poisson_2d_periodic_fftw

!> Solve Poisson equation on 2D mesh with periodic boundary conditions. 
!> return electric fields.
subroutine solve_e_fields_poisson_2d_periodic_fftw(self,e_x,e_y,rho,nrj)

   type(poisson_2d_periodic),intent(inout)  :: self
   sll_real64, dimension(:,:), intent(inout) :: rho
   sll_real64, dimension(:,:), intent(out)   :: e_x
   sll_real64, dimension(:,:), intent(out)   :: e_y
   sll_real64, optional                      :: nrj
   sll_int32  :: nc_x, nc_y
   sll_real64 :: dx, dy

   nc_x = self%nc_x
   nc_y = self%nc_y

   call fftw_execute_dft_r2c(self%fw, rho(1:nc_x,1:nc_y), self%rhot)

   self%ext(1,1) = 0.0_f64
   self%eyt(1,1) = 0.0_f64
   self%ext = -cmplx(0.0_f64,self%kx,kind=f64)*self%rhot
   self%eyt = -cmplx(0.0_f64,self%ky,kind=f64)*self%rhot

   call fftw_execute_dft_c2r(self%bw, self%ext, e_x(1:nc_x,1:nc_y))
   call fftw_execute_dft_c2r(self%bw, self%eyt, e_y(1:nc_x,1:nc_y))

   e_x = e_x / (nc_x*nc_y)
   e_y = e_y / (nc_x*nc_y)

   e_x(nc_x+1,:) = e_x(1,:)
   e_x(:,nc_y+1) = e_x(:,1)
   e_y(nc_x+1,:) = e_y(1,:)
   e_y(:,nc_y+1) = e_y(:,1)

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

end subroutine solve_e_fields_poisson_2d_periodic_fftw

subroutine free_poisson_2d_periodic_fftw(self)
type(poisson_2d_periodic) :: self
call fftw_free(self%p_rhot)
if (c_associated(self%p_ext)) call fftw_free(self%p_ext)
if (c_associated(self%p_eyt)) call fftw_free(self%p_eyt)
call fftw_destroy_plan(self%fw)
call fftw_destroy_plan(self%bw)
!if (nthreads > 1) then
   !call dfftw_cleanup_threads(error)
!end if

end subroutine

end module sll_poisson_2D_periodic
