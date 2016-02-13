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

#ifndef DOXYGEN_SHOULD_SKIP_THIS

!> @ingroup poisson_solvers
!> @brief  
!> Implements the Poisson solver in 2D with periodic boundary conditions
!> @details
!> This module uses fftw library
module sll_m_poisson_2d_periodic_fftw

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_fftw.h"

  use iso_c_binding, only: &
    c_associated, &
    c_double, &
    c_double_complex, &
    c_f_pointer, &
    c_ptr, &
    c_size_t

  use sll_m_constants, only: &
    sll_p_pi

#ifdef FFTW_F2003
  use sll_m_fftw3, only: &
    fftw_alloc_complex, &
    fftw_alloc_real, &
    fftw_destroy_plan, &
    fftw_estimate, &
    fftw_execute_dft_c2r, &
    fftw_execute_dft_r2c, &
    fftw_free, &
    fftw_plan_dft_c2r_2d, &
    fftw_plan_dft_r2c_2d
#else
  use sll_m_fftw3, only : &
       fftw_estimate
#endif

  implicit none

  public :: &
    sll_o_initialize, &
    sll_o_new, &
    sll_t_poisson_2d_periodic_fftw, &
    sll_o_solve

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!> Create a sll_o_new poisson solver on 1d mesh
interface sll_o_new
  module procedure new_poisson_2d_periodic_fftw
end interface

!> sll_o_initialize the Poisson solver using fftw library
interface sll_o_initialize
  module procedure initialize_poisson_2d_periodic_fftw
end interface

!> get the potential or electric fields
interface sll_o_solve
   module procedure solve_potential_poisson_2d_periodic_fftw
   module procedure solve_e_fields_poisson_2d_periodic_fftw
end interface

!> Delete the Poisson solver object
interface delete
   module procedure delete_poisson_2d_periodic_fftw
end interface

!> derived type to sll_o_solve the Poisson equation on 2d regular cartesian mesh 
!> with periodic boundary conditions on both sides
type :: sll_t_poisson_2d_periodic_fftw

   private
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
   fftw_comp                , pointer  :: exy(:,:) !< fft(ex and ey)
   fftw_plan                           :: p_rho    !< C array pointer
   fftw_plan                           :: p_exy    !< C array pointer
   fftw_plan                           :: p_tmp    !< C array pointer
   fftw_real                , pointer  :: tmp(:,:)

end type sll_t_poisson_2d_periodic_fftw


contains


  !> Create a sll_o_new solver
  !> @return a pointer to the solver derived type
function new_poisson_2d_periodic_fftw( &
   x_min,                              &
   x_max,                              &
   nc_x,                               &
   y_min,                              &
   y_max,                              &
   nc_y,                               &
   error)                              &
   result(self)

  type(sll_t_poisson_2d_periodic_fftw),pointer :: self   !< self object
  sll_int32,  intent(in)    :: nc_x   !< number of cells direction x
  sll_int32,  intent(in)    :: nc_y   !< number of cells direction y
  sll_real64, intent(in)    :: x_min  !< left corner direction x
  sll_real64, intent(in)    :: x_max  !< right corner direction x
  sll_real64, intent(in)    :: y_min  !< left corner direction y
  sll_real64, intent(in)    :: y_max  !< right corner direction y
  sll_int32,  intent(out)   :: error  !< error code

  SLL_ALLOCATE(self, error)
  call initialize_poisson_2d_periodic_fftw( &
          self, x_min, x_max, nc_x, y_min, y_max, nc_y, error )

end function new_poisson_2d_periodic_fftw 


!> sll_o_initialize the Poisson solver
subroutine initialize_poisson_2d_periodic_fftw(self, &
  x_min,                                             &
  x_max,                                             &
  nc_x,                                              &
  y_min,                                             &
  y_max,                                             &
  nc_y,                                              &
  error )

  type(sll_t_poisson_2d_periodic_fftw) :: self   !< Self data object
  sll_real64, intent(in)         :: x_min  !< left corner direction x
  sll_real64, intent(in)         :: x_max  !< right corner direction x
  sll_real64, intent(in)         :: y_min  !< left corner direction y
  sll_real64, intent(in)         :: y_max  !< right corner direction y
  sll_int32                      :: error  !< error code
  sll_int32                      :: nc_x   !< number of cells direction x
  sll_int32                      :: nc_y   !< number of cells direction y
  sll_int32                      :: ik, jk
  sll_real64                     :: kx1, kx0, ky0

  fftw_int                       :: sz


  self%nc_x = nc_x
  self%nc_y = nc_y

  self%dx   = (x_max-x_min) / real( nc_x, f64)
  self%dy   = (y_max-y_min) / real( nc_y, f64)

#ifdef DEBUG
  print*, " FFTW version of poisson 2d periodic solver "
#endif

#ifdef FFTW_F2003

  sz = int((nc_x/2+1)*nc_y,C_SIZE_T)

  self%p_rho = fftw_alloc_complex(sz)
  call c_f_pointer(self%p_rho, self%rht, [nc_x/2+1,nc_y])

  self%p_exy = fftw_alloc_complex(sz)
  call c_f_pointer(self%p_exy, self%exy, [nc_x/2+1,nc_y])

  self%p_tmp = fftw_alloc_real(int(nc_x*nc_y,C_SIZE_T))
  call c_f_pointer(self%p_tmp, self%tmp, [nc_x,nc_y])

  self%fw = fftw_plan_dft_r2c_2d(nc_y,nc_x,self%tmp,self%exy,FFTW_ESTIMATE)
  self%bw = fftw_plan_dft_c2r_2d(nc_y,nc_x,self%exy,self%tmp,FFTW_ESTIMATE)

#else

  SLL_ALLOCATE(self%rht(1:nc_x/2+1,1:nc_y),error)
  SLL_ALLOCATE(self%exy(1:nc_x/2+1,1:nc_y),error)
  SLL_ALLOCATE(self%tmp(nc_x,nc_y),error)
  call dfftw_plan_dft_r2c_2d(self%fw,nc_x,nc_y,self%tmp,self%rht,FFTW_ESTIMATE)
  call dfftw_plan_dft_c2r_2d(self%bw,nc_x,nc_y,self%rht,self%tmp,FFTW_ESTIMATE)
   
#endif

  SLL_ALLOCATE(self%kx(nc_x/2+1,nc_y), error)
  SLL_ALLOCATE(self%ky(nc_x/2+1,nc_y), error)
  SLL_ALLOCATE(self%k2(nc_x/2+1,nc_y), error)

  kx0 = 2._f64*sll_p_pi/(x_max-x_min)
  ky0 = 2._f64*sll_p_pi/(y_max-y_min)
  
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

!> sll_o_solve Poisson equation on 2D mesh with periodic boundary conditions. 
!> return potential.
subroutine solve_potential_poisson_2d_periodic_fftw(self, phi, rho)

  type(sll_t_poisson_2d_periodic_fftw)          :: self !< self data object
  sll_real64, dimension(:,:), intent(in)  :: rho  !< charge density
  sll_real64, dimension(:,:), intent(out) :: phi  !< electric potential
  sll_int32                               :: nc_x !< number of cells direction x
  sll_int32                               :: nc_y !< number of cells direction y

  nc_x = self%nc_x
  nc_y = self%nc_y

  self%tmp = rho(1:nc_x,1:nc_y)
  call fftw_execute_dft_r2c(self%fw, self%tmp, self%rht)

  self%rht = self%rht / self%k2

  call fftw_execute_dft_c2r(self%bw, self%rht, self%tmp)

  phi(1:nc_x,1:nc_y) = self%tmp / real(nc_x*nc_y, f64)   
  
  !Node centered case
  if(size(phi,1) == nc_x+1) phi(nc_x+1,:) = phi(1,:)
  if(size(phi,2) == nc_y+1) phi(:,nc_y+1) = phi(:,1)

end subroutine solve_potential_poisson_2d_periodic_fftw

!> sll_o_solve Poisson equation on 2D mesh with periodic boundary conditions. 
!> return electric fields.
subroutine solve_e_fields_poisson_2d_periodic_fftw(self,e_x,e_y,rho,nrj)

  type(sll_t_poisson_2d_periodic_fftw),intent(inout) :: self !< Self data object
  sll_real64, dimension(:,:),    intent(in)    :: rho  !< Charge density
  sll_real64, dimension(:,:),    intent(out)   :: e_x  !< Electric field x
  sll_real64, dimension(:,:),    intent(out)   :: e_y  !< Electric field y
  sll_real64, optional                         :: nrj  !< Energy 

  sll_int32  :: nc_x, nc_y
  sll_real64 :: dx, dy

  nc_x = self%nc_x
  nc_y = self%nc_y

  self%tmp = rho(1:nc_x,1:nc_y)
  call fftw_execute_dft_r2c(self%fw, self%tmp, self%rht)

  self%exy(1,1) = (0.0_f64,0.0_f64)
  self%exy = -cmplx(0.0_f64,self%kx,kind=f64)*self%rht
  call fftw_execute_dft_c2r(self%bw, self%exy, self%tmp)
  e_x(1:nc_x,1:nc_y) = self%tmp / real(nc_x*nc_y, f64)

  self%exy(1,1) = (0.0_f64,0.0_f64)
  self%exy = -cmplx(0.0_f64,self%ky,kind=f64)*self%rht
  call fftw_execute_dft_c2r(self%bw, self%exy, self%tmp)

  e_y(1:nc_x,1:nc_y) = self%tmp / real(nc_x*nc_y, f64)

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
  end if

end subroutine solve_e_fields_poisson_2d_periodic_fftw

!> Delete the Poisson object
subroutine delete_poisson_2d_periodic_fftw(self)

  type(sll_t_poisson_2d_periodic_fftw) :: self

#ifdef FFTW_F2003
  call fftw_free(self%p_rho)
  if (c_associated(self%p_exy)) call fftw_free(self%p_exy)
  if (c_associated(self%p_tmp)) call fftw_free(self%p_tmp)
  if (c_associated(self%p_rho)) call fftw_free(self%p_rho)
#endif

  call fftw_destroy_plan(self%fw)
  call fftw_destroy_plan(self%bw)

end subroutine delete_poisson_2d_periodic_fftw

end module sll_m_poisson_2d_periodic_fftw
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
