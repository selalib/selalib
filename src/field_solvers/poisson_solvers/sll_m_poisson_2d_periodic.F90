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

!> @ingroup poisson_solvers
!> @brief
!> Regular cartesian two dimensional mesh with periodic bounday conditions.
!> @details
!> Numerical method uses Fast Fourier Transform and periodic
!> boundary conditions.
!> @snippet poisson_solvers/test_poisson_2d_fft.F90 example
module sll_m_poisson_2d_periodic
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_assert.h"

use sll_m_poisson_2d_base, only: sll_c_poisson_2d_base
use sll_m_constants, only: sll_p_pi
use sll_m_fft, only: sll_t_fft, &
  sll_s_fft_init_c2r_2d,        &
  sll_s_fft_init_r2c_2d,        &
  sll_s_fft_exec_c2r_2d,        &
  sll_s_fft_exec_r2c_2d,        &
  sll_s_fft_free
implicit none

private

public :: sll_f_new_poisson_2d_periodic, sll_t_poisson_2d_periodic, &
          sll_t_poisson_2d_periodic_fft, sll_o_solve, &
          sll_o_initialize

!> derived type to sll_o_solve the Poisson equation on 2d regular 
!> cartesian mesh with periodic boundary conditions on both sides
type :: sll_t_poisson_2d_periodic_fft

   private
   sll_real64, dimension(:,:), pointer :: kx       !< wave number in x
   sll_real64, dimension(:,:), pointer :: ky       !< wave number in y
   sll_real64, dimension(:,:), pointer :: k2       !< \f[ k_x^2 + k_y^2 \f]
   sll_int32                           :: nc_x     !< cells number in x
   sll_int32                           :: nc_y     !< cells number in y
   sll_real64                          :: dx       !< x step size
   sll_real64                          :: dy       !< y step size
   sll_comp64               , pointer  :: rht(:,:) !< fft(rho)
   sll_comp64               , pointer  :: exy(:,:) !< fft(ex and ey)
   type(sll_t_fft)                     :: fw       !< forward fft plan
   type(sll_t_fft)                     :: bw       !< backward fft plan
   type(sll_t_fft)                     :: p_rho    !< C array pointer
   type(sll_t_fft)                     :: p_exy    !< C array pointer
   type(sll_t_fft)                     :: p_tmp    !< C array pointer
   sll_real64               , pointer  :: tmp(:,:)

end type sll_t_poisson_2d_periodic_fft

interface sll_o_initialize
  module procedure initialize_poisson_2d_periodic_fft
end interface

interface sll_o_solve
   module procedure solve_potential_poisson_2d_periodic_fft
   module procedure solve_e_fields_poisson_2d_periodic_fft
end interface


type, extends(sll_c_poisson_2d_base) :: sll_t_poisson_2d_periodic

  type(sll_t_poisson_2d_periodic_fft), private, pointer :: solver

contains

  !> Create the Poisson solver
  procedure, public, pass(poisson) :: initialize => &
    initialize_poisson_2d_periodic
  !> Compute potential solving the Poisson equation
  procedure, public, pass(poisson) :: compute_phi_from_rho => &
    compute_phi_from_rho_2d_fft
  !> Compute electric fields solving the Poisson equation
  procedure, public, pass(poisson) :: compute_E_from_rho => &
    compute_E_from_rho_2d_fft
    
end type sll_t_poisson_2d_periodic

contains

  !> @returns a pointer to the derived type sll_t_poisson_2d_periodic.
  function sll_f_new_poisson_2d_periodic( &
    eta1_min, &
    eta1_max, &
    nc_eta1, &
    eta2_min, &
    eta2_max, &
    nc_eta2) &     
    result(poisson)
      
    type(sll_t_poisson_2d_periodic),pointer :: poisson
    sll_real64 :: eta1_min
    sll_real64 :: eta1_max
    sll_int32 :: nc_eta1
    sll_real64 :: eta2_min
    sll_real64 :: eta2_max
    sll_int32 :: nc_eta2
    sll_int32 :: ierr
      
    SLL_ALLOCATE(poisson,ierr)
    call initialize_poisson_2d_periodic( &
    poisson, &
    eta1_min, &
    eta1_max, &
    nc_eta1, &
    eta2_min, &
    eta2_max, &
    nc_eta2)     
    
  end function sll_f_new_poisson_2d_periodic
  
  subroutine initialize_poisson_2d_periodic( &
    poisson, &
    eta1_min, &
    eta1_max, &
    nc_eta1, &
    eta2_min, &
    eta2_max, &
    nc_eta2)     
    class(sll_t_poisson_2d_periodic) :: poisson
    sll_real64 :: eta1_min
    sll_real64 :: eta1_max
    sll_int32 :: nc_eta1
    sll_real64 :: eta2_min
    sll_real64 :: eta2_max
    sll_int32 :: nc_eta2
    sll_int32 :: ierr
    
    SLL_ALLOCATE(poisson%solver,ierr)
    
    call sll_o_initialize( &
      poisson%solver, &
      eta1_min, &
      eta1_max, &
      nc_eta1, &
      eta2_min, &
      eta2_max, &
      nc_eta2, &
      ierr) 

  end subroutine initialize_poisson_2d_periodic
  
  !> solves \f$ -\Delta phi(x,y) = rho (x,y) \f$
  subroutine compute_phi_from_rho_2d_fft( poisson, phi, rho )
    class(sll_t_poisson_2d_periodic), target :: poisson
    sll_real64,dimension(:,:),intent(in) :: rho
    sll_real64,dimension(:,:),intent(out) :: phi
    
    call sll_o_solve( poisson%solver, phi, rho)
    
  end subroutine compute_phi_from_rho_2d_fft

  !> @brief
  !> sll_o_solve Poisson equation to compute electric fields
  !> @details
  !> solves 
  !> \f[ 
  !> E(x,y) = -\nabla \phi(x,y) \\
  !> -\Delta \phi(x,y) = \rho(x,y)
  !> \f]
  subroutine compute_E_from_rho_2d_fft( poisson, E1, E2, rho )
    class(sll_t_poisson_2d_periodic) :: poisson
    sll_real64,dimension(:,:),intent(in) :: rho
    sll_real64,dimension(:,:),intent(out) :: E1
    sll_real64,dimension(:,:),intent(out) :: E2
      
    call sll_o_solve( poisson%solver, E1, E2, rho)
      
  end subroutine compute_E_from_rho_2d_fft

!> Create a sll_o_new solver
!> @return a pointer to the solver derived type
function new_poisson_2d_periodic_fft( &
   x_min,                              &
   x_max,                              &
   nc_x,                               &
   y_min,                              &
   y_max,                              &
   nc_y,                               &
   error)                              &
   result(self)

  type(sll_t_poisson_2d_periodic_fft),pointer :: self   !< self object
  sll_int32,  intent(in)    :: nc_x   !< number of cells direction x
  sll_int32,  intent(in)    :: nc_y   !< number of cells direction y
  sll_real64, intent(in)    :: x_min  !< left corner direction x
  sll_real64, intent(in)    :: x_max  !< right corner direction x
  sll_real64, intent(in)    :: y_min  !< left corner direction y
  sll_real64, intent(in)    :: y_max  !< right corner direction y
  sll_int32,  intent(out)   :: error  !< error code

  SLL_ALLOCATE(self, error)
  call initialize_poisson_2d_periodic_fft( &
          self, x_min, x_max, nc_x, y_min, y_max, nc_y, error )

end function new_poisson_2d_periodic_fft 

!> sll_o_initialize the Poisson solver
subroutine initialize_poisson_2d_periodic_fft(self, &
  x_min,                                             &
  x_max,                                             &
  nc_x,                                              &
  y_min,                                             &
  y_max,                                             &
  nc_y,                                              &
  error )

  type(sll_t_poisson_2d_periodic_fft) :: self   !< Self data object
  sll_real64, intent(in)              :: x_min  !< left corner direction x
  sll_real64, intent(in)              :: x_max  !< right corner direction x
  sll_real64, intent(in)              :: y_min  !< left corner direction y
  sll_real64, intent(in)              :: y_max  !< right corner direction y
  sll_int32                           :: error  !< error code
  sll_int32                           :: nc_x   !< number of cells direction x
  sll_int32                           :: nc_y   !< number of cells direction y
  sll_int32                           :: ik, jk
  sll_real64                          :: kx1, kx0, ky0

  self%nc_x = nc_x
  self%nc_y = nc_y

  self%dx   = (x_max-x_min) / real( nc_x, f64)
  self%dy   = (y_max-y_min) / real( nc_y, f64)

#ifdef DEBUG
  print*, " FFT version of poisson 2d periodic solver "
#endif

  SLL_ALLOCATE(self%rht(1:nc_x/2+1,1:nc_y),error)
  SLL_ALLOCATE(self%exy(1:nc_x/2+1,1:nc_y),error)
  SLL_CLEAR_ALLOCATE(self%tmp(1:nc_x,1:nc_y),error)
  call sll_s_fft_init_r2c_2d(self%fw,nc_x,nc_y,self%tmp,self%rht)
  call sll_s_fft_init_c2r_2d(self%bw,nc_x,nc_y,self%rht,self%tmp)
   
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

end subroutine initialize_poisson_2d_periodic_fft

!> sll_o_solve Poisson equation on 2D mesh with periodic boundary conditions. 
!> return potential.
subroutine solve_potential_poisson_2d_periodic_fft(self, phi, rho)

  type(sll_t_poisson_2d_periodic_fft)          :: self !< self data object
  sll_real64, dimension(:,:), intent(in)  :: rho  !< charge density
  sll_real64, dimension(:,:), intent(out) :: phi  !< electric potential
  sll_int32                               :: nc_x !< number of cells direction x
  sll_int32                               :: nc_y !< number of cells direction y

  nc_x = self%nc_x
  nc_y = self%nc_y

  self%tmp = rho(1:nc_x,1:nc_y)
  call sll_s_fft_exec_r2c_2d(self%fw, self%tmp, self%rht)

  self%rht = self%rht / self%k2

  call sll_s_fft_exec_c2r_2d(self%bw, self%rht, self%tmp)

  phi(1:nc_x,1:nc_y) = self%tmp / real(nc_x*nc_y, f64)   
  
  !Node centered case
  if(size(phi,1) == nc_x+1) phi(nc_x+1,:) = phi(1,:)
  if(size(phi,2) == nc_y+1) phi(:,nc_y+1) = phi(:,1)

end subroutine solve_potential_poisson_2d_periodic_fft

!> sll_o_solve Poisson equation on 2D mesh with periodic boundary conditions. 
!> return electric fields.
subroutine solve_e_fields_poisson_2d_periodic_fft(self,e_x,e_y,rho,nrj)

  type(sll_t_poisson_2d_periodic_fft),intent(inout) :: self !< Self data object
  sll_real64, dimension(:,:),    intent(in)    :: rho  !< Charge density
  sll_real64, dimension(:,:),    intent(out)   :: e_x  !< Electric field x
  sll_real64, dimension(:,:),    intent(out)   :: e_y  !< Electric field y
  sll_real64, optional                         :: nrj  !< Energy 

  sll_int32  :: nc_x, nc_y
  sll_real64 :: dx, dy

  nc_x = self%nc_x
  nc_y = self%nc_y

  self%tmp = rho(1:nc_x,1:nc_y)
  call sll_s_fft_exec_r2c_2d(self%fw, self%tmp, self%rht)

  self%exy(1,1) = (0.0_f64,0.0_f64)
  self%exy = -cmplx(0.0_f64,self%kx,kind=f64)*self%rht
  call sll_s_fft_exec_c2r_2d(self%bw, self%exy, self%tmp)
  e_x(1:nc_x,1:nc_y) = self%tmp / real(nc_x*nc_y, f64)

  self%exy(1,1) = (0.0_f64,0.0_f64)
  self%exy = -cmplx(0.0_f64,self%ky,kind=f64)*self%rht
  call sll_s_fft_exec_c2r_2d(self%bw, self%exy, self%tmp)

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

end subroutine solve_e_fields_poisson_2d_periodic_fft

!> Delete the Poisson object
subroutine delete_poisson_2d_periodic_fft(self)

  type(sll_t_poisson_2d_periodic_fft) :: self

  call sll_s_fft_free(self%fw)
  call sll_s_fft_free(self%bw)

end subroutine delete_poisson_2d_periodic_fft

end module sll_m_poisson_2d_periodic
