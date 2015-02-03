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

!> @brief Spectral 1d advection
!> @details We are using FFTW to compute ffts. Boundary conditions
!> must be periodic.
module sll_module_advection_1d_spectral
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_fftw.h"
#include "sll_utilities.h"
#include "sll_constants.h"

use sll_boundary_condition_descriptors
use sll_module_advection_1d_base
use sll_module_interpolators_1d_base
use, intrinsic :: iso_c_binding
use fftw3

implicit none

private

type,extends(sll_advection_1d_base), public :: spectral_1d_advector
  
  sll_int32                          :: num_cells
  sll_int32                          :: num_points
  sll_real64                         :: eta_min
  sll_real64                         :: eta_max
  sll_int32                          :: bc
  sll_real64                         :: delta_eta
  sll_real64, dimension(:),  pointer :: d_dx
  sll_real64, dimension(:),  pointer :: kx
  fftw_plan                          :: fwx
  fftw_plan                          :: bwx
  fftw_plan                          :: p_tmp_x
  fftw_comp,  dimension(:),  pointer :: tmp_x

contains

  procedure, pass(adv) :: initialize 
  procedure, pass(adv) :: advect_1d 
  procedure, pass(adv) :: advect_1d_constant 
  procedure, pass(adv) :: delete

end type spectral_1d_advector

public ::  new_spectral_1d_advector

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function new_spectral_1d_advector( npts,     &
                                   eta_min,  &
                                   eta_max,  &
                                   bc ) result(adv)      

  type(spectral_1d_advector), pointer :: adv

  sll_int32,  intent(in)              :: npts
  sll_real64, intent(in), optional    :: eta_min
  sll_real64, intent(in), optional    :: eta_max
  sll_int32,  intent(in)              :: bc

  sll_int32                           :: ierr
  
  SLL_ALLOCATE(adv,ierr)
      
  call initialize( adv,        &
                   npts,       &
                   eta_min,    &
                   eta_max,    &
                   bc)    
  
end function new_spectral_1d_advector

subroutine initialize( adv, npts, eta_min, eta_max, bc)    

  class(spectral_1d_advector), intent(inout) :: adv
  sll_int32,                   intent(in)    :: npts
  sll_real64,        optional, intent(in)    :: eta_min
  sll_real64,        optional, intent(in)    :: eta_max
  sll_int32                                  :: bc

  sll_int32      :: i
  sll_int32      :: ierr
  sll_real64     :: kx0
  fftw_int       :: sz_tmp_x
  
  adv%num_points = npts
  adv%eta_min    = eta_min
  adv%eta_max    = eta_max
  adv%num_cells  = npts-1
  adv%delta_eta  = (eta_max-eta_min) / (npts-1)

  if( bc /= SLL_PERIODIC) then
    SLL_ERROR ('Boundary conditions for spectral advection must be periodic')
  else
    adv%bc = bc
  endif     

  FFTW_ALLOCATE(adv%tmp_x,adv%num_cells/2+1,sz_tmp_x,adv%p_tmp_x)
  SLL_CLEAR_ALLOCATE(adv%d_dx(1:adv%num_cells), ierr)

  NEW_FFTW_PLAN_R2C_1D(adv%fwx, adv%num_cells, adv%d_dx,  adv%tmp_x)
  NEW_FFTW_PLAN_C2R_1D(adv%bwx, adv%num_cells, adv%tmp_x, adv%d_dx)

  SLL_CLEAR_ALLOCATE(adv%kx(1:adv%num_cells/2+1), ierr)
   
  kx0 = 2._f64*sll_pi/(adv%num_cells*adv%delta_eta)

  adv%kx(1) = 1.0_f64
  do i=2,adv%num_cells/2+1
     adv%kx(i) = (i-1)*kx0
  end do
        
end subroutine initialize

subroutine advect_1d( adv, a, dt, input, output )

  class(spectral_1d_advector)           :: adv
  sll_real64, dimension(:), intent(in)  :: a
  sll_real64,               intent(in)  :: dt 
  sll_real64, dimension(:), intent(in)  :: input
  sll_real64, dimension(:), intent(out) :: output      

  sll_int32  :: nx
    
  nx = adv%num_cells

  adv%d_dx = input(1:nx)
  call fftw_execute_dft_r2c(adv%fwx, adv%d_dx, adv%tmp_x)

  !f = f^n exp(-i kx vx dt)

  adv%tmp_x = adv%tmp_x*cmplx(cos(adv%kx*a*dt),-sin(adv%kx*a*dt),kind=f64)

  call fftw_execute_dft_c2r(adv%bwx, adv%tmp_x, adv%d_dx)
  output(1:nx)= adv%d_dx / nx

  output(nx+1) = output(1)

end subroutine advect_1d

subroutine advect_1d_constant( adv, a, dt, input, output )

  class(spectral_1d_advector)              :: adv
  sll_real64,                  intent(in)  :: a
  sll_real64,                  intent(in)  :: dt 
  sll_real64, dimension(:),    intent(in)  :: input
  sll_real64, dimension(:),    intent(out) :: output      

  sll_int32  :: nx
    
  nx = adv%num_cells

  adv%d_dx = input(1:nx)
  call fftw_execute_dft_r2c(adv%fwx, adv%d_dx, adv%tmp_x)

  !f = f^n exp(-i kx vx dt)

  adv%tmp_x = adv%tmp_x*cmplx(cos(adv%kx*a*dt),-sin(adv%kx*a*dt),kind=f64)

  call fftw_execute_dft_c2r(adv%bwx, adv%tmp_x, adv%d_dx)
  output(1:nx)= adv%d_dx / nx

  output(nx+1) = output(1)

end subroutine advect_1d_constant

subroutine delete(adv)

  class(spectral_1d_advector), intent(inout) :: adv

#ifdef FFTW_F2003
  if (c_associated(adv%p_tmp_x)) call fftw_free(adv%p_tmp_x)
#endif

  call fftw_destroy_plan(adv%fwx)
  call fftw_destroy_plan(adv%bwx)

end subroutine delete

end module sll_module_advection_1d_spectral
