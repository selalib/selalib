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
#include "sll_utilities.h"
#include "sll_constants.h"

use sll_boundary_condition_descriptors
use sll_module_advection_1d_base
use sll_module_interpolators_1d_base
use sll_fft

implicit none

private

type,extends(sll_advection_1d_base), public :: spectral_1d_advector
  
  sll_int32                         :: num_cells
  sll_real64                        :: eta_min
  sll_real64                        :: eta_max
  sll_real64                        :: delta_eta
  sll_real64, dimension(:), pointer :: d_dx
  sll_real64, dimension(:), pointer :: kx
  type(sll_fft_plan),       pointer :: fwx
  type(sll_fft_plan),       pointer :: bwx
  sll_comp64, dimension(:), pointer :: tmp_x

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

function new_spectral_1d_advector( num_cells, &
                                   eta_min,   &
                                   eta_max    &
                                   ) result(adv)      

  type(spectral_1d_advector), pointer :: adv

  sll_int32,  intent(in)              :: num_cells
  sll_real64, intent(in), optional    :: eta_min
  sll_real64, intent(in), optional    :: eta_max

  sll_int32                           :: error
  
  SLL_ALLOCATE(adv,error)
      
  call initialize( adv, num_cells, eta_min, eta_max)    
  
end function new_spectral_1d_advector

subroutine initialize( adv, num_cells, eta_min, eta_max)    

  class(spectral_1d_advector), intent(inout) :: adv
  sll_int32,                   intent(in)    :: num_cells
  sll_real64,        optional, intent(in)    :: eta_min
  sll_real64,        optional, intent(in)    :: eta_max

  sll_int32      :: i
  sll_int32      :: error
  sll_real64     :: kx0
  
  adv%eta_min    = eta_min
  adv%eta_max    = eta_max
  adv%num_cells  = num_cells
  adv%delta_eta  = (eta_max-eta_min) / num_cells

  SLL_CLEAR_ALLOCATE(adv%d_dx(1:num_cells), error)
  SLL_CLEAR_ALLOCATE(adv%tmp_x(1:num_cells/2+1), error)

  adv%fwx => fft_new_plan(num_cells, adv%d_dx,  adv%tmp_x)
  adv%bwx => fft_new_plan(num_cells, adv%tmp_x, adv%d_dx)

  SLL_CLEAR_ALLOCATE(adv%kx(1:num_cells/2+1), error)
   
  kx0 = 2._f64*sll_pi/(eta_max-eta_min)

  adv%kx(1) = 1.0_f64
  do i=1,num_cells/2-1
     adv%kx(i+1)       = i * kx0
  end do
  adv%kx(num_cells/2+1) = sll_pi / adv%delta_eta

        
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
  call fft_apply_plan(adv%fwx, adv%d_dx, adv%tmp_x)

  !f = f^n exp(-i kx vx dt)

  adv%tmp_x = adv%tmp_x*cmplx(cos(adv%kx*a*dt),sin(adv%kx*a*dt),kind=f64)

  call fft_apply_plan(adv%bwx, adv%tmp_x, adv%d_dx)
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
  call fft_apply_plan(adv%fwx, adv%d_dx, adv%tmp_x)

  !f = f^n exp(-i kx vx dt)

  adv%tmp_x = adv%tmp_x* exp(-cmplx(0.0_f64,adv%kx*a*dt,kind=f64))



  call fft_apply_plan(adv%bwx, adv%tmp_x, adv%d_dx)
  output(1:nx)= adv%d_dx / nx

  output(nx+1) = output(1)

end subroutine advect_1d_constant

subroutine delete(adv)

  class(spectral_1d_advector), intent(inout) :: adv

  call fft_delete_plan(adv%fwx)
  call fft_delete_plan(adv%bwx)

end subroutine delete

end module sll_module_advection_1d_spectral
