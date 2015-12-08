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
module sll_m_advection_1d_ampere
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_errors.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_advection_1d_base, only: &
    sll_advection_1d_base

  use sll_m_constants, only: &
    sll_pi

  use sll_m_fft, only: &
    fft_apply_plan_c2r_1d, &
    fft_apply_plan_r2c_1d, &
    fft_delete_plan, &
    fft_new_plan_c2r_1d, &
    fft_new_plan_r2c_1d, &
    sll_fft_plan

  implicit none

  public :: &
    ampere_1d_advector_ptr, &
    new_ampere_1d_advector

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

type,extends(sll_advection_1d_base) :: ampere_1d_advector
  
  sll_int32                         :: nc_eta1
  sll_real64                        :: eta1_min
  sll_real64                        :: eta1_max
  sll_real64                        :: delta_eta1
  sll_real64, dimension(:), pointer :: d_dx
  sll_real64, dimension(:), pointer :: kx
  type(sll_fft_plan),       pointer :: fwx
  type(sll_fft_plan),       pointer :: bwx
  sll_comp64, dimension(:), pointer :: fk
  sll_comp64, dimension(:), pointer :: r0
  sll_comp64, dimension(:), pointer :: r1
  sll_comp64, dimension(:), pointer :: ek

contains

  procedure, pass(adv) :: initialize 
  procedure, pass(adv) :: advect_1d 
  procedure, pass(adv) :: advect_1d_constant 
  procedure, pass(adv) :: delete

end type ampere_1d_advector

type :: ampere_1d_advector_ptr 
  class(ampere_1d_advector), pointer :: ptr
end type ampere_1d_advector_ptr


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function new_ampere_1d_advector( nc_eta1,  &
                                 eta1_min, &
                                 eta1_max  ) result(adv)      

  type(ampere_1d_advector), pointer :: adv

  sll_int32,  intent(in) :: nc_eta1
  sll_real64, intent(in) :: eta1_min
  sll_real64, intent(in) :: eta1_max

  sll_int32              :: error
  
  SLL_ALLOCATE(adv,error)
      
  call initialize( adv, nc_eta1, eta1_min, eta1_max  )
  
end function new_ampere_1d_advector

subroutine initialize( adv, nc_eta1, eta1_min, eta1_max )

  class(ampere_1d_advector), intent(inout) :: adv
  sll_int32,                 intent(in)    :: nc_eta1
  sll_real64,                intent(in)    :: eta1_min
  sll_real64,                intent(in)    :: eta1_max

  sll_int32     :: i
  sll_int32     :: error
  sll_real64    :: kx0
  
  adv%eta1_min   = eta1_min
  adv%eta1_max   = eta1_max
  adv%nc_eta1    = nc_eta1
  adv%delta_eta1 = (eta1_max-eta1_min) / nc_eta1

  SLL_CLEAR_ALLOCATE(adv%d_dx(1:nc_eta1), error)
  SLL_ALLOCATE(adv%fk(1:nc_eta1/2+1), error)
  SLL_ALLOCATE(adv%ek(1:nc_eta1/2+1), error)
  SLL_ALLOCATE(adv%r0(1:nc_eta1/2+1), error)
  SLL_ALLOCATE(adv%r1(1:nc_eta1/2+1), error)
  adv%fk = (0.0_f64, 0.0_f64)
  adv%ek = (0.0_f64, 0.0_f64)
  adv%r0 = (0.0_f64, 0.0_f64)
  adv%r1 = (0.0_f64, 0.0_f64)

  adv%fwx => fft_new_plan_r2c_1d(nc_eta1, adv%d_dx,  adv%fk)
  adv%bwx => fft_new_plan_c2r_1d(nc_eta1, adv%fk, adv%d_dx)

  SLL_CLEAR_ALLOCATE(adv%kx(1:nc_eta1/2+1), error)
   
  kx0 = 2._f64*sll_pi/(eta1_max-eta1_min)

  adv%kx(1) = 1.0_f64
  do i=2,nc_eta1/2+1
     adv%kx(i) = (i-1)*kx0
  end do
        
end subroutine initialize

subroutine advect_1d( adv, a, dt, input, output )

  class(ampere_1d_advector)             :: adv
  sll_real64, dimension(:), intent(in)  :: a
  sll_real64,               intent(in)  :: dt 
  sll_real64, dimension(:), intent(in)  :: input
  sll_real64, dimension(:), intent(out) :: output      

  character(len=*), parameter :: this_sub_name = 'advect_1d'
  sll_int32                   :: num_cells

  num_cells = adv%nc_eta1
  output(1:num_cells) = a(1:num_cells) * input(1:num_cells) * dt
  SLL_ERROR( this_sub_name, "ampere advect_1d not implemented." )

end subroutine advect_1d

subroutine delete(adv)

  class(ampere_1d_advector), intent(inout) :: adv

  call fft_delete_plan(adv%fwx)
  call fft_delete_plan(adv%bwx)

end subroutine delete

subroutine advect_1d_constant( adv, a, dt, input, output )

  class(ampere_1d_advector)             :: adv
  sll_real64,               intent(in)  :: a 
  sll_real64,               intent(in)  :: dt 
  sll_real64, dimension(:), intent(in)  :: input
  sll_real64, dimension(:), intent(out) :: output      

  sll_int32  :: i, nc_x
    
  nc_x = adv%nc_eta1

  adv%d_dx = input(1:nc_x)
  call fft_apply_plan_r2c_1d(adv%fwx, adv%d_dx, adv%fk)
  do i = 2, nc_x/2+1
    adv%fk(i) = adv%fk(i)*cmplx(cos(adv%kx(i)*a*dt),-sin(adv%kx(i)*a*dt),kind=f64)
  end do
  call fft_apply_plan_c2r_1d(adv%bwx, adv%fk, adv%d_dx)

  output(1:nc_x) = adv%d_dx / nc_x
  output(nc_x+1) = output(1)

end subroutine advect_1d_constant

end module sll_m_advection_1d_ampere

