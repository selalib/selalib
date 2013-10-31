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

! for the moment mimic of sll_periodic_interpolator_1d.F90


module sll_module_advection_1d_periodic
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
use sll_boundary_condition_descriptors
use sll_module_advection_1d_base
use periodic_interp_module

implicit none

  type,extends(sll_advection_1d_base) :: periodic_1d_advector

     sll_int32                            :: num_cells
     sll_real64                           :: xmin
     sll_real64                           :: xmax
     type(periodic_interp_work), pointer  :: per_interp

  contains
    procedure, pass(adv) :: initialize => &
      initialize_periodic_1d_advector
    procedure, pass(adv) :: advect_1d_constant => &
      periodic_advect_1d_constant
  
  end type periodic_1d_advector
   
contains
  

  function new_periodic_1d_advector(&
    num_cells, &
    xmin, &
    xmax, &
    type, &
    order) &
    result(adv)      
    type(periodic_1d_advector), pointer :: adv
    sll_int32,  intent(in)               :: num_cells
    sll_real64, intent(in)               :: xmin
    sll_real64, intent(in)               :: xmax
    sll_int32,  intent(in)               :: type
    sll_int32,  intent(in)               :: order
    sll_int32 :: ierr

    SLL_ALLOCATE(adv,ierr)
    call initialize_periodic_1d_advector(&
      adv, &
      num_cells, &
      xmin, &
      xmax, &
      type, &
      order)
    
  end function new_periodic_1d_advector

  
  subroutine initialize_periodic_1d_advector(&
      adv, &
      num_cells, &
      xmin, &
      xmax, &
      type, &
      order)
      
    class(periodic_1d_advector) :: adv
    sll_int32,  intent(in)               :: num_cells
    sll_real64, intent(in)               :: xmin
    sll_real64, intent(in)               :: xmax
    sll_int32,  intent(in)               :: type
    sll_int32,  intent(in)               :: order

    call initialize_periodic_interp( &
      adv%per_interp, &
      num_cells, &
      type, &
      order)

    adv%num_cells = num_cells 
    adv%xmin = xmin 
    adv%xmax = xmax 
 
  
  end subroutine initialize_periodic_1d_advector   








  subroutine periodic_advect_1d_constant(&
    adv, &
    A, &
    dt, &
    input, &
    output)
    class(periodic_1d_advector) :: adv
    sll_real64, intent(in) :: A
    sll_real64, intent(in) :: dt 
    sll_real64, dimension(:), intent(in) :: input
    sll_real64, dimension(:), intent(out) :: output      
    sll_real64 :: shift
    sll_real64 :: xmin
    sll_real64 :: xmax
    sll_int32  :: num_cells
      
    num_cells = adv%num_cells
    xmin = adv%xmin
    xmax = adv%xmax
    shift = A*dt/(xmax-xmin)*real(num_cells,f64)
      
    call periodic_interp( &
      adv%per_interp, &
      output, &
      input, &
      shift)
    ! complete by periodicity
    output(num_cells+1) = output(1)
      
  end subroutine periodic_advect_1d_constant





end module sll_module_advection_1d_periodic