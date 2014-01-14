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

! in development
! use of CSL1D in 2D
! begin with a specific example (for characteristics and interpolation)


module sll_module_advection_2d_tensor_product
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
use sll_boundary_condition_descriptors
use sll_module_advection_2d_base
use sll_module_advection_1d_base
use sll_module_characteristics_2d_base
use sll_module_interpolators_2d_base
implicit none

  type,extends(sll_advection_2d_base) :: tensor_product_2d_advector
    
    class(sll_advection_1d_base), pointer  :: advect_x1
    class(sll_advection_1d_base), pointer  :: advect_x2
    sll_int32 :: Npts1
    sll_int32 :: Npts2
    sll_real64, dimension(:), pointer :: buf1d 
  contains
     procedure, pass(adv) :: initialize => &
       initialize_tensor_product_2d_advector
    procedure, pass(adv) :: advect_2d => &
      tensor_product_advect_2d
  
  end type tensor_product_2d_advector
   




contains
  function new_tensor_product_2d_advector( &
    advect_x1, &
    advect_x2, &
    Npts1, &
    Npts2 ) &
    result(adv)      
    type(tensor_product_2d_advector), pointer :: adv
    class(sll_advection_1d_base), pointer :: advect_x1
    class(sll_advection_1d_base), pointer :: advect_x2
    sll_int32, intent(in) :: Npts1 
    sll_int32, intent(in) :: Npts2 
    sll_int32 :: ierr
    
    SLL_ALLOCATE(adv,ierr)
        
    call initialize_tensor_product_2d_advector(&
      adv, &
      advect_x1, &
      advect_x2, &
      Npts1, &
      Npts2)
    
  end function  new_tensor_product_2d_advector


  subroutine initialize_tensor_product_2d_advector(&
    adv, &
    advect_x1, &
    advect_x2, &
    Npts1, &
    Npts2)
    class(tensor_product_2d_advector), intent(inout) :: adv
    class(sll_advection_1d_base), pointer :: advect_x1
    class(sll_advection_1d_base), pointer :: advect_x2
    sll_int32, intent(in) :: Npts1
    sll_int32, intent(in) :: Npts2
    sll_int32 :: ierr
    adv%advect_x1 => advect_x1
    adv%advect_x2 => advect_x2
    
    adv%Npts1 = Npts1
    adv%Npts2 = Npts2
    
    SLL_ALLOCATE(adv%buf1d(max(Npts1,Npts2)),ierr)
      
  end subroutine initialize_tensor_product_2d_advector

  subroutine tensor_product_advect_2d(&
    adv, &
    A1, &
    A2, &
    dt, &
    input, &
    output)
    class(tensor_product_2d_advector) :: adv
    sll_real64, dimension(:,:), intent(in) :: A1
    sll_real64, dimension(:,:), intent(in) :: A2
    sll_real64, intent(in) :: dt 
    sll_real64, dimension(:,:), intent(in) :: input
    sll_real64, dimension(:,:), intent(out) :: output      
    sll_int32 :: i1
    sll_int32 :: i2
    sll_int32 :: Npts1
    sll_int32 :: Npts2
    
    Npts1 = adv%Npts1
    Npts2 = adv%Npts2
    
    do i2=1,Npts2
      adv%buf1d(1:Npts1) = input(1:Npts1,i2)
      call adv%advect_x1%advect_1d( &
        A1(1:Npts1,i2), &
        0.5_f64*dt, &
        adv%buf1d(1:Npts1), &
        output(1:Npts1,i2))
    enddo

    do i1=1,Npts1
      adv%buf1d(1:Npts2) = output(i1,1:Npts2)
      call adv%advect_x2%advect_1d( &
        A2(i1,1:Npts2), &
        dt, &
        adv%buf1d(1:Npts2), &
        output(i1,1:Npts2))
    enddo

    do i2=1,Npts2
      adv%buf1d(1:Npts1) = input(1:Npts1,i2)
      call adv%advect_x1%advect_1d( &
        A1(1:Npts1,i2), &
        0.5_f64*dt, &
        adv%buf1d(1:Npts1), &
        output(1:Npts1,i2))
    enddo
          
  end subroutine tensor_product_advect_2d





end module sll_module_advection_2d_tensor_product
