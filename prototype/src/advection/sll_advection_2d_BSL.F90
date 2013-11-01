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

! in development; should be at least cubic splines
! attached with computation of characteristics


module sll_module_advection_2d_BSL
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
use sll_boundary_condition_descriptors
use sll_module_advection_2d_base
use sll_module_characteristics_2d_base
use sll_module_interpolators_2d_base
implicit none

  type,extends(sll_advection_2d_base) :: BSL_2d_advector
  
    class(sll_interpolator_2d_base), pointer  :: interp
    class(sll_characteristics_2d_base), pointer  :: charac
    sll_real64, dimension(:), pointer :: eta1_coords
    sll_real64, dimension(:), pointer :: eta2_coords
    sll_real64, dimension(:,:), pointer :: charac_foot1
    sll_real64, dimension(:,:), pointer :: charac_foot2
    sll_int32 :: Npts1
    sll_int32 :: Npts2  
  contains
     procedure, pass(adv) :: initialize => &
       initialize_BSL_2d_advector
    procedure, pass(adv) :: advect_2d => &
      BSL_advect_2d
  
  end type BSL_2d_advector
   




contains
  function new_BSL_2d_advector( &
    interp, &
    charac, &
    eta1_coords, &
    eta2_coords, &
    Npts1, &
    Npts2) &
    result(adv)      
    type(BSL_2d_advector), pointer :: adv
    class(sll_interpolator_2d_base), pointer :: interp
    class(sll_characteristics_2d_base), pointer  :: charac
    sll_int32 :: ierr
    sll_real64, dimension(:), pointer :: eta1_coords
    sll_real64, dimension(:), pointer :: eta2_coords
    sll_int32, intent(in) :: Npts1
    sll_int32, intent(in) :: Npts2
    
    SLL_ALLOCATE(adv,ierr)
    
    call initialize_BSL_2d_advector(&
      adv, &
      interp, &
      charac, &
      eta1_coords, &
      eta2_coords, &
      Npts1, &
      Npts2)
    
  end function  new_BSL_2d_advector


  subroutine initialize_BSL_2d_advector(&
    adv, &
    interp, &
    charac, &
    eta1_coords, &
    eta2_coords, &
    Npts1, &
    Npts2)
    class(BSL_2d_advector), intent(inout) :: adv
    class(sll_interpolator_2d_base), pointer :: interp
    class(sll_characteristics_2d_base), pointer  :: charac
    sll_real64, dimension(:), pointer :: eta1_coords
    sll_real64, dimension(:), pointer :: eta2_coords
    sll_int32, intent(in) :: Npts1
    sll_int32, intent(in) :: Npts2
    sll_int32 :: ierr
    
    
    adv%Npts1 = Npts1
    adv%Npts2 = Npts2
    adv%interp => interp
    adv%charac => charac
    !SLL_ALLOCATE(adv%x1_mesh(Npts1),ierr)
    !SLL_ALLOCATE(adv%x2_mesh(Npts2),ierr)

    SLL_ALLOCATE(adv%charac_foot1(Npts1,Npts2),ierr)
    SLL_ALLOCATE(adv%charac_foot2(Npts1,Npts2),ierr)

    
    if(size(eta1_coords,1)<Npts1)then
      print *,'#bad size for eta1_coords in initialize_BSL_2d_advector'
      stop
    endif
    if(size(eta2_coords,1)<Npts2)then
      print *,'#bad size for eta2_coords in initialize_BSL_2d_advector'
      stop
    endif
    
    adv%eta1_coords => eta1_coords
    adv%eta2_coords => eta2_coords
    
      
  end subroutine initialize_BSL_2d_advector

  subroutine BSL_advect_2d(&
    adv, &
    A1, &
    A2, &
    dt, &
    input, &
    output)
    class(BSL_2d_advector) :: adv
    sll_real64, dimension(:,:), intent(in) :: A1
    sll_real64, dimension(:,:), intent(in) :: A2
    sll_real64, intent(in) :: dt 
    sll_real64, dimension(:,:), intent(in) :: input
    sll_real64, dimension(:,:), intent(out) :: output      
    
    call adv%charac%compute_characteristics( &
      A1, &
      A2, &
      dt, &
      adv%eta1_coords, &
      adv%eta2_coords, &
      adv%charac_foot1, &
      adv%charac_foot2)
    
    
!    call adv%interp%compute_interpolants( &
!      input, &
!      adv%eta1_coords, &
!      adv%Npts1, &
!      adv%eta2_coords, &
!      adv%Npts2 )

    output = adv%interp%interpolate_array( &
      adv%Npts1, &
      adv%Npts2, &
      input, &
      adv%charac_foot1, &
      adv%charac_foot2)      
          
  end subroutine BSL_advect_2d





end module sll_module_advection_2d_BSL
