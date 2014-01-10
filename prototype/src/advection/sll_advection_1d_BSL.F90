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


module sll_module_advection_1d_BSL
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
use sll_boundary_condition_descriptors
use sll_module_advection_1d_base
use sll_module_characteristics_1d_base
use sll_module_interpolators_1d_base
implicit none

  type,extends(sll_advection_1d_base) :: BSL_1d_advector
  
    class(sll_interpolator_1d_base), pointer  :: interp
    class(sll_characteristics_1d_base), pointer  :: charac
    sll_real64, dimension(:), pointer :: eta_coords
    sll_real64, dimension(:), pointer :: charac_feet
    sll_int32 :: Npts
  contains
    procedure, pass(adv) :: initialize => &
       initialize_BSL_1d_advector
    procedure, pass(adv) :: advect_1d => &
      BSL_advect_1d
    procedure, pass(adv) :: advect_1d_constant => &
      BSL_advect_1d_constant
  
  end type BSL_1d_advector
   




contains
  function new_BSL_1d_advector( &
    interp, &
    charac, &
    Npts, &
    eta_min, &
    eta_max, &
    eta_coords) &  
    result(adv)      
    type(BSL_1d_advector), pointer :: adv
    class(sll_interpolator_1d_base), pointer :: interp
    class(sll_characteristics_1d_base), pointer  :: charac
    sll_int32, intent(in) :: Npts
    sll_real64, intent(in), optional :: eta_min
    sll_real64, intent(in), optional :: eta_max
    sll_real64, dimension(:), pointer, optional :: eta_coords
    sll_int32 :: ierr
    
    SLL_ALLOCATE(adv,ierr)
        
    call initialize_BSL_1d_advector(&
      adv, &
      interp, &
      charac, &
      Npts, &
      eta_min, &
      eta_max, &
      eta_coords)    
    
  end function  new_BSL_1d_advector


  subroutine initialize_BSL_1d_advector(&
    adv, &
    interp, &
    charac, &
    Npts, &
    eta_min, &
    eta_max, &
    eta_coords)    
    class(BSL_1d_advector), intent(inout) :: adv
    class(sll_interpolator_1d_base), pointer :: interp
    class(sll_characteristics_1d_base), pointer  :: charac
    sll_int32, intent(in) :: Npts
    sll_real64, intent(in), optional :: eta_min
    sll_real64, intent(in), optional :: eta_max
    sll_real64, dimension(:), pointer, optional :: eta_coords
    sll_int32 :: ierr
    sll_int32 :: i
    sll_real64 :: delta_eta
    
    
    adv%Npts = Npts
    adv%interp => interp
    adv%charac => charac
    SLL_ALLOCATE(adv%eta_coords(Npts),ierr)

    SLL_ALLOCATE(adv%charac_feet(Npts),ierr)

    if(present(eta_min).and.present(eta_max))then
      if(present(eta_coords))then
        print *,'#provide either eta_coords or eta_min and eta_max'
        print *,'#and not both in subroutine initialize_BSL_1d_advector'
        stop
      else
        delta_eta = (eta_max-eta_min)/real(Npts-1,f64)
        do i=1,Npts
          adv%eta_coords(i) = eta_min+real(i-1,f64)*delta_eta
        enddo
      endif
    else if(present(eta_coords))then
      if(size(eta_coords)<Npts)then
        print *,'#bad size for eta_coords in initialize_BSL_1d_advector'
        stop
      else
        adv%eta_coords(1:Npts) = eta_coords(1:Npts)
      endif     
    else
      print *,'#Warning, we assume eta_min = 0._f64 eta_max = 1._f64'
      delta_eta = 1._f64/real(Npts-1,f64)
      do i=1,Npts
          adv%eta_coords(i) = real(i-1,f64)*delta_eta
      enddo                      
    endif
          
  end subroutine initialize_BSL_1d_advector

  subroutine BSL_advect_1d(&
    adv, &
    A, &
    dt, &
    input, &
    output)
    class(BSL_1d_advector) :: adv
    sll_real64, dimension(:), intent(in) :: A
    sll_real64, intent(in) :: dt 
    sll_real64, dimension(:), intent(in) :: input
    sll_real64, dimension(:), intent(out) :: output      
    
    call adv%charac%compute_characteristics( &
      A, &
      dt, &
      adv%eta_coords, &
      adv%charac_feet)
    
    
!    call adv%interp%compute_interpolants( &
!      input, &
!      adv%eta1_coords, &
!      adv%Npts1, &
!      adv%eta2_coords, &
!      adv%Npts2 )

    output = adv%interp%interpolate_array( &
      adv%Npts, &
      input, &
      adv%charac_feet)      
          
  end subroutine BSL_advect_1d


  subroutine BSL_advect_1d_constant(&
    adv, &
    A, &
    dt, &
    input, &
    output)
    class(BSL_1d_advector) :: adv
    sll_real64, intent(in) :: A
    sll_real64, intent(in) :: dt 
    sll_real64, dimension(:), intent(in) :: input
    sll_real64, dimension(:), intent(out) :: output      
    sll_real64, dimension(:), allocatable :: A1
    sll_int32 :: ierr
    
    !this version is not optimized
    
    SLL_ALLOCATE(A1(adv%Npts),ierr)
    
    A1 = A
    
    call adv%charac%compute_characteristics( &
      A1, &
      dt, &
      adv%eta_coords, &
      adv%charac_feet)

!    call adv%interp%compute_interpolants( &
!      input, &
!      adv%eta1_coords, &
!      adv%Npts1, &
!      adv%eta2_coords, &
!      adv%Npts2 )

    output = adv%interp%interpolate_array( &
      adv%Npts, &
      input, &
      adv%charac_feet)      

    SLL_DEALLOCATE_ARRAY(A1,ierr)

          
  end subroutine BSL_advect_1d_constant






end module sll_module_advection_1d_BSL
