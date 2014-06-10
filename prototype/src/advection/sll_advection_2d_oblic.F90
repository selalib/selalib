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

! in development; use oblic interpolation
! data are on uniform (fine) gird in x1
! 

module sll_module_advection_2d_oblic
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
use sll_boundary_condition_descriptors
!use sll_module_advection_2d_base
use sll_module_characteristics_2d_base
use sll_module_interpolators_2d_base
use sll_fcisl_module
implicit none

  type  :: oblic_2d_advector
  
    class(sll_advection_1d_base), pointer :: adv_x1
    class(sll_advection_1d_base), pointer :: adv_aligned
    class(sll_interpolator_2d_base), pointer  :: interp
    class(sll_characteristics_2d_base), pointer  :: charac
    sll_real64, dimension(:), pointer :: eta1_coords
    sll_real64, dimension(:), pointer :: eta2_coords
    sll_real64, dimension(:,:), pointer :: charac_feet1
    sll_real64, dimension(:,:), pointer :: charac_feet2
    sll_real64, dimension(:,:), pointer :: phi_at_aligned
    sll_int32, dimension(:), pointer :: spaghetti_index
    sll_real64, dimension(:), pointer :: phi_at_aligned_1d
    sll_real64, dimension(:), pointer :: charac_feet_aligned_1d
    sll_int32 :: spaghetti_size
    sll_int32 :: Npts1
    sll_int32 :: Npts2
    sll_int32 :: shift  
  contains
  !   procedure, pass(adv) :: initialize => &
  !     initialize_oblic_2d_advector
  !  procedure, pass(adv) :: advect_2d => &
  !    oblic_advect_2d
  
  end type oblic_2d_advector
   




contains
  function new_oblic_2d_advector( &
    adv_x1, &
    adv_aligned, &
    interp, &
    charac, &
    Npts1, &
    Npts2, &
    eta1_min, &
    eta1_max, &
    eta2_min, &
    eta2_max, &
    eta1_coords, &
    eta2_coords) &  
    result(adv)      
    type(oblic_2d_advector), pointer :: adv
    class(sll_advection_1d_base), pointer :: adv_x1
    class(sll_advection_1d_base), pointer :: adv_aligned
    class(sll_interpolator_2d_base), pointer :: interp
    class(sll_characteristics_2d_base), pointer  :: charac
    sll_int32, intent(in) :: Npts1
    sll_int32, intent(in) :: Npts2
    sll_real64, intent(in), optional :: eta1_min
    sll_real64, intent(in), optional :: eta1_max
    sll_real64, intent(in), optional :: eta2_min
    sll_real64, intent(in), optional :: eta2_max
    sll_real64, dimension(:), pointer, optional :: eta1_coords
    sll_real64, dimension(:), pointer, optional :: eta2_coords
    sll_int32 :: ierr
    
    SLL_ALLOCATE(adv,ierr)
        
    call initialize_oblic_2d_advector(&
      adv, &
      adv_x1, &
      adv_aligned, &
      interp, &
      charac, &
      Npts1, &
      Npts2, &
      eta1_min, &
      eta1_max, &
      eta2_min, &
      eta2_max, &
      eta1_coords, &
      eta2_coords)    
    
  end function  new_oblic_2d_advector


  subroutine initialize_oblic_2d_advector(&
    adv, &
    adv_x1, &
    adv_aligned, &
    interp, &
    charac, &
    Npts1, &
    Npts2, &
    eta1_min, &
    eta1_max, &
    eta2_min, &
    eta2_max, &
    eta1_coords, &
    eta2_coords)    
    type(oblic_2d_advector), intent(inout) :: adv
    class(sll_advection_1d_base), pointer :: adv_x1
    class(sll_advection_1d_base), pointer :: adv_aligned
    class(sll_interpolator_2d_base), pointer :: interp
    class(sll_characteristics_2d_base), pointer  :: charac
    sll_int32, intent(in) :: Npts1
    sll_int32, intent(in) :: Npts2
    sll_real64, intent(in), optional :: eta1_min
    sll_real64, intent(in), optional :: eta1_max
    sll_real64, intent(in), optional :: eta2_min
    sll_real64, intent(in), optional :: eta2_max
    sll_real64, dimension(:), pointer, optional :: eta1_coords
    sll_real64, dimension(:), pointer, optional :: eta2_coords
    sll_int32 :: ierr
    sll_int32 :: i
    sll_real64 :: delta_eta1
    sll_real64 :: delta_eta2
    
    adv%shift = 0
    adv%Npts1 = Npts1
    adv%Npts2 = Npts2
    adv%adv_x1 => adv_x1
    adv%adv_aligned => adv_aligned
    adv%interp => interp
    adv%charac => charac
    !SLL_ALLOCATE(adv%x1_mesh(Npts1),ierr)
    !SLL_ALLOCATE(adv%x2_mesh(Npts2),ierr)
    SLL_ALLOCATE(adv%eta1_coords(Npts1),ierr)
    SLL_ALLOCATE(adv%eta2_coords(Npts2),ierr)

    SLL_ALLOCATE(adv%charac_feet1(Npts1,Npts2),ierr)
    SLL_ALLOCATE(adv%charac_feet2(Npts1,Npts2),ierr)

    SLL_ALLOCATE(adv%phi_at_aligned(Npts1,Npts2),ierr)
    SLL_ALLOCATE(adv%spaghetti_index(Npts1),ierr)
    SLL_ALLOCATE(adv%phi_at_aligned_1d(Npts1*Npts2),ierr)
    SLL_ALLOCATE(adv%charac_feet_aligned_1d(Npts1*Npts2),ierr)


    if(present(eta1_min).and.present(eta1_max))then
      if(present(eta1_coords))then
        print *,'#provide either eta1_coords or eta1_min and eta1_max'
        print *,'#and not both in subroutine initialize_BSL_2d_advector'
        stop
      else
        delta_eta1 = (eta1_max-eta1_min)/real(Npts1-1,f64)
        do i=1,Npts1
          adv%eta1_coords(i) = eta1_min+real(i-1,f64)*delta_eta1
        enddo
      endif
    else if(present(eta1_coords))then
      if(size(eta1_coords,1)<Npts1)then
        print *,'#bad size for eta1_coords in initialize_BSL_2d_advector'
        stop
      else
        adv%eta1_coords(1:Npts1) = eta1_coords(1:Npts1)
      endif     
    else
      print *,'#Warning, we assume eta1_min = 0._f64 eta1_max = 1._f64'
      delta_eta1 = 1._f64/real(Npts1-1,f64)
      do i=1,Npts1
          adv%eta1_coords(i) = real(i-1,f64)*delta_eta1
      enddo                      
    endif


    if(present(eta2_min).and.present(eta2_max))then
      if(present(eta2_coords))then
        print *,'#provide either eta2_coords or eta2_min and eta2_max'
        print *,'#and not both in subroutine initialize_BSL_2d_advector'
        stop
      else
        delta_eta2 = (eta2_max-eta2_min)/real(Npts2-1,f64)
        do i=1,Npts2
          adv%eta2_coords(i) = eta2_min+real(i-1,f64)*delta_eta2
        enddo
      endif
    else if(present(eta2_coords))then
      if(size(eta2_coords,1)<Npts2)then
        print *,'#bad size for eta2_coords in initialize_BSL_2d_advector'
        stop
      else
        adv%eta2_coords(1:Npts2) = eta2_coords(1:Npts2)
      endif     
    else
      print *,'#Warning, we assume eta2_min = 0._f64 eta2_max = 1._f64'
      delta_eta2 = 1._f64/real(Npts2-1,f64)
      do i=1,Npts2
          adv%eta2_coords(i) = real(i-1,f64)*delta_eta2
      enddo                      
    endif
    
      
  end subroutine initialize_oblic_2d_advector


!  subroutine set_shift(adv, shift)
!    type(oblic_2d_advector) :: adv  
!    sll_int32, intent(in) :: shift
!    adv%shift = shift
!  end subroutine set_shift
!
!
!  function get_shift(adv) result(shift)
!    type(oblic_2d_advector) :: adv  
!    sll_int32 :: shift
!    shift = adv%shift
!  end function get_shift


  subroutine oblic_advect_2d(&
    adv, &
    phi, &
    shift, &
    dt, &
    input, &
    output)
    type(oblic_2d_advector) :: adv
    sll_real64, dimension(:,:), intent(in) :: phi
    sll_int32, intent(in) :: shift
    sll_real64, intent(in) :: dt 
    sll_real64, dimension(:,:), intent(in) :: input
    sll_real64, dimension(:,:), intent(out) :: output      
    sll_int32 :: Npts1
    sll_int32 :: Npts2
    sll_real64 :: A
    sll_int32 :: i
    sll_real64 :: delta_x1
    print *,'#not implemented for the moment'
    
    Npts1 = adv%Npts1
    Npts2 = adv%Npts2
    A = real(shift,f64)/real(Npts2-1,f64)
    delta_x1 = (adv%eta1_coords(Npts1)-adv%eta1_coords(1))/real(Npts1-1,f64)

    !first compute phi at aligned points
    do i=1,Npts2
      call adv%adv_x1%advect_1d_constant( &
        A, &
        real(i-1,f64)*delta_x1, &
        phi(1:Npts1,i), &
        adv%phi_at_aligned(1:Npts1,i))      
    enddo

    call compute_spaghetti( &
      Npts1-1, &
      Npts2-1, &
      shift, &
      adv%spaghetti_index, &
      adv%spaghetti_size)
    
    call load_spaghetti( &
      adv%phi_at_aligned, &
      adv%phi_at_aligned_1d, &
      adv%spaghetti_index, &
      Npts1, &
      Npts2) 

!  do i=1,num_spaghetti
!    call adv%advect_1d_constant( &
!        tau, &
!        dt, &
!        buf1d_spaghetti(s:s+spaghetti_size), &
!        buf1d_spaghetto(1:spaghetti_size))      
!
!    buf1d_spaghetti(s:s+spaghetti_size)=buf1d_spaghetto(1:spaghetti_size) 
!    
!    s = s+spaghetti_size
!  enddo   

     
!    call adv%charac%compute_characteristics( &
!      A1, &
!      A2, &
!      dt, &
!      adv%eta1_coords, &
!      adv%eta2_coords, &
!      adv%charac_feet1, &
!      adv%charac_feet2)
    
    
!    call adv%interp%compute_interpolants( &
!      input, &
!      adv%eta1_coords, &
!      adv%Npts1, &
!      adv%eta2_coords, &
!      adv%Npts2 )

!    output = adv%interp%interpolate_array( &
!      adv%Npts1, &
!      adv%Npts2, &
!      input, &
!      adv%charac_feet1, &
!      adv%charac_feet2)      
          
  end subroutine oblic_advect_2d





end module sll_module_advection_2d_oblic
