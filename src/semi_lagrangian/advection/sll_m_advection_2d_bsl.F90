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

module sll_m_advection_2d_bsl
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_errors.h"

use sll_m_advection_2d_base,       only: sll_c_advector_2d
use sll_m_characteristics_2d_base, only: sll_c_characteristics_2d_base
use sll_m_interpolators_2d_base,   only: sll_c_interpolator_2d

implicit none

public :: sll_f_new_advector_2d_bsl, &
          sll_t_advector_2d_bsl,     &
          sll_s_advector_2d_bsl_init

private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

type,extends(sll_c_advector_2d) :: sll_t_advector_2d_bsl

  class(sll_c_interpolator_2d),         pointer :: interp
  class(sll_c_characteristics_2d_base), pointer :: charac
  sll_real64, dimension(:),             pointer :: eta1_coords
  sll_real64, dimension(:),             pointer :: eta2_coords
  sll_real64, dimension(:,:),           pointer :: charac_feet1
  sll_real64, dimension(:,:),           pointer :: charac_feet2
  sll_int32                                     :: npts1
  sll_int32                                     :: npts2  

contains

  procedure, pass(adv) :: init => sll_s_advector_2d_bsl_init
  procedure, pass(adv) :: advect_2d => bsl_advect_2d

end type sll_t_advector_2d_bsl

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function sll_f_new_advector_2d_bsl( interp, charac, npts1, npts2, &
  eta1_min, eta1_max, eta2_min, eta2_max, eta1_coords, eta2_coords) &  
  result(adv)      

  type(sll_t_advector_2d_bsl),          pointer :: adv
  class(sll_c_interpolator_2d),         pointer :: interp
  class(sll_c_characteristics_2d_base), pointer :: charac
  sll_int32,  intent(in)                        :: npts1
  sll_int32,  intent(in)                        :: npts2
  sll_real64, intent(in),              optional :: eta1_min
  sll_real64, intent(in),              optional :: eta1_max
  sll_real64, intent(in),              optional :: eta2_min
  sll_real64, intent(in),              optional :: eta2_max
  sll_real64, dimension(:), pointer,   optional :: eta1_coords
  sll_real64, dimension(:), pointer,   optional :: eta2_coords
  sll_int32 :: ierr
  
  SLL_ALLOCATE(adv,ierr)
      
  call adv%init( interp, charac, npts1, npts2, eta1_min, eta1_max, &
    eta2_min, eta2_max, eta1_coords, eta2_coords)    
  
end function sll_f_new_advector_2d_bsl

subroutine sll_s_advector_2d_bsl_init( adv, interp, charac, npts1, npts2, &
  eta1_min, eta1_max, eta2_min, eta2_max, eta1_coords, eta2_coords)    

  class(sll_t_advector_2d_bsl)                       :: adv
  class(sll_c_interpolator_2d),               target :: interp
  class(sll_c_characteristics_2d_base),       target :: charac
  sll_int32,  intent(in)                             :: npts1
  sll_int32,  intent(in)                             :: npts2
  sll_real64, intent(in),                   optional :: eta1_min
  sll_real64, intent(in),                   optional :: eta1_max
  sll_real64, intent(in),                   optional :: eta2_min
  sll_real64, intent(in),                   optional :: eta2_max
  sll_real64, dimension(:),   pointer    ,  optional :: eta1_coords
  sll_real64, dimension(:),   pointer    ,  optional :: eta2_coords
  sll_int32                                          :: ierr
  sll_int32                                          :: i
  sll_real64                                         :: delta_eta1
  sll_real64                                         :: delta_eta2

  adv%npts1 = npts1
  adv%npts2 = npts2
  adv%interp => interp
  adv%charac => charac
  
  SLL_ALLOCATE(adv%eta1_coords(npts1),ierr)
  SLL_ALLOCATE(adv%eta2_coords(npts2),ierr)

  SLL_ALLOCATE(adv%charac_feet1(npts1,npts2),ierr)
  SLL_ALLOCATE(adv%charac_feet2(npts1,npts2),ierr)

  if(present(eta1_min).and.present(eta1_max))then
    if(present(eta1_coords))then
      SLL_ERROR("initialize_advector_2d_bsl","provide either eta1_coords or eta1_min and eta1_max and not both")
    else
      delta_eta1 = (eta1_max-eta1_min)/real(npts1-1,f64)
      do i=1,npts1
        adv%eta1_coords(i) = eta1_min+real(i-1,f64)*delta_eta1
      enddo
    endif
  else if(present(eta1_coords))then
    if(size(eta1_coords,1)<npts1)then
      SLL_ERROR("initialize_advector_2d_bsl","bad size for eta1_coords")
    else
      adv%eta1_coords(1:npts1) = eta1_coords(1:npts1)
    endif     
  else
    print *,'#Warning, we assume eta1_min = 0._f64 eta1_max = 1._f64'
    delta_eta1 = 1._f64/real(npts1-1,f64)
    do i=1,npts1
        adv%eta1_coords(i) = real(i-1,f64)*delta_eta1
    enddo                      
  endif


  if(present(eta2_min).and.present(eta2_max))then
    if(present(eta2_coords))then
      SLL_ERROR("initialize_advector_2d_bsl","provide either eta2_coords or eta2_min and eta2_max and not both")
    else
      delta_eta2 = (eta2_max-eta2_min)/real(npts2-1,f64)
      do i=1,npts2
        adv%eta2_coords(i) = eta2_min+real(i-1,f64)*delta_eta2
      enddo
    endif
  else if(present(eta2_coords))then
    if(size(eta2_coords,1)<npts2)then
      SLL_ERROR("initialize_advector_2d_bsl","bad size for eta2_coords")
    else
      adv%eta2_coords(1:npts2) = eta2_coords(1:npts2)
    endif     
  else
    print *,'#Warning, we assume eta2_min = 0._f64 eta2_max = 1._f64'
    delta_eta2 = 1._f64/real(npts2-1,f64)
    do i=1,npts2
        adv%eta2_coords(i) = real(i-1,f64)*delta_eta2
    enddo                      
  endif
  
end subroutine sll_s_advector_2d_bsl_init

subroutine bsl_advect_2d( adv, a1, a2, dt, input, output)

  class(sll_t_advector_2d_bsl)            :: adv
  sll_real64, dimension(:,:), intent(in)  :: a1
  sll_real64, dimension(:,:), intent(in)  :: a2
  sll_real64,                 intent(in)  :: dt 
  sll_real64, dimension(:,:), intent(in)  :: input
  sll_real64, dimension(:,:), intent(out) :: output      
  
  call adv%charac%compute_characteristics( &
    a1,                                    &
    a2,                                    &
    dt,                                    &
    adv%eta1_coords,                       &
    adv%eta2_coords,                       &
    adv%charac_feet1,                      &
    adv%charac_feet2)
  
  !call adv%interp%compute_interpolants( &
  !   input,                             &
  !   adv%eta1_coords,                   &
  !   adv%npts1,                         &
  !   adv%eta2_coords,                   &
  !   adv%npts2 )

  call adv%interp%interpolate_array( &
    adv%npts1,                       &
    adv%npts2,                       &
    input,                           &
    adv%charac_feet1,                &
    adv%charac_feet2,                &
    output)      
        
end subroutine bsl_advect_2d

end module sll_m_advection_2d_bsl
