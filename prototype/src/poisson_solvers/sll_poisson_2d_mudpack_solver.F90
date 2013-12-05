!**************************************************************
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


!solves \sum_{i,j=1}^2 A_{i,j}\partial_{i,j} phi
!       +\sum_{i=1}^2B_i\partial_i phi
!       +C \phi = rho
!in polar coordinates
!this leads when A_{1,2}=A_{2,1}=0 and B_2 = 0
! A_11\partial_{1,1}\hat{phi}+B_1\partial_{1}\hat{phi}+(C+A_{2,2}k^2)\hat{phi} = \hat{rho}


module sll_module_poisson_2d_mudpack_solver
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
!use sll_boundary_condition_descriptors
use sll_module_poisson_2d_base
!use sll_poisson_2d_polar
implicit none

  integer, parameter :: SLL_SEPARABLE  = 1    !< type of equation
  integer, parameter :: SLL_NON_SEPARABLE_WITHOUT_CROSS_TERMS = 2    !< type of equation
  integer, parameter :: SLL_NON_SEPARABLE_WITH_CROSS_TERMS = 3    !< type of equation
  

  type,extends(sll_poisson_2d_base) :: poisson_2d_mudpack_solver     
  
  !type(sll_plan_poisson_polar), pointer                   :: poiss
  sll_real64, dimension(:,:), pointer :: cxx_2d
  sll_real64, dimension(:,:), pointer :: cxy_2d
  sll_real64, dimension(:,:), pointer :: cyy_2d
  sll_real64, dimension(:,:), pointer :: cx_2d
  sll_real64, dimension(:,:), pointer :: cy_2d
  sll_real64, dimension(:,:), pointer :: ce_2d
  sll_real64, dimension(:), pointer :: cxx_1d
  sll_real64, dimension(:), pointer :: cyy_1d
  sll_real64, dimension(:), pointer :: cx_1d
  sll_real64, dimension(:), pointer :: cy_1d
  sll_real64, dimension(:), pointer :: cex_1d
  sll_real64, dimension(:), pointer :: cey_1d
  sll_real64 :: cxx
  sll_real64 :: cyy
  sll_real64 :: cx
  sll_real64 :: cy
  sll_real64 :: ce
  sll_int32  :: mudpack_case
  
  
  contains
    procedure, pass(poisson) :: initialize => &
      initialize_poisson_2d_mudpack_solver
    procedure, pass(poisson) :: compute_phi_from_rho => &
      compute_phi_from_rho_2d_mudpack
    procedure, pass(poisson) :: compute_E_from_rho => &
      compute_E_from_rho_2d_mudpack
!    procedure, pass(poisson) :: compute_E_from_phi => &
!      compute_E_from_phi_2d_polar
      
  end type poisson_2d_mudpack_solver

contains
  function new_poisson_2d_mudpack_solver( &
    eta1_min, &
    eta1_max, &
    nc_eta1, &
    eta2_min, &
    eta2_max, &
    nc_eta2, &    
    bc_eta1_left, &
    bc_eta1_right, &
    bc_eta2_left, &
    bc_eta2_right, &
    mudpack_case, &
    cxx_2d, &
    cxy_2d, &
    cyy_2d, &
    cx_2d, &
    cy_2d, &
    ce_2d, &
    cxx_1d, &
    cyy_1d, &
    cx_1d, &
    cy_1d, &
    cex_1d, &
    cey_1d, &
    cxx, &
    cyy, &
    cx, &
    cy, &
    ce) &
    result(poisson)
      
    type(poisson_2d_mudpack_solver),pointer :: poisson
    sll_real64, intent(in) :: eta1_min
    sll_real64, intent(in) :: eta1_max
    sll_int32, intent(in) :: nc_eta1
    sll_real64, intent(in) :: eta2_min
    sll_real64, intent(in) :: eta2_max
    sll_int32, intent(in) :: nc_eta2
    sll_int32, intent(in) :: bc_eta1_left
    sll_int32, intent(in) :: bc_eta1_right
    sll_int32, intent(in) :: bc_eta2_left
    sll_int32, intent(in) :: bc_eta2_right
    sll_int32, intent(in) :: mudpack_case
    sll_real64, dimension(:,:), intent(in), optional :: cxx_2d
    sll_real64, dimension(:,:), intent(in), optional :: cxy_2d
    sll_real64, dimension(:,:), intent(in), optional :: cyy_2d
    sll_real64, dimension(:,:), intent(in), optional :: cx_2d
    sll_real64, dimension(:,:), intent(in), optional :: cy_2d
    sll_real64, dimension(:,:), intent(in), optional :: ce_2d
    sll_real64, dimension(:), intent(in), optional :: cxx_1d
    sll_real64, dimension(:), intent(in), optional :: cyy_1d
    sll_real64, dimension(:), intent(in), optional :: cx_1d
    sll_real64, dimension(:), intent(in), optional :: cy_1d
    sll_real64, dimension(:), intent(in), optional :: cex_1d
    sll_real64, dimension(:), intent(in), optional :: cey_1d
    sll_real64, intent(in), optional  :: cxx
    sll_real64, intent(in), optional  :: cyy
    sll_real64, intent(in), optional  :: cx
    sll_real64, intent(in), optional  :: cy
    sll_real64, intent(in), optional  :: ce

    sll_int32 :: ierr
      
    SLL_ALLOCATE(poisson,ierr)
    call initialize_poisson_2d_mudpack_solver( &
      poisson, &
      eta1_min, &
      eta1_max, &
      nc_eta1, &
      eta2_min, &
      eta2_max, &
      nc_eta2, &    
      bc_eta1_left, &
      bc_eta1_right, &
      bc_eta2_left, &
      bc_eta2_right, &
      mudpack_case, &
      cxx_2d, &
      cxy_2d, &
      cyy_2d, &
      cx_2d, &
      cy_2d, &
      ce_2d, &
      cxx_1d, &
      cyy_1d, &
      cx_1d, &
      cy_1d, &
      cex_1d, &
      cey_1d, &
      cxx, &
      cyy, &
      cx, &
      cy, &
      ce)
    
  end function new_poisson_2d_mudpack_solver
  
  
  subroutine initialize_poisson_2d_mudpack_solver( &
    poisson, &
    eta1_min, &
    eta1_max, &
    nc_eta1, &
    eta2_min, &
    eta2_max, &
    nc_eta2, &    
    bc_eta1_left, &
    bc_eta1_right, &
    bc_eta2_left, &
    bc_eta2_right, &
    mudpack_case, &
    cxx_2d, &
    cxy_2d, &
    cyy_2d, &
    cx_2d, &
    cy_2d, &
    ce_2d, &
    cxx_1d, &
    cyy_1d, &
    cx_1d, &
    cy_1d, &
    cex_1d, &
    cey_1d, &
    cxx, &
    cxy, &
    cyy, &
    cx, &
    cy, &
    ce)
    class(poisson_2d_mudpack_solver) :: poisson
    sll_real64, intent(in) :: eta1_min
    sll_real64, intent(in) :: eta1_max
    sll_int32, intent(in) :: nc_eta1
    sll_real64, intent(in) :: eta2_min
    sll_real64, intent(in) :: eta2_max
    sll_int32, intent(in) :: nc_eta2
    sll_int32, intent(in) :: bc_eta1_left
    sll_int32, intent(in) :: bc_eta1_right
    sll_int32, intent(in) :: bc_eta2_left
    sll_int32, intent(in) :: bc_eta2_right
    sll_int32, intent(in) :: mudpack_case
    sll_real64, dimension(:,:), intent(in), optional :: cxx_2d
    sll_real64, dimension(:,:), intent(in), optional :: cxy_2d
    sll_real64, dimension(:,:), intent(in), optional :: cyy_2d
    sll_real64, dimension(:,:), intent(in), optional :: cx_2d
    sll_real64, dimension(:,:), intent(in), optional :: cy_2d
    sll_real64, dimension(:,:), intent(in), optional :: ce_2d
    sll_real64, dimension(:), intent(in), optional :: cxx_1d
    sll_real64, dimension(:), intent(in), optional :: cyy_1d
    sll_real64, dimension(:), intent(in), optional :: cx_1d
    sll_real64, dimension(:), intent(in), optional :: cy_1d
    sll_real64, dimension(:), intent(in), optional :: cex_1d
    sll_real64, dimension(:), intent(in), optional :: cey_1d
    sll_real64, intent(in), optional  :: cxx
    sll_real64, intent(in), optional  :: cxy
    sll_real64, intent(in), optional  :: cyy
    sll_real64, intent(in), optional  :: cx
    sll_real64, intent(in), optional  :: cy
    sll_real64, intent(in), optional  :: ce
    sll_int32 :: ierr
        
    poisson%mudpack_case = mudpack_case 

    select case (mudpack_case)
      case (SLL_SEPARABLE)
        if(present(cxx_2d).or.present(cxy_2d).or.present(cyy_2d)&
          .or.present(cx_2d).or.present(cy_2d).or.present(ce_2d)) then
          print *,'#2d arrays should not be here'
          print *,'#in subroutine initialize_poisson_2d_mudpack_solver'
          stop
        endif
        
        if((.not.(present(cxx_1d))).and.(.not.(present(cxx)))) then
          print *,'#1d/0d array should be here for cxx'
          print *,'#in subroutine initialize_poisson_2d_mudpack_solver'
          stop        
        endif
        if(present(cxx_1d).and.present(cxx))then
          print *,'#please choose between 1d or 0d array for cxx'
          print *,'#in subroutine initialize_poisson_2d_mudpack_solver'
          stop        
        endif
        if(present(cxx_1d))then
          if(size(cxx_1d)<nc_eta1+1)then
            print *,'#Bad size for cxx_1d',size(cxx_1d),nc_eta1+1
            stop
          endif
          SLL_ALLOCATE(poisson%cxx_1d(nc_eta1+1),ierr)
          poisson%cxx_1d(1:nc_eta1+1)=cxx_1d(1:nc_eta1+1)
        endif
        if(present(cxx))then
          SLL_ALLOCATE(poisson%cxx_1d(nc_eta1+1),ierr)
          poisson%cxx_1d(1:nc_eta1+1)=cxx
        endif

        if((.not.(present(cyy_1d))).and.(.not.(present(cyy)))) then
          print *,'#1d/0d array should be here for cyy'
          print *,'#in subroutine initialize_poisson_2d_mudpack_solver'
          stop        
        endif
        if(present(cyy_1d).and.present(cyy))then
          print *,'#please choose between 1d or 0d array for cyy'
          print *,'#in subroutine initialize_poisson_2d_mudpack_solver'
          stop        
        endif
        if(present(cyy_1d))then
          if(size(cyy_1d)<nc_eta2+1)then
            print *,'#Bad size for cyy_1d',size(cyy_1d),nc_eta2+1
            stop
          endif
          SLL_ALLOCATE(poisson%cyy_1d(nc_eta2+1),ierr)
          poisson%cyy_1d(1:nc_eta2+1)=cyy_1d(1:nc_eta2+1)
        endif
        if(present(cyy))then
          SLL_ALLOCATE(poisson%cyy_1d(nc_eta2+1),ierr)
          poisson%cyy_1d(1:nc_eta2+1)=cyy
        endif

        if((.not.(present(cx_1d))).and.(.not.(present(cx)))) then
          SLL_ALLOCATE(poisson%cx_1d(nc_eta1+1),ierr)
          poisson%cx_1d(1:nc_eta1+1)=0._f64
        endif
        if(present(cx_1d).and.present(cx))then
          print *,'#please choose between 1d or 0d array for cx'
          print *,'#in subroutine initialize_poisson_2d_mudpack_solver'
          stop        
        endif
        if(present(cx_1d))then
          if(size(cx_1d)<nc_eta1+1)then
            print *,'#Bad size for cx_1d',size(cx_1d),nc_eta1+1
            stop
          endif
          SLL_ALLOCATE(poisson%cx_1d(nc_eta1+1),ierr)
          poisson%cx_1d(1:nc_eta1+1)=cx_1d(1:nc_eta1+1)
        endif
        if(present(cx))then
          SLL_ALLOCATE(poisson%cx_1d(nc_eta1+1),ierr)
          poisson%cx_1d(1:nc_eta1+1)=cx
        endif


        if((.not.(present(cy_1d))).and.(.not.(present(cy)))) then
          SLL_ALLOCATE(poisson%cy_1d(nc_eta2+1),ierr)
          poisson%cy_1d(1:nc_eta2+1)=0._f64
        endif
        if(present(cy_1d).and.present(cy))then
          print *,'#please choose between 1d or 0d array for cy'
          print *,'#in subroutine initialize_poisson_2d_mudpack_solver'
          stop        
        endif
        if(present(cy_1d))then
          if(size(cy_1d)<nc_eta2+1)then
            print *,'#Bad size for cy_1d',size(cy_1d),nc_eta2+1
            stop
          endif
          SLL_ALLOCATE(poisson%cy_1d(nc_eta2+1),ierr)
          poisson%cy_1d(1:nc_eta2+1)=cy_1d(1:nc_eta2+1)
        endif
        if(present(cy))then
          SLL_ALLOCATE(poisson%cy_1d(nc_eta2+1),ierr)
          poisson%cy_1d(1:nc_eta2+1)=cy
        endif

        if((.not.(present(cex_1d))).and.(.not.(present(ce)))) then
          SLL_ALLOCATE(poisson%cex_1d(nc_eta1+1),ierr)
          poisson%cex_1d = 0._f64          
        endif
        if(present(cex_1d).and.present(ce))then
          print *,'#please choose between 1d or 0d array for cex'
          print *,'#in subroutine initialize_poisson_2d_mudpack_solver'
          stop        
        endif
        if(present(cex_1d))then
          if(size(cex_1d)<nc_eta1+1)then
            print *,'#Bad size for cex_1d',size(cex_1d),nc_eta1+1
            stop
          endif
          SLL_ALLOCATE(poisson%cex_1d(nc_eta1+1),ierr)
          poisson%cex_1d(1:nc_eta1+1)=cex_1d(1:nc_eta1+1)
        endif
        if(present(ce))then
          SLL_ALLOCATE(poisson%cex_1d(nc_eta1+1),ierr)
          poisson%cex_1d(1:nc_eta1+1)=0.5_f64*ce
        endif

        if((.not.(present(cey_1d))).and.(.not.(present(ce)))) then
          SLL_ALLOCATE(poisson%cey_1d(nc_eta2+1),ierr)
          poisson%cey_1d(1:nc_eta2+1)=0._f64
        endif
        if(present(cey_1d).and.present(ce))then
          print *,'#please choose between 1d or 0d array for cey'
          print *,'#in subroutine initialize_poisson_2d_mudpack_solver'
          stop        
        endif
        if(present(cey_1d))then
          if(size(cey_1d)<nc_eta2+1)then
            print *,'#Bad size for cey_1d',size(cey_1d),nc_eta2+1
            stop
          endif
          SLL_ALLOCATE(poisson%cey_1d(nc_eta2+1),ierr)
          poisson%cey_1d(1:nc_eta2+1)=cey_1d(1:nc_eta2+1)
        endif
        if(present(ce))then
          SLL_ALLOCATE(poisson%cey_1d(nc_eta2+1),ierr)
          poisson%cey_1d(1:nc_eta2+1)=0.5_f64*ce
        endif









      
      
        
      case (SLL_NON_SEPARABLE_WITHOUT_CROSS_TERMS)
      case (SLL_NON_SEPARABLE_WITH_CROSS_TERMS)
      case default
        print *,'#bad mudpack_case',mudpack_case
        print *,'#in subroutine initialize_poisson_2d_mudpack_solver'
        stop 
    end select

        
  end subroutine initialize_poisson_2d_mudpack_solver
  
  ! solves -\Delta phi = rho in 2d
  subroutine compute_phi_from_rho_2d_mudpack( poisson, phi, rho )
    class(poisson_2d_mudpack_solver) :: poisson
    sll_real64,dimension(:,:),intent(in) :: rho
    sll_real64,dimension(:,:),intent(out) :: phi
    
    !call solve( poisson%poiss, rho, phi)
    
  end subroutine compute_phi_from_rho_2d_mudpack

    ! solves E = -\nabla Phi in 2d
!    subroutine compute_E_from_phi_2d_fft( poisson, phi, E1, E2 )
!      class(poisson_2d_fft_solver) :: poisson
!      sll_real64,dimension(:,:),intent(in) :: phi
!      sll_real64,dimension(:,:),intent(out) :: E1
!      sll_real64,dimension(:,:),intent(out) :: E2
!    end subroutine compute_E_from_phi_2d_fft

    ! solves E = -\nabla Phi with -\Delta phi = rho in 2d 
    subroutine compute_E_from_rho_2d_mudpack( poisson, rho, E1, E2 )
      class(poisson_2d_mudpack_solver) :: poisson
      sll_real64,dimension(:,:),intent(in) :: rho
      sll_real64,dimension(:,:),intent(out) :: E1
      sll_real64,dimension(:,:),intent(out) :: E2
      
      print *,'#compute_E_from_rho_2d_mudpack'      
      print *,'#not implemented for the moment'
      stop
      
      !call solve( poisson%poiss, E1, E2, rho)
      
    end subroutine compute_E_from_rho_2d_mudpack
  
  
  
  
end module sll_module_poisson_2d_mudpack_solver
