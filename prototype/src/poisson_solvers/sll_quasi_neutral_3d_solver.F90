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

module sll_module_quasi_neutral_3d_solver
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
!use sll_boundary_condition_descriptors
use sll_module_poisson_3d_base
use sll_module_poisson_2d_base
use sll_remapper
use sll_collective

implicit none

  type,extends(sll_poisson_3d_base) :: quasi_neutral_3d_solver     
  
  class(sll_poisson_2d_base), pointer                   :: poiss
  type(layout_3D),  pointer           :: layout_seq_x1x2 !< layout sequential in x
  type(layout_3D),  pointer           :: layout_seq_x3   !< layout sequential in y

  contains
    procedure, pass(poisson) :: initialize => &
      initialize_quasi_neutral_3d_solver
    procedure, pass(poisson) :: compute_phi_from_rho => &
      compute_phi_from_rho_3d_quasi_neutral
    procedure, pass(poisson) :: compute_E_from_rho => &
      compute_E_from_rho_3d_quasi_neutral
!    procedure, pass(poisson) :: compute_E_from_phi => &
!      compute_E_from_phi_2d_polar
      
  end type quasi_neutral_3d_solver

contains
  function new_quasi_neutral_3d_solver( &
    nc_eta1, &
    nc_eta2, &
    nc_eta3, &
    start_layout) &
    result(poisson)
      
    type(quasi_neutral_3d_solver),pointer :: poisson
    sll_int32, intent(in) :: nc_eta1
    sll_int32, intent(in) :: nc_eta2
    sll_int32, intent(in) :: nc_eta3
    type(layout_3D),  pointer           :: start_layout
    
    sll_int32 :: ierr
      
    SLL_ALLOCATE(poisson,ierr)
    call initialize_quasi_neutral_3d_solver( &
      poisson, &
      nc_eta1, &
      nc_eta2, &
      nc_eta3, &
      start_layout)    
  end function new_quasi_neutral_3d_solver
  
  
  subroutine initialize_quasi_neutral_3d_solver( &
    poisson, &
    nc_eta1, &
    nc_eta2, &
    nc_eta3, &
    start_layout )   
    class(quasi_neutral_3d_solver) :: poisson
    sll_int32, intent(in) :: nc_eta1
    sll_int32, intent(in) :: nc_eta2
    sll_int32, intent(in) :: nc_eta3
    type(layout_3D),  pointer           :: start_layout
    !sll_int32, intent(in) :: bc(2)
    !sll_int32 :: ierr
    sll_int64                        :: colsz ! collective size
    type(sll_collective_t), pointer  :: collective
    sll_int32                        :: loc_sz_x1
    sll_int32                        :: loc_sz_x2
    sll_int32                        :: loc_sz_x3
    sll_int32                        :: power2 
    sll_int32                        :: world_size 
    sll_int32                        :: nproc_x1 
    sll_int32                        :: nproc_x2
    sll_int32                        :: nproc_x3 

    ! The collective to be used is the one that comes with the given layout.
    collective => get_layout_collective( start_layout )
    colsz      = sll_get_collective_size( collective )
    world_size = sll_get_collective_size(sll_world_collective)
    power2 = int(log(real(world_size))/log(2.0))

    call compute_local_sizes_3d( &
      start_layout, &
      loc_sz_x1, &
      loc_sz_x2, &
      loc_sz_x3 )
    
    if((loc_sz_x1==nc_eta1+1).and.(loc_sz_x2==nc_eta2+1))then
      poisson%layout_seq_x1x2 => start_layout
      if(is_even(power2)) then
        nproc_x1 = 2**(power2/2)
        nproc_x2 = 2**(power2/2)
      else 
        nproc_x1 = 2**((power2-1)/2)
        nproc_x2 = 2**((power2+1)/2)
      end if
      nproc_x3 = 1
      call initialize_layout_with_distributed_3D_array( &
        nc_eta1+1, & 
        nc_eta2+1, & 
        nc_eta3+1, &
        nproc_x1, &
        nproc_x2, &
        nproc_x3, &
        poisson%layout_seq_x3 )
    else if (loc_sz_x3==nc_eta3+1) then
      poisson%layout_seq_x3 => start_layout
      poisson%layout_seq_x1x2  => new_layout_3D( sll_world_collective )
      nproc_x1 = 1
      nproc_x2 = 1
      nproc_x3 = world_size 
      call initialize_layout_with_distributed_3D_array( &
        nc_eta1+1, & 
        nc_eta2+1, & 
        nc_eta3+1, &
        nproc_x1, &
        nproc_x2, &
        nproc_x3, &
        poisson%layout_seq_x1x2 )
    else
      print *,'#Bad start_layout in initialize_quasi_neutral_3d_solver'
      print *,'#such parallelization not implemented'
      stop
    endif  


    
    
!    poisson%poiss => new_plan_poisson_polar( &
!      delta_eta,& 
!      eta1_min, &
!      nc_eta1, &
!      nc_eta2, &
!      bc, &
!      dlog_density, &
!      inv_Te)

        
  end subroutine initialize_quasi_neutral_3d_solver
  
  ! solves -\Delta phi = rho in 2d
  subroutine compute_phi_from_rho_3d_quasi_neutral( poisson, phi, rho )
    class(quasi_neutral_3d_solver) :: poisson
    sll_real64,dimension(:,:,:),intent(in) :: rho
    sll_real64,dimension(:,:,:),intent(out) :: phi

      phi = 0._f64
      print *,maxval(rho)
      
      if(.not.(associated(poisson%poiss)))then
        print *,'#poisson%poiss is not associated'
      endif

    
    !call solve( poisson%poiss, rho, phi)
    
  end subroutine compute_phi_from_rho_3d_quasi_neutral

    ! solves E = -\nabla Phi in 2d
!    subroutine compute_E_from_phi_2d_fft( poisson, phi, E1, E2 )
!      class(poisson_2d_fft_solver) :: poisson
!      sll_real64,dimension(:,:),intent(in) :: phi
!      sll_real64,dimension(:,:),intent(out) :: E1
!      sll_real64,dimension(:,:),intent(out) :: E2
!    end subroutine compute_E_from_phi_2d_fft

    ! solves E = -\nabla Phi with -\Delta phi = rho in 2d 
    subroutine compute_E_from_rho_3d_quasi_neutral( poisson, rho, E1, E2, E3 )
      class(quasi_neutral_3d_solver) :: poisson
      sll_real64,dimension(:,:,:),intent(in) :: rho
      sll_real64,dimension(:,:,:),intent(out) :: E1
      sll_real64,dimension(:,:,:),intent(out) :: E2
      sll_real64,dimension(:,:,:),intent(out) :: E3
      
      print *,'#compute_E_from_rho_3d_quasi_neutral'      
      print *,'#not implemented for the moment'
      
      
      E1 = 0._f64
      E2 = 0._f64
      E3 = 0._f64
      print *,maxval(rho)
      
      if(.not.(associated(poisson%poiss)))then
        print *,'#poisson%poiss is not associated'
      endif

      
      
      stop
      
      !call solve( poisson%poiss, E1, E2, rho)
      
    end subroutine compute_E_from_rho_3d_quasi_neutral
  
  
  
  
end module sll_module_quasi_neutral_3d_solver
