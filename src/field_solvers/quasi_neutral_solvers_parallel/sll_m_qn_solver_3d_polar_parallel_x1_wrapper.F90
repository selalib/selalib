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

!> wrapper around sll_m_qn_solver_3d_polar_parallel_x1
!> for abstract class flexibility purposes
!> here it is a specification of 3d poisson solver

module sll_m_qn_solver_3d_polar_parallel_x1_wrapper
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_poisson_3d_base, only: &
    sll_c_poisson_3d_base

  use sll_m_qn_solver_3d_polar_parallel_x1, only: &
    sll_o_new, &
    sll_t_qn_solver_3d_polar_parallel_x1, &
    sll_s_solve_qns3d_polar

  use sll_m_remapper, only: &
    sll_t_layout_2d

  implicit none

  public :: &
    sll_f_new_qn_solver_3d_polar_parallel_x1_wrapper

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  type,extends(sll_c_poisson_3d_base) :: qn_solver_3d_polar_parallel_x1_wrapper     
  
  type(sll_t_qn_solver_3d_polar_parallel_x1), pointer :: poiss


  contains
   ! procedure, pass(poisson) :: initialize => &
   !   initialize_qn_solver_3d_polar_parallel_x1_abstract
    procedure, pass(poisson) :: compute_phi_from_rho => &
      compute_phi_from_rho_3d_qns_polar_par_x1
    procedure, pass(poisson) :: compute_E_from_rho => &
      compute_E_from_rho_3d_qns_polar_par_x1
!    procedure, pass(poisson) :: compute_E_from_phi => &
!      compute_E_from_phi_2d_polar
      
  end type qn_solver_3d_polar_parallel_x1_wrapper

contains
  function sll_f_new_qn_solver_3d_polar_parallel_x1_wrapper( &
    layout_x1, &
    layout_x2, &
    x1_min, &
    x1_max, &
    num_cells_x1, &
    num_cells_x2, &
    num_cells_x3, &
    bc_rmin, &
    bc_rmax, &
    dlog_density, &
    inv_Te) &
    result(poisson)
    
    type(sll_t_layout_2d), pointer         :: layout_x1 !< sequential in x1 direction
    type(sll_t_layout_2d), pointer         :: layout_x2 !< sequential in x2 direction
    sll_real64          , intent(in) :: x1_min    !< rmin
    sll_real64          , intent(in) :: x1_max    !< rmax
    sll_int32           , intent(in) :: num_cells_x1      !< number of cells radial
    sll_int32           , intent(in) :: num_cells_x2  !< number of cells angular
    sll_int32           , intent(in) :: num_cells_x3  !< number of cells in x3 direction
    sll_int32,  optional, intent(in) :: bc_rmin !< radial boundary conditions
    sll_int32,  optional, intent(in) :: bc_rmax !< radial boundary conditions
    sll_real64, optional, intent(in) :: dlog_density(:) !< for quasi neutral solver
    sll_real64, optional, intent(in) :: inv_Te(:)    !< for quasi neutral solver
    
    type(qn_solver_3d_polar_parallel_x1_wrapper), pointer :: poisson

    sll_int32 :: ierr

    SLL_ALLOCATE(poisson,ierr)
    call initialize_qn_solver_3d_polar_parallel_x1_wrapper( &
      poisson, &
      layout_x1, &
      layout_x2, &
      x1_min, &
      x1_max, &
      num_cells_x1, &
      num_cells_x2, &
      num_cells_x3, &
      bc_rmin, &
      bc_rmax, &
      dlog_density, &
      inv_Te)

  end function sll_f_new_qn_solver_3d_polar_parallel_x1_wrapper
  
  ! TODO: check if this can be made a type-bound method
  subroutine initialize_qn_solver_3d_polar_parallel_x1_wrapper( &
    poisson, &
    layout_x1, &
    layout_x2, &
    x1_min, &
    x1_max, &
    num_cells_x1, &
    num_cells_x2, &
    num_cells_x3, &
    bc_rmin, &
    bc_rmax, &
    dlog_density, &
    inv_Te)

    type(qn_solver_3d_polar_parallel_x1_wrapper), intent(inout)  :: poisson     !< Poisson solver class
    type(sll_t_layout_2d), pointer         :: layout_x1 !< sequential in x1 direction
    type(sll_t_layout_2d), pointer         :: layout_x2 !< sequential in x2 direction
    sll_real64          , intent(in) :: x1_min    !< rmin
    sll_real64          , intent(in) :: x1_max    !< rmax
    sll_int32           , intent(in) :: num_cells_x1      !< number of cells radial
    sll_int32           , intent(in) :: num_cells_x2  !< number of cells angular
    sll_int32           , intent(in) :: num_cells_x3  !< number of cells in x3 direction
    sll_int32,  optional, intent(in) :: bc_rmin !< radial boundary conditions
    sll_int32,  optional, intent(in) :: bc_rmax !< radial boundary conditions
    sll_real64, optional, intent(in) :: dlog_density(:) !< for quasi neutral solver
    sll_real64, optional, intent(in) :: inv_Te(:)       !< for quasi neutral solver
    !local variables
    
    poisson%poiss => sll_o_new( &
      layout_x1, &
      layout_x2, &
      x1_min, &
      x1_max, &
      num_cells_x1, &
      num_cells_x2, &
      num_cells_x3, &
      bc_rmin, &
      bc_rmax, &
      dlog_density, &
      inv_Te)
          
  end subroutine initialize_qn_solver_3d_polar_parallel_x1_wrapper
  
  subroutine compute_phi_from_rho_3d_qns_polar_par_x1( poisson, phi, rho )
    class(qn_solver_3d_polar_parallel_x1_wrapper), target :: poisson
    sll_real64,dimension(:,:,:),intent(in) :: rho
    sll_real64,dimension(:,:,:),intent(out) :: phi
    
    call sll_s_solve_qns3d_polar(poisson%poiss,rho,phi)
    
  end subroutine compute_phi_from_rho_3d_qns_polar_par_x1
    
  subroutine compute_E_from_rho_3d_qns_polar_par_x1( poisson, E1, E2, E3, rho)
    class(qn_solver_3d_polar_parallel_x1_wrapper) :: poisson
    sll_real64,dimension(:,:,:),intent(in) :: rho
    sll_real64,dimension(:,:,:),intent(out) :: E1
    sll_real64,dimension(:,:,:),intent(out) :: E2
    sll_real64,dimension(:,:,:),intent(out) :: E3
      
    print *,'#compute_E_from_rho_3d_qns_polar_par_x1'      
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
  end subroutine compute_E_from_rho_3d_qns_polar_par_x1









end module sll_m_qn_solver_3d_polar_parallel_x1_wrapper
