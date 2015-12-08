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

!> @ingroup poisson_solvers
!> @brief
!> Solve Poisson equation on 2d mesh with polar coordinates.
!> @details
!> solves 
!> \f[
!>      \sum_{i,j=1}^2 A_{i,j}\partial_{i,j} phi
!>       +\sum_{i=1}^2B_i\partial_i \phi
!>       +C \phi = \rho
!> \f]
!> in polar coordinates
!> this leads when \f$ A_{1,2}=A_{2,1}=0 \f$ 
!> and \f$ B_2 = 0 \f$
!> \f[
!>    A_11\partial_{1,e1}\hat{\phi}+B_1\partial_{1}\hat{\phi}+(C+A_{2,2}k^2)\hat{\phi} 
!>    = \hat{\rho}
!> \f]
!> <b> Example </b>
!> @snippet poisson_solvers/test_poisson_2d_polar.F90 example
module sll_m_poisson_2d_polar_wrapper
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_poisson_2d_base, only: &
    sll_poisson_2d_base

  use sll_m_poisson_2d_polar, only: &
    new_plan_poisson_polar, &
    poisson_solve_polar, &
    sll_plan_poisson_polar, &
    solve_poisson_polar

  implicit none

  public :: &
    new_poisson_2d_polar, &
    sll_poisson_drift_kinetic

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Classic Poisson solver
  sll_int32, parameter :: SLL_POISSON_CLASSIC = 0
  !> Poisson solver for drift kinetic simulation
  sll_int32, parameter :: SLL_POISSON_DRIFT_KINETIC = 1

  !> Poisson solver in polar coordinates
  type, extends(sll_poisson_2d_base) :: poisson_2d_polar_solver     

    !> PLEASE ADD DOCUMENTATION
    type(sll_plan_poisson_polar), pointer :: solver
    !> PLEASE ADD DOCUMENTATION
    sll_int32                             :: poisson_case
    !> PLEASE ADD DOCUMENTATION
    sll_real64,dimension(:), pointer      :: dlog_density
    !> PLEASE ADD DOCUMENTATION
    sll_real64,dimension(:), pointer      :: inv_Te

  contains

    !> PLEASE ADD DOCUMENTATION
    procedure, pass(poisson) :: initialize => initialize_poisson_2d_polar_solver
    !> solves \f$ -\Delta phi(x,y) = rho(x,y) \f$
    procedure, pass(poisson) :: compute_phi_from_rho => compute_phi_from_rho_2d_polar
    !> solves \f$ -\Delta phi(x,y) = rho(x,y) \f$ and \f$ E = \nabla  \phi \f$
    procedure, pass(poisson) :: compute_E_from_rho => compute_E_from_rho_2d_polar

  end type poisson_2d_polar_solver


contains

  !> Allocate a new Poisson solver in polar coordinates
  !> @returns a pointer to the derived type
  function new_poisson_2d_polar( &
    eta1_min, &
    eta1_max, &
    nc_eta1, &
    nc_eta2, &
    bc, &
    dlog_density, &
    inv_Te, &
    poisson_case) &     
    result(poisson)
      
    type(poisson_2d_polar_solver), pointer :: poisson
    sll_real64, intent(in)                 :: eta1_min
    sll_real64, intent(in)                 :: eta1_max
    sll_int32,  intent(in)                 :: nc_eta1
    sll_int32,  intent(in)                 :: nc_eta2
    sll_int32,  intent(in)                 :: bc(2)
    sll_real64, intent(in), optional       :: dlog_density(:)
    sll_real64, intent(in), optional       :: inv_Te(:)
    sll_int32,  optional                   :: poisson_case
    sll_int32                              :: ierr
      
    SLL_ALLOCATE(poisson,ierr)
    call initialize_poisson_2d_polar_solver( &
      poisson, &
      eta1_min, &
      eta1_max, &
      nc_eta1, &
      nc_eta2, &
      bc, &
      dlog_density, &
      inv_Te, &
      poisson_case)
    
  end function new_poisson_2d_polar
  
  subroutine initialize_poisson_2d_polar_solver( &
    poisson, &
    eta1_min, &
    eta1_max, &
    nc_eta1, &
    nc_eta2, &
    bc, &
    dlog_density, &
    inv_Te, &
    poisson_case)

    class(poisson_2d_polar_solver) :: poisson

    sll_real64, intent(in) :: eta1_min
    sll_real64, intent(in) :: eta1_max
    sll_int32,  intent(in) :: nc_eta1
    sll_int32,  intent(in) :: nc_eta2
    sll_int32,  intent(in) :: bc(2)
    sll_real64, optional   :: dlog_density(:)
    sll_real64, optional   :: inv_Te(:)
    sll_int32,  optional   :: poisson_case

    sll_int32  :: ierr
    sll_real64 :: delta_eta
    
    delta_eta = (eta1_max-eta1_min)/real(nc_eta1,f64)
    
    if(present(poisson_case)) then
      poisson%poisson_case = poisson_case  
    else   
      poisson%poisson_case = SLL_POISSON_CLASSIC
    endif
    
    select case(poisson%poisson_case)
      case (SLL_POISSON_CLASSIC)
         poisson%solver => new_plan_poisson_polar( &
          delta_eta,& 
          eta1_min, &
          nc_eta1, &
          nc_eta2, &
          bc)
     case (SLL_POISSON_DRIFT_KINETIC)    
        SLL_ALLOCATE(poisson%dlog_density(nc_eta1+1),ierr)
        SLL_ALLOCATE(poisson%inv_Te(nc_eta1+1),ierr)
        if(.not.(present(dlog_density)))then
          print *,'#dlog_density should be present in initialize_poisson_2d_polar_solver'
          stop
        endif
        if(size(dlog_density)<nc_eta1+1)then
          print *,'#Bad size for dlog_density',size(dlog_density)
          stop
        endif

        if(.not.(present(inv_Te)))then
          print *,'#dlog_density should be present in initialize_poisson_2d_polar_solver'
          stop
        endif
        if(size(inv_Te)<nc_eta1+1)then
          print *,'#Bad size for dlog_density',size(inv_Te)
          stop
        endif
        poisson%dlog_density(1:nc_eta1+1)=dlog_density(1:nc_eta1+1)
        poisson%inv_Te(1:nc_eta1+1)=inv_Te(1:nc_eta1+1)
        poisson%solver => new_plan_poisson_polar( &
          delta_eta,& 
          eta1_min, &
          nc_eta1, &
          nc_eta2, &
          bc, &
          poisson%dlog_density, &
          poisson%inv_Te)
      case default
        print *,'#bad value of poisson_case=', poisson%poisson_case
        print *,'#not implemented'
        print *,'#in initialize_poisson_2d_polar_solver'
        stop
     end select   
  end subroutine initialize_poisson_2d_polar_solver
  
  subroutine compute_phi_from_rho_2d_polar( poisson, phi, rho )
    class(poisson_2d_polar_solver), target :: poisson
    sll_real64,dimension(:,:),intent(in) :: rho
    sll_real64,dimension(:,:),intent(out) :: phi
    
    select case(poisson%poisson_case)
      case (SLL_POISSON_CLASSIC)
        call poisson_solve_polar(poisson%solver,rho,phi)            
      case (SLL_POISSON_DRIFT_KINETIC)    
        call solve_poisson_polar(poisson%solver,rho,phi)
      case default
        print *,'#bad value of poisson_case=', poisson%poisson_case
        print *,'#not implemented'
        print *,'in compute_phi_from_rho_2d_polar'
        stop
     end select   

  end subroutine compute_phi_from_rho_2d_polar

    !> solves \f$ \vec{E} = -\nabla \phi \f$ with \f$ -\Delta \phi(x,y) = rho(x,y) \f$.
  subroutine compute_E_from_rho_2d_polar( poisson, E1, E2, rho )
    class(poisson_2d_polar_solver)        :: poisson
    sll_real64,dimension(:,:),intent(in)  :: rho
    sll_real64,dimension(:,:),intent(out) :: E1
    sll_real64,dimension(:,:),intent(out) :: E2
      
    print *,'#compute_E_from_rho_2d_polar'      
    print *,'#not implemented for the moment'
      
    E1 = 0._f64
    E2 = 0._f64
    print *,maxval(rho)
      
    if(.not.(associated(poisson%solver)))then
      print *,'#poisson%poiss is not associated'
    endif

!PN WTF ???
    stop
      
  end subroutine compute_E_from_rho_2d_polar
  
end module sll_m_poisson_2d_polar_wrapper
