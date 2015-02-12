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


!solves axisymmetric poisson

module sll_module_poisson_1d_polar_solver
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
!use sll_boundary_condition_descriptors
use sll_module_poisson_1d_base
!use sll_poisson_2d_polar
implicit none

  type,extends(sll_poisson_1d_base) :: poisson_1d_polar_solver     
    sll_real64 :: length
    sll_int32 :: nc_eta1
  !type(sll_plan_poisson_polar), pointer                   :: poiss
  
  
  contains
    procedure, pass(poisson) :: initialize => &
      initialize_poisson_1d_polar_solver
    procedure, pass(poisson) :: compute_phi_from_rho => &
      compute_phi_from_rho_1d_polar
    procedure, pass(poisson) :: compute_E_from_rho => &
      compute_E_from_rho_1d_polar
!    procedure, pass(poisson) :: compute_E_from_phi => &
!      compute_E_from_phi_2d_polar
      
  end type poisson_1d_polar_solver

contains
  function new_poisson_1d_polar_solver( &
    eta1_min, &
    eta1_max, &
    nc_eta1, &
    bc) &    
    result(poisson)
      
    type(poisson_1d_polar_solver),pointer :: poisson
    sll_real64, intent(in) :: eta1_min
    sll_real64, intent(in) :: eta1_max
    sll_int32, intent(in) :: nc_eta1
    sll_int32, intent(in), optional :: bc
    sll_int32 :: ierr
      
    SLL_ALLOCATE(poisson,ierr)
    call initialize_poisson_1d_polar_solver( &
      poisson, &
      eta1_min, &
      eta1_max, &
      nc_eta1, &
      bc )
    
  end function new_poisson_1d_polar_solver
  
  
  subroutine initialize_poisson_1d_polar_solver( &
    poisson, &
    eta1_min, &
    eta1_max, &
    nc_eta1, &
    bc )
    class(poisson_1d_polar_solver) :: poisson
    sll_real64, intent(in) :: eta1_min
    sll_real64, intent(in) :: eta1_max
    sll_int32, intent(in) :: nc_eta1
    sll_int32, intent(in), optional :: bc
    !sll_int32 :: ierr
    
    if(present(bc))then
      print *,'#Warning bc=',bc,'present but not used'
    endif    
    poisson%length = eta1_max-eta1_min
    poisson%nc_eta1 = nc_eta1
    
  end subroutine initialize_poisson_1d_polar_solver
  
  ! solves -\Delta phi = rho in 1d
  subroutine compute_phi_from_rho_1d_polar( poisson, phi, rho )
    class(poisson_1d_polar_solver), target :: poisson
    sll_real64,dimension(:),intent(in) :: rho
    sll_real64,dimension(:),intent(out) :: phi

    print *,'#compute_phi_from_rho_1d_polar'
    print *,'#not implemented yet'
    phi = 0._f64    
    print *,poisson%nc_eta1
    print *,maxval(rho)  
    stop
    stop
    
    
  end subroutine compute_phi_from_rho_1d_polar

    ! solves E = -\nabla Phi in 2d
!    subroutine compute_E_from_phi_2d_fft( poisson, phi, E1, E2 )
!      class(poisson_2d_fft_solver) :: poisson
!      sll_real64,dimension(:,:),intent(in) :: phi
!      sll_real64,dimension(:,:),intent(out) :: E1
!      sll_real64,dimension(:,:),intent(out) :: E2
!    end subroutine compute_E_from_phi_2d_fft

    ! solves E = -\nabla Phi with -\Delta phi = rho in 1d 
  subroutine compute_E_from_rho_1d_polar( poisson, E, rho )
    class(poisson_1d_polar_solver) :: poisson
    sll_real64,dimension(:),intent(in) :: rho
    sll_real64,dimension(:),intent(out) :: E
    sll_int32 :: N
    sll_real64 :: L
    
    N = poisson%nc_eta1
    L = poisson%length
    
    E(1:N+1) = rho(1:N+1)
    
    call poisson1dpolar(E, L, N)
      
      
  end subroutine compute_E_from_rho_1d_polar
  

  subroutine poisson1dpolar(E,L,N)
    integer,intent(in)::N
    !sll_real64,dimension(0:N),intent(inout)::E
    sll_real64, dimension(:), intent(inout) :: E
    sll_real64, intent(in) :: L
    integer :: i
    sll_real64 :: eold
    sll_real64 :: enew
    sll_real64 :: dx
    sll_real64 :: tmp
    
    !dx = L/real(N,f64)
    dx = L/(2._f64*real(N,f64))
    
    eold = E(1+N/2)*dx
    tmp = 0._f64
    E(1+N/2) = 0._f64
    do i=1,N/2
      enew = E(1+N/2+i)*dx
      tmp = (tmp-eold)*(1._f64-1._f64/real(i,f64))-enew
      eold = enew
      E(1+N/2+i)=tmp
      E(1+N/2-i)=-tmp
    enddo    

    
  end subroutine poisson1dpolar
  
  
  
end module sll_module_poisson_1d_polar_solver
