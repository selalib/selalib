module sll_module_poisson_2d_periodic_solver
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
!use sll_boundary_condition_descriptors
use sll_module_poisson_2d_base
#ifdef FFTW
use sll_poisson_2d_periodic_fftw
#else
use sll_poisson_2d_periodic_fftpack
#endif
implicit none


  type,extends(sll_poisson_2d_base) :: poisson_2d_periodic_solver     
  
#ifdef FFTW
  type(poisson_2d_periodic_fftw), pointer                   :: poiss
#else
  type(poisson_2d_periodic_fftpack), pointer                   :: poiss
#endif
  
  
  contains
    procedure, pass(poisson) :: initialize => &
      initialize_poisson_2d_periodic_solver
    procedure, pass(poisson) :: compute_phi_from_rho => &
      compute_phi_from_rho_2d_periodic
    procedure, pass(poisson) :: compute_E_from_rho => &
      compute_E_from_rho_2d_periodic
!    procedure, pass(poisson) :: compute_E_from_phi => &
!      compute_E_from_phi_2d_polar
      
  end type poisson_2d_periodic_solver

contains

  function new_poisson_2d_periodic_solver( &
    eta1_min, &
    eta1_max, &
    nc_eta1, &
    eta2_min, &
    eta2_max, &
    nc_eta2) &     
    result(poisson)
      
    type(poisson_2d_periodic_solver),pointer :: poisson
    sll_real64, intent(in) :: eta1_min
    sll_real64, intent(in) :: eta1_max
    sll_int32, intent(in) :: nc_eta1
    sll_real64, intent(in) :: eta2_min
    sll_real64, intent(in) :: eta2_max
    sll_int32, intent(in) :: nc_eta2
    sll_int32 :: ierr
      
    SLL_ALLOCATE(poisson,ierr)
    call initialize_poisson_2d_periodic_solver( &
      poisson, &
      eta1_min, &
      eta1_max, &
      nc_eta1, &
      eta2_min, &
      eta2_max, &
      nc_eta2)     
    
  end function new_poisson_2d_periodic_solver
  
  
  subroutine initialize_poisson_2d_periodic_solver( &
    poisson, &
    eta1_min, &
    eta1_max, &
    nc_eta1, &
    eta2_min, &
    eta2_max, &
    nc_eta2)     
    class(poisson_2d_periodic_solver) :: poisson
    sll_real64, intent(in) :: eta1_min
    sll_real64, intent(in) :: eta1_max
    sll_int32, intent(in) :: nc_eta1
    sll_real64, intent(in) :: eta2_min
    sll_real64, intent(in) :: eta2_max
    sll_int32, intent(in) :: nc_eta2
    sll_int32 :: ierr

    
    poisson%poiss => new( &
      eta1_min, &
      eta1_max, &
      nc_eta1, &
      eta2_min, &
      eta2_max, &
      nc_eta2, &
      ierr)


    
  end subroutine initialize_poisson_2d_periodic_solver
  
  ! solves -\Delta phi = rho in 2d
  subroutine compute_phi_from_rho_2d_periodic( poisson, phi, rho )
    class(poisson_2d_periodic_solver), target :: poisson
    sll_real64,dimension(:,:),intent(in) :: rho
    sll_real64,dimension(:,:),intent(out) :: phi
    
    call solve(poisson%poiss, phi, rho)
    
    
  end subroutine compute_phi_from_rho_2d_periodic

    ! solves E = -\nabla Phi in 2d
!    subroutine compute_E_from_phi_2d_fft( poisson, phi, E1, E2 )
!      class(poisson_2d_fft_solver) :: poisson
!      sll_real64,dimension(:,:),intent(in) :: phi
!      sll_real64,dimension(:,:),intent(out) :: E1
!      sll_real64,dimension(:,:),intent(out) :: E2
!    end subroutine compute_E_from_phi_2d_fft

    ! solves E = -\nabla Phi with -\Delta phi = rho in 2d 
  subroutine compute_E_from_rho_2d_periodic( poisson, E1, E2, rho )
    class(poisson_2d_periodic_solver) :: poisson
    sll_real64,dimension(:,:),intent(in) :: rho
    sll_real64,dimension(:,:),intent(out) :: E1
    sll_real64,dimension(:,:),intent(out) :: E2
      
    call solve(poisson%poiss, E1, E2, rho)
           
  end subroutine compute_E_from_rho_2d_periodic
  
  
  
  
end module sll_module_poisson_2d_periodic_solver
