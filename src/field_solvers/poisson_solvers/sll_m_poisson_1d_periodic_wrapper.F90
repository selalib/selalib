#ifndef DOXYGEN_SHOULD_SKIP_THIS
!> @ingroup poisson_solvers
module sll_m_poisson_1d_periodic_wrapper
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_poisson_1d_base, only: &
    sll_c_poisson_1d_base

  use sll_m_poisson_1d_periodic, only: &
    sll_o_new, &
    sll_t_poisson_1d_periodic, &
    sll_o_solve

  implicit none

  public :: &
    sll_f_new_poisson_1d_periodic_wrapper

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type,extends(sll_c_poisson_1d_base) :: poisson_1d_periodic_wrapper     
  
  type(sll_t_poisson_1d_periodic), pointer                   :: poiss
  
  
  contains
    procedure, pass(poisson) :: sll_o_initialize => &
      initialize_poisson_1d_periodic_wrapper
    procedure, pass(poisson) :: compute_phi_from_rho => &
      compute_phi_from_rho_1d_periodic
    procedure, pass(poisson) :: compute_E_from_rho => &
      compute_E_from_rho_1d_periodic
      
  end type poisson_1d_periodic_wrapper

contains
  function sll_f_new_poisson_1d_periodic_wrapper( &
    eta1_min, &
    eta1_max, &
    nc_eta1) &
    result(poisson)
      
    type(poisson_1d_periodic_wrapper),pointer :: poisson
    sll_real64, intent(in) :: eta1_min
    sll_real64, intent(in) :: eta1_max
    sll_int32, intent(in) :: nc_eta1
    sll_int32 :: ierr
      
    SLL_ALLOCATE(poisson,ierr)
    call initialize_poisson_1d_periodic_wrapper( &
      poisson, &
      eta1_min, &
      eta1_max, &
      nc_eta1)    
  end function sll_f_new_poisson_1d_periodic_wrapper
  
  
  subroutine initialize_poisson_1d_periodic_wrapper( &
    poisson, &
    eta1_min, &
    eta1_max, &
    nc_eta1)
    class(poisson_1d_periodic_wrapper) :: poisson
    sll_real64, intent(in) :: eta1_min
    sll_real64, intent(in) :: eta1_max
    sll_int32, intent(in) :: nc_eta1
    sll_int32 :: ierr

    
    poisson%poiss => sll_o_new( &
      eta1_min, &
      eta1_max, &
      nc_eta1, &
      ierr)


    
  end subroutine initialize_poisson_1d_periodic_wrapper
  
  ! solves -\Delta phi = rho in 2d
  subroutine compute_phi_from_rho_1d_periodic( poisson, phi, rho )
    class(poisson_1d_periodic_wrapper), target :: poisson
    sll_real64,dimension(:),intent(in) :: rho
    sll_real64,dimension(:),intent(out) :: phi
    
    print *,'#compute_phi_from_rho_1d_periodic'
    print *,'#not implemented yet'
    phi = 0._f64
    if(.not.(associated(poisson%poiss)))then
      print *,'#poisson%poiss not associated'
    endif
    print *,maxval(rho)  
    stop
    
  end subroutine compute_phi_from_rho_1d_periodic

  ! solves E = -\nabla Phi with -\Delta phi = rho in 2d 
  subroutine compute_E_from_rho_1d_periodic( poisson, E, rho )
    class(poisson_1d_periodic_wrapper) :: poisson
    sll_real64,dimension(:),intent(in) :: rho
    sll_real64,dimension(:),intent(out) :: E
      
    call sll_o_solve(poisson%poiss, E, rho)
           
  end subroutine compute_E_from_rho_1d_periodic
  
end module sll_m_poisson_1d_periodic_wrapper
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
