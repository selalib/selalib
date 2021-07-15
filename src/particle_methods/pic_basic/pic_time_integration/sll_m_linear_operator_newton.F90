
module sll_m_linear_operator_newton
#include "sll_working_precision.h"
  use sll_m_linear_operator_abstract

  use sll_m_time_propagator_pic_vm_1d2v_helper, only: &
       sll_t_time_propagator_pic_vm_1d2v_helper

  implicit none

  public :: sll_t_linear_operator_newton


  private

  type, extends(sll_t_linear_operator_abstract) :: sll_t_linear_operator_newton
     type( sll_t_time_propagator_pic_vm_1d2v_helper ), pointer  :: propagator => null()
     sll_real64 :: dt
     sll_real64 :: epsilon
     sll_real64, allocatable :: g0(:)
     
     
   contains
     procedure :: create => create_newton
     procedure :: free => free_newton
     procedure :: dot => dot_mono_r2r_newton
     procedure :: print_info => print_info_newton

  end type sll_t_linear_operator_newton


contains

  subroutine create_newton( self, propagator )
    class(sll_t_linear_operator_newton), intent( inout ) :: self
    type( sll_t_time_propagator_pic_vm_1d2v_helper ), target, intent( in )  :: propagator

    self%propagator => propagator


    self%n_rows = self%propagator%kernel_smoother_0%n_dofs*2
    self%n_cols = self%propagator%kernel_smoother_0%n_dofs*2

    self%n_global_rows = self%n_rows
    self%n_global_cols = self%n_cols

    allocate( self%g0(  self%propagator%kernel_smoother_0%n_dofs ) )
    
  end subroutine create_newton

  subroutine free_newton( self )
    class( sll_t_linear_operator_newton), intent( inout ) :: self

    self%propagator => null()
    
  end subroutine free_newton

  subroutine dot_mono_r2r_newton( self, x, y )
    class(sll_t_linear_operator_newton), intent( in ) :: self
    sll_real64, intent( in    ) :: x(:)
    sll_real64, intent(   out ) :: y(:)

    sll_int32 :: ndofs
    sll_real64 :: scratch(2*self%propagator%kernel_smoother_0%n_dofs)
    
    ndofs = self%propagator%kernel_smoother_0%n_dofs

    scratch = x*self%epsilon
    call self%propagator%evaluate_residual( self%dt, scratch, y )

    y = (y - self%g0)/self%epsilon
    
    
  end subroutine dot_mono_r2r_newton

    
  subroutine print_info_newton( self )
    class(sll_t_linear_operator_newton), intent(in) :: self 
  end subroutine print_info_newton

end module sll_m_linear_operator_newton
