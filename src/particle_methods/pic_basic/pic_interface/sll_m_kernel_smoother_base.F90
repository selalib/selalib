!> @ingroup pic_interface
!> @author Katharina Kormann, IPP
!> @brief Base class for kernel smoothers for accumulation and field evaluation in PIC.
!> @details This base class gives an abstract interface to the basic functions for accumulation of charge and current densities as well as the evaluation of a function at particle positions.
module sll_m_kernel_smoother_base

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  implicit none

  public :: &
    sll_p_collocation, &
    sll_p_galerkin, &
    sll_c_kernel_smoother

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  ! Define parameters to set if Galerkin or collocation scaling should be used in accumulation routines
  sll_int32, parameter :: sll_p_galerkin = 0
  sll_int32, parameter :: sll_p_collocation = 1

  !> Basic type of a kernel smoother used for PIC simulations
  type, abstract :: sll_c_kernel_smoother
     sll_int32              :: dim
     sll_int32              :: n_dofs  !< Number of degrees of freedom of the smoothing kernels.
     sll_int32, allocatable :: n_grid(:) !< Number of grid points per dimension for use on tensor product grid based smoothing kernels.
     
   contains

     procedure(add_single), deferred        :: add_charge !> Add the contribution of one particle to the charge density
     procedure(add_update) , deferred       :: add_current_update_v !> Add contribution of pne particle to the current density and update velocity 
!     procedure(add_single), deferred        :: add_current_single !> Add the contribution of one particle to the charge density
     procedure(eval_single), deferred       :: evaluate
     procedure(eval_multiple), deferred     :: evaluate_multiple


  end type sll_c_kernel_smoother


  !---------------------------------------------------------------------------!
  abstract interface
     subroutine add_single(this, position, weight, rho_dofs) 
       use sll_m_working_precision
       import sll_c_kernel_smoother
       class (sll_c_kernel_smoother), intent( inout ) :: this !< Kernel smoother object
       sll_real64,                intent( in ) :: position(this%dim) !< Position of the particle
       sll_real64,                intent( in ) :: weight !< Weight of the particle
       sll_real64,                    intent( inout ) :: rho_dofs(this%n_dofs) !< Coefficient vector of the charge distribution

     end subroutine add_single
  end interface
  
  !---------------------------------------------------------------------------!
  abstract interface
     subroutine add_update (this, position_old, position_new, weight, qoverm, bfield_dofs, vi, j_dofs)
       use sll_m_working_precision
       import sll_c_kernel_smoother
       class(sll_c_kernel_smoother), intent(inout) :: this !< kernel smoother object
       sll_real64, intent(in) :: position_old(this%dim)
       sll_real64, intent(in) :: position_new(this%dim)
       sll_real64, intent(in) :: weight
       sll_real64, intent(in) :: qoverm
       sll_real64, intent(in) :: bfield_dofs(this%n_dofs)
       sll_real64, intent(inout) :: vi(:)
       sll_real64, intent(inout) :: j_dofs(this%n_dofs)

     end subroutine add_update
  end interface



  !---------------------------------------------------------------------------!
  abstract interface
     subroutine eval_single(this, position, field_dofs, field_value)
       use sll_m_working_precision
       import sll_c_kernel_smoother
       class (sll_c_kernel_smoother), intent( inout ) :: this !< Kernel smoother object 
       sll_real64,                intent( in ) :: position(this%dim) !< Position of the particle
       sll_real64,                    intent( in ) :: field_dofs(this%n_dofs) !< Coefficient vector for the field DoFs
       sll_real64, intent(out) :: field_value !< Value(s) of the electric fields at given position
     end subroutine eval_single
  end interface

  !---------------------------------------------------------------------------!
  abstract interface
     subroutine eval_multiple(this, position, components, field_dofs, field_value)
       use sll_m_working_precision
       import sll_c_kernel_smoother
       class (sll_c_kernel_smoother), intent( inout ) :: this !< Kernel smoother object 
       sll_real64,                intent( in ) :: position(this%dim) !< Position of the particle
       sll_int32, intent(in) :: components(:)
       sll_real64,                    intent( in ) :: field_dofs(:,:) !< Coefficient vector for the field DoFs
       sll_real64, intent(out)                             :: field_value(:) !< Value(s) of the electric fields at given position
     end subroutine eval_multiple
  end interface




end module sll_m_kernel_smoother_base
