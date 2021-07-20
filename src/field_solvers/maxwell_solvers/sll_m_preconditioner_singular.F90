!> @ingroup maxwell_solvers
!> @brief
!> Module interface to solve Maxwell's equations with coordinate transformation in 3D
!> The linear systems are solved using PLAF
!> @details
!> 
!> @author
!> Benedikt Perse

module sll_m_preconditioner_singular
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_linear_solver_block, only : &
       sll_t_linear_solver_block

  use sll_m_linear_solver_abstract, only : &
       sll_t_linear_solver_abstract

  implicit none

  public :: &
       sll_t_preconditioner_singular

  private
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type, extends(sll_t_linear_solver_abstract) :: sll_t_preconditioner_singular
     type(sll_t_linear_solver_block), pointer   :: inverse_mass_fft !< block matrix solver
     sll_int32  :: n_total        !< number of Degrees of Freedom
     sll_real64, allocatable :: lumped_mass(:) !< lumped masslines 


   contains

     procedure :: create => create_3d_trafo
     procedure :: free => free_3d_trafo
     procedure :: solve_real => solve_3d_trafo
     procedure :: set_verbose
     procedure :: print_info
     procedure :: read_from_file

  end type sll_t_preconditioner_singular

contains


  subroutine create_3d_trafo( self, inverse_mass_fft, lumped_mass, n_total)
    class(sll_t_preconditioner_singular), intent( inout ) :: self      !< preconditioner
    type(sll_t_linear_solver_block), target    :: inverse_mass_fft !< block matrix solver
    sll_real64  :: lumped_mass(:) !< lumped masslines
    sll_int32, intent( in )  :: n_total !< number of Degrees of Freedom
  
    self%inverse_mass_fft => inverse_mass_fft
    self%n_total = n_total

    allocate( self%lumped_mass(n_total) )
    self%lumped_mass = lumped_mass

  end subroutine create_3d_trafo

  subroutine solve_3d_trafo(self, rhs, unknown)
    class(sll_t_preconditioner_singular), intent( inout ) :: self !< preconditioner
    sll_real64, intent(in   ) :: rhs(:) !< right hand side
    sll_real64, intent(  out) :: unknown(:) !< result
    !local variables
    sll_real64 :: scratch(self%n_total)

    scratch = rhs*self%lumped_mass
    call self%inverse_mass_fft%solve( scratch, unknown )
    unknown = unknown*self%lumped_mass

  end subroutine solve_3d_trafo

  subroutine set_verbose(self, verbose)
    class(sll_t_preconditioner_singular), intent( inout ) :: self   !< preconditioner
    logical, intent( in ) :: verbose !< logical for solver information

    self%verbose = verbose
  end subroutine set_verbose

  subroutine print_info(self)
    class(sll_t_preconditioner_singular), intent( in ) :: self !< preconditioner

  end subroutine print_info

  subroutine read_from_file(self, filename)
    class(sll_t_preconditioner_singular), intent( inout ) :: self   !< preconditioner
    character(len=*), intent( in ) :: filename !< filename
    
  end subroutine read_from_file


  subroutine free_3d_trafo( self )
    class(sll_t_preconditioner_singular), intent(inout) :: self !< preconditioner

    self%inverse_mass_fft => null()

  end subroutine free_3d_trafo


end module sll_m_preconditioner_singular
