!> @ingroup maxwell_solvers
!> @brief
!> Module interface to solve Maxwell's equations 
!> @details
!> 
!> @author
!> Benedikt Perse

module sll_m_preconditioner_jacobi
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_linear_operator_abstract, only : &
       sll_t_linear_operator_abstract

  use sll_m_linear_solver_abstract, only : &
       sll_t_linear_solver_abstract

  implicit none

  public :: &
       sll_t_preconditioner_jacobi

  private
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type, extends(sll_t_linear_solver_abstract) :: sll_t_preconditioner_jacobi
     sll_real64, allocatable :: diag(:) !< operator diagonal  

   contains

     procedure :: create => create_3d_trafo
     procedure :: free => free_3d_trafo
     procedure :: solve_real => solve_3d_trafo
     procedure :: set_verbose
     procedure :: print_info
     procedure :: read_from_file

  end type sll_t_preconditioner_jacobi

contains


  subroutine create_3d_trafo( self, linop)
    class(sll_t_preconditioner_jacobi), intent( inout ) :: self      !< preconditioner
    class(sll_t_linear_operator_abstract)    :: linop !< linear operator
    !local variables
    sll_int32 :: i
    sll_real64, allocatable :: scratch1(:), scratch2(:)
            
    allocate( self%diag(1:linop%n_rows) )
    allocate( scratch1(1:linop%n_rows) )
    allocate( scratch2(1:linop%n_rows) )
    
    do i = 1, linop%n_rows
       scratch1 = 0._f64
       scratch1(i) = 1._f64
       call linop%dot( scratch1, scratch2 )
       self%diag(i) = scratch2(i)
    end do
    
  end subroutine create_3d_trafo

  subroutine solve_3d_trafo(self, rhs, unknown)
    class(sll_t_preconditioner_jacobi), intent( inout ) :: self !< preconditioner
    sll_real64, intent(in   ) :: rhs(:) !< right hand side
    sll_real64, intent(  out) :: unknown(:) !< result
     
    unknown = rhs/self%diag

  end subroutine solve_3d_trafo

  subroutine set_verbose(self, verbose)
    class(sll_t_preconditioner_jacobi), intent( inout ) :: self   !< preconditioner
    logical, intent( in ) :: verbose !< logical for solver information

    self%verbose = verbose
  end subroutine set_verbose

  subroutine print_info(self)
    class(sll_t_preconditioner_jacobi), intent( in ) :: self !< preconditioner

  end subroutine print_info

  subroutine read_from_file(self, filename)
    class(sll_t_preconditioner_jacobi), intent( inout ) :: self   !< preconditioner
    character(len=*), intent( in ) :: filename !< filename
    
  end subroutine read_from_file


  subroutine free_3d_trafo( self )
    class(sll_t_preconditioner_jacobi), intent(inout) :: self !< preconditioner

  end subroutine free_3d_trafo


end module sll_m_preconditioner_jacobi
