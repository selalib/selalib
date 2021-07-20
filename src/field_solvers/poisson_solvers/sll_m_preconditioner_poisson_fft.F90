  !> @ingroup poisson_solvers
  !> @brief 
  !> This module is a wrapper around the spline FEM Poisson solver for the uniform grid with periodic boundary condtions inverted with FFTs
  !> The wrapper as a type sll_t_linear_solver_abstract allows for the use of the solver as a preconditioner to a linear solver
  !> @author
  !> Katharina Kormann

module sll_m_preconditioner_poisson_fft
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_linear_solver_abstract, only : &
       sll_t_linear_solver_abstract

  use sll_m_poisson_3d_fem_fft, only : &
       sll_t_poisson_3d_fem_fft

  implicit none

  public :: &
       sll_t_preconditioner_poisson_fft

  private
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type, extends(sll_t_linear_solver_abstract) :: sll_t_preconditioner_poisson_fft
     type(sll_t_poisson_3d_fem_fft) :: poisson_solver
     sll_real64, allocatable :: jacobi_precond(:)

     sll_int32 :: n_dofs
     
   contains

     procedure :: create => create_poisson_fft
     procedure :: free => free_poisson_fft
     procedure :: solve_real => solve_poisson_fft
     procedure :: set_verbose => set_verbose_poisson_fft
     procedure :: print_info => print_info_poisson_fft
     procedure :: read_from_file => read_from_file_poisson_fft

  end type sll_t_preconditioner_poisson_fft

contains

  subroutine create_poisson_fft( self, n_dofs, degree, delta_x, jacobi_in )
    class(sll_t_preconditioner_poisson_fft), intent( inout ) :: self !< preconditioner
    sll_int32, intent( in ) :: n_dofs(3)     !< No. of degrees of freedom in each direction
    sll_int32, intent( in ) :: degree(3)     !< spline degree in each direction
    sll_real64, intent( in ) :: delta_x(3)   !< grid size in each direction
    sll_real64, intent( in ) :: jacobi_in(:)

    call self%poisson_solver%init( n_dofs, degree, delta_x  )

    allocate( self%jacobi_precond( product(n_dofs)) )
    self%jacobi_precond = jacobi_in
    self%n_dofs = product(n_dofs)
    
  end subroutine create_poisson_fft

  subroutine free_poisson_fft( self )
    class(sll_t_preconditioner_poisson_fft), intent( inout ) :: self !< preconditioner

    call self%poisson_solver%free()
    
  end subroutine free_poisson_fft

  subroutine solve_poisson_fft(self,  rhs, unknown )
    class(sll_t_preconditioner_poisson_fft), intent( inout ) :: self !< preconditioner
    sll_real64, intent(in   ) :: rhs(:) !< right hand side
    sll_real64, intent(  out) :: unknown(:) !< result

    sll_real64, allocatable :: tmp(:)

    allocate( tmp(self%n_dofs) )
    unknown = self%jacobi_precond*rhs
   ! tmp = unknown
    call self%poisson_solver%compute_phi_from_rho( unknown, tmp )
    unknown = self%jacobi_precond*tmp
   
    
  end subroutine solve_poisson_fft

  subroutine read_from_file_poisson_fft(self, filename)
    class(sll_t_preconditioner_poisson_fft), intent( inout ) :: self   !< preconditioner
    character(len=*), intent( in ) :: filename !< filename
    
  end subroutine read_from_file_poisson_fft

  subroutine set_verbose_poisson_fft(self, verbose)
    class(sll_t_preconditioner_poisson_fft), intent( inout ) :: self   !< preconditioner
    logical, intent( in ) :: verbose !< logical for solver information

    self%verbose = verbose
  end subroutine set_verbose_poisson_fft

  subroutine print_info_poisson_fft(self)
    class(sll_t_preconditioner_poisson_fft), intent( in ) :: self !< preconditioner

  end subroutine print_info_poisson_fft

  
end module sll_m_preconditioner_poisson_fft
