!> @brief 
!> module for abstract iterative linear solvers 
!> @details
!> 
!>
!>
!> Maintainer   ARA
!> Modified by Benedikt Perse	
!> Stability	stable

module sll_m_linear_solver_iter_abstract
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"
  
  use sll_m_linear_solver_abstract, only: &
       sll_t_linear_solver_abstract
  
  use sll_m_linear_operator_abstract,  only: &
       sll_t_linear_operator_abstract
  
  use sll_m_linear_operator_penalized, only: &
       sll_t_linear_operator_penalized
  
  implicit none

  public :: &
       sll_t_linear_solver_iter_abstract

  private
  ! ..................................................
  !> @brief 
  !> class for abstract iterative linear solver 
  type, abstract, extends(sll_t_linear_solver_abstract) :: sll_t_linear_solver_iter_abstract
    
    sll_int32           :: n_maxiter  = 2000    !< maximum number of iterations
    logical           :: null_space = .false. !< true if singular linear operator:w

    sll_real64 :: atol       = 1.0d-9 !< absolute tolerance
   
    sll_real64,dimension(:), allocatable :: x_0 !< for the initialization

    class(sll_t_linear_operator_abstract), pointer :: ptr_linear_operator => null() !< pointer to the used linear operator
    class(sll_t_linear_operator_abstract), allocatable :: p_linear_operator !< used for nullspace
    class(sll_t_linear_solver_abstract), pointer :: ptr_pc_left  => null() !< pointer to a left pc
  contains
    
    procedure(sll_p_set_guess_linear_solver_iter_abstract)        , deferred :: set_guess 
    procedure(sll_p_check_convergence_linear_solver_iter_abstract), deferred :: check_convergence 
    
    procedure :: compute_residual_error => compute_residual_error_linear_solver_iter_abstract
    procedure :: set_linear_operator    => set_linop_linear_solver_iter_abstract
    procedure :: set_tolerance          => set_tolerance_linear_solver_iter_abstract
    !procedure :: free_abstract          => free_abstract_linear_solver_iter_abstract
  end type sll_t_linear_solver_iter_abstract
  ! ..................................................

 
  ! ..................................................
  abstract interface
     subroutine sll_p_set_guess_linear_solver_iter_abstract(self, x_0)
       use sll_m_working_precision
       import sll_t_linear_solver_iter_abstract

       class(sll_t_linear_solver_iter_abstract), intent(inout)  :: self
       sll_real64, dimension(:), intent(in)  :: x_0 
     end subroutine sll_p_set_guess_linear_solver_iter_abstract
  end interface
  ! ..................................................

  ! ..................................................
  abstract interface
     subroutine sll_p_check_convergence_linear_solver_iter_abstract(self, i_iteration, flag, r_err, arr_err)
       use sll_m_working_precision
       import sll_t_linear_solver_iter_abstract

       class(sll_t_linear_solver_iter_abstract), intent(in)     :: self
       sll_int32,                                  intent(in)     :: i_iteration 
       logical,                                  intent(inout)  :: flag 
       sll_real64,                   optional, intent(in)     :: r_err 
       sll_real64, dimension(:),     optional, intent(in)     :: arr_err 
     end subroutine sll_p_check_convergence_linear_solver_iter_abstract
  end interface
  ! ..................................................
contains
  ! ..................................................
  !> @brief  computes the residual error for an iterative solver 
  !>
  !> @param[inout] self    the current object 
  !> @param[inout] unknown the solution 
  !> @param[in]    rhs     the right hand side 
  subroutine compute_residual_error_linear_solver_iter_abstract(self, unknown, rhs, r_err)
  implicit none
    class(sll_t_linear_solver_iter_abstract), intent(in)    :: self
    sll_real64, dimension(:),               intent(in)    :: unknown 
    sll_real64, dimension(:),               intent(in)    :: rhs 
    sll_real64,                             intent(inout) :: r_err 
    ! local
    sll_int32 :: n
    sll_real64, dimension(:), allocatable :: residu

    ! ...
    n = size(unknown, 1)
    allocate(residu(n))
    residu = 0.0_f64
    ! ...

    ! ...
    call self % ptr_linear_operator % dot(unknown, residu)
    residu = rhs - residu

    r_err = maxval(abs(residu))
    ! ...

  end subroutine compute_residual_error_linear_solver_iter_abstract
  ! ..................................................

  ! ..................................................
  !> @brief  sets a linear operator for an iterative solver 
  !>
  !> @param[inout] self                the current object 
  !> @param[in]    linear_operator     a linear operator 
  subroutine set_linop_linear_solver_iter_abstract(self, linear_operator)
  implicit none
    class(sll_t_linear_solver_iter_abstract), target, intent(inout) :: self
    class(sll_t_linear_operator_abstract),    target, intent(in)    :: linear_operator
    ! local

    ! ..............................................................
    ! creation / allocation of the nullspace linear_operator if used
    ! ..............................................................
    ! ...
    if (self % null_space) then
      ! ... we need to distinguish between two cases
      !     - first call  => allocation + creation + set ptr_linear_operator
      !     - other calls => set ptr_linear_operator 
      ! ...

      ! ...
      if (.not. allocated(self % p_linear_operator)) then
        ! ...
        allocate(sll_t_linear_operator_penalized::self % p_linear_operator)
        ! ...
        
        ! ...
        select type (p_linear_operator => self % p_linear_operator)
        class is (sll_t_linear_operator_penalized)
          !> \todo add nullspace arguments
          call p_linear_operator % create(linear_operator=linear_operator)
        end select
        ! ...
      end if 
      ! ...

      ! ...
      self % ptr_linear_operator => self % p_linear_operator
      ! ...
    else
      ! ...
      self % ptr_linear_operator => linear_operator
      ! ...
    end if 
    ! ...

  end subroutine set_linop_linear_solver_iter_abstract
  ! ..................................................
  
  ! ..................................................
  !> @brief  set absolute tolerance 
  !>
  !> @param[inout] self   the current object 
  !> @param[in]    atol   the absolute tolerance to set
  subroutine set_tolerance_linear_solver_iter_abstract(self, atol)
  implicit none
    class(sll_t_linear_solver_iter_abstract), intent(inout) :: self
    sll_real64,                        intent(in)    :: atol 
    
    ! ...
    self % atol = atol
    ! ...
    
  end subroutine set_tolerance_linear_solver_iter_abstract
  ! ..................................................

  ! ..................................................
  !> @brief  abstract free for an iterative solver 
  !>
  !> @param[inout] self                the current object 
  subroutine free_abstract_linear_solver_iter_abstract(self)
  implicit none
    class(sll_t_linear_solver_iter_abstract), intent(inout) :: self
    ! local

    ! ............................................
    ! desctruction of the linear operator 
    ! ............................................
    select type (p_linear_operator => self % p_linear_operator)
    class is (sll_t_linear_operator_penalized)
      call p_linear_operator % free()
    end select
    if (allocated(self % p_linear_operator)) then
      deallocate(self % p_linear_operator)
    end if
    ! ............................................

    ! ...
    self % ptr_linear_operator =>  null() 
    ! ...

    ! ...
!    self % is_allocated = .false.
    ! ...

  end subroutine free_abstract_linear_solver_iter_abstract
  ! ..................................................


end module sll_m_linear_solver_iter_abstract
