!> @brief 
!> module for abstract linear solver 
!> @details
!> 
!>
!>
!> Maintainer   ARA
!> Modified by Benedikt Perse	
!> Stability	stable

module sll_m_linear_solver_abstract
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use sll_m_linear_operator_abstract, only: &
       sll_t_linear_operator_abstract

  use sll_m_linear_operator_block, only: &
       sll_t_linear_operator_block

  implicit none

  public :: &
       sll_t_linear_solver_abstract

  private
  ! ..................................................
  !> @brief 
  !> class for abstract linear solver 
  type, abstract :: sll_t_linear_solver_abstract
     sll_int32 :: n_rows         = 0       !< number of rows
     sll_int32 :: n_cols         = 0       !< number of columns

     sll_int32  :: n_global_rows = 0  !< number of rows different from n_rows for distributed linear operator 
     sll_int32  :: n_global_cols = 0  !< number of columns different from n_cols for linear operator

     sll_int32  :: n_total_rows  = 0  !< n_global_rows * n_block_rows 
     sll_int32  :: n_total_cols  = 0  !< n_global_cols * n_block_cols 

     logical :: is_allocated = .false.
     logical :: verbose  = .false. !< details output for diagnostic if true

   contains

     generic :: solve => solve_real

     procedure :: initialize_abstract  => initialize_linear_solver_abstract
     procedure :: set_verbose_abstract => set_verbose_linear_solver_abstract 

     procedure(sll_p_read_from_file_linear_solver_abstract), deferred :: read_from_file
     procedure(sll_p_set_verbose_linear_solver_abstract)   , deferred :: set_verbose
     procedure(sll_p_solve_real_linear_solver_abstract)    , deferred :: solve_real 
     procedure(sll_p_print_info_linear_solver_abstract)    , deferred :: print_info 
     procedure(sll_p_free_linear_solver_abstract)          , deferred :: free 
  end type sll_t_linear_solver_abstract
  ! ..................................................

  ! ..................................................
  abstract interface
     subroutine sll_p_solve_real_linear_solver_abstract(self, rhs, unknown)
       use sll_m_working_precision
       import sll_t_linear_solver_abstract

       class(sll_t_linear_solver_abstract), intent(inout)  :: self
       sll_real64, dimension(:)         , intent(in   ) :: rhs
       sll_real64, dimension(:)         , intent(  out) :: unknown 
     end subroutine sll_p_solve_real_linear_solver_abstract
  end interface
  ! ..................................................

  ! ..................................................
  abstract interface
     subroutine sll_p_set_verbose_linear_solver_abstract(self, verbose)
       use sll_m_working_precision
       import sll_t_linear_solver_abstract

       class(sll_t_linear_solver_abstract), intent(inout) :: self
       logical                            , intent(in)    :: verbose 
     end subroutine sll_p_set_verbose_linear_solver_abstract
  end interface
  ! ..................................................

  ! ..................................................
  abstract interface
     subroutine sll_p_read_from_file_linear_solver_abstract(self, filename)
       use sll_m_working_precision
       import sll_t_linear_solver_abstract

       class(sll_t_linear_solver_abstract), intent(inout) :: self
       character(len=*)                   , intent(in)    :: filename
     end subroutine sll_p_read_from_file_linear_solver_abstract
  end interface
  ! ..................................................

  ! ..................................................
  abstract interface
     subroutine sll_p_print_info_linear_solver_abstract(self)
       import sll_t_linear_solver_abstract

       class(sll_t_linear_solver_abstract), intent(in) :: self
     end subroutine sll_p_print_info_linear_solver_abstract
  end interface
  ! ..................................................


  ! ..................................................
  abstract interface
     subroutine sll_p_free_linear_solver_abstract(self)
       import sll_t_linear_solver_abstract
       class(sll_t_linear_solver_abstract), intent(inout)  :: self
     end subroutine sll_p_free_linear_solver_abstract
  end interface
  ! ..................................................


contains

  ! ........................................................
  !> @brief   initialize linear solver from linear operator
  !>
  !> @param[inout] self             the current object 
  !> @param[in]    linear_operator  a linear operator 
  subroutine initialize_linear_solver_abstract(self, linear_operator)
    implicit none
    class(sll_t_linear_solver_abstract)          , intent(inout) :: self 
    class(sll_t_linear_operator_abstract), target, intent(in)    :: linear_operator 

    self % n_rows = linear_operator % n_rows
    self % n_cols = linear_operator % n_cols

    self % n_global_rows = linear_operator % n_global_rows
    self % n_global_cols = linear_operator % n_global_cols

    select type (linear_operator) 
    class is (sll_t_linear_operator_block)
       self % n_total_rows =  self % n_global_rows
       self % n_total_cols =  self % n_global_cols 
       class default
       self % n_total_rows =  self % n_global_rows * linear_operator % n_block_rows
       self % n_total_cols =  self % n_global_cols * linear_operator % n_block_cols
    end select

  end subroutine initialize_linear_solver_abstract
  ! ........................................................

  ! ............................................
  !> @brief     sets the verbose for the linear solver object
  !>
  !> @param[inout] self     the current object 
  !> @param[in]    verbose  verbose flag
  subroutine set_verbose_linear_solver_abstract(self, verbose)
    implicit none
    class(sll_t_linear_solver_abstract), intent(inout) :: self
    logical                            , intent(in) :: verbose

    ! ...
    self % verbose = verbose
    ! ...

  end subroutine set_verbose_linear_solver_abstract
  ! ..................................................

end module sll_m_linear_solver_abstract
