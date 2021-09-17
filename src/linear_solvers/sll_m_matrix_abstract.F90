!> @brief 
!> module for abstract matrix 
!> @details
!> abstract matrix object
!>
!>
!> Maintainer   ARA
!> Modified by Benedikt Perse	
!> Stability	stable

module sll_m_matrix_abstract
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use sll_m_linear_operator_abstract, only: &
       sll_t_linear_operator_abstract

  implicit none

  public :: &
       sll_t_matrix_abstract

  private
  ! ...................................................
  !> @brief 
  !> abstract class for matrix 
  type, abstract, extends(sll_t_linear_operator_abstract) :: sll_t_matrix_abstract 

   contains
     generic :: get_diagonal => get_diagonal_default, &
          &  get_diagonal_block

     procedure(sll_p_add_values_matrix_abstract)        , deferred :: add_values 
     procedure(sll_p_set_values_matrix_abstract)        , deferred :: set_values
     procedure(sll_p_get_diagonal_matrix_abstract)      , deferred :: get_diagonal_default
     procedure(sll_p_get_diagonal_block_matrix_abstract), deferred :: get_diagonal_block
     procedure(sll_p_multiply_matrix_abstract)          , deferred :: multiply

  end type sll_t_matrix_abstract
  ! ...................................................

  ! ..................................................
  abstract interface
     subroutine sll_p_add_values_matrix_abstract(self, i_row, i_col, arr_x)
       use sll_m_working_precision
       import sll_t_matrix_abstract

       class(sll_t_matrix_abstract)   , intent(inout) :: self
       sll_int32                        , intent(in)    :: i_row 
       sll_int32                        , intent(in)    :: i_col
       sll_real64, dimension(:), intent(in)    :: arr_x
     end subroutine sll_p_add_values_matrix_abstract
  end interface
  ! ..................................................

  ! ..................................................
  abstract interface
     subroutine sll_p_set_values_matrix_abstract(self, i_row, i_col, arr_x)
       use sll_m_working_precision
       import sll_t_matrix_abstract

       class(sll_t_matrix_abstract)   , intent(inout) :: self
       sll_int32                        , intent(in)    :: i_row 
       sll_int32                        , intent(in)    :: i_col
       sll_real64, dimension(:), intent(in)    :: arr_x
     end subroutine sll_p_set_values_matrix_abstract
  end interface
  ! ..................................................

  ! ..................................................
  abstract interface
     subroutine sll_p_get_diagonal_matrix_abstract(self, diag, i_diag)
       use sll_m_working_precision
       import sll_t_matrix_abstract

       class(sll_t_matrix_abstract) , intent(in)    :: self
       sll_real64, dimension(:)   , intent(inout) :: diag
       sll_int32, optional, intent(in) :: i_diag

     end subroutine sll_p_get_diagonal_matrix_abstract
  end interface
  ! ..................................................

  ! ..................................................
  abstract interface
     subroutine sll_p_get_diagonal_block_matrix_abstract(self, diag, i_diag)
       use sll_m_working_precision
       import sll_t_matrix_abstract

       class(sll_t_matrix_abstract) , intent(in)        :: self
       sll_real64, dimension(:,:,:)   , intent(inout) :: diag
       sll_int32, optional, intent(in)                    :: i_diag

     end subroutine sll_p_get_diagonal_block_matrix_abstract
  end interface
  ! ..................................................


  ! ..................................................
  abstract interface
     subroutine sll_p_multiply_matrix_abstract(self, mat_a, mat_b)
       use sll_m_working_precision
       import sll_t_matrix_abstract

       class(sll_t_matrix_abstract) , intent(in)    :: self
       class(sll_t_matrix_abstract) , intent(in)    :: mat_a
       class(sll_t_matrix_abstract) , intent(inout) :: mat_b
     end subroutine sll_p_multiply_matrix_abstract
  end interface
  ! ..................................................


end module sll_m_matrix_abstract
