!> @brief 
!> module for abstract linear operator 
!> @details
!> a linear operator, is an application that associates to the x vector, the y vector. How this y vector is constructed, may be
!> specified by the user, through the subroutine "dot" 
!> the "dot" subroutine must be implemented in any extension of the abstract linear operator
!>
!>
!> Maintainer   ARA
!> Modified by Benedikt Perse	
!> Stability	stable

module sll_m_linear_operator_abstract
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  implicit none

  public :: &
       sll_t_linear_operator_abstract

  private
  ! ..................................................
  !> @brief 
  !> class for abstract linear operator 
  type, abstract :: sll_t_linear_operator_abstract
     sll_int32 :: n_rows        = 0  !< number of rows local to the processor 
     sll_int32 :: n_cols        = 0  !< number of columns local to the processor
     sll_int32 :: n_global_rows = 0  !< global number of rows, different from n_rows for distributed linear operator 
     sll_int32 :: n_global_cols = 0  !< global number of columns, different from n_cols for distributed linear operator
     sll_int32 :: n_block_rows  = 1  !< number of rows blocks
     sll_int32 :: n_block_cols  = 1  !< number of columns blocks 
     sll_int32 :: n_dof         = 1  !< number of degrees of freedom per node
     logical   :: is_allocated = .false.

     sll_int32, dimension(:), allocatable :: id_rows !< array of rows-vertices handled by the current proc
     sll_int32, dimension(:), allocatable :: id_cols !< array of cols-vertices handled by the current proc

   contains

     procedure :: initialize_abstract => initialize_linear_operator_abstract
     procedure :: print_info_abstract => print_info_linear_operator_abstract

     procedure(sll_p_dot_linear_operator_abstract), deferred :: dot
     procedure(sll_p_print_info_linear_operator_abstract), deferred :: print_info
     procedure(sll_p_free_linear_operator_abstract), deferred :: free 

  end type sll_t_linear_operator_abstract
  ! ..................................................

  ! ..................................................
  abstract interface
     subroutine sll_p_dot_linear_operator_abstract(self, x, y)
       use sll_m_working_precision
       import sll_t_linear_operator_abstract

       class(sll_t_linear_operator_abstract), intent(in)    :: self
       sll_real64, dimension(:),            intent(in   )    :: x
       sll_real64, dimension(:),            intent(  out) :: y 
     end subroutine sll_p_dot_linear_operator_abstract
  end interface
  ! ..................................................

  ! ..................................................
  abstract interface
     subroutine sll_p_print_info_linear_operator_abstract(self)
       import sll_t_linear_operator_abstract

       class(sll_t_linear_operator_abstract), intent(in) :: self
     end subroutine sll_p_print_info_linear_operator_abstract
  end interface
  ! ..................................................

  ! ..................................................
  abstract interface
     subroutine sll_p_free_linear_operator_abstract(self)
       import sll_t_linear_operator_abstract
       class(sll_t_linear_operator_abstract), intent(inout)  :: self
     end subroutine sll_p_free_linear_operator_abstract
  end interface
  ! ..................................................

contains 
  ! ........................................................
  !> @brief   initialize linear solver from a linear operator or given attributs 
  !>
  !> @param[inout] self           the current object 
  !> @param[inout] other          another object 
  !> @param[in]    n_rows         [optional] number of rows 
  !> @param[in]    n_cols         [optional] number of cols 
  !> @param[in]    n_global_rows  [optional] number of global rows 
  !> @param[in]    n_global_cols  [optional] number of global cols
  !> @param[in]    n_block_rows   [optional] number of rows blocks, default value: 1
  !> @param[in]    n_block_cols   [optional] number of cols blocks, default value: 1
  !> @param[in]    id_rows        [optional] array of rows-vertices handled by the current proc
  !> @param[in]    id_cols        [optional] array of cols-vertices handled by the current proc
  subroutine initialize_linear_operator_abstract( self, &
       & other, &
       & n_rows, n_cols, &
       & n_global_rows, n_global_cols, &
       & n_block_rows, n_block_cols, &
       & id_rows, id_cols)
    implicit none
    class(sll_t_linear_operator_abstract),           intent(inout) :: self 
    class(sll_t_linear_operator_abstract), optional, intent(in)    :: other 
    sll_int32,                               optional, intent(in)    :: n_rows
    sll_int32,                               optional, intent(in)    :: n_cols
    sll_int32,                               optional, intent(in)    :: n_global_rows
    sll_int32,                               optional, intent(in)    :: n_global_cols
    sll_int32,                               optional, intent(in)    :: n_block_rows
    sll_int32,                               optional, intent(in)    :: n_block_cols
    sll_int32, dimension(:),                 optional, intent(in)    :: id_rows 
    sll_int32, dimension(:),                 optional, intent(in)    :: id_cols 

    ! ...
    if (present(other)) then
       self % n_rows = other % n_rows
       self % n_cols = other % n_cols 

       self % n_global_rows = other % n_global_rows
       self % n_global_cols = other % n_global_cols 

       self % n_block_rows = other % n_block_rows
       self % n_block_cols = other % n_block_cols 

    else
       ! ...
       if (present(n_rows)) then
          self % n_rows = n_rows
       end if

       if (present(n_global_rows)) then  
          self % n_global_rows = n_global_rows
       end if

       if ((present(n_rows)) .and. (.not. (present(n_global_rows)))) then
          self % n_global_rows = n_rows
       end if

       if ((present(n_global_rows)) .and. (.not. (present(n_rows)))) then
          self % n_rows = n_global_rows
       end if
       ! ...

       ! ...
       if (present(n_cols)) then
          self % n_cols = n_cols
       end if

       if (present(n_global_cols)) then  
          self % n_global_cols = n_global_cols
       end if

       if ((present(n_cols)) .and. (.not. (present(n_global_cols)))) then
          self % n_global_cols = n_cols
       end if

       if ((present(n_global_cols)) .and. (.not. (present(n_cols)))) then
          self % n_cols = n_global_cols
       end if
       ! ...

       ! ...
       if ((present(n_block_rows))) then
          self % n_block_rows = n_block_rows
       end if
       ! ...

       ! ...
       if ((present(n_block_cols))) then
          self % n_block_cols = n_block_cols
       end if
       ! ...
    end if
    ! ...

    if (present(id_rows).and.present(id_cols)) then
       self % n_rows = size(id_rows)
       self % n_cols = size(id_cols)
       allocate( self % id_rows(self % n_rows))
       allocate( self % id_cols(self % n_cols))

       self % id_rows = id_rows
       self % id_cols = id_cols
    endif

    ! ...
    self % n_dof = self % n_block_rows * self % n_block_cols 
    ! ...
  end subroutine initialize_linear_operator_abstract
  ! ........................................................

  ! ........................................................
  !> @brief      prints a linear operator 
  !>
  !> @param[inout] self the current object 
  subroutine print_info_linear_operator_abstract(self)
    implicit none
    class(sll_t_linear_operator_abstract), intent(in) :: self 
    ! local

    print *, "* n_bloc_rows   : ", self % n_block_rows
    print *, "* n_bloc_cols   : ", self % n_block_cols
    print *, "* n_rows        : ", self % n_rows
    print *, "* n_cols        : ", self % n_cols    
    print *, "* n_global_rows : ", self % n_global_rows
    print *, "* n_global_cols : ", self % n_global_cols
    print *, "* n_dof         : ", self % n_dof    

  end subroutine print_info_linear_operator_abstract
  ! ........................................................

end module sll_m_linear_operator_abstract
