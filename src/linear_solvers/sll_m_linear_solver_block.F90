!> @brief 
!> module for a block linear solver 
!> @details
!> Linear solvers that can be written as blocks of other linear solvers. 
!>
!>
!> Maintainer   ARA
!> Modified by Benedikt Perse	
!> Stability	stable

module sll_m_linear_solver_block
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"
#include "sll_errors.h"

  use sll_m_linear_solver_abstract, only: &
       sll_t_linear_solver_abstract

  implicit none

  public :: &
       sll_t_linear_solver_block

  private
  ! .........................................................
  !> @brief 
  !> class that contains a pointer to a linear solver 
  type :: jrk_t_linear_solver_pointer
     class(sll_t_linear_solver_abstract), pointer :: ptr_linear_solver => null() 
  end type jrk_t_linear_solver_pointer
  ! .........................................................

  ! ..................................................
  !> @brief 
  !> class for a linear solver_block 
  type, extends(sll_t_linear_solver_abstract) :: sll_t_linear_solver_block
     sll_int32, dimension(:), allocatable :: arr_n_rows  
     !< number of rows for every block
     sll_int32, dimension(:), allocatable :: arr_n_cols  
     !< number of columns for every block

     sll_real64, dimension(:,:), allocatable :: arr_coeffs  
     !< coefficients for every block (default: 1)

     sll_int32 :: n_block_rows  = 1  !< number of rows blocks
     sll_int32 :: n_block_cols  = 1  !< number of columns blocks 

     type(jrk_t_linear_solver_pointer), dimension(:,:), allocatable :: linear_solvers !< linear solvers
   contains
     procedure :: create            => create_linear_solver_block
     procedure :: set               => set_linear_solver_block
     procedure :: read_from_file    => read_from_file_linear_solver_block
     procedure :: set_verbose       => set_verbose_linear_solver_block
     procedure :: solve_real        => solve_real_linear_solver_block
     procedure :: print_info        => print_info_linear_solver_block
     procedure :: free              => free_linear_solver_block
  end type sll_t_linear_solver_block
  ! ..................................................

contains

  ! ............................................
  !> @brief     creates a linear solver_block 
  !>            you must sets the linop using set
  !>
  !> @param[inout] self          the current object 
  !> @param[in]    n_block_rows  number of row blocks   
  !> @param[in]    n_block_cols  number of col blocks   
  subroutine create_linear_solver_block(self, n_block_rows, n_block_cols)
    implicit none
    class(sll_t_linear_solver_block), intent(inout) :: self 
    sll_int32                           , intent(in) :: n_block_rows
    sll_int32                           , intent(in) :: n_block_cols 
    ! local
    sll_int32 :: i_block_row
    sll_int32 :: i_block_col 

    ! ...
    self % n_block_rows = n_block_rows 
    self % n_block_cols = n_block_cols
    ! ...

    ! ...
    !call self % initialize_abstract ()
    ! ...

    ! ...
    allocate(self % arr_n_rows(n_block_rows))
    allocate(self % arr_n_cols(n_block_cols))
    allocate(self % linear_solvers(n_block_rows, n_block_cols))
    allocate(self % arr_coeffs(n_block_rows, n_block_cols))
    ! ...

    ! ...
    self % arr_n_rows = 0
    self % arr_n_cols = 0
    self % arr_coeffs = 1.0_f64
    ! ...

    ! ...
    do i_block_col = 1, n_block_cols
       do i_block_row = 1, n_block_rows
          self % linear_solvers(i_block_row, i_block_col) % ptr_linear_solver => null() 
       end do
    end do
    ! ...

  end subroutine create_linear_solver_block
  ! ............................................



  ! ............................................
  !> @brief     sets a linear solver 
  !>
  !> @param[inout] self         the current object 
  !> @param[in]    i_block_row  the block row identifier 
  !> @param[in]    i_block_col  the block col identifier 
  !> @param[in]    linop        the associated linear solver 
  !> @param[in]    r_coeff      a real coefficient to associate to the entry block [optional] 
  subroutine set_linear_solver_block(self, i_block_row, i_block_col, linop, r_coeff)
    implicit none
    class(sll_t_linear_solver_block)           , intent(inout) :: self 
    sll_int32                                      , intent(in)    :: i_block_row
    sll_int32                                      , intent(in)    :: i_block_col
    class(sll_t_linear_solver_abstract), target, intent(in)    :: linop 
    sll_real64, optional,                        intent(in)    :: r_coeff
    ! local
    logical, parameter :: verbose = .false.
    !> \todo use linop % n_rows/cols instead of n_global_rows/cols

    ! ...
    if (present(r_coeff)) then
       self % arr_coeffs(i_block_row, i_block_col) = r_coeff 
    end if
    ! ...

    ! ...
    self % linear_solvers(i_block_row, i_block_col) % ptr_linear_solver => linop 
    ! ...

    ! ...
    self % arr_n_rows(i_block_row) = max(self % arr_n_rows(i_block_row), linop % n_global_rows)
    self % arr_n_cols(i_block_col) = max(self % arr_n_cols(i_block_col), linop % n_global_cols)
    ! ...

    ! ...
    self % n_rows        = sum(self % arr_n_rows)
    self % n_cols        = sum(self % arr_n_cols)
    self % n_global_rows = sum(self % arr_n_rows)
    self % n_global_cols = sum(self % arr_n_cols)
    ! ...

  end subroutine set_linear_solver_block
  ! ............................................

  subroutine read_from_file_linear_solver_block(self, filename)
    class(sll_t_linear_solver_block), intent( inout ) :: self
    character(len=*), intent( in ) :: filename

  end subroutine read_from_file_linear_solver_block

  subroutine set_verbose_linear_solver_block( self, verbose )
    class(sll_t_linear_solver_block), intent( inout ) :: self
    logical, intent( in ) :: verbose

    self%verbose = verbose

  end subroutine set_verbose_linear_solver_block

  ! ...................................................
  !> @brief     apply the dot operation 
  !>
  !> @param[inout] self   the current object 
  !> @param[in]    x      a real valued vector 
  !> @param[inout] y      a real valued vector 
  subroutine solve_real_linear_solver_block(self, rhs, unknown)
    implicit none
    class(sll_t_linear_solver_block), intent(in) :: self
    sll_real64,dimension(:), intent(in   ) :: rhs
    sll_real64,dimension(:), intent(  out) :: unknown  
    ! local 
    sll_int32 :: n_block_rows
    sll_int32 :: n_block_cols 
    sll_int32 :: i_block_row
    sll_int32 :: i_block_col 
    sll_int32 :: i_begin_row 
    sll_int32 :: i_begin_col
    sll_int32 :: i_end_row 
    sll_int32 :: i_end_col
    class(sll_t_linear_solver_abstract), pointer :: ptr_linop => null()
    sll_real64, dimension(:), allocatable :: z

    ! ... 
    allocate(z(self % n_global_rows))
    z = 0.0_f64
    ! ...

    ! ...
    n_block_rows = self % n_block_rows
    n_block_cols = self % n_block_cols
    ! ...

    ! ...
    unknown = 0.0_f64
    i_begin_row = 1 
    do i_block_row = 1, n_block_rows
       i_end_row = i_begin_row - 1 + self % arr_n_rows(i_block_row) 

       ! ...
       i_begin_col = 1 
       do i_block_col = 1, n_block_cols
          i_end_col = i_begin_col - 1 + self % arr_n_cols(i_block_col) 

          ! ...
          ptr_linop => self % linear_solvers(i_block_row, i_block_col) % ptr_linear_solver 

          if (associated(ptr_linop)) then
             call ptr_linop % solve(rhs(i_begin_col:i_end_col), z(i_begin_row:i_end_row))
             unknown(i_begin_row:i_end_row) = unknown(i_begin_row:i_end_row) &
                  & + z(i_begin_row:i_end_row) * self % arr_coeffs(i_block_row, i_block_col)
          end if

          i_begin_col = i_begin_col + self % arr_n_cols(i_block_col)  
          ! ...
       end do
       i_begin_row = i_begin_row + self % arr_n_rows(i_block_row)  
       ! ...
    end do
    ! ...

  end subroutine solve_real_linear_solver_block
  ! ...................................................

  ! ............................................
  !> @brief     prints the current object 
  !>
  !> @param[in] self   the current object 
  subroutine print_info_linear_solver_block(self)
    implicit none
    class(sll_t_linear_solver_block), intent(in) :: self 

  end subroutine print_info_linear_solver_block
  ! ............................................

  ! ............................................
  !> @brief     destroys the current object 
  !>
  !> @param[inout] self   the current object 
  subroutine free_linear_solver_block(self)
    implicit none
    class(sll_t_linear_solver_block), intent(inout) :: self 
    ! local

    ! ...
    deallocate(self % arr_n_rows)
    deallocate(self % arr_n_cols)
    deallocate(self % linear_solvers)
    deallocate(self % arr_coeffs)
    ! ...

  end subroutine free_linear_solver_block
  ! ............................................

end module sll_m_linear_solver_block
