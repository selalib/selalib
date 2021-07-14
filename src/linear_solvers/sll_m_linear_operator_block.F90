!> @brief 
!> module for a block linear operator 
!> @details
!> Linear operators that can be written as blocks of other linear operators. 
!>
!>
!> Maintainer   ARA
!> Modified by Benedikt Perse	
!> Stability	stable

module sll_m_linear_operator_block
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"
  
  use sll_m_linear_operator_abstract, only: &
       sll_t_linear_operator_abstract
  implicit none

  public :: &
       sll_t_linear_operator_block
  
  private
  ! .........................................................
  !> @brief 
  !> class that contains a pointer to a linear operator 
  type  :: jrk_t_linear_operator_pointer
    class(sll_t_linear_operator_abstract), pointer :: ptr_linear_operator => null() 
  end type jrk_t_linear_operator_pointer
  ! .........................................................

  ! ..................................................
  !> @brief 
  !> class for a linear operator_block 
  type, extends(sll_t_linear_operator_abstract) :: sll_t_linear_operator_block
    sll_int32, dimension(:), allocatable :: arr_n_rows  
    !< number of rows for every block
    sll_int32, dimension(:), allocatable :: arr_n_cols  
    !< number of columns for every block

    sll_real64, dimension(:,:), allocatable :: arr_coeffs  
    !< coefficients for every block (default: 1) 

    type(jrk_t_linear_operator_pointer), dimension(:,:), allocatable :: linear_operators !< linear operators
  contains
    procedure :: create       => create_linear_operator_block
    procedure :: free         => free_linear_operator_block
    procedure :: dot => dot_linear_operator_block 
    procedure :: set          => set_linear_operator_block 
    procedure :: print_info   => print_info_linear_operator_block
  end type sll_t_linear_operator_block
  ! ..................................................

contains

  ! ............................................
  !> @brief     creates a linear operator_block 
  !>            you must sets the linop using set
  !>
  !> @param[inout] self          the current object 
  !> @param[in]    n_block_rows  number of row blocks   
  !> @param[in]    n_block_cols  number of col blocks   
  subroutine create_linear_operator_block(self, n_block_rows, n_block_cols)
  implicit none
    class(sll_t_linear_operator_block), intent(inout) :: self 
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
    call self % initialize_abstract ()
    ! ...

    ! ...
    allocate(self % arr_n_rows(n_block_rows))
    allocate(self % arr_n_cols(n_block_cols))
    allocate(self % linear_operators(n_block_rows, n_block_cols))
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
        self % linear_operators(i_block_row, i_block_col) % ptr_linear_operator => null() 
      end do
    end do
    ! ...

  end subroutine create_linear_operator_block
  ! ............................................

  ! ............................................
  !> @brief     destroys the current object 
  !>
  !> @param[inout] self   the current object 
  subroutine free_linear_operator_block(self)
  implicit none
    class(sll_t_linear_operator_block), intent(inout) :: self 
    ! local

    ! ...
    deallocate(self % arr_n_rows)
    deallocate(self % arr_n_cols)
    deallocate(self % linear_operators)
    deallocate(self % arr_coeffs)
    ! ...

  end subroutine free_linear_operator_block 
  ! ............................................

  ! ............................................
  !> @brief     sets a linear operator 
  !>
  !> @param[inout] self         the current object 
  !> @param[in]    i_block_row  the block row identifier 
  !> @param[in]    i_block_col  the block col identifier 
  !> @param[in]    linop        the associated linear operator 
  !> @param[in]    r_coeff      a real coefficient to associate to the entry block [optional] 
  subroutine set_linear_operator_block(self, i_block_row, i_block_col, linop, r_coeff)
  implicit none
    class(sll_t_linear_operator_block)           , intent(inout) :: self 
    sll_int32                                      , intent(in)    :: i_block_row
    sll_int32                                      , intent(in)    :: i_block_col
    class(sll_t_linear_operator_abstract), target, intent(in)    :: linop 
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
    self % linear_operators(i_block_row, i_block_col) % ptr_linear_operator => linop 
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
    
  end subroutine set_linear_operator_block 
  ! ............................................

  ! ............................................
  !> @brief     prints the current object 
  !>
  !> @param[in] self   the current object 
  subroutine print_info_linear_operator_block(self)
  implicit none
    class(sll_t_linear_operator_block), intent(in) :: self 
    ! local

    print *, ">>> linear_operator_block"
    call self % print_info_abstract()

    print *, "* arr_n_rows    : ", self % arr_n_rows
    print *, "* arr_n_cols    : ", self % arr_n_cols
    print *, "<<<"

  end subroutine print_info_linear_operator_block 
  ! ............................................

  ! ...................................................
  !> @brief     apply the dot operation 
  !>
  !> @param[inout] self   the current object 
  !> @param[in]    x      a real valued vector 
  !> @param[inout] y      a real valued vector 
  subroutine dot_linear_operator_block(self, x, y)
  implicit none
    class(sll_t_linear_operator_block), intent(in) :: self
    sll_real64,dimension(:), intent(in   ) :: x
    sll_real64,dimension(:), intent(  out) :: y  
    ! local 
    sll_int32 :: i
    sll_int32 :: n_block_rows
    sll_int32 :: n_block_cols 
    sll_int32 :: i_block_row
    sll_int32 :: i_block_col 
    sll_int32 :: i_begin_row 
    sll_int32 :: i_begin_col
    sll_int32 :: i_end_row 
    sll_int32 :: i_end_col
    class(sll_t_linear_operator_abstract), pointer :: ptr_linop => null()
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
    y = 0.0_f64
    i_begin_row = 1 
    do i_block_row = 1, n_block_rows
      i_end_row = i_begin_row - 1 + self % arr_n_rows(i_block_row) 
      
      ! ...
      i_begin_col = 1 
      do i_block_col = 1, n_block_cols
        i_end_col = i_begin_col - 1 + self % arr_n_cols(i_block_col) 
        
        ! ...
        ptr_linop => self % linear_operators(i_block_row, i_block_col) % ptr_linear_operator 

        if (associated(ptr_linop)) then
          call ptr_linop % dot(x(i_begin_col:i_end_col), z(i_begin_row:i_end_row))
          y(i_begin_row:i_end_row) = y(i_begin_row:i_end_row) &
                                 & + z(i_begin_row:i_end_row) * self % arr_coeffs(i_block_row, i_block_col)
        end if
        
        i_begin_col = i_begin_col + self % arr_n_cols(i_block_col)  
        ! ...
      end do
      i_begin_row = i_begin_row + self % arr_n_rows(i_block_row)  
      ! ...
    end do
    ! ...
     
  end subroutine dot_linear_operator_block
  ! ...................................................

end module sll_m_linear_operator_block
