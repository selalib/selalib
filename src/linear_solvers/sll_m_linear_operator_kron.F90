!> @brief
!> module for a linear operator of kronecker solver
!> @details
!> Defines a linear operator so that his dot represents the action of
!> L x = kron(A, B) x
!> or
!> L x = kron(A, B, C) x
!>
!>
!> Maintainer   ARA
!> Modified by Benedikt Perse	
!> Stability	stable


module sll_m_linear_operator_kron
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use sll_m_linear_operator_abstract, only: &
       sll_t_linear_operator_abstract

  implicit none

  public :: &
       sll_t_linear_operator_kron, &
       sll_transpose, &
       sll_matrix_to_vector, &
       sll_vector_to_matrix

  private
  ! ..................................................
  !> @brief
  !> class for a linear operator
  type, extends(sll_t_linear_operator_abstract) :: sll_t_linear_operator_kron
     sll_int32 :: n_rows_a     = 0   !< number of rows for the first direction
     sll_int32 :: n_cols_a     = 0   !< number of cols for the first direction
     sll_int32 :: n_rows_b     = 0   !< number of rows for the second direction
     sll_int32 :: n_cols_b     = 0   !< number of cols for the second direction
     sll_int32 :: n_rows_c     = 0   !< number of rows for the third direction
     sll_int32 :: n_cols_c     = 0   !< number of cols for the third direction

     class(sll_t_linear_operator_abstract), pointer :: ptr_linop_a  => null() !< pointer to the first linear operator
     class(sll_t_linear_operator_abstract), pointer :: ptr_linop_b  => null() !< pointer to the second linear operator
     class(sll_t_linear_operator_abstract), pointer :: ptr_linop_c  => null() !< pointer to the third linear operator

   contains
     procedure :: create       => create_linear_operator_kron
     procedure :: free         => free_linear_operator_kron
     procedure :: dot => dot_linear_operator_kron
     procedure :: print_info   => print_info_linear_operator_kron
  end type sll_t_linear_operator_kron
  ! ..................................................

contains

  ! ............................................
  !> @brief     creates a linear operator
  !>
  !> @param[inout] self       the current object
  !> @param[in]    linop_a    a first linear_operator object (applied to second dimension)
  !> @param[in]    linop_b    a second linear_operator object (applied to second dimension)
  !> @param[in]    linop_c    a third linear_operator object (applied to third dimension) [optional]
  subroutine create_linear_operator_kron( self, &
       & linop_a, linop_b, linop_c)
    implicit none
    class(sll_t_linear_operator_kron)                      , intent(inout) :: self
    class(sll_t_linear_operator_abstract), target          , intent(in)    :: linop_a
    class(sll_t_linear_operator_abstract), target          , intent(in)    :: linop_b
    class(sll_t_linear_operator_abstract), target, optional, intent(in)    :: linop_c
    ! local
    ! ...
    self % n_rows_a     =  linop_a % n_rows
    self % n_cols_a     =  linop_a % n_cols
    self % ptr_linop_a  => linop_a

    self % n_rows_b     = linop_b % n_rows
    self % n_cols_b     = linop_b % n_cols
    self % ptr_linop_b  => linop_b

    self % n_rows  = linop_a % n_rows * linop_b % n_rows
    self % n_cols  = linop_a % n_cols * linop_b % n_cols

    self % n_block_rows  = linop_a % n_block_rows * linop_b % n_block_rows
    self % n_block_cols  = linop_a % n_block_cols * linop_b % n_block_cols

    self % n_global_rows  = linop_a % n_global_rows * linop_b % n_global_rows
    self % n_global_cols  = linop_a % n_global_cols * linop_b % n_global_cols

    if (present(linop_c)) then
       self % n_rows_c     = linop_c % n_rows
       self % n_cols_c     = linop_c % n_cols
       self % ptr_linop_c  => linop_c

       self % n_rows         = self % n_rows * linop_c % n_rows
       self % n_cols         = self % n_cols * linop_c % n_cols

       self % n_block_rows   = self % n_block_rows * linop_c % n_block_rows
       self % n_block_cols   = self % n_block_cols * linop_c % n_block_cols

       self % n_global_rows  = self % n_global_rows * linop_c % n_global_rows
       self % n_global_cols  = self % n_global_cols * linop_c % n_global_cols
    end if
    ! ...

    ! ...
    call self % initialize_abstract ()
    ! ...

  end subroutine create_linear_operator_kron
  ! ............................................

  ! ............................................
  !> @brief     destroys the current object
  !>
  !> @param[inout] self   the current object
  subroutine free_linear_operator_kron(self)
    implicit none
    class(sll_t_linear_operator_kron), intent(inout) :: self
    ! TODO
  end subroutine free_linear_operator_kron
  ! ............................................

  ! ...................................................
  !> @brief     apply the dot operation
  !>
  !> @param[inout] self   the current object
  !> @param[in]    x      a real valued vector
  !> @param[inout] y      a real valued vector
  subroutine dot_linear_operator_kron(self, x, y)
    implicit none
    class(sll_t_linear_operator_kron), intent(in) :: self
    sll_real64,dimension(:), intent(in   ) :: x
    sll_real64,dimension(:), intent(  out) :: y
    ! local

    ! ...
    if (.not. associated(self % ptr_linop_c)) then
       call dot_2_linear_operator_kron(self, self % ptr_linop_a, self % ptr_linop_b, x, y)
    else
       call dot_3_linear_operator_kron(self, self % ptr_linop_a, self % ptr_linop_b, self % ptr_linop_c, x, y)
    end if
    ! ...

  end subroutine dot_linear_operator_kron
  ! ...................................................

  ! ...................................................
  !> @brief     apply the dot operation
  !>
  !> @param[inout] self     the current object
  !> @param[in]    linop_1  first linear operator
  !> @param[in]    linop_2  second linear operator
  !> @param[in]    x        a real valued vector
  !> @param[inout] y        a real valued vector
  subroutine dot_2_linear_operator_kron(self, linop_1, linop_2, x, y)
    implicit none
    class(sll_t_linear_operator_kron),     intent(in)    :: self
    class(sll_t_linear_operator_abstract), intent(in)    :: linop_1
    class(sll_t_linear_operator_abstract), intent(in)    :: linop_2
    sll_real64, dimension(:),            intent(in)    :: x
    sll_real64, dimension(:),            intent(inout) :: y
    ! local
    sll_int32 :: i
    sll_int32 :: n_rows_a
    sll_int32 :: n_rows_b
    sll_int32 :: n_cols_a
    sll_int32 :: n_cols_b

    sll_real64,dimension(:,:), allocatable :: x_matrix
    sll_real64,dimension(:,:), allocatable :: x_matrix_trans
    sll_real64,dimension(:,:), allocatable :: x_prime
    sll_real64,dimension(:,:), allocatable :: x_prime_trans
    sll_real64,dimension(:,:), allocatable :: y_matrix

    ! ...
    n_cols_a = linop_1 % n_cols
    n_rows_a = linop_1 % n_rows
    n_cols_b = linop_2 % n_cols
    n_rows_b = linop_2 % n_rows
    ! ...

    ! ...
    allocate(x_matrix(n_cols_b, n_cols_a))
    allocate(x_matrix_trans(n_cols_a, n_cols_b))
    allocate(x_prime(n_cols_b, n_rows_a))
    allocate(x_prime_trans(n_rows_a, n_cols_b))
    allocate(y_matrix(n_rows_b, n_rows_a))

    x_matrix = 0.0_f64

    ! ...
    call sll_vector_to_matrix(x,x_matrix)

    call sll_transpose(x_matrix, x_matrix_trans)

    do i=1,n_cols_b
       call linop_1 % dot(x_matrix_trans(:,i), x_prime_trans(:,i))
    enddo

    call sll_transpose(x_prime_trans, x_prime)

    do i=1,n_rows_a
       call linop_2 % dot(x_prime(:,i), y_matrix(:,i))
    enddo

    call sll_matrix_to_vector(y_matrix, y)


    ! ...
    deallocate(y_matrix)
    deallocate(x_matrix_trans)
    deallocate(x_prime)
    deallocate(x_prime_trans)
    deallocate(x_matrix)
    ! ...

  end subroutine dot_2_linear_operator_kron
  ! ...................................................

  ! ...................................................
  !> @brief     apply the dot operation
  !>
  !> @param[inout] self     the current object
  !> @param[in]    linop_1  first linear operator
  !> @param[in]    linop_2  second linear operator
  !> @param[in]    linop_3  third linear operator
  !> @param[in]    x        a real valued vector
  !> @param[inout] y        a real valued vector
  subroutine dot_3_linear_operator_kron(self, linop_1, linop_2, linop_3, x, y)
    implicit none
    class(sll_t_linear_operator_kron),     intent(in)    :: self
    class(sll_t_linear_operator_abstract), intent(in)    :: linop_1
    class(sll_t_linear_operator_abstract), intent(in)    :: linop_2
    class(sll_t_linear_operator_abstract), intent(in)    :: linop_3
    sll_real64, dimension(:),            intent(in)    :: x
    sll_real64, dimension(:),            intent(inout) :: y
    ! local
    sll_int32 :: i
    sll_int32 :: n_rows_1
    sll_int32 :: n_rows_23
    sll_int32 :: n_cols_1
    sll_int32 :: n_cols_23
    sll_real64,dimension(:,:), allocatable :: x_matrix
    sll_real64,dimension(:,:), allocatable :: x_matrix_trans
    sll_real64,dimension(:,:), allocatable :: x_prime
    sll_real64,dimension(:,:), allocatable :: x_prime_trans
    sll_real64,dimension(:,:), allocatable :: y_matrix

    ! ...
    n_rows_1  =  linop_1 % n_rows
    n_rows_23 =  linop_2 % n_rows * linop_3 % n_rows
    ! ...

    ! ...
    n_cols_1  =  linop_1 % n_cols
    n_cols_23 =  linop_2 % n_cols * linop_3 % n_cols
    ! ...

    ! ...
    allocate(y_matrix(n_rows_23, n_rows_1))
    allocate(x_prime(n_cols_23, n_rows_1))
    allocate(x_prime_trans(n_rows_1, n_cols_23))
    allocate(x_matrix(n_cols_23, n_cols_1))
    allocate(x_matrix_trans(n_cols_1, n_cols_23))
    ! ...

    ! ...
    call sll_vector_to_matrix(x,x_matrix)
    call sll_transpose(x_matrix, x_matrix_trans)

    !print*, "This is the matrix from vector x", x_matrix
    do i=1,n_cols_23
       call linop_1 % dot(x_matrix_trans(:,i),x_prime_trans(:,i))
    enddo


    call sll_transpose(x_prime_trans, x_prime)

    !print*, "This is the x_prime matrix!", x_prime

    do i=1,n_rows_1
       call dot_2_linear_operator_kron( self, &
            & self % ptr_linop_b, &
            & self % ptr_linop_c, &
            & x_prime(:,i), &
            & y_matrix(:,i))
    enddo

    call sll_matrix_to_vector(y_matrix,y)
    ! ...

    ! ...
    deallocate(y_matrix)
    deallocate(x_prime_trans)
    deallocate(x_matrix)
    deallocate(x_matrix_trans)
    ! ...

  end subroutine dot_3_linear_operator_kron
  ! ...................................................


  ! ...................................................
  !> @brief      prints a linear operator
  !>
  !> @param[inout] self the current object
  subroutine print_info_linear_operator_kron(self)
    implicit none
    class(sll_t_linear_operator_kron), intent(in) :: self
    ! local

    print *, ">>>> linear_operator_kron"
    call self % print_info_abstract()

    print *, "* n_rows_a        : ", self % n_rows_a
    print *, "* n_cols_a        : ", self % n_cols_a
    print *, "* n_rows_b        : ", self % n_rows_b
    print *, "* n_cols_b        : ", self % n_cols_b
    print *, "<<<< "

  end subroutine print_info_linear_operator_kron
  ! ...................................................


  ! ..................................................
  !> @brief  returnes the transpose of a 2d array      
  !>
  !> @param[in]    x     array to transpose 
  !> @param[inout] x_t   transposed array 
  recursive subroutine sll_transpose(x, x_t)
    implicit none
    sll_real64, dimension(:,:), intent(in)  :: x
    sll_real64, dimension(:,:), intent(out) :: x_t
    ! local
    sll_int32 :: m 
    sll_int32 :: n 
    sll_int32 :: i
    sll_int32 :: j 

    ! ... x of size (m,n) and x_t of size (n,m)
    m = size(x, 1)
    n = size(x, 2)
    ! ...

    !> \todo to be optimized depending on the biggest size
    ! ...
    do i=1, n
       do j=1, m
          x_t(i,j) = x(j,i)
       end do
    end do
    ! ...

  end subroutine sll_transpose
  ! ..................................................

  ! ..................................................
  !> @brief  converts a vector to a matrix 
  !>
  !> @param[in]    vec   input vector 
  !> @param[inout] mat   tensor reprensentation of the vector 
  subroutine sll_vector_to_matrix(vec, mat)
    implicit none
    sll_real64, dimension(:)  , intent(in)    :: vec
    sll_real64, dimension(:,:), intent(inout) :: mat
    ! local
    sll_int32 :: n
    sll_int32 :: m
    sll_int32 :: j
    sll_int32 :: l

    ! ...
    n = size(mat, 1)
    m = size(mat, 2)
    l = size(vec)
    ! ...

    ! ...
    do j=1, m 
       mat(1:n, j) = vec(1+(j-1)*n:j*n)
    end do
    ! ...

  end subroutine sll_vector_to_matrix
  ! ..................................................

  ! ..................................................
  !> @brief  converts a matrix to a vector 
  !>
  !> @param[in]    mat   input matrix 
  !> @param[inout] vec   array of shape 1 representation of the matrix 
  subroutine sll_matrix_to_vector(mat, vec)
    implicit none 
    sll_real64, dimension(:,:), intent(in)    :: mat
    sll_real64, dimension(:)  , intent(inout) :: vec
    !local
    sll_int32          :: n
    sll_int32          :: m
    sll_int32          :: j

    ! ...
    n = size(mat, 1)
    m = size(mat, 2)
    ! ...

    ! ...
    do j=1, m 
       vec(1+(j-1)*n:j*n) = mat(1:n, j) 
    end do
    ! ...

  end subroutine sll_matrix_to_vector
  ! ..................................................


end module sll_m_linear_operator_kron


