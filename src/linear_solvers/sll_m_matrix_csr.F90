!> @brief 
!> module for Compressed Sparse Row Matrix (CSR)
!> @details

!> \todo must use is_allocated in the create of csr matrices (probably we won't create the graph etc ...)
!>
!>
!> Maintainer   ARA
!> Modified by Benedikt Perse	
!> Stability	stable

module sll_m_matrix_csr
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use sll_m_matrix_abstract, only: &
       sll_t_matrix_abstract

#ifdef OPENMP_ENABLED
  use omp_lib
#endif


  implicit none

  public :: &
       sll_t_matrix_csr

  private
  ! ..................................................
  !> @brief 
  !> class for a csr matrix 
  type, extends(sll_t_matrix_abstract) :: sll_t_matrix_csr

     sll_int32 :: n_nnz            = 0        !< number of non zero elements
     logical :: is_allocated_ia  = .false.  !< true if ia is allocated
     logical :: is_allocated_jaa = .false.  !< true if jaa is allocated

     sll_int32, dimension(:), allocatable  :: arr_ia  !< row indicies
     sll_int32, dimension(:), allocatable  :: arr_ja  !< column indicies
     sll_real64, dimension(:,:), allocatable  :: arr_a !< real entries

   contains
     procedure :: create                   => create_matrix_csr
     procedure :: free                     => free_matrix_csr
     procedure :: dot                      => dot_matrix_csr
     procedure :: add_values               => add_values_matrix_csr
     procedure :: set_values               => set_values_matrix_csr
     procedure :: get_diagonal_default     => get_diagonal_default_matrix_csr
     procedure :: get_diagonal_block       => get_diagonal_block_matrix_csr
     procedure :: print_info               => print_info_matrix_csr
     procedure :: multiply                 => multiply_matrix_csr


  end type sll_t_matrix_csr
  ! ...................................................

contains

  ! ............................................
  !> @brief     creates a csr matrix 
  !>
  !> @param[inout] self                       the current object 
  !> @param[in]    n_rows                     number of rows 
  !> @param[in]    n_cols                     number of cols
  !> @param[in]    n_nnz                      number of nonzero elements 
  !> @param[in]    n_block_rows   [optional]  number of rows blocks, default value: 1
  !> @param[in]    n_block_cols   [optional]  number of cols blocks, default value: 1
  subroutine create_matrix_csr( self, &
       & n_rows, n_cols, n_nnz, &
       & n_block_rows, n_block_cols)
    implicit none
    !> param[inout] self : csr matrix structure
    class(sll_t_matrix_csr), intent(inout) :: self
    !> param[in] n_cols : number of columns
    sll_int32, optional, intent(in) :: n_cols
    !> param[in] n_rows : number of rows
    sll_int32, optional, intent(in) :: n_rows
    !> param[in] n_nnz : number of non zero elements
    sll_int32, optional, intent(in) :: n_nnz
    sll_int32, optional, intent(in) :: n_block_rows
    sll_int32, optional, intent(in) :: n_block_cols

    ! ...
    if (     (present(n_rows))         &
         &  .or. (present(n_cols))         &
         &  .or. (present(n_nnz)) ) then

       ! ...
       call allocate_matrix_csr( self, &
            n_rows=n_rows, &
            n_cols=n_cols, &
            n_nnz=n_nnz, &
            n_block_rows=n_block_rows, &
            n_block_cols=n_block_cols)
       ! ...

       ! ...
       self % n_global_rows = self % n_rows
       self % n_global_cols = self % n_cols
       ! ...
    elseif ((present(n_block_rows)) .or. (present(n_block_cols))) then
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

       ! ...
       call self % initialize_abstract ()
       ! ...

    else
       stop "create_matrix_csr: wrong arguments"
    end if
    ! ...

  end subroutine create_matrix_csr


  ! ............................................
  !> @brief     allocates a csr matrix 
  !>
  !> @param[inout] self                       the current object 
  !> @param[in]    n_rows                     number of rows 
  !> @param[in]    n_cols                     number of cols
  !> @param[in]    n_nnz                      number of nonzero elements 
  subroutine allocate_matrix_csr( self, &
       & n_rows, n_cols, n_nnz,          &
       & n_block_rows, n_block_cols)
    implicit none
    !> param[inout] self : csr matrix structure
    class(sll_t_matrix_csr), intent(inout) :: self
    !> param[in] n_cols : number of columns
    sll_int32, optional, intent(in) :: n_cols
    !> param[in] n_rows : number of rows
    sll_int32, optional, intent(in) :: n_rows
    !> param[in] n_nnz : number of non zero elements
    sll_int32, optional, intent(in) :: n_nnz
    sll_int32, optional, intent(in) :: n_block_rows
    sll_int32, optional, intent(in) :: n_block_cols

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

    ! ...
    call self % initialize_abstract ()
    ! ...

    if (present(n_rows)) then
       self%n_rows   = n_rows

       allocate(self%arr_ia(self%n_rows+1))
       self % is_allocated_ia = .true.
    end if

    if (present(n_cols)) then
       self%n_cols   = n_cols
    end if

    if (present(n_nnz)) then
       self%n_nnz  = n_nnz

       allocate(self%arr_ja(self%n_nnz))
       allocate(self%arr_a(self % n_dof, self%n_nnz))

       self % is_allocated     = .true.
       self % is_allocated_jaa = .true.
       self%arr_a = 0.0_f64
    end if


  end subroutine allocate_matrix_csr
  ! ...................................................

  ! ............................................
  !> @brief     destroys the current object 
  !>
  !> @param[inout] self   the current object 
  subroutine free_matrix_csr(self)
    implicit none
    class(sll_t_matrix_csr), intent(inout) :: self

    ! ...
    if (self % is_allocated_ia) then
       deallocate(self%arr_ia)

       self % is_allocated    = .false.
       self % is_allocated_ia = .false.
    end if
    ! ...

    ! ...
    if (self % is_allocated_jaa) then
       deallocate(self%arr_ja)
       deallocate(self%arr_a)

       self % is_allocated     = .false.
       self % is_allocated_jaa = .false.
    end if
    ! ...

  end subroutine free_matrix_csr
  ! ...................................................

  ! ...................................................
  !> @brief     apply the dot operation of real vector 
  !>
  !> @param[inout] self   the current object 
  !> @param[in]    x      a real valued vector 
  !> @param[out]   y      a real valued vector 
  subroutine dot_matrix_csr(self, x, y)
    implicit none
    class(sll_t_matrix_csr)  , intent(in)    :: self
    sll_real64,dimension(:), intent(in   ) :: x
    sll_real64,dimension(:), intent(  out) :: y

    call dot_vector_default(self, x, y)

  end subroutine dot_matrix_csr
  ! ...................................................


  ! ...................................................
  !> @brief     sequential dot operation 
  !>
  !> @param[inout] self   the current object 
  !> @param[in]    x      a real valued vector 
  !> @param[inout] y      a real valued vector 
  subroutine dot_vector_default(self, x, y)
    implicit none
    class(sll_t_matrix_csr)  , intent(in)  :: self
    sll_real64,dimension(:), intent(in   ) :: x
    sll_real64,dimension(:), intent(  out) :: y
    !local var
    sll_int32 :: i
    sll_int32 :: i_row
    sll_int32 :: i_col
    sll_int32 :: k_min
    sll_int32 :: k_max
    sll_int32 :: i_dof
    sll_int32 :: i_block_row
    sll_int32 :: i_block_col
    sll_int32 :: n_rows
    sll_int32 :: n_cols
    sll_int32 :: n_block_rows
    sll_int32 :: n_block_cols
    sll_int32 :: n_global_rows
    sll_int32 :: n_global_cols
    sll_real64, dimension(:,:), allocatable :: arr_x
    sll_int32 :: id
    sll_real64 :: seconds
    sll_int32 :: proc_num

    ! ...
    n_rows        = self % n_rows
    n_cols        = self % n_cols
    n_block_rows  = self % n_block_rows
    n_block_cols  = self % n_block_cols
    n_global_rows = self % n_global_rows
    n_global_cols = self % n_global_cols
    ! ...

    ! ... 
    allocate(arr_x(n_block_cols, self % n_global_cols))
    ! ...

    ! ...
    i = 0
    do i_col = 1, n_global_cols
       do i_block_col = 1, n_block_cols
          i = i + 1

          arr_x(i_block_col, i_col) = x(i)
       end do
    end do
    ! ...

    y = 0.0_f64

    ! ...
!!$#ifdef OPENMP_ENABLED
!!$      !export OMP_NUM_THREADS=2
!!$      !call omp_set_num_threads(4)
!!$      
!!$      seconds = omp_get_wtime()
!!$#endif
    !$OMP PARALLEL SHARED( y, n_block_rows, n_block_cols, arr_x, self ) PRIVATE( k_min, k_max, i_dof, i_block_row, i_block_col, i)
    !write(*,*)'Total Thread number: ', OMP_GET_NUM_THREADS()
    !write(*,*)'Thread rank: ', OMP_GET_THREAD_NUM()
    !$OMP DO 
    do i_row = 1, n_rows
       ! ...
       k_min = self % arr_ia ( i_row )
       k_max = self % arr_ia ( i_row + 1 ) - 1
       ! ...

       ! ...
       i_dof = 0
       do i_block_row = 1, n_block_rows
          do i_block_col = 1, n_block_cols
             i_dof = i_dof + 1

             i = i_block_row + (i_row - 1) * n_block_rows
             y(i) = y(i) + dot_product ( self % arr_a (i_dof, k_min : k_max ), &
                  arr_x (i_block_col, self%arr_ja ( k_min : k_max ) ) )

          end do
       end do
       ! ...
    end do
    !$OMP END DO
    !$OMP END PARALLEL
!!$#ifdef OPENMP_ENABLED
!!$      seconds = omp_get_wtime()-seconds
!!$      print*, 'time elapsed', seconds
!!$#endif
    ! ... 

  end subroutine dot_vector_default



  ! ............................................
  !> @brief     adds values to the matrix 
  !>
  !> @param[inout] self     the current object 
  !> @param[in]    i_row    row index
  !> @param[in]    i_col    col index
  !> @param[in]    arr_x    values to add
  subroutine add_values_matrix_csr(self, i_row, i_col, arr_x)
    implicit none
    class(sll_t_matrix_csr)        , intent(inout) :: self
    sll_int32                        , intent(in)    :: i_row
    sll_int32                        , intent(in)    :: i_col
    sll_real64, dimension(:), intent(in)    :: arr_x
    ! local
    sll_int32 :: i
    sll_int32 :: j
    sll_int32 :: k

    ! the current line is self%arr_ia(i_row)
    do k = self % arr_ia(i_row), self % arr_ia(i_row+1)-1
       j = self % arr_ja(k)
       if ( j == i_col ) then
          self % arr_a(:, k) = self % arr_a(:, k) + arr_x(:)
          exit
       end if
    end do

  end subroutine add_values_matrix_csr
  ! ...................................................

  ! ............................................
  !> @brief     sets values to the matrix 
  !>
  !> @param[inout] self     the current object 
  !> @param[in]    i_row    row index
  !> @param[in]    i_col    col index
  !> @param[in]    arr_x    values to set
  subroutine set_values_matrix_csr(self, i_row, i_col, arr_x)
    implicit none
    class(sll_t_matrix_csr)        , intent(inout) :: self
    sll_int32                        , intent(in)    :: i_row
    sll_int32                        , intent(in)    :: i_col
    sll_real64, dimension(:), intent(in)    :: arr_x
    ! local
    sll_int32 :: i
    sll_int32 :: j
    sll_int32 :: k

    ! the current line is self%arr_ia(i_row)
    do k = self % arr_ia(i_row), self % arr_ia(i_row+1)-1
       j = self % arr_ja(k)
       if ( j == i_col ) then
          self % arr_a(:, k) = arr_x(:)
          exit
       end if
    end do

  end subroutine set_values_matrix_csr
  ! ...................................................

  ! ............................................
  !> @brief      extracts a specified diagonal from a matrix
  !>        IMPORTANT: diag is a 1-rank array
  !> @param[in]     self    the current object 
  !> @param[inout]  diag    the  extracted diagonal (1-rank array)
  !> @param[in]     i_diag  diagonal index
  subroutine get_diagonal_default_matrix_csr(self, diag, i_diag)
    implicit none
    class(sll_t_matrix_csr)   , intent(in)    :: self
    sll_real64, dimension(:), intent(inout) :: diag
    sll_int32, optional, intent(in) :: i_diag
    !local
    sll_int32 :: i_row
    sll_int32 :: n_rows
    sll_int32 :: n_cols
    sll_int32 :: i
    sll_int32 :: i_dof
    sll_int32 :: i_block_row
    sll_int32 :: i_block_col
    sll_int32 :: n_block_rows
    sll_int32 :: n_block_cols
    sll_real64, dimension(:,:,:), allocatable :: block_diag

    ! ...
    n_block_rows = self % n_block_rows * self %n_block_rows
    n_block_cols = self % n_block_cols * self %n_block_cols
    n_rows       = self % n_rows * self %n_block_rows
    n_cols       = self % n_cols* self %n_block_cols

    ! ...

    ! ...
    diag = 0.0_f64
    ! ...

    ! ...
    allocate(block_diag(n_block_rows, n_block_cols, min(n_rows, n_cols)))
    block_diag = 0.0_f64
    ! ...

    ! ...
    call self % get_diagonal_block(block_diag, i_diag=i_diag)
    ! ...

    ! ...
    do i_row = 1, n_rows
       ! ...
       i_dof = 0
       do i_block_row = 1, n_block_rows
          do i_block_col = 1, n_block_cols
             i_dof = i_dof + 1

             if (i_block_row == i_block_col) then
                i = i_block_row + (i_row - 1) * n_block_rows
                diag(i) = block_diag(i_block_row, i_block_row, i_row)
             end if

          end do
       end do
       ! ...
    end do
    ! ...

  end subroutine get_diagonal_default_matrix_csr
  ! ...................................................

  ! ............................................
  !> @brief      extracts a specified block diagonal from a matrix
  !>        IMPORTANT: diag is a 3-rank array
  !> @param[in]     self    the current object 
  !> @param[inout]  diag    the  extracted diagonal (3-rank array)
  !> @param[in]     i_diag  diagonal index
  subroutine get_diagonal_block_matrix_csr(self, diag, i_diag)
    implicit none
    class(sll_t_matrix_csr)   , intent(in)    :: self
    sll_real64, dimension(:,:,:), intent(inout) :: diag
    sll_int32, optional, intent(in) :: i_diag
    !local
    sll_int32, parameter :: job = 0 ! don't altered the scr matrix
    sll_int32 :: n_nnz_diag
    ! poisitions in array arr_ja of the diagonal elements
    sll_int32 , dimension(:), allocatable  :: arr_ind_diag
    sll_int32 :: l_i_diag
    sll_int32 :: i_dof
    sll_int32 :: n_rows
    sll_int32 :: n_cols
    sll_int32 :: i_block_row
    sll_int32 :: i_block_col
    sll_int32 :: n_block_rows
    sll_int32 :: n_block_cols

    ! ... 0 => extart the main diagonal
    l_i_diag = 0
    if (present(i_diag)) then
       l_i_diag = i_diag
    end if
    ! ...

    ! ...
    n_block_rows = self % n_block_rows
    n_block_cols = self % n_block_cols
    n_rows       = self % n_rows
    n_cols       = self % n_cols
    ! ...

    ! ...
    diag = 0.0_f64
    ! ...

    ! ...
    allocate(arr_ind_diag(min(n_rows, n_cols)))
    ! ...

    ! ...
    i_dof = 0
    do i_block_row = 1, n_block_rows
       do i_block_col = 1, n_block_cols
          i_dof = i_dof + 1

          call getdia(n_rows, n_cols, job, &
               & self % arr_a(i_dof, :), self % arr_ja, self % arr_ia, &
               & n_nnz_diag, diag(i_block_row,i_block_col,:), arr_ind_diag, l_i_diag)

       end do
    end do
    ! ...

    ! ...
    deallocate(arr_ind_diag)
    ! ...

  end subroutine get_diagonal_block_matrix_csr
  ! ...................................................




  ! ............................................
  !> @brief     prints the current object 
  !>
  !> @param[in] self   the current object 
  subroutine print_info_matrix_csr(self)
    implicit none
    class(sll_t_matrix_csr), intent(in) :: self 
    ! local

    print *, ">>>>> matrix_csr"
    print *, "* n_bloc_rows   : ", self % n_block_rows
    print *, "* n_bloc_cols   : ", self % n_block_cols
    print *, "* n_rows        : ", self % n_rows
    print *, "* n_cols        : ", self % n_cols    
    print *, "* n_global_rows : ", self % n_global_rows
    print *, "* n_global_cols : ", self % n_global_cols
    print *, "* n_nnz         : ", self % n_nnz    
    print *, "<<<<<"

  end subroutine print_info_matrix_csr
  ! ............................................


  ! ............................................
  !> @brief performs the matrix by matrix product 
  !   mat_b = self *  mat_a 
  !> @param[in]    self   the current object 
  !> @param[in]    mat_a  csr matrix (A)
  !> @param[inout] mat_b  csr matrix (B)
  subroutine multiply_matrix_csr(self, mat_a, mat_b)
    implicit none
    class(sll_t_matrix_csr),      intent(in)    :: self
    class(sll_t_matrix_abstract), intent(in)    :: mat_a
    class(sll_t_matrix_abstract), intent(inout) :: mat_b
    ! local
    sll_real64, dimension(:), allocatable  :: arr_a
    sll_int32,      dimension(:), allocatable  :: arr_ia
    sll_int32,      dimension(:), allocatable  :: arr_ja
    sll_int32,      dimension(:), allocatable  :: arr_iw
    sll_int32 :: i
    sll_int32 :: n_rows
    sll_int32 :: n_cols
    sll_int32 :: n_nnz
    sll_int32 :: i_dof
    sll_int32 :: i_err
    sll_int32 :: n_nnz_max
    sll_int32, parameter :: i_job = 1
    sll_int32, parameter :: i_current_dof = 1

    ! ...
    if((self % n_dof)*(mat_a % n_dof).ne.1) then
       stop "multiply_matrix_csr: not implemented yet if n_dof > 1."
    endif
    ! ... 

    ! ... 
    n_rows    = self % n_rows
    n_cols    = mat_a % n_cols
    n_nnz_max = n_rows * n_cols
    ! ...

    ! ...
    allocate(arr_iw(n_cols)) 
    allocate(arr_ia(n_rows +1))
    allocate(arr_ja(n_nnz_max))
    allocate(arr_a (n_nnz_max))
    ! ...

    ! ...
    select type(mat_a)
    class is(sll_t_matrix_csr)

       call amub ( n_rows, n_cols, i_job, &
            & self % arr_a(i_current_dof,:), &
            & self % arr_ja, self % arr_ia, &
            & mat_a % arr_a(i_current_dof,:), &
            & mat_a % arr_ja, mat_a % arr_ia, &
            & arr_a, arr_ja, arr_ia, &
            & n_nnz_max, arr_iw, i_err) 

       class default
       stop "Error in multiply_matrix_csr: bad matrix_a type."
    end select
    ! ...

    ! ...
    n_nnz = arr_ia(n_rows + 1) -1

    select type(mat_b)
    class is(sll_t_matrix_csr)

       call mat_b % create( n_rows=n_rows, &
            & n_cols=n_cols, &
            & n_nnz=n_nnz)

       mat_b % arr_ia = arr_ia

       do i = 1, n_nnz
          mat_b % arr_ja(i) = arr_ja(i)
          mat_b % arr_a(i_current_dof, i)  = arr_a(i)
       end do
       class default
       stop "Error in multiply_matrix_csr: bad matrix_b type."
    end select
    ! ...

    ! ...
    deallocate(arr_iw)
    deallocate(arr_ja)
    deallocate(arr_ia)
    deallocate(arr_a)
    ! ...

  end subroutine multiply_matrix_csr
  ! ...................................................



end module sll_m_matrix_csr

