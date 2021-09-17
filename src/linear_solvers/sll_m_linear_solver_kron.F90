!> @brief 
!> module for kronecker linear solver 
!> @details
!>  
!>
!>
!> Maintainer   ARA
!> Modified by Benedikt Perse	
!> Stability	stable

module sll_m_linear_solver_kron
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"
  
  use sll_m_linear_operator_abstract, only: &
       sll_t_linear_operator_abstract
  
  use sll_m_linear_solver_abstract, only: &
       sll_t_linear_solver_abstract
  
  use sll_m_linear_solver_mgmres, only: &
       sll_t_linear_solver_mgmres
  
  use sll_m_linear_solver_cg, only: &
       sll_t_linear_solver_cg
  
  use sll_m_linear_operator_kron, only: &
       sll_t_linear_operator_kron, &
       sll_transpose, &
       sll_matrix_to_vector, &
       sll_vector_to_matrix
  
  use sll_m_matrix_abstract, only: &
       sll_t_matrix_abstract

  implicit none

  public :: &
       sll_t_linear_solver_kron

  private
  ! ..................................................
  !> @brief 
  !> class for the kronecker linear solver
  type, extends(sll_t_linear_solver_abstract) :: sll_t_linear_solver_kron
     sll_int32 :: nrowcol(3,2) !< number of rows and columns in each direction. only used in the 3d case

     class(sll_t_linear_solver_abstract), pointer :: ptr_linear_solver_a => null() !< pointer to the first solver
     class(sll_t_linear_solver_abstract), pointer :: ptr_linear_solver_b => null() !< pointer to the second solver
     class(sll_t_linear_solver_abstract), pointer :: ptr_linear_solver_c => null() !< pointer to the third solver, if given

     class(sll_t_linear_solver_abstract), private, allocatable :: p_linear_solver_a !< first solver, if created with a linear  operator
     class(sll_t_linear_solver_abstract), private, allocatable :: p_linear_solver_b !< second solver, if created with a linear  operator 
     class(sll_t_linear_solver_abstract), private, allocatable :: p_linear_solver_c !< third solver, if created with a linear  operator 
   contains
     procedure :: create        => create_linear_solver_kron
     procedure :: read_from_file => read_from_file_linear_solver_kron
     procedure :: set_verbose   => set_verbose_linear_solver_kron
     procedure :: solve_real    => solve_real_linear_solver_kron 
     procedure :: print_info     => print_info_linear_solver_kron 
     procedure :: free          => free_linear_solver_kron
  end type sll_t_linear_solver_kron
  ! ..................................................

contains

  ! ..................................................
  !> @brief     creates a kronecker linear solver 
  !> the user must provide
  !>    - either 2/3 linear solvers
  !>    - or 1 filename and 2/3 linear operators
  !>
  !> @param[inout] self               the current object 
  !> @param[in]    filnemae           file name containing the params namelist [optional] 
  !> @param[in]    linear_operator_a  1st linear operator [optional] 
  !> @param[in]    linear_operator_b  2nd linear operator [optional]
  !> @param[in]    linear_operator_c  3rd linear operator [optional] 
  subroutine create_linear_solver_kron( self, &
       & linear_operator_a, &
       & linear_operator_b, &
       & linear_operator_c, &
       & linear_solver_a, &
       & linear_solver_b, &
       & linear_solver_c, &
       & filename)
    implicit none
    class(sll_t_linear_solver_kron), target, intent(inout) :: self 
    class(sll_t_linear_operator_abstract), target, optional, intent(in) :: linear_operator_a
    class(sll_t_linear_operator_abstract), target, optional, intent(in) :: linear_operator_b
    class(sll_t_linear_operator_abstract), target, optional, intent(in) :: linear_operator_c
    class(sll_t_linear_solver_abstract), target, optional, intent(in) :: linear_solver_a
    class(sll_t_linear_solver_abstract), target, optional, intent(in) :: linear_solver_b
    class(sll_t_linear_solver_abstract), target, optional, intent(in) :: linear_solver_c
    character(len=*), optional, intent(in) :: filename
    ! local
    sll_int32, parameter :: input_file_id = 111
    character(len=256) :: solver_filename_a = "" 
    character(len=20)  :: solver_name_a = "" 
    character(len=256) :: solver_filename_b = "" 
    character(len=20)  :: solver_name_b = "" 
    character(len=256) :: solver_filename_c = "" 
    character(len=20)  :: solver_name_c = "" 
    logical :: flag_a 
    logical :: flag_b 
    logical :: flag_c 

    namelist /linear_solver/ solver_name_a, &
         & solver_filename_a, & 
         & solver_name_b, & 
         & solver_filename_b, & 
         & solver_name_c, & 
         & solver_filename_c

    ! ... the linear solvers are either given in the parameters namelist
    !     otherwise we create the lapack_lu linear solver
    if (present(filename)) then
       ! ... in opposition to the concrete solvers, 
       !     reading the file is mandatory in order 
       !     to allocate the right type for p_linear_solver 
       open(unit=input_file_id, file=trim(filename))
       read( input_file_id, linear_solver)
       close(input_file_id)
       ! ... 
    else
       solver_name_a = "lapack_lu"
       solver_name_b = "lapack_lu"
       solver_name_c = "lapack_lu"
    end if
    ! ...

    ! ...
    if (present(linear_solver_a)) then
       self % ptr_linear_solver_a => linear_solver_a 
    else
       ! ...
       if (present(linear_operator_a)) then
          ! ...
          call sll_create_linear_solver_with_name( solver_name_a, &
               & linear_operator_a, &
               & self % p_linear_solver_a, &
               & flag_a)

          if (flag_a) then
             self % ptr_linear_solver_a => self % p_linear_solver_a 
          end if
          ! ...
       else 
          stop "create_linear_solver_kron: linear_operator_a is expected"
       end if
       ! ...
    end if
    ! ...

    ! ...
    if (present(linear_solver_b)) then
       self % ptr_linear_solver_b => linear_solver_b 
    else
       ! ...
       if (present(linear_operator_b)) then
          ! ...
          call sll_create_linear_solver_with_name( solver_name_b, &
               & linear_operator_b, &
               & self % p_linear_solver_b, &
               & flag_b)

          if (flag_b) then
             self % ptr_linear_solver_b => self % p_linear_solver_b 
          end if
          ! ...
       else 
          stop "create_linear_solver_kron: linear_operator_b is expected"
       end if
       ! ...
    end if
    ! ...

    ! ...
    if (present(linear_solver_c)) then
       self % ptr_linear_solver_c => linear_solver_c 
    else
       ! ...
       if (present(linear_operator_c)) then
          ! ...
          call sll_create_linear_solver_with_name( solver_name_c, &
               & linear_operator_c, &
               & self % p_linear_solver_c, &
               & flag_c)

          if (flag_c) then
             self % ptr_linear_solver_c => self % p_linear_solver_c 
          end if
          ! ...
       end if
       ! ...
    end if
    ! ...

    ! ...
    self % nrowcol = 0 
    ! ...

    ! ... sets number of rows ans columns for the first direction
    self % nrowcol(1, 1) = self % ptr_linear_solver_a % n_global_rows 
    self % nrowcol(1, 2) = self % ptr_linear_solver_a % n_global_cols 
    ! ...

    ! ... sets number of rows ans columns for the second direction
    self % nrowcol(2, 1) = self % ptr_linear_solver_b % n_global_rows 
    self % nrowcol(2, 2) = self % ptr_linear_solver_b % n_global_cols 
    ! ...

    ! ... sets number of rows ans columns for the third direction
    if (associated(self % ptr_linear_solver_c)) then
       self % nrowcol(3, 1) = self % ptr_linear_solver_c % n_global_rows 
       self % nrowcol(3, 2) = self % ptr_linear_solver_c % n_global_cols 
    end if
    ! ...

  end subroutine create_linear_solver_kron
  ! ..................................................
 
  ! ..................................................
  !> @brief      create an linear solver using its name and a linear operator
  !>
  !> @param[in]    solver_name       string describing the linear solver to allocate 
  !> @param[in]    linear_operator   a linear operator 
  !> @param[inout] linear_solver     linear solver to be created
  !> @param[out]   flag              true if success, false otherwise 
  !> @param[in]    pc_left           a left preconditioner  [optional]
  subroutine sll_create_linear_solver_with_name( solver_name, & 
       & linear_operator, &
       & linear_solver, &
       & flag, &
       & pc_left)
    implicit none
    character(len=*)                                  , intent(in)    :: solver_name
    class(sll_t_linear_operator_abstract)  ,target           , intent(in)    :: linear_operator 
    class(sll_t_linear_solver_abstract)  , allocatable, intent(inout) :: linear_solver 
    logical                                           , intent(out)   :: flag
    class(sll_t_linear_solver_abstract)     , optional, intent(in)    :: pc_left 
    ! local
    class(sll_t_matrix_abstract), pointer :: ptr_matrix => null() 

    ! ...
    flag = .true.
    ! ...

    ! ...
    if (.not. allocated(linear_solver)) then
       call sll_allocate_linear_solver(solver_name, linear_solver, flag)
    end if
    ! ...

    ! ...
    select type (matrix => linear_operator) 
    class is (sll_t_matrix_abstract)
       ptr_matrix => matrix
    end select
    ! ...

    ! ...
    if (flag) then
       flag = .false.
       select type (solver => linear_solver)
       class is (sll_t_linear_solver_mgmres)
          call solver % create( linear_operator, &
               & pc_left=pc_left)
          flag = .true.
       class is (sll_t_linear_solver_cg)
          call solver % create( linear_operator, &
               & pc_left=pc_left)
          flag = .true.
       end select
       ! ...
    else
       stop 'sll_create_linear_solver_with_name: could not allocate the linear solver'
    end if
    ! ...

  end subroutine sll_create_linear_solver_with_name

  ! ............................................
  !> @brief     allocate the linear solver object
  !>
  !> @param[in]    solver_name   string describing the linear solver to allocate 
  !> @param[inout] linear_solver abstract object to be allocated 
  !> @param[out]   flag          true if success, false otherwise 
  subroutine sll_allocate_linear_solver(solver_name, linear_solver, flag)
    implicit none
    character(len=*)                                , intent(in)    :: solver_name
    class(sll_t_linear_solver_abstract), allocatable, intent(inout) :: linear_solver 
    logical                                         , intent(out)   :: flag 

    ! ...
    flag = .false.
    ! ...

    ! ...
    select case (solver_name)
    case ("mgmres")
       allocate(sll_t_linear_solver_mgmres::linear_solver)
       flag = .true.
    case ("cg")
       allocate(sll_t_linear_solver_cg::linear_solver)
       flag = .true.
    end select
    ! ...

  end subroutine sll_allocate_linear_solver

  ! ............................................
  !> @brief     sets the verbose for the linear solver object
  !>
  !> @param[inout] self     the current object 
  !> @param[in]    verbose  verbose flag
  subroutine set_verbose_linear_solver_kron(self, verbose)
    implicit none
    class(sll_t_linear_solver_kron), intent(inout) :: self
    logical                   , intent(in) :: verbose

    ! ...
    call self % set_verbose_abstract(verbose)
    ! ...

  end subroutine set_verbose_linear_solver_kron
  ! ..................................................

  ! ..................................................
  !> @brief     read from file 
  !>
  !> @param[inout] self                      the current object 
  !> @param[in]    filename       [optional] name of the output file
  subroutine read_from_file_linear_solver_kron(self, filename)
    implicit none
    class(sll_t_linear_solver_kron), intent(inout) :: self
    character(len=*)                , intent(in) :: filename
    ! local

  end subroutine read_from_file_linear_solver_kron
  ! ..................................................

  ! ..................................................
  !> @brief  solves the linear system with real vectors 
  !>
  !> @param[inout] self    the current object 
  !> @param[in]    rhs     the right hand side 
  !> @param[inout] unknown the solution 
  subroutine solve_real_linear_solver_kron(self, rhs, unknown)
    implicit none
    class(sll_t_linear_solver_kron), intent(inout)    :: self 
    sll_real64, dimension(:), intent(in   ) :: rhs 
    sll_real64, dimension(:), intent(  out) :: unknown
    ! local

    ! ...
    if (.not. associated(self % ptr_linear_solver_c)) then
       call solve_2_real_linear_solver_kron( self, &
            & self % ptr_linear_solver_a, &
            & self % ptr_linear_solver_b, &
            & rhs, unknown)
    else
       call solve_3_real_linear_solver_kron(self, rhs, unknown) 
    end if
    ! ...

  end subroutine solve_real_linear_solver_kron
  ! ..................................................

  ! ..................................................
  !> @brief  solves the linear system with real vectors 
  !>
  !> @param[inout] self    the current object 
  !> @param[in]    rhs     the right hand side 
  !> @param[inout] unknown the solution 
  subroutine solve_2_real_linear_solver_kron( self, &
       & linear_solver_a, linear_solver_b, &
       & rhs, unknown)
    implicit none
    class(sll_t_linear_solver_kron), intent(inout)    :: self 
    class(sll_t_linear_solver_abstract), intent(inout) :: linear_solver_a
    class(sll_t_linear_solver_abstract), intent(inout) :: linear_solver_b
    sll_real64, dimension(:), intent(in   ) :: rhs 
    sll_real64, dimension(:), intent(  out) :: unknown
    ! local
    sll_real64, dimension(:,:), allocatable :: y
    sll_real64, dimension(:,:), allocatable :: xprim
    sll_real64, dimension(:,:), allocatable :: xprim_trans
    sll_real64, dimension(:,:), allocatable :: x_sol_trans
    sll_real64, dimension(:,:), allocatable :: x_sol
    sll_int32 :: n_rows_a
    sll_int32 :: n_rows_b
    sll_int32 :: n_cols_a
    sll_int32 :: n_cols_b
    sll_int32 :: i_col_b

    ! ...
    n_cols_a = linear_solver_a % n_global_cols
    n_rows_a = linear_solver_a % n_global_rows
    n_cols_b = linear_solver_b % n_global_cols
    n_rows_b = linear_solver_b % n_global_rows
    ! ...

    ! ... TODO to be optimized, we don't need all these arrays
    allocate(y(n_rows_b, n_rows_a))
    allocate(xprim(n_cols_b, n_rows_a))
    allocate(xprim_trans(n_rows_a, n_cols_b))
    allocate(x_sol_trans(n_rows_a, n_cols_b))
    allocate(x_sol(n_cols_b, n_rows_a))
    ! ...

    ! ...
    y = 0.0_f64

    call sll_vector_to_matrix(rhs, y)
    ! ...

    ! ... solve every direction 
    xprim = 0.0_f64
    do i_col_b = 1, n_cols_b
       call linear_solver_b % solve_real( y(:, i_col_b), &
            & xprim(:, i_col_b))
    end do
    ! ...

    ! ...
    call sll_transpose(xprim, xprim_trans) 
    ! ...

    ! ...
    x_sol_trans = 0.0_f64
    do i_col_b = 1, n_cols_b
       call linear_solver_a % solve_real( xprim_trans(:, i_col_b), &
            & x_sol_trans(:, i_col_b))
    end do
    ! ...

    ! ...
    call sll_transpose(x_sol_trans, x_sol) 
    ! ...

    ! ...
    unknown = 0.0_f64
    call sll_matrix_to_vector(x_sol, unknown) 
    ! ...

    !> TODO 
    if(self%verbose) then
       print*, 'solve_2_real_linear_solver_kron, verbose = true: not yet implemented'
    end if
    !...

  end subroutine solve_2_real_linear_solver_kron
  ! ..................................................

  ! ..................................................
  !> @brief  apply the solve operation with real vectors
  !>
  !> @param[inout] self    the current object 
  !> @param[in]    rhs     the right hand side 
  !> @param[inout] unknown the solution 
  subroutine solve_3_real_linear_solver_kron(self, rhs, unknown)
    implicit none
    class(sll_t_linear_solver_kron), intent(in)    :: self !< Solver object
    sll_real64, dimension(:), intent(in   )    :: rhs !< right-hand-side vector
    sll_real64, dimension(:), intent(  out) :: unknown !< solution of the linear system

    !local
    sll_int32 :: i,j,k,istart,iend
    sll_int32 :: stride, counter, ind3d    
    ! scratch data
    sll_real64  :: scratch_nrow_a( self%nrowcol(1,1) )
    sll_real64  :: scratch_nrow_b( self%nrowcol(2,1) )
    sll_real64  :: scratch_nrow_c( self%nrowcol(3,1) )
    sll_real64  :: scratch_2_13( self%nrowcol(2,2), self%nrowcol(1,1)*self%nrowcol(3,2) )
    sll_real64  :: scratch_3_21(  self%nrowcol(3,2), self%nrowcol(2,1)*self%nrowcol(1,1) )

    stride = self%nrowcol(1,1)
    istart = 1
    iend = self%nrowcol(1,2)
    do k=1, self%nrowcol(3,2)
       do j=1, self%nrowcol(2,2)
          call self % ptr_linear_solver_a%solve ( rhs(istart:iend), scratch_nrow_a )
          do i=1, self%nrowcol(1,1)
             scratch_2_13 (j,(k-1)*stride+i) = scratch_nrow_a(i)
          end do
          istart = iend + 1
          iend = iend + self%nrowcol(1,2)
       end do
    end do

    stride = self%nrowcol(2,1)
    counter = 1
    do k=1, self%nrowcol(3,2)
       do i=1, self%nrowcol(1,1)
          call self % ptr_linear_solver_b%solve ( scratch_2_13 (:, counter), scratch_nrow_b )
          do j=1, self%nrowcol(2,1)
             scratch_3_21 (k,(i-1)*stride+j) = scratch_nrow_b(j)
          end do
          counter = counter +1
       end do
    end do


    stride = self%nrowcol(3,1)
    counter = 1
    do i=1, self%nrowcol(1,1)
       do j=1, self%nrowcol(2,1)
          call self % ptr_linear_solver_c%solve ( scratch_3_21 (:, counter), scratch_nrow_c )
          ind3d = i+(j-1)*self%nrowcol(1,1)
          do k=1, self%nrowcol(3,1)
             unknown (ind3d) = scratch_nrow_c(k)
             ind3d = ind3d + self%nrowcol(1,1)*self%nrowcol(2,1)
          end do
          counter = counter +1
       end do
    end do


    !>\todo 
    if(self%verbose) then
       print*, 'solve_3_real_linear_solver_kron, verbose = true: not yet implemented'
    end if
    !...


  end subroutine solve_3_real_linear_solver_kron
  ! ..................................................

  ! ............................................
  !> @brief      destroys a finite element cell 
  !>
  !> @param[inout] self the current object 
  subroutine print_info_linear_solver_kron(self)
    implicit none
    class(sll_t_linear_solver_kron), intent(in) :: self 
    ! local

    print *, ">>>> linear_solver_kron"
    print *, "* verbose    : ", self % verbose 
    print *, "<<<< "

  end subroutine print_info_linear_solver_kron
  ! ............................................

  ! ..................................................
  !> @brief     destroys the current object 
  !>
  !> @param[inout] self   the current object
  subroutine free_linear_solver_kron(self)
    implicit none
    class(sll_t_linear_solver_kron), intent(inout) :: self 
    ! local

    ! ...
    if (allocated(self % p_linear_solver_a)) then
       call self % p_linear_solver_a % free()
       deallocate(self % p_linear_solver_a)
       self % ptr_linear_solver_a => null()
    end if
    ! ...

    ! ...
    if (allocated(self % p_linear_solver_b)) then
       call self % p_linear_solver_b % free()
       deallocate(self % p_linear_solver_b)
       self % ptr_linear_solver_b => null()
    end if
    ! ...

    ! ...
    if (allocated(self % p_linear_solver_c)) then
       call self % p_linear_solver_c % free()
       deallocate(self % p_linear_solver_c)
       self % ptr_linear_solver_c => null()
    end if
    ! ...

  end subroutine free_linear_solver_kron
  ! ..................................................

end module sll_m_linear_solver_kron
