!> @brief 
!> module for a sequential gmres
!> @details
!>  
!>
!>
!> Maintainer   ARA
!> Modified by Benedikt Perse	
!> Stability	stable

module sll_m_linear_solver_mgmres
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use sll_m_linear_solver_iter_abstract, only: &
       sll_t_linear_solver_iter_abstract

  use sll_m_linear_operator_abstract, only: &
       sll_t_linear_operator_abstract

  use sll_m_linear_solver_abstract, only: &
       sll_t_linear_solver_abstract

  implicit none

  public :: &
       sll_t_linear_solver_mgmres

  private
  ! ..................................................
  sll_int32, parameter  :: sll_solver_bool_false              = 0               !< code id for False
  sll_int32, parameter  :: sll_solver_maxiter                 = 1000            !< default maximum number of iterations for
  sll_real64, parameter :: sll_solver_tolerance               = 1.0d-9          !< default tolerance for iterative solvers
  sll_int32, parameter  :: sll_solver_restart                 = 30              !< default number of restarts for gmres
  ! ..................................................
  !> @brief 
  !> class for a sequential gmres linear solver
  type, extends(sll_t_linear_solver_iter_abstract) :: sll_t_linear_solver_mgmres
     sll_int32 :: n_mr           = 2     !< number of restarts for gmres
     sll_real64 :: rtol = 1.0d-14 !< relative tolerance

   contains
     procedure :: create            => create_linear_solver_mgmres 
     procedure :: initialize        => initialize_linear_solver_mgmres
     procedure :: set_guess         => set_guess_linear_solver_mgmres
     procedure :: check_convergence => check_convergence_linear_solver_mgmres

     procedure :: read_from_file    => read_from_file_linear_solver_mgmres
     procedure :: set_verbose       => set_verbose_linear_solver_mgmres
     procedure :: solve_real        => solve_real_linear_solver_mgmres 
     procedure :: print_info        => print_info_linear_solver_mgmres 
     procedure :: free              => free_linear_solver_mgmres 

  end type sll_t_linear_solver_mgmres
  ! ..................................................

contains

  ! ..................................................
  !> @brief     creates a linear solver 
  !>
  !> @param[inout] self             the current object 
  !> @param[in]    linear_operator  a linear operator 
  !> @param[in]    pc_left          a left preconditioner. This should be a solver too. [optional] 
  !> @param[in]    pc_right          a right preconditioner. This should be a solver too. [optional] 
  !> @param[in]    filename         parameter filename [optional]
  subroutine create_linear_solver_mgmres(self, linear_operator, pc_left, filename)
    implicit none
    class(sll_t_linear_solver_mgmres)                      , intent(inout) :: self 
    class(sll_t_linear_operator_abstract)          , target, intent(in)    :: linear_operator 
    class(sll_t_linear_solver_abstract)  , optional, target, intent(in)    :: pc_left 
    character(len=*), optional,                              intent(in)    :: filename

    ! ...
    if (present(pc_left)) then
       self%ptr_pc_left => pc_left 
    end if
    ! ...

    ! ...
    call self%initialize(linear_operator)  
    ! ...

    ! ...
    allocate(self%x_0(self%n_total_rows))
    self%x_0 = 0.0_f64
    ! ...

    ! ...
    if (present(filename)) then
       call self%read_from_file(filename)
    end if
    ! ...

  end subroutine create_linear_solver_mgmres
  ! ..................................................

  ! ..................................................
  !> @brief     initializes the linear solver 
  !>
  !> @param[inout] self               the current object 
  !> @param[in]    linear_operator    a linear operator
  !> @param[in]    x_0 [optional] the initial guess, defalut value: 0
  subroutine initialize_linear_solver_mgmres(self, linear_operator, x_0)
    implicit none
    class(sll_t_linear_solver_mgmres), intent(inout) :: self 
    class(sll_t_linear_operator_abstract), target, intent(in) :: linear_operator 
    sll_real64,dimension(:), optional, intent(in) :: x_0

    ! ...
    if (present(x_0)) then
       self%x_0 = x_0
    end if
    ! ...

    ! ... linear solver abstract initialization
    call self%initialize_abstract(linear_operator)  
    ! ...

    ! ... iterative linear solver abstract initialization
    call self%set_linear_operator(linear_operator)  
    ! ...

  end subroutine initialize_linear_solver_mgmres
  ! ..................................................

  ! ..................................................
  !> @brief     sets the initial guess
  !>
  !> @param[inout] self     the current object 
  !> @param[in]    x_0  the initial guess
  subroutine set_guess_linear_solver_mgmres(self, x_0)
    implicit none
    class(sll_t_linear_solver_mgmres), intent(inout) :: self 
    sll_real64,dimension(:), intent(in) :: x_0

    self%x_0 = x_0

  end subroutine set_guess_linear_solver_mgmres
  ! ..................................................

  ! ..................................................
  !> @brief     check the convergence of the current linear solver 
  !>
  !> @param[inout] self                  the current object 
  !> @param[in]    i_iteration           the checking iteration 
  !> @param[in]    flag                  the verbose flag
  !> @param[in]    r_err      [optional] the error
  !> @param[in]    arr_err    [optional] the array of error
  ! ............................................
  subroutine check_convergence_linear_solver_mgmres(self, i_iteration, flag, r_err, arr_err)
    implicit none
    class(sll_t_linear_solver_mgmres),      intent(in)     :: self 
    sll_int32,                                  intent(in)     :: i_iteration 
    logical,                                  intent(inout)  :: flag 
    sll_real64,                   optional, intent(in)     :: r_err 
    sll_real64, dimension(:),     optional, intent(in)     :: arr_err 

    ! ... 
    if (self%verbose) then
       if (present(r_err)) then
          if (i_iteration <= self%n_maxiter) then
             print*, '* mgmres:  convegence after', i_iteration, 'iterations. Error ', r_err
          else
             print *, '* mgmres: Warning - max iterations', self%n_maxiter, 'achieved without convergence. Error', r_err
          end if
       end if
    end if
    !...

  end subroutine check_convergence_linear_solver_mgmres
  ! ............................................

  ! ..................................................
  !> @brief     read from file 
  !>
  !> @param[inout] self                      the current object 
  !> @param[in]    filename       [optional] name of the output file
  subroutine read_from_file_linear_solver_mgmres(self, filename)
    implicit none
    class(sll_t_linear_solver_mgmres), intent(inout) :: self
    character(len=*)                , intent(in) :: filename
    ! local
    sll_int32, parameter :: input_file_id = 111
    sll_int32           :: int_maxiter    = sll_solver_maxiter                 
    sll_int32           :: int_restart    = sll_solver_restart                 
    sll_int32           :: int_null_space = sll_solver_bool_false 
    sll_int32           :: int_verbose    = sll_solver_bool_false 
    sll_real64 :: real_rtol      = sll_solver_tolerance
    sll_real64 :: real_atol      = sll_solver_tolerance

    ! ... 
    namelist /linear_solver/  int_null_space, &
         &  int_maxiter, &
         &  int_restart, &
         &  real_rtol, &
         &  real_atol, &
         &  int_verbose
    ! ... 

    ! ... 
    open(unit=input_file_id, file=trim(filename))
    read( input_file_id, linear_solver)
    ! ... 

    ! ...  
    self%n_maxiter        = int_maxiter
    self%n_mr             = int_restart
    self%rtol             = real_rtol
    self%atol             = real_atol
    ! ... 

    ! ... 
    self%null_space       = .false. 
    if (int_null_space == 1) then
       self%null_space     = .true. 
    end if
    ! ... 

    ! ... 
    self%verbose       = .false. 
    if (int_verbose == 1) then
       self%verbose     = .true. 
    end if
    ! ... 

    ! ... 
    close(input_file_id)
    ! ... 

  end subroutine read_from_file_linear_solver_mgmres
  ! ..................................................

  ! ............................................
  !> @brief     sets the verbose for the linear solver object
  !>
  !> @param[inout] self     the current object 
  !> @param[in]    verbose  verbose flag
  subroutine set_verbose_linear_solver_mgmres(self, verbose)
    implicit none
    class(sll_t_linear_solver_mgmres), intent(inout) :: self
    logical                   , intent(in) :: verbose

    ! ...
    call self%set_verbose_abstract(verbose)
    ! ...

  end subroutine set_verbose_linear_solver_mgmres
  ! ..................................................

  ! ..................................................
  !> @brief  solves the linear system with real vectors 
  !>
  !> @param[inout] self    the current object 
  !> @param[in]    rhs     the right hand side 
  !> @param[inout] unknown the solution 
  subroutine solve_real_linear_solver_mgmres(self, rhs, unknown)
    implicit none
    class(sll_t_linear_solver_mgmres), intent(inout) :: self 
    sll_real64, dimension(:), intent(in   )  :: rhs 
    sll_real64, dimension(:), intent(  out)  :: unknown
    ! local
    sll_real64, dimension(:), allocatable  :: l_rhs
    sll_int32           :: itr_used
    sll_real64 :: res
    logical :: flag

    ! ...
    allocate(l_rhs(self%n_total_rows)) 
    l_rhs = rhs
    ! ...

    ! ...
    if (associated(self%ptr_pc_left)) then
       call self%ptr_pc_left % solve(rhs, l_rhs) 
    end if
    ! ...

    ! ...
    unknown = self%x_0
    ! ...

    ! ...
    call mgmres_linear_solver ( self, unknown, l_rhs, &
         & self%n_maxiter, self%n_mr, & 
         & self%atol, self%rtol, &
         & itr_used, res)
    ! ...

    ! ...
    deallocate(l_rhs)
    ! ...

    ! ...
    call self%check_convergence( i_iteration=itr_used, &
         & flag=flag, &
         & r_err=res)
    ! ...

  end subroutine solve_real_linear_solver_mgmres
  ! ..................................................


  ! ...................................................
  subroutine mgmres_linear_solver( self, x, rhs, itr_max, mr,  &
       &  tol_abs, tol_rel, itr_used, rho )

    !*****************************************************************************80
    !
    !! mgmres_csr applies restarted gmres to a sparse triplet matrix.
    !
    !  discussion:
    !
    !    the linear system a*x=b is solved iteratively.
    !
    !    the matrix a is assumed to be stored in sparse triplet form.  only
    !    the nonzero entries of a are stored.  for instance, the k-th nonzero
    !    entry in the matrix is stored by:
    !
    !      a(k) = value of entry,
    !      ia(k) = row of entry,
    !      ja(k) = column of entry.
    !
    !    thanks to jesus pueblas sanchez-guerra for supplying two
    !    corrections to the code on 31 may 2007.
    !
    !  licensing:
    !
    !    this code is distributed under the gnu lgpl license.
    !
    !  modified:
    !
    !    13 july 2007
    !
    !  author:
    !
    !    original c version by lili ju.
    !    fortran90 version by john burkardt.
    !
    !  reference:
    !
    !    richard barrett, michael berry, tony chan, james demmel,
    !    june donato, jack dongarra, victor eijkhout, roidan pozo,
    !    charles romine, henk van der vorst,
    !    templates for the solution of linear systems:
    !    building blocks for iterative methods,
    !    siam, 1994.
    !    isbn: 0898714710,
    !    lc: qa297.8.t45.
    !
    !    tim kelley,
    !    iterative methods for linear and nonlinear equations,
    !    siam, 2004,
    !    isbn: 0898713528,
    !    lc: qa297.8.k45.
    !
    !    yousef saad,
    !    iterative methods for sparse linear systems,
    !    second edition,
    !    siam, 2003,
    !    isbn: 0898715342,
    !    lc: qa188.s17.
    !
    !  parameters:
    !
    !    input, sll_int32 ( kind = 4 ) n, the order of the linear system.
    !
    !    input, sll_int32 ( kind = 4 ) nz_num, the number of nonzero matrix values.
    !
    !    input, sll_int32 ( kind = 4 ) ia(nz_num), ja(nz_num), the row and column
    !    indices of the matrix values.
    !
    !    input, real ( kind = 8 ) a(nz_num), the matrix values.
    !
    !    input/output, real ( kind = 8 ) x(n); on input, an approximation to
    !    the solution.  on output, an improved approximation.
    !
    !    input, real ( kind = 8 ) rhs(n), the right hand side of the linear system.
    !
    !    input, sll_int32 ( kind = 4 ) itr_max, the maximum number of (outer)
    !    iterations to take.
    !
    !    input, sll_int32 ( kind = 4 ) mr, the maximum number of (inner) iterations
    !    to take.  0 < mr <= n.
    !
    !    input, real ( kind = 8 ) tol_abs, an absolute tolerance applied to the
    !    current residual.
    !
    !    input, real ( kind = 8 ) tol_rel, a relative tolerance comparing the
    !    current residual to the initial residual.
    !
    implicit none
    class(sll_t_linear_solver_mgmres), intent(in) :: self 

    sll_int32 :: mr

    sll_real64 :: av
    sll_real64 :: c(1:mr)
    sll_real64, parameter :: delta = 1.0d-03
    sll_real64 :: g(1:mr+1)
    sll_real64 :: h(1:mr+1,1:mr)
    sll_real64 :: htmp
    sll_int32 :: i
    sll_int32 :: itr
    sll_int32 :: itr_max
    sll_int32 :: itr_used
    sll_int32 :: j
    sll_int32 :: k
    sll_int32 :: k_copy
    sll_real64 :: mu
    sll_real64 :: r(1:self%n_total_rows)
    sll_real64 :: p(1:self%n_total_rows)
    sll_real64 :: rho
    sll_real64 :: rho_tol
    sll_real64 :: rhs(1:self%n_total_rows)
    sll_real64 :: s(1:mr)
    sll_real64 :: tol_abs
    sll_real64 :: tol_rel
    sll_real64 :: v(1:self%n_total_rows,1:mr+1)
    ! logical, parameter :: verbose = .false.
    sll_real64 :: x(1:self%n_total_rows)
    sll_real64 :: y(1:mr+1)

    itr_used = 0

    if ( self%verbose .and. ( self%n_total_rows < mr ) ) then
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'mgmres_csr - fatal error!'
       write ( *, '(a)' ) '  n < mr.'
       write ( *, '(a,i8)' ) '  n = ', self%n_total_rows
       write ( *, '(a,i8)' ) '  mr = ', mr
       stop
    end if

    r =0.0_f64

    do itr = 1, itr_max

       call self%ptr_linear_operator % dot(x, r)

       ! ...
       if (associated(self%ptr_pc_left)) then
          call self%ptr_pc_left % solve(r, p)
          r=p
       end if
       ! ...

       r(1:self%n_total_rows) = rhs(1:self%n_total_rows) - r(1:self%n_total_rows)

       rho = sqrt ( dot_product ( r(1:self%n_total_rows), r(1:self%n_total_rows) ) )

       !      if ( self%verbose ) then
       !        write ( *, '(a,i8,a,g14.6)' ) '  itr = ', itr, '  residual = ', rho
       !      end if

       if ( itr == 1 ) then
          rho_tol = rho * tol_rel
       end if

       if ( rho <= rho_tol .or. rho <= tol_abs ) then
          exit
       end if
       v(1:self%n_total_rows,1) = r(1:self%n_total_rows) / rho


       g(1) = rho
       g(2:mr+1) = 0.0d+00

       h(1:mr+1,1:mr) = 0.0d+00

       do k = 1, mr

          k_copy = k

          call self%ptr_linear_operator % dot(v(1:self%n_total_rows,k), v(1:self%n_total_rows,k+1))  

          ! ...
          if (associated(self%ptr_pc_left)) then
             call self%ptr_pc_left % solve(v(1:self%n_total_rows,k+1), p)
             v(1:self%n_total_rows,k+1) = p
          end if
          ! ...

          av = sqrt ( dot_product ( v(1:self%n_total_rows,k+1), v(1:self%n_total_rows,k+1) ) )


          do j = 1, k
             h(j,k) = dot_product ( v(1:self%n_total_rows,k+1), v(1:self%n_total_rows,j) )
             v(1:self%n_total_rows,k+1) = v(1:self%n_total_rows,k+1) - h(j,k) * v(1:self%n_total_rows,j)
          end do

          h(k+1,k) = sqrt ( dot_product ( v(1:self%n_total_rows,k+1), v(1:self%n_total_rows,k+1) ) )

          if ( av + delta * h(k+1,k) == av ) then

             do j = 1, k
                htmp = dot_product ( v(1:self%n_total_rows,k+1), v(1:self%n_total_rows,j) )
                h(j,k) = h(j,k) + htmp
                v(1:self%n_total_rows,k+1) = v(1:self%n_total_rows,k+1) - htmp * v(1:self%n_total_rows,j)
             end do

             h(k+1,k) = sqrt ( dot_product ( v(1:self%n_total_rows,k+1), v(1:self%n_total_rows,k+1) ) )

          end if

          if ( h(k+1,k) /= 0.0d+00 ) then
             v(1:self%n_total_rows,k+1) = v(1:self%n_total_rows,k+1) / h(k+1,k)
          end if

          if ( 1 < k ) then

             y(1:k+1) = h(1:k+1,k)

             do j = 1, k - 1
                call mult_givens ( c(j), s(j), j, y(1:k+1) )
             end do

             h(1:k+1,k) = y(1:k+1)

          end if

          mu = sqrt ( h(k,k)**2 + h(k+1,k)**2 )
          c(k) = h(k,k) / mu
          s(k) = -h(k+1,k) / mu
          h(k,k) = c(k) * h(k,k) - s(k) * h(k+1,k)
          h(k+1,k) = 0.0d+00
          call mult_givens ( c(k), s(k), k, g(1:k+1) )
          rho = abs ( g(k+1) )

          itr_used = itr_used + 1

          !        if ( self%verbose ) then
          !          write ( *, '(a,i8,a,g14.6)' ) '  k =   ', k, '  residual = ', rho
          !        end if

          if ( rho <= rho_tol .or. rho <= tol_abs ) then
             exit
          end if

       end do

       k = k_copy - 1

       y(k+1) = g(k+1) / h(k+1,k+1)

       do i = k, 1, -1
          y(i) = ( g(i) - dot_product ( h(i,i+1:k+1), y(i+1:k+1) ) ) / h(i,i)
       end do

       do i = 1, self%n_total_rows 
          x(i) = x(i) + dot_product ( v(i,1:k+1), y(1:k+1) )
       end do

       if ( rho <= rho_tol .or. rho <= tol_abs ) then
          exit
       end if

    end do

    if ( self%verbose ) then
       write ( *, '(a)'       ) ' '
       write ( *, '(a)'       ) 'mgmres_csr:'
       write ( *, '(a,i8)'    ) '  iterations = ', itr_used
       write ( *, '(a,g14.6)' ) '  final residual = ', rho
    end if
    return
  end subroutine mgmres_linear_solver
  ! ...................................................

  ! ...................................................
  subroutine mult_givens ( c, s, k, g )

    !*****************************************************************************80
    !
    !! mult_givens applies a givens rotation to two successive entries of a vector.
    !
    !  discussion:
    !
    !    in order to make it easier to compare this code with the original c,
    !    the vector indexing is 0-based.
    !
    !  licensing:
    !
    !    this code is distributed under the gnu lgpl license.
    !
    !  modified:
    !
    !    08 august 2006
    !
    !  author:
    !
    !    original c version by lili ju.
    !    fortran90 version by john burkardt.
    !
    !  reference:
    !
    !    richard barrett, michael berry, tony chan, james demmel,
    !    june donato, jack dongarra, victor eijkhout, roidan pozo,
    !    charles romine, henk van der vorst,
    !    templates for the solution of linear systems:
    !    building blocks for iterative methods,
    !    siam, 1994.
    !    isbn: 0898714710,
    !    lc: qa297.8.t45.
    !
    !    tim kelley,
    !    iterative methods for linear and nonlinear equations,
    !    siam, 2004,
    !    isbn: 0898713528,
    !    lc: qa297.8.k45.
    !
    !    yousef saad,
    !    iterative methods for sparse linear systems,
    !    second edition,
    !    siam, 2003,
    !    isbn: 0898715342,
    !    lc: qa188.s17.
    !
    !  parameters:
    !
    !    input, real ( kind = 8 ) c, s, the cosine and sine of a givens
    !    rotation.
    !
    !    input, sll_int32 ( kind = 4 ) k, indicates the location of the first
    !    vector entry.
    !
    !    input/output, real ( kind = 8 ) g(1:k+1), the vector to be modified.
    !    on output, the givens rotation has been applied to entries g(k) and g(k+1).
    !
    implicit none

    sll_int32 :: k

    sll_real64 :: c
    sll_real64 :: g(1:k+1)
    sll_real64 :: g1
    sll_real64 :: g2
    sll_real64 :: s

    g1 = c * g(k) - s * g(k+1)
    g2 = s * g(k) + c * g(k+1)

    g(k)   = g1
    g(k+1) = g2

    return
  end subroutine mult_givens
  ! ...................................................

  ! ............................................
  !> @brief      destroys a finite element cell 
  !>
  !> @param[inout] self the current object 
  subroutine print_info_linear_solver_mgmres(self)
    implicit none
    class(sll_t_linear_solver_mgmres), intent(in) :: self 
    ! local

    print *, ">>>> linear_solver_mgmres"
    print *, "* verbose    : ", self%verbose 
    print *, "* n_mr       : ", self%n_mr 
    print *, "* rtol       : ", self%rtol 
    print *, "* atol       : ", self%atol 
    print *, "* null_space : ", self%null_space 
    print *, "* n_maxiter  : ", self%n_maxiter 
    print *, "<<<< "

  end subroutine print_info_linear_solver_mgmres
  ! ............................................

  ! ..................................................
  !> @brief     destroys the current object 
  !>
  !> @param[inout] self   the current object 
  subroutine free_linear_solver_mgmres(self)
    implicit none
    class(sll_t_linear_solver_mgmres), intent(inout) :: self 

    deallocate (self%x_0)

  end subroutine free_linear_solver_mgmres
  ! ..................................................

end module sll_m_linear_solver_mgmres
