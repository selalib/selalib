!> @brief 
!> module for conjugate gradient method in pure form  
!> @details
!>  a linear solver using the conjugate gradient method 
!>  This is a modification of linear_solver_cg 
!>
!> Maintainer   ARA, Katharina Kormann
!> Modified by Benedikt Perse	
!> Stability	stable

module sll_m_linear_solver_cg
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
       sll_t_linear_solver_cg

  private
  ! ..................................................
  sll_int32, parameter  :: sll_solver_bool_false              = 0               !< code id for False
  sll_int32, parameter  :: sll_solver_maxiter                 = 1000            !< default maximum number of iterations for
  sll_real64, parameter :: sll_solver_tolerance               = 1.0d-9          !< default tolerance for iterative solvers

  ! ..................................................
  !> @brief 
  !> class for the cg linear solver
  type, extends(sll_t_linear_solver_iter_abstract) :: sll_t_linear_solver_cg
   contains
     procedure :: create            => create_linear_solver_cg
     procedure :: initialize        => initialize_linear_solver_cg
     procedure :: set_guess         => set_guess_linear_solver_cg
     procedure :: check_convergence => check_convergence_linear_solver_cg
     procedure :: read_from_file    => read_from_file_linear_solver_cg 
     procedure :: solve_real        => solve_real_linear_solver_cg 
     procedure :: set_verbose       => set_verbose_linear_solver_cg
     procedure :: print_info        => print_info_linear_solver_cg
     procedure :: free              => free_linear_solver_cg
  end type sll_t_linear_solver_cg
  ! ..................................................

contains

  ! ..................................................
  !> @brief     creates a linear solver 
  !>
  !> @param[inout] self             the current object 
  !> @param[in]    linear_operator  a linear operator
  !> @param[in]    filename         parameter filename [optional]
  subroutine create_linear_solver_cg(self, linear_operator, pc_left, filename)
    implicit none
    class(sll_t_linear_solver_cg),                         intent(inout) :: self 
    class(sll_t_linear_operator_abstract),         target, intent(in)    :: linear_operator
    class(sll_t_linear_solver_abstract), optional, target, intent(in)    :: pc_left 
    character(len=*), optional,                            intent(in)    :: filename

    ! ...
    if( present(pc_left) ) then
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

  end subroutine create_linear_solver_cg
  ! ..................................................

  ! ..................................................
  !> @brief     initializes the linear solver 
  !>
  !> @param[inout] self               the current object 
  !> @param[in]    linear_operator    a linear operator
  !> @param[in]    x_0               [optional] the init guess, defalut value: 0
  subroutine initialize_linear_solver_cg(self, linear_operator, x_0)
    implicit none
    class(sll_t_linear_solver_cg), intent(inout) :: self 
    class(sll_t_linear_operator_abstract), target, intent(in) :: linear_operator 
    sll_real64,dimension(:), optional, intent(in) :: x_0
    ! local

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

  end subroutine initialize_linear_solver_cg
  ! ..................................................

  ! ..................................................
  !> @brief     sets the initial guess 
  !>
  !> @param[inout] self     the current object 
  !> @param[in]    x_guess  initial guess
  subroutine set_guess_linear_solver_cg(self, x_0)
    implicit none
    class(sll_t_linear_solver_cg), intent(inout) :: self 
    sll_real64,dimension(:), intent(in) :: x_0

    self%x_0 = x_0

  end subroutine set_guess_linear_solver_cg
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
  subroutine check_convergence_linear_solver_cg(self, i_iteration, flag, r_err, arr_err)
    implicit none
    class(sll_t_linear_solver_cg),      intent(in)     :: self 
    ! local
    sll_int32,                                  intent(in)     :: i_iteration 
    logical,                                  intent(inout)  :: flag 
    sll_real64,                   optional, intent(in)     :: r_err 
    sll_real64, dimension(:),     optional, intent(in)     :: arr_err 

    ! ... 
    if (self%verbose) then
       if (present(r_err)) then
          if (i_iteration <= self%n_maxiter) then
             print*, '* cg:  convergence after', i_iteration, 'iterations. Error ', r_err
          else
             print *, '* cg: Warning - max iterations', self%n_maxiter, 'achieved without convergence. Error', r_err
          end if
       end if
    end if
    !...

  end subroutine check_convergence_linear_solver_cg
  ! ............................................

  ! ..................................................
  !> @brief     read from file 
  !>
  !> @param[inout] self                      the current object 
  !> @param[in]    filename       [optional] name of the output file
  subroutine read_from_file_linear_solver_cg(self, filename)
    implicit none
    class(sll_t_linear_solver_cg), intent(inout) :: self
    character(len=*)                , intent(in) :: filename
    ! local
    sll_int32, parameter :: input_file_id = 111
    sll_int32           :: int_maxiter        = sll_solver_maxiter                 
    sll_int32           :: int_null_space     = sll_solver_bool_false
    sll_int32           :: int_verbose        = sll_solver_bool_false 
    sll_real64 :: real_atol          = sll_solver_tolerance
    ! ... 
    namelist /linear_solver/  int_null_space, &
         &  int_maxiter, &
         &  real_atol, &
         &  int_verbose
    ! ... 

    ! ... 
    open(unit=input_file_id, file=trim(filename))
    read( input_file_id, linear_solver)
    ! ... 

    ! ...  
    self%n_maxiter        = int_maxiter
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

  end subroutine read_from_file_linear_solver_cg
  ! ..................................................

  ! ..................................................
  !> @b af  solves the linear system with real vectors 
  !>
  !> @param[inout] self    the current object 
  !> @param[in]    rhs     the right hand side 
  !> @param[inout] unknown the solution 
  subroutine solve_real_linear_solver_cg(self, rhs, unknown)
    implicit none
    class(sll_t_linear_solver_cg), intent(inout) :: self 
    sll_real64, dimension(:), intent(in   )  :: rhs 
    sll_real64, dimension(:), intent(  out)  :: unknown
    ! local
    sll_real64, dimension(:), allocatable  :: l_rhs
    sll_int32, parameter :: current_dof = 1
    sll_int32           :: itr_used
    sll_real64 :: res
    logical :: flag

    ! ...
    allocate(l_rhs(self%n_total_rows)) 
    l_rhs = rhs
    ! ...


    ! ...
    unknown = self%x_0
    ! ...

    ! ...
    call cg_linear_solver(self, unknown, l_rhs, itr_used, res)
    ! ...

    ! ...
    deallocate(l_rhs)
    ! ...

    ! ...
    call self%check_convergence( i_iteration=itr_used, &
         & flag=flag, &
         & r_err=res)
    ! ...

  end subroutine solve_real_linear_solver_cg
  ! ..................................................

  ! ..................................................
  !> @brief  cg interface 
  !>
  !> @param[inout] self    the current object 
  !> @param[in]    b       the right hand side 
  !> @param[out]   u       the solution 
  !> @param[out]   niterx  number of iteration for convergence 
  !> @param[out]   res     final residual error
  subroutine cg_linear_solver(self, x, b, niterx, res)
    implicit none

    class(sll_t_linear_solver_cg), intent(in)  :: self
    sll_real64, dimension (:),   intent(inout) :: x
    sll_real64, dimension (:),   intent(in)  :: b
    sll_int32,                   intent(out) :: niterx
    sll_real64,                  intent(out) :: res
    ! local
    sll_int32 :: niterxmax
    sll_int32 :: k, l
    sll_int32 :: nb
    sll_real64 :: epsx
    sll_real64 :: lambda
    sll_real64 :: alpha
    sll_real64 :: w1
    sll_real64, dimension (self%n_total_rows) :: p
    sll_real64, dimension (self%n_total_rows) :: r
    sll_real64, dimension (self%n_total_rows) :: v
    sll_real64, dimension (self%n_total_rows) :: z

    ! ...
    niterxmax = self%n_maxiter
    nb        = self%n_total_rows
    epsx      = self%atol
    ! ...

    v = 0.0_f64
    call self%ptr_linear_operator%dot (x, v)
    !r_0=b-A x_0
    r = b - v
    !p_0=M^{-1}r_0
    if( associated( self%ptr_pc_left) )then
       call self%ptr_pc_left%solve(r, p)
    else
       p = r
    end if
    !alpha_0= <r_0,p_0>
    alpha = sum(r*p)

    niterx = 1
    do k = 1, niterxmax

       v = 0.0_f64
       !v_j=A p_j
       call self%ptr_linear_operator%dot (p, v)

       !lambda_j=<r_j,z_j>/<A p_j,p_j>
       lambda = alpha / (sum(p*v)+1.d-25)

       w1 = 0.0_f64
       !x_j+1=x_j+lambda_j p_j
       !r_j+r=r_j-lambda_j Ap_j
       do l=1,nb
          x(l)      = x(l) + lambda*p(l)
          r(l)      = r(l) - lambda*v(l)
          w1 = w1 + r(l) * r(l)
       end do

       res=sqrt(w1/real(nb,f64))
       if (res <= epsx) exit
       if( associated( self%ptr_pc_left) )then
          !z_j+1=M^{-1}r_j+1
          call self%ptr_pc_left%solve(r, z)
          !beta_j= <r_j+1,z_j+1>/<r_j,z_j>
          w1 = sum(r*z)
       else
          z = r
       end if
       !p_j+1=z_j+1+beta_j p_j
       p = z + w1/alpha * p
       alpha = w1

       niterx = k + 1

    end do

    return

  end subroutine cg_linear_solver
  ! ...................................................

  ! ............................................
  !> @brief     sets the verbose for the linear solver object
  !>
  !> @param[inout] self     the current object 
  !> @param[in]    verbose  verbose flag
  subroutine set_verbose_linear_solver_cg(self, verbose)
    implicit none
    class(sll_t_linear_solver_cg), intent(inout) :: self
    logical                   , intent(in) :: verbose

    ! ...
    call self%set_verbose_abstract(verbose)
    ! ...

  end subroutine set_verbose_linear_solver_cg
  ! ..................................................

  ! ............................................
  !> @brief      destroys a finite element cell 
  !>
  !> @param[inout] self the current object 
  subroutine print_info_linear_solver_cg(self)
    implicit none
    class(sll_t_linear_solver_cg), intent(in) :: self 
    ! local

    print *, ">>>> linear_solver_cg"
    print *, "* verbose    : ", self%verbose 
    print *, "* atol       : ", self%atol 
    print *, "* null_space : ", self%null_space 
    print *, "* n_maxiter  : ", self%n_maxiter 
    print *, "<<<< "

  end subroutine print_info_linear_solver_cg
  ! ............................................

  ! ..................................................
  !> @brief     destroys the current object 
  !>
  !> @param[inout] self   the current object
  subroutine free_linear_solver_cg(self)
    implicit none
    class(sll_t_linear_solver_cg), intent(inout) :: self 

    deallocate (self%x_0)

  end subroutine free_linear_solver_cg
  ! ..................................................

end module sll_m_linear_solver_cg
