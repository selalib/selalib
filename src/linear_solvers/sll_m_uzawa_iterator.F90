module sll_m_uzawa_iterator
#include "sll_working_precision.h"

  use sll_m_linear_solver_abstract, only: &
       sll_t_linear_solver_abstract

  use sll_m_linear_solver_iter_abstract, only: &
       sll_t_linear_solver_iter_abstract

  use sll_m_linear_operator_abstract,  only: &
       sll_t_linear_operator_abstract
       
  implicit none

  public :: sll_t_uzawa_iterator
  
  private
  type, extends(sll_t_linear_solver_iter_abstract) :: sll_t_uzawa_iterator
     class(sll_t_linear_solver_abstract), pointer :: solver_a
     class(sll_t_linear_operator_abstract), pointer :: operator_b
     class(sll_t_linear_operator_abstract), pointer :: operator_bt
     sll_int32 :: n_total0 !< product of number of degrees of freedom
     sll_int32 :: n_total1 !< product of number of degrees of freedom
 

   contains
     procedure :: create => create_uzawa_iterator
     procedure :: set_guess         => set_guess_uzawa_iterator
     procedure :: check_convergence => check_convergence_uzawa_iterator
     procedure :: read_from_file => read_from_file_uzawa_iterator
     procedure :: set_verbose => set_verbose_uzawa_iterator
     procedure :: solve_real => solve_uzawa_iterator
     procedure :: print_info => print_info_uzawa_iterator
     procedure :: free => free_uzawa_iterator

  end type sll_t_uzawa_iterator


contains

  subroutine create_uzawa_iterator( self, solver_a, operator_b, operator_bt )
    class(sll_t_uzawa_iterator), intent( inout ) :: self !< Uzawa iterator
    class(sll_t_linear_solver_abstract), target :: solver_a
    class(sll_t_linear_operator_abstract), target :: operator_b
    class(sll_t_linear_operator_abstract), target :: operator_bt
    
    self%solver_a => solver_a
    self%operator_b => operator_b
    self%operator_bt => operator_bt
  
    self%n_total0 = operator_b%n_global_cols
    self%n_total1 = operator_b%n_global_rows

    self%n_rows = solver_a%n_rows
    self%n_cols = solver_a%n_cols

    allocate(self%x_0(self%n_total0))
    self%x_0 = 0.0_f64
    
  end subroutine create_uzawa_iterator

  subroutine free_uzawa_iterator( self )
    class(sll_t_uzawa_iterator), intent( inout ) :: self !< Uzawa iterator

    deallocate (self%x_0)
    self%solver_a => null()
    self%operator_b => null()
    self%operator_bt => null()
    
  end subroutine free_uzawa_iterator
  
  
  subroutine solve_uzawa_iterator(self, rhs, unknown)
    class(sll_t_uzawa_iterator), intent( inout ) :: self !< Uzawa iterator
    sll_real64, intent( in    ) :: rhs(:) !< Inputvariable
    sll_real64, intent(   out ) :: unknown(:) !< Outputvariable
    !local variables
    sll_real64 :: rhs1(self%n_total1)
    sll_real64 :: x0(self%n_total0)
    sll_int32  :: itr_used
    sll_real64 :: res
    logical :: flag
    
    x0 = self%x_0
    call self%operator_b%dot(x0, rhs1)
    rhs1 = rhs - rhs1
    call self%solver_a%solve(rhs1, unknown)
    
    call uzawa_iterator(self, unknown, x0, itr_used, res)

    self%x_0 = x0

    call self%check_convergence( i_iteration=itr_used, &
         & flag=flag, &
         & r_err=res)

  end subroutine solve_uzawa_iterator


  subroutine uzawa_iterator(self, x1, x0, niterx, res) 
    class(sll_t_uzawa_iterator), intent( in ) :: self !< Uzawa iterator
    sll_real64,                  intent(inout) :: x1(:)
    sll_real64,                  intent(inout) :: x0(:)
    sll_int32,                   intent(out) :: niterx
    sll_real64,                  intent(out) :: res
    !local variables
    sll_int32 :: k
    sll_real64 :: alpha, beta
    sll_real64 :: a0(self%n_total0), p0(self%n_total0), r0(self%n_total0)
    sll_real64 :: a1(self%n_total1), p1(self%n_total1)
  

    call self%operator_bt%dot(x1, r0)
    p0=r0
    niterx = 1
    do k = 1, self%n_maxiter
       
       
       call self%operator_b%dot(p0, a1)
       call self%solver_a%solve(a1, p1)
       call self%operator_bt%dot(p1, a0)
       
       alpha = sum(p0*a0)/sum(p0*r0)
       
       x0 = x0 + alpha * p0
       r0 = r0 - alpha * a0
       x1 = x1 - alpha * p1
       
       res = sqrt(sum(r0*r0)/real(self%n_total0,f64))
       !print*, 'residuum', res
       if( res <= self%atol ) exit

       if(sqrt(sum(p0*p0)/real(self%n_total0,f64))<= self%atol) then
          print*, 'error uzawa iterator: krylov subspace exhausted'
          exit
          !p0 = r0
       else
          beta = sum(r0*a0)/sum(p0*a0)
          p0 = r0 - beta*p0
       end if
       
       niterx = k + 1
    end do
    
    return
    
  end subroutine uzawa_iterator

  subroutine uzawa_iterator2(self, x1, x0, niterx, res) 
    class(sll_t_uzawa_iterator), intent( in ) :: self !< Uzawa iterator
    sll_real64,                  intent(inout) :: x1(:)
    sll_real64,                  intent(inout) :: x0(:)
    sll_int32,                   intent(out) :: niterx
    sll_real64,                  intent(out) :: res
    !local variables
    sll_int32 :: k
    sll_real64 :: alpha, beta
    sll_real64 :: a0(self%n_total0), q0(self%n_total0), r0(self%n_total0), d0(self%n_total0)
    sll_real64 :: a1(self%n_total1), q1(self%n_total1)
  

    call self%operator_bt%dot(x1, q0)
    niterx = 1
    do k = 1, self%n_maxiter
       call self%operator_bt%dot(x1, d0)
       
       call self%operator_b%dot(q0, a1)
       call self%solver_a%solve(a1, q1)
       call self%operator_bt%dot(q1, a0)
       
       alpha = sum(d0*d0)/sum(d0*a0)
       
       x0 = x0 + alpha * q0
       x1 = x1 + alpha * q1

       call self%operator_bt%dot(x1, r0)
       res = sqrt(sum(r0*r0)/real(self%n_total0,f64))
       print*, 'residuum', res
       if( res <= self%atol ) exit

       if(sqrt(sum(q0*q0)/real(self%n_total0,f64))<= self%atol) then
          print*, 'error uzawa iterator: krylov subspace exhausted'
          exit
       else
          beta = sum(r0*r0)/sum(d0*d0)
          q0 = r0 + beta*q0
       end if
       
       niterx = k + 1
    end do
    
    return
    
  end subroutine uzawa_iterator2

  
  subroutine print_info_uzawa_iterator( self )
    class(sll_t_uzawa_iterator), intent(in) :: self !< Uzawa iterator

    print *, ">>>> uzawa_iterator"
    print *, "* verbose    : ", self%verbose 
    print *, "* atol       : ", self%atol 
    print *, "* n_maxiter  : ", self%n_maxiter 
    print *, "<<<< "
    
  end subroutine print_info_uzawa_iterator

   subroutine set_guess_uzawa_iterator(self, x_0)
     class(sll_t_uzawa_iterator), intent(inout) :: self !< Uzawa iterator
    sll_real64,dimension(:), intent(in) :: x_0

    self%x_0 = x_0

  end subroutine set_guess_uzawa_iterator

  ! ............................................
  subroutine check_convergence_uzawa_iterator(self, i_iteration, flag, r_err, arr_err)
    class(sll_t_uzawa_iterator),      intent(in)     :: self !< Uzawa iterator
    ! local
    sll_int32,                                  intent(in)     :: i_iteration 
    logical,                                  intent(inout)  :: flag 
    sll_real64,                   optional, intent(in)     :: r_err 
    sll_real64, dimension(:),     optional, intent(in)     :: arr_err 

    ! ... 
    if (self%verbose) then
       if (present(r_err)) then
          if (i_iteration <= self%n_maxiter) then
             print*, '* uzawa_iterator:  convergence after', i_iteration, 'iterations. Error ', r_err
          else
             print *, '* uzawa_iterator: Warning - max iterations', self%n_maxiter, 'achieved without convergence. Error', r_err
          end if
       end if
    end if
    !...

  end subroutine check_convergence_uzawa_iterator

  subroutine read_from_file_uzawa_iterator( self, filename )
    class(sll_t_uzawa_iterator), intent(inout) :: self !< Uzawa iterator
    character(len=*)                   , intent(in)    :: filename

  end subroutine read_from_file_uzawa_iterator

  subroutine set_verbose_uzawa_iterator( self, verbose )
    class(sll_t_uzawa_iterator), intent(inout) :: self !< Uzawa iterator
    logical                            , intent(in)    :: verbose 

    call self%set_verbose_abstract(verbose)
    
  end subroutine set_verbose_uzawa_iterator


end module sll_m_uzawa_iterator
