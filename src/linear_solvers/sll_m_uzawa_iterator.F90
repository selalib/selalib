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
     sll_int32 :: n_total !< product of number of degrees of freedom
 

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
  
    self%n_total = operator_b%n_global_cols

    self%n_rows = solver_a%n_rows
    self%n_cols = solver_a%n_cols

    allocate(self%x_0(self%n_total))
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
    sll_real64 :: rhs1(3*self%n_total)
    sll_real64 :: x2(self%n_total)
    sll_int32  :: itr_used
    sll_real64 :: res
    logical :: flag
    
    x2 = self%x_0
    call self%operator_b%dot(x2, rhs1)
    rhs1 = rhs - rhs1
    call self%solver_a%solve(rhs1, unknown)
    
    call uzawa_iterator(self, unknown, x2, itr_used, res)

    call self%check_convergence( i_iteration=itr_used, &
         & flag=flag, &
         & r_err=res)

  end subroutine solve_uzawa_iterator


  subroutine uzawa_iterator(self, x1, x2, niterx, res) 
    class(sll_t_uzawa_iterator), intent( in ) :: self !< Uzawa iterator
    sll_real64,                  intent(inout) :: x1(:)
    sll_real64,                  intent(inout) :: x2(:)
    sll_int32,                   intent(out) :: niterx
    sll_real64,                  intent(out) :: res
    !local variables
    sll_int32 :: k
    sll_real64 :: alpha, beta
    sll_real64 :: a1(3*self%n_total), p1(3*self%n_total)
    sll_real64 :: a2(self%n_total), p2(self%n_total), r2(self%n_total)

    call self%operator_bt%dot(x1, r2)
    p2=r2
    
    do k = 1, self%n_maxiter
       niterx    = k
       
       call self%operator_b%dot(p2, a1)
       call self%solver_a%solve(a1, p1)
       call self%operator_bt%dot(p1, a2)
       
       alpha = sum(p2*a2)/sum(p2*r2)
       
       x2 = x2 + alpha * p2
       r2 = r2 - alpha * a2
       x1 = x1 - alpha * p1
       
       res = sqrt(sum(r2*r2)/self%n_total)
       if( res <= self%atol ) exit

       beta = sum(r2*a2)/sum(p2*a2)
       p2 = r2 - beta*p2

    end do
    
    return
    
  end subroutine uzawa_iterator
  

  
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
             print*, '* uzawa_iterator:  convergence after', i_iteration, ' iterations. Error ', r_err
          else
             print *, '* uzawa_iterator: Warning - max iterations achieved without convergence. Error', r_err
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
