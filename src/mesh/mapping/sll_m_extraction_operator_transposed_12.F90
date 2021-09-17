module sll_m_extraction_operator_transposed_12
#include "sll_working_precision.h"
#include "sll_assert.h"
  
  use sll_m_linear_operator_abstract


  implicit none

  public :: sll_t_extraction_operator_transposed_12

  private
  type, extends(sll_t_linear_operator_abstract) :: sll_t_extraction_operator_transposed_12
     sll_int32 :: n_dofs(2)
     sll_int32 :: n0
     sll_int32 :: nk
     sll_real64, allocatable :: conp(:,:)
     sll_real64 :: tau
     sll_real64 :: tl(3, 3)
     

   contains
     procedure :: create => create_extraction_operator_transposed_12
     procedure :: free => free_extraction_operator_transposed_12
     procedure :: dot => dot_extraction_operator_transposed_12
     procedure :: print_info => print_info_extraction_operator_transposed_12

  end type sll_t_extraction_operator_transposed_12


contains

  subroutine create_extraction_operator_transposed_12( self, n_dofs, conp, k, tau )
    class(sll_t_extraction_operator_transposed_12), intent( inout ) :: self
    sll_int32, intent( in )  :: n_dofs(2) !< number of degrees of freedom
    sll_real64, intent( in ) :: conp(:,:) !< dofs of discrete spline mapping
    sll_int32, intent( in ) :: k !< C^k smooth polar splines
    sll_real64, intent( in ) :: tau !< diameter parameter for first control points of discrete spline mapping

    allocate( self%conp(product(n_dofs), 2) )
    self%n_dofs = n_dofs
    self%conp = conp
    self%tau = tau
    self%nk = 3!(k+1)*(k+2)/2
    self%n0 = 2*n_dofs(2)!(k+1)*n_dofs(2)

    self%tl(:,1) = 1._f64/3._f64
    self%tl(1,2) = 2._f64/(3._f64*self%tau)
    self%tl(2,2) = -1._f64/(3._f64*self%tau)
    self%tl(3,2) = -1._f64/(3._f64*self%tau)
    self%tl(1,3) = 0._f64
    self%tl(2,3) = 1._f64/(sqrt(3._f64)*self%tau)
    self%tl(3,3) = -1._f64/(sqrt(3._f64)*self%tau)
        
    self%n_rows = product(self%n_dofs)
    self%n_cols = product(self%n_dofs) - self%n0 + self%nk

    self%n_global_rows = self%n_rows
    self%n_global_cols = self%n_cols

  end subroutine create_extraction_operator_transposed_12

  subroutine free_extraction_operator_transposed_12( self )
    class(sll_t_extraction_operator_transposed_12), intent( inout ) :: self

  end subroutine free_extraction_operator_transposed_12


  subroutine dot_extraction_operator_transposed_12 ( self, x, y )
    class(sll_t_extraction_operator_transposed_12), intent( in ) :: self
    sll_real64, intent( in    ) :: x(:)
    sll_real64, intent(   out ) :: y(:)
    !local variables
    sll_int32 :: i, j, ind
    y=0.0_f64

    do i = 1, self%nk
       y(2) = y(2) + &
            ( (self%conp(1, 1)-self%conp(self%n_dofs(2), 1)) * self%tl(i,2) + &
            (self%conp(1, 2)-self%conp(self%n_dofs(2), 2)) * self%tl(i,3) )*&
            x(i)
    end do
    do j = 2, self%n_dofs(2)
       do i = 1, self%nk
          y(2+(j-1)*self%n_dofs(1)) = y(2+(j-1)*self%n_dofs(1)) + &
               ( (self%conp(j, 1)-self%conp(j-1, 1)) * self%tl(i,2) + &
               (self%conp(j, 2)-self%conp(j-1, 2)) * self%tl(i,3) )*&
               x(i)
       end do
    end do

    ind = self%nk+1
    do j = 1, self%n_dofs(2)
       do i = 3, self%n_dofs(1)
          y(i+(j-1)*self%n_dofs(1)) =  x(ind)
          ind = ind+1
       end do
    end do
    SLL_ASSERT(ind == self%n_cols+1)

  end subroutine dot_extraction_operator_transposed_12


  subroutine print_info_extraction_operator_transposed_12( self )
    class(sll_t_extraction_operator_transposed_12), intent(in) :: self

  end subroutine print_info_extraction_operator_transposed_12


end module sll_m_extraction_operator_transposed_12
