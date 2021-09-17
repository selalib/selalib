!> @brief 
!> module for a penalized linear operator 
!> @details
!> a penalized linear operator
!>
!>
!> Maintainer   ARA
!> Modified by Benedikt Perse	
!> Stability	stable

module sll_m_linear_operator_penalized
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use sll_m_linear_operator_abstract, only: &
       sll_t_linear_operator_abstract

  implicit none

  public :: &
       sll_t_linear_operator_penalized

  private
  ! ..................................................
  !> @brief 
  !> class for a linear operator_penalized 
  type, extends(sll_t_linear_operator_abstract) :: sll_t_linear_operator_penalized
     sll_int32 :: n_dim_nullspace = 1  !< dimension of the null space

     sll_real64, dimension(:,:), allocatable   :: vecs !< list of vectors that have a non vanishing compoenent on the kernel 

     class(sll_t_linear_operator_abstract), pointer :: ptr_linear_operator => null() !< pointer to the initial linear operator
   contains
     procedure :: create       => create_linear_operator_penalized
     procedure :: free         => free_linear_operator_penalized
     procedure :: dot => dot_linear_operator_penalized 
     procedure :: print_info   => print_info_linear_operator_penalized
  end type sll_t_linear_operator_penalized
  ! ..................................................

contains

  ! ............................................
  !> @brief     creates a linear operator_penalized 
  !>
  !> @param[inout] self                         the current object 
  !> @param[in]    linear_operator  [optional]  a linear_operator object 
  !> @param[in]    vecs             [optional]  vectors that spans the null space. dimension (n_dim_nullspace, n_cols).
  !> @param[in]    n_dim_nullspace  [optional]  dimension of the null space / kernel 
  subroutine create_linear_operator_penalized( self, &
       & linear_operator, &
       & vecs, n_dim_nullspace)
    implicit none
    class(sll_t_linear_operator_penalized)                 , intent(inout) :: self 
    class(sll_t_linear_operator_abstract), target, optional, intent(in)    :: linear_operator 
    sll_real64,dimension(:,:)                  , optional, intent(in)    :: vecs  
    sll_int32                                      , optional, intent(in)    :: n_dim_nullspace
    ! local

    if (present(linear_operator)) then
       self % ptr_linear_operator => linear_operator
    end if
    ! ...

    ! ...
    call self % initialize_abstract (other=linear_operator)
    ! ...

    ! ...
    self % n_dim_nullspace = 1  
    if (present(n_dim_nullspace)) then
       self % n_dim_nullspace = n_dim_nullspace
    end if
    ! ...

    ! ...
    allocate(self % vecs(self % n_dim_nullspace, self % n_global_cols)) 
    self % vecs = 0.0_f64
    ! ...

    ! ... \todo add assertions
    if (present(vecs)) then
       self % vecs = vecs 
    else
       if (self % n_dim_nullspace == 1) then
          self % vecs = 1.0_f64
       else
          stop "create_linear_operator_penalized: wrong arguments for the vecs attribut"
       end if
    end if
    ! ...

  end subroutine create_linear_operator_penalized
  ! ............................................

  ! ............................................
  !> @brief     destroys the current object 
  !>
  !> @param[inout] self   the current object 
  subroutine free_linear_operator_penalized(self)
    implicit none
    class(sll_t_linear_operator_penalized), intent(inout) :: self 
    ! local

    ! ...
    deallocate(self % vecs)
    ! ...

  end subroutine free_linear_operator_penalized
  ! ............................................

  ! ...................................................
  !> @brief     apply the dot operation 
  !>
  !> @param[inout] self   the current object 
  !> @param[in]    x      a real valued vector 
  !> @param[inout] y      a real valued vector 
  subroutine dot_linear_operator_penalized(self, x, y)
    implicit none
    class(sll_t_linear_operator_penalized), intent(in) :: self
    sll_real64,dimension(:), intent(in   ) :: x
    sll_real64,dimension(:), intent(  out) :: y  
    ! local 
    sll_int32 :: i

    ! ...
    if ( associated(self % ptr_linear_operator) .eqv. .true.) then
       call self % ptr_linear_operator % dot(x, y)
    else
       stop "dot_real_linear_operator_penalized: linear operator must be initialized"
    end if
    ! ...

    ! ... add penalization term
    do i = 1, self % n_dim_nullspace
       y(:) = y(:) + dot_product(self % vecs(i, :), x(:)) / real(self % n_global_cols, f64)
    end do
    ! ...

  end subroutine dot_linear_operator_penalized
  ! ...................................................

  ! ...................................................
  !> @brief      prints a linear operator 
  !>
  !> @param[inout] self the current object 
  subroutine print_info_linear_operator_penalized(self)
    implicit none
    class(sll_t_linear_operator_penalized), intent(in) :: self 
    ! local

    print *, ">>>> linear_operator_penalized"
    call self % print_info_abstract()

    print *, "* n_dim_nullspace        : ", self % n_dim_nullspace
    print *, "<<<< "

  end subroutine print_info_linear_operator_penalized
  ! ...................................................

end module sll_m_linear_operator_penalized
