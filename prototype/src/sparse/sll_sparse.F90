module sll_sparse
#include "sll_working_precision.h"
#include "sll_memory.h"

implicit none
 
private

!> Matrix in CSC format
type, public :: sll_csc_matrix
    sll_int32           :: n             !< Matrix dimension
    sll_int32           :: nnzeros       !< Number of non-zeros
    sll_int32, pointer  :: colptr(:)     !< CSC format colum pointer
    sll_int32, pointer  :: row(:)        !< CSC format row indices
    sll_real64, pointer :: avals(:)      !< non zeros values
end type sll_csc_matrix

!> Matrix in CSR format
type, public :: sll_csr_matrix
    sll_int32           :: n             !< Matrix dimension
    sll_int32           :: nnzeros       !< Number of non-zeros
    sll_int32, pointer  :: rowptr(:)     !< CSR format colum pointer
    sll_int32, pointer  :: col(:)        !< CSR format row indices
    sll_real64, pointer :: avals(:)      !< non zeros values
end type sll_csr_matrix

interface initialize
   module procedure initialize_csc_matrix
   module procedure initialize_csr_matrix
end interface initialize

interface todense
   module procedure csc_todense
   module procedure csr_todense
end interface todense

public initialize, todense

sll_int32 :: i, j, k, l

contains

subroutine initialize_csc_matrix(self, n, nnzeros, error)

   type(sll_csc_matrix)   :: self
   sll_int32, intent(in)  :: n         !< Matrix dimension
   sll_int32, intent(in)  :: nnzeros   !< Number of non-zeros
   sll_int32, intent(out) :: error

   self%n       = n
   self%nnzeros = nnzeros
   SLL_ALLOCATE(self%colptr(n+1), error)
   SLL_ALLOCATE(self%row(nnzeros), error)
   SLL_CLEAR_ALLOCATE(self%avals(1:nnzeros), error)

end subroutine initialize_csc_matrix

subroutine initialize_csr_matrix(this, n, nnzeros, error)

   type(sll_csr_matrix)  :: this
   sll_int32, intent(in) :: n         !< Matrix dimension
   sll_int32, intent(in) :: nnzeros   !< Number of non-zeros
   sll_int32             :: error

   this%n       = n
   this%nnzeros = nnzeros
   SLL_ALLOCATE(this%rowptr(n+1), error)
   SLL_ALLOCATE(this%col(nnzeros), error)
   SLL_ALLOCATE(this%avals(nnzeros), error)

end subroutine initialize_csr_matrix

subroutine csc_todense( this, dense_matrix)

   type(sll_csc_matrix)       :: this
   sll_real64, dimension(:,:) :: dense_matrix

   l = 0
   do j = 1, this%n 
      do k = this%colptr(j),this%colptr(j+1)-1 
         l = l + 1
         i = this%row(l)
         dense_matrix(i,j) = this%avals(l)
      end do
   end do

end subroutine csc_todense

subroutine csr_todense( this, dense_matrix)

   type(sll_csr_matrix)       :: this
   sll_real64, dimension(:,:) :: dense_matrix

   l = 0
   do i = 1, this%n 
      do k = this%rowptr(i),this%rowptr(i+1)-1 
         l = l + 1
         j = this%col(l)
         dense_matrix(i,j) = this%avals(l)
      end do
   end do

end subroutine csr_todense

end module sll_sparse
