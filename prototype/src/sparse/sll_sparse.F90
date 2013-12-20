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
end interface initialize

public initialize

contains

subroutine initialize_csc_matrix(this, n, nnzeros, error)

   type(sll_csc_matrix)  :: this
   sll_int32, intent(in) :: n         !< Matrix dimension
   sll_int32             :: nnzeros   !< Number of non-zeros
   sll_int32             :: error

   this%n       = n
   this%nnzeros = nnzeros
   SLL_ALLOCATE(this%colptr(n+1), error)
   SLL_ALLOCATE(this%row(nnzeros), error)
   SLL_ALLOCATE(this%avals(nnzeros), error)

end subroutine initialize_csc_matrix

subroutine initialize_csr_matrix(this, n, nnzeros, error)

   type(sll_csr_matrix)  :: this
   sll_int32, intent(in) :: n         !< Matrix dimension
   sll_int32             :: nnzeros   !< Number of non-zeros
   sll_int32             :: error

   this%n       = n
   this%nnzeros = nnzeros
   SLL_ALLOCATE(this%rowptr(n+1), error)
   SLL_ALLOCATE(this%col(nnzeros), error)
   SLL_ALLOCATE(this%avals(nnzeros), error)

end subroutine initialize_csr_matrix

end module sll_sparse
