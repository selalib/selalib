program test_sparse
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_constants.h"
#include "sll_utilities.h"

use sll_sparse

implicit none

type(sll_csc_matrix) :: csc_matrix
type(sll_csr_matrix) :: csr_matrix
sll_real64, dimension(:,:), allocatable :: dense_matrix

sll_int32 :: n, nnzeros
sll_int32 :: i, j
sll_int32 :: error

n = 6
SLL_CLEAR_ALLOCATE(dense_matrix(1:n,1:n), error)

nnzeros = 3*n+2

call initialize(csc_matrix, n, nnzeros, error)

!Lower matrix only
j=1
do i = 1, n
   csc_matrix%colptr(i) = j
   csc_matrix%row(j)    = i
   csc_matrix%avals(j) = 2
   j=j+1
   if (i /= n) then
      csc_matrix%row(j)   = i+1
      csc_matrix%avals(j) = -1.
      j=j+1
   end if
end do
csc_matrix%colptr(n+1) = j

call sll_display(csc_matrix%colptr, "i4")
call sll_display(csc_matrix%row, "i4")
call sll_display(csc_matrix%avals, "f8.3")
call todense(csc_matrix, dense_matrix)
call sll_display(dense_matrix, "f8.3")

call initialize(csr_matrix, n, nnzeros, error)

!Full symetrix matrix
i=1
do j = 1, n
   csr_matrix%rowptr(j) = i
   if (j > 1) then
      csr_matrix%avals(i)  = -1
      csr_matrix%col(i)    = j-1
      i = i+1
   end if
   csr_matrix%col(i) = j
   csr_matrix%avals(i) = 2
   i=i+1
   if (j < n) then
      csr_matrix%col(i)   = j+1
      csr_matrix%avals(i) = -1.
      i=i+1
   end if
end do
csr_matrix%rowptr(n+1) = i

call sll_display(csr_matrix%rowptr, "i4")
call sll_display(csr_matrix%col, "i4")
call sll_display(csr_matrix%avals, "f8.3")
call todense(csr_matrix, dense_matrix)
call sll_display(dense_matrix, "f8.3")

end program test_sparse
