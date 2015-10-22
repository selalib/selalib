program example_tridiag
  use sll_tridiagonal
#include "sll_memory.h"
#include "sll_working_precision.h"
  implicit none

  !Declaration
  sll_real64, dimension(:), pointer     :: a
  sll_real64, allocatable, dimension(:) :: x
  sll_real64, allocatable, dimension(:) :: b
  sll_real64, dimension(:), pointer     :: cts
  sll_int32,  dimension(:), pointer     :: ipiv
  sll_int32                             :: n, ierr

  !Define the size of the problem
  !We take A=3x3
  n = 3
  
  !Allocate the memory
  SLL_ALLOCATE(a(3*n),ierr)
  SLL_ALLOCATE(b(n),ierr)
  SLL_ALLOCATE(x(n),ierr)
  SLL_ALLOCATE(cts(7*n),ierr)
  SLL_ALLOCATE(ipiv(n),ierr)

  !Fill a
  a(2:3*n-1) = 1.0
  !Fill b
  b(:) = 1.0
  
  ! Solve ax=b and put the result in x
  ! You can change x by b for use only one vector.
  call setup_cyclic_tridiag( a, n, cts, ipiv ) !Compute the factoriazation of A=LU
  call solve_cyclic_tridiag( cts, ipiv, b, n, x ) !Solve LUx=b

  print *, 'We solve ax=b with:'
  print *, 'a=(1 1 0)'
  print *, '  (1 1 1)'
  print *, '  (0 1 1)'
  print *, 'b=(1)'
  print *, '  (1)'
  print *, '  (1)'
  print *, 'We obtain x=(',INT(x),')'
  
  SLL_DEALLOCATE(a,ierr)
  SLL_DEALLOCATE(cts,ierr)  
  SLL_DEALLOCATE(ipiv,ierr)  
  SLL_DEALLOCATE_ARRAY(b,ierr)  
  SLL_DEALLOCATE_ARRAY(x,ierr)  

end program example_tridiag

