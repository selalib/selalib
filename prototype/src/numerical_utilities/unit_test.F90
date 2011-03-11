program test_tridiag
  use sll_tridiagonal
#include "sll_memory.h"
#include "sll_working_precision.h"
  implicit none
  
  sll_real64, dimension(:), pointer :: a
  sll_real64, dimension(:), pointer :: x
  sll_real64, dimension(:), pointer :: b
  sll_real64, dimension(:), pointer :: cts
  sll_int32,  dimension(:), pointer  :: ipiv
  sll_real64 :: resid
  sll_int32 :: i, n, t, min_n, max_n, n_test
  sll_int32 :: err
  sll_int32 :: status=0
  
  min_n = 5
  max_n = 100
  n_test = 10
  
  print *, 'Tridiagonal setup/solve tester'
  print *, 'Allocating arrays'

  SLL_ALLOCATE(a(3*max_n),err)
  SLL_ALLOCATE(x(max_n),err)
  SLL_ALLOCATE(b(max_n),err)
  SLL_ALLOCATE(cts(7*max_n),err)
  SLL_ALLOCATE(ipiv(max_n),err)

  print *, 'entering main cycle...'  
  do t=1,n_test
     do n=min_n,max_n
        call random_number(a)
        call random_number(b)
 !       print *, 'array a:'
 !       print *, a(:)
 !       print *, 'array b:'
 !       print *, b(:)
        call setup_cyclic_tridiag( a, n, cts, ipiv )
 !       print *, 'cts:'
 !       print *, cts(:)
 !       print *, 'ipiv:'
 !       print *, ipiv(:)
        call solve_cyclic_tridiag( cts, ipiv, b, n, x )
 !       print *, 'x:'
 !       print *, x(:)
        ! compute the residual
        i = 1
        b(i) = b(i) - (a(i)*x(n) + a(i+1)*x(i) + a(i+2)*x(i+1))
        do i=2,n-1
           b(i) = b(i) - (a(3*i)*x(i+1) + a(3*i-1)*x(i) + a(3*i-2)*x(i-1))
        end do
        i = n
        b(i) = b(i) - (a(3*i)*x(1) + a(3*i-1)*x(i) + a(3*i-2)*x(i-1))

        resid = 0.0
        do i=1,n
           resid = resid + b(i)*b(i)
        end do
        resid = sqrt(resid/n)
        write (*,'(a, i12, 20es20.10)') 'average value of residual ', n, resid
        if( resid > 1.0e-10 ) then
           status=1
           write (*,'(a, i12, 20es20.10)') 'TOO HIGH residual value ', n, resid
        end if
     end do
  end do
  SLL_DEALLOCATE(a,    err)
  SLL_DEALLOCATE(x,    err)
  SLL_DEALLOCATE(b,    err)
  SLL_DEALLOCATE(cts,  err)
  SLL_DEALLOCATE(ipiv, err)
  if(status==0) then
     print *, 'setup/solve_tridiag: PASS'
  else
     print *, 'setup/solve_tridag: FAIL'
  end if
end program test_tridiag

