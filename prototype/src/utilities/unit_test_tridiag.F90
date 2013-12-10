program test_tridiag
  use sll_tridiagonal
#include "sll_memory.h"
#include "sll_working_precision.h"
  implicit none
  
  sll_real64, dimension(:), pointer :: a
  sll_real64, allocatable, dimension(:) :: x
  sll_real64, allocatable :: b(:)
  sll_comp64, allocatable, dimension(:) :: x_complex
  sll_comp64, allocatable :: b_complex(:)
  sll_real64, dimension(:), pointer :: cts
  sll_int32,  dimension(:), pointer  :: ipiv
  sll_real64 :: resid
  sll_comp64 :: resid_complex
  sll_int32 :: i, n, t, min_n, max_n, n_test
  sll_int32 :: err
  sll_int32 :: status=0
  
  min_n = 3
  max_n = 100
  n_test = 10
  
  print *, 'Tridiagonal setup/solve tester'
  print *, 'Allocating arrays'
  
  
  SLL_ALLOCATE(a(3*max_n),err)
  SLL_ALLOCATE(x(max_n),err)
  SLL_ALLOCATE(b(max_n),err)
  SLL_ALLOCATE(cts(7*max_n),err)
  SLL_ALLOCATE(ipiv(max_n),err)

  SLL_ALLOCATE(x_complex(max_n),err)
  SLL_ALLOCATE(b_complex(max_n),err)


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
        write (*,'(a, i12, a, 20es20.10)') '     for tridiagonal system of size ', n, ' the average value of the residual per point is ', resid
        if( abs(resid) > 1.0e-11 ) then
           status=1
           write (*,'(a, i2, a, 20es20.10)') 'TOO HIGH residual value found for size ', n, ' :', resid
        end if
     end do
  end do

!test for complex arrays

  do t=1,n_test
     do n=min_n,max_n
        call random_number(a)
        call random_number(b)
        do i=1,max_n
          b_complex(i) = cmplx(b(i),1._f64)
        enddo
 !       print *, 'array a:'
 !       print *, a(:)
 !       print *, 'array b:'
 !       print *, b(:)
        call setup_cyclic_tridiag( a, n, cts, ipiv )
 !       print *, 'cts:'
 !       print *, cts(:)
 !       print *, 'ipiv:'
 !       print *, ipiv(:)
        call solve_cyclic_tridiag( cts, ipiv, b_complex, n, x_complex )
 !       print *, 'x:'
 !       print *, x(:)
        ! compute the residual
        i = 1
        b_complex(i) = b_complex(i) - (a(i)*x_complex(n) + a(i+1)*x_complex(i) + a(i+2)*x_complex(i+1))
        do i=2,n-1
           b_complex(i) = b_complex(i) - (a(3*i)*x_complex(i+1) + a(3*i-1)*x_complex(i) + a(3*i-2)*x_complex(i-1))
        end do
        i = n
        b_complex(i) = b_complex(i) - (a(3*i)*x_complex(1) + a(3*i-1)*x_complex(i) + a(3*i-2)*x_complex(i-1))

        resid_complex = 0.0
        do i=1,n
           resid_complex = resid_complex + b_complex(i)*conjg(b_complex(i))
        end do
        resid_complex = sqrt(resid_complex/n)
        write (*,'(a, i12, a, 20es20.10)') '     for tridiagonal system of size ', n, ' the average value of the residual per point is ', resid_complex
        if( abs(resid_complex) > 1.0e-11 ) then
           status=1
           write (*,'(a, i2, a, 20es20.10)') 'TOO HIGH residual value found for size ', n, ' :', resid_complex
        end if
     end do
  end do









  SLL_DEALLOCATE(a,    err)
  SLL_DEALLOCATE(cts,  err)
  SLL_DEALLOCATE(ipiv, err)
  print *, ' '
  if(status==0) then
     print *, 'setup/solve_tridiag: PASS'
  else
     print *, 'setup/solve_tridag: FAIL'
  end if
end program test_tridiag

