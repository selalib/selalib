program gauss_legendre_tester
#include "sll_working_precision.h"
  use sll_constants
  use gauss_legendre_integration
  use test_function_module
  use sll_gausslobatto
  implicit none
  intrinsic :: dsin
  integer :: i,j
  type(gausslobatto1D) :: gausslob
  do i=2,10
     write (*,'(a, i8, a, f20.12)') 'case n = ', i, ': ', &
          gauss_legendre_integrate_1D( test_func, 0.0_f64, sll_pi/2.0, i)
  end do
  print *, 'Exact value: '
  write (*,'(f22.15)') 0.4674011002723395

  do i=2,10
     write (*,'(a, i8, a, f20.12)') 'case n = ', i, ': ', &
          gauss_legendre_integrate_1D( one, 0.0_f64, 1.0_f64, i)
  end do
  print *, 'Exact value: '
  write (*,'(f22.15)') 1.00000

  print *, 'Test gauss_points()'
  print *, gauss_points(5,-1.0_f64,1.0_f64)

  print*,'Test Gauss-Lobatto'
  do i=2,10
     call init_gausslobatto_1d(i,gausslob)
     !don't to it in real program, use transformation between real mesh and reference element
     !here it is done for simplicity
     gausslob%node(:)=(gausslob%node(:)+1.0d0)*(sll_pi/2.0d0)/2.0d0
     write (*,'(a, i8, a, e20.12)') 'case n = ', i, ': ', &
          & sum((/ (gausslob%weigh(j)*test_func(gausslob%node(j))*sll_pi/4.0d0,j=1,i) /))
     call delete_gausslobatto_1d(gausslob)
  end do

  print*,'Test Gauss-lobatto points and weight (5 points)'
  call init_gausslobatto_1d(5,gausslob)
  print*,gausslob%node
  print*,gausslob%weigh
  call delete_gausslobatto_1d(gausslob)

end program gauss_legendre_tester
