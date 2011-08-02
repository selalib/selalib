program gauss_legendre_tester
#include "sll_working_precision.h"
  use numeric_constants
  use gauss_legendre_integration
  use test_function_module
  implicit none
  intrinsic :: dsin
  integer :: i
  do i=2,10
     write (*,'(a, i8, a, e20.12)') 'case n = ', i, ': ', &
          gauss_legendre_integrate_1D( test_func, 0.0_f64, sll_pi/2.0, i)
  end do
  print *, 'Exact value: '
  write (*,'(e20.15)') 0.4674011002723395


end program gauss_legendre_tester
