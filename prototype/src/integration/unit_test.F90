program integration_tester
#include "sll_working_precision.h"
#include "sll_constants.h"
  use gauss_legendre_integration
  use gauss_lobatto_integration
  use test_function_module
  implicit none
  intrinsic :: dsin
  integer :: i,j,n
  sll_real64, dimension(10) :: x, w
  sll_real64, dimension(:,:), allocatable :: d, dlag
  character(len=18) :: string

  write (*,'(9x, a10, 9x, a10, 5x, a20 )') ' legendre ',' lobatto ', 'Exact value: '
  do i=2,10
     write (*,'(a, i2, a, 3e20.12)') 'n = ', i, ': ', &
          gauss_legendre_integrate_1D( test_func, 0._f64, sll_pi/2.0, i), &
          gauss_lobatto_integrate_1D(  test_func, 0._f64, sll_pi/2.0, i), &
          0.4674011002723395
  end do

  do i=2,10
     write (*,'(a, i8, a, f20.12)') 'case n = ', i, ': ', &
          gauss_legendre_integrate_1D( one, 0.0_f64, 1.0_f64, i)
  end do
  print *, 'Exact value: '
  write (*,'(f22.15)') 1.00000

  print *, 'Test gauss_points()'
  print *, gauss_points(5,-1.0_f64,1.0_f64)

  x = gauss_lobatto_points( 10, -1._f64, 1._f64)
  w = gauss_lobatto_weights(10)

  do i = 1, 10
     write(*,"(2f20.15)") &
     x(i), w(i)
  end do

! sage: x = PolynomialRing(RealField(200),'x').gen()
! sage: n = 10
! sage: P=legendre_P(n-1,x)
! sage: L=P.derivative()
! sage: proots = L.roots()
! sage: xk = [-1]+[ proots[i][0] for i in range(n-2)]+[1]
! sage: wk =[2/(n*(n-1))]+[2/(n*(n-1)*(legendre_P(n-1,xk[i]))^2) for i in range(1,n-1)]+[2/(n*(n-1))]
! sage: for i in range(10):
! sage:    print " %28.15f %28.15f " % (xk[i], wk[i])

  print*, " ** exact values with sage"

  write(*,*) " -1.000000000000000  0.022222222222222  "
  write(*,*) " -0.919533908166459  0.133305990851070  "
  write(*,*) " -0.738773865105505  0.224889342063126  "
  write(*,*) " -0.477924949810444  0.292042683679684  "
  write(*,*) " -0.165278957666387  0.327539761183897  "
  write(*,*) "  0.165278957666387  0.327539761183897  "
  write(*,*) "  0.477924949810444  0.292042683679684  "
  write(*,*) "  0.738773865105505  0.224889342063126  "
  write(*,*) "  0.919533908166459  0.133305990851070  "
  write(*,*) "  1.000000000000000  0.022222222222222  "

  write(*,"(/,a,/)") "Matrix of derivatives"
  n = 4
  allocate(d(n,n))
  d = gauss_lobatto_derivative_matrix(n, -1._f64, 1._f64) 

  write (string, '( "(",I2,"f20.15)" )' )  n
  do i = 1, n
     write(*,string) ( d(i,j), j = 1, n)
  end do 

  allocate(dlag(4,4))
  dlag(1,1) = -0.3000000000000000000000000000000000000000D1
  dlag(1,2) = -0.8090169943749474241022934171828190588602D0
  dlag(1,3) = 0.3090169943749474241022934171828190588602D0
  dlag(1,4) = -0.5000000000000000000000000000000000000000D0
  dlag(2,1) = 0.4045084971874737120511467085914095294301D1
  dlag(2,2) = -0.6D-39
  dlag(2,3) = -0.1118033988749894848204586834365638117720D1
  dlag(2,4) = 0.1545084971874737120511467085914095294300D1
  dlag(3,1) = -0.1545084971874737120511467085914095294300D1
  dlag(3,2) = 0.1118033988749894848204586834365638117720D1
  dlag(3,3) = 0.1D-38
  dlag(3,4) = -0.4045084971874737120511467085914095294301D1
  dlag(4,1) = 0.5000000000000000000000000000000000000000D0
  dlag(4,2) = -0.3090169943749474241022934171828190588602D0
  dlag(4,3) = 0.8090169943749474241022934171828190588602D0
  dlag(4,4) = 0.3000000000000000000000000000000000000000D1

  write(*,"(/,a,/)") " Exact values with maple "

  do i = 1, n
     write(*,string) ( dlag(i,j), j = 1, n)
  end do 

  print*, 'PASSED'

end program integration_tester
