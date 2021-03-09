program test_gauss_lobatto

#include "sll_memory.h"
#include "sll_working_precision.h"

   use sll_m_constants, only: sll_p_pi

   use sll_m_gauss_lobatto_integration, only: &
      sll_f_gauss_lobatto_derivative_matrix, &
      sll_o_gauss_lobatto_integrate_1d, &
      sll_f_gauss_lobatto_points, &
      sll_f_gauss_lobatto_weights

   use test_function_module

   implicit none

   sll_int32  :: i, j, n
   sll_real64 :: s

   sll_real64, dimension(10) :: x, w
   sll_real64, dimension(:, :), allocatable :: d
   sll_real64, dimension(:, :), allocatable :: dlag

   character(len=18) :: string

   write (*, '(5x, 2a16 )') 'lobatto', 'Exact value'
   do i = 2, 10
      do j = 1, i
         x(j) = real(j - 1, f64)*0.5_f64*sll_p_pi/real(i - 1, f64)
      end do
      write (*, '(a, i2, a, 2f16.12)') 'n = ', i, ': ', &
         sll_o_gauss_lobatto_integrate_1d(test_func, 0._f64, sll_p_pi/2._f64, i), &
         0.4674011002723395_f64
   end do

   write (*, *) " ------- "

   do i = 2, 10
      do j = 1, i
         x(j) = real(j - 1, f64)*1.0_f64/real(i - 1, f64)
      end do
      write (*, '(a, i2, a, 1f16.12)') 'n = ', i, ': ', &
         sll_o_gauss_lobatto_integrate_1d(one, 0.0_f64, 1.0_f64, i)
   end do
   print *, 'Exact value: '
   write (*, '(f22.15)') 1.00000

   x = sll_f_gauss_lobatto_points(10, -1._f64, 1._f64)
   w = sll_f_gauss_lobatto_weights(10)

   print *, 'Gauss Lobatto points and weights in [-1:1]'
   do i = 1, 10
      write (*, "(2f20.15)") x(i), w(i)
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

   print *, " ** exact values with sage"

   write (*, *) " -1.000000000000000  0.022222222222222  "
   write (*, *) " -0.919533908166459  0.133305990851070  "
   write (*, *) " -0.738773865105505  0.224889342063126  "
   write (*, *) " -0.477924949810444  0.292042683679684  "
   write (*, *) " -0.165278957666387  0.327539761183897  "
   write (*, *) "  0.165278957666387  0.327539761183897  "
   write (*, *) "  0.477924949810444  0.292042683679684  "
   write (*, *) "  0.738773865105505  0.224889342063126  "
   write (*, *) "  0.919533908166459  0.133305990851070  "
   write (*, *) "  1.000000000000000  0.022222222222222  "

   write (*, "(/,a,/)") "Matrix of derivatives"
   n = 4
   allocate (d(n, n))
   d = sll_f_gauss_lobatto_derivative_matrix(n)

   write (string, '( "(",I2,"f20.15)" )') n
   do i = 1, n
      write (*, string) (d(i, j), j=1, n)
   end do

   allocate (dlag(4, 4))
   dlag(1, 1) = -3.000000000000000_f64
   dlag(1, 2) = -0.809016994374947_f64
   dlag(1, 3) = 0.309016994374947_f64
   dlag(1, 4) = -0.500000000000000_f64
   dlag(2, 1) = 4.045084971874737_f64
   dlag(2, 2) = -0.6d-39
   dlag(2, 3) = -1.118033988749894_f64
   dlag(2, 4) = 1.545084971874737_f64
   dlag(3, 1) = -1.545084971874737_f64
   dlag(3, 2) = 1.118033988749894_f64
   dlag(3, 3) = 0.1d-38
   dlag(3, 4) = -4.045084971874737_f64
   dlag(4, 1) = 0.500000000000000_f64
   dlag(4, 2) = -0.309016994374947_f64
   dlag(4, 3) = 0.809016994374947_f64
   dlag(4, 4) = 3.000000000000000_f64

   write (*, "(/,a,/)") " Exact values with maple "

   do i = 1, n
      write (*, string) (dlag(j, i), j=1, n)
   end do

   s = sum(abs(transpose(dlag) - d))
   if (s > 1d-7) then
      print *, 'FAILED'
   else
      print *, 'PASSED'
   end if

end program test_gauss_lobatto
