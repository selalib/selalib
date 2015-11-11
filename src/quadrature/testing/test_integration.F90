program test_integration

#include "sll_memory.h"
#include "sll_working_precision.h"

use sll_m_rectangle_integration
use sll_m_trapz_integration
use sll_m_gauss_legendre_integration
use sll_m_gauss_lobatto_integration
use sll_m_fekete_integration
use sll_m_box_splines, only: &
     write_connectivity
use test_function_module, only: &
     one, &
     test_func, &
     one_2D, &
     test_func_2D
use sll_m_constants, only : &
     sll_pi
implicit none

integer :: i,j,n
sll_int32 :: ierr
!sll_int32 :: degree

sll_real64, dimension(10) :: x, w
sll_real64, dimension(:,:), allocatable :: d
sll_real64, dimension(:,:), allocatable :: dlag

character(len=18) :: string

! For the fekete quadrature:
sll_real64, dimension(2, 3) :: pxy1
sll_real64, dimension(2, 3) :: pxy2
sll_real64, dimension(:,:), allocatable :: xyw
!sll_real64 :: app_res
sll_int32  :: rule
!sll_int32  :: num_cells
!type(sll_hex_mesh_2d), pointer :: mesh
!sll_real64, dimension(:,:), allocatable     :: knots
!sll_int32,  dimension(:,:), allocatable     :: LM


write (*,'(5x, 5a16 )') &
'rectangle','trapz','legendre','lobatto', 'Exact value'
do i=2,10
  do j = 1, i
    x(j) = (j-1)*0.5_f64*sll_pi/(i-1)
  end do
  write (*,'(a, i2, a, 5f16.12)') 'n = ', i, ': ', &
   rectangle_integrate_1d( test_func, x, i), &
   trapz_integrate_1d( test_func, x, i), &
   gauss_legendre_integrate_1d( test_func, 0._f64, sll_pi/2._f64, i), &
   gauss_lobatto_integrate_1d(  test_func, 0._f64, sll_pi/2._f64, i), &
   0.4674011002723395
end do

write(*,*)" ------- "

do i=2,10
   do j = 1, i
     x(j) = (j-1)*1.0_f64/(i-1)
   end do
   write (*,'(a, i2, a, 4f16.12)') 'n = ', i, ': ', &
        rectangle_integrate_1d( one, x, i), &
        trapz_integrate_1d( one, x, i), &
        gauss_legendre_integrate_1d( one, 0.0_f64, 1.0_f64, i), &
        gauss_lobatto_integrate_1d( one, 0.0_f64, 1.0_f64, i)
end do
print *, 'Exact value: '
write (*,'(f22.15)') 1.00000

print *, 'Test gauss_points()'
print *, gauss_legendre_points_and_weights(5,-1.0_f64,1.0_f64)

x = gauss_lobatto_points( 10, -1._f64, 1._f64)
w = gauss_lobatto_weights(10)

print*, 'Gauss Lobatto points and weights in [-1:1]'
do i = 1, 10
   write(*,"(2f20.15)") x(i), w(i)
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
d = gauss_lobatto_derivative_matrix(n) 

write (string, '( "(",I2,"f20.15)" )' )  n
do i = 1, n
   write(*,string) ( d(i,j), j = 1, n)
end do

allocate(dlag(4,4))
dlag(1,1) = -0.3000000000000000000000000000000000000000e1_f64
dlag(1,2) = -0.8090169943749474241022934171828190588602e0_f64
dlag(1,3) = 0.3090169943749474241022934171828190588602e0_f64
dlag(1,4) = -0.5000000000000000000000000000000000000000e0_f64
dlag(2,1) = 0.4045084971874737120511467085914095294301e1_f64
dlag(2,2) = -0.6e-39_f64
dlag(2,3) = -0.1118033988749894848204586834365638117720e1_f64
dlag(2,4) = 0.1545084971874737120511467085914095294300e1_f64
dlag(3,1) = -0.1545084971874737120511467085914095294300e1_f64
dlag(3,2) = 0.1118033988749894848204586834365638117720e1_f64
dlag(3,3) = 0.1e-38_f64
dlag(3,4) = -0.4045084971874737120511467085914095294301e1_f64
dlag(4,1) = 0.5000000000000000000000000000000000000000e0_f64
dlag(4,2) = -0.3090169943749474241022934171828190588602e0_f64
dlag(4,3) = 0.8090169943749474241022934171828190588602e0_f64
dlag(4,4) = 0.3000000000000000000000000000000000000000e1_f64

write(*,"(/,a,/)") " Exact values with maple "

do i = 1, n
   write(*,string) ( dlag(i,j), j = 1, n)
end do

write(*,"(/,a)") "*********************************** "
write(*,"(a)") "       FEKETE QUAD TEST       "
write(*,"(a)") "*********************************** "

!Definition of first triangle
pxy1(:,1) = (/ 0._f64, 0._f64 /)
pxy1(:,2) = (/ 1._f64, 0._f64 /)
pxy1(:,3) = (/ 0._f64, 1._f64 /)

!Definition of first triangle
pxy2(:,1) = (/ 1._f64, 0._f64 /)
pxy2(:,2) = (/ 1._f64, 1._f64 /)
pxy2(:,3) = (/ 0._f64, 1._f64 /)

rule = 2
call fekete_order_num ( rule, n )
SLL_ALLOCATE(xyw(1:3, 1:n), ierr)

write(*,"(a)") " Computing Fekete points and weights on reference triangle "
write(*,"(/,a)") "           x                   y                    w"
xyw = fekete_points_and_weights(pxy1, rule)

do j = 1, n
   write(*, string) (xyw(i,j), i = 1, 3)
end do

print *, "sum weights = ", SUM(xyw(3,:))
! write(*,"(/,a)") " --Test for a constant real function (=1) "
! write(*,"(a)") "    on the squared domain [0,1]^2 divided on 2 triangles "

! app_res = fekete_integral(one_2D, pxy1) + fekete_integral(one_2D, pxy2)

! write (*,"(a, f20.12, a ,/)") " aprox = ", app_res, " (expected = 1.)"


print*, 'PASSED'

end program test_integration
