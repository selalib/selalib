!   This test  produces the b-spline coeff.s  bcoef  of the spline of order
!   k  with knots  t(i), i=1,..., n + k , which takes on the value
!   gtau(i) at  tau(i), i=1,..., n .
!
!******  Given data  ******
!  tau.....array of length  n , containing data point abscissae.
!    a s s u m p t i o n . . .  tau  is strictly increasing
!  gtau.....corresponding array of length  n , containing data point or-
!        dinates
!  t.....knot sequence, of length  n+k
!  n.....number of data points and dimension of spline space  s(k,t)
!  k.....order of spline
!
!******  Computed  ******
!  q.....array of size  (2*k-1)*n , containing the triangular factoriz-
!        ation of the coefficient matrix of the linear system for the b-
!        coefficients of the spline interpolant.
!           the b-coeffs for the interpolant of an additional data set
!        (tau(i),htau(i)), i=1,...,n  with the same data abscissae can
!        be obtained without going through all the calculations in this
!        routine, simply by loading  htau  into  bcoef  and then execut-
!        ing the    call banslv ( q, 2*k-1, n, k-1, k-1, bcoef )
!  bcoef.....the b-coefficients of the interpolant, of length  n
!  iflag.....an integer indicating success (= 1)  or failure (= 2)
!        the linear system to be solved is (theoretically) invertible if
!        and only if
!              t(i) .lt. tau(i) .lt. t(i+k),    all i.
!        violation of this condition is certain to lead to  iflag = 2 .
!
!******  m e t h o d  ******
!     the i-th equation of the linear system  a*bcoef = b  for the b-co-
!  effs of the interpolant enforces interpolation at  tau(i), i=1,...,n.
!  hence,  b(i) = gtau(i), all i, and  a  is a band matrix with  2k-1
!   bands (if it is invertible).
!     the matrix  a  is generated row by row and stored, diagonal by di-
!  agonal, in the  r o w s  of the array  q , with the main diagonal go-
!  ing into row  k .  see comments in the program below.
!     the banded system is then solved by a call to  banfac (which con-
!  structs the triangular factorization for  a  and stores it again in
!   q ), followed by a call to  banslv (which then obtains the solution
!   bcoef  by substitution).
!     banfac  does no pivoting, since the total positivity of the matrix
!  a  makes this unnecessary.
!

program test_pppack
implicit none

integer, parameter :: n = 11
integer, parameter :: k = 3
real(8) :: bcoef(n)
real(8) :: gtau(n)
real(8) :: tau(n)
real(8) :: z(n)
real(8) :: taui
real(8) :: q( 2*k-1,n)
real(8) :: t(n+k)
real(8) :: a(n,n)

real(8), parameter :: tau_min = -0.0_8
real(8), parameter :: tau_max = +1.0_8

integer :: iflag
integer :: i
integer :: ilp1mx
integer :: j
integer :: jj
integer :: km1
integer :: kpkm2
integer :: left
integer :: lenq
integer :: np1
real(8) :: fact
real(8) :: alpha, beta, gamma

do i = 1, n
  tau(i)  = tau_min + (i-1)*(tau_max-tau_min)/(n-1)
  gtau(i) = tau(i)
end do

t(1:k) = tau_min
if ( mod(k,2) == 0 ) then
  do i = k+1,n
    t(i) = tau(i-k/2) 
  end do
else
  do i = k+1, n
    t(i) = 0.5*(tau(i-(k-1)/2)+tau(i-1-(k-1)/2))
  end do
end if
t(n+1:n+k) = tau_max
write(*,"(' tau = ', 11f7.3)") tau
write(*,"(' t   = ', 14f7.3)") t

np1 = n + 1
km1 = k - 1
kpkm2 = 2*km1
left = k
lenq = n*(k+km1)

q = 0.0_8

!  ***   loop over i to construct the  n  interpolation equations
do i=1,n

  taui = tau(i)
  ilp1mx = min(i+k,n+1)
  
  !  Find LEFT in the closed interval (I,I+K-1) such that
  !    T(LEFT) <= TAU(I) < T(LEFT+1)
  !  The matrix is singular if this is not possible.
  
  left = max(left,i)
  
  if ( taui < t(left) ) stop ' The linear system is not invertible!'

  do while ( t(left+1) <= taui )
    left = left + 1
    if ( left < ilp1mx ) cycle
    left = left - 1
    if ( t(left+1) < taui ) stop ' The linear system is not invertible!'
    exit
  end do

  write(*,"(i3,3f7.3)") left, t(left), taui, t(left+1)
  ! *** the i-th equation enforces interpolation at taui, hence
  ! a(i,j) = b(j,k,t)(taui), all j. only the  k  entries with  j =
  ! left-k+1,...,left actually might be nonzero. these  k  numbers
  ! are returned, in  bcoef (used for temp.storage here), by the
  ! following
  call bsplvb ( t, k, 1, taui, left, bcoef )

  ! we therefore want  bcoef(j) = b(left-k+j)(taui) to go into
  ! a(i,left-k+j), i.e., into  q(i-(left+j)+2*k,(left+j)-k) since
  ! a(i+j,j)  is to go into  q(i+k,j), all i,j,  if we consider  q
  ! as a two-dim. array , with  2*k-1  rows (see comments in
  ! banfac). in the present program, we treat  q  as an equivalent
  ! one-dimensional array (because of fortran restrictions on
  ! dimension statements) . we therefore want  bcoef(j) to go into
  ! entry
  !     i -(left+j) + 2*k + ((left+j) - k-1)*(2*k-1)
  !            =  i-left+1 + (left -k)*(2*k-1) + (2*k-2)*j
  ! of  q .
  jj = i-left+1 + (left-k)*(k+km1)
  do j=1,k
    !jj = jj+kpkm2
    !q(jj) = bcoef(j)
    q(i-(left+j)+2*k,(left+j)-k) = bcoef(j)
  end do

end do

a = 0.0_8
do j = 1, n
  do i = -k+1,k-1
    if( i+j >=1 .and. i+j <=n) a(i+j,j) = q(i+k,j)
  end do
  write(*,"(11f7.3)") a(:,j)
end do

alpha        = q(2*k-2,n)
beta         = q(2,1)
gamma        = - q(1,1) 

q(1,1)      = q(1,1) - gamma 
q(2*k-1,n)  = q(2*k-1,n)-alpha*beta/gamma 

call banfac ( q, k+km1, n, km1, km1, iflag )
bcoef = gtau
call banslv ( q, k+km1, n, km1, km1, bcoef )

z    = 0.0_8    
z(1) = gamma
z(n) = alpha

call banslv ( q, k+km1, n, km1, km1, z )

print*, "gamma=",gamma
fact=(bcoef(1)+beta*bcoef(n)/gamma)/(1.0+z(1)+beta*z(n)/gamma)

do i=1,n 
  bcoef(i) = bcoef(i) - fact*z(i) 
end do

write(*,"(' bcoef = ', 11f7.3)") bcoef


end program test_pppack
