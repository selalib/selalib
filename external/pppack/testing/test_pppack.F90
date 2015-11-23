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

integer, parameter :: n = 8
integer, parameter :: k = 3
integer :: iflag, i,ilp1mx,j,jj,km1,kpkm2,left,lenq,np1
real(8) :: bcoef(n),gtau(n),tau(n),   taui
real(8) :: q((2*k-1)*n), t(n+k)

np1 = n + 1
km1 = k - 1
kpkm2 = 2*km1
left = k
lenq = n*(k+km1)

q = 0.0_8

!  ***   loop over i to construct the  n  interpolation equations
do i=1,n

  taui = tau(i)
  ilp1mx = min0(i+k,np1)
!        *** find  left  in the closed interval (i,i+k-1) such that
!                t(left) .le. tau(i) .lt. t(left+1)
!        matrix is singular if this is not possible
  left = max(left,i)
  do while (left .lt. ilp1mx .and. taui < t(left+1)) 
    left = left+1
  end do

  if (taui .lt. t(left+1)) then 
    ! *** the i-th equation enforces interpolation at taui, hence
    ! a(i,j) = b(j,k,t)(taui), all j. only the  k  entries with  j =
    ! left-k+1,...,left actually might be nonzero. these  k  numbers
    ! are returned, in  bcoef (used for temp.storage here), by the
    ! following
    call bsplvb ( t, k, 1, taui, left, bcoef )

  else 

    stop 'linear system is not invertible'

  end if
     
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
    jj = jj+kpkm2
    q(jj) = bcoef(j)
  end do

end do

! ***obtain factorization of  a  , stored again in  q.

call banfac ( q, k+km1, n, km1, km1, iflag )

! *** solve  a*bcoef = gtau  by backsubstitution
bcoef = gtau

call banslv ( q, k+km1, n, km1, km1, bcoef )

end program test_pppack
