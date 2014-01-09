function bvalue ( t, bcoef, n, k, x, jderiv )

!*****************************************************************************80
!
!! BVALUE evaluates a derivative of a spline from its B-spline representation.
!
!  Discussion:
!
!    The spline is taken to be continuous from the right.
!
!    The nontrivial knot interval (T(I),T(I+1)) containing X is
!    located with the aid of INTERV.  The K B-spline coefficients
!    of F relevant for this interval are then obtained from BCOEF,
!    or are taken to be zero if not explicitly available, and are
!    then differenced JDERIV times to obtain the B-spline
!    coefficients of (D**JDERIV)F relevant for that interval.
!
!    Precisely, with J = JDERIV, we have from X.(12) of the text that:
!
!      (D**J)F = sum ( BCOEF(.,J)*B(.,K-J,T) )
!
!    where
!                      / BCOEF(.),                    if J == 0
!                     /
!       BCOEF(.,J) = / BCOEF(.,J-1) - BCOEF(.-1,J-1)
!                   / -----------------------------,  if 0 < J
!                  /    (T(.+K-J) - T(.))/(K-J)
!
!    Then, we use repeatedly the fact that
!
!      sum ( A(.) * B(.,M,T)(X) ) = sum ( A(.,X) * B(.,M-1,T)(X) )
!
!    with
!                   (X - T(.))*A(.) + (T(.+M-1) - X)*A(.-1)
!      A(.,X) =   ---------------------------------------
!                   (X - T(.))      + (T(.+M-1) - X)
!
!    to write (D**J)F(X) eventually as a linear combination of
!    B-splines of order 1, and the coefficient for B(I,1,T)(X)
!    must then be the desired number (D**J)F(X).
!    See Chapter X, (17)-(19) of text.
!
!  Modified:
!
!    14 February 2007
!
!  Author:
!
!    Carl DeBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(N+K), the knot sequence.  T is assumed
!    to be nondecreasing.
!
!    Input, real ( kind = 8 ) BCOEF(N), B-spline coefficient sequence.
!
!    Input, integer N, the length of BCOEF.
!
!    Input, integer K, the order of the spline.
!
!    Input, real ( kind = 8 ) X, the point at which to evaluate.
!
!    Input, integer JDERIV, the order of the derivative to
!    be evaluated.  JDERIV is assumed to be zero or positive.
!
!    Output, real ( kind = 8 ) BVALUE, the value of the (JDERIV)-th
!    derivative of the spline at X.
!
  implicit none

  integer k
  integer n

  real ( kind = 8 ),dimension(:),pointer:: aj !(k)
  real ( kind = 8 ),dimension(:),pointer:: bcoef!(n)
  real ( kind = 8 ) bvalue
  real ( kind = 8 ):: tmp_value
  real ( kind = 8 ),dimension(:),pointer:: dl!(k)
  real ( kind = 8 ),dimension(:),pointer:: dr!(k)
  integer i
  integer ilo
  integer j
  integer jc
  integer jcmax
  integer jcmin
  integer jderiv
  integer jj
  integer mflag
  real ( kind = 8 ),dimension(:),pointer:: t!(n+k)
  real ( kind = 8 ) x

  bvalue = 0.0_8

  allocate(aj(k))
  allocate(dl(k))
  allocate(dr(k))

  aj(:)=0.0_8
  dl(:)=0.0_8
  dr(:)=0.0_8
  
  if ( k <= jderiv ) then
    return
  end if
!
!  Find I so that 1 <= I < N+K and T(I) < T(I+1) and T(I) <= X < T(I+1).
!
!  If no such I can be found, X lies outside the support of the
!  spline F and  BVALUE = 0.  The asymmetry in this choice of I makes F
!  right continuous.
!
  call interv ( t, n+k, x, i, mflag )

  if ( mflag /= 0 ) then
    return
  end if
!
!  If K = 1 (and JDERIV = 0), BVALUE = BCOEF(I).
!
  if ( k <= 1 ) then
    bvalue = bcoef(i)
    return
  end if
!
!  Store the K B-spline coefficients relevant for the knot interval
!  ( T(I),T(I+1) ) in AJ(1),...,AJ(K) and compute DL(J) = X - T(I+1-J),
!  DR(J) = T(I+J)-X, J=1,...,K-1.  Set any of the AJ not obtainable
!  from input to zero.
!
!  Set any T's not obtainable equal to T(1) or to T(N+K) appropriately.
!
  jcmin = 1

  if ( k <= i ) then

    do j = 1, k-1
      dl(j) = x - t(i+1-j)
    end do

  else

    jcmin = 1 - ( i - k )

    do j = 1, i
      dl(j) = x - t(i+1-j)
    end do

    do j = i, k-1
      aj(k-j) = 0.0_8
      dl(j) = dl(i)
    end do

  end if

  jcmax = k

  if ( n < i ) then

    jcmax = k + n - i
    do j = 1, k + n - i
      dr(j) = t(i+j) - x
    end do

    do j = k+n-i, k-1
      aj(j+1) = 0.0_8
      dr(j) = dr(k+n-i)
    end do

  else

    do j = 1, k-1
      dr(j) = t(i+j) - x
    end do

  end if

  do jc = jcmin, jcmax
    aj(jc) = bcoef(i-k+jc)
  end do
!
!  Difference the coefficients JDERIV times.
!
  do j = 1, jderiv

    ilo = k - j
    do jj = 1, k - j
      aj(jj) = ( ( aj(jj+1) - aj(jj) ) / ( dl(ilo) + dr(jj) ) ) &
        * real ( k - j, kind = 8 )
      ilo = ilo - 1
    end do

  end do
!
!  Compute value at X in (T(I),T(I+1)) of JDERIV-th derivative,
!  given its relevant B-spline coefficients in AJ(1),...,AJ(K-JDERIV).
!
  do j = jderiv+1, k-1
    ilo = k-j
    do jj = 1, k-j
      aj(jj) = ( aj(jj+1) * dl(ilo) + aj(jj) * dr(jj) ) &
        / ( dl(ilo) + dr(jj) )
      ilo = ilo - 1
    end do
  end do

  bvalue = aj(1)

  return
  deallocate(aj)
  deallocate(dl)
  deallocate(dr)


end function bvalue
