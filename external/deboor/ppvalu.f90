function ppvalu ( break, coef, l, k, x, jderiv )

!*****************************************************************************80
!
!! PPVALU evaluates a piecewise polynomial function or its derivative.
!
!  Discussion:
!
!    PPVALU calculates the value at X of the JDERIV-th derivative of
!    the piecewise polynomial function F from its piecewise
!    polynomial representation.
!
!    The interval index I, appropriate for X, is found through a
!    call to INTERV.  The formula for the JDERIV-th derivative
!    of F is then evaluated by nested multiplication.
!
!    The J-th derivative of F is given by:
!      (d^J) F(X) = 
!        COEF(J+1,I) + H * (
!        COEF(J+2,I) + H * (
!        ...
!        COEF(K-1,I) + H * (
!        COEF(K,  I) / (K-J-1) ) / (K-J-2) ... ) / 2 ) / 1
!    with
!      H = X - BREAK(I)
!    and
!      I = max ( 1, max ( J, BREAK(J) <= X, 1 <= J <= L ) ).
!
!  Modified:
!
!    16 February 2007
!
!  Author:
!
!    Carl de Boor
!
!  Reference:
!
!    Carl de Boor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663,
!    LC: QA1.A647.v27.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) BREAK(L+1), real COEF(*), integer L, the
!    piecewise polynomial representation of the function F to be evaluated.
!
!    Input, integer ( kind = 4 ) K, the order of the polynomial pieces that 
!    make up the function F.  The usual value for K is 4, signifying a 
!    piecewise cubic polynomial.
!
!    Input, real ( kind = 8 ) X, the point at which to evaluate F or
!    of its derivatives.
!
!    Input, integer ( kind = 4 ) JDERIV, the order of the derivative to be
!    evaluated.  If JDERIV is 0, then F itself is evaluated,
!    which is actually the most common case.  It is assumed
!    that JDERIV is zero or positive.
!
!    Output, real ( kind = 8 ) PPVALU, the value of the JDERIV-th
!    derivative of F at X.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) l

  real ( kind = 8 ) break(l+1)
  real ( kind = 8 ) coef(k,l)
  real ( kind = 8 ) fmmjdr
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) jderiv
  integer ( kind = 4 ) m
  integer ( kind = 4 ) ndummy
  real ( kind = 8 ) ppvalu
  real ( kind = 8 ) value
  real ( kind = 8 ) x

  value = 0.0D+00

  fmmjdr = k - jderiv
!
!  Derivatives of order K or higher are identically zero.
!
  if ( k <= jderiv ) then
    return
  end if
!
!  Find the index I of the largest breakpoint to the left of X.
!
  call interv ( break, l+1, x, i, ndummy )
!
!  Evaluate the JDERIV-th derivative of the I-th polynomial piece at X.
!
  h = x - break(i)
  m = k
 
  do

    value = ( value / fmmjdr ) * h + coef(m,i)
    m = m - 1
    fmmjdr = fmmjdr - 1.0D+00

    if ( fmmjdr <= 0.0D+00 ) then
      exit
    end if

  end do

  ppvalu = value
 
  return
end function ppvalu

