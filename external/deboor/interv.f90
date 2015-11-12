subroutine interv ( xt, lxt, x, left, mflag )

!*****************************************************************************80
!
!! INTERV brackets a real value in an ascending vector of values.
!
!  Discussion:
!
!    The XT array is a set of increasing values.  The goal of the routine
!    is to determine the largest index I so that 
!
!      XT(I) < XT(LXT)  and  XT(I) <= X.
!
!    The routine is designed to be efficient in the common situation
!    that it is called repeatedly, with X taken from an increasing
!    or decreasing sequence.
!
!    This will happen when a piecewise polynomial is to be graphed.
!    The first guess for LEFT is therefore taken to be the value
!    returned at the previous call and stored in the local variable ILO.
!
!    A first check ascertains that ILO < LXT.  This is necessary
!    since the present call may have nothing to do with the previous
!    call.  Then, if 
!      XT(ILO) <= X < XT(ILO+1), 
!    we set LEFT = ILO and are done after just three comparisons.
!
!    Otherwise, we repeatedly double the difference ISTEP = IHI - ILO
!    while also moving ILO and IHI in the direction of X, until
!      XT(ILO) <= X < XT(IHI)
!    after which we use bisection to get, in addition, ILO + 1 = IHI.
!    The value LEFT = ILO is then returned.
!
!    Thanks to Daniel Gloger for pointing out an important modification
!    to the routine, so that the piecewise polynomial in B-form is
!    left-continuous at the right endpoint of the basic interval,
!    17 April 2014.
!
!  Modified:
!
!    17 April 2014
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
!    Input, real ( kind = 8 ) XT(LXT), a nondecreasing sequence of values.
!
!    Input, integer ( kind = 4 ) LXT, the dimension of XT.
!
!    Input, real ( kind = 8 ) X, the point whose location with 
!    respect to the sequence XT is to be determined.
!
!    Output, integer ( kind = 4 ) LEFT, the index of the bracketing value:
!      1     if             X  <  XT(1)
!      I     if   XT(I)  <= X  < XT(I+1)
!      I     if   XT(I)  <  X == XT(I+1) == XT(LXT)
!
!    Output, integer ( kind = 4 ) MFLAG, indicates whether X lies within the
!    range of the data.
!    -1:            X  <  XT(1)
!     0: XT(I)   <= X  < XT(I+1)
!    +1: XT(LXT) <  X
!
  implicit none

  integer ( kind = 4 ) lxt

  integer ( kind = 4 ) left
  integer ( kind = 4 ) mflag
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ), save :: ilo = 1
  integer ( kind = 4 ) istep
  integer ( kind = 4 ) middle
  real ( kind = 8 ) x
  real ( kind = 8 ) xt(lxt)

  ihi = ilo + 1

  if ( lxt <= ihi ) then

    if ( xt(lxt) <= x ) then
      go to 110
    end if

    if ( lxt <= 1 ) then
      mflag = -1
      left = 1
      return
    end if

    ilo = lxt - 1
    ihi = lxt

  end if

  if ( xt(ihi) <= x ) then
    go to 20
  end if

  if ( xt(ilo) <= x ) then
    mflag = 0
    left = ilo
    return
  end if
!
!  Now X < XT(ILO).  Decrease ILO to capture X.
!
  istep = 1

10 continue

  ihi = ilo
  ilo = ihi - istep

  if ( 1 < ilo ) then
    if ( xt(ilo) <= x ) then
      go to 50
    end if
    istep = istep * 2
    go to 10
  end if

  ilo = 1

  if ( x < xt(1) ) then
    mflag = -1
    left = 1
    return
  end if

  go to 50
!
!  Now XT(IHI) <= X.  Increase IHI to capture X.
!
20 continue

  istep = 1

30 continue

  ilo = ihi
  ihi = ilo + istep

  if ( ihi < lxt ) then

    if ( x < xt(ihi) ) then
      go to 50
    end if

    istep = istep * 2
    go to 30

  end if

  if ( xt(lxt) <= x ) then
    go to 110
  end if
!
!  Now XT(ILO) < = X < XT(IHI).  Narrow the interval.
!
  ihi = lxt

50 continue

  do

    middle = ( ilo + ihi ) / 2

    if ( middle == ilo ) then
      mflag = 0
      left = ilo
      return
    end if
!
!  It is assumed that MIDDLE = ILO in case IHI = ILO+1.
!
    if ( xt(middle) <= x ) then
      ilo = middle
    else
      ihi = middle
    end if

  end do
!
!  Set output and return.
!
110 continue

  mflag = 1

  if ( x == xt(lxt) ) then
    mflag = 0
  end if

  do left = lxt - 1, 1, -1
    if ( xt(left) < xt(lxt) ) then
      return
    end if
  end do

  return
end subroutine
