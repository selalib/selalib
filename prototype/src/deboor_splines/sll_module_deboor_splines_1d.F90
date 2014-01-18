module sll_module_deboor_splines_1d

#include "sll_memory.h"
#include "sll_working_precision.h"
  implicit none 
  
  
contains
  
   subroutine interv( xt, lxt, x, left, mflag )
    
    !*************************************************************************
    !
    !! INTERV brackets a real value in an ascending vector of values.
    !
    !  Discussion:
    !
    !    The XT array is a set of increasing values.  The goal of the routine
    !    is to determine the largest index I so that XT(I) <= X.
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
    !
    !      XT(ILO) <= X < XT(ILO+1),
    !
    !    we set LEFT = ILO and are done after just three comparisons.
    !
    !    Otherwise, we repeatedly double the difference ISTEP = IHI - ILO
    !    while also moving ILO and IHI in the direction of X, until
    !
    !      XT(ILO) <= X < XT(IHI)
    !
    !    after which we use bisection to get, in addition, ILO + 1 = IHI.
    !    The value LEFT = ILO is then returned.
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
    !    Input, real ( kind = 8 ) XT(LXT), a nondecreasing sequence of values.
    !
    !    Input, integer LXT, the dimension of XT.
    !
    !    Input, real ( kind = 8 ) X, the point whose location with
    !    respect to the sequence XT is to be determined.
    !
    !    Output, integer LEFT, the index of the bracketing value:
    !      1     if             X  <  XT(1)
    !      I     if   XT(I)  <= X  < XT(I+1)
    !      LXT   if  XT(LXT) <= X
    !
    !    Output, integer MFLAG, indicates whether X lies within the
    !    range of the data.
    !    -1:            X  <  XT(1)
    !     0: XT(I)   <= X  < XT(I+1)
    !    +1: XT(LXT) <= X
    !
    implicit none
    
    sll_int32,intent(in):: lxt
    sll_int32,intent(out):: left
    sll_int32,intent(out):: mflag
    sll_int32:: ihi
    sll_int32, save :: ilo = 1
    sll_int32:: istep
    sll_int32:: middle
    sll_real64,intent(in) ::x
    sll_real64,dimension(:),pointer:: xt!(lxt)

    
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
    
10  continue
    
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
20  continue
    
    istep = 1
    
30  continue
    
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
    
50  continue
    
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
    
    do left = lxt, 1, -1
       if ( xt(left) < xt(lxt) ) then
          return
       end if
    end do
    
    return

  end subroutine interv


  subroutine bsplvb ( t, jhigh, index, x, left, biatx )

    !***********************************************************************
    !
    !! BSPLVB evaluates B-splines at a point X with a given knot sequence.
    !
    !  Discusion:
    !
    !    BSPLVB evaluates all possibly nonzero B-splines at X of order
    !
    !      JOUT = MAX ( JHIGH, (J+1)*(INDEX-1) )
    !
    !    with knot sequence T.
    !
    !    The recurrence relation
    !
    !                     X - T(I)               T(I+J+1) - X
    !    B(I,J+1)(X) = ----------- * B(I,J)(X) + --------------- * B(I+1,J)(X)
    !                  T(I+J)-T(I)               T(I+J+1)-T(I+1)
    !
    !    is used to generate B(LEFT-J:LEFT,J+1)(X) from B(LEFT-J+1:LEFT,J)(X)
    !    storing the new values in BIATX over the old.
    !
    !    The facts that
    !
    !      B(I,1)(X) = 1  if  T(I) <= X < T(I+1)
    !
    !    and that
    !
    !      B(I,J)(X) = 0  unless  T(I) <= X < T(I+J)
    !
    !    are used.
    !
    !    The particular organization of the calculations follows
    !    algorithm 8 in chapter X of the text.
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
    !Input, real ( kind = 8 ) T(LEFT+JOUT), the knot sequence.  T is assumed to
    !    be nondecreasing, and also, T(LEFT) must be strictly less than
    !    T(LEFT+1).
    !
    !    Input, integer JHIGH, INDEX, determine the order
    !    JOUT = max ( JHIGH, (J+1)*(INDEX-1) )
    !    of the B-splines whose values at X are to be returned.
    !    INDEX is used to avoid recalculations when several
    !    columns of the triangular array of B-spline values are
    !    needed, for example, in BVALUE or in BSPLVD.
    !    If INDEX = 1, the calculation starts from scratch and the entire
    !    triangular array of B-spline values of orders
    !    1, 2, ...,JHIGH is generated order by order, that is,
    !    column by column.
    !    If INDEX = 2, only the B-spline values of order J+1, J+2, ..., JOUT
    !    are generated, the assumption being that BIATX, J,
    !    DELTAL, DELTAR are, on entry, as they were on exit
    !    at the previous call.  In particular, if JHIGH = 0,
    !    then JOUT = J+1, that is, just the next column of B-spline
    !    values is generated.
    !    Warning: the restriction  JOUT <= JMAX (= 20) is
    !    imposed arbitrarily by the dimension statement for DELTAL
    !    and DELTAR, but is nowhere checked for.
    !
    !    Input, real ( kind = 8 ) X, the point at which the B-splines
    !    are to be evaluated.
    !
    !    Input, integer LEFT, an integer chosen so that
    !    T(LEFT) <= X <= T(LEFT+1).
    !
    !    Output, real ( kind = 8 ) BIATX(JOUT), with BIATX(I) containing the
    !    value at X of the polynomial of order JOUT which agrees
    !    with the B-spline B(LEFT-JOUT+I,JOUT,T) on the interval
    !    (T(LEFT),T(LEFT+1)).
    !
    implicit none
    
    sll_int32, parameter :: jmax = 20
    
    sll_int32:: jhigh
    
    sll_real64,dimension(jhigh):: biatx !(jhigh)
    sll_real64, save, dimension ( jmax ) :: deltal
    sll_real64, save, dimension ( jmax ) :: deltar
    sll_int32:: i
    sll_int32:: index
    sll_int32, save :: j = 1
    sll_int32:: left
    sll_real64:: saved
    sll_real64,dimension(left+jhigh):: t!() left+jhigh
    sll_real64:: term
    sll_real64:: x
    
    if ( index == 1 ) then
       j = 1
       biatx(1) = 1.0_8
       if ( jhigh <= j ) then
          return
       end if
    end if
    
    if ( t(left+1) <= t(left) ) then
       print*,'x=',x
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'BSPLVB - Fatal error!'
       write ( *, '(a)' ) '  It is required that T(LEFT) < T(LEFT+1).'
       write ( *, '(a,i8)' ) '  But LEFT = ', left
       write ( *, '(a,g14.6)' ) '  T(LEFT) =   ', t(left)
       write ( *, '(a,g14.6)' ) '  T(LEFT+1) = ', t(left+1)
       stop
    end if
    
    do
       
       deltar(j) = t(left+j) - x
       deltal(j) = x - t(left+1-j)
       
       saved = 0.0_f64
       do i = 1, j
          term = biatx(i) / ( deltar(i) + deltal(j+1-i) )
          biatx(i) = saved + deltar(i) * term
          saved = deltal(j+1-i) * term
       end do
    
       biatx(j+1) = saved
       j = j + 1
       
       if ( jhigh <= j ) then
          
          exit
       end if
    
    end do
    
    return
  end subroutine bsplvb
  


  subroutine bsplvd ( t, k, x, left, a, dbiatx, nderiv )

    !*************************************************************************
    !
    !! BSPLVD calculates the nonvanishing B-splines and derivatives at X.
    !
    !  Discussion:
    !
    !    Values at X of all the relevant B-splines of order K:K+1-NDERIV
    !    are generated via BSPLVB and stored temporarily in DBIATX.
    !
    !    Then the B-spline coefficients of the required derivatives
    !    of the B-splines of interest are generated by differencing,
    !    each from the preceding one of lower order, and combined with
    !    the values of B-splines of corresponding order in DBIATX
    !    to produce the desired values.
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
    !Input, real ( kind = 8 ) T(LEFT+K), the knot sequence.  It is assumed that
    !    T(LEFT) < T(LEFT+1).  Also, the output is correct only if
    !    T(LEFT) <= X <= T(LEFT+1).
    !
    !    Input, integer K, the order of the B-splines to be evaluated.
    !
    !    Input, real ( kind = 8 ) X, the point at which these values are sought.
    !
    !    Input, integer LEFT, indicates the left endpoint of the interval of
    !    interest.  The K B-splines whose support contains the interval
    !    ( T(LEFT), T(LEFT+1) ) are to be considered.
    !
    !    Workspace, real ( kind = 8 ) A(K,K).
    !
    !    Output, real ( kind = 8 ) DBIATX(K,NDERIV).  DBIATX(I,M) contains
    !    the value of the (M-1)st derivative of the (LEFT-K+I)-th B-spline
    !    of order K for knot sequence T, I=M,...,K, M=1,...,NDERIV.
    !
    !    Input, integer NDERIV, indicates that values of B-splines and their
    !    derivatives up to but not including the NDERIV-th are asked for.
    !
    implicit none
    
    sll_int32 :: k
    sll_int32 :: left
    sll_int32 :: nderiv
    
    sll_real64, dimension(k,k):: a!(k,k)
    sll_real64,dimension(k,nderiv), intent(out) :: dbiatx!(k,nderiv)
    sll_real64:: factor
    sll_real64:: fkp1mm
    sll_int32 :: i
    sll_int32 :: ideriv
    sll_int32 :: il
    sll_int32 :: j
    sll_int32 :: jlow
    sll_int32 :: jp1mid
    sll_int32 :: ldummy
    sll_int32 :: m
    sll_int32 :: mhigh
    !  sll_real64 sum1  ! this one is not used...
    sll_real64,dimension(left+k):: t ! (left+k)
    sll_real64:: x
    
    
    mhigh = max ( min ( nderiv, k ), 1 )
    !
    !  MHIGH is usually equal to NDERIV.
    !
    call bsplvb ( t, k+1-mhigh, 1, x, left, dbiatx )
    
    if ( mhigh == 1 ) then
       return
    end if
    !
    !  The first column of DBIATX always contains the B-spline values
    !  for the current order.  These are stored in column K+1-current
    !  order before BSPLVB is called to put values for the next
    !  higher order on top of it.
    !
    ideriv = mhigh
    do m = 2, mhigh
       jp1mid = 1
       do j = ideriv, k
          dbiatx(j,ideriv) = dbiatx(jp1mid,1)
          jp1mid = jp1mid + 1
          
       end do
       ideriv = ideriv - 1
       
       call bsplvb ( t, k+1-ideriv, 2, x, left, dbiatx )
       
    end do
    !
    !  At this point, B(LEFT-K+I, K+1-J)(X) is in DBIATX(I,J) for
    !  I=J,...,K and J=1,...,MHIGH ('=' NDERIV).
    !
    !  In particular, the first column of DBIATX is already in final form.
    !
    !  To obtain corresponding derivatives of B-splines in subsequent columns,
    !  generate their B-representation by differencing, then evaluate at X.
    !
    jlow = 1
    do i = 1, k
       a(jlow:k,i) = 0.0D+00
       jlow = i
       a(i,i) = 1.0D+00
    end do
    !
    !  At this point, A(.,J) contains the B-coefficients for the J-th of the
    !  K B-splines of interest here.
    !
    do m = 2, mhigh
       
       fkp1mm = real ( k + 1 - m, kind = 8 )
       il = left
       i = k
       !
       !  For J = 1,...,K, construct B-coefficients of (M-1)st derivative of
       !  B-splines from those for preceding derivative by differencing
       !  and store again in  A(.,J).  The fact that  A(I,J) = 0 for
       !  I < J is used.
       !
       do ldummy = 1, k+1-m
          
          factor = fkp1mm / ( t(il+k+1-m) - t(il) )
          !
          !  The assumption that T(LEFT) < T(LEFT+1) makes denominator
          !  in FACTOR nonzero.
          !
          a(i,1:i) = ( a(i,1:i) - a(i-1,1:i) ) * factor
          
          il = il - 1
          i = i - 1
       
       end do
       !
       !  For I = 1,...,K, combine B-coefficients A(.,I) with B-spline values
       !  stored in DBIATX(.,M) to get value of (M-1)st derivative of
       !  I-th B-spline (of interest here) at X, and store in DBIATX(I,M).
       !
       !  Storage of this value over the value of a B-spline
       !  of order M there is safe since the remaining B-spline derivatives
       !  of the same order do not use this value due to the fact
       !  that  A(J,I) = 0  for J < I.
       !
       do i = 1, k
          
          jlow = max ( i, m )
          
          dbiatx(i,m) = dot_product ( a(jlow:k,i), dbiatx(jlow:k,m) )
          
       end do
    
    end do
    return
  end subroutine bsplvd

 
  
  
  function bvalue( t, bcoef, n, k, x, jderiv ) result(res)
    
    !*********************************************************************
    !
    !! BVALUE evaluates a derivative of a spline from
    ! its B-spline representation.
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
    
    sll_int32 :: k
    sll_int32 :: n
    
    sll_real64,dimension(:),pointer:: aj !(k)
    sll_real64,dimension(:),pointer:: bcoef!(n)
    sll_real64:: res
    !sll_real64:: tmp_value
    sll_real64,dimension(:),pointer:: dl!(k)
    sll_real64,dimension(:),pointer:: dr!(k)
    sll_int32 :: i
    sll_int32 :: ilo
    sll_int32 :: j
    sll_int32 :: jc
    sll_int32 :: jcmax
    sll_int32 :: jcmin
    sll_int32 :: jderiv
    sll_int32 :: jj
    sll_int32 :: mflag
    sll_real64,dimension(:),pointer:: t!(n+k)
    sll_real64:: x
    sll_int32 :: ierr
    
    res = 0.0_8
    
    SLL_ALLOCATE(aj(k),ierr)
    SLL_ALLOCATE(dl(k),ierr)
    SLL_ALLOCATE(dr(k),ierr)
    
    aj(:)=0.0_8
    dl(:)=0.0_8
    dr(:)=0.0_8
    
    if ( k <= jderiv ) then
       SLL_DEALLOCATE(aj,ierr)
       SLL_DEALLOCATE(dl,ierr)
       SLL_DEALLOCATE(dr,ierr)
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
       SLL_DEALLOCATE(aj,ierr)
       SLL_DEALLOCATE(dl,ierr)
       SLL_DEALLOCATE(dr,ierr)
       return
    end if
    !
    !  If K = 1 (and JDERIV = 0), BVALUE = BCOEF(I).
    !
    if ( k <= 1 ) then
       res = bcoef(i)
       SLL_DEALLOCATE(aj,ierr)
       SLL_DEALLOCATE(dl,ierr)
       SLL_DEALLOCATE(dr,ierr)
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

    res = aj(1)
    
    SLL_DEALLOCATE(aj,ierr)
    SLL_DEALLOCATE(dl,ierr)
    SLL_DEALLOCATE(dr,ierr)
    
    return
    
    
    
  end function bvalue
  

  function dvalue1d(ar_x,&
       ai_nx,&
       ai_kx,&
       apr_Bcoef,&
       apr_tx,&
       deriv1) result(res)
    implicit none
    ! INPUT
    sll_real64 :: ar_x
    sll_real64 :: res
    sll_int32  :: ai_nx, ai_kx
    sll_int32  :: deriv1
    sll_real64, dimension ( : ),pointer :: apr_tx ! ai_nx + ai_kx
    sll_real64, dimension ( :),pointer :: apr_Bcoef ! ai_nx 			
    
    
    
    
    res = bvalue(apr_tx,&
         apr_Bcoef,&
         ai_nx,&
         ai_kx,&
         ar_x,&
         deriv1 )
    
  end function dvalue1d
  

  subroutine spli1d_dir ( &
     ai_nx,&
     ai_kx,&
     apr_taux,&
     apr_g,&
     apr_Bcoef,&
     apr_tx )
    implicit none
    ! INPUT
    sll_int32  :: ai_nx, ai_kx
    sll_real64, dimension ( ai_nx) :: apr_taux
    sll_real64, dimension ( ai_nx) :: apr_g
    ! OUTPUT
    sll_real64, dimension ( ai_nx  ) :: apr_Bcoef
    sll_real64, dimension ( ai_nx + ai_kx ) :: apr_tx
    ! LOCAL VARIABLES		
    sll_real64, dimension ( ai_nx ) :: lpr_work1
    !sll_real64, dimension ( ai_nx         ) :: lpr_work2
    sll_real64, dimension ( ai_nx *( 2*ai_kx-1) ) :: lpr_work31
    !sll_real64, dimension ( (ai_nx-ai_kx)*(2*ai_kx+3)+5*ai_kx+3 ) :: scrtch
    !sll_real64, dimension ( ai_nx + ai_kx ) :: t 
    sll_int32  :: li_i, li_iflag
    
    lpr_work1(:) = 0.0
    
    ! *** set up knots
  !     interpolate between knots
    ! x
    !  if (ai_kx <= 2) then 
    apr_tx ( 1 : ai_kx ) = apr_taux ( 1 )
    apr_tx ( ai_nx + 1 : ai_nx + ai_kx ) = apr_taux ( ai_nx )
  
    if ( mod(ai_kx,2) == 0 ) then
       do li_i = ai_kx + 1, ai_nx
          apr_tx ( li_i ) = apr_taux ( li_i - ai_kx/2 ) 
          
       end do
    else
       
       do li_i = ai_kx + 1, ai_nx
          apr_tx ( li_i ) = &
               0.5*( apr_taux ( li_i - (ai_kx-1)/2 ) + &
               apr_taux ( li_i -1 - (ai_kx-1)/2 ) )
          
       end do
       
    end if
    
    apr_Bcoef = 0.0_8
    do li_i = 1, ai_nx
       apr_Bcoef ( li_i ) = apr_g ( li_i )
    end do
    
    !  *** construct b-coefficients of interpolant
    !
    call splint ( &
         apr_taux,&
         apr_g,&
         apr_tx, &
         ai_nx,&
         ai_kx, &
         lpr_work31,&
         apr_Bcoef, &
         li_iflag )
  
  end subroutine spli1d_dir


  subroutine spli1d_per(&
       ar_L,&
       ai_nx,&
       ai_kx,&
       apr_taux,&
       apr_g,&
       apr_Bcoef,&
       apr_tx)
    ! CALLED WHEN WE WANT TO INTERPOL WITH A PERIODIC 
    ! FIRST PARAM WITH A PERIOD = ar_L
    implicit none
    ! INPUT
    sll_real64 :: ar_L 
    sll_int32  :: ai_nx, ai_kx
    sll_real64, dimension ( ai_nx) :: apr_taux
    sll_real64, dimension ( ai_nx) :: apr_g
    ! OUTPUT
    sll_real64, dimension ( ai_nx ) :: apr_Bcoef
    sll_real64, dimension ( ai_nx + ai_kx) :: apr_tx
    ! LOCAL VARIABLES		
    sll_real64, dimension ( ai_nx) :: lpr_taux
    sll_real64, dimension ( ai_nx) :: lpr_g
    sll_real64, dimension ( ai_nx *( 2*ai_kx-1) ) :: lpr_work
    !sll_real64, dimension ( (ai_nx-ai_kx)*(2*ai_kx+3)+5*ai_kx+3 ) :: scrtch
    sll_int32 :: iflag
    sll_int32 :: li_i
    if ( ar_L == 0.0_8 ) then
       print*,'Error spli1d_per : called with a period = 0 '
       stop
    end if
  
    
    lpr_taux ( 1 : ai_nx - 1 ) = apr_taux ( 1 : ai_nx - 1 )
    lpr_taux ( ai_nx ) = apr_taux ( ai_nx ) !apr_taux ( 1 ) + ar_L
    
    lpr_g ( 1 : ai_nx - 1  ) = apr_g ( 1 : ai_nx - 1 )
    lpr_g ( ai_nx) = apr_g ( ai_nx ) !apr_g ( 1 )		
    
    
    apr_tx ( 1 : ai_kx ) = lpr_taux ( 1 )
    apr_tx ( ai_nx + 1 : ai_nx + ai_kx ) = lpr_taux ( ai_nx )
  
    
    if ( mod(ai_kx,2) == 0 ) then
       do li_i = ai_kx + 1, ai_nx
          apr_tx ( li_i ) = lpr_taux ( li_i - ai_kx/2 ) 
          
       end do
    else
       
       do li_i = ai_kx + 1, ai_nx
          apr_tx ( li_i ) = &
               0.5*( lpr_taux ( li_i - (ai_kx-1)/2 ) + &
               lpr_taux ( li_i -1 - (ai_kx-1)/2 ) )
          
       end do
       
    end if

    
    call splint ( &
         lpr_taux,&
         lpr_g,&
         apr_tx,&
         ai_nx, &
         ai_kx, &
         lpr_work,&
         apr_Bcoef,&
         iflag)
    
  end subroutine spli1d_per
  


  subroutine splint ( tau, gtau, t, n, k, q, bcoef, iflag )
    
    !*************************************************************************
    !
    !! SPLINT produces the B-spline coefficients BCOEF of an 
    ! interpolating spline.
    !
    !  Discussion:
    !
    !    The spline is of order K with knots T(1:N+K), and takes on the 
    !    value GTAU(I) at TAU(I), for I = 1 to N.
    !
    !    The I-th equation of the linear system 
    !
    !      A * BCOEF = B 
    !
    !    for the B-spline coefficients of the interpolant enforces interpolation
    !    at TAU(1:N).
    !
    !    Hence, B(I) = GTAU(I), for all I, and A is a band matrix with 2*K-1
    !    bands, if it is invertible.
    !
    !    The matrix A is generated row by row and stored, diagonal by diagonal,
    !    in the rows of the array Q, with the main diagonal going
    !    into row K.  See comments in the program.
    !
    !    The banded system is then solved by a call to BANFAC, which 
    !    constructs the triangular factorization for A and stores it again in
    !    Q, followed by a call to BANSLV, which then obtains the solution
    !    BCOEF by substitution.
    !
    !    BANFAC does no pivoting, since the total positivity of the matrix
    !    A makes this unnecessary.
    !
    !    The linear system to be solved is (theoretically) invertible if
    !    and only if
    !      T(I) < TAU(I) < TAU(I+K), for all I.
    !    Violation of this condition is certain to lead to IFLAG = 2.
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
    !    ISBN: 0387953663,
    !    LC: QA1.A647.v27.
    !
    !  Parameters:
    !
    ! Input, real ( kind = 8 ) TAU(N), the data point abscissas.The entries in
    !    TAU should be strictly increasing.
    !
    !    Input, real ( kind = 8 ) GTAU(N), the data ordinates.
    !
    !    Input, real ( kind = 8 ) T(N+K), the knot sequence.
    !
    !    Input, integer ( kind = 4 ) N, the number of data points.
    !
    !    Input, integer ( kind = 4 ) K, the order of the spline.
    !
    !    Output, real ( kind = 8 ) Q((2*K-1)*N), the triangular factorization
    !    of the coefficient matrix of the linear system for the B-coefficients 
    !    of the spline interpolant.  The B-coefficients for the interpolant 
    !    of an additional data set can be obtained without going through all 
    !    the calculations in this routine, simply by loading HTAU into BCOEF 
    !    and then executing the call:
    !      call banslv ( q, 2*k-1, n, k-1, k-1, bcoef )
    !
    !    Output, real ( kind = 8 ) BCOEF(N), the B-spline coefficients of 
    !    the interpolant.
    !
    !    Output, integer ( kind = 4 ) IFLAG, error flag.
    !    1, = success.
    !    2, = failure.
    !
    implicit none
    
    sll_int32 ::  n
    
    sll_real64,dimension(n):: bcoef ! (n)
    sll_real64,dimension(n)::  gtau ! (n)
    sll_int32 :: i
    sll_int32 :: iflag
    sll_int32 :: ilp1mx
    sll_int32 :: j
    sll_int32 :: jj
    sll_int32 :: k
    sll_int32 :: kpkm2
    sll_int32 :: left
    sll_real64, dimension((2*k-1)*n) :: q!((2*k-1)*n)
    sll_real64,dimension(n+k) ::  t!(n+k)
    sll_real64,dimension(n) ::  tau!!(n)
    sll_real64:: taui
  
    kpkm2 = 2 * ( k - 1 )
    left = k
    q(1:(2*k-1)*n) = 0.0_f64
    !
    !  Loop over I to construct the N interpolation equations.
    !
    do i = 1, n
       
       taui = tau(i)
       ilp1mx = min ( i + k, n + 1 )
       !
       !  Find LEFT in the closed interval (I,I+K-1) such that
       !
       !    T(LEFT) <= TAU(I) < T(LEFT+1)
       !
       !  The matrix is singular if this is not possible.
       !
       left = max ( left, i )
       
       if ( taui < t(left) ) then
          iflag = 2
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'SPLINT - Fatal Error!'
          write ( *, '(a)' ) '  The linear system is not invertible!'
          return
       end if
   
       do while ( t(left+1) <= taui )
          
          left = left + 1
          
          if ( left < ilp1mx ) then
             cycle
          end if
          
          left = left - 1
          
          if ( t(left+1) < taui ) then
             iflag = 2
             write ( *, '(a)' ) ' '
             write ( *, '(a)' ) 'SPLINT - Fatal Error!'
             write ( *, '(a)' ) '  The linear system is not invertible!'
             return
          end if
          
          exit
          
       end do
       !
       !  The I-th equation enforces interpolation at TAUI, hence for all J,
       !
       !    A(I,J) = B(J,K,T)(TAUI).
       !
       !Only the K entries with J = LEFT-K+1,...,LEFT actually might be nonzero.
       !
       !These K numbers are returned, in BCOEF 
       ! (used for temporary storage here),
       !  by the following.
       !
       call bsplvb ( t, k, 1, taui, left, bcoef )
       !
       !  We therefore want BCOEF(J) = B(LEFT-K+J)(TAUI) to go into
       !  A(I,LEFT-K+J), that is, into Q(I-(LEFT+J)+2*K,(LEFT+J)-K) since
       !  A(I+J,J) is to go into Q(I+K,J), for all I, J, if we consider Q
       !  as a two-dimensional array, with  2*K-1 rows.  See comments in
       !  BANFAC.
       !
       !  In the present program, we treat Q as an equivalent
       !  one-dimensional array, because of fortran restrictions on
       !  dimension statements.
       !
       !  We therefore want  BCOEF(J) to go into the entry of Q with index:
       !
       !    I -(LEFT+J)+2*K + ((LEFT+J)-K-1)*(2*K-1)
       !   = I-LEFT+1+(LEFT -K)*(2*K-1) + (2*K-2)*J
       !
       jj = i - left + 1 + ( left - k ) * ( k + k - 1 )
       
       do j = 1, k
          jj = jj + kpkm2
          q(jj) = bcoef(j)
       end do
       
    end do
    !
    !  Obtain factorization of A, stored again in Q.
    !
    call banfac ( q, k+k-1, n, k-1, k-1, iflag )
    
    if ( iflag == 2 ) then
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'SPLINT - Fatal Error!'
       write ( *, '(a)' ) '  The linear system is not invertible!'
       return
    end if
    !
    !  Solve 
    !
    !    A * BCOEF = GTAU
    !
    !  by back substitution.
    !
    bcoef(1:n) = gtau(1:n)
    
    call banslv ( q, k+k-1, n, k-1, k-1, bcoef )
    
    return
  end subroutine splint

  
  
end module sll_module_deboor_splines_1d
