module sll_module_deboor_splines_2d

  

  implicit none 

  
contains
  
  
  

  subroutine interv( xt, lxt, x, left, mflag )
  
    !**********************************************************************
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
    sLL_real64,intent(in) ::x
    sLL_real64,dimension(:),pointer:: xt!(lxt)
  
    print*,lxt,x
    print*, xt(1)
    
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
 
    do left = lxt, 1, -1
       if ( xt(left) < xt(lxt) ) then
          return
       end if
    end do
    
    return
    print*, 'FFFFF'
end subroutine interv

  
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
    
    sll_int32 :: k
    sll_int32 :: n
    
    sll_real64,dimension(:),pointer:: aj !(k)
    sll_real64,dimension(:),pointer:: bcoef!(n)
    sll_real64:: bvalue
    sll_real64:: tmp_value
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
    
    bvalue = 0.0_8
    
    SLL_ALLOCATE(aj(k),ierr)
    SLL_ALLOCATE(dl(k),ierr)
    SLL_ALLOCATE(dr(k),ierr)
    
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
    SLL_DEALLOCATE(aj,ierr)
    SLL_DEALLOCATE(dl,ierr)
    SLL_DEALLOCATE(dr,ierr)
    
    
  end function bvalue

  

  
  
  ! Just for the record: what is ar_x, ar_y, ai_nx, etc., etc.? What is this
  ! function supposed to do? This is a deeply frustrating file.
  subroutine bvalue2d(&
       ar_x,&
       ar_y,&
       ai_nx,&
       ai_kx,&
       ai_ny,&
       ai_ky,&
       apr_Bcoef,&
       apr_tx,&
       sz_tx,&
       apr_ty,&
       sz_ty,&
       val )
    implicit none
    ! INPUT
    sll_real64 :: ar_x, ar_y
    sll_real64 :: bvalue,val
    sll_int32  :: ai_nx, ai_kx, ai_ny, ai_ky,sz_tx,sz_ty
    sll_real64, dimension(:), pointer :: apr_tx !  ai_nx + ai_kx 
    sll_real64, dimension(:), pointer :: apr_ty !  ai_ny + ai_ky	
    sll_real64, dimension(:,:),pointer :: apr_Bcoef!( ai_nx,ai_ny)			
    ! LOCAL VARIABLES
    sll_int32  :: li_i, li_j, li_mflag, li_lefty
    sll_real64, dimension(:),pointer :: lpr_coef ! ai_ky
    sll_real64, dimension(:),pointer :: tmp_tab
    sll_real64, dimension(:),pointer :: tmp_ty
    sll_int32 :: ierr
    
    
    SLL_ALLOCATE(lpr_coef(ai_ky),ierr)
    SLL_ALLOCATE(tmp_tab(ai_nx),ierr)
    SLL_ALLOCATE(tmp_ty( 2*ai_ky ),ierr)
    
    call interv ( apr_ty, sz_ty, ar_y, li_lefty, li_mflag )
    
    if ( li_mflag .NE. 0 ) then
       val = 0.0_8
       return 
    end if
    
    do li_j = 1, ai_ky
       
       
       tmp_tab = apr_bcoef ( 1:ai_nx , li_lefty - ai_ky + li_j )
       
       lpr_coef ( li_j ) = bvalue(&
            apr_tx,&
            tmp_tab,&
            ai_nx,&
            ai_kx,&
            ar_x,&
            0 )
       
       
    end do
    
    tmp_ty =  apr_ty ( li_lefty - ai_ky + 1 : li_lefty + ai_ky)
    val = bvalue(&
         tmp_ty,&
         lpr_coef,&
         ai_ky,&
         ai_ky,&
         ar_y,&
         0 )
    
    SLL_DEALLOCATE(lpr_coef,ierr)
    SLL_DEALLOCATE(tmp_tab,ierr)
    SLL_DEALLOCATE(tmp_ty,ierr)
  end subroutine bvalue2d


  function dvalue2d(&
       ar_x,&
       ar_y,&
       ai_nx,&
       ai_kx,&
       ai_ny,&
       ai_ky,&
       apr_Bcoef,&
       apr_tx,&
       apr_ty,deriv1,deriv2 ) result(res)
    implicit none
    ! INPUT
    sll_real64 :: ar_x, ar_y
    sll_real64 :: bvalue,res
    sll_int32  :: ai_nx, ai_kx, ai_ny, ai_ky
    sll_int32  :: deriv1,deriv2
    sll_real64, dimension ( : ),pointer :: apr_tx ! ai_nx + ai_kx
    sll_real64, dimension ( : ),pointer :: apr_ty ! ai_ny + ai_ky
    sll_real64, dimension ( : , : ),pointer :: apr_Bcoef !(ai_nx,ai_ny)
    ! LOCAL VARIABLES
    sll_int32  :: li_i, li_j, li_mflag, li_lefty
    sll_real64, dimension ( :) :: lpr_coef ! ai_ky			
    

    SLL_ALLOCATE(lpr_coef(ai_ky),ierr)
    SLL_ALLOCATE(tmp_coef(ai_nx),ierr)
    SLL_ALLOCATE(tmp_ty(2*ai_ky),ierr)
    
    call interv ( apr_ty, ai_ny + ai_ky, ar_y, li_lefty, li_mflag )
    
    if ( li_mflag .NE. 0 ) then
       res = 0.0_8
       return 
    end if
    
    do li_j = 1, ai_ky
       
       tmp_coef = apr_bcoef ( 1:ai_nx , li_lefty - ai_ky + li_j )
       lpr_coef ( li_j ) = bvalue(&
            apr_tx,&
            tmp_coef,&
            ai_nx,&
            ai_kx,&
            ar_x,&
            deriv1 )
       
    end do
    tmp_ty =  apr_ty ( li_lefty - ai_ky + 1 : li_lefty + ai_ky)
    res = bvalue(&
         tmp_ty,&
         lpr_coef,&
         ai_ky,&
         ai_ky,&
         ar_y,&
         deriv2 )
    
  end function dvalue2d
  


  subroutine spli2d ( tau, gtau, t, n, k, m, work, q, bcoef, iflag )

    !*****************************************************************************80
    !
    !! SPLI2D produces a interpolatory tensor product spline.
    !
    !  Discussion:
    !
    !    SPLI2D is an extended version of SPLINT.
    !
    !    SPLI2D produces the B-spline coefficients BCOEF(J,.) of the
    !    spline of order K with knots T(1:N+K), which takes on
    !    the value GTAU(I,J) at TAU(I), I=1,..., N, J=1,...,M.
    !
    !    The I-th equation of the linear system
    !
    !      A * BCOEF = B
    !
    !    for the B-spline coefficients of the interpolant enforces
    !    interpolation at TAU(I), I=1,...,N.  Hence,  B(I) = GTAU(I),
    !    for all I, and A is a band matrix with 2*K-1 bands, if it is
    !    invertible.
    !
    !    The matrix A is generated row by row and stored, diagonal by
    !    diagonal, in the rows of the array Q, with the main diagonal
    !    going into row K.
    !
    !    The banded system is then solved by a call to BANFAC, which
    !    constructs the triangular factorization for A and stores it
    !    again in Q, followed by a call to BANSLV, which then obtains
    !    the solution BCOEF by substitution.
    !
    !     The linear system to be solved is theoretically invertible if
    !     and only if
    !
    !       T(I) < TAU(I) < TAU(I+K), for all I.
    !
    !     Violation of this condition is certain to lead to IFLAG = 2.
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
    !    Input, real ( kind = 8 ) TAU(N), contains the data point abscissas.
    !    TAU must be strictly increasing
    !
    !    Input, real ( kind = 8 ) GTAU(N,M), contains the data point ordinates.
    !
    !    Input, real ( kind = 8 ) T(N+K), the knot sequence.
    !
    !    Input, integer N, the number of data points and the
    !    dimension of the spline space SPLINE(K,T)
    !
    !    Input, integer K, the order of the spline.
    !
    !    Input, integer M, the number of data sets.
    !
    !    Work space, real ( kind = 8 ) WORK(N).
    !
    !    Output, real ( kind = 8 ) Q(2*K-1)*N, the triangular
    !    factorization of the coefficient matrix of the linear
    !    system for the B-spline coefficients of the spline interpolant.
    !    The B-spline coefficients for the interpolant of an additional
    !    data set ( TAU(I), HTAU(I) ), I=1,...,N  with the same data
    !    abscissae can be obtained without going through all the
    !    calculations in this routine, simply by loading HTAU into
    !    BCOEF and then using the statement
    !      CALL BANSLV ( Q, 2*K-1, N, K-1, K-1, BCOEF )
    !
    !    Output, real ( kind = 8 ) BCOEF(N), the B-spline coefficients of
    !    the interpolant.
    !
    !    Output, integer IFLAG, error indicator.
    !    1, no error.
    !    2, an error occurred, which may have been caused by
    !       singularity of the linear system.
    !
    implicit none
    
    sll_int32 :: m
    sll_int32 :: n
    
    sll_real64,dimension(:,:),pointer:: bcoef !(m,n)
    sll_real64,dimension(:,:),pointer:: gtau  !(n,m)
    sll_int32 :: i
    sll_int32 :: iflag
    sll_int32 :: ilp1mx
    sll_int32 :: j
    sll_int32 :: jj
    sll_int32 :: k
    sll_int32 :: left
    sll_real64,dimension((2*k-1)*n):: q!((2*k-1)*n)
    sll_real64,dimension(:,:),pointer:: t!(n+k)
    sll_real64,dimension(:,:),pointer:: tau!(n)
    sll_real64:: taui
    sll_real64,dimension(n):: work!(n)
    
    left = k
    
    !print*, t
    q(1:(2*k-1)*n) = 0.0_f64
    !
    !  Construct the N interpolation equations.
    !
    
    do i = 1, n
       
       taui = tau(i)
       ilp1mx = min ( i + k, n + 1 )
       !
       !  Find the index LEFT in the closed interval (I,I+K-1) such that:
       !
       !    T(LEFT) < = TAU(I) < T(LEFT+1)
       !
       !  The matrix will be singular if this is not possible.
       !
       left = max ( left, i )
       
       if ( taui < t(left) ) then
          iflag = 2
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'SPLI2D - Fatal error!'
          write ( *, '(a)' ) '  The TAU array is not strictly increasing .'
          !print*, taui, t(left),left
          stop
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
             write ( *, '(a)' ) 'SPLI2D - Fatal error!'
             write ( *, '(a)' ) '  The TAU array is not strictly increasing.'
             !print*, taui, t(left+1),left
             stop
          end if
          
          exit

       end do
       !
       !  The I-th equation enforces interpolation at TAUI, hence
       !
       !    A(I,J) = B(J,K,T)(TAUI), for all J.
       !
       !  Only the K entries with J = LEFT-K+1, ..., LEFT actually might be
       !  nonzero.  These K numbers are returned, in WORK (used for
       !  temporary storage here), by the following call:
       !
       call bsplvb ( t, k, 1, taui, left, work )
       !print*, 'achtung',taui
       ! print*, 'work', work(1:k)
       !
       !  We therefore want
       !
       !    WORK(J) = B(LEFT-K+J)(TAUI)
       !
       !  to go into
       !
       !    A(I,LEFT-K+J),
       !
       !  that is, into  Q(I-(LEFT+J)+2*K,(LEFT+J)-K) since
       !  A(I+J,J) is to go into Q(I+K,J), for all I, J, if we consider Q
       !  as a two-dimensional array, with  2*K-1 rows.  See comments in
       !  BANFAC.
       !
       !  In the present program, we treat Q as an equivalent one-dimensional
       !  array, because of fortran restrictions on dimension statements.
       !
       !  We therefore want WORK(J) to go into the entry of Q with index:
       !    I -(LEFT+J)+2*K + ((LEFT+J)-K-1)*(2*K-1)
       !    = I-LEFT+1+(LEFT -K)*(2*K-1) + (2*K-2)*J
       !
       jj = i - left + 1 + ( left - k ) * ( k + k - 1 )
       
       do j = 1, k
          jj = jj + k + k - 2
          q(jj) = work(j)
       end do
       
    end do
    
    !
    !  Factor A, stored again in Q.
    !
    call banfac ( q, k+k-1, n, k-1, k-1, iflag )
    
    if ( iflag == 2 ) then
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'SPLI2D - Fatal error!'
       write ( *, '(a)' ) '  BANFAC reports that the matrix is singular.'
       stop
    end if
    !  Solve
    !
    !    A * BCOEF = GTAU
    !
    !  by back substitution.
    
    do j = 1, m
       
       work(1:n) = gtau(1:n,j)
       
       
       call banslv ( q, k+k-1, n, k-1, k-1, work )
       
       bcoef(j,1:n) = work(1:n)
       
    end do
   
    
    return
  end subroutine spli2d

  subroutine spli2d_custom ( &
     ai_nx,&
     ai_kx,&
     apr_taux,&
     ai_ny,&
     ai_ky,&
     apr_tauy,&
     apr_g,&
     apr_Bcoef,&
     apr_tx,&
     apr_ty )
    implicit none
    ! INPUT
    sll_int32  :: ai_nx, ai_kx, ai_ny, ai_ky
    sll_real64, dimension (:),pointer :: apr_taux !!ai_nx
    sll_real64, dimension (:),pointer :: apr_tauy	!! ai_ny	
    sll_real64, dimension (:,:),pointer :: apr_g    ! ai_nx,ai_ny	
    ! OUTPUT
    sll_real64, dimension (:,:),pointer :: apr_Bcoef !ai_nx , ai_ny 
    sll_real64, dimension ( : ),pointer :: apr_tx ! ai_nx + ai_kx
    sll_real64, dimension ( : ),pointer :: apr_ty ! ai_ny + ai_ky 
    ! LOCAL VARIABLES		
    sll_real64, dimension ( ai_nx , ai_ny ) :: lpr_work1
    sll_real64, dimension ( ai_nx         ) :: lpr_work2
    sll_real64, dimension ( ai_nx * ai_ny ) :: lpr_work3
    sll_real64, dimension ( ai_nx *( 2*ai_kx-1) ) :: lpr_work31
    sll_real64, dimension (( 2*ai_ky-1) * ai_ny ) :: lpr_work32
    sll_real64, dimension ( ai_ny         ) :: lpr_work4
    sll_real64, dimension ( ai_ny , ai_nx ) :: lpr_work5
    sll_real64, dimension ( (ai_nx-ai_kx)*(2*ai_kx+3)+5*ai_kx+3 ) :: scrtch
    sll_real64, dimension ( (ai_ny-ai_ky)*(2*ai_ky+3)+5*ai_ky+3 ) :: scrtch1
    sll_real64, dimension ( ai_nx + ai_kx ) :: t 
    sll_real64, dimension ( ai_ny ) :: apr_ty_bis
    sll_int32  :: li_i, li_j, li_iflag,iflag,iflag1
    sll_int32 :: o
    
    
    lpr_work1(:,:) = 0.0
    
    ! *** set up knots
    !     interpolate between knots
    
    apr_tx ( 1 : ai_kx ) = apr_taux ( 1 )
    apr_tx ( ai_nx + 1 : ai_nx + ai_kx ) = apr_taux ( ai_nx )
    
    if ( mod(ai_kx,2) == 0 ) then
       do li_i = ai_kx + 1, ai_nx
          apr_tx ( li_i ) = apr_taux ( li_i - ai_kx/2 ) 
          
       end do
     else
        
        do li_i = ai_kx + 1, ai_nx
           apr_tx ( li_i ) = &
                0.5*( apr_taux ( li_i - (ai_kx-1)/2 ) + apr_taux ( li_i -1 - (ai_kx-1)/2 ) )
           
        end do
     
     end if
     apr_Bcoef = 0.0_8
     do li_i = 1, ai_nx
        do li_j = 1, ai_ny
           apr_Bcoef ( li_i, li_j ) = apr_g ( li_i, li_j )
        end do
     end do
     !  *** construct b-coefficients of interpolant
     !
     apr_ty = 0.0_8		
     
     if ( mod(ai_ky,2) == 0 ) then
        do li_i = ai_ky + 1, ai_ny
           apr_ty ( li_i ) = apr_tauy ( li_i - ai_ky/2 ) 
           
        end do
     else
        
        do li_i = ai_ky + 1, ai_ny
           apr_ty ( li_i ) = &
                0.5*( apr_tauy ( li_i - (ai_ky-1)/2 ) + apr_tauy ( li_i -1 - (ai_ky-1)/2 ) )
           
        end do
        
     end if
     apr_ty ( 1 : ai_ky ) = apr_tauy ( 1 )
     apr_ty ( ai_ny + 1 : ai_ny + ai_ky ) = apr_tauy ( ai_ny )
     
     apr_ty_bis = apr_tauy(1:ai_ny)
     
     call spli2d ( &
          apr_taux,&
          apr_Bcoef,&
          apr_tx, &
          ai_nx,&
          ai_kx, &
          ai_ny, &
          lpr_work2,&
          lpr_work31,&
          lpr_work5, &
          li_iflag )
    
     apr_bcoef(:,:) =0.0_8
     lpr_work4 = 0.0_8
     lpr_work3 = 0.0_8
     lpr_work32= 0.0_8
     
     call spli2d ( &
          apr_ty_bis,&
          lpr_work5,&
          apr_ty,&
          ai_ny, &
          ai_ky, &
          ai_nx, &
          lpr_work4, &
          lpr_work32,&
          apr_bcoef, &
          li_iflag )
     
  
   end subroutine spli2d_custom
   

   
   subroutine spli2d_perdir (&
        ar_L,&
        ai_nx,&
        ai_kx,&
        apr_taux,&
        ai_ny,&
        ai_ky,&
        apr_tauy,&
        apr_g,&
        apr_Bcoef,&
        apr_tx,&
        apr_ty )
     ! CALLED WHEN WE WANT TO INTERPOL WITH A PERIODIC FIRST PARAM WITH A PERIOD = ar_L
     implicit none
     ! INPUT
     sll_real64 :: ar_L 
     sll_int32  :: ai_nx, ai_kx, ai_ny, ai_ky
     sll_real64, dimension ( :),pointer :: apr_taux ! ai_nx- 1
     sll_real64, dimension ( ai_ny	),pointer :: apr_tauy ! ai_ny		
     sll_real64, dimension ( :,:) :: apr_g	!ai_nx - 1, ai_ny
     ! OUTPUT
     sll_real64, dimension ( ai_nx , ai_ny	) :: apr_Bcoef
     sll_real64, dimension ( ai_nx + ai_kx	) :: apr_tx
     sll_real64, dimension ( ai_ny + ai_ky ) :: apr_ty
     ! LOCAL VARIABLES		
     sll_real64, dimension ( ai_nx	) :: lpr_taux		
     sll_real64, dimension ( ai_nx ,ai_ny) :: lpr_g	
     
     if ( ar_L == 0.0_8 ) then
        print*,'Error spli2d_per : called with a period = 0 '
        stop
     end if
     
    
     lpr_taux ( 1 : ai_nx - 1 ) = apr_taux ( 1 : ai_nx - 1 )				
     lpr_taux ( ai_nx ) = apr_taux ( 1 ) + ar_L						
     
     lpr_g ( 1 : ai_nx - 1 , 1 : ai_ny ) = apr_g ( 1 : ai_nx - 1 , 1 : ai_ny )
     lpr_g ( ai_nx , 1 : ai_ny ) = apr_g ( 1 , 1 : ai_ny )		
     
     
     call spli2d_custom ( &
          ai_nx, &
          ai_kx, &
          lpr_taux,&
          ai_ny,&
          ai_ky, &
          apr_tauy, &
          lpr_g,&
          apr_Bcoef,&
          apr_tx,&
          apr_ty )
     
   end subroutine spli2d_perdir
   
   subroutine spli2d_dirper (&
        ai_nx,&
        ai_kx,&
        apr_taux,&
        ar_L, &
        ai_ny,&
        ai_ky, &
        apr_tauy,&
        apr_g,&
        apr_Bcoef,&
        apr_tx,&
        apr_ty )
     ! CALLED WHEN WE WANT TO INTERPOL WITH A PERIODIC second PARAM WITH A PERIOD = ar_L
     implicit none
     ! INPUT
     sll_real64 :: ar_L
     sll_int32  :: ai_nx, ai_kx, ai_ny, ai_ky
     sll_real64, dimension ( :),pointer :: apr_taux ! ai_nx
     sll_real64, dimension (:),pointer :: apr_tauy !  ai_ny -1
     sll_real64, dimension ( :,:) :: apr_g ! ai_nx , ai_ny-1
     ! OUTPUT
     sll_real64, dimension ( ai_nx , ai_ny	) :: apr_Bcoef
     sll_real64, dimension ( ai_nx + ai_kx	) :: apr_tx
     sll_real64, dimension ( ai_ny + ai_ky ) :: apr_ty
     ! LOCAL VARIABLES
     sll_real64, dimension ( ai_nx	) :: lpr_tauy
     sll_real64, dimension ( ai_nx ,ai_ny) :: lpr_g
     
     if ( ar_L == 0.0_8 ) then
        print*,'Error spli2d_per : called with a period = 0 '
        stop
     end if
     
     
     lpr_tauy ( 1 : ai_ny - 1 ) = apr_tauy ( 1 : ai_ny - 1 )
     lpr_tauy ( ai_ny ) = apr_tauy ( 1 ) + ar_L
     
     lpr_g ( 1 : ai_nx , 1 : ai_ny -1 ) = apr_g ( 1 : ai_nx , 1 : ai_ny -1)
     lpr_g (1: ai_nx , ai_ny ) = apr_g ( 1 : ai_nx, 1 )
     
     call spli2d_custom (&
          ai_nx,&
          ai_kx,&
          apr_taux,&
          ai_ny, &
          ai_ky,&
          lpr_tauy, &
          lpr_g, &
          apr_Bcoef,&
          apr_tx,&
          apr_ty )
  
   end subroutine spli2d_dirper
   

   subroutine spli2d_perper(&
     ar_Lx,&
     ai_nx,&
     ai_kx,&
     apr_taux,&
     ar_Ly,&
     ai_ny,&
     ai_ky,&
     apr_tauy,&
     apr_g,&
     apr_Bcoef,&
     apr_tx,&
     apr_ty )
  ! CALLED WHEN WE WANT TO INTERPOL WITH A PERIODIC FIRST PARAM 
     !WITH A PERIOD = ar_L
     implicit none
     ! INPUT
     sll_real64 :: ar_Lx
     sll_real64 :: ar_Ly
     sll_int32  :: ai_nx, ai_kx, ai_ny, ai_ky
     sll_real64, dimension ( :),pointer :: apr_taux ! ai_nx -1
     sll_real64, dimension ( : ),pointer :: apr_tauy ! ai_ny - 1
     sll_real64, dimension (:,:),pointer :: apr_g !  ai_nx  - 1, ai_ny - 1
     ! OUTPUT
     sll_real64, dimension ( ai_nx , ai_ny) :: apr_Bcoef
     sll_real64, dimension ( ai_nx + ai_kx) :: apr_tx
     sll_real64, dimension ( ai_ny + ai_ky) :: apr_ty
     ! LOCAL VARIABLES
     sll_real64, dimension ( ai_nx) :: lpr_taux
     sll_real64, dimension ( ai_ny) :: lpr_tauy
     sll_real64, dimension ( ai_nx ,ai_ny) :: lpr_g
     
     if ( ar_Lx == 0.0_8 ) then
        print*,'Error spli2d_perper : called with a period = 0 '
        stop
     end if
     if ( ar_Ly == 0.0_8 ) then
        print*,'Error spli2d_perper : called with a period = 0 '
        stop
     end if
  
     lpr_taux ( 1 : ai_nx - 1 ) = apr_taux ( 1 : ai_nx - 1 )
     lpr_taux ( ai_nx ) = apr_taux ( 1 ) + ar_Lx
     
     
     lpr_tauy ( 1 : ai_ny - 1 ) = apr_tauy ( 1 : ai_ny - 1 )
     lpr_tauy ( ai_ny ) = apr_tauy ( 1 ) + ar_Ly
     
     
     lpr_g ( 1 : ai_nx - 1 , 1 : ai_ny - 1 ) = &
          apr_g ( 1 : ai_nx - 1 , 1 : ai_ny -1 )
     lpr_g ( ai_nx , 1 : ai_ny -1 ) = apr_g ( 1 , 1 : ai_ny -1 )
     lpr_g ( 1 : ai_nx -1 , ai_ny ) = apr_g ( 1 : ai_nx -1, 1 )
     lpr_g ( ai_nx , ai_ny ) = apr_g ( 1 , 1 )

     
     call spli2d_custom ( &
          ai_nx,&
          ai_kx,&
          lpr_taux,&
          ai_ny,&
          ai_ky,&
          lpr_tauy,&
          lpr_g,&
          apr_Bcoef,&
          apr_tx,&
          apr_ty )
     
end subroutine spli2d_perper

   
 end module sll_module_deboor_splines_2d
