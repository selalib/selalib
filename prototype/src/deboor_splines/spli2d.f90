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

  integer m
  integer n

  real ( kind = 8 ) bcoef(m,n)
  real ( kind = 8 ) gtau(n,m)
  integer i
  integer iflag
  integer ilp1mx
  integer j
  integer jj
  integer k
  integer left
  real ( kind = 8 ) q((2*k-1)*n)
  real ( kind = 8 ) t(n+k)
  real ( kind = 8 ) tau(n)
  real ( kind = 8 ) taui
  real ( kind = 8 ) work(n)

  left = k

  !print*, t
  q(1:(2*k-1)*n) = 0.0D+00
!
!  Construct the N interpolation equations.
!

  print*, n, tau(1:n)
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
  ! print*, 'qqqq_spli2d', q
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
 ! print*, 'rrttt',q
!
!  Solve
!
!    A * BCOEF = GTAU
!
!  by back substitution.

  !print*, 'za',gtau(:,2)
  !print*, "gt",size(q,1)!,size(q,2)
  do j = 1, m

    work(1:n) = gtau(1:n,j)
    
    
    call banslv ( q, k+k-1, n, k-1, k-1, work )

    bcoef(j,1:n) = work(1:n)
   ! print*, 'uyt',work(1:n)
  end do
  !print*,  bcoef(2,1:n)

  return
end subroutine spli2d
