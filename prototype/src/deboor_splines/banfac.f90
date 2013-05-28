subroutine banfac ( w, nroww, nrow, nbandl, nbandu, iflag )

!*****************************************************************************80
!
!! BANFAC factors a banded matrix without pivoting.
!
!  Discussion:
!
!    BANFAC returns in W the LU-factorization, without pivoting, of
!    the banded matrix A of order NROW with (NBANDL+1+NBANDU) bands
!    or diagonals in the work array W.
!
!    Gauss elimination without pivoting is used.  The routine is
!    intended for use with matrices A which do not require row
!    interchanges during factorization, especially for the totally
!    positive matrices which occur in spline calculations.
!
!    The matrix storage mode used is the same one used by LINPACK
!    and LAPACK, and results in efficient innermost loops.
!
!    Explicitly, A has
!
!      NBANDL bands below the diagonal
!      1     main diagonal
!      NBANDU bands above the diagonal
!
!    and thus, with MIDDLE=NBANDU+1,
!    A(I+J,J) is in W(I+MIDDLE,J) for I=-NBANDU,...,NBANDL, J=1,...,NROW.
!
!    For example, the interesting entries of a banded matrix
!    matrix of order 9, with NBANDL=1, NBANDU=2:
!
!      11 12 13  0  0  0  0  0  0
!      21 22 23 24  0  0  0  0  0
!       0 32 33 34 35  0  0  0  0
!       0  0 43 44 45 46  0  0  0
!       0  0  0 54 55 56 57  0  0
!       0  0  0  0 65 66 67 68  0
!       0  0  0  0  0 76 77 78 79
!       0  0  0  0  0  0 87 88 89
!       0  0  0  0  0  0  0 98 99
!
!    would appear in the first 1+1+2=4 rows of W as follows:
!
!       0  0 13 24 35 46 57 68 79
!       0 12 23 34 45 56 67 78 89
!      11 22 33 44 55 66 77 88 99
!      21 32 43 54 65 76 87 98  0
!
!    All other entries of W not identified in this way with an
!    entry of A are never referenced.
!
!    This routine makes it possible to solve any particular linear system
!    A*X=B for X by the call
!
!      call banslv ( w, nroww, nrow, nbandl, nbandu, b )
!
!    with the solution X contained in B on return.
!
!    If IFLAG=2, then one of NROW-1, NBANDL, NBANDU failed to be nonnegative,
!    or else one of the potential pivots was found to be zero
!    indicating that A does not have an LU-factorization.  This
!    implies that A is singular in case it is totally positive.
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
!    Input/output, real ( kind = 8 ) W(NROWW,NROW).
!    On input, W contains the "interesting" part of a banded
!    matrix A, with the diagonals or bands of A stored in the
!    rows of W, while columns of A correspond to columns of W.
!    On output, W contains the LU-factorization of A into a unit
!    lower triangular matrix L and an upper triangular matrix U
!    (both banded) and stored in customary fashion over the
!    corresponding entries of A.
!
!    Input, integer NROWW, the row dimension of the work array W.
!    NROWW must be at least NBANDL+1 + NBANDU.
!
!    Input, integer NROW, the number of rows in A.
!
!    Input, integer NBANDL, the number of bands of A below the main diagonal.
!
!    Input, integer NBANDU, the number of bands of A above the main diagonal.
!
!    Output, integer IFLAG, error flag.
!    1, success.
!    2, failure, the matrix was not factored.
!
  implicit none

  integer nrow
  integer nroww

  real ( kind = 8 ) factor
  integer i
  integer iflag
  integer j
  integer k
  integer middle
  integer nbandl
  integer nbandu
  real ( kind = 8 ) pivot
  real ( kind = 8 ) w(nroww,nrow)

  iflag = 1

  if ( nrow < 1 ) then
    iflag = 2
    return
  end if
!
!  W(MIDDLE,*) contains the main diagonal of A.
!
  middle = nbandu + 1

  if ( nrow == 1 ) then
    if ( w(middle,nrow) == 0.0D+00 ) then
      iflag = 2
    end if
    return
  end if
!
!  A is upper triangular.  Check that the diagonal is nonzero.
!
  if ( nbandl <= 0 ) then

    do i = 1, nrow-1
      if ( w(middle,i) == 0.0D+00 ) then
        iflag = 2
        return
      end if
    end do

    if ( w(middle,nrow) == 0.0D+00 ) then
      iflag = 2
    end if

    return
!
!  A is lower triangular.  Check that the diagonal is nonzero and
!  divide each column by its diagonal.
!
  else if ( nbandu <= 0 ) then

    do i = 1, nrow-1

      pivot = w(middle,i)

      if ( pivot == 0.0D+00 ) then
        iflag = 2
        return
      end if

      do j = 1, min ( nbandl, nrow-i )
        w(middle+j,i) = w(middle+j,i) / pivot
      end do

    end do

    return

  end if
!
!  A is not just a triangular matrix.
!  Construct the LU factorization.
!
  do i = 1, nrow-1
!
!  W(MIDDLE,I) is the pivot for the I-th step.
!
    if ( w(middle,i) == 0.0D+00 ) then
      iflag = 2
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BANFAC - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero pivot encountered in column ', i
      stop
    end if
!
!  Divide each entry in column I below the diagonal by PIVOT.
!
    do j = 1, min ( nbandl, nrow-i )
      w(middle+j,i) = w(middle+j,i) / w(middle,i)
    end do
!
!  Subtract A(I,I+K)*(I-th column) from (I+K)-th column (below row I).
!
    do k = 1, min ( nbandu, nrow-i )
      factor = w(middle-k,i+k)
      do j = 1, min ( nbandl, nrow-i )
        w(middle-k+j,i+k) = w(middle-k+j,i+k) - w(middle+j,i) * factor
      end do
    end do

  end do
!
!  Check the last diagonal entry.
!
  if ( w(middle,nrow) == 0.0D+00 ) then
    iflag = 2
  end if

  return
end subroutine banfac