module banded_linear_system_solver
implicit none

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

contains

subroutine banfac ( w, nroww, nrow, nbandl, nbandu, iflag )


  integer, intent(in)    :: nrow
  integer, intent(in)    :: nroww
  real(8), intent(inout) :: w(nroww,nrow)
  integer, intent(in)    :: nbandl
  integer, intent(in)    :: nbandu
  integer, intent(out)   :: iflag

  real ( kind = 8 ) factor
  integer i
  integer j
  integer k
  integer middle
  real ( kind = 8 ) pivot

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
       !do j = 1, nrow
       !   print*, 'col',j, 'values', w(:,j)
       !end do
       !print*, middle
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

!*****************************************************************************80
!
!! BANSLV solves a banded linear system A * X = B factored by BANFAC.
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
!    Input, real ( kind = 8 ) W(NROWW,NROW).  W contains the banded matrix,
!    after it has been factored by BANFAC.
!
!    Input, integer NROWW, the row dimension of the work array W.
!    NROWW must be at least NBANDL+1 + NBANDU.
!
!    Input, integer NROW, the number of rows in A.
!
!    Input, integer NBANDL, the number of bands of A below the
!    main diagonal.
!
!    Input, integer NBANDU, the number of bands of A above the
!    main diagonal.
!
!    Input/output, real ( kind = 8 ) B(NROW).
!    On input, B contains the right hand side of the system to be solved.
!    On output, B contains the solution, X.
!
subroutine banslv ( w, nroww, nrow, nbandl, nbandu, b )

  implicit none

  integer, intent(in)    :: nroww
  integer, intent(in)    :: nrow
  real(8), intent(in)    :: w(nroww,nrow)
  integer, intent(in)    :: nbandl
  integer, intent(in)    :: nbandu
  real(8), intent(inout) :: b(nrow)

  integer :: i
  integer :: j
  integer :: jmax
  integer :: middle

  middle = nbandu + 1

  if ( nrow == 1 ) then
    b(1) = b(1) / w(middle,1)
    return
  end if
!
!  Forward pass:
!
!  For I = 1, 2, ..., NROW-1, subtract RHS(I)*(I-th column of L)
!  from the right hand side, below the I-th row.
!
  if ( 0 < nbandl ) then
    do i = 1, nrow-1
      jmax = min ( nbandl, nrow-i )
      do j = 1, jmax
        b(i+j) = b(i+j) - b(i) * w(middle+j,i)
      end do
    end do
  end if
!
!  Backward pass:
!
!  For I=NROW, NROW-1,...,1, divide RHS(I) by
!  the I-th diagonal entry of U, then subtract
!  RHS(I)*(I-th column of U) from right hand side, above the I-th row.
!
  do i = nrow, 2, -1

    b(i) = b(i) / w(middle,i)

    do j = 1, min ( nbandu, i-1 )
      b(i-j) = b(i-j) - b(i) * w(middle-j,i)
    end do

  end do

  b(1) = b(1) / w(middle,1)

  return
end subroutine banslv

end module banded_linear_system_solver
  
