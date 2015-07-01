subroutine banslv ( w, nroww, nrow, nbandl, nbandu, b )

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
  implicit none

  integer nrow
  integer nroww

  real ( kind = 8 ) b(nrow)
  integer i
  integer j
  integer jmax
  integer middle
  integer nbandl
  integer nbandu
  real ( kind = 8 ) w(nroww,nrow)

  middle = nbandu + 1

  !print*, "ze",size(w,1),size(w,2)

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

  !print*, 'banslv',b
  return
end subroutine banslv
