subroutine bsplpp ( t, bcoef, n, k, scrtch, break, coef, l )

!*****************************************************************************80
!
!! BSPLPP converts from B-spline to piecewise polynomial form.
!
!  Discussion:
!
!    The B-spline representation of a spline is 
!      ( T, BCOEF, N, K ),
!    while the piecewise polynomial representation is 
!      ( BREAK, COEF, L, K ).
!
!    For each breakpoint interval, the K relevant B-spline coefficients 
!    of the spline are found and then differenced repeatedly to get the 
!    B-spline coefficients of all the derivatives of the spline on that 
!    interval. 
!
!    The spline and its first K-1 derivatives are then evaluated at the 
!    left end point of that interval, using BSPLVB repeatedly to obtain 
!    the values of all B-splines of the appropriate order at that point.
!
!  Modified:
!
!    14 February 2007
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
!    Input, real ( kind = 8 ) T(N+K), the knot sequence.
! 
!    Input, real ( kind = 8 ) BCOEF(N), the B spline coefficient sequence.
! 
!    Input, integer ( kind = 4 ) N, the number of B spline coefficients.
! 
!    Input, integer ( kind = 4 ) K, the order of the spline.
! 
!    Work array, real ( kind = 8 ) SCRTCH(K,K).
! 
!    Output, real ( kind = 8 ) BREAK(L+1), the piecewise polynomial breakpoint 
!    sequence.  BREAK contains the distinct points in the sequence T(K:N+1)
! 
!    Output, real ( kind = 8 ) COEF(K,N), with COEF(I,J) = (I-1)st derivative 
!    of the spline at BREAK(J) from the right.
! 
!    Output, integer ( kind = 4 ) L, the number of polynomial pieces which 
!    make up the spline in the interval ( T(K), T(N+1) ).
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) n

  real ( kind = 8 ) bcoef(n)
  real ( kind = 8 ) biatx(k)
  real ( kind = 8 ) break(*)
  real ( kind = 8 ) coef(k,n)
  real ( kind = 8 ) diff
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jp1
  integer ( kind = 4 ) left
  integer ( kind = 4 ) lsofar
  real ( kind = 8 ) scrtch(k,k)      
  real ( kind = 8 ) sum1
  real ( kind = 8 ) t(n+k)

  lsofar = 0
  break(1) = t(k)
  
  do left = k, n
!
!  Find the next nontrivial knot interval.
!
    if ( t(left+1) == t(left) ) then
      cycle
    end if

    lsofar = lsofar + 1
    break(lsofar+1) = t(left+1)

    if ( k <= 1 ) then
      coef(1,lsofar) = bcoef(left)
      cycle
    end if
!
!  Store the K B-spline coefficients relevant to current knot 
!  interval in SCRTCH(*,1).
!
    do i = 1, k
      scrtch(i,1) = bcoef(left-k+i)
    end do
!
!  For J=1,...,K-1, compute the  K-J  B-spline coefficients relevant to
!  the current knot interval for the J-th derivative by differencing
!  those for the (J-1)st derivative, and store in SCRTCH(.,J+1).
!
    do jp1 = 2, k
      j = jp1 - 1
      do i = 1, k - j
        diff = t(left+i) - t(left+i-(k-j))
        if ( 0.0D+00 < diff ) then
          scrtch(i,jp1) = ( ( scrtch(i+1,j) - scrtch(i,j) ) / diff ) &
            * real ( k - j, kind = 8 )
        end if
      end do
    end do
!
!  For J=0, ..., K-1, find the values at T(left) of the J+1
!  B-splines of order J+1 whose support contains the current
!  knot interval from those of order J (in  BIATX ), then combine
!  with the B-spline coefficients (in SCRTCH(.,K-J) ) found earlier
!  to compute the (K-J-1)st derivative at  T(LEFT) of the given
!  spline.
!
    call bsplvb ( t, 1, 1, t(left), left, biatx )

    coef(k,lsofar) = scrtch(1,k)
    
    do jp1 = 2, k
    
      call bsplvb ( t, jp1, 2, t(left), left, biatx )

      coef(k+1-jp1,lsofar) = dot_product ( biatx(1:jp1), scrtch(1:jp1,k+1-jp1) )
      
    end do

  end do
   
  l = lsofar

  return
end subroutine

