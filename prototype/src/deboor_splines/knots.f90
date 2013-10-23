subroutine knots ( break, l, kpm, m, t, n )
  
!*****************************************************************************80
  !
  !! KNOTS is to be called in COLLOC.
  !
  !  Discussion:
  !
  !    Note that the FORTRAN77 calling sequence has been modified, by
  !    adding the variable M.
  !
!    From the given breakpoint sequence BREAK, this routine constructs the 
  !    knot sequence T so that
!
!      SPLINE(K+M,T) = PP(K+M,BREAK) 
!
!    with M-1 continuous derivatives.
!
!    This means that T(1:N+KPM) is equal to BREAK(1) KPM times, then 
!    BREAK(2) through BREAK(L) each K times, then, finally, BREAK(L+1) 
!    KPM times.
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
!    Input, real ( kind = 8 ) BREAK(L+1), the breakpoint sequence.
!
!    Input, integer ( kind = 4 ) L, the number of intervals or pieces.
!
!    Input, integer ( kind = 4 ) KPM, = K+M, the order of the piecewise
!    polynomial function or spline.
!
!    Input, integer ( kind = 4 ) M, the order of the differential equation.
!
!    Output, real ( kind = 8 ) T(N+KPM), the knot sequence.
!
!    Output, integer ( kind = 4 ) N, = L*K+M = the dimension of SPLINE(K+M,T).
!
  implicit none

  integer ( kind = 4 ) kpm
  integer ( kind = 4 ) l
  integer ( kind = 4 ) n

  real ( kind = 8 ) break(l+1)
  integer ( kind = 4 ) iside
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) jjj
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ll
  integer ( kind = 4 ) m
  real ( kind = 8 ) t(n+kpm)
  real ( kind = 8 ) xside

  k = kpm - m
  n = l * k + m
  jj = n + kpm
  jjj = l + 1
  
  do ll = 1, kpm
    t(jj) = break(jjj)
    jj = jj - 1
  end do
  
  do j = 1, l
    jjj = jjj - 1
    do ll = 1, k
      t(jj) = break(jjj)
      jj = jj - 1
    end do
  end do
   
  t(1:kpm) = break(1)

  return
end
