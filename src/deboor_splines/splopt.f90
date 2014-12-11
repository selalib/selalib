subroutine splopt ( tau, n, k, scrtch, t, iflag )
  
!*****************************************************************************80
  !
!! SPLOPT computes the knots for an optimal recovery scheme. 
  !
  !  Discussion:
  !
  !    The optimal recovery scheme is of order K for data at TAU(1:N).
  !
  !    The interior knots T(K+1:N) are determined by Newton's method in 
  !    such a way that the signum function which changes sign at
  !      T(K+1:N)  and nowhere else in ( TAU(1), TAU(N) ) is 
  !    orthogonal to the spline space SPLINE ( K, TAU ) on that interval.
  !
  !    Let XI(J) be the current guess for T(K+J), J=1,...,N-K.  Then
  !    the next Newton iterate is of the form
  !
  !      XI(J) + (-)**(N-K-J)*X(J),  J=1,...,N-K,
  !
  !    with X the solution of the linear system
  !
  !      C * X = D.
  !
  !    Here, for all J,
  !
  !      C(I,J) = B(I)(XI(J)), 
  !
  !    with B(I) the I-th B-spline of order K for the knot sequence TAU, 
  !    for all I, and D is the vector given, for each I, by
  !
  !      D(I) = sum ( -A(J), J=I,...,N ) * ( TAU(I+K) - TAU(I) ) / K,
  !
  !    with, for I = 1 to N-1:  
  !
  !      A(I) = sum ( (-)**(N-K-J)*B(I,K+1,TAU)(XI(J)), J=1,...,N-K )
  !
  !    and  
  !
  !      A(N) = -0.5.
  !
  !    See Chapter XIII of text and references there for a derivation.
  !
  !    The first guess for T(K+J) is sum ( TAU(J+1:J+K-1) ) / ( K - 1 ).
  !
  !    The iteration terminates if max ( abs ( X(J) ) ) < TOL, with
  !
  !      TOL = TOLRTE * ( TAU(N) - TAU(1) ) / ( N - K ),
  !
  !    or else after NEWTMX iterations.
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
  !    Input, real ( kind = 8 ) TAU(N), the interpolation points.
  !    assumed to be nondecreasing, with TAU(I) < TAU(I+K), for all I.
  !
  !    Input, integer ( kind = 4 ) N, the number of data points.
  !
  !    Input, integer ( kind = 4 ) K, the order of the optimal recovery scheme 
  !    to be used.
  !
  !    Workspace, real ( kind = 8 ) SCRTCH((N-K)*(2*K+3)+5*K+3).  The various
  !    contents are specified in the text below.
  !
  !    Output, real ( kind = 8 ) T(N+K), the optimal knots ready for
  !    use in optimal recovery.  Specifically, T(1:K) = TAU(1), 
  !    T(N+1:N+K) = TAU(N), while the N - K interior knots T(K+1:N) 
  !    are calculated.
  !
  !    Output, integer ( kind = 4 ) IFLAG, error indicator.
  !    = 1, success.  T contains the optimal knots.
  !    = 2, failure.  K < 3 or N < K or the linear system was singular.
  !
  implicit none
  
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n
  
  real ( kind = 8 ) del
  real ( kind = 8 ) delmax
  real ( kind = 8 ) floatk
  integer ( kind = 4 ) i
  integer ( kind = 4 ) id
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) index
  integer ( kind = 4 ) j
  integer ( kind = 4 ) kp1
  integer ( kind = 4 ) kpkm1
  integer ( kind = 4 ) kpn
  integer ( kind = 4 ) left
  integer ( kind = 4 ) leftmk
  integer ( kind = 4 ) lenw
  integer ( kind = 4 ) ll
  integer ( kind = 4 ) llmax
  integer ( kind = 4 ) llmin
  integer ( kind = 4 ) na
  integer ( kind = 4 ) nb
  integer ( kind = 4 ) nc
  integer ( kind = 4 ) nd
  integer ( kind = 4 ), parameter :: newtmx = 10
  integer ( kind = 4 ) newton
  integer ( kind = 4 ) nmk
  integer ( kind = 4 ) nx
  real ( kind = 8 ) scrtch((n-k)*(2*k+3)+5*k+3)
  real ( kind = 8 ) t(n+k)
  real ( kind = 8 ) tau(n)
  real ( kind = 8 ) sign
  real ( kind = 8 ) signst
  real ( kind = 8 ) tol
  real ( kind = 8 ), parameter :: tolrte = 0.000001D+00
  real ( kind = 8 ) xij

  nmk = n - k
  
  if ( n < k ) then
     write ( *, '(a)' ) ' '
     write ( *, '(a)' ) 'SPLOPT - Fatal error!'
     write ( *, '(a)' ) '  N < K.'
     iflag = 2
     return
  end if
  
  if ( n == k ) then
     t(1:k) = tau(1)
     t(n+1:n+k) = tau(n)
     return
  end if
  
  if ( k <= 2 ) then
     write ( *, '(a)' ) ' '
     write ( *, '(a)' ) 'SPLOPT - Fatal error!'
     write ( *, '(a)' ) '  K < 2.'
     iflag = 2
     stop
  end if
 
  floatk = k
  kp1 = k + 1
  kpkm1 = k + k - 1
  kpn = k + n
  
  signst = -1.0D+00
  if ( ( nmk / 2 ) * 2 < nmk ) then
     signst = 1.0D+00
  end if
  !
  !  SCRTCH(I) = TAU-EXTENDED(I), I=1,...,N+K+K
  !
  nx = n + k + k + 1
  !
  !  SCRTCH(I+NX) = XI(I), I=0,...,N-K+1
  !
  na = nx + nmk + 1
  !
  !  SCRTCH(I+NA) = - A(I), I=1,...,N
  !
  nd = na + n
  !
  !  SCRTCH(I+ND) = X(I) or D(I), I=1,...,N-K
  !
  nb = nd + nmk
  !
  !  SCRTCH(I+NB) = BIATX(I), I=1,...,K+1
  !
  nc = nb + kp1
  !
  !  SCRTCH(I+(J-1)*(2K-1)+NC) = W(I,J) = C(I-K+J,J), I=J-K,...,J+K,
  !                                                     J=1,...,N-K.
  !
  lenw = kpkm1 * nmk
  !
  !  Extend TAU to a knot sequence and store in SCRTCH.
  !
  scrtch(1:k) = tau(1)
  scrtch(k+1:k+n) = tau(1:n)
  scrtch(kpn+1:kpn+k) = tau(n)
  !
  !  First guess for SCRTCH (.+NX) = XI.
  !
  scrtch(nx) = tau(1)
  scrtch(nmk+1+nx) = tau(n)
 
  do j = 1, nmk 
     scrtch(j+nx) = sum ( tau(j+1:j+k-1) ) / real ( k - 1, kind = 8 )
  end do
  !
  !  Last entry of SCRTCH (.+NA) = -A  is always ...
  !
  scrtch(n+na) = 0.5D+00
  !
  !  Start the Newton iteration.
  !
  newton = 1
  tol = tolrte * ( tau(n) - tau(1) ) / real ( nmk, kind = 8 )
  !
  !  Start the Newton step.
  !  Compute the 2*K-1 bands of the matrix C and store in SCRTCH(.+NC),
  !  and compute the vector SCRTCH(.+NA) = -A.
  !
  do newton = 1, newtmx
     
    scrtch(nc+1:nc+lenw) = 0.0D+00
    scrtch(na+1:na+n-1) = 0.0D+00
    
    sign = signst
    left = kp1
    
    do j = 1, nmk
       
       xij = scrtch(j+nx)
       
      do
         
         if ( xij < scrtch(left+1) ) then
            exit
         end if

         left = left + 1
         
         if ( kpn <= left ) then
            left = left - 1
            exit
         end if

      end do
      
      call bsplvb ( scrtch, k, 1, xij, left, scrtch(1+nb) )
      !
      !  The TAU sequence in SCRTCH is preceded by K additional knots.
      !
      !  Therefore, SCRTCH(LL+NB) now contains B(LEFT-2K+LL)(XIJ)
      !  which is destined for C(LEFT-2K+LL,J), and therefore for
      !
      !    W(LEFT-K-J+LL,J)= SCRTCH(LEFT-K-J+LL+(J-1)*KPKM1 + NC)
      !
      !  since we store the 2*K-1 bands of C in the 2*K-1 rows of
      !  the work array W, and W in turn is stored in SCRTCH,
      !  with W(1,1) = SCRTCH(1+NC).
      !
      !  Also, C being of order N - K, we would want  
      !    1 <= LEFT-2K+LL <= N - K or
      !    LLMIN=2K-LEFT  <=  LL <= N-LEFT+K = LLMAX.
      !
      leftmk = left - k
      index = leftmk - j + ( j - 1 ) * kpkm1 + nc
      llmin = max ( 1, k - leftmk )
      llmax = min ( k, n - leftmk )
      do ll = llmin, llmax
         scrtch(ll+index) = scrtch(ll+nb)
      end do
      
      call bsplvb ( scrtch, kp1, 2, xij, left, scrtch(1+nb) )
      
      id = max ( 0, leftmk - kp1 )
      llmin = 1 - min ( 0, leftmk - kp1 )
      do ll = llmin, kp1
         id = id + 1
         scrtch(id+na) = scrtch(id+na) - sign * scrtch(ll+nb)
      end do
     
      sign = - sign
      
   end do

   call banfac ( scrtch(1+nc), kpkm1, nmk, k-1, k-1, iflag )
   
   if ( iflag == 2 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPLOPT - Fatal error!'
      write ( *, '(a)' ) '  Matrix C is not invertible.'
      stop
   end if
   !
   !  Compute SCRTCH(.+ND) = D from SCRTCH(.+NA) = -A.
   !
   do i = n, 2, -1
      scrtch(i-1+na) = scrtch(i-1+na) + scrtch(i+na)
   end do
   
   do i = 1, nmk
      scrtch(i+nd) = scrtch(i+na) * ( tau(i+k) - tau(i) ) / floatk
   end do
   !
!  Compute SCRTCH(.+ND)= X.
   !
   call banslv ( scrtch(1+nc), kpkm1, nmk, k-1, k-1, scrtch(1+nd) )
   !
   !  Compute SCRTCH(.+ND) = change in XI.  Modify, if necessary, to
   !  prevent new XI from moving more than 1/3 of the way to its
   !  neighbors.  Then add to XI to obtain new XI in SCRTCH(.+NX).
   !
   delmax = 0.0D+00
   sign = signst
   
   do i = 1, nmk

      del = sign * scrtch(i+nd)
      delmax = max ( delmax, abs ( del ) )
      
      if ( 0.0D+00 < del ) then
         del = min ( del, ( scrtch(i+1+nx) - scrtch(i+nx) ) / 3.0D+00 )
      else
         del = max ( del, ( scrtch(i-1+nx) - scrtch(i+nx) ) / 3.0D+00 )
      end if
     
      sign = - sign
      scrtch(i+nx) = scrtch(i+nx) + del
      
   end do
   !
   !  Call it a day in case change in XI was small enough or too many
   !  steps were taken.
   !
   if ( delmax < tol ) then
      exit
   end if
   
end do

if ( tol <= delmax ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLOPT - Warning!'
    write ( *, '(a)' ) '  The Newton iteration did not converge.'
 end if
 
  t(1:k) = tau(1)
  t(k+1:n) = scrtch(nx+1:nx+n-k)  
  t(n+1:n+k) = tau(n)
  
  return
end subroutine splopt
