module sll_module_deboor_splines_2d

#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_utilities.h"
use sll_module_deboor_splines_1d
  
implicit none 
private

public dvalue2d, spli2d_custom
public bvalue2d, spli2d_custom_derder
  
contains
  
subroutine bvalue2d( x,       &
                     y,       &
                     nx,      &
                     kx,      &
                     ny,      &
                     ky,      &
                     bcoef,   &
                     tx,      &
                     ty,      &
                     val )

  sll_real64,                 intent(in)  :: x
  sll_real64,                 intent(in)  :: y
  sll_int32,                  intent(in)  :: nx
  sll_int32,                  intent(in)  :: kx
  sll_int32,                  intent(in)  :: ny
  sll_int32,                  intent(in)  :: ky
  sll_real64, dimension(:),   intent(in)  :: tx     !nx + kx 
  sll_real64, dimension(:),   intent(in)  :: ty     !ny + ky	
  sll_real64, dimension(:,:), intent(in)  :: bcoef  !( nx,ny)
  sll_real64,                 intent(out) :: val
  
  sll_int32  :: j
  sll_int32  :: lflag
  sll_int32  :: lefty
  sll_real64, dimension(1:ky) :: lpr_coef 
   
  call interv ( ty, ny+ky, y, lefty, lflag )

  if ( lflag .NE. 0 ) then
    val = 0.0_8
    return 
  end if
    
  do j = 1, ky
       
    lpr_coef(j) = bvalue(tx,                     &
                         bcoef(1:nx,lefty-ky+j), &
                               nx,               &
                               kx,               &
                               x,                &
                               0 )
       
  end do
   
  val = bvalue(ty(lefty-ky+1:lefty+ky), &
               lpr_coef,                &
               ky,                      &
               ky,                      &
               y,                       &
               0 )

end subroutine bvalue2d

function dvalue2d( x,      &
                   y,      &
                   nx,     &
                   kx,     &
                   ny,     &
                   ky,     &
                   bcoef,  &
                   tx,     &
                   ty,     &
                   deriv1, &
                   deriv2 ) 

  sll_real64, intent(in) :: x
  sll_real64, intent(in) :: y
  sll_int32,  intent(in) :: nx
  sll_int32,  intent(in) :: kx
  sll_int32,  intent(in) :: ny
  sll_int32,  intent(in) :: ky
  sll_real64, intent(in) :: bcoef(:,:) ! dimension(nx,ny)
  sll_real64, intent(in) :: tx(:)      ! dimension(nx+kx)
  sll_real64, intent(in) :: ty(:)      ! dimension(ny+ky)
  sll_int32,  intent(in) :: deriv1
  sll_int32,  intent(in) :: deriv2

  sll_real64 :: dvalue2d

  sll_int32  :: j
  sll_int32  :: lflag
  sll_int32  :: lefty

  sll_real64, dimension (1:ky),target:: lpr_coef ! ky
    
  call interv( ty, ny + ky, y, lefty, lflag)
    
  if ( lflag .NE. 0 ) then
    dvalue2d = 0.0_8
    return 
  end if
    
  do j = 1, ky
       
    lpr_coef(j) = bvalue( tx,                     &
                          bcoef(1:nx,lefty-ky+j), &
                          nx,                     &
                          kx,                     &
                          x,                      &
                          deriv1 )
       
  end do

  dvalue2d = bvalue( ty(lefty-ky+1:lefty+ky), &
                     lpr_coef,                &
                     ky,                      &
                     ky,                      &
                     y,                       &
                     deriv2 )

end function dvalue2d
  

!*****************************************************************************80
!
! SPLI2D produces a interpolatory tensor product spline.
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
subroutine spli2d ( tau, gtau, t, n, k, m, bcoef)
    
  sll_real64, dimension(:),   intent(in) :: tau   !(n)
  sll_real64, dimension(:,:), intent(in) :: gtau  !(n,m)
  sll_real64, dimension(:),   intent(in) :: t     !(n+k)
  sll_int32,                  intent(in) :: n
  sll_int32,                  intent(in) :: k
  sll_int32,                  intent(in) :: m

  sll_real64, dimension(:,:), intent(inout) :: bcoef !(m,n)

  sll_real64, dimension(n)         :: work  !(n)
  sll_real64, dimension((2*k-1),n) :: q     !((2*k-1)*n)

  sll_int32   :: i
  sll_int32   :: ilp1mx
  sll_int32   :: j
  sll_int32   :: jj
  sll_int32   :: left
  sll_int32   :: iflag
  sll_real64  :: taui
    
  left = k
   
  q = 0.0_f64
  !
  !  Construct the N interpolation equations.
  !
    
  do i = 1, n
       
    taui = tau(i)
    ilp1mx = min ( i+k, n+1 )
    !
    !  Find the index LEFT in the closed interval (I,I+K-1) such that:
    !
    !    T(LEFT) < = TAU(I) < T(LEFT+1)
    !
    !  The matrix will be singular if this is not possible.
    !
    left = max(left, i)
       
    SLL_ASSERT(taui >= t(left))
    ! Check The TAU array is strictly increasing.
       
    do while ( t(left+1) <= taui )
          
      left = left + 1
         
      if ( left < ilp1mx ) cycle
          
      left = left - 1
          
      SLL_ASSERT(taui >= t(left+1))
      !Check The TAU array is strictly increasing.
          
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
      !q(jj) = work(j)
      q(i-(left+j)+2*k,left+j-k) = work(j)
    end do
       
  end do
    
  !
  !  Factor A, stored again in Q.
  !
  call banfac ( q, k+k-1, n, k-1, k-1, iflag )
    
  if ( iflag == 2 ) then
    SLL_ERROR('BANFAC reports that the matrix is singular.')
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
   
end subroutine spli2d

subroutine spli2d_custom ( nx, kx, taux, ny, ky, tauy, g, bcoef, tx, ty)

  sll_int32,  intent(in)    :: nx
  sll_int32,  intent(in)    :: kx
  sll_real64, intent(in)    :: taux(:)    !nx
  sll_int32,  intent(in)    :: ny
  sll_int32,  intent(in)    :: ky
  sll_real64, intent(in)    :: tauy(:)    !ny	
  sll_real64, intent(in)    :: g(:,:)     !nx,ny	
  sll_real64, intent(inout) :: bcoef(:,:) !nx , ny 

  sll_real64, intent(in) :: tx(:) 
  sll_real64, intent(in) :: ty(:) 
  sll_real64 :: tmp(nx,ny)
    
  call spli2d(taux,   g, tx, nx, kx, ny, tmp)
  call spli2d(tauy, tmp, ty, ny, ky, nx, bcoef)
     
end subroutine spli2d_custom

!
!  Parameters:
!
!    Input, real ( kind = 8 ) TAU(N), contains the data point abscissas.
!    TAU must be strictly increasing
!
!    Input, real ( kind = 8 ) GTAU(N,M), contains the data point ordinates.
!    Input, integer (kind= 8 )TAU_DER(Np), contains the data point abscissas.
!    TAU must be strictly increasing
!
!    Input, real ( kind = 8 ) GTAU_DER(Np,M),contains the data point ordinates.
!
!    Input, real ( kind = 8 ) T(N+Np+K), the knot sequence.
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
subroutine spli2d_custom_derder ( nx,       &
                                  nx_der,   &
                                  kx,       &
                                  taux,     &
                                  taux_der, &
                                  ny,       &
                                  ny_der,   &
                                  ky,       &
                                  tauy,     &
                                  tauy_der, &
                                  g,        &
                                  g_der1,   &
                                  g_der2,   &
                                  bcoef,    &
                                  tx,       &
                                  ty )

  sll_int32,  intent(in)    :: nx
  sll_int32,  intent(in)    :: nx_der
  sll_int32,  intent(in)    :: kx
  sll_real64, intent(in)    :: taux(:)
  sll_int32,  intent(in)    :: tauy_der(:)
  sll_int32,  intent(in)    :: ny
  sll_int32,  intent(in)    :: ny_der
  sll_int32,  intent(in)    :: ky
  sll_real64, intent(in)    :: tauy(:)
  sll_int32,  intent(in)    :: taux_der(:)
  sll_real64, intent(in)    :: g(:,:)    
  sll_real64, intent(in)    :: g_der1(:,:)
  sll_real64, intent(in)    :: g_der2(:,:)

  sll_real64, intent(inout) :: bcoef(:,:)
  sll_real64, intent(inout) :: tx(:) 
  sll_real64, intent(inout) :: ty(:) 

  sll_real64, dimension(1:ny,1:nx+nx_der) :: tmp
  sll_int32 :: i
  sll_int32 :: j
   
  SLL_ASSERT(nx+nx_der+kx == nx+2*(kx-1))
  SLL_ASSERT(ny+ny_der+ky == ny+2*(ky-1))
    
  do j = 1, ny
       
    call splint_der(taux,g(:,j),taux_der,g_der1(:,j),tx,nx,nx_der,kx,tmp(j,:))
       
  end do

  do i = 1, nx+nx_der
       
    call splint_der(tauy,tmp(:,i),tauy_der,g_der2(:,i),ty,ny,ny_der,ky,bcoef(i,:))
       
  end do
    
end subroutine spli2d_custom_derder

end module sll_module_deboor_splines_2d
