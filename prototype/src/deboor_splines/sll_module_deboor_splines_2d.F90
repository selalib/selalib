module sll_module_deboor_splines_2d

#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_utilities.h"
use sll_module_deboor_splines_1d
  
implicit none 
private

public dvalue2d, spli2d_dirper, spli2d_perdir, spli2d_custom
public bvalue2d, spli2d_custom_derder, spli2d_perper
  
contains
  
! Just for the record: what is ar_x, ar_y, nx, etc., etc.? What is this
! function supposed to do? This is a deeply frustrating file.
subroutine bvalue2d( ar_x,        &
                     ar_y,        &
                     nx,       &
                     kx,       &
                     ny,       &
                     ky,       &
                     bcoef,   &
                     apr_tx,      &
                     apr_ty,      &
                     val )

  sll_real64,                 intent(in)  :: ar_x
  sll_real64,                 intent(in)  :: ar_y
  sll_int32,                  intent(in)  :: nx
  sll_int32,                  intent(in)  :: kx
  sll_int32,                  intent(in)  :: ny
  sll_int32,                  intent(in)  :: ky
  sll_real64, dimension(:),   intent(in)  :: apr_tx !nx + kx 
  sll_real64, dimension(:),   intent(in)  :: apr_ty !ny + ky	
  sll_real64, dimension(:,:), intent(in)  :: bcoef!( nx,ny)
  sll_real64,                 intent(out) :: val
  
  sll_int32  :: li_j, li_mflag, li_lefty
  sll_real64, dimension(1:ky) :: lpr_coef ! ky
   
  call interv ( apr_ty,ny + ky, ar_y, li_lefty, li_mflag )

  if ( li_mflag .NE. 0 ) then
    val = 0.0_8
    return 
  end if
    
  do li_j = 1, ky
       
    lpr_coef ( li_j ) = bvalue(apr_tx,                                 &
                               bcoef(1:nx,li_lefty-ky+li_j), &
                               nx,                                  &
                               kx,                                  &
                               ar_x,                                   &
                               0 )
       
  end do
   
  val = bvalue(apr_ty(li_lefty-ky+1:li_lefty+ky), &
               lpr_coef,                                &
               ky,                                   &
               ky,                                   &
               ar_y,                                    &
               0 )

end subroutine bvalue2d

function dvalue2d( ar_x,      &
                   ar_y,      &
                   nx,     &
                   kx,     &
                   ny,     &
                   ky,     &
                   bcoef, &
                   apr_tx,    &
                   apr_ty,    &
                   deriv1,    &
                   deriv2 ) result(res)

  sll_real64 :: ar_x
  sll_real64 :: ar_y
  sll_real64 :: res
  sll_int32  :: nx
  sll_int32  :: kx
  sll_int32  :: ny
  sll_int32  :: ky
  sll_int32  :: deriv1
  sll_int32  :: deriv2

  sll_real64, dimension(:)   :: apr_tx    ! nx + kx
  sll_real64, dimension(:)   :: apr_ty    ! ny + ky
  sll_real64, dimension(:,:) :: bcoef !(nx,ny)

  sll_int32  :: li_j
  sll_int32  :: li_mflag
  sll_int32  :: li_lefty

  sll_real64, dimension (1:ky),target:: lpr_coef ! ky
    
  call interv( apr_ty, ny + ky, ar_y, li_lefty, li_mflag)
    
  if ( li_mflag .NE. 0 ) then
    res = 0.0_8
    return 
  end if
    
  do li_j = 1, ky
       
    lpr_coef(li_j) = bvalue( apr_tx,                       &
                             bcoef(1:nx,li_lefty-ky+li_j), &
                             nx,                           &
                             kx,                           &
                             ar_x,                         &
                             deriv1 )
       
  end do

  res = bvalue( apr_ty(li_lefty-ky+1:li_lefty+ky), &
                lpr_coef,                          &
                ky,                                &
                ky,                                &
                ar_y,                              &
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
subroutine spli2d ( tau, gtau, t, n, k, m, bcoef, iflag )
    
  sll_real64, dimension(:),   intent(in)    :: tau   !(n)
  sll_real64, dimension(:,:), intent(in)    :: gtau  !(n,m)
  sll_real64, dimension(:),   intent(in)    :: t     !(n+k)
  sll_int32,                  intent(in)    :: n
  sll_int32,                  intent(in)    :: k
  sll_int32,                  intent(in)    :: m

  sll_real64, dimension(:,:), intent(inout) :: bcoef !(m,n)
  sll_int32,                  intent(out)   :: iflag

  sll_real64, dimension(n)         :: work  !(n)
  sll_real64, dimension((2*k-1)*n) :: q     !((2*k-1)*n)

  sll_int32                :: i
  sll_int32                :: ilp1mx
  sll_int32                :: j
  sll_int32                :: jj
  sll_int32                :: left


  sll_real64:: taui
    
  left = k
   
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
      q(jj) = work(j)
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

subroutine spli2d_custom ( nx,     &
                           kx,     &
                           taux,  &
                           ny,     &
                           ky,     &
                           tauy,  &
                           apr_g,     &
                           bcoef, &
                           apr_tx,    &
                           apr_ty )

  sll_int32,                 intent(in)    :: nx
  sll_int32,                 intent(in)    :: kx
  sll_real64, dimension(:),  intent(in)    :: taux  !nx
  sll_int32,                 intent(in)    :: ny
  sll_int32,                 intent(in)    :: ky
  sll_real64, dimension(:),  intent(in)    :: tauy  !ny	
  sll_real64, dimension(:,:),intent(inout) :: apr_g     !nx,ny	
  sll_real64, dimension(:,:),intent(inout) :: bcoef !nx , ny 
  sll_real64, dimension(:)  ,intent(inout) :: apr_tx    !nx + kx
  sll_real64, dimension(:)  ,intent(inout) :: apr_ty    !ny + ky 

  sll_real64, dimension(nx,ny) :: lpr_work1
  sll_real64, dimension(nx*ny) :: lpr_work3
  sll_real64, dimension(nx*(2*kx-1) ) :: q
  sll_real64, dimension(( 2*ky-1) * ny ) :: lpr_work32
  sll_real64, dimension( ny         ) :: lpr_work4
  sll_real64, dimension(1:ny,1:nx),target :: lpr_work5 !  ny , nx 
  sll_real64, dimension(:,:),pointer :: lpr_work5_ptr !  ny , nx 
  sll_real64, dimension(1:ny),target:: apr_ty_bis
  sll_real64, dimension(:),pointer:: apr_ty_bis_ptr

  sll_int32  :: i, j, iflag
    
  lpr_work1(:,:) = 0.0
    
  ! *** set up knots
  !     interpolate between knots
    
  apr_tx ( 1 : kx ) = taux ( 1 )
  apr_tx ( nx + 1 : nx + kx ) = taux ( nx )
    
  if ( mod(kx,2) == 0 ) then
    do i = kx + 1, nx
      apr_tx ( i ) = taux ( i - kx/2 ) 
    end do
  else
    do i = kx + 1, nx
      apr_tx ( i ) = 0.5*( taux ( i - (kx-1)/2 ) + &
                              taux ( i -1 - (kx-1)/2 ) )
    end do
  end if
  !  *** construct b-coefficients of interpolant
  apr_ty = 0.0_f64
  if ( mod(ky,2) == 0 ) then
    do i = ky + 1, ny
      apr_ty ( i ) = tauy ( i - ky/2 ) 
    end do
  else
    do i = ky + 1, ny
      apr_ty ( i ) = 0.5*( tauy ( i - (ky-1)/2 ) + &
                              tauy ( i -1 - (ky-1)/2 ) )
    end do
  end if
  apr_ty ( 1 : ky ) = tauy ( 1 )
  apr_ty ( ny + 1 : ny + ky ) = tauy ( ny )
  apr_ty_bis = tauy(1:ny)
  lpr_work5_ptr => lpr_work5

  bcoef(1:nx,1:ny) = apr_g

  call spli2d ( taux,      &
                bcoef,     &
                apr_tx,        &
                nx,         &
                kx,         &
                ny,         &
                apr_g,    &
                iflag )
    
  bcoef(:,:) = 0.0_8
  lpr_work4  = 0.0_8
  lpr_work3  = 0.0_8
  lpr_work32 = 0.0_8
     
  apr_ty_bis_ptr => apr_ty_bis
  call spli2d ( apr_ty_bis_ptr, &
                apr_g,  &
                apr_ty,         &
                ny,          &
                ky,          &
                nx,          &
                bcoef,      &
                iflag )
     
end subroutine spli2d_custom

subroutine spli2d_custom_derder ( nx,   &
                                  nx_der,&
                                  kx,&
                                  taux,&
                                  taux_der,&
                                  ny,&
                                  ny_der,&
                                  ky,&
                                  tauy,&
                                  tauy_der,&
                                  apr_g,&
                                  apr_g_der1,&
                                  apr_g_der2,&
                                  bcoef,&
                                  apr_tx,&
                                  apr_ty )

  sll_int32  :: nx, kx, ny, ky
  sll_int32  :: nx_der,ny_der
  sll_real64, dimension(:),pointer :: taux !!nx
  sll_real64, dimension(:),pointer :: tauy !! ny
  sll_int32,  dimension(:),pointer :: taux_der !!nx_der
  sll_int32,  dimension(:),pointer :: tauy_der !!ny_der
  sll_real64, dimension(:,:),pointer :: apr_g    ! nx,ny
  sll_real64, dimension(:,:),pointer :: apr_g_der1 ! nx_der,ny
  sll_real64, dimension(:,:),pointer :: apr_g_der2 !ny_der,nx + nx_der

  sll_real64, dimension(:,:),pointer::bcoef!nx + nx_der,ny+ ny_der 
  sll_real64, dimension( : ),pointer:: apr_tx ! nx + kx + nx_der
  sll_real64, dimension( : ),pointer:: apr_ty ! ny + ky + ny_der

  sll_real64, dimension ( nx + nx_der , ny + ny_der) :: lpr_work1
  sll_real64, dimension ( nx + nx_der ) :: lpr_work2
  sll_real64, dimension ( (nx + nx_der)* (ny+ny_der) ) :: lpr_work3
  sll_real64, dimension ( (nx+nx_der) *( 2*kx-1) ) :: lpr_work31
  sll_real64, dimension (( 2*ky-1) * (ny+ny_der) ) :: lpr_work32
  sll_real64, dimension ( ny +ny_der) :: lpr_work4
  sll_real64, dimension (1:ny,1:nx+nx_der),target:: lpr_work5 !  ny , nx
  sll_real64, dimension (:,:),pointer :: lpr_work5_ptr 
  sll_int32  :: li_iflag
   
  lpr_work1(:,:) = 0.0
  lpr_work5(:,:) = 0.0
    
  ! *** set up knots
  !     interpolate between knots
    
  apr_tx = 0.0_f64
  apr_tx ( 1 : kx ) = taux ( 1 )
  apr_tx ( nx+ nx_der + 1: nx + nx_der + kx ) = taux ( nx )
    
  if (nx + nx_der + kx == nx + 2*(kx-1)) then
    apr_tx (kx+1: nx+ nx_der) = taux(2:nx-1)
  else
    print*, 'problem with construction of knots' 
  end if
    
  !  *** construct b-coefficients of interpolant
  apr_ty = 0.0_f64
  apr_ty ( 1 : ky ) = tauy ( 1 )
  apr_ty ( ny+ ny_der + 1: ny + ny_der + ky ) = tauy ( ny )
    
  if (ny + ny_der + ky == ny + 2*(ky-1)) then
     apr_ty (ky+1: ny+ ny_der) = tauy(2:ny-1)
  else
     print*, 'problem with construction of knots' 
  end if
    
  lpr_work5_ptr => lpr_work5

  call spli2d_der ( taux,      &
                    apr_g,         &
                    taux_der,  &
                    apr_g_der1,    &
                    apr_tx,        &
                    nx,         &
                    nx_der,     &
                    kx,         &
                    ny,         &
                    lpr_work2,     &
                    lpr_work31,    &
                    lpr_work5_ptr, &
                    li_iflag )
    
   bcoef(:,:) =0.0_8
   lpr_work4 = 0.0_8
   lpr_work3 = 0.0_8
   lpr_work32= 0.0_8
     
   call spli2d_der ( tauy,        &
                     lpr_work5_ptr,   &
                     tauy_der,    &
                     apr_g_der2,      &
                     apr_ty,          &
                     ny,           &
                     ny_der,       &
                     ky,           &
                     nx+nx_der, &
                     lpr_work4,       &
                     lpr_work32,      &
                     bcoef,       &
                     li_iflag )

end subroutine spli2d_custom_derder

subroutine spli2d_perdir ( ar_L,      &
                           nx,     &
                           kx,     &
                           taux,  &
                           ny,     &
                           ky,     &
                           tauy,  &
                           apr_g,     &
                           bcoef, &
                           apr_tx,    &
                           apr_ty )

  ! CALLED WHEN WE WANT TO INTERPOL WITH A PERIODIC FIRST PARAM WITH A PERIOD = ar_L
  sll_real64 :: ar_L 
  sll_int32  :: nx, kx, ny, ky
  sll_real64, dimension ( :),pointer :: taux ! nx- 1
  sll_real64, dimension ( :),pointer :: tauy ! ny		
  sll_real64, dimension ( :,:) :: apr_g !nx - 1, ny

  sll_real64, dimension (:,:),pointer :: bcoef !  nx , ny	
  sll_real64, dimension (:),pointer :: apr_tx !  nx + kx
  sll_real64, dimension (:),pointer :: apr_ty ! ny + ky

  sll_real64, dimension (1:nx),target :: lpr_taux !  nx
  sll_real64, dimension (:),pointer :: lpr_taux_ptr
  sll_real64, dimension (1:nx,1:ny),target :: lpr_g !  nx ,ny
  sll_real64, dimension (:,:),pointer :: lpr_g_ptr

  if ( ar_L == 0.0_8 ) then
    print*,'Error spli2d_per : called with a period = 0 '
    stop
  end if
    
  lpr_taux ( 1 : nx - 1 ) = taux ( 1 : nx-1)
  lpr_taux ( nx ) = taux ( 1 ) + ar_L

  lpr_g ( 1 : nx - 1 , 1 : ny ) = apr_g ( 1 : nx - 1 , 1 : ny )
  lpr_g ( nx , 1 : ny ) = apr_g ( 1 , 1 : ny )

  lpr_taux_ptr => lpr_taux
  lpr_g_ptr => lpr_g
     
  call spli2d_custom ( nx,        & 
                       kx,        &
                       lpr_taux_ptr, &
                       ny,        &
                       ky,        &
                       tauy,     &
                       lpr_g_ptr,    &
                       bcoef,    &
                       apr_tx,       &
                       apr_ty )

     
end subroutine spli2d_perdir

subroutine spli2d_dirper (nx,      &
                          kx,      &
                          taux,   &
                          ar_L,       &
                          ny,      &
                          ky,      &  
                          tauy,   &
                          apr_g,      &
                          bcoef,  &
                          apr_tx,     &
                          apr_ty )

     ! CALLED WHEN WE WANT TO INTERPOL WITH A PERIODIC second PARAM WITH A PERIOD = ar_L

  sll_real64 :: ar_L
  sll_int32  :: nx, kx, ny, ky
  sll_real64, dimension ( :),pointer :: taux ! nx
  sll_real64, dimension (:),pointer :: tauy !  ny -1
  sll_real64, dimension ( :,:) :: apr_g ! nx , ny-1

  sll_real64, dimension (:,:),pointer :: bcoef !  nx , ny
  sll_real64, dimension ( :),pointer :: apr_tx ! nx + kx	
  sll_real64, dimension (:),pointer :: apr_ty ! ny + ky 

  sll_real64, dimension (1:ny),target :: lpr_tauy ! ny	
  sll_real64, dimension (1:nx,1:ny),target :: lpr_g  !  nx ,ny
  sll_real64, dimension (:),pointer :: lpr_tauy_ptr ! ny	
  sll_real64, dimension (:,:),pointer :: lpr_g_ptr
     
  if ( ar_L == 0.0_8 ) then
    print*,'Error spli2d_per : called with a period = 0 '
    stop
  end if
     
  lpr_tauy ( 1 : ny - 1 ) = tauy ( 1 : ny - 1 )
  lpr_tauy ( ny ) = tauy ( 1 ) + ar_L
     
  lpr_g ( 1 : nx , 1 : ny -1 ) = apr_g ( 1 : nx , 1 : ny -1)
  lpr_g (1: nx , ny ) = apr_g ( 1 : nx, 1 )
     
  lpr_tauy_ptr => lpr_tauy
  lpr_g_ptr => lpr_g

  call spli2d_custom (nx,        &
                      kx,        &
                      taux,     &
                      ny,        &
                      ky,        &
                      lpr_tauy_ptr, &
                      lpr_g_ptr,    &
                      bcoef,    &
                      apr_tx,       &
                      apr_ty )
  
 end subroutine spli2d_dirper
   
 subroutine spli2d_perper( Lx,     &
                           nx,     &
                           kx,     &
                           taux,  &
                           Ly,     &
                           ny,     &
                           ky,     &
                           tauy,  &
                           g,     &
                           bcoef, &
                           apr_tx,    &
                           apr_ty )

  sll_real64, intent(in) :: Lx
  sll_real64, intent(in) :: Ly
  sll_int32,  intent(in) :: nx
  sll_int32,  intent(in) :: kx
  sll_int32,  intent(in) :: ny
  sll_int32,  intent(in) :: ky

  sll_real64, dimension(:)   :: taux ! nx -1
  sll_real64, dimension(:)   :: tauy ! ny - 1
  sll_real64, dimension(:,:) :: g !  nx  - 1, ny - 1

  sll_real64, dimension(:,:) :: bcoef ! nx , ny
  sll_real64, dimension(:)   :: apr_tx !  nx + kx
  sll_real64, dimension(:)   :: apr_ty ! ny + ky

  sll_real64, dimension(1:nx),     target  :: lpr_taux ! tmp_ty
  sll_real64, dimension(1:ny),     target  :: lpr_tauy
  sll_real64, dimension(1:nx,1:ny),target  :: lpr_g !  ( nx ,ny)
  sll_real64, dimension(:),        pointer :: lpr_taux_ptr ! tmp_ty
  sll_real64, dimension(:),        pointer :: lpr_tauy_ptr
  sll_real64, dimension(:,:),      pointer :: lpr_g_ptr
     
  SLL_ASSERT( Lx /= 0.0_8 )
  SLL_ASSERT( Ly /= 0.0_8 )
  
  lpr_taux(1:nx-1) = taux(1:nx-1)
  lpr_taux(nx)     = taux(1) + Lx
     
  lpr_tauy(1:ny-1) = tauy(1:ny-1)
  lpr_tauy(ny)     = tauy(1)+Ly
     
  lpr_g(1:nx-1,1:ny-1) = g(1:nx-1,1:ny-1)
  lpr_g(    nx,1:ny-1) = g(     1,1:ny-1)
  lpr_g(1:nx-1,    ny) = g(1:nx-1,     1)
  lpr_g(    nx,    ny) = g(     1,     1)

  lpr_taux_ptr => lpr_taux
  lpr_tauy_ptr => lpr_tauy
  lpr_g_ptr    => lpr_g

  call spli2d_custom (  nx,           &
                        kx,           &
                        lpr_taux_ptr, & 
                        ny,           &
                        ky,           &
                        lpr_tauy_ptr, &
                        lpr_g_ptr,    &
                        bcoef,        &
                        apr_tx,       &
                        apr_ty )
     
end subroutine spli2d_perper

subroutine spli2d_der( tau,           &
                       gtau,          &
                       tau_der,       &
                       gtau_der,      &
                       t,             &
                       n,             &
                       np,            &
                       k,             &
                       m,             &
                       work,          &
                       q,             &
                       bcoef,         &
                       iflag )
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
  
  sll_int32 :: m
  sll_int32 :: n
  sll_int32 :: np
  sll_real64,dimension(:,:),pointer:: bcoef !(m,n+np)
  sll_real64,dimension(:,:),pointer:: gtau  !(n,m)
  sll_real64,dimension(:,:),pointer:: gtau_der!(np,n)
  sll_int32 :: iflag
  sll_int32 :: j
  sll_int32 :: k
  sll_real64,dimension((2*k-1)*n):: q!((2*k-1)*n)
  sll_real64,dimension(:),pointer:: t!(n+np+k)
  sll_real64,dimension(:),pointer:: tau!(n)
  sll_int32,dimension(:),pointer:: tau_der!np
  sll_real64,dimension(n):: work!(n)
  sll_real64,dimension(np):: work_der
  sll_real64,dimension(n+np):: work_result

  work_result = 0.0_f64

  do j = 1, m
       
    work(1:n) = gtau(:,j)
    work_der(1:np) = gtau_der(1:np,j)

    call splint_der( tau,         &
                     work,        &
                     tau_der,     &
                     work_der,    &
                     t,           &
                     n,           &
                     np,          &
                     k,           &
                     q,           &
                     work_result, &
                     iflag )
       
     bcoef(j,1:n+np) = work_result(1:n+np)
       
  end do
    
    
  return
end subroutine spli2d_der

end module sll_module_deboor_splines_2d
