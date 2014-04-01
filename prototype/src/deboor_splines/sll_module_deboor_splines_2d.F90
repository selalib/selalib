module sll_module_deboor_splines_2d

#include "sll_memory.h"
#include "sll_working_precision.h"
use sll_module_deboor_splines_1d
  
  implicit none 
  
  
contains
  
  
  
  
  
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
       apr_ty,&
       val )
    implicit none
    ! INPUT
    sll_real64 :: ar_x, ar_y
    sll_real64 ::val
    sll_int32  :: ai_nx, ai_kx, ai_ny, ai_ky
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
   
    call interv ( apr_ty,ai_ny + ai_ky, ar_y, li_lefty, li_mflag )

    if ( li_mflag .NE. 0 ) then
       val = 0.0_8
       SLL_DEALLOCATE(lpr_coef,ierr)
       SLL_DEALLOCATE(tmp_tab,ierr)
       SLL_DEALLOCATE(tmp_ty,ierr)
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
    sll_real64 :: res
    sll_int32  :: ai_nx, ai_kx, ai_ny, ai_ky
    sll_int32  :: deriv1,deriv2
    sll_real64, dimension ( : ),pointer :: apr_tx ! ai_nx + ai_kx
    sll_real64, dimension ( : ),pointer :: apr_ty ! ai_ny + ai_ky
    sll_real64, dimension ( : , : ),pointer :: apr_Bcoef !(ai_nx,ai_ny)
    ! LOCAL VARIABLES
    sll_int32  :: li_i, li_j, li_mflag, li_lefty
    sll_real64, dimension ( :),pointer :: lpr_coef ! ai_ky
    sll_real64, dimension ( :),pointer :: tmp_coef,tmp_ty
    sll_int32:: ierr
    
    
    SLL_ALLOCATE(lpr_coef(ai_ky),ierr)
    SLL_ALLOCATE(tmp_coef(ai_nx),ierr)
    SLL_ALLOCATE(tmp_ty(2*ai_ky),ierr)
    
    call interv ( apr_ty, ai_ny + ai_ky, ar_y, li_lefty, li_mflag )
    
    if ( li_mflag .NE. 0 ) then
       res = 0.0_8
       SLL_DEALLOCATE(lpr_coef,ierr)
       SLL_DEALLOCATE(tmp_coef,ierr)
       SLL_DEALLOCATE(tmp_ty,ierr)
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

    SLL_DEALLOCATE(lpr_coef,ierr)
    SLL_DEALLOCATE(tmp_coef,ierr)
    SLL_DEALLOCATE(tmp_ty,ierr)
    
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
    sll_real64,dimension(:),pointer:: t!(n+k)
    sll_real64,dimension(:),pointer:: tau!(n)
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
    sll_real64, dimension (:,:),pointer :: lpr_work5 !  ai_ny , ai_nx 
    sll_real64, dimension ( (ai_nx-ai_kx)*(2*ai_kx+3)+5*ai_kx+3 ) :: scrtch
    sll_real64, dimension ( (ai_ny-ai_ky)*(2*ai_ky+3)+5*ai_ky+3 ) :: scrtch1
    sll_real64, dimension ( ai_nx + ai_kx ) :: t 
    sll_real64, dimension ( : ),pointer :: apr_ty_bis
    sll_int32  :: li_i, li_j, li_iflag,iflag,iflag1
    sll_int32 :: o,ierr
    
    SLL_ALLOCATE(apr_ty_bis( ai_ny),ierr)
    SLL_ALLOCATE(lpr_work5(ai_ny , ai_nx ),ierr)
    
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
                0.5*( apr_taux ( li_i - (ai_kx-1)/2 ) + &
                apr_taux ( li_i -1 - (ai_kx-1)/2 ) )
           
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
                0.5*( apr_tauy ( li_i - (ai_ky-1)/2 ) + &
                apr_tauy ( li_i -1 - (ai_ky-1)/2 ) )
           
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
     
     SLL_DEALLOCATE(apr_ty_bis,ierr)
     SLL_DEALLOCATE(lpr_work5,ierr)
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
     sll_real64, dimension ( :),pointer :: apr_tauy ! ai_ny		
     sll_real64, dimension ( :,:) :: apr_g	!ai_nx - 1, ai_ny
     ! OUTPUT
     sll_real64, dimension (:,:),pointer :: apr_Bcoef !  ai_nx , ai_ny	
     sll_real64, dimension (:),pointer :: apr_tx !  ai_nx + ai_kx
     sll_real64, dimension (:),pointer :: apr_ty ! ai_ny + ai_ky
     ! LOCAL VARIABLES		
     sll_real64, dimension (:),pointer :: lpr_taux !  ai_nx		
     sll_real64, dimension (:,:),pointer :: lpr_g !  ai_nx ,ai_ny
     sll_int32 :: ierr
     
     SLL_ALLOCATE(lpr_taux(ai_nx), ierr)
     SLL_ALLOCATE(lpr_g(ai_nx ,ai_ny),ierr)
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

     SLL_DEALLOCATE(lpr_taux, ierr)
     SLL_DEALLOCATE(lpr_g,ierr)
     
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
     sll_real64, dimension (:,:),pointer :: apr_Bcoef !  ai_nx , ai_ny
     sll_real64, dimension ( :),pointer :: apr_tx ! ai_nx + ai_kx	
     sll_real64, dimension (:),pointer :: apr_ty ! ai_ny + ai_ky 
     ! LOCAL VARIABLES
     sll_real64, dimension ( :),pointer :: lpr_tauy ! ai_ny	
     sll_real64, dimension (:,:),pointer :: lpr_g  !  ai_nx ,ai_ny
     sll_int32 :: ierr
     
     SLL_ALLOCATE(lpr_tauy(ai_ny), ierr)
     SLL_ALLOCATE(lpr_g(ai_nx ,ai_ny),ierr)
     
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

     SLL_DEALLOCATE(lpr_tauy, ierr)
     SLL_DEALLOCATE(lpr_g,ierr)
  
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
     sll_real64, dimension (:),pointer :: apr_taux ! ai_nx -1
     sll_real64, dimension (:),pointer :: apr_tauy ! ai_ny - 1
     sll_real64, dimension (:,:),pointer :: apr_g !  ai_nx  - 1, ai_ny - 1
     ! OUTPUT
     sll_real64, dimension ( :,:),pointer :: apr_Bcoef ! ai_nx , ai_ny
     sll_real64, dimension (:),pointer :: apr_tx !  ai_nx + ai_kx
     sll_real64, dimension ( :),pointer :: apr_ty ! ai_ny + ai_ky
     ! LOCAL VARIABLES
     sll_real64, dimension ( :),pointer :: lpr_taux ! tmp_ty
     sll_real64, dimension ( :),pointer :: lpr_tauy
     sll_real64, dimension(:,:),pointer :: lpr_g !  ( ai_nx ,ai_ny)
     sll_int32 :: ierr
     
     SLL_ALLOCATE(lpr_taux(ai_nx),ierr)
     SLL_ALLOCATE(lpr_tauy(ai_ny),ierr)
     SLL_ALLOCATE( lpr_g( ai_nx ,ai_ny),ierr)
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
     

     SLL_DEALLOCATE(lpr_taux,ierr)
     SLL_DEALLOCATE(lpr_tauy,ierr)
     SLL_DEALLOCATE( lpr_g,ierr)
   end subroutine spli2d_perper
   




   
   
 end module sll_module_deboor_splines_2d
