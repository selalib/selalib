subroutine spli1d_dir ( &
     ai_nx,&
     ai_kx,&
     apr_taux,&
     apr_g,&
     apr_Bcoef,&
     apr_tx )
  implicit none
  ! INPUT
  integer  :: ai_nx, ai_kx
  real(8), dimension ( ai_nx) :: apr_taux		
  real(8), dimension ( ai_nx) :: apr_g	
  ! OUTPUT
  real(8), dimension ( ai_nx  ) :: apr_Bcoef
  real(8), dimension ( ai_nx + ai_kx ) :: apr_tx
  ! LOCAL VARIABLES		
  real(8), dimension ( ai_nx ) :: lpr_work1
  real(8), dimension ( ai_nx         ) :: lpr_work2
  real(8), dimension ( ai_nx *( 2*ai_kx-1) ) :: lpr_work31
  real(8), dimension ( (ai_nx-ai_kx)*(2*ai_kx+3)+5*ai_kx+3 ) :: scrtch
  real(8), dimension ( ai_nx + ai_kx ) :: t 
  integer  :: li_i, li_j, li_iflag,iflag,iflag1
  
  lpr_work1(:) = 0.0
  
  ! *** set up knots
  !     interpolate between knots
  ! x
  !  if (ai_kx <= 2) then 
  apr_tx ( 1 : ai_kx ) = apr_taux ( 1 )
  apr_tx ( ai_nx + 1 : ai_nx + ai_kx ) = apr_taux ( ai_nx )
  
!!$  do li_i = ai_kx + 1, ai_nx
!!$     apr_tx ( li_i ) = apr_taux ( 2 ) + &
!!$          (li_i-(ai_kx + 1))*&
!!$          ( apr_taux ( ai_nx-1 ) - apr_taux ( 2 ) ) / (ai_nx-(ai_kx + 1))
!!$     
!!$  end do

  if ( mod(ai_kx,2) == 0 ) then
     do li_i = ai_kx + 1, ai_nx
        apr_tx ( li_i ) = apr_taux ( li_i - ai_kx/2 ) 
        
     end do
  else
     
     do li_i = ai_kx + 1, ai_nx
        apr_tx ( li_i ) = 0.5*( apr_taux ( li_i - (ai_kx-1)/2 ) + apr_taux ( li_i -1 - (ai_kx-1)/2 ) )
        
     end do
     
  end if
  ! else
  !   scrtch(:) = 0.0_8
  !  print*, 'TEST SPLOT', size(apr_taux), ai_nx,ai_kx,size(apr_tx),size(scrtch)
  !   call splopt ( apr_taux,ai_nx, ai_kx, scrtch, apr_tx, iflag )
    ! print*, 'result', apr_tx
  
  !end if
 
 ! print*, 'knot',apr_tx
  
  apr_Bcoef = 0.0_8
  do li_i = 1, ai_nx
     apr_Bcoef ( li_i ) = apr_g ( li_i )
  end do
  !print*, 'coef',apr_Bcoef
     !  *** construct b-coefficients of interpolant
  !
  call splint ( &
       apr_taux,&
       apr_g,&
       apr_tx, &
       ai_nx,&
       ai_kx, &
       lpr_work31,&
       apr_Bcoef, &
       li_iflag )
  
!!$  print*, 'ok'
end subroutine spli1d_dir
