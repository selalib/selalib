subroutine spli1d_per (&
     ar_L,&
     ai_nx,&
     ai_kx,&
     apr_taux,&
     apr_g,&
     apr_Bcoef,&
     apr_tx)
  ! CALLED WHEN WE WANT TO INTERPOL WITH A PERIODIC FIRST PARAM WITH A PERIOD = ar_L
  implicit none
  ! INPUT
  real(8) :: ar_L 
  integer  :: ai_nx, ai_kx
  real(8), dimension ( ai_nx) :: apr_taux		
  real(8), dimension ( ai_nx) :: apr_g	
  ! OUTPUT
  real(8), dimension ( ai_nx ) :: apr_Bcoef
  real(8), dimension ( ai_nx + ai_kx) :: apr_tx
  ! LOCAL VARIABLES		
  real(8), dimension ( ai_nx) :: lpr_taux
  real(8), dimension ( ai_nx) :: lpr_g	
  real(8), dimension ( ai_nx *( 2*ai_kx-1) ) :: lpr_work
  real(8), dimension ( (ai_nx-ai_kx)*(2*ai_kx+3)+5*ai_kx+3 ) :: scrtch
  integer :: iflag
  integer :: li_i
  if ( ar_L == 0.0_8 ) then
     print*,'Error spli1d_per : called with a period = 0 '
     stop
  end if
  
  !print*, 'rer'
  lpr_taux ( 1 : ai_nx - 1 ) = apr_taux ( 1 : ai_nx - 1 )
  lpr_taux ( ai_nx ) = apr_taux ( ai_nx )! apr_taux ( 1 ) + ar_L
  
  lpr_g ( 1 : ai_nx - 1  ) = apr_g ( 1 : ai_nx - 1 )
  lpr_g ( ai_nx) = apr_g ( ai_nx )	!apr_g ( 1 )		
  

  apr_tx ( 1 : ai_kx ) = lpr_taux ( 1 )
  apr_tx ( ai_nx + 1 : ai_nx + ai_kx ) = lpr_taux ( ai_nx )
  

!!$  do li_i = ai_kx + 1, ai_nx
!!$     apr_tx ( li_i ) = lpr_taux ( 2 ) + &
!!$          (li_i-(ai_kx + 1))*&
!!$          ( lpr_taux ( ai_nx-1 ) - lpr_taux ( 2 ) ) / (ai_nx-(ai_kx + 1))
!!$        
!!$  end do


  if ( mod(ai_kx,2) == 0 ) then
     do li_i = ai_kx + 1, ai_nx
        apr_tx ( li_i ) = lpr_taux ( li_i - ai_kx/2 ) 
        
     end do
  else

      do li_i = ai_kx + 1, ai_nx
        apr_tx ( li_i ) = 0.5*( lpr_taux ( li_i - (ai_kx-1)/2 ) + lpr_taux ( li_i -1 - (ai_kx-1)/2 ) )
        
     end do
     
  end if

 ! call splopt ( apr_taux,ai_nx, ai_kx, scrtch, apr_tx, iflag )
 ! print*, 'coef', lpr_g
 ! print*, 'taux x',lpr_taux
  call splint ( &
       lpr_taux,&
       lpr_g,&
       apr_tx,&
       ai_nx, &
       ai_kx, &
       lpr_work,&
       apr_Bcoef,&
       iflag)

end subroutine spli1d_per
