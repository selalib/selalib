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
  real(8), dimension ( ai_nx- 1) :: apr_taux		
  real(8), dimension ( ai_nx- 1) :: apr_g	
  ! OUTPUT
  real(8), dimension ( ai_nx ) :: apr_Bcoef
  real(8), dimension ( ai_nx + ai_kx) :: apr_tx
  ! LOCAL VARIABLES		
  real(8), dimension ( ai_nx) :: lpr_taux		
  real(8), dimension ( ai_nx) :: lpr_g	
  
  if ( ar_L == 0.0_8 ) then
     print*,'Error spli1d_per : called with a period = 0 '
     stop
  end if
  
  !print*, 'rer'
  lpr_taux ( 1 : ai_nx - 1 ) = apr_taux ( 1 : ai_nx - 1 )
  lpr_taux ( ai_nx ) = apr_taux ( 1 ) + ar_L
  
  lpr_g ( 1 : ai_nx - 1  ) = apr_g ( 1 : ai_nx - 1 )
  lpr_g ( ai_nx) = apr_g ( 1 )		
  
  !print*,  lpr_g
  call splint ( &
       ai_nx, &
       ai_kx, &
       lpr_taux,&
       lpr_g,&
       apr_Bcoef,&
       apr_tx)
  ! print*, 'hello'
end subroutine spli1d_per
