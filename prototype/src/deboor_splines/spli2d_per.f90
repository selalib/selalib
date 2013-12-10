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
  real(8) :: ar_L 
  integer  :: ai_nx, ai_kx, ai_ny, ai_ky
  real(8), dimension ( ai_nx	- 1			) :: apr_taux
  real(8), dimension ( ai_ny				) :: apr_tauy		
  real(8), dimension ( ai_nx - 1, ai_ny	) :: apr_g	
  ! OUTPUT
  real(8), dimension ( ai_nx , ai_ny		) :: apr_Bcoef
  real(8), dimension ( ai_nx + ai_kx		) :: apr_tx
  real(8), dimension ( ai_ny + ai_ky     ) :: apr_ty
  ! LOCAL VARIABLES		
  real(8), dimension ( ai_nx				) :: lpr_taux		
  real(8), dimension ( ai_nx ,ai_ny		) :: lpr_g	
  
  if ( ar_L == 0.0_8 ) then
     print*,'Error spli2d_per : called with a period = 0 '
     stop
  end if
  
  !print*, 'rer'
  lpr_taux ( 1 : ai_nx - 1 ) = apr_taux ( 1 : ai_nx - 1 )						
  lpr_taux ( ai_nx ) = apr_taux ( 1 ) + ar_L						
  
  !print*, 'tx',lpr_taux
  !print*, 'ty',apr_tauy
  lpr_g ( 1 : ai_nx - 1 , 1 : ai_ny ) = apr_g ( 1 : ai_nx - 1 , 1 : ai_ny )
  lpr_g ( ai_nx , 1 : ai_ny ) = apr_g ( 1 , 1 : ai_ny )		
  
  !print*,  lpr_g
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
  real(8) :: ar_L
  integer  :: ai_nx, ai_kx, ai_ny, ai_ky
  real(8), dimension ( ai_nx	             	) :: apr_taux
  real(8), dimension ( ai_ny -1				) :: apr_tauy
  real(8), dimension ( ai_nx , ai_ny-1     	) :: apr_g
  ! OUTPUT
  real(8), dimension ( ai_nx , ai_ny		) :: apr_Bcoef
  real(8), dimension ( ai_nx + ai_kx		) :: apr_tx
  real(8), dimension ( ai_ny + ai_ky     ) :: apr_ty
  ! LOCAL VARIABLES
  real(8), dimension ( ai_nx				) :: lpr_tauy
  real(8), dimension ( ai_nx ,ai_ny		) :: lpr_g
  
  if ( ar_L == 0.0_8 ) then
     print*,'Error spli2d_per : called with a period = 0 '
     stop
  end if
  
  !print*, 'rer'
  lpr_tauy ( 1 : ai_ny - 1 ) = apr_tauy ( 1 : ai_ny - 1 )
  lpr_tauy ( ai_ny ) = apr_tauy ( 1 ) + ar_L
  
  lpr_g ( 1 : ai_nx , 1 : ai_ny -1 ) = apr_g ( 1 : ai_nx , 1 : ai_ny -1)
  lpr_g (1: ai_nx , ai_ny ) = apr_g ( 1 : ai_nx, 1 )
  !print*, 'rer1'
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
  
end subroutine spli2d_dirper
