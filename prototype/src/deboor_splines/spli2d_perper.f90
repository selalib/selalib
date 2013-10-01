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
  real(8) :: ar_Lx
  real(8) :: ar_Ly
  integer  :: ai_nx, ai_kx, ai_ny, ai_ky
  real(8), dimension ( ai_nx -1) :: apr_taux
  real(8), dimension ( ai_ny - 1 )   :: apr_tauy
  real(8), dimension ( ai_nx  - 1, ai_ny - 1) :: apr_g
  ! OUTPUT
  real(8), dimension ( ai_nx , ai_ny) :: apr_Bcoef
  real(8), dimension ( ai_nx + ai_kx) :: apr_tx
  real(8), dimension ( ai_ny + ai_ky) :: apr_ty
  ! LOCAL VARIABLES
  real(8), dimension ( ai_nx) :: lpr_taux
  real(8), dimension ( ai_ny) :: lpr_tauy
  real(8), dimension ( ai_nx ,ai_ny) :: lpr_g
  
  if ( ar_Lx == 0.0_8 ) then
     print*,'Error spli2d_perper : called with a period = 0 '
     stop
  end if
  if ( ar_Ly == 0.0_8 ) then
     print*,'Error spli2d_perper : called with a period = 0 '
     stop
  end if
  
!!$  print*,'argument', ar_Lx,&
!!$     ai_nx,&
!!$     ai_kx,&
!!$     apr_taux,&
!!$     ar_Ly,&
!!$     ai_ny,&
!!$     ai_ky,&
!!$     apr_tauy
 ! print*, 'copie size ? ', size(apr_g(:,1))
 ! print*, 'copie delete ? ', apr_g(:,3)
  lpr_taux ( 1 : ai_nx - 1 ) = apr_taux ( 1 : ai_nx - 1 )
  lpr_taux ( ai_nx ) = apr_taux ( 1 ) + ar_Lx
  
 
  lpr_tauy ( 1 : ai_ny - 1 ) = apr_tauy ( 1 : ai_ny - 1 )
  lpr_tauy ( ai_ny ) = apr_tauy ( 1 ) + ar_Ly


  lpr_g ( 1 : ai_nx - 1 , 1 : ai_ny - 1 ) = &
       apr_g ( 1 : ai_nx - 1 , 1 : ai_ny -1 )
  lpr_g ( ai_nx , 1 : ai_ny -1 ) = apr_g ( 1 , 1 : ai_ny -1 )
  lpr_g ( 1 : ai_nx -1 , ai_ny ) = apr_g ( 1 : ai_nx -1, 1 )
  lpr_g ( ai_nx , ai_ny ) = apr_g ( 1 , 1 )

!!  print*, '&&&&&&&&&&&&&'
!!$  print*, apr_g(1,1 : ai_ny -1)
!!$  print*, apr_g(2,1 : ai_ny -1)
!!$  print*, '&&&&&&&&&&&&&'
  !print*, 'in peper',lpr_taux


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

  !print*, 'test',apr_tx
  !print*, 'test',apr_ty
 ! print*, apr_Bcoef
  
end subroutine spli2d_perper
