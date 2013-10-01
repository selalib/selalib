real(8) function dvalue2d(&
     ar_x,&
     ar_y,&
     ai_nx,&
     ai_kx,&
     ai_ny,&
     ai_ky,&
     apr_Bcoef,&
     apr_tx,&
     apr_ty,deriv1,deriv2 )
  implicit none
  ! INPUT
  real(8) :: ar_x, ar_y
  real(8) :: bvalue
  integer  :: ai_nx, ai_kx, ai_ny, ai_ky
  integer  :: deriv1,deriv2
  real(8), dimension ( ai_nx + ai_kx ) :: apr_tx
  real(8), dimension ( ai_ny + ai_ky ) :: apr_ty
  real(8), dimension ( ai_nx , ai_ny ) :: apr_Bcoef
  ! LOCAL VARIABLES
  integer  :: li_i, li_j, li_mflag, li_lefty
  real(8), dimension ( ai_ky) :: lpr_coef			
  
  call interv ( apr_ty, ai_ny + 1, ar_y, li_lefty, li_mflag )
  
  if ( li_mflag .NE. 0 ) then
     dvalue2d = 0.0_8
     return 
  end if
  
  do li_j = 1, ai_ky
     
     lpr_coef ( li_j ) = bvalue(&
          apr_tx,&
          apr_bcoef ( 1 , li_lefty - ai_ky + li_j ),&
          ai_nx,&
          ai_kx,&
          ar_x,&
          deriv1 )
     
  end do
  !print*, lpr_coef
  dvalue2d = bvalue(&
       apr_ty ( li_lefty - ai_ky + 1 ),&
       lpr_coef,&
       ai_ky,&
       ai_ky,&
       ar_y,&
       deriv2 )
  
end function dvalue2d
