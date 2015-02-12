real(8) function dvalue1d(&
     ar_x,&
     ai_nx,&
     ai_kx,&
     apr_Bcoef,&
     apr_tx,&
     deriv1)
  implicit none
  ! INPUT
  real(8) :: ar_x
  real(8) :: bvalue
  integer  :: ai_nx, ai_kx
  integer  :: deriv1
  real(8), dimension ( ai_nx + ai_kx ) :: apr_tx
  real(8), dimension ( ai_nx ) :: apr_Bcoef			
  
 
  
  !print*, lpr_coef
  dvalue1d = bvalue(&
       apr_tx,&
       apr_Bcoef,&
       ai_nx,&
       ai_kx,&
       ar_x,&
       deriv1 )
  
end function dvalue1d
