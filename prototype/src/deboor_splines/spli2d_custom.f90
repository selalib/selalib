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
  integer  :: ai_nx, ai_kx, ai_ny, ai_ky
  real(8), dimension ( ai_nx) :: apr_taux
  real(8), dimension ( ai_ny) :: apr_tauy		
  real(8), dimension ( ai_nx,ai_ny) :: apr_g	
  ! OUTPUT
  real(8), dimension ( ai_nx , ai_ny ) :: apr_Bcoef
  real(8), dimension ( ai_nx + ai_kx ) :: apr_tx
  real(8), dimension ( ai_ny + ai_ky ) :: apr_ty
  ! LOCAL VARIABLES		
  real(8), dimension ( ai_nx , ai_ny ) :: lpr_work1
  real(8), dimension ( ai_nx         ) :: lpr_work2
  real(8), dimension ( ai_nx * ai_ny ) :: lpr_work3
  real(8), dimension ( ai_nx *( 2*ai_kx-1) ) :: lpr_work31
  real(8), dimension (( 2*ai_ky-1) * ai_ny ) :: lpr_work32
  real(8), dimension ( ai_ny         ) :: lpr_work4
  real(8), dimension ( ai_ny , ai_nx ) :: lpr_work5
  real(8), dimension ( (ai_nx-ai_kx)*(2*ai_kx+3)+5*ai_kx+3 ) :: scrtch
  real(8), dimension ( (ai_ny-ai_ky)*(2*ai_ky+3)+5*ai_ky+3 ) :: scrtch1
  real(8), dimension ( ai_nx + ai_kx ) :: t 
  integer  :: li_i, li_j, li_iflag,iflag,iflag1
  integer :: o
  

  lpr_work1(:,:) = 0.0

  ! *** set up knots
  !     interpolate between knots
  ! x
 ! if (ai_kx <= 4) then 
     apr_tx ( 1 : ai_kx ) = apr_taux ( 1 )
     apr_tx ( ai_nx + 1 : ai_nx + ai_kx ) = apr_taux ( ai_nx )
     
     do li_i = ai_kx + 1, ai_nx
        apr_tx ( li_i ) = apr_taux ( 2 ) + &
             (li_i-(ai_kx + 1))*&
             ( apr_taux ( ai_nx-1 ) - apr_taux ( 2 ) ) / (ai_nx-(ai_kx + 1))
        
     end do

  apr_Bcoef = 0.0_8
  do li_i = 1, ai_nx
     do li_j = 1, ai_ny
        apr_Bcoef ( li_i, li_j ) = apr_g ( li_i, li_j )
     end do
  end do
  !print*, 'coef',apr_Bcoef
  !  *** construct b-coefficients of interpolant
  !
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

 ! if (ai_ky <= 4) then 
     apr_ty ( 1 : ai_ky ) = apr_tauy ( 1 )
     apr_ty ( ai_ny + 1 : ai_ny + ai_ky ) = apr_tauy ( ai_ny )		
     

     do li_j = ai_ky + 1, ai_ny
        apr_ty ( li_j ) = apr_tauy ( 2 ) +&
             (li_j-(ai_ky + 1))*&
             ( apr_tauy ( ai_ny -1) - apr_tauy ( 2 ) ) / (ai_ny-(ai_ky + 1))
     end do
 ! else 
  !   call splopt ( apr_tauy,ai_ny, ai_ky, scrtch1, apr_ty, iflag1 )
 ! end if

 
  apr_bcoef(:,:) =0.0_8
  lpr_work4 = 0.0_8
  lpr_work3 = 0.0_8
  call spli2d ( &
       apr_tauy,&
       lpr_work5,&
       apr_ty,&
       ai_ny, &
       ai_ky, &
       ai_nx, &
       lpr_work4, &
       lpr_work32,&
       apr_bcoef, &
       li_iflag )


end subroutine spli2d_custom
