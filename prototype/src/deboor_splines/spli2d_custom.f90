subroutine spli2d_custom ( ai_nx, ai_kx, apr_taux, ai_ny, ai_ky, apr_tauy, apr_g, apr_Bcoef, apr_tx, apr_ty )
	implicit none
		! INPUT
		integer  :: ai_nx, ai_kx, ai_ny, ai_ky
		real(8), dimension ( ai_nx			) :: apr_taux
		real(8), dimension ( ai_ny			) :: apr_tauy		
		real(8), dimension ( ai_nx , ai_ny ) :: apr_g	
		! OUTPUT
		real(8), dimension ( ai_nx , ai_ny ) :: apr_Bcoef
		real(8), dimension ( ai_nx + ai_kx ) :: apr_tx
		real(8), dimension ( ai_ny + ai_ky ) :: apr_ty
		! LOCAL VARIABLES		
		real(8), dimension ( ai_nx , ai_ny ) :: lpr_work1
		real(8), dimension ( ai_nx         ) :: lpr_work2
		real(8), dimension ( ai_nx * ai_ny ) :: lpr_work3
		integer  :: li_i, li_j, li_iflag
									
									
! *** set up knots
!     interpolate between knots
		! x
		apr_tx ( 1 : ai_kx ) = apr_taux ( 1 )		
		apr_tx ( ai_nx + 1 : ai_nx + ai_kx ) = apr_taux ( ai_nx )		
		
		do li_i = ai_kx + 1, ai_nx
			apr_tx ( li_i ) = apr_taux ( 2 ) + (li_i-(ai_kx + 1))*( apr_taux ( ai_nx-1 ) - apr_taux ( 2 ) ) / (ai_nx-(ai_kx + 1))
 !( apr_taux ( li_i - ai_kx + 1 ) + apr_taux ( li_i - ai_kx + 2 ) ) / 2.0_8!
		end do

		apr_ty ( 1 : ai_ky ) = apr_tauy ( 1 )
		apr_ty ( ai_ny + 1 : ai_ny + ai_ky ) = apr_tauy ( ai_ny )		
		
		do li_j = ai_ky + 1, ai_ny
			apr_ty ( li_j ) = apr_tauy ( 2 ) + (li_j-(ai_ky + 1))*( apr_tauy ( ai_ny -1) - apr_tauy ( 2 ) ) / (ai_ny-(ai_ky + 1))
!( apr_tauy ( li_j - ai_ky + 1 ) + apr_tauy ( li_j - ai_ky + 2 ) ) / 2.0_8!
		end do

		do li_i = 1, ai_nx
			do li_j = 1, ai_ny
				apr_Bcoef ( li_i, li_j ) = apr_g ( li_i, li_j )
			end do
		end do
        !print*, 'ert',apr_tx
!  *** construct b-coefficients of interpolant
!
		call spli2d ( apr_taux, apr_bcoef, apr_tx, ai_nx, ai_kx, ai_ny, lpr_work2, lpr_work3, lpr_work1, li_iflag )
		call spli2d ( apr_tauy, lpr_work1, apr_ty, ai_ny, ai_ky, ai_nx, lpr_work2, lpr_work3, apr_bcoef, li_iflag )

	end subroutine spli2d_custom
