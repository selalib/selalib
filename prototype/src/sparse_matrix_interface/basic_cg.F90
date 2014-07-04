MODULE BASIC_CG_MODULE
USE SPM_DEF
USE SPM_TOOLS

CONTAINS
!*******************************************************!
!                                                       !
!          METHODE DU GRADIENT CONJUGUE SIMPLE          !
!                                                       !
!*******************************************************!

SUBROUTINE BASIC_CG ( ai_nR, apr_a, api_indices, api_indptr, ai_maxIter, ar_eps, apr_U, apr_B, ai_niter, ar_resnorm, apr_err )
IMPLICIT NONE
INTEGER(KIND=SPM_INTS_KIND)                        :: ai_nR
REAL(KIND=SPM_COEF_KIND)   , DIMENSION(:), POINTER :: apr_a 
INTEGER(KIND=SPM_INTL_KIND), DIMENSION(:), POINTER :: api_indices
INTEGER(KIND=SPM_INTS_KIND), DIMENSION(:), POINTER :: api_indptr  
INTEGER(KIND=SPM_INTS_KIND)                        :: ai_maxIter
REAL(KIND=SPM_COEF_KIND)                           :: ar_eps
REAL(KIND=SPM_COEF_KIND), DIMENSION(:), POINTER    :: apr_B
REAL(KIND=SPM_COEF_KIND), DIMENSION(:), POINTER    :: apr_U
INTEGER(KIND=SPM_INTS_KIND)                        :: ai_nIter
REAL(KIND=SPM_COEF_KIND)                           :: ar_resnorm
REAL(KIND=SPM_COEF_KIND), DIMENSION(:), POINTER    :: apr_err
! LOCAL
REAL(KIND=SPM_COEF_KIND), DIMENSION(:), POINTER :: lpr_Ad
REAL(KIND=SPM_COEF_KIND), DIMENSION(:), POINTER :: lpr_r
REAL(KIND=SPM_COEF_KIND), DIMENSION(:), POINTER :: lpr_d
REAL(KIND=SPM_COEF_KIND), DIMENSION(:), POINTER :: lpr_Ux
REAL(KIND=SPM_COEF_KIND) :: lr_Norm2r1
REAL(KIND=SPM_COEF_KIND) :: lr_Norm2r0
REAL(KIND=SPM_COEF_KIND) :: lr_NormInfb
REAL(KIND=SPM_COEF_KIND) :: lr_NormInfr
REAL(KIND=SPM_COEF_KIND) :: lr_ps
REAL(KIND=SPM_COEF_KIND) :: lr_beta
REAL(KIND=SPM_COEF_KIND) :: lr_alpha
LOGICAL  :: ll_continue
INTEGER(KIND=SPM_INTS_KIND)  :: li_iter
INTEGER(KIND=SPM_INTS_KIND)  :: li_err

apr_err  = 0.0
ai_niter = 0

IF ( ( dabs ( MAXVAL ( apr_B ) ) < ar_eps ) .AND. ( dabs ( MINVAL ( apr_B ) ) < ar_eps ) ) then
        apr_U = 0.0
        return
END IF

ALLOCATE(lpr_Ad(ai_nR))
ALLOCATE(lpr_r(ai_nR))
ALLOCATE(lpr_d(ai_nR))
ALLOCATE(lpr_Ux(ai_nR))

!================!
! initialisation !
!================!
lpr_Ux = apr_U
li_iter = 0
CALL CSR_GEMV(0, ai_nR, apr_a, api_indices, api_indptr, lpr_Ux, lpr_Ad, li_err)

!-------------------!
! calcul des normes !
!-------------------!
lpr_r = apr_B - lpr_Ad
lr_Norm2r0  = DOT_PRODUCT( lpr_r , lpr_r )
lr_NormInfb = maxval( dabs( apr_B ) )

lpr_d = lpr_r
!================!
ll_continue=.true.
DO while(ll_continue)
li_iter = li_iter + 1
!--------------------------------------!
! calcul du ak parametre optimal local !
!--------------------------------------!

CALL CSR_GEMV(0, ai_nR, apr_a, api_indices, api_indptr, lpr_d, lpr_Ad, li_err)

lr_ps = DOT_PRODUCT( lpr_Ad , lpr_d )
lr_alpha = lr_Norm2r0 / lr_ps

!==================================================!
! calcul de l'approximation Xk+1 et du residu Rk+1 !
!==================================================!
! calcul des composantes residuelles
!-----------------------------------
lpr_r = lpr_r - lr_alpha * lpr_Ad

!----------------------------------------!
! approximations ponctuelles au rang k+1 !
!----------------------------------------!
lpr_Ux = lpr_Ux + lr_alpha * lpr_d

!-------------------------------------------------------!
! (a) extraction de la norme infinie du residu          !
!     pour le test d'arret                              !
! (b) extraction de la norme euclidienne du residu rk+1 !
!-------------------------------------------------------!
lr_NormInfr = maxval(dabs( lpr_r ))
lr_Norm2r1 = DOT_PRODUCT( lpr_r , lpr_r )

!==================================================!
! calcul de la nouvelle direction de descente dk+1 !
!==================================================!
lr_beta = lr_Norm2r1 / lr_Norm2r0
lr_Norm2r0 = lr_Norm2r1
lpr_d = lpr_r + lr_beta * lpr_d

!-------------------!
! boucle suivante ? !
!-------------------!
ll_continue=( ( lr_NormInfr / lr_NormInfb ) >= ar_eps ) .AND. ( li_iter < ai_maxIter )

apr_err (li_iter) = lr_NormInfr / lr_NormInfb

END DO

apr_U = lpr_Ux

IF ( li_iter == ai_maxIter ) then
        PRINT*,'Warning CG : li_iter == ai_maxIter'
        PRINT*,'Error after CG =',( lr_NormInfr / lr_NormInfb )
        PRINT*,'Niter = ', li_iter
END IF

ai_niter = li_iter
ar_resnorm = lr_NormInfr

DEALLOCATE(lpr_Ad)
DEALLOCATE(lpr_d)
DEALLOCATE(lpr_r)
DEALLOCATE(lpr_Ux)

END SUBROUTINE BASIC_CG

END MODULE BASIC_CG_MODULE
