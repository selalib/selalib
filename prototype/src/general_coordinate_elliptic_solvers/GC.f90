module GC
  use tracelog_Module
  
  implicit none 

contains



	!*******************************************************!
	!                                                       !
	!          METHODE DU GRADIENT CONJUGUE SIMPLE          !
	!                                                       !
	!*******************************************************!

subroutine Gradient_conj ( this, apr_B,apr_U, ai_maxIter, ar_eps )
    use SparseMatrix_Module
    implicit none
    type(csr_matrix) :: this
    real(8), dimension(:) :: apr_U
    real(8), dimension(:) :: apr_B
    integer  :: ai_maxIter
    real(8) :: ar_eps
    !local var
    real(8), dimension(:), pointer :: lpr_Ad
    real(8), dimension(:), pointer :: lpr_r
    real(8), dimension(:), pointer :: lpr_d
    real(8), dimension(:), pointer :: lpr_Ux
    real(8) :: lr_Norm2r1
    real(8) :: lr_Norm2r0
    real(8) :: lr_NormInfb
    real(8) :: lr_NormInfr
    real(8) :: lr_ps
    real(8) :: lr_beta
    real(8) :: lr_alpha
    logical  :: ll_continue
    integer  :: li_iter
    integer  :: li_err
    integer  :: li_flag
		
    if ( this%oi_nR /= this%oi_nC ) then
            PRINT*,'ERROR Gradient_conj: The matrix must be square'
            stop
    end if

    if ( ( dabs ( MAXVAL ( apr_B ) ) < ar_eps ) .AND. ( dabs ( MINVAL ( apr_B ) ) < ar_eps ) ) then
            apr_U = 0.0_8
            return
    end if

    allocate(lpr_Ad(this%oi_nR),stat=li_err)
    if (li_err.ne.0) li_flag=10
    allocate(lpr_r(this%oi_nR),stat=li_err)
    if (li_err.ne.0) li_flag=20
    allocate(lpr_d(this%oi_nR),stat=li_err)
    if (li_err.ne.0) li_flag=30
    allocate(lpr_Ux(this%oi_nR),stat=li_err)
    if (li_err.ne.0) li_flag=40
    !================!
    ! initialisation !
    !================!
    apr_U(:)  = 0.0_8
    lpr_Ux(:) = apr_U(:)
    li_iter = 0
    call Mult_CSR_Matrix_Vector( this , lpr_Ux , lpr_Ad )
    !-------------------!
    ! calcul des normes !
    !-------------------!
    lpr_r       = apr_B - lpr_Ad
    lr_Norm2r0  = DOT_PRODUCT( lpr_r , lpr_r )
    lr_NormInfb = maxval( dabs( apr_B ) )

    lpr_d = lpr_r
    !================!

!    print *,'%%%%'
    ll_continue=.true.
    do while(ll_continue)
            li_iter = li_iter + 1
            !--------------------------------------!
            ! calcul du ak parametre optimal local !
            !--------------------------------------!

            call Mult_CSR_Matrix_Vector( this , lpr_d , lpr_Ad )
            
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
            !print*, 'norme infr = ',lr_NormInfr, 'norme infb=',  lr_NormInfb
            
    end do
    apr_U = lpr_Ux
    
    if ( li_iter == ai_maxIter ) then
            print*,'Warning Gradient_conj : li_iter == ai_maxIter'
            print*,'Error after CG =',( lr_NormInfr / lr_NormInfb )
    end if

!VG!    call printlog("Gradient_conj : duree du calcul =",ai_dtllevel=1)
!VG!    
!VG!    call printcputime ( )
    
    deallocate(lpr_Ad)
    deallocate(lpr_d)
    deallocate(lpr_r)
    deallocate(lpr_Ux)
end subroutine Gradient_conj


subroutine Gradient_conj_perper ( this, apr_B,apr_U, ai_maxIter, ar_eps )
    use SparseMatrix_Module
    implicit none
    type(csr_matrix) :: this
    real(8), dimension(:) :: apr_U
    real(8), dimension(:) :: apr_B
    integer  :: ai_maxIter
    integer  :: nb_spline_tot
    real(8) :: ar_eps
    !local var
    real(8), dimension(:), pointer :: lpr_Ad
!    real(8), dimension(:), pointer :: masse_fois_U
    real(8), dimension(:), pointer :: lpr_r
    real(8), dimension(:), pointer :: lpr_d
    real(8), dimension(:), pointer :: lpr_Ux
    real(8) :: lr_Norm2r1
    real(8) :: lr_Norm2r0
    real(8) :: lr_NormInfb
    real(8) :: lr_NormInfr
    real(8) :: lr_ps
    real(8) :: lr_beta
    real(8) :: lr_alpha
    logical  :: ll_continue
    integer  :: li_iter
    integer  :: li_err
    integer  :: li_flag
		
    if ( this%oi_nR /= this%oi_nC ) then
            PRINT*,'ERROR Gradient_conj: The matrix must be square'
            stop
    end if

    if ( ( dabs ( MAXVAL ( apr_B ) ) < ar_eps ) .AND. ( dabs ( MINVAL ( apr_B ) ) < ar_eps ) ) then
            apr_U = 0.0_8
            return
    end if

    allocate(lpr_Ad(this%oi_nR),stat=li_err)
    if (li_err.ne.0) li_flag=10
    allocate(lpr_r(this%oi_nR),stat=li_err)
    if (li_err.ne.0) li_flag=20
    allocate(lpr_d(this%oi_nR),stat=li_err)
    if (li_err.ne.0) li_flag=30
    allocate(lpr_Ux(this%oi_nR),stat=li_err)
    if (li_err.ne.0) li_flag=40
    !================!
    ! initialisation !
    !================!
    apr_U(:)  = 0.0_8
    lpr_Ux(:) = apr_U(:)
    li_iter = 0
    nb_spline_tot = size(apr_U)
    call Mult_CSR_Matrix_Vector( this , lpr_Ux , lpr_Ad )
    !-------------------!
    ! calcul des normes !
    !-------------------!
    lpr_r       = apr_B - lpr_Ad
    lr_Norm2r0  = DOT_PRODUCT( lpr_r , lpr_r )
    lr_NormInfb = maxval( dabs( apr_B ) )

    lpr_d = lpr_r
    !================!
!    print *,'%%%%'


    ll_continue=.true.
    do while(ll_continue)
            li_iter = li_iter + 1
            !--------------------------------------!
            ! calcul du ak parametre optimal local !
            !--------------------------------------!

            call Mult_CSR_Matrix_Vector( this , lpr_d , lpr_Ad )
            
            lr_ps = DOT_PRODUCT( lpr_Ad , lpr_d )
            lr_alpha = lr_Norm2r0 / lr_ps
            
            !==================================================!
            ! calcul de l'approximation Xk+1 et du residu Rk+1 !
            !==================================================!
            ! calcul des composantes residuelles
            !-----------------------------------
            lpr_r = lpr_r - lr_alpha * lpr_Ad
            lpr_r = lpr_r - sum(lpr_r)/nb_spline_tot 
            !----------------------------------------!
            ! approximations ponctuelles au rang k+1 !
            !----------------------------------------!
            lpr_Ux = lpr_Ux + lr_alpha * lpr_d

            !print *, 'a = ', mat_moyenne % opr_a(1:mat_moyenne%opi_ia(2)-1)
            !call Mult_CSR_Matrix_Vector(mat_moyenne, lpr_Ux, masse_fois_U)
            !print*, 'retereter',  sum(lpr_Ux)/nb_spline_tot
            !lpr_Ux = lpr_Ux - sum(lpr_Ux)/nb_spline_tot
            !print*, 'retereter', sum(lpr_Ux)/nb_spline_tot
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
            !print*, 'norme infr = ',lr_NormInfr, 'norme infb=',  lr_NormInfb
            
    end do
    apr_U = lpr_Ux
    
    if ( li_iter == ai_maxIter ) then
            print*,'Warning Gradient_conj : li_iter == ai_maxIter'
            print*,'Error after CG =',( lr_NormInfr / lr_NormInfb )
    end if

    call printlog("Gradient_conj : duree du calcul =",ai_dtllevel=1)
    
    call printcputime ( )
    
    deallocate(lpr_Ad)
    deallocate(lpr_d)
    deallocate(lpr_r)
    deallocate(lpr_Ux)
end subroutine Gradient_conj_perper


subroutine Gradient_conj_variante ( this, apr_B,apr_U, ai_maxIter, ar_eps,apr_U0 )
    use SparseMatrix_Module
    implicit none
    type(csr_matrix) :: this
    real(8), dimension(:) :: apr_U
    real(8), dimension(:),optional :: apr_U0
    real(8), dimension(:) :: apr_B
    integer  :: ai_maxIter
    real(8) :: ar_eps
    !local var
    real(8), dimension(:), pointer :: lpr_Ad
    real(8), dimension(:), pointer :: lpr_r
    real(8), dimension(:), pointer :: lpr_d
    real(8), dimension(:), pointer :: lpr_Ux
    real(8) :: lr_Norm2r1
    real(8) :: lr_Norm2r0
    real(8) :: lr_NormInfb
    real(8) :: lr_NormInfr
    real(8) :: lr_ps
    real(8) :: lr_beta
    real(8) :: lr_alpha
    logical  :: ll_continue
    integer  :: li_iter
    integer  :: li_err
    integer  :: li_flag
		
    if ( this%oi_nR /= this%oi_nC ) then
            PRINT*,'ERROR Gradient_conj: The matrix must be square'
            stop
    end if

    if ( ( dabs ( MAXVAL ( apr_B ) ) < ar_eps ) .AND. ( dabs ( MINVAL ( apr_B ) ) < ar_eps ) ) then
            apr_U = 0.0_8
            return
    end if

    allocate(lpr_Ad(this%oi_nR),stat=li_err)
    if (li_err.ne.0) li_flag=10
    allocate(lpr_r(this%oi_nR),stat=li_err)
    if (li_err.ne.0) li_flag=20
    allocate(lpr_d(this%oi_nR),stat=li_err)
    if (li_err.ne.0) li_flag=30
    allocate(lpr_Ux(this%oi_nR),stat=li_err)
    if (li_err.ne.0) li_flag=40
    !================!
    ! initialisation !
    !================!
    if ( present(apr_U0)) then 
       apr_U(:)  = apr_U0(:)
    else 
       apr_U(:)  = 0.0
    end if
    lpr_Ux(:) = apr_U(:)
    li_iter = 0
    call Mult_CSR_Matrix_Vector( this , lpr_Ux , lpr_Ad )
    !-------------------!
    ! calcul des normes !
    !-------------------!
    lpr_r       = apr_B - lpr_Ad
    lr_Norm2r0  = DOT_PRODUCT( lpr_r , lpr_r )
    lr_NormInfb = maxval( dabs( apr_B ) )
    
    lpr_d = lpr_r
    !================!
    
    !    print *,'%%%%'
    ll_continue=.true.
    do while(ll_continue)
       li_iter = li_iter + 1
       !--------------------------------------!
       ! calcul du ak parametre optimal local !
       !--------------------------------------!
       
       call Mult_CSR_Matrix_Vector( this , lpr_d , lpr_Ad )
       
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
       !print*, 'norme infr = ',lr_NormInfr, 'norme infb=',  lr_NormInfb
       
    end do
    apr_U = lpr_Ux
    
    if ( li_iter == ai_maxIter ) then
       print*,'Warning Gradient_conj : li_iter == ai_maxIter'
       print*,'Error after CG =',( lr_NormInfr / lr_NormInfb )
    end if

    call printlog("Gradient_conj : duree du calcul =",ai_dtllevel=1)
    
    call printcputime ( )
    
    deallocate(lpr_Ad)
    deallocate(lpr_d)
    deallocate(lpr_r)
    deallocate(lpr_Ux)
  end subroutine Gradient_conj_variante

subroutine Gradient_conj_perper_variante ( this, apr_B,apr_U, ai_maxIter, ar_eps,apr_U0 )
    use SparseMatrix_Module
    implicit none
    type(csr_matrix) :: this
    real(8), dimension(:) :: apr_U
    real(8), dimension(:),optional :: apr_U0
    real(8), dimension(:) :: apr_B
    integer  :: ai_maxIter
    integer  :: nb_spline_tot
    real(8) :: ar_eps
    !local var
    real(8), dimension(:), pointer :: lpr_Ad
    !    real(8), dimension(:), pointer :: masse_fois_U
    real(8), dimension(:), pointer :: lpr_r
    real(8), dimension(:), pointer :: lpr_d
    real(8), dimension(:), pointer :: lpr_Ux
    real(8) :: lr_Norm2r1
    real(8) :: lr_Norm2r0
    real(8) :: lr_NormInfb
    real(8) :: lr_NormInfr
    real(8) :: lr_ps
    real(8) :: lr_beta
    real(8) :: lr_alpha
    logical  :: ll_continue
    integer  :: li_iter
    integer  :: li_err
    integer  :: li_flag
		
    if ( this%oi_nR /= this%oi_nC ) then
            PRINT*,'ERROR Gradient_conj: The matrix must be square'
            stop
    end if

    if ( ( dabs ( MAXVAL ( apr_B ) ) < ar_eps ) .AND. ( dabs ( MINVAL ( apr_B ) ) < ar_eps ) ) then
            apr_U = 0.0_8
            return
    end if

    allocate(lpr_Ad(this%oi_nR),stat=li_err)
    if (li_err.ne.0) li_flag=10
    allocate(lpr_r(this%oi_nR),stat=li_err)
    if (li_err.ne.0) li_flag=20
    allocate(lpr_d(this%oi_nR),stat=li_err)
    if (li_err.ne.0) li_flag=30
    allocate(lpr_Ux(this%oi_nR),stat=li_err)
    if (li_err.ne.0) li_flag=40
    !================!
    ! initialisation !
    !================!
    if (present(apr_U0)) then
       apr_U(:)  = apr_U0(:)
    else
       apr_U(:)  = 0.0_8
    end if
    lpr_Ux(:) = apr_U(:)
    li_iter = 0
    nb_spline_tot = size(apr_U)
    call Mult_CSR_Matrix_Vector( this , lpr_Ux , lpr_Ad )
    !-------------------!
    ! calcul des normes !
    !-------------------!
    lpr_r       = apr_B - lpr_Ad
    lr_Norm2r0  = DOT_PRODUCT( lpr_r , lpr_r )
    lr_NormInfb = maxval( dabs( apr_B ) )

    lpr_d = lpr_r
    !================!
!    print *,'%%%%'


    ll_continue=.true.
    do while(ll_continue)
            li_iter = li_iter + 1
            !--------------------------------------!
            ! calcul du ak parametre optimal local !
            !--------------------------------------!

            call Mult_CSR_Matrix_Vector( this , lpr_d , lpr_Ad )
            
            lr_ps = DOT_PRODUCT( lpr_Ad , lpr_d )
            lr_alpha = lr_Norm2r0 / lr_ps
            
            !==================================================!
            ! calcul de l'approximation Xk+1 et du residu Rk+1 !
            !==================================================!
            ! calcul des composantes residuelles
            !-----------------------------------
            lpr_r = lpr_r - lr_alpha * lpr_Ad
            lpr_r = lpr_r - sum(lpr_r)/nb_spline_tot 
            !----------------------------------------!
            ! approximations ponctuelles au rang k+1 !
            !----------------------------------------!
            lpr_Ux = lpr_Ux + lr_alpha * lpr_d

            !print *, 'a = ', mat_moyenne % opr_a(1:mat_moyenne%opi_ia(2)-1)
            !call Mult_CSR_Matrix_Vector(mat_moyenne, lpr_Ux, masse_fois_U)
            !print*, 'retereter',  sum(lpr_Ux)/nb_spline_tot
            !lpr_Ux = lpr_Ux - sum(lpr_Ux)/nb_spline_tot
            !print*, 'retereter', sum(lpr_Ux)/nb_spline_tot
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
            !print*, 'norme infr = ',lr_NormInfr, 'norme infb=',  lr_NormInfb
            
    end do
    apr_U = lpr_Ux
    
    if ( li_iter == ai_maxIter ) then
            print*,'Warning Gradient_conj : li_iter == ai_maxIter'
            print*,'Error after CG =',( lr_NormInfr / lr_NormInfb )
    end if

    call printlog("Gradient_conj : duree du calcul =",ai_dtllevel=1)
    
    call printcputime ( )
    
    deallocate(lpr_Ad)
    deallocate(lpr_d)
    deallocate(lpr_r)
    deallocate(lpr_Ux)
end subroutine Gradient_conj_perper_variante


subroutine Gradient_conj_adjusted ( this, apr_B,apr_U,Masse_tot, ai_maxIter, ar_eps )
    use SparseMatrix_Module
    implicit none
    type(csr_matrix) :: this
    real(8), dimension(:) :: apr_U
    real(8), dimension(:) :: apr_B
    integer  :: ai_maxIter
    real(8) :: ar_eps
    !local var
    real(8), dimension(:), pointer :: lpr_Ad
    real(8), dimension(:), pointer :: lpr_r
    real(8), dimension(:), pointer :: lpr_d
    real(8), dimension(:), pointer :: lpr_Ux
    real(8), dimension(:), pointer :: Masse_tot
    real(8) :: lr_Norm2r1
    real(8) :: lr_Norm2r0
    real(8) :: lr_NormInfb
    real(8) :: lr_NormInfr
    real(8) :: lr_ps
    real(8) :: lr_beta
    real(8) :: lr_alpha
    logical  :: ll_continue
    integer  :: li_iter
    integer  :: li_err
    integer  :: li_flag
		
    if ( this%oi_nR /= this%oi_nC ) then
            PRINT*,'ERROR Gradient_conj: The matrix must be square'
            stop
    end if

    if ( ( dabs ( MAXVAL ( apr_B ) ) < ar_eps ) .AND. ( dabs ( MINVAL ( apr_B ) ) < ar_eps ) ) then
            apr_U = 0.0_8
            return
    end if

    allocate(lpr_Ad(this%oi_nR),stat=li_err)
    if (li_err.ne.0) li_flag=10
    allocate(lpr_r(this%oi_nR),stat=li_err)
    if (li_err.ne.0) li_flag=20
    allocate(lpr_d(this%oi_nR),stat=li_err)
    if (li_err.ne.0) li_flag=30
    allocate(lpr_Ux(this%oi_nR),stat=li_err)
    if (li_err.ne.0) li_flag=40
    !================!
    ! initialisation !
    !================!
    apr_U(:)  = 0.0_8
    lpr_Ux(:) = apr_U(:)
    li_iter = 0
    call Mult_CSR_Matrix_Vector( this , lpr_Ux , lpr_Ad )

    lpr_Ad = lpr_Ad - dot_product(Masse_tot, lpr_Ux)
    !print*, dot_product(Masse_tot, lpr_Ux)
    !-----------------
    !-------------------!
    ! calcul des normes !
    !-------------------!
    lpr_r       = apr_B - lpr_Ad
    lr_Norm2r0  = DOT_PRODUCT( lpr_r , lpr_r )
    lr_NormInfb = maxval( dabs( apr_B ) )

    lpr_d = lpr_r
    !================!

!    print *,'%%%%'
    ll_continue=.true.
    do while(ll_continue)
            li_iter = li_iter + 1
            !--------------------------------------!
            ! calcul du ak parametre optimal local !
            !--------------------------------------!

            call Mult_CSR_Matrix_Vector( this , lpr_d , lpr_Ad )
            
            lpr_Ad = lpr_Ad - dot_product(Masse_tot, lpr_d)
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
            !lpr_d(1) =  apr_B(1)
            !-------------------!
            ! boucle suivante ? !
            !-------------------!
            ll_continue=( ( lr_NormInfr / lr_NormInfb ) >= ar_eps ) .AND. ( li_iter < ai_maxIter )
            !print*, 'norme infr = ',lr_NormInfr, 'norme infb=',  lr_NormInfb
            
    end do
    apr_U = lpr_Ux
    
    if ( li_iter == ai_maxIter ) then
            print*,'Warning Gradient_conj : li_iter == ai_maxIter'
            print*,'Error after CG =',( lr_NormInfr / lr_NormInfb )
    end if

    call printlog("Gradient_conj : duree du calcul =",ai_dtllevel=1)
    
    call printcputime ( )
    
    deallocate(lpr_Ad)
    deallocate(lpr_d)
    deallocate(lpr_r)
    deallocate(lpr_Ux)
end subroutine Gradient_conj_adjusted
end module GC

