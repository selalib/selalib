!**************************************************************
!  Copyright Euratom-CEA
!  Authors : 
!     Virginie Grandgirard (virginie.grandgirard@cea.fr)
!     Chantal Passeron (chantal.passeron@cea.fr)
!     Guillaume Latu (guillaume.latu@cea.fr)
!     Xavier Garbet (xavier.garbet@cea.fr)
!     Philippe Ghendrih (philippe.ghendrih@cea.fr)
!     Yanick Sarazin (yanick.sarazin@cea.fr)
!  
!  This code GYSELA (for GYrokinetic SEmi-LAgrangian) 
!  is a 5D gyrokinetic global full-f code for simulating 
!  the plasma turbulence in a tokamak.
!  
!  This software is governed by the CeCILL-B license 
!  under French law and abiding by the rules of distribution 
!  of free software.  You can  use, modify and redistribute 
!  the software under the terms of the CeCILL-B license as 
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info". 
!**************************************************************
      
!---------------------------------------------
! file : spline1d
! date : 15/05/2000
! cubic spline interpolation for 
! 1D no-periodic and periodic function
!---------------------------------------------
module spline1d_class
  use prec_const
  use mem_alloc_module
  use utils_module, only : locate
  
  implicit none
  
  !*** 1D natural spline type ***
  type :: nspline1d
    integer     :: n	! dimension
    real(RKIND) :: h	! discretis. step  
    ! 2x2 matrix coefficients
    real(RKIND) :: eta1b, eta2b, eta3b, eta4b
    ! 2x2 matrix determinant
    real(RKIND) :: det_deltab
    ! index in the array
    integer    , dimension(:)  , pointer :: ipos1d
    ! diagonal  terms
    real(RKIND), dimension(:)  , pointer :: Adiag   
    ! off-diagonal terms
    real(RKIND), dimension(:)  , pointer :: Aodiag  
    ! A^-1*gamma
    real(RKIND), dimension(:,:), pointer :: Am1gamma	
    ! spline coefficients	
    real(RKIND), dimension(:)  , pointer :: scoef
    ! temporary arrays
    real(RKIND), dimension(:)  , pointer :: xprime
    real(RKIND), dimension(:)  , pointer :: rhs
    ! array used for integration
    integer    , dimension(1:6) :: indx
    real(RKIND), dimension(1:6) :: factor_int
  end type nspline1d
  
  !*** 1D periodic spline type ***
  type :: pspline1d
    integer     :: n	! dimension
    real(RKIND) :: h	! discretis. step  
    ! extrema of parallel region
    real(RKIND) :: eta1b, eta2b, eta3b, eta4b
    ! 2x2 matrix determinant
    real(RKIND) :: det_deltab
    ! index in the array
    integer    , dimension(:)  , pointer :: ipos1d    		
    ! diagonal  terms
    real(RKIND), dimension(:)  , pointer :: Adiag   
    ! off-diagonal terms
    real(RKIND), dimension(:)  , pointer :: Aodiag  
    ! A^-1*gamma
    real(RKIND), dimension(:,:), pointer :: Am1gamma	
    ! spline coefficients	
    real(RKIND), dimension(:)  , pointer :: scoef	
    ! temporary arrays
    real(RKIND), dimension(:)  , pointer :: xprime
    real(RKIND), dimension(:)  , pointer :: rhs
    ! array used for integration
    integer    , dimension(1:6) :: indx
    real(RKIND), dimension(1:6) :: factor_int
  end type pspline1d
      
  !*** computation of the cubic splines basis ***
  interface spline_basis
    module procedure spline_basis_point, spline_basis_knot
  end interface
  
  !*** computation of the first derivate of the ***
  !***  cubic splines basis                     ***
  interface spline_basisderiv
    module procedure spline_basisderiv_point, &
      spline_basisderiv_knot
  end interface
  
  !*** computation of the second derivate of the ***
  !***  cubic splines basis                      ***
  interface spline_basis2deriv
    module procedure spline_basis2deriv_point, &
      spline_basis2deriv_knot
  end interface
  
      
  !******************************
  contains
  !******************************
      
  !*******************************************************
  ! COMMON FUNCTIONS FOR NATURAL AND PERIODIC CONDITIONS *
  !*******************************************************
      
  !--------------------------------------------- 
  ! computation of the cubic spline basis 
  !   forall point of the space
  ! (available for periodic or natural function)
  !---------------------------------------------
  subroutine spline_basis_point(xprev,xstar,xnext,h,sbase)
    real(RKIND)                 , intent(in)    :: xprev,xstar
    real(RKIND)                 , intent(in)    :: xnext,h
    real(RKIND), dimension(-1:2), intent(inout) :: sbase
  
    real(RKIND) :: coeff, d_prev, d_next
  
    coeff     = 1._RKIND/(h*h*h)
    d_prev    = xstar-xprev
    d_next    = xnext-xstar
    sbase(2)  = coeff*d_prev*d_prev*d_prev
    sbase(1)  = 1._RKIND+coeff*(3._RKIND*d_prev* &
      (h*h+h*d_prev-d_prev*d_prev))
    sbase(0)  = 1._RKIND+coeff*(3._RKIND*d_next* &
      (h*h+h*d_next-d_next*d_next))
    sbase(-1) = coeff*d_next*d_next*d_next  
  end subroutine spline_basis_point
  
  
  !--------------------------------------------- 
  ! computation of the cubic spline basis 
  !  for a knot of the mesh
  ! (available for periodic or natural function)
  !---------------------------------------------
  subroutine spline_basis_knot(num,nb_knots,sbase)
    integer                     , intent(in)    :: num
    integer                     , intent(in)    :: nb_knots
    real(RKIND), dimension(-1:2), intent(INOUT) :: sbase
  
    if (num.lt.nb_knots) then
      sbase(-1) = 1._RKIND
      sbase(0)  = 4._RKIND
      sbase(1)  = 1._RKIND
      sbase(2)  = 0._RKIND
    else
      sbase(-1) = 0._RKIND
      sbase(0)  = 1._RKIND
      sbase(1)  = 4._RKIND
      sbase(2)  = 1._RKIND
    endif
  end subroutine spline_basis_knot
  
  
  !------------------------------------------------------------
  !  computation of the second derivative of the cubic spline 
  !    basis forall point of the space
  !------------------------------------------------------------ 
  subroutine spline_basisderiv_point(xprev,xstar,xnext,h,sbase)
    real(RKIND)                 , intent(in)   :: xprev
    real(RKIND)                 , intent(in)   :: xstar
    real(RKIND)                 , intent(in)   :: xnext
    real(RKIND)                 , intent(in)   :: h
    real(RKIND), dimension(-1:2), intent(inout):: sbase
    
    real(RKIND) :: coeff, d_prev, d_next
    
    coeff     = 3._RKIND/(h*h*h)
    d_prev    = xstar-xprev
    d_next    = xnext-xstar
    sbase(2)  = coeff*d_prev*d_prev
    sbase(1)  = coeff * &
      (h*h+2._RKIND*h*d_prev-3._RKIND*d_prev*d_prev)
    sbase(0)  = -1._RKIND*coeff * &
      (h*h+2._RKIND*h*d_next-3._RKIND*d_next*d_next)
    sbase(-1) = -1._RKIND*coeff*d_next*d_next
  end subroutine spline_basisderiv_point
      
  !--------------------------------------------- 
  ! computation of the first derivate of 
  !  the cubic spline basis for a knot of the mesh
  ! (available for periodic or natural function)
  !---------------------------------------------
  subroutine spline_basisderiv_knot(num,nb_knots,h,sbase_prime)
    integer                     , intent(in)    :: num
    integer                     , intent(in)    :: nb_knots
    real(RKIND)                 , intent(in)    :: h
    real(RKIND), dimension(-1:2), intent(INOUT) :: sbase_prime
  
    if (num.lt.nb_knots) then
      sbase_prime(-1) = -3._RKIND/h
      sbase_prime(0)  = 0._RKIND
      sbase_prime(1)  = 3._RKIND/h
      sbase_prime(2)  = 0._RKIND
    else
      sbase_prime(-1) = 0._RKIND
      sbase_prime(0)  = -3._RKIND/h
      sbase_prime(1)  = 0._RKIND
      sbase_prime(2)  = 3._RKIND/h
    endif
  end subroutine spline_basisderiv_knot
      
  !------------------------------------------------------------
  !  computation of the second derivative of the cubic spline 
  !    basis forall point of the space
  !------------------------------------------------------------ 
  subroutine spline_basis2deriv_point(xprev,xstar,xnext,h, &
    sbase_second)
    real(RKIND)                 , intent(in)   :: xprev
    real(RKIND)                 , intent(in)   :: xstar
    real(RKIND)                 , intent(in)   :: xnext
    real(RKIND)                 , intent(in)   :: h
    real(RKIND), dimension(-1:2), intent(inout):: sbase_second
    
    real(RKIND) :: coeff, d_prev, d_next
    
    coeff            = 6._RKIND/(h*h*h)
    d_prev           = xstar-xprev
    d_next           = xnext-xstar
    sbase_second(2)  = coeff*d_prev
    sbase_second(1)  = coeff*(h-3._RKIND*d_prev)
    sbase_second(0)  = coeff*(h-3._RKIND*d_next)
    sbase_second(-1) = coeff*d_next
  end subroutine spline_basis2deriv_point
  
  
  !--------------------------------------------- 
  ! computation of the second derivate of the
  !  cubic spline basis for a knot of the mesh
  ! (available for periodic or natural function)
  !---------------------------------------------
  subroutine spline_basis2deriv_knot(num,nb_knots,h,sbase_second)
    integer                     , intent(in)    :: num
    integer                     , intent(in)    :: nb_knots
    real(RKIND)                 , intent(in)    :: h
    real(RKIND), dimension(-1:2), intent(INOUT) :: sbase_second
  
    real(RKIND) :: h2
      
    h2 = h*h
    if (num.lt.nb_knots) then
      sbase_second(-1) = 6._RKIND/h2
      sbase_second(0)  = -12._RKIND/h2
      sbase_second(1)  = 6._RKIND/h2
      sbase_second(2)  = 0._RKIND
    else
      sbase_second(-1) = 0._RKIND
      sbase_second(0)  = 6._RKIND/h2
      sbase_second(1)  = -12._RKIND/h2
      sbase_second(2)  = 6._RKIND/h2
    endif
  end subroutine spline_basis2deriv_knot
      
  !---------------------------------------------------
  ! computation of the coefficients (used for
  !  the integration calculation) at the boundaries, 
  !  because the spline -1,0,1,n-1,n,n+1 are not 
  !  totally in the integration domain
  !---------------------------------------------------
  subroutine integration_coef_boundaries(n,h,indx,factor)
    integer                    , intent(in) :: n
    real(RKIND)                , intent(in) :: h
    integer    , dimension(1:6), intent(out) :: indx
    real(RKIND), dimension(1:6), intent(out) :: factor     
      
    !*** array initialisation for the boundary terms ***
    indx(1) = -1   ; factor(1) = h/4._RKIND
    indx(2) = 0    ; factor(2) = 3._RKIND*h
    indx(3) = 1    ; factor(3) = (23._RKIND/4._RKIND)*h
    indx(4) = n-1  ; factor(4) = (23._RKIND/4._RKIND)*h
    indx(5) = n    ; factor(5) = 3._RKIND*h
    indx(6) = n+1  ; factor(6) = h/4._RKIND
  end subroutine integration_coef_boundaries
      
  !********************************
  ! FUNCTIONS FOR NATURAL SPLINES *
  !********************************
  !----------------------------------------------------- 
  ! 1D spline initialisation for non-periodic function,
  !  with two possibilities for the boundary conditions : 
  !   - dirichlet conditions 
  !   - information on the first derivative
  ! In the case of Dirichlet conditions, an Lagrange 
  ! interpolation is used for the both extremities
  !-----------------------------------------------------   
  subroutine new_spline1d_natural(nsthis,n,h)
    use globals, only : memory_test
    type(nspline1d), intent(out) :: nsthis
    integer        , intent(in)  :: n	! problem dimension
    real(RKIND)    , intent(in)  :: h	! discretisation step
    
    integer :: i, err
      
    nsthis%n = n
    nsthis%h = h
        
    !*** memory allocation ***
    call glob_allocate(nsthis%xprime,0,n,'xprime')
    call glob_allocate(nsthis%rhs,0,n,'rhs')
    call glob_allocate(nsthis%ipos1d,0,n,'ipos1d') 
    call glob_allocate(nsthis%Adiag,0,n,'Adiag')
    call glob_allocate(nsthis%Aodiag,0,n-1,'Aodiag')
    call glob_allocate(nsthis%Am1gamma,0,n,0,1,'Am1gamma')
    call glob_allocate(nsthis%scoef,-1,n+1,'scoef')
    
    if (.not.memory_test) then
      !*** ipos1d array initialization ***
      do i = 0,n-1
        nsthis%ipos1d(i) = i
      enddo
      nsthis%ipos1d(n) = n-1
      
      !*** A factorisation ***
      nsthis%Adiag  = 4._RKIND
      nsthis%Aodiag = 1._RKIND
      call dpttrf(n+1,nsthis%Adiag,nsthis%Aodiag,err)
      if (err.ne.0) then
        write(6,*) 'Natural case : problem in the A factorisation'
        stop
      end if
      
      !*** Computation of A-1.gamma ***
      nsthis%Am1gamma      = 0._RKIND
      nsthis%Am1gamma(0,1) = 1._RKIND 
      nsthis%Am1gamma(n,0) = 1._RKIND
      call dpttrs(n+1,2,nsthis%Adiag, &
        nsthis%Aodiag,nsthis%Am1gamma,n+1,err)
      if (err.ne.0) then
        write(6,*) 'Natural case : problem in A-1.gamma solving'
        stop
      end if
      
      !**************************************** 
      !*   2x2 matrice (deltab) computation : *
      !*    deltab = delta-lambda.A-1.gamma   *
      !*     where delta = 1  0               *
      !*                   0 -1               *
      !****************************************     
      nsthis%eta1b = 1._RKIND+nsthis%Am1gamma(n-1,0)
      nsthis%eta2b = nsthis%Am1gamma(n-1,1)
      nsthis%eta3b = -nsthis%Am1gamma(1,0)
      nsthis%eta4b = -(1._RKIND+nsthis%Am1gamma(1,1))
      
      ! Computation of deltab determinant
      nsthis%det_deltab = nsthis%eta1b*nsthis%eta4b - &
        nsthis%eta2b*nsthis%eta3b
      if (nsthis%det_deltab.eq.0._RKIND) then
        write(6,*) &
          'Natural case : problem with determinant equal to'
        write(6,*) nsthis%det_deltab
        stop
      end if
      
      ! Initialisation of the coefficients required for integration
      call integration_coef_boundaries(nsthis%n,nsthis%h, &
        nsthis%indx,nsthis%factor_int)
    end if
  end subroutine new_spline1d_natural
  
      
  !----------------------------------------------------- 
  ! 1D natural spline destructor
  !-----------------------------------------------------   
  subroutine del_spline1d_natural(nsthis)
    type(nspline1d), intent(inout) :: nsthis
    
    !*** memory allocation ***
    call glob_deallocate(nsthis%xprime)
    call glob_deallocate(nsthis%rhs)
    call glob_deallocate(nsthis%ipos1d)
    call glob_deallocate(nsthis%Adiag)
    call glob_deallocate(nsthis%Aodiag)
    call glob_deallocate(nsthis%Am1gamma)
    call glob_deallocate(nsthis%scoef)
  end subroutine del_spline1d_natural
      
  
  !----------------------------------------------------------
  ! Compute boundary condition in the case of non-periodic
  !  function, based on the approximation of the second
  !  derivative by Lagrange polynomials
  !   (used for the cubic spline computation)
  ! Input  : BC_left and BC_right 
  !   for the boundary conditions specifications
  !   = 0 if Dirichlet condition (func=0)
  !   = 1 if Neumann (func'=0)
  !   = 2 if Hermite conditions 
  !     (func'=approximation of the derivate)
  ! Output : deriv_funcBC(0:1) with
  !   - deriv_funcBC(0) = func'(n) 
  !       = 11/6h*func(n) - 
  !         3/h*func(n-1)+3/2h*func(n-2)-1/3h*func(n-3)  
  !   - deriv_funcBC(1) = func'(0)
  !       = -11/6h*func(0) + 
  !         3/h*func(1)-3/2h*func(2)+1/3h*func(3)
  !----------------------------------------------------------
  subroutine compute_natural_BC(func,n,h, &
    BC_left,BC_right,deriv_funcBC)
    use globals, only : pglobal_id
    real(RKIND), dimension(0:), intent(in)    :: func
    integer                   , intent(in)    :: n
    real(RKIND)               , intent(in)    :: h
    integer                   , intent(in)    :: BC_left, BC_right
    real(RKIND), dimension(0:), intent(inout) :: deriv_funcBC
    
    !-> computation of the approximate derivative 
    !   at the right hand side
    select case (BC_right)
    case(0)
      deriv_funcBC(0) = (-3._RKIND*func(n-1) + &
        1.5_RKIND*func(n-2) - &
        1._RKIND*func(n-3)/3._RKIND)/h
    case(1)
      deriv_funcBC(0) = 0._RKIND
    case(2)
      deriv_funcBC(0) = (11._RKIND*func(n)/6._RKIND - &
        3._RKIND*func(n-1) + 1.5_RKIND*func(n-2) - &
        1._RKIND*func(n-3)/3._RKIND)/h
    case default
      if (pglobal_id.eq.0) then
        print*, 'BC_right = ', BC_right, ' must be 0, 1 or 2'
        stop
      end if
    end select
    !-> computation of the approximate derivative 
    !   at the left hand side
    select case (BC_left)
    case(0)
      deriv_funcBC(1) = (3._RKIND*func(1) - 1.5_RKIND*func(2) + &
        1._RKIND*func(3)/3._RKIND)/h
    case(1)
      deriv_funcBC(1) = 0._RKIND
    case(2)
      deriv_funcBC(1) = (-11._RKIND*func(0)/6._RKIND + &
        3._RKIND*func(1) - 1.5_RKIND*func(2) + &
        1._RKIND*func(3)/3._RKIND)/h
    case default
      if (pglobal_id.eq.0) then
        print*, 'BC_left = ', BC_left, ' must be 0, 1 or 2'
        stop
      end if
    end select
  end subroutine compute_natural_BC
      
  !----------------------------------------------------
  ! 1D spline coefficient computation for 
  !  non-periodic function, i.e solving of
  !      Atilde| x | = | u |
  !            | y | = | v |
  ! where x = (c0,...,cn)t
  !       y = (cn+1,c-1)t
  !       u = (rhs(0),...,rhs(n))t
  ! and
  !       v = (rhs'(n),rhs'(0))t
  !       v = (sigma1,sigma2)t
  ! where sigma1 = h/3*rhs'(n) and
  !       sigma2 = h/3*rhs'(0)
  ! with the equivalent of a Shur complement method
  !  Rk1 : rhs'(0) and rhs'(n) are approximated
  !       with Lagrange polynomials 
  !  Rk2 : The boundary conditions are given by
  !    BC_left and BC_right which are equal to
  !     . 0 if Dirichlet condition (func=0)
  !     . 1 if Neumann (func'=0)
  !     . 2 if Hermite conditions 
  !    (func'=approximation of the derivate)
  !----------------------------------------------------
  subroutine natural_spline_coef(nsthis,rhs, &
    BC_left,BC_right,Sderiv_rhs)
    type(nspline1d)            , intent(inout)        :: nsthis
    real(RKIND), dimension(0:) , intent(in)           :: rhs
    integer                    , intent(in)           :: BC_left
    integer                    , intent(in)           :: BC_right
    ! deriv_rhs(0) = rhs'(n) and deriv_rhs(1) = rhs'(0)
    real(RKIND), dimension(0:) , intent(in), optional :: Sderiv_rhs 
    
    integer                     :: i, n, err
    real(RKIND)                 :: coeff
    real(RKIND), dimension(0:1) :: deriv_rhs, rhs2d, yprime 
     
    n = nsthis%n
      
    !*** Solving of x'=A-1.u ***
    do i = 0,n
      nsthis%xprime(i) = rhs(i)
    enddo
    !-> if dirichlet condition at r=rmin
    if (BC_left.eq.0) then
      nsthis%xprime(0) = 0._RKIND
    end if
    !-> if dirichlet condition at r=rmax
    if (BC_right.eq.0) then
      nsthis%xprime(n) = 0._RKIND
    end if
    call dpttrs(n+1,1,nsthis%Adiag,nsthis%Aodiag, &
      nsthis%xprime,n+1,err)
    if (err.ne.0) then
      write(6,*) 'Natural case : problem in x''=A-1.u solving'
      stop
    end if
      
    !*** computation of rhs'(0) and rhs'(n)  ***
    !-> deriv_rhs(0) = rhs'(n) and deriv_rhs(1) = rhs'(0)
    if (.not.present(Sderiv_rhs)) then
      call compute_natural_BC(rhs,nsthis%n,nsthis%h, &
        BC_left,BC_right,deriv_rhs)
    else
      deriv_rhs(0) = Sderiv_rhs(0)
      deriv_rhs(1) = Sderiv_rhs(1)
    end if
      
    !*** v-lamda.A-1u assembling ***
    coeff    = nsthis%h/3._RKIND
    rhs2d(0) = coeff*deriv_rhs(0)+nsthis%xprime(n-1)
    rhs2d(1) = coeff*deriv_rhs(1)-nsthis%xprime(1)
      
    !*** Solving of the 2X2 system :      ***
    !***   deltab.yprime = v-lamda.A-1u   ***
    yprime(0) = (1._RKIND/nsthis%det_deltab)* &
      (nsthis%eta4b*rhs2d(0)-nsthis%eta2b*rhs2d(1))
    yprime(1) = (1._RKIND/nsthis%det_deltab)* &
      (-nsthis%eta3b*rhs2d(0)+nsthis%eta1b*rhs2d(1))
		
    !*** Computation of x = x'-A-1.gamma.y ***
    do i = 0,n
      nsthis%scoef(i) = nsthis%xprime(i) - &
        nsthis%Am1gamma(i,0)*yprime(0) &
        -nsthis%Am1gamma(i,1)*yprime(1)
    enddo
    nsthis%scoef(-1)  = yprime(1)
    nsthis%scoef(n+1) = yprime(0)
  end subroutine natural_spline_coef
    
      
  !-----------------------------------------------------
  ! 1D spline coefficient computation for 
  !  non-periodic function, i.e solving of
  !      Atilde| x | = | u |
  !            | y | = | v |
  ! where x = (c0,...,cn)t
  !       y = (cn+1,c-1)t
  !       u = (rhs(0),...,rhs(n))t
  ! and
  !       v = (rhs'(n),rhs'(0))t
  !       v = (sigma1,sigma2)t
  ! where sigma1 = h/3*rhs'(n) and
  !       sigma2 = h/3*rhs'(0)
  ! with the equivalent of a Shur complement method
  !  Rk1 : rhs'(0) and rhs'(n) are approximated
  !       with Lagrange polynomials 
  !  Rk2 : The boundary conditions are given by
  !    BC_left and BC_right which are equal to
  !     . 0 if Dirichlet condition (func=0)
  !     . 1 if Neumann (func'=0)
  !     . 2 if Hermite conditions 
  !        (func'=approximation of the derivate)
  !-----------------------------------------------------
  subroutine natural_omp_spline_coef(nsthis,rhs, &
    BC_left,BC_right,scoef,Sderiv_rhs)
    type(nspline1d)            , intent(in)           :: nsthis
    real(RKIND), dimension(:)  , pointer              :: rhs
    integer                    , intent(in)           :: BC_left
    integer                    , intent(in)           :: BC_right
    ! deriv_rhs(0) = rhs'(n) and deriv_rhs(1) = rhs'(0)
    real(RKIND), dimension(:)  , pointer              :: scoef
    real(RKIND), dimension(0:) , intent(in), optional :: Sderiv_rhs 
      
    integer                            :: i, n, err
    real(RKIND)                        :: coeff
    real(RKIND), dimension(0:1)        :: deriv_rhs, rhs2d, yprime 
    real(RKIND), dimension(:), pointer :: xprime
     
    xprime => scoef
    n = nsthis%n
    do i = 0,n
      xprime(i) = rhs(i)
    end do
      
    !*** Solving of x'=A-1.u ***
    do i = 0,n
      xprime(i) = rhs(i)
    enddo
    !-> if dirichlet condition at r=rmin
    if (BC_left.eq.0) then
      xprime(0) = 0._RKIND
    end if
    !-> if dirichlet condition at r=rmax
    if (BC_right.eq.0) then
      xprime(n) = 0._RKIND
    end if
    call dpttrs(n+1,1,nsthis%Adiag,nsthis%Aodiag,xprime(0),n+1,err)
    if (err.ne.0) then
      write(6,*) 'Natural case : problem in x''=A-1.u solving'
      stop
    end if
    
    !*** computation of rhs'(0) and rhs'(n)  ***
    !-> deriv_rhs(0) = rhs'(n) and deriv_rhs(1) = rhs'(0)
    if (.not.present(Sderiv_rhs)) then
      call compute_natural_BC(rhs,nsthis%n,nsthis%h, &
        BC_left,BC_right,deriv_rhs)
    else
      deriv_rhs(0) = Sderiv_rhs(0)
      deriv_rhs(1) = Sderiv_rhs(1)
    end if
      
    !*** v-lamda.A-1u assembling ***
    coeff    = nsthis%h/3._RKIND
    rhs2d(0) = coeff*deriv_rhs(0)+xprime(n-1)
    rhs2d(1) = coeff*deriv_rhs(1)-xprime(1)
      
    !*** Solving of the 2X2 system :    ***
    !***   deltab.yprime = v-lamda.A-1u ***
    yprime(0) = (1._RKIND/nsthis%det_deltab)* &
      (nsthis%eta4b*rhs2d(0)-nsthis%eta2b*rhs2d(1))
    yprime(1) = (1._RKIND/nsthis%det_deltab)* &
      (-nsthis%eta3b*rhs2d(0)+nsthis%eta1b*rhs2d(1))
		
    !*** Computation of x = x'-A-1.gamma.y ***
    do i = 0,n
      scoef(i) = xprime(i)-nsthis%Am1gamma(i,0)*yprime(0) &
        -nsthis%Am1gamma(i,1)*yprime(1)
    enddo
    scoef(-1)  = yprime(1)
    scoef(n+1) = yprime(0)
  end subroutine natural_omp_spline_coef
 
      
  !--------------------------------------------- 
  !  natural function interpolation by using the 
  !  cubic spline coefficients calculated 
  !---------------------------------------------
  subroutine natural_interpol1d(nsthis,x,fgrid,BC_left,BC_right, &
    nbstar,xstar,finterpol)
    type(nspline1d)           , intent(inout) :: nsthis
    real(RKIND), dimension(0:), intent(in)    :: x, fgrid
    integer                   , intent(in)    :: BC_left
    integer                   , intent(in)    :: BC_right
    integer                   , intent(in)    :: nbstar
    real(RKIND), dimension(0:), intent(in)    :: xstar
    real(RKIND), dimension(0:), intent(inout) :: finterpol
      
    integer                      :: i, ipos, k    
    real(RKIND), dimension(-1:2) :: sbase
    
    !*** spline coefficient computation *** 
    call natural_spline_coef(nsthis,fgrid,BC_left,BC_right) 
      
    do i = 0, nbstar
      !*** array position location ***
      call locate(xstar(i),x,nsthis%n,nsthis%h, &
        " natural_interpol1d "//char(0),ipos)
    
      !*** Calculation of cubic spline basis ***
      call spline_basis(x(ipos),xstar(i),x(ipos+1),nsthis%h,sbase)
      
      !*** Computes f(x*) ***
      finterpol(i) = 0._RKIND
      do k=-1,2 	
        finterpol(i) = finterpol(i)+nsthis%scoef(ipos+k)*sbase(k)
      end do            
    enddo
  end subroutine natural_interpol1d
  
      
  !----------------------------------------------------------- 
  ! Compute the first derivative of a non-periodic function
  !-----------------------------------------------------------
  subroutine natural_deriv1d(nsthis,n,h, &
    BC_left,BC_right,func,dfuncdx)
    type(nspline1d)            , intent(inout) :: nsthis
    integer                    , intent(in)    :: n
    real(RKIND)                , intent(in)    :: h
    integer                    , intent(in)    :: BC_left
    integer                    , intent(in)    :: BC_right
    real(RKIND), dimension(0:n), intent(in)    :: func
    real(RKIND), dimension(0:n), intent(inout) :: dfuncdx
      
    integer                      :: ix, i
    real(RKIND), dimension(-1:2) :: sbase_prime
    real(RKIND)                  :: dfuncdx_tmp
      
    !*** computation of the spline coefficients of the function ***
    call natural_spline_coef(nsthis,func,BC_left,BC_right)
    
    !*** computation of the first derivative of the function ***
    do ix = 0,n
      call spline_basisderiv(ix,n,h,sbase_prime)
      dfuncdx_tmp = 0._RKIND	  
      do i = -1,2
        dfuncdx_tmp = dfuncdx_tmp + & 
          sbase_prime(i) * nsthis%scoef(nsthis%ipos1d(ix)+i) 
      enddo
      dfuncdx(ix) = dfuncdx_tmp 
    enddo
  end subroutine natural_deriv1d
      
  !----------------------------------------------------------- 
  ! Compute the first derivative of a non-periodic function
  !-----------------------------------------------------------
  subroutine natural_omp_deriv1d(nsthis,n,h,&
    BC_left,BC_right,rhs,scoef,drhs_dr)
    type(nspline1d)          , intent(in) :: nsthis
    integer                  , intent(in) :: n
    real(RKIND)              , intent(in) :: h
    integer                  , intent(in) :: BC_left
    integer                  , intent(in) :: BC_right
    real(RKIND), dimension(:), pointer    :: rhs
    real(RKIND), dimension(:), pointer    :: scoef
    real(RKIND), dimension(:), pointer    :: drhs_dr
      
    integer                      :: ix, i
    real(RKIND), dimension(-1:2) :: sbase_prime
    real(RKIND)                  :: dfuncdx_tmp
      
    !*** computation of the spline coefficients of the function ***
    call natural_omp_spline_coef(nsthis,rhs,BC_left,BC_right,scoef)
      
    !*** computation of the first derivative of the function ***
    do ix = 0,n
      call spline_basisderiv(ix,n,h,sbase_prime)
      dfuncdx_tmp = 0._RKIND	  
      do i = -1,2
        dfuncdx_tmp = dfuncdx_tmp + & 
          sbase_prime(i) * scoef(nsthis%ipos1d(ix)+i) 
      enddo
      drhs_dr(ix) = dfuncdx_tmp 
    enddo
  end subroutine natural_omp_deriv1d
      
  !*********************************
  ! FUNCTIONS FOR PERIODIC SPLINES *
  !*********************************
      
  !----------------------------------------------------- 
  ! 1D spline initialisation for periodic function
  !-----------------------------------------------------   
  subroutine new_spline1d_period(psthis,n,h)
    use globals, only : memory_test
    type(pspline1d), intent(out) :: psthis
    real(RKIND)    , intent(in)  :: h
    integer        , intent(in)  :: n
    
    integer     :: i, err
    real(RKIND) :: coeff1, coeff2
      
    psthis%n = n
    psthis%h = h
      
    !*** memory allocation ***
    call glob_allocate(psthis%xprime,0,n,'xprime')
    call glob_allocate(psthis%rhs,0,n,'rhs')
    call glob_allocate(psthis%ipos1d,0,n,'ipos1d')
    call glob_allocate(psthis%Adiag,0,n,'Adiag')
    call glob_allocate(psthis%Aodiag,0,n-1,'Aodiag')
    call glob_allocate(psthis%Am1gamma,0,n,0,1,'Am1gamma')
    call glob_allocate(psthis%scoef,-1,n+1,'scoef')
      
    if (.not.memory_test) then
      !*** A factorisation ***
      psthis%Adiag  = 4._RKIND
      psthis%Aodiag = 1._RKIND
      call dpttrf(n+1,psthis%Adiag,psthis%Aodiag,err)
      if (err.ne.0) then
        write(6,*) 'Periodic case : problem in the A factorisation'
        stop
      end if
      
      !*** ipos1d array initialization ***
      do i = 0,n-1
        psthis%ipos1d(i) = i
      enddo
      psthis%ipos1d(n) = n-1
      
      !*** Computation of A-1.gamma ***
      psthis%Am1gamma      = 0._RKIND
      psthis%Am1gamma(0,1) = 1._RKIND 
      psthis%Am1gamma(n,0) = 1._RKIND
      call dpttrs(n+1,2,psthis%Adiag,psthis%Aodiag, &
        psthis%Am1gamma,n+1,err)
      if (err.ne.0) then
        write(6,*) 'Periodic case : problem in A-1.gamma solving'
        stop
      end if
      
      !**************************************** 
      !*   2x2 matrice (deltab) computation : *
      !*    deltab = delta-lambda.A-1.gamma   *
      !*     where delta = 4 0                *
      !*                   1 0                *
      !****************************************     
      coeff1 = 3._RKIND/h
      coeff2 = 6._RKIND/(h**2)
      psthis%eta1b = coeff1*(-1._RKIND-psthis%Am1gamma(1,0)- &
        psthis%Am1gamma(n-1,0))
      psthis%eta2b = coeff1*(-1._RKIND-psthis%Am1gamma(1,1)- &
        psthis%Am1gamma(n-1,1))
      psthis%eta3b = coeff2 * &
        (-1._RKIND+2._RKIND*psthis%Am1gamma(0,0) &
        -psthis%Am1gamma(1,0)+psthis%Am1gamma(n-1,0)- &
        2._RKIND*psthis%Am1gamma(n,0))
      psthis%eta4b = coeff2 * &
        (1._RKIND+2._RKIND*psthis%Am1gamma(0,1) &
        -psthis%Am1gamma(1,1)+psthis%Am1gamma(n-1,1)- &
        2._RKIND*psthis%Am1gamma(n,1))
      
      ! Computation of deltab determinant
      psthis%det_deltab = psthis%eta1b*psthis%eta4b - &
        psthis%eta2b*psthis%eta3b
      if (psthis%det_deltab.eq.0._RKIND) then
        write(6,*) &
          'Periodic case : problem with determinant equal to'
        write(6,*) psthis%det_deltab
        stop
      end if
      
      ! Initialisation of the coefficients required for integration
      call integration_coef_boundaries(psthis%n,psthis%h, &
        psthis%indx,psthis%factor_int)         
    end if
  end subroutine new_spline1d_period
      
  !----------------------------------------------------- 
  ! 1D periodic spline destructor
  !-----------------------------------------------------   
  subroutine del_spline1d_period(psthis)
    type(pspline1d), intent(out) :: psthis
    
    !*** memory allocation ***
    call glob_deallocate(psthis%xprime)
    call glob_deallocate(psthis%rhs)
    call glob_deallocate(psthis%ipos1d)
    call glob_deallocate(psthis%Adiag)
    call glob_deallocate(psthis%Aodiag)
    call glob_deallocate(psthis%Am1gamma)
    call glob_deallocate(psthis%scoef)
  end subroutine del_spline1d_period
      
  !------------------------------------------------ 
  ! 1D spline coefficient computation for 
  !  periodic function, i.e solving of
  !      Atilde| x | = | u |
  !            | y | = | v |
  ! where x = (c0,...,cn-2)t
  !       y = cn+1
  !       u = (rhs(0),...,rhs(n-2))t
  !       v = rhs(n-1)
  ! with the equivalent of a Shur complement method
  !------------------------------------------------
   subroutine period_spline_coef(psthis,rhs)
    type(pspline1d)            , intent(inout) :: psthis
    real(RKIND), dimension(0:) , intent(in)    :: rhs
    
    integer                     :: i, n, err
    real(RKIND), dimension(0:1) :: rhs2d, yprime  
    real(RKIND)                 :: coeff1, coeff2    
      
    n = psthis%n    
    !*** Solving of x'=A-1.u ***
    do i = 0,n
      psthis%xprime(i) = rhs(i)
    enddo
    call dpttrs(n+1,1,psthis%Adiag,psthis%Aodiag, &
      psthis%xprime,n+1,err)
    if (err.ne.0) then
      write(6,*) 'Periodic case : problem in x''=A-1.u solving'
      stop
    end if
      
    !*** v-lamda.A-1u assembling ***
    coeff1   = 3._RKIND/psthis%h
    coeff2   = 6._RKIND/(psthis%h**2)
    rhs2d(0) = -coeff1*(psthis%xprime(1)+psthis%xprime(n-1))
    rhs2d(1) = -coeff2 * &
      (-2._RKIND*psthis%xprime(0) + &
      psthis%xprime(1)-psthis%xprime(n-1)+ &
      2._RKIND*psthis%xprime(n))
    
    !*** Solving of the 2X2 system :     ***
    !***    deltab.yprime = v-lamda.A-1u ***
    yprime(0) = (1._RKIND/psthis%det_deltab)* &
                (psthis%eta4b*rhs2d(0)-psthis%eta2b*rhs2d(1))
    yprime(1) = (1._RKIND/psthis%det_deltab)* &
                (-psthis%eta3b*rhs2d(0)+psthis%eta1b*rhs2d(1))
		
    !*** Computation of x = x'-A-1.gamma.y ***
    do i = 0,n
      psthis%scoef(i) = psthis%xprime(i) - &
        psthis%Am1gamma(i,0)*yprime(0) &
        -psthis%Am1gamma(i,1)*yprime(1)
    enddo
    psthis%scoef(n+1) = yprime(0)
    psthis%scoef(-1)  = yprime(1)
  end subroutine period_spline_coef
      
  !------------------------------------------------ 
  ! 1D spline coefficient computation for 
  !  periodic function, i.e solving of
  !      Atilde| x | = | u |
  !            | y | = | v |
  ! where x = (c0,...,cn-2)t
  !       y = cn+1
  !       u = (rhs(0),...,rhs(n-2))t
  !       v = rhs(n-1)
  ! with the equivalent of a Shur complement method
  !------------------------------------------------
  subroutine period_omp_spline_coef(psthis,rhs,scoef)
    type(pspline1d)            , intent(in) :: psthis
    real(RKIND), dimension(:)  , pointer    :: rhs
    real(RKIND), dimension(:)  , pointer    :: scoef
      
    integer                            :: tid_tmp
    integer                            :: i, n, err
    real(RKIND), dimension(0:1)        :: rhs2d, yprime  
    real(RKIND)                        :: coeff1, coeff2
    real(RKIND), dimension(:), pointer :: xprime 
      
    !*** Solving of x'=A-1.u ***
    xprime => scoef
    n = psthis%n
    do i = 0,n
      xprime(i) = rhs(i)
    end do
    call dpttrs(n+1,1,psthis%Adiag,psthis%Aodiag,xprime(0),n+1,err)
    if (err.ne.0) then
      write(6,*) 'Periodic case : problem in x''=A-1.u solving'
      stop
    end if
      
    !*** v-lamda.A-1u assembling ***
    coeff1   = 3._RKIND/psthis%h
    coeff2   = 6._RKIND/(psthis%h**2)
    rhs2d(0) = -coeff1*(xprime(1)+xprime(n-1))
    rhs2d(1) = -coeff2*(-2._RKIND*xprime(0)+xprime(1)-xprime(n-1)+ &
      2._RKIND*xprime(n))
    
    !*** Solving of the 2X2 system :   ***
    !***  deltab.yprime = v-lamda.A-1u ***
    yprime(0) = (1._RKIND/psthis%det_deltab)* &
                (psthis%eta4b*rhs2d(0)-psthis%eta2b*rhs2d(1))
    yprime(1) = (1._RKIND/psthis%det_deltab)* &
                (-psthis%eta3b*rhs2d(0)+psthis%eta1b*rhs2d(1))
		
    !*** Computation of x = x'-A-1.gamma.y ***
    do i = 0,n
      scoef(i) = xprime(i)-psthis%Am1gamma(i,0)*yprime(0) &
                      -psthis%Am1gamma(i,1)*yprime(1)
    enddo
    scoef(n+1) = yprime(0)
    scoef(-1)  = yprime(1)
  end subroutine period_omp_spline_coef
      
  !--------------------------------------------- 
  !  periodic function interpolation by using the 
  !  cubic spline coefficients calculated 
  !---------------------------------------------  
  subroutine period_interpol1d(psthis,x,fgrid,nbstar, &
    xstar,finterpol)
    type(pspline1d)                   , intent(inout) :: psthis
    real(RKIND), dimension(0:psthis%n), intent(in)    :: x, fgrid
    integer                           , intent(in)    :: nbstar
    real(RKIND), dimension(0:nbstar)  , intent(in)    :: xstar
    real(RKIND), dimension(0:nbstar)  , intent(inout) :: finterpol
      
    integer                      :: i, ipos, k    
    real(RKIND), dimension(-1:2) :: sbase
    
    !*** spline coefficient computation *** 
    call period_spline_coef(psthis,fgrid) 
      
    do i = 0, nbstar
      !*** array position location ***
      call locate(xstar(i),x,psthis%n,psthis%h, &
        " period_interpol1d "//char(0),ipos)
    
      !*** Calculation of cubic spline basis ***
      call spline_basis(x(ipos),xstar(i),x(ipos+1),psthis%h,sbase)
      
      !*** Computes f(x*) ***
      finterpol(i) = 0._RKIND
      do k=-1,2 	
        finterpol(i) = finterpol(i)+psthis%scoef(ipos+k)*sbase(k)
      end do            
    enddo
  end subroutine period_interpol1d
      
  !----------------------------------------------------------- 
  ! Compute the first derivative of a periodic function
  !-----------------------------------------------------------
  subroutine period_deriv1d(psthis,n,h,func,dfuncdx)
    type(pspline1d)            , intent(inout) :: psthis
    integer                    , intent(in)    :: n
    real(RKIND)                , intent(in)    :: h
    real(RKIND), dimension(0:n), intent(in)    :: func
    real(RKIND), dimension(0:n), intent(inout) :: dfuncdx
      
    integer                            :: ix, i 
    real(RKIND), dimension(-1:2)       :: sbase_prime
    real(RKIND), dimension(0:1)        :: deriv
    real(RKIND)                        :: dfuncdx_tmp
      
    !*** computation of the spline coefficients of the function ***
    psthis%rhs(0:n) = func(0:n)
    call period_spline_coef(psthis,psthis%rhs)
      
    !*** computation of the first derivative of the function ***
    do ix = 0,n
      call spline_basisderiv(ix,n,h,sbase_prime)
      dfuncdx_tmp = 0._RKIND	  
      do i = -1,2
        dfuncdx_tmp = dfuncdx_tmp + & 
          sbase_prime(i) * psthis%scoef(psthis%ipos1d(ix)+i) 
      enddo
      dfuncdx(ix) = dfuncdx_tmp 
    enddo
  end subroutine period_deriv1d
      
  !----------------------------------------------------------- 
  ! Compute the first derivative of a periodic function
  !-----------------------------------------------------------
  subroutine period_omp_deriv1d(psthis,n,h,rhs,scoef)
    type(pspline1d)          , intent(in) :: psthis
    integer                  , intent(in) :: n
    real(RKIND)              , intent(in) :: h
    real(RKIND), dimension(:), pointer    :: rhs
    real(RKIND), dimension(:), pointer    :: scoef
      
    integer                      :: ix, i 
    real(RKIND), dimension(-1:2) :: sbase_prime
    real(RKIND), dimension(0:1)  :: deriv
    real(RKIND)                  :: dfuncdx_tmp
      
    !*** computation of the spline coefficients of the function ***
    call period_omp_spline_coef(psthis,rhs,scoef)
      
    !*** computation of the first derivative of the function ***
    do ix = 0,n
      call spline_basisderiv(ix,n,h,sbase_prime)
      dfuncdx_tmp = 0._RKIND	  
      do i = -1,2
        dfuncdx_tmp = dfuncdx_tmp + & 
          sbase_prime(i) * scoef(psthis%ipos1d(ix)+i) 
      enddo
      rhs(ix) = dfuncdx_tmp 
    enddo
  end subroutine period_omp_deriv1d
      
  !************************************************
  ! COEFFICIENTS USED FOR INTEGRAL COMPUTATION
  !************************************************
  !----------------------------------------------------------- 
  ! Compute 2D cubic spline coefficients for a 2D array
  !-----------------------------------------------------------
  subroutine compute_scoef2D(geom,nspline1d_r,BCr_left,BCr_right,&
    pspline1d_theta,array2D,scoef2D)
    use globals, only : Rarray1_Nr, Rarray1_Ntheta
    use geometry_class
    type(geometry)       , intent(in)    :: geom
    type(nspline1d)      , intent(inout) :: nspline1d_r
    integer              , intent(in)    :: BCr_left
    integer              , intent(in)    :: BCr_right
    type(pspline1d)      , intent(inout) :: pspline1d_theta
    real(RKIND), &
      dimension(0:,0:)   , intent(in)    :: array2D
    real(RKIND), &
      dimension(-1:,-1:) , intent(out)   :: scoef2D
      
    integer                            :: itheta, ir, Nr, Ntheta
    integer                            :: i1, i2
    real(RKIND), dimension(:), pointer :: rhs1, rhs2
      
    Nr     = geom%Nr
    Ntheta = geom%Ntheta
      
    rhs1 => Rarray1_Nr
    rhs2 => Rarray1_Ntheta
      
    do itheta = 0,Ntheta-1
      do ir = 0,Nr
        rhs1(ir) = array2D(ir,itheta)
      enddo
      !*** Cubic Lagrangian interpolation for the radial ***
      !*** derivatives at the boundaries
      call natural_spline_coef(nspline1d_r,rhs1,BCr_left,BCr_right)
      do i1 = -1,Nr+1  
        scoef2D(i1,itheta) = nspline1d_r%scoef(i1)
      enddo
    enddo
    !*** Solving in the second direction ***
    do i1 = -1,Nr+1
      do itheta = 0,Ntheta-1
        rhs2(itheta) = scoef2D(i1,itheta)
      enddo
      rhs2(Ntheta) = rhs2(0)
      call period_spline_coef(pspline1d_theta,rhs2)
      do i2 = -1,Ntheta+1
        !*** cubic splines in r and theta directions ***
        scoef2D(i1,i2) = pspline1d_theta%scoef(i2)
      end do
    enddo
  end subroutine compute_scoef2D
end module spline1d_class
