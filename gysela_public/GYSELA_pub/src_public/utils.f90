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
      
!-----------------------------------------------------
! file : utils.f90
! date : 04/05/2001
!  routines used in efield module and in poisson class
!-----------------------------------------------------
module utils_module
  use mem_alloc_module
      
  implicit none
      
  interface tridiag
     module procedure tridiag_r, tridiag_c
  end interface
      
  !******************************
  contains
  !******************************
      
  !*******************************************************
  !            LOCALISATION IN A 1D ARRAY
  !*******************************************************
  !------------------------------------------------ 
  ! Returns the position in the array of the 
  ! first knot at the left of xstar 
  ! -> in the case of an equidistant mesh where
  !     x(i) = x(0) + i*h
  !
  ! Rk : Becareful, this subroutine is only valid
  !      for the case of x beginning at index 0
  !------------------------------------------------
  subroutine locate(xstar,x,n,h,subp,ipos_left)
    real(RKIND)                , intent(in)    :: xstar
    integer                    , intent(in)    :: n
    real(RKIND), dimension(0:n), intent(in)    :: x
    real(RKIND)                , intent(in)    :: h
    character (LEN=50)         , intent(in)    :: subp
    integer                    , intent(inout) :: ipos_left
    
    if (h.eq.0._RKIND) then
      ipos_left = 0
    else
      !*** Searchs ipos_left such that xstar(i) belongs to ***
      !***   [x(ipos_left),x(ipos_left+1)[                 ***
      ipos_left = floor( (xstar-x(0))/h )
      if ( (ipos_left<0) .or. (ipos_left>n) ) then
        write(6,*) 'Error : ', xstar, &
          ' has not be located between ', &
          x(0), 'and', x(n), ' => ipos_left = ', ipos_left, &
          " called by ",subp
        stop
      end if
      if (ipos_left.eq.n) then
        ipos_left = ipos_left-1
      end if
    end if
  end subroutine locate
      
  !--------------------------------------------- 
  ! Returns the position in the array of the 
  ! first knot at the left of xstar
  !  -> in the case of a non-equidistant mesh x
  ! (using of the dichotomie method)
  !
  ! Rk : Becareful, this subroutine is only valid
  !      for the case of x beginning at index 0
  !---------------------------------------------
  subroutine locate_nonequidistant(xstar,x,n,subp,ipos_left)
    real(RKIND)               , intent(in)  :: xstar
    integer                   , intent(in)  :: n
    real(RKIND), dimension(0:), intent(in)  :: x
    character (LEN=50)        , intent(in)  :: subp
    integer                   , intent(out) :: ipos_left 
    
    integer :: first_indx
    integer :: indx_L, indx_R    ! for left and right indexes
    integer :: indx_M            ! for middle index
      
    first_indx = 0
    ipos_left  = 0
    if ( (xstar.lt.x(first_indx)) .or. (xstar.gt.x(n)) ) then
      write(6,*) 'Error : ', xstar, &
        ' has not be located between ', &
        x(first_indx), 'and', x(n), &
        ' , called by ',subp
      stop
    endif
    indx_L = first_indx
    indx_R = n
    indx_M = floor(0.5_RKIND*(indx_L+indx_R))
    do while (indx_M.ne.indx_L)
      if (xstar.lt.x(indx_M)) then
        indx_R = indx_M
      else
        indx_L = indx_M
      end if
      indx_M = floor(0.5_RKIND*(indx_L+indx_R))
    end do
    if (.not.((xstar.ge.x(indx_L)).and.(xstar.lt.x(indx_R))) .and. (n .ne. 0)) then 
      write(6,'(A,e12.5,A,I4,A,E12.5,A,I4,A,E12.5,A,A)') &
        ' Error : ', xstar, &
        ' has not be located; indx_L ',indx_L, &
        ' x(indx_L) ',x(indx_L),' indx_R ',indx_R, &
        ' x(indx_R) ',x(indx_R),' , called by ',subp
      stop
    END if
    ipos_left = indx_L
  end subroutine locate_nonequidistant
    
      
  !*******************************************************
  !            NUMERICAL DERIVATIVES
  !*******************************************************
  !--------------------------------------------------------------
  ! first derivative computed in finite differences of 4th order
  !   input:  - U   : function to be derived
  !           - dx  : step of the mesh
  !           - nx  : number of grid points
  !           - BC  : "Boundary Conditions": 
  !                    0=non-periodic, 1=periodic
  !   output: - dxU : first derivative in x of U 
  !--------------------------------------------------------------
  subroutine deriv1(U,dxU,nx,dx,BC)
    use prec_const
      
    implicit none
    integer                      :: nx, i, BC
    real(RKIND)                  :: dx, cdx, c1a, c1b, c2a, c2b
    real(RKIND), dimension(0:nx) :: U, dxU
      
    cdx = TW/(TH*dx)
    c1a = 25._RKIND/12._RKIND
    c1b = FO/TH
    c2a = 5._RKIND/6._RKIND
    c2b = TH/TW
      
    do i=2,nx-2
      dxU(i) = cdx*(U(i+1)-U(i-1) + (U(i-2)-U(i+2))/8._RKIND)
    end do
    IF (BC.eq.0) THEN
      dxU(0)    = (-c1a*U(0)+FO*U(1)-TH*U(2)+c1b*U(3)-U(4)/FO)/dx
      dxU(1)    = (-U(0)/FO-c2a*U(1)+c2b*U(2)-HF*U(3)+U(4)/ &
        12._RKIND)/dx
      dxU(nx)   = (c1a*U(nx)-FO*U(nx-1)+TH*U(nx-2)- &
        c1b*U(nx-3)+U(nx-4)/FO)/dx
      dxU(nx-1) = (U(nx)/FO+c2a*U(nx-1)-c2b*U(nx-2) + &
        HF*U(nx-3)-U(nx-4)/12._RKIND)/dx
    ELSEIF (BC.eq.1) THEN
      dxU(0)    = cdx*( U(1)-U(nx)+(U(nx-1)-U(2))/8._RKIND )
      dxU(1)    = cdx*( U(2)-U(0)+(U(nx)-U(3))/8._RKIND )
      dxU(nx)   = cdx*( U(0)-U(nx-1)+(U(nx-2)-U(1))/8._RKIND )
      dxU(nx-1) = cdx*( U(nx)-U(nx-2)+(U(nx-3)-U(0))/8._RKIND )
    ELSE
      print*,'Error in the choice of BC'
      stop
    END IF
  end subroutine deriv1
      
  !-----------------------------------------------------------
  ! Solves for a vector U of length Nx the tridiagonal 
  !  linear set given by the equation M*U = R , 
  !  where A,B,C and R are input vectors and are not modified
  !-----------------------------------------------------------
  subroutine tridiag_r(A,B,C,R,U,Nx)
    use prec_const
      
    implicit none
    integer                        :: i, Nx
    real(RKIND)                    :: bet
    real(RKIND), dimension(0:Nx+1) :: gam
    real(RKIND), dimension(0:Nx)   :: A, B, C, R, U
      
    if (B(0).eq.0.) then
      write(6,*) 'error in tridiag_B(0)',B(0),B(1),B(Nx)
      stop
    endif
    bet  = B(0)
    U(0) = R(0) / bet
    do i = 1,Nx
      gam(i) = C(i-1)/bet
      bet    = B(i) - A(i) * gam(i)
      if (bet.eq.0.) write(6,*) 'error in tridiag_bet'
      U(i)   = (R(i) - A(i) * U(i-1)) / bet
    end do
    do i = Nx-1,0,-1
      U(i) = U(i) - gam(i+1) * U(i+1)
    end do
  end subroutine tridiag_r
      
  !------------------------------------------------
  ! Same as TRIDIAG when R and U are COMPLEX data
  !------------------------------------------------
  subroutine tridiag_c(A,B,C,R,U,Nx)
    use prec_const
      
    implicit none
    integer                           :: i, Nx
    real(RKIND)                       :: bet
    real(RKIND),    dimension(0:Nx+1) :: gam
    real(RKIND),    dimension(0:Nx)   :: A, B, C
    complex(RKIND), dimension(0:Nx)   :: R, U
      
    if (B(0).eq.0.) then 
      write(6,*) 'error in tridiag_c'
      stop
    end if
    bet  = B(0)
    U(0) = R(0) / bet
    do i = 1,Nx
      gam(i) = C(i-1)/bet
      bet    = B(i) - A(i) * gam(i)
      if (bet.eq.0.) write(6,*) 'error in tridiag_c'
      U(i)   = (R(i) - A(i) * U(i-1)) / bet
    end do
    do i = Nx-1,0,-1
      U(i) = U(i) - gam(i+1) * U(i+1)
    end do
  end subroutine tridiag_c
      
  !---------------------------------------------------------
  ! date : 10/03/2006 
  ! generator of noise of uniform density
  ! seed value : G = 1431655765._RKIND
  !  
  !   input:  - G previous random value generated
  !  
  !   output: - noise = real -1 or 1 
  !---------------------------------------------------------
  subroutine binNOISE(G, noise)
    real(RKIND), intent(inout) :: G
    real(RKIND), intent(out)   :: noise
      
    real(RKIND) :: UNIF
      
    G    = mod(G*9228907._RKIND,4294967296._RKIND) 
    UNIF = G / 4294967296._RKIND            ! real in [0,1] 
    noise = sign(1._RKIND,UNIF-0.5_RKIND)   ! real -1 or 1 
     
  end subroutine binNOISE
      
  !*******************************************************
  !            NUMERICAL INTEGRATIONS
  !*******************************************************
  !---------------------------------------------------------
  ! Computation of the coefficients required for an 
  !  integration by the trapeze method 
  !       (method of second order)
  !---------------------------------------------------------
  subroutine compute_trapeze_coeff(Nx,dx,periodic,integr_coeff)
    integer                   , intent(in)    :: Nx
    real(RKIND)               , intent(in)    :: dx
    logical                   , intent(in)    :: periodic
    real(RKIND), dimension(0:), intent(inout) :: integr_coeff
      
    integer :: ix
      
    integr_coeff = 0._RKIND
    do ix = 0,Nx-1
      integr_coeff(ix)   = integr_coeff(ix)   + 1._RKIND
      integr_coeff(ix+1) = integr_coeff(ix+1) + 1._RKIND
    end do
    integr_coeff = dx*integr_coeff/2._RKIND
    if (periodic) then
      integr_coeff(Nx) = 0._RKIND
      integr_coeff(0)  = 2._RKIND * &
        integr_coeff(0)
    end if
  end subroutine compute_trapeze_coeff
      
  !---------------------------------------------------------
  ! Computation of the coefficients required for an 
  !  integration by the Simpson method 
  !       (method of fourth order)
  !---------------------------------------------------------
  subroutine compute_simpson_coeff(Nx,dx,periodic,integr_coeff)
    integer                   , intent(in)    :: Nx
    real(RKIND)               , intent(in)    :: dx
    logical                   , intent(in)    :: periodic
    real(RKIND), dimension(0:), intent(inout) :: integr_coeff    
      
    integer :: ix
      
    if (.not.periodic) then
      integr_coeff(0)  = 1._RKIND
      integr_coeff(Nx) = 1._RKIND
    else
      integr_coeff(0)  = 2._RKIND
      integr_coeff(Nx) = 0._RKIND
    end if
    do ix = 1,Nx-1
      if (mod(ix,2).eq.0) then
        integr_coeff(ix) = 2._RKIND
      else
        integr_coeff(ix) = 4._RKIND
      end if
    end do
    integr_coeff = dx*integr_coeff/3._RKIND
  end subroutine compute_simpson_coeff
      
  !---------------------------------------------------------
  ! Computation of the coefficients required for an 
  !  integration by the Villarceau method
  !       (method of order > simpson)
  !---------------------------------------------------------
  subroutine compute_villarceau_coeff(Nx,dx,integr_coeff)
    integer                   , intent(in)    :: Nx
    real(RKIND)               , intent(in)    :: dx
    real(RKIND), dimension(0:), intent(inout) :: integr_coeff
      
    integer :: ix
      
    integr_coeff = 0._RKIND
    do ix = 0,Nx-4,4
      integr_coeff(ix)   = integr_coeff(ix)   + 14._RKIND
      integr_coeff(ix+1) = integr_coeff(ix+1) + 64._RKIND
      integr_coeff(ix+2) = integr_coeff(ix+2) + 24._RKIND
      integr_coeff(ix+3) = integr_coeff(ix+3) + 64._RKIND
      integr_coeff(ix+4) = integr_coeff(ix+4) + 14._RKIND
    end do
    integr_coeff = dx*integr_coeff/45._RKIND
  end subroutine compute_villarceau_coeff
      
  !---------------------------------------------------------
  ! Computation of the coefficients required for an 
  !  integration by the Hardy method 
  !       (method of order > Villarceau)
  !---------------------------------------------------------
  subroutine compute_hardy_coeff(Nx,dx,integr_coeff)
    integer                   , intent(in)    :: Nx
    real(RKIND)               , intent(in)    :: dx
    real(RKIND), dimension(0:), intent(inout) :: integr_coeff
      
    integer :: ix
      
    integr_coeff = 0._RKIND
    do ix = 0,Nx-6,6
      integr_coeff(ix)   = integr_coeff(ix)   + 41._RKIND
      integr_coeff(ix+1) = integr_coeff(ix+1) + 216._RKIND
      integr_coeff(ix+2) = integr_coeff(ix+2) + 27._RKIND
      integr_coeff(ix+3) = integr_coeff(ix+3) + 272._RKIND
      integr_coeff(ix+4) = integr_coeff(ix+4) + 27._RKIND
      integr_coeff(ix+5) = integr_coeff(ix+5) + 216._RKIND
      integr_coeff(ix+6) = integr_coeff(ix+6) + 41._RKIND
    end do
    integr_coeff = dx*integr_coeff/140._RKIND
  end subroutine compute_hardy_coeff
      
  !------------------------------------------------------- 
  ! Initialisation of the random seed
  !------------------------------------------------------- 
  subroutine init_random_seed()
    integer :: i
    integer :: n     = 32
    integer :: clock = 1
      
    integer, dimension(:), allocatable :: seed
      
    call RANDOM_SEED(size = n)
    allocate(seed(n))
    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    call RANDOM_SEED(PUT = seed)
    deallocate(seed)
  end subroutine init_random_seed
end module utils_module
