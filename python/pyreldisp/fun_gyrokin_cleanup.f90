!-----------------------------------------------------------------------
!**BEGIN PROLOGUE Function_Input_Module 
!**DATE WRITTEN   1999/06/01  (year/month/day)
!**AUTHORS
!
!  Peter Kravanja, Marc Van Barel
!
!    Department of Computer Science
!    Katholieke Universiteit Leuven
!    Celestijnenlaan 200 A
!    B-3001 Heverlee
!    BELGIUM
!    E-mail: Peter.Kravanja@na-net.ornl.gov
!            Marc.VanBarel@cs.kuleuven.ac.be
!
!  Omiros Ragos, Michael N. Vrahatis, and Filareti A. Zafiropoulos
!
!    Department of Mathematics
!    University of Patras
!    GR-261.10 Patras
!    GREECE
!    E-mail: {ragos,vrahatis,phikapa}@math.upatras.gr
!
!**DESCRIPTION
!
!  The user can specify the function (as well as its first derivative)
!  that is to be considered, the region where this function is analytic
!  and the values of the parameters NEWTONZ and NEWTONF, which are used
!  in the iterative refinement procedure that improves the computed 
!  approximations for the zeros. 
!
!**END PROLOGUE Function_Input_Module
!-----------------------------------------------------------------------
!
!  Modified to add plasma dispersion functions for 
!    4D drift-kinetic cylindrical case
!   based on the complex error function. 
!
!  - Eric Sonnendrucker: 2008/12/11
!  - Modified by Virginie Grandgirard: 2013/12/11
!-----------------------------------------------------------------------

MODULE Function_Input_Module

  use Precision_Module

  implicit none

  integer, parameter :: max=61
  double precision :: ordre_grandeur = 1
  complex (kind=dp), allocatable :: vector(:,:)
  !---------------------------------------------------------------------
  !**ACCESSIBILITY
  !
  public

  !---------------------------------------------------------------------
  !  If ICON = 3 or 4 (as specified in Zeal_Input_Module), then specify 
  !  the values of NEWTONZ and NEWTONF.
  !  These variables are used as follows. The modified Newton's method,
  !  which takes into account the multiplicity of a zero and converges
  !  quadratically, is used to refine the calculated approximations for
  !  the zeros. The iteration stops if
  !    - the relative distance between two successive approximations is
  !      at most NEWTONZ
  !  or
  !    - the absolute value of the function at the last approximation is
  !      at most NEWTONF
  !  or if a maximum number of iterations is exceeded.
  !
  real (kind=dp) :: NEWTONZ = 5.0E-08_DP
  real (kind=dp) :: NEWTONF = 1.0E-14_DP
  real (kind=dp) :: PI = 3.1415926535897931_DP

  real (kind=dp) :: kpar, ktheta
  integer        :: NNr
  real (kind=dp) :: dr
  real (kind=dp) :: B0

  real (kind=dp), allocatable :: Ti(:)
  real (kind=dp), allocatable :: dlogTi(:)
  real (kind=dp), allocatable :: Te(:)
  real (kind=dp), allocatable :: dlogn0(:)
  real (kind=dp), allocatable :: ddlogn0(:)
  real (kind=dp), allocatable :: rmesh(:)

  real (kind=dp), allocatable :: btheta(:)
  real (kind=dp), allocatable :: bz(:)

  ! matrix of operator
  complex (kind=dp), allocatable :: A(:,:) 
  ! derivative of A (diagonal)
  complex (kind=dp), allocatable :: Ap(:)   

contains

  subroutine FDF(omega,F,DF)
    !-------------------------------------------------------------
    !**PURPOSE
    !  Define the function F and its derivative DF.
    !-------------------------------------------------------------
    complex (kind=dp), intent(in   ) :: omega
    complex (kind=dp), intent(  out) :: F, DF
    complex (kind=dp)                :: res1, res2

    call Matrice(omega)
    call Dmatrice(omega)
    call det2(res1, res2)
    F  = res1/ordre_grandeur
    DF = res2/ordre_grandeur

  end subroutine FDF


  !-------------------------------------------------------------------- 
  subroutine matrice(omega)
    ! Fill in terms of matrix A for given value of omega
    ! A is tridiagonal and band stored : 
    !       A(1,j) is the upper diagonal term in column j
    !       A(2,j) is the diagonal term in column j
    !       A(3,j) is the lower diagonal term in column j
    !-----------------------------------------------------
    complex (kind=dp), intent(in) :: omega
    complex (kind=dp)             :: zeta, dzeta, rho
    real    (kind=dp)             :: kpar_modif
    integer                       :: i,j
    
    do i = 1, NNr 
      kpar_modif = (kpar*bz(i)+(ktheta*btheta(i)/rmesh(i)))*sqrt(2.0_dp*Ti(i))
      rho = omega/(kpar_modif)
      call FriedConte(rho,zeta,dzeta)
      A(2,i) = 2._dp+ dr*dr*( (ktheta**2-0.25_dp)/rmesh(i)**2 &
        !+ 1._dp/(Zi*Te(i)) + 0.5_dp*(0.5_dp*dlogn0(i)**2 &
        + 1._dp/Te(i) + 0.5_dp*(0.5_dp*dlogn0(i)**2 &
        + dlogn0(i)/rmesh(i) + ddlogn0(i)) &
        + ( (1+rho*zeta)/Ti(i) &
        - ktheta/(rmesh(i)*B0*kpar_modif) &
        *((dlogn0(i)-0.5_DP*dlogTi(i))*zeta &
        + dlogTi(i)*rho*(1+rho*zeta))))
      A(1,i) = -1.0_dp
      A(3,i) = -1.0_dp
    end do
  end subroutine matrice


  !---------------------------------------------------------
  subroutine Dmatrice(omega)
    ! Compute Ap the term to term derivative of matrice A
    ! Ap is a diagonal matrix
    !----------------------------------------
    complex (kind=dp), intent(in) :: omega

    integer  i, j
    complex (kind=dp) :: zeta, dzeta, rho
    real (kind=dp) :: kpar_modif
    do i = 1, NNr  
      kpar_modif = (kpar*bz(i)+(ktheta*btheta(i)/rmesh(i)))*sqrt(2.0_dp*Ti(i))
      rho = omega/(kpar_modif)
      call FriedConte(rho,zeta,dzeta)
      Ap(i) = dr*dr* &
        ((zeta+rho*dzeta)/Ti(i) - &
        ktheta/(B0*rmesh(i)*kpar_modif) * &
        (dlogn0(i)*dzeta + dlogTi(i) * &
        (1._dp+2._dp*rho*zeta+(rho**2-0.5_dp)*dzeta))) / &
        (kpar_modif)
    end do
  end subroutine Dmatrice


  !---------------------------------------------
!  subroutine det(res)
!    complex (kind=dp) ::res
!    complex (kind=dp) :: b, c,temp
!    integer :: i, j
!
!    b = 1
!    c = A(2,1)
!    do i = 2, NNr-1
!       temp = c
!       c = (A(2,i))*c - A(1,i)*A(3,i-1)*b
!       b = temp
!    end do
!
!    res = c
!    return
!  end subroutine det


   !---------------------------------------------------------
!  subroutine det2(res1, res2)
!    ! Ap doit être la matrice des dérivées des termes de A
!    ! Ap est diagonale dans notre cas
!
!    complex (kind=dp) :: res1, res2
!    complex (kind=dp) :: b1, c1,temp1, b2, c2, temp2, c3, d1
!    integer :: i, j
!
!    b1 = 1
!    c1 = A(2,1)
!    b2 = 0
!    c2 = Ap(1)
!    do i = 2, NNr-1
!       temp1 = c1
!       temp2 = c2
!       c1 = A(2,i)*c1 - A(1,i)*A(3,i-1)*b1
!       ! Ap est diagonale
!       c2 = Ap(i)*temp1 + A(2,i)*c2 - A(1,i)*A(3,i-1)*b2
!  
!       d1 = b1
!       b1 = temp1
!       b2 = temp2
!
!    end do
!    res2 = c2
!    res1 = c1
!    return
!  end subroutine det2


  ! Implement equations (36) and (37) of
  ! D. Coulette, N. Besse / Journal of Computational Physics 248 (2013) 1-32
  subroutine det2(res1, res2)

    complex (kind=dp) :: res1, res2
    complex (kind=dp) :: D, Dp, D_tmp1, D_tmp2, Dp_tmp1, Dp_tmp2
    integer :: i

    D_tmp1=A(2,1)
    D_tmp2=1
    Dp_tmp1=Ap(1)
    Dp_tmp2=0
    do i=2,NNr-1
      D=A(2,i)*D_tmp1-D_tmp2
      Dp=Ap(i)*D_tmp1+A(2,i)*Dp_tmp1-Dp_tmp2

      D_tmp2=D_tmp1
      D_tmp1=D
      Dp_tmp2=Dp_tmp1
      Dp_tmp1=Dp
    enddo

    res1 = D
    res2 = Dp
    return
  end subroutine det2


  !---------------------------------------------------------
  subroutine init
   !complex (kind=dp) :: res, omega
    complex (kind=dp) :: res1, res2, omega
    real (kind=dp) :: tab(20)
    integer :: i

    if (.not.allocated(A)) allocate(A(3,NNr))
    if (.not.allocated(Ap)) allocate(Ap(NNr))

    do i = 1, 20
       omega = cmplx(-1+i*0.1,0.0)
       call matrice(omega)
       call det2(res1,res2)
       tab(i) = real(res1)
    end do
    ordre_grandeur = maxval(tab)
  end subroutine init


  !---------------------------------------------------
  subroutine kernel(omega)
    complex (kind=dp), intent(in) :: omega
    complex (kind=dp)             :: res1,res2
    complex (kind=dp)             :: pt    (NNr,NNr)
    complex (kind=dp)             :: Q     (NNr,NNr)
    complex (kind=dp)             :: C     (NNr,NNr)
    complex (kind=dp)             :: work  (NNr    )
    real    (kind=dp)             :: D     (NNr    )
    real    (kind=dp)             :: E     (NNr-1  )
    real    (kind=dp)             :: rwork (NNr    )
    real    (kind=dp)             :: rrwork(4*NNr  )
    real    (kind=dp)             :: EE    (NNr    )
    integer                       :: i, infosub, j, ii
    
    call Matrice(omega)

    call zgbbrd('B',NNr,NNr,0,1,1,A,3,D,E,Q,NNr,PT,NNr,C,NNr,work,rwork,infosub)
    if (infosub .NE. 0) then
       print*, 'zeal failed because there is an error in the compute &
               &of singular value decomposition'
       stop
    end if

    do i = 1 , NNr-1
       EE(i) = E(i)
    end do

    call zbdsqr('U',NNr,NNr,NNr,0,D,EE,PT,NNr,Q,NNr,C,1,rrwork,infosub)
    if (infosub .NE. 0) then
       print*, 'zeal failed because there is an error in the compute &
               &of singular value decomposition'
       stop
    end if
    ! Singular values are given in D in decreasing order, so
    ! the kernel corresponds to the last singular values.
    ii = 0
    do while(D(NNr-ii).lt.1.e-10)
       ii=ii+1
    end do
    print '(a,i0,/)', 'Dimension of kernel = ', ii
    print '(a,f16.14,a,f16.14,a)', &
          'omega = ', real(omega), '+', imag(omega), 'j'
   !print*, 'omega =', omega
    if (ii.eq. 0) then
       print*, 'zealpy failed because A(', omega, ') has 0-dimensional kernel'
    else
       if (.not.allocated(vector)) allocate(vector(NNr, ii))

       do i = 1, NNr
          do j = 1, ii
             vector(i,j) = PT(NNr - ii + j,i)
             ! print*,'vect propre ',i,vector(i,j)
          end do
       end do
    end if
  end subroutine kernel


  !-------------------------------------------------------------------
  FUNCTION VALREG(LV,H)
    !------------------------------------------------------------------
    !**PURPOSE
    !  Given a rectangular region specified by LV and H (cf. the module
    !   Zeal_Input_Module), decide whether the function is analytic 
    !   inside this region or not.
    !-------------------------------------------------------------------
    LOGICAL VALREG
    REAL(KIND=DP), INTENT(IN) :: LV(2), H(2)

    VALREG = .TRUE.
    !  The following statement can be used for functions that have a
    !  branch cut along the non-positive real axis.
    !
    !   VALREG = .NOT. ( LV(2)*(LV(2)+H(2)) <= ZERO .AND. LV(1) <= ZERO )
  END FUNCTION VALREG

END MODULE Function_Input_Module

