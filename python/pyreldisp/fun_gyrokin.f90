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
  USE Precision_Module

  IMPLICIT NONE
  integer, parameter :: max=61
  double precision :: ordre_grandeur = 1
  complex(kind=dp), dimension(:,:), allocatable :: vecteur
  !---------------------------------------------------------------------
  !**ACCESSIBILITY
  !
  PUBLIC

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
  REAL(KIND=DP) :: NEWTONZ = 5.0E-08_DP
  REAL(KIND=DP) :: NEWTONF = 1.0E-14_DP
  REAL(KIND=DP) :: PI = 3.1415926535897931_DP

  INTEGER       :: Zi
  REAL(KIND=DP) :: kpar, ktheta
  integer       :: NNr
  REAL(KIND=DP) :: dr

  REAL(KIND=DP), dimension(:), allocatable :: Ti
  REAL(KIND=DP), dimension(:), allocatable :: dlogTi
  REAL(KIND=DP), dimension(:), allocatable :: Te
  REAL(KIND=DP), dimension(:), allocatable :: n0
  REAL(KIND=DP), dimension(:), allocatable :: dlogn0
  REAL(KIND=DP), dimension(:), allocatable :: ddlogn0
  REAL(KIND=DP), dimension(:), allocatable :: rmesh

  ! matrix of operator
  COMPLEX (KIND=DP), dimension(:,:), allocatable :: A 
  ! derivative of A (diagonal)
  COMPLEX (KIND=DP), dimension(:), allocatable :: Ap   

CONTAINS

  SUBROUTINE FDF(omega,F,DF)
    !-------------------------------------------------------------
    !**PURPOSE
    !  Define the function F and its derivative DF.
    !-------------------------------------------------------------
    COMPLEX(KIND=DP), INTENT(IN)   :: omega
    COMPLEX(KIND=DP), INTENT(OUT)  :: F, DF

    ! local variables

    complex (kind=dp) :: res1, res2

       call Matrice(omega)
       call Dmatrice(omega)
       call det2(res1, res2)
       F = res1/ordre_grandeur
       DF = res2/ordre_grandeur
  END SUBROUTINE FDF


  !-------------------------------------------------------------------- 
  subroutine matrice(omega)
    ! Fill in terms of matrix A for given value of omega
    ! A is tridiagonal and band stored : 
    !       A(1,j) is the upper diagonal term in column j
    !       A(2,j) is the diagonal term in column j
    !       A(3,j) is the lower diagonal term in column j
    !-----------------------------------------------------
    COMPLEX(KIND=DP), intent(in) :: omega
    integer i, j
    COMPLEX(KIND=DP) :: zeta, dzeta, rho
    
    do i = 1, NNr 
       rho = omega/(kpar*sqrt(2.0_dp*Ti(i)))
       call FriedConte(rho,zeta,dzeta)
       A(2,i) = 2+ dr*dr*( (ktheta**2-0.25_dp)/rmesh(i)**2 &
         + 1/(Zi*Te(i)) + 0.5_dp*(0.5_dp*dlogn0(i)**2 &
         + dlogn0(i)/rmesh(i) + ddlogn0(i)) &
         + (1+rho*zeta)/Ti(i) - ktheta/(rmesh(i)*kpar*sqrt(2*Ti(i))) &
         *((dlogn0(i)-0.5_DP*dlogTi(i))*zeta &
         + dlogTi(i)*rho*(1+rho*zeta)))
       A(1,i) = -1.0_dp
       A(3,i) = -1.0_dp
    end do
  end subroutine matrice


  !---------------------------------------------------------
  subroutine Dmatrice(omega)
    ! Compute Ap the term to term derivative of matrice A
    ! Ap is a diagonal matrix
    !----------------------------------------
    COMPLEX(KIND=DP), intent(in) :: omega

    integer  i, j
    COMPLEX(KIND=DP) :: zeta, dzeta, rho

    do i = 1, NNr  
       rho = omega/(kpar*sqrt(2.0_dp*Ti(i)))
       call FriedConte(rho,zeta,dzeta)
       Ap(i) = dr*dr*((zeta+rho*dzeta)/Ti(i) - &
         ktheta/(rmesh(i)*kpar*sqrt(2*Ti(i))) * &
         (dlogn0(i)*dzeta + dlogTi(i) * &
         (1+2*rho*zeta+(rho**2-0.5_dp)*dzeta))) / &
         (kpar*sqrt(2.0_dp*Ti(i)))
    end do
  end subroutine Dmatrice


  !---------------------------------------------
  subroutine det(res)
    complex (kind=dp) ::res
    complex (kind=dp) :: b, c,temp
    integer :: i, j

    b = 1
    c = A(2,1)
    do i = 2, NNr-1
       temp = c
       c = (A(2,i))*c - A(1,i)*A(3,i-1)*b
       b = temp
    end do

    res = c
    return
  end subroutine det


  !---------------------------------------------------------
  subroutine det2(res1, res2)
    ! Ap doit être la matrice des dérivées des termes de A
    ! Ap est diagonale dans notre cas

    complex (kind=dp) :: res1, res2
    complex (kind=dp) :: b1, c1,temp1, b2, c2, temp2, c3, d1
    integer :: i, j

    b1 = 1
    c1 = A(2,1)
    b2 = 0
    c2 = Ap(1)
    do i = 2, NNr-1
       temp1 = c1
       temp2 = c2
       c1 = A(2,i)*c1 - A(1,i)*A(3,i-1)*b1
       ! Ap est diagonale
       c2 = Ap(i)*temp1 + A(2,i)*c2 - A(1,i)*A(3,i-1)*b2
  
       d1 = b1
       b1 = temp1
       b2 = temp2

    end do
    res2 = c2
    res1 = c1
    return
  end subroutine det2

  !---------------------------------------------------------
  subroutine init
    COMPLEX(KIND = DP) :: res,  omega
    real(kind = dp), dimension(20) :: tab
    integer :: i

    if (.not.allocated(A)) allocate(A(3,NNr))
    if (.not.allocated(Ap)) allocate(Ap(NNr))

    do i = 1, 20
       omega = cmplx(-1+i*0.1,0.0)
       call matrice(omega)
       call det(res)
       tab(i) = real(res)
    end do
    ordre_grandeur = maxval(tab)
  end subroutine init


  !---------------------------------------------------
  subroutine kernel(omega)
    complex(kind=dp), intent(in) :: omega
    complex(kind=dp) :: res1,res2
    complex(kind=dp), dimension(NNr,NNr) :: pt
    complex(kind=dp), dimension(NNr,NNr) :: Q
    complex(kind=dp), dimension(NNr,NNr) :: C
    complex(kind=dp), dimension(NNr)     :: WORK
    real(kind=dp)   , dimension(NNr)     :: D
    real(kind=dp)   , dimension(NNr-1)   :: E
    real(kind=dp)   , dimension(NNr)     :: RWORK
    real(kind=dp)   , dimension(4*NNr)   :: RRWORK
    real(kind=dp)   , dimension(NNr)     :: EE
    integer :: i, infosub, j, ii
    
    call Matrice(omega)

    call zgbbrd('B', NNr, NNr, 0, 1, 1, A, 3, D, E, Q, NNr, PT, NNr, C, NNr, WORK, RWORK,INFOSUB)
    if (infosub .NE. 0) then
       print*, 'zeal failed because there is an error in the compute of singular value decomposition'
       stop
    end if

    do i = 1 , NNr-1
       EE(i) = E(i)
    end do

    call zbdsqr('U', NNr, NNr, NNr, 0, D, EE , PT, NNr , Q, NNr, C, 1, rrwork, infosub)
    if (infosub .NE. 0) then
       print*, 'zeal failed because there is an error in the compute of singular value decomposition'
       stop
    end if
    ! singular values are given in D in decreasing order. So the kernel corresponds to the last
    ! singular values
    ii = 0
    do while(D(NNr-ii).lt.1.e-10)
       ii=ii+1
    end do
    print*, 'dim of kernel=',ii, 'omega= ',omega
    if (ii.eq. 0) then
       print*, 'zealpy failed because A(', omega, ') has 0-dimensional kernel'
    else
       if (.not.allocated(vecteur)) allocate(vecteur(NNr, ii))

       do i = 1, NNr
          do j = 1, ii
             vecteur(i,j) = PT(NNr - ii + j,i)
             ! print*,'vect propre ',i,vecteur(i,j)
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

