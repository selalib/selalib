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
!     vlasov-poisson 1D, Landau case
!   based on the complex error function. 
!
!  - Eric Sonnendrucker: 2008/12/11
!  - Modified by Virginie Grandgirard: 2013/12/11
!-----------------------------------------------------------------------

MODULE Function_Input_Module
  
  use Precision_Module
  use plasma_functions, only: FriedConte

  IMPLICIT NONE
  integer, parameter :: max=61
  double precision :: order_magnitude = 1.0_dp
  complex(kind=dp), dimension(:,:), allocatable :: vector
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

  REAL(KIND=DP) :: kmode

CONTAINS

  SUBROUTINE FDF(omega,F,DF)
    !---------------------------------------------------
    !**PURPOSE
    !  Define the function F and its derivative DF.
    !---------------------------------------------------
    COMPLEX(KIND=DP), INTENT(IN)   :: omega
    COMPLEX(KIND=DP), INTENT(OUT)  :: F, DF

    ! local variables

    COMPLEX(KIND=DP) :: zeta, dzeta

    call FriedConte(omega/(sqrt(2.0_dp)*kmode),zeta,dzeta) 
    F= 1.0_dp-(1.0_dp/(2.0_dp*kmode**2))*dzeta
    DF= 1.0_dp/(sqrt(2.0_dp)*kmode**3) * &
      (zeta+omega/(sqrt(2.0_dp)*kmode)*dzeta)
  END SUBROUTINE FDF


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

