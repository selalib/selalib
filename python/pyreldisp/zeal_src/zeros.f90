!-----------------------------------------------------------------------
!**BEGIN PROLOGUE Zeros_Module 
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
!  This modules is basically composed of the subroutine APPROXIMATE,
!  which is used to approximate the zeros that exist inside a given box
!  as well as to find their multiplicities. It also conatains the
!  subroutine INPROD which calculates inner products of FOPs and the
!  function FQUAD which is an interface to QUADPACK.
!
!**END PROLOGUE Zeros_Module 
!-----------------------------------------------------------------------

MODULE Zeros_Module 

     USE Precision_Module
     USE Integration_Input_Module
     USE Function_Input_Module
     USE Zeal_Input_Module
     USE Error_Module

     IMPLICIT NONE

!-----------------------------------------------------------------------
!**ACCESSIBILITY
!
     PUBLIC :: APPROXIMATE 

!-----------------------------------------------------------------------
!**GLOBAL MODULE VARIABLES
!
     LOGICAL                          :: QUAD_REAL_FLAG
     INTEGER                          :: QUAD_SIDE, QUAD_DEGREE
     REAL(KIND=DP), DIMENSION(2)      :: QUAD_POINT, QUAD_STEP
     COMPLEX(KIND=DP), DIMENSION(2*M) :: QUAD_POLY


CONTAINS


     SUBROUTINE APPROXIMATE(POINT,STEP,NUMBER,ZEROS,FZEROS,  &
                          MULTIPLICITIES,NUMBERDISTINCT)
!-----------------------------------------------------------------------
!**PURPOSE
!  Given a rectangular region specified by the variables POINT and STEP, 
!  and given the total number NUMBER of zeros that lie inside this region,
!  this subroutine computes the number NUMBERDISTINCT of mutually distinct 
!  zeros that lie inside the region, approximations ZEROS for the zeros 
!  themselves, the values FZEROS of the function at these approximations,
!  and the corresponding multiplicities MULTIPLICITIES.
!
!-----------------------------------------------------------------------
!  Parameters
!
     INTEGER, INTENT(IN)  :: NUMBER
     INTEGER, INTENT(OUT) :: NUMBERDISTINCT
     INTEGER, DIMENSION(M), INTENT(OUT) :: MULTIPLICITIES 
     REAL(KIND=DP), DIMENSION(2), INTENT(IN) :: POINT, STEP
     COMPLEX(KIND=DP), DIMENSION(M), INTENT(OUT) :: ZEROS, FZEROS

!-----------------------------------------------------------------------
!  Variables used by ZGEGV
!
     INTEGER :: EIGINF
     REAL(KIND=DP), DIMENSION(8*M)    :: RWORK
     COMPLEX(KIND=DP), DIMENSION(M)   :: ALPHA, BETA
     COMPLEX(KIND=DP), DIMENSION(2*M) :: WORK
     COMPLEX(KIND=DP), DIMENSION(1,1) :: VL, VR
     COMPLEX(KIND=DP), DIMENSION(M,M) :: A, B 

!-----------------------------------------------------------------------
!  Variables used by ZGETRF and ZGETRS 
!
     INTEGER :: INFOFACT, INFOSUB
     INTEGER, DIMENSION(:), ALLOCATABLE :: IPIV
     COMPLEX(KIND=DP), DIMENSION(:), ALLOCATABLE :: COMPLEX_MULT
     COMPLEX(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: V

!-----------------------------------------------------------------------
!  Other local variables
!
     INTEGER :: DEG1, DEG2, R, T, J, K, Q, QSTART, L, TT 
     INTEGER :: ROUNDM, NDISTINCT
     LOGICAL :: TAKE_REGULAR, ALLSMALL
     REAL(KIND=DP) :: ACCPARAM, ABSIP
     REAL(KIND=DP) :: LEFTR, LEFTI, RIGHTR, RIGHTI
     COMPLEX(KIND=DP) :: CENTER, MU, NORM2, IP, Z, F, DF
     COMPLEX(KIND=DP), DIMENSION(M) :: POLY1, POLY2, ZEROS_TMP 
     COMPLEX(KIND=DP), DIMENSION(:), ALLOCATABLE :: ORDMOM
     COMPLEX(KIND=DP), DIMENSION(:,:), ALLOCATABLE :: FOPS, G, G1

!-----------------------------------------------------------------------
     INTRINSIC CMPLX, MAX, NINT, REAL
     EXTERNAL ZGEGV, ZGETRF, ZGETRS

!
!  Calculate the arithmetic mean MU of the zeros.
!
     CENTER = CMPLX(POINT(1)+0.5_DP*STEP(1),POINT(2)+0.5_DP*STEP(2),DP)
     DEG1 = 1
     POLY1(1) = CENTER
     DEG2 = 0
     CALL INPROD(POINT,STEP,POLY1,DEG1,POLY2,DEG2,MU,ACCPARAM)
     IF ( INFO .EQ. 4 ) RETURN
     MU = MU/NUMBER + CENTER
!
! Initialize FOPS, G and G1 
!
     ALLOCATE(FOPS(NUMBER,NUMBER))
     ALLOCATE(G(NUMBER,NUMBER))
     ALLOCATE(G1(NUMBER,NUMBER))
     FOPS = ZERO; G = ZERO; G1 = ZERO
     FOPS(1,1) = MU; G(1,1) = NUMBER

     R = 1; T = 0
     TAKE_REGULAR = .TRUE.

     DEG1 = R+T
     DEG2 = R+T
       POLY1(1:DEG1) = FOPS(1:DEG1,DEG1)
       POLY2(1:DEG2) = FOPS(1:DEG2,DEG2)
     CALL INPROD(POINT,STEP,POLY1,DEG1,POLY2,DEG2,NORM2,ACCPARAM)
     IF ( INFO .EQ. 4 ) RETURN
!
!  Stopping criterion
!
     ALLSMALL = .TRUE.
     TT = 0

     DO
       IF ( .NOT.(ALLSMALL .AND. TT .LE. NUMBER-2) ) EXIT

       IF ( TT .EQ. 0 ) THEN
         DEG1 = 1
         DEG2 = 1
         POLY1(1) = FOPS(1,1)
         POLY2(1) = FOPS(1,1)
         CALL INPROD(POINT,STEP,POLY1,DEG1,POLY2,DEG2,IP,ACCPARAM)
         IF ( INFO .EQ. 4 ) RETURN
         NORM2 = IP
       ELSE
         FOPS(1,TT+1) = FOPS(1,1)
         DO J = 1, TT
           FOPS(J+1,TT+1) = FOPS(J,TT)
         END DO
         DEG1 = TT+1
         DEG2 = 1
         POLY1(1:DEG1) = FOPS(1:DEG1,DEG1)
         POLY2(1) = FOPS(1,1)
         CALL INPROD(POINT,STEP,POLY1,DEG1,POLY2,DEG2,IP,ACCPARAM)
         IF ( INFO .EQ. 4 ) RETURN
       END IF
       ABSIP = ABS(IP)

       ALLSMALL = (ABSIP/ACCPARAM .LT. EPS_STOP)
       TT = TT + 1
     END DO

     IF ( ALLSMALL ) THEN
       NUMBERDISTINCT = 1
       ZEROS(1) = FOPS(1,1)
     ELSE
!
!  Main loop.
!
       EXTDO: DO
         IF ( R+T .GE. NUMBER ) EXIT EXTDO

         K = T 
!
!  Update the matrix G
!
         QSTART = MAX(R-K,0)

         DEG1 = R+K

         DO Q = QSTART, R+K-1
           POLY1(1:DEG1) = FOPS(1:DEG1,DEG1)
           IF ( Q .EQ. 0 ) THEN
             DEG2 = 0
             CALL INPROD(POINT,STEP,POLY1,DEG1,POLY2,DEG2,IP,ACCPARAM)
           ELSE
             DEG2 = Q
             POLY2(1:DEG2) = FOPS(1:DEG2,Q)
             CALL INPROD(POINT,STEP,POLY1,DEG1,POLY2,DEG2,IP,ACCPARAM)
           END IF 
           IF ( INFO .EQ. 4 ) RETURN
           G(DEG1+1,Q+1) = IP
         END DO

         IF ( .NOT. TAKE_REGULAR ) THEN
           DEG2 = R+K
           POLY1(1:DEG1) = FOPS(1:DEG1,DEG1)
           POLY2(1:DEG2) = FOPS(1:DEG2,DEG2)
           CALL INPROD(POINT,STEP,POLY1,DEG1,POLY2,DEG2,NORM2,ACCPARAM)
           IF ( INFO .EQ. 4 ) RETURN
         END IF
         G(DEG1+1,DEG1+1) = NORM2
       
         DO J = QSTART+1, DEG1
           G(J,DEG1+1) = G(DEG1+1,J)
         END DO
!
!  Update the matrix G1
!
         QSTART = MAX(R-K-1,0)

         DO Q = QSTART, DEG1
           POLY1(1:DEG1) = FOPS(1:DEG1,DEG1)
           IF ( Q .EQ. 0 ) THEN
             DEG2 = 1
             POLY2(1) = MU
             CALL INPROD(POINT,STEP,POLY1,DEG1,POLY2,DEG2,IP,ACCPARAM)
           ELSE
             DEG2 = Q+1
             POLY2(1:Q) = FOPS(1:Q,Q)
             POLY2(Q+1) = MU
             CALL INPROD(POINT,STEP,POLY1,DEG1,POLY2,DEG2,IP,ACCPARAM)
           END IF
           IF ( INFO .EQ. 4 ) RETURN
           G1(DEG1+1,Q+1) = IP
         END DO

         DO J = QSTART+1, DEG1
           G1(J,DEG1+1) = G1(DEG1+1,J)
         END DO
!
!  Solve the generalized eigenvalue problem G1 - \lambda G
!
         A(1:R+T+1,1:R+T+1) = G1(1:R+T+1,1:R+T+1)
         B(1:R+T+1,1:R+T+1) = G(1:R+T+1,1:R+T+1)

         CALL ZGEGV('N','N',R+T+1,A,M,B,M,ALPHA,BETA,    &
                    VL,1,VR,1,WORK,2*M,RWORK,EIGINF)

         IF ( EIGINF .NE. 0 ) INFO = 4
         IF ( INFO .EQ. 4 ) RETURN

       FOPS(1:R+T+1,R+T+1) = MU + ALPHA(1:R+T+1)/BETA(1:R+T+1)
!
!  Do we take a regular step or an inner step?
!
         LEFTR = POINT(1) - 0.1_DP*STEP(1)
         RIGHTR = POINT(1) + 1.1_DP*STEP(1)
         LEFTI = POINT(2) - 0.1_DP*STEP(2)
         RIGHTI = POINT(2) + 1.1_DP*STEP(2)

         TAKE_REGULAR = .TRUE.
         DO J = 1, R+T+1
           TAKE_REGULAR = TAKE_REGULAR         .AND. &
             LEFTR .LE. REAL(FOPS(J,R+T+1),DP) .AND. &
             REAL(FOPS(J,R+T+1),DP) .LE. RIGHTR
         END DO
         DO J = 1, R+T+1
           TAKE_REGULAR = TAKE_REGULAR       .AND. &
             LEFTI .LE. AIMAG(FOPS(J,R+T+1)) .AND. &
             AIMAG(FOPS(J,R+T+1)) .LE. RIGHTI
         END DO
!
!  Regular step
!
         IF ( TAKE_REGULAR ) THEN
         
           R = R+T+1
           T = 0
           ALLSMALL = .TRUE.
           TT = 0

           DO
             IF ( .NOT.(ALLSMALL .AND. TT .LE. NUMBER-1-R) ) EXIT
 
             IF ( TT .EQ. 0 ) THEN
               DEG1 = R
               DEG2 = R
               POLY1(1:DEG1) = FOPS(1:DEG1,DEG1)
               POLY2(1:DEG2) = FOPS(1:DEG2,DEG2) 
               CALL INPROD(POINT,STEP,POLY1,DEG1,POLY2,DEG2,IP,ACCPARAM)
               IF ( INFO .EQ. 4 ) RETURN
               NORM2 = IP
             ELSE
               FOPS(1:R,R+TT) = FOPS(1:R,R)
               FOPS(R+1:R+TT,R+TT) = FOPS(1:TT,TT)
               DEG1 = R+TT
               DEG2 = R
               POLY1(1:DEG1) = FOPS(1:DEG1,DEG1)
               POLY2(1:DEG2) = FOPS(1:DEG2,DEG2)
               CALL INPROD(POINT,STEP,POLY1,DEG1,POLY2,DEG2,IP,ACCPARAM)
               IF ( INFO .EQ. 4 ) RETURN
             END IF

             ABSIP = ABS(IP)

             ALLSMALL = (ABSIP/ACCPARAM .LT. EPS_STOP)
             TT = TT + 1
           END DO

           IF ( ALLSMALL ) THEN
             NUMBERDISTINCT = R
             ZEROS(1:R) = FOPS(1:R,R)
             EXIT EXTDO
           END IF

         ELSE
!
!  Inner step
!
           T = T + 1

           FOPS(1:R,R+T) = FOPS(1:R,R)
           FOPS(R+1:R+T,R+T) = FOPS(1:T,T)

         END IF

       END DO EXTDO

     END IF

     IF ( .NOT. ALLSMALL ) THEN
       NUMBERDISTINCT = NUMBER
       ZEROS(1:NUMBER) = FOPS(1:NUMBER,NUMBER)
     END IF
!
!  Compute the multiplicities. 
!
!  We start by calculating the first NUMBERDISTINCT ordinary moments.
!
     ALLOCATE(ORDMOM(NUMBERDISTINCT))
     ORDMOM(1) = NUMBER
     IF ( NUMBERDISTINCT .GT. 1 ) THEN
       ORDMOM(2) = NUMBER*MU
       DO J = 3, NUMBERDISTINCT
         DEG1 = J-1
         POLY1 = ZERO
         DEG2 = 0
         CALL INPROD(POINT,STEP,POLY1,DEG1,POLY2,DEG2,IP,ACCPARAM)
         IF ( INFO .EQ. 4 ) RETURN
         ORDMOM(J) = IP
       END DO     
     END IF
!
!  Then we solve a Vandermonde system.
!
     ALLOCATE(V(NUMBERDISTINCT,NUMBERDISTINCT))
     ALLOCATE(IPIV(NUMBERDISTINCT)) 
    
     DO J = 1, NUMBERDISTINCT 
       DO L = 1, NUMBERDISTINCT
         IF ( J .EQ. 1 ) THEN
           V(J,L) = CMPLX(ONE,ZERO,DP)
         ELSE
           V(J,L) = ZEROS(L)**(J-1)
         END IF
       END DO 
     END DO

     CALL ZGETRF(NUMBERDISTINCT,NUMBERDISTINCT,V,NUMBERDISTINCT, &
                 IPIV,INFOFACT)
     IF ( INFOFACT .NE. 0 ) INFO = 4
     IF ( INFO .EQ. 4 ) RETURN 

     ALLOCATE(COMPLEX_MULT(NUMBERDISTINCT))
     COMPLEX_MULT = ORDMOM
     CALL ZGETRS('N',NUMBERDISTINCT,1,V,NUMBERDISTINCT,IPIV,     &
                 COMPLEX_MULT,NUMBERDISTINCT,INFOSUB)
     IF ( INFOSUB .NE. 0 ) INFO = 4
     IF ( INFO .EQ. 4 ) RETURN
!
!  The computed multiplicities are rounded to the nearest integers.
!  Multiplicities that are equal to zero correspond to spurious zeros
!  and are thrown away. In the process, we calculate the final value
!  of NUMBERDISTINCT.
!  
     NDISTINCT = 0
     DO J = 1, NUMBERDISTINCT
       ROUNDM = NINT(REAL(COMPLEX_MULT(J),DP))
       IF ( ROUNDM .LT. 0 ) THEN
         INFO = 4
         RETURN
       END IF

       IF ( ROUNDM .GE. 1 ) THEN
         NDISTINCT = NDISTINCT + 1
         MULTIPLICITIES(NDISTINCT) = ROUNDM
         ZEROS_TMP(NDISTINCT) = ZEROS(J)
       END IF
     END DO

     IF ( NDISTINCT .EQ. 0 ) THEN
       INFO = 4
       RETURN
     END IF

     NUMBERDISTINCT = NDISTINCT
     ZEROS = ZEROS_TMP
!
!  Calculate the function values
!
     DO J = 1, NUMBERDISTINCT
       Z = ZEROS(J)
       CALL FDF(Z,F,DF)
       FZEROS(J) = F
     END DO 
!
!  Last statements of APPROXIMATE.
!
     DEALLOCATE(FOPS)
     DEALLOCATE(G)
     DEALLOCATE(G1)
     DEALLOCATE(ORDMOM)
     DEALLOCATE(V)
     DEALLOCATE(IPIV)
     DEALLOCATE(COMPLEX_MULT)

     END SUBROUTINE APPROXIMATE 


!-----------------------------------------------------------------------

     
     SUBROUTINE INPROD(POINT,STEP,POLY1,DEG1,POLY2,DEG2,RESULT,ACCPARAM) 
!-----------------------------------------------------------------------
!  Parameters
!
     INTEGER, INTENT(IN)                         :: DEG1, DEG2
     REAL(KIND=DP), INTENT(OUT)                  :: ACCPARAM
     REAL(KIND=DP), DIMENSION(2), INTENT(IN)     :: POINT, STEP
     COMPLEX(KIND=DP), INTENT(OUT)               :: RESULT
     COMPLEX(KIND=DP), DIMENSION(M), INTENT(IN)  :: POLY1, POLY2

!-----------------------------------------------------------------------
!  Local variables needed by DQAGX
!
     INTEGER                        :: KEY, IER, NEVAL, LAST
     INTEGER, PARAMETER             :: LIMIT = 5000
     INTEGER, PARAMETER             :: LENW = 4*LIMIT
     INTEGER, DIMENSION(LIMIT)      :: IWORK
     REAL(KIND=DP)                  :: RES, ACCPAR, ABSERR
     REAL(KIND=DP), DIMENSION(LENW) :: WORK

!-----------------------------------------------------------------------
!  Local variables
!
     INTEGER :: SIDE, J, K
     REAL(KIND=DP) :: RESULT_R, RESULT_I, ACC_R, ACC_I
     REAL(KIND=DP), DIMENSION(4) :: REALS, IMAGS, AP_REAL, AP_IMAG

!-----------------------------------------------------------------------
     INTRINSIC CMPLX, ABS, MAX
     EXTERNAL DQAGX


     QUAD_POINT = POINT
     QUAD_STEP = STEP

     QUAD_DEGREE = DEG1 + DEG2
     QUAD_POLY = 0.0_DP
     DO J = 1, DEG1
       QUAD_POLY(J) = POLY1(J)
     END DO
     DO J = 1, DEG2
       QUAD_POLY(DEG1+J) = POLY2(J)
     END DO
!
!  Specify which quadrature formula is to be used.
!
     KEY = 3
!
!  We consider the four sides of the box specified by POINT and STEP.
!  For each side, we calculate the real part and the imaginary part
!  of the integral.
!
     DO SIDE = 1, 4
       QUAD_SIDE = SIDE

       QUAD_REAL_FLAG = .TRUE.
       CALL DQAGX(FQUAD,ZERO,ONE,INTABS,INTREL,KEY,RES,ACCPAR, &
                  ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)
       IF ( IER .EQ. 1 .OR. IER .EQ. 3 ) THEN
         INFO = 4 
         RETURN
       END IF
       REALS(SIDE) = RES
       AP_REAL(SIDE) = ACCPAR

       QUAD_REAL_FLAG = .FALSE.        
       CALL DQAGX(FQUAD,ZERO,ONE,INTABS,INTREL,KEY,RES,ACCPAR, &
                  ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)
       IF ( IER .EQ. 1 .OR. IER .EQ. 3 ) THEN
         INFO = 4
         RETURN
       END IF
       IMAGS(SIDE) = RES
       AP_IMAG(SIDE) = ACCPAR
     END DO
  
     RESULT_R = ZERO
     RESULT_I = ZERO
     ACC_R = ZERO
     ACC_I = ZERO
     DO K = 1, 4
       RESULT_R = RESULT_R + REALS(K)
       RESULT_I = RESULT_I + IMAGS(K)
       IF ( ABS(RESULT_R) .GT. ACC_R ) ACC_R = ABS(RESULT_R)
       IF ( ABS(RESULT_I) .GT. ACC_I ) ACC_I = ABS(RESULT_I)
     END DO

     RESULT = CMPLX(RESULT_R,RESULT_I,DP)

     ACC_R = MAX(ACC_R,AP_REAL(1),AP_REAL(2),AP_REAL(3),AP_REAL(4))
     ACC_I = MAX(ACC_I,AP_IMAG(1),AP_IMAG(2),AP_IMAG(3),AP_IMAG(4))

     ACCPARAM = ABS(CMPLX(ACC_R,ACC_I,DP)) 
 
     END SUBROUTINE INPROD


!-----------------------------------------------------------------------


     FUNCTION FQUAD(T)
!-----------------------------------------------------------------------
     INTRINSIC ATAN, CMPLX, REAL, AIMAG

     INTEGER          :: K
     REAL(KIND=DP)    :: T, FQUAD
     COMPLEX(KIND=DP) :: TWOPII, DELTA, Z0, Z, F, DF, POLY, RESULT
     

     TWOPII = 8.0_DP*ATAN(ONE)*I

     IF ( QUAD_SIDE .EQ. 1 ) DELTA = QUAD_STEP(1)
     IF ( QUAD_SIDE .EQ. 2 ) DELTA = I*QUAD_STEP(2)
     IF ( QUAD_SIDE .EQ. 3 ) DELTA = -QUAD_STEP(1)
     IF ( QUAD_SIDE .EQ. 4 ) DELTA = -I*QUAD_STEP(2) 
    
     Z0 = CMPLX(QUAD_POINT(1),QUAD_POINT(2),DP)
     IF ( QUAD_SIDE .EQ. 1 ) Z = Z0 + T*QUAD_STEP(1)
     IF ( QUAD_SIDE .EQ. 2 ) Z = Z0 + QUAD_STEP(1) + I*T*QUAD_STEP(2)
     IF ( QUAD_SIDE .EQ. 3 ) &
       Z = Z0 + CMPLX(QUAD_STEP(1),QUAD_STEP(2),DP)-T*QUAD_STEP(1)
     IF ( QUAD_SIDE .EQ. 4 ) Z = Z0 + I*QUAD_STEP(2) - I*T*QUAD_STEP(2)

     CALL FDF(Z,F,DF)

     POLY = ONE 
     DO K = 1, QUAD_DEGREE 
       POLY = POLY*(Z - QUAD_POLY(K))
     END DO

     RESULT = POLY*(DF/F)*(DELTA/TWOPII)

     IF ( QUAD_REAL_FLAG ) THEN
       FQUAD = REAL(RESULT,DP)
     ELSE
       FQUAD = AIMAG(RESULT)
     END IF

     END FUNCTION FQUAD


END MODULE Zeros_Module 
