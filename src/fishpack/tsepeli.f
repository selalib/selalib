c
c     file tsepeli.f
c
c  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c  .                                                             .
c  .                  copyright (c) 2004 by UCAR                 .
c  .                                                             .
c  .       UNIVERSITY CORPORATION for ATMOSPHERIC RESEARCH       .
c  .                                                             .
c  .                      all rights reserved                    .
c  .                                                             .
c  .                                                             .
c  .                      FISHPACK version 5.0                   .
c  .                                                             .
c  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C     *                                                               *
C     *                        F I S H P A C K                        *
C     *                                                               *
C     *                                                               *
C     *     A PACKAGE OF FORTRAN SUBPROGRAMS FOR THE SOLUTION OF      *
C     *                                                               *
C     *      SEPARABLE ELLIPTIC PARTIAL DIFFERENTIAL EQUATIONS        *
C     *                                                               *
C     *                  (Version 5.0 , JUNE 2004)                    *
C     *                                                               *
C     *                             BY                                *
C     *                                                               *
C     *        JOHN ADAMS, PAUL SWARZTRAUBER AND ROLAND SWEET         *
C     *                                                               *
C     *                             OF                                *
C     *                                                               *
C     *         THE NATIONAL CENTER FOR ATMOSPHERIC RESEARCH          *
C     *                                                               *
C     *                BOULDER, COLORADO  (80307)  U.S.A.             *
C     *                                                               *
C     *                   WHICH IS SPONSORED BY                       *
C     *                                                               *
C     *              THE NATIONAL SCIENCE FOUNDATION                  *
C     *                                                               *
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
 
      PROGRAM TSEPELI
      USE fish
      implicit none
      TYPE (fishworkspace) :: w
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER::M,N,NX,NY,I,J,MBDCND,NBDCND,IDMN,INTL,IORDER,IERROR
      REAL , DIMENSION(33,33) :: USOL, GRHS
      REAL , DIMENSION(33) :: BDA, BDB
      REAL :: A, B, C, D, DLX, DLY, X, AF, BF, CF, Y, DF, EF, FF, ALPHA
     1   , BETA, DUM, PERTRB, ERR, ERR2, ERR4
C-----------------------------------------------
C
C     DECLARE COEFFICIENT SUBROUTINES EXTERNAL
C
      external cofx,cofy
C
C     DEFINE ARITHMETIC FUNCTIONS GIVING EXACT SOLUTION
C
C
C     SET LIMITS ON REGION
C
      A = 0.0
      B = 1.0
      C = 0.0
      D = 1.0
C
C     SET GRID SIZE
C
      M = 32
      N = 32
      DLX = (B - A)/FLOAT(M)
      DLY = (D - C)/FLOAT(N)
      NX = M + 1
      NY = N + 1
      DO I = 1, NX
         X = A + FLOAT(I - 1)*DLX
C
C     SET SPECIFIED BOUNDARY CONDITIONS AT Y=C,D
C
         USOL(I,1) = UE(X,C)
         USOL(I,NY) = UE(X,D)
         CALL COFX (X, AF, BF, CF)
         DO J = 1, NY
            Y = C + FLOAT(J - 1)*DLY
            CALL COFY (Y, DF, EF, FF)
C
C     SET RIGHT HAND SIDE
C
            GRHS(I,J) = AF*UXXE(X,Y) + BF*UXE(X,Y) + CF*UE(X,Y) + DF*
     1         UYYE(X,Y) + EF*UYE(X,Y) + FF*UE(X,Y)
         END DO
      END DO
C
C     SET MIXED BOUNDARY CONDITIONS AT X=A,B
C
      ALPHA = 1.0
      BETA = 1.0
      DO J = 1, NY
         Y = C + FLOAT(J - 1)*DLY
         BDA(J) = UXE(A,Y) + ALPHA*UE(A,Y)
         BDB(J) = UXE(B,Y) + BETA*UE(B,Y)
      END DO
C
C     SET BOUNDARY SWITHCES
C
      MBDCND = 3
      NBDCND = 1
C
C     SET FIRST DIMENSION OF USOL,GRHS
C
      IDMN = 33
!     set for initialization of sepeli
      INTL = 0
C
C     OBTAIN SECOND ORDER APPROXIMATION
C
      IORDER = 2
      CALL SEPELI (INTL, IORDER, A, B, M, MBDCND, BDA, ALPHA, BDB, BETA
     1   , C, D, N, NBDCND, DUM, DUM, DUM, DUM, COFX, COFY, GRHS, USOL, 
     2   IDMN, W, PERTRB, IERROR)
      ERR = 0.0
      DO I = 1, NX
         X = A + FLOAT(I - 1)*DLX
         DO J = 1, NY
            Y = C + FLOAT(J - 1)*DLY
            ERR = AMAX1(ERR,ABS((USOL(I,J)-UE(X,Y))/UE(X,Y)))
         END DO
      END DO
      ERR2 = ERR
C
C     OBTAIN FOURTH ORDER APPROXIMATION
C
      IORDER = 4
C
C     NON-INITIAL CALL
C
      INTL = 1
      CALL SEPELI (INTL, IORDER, A, B, M, MBDCND, BDA, ALPHA, BDB, BETA
     1   , C, D, N, NBDCND, DUM, DUM, DUM, DUM, COFX, COFY, GRHS, USOL, 
     2   IDMN, W, PERTRB, IERROR)
C
C     COMPUTE DISCRETIZATION ERROR
C
      ERR = 0.0
      DO J = 1, NY
         Y = C + FLOAT(J - 1)*DLY
         DO I = 1, NX
            X = A + FLOAT(I - 1)*DLX
            ERR = AMAX1(ERR,ABS((USOL(I,J)-UE(X,Y))/UE(X,Y)))
         END DO
      END DO
      ERR4 = ERR
!     Print earlier output from platforms with 32 and 64 bit floating point
!     arithemtic followed by the output from this computer
      WRITE (*, *) '    SEPELI TEST RUN *** '
      WRITE (*, *) 
     1   '    Previous 64 bit floating point arithmetic result '
      WRITE (*, *) '    IERROR = 0'
      WRITE (*, *) '    Second Order Discretization Error = 9.7891E-5'
      WRITE (*, *) '    Fourth Order Discretization Error = 1.4735E-6'
      WRITE (*, *) 
     1   '    Previous 32 bit floating point arithmetic result '
      WRITE (*, *) '    IERROR = 0'
      WRITE (*, *) '    Second Order Discretization Error = 1.2708E-4'
      WRITE (*, *) '    Fourth Order Discretization Error = 3.1948E-5'
      WRITE (*, *) '    The output from your computer is: '
      WRITE (*, *) '    IERROR =', IERROR
      WRITE (*, *) '    Second Order Discretization Error =', ERR2
      WRITE (*, *) '    Fourth Order Discretization Error =', ERR4
!     release dynamically allocated real and complex work space
      CALL FISHFIN (W)
      STOP 
      CONTAINS


      REAL FUNCTION UE (S, T)
      REAL, INTENT(IN) :: S
      REAL, INTENT(IN) :: T
      UE = (S*T)**3 + 1.0
      RETURN 
      END FUNCTION UE


      REAL FUNCTION UXE (S, T)
      REAL, INTENT(IN) :: S
      REAL, INTENT(IN) :: T
      UXE = 3.0*S**2*T**3
      RETURN 
      END FUNCTION UXE


      REAL FUNCTION UXXE (S, T)
      REAL, INTENT(IN) :: S
      REAL, INTENT(IN) :: T
      UXXE = 6.0*S*T**3
      RETURN 
      END FUNCTION UXXE


      REAL FUNCTION UYE (S, T)
      REAL, INTENT(IN) :: S
      REAL, INTENT(IN) :: T
      UYE = 3.0*S**3*T**2
      RETURN 
      END FUNCTION UYE


      REAL FUNCTION UYYE (S, T)
      REAL, INTENT(IN) :: S
      REAL, INTENT(IN) :: T
      UYYE = 6.0*S**3*T
      RETURN 
      END FUNCTION UYYE
      END PROGRAM TSEPELI


      SUBROUTINE COFX(X, AF, BF, CF)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL , INTENT(IN) :: X
      REAL , INTENT(OUT) :: AF
      REAL , INTENT(OUT) :: BF
      REAL , INTENT(OUT) :: CF
C-----------------------------------------------
C
C     SET COEFFICIENTS IN THE X-DIRECTION.
C
      AF = (X + 1.)**2
      BF = 2.0*(X + 1.)
      CF = -X
      RETURN 
      END SUBROUTINE COFX


      SUBROUTINE COFY(Y, DF, EF, FF)
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL , INTENT(IN) :: Y
      REAL , INTENT(OUT) :: DF
      REAL , INTENT(OUT) :: EF
      REAL , INTENT(OUT) :: FF
C-----------------------------------------------
C
C     SET COEFFICIENTS IN Y DIRECTION
C
      DF = EXP(Y)
      EF = 0.0
      FF = -Y
      RETURN 
      END SUBROUTINE COFY
