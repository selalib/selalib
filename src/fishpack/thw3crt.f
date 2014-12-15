c
c     file thw3crt.f
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
      PROGRAM THW3CRT
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: LBDCND, MBDCND, NBDCND, L, M, N, LDIMF, MDIMF, LP1, I, 
     1   MP1, J, NP1, K, IERROR
      REAL , DIMENSION(11,41,16) :: F
      REAL , DIMENSION(11,41) :: BDZF
      REAL , DIMENSION(11) :: X
      REAL , DIMENSION(41) :: Y
      REAL , DIMENSION(16) :: Z
      REAL :: ELMBDA, XS, XF, YS, PI, YF, ZS, ZF, DX, DY, DZ, BDXS, BDXF
     1   , BDYS, BDYF, BDZS, PERTRB, ERR, T
C-----------------------------------------------
 
C
C        FROM THE DESCRIPTION OF THE PROBLEM GIVEN ABOVE, WE DEFINE
C     THE FOLLOWING QUANTITIES
C
      ELMBDA = -3.
      XS = 0.
      XF = 1.
      LBDCND = 1
      YS = 0.
      PI = 4.0*ATAN(1.0)
      YF = 2.*PI
      MBDCND = 0
      ZS = 0.
      ZF = PI/2.
      NBDCND = 2
      L = 10
      M = 40
      N = 15
C
C     FROM THE DIMENSION STATEMENT ABOVE WE DEFINE
C
      LDIMF = 11
      MDIMF = 41
C
C     WE DEFINE THE GRID POINTS FOR LATER USE.
C
      LP1 = L + 1
      DX = (XF - XS)/FLOAT(L)
      DO I = 1, LP1
         X(I) = XS + FLOAT(I - 1)*DX
      END DO
      MP1 = M + 1
      DY = (YF - YS)/FLOAT(M)
      DO J = 1, MP1
         Y(J) = YS + FLOAT(J - 1)*DY
      END DO
      NP1 = N + 1
      DZ = (ZF - ZS)/FLOAT(N)
      DO K = 1, NP1
         Z(K) = ZS + FLOAT(K - 1)*DZ
      END DO
C
C     WE DEFINE THE ARRAY OF DERIVATIVE BOUNDARY VALUES.
C
      DO I = 1, LP1
         DO J = 1, MP1
            BDZF(I,J) = -X(I)**4*SIN(Y(J))
         END DO
      END DO
C
C     NOTE THAT FOR THIS EXAMPLE ALL OTHER BOUNDARY ARRAYS ARE
C     DUMMY VARIABLES.
C     WE DEFINE THE FUNCTION BOUNDARY VALUES IN THE F ARRAY.
C
      DO J = 1, MP1
         DO K = 1, NP1
            F(1,J,K) = 0.
            F(LP1,J,K) = SIN(Y(J))*COS(Z(K))
         END DO
      END DO
      DO I = 1, LP1
         DO J = 1, MP1
            F(I,J,1) = X(I)**4*SIN(Y(J))
         END DO
      END DO
C
C     WE NOW DEFINE THE VALUES OF THE RIGHT SIDE OF THE HELMHOLTZ
C     EQUATION.
C
      DO I = 2, L
         DO J = 1, MP1
            DO K = 2, NP1
               F(I,J,K) = 4.*X(I)**2*(3. - X(I)**2)*SIN(Y(J))*COS(Z(K))
            END DO
         END DO
      END DO
C
C     CALL HW3CRT TO GENERATE AND SOLVE THE FINITE DIFFERENCE EQUATION.
C
      CALL HW3CRT (XS, XF, L, LBDCND, BDXS, BDXF, YS, YF, M, MBDCND, 
     1   BDYS, BDYF, ZS, ZF, N, NBDCND, BDZS, BDZF, ELMBDA, LDIMF, MDIMF
     2   , F, PERTRB, IERROR)
C
C     COMPUTE DISCRETIZATION ERROR.  THE EXACT SOLUTION TO THE
C     PROBLEM IS
C
C        U(X,Y,Z) = X**4*SIN(Y)*COS(Z)
C
      ERR = 0.
      DO I = 1, LP1
         DO J = 1, MP1
            DO K = 1, NP1
               T = ABS(F(I,J,K)-X(I)**4*SIN(Y(J))*COS(Z(K)))
               ERR = AMAX1(T,ERR)
            END DO
         END DO
      END DO
!     Print earlier output from platforms with 32 and 64 bit floating point
!     arithemtic followed by the output from this computer
      WRITE (*, *) '    HW3CRT TEST RUN *** '
      WRITE (*, *) 
     1   '    Previous 64 bit floating point arithmetic result '
      WRITE (*, *) '    IERROR = 0,  Discretization Error = 9.6480E-3'
      WRITE (*, *) 
     1   '    Previous 32 bit floating point arithmetic result '
      WRITE (*, *) '    IERROR = 0,  Discretization Error = 9.6480E-3'
      WRITE (*, *) '    The output from your computer is: '
      WRITE (*, *) '    IERROR =', IERROR, ' Discretization Error = ', 
     1   ERR
      STOP 
      END PROGRAM THW3CRT
