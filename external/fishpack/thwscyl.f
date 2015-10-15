c
c     file thwscyl.f
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
      PROGRAM THWSCYL
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: IDIMF, M, MBDCND, N, NBDCND, MP1, NP1, I, J, IERROR
      REAL , DIMENSION(75,105) :: F
      REAL , DIMENSION(101) :: BDA, BDB
      REAL , DIMENSION(51) :: BDC, BDD, R
      REAL , DIMENSION(101) :: Z
      REAL :: A, B, C, D, ELMBDA, PERTRB, X, ERR
C-----------------------------------------------
C
C     FROM DIMENSION STATEMENT WE GET VALUE OF IDIMF.
C
      IDIMF = 75
      A = 0.
      B = 1.
      M = 50
      MBDCND = 6
      C = 0.
      D = 1.
      N = 100
      NBDCND = 3
      ELMBDA = 0.
C
C     AUXILIARY QUANTITIES.
C
      MP1 = M + 1
      NP1 = N + 1
C
C     GENERATE AND STORE GRID POINTS FOR THE PURPOSE OF COMPUTING
C     BOUNDARY DATA AND THE RIGHT SIDE OF THE POISSON EQUATION.
C
      DO I = 1, MP1
         R(I) = FLOAT(I - 1)/50.
      END DO
      DO J = 1, NP1
         Z(J) = FLOAT(J - 1)/100.
      END DO
C
C     GENERATE BOUNDARY DATA.
C
      BDB(:NP1) = 4.*Z(:NP1)**4
C
C     GENERATE BOUNDARY DATA.
C
      BDC(:MP1) = 0.
      BDD(:MP1) = 4.*R(:MP1)**4
C
C     BDA IS A DUMMY VARIABLE.
C
C
C     GENERATE RIGHT SIDE OF EQUATION.
C
      DO I = 1, MP1
         F(I,:NP1) = 4.*R(I)**2*Z(:NP1)**2*(4.*Z(:NP1)**2+3.*R(I)**2)
      END DO
      CALL HWSCYL (A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND, BDC, BDD
     1   , ELMBDA, F, IDIMF, PERTRB, IERROR)
C
C     COMPUTE DISCRETIZATION ERROR BY MINIMIZING OVER ALL A THE FUNCTION
C     NORM(F(I,J) - A*1 - U(R(I),Z(J))).  THE EXACT SOLUTION IS
C                U(R,Z) = (R*Z)**4 + ARBITRARY CONSTANT.
C
      X = 0.
      DO I = 1, MP1
         X = X + SUM(F(I,:NP1)-(R(I)*Z(:NP1))**4)
      END DO
      X = X/FLOAT(NP1*MP1)
      F(:MP1,:NP1) = F(:MP1,:NP1) - X
      ERR = 0.
      DO I = 1, MP1
         DO J = 1, NP1
            X = ABS(F(I,J)-(R(I)*Z(J))**4)
            ERR = AMAX1(X,ERR)
         END DO
      END DO
!     Print earlier output from platforms with 32 and 64 bit floating point
!     arithemtic followed by the output from this computer
      WRITE (*, *) '    HWSCYL TEST RUN *** '
      WRITE (*, *) 
     1   '    Previous 64 bit floating point arithmetic result '
      WRITE (*, *) '    IERROR = 0,  PERTRB = 2.2674E-4'
      WRITE (*, *) '    Discretization Error = 3.7367E-4 '
      WRITE (*, *) 
     1   '    Previous 32 bit floating point arithmetic result '
      WRITE (*, *) '    IERROR = 0,  PERTRB = 2.26976-4'
      WRITE (*, *) '    Discretization Error = 3.5554E-4'
      WRITE (*, *) '    The output from your computer is: '
      WRITE (*, *) '    IERROR =', IERROR, ' PERTRB = ', PERTRB
      WRITE (*, *) '    Discretization Error = ', ERR
      STOP 
      END PROGRAM THWSCYL
