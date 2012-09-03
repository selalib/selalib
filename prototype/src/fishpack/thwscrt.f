c
c     file thwscrt.f
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
      PROGRAM THWSCRT
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: IDIMF, M, MBDCND, N, NBDCND, MP1, NP1, I, J, IERROR
      REAL , DIMENSION(45,82) :: F
      REAL , DIMENSION(81) :: BDB, Y
      REAL , DIMENSION(41) :: X
      REAL::A,B,C,D,ELMBDA,PI,DUM,PIBY2,PISQ,BDA,BDC,BDD,PERTRB,ERR,Z
C-----------------------------------------------
C
C     FROM DIMENSION STATEMENT WE GET VALUE OF IDIMF.
C
      IDIMF = 45
      A = 0.
      B = 2.
      M = 40
      MBDCND = 2
      C = -1.
      D = 3.
      N = 80
      NBDCND = 0
      ELMBDA = -4.
C
C     AUXILIARY QUANTITIES.
C
      PI = 4.0*ATAN(1.0)
      PIBY2 = PI/2.
      PISQ = PI**2
      MP1 = M + 1
      NP1 = N + 1
C
C     GENERATE AND STORE GRID POINTS FOR THE PURPOSE OF COMPUTING
C     BOUNDARY DATA AND THE RIGHT SIDE OF THE HELMHOLTZ EQUATION.
C
      DO I = 1, MP1
         X(I) = FLOAT(I - 1)/20.
      END DO
      DO J = 1, NP1
         Y(J) = (-1.) + FLOAT(J - 1)/20.
      END DO
C
C     GENERATE BOUNDARY DATA.
C
      DO J = 1, NP1
         BDB(J) = 4.*COS((Y(J)+1.)*PIBY2)
      END DO
C
C     BDA, BDC, AND BDD ARE DUMMY VARIABLES.
C
      F(1,:NP1) = 0.
C
C     GENERATE RIGHT SIDE OF EQUATION.
C
      DO I = 2, MP1
         DO J = 1, NP1
            F(I,J) = (2. - (4. + PISQ/4.)*X(I)**2)*COS((Y(J)+1.)*PIBY2)
         END DO
      END DO
      CALL HWSCRT (A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND, BDC, BDD
     1   , ELMBDA, F, IDIMF, PERTRB, IERROR)
C
C     COMPUTE DISCRETIZATION ERROR.  THE EXACT SOLUTION IS
C                U(X,Y) = X**2*COS((Y+1)*PIBY2)
C
      ERR = 0.
      DO I = 1, MP1
         DO J = 1, NP1
            Z = ABS(F(I,J)-X(I)**2*COS((Y(J)+1.)*PIBY2))
            ERR = AMAX1(Z,ERR)
         END DO
      END DO
!     Print earlier output from platforms with 32 and 64 bit floating point
!     arithemtic followed by the output from this computer
      WRITE (*, *) '    HWSCRT TEST RUN *** '
      WRITE (*, *) 
     1   '    Previous 64 bit floating point arithmetic result '
      WRITE (*, *) '    IERROR = 0,  Discretization Error = 5.36508-4'
      WRITE (*, *) 
     1   '    Previous 32 bit floating point arithmetic result '
      WRITE (*, *) '    IERROR = 0,  Discretization Error = 4.9305E-4'
      WRITE (*, *) '    The output from your computer is: '
      WRITE (*, *) '    IERROR =', IERROR, ' Discretization Error = ', 
     1   ERR
      STOP 
      END PROGRAM THWSCRT
