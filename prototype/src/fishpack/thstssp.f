c
c     file thstssp.f
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
      PROGRAM THSTSSP
      implicit none
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: M, MBDCND, N, NBDCND, IDIMF, I, J, IERROR
      REAL , DIMENSION(18,72) :: F
      REAL , DIMENSION(72) :: BDB
      REAL , DIMENSION(18) :: SINT
      REAL , DIMENSION(72) :: SINP
      REAL::PI,A,B,C,D,ELMBDA,DTHETA,DPHI,BDA,BDC,BDD,PERTRB,ERR,Z
C-----------------------------------------------
C
C     THE VALUE OF IDIMF IS THE FIRST DIMENSION OF F.
C
      PI = 4.0*ATAN(1.0)
      A = 0.
      B = PI/2.
      M = 18
      MBDCND = 6
      C = 0.
      D = 2.*PI
      N = 72
      NBDCND = 0
      ELMBDA = 0.
      IDIMF = 18
C
C     GENERATE SINES FOR USE IN SUBSEQUENT COMPUTATIONS
C
      DTHETA = B/FLOAT(M)
      DO I = 1, M
         SINT(I) = SIN((FLOAT(I) - 0.5)*DTHETA)
      END DO
      DPHI = D/FLOAT(N)
      DO J = 1, N
         SINP(J) = SIN((FLOAT(J) - 0.5)*DPHI)
      END DO
C
C     COMPUTE RIGHT SIDE OF EQUATION AND STORE IN F
C
      DO J = 1, N
         F(:M,J) = 2. - 6.*(SINT(:M)*SINP(J))**2
      END DO
C
C     STORE DERIVATIVE DATA AT THE EQUATOR
C
      BDB(:N) = 0.
C
C     BDA, BDC, AND BDD ARE DUMMY VARIABLES.
C
      CALL HSTSSP (A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND, BDC, BDD
     1   , ELMBDA, F, IDIMF, PERTRB, IERROR)
C
C     COMPUTE DISCRETIZATION ERROR. SINCE PROBLEM IS SINGULAR, THE
C     SOLUTION MUST BE NORMALIZED.
C
      ERR = 0.
      DO J = 1, N
         DO I = 1, M
            Z = ABS(F(I,J)-(SINT(I)*SINP(J))**2-F(1,1))
            ERR = AMAX1(Z,ERR)
         END DO
      END DO
!     Print earlier output from platforms with 32 and 64 bit floating point
!     arithemtic followed by the output from this computer
      WRITE (*, *) '    HSTSSP TEST RUN *** '
      WRITE (*, *) 
     1   '    Previous 64 bit floating point arithmetic result '
      WRITE (*, *) '    IERROR = 0,  PERTRB = 6.35830E-4'
      WRITE (*, *) '    discretization error = 3.37523E-3'
      WRITE (*, *) 
     1   '    Previous 32 bit floating point arithmetic result '
      WRITE (*, *) '    IERROR = 0,  PERTRB = 6.35919E-4'
      WRITE (*, *) '    discretization error = 3.38144E-3'
      WRITE (*, *) '    The output from your computer is: '
      WRITE (*, *) '    IERROR =', IERROR, ' PERTRB = ', PERTRB
      WRITE (*, *) '    discretization error = ', ERR
      STOP 
      END PROGRAM THSTSSP
