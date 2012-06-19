PROGRAM POISSON
IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
INTEGER :: IDIMF, M, MBDCND, N, NBDCND, MP1, NP1, I, J, IERROR
REAL(8) , DIMENSION(45,82) :: F
REAL(8) , DIMENSION(81) :: BDB, Y
REAL(8) , DIMENSION(41) :: X
REAL(8)::A,B,C,D,ELMBDA,PI,DUM,PIBY2,PISQ,BDA,BDC,BDD,PERTRB,ERR,Z
!-----------------------------------------------
!
!     FROM DIMENSION STATEMENT WE GET VALUE OF IDIMF.
!
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

!     AUXILIARY QUANTITIES.

PI = 4.0*ATAN(1.0)
PIBY2 = PI/2.
PISQ = PI**2
MP1 = M + 1
NP1 = N + 1

!     GENERATE AND STORE GRID POINTS FOR THE PURPOSE OF COMPUTING
!     BOUNDARY DATA AND THE RIGHT SIDE OF THE HELMHOLTZ EQUATION.

DO I = 1, MP1
   X(I) = A+(B-A)*FLOAT(I-1)/M
END DO

DO J = 1, NP1
   Y(J) = C+(D-C)*FLOAT(J-1)/N
END DO

!     GENERATE BOUNDARY DATA.

DO J = 1, NP1
   BDB(J) = 4.*COS((Y(J)+1.)*PIBY2)
END DO

!     BDA, BDC, AND BDD ARE DUMMY VARIABLES.

F(1,:NP1) = 0.

!     GENERATE RIGHT SIDE OF EQUATION.

DO I = 2, MP1
   DO J = 1, NP1
      F(I,J) = (2. - (4. + PISQ/4.)*X(I)**2)*COS((Y(J)+1.)*PIBY2)
   END DO
END DO
CALL HWSCRT (A, B, M, MBDCND, BDA, BDB, C, D, N, NBDCND, BDC, BDD &
      , ELMBDA, F, IDIMF, PERTRB, IERROR)

!     COMPUTE DISCRETIZATION ERROR.  THE EXACT SOLUTION IS
!                U(X,Y) = X**2*COS((Y+1)*PIBY2)

ERR = 0.
DO I = 1, MP1
   DO J = 1, NP1
      Z = ABS(F(I,J)-X(I)**2*COS((Y(J)+1.)*PIBY2))
      ERR = MAX(Z,ERR)
   END DO
END DO

!     Print earlier output from platforms with 32 and 64 bit floating point
!     arithemtic followed by the output from this computer

WRITE (*, *) '    HWSCRT TEST RUN *** '
WRITE (*, *) '    Previous 64 bit floating point arithmetic result '
WRITE (*, *) '    IERROR = 0,  Discretization Error = 5.36508-4'
WRITE (*, *) '    Previous 32 bit floating point arithmetic result '
WRITE (*, *) '    IERROR = 0,  Discretization Error = 4.9305E-4'
WRITE (*, *) '    The output from your computer is: '
WRITE (*, *) '    IERROR =', IERROR, ' Discretization Error = ', ERR
STOP 
END PROGRAM POISSON
