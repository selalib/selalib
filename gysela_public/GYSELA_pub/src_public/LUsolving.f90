!**************************************************************
!  Copyright Euratom-CEA
!  Authors : 
!     Virginie Grandgirard (virginie.grandgirard@cea.fr)
!     Chantal Passeron (chantal.passeron@cea.fr)
!     Guillaume Latu (guillaume.latu@cea.fr)
!     Xavier Garbet (xavier.garbet@cea.fr)
!     Philippe Ghendrih (philippe.ghendrih@cea.fr)
!     Yanick Sarazin (yanick.sarazin@cea.fr)
!  
!  This code GYSELA (for GYrokinetic SEmi-LAgrangian) 
!  is a 5D gyrokinetic global full-f code for simulating 
!  the plasma turbulence in a tokamak.
!  
!  This software is governed by the CeCILL-B license 
!  under French law and abiding by the rules of distribution 
!  of free software.  You can  use, modify and redistribute 
!  the software under the terms of the CeCILL-B license as 
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info". 
!**************************************************************
     SUBROUTINE ZGTTRF( N, DL, D, DU, DU2, IPIV, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      INTEGER            INFO, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         D( * ), DL( * ), DU( * ), DU2( * )
!     ..
!
!  Purpose
!  =======
!
!  ZGTTRF computes an LU factorization of a complex 
!   tridiagonal matrix A using elimination with partial 
!   pivoting and row interchanges.
!
!  The factorization has the form
!     A = L * U
!  where L is a product of permutation and unit lower bidiagonal
!  matrices and U is upper triangular with nonzeros in only the main
!  diagonal and first two superdiagonals.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the matrix A.
!
!  DL      (input/output) COMPLEX*16 array, dimension (N-1)
!          On entry, DL must contain the (n-1) sub-diagonal 
!          elements of A.
!
!          On exit, DL is overwritten by the (n-1) multipliers that
!          define the matrix L from the LU factorization of A.
!
!  D       (input/output) COMPLEX*16 array, dimension (N)
!          On entry, D must contain the diagonal elements of A.
!
!          On exit, D is overwritten by the n diagonal 
!          elements of the upper triangular matrix U 
!          from the LU factorization of A.
!
!  DU      (input/output) COMPLEX*16 array, dimension (N-1)
!          On entry, DU must contain the (n-1) 
!          super-diagonal elements of A.
!
!          On exit, DU is overwritten by the (n-1) elements 
!          of the first super-diagonal of U.
!
!  DU2     (output) COMPLEX*16 array, dimension (N-2)
!          On exit, DU2 is overwritten by the (n-2) elements of the
!          second super-diagonal of U.
!
!  IPIV    (output) INTEGER array, dimension (N)
!          The pivot indices; for 1 <= i <= n, row i of 
!          the matrix was interchanged with row IPIV(i).  
!          IPIV(i) will always be either i or i+1; 
!          IPIV(i) = i indicates a row interchange was not
!          required.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -k, the k-th argument 
!          had an illegal value
!          > 0:  if INFO = k, U(k,k) is exactly zero. 
!          The factorization has been completed, 
!          but the factor U is exactly singular, and division 
!          by zero will occur if it is used to solve a system 
!          of equations.
!
!  ===============================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I
      COMPLEX*16         FACT, TEMP, ZDUM
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DIMAG
!     ..
!     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
!     ..
!     .. Statement Function definitions ..
      CABS1( ZDUM ) = ABS( DBLE( ZDUM ) ) + ABS( DIMAG( ZDUM ) )
!     ..
!     .. Executable Statements ..
!
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
         WRITE(6,*)'ZGTTRF', -INFO
         STOP
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) &
        RETURN
!
!     Initialize IPIV(i) = i and DU2(i) = 0
!
      DO 10 I = 1, N
         IPIV( I ) = I
   10 CONTINUE
      DO 20 I = 1, N - 2
         DU2( I ) = ZERO
   20 CONTINUE
!
      DO 30 I = 1, N - 2
         IF( CABS1( D( I ) ).GE.CABS1( DL( I ) ) ) THEN
!
!           No row interchange required, eliminate DL(I)
!
            IF( CABS1( D( I ) ).NE.ZERO ) THEN
               FACT = DL( I ) / D( I )
               DL( I ) = FACT
               D( I+1 ) = D( I+1 ) - FACT*DU( I )
            END IF
         ELSE
!
!           Interchange rows I and I+1, eliminate DL(I)
!
            FACT = D( I ) / DL( I )
            D( I ) = DL( I )
            DL( I ) = FACT
            TEMP = DU( I )
            DU( I ) = D( I+1 )
            D( I+1 ) = TEMP - FACT*D( I+1 )
            DU2( I ) = DU( I+1 )
            DU( I+1 ) = -FACT*DU( I+1 )
            IPIV( I ) = I + 1
         END IF
   30 CONTINUE
      IF( N.GT.1 ) THEN
         I = N - 1
         IF( CABS1( D( I ) ).GE.CABS1( DL( I ) ) ) THEN
            IF( CABS1( D( I ) ).NE.ZERO ) THEN
               FACT = DL( I ) / D( I )
               DL( I ) = FACT
               D( I+1 ) = D( I+1 ) - FACT*DU( I )
            END IF
         ELSE
            FACT = D( I ) / DL( I )
            D( I ) = DL( I )
            DL( I ) = FACT
            TEMP = DU( I )
            DU( I ) = D( I+1 )
            D( I+1 ) = TEMP - FACT*D( I+1 )
            IPIV( I ) = I + 1
         END IF
      END IF
!
!     Check for a zero on the diagonal of U.
!
      DO 40 I = 1, N
         IF( CABS1( D( I ) ).EQ.ZERO ) THEN
            INFO = I
            GO TO 50
         END IF
   40 CONTINUE
   50 CONTINUE
!
      RETURN
!
!     End of ZGTTRF
!
      END
      
      SUBROUTINE ZGTTRS( TRANS, N, NRHS, DL, D, DU, DU2, IPIV, B, LDB, &
                        INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         B( LDB, * ), D( * ), DL( * ), DU( * ), DU2( * )
!     ..
!
!  Purpose
!  =======
!
!  ZGTTRS solves one of the systems of equations
!     A * X = B,  A**T * X = B,  or  A**H * X = B,
!  with a tridiagonal matrix A using the LU factorization computed
!  by ZGTTRF.
!
!  Arguments
!  =========
!
!  TRANS   (input) CHARACTER
!          Specifies the form of the system of equations.
!          = 'N':  A * X = B     (No transpose)
!          = 'T':  A**T * X = B  (Transpose)
!          = 'C':  A**H * X = B  (Conjugate transpose)
!
!  N       (input) INTEGER
!          The order of the matrix A.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  DL      (input) COMPLEX*16 array, dimension (N-1)
!          The (n-1) multipliers that define the matrix L from the
!          LU factorization of A.
!
!  D       (input) COMPLEX*16 array, dimension (N)
!          The n diagonal elements of the upper triangular matrix U from
!          the LU factorization of A.
!
!  DU      (input) COMPLEX*16 array, dimension (N-1)
!          The (n-1) elements of the first super-diagonal of U.
!
!  DU2     (input) COMPLEX*16 array, dimension (N-2)
!          The (n-2) elements of the second super-diagonal of U.
!
!  IPIV    (input) INTEGER array, dimension (N)
!          The pivot indices; for 1 <= i <= n, row i of the matrix was
!          interchanged with row IPIV(i).  IPIV(i) will always be either
!          i or i+1; IPIV(i) = i indicates a row interchange was not
!          required.
!
!  B       (input/output) COMPLEX*16 array, dimension (LDB,NRHS)
!          On entry, the matrix of right hand side vectors B.
!          On exit, B is overwritten by the solution vectors X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -k, the k-th argument had an illegal value
!
!  =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            NOTRAN
      INTEGER            ITRANS, J, JB, NB
!     ..
!     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA, ZGTTS2
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
      INFO = 0
      NOTRAN = ( TRANS.EQ.'N' .OR. TRANS.EQ.'n' )
      IF( .NOT.NOTRAN .AND. .NOT.( TRANS.EQ.'T' .OR. TRANS.EQ. &
         't' ) .AND. .NOT.( TRANS.EQ.'C' .OR. TRANS.EQ.'c' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDB.LT.MAX( N, 1 ) ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         WRITE(6,*)'ZGTTRF', INFO
         STOP
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 .OR. NRHS.EQ.0 ) &
        RETURN
!
!     Decode TRANS
!
      IF( NOTRAN ) THEN
         ITRANS = 0
      ELSE IF( TRANS.EQ.'T' .OR. TRANS.EQ.'t' ) THEN
         ITRANS = 1
      ELSE
         ITRANS = 2
      END IF
!
!     Determine the number of right-hand sides to solve at a time.
!
      IF( NRHS.EQ.1 ) THEN
         NB = 1
      ELSE
         !NB = MAX(1,ILAENV( 1,'ZGTTRS',TRANS,N,NRHS,-1,-1))
         NB = 1
      END IF
!
      IF( NB.GE.NRHS ) THEN
         CALL ZGTTS2( ITRANS, N, NRHS, DL, D, DU, DU2, IPIV, &
           B, LDB )
      ELSE
         DO 10 J = 1, NRHS, NB
            JB = MIN( NRHS-J+1, NB )
            CALL ZGTTS2( ITRANS, N, JB, DL, D, DU, DU2, IPIV, &
              B( 1, J ), LDB )
   10    CONTINUE
      END IF
!
!     End of ZGTTRS
!
      END
      
      SUBROUTINE ZGTTS2( ITRANS, N, NRHS, DL, D, DU, &
        DU2, IPIV, B, LDB )
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      INTEGER            ITRANS, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         B( LDB, * ), D( * ), DL( * )
      COMPLEX*16         DU( * ), DU2( * )
!     ..
!
!  Purpose
!  =======
!
!  ZGTTS2 solves one of the systems of equations
!     A * X = B,  A**T * X = B,  or  A**H * X = B,
!  with a tridiagonal matrix A using the LU factorization computed
!  by ZGTTRF.
!
!  Arguments
!  =========
!
!  ITRANS  (input) INTEGER
!          Specifies the form of the system of equations.
!          = 0:  A * X = B     (No transpose)
!          = 1:  A**T * X = B  (Transpose)
!          = 2:  A**H * X = B  (Conjugate transpose)
!
!  N       (input) INTEGER
!          The order of the matrix A.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., 
!          the number of columns of the matrix B.  NRHS >= 0.
!
!  DL      (input) COMPLEX*16 array, dimension (N-1)
!          The (n-1) multipliers that define the matrix L from the
!          LU factorization of A.
!
!  D       (input) COMPLEX*16 array, dimension (N)
!          The n diagonal elements of the upper triangular 
!          matrix U from the LU factorization of A.
!
!  DU      (input) COMPLEX*16 array, dimension (N-1)
!          The (n-1) elements of the first super-diagonal of U.
!
!  DU2     (input) COMPLEX*16 array, dimension (N-2)
!          The (n-2) elements of the second super-diagonal of U.
!
!  IPIV    (input) INTEGER array, dimension (N)
!          The pivot indices; for 1 <= i <= n, 
!          row i of the matrix was interchanged with row IPIV(i).
!          IPIV(i) will always be either i or i+1; 
!          IPIV(i) = i indicates a row interchange was not
!          required.
!
!  B       (input/output) COMPLEX*16 array, dimension (LDB,NRHS)
!          On entry, the matrix of right hand side vectors B.
!          On exit, B is overwritten by the solution vectors X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  ===============================================================
!
!     .. Local Scalars ..
      INTEGER            I, J
      COMPLEX*16         TEMP
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DCONJG
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF( N.EQ.0 .OR. NRHS.EQ.0 ) &
        RETURN
!
      IF( ITRANS.EQ.0 ) THEN
!
!        Solve A*X = B using the LU factorization of A,
!        overwriting each right hand side vector with its solution.
!
         IF( NRHS.LE.1 ) THEN
            J = 1
   10       CONTINUE
!
!           Solve L*x = b.
!
            DO 20 I = 1, N - 1
               IF( IPIV( I ).EQ.I ) THEN
                  B( I+1, J ) = B( I+1, J ) - DL( I )*B( I, J )
               ELSE
                  TEMP = B( I, J )
                  B( I, J ) = B( I+1, J )
                  B( I+1, J ) = TEMP - DL( I )*B( I, J )
               END IF
   20       CONTINUE
!
!           Solve U*x = b.
!
            B( N, J ) = B( N, J ) / D( N )
            IF( N.GT.1 ) &
              B( N-1, J ) = ( B( N-1, J )-DU( N-1 )*B( N, J ) ) / &
                            D( N-1 )
            DO 30 I = N - 2, 1, -1
               B( I, J ) = ( B( I, J )-DU( I )* &
                 B( I+1, J )-DU2( I )*B( I+2, J ) ) / D( I )
   30       CONTINUE
            IF( J.LT.NRHS ) THEN
               J = J + 1
               GO TO 10
            END IF
         ELSE
            DO 60 J = 1, NRHS
!
!           Solve L*x = b.
!
               DO 40 I = 1, N - 1
                  IF( IPIV( I ).EQ.I ) THEN
                     B( I+1, J ) = B( I+1, J ) - DL( I )*B( I, J )
                  ELSE
                     TEMP = B( I, J )
                     B( I, J ) = B( I+1, J )
                     B( I+1, J ) = TEMP - DL( I )*B( I, J )
                  END IF
   40          CONTINUE
!
!           Solve U*x = b.
!
               B( N, J ) = B( N, J ) / D( N )
               IF( N.GT.1 ) &
                 B( N-1, J ) = ( B( N-1, J )-DU( N-1 ) * &
                 B( N, J ) ) / D( N-1 )
               DO 50 I = N - 2, 1, -1
                  B( I, J ) = ( B( I, J )-DU( I ) * &
                    B( I+1, J )-DU2( I )*B( I+2, J ) ) / D( I )
   50          CONTINUE
   60       CONTINUE
         END IF
      ELSE IF( ITRANS.EQ.1 ) THEN
!
!        Solve A**T * X = B.
!
         IF( NRHS.LE.1 ) THEN
            J = 1
   70       CONTINUE
!
!           Solve U**T * x = b.
!
            B( 1, J ) = B( 1, J ) / D( 1 )
            IF( N.GT.1 ) &
              B( 2, J ) = ( B( 2, J )-DU( 1 )*B( 1, J ) ) / D( 2 )
            DO 80 I = 3, N
               B( I, J ) = ( B( I, J )-DU( I-1 ) * &
                 B( I-1, J )-DU2( I-2 )* B( I-2, J ) ) / D( I )
   80       CONTINUE
!
!           Solve L**T * x = b.
!
            DO 90 I = N - 1, 1, -1
               IF( IPIV( I ).EQ.I ) THEN
                  B( I, J ) = B( I, J ) - DL( I )*B( I+1, J )
               ELSE
                  TEMP = B( I+1, J )
                  B( I+1, J ) = B( I, J ) - DL( I )*TEMP
                  B( I, J ) = TEMP
               END IF
   90       CONTINUE
            IF( J.LT.NRHS ) THEN
               J = J + 1
               GO TO 70
            END IF
         ELSE
            DO 120 J = 1, NRHS
!
!           Solve U**T * x = b.
!
               B( 1, J ) = B( 1, J ) / D( 1 )
               IF( N.GT.1 ) &
                 B( 2, J ) = ( B( 2, J )-DU( 1 ) * &
                 B( 1, J ) ) / D( 2 )
               DO 100 I = 3, N
                  B( I, J ) = ( B( I, J )-DU( I-1 ) * &
                    B( I-1, J )-DU2( I-2 )*B( I-2, J ) ) / D( I )
  100          CONTINUE
!
!           Solve L**T * x = b.
!
               DO 110 I = N - 1, 1, -1
                  IF( IPIV( I ).EQ.I ) THEN
                     B( I, J ) = B( I, J ) - DL( I )*B( I+1, J )
                  ELSE
                     TEMP = B( I+1, J )
                     B( I+1, J ) = B( I, J ) - DL( I )*TEMP
                     B( I, J ) = TEMP
                  END IF
  110          CONTINUE
  120       CONTINUE
         END IF
      ELSE
!
!        Solve A**H * X = B.
!
         IF( NRHS.LE.1 ) THEN
            J = 1
  130       CONTINUE
!
!           Solve U**H * x = b.
!
            B( 1, J ) = B( 1, J ) / DCONJG( D( 1 ) )
            IF( N.GT.1 ) &
              B( 2, J ) = ( B( 2, J )-DCONJG( DU( 1 ) ) * &
              B( 1, J ) ) / DCONJG( D( 2 ) )
            DO 140 I = 3, N
               B( I, J ) = ( B( I, J )-DCONJG( DU( I-1 ) ) * &
                 B( I-1, J )- DCONJG( DU2( I-2 ) ) * &
                 B( I-2, J ) ) / DCONJG( D( I ) )
  140       CONTINUE
!
!           Solve L**H * x = b.
!
            DO 150 I = N - 1, 1, -1
               IF( IPIV( I ).EQ.I ) THEN
                  B( I, J ) = B( I, J ) - &
                    DCONJG( DL( I ) )*B( I+1, J )
               ELSE
                  TEMP = B( I+1, J )
                  B( I+1, J ) = B( I, J ) - DCONJG( DL( I ) )*TEMP
                  B( I, J ) = TEMP
               END IF
  150       CONTINUE
            IF( J.LT.NRHS ) THEN
               J = J + 1
               GO TO 130
            END IF
         ELSE
            DO 180 J = 1, NRHS
!
!           Solve U**H * x = b.
!
               B( 1, J ) = B( 1, J ) / DCONJG( D( 1 ) )
               IF( N.GT.1 ) &
                 B( 2, J ) = ( B( 2, J ) - &
                 DCONJG( DU( 1 ) )*B( 1, J ) ) / DCONJG( D( 2 ) )
               DO 160 I = 3, N
                  B( I, J ) = ( B( I, J )-DCONJG( DU( I-1 ) ) * &
                             B( I-1, J )-DCONJG( DU2( I-2 ) ) * &
                             B( I-2, J ) ) / DCONJG( D( I ) )
  160          CONTINUE
!
!           Solve L**H * x = b.
!
               DO 170 I = N - 1, 1, -1
                  IF( IPIV( I ).EQ.I ) THEN
                     B( I, J ) = B( I, J ) - DCONJG( DL( I ) ) * &
                                B( I+1, J )
                  ELSE
                     TEMP = B( I+1, J )
                     B( I+1, J ) = B( I, J ) - &
                       DCONJG( DL( I ) )*TEMP
                     B( I, J ) = TEMP
                  END IF
  170          CONTINUE
  180       CONTINUE
         END IF
      END IF
!
!     End of ZGTTS2
!
      END
      
      SUBROUTINE DPTTRF( N, D, E, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      INTEGER            INFO, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * )
!     ..
!
!  Purpose
!  =======
!
!  DPTTRF computes the L*D*L' factorization of a real symmetric
!  positive definite tridiagonal matrix A.  
!  The factorization may also be regarded as having 
!  the form A = U'*D*U.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  D       (input/output) DOUBLE PRECISION array, dimension (N)
!          On entry, the n diagonal elements of the 
!          tridiagonal matrix A.  On exit, the n diagonal 
!          elements of the diagonal matrix D from the 
!          L*D*L' factorization of A.
!
!  E       (input/output) DOUBLE PRECISION array, dimension (N-1)
!          On entry, the (n-1) subdiagonal elements 
!          of the tridiagonal matrix A.  
!          On exit, the (n-1) subdiagonal elements of the
!          unit bidiagonal factor L from the L*D*L' 
!          factorization of A.
!          E can also be regarded as the superdiagonal of the unit
!          bidiagonal factor U from the U'*D*U factorization of A.
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had 
!          an illegal value
!          > 0: if INFO = k, the leading minor of order k is not
!               positive definite; if k < N, 
!               the factorization could not be completed, 
!               while if k = N, the factorization was
!               completed, but D(N) = 0.
!
!  ===============================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, I4
      DOUBLE PRECISION   EI
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MOD
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
         WRITE(6,*)'DPTTRF', -INFO
         STOP
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) &
        RETURN
!
!     Compute the L*D*L' (or U'*D*U) factorization of A.
!
      I4 = MOD( N-1, 4 )
      DO 10 I = 1, I4
         IF( D( I ).LE.ZERO ) THEN
            INFO = I
            GO TO 30
         END IF
         EI = E( I )
         E( I ) = EI / D( I )
         D( I+1 ) = D( I+1 ) - E( I )*EI
   10 CONTINUE
!
      DO 20 I = I4 + 1, N - 4, 4
!
!        Drop out of the loop if d(i) <= 0: the matrix is 
!        not positive definite.
!
         IF( D( I ).LE.ZERO ) THEN
            INFO = I
            GO TO 30
         END IF
!
!        Solve for e(i) and d(i+1).
!
         EI = E( I )
         E( I ) = EI / D( I )
         D( I+1 ) = D( I+1 ) - E( I )*EI
!
         IF( D( I+1 ).LE.ZERO ) THEN
            INFO = I + 1
            GO TO 30
         END IF
!
!        Solve for e(i+1) and d(i+2).
!
         EI = E( I+1 )
         E( I+1 ) = EI / D( I+1 )
         D( I+2 ) = D( I+2 ) - E( I+1 )*EI
!
         IF( D( I+2 ).LE.ZERO ) THEN
            INFO = I + 2
            GO TO 30
         END IF
!
!        Solve for e(i+2) and d(i+3).
!
         EI = E( I+2 )
         E( I+2 ) = EI / D( I+2 )
         D( I+3 ) = D( I+3 ) - E( I+2 )*EI
!
         IF( D( I+3 ).LE.ZERO ) THEN
            INFO = I + 3
            GO TO 30
         END IF
!
!        Solve for e(i+3) and d(i+4).
!
         EI = E( I+3 )
         E( I+3 ) = EI / D( I+3 )
         D( I+4 ) = D( I+4 ) - E( I+3 )*EI
   20 CONTINUE
!
!     Check d(n) for positive definiteness.
!
      IF( D( N ).LE.ZERO ) &
        INFO = N
!
   30 CONTINUE
      RETURN
!
!     End of DPTTRF
!
      END
      
      SUBROUTINE DPTTRS( N, NRHS, D, E, B, LDB, INFO )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   B( LDB, * ), D( * ), E( * )
!     ..
!
!  Purpose
!  =======
!
!  DPTTRS solves a tridiagonal system of the form
!     A * X = B
!  using the L*D*L' factorization of A computed by DPTTRF.  D is a
!  diagonal matrix specified in the vector D, L is a unit 
!  bidiagonal matrix whose subdiagonal is specified in the 
!  vector E, and X and B are N by NRHS matrices.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the tridiagonal matrix A.  N >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., 
!          the number of columns of the matrix B.  NRHS >= 0.
!
!  D       (input) DOUBLE PRECISION array, dimension (N)
!          The n diagonal elements of the diagonal matrix D 
!          from the L*D*L' factorization of A.
!
!  E       (input) DOUBLE PRECISION array, dimension (N-1)
!          The (n-1) subdiagonal elements of the unit 
!          bidiagonal factor L from the L*D*L' factorization of A.
!          E can also be regarded as the superdiagonal of 
!          the unit bidiagonal factor U from the factorization 
!          A = U'*D*U.
!
!  B       (input/output) DOUBLE PRECISION array, 
!          dimension (LDB,NRHS)
!          On entry, the right hand side vectors B for the 
!          system of linear equations.
!          On exit, the solution vectors, X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!
!  ================================================================
!
!     .. Local Scalars ..
      INTEGER            J, JB, NB
!     ..
!     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
!     ..
!     .. External Subroutines ..
      EXTERNAL           DPTTS2, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments.
!
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         WRITE(6,*)'DPTTRS', -INFO
         STOP
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 .OR. NRHS.EQ.0 ) &
        RETURN
!
!     Determine the number of right-hand sides to solve at a time.
!
      IF( NRHS.EQ.1 ) THEN
         NB = 1
      ELSE
         !NB = MAX( 1, ILAENV( 1, 'DPTTRS', ' ', N, NRHS, -1, -1 ) )
         NB = 1
      END IF
!
      IF( NB.GE.NRHS ) THEN
         CALL DPTTS2( N, NRHS, D, E, B, LDB )
      ELSE
         DO 10 J = 1, NRHS, NB
            JB = MIN( NRHS-J+1, NB )
            CALL DPTTS2( N, JB, D, E, B( 1, J ), LDB )
   10    CONTINUE
      END IF
!
      RETURN
!
!     End of DPTTRS
!
      END
      
      SUBROUTINE DPTTS2( N, NRHS, D, E, B, LDB )
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!
!     .. Scalar Arguments ..
      INTEGER            LDB, N, NRHS
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   B( LDB, * ), D( * ), E( * )
!     ..
!
!  Purpose
!  =======
!
!  DPTTS2 solves a tridiagonal system of the form
!     A * X = B
!  using the L*D*L' factorization of A computed by DPTTRF.  D is a
!  diagonal matrix specified in the vector D, L is a unit bidiagonal
!  matrix whose subdiagonal is specified in the vector E, 
!  and X and B are N by NRHS matrices.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the tridiagonal matrix A.  N >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., 
!          the number of columns of the matrix B.  NRHS >= 0.
!
!  D       (input) DOUBLE PRECISION array, dimension (N)
!          The n diagonal elements of the diagonal matrix D 
!          from the L*D*L' factorization of A.
!
!  E       (input) DOUBLE PRECISION array, dimension (N-1)
!          The (n-1) subdiagonal elements of the unit 
!          bidiagonal factor L from the L*D*L' factorization of A.
!          E can also be regarded as the superdiagonal of 
!          the unit bidiagonal factor U from the factorization 
!          A = U'*D*U.
!
!  B       (input/output) DOUBLE PRECISION array, 
!          dimension (LDB,NRHS)
!          On entry, the right hand side vectors B for 
!          the system of linear equations.
!          On exit, the solution vectors, X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  ================================================================
!
!     .. Local Scalars ..
      INTEGER            I, J
!     ..
!     .. External Subroutines ..
      EXTERNAL           DSCAL
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF( N.LE.1 ) THEN
         IF( N.EQ.1 ) &
           CALL DSCAL( NRHS, 1.D0 / D( 1 ), B, LDB )
         RETURN
      END IF
!
!     Solve A * X = B using the factorization A = L*D*L',
!     overwriting each right hand side vector with its solution.
!
      DO 30 J = 1, NRHS
!
!           Solve L * x = b.
!
         DO 10 I = 2, N
            B( I, J ) = B( I, J ) - B( I-1, J )*E( I-1 )
   10    CONTINUE
!
!           Solve D * L' * x = b.
!
         B( N, J ) = B( N, J ) / D( N )
         DO 20 I = N - 1, 1, -1
            B( I, J ) = B( I, J ) / D( I ) - B( I+1, J )*E( I )
   20    CONTINUE
   30 CONTINUE
!
      RETURN
!
!     End of DPTTS2
!
      END
      
      subroutine  dscal(n,da,dx,incx)
!
!     scales a vector by a constant.
!     uses unrolled loops for increment equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      double precision da,dx(*)
      integer i,incx,m,mp1,n,nincx
!
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
!
!        code for increment not equal to 1
!
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dx(i) = da*dx(i)
   10 continue
      return
!
!        code for increment equal to 1
!
!
!        clean-up loop
!
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dx(i) = da*dx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 continue
      return
      end
