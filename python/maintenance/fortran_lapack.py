data_types = [ \
'S',  # real
'D',  # double precision
'C',  # complex
'Z',  # complex*16 or double complex
]

matrix_types = [ \
'BD',  # bidiagonal
'DI',  # diagonal
'GB',  # general band
'GE',  # general (i.e., unsymmetric, in some cases rectangular)
'GG',  # general matrices, generalized problem (i.e., a pair of general matrices)
'GT',  # general tridiagonal
'HB',  # (complex) Hermitian band
'HE',  # (complex) Hermitian
'HG',  # upper Hessenberg matrix, generalized problem (i.e a Hessenberg and a
       # triangular matrix)
'HP',  # (complex) Hermitian, packed storage
'HS',  # upper Hessenberg
'OP',  # (real) orthogonal, packed storage
'OR',  # (real) orthogonal
'PB',  # symmetric or Hermitian positive definite band
'PO',  # symmetric or Hermitian positive definite
'PP',  # symmetric or Hermitian positive definite, packed storage
'PT',  # symmetric or Hermitian positive definite tridiagonal
'SB',  # (real) symmetric band
'SP',  # symmetric, packed storage
'ST',  # (real) symmetric tridiagonal
'SY',  # symmetric
'TB',  # triangular band
'TG',  # triangular matrices, generalized problem (i.e., a pair of triangular matrices)
'TP',  # triangular, packed storage
'TR',  # triangular (or in some cases quasi-triangular)
'TZ',  # trapezoidal
'UN',  # (complex) unitary
'UP',  # (complex) unitary, packed storage
]

computations = [ \
'BAK',  # back transformation of eigenvectors after balancing 
'BAL',  # permute and/or balance to isolate eigenvalues 
'BRD',  # reduce to bidiagonal form by orthogonal transformations 
'CON',  # estimate condition number
'EBZ',  # compute selected eigenvalues by bisection 
'EDC',  # compute eigenvectors using divide and conquer 
'EIN',  # compute selected eigenvectors by inverse iteration 
'EQR',  # compute eigenvalues and/or the Schur form using the QR algorithm 
'EGR',  # compute selected eigenvalues, and optionally, eigenvectors 
        # using Relatively Robust Representations 
'EQU',  # equilibrate a matrix to reduce its condition number 
'EQZ',  # compute generalized eigenvalues and/or generalized Schur form by QZ method 
'ERF',  # compute eigenvectors using the Pal-Walker-Kahan variant of the QL or QR 
        # algorithm 
'EVC',  # compute eigenvectors from Schur factorization 
'EXC',  # swap adjacent diagonal blocks in a quasi-upper triangular matrix 
'GBR',  # generate the orthogonal/unitary matrix from xGEBRD 
'GHR',  # generate the orthogonal/unitary matrix from xGEHRD 
'GLQ',  # generate the orthogonal/unitary matrix from xGELQF 
'GQL',  # generate the orthogonal/unitary matrix from xGEQLF 
'GQR',  # generate the orthogonal/unitary matrix from xGEQRF 
'GRQ',  # generate the orthogonal/unitary matrix from xGERQF 
'GST',  # reduce a symmetric-definite generalized eigenvalue problem to standard form 
'GTR',  # generate the orthogonal/unitary matrix from xxxTRD 
'HRD',  # reduce to upper Hessenberg form by orthogonal transformations 
'LQF',  # compute an LQ factorization without pivoting 
'MBR',  # multiply by the orthogonal/unitary matrix from xGEBRD 
'MHR',  # multiply by the orthogonal/unitary matrix from xGEHRD 
'MLQ',  # multiply by the orthogonal/unitary matrix from xGELQF 
'MQL',  # multiply by the orthogonal/unitary matrix from xGEQLF 
'MQR',  # multiply by the orthogonal/unitary matrix from xGEQRF 
'MRQ',  # multiply by the orthogonal/unitary matrix from xGERQF 
'MRZ',  # multiply by the orthogonal/unitary matrix from xTZRZF 
'MTR',  # multiply by the orthogonal/unitary matrix from xxxTRD 
'QLF',  # compute a QL factorization without pivoting 
'QPF',  # compute a QR factorization with column pivoting
'QRF',  # compute a QR factorization without pivoting 
'RFS',  # refine initial solution returned by TRS routines 
'RQF',  # compute an RQ factorization without pivoting 
'RZF',  # compute an RZ factorization without pivoting 
'SDC',  #compute using divide and conquer SVD 
'SEN',  # compute a basis and/or reciprocal condition number (sensitivity) of an 
        # invariant subspace 
'SJA',  # obtain singular values, and optionally vectors, using Jacobi's method 
'SNA',  # estimate reciprocal condition numbers of eigenvalue/-vector pairs
'SQR',  # compute singular values and/or singular vectors using the QR algorithm 
'SVP',  # preprocessing for GSVD 
'SYL',  # solve the Sylvester matrix equation 
'TRD',  # reduce a symmetric matrix to real symmetric tridiagonal form 
'TRF',  # compute a triangular factorization (LU, Cholesky, etc.) 
'TRI',  # compute inverse (based on triangular factorization)
'TRS',  # solve systems of linear equations (based on triangular factorization)
]

drivers = [\
'SV' ,  # factor the matrix and solve a system of equations
'SVX',  # equilibrate, factor, solve, compute error bounds and do iterative refinement, and estimate the condition number
'GLM',  #  solves the generalized linear regression model
'LS' ,  # solve over- or underdetermined linear system using orthogonal factorizations
'LSE',  #  solves the constrained linear least squares problem
'LSX',  # compute a minimum-norm solution using a complete orthogonal factorization (using QR with column pivoting xGEQPF, xTZRQF, and xLATZM)
'LSY',  # compute a minimum-norm solution using a complete orthogonal factorization (using QR with column pivoting xGEQP3, xTZRZF, and xORMRZ)
'LSS',  # solve least squares problem using the SVD
'LSD',  # solve least squares problem using the divide and conquer SVD
'EV' ,  # compute all eigenvalues and/or eigenvectors
'EVD',  # compute all eigenvalues and/or eigenvectors; if eigenvectors are desired, it uses a divide and conquer algorithm.
'EVX',  # compute selected eigenvalues and eigenvectors
'EVR',  # compute selected eigenvalues, and optionally, eigenvectors using the Relatively Robust Representation
'ES' ,  # compute all eigenvalues, Schur form, and/or Schur vectors
'ESX',  # compute all eigenvalues, Schur form, and/or Schur vectors and the conditioning # of selected eigenvalues or eigenvectors
'GV' ,  # compute generalized eigenvalues and/or generalized eigenvectors
'GVD',  # compute generalized eigenvalues, and optionally, generalized eigenvectors using a divide and conquer algorithm
'GVX',  # compute selected generalized eigenvalues, and optionally, generalized eigenvectors
'GS' ,  # compute generalized eigenvalues, Schur form, and/or Schur vectors
'SDD',  # compute the divide-and-conquer SVD
'SVD',  # compute the SVD and/or singular vectors
]

pattern = r"\b([{}](?:{})(?:{}))\b".format( \
        '' .join( data_types ),
        '|'.join( matrix_types ),
        '|'.join( computations + drivers ),
        )
