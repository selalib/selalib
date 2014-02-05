MODULE SPM
USE SPM_DEF
USE SPM_TOOLS
USE BASIC_CG_MODULE
USE GAUSS_SEIDEL_MODULE
IMPLICIT NONE 
    
    !**************************************************
    !		SparseMatrix
    !**************************************************
    TYPE, PRIVATE :: MATRIX
        INTEGER(KIND=SPM_INTS_KIND) :: oi_fmt !the current format of the matrix 
        INTEGER(KIND=SPM_INTS_KIND) :: oi_nR !NUMBER OF ROWS
        INTEGER(KIND=SPM_INTS_KIND) :: oi_nC !NUMBER OF COLUMNS
        INTEGER(KIND=SPM_INTL_KIND) :: oi_nnz !NUMBER OF NON ZERO ELTS
        INTEGER(KIND=SPM_INTL_KIND) :: oi_iassembly !Current index for the assembling process

        INTEGER(KIND=SPM_INTS_KIND), DIMENSION(:), POINTER :: opi_indptr   !CSC/CSR format
        INTEGER(KIND=SPM_INTL_KIND), DIMENSION(:), POINTER :: opi_indices !CSC/CSR format

        REAL(KIND=SPM_COEF_KIND)   , DIMENSION(:), POINTER :: opr_a !values, for CSC/CSR

        INTEGER(KIND=SPM_INTL_KIND), DIMENSION(:), POINTER :: opi_I !IJV format, unsorted
        INTEGER(KIND=SPM_INTL_KIND), DIMENSION(:), POINTER :: opi_J !IJV format, unsorted
        REAL(KIND=SPM_COEF_KIND)   , DIMENSION(:), POINTER :: opr_V !values IJV format, unsorted

        LOGICAL :: ol_assolver ! true if we use the current matrix to solve a linear system

        ! ... 
        ! Upper and Lower decompositions of the matrix
        ! ...
        INTEGER(KIND=SPM_INTL_KIND) :: oi_NNZ_U
        INTEGER(KIND=SPM_INTL_KIND) :: oi_NNZ_L

        REAL(KIND=SPM_COEF_KIND)   , DIMENSION(:), POINTER :: opr_U 
        INTEGER(KIND=SPM_INTS_KIND), DIMENSION(:), POINTER :: opi_IU
        INTEGER(KIND=SPM_INTL_KIND), DIMENSION(:), POINTER :: opi_JU

        REAL(KIND=SPM_COEF_KIND)   , DIMENSION(:), POINTER :: opr_L 
        INTEGER(KIND=SPM_INTS_KIND), DIMENSION(:), POINTER :: opi_IL
        INTEGER(KIND=SPM_INTL_KIND), DIMENSION(:), POINTER :: opi_JL

        REAL(KIND=SPM_COEF_KIND)   , DIMENSION(:), POINTER  :: opr_DIAGA         
        ! ...
    END TYPE MATRIX
    !***************************************************

    INTEGER, PRIVATE                                 :: mi_nmatrices
    TYPE(MATRIX), DIMENSION(:), POINTER, PRIVATE     :: mpo_M
    INTEGER(KIND=SPM_INTS_KIND), DIMENSION(:,:), POINTER, PRIVATE     :: mpi_iparam_solver
    REAL(KIND=SPM_COEF_KIND)   , DIMENSION(:,:), POINTER, PRIVATE     :: mpr_rparam_solver

CONTAINS

! .............................................
SUBROUTINE SPM_INITIALIZE(nmatrices, IERROR)
IMPLICIT NONE
INTEGER(KIND=SPM_INTS_KIND),      INTENT(IN)  :: nmatrices
INTEGER(KIND=SPM_INTS_KIND),      INTENT(OUT) :: IERROR
! LOCAL
INTEGER :: ID

mi_nmatrices = nmatrices
ALLOCATE ( mpo_M(0:nmatrices - 1) )
ALLOCATE ( mpi_iparam_solver(0:nmatrices - 1,SPM_NIPARAM_SOLVER ) )
ALLOCATE ( mpr_rparam_solver(0:nmatrices - 1,SPM_NRPARAM_SOLVER ) )

DO ID = 0, nmatrices-1
mpo_M (ID) % oi_fmt = SPM_NULL_FMT
mpo_M (ID) % ol_assolver = .FALSE.
END DO


IERROR = SPM_SUCCESS

END SUBROUTINE SPM_INITIALIZE
! .............................................

! .............................................
SUBROUTINE SPM_FINALIZE(IERROR)
IMPLICIT NONE
INTEGER(KIND=SPM_INTS_KIND),      INTENT(OUT) :: IERROR
! LOCAL
INTEGER(KIND=SPM_INTS_KIND) :: ID

DEALLOCATE ( mpo_M )
IERROR = SPM_SUCCESS 

END SUBROUTINE SPM_FINALIZE
! .............................................

! .............................................
SUBROUTINE SPM_SetDefaultOptions(ID, stratnum, IERROR)
IMPLICIT NONE
INTEGER(KIND=SPM_INTS_KIND),      INTENT(IN)  :: ID
INTEGER(KIND=SPM_INTS_KIND),      INTENT(IN)  :: stratnum
INTEGER(KIND=SPM_INTS_KIND),      INTENT(OUT) :: IERROR

PRINT *, 'SPM_SetDefaultOptions : Not Yet Implemented.'
IERROR = SPM_SUCCESS
END SUBROUTINE SPM_SetDefaultOptions
! .............................................

! .............................................
SUBROUTINE SPM_SetOptionINT(ID, NUMBER, VALUE, IERROR)
IMPLICIT NONE
INTEGER(KIND=SPM_INTS_KIND),      INTENT(IN)  :: ID
INTEGER(KIND=SPM_INTS_KIND),      INTENT(IN)  :: NUMBER 
INTEGER(KIND=SPM_INTS_KIND),      INTENT(IN)  :: VALUE
INTEGER(KIND=SPM_INTS_KIND),      INTENT(OUT) :: IERROR

mpi_iparam_solver (ID, NUMBER) = VALUE

IERROR = SPM_SUCCESS
END SUBROUTINE SPM_SetOptionINT
! .............................................

! .............................................
SUBROUTINE SPM_SetOptionREAL(ID, NUMBER, VALUE, IERROR)
IMPLICIT NONE
INTEGER(KIND=SPM_INTS_KIND),      INTENT(IN)  :: ID
INTEGER(KIND=SPM_INTS_KIND),      INTENT(IN)  :: NUMBER 
REAL(KIND=SPM_COEF_KIND)   ,      INTENT(IN)  :: VALUE
INTEGER(KIND=SPM_INTS_KIND),      INTENT(OUT) :: IERROR

mpr_rparam_solver (ID, NUMBER) = VALUE

IERROR = SPM_SUCCESS
END SUBROUTINE SPM_SetOptionREAL
! .............................................

! .............................................
SUBROUTINE SPM_GetOptionINT(ID, NUMBER, VALUE, IERROR)
IMPLICIT NONE
INTEGER(KIND=SPM_INTS_KIND),      INTENT(IN)  :: ID
INTEGER(KIND=SPM_INTS_KIND),      INTENT(IN)  :: NUMBER 
INTEGER(KIND=SPM_INTS_KIND),      INTENT(OUT)  :: VALUE
INTEGER(KIND=SPM_INTS_KIND),      INTENT(OUT) :: IERROR

VALUE = mpi_iparam_solver (ID, NUMBER)

IERROR = SPM_SUCCESS
END SUBROUTINE SPM_GetOptionINT
! .............................................

! .............................................
SUBROUTINE SPM_GetOptionREAL(ID, NUMBER, VALUE, IERROR)
IMPLICIT NONE
INTEGER(KIND=SPM_INTS_KIND),      INTENT(IN)  :: ID
INTEGER(KIND=SPM_INTS_KIND),      INTENT(IN)  :: NUMBER 
REAL(KIND=SPM_COEF_KIND)   ,      INTENT(OUT)  :: VALUE
INTEGER(KIND=SPM_INTS_KIND),      INTENT(OUT) :: IERROR

VALUE = mpr_rparam_solver (ID, NUMBER)

IERROR = SPM_SUCCESS
END SUBROUTINE SPM_GetOptionREAL
! .............................................

! .............................................
SUBROUTINE SPM_GRAPHBEGIN(ID, N, EDGENBR, IERROR)
IMPLICIT NONE
INTEGER(KIND=SPM_INTS_KIND),      INTENT(IN)  :: ID
INTEGER(KIND=SPM_INTS_KIND),      INTENT(IN)  :: N
INTEGER(KIND=SPM_INTL_KIND),      INTENT(IN)  :: EDGENBR
INTEGER(KIND=SPM_INTS_KIND),      INTENT(OUT) :: IERROR

PRINT *, 'SPM_GRAPHBEGIN : Not Yet Implemented.'
IERROR = SPM_SUCCESS
END SUBROUTINE SPM_GRAPHBEGIN
! .............................................

! .............................................
SUBROUTINE SPM_GRAPHEDGE(ID, ROW, COL, IERROR)
IMPLICIT NONE
INTEGER(KIND=SPM_INTS_KIND),      INTENT(IN)  :: ID
INTEGER(KIND=SPM_INTS_KIND),      INTENT(IN)  :: ROW 
INTEGER(KIND=SPM_INTS_KIND),      INTENT(IN)  :: COL
INTEGER(KIND=SPM_INTS_KIND),      INTENT(OUT) :: IERROR

PRINT *, 'SPM_GRAPHEDGE : Not Yet Implemented.'
IERROR = SPM_SUCCESS
END SUBROUTINE SPM_GRAPHEDGE
! .............................................

! .............................................
SUBROUTINE SPM_GRAPHEND(ID, IERROR)
IMPLICIT NONE
INTEGER(KIND=SPM_INTS_KIND),      INTENT(IN)  :: ID
INTEGER(KIND=SPM_INTS_KIND),      INTENT(OUT) :: IERROR

PRINT *, ' SPM_GRAPHEND: Not Yet Implemented.'
IERROR = SPM_SUCCESS
END SUBROUTINE SPM_GRAPHEND
! .............................................

! .............................................
SUBROUTINE SPM_GraphGlobalCSC(ID, N, INDICES, INDPTR, ROOT, IERROR)
IMPLICIT NONE
!id	Solver instance identIFication number.
!N	Global number of columns
!INDICES	Index of the first element of each column in INDPTR array.
!INDPTR	Global row number array.
!root	Root processor : this processor enter the global data (-1 for all processors).
INTEGER(KIND=SPM_INTS_KIND),                    INTENT(IN)  :: ID
INTEGER(KIND=SPM_INTS_KIND),                    INTENT(IN)  :: N
INTEGER(KIND=SPM_INTL_KIND), DIMENSION(:),      INTENT(IN)  :: INDICES
INTEGER(KIND=SPM_INTS_KIND), DIMENSION(:),      INTENT(IN)  :: INDPTR
INTEGER(KIND=SPM_INTS_KIND),                    INTENT(IN)  :: ROOT 
INTEGER(KIND=SPM_INTS_KIND),                    INTENT(OUT) :: IERROR
! LOCAL
TYPE(MATRIX), POINTER :: M
INTEGER :: NR, NC, NNZ
PRINT *, 'SPM_GraphGlobalCSC: not yet implemented'
!M => mpo_M ( ID )

NNZ     = SIZE (INDICES, 1)
NC      = SIZE (INDPTR, 1) - 1 
NR      = N

mpo_M ( ID ) % oi_fmt = SPM_CSC_FMT
mpo_M ( ID ) % oi_nnz = NNZ 
mpo_M ( ID ) % oi_nR  = NR
mpo_M ( ID ) % oi_nC  = NC

ALLOCATE ( mpo_M ( ID ) % opi_indices(1:NNZ) )
ALLOCATE ( mpo_M ( ID ) % opi_indptr  (1:NC+1) )
ALLOCATE ( mpo_M ( ID ) % opr_a     (1:NNZ) )

mpo_M ( ID ) % opi_indices (1:NNZ)     = INDICES (1:NNZ) 
mpo_M ( ID ) % opi_indptr   (1:NR+1)    = INDPTR   (1:NR+1) 
mpo_M ( ID ) % opr_a                  = 0.0

IERROR = SPM_SUCCESS
END SUBROUTINE SPM_GraphGlobalCSC
! .............................................

! .............................................
SUBROUTINE SPM_GraphGlobalCSR(ID, N, INDICES, INDPTR, ROOT, IERROR)
IMPLICIT NONE
INTEGER(KIND=SPM_INTS_KIND),                    INTENT(IN)  :: ID
INTEGER(KIND=SPM_INTS_KIND),                    INTENT(IN)  :: N
INTEGER(KIND=SPM_INTL_KIND), DIMENSION(:),      INTENT(IN)  :: INDICES
INTEGER(KIND=SPM_INTS_KIND), DIMENSION(:),      INTENT(IN)  :: INDPTR
INTEGER(KIND=SPM_INTS_KIND),                    INTENT(IN)  :: ROOT 
INTEGER(KIND=SPM_INTS_KIND),                    INTENT(OUT) :: IERROR
! LOCAL
TYPE(MATRIX), POINTER :: M
INTEGER :: NR, NC, NNZ

!M => mpo_M ( ID )

NNZ     = SIZE (INDICES, 1)
NR      = SIZE (INDPTR, 1) - 1 
NC      = N

mpo_M ( ID ) % oi_fmt = SPM_CSR_FMT
mpo_M ( ID ) % oi_nnz = NNZ 
mpo_M ( ID ) % oi_nR  = NR
mpo_M ( ID ) % oi_nC  = NC

ALLOCATE ( mpo_M ( ID ) % opi_indices(1:NNZ) )
ALLOCATE ( mpo_M ( ID ) % opi_indptr  (1:NR+1) )
ALLOCATE ( mpo_M ( ID ) % opr_a     (1:NNZ) )

mpo_M ( ID ) % opi_indices (1:NNZ)     = INDICES (1:NNZ) 
mpo_M ( ID ) % opi_indptr   (1:NR+1)    = INDPTR   (1:NR+1) 
mpo_M ( ID ) % opr_a                  = 0.0

IERROR = SPM_SUCCESS
END SUBROUTINE SPM_GraphGlobalCSR
! .............................................

! .............................................
SUBROUTINE SPM_GraphGlobalIJV(ID, NR, NC, ROOT, IERROR)
IMPLICIT NONE
INTEGER(KIND=SPM_INTS_KIND),                    INTENT(IN)  :: ID
INTEGER(KIND=SPM_INTS_KIND),                    INTENT(IN)  :: NR
INTEGER(KIND=SPM_INTS_KIND),                    INTENT(IN)  :: NC
INTEGER(KIND=SPM_INTS_KIND),                    INTENT(IN)  :: ROOT 
INTEGER(KIND=SPM_INTS_KIND),                    INTENT(OUT) :: IERROR
! LOCAL
TYPE(MATRIX), POINTER :: M

M => mpo_M ( ID )

M % oi_fmt = SPM_IJV_FMT
M % oi_nnz = -1 
M % oi_nR  = NR
M % oi_nC  = NC

! opi_I(0) = length of opi_I
ALLOCATE ( M % opi_I (0:SPM_IJV_MAXSIZE) )
ALLOCATE ( M % opi_J (0:SPM_IJV_MAXSIZE) )
ALLOCATE ( M % opr_V (0:SPM_IJV_MAXSIZE) )

M % opi_I       = 0 
M % opi_J       = 0 
M % opr_V       = 0.0

IERROR = SPM_SUCCESS
END SUBROUTINE SPM_GraphGlobalIJV
! .............................................

! .............................................
SUBROUTINE SPM_CLEAN(ID, IERROR)
IMPLICIT NONE
INTEGER(KIND=SPM_INTS_KIND),      INTENT(IN)  :: ID
INTEGER(KIND=SPM_INTS_KIND),      INTENT(OUT) :: IERROR
! LOCAL
TYPE(MATRIX), POINTER :: M

M => mpo_M ( ID )

! we free the solver if used
IF ( M % ol_assolver ) THEN
CALL SPM_FREE_SOLVE(ID, IERROR)
END IF

IF (M % oi_fmt==SPM_COO_FMT) THEN
        PRINT *, 'SPM_CLEAN - SPM_COO_FMT : Not done yet'
END IF

IF ( (M % oi_fmt==SPM_CSC_FMT) &
.OR. (M % oi_fmt==SPM_CSR_FMT) )THEN
DEALLOCATE ( M % opi_indices )
DEALLOCATE ( M % opi_indptr )
DEALLOCATE ( M % opr_a )
END IF

IF (M % oi_fmt==SPM_IJV_FMT) THEN
DEALLOCATE ( M % opi_I )
DEALLOCATE ( M % opi_J )
DEALLOCATE ( M % opr_V )
END IF

IERROR = SPM_SUCCESS
END SUBROUTINE SPM_CLEAN
! .............................................

! .............................................
SUBROUTINE SPM_CLEANALL(IERROR)
IMPLICIT NONE
INTEGER(KIND=SPM_INTS_KIND),      INTENT(OUT) :: IERROR
! LOCAL
INTEGER  :: ID


DO ID = 0, mi_nmatrices-1
CALL SPM_CLEAN(ID, IERROR)
END DO

IERROR = SPM_SUCCESS
END SUBROUTINE SPM_CLEANALL
! .............................................

! .............................................
SUBROUTINE SPM_MATRIXRESET(ID, IERROR)
!(ID, IERROR)
IMPLICIT NONE
INTEGER(KIND=SPM_INTS_KIND),      INTENT(IN)  :: ID
INTEGER(KIND=SPM_INTS_KIND),      INTENT(OUT) :: IERROR
! LOCAL
TYPE(MATRIX), POINTER :: M

M => mpo_M ( ID )

IF (M % oi_fmt==SPM_COO_FMT) THEN
        PRINT *, 'SPM_CLEAN - SPM_COO_FMT : Not done yet'
END IF

IF ( (M % oi_fmt==SPM_CSC_FMT) &
.OR. (M % oi_fmt==SPM_CSR_FMT) )THEN
M % opr_a = 0.0
END IF

IF (M % oi_fmt==SPM_IJV_FMT) THEN
M % opi_I = 0 
M % opi_J = 0
M % opr_V = 0.0
END IF

IERROR = SPM_SUCCESS
END SUBROUTINE SPM_MATRIXRESET
! .............................................

! .............................................
SUBROUTINE SPM_ASSEMBLYBEGIN(ID, COEFNBR, OP, OP2, MODE, SYM, IERROR)
! in the case of ijv unsorted. coefnbr is not used and can be set to 0
IMPLICIT NONE
INTEGER(KIND=SPM_INTS_KIND),      INTENT(IN)  :: ID 
INTEGER(KIND=SPM_INTS_KIND),      INTENT(IN)  :: OP
INTEGER(KIND=SPM_INTS_KIND),      INTENT(IN)  :: OP2
INTEGER(KIND=SPM_INTS_KIND),      INTENT(IN)  :: MODE
INTEGER(KIND=SPM_INTS_KIND),      INTENT(IN)  :: SYM
INTEGER(KIND=SPM_INTL_KIND),      INTENT(IN)  :: COEFNBR
INTEGER(KIND=SPM_INTS_KIND),      INTENT(OUT) :: IERROR
! LOCAL
TYPE(MATRIX), POINTER :: M

M => mpo_M ( ID )

IERROR = SPM_SUCCESS
END SUBROUTINE SPM_ASSEMBLYBEGIN
! .............................................

! .............................................
SUBROUTINE SPM_ASSEMBLYSETVALUE(ID, ROW, COL, VALUE, IERROR)
INTEGER(KIND=SPM_INTS_KIND),      INTENT(IN)  :: ID
INTEGER(KIND=SPM_INTS_KIND),      INTENT(IN)  :: ROW
INTEGER(KIND=SPM_INTS_KIND),      INTENT(IN)  :: COL
REAL(KIND=SPM_COEF_KIND)   ,      INTENT(IN)  :: VALUE
INTEGER(KIND=SPM_INTS_KIND),      INTENT(OUT) :: IERROR
! LOCAL
TYPE(MATRIX), POINTER :: M

M => mpo_M ( ID )

IF (M % oi_fmt==SPM_COO_FMT) THEN
PRINT *, 'SPM_ASSEMBLYSETVALUE : Not Yet Implemented.'
END IF

IF (M % oi_fmt==SPM_CSC_FMT) THEN
CALL CSC_ADD_TO_MATRIX(M % opi_indices, M % opi_indptr, M % opr_a, VALUE, ROW, COL, IERROR)
END IF

IF (M % oi_fmt==SPM_CSR_FMT) THEN
CALL CSR_ADD_TO_MATRIX(M % opi_indices, M % opi_indptr, M % opr_a, VALUE, ROW, COL, IERROR)
END IF

IF (M % oi_fmt==SPM_IJV_FMT) THEN
CALL IJV_ADD_TO_MATRIX(M % opi_I, M % opi_J, M % opr_V, VALUE, ROW, COL, IERROR)
END IF

IERROR = SPM_SUCCESS
END SUBROUTINE SPM_ASSEMBLYSETVALUE
! .............................................

! .............................................
SUBROUTINE SPM_ASSEMBLYEND(ID, IERROR)
!(ID, IERROR)
IMPLICIT NONE
INTEGER(KIND=SPM_INTS_KIND),      INTENT(IN)  :: ID
INTEGER(KIND=SPM_INTS_KIND),      INTENT(OUT) :: IERROR

! convert the matrix to a specIFic format IF needed
IERROR = SPM_SUCCESS
END SUBROUTINE SPM_ASSEMBLYEND
! .............................................

! .............................................
SUBROUTINE SPM_GetnR(ID, VALUE, IERROR)
IMPLICIT NONE
INTEGER(KIND=SPM_INTS_KIND),      INTENT(IN)  :: ID
INTEGER(KIND=SPM_INTS_KIND),      INTENT(OUT) :: VALUE
INTEGER(KIND=SPM_INTS_KIND),      INTENT(OUT) :: IERROR

VALUE = mpo_M (ID) % oi_nR

IERROR = SPM_SUCCESS
END SUBROUTINE SPM_GetnR
! .............................................

! .............................................
SUBROUTINE SPM_GetnC(ID, VALUE, IERROR)
IMPLICIT NONE
INTEGER(KIND=SPM_INTS_KIND),      INTENT(IN)  :: ID
INTEGER(KIND=SPM_INTS_KIND),      INTENT(OUT) :: VALUE
INTEGER(KIND=SPM_INTS_KIND),      INTENT(OUT) :: IERROR

VALUE = mpo_M (ID) % oi_nC

IERROR = SPM_SUCCESS
END SUBROUTINE SPM_GetnC
! .............................................

! .............................................
SUBROUTINE SPM_Getnnz(ID, VALUE, IERROR)
IMPLICIT NONE
INTEGER(KIND=SPM_INTS_KIND),      INTENT(IN)  :: ID
INTEGER(KIND=SPM_INTL_KIND),      INTENT(OUT) :: VALUE
INTEGER(KIND=SPM_INTS_KIND),      INTENT(OUT) :: IERROR

VALUE = mpo_M (ID) % oi_nnz

IERROR = SPM_SUCCESS
END SUBROUTINE SPM_Getnnz
! .............................................

! .............................................
SUBROUTINE SPM_GetINDPTR(ID, SIZE, ARRAY, IERROR)
IMPLICIT NONE
INTEGER(KIND=SPM_INTS_KIND)                     , INTENT(IN)  :: ID
INTEGER(KIND=SPM_INTS_KIND)                     , INTENT(IN)  :: SIZE
INTEGER(KIND=SPM_INTS_KIND)    , DIMENSION(SIZE), INTENT(OUT) :: ARRAY
INTEGER(KIND=SPM_INTS_KIND)                     , INTENT(OUT) :: IERROR
! LOCAL
TYPE(MATRIX), POINTER :: M

M => mpo_M ( ID )

IF ( (M % oi_fmt/=SPM_CSC_FMT) &
.AND. (M % oi_fmt/=SPM_CSR_FMT) )THEN
RETURN
END IF

ARRAY = M % opi_indptr

IERROR = SPM_SUCCESS
END SUBROUTINE SPM_GetINDPTR
! .............................................

! .............................................
SUBROUTINE SPM_GetINDICES(ID, SIZE, ARRAY, IERROR)
IMPLICIT NONE
INTEGER(KIND=SPM_INTS_KIND)                     , INTENT(IN)  :: ID
INTEGER(KIND=SPM_INTL_KIND)                     , INTENT(IN)  :: SIZE
INTEGER(KIND=SPM_INTL_KIND)    , DIMENSION(SIZE), INTENT(OUT) :: ARRAY
INTEGER(KIND=SPM_INTS_KIND)                     , INTENT(OUT) :: IERROR
! LOCAL
TYPE(MATRIX), POINTER :: M

M => mpo_M ( ID )

IF ( (M % oi_fmt/=SPM_CSC_FMT) &
.AND. (M % oi_fmt/=SPM_CSR_FMT) )THEN
RETURN
END IF

ARRAY = M % opi_indices

IERROR = SPM_SUCCESS
END SUBROUTINE SPM_GetINDICES
! .............................................

! .............................................
SUBROUTINE SPM_GetDATA(ID, SIZE, ARRAY, IERROR)
IMPLICIT NONE
INTEGER(KIND=SPM_INTS_KIND)                     , INTENT(IN)  :: ID
INTEGER(KIND=SPM_INTL_KIND)                     , INTENT(IN)  :: SIZE
REAL(KIND=SPM_COEF_KIND)       , DIMENSION(SIZE), INTENT(OUT) :: ARRAY
INTEGER(KIND=SPM_INTS_KIND)                     , INTENT(OUT) :: IERROR
! LOCAL
TYPE(MATRIX), POINTER :: M

M => mpo_M ( ID )

IF ( (M % oi_fmt/=SPM_CSC_FMT) &
.AND. (M % oi_fmt/=SPM_CSR_FMT) )THEN
RETURN
END IF

ARRAY = M % opr_a

IERROR = SPM_SUCCESS
END SUBROUTINE SPM_GetDATA
! .............................................

! .............................................
SUBROUTINE SPM_SetDATA(ID, SIZE, ARRAY, IERROR)
IMPLICIT NONE
INTEGER(KIND=SPM_INTS_KIND)                     , INTENT(IN)  :: ID
INTEGER(KIND=SPM_INTL_KIND)                     , INTENT(IN)  :: SIZE
REAL(KIND=SPM_COEF_KIND)       , DIMENSION(SIZE), INTENT(IN)  :: ARRAY
INTEGER(KIND=SPM_INTS_KIND)                     , INTENT(OUT) :: IERROR
! LOCAL
TYPE(MATRIX), POINTER :: M

M => mpo_M ( ID )

IF ( (M % oi_fmt/=SPM_CSC_FMT) &
.AND. (M % oi_fmt/=SPM_CSR_FMT) )THEN
RETURN
END IF

M % opr_a = ARRAY

IERROR = SPM_SUCCESS
END SUBROUTINE SPM_SetDATA
! .............................................

! .............................................
SUBROUTINE SPM_GetDATAIJV(ID, SIZE, ARRAY_I, ARRAY_J, ARRAY_V, IERROR)
IMPLICIT NONE
INTEGER(KIND=SPM_INTS_KIND)                     , INTENT(IN)    :: ID
INTEGER(KIND=SPM_INTL_KIND)                 , INTENT(IN) :: SIZE
INTEGER(KIND=SPM_INTL_KIND), DIMENSION(SIZE), INTENT(INOUT) :: ARRAY_I
INTEGER(KIND=SPM_INTL_KIND), DIMENSION(SIZE), INTENT(INOUT) :: ARRAY_J
REAL(KIND=SPM_COEF_KIND)   , DIMENSION(SIZE), INTENT(INOUT) :: ARRAY_V
INTEGER(KIND=SPM_INTS_KIND)                     , INTENT(OUT)   :: IERROR
! LOCAL
TYPE(MATRIX), POINTER :: M

M => mpo_M ( ID )

IF (M % oi_fmt/=SPM_IJV_FMT) THEN
RETURN
END IF

IF (SIZE /= M % opi_I(0)) THEN
PRINT *, "SPM_DATAIJV error: please give the correct value for SIZE"
END IF

ARRAY_I (1:SIZE) = M % opi_I (1:SIZE)
ARRAY_J (1:SIZE) = M % opi_J (1:SIZE)
ARRAY_V (1:SIZE) = M % opr_V (1:SIZE)

IERROR = SPM_SUCCESS
END SUBROUTINE SPM_GetDATAIJV
! .............................................

! .............................................
SUBROUTINE SPM_GetSIZEIJV(ID, SIZE, IERROR)
IMPLICIT NONE
INTEGER(KIND=SPM_INTS_KIND)                     , INTENT(IN)    :: ID
INTEGER(KIND=SPM_INTL_KIND)                     , INTENT(INOUT) :: SIZE
INTEGER(KIND=SPM_INTS_KIND)                     , INTENT(OUT)   :: IERROR
! LOCAL
TYPE(MATRIX), POINTER :: M

M => mpo_M ( ID )

IF (M % oi_fmt/=SPM_IJV_FMT) THEN
RETURN
END IF

SIZE = M % opi_I(0)

IERROR = SPM_SUCCESS
END SUBROUTINE SPM_GetSIZEIJV
! .............................................

! .............................................
SUBROUTINE SPM_GEMV(ID, TRANSA, X, Y, IERROR)
INTEGER(KIND=SPM_INTS_KIND),                             INTENT(IN)     :: ID
INTEGER(KIND=SPM_INTS_KIND),                             INTENT(IN)     :: TRANSA 
REAL(KIND=SPM_COEF_KIND)   , DIMENSION(:),      INTENT(IN)     :: X
REAL(KIND=SPM_COEF_KIND)   , DIMENSION(:),      INTENT(INOUT)  :: Y
INTEGER(KIND=SPM_INTS_KIND),                             INTENT(OUT)    :: IERROR
! LOCAL
TYPE(MATRIX), POINTER :: M

M => mpo_M ( ID )

IF (M % oi_fmt==SPM_COO_FMT) THEN
PRINT *, 'SPM_GEMV: SPM_COO_FMT Not Yet Implemented.'
END IF

IF (M % oi_fmt==SPM_CSC_FMT) THEN
PRINT *, 'SPM_GEMV: SPM_CSC_FMT Not Yet Implemented.'
END IF

IF (M % oi_fmt==SPM_CSR_FMT) THEN
CALL CSR_GEMV(TRANSA, M % oi_nR, M % opr_a, M % opi_indices, M % opi_indptr, X, Y, IERROR)
END IF

IERROR = SPM_SUCCESS
END SUBROUTINE SPM_GEMV
! .............................................

! .............................................
SUBROUTINE SPM_MV(ID, TRANSA, ALPHA, X, BETA, Y, IERROR)
INTEGER(KIND=SPM_INTS_KIND),                             INTENT(IN)     :: ID
INTEGER(KIND=SPM_INTS_KIND),                             INTENT(IN)     :: TRANSA 
REAL(KIND=SPM_COEF_KIND)   ,                             INTENT(IN)     :: ALPHA
REAL(KIND=SPM_COEF_KIND)   , DIMENSION(:),      INTENT(IN)     :: X
REAL(KIND=SPM_COEF_KIND)   ,                             INTENT(IN)     :: BETA
REAL(KIND=SPM_COEF_KIND)   , DIMENSION(:),      INTENT(INOUT)  :: Y
INTEGER(KIND=SPM_INTS_KIND),                             INTENT(OUT)    :: IERROR
! LOCAL
TYPE(MATRIX), POINTER :: M

M => mpo_M ( ID )

IF (M % oi_fmt==SPM_COO_FMT) THEN
PRINT *, 'SPM_MV: Not Yet Implemented.'
END IF

IF (M % oi_fmt==SPM_CSC_FMT) THEN
CALL CSC_MV(TRANSA, M % oi_nR, M % oi_nC, ALPHA, M % opr_a, M % opi_indices, M % opi_indptr, X, BETA, Y, IERROR)
END IF

IF (M % oi_fmt==SPM_CSR_FMT) THEN
PRINT *, 'SPM_MV: Not Yet Implemented.'
END IF

IERROR = SPM_SUCCESS
END SUBROUTINE SPM_MV
! .............................................


! .............................................
SUBROUTINE SPM_INITIALIZE_SOLVE(ID, IERROR)
INTEGER(KIND=SPM_INTS_KIND),                             INTENT(IN)     :: ID
INTEGER(KIND=SPM_INTS_KIND),                             INTENT(OUT)    :: IERROR
! LOCAL
TYPE(MATRIX), POINTER :: M

M => mpo_M ( ID )

M % ol_assolver = .TRUE.

IF (M % oi_fmt==SPM_CSR_FMT) THEN
! *******************************************
!               BASIC-GS GAUSS SEIDEL 
! *******************************************
IF ( mpi_iparam_solver (ID, SPM_IPARAM_SOLVER) == SPM_SOLVER_BASIC_GS ) THEN
CALL GS_INITIALIZE( M % oi_nR, M % oi_nnz, M % opr_a, M % opi_indices, M % opi_indptr &
                  , M % oi_NNZ_U, M % opr_U, M % opi_JU, M % opi_IU &
                  , M % oi_NNZ_L, M % opr_L, M % opi_JL, M % opi_IL &
                  , M % opr_DIAGA)
END IF
! *******************************************

END IF

IF (M % oi_fmt==SPM_COO_FMT) THEN
PRINT *, 'SPM_SOLVE: COO Not Yet Implemented.'
END IF

IF (M % oi_fmt==SPM_CSC_FMT) THEN
PRINT *, 'SPM_SOLVE: CSC Not Yet Implemented.'
END IF


IERROR = SPM_SUCCESS
END SUBROUTINE SPM_INITIALIZE_SOLVE
! .............................................


! .............................................
SUBROUTINE SPM_FREE_SOLVE(ID, IERROR)
INTEGER(KIND=SPM_INTS_KIND),                             INTENT(IN)     :: ID
INTEGER(KIND=SPM_INTS_KIND),                             INTENT(OUT)    :: IERROR
! LOCAL
TYPE(MATRIX), POINTER :: M

M => mpo_M ( ID )

IF (M % oi_fmt==SPM_CSR_FMT) THEN
! *******************************************
!               BASIC-GS GAUSS SEIDEL 
! *******************************************
IF ( mpi_iparam_solver (ID, SPM_IPARAM_SOLVER) == SPM_SOLVER_BASIC_GS ) THEN
CALL GS_FREE( M % opr_U, M % opi_JU, M % opi_IU &
            , M % opr_L, M % opi_JL, M % opi_IL &
            , M % opr_DIAGA)
END IF
! *******************************************

END IF

IF (M % oi_fmt==SPM_COO_FMT) THEN
PRINT *, 'SPM_SOLVE: COO Not Yet Implemented.'
END IF

IF (M % oi_fmt==SPM_CSC_FMT) THEN
PRINT *, 'SPM_SOLVE: CSC Not Yet Implemented.'
END IF


IERROR = SPM_SUCCESS
END SUBROUTINE SPM_FREE_SOLVE
! .............................................


! .............................................
SUBROUTINE SPM_SOLVE(ID, X, Y, IERROR)
INTEGER(KIND=SPM_INTS_KIND),                             INTENT(IN)     :: ID
REAL(KIND=SPM_COEF_KIND)   , DIMENSION(:), POINTER,      INTENT(IN)     :: X
REAL(KIND=SPM_COEF_KIND)   , DIMENSION(:), POINTER,      INTENT(INOUT)  :: Y
INTEGER(KIND=SPM_INTS_KIND),                             INTENT(OUT)    :: IERROR
! LOCAL
TYPE(MATRIX), POINTER :: M
INTEGER(KIND=SPM_INTS_KIND), DIMENSION(10) :: lpi_param
REAL(KIND=SPM_COEF_KIND), DIMENSION(1) :: lpr_param
REAL(KIND=SPM_COEF_KIND) :: lr_tol
REAL(KIND=SPM_COEF_KIND) :: lr_resnorm
INTEGER(KIND=SPM_INTS_KIND) :: li_maxiter
INTEGER(KIND=SPM_INTS_KIND) :: li_niter
REAL(KIND=SPM_COEF_KIND)   , DIMENSION(:), POINTER  :: lpr_err

M => mpo_M ( ID )

li_maxiter      = mpi_iparam_solver (ID, SPM_IPARAM_SOLVER_MAXITER)
lr_tol          = mpr_rparam_solver (ID, SPM_RPARAM_SOLVER_RTOL)

IF (M % oi_fmt==SPM_CSR_FMT) THEN

select case ( mpi_iparam_solver (ID, SPM_IPARAM_SOLVER) )
! *******************************************
!               CUSP-CG
! *******************************************
case ( SPM_SOLVER_CUSP_CG )
#ifdef _CUSP
lpi_param(1) = M % oi_nR 
lpi_param(2) = M % oi_nC
lpi_param(3) = M % oi_nnz 
lpi_param(4) = li_maxiter  ! niterations 
lpi_param(5) = mpi_iparam_solver (ID, SPM_IPARAM_SOLVER_PRECOND) !preconditioner

lpr_param(1) = lr_tol  ! relative tolerance

!CALL cusp_cg_precond(M % opi_indptr, M % opi_indices, M % opr_a, lpi_param, lpr_param, X, Y)
CALL cusp_cg(M % opi_indptr, M % opi_indices, M % opr_a, lpi_param, lpr_param, X, Y)
#endif
! *******************************************

! *******************************************
!               CUSP-GMRES
! *******************************************
case ( SPM_SOLVER_CUSP_GMRES )
#ifdef _CUSP
lpi_param(1) = M % oi_nR 
lpi_param(2) = M % oi_nC
lpi_param(3) = M % oi_nnz 
lpi_param(4) = mpi_iparam_solver (ID, SPM_IPARAM_SOLVER_MAXITER)  ! niterations 
lpi_param(5) = mpi_iparam_solver (ID, SPM_IPARAM_SOLVER_RESTART)  ! n restart

lpr_param(1) = mpr_rparam_solver (ID, SPM_RPARAM_SOLVER_RTOL)  ! relative tolerance

CALL cusp_gmres(M % opi_indptr, M % opi_indices, M % opr_a, lpi_param, lpr_param, X, Y)
#endif
! *******************************************

! *******************************************
!               CUSP-BICGSTAB
! *******************************************
case ( SPM_SOLVER_CUSP_BICGSTAB )
#ifdef _CUSP
lpi_param(1) = M % oi_nR 
lpi_param(2) = M % oi_nC
lpi_param(3) = M % oi_nnz 
lpi_param(4) = mpi_iparam_solver (ID, SPM_IPARAM_SOLVER_MAXITER)  ! niterations 

lpr_param(1) = mpr_rparam_solver (ID, SPM_RPARAM_SOLVER_RTOL)  ! relative tolerance

CALL cusp_bicgstab(M % opi_indptr, M % opi_indices, M % opr_a, lpi_param, lpr_param, X, Y)
#endif
! *******************************************

! *******************************************
!               CUSP-JACOBI
! *******************************************
!case ( SPM_SOLVER_CUSP_JACOBI )
!#ifdef _CUSP
!lpi_param(1) = M % oi_nR 
!lpi_param(2) = M % oi_nC
!lpi_param(3) = M % oi_nnz 
!
!lpr_param(1) = mpr_rparam_solver (ID, SPM_RPARAM_SOLVER_OMEGA)  ! relative tolerance
!
!CALL cusp_jacobi(M % opi_indptr, M % opi_indices, M % opr_a, lpi_param, lpr_param, X, Y)
!#endif
! *******************************************

! *******************************************
!               BASIC-CG CONJUGATE GRADIENT
! *******************************************
case ( SPM_SOLVER_BASIC_CG )
ALLOCATE(lpr_err(li_maxiter))
CALL BASIC_CG ( M % oi_nR, M % opr_a, M % opi_indices, M % opi_indptr &
, li_maxiter, lr_tol, X, Y, li_niter, lr_resnorm, lpr_err )
mpi_iparam_solver (ID, SPM_IPARAM_SOLVER_NITER) = li_niter
mpr_rparam_solver (ID, SPM_RPARAM_SOLVER_RESNORM) = lr_resnorm
DEALLOCATE(lpr_err)
! *******************************************

! *******************************************
!               BASIC-GS GAUSS SEIDEL 
! *******************************************
case ( SPM_SOLVER_BASIC_GS )
ALLOCATE(lpr_err(li_maxiter))
CALL GS_SOLVE ( M % oi_nR, M % opr_a, M % opi_indices, M % opi_indptr &
              , M % opr_U, M % opi_JU, M % opi_IU &
              , M % opr_L, M % opi_JL, M % opi_IL &
              , M % opr_DIAGA                     &
              , X, Y, lr_tol, li_maxiter, li_niter, lr_resnorm, lpr_err)
mpi_iparam_solver (ID, SPM_IPARAM_SOLVER_NITER)   = li_niter
mpr_rparam_solver (ID, SPM_RPARAM_SOLVER_RESNORM) = lr_resnorm
DEALLOCATE(lpr_err)
! *******************************************
case Default
PRINT *, "SPM_SOLVE-CSR: Type Solver Not Yet implemented"
end select

END IF


IF (M % oi_fmt==SPM_COO_FMT) THEN
PRINT *, 'SPM_SOLVE-COO: Not Yet Implemented.'
END IF

IF (M % oi_fmt==SPM_CSC_FMT) THEN
PRINT *, 'SPM_SOLVE-CSC: Not Yet Implemented.'
END IF


IERROR = SPM_SUCCESS
END SUBROUTINE SPM_SOLVE
! .............................................

END MODULE SPM
