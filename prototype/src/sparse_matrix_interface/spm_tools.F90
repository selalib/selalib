MODULE SPM_TOOLS
USE SPM_DEF
IMPLICIT NONE

CONTAINS
! .............................................
SUBROUTINE CSC_ADD_TO_MATRIX(api_indices, api_indptr, apr_a, ar_value, ai_A, ai_Aprime, IERROR)
IMPLICIT NONE
INTEGER(KIND=SPM_INTL_KIND), DIMENSION(:), POINTER,      INTENT(IN)     :: api_indices
INTEGER(KIND=SPM_INTS_KIND), DIMENSION(:), POINTER,      INTENT(IN)     :: api_indptr
REAL(KIND=SPM_COEF_KIND)   , DIMENSION(:), POINTER,      INTENT(INOUT)  :: apr_a
REAL(KIND=SPM_COEF_KIND)                          ,      INTENT(IN)  :: ar_value
INTEGER(KIND=SPM_INTS_KIND)                       ,      INTENT(IN)  :: ai_A
INTEGER(KIND=SPM_INTS_KIND)                       ,      INTENT(IN)  :: ai_Aprime
INTEGER(KIND=SPM_INTS_KIND)                       ,      INTENT(OUT) :: IERROR
! LOCAL
INTEGER:: li_i
INTEGER:: li_j
INTEGER:: li_k
PRINT *, 'CSC_ADD_TO_MATRIX: not yet implemented'
DO li_k = api_indptr ( ai_Aprime ), api_indptr ( ai_Aprime + 1 ) - 1
li_i = api_indices ( li_k )

IF (li_i == ai_A) THEN
apr_a ( li_k ) = apr_a ( li_k ) + ar_value
EXIT
END IF

END DO
IERROR = SPM_SUCCESS

END SUBROUTINE CSC_ADD_TO_MATRIX
! .............................................
SUBROUTINE CSR_ADD_TO_MATRIX(api_indices, api_indptr, apr_a, ar_value, ai_A, ai_Aprime, IERROR)
IMPLICIT NONE
INTEGER(KIND=SPM_INTL_KIND), DIMENSION(:), POINTER,      INTENT(IN)     :: api_indices
INTEGER(KIND=SPM_INTS_KIND), DIMENSION(:), POINTER,      INTENT(IN)     :: api_indptr
REAL(KIND=SPM_COEF_KIND)   , DIMENSION(:), POINTER,      INTENT(INOUT)  :: apr_a
REAL(KIND=SPM_COEF_KIND)                          ,      INTENT(IN)  :: ar_value
INTEGER(KIND=SPM_INTS_KIND)                       ,      INTENT(IN)  :: ai_A
INTEGER(KIND=SPM_INTS_KIND)                       ,      INTENT(IN)  :: ai_Aprime
INTEGER(KIND=SPM_INTS_KIND)                       ,      INTENT(OUT) :: IERROR
! LOCAL
INTEGER:: li_i
INTEGER:: li_j
INTEGER:: li_k

! THE CURRENT LINE IS api_indptr(ai_A)
DO li_k = api_indptr ( ai_A ), api_indptr ( ai_A + 1 ) - 1
li_j = api_indices ( li_k )

IF (li_j == ai_Aprime) THEN
apr_a ( li_k ) = apr_a ( li_k ) + ar_value
EXIT
END IF

END DO
IERROR = SPM_SUCCESS

END SUBROUTINE CSR_ADD_TO_MATRIX
! .............................................
SUBROUTINE IJV_ADD_TO_MATRIX(api_i, api_j, apr_v, ar_value, ai_A, ai_Aprime, IERROR)
IMPLICIT NONE
INTEGER(KIND=SPM_INTL_KIND), DIMENSION(:), POINTER,      INTENT(INOUT)  :: api_i
INTEGER(KIND=SPM_INTL_KIND), DIMENSION(:), POINTER,      INTENT(INOUT)  :: api_j
REAL(KIND=SPM_COEF_KIND)   , DIMENSION(:), POINTER,      INTENT(INOUT)  :: apr_v
REAL(KIND=SPM_COEF_KIND)                          ,      INTENT(IN)     :: ar_value
INTEGER(KIND=SPM_INTS_KIND)                       ,      INTENT(IN)     :: ai_A
INTEGER(KIND=SPM_INTS_KIND)                       ,      INTENT(IN)     :: ai_Aprime
INTEGER(KIND=SPM_INTS_KIND)                       ,      INTENT(OUT)    :: IERROR
! LOCAL
INTEGER(KIND=SPM_INTL_KIND), DIMENSION(:), POINTER :: lpi_i
INTEGER(KIND=SPM_INTS_KIND), DIMENSION(:), POINTER :: lpi_j
REAL(KIND=SPM_COEF_KIND)   , DIMENSION(:), POINTER :: lpr_v
INTEGER:: li_n
INTEGER:: li_size

li_n = api_I (0)

! reallocate memory if necessary
li_size = SIZE(api_I, 1) - 1 ! minus 1 because we start from 0

IF (li_n==li_size) then

ALLOCATE(lpi_I(0:li_n)) ; ALLOCATE(lpi_J(0:li_n)) ; ALLOCATE(lpr_V(0:li_n))
lpi_I(0:li_n) = api_I(0:li_n) ; lpi_J(0:li_n) = api_J(0:li_n) ; lpr_V(0:li_n) = apr_V(0:li_n)
DEALLOCATE(api_I, api_J, apr_V)

ALLOCATE(api_I(0:2 * li_size)) ; ALLOCATE(api_J(0:2 * li_size)) ; ALLOCATE(apr_V(0:2 * li_size))
api_I(0:li_n) = lpi_I(0:li_n) ; api_J(0:li_n) = lpi_J(0:li_n) ; apr_V(0:li_n) = lpr_V(0:li_n)
DEALLOCATE(lpi_I, lpi_J, lpr_V)

end if

! insert the new ijv
li_n = li_n + 1
api_I(li_n) = ai_A
api_J(li_n) = ai_Aprime 
apr_V(li_n) = ar_value 
api_I (0) = li_n
api_J (0) = li_n
apr_V (0) = li_n

IERROR = SPM_SUCCESS
END SUBROUTINE IJV_ADD_TO_MATRIX
! .............................................
SUBROUTINE CSC_GEMV(ai_transa, ai_nC, apr_a, api_indices, api_indptr, apr_x, apr_y, IERROR)
! implments y := A*x or y := A'*x
IMPLICIT NONE
INTEGER(KIND=SPM_INTS_KIND),                             INTENT(IN)     :: ai_transa
INTEGER(KIND=SPM_INTS_KIND),                             INTENT(IN)     :: ai_nC
REAL(KIND=SPM_COEF_KIND)   , DIMENSION(:), POINTER,      INTENT(IN)     :: apr_a
INTEGER(KIND=SPM_INTL_KIND), DIMENSION(:), POINTER,      INTENT(IN)     :: api_indices
INTEGER(KIND=SPM_INTS_KIND), DIMENSION(:), POINTER,      INTENT(IN)     :: api_indptr
REAL(KIND=SPM_COEF_KIND)   , DIMENSION(:),      INTENT(IN)     :: apr_x
REAL(KIND=SPM_COEF_KIND)   , DIMENSION(:),      INTENT(INOUT)  :: apr_y
INTEGER(KIND=SPM_INTS_KIND)                       ,      INTENT(OUT) :: IERROR
!local var
INTEGER :: li_j, li_k_1, li_k_2

!IF (ai_transa==1) THEN
!PRINT *, 'CSC_GEMV: not yet implemented'
!STOP
!END IF
!!print *, SIZE(apr_a, 1)
!!print *, api_indptr(ai_nR)
!DO li_j = 1, ai_nC
!
!li_k_1 = api_indptr(li_i)
!li_k_2 = api_indptr(li_i + 1) - 1
!
!apr_y(li_i) = DOT_PRODUCT(apr_a(li_k_1: li_k_2), apr_x(api_indices(li_k_1: li_k_2)))
!    
!END DO
PRINT *, 'CSC_GEMV: not yet implemented'
IERROR = SPM_SUCCESS
END SUBROUTINE CSC_GEMV
! .............................................
SUBROUTINE CSR_GEMV(ai_transa, ai_nR, apr_a, api_indices, api_indptr, apr_x, apr_y, IERROR)
! implments y := A*x or y := A'*x
IMPLICIT NONE
INTEGER(KIND=SPM_INTS_KIND),                             INTENT(IN)     :: ai_transa
INTEGER(KIND=SPM_INTS_KIND),                             INTENT(IN)     :: ai_nR
REAL(KIND=SPM_COEF_KIND)   , DIMENSION(:), POINTER,      INTENT(IN)     :: apr_a
INTEGER(KIND=SPM_INTL_KIND), DIMENSION(:), POINTER,      INTENT(IN)     :: api_indices
INTEGER(KIND=SPM_INTS_KIND), DIMENSION(:), POINTER,      INTENT(IN)     :: api_indptr
REAL(KIND=SPM_COEF_KIND)   , DIMENSION(:),      INTENT(IN)     :: apr_x
REAL(KIND=SPM_COEF_KIND)   , DIMENSION(:),      INTENT(INOUT)  :: apr_y
INTEGER(KIND=SPM_INTS_KIND)                       ,      INTENT(OUT) :: IERROR
!local var
INTEGER :: li_i, li_k_1, li_k_2
INTEGER, DIMENSION(3) :: lpi_param

IF (ai_transa==1) THEN
PRINT *, 'CSR_GEMV: not yet implemented'
STOP
END IF
!print *, SIZE(apr_a, 1)
!print *, api_indptr(ai_nR)

#ifdef _CUSP
PRINT *, 'CSR_GEMV : MUST GIVE THE CORRECT VALUES FOR nC and nnz'
STOP
lpi_param(1) = ai_nR 
!lpi_param(2) = ai_nC
!lpi_param(3) = ai_nnz 
api_indptr = api_indptr - 1
api_indices= api_indices- 1
CALL cusp_csr_mv(api_indptr, api_indices, apr_a, lpi_param, apr_x, apr_y)
api_indptr = api_indptr + 1
api_indices= api_indices+ 1
#elif defined _MKL

CALL mkl_dcsrgemv('n', ai_nR &
, apr_a  &
, api_indptr &
, api_indices  &
, apr_x, apr_y)

#else    

DO li_i = 1, ai_nR

li_k_1 = api_indptr(li_i)
li_k_2 = api_indptr(li_i + 1) - 1

apr_y(li_i) = DOT_PRODUCT(apr_a(li_k_1: li_k_2), apr_x(api_indices(li_k_1: li_k_2)))
    
END DO
#endif

IERROR = SPM_SUCCESS
END SUBROUTINE CSR_GEMV
! .............................................
SUBROUTINE CSC_MV(ai_transa, ai_nR, ai_nC, ar_alpha, apr_a, api_indices, api_indptr, apr_x, ar_beta, apr_y, IERROR)
! implments y := alpha*A*x + beta*y
! or y := alpha*A'*x + beta*y,
IMPLICIT NONE
INTEGER(KIND=SPM_INTS_KIND),                             INTENT(IN)     :: ai_transa
INTEGER(KIND=SPM_INTS_KIND),                             INTENT(IN)     :: ai_nR
INTEGER(KIND=SPM_INTS_KIND),                             INTENT(IN)     :: ai_nC
REAL(KIND=SPM_COEF_KIND)   ,                             INTENT(IN)     :: ar_alpha
REAL(KIND=SPM_COEF_KIND)   , DIMENSION(:), POINTER,      INTENT(IN)     :: apr_a
INTEGER(KIND=SPM_INTL_KIND), DIMENSION(:), POINTER,      INTENT(IN)     :: api_indices
INTEGER(KIND=SPM_INTS_KIND), DIMENSION(:), POINTER,      INTENT(IN)     :: api_indptr
REAL(KIND=SPM_COEF_KIND)   , DIMENSION(:),      INTENT(IN)     :: apr_x
REAL(KIND=SPM_COEF_KIND)   ,                             INTENT(IN)     :: ar_beta
REAL(KIND=SPM_COEF_KIND)   , DIMENSION(:),      INTENT(INOUT)  :: apr_y
INTEGER(KIND=SPM_INTS_KIND)                       ,      INTENT(OUT) :: IERROR
!local var
INTEGER :: li_i, li_k_1, li_k_2

IF (ai_transa==1) THEN
PRINT *, 'CSC_MV: not yet implemented'
STOP
END IF
PRINT *, 'CSC_MV: not yet implemented'
!DO li_i = 1, ai_nR
!
!li_k_1 = api_indptr(li_i)
!li_k_2 = api_indptr(li_i + 1) - 1
!
!apr_y(li_i) = ar_alpha * DOT_PRODUCT(apr_a(li_k_1: li_k_2), apr_x(api_indices(li_k_1: li_k_2))) &
!+ ar_beta * apr_y(li_i)
!    
!END DO

IERROR = SPM_SUCCESS
END SUBROUTINE CSC_MV
! .............................................
SUBROUTINE CSR_MV(ai_transa, ai_nR, ai_nC, ar_alpha, apr_a, api_indices, api_indptr, apr_x, ar_beta, apr_y, IERROR)
! implments y := alpha*A*x + beta*y
! or y := alpha*A'*x + beta*y,
IMPLICIT NONE
INTEGER(KIND=SPM_INTS_KIND),                             INTENT(IN)     :: ai_transa
INTEGER(KIND=SPM_INTS_KIND),                             INTENT(IN)     :: ai_nR
INTEGER(KIND=SPM_INTS_KIND),                             INTENT(IN)     :: ai_nC
REAL(KIND=SPM_COEF_KIND)   ,                             INTENT(IN)     :: ar_alpha
REAL(KIND=SPM_COEF_KIND)   , DIMENSION(:), POINTER,      INTENT(IN)     :: apr_a
INTEGER(KIND=SPM_INTL_KIND), DIMENSION(:), POINTER,      INTENT(IN)     :: api_indices
INTEGER(KIND=SPM_INTS_KIND), DIMENSION(:), POINTER,      INTENT(IN)     :: api_indptr
REAL(KIND=SPM_COEF_KIND)   , DIMENSION(:),      INTENT(IN)     :: apr_x
REAL(KIND=SPM_COEF_KIND)   ,                             INTENT(IN)     :: ar_beta
REAL(KIND=SPM_COEF_KIND)   , DIMENSION(:),      INTENT(INOUT)  :: apr_y
INTEGER(KIND=SPM_INTS_KIND)                       ,      INTENT(OUT) :: IERROR
!local var
INTEGER :: li_i, li_k_1, li_k_2

IF (ai_transa==1) THEN
PRINT *, 'CSR_MV: not yet implemented'
STOP
END IF

DO li_i = 1, ai_nR

li_k_1 = api_indptr(li_i)
li_k_2 = api_indptr(li_i + 1) - 1

apr_y(li_i) = ar_alpha * DOT_PRODUCT(apr_a(li_k_1: li_k_2), apr_x(api_indices(li_k_1: li_k_2))) &
+ ar_beta * apr_y(li_i)
    
END DO

IERROR = SPM_SUCCESS
END SUBROUTINE CSR_MV

END MODULE SPM_TOOLS
