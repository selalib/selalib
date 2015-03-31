MODULE SOLVER
  USE SPM_DEF
  USE SPM
  USE MODEL_DEF
  USE JOREK_PARAM_DEF
  USE JOREK_PARAM
  IMPLICIT NONE
  
  INTEGER, PARAMETER, PRIVATE :: mi_dtllevel_base = 0


  INTEGER, PRIVATE :: mi_maxiter
  REAL(KIND=RK) , PRIVATE :: mr_rtol

  INTEGER, PRIVATE :: mi_solver_id
  INTEGER, PRIVATE :: mi_solver_ksp, mi_solver_pc_ksp
  INTEGER, PRIVATE :: mi_null_space

CONTAINS
! ----------------------------------------------------------------------
! Initialize the SPM solver module
! user must specify the solver type Petsc (LU, ...), ...
SUBROUTINE INITIALIZE_SOLVER( )
  IMPLICIT NONE
     INTEGER          :: ierr

     ! ... get solver parameters from JOREK_PARAM
     CALL JOREK_Param_GETInt(INT_LINEAR_SOLVER_ID_ID    , mi_solver_id , ierr) 
     CALL JOREK_Param_GETInt(INT_LINEAR_SOLVER_KRYLOV_ID, mi_solver_ksp, ierr) 
     CALL JOREK_Param_GETInt(INT_LINEAR_SOLVER_KRYLOV_PC_ID, mi_solver_pc_ksp, ierr)
     CALL JOREK_Param_GETInt(INT_SOLVER_MAXITER_ID      , mi_maxiter   , ierr) 
     CALL JOREK_Param_GETREAL(REAL_SOLVER_TOLERANCE_ID  , mr_rtol      , ierr)
     ! ...

     ! ... set solver parameters into SPM
     CALL SPM_SetOptionINT(Matrix_ID, SPM_IPARAM_SOLVER_MAXITER, mi_maxiter, ierr)
     CALL SPM_SetOptionREAL(Matrix_ID, SPM_RPARAM_SOLVER_RTOL, mr_rtol, ierr)

     CALL SPM_SetOptionINT(Matrix_ID, SPM_IPARAM_SOLVER, mi_solver_id, ierr)
     IF (mi_solver_id == SPM_SOLVER_PETSC_KSP) THEN
        CALL SPM_SetOptionINT(Matrix_ID, SPM_IPARAM_SOLVER_KSP_TYPE, mi_solver_ksp, ierr)
        CALL SPM_SetOptionINT(Matrix_ID, SPM_IPARAM_SOLVER_KSP_PC_TYPE, mi_solver_pc_ksp, ierr)
     END IF

     IF (mi_null_space == 1) THEN
        CALL SPM_SetOptionINT(Matrix_ID, SPM_IPARAM_SOLVER_NULL_SPACE, SPM_SOLVER_WITH_NULL_SPACE, ierr)
     END IF
     ! ...

     ! ... SPM Solver Initialization
     CALL SPM_INITIALIZE_SOLVE(Matrix_ID, ierr)
     ! ...

END SUBROUTINE INITIALIZE_SOLVER
! ----------------------------------------------------------------------
! deallocate memory used by the SPM solver module
SUBROUTINE FREE_SOLVER( )
  IMPLICIT NONE
     INTEGER          :: ierr

     CALL SPM_FREE_SOLVE(Matrix_ID, ierr)

END SUBROUTINE FREE_SOLVER
! ----------------------------------------------------------------------
! write here the kernel solver
! the kernel solver is the linear system to solve at each time step 
! (if evolution problem)
! solves My = x
SUBROUTINE KERNEL_SOLVER(X, Y, ierr)
  IMPLICIT NONE
     REAL(KIND=SPM_COEF_KIND)   , DIMENSION(:), INTENT(IN)  :: X 
     REAL(KIND=SPM_COEF_KIND)   , DIMENSION(:), INTENT(OUT) :: Y
     ! LOCAL
     INTEGER          :: ierr    

     ! solve X = M^-1 Y
!     CALL SPM_SOLVE(Matrix_ID, X, Y, ierr)

     CALL SPM_SETLOCALRHS(Matrix_ID, X, SPM_ASSEMBLY_OVW, SPM_ASSEMBLY_OVW, ierr)
     CALL SPM_SetGlobalRHS(Matrix_ID, X, SPM_ASSEMBLY_OVW, 0, ierr)
     CALL SPM_GetGlobalSolution(Matrix_ID, Y, 0, ierr)
!     CALL SPM_GetLocalSolution(Matrix_ID, Y, ierr)

END SUBROUTINE KERNEL_SOLVER
! ----------------------------------------------------------------------
! write here the full solve procedure by calling the kernel solver
SUBROUTINE SOLVE(X, Y)
  IMPLICIT NONE
     REAL(KIND=SPM_COEF_KIND)   , DIMENSION(:), INTENT(IN)  :: X 
     REAL(KIND=SPM_COEF_KIND)   , DIMENSION(:), INTENT(OUT) :: Y
     ! LOCAL
     INTEGER          :: ierr    

     CALL KERNEL_SOLVER(X, Y, ierr)

END SUBROUTINE SOLVE

END MODULE SOLVER
