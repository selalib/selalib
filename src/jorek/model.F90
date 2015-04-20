MODULE MODEL
  USE MODEL_DEF
  USE TypeDef        
  USE SPM_DEF
  USE SPM
  USE Assembly
  USE Evaluator 
  USE OUTPUT 
  USE DIRICHLET_MOD
  USE ELEMENT_CONTRIBUTION
  USE TESTCASE
  USE JOREK_GLOB
  USE GRAPH
  USE NUMBERING
  USE DIAGNOSTICS_COMPUTATION
  USE MESH_DEF
  USE MESH
  USE QUADRATURES_DEF
  USE QUADRATURES
  USE FIELD_DEF
  USE FIELD
  USE MATRIX_DEF
  USE MATRIX
  USE LINEAR_SOLVER_DEF
  USE LINEAR_SOLVER
  USE SPACE_DEF
  USE SPACE
  USE FEM_DEF
  USE FEM
  IMPLICIT NONE

CONTAINS
  ! ..................................................
  SUBROUTINE DEFINE_MODEL( )
  IMPLICIT NONE
    INTEGER, PARAMETER          :: N_DIM = 2
    INTEGER :: li_poloidal_basis
    INTEGER :: li_err
    CHARACTER(LEN = 1024)           :: dirname
    CHARACTER(LEN = 1024)           :: argname

    current_model       = 1

    n_var_unknown       = 1
    n_var_sys           = 1
    
    i_Vp_Rho            = 1

    nmatrices           = 1
    i_Vu_Rho            = 1 

    CALL JOREK_Param_GETInt(INT_MODES_M1_ID, mi_mode_m1, li_err)
    CALL JOREK_Param_GETInt(INT_MODES_N1_ID, mi_mode_n1, li_err)

    CALL JOREK_Param_GETReal(REAL_RGEO_ID,mr_R0,li_err)
    CALL JOREK_Param_GETReal(REAL_ZGEO_ID,mr_Z0,li_err)
    CALL JOREK_Param_GETReal(REAL_AMIN_ID,mr_a,li_err) 
    CALL JOREK_Param_GETReal(REAL_ACENTER_ID,mr_acenter,li_err)
    CALL JOREK_Param_GETInt(INT_NSTEP_MAX_ID,mi_nstep_max,li_err)

    ! ...
    CALL SPACE_CREATE(space_trial, fem_mesh)
    ptr_space_trial => space_trial
    CALL SPACE_CREATE(space_test, fem_mesh)
    ptr_space_test => space_test
    ! ...

    ! ...
    argname = "--geometry"
    CALL JOREK_GET_ARGUMENTS(argname, dirname, li_err)
    CALL CREATE_MESH(fem_mesh, dirname)
    ! ...

    CALL MODEL_CREATE(fem_model, fem_mesh, ptr_space_trial, ptr_space_test, n_var_unknown, n_var_sys)

    CALL FIELD_CREATE(field_U, fem_model % space_trial, field_name="density")
    ptr_field => field_U
    CALL MODEL_APPEND_FIELD(fem_model, ptr_field)

    CALL MATRIX_CREATE(matrix_Stiffnes)
    ptr_system => matrix_Stiffnes

    CALL MATRIX_APPEND_UNKNOWN_FIELD(ptr_system, ptr_field)
    CALL MODEL_APPEND_MATRIX(fem_model, ptr_system)

    CALL MODEL_INITIALIZE(fem_model)    

  END SUBROUTINE DEFINE_MODEL
  ! ..................................................

  ! ..................................................
  SUBROUTINE FREE_MODEL( )
  IMPLICIT NONE

    CALL MODEL_FREE(fem_model)

  END SUBROUTINE FREE_MODEL
  ! ..................................................

  ! ..................................................
  SUBROUTINE RUN_MODEL( )
  IMPLICIT NONE
     REAL(KIND=RK)    :: T_fin, T_deb
     CHARACTER(LEN = 1024)           :: filename
     CHARACTER(LEN = 1024)           :: argname
     INTEGER(KIND=SPM_INTS_KIND) :: nRows ! Number of Rows
     INTEGER                     :: li_err
     INTEGER                     :: li_myRank
     character(len=20) :: ls_stamp_default
     character(len=20)   :: ls_msg
     REAL(KIND=RK), DIMENSION(:), ALLOCATABLE  :: lpr_Y
     REAL(KIND=RK), DIMENSION(1:1,1:2) :: lpr_x
     REAL(KIND=RK), DIMENSION(1:1,1:2) :: lpr_v
     INTEGER                     :: li_i   

!     Find out the rank of the current proc
     li_myRank = 0
#ifdef MPI_ENABLED
      CALL MPI_Comm_rank ( MPI_COMM_WORLD, li_myRank, li_err)
#endif

     ! ... define weak formulation 
     ptr_system % ptr_matrix_contribution => Matrix_for_Vi_Vj
     ptr_system % ptr_rhs_contribution    => RHS_for_Vi
     ! ...

     ! ... Loop over elements 
     CALL CPU_TIME(T_deb)
     CALL MODEL_ASSEMBLY(fem_model, ptr_system)
     CALL CPU_TIME(T_fin)
     ! ...

!     print *, "RHS ", ptr_system % opr_global_rhs(1:6)

     ! ... example of solver calls
     argname = "--solver"
     CALL JOREK_GET_ARGUMENTS(argname, filename, li_err)
     CALL LINEAR_SOLVER_CREATE(solver, ptr_system, filename)
     CALL LINEAR_SOLVER_SOLVE(solver, ptr_system % opr_global_rhs, ptr_system % opr_global_unknown) 
     CALL LINEAR_SOLVER_FREE(solver)
     ! ...

     ! ... update related fields to the linear system
     CALL UPDATE_FIELDS(fem_model, ptr_system)
!     print *, ptr_field % opr_coeffs
     ! ...

     ! ... evaluate field
     lpr_x = 1.0
     lpr_v = 0.0
     CALL FIELD_EVALUATE(field_U, 1, lpr_x, lpr_v)
     print *, "field value at ", lpr_x, " is ", lpr_v
     ! ...

     CALL GET_NR_MATRIX(ptr_system, nRows)

     write(ls_msg,*) li_myRank
     ls_stamp_default = "-proc_" // TRIM ( ADJUSTL ( ls_msg ) )
     ls_stamp_default = TRIM ( ADJUSTL ( ADJUSTR ( ls_stamp_default ) ) )

     CALL getarg(0, filename)

     ALLOCATE(lpr_Y(nRows))
     lpr_Y=0.d0 
     CALL SPM_MATMULT(Matrix_A_ID, ptr_system % opr_global_unknown, lpr_Y, li_err)

     PRINT *, MAXVAL(lpr_Y-ptr_system % opr_global_rhs), MINVAL(lpr_Y-ptr_system % opr_global_rhs) 
!     print *, "UNKNOWN ", ptr_system % opr_global_unknown

     ! ... Evaluate unknowns on vertecies and compte model norms 
     CALL MODEL_DIAGNOSTICS(fem_model, ANALYTICAL_MODEL, Assembly_Diagnostics, Plot_Diagnostics, ptr_system)

      filename = TRIM(filename) 
     CALL MODEL_SAVE(fem_model, ANALYTICAL_MODEL, filename, 0)

     ! ...
     OPEN(UNIT=12, FILE=TRIM ("RHS" // ADJUSTL ( ADJUSTR ( ls_stamp_default ) ) )  // ".txt"&
             & , ACTION="write", STATUS="replace")
     DO li_i=1,nRows
        WRITE(12,*) ptr_system % opr_global_rhs(li_i) 
     END DO
     CLOSE(12)
     OPEN(UNIT=13, FILE=TRIM ("UNKNOWN" // ADJUSTL ( ADJUSTR ( ls_stamp_default ) ) )  // ".txt"&
             & , ACTION="write", STATUS="replace")
     DO li_i=1,nRows
        WRITE(13,*) ptr_system % opr_global_unknown(li_i) 
     END DO
     CLOSE(13)
     OPEN(UNIT=14, FILE=TRIM ("VAR" // ADJUSTL ( ADJUSTR ( ls_stamp_default ) ) )  // ".txt"&
             & , ACTION="write", STATUS="replace")
     DO li_i=1, fem_model % ptr_mesh % oi_n_nodes
        WRITE(14,*) fem_model % opr_global_var(:,li_i) 
     END DO
     CLOSE(14)
     ! ...

#ifdef DEBUG_TRACE
    CALL concatmsg(" Min of VAR ", ai_dtllevel = 0)
    CALL concatmsg(MINVAL(fem_model % opr_global_var(1,:)), ai_dtllevel = 0)
    CALL concatmsg(" Max of VAR ", ai_dtllevel = 0)
    CALL concatmsg(MAXVAL(fem_model % opr_global_var(1,:)), ai_dtllevel = 0)
    CALL printmsg(ai_dtllevel = 0)
#endif

    

  END SUBROUTINE RUN_MODEL
  ! ..................................................

END MODULE MODEL
