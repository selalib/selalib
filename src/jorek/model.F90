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
  USE MATRIX
  IMPLICIT NONE

CONTAINS
  ! ..................................................
  SUBROUTINE INITIALIZE_MODEL( )
  IMPLICIT NONE
    INTEGER                     :: ierr
    INTEGER(KIND=SPM_INTS_KIND) :: nRows ! Number of Rows
    INTEGER(KIND=SPM_INTS_KIND) :: nCols ! Number of Columns
    INTEGER(KIND=SPM_INTS_KIND) :: li_locsize
    INTEGER                     :: li_ncolors
    INTEGER                     :: li_i
    INTEGER, DIMENSION(:), ALLOCATABLE  :: lpi_rvars
    INTEGER, DIMENSION(:), ALLOCATABLE  :: lpi_cvars
    INTEGER                     :: li_n_nodes_global
    INTEGER, PARAMETER          :: N_DIM = 2

    ALLOCATE(lpi_rvars(0:n_var_sys))
    ALLOCATE(lpi_cvars(0:n_var_sys))

    lpi_rvars(0) = n_var_sys
    DO li_i = 1, n_var_sys
       lpi_rvars(li_i) = li_i
    END DO

    lpi_cvars(0) = n_var_sys
    DO li_i = 1, n_var_sys
       lpi_cvars(li_i) = li_i
    END DO

    CALL INITIALIZE_MATRIX(Matrix_A_ID, Mesh2D, lpi_rvars, lpi_cvars)

    ! ... Get the new size of the matrix 
    CALL SPM_GetnR(Matrix_A_ID, nRows, ierr)
    CALL SPM_GetnC(Matrix_A_ID, nCols, ierr)
    ! ...

    li_n_nodes_global = Mesh2D % oi_n_nodes

    ALLOCATE(Global_Rhs(nCols))
    ALLOCATE(Global_Unknown(nCols))

    ALLOCATE(Var(Mesh2D% oi_n_max_order * n_var_unknown, li_n_nodes_global) )
    Allocate(Diagnostics(N_diag,mi_nstep_max+1))

    Global_Rhs     = 0.0
    Global_Unknown = 0.0
    Var            = 0.0
    Diagnostics    = 0.0

    CALL CREATE_GREENBOX(GBox2D, &
            & n_var_unknown, &
            & n_var_sys, &
            & Mesh2D % ptr_quad % oi_n_points)

  END SUBROUTINE INITIALIZE_MODEL
  ! ..................................................

  ! ..................................................
  SUBROUTINE FREE_MODEL( )
  IMPLICIT NONE
    DEALLOCATE(Global_Rhs)
    DEALLOCATE(Global_Unknown)
    DEALLOCATE(Var)
    DEALLOCATE(Diagnostics)

    CALL FREE_GREENBOX(GBox2D)
  END SUBROUTINE FREE_MODEL
  ! ..................................................

  ! ..................................................
  SUBROUTINE RUN_MODEL( )
  IMPLICIT NONE
     REAL(KIND=RK)    :: T_fin, T_deb
     CHARACTER(LEN = 1024)           :: filename
     INTEGER(KIND=SPM_INTS_KIND) :: nRows ! Number of Rows
     INTEGER                     :: ierr
     INTEGER                     :: li_myRank
     character(len=20) :: ls_stamp_default
     character(len=20)   :: ls_msg
     REAL(KIND=RK), DIMENSION(:), ALLOCATABLE  :: lpr_Y
     INTEGER, DIMENSION(:), ALLOCATABLE  :: lpi_rvars
     INTEGER, DIMENSION(:), ALLOCATABLE  :: lpi_cvars
     INTEGER                     :: li_i   

     ALLOCATE(lpi_rvars(0:n_var_sys))
     ALLOCATE(lpi_cvars(0:n_var_sys))
   
     lpi_rvars(0) = n_var_sys
     DO li_i = 1, n_var_sys
        lpi_rvars(li_i) = li_i
     END DO
   
     lpi_cvars(0) = n_var_sys
     DO li_i = 1, n_var_sys
        lpi_cvars(li_i) = li_i
     END DO

!     Find out the rank of the current proc
     li_myRank = 0
#ifdef MPI_ENABLED
      CALL MPI_Comm_rank ( MPI_COMM_WORLD, li_myRank, ierr)
#endif

     CALL SPM_GetnR(Matrix_A_ID, nRows, ierr)

     write(ls_msg,*) li_myRank
     ls_stamp_default = "-proc_" // TRIM ( ADJUSTL ( ls_msg ) )
     ls_stamp_default = TRIM ( ADJUSTL ( ADJUSTR ( ls_stamp_default ) ) )

     CALL getarg(0, filename)

     ! ... Loop over elements 
     CALL CPU_TIME(T_deb)
     CALL Loop_On_Elmts(Matrix_A_ID, li_myRank, &
             & Mesh2D, GBox2D, &
             & RHS_for_Vi, Matrix_for_Vi_Vj, &
             & Global_Unknown, Var, Global_Rhs, &
             & lpi_rvars, lpi_cvars)
     CALL CPU_TIME(T_fin)
     ! ...

!     print *, "RHS ", Global_Rhs(1:6)

     ! ... example of solver calls
     CALL LINEAR_SOLVER_NEW_WITH_MATRIX_ID(Solver, Matrix_A_ID)
     CALL LINEAR_SOLVER_SOLVE(Solver, Global_Rhs, Global_Unknown) 
     CALL LINEAR_SOLVER_FREE(Solver)
     ! ...

     ALLOCATE(lpr_Y(nRows))
     lpr_Y=0.d0 
     CALL SPM_MATMULT(Matrix_A_ID, Global_Unknown, lpr_Y, ierr)

     PRINT *, MAXVAL(lpr_Y-Global_Rhs), MINVAL(lpr_Y-Global_Rhs) 

!     print *, "UNKNOWN ", Global_Unknown

     ! ... Evaluate unknowns on vertecies and compte model norms 
     CALL Evaluate_On_Elmts(Matrix_A_ID, li_myRank, &
             & Mesh2D, GBox2D, &
             & Global_Unknown, Var, ANALYTICAL_MODEL, &
             & Assembly_Diagnostics,Plot_Diagnostics, &
             & lpi_rvars, lpi_cvars,0)

      filename = TRIM(filename) 
     CALL SaveMeshes(Mesh2D, GBox2D, Var, ANALYTICAL_MODEL, filename,0)

     ! ...
     OPEN(UNIT=12, FILE=TRIM ("RHS" // ADJUSTL ( ADJUSTR ( ls_stamp_default ) ) )  // ".txt"&
             & , ACTION="write", STATUS="replace")
     DO li_i=1,nRows
        WRITE(12,*) Global_Rhs(li_i) 
     END DO
     CLOSE(12)
     OPEN(UNIT=13, FILE=TRIM ("UNKNOWN" // ADJUSTL ( ADJUSTR ( ls_stamp_default ) ) )  // ".txt"&
             & , ACTION="write", STATUS="replace")
     DO li_i=1,nRows
        WRITE(13,*) Global_Unknown(li_i) 
     END DO
     CLOSE(13)
     OPEN(UNIT=14, FILE=TRIM ("VAR" // ADJUSTL ( ADJUSTR ( ls_stamp_default ) ) )  // ".txt"&
             & , ACTION="write", STATUS="replace")
     DO li_i=1, Mesh2D % oi_n_nodes
        WRITE(14,*) Var(:,li_i) 
     END DO
     CLOSE(14)
     ! ...

#ifdef DEBUG_TRACE
    CALL concatmsg(" Min of Variables ", ai_dtllevel = 0)
    CALL concatmsg(MINVAL(VAR(1,:)), ai_dtllevel = 0)
    CALL concatmsg(" Max of Variables ", ai_dtllevel = 0)
    CALL concatmsg(MAXVAL(VAR(1,:)), ai_dtllevel = 0)
    CALL printmsg(ai_dtllevel = 0)
#endif

    

  END SUBROUTINE RUN_MODEL
  ! ..................................................

END MODULE MODEL
