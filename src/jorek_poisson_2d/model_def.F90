MODULE MODEL_DEF
  USE TypeDef
  USE SPM_DEF
  USE INDICES_DEF
  USE JOREK_CONTEXT_DEF
  USE JOREK_CONTEXT
  USE JOREK_PARAM_DEF
  USE JOREK_PARAM
  IMPLICIT NONE

  TYPE(Def_Node)              :: TheNodes
  TYPE(Def_Elmt2D)            :: Elmts2D

  REAL(KIND=RK), DIMENSION(:), ALLOCATABLE  :: Global_Rhs 
  REAL(KIND=RK), DIMENSION(:), ALLOCATABLE  :: Global_Unknown
  REAL(KIND=RK), DIMENSION(:,:), ALLOCATABLE       :: Var

  INTEGER(KIND=JOREK_INTS_KIND) :: mi_nstep_max 
  
!  INTEGER, PARAMETER          :: Matrix_ID              = 0
  INTEGER, PARAMETER, PRIVATE :: mi_dtllevel_base       = 0

  INTEGER, PARAMETER	      :: mi_nvar	        = 1

  INTEGER, PARAMETER          :: Matrix_A_ID              = 0

  TYPE(CONTEXT_DEF) :: CONTEXT
  INTEGER :: mi_mode_m1
  INTEGER :: mi_mode_n1

  REAL(KIND=JOREK_COEF_KIND)  :: mr_a
  REAL(KIND=JOREK_COEF_KIND)  :: mr_acenter
  REAL(KIND=JOREK_COEF_KIND)  :: mr_R0
  REAL(KIND=JOREK_COEF_KIND)  :: mr_Z0

CONTAINS
  ! ..................................................
  SUBROUTINE DEFINE_MODEL( )
  IMPLICIT NONE
    INTEGER :: ierr

    n_dim               = 2

    current_model       = 1

    n_var_unknown       = 1
    n_var_sys           = 1
    
    i_Vp_Rho            = 1
    

    nmatrices           = 1

    i_Vu_Rho            = 1 
    ALLOCATE(NamesVarU(n_var_unknown)) 
    NamesVarU(i_Vu_Rho)  = "Density"

    ALLOCATE(NamesVarP(n_var_unknown)) 
    NamesVarP(i_Vp_Rho)  = "Density"

    CALL JOREK_Param_GETInt(INT_MODES_M1_ID, mi_mode_m1, ierr)
    CALL JOREK_Param_GETInt(INT_MODES_N1_ID, mi_mode_n1, ierr)

    CALL JOREK_Param_GETReal(REAL_RGEO_ID,mr_R0,ierr)
    CALL JOREK_Param_GETReal(REAL_ZGEO_ID,mr_Z0,ierr)
    CALL JOREK_Param_GETReal(REAL_AMIN_ID,mr_a,ierr) 
    CALL JOREK_Param_GETReal(REAL_ACENTER_ID,mr_acenter,ierr)
    CALL JOREK_Param_GETInt(INT_NSTEP_MAX_ID,mi_nstep_max,ierr)

  END SUBROUTINE DEFINE_MODEL
  ! ..................................................

END MODULE MODEL_DEF
