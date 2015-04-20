MODULE MODEL_DEF
  USE TypeDef
  USE SPM_DEF
  USE INDICES_DEF
  USE JOREK_PARAM_DEF
  USE MESH_DEF
  USE BLACKBOX_DEF
  USE GREENBOX_DEF
  USE QUADRATURES_DEF 
  USE LINEAR_SOLVER_DEF
  USE BASIS_DEF
  USE FEBasis
  USE FEBasis2D_Bezier
!  USE FEBasis2D_Splines
!  USE FEBasis2D_BoxSplines
  USE SPACE_DEF
  USE FIELD_DEF
  USE MATRIX_DEF
  USE FEM_DEF

  IMPLICIT NONE

  ! ... quad mesh
  TYPE(DEF_FEM_QUAD_2D)  , TARGET :: fem_model 

  TYPE(DEF_SPACE_QUAD_2D), TARGET :: space_trial
  CLASS(DEF_SPACE_QUAD_2D), POINTER :: ptr_space_trial

  TYPE(DEF_SPACE_QUAD_2D), TARGET :: space_test
  CLASS(DEF_SPACE_QUAD_2D), POINTER :: ptr_space_test
  ! ...

!  ! ... triangle mesh
!  TYPE(DEF_FEM_TRIANGLE_2D)  , TARGET :: fem_model 
!
!  TYPE(DEF_SPACE_TRIANGLE_2D), TARGET :: space_trial
!  CLASS(DEF_SPACE_TRIANGLE_2D), POINTER :: ptr_space_trial
!
!  TYPE(DEF_SPACE_TRIANGLE_2D), TARGET :: space_test
!  CLASS(DEF_SPACE_TRIANGLE_2D), POINTER :: ptr_space_test
!  ! ...

  TYPE(DEF_MESH_2D), TARGET       :: fem_mesh

  TYPE(DEF_0_FORM_2D) , TARGET  :: field_U
  CLASS(DEF_FIELD_2D), POINTER :: ptr_field

  TYPE(DEF_MATRIX_2D) , TARGET  :: matrix_Stiffnes
  CLASS(DEF_MATRIX_2D), POINTER :: ptr_system

  TYPE(DEF_LINEAR_SOLVER_2D)                 :: solver

  INTEGER :: n_var_sys 
  INTEGER :: n_var_unknown

  REAL(KIND=RK), DIMENSION(:,:), POINTER  :: Diagnostics

  INTEGER(KIND=JOREK_INTS_KIND) :: mi_nstep_max 
  
  INTEGER, PARAMETER, PRIVATE :: mi_dtllevel_base       = 0

  INTEGER, PARAMETER	      :: mi_nvar	        = 1

  INTEGER, PARAMETER          :: Matrix_A_ID              = 0

  INTEGER :: mi_mode_m1
  INTEGER :: mi_mode_n1

  REAL(KIND=JOREK_COEF_KIND)  :: mr_a
  REAL(KIND=JOREK_COEF_KIND)  :: mr_acenter
  REAL(KIND=JOREK_COEF_KIND)  :: mr_R0
  REAL(KIND=JOREK_COEF_KIND)  :: mr_Z0

END MODULE MODEL_DEF
