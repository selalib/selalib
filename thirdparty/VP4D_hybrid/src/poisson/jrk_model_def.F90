MODULE JRK_MODEL_DEF
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
  USE SPACE_DEF
  USE FIELD_DEF
  USE MATRIX_DEF
  USE FEM_DEF

  IMPLICIT NONE

  REAL(KIND=RK), DIMENSION(:,:), POINTER :: g_diagnostics 

  ! ... quad mesh
  TYPE, PUBLIC :: JRK_POISSON_2D 
     TYPE(DEF_FEM_QUAD_2D)      :: fem_model 
     TYPE(DEF_SPACE_QUAD_2D)    :: space_trial
     TYPE(DEF_SPACE_QUAD_2D)    :: space_test
     TYPE(DEF_MESH_2D)          :: fem_mesh
     TYPE(DEF_0_FORM_2D)        :: field_U
     TYPE(DEF_MATRIX_2D)        :: matrix_stiffness
     TYPE(DEF_LINEAR_SOLVER_2D) :: solver
     LOGICAL                    :: use_mass_matrix=.FALSE.
     TYPE(DEF_MATRIX_2D)        :: matrix_mass
     TYPE(DEF_LINEAR_SOLVER_2D) :: mass_solver
  END TYPE JRK_POISSON_2D
  ! ...

  ! ...
  CLASS(DEF_FEM_QUAD_2D), POINTER           :: ptr_fem_model
  CLASS(DEF_SPACE_QUAD_2D), POINTER         :: ptr_space_trial
  CLASS(DEF_SPACE_QUAD_2D), POINTER         :: ptr_space_test
  CLASS(DEF_FIELD_2D), POINTER              :: ptr_field
  CLASS(DEF_MATRIX_2D), POINTER             :: ptr_system
  CLASS(DEF_MATRIX_2D), POINTER             :: ptr_mass_system
  ! ...

  INTEGER(KIND=JOREK_INTS_KIND)             :: mi_nstep_max 
  
  INTEGER, PARAMETER, PRIVATE               :: mi_dtllevel_base       = 0

  INTEGER                                   :: mi_mesh_type

  INTEGER                                   :: mi_mode_m1
  INTEGER                                   :: mi_mode_n1
  INTEGER                                   :: mi_freqsave

  REAL(KIND=JOREK_COEF_KIND)                :: mr_a
  REAL(KIND=JOREK_COEF_KIND)                :: mr_acenter
  REAL(KIND=JOREK_COEF_KIND)                :: mr_R0
  REAL(KIND=JOREK_COEF_KIND)                :: mr_Z0
  REAL(KIND=JOREK_COEF_KIND)                :: mr_e_l2_norm
  REAL(KIND=JOREK_COEF_KIND)                :: mr_rscale
  REAL(KIND=JOREK_COEF_KIND)                :: mr_zscale
END MODULE JRK_MODEL_DEF
