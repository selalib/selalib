MODULE JRK_MODEL
  ! ------------------------------ Modules used
  USE JRK_MODEL_DEF
  USE TypeDef
  USE SPM_DEF
  USE SPM
  USE Assembly
  USE Evaluator
  USE OUTPUT
  USE DIRICHLET_MOD
  USE JRK_ELEMENT_CONTRIBUTION
  USE JOREK_GLOB
  USE GRAPH
  USE NUMBERING
  USE JRK_DIAGNOSTICS_COMPUTATION
  USE JRK_TESTCASE
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
  USE JRK_TESTCASE
  ! ------------------------------
  IMPLICIT NONE

  ! .............................................
  ! ... TODO : to be removed after renaming define_model => initialize_model
  interface initialize_model
     module procedure  define_model
  end interface initialize_model
  ! .............................................

CONTAINS
  ! ..................................................

  ! ..................................................
  SUBROUTINE DEFINE_MODEL(self, use_mass_matrix, testcase_id, mpi_communicator)
    IMPLICIT NONE
    TYPE(JRK_POISSON_2D), TARGET, INTENT(INOUT) :: self
    LOGICAL, OPTIONAL, INTENT(IN) :: use_mass_matrix
    INTEGER, INTENT(IN) :: testcase_id
    INTEGER, OPTIONAL, INTENT(IN) :: mpi_communicator

    ! LOCAL
    ! ------------------------------ Local variables >> Integers
    INTEGER, PARAMETER              :: N_DIM = 2
    INTEGER                         :: li_poloidal_basis
    INTEGER                         :: li_err
    INTEGER                         :: tracelogoutput
    INTEGER                         :: tracelogdetail
    INTEGER                         :: li_nargs
    ! ------------------------------ Local variables >> Reals
    REAL(KIND=RK)                   :: T_fin
    REAL(KIND=RK)                   :: T_deb
    REAL(KIND=RK), DIMENSION(1:2)   :: lpr_scale
    ! ------------------------------ Local variables >> Logical
    LOGICAL                         :: ll_stdoutput
    LOGICAL                         :: ll_use_mass_matrix
    ! ------------------------------ Local variables >> Characters
    CHARACTER(LEN = 1024)           :: dirname
    CHARACTER(LEN = 1024)           :: argname
    CHARACTER(LEN = 1024)           :: django_parameters
    CHARACTER(LEN = 1024)           :: filename_parameter
    CHARACTER(LEN = 1024)           :: progname
    CHARACTER(LEN = 1024)           :: rootname
    CHARACTER(LEN = 1024)           :: filename
    CHARACTER(LEN = 1024)           :: exec
    ! ..................................................
    !  print *, "DEFINE_MODEL : BEGIN"

    ll_use_mass_matrix = .FALSE.
    IF (PRESENT(use_mass_matrix)) THEN
       ll_use_mass_matrix = use_mass_matrix
    END IF

    argname = "--django-params"

    CALL JOREK_GET_ARGUMENTS(argname, django_parameters, li_err)

    if (li_err .ne. 0) then
       print *, "Error in comand line:" // &
            "I was expecting an argument '--django-params $DJANGO_PARAMS' " // &
            " where $DJANGO_PARAMS contains the files parameter.txt " // &
            " and solver.txt in order to initialize the DJANGO solver"
    end if
    filename_parameter = trim(django_parameters) // "/parameter.txt"

    CALL INITIALIZE_JOREK_PARAMETERS(filename_parameter)

#ifdef DEBUG_TRACE
    CALL JOREK_PARAM_GETINT(INT_TRACE_DETAIL_ID,TraceLogDetail,li_err)
    CALL JOREK_PARAM_GETINT(INT_TRACE_OUTPUT_ID,TraceLogOutput,li_err)
    ll_stdoutput = (TraceLogOutput==1)
    CALL opentracelog(al_stdoutput = ll_stdoutput, ai_dtllevel = TraceLogDetail)

    CALL printlog("=====================================================", ai_dtllevel = 0)
#ifdef PETSC_ENABLED
    CALL printlog("== PETSC ENABLED                                   ==", ai_dtllevel = 0)
#else
    CALL printlog("== PETSC DISABLED                                  ==", ai_dtllevel = 0)
#endif

#ifdef DEBUG_ELEMENT
    CALL printlog("== Jorek is working in Debug Element Assembly mode ==", ai_dtllevel = 0)
#endif
    CALL printlog("=====================================================", ai_dtllevel = 0)
#endif

    CALL getarg(0, exec)
    rootname = "output"

    PRINT *, "RunName : ", Trim(exec)

    IF (PRESENT(mpi_communicator)) THEN
       call SPM_INITIALIZE(li_err, mpi_communicator=mpi_communicator)
    ELSE
       call SPM_INITIALIZE(li_err)
    END IF

    ! ... define model parameters
    ! ..................................................
    CALL JOREK_Param_GETInt(INT_MODES_M1_ID , mi_mode_m1  , li_err)
    CALL JOREK_Param_GETInt(INT_MODES_N1_ID , mi_mode_n1  , li_err)
    CALL JOREK_Param_GETInt(INT_NSTEP_MAX_ID, mi_nstep_max, li_err)
    CALL JOREK_Param_GETInt(INT_TYPEMESH_ID , mi_mesh_type, li_err)
    CALL JOREK_Param_GETInt(INT_DIAGNOSTICS_FREQSAVE_ID, mi_freqsave, li_err)

    CALL JOREK_Param_GETReal(REAL_RGEO_ID    , mr_R0       , li_err)
    CALL JOREK_Param_GETReal(REAL_ZGEO_ID    , mr_Z0       , li_err)
    CALL JOREK_Param_GETReal(REAL_AMIN_ID    , mr_a        , li_err) 
    CALL JOREK_Param_GETReal(REAL_ACENTER_ID , mr_acenter  , li_err)
    CALL JOREK_Param_GETReal(REAL_RSCALE_ID  , mr_rscale       , li_err)
    CALL JOREK_Param_GETReal(REAL_ZSCALE_ID  , mr_zscale       , li_err)
    ! ...

    self % use_mass_matrix = ll_use_mass_matrix

    !    if the basis is given in a direcroty dirname given by the argument argname
    IF (mi_mesh_type == INT_MESH_BEZIER_DESCRIPTION) THEN
       argname = "--basis"
       CALL JOREK_GET_ARGUMENTS(argname, dirname, li_err)
       CALL SPACE_CREATE(self % space_trial, self % fem_mesh, dirname=dirname)
       ptr_space_trial => self % space_trial
       CALL SPACE_CREATE(self % space_test, self % fem_mesh, dirname=dirname)
       ptr_space_test => self % space_test
    ELSE
       !    if the basis is defined internally (hermite-bezier for example)
       CALL SPACE_CREATE(self % space_trial, self % fem_mesh)
       ptr_space_trial => self % space_trial
       CALL SPACE_CREATE(self % space_test, self % fem_mesh)
       ptr_space_test => self % space_test
       ! ...
    END IF

    ! ...
    argname = "--geometry"
    CALL JOREK_GET_ARGUMENTS(argname, dirname, li_err)

    CALL CREATE_MESH(self % fem_mesh, dirname=dirname)
    ! ...
    lpr_scale(1) = mr_rscale
    lpr_scale(2) = mr_zscale 



    CALL MESH_SCALE(self % fem_mesh,lpr_scale)

    CALL MODEL_CREATE(self % fem_model, self % fem_mesh, ptr_space_trial, ptr_space_test)
    ptr_fem_model => self % fem_model

    CALL FIELD_CREATE(self % field_U, self % fem_model % space_trial, field_name="density")
    ptr_field => self % field_U
    CALL MODEL_APPEND_FIELD(self % fem_model, ptr_field)

    ! ... STIFFNESS MATRIX
    CALL MATRIX_CREATE(self % matrix_stiffness)
    ptr_system => self % matrix_stiffness

    CALL MATRIX_APPEND_UNKNOWN_FIELD(ptr_system, ptr_field)
    CALL MODEL_APPEND_MATRIX(self % fem_model, ptr_system)

    ! ... define weak formulation 
    ptr_system % ptr_matrix_contribution => Stiffness_Matrix_for_Vi_Vj

    ! .....................................
    !        If using Mass Matrix
    ! .....................................
    IF (self % use_mass_matrix) THEN

       ! ... MASS MATRIX
       CALL MATRIX_CREATE(self % matrix_mass)
       ptr_mass_system => self % matrix_mass

       CALL MATRIX_APPEND_UNKNOWN_FIELD(ptr_mass_system, ptr_field)
       CALL MODEL_APPEND_MATRIX(self % fem_model, ptr_mass_system)

       ! ... define weak formulation
       ptr_mass_system % ptr_matrix_contribution => Mass_Matrix_for_Vi_Vj
    END IF
    ! .....................................
    CALL MODEL_INITIALIZE(self % fem_model)
    ! ... define testcase
    CALL SET_TESTCASE_MODEL(self, testcase_id)

    !    ! ... Assemble Stiffness matrix 
    !    ptr_system => self % matrix_stiffness
    !    CALL CPU_TIME(T_deb)
    !    CALL MODEL_ASSEMBLY(self % fem_model, ptr_system)
    !    CALL CPU_TIME(T_fin)
    !    ! ...
    !
    !    ! .....................................
    !    !        If using Mass Matrix
    !    ! .....................................
    !    IF (self % use_mass_matrix) THEN
    !       ! ... Assemble Mass matrix 
    !       ptr_mass_system => self % matrix_mass
    !       CALL CPU_TIME(T_deb)
    !       CALL MODEL_ASSEMBLY(self % fem_model, ptr_mass_system)
    !       CALL CPU_TIME(T_fin)
    !       ! ...
    !    END IF
    !    ! .....................................

    ! ... create solver for stiffness matrix 
    ptr_system => self % matrix_stiffness
    filename = trim(django_parameters) // "/solver.txt"
    CALL LINEAR_SOLVER_CREATE(self % solver, ptr_system, filename)

    ! ...

    ! .....................................
    !        If using Mass Matrix
    ! .....................................
    IF (self % use_mass_matrix) THEN
       ! ... create solver for mass matrix 
       ptr_mass_system => self % matrix_mass
       CALL LINEAR_SOLVER_CREATE(self % mass_solver, ptr_mass_system, filename)
       ! ...
    END IF
    ! .....................................

    self % fem_model % opr_diagnostics(:, :) = 0.0
    g_diagnostics => self % fem_model % opr_diagnostics

    !    print *, "DEFINE_MODEL : END"
  END SUBROUTINE DEFINE_MODEL
  ! ..................................................

  ! ..................................................
  SUBROUTINE FREE_MODEL(self )
    IMPLICIT NONE
    TYPE(JRK_POISSON_2D), INTENT(INOUT) :: self 

    CALL MODEL_FREE(self % fem_model)
    CALL LINEAR_SOLVER_FREE(self % solver)
    IF (self % use_mass_matrix) THEN
       CALL LINEAR_SOLVER_FREE(self % mass_solver)
    END IF
  END SUBROUTINE FREE_MODEL
  ! ..................................................

  ! ..................................................
  SUBROUTINE ASSEMBLE_MODEL(self)
    IMPLICIT NONE
    TYPE(JRK_POISSON_2D), TARGET, INTENT(INOUT) :: self 
    REAL(KIND=RK)                   :: T_fin
    REAL(KIND=RK)                   :: T_deb

    ! ... Assemble Stiffness matrix 
    ptr_system => self % matrix_stiffness
    CALL CPU_TIME(T_deb)
    CALL MODEL_ASSEMBLY(self % fem_model, ptr_system)
    CALL CPU_TIME(T_fin)

    ! ...

    ! .....................................
    !        If using Mass Matrix
    ! .....................................
    IF (self % use_mass_matrix) THEN
       ! ... Assemble Mass matrix 
       ptr_mass_system => self % matrix_mass
       CALL CPU_TIME(T_deb)
       CALL MODEL_ASSEMBLY(self % fem_model, ptr_mass_system)
       CALL CPU_TIME(T_fin)
       ! ...
    END IF
    ! .....................................
  END SUBROUTINE ASSEMBLE_MODEL
  ! ..................................................

  ! ..................................................
  SUBROUTINE SET_TESTCASE_MODEL(self, testcase_id)
    IMPLICIT NONE
    TYPE(JRK_POISSON_2D), INTENT(INOUT) :: self 
    INTEGER, INTENT(IN) :: testcase_id
    ! LOCAL
    INTEGER :: ierr

    ! ... Stiffness matrix testcases
    SELECT CASE (testcase_id)

    CASE (1)
       ptr_system % ptr_rhs_contribution => RHS_for_Vi_1
       ptr_field  % ptr_func_analytical  => ANALYTICAL_MODEL_1
    CASE (2)
       ptr_system % ptr_rhs_contribution => RHS_for_Vi_2
       ptr_field  % ptr_func_analytical  => ANALYTICAL_MODEL_2
    CASE (3)
       ptr_system % ptr_rhs_contribution => RHS_for_Vi_3
       ptr_field  % ptr_func_analytical  => ANALYTICAL_MODEL_3
    CASE (4)
       ptr_system % ptr_rhs_contribution => RHS_for_Vi_4
       ptr_field  % ptr_func_analytical  => ANALYTICAL_MODEL_4
    CASE (5)
       ptr_system % ptr_rhs_contribution => RHS_for_Vi_5
       ptr_field  % ptr_func_analytical  => ANALYTICAL_MODEL_5
    CASE (6)
       ptr_system % ptr_rhs_contribution => RHS_for_Vi_6
       ptr_field  % ptr_func_analytical  => ANALYTICAL_MODEL_6
    CASE (100)
       ptr_system % ptr_rhs_contribution => RHS_for_Vi_100
       ptr_field  % ptr_func_analytical  => ANALYTICAL_MODEL_100
    CASE DEFAULT

       ptr_system % ptr_rhs_contribution => RHS_for_Vi_ONE
       ptr_field  % ptr_func_analytical  => ANALYTICAL_MODEL_ZERO
    END SELECT
    ! ...

    ! .....................................
    !        If using Mass Matrix
    ! .....................................
    IF (self % use_mass_matrix) THEN

       SELECT CASE (testcase_id)
       CASE (10)
          ptr_mass_system % ptr_rhs_contribution => NULL()
          ptr_field  % ptr_func_analytical       => ANALYTICAL_MODEL_MASS_1
          !          ptr_field  % ptr_func_analytical       => ANALYTICAL_MODEL_ZERO
       CASE (11)
          ptr_mass_system % ptr_rhs_contribution => NULL()
          ptr_field  % ptr_func_analytical       => ANALYTICAL_MODEL_MASS_2
       CASE (12)
          ptr_mass_system % ptr_rhs_contribution => NULL()
          ptr_field  % ptr_func_analytical       => ANALYTICAL_MODEL_MASS_3
       CASE DEFAULT
          ptr_mass_system % ptr_rhs_contribution => NULL() 
          ptr_field  % ptr_func_analytical       => ANALYTICAL_MODEL_ZERO
       END SELECT
    END IF
    ! .....................................

    CALL JOREK_Param_SETInt(INT_TESTCASE_ID , testcase_id , ierr)

  END SUBROUTINE SET_TESTCASE_MODEL
  ! ..................................................

  ! ..................................................
  SUBROUTINE DIAGNOSTICS_MODEL(self, nstep )
    IMPLICIT NONE
    TYPE(JRK_POISSON_2D), INTENT(INOUT) :: self 

    INTEGER                    :: nstep
    !LOCAL
    INTEGER(KIND=SPM_INTS_KIND)     :: nRows ! Number of Rows
    INTEGER                         :: li_err
    INTEGER                         :: li_i   
    INTEGER                         :: li_elmt
    REAL(KIND=RK)                   :: T_fin, T_deb
    REAL(KIND=RK), DIMENSION(1:1,1:2) :: lpr_x
    REAL(KIND=RK), DIMENSION(1:1,1:2) :: lpr_v
    REAL(KIND=RK), DIMENSION(:)  , POINTER  :: lpr_Uu
    CHARACTER(LEN = 1024)           :: filename
    CHARACTER(LEN = 1024)           :: argname
    CHARACTER(len=20)               :: ls_stamp_default
    CHARACTER(len=20)               :: ls_msg
    CHARACTER(LEN = 1024)           :: exec
    CLASS(DEF_MATRIX_2D), POINTER   :: ptr_current_system

    ! ...

    ptr_current_system => ptr_system 
    IF (self % use_mass_matrix) THEN  
       ptr_current_system => ptr_mass_system 
    END IF
    ! ...

    !print *, "before diags", mpo_M(0) % oi_nR, "virtual",mpo_M(0) % oi_nR_virtual

    ptr_current_system % ptr_assembly_diagnostics => Assembly_Diagnostics!_Field_Evaluate
       
    CALL MODEL_DIAGNOSTICS(self % fem_model, Plot_Diagnostics, ptr_current_system, nstep)

    !TODO : change filename to a name that wont be erased
    CALL getarg(0, filename)
    filename = TRIM(filename) 
    CALL MODEL_SAVE(self % fem_model, ptr_current_system, filename, 0)
    ! ...

  END SUBROUTINE DIAGNOSTICS_MODEL
  ! ..................................................

  ! ..................................................
  SUBROUTINE GET_SHAPE_MODEL(self, arr_shape )
    IMPLICIT NONE
    TYPE(JRK_POISSON_2D), INTENT(IN) :: self 
    INTEGER, DIMENSION(:), INTENT(INOUT) :: arr_shape
    !LOCAL
    INTEGER :: n_rows
    INTEGER :: n_cols

    CALL GET_NR_MATRIX(ptr_system, n_rows)
    CALL GET_NC_MATRIX(ptr_system, n_cols)

    arr_shape(1) = n_rows
    arr_shape(2) = n_cols

  END SUBROUTINE GET_SHAPE_MODEL
  ! ..................................................

  ! ..................................................
  SUBROUTINE GET_RHS_MODEL(self, arr_rhs, use_mass_matrix)
    IMPLICIT NONE
    TYPE(JRK_POISSON_2D), INTENT(IN) :: self 
    REAL(KIND=RK), DIMENSION(:), INTENT(INOUT) :: arr_rhs
    LOGICAL, OPTIONAL, INTENT(IN) :: use_mass_matrix
    !LOCAL
    LOGICAL                         :: ll_use_mass_matrix

    ll_use_mass_matrix = .FALSE.
    IF (PRESENT(use_mass_matrix)) THEN
       ll_use_mass_matrix = use_mass_matrix   
    END IF

    IF (ll_use_mass_matrix) THEN
       arr_rhs = ptr_mass_system % opr_global_rhs
    ELSE
       arr_rhs = ptr_system % opr_global_rhs
    END IF
  END SUBROUTINE GET_RHS_MODEL
  ! ..................................................

  ! ..................................................
  SUBROUTINE SET_RHS_MODEL(self, arr_rhs, use_mass_matrix)
    IMPLICIT NONE
    TYPE(JRK_POISSON_2D), INTENT( INOUT )    :: self
    REAL(KIND=RK), DIMENSION(:), INTENT(IN) :: arr_rhs
    LOGICAL, OPTIONAL, INTENT(IN) :: use_mass_matrix
    !LOCAL
    LOGICAL                         :: ll_use_mass_matrix

    ll_use_mass_matrix = .FALSE.
    IF (PRESENT(use_mass_matrix)) THEN
       ll_use_mass_matrix = use_mass_matrix   
    END IF

    IF (ll_use_mass_matrix) THEN
       ptr_mass_system % opr_global_rhs = arr_rhs 
    ELSE
       ptr_system % opr_global_rhs = arr_rhs 
    END IF


  END SUBROUTINE SET_RHS_MODEL
  ! ..................................................

  ! ..................................................
  subroutine solve_model(self)
    IMPLICIT NONE
    TYPE(JRK_POISSON_2D)              , INTENT( INOUT )    :: self
    !    print *, "hello from solve_model"
    ! ...  Solver call
    CALL LINEAR_SOLVER_SOLVE(self % solver, &
         & ptr_system % opr_global_rhs,  &
         & ptr_system % opr_global_unknown) 

    ! print *,  "MIN Rho : ", MINVAL(ptr_system % opr_global_rhs), &
    !         & "MAX Rho : ", MAXVAL(ptr_system % opr_global_rhs)

    ! print *,  "MIN Phi : ", MINVAL(ptr_system % opr_global_unknown), &
    !         & "MAX Phi : ", MAXVAL(ptr_system % opr_global_unknown)
    ! ! ...
    ! print *, "linear_solver_solve ok"
    ! ... update related fields to the linear system
    CALL UPDATE_FIELDS(self % fem_model, ptr_system)
    !   print *, "fields updated"

    ! ...
  end subroutine solve_model
  ! ..................................................

  ! ..................................................
  subroutine project_model(self)
    IMPLICIT NONE
    TYPE(JRK_POISSON_2D)              , INTENT( INOUT )    :: self

    ! .....................................
    !        If using Mass Matrix
    ! .....................................
    IF (self % use_mass_matrix) THEN
       ! ...  Solver call
       CALL LINEAR_SOLVER_SOLVE(self % mass_solver, &
            & ptr_mass_system % opr_global_rhs,  &
            & ptr_mass_system % opr_global_unknown) 

       print *,  "MIN Y : ", MINVAL(ptr_mass_system % opr_global_rhs), &
            & "MAX Y : ", MAXVAL(ptr_mass_system % opr_global_rhs)

       print *,  "MIN X : ", MINVAL(ptr_mass_system % opr_global_unknown), &
            & "MAX X : ", MAXVAL(ptr_mass_system % opr_global_unknown)
       ! ...

       ! ... update related fields to the linear system
       CALL UPDATE_FIELDS(self % fem_model, ptr_mass_system)

       ! ...
    END IF
    ! .....................................
  end subroutine project_model
  ! ..................................................

  ! ..................................................
  subroutine add_contribution_field_model(self, element_id, logical_position, weight, contributions)
    IMPLICIT NONE
    TYPE(JRK_POISSON_2D), INTENT( INOUT )         :: self
    INTEGER                  , INTENT(IN)         :: element_id 
    REAL(KIND=RK), DIMENSION(:,:), INTENT(IN)     :: logical_position
    REAL(KIND=RK), INTENT(IN)                     :: weight
    REAL(KIND=RK), DIMENSION(:) , INTENT(INOUT)   :: contributions
    ! LOCAL

    ! adding contribution in the table "lpr_contrib"
    ! routine in jorek/src/fem
    CALL ADD_CONTRIBUTION_FIELD(ptr_fem_model, &
         & ptr_system,                         &
         & element_id,                         &
         & logical_position,                   &
         & weight,                             &
         & contributions)

  end subroutine add_contribution_field_model
  ! ..................................................

  ! ..................................................
  subroutine evaluate_field_model_logical(self, element_id, logical_position, values)
    IMPLICIT NONE
    TYPE(JRK_POISSON_2D)      , INTENT( INOUT )    :: self
    INTEGER                  , INTENT(IN)         :: element_id 
    REAL(KIND=RK), DIMENSION(:,:) , INTENT(IN)    :: logical_position
    REAL(KIND=RK), DIMENSION(:,:) , INTENT(OUT)   :: values

    CALL FIELD_EVALUATE_LOGICAL(self%field_U, element_id, logical_position, values)

  end subroutine evaluate_field_model_logical
  ! ..................................................
  ! ..................................................
  subroutine evaluate_field_model(self, element_id, logical_position, values)
    IMPLICIT NONE
    TYPE(JRK_POISSON_2D)      , INTENT( INOUT )    :: self
    INTEGER                  , INTENT(IN)         :: element_id 
    REAL(KIND=RK), DIMENSION(:,:) , INTENT(IN)    :: logical_position
    REAL(KIND=RK), DIMENSION(:,:) , INTENT(OUT)   :: values

    CALL FIELD_EVALUATE(self%field_U, element_id, logical_position, values)

  end subroutine evaluate_field_model
  ! ..................................................

  ! ..................................................
  SUBROUTINE ANALYTICAL_MODEL_ZERO(self, apr_x, apr_v,apr_info,api_info,n_variable,n_dimension)
    ! .. TODO to remove from here, must have a default function set to the zero function for every field
    IMPLICIT NONE
    CLASS(DEF_FIELD_2D), POINTER                       :: self
    INTEGER                                            :: n_variable
    INTEGER                                            :: n_dimension
    REAL(KIND=8), DIMENSION(n_dimension)               :: apr_x
    REAL(KIND=8), DIMENSION(N_variable, n_dimension+1) :: apr_v
    REAL(KIND=8), DIMENSION(10)                        :: apr_info
    INTEGER, DIMENSION(10)                             :: api_info
    ! LOCAL

    apr_v = 0.0

  END SUBROUTINE ANALYTICAL_MODEL_ZERO
  ! ..................................................

END MODULE JRK_MODEL
