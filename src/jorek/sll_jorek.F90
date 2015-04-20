MODULE jorek_model
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

END MODULE jorek_model

module sll_jorek

  ! ------------------------------ Modules used
  USE TypeDef
  USE FEBasis
  USE tracelog_module
  USE SPM_DEF
  USE JOREK_MODEL
  USE JOREK_PARAM
  USE JOREK_PARAM_DEF
  USE COORDINATES_DEF
  USE COORDINATES
  USE FIELD_DEF
  USE JOREK_GLOB_DEF
  USE INDICES_DEF
  USE MODEL_DEF
  use fem_def
  use fem

  ! ------------------------------ 
  IMPLICIT NONE

  type, public :: sll_jorek_solver

  end type sll_jorek_solver

  interface sll_create
  module procedure initialize_jorek
  end interface sll_create

contains

  subroutine initialize_jorek(jorek)

  type(sll_jorek_solver) :: jorek

  ! ------------------------------ Local variables >> Integers
  INTEGER          :: nb_args
  INTEGER          :: is, ie, il, i1, i2, is1, is2
  INTEGER          :: imesh=1, Nrefine
  INTEGER          :: iv, iv_Pol, ipol, Ns_3D
  INTEGER          :: ierr
  INTEGER          :: TraceLogDetail 
  INTEGER          :: TraceLogOutput

  INTEGER, DIMENSION(:,:), POINTER     :: Nu
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: Ok
  INTEGER(KIND=SPM_INTS_KIND) :: nRows !NUMBER OF ROWS
  INTEGER(KIND=SPM_INTS_KIND) :: nCols !NUMBER OF COLUMNS

  ! ------------------------------ Local variables >> Reals

  REAL(KIND=RK), DIMENSION(:,:), POINTER :: Coor  

  ! ------------------------------ Local variables >> Characters
  CHARACTER(LEN = 1024)         :: exec
  CHARACTER(LEN = 1024)         :: rootname
  CHARACTER(LEN = *), PARAMETER :: mod_name = "Main"

!  ! ------------------------------ Local variables >> Structures
  LOGICAL :: ll_stdoutput
  INTEGER :: li_assembly_proc
  INTEGER :: li_toroidal_basis
  character filename_parameter*1024
  INTEGER :: li_n_Gauss_Rp, li_n_Gauss_Zp

    INTEGER, PARAMETER          :: N_DIM = 2
    INTEGER :: li_poloidal_basis
    INTEGER :: li_err
    CHARACTER(LEN = 1024)           :: dirname

     REAL(KIND=RK)    :: T_fin, T_deb
     CHARACTER(LEN = 1024)           :: filename
     CHARACTER(LEN = 1024)           :: argname
     INTEGER                     :: li_myRank
     character(len=20) :: ls_stamp_default
     character(len=20)   :: ls_msg
     REAL(KIND=RK), DIMENSION(:), ALLOCATABLE  :: lpr_Y
     REAL(KIND=RK), DIMENSION(1:1,1:2) :: lpr_x
     REAL(KIND=RK), DIMENSION(1:1,1:2) :: lpr_v
     INTEGER                     :: li_i   

  argname = "--parameters"
  CALL JOREK_GET_ARGUMENTS(argname, filename_parameter, ierr)

!  CALL getarg( 1, filename_parameter )
!  PRINT *, "XXXXXXXXXXXXXXXXXXXXXXXXXXXX ", filename_parameter

  CALL INITIALIZE_JOREK_PARAMETERS(filename_parameter)

#ifdef DEBUG_TRACE   
  ! ... Only used in 3D case
  CALL JOREK_Param_GETInt(INT_TRACE_DETAIL_ID,TraceLogDetail,ierr)
  CALL JOREK_Param_GETInt(INT_TRACE_OUTPUT_ID,TraceLogOutput,ierr)
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

  nb_args = iargc()

!  IF (nb_args > 1) THEN
!          PRINT *, "ERROR: Jorek runs without argument. & 
!                  &Input data must be provided in a text file.&
!                  &Please take a look at Jorek Documentation."
!     CALL delegate_stop()
!  END IF

  PRINT *, "RunName : ", Trim(exec)
  
  ! ... Initialize SPM, MPI 
  CALL SPM_INITIALIZE(ierr)
  ! ...

  ! ... Define Model Parameters 
  ! ..................................................

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

  ! ..................................................
  ! ...

  ! -----------------------------------------------------------------
  
  ! output file
  open(unit=li_file_stream_norm, file='output_Var_diag.dat', status='unknown')
  open(unit=li_file_stream_visu, file='output_Var_visu.dat', status='unknown')

  !RUN

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

    

  ! ..................................................
  !END RUN
  
    CALL MODEL_FREE(fem_model)

  close(li_file_stream_norm)
  close(li_file_stream_visu)
  
  CALL SPM_CLEANALL(ierr)
  CALL SPM_FINALIZE(ierr)

#ifdef DEBUG_TRACE   
  CALL closetracelog()
#endif

end subroutine initialize_jorek






    ! ............................................................... 
    SUBROUTINE  RHS_for_Vi(ptr_matrix, ao_BBox2Di, ao_GBox2D)
    IMPLICIT NONE
       CLASS(DEF_MATRIX_2D), POINTER :: ptr_matrix
       TYPE(DEF_BLACKBOX_2D) :: ao_BBox2Di
       TYPE(DEF_GREENBOX_2D)                :: ao_GBox2D
       ! LOCAL
       REAL(KIND=RK) :: f_rhs
       REAL(KIND=RK) :: contribution
       INTEGER       :: ijg
       REAL(KIND=RK) :: wVol
       REAL(KIND=RK) :: Vi_0
       REAL(KIND=RK) :: Vi_R
       REAL(KIND=RK) :: Vi_Z
       REAL(KIND=RK) :: Vi_RR
       REAL(KIND=RK) :: Vi_RZ
       REAL(KIND=RK) :: Vi_ZZ
       REAL(KIND=RK), DIMENSION(:), POINTER :: Rhs_Contribution

       Rhs_Contribution => ao_GBox2D % Rhs_Contribution

       ijg   = ao_BBox2Di % ijg
       wVol  = ao_BBox2Di % wVol(ijg)
       Vi_0  = ao_BBox2Di % B_0(ijg)
       Vi_R  = ao_BBox2Di % B_x1(ijg)
       Vi_Z  = ao_BBox2Di % B_x2(ijg)
       Vi_RR = ao_BBox2Di % B_x1x1(ijg)
       Vi_RZ = ao_BBox2Di % B_x1x2(ijg)
       Vi_ZZ = ao_BBox2Di % B_x2x2(ijg)
  
       ! ... ADD L2 contribution
       f_rhs       = ANALYTICAL_RHS(ao_BBox2Di)

       contribution                = Vi_0 *wVol*f_rhs
       Rhs_Contribution(i_Vu_Rho)  =  contribution 
       ! ...

    END SUBROUTINE RHS_for_Vi
    ! ............................................................... 
  
    ! ............................................................... 
    SUBROUTINE Matrix_for_Vi_Vj(ptr_matrix, ao_BBox2Di, ao_BBox2Dj, ao_GBox2D)
    IMPLICIT NONE
       CLASS(DEF_MATRIX_2D), POINTER :: ptr_matrix
       TYPE(DEF_BLACKBOX_2D) :: ao_BBox2Di
       TYPE(DEF_BLACKBOX_2D) :: ao_BBox2Dj
       TYPE(DEF_GREENBOX_2D)                :: ao_GBox2D
       ! LOCAL
       REAL(KIND=RK) :: contribution
       INTEGER       :: ijg
       REAL(KIND=RK) :: wVol
       REAL(KIND=RK) :: Vi_0
       REAL(KIND=RK) :: Vi_R
       REAL(KIND=RK) :: Vi_Z
       REAL(KIND=RK) :: Vi_RR
       REAL(KIND=RK) :: Vi_RZ
       REAL(KIND=RK) :: Vi_ZZ
       REAL(KIND=RK) :: Vj_0
       REAL(KIND=RK) :: Vj_R
       REAL(KIND=RK) :: Vj_Z
       REAL(KIND=RK) :: Vj_RR
       REAL(KIND=RK) :: Vj_RZ
       REAL(KIND=RK) :: Vj_ZZ
       REAL(KIND=RK), DIMENSION(:,:), POINTER :: Matrix_Contribution

       Matrix_Contribution => ao_GBox2D % Matrix_Contribution

       ijg   = ao_BBox2Di % ijg
       wVol  = ao_BBox2Di % wVol(ijg)
       Vi_0  = ao_BBox2Di % B_0(ijg)    ; Vj_0  = ao_BBox2Dj % B_0(ijg)
       Vi_R  = ao_BBox2Di % B_x1(ijg)   ; Vj_R  = ao_BBox2Dj % B_x1(ijg)
       Vi_Z  = ao_BBox2Di % B_x2(ijg)   ; Vj_Z  = ao_BBox2Dj % B_x2(ijg)
       Vi_RR = ao_BBox2Di % B_x1x1(ijg) ; Vj_RR = ao_BBox2Dj % B_x1x1(ijg)
       Vi_RZ = ao_BBox2Di % B_x1x2(ijg) ; Vj_RZ = ao_BBox2Dj % B_x1x2(ijg)
       Vi_ZZ = ao_BBox2Di % B_x2x2(ijg) ; Vj_ZZ = ao_BBox2Dj % B_x2x2(ijg)

       ! ... ADD STIFFNESS CONTRIBUTION
       contribution = ( Vi_R*Vj_R + Vi_Z*Vj_Z ) * wVol

       Matrix_Contribution(i_Vu_Rho, i_Vu_Rho) =  contribution 
       ! ...

    END SUBROUTINE Matrix_for_Vi_Vj
    ! ............................................................... 


 ! .........................................
SUBROUTINE Assembly_Diagnostics(ao_BBox2D, ao_GBox2D,nstep)
IMPLICIT NONE
TYPE(DEF_BLACKBOX_2D) :: ao_BBox2D
TYPE(DEF_GREENBOX_2D)                :: ao_GBox2D
INTEGER :: mi_COORDINATES_POLOIDAL,ierr, ijg, kg  , nstep
REAL(KIND=RK) ::wVol

  IF(nstep .ge. 0) THEN

       ijg   = ao_BBox2D % ijg
       wVol  = ao_BBox2D % wVol(ijg)
 
       fem_model % opr_diagnostics(1,nstep+1) = fem_model % opr_diagnostics(1,nstep+1) + &
		  & ao_GBox2D%VarN_0(1, ijg) * &
		  & wVol

       fem_model % opr_diagnostics(2,nstep+1) = fem_model % opr_diagnostics(2,nstep+1) + &
		  & ao_GBox2D%VarN_0(1, ijg) * &
		  & ao_GBox2D%VarN_0(1, ijg) * &
		  & wVol

       fem_model % opr_diagnostics(3,nstep+1) = fem_model % opr_diagnostics(3,nstep+1) + &
                  & ao_GBox2D%VarN_x1(1, ijg) * &
		  & ao_GBox2D%VarN_x1(1, ijg) * &
		  & wVol  + &
		  & ao_GBox2D%VarN_x2(1, ijg) * &
		  & ao_GBox2D%VarN_x2(1, ijg) * &
		  & wVol 

       !... Diff norms
       fem_model % opr_diagnostics(4,nstep+1) = fem_model % opr_diagnostics(4,nstep+1) + &
                  & ( ao_GBox2D%VarN_0(1, ijg) - ao_GBox2D%Sol_analytical(ijg, 1, 1) ) * &
		  & ( ao_GBox2D%VarN_0(1, ijg) - ao_GBox2D%Sol_analytical(ijg, 1, 1) ) * &
		  & wVol 

       fem_model % opr_diagnostics(5,nstep+1) = fem_model % opr_diagnostics(5,nstep+1) + &
		  & ( ao_GBox2D%VarN_x1(1, ijg) - ao_GBox2D%Sol_analytical(ijg, 1, 2) ) * &
		  & ( ao_GBox2D%VarN_x1(1, ijg) - ao_GBox2D%Sol_analytical(ijg, 1, 2) ) * &
		  & wVol + &
		  & ( ao_GBox2D%VarN_x2(1, ijg) - ao_GBox2D%Sol_analytical(ijg, 1, 3) ) * &
		  & ( ao_GBox2D%VarN_x2(1, ijg) - ao_GBox2D%Sol_analytical(ijg, 1, 3) ) * &
		  & wVol 

  END IF
  
  RETURN 

END SUBROUTINE Assembly_Diagnostics
! .........................................

SUBROUTINE Plot_Diagnostics(nstep)
  IMPLICIT NONE
  INTEGER ::mi_COORDINATES_POLOIDAL,ierr,nstep
  REAL(KIND=RK) :: lr_dt, lr_time
  
  CALL JOREK_Param_GETInt(INT_COORDINATES_POLOIDAL_ID,mi_COORDINATES_POLOIDAL,ierr)
  CALL JOREK_Param_GETReal(REAL_DT_ID,lr_dt,ierr)

    IF(nstep .ge. 0) THEN
  
     fem_model % opr_diagnostics(2,nstep+1) = SQRT(fem_model % opr_diagnostics(2,nstep+1))
     fem_model % opr_diagnostics(3,nstep+1) = SQRT(fem_model % opr_diagnostics(3,nstep+1))

     fem_model % opr_diagnostics(4,nstep+1) = SQRT(fem_model % opr_diagnostics(4,nstep+1))
     fem_model % opr_diagnostics(5,nstep+1) = SQRT(fem_model % opr_diagnostics(5,nstep+1))
     
     PRINT *, "======      Masse ======"
     PRINT *,'Masse :',fem_model % opr_diagnostics(1,nstep+1)
     PRINT *, "======      NORMS ======"
     PRINT *,'Norm L2 :',fem_model % opr_diagnostics(2,nstep+1)
     PRINT *,'Norm H1 :',fem_model % opr_diagnostics(3,nstep+1)
     PRINT *, "====== DIFF-NORMS ======"
     PRINT *,'Error L2 :',fem_model % opr_diagnostics(4,nstep+1)/fem_model % opr_diagnostics(2,nstep+1)
     PRINT *,'Error H1 :',fem_model % opr_diagnostics(5,nstep+1)/fem_model % opr_diagnostics(3,nstep+1)
     PRINT *, "========================"
     
   
        IF(nstep .eq. 0) THEN
           write(li_file_stream_norm,*) '# ,time, Var_id, masse,  norm_L2, semi_norm_H1, diff_norm_L2, diff_semi_norm_H1'
        END IF
        lr_time=nstep*lr_dt
        write(li_file_stream_norm,*) lr_time, 1, fem_model % opr_diagnostics(1,nstep+1), fem_model % opr_diagnostics(2,nstep+1), &
             fem_model % opr_diagnostics(3,nstep+1), fem_model % opr_diagnostics(4,nstep+1),fem_model % opr_diagnostics(5,nstep+1)
     END IF

  
  RETURN 
END SUBROUTINE Plot_Diagnostics
! .........................................


  !! Test case 1 :
  !!             - mesh: square, collela or square-periodic
  !!             - R0=Z0=1 (lenght square)
  !!             - Cylindrical or Cartesian coordindate
  !!             - Dirichet or periodic boundary condition
  !!             - f(x,y)= 8pi**2 * SIN(2pi*R)*SIN(2pi*Z)

  !! Test case 2 : validate this case

  !! Test case 3 : validate this case
  
  !! Test case 4 :
  !!             - mesh: square-periodic
  !!             - R0=Z0=1 (lenght square)
  !!             - Cylindrical or Cartesian coordindate
  !!             - Periodic boundary condition
  !!             - f(x,y)= 8pi**2 * cos(2pi*R)*cos(2pi*Z)

  
  FUNCTION ANALYTICAL_RHS(ao_BBox2D)
  IMPLICIT NONE
     REAL(KIND=RK) :: ANALYTICAL_RHS 
     TYPE(DEF_BLACKBOX_2D) :: ao_BBox2D
     REAL(KIND=RK) :: lr_R, lr_Z
     ! LOCAL
     REAL(KIND=RK) :: lr_k1
     REAL(KIND=RK) :: lr_k2
     INTEGER       :: li_testcase
     INTEGER       :: ierr
     REAL(KIND=RK) :: lr_R0
     REAL(KIND=RK) :: lr_a
     REAL(KIND=RK) :: lr_acenter
     INTEGER       :: ijg

     ijg = ao_BBox2D%ijg

     lr_R   = ao_BBox2D%Xp_0(1,ijg)
     lr_Z   = ao_BBox2D%Xp_0(2,ijg)

     lr_R0         = mr_R0 
     lr_a          = mr_a
     lr_acenter    = mr_acenter    
    
     lr_k1 = 2.0 * PI * FLOAT(mi_mode_m1)
     lr_k2 = 2.0 * PI * FLOAT(mi_mode_n1)
    
     CALL JOREK_Param_GETInt(INT_TESTCASE_ID, li_testcase, ierr)
     IF(li_testcase .eq. 1) THEN
        ANALYTICAL_RHS = (lr_k1**2 + lr_k2**2) * SIN(lr_k1*lr_R)*SIN(lr_k2*lr_Z)
     ENDIF    
     IF(li_testcase .eq. 2) THEN
        ANALYTICAL_RHS = 4.0 * ( lr_R**2 + lr_Z**2 ) * SIN ( 1.0 - lr_R**2 - lr_Z**2 ) &
		& + 4.0 * COS ( 1.0 - lr_R**2 - lr_Z**2 ) 
     ENDIF
     IF(li_testcase .eq. 3) THEN
        ANALYTICAL_RHS =8*lr_Z**2/lr_a**2 - 4 + 2*(2*lr_R - 2*lr_R0)**2/lr_a**2 &
		&+ 4*(lr_Z**2 + (lr_R - lr_R0)**2)/lr_a**2 &
		& + 4*(lr_Z**2 - lr_acenter**2 + (lr_R - lr_R0)**2)/lr_a**2 
     ENDIF

     ! periodic but not dirichet homogeneous
      IF(li_testcase .eq. 4) THEN
        ANALYTICAL_RHS = (lr_k1**2 + lr_k2**2) * cos(lr_k1*lr_R)*cos(lr_k2*lr_Z)
     ENDIF 
     
     
  
  END FUNCTION ANALYTICAL_RHS
  ! ..................................................
  
  ! ..................................................
  SUBROUTINE ANALYTICAL_MODEL(apr_x, apr_v,apr_info,api_info,n_variable,n_dimension)
  IMPLICIT NONE
     INTEGER                            :: n_variable,n_dimension
     REAL(KIND=8), DIMENSION(n_dimension), INTENT(IN)  :: apr_x
     REAL(KIND=8), DIMENSION(N_variable, n_dimension+1), INTENT(OUT) :: apr_v
     REAL(KIND=8), DIMENSION(10), INTENT(IN)  :: apr_info
     INTEGER, DIMENSION(10), INTENT(IN) :: api_info
     ! LOCAL
     REAL(KIND=RK) :: lr_R
     REAL(KIND=RK) :: lr_Z
     REAL(KIND=RK) :: lr_k1
     REAL(KIND=RK) :: lr_k2
     INTEGER       :: li_testcase
     INTEGER       :: ierr
     REAL(KIND=RK) :: lr_R0
     REAL(KIND=RK) :: lr_a
     REAL(KIND=RK) :: lr_acenter

     lr_R0         = mr_R0 
     lr_a          = mr_a
     lr_acenter    = mr_acenter

     lr_R = apr_x(1)
     lr_Z = apr_x(2)

     lr_k1 = 2.0 * PI * FLOAT(mi_mode_m1)
     lr_k2 = 2.0 * PI * FLOAT(mi_mode_n1)

     CALL JOREK_Param_GETInt(INT_TESTCASE_ID, li_testcase, ierr)
     IF(li_testcase .eq. 1) THEN
        ! ... u
        apr_v(1, 1) = SIN(lr_k1*lr_R)*SIN(lr_k2*lr_Z)
        ! ... u_R
        apr_v(1, 2) = lr_k1*COS(lr_k1*lr_R)*SIN(lr_k2*lr_Z)
        ! ... u_Z
        apr_v(1, 3) = lr_k2*SIN(lr_k1*lr_R)*COS(lr_k2*lr_Z)
        ! ...
     ENDIF    
     IF(li_testcase .eq. 2) THEN
        ! ... u
        apr_v(1, 1) = SIN ( 1.0 - lr_R**2 - lr_Z**2 )
        ! ... u_R
        apr_v(1, 2) = - 2.0 * lr_R * COS ( 1.0 - lr_R**2 - lr_Z**2 )
        ! ... u_Z
        apr_v(1, 3) = - 2.0 * lr_Z * COS ( 1.0 - lr_R**2 - lr_Z**2 )
        ! ...
     ENDIF
     IF(li_testcase .eq. 3) THEN
        ! ... u
        apr_v(1, 1) = (1 - (lr_Z**2 + (lr_R - lr_R0)**2)/lr_a**2)*(lr_Z**2 - lr_acenter**2 + (lr_R - lr_R0)**2)
        ! ... u_R
        apr_v(1, 2) = (1 - (lr_Z**2 + (lr_R - lr_R0)**2)/lr_a**2)*(2*lr_R - 2*lr_R0) &
		& - (2*lr_R - 2*lr_R0)*(lr_Z**2 - lr_acenter**2 + (lr_R - lr_R0)**2)/lr_a**2
        ! ... u_Z
        apr_v(1, 3) = 2*lr_Z*(1 - (lr_Z**2 + (lr_R - lr_R0)**2)/lr_a**2) &
		& - 2*lr_Z*(lr_Z**2 - lr_acenter**2 + (lr_R - lr_R0)**2)/lr_a**2
        ! ...
     ENDIF

     ! periodic but not dirichet homogeneous
     IF(li_testcase .eq. 4) THEN
        ! ... u
        apr_v(1, 1) = cos(lr_k1*lr_R)*cos(lr_k2*lr_Z)
        ! ... u_R
        apr_v(1, 2) = -lr_k1*sin(lr_k1*lr_R)*cos(lr_k2*lr_Z)
        ! ... u_Z
        apr_v(1, 3) = -lr_k2*cos(lr_k1*lr_R)*sin(lr_k2*lr_Z)
        ! ...
     ENDIF    
     

  END SUBROUTINE ANALYTICAL_MODEL
  ! ..................................................
end module sll_jorek
