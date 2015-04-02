module sll_jorek
#include "sll_working_precision.h"
#include "sll_jorek.h"

implicit none

  TYPE(DEF_MESH_2D)              :: Mesh2D
  TYPE(DEF_GREENBOX_2D)   :: GBox2D
  TYPE(DEF_QUADRATURE_SQUARE), TARGET    :: Quad
  TYPE(DEF_LINEAR_SOLVER)   :: Solver

  INTEGER :: n_var_sys 
  INTEGER :: n_var_unknown

  REAL(KIND=RK), DIMENSION(:), ALLOCATABLE  :: Global_Rhs 
  REAL(KIND=RK), DIMENSION(:), ALLOCATABLE  :: Global_Unknown
  REAL(KIND=RK), DIMENSION(:,:), ALLOCATABLE       :: Var
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine initialize_jorek()

  sll_int32 :: nb_args
  sll_int32 :: imesh=1
  sll_int32 :: ierr

  character(len=1024) :: exec
  character(len=1024) :: rootname
  character(len=1024) :: flag
  character(len=1024) :: filename_parameter*1024

  character(len = *), parameter :: mod_name = "poisson_2d"
  sll_int32 :: li_n_gauss_rp, li_n_gauss_zp

  flag = "--parameters"
  call jorek_get_arguments(flag, filename_parameter, ierr)

  call initialize_jorek_parameters(filename_parameter)

  call getarg(0, exec)
  rootname = "output"

  nb_args = iargc()

  print*, "runname : ", trim(exec)
 
  ! ... define model parameters 
  call define_model( )

  ! ... define basis functions 
  call initbasis( mesh2d%oi_n_max_order,         &
                  mesh2d%oi_n_max_order,         &
                  mesh2d%oi_n_max_vtex_per_elmt, &
                  mesh2d%oi_n_max_vtex_per_elmt, &
                  li_n_gauss_rp,                 &
                  li_n_gauss_zp,                 &
                  mesh2d%ptr_quad%oi_n_points)

  call jorek_param_getint(int_typemesh_id,imesh,ierr)

  call initgrid(mesh2d,                          &
                imesh,                           &
                mesh2d%oi_n_nodes,               &
                mesh2d%oi_n_elmts,               &
                mesh2d%oi_n_max_vtex_per_elmt,   &
                mesh2d%oi_n_max_order,           &
                mesh2d%oi_n_max_vtex_per_elmt,   &
                mesh2d%oi_n_max_nbnet,           &
                n_var_unknown,                   &
                mesh2d%oi_n_max_order,           &
                mesh2d%oi_n_nzero_bloc,          &
                mesh2d%ptr_quad%oi_n_points      )

  call spm_initialize(nmatrices, ierr)

  call initialize_model()

  call coordinates_initialize(2, mesh2d%ptr_quad%oi_n_points, ierr)

  call run_model()
  
end subroutine initialize_jorek

subroutine delete_jorek()
  sll_int32 :: ierr

  call free_model()
  call spm_cleanall(ierr)
  call spm_finalize(ierr)

end subroutine delete_jorek

  SUBROUTINE DEFINE_MODEL( )
  IMPLICIT NONE
    INTEGER, PARAMETER          :: N_DIM = 2
    INTEGER :: ierr

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

    CALL CREATE_QUADRATURE(Quad, N_DIM, 0, 3, 3)

    Mesh2D % ptr_quad => Quad
  END SUBROUTINE DEFINE_MODEL

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

  SUBROUTINE FREE_MODEL( )
  IMPLICIT NONE
    DEALLOCATE(Global_Rhs)
    DEALLOCATE(Global_Unknown)
    DEALLOCATE(Var)
    DEALLOCATE(Diagnostics)

    CALL FREE_GREENBOX(GBox2D)
  END SUBROUTINE FREE_MODEL

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

  END SUBROUTINE RUN_MODEL


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

SUBROUTINE Assembly_Diagnostics(ao_BBox2D, ao_GBox2D,nstep)
IMPLICIT NONE
TYPE(DEF_BLACKBOX_2D) :: ao_BBox2D
TYPE(DEF_GREENBOX_2D)                :: ao_GBox2D
INTEGER :: mi_COORDINATES_POLOIDAL,ierr, ijg, kg  , nstep
REAL(KIND=RK) ::wVol

  IF(nstep .ge. 0) THEN

       ijg   = ao_BBox2D % ijg
       wVol  = ao_BBox2D % wVol(ijg)
 
       Diagnostics(1,nstep+1) = Diagnostics(1,nstep+1) + &
		  & ao_GBox2D%VarN_0(1, ijg) * &
		  & wVol

       Diagnostics(2,nstep+1) = Diagnostics(2,nstep+1) + &
		  & ao_GBox2D%VarN_0(1, ijg) * &
		  & ao_GBox2D%VarN_0(1, ijg) * &
		  & wVol

       Diagnostics(3,nstep+1) = Diagnostics(3,nstep+1) + &
                  & ao_GBox2D%VarN_x1(1, ijg) * &
		  & ao_GBox2D%VarN_x1(1, ijg) * &
		  & wVol  + &
		  & ao_GBox2D%VarN_x2(1, ijg) * &
		  & ao_GBox2D%VarN_x2(1, ijg) * &
		  & wVol 

       !... Diff norms
       Diagnostics(4,nstep+1) = Diagnostics(4,nstep+1) + &
                  & ( ao_GBox2D%VarN_0(1, ijg) - ao_GBox2D%Sol_analytical(ijg, 1, 1) ) * &
		  & ( ao_GBox2D%VarN_0(1, ijg) - ao_GBox2D%Sol_analytical(ijg, 1, 1) ) * &
		  & wVol 

       Diagnostics(5,nstep+1) = Diagnostics(5,nstep+1) + &
		  & ( ao_GBox2D%VarN_x1(1, ijg) - ao_GBox2D%Sol_analytical(ijg, 1, 2) ) * &
		  & ( ao_GBox2D%VarN_x1(1, ijg) - ao_GBox2D%Sol_analytical(ijg, 1, 2) ) * &
		  & wVol + &
		  & ( ao_GBox2D%VarN_x2(1, ijg) - ao_GBox2D%Sol_analytical(ijg, 1, 3) ) * &
		  & ( ao_GBox2D%VarN_x2(1, ijg) - ao_GBox2D%Sol_analytical(ijg, 1, 3) ) * &
		  & wVol 

  END IF
  
  RETURN 

END SUBROUTINE Assembly_Diagnostics

SUBROUTINE Plot_Diagnostics(nstep)
  IMPLICIT NONE
  INTEGER ::mi_COORDINATES_POLOIDAL,ierr,nstep
  REAL(KIND=RK) :: lr_dt, lr_time
  
  CALL JOREK_Param_GETInt(INT_COORDINATES_POLOIDAL_ID,mi_COORDINATES_POLOIDAL,ierr)
  CALL JOREK_Param_GETReal(REAL_DT_ID,lr_dt,ierr)

    IF(nstep .ge. 0) THEN
  
     Diagnostics(2,nstep+1) = SQRT(Diagnostics(2,nstep+1))
     Diagnostics(3,nstep+1) = SQRT(Diagnostics(3,nstep+1))

     Diagnostics(4,nstep+1) = SQRT(Diagnostics(4,nstep+1))
     Diagnostics(5,nstep+1) = SQRT(Diagnostics(5,nstep+1))
     
     PRINT *, "======      Masse ======"
     PRINT *,'Masse :',Diagnostics(1,nstep+1)
     PRINT *, "======      NORMS ======"
     PRINT *,'Norm L2 :',Diagnostics(2,nstep+1)
     PRINT *,'Norm H1 :',Diagnostics(3,nstep+1)
     PRINT *, "====== DIFF-NORMS ======"
     PRINT *,'Error L2 :',Diagnostics(4,nstep+1)/Diagnostics(2,nstep+1)
     PRINT *,'Error H1 :',Diagnostics(5,nstep+1)/Diagnostics(3,nstep+1)
     PRINT *, "========================"
     
   
        IF(nstep .eq. 0) THEN
           write(li_file_stream_norm,*) '# ,time, Var_id, masse,  norm_L2, semi_norm_H1, diff_norm_L2, diff_semi_norm_H1'
        END IF
        lr_time=nstep*lr_dt
        write(li_file_stream_norm,*) lr_time, 1, Diagnostics(1,nstep+1), Diagnostics(2,nstep+1), &
             Diagnostics(3,nstep+1), Diagnostics(4,nstep+1),Diagnostics(5,nstep+1)
     END IF

  
  RETURN 
END SUBROUTINE Plot_Diagnostics

    SUBROUTINE  RHS_for_Vi(ao_BBox2Di, ao_GBox2D)
    IMPLICIT NONE
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
  
    SUBROUTINE Matrix_for_Vi_Vj(ao_BBox2Di, ao_BBox2Dj, ao_GBox2D)
    IMPLICIT NONE
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

    END SUBROUTINE Matrix_for_Vi_Vj

end module sll_jorek
