PROGRAM Main

  ! ------------------------------ Modules used
  USE TypeDef
  USE MeshGen
  USE FEBasis
  USE tracelog_module
  USE SPM_DEF
  USE MODEL
  USE JOREK_PARAM
  USE JOREK_PARAM_DEF
  USE COORDINATES_DEF
  USE COORDINATES

  ! ------------------------------ 
  IMPLICIT NONE

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
  REAL(KIND=RK)    :: T_fin, T_deb

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
  CHARACTER(LEN = 1024)           :: argname
  INTEGER :: li_n_Gauss_Rp, li_n_Gauss_Zp

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
 
  ! ... Define Model Parameters 
  CALL DEFINE_MODEL( )
  ! ...

  ! ... Define Basis Functions 
  CALL InitBasis(&
               & Mesh2D%oi_n_max_order, &
               & Mesh2D%oi_n_max_order, &
               & Mesh2D%oi_n_max_vtex_per_elmt, &
               & Mesh2D%oi_n_max_vtex_per_elmt, &
               & li_n_Gauss_Rp, &
               & li_n_Gauss_Zp, &
               & Mesh2D % ptr_quad % oi_n_points)
  ! ...

  CALL JOREK_Param_GETInt(INT_TYPEMESH_ID,imesh,ierr)
!  PRINT *, ">>> IMESH ", imesh

  CALL InitGrid(Mesh2D,imesh, &
              & Mesh2D%oi_n_Nodes, &
              & Mesh2D%oi_n_elmts, &
              & Mesh2D%oi_n_max_vtex_per_elmt, &
              & Mesh2D%oi_n_max_order, &
              & Mesh2D%oi_n_max_vtex_per_elmt, &
              & Mesh2D%oi_n_max_nbnet, &
              & n_var_unknown, &
              & Mesh2D%oi_n_max_order, &
              & Mesh2D%oi_n_nzero_Bloc, &
              & Mesh2D % ptr_quad % oi_n_points&
              & )

  ! -----------------------------------------------------------------
  
  ! output file
  open(unit=li_file_stream_norm, file='output_Var_diag.dat', status='unknown')
  open(unit=li_file_stream_visu, file='output_Var_visu.dat', status='unknown')
 
  CALL SPM_INITIALIZE(nmatrices, ierr)

  CALL INITIALIZE_MODEL()

  CALL COORDINATES_INITIALIZE(2, Mesh2D % ptr_quad % oi_n_points, ierr)

  CALL RUN_MODEL()
  
  CALL FREE_MODEL()

  close(li_file_stream_norm)
  close(li_file_stream_visu)
  
  CALL SPM_CLEANALL(ierr)
  CALL SPM_FINALIZE(ierr)

#ifdef DEBUG_TRACE   
  CALL closetracelog()
#endif

END PROGRAM Main
