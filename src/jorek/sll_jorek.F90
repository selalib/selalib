module sll_jorek
#include "sll_working_precision.h"

  use TypeDef
  use MeshGen
  use FEBasis
  use SPM_DEF
  use MODEL
  use JOREK_PARAM
  use JOREK_PARAM_DEF
  use COORDINATES_DEF
  use COORDINATES

  ! ------------------------------ 
  IMPLICIT NONE

contains


subroutine initialize_jorek()

  sll_int32 :: nb_args
  sll_int32 :: is, ie, il, i1, i2, is1, is2
  sll_int32 :: imesh=1, Nrefine
  sll_int32 :: iv, iv_Pol, ipol, Ns_3D
  sll_int32 :: ierr
  sll_int32 :: TraceLogDetail 
  sll_int32 :: TraceLogOutput

  sll_int32, dimension(:,:), pointer     :: Nu
  sll_int32, dimension(:,:), allocatable :: Ok
  sll_int32(kind=SPM_INTS_KIND)          :: nRows !NUMBER OF ROWS
  sll_int32(kind=SPM_INTS_KIND)          :: nCols !NUMBER OF COLUMNS

  real(kind=RK)    :: T_fin, T_deb

  real(kind=RK), dimension(:,:), pointer :: Coor  

  character(len = 1024)         :: exec
  character(len = 1024)         :: rootname
  character(len = *), parameter :: mod_name = "poisson_2d"

  LOGICAL :: ll_stdoutput
  sll_int32 :: li_assembly_proc
  sll_int32 :: li_toroidal_basis
  character filename_parameter*1024
  character(len = 1024)           :: argname
  sll_int32 :: li_n_Gauss_Rp, li_n_Gauss_Zp

  argname = "--parameters"
  call jorek_get_arguments(argname, filename_parameter, ierr)

  call initialize_jorek_parameters(filename_parameter)

  call getarg(0, exec)
  rootname = "output"

  nb_args = iargc()

  print*, "RunName : ", trim(exec)
 
  ! ... Define Model Parameters 
  call define_model( )

  ! ... Define Basis Functions 
  call InitBasis( mesh2d%oi_n_max_order,         &
                  mesh2d%oi_n_max_order,         &
                  mesh2d%oi_n_max_vtex_per_elmt, &
                  mesh2d%oi_n_max_vtex_per_elmt, &
                  li_n_Gauss_Rp,                 &
                  li_n_Gauss_Zp,                 &
                  mesh2d%ptr_quad%oi_n_points)

  call JOREK_Param_GETInt(INT_TYPEMESH_ID,imesh,ierr)

  call InitGrid(mesh2d,                          &
                imesh,                           &
                mesh2d%oi_n_Nodes,               &
                mesh2d%oi_n_elmts,               &
                mesh2d%oi_n_max_vtex_per_elmt,   &
                mesh2d%oi_n_max_order,           &
                mesh2d%oi_n_max_vtex_per_elmt,   &
                mesh2d%oi_n_max_nbnet,           &
                n_var_unknown,                   &
                mesh2d%oi_n_max_order,           &
                mesh2d%oi_n_nzero_Bloc,          &
                mesh2d%ptr_quad%oi_n_points      )

  ! output file
  open(unit=li_file_stream_norm, file='output_Var_diag.dat', status='unknown')
  open(unit=li_file_stream_visu, file='output_Var_visu.dat', status='unknown')
 
  call SPM_INITIALIZE(nmatrices, ierr)

  call INITIALIZE_MODEL()

  call COORDINATES_INITIALIZE(2, mesh2d % ptr_quad % oi_n_points, ierr)

  call RUN_MODEL()
  
  call FREE_MODEL()

  close(li_file_stream_norm)
  close(li_file_stream_visu)
  
  call SPM_CLEANALL(ierr)
  call SPM_FINALIZE(ierr)

end subroutine initialize_jorek

end module sll_jorek
