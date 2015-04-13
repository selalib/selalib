module sll_jorek
#include "sll_working_precision.h"

  use typedef
  use meshgen
  use febasis
  use spm_def
  use model
  use jorek_param
  use jorek_param_def
  use coordinates_def
  use coordinates

  ! ------------------------------ 
  IMPLICIT NONE

contains


subroutine initialize_jorek()

  sll_int32 :: nb_args
  sll_int32 :: imesh=1
  sll_int32 :: ierr

  character(len = 1024) :: exec
  character(len = 1024) :: rootname
  character(len = 1024) :: flag
  character(len = *), parameter :: mod_name = "poisson_2d"
  character filename_parameter*1024
  sll_int32 :: li_n_Gauss_Rp, li_n_Gauss_Zp

  flag = "--parameters"
  call jorek_get_arguments(flag, filename_parameter, ierr)

  call initialize_jorek_parameters(filename_parameter)

  call getarg(0, exec)
  rootname = "output"

  nb_args = iargc()

  print*, "RunName : ", trim(exec)
 
  ! ... Define Model Parameters 
  call define_model( )

  ! ... Define Basis Functions 
  call initbasis( mesh2d%oi_n_max_order,         &
                  mesh2d%oi_n_max_order,         &
                  mesh2d%oi_n_max_vtex_per_elmt, &
                  mesh2d%oi_n_max_vtex_per_elmt, &
                  li_n_Gauss_Rp,                 &
                  li_n_Gauss_Zp,                 &
                  mesh2d%ptr_quad%oi_n_points)

  call jorek_paral_getint(int_typemesh_id,imesh,ierr)

  call initgrid(mesh2d,                          &
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

end module sll_jorek
