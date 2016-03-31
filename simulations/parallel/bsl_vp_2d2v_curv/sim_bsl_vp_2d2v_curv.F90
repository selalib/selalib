! Sample computation with the following characteristics:
! - vlasov-poisson
! - 4D: x, y, vx, vy (or x1, x2, x3, x4) with arbitrary coordinate 
!   transformation
!   in the x,y variables.
! - parallel

program sim_bsl_vp_2d2v_curv
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use sll_m_cartesian_meshes, only: &
    sll_f_new_cartesian_mesh_2d, &
    sll_t_cartesian_mesh_2d, &
    sll_o_delete

  use sll_m_collective, only: &
    sll_s_boot_collective, &
    sll_s_halt_collective

  use sll_m_common_array_initializers, only: &
    sll_f_landau_initializer_4d

  use sll_m_common_coordinate_transformations, only: &
    sll_f_identity_jac11, &
    sll_f_identity_jac12, &
    sll_f_identity_jac21, &
    sll_f_identity_jac22, &
    sll_f_identity_x1, &
    sll_f_identity_x2

  use sll_m_constants, only: &
    sll_p_pi

  use sll_m_coordinate_transformation_2d_base, only: &
    sll_c_coordinate_transformation_2d_base

  use sll_m_coordinate_transformations_2d, only: &
    sll_f_new_coordinate_transformation_2d_analytic

  use sll_m_sim_bsl_vp_2d2v_curv, only: &
    sll_s_initialize_vp4d_general, &
    sll_t_simulation_4d_vp_general, &
    sll_o_delete

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  character(len=256) :: filename
  character(len=256) :: filename_local
  type(sll_t_simulation_4d_vp_general)      :: simulation
  type(sll_t_cartesian_mesh_2d), pointer      :: mx
  type(sll_t_cartesian_mesh_2d), pointer      :: mv
  class(sll_c_coordinate_transformation_2d_base), pointer :: transformation_x
  sll_real64, dimension(1:5) :: landau_params

  print *, 'Booting parallel environment...'
  call sll_s_boot_collective() ! Wrap this up somewhere else

  ! In this test, the name of the file to open is provided as a command line
  ! argument.
  call get_command_argument(1, filename)
  filename_local = trim(filename)

  ! To initialize the simulation type, there should be two options. One is to
  ! initialize from a file:
  print *, 'executing unit_test_4d_vp_general'
  
  call simulation%init_from_file(filename_local)
  print *, 'finished initialization'
  
  ! The second is to initialize 'manually' with a routine whose parameters
  ! allow to configure the different types of objects in the simulation. For
  ! instance, the type of coordinate mapping. Here we use both methods while
  ! we develop and sort out the interfaces.
  ! Eventually, when using the module, one should only need to use one 
  ! way to initialize the simulation object, in development we are using them
  ! both...

! hardwired, this should be consistent with whatever is read from a file
#define NPTS1 32
#define NPTS2 32
#define NPTS3 32
#define NPTS4 32

  ! logical mesh for space coordinates
  mx => sll_f_new_cartesian_mesh_2d( NPTS1, NPTS2,       & 
       eta1_min=.0_f64, eta1_max=4.0_f64*sll_p_pi)

  ! logical mesh for velocity coordinates
  mv => sll_f_new_cartesian_mesh_2d( NPTS3, NPTS4, &
       eta1_min=-6.0_f64, eta1_max=6.0_f64, &
       eta2_min=-6.0_f64, eta2_max=6.0_f64)
  print *, 'allocated logical meshes'
!  ! logical mesh for space coordinates
!  mx => sll_f_new_cartesian_mesh_2d( NPTS1, NPTS2)
!
!  ! logical mesh for velocity coordinates
!  mv => sll_f_new_cartesian_mesh_2d( NPTS1, NPTS2, &
!       eta1_min=-6.0_f64, eta1_max=6.0_f64, &
!       eta2_min=-6.0_f64, eta2_max=6.0_f64)
!
  ! coordinate transformation associated with space coordinates
  transformation_x => sll_f_new_coordinate_transformation_2d_analytic( &
       "analytic_identity_transformation", &
       mx, &
       sll_f_identity_x1, &
       sll_f_identity_x2, &
       sll_f_identity_jac11, &
       sll_f_identity_jac12, &
       sll_f_identity_jac21, &
       sll_f_identity_jac22, &
       (/ 0.0_f64 /) )
 print *, 'allocated coordinate transformation'

!  transformation_x => sll_f_new_coordinate_transformation_2d_analytic( &
!       "analytic_sinprod_transformation", &
!       mx, &
!       sll_f_sinprod_x1, &
!       sll_f_sinprod_x2, &
!       sll_f_sinprod_jac11, &
!       sll_f_sinprod_jac12, &
!       sll_f_sinprod_jac21, &
!       sll_f_sinprod_jac22 )

  ! define the values of the parameters for the landau initializer

!!$  gaussian_params(1) = 2.0*sll_p_pi !xc
!!$  gaussian_params(2) = 2.0*sll_p_pi !yc
!!$  gaussian_params(3) = 0.0        !vxc
!!$  gaussian_params(4) = 0.0        !vyc
!!$  gaussian_params(5) = 1.0        !vxc
!!$  gaussian_params(6) = 0.0        !vyc

  landau_params(1) = 0.0_f64      !eta1_min
  landau_params(2) = mx%eta1_max
  landau_params(3) = 0.0_f64      !eta2_min
  landau_params(4) = mx%eta2_max
  landau_params(5) = 0.05_f64!0.01     !eps

  ! initialize simulation object with the above parameters
  call sll_s_initialize_vp4d_general( &
       simulation, &
       mx, &
       mv, &
       transformation_x, &
       sll_f_landau_initializer_4d, &
       landau_params )
  print *, 'initialized simulation object'
!  ! define the values of the parameters for the landau initializer
!  gaussian_params(1) = 3.0*sll_p_pi !xc
!  gaussian_params(2) = 2.0*sll_p_pi !yc
!  gaussian_params(3) = 0.0        !vxc
!  gaussian_params(4) = 0.0        !vyc
!
!  ! initialize simulation object with the above parameters
!  call sll_s_initialize_vp4d_general( &
!       simulation, &
!       mx, &
!       mv, &
!       transformation_x, &
!       sll_gaussian_initializer_4d, &
!       gaussian_params )
  print *, ' f initialized '

  call simulation%run( )
  call sll_o_delete(simulation)
  print *, 'reached end of vp4d test'
  print *, 'PASSED'

  call sll_s_halt_collective()


end program sim_bsl_vp_2d2v_curv


