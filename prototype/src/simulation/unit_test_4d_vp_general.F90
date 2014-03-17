! Sample computation with the following characteristics:
! - vlasov-poisson
! - 4D: x, y, vx, vy (or x1, x2, x3, x4) with arbitrary coordinate 
!   transformation
!   in the x,y variables.
! - parallel

program vlasov_poisson_4d_general
#include "sll_working_precision.h"
  use sll_simulation_4d_vlasov_poisson_general
  use sll_collective
  use sll_constants
  use sll_logical_meshes
  use sll_module_coordinate_transformations_2d
  use sll_common_coordinate_transformations
  use sll_common_array_initializers_module
  implicit none

  character(len=256) :: filename
  character(len=256) :: filename_local
  type(sll_simulation_4d_vp_general)      :: simulation
  type(sll_logical_mesh_2d), pointer      :: mx
  type(sll_logical_mesh_2d), pointer      :: mv
  class(sll_coordinate_transformation_2d_base), pointer :: transformation_x
  sll_real64, dimension(1:5) :: landau_params
  sll_real64, dimension(1:6) :: gaussian_params

  print *, 'Booting parallel environment...'
  call sll_boot_collective() ! Wrap this up somewhere else

  ! In this test, the name of the file to open is provided as a command line
  ! argument.
  call getarg(1, filename)
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
  mx => new_logical_mesh_2d( NPTS1, NPTS2,       & 
       eta1_min=.0_f64, eta1_max=4.0_f64*sll_pi)

  ! logical mesh for velocity coordinates
  mv => new_logical_mesh_2d( NPTS3, NPTS4, &
       eta1_min=-6.0_f64, eta1_max=6.0_f64, &
       eta2_min=-6.0_f64, eta2_max=6.0_f64)
  print *, 'allocated logical meshes'
!  ! logical mesh for space coordinates
!  mx => new_logical_mesh_2d( NPTS1, NPTS2)
!
!  ! logical mesh for velocity coordinates
!  mv => new_logical_mesh_2d( NPTS1, NPTS2, &
!       eta1_min=-6.0_f64, eta1_max=6.0_f64, &
!       eta2_min=-6.0_f64, eta2_max=6.0_f64)
!
  ! coordinate transformation associated with space coordinates
  transformation_x => new_coordinate_transformation_2d_analytic( &
       "analytic_identity_transformation", &
       mx, &
       identity_x1, &
       identity_x2, &
       identity_jac11, &
       identity_jac12, &
       identity_jac21, &
       identity_jac22, &
       (/ 0.0_f64 /) )
 print *, 'allocated coordinate transformation'

!  transformation_x => new_coordinate_transformation_2d_analytic( &
!       "analytic_sinprod_transformation", &
!       mx, &
!       sinprod_x1, &
!       sinprod_x2, &
!       sinprod_jac11, &
!       sinprod_jac12, &
!       sinprod_jac21, &
!       sinprod_jac22 )

  ! define the values of the parameters for the landau initializer

!!$  gaussian_params(1) = 2.0*sll_pi !xc
!!$  gaussian_params(2) = 2.0*sll_pi !yc
!!$  gaussian_params(3) = 0.0        !vxc
!!$  gaussian_params(4) = 0.0        !vyc
!!$  gaussian_params(5) = 1.0        !vxc
!!$  gaussian_params(6) = 0.0        !vyc

  landau_params(1) = 0.0      !eta1_min
  landau_params(2) = mx%eta1_max
  landau_params(3) = 0.0      !eta2_min
  landau_params(4) = mx%eta2_max
  landau_params(5) = 0.05!0.01     !eps

  ! initialize simulation object with the above parameters
  call initialize_vp4d_general( &
       simulation, &
       mx, &
       mv, &
       transformation_x, &
       sll_landau_initializer_4d, &
       landau_params )
  print *, 'initialized simulation object'
!  ! define the values of the parameters for the landau initializer
!  gaussian_params(1) = 3.0*sll_pi !xc
!  gaussian_params(2) = 2.0*sll_pi !yc
!  gaussian_params(3) = 0.0        !vxc
!  gaussian_params(4) = 0.0        !vyc
!
!  ! initialize simulation object with the above parameters
!  call initialize_vp4d_general( &
!       simulation, &
!       mx, &
!       mv, &
!       transformation_x, &
!       sll_gaussian_initializer_4d, &
!       gaussian_params )
  print *, ' f initialized '

  call simulation%run( )
  call sll_delete(simulation)
  print *, 'reached end of vp4d test'
  print *, 'PASSED'

  call sll_halt_collective()


end program vlasov_poisson_4d_general


