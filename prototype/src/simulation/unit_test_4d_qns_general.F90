! Sample computation with the following characteristics:
! - vlasov-poisson
! - 4D: x, y, vx, vy (or x1, x2, x3, x4) with arbitrary coordinate 
!   transformation
!   in the x,y variables.
! - parallel

program qns_4d_general
#include "sll_working_precision.h"
  use sll_simulation_4d_qns_general_module
  use sll_collective
  use sll_constants
  use sll_logical_meshes
  use sll_module_coordinate_transformations_2d
  use sll_common_coordinate_transformations
  use sll_module_coordinate_transformations_2d_nurbs
  use sll_common_array_initializers_module
  implicit none

  character(len=256) :: filename
  character(len=256) :: filename_local
  type(sll_simulation_4d_qns_general)      :: simulation
  type(sll_logical_mesh_2d), pointer      :: mx
  type(sll_logical_mesh_2d), pointer      :: mv
  class(sll_coordinate_transformation_2d_base), pointer :: transformation_x
  !class(sll_coordinate_transformation_2d_nurbs), pointer :: transformation_x
  sll_real64, dimension(1:8) :: landau_params
  sll_real64, dimension(1:6) :: gaussian_params
  sll_real64, dimension(1:2) :: elec_field_ext_params
  sll_real64, external :: func_zero, func_one, func_minus_one,func_epsi
  sll_real64, dimension(1) :: f_zero_params
  sll_real64, dimension(1) :: f_one_params
  sll_real64, dimension(1) :: f_minus_one_params
  sll_real64, dimension(1) :: f_epsi_params
  sll_real64, external :: electric_field_ext_1
  sll_real64, external :: electric_field_ext_2


  print *, 'Booting parallel environment...'
  call sll_boot_collective() ! Wrap this up somewhere else

  ! In this test, the name of the file to open is provided as a command line
  ! argument.
  call getarg(1, filename)
  filename_local = trim(filename)
  
  f_zero_params(:) = (/0.0_f64/)
  f_one_params(:) = (/0.0_f64/)
  f_minus_one_params(:) = (/0.0_f64/)
  f_epsi_params(:) = (/0.0_f64/)
  elec_field_ext_params(:) = (/1.25_f64,1.25_f64/)

  ! To initialize the simulation type, there should be two options. One is to
  ! initialize from a file:


  
  call simulation%init_from_file(filename_local)
  
  ! The second is to initialize 'manually' with a routine whose parameters
  ! allow to configure the different types of objects in the simulation. For
  ! instance, the type of coordinate mapping. Here we use both methods while
  ! we develop and sort out the interfaces.
  ! Eventually, when using the module, one should only need to use one 
  ! way to initialize the simulation object, in development we are using them
  ! both...

! hardwired, this should be consistent with whatever is read from a file
#define NPTS1 64
#define NPTS2 64
#define NPTS3 64
#define NPTS4 64
#define SPL_DEG_ETA1 3 
#define SPL_DEG_ETA2 3
#define SPL_DEG_VX 3 
#define SPL_DEG_VY 3

  ! ---------------------------------------------------------------------
  ! logical mesh for space coordinates
  ! NPTS1 and NPTS2 is the number of cells
  ! ---------------------------------------------------------------------


  ! mesh for the test case 
  !  sll_periodic_periodic_gaussian2009_initializer_4d   
  !  sll_periodic_periodic_gaussian2002_initializer_4d
  ! sll_landau_initializer_4d

!!$  mx => new_logical_mesh_2d( NPTS1, NPTS2,       & 
!!$       eta1_min=0.0_f64, eta1_max= 4.0_f64*sll_pi, &
!!$       eta2_min=0.0_f64, eta2_max= 4.0_f64*sll_pi )
!!$  
!!$  ! logical mesh for velocity coordinates
!!$  mv => new_logical_mesh_2d( NPTS3, NPTS4, &
!!$       eta1_min=-6.0_f64, eta1_max=6.0_f64, &
!!$       eta2_min=-6.0_f64, eta2_max=6.0_f64)


  ! mesh for the test case 
  ! sll_gaussian_beam_initializer_4d
  
  mx => new_logical_mesh_2d( NPTS1, NPTS2,  & 
       eta1_min=-9.0_f64, eta1_max= 9.0_f64, &
       eta2_min=-9.0_f64, eta2_max= 9.0_f64 )
  
  ! logical mesh for velocity coordinates
  mv => new_logical_mesh_2d( NPTS3, NPTS4, &
       eta1_min=-9.0_f64, eta1_max=9.0_f64, &
       eta2_min=-9.0_f64, eta2_max=9.0_f64)

  ! ---------------------------------------------------------------------
  ! coordinate transformation associated with space coordinates
  ! ---------------------------------------------------------------------

  ! identity transformation

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

  ! colella transformation
  
!!$  transformation_x => new_coordinate_transformation_2d_analytic( &
!!$       "analytic_colela_transformation", &
!!$       mx, &
!!$       sinprod_x1, &
!!$       sinprod_x2, &
!!$       sinprod_jac11, &
!!$       sinprod_jac12, &
!!$       sinprod_jac21, &
!!$       sinprod_jac22, &
!!$       (/ 0.1_f64,0.1_f64,4.0_f64*sll_pi,4.0_f64*sll_pi /) )


!!$   transformation_x => new_nurbs_2d_transformation_from_file("../src/coordinate_transformations/n31x31p3x3_patch0.nml")
!!$   transformation_x%mesh => mx
!!$   print*, 'transformation ok'
  ! ---------------------------------------------------------------------
  ! define the values of the parameters for the landau initializer
  ! ---------------------------------------------------------------------

!!$  ! sll_landau_initializer_4d
!!$  landau_params(1) = mx%eta1_min      !eta1_min
!!$  landau_params(2) = mx%eta1_max
!!$  landau_params(3) = mx%eta2_min      !eta2_min
!!$  landau_params(4) = mx%eta2_max
!!$  landau_params(5) = 0.05     !eps
!!$

!!$  ! 2002 et 2009 sll_periodic_periodic_gaussian2009_initializer_4d   
!!$  !  sll_periodic_periodic_gaussian2002_initializer_4d
!!$  landau_params(1) = mx%eta1_min      !eta1_min
!!$  landau_params(2) = mx%eta1_max
!!$  landau_params(3) = mx%eta2_min      !eta2_min
!!$  landau_params(4) = mx%eta2_max
!!$  landau_params(5) = 0.0_f64     !eps
!!$  landau_params(6) = 0.0_f64     !eps
!!$  landau_params(7) = 0.05_f64     !eps
!!$  landau_params(8) = 1._f64     !eps


  ! sll_gaussian_beam_initializer_4d parameters

  landau_params(1) = 1.0_f64!vth
  landau_params(2) = 1.0_f64!xth
  landau_params(3) = 1.0_f64!sigma_x
  landau_params(4) = 1.0_f64!sigma_v
  landau_params(5) = 0.0_f64!vxc
  landau_params(6) = 0.0_f64!vyc
  landau_params(7) = 0.0_f64!xc
  landau_params(8) = 0.0_f64!yc
  landau_params(9) = 1.0_f64!n0
  

  ! ---------------------------------------------------------------------
  ! initialize simulation object with the above parameters
  ! ---------------------------------------------------------------------
  call initialize_4d_qns_general( &
       simulation, &
       mx, &
       mv, &
       transformation_x, &
       sll_gaussian_beam_initializer_4d, &
       landau_params, &
       func_one,  &  ! a11
       f_one_params, &
       func_zero, &   !a12
       f_zero_params, &
       func_zero, &   !a21
       f_zero_params, &
       func_one,  &   !a22
       f_one_params, &
       func_zero, &   !b1
       f_zero_params, &
       func_zero, &   !der1 b1
       func_zero, &   !der2 b1
       func_zero, &   ! b2
       f_zero_params, &
       func_zero, &   !der1 b2
       func_zero, &   !der2 b2
       func_zero, &   ! c
       f_zero_params, &
       SPL_DEG_ETA1, & 
       SPL_DEG_ETA2, & 
       SPL_DEG_VX, & 
       SPL_DEG_VY, & 
       SLL_DIRICHLET,&!SLL_PERIODIC, &
       SLL_DIRICHLET,&!SLL_PERIODIC, &
       SLL_DIRICHLET,&!SLL_PERIODIC, &
       SLL_DIRICHLET,&!SLL_PERIODIC, &
       SLL_DIRICHLET,&!SLL_PERIODIC, &
       SLL_DIRICHLET,&!SLL_PERIODIC, &
       SLL_DIRICHLET,&!SLL_PERIODIC, &
       SLL_DIRICHLET,&!SLL_PERIODIC, &
       electric_field_ext_1,&
       electric_field_ext_2,&
       elec_field_ext_params)


  print *, ' f initialized '

  call simulation%run( )
 ! call delete(simulation)
  print *, 'reached end of vp4d test'
  print *, 'PASSED'

  call sll_halt_collective()

  !call delete(simulation)
  !call delete(transformation_x)

end program qns_4d_general

! External functions used as parameters in the above unit test:


function func_one( eta1, eta2, params ) result(res)
  real(8), intent(in) :: eta1
  real(8), intent(in) :: eta2
  real(8), dimension(:), intent(in) :: params
  real(8) :: res
  res = 1.0_8
end function func_one

function func_minus_one( eta1, eta2, params ) result(res)
  real(8), intent(in) :: eta1
  real(8), intent(in) :: eta2
  real(8), dimension(:), intent(in) :: params
  real(8) :: res
  res = -1.0_8
end function func_minus_one

function func_zero( eta1, eta2, params ) result(res)
  real(8), intent(in) :: eta1
  real(8), intent(in) :: eta2
  real(8), dimension(:), intent(in) :: params
  real(8) :: res
  res = 0.0_8
end function func_zero


function func_epsi( eta1, eta2, params ) result(res)
  real(8), intent(in) :: eta1
  real(8), intent(in) :: eta2
  real(8), dimension(:), intent(in) :: params
  real(8) :: res
  res = 0.00001_8
end function func_epsi

function electric_field_ext_1(x,y,params) result(res)
  real(8), intent(in) :: x
  real(8), intent(in) :: y
  real(8), dimension(:), intent(in) :: params
  real(8) :: res
  res = params(1)*x
end function electric_field_ext_1


function electric_field_ext_2(x,y,params) result(res)
  real(8), intent(in) :: x
  real(8), intent(in) :: y
  real(8), dimension(:), intent(in) :: params
  real(8) :: res
  res = params(2)*y
end function electric_field_ext_2



  
