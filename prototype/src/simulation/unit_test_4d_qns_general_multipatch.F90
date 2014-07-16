! Sample computation with the following characteristics:
! - vlasov-poisson
! - 4D: x, y, vx, vy (or x1, x2, x3, x4) with an arbitrary coordinate 
!   transformation in the x,y variables.
! - parallel
! - The coordinate transformation is defined by patches.

program qns_4d_general_multipatch
#include "sll_working_precision.h"
  use sll_coordinate_transformation_multipatch_module
  use sll_simulation_4d_qns_general_multipatch_module, only: &
     sll_simulation_4d_qns_general_multipatch, &
     initialize_4d_qns_gen_mp, &
     run_4d_qns_general_mp
  use sll_collective
!  use sll_constants
  use sll_logical_meshes
  use sll_module_coordinate_transformations_2d
  use sll_common_coordinate_transformations
!  use sll_module_coordinate_transformations_2d_nurbs
  use sll_common_array_initializers_module
  use sll_module_poisson_2d_elliptic_solver, &
     only: es_gauss_legendre
  use sll_coordinate_transformation_multipatch_module, only: &
     sll_coordinate_transformation_multipatch_2d
!  use sll_module_scalar_field_2d_multipatch
  implicit none


  character(len=256) :: filename
  character(len=256) :: filename_local
  type(sll_simulation_4d_qns_general_multipatch)     :: simulation
  type(sll_logical_mesh_2d), pointer      :: mv ! delete
  type(sll_coordinate_transformation_multipatch_2d) :: mp

  sll_real64, dimension(1:8) :: gaussian_beam_params
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
  elec_field_ext_params(:) =  (/64.0_f64,64.0_f64/)!(/36.0_f64,36.0_f64/)

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
#define NPTS3 64
#define NPTS4 64
#define SPL_DEG_VX 3 
#define SPL_DEG_VY 3


  ! logical mesh for velocity coordinates
  mv => new_logical_mesh_2d( NPTS3, NPTS4, &
       eta1_min=-9.0_f64, eta1_max=9.0_f64, &
       eta2_min=-9.0_f64, eta2_max=9.0_f64)

  ! ---------------------------------------------------------------------
  ! coordinate transformation associated with space coordinates
  ! ---------------------------------------------------------------------

  call mp%read_from_file("identity_mp2")

  ! ---------------------------------------------------------------------
  ! define the values of the parameters for the gaussian beam initializer
  ! ---------------------------------------------------------------------

  ! sll_gaussian_beam_initializer_4d parameters
  !TODO : fix parameteres

  gaussian_beam_params(1) = 1.0_f64!3.0_f64/2.0_f64!vth
  gaussian_beam_params(2) = 0.05_f64!1.0_f64!xth
  gaussian_beam_params(3) = 0.0_f64!vxc
  gaussian_beam_params(4) = 0.0_f64!vyc
  gaussian_beam_params(5) = 0.5_f64!xc
  gaussian_beam_params(6) = 0.5_f64!yc
  gaussian_beam_params(7) = sll_pi*15*8!sll_pi*15*18!n0
  gaussian_beam_params(8) = 1.0_f64 !radius
  

  ! ---------------------------------------------------------------------
  ! initialize simulation object with the above parameters
  ! ---------------------------------------------------------------------
  call simulation%initialize( &
       mv, &
       mp, &
       sll_gaussian_beam_initializer_4d, &
       gaussian_beam_params, &
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
       SPL_DEG_VX, & 
       SPL_DEG_VY, & 
       ES_GAUSS_LEGENDRE,&
       ES_GAUSS_LEGENDRE,&
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
       elec_field_ext_params,&
       100)


  print *, ' f initialized '

  call simulation%run( )
 ! call delete(simulation)
  print *, 'reached end of vp4d test'
  print *, 'PASSED'

  call sll_halt_collective()

  !call delete(simulation)
  !call delete(transformation_x)

end program qns_4d_general_multipatch

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



  
