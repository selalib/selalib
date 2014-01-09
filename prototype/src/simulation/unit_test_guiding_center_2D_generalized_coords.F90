! Sample computation with the following characteristics:
! - guiding center-poisson
! - 2D: x, y with arbitrary coordinate 
!   transformation
!   in the x,y variables.


program gc_2d_general
#include "sll_working_precision.h"
  use sll_simulation_2d_guiding_center_generalized_coords_module
  use sll_constants
  use sll_logical_meshes
  use sll_module_coordinate_transformations_2d
  use sll_common_coordinate_transformations
  use sll_common_array_initializers_module
  use sll_module_scalar_field_2d_alternative
  use sll_module_scalar_field_2d_base
  use sll_arbitrary_degree_spline_interpolator_2d_module
  use sll_module_interpolators_2d_base
  implicit none

  character(len=256) :: filename
  character(len=256) :: filename_local
  type(sll_simulation_2d_guiding_center_generalized)      :: simulation
  type(sll_logical_mesh_2d), pointer      :: M
  class(sll_coordinate_transformation_2d_base), pointer :: transformation
  
  sll_real64, external ::     func_zero, func_one, func_minus_one
  
  sll_real64, external ::     func_khp1,func_khp2
  sll_real64, external ::     sinprod2_x1, sinprod2_x2
  sll_real64, external ::     sinprod2_jac11, sinprod2_jac12
  sll_real64, external ::     sinprod2_jac21, sinprod2_jac22 
  
  sll_real64 :: eta1_min
  sll_real64 :: eta2_min  
  sll_real64 :: eta1_max
  sll_real64 :: eta2_max
  sll_real64, dimension(1:2) :: landau_params
  sll_real64, dimension(1)   :: params_field 
  sll_real64, dimension(4)   :: params_mesh
  
  
  print *, 'Start...'
  

  ! In this test, the name of the file to open is provided as a command line
  ! argument.
  
  !call getarg(1, filename)
  !filename_local = trim(filename)

  ! To initialize the simulation type, there should be two options. One is to
  ! initialize from a file:
  
  
  
  ! hardwired, this should be consistent with whatever is read from a file
#define NCELL1 64
#define NCELL2 64
#define SPL_DEG1 2
#define SPL_DEG2 2

    landau_params(1)=0.5_f64 ! mode
    landau_params(2)=0.015_f64
    params_field = (/1.0_f64/)
     
    eta1_min = 0._f64
    eta1_max = 2._f64  * sll_pi/landau_params(1)
    eta2_min = 0._f64 
    eta2_max = 2._f64  * sll_pi

    !  In collela  mesh params_mesh =( alpha1, alpha2, L1, L2 ) such that :
    !  x1= eta1 + alpha1*sin(2*pi*eta1/L1)*sin(2*pi*eta2/L2)
    params_mesh = (/ 0.1_f64, 0.1_f64, 1.0_f64, 1.0_f64/)
 
  
  !call simulation%init_from_file(filename_local)
    simulation%dt = 0.1_f64
    simulation%num_iterations = 5 !600
    simulation%carac_case = 3
    simulation%time_scheme = 2
    simulation%visu_step = 100
    simulation%nc_x1 = NCELL1
    simulation%nc_x2 = NCELL2
    
  
  ! The second is to initialize 'manually' with a routine whose parameters
  ! allow to configure the different types of objects in the simulation. For
  ! instance, the type of coordinate mapping. Here we use both methods while
  ! we develop and sort out the interfaces.
  ! Eventually, when using the module, one should only need to use one 
  ! way to initialize the simulation object, in development we are using them
  ! both...

print*,'pass -1'
  ! logical mesh for space coordinates
  M => new_logical_mesh_2d( NCELL1, NCELL2,       & 
       eta1_min, eta1_max,eta2_min, eta2_max)

 

!  ! logical mesh for space coordinates
!  mx => new_logical_mesh_2d( NPTS1, NPTS2)
print*,'pass 0'
  ! coordinate transformation associated with space coordinates
transformation => new_coordinate_transformation_2d_analytic( &
       "analytic_identity_transformation", &
       M, &
       identity_x1, &
       identity_x2, &
       identity_jac11, &
       identity_jac12, &
       identity_jac21, &
       identity_jac22, &
       params_mesh   )
! transformation => new_coordinate_transformation_2d_analytic( &
!       "analytic_polar_transformation", &
!       M, &
!       polar_x1, &
!       polar_x2, &
!       polar_jac11, &
!       polar_jac12, &
!       polar_jac21, &
!       polar_jac22, &
!       params_mesh  )     

! transformation => new_coordinate_transformation_2d_analytic( &
!       "analytic_collela_transformation", &
!       M, &
!       sinprod_x1, &
!       sinprod_x2, &
!       sinprod_jac11, &
!       sinprod_jac12, &
!       sinprod_jac21, &
!       sinprod_jac22, &
!       params_mesh  )

 print*,'pass 1'
  ! initialize simulation object with the above parameters
  call initialize_2d_gc_general( &
       simulation, &
       M, &
       transformation, &
       func_khp1, &
       landau_params, &
       params_field, &
       func_one, &
       func_zero, &
       func_zero, &
       func_one, &
       func_zero, &
       SPL_DEG1, & 
       SPL_DEG2, & 
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC )
       
  

  print *, ' f initialized '

  call simulation%run( )
  !call delete(simulation)
  call delete_2d_gc_general(simulation)
  call transformation%delete()
  
  print *, 'reached end of gc test'
  print *, 'PASSED'




end program gc_2d_general

! External functions used as parameters in the above unit test:


function func_one( eta1, eta2, params ) result(res)
  real(8), intent(in) :: eta1
  real(8), intent(in) :: eta2
  real(8), dimension(:), intent(in), optional :: params
  real(8) :: res
  res = 1.0_8
end function func_one

function func_minus_one( eta1, eta2, params ) result(res)
  real(8), intent(in) :: eta1
  real(8), intent(in) :: eta2
  real(8), dimension(:), intent(in), optional :: params
  real(8) :: res
  res = -1.0_8
end function func_minus_one

function func_zero( eta1, eta2, params ) result(res)
  real(8), intent(in) :: eta1
  real(8), intent(in) :: eta2
  real(8), dimension(:), intent(in), optional :: params
  real(8) :: res
  res = 0.0_8
end function func_zero


function func_khp1( x, y, params ) result(res)
  use sll_constants
  real(8), intent(in) :: x
  real(8), intent(in) :: y
  real(8) :: landau_mode,landau_alpha
  real(8), dimension(:), intent(in):: params
  real(8) :: res
  landau_mode = params(1)
  landau_alpha = params(2)
  res = sin(y)+ landau_alpha*cos(landau_mode*x)   
  !res = 0.001*cos(2*sll_pi*x)
 end function func_khp1
 
 function func_khp2( x, y, params ) result(res)   
  real(8), intent(in) :: x
  real(8), intent(in) :: y
  real(8) :: landau_mode,landau_alpha
  real(8), dimension(:), intent(in) :: params
  real(8) :: res
  landau_mode = params(1)
  landau_alpha = params(2) 
  res = sin(x)*sin(y)+ landau_alpha*sin(2*landau_mode*x)*sin(2*y)      
 end function func_khp2
  
 ! **************************************************************************
  !
  ! "Colella transformation";
  ! sinusoidal product (see P. Colella et al. JCP 230 (2011) formula 
  ! (102) p 2968):
  !
  !        x1 = eta1 + alpha * sin(eta1) * sin(eta2)
  !        x2 = eta2 + alpha * sin(eta1) * sin(eta2)
  !
  ! **************************************************************************

  ! direct mapping
  function sinprod2_x1 ( eta1, eta2,alpha )
    real(8)  :: sinprod2_x1
    real(8), intent(in)   :: eta1
    real(8), intent(in)   :: eta2
    real(8), intent(inout),optional :: alpha
    alpha = 0.1
    sinprod2_x1 = eta1 + alpha* sin(eta1) * sin(eta2)
  end function sinprod2_x1

  function sinprod2_x2 ( eta1, eta2,alpha )
    real(8)  :: sinprod2_x2
    real(8), intent(in)   :: eta1
    real(8), intent(in)   :: eta2
    real(8), intent(inout),optional :: alpha
    alpha = 0.1
    sinprod2_x2 = eta2 + alpha * sin(eta1) * sin(eta2)
  end function sinprod2_x2

  ! jacobian matrix
  function sinprod2_jac11 ( eta1, eta2,alpha )
    real(8)  :: sinprod2_jac11
    real(8), intent(in)   :: eta1
    real(8), intent(in)   :: eta2
    real(8), intent(inout),optional :: alpha
    alpha = 0.1
    sinprod2_jac11 = 1.0_8 + alpha * cos(eta1) * sin(eta2)
  end function sinprod2_jac11

    function sinprod2_jac12 ( eta1, eta2,alpha )
    real(8)  :: sinprod2_jac12
    real(8), intent(in)   :: eta1
    real(8), intent(in)   :: eta2
    real(8), intent(inout),optional :: alpha
    alpha = 0.1
    sinprod2_jac12 = alpha* sin(eta1) * cos(eta2)
  end function sinprod2_jac12

  function sinprod2_jac21 ( eta1, eta2,alpha )
    real(8)  :: sinprod2_jac21
    real(8), intent(in)   :: eta1
    real(8), intent(in)   :: eta2
    real(8), intent(inout),optional :: alpha
    alpha = 0.1
    sinprod2_jac21 = alpha * cos(eta1) * sin(eta2)
  end function sinprod2_jac21

  function sinprod2_jac22 ( eta1, eta2,alpha )
    real(8)  :: sinprod2_jac22
    real(8), intent(in)   :: eta1
    real(8), intent(in)   :: eta2
    real(8), intent(inout),optional :: alpha
    alpha = 0.1
    sinprod2_jac22 = 1.0_8 + alpha * sin(eta1)*cos(eta2)
  end function sinprod2_jac22

   ! jacobian ie determinant of jacobian matrix
  function sinprod2_jac ( eta1, eta2,alpha )
    real(8)  :: sinprod2_jac
    real(8), intent(in)   :: eta1
    real(8), intent(in)   :: eta2
    real(8), intent(inout),optional  :: alpha
    alpha = 0.1
    sinprod2_jac = 1.0_8 + alpha * sin(eta1+eta2) 
  end function sinprod2_jac
