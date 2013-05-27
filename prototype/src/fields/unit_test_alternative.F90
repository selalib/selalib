program unit_test_alternative
#include "sll_working_precision.h"
#include "sll_memory.h"
  use sll_logical_meshes
  use sll_constants
  use sll_module_scalar_field_2d_alternative
  use sll_module_coordinate_transformations_2d
  use sll_common_coordinate_transformations
  implicit none
  
  type(sll_logical_mesh_2d), pointer               :: mesh_2d
  class(sll_coordinate_transformation_2d_base), pointer :: T
  class(sll_scalar_field_2d_base), pointer              :: field_2d
  real(8), external :: test_function
  sll_int32 :: nc1, nc2, iplot
 ! procedure(polar_x1), pointer :: px1, px2, pjac11, pjac12, pjac21, pjac22
 ! type(init_landau_2d), target :: init_landau
 ! class(scalar_field_2d_initializer_base), pointer    :: pfinit
 ! type(cubic_spline_1d_interpolator), target  :: interp_eta1
 ! type(cubic_spline_1d_interpolator), target  :: interp_eta2
 ! class(sll_interpolator_1d_base), pointer :: interp_eta1_ptr
 ! class(sll_interpolator_1d_base), pointer :: interp_eta2_ptr

  ! logical mesh
  nc1 = 32
  nc2 = 32

  mesh_2d => new_logical_mesh_2d( nc1, nc2, &
       0.0_f64, 2.0*sll_pi, 0.0_f64,2.0*sll_pi )

  print *, 'initialized mesh 2D'

  ! coordinate transformation
  T => new_coordinate_transformation_2d_analytic( &
       "analytic", &
       mesh_2d, &
       identity_x1, &
       identity_x2, &
       identity_jac11, &
       identity_jac12, &
       identity_jac21, &
       identity_jac22 )
  print *, 'initialized transformation'

  field_2d => new_scalar_field_2d_analytic_alt( &
       test_function, &
       'doubly_periodic', &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC )
  print *, 'initialized field 2d'

  print *, 'field value at 0,0 = ', field_2d%value_at_point(0.0_f64,0.0_f64)
  print *, 'field value at indices 1,1 = ', &
       field_2d%value_at_indices(1,1)



!!$  do iplot = 1, 10
!!$     field%data = exp(-(mesh%x1_node**2+mesh%x2_node**2)*iplot*0.1)
!!$     call write_scalar_field_2d( field, multiply_by_jacobian=.true. )
!!$     call write_scalar_field_2d( field, multiply_by_jacobian=.true., output_file_name="field" )
!!$  end do
!!$
!!$  call mesh%write_to_file()

  print *, 'PASSED'
  
end program unit_test_alternative


function test_function( eta1, eta2, params ) result(res)
  intrinsic :: cos
  real(8) :: res
  real(8), intent(in) :: eta1
  real(8), intent(in) :: eta2
  real(8), dimension(:), intent(in), optional :: params
  res = 2.0*cos(eta1)*cos(eta2)
end function test_function


