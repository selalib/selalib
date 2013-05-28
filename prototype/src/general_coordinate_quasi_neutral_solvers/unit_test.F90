program test_general_qns
  use sll_logical_meshes
  use sll_module_coordinate_transformations_2d
  use sll_common_coordinate_transformations
  use sll_general_coordinate_qn_solver_module
  use sll_module_scalar_field_2d_alternative
  use sll_constants
  use sll_arbitrary_degree_spline_interpolator_2d_module
#include "sll_memory.h"
#include "sll_working_precision.h"
  implicit none

#define SPLINE_DEG1 1
#define SPLINE_DEG2 1
#define NUM_CELLS1  64
#define NUM_CELLS2  64
#define ETA1MIN  0.0_f64
#define ETA1MAX  (2.0*sll_pi)
#define ETA2MIN  0.0_f64
#define ETA2MAX  (2.0*sll_pi)

  type(sll_logical_mesh_2d), pointer          :: mesh_2d
  class(sll_coordinate_transformation_2d_base), pointer :: T
  type(general_coordinate_qn_solver), pointer :: qns
  type(arb_deg_2d_interpolator), target :: interp_2d
  class(sll_interpolator_2d_base), pointer :: interp_2d_ptr
!  class(sll_scalar_field_2d_analytic_alt), dimension(2,2) :: a_field_mat
  type(sll_scalar_field_2d_base_ptr), dimension(2,2) :: a_field_mat
  class(sll_scalar_field_2d_base), pointer    :: c_field
  class(sll_scalar_field_2d_base), pointer    :: rho
  class(sll_scalar_field_2d_base), pointer    :: phi
  real(8), external :: func_zero, func_one, source_term
  sll_real64, dimension(:,:), allocatable :: values
  sll_int32 :: ierr

  ! Table to represent the node values of phi
  SLL_ALLOCATE(values(NUM_CELLS1+1,NUM_CELLS2+1),ierr)

  ! First thing, initialize the logical mesh associated with this problem.
  mesh_2d => new_logical_mesh_2d( NUM_CELLS1, NUM_CELLS2, &
       0.0_f64, 2.0*sll_pi, 0.0_f64,2.0*sll_pi )

  ! Second, initialize the coordinate transformation associated with this 
  ! problem.
  T => new_coordinate_transformation_2d_analytic( &
       "analytic", &
       mesh_2d, &
       identity_x1, &
       identity_x2, &
       identity_jac11, &
       identity_jac12, &
       identity_jac21, &
       identity_jac22 )
  print *, 'initialized coordinate transformation'

  ! Thirdly, each field object must be initialized using the same logical
  ! mesh and coordinate transformation.
  a_field_mat(1,1)%base => new_scalar_field_2d_analytic_alt( &
       func_one, &
       "a11", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC ) 

  a_field_mat(1,2)%base => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "a12", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC ) 

  a_field_mat(2,1)%base => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "a21", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC ) 

  a_field_mat(2,2)%base => new_scalar_field_2d_analytic_alt( &
       func_one, &
       "a22", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC ) 

!!$  call initialize_scalar_field_2d_analytic_alt( &
!!$       a_field_mat(1,1), &
!!$       func_one, &
!!$       T, &
!!$       SLL_PERIODIC, &
!!$       SLL_PERIODIC, &
!!$       SLL_PERIODIC, &
!!$       SLL_PERIODIC )

!!$  call initialize_scalar_field_2d_analytic_alt( &
!!$       a_field_mat(1,2), &
!!$       func_zero, &
!!$       T, &
!!$       SLL_PERIODIC, &
!!$       SLL_PERIODIC, &
!!$       SLL_PERIODIC, &
!!$       SLL_PERIODIC )
!!$
!!$  call initialize_scalar_field_2d_analytic_alt( &
!!$       a_field_mat(2,1), &
!!$       func_zero, &
!!$       T, &
!!$       SLL_PERIODIC, &
!!$       SLL_PERIODIC, &
!!$       SLL_PERIODIC, &
!!$       SLL_PERIODIC )
!!$
!!$  call initialize_scalar_field_2d_analytic_alt( &
!!$       a_field_mat(2,2), &
!!$       func_one, &
!!$       T, &
!!$       SLL_PERIODIC, &
!!$       SLL_PERIODIC, &
!!$       SLL_PERIODIC, &
!!$       SLL_PERIODIC )

  c_field => new_scalar_field_2d_analytic_alt( &
       func_one, &
       "c_field", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC )

  rho => new_scalar_field_2d_analytic_alt( &
       source_term, &
       "rho", &     
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC )

  call initialize_ad2d_interpolator( &
       interp_2d, &
       NUM_CELLS1+1, &
       NUM_CELLS2+1, &
       ETA1MIN, &
       ETA1MAX, &
       ETA2MIN, &
       ETA2MAX, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SPLINE_DEG1, &
       SPLINE_DEG2 )

  interp_2d_ptr => interp_2d

  phi => new_scalar_field_2d_discrete_alt( &
       values, &
       "phi", &
       interp_2d_ptr, &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC )

  print *, 'initialized fields...'

  qns => new_general_qn_solver( &
       SPLINE_DEG1, &
       SPLINE_DEG2, &
       NUM_CELLS1, &
       NUM_CELLS2, &
       QNS_GAUSS_LEGENDRE, &
       QNS_GAUSS_LEGENDRE, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       ETA1MIN, &
       ETA1MAX, &
       ETA2MIN, &
       ETA2MAX )

  print *, 'Initialized QNS object'

  ! solve the field
  call solve_quasi_neutral_eq_general_coords( &
       qns, &
       a_field_mat, &
       c_field, &
       rho, &
       phi )

  print *, 'Completed solution'


  ! delete things...
  call delete(qns)
  print *, 'PASSED'
end program test_general_qns

function func_one( eta1, eta2, params ) result(res)
  real(8), intent(in) :: eta1
  real(8), intent(in) :: eta2
  real(8), dimension(:), intent(in), optional :: params
  real(8) :: res
  res = 1.0_8
end function func_one

function func_zero( eta1, eta2, params ) result(res)
  real(8), intent(in) :: eta1
  real(8), intent(in) :: eta2
  real(8), dimension(:), intent(in), optional :: params
  real(8) :: res
  res = 0.0_8
end function func_zero

function source_term( eta1, eta2, params ) result(res)
  intrinsic :: cos
  real(8), intent(in) :: eta1
  real(8), intent(in) :: eta2
  real(8), dimension(:), intent(in), optional :: params
  real(8) :: res
  res = 2.0*cos(eta1)*cos(eta2)
end function source_term
