program test_general_elliptic_solver
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_utilities.h"
#include "sll_file_io.h"

use sll_cartesian_meshes
use sll_module_coordinate_transformations_2d
use sll_common_coordinate_transformations
use sll_module_scalar_field_2d
use sll_constants
use sll_module_arbitrary_degree_spline_interpolator_2d
use sll_timer
use sll_module_deboor_splines_2d

#ifdef _UMFPACK
  use sll_general_coordinate_elliptic_solver_module_umfpack
#else
  use sll_general_coordinate_elliptic_solver_module
#endif

implicit none

#define SPLINE_DEG1       3
#define SPLINE_DEG2       3
#define NUM_CELLS1        4
#define NUM_CELLS2        4
#define ETA1MIN           0.0_f64
#define ETA1MAX           1.0_f64
#define ETA2MIN           0.0_f64
#define ETA2MAX           1.0_f64
#define PRINT_COMPARISON  .false.

type(sll_cartesian_mesh_2d), pointer                      :: mesh_2d
class(sll_coordinate_transformation_2d_base), pointer     :: T
type(general_coordinate_elliptic_solver)                  :: es
type(sll_arbitrary_degree_spline_interpolator_2d), target :: interp_2d
type(sll_arbitrary_degree_spline_interpolator_2d), target :: interp_2d_rhs
class(sll_scalar_field_2d_base), pointer                  :: a11_field_mat
class(sll_scalar_field_2d_base), pointer                  :: a12_field_mat
class(sll_scalar_field_2d_base), pointer                  :: a21_field_mat
class(sll_scalar_field_2d_base), pointer                  :: a22_field_mat
class(sll_scalar_field_2d_base), pointer                  :: b1_field_vect
class(sll_scalar_field_2d_base), pointer                  :: b2_field_vect
class(sll_scalar_field_2d_base), pointer                  :: c_field
class(sll_scalar_field_2d_base), pointer                  :: rho
type(sll_scalar_field_2d_discrete), pointer               :: phi
type(sll_time_mark)                                       :: t_reference

sll_real64 :: ti
sll_real64 :: te

sll_real64 :: acc
sll_real64 :: normL2
sll_real64 :: normH1

sll_real64, dimension(NUM_CELLS1+1,NUM_CELLS2+1) :: values
sll_real64, dimension(NUM_CELLS1+1,NUM_CELLS2+1) :: calculated
sll_real64, dimension(NUM_CELLS1+1,NUM_CELLS2+1) :: reference

sll_int32  :: i, j
sll_real64 :: h1,h2,node_val,ref
sll_real64 :: eta1(NUM_CELLS1+1)
sll_real64 :: eta2(NUM_CELLS2+1)
sll_int32  :: npts1,npts2

real(8) :: integral_solution
real(8) :: integral_exact_solution

sll_real64 :: grad1_node_val,grad2_node_val,grad1ref,grad2ref

mesh_2d => new_cartesian_mesh_2d( NUM_CELLS1, &
                                  NUM_CELLS2, &
                                  ETA1MIN,    &
                                  ETA1MAX,    &
                                  ETA2MIN,    &
                                  ETA2MAX )
acc   = 0.0_f64
npts1 =  NUM_CELLS1 + 1
npts2 =  NUM_CELLS2 + 1
h1    = (ETA1MAX-ETA1MIN)/real(NPTS1-1,f64)
h2    = (ETA2MAX-ETA2MIN)/real(NPTS2-1,f64)

do j=1,npts2
   do i=1,npts1
      eta1(i)  = (i-1)*h1 + ETA1MIN
      eta2(j)  = (j-1)*h2 + ETA2MIN
   end do
end do

normL2 = 0.0_f64
normH1 = 0.0_f64

print*, "-------------------------------------------------------------"
print*, " test case witout change of coordinates"
print*, " dirichlet-dirichlet boundary conditions"
print*, "-------------------------------------------------------------"

T => new_coordinate_transformation_2d_analytic( &
     "analytic",                                &
     mesh_2d,                                   &
     identity_x1,                               &
     identity_x2,                               &
     identity_jac11,                            &
     identity_jac12,                            &
     identity_jac21,                            &
     identity_jac22,                            &
     [0.0_f64]                                  )

call initialize_fields( SLL_DIRICHLET, SLL_DIRICHLET, SLL_DIRICHLET, SLL_DIRICHLET)

rho => new_scalar_field_2d_analytic( &
     func_zero,                      &
     "rho",                          &     
     T,                              &
     SLL_DIRICHLET,                  &
     SLL_DIRICHLET,                  &
     SLL_DIRICHLET,                  &
     SLL_DIRICHLET,                  &
     [0.0_f64]                       )

call solve_fields( SLL_DIRICHLET, SLL_DIRICHLET, &
                   SLL_DIRICHLET, SLL_DIRICHLET, ti, te)

call phi%write_to_file(0)

do j=1,npts2
do i=1,npts1
  node_val = calculated(i,j)
  grad1_node_val  = phi%first_deriv_eta1_value_at_point(eta1(i), eta2(j))
  grad2_node_val  = phi%first_deriv_eta2_value_at_point(eta1(i), eta2(j))
  ref             = 0.0_f64
  grad1ref        = 0.0_f64
  grad2ref        = 0.0_f64
  reference(i,j)  = ref
  normL2          = normL2 + (node_val-ref)**2*h1*h2
  normH1          = normH1 + ((grad1_node_val-grad1ref)**2+(grad2_node_val-grad2ref)**2)*h1*h2
  if (PRINT_COMPARISON) call printout_comparison()
end do
end do

integral_solution = sum(calculated(1:NUM_CELLS1,1:NUM_CELLS2))*h1*h2
integral_exact_solution = sum(reference(1:NUM_CELLS1,1:NUM_CELLS2))*h1*h2

call delete_things()
call check_error()

print*, 'PASSED'

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine printout_comparison()

print*,'(eta1,eta2) = ',eta1, eta2(j), 'calculated = ', node_val, &
       'theoretical = ',ref, 'difference=', node_val-ref
print*,'(eta1,eta2) = ',eta1, eta2(j), 'calculated = ', grad1_node_val, &
       'theoretical = ',grad1ref, 'difference=',grad1ref-grad1_node_val
print*,'(eta1,eta2) = ',eta1, eta2(j), 'calculated = ', grad2_node_val, &
       'theoretical = ',grad2ref, 'difference=',grad2ref-grad2_node_val

end subroutine printout_comparison

! Each field object must be initialized using the same logical
! mesh and coordinate transformation.
subroutine initialize_fields( bc1_min, bc1_max, bc2_min, bc2_max)

sll_int32, intent(in) :: bc1_min
sll_int32, intent(in) :: bc2_min
sll_int32, intent(in) :: bc1_max
sll_int32, intent(in) :: bc2_max

a11_field_mat => new_scalar_field_2d_analytic( &
  func_one,                                    &
  "a11",                                       &
  T,                                           &
  bc1_min,                                     &
  bc1_max,                                     &
  bc2_min,                                     &
  bc2_max,                                     &
  [0.0_f64]  ) 

a12_field_mat => new_scalar_field_2d_analytic( &
  func_zero,                                   &
  "a12",                                       &
  T,                                           &
  bc1_min,                                     &
  bc1_max,                                     &
  bc2_min,                                     &
  bc2_max,                                     &
  [0.0_f64] )

a21_field_mat => new_scalar_field_2d_analytic( &
  func_zero,                                   &
  "a21",                                       &
  T,                                           &
  bc1_min,                                     &
  bc1_max,                                     &
  bc2_min,                                     &
  bc2_max,                                     &
  [0.0_f64] ) 

a22_field_mat => new_scalar_field_2d_analytic( &
  func_one,                                    &
  "a22",                                       &
  T,                                           &
  bc1_min,                                     &
  bc1_max,                                     &
  bc2_min,                                     &
  bc2_max,                                     &
  [0.0_f64])

b1_field_vect => new_scalar_field_2d_analytic( &
  func_zero,                                   &
  "b1",                                        &
  T,                                           &
  bc1_min,                                     &
  bc1_max,                                     &
  bc2_min,                                     &
  bc2_max,                                     &
  [0.0_f64],                                   & 
  first_deriv_eta1 = func_zero,                &
  first_deriv_eta2 = func_zero) 

b2_field_vect => new_scalar_field_2d_analytic( &
  func_zero,                                   &
  "b2",                                        &
  T,                                           &
  bc1_min,                                     &
  bc1_max,                                     &
  bc2_min,                                     &
  bc2_max,                                     &
  [0.0_f64],                                   &
  first_deriv_eta1 = func_zero,                &
  first_deriv_eta2 = func_zero)

c_field => new_scalar_field_2d_analytic(       &
  func_zero,                                   &
  "c_field",                                   &
  T,                                           &
  bc1_min,                                     &
  bc1_max,                                     &
  bc2_min,                                     &
  bc2_max,                                     &
  [0.0_f64]  )

call initialize_ad2d_interpolator(             &
  interp_2d,                                   &
  NUM_CELLS1+1,                                &
  NUM_CELLS2+1,                                &
  ETA1MIN,                                     &
  ETA1MAX,                                     &
  ETA2MIN,                                     &
  ETA2MAX,                                     &
  bc1_min,                                     &
  bc1_max,                                     &
  bc2_min,                                     &
  bc2_max,                                     &
  SPLINE_DEG1,                                 &
  SPLINE_DEG2 )

call initialize_ad2d_interpolator(             &
  interp_2d_rhs,                               &
  NUM_CELLS1+1,                                &
  NUM_CELLS2+1,                                &
  ETA1MIN,                                     &
  ETA1MAX,                                     &
  ETA2MIN,                                     &
  ETA2MAX,                                     &
  bc1_min,                                     &
  bc1_max,                                     &
  bc2_min,                                     &
  bc2_max,                                     &
  SPLINE_DEG1,                                 &
  SPLINE_DEG2 )

phi => new_scalar_field_2d_discrete(           &
  "phi",                                       &
  interp_2d,                                   &
  T,                                           &
  bc1_min,                                     &
  bc1_max,                                     &
  bc2_min,                                     &
  bc2_max )

end subroutine initialize_fields

subroutine delete_things()

  call sll_delete(es)
  call rho%delete()
  call c_field%delete()
  call phi%delete()
  call a11_field_mat%delete()
  call a12_field_mat%delete()
  call a21_field_mat%delete()
  call b1_field_vect%delete()
  call b2_field_vect%delete()
  call a22_field_mat%delete()
  call T%delete()

end subroutine delete_things

subroutine solve_fields( bc1_min, &
                         bc1_max, &
                         bc2_min, &
                         bc2_max, &
                         ti,      &
                         te       )

sll_int32,  intent(in)  :: bc1_min
sll_int32,  intent(in)  :: bc2_min
sll_int32,  intent(in)  :: bc1_max
sll_int32,  intent(in)  :: bc2_max
sll_real64, intent(out) :: ti
sll_real64, intent(out) :: te

integral_solution = 0.0_f64
integral_exact_solution = 0.0_f64

call sll_set_time_mark(t_reference)

call sll_create(       &
  es,                  &
  SPLINE_DEG1,         &
  SPLINE_DEG2,         &
  NUM_CELLS1,          &
  NUM_CELLS2,          &
  ES_GAUSS_LEGENDRE,   &
  ES_GAUSS_LEGENDRE,   &
  bc1_min,             &
  bc1_max,             &
  bc2_min,             &
  bc2_max,             &
  ETA1MIN,             &
  ETA1MAX,             &
  ETA2MIN,             &
  ETA2MAX              )
 
call factorize_mat_es( &
  es,                  &
  a11_field_mat,       &
  a12_field_mat,       &
  a21_field_mat,       &
  a22_field_mat,       &
  b1_field_vect,       &
  b2_field_vect,       &
  c_field              )

ti = sll_time_elapsed_since(t_reference)

call sll_set_time_mark(t_reference)

values = 0.0_f64
call phi%set_field_data(values)
call phi%update_interpolation_coefficients()

call sll_solve( es, rho, phi)


te = sll_time_elapsed_since(t_reference)

do j=1,npts2
  do i=1,npts1
    calculated(i,j) = phi%value_at_point(eta1(i),eta2(j))
  end do
end do

integral_solution       = 0.0_f64
integral_exact_solution = 0.0_f64

end subroutine solve_fields

subroutine check_error()

print"('integral solution       =',g15.3)", integral_solution
print"('integral exact solution =',g15.3)", integral_exact_solution
acc = sum(abs(calculated-reference))/(npts1*npts2)
if (sqrt(normL2) <= h1**(SPLINE_DEG1-1)  .and. &
    sqrt(normH1) <= h1**(SPLINE_DEG1-2)) then     
   print"('error=',g15.3, 4x, 'OK' )", acc
else
  print*, ' L2 norm :', sqrt(normL2), h1**(SPLINE_DEG1-1)
  print*, ' H1 norm :', sqrt(normH1), h1**(SPLINE_DEG1-2)
  stop 'FAILED'
end if

end subroutine check_error

function func_one( eta1, eta2, params ) result(res)
real(8), intent(in) :: eta1
real(8), intent(in) :: eta2
real(8), dimension(:), intent(in) :: params
real(8) :: res
res = 1.0_f64
end function func_one

function func_zero( eta1, eta2, params ) result(res)
real(8), intent(in) :: eta1
real(8), intent(in) :: eta2
real(8), dimension(:), intent(in) :: params
real(8) :: res
res = 0.0_f64
end function func_zero

end program test_general_elliptic_solver
