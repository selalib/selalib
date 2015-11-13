!Old test for general elliptic solver made by Aurore
program test_general_elliptic_solver
#include "sll_memory.h"
#include "sll_working_precision.h"

use sll_m_cartesian_meshes
use sll_m_coordinate_transformations_2d
use sll_m_common_coordinate_transformations
use sll_m_scalar_field_2d
use sll_m_constants
use sll_m_arbitrary_degree_spline_interpolator_2d
use sll_m_timer
use sll_m_deboor_splines_2d
use m_init_functions, only: &
  func_zero, &
  func_four, &
  func_one, &
  func_epsi, &
  source_term_perper, &
  source_term_perdir, &
  source_term_dirper, &
  source_term_chgt_perper, &
  source_term_chgt_perdir, &
  source_term_chgt_dirper, &
  source_term_chgt_dirdir, &
  f_sin, u_sin, u_sin_der1, u_sin_der2, &
  f_cos, u_cos, u_cos_der1, u_cos_der2, &
  sol_exacte_perper, &
  sol_exacte_perper_der1, &
  sol_exacte_perper_der2, &
  sol_exacte_perdir  , &
  sol_exacte_perdir_der1, &
  sol_exacte_perdir_der2, &
  sol_exacte_dirper, &
  sol_exacte_dirper_der1, &
  sol_exacte_dirper_der2, &
  sol_exacte_chgt_perper, &
  sol_exacte_chgt_perper_der1, &
  sol_exacte_chgt_perper_der2, &
  sol_exacte_chgt_perdir  , &
  sol_exacte_chgt_perdir_der1, &
  sol_exacte_chgt_perdir_der2, &
  sol_exacte_chgt_dirper, &
  sol_exacte_chgt_dirper_der1, &
  sol_exacte_chgt_dirper_der2, &
  sol_exacte_chgt_dirdir, &
  sol_exacte_chgt_dirdir_der1, &
  sol_exacte_chgt_dirdir_der2, &
  adimension_chgt_x, &
  adimension_chgt_y, &
  jac11_adimension_chgt, &
  jac12_adimension_chgt, &
  jac21_adimension_chgt, &
  jac22_adimension_chgt, &
  sol_exacte_chgt_adim, &
  source_term_chgt_adim



#ifdef _UMFPACK
  use sll_general_coordinate_elliptic_solver_module_umfpack
#else
  use sll_m_general_coordinate_elliptic_solver
#endif

implicit none

#define SPLINE_DEG1       3
#define SPLINE_DEG2       3
#define NUM_CELLS1        64
#define NUM_CELLS2        64
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
class(sll_interpolator_2d_base), pointer                  :: rhs_interp
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

sll_real64 :: ti(16), te(16)
sll_real64 :: acc(16)    
sll_real64 :: normL2(16)
sll_real64 :: normH1(16)

sll_real64, dimension(NUM_CELLS1+1,NUM_CELLS2+1) :: values
sll_real64, dimension(NUM_CELLS1+1,NUM_CELLS2+1) :: calculated
sll_real64, dimension(NUM_CELLS1+1,NUM_CELLS2+1) :: reference
sll_real64, dimension(NUM_CELLS1+1,NUM_CELLS2+1) :: tab_rho

sll_real64 :: val_jac
sll_int32  :: i, j, k
sll_real64 :: h1,h2,node_val,ref
sll_real64 :: eta1(NUM_CELLS1+1)
sll_real64 :: eta2(NUM_CELLS2+1)
sll_int32  :: npts1,npts2
!sll_int32  :: ierr

real(8) :: integral_solution
real(8) :: integral_exact_solution

character(len=10) :: cmd
integer           :: itest1
integer           :: itest2
character(len=4)  :: ccase
!sll_int32         :: file_id


sll_real64 :: grad1_node_val,grad2_node_val,grad1ref,grad2ref
sll_real64, dimension(1) :: whatever  ! dummy params array
character(len=49) :: case_name(16)

case_name = ['(per-per) with identity and source term analytic', &
             '(per-dir) with identity and source term analytic', &
             '(dir-dir) with identity and source term analytic', &
             '(dir-per) with identity and source term analytic', &
             '(per-per) with colella  and source term analytic', &
             '(per-dir) with colella  and source term analytic', &
             '(dir-dir) with colella  and source term analytic', &
             '(dir-per) with colella  and source term analytic', &
             '(per-per) with identity and source term discrete', &
             '(per-per) with colella  and source term discrete', &
             '(per-dir) with colella  and source term discrete', &
             '(dir-dir) with colella  and source term discrete', &
             '(dir-per) with colella  and source term discrete', &
             '(dir-per) with polar    and source term analytic', &
             '(dir-per) with polar    and source term analytic', &
             '(dir-dir) with identity and source term analytic']

! First thing, initialize the logical mesh associated with this problem. 
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

call get_command_argument(1,cmd)
read (cmd,'(I2)') itest1
call get_command_argument(2,cmd)
read (cmd,'(I2)') itest2
print *, itest1, itest2

if (itest1 ==0) itest1 = 01
if (itest2 ==0) itest2 = 13

do k = itest1, itest2

  call int2string(k, ccase)

  select case(k)

  case(1)

  print*, "-------------------------------------------------------------"
  print*, "1 test case witout change of coordinates"
  print*, " periodic-periodic boundary conditions "
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
       (/ 0.0_f64 /) )

  call initialize_fields( SLL_PERIODIC, SLL_PERIODIC, SLL_PERIODIC, SLL_PERIODIC)

  rho => new_scalar_field_2d_analytic( &
       source_term_perper,             &
       "rho"//ccase,                   &     
       T,                              &
       SLL_PERIODIC,                   &
       SLL_PERIODIC,                   &
       SLL_PERIODIC,                   &
       SLL_PERIODIC,                   &
       whatever  )
  
  call solve_fields( SLL_PERIODIC, SLL_PERIODIC, &
                     SLL_PERIODIC, SLL_PERIODIC, ti(k), te(k))

  do j=1,npts2
  do i=1,npts1
    node_val        = calculated(i,j)
    grad1_node_val  = phi%first_deriv_eta1_value_at_point(eta1(i),eta2(j))
    grad2_node_val  = phi%first_deriv_eta2_value_at_point(eta1(i),eta2(j))
    ref             = sol_exacte_perper(eta1(i),eta2(j))
    grad1ref        = sol_exacte_perper_der1(eta1(i),eta2(j))
    grad2ref        = sol_exacte_perper_der2(eta1(i),eta2(j))
    reference( i,j) = ref
    normL2(k)       = normL2(k)+(node_val-ref)**2*h1*h2
    normH1(k)       = normH1(k)+((grad1_node_val-grad1ref)**2 &
                                +(grad2_node_val-grad2ref)**2)*h1*h2
  end do
  end do

  integral_solution = sum(calculated(1:NUM_CELLS1,1:NUM_CELLS2))*h1*h2
  integral_exact_solution = sum(reference(1:NUM_CELLS1,1:NUM_CELLS2))*h1*h2

  call delete_things()
  
  call check_error(k)

  case(2)
  print*, "-------------------------------------------------------------"
  print*, " 2 test case witout change of coordinates"
  print*, " periodic-dirichlet boundary conditions"
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
       (/0.0_f64/) )

  call initialize_fields( SLL_PERIODIC, SLL_PERIODIC, SLL_DIRICHLET, SLL_DIRICHLET)

  rho => new_scalar_field_2d_analytic( &
       source_term_perdir,             &
       "rho"//ccase,                   &     
       T,                              &
       SLL_PERIODIC,                   &
       SLL_PERIODIC,                   &
       SLL_DIRICHLET,                  &
       SLL_DIRICHLET,                  &
       whatever )

  call solve_fields( SLL_PERIODIC,  SLL_PERIODIC,  &
                     SLL_DIRICHLET, SLL_DIRICHLET, &
                     ti(k), te(k))
  
  do j=1,npts2
  do i=1,npts1
    node_val       = calculated(i,j)
    grad1_node_val = phi%first_deriv_eta1_value_at_point(eta1(i), eta2(j))
    grad2_node_val = phi%first_deriv_eta2_value_at_point(eta1(i), eta2(j))
    ref            = sol_exacte_perdir(eta1(i),eta2(j))
    grad1ref       = sol_exacte_perdir_der1(eta1(i),eta2(j))
    grad2ref       = sol_exacte_perdir_der2(eta1(i),eta2(j))
    reference(i,j) = ref
    normL2(k)      = normL2(k)+(node_val-ref)**2*h1*h2
    normH1(k)      = normH1(k)+((grad1_node_val-grad1ref)**2+(grad2_node_val-grad2ref)**2)*h1*h2
    if(PRINT_COMPARISON) call printout_comparison()
  end do
  end do
  
  integral_solution = sum(calculated(1:NUM_CELLS1,1:NUM_CELLS2))*h1*h2
  integral_exact_solution = sum(reference(1:NUM_CELLS1,1:NUM_CELLS2))*h1*h2

  call delete_things()
  call check_error(k)

  case(3)
  print*, "-------------------------------------------------------------"
  print*, " 3 test case witout change of coordinates"
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
       (/0.0_f64/) )
  
  call initialize_fields( SLL_DIRICHLET, SLL_DIRICHLET, SLL_DIRICHLET, SLL_DIRICHLET)

  rho => new_scalar_field_2d_analytic( &
       source_term_perdir,             &
       "rho"//ccase,                   &     
       T,                              &
       SLL_DIRICHLET,                  &
       SLL_DIRICHLET,                  &
       SLL_DIRICHLET,                  &
       SLL_DIRICHLET,                  &
       whatever )
  
  call solve_fields( SLL_DIRICHLET, SLL_DIRICHLET, &
                     SLL_DIRICHLET, SLL_DIRICHLET, ti(k), te(k))

  do j=1,npts2
  do i=1,npts1
    node_val = calculated(i,j)
    grad1_node_val  = phi%first_deriv_eta1_value_at_point(eta1(i), eta2(j))
    grad2_node_val  = phi%first_deriv_eta2_value_at_point(eta1(i), eta2(j))
    ref             = sol_exacte_perdir(eta1(i),eta2(j))
    grad1ref        = sol_exacte_perdir_der1(eta1(i),eta2(j))
    grad2ref        = sol_exacte_perdir_der2(eta1(i),eta2(j))
    reference(i,j)  = ref
    normL2(k)       = normL2(k) + (node_val-ref)**2*h1*h2
    normH1(k)       = normH1(k) + ((grad1_node_val-grad1ref)**2+(grad2_node_val-grad2ref)**2)*h1*h2
    if (PRINT_COMPARISON) call printout_comparison()
  end do
  end do
  
  integral_solution = sum(calculated(1:NUM_CELLS1,1:NUM_CELLS2))*h1*h2
  integral_exact_solution = sum(reference(1:NUM_CELLS1,1:NUM_CELLS2))*h1*h2

  call delete_things()
  call check_error(k)

  case(4)
  print*, "-------------------------------------------------------------"
  print*, " 4 test case witout change of coordinates"
  print*, " dirichlet-periodic boundary conditions"
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
       (/0.0_f64/) )
  
  call initialize_fields( SLL_DIRICHLET, SLL_DIRICHLET, &
                          SLL_PERIODIC,  SLL_PERIODIC)
  
  rho => new_scalar_field_2d_analytic( &
       source_term_dirper,             &
       "rho"//ccase,                   &     
       T,                              &
       SLL_DIRICHLET,                  &
       SLL_DIRICHLET,                  &
       SLL_PERIODIC,                   &
       SLL_PERIODIC,                   &
       whatever )
  
  call solve_fields( SLL_DIRICHLET, SLL_DIRICHLET, &
                     SLL_PERIODIC, SLL_PERIODIC, ti(k), te(k))
  
  do j=1,npts2
  do i=1,npts1
    node_val = calculated(i,j)
    node_val        = phi%value_at_point(eta1(i),eta2(j))
    grad1_node_val  = phi%first_deriv_eta1_value_at_point(eta1(i), eta2(j))
    grad2_node_val  = phi%first_deriv_eta2_value_at_point(eta1(i), eta2(j))
    ref             = sol_exacte_dirper(eta1(i),eta2(j))
    grad1ref        = sol_exacte_dirper_der1(eta1(i),eta2(j))
    grad2ref        = sol_exacte_dirper_der2(eta1(i),eta2(j))
    reference(i,j)  = ref
    normL2(k)       = normL2(k) + (node_val-ref)**2*h1*h2
    normH1(k)       = normH1(k) + &
    ((grad1_node_val-grad1ref)**2+(grad2_node_val-grad2ref)**2)*h1*h2
    if(PRINT_COMPARISON) call printout_comparison()
  end do
  end do

  integral_solution = sum(calculated(1:NUM_CELLS1,1:NUM_CELLS2))*h1*h2
  integral_exact_solution = sum(reference(1:NUM_CELLS1,1:NUM_CELLS2))*h1*h2

  call delete_things()
  
  call check_error(k)

  case(5)
  print*, "-----------------------------------------------"
  print*, " 5 test case with colella change of coordinates"
  print*, " periodic-periodic boundary conditions         "
  print*, "-----------------------------------------------"
  
  T => new_coordinate_transformation_2d_analytic( &
       "analytic",                                &
       mesh_2d,                                   &
       sinprod_x1,                                &
       sinprod_x2,                                &
       sinprod_jac11,                             &
       sinprod_jac12,                             &
       sinprod_jac21,                             &
       sinprod_jac22,                             &
       (/ 0.1_f64, 0.1_f64, 1.0_f64, 1.0_f64/) )
  
  call initialize_fields( SLL_PERIODIC, SLL_PERIODIC, &
                          SLL_PERIODIC, SLL_PERIODIC)

  rho => new_scalar_field_2d_analytic( &
       source_term_chgt_perper,        &
       "rho"//ccase,                   &     
       T,                              &
       SLL_PERIODIC,                   &
       SLL_PERIODIC,                   &
       SLL_PERIODIC,                   &
       SLL_PERIODIC,                   &
       whatever )
  
  call solve_fields( SLL_PERIODIC, SLL_PERIODIC, &
                     SLL_PERIODIC, SLL_PERIODIC, ti(k), te(k))
  
  do j=1,npts2
  do i=1,npts1
    node_val = calculated(i,j)
    grad1_node_val = phi%first_deriv_eta1_value_at_point(eta1(i), eta2(j))
    grad2_node_val = phi%first_deriv_eta2_value_at_point(eta1(i), eta2(j))
    ref        = sol_exacte_chgt_perper(eta1(i),eta2(j))
    grad1ref   = sol_exacte_chgt_perper_der1(eta1(i),eta2(j))
    grad2ref   = sol_exacte_chgt_perper_der2(eta1(i),eta2(j))
    reference(i,j) = ref
    val_jac = sinprod_jac11(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))*&
              sinprod_jac22(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))-&
              sinprod_jac12(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))*&
              sinprod_jac21(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))
    if(PRINT_COMPARISON) call printout_comparison()
    if ( i < NUM_CELLS1 .and. j < NUM_CELLS2) then
       integral_solution = integral_solution + node_val*val_jac * h1*h2
       integral_exact_solution = integral_exact_solution + ref*val_jac * h1*h2
       normL2(k)    = normL2(k) + (node_val-ref)**2*h1*h2*val_jac
       normH1(k)    = normH1(k) + &
      ((grad1_node_val-grad1ref)**2+(grad2_node_val-grad2ref)**2)*h1*h2*val_jac
    end if
  end do
  end do
  
  integral_solution = sum(calculated(1:NUM_CELLS1,1:NUM_CELLS2))*h1*h2
  integral_exact_solution = sum(reference(1:NUM_CELLS1,1:NUM_CELLS2))*h1*h2

  call delete_things()
  
  call check_error(k)

  case(6)
  print*, "-------------------------------------------------------------"
  print*, " 6 test case with colella change of coordinates"
  print*, " periodic-dirichlet boundary conditions"
  print*, "-------------------------------------------------------------"
  
  T => new_coordinate_transformation_2d_analytic( &
       "analytic",                                &
       mesh_2d,                                   &
       sinprod_x1,                                &
       sinprod_x2,                                &
       sinprod_jac11,                             &
       sinprod_jac12,                             &
       sinprod_jac21,                             &
       sinprod_jac22,                             &
       (/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))
  
  call initialize_fields( SLL_PERIODIC,  SLL_PERIODIC, SLL_DIRICHLET, SLL_DIRICHLET)
  
  rho => new_scalar_field_2d_analytic( &
       source_term_chgt_perdir,        &
       "rho"//ccase,                   &     
       T,                              &
       SLL_PERIODIC,                   &
       SLL_PERIODIC,                   &
       SLL_DIRICHLET,                  &
       SLL_DIRICHLET,                  &
       whatever )
  
  call solve_fields( SLL_PERIODIC, SLL_PERIODIC, &
                     SLL_DIRICHLET, SLL_DIRICHLET, ti(k), te(k))

  do j=1,npts2
  do i=1,npts1
    node_val = calculated(i,j)
    grad1_node_val = phi%first_deriv_eta1_value_at_point(eta1(i), eta2(j))
    grad2_node_val = phi%first_deriv_eta2_value_at_point(eta1(i), eta2(j))
    ref        = sol_exacte_chgt_perdir(eta1(i),eta2(j))
    grad1ref   = sol_exacte_chgt_perdir_der1(eta1(i),eta2(j))
    grad2ref   = sol_exacte_chgt_perdir_der2(eta1(i),eta2(j))
    reference(i,j) = ref
    if(PRINT_COMPARISON) call printout_comparison()
    
    val_jac = sinprod_jac11(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))*&
              sinprod_jac22(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))-&
              sinprod_jac12(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))*&
              sinprod_jac21(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))
    if ( i < NUM_CELLS1 .and. j < NUM_CELLS2 ) then
       integral_solution = integral_solution + node_val*val_jac* h1*h2
       integral_exact_solution = integral_exact_solution + ref*val_jac* h1*h2
       normL2(k)    = normL2(k) + (node_val-ref)**2*h1*h2
       normH1(k)    = normH1(k) + &
      ((grad1_node_val-grad1ref)**2+(grad2_node_val-grad2ref)**2)*h1*h2*val_jac
       
    end if
  end do
  end do
  
  call delete_things()
  
  case(7)
  print*, "-------------------------------------------------------------"
  print*, " 7 test case with colella change of coordinates"
  print*, " dirichlet-dirichlet boundary conditions"
  print*, "-------------------------------------------------------------"
  
  T => new_coordinate_transformation_2d_analytic( &
       "analytic",    &
       mesh_2d,       &
       sinprod_x1,    &
       sinprod_x2,    &
       sinprod_jac11, &
       sinprod_jac12, &
       sinprod_jac21, &
       sinprod_jac22, &
       (/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))
  
  call initialize_fields( SLL_DIRICHLET, SLL_DIRICHLET, &
                          SLL_DIRICHLET, SLL_DIRICHLET)

  rho => new_scalar_field_2d_analytic( &
       source_term_chgt_dirdir,        &
       "rho"//ccase,                   &     
       T,                              &
       SLL_DIRICHLET,                  &
       SLL_DIRICHLET,                  &
       SLL_DIRICHLET,                  &
       SLL_DIRICHLET,                  &
       whatever )
  
  call solve_fields( SLL_DIRICHLET, SLL_DIRICHLET, &
                     SLL_DIRICHLET, SLL_DIRICHLET, ti(k), te(k))
  
  do j=1,npts2
     do i=1,npts1
        
        node_val = calculated(i,j)
        grad1_node_val = phi%first_deriv_eta1_value_at_point(eta1(i), eta2(j))
        grad2_node_val = phi%first_deriv_eta2_value_at_point(eta1(i), eta2(j))
        ref        = sol_exacte_chgt_dirdir(eta1(i),eta2(j))
        grad1ref   = sol_exacte_chgt_dirdir_der1(eta1(i),eta2(j))
        grad2ref   = sol_exacte_chgt_dirdir_der2(eta1(i),eta2(j))
        reference(i,j) = ref
        if(PRINT_COMPARISON) call printout_comparison()
        val_jac = sinprod_jac11(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))*&
                  sinprod_jac22(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))-&
                  sinprod_jac12(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))*&
                  sinprod_jac21(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))
        if ( i < NUM_CELLS1 .and. j < NUM_CELLS2) then
           integral_solution = integral_solution + node_val*val_jac * h1*h2
           integral_exact_solution = integral_exact_solution + ref*val_jac * h1*h2
           normL2(k)    = normL2(k) + (node_val-ref)**2*h1*h2*val_jac
           normH1(k)    = normH1(k) + &
         ((grad1_node_val-grad1ref)**2+(grad2_node_val-grad2ref)**2)*h1*h2*val_jac
           
        end if
     end do
  end do
  call delete_things()

  call check_error(k)

  case(8)
  print*, "---------------------"
  print*, " 8 test case with colella change of coordinates"
  print*, " dirichlet-periodic boundary conditions"
  print*, "---------------------"
  
  T => new_coordinate_transformation_2d_analytic( &
       "analytic",                                &
       mesh_2d,                                   &
       sinprod_x1,                                &
       sinprod_x2,                                &
       sinprod_jac11,                             &
       sinprod_jac12,                             &
       sinprod_jac21,                             &
       sinprod_jac22,                             &
       (/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))

  call initialize_fields( SLL_DIRICHLET, SLL_DIRICHLET, &
                          SLL_PERIODIC,  SLL_PERIODIC)
  
  rho => new_scalar_field_2d_analytic( &
       source_term_chgt_dirper,        &
       "rho"//ccase,                   &     
       T,                              &
       SLL_DIRICHLET,                  &
       SLL_DIRICHLET,                  &
       SLL_PERIODIC,                   &
       SLL_PERIODIC,                   &
       whatever)
  
  call solve_fields( SLL_DIRICHLET, SLL_DIRICHLET, &
                     SLL_PERIODIC,  SLL_PERIODIC, ti(k), te(k))
  
  do j=1,npts2
     do i=1,npts1
        
        node_val = calculated(i,j)
        grad1_node_val = phi%first_deriv_eta1_value_at_point(eta1(i), eta2(j))
        grad2_node_val = phi%first_deriv_eta2_value_at_point(eta1(i), eta2(j))
        ref        = sol_exacte_chgt_dirper(eta1(i),eta2(j))
        grad1ref   = sol_exacte_chgt_dirper_der1(eta1(i),eta2(j))
        grad2ref   = sol_exacte_chgt_dirper_der2(eta1(i),eta2(j))
        reference(i,j) = ref
        if(PRINT_COMPARISON) call printout_comparison()
        val_jac = sinprod_jac11(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))*&
                  sinprod_jac22(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))-&
                  sinprod_jac12(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))*&
                  sinprod_jac21(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))
        if ( i < NUM_CELLS1 .and. j < NUM_CELLS2) then
           integral_solution = integral_solution + node_val*val_jac * h1*h2
           integral_exact_solution = integral_exact_solution + ref*val_jac * h1*h2
           normL2(k)    = normL2(k) + (node_val-ref)**2*h1*h2*val_jac
           normH1(k)    = normH1(k) + &
          ((grad1_node_val-grad1ref)**2+(grad2_node_val-grad2ref)**2)*h1*h2*val_jac
        end if
     end do
  end do

  call delete_things()

  call check_error(k)

  case(9)
  print*, "---------------------"
  print*, " 9 test case without change of coordinates"
  print*, " periodic-periodic boundary conditions"
  print*, " with non analytic source term " 
  print*, "---------------------"
  
  T => new_coordinate_transformation_2d_analytic( &
       "analytic",                                &
       mesh_2d,                                   &
       identity_x1,                               &
       identity_x2,                               &
       identity_jac11,                            &
       identity_jac12,                            &
       identity_jac21,                            &
       identity_jac22,                            &
       (/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))
  
  call initialize_fields( SLL_PERIODIC, SLL_PERIODIC, &
                          SLL_PERIODIC, SLL_PERIODIC)

  do j=1,npts2
    do i=1,npts1
      tab_rho(i,j) = source_term_perper( eta1(i), eta2(j), whatever )
    end do
  end do
  
  rho => new_scalar_field_2d_discrete( &
       "rho"//ccase,                   &
       interp_2d_rhs,                  &
       T,                              &
       SLL_PERIODIC,                   &
       SLL_PERIODIC,                   &
       SLL_PERIODIC,                   &
       SLL_PERIODIC,                   &
       eta1,                           &
       NUM_CELLS1,                     &
       eta2,                           &
       NUM_CELLS2)  

  call rho%set_field_data(tab_rho)
  call rho%update_interpolation_coefficients()

  call solve_fields( SLL_PERIODIC, SLL_PERIODIC, &
                     SLL_PERIODIC, SLL_PERIODIC, ti(k), te(k))
  
  do j=1,npts2
  do i=1,npts1
     
    node_val = calculated(i,j)
    grad1_node_val = phi%first_deriv_eta1_value_at_point(eta1(i), eta2(j))
    grad2_node_val = phi%first_deriv_eta2_value_at_point(eta1(i), eta2(j))
    ref        = sol_exacte_perper(eta1(i),eta2(j))
    grad1ref   = sol_exacte_perper_der1(eta1(i),eta2(j))
    grad2ref   = sol_exacte_perper_der2(eta1(i),eta2(j))
    reference(i,j) = ref
    val_jac = 1.0_f64
    if(PRINT_COMPARISON) call printout_comparison()
    if ( i < NUM_CELLS1 .and. j < NUM_CELLS2) then
       integral_solution = integral_solution + node_val * h1*h2
       integral_exact_solution = integral_exact_solution + ref * h1*h2
       normL2(k)    = normL2(k) + (node_val-ref)**2*h1*h2*val_jac
       normH1(k)    = normH1(k) + &
      ((grad1_node_val-grad1ref)**2+(grad2_node_val-grad2ref)**2)*h1*h2*val_jac
    end if

  end do
  end do

  call delete_things()

  call check_error(k)

  case(10)
  print*, "------------------------------------------------"
  print*, " 10 test case with colella change of coordinates"
  print*, " periodic-periodic boundary conditions          "
  print*, " with non analytic source term                  " 
  print*, "------------------------------------------------"
  
  T => new_coordinate_transformation_2d_analytic( &
       "analytic",                                &
       mesh_2d,                                   &
       sinprod_x1,                                &
       sinprod_x2,                                &
       sinprod_jac11,                             &
       sinprod_jac12,                             &
       sinprod_jac21,                             &
       sinprod_jac22,                             &
       (/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))
  
  call initialize_fields( SLL_PERIODIC, SLL_PERIODIC, &
                          SLL_PERIODIC, SLL_PERIODIC)
  
  do j=1,npts2
  do i=1,npts1
    tab_rho(i,j) = source_term_chgt_perper( eta1(i), eta2(j), whatever )
  end do
  end do
  
  rhs_interp => interp_2d_rhs
  tab_rho(:,:) = tab_rho - sum(tab_rho)/real(NUM_CELLS1*NUM_CELLS2,f64)

  rho => new_scalar_field_2d_discrete( &
       "rho"//ccase,                   &
       rhs_interp,                     &
       T,                              &
       SLL_PERIODIC,                   &
       SLL_PERIODIC,                   &
       SLL_PERIODIC,                   &
       SLL_PERIODIC ,                  &
       eta1,                           &
       NUM_CELLS1,                     &
       eta2,                           &
       NUM_CELLS2)

  call rho%set_field_data(tab_rho)
  call rho%update_interpolation_coefficients()
  
  call solve_fields( SLL_PERIODIC, SLL_PERIODIC, &
                     SLL_PERIODIC, SLL_PERIODIC, ti(k), te(k))
  
  do j=1,npts2
  do i=1,npts1
     
    node_val = calculated(i,j)
    grad1_node_val = phi%first_deriv_eta1_value_at_point(eta1(i), eta2(j))
    grad2_node_val = phi%first_deriv_eta2_value_at_point(eta1(i), eta2(j))
    ref        = sol_exacte_chgt_perper(eta1(i),eta2(j))
    grad1ref   = sol_exacte_chgt_perper_der1(eta1(i),eta2(j))
    grad2ref   = sol_exacte_chgt_perper_der2(eta1(i),eta2(j))
    reference(i,j) = ref
    if(PRINT_COMPARISON) call printout_comparison()
    val_jac = sinprod_jac11(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))*&
              sinprod_jac22(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))-&
              sinprod_jac12(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))*&
              sinprod_jac21(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))
    
    if ( i < NUM_CELLS1 .and. j < NUM_CELLS2) then
       integral_solution = integral_solution + node_val*val_jac * h1*h2
       integral_exact_solution = integral_exact_solution + ref*val_jac * h1*h2
       normL2(k)    = normL2(k) + (node_val-ref)**2*h1*h2*val_jac
       normH1(k)    = normH1(k) + &
      ((grad1_node_val-grad1ref)**2+(grad2_node_val-grad2ref)**2)*h1*h2*val_jac
    end if

  end do
  end do
  
  call delete_things()

  call check_error(k)

  case(11)
  print*, "------------------------------------------------"
  print*, " 11 test case with colella change of coordinates"
  print*, " periodic-dirichlet boundary conditions         "
  print*, " with non analytic source term                  " 
  print*, "------------------------------------------------"
  
  T => new_coordinate_transformation_2d_analytic( &
       "analytic",                                &
       mesh_2d,                                   &
       sinprod_x1,                                &
       sinprod_x2,                                &
       sinprod_jac11,                             &
       sinprod_jac12,                             &
       sinprod_jac21,                             &
       sinprod_jac22,                             &
       (/0.1_f64, 0.1_f64, 1.0_f64, 1.0_f64/)) 

  call initialize_fields( SLL_PERIODIC, SLL_PERIODIC, &
                          SLL_DIRICHLET, SLL_DIRICHLET)

  do j=1,npts2
  do i=1,npts1
    tab_rho(i,j) = source_term_chgt_perdir( eta1(i), eta2(j), whatever )
  end do
  end do

  rhs_interp => interp_2d_rhs

  rho => new_scalar_field_2d_discrete( &
       "rho"//ccase,                   &
       rhs_interp,                     &
       T,                              &
       SLL_PERIODIC,                   &
       SLL_PERIODIC,                   &
       SLL_DIRICHLET,                  &
       SLL_DIRICHLET,                  &
       eta1,                           &
       NUM_CELLS1,                     &
       eta2,                           &
       npts2)

  call rho%set_field_data(tab_rho)
  call rho%update_interpolation_coefficients()
  
  call solve_fields( SLL_PERIODIC, SLL_PERIODIC, &
                     SLL_DIRICHLET, SLL_DIRICHLET, ti(k), te(k))
  
  do j=1,npts2
  do i=1,npts1
    
    node_val = calculated(i,j)
    grad1_node_val = phi%first_deriv_eta1_value_at_point(eta1(i), eta2(j))
    grad2_node_val = phi%first_deriv_eta2_value_at_point(eta1(i), eta2(j))
    ref        = sol_exacte_chgt_perdir(eta1(i),eta2(j))
    grad1ref   = sol_exacte_chgt_perdir_der1(eta1(i),eta2(j))
    grad2ref   = sol_exacte_chgt_perdir_der2(eta1(i),eta2(j))
    reference(i,j) = ref
    if(PRINT_COMPARISON) call printout_comparison()
    val_jac = sinprod_jac11(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))*&
              sinprod_jac22(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))-&
              sinprod_jac12(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))*&
              sinprod_jac21(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))
    
    if ( i < NUM_CELLS1 .and. j < NUM_CELLS2) then
       integral_solution = integral_solution + node_val*val_jac * h1*h2
       integral_exact_solution = integral_exact_solution + &
            ref*val_jac * h1*h2
       normL2(k)    = normL2(k) + (node_val-ref)**2*h1*h2*val_jac
       normH1(k)   = normH1(k) + ((grad1_node_val-grad1ref)**2+&
            (grad2_node_val-grad2ref)**2)*h1*h2*val_jac
       
    end if
  end do
  end do

  call delete_things()
  
  call check_error(k)

  case(12)

  print*, "------------------------------------------------"
  print*, " 12 test case with colella change of coordinates"
  print*, " dirichlet-dirichlet boundary conditions        "
  print*, " with non analytic source term                  " 
  print*, "------------------------------------------------"
  
   T => new_coordinate_transformation_2d_analytic( &
       "analytic",                                 &
       mesh_2d,                                    &
       sinprod_x1,                                 &
       sinprod_x2,                                 &
       sinprod_jac11,                              &
       sinprod_jac12,                              &
       sinprod_jac21,                              &
       sinprod_jac22,                              &
       (/0.1_f64, 0.1_f64, 1.0_f64, 1.0_f64/) )

  call initialize_fields( SLL_DIRICHLET, SLL_DIRICHLET, &
                          SLL_DIRICHLET, SLL_DIRICHLET)

  do j=1,npts2
  do i=1,npts1
    tab_rho(i,j)  = source_term_chgt_dirdir( eta1(i), eta2(j), whatever )
  end do
  end do

  rhs_interp => interp_2d_rhs

  rho => new_scalar_field_2d_discrete( &
       "rho"//ccase,                   &
       rhs_interp,                     &
       T,                              &
       SLL_DIRICHLET,                  &
       SLL_DIRICHLET,                  &
       SLL_DIRICHLET,                  &
       SLL_DIRICHLET,                  &
       eta1,                           &
       npts1,                          &
       eta2,                           &
       npts2)

  call rho%set_field_data(tab_rho)
  call rho%update_interpolation_coefficients()
  
  call solve_fields( SLL_DIRICHLET, SLL_DIRICHLET, &
                     SLL_DIRICHLET, SLL_DIRICHLET, ti(k), te(k))
  
  do j=1,npts2
  do i=1,npts1
        
    node_val       = calculated(i,j)
    grad1_node_val = phi%first_deriv_eta1_value_at_point(eta1(i), eta2(j))
    grad2_node_val = phi%first_deriv_eta2_value_at_point(eta1(i), eta2(j))
    ref            = sol_exacte_chgt_dirdir(eta1(i),eta2(j))
    grad1ref       = sol_exacte_chgt_dirdir_der1(eta1(i),eta2(j))
    grad2ref       = sol_exacte_chgt_dirdir_der2(eta1(i),eta2(j))
    reference(i,j) = ref

    if (PRINT_COMPARISON) call printout_comparison()
    val_jac = sinprod_jac11(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))*&
              sinprod_jac22(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))-&
              sinprod_jac12(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))*&
              sinprod_jac21(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))
    
    if ( i < NUM_CELLS1 .and. j < NUM_CELLS2) then
       integral_solution = integral_solution + node_val*val_jac * h1*h2
       integral_exact_solution = integral_exact_solution + &
            ref*val_jac * h1*h2
       normL2(k)    = normL2(k) + (node_val-ref)**2*h1*h2*val_jac
       normH1(k)    = normH1(k) + ((grad1_node_val-grad1ref)**2+&
            (grad2_node_val-grad2ref)**2)*h1*h2*val_jac
    end if

  end do
  end do

  call delete_things()

  call check_error(k)

  case(13)

  print*, "------------------------------------------------"
  print*, " 13 test case with colella change of coordinates"
  print*, " dirichlet-periodic  boundary conditions        "
  print*, " with non analytic source term                  " 
  print*, "------------------------------------------------"
  
  T => new_coordinate_transformation_2d_analytic( &
       "analytic",                                &
       mesh_2d,                                   &
       sinprod_x1,                                &
       sinprod_x2,                                &
       sinprod_jac11,                             &
       sinprod_jac12,                             &
       sinprod_jac21,                             &
       sinprod_jac22,                             &
       (/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))
  
  call initialize_fields( SLL_DIRICHLET, SLL_DIRICHLET, &
                          SLL_PERIODIC,  SLL_PERIODIC)

  do j=1,npts2
     do i=1,npts1
        tab_rho(i,j)  = source_term_chgt_dirper( eta1(i), eta2(j), whatever )
     end do
  end do
  
  rhs_interp => interp_2d_rhs

  rho => new_scalar_field_2d_discrete( &
       "rho"//ccase,                   &
       rhs_interp,                     &
       T,                              &
       SLL_DIRICHLET,                  &
       SLL_DIRICHLET,                  &
       SLL_PERIODIC,                   &
       SLL_PERIODIC,                   &
       eta1,                           &
       npts1,                          &
       eta2,                           &
       NUM_CELLS2)

  call rho%set_field_data(tab_rho)
  call rho%update_interpolation_coefficients()

  call solve_fields( SLL_DIRICHLET, SLL_DIRICHLET,  &
                     SLL_PERIODIC,  SLL_PERIODIC, ti(k), te(k))
   
  do j=1,npts2
  do i=1,npts1
     
    node_val = calculated(i,j)
    grad1_node_val = phi%first_deriv_eta1_value_at_point(eta1(i), eta2(j))
    grad2_node_val = phi%first_deriv_eta2_value_at_point(eta1(i), eta2(j))
    ref        = sol_exacte_chgt_dirper(eta1(i),eta2(j))
    grad1ref   = sol_exacte_chgt_dirper_der1(eta1(i),eta2(j))
    grad2ref   = sol_exacte_chgt_dirper_der2(eta1(i),eta2(j))
    reference(i,j)  = ref
    if (PRINT_COMPARISON) call printout_comparison()
    val_jac = sinprod_jac11(eta1(i),eta2(j),(/.1_f64,.1_f64,1._f64,1._f64/))*&
              sinprod_jac22(eta1(i),eta2(j),(/.1_f64,.1_f64,1._f64,1._f64/))-&
              sinprod_jac12(eta1(i),eta2(j),(/.1_f64,.1_f64,1._f64,1._f64/))*&
              sinprod_jac21(eta1(i),eta2(j),(/.1_f64,.1_f64,1._f64,1._f64/))
    
    if ( i < NUM_CELLS1 .and. j < NUM_CELLS2) then
       integral_solution = integral_solution + node_val*val_jac * h1*h2
       integral_exact_solution = integral_exact_solution + &
            ref*val_jac * h1*h2
       normL2(k)    = normL2(k) + (node_val-ref)**2*h1*h2*val_jac
       normH1(k)    = normH1(k) + ((grad1_node_val-grad1ref)**2+&
            (grad2_node_val-grad2ref)**2)*h1*h2*val_jac
       
    end if
  end do
  end do

  call delete_things()
  call check_error(k)

  case(14)
  print*, "--------------------------------------------------"
  print*, " 14 test case with polar change of coordinates    "
  print*, " dirichlet-periodic  boundary conditions          "
  print*, " with analytic source term                        " 
  print*, "--------------------------------------------------"
  
  T => new_coordinate_transformation_2d_analytic( &
       "polar",                                   &
       mesh_2d,                                   &
       x1_polar_f,                                &
       x2_polar_f,                                &
       deriv_x1_polar_f_eta1,                     & 
       deriv_x1_polar_f_eta2,                     &
       deriv_x2_polar_f_eta1,                     &
       deriv_x2_polar_f_eta2,                     &
       [1.0_f64,2.0_f64] )
  
  call initialize_fields( SLL_DIRICHLET, SLL_DIRICHLET, &
                          SLL_PERIODIC,  SLL_PERIODIC)

  rho => new_scalar_field_2d_analytic( &
       f_sin,                          &
       "fsin",                         &
       T,                              &
       SLL_DIRICHLET,                  &
       SLL_DIRICHLET,                  &
       SLL_PERIODIC,                   &
       SLL_PERIODIC,                   &
       whatever)

  do j = 1, npts2
  do i = 1, npts1
    values(i,j) = u_sin( eta1(i), eta2(j) )
  end do
  end do

  call phi%set_field_data(values)

  call solve_fields( SLL_DIRICHLET, SLL_DIRICHLET, &
                     SLL_PERIODIC,  SLL_PERIODIC, ti(k), te(k))
   
  do j=1,npts2
  do i=1,npts1
        
    node_val       = calculated(i,j)
    grad1_node_val = phi%first_deriv_eta1_value_at_point(eta1(i), eta2(j))
    grad2_node_val = phi%first_deriv_eta2_value_at_point(eta1(i), eta2(j))
    ref            = u_sin(eta1(i),eta2(j))
    grad1ref       = u_sin_der1(eta1(i),eta2(j))
    grad2ref       = u_sin_der2(eta1(i),eta2(j))
    reference(i,j) = ref

    if (PRINT_COMPARISON) call printout_comparison()

    val_jac = deriv_x1_polar_f_eta1(eta1(i),eta2(j),[1.0_f64,2.0_f64])*&
              deriv_x2_polar_f_eta2(eta1(i),eta2(j),[1.0_f64,2.0_f64])-&
              deriv_x1_polar_f_eta2(eta1(i),eta2(j),[1.0_f64,2.0_f64])*&
              deriv_x2_polar_f_eta1(eta1(i),eta2(j),[1.0_f64,2.0_f64])
      
    if ( i < NUM_CELLS1 .and. j < NUM_CELLS2) then
      integral_solution = integral_solution + node_val*val_jac * h1*h2
      integral_exact_solution = integral_exact_solution + ref*val_jac * h1*h2
      normL2(k)    = normL2(k) + (node_val-ref)**2*h1*h2*val_jac

      !PN error on derivative is disabled because grad1ref and grad2ref are not 
      !PN computed correctly
      !normH1(k) = normH1(k) + ((grad1_node_val-grad1ref)**2+&
      !                         (grad2_node_val-grad2ref)**2)*h1*h2*val_jac

    end if

  end do
  end do

  call delete_things()
  call check_error(k)

  case(15)

  print*, "--------------------------------------------------"
  print*, " 15 test case with polar change of coordinates    "
  print*, " dirichlet-periodic  boundary conditions          "
  print*, " with analytic source term                        " 
  print*, "--------------------------------------------------"
  
  T => new_coordinate_transformation_2d_analytic( &
       "polar",                                   &
       mesh_2d,                                   &
       x1_polar_f,                                &
       x2_polar_f,                                &
       deriv_x1_polar_f_eta1,                     &
       deriv_x1_polar_f_eta2,                     &
       deriv_x2_polar_f_eta1,                     &
       deriv_x2_polar_f_eta2,                     &
       [1.0_f64,2.0_f64] )
  
  call initialize_fields( SLL_DIRICHLET, SLL_DIRICHLET, SLL_PERIODIC,  SLL_PERIODIC)

  rho => new_scalar_field_2d_analytic( &
       f_cos,                          & 
       "f_cos",                        &
       T,                              &
       SLL_DIRICHLET,                  &
       SLL_DIRICHLET,                  &
       SLL_PERIODIC,                   &
       SLL_PERIODIC,                   &
       whatever)


  call solve_fields( SLL_DIRICHLET, SLL_DIRICHLET, &
                     SLL_PERIODIC,  SLL_PERIODIC, ti(k), te(k))
   
  do j=1,npts2
  do i=1,npts1
        
    node_val       = calculated(i,j)
    grad1_node_val = phi%first_deriv_eta1_value_at_point(eta1(i), eta2(j))
    grad2_node_val = phi%first_deriv_eta2_value_at_point(eta1(i), eta2(j))
    ref            = u_cos(eta1(i),eta2(j))
    grad1ref       = u_cos_der1(eta1(i),eta2(j))
    grad2ref       = u_cos_der2(eta1(i),eta2(j))
    reference(i,j) = ref

    if (PRINT_COMPARISON) call printout_comparison()

    val_jac = deriv_x1_polar_f_eta1(eta1(i),eta2(j),[1.0_f64,2.0_f64])*&
              deriv_x2_polar_f_eta2(eta1(i),eta2(j),[1.0_f64,2.0_f64])-&
              deriv_x1_polar_f_eta2(eta1(i),eta2(j),[1.0_f64,2.0_f64])*&
              deriv_x2_polar_f_eta1(eta1(i),eta2(j),[1.0_f64,2.0_f64])
      
    if ( i < NUM_CELLS1 .and. j < NUM_CELLS2) then
      integral_solution       = integral_solution + node_val*val_jac * h1*h2
      integral_exact_solution = integral_exact_solution + ref*val_jac * h1*h2
      normL2(k) = normL2(k) + (node_val-ref)**2*h1*h2*val_jac
!PN error on derivative is disabled because grad1ref and grad2ref are not 
!PN computed correctly
!      normH1(k) = normH1(k) + ((grad1_node_val-grad1ref)**2+&
!                               (grad2_node_val-grad2ref)**2)*h1*h2*val_jac
    end if

  end do
  end do

  call delete_things()

  call check_error(k)

  case(16)
  print*, "--------------------------------------------------"
  print*, " 16 test case with                                "
  print*, " dirichlet-dirichlet  boundary conditions         "
  print*, " with source term = 4                             " 
  print*, "--------------------------------------------------"
  
  T => new_coordinate_transformation_2d_analytic( &
       "analytic",                                &
       mesh_2d,                                   &
       identity_x1,                               &
       identity_x2,                               &
       identity_jac11,                            &
       identity_jac12,                            &
       identity_jac21,                            &
       identity_jac22,                            &
       (/0.0_f64/) )
  
  call initialize_fields( SLL_DIRICHLET, SLL_DIRICHLET, &
                          SLL_DIRICHLET, SLL_DIRICHLET)

  rho => new_scalar_field_2d_analytic( &
       func_four,                      & 
       "func_four",                    &
       T,                              &
       SLL_DIRICHLET,                  &
       SLL_DIRICHLET,                  &
       SLL_DIRICHLET,                  &
       SLL_DIRICHLET,                  &
       whatever)

  do j = 1, npts2
  do i = 1, npts1
    values(i,j) = eta1(i)*eta1(i) + eta2(j)*eta2(j)
  end do
  end do

!  call phi%interp_2d%set_values_at_boundary(values(1,:),     &
!                                            values(npts1,:), &
!                                            values(:,1),     &
!                                            values(:,npts2))

  call solve_fields( SLL_DIRICHLET, SLL_DIRICHLET, &
                     SLL_DIRICHLET, SLL_DIRICHLET, ti(k), te(k))
   
  do j=1,npts2
  do i=1,npts1
        
    node_val       = calculated(i,j)
    grad1_node_val = phi%first_deriv_eta1_value_at_point(eta1(i), eta2(j))
    grad2_node_val = phi%first_deriv_eta2_value_at_point(eta1(i), eta2(j))
    ref            = eta1(i)*eta1(i)+eta2(j)*eta2(j)
    grad1ref       = 2.0_f64*eta1(i)
    grad2ref       = 2.0_f64*eta2(j)
    reference(i,j) = ref

    if (PRINT_COMPARISON) call printout_comparison()

    integral_solution       = integral_solution + node_val*h1*h2
    integral_exact_solution = integral_exact_solution + ref*h1*h2
    normL2(k) = normL2(k) + (node_val-ref)**2*h1*h2
    normH1(k) = normH1(k) + ((grad1_node_val-grad1ref)**2+&
                             (grad2_node_val-grad2ref)**2)*h1*h2
  end do
  end do

  call delete_things()

  call check_error(k)

  end select

end do

!call sll_ascii_file_create("solutions_gces.gnu",file_id,ierr)
do k = itest1, itest2
  write(*,"(a)") case_name(k)
  print"('test',i2,' : ','norm L2=',g15.3,' norm H1=',g15.3,' times=',2g15.3)" &
    ,k,normL2(k),normH1(k),ti(k),te(k)
!  call int2string(k, ccase)
!  write(file_id,"(a)") "set title '"//case_name(k)//"'"
!  write(file_id,"(a)") "load 'phi_"//ccase//".gnu'"
!  write(file_id,"(a)") " pause -1"
end do

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
subroutine initialize_fields( bc_eta1_min, bc_eta1_max, bc_eta2_min, bc_eta2_max)

  sll_int32, intent(in) :: bc_eta1_min
  sll_int32, intent(in) :: bc_eta2_min
  sll_int32, intent(in) :: bc_eta1_max
  sll_int32, intent(in) :: bc_eta2_max

  a11_field_mat => new_scalar_field_2d_analytic(     &
    func_one,                                        &
    "a11",                                           &
    T,                                               &
    bc_eta1_min,                                     &
    bc_eta1_max,                                     &
    bc_eta2_min,                                     &
    bc_eta2_max,                                     &
    whatever  ) 
  
  a12_field_mat => new_scalar_field_2d_analytic(     &
    func_zero,                                       &
    "a12",                                           &
    T,                                               &
    bc_eta1_min,                                     &
    bc_eta1_max,                                     &
    bc_eta2_min,                                     &
    bc_eta2_max,                                     &
    whatever )
  
  a21_field_mat => new_scalar_field_2d_analytic(     &
    func_zero,                                       &
    "a21",                                           &
    T,                                               &
    bc_eta1_min,                                     &
    bc_eta1_max,                                     &
    bc_eta2_min,                                     &
    bc_eta2_max,                                     &
    whatever ) 
  
  a22_field_mat => new_scalar_field_2d_analytic(     &
    func_one,                                        &
    "a22",                                           &
    T,                                               &
    bc_eta1_min,                                     &
    bc_eta1_max,                                     &
    bc_eta2_min,                                     &
    bc_eta2_max,                                     &
    whatever)

  b1_field_vect => new_scalar_field_2d_analytic(     &
    func_zero,                                       &
    "b1",                                            &
    T,                                               &
    bc_eta1_min,                                     &
    bc_eta1_max,                                     &
    bc_eta2_min,                                     &
    bc_eta2_max,                                     &
    whatever,                                        & 
    first_deriv_eta1 = func_zero,                    &
    first_deriv_eta2 = func_zero) 

  b2_field_vect => new_scalar_field_2d_analytic(     &
    func_zero,                                       &
    "b2",                                            &
    T,                                               &
    bc_eta1_min,                                     &
    bc_eta1_max,                                     &
    bc_eta2_min,                                     &
    bc_eta2_max,                                     &
    whatever,                                        &
    first_deriv_eta1 = func_zero,                    &
    first_deriv_eta2 = func_zero)

  c_field => new_scalar_field_2d_analytic(           &
    func_zero,                                       &
    "c_field",                                       &
    T,                                               &
    bc_eta1_min,                                     &
    bc_eta1_max,                                     &
    bc_eta2_min,                                     &
    bc_eta2_max,                                     &
    whatever  )

  call initialize_ad2d_interpolator(                 &
    interp_2d,                                       &
    NUM_CELLS1+1,                                    &
    NUM_CELLS2+1,                                    &
    ETA1MIN,                                         &
    ETA1MAX,                                         &
    ETA2MIN,                                         &
    ETA2MAX,                                         &
    bc_eta1_min,                                     &
    bc_eta1_max,                                     &
    bc_eta2_min,                                     &
    bc_eta2_max,                                     &
    SPLINE_DEG1,                                     &
    SPLINE_DEG2 )

  call initialize_ad2d_interpolator(                 &
    interp_2d_rhs,                                   &
    NUM_CELLS1+1,                                    &
    NUM_CELLS2+1,                                    &
    ETA1MIN,                                         &
    ETA1MAX,                                         &
    ETA2MIN,                                         &
    ETA2MAX,                                         &
    bc_eta1_min,                                     &
    bc_eta1_max,                                     &
    bc_eta2_min,                                     &
    bc_eta2_max,                                     &
    SPLINE_DEG1,                                     &
    SPLINE_DEG2 )

  phi => new_scalar_field_2d_discrete(               &
    "phi_"//ccase,                                   &
    interp_2d,                                       &
    T,                                               &
    bc_eta1_min,                                     &
    bc_eta1_max,                                     &
    bc_eta2_min,                                     &
    bc_eta2_max )

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

subroutine solve_fields( bc_eta1_min, &
                         bc_eta1_max, &
                         bc_eta2_min, &
                         bc_eta2_max, &
                         ti,          &
                         te           )

sll_int32,  intent(in)  :: bc_eta1_min
sll_int32,  intent(in)  :: bc_eta2_min
sll_int32,  intent(in)  :: bc_eta1_max
sll_int32,  intent(in)  :: bc_eta2_max
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
  bc_eta1_min,         &
  bc_eta1_max,         &
  bc_eta2_min,         &
  bc_eta2_max,         &
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

!call phi%write_to_file(0)

te = sll_time_elapsed_since(t_reference)

do j=1,npts2
  do i=1,npts1
    calculated(i,j) = phi%value_at_point(eta1(i),eta2(j))
  end do
end do

integral_solution       = 0.0_f64
integral_exact_solution = 0.0_f64

end subroutine solve_fields

subroutine check_error(icase)

integer, intent(in) :: icase
print"('integral solution       =',g15.3)", integral_solution
print"('integral exact solution =',g15.3)", integral_exact_solution
acc(icase) = sum(abs(calculated-reference))/real(npts1*npts2,f64)
if ((sqrt(normL2(icase)) <= h1**(SPLINE_DEG1-1))   .AND. &
    (sqrt(normH1(icase)) <= h1**(SPLINE_DEG1-2))) then     
   print"('test:',i2,4x,'error=',g15.3, 4x, 'OK' )", icase, acc(icase)
else
  print*, ' L2 norm :', sqrt(normL2(icase)), h1**(SPLINE_DEG1-1)
  print*, ' H1 norm :', sqrt(normH1(icase)), h1**(SPLINE_DEG1-2)
  stop 'FAILED'
end if

end subroutine check_error

end program test_general_elliptic_solver
