program test_scalar_field_1d
#include "sll_working_precision.h"
#include "sll_memory.h"
use sll_m_cartesian_meshes
use sll_m_constants
use sll_m_scalar_field_1d
implicit none
  
#define SPLINE_DEG1 3
#define NUM_CELLS1  64
#define ETA1MIN  0.0_f64
#define ETA1MAX  1.0_f64
#define PRINT_COMPARISON .false.
  
type(sll_cartesian_mesh_1d), pointer                      :: mesh_1d
type(sll_arbitrary_degree_spline_interpolator_1d), target :: interp_1d

class(sll_scalar_field_1d_base), pointer :: periodic_analytic
class(sll_scalar_field_1d_base), pointer :: dirichlet_analytic
class(sll_scalar_field_1d_base), pointer :: periodic_discrete
class(sll_scalar_field_1d_base), pointer :: dirichlet_discrete

sll_int32                                :: nc1!, iplot
sll_real64                               :: grad1_node_val,grad1ref
sll_real64, dimension(:), allocatable    :: tab_values
sll_real64                               :: node_val,ref
sll_real64, dimension(:),   allocatable  :: point1
sll_real64                               :: eta1
sll_real64                               :: h1
sll_int32                                :: i


sll_real64 :: normL2_1,normL2_2,normL2_3,normL2_4
sll_real64 :: normH1_1,normH1_2,normH1_3,normH1_4

nc1 = NUM_CELLS1
h1 = (ETA1MAX-ETA1MIN)/real(nc1,f64)
  
! First thing, initialize the logical mesh associated with this problem.        
mesh_1d => new_cartesian_mesh_1d( NUM_CELLS1,ETA1MIN, ETA1MAX)
  
! --------------------------------------------------------------------------
!   Test case periodic analytic
!----------------------------------------------------------------------------
  
! ----> initialization of the field
periodic_analytic  => new_scalar_field_1d_analytic( &
  test_function_per, &
   "periodic_analytic",  &
  SLL_PERIODIC,      &
  SLL_PERIODIC,      &
  mesh_1d,           &
  first_derivative=test_function_per_der1)
  
! -------> compute error norm L2 and H1
normL2_1 = 0.0_f64
normH1_1 = 0.0_f64
do i=1,nc1+1
  eta1       = real(i-1,f64)*h1 + ETA1MIN
  node_val   = periodic_analytic%value_at_point(eta1)
  grad1_node_val = periodic_analytic%derivative_value_at_point(eta1)
  ref        = test_function_per(eta1)
  grad1ref   = test_function_per_der1(eta1)
  if(PRINT_COMPARISON) then
    print *, 'eta1 = ', eta1 , 'calculated = ', node_val, &
           'theoretical = ', ref, 'difference=', node_val-ref
    print *, 'eta1 = ', eta1, 'calculated = ', grad1_node_val, &
           'theoretical = ', grad1ref, 'difference=',grad1ref-grad1_node_val
        
  end if
     
  normL2_1    = normL2_1 + (node_val-ref)**2*h1
  normH1_1    = normH1_1 + ((grad1_node_val-grad1ref)**2)*h1
     
end do
print*, 'PASSED'
call periodic_analytic%delete()
  
! --------------------------------------------------------------------------
!   Test case dirichlet analytic
!----------------------------------------------------------------------------
  
dirichlet_analytic  => new_scalar_field_1d_analytic( &
  test_function_dir,                                 &
  "dirichlet_analytic",                                  &
  SLL_PERIODIC,                                      &
  SLL_PERIODIC,                                      &
  mesh_1d,                                           &
  first_derivative=test_function_dir_der1)
  
normL2_2 = 0.0_f64
normH1_2 = 0.0_f64
do i=1,nc1+1
  eta1           = real(i-1,f64)*h1 + ETA1MIN
  node_val       = dirichlet_analytic%value_at_point(eta1)
  grad1_node_val = dirichlet_analytic%derivative_value_at_point(eta1)
  ref            = test_function_dir(eta1)
  grad1ref   = test_function_dir_der1(eta1)
  if(PRINT_COMPARISON) then
    print *, 'eta1 = ', eta1, 'calculated = ', node_val, &
          'theoretical = ', ref, 'difference=', node_val-ref
    print *, 'eta1 = ', eta1, 'calculated = ', grad1_node_val, &
             'theoretical = ', grad1ref, 'difference=',grad1ref-grad1_node_val
        
  end if
        
  normL2_2 = normL2_2 + (node_val-ref)**2*h1
  normH1_2 = normH1_2 + ((grad1_node_val-grad1ref)**2)*h1
end do

call dirichlet_analytic%delete()
  
! --------------------------------------------------------------------------
!   Test case periodic non analytic
!----------------------------------------------------------------------------
allocate(point1(nc1 + 1))
allocate(tab_values(nc1 + 1))
do i=1,nc1 + 1
  point1(i)     = (i-1)*h1 + ETA1MIN 
  tab_values(i) = test_function_per(point1(i))
end do
  
call initialize_ad1d_interpolator( interp_1d,    &
                                   NUM_CELLS1+1, &
                                   ETA1MIN,      &
                                   ETA1MAX,      &
                                   SLL_PERIODIC, &
                                   SLL_PERIODIC, &
                                   SPLINE_DEG1)
  
periodic_discrete => new_scalar_field_1d_discrete( &
       "periodic_discrete", &
       interp_1d, &
       SLL_PERIODIC, &
       SLL_PERIODIC,&
       mesh_1d,&
       point1,&
       nc1+1)

! ------- > allocation values of field
call periodic_discrete%set_field_data(tab_values)
! --------> Compute coefficients of the field
call periodic_discrete%update_interpolation_coefficients( )

! -------> compute error norm L2 and H1
normL2_3 = 0.0_f64
normH1_3 = 0.0_f64
do i=1,nc1 + 1 
  eta1 = real(i-1,f64)*h1 + ETA1MIN 
  node_val       = periodic_discrete%value_at_point(eta1)
  grad1_node_val = periodic_discrete%derivative_value_at_point(eta1)
  ref        = test_function_per(eta1)
  grad1ref   = test_function_per_der1(eta1)
  if(PRINT_COMPARISON) then
     print *, 'eta1 = ', eta1, 'calculated = ', node_val, &
          'theoretical = ', ref, 'difference=', node_val-ref
     print *, 'eta1 = ', eta1, 'calculated = ', grad1_node_val, &
          'theoretical = ', grad1ref, 'difference=',grad1ref-grad1_node_val
  end if
     
  normL2_3    = normL2_3 + (node_val-ref)**2*h1
  normH1_3    = normH1_3 + ((grad1_node_val-grad1ref)**2)*h1
     
end do

call periodic_discrete%delete()

DEALLOCATE(point1)
DEALLOCATE(tab_values)
  
! --------------------------------------------------------------------------
!   Test case dirichlet non analytic
!----------------------------------------------------------------------------
  
allocate(point1(nc1 + 1))
allocate(tab_values(nc1 + 1))
do i=1,nc1 + 1
  point1(i)       = real(i-1,f64)*h1 + ETA1MIN 
  tab_values(i)   = test_function_dir(point1(i) )
end do

call initialize_ad1d_interpolator( &
       interp_1d, &
       NUM_CELLS1+1, &
       ETA1MIN, &
       ETA1MAX, &
       SLL_DIRICHLET,&
       SLL_DIRICHLET,&
       SPLINE_DEG1)
  
dirichlet_discrete => new_scalar_field_1d_discrete( &
       "dirichlet_discrete", &
       interp_1d, &
       SLL_DIRICHLET,&
       SLL_DIRICHLET,&
       mesh_1d,&
       point1,&
       nc1+1)
  
call dirichlet_discrete%set_field_data(tab_values)
call dirichlet_discrete%update_interpolation_coefficients( )
  
normL2_4 = 0.0_f64
normH1_4 = 0.0_f64
  
do i=1,nc1 + 1 
  eta1 = real(i-1,f64)*h1 + ETA1MIN 
  node_val       = dirichlet_discrete%value_at_point(eta1)
  grad1_node_val = dirichlet_discrete%derivative_value_at_point(eta1)
  ref        = test_function_dir(eta1)
  grad1ref   = test_function_dir_der1(eta1)
  if(PRINT_COMPARISON) then
     print *, 'eta1 = ', eta1, 'calculated = ', node_val, &
          'theoretical = ', ref, 'difference=', node_val-ref
     print *, 'eta1 = ', eta1, 'calculated = ', grad1_node_val, &
          'theoretical = ', grad1ref, 'difference=',grad1ref-grad1_node_val
     
  end if
  
  normL2_4    = normL2_4 + (node_val-ref)**2*h1
  normH1_4    = normH1_4 + ((grad1_node_val-grad1ref)**2)*h1
     
end do

call dirichlet_discrete%delete()

deallocate(point1)
deallocate(tab_values)

! **********************  TESTS **************************************
  
print*, '-------------------------------------------------------'
print*, ' PERIODIC ANALYTIC' 
print*, '-------------------------------------------------------'
print *,'Norm L2',sqrt(normL2_1),'Norm H1',sqrt(normH1_1),&
     h1**(SPLINE_DEG1),  h1**(SPLINE_DEG1-1)

print*, '-------------------------------------------------------'
print*, ' DIRICHLET ANALYTIC' 
print*, '-------------------------------------------------------'
print *,'Norm L2',sqrt(normL2_2),'Norm H1',sqrt(normH1_2),&
     h1**(SPLINE_DEG1),  h1**(SPLINE_DEG1-1)

print*, '-------------------------------------------------------'
print*, ' PERIODIC NON ANALYTIC' 
print*, '-------------------------------------------------------'
print *,'Norm L2',sqrt(normL2_3),'Norm H1',sqrt(normH1_3),&
     h1**(SPLINE_DEG1),  h1**(SPLINE_DEG1-1)

print*, '-------------------------------------------------------'
print*, ' DIRICHLET NON ANALYTIC' 
print*, '-------------------------------------------------------'
print *,'Norm L2',sqrt(normL2_4),'Norm H1',sqrt(normH1_4),&
     h1**(SPLINE_DEG1),  h1**(SPLINE_DEG1-1)

if ( ( sqrt(normL2_1) <= h1**(SPLINE_DEG1))   .AND. &
     ( sqrt(normL2_2) <= h1**(SPLINE_DEG1))   .AND. &
     ( sqrt(normL2_3) <= h1**(SPLINE_DEG1))   .AND. &
     ( sqrt(normL2_4) <= h1**(SPLINE_DEG1))   .AND. &
     ( sqrt(normH1_1) <= h1**(SPLINE_DEG1-1)) .AND. &
     ( sqrt(normH1_2) <= h1**(SPLINE_DEG1-1)) .AND. &
     ( sqrt(normH1_3) <= h1**(SPLINE_DEG1-1)) .AND. &
     ( sqrt(normH1_4) <= h1**(SPLINE_DEG1-1))) then
   print *, 'PASSED'
end if

contains

function test_function_per( eta1, params ) result(res)
   sll_real64 :: res
   sll_real64, intent(in) :: eta1
   sll_real64, dimension(:), intent(in), optional :: params
  intrinsic :: cos
  res = cos(2*sll_pi*eta1)
end function test_function_per


function test_function_per_der1( eta1, params) result(res)
  intrinsic :: cos,sin
  sll_real64 :: res
  sll_real64, intent(in) :: eta1
  sll_real64, dimension(:), intent(in), optional :: params

  res = -2*sll_pi*sin(2*sll_pi*eta1)
end function test_function_per_der1


function test_function_dir( eta1, params) result(res)
  intrinsic :: cos,sin
  sll_real64 :: res
  sll_real64, intent(in) :: eta1
   sll_real64, dimension(:), intent(in), optional :: params

  res = sin(2*sll_pi*eta1)
end function test_function_dir

function test_function_dir_der1( eta1, params) result(res)
  intrinsic :: cos
  sll_real64 :: res
  sll_real64, intent(in) :: eta1
   sll_real64, dimension(:), intent(in), optional :: params

  res = 2.0*sll_pi*cos(2*sll_pi*eta1)
end function test_function_dir_der1

end program test_scalar_field_1d

