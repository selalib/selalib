program unit_test_1d_alternative
#include "sll_working_precision.h"
#include "sll_memory.h"
  use sll_logical_meshes
  use sll_constants
  use sll_module_scalar_field_1d_alternative
  implicit none
  
#define SPLINE_DEG1 3
#define NUM_CELLS1  64
#define ETA1MIN  0.0_f64
#define ETA1MAX  1.0_f64
#define PRINT_COMPARISON .false.
  
  type(sll_logical_mesh_1d), pointer               :: mesh_1d
  ! either of these type declarations can be used to work. Initialization is
  ! different.
  class(sll_scalar_field_1d_base), pointer              :: field_1d

  class(sll_scalar_field_1d_base), pointer              :: periodic_anal
  class(sll_scalar_field_1d_base), pointer              :: dirichlet_anal
  class(sll_scalar_field_1d_base), pointer              :: periodic_discrete
  class(sll_scalar_field_1d_base), pointer              :: dirichlet_discrete
  sll_int32 :: nc1, iplot
  sll_real64 :: grad1_node_val,grad1ref
  sll_real64, dimension(:), pointer :: tab_values
  type(arb_deg_1d_interpolator), target                 :: interp_1d
  sll_real64 :: node_val,ref

  sll_real64, dimension(:),   allocatable    :: point1
  sll_real64 :: eta1
  sll_real64  :: h1
  sll_int32 :: npts1,npts2
  sll_int32 :: i,j
  sll_int32 :: ierr
  real(8), external :: test_function_per
  real(8), external :: test_function_per_der1
  real(8), external :: test_function_dir
  real(8), external :: test_function_dir_der1

  sll_real64 :: normL2_1,normL2_2,normL2_3,normL2_4
  sll_real64 :: normH1_1,normH1_2,normH1_3,normH1_4
  ! logical mesh
  nc1 = NUM_CELLS1
  h1 = (ETA1MAX-ETA1MIN)/real(nc1,f64)
  print *, 'h1 = ', h1
  


  ! First thing, initialize the logical mesh associated with this problem.        
  mesh_1d => new_logical_mesh_1d( NUM_CELLS1,ETA1MIN, ETA1MAX)
  
  print *, 'initialized mesh 1D'
  
  
  
  ! ******************************************************************
  ! ------------------ TEST ANALYTIC ------------------------------  
  ! ******************************************************************
  
!!! --------------------------------------------------------------------------
  !   Test case periodic analytic
  !----------------------------------------------------------------------------
  
  
  ! ----> initialization of the field
  periodic_anal  => new_scalar_field_1d_analytic_alt( &
       test_function_per, &
       "periodic_anal", &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       mesh_1d, &
       first_derivative=test_function_per_der1)
  
  print *, 'initialized field 1d'
  
  ! -------> compute error norm L2 and H1
  normL2_1 = 0.0_f64
  normH1_1 = 0.0_f64
  do i=1,nc1+1
     eta1       = real(i-1,f64)*h1 + ETA1MIN
     node_val   = periodic_anal%value_at_point(eta1)
    
     grad1_node_val = periodic_anal%derivative_value_at_point(eta1)
     
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
   !print*, normL2_1
    print*, 'passed'
   ! -------> field visualization 
   !call periodic_anal%write_to_file(0)
   
  ! -------> delete field
  call periodic_anal%delete()
  
!!! --------------------------------------------------------------------------
  !   Test case dirichlet analytic
  !----------------------------------------------------------------------------
  
  ! ----> initialization of the field
  dirichlet_anal  => new_scalar_field_1d_analytic_alt( &
       test_function_dir, &
       "dirichlet_anal", &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       mesh_1d, &
       first_derivative=test_function_dir_der1)
  print *, 'initialized field 1d'
  
  ! -------> compute error norm L2 and H1
  normL2_2 = 0.0_f64
  normH1_2 = 0.0_f64
  do i=1,nc1+1
     eta1           = real(i-1,f64)*h1 + ETA1MIN
     node_val       = dirichlet_anal%value_at_point(eta1)
     grad1_node_val = dirichlet_anal%derivative_value_at_point(eta1)
     ref            = test_function_dir(eta1)
     grad1ref   = test_function_dir_der1(eta1)
     if(PRINT_COMPARISON) then
        print *, 'eta1 = ', eta1, 'calculated = ', node_val, &
             'theoretical = ', ref, 'difference=', node_val-ref
        print *, 'eta1 = ', eta1, 'calculated = ', grad1_node_val, &
             'theoretical = ', grad1ref, 'difference=',grad1ref-grad1_node_val
        
     end if
        
     normL2_2    = normL2_2 + (node_val-ref)**2*h1
     normH1_2    = normH1_2 + ((grad1_node_val-grad1ref)**2)*h1
  end do

  ! -------> field visualization
  !call dirichlet_anal%write_to_file(0)
  
  ! -------> delete field
  call dirichlet_anal%delete()
  
  
  
  ! ******************************************************************  
  ! ------------------ TEST NON ANALYTIC ------------------------------
! ******************************************************************
  
!!! --------------------------------------------------------------------------
  !   Test case periodic non analytic
  !----------------------------------------------------------------------------
  
  
  ! ----> allocation table for field
  allocate(point1(nc1 + 1))
  allocate(tab_values(nc1 + 1))
  do i=1,nc1 + 1
     point1(i)       = real(i-1,f64)*h1 + ETA1MIN 
     tab_values(i)  = test_function_per(point1(i))
  end do
  
  
  ! ----> initializatio of the interpolator for the field
  
  call initialize_ad1d_interpolator( &
       interp_1d, &
       NUM_CELLS1+1, &
       ETA1MIN, &
       ETA1MAX, &
       SLL_PERIODIC, &
       SLL_PERIODIC,&
       SPLINE_DEG1)
  
  ! ----> initialization of the field
  
  periodic_discrete => new_scalar_field_1d_discrete_alt( &
       "periodic_discrete", &
       interp_1d, &
       SLL_PERIODIC, &
       SLL_PERIODIC,&
       mesh_1d,&
       point1,&
       nc1+1)

  print*, 'pasesd'
  ! ------- > allocation values of field
  call periodic_discrete%set_field_data(tab_values)
  print*, 'pasesd'
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
  !print*, normH1_5
  ! -------> field visualization 
!!$  call periodic_discrete%write_to_file(0)


  ! -------> delete field
  call periodic_discrete%delete()

  ! -------> delete table
  DEALLOCATE(point1)
  DEALLOCATE(tab_values)
  
!!! --------------------------------------------------------------------------
  !   Test case dirichlet non analytic
  !----------------------------------------------------------------------------
  
  
   ! ----> allocation table for field
  allocate(point1(nc1 + 1))
  allocate(tab_values(nc1 + 1))
  do i=1,nc1 + 1
     point1(i)       = real(i-1,f64)*h1 + ETA1MIN 
     tab_values(i)   = test_function_dir(point1(i) )
  end do

  

  ! ----> initializatio of the interpolator for the field
  
  call initialize_ad1d_interpolator( &
       interp_1d, &
       NUM_CELLS1+1, &
       ETA1MIN, &
       ETA1MAX, &
       SLL_DIRICHLET,&
       SLL_DIRICHLET,&
       SPLINE_DEG1)
  
  ! ----> initialization of the field
  
  dirichlet_discrete => new_scalar_field_1d_discrete_alt( &
       "dirichlet_discrete", &
       interp_1d, &
       SLL_DIRICHLET,&
       SLL_DIRICHLET,&
       mesh_1d,&
       point1,&
       nc1+1)
  
  ! ------- > allocation values of field
  call dirichlet_discrete%set_field_data(tab_values)
  ! --------> Compute coefficients of the field
  call dirichlet_discrete%update_interpolation_coefficients( )
  

  ! -------> compute error norm L2 and H1
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
  ! -------> field visualization 
!!$  call dirichlet_discrete%write_to_file(0)


  ! -------> delete field
  call dirichlet_discrete%delete()

  ! -------> delete table
  DEALLOCATE(point1)
  DEALLOCATE(tab_values)

  

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
end program unit_test_1d_alternative

 ! ------------->FUNCTION PERIODIC
function test_function_per( eta1) result(res)
   use sll_constants
   real(8) :: res
   real(8), intent(in) :: eta1
  intrinsic :: cos
  !real(8), dimension(:), intent(in), optional :: params
  res = cos(2*sll_pi*eta1)
end function test_function_per

function test_function_per_der1( eta1) result(res)
  use sll_constants
  intrinsic :: cos,sin
  real(8) :: res
  real(8), intent(in) :: eta1
 
  res = -2*sll_pi*sin(2*sll_pi*eta1)
end function test_function_per_der1


! ------------->FUNCTION DIRICHLET
function test_function_dir( eta1) result(res)
  use sll_constants
  intrinsic :: cos,sin
  real(8) :: res
  real(8), intent(in) :: eta1
  
  
  res = sin(2*sll_pi*eta1)
end function test_function_dir

function test_function_dir_der1( eta1) result(res)
  use sll_constants
  intrinsic :: cos
  real(8) :: res
  real(8), intent(in) :: eta1
  res = 2.0*sll_pi*cos(2*sll_pi*eta1)
end function test_function_dir_der1
