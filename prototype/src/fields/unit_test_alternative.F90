program unit_test_alternative
#include "sll_working_precision.h"
#include "sll_memory.h"
  use sll_logical_meshes
  use sll_constants
  use sll_module_scalar_field_2d_alternative
  use sll_module_coordinate_transformations_2d
  use sll_common_coordinate_transformations
  use helper_functions
  implicit none

#define SPLINE_DEG1 3
#define SPLINE_DEG2 3
#define NUM_CELLS1  64
#define NUM_CELLS2  64
#define ETA1MIN  0.0_f64
#define ETA1MAX  1.0_f64
#define ETA2MIN  0.0_f64
#define ETA2MAX  1.0_f64
#define PRINT_COMPARISON .false.
  
  type(sll_logical_mesh_2d), pointer               :: mesh_2d
  class(sll_coordinate_transformation_2d_base), pointer :: T
  ! either of these type declarations can be used to work. Initialization is
  ! different.
  class(sll_scalar_field_2d_base), pointer         :: field_2d

  class(sll_scalar_field_2d_base), pointer         :: doubly_periodic_anal
  class(sll_scalar_field_2d_base), pointer         :: periodique_dirichlet_anal
  class(sll_scalar_field_2d_base), pointer         :: dirichlet_dirichlet_anal
  class(sll_scalar_field_2d_base), pointer      :: dirichlet_periodique_anal
  class(sll_scalar_field_2d_base), pointer      :: doubly_periodic_discrete
  class(sll_scalar_field_2d_base), pointer      :: periodique_dirichlet_discrete
  class(sll_scalar_field_2d_base), pointer      :: dirichlet_dirichlet_discrete
  class(sll_scalar_field_2d_base), pointer      :: dirichlet_periodique_discrete
  sll_int32 :: nc1, nc2, iplot
  sll_real64 :: grad1_node_val,grad2_node_val,grad1ref,grad2ref
  sll_real64, dimension(:,:), pointer :: tab_values
  type(arb_deg_2d_interpolator), target                 :: interp_2d
  sll_real64 :: node_val,ref
  ! procedure(polar_x1), pointer :: px1, px2, pjac11, pjac12, pjac21, pjac22
  ! type(init_landau_2d), target :: init_landau
  ! class(scalar_field_2d_initializer_base), pointer    :: pfinit
  ! type(cubic_spline_1d_interpolator), target  :: interp_eta1
  ! type(cubic_spline_1d_interpolator), target  :: interp_eta2
  ! class(sll_interpolator_1d_base), pointer :: interp_eta1_ptr
  ! class(sll_interpolator_1d_base), pointer :: interp_eta2_ptr
  type(arb_deg_2d_interpolator), target                 :: interp_2d_term_source
  class(sll_scalar_field_2d_base), pointer              :: rho
  sll_real64, dimension(:,:), allocatable    :: calculated
  sll_real64, dimension(:,:), allocatable    :: difference
  sll_real64, dimension(:,:), pointer    :: tab_rho
  sll_real64, dimension(:),   allocatable    :: point1
  sll_real64, dimension(:),   allocatable    :: point2
  sll_real64 :: eta1,eta2
  sll_real64  :: h1,h2
  sll_int32 :: npts1,npts2
  sll_int32 :: i,j
  sll_int32 :: ierr
!!$  real(8), external :: test_function_perper
!!$  real(8), external :: test_function_perper_der1
!!$  real(8), external :: test_function_perper_der2
!!$  real(8), external :: test_function_perdir
!!$  real(8), external :: test_function_perdir_der1
!!$  real(8), external :: test_function_perdir_der2
!!$  real(8), external :: test_function_dirper
!!$  real(8), external :: test_function_dirper_der1
!!$  real(8), external :: test_function_dirper_der2
!!$  real(8), external :: test_function_dirdir
!!$  real(8), external :: test_function_dirdir_der1
!!$  real(8), external :: test_function_dirdir_der2
  
  sll_real64 :: normL2_1,normL2_2,normL2_3,normL2_4
  sll_real64 :: normL2_5,normL2_6,normL2_7,normL2_8
  sll_real64 :: normH1_1,normH1_2,normH1_3,normH1_4
  sll_real64 :: normH1_5,normH1_6,normH1_7,normH1_8
  sll_real64, dimension(1) :: params_identity

  params_identity(:) = (/ 0.0_f64 /)

  ! logical mesh
  nc1 = NUM_CELLS1
  nc2 = NUM_CELLS1
  h1 = (ETA1MAX-ETA1MIN)/real(nc1,f64)
  h2 = (ETA2MAX-ETA2MIN)/real(nc2,f64)
  print *, 'h1 = ', h1
  print *, 'h2 = ', h2
  


  ! First thing, initialize the logical mesh associated with this problem.        
  mesh_2d => new_logical_mesh_2d( NUM_CELLS1, NUM_CELLS2, &
       ETA1MIN, ETA1MAX, ETA2MIN,ETA2MAX )
  
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
       identity_jac22, &
       params_identity )
  print *, 'initialized transformation'



  ! ******************************************************************
  ! ------------------ TEST ANALYTIC ------------------------------  
  ! ******************************************************************
  
!!! --------------------------------------------------------------------------
  !   Test case periodic-periodic analytic
  !----------------------------------------------------------------------------
  
  
  
!!$  call field_2d_a%initialize( &
!!$       test_function_perper, &
!!$       'doubly_periodic', &
!!$       T, &
!!$       SLL_PERIODIC, &
!!$       SLL_PERIODIC, &
!!$       SLL_PERIODIC, &
!!$       SLL_PERIODIC )
  
  ! ----> initialization of the field
  doubly_periodic_anal  => new_scalar_field_2d_analytic_alt( &
       test_function_perper, &
       "doubly_periodic_anal", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC,&
       (/0.0_f64/), & ! could be anything in this case
       first_deriv_eta1=test_function_perper_der1,&
       first_deriv_eta2=test_function_perper_der2)
  print *, 'initialized field 2d'
  
  ! -------> compute error norm L2 and H1
  normL2_1 = 0.0_f64
  normH1_1 = 0.0_f64
  do j=1,nc2+1
     do i=1,nc1+1
        eta1       = real(i-1,f64)*h1 + ETA1MIN
        eta2       = real(j-1,f64)*h2 + ETA2MIN
        node_val   = doubly_periodic_anal%value_at_point(eta1,eta2)
        grad1_node_val = doubly_periodic_anal%first_deriv_eta1_value_at_point(eta1, eta2)
        grad2_node_val = doubly_periodic_anal%first_deriv_eta2_value_at_point(eta1, eta2)
        ref        = test_function_perper(eta1,eta2,(/0.0_f64/))
        grad1ref   = test_function_perper_der1(eta1,eta2,(/0.0_f64/))
        grad2ref   = test_function_perper_der2(eta1,eta2,(/0.0_f64/))
        
        if(PRINT_COMPARISON) then
           print *, '(eta1,eta2) = ', eta1, eta2, 'calculated = ', node_val, &
                'theoretical = ', ref, 'difference=', node_val-ref
           print *, '(eta1,eta2) = ', eta1, eta2, 'calculated = ', grad1_node_val, &
                'theoretical = ', grad1ref, 'difference=',grad1ref-grad1_node_val
           
        end if
        
        normL2_1    = normL2_1 + (node_val-ref)**2*h1*h2
        normH1_1    = normH1_1 + ((grad1_node_val-grad1ref)**2+&
             (grad2_node_val-grad2ref)**2)*h1*h2
        
     end do
  end do
  !print*, normL2_1
  
  ! -------> field visualization 
  call doubly_periodic_anal%write_to_file(0)
  
  ! the following call can also be made as:
  ! call field_2d_a%delete()
  ! we leave this as follows to test if any compilers complain about this
  ! syntax.
  !call delete(field_2d_a)
  
  ! -------> delete field
  call doubly_periodic_anal%delete()
  
!!! --------------------------------------------------------------------------
  !   Test case periodic-dirichlet analytic
  !----------------------------------------------------------------------------
  
  ! ----> initialization of the field
  periodique_dirichlet_anal  => new_scalar_field_2d_analytic_alt( &
       test_function_perdir, &
       "periodique_dirichlet_anal", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET,&
       (/0.0_f64/), & ! could be anything
       first_deriv_eta1=test_function_perdir_der1,&
       first_deriv_eta2=test_function_perdir_der2) 
  print *, 'initialized field 2d'
  
  ! -------> compute error norm L2 and H1
  normL2_2 = 0.0_f64
  normH1_2 = 0.0_f64
  do j=1,nc2+1
     do i=1,nc1+1
        eta1       = real(i-1,f64)*h1 + ETA1MIN
        eta2       = real(j-1,f64)*h2 + ETA2MIN
        node_val   = periodique_dirichlet_anal%value_at_point(eta1,eta2)
        grad1_node_val = periodique_dirichlet_anal%first_deriv_eta1_value_at_point(&
                                                                              eta1, eta2)
        grad2_node_val = &
             periodique_dirichlet_anal%first_deriv_eta2_value_at_point(&
             eta1, eta2)
        ref        = test_function_perdir(eta1,eta2,(/0.0_f64/))
        grad1ref   = test_function_perdir_der1(eta1,eta2,(/0.0_f64/))
        grad2ref   = test_function_perdir_der2(eta1,eta2,(/0.0_f64/))
        
        if(PRINT_COMPARISON) then
           print *, '(eta1,eta2) = ', eta1, eta2, 'calculated = ', node_val, &
                'theoretical = ', ref, 'difference=', node_val-ref
           print *, '(eta1,eta2) = ', eta1, eta2, 'calculated = ', grad1_node_val, &
                'theoretical = ', grad1ref, 'difference=',grad1ref-grad1_node_val
           
        end if
        
        normL2_2    = normL2_2 + (node_val-ref)**2*h1*h2
        normH1_2    = normH1_2 + ((grad1_node_val-grad1ref)**2+&
             (grad2_node_val-grad2ref)**2)*h1*h2
        
     end do
  end do

  ! -------> field visualization
   call periodique_dirichlet_anal%write_to_file(0)
  
  ! the following call can also be made as:
  ! call field_2d_a%delete()
  ! we leave this as follows to test if any compilers complain about this
  ! syntax.
  !call delete(field_2d_a)
  
   ! -------> delete field
  call periodique_dirichlet_anal%delete()
  
!!! --------------------------------------------------------------------------
  !   Test case dirichlet-periodic analytic
  !----------------------------------------------------------------------------
  


  ! ----> initialization of the field
  dirichlet_periodique_anal  => new_scalar_field_2d_analytic_alt( &
       test_function_dirper, &
       "dirichlet_periodique_anal", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET,&
       SLL_PERIODIC, &
       SLL_PERIODIC,&
       (/0.0_f64/), &
       first_deriv_eta1=test_function_dirper_der1,&
       first_deriv_eta2=test_function_dirper_der2)

  print *, 'initialized field 2d'
  
  ! -------> compute error norm L2 and H1
  normL2_3 = 0.0_f64
  normH1_3 = 0.0_f64
  do j=1,nc2+1
     do i=1,nc1+1
        eta1       = real(i-1,f64)*h1 + ETA1MIN
        eta2       = real(j-1,f64)*h2 + ETA2MIN
        node_val   = dirichlet_periodique_anal%value_at_point(eta1,eta2)
        grad1_node_val = dirichlet_periodique_anal%first_deriv_eta1_value_at_point(&
                                                                              eta1, eta2)
        grad2_node_val = &
             dirichlet_periodique_anal%first_deriv_eta2_value_at_point(&
             eta1, eta2)
        ref        = test_function_dirper(eta1,eta2,(/0.0_f64/))
        grad1ref   = test_function_dirper_der1(eta1,eta2,(/0.0_f64/))
        grad2ref   = test_function_dirper_der2(eta1,eta2,(/0.0_f64/))
        
        if(PRINT_COMPARISON) then
           print *, '(eta1,eta2) = ', eta1, eta2, 'calculated = ', node_val, &
                'theoretical = ', ref, 'difference=', node_val-ref
           print *, '(eta1,eta2) = ', eta1, eta2, 'calculated = ', grad1_node_val, &
                'theoretical = ', grad1ref, 'difference=',grad1ref-grad1_node_val
           
        end if
        
        normL2_3    = normL2_3 + (node_val-ref)**2*h1*h2
        normH1_3    = normH1_3 + ((grad1_node_val-grad1ref)**2+&
             (grad2_node_val-grad2ref)**2)*h1*h2
        
     end do
  end do

  ! -------> field visualization 
   call dirichlet_periodique_anal%write_to_file(0)
  
  ! the following call can also be made as:
  ! call field_2d_a%delete()
  ! we leave this as follows to test if any compilers complain about this
  ! syntax.
  !call delete(field_2d_a)
  
   ! -------> delete field
  call dirichlet_periodique_anal%delete()


!!! --------------------------------------------------------------------------
  !   Test case dirichlet-dirichlet analytic
  !----------------------------------------------------------------------------
  
  ! ----> initialization of the field
  dirichlet_dirichlet_anal  => new_scalar_field_2d_analytic_alt( &
       test_function_dirdir, &
       "dirichlet_dirichlet_anal", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET,&
       SLL_DIRICHLET, &
       SLL_DIRICHLET,&
       (/0.0_f64/), &
       first_deriv_eta1= test_function_dirdir_der1,&
       first_deriv_eta2=test_function_dirdir_der2)
  
  print *, 'initialized field 2d'
  
  ! -------> compute error norm L2 and H1
  normL2_4 = 0.0_f64
  normH1_4 = 0.0_f64
  do j=1,nc2+1
     do i=1,nc1+1
        eta1       = real(i-1,f64)*h1 + ETA1MIN
        eta2       = real(j-1,f64)*h2 + ETA2MIN
        node_val   = dirichlet_dirichlet_anal%value_at_point(eta1,eta2)
        grad1_node_val = dirichlet_dirichlet_anal%first_deriv_eta1_value_at_point(&
                                                                              eta1, eta2)
        grad2_node_val = dirichlet_dirichlet_anal%first_deriv_eta2_value_at_point(&
                                                                              eta1, eta2)
        ref        = test_function_dirdir(eta1,eta2,(/0.0_f64/))
        grad1ref   = test_function_dirdir_der1(eta1,eta2,(/0.0_f64/))
        grad2ref   = test_function_dirdir_der2(eta1,eta2,(/0.0_f64/))
        
        if(PRINT_COMPARISON) then
           print *, '(eta1,eta2) = ', eta1, eta2, 'calculated = ', node_val, &
                'theoretical = ', ref, 'difference=', node_val-ref
           print *, '(eta1,eta2) = ', eta1, eta2, 'calculated = ', grad1_node_val, &
                'theoretical = ', grad1ref, 'difference=',grad1ref-grad1_node_val
           
        end if
        
        normL2_4    = normL2_4 + (node_val-ref)**2*h1*h2
        normH1_4    = normH1_4 + ((grad1_node_val-grad1ref)**2+&
             (grad2_node_val-grad2ref)**2)*h1*h2
        
     end do
  end do
  ! -------> field visualization 
   call dirichlet_dirichlet_anal%write_to_file(0)
  
  ! the following call can also be made as:
  ! call field_2d_a%delete()
  ! we leave this as follows to test if any compilers complain about this
  ! syntax.
  !call delete(field_2d_a)
  
   ! -------> delete field   
  call dirichlet_dirichlet_anal%delete()
  

  
  
! ******************************************************************  
  ! ------------------ TEST NON ANALYTIC ------------------------------
! ******************************************************************
  
!!! --------------------------------------------------------------------------
  !   Test case periodic-periodic non analytic
  !----------------------------------------------------------------------------
  

  ! ----> allocation table for field
  allocate(point1(nc1 + 1))
  allocate(point2(nc2 + 1))
  allocate(tab_values(nc1 + 1,nc2 + 1))
  do j=1,nc2 + 1
     do i=1,nc1 + 1
        point1(i)       = real(i-1,f64)*h1 + ETA1MIN 
        point2(j)       = real(j-1,f64)*h2 + ETA2MIN 
        tab_values(i,j)  = test_function_perper(point1(i),point2(j),(/0.0_f64/))
     end do
  end do
  

  ! ----> initializatio of the interpolator for the field
  
  call initialize_ad2d_interpolator( &
       interp_2d, &
       NUM_CELLS1+1, &
       NUM_CELLS2+1, &
       ETA1MIN, &
       ETA1MAX, &
       ETA2MIN, &
       ETA2MAX, &
       SLL_PERIODIC, &
       SLL_PERIODIC,&
       SLL_PERIODIC,&
       SLL_PERIODIC,&
       SPLINE_DEG1, &
       SPLINE_DEG2)
  
  ! ----> initialization of the field
  
  doubly_periodic_discrete => new_scalar_field_2d_discrete_alt( &
       "doubly_periodic_discrete", &
       interp_2d, &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC,&
       SLL_PERIODIC,&
       SLL_PERIODIC,&
       point1,&
       nc1+1,&
       point2,&
       nc2+1)
  
  ! ------- > allocation values of field
  call doubly_periodic_discrete%set_field_data(tab_values)
  ! --------> Compute coefficients of the field
  call doubly_periodic_discrete%update_interpolation_coefficients( )
  

  ! -------> compute error norm L2 and H1
  normL2_5 = 0.0_f64
  normH1_5 = 0.0_f64
  do j=1,nc2 + 1
     do i=1,nc1 + 1 
        eta1 = real(i-1,f64)*h1 + ETA1MIN 
        eta2 = real(j-1,f64)*h2 + ETA2MIN
        node_val   = doubly_periodic_discrete%value_at_point(eta1,eta2)
        grad1_node_val = doubly_periodic_discrete%first_deriv_eta1_value_at_point(&
                                                                              eta1, eta2)
        grad2_node_val = doubly_periodic_discrete%first_deriv_eta2_value_at_point(&
                                                                              eta1, eta2)
        ref        = test_function_perper(eta1,eta2,(/0.0_f64/))
        grad1ref   = test_function_perper_der1(eta1,eta2,(/0.0_f64/))
        grad2ref   = test_function_perper_der2(eta1,eta2,(/0.0_f64/))
        if(PRINT_COMPARISON) then
           print *, '(eta1,eta2) = ', eta1, eta2, 'calculated = ', node_val, &
                'theoretical = ', ref, 'difference=', node_val-ref
           print *, '(eta1,eta2) = ', eta1, eta2, 'calculated = ', grad1_node_val, &
                'theoretical = ', grad1ref, 'difference=',grad1ref-grad1_node_val
           print *, '(eta1,eta2) = ', eta1, eta2, 'calculated = ', grad2_node_val, &
                'theoretical = ', grad2ref, 'difference=',grad2ref-grad2_node_val
        end if
        
        normL2_5    = normL2_5 + (node_val-ref)**2*h1*h2
        normH1_5    = normH1_5 + ((grad1_node_val-grad1ref)**2+&
             (grad2_node_val-grad2ref)**2)*h1*h2
        
     end do
  end do
  !print*, normH1_5
  ! -------> field visualization 
  call doubly_periodic_discrete%write_to_file(0)


  ! -------> delete field
  call doubly_periodic_discrete%delete()

  ! -------> delete table
  DEALLOCATE(point1)
  DEALLOCATE(point2)
  DEALLOCATE(tab_values)
  
!!! --------------------------------------------------------------------------
  !   Test case periodic-dirichlet non analytic
!----------------------------------------------------------------------------


   ! ----> allocation table for field
  allocate(point1(nc1 + 1))
  allocate(point2(nc2 + 1))
  allocate(tab_values(nc1 + 1,nc2 + 1))
  do j=1,nc2 + 1
     do i=1,nc1 + 1
        point1(i)       = real(i-1,f64)*h1 + ETA1MIN 
        point2(j)       = real(j-1,f64)*h2 + ETA2MIN 
        tab_values(i,j)  = test_function_perdir(point1(i),point2(j),(/0.0_f64/))
     end do
  end do
  

  ! ----> initializatio of the interpolator for the field
  
  call initialize_ad2d_interpolator( &
       interp_2d, &
       NUM_CELLS1+1, &
       NUM_CELLS2+1, &
       ETA1MIN, &
       ETA1MAX, &
       ETA2MIN, &
       ETA2MAX, &
       SLL_PERIODIC, &
       SLL_PERIODIC,&
       SLL_DIRICHLET,&
       SLL_DIRICHLET,&
       SPLINE_DEG1, &
       SPLINE_DEG2)
  
  ! ----> initialization of the field
  
  periodique_dirichlet_discrete => new_scalar_field_2d_discrete_alt( &
       "periodique_dirichlet_discrete", &
       interp_2d, &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC,&
       SLL_DIRICHLET,&
       SLL_DIRICHLET,&
       point1,&
       nc1+1,&
       point2,&
       nc2+1)
  
  ! ------- > allocation values of field
  call periodique_dirichlet_discrete%set_field_data(tab_values)
  ! --------> Compute coefficients of the field
  call periodique_dirichlet_discrete%update_interpolation_coefficients( )
  

  ! -------> compute error norm L2 and H1
  normL2_6 = 0.0_f64
  normH1_6 = 0.0_f64
  do j=1,nc2 + 1
     do i=1,nc1 + 1 
        eta1 = real(i-1,f64)*h1 + ETA1MIN 
        eta2 = real(j-1,f64)*h2 + ETA2MIN
        node_val   = periodique_dirichlet_discrete%value_at_point(eta1,eta2)
        grad1_node_val = &
             periodique_dirichlet_discrete%first_deriv_eta1_value_at_point(&
             eta1, eta2)
        grad2_node_val = &
             periodique_dirichlet_discrete%first_deriv_eta2_value_at_point(&
             eta1, eta2)
        ref        = test_function_perdir(eta1,eta2,(/0.0_f64/))
        grad1ref   = test_function_perdir_der1(eta1,eta2,(/0.0_f64/))
        grad2ref   = test_function_perdir_der2(eta1,eta2,(/0.0_f64/))
        if(PRINT_COMPARISON) then
           print *, '(eta1,eta2) = ', eta1, eta2, 'calculated = ', node_val, &
                'theoretical = ', ref, 'difference=', node_val-ref
           print *, '(eta1,eta2) = ', eta1, eta2, 'calculated = ', grad1_node_val, &
                'theoretical = ', grad1ref, 'difference=',grad1ref-grad1_node_val
           
        end if
        
        normL2_6    = normL2_6 + (node_val-ref)**2*h1*h2
        normH1_6    = normH1_6 + ((grad1_node_val-grad1ref)**2+&
             (grad2_node_val-grad2ref)**2)*h1*h2
        
     end do
  end do
  
  ! -------> field visualization 
  call periodique_dirichlet_discrete%write_to_file(0)


  ! -------> delete field
  call periodique_dirichlet_discrete%delete()

  ! -------> delete table
  DEALLOCATE(point1)
  DEALLOCATE(point2)
  DEALLOCATE(tab_values)




!!! --------------------------------------------------------------------------
  !   Test case dirichlet-periodic non analytic
  !----------------------------------------------------------------------------



  

   ! ----> allocation table for field
  allocate(point1(nc1 + 1))
  allocate(point2(nc2 + 1))
  allocate(tab_values(nc1 + 1,nc2 + 1))
  do j=1,nc2 + 1
     do i=1,nc1 + 1
        point1(i)       = real(i-1,f64)*h1 + ETA1MIN 
        point2(j)       = real(j-1,f64)*h2 + ETA2MIN 
        tab_values(i,j)  = test_function_dirper(point1(i),point2(j),(/0.0_f64/))
     end do
  end do
  

  ! ----> initializatio of the interpolator for the field
  
  call initialize_ad2d_interpolator( &
       interp_2d, &
       NUM_CELLS1+1, &
       NUM_CELLS2+1, &
       ETA1MIN, &
       ETA1MAX, &
       ETA2MIN, &
       ETA2MAX, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET,&
       SLL_PERIODIC,&
       SLL_PERIODIC,&
       SPLINE_DEG1, &
       SPLINE_DEG2)
  
  ! ----> initialization of the field
  
  dirichlet_periodique_discrete => new_scalar_field_2d_discrete_alt( &
       "dirichlet_periodique_discrete", &
       interp_2d, &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET,&
       SLL_PERIODIC,&
       SLL_PERIODIC,&
       point1,&
       nc1+1,&
       point2,&
       nc2+1)
  
  ! ------- > allocation values of field
  call dirichlet_periodique_discrete%set_field_data(tab_values)
  ! --------> Compute coefficients of the field
  call dirichlet_periodique_discrete%update_interpolation_coefficients( )
  

  ! -------> compute error norm L2 and H1
  normL2_7 = 0.0_f64
  normH1_7 = 0.0_f64
  do j=1,nc2 + 1
     do i=1,nc1 + 1 
        eta1 = real(i-1,f64)*h1 + ETA1MIN 
        eta2 = real(j-1,f64)*h2 + ETA2MIN
        node_val   = dirichlet_periodique_discrete%value_at_point(eta1,eta2)
        grad1_node_val = dirichlet_periodique_discrete%first_deriv_eta1_value_at_point(&
                                                                              eta1, eta2)
        grad2_node_val = dirichlet_periodique_discrete%first_deriv_eta2_value_at_point(&
                                                                              eta1, eta2)
        ref        = test_function_dirper(eta1,eta2,(/0.0_f64/))
        grad1ref   = test_function_dirper_der1(eta1,eta2,(/0.0_f64/))
        grad2ref   = test_function_dirper_der2(eta1,eta2,(/0.0_f64/))
        if(PRINT_COMPARISON) then
           print *, '(eta1,eta2) = ', eta1, eta2, 'calculated = ', node_val, &
                'theoretical = ', ref, 'difference=', node_val-ref
           print *, '(eta1,eta2) = ', eta1, eta2, 'calculated = ', grad1_node_val, &
                'theoretical = ', grad1ref, 'difference=',grad1ref-grad1_node_val
           
        end if
        
        normL2_7    = normL2_7 + (node_val-ref)**2*h1*h2
        normH1_7    = normH1_7 + ((grad1_node_val-grad1ref)**2+&
             (grad2_node_val-grad2ref)**2)*h1*h2
        
     end do
  end do
  
  ! -------> field visualization 
  call dirichlet_periodique_discrete%write_to_file(0)


  ! -------> delete field
  call dirichlet_periodique_discrete%delete()

  ! -------> delete table
  DEALLOCATE(point1)
  DEALLOCATE(point2)
  DEALLOCATE(tab_values)



!!! --------------------------------------------------------------------------
!   Test case dirichlet-dirichlet non analytic
  !----------------------------------------------------------------------------

  
  

   ! ----> allocation table for field
  allocate(point1(nc1 + 1))
  allocate(point2(nc2 + 1))
  allocate(tab_values(nc1 + 1,nc2 + 1))
  do j=1,nc2 + 1
     do i=1,nc1 + 1
        point1(i)       = real(i-1,f64)*h1 + ETA1MIN 
        point2(j)       = real(j-1,f64)*h2 + ETA2MIN 
        tab_values(i,j)  = test_function_dirdir(point1(i),point2(j),(/0.0_f64/))
     end do
  end do
  

  ! ----> initializatio of the interpolator for the field
  
  call initialize_ad2d_interpolator( &
       interp_2d, &
       NUM_CELLS1+1, &
       NUM_CELLS2+1, &
       ETA1MIN, &
       ETA1MAX, &
       ETA2MIN, &
       ETA2MAX, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET,&
       SLL_DIRICHLET,&
       SLL_DIRICHLET,&
       SPLINE_DEG1, &
       SPLINE_DEG2)
  
  ! ----> initialization of the field
  
  dirichlet_dirichlet_discrete => new_scalar_field_2d_discrete_alt( &
       "dirichlet_dirichlet_discrete", &
       interp_2d, &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET,&
       SLL_DIRICHLET,&
       SLL_DIRICHLET,&
       point1,&
       nc1+1,&
       point2,&
       nc2+1)
  
  ! ------- > allocation values of field
  call dirichlet_dirichlet_discrete%set_field_data(tab_values)
  ! --------> Compute coefficients of the field
  call dirichlet_dirichlet_discrete%update_interpolation_coefficients( )
  

  ! -------> compute error norm L2 and H1
  normL2_8 = 0.0_f64
  normH1_8 = 0.0_f64
  do j=1,nc2 + 1
     do i=1,nc1 + 1 
        eta1 = real(i-1,f64)*h1 + ETA1MIN 
        eta2 = real(j-1,f64)*h2 + ETA2MIN
        node_val   = dirichlet_dirichlet_discrete%value_at_point(eta1,eta2)
        grad1_node_val = &
             dirichlet_dirichlet_discrete%first_deriv_eta1_value_at_point(&
             eta1, eta2)
        grad2_node_val = &
             dirichlet_dirichlet_discrete%first_deriv_eta2_value_at_point(&
             eta1, eta2)
        ref        = test_function_dirdir(eta1,eta2,(/0.0_f64/))
        grad1ref   = test_function_dirdir_der1(eta1,eta2,(/0.0_f64/))
        grad2ref   = test_function_dirdir_der2(eta1,eta2,(/0.0_f64/))
        if(PRINT_COMPARISON) then
           print *, '(eta1,eta2) = ', eta1, eta2, 'calculated = ', node_val, &
                'theoretical = ', ref, 'difference=', node_val-ref
           print *, '(eta1,eta2) = ', eta1, eta2, 'calculated = ', grad1_node_val, &
                'theoretical = ', grad1ref, 'difference=',grad1ref-grad1_node_val
           
        end if
        
        normL2_8    = normL2_8 + (node_val-ref)**2*h1*h2
        normH1_8    = normH1_8 + ((grad1_node_val-grad1ref)**2+&
             (grad2_node_val-grad2ref)**2)*h1*h2
        

!!$   do j=0,npts2-2
!!$     do i=0,npts1-2
!!$        eta1 = real(i,f64)*(ETA1MAX-ETA1MIN)/(2*(npts1-1)) + ETA1MIN 
!!$        eta2 = real(j,f64)*(ETA2MAX-ETA2MIN)/(2*(npts2-1)) + ETA2MIN
!!$       calculated(i+1,j+1) = rho%value_at_point(eta1,eta2)
!!$       difference(i+1,j+1) = calculated(i+1,j+1)-cos(2.0_f64*sll_pi*eta2)*cos(2.0_f64*sll_pi*eta1)
!       print*, 'point=',eta1,eta2,'difference=', difference(i+1,j+1), calculated(i+1,j+1),cos(2.0_f64*sll_pi*eta2)*cos(2.0_f64*sll_pi*eta1)

     end do
  end do
  
  ! -------> field visualization 
  call dirichlet_dirichlet_discrete%write_to_file(0)


  ! -------> delete field
  call dirichlet_dirichlet_discrete%delete()

  ! -------> delete table
  DEALLOCATE(point1)
  DEALLOCATE(point2)
  DEALLOCATE(tab_values)




  ! **********************  TESTS **************************************
  
  print*, '-------------------------------------------------------'
  print*, ' PERIODIC-PERIODIC ANALYTIC' 
  print*, '-------------------------------------------------------'
  print *,'Norm L2',sqrt(normL2_1),'Norm H1',sqrt(normH1_1),&
       h1**(SPLINE_DEG1),  h1**(SPLINE_DEG1-1)
  
  print*, '-------------------------------------------------------'
  print*, ' PERIODIC-DIRICHLET ANALYTIC' 
  print*, '-------------------------------------------------------'
  print *,'Norm L2',sqrt(normL2_2),'Norm H1',sqrt(normH1_2),&
       h1**(SPLINE_DEG1),  h1**(SPLINE_DEG1-1)
  
  print*, '-------------------------------------------------------'
  print*, ' DIRICHLET-PERIODIC ANALYTIC' 
  print*, '-------------------------------------------------------'
  print *,'Norm L2',sqrt(normL2_3),'Norm H1',sqrt(normH1_3),&
       h1**(SPLINE_DEG1),  h1**(SPLINE_DEG1-1)
  
  print*, '-------------------------------------------------------'
  print*, ' DIRICHLET-DIRCHLET ANALYTIC' 
  print*, '-------------------------------------------------------'
  print *,'Norm L2',sqrt(normL2_4),'Norm H1',sqrt(normH1_4),&
       h1**(SPLINE_DEG1),  h1**(SPLINE_DEG1-1)

  
  print*, '-------------------------------------------------------'
  print*, ' PERIODIC-PERIODIC NON ANALYTIC' 
  print*, '-------------------------------------------------------'
  print *,'Norm L2',sqrt(normL2_5),'Norm H1',sqrt(normH1_5),&
       h1**(SPLINE_DEG1),  h1**(SPLINE_DEG1-1)
  
  print*, '-------------------------------------------------------'
  print*, ' PERIODIC-DIRICHLET NON ANALYTIC' 
  print*, '-------------------------------------------------------'
  print *,'Norm L2',sqrt(normL2_6),'Norm H1',sqrt(normH1_6),&
       h1**(SPLINE_DEG1),  h1**(SPLINE_DEG1-1)

  print*, '-------------------------------------------------------'
  print*, ' DIRICHLET-PERIODIC NON ANALYTIC' 
  print*, '-------------------------------------------------------'
  print *,'Norm L2',sqrt(normL2_7),'Norm H1',sqrt(normH1_7),&
       h1**(SPLINE_DEG1),  h1**(SPLINE_DEG1-1)
  
  print*, '-------------------------------------------------------'
  print*, ' DIRICHLET-DIRCHLET NON ANALYTIC' 
  print*, '-------------------------------------------------------'
  print *,'Norm L2',sqrt(normL2_8),'Norm H1',sqrt(normH1_8),&
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
end program unit_test_alternative

