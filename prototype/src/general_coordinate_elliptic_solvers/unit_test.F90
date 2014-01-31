program test_general_elliptic_solver
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_logical_meshes
  use sll_module_coordinate_transformations_2d
  use sll_common_coordinate_transformations
  use sll_module_scalar_field_2d_alternative
  use sll_constants
  use sll_arbitrary_degree_spline_interpolator_2d_module
  use sll_timer


#ifdef _UMFPACK
  use sll_general_coordinate_elliptic_solver_module_umfpack
#else
  use sll_general_coordinate_elliptic_solver_module
#endif

  implicit none

#define SPLINE_DEG1 3
#define SPLINE_DEG2 3
#define NUM_CELLS1  10
#define NUM_CELLS2  10
#define ETA1MIN  0.0_f64
#define ETA1MAX  1.0_f64
#define ETA2MIN  0.0_f64
#define ETA2MAX  1.0_f64
#define PRINT_COMPARISON .false.

  type(sll_logical_mesh_2d), pointer                    :: mesh_2d
  class(sll_coordinate_transformation_2d_base), pointer :: T
  type(general_coordinate_elliptic_solver)              :: es
  type(arb_deg_2d_interpolator), target                 :: interp_2d
  type(arb_deg_2d_interpolator), target                 :: interp_2d_term_source
 ! class(sll_interpolator_2d_base), pointer              :: interp_2d_ptr
  class(sll_interpolator_2d_base), pointer              :: terme_source_interp
  class(sll_scalar_field_2d_base), pointer              :: a11_field_mat
  class(sll_scalar_field_2d_base), pointer              :: a12_field_mat
  class(sll_scalar_field_2d_base), pointer              :: a21_field_mat
  class(sll_scalar_field_2d_base), pointer              :: a22_field_mat
  class(sll_scalar_field_2d_base), pointer              :: b1_field_vect
  class(sll_scalar_field_2d_base), pointer              :: b2_field_vect
  class(sll_scalar_field_2d_base), pointer              :: c_field
  class(sll_scalar_field_2d_base), pointer              :: rho
  type(sll_scalar_field_2d_discrete_alt), pointer       :: phi
  type(sll_time_mark) :: t_reference
  sll_real64 :: t1i, t1e, t2i, t2e, t3i, t3e, t4i, t4e, t5i, t5e, t6i, t6e, &
       t7i, t7e, t8i, t8e, t9i,t9e,t10i,t10e,t11i,t11e,t12i,t12e,t95e,t95i
  ! sll_real64 :: t105e,t105i,t115e,t115i,t125i,t125e
  real(8), external :: func_zero
  real(8), external :: func_one
  real(8), external :: func_epsi
  real(8), external :: source_term_perper
  real(8), external :: source_term_perdir
  real(8), external :: source_term_dirper
  real(8), external :: source_term_chgt_perper
  real(8), external :: source_term_chgt_perdir
  real(8), external :: source_term_chgt_dirper
  real(8), external :: source_term_chgt_dirdir
  sll_real64, dimension(:,:), pointer :: values
  sll_real64 :: acc1,acc2,acc3,acc4,acc5,acc6,acc7,acc8,acc9
  sll_real64 :: acc10,acc11,acc12,acc95
  sll_real64 :: normL2_1,normL2_2,normL2_3,normL2_4,normL2_5,normL2_6
  sll_real64 :: normL2_7,normL2_8,normL2_9,normL2_95,normL2_10,normL2_11,normL2_12
  sll_real64 :: normH1_1,normH1_2,normH1_3,normH1_4,normH1_5,normH1_6
  sll_real64 :: normH1_7,normH1_8,normH1_9,normH1_95,normH1_10,normH1_11,normH1_12
 ! sll_real64 :: acc105,acc115,acc125
  sll_real64, dimension(:,:), allocatable    :: calculated
  sll_real64, dimension(:,:), allocatable    :: difference
  sll_real64, dimension(:,:), allocatable    :: reference
  sll_real64, dimension(:,:), allocatable    :: tab_rho
  sll_real64, dimension(:),   allocatable    :: point1
  sll_real64, dimension(:),   allocatable    :: point2
!  sll_real64, dimension(:,:), pointer        :: test_coeff
  sll_real64 :: val_jac
  sll_int32 :: ierr
  sll_int32  :: i, j
  sll_real64 :: h1,h2,eta1,eta2,node_val,ref
  sll_int32 :: npts1,npts2
  real(8), external :: sol_exacte_perper
  real(8), external :: sol_exacte_perper_der1
  real(8), external :: sol_exacte_perper_der2
  real(8), external :: sol_exacte_perdir  
  real(8), external :: sol_exacte_perdir_der1
  real(8), external :: sol_exacte_perdir_der2
  real(8), external :: sol_exacte_dirper
  real(8), external :: sol_exacte_dirper_der1
  real(8), external :: sol_exacte_dirper_der2
  real(8), external :: sol_exacte_chgt_perper
  real(8), external :: sol_exacte_chgt_perper_der1
  real(8), external :: sol_exacte_chgt_perper_der2
  real(8), external :: sol_exacte_chgt_perdir  
  real(8), external :: sol_exacte_chgt_perdir_der1
  real(8), external :: sol_exacte_chgt_perdir_der2
  real(8), external :: sol_exacte_chgt_dirper
  real(8), external :: sol_exacte_chgt_dirper_der1
  real(8), external :: sol_exacte_chgt_dirper_der2
  real(8), external :: sol_exacte_chgt_dirdir
  real(8), external :: sol_exacte_chgt_dirdir_der1
  real(8), external :: sol_exacte_chgt_dirdir_der2
  real(8), external :: adimension_chgt_x
  real(8), external :: adimension_chgt_y
  real(8), external :: jac11_adimension_chgt
  real(8), external :: jac12_adimension_chgt
  real(8), external :: jac21_adimension_chgt
  real(8), external :: jac22_adimension_chgt
  real(8), external :: sol_exacte_chgt_adim
  real(8), external :: source_term_chgt_adim
  real(8) :: integrale_solution
  real(8) :: integrale_solution_exacte
 ! sll_real64 :: node_val1
  !sll_real64 :: epsi
  !sll_real64 :: epsi1
  sll_real64 :: grad1_node_val,grad2_node_val,grad1ref,grad2ref
 ! epsi  =  0.000_f64
 ! epsi1 =  0.000_f64 ! penalization method
  sll_real64, dimension(1) :: whatever  ! dummy params array
  
  
!!$  !*******************************************************************
!!$  !        WHITHOUT CHANGE OF COORDINATES AND ANALYTIC DATA
!!$  !*******************************************************************
!!$  !--------------------------------------------------------------------
!!$  !     1 test case without change of coordinates 
!!$  !      periodic-periodic boundary conditions
!!$  !--------------------------------------------------------------------
!!$  
  print*, "-------------------------------------------------------------"
  print*, "1 test case witout change of coordinates"
  print*, "-------------------------------------------------------------"
  npts1 =  NUM_CELLS1 + 1
  npts2 =  NUM_CELLS2 + 1
  h1 = (ETA1MAX-ETA1MIN)/real(NPTS1-1,f64)
  h2 = (ETA2MAX-ETA2MIN)/real(NPTS2-1,f64)
  print *, 'h1 = ', h1
  print *, 'h2 = ', h2

  ! Table to represent the node values of phi
  SLL_ALLOCATE(values(NUM_CELLS1+1,NUM_CELLS2+1),ierr)
  SLL_ALLOCATE(calculated(npts1,npts2),ierr)
  SLL_ALLOCATE(difference(npts1,npts2),ierr)
  SLL_ALLOCATE(reference(npts1,npts2),ierr)
  values(:,:) = 0.0_f64

  ! First thing, initialize the logical mesh associated with this problem.        
  mesh_2d => new_logical_mesh_2d( NUM_CELLS1, NUM_CELLS2, &
       ETA1MIN, ETA1MAX, ETA2MIN,ETA2MAX )

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
       identity_jac22, &
       (/ 0.0_f64 /) )
  print *, 'initialized coordinate transformation'

  ! Thirdly, each field object must be initialized using the same logical
  ! mesh and coordinate transformation.
  a11_field_mat => new_scalar_field_2d_analytic_alt( &
       func_one, &
       "a11", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       whatever ) 

 
  a12_field_mat => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "a12", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       whatever ) 

  a21_field_mat => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "a21", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       whatever  ) 

  a22_field_mat => new_scalar_field_2d_analytic_alt( &
       func_one, &
       "a22", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       whatever  ) 


  b1_field_vect => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "b1", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       whatever, & 
       first_deriv_eta1 = func_zero, &
       first_deriv_eta2 = func_zero) 

  b2_field_vect => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "b2", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       whatever, &
       first_deriv_eta1 = func_zero, &
       first_deriv_eta2 = func_zero)


  c_field => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "c_field", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       whatever  )

  rho => new_scalar_field_2d_analytic_alt( &
       source_term_perper, &
       "rho1", &     
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       whatever  )

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

 ! interp_2d_ptr => interp_2d

  phi => new_scalar_field_2d_discrete_alt( &
       "phi1", &
       interp_2d, &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC )
  
  call phi%set_field_data( values )
  call phi%update_interpolation_coefficients( )

  print *, 'initialized fields...'

  call sll_set_time_mark(t_reference)

  call initialize_general_elliptic_solver( &
       es, &
       SPLINE_DEG1, &
       SPLINE_DEG2, &
       NUM_CELLS1, &
       NUM_CELLS2, &
       ES_GAUSS_LEGENDRE, &
       ES_GAUSS_LEGENDRE, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       ETA1MIN, &
       ETA1MAX, &
       ETA2MIN, &
       ETA2MAX)
 
  t1i = sll_time_elapsed_since(t_reference)
 
  print *, 'Initialized ES object'

  call sll_set_time_mark(t_reference)

  ! compute matrix the field
  call factorize_mat_es(&
       es, &
       a11_field_mat, &
       a12_field_mat,&
       a21_field_mat,&
       a22_field_mat,&
       b1_field_vect,&
       b2_field_vect,&
       c_field)!, &
      ! rho)

  ! solve the field
  call solve_general_coordinates_elliptic_eq(&
       es,&
       rho,&
       phi)
  
  t1e = sll_time_elapsed_since(t_reference)

  !print *, 'Completed solution',es%phi_vec
!  print*, 'reorganizaton of splines coefficients of solution'

    print *, 'Compare the values of the transformation at the nodes: '
!!$  
    
  acc1 = 0.0_f64
  normL2_1 = 0.0_f64
  normH1_1 = 0.0_f64
  integrale_solution = 0.0_f64
  integrale_solution_exacte = 0.0_f64
  do j=0,npts2-1
     do i=0,npts1-1
        eta1       = real(i,f64)*h1 + ETA1MIN
        eta2       = real(j,f64)*h2 + ETA2MIN
    
        node_val   = phi%value_at_point(eta1,eta2)

        grad1_node_val = phi%first_deriv_eta1_value_at_point(eta1, eta2)
        grad2_node_val = phi%first_deriv_eta2_value_at_point(eta1, eta2)
        ref        = sol_exacte_perper(eta1,eta2)
        grad1ref   = sol_exacte_perper_der1(eta1,eta2)
        grad2ref   = sol_exacte_perper_der2(eta1,eta2)
        calculated(i+1,j+1) = node_val
        difference(i+1,j+1) = ref-node_val
        reference( i+1,j+1) = ref
        if(PRINT_COMPARISON) then
           print *, '(eta1,eta2) = ', eta1, eta2, 'calculated = ', node_val, &
                'theoretical = ', ref, 'difference=', node_val-ref
           print *, '(eta1,eta2) = ', eta1, eta2, 'calculated = ', grad1_node_val, &
                'theoretical = ', grad1ref, 'difference=',grad1ref-grad1_node_val
           print *, '(eta1,eta2) = ', eta1, eta2, 'calculated = ', grad2_node_val, &
                'theoretical = ', grad2ref, 'difference=',grad2ref-grad2_node_val
           
        end if
        acc1        = acc1 + abs(node_val-ref)
        normL2_1    = normL2_1 + (node_val-ref)**2*h1*h2
        normH1_1    = normH1_1 + ((grad1_node_val-grad1ref)**2+(grad2_node_val-grad2ref)**2)*h1*h2
        integrale_solution = integrale_solution + ref*h1*h2!node_val
     end do
  end do

  print*, 'integrale de la solution =', sum(calculated(1:npts1-1,1:npts2-1))*h1*h2,&
       'integrale de la solution exacte=',sum(reference(1:npts1-1,1:npts2-1))*h1*h2
  call phi%write_to_file(0)
  
 
  ! delete things...
  call delete(es)
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
  
  SLL_DEALLOCATE(values, ierr)
  SLL_DEALLOCATE_ARRAY(calculated,ierr)
  SLL_DEALLOCATE_ARRAY(difference,ierr)
  SLL_DEALLOCATE_ARRAY(reference,ierr)
  
!!$  !--------------------------------------------------------------------
!!$  
!!$  !     2 test case without chane of coordinates 
!!$  !      periodic-dirichlet boundary conditions
!!$  
!!$  !--------------------------------------------------------------------
!!$  
!!$  
!!$  
  print*, "-------------------------------------------------------------"
  print*, " 2 test case witout change of coordinates"
  print*, " periodic-dirichlet boundary conditions"
  print*, "-------------------------------------------------------------"
  npts1 =  NUM_CELLS1 + 1
  npts2 =  NUM_CELLS2 + 1
  h1 = (ETA1MAX-ETA1MIN)/real(NPTS1-1,f64)
  h2 = (ETA2MAX-ETA2MIN)/real(NPTS2-1,f64)
  print *, 'h1 = ', h1
  print *, 'h2 = ', h2
  ! Table to represent the node values of phi
  SLL_ALLOCATE(values(NUM_CELLS1+1,NUM_CELLS2+1),ierr)
  SLL_ALLOCATE(calculated(npts1,npts2),ierr)
  SLL_ALLOCATE(difference(npts1,npts2),ierr)
  SLL_ALLOCATE(reference(npts1,npts2),ierr)
  values(:,:) = 0.0_f64
  
  ! First thing, initialize the logical mesh associated with this problem.        
  mesh_2d => new_logical_mesh_2d( NUM_CELLS1, NUM_CELLS2, &
       ETA1MIN, ETA1MAX, ETA2MIN,ETA2MAX )
  
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
       identity_jac22, &
       (/0.0_f64/) )
  print *, 'initialized coordinate transformation'

  ! Thirdly, each field object must be initialized using the same logical
  ! mesh and coordinate transformation.
  a11_field_mat => new_scalar_field_2d_analytic_alt( &
       func_one, &
       "a11", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       whatever  ) 

  a12_field_mat => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "a12", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       whatever  ) 

  a21_field_mat => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "a21", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       whatever  ) 
  
  a22_field_mat  => new_scalar_field_2d_analytic_alt( &
       func_one, &
       "a22", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       whatever  ) 
  

  b1_field_vect => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "b1", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       whatever, &
       first_deriv_eta1 = func_zero, &
       first_deriv_eta2 = func_zero)  

  b2_field_vect => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "b2", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       whatever, &
       first_deriv_eta1 = func_zero, &
       first_deriv_eta2 = func_zero)
  
  c_field => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "c_field", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       whatever )

  rho => new_scalar_field_2d_analytic_alt( &
       source_term_perdir, &
       "rho2", &     
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       whatever )

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
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SPLINE_DEG1, &
       SPLINE_DEG2 )

  
  phi => new_scalar_field_2d_discrete_alt( &
       "phi2", &
       interp_2d, &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET)
  call phi%set_field_data(values)
  call phi%update_interpolation_coefficients()
  
  print *, 'initialized fields...'
!  print *, 'a = ', es%csr_mat%opr_a

  call sll_set_time_mark(t_reference)

  call initialize_general_elliptic_solver( &
       es, &
       SPLINE_DEG1, &
       SPLINE_DEG2, &
       NUM_CELLS1, &
       NUM_CELLS2, &
       ES_GAUSS_LEGENDRE, &
       ES_GAUSS_LEGENDRE, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       ETA1MIN, &
       ETA1MAX, &
       ETA2MIN, &
       ETA2MAX)
  t2i = sll_time_elapsed_since(t_reference) 
  print *, 'Initialized ES object'
  
  call sll_set_time_mark(t_reference)

 call factorize_mat_es(&
       es, &
       a11_field_mat, &
       a12_field_mat,&
       a21_field_mat,&
       a22_field_mat,&
       b1_field_vect,&
       b2_field_vect,&
       c_field)!, &
       !rho)
  
!!$
  call solve_general_coordinates_elliptic_eq(&
       es,&
       rho,&
       phi)
  
  !print *, 'Completed solution',es%phi_vec
  
  t2e = sll_time_elapsed_since(t_reference)
  
!  print *, 'Compare the values of the transformation at the nodes: '
  
  acc2 = 0.0_f64
  normL2_2 = 0.0_f64
  normH1_2 = 0.0_f64
  do j=0,npts2-1
     do i=0,npts1-1
        eta1       = real(i,f64)*h1 + ETA1MIN
        eta2       = real(j,f64)*h2 + ETA2MIN
        node_val   = phi%value_at_point(eta1,eta2)
        grad1_node_val = phi%first_deriv_eta1_value_at_point(eta1, eta2)
        grad2_node_val = phi%first_deriv_eta2_value_at_point(eta1, eta2)
        ref        = sol_exacte_perdir(eta1,eta2)
        grad1ref   = sol_exacte_perdir_der1(eta1,eta2)
        grad2ref   = sol_exacte_perdir_der2(eta1,eta2)
        !        print*,sin(2*sll_pi*eta1)*cos(2*sll_pi*eta1)
        calculated(i+1,j+1) = node_val
        difference(i+1,j+1) = ref-node_val
        reference( i+1,j+1) = ref
        if(PRINT_COMPARISON) then
           print *, '(eta1,eta2) = ', eta1, eta2, 'calculated = ', node_val, &
                'theoretical = ', ref,'difference=',ref-node_val
        end if
        acc2        = acc2 + abs(node_val-ref)
        normL2_2    = normL2_2 + (node_val-ref)**2*h1*h2
        normH1_2    = normH1_2 + ((grad1_node_val-grad1ref)**2+(grad2_node_val-grad2ref)**2)*h1*h2
     end do
  end do
  
  print*, 'integrale de la solution =', sum(calculated(1:npts1-1,1:npts2-1))*h1*h2,&
       'integrale de la solution exacte=', sum(reference(1:npts1-1,1:npts2-1))*h1*h2
  call phi%write_to_file(0)


  ! delete things...
  call delete(es)
  call rho%delete()
  call c_field%delete()
  call phi%delete()
  call a11_field_mat%delete()
  call a12_field_mat%delete()
  call a21_field_mat%delete()
  call a22_field_mat%delete()
  call b1_field_vect%delete()
  call b2_field_vect%delete()
  call T%delete()
  
  SLL_DEALLOCATE(values, ierr)
  SLL_DEALLOCATE_ARRAY(calculated,ierr)
  SLL_DEALLOCATE_ARRAY(difference,ierr)
  SLL_DEALLOCATE_ARRAY(reference,ierr)
  
!!$  !--------------------------------------------------------------------
!!$  
!!$  !     3 test case without change of coordinates 
!!$  !      dirichlet-dirichlet boundary conditions
!!$  
!!$  !--------------------------------------------------------------------
!!$  
  
  
  print*, "-------------------------------------------------------------"
  print*, " 3 test case witout change of coordinates"
  print*, " dirichlet-dirichlet boundary conditions"
  print*, "-------------------------------------------------------------"
  npts1 =  NUM_CELLS1 + 1
  npts2 =  NUM_CELLS2 + 1
  h1 = (ETA1MAX-ETA1MIN)/real(NPTS1-1,f64)
  h2 = (ETA2MAX-ETA2MIN)/real(NPTS2-1,f64)
  print *, 'h1 = ', h1
  print *, 'h2 = ', h2
  
  ! Table to represent the node values of phi
  SLL_ALLOCATE(values(NUM_CELLS1+1,NUM_CELLS2+1),ierr)
  SLL_ALLOCATE(calculated(npts1,npts2),ierr)
  SLL_ALLOCATE(difference(npts1,npts2),ierr)
  SLL_ALLOCATE(reference(npts1,npts2),ierr)
  values(:,:) = 0.0_f64
  
  ! First thing, initialize the logical mesh associated with this problem.    
  mesh_2d => new_logical_mesh_2d( NUM_CELLS1, NUM_CELLS2, &
       ETA1MIN, ETA1MAX, ETA2MIN,ETA2MAX )
  
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
       identity_jac22, &
       (/0.0_f64/) )
  print *, 'initialized coordinate transformation'
  
  ! Thirdly, each field object must be initialized using the same logical
  ! mesh and coordinate transformation.
  a11_field_mat => new_scalar_field_2d_analytic_alt( &
       func_one, &
       "a11", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       whatever ) 
  
  a12_field_mat => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "a12", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       whatever ) 
  
  a21_field_mat => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "a21", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       whatever ) 
  
  a22_field_mat => new_scalar_field_2d_analytic_alt( &
       func_one, &
       "a22", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       whatever) 
  

  b1_field_vect => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "b1", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       whatever, &
       first_deriv_eta1 = func_zero, &
       first_deriv_eta2 = func_zero) 
  
  b2_field_vect => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "b2", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       whatever, &
       first_deriv_eta1 = func_zero, &
       first_deriv_eta2 = func_zero)
  
  c_field => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "c_field", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       whatever)

  rho => new_scalar_field_2d_analytic_alt( &
       source_term_perdir, &
       "rho3", &     
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       whatever )
  
  call initialize_ad2d_interpolator( &
       interp_2d, &
       NUM_CELLS1+1, &
       NUM_CELLS2+1, &
       ETA1MIN, &
       ETA1MAX, &
       ETA2MIN, &
       ETA2MAX, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SPLINE_DEG1, &
       SPLINE_DEG2 )
 
  
  phi => new_scalar_field_2d_discrete_alt( &
       "phi3", &
       interp_2d, &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET)
  call phi%set_field_data(values)
  call phi%update_interpolation_coefficients()
  
  print *, 'initialized fields...'
  
  call sll_set_time_mark(t_reference)

  call initialize_general_elliptic_solver( &
       es, &
       SPLINE_DEG1, &
       SPLINE_DEG2, &
       NUM_CELLS1, &
       NUM_CELLS2, &
       ES_GAUSS_LEGENDRE, &
       ES_GAUSS_LEGENDRE, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       ETA1MIN, &
       ETA1MAX, &
       ETA2MIN, &
       ETA2MAX)
  
  t3i = sll_time_elapsed_since(t_reference) 

  print *, 'Initialized ES object'

  call sll_set_time_mark(t_reference)

  call factorize_mat_es(&
       es, &
       a11_field_mat, &
       a12_field_mat,&
       a21_field_mat,&
       a22_field_mat,&
       b1_field_vect,&
       b2_field_vect,&
       c_field)!, &
      ! rho)

  ! solve the field
  call solve_general_coordinates_elliptic_eq(&
       es,&
       rho,&
       phi)

  !call interp_2d%set_coefficients( es%phi_vec)
 
  
  t3e = sll_time_elapsed_since(t_reference)


  acc3 = 0.0_f64
  normL2_3 = 0.0_f64
  normH1_3 = 0.0_f64
  do j=0,npts2-1
     do i=0,npts1-1
        eta1       = real(i,f64)*h1 + ETA1MIN
        eta2       = real(j,f64)*h2 + ETA2MIN

        node_val   =phi%value_at_point(eta1,eta2)
        grad1_node_val = phi%first_deriv_eta1_value_at_point(eta1, eta2)
        grad2_node_val = phi%first_deriv_eta2_value_at_point(eta1, eta2)

        ref        = sol_exacte_perdir(eta1,eta2)
        grad1ref   = sol_exacte_perdir_der1(eta1,eta2)
        grad2ref   = sol_exacte_perdir_der2(eta1,eta2)
        calculated(i+1,j+1) = node_val
        difference(i+1,j+1) = ref-node_val
        reference(i+1,j+1) = ref
        if(PRINT_COMPARISON) then
           print *, '(eta1,eta2) = ', eta1, eta2, 'calculated = ', node_val, &
                'theoretical = ', ref,'difference=', ref-node_val
        end if
        acc3        = acc3 + abs(node_val-ref)
        normL2_3    = normL2_3 + (node_val-ref)**2*h1*h2
        normH1_3    = normH1_3 + ((grad1_node_val-grad1ref)**2+(grad2_node_val-grad2ref)**2)*h1*h2
     end do
  end do
  print*, 'integrale de la solution =', sum(calculated(1:npts1-1,1:npts2-1))*h1*h2,&
       ' integrale de la solution exacte =', sum(reference(1:npts1-1,1:npts2-1))*h1*h2

  
  call phi%write_to_file(0)
  ! delete things...
  call delete(es)
  call rho%delete()
  call c_field%delete()
  call phi%delete()
  call a11_field_mat%delete()
  call a12_field_mat%delete()
  call a21_field_mat%delete()
  call a22_field_mat%delete()
  call b1_field_vect%delete()
  call b2_field_vect%delete()
  call T%delete()
  
  SLL_DEALLOCATE(values, ierr)
  SLL_DEALLOCATE_ARRAY(calculated,ierr)
  SLL_DEALLOCATE_ARRAY(difference,ierr)
  SLL_DEALLOCATE_ARRAY(reference,ierr)

!!$  !--------------------------------------------------------------------
!!$  
!!$  !     4 test case without change of coordinates 
!!$  !     dirichlet-periodic boundary conditions
!!$  
!!$  !--------------------------------------------------------------------
!!$  
  print*, "-------------------------------------------------------------"
  print*, " 4 test case witout change of coordinates"
  print*, " dirichlet-periodic boundary conditions"
  print*, "-------------------------------------------------------------"
  npts1 =  NUM_CELLS1 + 1
  npts2 =  NUM_CELLS2 + 1
  h1 = 1.0_f64/real(NPTS1-1,f64)
  h2 = 1.0_f64/real(NPTS2-1,f64)
  print *, 'h1 = ', h1
  print *, 'h2 = ', h2
  
  ! Table to represent the node values of phi
  SLL_ALLOCATE(values(NUM_CELLS1+1,NUM_CELLS2+1),ierr)
  SLL_ALLOCATE(calculated(npts1,npts2),ierr)
  SLL_ALLOCATE(difference(npts1,npts2),ierr)
  SLL_ALLOCATE(reference(npts1,npts2),ierr)
  values(:,:) = 0.0_f64
  
  ! First thing, initialize the logical mesh associated with this problem.        
  mesh_2d => new_logical_mesh_2d( NUM_CELLS1, NUM_CELLS2, &
       ETA1MIN, ETA1MAX, ETA2MIN,ETA2MAX )
  
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
       identity_jac22, &
       (/0.0_f64/) )
  print *, 'initialized coordinate transformation'
  
  ! Thirdly, each field object must be initialized using the same logical
  ! mesh and coordinate transformation.
  a11_field_mat => new_scalar_field_2d_analytic_alt( &
       func_one, &
       "a11", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       whatever ) 
  
  a12_field_mat => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "a12", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       whatever ) 
  
  a21_field_mat => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "a21", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       whatever ) 
  
  a22_field_mat => new_scalar_field_2d_analytic_alt( &
       func_one, &
       "a22", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       whatever ) 
  
  b1_field_vect => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "b1", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       whatever, &
       first_deriv_eta1 = func_zero, &
       first_deriv_eta2 = func_zero)

  b2_field_vect => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "b2", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       whatever, &
       first_deriv_eta1 = func_zero, &
       first_deriv_eta2 = func_zero)

  c_field => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "c_field", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       whatever)
  
  rho => new_scalar_field_2d_analytic_alt( &
       source_term_dirper, &
       "rho4", &     
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       whatever )
  
  call initialize_ad2d_interpolator( &
       interp_2d, &
       NUM_CELLS1+1, &
       NUM_CELLS2+1, &
       ETA1MIN, &
       ETA1MAX, &
       ETA2MIN, &
       ETA2MAX, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SPLINE_DEG1, &
       SPLINE_DEG2 )
  
  phi => new_scalar_field_2d_discrete_alt( &
       "phi4", &
       interp_2d, &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_PERIODIC, &
       SLL_PERIODIC)
  call phi%set_field_data(values)
  call phi%update_interpolation_coefficients()
  
  print *, 'initialized fields...'
  
  call sll_set_time_mark(t_reference)
  
  call initialize_general_elliptic_solver( &
       es, &
       SPLINE_DEG1, &
       SPLINE_DEG2, &
       NUM_CELLS1, &
       NUM_CELLS2, &
       ES_GAUSS_LEGENDRE, &
       ES_GAUSS_LEGENDRE, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       ETA1MIN, &
       ETA1MAX, &
       ETA2MIN, &
       ETA2MAX)
  
  t4i = sll_time_elapsed_since(t_reference) 
  print *, 'Initialized ES object'
  
  call sll_set_time_mark(t_reference)
  
  call factorize_mat_es(&
       es, &
       a11_field_mat, &
       a12_field_mat,&
       a21_field_mat,&
       a22_field_mat,&
       b1_field_vect,&
       b2_field_vect,&
       c_field)!, &
  
  ! solve the field
  call solve_general_coordinates_elliptic_eq(&
       es,&
       rho,&
       phi)
  t4e = sll_time_elapsed_since(t_reference) 
  
  
  acc4 = 0.0_f64
  normL2_4 = 0.0_f64
  normH1_4 = 0.0_f64
  do j=0,npts2-1
     do i=0,npts1-1
        eta1       = real(i,f64)*h1 + ETA1MIN
        eta2       = real(j,f64)*h2 + ETA2MIN
        node_val   = phi%value_at_point(eta1,eta2)
        grad1_node_val = phi%first_deriv_eta1_value_at_point(eta1, eta2)
        grad2_node_val = phi%first_deriv_eta2_value_at_point(eta1, eta2)
        ref        = sol_exacte_dirper(eta1,eta2)
        grad1ref   = sol_exacte_dirper_der1(eta1,eta2)
        grad2ref   = sol_exacte_dirper_der2(eta1,eta2)
        calculated(i+1,j+1) = node_val
        difference(i+1,j+1) = ref-node_val
        reference(i+1,j+1) = ref
        if(PRINT_COMPARISON) then
           print *, '(eta1,eta2) = ', eta1, eta2, 'calculated = ', node_val, &
                'theoretical = ', ref,'difference=',ref-node_val
        end if
        acc4        = acc4 + abs(node_val-ref)
        normL2_4    = normL2_4 + (node_val-ref)**2*h1*h2
        normH1_4    = normH1_4 + ((grad1_node_val-grad1ref)**2+(grad2_node_val-grad2ref)**2)*h1*h2
     end do
  end do
  print*, 'integrale de la solution =', sum(calculated(1:npts1-1,1:npts2-1))*h1*h2,&
       ' integrale de la solution exacte=',sum(reference(1:npts1-1,1:npts2-1))*h1*h2
  call phi%write_to_file(0)
  ! delete things...
  call delete(es)
  call rho%delete()
  call c_field%delete()
  call phi%delete()
  call a11_field_mat%delete()
  call a12_field_mat%delete()
  call a21_field_mat%delete()
  call a22_field_mat%delete()
  call b1_field_vect%delete()
  call b2_field_vect%delete()
  call T%delete()
  
  SLL_DEALLOCATE(values, ierr)
  SLL_DEALLOCATE_ARRAY(calculated,ierr)
  SLL_DEALLOCATE_ARRAY(difference,ierr)
  SLL_DEALLOCATE_ARRAY(reference,ierr)

  !****************************************************************************
  !        WHITH CHANGE OF COORDINATES AND ANALYTIC DATA
  !*****************************************************************************
!!$  !--------------------------------------------------------------------
!!$  
!!$  !     5 test case with colella change of coordinates 
!!$  !     periodic-periodic boundary conditions
!!$  
!!$  !--------------------------------------------------------------------
!!$  
  print*, "---------------------"
  print*, " 5 test case with colella change of coordinates"
  print*, " periodic-periodic boundary conditions"
  print*, "---------------------"
  npts1 =  NUM_CELLS1 + 1
  npts2 =  NUM_CELLS2 + 1
  h1 = (ETA1MAX - ETA1MIN)/real(NPTS1-1,f64)
  h2 = (ETA2MAX - ETA2MIN)/real(NPTS2-1,f64)
  print *, 'h1 = ', h1
  print *, 'h2 = ', h2
  
  ! Table to represent the node values of phi
  SLL_ALLOCATE(values(NUM_CELLS1+1,NUM_CELLS2+1),ierr)
  SLL_ALLOCATE(calculated(npts1,npts2),ierr)
  SLL_ALLOCATE(difference(npts1,npts2),ierr)
  SLL_ALLOCATE(reference(npts1,npts2),ierr)
  values(:,:) = 0.0_f64
  
  ! First thing, initialize the logical mesh associated with this problem.  
  mesh_2d => new_logical_mesh_2d( NUM_CELLS1, NUM_CELLS2, &
       ETA1MIN, ETA1MAX, ETA2MIN,ETA2MAX )
  
  ! Second, initialize the coordinate transformation associated with this 
  ! problem.
  T => new_coordinate_transformation_2d_analytic( &
       "analytic", &
       mesh_2d, &
       sinprod_x1, &
       sinprod_x2, &
       sinprod_jac11, &
       sinprod_jac12, &
       sinprod_jac21, &
       sinprod_jac22, &
       (/ 0.1_f64, 0.1_f64, 1.0_f64, 1.0_f64/) )
  print *, 'initialized coordinate transformation'
  
  ! Thirdly, each field object must be initialized using the same logical
  ! mesh and coordinate transformation.
  a11_field_mat => new_scalar_field_2d_analytic_alt( &
       func_one, &
       "a11", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       whatever ) 
  
  a12_field_mat => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "a12", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       whatever ) 
  
  a21_field_mat => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "a21", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       whatever ) 
  
  a22_field_mat => new_scalar_field_2d_analytic_alt( &
       func_one, &
       "a22", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       whatever ) 
  

  b1_field_vect => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "b1", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       whatever, &
       first_deriv_eta1 = func_zero, &
       first_deriv_eta2 = func_zero)

  b2_field_vect => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "b2", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       whatever, &
       first_deriv_eta1 = func_zero, &
       first_deriv_eta2 = func_zero)


  c_field => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "c_field", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       whatever)
  
  rho => new_scalar_field_2d_analytic_alt( &
       source_term_chgt_perper, &
       "rho5", &     
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       whatever )
  
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
  
  
  phi => new_scalar_field_2d_discrete_alt( &
       "phi5", &
       interp_2d, &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC)
  call phi%set_field_data(values)
  call phi%update_interpolation_coefficients()
  
  print *, 'initialized fields...'
  
  call sll_set_time_mark(t_reference)
  
  call initialize_general_elliptic_solver( &
       es, &
       SPLINE_DEG1, &
       SPLINE_DEG2, &
       NUM_CELLS1, &
       NUM_CELLS2, &
       ES_GAUSS_LEGENDRE, &
       ES_GAUSS_LEGENDRE, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       ETA1MIN, &
       ETA1MAX, &
       ETA2MIN, &
       ETA2MAX)
  
  t5i = sll_time_elapsed_since(t_reference) 
  print *, 'Initialized ES object'
  
  call sll_set_time_mark(t_reference)  

  call factorize_mat_es(&
       es, &
       a11_field_mat, &
       a12_field_mat,&
       a21_field_mat,&
       a22_field_mat,&
       b1_field_vect,&
       b2_field_vect,&
       c_field)!, &
  
  ! solve the field
  call solve_general_coordinates_elliptic_eq(&
       es,&
       rho,&
       phi)
  !print *, 'Completed solution',es%phi_vec
  
  t5e = sll_time_elapsed_since(t_reference)  
  
  !  print *, 'Compare the values of the transformation at the nodes: '
  
  acc5 = 0.0_f64
  integrale_solution = 0.0
  integrale_solution_exacte = 0.0
  normH1_5 =  0.0
  normL2_5 =  0.0
  do j=0,npts2-1
     do i=0,npts1-1
        eta1       = real(i,f64)*h1 + ETA1MIN
        eta2       = real(j,f64)*h2 + ETA2MIN
        node_val   = phi%value_at_point(eta1,eta2)
        grad1_node_val = phi%first_deriv_eta1_value_at_point(eta1, eta2)
        grad2_node_val = phi%first_deriv_eta2_value_at_point(eta1, eta2)
        !print*, 'rer'
        ref        = sol_exacte_chgt_perper(eta1,eta2)
        grad1ref   = sol_exacte_chgt_perper_der1(eta1,eta2)
        grad2ref   = sol_exacte_chgt_perper_der2(eta1,eta2)
        calculated(i+1,j+1) = node_val
        difference(i+1,j+1) = ref-node_val
        reference(i+1,j+1) = ref
        ! jac_mat(:,:) = mesh_2d%get_jacobian_matrix(eta1,eta2)
        val_jac = sinprod_jac11(eta1,eta2,(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))*&
                  sinprod_jac22(eta1,eta2,(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))-&
                  sinprod_jac12(eta1,eta2,(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))*&
                  sinprod_jac21(eta1,eta2,(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))
        if(PRINT_COMPARISON) then! k = 5 0.00000001 k = 3 0.00001, k=4 0.0000001
           print *, '(eta1,eta2) = ', eta1, eta2, 'calculated = ', node_val, &
                'theoretical = ', ref, 'difference = ', ref-node_val
           
        end if
        acc5        = acc5 + abs(node_val-ref)
        if ( i < npts1-1 .and. j < npts2-1) then
           integrale_solution = integrale_solution + node_val*val_jac * h1*h2
           integrale_solution_exacte = integrale_solution_exacte + ref*val_jac * h1*h2
           normL2_5    = normL2_5 + (node_val-ref)**2*h1*h2*val_jac
           normH1_5    = normH1_5 + ((grad1_node_val-grad1ref)**2+(grad2_node_val-grad2ref)**2)*h1*h2*val_jac
        end if
     end do
  end do
  
  print*, ' integrale solution', integrale_solution, &
       'integrale de la solution exacte=', integrale_solution_exacte
  call phi%write_to_file(0)
  ! delete things...
  call delete(es)
  call rho%delete()
  call c_field%delete()
  call phi%delete()
  call a11_field_mat%delete()
  call a12_field_mat%delete()
  call a21_field_mat%delete()
  call a22_field_mat%delete()
  call b1_field_vect%delete()
  call b2_field_vect%delete()
  call T%delete()
  
  SLL_DEALLOCATE(values, ierr)
  SLL_DEALLOCATE_ARRAY(calculated,ierr)
  SLL_DEALLOCATE_ARRAY(difference,ierr)
  SLL_DEALLOCATE_ARRAY(reference,ierr)



!!$
!!$  !--------------------------------------------------------------------
!!$  
!!$  !     6 test case with colella change of coordinates 
!!$  !     periodic-dirichlet boundary conditions
!!$  
!!$  !--------------------------------------------------------------------
!!$  
!!$  
  print*, "-------------------------------------------------------------"
  print*, " 6 test case with colella change of coordinates"
  print*, " periodic-dirichlet boundary conditions"
  print*, "-------------------------------------------------------------"
  npts1 =  NUM_CELLS1 + 1
  npts2 =  NUM_CELLS2 + 1
  h1 = (ETA1MAX - ETA1MIN)/real(NPTS1-1,f64)
  h2 = (ETA2MAX - ETA2MIN)/real(NPTS2-1,f64)
  print *, 'h1 = ', h1
  print *, 'h2 = ', h2
  
  ! Table to represent the node values of phi
  SLL_ALLOCATE(values(NUM_CELLS1+1,NUM_CELLS2+1),ierr)
  SLL_ALLOCATE(calculated(npts1,npts2),ierr)
  SLL_ALLOCATE(difference(npts1,npts2),ierr)
  SLL_ALLOCATE(reference(npts1,npts2),ierr)

  values(:,:) = 0.0_f64
  
  ! First thing, initialize the logical mesh associated with this problem.        
  mesh_2d => new_logical_mesh_2d( NUM_CELLS1, NUM_CELLS2, &
       ETA1MIN, ETA1MAX, ETA2MIN,ETA2MAX )
  
  ! Second, initialize the coordinate transformation associated with this 
  ! problem.
  T => new_coordinate_transformation_2d_analytic( &
       "analytic", &
       mesh_2d, &
       sinprod_x1, &
       sinprod_x2, &
       sinprod_jac11, &
       sinprod_jac12, &
       sinprod_jac21, &
       sinprod_jac22, &
       (/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))
  print *, 'initialized coordinate transformation'
  
  ! Thirdly, each field object must be initialized using the same logical
  ! mesh and coordinate transformation.
  a11_field_mat => new_scalar_field_2d_analytic_alt( &
       func_one, &
       "a11", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       whatever ) 
  
  a12_field_mat => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "a12", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       whatever ) 
  
  a21_field_mat => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "a21", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       whatever ) 
  
  a22_field_mat => new_scalar_field_2d_analytic_alt( &
       func_one, &
       "a22", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       whatever ) 
  
  b1_field_vect => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "b1", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       whatever, &
       first_deriv_eta1 = func_zero, &
       first_deriv_eta2 = func_zero)
  
  b2_field_vect => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "b2", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       whatever, &
       first_deriv_eta1 = func_zero, &
       first_deriv_eta2 = func_zero)
  
  c_field => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "c_field", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       whatever )

  rho => new_scalar_field_2d_analytic_alt( &
       source_term_chgt_perdir, &
       "rho6", &     
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       whatever )
  
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
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SPLINE_DEG1, &
       SPLINE_DEG2 )
  
  phi => new_scalar_field_2d_discrete_alt( &
       "phi6", &
       interp_2d, &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET)
  call phi%set_field_data(values)
  call phi%update_interpolation_coefficients()
  
  print *, 'initialized fields...'
  
  !  print *, 'a = ', es%csr_mat%opr_a
  
  call sll_set_time_mark(t_reference)

  call initialize_general_elliptic_solver( &
       es, &
       SPLINE_DEG1, &
       SPLINE_DEG2, &
       NUM_CELLS1, &
       NUM_CELLS2, &
       ES_GAUSS_LEGENDRE, &
       ES_GAUSS_LEGENDRE, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       ETA1MIN, &
       ETA1MAX, &
       ETA2MIN, &
       ETA2MAX)

  t6i = sll_time_elapsed_since(t_reference) 
  
  print *, 'Initialized ES object'
  
  call sll_set_time_mark(t_reference)

  call factorize_mat_es(&
       es, &
       a11_field_mat, &
       a12_field_mat,&
       a21_field_mat,&
       a22_field_mat,&
       b1_field_vect,&
       b2_field_vect,&
       c_field)!, &
  
  ! solve the field
  call solve_general_coordinates_elliptic_eq(&
       es,&
       rho,&
       phi)
  t6e = sll_time_elapsed_since(t_reference)
  
  
  !  print *, 'Compare the values of the transformation at the nodes: '
  
  acc6 = 0.0_f64
  integrale_solution = 0.0_f64
  integrale_solution_exacte = 0.0_f64
  normL2_6 = 0.0_f64
  normH1_6 = 0.0_f64
  do j=0,npts2-1
     do i=0,npts1-1
        eta1       = real(i,f64)*h1 + ETA1MIN
        eta2       = real(j,f64)*h2 + ETA2MIN
        
        node_val   = phi%value_at_point(eta1,eta2)
        grad1_node_val = phi%first_deriv_eta1_value_at_point(eta1, eta2)
        grad2_node_val = phi%first_deriv_eta2_value_at_point(eta1, eta2)
        !print*, 'rer'
        ref        = sol_exacte_chgt_perdir(eta1,eta2)
        grad1ref   = sol_exacte_chgt_perdir_der1(eta1,eta2)
        grad2ref   = sol_exacte_chgt_perdir_der2(eta1,eta2)
        calculated(i+1,j+1) = node_val
        difference(i+1,j+1) = ref-node_val
        reference(i+1,j+1) = ref
        
        if(PRINT_COMPARISON) then
           print *, '(eta1,eta2) = ', eta1, eta2, 'calculated = ', node_val, &
                'theoretical = ', ref,'difference=', ref-node_val
        end if
        acc6        = acc6 + abs(node_val-ref)
        
        val_jac = sinprod_jac11(eta1,eta2,(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))*&
                  sinprod_jac22(eta1,eta2,(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))-&
                  sinprod_jac12(eta1,eta2,(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))*&
                  sinprod_jac21(eta1,eta2,(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))
        if ( i < npts1-1 .and. j < npts2-1 ) then
           integrale_solution = integrale_solution + node_val*val_jac* h1*h2
           integrale_solution_exacte = integrale_solution_exacte + ref*val_jac* h1*h2
           normL2_6    = normL2_6 + (node_val-ref)**2*h1*h2
           normH1_6    = normH1_6 + ((grad1_node_val-grad1ref)**2+(grad2_node_val-grad2ref)**2)*h1*h2*val_jac
           
        end if
     end do
  end do
  
  print*, 'integrale solution=',integrale_solution,&
       'integrale de la solution exacte=', integrale_solution_exacte
  call phi%write_to_file(0)
  ! delete things...
  call delete(es)
  call rho%delete()
  call c_field%delete()
  call phi%delete()
  call a11_field_mat%delete()
  call a12_field_mat%delete()
  call a21_field_mat%delete()
  call a22_field_mat%delete()
  call b1_field_vect%delete()
  call b2_field_vect%delete()
  call T%delete()
  
  SLL_DEALLOCATE(values, ierr)
  SLL_DEALLOCATE_ARRAY(calculated,ierr)
  SLL_DEALLOCATE_ARRAY(difference,ierr)
  SLL_DEALLOCATE_ARRAY(reference,ierr)


  
  !--------------------------------------------------------------------
  
  !     7 test case with colella change of coordinates 
  !     dirichlet-dirichlet boundary conditions
  
  !--------------------------------------------------------------------
  
  print*, "-------------------------------------------------------------"
  print*, " 7 test case with colella change of coordinates"
  print*, " dirichlet-dirichlet boundary conditions"
  print*, "-------------------------------------------------------------"
  npts1 =  NUM_CELLS1 + 1
  npts2 =  NUM_CELLS2 + 1
  h1 = (ETA1MAX - ETA1MIN)/real(NPTS1-1,f64)
  h2 = (ETA2MAX - ETA2MIN)/real(NPTS2-1,f64)
  print *, 'h1 = ', h1
  print *, 'h2 = ', h2
  
  ! Table to represent the node values of phi
  SLL_ALLOCATE(values(NUM_CELLS1+1,NUM_CELLS2+1),ierr)
  SLL_ALLOCATE(calculated(npts1,npts2),ierr)
  SLL_ALLOCATE(difference(npts1,npts2),ierr)
  SLL_ALLOCATE(reference(npts1,npts2),ierr)
  values(:,:) = 0.0_f64

  ! First thing, initialize the logical mesh associated with this problem.        
  mesh_2d => new_logical_mesh_2d( NUM_CELLS1, NUM_CELLS2, &
       ETA1MIN, ETA1MAX, ETA2MIN,ETA2MAX )
  
  ! Second, initialize the coordinate transformation associated with this 
  ! problem.
  T => new_coordinate_transformation_2d_analytic( &
       "analytic", &
       mesh_2d, &
       sinprod_x1, &
       sinprod_x2, &
       sinprod_jac11, &
       sinprod_jac12, &
       sinprod_jac21, &
       sinprod_jac22, &
       (/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))
  print *, 'initialized coordinate transformation'
  
   ! Thirdly, each field object must be initialized using the same logical
  ! mesh and coordinate transformation.
  a11_field_mat => new_scalar_field_2d_analytic_alt( &
       func_one, &
       "a11", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       whatever ) 
  
  a12_field_mat => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "a12", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       whatever) 
  
  a21_field_mat => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "a21", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       whatever ) 
  
  a22_field_mat => new_scalar_field_2d_analytic_alt( &
       func_one, &
       "a22", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       whatever) 
  
  b1_field_vect => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "b1", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       whatever, &
       first_deriv_eta1 = func_zero, &
       first_deriv_eta2 = func_zero)

  b2_field_vect => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "b2", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       whatever, &
       first_deriv_eta1 = func_zero, &
       first_deriv_eta2 = func_zero)


  c_field => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "c_field", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       whatever )
  
  rho => new_scalar_field_2d_analytic_alt( &
       source_term_chgt_dirdir, &
       "rho7", &     
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       whatever )
  
  call initialize_ad2d_interpolator( &
       interp_2d, &
       NUM_CELLS1+1, &
       NUM_CELLS2+1, &
       ETA1MIN, &
       ETA1MAX, &
       ETA2MIN, &
       ETA2MAX, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SPLINE_DEG1, &
       SPLINE_DEG2 )
  
  phi => new_scalar_field_2d_discrete_alt( &
       "phi7", &
       interp_2d, &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET)
  call phi%set_field_data(values)
  call phi%update_interpolation_coefficients()
  
  print *, 'initialized fields...'
  
  call sll_set_time_mark(t_reference)
  call initialize_general_elliptic_solver( &
       es, &
       SPLINE_DEG1, &
       SPLINE_DEG2, &
       NUM_CELLS1, &
       NUM_CELLS2, &
       ES_GAUSS_LEGENDRE, &
       ES_GAUSS_LEGENDRE, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       ETA1MIN, &
       ETA1MAX, &
       ETA2MIN, &
       ETA2MAX)
  
  t7i = sll_time_elapsed_since(t_reference) 
  
  print *, 'Initialized ES object'

  call sll_set_time_mark(t_reference)

  call factorize_mat_es(&
       es, &
       a11_field_mat, &
       a12_field_mat,&
       a21_field_mat,&
       a22_field_mat,&
       b1_field_vect,&
       b2_field_vect,&
       c_field)!, &
  
  ! solve the field
  call solve_general_coordinates_elliptic_eq(&
       es,&
       rho,&
       phi)
!!$ !!$
  t7e = sll_time_elapsed_since(t_reference)
  
  
!  print *, 'Compare the values of the transformation at the nodes: '
  
  acc7 = 0.0_f64
  integrale_solution = 0.0_f64
  integrale_solution_exacte = 0.0_f64
  normL2_7 = 0.0_f64
  normH1_7 = 0.0_f64
  do j=0,npts2-1
     do i=0,npts1-1
        eta1       = real(i,f64)*h1 + ETA1MIN
        eta2       = real(j,f64)*h2 + ETA2MIN
        
        node_val   = phi%value_at_point(eta1,eta2)
        grad1_node_val = phi%first_deriv_eta1_value_at_point(eta1, eta2)
        grad2_node_val = phi%first_deriv_eta2_value_at_point(eta1, eta2)
        !print*, 'rer'
        ref        = sol_exacte_chgt_dirdir(eta1,eta2)
        grad1ref   = sol_exacte_chgt_dirdir_der1(eta1,eta2)
        grad2ref   = sol_exacte_chgt_dirdir_der2(eta1,eta2)
        calculated(i+1,j+1) = node_val
        difference(i+1,j+1) = ref-node_val
        reference(i+1,j+1) = ref
        if(PRINT_COMPARISON) then
           print *, '(eta1,eta2) = ', eta1, eta2, 'calculated = ', node_val, &
                'theoretical = ', ref,'difference=',ref-node_val
        end if
        acc7        = acc7 + abs(node_val-ref)
        val_jac = sinprod_jac11(eta1,eta2,(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))*&
                  sinprod_jac22(eta1,eta2,(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))-&
                  sinprod_jac12(eta1,eta2,(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))*&
                  sinprod_jac21(eta1,eta2,(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))
        if ( i < npts1-1 .and. j < npts2-1) then
           integrale_solution = integrale_solution + node_val*val_jac * h1*h2
           integrale_solution_exacte = integrale_solution_exacte + ref*val_jac * h1*h2
           normL2_7    = normL2_7 + (node_val-ref)**2*h1*h2*val_jac
           normH1_7    = normH1_7 + ((grad1_node_val-grad1ref)**2+(grad2_node_val-grad2ref)**2)*h1*h2*val_jac
           
        end if
     end do
  end do
  print*, 'integrale solution=',integrale_solution,&
       'integrale de la solution excate=', integrale_solution_exacte
  call phi%write_to_file(0)
  ! delete things...
  call delete(es)
  call rho%delete()
  call c_field%delete()
  call phi%delete()
  call a11_field_mat%delete()
  call a12_field_mat%delete()
  call a21_field_mat%delete()
  call a22_field_mat%delete()
  call b1_field_vect%delete()
  call b2_field_vect%delete()
  call T%delete()

  SLL_DEALLOCATE(values, ierr)
  SLL_DEALLOCATE_ARRAY(calculated,ierr)
  SLL_DEALLOCATE_ARRAY(difference,ierr)
  SLL_DEALLOCATE_ARRAY(reference,ierr)

  

!!$  !--------------------------------------------------------------------
!!$  
!!$  !     8  test case with colella change of coordinates 
!!$  !     dirichlet-periodic boundary conditions
!!$  
!!$  !--------------------------------------------------------------------
!!$  
!!$  
  print*, "---------------------"
  print*, " 8 test case with colella change of coordinates"
  print*, " dirichlet-periodic boundary conditions"
  print*, "---------------------"
  npts1 =  NUM_CELLS1 + 1
  npts2 =  NUM_CELLS2 + 1
  h1 = (ETA1MAX - ETA1MIN)/real(NPTS1-1,f64)
  h2 = (ETA2MAX - ETA2MIN)/real(NPTS2-1,f64)
  print *, 'h1 = ', h1
  print *, 'h2 = ', h2
  
  ! Table to represent the node values of phi
  SLL_ALLOCATE(values(NUM_CELLS1+1,NUM_CELLS2+1),ierr)
  SLL_ALLOCATE(calculated(npts1,npts2),ierr)
  SLL_ALLOCATE(difference(npts1,npts2),ierr)
  SLL_ALLOCATE(reference(npts1,npts2),ierr)
  values(:,:) = 0.0_f64
  
  ! First thing, initialize the logical mesh associated with this problem.        
  mesh_2d => new_logical_mesh_2d( NUM_CELLS1, NUM_CELLS2, &
       ETA1MIN, ETA1MAX, ETA2MIN,ETA2MAX )
  
  ! Second, initialize the coordinate transformation associated with this 
  ! problem.
  T => new_coordinate_transformation_2d_analytic( &
       "analytic", &
       mesh_2d, &
       sinprod_x1, &
       sinprod_x2, &
       sinprod_jac11, &
       sinprod_jac12, &
       sinprod_jac21, &
       sinprod_jac22, &
       (/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))
  print *, 'initialized coordinate transformation'
  
  ! Thirdly, each field object must be initialized using the same logical
  ! mesh and coordinate transformation.
  a11_field_mat => new_scalar_field_2d_analytic_alt( &
       func_one, &
       "a11", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       whatever ) 
  
  a12_field_mat => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "a12", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET,&
       SLL_PERIODIC,&
       SLL_PERIODIC, &
       whatever )
  
  a21_field_mat => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "a21", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET,&
       SLL_PERIODIC,&
       SLL_PERIODIC, &
       whatever ) 
  
  a22_field_mat => new_scalar_field_2d_analytic_alt( &
       func_one, &
       "a22", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET,&
       SLL_PERIODIC,&
       SLL_PERIODIC,&
       whatever) 
  
  b1_field_vect => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "b1", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       whatever, &
       first_deriv_eta1 = func_zero, &
       first_deriv_eta2 = func_zero)

  b2_field_vect => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "b2", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       whatever, &
       first_deriv_eta1 = func_zero, &
       first_deriv_eta2 = func_zero)


  c_field => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "c_field", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       whatever)
  
  rho => new_scalar_field_2d_analytic_alt( &
       source_term_chgt_dirper, &
       "rho8", &     
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET,&
       SLL_PERIODIC,&
       SLL_PERIODIC, &
       whatever)
  
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
       SPLINE_DEG2 )
  
  
  phi => new_scalar_field_2d_discrete_alt( &
       "phi8", &
       interp_2d, &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET,&
       SLL_PERIODIC, &
       SLL_PERIODIC)
  call phi%set_field_data(values)
  call phi%update_interpolation_coefficients()
  
  print *, 'initialized fields...'
  
  call sll_set_time_mark(t_reference)
  
  call initialize_general_elliptic_solver( &
       es, &
       SPLINE_DEG1, &
       SPLINE_DEG2, &
       NUM_CELLS1, &
       NUM_CELLS2, &
       ES_GAUSS_LEGENDRE, &
       ES_GAUSS_LEGENDRE, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET,&
       SLL_PERIODIC,&
       SLL_PERIODIC,&
       ETA1MIN, &
       ETA1MAX, &
       ETA2MIN, &
       ETA2MAX)
  
  t8i = sll_time_elapsed_since(t_reference) 
  
  print *, 'Initialized ES object'
  call sll_set_time_mark(t_reference)
  
  call factorize_mat_es(&
       es, &
       a11_field_mat, &
       a12_field_mat,&
       a21_field_mat,&
       a22_field_mat,&
       b1_field_vect,&
       b2_field_vect,&
       c_field)!, &
  
  ! solve the field
  call solve_general_coordinates_elliptic_eq(&
       es,&
       rho,&
       phi)
  
  !print *, 'Completed solution',es%phi_vec
  
  t8e = sll_time_elapsed_since(t_reference)
  
  
  
  !  print *, 'Compare the values of the transformation at the nodes: '
  
  acc8 = 0.0_f64
  integrale_solution = 0.0_f64
  integrale_solution_exacte = 0.0_f64
  normH1_8 =  0.0
  normL2_8 =  0.0
  do j=0,npts2-1
     do i=0,npts1-1
        eta1       = real(i,f64)*h1 + ETA1MIN
        eta2       = real(j,f64)*h2 + ETA2MIN
        
        node_val   =phi%value_at_point(eta1,eta2)
        grad1_node_val = phi%first_deriv_eta1_value_at_point(eta1, eta2)
        grad2_node_val = phi%first_deriv_eta2_value_at_point(eta1, eta2)
        !print*, 'rer'
        ref        = sol_exacte_chgt_dirper(eta1,eta2)
        grad1ref   = sol_exacte_chgt_dirper_der1(eta1,eta2)
        grad2ref   = sol_exacte_chgt_dirper_der2(eta1,eta2)
        calculated(i+1,j+1) = node_val
        difference(i+1,j+1) = ref-node_val
        reference(i+1,j+1) = ref
        if(PRINT_COMPARISON) then
           print *, '(eta1,eta2) = ', eta1, eta2, 'calculated = ', node_val, &
                'theoretical = ', ref,'difference=',ref-node_val
        end if
        acc8        = acc8 + abs(node_val-ref)
        val_jac = sinprod_jac11(eta1,eta2,(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))*&
                  sinprod_jac22(eta1,eta2,(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))-&
                  sinprod_jac12(eta1,eta2,(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))*&
                  sinprod_jac21(eta1,eta2,(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))
        if ( i < npts1-1 .and. j < npts2-1) then
           integrale_solution = integrale_solution + node_val*val_jac * h1*h2
           integrale_solution_exacte = integrale_solution_exacte + ref*val_jac * h1*h2
           normL2_8    = normL2_8 + (node_val-ref)**2*h1*h2*val_jac
           normH1_8    = normH1_8 + ((grad1_node_val-grad1ref)**2+(grad2_node_val-grad2ref)**2)*h1*h2*val_jac
        end if
     end do
  end do
  print*, 'integrale solution=',integrale_solution,&
       'integrale de la solution exacte=', integrale_solution_exacte
  call phi%write_to_file(0)
  ! delete things...
  call delete(es)
  call rho%delete()
  call c_field%delete()
  call phi%delete()
  call a11_field_mat%delete()
  call a12_field_mat%delete()
  call a21_field_mat%delete()
  call a22_field_mat%delete()
  call b1_field_vect%delete()
  call b2_field_vect%delete()
  call T%delete()
  
  SLL_DEALLOCATE(values, ierr)
  SLL_DEALLOCATE_ARRAY(calculated,ierr)
  SLL_DEALLOCATE_ARRAY(difference,ierr)
  SLL_DEALLOCATE_ARRAY(reference,ierr)

  



 !******************************************************************************
  !        WHITH CHANGE OF COORDINATES AND NON-ANALYTIC SOURCE TERM
  !****************************************************************************
!!$  !--------------------------------------------------------------------
!!$  
!!$  !     95  test case without change of coordinates 
!!$  !      periodic-periodic boundary conditions
!!$  !      and with a non analytic source term
!!$  !      with separate factorization of solve 
!!$  
!!$  !--------------------------------------------------------------------
!!$  
!!$  
!!$  
  print*, "---------------------"
  print*, " 95 test case without change of coordinates"
  print*, " periodic-periodic boundary conditions"
  print*, " with non analytic source term " 
  print*, "---------------------"
  npts1 =  NUM_CELLS1 + 1
  npts2 =  NUM_CELLS2 + 1
  h1 = (ETA1MAX - ETA1MIN)/real(NPTS1-1,f64)
  h2 = (ETA2MAX - ETA2MIN)/real(NPTS2-1,f64)
  print *, 'h1 = ', h1
  print *, 'h2 = ', h2
  
  ! Table to represent the node values of phi
  SLL_ALLOCATE(values(NUM_CELLS1+1,NUM_CELLS2+1),ierr)
  SLL_ALLOCATE(calculated(npts1,npts2),ierr)
  SLL_ALLOCATE(difference(npts1,npts2),ierr)
  SLL_ALLOCATE(reference(npts1,npts2),ierr)
  values(:,:) = 0.0_f64
  
  ! First thing, initialize the logical mesh associated with this problem.        
  mesh_2d => new_logical_mesh_2d( NUM_CELLS1, NUM_CELLS2, &
       ETA1MIN, ETA1MAX, ETA2MIN,ETA2MAX )
  
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
       identity_jac22, &
       (/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))
  print *, 'initialized coordinate transformation'
  
  ! Thirdly, each field object must be initialized using the same logical
  ! mesh and coordinate transformation.
  a11_field_mat => new_scalar_field_2d_analytic_alt( &
       func_one, &
       "a11", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       whatever ) 
  
  a12_field_mat => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "a12", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC,&
       SLL_PERIODIC,&
       SLL_PERIODIC, &
       whatever )
  
  a21_field_mat => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "a21", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC,&
       SLL_PERIODIC,&
       SLL_PERIODIC, &
       whatever ) 
  
  a22_field_mat => new_scalar_field_2d_analytic_alt( &
       func_one, &
       "a22", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC,&
       SLL_PERIODIC,&
       SLL_PERIODIC,&
       whatever)
    
  b1_field_vect => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "b1", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       whatever, &
       first_deriv_eta1 = func_zero, &
       first_deriv_eta2 = func_zero)

  b2_field_vect => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "b2", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       whatever, &
       first_deriv_eta1 = func_zero, &
       first_deriv_eta2 = func_zero)
  
  c_field => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "c_field", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       whatever)
  
  
  allocate(point1(npts1))
  allocate(point2(npts2))
  allocate(tab_rho(npts1,npts2))
  do j=0,npts2-1
     do i=0,npts1-1
        point1(i+1)       = real(i,f64)*(ETA1MAX-ETA1MIN)/(npts1-1) + ETA1MIN 
        point2(j+1)       = real(j,f64)*(ETA2MAX-ETA2MIN)/(npts2-1) + ETA2MIN 
        tab_rho(i+1,j+1)  = source_term_perper(point1(i+1),point2(j+1))
        !print*, tab_rho(i+1,j+1), -2*(2.0*sll_pi)**2*cos(2.0*sll_pi* point1(i+1))*cos(2.0*sll_pi* point2(j+1))
     end do
  end do
  

  print*, 'sum',sum(tab_rho)*h2*h2
  call initialize_ad2d_interpolator( &
       interp_2d_term_source, &
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
  
  ! terme_source_interp => interp_2d_term_source
  
  rho => new_scalar_field_2d_discrete_alt( &
       "rho95", &
       interp_2d_term_source, &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC,&
       SLL_PERIODIC,&
       SLL_PERIODIC,&
       point1,&
       npts1-1,&
       point2,&
       npts2-1)  
  call rho%set_field_data(tab_rho)
  call rho%update_interpolation_coefficients()

 call rho%write_to_file(0)
 
 
!!$  
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
       SPLINE_DEG2 )
!!$  
  
  phi => new_scalar_field_2d_discrete_alt( &
       "phi95", &
       interp_2d, &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC,&
       SLL_PERIODIC, &
       SLL_PERIODIC)
  call phi%set_field_data(values)
  call phi%update_interpolation_coefficients()
  
  print *, 'initialized fields...'
  
  call sll_set_time_mark(t_reference)
  
  call initialize_general_elliptic_solver( &
       es, &
       SPLINE_DEG1, &
       SPLINE_DEG2, &
       NUM_CELLS1, &
       NUM_CELLS2, &
       ES_GAUSS_LEGENDRE, &
       ES_GAUSS_LEGENDRE, &
       SLL_PERIODIC, &
       SLL_PERIODIC,&
       SLL_PERIODIC,&
       SLL_PERIODIC,&
       ETA1MIN, &
       ETA1MAX, &
       ETA2MIN, &
       ETA2MAX)
  
  t95i = sll_time_elapsed_since(t_reference) 
  
  print *, 'Initialized ES object'
  call sll_set_time_mark(t_reference)
  
  call factorize_mat_es(&
       es, &
       a11_field_mat, &
       a12_field_mat,&
       a21_field_mat,&
       a22_field_mat,&
       b1_field_vect,&
       b2_field_vect,&
       c_field)!, &

  ! solve the field
  call solve_general_coordinates_elliptic_eq(&
       es,&
       rho,&
       phi)
  
  !print *, 'Completed solution',es%phi_vec
  
  t95e = sll_time_elapsed_since(t_reference)
  
  print *, 'Compare the values of the transformation at the nodes: '
  
  acc95 = 0.0_f64
  integrale_solution = 0.0_f64
  integrale_solution_exacte = 0.0_f64
  normL2_95 = 0.0
  normH1_95 = 0.0
  do j=0,npts2-1
     do i=0,npts1-1
        eta1       = real(i,f64)*h1 + ETA1MIN
        eta2       = real(j,f64)*h2 + ETA2MIN
        
        node_val   =     phi%value_at_point(eta1,eta2)
        grad1_node_val = phi%first_deriv_eta1_value_at_point(eta1, eta2)
        grad2_node_val = phi%first_deriv_eta2_value_at_point(eta1, eta2)
        !print*, 'value at node', node_val
        !print*, 'rer'
        ref        = sol_exacte_perper(eta1,eta2)
        grad1ref   = sol_exacte_perper_der1(eta1,eta2)
        grad2ref   = sol_exacte_perper_der2(eta1,eta2)
        calculated(i+1,j+1) = node_val
        difference(i+1,j+1) = ref-node_val
        reference(i+1,j+1) = ref
        val_jac = 1.0
        if(PRINT_COMPARISON) then
           print *, '(eta1,eta2) = ', eta1, eta2, 'calculated = ', node_val, &
                'theoretical = ', ref, 'difference=',ref-node_val
           !end if
         !  print *, '(eta1,eta2) = ', eta1, eta2, 'calculated = ', grad1_node_val, &
          !      'theoretical = ', grad1ref, 'difference=',grad1ref-grad1_node_val
        end if
        acc95        = acc95 + abs(node_val-ref)
        if ( i < npts1-1 .and. j < npts2-1) then
           integrale_solution = integrale_solution + node_val * h1*h2
           integrale_solution_exacte = integrale_solution_exacte + ref * h1*h2
           normL2_95    = normL2_95 + (node_val-ref)**2*h1*h2*val_jac
           normH1_95    = normH1_95 + ((grad1_node_val-grad1ref)**2+(grad2_node_val-grad2ref)**2)*h1*h2*val_jac
          ! print*, normH1_95,normL2_95,grad1_node_val,grad2_node_val,grad1ref,grad2ref,val_jac
        end if
     end do
  end do



  print*, 'integrale solution=', integrale_solution,&
       'integrale de la solution exacte=', integrale_solution_exacte
  call phi%write_to_file(0)
  ! delete things...
  call delete(es)
  call rho%delete()
  call c_field%delete()
  call phi%delete()
  call a11_field_mat%delete()
  call a12_field_mat%delete()
  call a21_field_mat%delete()
  call a22_field_mat%delete()
  call b1_field_vect%delete()
  call b2_field_vect%delete()
  call T%delete()
  
  SLL_DEALLOCATE(values, ierr)
  SLL_DEALLOCATE_ARRAY(calculated,ierr)
  SLL_DEALLOCATE_ARRAY(difference,ierr)
  SLL_DEALLOCATE_ARRAY(reference,ierr)
  SLL_DEALLOCATE_ARRAY(point1,ierr)
  SLL_DEALLOCATE_ARRAY(point2,ierr)
  SLL_DEALLOCATE_ARRAY(tab_rho,ierr)



  
  

  !--------------------------------------------------------------------
  
  !     9  test case with colella change of coordinates 
  !      periodic-periodic boundary conditions
  !      and with a non analytic source term
  
  !--------------------------------------------------------------------
!!$  
!!$  
!!$  
  print*, "---------------------"
  print*, " 9 test case with colella change of coordinates"
  print*, " periodic-periodic boundary conditions"
  print*, " with non analytic source term " 
  print*, "---------------------"
  npts1 =  NUM_CELLS1 + 1
  npts2 =  NUM_CELLS2 + 1
  h1 = (ETA1MAX - ETA1MIN)/real(NPTS1-1,f64)
  h2 = (ETA2MAX - ETA2MIN)/real(NPTS2-1,f64)
  print *, 'h1 = ', h1
  print *, 'h2 = ', h2
  
  ! Table to represent the node values of phi
  SLL_ALLOCATE(values(NUM_CELLS1+1,NUM_CELLS2+1),ierr)
  SLL_ALLOCATE(calculated(npts1,npts2),ierr)
  SLL_ALLOCATE(difference(npts1,npts2),ierr)
  SLL_ALLOCATE(reference(npts1,npts2),ierr)
  values(:,:) = 0.0_f64
  
  ! First thing, initialize the logical mesh associated with this problem.        
  mesh_2d => new_logical_mesh_2d( NUM_CELLS1, NUM_CELLS2, &
       ETA1MIN, ETA1MAX, ETA2MIN,ETA2MAX )
  
  ! Second, initialize the coordinate transformation associated with this 
  ! problem.
  T => new_coordinate_transformation_2d_analytic( &
       "analytic", &
       mesh_2d, &
       sinprod_x1, &
       sinprod_x2, &
       sinprod_jac11, &
       sinprod_jac12, &
       sinprod_jac21, &
       sinprod_jac22, &
       (/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))
  print *, 'initialized coordinate transformation'
  
  ! Thirdly, each field object must be initialized using the same logical
  ! mesh and coordinate transformation.
  a11_field_mat => new_scalar_field_2d_analytic_alt( &
       func_one, &
       "a11", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       whatever  ) 
  
  a12_field_mat => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "a12", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC,&
       SLL_PERIODIC,&
       SLL_PERIODIC, &
       whatever )
  
  a21_field_mat => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "a21", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC,&
       SLL_PERIODIC,&
       SLL_PERIODIC, &
       whatever ) 
  
  a22_field_mat => new_scalar_field_2d_analytic_alt( &
       func_one, &
       "a22", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC,&
       SLL_PERIODIC,&
       SLL_PERIODIC, &
       whatever)
  
  
  b1_field_vect => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "b1", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       whatever, &
       first_deriv_eta1 = func_zero, &
       first_deriv_eta2 = func_zero)

  b2_field_vect => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "b2", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       whatever, &
       first_deriv_eta1 = func_zero, &
       first_deriv_eta2 = func_zero)
  
  c_field => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "c_field", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       whatever)
  
  

  allocate(point1(npts1))
  allocate(point2(npts2))
  allocate(tab_rho(npts1,npts2))
  
  do j=0,npts2-1
     do i=0,npts1-1
        point1(i+1)       = real(i,f64)*h1 + ETA1MIN
        point2(j+1)       = real(j,f64)*h2 + ETA2MIN
        tab_rho(i+1,j+1)  = source_term_chgt_perper(point1(i+1),point2(j+1))
     end do
  end do
  
  
  call initialize_ad2d_interpolator( &
       interp_2d_term_source, &
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
       SPLINE_DEG2 )
  
  terme_source_interp => interp_2d_term_source
  
  tab_rho(:,:) = tab_rho - sum(tab_rho)/((npts1-1)*(npts2-1))
  print*,'moyenne', sum(tab_rho)
  rho => new_scalar_field_2d_discrete_alt( &
       "rho9", &
       terme_source_interp, &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC,&
       SLL_PERIODIC,&
       SLL_PERIODIC ,&
       point1,&
       npts1-1,&
       point2,&
       npts2-1)
  call rho%set_field_data(tab_rho)
  call rho%update_interpolation_coefficients()
  
!!$  
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
       SPLINE_DEG2 )
!!$  
  
  phi => new_scalar_field_2d_discrete_alt( &
       "phi9", &
       interp_2d, &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC,&
       SLL_PERIODIC, &
       SLL_PERIODIC)
  call phi%set_field_data(values)
  call phi%update_interpolation_coefficients()
  
  print *, 'initialized fields...'
  
  call sll_set_time_mark(t_reference)
  
  call initialize_general_elliptic_solver( &
       es, &
       SPLINE_DEG1, &
       SPLINE_DEG2, &
       NUM_CELLS1, &
       NUM_CELLS2, &
       ES_GAUSS_LEGENDRE, &
       ES_GAUSS_LEGENDRE, &
       SLL_PERIODIC, &
       SLL_PERIODIC,&
       SLL_PERIODIC,&
       SLL_PERIODIC,&
       ETA1MIN, &
       ETA1MAX, &
       ETA2MIN, &
       ETA2MAX)
  
  t9i = sll_time_elapsed_since(t_reference) 
  
  print *, 'Initialized ES object'
  call sll_set_time_mark(t_reference)
  

  call factorize_mat_es(&
       es, &
       a11_field_mat, &
       a12_field_mat,&
       a21_field_mat,&
       a22_field_mat,&
       b1_field_vect,&
       b2_field_vect,&
       c_field)!, &
  
  ! solve the field
  call solve_general_coordinates_elliptic_eq(&
       es,&
       rho,&
       phi)
  !print *, 'Completed solution',es%phi_vec
  
  t9e = sll_time_elapsed_since(t_reference)
  
  print *, 'Compare the values of the transformation at the nodes: '
  
  acc9 = 0.0_f64
  integrale_solution = 0.0_f64
  integrale_solution_exacte = 0.0_f64
  normL2_9 = 0.0
  normH1_9 = 0.0
  do j=0,npts2-1
     do i=0,npts1-1
        eta1       = real(i,f64)*h1 + ETA1MIN
        eta2       = real(j,f64)*h2 + ETA2MIN
        
        node_val   =phi%value_at_point(eta1,eta2)
        grad1_node_val = phi%first_deriv_eta1_value_at_point(eta1, eta2)
        grad2_node_val = phi%first_deriv_eta2_value_at_point(eta1, eta2)
        !print*, 'value at node', node_val
        !print*, 'rer'
        ref        = sol_exacte_chgt_perper(eta1,eta2)
        grad1ref   = sol_exacte_chgt_perper_der1(eta1,eta2)
        grad2ref   = sol_exacte_chgt_perper_der2(eta1,eta2)
        calculated(i+1,j+1) = node_val
        difference(i+1,j+1) = ref-node_val
        reference(i+1,j+1) = ref
        if(PRINT_COMPARISON) then
           print *, '(eta1,eta2) = ', eta1, eta2, 'calculated = ', node_val, &
                'theoretical = ', ref,'difference=',ref-node_val
        end if
        acc9        = acc9 + abs(node_val-ref)
        
        val_jac = sinprod_jac11(eta1,eta2,(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))*&
                  sinprod_jac22(eta1,eta2,(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))-&
                  sinprod_jac12(eta1,eta2,(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))*&
                  sinprod_jac21(eta1,eta2,(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))
        
        if ( i < npts1-1 .and. j < npts2-1) then
           integrale_solution = integrale_solution + node_val*val_jac * h1*h2
           integrale_solution_exacte = integrale_solution_exacte + ref*val_jac * h1*h2
           normL2_9    = normL2_9 + (node_val-ref)**2*h1*h2*val_jac
           normH1_9    = normH1_9 + ((grad1_node_val-grad1ref)**2+(grad2_node_val-grad2ref)**2)*h1*h2*val_jac
           
        end if
     end do
  end do
  
  print*, 'integrale solution', integrale_solution,&
       ' integrale de la solution exacte=', integrale_solution_exacte
  call phi%write_to_file(0)
  call rho%write_to_file(0)
  ! delete things...
  call delete(es)
  call rho%delete()
  call c_field%delete()
  call phi%delete()
  call a11_field_mat%delete()
  call a12_field_mat%delete()
  call a21_field_mat%delete()
  call a22_field_mat%delete()
  call b1_field_vect%delete()
  call b2_field_vect%delete()
  call T%delete()
  
  SLL_DEALLOCATE(values, ierr)
  SLL_DEALLOCATE_ARRAY(calculated,ierr)
  SLL_DEALLOCATE_ARRAY(difference,ierr)
  SLL_DEALLOCATE_ARRAY(reference,ierr)
  SLL_DEALLOCATE_ARRAY(point1,ierr)
  SLL_DEALLOCATE_ARRAY(point2,ierr)
  SLL_DEALLOCATE_ARRAY(tab_rho,ierr)
  
  
  
   !--------------------------------------------------------------------
  
  !     10  test case with colella change of coordinates 
  !      periodic-dirichlet boundary conditions
  !     and non analytic source term
  !--------------------------------------------------------------------
!!$  
  
  
  print*, "---------------------"
  print*, " 10 test case with colella change of coordinates"
  print*, " periodic-dirichlet boundary conditions"
  print*, " with non analytic source term " 
  print*, "---------------------"
  npts1 =  NUM_CELLS1 + 1
  npts2 =  NUM_CELLS2 + 1
  h1 = (ETA1MAX - ETA1MIN)/real(NPTS1-1,f64)
  h2 = (ETA2MAX - ETA2MIN)/real(NPTS2-1,f64)
  print *, 'h1 = ', h1
  print *, 'h2 = ', h2
  
  ! Table to represent the node values of phi
  SLL_ALLOCATE(values(NUM_CELLS1+1,NUM_CELLS2+1),ierr)
  SLL_ALLOCATE(calculated(npts1,npts2),ierr)
  SLL_ALLOCATE(difference(npts1,npts2),ierr)
   SLL_ALLOCATE(reference(npts1,npts2),ierr)
  values(:,:) = 0.0_f64

  ! First thing, initialize the logical mesh associated with this problem.        
  mesh_2d => new_logical_mesh_2d( NUM_CELLS1, NUM_CELLS2, &
       ETA1MIN, ETA1MAX, ETA2MIN,ETA2MAX )
  
  ! Second, initialize the coordinate transformation associated with this 
  ! problem.
   T => new_coordinate_transformation_2d_analytic( &
       "analytic", &
       mesh_2d, &
       sinprod_x1, &
       sinprod_x2, &
       sinprod_jac11, &
       sinprod_jac12, &
       sinprod_jac21, &
       sinprod_jac22, &
       (/0.1_f64, 0.1_f64, 1.0_f64, 1.0_f64/)) 

  print *, 'initialized coordinate transformation'
  
  ! Thirdly, each field object must be initialized using the same logical
  ! mesh and coordinate transformation.
  a11_field_mat => new_scalar_field_2d_analytic_alt( &
       func_one, &
       "a11", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       whatever ) 
  
  a12_field_mat => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "a12", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC,&
       SLL_DIRICHLET,&
       SLL_DIRICHLET, &
       whatever )
  
  a21_field_mat => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "a21", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC,&
       SLL_DIRICHLET,&
       SLL_DIRICHLET, &
       whatever ) 
  
  a22_field_mat => new_scalar_field_2d_analytic_alt( &
       func_one, &
       "a22", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC,&
       SLL_DIRICHLET,&
       SLL_DIRICHLET, &
       whatever)
  

  b1_field_vect => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "b1", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       whatever, &
       first_deriv_eta1 = func_zero, &
       first_deriv_eta2 = func_zero)
  
  b2_field_vect => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "b2", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       whatever, &
       first_deriv_eta1 = func_zero, &
       first_deriv_eta2 = func_zero)
 
  
  c_field => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "c_field", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       whatever)

  

  allocate(point1(npts1))
  allocate(point2(npts2))
  allocate(tab_rho(npts1,npts2))
  do j=0,npts2-1
     do i=0,npts1-1
        point1(i+1)       = real(i,f64)*h1 + ETA1MIN
        point2(j+1)       = real(j,f64)*h2 + ETA2MIN
        tab_rho(i+1,j+1)  = source_term_chgt_perdir(point1(i+1),point2(j+1))
     end do
  end do

  call initialize_ad2d_interpolator( &
       interp_2d_term_source, &
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
       SPLINE_DEG2 )
  
  terme_source_interp => interp_2d_term_source


  rho => new_scalar_field_2d_discrete_alt( &
       "rho10", &
       terme_source_interp, &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC,&
       SLL_DIRICHLET,&
       SLL_DIRICHLET,&
       point1,&
       npts1-1,&
       point2,&
       npts2)
  call rho%set_field_data(tab_rho)
  call rho%update_interpolation_coefficients()

!!$  
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
       SPLINE_DEG2 )
!!$  
  
  phi => new_scalar_field_2d_discrete_alt( &
       "phi10", &
       interp_2d, &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC,&
       SLL_DIRICHLET, &
       SLL_DIRICHLET)
  call phi%set_field_data(values)
  call phi%update_interpolation_coefficients()
  
  print *, 'initialized fields...'
  
  call sll_set_time_mark(t_reference)
  
  call initialize_general_elliptic_solver( &
       es, &
       SPLINE_DEG1, &
       SPLINE_DEG2, &
       NUM_CELLS1, &
       NUM_CELLS2, &
       ES_GAUSS_LEGENDRE, &
       ES_GAUSS_LEGENDRE, &
       SLL_PERIODIC, &
       SLL_PERIODIC,&
       SLL_DIRICHLET,&
       SLL_DIRICHLET,&
       ETA1MIN, &
       ETA1MAX, &
       ETA2MIN, &
       ETA2MAX)

  t10i = sll_time_elapsed_since(t_reference) 
  
  print *, 'Initialized ES object'
  call sll_set_time_mark(t_reference)
  

  call factorize_mat_es(&
       es, &
       a11_field_mat, &
       a12_field_mat,&
       a21_field_mat,&
       a22_field_mat,&
       b1_field_vect,&
       b2_field_vect,&
       c_field)!, &
  
  ! solve the field
  call solve_general_coordinates_elliptic_eq(&
       es,&
       rho,&
       phi)
!!$  
  t10e = sll_time_elapsed_since(t_reference)
  
  print *, 'Compare the values of the transformation at the nodes: '
  
  acc10 = 0.0_f64
  integrale_solution = 0.0 
  integrale_solution_exacte = 0.0
  normH1_10 =  0.0
  normL2_10 =  0.0
  do j=0,npts2-1
     do i=0,npts1-1
        eta1       = real(i,f64)*h1 + ETA1MIN
        eta2       = real(j,f64)*h2 + ETA2MIN
        
        node_val   =phi%value_at_point(eta1,eta2)
        grad1_node_val = phi%first_deriv_eta1_value_at_point(eta1, eta2)
        grad2_node_val = phi%first_deriv_eta2_value_at_point(eta1, eta2)
        !print*, 'rer'
        ref        = sol_exacte_chgt_perdir(eta1,eta2)
        grad1ref   = sol_exacte_chgt_perdir_der1(eta1,eta2)
        grad2ref   = sol_exacte_chgt_perdir_der2(eta1,eta2)
        calculated(i+1,j+1) = node_val
        difference(i+1,j+1) = ref-node_val
        reference(i+1,j+1) = ref
        if(PRINT_COMPARISON) then
           print *, '(eta1,eta2) = ', eta1, eta2, 'calculated = ', node_val, &
                'theoretical = ', ref,'difference=',ref-node_val
        end if
        acc10        = acc10 + abs(node_val-ref)
        val_jac = sinprod_jac11(eta1,eta2,(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))*&
                  sinprod_jac22(eta1,eta2,(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))-&
                  sinprod_jac12(eta1,eta2,(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))*&
                  sinprod_jac21(eta1,eta2,(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))
        
        if ( i < npts1-1 .and. j < npts2-1) then
           integrale_solution = integrale_solution + node_val*val_jac * h1*h2
           integrale_solution_exacte = integrale_solution_exacte + &
                ref*val_jac * h1*h2
           normL2_10    = normL2_10 + (node_val-ref)**2*h1*h2*val_jac
           normH1_10   = normH1_10 + ((grad1_node_val-grad1ref)**2+&
                (grad2_node_val-grad2ref)**2)*h1*h2*val_jac
           
        end if
     end do
  end do
  print*, 'integrale de la solution=',integrale_solution,&
       'integrale de la solution exacte=',integrale_solution_exacte

  call phi%write_to_file(0)
  ! delete things...
  call delete(es)
  call rho%delete()
  call c_field%delete()
  call phi%delete()
  call a11_field_mat%delete()
  call a12_field_mat%delete()
  call a21_field_mat%delete()
  call a22_field_mat%delete()
  call b1_field_vect%delete()
  call b2_field_vect%delete()
  call T%delete()
  
  SLL_DEALLOCATE(values, ierr)
  SLL_DEALLOCATE_ARRAY(calculated,ierr)
  SLL_DEALLOCATE_ARRAY(difference,ierr)
  SLL_DEALLOCATE_ARRAY(reference,ierr)
  SLL_DEALLOCATE_ARRAY(point1,ierr)
  SLL_DEALLOCATE_ARRAY(point2,ierr)
  SLL_DEALLOCATE_ARRAY(tab_rho,ierr)

!!$
!!$  !--------------------------------------------------------------------
!!$  
!!$  !     11  test case with colella change of coordinates 
!!$  !      dirichlet-dirichlet boundary conditions
!!$  !     and non analytic source term
!!$  !--------------------------------------------------------------------
!!$  
!!$  
!!$  
  print*, "---------------------"
  print*, " 11 test case with colella change of coordinates"
  print*, " dirichlet-dirichlet boundary conditions"
  print*, " with non analytic source term " 
  print*, "---------------------"
  npts1 =  NUM_CELLS1 + 1
  npts2 =  NUM_CELLS2 + 1
  h1 = (ETA1MAX - ETA1MIN)/real(NPTS1-1,f64)
  h2 = (ETA2MAX - ETA2MIN)/real(NPTS2-1,f64)
  print *, 'h1 = ', h1
  print *, 'h2 = ', h2
  
  ! Table to represent the node values of phi
  SLL_ALLOCATE(values(NUM_CELLS1+1,NUM_CELLS2+1),ierr)
  SLL_ALLOCATE(calculated(npts1,npts2),ierr)
  SLL_ALLOCATE(difference(npts1,npts2),ierr)
  SLL_ALLOCATE(reference(npts1,npts2),ierr)

  values(:,:) = 0.0_f64
  
  ! First thing, initialize the logical mesh associated with this problem.        
  mesh_2d => new_logical_mesh_2d( NUM_CELLS1, NUM_CELLS2, &
       ETA1MIN, ETA1MAX, ETA2MIN,ETA2MAX )
  
  ! Second, initialize the coordinate transformation associated with this 
  ! problem.
   T => new_coordinate_transformation_2d_analytic( &
       "analytic", &
       mesh_2d, &
       sinprod_x1, &
       sinprod_x2, &
       sinprod_jac11, &
       sinprod_jac12, &
       sinprod_jac21, &
       sinprod_jac22, &
       (/0.1_f64, 0.1_f64, 1.0_f64, 1.0_f64/) )

  print *, 'initialized coordinate transformation'
  
  ! Thirdly, each field object must be initialized using the same logical
  ! mesh and coordinate transformation.
  a11_field_mat => new_scalar_field_2d_analytic_alt( &
       func_one, &
       "a11", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       whatever ) 
  
  a12_field_mat => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "a12", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET,&
       SLL_DIRICHLET,&
       SLL_DIRICHLET, &
       whatever )
  
  a21_field_mat => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "a21", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET,&
       SLL_DIRICHLET,&
       SLL_DIRICHLET, &
       whatever ) 
  
  a22_field_mat => new_scalar_field_2d_analytic_alt( &
       func_one, &
       "a22", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET,&
       SLL_DIRICHLET,&
       SLL_DIRICHLET,&
       whatever) 
  
  b1_field_vect => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "b1", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       whatever, &
       first_deriv_eta1 = func_zero, &
       first_deriv_eta2 = func_zero)
  
  b2_field_vect => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "b2", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       whatever, &
       first_deriv_eta1 = func_zero, &
       first_deriv_eta2 = func_zero)
  
  c_field => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "c_field", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       whatever)



  allocate(point1(npts1))
  allocate(point2(npts2))
  allocate(tab_rho(npts1,npts2))
  do j=0,npts2-1
     do i=0,npts1-1
        point1(i+1)       = real(i,f64)*h1 + ETA1MIN
        point2(j+1)       = real(j,f64)*h2 + ETA2MIN
        tab_rho(i+1,j+1)  = source_term_chgt_dirdir(point1(i+1),point2(j+1))
     end do
  end do

  call initialize_ad2d_interpolator( &
       interp_2d_term_source, &
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
       SPLINE_DEG2 )
  
  terme_source_interp => interp_2d_term_source



  rho => new_scalar_field_2d_discrete_alt( &
       "rho11", &
       terme_source_interp, &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET,&
       SLL_DIRICHLET,&
       SLL_DIRICHLET,&
       point1,&
       npts1,&
       point2,&
       npts2)
  call rho%set_field_data(tab_rho)
  call rho%update_interpolation_coefficients()

!!$  
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
       SPLINE_DEG2 )
!!$    
  phi => new_scalar_field_2d_discrete_alt( &
       "phi11", &
       interp_2d, &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET,&
       SLL_DIRICHLET, &
       SLL_DIRICHLET)
  call phi%set_field_data(values)
  call phi%update_interpolation_coefficients()
  
  print *, 'initialized fields...'
  call sll_set_time_mark(t_reference)
  
  call initialize_general_elliptic_solver( &
       es, &
       SPLINE_DEG1, &
       SPLINE_DEG2, &
       NUM_CELLS1, &
       NUM_CELLS2, &
       ES_GAUSS_LEGENDRE, &
       ES_GAUSS_LEGENDRE, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET,&
       SLL_DIRICHLET,&
       SLL_DIRICHLET,&
       ETA1MIN, &
       ETA1MAX, &
       ETA2MIN, &
       ETA2MAX)
  
  t11i = sll_time_elapsed_since(t_reference) 

  !print *, 'Initialized ES object'
  call sll_set_time_mark(t_reference)

  call factorize_mat_es(&
       es, &
       a11_field_mat, &
       a12_field_mat,&
       a21_field_mat,&
       a22_field_mat,&
       b1_field_vect,&
       b2_field_vect,&
       c_field)!, &

  !print*, 'solve'
  ! solve the field
  call solve_general_coordinates_elliptic_eq(&
       es,&
       rho,&
       phi)
 ! print*, 'aye'
!!$  
  t11e = sll_time_elapsed_since(t_reference)
  
  print *, 'Compare the values of the transformation at the nodes: '
  
  acc11 = 0.0_f64
  integrale_solution_exacte = 0.0
  integrale_solution = 0.0
  normL2_11 = 0.0_f64
  normH1_11 = 0.0_f64
  do j=0,npts2-1
     do i=0,npts1-1
        eta1       = real(i,f64)*h1 + ETA1MIN
        eta2       = real(j,f64)*h2 + ETA2MIN
        
        node_val   =phi%value_at_point(eta1,eta2)
        grad1_node_val = phi%first_deriv_eta1_value_at_point(eta1, eta2)
        grad2_node_val = phi%first_deriv_eta2_value_at_point(eta1, eta2)
        !print*, 'rer'
        ref        = sol_exacte_chgt_dirdir(eta1,eta2)
        grad1ref   = sol_exacte_chgt_dirdir_der1(eta1,eta2)
        grad2ref   = sol_exacte_chgt_dirdir_der2(eta1,eta2)
        calculated(i+1,j+1) = node_val
        difference(i+1,j+1) = ref-node_val
        reference(i+1,j+1)  = ref

        if(PRINT_COMPARISON) then
           print *, '(eta1,eta2) = ', eta1, eta2, 'calculated = ', node_val, &
                'theoretical = ', ref, 'difference=', ref-node_val
        end if
        acc11        = acc11 + abs(node_val-ref)
        val_jac = sinprod_jac11(eta1,eta2,(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))*&
                  sinprod_jac22(eta1,eta2,(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))-&
                  sinprod_jac12(eta1,eta2,(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))*&
                  sinprod_jac21(eta1,eta2,(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))
        
        if ( i < npts1-1 .and. j < npts2-1) then
           integrale_solution = integrale_solution + node_val*val_jac * h1*h2
           integrale_solution_exacte = integrale_solution_exacte + &
                ref*val_jac * h1*h2
           normL2_11    = normL2_11 + (node_val-ref)**2*h1*h2*val_jac
           normH1_11    = normH1_11 + ((grad1_node_val-grad1ref)**2+&
                (grad2_node_val-grad2ref)**2)*h1*h2*val_jac
        end if
     end do
  end do
  print*, 'integrale de la solution=', integrale_solution,&
       'integrale de la solution exacte=', integrale_solution_exacte
  
  call phi%write_to_file(0)
  ! delete things...
  call delete(es)
  call rho%delete()
  call c_field%delete()
  call phi%delete()
  call a11_field_mat%delete()
  call a12_field_mat%delete()
  call a21_field_mat%delete()
  call a22_field_mat%delete()

  call b1_field_vect%delete()
  call b2_field_vect%delete()
  
  call T%delete()

  SLL_DEALLOCATE(values, ierr)
  SLL_DEALLOCATE_ARRAY(calculated,ierr)
  SLL_DEALLOCATE_ARRAY(difference,ierr)
  SLL_DEALLOCATE_ARRAY(reference,ierr)
  SLL_DEALLOCATE_ARRAY(point1,ierr)
  SLL_DEALLOCATE_ARRAY(point2,ierr)
  SLL_DEALLOCATE_ARRAY(tab_rho,ierr)

!!$  !--------------------------------------------------------------------
!!$  
!!$  !     12  test case with colella change of coordinates 
!!$  !      dirichlet-periodic boundary conditions
!!$  !     and non analytic source term
!!$  !--------------------------------------------------------------------
!!$  
!!$  
  
  print*, "---------------------"
  print*, " 12 test case with colella change of coordinates"
  print*, " dirichlet-periodic  boundary conditions"
  print*, " with non analytic source term " 
  print*, "---------------------"
  npts1 =  NUM_CELLS1 + 1
  npts2 =  NUM_CELLS2 + 1
  h1 = (ETA1MAX - ETA1MIN)/real(NPTS1-1,f64)
  h2 = (ETA2MAX - ETA2MIN)/real(NPTS2-1,f64)
  print *, 'h1 = ', h1
  print *, 'h2 = ', h2
  
  ! Table to represent the node values of phi
  SLL_ALLOCATE(values(NUM_CELLS1+1,NUM_CELLS2+1),ierr)
  SLL_ALLOCATE(calculated(npts1,npts2),ierr)
  SLL_ALLOCATE(difference(npts1,npts2),ierr)
  SLL_ALLOCATE(reference(npts1,npts2),ierr)
  values(:,:) = 0.0_f64
  
  ! First thing, initialize the logical mesh associated with this problem.        
  mesh_2d => new_logical_mesh_2d( NUM_CELLS1, NUM_CELLS2, &
       ETA1MIN, ETA1MAX, ETA2MIN,ETA2MAX )
  
  ! Second, initialize the coordinate transformation associated with this 
  ! problem.
  T => new_coordinate_transformation_2d_analytic( &
       "analytic", &
       mesh_2d, &
       sinprod_x1, &
       sinprod_x2, &
       sinprod_jac11, &
       sinprod_jac12, &
       sinprod_jac21, &
       sinprod_jac22, &
       (/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))
  print *, 'initialized coordinate transformation'
  
  ! Thirdly, each field object must be initialized using the same logical
  ! mesh and coordinate transformation.
  a11_field_mat => new_scalar_field_2d_analytic_alt( &
       func_one, &
       "a11", &
       T, &
       SLL_DIRICHLET,&
       SLL_DIRICHLET,&
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       whatever ) 
  
  a12_field_mat => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "a12", &
       T, &
       SLL_DIRICHLET,&
       SLL_DIRICHLET,&
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       whatever )
  
  a21_field_mat => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "a21", &
       T, &
       SLL_DIRICHLET,&
       SLL_DIRICHLET,&
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       whatever ) 
  
  a22_field_mat => new_scalar_field_2d_analytic_alt( &
       func_one, &
       "a22", &
       T, &
       SLL_DIRICHLET,&
       SLL_DIRICHLET,&
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       whatever) 
  

  b1_field_vect => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "b1", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       whatever, &
       first_deriv_eta1 = func_zero, &
       first_deriv_eta2 = func_zero)
  
  b2_field_vect => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "b2", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       whatever, &
       first_deriv_eta1 = func_zero, &
       first_deriv_eta2 = func_zero) 
  
  c_field => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "c_field", &
       T, &
       SLL_DIRICHLET,&
       SLL_DIRICHLET,&
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       whatever)
  
  
  
  SLL_ALLOCATE(point1(npts1),ierr)
  SLL_ALLOCATE(point2(npts2),ierr)
  SLL_ALLOCATE(tab_rho(npts1,npts2),ierr)
  do j=0,npts2-1
     do i=0,npts1-1
        point1(i+1)       = real(i,f64)*h1 + ETA1MIN
        point2(j+1)       = real(j,f64)*h2 + ETA2MIN
        tab_rho(i+1,j+1)  = source_term_chgt_dirper(point1(i+1),point2(j+1))
     end do
  end do
  
  call initialize_ad2d_interpolator( &
       interp_2d_term_source, &
       NUM_CELLS1+1, &
       NUM_CELLS2+1, &
       ETA1MIN, &
       ETA1MAX, &
       ETA2MIN, &
       ETA2MAX, &
       SLL_DIRICHLET,&
       SLL_DIRICHLET,&
       SLL_PERIODIC, &
       SLL_PERIODIC,&
       SPLINE_DEG1, &
       SPLINE_DEG2 )
  
  terme_source_interp => interp_2d_term_source


  rho => new_scalar_field_2d_discrete_alt( &
       "rho12", &
       terme_source_interp, &
       T, &
       SLL_DIRICHLET,&
       SLL_DIRICHLET,&
       SLL_PERIODIC, &
       SLL_PERIODIC,&
       point1,&
       npts1,&
       point2,&
       npts2-1)
  call rho%set_field_data(tab_rho)
  call rho%update_interpolation_coefficients()

!!$  
  call initialize_ad2d_interpolator( &
       interp_2d, &
       NUM_CELLS1+1, &
       NUM_CELLS2+1, &
       ETA1MIN, &
       ETA1MAX, &
       ETA2MIN, &
       ETA2MAX, &
       SLL_DIRICHLET,&
       SLL_DIRICHLET,&
       SLL_PERIODIC, &
       SLL_PERIODIC,&
       SPLINE_DEG1, &
       SPLINE_DEG2 )
!!$  
  phi => new_scalar_field_2d_discrete_alt( &
       "phi12", &
       interp_2d, &
       T, &
       SLL_DIRICHLET,&
       SLL_DIRICHLET,&
       SLL_PERIODIC, &
       SLL_PERIODIC)
  call phi%set_field_data(values)
  call phi%update_interpolation_coefficients()
   
  print *, 'initialized fields...'
  
  call sll_set_time_mark(t_reference)

  call initialize_general_elliptic_solver( &
       es, &
       SPLINE_DEG1, &
       SPLINE_DEG2, &
       NUM_CELLS1, &
       NUM_CELLS2, &
       ES_GAUSS_LEGENDRE, &
       ES_GAUSS_LEGENDRE, &
       SLL_DIRICHLET,&
       SLL_DIRICHLET,&
       SLL_PERIODIC, &
       SLL_PERIODIC,&
       ETA1MIN, &
       ETA1MAX, &
       ETA2MIN, &
       ETA2MAX)

  t12i = sll_time_elapsed_since(t_reference)

  call factorize_mat_es(&
       es, &
       a11_field_mat, &
       a12_field_mat,&
       a21_field_mat,&
       a22_field_mat,&
       b1_field_vect,&
       b2_field_vect,&
       c_field)!, &

  ! solve the field
  call solve_general_coordinates_elliptic_eq(&
       es,&
       rho,&
       phi)
  
  t12e = sll_time_elapsed_since(t_reference)
  
  print *, 'Compare the values of the transformation at the nodes: '
  
  acc12 = 0.0_f64
  integrale_solution = 0.0
  integrale_solution_exacte = 0.0
  normL2_12 = 0.0_f64
  normH1_12 = 0.0_f64
  
  do j=0,npts2-1
     do i=0,npts1-1
        eta1       = real(i,f64)*h1 + ETA1MIN
        eta2       = real(j,f64)*h2 + ETA2MIN
        
        node_val   = phi%value_at_point(eta1,eta2)
        grad1_node_val = phi%first_deriv_eta1_value_at_point(eta1, eta2)
        grad2_node_val = phi%first_deriv_eta2_value_at_point(eta1, eta2)
        !print*, 'rer'
        ref        = sol_exacte_chgt_dirper(eta1,eta2)
        grad1ref   = sol_exacte_chgt_dirper_der1(eta1,eta2)
        grad2ref   = sol_exacte_chgt_dirper_der2(eta1,eta2)
        calculated(i+1,j+1) = node_val
        difference(i+1,j+1) = ref-node_val
        reference(i+1,j+1)  = ref
        if(PRINT_COMPARISON) then
           print *, '(eta1,eta2) = ', eta1, eta2, 'calculated = ', node_val, &
                'theoretical = ', ref,'difference=', ref-node_val
        end if
        acc12        = acc12 + abs(node_val-ref)
        val_jac = sinprod_jac11(eta1,eta2,(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))*&
                  sinprod_jac22(eta1,eta2,(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))-&
                  sinprod_jac12(eta1,eta2,(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))*&
                  sinprod_jac21(eta1,eta2,(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))
        
        if ( i < npts1-1 .and. j < npts2-1) then
           integrale_solution = integrale_solution + node_val*val_jac * h1*h2
           integrale_solution_exacte = integrale_solution_exacte + &
                ref*val_jac * h1*h2
           normL2_12    = normL2_12 + (node_val-ref)**2*h1*h2*val_jac
           normH1_12    = normH1_12 + ((grad1_node_val-grad1ref)**2+&
                (grad2_node_val-grad2ref)**2)*h1*h2*val_jac
           
        end if
     end do
  end do
  print*, 'integrale de la solution =', integrale_solution,&
       'integrale de la solution exacte=', integrale_solution_exacte
  call phi%write_to_file(0)
  ! delete things...
  call delete(es)
  call rho%delete()
  call c_field%delete()
  call phi%delete()
  call a11_field_mat%delete()
  call a12_field_mat%delete()
  call a21_field_mat%delete()
  call a22_field_mat%delete()
  call b1_field_vect%delete()
  call b2_field_vect%delete()
  
  call T%delete()

  SLL_DEALLOCATE(values, ierr)
  SLL_DEALLOCATE_ARRAY(calculated,ierr)
  SLL_DEALLOCATE_ARRAY(difference,ierr)
  SLL_DEALLOCATE_ARRAY(reference,ierr)
  SLL_DEALLOCATE_ARRAY(point1,ierr)
  SLL_DEALLOCATE_ARRAY(point2,ierr)
  SLL_DEALLOCATE_ARRAY(tab_rho,ierr)


  print*, '------------------------------------------------------'
  print*, ' WITHOUT CHANGE OF COORDINATES AND ANALYTIC DATA' 
  print*, '-----------------------------------------------------'

  print *,'Average error in nodes (per-per) without change of coordinates='&
       ,acc1/(npts1*npts2), ',  initialization time (s): ', t1i, &
       ',  solution time (s): ', t1e,'Norm L2',sqrt(normL2_1),'Norm H1',sqrt(normH1_1)
  print *,'Average error in nodes (per-dir) without change of coordinates='&
       ,acc2/(npts1*npts2), ',  initialization time (s): ', t2i, &
       ',  solution time (s): ', t2e,'Norm L2',sqrt(normL2_2),'Norm H1',sqrt(normH1_2)
  print *,'Average error in nodes (dir-dir) without change of coordinates='&
       ,acc3/(npts1*npts2), ',  initialization time (s): ', t3i, &
       ',  solution time (s): ', t3e,'Norm L2',sqrt(normL2_3),'Norm H1',sqrt(normH1_3)
  print *,'Average error in nodes (dir-per) without change of coordinates='&
       ,acc4/(npts1*npts2), ',  initialization time (s): ', t4i, &
       ',  solution time (s): ', t4e,'Norm L2',sqrt(normL2_4),'Norm H1',sqrt(normH1_4)
!!$
  print*, '-------------------------------------------------------'
  print*, ' COLELLA CHANGE OF COORDINATES AND ANALYTIC DATA' 
  print*, '-------------------------------------------------------'
  print *,'Average error in nodes (per-per) '
  print*, 'with colella change of coordinates='&
       ,acc5/(npts1*npts2), ',  initialization time (s): ', t5i, &
       ',  solution time (s): ', t5e,'Norm L2',sqrt(normL2_5),'Norm H1',sqrt(normH1_5)
  print *,'Average error in nodes (per-dir) '
  print*, 'with colella change of coordinates='&
       ,acc6/(npts1*npts2), ',  initialization time (s): ', t6i, &
       ',  solution time (s): ', t6e,'Norm L2',sqrt(normL2_6),'Norm H1',sqrt(normH1_6)
  print *,'Average error in nodes (dir-dir) '
  print*, 'with colella change of coordinates='&
       ,acc7/(npts1*npts2), ',  initialization time (s): ', t7i, &
       ',  solution time (s): ', t7e,'Norm L2',sqrt(normL2_7),'Norm H1',sqrt(normH1_7)
  print *,'Average error in nodes (dir-per) '
  print*, 'with colella change of coordinates='&
       ,acc8/(npts1*npts2), ',  initialization time (s): ', t8i, &
       ',  solution time (s): ', t8e,'Norm L2',sqrt(normL2_8),'Norm H1',sqrt(normH1_8)
!!$  
  print*, '-------------------------------------------------------'
  print*, ' WITHOUT CHANGE OF COORDINATES AND WITH A SOURCE TERM NON-ANALYTIC' 
  print*, '-------------------------------------------------------'
  print *,'Average error in nodes (per-per) '
  print*, 'without change of coordinates='&
       ,acc95/(npts1*npts2), ',  initialization time (s): ', t95i, &
       ',  solution time (s): ', t95e,'Norm L2',sqrt(normL2_95),'Norm H1',sqrt(normH1_95),&
       h1**(SPLINE_DEG1),  h1**(SPLINE_DEG1-1)
!!$
!!$
  print*, '-------------------------------------------------------'
  print*, ' COLELLA CHANGE OF COORDINATES AND WITH A SOURCE TERM NON-ANALYTIC' 
  print*, '-------------------------------------------------------'
  print *,'Average error in nodes (per-per) '
  print*, 'with colella change of coordinates and source term non analytic='&
       ,acc9/(npts1*npts2), ',  initialization time (s): ', t9i, &
       ',  solution time (s): ', t9e,'Norm L2',sqrt(normL2_9),'Norm H1',sqrt(normH1_9)
  print *,'Average error in nodes (per-dir) '
  print*, 'with colella change of coordinates and source term non analytic='&
       ,acc10/(npts1*npts2), ',  initialization time (s): ', t10i, &
       ',  solution time (s): ', t10e,'Norm L2',sqrt(normL2_10),'Norm H1',sqrt(normH1_10)
  print *,'Average error in nodes (dir-dir)'
  print*, 'with colella change of coordinates and source term non analytic='&
       ,acc11/(npts1*npts2), ',  initialization time (s): ', t11i, &
       ',  solution time (s): ', t11e,'Norm L2',sqrt(normL2_11),'Norm H1',sqrt(normH1_11)
  print *,'Average error in nodes (dir-per) '
  print*, 'with colella change of coordinates and source term non analytic='&
       ,acc12/(npts1*npts2), ',  initialization time (s): ', t12i, &
       ',  solution time (s): ', t12e,'Norm L2',sqrt(normL2_12),'Norm H1',sqrt(normH1_12)
  
!!$  
!print*,h1**(SPLINE_DEG1-2)*h2**(SPLINE_DEG2-2)
!borne_L2 = 1.8*sll_pi**2*h1**(SPLINE_DEG1-1)*h2**(SPLINE_DEG2-1)
  if ( ( sqrt(normL2_1) <= h1**(SPLINE_DEG1-1))   .AND. &
       ( sqrt(normL2_2) <= h1**(SPLINE_DEG1-1))   .AND. &
       ( sqrt(normL2_3) <= h1**(SPLINE_DEG1 -1))   .AND. &
       ( sqrt(normL2_4) <= h1**(SPLINE_DEG1 -1))   .AND. &
       ( sqrt(normL2_5) <= h1**(SPLINE_DEG1 -1))   .AND. &
       ( sqrt(normL2_6) <= h1**(SPLINE_DEG1 -1))   .AND. &
       ( sqrt(normL2_7) <= h1**(SPLINE_DEG1 -1))   .AND. &
       ( sqrt(normL2_8) <= h1**(SPLINE_DEG1 -1))   .AND. &
       ( sqrt(normL2_95)<= h1**(SPLINE_DEG1 -1))   .AND. &
       ( sqrt(normL2_9) <= h1**(SPLINE_DEG1 -1))   .AND. &
       ( sqrt(normL2_10)<= h1**(SPLINE_DEG1 -1))   .AND. &
       ( sqrt(normL2_11)<= h1**(SPLINE_DEG1 -1))   .AND. &
       ( sqrt(normL2_12)<= h1**(SPLINE_DEG1 -1))   .AND. &
       ( sqrt(normH1_1) <= h1**(SPLINE_DEG1 -1-1))   .AND. &
       ( sqrt(normH1_2) <= h1**(SPLINE_DEG1 -1-1))   .AND. &
       ( sqrt(normH1_3) <= h1**(SPLINE_DEG1 -1-1))   .AND. &
       ( sqrt(normH1_4) <= h1**(SPLINE_DEG1 -1-1))   .AND. &
       ( sqrt(normH1_5) <= h1**(SPLINE_DEG1 -1-1))   .AND. &
       ( sqrt(normH1_6) <= h1**(SPLINE_DEG1 -1-1))   .AND. &
       ( sqrt(normH1_7) <= h1**(SPLINE_DEG1 -1-1))   .AND. &
       ( sqrt(normH1_8) <= h1**(SPLINE_DEG1 -1-1))   .AND. &
       ( sqrt(normH1_95)<= h1**(SPLINE_DEG1 -1-1))   .AND. &
       ( sqrt(normH1_9) <= h1**(SPLINE_DEG1 -1-1))   .AND. &
       ( sqrt(normH1_10)<= h1**(SPLINE_DEG1 -1-1))   .AND. &
       ( sqrt(normH1_11)<= h1**(SPLINE_DEG1 -1-1))   .AND. &
       ( sqrt(normH1_12)<= h1**(SPLINE_DEG1 -1-1)) ) then
     
       
     print *, 'PASSED'
  end if
end program test_general_elliptic_solver




! External functions used as parameters in the above unit test:


function func_one( eta1, eta2, params ) result(res)
  real(8), intent(in) :: eta1
  real(8), intent(in) :: eta2
  real(8), dimension(:), intent(in) :: params
  real(8) :: res
  res = 1.0_8
end function func_one

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

  res = 0.0_8
end function func_epsi


!----------------------------------------------------------
!  Solution for a identity change of coordinates 
!   and periodic-periodic conditions
!   the matrix A is equal to identity 
!   the scalar c is equal to zero 
!-------------------------------------------------------------
function source_term_perper( eta1, eta2) result(res)
  use sll_constants
  intrinsic :: cos

  real(8), intent(in) :: eta1
  real(8), intent(in) :: eta2
  ! real(8), dimension(:), intent(in), optional :: params
  real(8) :: res

  res =  0.001*cos(2*sll_pi*eta1)
  !!-2*(2.0*sll_pi)**2*cos(2.0*sll_pi*eta1)*cos(2.0*sll_pi*eta2)! 0.001*cos(2*sll_pi*eta1)!
end function source_term_perper

real(8) function sol_exacte_perper(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  
  !real(8), dimension(:), intent(in), optional :: params
  sol_exacte_perper = -0.001/((2*sll_pi)**2)*cos(2*sll_pi*eta1)!cos(2.0*sll_pi*eta1)*cos(2.0*sll_pi*eta2)!-0.001/((2*sll_pi)**2)*cos(2*sll_pi*eta1)
end function sol_exacte_perper

real(8) function sol_exacte_perper_der1(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  
  !real(8), dimension(:), intent(in), optional :: params
  sol_exacte_perper_der1 = 0.001/(2*sll_pi)*sin(2*sll_pi*eta1) !-2.0*sll_pi*sin(2.0*sll_pi*eta1)*cos(2.0*sll_pi*eta2)
end function sol_exacte_perper_der1
real(8) function sol_exacte_perper_der2(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  
  !real(8), dimension(:), intent(in), optional :: params
  sol_exacte_perper_der2 = 0.0_f64!-2.0*sll_pi*cos(2.0*sll_pi*eta1)*sin(2.0*sll_pi*eta2)
end function sol_exacte_perper_der2

!----------------------------------------------------------
!  Solution for a identity change of coordinates 
!   and periodic-dirichlet conditions
!   and also dirichlet-dirichlet conditons
!   the matrix A is equal to identity 
!   the scalar c is equal to zero 
!-------------------------------------------------------------

real(8) function source_term_perdir(eta1,eta2,params) ! in the path
  use sll_constants
  intrinsic :: cos
  intrinsic :: sin 
  real(8),intent(in) :: eta1,eta2
  real(8), dimension(:), intent(in), optional :: params

  
  source_term_perdir = -2*(2*sll_pi)**2* sin(2*sll_pi*eta1)*sin(2*sll_pi*eta2)
      ! -(16.0*sll_pi**2*eta2**4 &
      ! - 16.0*sll_pi**2*eta2**2 &
      ! - 12.0*eta2**2 + 2.0)*cos(2*sll_pi*eta1)*sin(2*sll_pi*eta1)
  
end function source_term_perdir


real(8) function sol_exacte_perdir(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  intrinsic :: cos
  intrinsic :: sin
  !real(8), dimension(:), intent(in), optional :: params
  sol_exacte_perdir = sin(2.0*sll_pi*eta1)*sin(2.0*sll_pi*eta2)!eta2 ** 2 * (eta2**2-1)&
      ! * cos(2.0*sll_pi*eta1)*sin(2.0*sll_pi*eta1)
  
  !print*, 'heho'
end function sol_exacte_perdir


real(8) function sol_exacte_perdir_der1(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  intrinsic :: cos
  intrinsic :: sin
  !real(8), dimension(:), intent(in), optional :: params
  sol_exacte_perdir_der1 = 2.0*sll_pi*cos(2.0*sll_pi*eta1)*sin(2.0*sll_pi*eta2)
end function sol_exacte_perdir_der1


real(8) function sol_exacte_perdir_der2(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  intrinsic :: cos
  intrinsic :: sin
  !real(8), dimension(:), intent(in), optional :: params
  sol_exacte_perdir_der2 = 2.0*sll_pi*sin(2.0*sll_pi*eta1)*cos(2.0*sll_pi*eta2)
end function sol_exacte_perdir_der2

!  Solution for a identity change of coordinates 
!   and also dirichlet-periodicconditons
!   the matrix A is equal to identity 
!   the scalar c is equal to zero 
!-------------------------------------------------------------

real(8) function source_term_dirper(eta1,eta2,params) ! in the path
  use sll_constants
  real(8),intent(in) :: eta1,eta2
  real(8), dimension(:), intent(in), optional :: params

  source_term_dirper = -2*(2*sll_pi)**2* sin(2*sll_pi*eta1)*cos(2*sll_pi*eta2)
     ! -(16.0*sll_pi**2*eta1**4 &
     ! - 16.0*sll_pi**2*eta1**2 &
     ! - 12.0*eta1**2 + 2.0)*sin(2*sll_pi*eta2)*cos(2*sll_pi*eta2)
end function source_term_dirper


real(8) function sol_exacte_dirper(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  !real(8), dimension(:), intent(in), optional :: params
  sol_exacte_dirper = sin(2.0*sll_pi*eta1)*cos(2.0*sll_pi*eta2)
       !eta1 ** 2 * (eta1**2-1)* cos(2*sll_pi*eta2)*sin(2*sll_pi*eta2)
  
  
end function sol_exacte_dirper

real(8) function sol_exacte_dirper_der1(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  !real(8), dimension(:), intent(in), optional :: params
  sol_exacte_dirper_der1 = 2*sll_pi*cos(2.0*sll_pi*eta1)*cos(2.0*sll_pi*eta2)
       !eta1 ** 2 * (eta1**2-1)* cos(2*sll_pi*eta2)*sin(2*sll_pi*eta2)
  
  
end function sol_exacte_dirper_der1

real(8) function sol_exacte_dirper_der2(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  !real(8), dimension(:), intent(in), optional :: params
  sol_exacte_dirper_der2 = -2.0*sll_pi*sin(2.0*sll_pi*eta1)*sin(2.0*sll_pi*eta2)
       !eta1 ** 2 * (eta1**2-1)* cos(2*sll_pi*eta2)*sin(2*sll_pi*eta2)
  
  
end function sol_exacte_dirper_der2

!----------------------------------------------------------
!  Solution for a r theta change of coordinates 
!   and periodic-dirichlet conditions
!   and also dirivhlet-dirichlet conditons
!   the matrix A is equal to identity 
!   the scalar c is equal to zero 
!-------------------------------------------------------------


real(8) function rho_rtheta(eta1,eta2,params) ! in the path
  use sll_constants
  intrinsic :: cos
  intrinsic :: sin
  real(8),intent(in) :: eta1,eta2
  real(8) :: x, y
  real(8), dimension(:), intent(in), optional :: params
  
  if (present(params)) print*, params

  x = eta2*cos(2*sll_pi*eta1)
  y = eta2*sin(2*sll_pi*eta1)
  
  rho_rtheta = x*y*(-32.0*x**2 - 32.0*y**2 + 15.0)  
  
  
end function rho_rtheta


real(8) function sol_exacte_rtheta(eta1,eta2,params) ! in the path
  use sll_constants
  real(8),intent(in) :: eta1,eta2
  intrinsic :: cos
  intrinsic :: sin
  real(8), dimension(:), intent(in), optional :: params
  
  if (present(params)) print*, params
  
  sol_exacte_rtheta = ( eta2**2-1)*(eta2**2-0.5**2)*eta2**2&
       *cos(2*sll_pi*eta1)*sin(2*sll_pi*eta1)
  
    
end function sol_exacte_rtheta


!----------------------------------------------------------
!  Solution for a colella change of coordinates 
!   and periodic-periodic conditons
!   the matrix A is equal to identity 
!   the scalar c is equal to zero 
!-------------------------------------------------------------




real(8) function source_term_chgt_perper(eta1,eta2) ! in the path
  use sll_constants
  intrinsic :: cos
  intrinsic :: sin
  real(8):: eta1,eta2
  real(8) :: x, y
  ! real(8), dimension(:), intent(in), optional :: params
  
  x =   eta1 + 0.1_f64*sin(2*sll_pi*eta1) * sin(2*sll_pi*eta2)
  y =   eta2 + 0.1_f64*sin(2*sll_pi*eta1) * sin(2*sll_pi*eta2)
  
  source_term_chgt_perper = -8.0*sll_pi**2*cos(2*sll_pi*x)*cos(2*sll_pi*y) 
  
end function source_term_chgt_perper




real(8) function sol_exacte_chgt_perper(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  real(8) :: x,y
  intrinsic :: cos
  intrinsic :: sin
  
  x =   eta1 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  y =   eta2 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  
  sol_exacte_chgt_perper = cos(2*sll_pi*x)*cos(2*sll_pi*y)
  
  
end function sol_exacte_chgt_perper

real(8) function sol_exacte_chgt_perper_der1(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  real(8) :: x,y
  intrinsic :: cos
  intrinsic :: sin
  
  x =   eta1 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  y =   eta2 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  
  sol_exacte_chgt_perper_der1 = -2*sll_pi*sin(2*sll_pi*x)*cos(2*sll_pi*y)&
       * ( 1.0_f64 + 0.1*2*sll_pi*cos(2* sll_pi*eta1) * sin(2*sll_pi*eta2) )&
       -2*sll_pi*cos(2*sll_pi*x)*sin(2*sll_pi*y)&
       * ( 0.1*2*sll_pi*cos(2* sll_pi*eta1) * sin(2*sll_pi*eta2) )
  
  
end function sol_exacte_chgt_perper_der1

real(8) function sol_exacte_chgt_perper_der2(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  real(8) :: x,y
  intrinsic :: cos
  intrinsic :: sin
  
  x =   eta1 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  y =   eta2 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  
  !sol_exacte_chgt_perper_der2 = -2*sll_pi*cos(2*sll_pi*x)*sin(2*sll_pi*y)
  
  sol_exacte_chgt_perper_der2 = -2*sll_pi*sin(2*sll_pi*x)*cos(2*sll_pi*y)&
       * ( 0.1*2*sll_pi*sin(2* sll_pi*eta1) * cos(2*sll_pi*eta2) )&
       -2*sll_pi*cos(2*sll_pi*x)*sin(2*sll_pi*y)&
       * ( 1.0_f64 + 0.1*2*sll_pi*sin(2* sll_pi*eta1)*cos(2*sll_pi*eta2) )
end function sol_exacte_chgt_perper_der2


!----------------------------------------------------------
!  Solution for a colella change of coordinates 
!   and periodic-dirichlet conditons
!   and dircihlet-diichlet conditions
!   the matrix A is equal to identity 
!   the scalar c is equal to zero 
!-------------------------------------------------------------




real(8) function source_term_chgt_perdir(eta1,eta2) ! in the path
  use sll_constants
  real(8),intent(in) :: eta1,eta2
  real(8) :: x, y
  intrinsic :: cos
  intrinsic :: sin
  !real(8), dimension(:), intent(in), optional :: params
    
  x =   eta1 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  y =   eta2 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  
  source_term_chgt_perdir= -2*(2*sll_pi)**2 * sin(2*sll_pi*y)*cos(2*sll_pi*x)
  
  
end function source_term_chgt_perdir



real(8) function sol_exacte_chgt_perdir(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  real(8) :: x,y
  intrinsic :: cos
  intrinsic :: sin
 ! real(8), dimension(:), intent(in), optional :: params
  x =   eta1 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  y =   eta2 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  
  sol_exacte_chgt_perdir = cos(2*sll_pi*x)*sin(2*sll_pi*y)
  
  
end function sol_exacte_chgt_perdir

real(8) function sol_exacte_chgt_perdir_der1(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  real(8) :: x,y
  intrinsic :: cos
  intrinsic :: sin
 ! real(8), dimension(:), intent(in), optional :: params
  x =   eta1 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  y =   eta2 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  
  sol_exacte_chgt_perdir_der1 = -2*sll_pi*sin(2*sll_pi*x)*sin(2*sll_pi*y)&
       * ( 1.0_f64 + 0.1*2*sll_pi*cos(2*sll_pi*eta1) * sin(2*sll_pi*eta2) )&
       + 2*sll_pi*cos(2*sll_pi*x)*cos(2*sll_pi*y)&
       * ( 2*sll_pi*0.1*cos(2* sll_pi*eta1) * sin(2*sll_pi*eta2) ) 
  
  
end function sol_exacte_chgt_perdir_der1

real(8) function sol_exacte_chgt_perdir_der2(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  real(8) :: x,y
  intrinsic :: cos
  intrinsic :: sin
 ! real(8), dimensi@on(:), intent(in), optional :: params
  x =   eta1 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  y =   eta2 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  
  sol_exacte_chgt_perdir_der2 = -2*sll_pi*sin(2*sll_pi*x)*sin(2*sll_pi*y)&
       * ( 0.1*2*sll_pi*sin(2*sll_pi*eta1) * cos(2*sll_pi*eta2) ) &
       + 2*sll_pi*cos(2*sll_pi*x)*cos(2*sll_pi*y)&
       * ( 1.0_f64 + 2*sll_pi*0.1*sin(2* sll_pi*eta1) *cos(2*sll_pi*eta2) ) 
  
  
end function sol_exacte_chgt_perdir_der2

!----------------------------------------------------------
!  Solution for a colella change of coordinates 
!   and dirchlet-periodic conditions
!   the matrix A is equal to identity 
!   the scalar c is equal to zero 
!-------------------------------------------------------------



real(8) function source_term_chgt_dirdir(eta1,eta2) ! in the path
  use sll_constants
  real(8),intent(in) :: eta1,eta2
  real(8) :: x, y
  intrinsic :: cos
  intrinsic :: sin
  ! -------------------------------------------------
  ! In the case without change of coordinates
  ! -------------------------------------------------
  x =   eta1 + 0.1_8*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  y =   eta2 + 0.1_8*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  
  
  source_term_chgt_dirdir = &
       -2*(2.0*sll_pi)**2*sin(2*sll_pi*x)*sin(2*sll_pi*y)
  
end function source_term_chgt_dirdir




real(8) function sol_exacte_chgt_dirdir(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  real(8) :: x,y
  intrinsic :: cos
  intrinsic :: sin
  
  x =   eta1 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  y =   eta2 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  
  
  sol_exacte_chgt_dirdir = sin(2* sll_pi*y)*sin(2* sll_pi*x)
  
end function sol_exacte_chgt_dirdir


real(8) function sol_exacte_chgt_dirdir_der1(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  real(8) :: x,y
  intrinsic :: cos
  intrinsic :: sin
  
  x =   eta1 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  y =   eta2 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  
  
  sol_exacte_chgt_dirdir_der1 = 2*sll_pi*cos(2* sll_pi*x)*sin(2* sll_pi*y)&
       * ( 1.0_f64 + 0.1*2*sll_pi*cos(2*sll_pi*eta1) * sin(2*sll_pi*eta2) )&
       + 2*sll_pi*sin(2* sll_pi*x)*cos(2* sll_pi*y) &
       * ( 2*sll_pi*0.1*cos(2* sll_pi*eta1) * sin(2*sll_pi*eta2) )
end function sol_exacte_chgt_dirdir_der1


real(8) function sol_exacte_chgt_dirdir_der2(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  real(8) :: x,y
  intrinsic :: cos
  intrinsic :: sin
  
  x =   eta1 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  y =   eta2 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  
  
  sol_exacte_chgt_dirdir_der2 =  2*sll_pi*cos(2* sll_pi*x)*sin(2* sll_pi*y)&
       * ( 0.1*2*sll_pi*sin(2*sll_pi*eta1) * cos(2*sll_pi*eta2)  )&
       + 2*sll_pi*sin(2* sll_pi*x)*cos(2* sll_pi*y) &
       * ( 1.0_f64 + 2*sll_pi*0.1*sin(2* sll_pi*eta1) *cos(2*sll_pi*eta2) )
  
end function sol_exacte_chgt_dirdir_der2





!----------------------------------------------------------
!  Solution for a colella change of coordinates 
!   and dirchlet-periodic conditions
!   the matrix A is equal to identity 
!   the scalar c is equal to zero 
!-------------------------------------------------------------



real(8) function source_term_chgt_dirper(eta1,eta2) ! in the path
  use sll_constants
  real(8),intent(in) :: eta1,eta2
  real(8) :: x, y
  intrinsic :: cos
  intrinsic :: sin
  !real(8), dimension(:), intent(in), optional :: params
  ! -------------------------------------------------
  ! In the case without change of coordinates
  ! -------------------------------------------------
  x =   eta1 + 0.1_8*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  y =   eta2 + 0.1_8*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  
  
  source_term_chgt_dirper = -2*(2*sll_pi)**2*sin(2*sll_pi*x)*cos(2*sll_pi*y)
  
end function source_term_chgt_dirper




real(8) function sol_exacte_chgt_dirper(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  real(8) :: x,y
  intrinsic :: cos
  intrinsic :: sin

  x =   eta1 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  y =   eta2 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  
  
  sol_exacte_chgt_dirper = sin(2* sll_pi*x)*cos(2* sll_pi*y)
  
end function sol_exacte_chgt_dirper


real(8) function sol_exacte_chgt_dirper_der1(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  real(8) :: x,y
  intrinsic :: cos
  intrinsic :: sin
  
  x =   eta1 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  y =   eta2 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  
  
  sol_exacte_chgt_dirper_der1 = 2*sll_pi*cos(2* sll_pi*x)*cos(2* sll_pi*y) &
       * ( 1.0_f64 + 0.1*2*sll_pi*cos(2*sll_pi*eta1) * sin(2*sll_pi*eta2) )&
       - 2*sll_pi*sin(2* sll_pi*x)*sin(2* sll_pi*y)&
       * ( 2*sll_pi*0.1*cos(2* sll_pi*eta1) * sin(2*sll_pi*eta2) ) 
end function sol_exacte_chgt_dirper_der1


real(8) function sol_exacte_chgt_dirper_der2(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  real(8) :: x,y
  intrinsic :: cos
  intrinsic :: sin
  
  x =   eta1 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  y =   eta2 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  
  
  sol_exacte_chgt_dirper_der2 = 2*sll_pi*cos(2* sll_pi*x)*cos(2* sll_pi*y) &
       * ( 0.1*2*sll_pi*sin(2*sll_pi*eta1) * cos(2*sll_pi*eta2)  )&
       - 2*sll_pi*sin(2* sll_pi*x)*sin(2* sll_pi*y)&
       * (1.0_f64 + 2*sll_pi*0.1*sin(2* sll_pi*eta1) *cos(2*sll_pi*eta2) ) 
  
end function sol_exacte_chgt_dirper_der2


!!!!!! test case with F(theta,phi) = (2pi theta , 2pi phi)

real(8) function adimension_chgt_x(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  print*, eta1, eta2
  adimension_chgt_x = 2*sll_pi*eta1 !+ eta2)
end function adimension_chgt_x

real(8) function adimension_chgt_y(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  print*, eta1, eta2
  adimension_chgt_y = 2*sll_pi*eta2
end function adimension_chgt_y


real(8) function jac11_adimension_chgt(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  print*, eta1, eta2
  jac11_adimension_chgt = 2*sll_pi
end function jac11_adimension_chgt

real(8) function jac12_adimension_chgt(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  print*, eta1, eta2
  jac12_adimension_chgt = 0.0!sll_pi
end function jac12_adimension_chgt

real(8) function jac21_adimension_chgt(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  print*, eta1, eta2
  jac21_adimension_chgt = 0.0!2*sll_pi!0.0
end function jac21_adimension_chgt

real(8) function jac22_adimension_chgt(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  print*, eta1, eta2
  jac22_adimension_chgt = 2*sll_pi
end function jac22_adimension_chgt



real(8) function sol_exacte_chgt_adim(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  real(8) :: x,y
  intrinsic :: cos
  intrinsic :: sin
  
  x =   2*sll_pi*eta1!+eta2)
  y =   2* sll_pi*eta2
  
  
  sol_exacte_chgt_adim = cos(x)*cos(y)
  
end function sol_exacte_chgt_adim


real(8) function source_term_chgt_adim(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  real(8) :: x,y
  intrinsic :: cos
  intrinsic :: sin
  
  x =   2*sll_pi*eta1 !+eta2)
  y =   2* sll_pi*eta2
  
  
  source_term_chgt_adim = -2*cos(x)*cos(y)
  
end function source_term_chgt_adim
