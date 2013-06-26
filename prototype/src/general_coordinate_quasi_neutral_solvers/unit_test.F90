program test_general_qns
  use sll_logical_meshes
  use sll_module_coordinate_transformations_2d
  use sll_common_coordinate_transformations
  use sll_general_coordinate_qn_solver_module
  use sll_module_scalar_field_2d_alternative
  use sll_constants
  use sll_arbitrary_degree_spline_interpolator_2d_module
  use sll_timer
#include "sll_memory.h"
#include "sll_working_precision.h"
  implicit none

#define SPLINE_DEG1 3
#define SPLINE_DEG2 3
#define NUM_CELLS1  30
#define NUM_CELLS2  30
#define ETA1MIN  0.0_f64
#define ETA1MAX  1.0_f64
#define ETA2MIN  0.0_f64
#define ETA2MAX  1.0_f64
#define PRINT_COMPARISON .false.

  type(sll_logical_mesh_2d), pointer                    :: mesh_2d
  class(sll_coordinate_transformation_2d_base), pointer :: T
  type(general_coordinate_qn_solver)                    :: qns
  type(arb_deg_2d_interpolator), target                 :: interp_2d
  type(arb_deg_2d_interpolator), target                 :: interp_2d_term_source
 ! class(sll_interpolator_2d_base), pointer              :: interp_2d_ptr
  class(sll_interpolator_2d_base), pointer              :: terme_source_interp
!  class(sll_scalar_field_2d_analytic_alt), dimension(2,2) :: a_field_mat
  type(sll_scalar_field_2d_base_ptr), dimension(2,2)    :: a_field_mat
  class(sll_scalar_field_2d_base), pointer              :: c_field
  class(sll_scalar_field_2d_base), pointer              :: rho
  type(sll_scalar_field_2d_discrete_alt), pointer       :: phi
  type(sll_time_mark) :: t_reference
  sll_real64 :: t1i, t1e, t2i, t2e, t3i, t3e, t4i, t4e, t5i, t5e, t6i, t6e, &
       t7i, t7e, t8i, t8e, t9i,t9e,t10i,t10e,t11i,t11e,t12i,t12e,t95e,t95i
 ! sll_real64 :: t105e,t105i,t115e,t115i,t125i,t125e
  real(8), external :: func_zero
  real(8), external :: func_one
  real(8), external :: source_term_perper
  real(8), external :: source_term_perdir
  real(8), external :: source_term_dirper
  real(8), external :: source_term_chgt_perper
  real(8), external :: source_term_chgt_perdir
  real(8), external :: source_term_chgt_dirper
  real(8), external :: source_term_chgt_dirdir
  sll_real64, dimension(:,:), allocatable :: values
  sll_real64 :: acc1,acc2,acc3,acc4,acc5,acc6,acc7,acc8,acc9
  sll_real64 :: acc10,acc11,acc12,acc95
 ! sll_real64 :: acc105,acc115,acc125
  sll_real64, dimension(:,:), allocatable    :: calculated
  sll_real64, dimension(:,:), allocatable    :: difference
  sll_real64, dimension(:,:), allocatable    :: tab_rho
  sll_real64, dimension(:),   allocatable    :: point1
  sll_real64, dimension(:),   allocatable    :: point2
!  sll_real64, dimension(:,:), pointer        :: test_coeff
  sll_int32 :: ierr
  sll_int32  :: i, j
  sll_real64 :: h1,h2,eta1,eta2,node_val,ref
  sll_int32 :: npts1,npts2
  real(8), external :: sol_exacte_perper
  real(8), external :: sol_exacte_perdir  
  real(8), external :: sol_exacte_dirper
  real(8), external :: sol_exacte_chgt_perper
  real(8), external :: sol_exacte_chgt_perdir  
  real(8), external :: sol_exacte_chgt_dirper
  real(8), external :: sol_exacte_chgt_dirdir
 ! sll_real64 :: node_val1
  !sll_real64 :: epsi
  !sll_real64 :: epsi1

 ! epsi  =  0.000_f64
 ! epsi1 =  0.000_f64 ! penalization method
  
  
  
    
  !--------------------------------------------------------------------
  !     1 test case without chane of coordinates 
  !      periodic-periodic boundary conditions
  !--------------------------------------------------------------------
  
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


  c_field => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "c_field", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC )

  rho => new_scalar_field_2d_analytic_alt( &
       source_term_perper, &
       "rho1", &     
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

 ! interp_2d_ptr => interp_2d

  phi => new_scalar_field_2d_discrete_alt( &
       values, &
       "phi1", &
       interp_2d, &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC )

  print *, 'initialized fields...'

  call set_time_mark(t_reference)

  call initialize_general_qn_solver( &
       qns, &
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
       ETA2MAX)
 
  t1i = time_elapsed_since(t_reference)
 
  print *, 'Initialized QNS object'

  call set_time_mark(t_reference)

  ! solve the field
  call solve_quasi_neutral_eq_general_coords( &
       qns, &
       a_field_mat, &
       c_field, &
       rho, &
       phi )

  t1e = time_elapsed_since(t_reference)

  !print *, 'Completed solution',qns%phi_vec
  print*, 'reorganizaton of splines coefficients of solution'
  
!  print *, 'Compare the values of the transformation at the nodes: '
  
  acc1 = 0.0_f64
  
  do j=0,npts2-1
     do i=0,npts1-1
        eta1       = real(i,f64)*h1 + ETA1MIN
        eta2       = real(j,f64)*h2 + ETA2MIN
        node_val   = phi%value_at_point(eta1,eta2)
        ref        = sol_exacte_perper(eta1,eta2)
        calculated(i+1,j+1) = node_val
        difference(i+1,j+1) = ref-node_val
        if(PRINT_COMPARISON) then
        print *, '(eta1,eta2) = ', eta1, eta2, 'calculated = ', node_val, &
             'theoretical = ', ref
        end if
        acc1        = acc1 + abs(node_val-ref)
     end do
  end do

  call phi%write_to_file(0)


  ! delete things...
  call delete(qns)
  call rho%delete()
  call c_field%delete()
  call phi%delete()

  call a_field_mat(1,1)%base%delete()
  call a_field_mat(1,2)%base%delete()
  call a_field_mat(2,1)%base%delete()
  call a_field_mat(2,2)%base%delete()

  call T%delete()
  
  DEALLOCATE(values)
  DEALLOCATE(calculated)
  DEALLOCATE(difference)

  !--------------------------------------------------------------------
  
  !     2 test case without chane of coordinates 
  !      periodic-dirichlet boundary conditions
  
  !--------------------------------------------------------------------
  
  
  
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
       SLL_DIRICHLET, &
       SLL_DIRICHLET ) 

  a_field_mat(1,2)%base => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "a12", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET ) 

  a_field_mat(2,1)%base => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "a21", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET ) 
  
  a_field_mat(2,2)%base => new_scalar_field_2d_analytic_alt( &
       func_one, &
       "a22", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET ) 
  
  
  c_field => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "c_field", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET )

  rho => new_scalar_field_2d_analytic_alt( &
       source_term_perdir, &
       "rho2", &     
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET )

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
       values, &
       "phi2", &
       interp_2d, &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET)
  
  print *, 'initialized fields...'
!  print *, 'a = ', qns%csr_mat%opr_a

  call set_time_mark(t_reference)

  call initialize_general_qn_solver( &
       qns, &
       SPLINE_DEG1, &
       SPLINE_DEG2, &
       NUM_CELLS1, &
       NUM_CELLS2, &
       QNS_GAUSS_LEGENDRE, &
       QNS_GAUSS_LEGENDRE, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       ETA1MIN, &
       ETA1MAX, &
       ETA2MIN, &
       ETA2MAX)
  t2i = time_elapsed_since(t_reference) 
  print *, 'Initialized QNS object'
  
  call set_time_mark(t_reference)

  call solve_quasi_neutral_eq_general_coords( &
       qns, &
       a_field_mat, &
       c_field, &
       rho, &
       phi )

  !print *, 'Completed solution',qns%phi_vec

  t2e = time_elapsed_since(t_reference)
  
!  print *, 'Compare the values of the transformation at the nodes: '
  
  acc2 = 0.0_f64
  
  do j=0,npts2-1
     do i=0,npts1-1
        eta1       = real(i,f64)*h1 + ETA1MIN
        eta2       = real(j,f64)*h2 + ETA2MIN
        node_val   = phi%value_at_point(eta1,eta2)
        ref        = sol_exacte_perdir(eta1,eta2)
!        print*,sin(2*sll_pi*eta1)*cos(2*sll_pi*eta1)
        calculated(i+1,j+1) = node_val
        difference(i+1,j+1) = ref-node_val
        if(PRINT_COMPARISON) then
           print *, '(eta1,eta2) = ', eta1, eta2, 'calculated = ', node_val, &
                'theoretical = ', ref
        end if
        acc2        = acc2 + abs(node_val-ref)
     end do
  end do

  call phi%write_to_file(0)
    
  ! delete things...
  call delete(qns)
  call rho%delete()
  call c_field%delete()
  call phi%delete()
  call a_field_mat(1,1)%base%delete()
  call a_field_mat(1,2)%base%delete()
  call a_field_mat(2,1)%base%delete()
  call a_field_mat(2,2)%base%delete()

  call T%delete()
  
  DEALLOCATE(values)
  DEALLOCATE(calculated)
  DEALLOCATE(difference)
  
  !--------------------------------------------------------------------
  
  !     3 test case without chane of coordinates 
  !      dirichlet-dirichlet boundary conditions
  
  !--------------------------------------------------------------------
  
  
  
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
       identity_jac22 )
  print *, 'initialized coordinate transformation'
  
  ! Thirdly, each field object must be initialized using the same logical
  ! mesh and coordinate transformation.
  a_field_mat(1,1)%base => new_scalar_field_2d_analytic_alt( &
       func_one, &
       "a11", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET ) 
  
  a_field_mat(1,2)%base => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "a12", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET ) 
  
  a_field_mat(2,1)%base => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "a21", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET ) 
  
  a_field_mat(2,2)%base => new_scalar_field_2d_analytic_alt( &
       func_one, &
       "a22", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET ) 
  
  
  c_field => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "c_field", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET )

  rho => new_scalar_field_2d_analytic_alt( &
       source_term_perdir, &
       "rho3", &     
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET )
  
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
       values, &
       "phi3", &
       interp_2d, &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET)
  
  print *, 'initialized fields...'
  
  call set_time_mark(t_reference)

  call initialize_general_qn_solver( &
       qns, &
       SPLINE_DEG1, &
       SPLINE_DEG2, &
       NUM_CELLS1, &
       NUM_CELLS2, &
       QNS_GAUSS_LEGENDRE, &
       QNS_GAUSS_LEGENDRE, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       ETA1MIN, &
       ETA1MAX, &
       ETA2MIN, &
       ETA2MAX)
  
  t3i = time_elapsed_since(t_reference) 

  print *, 'Initialized QNS object'

  call set_time_mark(t_reference)

  ! solve the field
  call solve_quasi_neutral_eq_general_coords( &
       qns, &
       a_field_mat, &
       c_field, &
       rho, &
       phi )

  call interp_2d%set_coefficients( qns%phi_vec)
 
  
  t3e = time_elapsed_since(t_reference)


  acc3 = 0.0_f64
  
  do j=0,npts2-1
     do i=0,npts1-1
        eta1       = real(i,f64)*h1 + ETA1MIN
        eta2       = real(j,f64)*h2 + ETA2MIN

        node_val   =phi%value_at_point(eta1,eta2)

        ref        = sol_exacte_perdir(eta1,eta2)
        calculated(i+1,j+1) = node_val
        difference(i+1,j+1) = ref-node_val
        if(PRINT_COMPARISON) then
           print *, '(eta1,eta2) = ', eta1, eta2, 'calculated = ', node_val, &
                'theoretical = ', ref
        end if
        acc3        = acc3 + abs(node_val-ref)
     end do
  end do
  
  call phi%write_to_file(0)
  ! delete things...
  call delete(qns)
  call rho%delete()
  call c_field%delete()
  call phi%delete()
  call a_field_mat(1,1)%base%delete()
  call a_field_mat(1,2)%base%delete()
  call a_field_mat(2,1)%base%delete()
  call a_field_mat(2,2)%base%delete()

  call T%delete()
  
  DEALLOCATE(values)
  DEALLOCATE(calculated)
  DEALLOCATE(difference)


  
  !--------------------------------------------------------------------
  
  !     4 test case without change of coordinates 
  !     dirichlet-periodic boundary conditions
  
  !--------------------------------------------------------------------
    
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
       identity_jac22 )
  print *, 'initialized coordinate transformation'
  
  ! Thirdly, each field object must be initialized using the same logical
  ! mesh and coordinate transformation.
  a_field_mat(1,1)%base => new_scalar_field_2d_analytic_alt( &
       func_one, &
       "a11", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_PERIODIC, &
       SLL_PERIODIC ) 
  
  a_field_mat(1,2)%base => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "a12", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_PERIODIC, &
       SLL_PERIODIC ) 
  
  a_field_mat(2,1)%base => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "a21", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_PERIODIC, &
       SLL_PERIODIC ) 
  
  a_field_mat(2,2)%base => new_scalar_field_2d_analytic_alt( &
       func_one, &
       "a22", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_PERIODIC, &
       SLL_PERIODIC ) 
  
  c_field => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "c_field", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_PERIODIC, &
       SLL_PERIODIC )

  rho => new_scalar_field_2d_analytic_alt( &
       source_term_dirper, &
       "rho4", &     
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
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
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SPLINE_DEG1, &
       SPLINE_DEG2 )

  phi => new_scalar_field_2d_discrete_alt( &
       values, &
       "phi4", &
       interp_2d, &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_PERIODIC, &
       SLL_PERIODIC)
  
  print *, 'initialized fields...'

  call set_time_mark(t_reference)

  call initialize_general_qn_solver( &
       qns, &
       SPLINE_DEG1, &
       SPLINE_DEG2, &
       NUM_CELLS1, &
       NUM_CELLS2, &
       QNS_GAUSS_LEGENDRE, &
       QNS_GAUSS_LEGENDRE, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       ETA1MIN, &
       ETA1MAX, &
       ETA2MIN, &
       ETA2MAX)

  t4i = time_elapsed_since(t_reference) 
  print *, 'Initialized QNS object'
  
  call set_time_mark(t_reference)

  ! solve the field
  call solve_quasi_neutral_eq_general_coords( &
       qns, &
       a_field_mat, &
       c_field, &
       rho, &
       phi )

 
  t4e = time_elapsed_since(t_reference) 
  
  
  acc4 = 0.0_f64
  
  do j=0,npts2-1
     do i=0,npts1-1
        eta1       = real(i,f64)*h1 + ETA1MIN
        eta2       = real(j,f64)*h2 + ETA2MIN
        node_val   = phi%value_at_point(eta1,eta2)
        ref        = sol_exacte_dirper(eta1,eta2)
        calculated(i+1,j+1) = node_val
        difference(i+1,j+1) = ref-node_val
        if(PRINT_COMPARISON) then
           print *, '(eta1,eta2) = ', eta1, eta2, 'calculated = ', node_val, &
                'theoretical = ', ref
        end if
        acc4        = acc4 + abs(node_val-ref)
     end do
  end do

  call phi%write_to_file(0)
  ! delete things...
  call delete(qns)
  call rho%delete()
  call c_field%delete()
  call phi%delete()
  call a_field_mat(1,1)%base%delete()
  call a_field_mat(1,2)%base%delete()
  call a_field_mat(2,1)%base%delete()
  call a_field_mat(2,2)%base%delete()

  call T%delete()
  
  DEALLOCATE(values)
  DEALLOCATE(calculated)
  DEALLOCATE(difference)


  !--------------------------------------------------------------------
  
  !     5 test case with colella change of coordinates 
  !     periodic-periodic boundary conditions
  
  !--------------------------------------------------------------------
  
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
       sinprod_jac22)
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
  
  c_field => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "c_field", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC )
  
  rho => new_scalar_field_2d_analytic_alt( &
       source_term_chgt_perper, &
       "rho5", &     
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
  
  
  phi => new_scalar_field_2d_discrete_alt( &
       values, &
       "phi5", &
       interp_2d, &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC)
  
  print *, 'initialized fields...'
  
  call set_time_mark(t_reference)

  call initialize_general_qn_solver( &
       qns, &
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
       ETA2MAX)

  t5i = time_elapsed_since(t_reference) 
  print *, 'Initialized QNS object'

  call set_time_mark(t_reference)  
  ! solve the field
  call solve_quasi_neutral_eq_general_coords( &
       qns, &
       a_field_mat, &
       c_field, &
       rho, &
       phi )

  !print *, 'Completed solution',qns%phi_vec
  
  t5e = time_elapsed_since(t_reference)  
  
!  print *, 'Compare the values of the transformation at the nodes: '
  
  acc5 = 0.0_f64
  
  do j=0,npts2-1
     do i=0,npts1-1
        eta1       = real(i,f64)*h1 + ETA1MIN
        eta2       = real(j,f64)*h2 + ETA2MIN
        node_val   = phi%value_at_point(eta1,eta2)
        !print*, 'rer'
        ref        = sol_exacte_chgt_perper(eta1,eta2)
        calculated(i+1,j+1) = node_val
        difference(i+1,j+1) = ref-node_val
        if(PRINT_COMPARISON) then
           print *, '(eta1,eta2) = ', eta1, eta2, 'calculated = ', node_val, &
                'theoretical = ', ref
        end if
        acc5        = acc5 + abs(node_val-ref)
     end do
  end do

  call phi%write_to_file(0)
  ! delete things...
  call delete(qns)
  call rho%delete()
  call c_field%delete()
  call phi%delete()
  call a_field_mat(1,1)%base%delete()
  call a_field_mat(1,2)%base%delete()
  call a_field_mat(2,1)%base%delete()
  call a_field_mat(2,2)%base%delete()

  call T%delete()
  
  DEALLOCATE(values)
  DEALLOCATE(calculated)
  DEALLOCATE(difference)

  !--------------------------------------------------------------------
  
  !     6 test case with colella change of coordinates 
  !     periodic-dirichlet boundary conditions
  
  !--------------------------------------------------------------------
  
  
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
       sinprod_jac22)
  print *, 'initialized coordinate transformation'
  
  ! Thirdly, each field object must be initialized using the same logical
  ! mesh and coordinate transformation.
  a_field_mat(1,1)%base => new_scalar_field_2d_analytic_alt( &
       func_one, &
       "a11", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET ) 
  
  a_field_mat(1,2)%base => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "a12", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET) 
  
  a_field_mat(2,1)%base => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "a21", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET ) 
  
  a_field_mat(2,2)%base => new_scalar_field_2d_analytic_alt( &
       func_one, &
       "a22", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET ) 
  
  
  c_field => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "c_field", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET )
  
  rho => new_scalar_field_2d_analytic_alt( &
       source_term_chgt_perdir, &
       "rho6", &     
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET )
  
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
       values, &
       "phi6", &
       interp_2d, &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET)
  
  print *, 'initialized fields...'

!  print *, 'a = ', qns%csr_mat%opr_a

  call set_time_mark(t_reference)

  call initialize_general_qn_solver( &
       qns, &
       SPLINE_DEG1, &
       SPLINE_DEG2, &
       NUM_CELLS1, &
       NUM_CELLS2, &
       QNS_GAUSS_LEGENDRE, &
       QNS_GAUSS_LEGENDRE, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       ETA1MIN, &
       ETA1MAX, &
       ETA2MIN, &
       ETA2MAX)

  t6i = time_elapsed_since(t_reference) 

  print *, 'Initialized QNS object'

  call set_time_mark(t_reference)
  
  ! solve the field
  call solve_quasi_neutral_eq_general_coords( &
       qns, &
       a_field_mat, &
       c_field, &
       rho, &
       phi )

  t6e = time_elapsed_since(t_reference)

  
!  print *, 'Compare the values of the transformation at the nodes: '
  
  acc6 = 0.0_f64
  
  do j=0,npts2-1
     do i=0,npts1-1
        eta1       = real(i,f64)*h1 + ETA1MIN
        eta2       = real(j,f64)*h2 + ETA2MIN
        
        node_val   = phi%value_at_point(eta1,eta2)
        !print*, 'rer'
        ref        = sol_exacte_chgt_perdir(eta1,eta2)
        calculated(i+1,j+1) = node_val
        difference(i+1,j+1) = ref-node_val

        if(PRINT_COMPARISON) then
           print *, '(eta1,eta2) = ', eta1, eta2, 'calculated = ', node_val, &
                'theoretical = ', ref
        end if
         acc6        = acc6 + abs(node_val-ref)
     end do
  end do
  call phi%write_to_file(0)
  ! delete things...
  call delete(qns)
  call rho%delete()
  call c_field%delete()
  call phi%delete()
  call a_field_mat(1,1)%base%delete()
  call a_field_mat(1,2)%base%delete()
  call a_field_mat(2,1)%base%delete()
  call a_field_mat(2,2)%base%delete()

  call T%delete()
  
  DEALLOCATE(values)
  DEALLOCATE(calculated)
  DEALLOCATE(difference)



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
       sinprod_jac22)
  print *, 'initialized coordinate transformation'
  
  ! Thirdly, each field object must be initialized using the same logical
  ! mesh and coordinate transformation.
  a_field_mat(1,1)%base => new_scalar_field_2d_analytic_alt( &
       func_one, &
       "a11", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET ) 
  
  a_field_mat(1,2)%base => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "a12", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET) 
  
  a_field_mat(2,1)%base => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "a21", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET ) 
  
  a_field_mat(2,2)%base => new_scalar_field_2d_analytic_alt( &
       func_one, &
       "a22", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET ) 
  
  
  c_field => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "c_field", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET )
  
  rho => new_scalar_field_2d_analytic_alt( &
       source_term_chgt_dirdir, &
       "rho7", &     
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET )
  
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
       values, &
       "phi7", &
       interp_2d, &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET)
  
  print *, 'initialized fields...'

  call set_time_mark(t_reference)
  call initialize_general_qn_solver( &
       qns, &
       SPLINE_DEG1, &
       SPLINE_DEG2, &
       NUM_CELLS1, &
       NUM_CELLS2, &
       QNS_GAUSS_LEGENDRE, &
       QNS_GAUSS_LEGENDRE, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       ETA1MIN, &
       ETA1MAX, &
       ETA2MIN, &
       ETA2MAX)

  t7i = time_elapsed_since(t_reference) 

  print *, 'Initialized QNS object'

  call set_time_mark(t_reference)
  
  ! solve the field
  call solve_quasi_neutral_eq_general_coords( &
       qns, &
       a_field_mat, &
       c_field, &
       rho, &
       phi )
  

  t7e = time_elapsed_since(t_reference)

  
!  print *, 'Compare the values of the transformation at the nodes: '
  
  acc7 = 0.0_f64
  
  do j=0,npts2-1
     do i=0,npts1-1
        eta1       = real(i,f64)*h1 + ETA1MIN
        eta2       = real(j,f64)*h2 + ETA2MIN
        
        node_val   = phi%value_at_point(eta1,eta2)
        !print*, 'rer'
        ref        = sol_exacte_chgt_dirdir(eta1,eta2)
        calculated(i+1,j+1) = node_val
        difference(i+1,j+1) = ref-node_val
        if(PRINT_COMPARISON) then
           print *, '(eta1,eta2) = ', eta1, eta2, 'calculated = ', node_val, &
                'theoretical = ', ref
        end if
        acc7        = acc7 + abs(node_val-ref)
     end do
  end do
  call phi%write_to_file(0)
  ! delete things...
  call delete(qns)
  call rho%delete()
  call c_field%delete()
  call phi%delete()
  call a_field_mat(1,1)%base%delete()
  call a_field_mat(1,2)%base%delete()
  call a_field_mat(2,1)%base%delete()
  call a_field_mat(2,2)%base%delete()

  call T%delete()
  
  DEALLOCATE(values)
  DEALLOCATE(calculated)
  DEALLOCATE(difference)


  !--------------------------------------------------------------------
  
  !     8  test case with colella change of coordinates 
  !     dirichlet-dirichlet boundary conditions
  
  !--------------------------------------------------------------------
  
  
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
       sinprod_jac22)
  print *, 'initialized coordinate transformation'
  
  ! Thirdly, each field object must be initialized using the same logical
  ! mesh and coordinate transformation.
  a_field_mat(1,1)%base => new_scalar_field_2d_analytic_alt( &
       func_one, &
       "a11", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_PERIODIC, &
       SLL_PERIODIC ) 
  
  a_field_mat(1,2)%base => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "a12", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET,&
       SLL_PERIODIC,&
       SLL_PERIODIC)
  
  a_field_mat(2,1)%base => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "a21", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET,&
       SLL_PERIODIC,&
       SLL_PERIODIC) 
  
  a_field_mat(2,2)%base => new_scalar_field_2d_analytic_alt( &
       func_one, &
       "a22", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET,&
       SLL_PERIODIC,&
       SLL_PERIODIC) 
  
  
  c_field => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "c_field", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_PERIODIC, &
       SLL_PERIODIC )
  
  rho => new_scalar_field_2d_analytic_alt( &
       source_term_chgt_dirper, &
       "rho8", &     
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET,&
       SLL_PERIODIC,&
       SLL_PERIODIC)
  
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
       values, &
       "phi8", &
       interp_2d, &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET,&
       SLL_PERIODIC, &
       SLL_PERIODIC)
  
  print *, 'initialized fields...'

  call set_time_mark(t_reference)
  
  call initialize_general_qn_solver( &
       qns, &
       SPLINE_DEG1, &
       SPLINE_DEG2, &
       NUM_CELLS1, &
       NUM_CELLS2, &
       QNS_GAUSS_LEGENDRE, &
       QNS_GAUSS_LEGENDRE, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET,&
       SLL_PERIODIC,&
       SLL_PERIODIC,&
       ETA1MIN, &
       ETA1MAX, &
       ETA2MIN, &
       ETA2MAX)
  
  t8i = time_elapsed_since(t_reference) 

  print *, 'Initialized QNS object'
  call set_time_mark(t_reference)
  
  ! solve the field
  call solve_quasi_neutral_eq_general_coords( &
       qns, &
       a_field_mat, &
       c_field, &
       rho, &
       phi )
  
  !print *, 'Completed solution',qns%phi_vec
    
  t8e = time_elapsed_since(t_reference)


  
!  print *, 'Compare the values of the transformation at the nodes: '
  
  acc8 = 0.0_f64
  
  do j=0,npts2-1
     do i=0,npts1-1
        eta1       = real(i,f64)*h1 + ETA1MIN
        eta2       = real(j,f64)*h2 + ETA2MIN
        
        node_val   =phi%value_at_point(eta1,eta2)
        !print*, 'rer'
        ref        = sol_exacte_chgt_dirper(eta1,eta2)
        calculated(i+1,j+1) = node_val
        difference(i+1,j+1) = ref-node_val

        if(PRINT_COMPARISON) then
           print *, '(eta1,eta2) = ', eta1, eta2, 'calculated = ', node_val, &
                'theoretical = ', ref
        end if
        acc8        = acc8 + abs(node_val-ref)
     end do
  end do
  call phi%write_to_file(0)
  ! delete things...
  call delete(qns)
  call rho%delete()
  call c_field%delete()
  call phi%delete()
  call a_field_mat(1,1)%base%delete()
  call a_field_mat(1,2)%base%delete()
  call a_field_mat(2,1)%base%delete()
  call a_field_mat(2,2)%base%delete()

  call T%delete()
  
  DEALLOCATE(values)
  DEALLOCATE(calculated)
  DEALLOCATE(difference)



 !--------------------------------------------------------------------
  
  !     95  test case without change of coordinates 
  !      periodic-periodic boundary conditions
  !      and with a non analytic source term
  
  !--------------------------------------------------------------------
  
  
  
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
       identity_jac22)
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
       SLL_PERIODIC,&
       SLL_PERIODIC,&
       SLL_PERIODIC)
  
  a_field_mat(2,1)%base => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "a21", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC,&
       SLL_PERIODIC,&
       SLL_PERIODIC) 
  
  a_field_mat(2,2)%base => new_scalar_field_2d_analytic_alt( &
       func_one, &
       "a22", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC,&
       SLL_PERIODIC,&
       SLL_PERIODIC) 
  
  
  c_field => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "c_field", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC )

  

  allocate(point1(2*npts1-1))
  allocate(point2(2*npts2-1))
  allocate(tab_rho(2*npts1-1,2*npts2-1))
  do j=0,2*npts2-2
     do i=0,2*npts1-2
        point1(i+1)       = real(i,f64)*(ETA1MAX-ETA1MIN)/(2*npts1-1) + ETA1MIN 
        point2(j+1)       = real(j,f64)*(ETA2MAX-ETA2MIN)/(2*npts2-1) + ETA2MIN 
        tab_rho(i+1,j+1)  = source_term_perper(point1(i+1),point2(j+1))
     end do
  end do

  !print*, point1,point2
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
       tab_rho, &
       "rho95", &
       interp_2d_term_source, &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC,&
       SLL_PERIODIC,&
       SLL_PERIODIC,&
       point1,&
       2*npts1-1,&
       point2,&
       2*npts2-1)



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
       values, &
       "phi95", &
       interp_2d, &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC,&
       SLL_PERIODIC, &
       SLL_PERIODIC)
  
  print *, 'initialized fields...'
  
  call set_time_mark(t_reference)
  
  call initialize_general_qn_solver( &
       qns, &
       SPLINE_DEG1, &
       SPLINE_DEG2, &
       NUM_CELLS1, &
       NUM_CELLS2, &
       QNS_GAUSS_LEGENDRE, &
       QNS_GAUSS_LEGENDRE, &
       SLL_PERIODIC, &
       SLL_PERIODIC,&
       SLL_PERIODIC,&
       SLL_PERIODIC,&
       ETA1MIN, &
       ETA1MAX, &
       ETA2MIN, &
       ETA2MAX)
  
  t95i = time_elapsed_since(t_reference) 
  
  print *, 'Initialized QNS object'
  call set_time_mark(t_reference)
  
  ! solve the field
  call solve_quasi_neutral_eq_general_coords( &
       qns, &
       a_field_mat, &
       c_field, &
       rho, &
       phi )
  
  !print *, 'Completed solution',qns%phi_vec
    
  t95e = time_elapsed_since(t_reference)
  
  print *, 'Compare the values of the transformation at the nodes: '
  
  acc95 = 0.0_f64
  
  do j=0,npts2-1
     do i=0,npts1-1
        eta1       = real(i,f64)*h1 + ETA1MIN
        eta2       = real(j,f64)*h2 + ETA2MIN
        
        node_val   =phi%value_at_point(eta1,eta2)
        !print*, 'value at node', node_val
        !print*, 'rer'
        ref        = sol_exacte_perper(eta1,eta2)
        calculated(i+1,j+1) = node_val
        difference(i+1,j+1) = ref-node_val
        if(PRINT_COMPARISON) then
           print *, '(eta1,eta2) = ', eta1, eta2, 'calculated = ', node_val, &
                'theoretical = ', ref
        end if
        acc95        = acc95 + abs(node_val-ref)
     end do
  end do
  call phi%write_to_file(0)
  ! delete things...
  call delete(qns)
  call rho%delete()
  call c_field%delete()
  call phi%delete()
  call a_field_mat(1,1)%base%delete()
  call a_field_mat(1,2)%base%delete()
  call a_field_mat(2,1)%base%delete()
  call a_field_mat(2,2)%base%delete()

  call T%delete()
  
  DEALLOCATE(values)
  DEALLOCATE(calculated)
  DEALLOCATE(difference)
  DEALLOCATE(point1)
  DEALLOCATE(point2)
  DEALLOCATE(tab_rho)




   !--------------------------------------------------------------------
  
  !     9  test case with colella change of coordinates 
  !      periodic-periodic boundary conditions
  !      and with a non analytic source term
  
  !--------------------------------------------------------------------
  
  
  
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
       sinprod_jac22)
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
       SLL_PERIODIC,&
       SLL_PERIODIC,&
       SLL_PERIODIC)
  
  a_field_mat(2,1)%base => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "a21", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC,&
       SLL_PERIODIC,&
       SLL_PERIODIC) 
  
  a_field_mat(2,2)%base => new_scalar_field_2d_analytic_alt( &
       func_one, &
       "a22", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC,&
       SLL_PERIODIC,&
       SLL_PERIODIC) 
  
  
  c_field => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "c_field", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC )

  

  allocate(point1(npts1-1))
  allocate(point2(npts2-1))
  allocate(tab_rho(npts1-1,npts2-1))
  
  do j=0,npts2-2
     do i=0,npts1-2
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
       tab_rho, &
       "rho9", &
       terme_source_interp, &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC,&
       SLL_PERIODIC,&
       SLL_PERIODIC,&
       point1,&
       npts1-1,&
       point2,&
       npts2-1)


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
       values, &
       "phi9", &
       interp_2d, &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC,&
       SLL_PERIODIC, &
       SLL_PERIODIC)
  
  print *, 'initialized fields...'
  
  call set_time_mark(t_reference)
  
  call initialize_general_qn_solver( &
       qns, &
       SPLINE_DEG1, &
       SPLINE_DEG2, &
       NUM_CELLS1, &
       NUM_CELLS2, &
       QNS_GAUSS_LEGENDRE, &
       QNS_GAUSS_LEGENDRE, &
       SLL_PERIODIC, &
       SLL_PERIODIC,&
       SLL_PERIODIC,&
       SLL_PERIODIC,&
       ETA1MIN, &
       ETA1MAX, &
       ETA2MIN, &
       ETA2MAX)
  
  t9i = time_elapsed_since(t_reference) 
  
  print *, 'Initialized QNS object'
  call set_time_mark(t_reference)
  
  ! solve the field
  call solve_quasi_neutral_eq_general_coords( &
       qns, &
       a_field_mat, &
       c_field, &
       rho, &
       phi )
  
  !print *, 'Completed solution',qns%phi_vec
    
  t9e = time_elapsed_since(t_reference)
  
  print *, 'Compare the values of the transformation at the nodes: '
  
  acc9 = 0.0_f64
  
  do j=0,npts2-1
     do i=0,npts1-1
        eta1       = real(i,f64)*h1 + ETA1MIN
        eta2       = real(j,f64)*h2 + ETA2MIN
        
        node_val   =phi%value_at_point(eta1,eta2)
        !print*, 'value at node', node_val
        !print*, 'rer'
        ref        = sol_exacte_chgt_perper(eta1,eta2)
        calculated(i+1,j+1) = node_val
        difference(i+1,j+1) = ref-node_val
        if(PRINT_COMPARISON) then
           print *, '(eta1,eta2) = ', eta1, eta2, 'calculated = ', node_val, &
                'theoretical = ', ref
        end if
        acc9        = acc9 + abs(node_val-ref)
     end do
  end do
  call phi%write_to_file(0)
  call rho%write_to_file(0)
  ! delete things...
  call delete(qns)
  call rho%delete()
  call c_field%delete()
  call phi%delete()
  call a_field_mat(1,1)%base%delete()
  call a_field_mat(1,2)%base%delete()
  call a_field_mat(2,1)%base%delete()
  call a_field_mat(2,2)%base%delete()

  call T%delete()
  
  DEALLOCATE(values)
  DEALLOCATE(calculated)
  DEALLOCATE(difference)
  DEALLOCATE(point1)
  DEALLOCATE(point2)
  DEALLOCATE(tab_rho)

!!$
!!$
!!$   !--------------------------------------------------------------------
!!$  
!!$  !     10  test case with colella change of coordinates 
!!$  !      periodic-dirichlet boundary conditions
!!$  !     and non analytic source term
!!$  !--------------------------------------------------------------------
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
       sinprod_jac22)
  print *, 'initialized coordinate transformation'
  
  ! Thirdly, each field object must be initialized using the same logical
  ! mesh and coordinate transformation.
  a_field_mat(1,1)%base => new_scalar_field_2d_analytic_alt( &
       func_one, &
       "a11", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET) 
  
  a_field_mat(1,2)%base => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "a12", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC,&
       SLL_DIRICHLET,&
       SLL_DIRICHLET)
  
  a_field_mat(2,1)%base => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "a21", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC,&
       SLL_DIRICHLET,&
       SLL_DIRICHLET) 
  
  a_field_mat(2,2)%base => new_scalar_field_2d_analytic_alt( &
       func_one, &
       "a22", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC,&
       SLL_DIRICHLET,&
       SLL_DIRICHLET) 
  
 
  
  c_field => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "c_field", &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET )

  

  allocate(point1(npts1-1))
  allocate(point2(npts2))
  allocate(tab_rho(npts1-1,npts2))
  do j=0,npts2-1
     do i=0,npts1-2
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
       tab_rho, &
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
       values, &
       "phi10", &
       interp_2d, &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC,&
       SLL_DIRICHLET, &
       SLL_DIRICHLET)
  
  print *, 'initialized fields...'
  
  call set_time_mark(t_reference)
  
  call initialize_general_qn_solver( &
       qns, &
       SPLINE_DEG1, &
       SPLINE_DEG2, &
       NUM_CELLS1, &
       NUM_CELLS2, &
       QNS_GAUSS_LEGENDRE, &
       QNS_GAUSS_LEGENDRE, &
       SLL_PERIODIC, &
       SLL_PERIODIC,&
       SLL_DIRICHLET,&
       SLL_DIRICHLET,&
       ETA1MIN, &
       ETA1MAX, &
       ETA2MIN, &
       ETA2MAX)

  t10i = time_elapsed_since(t_reference) 
  
  print *, 'Initialized QNS object'
  call set_time_mark(t_reference)
  
  ! solve the field
  call solve_quasi_neutral_eq_general_coords( &
       qns, &
       a_field_mat, &
       c_field, &
       rho, &
       phi )
 
  
  t10e = time_elapsed_since(t_reference)
  
  print *, 'Compare the values of the transformation at the nodes: '
  
  acc10 = 0.0_f64
  
  do j=0,npts2-1
     do i=0,npts1-1
        eta1       = real(i,f64)*h1 + ETA1MIN
        eta2       = real(j,f64)*h2 + ETA2MIN
        
        node_val   =phi%value_at_point(eta1,eta2)
        !print*, 'rer'
        ref        = sol_exacte_chgt_perdir(eta1,eta2)
        calculated(i+1,j+1) = node_val
        difference(i+1,j+1) = ref-node_val
        if(PRINT_COMPARISON) then
           print *, '(eta1,eta2) = ', eta1, eta2, 'calculated = ', node_val, &
                'theoretical = ', ref
        end if
        acc10        = acc10 + abs(node_val-ref)
     end do
  end do

  call phi%write_to_file(0)
  ! delete things...
  call delete(qns)
  call rho%delete()
  call c_field%delete()
  call phi%delete()
  call a_field_mat(1,1)%base%delete()
  call a_field_mat(1,2)%base%delete()
  call a_field_mat(2,1)%base%delete()
  call a_field_mat(2,2)%base%delete()

  call T%delete()
  
  DEALLOCATE(values)
  DEALLOCATE(calculated)
  DEALLOCATE(difference)
  DEALLOCATE(point1)
  DEALLOCATE(point2)
  DEALLOCATE(tab_rho)
!!$

  !--------------------------------------------------------------------
  
  !     11  test case with colella change of coordinates 
  !      dirichlet-dirichlet boundary conditions
  !     and non analytic source term
  !--------------------------------------------------------------------
  
  
  
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
       sinprod_jac22)
  print *, 'initialized coordinate transformation'
  
  ! Thirdly, each field object must be initialized using the same logical
  ! mesh and coordinate transformation.
  a_field_mat(1,1)%base => new_scalar_field_2d_analytic_alt( &
       func_one, &
       "a11", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET) 
  
  a_field_mat(1,2)%base => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "a12", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET,&
       SLL_DIRICHLET,&
       SLL_DIRICHLET)
  
  a_field_mat(2,1)%base => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "a21", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET,&
       SLL_DIRICHLET,&
       SLL_DIRICHLET) 
  
  a_field_mat(2,2)%base => new_scalar_field_2d_analytic_alt( &
       func_one, &
       "a22", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET,&
       SLL_DIRICHLET,&
       SLL_DIRICHLET) 
  
  
  c_field => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "c_field", &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET )



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
       tab_rho, &
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
       values, &
       "phi11", &
       interp_2d, &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET,&
       SLL_DIRICHLET, &
       SLL_DIRICHLET)
  
  print *, 'initialized fields...'
  call set_time_mark(t_reference)
  
  call initialize_general_qn_solver( &
       qns, &
       SPLINE_DEG1, &
       SPLINE_DEG2, &
       NUM_CELLS1, &
       NUM_CELLS2, &
       QNS_GAUSS_LEGENDRE, &
       QNS_GAUSS_LEGENDRE, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET,&
       SLL_DIRICHLET,&
       SLL_DIRICHLET,&
       ETA1MIN, &
       ETA1MAX, &
       ETA2MIN, &
       ETA2MAX)
  
  t11i = time_elapsed_since(t_reference) 

  print *, 'Initialized QNS object'
  call set_time_mark(t_reference)
  
  ! solve the field
  call solve_quasi_neutral_eq_general_coords( &
       qns, &
       a_field_mat, &
       c_field, &
       rho, &
       phi )
 
  
  t11e = time_elapsed_since(t_reference)
  
  print *, 'Compare the values of the transformation at the nodes: '
  
  acc11 = 0.0_f64
  
  do j=0,npts2-1
     do i=0,npts1-1
        eta1       = real(i,f64)*h1 + ETA1MIN
        eta2       = real(j,f64)*h2 + ETA2MIN
        
        node_val   =phi%value_at_point(eta1,eta2)
        !print*, 'rer'
        ref        = sol_exacte_chgt_dirdir(eta1,eta2)
        calculated(i+1,j+1) = node_val
        difference(i+1,j+1) = ref-node_val
        if(PRINT_COMPARISON) then
           print *, '(eta1,eta2) = ', eta1, eta2, 'calculated = ', node_val, &
                'theoretical = ', ref
        end if
        acc11        = acc11 + abs(node_val-ref)
     end do
  end do
  call phi%write_to_file(0)
  ! delete things...
  call delete(qns)
  call rho%delete()
  call c_field%delete()
  call phi%delete()
  call a_field_mat(1,1)%base%delete()
  call a_field_mat(1,2)%base%delete()
  call a_field_mat(2,1)%base%delete()
  call a_field_mat(2,2)%base%delete()

  call T%delete()
  
  DEALLOCATE(values)
  DEALLOCATE(calculated)
  DEALLOCATE(difference)
  DEALLOCATE(point1)
  DEALLOCATE(point2)
  DEALLOCATE(tab_rho)

  !--------------------------------------------------------------------
  
  !     12  test case with colella change of coordinates 
  !      dirichlet-periodic boundary conditions
  !     and non analytic source term
  !--------------------------------------------------------------------
  
  
  
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
       sinprod_jac22)
  print *, 'initialized coordinate transformation'
  
  ! Thirdly, each field object must be initialized using the same logical
  ! mesh and coordinate transformation.
  a_field_mat(1,1)%base => new_scalar_field_2d_analytic_alt( &
       func_one, &
       "a11", &
       T, &
       SLL_DIRICHLET,&
       SLL_DIRICHLET,&
       SLL_PERIODIC, &
       SLL_PERIODIC) 
  
  a_field_mat(1,2)%base => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "a12", &
       T, &
       SLL_DIRICHLET,&
       SLL_DIRICHLET,&
       SLL_PERIODIC, &
       SLL_PERIODIC)
  
  a_field_mat(2,1)%base => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "a21", &
       T, &
       SLL_DIRICHLET,&
       SLL_DIRICHLET,&
       SLL_PERIODIC, &
       SLL_PERIODIC) 
  
  a_field_mat(2,2)%base => new_scalar_field_2d_analytic_alt( &
       func_one, &
       "a22", &
       T, &
       SLL_DIRICHLET,&
       SLL_DIRICHLET,&
       SLL_PERIODIC, &
       SLL_PERIODIC) 
  
  
  c_field => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "c_field", &
       T, &
       SLL_DIRICHLET,&
       SLL_DIRICHLET,&
       SLL_PERIODIC, &
       SLL_PERIODIC )
  


  allocate(point1(npts1))
  allocate(point2(npts2-1))
  allocate(tab_rho(npts1,npts2-1))
  do j=0,npts2-2
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
       tab_rho, &
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
       values, &
       "phi12", &
       interp_2d, &
       T, &
       SLL_DIRICHLET,&
       SLL_DIRICHLET,&
       SLL_PERIODIC, &
       SLL_PERIODIC)
   
  print *, 'initialized fields...'
  
  call set_time_mark(t_reference)

  call initialize_general_qn_solver( &
       qns, &
       SPLINE_DEG1, &
       SPLINE_DEG2, &
       NUM_CELLS1, &
       NUM_CELLS2, &
       QNS_GAUSS_LEGENDRE, &
       QNS_GAUSS_LEGENDRE, &
       SLL_DIRICHLET,&
       SLL_DIRICHLET,&
       SLL_PERIODIC, &
       SLL_PERIODIC,&
       ETA1MIN, &
       ETA1MAX, &
       ETA2MIN, &
       ETA2MAX)
  
  t12i = time_elapsed_since(t_reference) 

  print *, 'Initialized QNS object'
  call set_time_mark(t_reference)
  
  ! solve the field
  call solve_quasi_neutral_eq_general_coords( &
       qns, &
       a_field_mat, &
       c_field, &
       rho, &
       phi )
  
  
  t12e = time_elapsed_since(t_reference)
  
  print *, 'Compare the values of the transformation at the nodes: '
  
  acc12 = 0.0_f64
  
  do j=0,npts2-1
     do i=0,npts1-1
        eta1       = real(i,f64)*h1 + ETA1MIN
        eta2       = real(j,f64)*h2 + ETA2MIN
        
        node_val   = phi%value_at_point(eta1,eta2)
        !print*, 'rer'
        ref        = sol_exacte_chgt_dirper(eta1,eta2)
        calculated(i+1,j+1) = node_val
        difference(i+1,j+1) = ref-node_val
        if(PRINT_COMPARISON) then
           print *, '(eta1,eta2) = ', eta1, eta2, 'calculated = ', node_val, &
                'theoretical = ', ref
        end if
        acc12        = acc12 + abs(node_val-ref)
     end do
  end do
  call phi%write_to_file(0)
  ! delete things...
  call delete(qns)
  call rho%delete()
  call c_field%delete()
  call phi%delete()
  call a_field_mat(1,1)%base%delete()
  call a_field_mat(1,2)%base%delete()
  call a_field_mat(2,1)%base%delete()
  call a_field_mat(2,2)%base%delete()

  call T%delete()
  
  DEALLOCATE(values)
  DEALLOCATE(calculated)
  DEALLOCATE(difference)
  DEALLOCATE(point1)
  DEALLOCATE(point2)
  DEALLOCATE(tab_rho)
!!$
!!$
  print*, '------------------------------------------------------'
  print*, ' WITHOUT CHANGE OF COORDINATES AND ANALYTIC DATA' 
  print*, '-----------------------------------------------------'

  print *,'Average error in nodes (per-per) without change of coordinates='&
       ,acc1/(npts1*npts2), ',  initialization time (s): ', t1i, &
       ',  solution time (s): ', t1e
  print *,'Average error in nodes (per-dir) without change of coordinates='&
       ,acc2/(npts1*npts2), ',  initialization time (s): ', t2i, &
       ',  solution time (s): ', t2e
  print *,'Average error in nodes (dir-dir) without change of coordinates='&
       ,acc3/(npts1*npts2), ',  initialization time (s): ', t3i, &
       ',  solution time (s): ', t3e
  print *,'Average error in nodes (dir-per) without change of coordinates='&
       ,acc4/(npts1*npts2), ',  initialization time (s): ', t4i, &
       ',  solution time (s): ', t4e

  print*, '-------------------------------------------------------'
  print*, ' COLELLA CHANGE OF COORDINATES AND ANALYTIC DATA' 
  print*, '-------------------------------------------------------'
  print *,'Average error in nodes (per-per) '
  print*, 'with colella change of coordinates='&
       ,acc5/(npts1*npts2), ',  initialization time (s): ', t5i, &
       ',  solution time (s): ', t5e
  print *,'Average error in nodes (per-dir) '
  print*, 'with colella change of coordinates='&
       ,acc6/(npts1*npts2), ',  initialization time (s): ', t6i, &
       ',  solution time (s): ', t6e
  print *,'Average error in nodes (dir-dir) '
  print*, 'with colella change of coordinates='&
       ,acc7/(npts1*npts2), ',  initialization time (s): ', t7i, &
       ',  solution time (s): ', t7e
  print *,'Average error in nodes (dir-per) '
  print*, 'with colella change of coordinates='&
       ,acc8/(npts1*npts2), ',  initialization time (s): ', t8i, &
       ',  solution time (s): ', t8e
  
  print*, '-------------------------------------------------------'
  print*, ' WITHOUT CHANGE OF COORDINATES AND ANALYTIC DATA' 
  print*, '-------------------------------------------------------'
  print *,'Average error in nodes (per-per) '
  print*, 'without change of coordinates='&
       ,acc95/(npts1*npts2), ',  initialization time (s): ', t95i, &
       ',  solution time (s): ', t95e


  print*, '-------------------------------------------------------'
  print*, ' COLELLA CHANGE OF COORDINATES AND WITH A SOURCE TERM NON-ANALYTIC' 
  print*, '-------------------------------------------------------'
  print *,'Average error in nodes (per-per) '
  print*, 'with colella change of coordinates and source term non analytic='&
       ,acc9/(npts1*npts2), ',  initialization time (s): ', t9i, &
       ',  solution time (s): ', t9e
  print *,'Average error in nodes (per-dir) '
  print*, 'with colella change of coordinates and source term non analytic='&
       ,acc10/(npts1*npts2), ',  initialization time (s): ', t10i, &
       ',  solution time (s): ', t10e
  print *,'Average error in nodes (dir-dir)'
  print*, 'with colella change of coordinates and source term non analytic='&
       ,acc11/(npts1*npts2), ',  initialization time (s): ', t11i, &
       ',  solution time (s): ', t11e
  print *,'Average error in nodes (dir-per) '
  print*, 'with colella change of coordinates and source term non analytic='&
       ,acc12/(npts1*npts2), ',  initialization time (s): ', t12i, &
       ',  solution time (s): ', t12e

  
  print *, 'PASSED'
end program test_general_qns




! External functions used as parameters in the above unit test:


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
  res = 2*(2.0*sll_pi)**2*cos(2.0*sll_pi*eta1)*cos(2.0*sll_pi*eta2)
end function source_term_perper

real(8) function sol_exacte_perper(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  !real(8), dimension(:), intent(in), optional :: params
  sol_exacte_perper = cos(2.0*sll_pi*eta1)*cos(2.0*sll_pi*eta2)
end function sol_exacte_perper


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
  source_term_perdir = &
       (16.0*sll_pi**2*eta2**4 &
       - 16.0*sll_pi**2*eta2**2 &
       - 12.0*eta2**2 + 2.0)*cos(2*sll_pi*eta1)*sin(2*sll_pi*eta1)

end function source_term_perdir


real(8) function sol_exacte_perdir(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  intrinsic :: cos
  intrinsic :: sin
  !real(8), dimension(:), intent(in), optional :: params
  sol_exacte_perdir = eta2 ** 2 * (eta2**2-1)&
       * cos(2.0*sll_pi*eta1)*sin(2.0*sll_pi*eta1)
  
  !print*, 'heho'
end function sol_exacte_perdir


!  Solution for a identity change of coordinates 
!   and also dirichlet-periodicconditons
!   the matrix A is equal to identity 
!   the scalar c is equal to zero 
!-------------------------------------------------------------

real(8) function source_term_dirper(eta1,eta2,params) ! in the path
  use sll_constants
  real(8),intent(in) :: eta1,eta2
  real(8), dimension(:), intent(in), optional :: params
  source_term_dirper = &
      (16.0*sll_pi**2*eta1**4 &
      - 16.0*sll_pi**2*eta1**2 &
      - 12.0*eta1**2 + 2.0)*sin(2*sll_pi*eta2)*cos(2*sll_pi*eta2)
end function source_term_dirper


real(8) function sol_exacte_dirper(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  !real(8), dimension(:), intent(in), optional :: params
  sol_exacte_dirper = &
       eta1 ** 2 * (eta1**2-1)* cos(2*sll_pi*eta2)*sin(2*sll_pi*eta2)
  
  
end function sol_exacte_dirper




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
  
  x =   eta1 + 0.1_8*sin(2*sll_pi*eta1) * sin(2*sll_pi*eta2)
  y =   eta2 + 0.1_8*sin(2*sll_pi*eta1) * sin(2*sll_pi*eta2)
  
  source_term_chgt_perper = 8.0*sll_pi**2*cos(2*sll_pi*x)*cos(2*sll_pi*y) 
  
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
  
  source_term_chgt_perdir= 2*(2*sll_pi)**2 * sin(2*sll_pi*y)*cos(2*sll_pi*x)
  
  
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
       2*(2.0*sll_pi)**2*sin(2*sll_pi*x)*sin(2*sll_pi*y)
  
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
  
  
  source_term_chgt_dirper = 2*(2*sll_pi)**2*sin(2*sll_pi*x)*cos(2*sll_pi*y)
  
end function source_term_chgt_dirper




real(8) function sol_exacte_chgt_dirper(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  real(8) :: x,y
  intrinsic :: cos
  intrinsic :: sin

  x =   eta1 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  y =   eta2 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  
  
  sol_exacte_chgt_dirper = cos(2* sll_pi*y)*sin(2* sll_pi*x)
  
end function sol_exacte_chgt_dirper

