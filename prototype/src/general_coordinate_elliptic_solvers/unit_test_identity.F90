program test_general_elliptic_solver
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_logical_meshes
  use sll_module_coordinate_transformations_2d
  use sll_common_coordinate_transformations
  use sll_module_scalar_field_2d_alternative
  use sll_constants
  use sll_module_arbitrary_degree_spline_interpolator_2d
  use sll_timer


#ifdef _UMFPACK
  use sll_general_coordinate_elliptic_solver_module_umfpack
#else
  use sll_general_coordinate_elliptic_solver_module
#endif

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

  type(sll_logical_mesh_2d), pointer                    :: mesh_2d
  class(sll_coordinate_transformation_2d_base), pointer :: T
  type(general_coordinate_elliptic_solver)              :: es
  type(sll_arbitrary_degree_spline_interpolator_2d), target                 :: interp_2d
  type(sll_arbitrary_degree_spline_interpolator_2d), target                 :: interp_2d_term_source
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

  !*******************************************************************
  !        WHITHOUT CHANGE OF COORDINATES AND ANALYTIC DATA
  !*******************************************************************
  !--------------------------------------------------------------------
  !     1 test case without change of coordinates 
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

  call sll_create( &
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
  call sll_solve(&
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

  print"('integrale de la solution =',g15.3)", &
     sum(calculated(1:npts1-1,1:npts2-1))*h1*h2
  print"('integrale de la solution exacte =',g15.3)", &
     sum(reference(1:npts1-1,1:npts2-1))*h1*h2

  call phi%write_to_file(0)
  
  ! delete things...
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
  
  SLL_DEALLOCATE(values, ierr)
  SLL_DEALLOCATE_ARRAY(calculated,ierr)
  SLL_DEALLOCATE_ARRAY(difference,ierr)
  SLL_DEALLOCATE_ARRAY(reference,ierr)

  print*, '------------------------------------------------------'
  print*, ' WITHOUT CHANGE OF COORDINATES AND ANALYTIC DATA' 
  print*, '-----------------------------------------------------'

  print *,'Average error in nodes (per-per) without change of coordinates='&
       ,acc1/(npts1*npts2), ',  initialization time (s): ', t1i, &
       ',  solution time (s): ', t1e,'Norm L2',sqrt(normL2_1),'Norm H1',sqrt(normH1_1)
  
!!$  
!print*,h1**(SPLINE_DEG1-2)*h2**(SPLINE_DEG2-2)
!borne_L2 = 1.8*sll_pi**2*h1**(SPLINE_DEG1-1)*h2**(SPLINE_DEG2-1)
  if ((sqrt(normL2_1) <= h1**(SPLINE_DEG1-1))  .AND. &
      (sqrt(normH1_1) <= h1**(SPLINE_DEG1-1-1)) ) then
     print *, 'PASSED'
  else
     print *, 'FAILED'
  end if

end program test_general_elliptic_solver




