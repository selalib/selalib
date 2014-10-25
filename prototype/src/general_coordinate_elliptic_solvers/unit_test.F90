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
  use sll_module_deboor_splines_2d

#ifdef _UMFPACK
  use sll_general_coordinate_elliptic_solver_module_umfpack
#else
  use sll_general_coordinate_elliptic_solver_module
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

  type(sll_logical_mesh_2d), pointer                        :: mesh_2d
  class(sll_coordinate_transformation_2d_base), pointer     :: T
  type(general_coordinate_elliptic_solver)                  :: es
  type(sll_arbitrary_degree_spline_interpolator_2d), target :: interp_2d
  type(sll_arbitrary_degree_spline_interpolator_2d), target :: interp_2d_term_source
  class(sll_interpolator_2d_base), pointer                  :: terme_source_interp
  class(sll_scalar_field_2d_base), pointer                  :: a11_field_mat
  class(sll_scalar_field_2d_base), pointer                  :: a12_field_mat
  class(sll_scalar_field_2d_base), pointer                  :: a21_field_mat
  class(sll_scalar_field_2d_base), pointer                  :: a22_field_mat
  class(sll_scalar_field_2d_base), pointer                  :: b1_field_vect
  class(sll_scalar_field_2d_base), pointer                  :: b2_field_vect
  class(sll_scalar_field_2d_base), pointer                  :: c_field
  class(sll_scalar_field_2d_base), pointer                  :: rho
  type(sll_scalar_field_2d_discrete_alt), pointer           :: phi
  type(sll_time_mark)                                       :: t_reference

  sll_real64 :: ti(13), te(13)

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
  sll_real64 :: acc(13)
  sll_real64 :: normL2(13)
  sll_real64 :: normH1(13)

  sll_real64, dimension(:,:), allocatable    :: calculated
  sll_real64, dimension(:,:), allocatable    :: difference
  sll_real64, dimension(:,:), allocatable    :: reference
  sll_real64, dimension(:,:), allocatable    :: tab_rho

  sll_real64 :: val_jac
  sll_int32  :: ierr
  sll_int32  :: i, j, k
  sll_real64 :: h1,h2,node_val,ref
  sll_real64 :: eta1(NUM_CELLS1+1)
  sll_real64 :: eta2(NUM_CELLS2+1)
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

  sll_real64 :: grad1_node_val,grad2_node_val,grad1ref,grad2ref
  sll_real64, dimension(1) :: whatever  ! dummy params array

  ! First thing, initialize the logical mesh associated with this problem. 
  mesh_2d => new_logical_mesh_2d( NUM_CELLS1, &
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

  k = 1
  print*, "-------------------------------------------------------------"
  print*, "1 test case witout change of coordinates"
  print*, " periodic-periodic boundary conditions "
  print*, "-------------------------------------------------------------"

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

  call initialize_fields( SLL_PERIODIC, SLL_PERIODIC, SLL_PERIODIC, SLL_PERIODIC)

  rho => new_scalar_field_2d_analytic_alt( &
       source_term_perper, &
       "rho1", &     
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       whatever  )
  
  call phi%set_field_data( values )
  call phi%update_interpolation_coefficients( )
  call solve_fields( SLL_PERIODIC, SLL_PERIODIC, SLL_PERIODIC, SLL_PERIODIC, ti(k), te(k))

  normL2(k)= 0.0_f64
  normH1_1 = 0.0_f64

  do j=1,npts2
     do i=1,npts1
        node_val   = phi%value_at_point(eta1(i),eta2(j))
        grad1_node_val = phi%first_deriv_eta1_value_at_point(eta1(i),eta2(j))
        grad2_node_val = phi%first_deriv_eta2_value_at_point(eta1(i),eta2(j))
        ref        = sol_exacte_perper(eta1(i),eta2(j))
        grad1ref   = sol_exacte_perper_der1(eta1(i),eta2(j))
        grad2ref   = sol_exacte_perper_der2(eta1(i),eta2(j))
        calculated(i,j) = node_val
        difference(i,j) = ref-node_val
        reference( i,j) = ref
        acc(k)     = acc(k) + abs(node_val-ref)
        normL2(k)  = normL2(k)+ (node_val-ref)**2*h1*h2
        normH1_1   = normH1_1 + ((grad1_node_val-grad1ref)**2+(grad2_node_val-grad2ref)**2)*h1*h2
     end do
  end do


  call phi%write_to_file(0)
  
  call delete_things()
  
  print*, 'TEST 1'
  if ( ( sqrt(normL2(k) <= h1**(SPLINE_DEG1-1))   .AND. &
       ( sqrt(normH1_1) <= h1**(SPLINE_DEG1-1-1))) then
     print *, 'PASSED'
  else
     stop     'FAILED'
  end if

  print*, "-------------------------------------------------------------"
  print*, " 2 test case witout change of coordinates"
  print*, " periodic-dirichlet boundary conditions"
  print*, "-------------------------------------------------------------"
  k = k+1
  
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

  call initialize_fields( SLL_PERIODIC, SLL_PERIODIC, SLL_DIRICHLET, SLL_DIRICHLET)

  rho => new_scalar_field_2d_analytic_alt( &
       source_term_perdir, &
       "rho2", &     
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       whatever )

  call phi%set_field_data(values)
  call phi%update_interpolation_coefficients()
  call solve_fields( SLL_PERIODIC, SLL_PERIODIC, SLL_DIRICHLET, SLL_DIRICHLET, ti(k), te(k))
  
  normL2(k) = 0.0_f64
  normH1_2 = 0.0_f64
  do j=1,npts2
     do i=1,npts1
        node_val   = phi%value_at_point(eta1(i),eta2(j))
        grad1_node_val = phi%first_deriv_eta1_value_at_point(eta1(i), eta2(j))
        grad2_node_val = phi%first_deriv_eta2_value_at_point(eta1(i), eta2(j))
        ref        = sol_exacte_perdir(eta1(i),eta2(j))
        grad1ref   = sol_exacte_perdir_der1(eta1(i),eta2(j))
        grad2ref   = sol_exacte_perdir_der2(eta1(i),eta2(j))
        !        print*,sin(2*sll_pi*eta1)*cos(2*sll_pi*eta1)
        calculated(i,j) = node_val
        difference(i,j) = ref-node_val
        reference( i,j) = ref
        if(PRINT_COMPARISON) call printout_comparison()
        acc(k)      = acc(k) + abs(node_val-ref)
        normL2(k)    = normL2(k) + (node_val-ref)**2*h1*h2
        normH1_2    = normH1_2 + ((grad1_node_val-grad1ref)**2+(grad2_node_val-grad2ref)**2)*h1*h2
     end do
  end do
  
  call phi%write_to_file(0)
  call delete_things()
  
  if ( ( sqrt(normL2(k)) <= h1**(SPLINE_DEG1-1))   .AND. &
       ( sqrt(normH1_2) <= h1**(SPLINE_DEG1-1-1))) then     
     print *, 'PASSED'
  else
     stop     'FAILED'
  end if

  print*, "-------------------------------------------------------------"
  print*, " 3 test case witout change of coordinates"
  print*, " dirichlet-dirichlet boundary conditions"
  print*, "-------------------------------------------------------------"
  k = k+1
  
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
  
  call initialize_fields( SLL_DIRICHLET, SLL_DIRICHLET, SLL_DIRICHLET, SLL_DIRICHLET)

  rho => new_scalar_field_2d_analytic_alt( &
       source_term_perdir, &
       "rho3", &     
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       whatever )
  
  call phi%set_field_data(values)
  call phi%update_interpolation_coefficients()
  
  call solve_fields( SLL_DIRICHLET, SLL_DIRICHLET, SLL_DIRICHLET, SLL_DIRICHLET, ti(k), te(k))

  normL2(k) = 0.0_f64
  normH1_3 = 0.0_f64
  do j=1,npts2
     do i=1,npts1

        node_val   =phi%value_at_point(eta1(i),eta2(j))
        grad1_node_val = phi%first_deriv_eta1_value_at_point(eta1(i), eta2(j))
        grad2_node_val = phi%first_deriv_eta2_value_at_point(eta1(i), eta2(j))

        ref        = sol_exacte_perdir(eta1(i),eta2(j))
        grad1ref   = sol_exacte_perdir_der1(eta1(i),eta2(j))
        grad2ref   = sol_exacte_perdir_der2(eta1(i),eta2(j))
        calculated(i,j) = node_val
        difference(i,j) = ref-node_val
        reference(i,j) = ref
        if (PRINT_COMPARISON) call printout_comparison()
        acc(k)      = acc(k) + abs(node_val-ref)
        normL2(k)    = normL2(k) + (node_val-ref)**2*h1*h2
        normH1_3    = normH1_3 + ((grad1_node_val-grad1ref)**2+(grad2_node_val-grad2ref)**2)*h1*h2
     end do
  end do
  
  call phi%write_to_file(0)
  call delete_things()
  
  if ( ( sqrt(normL2(k)) <= h1**(SPLINE_DEG1-1))   .AND. &
       ( sqrt(normH1_3) <= h1**(SPLINE_DEG1-1-1))) then     
     print *, 'PASSED'
  else
     stop     'FAILED'
  end if

  print*, "-------------------------------------------------------------"
  print*, " 4 test case witout change of coordinates"
  print*, " dirichlet-periodic boundary conditions"
  print*, "-------------------------------------------------------------"
  k = k+1
  
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
  
  call initialize_fields( SLL_DIRICHLET, SLL_DIRICHLET, &
                          SLL_PERIODIC, SLL_PERIODIC)
  
  rho => new_scalar_field_2d_analytic_alt( &
       source_term_dirper, &
       "rho4", &     
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       whatever )
  
  call phi%set_field_data(values)
  call phi%update_interpolation_coefficients()
  
  call solve_fields( SLL_DIRICHLET, SLL_DIRICHLET, SLL_PERIODIC, SLL_PERIODIC, ti(k), te(k))
  
  normL2(k) = 0.0_f64
  normH1_4 = 0.0_f64
  do j=1,npts2
     do i=1,npts1
        node_val   = phi%value_at_point(eta1(i),eta2(j))
        grad1_node_val = phi%first_deriv_eta1_value_at_point(eta1(i), eta2(j))
        grad2_node_val = phi%first_deriv_eta2_value_at_point(eta1(i), eta2(j))
        ref        = sol_exacte_dirper(eta1(i),eta2(j))
        grad1ref   = sol_exacte_dirper_der1(eta1(i),eta2(j))
        grad2ref   = sol_exacte_dirper_der2(eta1(i),eta2(j))
        calculated(i,j) = node_val
        difference(i,j) = ref-node_val
        reference(i,j) = ref
        if(PRINT_COMPARISON) call printout_comparison()
        acc(k)        = acc(k) + abs(node_val-ref)
        normL2(k)    = normL2(k)+ (node_val-ref)**2*h1*h2
        normH1_4    = normH1_4 + ((grad1_node_val-grad1ref)**2+(grad2_node_val-grad2ref)**2)*h1*h2
     end do
  end do
  call phi%write_to_file(0)

  call delete_things()
  
  if ( ( sqrt(normL2(k) <= h1**(SPLINE_DEG1-1))   .AND. &
       ( sqrt(normH1_4) <= h1**(SPLINE_DEG1-1-1))) then     
     print *, 'PASSED'
  else
     stop     'FAILED'
  end if

  print*, "---------------------"
  print*, " 5 test case with colella change of coordinates"
  print*, " periodic-periodic boundary conditions"
  print*, "---------------------"
  k = k+1
  
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
  
  call initialize_fields( SLL_PERIODIC, SLL_PERIODIC, &
                          SLL_PERIODIC, SLL_PERIODIC)

  rho => new_scalar_field_2d_analytic_alt( &
       source_term_chgt_perper, &
       "rho5", &     
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       whatever )
  
  call phi%set_field_data(values)
  call phi%update_interpolation_coefficients()
  
  call solve_fields( SLL_PERIODIC, SLL_PERIODIC, SLL_PERIODIC, SLL_PERIODIC, ti(k), te(k))
  
  normH1_5 =  0.0
  normL2(k)=  0.0
  do j=1,npts2
     do i=1,npts1
        node_val   = phi%value_at_point(eta1(i),eta2(j))
        grad1_node_val = phi%first_deriv_eta1_value_at_point(eta1(i), eta2(j))
        grad2_node_val = phi%first_deriv_eta2_value_at_point(eta1(i), eta2(j))
        ref        = sol_exacte_chgt_perper(eta1(i),eta2(j))
        grad1ref   = sol_exacte_chgt_perper_der1(eta1(i),eta2(j))
        grad2ref   = sol_exacte_chgt_perper_der2(eta1(i),eta2(j))
        calculated(i,j) = node_val
        difference(i,j) = ref-node_val
        reference(i,j) = ref
        val_jac = sinprod_jac11(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))*&
                  sinprod_jac22(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))-&
                  sinprod_jac12(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))*&
                  sinprod_jac21(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))
        if(PRINT_COMPARISON) call printout_comparison()
        acc(k)        = acc(k) + abs(node_val-ref)
        if ( i < npts1-1 .and. j < npts2-1) then
           integrale_solution = integrale_solution + node_val*val_jac * h1*h2
           integrale_solution_exacte = integrale_solution_exacte + ref*val_jac * h1*h2
           normL2(k)   = normL2(k)+ (node_val-ref)**2*h1*h2*val_jac
           normH1_5    = normH1_5 + ((grad1_node_val-grad1ref)**2+(grad2_node_val-grad2ref)**2)*h1*h2*val_jac
        end if
     end do
  end do
  
  print*, ' integrale solution', integrale_solution
  print*, ' integrale de la solution exacte=', integrale_solution_exacte
  call phi%write_to_file(0)

  call delete_things()
  
  if ( ( sqrt(normL2(k) <= h1**(SPLINE_DEG1-1))   .AND. &
       ( sqrt(normH1_5) <= h1**(SPLINE_DEG1-1-1))) then     
     print *, 'PASSED'
  else
     stop     'FAILED'
  end if

  print*, "-------------------------------------------------------------"
  print*, " 6 test case with colella change of coordinates"
  print*, " periodic-dirichlet boundary conditions"
  print*, "-------------------------------------------------------------"
  k = k+1
  
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
  
  call initialize_fields( SLL_PERIODIC,  SLL_PERIODIC, &
                          SLL_DIRICHLET, SLL_DIRICHLET)
  
  rho => new_scalar_field_2d_analytic_alt( &
       source_term_chgt_perdir, &
       "rho6", &     
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       whatever )
  
  call phi%set_field_data(values)
  call phi%update_interpolation_coefficients()
  
  call sll_set_time_mark(t_reference)
  call solve_fields( SLL_PERIODIC, SLL_PERIODIC, SLL_DIRICHLET, SLL_DIRICHLET, ti(k), te(k))

  integrale_solution = 0.0_f64
  integrale_solution_exacte = 0.0_f64
  normL2(k)= 0.0_f64
  normH1_6 = 0.0_f64
  do j=1,npts2
     do i=1,npts1
        
        node_val   = phi%value_at_point(eta1(i),eta2(j))
        grad1_node_val = phi%first_deriv_eta1_value_at_point(eta1(i), eta2(j))
        grad2_node_val = phi%first_deriv_eta2_value_at_point(eta1(i), eta2(j))
        !print*, 'rer'
        ref        = sol_exacte_chgt_perdir(eta1(i),eta2(j))
        grad1ref   = sol_exacte_chgt_perdir_der1(eta1(i),eta2(j))
        grad2ref   = sol_exacte_chgt_perdir_der2(eta1(i),eta2(j))
        calculated(i,j) = node_val
        difference(i,j) = ref-node_val
        reference(i,j) = ref
        
        if(PRINT_COMPARISON) call printout_comparison()
        acc(k)        = acc(k) + abs(node_val-ref)
        
        val_jac = sinprod_jac11(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))*&
                  sinprod_jac22(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))-&
                  sinprod_jac12(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))*&
                  sinprod_jac21(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))
        if ( i < npts1-1 .and. j < npts2-1 ) then
           integrale_solution = integrale_solution + node_val*val_jac* h1*h2
           integrale_solution_exacte = integrale_solution_exacte + ref*val_jac* h1*h2
           normL2(k)   = normL2(k)+ (node_val-ref)**2*h1*h2
           normH1_6    = normH1_6 + ((grad1_node_val-grad1ref)**2+(grad2_node_val-grad2ref)**2)*h1*h2*val_jac
           
        end if
     end do
  end do
  
  print*, 'integrale solution=',integrale_solution,&
       'integrale de la solution exacte=', integrale_solution_exacte
  call phi%write_to_file(0)

  call delete_things()
  
  print*, 'TEST 6'
  if ( ( sqrt(normL2(k) <= h1**(SPLINE_DEG1-1))   .AND. &
       ( sqrt(normH1_6) <= h1**(SPLINE_DEG1-1-1))) then     
     print *, 'PASSED'
  else
     stop     'FAILED'
  end if

  print*, "-------------------------------------------------------------"
  print*, " 7 test case with colella change of coordinates"
  print*, " dirichlet-dirichlet boundary conditions"
  print*, "-------------------------------------------------------------"
  k = k+1
  
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
  
  call initialize_fields( SLL_DIRICHLET, SLL_DIRICHLET, &
                          SLL_DIRICHLET, SLL_DIRICHLET)

  rho => new_scalar_field_2d_analytic_alt( &
       source_term_chgt_dirdir, &
       "rho7", &     
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET, &
       whatever )
  
  call phi%set_field_data(values)
  call phi%update_interpolation_coefficients()
  
  call solve_fields( SLL_DIRICHLET, SLL_DIRICHLET, SLL_DIRICHLET, SLL_DIRICHLET, ti(k), te(k))
  
  integrale_solution = 0.0_f64
  integrale_solution_exacte = 0.0_f64
  normL2(k)= 0.0_f64
  normH1_7 = 0.0_f64
  do j=1,npts2
     do i=1,npts1
        
        node_val   = phi%value_at_point(eta1(i),eta2(j))
        grad1_node_val = phi%first_deriv_eta1_value_at_point(eta1(i), eta2(j))
        grad2_node_val = phi%first_deriv_eta2_value_at_point(eta1(i), eta2(j))
        ref        = sol_exacte_chgt_dirdir(eta1(i),eta2(j))
        grad1ref   = sol_exacte_chgt_dirdir_der1(eta1(i),eta2(j))
        grad2ref   = sol_exacte_chgt_dirdir_der2(eta1(i),eta2(j))
        calculated(i,j) = node_val
        difference(i,j) = ref-node_val
        reference(i,j) = ref
        if(PRINT_COMPARISON) call printout_comparison()
        acc(k)        = acc(k) + abs(node_val-ref)
        val_jac = sinprod_jac11(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))*&
                  sinprod_jac22(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))-&
                  sinprod_jac12(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))*&
                  sinprod_jac21(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))
        if ( i < npts1-1 .and. j < npts2-1) then
           integrale_solution = integrale_solution + node_val*val_jac * h1*h2
           integrale_solution_exacte = integrale_solution_exacte + ref*val_jac * h1*h2
           normL2(k)   = normL2(k)+ (node_val-ref)**2*h1*h2*val_jac
           normH1_7    = normH1_7 + ((grad1_node_val-grad1ref)**2+(grad2_node_val-grad2ref)**2)*h1*h2*val_jac
           
        end if
     end do
  end do
  print*, 'integrale solution=',integrale_solution,&
       'integrale de la solution excate=', integrale_solution_exacte
  call phi%write_to_file(0)

  call delete_things()

  if ( ( sqrt(normL2(k) <= h1**(SPLINE_DEG1-1))   .AND. &
       ( sqrt(normH1_7) <= h1**(SPLINE_DEG1-1-1))) then     
     print *, 'PASSED'
  else
     stop     'FAILED'
  end if

  print*, "---------------------"
  print*, " 8 test case with colella change of coordinates"
  print*, " dirichlet-periodic boundary conditions"
  print*, "---------------------"
  k = k+1
  
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

  call initialize_fields( SLL_DIRICHLET, SLL_DIRICHLET, &
                          SLL_PERIODIC,  SLL_PERIODIC)
  
  rho => new_scalar_field_2d_analytic_alt( &
       source_term_chgt_dirper, &
       "rho8", &     
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET,&
       SLL_PERIODIC,&
       SLL_PERIODIC, &
       whatever)
  
  call phi%set_field_data(values)
  call phi%update_interpolation_coefficients()
  
  call solve_fields( SLL_DIRICHLET, SLL_DIRICHLET, SLL_PERIODIC, SLL_PERIODIC, ti(k), te(k))
  
  integrale_solution = 0.0_f64
  integrale_solution_exacte = 0.0_f64
  normH1_8 =  0.0
  normL2(k)=  0.0
  do j=1,npts2
     do i=1,npts1
        
        node_val   =phi%value_at_point(eta1(i),eta2(j))
        grad1_node_val = phi%first_deriv_eta1_value_at_point(eta1(i), eta2(j))
        grad2_node_val = phi%first_deriv_eta2_value_at_point(eta1(i), eta2(j))
        ref        = sol_exacte_chgt_dirper(eta1(i),eta2(j))
        grad1ref   = sol_exacte_chgt_dirper_der1(eta1(i),eta2(j))
        grad2ref   = sol_exacte_chgt_dirper_der2(eta1(i),eta2(j))
        calculated(i,j) = node_val
        difference(i,j) = ref-node_val
        reference(i,j) = ref
        if(PRINT_COMPARISON) call printout_comparison()
        acc(k)        = acc(k) + abs(node_val-ref)
        val_jac = sinprod_jac11(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))*&
                  sinprod_jac22(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))-&
                  sinprod_jac12(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))*&
                  sinprod_jac21(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))
        if ( i < npts1-1 .and. j < npts2-1) then
           integrale_solution = integrale_solution + node_val*val_jac * h1*h2
           integrale_solution_exacte = integrale_solution_exacte + ref*val_jac * h1*h2
           normL2(k)   = normL2(k)+ (node_val-ref)**2*h1*h2*val_jac
           normH1_8    = normH1_8 + ((grad1_node_val-grad1ref)**2+(grad2_node_val-grad2ref)**2)*h1*h2*val_jac
        end if
     end do
  end do
  print*, 'integrale solution=',integrale_solution,&
       'integrale de la solution exacte=', integrale_solution_exacte
  call phi%write_to_file(0)

  call delete_things()

  if ( ( sqrt(normL2(k) <= h1**(SPLINE_DEG1-1))   .AND. &
       ( sqrt(normH1_8) <= h1**(SPLINE_DEG1-1-1))) then     
     print *, 'PASSED'
  else
     stop     'FAILED'
  end if

  print*, "---------------------"
  print*, " 95 test case without change of coordinates"
  print*, " periodic-periodic boundary conditions"
  print*, " with non analytic source term " 
  print*, "---------------------"
  k = k+1
  
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
  
  call initialize_fields( SLL_PERIODIC, SLL_PERIODIC, &
                          SLL_PERIODIC, SLL_PERIODIC)

  allocate(tab_rho(npts1,npts2))
  do j=1,npts2
     do i=1,npts1
        tab_rho(i,j)  = source_term_perper(eta1(i),eta2(j))
     end do
  end do
  
  rho => new_scalar_field_2d_discrete_alt( &
       "rho95", &
       interp_2d_term_source, &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC,&
       SLL_PERIODIC,&
       SLL_PERIODIC,&
       eta1,&
       npts1-1,&
       eta2,&
       npts2-1)  
  call rho%set_field_data(tab_rho)
  call rho%update_interpolation_coefficients()

  call rho%write_to_file(0)
 
  call phi%set_field_data(values)
  call phi%update_interpolation_coefficients()
  
  call solve_fields( SLL_PERIODIC, SLL_PERIODIC, SLL_PERIODIC, SLL_PERIODIC, ti(k), te(k))
  
  integrale_solution = 0.0_f64
  integrale_solution_exacte = 0.0_f64
  normL2(k) = 0.0
  normH1_95 = 0.0
  do j=1,npts2
     do i=1,npts1
        
        node_val   =     phi%value_at_point(eta1(i),eta2(j))
        grad1_node_val = phi%first_deriv_eta1_value_at_point(eta1(i), eta2(j))
        grad2_node_val = phi%first_deriv_eta2_value_at_point(eta1(i), eta2(j))
        ref        = sol_exacte_perper(eta1(i),eta2(j))
        grad1ref   = sol_exacte_perper_der1(eta1(i),eta2(j))
        grad2ref   = sol_exacte_perper_der2(eta1(i),eta2(j))
        calculated(i,j) = node_val
        difference(i,j) = ref-node_val
        reference(i,j) = ref
        val_jac = 1.0
        if(PRINT_COMPARISON) call printout_comparison()
        acc(k)        = acc(k) + abs(node_val-ref)
        if ( i < npts1-1 .and. j < npts2-1) then
           integrale_solution = integrale_solution + node_val * h1*h2
           integrale_solution_exacte = integrale_solution_exacte + ref * h1*h2
           normL2(k)    = normL2(k) + (node_val-ref)**2*h1*h2*val_jac
           normH1_95    = normH1_95 + ((grad1_node_val-grad1ref)**2+(grad2_node_val-grad2ref)**2)*h1*h2*val_jac
        end if
     end do
  end do

  print*, 'integrale solution=', integrale_solution,&
       'integrale de la solution exacte=', integrale_solution_exacte
  call phi%write_to_file(0)

  call delete_things()

  if ( ( sqrt(normL2(k))<= h1**(SPLINE_DEG1-1))   .AND. &
       ( sqrt(normH1_95)<= h1**(SPLINE_DEG1-1-1))) then     
     print *, 'PASSED'
  else
     stop     'FAILED'
  end if

  print*, "---------------------"
  print*, " 9 test case with colella change of coordinates"
  print*, " periodic-periodic boundary conditions"
  print*, " with non analytic source term " 
  print*, "---------------------"
  k = k+1
  
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
  
  call initialize_fields( SLL_PERIODIC, SLL_PERIODIC, SLL_PERIODIC, SLL_PERIODIC)
  
  do j=1,npts2
     do i=1,npts1
        tab_rho(i,j)  = source_term_chgt_perper(eta1(i),eta2(j))
     end do
  end do
  
  terme_source_interp => interp_2d_term_source
  
  tab_rho(:,:) = tab_rho - sum(tab_rho)/((npts1-1)*(npts2-1))
  rho => new_scalar_field_2d_discrete_alt( &
       "rho9", &
       terme_source_interp, &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC,&
       SLL_PERIODIC,&
       SLL_PERIODIC ,&
       eta1,&
       npts1-1,&
       eta2,&
       npts2-1)
  call rho%set_field_data(tab_rho)
  call rho%update_interpolation_coefficients()
  
  call phi%set_field_data(values)
  call phi%update_interpolation_coefficients()
  
  call solve_fields( SLL_PERIODIC, SLL_PERIODIC, SLL_PERIODIC, SLL_PERIODIC, ti(k), te(k))
  
  integrale_solution = 0.0_f64
  integrale_solution_exacte = 0.0_f64
  normL2(k)= 0.0
  normH1_9 = 0.0
  do j=1,npts2
     do i=1,npts1
        
        node_val   =phi%value_at_point(eta1(i),eta2(j))
        grad1_node_val = phi%first_deriv_eta1_value_at_point(eta1(i), eta2(j))
        grad2_node_val = phi%first_deriv_eta2_value_at_point(eta1(i), eta2(j))
        ref        = sol_exacte_chgt_perper(eta1(i),eta2(j))
        grad1ref   = sol_exacte_chgt_perper_der1(eta1(i),eta2(j))
        grad2ref   = sol_exacte_chgt_perper_der2(eta1(i),eta2(j))
        calculated(i,j) = node_val
        difference(i,j) = ref-node_val
        reference(i,j) = ref
        if(PRINT_COMPARISON) call printout_comparison()
        acc(k)        = acc(k) + abs(node_val-ref)
        
        val_jac = sinprod_jac11(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))*&
                  sinprod_jac22(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))-&
                  sinprod_jac12(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))*&
                  sinprod_jac21(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))
        
        if ( i < npts1-1 .and. j < npts2-1) then
           integrale_solution = integrale_solution + node_val*val_jac * h1*h2
           integrale_solution_exacte = integrale_solution_exacte + ref*val_jac * h1*h2
           normL2(k)   = normL2(k)+ (node_val-ref)**2*h1*h2*val_jac
           normH1_9    = normH1_9 + ((grad1_node_val-grad1ref)**2+(grad2_node_val-grad2ref)**2)*h1*h2*val_jac
           
        end if
     end do
  end do
  
  print*, 'integrale solution', integrale_solution
  print*, ' integrale de la solution exacte=', integrale_solution_exacte

  call phi%write_to_file(0)
  call rho%write_to_file(0)

  call delete_things()

  if ( ( sqrt(normL2(k) <= h1**(SPLINE_DEG1-1))   .AND. &
       ( sqrt(normH1_9) <= h1**(SPLINE_DEG1-1-1))) then     
     print *, 'PASSED'
  else
     stop     'FAILED'
  end if

  print*, "---------------------"
  print*, " 10 test case with colella change of coordinates"
  print*, " periodic-dirichlet boundary conditions"
  print*, " with non analytic source term " 
  print*, "---------------------"
  k = k+1
  
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

  call initialize_fields( SLL_PERIODIC,  SLL_PERIODIC, &
                          SLL_DIRICHLET, SLL_DIRICHLET)

  do j=1,npts2
     do i=1,npts1
        tab_rho(i,j)  = source_term_chgt_perdir(eta1(i),eta2(j))
     end do
  end do

  terme_source_interp => interp_2d_term_source

  rho => new_scalar_field_2d_discrete_alt( &
       "rho10", &
       terme_source_interp, &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC,&
       SLL_DIRICHLET,&
       SLL_DIRICHLET,&
       eta1,&
       npts1-1,&
       eta2,&
       npts2)

  call rho%set_field_data(tab_rho)
  call rho%update_interpolation_coefficients()

  call phi%set_field_data(values)
  call phi%update_interpolation_coefficients()
  
  call solve_fields( SLL_PERIODIC, SLL_PERIODIC, SLL_DIRICHLET, SLL_DIRICHLET, ti(k), te(k))
  
  integrale_solution = 0.0 
  integrale_solution_exacte = 0.0
  normH1_10 =  0.0
  normL2(k) =  0.0
  do j=1,npts2
     do i=1,npts1
        
        node_val   =phi%value_at_point(eta1(i),eta2(j))
        grad1_node_val = phi%first_deriv_eta1_value_at_point(eta1(i), eta2(j))
        grad2_node_val = phi%first_deriv_eta2_value_at_point(eta1(i), eta2(j))
        ref        = sol_exacte_chgt_perdir(eta1(i),eta2(j))
        grad1ref   = sol_exacte_chgt_perdir_der1(eta1(i),eta2(j))
        grad2ref   = sol_exacte_chgt_perdir_der2(eta1(i),eta2(j))
        calculated(i,j) = node_val
        difference(i,j) = ref-node_val
        reference(i,j) = ref
        if(PRINT_COMPARISON) call printout_comparison()
        acc(k)        = acc(k) + abs(node_val-ref)
        val_jac = sinprod_jac11(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))*&
                  sinprod_jac22(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))-&
                  sinprod_jac12(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))*&
                  sinprod_jac21(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))
        
        if ( i < npts1-1 .and. j < npts2-1) then
           integrale_solution = integrale_solution + node_val*val_jac * h1*h2
           integrale_solution_exacte = integrale_solution_exacte + &
                ref*val_jac * h1*h2
           normL2(k)    = normL2(k) + (node_val-ref)**2*h1*h2*val_jac
           normH1_10   = normH1_10 + ((grad1_node_val-grad1ref)**2+&
                (grad2_node_val-grad2ref)**2)*h1*h2*val_jac
           
        end if
     end do
  end do
  print*, 'integrale de la solution=',integrale_solution
  print*, 'integrale de la solution exacte=',integrale_solution_exacte

  call phi%write_to_file(0)
  call delete_things()
  
  if ( ( sqrt(normL2(k))<= h1**(SPLINE_DEG1-1))   .AND. &
       ( sqrt(normH1_10)<= h1**(SPLINE_DEG1-1-1))) then     
     print *, 'PASSED'
  else
     stop     'FAILED'
  end if

  print*, "---------------------"
  print*, " 11 test case with colella change of coordinates"
  print*, " dirichlet-dirichlet boundary conditions"
  print*, " with non analytic source term " 
  print*, "---------------------"
  k = k+1
  
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

  call initialize_fields( SLL_DIRICHLET, SLL_DIRICHLET, &
                          SLL_DIRICHLET, SLL_DIRICHLET)

  do j=1,npts2
     do i=1,npts1
        tab_rho(i,j)  = source_term_chgt_dirdir(eta1(i),eta2(j))
     end do
  end do

  terme_source_interp => interp_2d_term_source

  rho => new_scalar_field_2d_discrete_alt( &
       "rho11", &
       terme_source_interp, &
       T, &
       SLL_DIRICHLET, &
       SLL_DIRICHLET,&
       SLL_DIRICHLET,&
       SLL_DIRICHLET,&
       eta1,&
       npts1,&
       eta2,&
       npts2)

  call rho%set_field_data(tab_rho)
  call rho%update_interpolation_coefficients()
  call phi%set_field_data(values)
  call phi%update_interpolation_coefficients()
  
  call solve_fields( SLL_DIRICHLET, SLL_DIRICHLET, SLL_DIRICHLET, SLL_DIRICHLET, ti(k), te(k))
  
  integrale_solution_exacte = 0.0
  integrale_solution = 0.0
  normL2(k) = 0.0_f64
  normH1_11 = 0.0_f64
  do j=1,npts2
     do i=1,npts1
        
        node_val   =phi%value_at_point(eta1(i),eta2(j))
        grad1_node_val = phi%first_deriv_eta1_value_at_point(eta1(i), eta2(j))
        grad2_node_val = phi%first_deriv_eta2_value_at_point(eta1(i), eta2(j))
        ref        = sol_exacte_chgt_dirdir(eta1(i),eta2(j))
        grad1ref   = sol_exacte_chgt_dirdir_der1(eta1(i),eta2(j))
        grad2ref   = sol_exacte_chgt_dirdir_der2(eta1(i),eta2(j))
        calculated(i,j) = node_val
        difference(i,j) = ref-node_val
        reference(i,j)  = ref

        if(PRINT_COMPARISON) call printout_comparison()
        acc(k)        = acc(k) + abs(node_val-ref)
        val_jac = sinprod_jac11(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))*&
                  sinprod_jac22(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))-&
                  sinprod_jac12(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))*&
                  sinprod_jac21(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))
        
        if ( i < npts1-1 .and. j < npts2-1) then
           integrale_solution = integrale_solution + node_val*val_jac * h1*h2
           integrale_solution_exacte = integrale_solution_exacte + &
                ref*val_jac * h1*h2
           normL2(k)    = normL2(k) + (node_val-ref)**2*h1*h2*val_jac
           normH1_11    = normH1_11 + ((grad1_node_val-grad1ref)**2+&
                (grad2_node_val-grad2ref)**2)*h1*h2*val_jac
        end if
     end do
  end do

  print*, 'integrale de la solution=', integrale_solution
  print*, 'integrale de la solution exacte=', integrale_solution_exacte
  
  call phi%write_to_file(0)
  call delete_things()

  if ( ( sqrt(normL2(k))<= h1**(SPLINE_DEG1-1))   .AND. &
       ( sqrt(normH1_11)<= h1**(SPLINE_DEG1-1-1))) then     
     print *, 'PASSED'
  else
     stop     'FAILED'
  end if

  print*, "---------------------"
  print*, " 12 test case with colella change of coordinates"
  print*, " dirichlet-periodic  boundary conditions"
  print*, " with non analytic source term " 
  print*, "---------------------"
  k = k+1
  
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
  
  call initialize_fields( SLL_DIRICHLET, SLL_DIRICHLET, &
                          SLL_PERIODIC,  SLL_PERIODIC)

  do j=1,npts2
     do i=1,npts1
        tab_rho(i,j)  = source_term_chgt_dirper(eta1(i),eta2(j))
     end do
  end do
  
  terme_source_interp => interp_2d_term_source

  rho => new_scalar_field_2d_discrete_alt( &
       "rho12", &
       terme_source_interp, &
       T, &
       SLL_DIRICHLET,&
       SLL_DIRICHLET,&
       SLL_PERIODIC, &
       SLL_PERIODIC,&
       eta1,&
       npts1,&
       eta2,&
       npts2-1)

  call rho%set_field_data(tab_rho)
  call rho%update_interpolation_coefficients()
  call phi%set_field_data(values)
  call phi%update_interpolation_coefficients()

  call solve_fields( SLL_DIRICHLET, SLL_DIRICHLET, SLL_PERIODIC, SLL_PERIODIC, ti(k), te(k))
   
  integrale_solution = 0.0
  integrale_solution_exacte = 0.0
  normL2(k) = 0.0_f64
  normH1_12 = 0.0_f64
  
  do j=1,npts2
     do i=1,npts1
        
        node_val   = phi%value_at_point(eta1(i),eta2(j))
        grad1_node_val = phi%first_deriv_eta1_value_at_point(eta1(i), eta2(j))
        grad2_node_val = phi%first_deriv_eta2_value_at_point(eta1(i), eta2(j))
        ref        = sol_exacte_chgt_dirper(eta1(i),eta2(j))
        grad1ref   = sol_exacte_chgt_dirper_der1(eta1(i),eta2(j))
        grad2ref   = sol_exacte_chgt_dirper_der2(eta1(i),eta2(j))
        calculated(i,j) = node_val
        difference(i,j) = ref-node_val
        reference(i,j)  = ref
        if (PRINT_COMPARISON) call printout_comparison()
        acc(k)        = acc(k) + abs(node_val-ref)
        val_jac = sinprod_jac11(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))*&
                  sinprod_jac22(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))-&
                  sinprod_jac12(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))*&
                  sinprod_jac21(eta1(i),eta2(j),(/0.1_f64,0.1_f64,1.0_f64,1.0_f64/))
        
        if ( i < npts1-1 .and. j < npts2-1) then
           integrale_solution = integrale_solution + node_val*val_jac * h1*h2
           integrale_solution_exacte = integrale_solution_exacte + &
                ref*val_jac * h1*h2
           normL2(k)    = normL2(k) + (node_val-ref)**2*h1*h2*val_jac
           normH1_12    = normH1_12 + ((grad1_node_val-grad1ref)**2+&
                (grad2_node_val-grad2ref)**2)*h1*h2*val_jac
           
        end if
     end do
  end do
  print*, 'integrale de la solution =', integrale_solution
  print*, 'integrale de la solution exacte=', integrale_solution_exacte
  call phi%write_to_file(0)
  call delete_things()

  if ( ( sqrt(normL2(k))<= h1**(SPLINE_DEG1-1))   .AND. &
       ( sqrt(normH1_12)<= h1**(SPLINE_DEG1-1-1))) then     
     print *, 'PASSED'
  else
     stop     'FAILED'
  end if

  print*, '-----------------------------------------------------'
  print*, ' WITHOUT CHANGE OF COORDINATES AND ANALYTIC DATA' 
  print*, '-----------------------------------------------------'
  acc = acc/(npts1*npts2)
  print*,'Average error in nodes (per-per) without change of coordinates=',acc(1)
  print*,'Norm L2',sqrt(normL2(k),'Norm H1',sqrt(normH1_1)
  print*,'Average error in nodes (per-dir) without change of coordinates=',acc(2)
  print*,'Norm L2',sqrt(normL2(k),'Norm H1',sqrt(normH1_2)
  print*,'Average error in nodes (dir-dir) without change of coordinates=',acc(3)
  print*,'Norm L2',sqrt(normL2(k),'Norm H1',sqrt(normH1_3)
  print*,'Average error in nodes (dir-per) without change of coordinates=',acc(4)
  print*,'Norm L2',sqrt(normL2(k),'Norm H1',sqrt(normH1_4)
  print*,'-------------------------------------------------------'
  print*,' COLELLA CHANGE OF COORDINATES AND ANALYTIC DATA' 
  print*,'-------------------------------------------------------'
  print*,'Average error in nodes (per-per) with colella change of coordinates=',acc(5)
  print*,'Norm L2',sqrt(normL2(k),'Norm H1',sqrt(normH1_5)
  print*,'Average error in nodes (per-dir) with colella change of coordinates=',acc(6)
  print*,'Norm L2',sqrt(normL2(k),'Norm H1',sqrt(normH1_6)
  print*,'Average error in nodes (dir-dir) with colella change of coordinates=',acc(7)
  print*,'Norm L2',sqrt(normL2(k),'Norm H1',sqrt(normH1_7)
  print*,'Average error in nodes (dir-per) with colella change of coordinates=',acc(8)
  print*,'Norm L2',sqrt(normL2(k),'Norm H1',sqrt(normH1_8)
  print*,'------------------------------------------------------------------'
  print*,' WITHOUT CHANGE OF COORDINATES AND WITH A SOURCE TERM NON-ANALYTIC' 
  print*,'------------------------------------------------------------------'
  print*,'Average error in nodes (per-per) without change of coordinates=',acc(9)
  print*,'Norm L2',sqrt(normL2(k)),'Norm H1',sqrt(normH1_95)
  print*,'------------------------------------------------------------------'
  print*,' COLELLA CHANGE OF COORDINATES AND WITH A SOURCE TERM NON-ANALYTIC' 
  print*,'------------------------------------------------------------------'
  print*,'Average error in nodes (per-per) with colella and source term non analytic=',acc(10)
  print*,'Norm L2',sqrt(normL2(k),'Norm H1',sqrt(normH1_9)
  print*,'Average error in nodes (per-dir) with colella and source term non analytic=',acc(11)
  print*,'Norm L2',sqrt(normL2(k)),'Norm H1',sqrt(normH1_10)
  print*,'Average error in nodes (dir-dir) with colella and source term non analytic=',acc(12)
  print*,'Norm L2',sqrt(normL2(k)),'Norm H1',sqrt(normH1_11)
  print*,'Average error in nodes (dir-per) with colella and source term non analytic=',acc(13)
  print*,'Norm L2',sqrt(normL2(k)),'Norm H1',sqrt(normH1_12)
  
contains

  subroutine printout_comparison()
    print *, '(eta1,eta2) = ', eta1, eta2(j), 'calculated = ', node_val, &
             'theoretical = ', ref, 'difference=', node_val-ref
    print *, '(eta1,eta2) = ', eta1, eta2(j), 'calculated = ', grad1_node_val, &
             'theoretical = ', grad1ref, 'difference=',grad1ref-grad1_node_val
    print *, '(eta1,eta2) = ', eta1, eta2(j), 'calculated = ', grad2_node_val, &
             'theoretical = ', grad2ref, 'difference=',grad2ref-grad2_node_val
  end subroutine printout_comparison

  ! Each field object must be initialized using the same logical
  ! mesh and coordinate transformation.
  subroutine initialize_fields( bc_eta1_min, bc_eta1_max, bc_eta2_min, bc_eta2_max)
  sll_int32, intent(in) :: bc_eta1_min
  sll_int32, intent(in) :: bc_eta2_min
  sll_int32, intent(in) :: bc_eta1_max
  sll_int32, intent(in) :: bc_eta2_max

  SLL_ALLOCATE(values(npts1,npts2),ierr)
  SLL_ALLOCATE(calculated(npts1,npts2),ierr)
  SLL_ALLOCATE(difference(npts1,npts2),ierr)
  SLL_ALLOCATE(reference(npts1,npts2),ierr)
  values(:,:) = 0.0_f64

  a11_field_mat => new_scalar_field_2d_analytic_alt( &
       func_one, &
       "a11", &
       T, &
       bc_eta1_min, &
       bc_eta1_max, &
       bc_eta2_min, &
       bc_eta2_max, &
       whatever  ) 
  
  a12_field_mat => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "a12", &
       T, &
       bc_eta1_min, &
       bc_eta1_max,&
       bc_eta2_min,&
       bc_eta2_max, &
       whatever )
  
  a21_field_mat => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "a21", &
       T, &
       bc_eta1_min, &
       bc_eta1_max,&
       bc_eta2_min,&
       bc_eta2_max, &
       whatever ) 
  
  a22_field_mat => new_scalar_field_2d_analytic_alt( &
       func_one, &
       "a22", &
       T, &
       bc_eta1_min, &
       bc_eta1_max,&
       bc_eta2_min,&
       bc_eta2_max, &
       whatever)

  b1_field_vect => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "b1", &
       T, &
       bc_eta1_min, &
       bc_eta1_max, &
       bc_eta2_min, &
       bc_eta2_max, &
       whatever, & 
       first_deriv_eta1 = func_zero, &
       first_deriv_eta2 = func_zero) 

  b2_field_vect => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "b2", &
       T, &
       bc_eta1_min, &
       bc_eta1_max, &
       bc_eta2_min, &
       bc_eta2_max, &
       whatever, &
       first_deriv_eta1 = func_zero, &
       first_deriv_eta2 = func_zero)


  c_field => new_scalar_field_2d_analytic_alt( &
       func_zero, &
       "c_field", &
       T, &
       bc_eta1_min, &
       bc_eta1_max, &
       bc_eta2_min, &
       bc_eta2_max, &
       whatever  )

  call initialize_ad2d_interpolator( &
       interp_2d, &
       NUM_CELLS1+1, &
       NUM_CELLS2+1, &
       ETA1MIN, &
       ETA1MAX, &
       ETA2MIN, &
       ETA2MAX, &
       bc_eta1_min, &
       bc_eta1_max, &
       bc_eta2_min, &
       bc_eta2_max, &
       SPLINE_DEG1, &
       SPLINE_DEG2 )

  call initialize_ad2d_interpolator( &
       interp_2d_term_source, &
       NUM_CELLS1+1, &
       NUM_CELLS2+1, &
       ETA1MIN, &
       ETA1MAX, &
       ETA2MIN, &
       ETA2MAX, &
       bc_eta1_min, &
       bc_eta1_max,&
       bc_eta2_min,&
       bc_eta2_max,&
       SPLINE_DEG1, &
       SPLINE_DEG2 )

  phi => new_scalar_field_2d_discrete_alt( &
       "phi", &
       interp_2d, &
       T, &
       bc_eta1_min, &
       bc_eta1_max, &
       bc_eta2_min, &
       bc_eta2_max )

  end subroutine initialize_fields

  subroutine delete_things()
    print"('integrale de la solution =',g15.3)", sum(calculated(1:npts1-1,1:npts2-1))*h1*h2
    print"('integrale de la solution exacte =',g15.3)", sum(reference(1:npts1-1,1:npts2-1))*h1*h2
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
 end subroutine delete_things

  subroutine solve_fields( bc_eta1_min, bc_eta1_max, bc_eta2_min, bc_eta2_max, &
                           ti, te)

  sll_int32,  intent(in)  :: bc_eta1_min
  sll_int32,  intent(in)  :: bc_eta2_min
  sll_int32,  intent(in)  :: bc_eta1_max
  sll_int32,  intent(in)  :: bc_eta2_max
  sll_real64, intent(out) :: ti
  sll_real64, intent(out) :: te
  sll_int32               :: istep

  integrale_solution = 0.0_f64
  integrale_solution_exacte = 0.0_f64

  call sll_set_time_mark(t_reference)

  call sll_create( &
       es, &
       SPLINE_DEG1, &
       SPLINE_DEG2, &
       NUM_CELLS1, &
       NUM_CELLS2, &
       ES_GAUSS_LEGENDRE, &
       ES_GAUSS_LEGENDRE, &
       bc_eta1_min, &
       bc_eta1_max, &
       bc_eta2_min, &
       bc_eta2_max, &
       ETA1MIN, &
       ETA1MAX, &
       ETA2MIN, &
       ETA2MAX)
 
  ti = sll_time_elapsed_since(t_reference)
 
  call sll_set_time_mark(t_reference)

  call factorize_mat_es(&
       es, &
       a11_field_mat, &
       a12_field_mat,&
       a21_field_mat,&
       a22_field_mat,&
       b1_field_vect,&
       b2_field_vect,&
       c_field)

  do istep = 1, 10

    call sll_solve( es, rho, phi)

  end do

  te = sll_time_elapsed_since(t_reference)

  end subroutine solve_fields

end program test_general_elliptic_solver
