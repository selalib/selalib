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
  class(sll_interpolator_2d_base), pointer                  :: rhs_interp
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

  sll_real64 :: acc(13)
  sll_real64 :: normL2(13)
  sll_real64 :: normH1(13)

  sll_real64, dimension(NUM_CELLS1+1,NUM_CELLS2+1) :: values
  sll_real64, dimension(NUM_CELLS1+1,NUM_CELLS2+1) :: calculated
  sll_real64, dimension(NUM_CELLS1+1,NUM_CELLS2+1) :: difference
  sll_real64, dimension(NUM_CELLS1+1,NUM_CELLS2+1) :: reference
  sll_real64, dimension(NUM_CELLS1+1,NUM_CELLS2+1) :: tab_rho

  sll_real64 :: val_jac
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

  real(8) :: integral_solution
  real(8) :: integral_exact_solution

  CHARACTER(len=10) :: cmd

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

  normL2 = 0.0_f64
  normH1 = 0.0_f64

  call Get_command_argument(1,cmd)
  read (cmd,'(I2)') k
  print *, i

  do k = 1, 13
  select case(k)
  case(1)
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
  
  call solve_fields( SLL_PERIODIC, SLL_PERIODIC, SLL_PERIODIC, SLL_PERIODIC, ti(k), te(k))

  do j=1,npts2
     do i=1,npts1
        node_val = calculated(i,j)
        grad1_node_val = phi%first_deriv_eta1_value_at_point(eta1(i),eta2(j))
        grad2_node_val = phi%first_deriv_eta2_value_at_point(eta1(i),eta2(j))
        ref        = sol_exacte_perper(eta1(i),eta2(j))
        grad1ref   = sol_exacte_perper_der1(eta1(i),eta2(j))
        grad2ref   = sol_exacte_perper_der2(eta1(i),eta2(j))
        difference(i,j) = ref-node_val
        reference( i,j) = ref
        normL2(k)   = normL2(k) + (node_val-ref)**2*h1*h2
        normH1(k)   = normH1(k) + ((grad1_node_val-grad1ref)**2+(grad2_node_val-grad2ref)**2)*h1*h2
     end do
  end do

  integral_solution = sum(calculated(1:NUM_CELLS1,1:NUM_CELLS2))*h1*h2
  integral_exact_solution = sum(reference(1:NUM_CELLS1,1:NUM_CELLS2))*h1*h2

  call delete_things()
  
  call check_error()

  case(2)
  print*, "-------------------------------------------------------------"
  print*, " 2 test case witout change of coordinates"
  print*, " periodic-dirichlet boundary conditions"
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

  call solve_fields( SLL_PERIODIC, SLL_PERIODIC, SLL_DIRICHLET, SLL_DIRICHLET, ti(k), te(k))
  
  do j=1,npts2
     do i=1,npts1
        node_val = calculated(i,j)
        grad1_node_val = phi%first_deriv_eta1_value_at_point(eta1(i), eta2(j))
        grad2_node_val = phi%first_deriv_eta2_value_at_point(eta1(i), eta2(j))
        ref        = sol_exacte_perdir(eta1(i),eta2(j))
        grad1ref   = sol_exacte_perdir_der1(eta1(i),eta2(j))
        grad2ref   = sol_exacte_perdir_der2(eta1(i),eta2(j))
        difference(i,j) = ref-node_val
        reference( i,j) = ref
        if(PRINT_COMPARISON) call printout_comparison()
        normL2(k)    = normL2(k) + (node_val-ref)**2*h1*h2
        normH1(k)    = normH1(k) + ((grad1_node_val-grad1ref)**2+(grad2_node_val-grad2ref)**2)*h1*h2
     end do
  end do
  
  integral_solution = sum(calculated(1:NUM_CELLS1,1:NUM_CELLS2))*h1*h2
  integral_exact_solution = sum(reference(1:NUM_CELLS1,1:NUM_CELLS2))*h1*h2

  call delete_things()
  
  call check_error()

  case(3)
  print*, "-------------------------------------------------------------"
  print*, " 3 test case witout change of coordinates"
  print*, " dirichlet-dirichlet boundary conditions"
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
  
  call solve_fields( SLL_DIRICHLET, SLL_DIRICHLET, SLL_DIRICHLET, SLL_DIRICHLET, ti(k), te(k))

  do j=1,npts2
     do i=1,npts1
        node_val = calculated(i,j)
        grad1_node_val  = phi%first_deriv_eta1_value_at_point(eta1(i), eta2(j))
        grad2_node_val  = phi%first_deriv_eta2_value_at_point(eta1(i), eta2(j))
        ref             = sol_exacte_perdir(eta1(i),eta2(j))
        grad1ref        = sol_exacte_perdir_der1(eta1(i),eta2(j))
        grad2ref        = sol_exacte_perdir_der2(eta1(i),eta2(j))
        difference(i,j) = ref-node_val
        reference(i,j)  = ref
        normL2(k)  = normL2(k) + (node_val-ref)**2*h1*h2
        normH1(k)  = normH1(k) + ((grad1_node_val-grad1ref)**2+(grad2_node_val-grad2ref)**2)*h1*h2
        if (PRINT_COMPARISON) call printout_comparison()
     end do
  end do
  
  integral_solution = sum(calculated(1:NUM_CELLS1,1:NUM_CELLS2))*h1*h2
  integral_exact_solution = sum(reference(1:NUM_CELLS1,1:NUM_CELLS2))*h1*h2

  call delete_things()
  call check_error()

  case(4)
  print*, "-------------------------------------------------------------"
  print*, " 4 test case witout change of coordinates"
  print*, " dirichlet-periodic boundary conditions"
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
  
  call solve_fields( SLL_DIRICHLET, SLL_DIRICHLET, SLL_PERIODIC, SLL_PERIODIC, ti(k), te(k))
  
  do j=1,npts2
     do i=1,npts1
        node_val = calculated(i,j)
        node_val        = phi%value_at_point(eta1(i),eta2(j))
        grad1_node_val  = phi%first_deriv_eta1_value_at_point(eta1(i), eta2(j))
        grad2_node_val  = phi%first_deriv_eta2_value_at_point(eta1(i), eta2(j))
        ref             = sol_exacte_dirper(eta1(i),eta2(j))
        grad1ref        = sol_exacte_dirper_der1(eta1(i),eta2(j))
        grad2ref        = sol_exacte_dirper_der2(eta1(i),eta2(j))
        calculated(i,j) = node_val
        difference(i,j) = ref-node_val
        reference(i,j)  = ref
        normL2(k)       = normL2(k) + (node_val-ref)**2*h1*h2
        normH1(k)       = normH1(k) + ((grad1_node_val-grad1ref)**2+(grad2_node_val-grad2ref)**2)*h1*h2
        if(PRINT_COMPARISON) call printout_comparison()
     end do
  end do

  integral_solution = sum(calculated(1:NUM_CELLS1,1:NUM_CELLS2))*h1*h2
  integral_exact_solution = sum(reference(1:NUM_CELLS1,1:NUM_CELLS2))*h1*h2

  call delete_things()
  
  call check_error()

  case(5)
  print*, "---------------------"
  print*, " 5 test case with colella change of coordinates"
  print*, " periodic-periodic boundary conditions"
  print*, "---------------------"
  
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
  
  call solve_fields( SLL_PERIODIC, SLL_PERIODIC, SLL_PERIODIC, SLL_PERIODIC, ti(k), te(k))
  
  do j=1,npts2
     do i=1,npts1
        node_val = calculated(i,j)
        grad1_node_val = phi%first_deriv_eta1_value_at_point(eta1(i), eta2(j))
        grad2_node_val = phi%first_deriv_eta2_value_at_point(eta1(i), eta2(j))
        ref        = sol_exacte_chgt_perper(eta1(i),eta2(j))
        grad1ref   = sol_exacte_chgt_perper_der1(eta1(i),eta2(j))
        grad2ref   = sol_exacte_chgt_perper_der2(eta1(i),eta2(j))
        difference(i,j) = ref-node_val
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
           normH1(k)    = normH1(k) + ((grad1_node_val-grad1ref)**2+(grad2_node_val-grad2ref)**2)*h1*h2*val_jac
        end if
     end do
  end do
  
  integral_solution = sum(calculated(1:NUM_CELLS1,1:NUM_CELLS2))*h1*h2
  integral_exact_solution = sum(reference(1:NUM_CELLS1,1:NUM_CELLS2))*h1*h2

  call delete_things()
  
  call check_error()

  case(6)
  print*, "-------------------------------------------------------------"
  print*, " 6 test case with colella change of coordinates"
  print*, " periodic-dirichlet boundary conditions"
  print*, "-------------------------------------------------------------"
  
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
  
  
  call sll_set_time_mark(t_reference)
  call solve_fields( SLL_PERIODIC, SLL_PERIODIC, SLL_DIRICHLET, SLL_DIRICHLET, ti(k), te(k))

  do j=1,npts2
     do i=1,npts1
        
        node_val = calculated(i,j)
        grad1_node_val = phi%first_deriv_eta1_value_at_point(eta1(i), eta2(j))
        grad2_node_val = phi%first_deriv_eta2_value_at_point(eta1(i), eta2(j))
        !print*, 'rer'
        ref        = sol_exacte_chgt_perdir(eta1(i),eta2(j))
        grad1ref   = sol_exacte_chgt_perdir_der1(eta1(i),eta2(j))
        grad2ref   = sol_exacte_chgt_perdir_der2(eta1(i),eta2(j))
        difference(i,j) = ref-node_val
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
           normH1(k)    = normH1(k) + ((grad1_node_val-grad1ref)**2+(grad2_node_val-grad2ref)**2)*h1*h2*val_jac
           
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
  
  call solve_fields( SLL_DIRICHLET, SLL_DIRICHLET, SLL_DIRICHLET, SLL_DIRICHLET, ti(k), te(k))
  
  do j=1,npts2
     do i=1,npts1
        
        node_val = calculated(i,j)
        grad1_node_val = phi%first_deriv_eta1_value_at_point(eta1(i), eta2(j))
        grad2_node_val = phi%first_deriv_eta2_value_at_point(eta1(i), eta2(j))
        ref        = sol_exacte_chgt_dirdir(eta1(i),eta2(j))
        grad1ref   = sol_exacte_chgt_dirdir_der1(eta1(i),eta2(j))
        grad2ref   = sol_exacte_chgt_dirdir_der2(eta1(i),eta2(j))
        difference(i,j) = ref-node_val
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
           normH1(k)    = normH1(k) + ((grad1_node_val-grad1ref)**2+(grad2_node_val-grad2ref)**2)*h1*h2*val_jac
           
        end if
     end do
  end do
  call delete_things()

  call check_error()

  case(8)
  print*, "---------------------"
  print*, " 8 test case with colella change of coordinates"
  print*, " dirichlet-periodic boundary conditions"
  print*, "---------------------"
  
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
  
  call solve_fields( SLL_DIRICHLET, SLL_DIRICHLET, SLL_PERIODIC, SLL_PERIODIC, ti(k), te(k))
  
  do j=1,npts2
     do i=1,npts1
        
        node_val = calculated(i,j)
        grad1_node_val = phi%first_deriv_eta1_value_at_point(eta1(i), eta2(j))
        grad2_node_val = phi%first_deriv_eta2_value_at_point(eta1(i), eta2(j))
        ref        = sol_exacte_chgt_dirper(eta1(i),eta2(j))
        grad1ref   = sol_exacte_chgt_dirper_der1(eta1(i),eta2(j))
        grad2ref   = sol_exacte_chgt_dirper_der2(eta1(i),eta2(j))
        difference(i,j) = ref-node_val
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
           normH1(k)    = normH1(k) + ((grad1_node_val-grad1ref)**2+(grad2_node_val-grad2ref)**2)*h1*h2*val_jac
        end if
     end do
  end do

  call delete_things()

  call check_error()

  case(9)
  print*, "---------------------"
  print*, " 95 test case without change of coordinates"
  print*, " periodic-periodic boundary conditions"
  print*, " with non analytic source term " 
  print*, "---------------------"
  
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
       NUM_CELLS1,&
       eta2,&
       NUM_CELLS2)  
  call rho%set_field_data(tab_rho)
  call rho%update_interpolation_coefficients()

  call solve_fields( SLL_PERIODIC, SLL_PERIODIC, SLL_PERIODIC, SLL_PERIODIC, ti(k), te(k))
  
  do j=1,npts2
     do i=1,npts1
        
        node_val = calculated(i,j)
        grad1_node_val = phi%first_deriv_eta1_value_at_point(eta1(i), eta2(j))
        grad2_node_val = phi%first_deriv_eta2_value_at_point(eta1(i), eta2(j))
        ref        = sol_exacte_perper(eta1(i),eta2(j))
        grad1ref   = sol_exacte_perper_der1(eta1(i),eta2(j))
        grad2ref   = sol_exacte_perper_der2(eta1(i),eta2(j))
        difference(i,j) = ref-node_val
        reference(i,j) = ref
        val_jac = 1.0
        if(PRINT_COMPARISON) call printout_comparison()
        if ( i < NUM_CELLS1 .and. j < NUM_CELLS2) then
           integral_solution = integral_solution + node_val * h1*h2
           integral_exact_solution = integral_exact_solution + ref * h1*h2
           normL2(k)    = normL2(k) + (node_val-ref)**2*h1*h2*val_jac
           normH1(k)    = normH1(k) + ((grad1_node_val-grad1ref)**2+(grad2_node_val-grad2ref)**2)*h1*h2*val_jac
        end if
     end do
  end do

  call delete_things()

  call check_error()

  case(10)
  print*, "---------------------"
  print*, " 9 test case with colella change of coordinates"
  print*, " periodic-periodic boundary conditions"
  print*, " with non analytic source term " 
  print*, "---------------------"
  
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
  
  rhs_interp => interp_2d_term_source
  
  tab_rho(:,:) = tab_rho - sum(tab_rho)/(NUM_CELLS1*NUM_CELLS2)

  rho => new_scalar_field_2d_discrete_alt( &
       "rho9", &
       rhs_interp, &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC,&
       SLL_PERIODIC,&
       SLL_PERIODIC ,&
       eta1,&
       NUM_CELLS1,&
       eta2,&
       NUM_CELLS2)

  call rho%set_field_data(tab_rho)
  call rho%update_interpolation_coefficients()
  
  
  call solve_fields( SLL_PERIODIC, SLL_PERIODIC, SLL_PERIODIC, SLL_PERIODIC, ti(k), te(k))
  
  do j=1,npts2
     do i=1,npts1
        
        node_val = calculated(i,j)
        grad1_node_val = phi%first_deriv_eta1_value_at_point(eta1(i), eta2(j))
        grad2_node_val = phi%first_deriv_eta2_value_at_point(eta1(i), eta2(j))
        ref        = sol_exacte_chgt_perper(eta1(i),eta2(j))
        grad1ref   = sol_exacte_chgt_perper_der1(eta1(i),eta2(j))
        grad2ref   = sol_exacte_chgt_perper_der2(eta1(i),eta2(j))
        difference(i,j) = ref-node_val
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
           normH1(k)    = normH1(k) + ((grad1_node_val-grad1ref)**2+(grad2_node_val-grad2ref)**2)*h1*h2*val_jac
           
        end if
     end do
  end do
  
  call delete_things()

  call check_error()

  case(11)
  print*, "---------------------"
  print*, " 10 test case with colella change of coordinates"
  print*, " periodic-dirichlet boundary conditions"
  print*, " with non analytic source term " 
  print*, "---------------------"
  
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

  rhs_interp => interp_2d_term_source

  rho => new_scalar_field_2d_discrete_alt( &
       "rho10", &
       rhs_interp, &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC,&
       SLL_DIRICHLET,&
       SLL_DIRICHLET,&
       eta1,&
       NUM_CELLS1,&
       eta2,&
       npts2)

  call rho%set_field_data(tab_rho)
  call rho%update_interpolation_coefficients()

  
  call solve_fields( SLL_PERIODIC, SLL_PERIODIC, SLL_DIRICHLET, SLL_DIRICHLET, ti(k), te(k))
  
  do j=1,npts2
     do i=1,npts1
        
        node_val = calculated(i,j)
        grad1_node_val = phi%first_deriv_eta1_value_at_point(eta1(i), eta2(j))
        grad2_node_val = phi%first_deriv_eta2_value_at_point(eta1(i), eta2(j))
        ref        = sol_exacte_chgt_perdir(eta1(i),eta2(j))
        grad1ref   = sol_exacte_chgt_perdir_der1(eta1(i),eta2(j))
        grad2ref   = sol_exacte_chgt_perdir_der2(eta1(i),eta2(j))
        difference(i,j) = ref-node_val
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
  
  call check_error()

  case(12)
  print*, "---------------------"
  print*, " 11 test case with colella change of coordinates"
  print*, " dirichlet-dirichlet boundary conditions"
  print*, " with non analytic source term " 
  print*, "---------------------"
  
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

  rhs_interp => interp_2d_term_source

  rho => new_scalar_field_2d_discrete_alt( &
       "rho11", &
       rhs_interp, &
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
  
  call solve_fields( SLL_DIRICHLET, SLL_DIRICHLET, SLL_DIRICHLET, SLL_DIRICHLET, ti(k), te(k))
  
  do j=1,npts2
     do i=1,npts1
        
        node_val = calculated(i,j)
        grad1_node_val = phi%first_deriv_eta1_value_at_point(eta1(i), eta2(j))
        grad2_node_val = phi%first_deriv_eta2_value_at_point(eta1(i), eta2(j))
        ref        = sol_exacte_chgt_dirdir(eta1(i),eta2(j))
        grad1ref   = sol_exacte_chgt_dirdir_der1(eta1(i),eta2(j))
        grad2ref   = sol_exacte_chgt_dirdir_der2(eta1(i),eta2(j))
        difference(i,j) = ref-node_val
        reference(i,j)  = ref

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
           normH1(k)    = normH1(k) + ((grad1_node_val-grad1ref)**2+&
                (grad2_node_val-grad2ref)**2)*h1*h2*val_jac
        end if
     end do
  end do

  call delete_things()

  call check_error()

  case(13)
  print*, "---------------------"
  print*, " 12 test case with colella change of coordinates"
  print*, " dirichlet-periodic  boundary conditions"
  print*, " with non analytic source term " 
  print*, "---------------------"
  
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
  
  rhs_interp => interp_2d_term_source

  rho => new_scalar_field_2d_discrete_alt( &
       "rho12", &
       rhs_interp, &
       T, &
       SLL_DIRICHLET,&
       SLL_DIRICHLET,&
       SLL_PERIODIC, &
       SLL_PERIODIC,&
       eta1,&
       npts1,&
       eta2,&
       NUM_CELLS2)

  call rho%set_field_data(tab_rho)
  call rho%update_interpolation_coefficients()

  call solve_fields( SLL_DIRICHLET, SLL_DIRICHLET, SLL_PERIODIC, SLL_PERIODIC, ti(k), te(k))
   
  do j=1,npts2
     do i=1,npts1
        
        node_val = calculated(i,j)
        grad1_node_val = phi%first_deriv_eta1_value_at_point(eta1(i), eta2(j))
        grad2_node_val = phi%first_deriv_eta2_value_at_point(eta1(i), eta2(j))
        ref        = sol_exacte_chgt_dirper(eta1(i),eta2(j))
        grad1ref   = sol_exacte_chgt_dirper_der1(eta1(i),eta2(j))
        grad2ref   = sol_exacte_chgt_dirper_der2(eta1(i),eta2(j))
        difference(i,j) = ref-node_val
        reference(i,j)  = ref
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

  call check_error()

  end select
  end do

  acc = acc/(npts1*npts2)

  print*,'error (per-per) with identity   =',acc(1)
  print*,'error (per-dir) with identity   =',acc(2)
  print*,'error (dir-dir) with identity   =',acc(3)
  print*,'error (dir-per) with identity   =',acc(4)
  print*,'error (per-per) with colella    =',acc(5)
  print*,'error (per-dir) with colella    =',acc(6)
  print*,'error (dir-dir) with colella    =',acc(7)
  print*,'error (dir-per) with colella    =',acc(8)
  print*,'error (per-per) with identity and source term non analytic=',acc(9)
  print*,'error (per-per) with colella  and source term non analytic=',acc(10)
  print*,'error (per-dir) with colella  and source term non analytic=',acc(11)
  print*,'error (dir-dir) with colella  and source term non analytic=',acc(12)
  print*,'error (dir-per) with colella  and source term non analytic=',acc(13)

  do k = 1, 13
    print"('test',i2,' : ','norm L2=',g15.3,' norm H1=',g15.3,' times=',2g15.3)" &
      ,k,normL2(k),normH1(k),ti(k),te(k)
  end do
  
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

  subroutine solve_fields( bc_eta1_min, bc_eta1_max, bc_eta2_min, bc_eta2_max, &
                           ti, te)

  sll_int32,  intent(in)  :: bc_eta1_min
  sll_int32,  intent(in)  :: bc_eta2_min
  sll_int32,  intent(in)  :: bc_eta1_max
  sll_int32,  intent(in)  :: bc_eta2_max
  sll_real64, intent(out) :: ti
  sll_real64, intent(out) :: te
  sll_int32               :: istep

  integral_solution = 0.0_f64
  integral_exact_solution = 0.0_f64

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
 
 
  call factorize_mat_es(&
    es, &
    a11_field_mat, &
    a12_field_mat,&
    a21_field_mat,&
    a22_field_mat,&
    b1_field_vect,&
    b2_field_vect,&
    c_field)

  ti = sll_time_elapsed_since(t_reference)

  call sll_set_time_mark(t_reference)

  write(*,"(' ')",advance="no")
  do istep = 1, 20

    values = 0.0_f64
    call phi%set_field_data(values)
    call phi%update_interpolation_coefficients()

    call sll_solve( es, rho, phi)
    write(*,"(i3)",advance="no") istep

  end do
  write(*,*) ' steps'

  te = sll_time_elapsed_since(t_reference)

  do j=1,npts2
     do i=1,npts1
        calculated(i,j) = phi%value_at_point(eta1(i),eta2(j))
     end do
  end do

  integral_solution = 0.0_f64
  integral_exact_solution = 0.0_f64

  end subroutine solve_fields

  subroutine check_error()
  print"('integral de la solution =',g15.3)", integral_solution
  print"('integral de la solution exacte =',g15.3)", integral_exact_solution
  acc = sum(abs(difference))
  print*, 'TEST ',k
  if ( ( sqrt(normL2(k)) <= h1**(SPLINE_DEG1-1))   .AND. &
       ( sqrt(normH1(k)) <= h1**(SPLINE_DEG1-1-1))) then     
     print *, 'PASSED'
  else
     stop     'FAILED'
  end if
  end subroutine check_error

end program test_general_elliptic_solver
