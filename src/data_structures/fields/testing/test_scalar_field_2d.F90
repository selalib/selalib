program unit_test_2d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

   use helper_functions, only: &
      test_function_dirdir, &
      test_function_dirdir_der1, &
      test_function_dirdir_der2, &
      test_function_dirper, &
      test_function_dirper_der1, &
      test_function_dirper_der2, &
      test_function_perdir, &
      test_function_perdir_der1, &
      test_function_perdir_der2, &
      test_function_perper, &
      test_function_perper_der1, &
      test_function_perper_der2

   use sll_m_boundary_condition_descriptors, only: &
      sll_p_dirichlet, &
      sll_p_periodic

   use sll_m_cartesian_meshes

   use sll_m_common_coordinate_transformations, only: &
      sll_f_identity_jac11, &
      sll_f_identity_jac12, &
      sll_f_identity_jac21, &
      sll_f_identity_jac22, &
      sll_f_identity_x1, &
      sll_f_identity_x2

   use sll_m_coordinate_transformation_2d_base, only: &
      sll_c_coordinate_transformation_2d_base

   use sll_m_coordinate_transformations_2d

   use sll_m_scalar_field_2d, only: &
      sll_t_scalar_field_2d_analytic, &
      sll_t_scalar_field_2d_discrete

   use sll_m_scalar_field_2d_base, only: &
      sll_c_scalar_field_2d_base

   implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#define SPLINE_DEG1 3
#define SPLINE_DEG2 3
#define NUM_CELLS1  64
#define NUM_CELLS2  64
#define ETA1MIN  0.0_f64
#define ETA1MAX  1.0_f64
#define ETA2MIN  0.0_f64
#define ETA2MAX  1.0_f64
#define PRINT_COMPARISON .false.

   type(sll_t_cartesian_mesh_2d)                               :: mesh_2d
   class(sll_c_coordinate_transformation_2d_base), pointer :: T
   type(sll_t_coordinate_transformation_2d_analytic), target  :: T_analytic
   !type(sll_t_coordinate_transformation_2d_discrete),  target  :: T_discrete

   class(sll_c_scalar_field_2d_base), pointer :: scalar_field

   type(sll_t_scalar_field_2d_analytic), target :: doubly_periodic_analytic
   type(sll_t_scalar_field_2d_analytic), target :: periodic_dirichlet_analytic
   type(sll_t_scalar_field_2d_analytic), target :: dirichlet_dirichlet_analytic
   type(sll_t_scalar_field_2d_analytic), target :: dirichlet_periodic_analytic
   type(sll_t_scalar_field_2d_discrete), target :: doubly_periodic_discrete
   type(sll_t_scalar_field_2d_discrete), target :: periodic_dirichlet_discrete
   type(sll_t_scalar_field_2d_discrete), target :: dirichlet_dirichlet_discrete
   type(sll_t_scalar_field_2d_discrete), target :: dirichlet_periodic_discrete

   sll_int32 :: nc1, nc2!, iplot
   sll_real64 :: grad1_node_val, grad2_node_val, grad1ref, grad2ref
   sll_real64, dimension(:, :), pointer :: tab_values
   sll_real64 :: node_val, ref

   sll_real64, dimension(:), allocatable    :: point1
   sll_real64, dimension(:), allocatable    :: point2
   sll_real64 :: eta1, eta2
   sll_real64  :: h1, h2
   sll_int32 :: i, j

   sll_real64 :: normL2_1, normL2_2, normL2_3, normL2_4
   sll_real64 :: normL2_5, normL2_6, normL2_7, normL2_8
   sll_real64 :: normH1_1, normH1_2, normH1_3, normH1_4
   sll_real64 :: normH1_5, normH1_6, normH1_7, normH1_8

   sll_real64, dimension(1), parameter :: params_identity = [0.0_f64]

   ! logical mesh
   nc1 = NUM_CELLS1
   nc2 = NUM_CELLS1
   h1 = (ETA1MAX - ETA1MIN)/real(nc1, f64)
   h2 = (ETA2MAX - ETA2MIN)/real(nc2, f64)
   print *, 'h1 = ', h1
   print *, 'h2 = ', h2

   ! First thing, initialize the logical mesh associated with this problem.
   call sll_s_cartesian_mesh_2d_init(mesh_2d, NUM_CELLS1, NUM_CELLS2, &
                                     ETA1MIN, ETA1MAX, ETA2MIN, ETA2MAX)

   print *, 'initialized mesh 2D'

   ! coordinate transformation
   call sll_s_coordinate_transformation_2d_analytic_init( &
      T_analytic, &
      "analytic", &
      mesh_2d, &
      sll_f_identity_x1, &
      sll_f_identity_x2, &
      sll_f_identity_jac11, &
      sll_f_identity_jac12, &
      sll_f_identity_jac21, &
      sll_f_identity_jac22, &
      params_identity)
   print *, 'initialized transformation'
   T => T_analytic

   ! ******************************************************************
   ! ------------------ TEST ANALYTIC ------------------------------
   ! ******************************************************************

   ! --------------------------------------------------------------------------
   !   Test case periodic-periodic analytic
   !----------------------------------------------------------------------------

   ! ----> initialization of the field
   call doubly_periodic_analytic%init( &
      test_function_perper, &
      "doubly_periodic_analytic", &
      T, &
      sll_p_periodic, &
      sll_p_periodic, &
      sll_p_periodic, &
      sll_p_periodic, &
      (/0.0_f64/), & ! could be anything in this case
      first_deriv_eta1=test_function_perper_der1, &
      first_deriv_eta2=test_function_perper_der2)

   scalar_field => doubly_periodic_analytic

   print *, 'initialized field 2d'

   ! -------> compute error norm L2 and H1
   normL2_1 = 0.0_f64
   normH1_1 = 0.0_f64
   do j = 1, nc2 + 1
      do i = 1, nc1 + 1
         eta1 = real(i - 1, f64)*h1 + ETA1MIN
         eta2 = real(j - 1, f64)*h2 + ETA2MIN
         node_val = scalar_field%value_at_point(eta1, eta2)
         grad1_node_val = scalar_field%first_deriv_eta1_value_at_point(eta1, eta2)
         grad2_node_val = scalar_field%first_deriv_eta2_value_at_point(eta1, eta2)
         ref = test_function_perper(eta1, eta2, (/0.0_f64/))
         grad1ref = test_function_perper_der1(eta1, eta2, (/0.0_f64/))
         grad2ref = test_function_perper_der2(eta1, eta2, (/0.0_f64/))

         if (PRINT_COMPARISON) then
            print *, '(eta1,eta2) = ', eta1, eta2, 'calculated = ', node_val, &
               'theoretical = ', ref, 'difference=', node_val - ref
            print *, '(eta1,eta2) = ', eta1, eta2, 'calculated = ', grad1_node_val, &
               'theoretical = ', grad1ref, 'difference=', grad1ref - grad1_node_val

         end if

         normL2_1 = normL2_1 + (node_val - ref)**2*h1*h2
         normH1_1 = normH1_1 + ((grad1_node_val - grad1ref)**2 + &
                                (grad2_node_val - grad2ref)**2)*h1*h2

      end do
   end do

   ! -------> field visualization
   call scalar_field%write_to_file(1)

   ! -------> delete field
   call scalar_field%free()

   ! --------------------------------------------------------------------------
   !   Test case periodic-dirichlet analytic
   !----------------------------------------------------------------------------

   ! ----> initialization of the field
   call periodic_dirichlet_analytic%init( &
      test_function_perdir, &
      "periodic_dirichlet_analytic", &
      T, &
      sll_p_periodic, &
      sll_p_periodic, &
      sll_p_dirichlet, &
      sll_p_dirichlet, &
      (/0.0_f64/), & ! could be anything
      first_deriv_eta1=test_function_perdir_der1, &
      first_deriv_eta2=test_function_perdir_der2)

   scalar_field => periodic_dirichlet_analytic
   print *, 'initialized field 2d'

   ! -------> compute error norm L2 and H1
   normL2_2 = 0.0_f64
   normH1_2 = 0.0_f64
   do j = 1, nc2 + 1
      do i = 1, nc1 + 1
         eta1 = real(i - 1, f64)*h1 + ETA1MIN
         eta2 = real(j - 1, f64)*h2 + ETA2MIN
         node_val = scalar_field%value_at_point(eta1, eta2)
         grad1_node_val = scalar_field%first_deriv_eta1_value_at_point( &
                          eta1, eta2)
         grad2_node_val = &
            scalar_field%first_deriv_eta2_value_at_point( &
            eta1, eta2)
         ref = test_function_perdir(eta1, eta2, (/0.0_f64/))
         grad1ref = test_function_perdir_der1(eta1, eta2, (/0.0_f64/))
         grad2ref = test_function_perdir_der2(eta1, eta2, (/0.0_f64/))

         if (PRINT_COMPARISON) then
            print *, '(eta1,eta2) = ', eta1, eta2, 'calculated = ', node_val, &
               'theoretical = ', ref, 'difference=', node_val - ref
            print *, '(eta1,eta2) = ', eta1, eta2, 'calculated = ', grad1_node_val, &
               'theoretical = ', grad1ref, 'difference=', grad1ref - grad1_node_val

         end if

         normL2_2 = normL2_2 + (node_val - ref)**2*h1*h2
         normH1_2 = normH1_2 + ((grad1_node_val - grad1ref)**2 + &
                                (grad2_node_val - grad2ref)**2)*h1*h2

      end do
   end do

   ! -------> field visualization
   call scalar_field%write_to_file(1)

   ! -------> delete field
   call scalar_field%free()

   ! --------------------------------------------------------------------------
   !   Test case dirichlet-periodic analytic
   !----------------------------------------------------------------------------

   ! ----> initialization of the field
   call dirichlet_periodic_analytic%init( &
      test_function_dirper, &
      "dirichlet_periodic_analytic", &
      T, &
      sll_p_dirichlet, &
      sll_p_dirichlet, &
      sll_p_periodic, &
      sll_p_periodic, &
      (/0.0_f64/), &
      first_deriv_eta1=test_function_dirper_der1, &
      first_deriv_eta2=test_function_dirper_der2)

   print *, 'initialized field 2d'

   ! -------> compute error norm L2 and H1
   normL2_3 = 0.0_f64
   normH1_3 = 0.0_f64
   do j = 1, nc2 + 1
      do i = 1, nc1 + 1
         eta1 = real(i - 1, f64)*h1 + ETA1MIN
         eta2 = real(j - 1, f64)*h2 + ETA2MIN
         node_val = dirichlet_periodic_analytic%value_at_point(eta1, eta2)
         grad1_node_val = dirichlet_periodic_analytic%first_deriv_eta1_value_at_point( &
                          eta1, eta2)
         grad2_node_val = &
            dirichlet_periodic_analytic%first_deriv_eta2_value_at_point( &
            eta1, eta2)
         ref = test_function_dirper(eta1, eta2, (/0.0_f64/))
         grad1ref = test_function_dirper_der1(eta1, eta2, (/0.0_f64/))
         grad2ref = test_function_dirper_der2(eta1, eta2, (/0.0_f64/))

         if (PRINT_COMPARISON) then
            print *, '(eta1,eta2) = ', eta1, eta2, 'calculated = ', node_val, &
               'theoretical = ', ref, 'difference=', node_val - ref
            print *, '(eta1,eta2) = ', eta1, eta2, 'calculated = ', grad1_node_val, &
               'theoretical = ', grad1ref, 'difference=', grad1ref - grad1_node_val

         end if

         normL2_3 = normL2_3 + (node_val - ref)**2*h1*h2
         normH1_3 = normH1_3 + ((grad1_node_val - grad1ref)**2 + &
                                (grad2_node_val - grad2ref)**2)*h1*h2

      end do
   end do

   ! -------> field visualization
   call dirichlet_periodic_analytic%write_to_file(1)

   ! -------> delete field
   call dirichlet_periodic_analytic%free()

   ! --------------------------------------------------------------------------
   !   Test case dirichlet-dirichlet analytic
   !----------------------------------------------------------------------------

   ! ----> initialization of the field
   call dirichlet_dirichlet_analytic%init( &
      test_function_dirdir, &
      "dirichlet_dirichlet_analytic", &
      T, &
      sll_p_dirichlet, &
      sll_p_dirichlet, &
      sll_p_dirichlet, &
      sll_p_dirichlet, &
      (/0.0_f64/), &
      first_deriv_eta1=test_function_dirdir_der1, &
      first_deriv_eta2=test_function_dirdir_der2)

   print *, 'initialized field 2d'

   ! -------> compute error norm L2 and H1
   normL2_4 = 0.0_f64
   normH1_4 = 0.0_f64
   do j = 1, nc2 + 1
      do i = 1, nc1 + 1
         eta1 = real(i - 1, f64)*h1 + ETA1MIN
         eta2 = real(j - 1, f64)*h2 + ETA2MIN
         node_val = dirichlet_dirichlet_analytic%value_at_point(eta1, eta2)
         grad1_node_val = dirichlet_dirichlet_analytic%first_deriv_eta1_value_at_point( &
                          eta1, eta2)
         grad2_node_val = dirichlet_dirichlet_analytic%first_deriv_eta2_value_at_point( &
                          eta1, eta2)
         ref = test_function_dirdir(eta1, eta2, (/0.0_f64/))
         grad1ref = test_function_dirdir_der1(eta1, eta2, (/0.0_f64/))
         grad2ref = test_function_dirdir_der2(eta1, eta2, (/0.0_f64/))

         if (PRINT_COMPARISON) then
            print *, '(eta1,eta2) = ', eta1, eta2, 'calculated = ', node_val, &
               'theoretical = ', ref, 'difference=', node_val - ref
            print *, '(eta1,eta2) = ', eta1, eta2, 'calculated = ', grad1_node_val, &
               'theoretical = ', grad1ref, 'difference=', grad1ref - grad1_node_val

         end if

         normL2_4 = normL2_4 + (node_val - ref)**2*h1*h2
         normH1_4 = normH1_4 + ((grad1_node_val - grad1ref)**2 + &
                                (grad2_node_val - grad2ref)**2)*h1*h2

      end do
   end do
   ! -------> field visualization
   call dirichlet_dirichlet_analytic%write_to_file(1)

   ! -------> delete field
   call dirichlet_dirichlet_analytic%free()

   ! **********************  TESTS **************************************

   print *, '-------------------------------------------------------'
   print *, ' PERIODIC-PERIODIC ANALYTIC'
   print *, '-------------------------------------------------------'
   print *, 'Norm L2', sqrt(normL2_1), 'Norm H1', sqrt(normH1_1), &
      h1**(SPLINE_DEG1), h1**(SPLINE_DEG1 - 1)

   print *, '-------------------------------------------------------'
   print *, ' PERIODIC-DIRICHLET ANALYTIC'
   print *, '-------------------------------------------------------'
   print *, 'Norm L2', sqrt(normL2_2), 'Norm H1', sqrt(normH1_2), &
      h1**(SPLINE_DEG1), h1**(SPLINE_DEG1 - 1)

   print *, '-------------------------------------------------------'
   print *, ' DIRICHLET-PERIODIC ANALYTIC'
   print *, '-------------------------------------------------------'
   print *, 'Norm L2', sqrt(normL2_3), 'Norm H1', sqrt(normH1_3), &
      h1**(SPLINE_DEG1), h1**(SPLINE_DEG1 - 1)

   print *, '-------------------------------------------------------'
   print *, ' DIRICHLET-DIRICHLET ANALYTIC'
   print *, '-------------------------------------------------------'
   print *, 'Norm L2', sqrt(normL2_4), 'Norm H1', sqrt(normH1_4), &
      h1**(SPLINE_DEG1), h1**(SPLINE_DEG1 - 1)

   if ((sqrt(normL2_1) <= h1**(SPLINE_DEG1)) .AND. &
       (sqrt(normL2_2) <= h1**(SPLINE_DEG1)) .AND. &
       (sqrt(normL2_3) <= h1**(SPLINE_DEG1)) .AND. &
       (sqrt(normL2_4) <= h1**(SPLINE_DEG1)) .AND. &
       (sqrt(normH1_1) <= h1**(SPLINE_DEG1 - 1)) .AND. &
       (sqrt(normH1_2) <= h1**(SPLINE_DEG1 - 1)) .AND. &
       (sqrt(normH1_3) <= h1**(SPLINE_DEG1 - 1)) .AND. &
       (sqrt(normH1_4) <= h1**(SPLINE_DEG1 - 1))) then

      print *, 'PASSED'

   end if

end program unit_test_2d

