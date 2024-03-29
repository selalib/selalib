program test_scalar_field_1d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

   use sll_m_boundary_condition_descriptors, only: &
      sll_p_dirichlet, &
      sll_p_periodic

   use sll_m_cartesian_meshes, only: &
      sll_f_new_cartesian_mesh_1d, &
      sll_t_cartesian_mesh_1d

   use sll_m_constants, only: &
      sll_p_pi

   use sll_m_scalar_field_1d, only: &
      sll_t_scalar_field_1d_analytic, &
      sll_f_new_scalar_field_1d_analytic, &
      sll_f_new_scalar_field_1d_discrete

   use sll_m_scalar_field_1d_base, only: &
      sll_c_scalar_field_1d_base

   implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#define SPLINE_DEG1 3
#define NUM_CELLS1  64
#define ETA1MIN  0.0_f64
#define ETA1MAX  1.0_f64
#define PRINT_COMPARISON .false.

   type(sll_t_cartesian_mesh_1d), pointer                      :: mesh_1d

   class(sll_c_scalar_field_1d_base), pointer :: periodic_analytic
   class(sll_c_scalar_field_1d_base), pointer :: dirichlet_analytic
   class(sll_c_scalar_field_1d_base), pointer :: periodic_discrete
   class(sll_c_scalar_field_1d_base), pointer :: dirichlet_discrete

   sll_int32                                :: nc1!, iplot
   sll_real64                               :: grad1_node_val, grad1ref
   sll_real64, dimension(:), allocatable    :: tab_values
   sll_real64                               :: node_val, ref
   sll_real64, dimension(:), allocatable  :: point1
   sll_real64                               :: eta1
   sll_real64                               :: h1
   sll_int32                                :: i

   sll_real64 :: normL2_1, normL2_2, normL2_3, normL2_4
   sll_real64 :: normH1_1, normH1_2, normH1_3, normH1_4

   nc1 = NUM_CELLS1
   h1 = (ETA1MAX - ETA1MIN)/real(nc1, f64)

! First thing, initialize the logical mesh associated with this problem.
   mesh_1d => sll_f_new_cartesian_mesh_1d(NUM_CELLS1, ETA1MIN, ETA1MAX)

! --------------------------------------------------------------------------
!   Test case periodic analytic
!----------------------------------------------------------------------------

! ----> initialization of the field
   allocate (sll_t_scalar_field_1d_analytic :: periodic_analytic)
   select type (periodic_analytic)
   type is (sll_t_scalar_field_1d_analytic)
      call periodic_analytic%init( &
         test_function_per, &
         "periodic_analytic", &
         sll_p_periodic, &
         sll_p_periodic, &
         mesh_1d, &
         first_derivative=test_function_per_der1)
   end select

! -------> compute error norm L2 and H1
   normL2_1 = 0.0_f64
   normH1_1 = 0.0_f64
   do i = 1, nc1 + 1
      eta1 = real(i - 1, f64)*h1 + ETA1MIN
      node_val = periodic_analytic%value_at_point(eta1)
      grad1_node_val = periodic_analytic%derivative_value_at_point(eta1)
      ref = test_function_per(eta1)
      grad1ref = test_function_per_der1(eta1)
      if (PRINT_COMPARISON) then
         print *, 'eta1 = ', eta1, 'calculated = ', node_val, &
            'theoretical = ', ref, 'difference=', node_val - ref
         print *, 'eta1 = ', eta1, 'calculated = ', grad1_node_val, &
            'theoretical = ', grad1ref, 'difference=', grad1ref - grad1_node_val

      end if

      normL2_1 = normL2_1 + (node_val - ref)**2*h1
      normH1_1 = normH1_1 + ((grad1_node_val - grad1ref)**2)*h1

   end do
   print *, 'PASSED'
   call periodic_analytic%delete()

! --------------------------------------------------------------------------
!   Test case dirichlet analytic
!----------------------------------------------------------------------------

   dirichlet_analytic => sll_f_new_scalar_field_1d_analytic( &
                         test_function_dir, &
                         "dirichlet_analytic", &
                         sll_p_periodic, &
                         sll_p_periodic, &
                         mesh_1d, &
                         first_derivative=test_function_dir_der1)

   normL2_2 = 0.0_f64
   normH1_2 = 0.0_f64
   do i = 1, nc1 + 1
      eta1 = real(i - 1, f64)*h1 + ETA1MIN
      node_val = dirichlet_analytic%value_at_point(eta1)
      grad1_node_val = dirichlet_analytic%derivative_value_at_point(eta1)
      ref = test_function_dir(eta1)
      grad1ref = test_function_dir_der1(eta1)
      if (PRINT_COMPARISON) then
         print *, 'eta1 = ', eta1, 'calculated = ', node_val, &
            'theoretical = ', ref, 'difference=', node_val - ref
         print *, 'eta1 = ', eta1, 'calculated = ', grad1_node_val, &
            'theoretical = ', grad1ref, 'difference=', grad1ref - grad1_node_val

      end if

      normL2_2 = normL2_2 + (node_val - ref)**2*h1
      normH1_2 = normH1_2 + ((grad1_node_val - grad1ref)**2)*h1
   end do

   call dirichlet_analytic%delete()

! **********************  TESTS **************************************

   print *, '-------------------------------------------------------'
   print *, ' PERIODIC ANALYTIC'
   print *, '-------------------------------------------------------'
   print *, 'Norm L2', sqrt(normL2_1), 'Norm H1', sqrt(normH1_1), &
      h1**(SPLINE_DEG1), h1**(SPLINE_DEG1 - 1)

   print *, '-------------------------------------------------------'
   print *, ' DIRICHLET ANALYTIC'
   print *, '-------------------------------------------------------'
   print *, 'Norm L2', sqrt(normL2_2), 'Norm H1', sqrt(normH1_2), &
      h1**(SPLINE_DEG1), h1**(SPLINE_DEG1 - 1)

   if ((sqrt(normL2_1) <= h1**(SPLINE_DEG1)) .AND. &
       (sqrt(normL2_2) <= h1**(SPLINE_DEG1)) .AND. &
       (sqrt(normH1_1) <= h1**(SPLINE_DEG1 - 1)) .AND. &
       (sqrt(normH1_2) <= h1**(SPLINE_DEG1 - 1))) then
      print *, 'PASSED'
   end if

contains

   function test_function_per(eta1, params) result(res)
      sll_real64 :: res
      sll_real64, intent(in) :: eta1
      sll_real64, dimension(:), intent(in), optional :: params
      intrinsic :: cos
      res = cos(2*sll_p_pi*eta1)
   end function test_function_per

   function test_function_per_der1(eta1, params) result(res)
      intrinsic :: cos, sin
      sll_real64 :: res
      sll_real64, intent(in) :: eta1
      sll_real64, dimension(:), intent(in), optional :: params

      res = -2*sll_p_pi*sin(2*sll_p_pi*eta1)
   end function test_function_per_der1

   function test_function_dir(eta1, params) result(res)
      intrinsic :: cos, sin
      sll_real64 :: res
      sll_real64, intent(in) :: eta1
      sll_real64, dimension(:), intent(in), optional :: params

      res = sin(2*sll_p_pi*eta1)
   end function test_function_dir

   function test_function_dir_der1(eta1, params) result(res)
      intrinsic :: cos
      sll_real64 :: res
      sll_real64, intent(in) :: eta1
      sll_real64, dimension(:), intent(in), optional :: params

      res = 2.0*sll_p_pi*cos(2*sll_p_pi*eta1)
   end function test_function_dir_der1

end program test_scalar_field_1d

