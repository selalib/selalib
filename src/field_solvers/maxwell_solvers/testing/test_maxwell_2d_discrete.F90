program test_maxwell_2d_discrete
!--------------------------------------------------------------------------
!  test 2D Maxwell solver based on discontinuous galerkine on a mapped mesh
!  with discrete coordinate transformation
!--------------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_maxwell_solvers_macros.h"

   use m_maxwell_helper_functions, only: &
      sol_bz

   use sll_m_boundary_condition_descriptors, only: &
      sll_p_hermite, &
      sll_p_periodic

   use sll_m_cartesian_meshes, only: &
      sll_f_new_cartesian_mesh_2d, &
      sll_t_cartesian_mesh_2d

   use sll_m_common_coordinate_transformations, only: &
      sll_f_deriv1_jacobian_polar_f, &
      sll_f_deriv_x1_polar_f_eta1, &
      sll_f_deriv_x1_polar_f_eta2, &
      sll_f_deriv_x2_polar_f_eta1, &
      sll_f_deriv_x2_polar_f_eta2, &
      sll_f_jacobian_polar_f, &
      sll_f_x1_polar_f, &
      sll_f_x2_polar_f

   use sll_m_coordinate_transformation_2d_base, only: &
      sll_c_coordinate_transformation_2d_base

   use sll_m_coordinate_transformations_2d, only: &
      sll_f_new_coordinate_transformation_2d_analytic, &
      sll_f_new_coordinate_transformation_2d_discrete

   use sll_m_cubic_spline_interpolator_2d, only: &
      sll_t_cubic_spline_interpolator_2d

   use sll_m_dg_fields, only: &
      sll_t_dg_field_2d

   use sll_m_maxwell_2d_diga, only: &
      sll_t_maxwell_2d_diga, &
      sll_s_maxwell_2d_diga_init, &
      sll_s_solve_maxwell_2d_diga, &
      sll_p_uncentered, &
      sll_f_new_maxwell_2d_diga

   implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#define NPTS1 33
#define NPTS2 33

   type(sll_t_cartesian_mesh_2d), pointer :: mesh
   class(sll_c_coordinate_transformation_2d_base), pointer :: tau_d ! discrete transf
   class(sll_c_coordinate_transformation_2d_base), pointer :: tau_a ! analytic transf
! for the discrete case...
   type(sll_t_cubic_spline_interpolator_2d)  :: x1_interp
   type(sll_t_cubic_spline_interpolator_2d)  :: x2_interp
   type(sll_t_cubic_spline_interpolator_2d)  :: j_interp
   sll_real64, dimension(:, :), allocatable :: x1_tab
   sll_real64, dimension(:, :), allocatable :: x2_tab
   sll_real64, dimension(:), allocatable :: x1_eta1_min, x1_eta1_max
   sll_real64, dimension(:), allocatable :: x2_eta1_min, x2_eta1_max
   sll_real64, dimension(:, :), allocatable :: jacs
   sll_int32  :: i, j
   sll_real64 :: eta1, eta2, h1, h2
   sll_real64, dimension(2) :: params   ! for the polar transformation

   type(sll_t_maxwell_2d_diga), pointer  :: maxwell_d
   type(sll_t_maxwell_2d_diga), pointer  :: maxwell_a
   sll_int32, parameter                  :: degree = 2
   type(sll_t_dg_field_2d)      :: ex_a, dx_a
   type(sll_t_dg_field_2d)      :: ey_a, dy_a
   type(sll_t_dg_field_2d)      :: bz_a, dz_a
   type(sll_t_dg_field_2d)      :: ex_d, dx_d
   type(sll_t_dg_field_2d)      :: ey_d, dy_d
   type(sll_t_dg_field_2d)      :: bz_d, dz_d

#define RMIN 0.1_f64
#define RMAX 1.0_f64

   params(:) = (/RMIN, RMAX/)

   print *, 'filling out discrete arrays for x1 and x2 '
   print *, 'needed in the discrete case'

   h1 = 1.0_f64/real(NPTS1 - 1, f64)
   h2 = 1.0_f64/real(NPTS2 - 1, f64)
   print *, 'h1 = ', h1
   print *, 'h2 = ', h2
   allocate (x1_tab(NPTS1, NPTS2))
   allocate (x2_tab(NPTS1, NPTS2))
   allocate (x1_eta1_min(NPTS2))
   allocate (x1_eta1_max(NPTS2))
   allocate (x2_eta1_min(NPTS2))
   allocate (x2_eta1_max(NPTS2))
   allocate (jacs(NPTS1, NPTS2))

   mesh => sll_f_new_cartesian_mesh_2d(NPTS1 - 1, NPTS2 - 1)

   do j = 0, NPTS2 - 1
      do i = 0, NPTS1 - 1
         eta1 = real(i, f64)*h1
         eta2 = real(j, f64)*h2
         x1_tab(i + 1, j + 1) = sll_f_x1_polar_f(eta1, eta2, params)
         x2_tab(i + 1, j + 1) = sll_f_x2_polar_f(eta1, eta2, params)
         jacs(i + 1, j + 1) = sll_f_jacobian_polar_f(eta1, eta2, params)
      end do
   end do

   do j = 0, NPTS2 - 1
      eta1 = 0.0_f64
      eta2 = real(j, f64)*h2
      x1_eta1_min(j + 1) = sll_f_deriv_x1_polar_f_eta1(eta1, eta2, params)
      x2_eta1_min(j + 1) = sll_f_deriv_x2_polar_f_eta1(eta1, eta2, params)
      eta1 = 1.0_f64
      x1_eta1_max(j + 1) = sll_f_deriv_x1_polar_f_eta1(eta1, eta2, params)
      x2_eta1_max(j + 1) = sll_f_deriv_x2_polar_f_eta1(eta1, eta2, params)
   end do

   print *, 'initializing the interpolators: '

   call x1_interp%init( &
      NPTS1, &
      NPTS2, &
      0.0_f64, &
      1.0_f64, &
      0.0_f64, &
      1.0_f64, &
      sll_p_hermite, &
      sll_p_periodic, &
      eta1_min_slopes=x1_eta1_min, &
      eta1_max_slopes=x1_eta1_max)

   call x2_interp%init( &
      NPTS1, &
      NPTS2, &
      0.0_f64, &
      1.0_f64, &
      0.0_f64, &
      1.0_f64, &
      sll_p_hermite, &
      sll_p_periodic, &
      eta1_min_slopes=x2_eta1_min, &
      eta1_max_slopes=x2_eta1_max)

   call j_interp%init( &
      NPTS1, &
      NPTS2, &
      0.0_f64, &
      1.0_f64, &
      0.0_f64, &
      1.0_f64, &
      sll_p_hermite, &
      sll_p_periodic, &
      const_eta1_min_slope=sll_f_deriv1_jacobian_polar_f(0.0_f64, 0.0_f64, params), &
      const_eta1_max_slope=sll_f_deriv1_jacobian_polar_f(1.0_f64, 0.0_f64, params))

   print *, 'Initialized interpolators...'

   tau_d => sll_f_new_coordinate_transformation_2d_discrete( &
            mesh, &
            "polar_discrete", &
            x1_interp, &
            x2_interp, &
            j_interp, &
            x1_tab, &
            x2_tab, &
            jacobians_node=jacs)

   call tau_d%write_to_file()

   maxwell_d => sll_f_new_maxwell_2d_diga(tau_d, degree, TE_POLARIZATION, &
                                          sll_p_periodic, sll_p_periodic, &
                                          sll_p_periodic, sll_p_periodic, &
                                          sll_p_uncentered)

   call ex_d%init(degree, tau_d)
   call ey_d%init(degree, tau_d)
   call bz_d%init(degree, tau_d, sol_bz)
   call bz_d%write_to_file('bz_d')

   call dx_d%init(degree, tau_d)
   call dy_d%init(degree, tau_d)
   call dz_d%init(degree, tau_d)

   call sll_s_solve_maxwell_2d_diga(maxwell_d, ex_d, ey_d, bz_d, dx_d, dy_d, dz_d)
   call dx_d%write_to_file('dx_d')
   call dy_d%write_to_file('dy_d')

   tau_a => sll_f_new_coordinate_transformation_2d_analytic( &
            "polar_analytic", &
            mesh, &
            sll_f_x1_polar_f, &
            sll_f_x2_polar_f, &
            sll_f_deriv_x1_polar_f_eta1, &
            sll_f_deriv_x1_polar_f_eta2, &
            sll_f_deriv_x2_polar_f_eta1, &
            sll_f_deriv_x2_polar_f_eta2, &
            [RMIN, RMAX])

   call tau_a%write_to_file()

   maxwell_a => sll_f_new_maxwell_2d_diga(tau_a, degree, TE_POLARIZATION, &
                                          sll_p_periodic, sll_p_periodic, &
                                          sll_p_periodic, sll_p_periodic, &
                                          sll_p_uncentered)

   call ex_a%init(degree, tau_a)
   call ex_a%write_to_file('ex_a')
   call ey_a%init(degree, tau_a)
   call ey_a%write_to_file('ey_a')
   call bz_a%init(degree, tau_a, sol_bz)
   call bz_a%write_to_file('bz_a')

   call dx_a%init(degree, tau_a)
   call dy_a%init(degree, tau_a)
   call dz_a%init(degree, tau_a)

   call sll_s_solve_maxwell_2d_diga(maxwell_a, ex_a, ey_a, bz_a, dx_a, dy_a, dz_a)
   call dx_a%write_to_file('dx_a')
   call dy_a%write_to_file('dy_a')

   print"(a)", 'PASSED'

end program test_maxwell_2d_discrete
