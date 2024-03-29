! Sample computation with the following characteristics:
! - vlasov-poisson
! - 4D: x, y, vx, vy (or x1, x2, x3, x4) with arbitrary coordinate
!   transformation
!   in the x,y variables.
! - parallel

program sim_bsl_vp_2d2v_polar
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

   use sll_m_cartesian_meshes, only: &
      sll_f_new_cartesian_mesh_2d, &
      sll_t_cartesian_mesh_2d

   use sll_m_collective, only: &
      sll_s_boot_collective, &
      sll_s_halt_collective

   use sll_m_common_array_initializers, only: &
      sll_f_periodic_gaussian_initializer_4d

   use sll_m_common_coordinate_transformations, only: &
      sll_f_polar_jac11, &
      sll_f_polar_jac12, &
      sll_f_polar_jac21, &
      sll_f_polar_jac22, &
      sll_f_polar_x1, &
      sll_f_polar_x2

   use sll_m_constants, only: &
      sll_p_pi

   use sll_m_coordinate_transformation_2d_base, only: &
      sll_c_coordinate_transformation_2d_base

   use sll_m_coordinate_transformations_2d, only: &
      sll_f_new_coordinate_transformation_2d_analytic

   use sll_m_sim_bsl_vp_2d2v_polar, only: &
      sll_s_initialize_vp4d_polar, &
      sll_o_delete, &
      sll_t_simulation_4d_vp_polar

   implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   character(len=256)                  :: filename
   type(sll_t_simulation_4d_vp_polar)    :: simulation
   type(sll_t_cartesian_mesh_2d), pointer  :: mx
   type(sll_t_cartesian_mesh_2d), pointer  :: mv
   sll_real64                          :: params(6)
   sll_int32                           :: nargs

   class(sll_c_coordinate_transformation_2d_base), pointer :: transformation

   call sll_s_boot_collective() ! Wrap this up somewhere else

   nargs = command_argument_count()
   if (nargs > 0) then

      call get_command_argument(1, filename)
      call simulation%init_from_file(trim(filename))

   else

#define NPTS1 64
#define NPTS2 64
#define NPTS3 32
#define NPTS4 32

      ! logical mesh for space coordinates
      mx => sll_f_new_cartesian_mesh_2d(NPTS1, NPTS2, &
                                        eta1_min=2.0_f64, eta1_max=8.0_f64, &
                                        eta2_min=.0_f64, eta2_max=2.0_f64*sll_p_pi)

      ! logical mesh for velocity coordinates
      mv => sll_f_new_cartesian_mesh_2d(NPTS3, NPTS4, &
                                        eta1_min=-6.0_f64, eta1_max=6.0_f64, &
                                        eta2_min=-6.0_f64, eta2_max=6.0_f64)

   end if

   ! coordinate transformation associated with space coordinates
   transformation => sll_f_new_coordinate_transformation_2d_analytic( &
                     "analytic_polar_transformation", &
                     mx, &
                     sll_f_polar_x1, &
                     sll_f_polar_x2, &
                     sll_f_polar_jac11, &
                     sll_f_polar_jac12, &
                     sll_f_polar_jac21, &
                     sll_f_polar_jac22, &
                     (/0.0_f64/)) ! this particular transformation is not parametrizable

   ! sll_o_initialize simulation object with the above parameters
   call sll_s_initialize_vp4d_polar( &
      simulation, &
      mx, &
      mv, &
      transformation, &
      sll_f_periodic_gaussian_initializer_4d, &
      params)

!  function defined in  parallel_array_initializers/sll_m_common_array_initializers.F90
!  sll_f_periodic_gaussian_initializer_4d(x,y,vx,xy) =
!  val = alpha*exp(-0.5_f64*((x -xc )**2+(y -yc )**2)) + &
!        beta *exp(-0.5_f64*((vx-vxc)**2+(vy-vyc)**2))

   params(1) = 5.0_f64 !xc
   params(2) = 0.0_f64 !yc
   params(3) = 0.0_f64 !vxc
   params(4) = 0.0_f64 !vyc
   params(5) = 1.0_f64 !alpha
   params(6) = 0.0_f64 !beta

   simulation%num_iterations = 100 ! run 100 iterations
   simulation%dt = 0.1_f64 ! time step
   call simulation%run()
   call sll_o_delete(simulation)

   print *, 'PASSED'

   call sll_s_halt_collective()

end program sim_bsl_vp_2d2v_polar
