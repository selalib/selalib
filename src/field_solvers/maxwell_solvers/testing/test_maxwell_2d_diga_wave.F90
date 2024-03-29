program test_maxwell_2d_diga_wave
   !--------------------------------------------------------------------------
   !  test 2D Maxwell solver based on discontinuous galerkine on a mapped mesh
   !--------------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_maxwell_solvers_macros.h"

   use m_maxwell_helper_functions, only: &
      gaussian

   use sll_m_boundary_condition_descriptors, only: &
      sll_p_silver_muller

   use sll_m_cartesian_meshes, only: &
      sll_s_cartesian_mesh_2d_init, &
      sll_t_cartesian_mesh_2d

   use sll_m_common_coordinate_transformations, only: &
      sll_f_affine_jac11, &
      sll_f_affine_jac12, &
      sll_f_affine_jac21, &
      sll_f_affine_jac22, &
      sll_f_affine_x1, &
      sll_f_affine_x2, &
      sll_f_homography_jac11, &
      sll_f_homography_jac12, &
      sll_f_homography_jac21, &
      sll_f_homography_jac22, &
      sll_f_homography_x1, &
      sll_f_homography_x2, &
      sll_f_identity_jac11, &
      sll_f_identity_jac12, &
      sll_f_identity_jac21, &
      sll_f_identity_jac22, &
      sll_f_identity_x1, &
      sll_f_identity_x2, &
      sll_f_rubber_sheeting_jac11, &
      sll_f_rubber_sheeting_jac12, &
      sll_f_rubber_sheeting_jac21, &
      sll_f_rubber_sheeting_jac22, &
      sll_f_rubber_sheeting_x1, &
      sll_f_rubber_sheeting_x2, &
      sll_f_sinprod_jac11, &
      sll_f_sinprod_jac12, &
      sll_f_sinprod_jac21, &
      sll_f_sinprod_jac22, &
      sll_f_sinprod_x1, &
      sll_f_sinprod_x2

   use sll_m_constants, only: &
      sll_p_pi

   use sll_m_coordinate_transformation_2d_base, only: &
      sll_c_coordinate_transformation_2d_base, &
      sll_p_io_mtv

   use sll_m_coordinate_transformations_2d, only: &
      sll_t_coordinate_transformation_2d_analytic, &
      sll_s_coordinate_transformation_2d_analytic_init, &
      sll_f_new_coordinate_transformation_2d_analytic

   use sll_m_dg_fields, only: &
      sll_t_dg_field_2d, &
      sll_f_new_dg_field_2d

   use sll_m_maxwell_2d_diga, only: &
      sll_s_maxwell_2d_diga_init, &
      sll_t_maxwell_2d_diga, &
      sll_s_solve_maxwell_2d_diga, &
      sll_p_uncentered

   implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   !=====================================!
   ! Simulation parameters               !
   !=====================================!
   sll_int32, parameter :: nc_eta1 = 10  !
   sll_int32, parameter :: nc_eta2 = 10  !
   sll_int32, parameter :: degree = 4   !
   !=====================================!

   sll_int32  :: nstep
   sll_real64 :: eta1_max, eta1_min
   sll_real64 :: eta2_max, eta2_min
   sll_real64 :: delta_eta1, delta_eta2

   type(sll_t_cartesian_mesh_2d) :: mesh
   class(sll_c_coordinate_transformation_2d_base), pointer :: tau
   type(sll_t_coordinate_transformation_2d_analytic), target  :: tau_a

   type(sll_t_maxwell_2d_diga)   :: maxwell_TE

   type(sll_t_dg_field_2d) :: ex, ex0, dx, sx
   type(sll_t_dg_field_2d) :: ey, ey0, dy, sy
   type(sll_t_dg_field_2d) :: bz, bz0, dz, sz

   sll_real64  :: time
   sll_int32   :: istep
   sll_real64  :: dt
   sll_real64  :: cfl = 0.5_f64
   sll_int32   :: itest
   !character(len=4) :: cstep

   call sll_s_cartesian_mesh_2d_init(mesh, nc_eta1, nc_eta2, &
                                     eta1_min=-5._f64, eta1_max=5._f64, &
                                     eta2_min=-5._f64, eta2_max=5._f64)

   write (*, "(3f8.3,i4)") mesh%eta1_min, mesh%eta1_max, mesh%delta_eta1, mesh%num_cells1
   write (*, "(3f8.3,i4)") mesh%eta2_min, mesh%eta2_max, mesh%delta_eta2, mesh%num_cells2

   eta1_min = mesh%eta1_min
   eta1_max = mesh%eta1_max
   eta2_min = mesh%eta2_min
   eta2_max = mesh%eta2_max

   delta_eta1 = mesh%delta_eta1
   delta_eta2 = mesh%delta_eta2

   do itest = 1, 5

      select case (itest)

      case (1)

         print *, "Identity transformation"; 
         call sll_s_coordinate_transformation_2d_analytic_init( &
            tau_a, &
            "identity_transformation", &
            mesh, &
            sll_f_identity_x1, &
            sll_f_identity_x2, &
            sll_f_identity_jac11, &
            sll_f_identity_jac12, &
            sll_f_identity_jac21, &
            sll_f_identity_jac22, &
            SLL_NULL_REAL64)

      case (2)

         print *, "Affine transformation"; 
         ! x1 = (b-a)*(cos(e)*eta1-sin(e)*eta2) + a
         ! x2 = (d-c)*(sin(e)*eta1+cos(e)*eta2) + c

         call sll_s_coordinate_transformation_2d_analytic_init( &
            tau_a, &
            "affine_transformation", &
            mesh, &
            sll_f_affine_x1, &
            sll_f_affine_x2, &
            sll_f_affine_jac11, &
            sll_f_affine_jac12, &
            sll_f_affine_jac21, &
            sll_f_affine_jac22, &
            [0.0_f64, 1.0_f64, &
             0.0_f64, 1.0_f64, &
             0.25*sll_p_pi])

      case (3)

         print *, "Homography transformation"

         ! x1 = (a*eta1*b*eta2+c)/(g*eta1+h*eta2+1)
         ! x2 = (d*eta1*e*eta2+f)/(g*eta1+h*eta2+1)
         !
         !  a = fixed scale factor in x1 direction with scale x2 unchanged.
         !  b = scale factor in x1 direction proportional to x2 distance from origin.
         !  c = origin translation in x1 direction.
         !  d = scale factor in x2 direction proportional to x1 distance from origin.
         !  e = fixed scale factor in x2 direction with scale x1 unchanged.
         !  f = origin translation in x2 direction.
         !  g = proportional scale factors x1 and x2 in function of x1.
         !  h = proportional scale factors x1 and x2 in function of x2.

         call sll_s_coordinate_transformation_2d_analytic_init( &
            tau_a, &
            "homography_transformation", &
            mesh, &
            sll_f_homography_x1, &
            sll_f_homography_x2, &
            sll_f_homography_jac11, &
            sll_f_homography_jac12, &
            sll_f_homography_jac21, &
            sll_f_homography_jac22, &
            [01.0_f64, 00.2_f64, 00.0_f64, &
             -00.2_f64, 01.0_f64, 00.0_f64, &
             00.0_f64, 00.0_f64])

      case (4)

         print *, "Colella transformation"; 
         ! sinusoidal product (see P. Colella et al. JCP 230 (2011) formula
         ! (102) p 2968):
         !
         ! x1 = eta1 + 0.1 * sin(2*pi*eta1) * sin(2*pi*eta2)
         ! x2 = eta2 + 0.1 * sin(2*pi*eta1) * sin(2*pi*eta2)

         call sll_s_coordinate_transformation_2d_analytic_init( &
            tau_a, &
            "collela_transformation", &
            mesh, &
            sll_f_sinprod_x1, &
            sll_f_sinprod_x2, &
            sll_f_sinprod_jac11, &
            sll_f_sinprod_jac12, &
            sll_f_sinprod_jac21, &
            sll_f_sinprod_jac22, &
            [0.1_f64, 0.1_f64, eta1_max - eta1_min, eta2_max - eta2_min])

      case (5)

         print *, "Rubber-Sheeting transformation"

         ! x1 = a*eta1*eta2+b*eta1+c*eta2+d
         ! x2 = e*eta1*eta2+f*eta1+g*eta2+h

         call sll_s_coordinate_transformation_2d_analytic_init( &
            tau_a, &
            "rubber_sheeting_transformation", &
            mesh, &
            sll_f_rubber_sheeting_x1, &
            sll_f_rubber_sheeting_x2, &
            sll_f_rubber_sheeting_jac11, &
            sll_f_rubber_sheeting_jac12, &
            sll_f_rubber_sheeting_jac21, &
            sll_f_rubber_sheeting_jac22, &
            [00.0_f64, 01.0_f64, 00.2_f64, 00.0_f64, &
             00.0_f64, 00.2_f64, 01.0_f64, 00.0_f64])

      end select

      tau => tau_a

      call tau%write_to_file(sll_p_io_mtv)

      time = 0.0_f64

      call ex%init(degree, tau)
      call ey%init(degree, tau)
      call bz%init(degree, tau, gaussian)

      call ex0%init(degree, tau)
      call ey0%init(degree, tau)
      call bz0%init(degree, tau)

      call dx%init(degree, tau)
      call dy%init(degree, tau)
      call dz%init(degree, tau)

      call sx%init(degree, tau)
      call sy%init(degree, tau)
      call sz%init(degree, tau)

      dt = cfl/sqrt(1./(delta_eta1/(degree + 1))**2 + 1./(delta_eta2/(degree + 1))**2)
      nstep = ceiling(15.0/dt)

      call sll_s_maxwell_2d_diga_init(maxwell_TE, &
                                      tau, &
                                      degree, &
                                      TE_POLARIZATION, &
                                      sll_p_silver_muller, &  !South
                                      sll_p_silver_muller, &  !East
                                      sll_p_silver_muller, &  !North
                                      sll_p_silver_muller, &  !West
                                      sll_p_uncentered)      !Flux

      do istep = 1, nstep !*** Loop over time

         !write(*,"(' istep = ',I6,' time = ',g12.3,' sec')") istep, time

         call rksetup()

         call sll_s_solve_maxwell_2d_diga(maxwell_TE, ex, ey, bz, dx, dy, dz)

         call accumulate(1._f64/6._f64)
         call rkstage(0.5_f64)

         call sll_s_solve_maxwell_2d_diga(maxwell_TE, ex, ey, bz, dx, dy, dz)
         call accumulate(1._f64/3._f64)
         call rkstage(0.5_f64)

         call sll_s_solve_maxwell_2d_diga(maxwell_TE, ex, ey, bz, dx, dy, dz)
         call accumulate(1._f64/3._f64)
         call rkstage(1.0_f64)

         call sll_s_solve_maxwell_2d_diga(maxwell_TE, ex, ey, bz, dx, dy, dz)
         call accumulate(1._f64/6._f64)

         call rkstep()

         time = time + dt

         !call sll_s_int2string(istep, cstep)
         !call bz%write_to_file('bz')
         !call maxwell_TE%po%write_to_file('phi')
         !call bz%write_to_file('bz'//cstep, sll_p_io_xdmf, time)

      end do ! next time step

      if (sqrt(sum(bz%array*bz%array)) &
          /real(nc_eta1*nc_eta2*(degree + 1)*(degree + 1), f64) < 0.001) then
         print"(a)", 'PASSED'
      else
         stop 'FAILED'
      end if

      call tau%delete()

   end do

contains

   subroutine rksetup()

      sx%array = 0.0_f64
      sy%array = 0.0_f64
      sz%array = 0.0_f64

      ex0%array = ex%array
      ey0%array = ey%array
      bz0%array = bz%array

   end subroutine rksetup

   subroutine rkstage(coef)

      sll_real64, intent(in) :: coef

      ex%array = ex0%array + coef*dt*dx%array
      ey%array = ey0%array + coef*dt*dy%array
      bz%array = bz0%array + coef*dt*dz%array

   end subroutine rkstage

   subroutine accumulate(coef)

      sll_real64, intent(in) :: coef

      sx%array = sx%array + coef*dx%array
      sy%array = sy%array + coef*dy%array
      sz%array = sz%array + coef*dz%array

   end subroutine accumulate

   subroutine rkstep()

      ex%array = ex0%array + dt*sx%array
      ey%array = ey0%array + dt*sy%array
      bz%array = bz0%array + dt*sz%array

   end subroutine rkstep

end program test_maxwell_2d_diga_wave
