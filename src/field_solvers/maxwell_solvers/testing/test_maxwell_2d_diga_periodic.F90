! Here we solve Maxwell system without sources using
! discontnuous Galerkin numerical method on cartesian mesh
! Initiallly, Bz is set to cos product and boundary conditions are periodic.
! Ex and Ey are equals to analytic solutions.
! For time scheme we use classic RK4
!>$L_x,L_y$ domain dimensions and M,N are integers.
!>$
!>\omega = \sqrt{(\frac{M\pi}{L_x})^2+(\frac{N\pi}{L_y})^2}
!>$
!>$
!>B_z(x,y,t) =   - \cos(M \pi \frac{x}{L_x})  \cos(N \pi \frac{y}{L_y}) \cos(\omega t)
!>$
!>$
!>E_x(x,y,t) = \frac{c^2 N \pi }{\omega Ly} cos(M \pi \frac{x}{L_x}) \sin(N \pi  \frac{y}{L_y}) \sin(\omega t)
!>$
!>$
!>E_y(x,y,t) = - \frac{c^2 M \pi }{\omega Lx} \sin (M \pi \frac{x}{L_x}) \cos (N \pi  \frac{y}{L_y}) \sin(\omega t)
!>$
!>
program test_maxwell_2d_diga_periodic
!--------------------------------------------------------------------------
!  test 2D Maxwell solver based on discontinuous galerkine on a mapped mesh
!--------------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_maxwell_solvers_macros.h"

   use m_maxwell_helper_functions, only: &
      sol_bz, &
      sol_ex, &
      sol_ey

   use sll_m_boundary_condition_descriptors, only: &
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
      sll_c_coordinate_transformation_2d_base, &
      sll_p_io_mtv

   use sll_m_coordinate_transformations_2d

   use sll_m_dg_fields, only: &
      sll_t_dg_field_2d

   use sll_m_maxwell_2d_diga, only: &
      sll_s_maxwell_2d_diga_init, &
      sll_t_maxwell_2d_diga, &
      sll_s_solve_maxwell_2d_diga, &
      sll_p_uncentered

   use sll_m_maxwell_solvers_base, only: &
      sll_s_plot_two_fields

   implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!=====================================!
! Simulation parameters               !
!=====================================!
   sll_int32, parameter :: nc_eta1 = 10  !
   sll_int32, parameter :: nc_eta2 = 10  !
   sll_int32, parameter :: degree = 5   !
!=====================================!

   sll_int32  :: nstep
   sll_real64 :: eta1_max, eta1_min
   sll_real64 :: eta2_max, eta2_min
   sll_real64 :: delta_eta1, delta_eta2

   type(sll_t_cartesian_mesh_2d)                              :: mesh
   class(sll_c_coordinate_transformation_2d_base), pointer :: tau
   type(sll_t_coordinate_transformation_2d_analytic), target  :: tau_analytic

   type(sll_t_maxwell_2d_diga)   :: maxwell_TE

   type(sll_t_dg_field_2d) :: ex, ex0, dx, sx
   type(sll_t_dg_field_2d) :: ey, ey0, dy, sy
   type(sll_t_dg_field_2d) :: bz, bz0, dz, sz
   type(sll_t_dg_field_2d) :: exact

   sll_real64  :: time
   sll_int32   :: istep
   sll_real64  :: dt
   sll_real64  :: cfl = 0.1_f64
   sll_real64  :: error

#ifdef DEBUG
   sll_int32 :: i, j
   sll_real64, dimension(nc_eta1, nc_eta2) :: f1
   sll_real64, dimension(nc_eta1, nc_eta2) :: f2
#endif

   call sll_s_cartesian_mesh_2d_init(mesh, nc_eta1, nc_eta2, &
                                     eta1_min=0._f64, eta1_max=1._f64, &
                                     eta2_min=0._f64, eta2_max=1._f64)

   write (*, "(3f8.3,i4)") mesh%eta1_min, mesh%eta1_max, mesh%delta_eta1, mesh%num_cells1
   write (*, "(3f8.3,i4)") mesh%eta2_min, mesh%eta2_max, mesh%delta_eta2, mesh%num_cells2

   eta1_min = mesh%eta1_min
   eta1_max = mesh%eta1_max
   eta2_min = mesh%eta2_min
   eta2_max = mesh%eta2_max

   delta_eta1 = mesh%delta_eta1
   delta_eta2 = mesh%delta_eta2

! "Identity transformation";
   call sll_s_coordinate_transformation_2d_analytic_init( &
      tau_analytic, &
      "identity_transformation", &
      mesh, &
      sll_f_identity_x1, &
      sll_f_identity_x2, &
      sll_f_identity_jac11, &
      sll_f_identity_jac12, &
      sll_f_identity_jac21, &
      sll_f_identity_jac22, &
      SLL_NULL_REAL64)

   tau => tau_analytic

   call tau%write_to_file(sll_p_io_mtv)

   time = 0.0_f64

   call ex%init(degree, tau, sol_ex)
   call ey%init(degree, tau, sol_ey)
   call bz%init(degree, tau, sol_bz)

   call ex0%init(degree, tau)
   call ey0%init(degree, tau)
   call bz0%init(degree, tau)

   call dx%init(degree, tau)
   call dy%init(degree, tau)
   call dz%init(degree, tau)

   call sx%init(degree, tau)
   call sy%init(degree, tau)
   call sz%init(degree, tau)

   call exact%init(degree, tau)

   dt = cfl/sqrt(1./(delta_eta1/(degree + 1))**2 + 1./(delta_eta2/(degree + 1))**2)
   nstep = 100

   call sll_s_maxwell_2d_diga_init(maxwell_TE, tau, degree, TE_POLARIZATION, &
                                   sll_p_periodic, sll_p_periodic, sll_p_periodic, sll_p_periodic, &
                                   sll_p_uncentered)

   do istep = 1, nstep !*** Loop over time

      write (*, "(' istep = ',I6,' time = ',g12.3,' sec')") istep, time

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

#ifdef DEBUG

      call exact%set_value(sol_bz, time)
      do j = 1, nc_eta2
      do i = 1, nc_eta1
         f1(i, j) = bz%array(1, 1, i, j)
         f2(i, j) = exact%array(1, 1, i, j)
      end do
      end do

      call sll_s_plot_two_fields('bz', nc_eta1, nc_eta2, f1, f2, istep, time)

#endif

   end do ! next time step

   call check_error_ex(time)
   call check_error_ey(time)
   call check_error_bz(time)

   print"(a)", 'PASSED'

contains

   subroutine check_error_ex(time)

      sll_real64, intent(in) :: time

      call exact%set_value(sol_ex, time)

      error = maxval(abs(ex%array - exact%array))
      write (*, "(' EX erreur = ',g15.5)") error
      if (error > 0.001) STOP 'FAILED'

   end subroutine check_error_ex

   subroutine check_error_ey(time)

      sll_real64, intent(in) :: time

      call exact%set_value(sol_ey, time)

      error = maxval(abs(ey%array - exact%array))
      write (*, "(' EY erreur = ',g15.5)") error
      if (error > 0.001) STOP 'FAILED'

   end subroutine check_error_ey

   subroutine check_error_bz(time)

      sll_real64, intent(in) :: time

      call exact%set_value(sol_bz, time)

      error = maxval(abs(bz%array - exact%array))
      write (*, "(' BZ erreur = ',g15.5)") error
      if (error > 0.001) STOP 'FAILED'

   end subroutine check_error_bz

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

end program test_maxwell_2d_diga_periodic

