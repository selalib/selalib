program test_dg_fields
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

   use sll_m_cartesian_meshes, only: &
      sll_f_new_cartesian_mesh_2d, &
      sll_t_cartesian_mesh_2d

   use sll_m_common_coordinate_transformations, only: &
      sll_f_identity_jac11, &
      sll_f_identity_jac12, &
      sll_f_identity_jac21, &
      sll_f_identity_jac22, &
      sll_f_identity_x1, &
      sll_f_identity_x2, &
      sll_f_sinprod_jac11, &
      sll_f_sinprod_jac12, &
      sll_f_sinprod_jac21, &
      sll_f_sinprod_jac22, &
      sll_f_sinprod_x1, &
      sll_f_sinprod_x2

   use sll_m_coordinate_transformation_2d_base, only: &
      sll_c_coordinate_transformation_2d_base, &
      sll_p_io_gmsh, &
      sll_p_io_gnuplot, &
      sll_p_io_mtv, &
      sll_p_io_xdmf

   use sll_m_coordinate_transformations_2d, only: &
      sll_f_new_coordinate_transformation_2d_analytic

   use sll_m_dg_fields, only: sll_t_dg_field_2d

   implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!=====================================!
! Simulation parameters               !
!=====================================!
   sll_int32, parameter :: nc_eta1 = 10  !
   sll_int32, parameter :: nc_eta2 = 20  !
   sll_int32, parameter :: degree = 3   !
!=====================================!

   sll_real64 :: eta1_max, eta1_min
   sll_real64 :: eta2_max, eta2_min
   sll_real64 :: delta_eta1, delta_eta2

   type(sll_t_cartesian_mesh_2d), pointer :: mesh
   class(sll_c_coordinate_transformation_2d_base), pointer :: tau
   class(sll_c_coordinate_transformation_2d_base), pointer :: collela

   type(sll_t_dg_field_2d) :: ex
   type(sll_t_dg_field_2d) :: bz

   mesh => sll_f_new_cartesian_mesh_2d(nc_eta1, nc_eta2, &
                                       eta1_min=-1._f64, eta1_max=1._f64, &
                                       eta2_min=-1._f64, eta2_max=1._f64)

   write (*, "(3f8.3,i4)") mesh%eta1_min, mesh%eta1_max, mesh%delta_eta1, mesh%num_cells1
   write (*, "(3f8.3,i4)") mesh%eta2_min, mesh%eta2_max, mesh%delta_eta2, mesh%num_cells2

   eta1_min = mesh%eta1_min
   eta1_max = mesh%eta1_max
   eta2_min = mesh%eta2_min
   eta2_max = mesh%eta2_max

   delta_eta1 = mesh%delta_eta1
   delta_eta2 = mesh%delta_eta2

! "Colella transformation";
! sinusoidal product (see P. Colella et al. JCP 230 (2011) formula
! (102) p 2968):
!
! x1 = eta1 + 0.1 * sin(2*pi*eta1) * sin(2*pi*eta2)
! x2 = eta2 + 0.1 * sin(2*pi*eta1) * sin(2*pi*eta2)

   collela => sll_f_new_coordinate_transformation_2d_analytic( &
              "collela_transformation", &
              mesh, &
              sll_f_sinprod_x1, &
              sll_f_sinprod_x2, &
              sll_f_sinprod_jac11, &
              sll_f_sinprod_jac12, &
              sll_f_sinprod_jac21, &
              sll_f_sinprod_jac22, &
              (/0.1_f64, 0.1_f64, 2.0_f64, 2.0_f64/))

   call collela%write_to_file(sll_p_io_gnuplot)

   tau => sll_f_new_coordinate_transformation_2d_analytic( &
          "identity_transformation", &
          mesh, &
          sll_f_identity_x1, &
          sll_f_identity_x2, &
          sll_f_identity_jac11, &
          sll_f_identity_jac12, &
          sll_f_identity_jac21, &
          sll_f_identity_jac22, &
          SLL_NULL_REAL64)

   call tau%write_to_file(sll_p_io_mtv)

   call ex%init(degree, tau, add)
   call ex%write_to_file('ex', sll_p_io_gmsh)
   call ex%write_to_file('ex', sll_p_io_mtv)
   call ex%write_to_file('ex', sll_p_io_xdmf)

   call bz%init(degree, collela, gaussian)
   call bz%write_to_file('bz', sll_p_io_gmsh)
   call bz%write_to_file('bz', sll_p_io_mtv)
   call bz%write_to_file('bz', sll_p_io_xdmf)

   print *, 'PASSED'

contains

   sll_real64 function gaussian(x1, x2, time)

      use sll_m_constants
      implicit none

      sll_real64, intent(in) :: x1
      sll_real64, intent(in) :: x2
      sll_real64, intent(in) :: time

      gaussian = exp(-(x1*x1 + x2*x2))*cos(time)
      return

   end function gaussian

   sll_real64 function add(x1, x2, time)

      use sll_m_constants
      implicit none

      sll_real64, intent(in) :: x1
      sll_real64, intent(in) :: x2
      sll_real64, intent(in) :: time

      add = x1 + x2
      return
      print *, time

   end function add

end program test_dg_fields
