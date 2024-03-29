!> @ingroup quasi_neutral_solvers_parallel
!> @brief   Parallel 3D quasi-neutrality solver on "extruded" 2D polar mesh.
!> @authors Yaman Güçlü, IPP Garching
!> @authors Edoardo Zoni, IPP Garching
!>
!> @details
!> This module is a 3D wrapper around the 2D quasi-neutral solver in polar
!> coordinates, with simple cycle over the 3rd direction.

module sll_m_qn_solver_3d_polar_par
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_working_precision.h"

   use sll_m_remapper, only: &
      sll_t_layout_2d

   use sll_m_qn_solver_2d_polar_par, only: &
      sll_t_qn_solver_2d_polar_par, &
      sll_s_qn_solver_2d_polar_par_init, &
      sll_s_qn_solver_2d_polar_par_solve, &
      sll_s_qn_solver_2d_polar_par_free

   implicit none

   public :: &
      sll_t_qn_solver_3d_polar_par, &
      sll_s_qn_solver_3d_polar_par_init, &
      sll_s_qn_solver_3d_polar_par_solve, &
      sll_s_qn_solver_3d_polar_par_free

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   !> Class for the 3D quasi-neutral solver in polar coordinates
   !> It is basically a wrapper around a 2D solver, which is reused within
   !> a cycle over the x3 coordinate
   type sll_t_qn_solver_3d_polar_par

      type(sll_t_qn_solver_2d_polar_par), private :: solver_2d

   end type sll_t_qn_solver_3d_polar_par

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
contains
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

   !> Initialize the 3D (wrapper) solver
   subroutine sll_s_qn_solver_3d_polar_par_init(solver, &
                                                layout_r, &
                                                layout_a, &
                                                rmin, &
                                                rmax, &
                                                nr, &
                                                ntheta, &
                                                bc_rmin, &
                                                bc_rmax, &
                                                rho_m0, &
                                                b_magn, &
                                                lambda, &
                                                use_zonal_flow, &
                                                epsilon_0, &
                                                rgrid)

      type(sll_t_qn_solver_3d_polar_par), intent(inout) :: solver         !< Poisson solver class
      type(sll_t_layout_2d), pointer       :: layout_r       !< sequential in r direction
      type(sll_t_layout_2d), pointer       :: layout_a       !< sequential in theta direction
      real(f64), intent(in) :: rmin           !< rmin
      real(f64), intent(in) :: rmax           !< rmax
      integer(i32), intent(in) :: nr             !< number of cells radial
      integer(i32), intent(in) :: ntheta         !< number of cells angular
      integer(i32), intent(in) :: bc_rmin        !< boundary condition at r_min
      integer(i32), intent(in) :: bc_rmax        !< boundary condition at r_max
      real(f64), intent(in) :: rho_m0(:)      !< radial profile: total mass density of equilibrium
      real(f64), intent(in) :: b_magn(:)      !< radial profile: intensity of magnetic field
      real(f64), optional, intent(in) :: lambda(:)      !< radial profile: electron Debye length
      logical, optional, intent(in) :: use_zonal_flow !< if .false. set flux average to zero
      real(f64), optional, intent(in) :: epsilon_0      !< override default: vacuum permittivity
      real(f64), target, optional, intent(in) :: rgrid(:)       !< grid points along r

      ! Initialize the 2D solver
      call sll_s_qn_solver_2d_polar_par_init(solver%solver_2d, &
                                             layout_r, &
                                             layout_a, &
                                             rmin, &
                                             rmax, &
                                             nr, &
                                             ntheta, &
                                             bc_rmin, &
                                             bc_rmax, &
                                             rho_m0, &
                                             b_magn, &
                                             lambda, &
                                             use_zonal_flow, &
                                             epsilon_0, &
                                             rgrid)

   end subroutine sll_s_qn_solver_3d_polar_par_init

   !=============================================================================
   !> Solve the quasi-neutrality equation and get the electrostatic potential
   subroutine sll_s_qn_solver_3d_polar_par_solve(solver, rhs, phi)
      type(sll_t_qn_solver_3d_polar_par), intent(inout) :: solver     !< 3D solver
      real(f64), intent(in) :: rhs(:, :, :) !< Charge density
      real(f64), intent(out) :: phi(:, :, :) !< Potential

      integer(i32) :: i3

      SLL_ASSERT(all(shape(rhs) == shape(phi)))

      do i3 = lbound(rhs, 3), ubound(rhs, 3)
         call sll_s_qn_solver_2d_polar_par_solve(solver%solver_2d, rhs(:, :, i3), phi(:, :, i3))
      end do

   end subroutine sll_s_qn_solver_3d_polar_par_solve

   !=============================================================================
   !> Delete contents (local storage) of quasi-neutrality solver
   subroutine sll_s_qn_solver_3d_polar_par_free(solver)
      type(sll_t_qn_solver_3d_polar_par), intent(inout) :: solver

      call sll_s_qn_solver_2d_polar_par_free(solver%solver_2d)

   end subroutine sll_s_qn_solver_3d_polar_par_free

end module sll_m_qn_solver_3d_polar_par
