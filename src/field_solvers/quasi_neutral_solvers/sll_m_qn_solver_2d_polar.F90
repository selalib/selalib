!**************************************************************
!  Copyright INRIA
!  Authors :
!     CALVI project team
!
!  This code SeLaLib (for Semi-Lagrangian-Library)
!  is a parallel library for simulating the plasma turbulence
!  in a tokamak.
!
!  This software is governed by the CeCILL-B license
!  under French law and abiding by the rules of distribution
!  of free software.  You can  use, modify and redistribute
!  the software under the terms of the CeCILL-B license as
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info".
!**************************************************************

!> @ingroup quasi_neutral_solvers
!>
!> @brief
!> Serial quasi-neutrality solver on 2D polar mesh; uses FFT in theta and
!> 2nd-order FD in r.
!>
!> @authors Yaman Güçlü, IPP Garching
!> @authors Edoardo Zoni, IPP Garching
!>
!> @details
!>
!> #### Model equations ####
!>
!> This module solves the quasi-neutrality equation on a 2D polar mesh
!> \f$ (r,\theta) \in [r_\textrm{min},r_\textrm{max}]\times[0,2\pi) \f$.
!> The general quasi-neutrality equation is of the form
!> \f[
!> -\nabla_\perp\cdot
!> \bigg[\frac{\rho_{m0}}{\epsilon_0 B^2}\nabla_\perp\phi\bigg]
!> +\frac{1}{\lambda_D^2}\big[\phi-\chi\langle\phi\rangle_f\big]
!> = \frac{1}{\epsilon_0}\rho_{c1},
!> \f]
!> where \f$ \phi \f$ is the electrostatic potential,
!> \f$ \langle\phi\rangle_f \f$ is the flux-surface average of \f$ \phi \f$,
!> \f$ \rho_{c1} \f$ is the charge perturbation density due to the
!> kinetic species, \f$ \rho_{m0} \f$ is the equilibrium mass density,
!> \f$ B \f$ is the intensity of the equilibrium background magnetic field,
!> \f$ \epsilon_0 \f$ is the permittivity of free space and \f$ \lambda_D \f$ is
!> the electron Debye length
!> \f[
!> \lambda_D\equiv\sqrt{\frac{\epsilon_0\kappa T_e}{q_e^{2}n_{e0}}}.
!> \f]
!> The equation in polar coordinates reads
!> \f[
!> -\bigg[
!> \frac{g}{r}\frac{\partial}{\partial r}
!> +\frac{\partial}{\partial r}\bigg(g\frac{\partial}{\partial r}\bigg)
!> +\frac{g}{r^2}\frac{\partial^2}{\partial\theta^2}
!> \bigg]\phi(r,\theta)
!> +\frac{1}{\lambda_D^2}\big[\phi(r,\theta)-\chi\langle\phi\rangle_\theta(r)\big]
!> = \frac{1}{\epsilon_0}\rho_{c1}(r,\theta),
!> \f]
!> where \f$ g\equiv\rho_{m0}/\epsilon_0 B^{2} \f$ is assumed to be independent
!> of \f$ \theta \f$ and the flux-surface average is replaced by a simple average
!> over \f$ \theta \f$.
!>
!> The boundary conditions (BCs) on \f$ \phi(r,\theta) \f$ are set as follows.
!> \f$ \phi(r,\theta) \f$ is \f$ 2\pi\f$-periodic along \f$ \theta \f$ and the
!> BCs along \f$ r \f$ can be chosen among the following types
!> (\f$ \overline{r} = r_\textrm{min} \f$ or \f$ \overline{r} = r_\textrm{max} \f$):
!>  <table>
!>  <tr><th> \c Option name <th> Description <th> Condition on solution <th> Condition on Fourier modes
!>  <tr><td> \c sll_p_dirichlet <td> Homogeneous Dirichlet
!>      <td> \f$ \phi(\overline{r},\theta)\equiv 0\f$
!>      <td> \f$ \widehat{\phi}_k(\overline{r})= 0\quad \forall k \f$
!>  <tr><td> \c sll_p_neumann <td> Homogeneous Neumann
!>      <td> \f$ \partial_r\phi(\overline{r},\theta)\equiv 0\f$
!>      <td> \f$ \partial_r\widehat{\phi}_k(\overline{r})=0\quad \forall k \f$
!>  <tr><td> \c sll_p_neumann_mode_0 <td> Homogeneous Neumann mode 0
!>      <td> \f$ \partial_\theta\phi(\overline{r},\theta)\equiv 0 \f$ and
!>           \f$ \int_0^{2\pi} \partial_r\phi(\overline{r},\theta) d\theta = 0 \f$
!>      <td> \f$ \partial_r\widehat{\phi}_0(\overline{r})=0 \f$ and
!>           \f$ \widehat{\phi}_k(\overline{r})=0 \f$ for \f$ k\neq 0 \f$
!>  <tr><td> \c sll_p_polar_origin   <td> \f$ \overline{r}=r_\textrm{min}=0 \f$ is center of circular disk
!>           <td> \f$ \phi(-\Delta r/2,\theta) \equiv \phi(\Delta r/2,\theta+\pi) \f$
!>           <td> \f$ \widehat{\phi}_k(-\Delta r/2) = (-1)^k \widehat{\phi}_k(\Delta r/2)\quad \forall k \f$
!>  </table>
!>
!> Please note that the option \c sll_p_polar_origin can only be used at
!> \f$ \overline{r} = r_\textrm{min} \f$ when one sets \f$ r_\textrm{min} = 0 \f$.
!>
!> The following arguments are given by the user at initialization:
!> \f$ \rho_{m0} \f$, \f$ B \f$, \f$ \lambda_D \f$, \f$ \epsilon_0 \f$ and the
!> additional parameter \f$ \chi \f$ (default is \f$ \chi=1 \f$).
!> If \f$ \lambda_D \f$ is not given to the solver, we assume
!> \f$ 1/\lambda_D^2=0 \f$: the electrons form a kinetic species and their
!> contribution goes into \f$ \rho_{c1} \f$. If \f$ \epsilon_0 \f$ is not given
!> to the solver, we assume \f$ \epsilon_0=8.854187817\times10^{-12} [F/m]\f$
!> according to the parameter \c sll_p_epsilon_0 in module sll_m_constants.
!>
!> #### Numerical methods ####
!>
!> This module uses fast Fourier transforms (FFTs) in \f$ \theta \f$ and a
!> 2nd-order finite-difference method in \f$ r \f$.
!> Thanks to the linearity of the differential operator and the periodicity of
!> the domain, a discrete Fourier transform (DFT) in \f$ \theta \f$ is applied
!> to both sides of the above elliptic PDE. Then, each Fourier coefficient
!> \f$ \widehat{\phi}_k(r) \f$ solves an independent 1D boundary value problem
!> on \f$ [r_\textrm{min},r_\textrm{max}] \f$:
!> \f[
!> -\bigg[
!> \frac{g}{r}\frac{\partial}{\partial r}
!> +\frac{\partial}{\partial r}\bigg(g\frac{\partial}{\partial r}\bigg)
!> -\frac{k^2}{r^2}g
!> \bigg]\widehat{\phi}_k(r)
!> +\frac{1}{\lambda_D^2}(1-\chi\,\delta_{k0})\widehat{\phi}_k(r)
!> = \frac{1}{\epsilon_0}\widehat{\rho}_{c1,k}(r),
!> \f]
!> For each mode \f$ k \f$, the resulting ODE is solved with a 2nd-order
!> finite-difference collocation method.
!>
!> Since \f$ \phi(r,\theta) \f$ is real, we have
!> \f$ \widehat{\phi}_{-k}=\overline{\widehat{\phi}_k} \f$ for each
!> \f$ k=-N_\theta/2,-N_\theta/2+1,\dots,-1,0,1,\dots,N_\theta/2-1,N_\theta/2\f$,
!> where \f$ N_\theta\f$ is the number of cells along \f$ \theta \f$.
!> Therefore, it is enough to consider only \f$ N_\theta/2+1 \f$ coefficients
!> \f$ \widehat{\phi}_0, \dots, \widehat{\phi}_{N_\theta/2} \f$. For this purpose,
!> the \c r2c and \c c2r interfaces of FFTW are used for the forward and backward
!> FFTs, respectively.
!>
!> If the user does not provide a radial grid, a collocation grid with uniform
!> spacing \f$ \Delta r \f$ is calculated from the input parameters
!> \f$ r_\textrm{min} \f$, \f$ r_\textrm{max} \f$ and \f$N_r\f$
!> (number of radial cells).
!> In the case of a circular annulus (\f$ r_\textrm{min} > 0\f$),
!> the collocation grid has \f$ N_r+1 \f$ points, with first point
!> \f$ r_1 = r_\textrm{min} \f$ and last point \f$ r_{N_r+1} = r_\textrm{max} \f$.
!> In the case of a circular disc (\f$ r_\textrm{min} = 0 \f$),
!> the collocation grid has only \f$ N_r \f$ points, with first point
!> \f$ r_1 = \Delta r/2 \f$ and last point \f$ r_{N_r} = r_\textrm{max} \f$.
!>
!> For a given mode \f$k\f$, at each interior point in the collocation grid
!> (\f$ r_i \neq r_\textrm{min}, r_\textrm{max} \f$) the equation above
!> is approximated using centered finite-differences on a 3-point stencil,
!> valid on a general non-uniform grid with
!> \f$ h_i^{+}=r_{i+1}-r_i \f$ and \f$ h_i^{-}=r_{i}-r_{i-1} \f$:
!>
!> \f[
!> \frac{\partial}{\partial r} \widehat\phi_k(r_i) =
!> \frac{1}{h_i^+ + h_i^-} \left[
!>  -\frac{h_i^+}{h_i^-}\, \widehat{\phi}_k(r_{i-1})
!>  +\left( \frac{h_i^+}{h_i^-} - \frac{h_i^-}{h_i^+} \right)\, \widehat{\phi}_k(r_i)
!>  +\frac{h_i^-}{h_i^+}\, \widehat{\phi}_k(r_{i+1})
!> \right]
!>  \,+\, \left\{ -\frac{\partial^3}{\partial r^3} \widehat\phi_k(r_i) \frac{h_i^+ h_i^-}{6} \right\}
!>  \,+\, \textit{H.O.T.},
!> \f]
!>
!> \f[
!> \frac{\partial^2}{\partial r^2} \widehat\phi_k(r_i) =
!> \frac{1}{h_i^+ h_i^-} \left[
!>   \frac{2h_i^+}{h_i^++h_i^-}\, \widehat{\phi}_k(r_{i-1})
!>                             -2 \widehat{\phi}_k(r_i)
!>  +\frac{2h_i^-}{h_i^++h_i^-}\, \widehat{\phi}_k(r_{i+1})
!> \right]
!>  \,+\, \left\{ -\frac{\partial^3}{\partial r^3} \widehat\phi_k(r_i) \frac{h_i^+ - h_i^-}{3}
!>  \,            -\frac{\partial^4}{\partial r^4} \widehat\phi_k(r_i) \frac{(h_i^+)^2 - h_i^+ h_i^- +(h_i^-)^2}{12} \right\}
!>  \,+\, \textit{H.O.T.}.
!> \f]
!>
!> On the right-hand side of each equation, the first term is the
!> finite-difference approximation used in the code,
!> the terms in curly braces represent the leading truncation error,
!> and \f$ \textit{H.O.T} \f$ means "higher order terms".
!> The coefficients that multiply \f$\widehat{\phi}_k(r_{i-1})\f$, \f$\widehat{\phi}_k(r_i)\f$
!> and \f$\widehat{\phi}_k(r_{i+1})\f$ are the so-called "finite-difference coefficients".
!> Although both formulas are exact for a parabolic profile, for a general profile
!> the finite-difference approximation to the second derivative is only 1st-order
!> accurate when \f$ h_i^+ \neq h_i^- \f$.
!> In the future the finite-difference coefficients for the second derivative
!> can be modified to cancel the 1st-order error term: the Fourier-transformed
!> quasi-neutrality equation can be differentiated to obtain an expression for
!> \f$ \partial^3\widehat{\phi}_k / \partial r^3 \f$ that involves only
!> 1st and 2nd derivatives,
!> which are then approximated with the same finite-difference formulas.
!>
!> #### Usage example ####
!>
!> \code
!> use sll_m_qn_solver_2d_polar, only: &
!>   sll_t_qn_solver_2d_polar,         &
!>   sll_s_qn_solver_2d_polar_init,    &
!>   sll_s_qn_solver_2d_polar_solve,   &
!>   sll_s_qn_solver_2d_polar_free
!>
!> ...
!>
!> type(sll_t_qn_solver_2d_polar) :: solver
!> real(f64), allocatable         :: rho(:,:), phi(:,:)
!> real(f64), allocatable         :: rho_m0(:), b_magn(:), lambda(:)
!> real(f64)                      :: rmin, rmax
!> integer(i32)                   :: nr, ntheta, bc_rmin, bc_rmax
!>
!> ...
!>
!> call sll_s_qn_solver_2d_polar_init ( solver, rmin, rmax, nr, ntheta, bc_rmin, bc_rmax, rho_m0, b_magn, lambda, use_zonal_flow=.false., epsilon_0=1.0_f64 )
!> call sll_s_qn_solver_2d_polar_solve( solver, rho, phi )
!> call sll_s_qn_solver_2d_polar_free ( solver )
!> \endcode

module sll_m_qn_solver_2d_polar
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"
#include "sll_working_precision.h"

   use sll_m_constants, only: &
      sll_p_epsilon_0

   use sll_m_boundary_condition_descriptors, only: &
      sll_p_dirichlet, &
      sll_p_neumann, &
      sll_p_neumann_mode_0, &
      sll_p_polar_origin

   use sll_m_fft, only: &
      sll_t_fft, &
      sll_f_fft_allocate_aligned_complex, &
      sll_f_fft_allocate_aligned_real, &
      sll_s_fft_deallocate_aligned_complex, &
      sll_s_fft_deallocate_aligned_real, &
      sll_s_fft_init_r2c_1d, &
      sll_s_fft_init_c2r_1d, &
      sll_s_fft_exec_r2c_1d, &
      sll_s_fft_exec_c2r_1d, &
      sll_s_fft_free

   use sll_m_tridiagonal, only: &
      sll_s_setup_cyclic_tridiag, &
      sll_o_solve_cyclic_tridiag

   use sll_m_utilities, only: &
      sll_s_new_array_linspace

   implicit none

   public :: &
      sll_t_qn_solver_2d_polar, &
      sll_s_qn_solver_2d_polar_init, &
      sll_s_qn_solver_2d_polar_solve, &
      sll_s_qn_solver_2d_polar_free

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   !> Class for 2D Poisson solver in polar coordinates
   type sll_t_qn_solver_2d_polar

      real(f64)              :: rmin      !< Min value of r coordinate
      real(f64)              :: rmax      !< Max value of r coordinate
      integer(i32)              :: nr        !< Number of cells along r
      integer(i32)              :: nt        !< Number of cells along theta
      integer(i32)              :: bc(2)     !< Boundary conditions options
      real(f64)              :: epsilon_0 !< Vacuum permittivity (user may override)
      real(f64), allocatable :: g(:)     !< \f$ g(r) = \rho_{m,0}/(B^2\epsilon_0) \f$

      type(sll_t_fft)           :: fw        !< Forward FFT plan
      type(sll_t_fft)           :: bw        !< Inverse FFT plan
      complex(f64), allocatable :: z(:, :) !< 2D work array needed for transposition
      complex(f64), pointer     :: temp_c(:) !< 1D work array, complex
      real(f64), pointer     :: temp_r(:) !< 1D work array, real
      real(f64), allocatable :: mat(:, :) !< Tridiagonal matrix (one for each k)
      real(f64), allocatable :: cts(:)   !< Lapack coefficients
      integer(i32), allocatable :: ipiv(:)   !< Lapack pivot indices

      real(f64) :: bc_coeffs_rmin(2:3) ! needed for Neumann at r=r_min
      real(f64) :: bc_coeffs_rmax(-2:-1) ! needed for Neumann at r=r_max
      integer   :: skip0                 ! needed for full circle

   end type sll_t_qn_solver_2d_polar

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   !=============================================================================
   !> Initialize the solver
   subroutine sll_s_qn_solver_2d_polar_init(solver, &
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

      type(sll_t_qn_solver_2d_polar), intent(out) :: solver         !< solver object
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

      character(len=*), parameter :: this_sub_name = 'sll_s_qn_solver_2d_polar_init'

      real(f64) :: hp, hm
      real(f64) :: inv_r
      real(f64) :: d0_coeffs(3)
      real(f64) :: d1_coeffs(3)
      real(f64) :: d2_coeffs(3)
      real(f64) :: ddr_ln_g
      real(f64) :: c
      real(f64), allocatable :: r_nodes(:)

      integer(i32) :: i, k
      integer(i32) :: bck(2)
      integer(i32) :: last

      ! Consistency checks
      SLL_ASSERT_ALWAYS(rmin >= 0.0_f64)
      SLL_ASSERT_ALWAYS(rmin < rmax)
      SLL_ASSERT_ALWAYS(nr >= 1)
      SLL_ASSERT_ALWAYS(ntheta >= 1)

      if (bc_rmin == sll_p_polar_origin .and. rmin /= 0.0_f64) then
         SLL_ERROR(this_sub_name, "BC option 'sll_p_polar_origin' requires r_min = 0")
      end if

      ! Set boundary condition at r_min
      select case (bc_rmin)
      case (sll_p_dirichlet, sll_p_neumann, sll_p_neumann_mode_0, sll_p_polar_origin)
         solver%bc(1) = bc_rmin
      case default
         SLL_ERROR(this_sub_name, 'Unrecognized boundary condition at r_min')
      end select

      ! Set boundary condition at r_max
      select case (bc_rmax)
      case (sll_p_dirichlet, sll_p_neumann, sll_p_neumann_mode_0)
         solver%bc(2) = bc_rmax
      case default
         SLL_ERROR(this_sub_name, 'Unrecognized boundary condition at r_max')
      end select

      ! Override vacuum permittivity in SI units
      if (present(epsilon_0)) then
         solver%epsilon_0 = epsilon_0
      else
         solver%epsilon_0 = sll_p_epsilon_0
      end if

      ! Store information in solver
      solver%rmin = rmin
      solver%rmax = rmax
      solver%nr = nr
      solver%nt = ntheta

      ! Important: if full circle is simulated, center point is not solved for!
      !            Therefore, solution is calculated on nr+1-skip0 points.
      solver%skip0 = merge(1, 0, bc_rmin == sll_p_polar_origin)

      ! r grid (possibly non-uniform)
      allocate (r_nodes(nr + 1))

      if (present(rgrid)) then  !--> Create computational grid from user data

         SLL_ASSERT_ALWAYS(all(rgrid > 0.0_f64))
         if (bc_rmin == sll_p_polar_origin) then
            SLL_ASSERT_ALWAYS(size(rgrid) == nr)
            SLL_ASSERT_ALWAYS(rgrid(nr) == rmax)
            r_nodes(1) = -rgrid(1)
            r_nodes(2:) = rgrid(:)
         else
            SLL_ASSERT_ALWAYS(size(rgrid) == nr + 1)
            SLL_ASSERT_ALWAYS(rgrid(1) == rmin)
            SLL_ASSERT_ALWAYS(rgrid(nr + 1) == rmax)
            r_nodes(:) = rgrid(:)
         end if

      else  !-------------------------> Create uniform grid

         if (bc_rmin == sll_p_polar_origin) then
            associate (rmin => rmax/real(2*nr + 1, f64))
               r_nodes(1) = -rmin
               call sll_s_new_array_linspace(r_nodes(2:), rmin, rmax, endpoint=.true.)
            end associate
         else
            call sll_s_new_array_linspace(r_nodes, rmin, rmax, endpoint=.true.)
         end if

      end if

      ! Allocate arrays for solution of linear systems along r
      associate (sh => solver%skip0)
         allocate (solver%g(nr + 1))
         allocate (solver%z(nr + 1 - sh, 0:ntheta/2))
         allocate (solver%mat((nr - 1)*3, 0:ntheta/2)) ! for each k, matrix depends on r
         allocate (solver%cts((nr - 1)*7))
         allocate (solver%ipiv(nr - 1))
      end associate

      ! Allocate in ALIGNED fashion 1D work arrays for FFT
      solver%temp_r => sll_f_fft_allocate_aligned_real(ntheta)
      solver%temp_c => sll_f_fft_allocate_aligned_complex(ntheta/2 + 1)

      ! Initialize plans for forward and backward FFTs
      call sll_s_fft_init_r2c_1d(solver%fw, &
                                 ntheta, &
                                 solver%temp_r(:), &
                                 solver%temp_c(:), &
                                 aligned=.true., &
                                 normalized=.true.)

      call sll_s_fft_init_c2r_1d(solver%bw, &
                                 ntheta, &
                                 solver%temp_c(:), &
                                 solver%temp_r(:), &
                                 aligned=.true., &
                                 normalized=.false.)

      ! Store non-dimensional coefficient g(r) = \rho(r) / (B(r)^2 \epsilon_0)
      associate (sh => solver%skip0)
         solver%g(1 + sh:) = rho_m0(:)/(b_magn(:)**2*solver%epsilon_0)
         if (sh == 1) then
            solver%g(1) = solver%g(2) ! ghost point
         end if
      end associate

      ! Store matrix coefficients into solver%mat
      ! Cycle over k
      do k = 0, ntheta/2

         !--------------------------------------------
         ! Compute boundary conditions type for mode k
         !--------------------------------------------
         bck(:) = solver%bc(:)
         do i = 1, 2
            if (bck(i) == sll_p_neumann_mode_0) then
               if (k == 0) then
                  bck(i) = sll_p_neumann
               else
                  bck(i) = sll_p_dirichlet
               end if
            end if
         end do

         !--------------------------------------------
         ! Compute matrix coefficients for a given k_j
         !--------------------------------------------
         do i = 2, nr

            ! Finite difference coefficients for 0th (trivial), 1st and 2nd derivatives
            hp = r_nodes(i + 1) - r_nodes(i)
            hm = r_nodes(i) - r_nodes(i - 1)
            inv_r = 1.0_f64/r_nodes(i)
            d0_coeffs(:) = [0.0_f64, 1.0_f64, 0.0_f64]
            d1_coeffs(:) = [-hp/hm, (hp**2 - hm**2)/(hp*hm), hm/hp]/(hp + hm)
            d2_coeffs(:) = [2*hp/(hp + hm), -2.0_f64, 2*hm/(hp + hm)]/(hp*hm)

            ! Compute d/dr(ln(g)) as (1/g)*(dg/dr) using finite differences
            ddr_ln_g = dot_product(d1_coeffs(:), solver%g(i - 1:i + 1))/solver%g(i)

            ! Determine contribution from zonal flow
            if (present(lambda)) then
               associate (sh => solver%skip0)
                  c = 1.0_f64/(lambda(i - sh)**2*solver%g(i))
               end associate
               if (present(use_zonal_flow)) then
                  if (use_zonal_flow .and. k == 0) then
                     c = 0.0_f64
                  end if
               end if
            else
               c = 0.0_f64
            end if

            ! Fill in elements of i-th matrix row
            solver%mat(3*(i - 1) - 2:3*(i - 1), k) = &
               -d2_coeffs(:) &
               - d1_coeffs(:)*(inv_r + ddr_ln_g) &
               + d0_coeffs(:)*((k*inv_r)**2 + c)

         end do

         !--------------------------------------------
         ! Set boundary condition at rmin
         !--------------------------------------------
         if (bck(1) == sll_p_dirichlet) then ! Dirichlet
            solver%mat(1, k) = 0.0_f64

         else if (bck(1) == sll_p_neumann) then ! Neumann

            ! Coefficients of homogeneous boundary condition
            hp = r_nodes(3) - r_nodes(2)
            hm = r_nodes(2) - r_nodes(1)
            d1_coeffs(:) = [-2 - hp/hm, 2 + hp/hm + hm/hp, -hm/hp]
            solver%bc_coeffs_rmin(2:3) = -d1_coeffs(2:3)/d1_coeffs(1)

            ! Gaussian elimination: remove phi(1) variable using boundary condition
            solver%mat(3, k) = solver%mat(3, k) + solver%bc_coeffs_rmin(3)*solver%mat(1, k)
            solver%mat(2, k) = solver%mat(2, k) + solver%bc_coeffs_rmin(2)*solver%mat(1, k)
            solver%mat(1, k) = 0.0_f64

         else if (bck(1) == sll_p_polar_origin) then ! center of circular domain

            ! Gaussian elimination: phi(1) = (-1)^k phi(2)
            solver%mat(2, k) = solver%mat(2, k) + (-1)**k*solver%mat(1, k)
            solver%mat(1, k) = 0.0_f64

         end if

         !--------------------------------------------
         ! Set boundary condition at rmax
         !--------------------------------------------
         last = 3*(nr - 1)
         if (bck(2) == sll_p_dirichlet) then ! Dirichlet
            solver%mat(last, k) = 0.0_f64

         else if (bck(2) == sll_p_neumann) then ! Neumann

            ! Coefficients of homogeneous boundary condition
            hp = r_nodes(nr + 1) - r_nodes(nr)
            hm = r_nodes(nr) - r_nodes(nr - 1)
            d1_coeffs(:) = [hp/hm, -2 - hp/hm - hm/hp, 2 + hm/hp]
            solver%bc_coeffs_rmax(-2:-1) = -d1_coeffs(1:2)/d1_coeffs(3)

            ! Gaussian elimination: remove phi(last) variable using boundary condition
            solver%mat(last - 2, k) = solver%mat(last - 2, k) + solver%bc_coeffs_rmax(-2)*solver%mat(last, k)
            solver%mat(last - 1, k) = solver%mat(last - 1, k) + solver%bc_coeffs_rmax(-1)*solver%mat(last, k)
            solver%mat(last, k) = 0.0_f64

         end if

      end do

   end subroutine sll_s_qn_solver_2d_polar_init

   !=============================================================================
   !> Solve the quasi-neutrality equation and get the electrostatic potential
   subroutine sll_s_qn_solver_2d_polar_solve(solver, rho, phi)
      type(sll_t_qn_solver_2d_polar), intent(inout) :: solver   !< Solver object
      real(f64), intent(in) :: rho(:, :) !< Charge density
      real(f64), intent(out) :: phi(:, :) !< Potential

      integer(i32) :: nr, ntheta, bck(2)
      integer(i32) :: i, k
      integer(i32) :: nrpts
      integer(i32) :: sh

      nr = solver%nr
      ntheta = solver%nt

      ! Shift in radial grid indexing
      sh = solver%skip0 ! =1 if bc_rmin==sll_p_polar_origin, =0 otherwise

      ! Number of points in radial grid
      nrpts = nr + 1 - sh

      ! Consistency check: 'rho' and 'phi' have shape defined at initialization
      SLL_ASSERT_ALWAYS(all(shape(rho) == [nrpts, ntheta]))
      SLL_ASSERT_ALWAYS(all(shape(phi) == [nrpts, ntheta]))

      ! For each r_i, compute FFT of rho(r_i,theta) to obtain \hat{rho}(r_i,k)
      do i = 1, nrpts
         solver%temp_r(:) = rho(i, :)
         call sll_s_fft_exec_r2c_1d(solver%fw, solver%temp_r(:), solver%temp_c(:))
         solver%z(i, :) = solver%temp_c(:)
      end do

      ! Cycle over k
      do k = 0, ntheta/2

         ! rhok(r) is k-th Fourier mode of rho(r,theta)
         ! phik(r) is k-th Fourier mode of phi(r,theta)
         ! rhok is 1D contiguous slice (column) of solver%z
         ! we will overwrite rhok with phik
         associate (rhok => solver%z(:, k), phik => solver%z(:, k))

            rhok(:) = rhok(:)/(solver%g(1 + sh:)*solver%epsilon_0)

            ! Solve tridiagonal system to obtain \hat{phi}_{k_j}(r) at internal points
            call sll_s_setup_cyclic_tridiag(solver%mat(:, k), nr - 1, solver%cts, solver%ipiv)
            call sll_o_solve_cyclic_tridiag(solver%cts, solver%ipiv, &
                                            rhok(2 - sh:nrpts - 1), nr - 1, phik(2 - sh:nrpts - 1))

            ! Compute boundary conditions type for mode k
            bck(:) = solver%bc(:)
            do i = 1, 2
               if (bck(i) == sll_p_neumann_mode_0) then
                  if (k == 0) then
                     bck(i) = sll_p_neumann
                  else
                     bck(i) = sll_p_dirichlet
                  end if
               end if
            end do

            ! Boundary condition at rmin
            if (bck(1) == sll_p_dirichlet) then ! Dirichlet
               phik(1) = (0.0_f64, 0.0_f64)
            else if (bck(1) == sll_p_neumann) then ! Neumann
               associate (c => solver%bc_coeffs_rmin)
                  phik(1) = c(2)*phik(2) + c(3)*phik(3)
               end associate
            end if

            ! Boundary condition at rmax
            if (bck(2) == sll_p_dirichlet) then ! Dirichlet
               phik(nrpts) = (0.0_f64, 0.0_f64)
            else if (bck(2) == sll_p_neumann) then ! Neumann
               associate (c => solver%bc_coeffs_rmax)
                  phik(nrpts) = c(-2)*phik(nrpts - 2) + c(-1)*phik(nrpts - 1)
               end associate
            end if

         end associate

      end do

      ! For each r_i, compute inverse FFT of \hat{phi}(r_i,k) to obtain phi(r_i,theta)
      do i = 1, nrpts
         solver%temp_c(:) = solver%z(i, :)
         call sll_s_fft_exec_c2r_1d(solver%bw, solver%temp_c(:), solver%temp_r(:))
         phi(i, :) = solver%temp_r(:)
      end do

   end subroutine sll_s_qn_solver_2d_polar_solve

   !=============================================================================
   !> Delete contents (local storage) of quasi-neutrality solver
   subroutine sll_s_qn_solver_2d_polar_free(solver)
      type(sll_t_qn_solver_2d_polar), intent(inout) :: solver

      call sll_s_fft_free(solver%fw)
      call sll_s_fft_free(solver%bw)

      call sll_s_fft_deallocate_aligned_real(solver%temp_r)
      call sll_s_fft_deallocate_aligned_complex(solver%temp_c)

      deallocate (solver%z)
      deallocate (solver%mat)
      deallocate (solver%cts)
      deallocate (solver%ipiv)

   end subroutine sll_s_qn_solver_2d_polar_free

end module sll_m_qn_solver_2d_polar
