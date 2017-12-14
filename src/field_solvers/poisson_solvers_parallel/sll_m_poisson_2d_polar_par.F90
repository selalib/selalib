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

!> @ingroup poisson_solvers_parallel
!> @brief
!> Parallel Poisson solver on 2D polar mesh; uses FFT in theta and 2nd-order FD in r.
!>
!> @authors Yaman Güçlü, IPP Garching
!> @authors Edoardo Zoni, IPP Garching
!>
!> @details
!>
!> #### Model equations ####
!>
!> This module solves the Poisson equation \f$ -\nabla^2 \phi = \rho \f$
!> (permittivity of free space \f$ \epsilon_0 = 1 \f$) on a 2D polar mesh
!> \f$ (r,\theta) \in [r_\textrm{min},r_\textrm{max}]\times[0,2\pi) \f$.
!> The Poisson equation in polar coordinates reads
!> \f[
!> -\bigg(
!> \frac{\partial^2}{\partial r^{2}}
!> +\frac{1}{r}\frac{\partial}{\partial r}
!> +\frac{1}{r^{2}}\frac{\partial^{2}}{\partial\theta^{2}}
!> \bigg)\phi(r,\theta) = \rho(r,\theta).
!> \f]
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
!> -\bigg(
!> \frac{\partial^2}{\partial r^{2}}
!> +\frac{1}{r}\frac{\partial}{\partial r}
!> -\frac{k^{2}}{r^{2}}
!> \bigg)\widehat{\phi}_k(r) = \widehat{\rho}_k(r).
!> \f]
!> For each mode \f$ k \f$, the resulting ODE is solved with a 2nd-order
!> finite-difference collocation method.
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
!> Poisson equation is differentiated to obtain an expression for
!> \f$ \partial^3\widehat{\phi}_k / \partial r^3 \f$ that involves only
!> 1st and 2nd derivatives,
!>
!> \f[
!>                     \frac{\partial^3\widehat\phi_k}{\partial r^3} =
!>  -\frac{1}{r}       \frac{\partial^2\widehat\phi_k}{\partial r^2}
!>  +\frac{1+k^2}{r^2} \frac{\partial  \widehat\phi_k}{\partial r  }
!>  -\frac{2 k^2}{r^3}      {          \widehat\phi_k}
!>  -                  \frac{\partial  \widehat\rho_k}{\partial r  },
!> \f]
!>
!> which are then approximated with the same finite-difference formulas.
!>
!> #### Parallelization ####
!>
!> - The user must provide two compatible 2D layouts:
!>   - \c layout_a sequential in \f$ \theta \f$;
!>   - \c layout_r sequential in \f$ r \f$;
!> - Input and output data are given in \c layout_a;
!> - \c layout_r is only used internally;
!> - FFTs are computed on one row of \c layout_a;
!> - The linear solver (for each mode \f$ k \f$) works on one column of \c layout_r;
!> - The maximum number of parallel solve is \f$ N_\theta \f$,
!>   where \f$ N_\theta\f$ is the number of cells along \f$ \theta \f$.
!>   Therefore, instead of working with \f$ N_\theta/2+1 \f$ Fourier modes, we
!>   work with \f$ N_\theta \f$ real modes. For this purpose, the \c r2r interface
!>   of FFTW is used: this corresponds to DFTs of real input and complex-Hermitian
!>   output in halfcomplex format.
!>
!> #### Usage example ####
!>
!> \code
!> use sll_m_poisson_2d_polar_par, only: &
!>   sll_t_poisson_2d_polar_par,         &
!>   sll_s_poisson_2d_polar_par_init,    &
!>   sll_s_poisson_2d_polar_par_solve,   &
!>   sll_s_poisson_2d_polar_par_free
!>
!> use sll_m_remapper, only: &
!>   sll_t_layout_2d
!>
!> ...
!>
!> type(sll_t_poisson_2d_polar_par) :: solver
!> type(sll_t_layout_2d)            :: layout_a, layout_r
!> real(f64), allocatable           :: rho(:,:), phi(:,:)
!> real(f64)                        :: rmin, rmax
!> integer(i32)                     :: nr, ntheta, bc_rmin, bc_rmax
!>
!> ...
!>
!> call sll_s_poisson_2d_polar_par_init ( solver, layout_r, layout_a, rmin, rmax, nr, ntheta, bc_rmin, bc_rmax )
!> call sll_s_poisson_2d_polar_par_solve( solver, rho, phi )
!> call sll_s_poisson_2d_polar_par_free ( solver )
!> \endcode

module sll_m_poisson_2d_polar_par
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"
#include "sll_working_precision.h"

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_dirichlet, &
    sll_p_neumann, &
    sll_p_neumann_mode_0, &
    sll_p_polar_origin

  use sll_m_collective, only: &
    sll_f_get_collective_size, &
    sll_s_halt_collective, &
    sll_v_world_collective

  use sll_m_fft, only: &
    sll_t_fft, &
    sll_p_fft_forward, &
    sll_p_fft_backward, &
    sll_f_fft_allocate_aligned_real, &
    sll_s_fft_deallocate_aligned_real, &
    sll_s_fft_init_r2r_1d, &
    sll_s_fft_exec_r2r_1d, &
    sll_s_fft_get_k_list_r2r_1d, &
    sll_s_fft_free

  use sll_m_remapper, only: &
    sll_t_layout_2d, &
    sll_o_compute_local_sizes, &
    sll_o_get_layout_global_size_i, &
    sll_o_get_layout_global_size_j, &
    sll_o_local_to_global, &
    sll_t_remap_plan_2d_real64, &
    sll_o_new_remap_plan, &
    sll_o_apply_remap_2d

  use sll_m_tridiagonal, only: &
    sll_s_setup_cyclic_tridiag, &
    sll_o_solve_cyclic_tridiag

  use sll_m_utilities, only: &
    sll_s_new_array_linspace

  implicit none

  public :: &
    sll_t_poisson_2d_polar_par, &
    sll_s_poisson_2d_polar_par_init, &
    sll_s_poisson_2d_polar_par_solve, &
    sll_s_poisson_2d_polar_par_free

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Class for the Poisson solver in polar coordinate
  type sll_t_poisson_2d_polar_par

   real   (f64)              :: rmin       !< Min value of r coordinate
   real   (f64)              :: rmax       !< Max value of r coordinate
   integer(i32)              :: nr         !< Number of cells along r
   integer(i32)              :: ntheta     !< Number of cells along theta
   integer(i32)              :: bc(2)      !< Boundary conditions options

   type(sll_t_fft)           :: fw         !< Forward FFT plan
   type(sll_t_fft)           :: bw         !< Inverse FFT plan
   real   (f64), pointer     :: tmp (:)    !< 1D work array for FFT, real
   integer(i32), allocatable :: k_list(:)
   real   (f64), allocatable :: z_r (:,:)  !< 2D array sequential in r
   real   (f64), pointer     :: z_a (:,:)  !< 2D array sequential in theta
   real   (f64), allocatable :: mat (:,:)  !< Tridiagonal matrix (one for each k)
   real   (f64), allocatable :: cts (:)    !< Lapack coefficients
   integer(i32), allocatable :: ipiv(:)    !< Lapack pivot indices

   real(f64) :: bc_coeffs_rmin( 2: 3) ! needed for Neumann at r=r_min
   real(f64) :: bc_coeffs_rmax(-2:-1) ! needed for Neumann at r=r_max
   integer   :: skip0                 ! needed for full circle

   type(sll_t_layout_2d)           , pointer :: layout_r !< layout sequential in r
   type(sll_t_layout_2d)           , pointer :: layout_a !< layout sequential in theta
   type(sll_t_remap_plan_2d_real64), pointer :: rmp_ra   !< remap r->theta 
   type(sll_t_remap_plan_2d_real64), pointer :: rmp_ar   !< remap theta->r

  end type sll_t_poisson_2d_polar_par

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !=============================================================================
  !> sll_o_initialize the Poisson solver in polar coordinates
  subroutine sll_s_poisson_2d_polar_par_init( solver, &
      layout_r, &
      layout_a, &
      rmin    , &
      rmax    , &
      nr      , &
      ntheta  , &
      bc_rmin , &
      bc_rmax , &
      rgrid   )

    type(sll_t_poisson_2d_polar_par), intent(  out) :: solver   !< solver object
    type(sll_t_layout_2d)           , pointer       :: layout_r !< layout sequential in r direction
    type(sll_t_layout_2d)           , pointer       :: layout_a !< layout sequential in theta direction
    real(f64)                       , intent(in   ) :: rmin     !< rmin
    real(f64)                       , intent(in   ) :: rmax     !< rmax
    integer(i32)                    , intent(in   ) :: nr       !< number of cells radial
    integer(i32)                    , intent(in   ) :: ntheta   !< number of cells angular
    integer(i32)                    , intent(in   ) :: bc_rmin  !< boundary condition at r_min
    integer(i32)                    , intent(in   ) :: bc_rmax  !< boundary condition at r_max
    real(f64),      target, optional, intent(in   ) :: rgrid(:) !< grid points along r

    character(len=*), parameter :: this_sub_name = 'sll_s_poisson_2d_polar_par_init'

    real(f64) :: hp, hm
    real(f64) :: inv_r
    real(f64) :: d1_coeffs(3)
    real(f64) :: d2_coeffs(3)
    real(f64), allocatable :: r_nodes(:)

    integer(i32)  :: i, j, k
    integer(i32)  :: bck(2)
    integer(i32)  :: last
    integer(i32)  :: sh

    integer(i32)  :: loc_sz_r(2) ! sequential in r direction
    integer(i32)  :: loc_sz_a(2) ! sequential in theta direction
    integer(i32)  :: glob_idx(2)

    integer(i32), allocatable :: k_list_glob(:)

    ! Consistency checks
    SLL_ASSERT_ALWAYS( rmin >= 0.0_f64 )
    SLL_ASSERT_ALWAYS( rmin < rmax )
    SLL_ASSERT_ALWAYS( nr >= 1 )
    SLL_ASSERT_ALWAYS( ntheta >= 1 )

    if (bc_rmin == sll_p_polar_origin .and. rmin /= 0.0_f64) then
      SLL_ERROR( this_sub_name, "BC option 'sll_p_polar_origin' requires r_min = 0" )
    end if

    ! Set boundary condition at r_min
    select case( bc_rmin )
    case( sll_p_dirichlet, sll_p_neumann, sll_p_neumann_mode_0, sll_p_polar_origin )
      solver%bc(1) = bc_rmin
    case default
      SLL_ERROR( this_sub_name, 'Unrecognized boundary condition at r_min' )
    end select

    ! Set boundary condition at r_max
    select case( bc_rmax )
    case( sll_p_dirichlet, sll_p_neumann, sll_p_neumann_mode_0 )
      solver%bc(2) = bc_rmax
    case default
      SLL_ERROR( this_sub_name, 'Unrecognized boundary condition at r_max' )
    end select

    ! Important: if full circle is simulated, center point is not solved for!
    !            Therefore, solution is calculated on nr+1-skip0 points.
    sh = merge( 1, 0, bc_rmin==sll_p_polar_origin )

    ! Consistency check: global size of 2D layouts must be (nr+1,ntheta)
    SLL_ASSERT_ALWAYS( sll_o_get_layout_global_size_i( layout_r ) == nr+1-sh )
    SLL_ASSERT_ALWAYS( sll_o_get_layout_global_size_j( layout_r ) == ntheta  )
    !
    SLL_ASSERT_ALWAYS( sll_o_get_layout_global_size_i( layout_a ) == nr+1-sh )
    SLL_ASSERT_ALWAYS( sll_o_get_layout_global_size_j( layout_a ) == ntheta  )

    ! Compute local size of 2D arrays in the two layouts
    call sll_o_compute_local_sizes( layout_r, loc_sz_r(1), loc_sz_r(2) )
    call sll_o_compute_local_sizes( layout_a, loc_sz_a(1), loc_sz_a(2) )

    ! Consistency check: layout_r sequential in r, layout_a sequential in theta
    SLL_ASSERT_ALWAYS( loc_sz_r(1) == nr+1-sh )
    SLL_ASSERT_ALWAYS( loc_sz_a(2) == ntheta  )

    ! Store global information in solver
    solver%rmin     =  rmin
    solver%rmax     =  rmax
    solver%nr       =  nr
    solver%ntheta   =  ntheta
    solver%layout_a => layout_a
    solver%layout_r => layout_r
    solver%skip0    =  sh

    ! r grid (possibly non-uniform)
    allocate( r_nodes(nr+1) )

    if (present( rgrid )) then  !--> Create computational grid from user data

      SLL_ASSERT_ALWAYS( all( rgrid > 0.0_f64 ) )
      if (bc_rmin == sll_p_polar_origin) then
        SLL_ASSERT_ALWAYS( size(rgrid) == nr   )
        SLL_ASSERT_ALWAYS( rgrid(nr  ) == rmax )
        r_nodes(1 ) = -rgrid(1)
        r_nodes(2:) =  rgrid(:)
      else
        SLL_ASSERT_ALWAYS( size(rgrid) == nr+1 )
        SLL_ASSERT_ALWAYS( rgrid(   1) == rmin )
        SLL_ASSERT_ALWAYS( rgrid(nr+1) == rmax )
        r_nodes(:) = rgrid(:)
      end if

    else  !-------------------------> Create uniform grid

      if (bc_rmin == sll_p_polar_origin) then
        associate( rmin => rmax / real(2*nr+1,f64) )
          r_nodes(1) = -rmin
          call sll_s_new_array_linspace( r_nodes(2:), rmin, rmax, endpoint=.true. )
        end associate
      else
        call sll_s_new_array_linspace( r_nodes, rmin, rmax, endpoint=.true. )
      end if

    end if

    ! Allocate arrays global in r
    allocate( solver%z_r (nr+1+sh,loc_sz_r(2)) )
    allocate( solver%mat((nr-1)*3,loc_sz_r(2)) ) ! for each k, matrix depends on r
    allocate( solver%cts((nr-1)*7) )
    allocate( solver%ipiv(nr-1) )

    ! Remap objects between two layouts (for transposing between f_r and f_a)
    allocate( solver%z_a( loc_sz_a(1), ntheta ) )
    solver%rmp_ra => sll_o_new_remap_plan( solver%layout_r, solver%layout_a, solver%z_r )
    solver%rmp_ar => sll_o_new_remap_plan( solver%layout_a, solver%layout_r, solver%z_a )
    deallocate( solver%z_a )

    ! Allocate in ALIGNED fashion 1D work array for FFT
    solver%tmp => sll_f_fft_allocate_aligned_real( ntheta )

    ! Initialize plans for forward and backward FFTs
    call sll_s_fft_init_r2r_1d( solver%fw, &
      ntheta             , &
      solver%tmp(:)      , &
      solver%tmp(:)      , &
      sll_p_fft_forward  , &
      aligned    = .true., &
      normalized = .true. )

    call sll_s_fft_init_r2r_1d( solver%bw, &
      ntheta             , &
      solver%tmp(:)      , &
      solver%tmp(:)      , &
      sll_p_fft_backward , &
      aligned    = .true., &
      normalized = .false. )

    ! Determine global k_list
    allocate( k_list_glob( ntheta ) )
    call sll_s_fft_get_k_list_r2r_1d( solver%fw, k_list_glob )
    ! Extract local k list
    allocate( solver%k_list(loc_sz_r(2)) )
    do j = 1, loc_sz_r(2)
      glob_idx(:) = sll_o_local_to_global( solver%layout_r, [1,j] )
      solver%k_list(j) = k_list_glob(glob_idx(2))
    end do
    deallocate( k_list_glob )

    ! Store matrix coefficients into solver%mat
    ! Cycle over k_j
    do j = 1, loc_sz_r(2)

      ! Get value of k_j from precomputed list of local values
      k = solver%k_list(j)

      !----------------------------------------------
      ! Compute boundary conditions type for mode k_j
      !----------------------------------------------
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

      !----------------------------------------------
      ! Compute matrix coefficients for a given k_j
      !----------------------------------------------
      do i = 2, nr
        hp = r_nodes(i+1)-r_nodes(i)
        hm = r_nodes(i)  -r_nodes(i-1)
        inv_r = 1.0_f64 / r_nodes(i)
        d1_coeffs(:) = [-hp/hm, (hp**2-hm**2)/(hp*hm), hm/hp] / (hp+hm)
        d2_coeffs(:) = [2*hp/(hp+hm), -2.0_f64, 2*hm/(hp+hm)] / (hp*hm)

        solver%mat(3*(i-1)-2:3*(i-1), j) = &
          - d2_coeffs(:) - d1_coeffs(:)*inv_r + [0.0_f64, (k*inv_r)**2, 0.0_f64]
      end do

      !----------------------------------------------
      ! Set boundary condition at rmin
      !----------------------------------------------
      if (bck(1) == sll_p_dirichlet) then ! Dirichlet
        solver%mat(1,j) = 0.0_f64

      else if (bck(1) == sll_p_neumann) then ! Neumann

        ! Coefficients of homogeneous boundary condition
        hp = r_nodes(3)-r_nodes(2)
        hm = r_nodes(2)-r_nodes(1)
        d1_coeffs(:) = [-2-hp/hm, 2+hp/hm+hm/hp, -hm/hp]
        solver % bc_coeffs_rmin(2:3) = -d1_coeffs(2:3)/d1_coeffs(1)

        ! Gaussian elimination: remove phi(1) variable using boundary condition
        solver%mat(3,j) = solver%mat(3,j) + solver%bc_coeffs_rmin(3) * solver%mat(1,j)
        solver%mat(2,j) = solver%mat(2,j) + solver%bc_coeffs_rmin(2) * solver%mat(1,j)
        solver%mat(1,j) = 0.0_f64

      else if (bck(1) == sll_p_polar_origin) then ! center of circular domain

        ! Gaussian elimination: phi(1) = (-1)^k phi(2)
        solver%mat(2,j) = solver%mat(2,j) + (-1)**k * solver%mat(1,j)
        solver%mat(1,j) = 0.0_f64

      end if

      !----------------------------------------------
      ! Set boundary condition at rmax
      !----------------------------------------------
      last = 3*(nr-1)
      if (bck(2) == sll_p_dirichlet) then ! Dirichlet
        solver%mat(last,j) = 0.0_f64

      else if (bck(2) == sll_p_neumann) then ! Neumann

        ! Coefficients of homogeneous boundary condition
        hp = r_nodes(nr+1)-r_nodes(nr)
        hm = r_nodes(nr)-r_nodes(nr-1)
        d1_coeffs(:) = [hp/hm, -2-hp/hm-hm/hp, 2+hm/hp]
        solver % bc_coeffs_rmax(-2:-1) = -d1_coeffs(1:2)/d1_coeffs(3)

        ! Gaussian elimination: remove phi(last) variable using boundary condition
        solver%mat(last-2,j) = solver%mat(last-2,j) + solver%bc_coeffs_rmax(-2) * solver%mat(last,j)
        solver%mat(last-1,j) = solver%mat(last-1,j) + solver%bc_coeffs_rmax(-1) * solver%mat(last,j)
        solver%mat(last  ,j) = 0.0_f64

      end if

    end do

  end subroutine sll_s_poisson_2d_polar_par_init


  !=============================================================================
  !> Solve the Poisson equation and get the electrostatic potential
  subroutine sll_s_poisson_2d_polar_par_solve( solver, rho, phi )
    type(sll_t_poisson_2d_polar_par) , intent(inout) :: solver   !< Solver object
    real(f64)                        , intent(in   ) :: rho(:,:) !< Charge density
    real(f64), target                , intent(  out) :: phi(:,:) !< Potential

    integer(i32) :: nr, ntheta, bck(2)
    integer(i32) :: i, j, k
    integer(i32) :: nrpts
    integer(i32) :: sh

    nr     = solver%nr
    ntheta = solver%ntheta

    ! Shift in radial grid indexing
    sh = solver%skip0 ! =1 if bc_rmin==sll_p_polar_origin, =0 otherwise

    ! Number of points in radial grid
    nrpts = nr+1-sh

    ! Consistency check: rho and phi must be given in layout sequential in theta
    call verify_argument_sizes_par( solver%layout_a, rho, 'rho' )
    call verify_argument_sizes_par( solver%layout_a, phi, 'phi' )

    ! Use output array 'phi' as 2D work array sequential in theta
    solver%z_a => phi(:,:)

    ! For each r_i, compute FFT of rho(r_i,theta) to obtain \hat{rho}(r_i,k)
    do i = 1, ubound( rho, 1 )
      solver%tmp(:) = rho(i,:)
      call sll_s_fft_exec_r2r_1d( solver%fw, solver%tmp(:), solver%tmp(:) )
      solver%z_a(i,:) = solver%tmp(:)
    end do

    ! Remap \hat{rho}(k) to layout distributed in k (global in r) -> \hat{rho}_k(r)
    call sll_o_apply_remap_2d( solver%rmp_ar, solver%z_a, solver%z_r )

    ! Cycle over k_j
    do j = 1, ubound( solver%z_r, 2 )

      ! rhok(r) is k-th Fourier mode of rho(r,theta)
      ! phik(r) is k-th Fourier mode of phi(r,theta)
      ! rhok is 1D contiguous slice (column) of solver%z
      ! we will overwrite rhok with phik
      associate( rhok => solver%z_r(:,j), phik => solver%z_r(:,j) )

        ! Solve tridiagonal system to obtain \hat{phi}_{k_j}(r) at internal points
        call sll_s_setup_cyclic_tridiag( solver%mat(:,j), nr-1, solver%cts, solver%ipiv )
        call sll_o_solve_cyclic_tridiag( solver%cts, solver%ipiv, rhok(2-sh:nrpts-1), nr-1, phik(2-sh:nrpts-1) )

        ! Get value of k_j from precomputed list of local values
        k = solver%k_list(j)

        ! Compute boundary conditions type for mode k_j
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
          phik(1) = 0.0_f64
        else if (bck(1) == sll_p_neumann) then ! Neumann
          associate( c => solver % bc_coeffs_rmin )
            phik(1) = c(2)*phik(2) + c(3)*phik(3)
          end associate
        end if

        ! Boundary condition at rmax
        if (bck(2) == sll_p_dirichlet) then ! Dirichlet
          phik(nrpts) = 0.0_f64
        else if (bck(2) == sll_p_neumann) then ! Neumann
          associate( c => solver % bc_coeffs_rmax )
            phik(nrpts) = c(-2)*phik(nrpts-2) + c(-1)*phik(nrpts-1)
          end associate
        endif

      end associate

    end do

    ! Redistribute \hat{phi}(r_i,k_j) into layout global in k
    call sll_o_apply_remap_2d( solver%rmp_ra, solver%z_r, solver%z_a )

    ! For each r_i, compute inverse FFT of \hat{phi}(r_i,k) to obtain phi(r_i,theta)
    do i = 1, ubound( solver%z_a, 1 )
      solver%tmp(:) = solver%z_a(i,:)
      call sll_s_fft_exec_r2r_1d( solver%bw, solver%tmp(:), solver%tmp(:) )
      phi(i,:) = solver%tmp(:)
    end do

  end subroutine sll_s_poisson_2d_polar_par_solve


  !=============================================================================
  !> Delete contents (local storage) of Poisson's solver
  subroutine sll_s_poisson_2d_polar_par_free( solver )
    type(sll_t_poisson_2d_polar_par) , intent(inout) :: solver

    call sll_s_fft_free( solver%fw )
    call sll_s_fft_free( solver%bw )

    call sll_s_fft_deallocate_aligned_real( solver%tmp )

    deallocate( solver%z_r  )
    deallocate( solver%mat  )
    deallocate( solver%cts  )
    deallocate( solver%ipiv )

    deallocate( solver%rmp_ra )
    deallocate( solver%rmp_ar )

    solver%layout_r => null()
    solver%layout_a => null()
    solver%rmp_ra   => null()
    solver%rmp_ar   => null()

  end subroutine sll_s_poisson_2d_polar_par_free


  !=============================================================================
  !> Check if array sizes are compatible with the layout 
  subroutine verify_argument_sizes_par( layout, array, array_name )
    type(sll_t_layout_2d), pointer    :: layout
    real(f64)            , intent(in) :: array(:,:)
    character(len=*)     , intent(in) :: array_name

    integer(i32)                :: n(2) ! nx_loc, ny_loc
    character(len=*), parameter :: sfmt = "('['(i0)','(i0)']')"
    character(len=256)          :: err_fmt
    character(len=256)          :: err_msg

    call sll_o_compute_local_sizes( layout, n(1), n(2) )

    if ( .not. all( shape(array)==n(:) ) ) then
      ! Create error format string
      err_fmt = "(a,"//sfmt//",a,"//sfmt//",a)"
      ! Write error message to string
      write( err_msg, err_fmt ) &
        "shape("//array_name//") = ", shape(array), &
        " is not compatible with local shape ", n ," of 2D layout"
      ! Stop execution with error
      SLL_ERROR( "sll_s_poisson_2d_polar_par_solve", trim(err_msg) )
    end if

  end subroutine verify_argument_sizes_par

end module sll_m_poisson_2d_polar_par
