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

!> @ingroup quasi_neutral_solvers_parallel
!>
!> @brief
!> Parallel quasi-neutrality solver on 2D polar mesh; uses FFTs along theta.
!>
!> @authors Yaman Güçlü, IPP Garching
!> @authors Edoardo Zoni, IPP Garching
!>
!> @details
!> This module solves the quasi-neutrality equation on a 2D polar mesh
!> \f$ (r,\theta) \in [r_\textrm{min},r_\textrm{max}]\times[0,2\pi) \f$ 
!> using fast Fourier transforms (FFTs) in \f$ \theta \f$ and 2nd-order
!> finite-difference methods in \f$ r \f$.
!>
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
!> - Homogeneous Dirichlet: \f$ \phi(\overline{r},\theta)=0\f$;
!> - Homogeneous Neumann: \f$ \partial_r\phi(\overline{r},\theta)=0\f$;
!> - Homogeneous Neumann mode 0:
!>   \f$ \partial_r\widehat{\phi}_0(\overline{r})=0\f$ and
!>   \f$ \widehat{\phi}_k(\overline{r})=0 \f$ for \f$ k\neq 0 \f$.
!>
!> The following arguments are given by the user at initialization:
!> \f$ \rho_{m0} \f$, \f$ B \f$, \f$ \lambda_D \f$, \f$ \epsilon_0 \f$ and the
!> additional parameter \f$ \chi \f$ (default is \f$ \chi=1 \f$).
!> If \f$ \lambda_D \f$ is not given to the solver, we assume
!> \f$ 1/\lambda_D^2=0 \f$: the electrons form a kinetic species and their
!> contribution goes into \f$ \rho_{c1} \f$. If \f$ \epsilon_0 \f$ is not given
!> to the solver, we assume \f$ \epsilon_0=1 \f$.
!>
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
!> ##### Note on parallelization #####
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

module sll_m_qn_solver_2d_polar_par
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"
#include "sll_working_precision.h"

  use sll_m_constants, only: &
    sll_p_epsilon_0

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_dirichlet, &
    sll_p_neumann, &
    sll_p_neumann_mode_0

  use sll_m_collective, only: &
    sll_f_get_collective_size, &
    sll_s_halt_collective, &
    sll_v_world_collective

  use sll_m_fft, only: &
    sll_t_fft, &
    sll_p_fft_forward, &
    sll_p_fft_backward, &
    sll_f_fft_allocate_aligned_real, &
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

  implicit none

  public :: &
    sll_t_qn_solver_2d_polar_par, &
    sll_s_qn_solver_2d_polar_par_init, &
    sll_s_qn_solver_2d_polar_par_solve, &
    sll_s_qn_solver_2d_polar_par_free

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Class for the Poisson solver in polar coordinate
  type sll_t_qn_solver_2d_polar_par

   sll_real64              :: rmin      !< Min value of r coordinate
   sll_real64              :: rmax      !< Max value of r coordinate
   sll_int32               :: nr        !< Number of cells along r
   sll_int32               :: ntheta    !< Number of cells along theta
   sll_int32               :: bc(2)     !< Boundary conditions options
   sll_real64              :: epsilon_0 !< Vacuum permittivity (user may override)
   sll_real64, allocatable :: g(:)      !< \f$ g(r) = \rho_{m,0}/(B^2\epsilon_0) \f$

   type(sll_t_fft)         :: fw        !< Forward FFT plan
   type(sll_t_fft)         :: bw        !< Inverse FFT plan
   sll_real64, pointer     :: tmp (:)   !< 1D work array for FFT, real
   sll_int32 , allocatable :: k_list(:)
   sll_real64, allocatable :: z_r (:,:) !< 2D array sequential in r
   sll_real64, pointer     :: z_a (:,:) !< 2D array sequential in theta
   sll_real64, allocatable :: mat (:,:) !< Tridiagonal matrix (one for each k)
   sll_real64, allocatable :: cts (:)   !< Lapack coefficients
   sll_int32 , allocatable :: ipiv(:)   !< Lapack pivot indices

   type(sll_t_layout_2d)           , pointer :: layout_r !< layout sequential in r
   type(sll_t_layout_2d)           , pointer :: layout_a !< layout sequential in theta
   type(sll_t_remap_plan_2d_real64), pointer :: rmp_ra   !< remap r->theta 
   type(sll_t_remap_plan_2d_real64), pointer :: rmp_ar   !< remap theta->r

  end type sll_t_qn_solver_2d_polar_par

  ! Allowed boundary conditions
  sll_int32, parameter :: bc_opts(3) = &
    [sll_p_dirichlet, sll_p_neumann, sll_p_neumann_mode_0]

  ! Local parameters
  sll_real64, parameter ::  one_third  = 1.0_f64 / 3.0_f64
  sll_real64, parameter :: four_thirds = 4.0_f64 / 3.0_f64

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !=============================================================================
  !> Initialize the solver
  subroutine sll_s_qn_solver_2d_polar_par_init( solver, &
      layout_r      , &
      layout_a      , &
      rmin          , &
      rmax          , &
      nr            , &
      ntheta        , &
      bc_rmin       , &
      bc_rmax       , &
      rho_m0        , &
      b_magn        , &
      lambda        , &
      use_zonal_flow, &
      epsilon_0 )

    type(sll_t_qn_solver_2d_polar_par), intent(out) :: solver !< solver object
    type(sll_t_layout_2d), pointer    :: layout_r       !< layout sequential in r direction
    type(sll_t_layout_2d), pointer    :: layout_a       !< layout sequential in theta direction
    sll_real64           , intent(in) :: rmin           !< rmin
    sll_real64           , intent(in) :: rmax           !< rmax
    sll_int32            , intent(in) :: nr             !< number of cells radial
    sll_int32            , intent(in) :: ntheta         !< number of cells angular
    sll_int32            , intent(in) :: bc_rmin        !< boundary condition at r_min
    sll_int32            , intent(in) :: bc_rmax        !< boundary condition at r_max
    sll_real64           , intent(in) :: rho_m0(:)      !< radial profile: total mass density of equilibrium
    sll_real64           , intent(in) :: b_magn(:)      !< radial profile: intensity of magnetic field
    sll_real64,  optional, intent(in) :: lambda(:)      !< radial profile: electron Debye length
    logical   ,  optional, intent(in) :: use_zonal_flow !< if .false. set flux average to zero
    sll_real64,  optional, intent(in) :: epsilon_0      !< override default: vacuum permittivity

    character(len=*), parameter :: this_sub_name = 'sll_s_qn_solver_2d_polar_par_init'

    sll_real64 :: dr
    sll_real64 :: inv_r
    sll_real64 :: inv_dr
    sll_real64 :: c
    sll_real64 :: d1(-1:+1)
    sll_real64 :: d2(-1:+1)

    sll_int32  :: i, j, k
    sll_int32  :: bck(2)
    sll_int32  :: last

    sll_int32  :: loc_sz_r(2) ! local shape of layout_r
    sll_int32  :: loc_sz_a(2) ! local shape of layout_a
    sll_int32  :: glob_idx(2) ! global indices

    sll_int32, allocatable :: k_list_glob(:)

    ! Consistency check: boundary conditions must be one of three options
    if( any( bc_rmin == bc_opts ) ) then
      solver%bc(1) = bc_rmin
    else
      SLL_ERROR( this_sub_name, 'Unrecognized boundary condition at r_min' )
    end if
    !
    if( any( bc_rmax == bc_opts ) ) then
      solver%bc(2) = bc_rmax
    else
      SLL_ERROR( this_sub_name, 'Unrecognized boundary condition at r_max' )
    end if

    ! Consistency check: global size of 2D layouts must be (nr+1,ntheta)
    SLL_ASSERT( sll_o_get_layout_global_size_i( layout_r ) == nr+1   )
    SLL_ASSERT( sll_o_get_layout_global_size_j( layout_r ) == ntheta )
    !
    SLL_ASSERT( sll_o_get_layout_global_size_i( layout_a ) == nr+1   )
    SLL_ASSERT( sll_o_get_layout_global_size_j( layout_a ) == ntheta )

    ! Compute local size of 2D arrays in the two layouts
    call sll_o_compute_local_sizes( layout_r, loc_sz_r(1), loc_sz_r(2) )
    call sll_o_compute_local_sizes( layout_a, loc_sz_a(1), loc_sz_a(2) )

    ! Consistency check: layout_r sequential in r, layout_a sequential in theta
    SLL_ASSERT( loc_sz_r(1) == nr+1   )
    SLL_ASSERT( loc_sz_a(2) == ntheta )

    ! Override vacuum permittivity in SI units
    if (present( epsilon_0 )) then
      solver%epsilon_0 = epsilon_0
    else
      solver%epsilon_0 = sll_p_epsilon_0
    end if

    ! Store global information in solver
    solver%rmin     =  rmin
    solver%rmax     =  rmax
    solver%nr       =  nr
    solver%ntheta   =  ntheta
    solver%layout_a => layout_a
    solver%layout_r => layout_r

    ! Allocate arrays global in r
    allocate( solver%g   (nr+1) )
    allocate( solver%z_r (nr+1,   loc_sz_r(2)) )
    allocate( solver%mat((nr-1)*3,loc_sz_r(2)) ) ! for each k, matrix depends on r
    allocate( solver%cts((nr-1)*7) )
    allocate( solver%ipiv(nr-1) )

    ! Remap objects between two layouts (for transposing between f_r and f_a)
    allocate( solver%z_a( loc_sz_a(1), ntheta ) )
    solver%rmp_ra => sll_o_new_remap_plan( solver%layout_r, solver%layout_a, solver%z_r )
    solver%rmp_ar => sll_o_new_remap_plan( solver%layout_a, solver%layout_r, solver%z_a )
    deallocate( solver%z_a )

    ! Allocate in ALIGNED fashion 1D array for storing one row of f_a
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

    ! Store non-dimensional coefficient g(r) = \rho(r) / (B(r)^2 \epsilon_0)
    solver%g(:) = rho_m0(:) / (b_magn(:)**2 * solver%epsilon_0)

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

    ! Precompute convenient parameters
    dr     = (rmax-rmin)/nr
    inv_dr = 1.0_f64 / dr

    ! Store matrix coefficients into solver%mat
    ! Cycle over k_j
    do j = 1, loc_sz_r(2)

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

      ! Compute matrix coefficients for a given k_j
      do i = 2, nr

        associate( g => solver%g(:) )

        if (present( lambda )) then
          c = 1.0_f64 / (lambda(i)**2 * g(i))
          if (present( use_zonal_flow )) then
            if (use_zonal_flow .and. k == 0) then
              c = 0.0_f64
            end if
          end if
        else
          c = 0.0_f64
        end if

        inv_r    = 1.0_f64 / (rmin + (i-1)*dr)
        d1(-1:1) = [-0.5_f64, 0.0_f64, 0.5_f64] * inv_dr
        d2(-1:1) = [0.5_f64*(g(i-1)+g(i))/g(i), -2.0_f64, 0.5_f64*(g(i)+g(i+1))/g(i)] * inv_dr**2

        end associate

        solver%mat(3*(i-1)  ,j) = -d2( 1) -d1( 1)*inv_r
        solver%mat(3*(i-1)-1,j) = -d2( 0) -d1( 0)*inv_r  + (k*inv_r)**2 + c
        solver%mat(3*(i-1)-2,j) = -d2(-1) -d1(-1)*inv_r
      end do

      ! Set boundary condition at rmin
      if (bck(1) == sll_p_dirichlet) then ! Dirichlet
        solver%mat(1,j) = 0.0_f64
      else if (bck(1) == sll_p_neumann) then ! Neumann
        solver%mat(3,j) = solver%mat(3,j) -  one_third  * solver%mat(1,j)
        solver%mat(2,j) = solver%mat(2,j) + four_thirds * solver%mat(1,j)
        solver%mat(1,j) = 0.0_f64
      end if

      ! Set boundary condition at rmax
      last = 3*(nr-1)
      if (bck(2) == sll_p_dirichlet) then ! Dirichlet
        solver%mat(last,j) = 0.0_f64
      else if (bck(2) == sll_p_neumann) then ! Neumann
        solver%mat(last-2,j) = solver%mat(last-2,j) -  one_third  *solver%mat(last,j)
        solver%mat(last-1,j) = solver%mat(last-1,j) + four_thirds *solver%mat(last,j)
        solver%mat(last  ,j) = 0.0_f64
      end if

    end do

  end subroutine sll_s_qn_solver_2d_polar_par_init


  !=============================================================================
  !> Solve the quasi-neutrality equation and get the electrostatic potential
  subroutine sll_s_qn_solver_2d_polar_par_solve( solver, rho, phi )
    type(sll_t_qn_solver_2d_polar_par) , intent(inout) :: solver   !< Solver object
    sll_real64                         , intent(in   ) :: rho(:,:) !< Charge density
    sll_real64, target                 , intent(  out) :: phi(:,:) !< Potential

    sll_int32  :: nr, ntheta, bck(2)
    sll_int32  :: i, j, k

    nr     = solver%nr
    ntheta = solver%ntheta

    ! Consistency check: rho and phi must be given in layout sequential in theta
    call verify_argument_sizes_par( solver%layout_a, rho )
    call verify_argument_sizes_par( solver%layout_a, phi )

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
      associate( rhok => solver%z_r(:,j)/(solver%g(:)*solver%epsilon_0), &
                 phik => solver%z_r(:,j) )

        ! Solve tridiagonal system to obtain \hat{phi}_{k_j}(r) at internal points
        call sll_s_setup_cyclic_tridiag( solver%mat(:,j), nr-1, solver%cts, solver%ipiv )
        call sll_o_solve_cyclic_tridiag( solver%cts, solver%ipiv, rhok(2:nr), &
             nr-1, phik(2:nr) )

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
          phik(1) = four_thirds*phik(2) - one_third*phik(3)
        end if

        ! Boundary condition at rmax
        if (bck(2) == sll_p_dirichlet) then ! Dirichlet
          phik(nr+1) = 0.0_f64
        else if (bck(2) == sll_p_neumann) then ! Neumann
          phik(nr+1) = four_thirds*phik(nr) - one_third*phik(nr-1)
        end if

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

  end subroutine sll_s_qn_solver_2d_polar_par_solve


  !=============================================================================
  !> Delete contents (local storage) of quasi-neutrality solver
  subroutine sll_s_qn_solver_2d_polar_par_free( solver )
    type(sll_t_qn_solver_2d_polar_par), intent(inout) :: solver

    call sll_s_fft_free( solver%fw )
    call sll_s_fft_free( solver%bw )
    deallocate( solver%tmp  )

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

  end subroutine sll_s_qn_solver_2d_polar_par_free


  !=============================================================================
  !> Check if array sizes are compatble with the layout 
  subroutine verify_argument_sizes_par(layout, array)

    type(sll_t_layout_2d), pointer       :: layout
    sll_real64, dimension(:,:)     :: array
    sll_int32,  dimension(2)       :: n ! nx_loc, ny_loc
    sll_int32                      :: i

    ! Note that this checks for strict sizes, not an array being bigger
    ! than a certain size, but exactly a desired size... This may be a bit
    ! too stringent.
    call sll_o_compute_local_sizes( layout, n(1), n(2) )

    do i=1,2
       if ( (n(i)/=size(array,i))) then
          print*, 'ERROR: solve_poisson_polar_parallel()', &
               'size of either rho or phi does not match expected size. '
          if (i==1) then
             print*, 'solve_poisson_polar_parallel(): ', &
                  'mismatch in direction r'
          else if (i==2) then
             print*, 'solve_poisson_polar_parallel(): ', &
                  'mismatch in direction theta'
          endif
          print *, 'Exiting...'
          call sll_s_halt_collective()
          stop
       endif
    enddo

  end subroutine verify_argument_sizes_par

end module sll_m_qn_solver_2d_polar_par
