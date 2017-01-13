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
!> @author  Yaman Güçlü, IPP Garching
!> @author  Edoardo Zoni, IPP Garching
!>
!> @details
!> Module to solve the quasi-neutrality equation on a 2D polar mesh using
!> FFT transforms in theta and 2nd order finite differences in r
!>
!> The general equation is of the form
!>
!> -\nabla_\perp \cdot ( \rho_{m,0}/(B^2 \epsilon_0) \nabla_\perp \phi )
!> + (\phi - \chi \langle\phi\rangle) / \lambda^2
!> = \rho_{c,1} / \epsilon_0,
!>
!> where
!>
!> \phi                 is electrostatic potential 
!> \langle\phi\rangle   is flux-average of \phi
!> \rho_{c,1}           is charge density perturbation due to kinetic species

!> The following terms are given by the user at initialization:
!>
!> \rho_{m,0}           is total mass density of equilibrium
!> B                    is intensity of equilibrium magnetic field 
!> \lambda              is electron Debye length of equilibrium
!> \epsilon_0           is vacuum permittivity (= 4\pi if CGS is used)
!>
!> Note:
!>
!> \chi = 1 by default, but it can be set to 0
!> If \lambda is not given to the solver, we assume 1/lambda=0: the electrons
!> are treated kinetically and therefore their contribution is part of \rho;
!>
!> In 2D polar geometry we solve
!>
!> -\partial_r ( g \partial_r \phi ) - \frac{g}{r} \partial_r \phi
!> -\frac{g}{r^2} \partial_\theta \phi + (\phi-\langle\phi\rangle) / \lambda^2
!> = \rho_{c,1} / \epsilon_0,
!>
!> where $g = \rho_{m,0}/(B^2 \epsilon_0)$ and we assume that
!> $g(r)$ and $\lambda(r)$ do not depend on theta.

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
    sll_p_neumann_mode_0

  use sll_m_fft, only: &
    sll_t_fft, &
    sll_p_fft_forward, &
    sll_p_fft_backward, &
    sll_f_fft_allocate_aligned_complex, &
    sll_s_fft_init_c2c_1d, &
    sll_s_fft_exec_c2c_1d, &
    sll_s_fft_free

  use sll_m_tridiagonal, only: &
    sll_s_setup_cyclic_tridiag, &
    sll_o_solve_cyclic_tridiag

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

   sll_real64              :: rmin     !< Min value of r coordinate
   sll_real64              :: rmax     !< Max value of r coordinate
   sll_int32               :: nr       !< Number of cells along r
   sll_int32               :: nt       !< Number of cells along theta
   sll_int32               :: bc(2)    !< Boundary conditions options
   sll_real64              :: epsilon_0!< Vacuum permittivity (user may override)
   sll_real64, allocatable :: g   (:)  !< g(r) = \rho_{m,0}/(B^2\epsilon_0)

   type(sll_t_fft)         :: fw       !< Forward FFT plan
   type(sll_t_fft)         :: bw       !< Inverse FFT plan
   sll_comp64, allocatable :: z   (:,:)!< 2D work array needed for transposition
   sll_comp64, pointer     :: zrow(:)  !< 1D slice (one row) of z, ALIGNED
   sll_comp64, allocatable :: fk  (:)  !< k-th Fourier mode of rho
   sll_comp64, allocatable :: phik(:)  !< k-th Fourier mode of phi
   sll_real64, allocatable :: mat (:,:)!< Tridiagonal matrix (one for each k)
   sll_real64, allocatable :: cts (:)  !< Lapack coefficients
   sll_int32 , allocatable :: ipiv(:)  !< Lapack pivot indices

  end type sll_t_qn_solver_2d_polar

  ! Allowed boundary conditions
  sll_int32, parameter :: bc_opts(3) = &
    [sll_p_dirichlet, sll_p_neumann, sll_p_neumann_mode_0]

  ! Local parameters
  sll_real64, parameter ::  one_third  = 1.0_f64 / 3.0_f64
  sll_real64, parameter :: four_thirds = 4.0_f64 / 3.0_f64

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
contains
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  !> Initialize the solver
  subroutine sll_s_qn_solver_2d_polar_init( solver, &
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

    type(sll_t_qn_solver_2d_polar), intent(out) :: solver !< Poisson solver class
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

    character(len=*), parameter :: this_sub_name = 'sll_s_qn_solver_2d_polar_init'

    sll_real64 :: dr
    sll_real64 :: inv_r
    sll_real64 :: inv_dr
    sll_real64 :: c
    sll_real64 :: d1(-1:+1)
    sll_real64 :: d2(-1:+1)

    sll_int32  :: i, j, k
    sll_int32  :: bc(2)
    sll_int32  :: last

    ! Consistency check: boundary conditions must be one of three options
    if( any( bc_rmin == bc_opts ) ) then
      bc(1) = bc_rmin
    else
      SLL_ERROR( this_sub_name, 'Unrecognized boundary condition at r_min' )
    end if
    !
    if( any( bc_rmax == bc_opts ) ) then
      bc(2) = bc_rmax
    else
      SLL_ERROR( this_sub_name, 'Unrecognized boundary condition at r_max' )
    end if

    ! Override vacuum permittivity in SI units
    if (present( epsilon_0 )) then
      solver%epsilon_0 = epsilon_0
    else
      solver%epsilon_0 = sll_p_epsilon_0
    end if

    ! Store information in solver
    solver%rmin  = rmin
    solver%rmax  = rmax
    solver%nr    = nr
    solver%nt    = ntheta
    solver%bc(:) = bc(:)

    ! Allocate arrays
    allocate( solver%g   (nr+1) )
    allocate( solver%z   (nr+1,ntheta) )
    allocate( solver%fk  (nr+1) )
    allocate( solver%phik(nr+1) )
    allocate( solver%mat((nr-1)*3,ntheta) ) ! for each k, matrix depends on r
    allocate( solver%cts((nr-1)*7) )
    allocate( solver%ipiv(nr-1) )

    ! Allocate in ALIGNED fashion 1D array for storing one row of z
    solver%zrow => sll_f_fft_allocate_aligned_complex( ntheta )

    ! Initialize plans for forward and backward FFTs
    call sll_s_fft_init_c2c_1d( solver%fw, &
      ntheta             , &
      solver%zrow(:)     , &
      solver%zrow(:)     , &
      sll_p_fft_forward  , &
      aligned    = .true., &
      normalized = .true. )

    call sll_s_fft_init_c2c_1d( solver%bw, &
      ntheta             , &
      solver%zrow(:)     , &
      solver%zrow(:)     , &
      sll_p_fft_backward , &
      aligned    = .true., &
      normalized = .false. )

    ! Store non-dimensional coefficient g(r) = \rho(r) / (B(r)^2 \epsilon_0)
    solver%g(:) = rho_m0(:) / (b_magn(:)**2 * solver%epsilon_0)

    ! Precompute convenient parameters
    dr     = (rmax-rmin)/nr
    inv_dr = 1.0_f64 / dr

    ! Store matrix coefficients into solver%mat
    ! Cycle over k_j
    do j = 1, ntheta

      ! Determine value of k_j (careful: ordering is not trivial)
      if (j <= ntheta/2) then
        k = j-1
      else
        k = ntheta-(j-1)
      end if

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
      if (bc(1) == sll_p_dirichlet) then ! Dirichlet
        solver%mat(1,j) = 0.0_f64
      else if (bc(1) == sll_p_neumann) then ! Neumann
        solver%mat(3,j) = solver%mat(3,j) -  one_third  * solver%mat(1,j)
        solver%mat(2,j) = solver%mat(2,j) + four_thirds * solver%mat(1,j)
        solver%mat(1,j) = 0.0_f64
      else if (bc(1) == sll_p_neumann_mode_0) then 
        if (k == 0) then ! Neumann for mode zero
          solver%mat(3,j) = solver%mat(3,j) -  one_third  * solver%mat(1,j)
          solver%mat(2,j) = solver%mat(2,j) + four_thirds * solver%mat(1,j)
          solver%mat(1,j) = 0.0_f64
        else             ! Dirichlet for other modes
          solver%mat(1,j) = 0.0_f64
        endif
      endif

      ! Set boundary condition at rmax
      last = 3*(nr-1)
      if (bc(2) == sll_p_dirichlet) then ! Dirichlet
        solver%mat(last,j) = 0.0_f64
      else if (bc(2) == sll_p_neumann) then ! Neumann
        solver%mat(last-2,j) = solver%mat(last-2,j) -  one_third  *solver%mat(last,j)
        solver%mat(last-1,j) = solver%mat(last-1,j) + four_thirds *solver%mat(last,j)
        solver%mat(last  ,j) = 0.0_f64
      else if (bc(2) == sll_p_neumann_mode_0) then 
        if (k == 0) then ! Neumann for mode zero
          solver%mat(last-2,j) = solver%mat(last-2,j) -  one_third  *solver%mat(last,j)
          solver%mat(last-1,j) = solver%mat(last-1,j) + four_thirds *solver%mat(last,j)
          solver%mat(last  ,j) = 0.0_f64
        else             ! Dirichlet for other modes
          solver%mat(last,j) = 0.0_f64
        endif
      endif

    end do

  end subroutine sll_s_qn_solver_2d_polar_init


  !=============================================================================
  !> Solve the quasi-neutrality equation and get the electrostatic potential
  subroutine sll_s_qn_solver_2d_polar_solve( solver, rho, phi )
    type(sll_t_qn_solver_2d_polar) , intent(inout) :: solver   !< Solver object
    sll_real64                     , intent(in   ) :: rho(:,:) !< Charge density
    sll_real64                     , intent(  out) :: phi(:,:) !< Potential

    sll_int32  :: nr, ntheta, bc(2)
    sll_int32  :: i, j, k

    nr     = solver%nr
    ntheta = solver%nt
    bc     = solver%bc

    ! Consistency check: 'rho' and 'phi' have shape defined at initialization
    SLL_ASSERT( all( shape(rho) == [nr+1,ntheta] ) )
    SLL_ASSERT( all( shape(phi) == [nr+1,ntheta] ) )

    ! Initialize 2D complex array for FFTs along theta and FD solver along r
    solver%z(:,:) = cmplx( rho(:,:), 0.0_f64, kind=f64 )

    ! For each r_i, compute FFT of rho(r_i,theta) to obtain \hat{rho}(r_i,k)
    do i = 1, nr+1
      solver%zrow(:) = solver%z(i,:)
      call sll_s_fft_exec_c2c_1d( solver%fw, solver%zrow(:), solver%zrow(:) )
      solver%z(i,:) = solver%zrow(:)
    end do

    ! Cycle over k_j
    do j = 1, ntheta

      ! Determine value of k_j (careful: ordering is not trivial)
      if (k <= ntheta/2) then
        k = j-1
      else
        k = ntheta-(j-1)
      endif

      ! Copy 1D slice of \hat{rho} into separate array and divide it by (g*eps0)
      solver%fk(:) = solver%z(:,j) / (solver%g(:) * solver%epsilon_0)

      ! Solve tridiagonal system to obtain \hat{phi}_{k_j}(r) at internal points
      call sll_s_setup_cyclic_tridiag( solver%mat(:,j), nr-1, solver%cts, solver%ipiv )
      call sll_o_solve_cyclic_tridiag( solver%cts, solver%ipiv, solver%fk(2:nr), &
        nr-1, solver%phik(2:nr) )

      ! Boundary condition at rmin
      if (bc(1) == sll_p_dirichlet) then ! Dirichlet
        solver%phik(1) = (0.0_f64, 0.0_f64)
      else if (bc(1) == sll_p_neumann) then ! Neumann
        solver%phik(1) = four_thirds*solver%phik(2) - one_third*solver%phik(3)
      else if (bc(1) == sll_p_neumann_mode_0) then 
        if (k==0) then ! Neumann for mode zero
          solver%phik(1) = four_thirds*solver%phik(2) - one_third*solver%phik(3)
        else           ! Dirichlet for other modes
          solver%phik(1) = (0.0_f64, 0.0_f64)
        endif
      endif

      ! Boundary condition at rmax
      if (bc(2) == sll_p_dirichlet) then ! Dirichlet
        solver%phik(nr+1) = (0.0_f64, 0.0_f64)
      else if (bc(2) == sll_p_neumann) then ! Neumann
        solver%phik(nr+1) = four_thirds*solver%phik(nr) - one_third*solver%phik(nr-1)
      else if (bc(2) == sll_p_neumann_mode_0) then 
        if(k==0) then ! Neumann for mode zero
          solver%phik(nr+1) = four_thirds*solver%phik(nr) - one_third*solver%phik(nr-1)
        else          ! Dirichlet for other modes
          solver%phik(nr+1) = (0.0_f64, 0.0_f64)
        endif
      endif

      ! Store \hat{phi}_{k_j}(r) into 1D slice of 2D array
      solver%z(:,j) = solver%phik(:)

    end do

    ! For each r_i, compute inverse FFT of \hat{phi}(r_i,k) to obtain phi(r_i,theta)
    do i = 1, nr+1
      solver%zrow(:) = solver%z(i,:)
      call sll_s_fft_exec_c2c_1d( solver%bw, solver%zrow(:), solver%zrow(:) )
      phi(i,:) = real( solver%zrow(:) )
    end do

  end subroutine sll_s_qn_solver_2d_polar_solve


  !=============================================================================
  !> Delete contents (local storage) of quasi-neutrality solver
  subroutine sll_s_qn_solver_2d_polar_free( solver )
    type(sll_t_qn_solver_2d_polar), intent(inout) :: solver

    call sll_s_fft_free( solver%fw )
    call sll_s_fft_free( solver%bw )

    deallocate( solver%z    )
    deallocate( solver%zrow )
    deallocate( solver%fk   )
    deallocate( solver%phik )
    deallocate( solver%mat  )
    deallocate( solver%cts  )
    deallocate( solver%ipiv )

  end subroutine sll_s_qn_solver_2d_polar_free


end module sll_m_qn_solver_2d_polar
