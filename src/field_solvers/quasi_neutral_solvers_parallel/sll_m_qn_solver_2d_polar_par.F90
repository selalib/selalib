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
    sll_s_fft_init_c2c_1d, &
    sll_s_fft_exec_c2c_1d, &
    sll_s_fft_free

  use sll_m_remapper, only: &
    sll_t_layout_2d, &
    sll_o_compute_local_sizes, &
    sll_o_get_layout_global_size_i, &
    sll_o_get_layout_global_size_j, &
    sll_o_local_to_global, &
    sll_t_remap_plan_2d_comp64, &
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

   sll_real64              :: rmin     !< Min value of r coordinate
   sll_real64              :: rmax     !< Max value of r coordinate
   sll_int32               :: nr       !< Number of cells along r
   sll_int32               :: nt       !< Number of cells along theta
   sll_int32               :: bc(2)    !< Boundary conditions options
   sll_real64              :: epsilon_0!< Vacuum permittivity (user may override)
   sll_real64, allocatable :: g(:)     !< g(r) = \rho_{m,0}/(B^2\epsilon_0)

   type(sll_t_fft)         :: fw       !< Forward FFT plan
   type(sll_t_fft)         :: bw       !< Inverse FFT plan
   sll_comp64, allocatable :: f_r (:,:)!< 2D array sequential in r
   sll_comp64, allocatable :: f_a (:,:)!< 2D array sequential in theta
   sll_comp64, allocatable :: fk  (:)  !< k-th Fourier mode of rho
   sll_comp64, allocatable :: phik(:)  !< k-th Fourier mode of phi
   sll_real64, allocatable :: mat (:,:)!< Tridiagonal matrix (one for each k)
   sll_real64, allocatable :: cts (:)  !< Lapack coefficients
   sll_int32 , allocatable :: ipiv(:)  !< Lapack pivot indices

   type(sll_t_layout_2d)           , pointer :: layout_r !< layout sequential in r
   type(sll_t_layout_2d)           , pointer :: layout_a !< layout sequential in theta
   type(sll_t_remap_plan_2d_comp64), pointer :: rmp_ra   !< remap r->theta 
   type(sll_t_remap_plan_2d_comp64), pointer :: rmp_ar   !< remap theta->r

  end type sll_t_qn_solver_2d_polar_par

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

    sll_real64              :: dr
    sll_real64              :: inv_r
    sll_real64              :: inv_dr
    sll_real64              :: c
    sll_real64              :: d1(-1:+1)
    sll_real64              :: d2(-1:+1)

    sll_comp64, allocatable :: buf(:)
    sll_int32               :: i, j, k
    sll_int32               :: bc(2)
    sll_int32               :: last

    sll_int32               :: loc_sz_r(2) ! local shape of layout_r
    sll_int32               :: loc_sz_a(2) ! local shape of layout_a
    sll_int32               :: glob_idx(2) ! global indices

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
    solver%nt       =  ntheta
    solver%bc(:)    =  bc(:)
    solver%layout_a => layout_a
    solver%layout_r => layout_r

    ! Allocate arrays global in r
    allocate( solver%g   (nr+1) )
    allocate( solver%fk  (nr+1) )
    allocate( solver%phik(nr+1) )
    allocate( solver%mat((nr-1)*3,loc_sz_r(2)) ) ! for each k, matrix depends on r
    allocate( solver%cts((nr-1)*7) )
    allocate( solver%ipiv(nr-1) )

    ! Array for FFTs in theta-direction
    allocate( solver%f_a (loc_sz_a(1), loc_sz_a(2)) )
    solver%f_a(:,:) = (0.0_f64, 0.0_f64)

    ! Array for FD solver in r-direction
    allocate( solver%f_r (loc_sz_r(1), loc_sz_r(2)) )
    solver%f_r(:,:) = (0.0_f64, 0.0_f64)

    ! Remap objects between two layouts (for transposing between f_r and f_a)
    solver%rmp_ra => sll_o_new_remap_plan( solver%layout_r, solver%layout_a, solver%f_r )
    solver%rmp_ar => sll_o_new_remap_plan( solver%layout_a, solver%layout_r, solver%f_a )

    ! Initialize plans for forward and backward FFTs
    allocate( buf(ntheta) )
    call sll_s_fft_init_c2c_1d( solver%fw, ntheta, buf, buf, sll_p_fft_forward, normalized=.true. )
    call sll_s_fft_init_c2c_1d( solver%bw, ntheta, buf, buf, sll_p_fft_backward )
    deallocate( buf )

    ! Store non-dimensional coefficient g(r) = \rho(r) / (B(r)^2 \epsilon_0)
    solver%g(:) = rho_m0(:) / (b_magn(:)**2 * solver%epsilon_0)

    ! Precompute convenient parameters
    dr     = (rmax-rmin)/nr
    inv_dr = 1.0_f64 / dr

    ! Store matrix coefficients into solver%mat
    ! Cycle over k_j
    do j = 1, loc_sz_r(2)

      ! Determine value of k_j (careful: ordering is not trivial)
      glob_idx(:) = sll_o_local_to_global( solver%layout_r, [1,j] )
      k = glob_idx(2)
      if (k <= ntheta/2) then
        k = k-1
      else
        k = ntheta-(k-1)
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

  end subroutine sll_s_qn_solver_2d_polar_par_init


  !=============================================================================
  !> Solve the quasi-neutrality equation and get the electrostatic potential
  subroutine sll_s_qn_solver_2d_polar_par_solve( solver, rho, phi )
    type(sll_t_qn_solver_2d_polar_par) , intent(inout) :: solver   !< Solver object
    sll_real64                         , intent(in   ) :: rho(:,:) !< Charge density
    sll_real64                         , intent(  out) :: phi(:,:) !< Potential

    sll_int32  :: nr, ntheta, bc(2)
    sll_int32  :: i, j, k
    sll_int32  :: glob_idx(2)

    nr     = solver%nr
    ntheta = solver%nt
    bc     = solver%bc

    ! Consistency check: rho and phi must be given in layout sequential in theta
    call verify_argument_sizes_par( solver%layout_a, rho )
    call verify_argument_sizes_par( solver%layout_a, phi )

    ! Copy charge into 2D complex array
    solver%f_a(:,:) = cmplx( rho, 0.0_f64, kind=f64 )

    ! For each r_i, compute FFT of rho(r_i,theta) to obtain \hat{rho}(r_i,k)
    do i = 1, ubound( solver%f_a, 1 )
      call sll_s_fft_exec_c2c_1d( solver%fw, solver%f_a(i,:), solver%f_a(i,:) )
    end do

    ! Remap \hat{rho}(k) to layout distributed in k (global in r) -> \hat{rho}_k(r)
    call sll_o_apply_remap_2d( solver%rmp_ar, solver%f_a, solver%f_r )

    ! Cycle over k_j
    do j = 1, ubound( solver%f_r, 2 )

      ! Determine value of k_j (careful: ordering is not trivial)
      glob_idx(:) = sll_o_local_to_global( solver%layout_r, [1,j] )
      k = glob_idx(2)
      if (k <= ntheta/2) then
        k = k-1
      else
        k = ntheta-(k-1)
      endif

      ! Copy 1D slice of \hat{rho} into separate array and divide it by (g*eps0)
      solver%fk(:) = solver%f_r(:,j) / (solver%g(:) * solver%epsilon_0)

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

      ! Store \hat{phi}_{k_j}(r) into 1D slice of 2D array (ready for remap)
      solver%f_r(:,j) = solver%phik(:)

    end do

    ! Redistribute \hat{phi}(r_i,k_j) into layout global in k
    call sll_o_apply_remap_2d( solver%rmp_ra, solver%f_r, solver%f_a )
    call verify_argument_sizes_par( solver%layout_a, phi )
    
    ! For each r_i, compute inverse FFT of \hat{phi}(r_i,k) to obtain phi(r_i,theta)
    do i = 1, ubound( solver%f_a, 1 )
      call sll_s_fft_exec_c2c_1d( solver%bw, solver%f_a(i,:), solver%f_a(i,:) )
      phi(i,:) = real(solver%f_a(i,:))
    end do

  end subroutine sll_s_qn_solver_2d_polar_par_solve


  !=============================================================================
  !> Delete contents (local storage) of quasi-neutrality solver
  subroutine sll_s_qn_solver_2d_polar_par_free( solver )
    type(sll_t_qn_solver_2d_polar_par), intent(inout) :: solver

    call sll_s_fft_free( solver%fw )
    call sll_s_fft_free( solver%bw )

    deallocate( solver%f_r  )
    deallocate( solver%f_a  )
    deallocate( solver%fk   )
    deallocate( solver%phik )
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
