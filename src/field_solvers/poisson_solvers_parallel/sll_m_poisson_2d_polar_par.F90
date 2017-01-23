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

!> @ingroup poisson_solvers
!> @author  Yaman Güçlü, IPP Garching
!> @author  Edoardo Zoni, IPP Garching
!>
!> This module solves the Poisson equation \f$ -\nabla^2 \phi = \rho \f$
!> (permittivity of free space \f$ \epsilon_0 = 1 \f$) on a 2D polar mesh with
!> coordinates \f$ (r,\theta) \f$:
!> \f[
!> -\bigg(
!> \frac{\partial^2}{\partial r^{2}}
!> +\frac{1}{r}\frac{\partial}{\partial r}
!> +\frac{1}{r^{2}}\frac{\partial^{2}}{\partial\theta^{2}}
!> \bigg)\phi(r,\theta) = \rho(r,\theta).
!> \f]
!> The boundary conditions on \f$ \phi(r,\theta) \f$ are as follows:
!> \f$ 2\pi \f$-periodicity along \f$ \theta \f$;
!> Neumann mode 0 at \f$ r = r_\textrm{min} \f$:
!> \f$ \partial_r \widehat{\phi}_0(r_\textrm{min}) = 0 \f$
!> and \f$ \widehat{\phi}_k(r_\textrm{min})=0 \f$ for \f$ k\neq 0 \f$;
!> homogeneous Dirichlet at \f$ r = r_\textrm{max} \f$:
!> \f$ \phi(r_\textrm{max},\theta) = 0 \f$.
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

module sll_m_poisson_2d_polar_par
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"
#include "sll_working_precision.h"

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
    sll_t_poisson_2d_polar_par, &
    sll_s_poisson_2d_polar_par_init, &
    sll_s_poisson_2d_polar_par_solve, &
    sll_s_poisson_2d_polar_par_free

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Class for the Poisson solver in polar coordinate
  type sll_t_poisson_2d_polar_par

   sll_real64              :: rmin       !< Min value of r coordinate
   sll_real64              :: rmax       !< Max value of r coordinate
   sll_int32               :: nr         !< Number of cells along r
   sll_int32               :: ntheta     !< Number of cells along theta
   sll_int32               :: bc(2)      !< Boundary conditions options

   type(sll_t_fft)         :: fw         !< Forward FFT plan
   type(sll_t_fft)         :: bw         !< Inverse FFT plan
   sll_real64, pointer     :: tmp (:)    !< 1D work array for FFT, real
   sll_int32 , allocatable :: k_list(:)
   sll_real64, allocatable :: z_r (:,:)  !< 2D array sequential in r
   sll_real64, pointer     :: z_a (:,:)  !< 2D array sequential in theta
   sll_real64, allocatable :: mat (:,:)  !< Tridiagonal matrix (one for each k)
   sll_real64, allocatable :: cts (:)    !< Lapack coefficients
   sll_int32 , allocatable :: ipiv(:)    !< Lapack pivot indices

   type(sll_t_layout_2d)           , pointer :: layout_r !< layout sequential in r
   type(sll_t_layout_2d)           , pointer :: layout_a !< layout sequential in theta
   type(sll_t_remap_plan_2d_real64), pointer :: rmp_ra   !< remap r->theta 
   type(sll_t_remap_plan_2d_real64), pointer :: rmp_ar   !< remap theta->r

  end type sll_t_poisson_2d_polar_par

  ! Allowed boundary conditions
  sll_int32, parameter :: bc_opts(3) = &
    [sll_p_dirichlet, sll_p_neumann, sll_p_neumann_mode_0]

  ! Local parameters
  sll_real64, parameter ::  one_third  = 1.0_f64 / 3.0_f64
  sll_real64, parameter :: four_thirds = 4.0_f64 / 3.0_f64

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
contains
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

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
      bc_rmax )

    type(sll_t_poisson_2d_polar_par), intent(out) :: solver !< solver object
    type(sll_t_layout_2d), pointer    :: layout_r !< layout sequential in r direction
    type(sll_t_layout_2d), pointer    :: layout_a !< layout sequential in theta direction
    sll_real64           , intent(in) :: rmin     !< rmin
    sll_real64           , intent(in) :: rmax     !< rmax
    sll_int32            , intent(in) :: nr       !< number of cells radial
    sll_int32            , intent(in) :: ntheta   !< number of cells angular
    sll_int32            , intent(in) :: bc_rmin  !< boundary condition at r_min
    sll_int32            , intent(in) :: bc_rmax  !< boundary condition at r_max

    character(len=*), parameter :: this_sub_name = 'sll_s_poisson_2d_polar_par_init'

    sll_real64 :: dr
    sll_real64 :: inv_r
    sll_real64 :: inv_dr

    sll_int32  :: loc_sz_r(2) ! sequential in r direction
    sll_int32  :: loc_sz_a(2) ! sequential in theta direction
    sll_int32  :: i, j, k
    sll_int32  :: bck(2)
    sll_int32  :: last
    sll_int32  :: glob_idx(2)

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

    ! Store global information in solver
    solver%rmin     =  rmin
    solver%rmax     =  rmax
    solver%nr       =  nr
    solver%ntheta   =  ntheta
    solver%layout_a => layout_a
    solver%layout_r => layout_r

    ! Allocate arrays global in r
    allocate( solver%z_r (nr+1,   loc_sz_r(2)) )
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

    ! Precompute convenient parameters
    dr = (rmax-rmin)/nr
    inv_dr = 1.0_f64/dr

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
        inv_r = 1.0_f64 / (rmin + (i-1)*dr)
        solver%mat(3*(i-1)  ,j) =        -inv_dr**2 - 0.5_f64*inv_dr*inv_r
        solver%mat(3*(i-1)-1,j) = 2.0_f64*inv_dr**2 + (k*inv_r)**2
        solver%mat(3*(i-1)-2,j) =        -inv_dr**2 + 0.5_f64*inv_dr*inv_r
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

  end subroutine sll_s_poisson_2d_polar_par_init


  !=============================================================================
  !> Solve the Poisson equation and get the electrostatic potential
  subroutine sll_s_poisson_2d_polar_par_solve( solver, rho, phi )
    type(sll_t_poisson_2d_polar_par) , intent(inout) :: solver   !< Solver object
    sll_real64                       , intent(in   ) :: rho(:,:) !< Charge density
    sll_real64, target               , intent(  out) :: phi(:,:) !< Potential

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
      associate( rhok => solver%z_r(:,j), phik => solver%z_r(:,j) )

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

  end subroutine sll_s_poisson_2d_polar_par_free


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

end module sll_m_poisson_2d_polar_par
