program test_poisson_2d_fem_sps
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"

  use sll_m_working_precision, only: f64

  use sll_m_constants, only: &
    sll_p_pi, &
    sll_p_twopi

  use sll_m_utilities, only: sll_s_new_array_linspace

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_greville, &
    sll_p_periodic

  use sll_m_bsplines, only: &
    sll_c_bsplines, &
    sll_s_bsplines_new

  use sll_m_spline_interpolator_1d, only: sll_s_spline_1d_compute_num_cells

  use sll_m_spline_2d, only: sll_t_spline_2d

  use sll_m_spline_interpolator_2d, only: sll_t_spline_interpolator_2d

  use sll_m_polar_mapping_analytical, only: sll_c_polar_mapping_analytical

  use sll_m_polar_mapping_analytical_target, only: sll_t_polar_mapping_analytical_target

  use sll_m_polar_mapping_analytical_czarny, only: sll_t_polar_mapping_analytical_czarny

  use sll_m_polar_mapping_iga, only: sll_t_polar_mapping_iga

  use sll_m_poisson_2d_fem_sps_stencil, only: sll_t_poisson_2d_fem_sps_stencil

  use sll_m_poisson_2d_fem_sps_stencil_new, only: sll_t_poisson_2d_fem_sps_stencil_new

  use sll_m_boundary_condition_descriptors, only: sll_p_dirichlet

  use sll_m_timer, only: &
    sll_t_time_mark    , &
    sll_s_set_time_mark, &
    sll_f_time_elapsed_between

  use sll_m_hdf5_io_serial, only: &
    sll_t_hdf5_ser_handle     , &
    sll_s_hdf5_ser_file_create, &
    sll_s_hdf5_ser_file_close , &
    sll_o_hdf5_ser_write_array, &
    sll_o_hdf5_ser_write_attribute

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Working precision
  integer, parameter :: wp = f64

  ! To initialize B-splines (p1,p2 degrees)
  integer :: mm, n1, n2, p1, p2, ncells1, ncells2, maptype!, nn

  ! B-splines break points
  real(wp), allocatable :: breaks_eta1(:)
  real(wp), allocatable :: breaks_eta2(:)

  ! 1D B-splines
  class(sll_c_bsplines), allocatable :: bsplines_eta1
  class(sll_c_bsplines), allocatable :: bsplines_eta2

  ! Analytical and discrete mappings
  class(sll_c_polar_mapping_analytical), allocatable :: mapping_analytic
  type(sll_t_polar_mapping_iga) :: mapping_discrete
  real(wp) :: d0, e0, x0(2)

  ! 2D splines representing right hand side and solution
  type(sll_t_spline_2d) :: spline_2d_rhs
  type(sll_t_spline_2d) :: spline_2d_phi

  ! 2D spline interpolator
  type(sll_t_spline_interpolator_2d) :: spline_interp_2d

  ! Needed for 2D interpolation of right hand side
  real(wp), allocatable :: tau_eta1(:)
  real(wp), allocatable :: tau_eta2(:)
  real(wp), allocatable :: gtau(:,:)

  ! Poisson solver
!  type(sll_t_poisson_2d_fem_sps_stencil    ) :: solver
  type(sll_t_poisson_2d_fem_sps_stencil_new) :: solver

!  ! Stiffness/mass dense matrices and C1 projections
!  real(wp), allocatable :: A (:,:)
!  real(wp), allocatable :: M (:,:)
!  real(wp), allocatable :: Ap(:,:)
!  real(wp), allocatable :: Mp(:,:)

  ! Auxiliary variables
  integer  :: i1, i2
  real(wp) :: eta(2), x(2)

  ! Timing
  type(sll_t_time_mark) :: t0, t1
  real(wp) :: dt

  ! For hdf5 I/O
  type(sll_t_hdf5_ser_handle) :: file_id
  integer                     :: h5_error

  real(wp), allocatable :: phi_spl(:,:), err(:,:)
  real(wp) :: err_L2_norm, err_Linf_norm, dV, phi_exact

  !-----------------------------------------------------------------------------
  ! Initialize B-splines basis functions
  !-----------------------------------------------------------------------------

  call sll_s_set_time_mark( t0 )

  ! Number of degrees of freedom (control points) along s and theta
  mm = 128
  n1 = mm
  n2 = mm * 2
!  nn = 3 + (n1-2) * n2

  ! Spline degrees along s and theta
  p1 = 3
  p2 = 3

  ! Mapping type: 0 circle, 1 target (no manufactured solution for Czarny's mapping)
  maptype = 1

  ! Create HDF5 file
  call sll_s_hdf5_ser_file_create( 'poisson_2d_fem_sps.h5', file_id, h5_error )

  ! HDF5 I/O
  call sll_o_hdf5_ser_write_attribute( file_id, "/", "n1", n1, h5_error )
  call sll_o_hdf5_ser_write_attribute( file_id, "/", "n2", n2, h5_error )
  call sll_o_hdf5_ser_write_attribute( file_id, "/", "p1", p1, h5_error )
  call sll_o_hdf5_ser_write_attribute( file_id, "/", "p2", p2, h5_error )

  ! Compute number of cells from number of interpolation points along s
  call sll_s_spline_1d_compute_num_cells( &
    degree  = p1            , &
    bc_xmin = sll_p_greville, &
    bc_xmax = sll_p_greville, &
    nipts   = n1            , &
    ncells  = ncells1 )

  ! Construct break points along s to initialize non-uniform spline basis
  allocate( breaks_eta1( ncells1+1 ) )
  call sll_s_new_array_linspace( breaks_eta1, 0.0_wp, 1.0_wp, endpoint=.true. )

  ! Create 1D spline basis along s in [0,1]
  call sll_s_bsplines_new( &
    bsplines = bsplines_eta1, &
    degree   = p1           , &
    periodic = .false.      , &
    xmin     = 0.0_wp       , &
    xmax     = 1.0_wp       , &
    ncells   = ncells1      , &
    breaks   = breaks_eta1 )

  ! Compute number of cells from number of interpolation points along theta
  call sll_s_spline_1d_compute_num_cells( &
    degree  = p2            , &
    bc_xmin = sll_p_periodic, &
    bc_xmax = sll_p_periodic, &
    nipts   = n2            , &
    ncells  = ncells2 )

  ! Construct break points along theta
  allocate( breaks_eta2( ncells2+1 ) )
  call sll_s_new_array_linspace( breaks_eta2, 0.0_wp, sll_p_twopi, endpoint=.true. )

  ! Create 1D spline basis along theta in [0,2pi]
  call sll_s_bsplines_new( &
    bsplines = bsplines_eta2, &
    degree   = p2           , &
    periodic = .true.       , &
    xmin     = 0.0_wp       , &
    xmax     = sll_p_twopi  , &
    ncells   = ncells2 )

  !-----------------------------------------------------------------------------
  ! Initialize mapping and polar B-splines
  !-----------------------------------------------------------------------------

  ! Allocate analytical mapping
  allocate( sll_t_polar_mapping_analytical_target :: mapping_analytic )

  ! Initialize analytical mapping
  select type ( mapping_analytic )
    type is ( sll_t_polar_mapping_analytical_target )
      if ( maptype == 0 ) then
        x0 = (/ 0.0_wp, 0.0_wp /)
        d0 = 0.0_wp
        e0 = 0.0_wp
        call mapping_analytic % init( x0, d0, e0 )
      else if ( maptype == 1 ) then
        x0 = (/ 0.0_wp, 0.0_wp /)
        d0 = 0.2_wp
        e0 = 0.3_wp
        call mapping_analytic % init( x0, d0, e0 )
      end if
   end select

  ! Discrete mapping
  call mapping_discrete % init( bsplines_eta1, bsplines_eta2, mapping_analytic )

  ! Write mapping info
  call mapping_discrete % store_data( n1, n2, file_id )

  call sll_s_set_time_mark( t1 )

  dt = sll_f_time_elapsed_between( t0, t1 )
  write(*,'(/a,es8.1/)') " Time required for initialization of B-splines and mapping: ", dt

  !-----------------------------------------------------------------------------
  ! Interpolate right hand side in tensor-product space
  !-----------------------------------------------------------------------------

  call sll_s_set_time_mark( t0 )

  ! Initialize 2D spline for right hand side
  call spline_2d_rhs % init( bsplines_eta1, bsplines_eta2 )

  ! Initialize 2D spline interpolator
  call spline_interp_2d % init( &
    bsplines_eta1, &
    bsplines_eta2, &
    [sll_p_greville,sll_p_periodic], &
    [sll_p_greville,sll_p_periodic] )

  ! Get interpolation points from 2D spline interpolator
  call spline_interp_2d % get_interp_points( tau_eta1, tau_eta2 )

  ! Evaluate right hand side on interpolation points
  associate( nt1 => size( tau_eta1 ), nt2 => size( tau_eta2 ) )

    allocate( gtau( nt1, nt2 ) )

    if ( maptype == 0 ) then

      do i2 = 1, nt2
        do i1 = 1, nt1
          eta = (/ tau_eta1(i1), tau_eta2(i2) /)
          x   = mapping_discrete % eval( eta )
          gtau(i1,i2) = rhs_cart( x )
        end do
      end do

     else if ( maptype == 1 ) then

       do i2 = 1, nt2
         do i1 = 1, nt1
           eta = (/ tau_eta1(i1), tau_eta2(i2) /)
           gtau(i1,i2) = rhs_log( eta )
         end do
       end do

     end if

  end associate

  ! Compute interpolant spline
  call spline_interp_2d % compute_interpolant( spline_2d_rhs, gtau )

  call sll_s_set_time_mark( t1 )

  dt = sll_f_time_elapsed_between( t0, t1 )
  write(*,'(a,es8.1/)' ) " Time required for interpolation of right hand side: ", dt

  !-----------------------------------------------------------------------------
  ! Poisson solver
  !-----------------------------------------------------------------------------

  ! Initialize 2D spline for solution
  call spline_2d_phi % init( bsplines_eta1, bsplines_eta2 )

  call sll_s_set_time_mark( t0 )

  ! Initialize Poisson solver
  call solver % init( bsplines_eta1, bsplines_eta2, breaks_eta1, breaks_eta2, mapping_discrete )

  call sll_s_set_time_mark( t1 )

  dt = sll_f_time_elapsed_between( t0, t1 )
  write(*,'(a,es8.1/)' ) " Time required for initialization of Poisson solver: ", dt

!  ! Allocate stiffness/mass dense matrices and C1 projections
!  allocate( A( n1*n2, n1*n2 ) )
!  allocate( M( n1*n2, n1*n2 ) )
!  allocate( Ap( nn, nn ) )
!  allocate( Mp( nn, nn ) )
!
!  ! Convert stencil to dense
!  A = 0.0_wp
!  M = 0.0_wp
!  call solver % A_linop_stencil % to_array( A )
!  call solver % M_linop_stencil % to_array( M )
!
!  ! Convert C1 block to dense
!  Ap = 0.0_wp
!  Mp = 0.0_wp
!  call solver % Ap_linop_c1_block % to_array( Ap )
!  call solver % Mp_linop_c1_block % to_array( Mp )

  call solver % set_boundary_conditions( sll_p_dirichlet )

  call sll_s_set_time_mark( t0 )

  ! Solve
  call solver % reset_charge()
  ! TODO: test properly signature with callable function
  call solver % accumulate_charge( spline_2d_rhs ) ! Right hand side is 2D spline
  call solver % solve( spline_2d_phi )

  call sll_s_set_time_mark( t1 )

  dt = sll_f_time_elapsed_between( t0, t1 )
  write(*,'(a,es8.1/)' ) " Time required for solution of Poisson equation: ", dt

  ! Evaluate phi spline on logical grid and compute spatial L2-norm of error

  allocate( phi_spl( n1, n2 ) )
  allocate( err( n1, n2 ) )

  err_L2_norm   = 0.0_wp
  err_Linf_norm = 0.0_wp

  ! Integral volume in logical space: ds * dtheta
  dV = ( 1.0 / n1 ) * ( sll_p_twopi / n2 )

  do i2 = 1, n2
    do i1 = 1, n1

      eta(1) = real( i1-1, wp ) / real( n1-1, wp )
      eta(2) = real( i2-1, wp ) * sll_p_twopi / real( n2, wp )

      x = mapping_discrete % eval( eta )

      phi_spl(i1,i2) = spline_2d_phi % eval( eta(1), eta(2) )

      if ( maptype == 0 ) then
        phi_exact = ( 1.0_wp - eta(1)**2 ) * cos( sll_p_twopi * x(1) ) * sin( sll_p_twopi * x(2) )
      else if ( maptype == 1 ) then
        phi_exact = eta(1)**2 * ( 1.0_wp - eta(1)**2 ) * cos( eta(2) )
      end if

      err(i1,i2) = phi_spl(i1,i2) - phi_exact

      err_L2_norm = err_L2_norm + err(i1,i2)**2 * mapping_discrete % jdet( eta ) * dV

      err_Linf_norm = merge( abs( err(i1,i2) ), err_Linf_norm, abs( err(i1,i2) ) > err_Linf_norm )

    end do
  end do

  err_L2_norm = sqrt( err_L2_norm )

  !-----------------------------------------------------------------------------
  ! HDF5 I/O
  !-----------------------------------------------------------------------------

  call sll_s_set_time_mark( t0 )

!  ! Write stiffness matrix
!  call sll_o_hdf5_ser_write_array( file_id, A, "/A", h5_error )
!
!  ! Write mass matrix
!  call sll_o_hdf5_ser_write_array( file_id, M, "/M", h5_error )
!
!!  ! Write right hand side
!!  call sll_o_hdf5_ser_write_array( file_id, solver % b, "/b", h5_error )

  ! Write solution
  call sll_o_hdf5_ser_write_array( file_id, solver % x, "/x", h5_error )

!  ! Write C1 projection of stiffness matrix
!  call sll_o_hdf5_ser_write_array( file_id, Ap, "/Ap", h5_error )
!
!  ! Write C1 projection of mass matrix
!  call sll_o_hdf5_ser_write_array( file_id, Mp, "/Mp", h5_error )
!
!!  ! Write L matrix needed for projection
!!  call sll_o_hdf5_ser_write_array( file_id, solver % L, "/L", h5_error )

  ! Write reshaped solution
  call sll_o_hdf5_ser_write_array( file_id, spline_2d_phi % bcoef(1:n1,1:n2), "/phi", h5_error )

  ! Write phi spline
  call sll_o_hdf5_ser_write_array( file_id, phi_spl, "/phi_spl", h5_error )

  ! Write error
  call sll_o_hdf5_ser_write_array( file_id, err, "/error", h5_error )

  ! Close HDF5 file
  call sll_s_hdf5_ser_file_close ( file_id, h5_error )

  call sll_s_set_time_mark( t1 )

  dt = sll_f_time_elapsed_between( t0, t1 )
  write(*,'(a,es8.1/)') " Time required for writing HDF5 output: ", dt

  write(*,'(a,es8.2/)') " Spatial L2    norm of error: ", err_L2_norm

  write(*,'(a,es8.2/)') " Spatial L-inf norm of error: ", err_Linf_norm

  !-----------------------------------------------------------------------------
  ! Deallocations and free
  !-----------------------------------------------------------------------------

  deallocate( breaks_eta1 )
  deallocate( breaks_eta2 )

  deallocate( phi_spl )

!  deallocate( A )
!  deallocate( M )
!  deallocate( Ap )
!  deallocate( Mp )

  deallocate( mapping_analytic )

  call bsplines_eta1 % free()
  call bsplines_eta2 % free()

  call mapping_discrete % free()

  call spline_2d_rhs % free()
  call spline_2d_phi % free()

  call solver % free()

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  SLL_PURE function rhs_cart( x )
    real(wp), intent(in) :: x(2)
    real(wp) :: rhs_cart

    associate( sq => x(1)**2 + x(2)**2, xp => sll_p_twopi * x(1), yp => sll_p_twopi * x(2) )

      rhs_cart =  8.0_wp*sll_p_pi**2*(1.0_wp-sq)*cos(xp)*sin(yp)+4.0_wp*cos(xp)*sin(yp) &
                 -8.0_wp*sll_p_pi*(x(1)*sin(xp)*sin(yp)-x(2)*cos(xp)*cos(yp))

    end associate

  end function rhs_cart

  SLL_PURE function rhs_log( eta )
    real(wp), intent(in) :: eta(2)
    real(wp) :: rhs_log

    associate( s => eta(1), t => eta(2) )

      rhs_log = (8.0_wp*d0**3)*s**5+(-24.0_wp*d0**2*e0*cos(t)**3+20.0_wp*d0**2*e0*cos(t)                  &
                +24.0_wp*d0**2*cos(t)**3-20.0_wp*d0**2*cos(t))*s**4+(-8.0_wp*d0**3                        &
                -42.0_wp*d0*e0**2*cos(t)**2+16.0_wp*d0*e0**2-24.0_wp*d0*e0*(-cos(t)**2+1.0_wp)**2         &
                -12.0_wp*d0*e0*cos(t)**2-42.0_wp*d0*cos(t)**2+16.0_wp*d0)*s**3+(8.0_wp*d0**2*e0*cos(t)**3 &
                -12.0_wp*d0**2*e0*cos(t)-8.0_wp*d0**2*cos(t)**3+12.0_wp*d0**2*cos(t)-15.0_wp*e0**3*cos(t) &
                -12.0_wp*e0**2*cos(t)**3+9.0_wp*e0**2*cos(t)+12.0_wp*e0*cos(t)**3-9.0_wp*e0*cos(t)        &
                +15.0_wp*cos(t))*s**2+(10.0_wp*d0*e0**2*cos(t)**2-8.0_wp*d0*e0**2                         &
                -8.0_wp*d0*e0*(-cos(t)**2+1.0_wp)**2-20.0_wp*d0*e0*cos(t)**2+16.0_wp*d0*e0                &
                +10.0_wp*d0*cos(t)**2-8.0_wp*d0)*s-3.0_wp*cos(t)+3.0_wp*e0**3*cos(t)                      &
                +e0**2*(-4.0_wp*cos(t)**3+3.0_wp*cos(t))+e0*(4.0_wp*cos(t)**3-3.0_wp*cos(t))

      rhs_log = - rhs_log / ( (1.0_wp+e0)**2 * (2.0_wp*d0*s*cos(t)+e0-1.0_wp)**3 )

    end associate

  end function rhs_log

end program test_poisson_2d_fem_sps
