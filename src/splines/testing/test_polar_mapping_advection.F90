program test_polar_mapping_advection
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"

  use sll_m_working_precision, only: f64

  use sll_m_constants, only: &
    sll_p_pi, &
    sll_p_twopi

  use sll_m_utilities, only: sll_s_new_array_linspace

  use sll_m_polar_mapping_analytical, only: sll_c_polar_mapping_analytical

  use sll_m_polar_mapping_analytical_target, only: sll_t_polar_mapping_analytical_target

  use sll_m_polar_mapping_analytical_czarny, only: sll_t_polar_mapping_analytical_czarny

  use sll_m_polar_mapping_iga, only: sll_t_polar_mapping_iga

  use sll_m_jacobian_2d_pseudo_cartesian, only: sll_t_jacobian_2d_pseudo_cartesian

  use sll_m_polar_advector_base, only: sll_c_polar_advector

  use sll_m_polar_advector_constant, only: sll_t_polar_advector_constant

  use sll_m_polar_advector_rotating, only: sll_t_polar_advector_rotating

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_greville, &
    sll_p_periodic

  use sll_m_bsplines, only: &
    sll_c_bsplines, &
    sll_s_bsplines_new

  use sll_m_spline_2d, only: sll_t_spline_2d

  use sll_m_spline_interpolator_1d, only: sll_s_spline_1d_compute_num_cells 

  use sll_m_spline_interpolator_2d, only: sll_t_spline_interpolator_2d

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

  ! To initialize B-splines
  integer :: mm, n1, n2, deg1, deg2, ncells1, ncells2

  ! To initialize non-uniform 1D B-splines along eta1
  real(wp), allocatable :: breaks_eta1(:)

  ! 1D B-splines
  class(sll_c_bsplines), allocatable :: spline_basis_eta1
  class(sll_c_bsplines), allocatable :: spline_basis_eta2

  ! Analytical and discrete IGA mappings
  class(sll_c_polar_mapping_analytical), allocatable :: mapping_analytical
  type(sll_t_polar_mapping_iga) :: mapping_iga

  type(sll_t_jacobian_2d_pseudo_cartesian) :: jac_2d_pcart

  ! 2D tensor-product splines for advection fields A1,A2 and distribution function f
  type(sll_t_spline_2d) :: spline_2d_f
  type(sll_t_spline_2d) :: spline_2d_a1
  type(sll_t_spline_2d) :: spline_2d_a2

  ! 2D tensor-product spline interpolator
  type(sll_t_spline_interpolator_2d) :: spline_interpolator_2d

  ! Interpolation points and fields
  integer :: npts1, npts2
  real(wp), allocatable :: e1_node(:)
  real(wp), allocatable :: e2_node(:)
  real(wp), allocatable :: f(:,:), f_ex(:,:)
  real(wp), allocatable :: a1(:,:), a2(:,:)

  ! For 2D advection
  type(sll_t_polar_advector_rotating) :: advector

  ! Auxiliary variables
  integer  :: i1, i2, i, iterations, maxiter
  real(wp) :: eta(2), eta_new(2)
  real(wp) :: x(2), xi(2), a(2)
  real(wp) :: dt, dV, abs_tol, rel_tol
  real(wp) :: err, err_L2_norm_space, err_Linf_norm_space, err_Linf_norm_time_L2, err_Linf_norm_time_Li

  ! For hdf5 I/O
  type(sll_t_hdf5_ser_handle) :: file_id
  integer                     :: h5_error
  character(len=32)           :: attr_name

  !=============================================================================
  ! Initialize B-splines basis functions
  !=============================================================================

  mm = 64 * 16
  n1 = mm
  n2 = mm * 2

  deg1 = 3
  deg2 = 3

  ! Compute number of cells from number of interpolation points along eta1
  call sll_s_spline_1d_compute_num_cells( &
    degree  = deg1          , &
    bc_xmin = sll_p_greville, &
    bc_xmax = sll_p_greville, &
    nipts   = n1            , &
    ncells  = ncells1 )

  allocate( breaks_eta1( ncells1+1 ) )

  call sll_s_new_array_linspace( breaks_eta1, 0.0_wp, 1.0_wp, endpoint=.true. )

  ! Create 1D spline basis along eta1 in [0,1]
  call sll_s_bsplines_new( &
    bsplines = spline_basis_eta1, &
    degree   = deg1             , &
    periodic = .false.          , &
    xmin     = 0.0_wp           , &
    xmax     = 1.0_wp           , &
    ncells   = ncells1          , &
    breaks   = breaks_eta1 )

  deallocate( breaks_eta1 )

  ! Compute number of cells from number of interpolation points along eta2
  call sll_s_spline_1d_compute_num_cells( &
    degree  = deg2          , &
    bc_xmin = sll_p_periodic, &
    bc_xmax = sll_p_periodic, &
    nipts   = n2            , &
    ncells  = ncells2 )

  ! Create 1D spline basis along eta2 in [0,2pi]
  call sll_s_bsplines_new( &
    bsplines = spline_basis_eta2, &
    degree   = deg2             , &
    periodic = .true.           , &
    xmin     = 0.0_wp           , &
    xmax     = sll_p_twopi      , &
    ncells   = ncells2 )

  !=============================================================================
  ! Initialize standard 2D tensor-product splines and spline interpolator
  !=============================================================================

  ! Initialize 2D tensor-product splines
  call spline_2d_f  % init( spline_basis_eta1, spline_basis_eta2 )
  call spline_2d_a1 % init( spline_basis_eta1, spline_basis_eta2 )
  call spline_2d_a2 % init( spline_basis_eta1, spline_basis_eta2 )

  ! Initialize 2D tensor-product spline interpolator
  call spline_interpolator_2d % init( spline_basis_eta1, spline_basis_eta2, &
                                      [ sll_p_greville, sll_p_periodic ], &
                                      [ sll_p_greville, sll_p_periodic ] )

  ! Get interpolation points and allocate 2D array of values
  call spline_interpolator_2d % get_interp_points( e1_node, e2_node )

  npts1 = size( e1_node )
  npts2 = size( e2_node )
  allocate( f   ( npts1, npts2+1 ) ) ! repeated point along eta2
  allocate( f_ex( npts1, npts2+1 ) ) ! repeated point along eta2

  ! Create HDF5 file for output
  call sll_s_hdf5_ser_file_create( 'mapping_test_advection.h5', file_id, h5_error )

  !=============================================================================
  ! Initialize mapping
  !=============================================================================

!  allocate( sll_t_polar_mapping_analytical_target :: mapping_analytical )
  allocate( sll_t_polar_mapping_analytical_czarny :: mapping_analytical )

  ! Initialize analytica mapping
  select type ( mapping_analytical )
    type is ( sll_t_polar_mapping_analytical_target )
      call mapping_analytical % init() ! circular mapping
!      call mapping_analytical % init( x0=[0.0_wp,0.0_wp], d0=0.2_wp, e0=0.3_wp )
    type is ( sll_t_polar_mapping_analytical_czarny )
      call mapping_analytical % init( x0=[0.0_wp,0.0_wp], b =1.4_wp, e =0.3_wp )
  end select

  ! Initialize discrete IGA mapping
  call mapping_iga % init( spline_basis_eta1, spline_basis_eta2, mapping_analytical )

  call mapping_iga % store_data( npts1, npts2, file_id )

  !=============================================================================
  ! Set initial distribution function
  !=============================================================================

  ! Initial condition
  do i2 = 1, npts2
    do i1 = 1, npts1
      x = mapping_iga % eval( [ e1_node(i1), e2_node(i2) ] )
      f(i1,i2) = f_initial( x )
    end do
  end do

  ! Apply periodicity along eta2
  f(:,npts2+1) = f(:,1)

  ! Store initial solution
  call sll_o_hdf5_ser_write_array( file_id, f, "/f_0", h5_error )
  call sll_o_hdf5_ser_write_array( file_id, f, "/f_ex_0", h5_error )

  !=============================================================================
  ! Evolve f in time
  !=============================================================================

  ! Initialize advector
  call advector % init( xc=[0.25_wp,0.0_wp], omega=sll_p_twopi )
!  call advector % init( xc=[-0.15_wp,0.0_wp], omega=sll_p_twopi )

  ! Parameters for time cycle
  iterations = 10 * 16
  dt = 1.0e-01_wp / 16.0_wp

  call sll_o_hdf5_ser_write_attribute( file_id, "/", "iterations", iterations, h5_error )

  ! Compute interpolating splines for advection fields
  allocate( a1( npts1, npts2 ) )
  allocate( a2( npts1, npts2 ) )
  do i2 = 1, npts2
    do i1 = 1, npts1
      x = mapping_iga % eval( [ e1_node(i1), e2_node(i2) ] )
      call advector % velocity_field( x, a )
      a1(i1,i2) = a(1)
      a2(i1,i2) = a(2)
    end do
  end do
  call spline_interpolator_2d % compute_interpolant( spline_2d_a1, a1(:,:) )
  call spline_interpolator_2d % compute_interpolant( spline_2d_a2, a2(:,:) )

  ! Initialize Jacobian of pseudo-Cartesian coordinates
  call jac_2d_pcart % init( mapping_iga )

!  ! Set tolerance and maximum iterations for mapping inversion
!  tol = 1.0e-14_wp
!  maxiter = 100

  abs_tol = 1.0e-14_wp
  rel_tol = 1.0e-12_wp
  maxiter = 100

  ! Integral volume in logical space: ds * dtheta
  dV = ( 1.0 / n1 ) * ( sll_p_twopi / n2 )

  err_Linf_norm_time_L2 = 0.0_wp
  err_Linf_norm_time_Li = 0.0_wp

  ! Time cycle
  do i = 1, iterations

    ! Compute interpolating spline for f
    call spline_interpolator_2d % compute_interpolant( spline_2d_f, f(:,1:npts2) )

    err_L2_norm_space   = 0.0_wp
    err_Linf_norm_space = 0.0_wp

    ! 2D advection
    do i2 = 1, npts2
      do i1 = 1, npts1

        eta(1) = e1_node(i1)
        eta(2) = e2_node(i2)

        ! Map logical coordinates to Cartesian coordinates using discrete IGA mapping
        x = mapping_iga % eval( eta )

!        write(*,'(a,es21.14,a,es21.14,a)') " (e1,e2) = (", eta(1), ",", eta(2), ") at t"
!        write(*,'(a,es21.14,a,es21.14,a)') " (x1,x2) = (",   x(1), ",",   x(2), ") at t"

!        ! Advect point using analytical flow field
!        x_new   = advector % flow_field( x, -dt )
!        eta_new = mapping_iga % eval_inverse( x_new, eta, tol, maxiter )

!        ! Advect point using explicit Runge-Kutta 4th order (Cartesian coordinates)
!        eta_new = advector % advect_x1x2( eta, -dt, mapping_iga, spline_2d_a1, spline_2d_a2, tol, maxiter )

        ! Advect point using explicit Runge-Kutta 4th order (intermediate coordinates)
        eta_new = advector % advect_y1y2( &
          eta         , &
          -dt         , &
          jac_2d_pcart, &
          spline_2d_a1, &
          spline_2d_a2, &
          rel_tol     , &
          abs_tol     , &
          maxiter )

!        write(*,'(a,es21.14,a,es21.14,a)') " (x1,x2) = (", x_new(1), ",", x_new(2), ") at t+dt"
!        write(*,'(a,es21.14,a,es21.14,a)') " (e1,e2) = (", eta_new(1), ",", eta_new(2), ") at t+dt"

        ! Evaluate distribution function at origin of characteristics
        f(i1,i2) = spline_2d_f % eval( eta_new(1), eta_new(2) )

        ! Exact solution using method of characteristics
        xi = advector % flow_field( x, -dt*real(i,wp) )
        f_ex(i1,i2) = f_initial( xi )

        err = f(i1,i2) - f_ex(i1,i2)

        err_L2_norm_space = err_L2_norm_space + err**2 * abs( mapping_iga % jdet( eta ) ) * dV

        err = abs( f(i1,i2) - f_ex(i1,i2) )

        err_Linf_norm_space = merge( err, err_Linf_norm_space, err > err_Linf_norm_space )

      end do
    end do

    err_L2_norm_space = sqrt( err_L2_norm_space )

    ! Compute L-inf norm of error in time of L2 norm of error in space
    err_Linf_norm_time_L2 = merge( err_L2_norm_space, err_Linf_norm_time_L2, err_L2_norm_space > err_Linf_norm_time_L2 )

    ! Compute L-inf norm of error in time of L-inf norm of error in space
    err_Linf_norm_time_Li = merge( err_Linf_norm_space, err_Linf_norm_time_Li, err_Linf_norm_space > err_Linf_norm_time_Li )

    ! Apply periodicity along eta2
    f   (:,npts2+1) = f   (:,1)
    f_ex(:,npts2+1) = f_ex(:,1)

    ! Store solution
    write( attr_name, '(a,i0)' ) "/f_", i 
    call sll_o_hdf5_ser_write_array( file_id, f, trim(attr_name), h5_error )

    ! Store exact solution
    write( attr_name, '(a,i0)' ) "/f_ex_", i 
    call sll_o_hdf5_ser_write_array( file_id, f_ex, trim(attr_name), h5_error )

  end do

  write(*,*)
  write(*,'(a,es8.2)') " Maximum in time of spatial L2 norm of error = ", err_Linf_norm_time_L2
  write(*,*)

  write(*,*)
  write(*,'(a,es8.2)') " Maximum in time of spatial L-inf norm of error = ", err_Linf_norm_time_Li
  write(*,*)

  ! Close HDF5 file
  call sll_s_hdf5_ser_file_close ( file_id, h5_error )

  !=============================================================================
  ! Deallocate allocatables and free objects
  !=============================================================================

  deallocate( e1_node )
  deallocate( e2_node )
  deallocate( f    )
  deallocate( f_ex )
  deallocate( a1 )
  deallocate( a2 )

  deallocate( mapping_analytical )

  call mapping_iga  % free()
  call spline_2d_f  % free()
  call spline_2d_a1 % free()
  call spline_2d_a2 % free()
  call spline_interpolator_2d % free()
  call spline_basis_eta1 % free()
  call spline_basis_eta2 % free()
  call jac_2d_pcart % free()
  call advector % free()

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !-----------------------------------------------------------------------------
  ! Initial distribution function
  !-----------------------------------------------------------------------------
  SLL_PURE function f_initial( x ) result( f0 )
    real(wp), intent(in) :: x(2)
    real(wp) :: f0

!    ! Cosine bell
!
!    real(wp), parameter :: x0(2) = [0.0_wp,0.0_wp]
!    integer , parameter :: p = 4
!    real(wp) :: xmod
!
!    xmod = norm2( x-x0 )
!
!    if ( xmod < 0.25_wp ) then
!      f0 = cos( sll_p_twopi * (x(1)-x0(1)) )**p * cos( sll_p_twopi * (x(2)-x0(2)) )**p
!    else
!      f0 = 0.0_wp
!    end if

    ! Superposition of cosine bells with elliptical cross-sections

    real(wp), parameter :: x0(2) = [-0.15_wp,0.0_wp]
!    real(wp), parameter :: x0(2) = [0.25_wp,0.0_wp]
    integer , parameter :: p = 4
    real(wp), parameter :: a = 0.3_wp
    real(wp) :: r1, r2

    r1 = sqrt( ( x(1) - x0(1) )**2 + 8.0_wp * ( x(2) - x0(2) )**2 )
    r2 = sqrt( 8.0_wp * ( x(1) - x0(1) )**2 + ( x(2) - x0(2) )**2 )

    f0 = 0.0_wp

    if ( r1 < a ) f0 = 0.5_wp * cos( 0.5_wp * sll_p_pi * r1 / a )**p
    if ( r2 < a ) f0 = f0 + 0.5_wp * cos( 0.5_wp * sll_p_pi * r2 / a )**p

  end function f_initial

end program test_polar_mapping_advection
