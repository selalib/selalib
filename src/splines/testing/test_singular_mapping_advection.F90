program test_singular_mapping_advection
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"

  use sll_m_working_precision, only: f64

  use sll_m_constants, only: &
    sll_p_pi, &
    sll_p_twopi

  use sll_m_utilities, only: sll_s_new_array_linspace

  use sll_m_singular_mapping_analytic, only: sll_c_singular_mapping_analytic

  use sll_m_singular_mapping_analytic_target, only: sll_t_singular_mapping_analytic_target

  use sll_m_singular_mapping_analytic_czarny, only: sll_t_singular_mapping_analytic_czarny

  use sll_m_singular_mapping_discrete, only: sll_t_singular_mapping_discrete

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

  use sll_m_gauss_legendre_integration, only: &
    sll_f_gauss_legendre_points, &
    sll_f_gauss_legendre_weights

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

  ! Tolerance valid for n1=128, n2=256, p1=p2=3, dt=0.05 and iter=20 in input file
  ! (needs to be adjusted in case of different input values)
  real(wp), parameter :: tolerance = 1e-1_wp

  ! Integer variables
  integer :: mm, n1, n2, p1, p2, ncells1, ncells2, ntau1, ntau2, Nk1, Nk2, Nq1, Nq2, &
             i1, i2, k1, k2, q1, q2, i, iter, maxiter, maptype, file_er, file_unit

  !Real variables
  real(wp) :: dt, abs_tol, rel_tol, eta(2), eta_new(2), x(2), xi(2), a(2), &
              err, L2_norm_space, Linf_norm_space, L2_error, Linf_error, f_ex_quad

  ! Character variables
  character(len=32) :: attr_name
  character(len=:), allocatable :: input_file

  ! Logical variables
  logical :: success

  ! Real 1D allocatables
  real(wp), allocatable :: breaks_eta1(:), breaks_eta2(:), tau_eta1(:), tau_eta2(:)

  ! Real 2D allocatables
  real(wp), allocatable :: f(:,:), f_ex(:,:), Ax(:,:), Ay(:,:), quad_points_eta1(:,:), &
                           quad_points_eta2(:,:), quad_weights_eta1(:,:), quad_weights_eta2(:,:)

  ! Real 4D allocatables
  real(wp), allocatable :: volume(:,:,:,:)

  ! Abstract polymorphic types
  class(sll_c_bsplines)                 , allocatable :: bsplines_eta1, bsplines_eta2
  class(sll_c_singular_mapping_analytic), allocatable :: mapping_analytic

  ! Concrete types
  type(sll_t_singular_mapping_discrete)    :: mapping_discrete
  type(sll_t_spline_2d)                    :: spline_2d_f, spline_2d_Ax, spline_2d_Ay
  type(sll_t_spline_interpolator_2d)       :: spline_interp_2d
  type(sll_t_jacobian_2d_pseudo_cartesian) :: jac_2d_pcart
  type(sll_t_polar_advector_rotating)      :: advector
  type(sll_t_hdf5_ser_handle)              :: file_id

  ! Namelists
  namelist /splines/ &
    n1, &
    n2, &
    p1, &
    p2

  namelist /geometry/ &
    maptype

  namelist /time_integration/ &
    dt, &
    iter

  ! Read input file
  input_file = trim( "test_singular_mapping_advection.nml" )
  open( file=input_file, status='old', action='read', newunit=file_unit )
  read( file_unit, splines  )        ; rewind( file_unit )
  read( file_unit, geometry )        ; rewind( file_unit )
  read( file_unit, time_integration );  close( file_unit )

  !-----------------------------------------------------------------------------
  ! Initialize B-splines basis functions
  !-----------------------------------------------------------------------------

  ! Create HDF5 file for output
  call sll_s_hdf5_ser_file_create( 'mapping_test_advection.h5', file_id, file_er )

  ! HDF5 I/O
  call sll_o_hdf5_ser_write_attribute( file_id, "/", "n1", n1, file_er )
  call sll_o_hdf5_ser_write_attribute( file_id, "/", "n2", n2, file_er )
  call sll_o_hdf5_ser_write_attribute( file_id, "/", "p1", p1, file_er )
  call sll_o_hdf5_ser_write_attribute( file_id, "/", "p2", p2, file_er )

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

  ! Compute number of cells from number of interpolation points along eta2
  call sll_s_spline_1d_compute_num_cells( &
    degree  = p2          , &
    bc_xmin = sll_p_periodic, &
    bc_xmax = sll_p_periodic, &
    nipts   = n2            , &
    ncells  = ncells2 )

  ! Construct break points along theta
  allocate( breaks_eta2( ncells2+1 ) )
  call sll_s_new_array_linspace( breaks_eta2, 0.0_wp, sll_p_twopi, endpoint=.true. )

  ! Create 1D spline basis along eta2 in [0,2pi]
  call sll_s_bsplines_new( &
    bsplines = bsplines_eta2, &
    degree   = p2           , &
    periodic = .true.       , &
    xmin     = 0.0_wp       , &
    xmax     = sll_p_twopi  , &
    ncells   = ncells2 )

  !-----------------------------------------------------------------------------
  ! Initialize standard 2D tensor-product splines and spline interpolator
  !-----------------------------------------------------------------------------

  ! Initialize 2D tensor-product splines
  call spline_2d_f  % init( bsplines_eta1, bsplines_eta2 )
  call spline_2d_Ax % init( bsplines_eta1, bsplines_eta2 )
  call spline_2d_Ay % init( bsplines_eta1, bsplines_eta2 )

  ! Initialize 2D tensor-product spline interpolator
  call spline_interp_2d % init( bsplines_eta1, bsplines_eta2, &
                                      [ sll_p_greville, sll_p_periodic ], &
                                      [ sll_p_greville, sll_p_periodic ] )

  ! Get interpolation points and allocate 2D array of values
  call spline_interp_2d % get_interp_points( tau_eta1, tau_eta2 )

  ntau1 = size( tau_eta1 )
  ntau2 = size( tau_eta2 )

  ! Repeated point in theta
  allocate( f   ( ntau1, ntau2+1 ) )
  allocate( f_ex( ntau1, ntau2+1 ) )

  !-----------------------------------------------------------------------------
  ! Initialize mapping
  !-----------------------------------------------------------------------------

  ! Allocate analytical mapping
  if ( maptype == 0 .or. maptype == 1 ) then
    allocate( sll_t_singular_mapping_analytic_target :: mapping_analytic )
  else if ( maptype == 2 ) then
    allocate( sll_t_singular_mapping_analytic_czarny :: mapping_analytic )
  end if

  ! Initialize analytical mapping
  select type ( mapping_analytic )
    type is ( sll_t_singular_mapping_analytic_target )
      if ( maptype == 0 ) then
        call mapping_analytic % init( x0=[0.0_wp,0.0_wp], d0=0.0_wp, e0=0.0_wp )
      else if ( maptype == 1 ) then
        call mapping_analytic % init( x0=[0.0_wp,0.0_wp], d0=0.2_wp, e0=0.3_wp )
      end if
    type is ( sll_t_singular_mapping_analytic_czarny )
      call mapping_analytic % init( x0=[0.0_wp,0.0_wp], b=1.4_wp, e=0.3_wp )
  end select

  ! Initialize discrete mapping
  call mapping_discrete % init( bsplines_eta1, bsplines_eta2, mapping_analytic )

  ! Write mapping info
  call mapping_discrete % store_data( n1, n2, file_id )

  !-----------------------------------------------------------------------------
  ! Set initial distribution function
  !-----------------------------------------------------------------------------

  ! Initial condition
  do i2 = 1, ntau2
    do i1 = 1, ntau1
      eta(1) = tau_eta1(i1)
      eta(2) = tau_eta2(i2)
      x = mapping_discrete % eval( eta )
      f(i1,i2) = f_initial( x )
    end do
  end do
  ! Apply periodicity along eta2
  f(:,ntau2+1) = f(:,1)

  ! Store initial solution
  call sll_o_hdf5_ser_write_array( file_id, f, "/f_0", file_er )
  call sll_o_hdf5_ser_write_array( file_id, f, "/f_ex_0", file_er )

  ! Compute interpolating spline for f
  call spline_interp_2d % compute_interpolant( spline_2d_f, f(:,1:ntau2) )

  !-----------------------------------------------------------------------------
  ! Evolve f in time
  !-----------------------------------------------------------------------------

  ! Initialize advector
  call advector % init( xc=[0.25_wp,0.0_wp], omega=sll_p_twopi )

  call sll_o_hdf5_ser_write_attribute( file_id, "/", "iterations", iter, file_er )

  ! Compute interpolating splines for advection fields
  allocate( Ax( ntau1, ntau2 ) )
  allocate( Ay( ntau1, ntau2 ) )

  do i2 = 1, ntau2
    do i1 = 1, ntau1
      eta(1) = tau_eta1(i1)
      eta(2) = tau_eta2(i2)
      x = mapping_discrete % eval( eta )
      call advector % velocity_field( x, a )
      Ax(i1,i2) = a(1)
      Ay(i1,i2) = a(2)
    end do
  end do

  call spline_interp_2d % compute_interpolant( spline_2d_Ax, Ax(:,:) )
  call spline_interp_2d % compute_interpolant( spline_2d_Ay, Ay(:,:) )

  ! Initialize Jacobian of pseudo-Cartesian coordinates
  call jac_2d_pcart % init( mapping_discrete )
  call jac_2d_pcart % pole( tau_eta2 )

  ! Quadrature points to compute L2-norm of error
  Nk1 = ncells1
  Nk2 = ncells2
  Nq1 = 1 + p1
  Nq2 = 1 + p2

  allocate( quad_points_eta1 ( Nq1, Nk1 ) )
  allocate( quad_points_eta2 ( Nq2, Nk2 ) )
  allocate( quad_weights_eta1( Nq1, Nk1 ) )
  allocate( quad_weights_eta2( Nq2, Nk2 ) )

  allocate( volume( Nq1, Nq2, Nk1, Nk2 ) )

  ! Quadrature points and weights along s
  do k1 = 1, Nk1
    quad_points_eta1 (:,k1) = sll_f_gauss_legendre_points ( Nq1, breaks_eta1(k1), breaks_eta1(k1+1) )
    quad_weights_eta1(:,k1) = sll_f_gauss_legendre_weights( Nq1, breaks_eta1(k1), breaks_eta1(k1+1) )
  end do

  ! Quadrature points and weights along theta
  do k2 = 1, Nk2
    quad_points_eta2 (:,k2) = sll_f_gauss_legendre_points ( Nq2, breaks_eta2(k2), breaks_eta2(k2+1) )
    quad_weights_eta2(:,k2) = sll_f_gauss_legendre_weights( Nq2, breaks_eta2(k2), breaks_eta2(k2+1) )
  end do

  do k2 = 1, Nk2
    do k1 = 1, Nk1
      do q2 = 1, Nq2
        do q1 = 1, Nq1
          eta(1) = quad_points_eta1(q1,k1)
          eta(2) = quad_points_eta2(q2,k2)
          volume(q1,q2,k1,k2) = abs( mapping_discrete % jdet( eta ) ) &
                                * quad_weights_eta1(q1,k1) * quad_weights_eta2(q2,k2)
        end do
      end do
    end do
  end do

  L2_error   = 0.0_wp
  Linf_error = 0.0_wp

  ! Time cycle
  do i = 1, iter

    Linf_norm_space = 0.0_wp

    do i2 = 1, ntau2
      do i1 = 1, ntau1

        eta(1) = tau_eta1(i1)
        eta(2) = tau_eta2(i2)

        ! Map logical coordinates to Cartesian coordinates using discrete IGA mapping
        x = mapping_discrete % eval( eta )

!        ! Advect point using analytical flow field
!        x_new   = advector % flow_field( x, -dt )
!        eta_new = mapping_discrete % eval_inverse( x_new, eta, tol, maxiter )

!        ! Advect point using Cartesian coordinates (RK3)
!        eta_new = advector % advect_cart( eta, -dt, mapping_discrete, spline_2d_Ax, spline_2d_Ay )

        ! Advect point using pseudo-Cartesian coordinates (RK3)
        eta_new = advector % advect_pseudo_cart( eta, -dt, jac_2d_pcart, spline_2d_Ax, spline_2d_Ay )

        ! Evaluate distribution function at origin of characteristics
        f(i1,i2) = spline_2d_f % eval( eta_new(1), eta_new(2) )

        ! Exact solution using method of characteristics
        xi = advector % flow_field( x, -dt*real(i,wp) )
        f_ex(i1,i2) = f_initial( xi )

        err = f(i1,i2) - f_ex(i1,i2)

        Linf_norm_space = merge( abs( err ), Linf_norm_space, abs( err ) > Linf_norm_space )

      end do
    end do

    ! Compute interpolating spline for f
    call spline_interp_2d % compute_interpolant( spline_2d_f, f(:,1:ntau2) )

    ! Compute spatial L2-norm of error
    L2_norm_space   = 0.0_wp

    do k2 = 1, Nk2
      do k1 = 1, Nk1
        do q2 = 1, Nq2
          do q1 = 1, Nq1

            eta = (/ quad_points_eta1(q1,k1), quad_points_eta2(q2,k2) /)
            x   = mapping_discrete % eval( eta )

            xi = advector % flow_field( x, -dt*real(i,wp) )
            f_ex_quad = f_initial( xi )

            L2_norm_space = L2_norm_space + volume(q1,q2,k1,k2) &
                            * ( spline_2d_f % eval( eta(1), eta(2) ) - f_ex_quad )**2

          end do
        end do
      end do
    end do

    L2_norm_space = sqrt( L2_norm_space )

    ! Compute L-inf norm in time of spatial L2 norm of error
    L2_error = merge( L2_norm_space, L2_error, L2_norm_space > L2_error )

    ! Compute L-inf norm in time of spatial L-inf norm of error
    Linf_error = merge( Linf_norm_space, Linf_error, Linf_norm_space > Linf_error )

    ! Apply periodicity along eta2
    f   (:,ntau2+1) = f   (:,1)
    f_ex(:,ntau2+1) = f_ex(:,1)

    ! Store solution
    write( attr_name, '(a,i0)' ) "/f_", i 
    call sll_o_hdf5_ser_write_array( file_id, f, trim(attr_name), file_er )

    ! Store exact solution
    write( attr_name, '(a,i0)' ) "/f_ex_", i 
    call sll_o_hdf5_ser_write_array( file_id, f_ex, trim(attr_name), file_er )

  end do

  ! Close HDF5 file
  call sll_s_hdf5_ser_file_close ( file_id, file_er )

  write(*,'(/a,es8.2)')  " Maximum in time of spatial   L2-norm of error: ", L2_error
  write(*,'(/a,es8.2/)') " Maximum in time of spatial Linf-norm of error: ", Linf_error

  ! Check if test passed
  success = .true.
  if( max( L2_error, Linf_error ) > tolerance ) then
    success = .false.
    write(*,'(a,es7.1,a/)') "Test FAILED (tolerance ", tolerance, ")"
  end if
  if( success ) write(*,'(a/)') "Test PASSED"

  !-----------------------------------------------------------------------------
  ! Deallocate allocatables and free objects
  !-----------------------------------------------------------------------------

  deallocate( breaks_eta1, breaks_eta2, tau_eta1, tau_eta2, f, f_ex, Ax, Ay, &
              quad_points_eta1, quad_points_eta2, quad_weights_eta1, quad_weights_eta2, volume )

  deallocate( mapping_analytic )

  call mapping_discrete  % free()
  call bsplines_eta1 % free()
  call bsplines_eta2 % free()
  call spline_2d_f  % free()
  call spline_2d_Ax % free()
  call spline_2d_Ay % free()
  call spline_interp_2d % free()
  call jac_2d_pcart % free()
  call advector % free()

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Initial distribution function (two cosine bells with elliptical cross-sections)
  SLL_PURE function f_initial( x ) result( f0 )
    real(wp), intent(in) :: x(2)
    real(wp) :: f0

    real(wp), parameter :: x0(2) = [-0.15_wp,0.0_wp]
    integer , parameter :: p = 4
    real(wp), parameter :: a = 0.3_wp
    real(wp) :: r1, r2

    r1 = sqrt( ( x(1) - x0(1) )**2 + 8.0_wp * ( x(2) - x0(2) )**2 )
    r2 = sqrt( 8.0_wp * ( x(1) - x0(1) )**2 + ( x(2) - x0(2) )**2 )

    f0 = 0.0_wp

    if ( r1 < a ) f0 = 0.5_wp * cos( 0.5_wp * sll_p_pi * r1 / a )**p
    if ( r2 < a ) f0 = f0 + 0.5_wp * cos( 0.5_wp * sll_p_pi * r2 / a )**p

  end function f_initial

end program test_singular_mapping_advection
