program sim_bsl_gc_2d0v_smooth_polar_splines
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"

  use sll_m_working_precision, only: f64

  use sll_m_constants, only: &
    sll_p_pi, &
    sll_p_twopi

  use sll_m_utilities, only: sll_s_new_array_linspace

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_periodic, &
    sll_p_greville, &
    sll_p_dirichlet

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

  use sll_m_poisson_2d_fem_sps_stencil_new, only: sll_t_poisson_2d_fem_sps_stencil_new

  use sll_m_timer, only: &
    sll_t_time_mark    , &
    sll_s_set_time_mark, &
    sll_f_time_elapsed_between

  use sll_m_hdf5_io_serial, only: &
    sll_t_hdf5_ser_handle     , &
    sll_s_hdf5_ser_file_create, &
    sll_s_hdf5_ser_file_open  , &
    sll_s_hdf5_ser_file_close , &
    sll_o_hdf5_ser_write_array, &
    sll_o_hdf5_ser_write_attribute, &
    sll_o_hdf5_ser_read_array

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Working precision
  integer, parameter :: wp = f64

  ! Integer variables
  integer :: n1, n2, p1, p2, ncells1, ncells2, ntau1, ntau2, i1, i2, it, maptype, h5_error, file_unit

  ! Real variables
  real(wp) :: sigma, inf_norm_phi, residual, t_diff, t_iter, eta(2)

!  real(wp), parameter :: max_phi = 0.11268967746569288_wp
  real(wp), parameter :: max_phi = 0.113_wp

  ! Namelists

  namelist /splines/ &
    n1, &
    n2, &
    p1, &
    p2

  namelist /geometry/ &
    maptype

  real(wp), parameter :: tolerance = 1.0e-12_wp

  ! Character variables
  character(len=:), allocatable :: input_file
  character(len=32) :: attr_name

  ! Real 1D allocatables
  real(wp), allocatable :: breaks_eta1(:), breaks_eta2(:), tau_eta1(:), tau_eta2(:)

  ! Real 2D allocatables
  real(wp), allocatable :: rho(:,:), phi(:,:)

  ! Abstract polymorphic types
  class(sll_c_bsplines)                , allocatable :: bsplines_eta1, bsplines_eta2
  class(sll_c_polar_mapping_analytical), allocatable :: mapping_analytic

  ! Concrete types
  type(sll_t_polar_mapping_iga)              :: mapping_discrete
  type(sll_t_spline_2d)                      :: spline_2d_rho, spline_2d_phi
  type(sll_t_spline_interpolator_2d)         :: spline_interp_2d
  type(sll_t_poisson_2d_fem_sps_stencil_new) :: poisson_solver
  type(sll_t_time_mark)                      :: t0, t1
  type(sll_t_hdf5_ser_handle)                :: file_id

  ! Parse input argument
  call s_parse_command_arguments( input_file )

  ! Read input file
  open ( file=trim( input_file ), status='old', action='read', newunit=file_unit )
  read ( file_unit, splines  ); rewind( file_unit )
  read ( file_unit, geometry ); close ( file_unit )

  ! Create HDF5 file
  call sll_s_hdf5_ser_file_create( 'sim_bsl_gc_2d0v_smooth_polar_splines_equilibrium.h5', file_id, h5_error )

  ! Write data to HDF5 file
  call sll_o_hdf5_ser_write_attribute( file_id, "/", "n1", n1, h5_error )
  call sll_o_hdf5_ser_write_attribute( file_id, "/", "n2", n2, h5_error )
  call sll_o_hdf5_ser_write_attribute( file_id, "/", "p1", p1, h5_error )
  call sll_o_hdf5_ser_write_attribute( file_id, "/", "p2", p2, h5_error )

  write(*,'(/a)') " >> Initializing B-splines"

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

  write(*,'(/a)') " >> Initializing mapping"

  ! Allocate analytical mapping
  if ( maptype == 0 .or. maptype == 1 ) then
    allocate( sll_t_polar_mapping_analytical_target :: mapping_analytic )
  else if ( maptype == 2 ) then
    allocate( sll_t_polar_mapping_analytical_czarny :: mapping_analytic )
  end if

  ! Initialize analytical mapping
  select type ( mapping_analytic )
    type is ( sll_t_polar_mapping_analytical_target )
      if ( maptype == 0 ) then
        call mapping_analytic % init( x0=[0.0_wp,0.0_wp], d0=0.0_wp, e0=0.0_wp )
      else if ( maptype == 1 ) then
        call mapping_analytic % init( x0=[0.0_wp,0.0_wp], d0=0.2_wp, e0=0.3_wp )
      end if
    type is ( sll_t_polar_mapping_analytical_czarny )
      call mapping_analytic % init( x0=[0.0_wp,0.0_wp], b=1.4_wp, e=0.3_wp )
  end select

  ! Initialize discrete mapping
  call mapping_discrete % init( bsplines_eta1, bsplines_eta2, mapping_analytic )

  ! Write data to HDF5 file
  call mapping_discrete % store_data( n1, n2, file_id )

  write(*,'(/a)') " >> Initializing 2D splines and interpolator"

  ! Initialize 2D splines
  call spline_2d_rho % init( bsplines_eta1, bsplines_eta2 )
  call spline_2d_phi % init( bsplines_eta1, bsplines_eta2 )

  ! Initialize 2D spline interpolator
  call spline_interp_2d % init( &
    bsplines_eta1, &
    bsplines_eta2, &
    [sll_p_greville,sll_p_periodic], &
    [sll_p_greville,sll_p_periodic] )

  ! Get interpolation points from 2D spline interpolator
  call spline_interp_2d % get_interp_points( tau_eta1, tau_eta2 )

  ntau1 = size( tau_eta1 )
  ntau2 = size( tau_eta2 )

  ! Write data to HDF5 file
  call sll_o_hdf5_ser_write_attribute( file_id, "/", "ntau1", ntau1, h5_error )
  call sll_o_hdf5_ser_write_attribute( file_id, "/", "ntau2", ntau2, h5_error )

  write(*,'(/a)',advance='no') " >> Initializing Poisson solver"
 
  call sll_s_set_time_mark( t0 )

  ! Initialize Poisson solver
  call poisson_solver % init( &
    bsplines_eta1   , &
    bsplines_eta2   , &
    breaks_eta1     , &
    breaks_eta2     , &
    mapping_discrete )

  ! Set boundary conditions
  call poisson_solver % set_boundary_conditions( sll_p_dirichlet )

  call sll_s_set_time_mark( t1 )

  t_diff = sll_f_time_elapsed_between( t0, t1 )
  write(*,'(a,es7.1,a)') " ( ", t_diff, " s )"

  ! Repeated point along theta
  allocate( rho    ( ntau1, ntau2+1 ) )
  allocate( phi    ( ntau1, ntau2+1 ) )

  ! Initial guess for phi
  phi(:,1:ntau2) = 0.1_wp
  phi(:,ntau2+1) = phi(:,1)
  ! Write phi on interpolation points
  write( attr_name, '(a,i0)' ) "/phi_", 0
  call sll_o_hdf5_ser_write_array( file_id, phi, trim(attr_name), h5_error )

  ! Initial guess for sigma
  sigma = 1.0_wp

  ! Initial rho
  do i2 = 1, ntau2
    do i1 = 1, ntau1
      eta(1) = tau_eta1(i1)
      eta(2) = tau_eta2(i2)
      rho(i1,i2) = sigma * f( phi(i1,i2) )
    end do
  end do
  rho(:,ntau2+1) = rho(:,1)
  ! Write rho on interpolation points
  write( attr_name, '(a,i0)' ) "/rho_", 0
  call sll_o_hdf5_ser_write_array( file_id, rho, trim(attr_name), h5_error )

  ! Set residual to start Picard iteration
  residual  = 1.0_wp

  write(*,'(/a/)') " >> Entering Picard cycle"

  it = 0
  t_iter = 0.0_wp
  do while ( residual > tolerance )

    it = it+1
    write(*,'(a,i4)',advance='no') "    iteration ", it
 
    call sll_s_set_time_mark( t0 )

    ! rho_n
    call spline_interp_2d % compute_interpolant( spline_2d_rho, rho(:,1:ntau2) )

    ! phi_*
    call poisson_solver % reset_charge()
    call poisson_solver % accumulate_charge( spline_2d_rho )
    call poisson_solver % solve( spline_2d_phi )

    do i2 = 1, ntau2
      do i1 = 1, ntau1
        eta(1) = tau_eta1(i1)
        eta(2) = tau_eta2(i2)
        phi(i1,i2) = spline_2d_phi % eval( eta(1), eta(2) )
      end do
    end do
    phi(:,ntau2+1) = phi(:,1)

    ! || phi_* ||_inf
    inf_norm_phi = maxval( abs( phi(:,:) ) )

    ! phi_n = phi_* / || phi_* ||_inf
    spline_2d_phi % bcoef(:,:) = spline_2d_phi % bcoef(:,:) / inf_norm_phi * max_phi

    residual = sigma

    ! sigma_n = sigma_(n-1) / || phi_* ||_inf
    sigma = sigma / inf_norm_phi * max_phi

    residual = abs( sigma - residual )

    do i2 = 1, ntau2
      do i1 = 1, ntau1
        eta(1) = tau_eta1(i1)
        eta(2) = tau_eta2(i2)
        phi(i1,i2) = spline_2d_phi % eval( eta(1), eta(2) )
        rho(i1,i2) = sigma * f( phi(i1,i2) )
      end do
    end do
    phi(:,ntau2+1) = phi(:,1)
    rho(:,ntau2+1) = rho(:,1)

    ! Write rho on interpolation points
    write( attr_name, '(a,i0)' ) "/rho_", it
    call sll_o_hdf5_ser_write_array( file_id, rho, trim(attr_name), h5_error )

    ! Write phi on interpolation points
    write( attr_name, '(a,i0)' ) "/phi_", it
    call sll_o_hdf5_ser_write_array( file_id, phi, trim(attr_name), h5_error )

    call sll_s_set_time_mark( t1 )

    t_diff = sll_f_time_elapsed_between( t0, t1 )
    write(*,'(a,es7.1,a)',advance='no') " ( ", t_diff, " s )"

    t_iter = t_iter + t_diff

    write(*,'(a,es21.14)') " sigma = ", sigma

  end do

  call sll_o_hdf5_ser_write_attribute( file_id, "/", "iterations", it, h5_error )

  ! Average time per iteration
  if ( it /= 1 ) t_iter = t_iter / real(it-1,wp)
  write(*,'(/a,es7.1,a)') " >> Average time per iteration: ", t_iter, " s"

  ! Close HDF5 file
  call sll_s_hdf5_ser_file_close( file_id, h5_error )

  !-----------------------------------------------------------------------------
  ! Deallocations and free
  !-----------------------------------------------------------------------------

  ! Deallocate real 1D allocatables
  deallocate( breaks_eta1, breaks_eta2, tau_eta1, tau_eta2 )

  ! Deallocate real 2D allocatables
  deallocate( rho, phi )

  ! Free abstract polymorphic types
  call bsplines_eta1 % free()
  call bsplines_eta2 % free()

  ! Deallocate abstract polymorphic types
  deallocate( bsplines_eta1, bsplines_eta2, mapping_analytic )

  ! Free concrete types
  call mapping_discrete   % free()
  call spline_2d_rho      % free()
  call spline_2d_phi      % free()
  call spline_interp_2d   % free()
  call poisson_solver     % free()

  write(*,'(/a/)') " >> End of program"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine s_parse_command_arguments( input_file )
    character(len=:), allocatable, intent(  out) :: input_file

    character(len=*), parameter :: this_sub_name = "s_parse_command_arguments"

    integer :: argc
    character(len=256) :: string
    character(len=256) :: err_msg

    ! Count command line arguments
    argc = command_argument_count()

    ! Stop if there are no arguments
    if (argc == 0) then
      err_msg = "Input file not found"
      SLL_ERROR( this_sub_name, err_msg )
      return
    end if

    ! Read file name
    call get_command_argument( 1, string )
    string  = adjustl( string )
    input_file = trim( string )

  end subroutine s_parse_command_arguments

  !-----------------------------------------------------------------------------
  function f( psi )
    real(wp), intent(in) :: psi 
    real(wp) :: f

!    f = psi

    real(wp) :: b, m, q, L, k, s

    b = 1.0e+15_wp

    m = 8.8_wp
    q = -0.22_wp

    L = 1.0_wp
    k = 1.0e+02_wp

    f = log( 1.0_wp + b**( m*psi + q ) +  L / ( 1.0_wp + exp( -k * ( psi - 0.095_wp ) ) ) * &
             b**( 4.3_wp * ( -(max_phi-psi)**0.6_wp + (max_phi-0.025_wp)**0.6_wp ) ) ) / log( b )

  end function f

end program sim_bsl_gc_2d0v_smooth_polar_splines
