program sim_bsl_gc_2d0v_smooth_polar_splines
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"

  use sll_m_working_precision, only: f64

  use sll_m_constants, only: sll_p_twopi

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

  use sll_m_electric_field, only: sll_t_electric_field

  use sll_m_diagnostics, only: sll_t_diagnostics

  use sll_m_point_charge, only: sll_t_point_charge

  use sll_m_simulation_state, only: sll_t_simulation_state

  use sll_m_time_integrator, only: sll_t_time_integrator

  use sll_m_timer, only: &
    sll_t_time_mark    , &
    sll_s_set_time_mark, &
    sll_f_time_elapsed_between

  use sll_m_hdf5_io_serial, only: &
    sll_t_hdf5_ser_handle         , &
    sll_s_hdf5_ser_file_create    , &
    sll_s_hdf5_ser_file_open      , &
    sll_s_hdf5_ser_file_close     , &
    sll_o_hdf5_ser_write_array    , &
    sll_o_hdf5_ser_write_attribute, &
    sll_o_hdf5_ser_read_array     , &
    sll_o_hdf5_ser_read_attribute

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Working precision
  integer, parameter :: wp = f64

  ! Integer variables
  integer :: n1, n2, p1, p2, ncells1, ncells2, ntau1, ntau2, nx1, nx2, i1, i2, it, &
             iter, diag_freq, maxiter, l, p, maptype, h5_error, file_unit, nc, ic

  ! Real variables
  real(wp) :: dt, abs_tol, rel_tol, smin, smax, ampl, t_diff, t_iter, eta(2)

  ! Logical variables
  logical :: equil_num, success

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

  namelist /characteristics/ &
    abs_tol, &
    rel_tol, &
    maxiter

  namelist /initial_condition/ &
    l   , &
    p   , &
    smin, &
    smax, &
    ampl

  namelist /equilibrium/ &
    equil_num

  namelist /diagnostics/ &
    diag_freq, &
    nx1      , &
    nx2

  namelist /point_charges/ &
    nc

  namelist /point_charges_specs/ &
    intensity, &
    location

  ! Real parameters
  real(wp), parameter :: epsi = 1.0e-12_wp

  ! Character variables
  character(len=:), allocatable :: input_file
  character(len=32) :: file_name
  character(len=10) :: status
  character(len=10) :: position
  character(len=32) :: attr_name

  ! Real 1D allocatables
  real(wp), allocatable :: breaks_eta1(:), breaks_eta2(:), tau_eta1(:), tau_eta2(:), intensity(:)

  ! Real 2D allocatables
  real(wp), allocatable :: phi(:,:), location(:,:)

  ! Abstract polymorphic types
  class(sll_c_bsplines)                , allocatable :: bsplines_eta1, bsplines_eta2
  class(sll_c_polar_mapping_analytical), allocatable :: mapping_analytic

  ! Concrete types
  type(sll_t_polar_mapping_iga)              :: mapping_discrete
  type(sll_t_spline_2d)                      :: spline_2d_rho, spline_2d_phi
  type(sll_t_spline_interpolator_2d)         :: spline_interp_2d
  type(sll_t_poisson_2d_fem_sps_stencil_new) :: poisson_solver
  type(sll_t_electric_field)                 :: electric_field
  type(sll_t_simulation_state)               :: sim_state
  type(sll_t_time_integrator)                :: time_integrator
  type(sll_t_diagnostics)                    :: diag
  type(sll_t_time_mark)                      :: t0, t1
  type(sll_t_hdf5_ser_handle)                :: file_id, file_id_eq

  ! Parse input argument
  call s_parse_command_arguments( input_file )

  ! Read input file
  open( file=trim( input_file ), status='old', action='read', newunit=file_unit )
  read( file_unit, splines             ); rewind( file_unit )
  read( file_unit, geometry            ); rewind( file_unit )
  read( file_unit, time_integration    ); rewind( file_unit )
  read( file_unit, characteristics     ); rewind( file_unit )
  read( file_unit, initial_condition   ); rewind( file_unit )
  read( file_unit, equilibrium         ); rewind( file_unit )

  ! for point charges: reading number of point charges, intesities and positions
  read( file_unit, point_charges       ); rewind( file_unit )
  allocate( intensity(    nc ) )
  allocate( location ( 2, nc ) )
  read( file_unit, point_charges_specs ); rewind( file_unit )

  read( file_unit, diagnostics         ); close ( file_unit )

  ! Create HDF5 file
  call sll_s_hdf5_ser_file_create( 'sim_bsl_gc_2d0v_smooth_polar_splines.h5', file_id, h5_error )

  ! Write data to HDF5 file
  call sll_o_hdf5_ser_write_attribute( file_id, "/", "n1" , n1 , h5_error )
  call sll_o_hdf5_ser_write_attribute( file_id, "/", "n2" , n2 , h5_error )
  call sll_o_hdf5_ser_write_attribute( file_id, "/", "p1" , p1 , h5_error )
  call sll_o_hdf5_ser_write_attribute( file_id, "/", "p2" , p2 , h5_error )
  call sll_o_hdf5_ser_write_attribute( file_id, "/", "nx1", nx1, h5_error )
  call sll_o_hdf5_ser_write_attribute( file_id, "/", "nx2", nx2, h5_error )
  call sll_o_hdf5_ser_write_attribute( file_id, "/", "nc" , nc , h5_error )

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

  ! Initialize electric field
  call electric_field % init( mapping_discrete, spline_2d_phi )

  ! Repeated point along theta
  allocate( phi( ntau1, ntau2+1 ) )

  ! Initialize simulation state
  call sim_state % init( &
    ntau1        , &
    ntau2        , &
    nc           , &
    intensity    , &
    location     , &
    spline_2d_rho, &
    spline_2d_phi, &
    electric_field )

  ! Compute equilibrium density on interpolation points
  if ( equil_num ) then

    call sll_s_hdf5_ser_file_open( 'sim_bsl_gc_2d0v_smooth_polar_splines_equilibrium.h5', file_id_eq, h5_error )
    call sll_o_hdf5_ser_read_attribute( file_id_eq, "/", "iterations", it, h5_error )
    write( attr_name, '(a,i0)' ) "/rho_", it
    call sll_o_hdf5_ser_read_array( file_id_eq, sim_state % rho, trim(attr_name), h5_error )
    call sll_o_hdf5_ser_read_array( file_id_eq, phi, trim(attr_name), h5_error )
!    call sll_s_hdf5_ser_file_close( file_id_eq, h5_error )

    ! Write data to HDF5 file
    call sll_o_hdf5_ser_write_array( file_id, phi, "/phi_eq", h5_error )

  else

    do i2 = 1, ntau2
      do i1 = 1, ntau1
        eta(1) = tau_eta1(i1)
        eta(2) = tau_eta2(i2)
        sim_state % rho(i1,i2) = rho_equilibrium( eta )
      end do
    end do
    ! Apply periodicity along theta
    sim_state % rho(:,ntau2+1) = sim_state % rho(:,1)

  end if

  ! Write data to HDF5 file
  call sll_o_hdf5_ser_write_array( file_id, sim_state % rho, "/rho_eq", h5_error )

  ! Add perturbation to equilibrium
  do i2 = 1, ntau2
    do i1 = 1, ntau1
      eta(1) = tau_eta1(i1)
      eta(2) = tau_eta2(i2)
      sim_state % rho(i1,i2) = sim_state % rho(i1,i2) + rho_perturbation( eta )
    end do
  end do
  ! Apply periodicity along theta
  sim_state % rho(:,ntau2+1) = sim_state % rho(:,1)

  ! Initialize time integrator
   call time_integrator % init( &
     dt              , &
     tau_eta1        , &
     tau_eta2        , &
     mapping_discrete, &
     spline_interp_2d, &
     sim_state       , &
     poisson_solver  , &
     abs_tol         , &
     rel_tol         , &
     maxiter )

  ! Compute interpolant spline for initial density
  call spline_interp_2d % compute_interpolant( sim_state % spline_2d_rho, sim_state % rho(:,1:ntau2) )

  ! Solve Poisson equation
  call poisson_solver % reset_charge()
  call poisson_solver % accumulate_charge( sim_state % spline_2d_rho )
  if ( sim_state % point_charges_present ) then
    do ic = 1, nc
      associate( intensity => sim_state % point_charges(ic)%intensity, &
                 location  => sim_state % point_charges(ic)%location )
        call poisson_solver % accumulate_charge( intensity, location )
      end associate
    end do
  end if
  call poisson_solver % solve( sim_state % spline_2d_phi )

  ! Initialize scalar diagnostics
  call diag % init( &
    ncells1         , &
    ncells2         , &
    p1              , &
    p2              , &
    nx1             , &
    nx2             , &
    tau_eta1        , &
    tau_eta2        , &
    breaks_eta1     , &
    breaks_eta2     , &
    mapping_discrete, &
    mapping_analytic, &
    sim_state )

  file_name = "scalar_diagnostics.dat"
  status    = 'replace'
  position  = 'asis'
  open( file=file_name, newunit=file_unit, action='write', status=status, position=position )

  ! Write scalar diagnostics
  call diag % write_scalar_data( file_unit, 0.0_wp )

  ! Write 2D data on interpolation grid
  call diag % write_on_interpolation_grid( file_id, 0 )

  ! Write electric field on Cartesian grid
  call diag % write_on_cartesian_grid( file_id, 0 )

  ! Write point charges trajectory
  if ( sim_state % point_charges_present ) call diag % write_point_charges( file_id, 0 )

  ! HDF5 I/O
  call sll_o_hdf5_ser_write_attribute( file_id, "/", "time_step", dt, h5_error )
  call sll_o_hdf5_ser_write_attribute( file_id, "/", "diag_freq", diag_freq, h5_error )

  !-----------------------------------------------------------------------------
  ! Time cycle
  !-----------------------------------------------------------------------------

  if ( iter > 0 ) write(*,'(/a/)') " >> Entering time cycle"

  ! Average time per iteration
  t_iter = 0.0_wp

  do it = 1, iter

    write(*,'(a,i4)',advance='no') "    iteration ", it

    call sll_s_set_time_mark( t0 )

    call time_integrator % advance_in_time( success )
    if ( .not. success ) exit

    ! Write scalar diagnostics
    call diag % write_scalar_data( file_unit, it*dt )

    ! Write 2D diagnostics
    if ( mod( it, diag_freq ) == 0 ) then

      ! Write point charges trajectory
      if ( sim_state % point_charges_present ) call diag % write_point_charges( file_id, it )

      ! Write 2D data on interpolation grid
      call diag % write_on_interpolation_grid( file_id, it )

      ! Write electric field on Cartesian grid
      call diag % write_on_cartesian_grid( file_id, it )

    end if

    call sll_s_set_time_mark( t1 )

    t_diff = sll_f_time_elapsed_between( t0, t1 )
    write(*,'(a,es7.1,a)') " ( ", t_diff, " s )"

    t_iter = t_iter + t_diff

  end do

  call sll_o_hdf5_ser_write_attribute( file_id, "/", "iterations", it, h5_error )

  if ( iter > 0 ) then
    ! Average time per iteration
    if ( it /= 1 ) t_iter = t_iter / real(it-1,wp)
    write(*,'(/a,es7.1,a)') " >> Average time per iteration: ", t_iter, " s"
  end if

  ! Close HDF5 file
  call sll_s_hdf5_ser_file_close( file_id, h5_error )

  close( file_unit )

  !-----------------------------------------------------------------------------
  ! Deallocations and free
  !-----------------------------------------------------------------------------

  ! Deallocate real 1D allocatables
  deallocate( breaks_eta1, breaks_eta2, tau_eta1, tau_eta2, intensity )

  ! Deallocate real 2D allocatables
  deallocate( phi, location )

  ! Free abstract polymorphic types
  call bsplines_eta1 % free()
  call bsplines_eta2 % free()

  ! Deallocate abstract polymorphic types
  deallocate( bsplines_eta1, bsplines_eta2, mapping_analytic )

  ! Free concrete types
  call mapping_discrete % free()
  call spline_2d_rho    % free()
  call spline_2d_phi    % free()
  call spline_interp_2d % free()
  call poisson_solver   % free()
  call electric_field   % free()
  call diag             % free()
  call sim_state        % free()
  call time_integrator  % free()

  write(*,'(/a/)') " >> End of simulation"

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
  SLL_PURE function rho_equilibrium( eta )
    real(wp), intent(in) :: eta(2)
    real(wp) :: rho_equilibrium

!    ! Diocotron instability
!    associate( smid => 0.5_wp * ( smin + smax ), dist => 0.5_wp * ( smax - smin ) )
!
!      ! smoothing on initial condition suggested in doi:10.1140/epjd/e2014-50180-9 (eq. 10)
!      rho_equilibrium = exp( - ( (eta(1)-smid)/dist )**p )
!
!    end associate

    if ( eta(1) <= 0.8_wp ) then
      rho_equilibrium = 1.0_wp - 1.25_wp * eta(1)
    else
      rho_equilibrium = 0.0_wp
    end if

  end function rho_equilibrium

  !-----------------------------------------------------------------------------
  SLL_PURE function rho_perturbation( eta )
    real(wp), intent(in) :: eta(2)
    real(wp) :: rho_perturbation

!    ! Diocotron instability
!    associate( smid => 0.5_wp * ( smin + smax ), dist => 0.5_wp * ( smax - smin ) )
!
!      ! smoothing on initial condition suggested in doi:10.1140/epjd/e2014-50180-9 (eq. 10)
!      rho_perturbation = ampl * cos( l*eta(2) ) * exp( - ( (eta(1)-smid)/dist )**p )
!
!    end associate

!    ! Vortex merger
!    real(wp) :: s, x(2), x0(2), rho_p1, rho_p2
!
!    s  = 0.08_wp
!
!    x  = mapping_discrete % eval( eta )
!
!    x0 = (/ 0.08_wp, -0.14_wp /)
!    rho_p1 = ampl * exp( - ( 0.5_wp*((x(1)-x0(1))**2+(x(2)-x0(2))**2)/s**2 ) )
!
!    x0 = (/ -0.08_wp, 0.14_wp /)
!    rho_p2 = ampl * exp( - ( 0.5_wp*((x(1)-x0(1))**2+(x(2)-x0(2))**2)/s**2 ) )
!
!    rho_perturbation = rho_p1 + rho_p2

    ! Point-like vortexes
    rho_perturbation = 0.0_wp

  end function rho_perturbation

end program sim_bsl_gc_2d0v_smooth_polar_splines
