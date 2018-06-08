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
    sll_s_hdf5_ser_file_close , &
    sll_o_hdf5_ser_write_array, &
    sll_o_hdf5_ser_write_attribute

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Working precision
  integer, parameter :: wp = f64

  ! To initialize B-splines (p1,p2 degrees)
  integer :: mm, n1, n2, p1, p2, ncells1, ncells2, ntau1, ntau2

  ! B-splines break points
  real(wp), allocatable :: breaks_eta1(:)
  real(wp), allocatable :: breaks_eta2(:)

  ! 1D B-splines
  class(sll_c_bsplines), allocatable :: bsplines_eta1
  class(sll_c_bsplines), allocatable :: bsplines_eta2

  ! Analytical and discrete mappings
  class(sll_c_polar_mapping_analytical), allocatable :: mapping_analytical
  type(sll_t_polar_mapping_iga) :: mapping_iga

  ! 2D splines for rho and phi
  type(sll_t_spline_2d) :: spline_2d_rho
  type(sll_t_spline_2d) :: spline_2d_phi

  ! 2D spline interpolator
  type(sll_t_spline_interpolator_2d) :: spline_interp_2d

  ! Needed for 2D interpolation of right hand side
  real(wp), allocatable :: tau_eta1(:)
  real(wp), allocatable :: tau_eta2(:)
  real(wp), allocatable :: rho(:,:), rho_new(:,:)

  ! Poisson solver
  type(sll_t_poisson_2d_fem_sps_stencil_new) :: poisson_solver

  ! Auxiliary variables
  integer  :: i1, i2
  integer  :: i, iter
  real(wp) :: dt, half_dt
  real(wp) :: eta(2)

  ! Timing
  type(sll_t_time_mark) :: t0, t1
  real(wp) :: t_diff, t_iter

  ! ODE integrator
  real(wp), parameter :: abs_tol = 1.0e-14_wp
  real(wp), parameter :: rel_tol = 1.0e-12_wp
  integer , parameter :: maxiter = 100
  logical :: success

  ! HDF5 I/O
  type(sll_t_hdf5_ser_handle) :: file_id
  integer                     :: h5_error
  character(len=32)           :: attr_name

  !-----------------------------------------------------------------------------
  ! Initialize B-splines basis functions
  !-----------------------------------------------------------------------------

  ! Number of degrees of freedom (control points) along s and theta
  mm = 128
  n1 = mm
  n2 = mm * 2

  ! Spline degrees along s and theta
  p1 = 3
  p2 = 3

  ! Create HDF5 file
  call sll_s_hdf5_ser_file_create( 'sim_bsl_gc_2d0v_smooth_polar_splines.h5', file_id, h5_error )

  ! HDF5 I/O
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

  !-----------------------------------------------------------------------------
  ! Initialize mapping
  !-----------------------------------------------------------------------------

  write(*,'(/a)') " >> Initializing mapping"

  allocate( sll_t_polar_mapping_analytical_target :: mapping_analytical )
!  allocate( sll_t_polar_mapping_analytical_czarny :: mapping_analytical )

  ! Analytical mapping
  select type ( mapping_analytical )
    type is ( sll_t_polar_mapping_analytical_target )
      call mapping_analytical % init( x0=[0.0_wp,0.0_wp], d0=0.0_wp, e0=0.0_wp ) ! circular mapping
!      call mapping_analytical % init( x0=[0.0_wp,0.0_wp], d0=0.2_wp, e0=0.3_wp )
    type is ( sll_t_polar_mapping_analytical_czarny )
      call mapping_analytical % init( x0=[0.0_wp,0.0_wp], b =1.4_wp, e =0.3_wp )
  end select

  ! Discrete mapping
  call mapping_iga % init( bsplines_eta1, bsplines_eta2, mapping_analytical )

  ! HDF5 I/O
  call mapping_iga % store_data( n1, n2, file_id )

  !-----------------------------------------------------------------------------
  ! Initialize splines and spline interpolator
  !-----------------------------------------------------------------------------

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

  !-----------------------------------------------------------------------------
  ! Initialize Poisson solver
  !-----------------------------------------------------------------------------

  write(*,'(/a)',advance='no') " >> Initializing Poisson solver"
 
  call sll_s_set_time_mark( t0 )

  ! Initialize Poisson solver
  call poisson_solver % init( bsplines_eta1, bsplines_eta2, breaks_eta1, breaks_eta2, mapping_iga )

  ! Set boundary conditions
  call poisson_solver % set_boundary_conditions( sll_p_dirichlet )

  call sll_s_set_time_mark( t1 )

  t_diff = sll_f_time_elapsed_between( t0, t1 )
  write(*,'(a,es7.1,a)') " ( ", t_diff, " s )"

  !-----------------------------------------------------------------------------
  ! Compute interpolators
  !-----------------------------------------------------------------------------

  allocate( rho    ( ntau1, ntau2+1 ) ) ! repeated point along theta
  allocate( rho_new( ntau1, ntau2+1 ) ) ! repeated point along theta

  ! Evaluate right hand side on interpolation points
  do i2 = 1, ntau2
    do i1 = 1, ntau1
      eta(1) = tau_eta1(i1)
      eta(2) = tau_eta2(i2)
      rho(i1,i2) = rho_initial( eta )
    end do
  end do

  ! Apply periodicity along theta
  rho(:,ntau2+1) = rho(:,1)

  ! HDF5 I/O
  write( attr_name, '(a,i0)' ) "/rho_", 0
  call sll_o_hdf5_ser_write_array( file_id, rho, trim(attr_name), h5_error )

  ! Compute interpolant spline for initial density
  call spline_interp_2d % compute_interpolant( spline_2d_rho, rho(:,1:ntau2) )

  !-----------------------------------------------------------------------------
  ! Time cycle
  !-----------------------------------------------------------------------------

  ! Solve Poisson equation
  call poisson_solver % solve( spline_2d_rho, spline_2d_phi )

  ! Time step
  dt      = 0.1_wp
  half_dt = 0.5_wp * dt

  ! HDF5 I/O
  call sll_o_hdf5_ser_write_attribute( file_id, "/", "time_step", dt, h5_error )

  write(*,'(/a/)') " >> Entering time cycle"

  ! Average time per iteration
  t_iter = 0.0_wp

  iter = 1000
  do i = 1, iter

    write(*,'(a,i4)',advance='no') "    iteration ", i

    call sll_s_set_time_mark( t0 )

    ! Predictor
    call advect_pseudo_cartesian( &
      -half_dt          , &
      tau_eta1          , &
      tau_eta2          , &
      mapping_analytical, &
      spline_2d_phi     , &
      spline_2d_rho     , &
      abs_tol           , &
      rel_tol           , &
      maxiter           , &
      success           , &
      rho_new )

    ! If integration of characteristics did not converge, exit time cycle, close HDF5 file and end simulation
    if ( .not. success ) exit

    call spline_interp_2d % compute_interpolant( spline_2d_rho, rho_new(:,1:ntau2) )

    ! Solve Poisson equation
    call poisson_solver % solve ( spline_2d_rho, spline_2d_phi )

    call spline_interp_2d % compute_interpolant( spline_2d_rho, rho(:,1:ntau2) )

    ! Corrector
    call advect_pseudo_cartesian( &
      -dt               , &
      tau_eta1          , &
      tau_eta2          , &
      mapping_analytical, &
      spline_2d_phi     , &
      spline_2d_rho     , &
      abs_tol           , &
      rel_tol           , &
      maxiter           , &
      success           , &
      rho_new )

    ! If integration of characteristics did not converge, exit time cycle, close HDF5 file and end simulation
    if ( .not. success ) exit

    call spline_interp_2d % compute_interpolant( spline_2d_rho, rho_new(:,1:ntau2) )

    ! Solve Poisson equation
    call poisson_solver % solve ( spline_2d_rho, spline_2d_phi )

    rho = rho_new

    ! HDF5 I/O
    write( attr_name, '(a,i0)' ) "/rho_", i
    call sll_o_hdf5_ser_write_array( file_id, rho, trim(attr_name), h5_error )

    call sll_s_set_time_mark( t1 )

    t_diff = sll_f_time_elapsed_between( t0, t1 )
    write(*,'(a,es7.1,a)') " ( ", t_diff, " s )"

    t_iter = t_iter + t_diff

  end do

  ! HDF5 I/O
  call sll_o_hdf5_ser_write_attribute( file_id, "/", "iterations", i, h5_error )

  ! Average time per iteration
  if ( i /= 1 ) t_iter = t_iter / real(i-1,wp)
  write(*,'(/a,es7.1,a)') " >> Average time per iteration: ", t_iter, " s"

  ! Close HDF5 file
  call sll_s_hdf5_ser_file_close ( file_id, h5_error )

  !-----------------------------------------------------------------------------
  ! Deallocations and free
  !-----------------------------------------------------------------------------

  deallocate( breaks_eta1 )
  deallocate( breaks_eta2 )

  call bsplines_eta1 % free()
  call bsplines_eta2 % free()

  deallocate( mapping_analytical )
  call mapping_iga % free()

  call spline_2d_rho % free()
  call spline_2d_phi % free()

  call spline_interp_2d % free()

  deallocate( tau_eta1 )
  deallocate( tau_eta2 )
  deallocate( rho     )
  deallocate( rho_new )

  call poisson_solver % free()

  write(*,'(/a/)') " >> End of simulation"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  SLL_PURE function logical_to_pseudo_cartesian( eta ) result( x )
    real(wp), intent(in) :: eta(2)
    real(wp) :: x(2)

    x(1) = eta(1) * cos( eta(2) )
    x(2) = eta(1) * sin( eta(2) )

  end function logical_to_pseudo_cartesian

  !-----------------------------------------------------------------------------
  SLL_PURE function pseudo_cartesian_to_logical( x ) result( eta )
    real(wp), intent(in) :: x(2)
    real(wp) :: eta(2)

    eta(1) = norm2( x )
    eta(2) = modulo( atan2( x(2), x(1) ), sll_p_twopi ) ! atan2 returns theta in [-pi,pi)

    ! For characteristics going out of the domain
    if ( eta(1) > 1.0_wp ) eta(1) = 1.0_wp

  end function pseudo_cartesian_to_logical

  !-----------------------------------------------------------------------------
  SLL_PURE function rho_initial( eta )
    real(wp), intent(in) :: eta(2)
    real(wp) :: rho_initial

    integer  :: l
    real(wp) :: smin, smax, eps

    ! TODO: move this as input (read from input file in simulation)
    l    = 6
    smin = 0.45_wp
    smax = 0.55_wp
    eps  = 1.0e-01_wp

    ! see doi:10.1140/epjd/e2014-50180-9, equation (10)
    if ( smin <= eta(1) .and. eta(1) <= smax ) then
      rho_initial = ( 1.0_wp + eps*cos( l*eta(2) ) ) * exp( -( eta(1) - 0.5_wp )**2 / 0.002_wp )
!      rho_initial = 1.0_wp + eps * cos( l*eta(2) )
    else
      rho_initial = 0.0_wp
    end if

  end function rho_initial

  !-----------------------------------------------------------------------------
  function electric_field( eta, mapping, spline_2d_phi )
    real(wp)                             , intent(in) :: eta(2)
    class(sll_c_polar_mapping_analytical), intent(in) :: mapping
    type(sll_t_spline_2d)                , intent(in) :: spline_2d_phi
    real(wp) :: electric_field(2)

    ! NOTE: this parameter can not be arbitrarily small
    real(wp), parameter :: eps = 1.0e-02_wp

    real(wp) :: jdet, jmat(2,2)
    real(wp) :: d1, d2, d3, d4, d5, d6
    real(wp) :: th1, th2
    real(wp) :: ef_0(2), ef_eps(2)

    if ( eta(1) == 0.0_wp ) then

      th1 = 0.1_wp * sll_p_pi
      th2 = 0.3_wp * sll_p_pi

      d1 = spline_2d_phi % eval_deriv_x1( eta(1), th1 ) ! dphi/ds(0,theta_1)
      d2 = spline_2d_phi % eval_deriv_x1( eta(1), th2 ) ! dphi/ds(0,theta_2)

!      write(*,*)
!      write(*,'(a)',advance='no') "    d1  = " ; write(*,*) d1
!      write(*,'(a)',advance='no') "    d2  = " ; write(*,*) d2

      jmat = mapping % jmat( (/ 0.0_wp, th1 /) )
 
      d3 = jmat(1,1) ! dx/ds(0,theta_1)
      d4 = jmat(2,1) ! dy/ds(0,theta_1)

!      write(*,'(a)',advance='no') "    d3  = " ; write(*,*) d3
!      write(*,'(a)',advance='no') "    d4  = " ; write(*,*) d4

      jmat = mapping % jmat( (/ 0.0_wp, th2 /) )
 
      d5 = jmat(1,1) ! dx/ds(0,theta_2)
      d6 = jmat(2,1) ! dy/ds(0,theta_2)

!      write(*,'(a)',advance='no') "    d5  = " ; write(*,*) d5
!      write(*,'(a)',advance='no') "    d6  = " ; write(*,*) d6

      electric_field(1) = - ( d4*d2 - d1*d6 ) / ( d4*d5 - d3*d6 )
      electric_field(2) = - ( d1*d5 - d3*d2 ) / ( d4*d5 - d3*d6 )

    else if ( 0.0_wp < eta(1) .and. eta(1) < eps ) then

      th1 = 0.1_wp * sll_p_pi
      th2 = 0.3_wp * sll_p_pi

      d1 = spline_2d_phi % eval_deriv_x1( eta(1), th1 ) ! dphi/ds(0,theta_1)
      d2 = spline_2d_phi % eval_deriv_x1( eta(1), th2 ) ! dphi/ds(0,theta_2)

      jmat = mapping % jmat( (/ 0.0_wp, th1 /) )
 
      d3 = jmat(1,1) ! dx/ds(0,theta_1)
      d4 = jmat(2,1) ! dy/ds(0,theta_1)

      jmat = mapping % jmat( (/ 0.0_wp, th2 /) )
 
      d5 = jmat(1,1) ! dx/ds(0,theta_2)
      d6 = jmat(2,1) ! dy/ds(0,theta_2)

      ! E(0,theta)
      ef_0(1) = - ( d4*d2 - d1*d6 ) / ( d4*d5 - d3*d6 )
      ef_0(2) = - ( d1*d5 - d3*d2 ) / ( d4*d5 - d3*d6 )

      jmat = mapping % jmat( eta )
      jdet = mapping % jdet( eta )

      ! E(eps,theta): J^(-T) times 'logical' gradient
      ef_eps(1) = - (   jmat(2,2) * spline_2d_phi % eval_deriv_x1( eta(1), eta(2) ) &
                      - jmat(2,1) * spline_2d_phi % eval_deriv_x2( eta(1), eta(2) ) ) / jdet
      ef_eps(2) = - ( - jmat(1,2) * spline_2d_phi % eval_deriv_x1( eta(1), eta(2) ) &
                      + jmat(1,1) * spline_2d_phi % eval_deriv_x2( eta(1), eta(2) ) ) / jdet

      ! Linear interpolation between 0 and eps
      electric_field(1) = (1.0_wp-eta(1)/eps)*ef_0(1) + eta(1)/eps*ef_eps(1)
      electric_field(2) = (1.0_wp-eta(1)/eps)*ef_0(2) + eta(1)/eps*ef_eps(2)

    else if ( eps <= eta(1) ) then

      jmat = mapping % jmat( eta )
      jdet = mapping % jdet( eta )

      ! J^(-T) times 'logical' gradient
      electric_field(1) = - (   jmat(2,2) * spline_2d_phi % eval_deriv_x1( eta(1), eta(2) ) &
                              - jmat(2,1) * spline_2d_phi % eval_deriv_x2( eta(1), eta(2) ) ) / jdet
      electric_field(2) = - ( - jmat(1,2) * spline_2d_phi % eval_deriv_x1( eta(1), eta(2) ) &
                              + jmat(1,1) * spline_2d_phi % eval_deriv_x2( eta(1), eta(2) ) ) / jdet

    end if

  end function electric_field

  !-----------------------------------------------------------------------------
  subroutine advect_pseudo_cartesian( &
    h            , &
    tau_eta1     , &
    tau_eta2     , &
    mapping      , &
    spline_2d_phi, &
    spline_2d_rho, &
    abs_tol      , &
    rel_tol      , &
    maxiter      , &
    success      , &
    rho_new )
    real(wp)                             , intent(in   ) :: h
    real(wp)                             , intent(in   ) :: tau_eta1(:)
    real(wp)                             , intent(in   ) :: tau_eta2(:)
    class(sll_c_polar_mapping_analytical), intent(in   ) :: mapping
    type(sll_t_spline_2d)                , intent(in   ) :: spline_2d_phi
    type(sll_t_spline_2d)                , intent(in   ) :: spline_2d_rho
    real(wp)                             , intent(in   ) :: abs_tol
    real(wp)                             , intent(in   ) :: rel_tol
    integer                              , intent(in   ) :: maxiter
    logical                              , intent(  out) :: success
    real(wp)                             , intent(inout) :: rho_new(:,:)

    integer  :: i1, i2, j, ntau1, ntau2
    real(wp) :: tol_sqr
    real(wp) :: x0(2), x(2), dx(2), temp(2), k2(2), a0(2), E(2), eta(2)
    real(wp) :: jmat_comp(2,2)

    character(len=*), parameter :: this_sub_name = "advect_pseudo_cartesian"
    character(len=256) :: err_msg

    ntau1 = size( tau_eta1 )
    ntau2 = size( tau_eta2 )

    success = .true.

!    write(*,*)

    do i2 = 1, ntau2
      do i1 = 1, ntau1

!        write(*,*)
!        write(*,'(a,i3,a,i3)') "    i1 = ", i1, "    i2 = ", i2

        eta(1) = tau_eta1(i1)
        eta(2) = tau_eta2(i2)

        !-----------------------------------------------------------------------
        ! Trapezoidal rule for integrating characteristics
        !-----------------------------------------------------------------------

        x0 = logical_to_pseudo_cartesian( eta )

!        write(*,*)
!        write(*,'(a)',advance='no') "    x0  = " ; write(*,*) x0
!        write(*,'(a)',advance='no') "    eta = " ; write(*,*) eta

        tol_sqr = ( abs_tol + rel_tol * norm2( x0 ) )**2

!        write(*,*)
!        write(*,'(a)',advance='no') "    tol = " ; write(*,*) tol_sqr

        ! Cartesian components of electric field
        E = electric_field( eta, mapping, spline_2d_phi )

!        write(*,*)
!        write(*,'(a)',advance='no') "    E   = " ; write(*,*) E

        jmat_comp = mapping % jmat_comp( eta )

        ! Pseudo-Cartesian components of advection field
        a0(1) = jmat_comp(1,1) * (-E(2)) + jmat_comp(1,2) * E(1)
        a0(2) = jmat_comp(2,1) * (-E(2)) + jmat_comp(2,2) * E(1)

!        write(*,*)
!        write(*,'(a)',advance='no') "    a0  = " ; write(*,*) a0

        ! First iteration
        dx  = h * a0
        x   = x0 + dx
        eta = pseudo_cartesian_to_logical( x )

!        write(*,*)
!        write(*,'(a)',advance='no') "    x   = " ; write(*,*) x
!        write(*,'(a)',advance='no') "    eta = " ; write(*,*) eta

        ! Successive iterations if first iteration did not converge
        if ( dot_product( dx, dx ) > tol_sqr ) then

          success = .false.

          associate( k1 => a0, h_half => 0.5_wp*h, dx_old => temp, error => temp )

            do j = 2, maxiter

              jmat_comp = mapping % jmat_comp( eta )

              ! Cartesian components of advection field
              E = electric_field( eta, mapping, spline_2d_phi )

!              write(*,*)
!              write(*,'(a)',advance='no') "    E   = " ; write(*,*) E

              ! k2 = f(t,x_{i-1})
              k2(1) = jmat_comp(1,1) * (-E(2)) + jmat_comp(1,2) * E(1)
              k2(2) = jmat_comp(2,1) * (-E(2)) + jmat_comp(2,2) * E(1)

              dx_old = dx
              dx     = h_half*(k1+k2)
              error  = dx_old - dx
              x      = x0 + dx
              eta    = pseudo_cartesian_to_logical( x )

!              write(*,*)
!              write(*,'(a)',advance='no') "    x   = " ; write(*,*) x
!              write(*,'(a)',advance='no') "    eta = " ; write(*,*) eta

              if ( dot_product( error, error ) <= tol_sqr ) then
                success = .true. ; exit
              end if

            end do

          end associate

        end if

        ! Check if integrator converged
        if ( success ) then

          rho_new(i1,i2) = spline_2d_rho % eval( eta(1), eta(2) )

        else

          write( err_msg, '(a,i0,a)' ) "integration of characteristics did not converge after ", maxiter, " iterations"
          SLL_WARNING( this_sub_name, err_msg )

          return

        end if

      end do
    end do

    ! Apply periodicity along theta
    rho_new(:,ntau2+1) = rho_new(:,1)

  end subroutine advect_pseudo_cartesian

end program sim_bsl_gc_2d0v_smooth_polar_splines
