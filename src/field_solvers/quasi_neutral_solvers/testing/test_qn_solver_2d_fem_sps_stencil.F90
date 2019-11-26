program test_qn_solver_2d_fem_sps_stencil
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

  use sll_m_singular_mapping_analytic, only: sll_c_singular_mapping_analytic

  use sll_m_singular_mapping_analytic_target, only: sll_t_singular_mapping_analytic_target

  use sll_m_singular_mapping_analytic_czarny, only: sll_t_singular_mapping_analytic_czarny

  use sll_m_singular_mapping_discrete, only: sll_t_singular_mapping_discrete

  use sll_m_qn_solver_2d_fem_sps_stencil_new, only: sll_t_qn_solver_2d_fem_sps_stencil_new

  use sll_m_boundary_condition_descriptors, only: sll_p_dirichlet

  use sll_m_gauss_legendre_integration, only: &
    sll_f_gauss_legendre_points, &
    sll_f_gauss_legendre_weights

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

  ! Tolerance valid for n1=128, n2=256 and p1=p2=3 in input file
  ! (needs to be adjusted in case of different input values)
  real(wp), parameter :: tolerance = 1e-6_wp

  ! Integer variables
  integer :: mm, n1, n2, p1, p2, ncells1, ncells2, ntau1, ntau2, Nk1, Nk2, Nq1, Nq2, &
             i1, i2, k1, k2, q1, q2, maptype, file_er(2), file_unit!, nn

  ! Character variables
  character(len=:), allocatable :: input_file

  ! Logical variables
  logical :: success

  ! Real variables
  real(wp) :: dt, L2_norm, Linf_norm, eta(2), x(2)

  ! Real 1D allocatables
  real(wp), allocatable :: breaks_eta1(:), breaks_eta2(:), tau_eta1(:), tau_eta2(:)

  ! Real 2D allocatables
  real(wp), allocatable :: gtau(:,:), coeffs1(:,:), coeffs2(:,:), phi_spl(:,:), err(:,:), &
                           quad_points_eta1(:,:), quad_points_eta2(:,:), quad_weights_eta1(:,:), quad_weights_eta2(:,:)

  ! Real 4D allocatables
  real(wp), allocatable :: volume(:,:,:,:)

  ! Abstract polymorphic types
  class(sll_c_bsplines)                 , allocatable :: bsplines_eta1, bsplines_eta2
  class(sll_c_singular_mapping_analytic), allocatable :: mapping_analytic

  ! Concrete types
  type(sll_t_singular_mapping_discrete)        :: mapping_discrete
  type(sll_t_spline_2d)                        :: spline_2d_rhs, spline_2d_phi
  type(sll_t_spline_interpolator_2d)           :: spline_interp_2d
  type(sll_t_qn_solver_2d_fem_sps_stencil_new) :: solver
  type(sll_t_time_mark)                        :: t0, t1
  type(sll_t_hdf5_ser_handle)                  :: file_id(2)

!  ! Stiffness/mass dense matrices and C1 projections
!  real(wp), allocatable :: A (:,:)
!  real(wp), allocatable :: M (:,:)
!  real(wp), allocatable :: Ap(:,:)
!  real(wp), allocatable :: Mp(:,:)

  ! Namelists in input file
  namelist /splines/ &
    n1, &
    n2, &
    p1, &
    p2

  namelist /geometry/ &
    maptype

  ! Read input file
  input_file = trim( "test_qn_solver_2d_fem_sps_stencil.nml" )
  open( file=input_file, status='old', action='read', newunit=file_unit )
  read( file_unit, splines  ); rewind( file_unit )
  read( file_unit, geometry );  close( file_unit )

  !-----------------------------------------------------------------------------
  ! Initialize B-splines basis functions
  !-----------------------------------------------------------------------------

  ! Create HDF5 file
  call sll_s_hdf5_ser_file_create( 'test_qn_solver_2d_fem_sps_stencil.h5', file_id(1), file_er(1) )

  ! HDF5 I/O
  call sll_o_hdf5_ser_write_attribute( file_id(1), "/", "n1", n1, file_er(1) )
  call sll_o_hdf5_ser_write_attribute( file_id(1), "/", "n2", n2, file_er(1) )
  call sll_o_hdf5_ser_write_attribute( file_id(1), "/", "p1", p1, file_er(1) )
  call sll_o_hdf5_ser_write_attribute( file_id(1), "/", "p2", p2, file_er(1) )

  call sll_s_set_time_mark( t0 )

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
        call mapping_analytic % init( x0=[0.0_wp,0.0_wp], Delta=0.0_wp, kappa=0.0_wp )
      else if ( maptype == 1 ) then
        call mapping_analytic % init( x0=[0.0_wp,0.0_wp], Delta=0.2_wp, kappa=0.3_wp )
      end if
    type is ( sll_t_singular_mapping_analytic_czarny )
      call mapping_analytic % init( x0=[0.0_wp,0.0_wp], e=1.4_wp, eps=0.3_wp )
  end select

  ! Initialize discrete mapping
  call mapping_discrete % init( bsplines_eta1, bsplines_eta2, mapping_analytic )

  ! Write mapping info
  call mapping_discrete % store_data( n1, n2, file_id(1) )

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

  ntau1 = size( tau_eta1 )
  ntau2 = size( tau_eta2 )

  ! Evaluate right hand side on interpolation points
  allocate( gtau( ntau1, ntau2 ) )

  call sll_s_hdf5_ser_file_open( 'test_qn_solver_2d_fem_sps_stencil_rho.h5', file_id(2), file_er(2) )
  call sll_o_hdf5_ser_read_array( file_id(2), gtau, '/rho_grid', file_er(2) )
  !call sll_s_hdf5_ser_file_close( file_id(2), file_er(2) )

  ! Compute interpolant spline
  call spline_interp_2d % compute_interpolant( spline_2d_rhs, gtau )

  call sll_s_set_time_mark( t1 )

  dt = sll_f_time_elapsed_between( t0, t1 )
  write(*,'(a,es8.1/)' ) " Time required for interpolation of right hand side: ", dt

  ! Assign 2D coefficients
  allocate( coeffs1( ntau1, ntau2 ) )
  allocate( coeffs2( ntau1, ntau2 ) )

  do i2 = 1, ntau2
    do i1 = 1, ntau1
      eta = (/ tau_eta1(i1), tau_eta2(i2) /)
      coeffs1(i1,i2) = fun_coeffs1( eta )
      coeffs2(i1,i2) = fun_coeffs2( eta )
    end do
  end do

  !-----------------------------------------------------------------------------
  ! Quasi-neutrality solver
  !-----------------------------------------------------------------------------

  ! Initialize 2D spline for solution
  call spline_2d_phi % init( bsplines_eta1, bsplines_eta2 )

  call sll_s_set_time_mark( t0 )

  ! Initialize quasi-neutrality solver
  call solver % init( &
    bsplines_eta1   , &
    bsplines_eta2   , &
    breaks_eta1     , &
    breaks_eta2     , &
    mapping_discrete, &
    coeffs1         , &
    coeffs2 )

  call sll_s_set_time_mark( t1 )

  dt = sll_f_time_elapsed_between( t0, t1 )
  write(*,'(a,es8.1/)' ) " Time required for initialization of quasi-neutrality solver: ", dt

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

  ! Set boundary conditions
  call solver % set_boundary_conditions( sll_p_dirichlet )

  call sll_s_set_time_mark( t0 )

  ! Solve
  call solver % reset_charge()
  ! TODO: test properly signature with callable function
  call solver % accumulate_charge( spline_2d_rhs ) ! Right hand side is 2D spline
  call solver % solve( spline_2d_phi )

  call sll_s_set_time_mark( t1 )

  dt = sll_f_time_elapsed_between( t0, t1 )
  write(*,'(a,es8.1/)' ) " Time required for solution of quasi-neutrality equation: ", dt

  ! Evaluate phi spline on logical grid and compute spatial L2-norm of error

  ! Repeated point in theta
  allocate( phi_spl( ntau1, ntau2+1 ) )
  allocate( err( ntau1, ntau2+1 ) )

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

  ! Compute L2-norm error
  L2_norm = 0.0_wp

  do k2 = 1, Nk2
    do k1 = 1, Nk1
      do q2 = 1, Nq2
        do q1 = 1, Nq1

          eta = (/ quad_points_eta1(q1,k1), quad_points_eta2(q2,k2) /)
          x   = mapping_discrete % eval( eta )

          L2_norm = L2_norm + volume(q1,q2,k1,k2) &
                    * ( spline_2d_phi % eval( eta(1), eta(2) ) - phi_exact( eta, x ) )**2

        end do
      end do
    end do
  end do

  L2_norm = sqrt( L2_norm )

  ! Compute Linf-norm of error
  Linf_norm = 0.0_wp

  do i2 = 1, ntau2
    do i1 = 1, ntau1

      eta = (/ tau_eta1(i1), tau_eta2(i2) /)
      x   = mapping_discrete % eval( eta )

      phi_spl(i1,i2) = spline_2d_phi % eval( eta(1), eta(2) )
      phi_spl(:,ntau2+1) = phi_spl(:,1)

      err(i1,i2) = phi_spl(i1,i2) - phi_exact( eta, x )
      err(:,ntau2+1) = err(:,1)

      Linf_norm = merge( abs( err(i1,i2) ), Linf_norm, abs( err(i1,i2) ) > Linf_norm )

    end do
  end do

  !-----------------------------------------------------------------------------
  ! HDF5 I/O
  !-----------------------------------------------------------------------------

  call sll_s_set_time_mark( t0 )

  ! Write solution
  call sll_o_hdf5_ser_write_array( file_id(1), solver % x, "/x", file_er(1) )

  ! Write reshaped solution
  call sll_o_hdf5_ser_write_array( file_id(1), spline_2d_phi % bcoef(1:n1,1:n2), "/phi", file_er(1) )

  ! Write phi spline
  call sll_o_hdf5_ser_write_array( file_id(1), phi_spl, "/phi_spl", file_er(1) )

  ! Write error
  call sll_o_hdf5_ser_write_array( file_id(1), err, "/error", file_er(1) )

!  ! Write stiffness matrix
!  call sll_o_hdf5_ser_write_array( file_id, A, "/A", file_er(1) )
!
!  ! Write mass matrix
!  call sll_o_hdf5_ser_write_array( file_id, M, "/M", file_er(1) )
!
!!  ! Write right hand side
!!  call sll_o_hdf5_ser_write_array( file_id, solver % b, "/b", file_er(1) )
!
!  ! Write C1 projection of stiffness matrix
!  call sll_o_hdf5_ser_write_array( file_id, Ap, "/Ap", file_er(1) )
!
!  ! Write C1 projection of mass matrix
!  call sll_o_hdf5_ser_write_array( file_id, Mp, "/Mp", file_er(1) )
!
!!  ! Write L matrix needed for projection
!!  call sll_o_hdf5_ser_write_array( file_id, solver % L, "/L", file_er(1) )

  ! Close HDF5 file
  call sll_s_hdf5_ser_file_close ( file_id(1), file_er(1) )

  call sll_s_set_time_mark( t1 )

  dt = sll_f_time_elapsed_between( t0, t1 )

  write(*,'(a,es8.1/)') " Time required for writing HDF5 output: ", dt
  write(*,'(a,es8.2/)') "   L2-norm of error: ", L2_norm
  write(*,'(a,es8.2/)') " Linf-norm of error: ", Linf_norm

  ! Check if test passed
  success = .true.
  if( max( L2_norm, Linf_norm ) > tolerance ) then
    success = .false.
    write(*,'(a,es7.1,a/)') "Test FAILED (tolerance ", tolerance, ")"
  end if
  if( success ) write(*,'(a/)') "Test PASSED"

  !-----------------------------------------------------------------------------
  ! Deallocations and free
  !-----------------------------------------------------------------------------

  deallocate( breaks_eta1, breaks_eta2, tau_eta1, tau_eta2, gtau, coeffs1, coeffs2, phi_spl, err, &
              quad_points_eta1, quad_points_eta2, quad_weights_eta1, quad_weights_eta2, volume )

!  deallocate( A, M, Ap, Mp )

  deallocate( mapping_analytic )

  call mapping_discrete % free()
  call bsplines_eta1 % free()
  call bsplines_eta2 % free()
  call spline_2d_rhs % free()
  call spline_2d_phi % free()
  call spline_interp_2d % free()
  call solver % free()

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  SLL_PURE function fun_coeffs1( eta )
    real(wp), intent(in) :: eta(2)
    real(wp) :: fun_coeffs1

    associate( s => eta(1), t => eta(2) )

      fun_coeffs1 = exp( - tanh( ( s - 0.5_wp ) / 0.1_wp ) )

    end associate

  end function fun_coeffs1

  SLL_PURE function fun_coeffs2( eta )
    real(wp), intent(in) :: eta(2)
    real(wp) :: fun_coeffs2

    associate( s => eta(1), t => eta(2) )

      fun_coeffs2 = 1.0_wp / exp( - tanh( ( s - 0.5_wp ) / 0.2_wp ) )

    end associate

  end function fun_coeffs2

  SLL_PURE function phi_exact( eta, x )
    real(wp), intent(in) :: eta(2)
    real(wp), intent(in) :: x(2)
    real(wp) :: phi_exact

    phi_exact = ( 1.0_wp - eta(1)**2 ) * cos( sll_p_twopi * x(1) ) * sin( sll_p_twopi * x(2) )

  end function phi_exact

end program test_qn_solver_2d_fem_sps_stencil
