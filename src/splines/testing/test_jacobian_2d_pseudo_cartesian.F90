program test_jacobian_2d_pseudo_cartesian
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"

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

  use sll_m_spline_2d, only: sll_t_spline_2d

  use sll_m_spline_interpolator_1d, only: &
    sll_s_spline_1d_compute_num_cells, &
    sll_t_spline_interpolator_1d

  use sll_m_spline_interpolator_2d, only: sll_t_spline_interpolator_2d

  use sll_m_polar_mapping_analytical, only: sll_c_polar_mapping_analytical

  use sll_m_polar_mapping_analytical_target, only: sll_t_polar_mapping_analytical_target

  use sll_m_polar_mapping_analytical_czarny, only: sll_t_polar_mapping_analytical_czarny

  use sll_m_polar_mapping_iga, only: sll_t_polar_mapping_iga

  use sll_m_jacobian_2d_pseudo_cartesian, only: sll_t_jacobian_2d_pseudo_cartesian

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Working precision
  integer, parameter :: wp = f64

  ! To initialize B-splines
  integer :: n1, n2, deg1, deg2, ncells1, ncells2

  ! To initialize non-uniform 1D B-splines along eta1
  real(wp), allocatable :: breaks_eta1(:)

  ! 1D B-splines
  class(sll_c_bsplines), allocatable :: spline_basis_eta1
  class(sll_c_bsplines), allocatable :: spline_basis_eta2

  ! 1D/2D tensor-product spline interpolators
  type(sll_t_spline_interpolator_2d) :: spline_interpolator_2d

  ! Interpolation points and fields
  integer :: npts1, npts2
  real(wp), allocatable :: e1_node(:)
  real(wp), allocatable :: e2_node(:)

  ! Auxiliary variables
  integer  :: i2, mm
  real(wp) :: eta(2)
  real(wp) :: jmat(2,2), jmat_analytical(2,2)
  real(wp) :: maxerr, errmat(2,2)
  real(wp), parameter :: kappa = 0.3_wp
  real(wp), parameter :: delta = 0.2_wp
  real(wp), parameter :: epsil = 0.3_wp
  real(wp), parameter :: ellip = 1.4_wp

  ! Analytical and discrete mappings
  class(sll_c_polar_mapping_analytical), allocatable :: mapping_analytical
  type(sll_t_polar_mapping_iga) :: mapping_discrete

  ! Quasi-Cartesian Jacobian
  type(sll_t_jacobian_2d_pseudo_cartesian) :: jacobian

  mm = 16
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

  ! Initialize 1D/2D tensor-product spline interpolators
  call spline_interpolator_2d % init( spline_basis_eta1, spline_basis_eta2, &
                                      [ sll_p_greville, sll_p_periodic ], &
                                      [ sll_p_greville, sll_p_periodic ] )

  ! Get interpolation points and allocate 1D/2D array of values
  call spline_interpolator_2d % get_interp_points( e1_node, e2_node )

  npts1 = size( e1_node )
  npts2 = size( e2_node )

  allocate( sll_t_polar_mapping_analytical_target :: mapping_analytical )

  ! Initialize analytical mapping
  select type ( mapping_analytical )
    type is ( sll_t_polar_mapping_analytical_target )
      call mapping_analytical % init( x0=[0.0_wp,0.0_wp], d0=delta, e0=kappa )
    type is ( sll_t_polar_mapping_analytical_czarny )
      call mapping_analytical % init( x0=[0.0_wp,0.0_wp], b =ellip, e =epsil )
  end select

  ! Initialize discrete mapping
  call mapping_discrete % init( spline_basis_eta1, spline_basis_eta2, mapping_analytical )

  ! Initialize composite inverse Jacobian
  call jacobian % init( mapping_discrete )

  ! Compute analytical composite inverse Jacobian
  select type ( mapping_analytical )
    type is ( sll_t_polar_mapping_analytical_target )
      jmat_analytical(1,1) = 1.0_wp / ( 1.0_wp - kappa )
      jmat_analytical(1,2) = 0.0_wp
      jmat_analytical(2,1) = 0.0_wp
      jmat_analytical(2,2) = 1.0_wp / ( 1.0_wp + kappa )
    type is ( sll_t_polar_mapping_analytical_czarny )
      jmat_analytical(1,1) = - sqrt( 1.0_wp + epsil**2 )
      jmat_analytical(1,2) = 0.0_wp
      jmat_analytical(2,1) = 0.0_wp
      jmat_analytical(2,2) = ( 2.0_wp - sqrt( 1.0_wp + epsil**2 ) ) / &
                             ( ellip /  sqrt( 1.0_wp - epsil**2 * 0.25_wp ) )
  end select

  write(*,*)

  maxerr = 0.0_wp

  do i2 = 1, npts2

      eta(1) = 0.0_wp
      eta(2) = e2_node(i2)

      jmat   = jacobian % eval( eta )
      errmat = abs( jmat - jmat_analytical )

      write(*,'(a,f18.15)') 'theta  = ', eta(2)

      write(*,'(a,f18.15,a,es8.2)') 'J(1,1) = ', jmat(1,1), '  error = ', errmat(1,1)
      maxerr = merge( errmat(1,1), maxerr, errmat(1,1) > maxerr )

      write(*,'(a,f18.15,a,es8.2)') 'J(2,2) = ', jmat(2,2), '  error = ', errmat(2,2)
      maxerr = merge( errmat(2,2), maxerr, errmat(2,2) > maxerr )

      write(*,'(a,f18.15,a,es8.2)') 'J(1,2) = ', jmat(1,2), '  error = ', errmat(1,2)
      maxerr = merge( errmat(1,2), maxerr, errmat(1,2) > maxerr )

      write(*,'(a,f18.15,a,es8.2)') 'J(2,1) = ', jmat(2,1), '  error = ', errmat(2,1)
      maxerr = merge( errmat(2,1), maxerr, errmat(2,1) > maxerr )

      write(*,*)

  end do

  write(*,'(a,es8.2)') 'maximum error = ', maxerr
  write(*,*)

  deallocate( e1_node )
  deallocate( e2_node )

  call spline_interpolator_2d % free()
  call spline_basis_eta1 % free()
  call spline_basis_eta2 % free()
  call mapping_discrete % free()
  call jacobian % free()

end program test_jacobian_2d_pseudo_cartesian
