program test_polar_mapping
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  use sll_m_working_precision, only: f64

  use sll_m_constants, only: sll_p_twopi

  use sll_m_utilities, only: sll_s_new_array_linspace

  use sll_m_polar_mapping_analytical_target, only: sll_t_polar_mapping_analytical_target

  use sll_m_polar_mapping_analytical_czarny, only: sll_t_polar_mapping_analytical_czarny

  use sll_m_polar_mapping_iga, only: sll_t_polar_mapping_iga

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_greville, &
    sll_p_periodic

  use sll_m_bsplines, only: &
    sll_c_bsplines, &
    sll_s_bsplines_new

  use sll_m_spline_interpolator_1d, only: sll_s_spline_1d_compute_num_cells 

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Working precision
  integer, parameter :: wp = f64

  ! To initialize B-splines
  integer :: n1, n2, deg1, deg2, ncells1, ncells2
  ! To plot mesh
  integer :: npts1, npts2

  real(wp), allocatable :: breaks_eta1(:)

  ! 1D B-splines
  class(sll_c_bsplines), allocatable :: spline_basis_eta1
  class(sll_c_bsplines), allocatable :: spline_basis_eta2

  ! Note: Circle corresponds to Target or Czarny default mapping
  type(sll_t_polar_mapping_analytical_target) :: mapping_analytical_circle
  type(sll_t_polar_mapping_analytical_target) :: mapping_analytical_target
  type(sll_t_polar_mapping_analytical_czarny) :: mapping_analytical_czarny
  type(sll_t_polar_mapping_iga) :: mapping_iga

  npts1 = 10
  npts2 = 40

  !-----------------------------------------------------------------------------
  ! Analytical mappings
  !-----------------------------------------------------------------------------

  ! Circle (Target or Czarny default mapping)
  call mapping_analytical_circle % init()

  call mapping_analytical_circle % store_data( npts1, npts2, "mapping_analytical_circle")

  ! Target
  call mapping_analytical_target % init( x0=[0.2_wp,-0.1_wp], d0=0.2_wp, e0=0.3_wp )

  call mapping_analytical_target % store_data( npts1, npts2, "mapping_analytical_target")

  ! Czarny
  call mapping_analytical_czarny % init( x0=[0.2_wp,-0.1_wp], b=1.4_wp, e=0.3_wp )

  call mapping_analytical_czarny % store_data( npts1, npts2, "mapping_analytical_czarny")

  !-----------------------------------------------------------------------------
  ! Discrete IGA mappings
  !-----------------------------------------------------------------------------
  n1 = 7
  n2 = 20

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

  ! Compute number of cells from number of interpolation points along eta2
  call sll_s_spline_1d_compute_num_cells( &
    degree  = deg2          , &
    bc_xmin = sll_p_periodic, &
    bc_xmax = sll_p_periodic, &
    nipts   = n2            , &
    ncells  = ncells2 )

  ! Create 1D spline basis along eta1 in [0,1]
  call sll_s_bsplines_new( &
    bsplines = spline_basis_eta1, &
    degree   = deg1             , &
    periodic = .false.          , &
    xmin     = 0.0_wp           , &
    xmax     = 1.0_wp           , &
    ncells   = ncells1          , &
    breaks   = breaks_eta1 )

  ! Create 1D spline basis along eta2 in [0,2pi]
  call sll_s_bsplines_new( &
    bsplines = spline_basis_eta2, &
    degree   = deg2             , &
    periodic = .true.           , &
    xmin     = 0.0_wp           , &
    xmax     = sll_p_twopi      , &
    ncells   = ncells2 )

  ! Circle (Target or Czarny default mapping)
  call mapping_iga % init( spline_basis_eta1, spline_basis_eta2, mapping_analytical_circle )

  call mapping_iga % store_data( npts1, npts2, "mapping_iga_circle")

  call mapping_iga % free()

  ! Target
  call mapping_iga % init( spline_basis_eta1, spline_basis_eta2, mapping_analytical_target )

  call mapping_iga % store_data( npts1, npts2, "mapping_iga_target")

  call mapping_iga % free()

  ! Czarny
  call mapping_iga % init( spline_basis_eta1, spline_basis_eta2, mapping_analytical_czarny )

  call mapping_iga % store_data( npts1, npts2, "mapping_iga_czarny")

  call mapping_iga % free()

end program test_polar_mapping
