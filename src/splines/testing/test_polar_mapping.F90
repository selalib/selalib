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

  use sll_m_spline_2d, only: sll_t_spline_2d

  use sll_m_spline_interpolator_1d, only: sll_s_spline_1d_compute_num_cells 

  use sll_m_spline_interpolator_2d, only: sll_t_spline_interpolator_2d

  use sll_m_polar_bsplines_2d, only: sll_t_polar_bsplines_2d

  use m_analytical_profiles_2d_cos_cos, only: t_analytical_profile_2d_cos_cos

  use sll_m_hdf5_io_serial, only: &
    sll_t_hdf5_ser_handle     , &
    sll_s_hdf5_ser_file_create, &
    sll_s_hdf5_ser_file_close , &
    sll_o_hdf5_ser_write_array

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Working precision
  integer, parameter :: wp = f64

  ! To initialize B-splines
  integer :: n1, n2, deg1, deg2, ncells1, ncells2

  ! To plot mesh
  integer :: npts1, npts2

  ! To initialize mappings
  real(wp) :: x0(2), d0, e0, b, e

  ! To initialize non-uniform 1D B-splines along eta1
  real(wp), allocatable :: breaks_eta1(:)

  ! 1D B-splines
  class(sll_c_bsplines), allocatable :: spline_basis_eta1
  class(sll_c_bsplines), allocatable :: spline_basis_eta2

  ! Note: Circle corresponds to Target or Czarny default mapping
  type(sll_t_polar_mapping_analytical_target) :: mapping_analytical_circle
  type(sll_t_polar_mapping_analytical_target) :: mapping_analytical_target
  type(sll_t_polar_mapping_analytical_czarny) :: mapping_analytical_czarny
  type(sll_t_polar_mapping_iga) :: mapping_iga

  ! New polar B-splines
  type(sll_t_polar_bsplines_2d) :: polar_spline_basis

  integer :: i1, i2

  real(wp), allocatable :: y1_grid(:), y2_grid(:)     ! intermediate coordinates
  real(wp), allocatable :: e1_grid(:,:), e2_grid(:,:) ! logical coordinates
  real(wp), allocatable :: x1_grid(:,:), x2_grid(:,:) ! physical coordinates
  real(wp), allocatable :: jacobian(:,:)

  real(wp) :: eta(2), x(2)

  ! For hdf5 I/O
  type(sll_t_hdf5_ser_handle) :: file_id
  integer                     :: h5_error

  !-----------------------------------------------------------------------------
  ! For interpolation tests
  !-----------------------------------------------------------------------------

  ! 2D tensor-product spline
  type(sll_t_spline_2d) :: spline_2d

  ! 2D tensor-product spline interpolator
  type(sll_t_spline_interpolator_2d) :: spline_interpolator_2d

  ! Interpolation points and profile values
  real(wp), allocatable :: tau1(:)
  real(wp), allocatable :: tau2(:)
  real(wp), allocatable :: gtau(:,:)

  ! Profile and profile parameters
  type(t_analytical_profile_2d_cos_cos) :: profile
  integer  :: m1, m2
  real(wp) :: c1, c2

  real(wp) :: error, max_error
  real(wp), allocatable :: interp_funct(:,:), interp_error(:,:)

  !=============================================================================
  ! Initialize B-splines basis functions
  !=============================================================================

  n1 = 40
  n2 = 80

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

  ! Number of points on test grids
  npts1 = 10
  npts2 = 40

  !-----------------------------------------------------------------------------
  ! Circle (Target or Czarny default mapping)
  !-----------------------------------------------------------------------------

  call mapping_analytical_circle % init()

  call mapping_analytical_circle % store_data( npts1, npts2, "mapping_analytic_circle")

  call mapping_iga % init( spline_basis_eta1, spline_basis_eta2, mapping_analytical_circle )

  call mapping_iga % store_data( npts1, npts2, "mapping_discrete_circle")

  call polar_spline_basis % init( spline_basis_eta1, spline_basis_eta2, mapping_iga )

  call polar_spline_basis % free()

  call mapping_iga % free()

  !-----------------------------------------------------------------------------
  ! Target
  !-----------------------------------------------------------------------------

  x0 = [0.08_wp,0.0_wp]
  d0 = 0.2_wp
  e0 = 0.3_wp

  call mapping_analytical_target % init( x0=x0, d0=d0, e0=e0 )

  call mapping_analytical_target % store_data( npts1, npts2, "mapping_analytic_target")

  call mapping_iga % init( spline_basis_eta1, spline_basis_eta2, mapping_analytical_target )

  call mapping_iga % store_data( npts1, npts2, "mapping_discrete_target")

  call polar_spline_basis % init( spline_basis_eta1, spline_basis_eta2, mapping_iga )

  call polar_spline_basis % free()

  call mapping_iga % free()

  !-----------------------------------------------------------------------------
  ! Czarny
  !-----------------------------------------------------------------------------

  x0 = [0.08_wp,0.0_wp]
  b  = 1.4_wp
  e  = 0.3_wp

  call mapping_analytical_czarny % init( x0=x0, b=b, e=e )

  call mapping_analytical_czarny % store_data( npts1, npts2, "mapping_analytic_czarny")

  call mapping_iga % init( spline_basis_eta1, spline_basis_eta2, mapping_analytical_czarny )

  call mapping_iga % store_data( npts1, npts2, "mapping_discrete_czarny")

  call polar_spline_basis % init( spline_basis_eta1, spline_basis_eta2, mapping_iga )

  call polar_spline_basis % free()

  call mapping_iga % free()

  !=============================================================================
  ! Test intermediate mapping
  !=============================================================================

  npts1 = 11
  npts2 = 11

  allocate( y1_grid( npts1 ) )
  allocate( y2_grid( npts2 ) )

  ! Create uniform grids in X and Y: X=s*cos(t), Y=s*sin(t)
  call sll_s_new_array_linspace( y1_grid, -1.0_wp, 1.0_wp, endpoint=.true. )
  call sll_s_new_array_linspace( y2_grid, -1.0_wp, 1.0_wp, endpoint=.true. )

  allocate( e1_grid( npts1, npts2 ) )
  allocate( e2_grid( npts1, npts2 ) )

  ! Create corresponding grids in s and t (theta)
  do i2 = 1, npts2
    do i1 = 1, npts1
      e1_grid(i1,i2) = sqrt ( y1_grid(i1)**2 + y2_grid(i2)**2 )
      e2_grid(i1,i2) = atan2( y2_grid(i2), y1_grid(i1) )
    end do
  end do

  allocate( x1_grid ( npts1, npts2 ) )
  allocate( x2_grid ( npts1, npts2 ) )
  allocate( jacobian( npts1, npts2 ) )

  ! Create corresponding grids in x and y plus Jacobian determinant (circle)

  call sll_s_hdf5_ser_file_create( 'mapping_intermed_circle.h5', file_id, h5_error )

  do i2 = 1, npts2
    do i1 = 1, npts1
      eta(1) = e1_grid(i1,i2)
      eta(2) = e2_grid(i1,i2)
      x(:) = mapping_analytical_circle % eval( eta )
      x1_grid(i1,i2) = x(1)
      x2_grid(i1,i2) = x(2)
      jacobian(i1,i2) = mapping_analytical_circle % jdet( eta )
    end do
  end do

  ! Store mesh
  call sll_o_hdf5_ser_write_array( file_id, x1_grid, "/x1", h5_error )
  call sll_o_hdf5_ser_write_array( file_id, x2_grid, "/x2", h5_error )
  call sll_o_hdf5_ser_write_array( file_id, jacobian, "/jacobian", h5_error )

  call sll_s_hdf5_ser_file_close ( file_id, h5_error )

  ! Create corresponding grids in x and y plus Jacobian determinant (Target)

  call sll_s_hdf5_ser_file_create( 'mapping_intermed_target.h5', file_id, h5_error )

  do i2 = 1, npts2
    do i1 = 1, npts1
      eta(1) = e1_grid(i1,i2)
      eta(2) = e2_grid(i1,i2)
      x(:) = mapping_analytical_target % eval( eta )
      x1_grid(i1,i2) = x(1)
      x2_grid(i1,i2) = x(2)
      jacobian(i1,i2) = mapping_analytical_target % jdet( eta )
    end do
  end do

  ! Store mesh
  call sll_o_hdf5_ser_write_array( file_id, x1_grid, "/x1", h5_error )
  call sll_o_hdf5_ser_write_array( file_id, x2_grid, "/x2", h5_error )
  call sll_o_hdf5_ser_write_array( file_id, jacobian, "/jacobian", h5_error )

  call sll_s_hdf5_ser_file_close ( file_id, h5_error )

  ! Create corresponding grids in x and y plus Jacobian determinant (Czarny)

  call sll_s_hdf5_ser_file_create( 'mapping_intermed_czarny.h5', file_id, h5_error )

  do i2 = 1, npts2
    do i1 = 1, npts1
      eta(1) = e1_grid(i1,i2)
      eta(2) = e2_grid(i1,i2)
      x(:) = mapping_analytical_czarny % eval( eta )
      x1_grid(i1,i2) = x(1)
      x2_grid(i1,i2) = x(2)
      jacobian(i1,i2) = mapping_analytical_czarny % jdet( eta )
    end do
  end do

  ! Store mesh
  call sll_o_hdf5_ser_write_array( file_id, x1_grid, "/x1", h5_error )
  call sll_o_hdf5_ser_write_array( file_id, x2_grid, "/x2", h5_error )
  call sll_o_hdf5_ser_write_array( file_id, jacobian, "/jacobian", h5_error )

  call sll_s_hdf5_ser_file_close ( file_id, h5_error )

  deallocate( y1_grid  )
  deallocate( y2_grid  )
  deallocate( e1_grid  )
  deallocate( e2_grid  )
  deallocate( x1_grid  )
  deallocate( x2_grid  )
  deallocate( jacobian )

  !=============================================================================
  ! Test standard 2D interpolation
  !=============================================================================

  ! Number of points on test grids
  npts1 = 10
  npts2 = 40

  ! Initialize 2D profile
  m1 = 2
  m2 = 3
  c1 = 0.0_wp
  c2 = 0.0_wp

  call profile % init( m1, m2, c1, c2 )

  ! Initialize 2D tensor-product spline
  call spline_2d % init( spline_basis_eta1, spline_basis_eta2 )

  ! Initialize 2D tensor-product spline interpolator
  call spline_interpolator_2d % init( spline_basis_eta1, spline_basis_eta2, &
                                      [ sll_p_greville, sll_p_periodic ], &
                                      [ sll_p_greville, sll_p_periodic ] )

  ! Get interpolation points and allocate 2D array of values
  call spline_interpolator_2d % get_interp_points( tau1, tau2 )
  allocate( gtau(size(tau1),size(tau2)) )

  !-----------------------------------------------------------------------------
  ! Circle
  !-----------------------------------------------------------------------------

  call mapping_analytical_circle % init()

  call mapping_iga % init( spline_basis_eta1, spline_basis_eta2, mapping_analytical_circle )

  ! Evaluate analytical profile at interpolation points
  do i2 = 1, n2
    do i1 = 1, n1
      x = mapping_iga % eval( [tau1(i1),tau2(i2)] )
      gtau(i1,i2) = profile % eval( x(1), x(2) )
    end do
  end do

  ! Compute spline interpolating profile on grid
  call spline_interpolator_2d % compute_interpolant( spline_2d, gtau )

  write(*,*)
  write(*,'(a)') " =============="
  write(*,'(a)') " CIRCLE mapping"
  write(*,'(a)') " =============="

  write(*,*)
  write(*,'(a)') " ***************************************************************************"
  write(*,'(a)') " TEST #1: evaluate standard 2D tensor product spline at interpolation points"
  write(*,'(a)') " ***************************************************************************"

  ! Evaluate spline at interpolation points
  max_error = 0.0_wp
  do i2 = 1, n2
    do i1 = 1, n1
      error = gtau(i1,i2) - spline_2d % eval( tau1(i1), tau2(i2) )
      max_error = max( max_error, abs( error ) )
    end do
  end do

  write(*,*)
  write(*,'(a,es8.2)') " Maximum error = ", max_error

  write(*,*)
  write(*,'(a)') " ***************************************************************************"
  write(*,'(a)') " TEST #2: evaluate standard 2D tensor product spline on test grid           "
  write(*,'(a)') " ***************************************************************************"

  ! Allocate grid for interpolation error
  allocate( interp_funct( npts1, npts2+1 ) )
  allocate( interp_error( npts1, npts2+1 ) )

  call sll_s_hdf5_ser_file_create( 'mapping_interpol_circle.h5', file_id, h5_error )

  ! Evaluate spline at interpolation points
  max_error = 0.0_wp
  do i2 = 1, npts2+1
    do i1 = 1, npts1
      eta(1) = real( i1-1, wp ) / real( npts1-1, wp )
      eta(2) = real( i2-1, wp ) * sll_p_twopi / real( npts2, wp )
      x  (:) = mapping_iga % eval( eta )
      interp_funct(i1,i2) = profile % eval( x(1), x(2) )
      interp_error(i1,i2) = profile % eval( x(1), x(2) ) - spline_2d % eval( eta(1), eta(2) )
      max_error = max( max_error, abs( interp_error(i1,i2) ) )
    end do
  end do

  ! Store mesh
  call sll_o_hdf5_ser_write_array( file_id, interp_funct, "/interp_funct", h5_error )
  call sll_o_hdf5_ser_write_array( file_id, abs( interp_error ), "/interp_error", h5_error )

  call sll_s_hdf5_ser_file_close ( file_id, h5_error )

  write(*,*)
  write(*,'(a,es8.2)') " Maximum error = ", max_error
  write(*,*)

  call mapping_iga % free()

  !-----------------------------------------------------------------------------
  ! Target
  !-----------------------------------------------------------------------------

  ! Initialize mapping
  x0 = [0.08_wp,0.0_wp]
  d0 = 0.2_wp
  e0 = 0.3_wp

  call mapping_analytical_target % init( x0=x0, d0=d0, e0=e0 )

  call mapping_iga % init( spline_basis_eta1, spline_basis_eta2, mapping_analytical_target )

  ! Evaluate analytical profile at interpolation points
  do i2 = 1, n2
    do i1 = 1, n1
      x = mapping_iga % eval( [tau1(i1),tau2(i2)] )
      gtau(i1,i2) = profile % eval( x(1), x(2) )
    end do
  end do

  ! Compute spline interpolating profile on grid
  call spline_interpolator_2d % compute_interpolant( spline_2d, gtau )

  write(*,*)
  write(*,'(a)') " =============="
  write(*,'(a)') " TARGET mapping"
  write(*,'(a)') " =============="

  write(*,*)
  write(*,'(a)') " ***************************************************************************"
  write(*,'(a)') " TEST #1: evaluate standard 2D tensor product spline at interpolation points"
  write(*,'(a)') " ***************************************************************************"

  ! Evaluate spline at interpolation points
  max_error = 0.0_wp
  do i2 = 1, n2
    do i1 = 1, n1
      error = gtau(i1,i2) - spline_2d % eval( tau1(i1), tau2(i2) )
      max_error = max( max_error, abs( error ) )
    end do
  end do

  write(*,*)
  write(*,'(a,es8.2)') " Maximum error = ", max_error

  write(*,*)
  write(*,'(a)') " ***************************************************************************"
  write(*,'(a)') " TEST #2: evaluate standard 2D tensor product spline on test grid           "
  write(*,'(a)') " ***************************************************************************"

  call sll_s_hdf5_ser_file_create( 'mapping_interpol_target.h5', file_id, h5_error )

  ! Evaluate spline at interpolation points
  max_error = 0.0_wp
  do i2 = 1, npts2+1
    do i1 = 1, npts1
      eta(1) = real( i1-1, wp ) / real( npts1-1, wp )
      eta(2) = real( i2-1, wp ) * sll_p_twopi / real( npts2, wp )
      x  (:) = mapping_iga % eval( eta )
      interp_funct(i1,i2) = profile % eval( x(1), x(2) )
      interp_error(i1,i2) = profile % eval( x(1), x(2) ) - spline_2d % eval( eta(1), eta(2) )
      max_error = max( max_error, abs( interp_error(i1,i2) ) )
    end do
  end do

  ! Store mesh
  call sll_o_hdf5_ser_write_array( file_id, interp_funct, "/interp_funct", h5_error )
  call sll_o_hdf5_ser_write_array( file_id, abs( interp_error ), "/interp_error", h5_error )

  call sll_s_hdf5_ser_file_close ( file_id, h5_error )

  write(*,*)
  write(*,'(a,es8.2)') " Maximum error = ", max_error
  write(*,*)

  call mapping_iga % free()

  !-----------------------------------------------------------------------------
  ! Czarny
  !-----------------------------------------------------------------------------

  ! Initialize mapping
  x0 = [0.08_wp,0.0_wp]
  b  = 1.4_wp
  e  = 0.3_wp

  call mapping_analytical_czarny % init( x0=x0, b=b, e=e )

  call mapping_iga % init( spline_basis_eta1, spline_basis_eta2, mapping_analytical_czarny )

  ! Evaluate analytical profile at interpolation points
  do i2 = 1, n2
    do i1 = 1, n1
      x = mapping_iga % eval( [tau1(i1),tau2(i2)] )
      gtau(i1,i2) = profile % eval( x(1), x(2) )
    end do
  end do

  ! Compute spline interpolating profile on grid
  call spline_interpolator_2d % compute_interpolant( spline_2d, gtau )

  write(*,*)
  write(*,'(a)') " =============="
  write(*,'(a)') " CZARNY mapping"
  write(*,'(a)') " =============="

  write(*,*)
  write(*,'(a)') " ***************************************************************************"
  write(*,'(a)') " TEST #1: evaluate standard 2D tensor product spline at interpolation points"
  write(*,'(a)') " ***************************************************************************"

  ! Evaluate spline at interpolation points
  max_error = 0.0_wp
  do i2 = 1, n2
    do i1 = 1, n1
      error = gtau(i1,i2) - spline_2d % eval( tau1(i1), tau2(i2) )
      max_error = max( max_error, abs( error ) )
    end do
  end do

  write(*,*)
  write(*,'(a,es8.2)') " Maximum error = ", max_error

  write(*,*)
  write(*,'(a)') " ***************************************************************************"
  write(*,'(a)') " TEST #2: evaluate standard 2D tensor product spline on test grid           "
  write(*,'(a)') " ***************************************************************************"

  call sll_s_hdf5_ser_file_create( 'mapping_interpol_czarny.h5', file_id, h5_error )

  ! Evaluate spline at interpolation points
  max_error = 0.0_wp
  do i2 = 1, npts2+1
    do i1 = 1, npts1
      eta(1) = real( i1-1, wp ) / real( npts1-1, wp )
      eta(2) = real( i2-1, wp ) * sll_p_twopi / real( npts2, wp )
      x  (:) = mapping_iga % eval( eta )
      interp_funct(i1,i2) = profile % eval( x(1), x(2) )
      interp_error(i1,i2) = profile % eval( x(1), x(2) ) - spline_2d % eval( eta(1), eta(2) )
      max_error = max( max_error, abs( interp_error(i1,i2) ) )
    end do
  end do

  ! Store mesh
  call sll_o_hdf5_ser_write_array( file_id, interp_funct, "/interp_funct", h5_error )
  call sll_o_hdf5_ser_write_array( file_id, abs( interp_error ), "/interp_error", h5_error )

  call sll_s_hdf5_ser_file_close ( file_id, h5_error )

  write(*,*)
  write(*,'(a,es8.2)') " Maximum error = ", max_error
  write(*,*)

  call mapping_iga % free()

  !-----------------------------------------------------------------------------
  ! Deallocate and free remaining objects
  !-----------------------------------------------------------------------------

  deallocate( tau1 )
  deallocate( tau2 )
  deallocate( gtau )

  call spline_2d % free()
  call spline_interpolator_2d % free()

end program test_polar_mapping
