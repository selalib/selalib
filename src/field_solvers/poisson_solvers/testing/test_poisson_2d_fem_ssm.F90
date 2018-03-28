program test_poisson_2d_fem_ssm
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_errors.h"

  use sll_m_working_precision, only: f64

  use sll_m_constants, only: sll_p_twopi

  use sll_m_utilities, only: sll_s_new_array_linspace

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_greville, &
    sll_p_periodic

  use sll_m_bsplines, only: &
    sll_c_bsplines, &
    sll_s_bsplines_new

  use sll_m_spline_interpolator_1d, only: sll_s_spline_1d_compute_num_cells 

  use sll_m_polar_bsplines_2d, only: sll_t_polar_bsplines_2d

  use sll_m_polar_mapping_analytical, only: sll_c_polar_mapping_analytical

  use sll_m_polar_mapping_analytical_target, only: sll_t_polar_mapping_analytical_target

  use sll_m_polar_mapping_analytical_czarny, only: sll_t_polar_mapping_analytical_czarny

  use sll_m_polar_mapping_iga, only: sll_t_polar_mapping_iga

  use sll_m_poisson_2d_fem_ssm_weak_form, only: sll_t_poisson_2d_fem_ssm_weak_form

  use sll_m_poisson_2d_fem_ssm_assembler, only: sll_t_poisson_2d_fem_ssm_assembler

  use sll_m_poisson_2d_fem_ssm_projector, only: sll_t_poisson_2d_fem_ssm_projector

  use sll_m_vector_space_real_arrays, only: sll_t_vector_space_real_1d

  use sll_m_linear_operator_matrix_dense, only: sll_t_linear_operator_matrix_dense

  use sll_m_conjugate_gradient, only: sll_t_conjugate_gradient

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

  ! To initialize B-splines (p1,p2 degrees)
  integer :: mm, n1, n2, p1, p2, ncells1, ncells2

  ! B-splines break points
  real(wp), allocatable :: breaks_eta1(:)
  real(wp), allocatable :: breaks_eta2(:)

  ! 1D B-splines
  class(sll_c_bsplines), allocatable :: bsplines_eta1
  class(sll_c_bsplines), allocatable :: bsplines_eta2

  ! Analytical and discrete mappings
  class(sll_c_polar_mapping_analytical), allocatable :: mapping_analytical
  type(sll_t_polar_mapping_iga) :: mapping_iga

  ! Polar B-splines
  type(sll_t_polar_bsplines_2d) :: polar_bsplines

  ! Number of finite elements
  integer :: Nk1, Nk2

  ! Number of Gauss-Legendre quadrature points
  integer :: Nq1, Nq2

  ! Quadrature points
  real(wp), allocatable :: quad_points_eta1(:,:)
  real(wp), allocatable :: quad_points_eta2(:,:)

  ! Quadrature weights
  real(wp), allocatable :: quad_weights_eta1(:,:)
  real(wp), allocatable :: quad_weights_eta2(:,:)

  ! 1D data (B-splines values and derivatives)
  real(wp), allocatable :: data_1d_eta1(:,:,:,:)
  real(wp), allocatable :: data_1d_eta2(:,:,:,:)

  ! 2D data (inverse metric tensor, integral volume, right hand side)
  real(wp), allocatable :: int_volume(:,:,:,:)
  real(wp), allocatable :: inv_metric(:,:,:,:,:,:)
  real(wp), allocatable :: data_2d_rhs(:,:,:,:)

  ! Stiffness and mass matrices and C1 projections
  real(wp), allocatable :: A (:,:)
  real(wp), allocatable :: M (:,:)
  real(wp), allocatable :: Ap(:,:)
  real(wp), allocatable :: Ap_temp(:,:)
  real(wp), allocatable :: Mp(:,:)

  ! Matrix for barycentric coordinates
  real(wp), allocatable :: L(:,:)

  ! Right hand side vector and C1 projection
  real(wp), allocatable :: b (:)
  real(wp), allocatable :: bp(:)
  real(wp), allocatable :: bp_temp(:)

  ! Solution and C1 projection
  real(wp), allocatable :: x (:)
  real(wp), allocatable :: xp(:)

  ! Weak form
  type(sll_t_poisson_2d_fem_ssm_weak_form) :: weak_form

  ! Assembler
  type(sll_t_poisson_2d_fem_ssm_assembler) :: assembler

  ! C1 projector
  type(sll_t_poisson_2d_fem_ssm_projector) :: projector

  ! Linear solver
  type(sll_t_vector_space_real_1d)         :: bp_vecsp
  type(sll_t_vector_space_real_1d)         :: xp_vecsp
  type(sll_t_linear_operator_matrix_dense) :: Ap_linop
  type(sll_t_conjugate_gradient)           :: cjsolver
  real(wp), parameter :: tol = 1.0e-14_wp
  logical , parameter :: verbose = .false.

  ! Auxiliary/temporary variables
  integer  :: i, j, k, k1, k2, q1, q2
  real(wp) :: cx, cy, jdet, eta(2), jmat(2,2)

  ! For hdf5 I/O
  type(sll_t_hdf5_ser_handle) :: file_id
  integer                     :: h5_error

  !-----------------------------------------------------------------------------
  ! Initialize B-splines basis functions
  !-----------------------------------------------------------------------------

  ! Number of degrees of freedom (control points) along s and theta
  mm = 10
  n1 = mm * 1
  n2 = mm * 2

  ! Spline degrees along s and theta
  p1 = 3
  p2 = 5

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

!  allocate( sll_t_polar_mapping_analytical_target :: mapping_analytical )
  allocate( sll_t_polar_mapping_analytical_czarny :: mapping_analytical )

  ! Initialize analytical mapping
  select type ( mapping_analytical )
    type is ( sll_t_polar_mapping_analytical_target )
      call mapping_analytical % init() ! circle
!      call mapping_analytical % init( x0=[0.0_wp,0.0_wp], d0=0.2_wp, e0=0.3_wp )
    type is ( sll_t_polar_mapping_analytical_czarny )
      call mapping_analytical % init( x0=[0.0_wp,0.0_wp], b =1.4_wp, e =0.3_wp )
  end select

  ! Initialize discrete mapping
  call mapping_iga % init( bsplines_eta1, bsplines_eta2, mapping_analytical )

  !-----------------------------------------------------------------------------
  ! Initialize polar B-splines
  !-----------------------------------------------------------------------------

  call polar_bsplines % init( bsplines_eta1, bsplines_eta2, mapping_iga )

  !-----------------------------------------------------------------------------
  ! Allocate and initialize quadrature points and weights
  !-----------------------------------------------------------------------------

  ! Number of finite elements
  Nk1 = size( breaks_eta1 )-1
  Nk2 = size( breaks_eta2 )-1

  ! Number of Gauss-Legendre quadrature points
  Nq1 = p1 + 1
  Nq2 = p2 + 1

  ! Allocate quadrature points along s and theta
  allocate( quad_points_eta1( Nq1, Nk1 ) )
  allocate( quad_points_eta2( Nq2, Nk2 ) )

  ! Allocate quadrature weights along s and theta
  allocate( quad_weights_eta1( Nq1, Nk1 ) )
  allocate( quad_weights_eta2( Nq2, Nk2 ) )

  ! Initialize quadrature points and weights along s
  do k1 = 1, Nk1
    quad_points_eta1 (:,k1) = sll_f_gauss_legendre_points ( Nq1, breaks_eta1(k1), breaks_eta1(k1+1) )
    quad_weights_eta1(:,k1) = sll_f_gauss_legendre_weights( Nq1, breaks_eta1(k1), breaks_eta1(k1+1) )
  end do

  ! Initialize quadrature points and weights along theta
  do k2 = 1, Nk2
    quad_points_eta2 (:,k2) = sll_f_gauss_legendre_points ( Nq2, breaks_eta2(k2), breaks_eta2(k2+1) )
    quad_weights_eta2(:,k2) = sll_f_gauss_legendre_weights( Nq2, breaks_eta2(k2), breaks_eta2(k2+1) )
  end do

  !-----------------------------------------------------------------------------
  ! Allocate and fill 1D data
  !-----------------------------------------------------------------------------

  ! Allocate 1D data along s and theta
  allocate( data_1d_eta1( Nq1, 2, 1+p1, Nk1 ) )
  allocate( data_1d_eta2( Nq2, 2, 1+p2, Nk2 ) )

  ! Initialize 1D data along s
  do k1 = 1, Nk1
    do q1 = 1, Nq1
      call bsplines_eta1 % eval_basis_and_n_derivs( &
        x      = quad_points_eta1(q1,k1), &
        n      = 1, &
        derivs = data_1d_eta1(q1,:,:,k1), &
        jmin   = i )
    end do
  end do

  ! Initialize 1D data along theta
  do k2 = 1, Nk2
    do q2 = 1, Nq2
      call bsplines_eta2 % eval_basis_and_n_derivs( &
        x      = quad_points_eta2(q2,k2), &
        n      = 1, &
        derivs = data_1d_eta2(q2,:,:,k2), &
        jmin   = i )
    end do
  end do

  !-----------------------------------------------------------------------------
  ! Allocate and fill 2D data
  !-----------------------------------------------------------------------------

  allocate( int_volume( Nq1, Nq2, Nk1, Nk2 ) )
  allocate( inv_metric( Nq1, Nq2, Nk1, Nk2, 2, 2 ) )

  ! Initialize 2D data
  do k2 = 1, Nk2
    do k1 = 1, Nk1
      do q2 = 1, Nq2
        do q1 = 1, Nq1

          eta  = [ quad_points_eta1(q1,k1), quad_points_eta2(q2,k2) ]

          jdet = mapping_iga % jdet( eta )
          jmat = mapping_iga % jmat( eta )

          ! Area associated to each quadrature point
          int_volume(q1,q2,k1,k2) = abs( jdet ) * quad_weights_eta1(q1,k1) * quad_weights_eta2(q2,k2)

          ! Inverse metric tensor
          inv_metric(q1,q2,k1,k2,1,1) = jmat(1,2)**2 + jmat(2,2)**2
          inv_metric(q1,q2,k1,k2,1,2) = - jmat(1,1)*jmat(1,2) - jmat(2,1)*jmat(2,2)
          inv_metric(q1,q2,k1,k2,2,1) = inv_metric(q1,q2,k1,k2,1,2) ! symmetry
          inv_metric(q1,q2,k1,k2,2,2) = jmat(1,1)**2 + jmat(2,1)**2
          inv_metric(q1,q2,k1,k2,:,:) = inv_metric(q1,q2,k1,k2,:,:) / jdet**2

        end do
      end do
    end do
  end do

  allocate( data_2d_rhs( Nq1, Nq2, Nk1, Nk2 ) )

  ! Initialize 2D discrete RHS
  do k2 = 1, Nk2
    do k1 = 1, Nk1
      do q2 = 1, Nq2
        do q1 = 1, Nq1

          eta = [ quad_points_eta1(q1,k1), quad_points_eta2(q2,k2) ]

          data_2d_rhs(q1,q2,k1,k2) = rhs( eta, mapping_iga )

        end do
      end do
    end do
  end do

  !-----------------------------------------------------------------------------
  ! Matrix of barycentric coordinates needed for C1 projection
  !-----------------------------------------------------------------------------

  allocate( L( 2*n2, 3 ) )

  do i = 1, 2*n2

    ! Indices to access control points
    j = (i-1) / n2 + 1
    k = modulo( (i-1), n2 ) + 1

    cx = mapping_iga % spline_2d_x1 % bcoef(j,k)
    cy = mapping_iga % spline_2d_x2 % bcoef(j,k)

    L(i,1) = polar_bsplines % eval_l0( cx, cy )
    L(i,2) = polar_bsplines % eval_l1( cx, cy )
    L(i,3) = polar_bsplines % eval_l2( cx, cy )

  end do

  !-----------------------------------------------------------------------------
  ! Initialize assembler and projector
  !-----------------------------------------------------------------------------

  call assembler % init( n1, n2, weak_form )
  call projector % init( n1, n2, L )

  !-----------------------------------------------------------------------------
  ! Allocate and fill stiffness and mass matrices
  !-----------------------------------------------------------------------------

  allocate( A( n1*n2, n1*n2 ) )
  allocate( M( n1*n2, n1*n2 ) )

  A = 0.0_wp
  M = 0.0_wp

  ! Cycle over finite elements
  do k2 = 1, Nk2
    do k1 = 1, Nk1
      call assembler % add_element_mat( &
        k1          , &
        k2          , &
        data_1d_eta1, &
        data_1d_eta2, &
        int_volume  , &
        inv_metric  , &
        A           , &
        M )
    end do
  end do

  !-----------------------------------------------------------------------------
  ! Allocate and fill right hand side
  !-----------------------------------------------------------------------------

  allocate( b( n1*n2 ) )

  b = 0.0_wp

  ! Cycle over finite elements
  do k2 = 1, Nk2
    do k1 = 1, Nk1
      call assembler % add_element_rhs( &
        k1          , &
        k2          , &
        data_1d_eta1, &
        data_1d_eta2, &
        data_2d_rhs , &
        int_volume  , &
        b )
    end do
  end do

  !-----------------------------------------------------------------------------
  ! Allocate and fill C1 projections of stiffness and mass matrices
  !-----------------------------------------------------------------------------

  allocate( Ap( 3+(n1-2)*n2, 3+(n1-2)*n2 ) )
  allocate( Mp( 3+(n1-2)*n2, 3+(n1-2)*n2 ) )

  call projector % change_basis_matrix( A, Ap )
  call projector % change_basis_matrix( M, Mp )

  allocate( Ap_temp( size(Ap,1), size(Ap,2) ), source=Ap )

  associate( nn => 3+(n1-2)*n2, idx => 3+(n1-3)*n2 )

!    do j = idx+1, nn
!      do i = idx+1, nn
!        Ap_temp(i,j) = 0.0_wp
!      end do
!    end do
    do i = idx+1, nn
      Ap_temp(i,:) = 0.0_wp
    end do
    do j = idx+1, nn
      Ap_temp(:,j) = 0.0_wp
    end do

  end associate

  !-----------------------------------------------------------------------------
  ! Allocate and fill C1 projection of right hand side
  !-----------------------------------------------------------------------------

  allocate( bp( 3+(n1-2)*n2 ) )

  call projector % change_basis_vector( b, bp )

  allocate( bp_temp( size(bp) ), source=bp )

  associate( nn => 3+(n1-2)*n2, idx => 3+(n1-3)*n2 )

    do i = idx+1, nn
      bp_temp(i) = 0.0_wp
    end do

  end associate

  !-----------------------------------------------------------------------------
  ! Allocate C1 projection of solution
  !-----------------------------------------------------------------------------

  allocate( xp( 3+(n1-2)*n2 ) )

  xp = 0.0_wp

  !-----------------------------------------------------------------------------
  ! Initialize and solve linear system
  !-----------------------------------------------------------------------------

  ! Construct linear operator from matrix Ap_temp
  call Ap_linop % init( Ap_temp )

  ! Construct vector space from vector bp_temp
  call bp_vecsp % attach( bp_temp )

  ! Construct vector space for solution
  call xp_vecsp % attach( xp )

  ! Initialize conjugate gradient solver
  call cjsolver % init( tol=tol, verbose=verbose, template_vector=xp_vecsp )

  ! Solve linear system Ap*xp=bp
  call cjsolver % solve( A=Ap_linop, b=bp_vecsp, x=xp_vecsp )

  ! Copy solution into array xp
  xp = xp_vecsp % array

  !-----------------------------------------------------------------------------
  ! Allocate and compute solution in tensor-product space
  !-----------------------------------------------------------------------------

  allocate( x( n1*n2 ) )

  call projector % change_basis_vector_inverse( xp, x )

  !-----------------------------------------------------------------------------
  ! HDF5 I/O
  !-----------------------------------------------------------------------------

  ! Create HDF5 file for output
  call sll_s_hdf5_ser_file_create( 'poisson_2d_fem_ssm.h5', file_id, h5_error )

  ! Write stiffness matrix
  call sll_o_hdf5_ser_write_array( file_id, A, "/A", h5_error )

  ! Write mass matrix
  call sll_o_hdf5_ser_write_array( file_id, M, "/M", h5_error )

  ! Write L matrix needed for projection
  call sll_o_hdf5_ser_write_array( file_id, L, "/L", h5_error )

  ! Write C1 projection of stiffness matrix
  call sll_o_hdf5_ser_write_array( file_id, Ap, "/Ap", h5_error )
  call sll_o_hdf5_ser_write_array( file_id, Ap_temp, "/Ap_temp", h5_error )

  ! Write C1 projection of mass matrix
  call sll_o_hdf5_ser_write_array( file_id, Mp, "/Mp", h5_error )

  ! Write right hand side
  call sll_o_hdf5_ser_write_array( file_id, b, "/b", h5_error )

  ! Write C1 projection of right hand side
  call sll_o_hdf5_ser_write_array( file_id, bp, "/bp", h5_error )
  call sll_o_hdf5_ser_write_array( file_id, bp_temp, "/bp_temp", h5_error )

  ! Write solution
  call sll_o_hdf5_ser_write_array( file_id, x, "/x", h5_error )

  ! Write C1 projection of solution
  call sll_o_hdf5_ser_write_array( file_id, xp, "/xp", h5_error )

  ! Close HDF5 file
  call sll_s_hdf5_ser_file_close ( file_id, h5_error )

  !-----------------------------------------------------------------------------
  ! Deallocate allocatables and free objects
  !-----------------------------------------------------------------------------

  deallocate( breaks_eta1 )
  deallocate( breaks_eta2 )
  deallocate( quad_points_eta1 )
  deallocate( quad_points_eta2 )
  deallocate( quad_weights_eta1 )
  deallocate( quad_weights_eta2 )
  deallocate( data_1d_eta1 )
  deallocate( data_1d_eta2 )
  deallocate( int_volume )
  deallocate( inv_metric )
  deallocate( A  )
  deallocate( M  )
  deallocate( Ap )
  deallocate( Mp )
  deallocate( b  )
  deallocate( bp )
  deallocate( L  )

  deallocate( mapping_analytical )

  call bsplines_eta1  % free()
  call bsplines_eta2  % free()
  call mapping_iga    % free()
  call polar_bsplines % free()

  call projector % free()

  call cjsolver % free()
  call Ap_linop % free()
  call bp_vecsp % delete()
  call xp_vecsp % delete()

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  SLL_PURE function rhs( eta, mapping )
    real(wp)                     , intent(in) :: eta(2)
    type(sll_t_polar_mapping_iga), intent(in) :: mapping
    real(wp) :: rhs

    real(wp) :: x(2)

    x   = mapping % eval( eta )

    rhs = sin( sll_p_twopi * x(1) ) * cos( sll_p_twopi * x(2) )

    return

  end function rhs

end program test_poisson_2d_fem_ssm
