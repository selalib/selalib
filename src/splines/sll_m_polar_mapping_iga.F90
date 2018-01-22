module sll_m_polar_mapping_iga
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  use sll_m_working_precision, only: f64

  use sll_m_constants, only: sll_p_twopi

  use sll_m_polar_mapping_base, only: sll_c_polar_mapping

  use sll_m_polar_mapping_analytical, only: sll_c_polar_mapping_analytical

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_greville, &
    sll_p_periodic

  use sll_m_bsplines, only: sll_c_bsplines

  use sll_m_spline_2d, only: sll_t_spline_2d

  use sll_m_spline_interpolator_2d, only: sll_t_spline_interpolator_2d

  use sll_m_hdf5_io_serial, only: &
    sll_t_hdf5_ser_handle     , &
    sll_s_hdf5_ser_file_create, &
    sll_s_hdf5_ser_file_close , &
    sll_o_hdf5_ser_write_array

  public :: sll_t_polar_mapping_iga

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Working precision
  integer, parameter :: wp = f64

  !> Concrete type, iga polar mapping
  type, extends(sll_c_polar_mapping) :: sll_t_polar_mapping_iga

    ! Number of basis functions along eta1 and eta2
    integer :: nbasis_eta1
    integer :: nbasis_eta2

    ! Number of basis functions along eta1 and eta2
    integer :: degree_eta1
    integer :: degree_eta2

    ! Interpolation points
    real(wp), allocatable :: tau_eta1(:)
    real(wp), allocatable :: tau_eta2(:)

    real(wp), allocatable :: gtau_x1(:,:) 
    real(wp), allocatable :: gtau_x2(:,:) 

    ! 2D splines
    type(sll_t_spline_2d) :: spline_2d_x1
    type(sll_t_spline_2d) :: spline_2d_x2

    ! 2D spline interpolator
    type(sll_t_spline_interpolator_2d) :: spline_interp_eta12

  contains

    procedure :: init       => s_polar_mapping_iga__init
    procedure :: eval       => f_polar_mapping_iga__eval
    procedure :: jacobian   => f_polar_mapping_iga__jacobian ! Jacobian determinant
    procedure :: store_data => s_polar_mapping_iga__store_data
    procedure :: free       => s_polar_mapping_iga__free

  end type sll_t_polar_mapping_iga

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine s_polar_mapping_iga__init( self, spline_basis_eta1, spline_basis_eta2, mapping_analytical) 
    class(sll_t_polar_mapping_iga)       , intent(inout) :: self
    class(sll_c_bsplines)                , intent(in   ) :: spline_basis_eta1
    class(sll_c_bsplines)                , intent(in   ) :: spline_basis_eta2
    class(sll_c_polar_mapping_analytical), intent(in   ) :: mapping_analytical

    integer  :: i1, i2
    real(wp) :: eta(2), x(2)

    ! Store number of basis functions along eta1 and eta2
    self % nbasis_eta1 = spline_basis_eta1 % nbasis
    self % nbasis_eta2 = spline_basis_eta2 % nbasis

    ! Store degree of basis functions along eta1 and eta2
    self % degree_eta1 = spline_basis_eta1 % degree
    self % degree_eta2 = spline_basis_eta2 % degree

    ! Initialize 2D spline along x1
    call self % spline_2d_x1 % init( &
      bsplines_x1 = spline_basis_eta1, &
      bsplines_x2 = spline_basis_eta2 )

    ! Initialize 2D spline along x2
    call self % spline_2d_x2 % init( &
      bsplines_x1 = spline_basis_eta1, &
      bsplines_x2 = spline_basis_eta2 )

    ! Initialize 2D spline interpolator
    call self % spline_interp_eta12 % init( &
      bspl1 = spline_basis_eta1, &
      bspl2 = spline_basis_eta2, &
      bc_xmin = [sll_p_greville,sll_p_periodic], &
      bc_xmax = [sll_p_greville,sll_p_periodic] )

    ! Get interpolation points from spline interpolator
    call self % spline_interp_eta12 % get_interp_points( &
      tau1 = self % tau_eta1, &
      tau2 = self % tau_eta2 )

    allocate( self % gtau_x1( size( self % tau_eta1 ), size( self % tau_eta2 ) ) )
    allocate( self % gtau_x2( size( self % tau_eta1 ), size( self % tau_eta2 ) ) )

    associate( ntau1 => size( self % tau_eta1 ), ntau2 => size( self % tau_eta2 ) )

      ! Evaluate analytical mapping on interpolation points
      do i2 = 1, ntau2
        do i1 = 1, ntau1
          eta(1) = self % tau_eta1(i1)
          eta(2) = self % tau_eta2(i2)
          x  (:) = mapping_analytical % eval( eta )
          self % gtau_x1(i1,i2) = x(1)
          self % gtau_x2(i1,i2) = x(2)
        end do
      end do

    end associate

    ! Compute interpolant 2D spline along x1
    call self % spline_interp_eta12 % compute_interpolant( &
      spline = self % spline_2d_x1, &
      gtau   = self % gtau_x1 )

    ! Compute interpolant 2D spline along x2
    call self % spline_interp_eta12 % compute_interpolant( &
      spline = self % spline_2d_x2, &
      gtau   = self % gtau_x2 )

  end subroutine s_polar_mapping_iga__init

  !-----------------------------------------------------------------------------
  function f_polar_mapping_iga__eval( self, eta ) result( x )
    class(sll_t_polar_mapping_iga), intent(in) :: self
    real(wp)                      , intent(in) :: eta(2)
    real(wp) :: x(2)

    x(1) = self % spline_2d_x1 % eval( eta(1), eta(2) )
    x(2) = self % spline_2d_x2 % eval( eta(1), eta(2) )

  end function f_polar_mapping_iga__eval

  !-----------------------------------------------------------------------------
  function f_polar_mapping_iga__jacobian( self, eta ) result( jacobian )
    class(sll_t_polar_mapping_iga), intent(in) :: self
    real(wp)                      , intent(in) :: eta(2)
    real(wp) :: jacobian

    real(wp) :: d1, d2, d3, d4

    ! d1 = d(x1)/d(eta1)
    ! d2 = d(x1)/d(eta2)
    ! d3 = d(x2)/d(eta1)
    ! d4 = d(x2)/d(eta2)
    d1 = self % spline_2d_x1 % eval_deriv_x1( eta(1), eta(2) )
    d2 = self % spline_2d_x1 % eval_deriv_x2( eta(1), eta(2) )
    d3 = self % spline_2d_x2 % eval_deriv_x1( eta(1), eta(2) )
    d4 = self % spline_2d_x2 % eval_deriv_x2( eta(1), eta(2) )

    jacobian = d1*d4 - d2*d3

  end function f_polar_mapping_iga__jacobian

  !-----------------------------------------------------------------------------
  subroutine s_polar_mapping_iga__store_data( self, n1, n2, file_name )
    class(sll_t_polar_mapping_iga), intent(in) :: self
    integer                       , intent(in) :: n1
    integer                       , intent(in) :: n2
    character(len=*)              , intent(in) :: file_name

    integer  :: i1, i2
    real(wp) :: eta(2), x(2)

    real(wp), allocatable :: x1(:,:), x2(:,:), jacobian(:,:)

    ! For hdf5 I/O
    type(sll_t_hdf5_ser_handle) :: file_id
    integer                     :: error

    ! Allocate physical mesh
    allocate( x1( n1, n2+1 ) ) ! repeated point along eta2
    allocate( x2( n1, n2+1 ) ) ! repeated point along eta2

    ! Allocate Jacobian determinant
    allocate( jacobian( n1, n2+1 ) ) ! repeated point along eta2

    ! Create HDF5 file
    call sll_s_hdf5_ser_file_create( trim( file_name )//'.h5', file_id, error )

    ! Compute physical mesh and Jacobian determinant
    do i2 = 1, n2+1 ! repeated point along eta2
      do i1 = 1, n1
        eta(1) = real( i1-1, wp ) / real( n1-1, wp )
        eta(2) = real( i2-1, wp ) * sll_p_twopi / real( n2, wp )
        x  (:) = self % eval( eta )
        x1(i1,i2) = x(1)
        x2(i1,i2) = x(2)
        jacobian(i1,i2) = self % jacobian( eta )
      end do
    end do

    ! Store physical mesh
    call sll_o_hdf5_ser_write_array( file_id, x1, "/x1", error )
    call sll_o_hdf5_ser_write_array( file_id, x2, "/x2", error )

    ! Store Jacobian determinant
    call sll_o_hdf5_ser_write_array( file_id, jacobian, "/jacobian", error )

    ! Store control points
    associate( nb1 => self % nbasis_eta1, nb2 => self % nbasis_eta2 )
      call sll_o_hdf5_ser_write_array( file_id, self % spline_2d_x1 % bcoef( 1:nb1, 1:nb2 ), "/c_x1", error )
      call sll_o_hdf5_ser_write_array( file_id, self % spline_2d_x2 % bcoef( 1:nb1, 1:nb2 ), "/c_x2", error )
    end associate

    call sll_s_hdf5_ser_file_close ( file_id, error )

  end subroutine s_polar_mapping_iga__store_data

  !-----------------------------------------------------------------------------
  subroutine s_polar_mapping_iga__free( self )
    class(sll_t_polar_mapping_iga), intent(inout) :: self

    deallocate( self % tau_eta1 )
    deallocate( self % tau_eta2 )
    deallocate( self % gtau_x1  )
    deallocate( self % gtau_x2  )

    call self % spline_2d_x1 % free()
    call self % spline_2d_x2 % free()

    call self % spline_interp_eta12 % free()

  end subroutine s_polar_mapping_iga__free

end module sll_m_polar_mapping_iga
