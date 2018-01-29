module sll_m_polar_mapping_analytical
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  use sll_m_working_precision, only: f64

  use sll_m_constants, only: sll_p_twopi

  use sll_m_polar_mapping_base, only: sll_c_polar_mapping

  use sll_m_hdf5_io_serial, only: &
    sll_t_hdf5_ser_handle     , &
    sll_s_hdf5_ser_file_create, &
    sll_s_hdf5_ser_file_close , &
    sll_o_hdf5_ser_write_array

  implicit none

  public :: sll_c_polar_mapping_analytical

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Working precision
  integer, parameter :: wp = f64

  !> Abstract type, analytical polar mapping
  !  (may contain common components/methods)
  type, extends(sll_c_polar_mapping), abstract :: sll_c_polar_mapping_analytical

  contains

    procedure :: store_data => s_polar_mapping_analytical__store_data

  end type sll_c_polar_mapping_analytical

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine s_polar_mapping_analytical__store_data( self, n1, n2, file_name )
    class(sll_c_polar_mapping_analytical), intent(in) :: self
    integer                              , intent(in) :: n1
    integer                              , intent(in) :: n2
    character(len=*)                     , intent(in) :: file_name

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
        jacobian(i1,i2) = self % jdet( eta )
      end do
    end do

    ! Store physical mesh
    call sll_o_hdf5_ser_write_array( file_id, x1, "/x1", error )
    call sll_o_hdf5_ser_write_array( file_id, x2, "/x2", error )

    ! Store Jacobian determinant
    call sll_o_hdf5_ser_write_array( file_id, jacobian, "/jacobian", error )

    call sll_s_hdf5_ser_file_close ( file_id, error )

  end subroutine s_polar_mapping_analytical__store_data

end module sll_m_polar_mapping_analytical
