module sll_m_polar_mapping_base
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  use sll_m_working_precision, only: f64

  use sll_m_constants, only: sll_p_twopi

  use sll_m_hdf5_io_serial, only: &
    sll_t_hdf5_ser_handle     , &
    sll_s_hdf5_ser_file_create, &
    sll_s_hdf5_ser_file_close , &
    sll_o_hdf5_ser_write_array

  implicit none

  public :: sll_c_polar_mapping

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Working precision
  integer, parameter :: wp = f64

  !> Abstract  type, polar mapping
  type, abstract :: sll_c_polar_mapping

  contains

    ! Deferred procedures
    procedure(i_fun_eval), deferred :: eval

    ! Non-deferred procedures
    procedure :: write_mesh => s_polar_mapping__write_mesh

  end type sll_c_polar_mapping

  ! Interfaces for deferred procedures
  abstract interface

    function i_fun_eval( self, eta ) result( x )
      import sll_c_polar_mapping, wp
      class(sll_c_polar_mapping), intent(in) :: self
      real(wp)                  , intent(in) :: eta(2)
      real(wp) :: x(2)
    end function i_fun_eval

  end interface

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine s_polar_mapping__write_mesh( self, n1, n2, file_name )
    class(sll_c_polar_mapping), intent(inout) :: self
    integer                   , intent(in   ) :: n1
    integer                   , intent(in   ) :: n2
    character(len=*)          , intent(in   ) :: file_name

    integer  :: i1, i2
    real(wp) :: eta(2), x(2)

    real(wp), allocatable :: x1(:,:)
    real(wp), allocatable :: x2(:,:)

    ! For hdf5 I/O
    type(sll_t_hdf5_ser_handle) :: file_id
    integer                     :: error

    ! Allocate physical mesh
    allocate( x1( n1, n2+1 ) ) ! repeated point along eta2
    allocate( x2( n1, n2+1 ) ) ! repeated point along eta2

    ! Initialize physical mesh (assuming domain [0,1]x[0,2\pi] )
    do i2 = 1, n2+1 ! repeated point along eta2
      do i1 = 1, n1
        eta(1) = real( i1-1, wp ) / real( n1-1, wp )
        eta(2) = real( i2-1, wp ) * sll_p_twopi / real( n2, wp )
        x  (:) = self % eval( eta )
        x1(i1,i2) = x(1)
        x2(i1,i2) = x(2)
      end do
    end do

    call sll_s_hdf5_ser_file_create( trim( file_name )//".h5", file_id, error )
    call sll_o_hdf5_ser_write_array( file_id, x1, "/x1", error )
    call sll_o_hdf5_ser_write_array( file_id, x2, "/x2", error )
    call sll_s_hdf5_ser_file_close ( file_id, error ) 

  end subroutine s_polar_mapping__write_mesh

end module sll_m_polar_mapping_base
