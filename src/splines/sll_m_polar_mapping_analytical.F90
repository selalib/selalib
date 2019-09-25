module sll_m_singular_mapping_analytic
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"

  use sll_m_working_precision, only: f64

  use sll_m_constants, only: sll_p_twopi

  use sll_m_singular_mapping_base, only: sll_c_singular_mapping

  use sll_m_hdf5_io_serial, only: &
    sll_t_hdf5_ser_handle     , &
    sll_s_hdf5_ser_file_create, &
    sll_s_hdf5_ser_file_close , &
    sll_o_hdf5_ser_write_array

  implicit none

  public :: sll_c_singular_mapping_analytic

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Working precision
  integer, parameter :: wp = f64

  !> Abstract type, analytical singular mapping
  !  (may contain common components/methods)
  type, extends(sll_c_singular_mapping), abstract :: sll_c_singular_mapping_analytic

  contains

    ! Deferred procedures
    procedure(i_fun_jmat_comp), deferred :: jmat_comp

    ! Non-deferred procedures
    procedure :: store_data => s_singular_mapping_analytic__store_data

  end type sll_c_singular_mapping_analytic

  ! Interfaces for deferred procedures
  abstract interface

    SLL_PURE function i_fun_jmat_comp( self, eta ) result( jmat_comp )
      import sll_c_singular_mapping_analytic, wp
      class(sll_c_singular_mapping_analytic), intent(in) :: self
      real(wp)                              , intent(in) :: eta(2)
      real(wp) :: jmat_comp(2,2)
    end function i_fun_jmat_comp

  end interface

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine s_singular_mapping_analytic__store_data( self, n1, n2, file_id )
    class(sll_c_singular_mapping_analytic), intent(in) :: self
    integer                               , intent(in) :: n1
    integer                               , intent(in) :: n2
    type(sll_t_hdf5_ser_handle)           , intent(in) :: file_id

    integer  :: i1, i2
    real(wp) :: eta(2), x(2)

    real(wp), allocatable :: x1(:,:), x2(:,:), jacobian(:,:)

    ! For hdf5 I/O
    integer :: error

    ! Allocate physical mesh
    allocate( x1( n1, n2+1 ) ) ! repeated point along eta2
    allocate( x2( n1, n2+1 ) ) ! repeated point along eta2

    ! Allocate Jacobian determinant
    allocate( jacobian( n1, n2+1 ) ) ! repeated point along eta2

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

  end subroutine s_singular_mapping_analytic__store_data

end module sll_m_singular_mapping_analytic
