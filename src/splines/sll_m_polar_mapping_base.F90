module sll_m_polar_mapping_base
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  use sll_m_working_precision, only: f64

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
    procedure(i_fun_eval)      , deferred :: eval
    procedure(i_fun_jacobian)  , deferred :: jacobian ! Jacobian determinant
    procedure(i_sub_store_data), deferred :: store_data

  end type sll_c_polar_mapping

  ! Interfaces for deferred procedures
  abstract interface

    function i_fun_eval( self, eta ) result( x )
      import sll_c_polar_mapping, wp
      class(sll_c_polar_mapping), intent(in) :: self
      real(wp)                  , intent(in) :: eta(2)
      real(wp) :: x(2)
    end function i_fun_eval

    function i_fun_jacobian( self, eta ) result( jacobian )
      import sll_c_polar_mapping, wp
      class(sll_c_polar_mapping), intent(in) :: self
      real(wp)                  , intent(in) :: eta(2)
      real(wp) :: jacobian
    end function i_fun_jacobian

    subroutine i_sub_store_data( self, n1, n2, file_name )
      import sll_c_polar_mapping
      class(sll_c_polar_mapping), intent(in) :: self
      integer                   , intent(in) :: n1
      integer                   , intent(in) :: n2
      character(len=*)          , intent(in) :: file_name
    end subroutine i_sub_store_data

  end interface

end module sll_m_polar_mapping_base
