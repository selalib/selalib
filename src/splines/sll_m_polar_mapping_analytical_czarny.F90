module sll_m_polar_mapping_analytical_czarny
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  use sll_m_working_precision, only: f64

  use sll_m_constants, only: sll_p_twopi

  use sll_m_polar_mapping_analytical, only: sll_c_polar_mapping_analytical

  implicit none

  public :: sll_t_polar_mapping_analytical_czarny

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Working precision
  integer, parameter :: wp = f64

  !> Concrete type, analytical polar mapping
  type, extends(sll_c_polar_mapping_analytical) :: sll_t_polar_mapping_analytical_czarny

    ! Default values, can be overwritten from 'init' method
    real(wp) :: b = 1.4_wp
    real(wp) :: e = 0.3_wp

  contains

    procedure :: init => s_polar_mapping_analytical_czarny__init
    procedure :: eval => f_polar_mapping_analytical_czarny__eval

  end type sll_t_polar_mapping_analytical_czarny

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine s_polar_mapping_analytical_czarny__init( self, b, e )
    class(sll_t_polar_mapping_analytical_czarny), intent(inout) :: self
    real(wp), optional                          , intent(in   ) :: b
    real(wp), optional                          , intent(in   ) :: e

    ! Overwrite parameters
    if ( present( b ) ) self % b = b 
    if ( present( e ) ) self % e = e

  end subroutine s_polar_mapping_analytical_czarny__init

  !-----------------------------------------------------------------------------
  function f_polar_mapping_analytical_czarny__eval( self, eta ) result( x )
    class(sll_t_polar_mapping_analytical_czarny), intent(in) :: self
    real(wp)                                    , intent(in) :: eta(2)
    real(wp) :: x(2)

    associate( b => self % b, e => self % e )

      x(1) = (1.0_wp-sqrt( 1.0_wp + e*(e+2.0_wp*eta(1)*cos( eta(2) )) ))/e
      x(2) = b*eta(1)*sin( eta(2) )/((1.0_wp+e*x(1))*sqrt( 1.0_wp-e**2*0.25_wp ))

    end associate

  end function f_polar_mapping_analytical_czarny__eval

end module sll_m_polar_mapping_analytical_czarny
