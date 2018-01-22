module sll_m_polar_mapping_analytical
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  use sll_m_polar_mapping_base, only: sll_c_polar_mapping

  implicit none

  public :: sll_c_polar_mapping_analytical

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Abstract type, analytical polar mapping
  !  (may contain common components/methods)
  type, extends(sll_c_polar_mapping), abstract :: sll_c_polar_mapping_analytical
  end type sll_c_polar_mapping_analytical

end module sll_m_polar_mapping_analytical
