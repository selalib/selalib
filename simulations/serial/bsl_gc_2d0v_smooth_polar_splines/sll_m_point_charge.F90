module sll_m_point_charge
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   use sll_m_working_precision, only: f64

   implicit none

   public :: sll_t_point_charge

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   ! Working precision
   integer, parameter :: wp = f64

   type :: sll_t_point_charge

      real(wp) :: intensity
      real(wp) :: location(2)

   contains

      procedure :: init => s_point_charge__init

   end type sll_t_point_charge

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine s_point_charge__init(self, intensity, location)
      class(sll_t_point_charge), intent(inout) :: self
      real(wp), intent(in) :: intensity
      real(wp), intent(in) :: location(2)

      self%intensity = intensity
      self%location(:) = location(:)

   end subroutine s_point_charge__init

end module sll_m_point_charge
