!> @authors  Yaman Güçlü, IPP Garching
!> @authors Edoardo Zoni, IPP Garching
!>
!> @brief
!> $f(x) = \sum_{i=0}^n c_i x^i$

module m_analytical_profiles_1d_poly
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_errors.h"

   use sll_m_working_precision, only: &
      f64

   use sll_m_constants, only: &
      sll_p_twopi

   use m_analytical_profiles_1d_base, only: &
      t_profile_1d_info, &
      c_analytical_profile_1d

   implicit none

   public :: &
      t_analytical_profile_1d_poly

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   !> Working precision
   integer, parameter :: wp = f64

   !> Concrete type for 1D interpolation test-cases
   type, extends(c_analytical_profile_1d) :: t_analytical_profile_1d_poly

      ! Profile info (not overwritable by user)
      type(t_profile_1d_info), private :: info

      integer :: deg
      real(wp), allocatable :: coeffs(:)

   contains

      ! Constructor
      procedure :: init

      ! Abstract interface
      procedure :: get_info
      procedure :: eval
      procedure :: max_norm ! general formula not implemented for this profile

   end type t_analytical_profile_1d_poly

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   !-----------------------------------------------------------------------------
   !> Initialize profile
   !-----------------------------------------------------------------------------
   subroutine init(self, deg)
      class(t_analytical_profile_1d_poly), intent(out) :: self
      integer, intent(in) :: deg

      ! Initialize basic profile info
      self%info%xmin = -1.0_wp
      self%info%xmax = 1.0_wp
      self%info%periodic = .false.
      self%info%poly_order = deg

      ! Generate random polynomial coefficients
      if (allocated(self%coeffs)) deallocate (self%coeffs)
      allocate (self%coeffs(0:deg))
      self%deg = deg
      call random_number(self%coeffs)     ! 0 <= c < 1
      self%coeffs = 1.0_wp - self%coeffs  ! 0 < c <= 1

   end subroutine init

   !-----------------------------------------------------------------------------
   !> Get profile info
   !-----------------------------------------------------------------------------
   pure subroutine get_info(self, info)
      class(t_analytical_profile_1d_poly), intent(in) :: self
      type(t_profile_1d_info), intent(out) :: info

      info = self%info

   end subroutine get_info

   !-----------------------------------------------------------------------------
   !> Evaluate 1D profile (or one of its derivatives)
   !-----------------------------------------------------------------------------
   pure function eval(self, x, diff) result(f)
      class(t_analytical_profile_1d_poly), intent(in) :: self
      real(wp), intent(in) :: x
      integer, optional, intent(in) :: diff
      real(wp) :: f

      integer  :: d
      integer  :: i
      real(wp) :: c

      if (present(diff)) then; d = diff; else; d = 0; end if

      f = 0.0_wp
      if (d > self%deg) return

      do i = d, self%deg
         c = real(falling_factorial(i, d), kind=wp)*x**(i - d)
         f = f + self%coeffs(i)*c
      end do

   end function eval

   !-----------------------------------------------------------------------------
   !> Evaluate max norm of profile (or one of its derivatives) over domain
   !-----------------------------------------------------------------------------
   function max_norm(self, diff) result(norm)
      class(t_analytical_profile_1d_poly), intent(in) :: self
      integer, optional, intent(in) :: diff
      real(wp) :: norm

      integer  :: d

      if (present(diff)) then; d = diff; else; d = 0; end if

      if (self%info%xmax >= abs(self%info%xmin)) then
         ! For xmax >= |xmin|:
         ! max(|f^(d)(x)|) = f^(d)(xmax)
         norm = self%eval(self%info%xmax, d)
      else
         SLL_ERROR("t_analytical_profile_1d_poly % max_norm", "General formula not implemented")
      end if

   end function max_norm

   !-----------------------------------------------------------------------------
   !> Calculate falling factorial of x
   !> [ https://en.wikipedia.org/wiki/Falling_and_rising_factorials ]
   !-----------------------------------------------------------------------------
   pure function falling_factorial(x, n) result(c)
      integer, intent(in) :: x
      integer, intent(in) :: n
      integer :: c
      integer :: k
      c = 1
      do k = 0, n - 1
         c = c*(x - k)
      end do
   end function falling_factorial

end module m_analytical_profiles_1d_poly
