!> @authors Yaman Güçlü , IPP Garching
!> @authors Edoardo Zoni, IPP Garching
!>
!> @brief
!> f(x) = cos(2\pi x)

module m_analytical_profiles_1d_cos
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   use sll_m_working_precision, only: &
      f64

   use sll_m_constants, only: &
      sll_p_pi, &
      sll_p_twopi

   use m_analytical_profiles_1d_base, only: &
      t_profile_1d_info, &
      c_analytical_profile_1d

   implicit none

   public :: &
      t_analytical_profile_1d_cos

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   !> Working precision
   integer, parameter :: wp = f64

   !> Concrete type for 1D interpolation test-cases
   !> f(x) = cos( k*x + phi )
   type, extends(c_analytical_profile_1d) :: t_analytical_profile_1d_cos

      ! Default parameters
      real(wp) :: k = sll_p_twopi
      real(wp) :: phi = 0.0_wp

      ! Profile info (not overwritable by user)
      type(t_profile_1d_info), private :: info

   contains

      ! Constructor
      procedure :: init

      ! Abstract interface
      procedure :: get_info
      procedure :: eval
      procedure :: max_norm

   end type t_analytical_profile_1d_cos

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   !-----------------------------------------------------------------------------
   !> Initialize profile
   !-----------------------------------------------------------------------------
   subroutine init(self, n, c)
      class(t_analytical_profile_1d_cos), intent(out) :: self
      integer, optional, intent(in) :: n ! k   = 2*pi*n
      real(wp), optional, intent(in) :: c ! phi = 2*pi*c, c in [0,1)

      ! Initialize basic profile info
      self%info%xmin = 0.0_wp
      self%info%xmax = 1.0_wp
      self%info%periodic = .true.
      self%info%poly_order = -1

      ! Overwrite defaults with user parameters
      if (present(n)) self%k = sll_p_twopi*n
      if (present(c)) self%phi = sll_p_twopi*c

   end subroutine init

   !-----------------------------------------------------------------------------
   !> Get profile info
   !-----------------------------------------------------------------------------
   pure subroutine get_info(self, info)
      class(t_analytical_profile_1d_cos), intent(in) :: self
      type(t_profile_1d_info), intent(out) :: info

      info = self%info

   end subroutine get_info

   !-----------------------------------------------------------------------------
   !> Evaluate 1D profile (or one of its derivatives)
   !-----------------------------------------------------------------------------
   pure function eval(self, x, diff) result(f)
      class(t_analytical_profile_1d_cos), intent(in) :: self
      real(wp), intent(in) :: x
      integer, optional, intent(in) :: diff
      real(wp) :: f

      integer :: d

      if (present(diff)) then; d = diff; else; d = 0; end if

      f = self%k**d*cos(0.5_wp*sll_p_pi*d + self%k*x + self%phi)

   end function eval

   !-----------------------------------------------------------------------------
   !> Evaluate max norm of profile (or one of its derivatives) over domain
   !-----------------------------------------------------------------------------
   function max_norm(self, diff) result(norm)
      class(t_analytical_profile_1d_cos), intent(in) :: self
      integer, optional, intent(in) :: diff
      real(wp) :: norm

      integer :: d

      if (present(diff)) then; d = diff; else; d = 0; end if

      norm = self%k**d

   end function max_norm

end module m_analytical_profiles_1d_cos
