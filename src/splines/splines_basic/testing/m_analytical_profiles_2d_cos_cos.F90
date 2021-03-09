!> @authors Yaman Güçlü, IPP Garching
!>
!> @brief
!> f(x,y) = cos(2\pi x1) cos(2\pi x2)

module m_analytical_profiles_2d_cos_cos
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   use sll_m_working_precision, only: &
      f64

   use sll_m_constants, only: &
      sll_p_pi, &
      sll_p_twopi

   use m_analytical_profiles_2d_base, only: &
      t_profile_2d_info, &
      c_analytical_profile_2d

   implicit none

   public :: &
      t_analytical_profile_2d_cos_cos

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   !> Working precision
   integer, parameter :: wp = f64

   !> Concrete type for 2D interpolation test-cases
   !> f(x1,x2) = cos( k1*x1 + phi1 ) * cos( k2*x2 + phi2 )
   type, extends(c_analytical_profile_2d) :: t_analytical_profile_2d_cos_cos

      ! Default parameters
      real(wp) :: k1 = sll_p_twopi
      real(wp) :: k2 = sll_p_twopi
      real(wp) :: phi1 = 0.0_wp
      real(wp) :: phi2 = 0.0_wp

      ! Profile info (not overwritable by user)
      type(t_profile_2d_info), private :: info

   contains

      ! Constructor
      procedure :: init

      ! Abstract interface
      procedure :: get_info
      procedure :: eval
      procedure :: max_norm

   end type t_analytical_profile_2d_cos_cos

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   !-----------------------------------------------------------------------------
   !> Initialize profile
   !-----------------------------------------------------------------------------
   subroutine init(self, n1, n2, c1, c2)
      class(t_analytical_profile_2d_cos_cos), intent(out) :: self
      integer, optional, intent(in) :: n1 ! k1 = 2*pi*n1
      integer, optional, intent(in) :: n2 ! k2 = 2*pi*n2
      real(wp), optional, intent(in) :: c1 ! phi1 = 2*pi*c1, c1 in [0,1)
      real(wp), optional, intent(in) :: c2 ! phi2 = 2*pi*c2, c2 in [0,1)

      ! Initialize basic profile info
      self%info%x1_min = 0.0_wp
      self%info%x1_max = 1.0_wp
      self%info%x2_min = 0.0_wp
      self%info%x2_max = 1.0_wp
      self%info%x1_periodic = .true.
      self%info%x2_periodic = .true.
      self%info%x1_poly_order = -1
      self%info%x2_poly_order = -1

      ! Overwrite defaults with user parameters
      if (present(n1)) self%k1 = sll_p_twopi*n1
      if (present(n2)) self%k2 = sll_p_twopi*n2
      if (present(c1)) self%phi1 = sll_p_twopi*c1
      if (present(c2)) self%phi2 = sll_p_twopi*c2

   end subroutine init

   !-----------------------------------------------------------------------------
   !> Get profile info
   !-----------------------------------------------------------------------------
   pure subroutine get_info(self, info)
      class(t_analytical_profile_2d_cos_cos), intent(in) :: self
      type(t_profile_2d_info), intent(out) :: info

      info = self%info

   end subroutine get_info

   !-----------------------------------------------------------------------------
   !> Evaluate 2D profile (or one of its derivatives)
   !-----------------------------------------------------------------------------
   pure function eval(self, x1, x2, diff_x1, diff_x2) result(f)
      class(t_analytical_profile_2d_cos_cos), intent(in) :: self
      real(wp), intent(in) :: x1
      real(wp), intent(in) :: x2
      integer, optional, intent(in) :: diff_x1
      integer, optional, intent(in) :: diff_x2
      real(wp) :: f

      integer :: d1, d2

      if (present(diff_x1)) then; d1 = diff_x1; else; d1 = 0; end if
      if (present(diff_x2)) then; d2 = diff_x2; else; d2 = 0; end if

      f = self%k1**d1*cos(0.5_wp*sll_p_pi*d1 + self%k1*x1 + self%phi1) &
          *self%k2**d2*cos(0.5_wp*sll_p_pi*d2 + self%k2*x2 + self%phi2)

   end function eval

   !-----------------------------------------------------------------------------
   !> Evaluate max norm of profile (or one of its derivatives) over domain
   !-----------------------------------------------------------------------------
   function max_norm(self, diff_x1, diff_x2) result(norm)
      class(t_analytical_profile_2d_cos_cos), intent(in) :: self
      integer, optional, intent(in) :: diff_x1
      integer, optional, intent(in) :: diff_x2
      real(wp) :: norm

      integer :: d1, d2

      if (present(diff_x1)) then; d1 = diff_x1; else; d1 = 0; end if
      if (present(diff_x2)) then; d2 = diff_x2; else; d2 = 0; end if

      norm = self%k1**d1*self%k2**d2

   end function max_norm

end module m_analytical_profiles_2d_cos_cos
