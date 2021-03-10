!> @authors Yaman Güçlü, IPP Garching
!>
!> @brief
!> Define 1D functions on real interval to test interpolation
!>

module m_analytical_profiles_1d_base
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   use sll_m_working_precision, only: f64

   implicit none

   public :: &
      t_profile_1d_info, &
      c_analytical_profile_1d

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   !> Working precision
   integer, parameter :: wp = f64

   !-----------------------------------------------------------------------------
   !> Structure with basic info about 1D profile
   !-----------------------------------------------------------------------------
   type :: t_profile_1d_info

      real(wp) :: xmin
      real(wp) :: xmax
      logical  :: periodic    ! .true. if profile is periodic
      integer  :: poly_order  ! polynomial order (-1 if not a polynomial)

   end type t_profile_1d_info

   !-----------------------------------------------------------------------------
   !> Abstract type for 2D analytical profile to be used for testing interpolation
   !-----------------------------------------------------------------------------
   type, abstract :: c_analytical_profile_1d

   contains
      ! Get profile info
      procedure(i_subr_get_info), deferred :: get_info

      ! Evaluate 1D profile (or one of its derivatives)
      procedure(i_func_eval_profile), deferred :: eval

      ! Evaluate max norm of profile (or one of its derivatives) over domain
      procedure(i_func_max_norm), deferred :: max_norm

   end type c_analytical_profile_1d

   !-----------------------------------------------------------------------------
   abstract interface

      ! Get profile info
      pure subroutine i_subr_get_info(self, info)
         import c_analytical_profile_1d, t_profile_1d_info
         class(c_analytical_profile_1d), intent(in) :: self
         type(t_profile_1d_info), intent(out) :: info
      end subroutine i_subr_get_info

      ! Evaluate 1D profile (or one of its derivatives)
      pure function i_func_eval_profile(self, x, diff) result(f)
         import c_analytical_profile_1d, wp
         class(c_analytical_profile_1d), intent(in) :: self
         real(wp), intent(in) :: x
         integer, optional, intent(in) :: diff
         real(wp) :: f
      end function i_func_eval_profile

      ! Evaluate max norm of profile (or one of its derivatives) over domain
      function i_func_max_norm(self, diff) result(norm)
         import c_analytical_profile_1d, wp
         class(c_analytical_profile_1d), intent(in) :: self
         integer, optional, intent(in) :: diff
         real(wp) :: norm
      end function i_func_max_norm

   end interface

end module m_analytical_profiles_1d_base
