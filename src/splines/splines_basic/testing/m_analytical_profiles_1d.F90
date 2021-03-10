module m_analytical_profiles_1d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   use m_analytical_profiles_1d_base, only: &
      t_profile_1d_info, &
      c_analytical_profile_1d

   use m_analytical_profiles_1d_cos, only: &
      t_analytical_profile_1d_cos

   use m_analytical_profiles_1d_poly, only: &
      t_analytical_profile_1d_poly

   implicit none

   public :: &
      t_profile_1d_info, &
      c_analytical_profile_1d, &
      t_analytical_profile_1d_cos, &
      t_analytical_profile_1d_poly

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end module m_analytical_profiles_1d
