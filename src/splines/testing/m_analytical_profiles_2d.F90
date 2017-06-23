module m_analytical_profiles_2d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  use m_analytical_profiles_2d_base, only: &
    t_profile_2d_info, &
    c_analytical_profile_2d

  use m_analytical_profiles_2d_cos_cos, only: &
    t_analytical_profile_2d_cos_cos

  implicit none

  public :: &
    t_profile_2d_info, &
    c_analytical_profile_2d, &
    t_analytical_profile_2d_cos_cos

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end module m_analytical_profiles_2d
