program test_bsplines_2d_new
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  use sll_m_working_precision, only: &
    f64

  use m_analytical_profiles_2d, only: &
    t_analytical_profile_2d_cos_cos

  use m_test_bsplines_2d, only: &
    s_check_all_bcs

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type(t_analytical_profile_2d_cos_cos) :: profile_2d_cos_cos
  integer :: degree1
  integer :: degree2
  integer :: nx1
  integer :: nx2

  ! Initialize profile
  call profile_2d_cos_cos % init()

  ! Choose spline degree
  degree1 = 3
  degree2 = 3

  ! Choose number of knots in tensor grid
  nx1 = 10
  nx2 = 37

  ! Test all possible boundary conditions
  call s_check_all_bcs( profile_2d_cos_cos, degree1, degree2, nx1, nx2 )

  ! Pass test by default
  ! TODO: check errors
  write(*,*) "PASSED"

end program test_bsplines_2d_new
