program unit_test_2d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  use sll_m_coordinate_transformation_multipatch, only: &
    sll_t_coordinate_transformation_multipatch_2d

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type(sll_t_coordinate_transformation_multipatch_2d) :: mp

!  call mp%read_from_file("identity_mp_info.nml")
  call mp%read_from_file("circle_mp5_pts12")

  print *, 'connectivity patch 1, face 3', mp%get_connectivity(1, 3)
  print *, 'connectivity patch 3, face 1', mp%get_connectivity(3, 1)

  call mp%write_to_file()
  call mp%delete()

  print *, "PASSED"

end program unit_test_2d
