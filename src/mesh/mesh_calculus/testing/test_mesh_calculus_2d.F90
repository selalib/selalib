program test_mesh_calculus
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use sll_m_cartesian_meshes, only: &
    sll_f_new_cartesian_mesh_2d, &
    sll_t_cartesian_mesh_2d

  use sll_m_common_coordinate_transformations, only: &
    sll_f_deriv_x1_polar_f_eta1, &
    sll_f_deriv_x1_polar_f_eta2, &
    sll_f_deriv_x2_polar_f_eta1, &
    sll_f_deriv_x2_polar_f_eta2, &
    sll_f_identity_jac11, &
    sll_f_identity_jac12, &
    sll_f_identity_jac21, &
    sll_f_identity_jac22, &
    sll_f_identity_x1, &
    sll_f_identity_x2, &
    sll_f_x1_polar_f, &
    sll_f_x2_polar_f

  use sll_m_coordinate_transformation_2d_base, only: &
    sll_c_coordinate_transformation_2d_base

  use sll_m_coordinate_transformations_2d, only: &
    sll_f_new_coordinate_transformation_2d_analytic

  use sll_m_mesh_calculus_2d, only: &
    sll_f_cell_volume, &
    sll_f_edge_length_eta1_minus, &
    sll_f_edge_length_eta1_plus, &
    sll_f_edge_length_eta2_minus, &
    sll_f_edge_length_eta2_plus, &
    sll_f_normal_integral_eta1_plus

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#define NCELLS1 128
#define NCELLS2 128

  type(sll_t_cartesian_mesh_2d), pointer                    :: m
  class(sll_c_coordinate_transformation_2d_base), pointer :: Ti   ! identity 
  class(sll_c_coordinate_transformation_2d_base), pointer :: Tp   ! polar
  sll_real64 :: volume, length_east, length_west, length_north, length_south
  sll_real64, dimension(2) :: int_east
  sll_int32  :: i, j

  m => sll_f_new_cartesian_mesh_2d( &
       NCELLS1, &
       NCELLS2 )

  Ti => sll_f_new_coordinate_transformation_2d_analytic( &
       "identity_mesh_calculus", &
       m, &
       sll_f_identity_x1, &
       sll_f_identity_x2, &
       sll_f_identity_jac11, &
       sll_f_identity_jac12, &
       sll_f_identity_jac21, &
       sll_f_identity_jac22, &
       (/0.0_f64, 0.0_f64/) ) ! this doesn't matter for identity transformation

  Tp => sll_f_new_coordinate_transformation_2d_analytic( &
       "polar_mesh_calculus", &
       m, &
       sll_f_x1_polar_f, &
       sll_f_x2_polar_f, &
       sll_f_deriv_x1_polar_f_eta1, &
       sll_f_deriv_x1_polar_f_eta2, &
       sll_f_deriv_x2_polar_f_eta1, &
       sll_f_deriv_x2_polar_f_eta2, &
       (/0.1_f64, 1.0_f64/) )

  volume = 0.0_f64
  do j=1,NCELLS1
     do i=1,NCELLS2
        volume = volume + sll_f_cell_volume(Ti,i,j,3)
     end do
  end do
  print *, 'volume = ', volume

  length_east  = 0.0_f64
  length_west  = 0.0_f64
  length_north = 0.0_f64
  length_south = 0.0_f64

  do i=1,NCELLS1
     length_north = length_north + sll_f_edge_length_eta2_plus(Ti,i,1,3)
     length_south = length_south + sll_f_edge_length_eta2_minus(Ti,i,1,3)
  end do

  do j=1,NCELLS2
     length_east  = length_east + sll_f_edge_length_eta1_plus(Ti,1,j,3)
     length_west  = length_west + sll_f_edge_length_eta1_minus(Ti,1,j,3)
  end do

  print *, 'length east edge = ', length_east, length_west, length_north, &
       length_south

  print *, '---------------------------------------------------------------'
  print *, '                   TEST OF NORMAL INTEGRALS                    '
  print *, '---------------------------------------------------------------'

  int_east = sll_f_normal_integral_eta1_plus(Ti,1,1,3)
  print *, 'integral of the normal vector on east edge = ', int_east

  print *, 'PASSED'

end program test_mesh_calculus
