program test_mesh_calculus
#include "sll_working_precision.h" 
#include "sll_assert.h"
  use sll_module_coordinate_transformations_2d
  use sll_common_coordinate_transformations
  use sll_logical_meshes
  use sll_mesh_calculus_2d_module
  implicit none

#define NCELLS1 128
#define NCELLS2 128

  type(sll_logical_mesh_2d), pointer                    :: m
  class(sll_coordinate_transformation_2d_base), pointer :: Ti   ! identity 
  class(sll_coordinate_transformation_2d_base), pointer :: Tp   ! polar
  sll_real64 :: volume, length_east, length_west, length_north, length_south
  sll_real64, dimension(2) :: int_east
  sll_int32  :: i, j

  m => new_logical_mesh_2d( &
       NCELLS1, &
       NCELLS2 )

  Ti => new_coordinate_transformation_2d_analytic( &
       "identity_mesh_calculus", &
       m, &
       identity_x1, &
       identity_x2, &
       identity_jac11, &
       identity_jac12, &
       identity_jac21, &
       identity_jac22, &
       (/0.0_f64, 0.0_f64/) ) ! this doesn't matter for identity transformation

  Tp => new_coordinate_transformation_2d_analytic( &
       "polar_mesh_calculus", &
       m, &
       x1_polar_f, &
       x2_polar_f, &
       deriv_x1_polar_f_eta1, &
       deriv_x1_polar_f_eta2, &
       deriv_x2_polar_f_eta1, &
       deriv_x2_polar_f_eta2, &
       (/0.1_f64, 1.0_f64/) )

  volume = 0.0_f64
  do j=1,NCELLS1
     do i=1,NCELLS2
        volume = volume + cell_volume(Ti,i,j,3)
     end do
  end do
  print *, 'volume = ', volume

  length_east  = 0.0_f64
  length_west  = 0.0_f64
  length_north = 0.0_f64
  length_south = 0.0_f64

  do i=1,NCELLS1
     length_north = length_north + edge_length_eta2_plus(Ti,i,1,3)
     length_south = length_south + edge_length_eta2_minus(Ti,i,1,3)
  end do

  do j=1,NCELLS2
     length_east  = length_east + edge_length_eta1_plus(Ti,1,j,3)
     length_west  = length_west + edge_length_eta1_minus(Ti,1,j,3)
  end do

  print *, 'length east edge = ', length_east, length_west, length_north, &
       length_south

  print *, '---------------------------------------------------------------'
  print *, '                   TEST OF NORMAL INTEGRALS                    '
  print *, '---------------------------------------------------------------'

  int_east = normal_integral_eta1_plus(Ti,1,1,3)
  print *, 'integral of the normal vector on east edge = ', int_east

  print *, 'PASSED'

end program test_mesh_calculus
