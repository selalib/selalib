program unit_test_2d
#include "sll_working_precision.h"
  use sll_coordinate_transformation_multipatch_module
  use sll_distribution_function_4d_multipatch_module
  implicit none

  type(sll_coordinate_transformation_multipatch_2d) :: t_mp
  type(sll_distribution_function_4d_multipatch), pointer  :: f_mp
  sll_int32 :: npts_x3
  sll_int32 :: npts_x4

  npts_x3 = 32
  npts_x4 = 32

!  call mp%read_from_file("identity_mp_info.nml")
  call t_mp%read_from_file("square_4p_n10")

  f_mp => sll_new_distribution_function_4d_multipatch( t_mp, npts_x3, npts_x4 )


!  call mp%delete()

  print *, "PASSED"

end program unit_test_2d
