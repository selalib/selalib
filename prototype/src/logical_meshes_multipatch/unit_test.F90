program unit_test_logical_meshes_multipatch
  use sll_logical_meshes_multipatch
#include "sll_working_precision.h"
  implicit none

#define NUM_PATCHES 5
#define NUM_CELLS1 32
#define NUM_CELLS2 32

  type(sll_logical_mesh_multipatch_2d), pointer :: mp2d
  sll_int32 :: i

  mp2d => new_logical_mesh_multipatch_2d( NUM_PATCHES )

  print *, 'allocated mesh'

  do i=0,NUM_PATCHES-1
     call mp2d%initialize_patch( i, NUM_CELLS1, NUM_CELLS2 )
  end do

  print *, 'initialized patches'

  do i=0,NUM_PATCHES-1
     print *, mp2d%get_delta_eta1(i)
  end do

  call sll_delete(mp2d)

  print *, 'PASSED'

end program unit_test_logical_meshes_multipatch
