program unit_test
#include "sll_working_precision.h"
  use numeric_constants
  use sll_module_mapped_meshes_2d_cartesian
  implicit none

#define NPTS1 30
#define NPTS2 60 

  type(sll_mapped_mesh_2d_cartesian)    :: map_a    ! analytic map
  sll_int32  :: i, j
  sll_real64 :: eta1, eta2, h1, h2, delta, acc, acc1
  sll_real64 :: node, node_a, node_d, interp, val_a

  h1 = 1.0_f64/real(NPTS1-1,f64)
  h2 = 1.0_f64/real(NPTS2-1,f64)

  print *, '**********************************************************'
  print *, '              TESTING THE CARTESIAN MAP                    '
  print *, '**********************************************************'

  ! Need to do something about these variables being always on the stack...
!  map_a => new_mapped_mesh_2D_general( ANALYTIC_MAP )

  call map_a%initialize( &
       "map_cartesian", &
       NPTS1, &
       NPTS2 )

  print *, 'initialized analytic map'

  print *, 'jacobian_2d(map_a, 0.5, 0.5) = ', map_a%jacobian(0.5_f64,0.5_f64)

  acc  = 0.0_f64
  acc1 = 0.0_f64
  do j=0,NPTS2-1
     do i=0,NPTS1-1
        eta1    = real(i,f64)*h1
        eta2    = real(j,f64)*h2
        node_a  = map_a%x1_at_node(i+1,j+1)
        val_a   = map_a%x1(eta1,eta2)
        acc     = acc + abs(node_a-val_a)
        node_a  = map_a%x2_at_node(i+1,j+1)
        val_a   = map_a%x2(eta1,eta2)
        acc1    = acc1 + abs(node_a-val_a)
     end do
  end do
  print *, 'Average error in nodes, x1 transformation = ', acc/(NPTS1*NPTS2)
  print *, 'Average error in nodes, x2 transformation = ', acc1/(NPTS1*NPTS2)

  call map_a%write_to_file()

  PRINT*,'EXIT PROGRAM'
  call exit(1)
end program unit_test
