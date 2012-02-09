program unit_test
#include "sll_working_precision.h"
  use sll_coordinate_transformation
  use geometry_functions
  implicit none

#define NPTS1 32
#define NPTS2 32

  type(map_2D), pointer     :: map_a    ! analytic map
  type(map_2D), pointer     :: map_d    ! discrete map
  sll_real64, dimension(:,:), allocatable :: x1
  sll_real64, dimension(:,:), allocatable :: x2
  sll_real64, dimension(:,:), allocatable :: jacs
  sll_int32  :: i, j
  sll_real64 :: eta1, eta2, h1, h2, delta, acc, node, interp

  print *,  'filling out discrete arrays for x1 and x2 ', &
       'needed in the discrete case'
  h1 = 1.0_f64/real(NPTS1-1,f64)
  h2 = 1.0_f64/real(NPTS2-1,f64)
  allocate(x1(NPTS1,NPTS2))
  allocate(x2(NPTS1,NPTS2))
  allocate(jacs(NPTS1,NPTS2))
  do j=0,NPTS2-1
     do i=0,NPTS1-1
        eta1          = real(i,f64)*h1
        eta2          = real(j,f64)*h2
        x1(i+1,j+1)   = eta1*cos(eta2)
        x2(i+1,j+1)   = eta1*sin(eta2)
        jacs(i+1,j+1) = eta1
     end do
  end do

  print *, '**********************************************************'
  print *, '              TESTING THE ANALYTIC MAP                    '
  print *, '**********************************************************'

  map_a => new_map_2D( ANALYTIC_MAP )
  print *, 'allocated map'

  call initialize_map_2D( &
       map_a, &
       NPTS1, &
       NPTS2, &
       x1_func=polar_x1, &
       x2_func=polar_x2, &
       j11_func=polar_jac11, &
       j12_func=polar_jac12, &
       j21_func=polar_jac21, &
       j22_func=polar_jac22 )
  print *, 'initialized map'

  print *, 'jacobian_2d(map_a, 0.5, 0.5) = ', jacobian_2d(map_a,0.5_f64,0.5_f64)

  print *, '**********************************************************'
  print *, '              TESTING THE DISCRETE MAP                    '
  print *, '**********************************************************'

  map_d => new_map_2D( DISCRETE_MAP )
  print *, 'allocated discrete map'

  call initialize_map_2D( &
       map_d, &
       NPTS1, &
       NPTS2, &
       x1_node=x1, &
       x2_node=x2, &
       jacobians_node=jacs, &
       eta1_bc_type_x1=HERMITE_MAP_BC, &
       eta2_bc_type_x1=PERIODIC_MAP_BC, &
       eta1_bc_type_x2=HERMITE_MAP_BC, &
       eta2_bc_type_x2=PERIODIC_MAP_BC )

  print *, 'Compare the values of the jacobian at the nodes, resulting from ',&
       'calls to map_2d_jacobian_node() and jacobian_2D(map, eta1, eta2)'
  acc = 0.0_f64
  do j=0,NPTS1-1
     do i=0,NPTS2-1
        eta1   = real(i,f64)*h1
        eta2   = real(j,f64)*h2
        node   = map_2d_jacobian_node(map_d,i+1,j+1)
        interp = jacobian_2D(map_d,eta1,eta2) 
        delta  =  node - interp
             
        print *, 'NODE = ', node, ', INTERP = ', interp, &
             '. DIFFERENCE at (', i,j,') = ', delta
        acc = acc + abs(delta)
     end do
  end do

  print *, 'Average error = ', acc/real(NPTS1*NPTS2,f64)
  call delete(map_a)
  call delete(map_d)
  print *, 'deleted maps'
  print *, 'reached end of unit test'
end program unit_test
