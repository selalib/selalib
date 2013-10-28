program unit_test_1d
#include "sll_working_precision.h"
  use sll_constants
! The next line cause an error with the definition of 
! compute_interpolants in sll_cubic_spline_interpolator_1d
! use sll_module_mapped_meshes_2d
!  use geometry_functions
  use sll_common_coordinate_transformations
  use sll_cubic_spline_interpolator_1d
  use sll_module_mapped_meshes_1d
  implicit none

#define NPTS1 33

  type(sll_mapped_mesh_1d_analytic)     :: map_a    ! analytic map
  type(sll_mapped_mesh_1d_discrete)     :: map_d    ! discrete map
  ! for the discrete case...
  type(cubic_spline_1d_interpolator)    :: x1_interp
  type(cubic_spline_1d_interpolator)    :: j_interp
  sll_real64, dimension(:), allocatable :: x1
  sll_real64                            :: x1_eta1_min, x1_eta1_max
  sll_real64, dimension(:), allocatable :: jacs
  sll_int32  :: i
  sll_real64 :: eta1, h1, delta, acc, acc1
  sll_real64 :: node, node_a, node_d, interp, val_a
  sll_real64, dimension(2) :: params

  params(:) = (/-1.0_f64, 1.0_f64/)

  print *,  'filling out discrete arrays for x1 ', &
       'needed in the discrete case'
  h1 = 1.0_f64/real(NPTS1-1,f64)
  print *, 'h1 = ', h1
  allocate(x1(NPTS1))
  allocate(jacs(NPTS1))

  do i=0,NPTS1-1
     eta1      = real(i,f64)*h1
     x1(i+1)   = linear_map_f(eta1,params) 
     jacs(i+1) = linear_map_jac_f(eta1,params)
  end do
  
  ! Fill out the transformation's slopes at the borders
  eta1           = 0.0_f64
  ! In 1D the slope is the jacobian.
  x1_eta1_min = linear_map_jac_f(eta1,params)
  eta1           = 1.0_f64
  x1_eta1_max = linear_map_jac_f(eta1,params)

  print *, '**********************************************************'
  print *, '              TESTING THE ANALYTIC MAP                    '
  print *, '**********************************************************'

  ! Need to do something about these variables being always on the stack...
!  map_a => new_mapped_mesh_2D_general( ANALYTIC_MAP )

  !print *, x1_polar_f(1.0_f64,1.0_f64)
#ifdef STDF95
  call initialize( map_a, &
#else
  call map_a%initialize( &
#endif
       "map_a", &
       NPTS1, &
       linear_map_f, &
       linear_map_jac_f, &
       params )
  print *, 'initialized analytic map in 1D'

#ifdef STDF95
  print *, 'jacobian_1d(map_a, 0.5) = ', linear_map_jac_f(0.5_f64)
#else
  print *, 'jacobian_1d(map_a, 0.5) = ', map_a%jacobian(0.5_f64)
#endif

  acc  = 0.0_f64
  acc1 = 0.0_f64
  do i=0,NPTS1-1
     eta1    = real(i,f64)*h1
#ifdef STDF95
     node_a  = x1_at_node(map_a,i+1)
     val_a   = linear_map_f(eta1)
#else
     node_a  = map_a%x1_at_node(i+1)
     val_a   = map_a%x1(eta1)
#endif
     acc     = acc + abs(node_a-val_a)
  end do
  print *, 'Average error in nodes, x1 transformation = ', acc/(NPTS1)

#ifdef STDF95
  call write_to_file(map_a)
#else
  call map_a%write_to_file()
#endif


  print *, '**********************************************************'
  print *, '              TESTING THE DISCRETE MAP                    '
  print *, '**********************************************************'

  print *, 'initializing the interpolator: '

#ifdef STDF95
  call cubic_spline_1d_interpolator_initialize( x1_interp,&
#else
  call x1_interp%initialize( &
#endif
       NPTS1, &
       0.0_f64, &
       1.0_f64, &
       SLL_HERMITE, &
       slope_left=x1_eta1_min, &
       slope_right=x1_eta1_max )

#ifdef STDF95
  call cubic_spline_1d_interpolator_initialize( j_interp,&
#else
  call j_interp%initialize( &
#endif
       NPTS1, &
       0.0_f64, &
       1.0_f64, &
       SLL_HERMITE, &
       slope_left=linear_map_jac_f(0.0_f64,params), &
       slope_right=linear_map_jac_f(1.0_f64,params) )

#ifdef STDF95
  call initialize( map_d,&
#else
  call map_d%initialize( &
#endif
       "map_d", &
       NPTS1, &
       x1, &
       x1_interp, &
       j_interp, &
       jacobians_node=jacs )

 ! print *, 'x1: '
 ! print *, map_d%x1_node(:,:)

  print *, 'Compare the values of the transformation at the nodes: '
  acc  = 0.0_f64
  acc1 = 0.0_f64
  do i=1,NPTS1
#ifdef STDF95
     node_a   = x1_at_node(map_a,i)
     node_d   = x1_at_node(map_d,i)
#else
     node_a   = map_a%x1_at_node(i)
     node_d   = map_d%x1_at_node(i)
#endif
     acc      = acc + abs(node_a-node_d)
  end do
  print *, 'Average error in nodes, x1 transformation = ', acc/(NPTS1)

  print *, 'Compare the values of the jacobian at the nodes, resulting from ',&
       'calls to map_1d_jacobian_node() and jacobian_1D(map, eta1)'
  acc = 0.0_f64

  do i=0,NPTS1-1
     PRINT*,'i=',i
     eta1   = real(i,f64)*h1
#ifdef STDF95
     ! For analytic meshe the user give directly the jacobian
     ! In F95 standart qui can't write 
     !node   = mapped_meshes_jacobian(map_a,eta1)
     ! So call simply the user function
     node   = linear_map_jac_f(eta1)
     interp = jacobian(map_d,eta1)
#else
     !        print *, 'values: ', i, j, eta1, eta2
     !        print *, 'about to call map_a%jacobian(eta1,eta2)'
     node   = map_a%jacobian(eta1)
     !        node   = map_2d_jacobian_node(map_d,i+1,j+1)
     !        print *, 'about to call map_d%jacobian(eta1,eta2)'
     interp = map_d%jacobian(eta1)
#endif
     delta  =  node - interp
     print *, 'eta1 = ', eta1
     print *, '(',i+1,'): ANALYT = ', node, ', DISCR = ', interp, &
             '. DIFFERENCE  = ', delta
     !        print *, '(',i+1,j+1,'): NODE = ', node, ', ANALYT = ', jac_analyt, &
     !             '. DIFFERENCE  = ', delta2
     acc = acc + abs(delta)
  end do

#ifdef STDF95
  call write_to_file(map_d)
  call write_to_file(map_d)
#else
  call map_d%write_to_file()
  call map_d%write_to_file()
#endif
  print *, 'Average error = ', acc/real(NPTS1,f64)
  !  call delete(map_a)
  !  call delete(map_d)

  print *, 'deleted maps'
  print *, 'reached end of unit test'
end program unit_test_1d
