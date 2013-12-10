program unit_test
#include "sll_working_precision.h"
  use sll_constants
  use sll_module_mapped_meshes_2d
  use sll_common_coordinate_transformations
  use sll_cubic_spline_interpolator_2d
  implicit none

#define NPTS1 33
#define NPTS2 33 

  type(sll_mapped_mesh_2d_analytic)    :: map_a    ! analytic map
  type(sll_mapped_mesh_2d_discrete)    :: map_d    ! discrete map
  type(sll_mapped_mesh_2d_analytic), pointer :: map_a_ptr !test
  ! for the discrete case...
  type(cubic_spline_2d_interpolator)   :: x1_interp
  type(cubic_spline_2d_interpolator)   :: x2_interp
  type(cubic_spline_2d_interpolator)   :: j_interp
  sll_real64, dimension(:,:), allocatable :: x1_tab
  sll_real64, dimension(:,:), allocatable :: x2_tab
  sll_real64, dimension(:), allocatable   :: x1_eta1_min, x1_eta1_max
  sll_real64, dimension(:), allocatable   :: x2_eta1_min, x2_eta1_max
  sll_real64, dimension(:,:), allocatable :: jacs
  sll_int32  :: i, j
  sll_real64 :: eta1, eta2, h1, h2, delta, acc, acc1
  sll_real64 :: node, node_a, node_d, interp, val_a
  sll_real64, dimension(2) :: params

  params(:) = (/0.1_f64, 1.0_f64/)

  print *,  'filling out discrete arrays for x1 and x2 ', &
       'needed in the discrete case'
  h1 = 1.0_f64/real(NPTS1-1,f64)
  h2 = 1.0_f64/real(NPTS2-1,f64)
  print *, 'h1 = ', h1
  print *, 'h2 = ', h2
  allocate(x1_tab(NPTS1,NPTS2))
  allocate(x2_tab(NPTS1,NPTS2))
  allocate(x1_eta1_min(NPTS2))
  allocate(x1_eta1_max(NPTS2))
  allocate(x2_eta1_min(NPTS2))
  allocate(x2_eta1_max(NPTS2))
  allocate(jacs(NPTS1,NPTS2))

  do j=0,NPTS2-1
     do i=0,NPTS1-1
        eta1          = real(i,f64)*h1
        eta2          = real(j,f64)*h2
        x1_tab(i+1,j+1)   = x1_polar_f(eta1,eta2,params) 
        x2_tab(i+1,j+1)   = x2_polar_f(eta1,eta2,params) 
        jacs(i+1,j+1) = jacobian_polar_f(eta1,eta2,params)
     end do
  end do

  ! Fill out the transformation's slopes at the borders
  do j=0,NPTS2-1
     eta1           = 0.0_f64
     eta2           = real(j,f64)*h2
     x1_eta1_min(j+1) = deriv_x1_polar_f_eta1(eta1,eta2,params)
     x2_eta1_min(j+1) = deriv_x2_polar_f_eta1(eta1,eta2,params)
     eta1           = 1.0_f64
     x1_eta1_max(j+1) = deriv_x1_polar_f_eta1(eta1,eta2,params)
     x2_eta1_max(j+1) = deriv_x2_polar_f_eta1(eta1,eta2,params)
  end do

  print *, '**********************************************************'
  print *, '              TESTING THE ANALYTIC MAP                    '
  print *, '**********************************************************'

  ! Need to do something about these variables being always on the stack...
!  map_a => new_mapped_mesh_2D_general( ANALYTIC_MAP )

  print *, x1_polar_f(1.0_f64,1.0_f64,params)
#ifdef STDF95
  call initialize( map_a,&
#else
  call map_a%initialize( &
#endif
       "map_a", &
       NPTS1, &
       NPTS2, &
       x1_polar_f, &
       x2_polar_f, &
       deriv_x1_polar_f_eta1, &
       deriv_x1_polar_f_eta2, &
       deriv_x2_polar_f_eta1, &
       deriv_x2_polar_f_eta2, &
       params )
  print *, 'initialized analytic map'

  ! The following pointer is not used but wanted to test the 'new' function
  ! wrapper.
  map_a_ptr => new_mesh_2d_analytic( &
       "map_a", &
       NPTS1, &
       NPTS2, &
       x1_polar_f, &
       x2_polar_f, &
       deriv_x1_polar_f_eta1, &
       deriv_x1_polar_f_eta2, &
       deriv_x2_polar_f_eta1, &
       deriv_x2_polar_f_eta2, &
       params )

#ifdef STDF95
  print *, 'jacobian_2d(map_a, 0.5, 0.5) = ', &
       deriv_x1_polar_f_eta1(0.5_f64,0.5_f64,params)*&
       deriv_x2_polar_f_eta2(0.5_f64,0.5_f64,params)-&
       deriv_x1_polar_f_eta2(0.5_f64,0.5_f64,params)*&
       deriv_x2_polar_f_eta1(0.5_f64,0.5_f64,params)
#else
  print *, 'jacobian_2d(map_a, 0.5, 0.5) = ', map_a%jacobian(0.5_f64,0.5_f64)
#endif

  acc  = 0.0_f64
  acc1 = 0.0_f64
  do j=0,NPTS2-1
     do i=0,NPTS1-1
        eta1    = real(i,f64)*h1
        eta2    = real(j,f64)*h2
#ifdef STDF95
        node_a  = x1_at_node(map_a,i+1,j+1)
        val_a   = x1_polar_f(eta1,eta2,params)
#else
        node_a  = map_a%x1_at_node(i+1,j+1)
        val_a   = map_a%x1(eta1,eta2)
#endif
        acc     = acc + abs(node_a-val_a)
#ifdef STDF95
        node_a  = x2_at_node(map_a,i+1,j+1)
        val_a   = x2_polar_f(eta1,eta2,params)
#else
        node_a  = map_a%x2_at_node(i+1,j+1)
        val_a   = map_a%x2(eta1,eta2)
#endif
        acc1    = acc1 + abs(node_a-val_a)
     end do
  end do
  print *, 'Average error in nodes, x1 transformation = ', acc/(NPTS1*NPTS2)
  print *, 'Average error in nodes, x2 transformation = ', acc1/(NPTS1*NPTS2)

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
  call cubic_spline_2d_initialize( x1_interp,&
#else
  call x1_interp%initialize( &
#endif
       NPTS1, &
       NPTS2, &
       0.0_f64, &
       1.0_f64, &
       0.0_f64, &
       1.0_f64, &
       SLL_HERMITE, &
       SLL_PERIODIC, &
       eta1_min_slopes=x1_eta1_min, &
       eta1_max_slopes=x1_eta1_max )

#ifdef STDF95
  call cubic_spline_2d_initialize( x2_interp,&
#else
  call x2_interp%initialize( &
#endif
       NPTS1, &
       NPTS2, &
       0.0_f64, &
       1.0_f64, &
       0.0_f64, &
       1.0_f64, &
       SLL_HERMITE, &
       SLL_PERIODIC, &
       eta1_min_slopes=x2_eta1_min, &
       eta1_max_slopes=x2_eta1_max )

#ifdef STDF95
  call cubic_spline_2d_initialize( j_interp,&
#else
  call j_interp%initialize( &
#endif
       NPTS1, &
       NPTS2, &
       0.0_f64, &
       1.0_f64, &
       0.0_f64, &
       1.0_f64, &
       SLL_HERMITE, &
       SLL_PERIODIC, &
       const_eta1_min_slope=deriv1_jacobian_polar_f(0.0_f64,0.0_f64,params), &
       const_eta1_max_slope=deriv1_jacobian_polar_f(0.0_f64,0.0_f64,params) )

#ifdef STDF95
  call initialize( map_d,&
#else
  call map_d%initialize( &
#endif
       "map_d", &
       NPTS1, &
       NPTS2, &
       x1_tab, &
       x2_tab, &
       x1_interp, &
       x2_interp, &
       j_interp, &
       jacobians_node=jacs )

 ! print *, 'x1: '
 ! print *, map_d%x1_node(:,:)

  print *, 'Compare the values of the transformation at the nodes: '
  acc  = 0.0_f64
  acc1 = 0.0_f64
  do j=1,NPTS2
     do i=1,NPTS1
#ifdef STDF95
        node_a   = x1_at_node(map_a,i,j)
        node_d   = x1_at_node(map_d,i,j)
#else
        node_a   = map_a%x1_at_node(i,j)
        node_d   = map_d%x1_at_node(i,j)
#endif
        acc      = acc + abs(node_a-node_d)
#ifdef STDF95
        node_a   = x2_at_node(map_a,i,j)
        node_d   = x2_at_node(map_d,i,j)
#else
        node_a   = map_a%x2_at_node(i,j)
        node_d   = map_d%x2_at_node(i,j)
#endif
        acc1     = acc1 + abs(node_a-node_d)
     end do
  end do
  print *, 'Average error in nodes, x1 transformation = ', acc/(NPTS1*NPTS2)
  print *, 'Average error in nodes, x2 transformation = ', acc1/(NPTS1*NPTS2)

  print *, 'Compare the values of the jacobian at the nodes, resulting from ',&
       'calls to map_2d_jacobian_node() and jacobian_2D(map, eta1, eta2)'
  acc = 0.0_f64
  do j=0,NPTS2-1
     do i=0,NPTS1-1
        eta1   = real(i,f64)*h1
        eta2   = real(j,f64)*h2
#ifdef STDF95
        node   =  deriv_x1_polar_f_eta1(eta1,eta2,params)*&
                  deriv_x2_polar_f_eta2(eta1,eta2,params)-&
                  deriv_x1_polar_f_eta2(eta1,eta2,params)*&
                  deriv_x2_polar_f_eta1(eta1,eta2,params)
        interp = jacobian(map_d,eta1,eta2) 
#else
!        print *, 'values: ', i, j, eta1, eta2
!        print *, 'about to call map_a%jacobian(eta1,eta2)'
        node   = map_a%jacobian(eta1,eta2)
!        node   = map_2d_jacobian_node(map_d,i+1,j+1)
!        print *, 'about to call map_d%jacobian(eta1,eta2)'
        interp = map_d%jacobian(eta1,eta2)
#endif
        delta  =  node - interp
        print *, 'eta1 = ', eta1, 'eta2 = ', eta2
        print *, '(',i+1,j+1,'): ANALYT = ', node, ', DISCR = ', interp, &
             '. DIFFERENCE  = ', delta
!        print *, '(',i+1,j+1,'): NODE = ', node, ', ANALYT = ', jac_analyt, &
!             '. DIFFERENCE  = ', delta2
        acc = acc + abs(delta)
     end do
  end do

#ifdef STDF95
  call write_to_file(map_d)
  call write_to_file(map_d)
#else
  call map_d%write_to_file()
  call map_d%write_to_file()
#endif

  print *, 'Average error = ', acc/real(NPTS1*NPTS2,f64)
!  call delete(map_a)
!  call delete(map_d)

  print *, 'deleted maps'
  print *, 'reached end of unit test'

  deallocate(x1_eta1_min)
  deallocate(x1_eta1_max)
  deallocate(x2_eta1_min)
  deallocate(x2_eta1_max)

end program unit_test
