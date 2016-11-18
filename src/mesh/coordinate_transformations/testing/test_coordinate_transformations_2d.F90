program unit_test_2d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_hermite, &
    sll_p_periodic

  use sll_m_cartesian_meshes, only: &
    sll_f_new_cartesian_mesh_2d, &
    sll_t_cartesian_mesh_2d, &
    sll_o_delete

  use sll_m_common_coordinate_transformations, only: &
    sll_f_deriv1_jacobian_polar_f, &
    sll_f_deriv_x1_polar_f_eta1, &
    sll_f_deriv_x1_polar_f_eta2, &
    sll_f_deriv_x2_polar_f_eta1, &
    sll_f_deriv_x2_polar_f_eta2, &
    sll_f_jacobian_polar_f, &
    sll_f_sinprod_gen_jac, &
    sll_f_sinprod_gen_jac11, &
    sll_f_sinprod_gen_jac12, &
    sll_f_sinprod_gen_jac21, &
    sll_f_sinprod_gen_jac22, &
    sll_f_sinprod_gen_x1, &
    sll_f_sinprod_gen_x2, &
    sll_f_sinprod_jac, &
    sll_f_sinprod_jac11, &
    sll_f_sinprod_jac12, &
    sll_f_sinprod_jac21, &
    sll_f_sinprod_jac22, &
    sll_f_sinprod_x1, &
    sll_f_sinprod_x2, &
    sll_f_x1_polar_f, &
    sll_f_x2_polar_f

  use sll_m_coordinate_transformation_2d_base, only: &
    sll_p_io_mtv

  use sll_m_coordinate_transformations_2d, only: &
    sll_f_new_coordinate_transformation_2d_analytic, &
    sll_t_coordinate_transformation_2d_analytic, &
    sll_t_coordinate_transformation_2d_discrete, &
    sll_o_delete

  use sll_m_cubic_spline_interpolator_2d, only: &
    sll_t_cubic_spline_interpolator_2d

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#define NPTS1 33
#define NPTS2 33 
  type(sll_t_cartesian_mesh_2d), pointer :: mesh
  type(sll_t_coordinate_transformation_2d_analytic) :: t_a    ! analytic transf
  type(sll_t_coordinate_transformation_2d_discrete) :: t_d    ! discrete transf
  type(sll_t_coordinate_transformation_2d_analytic), pointer :: t_a_ptr !test
  ! for the discrete case...
  type(sll_t_cubic_spline_interpolator_2d)   :: x1_interp
  type(sll_t_cubic_spline_interpolator_2d)   :: x2_interp
  type(sll_t_cubic_spline_interpolator_2d)   :: j_interp
  sll_real64, dimension(:,:), allocatable :: x1_tab
  sll_real64, dimension(:,:), allocatable :: x2_tab
  sll_real64, dimension(:), allocatable   :: x1_eta1_min, x1_eta1_max
  sll_real64, dimension(:), allocatable   :: x2_eta1_min, x2_eta1_max
  sll_real64, dimension(:,:), allocatable :: jacs
  sll_int32  :: i, j
  sll_real64 :: eta1, eta2, h1, h2, delta, acc, acc1
  sll_real64 :: node, node_a, node_d, interp, val_a
  sll_real64, dimension(2) :: params   ! for the polar transformation

#define RMIN 0.1_f64
#define RMAX 1.0_f64

  params(:) = (/RMIN, RMAX/)

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
        x1_tab(i+1,j+1)   = sll_f_x1_polar_f(eta1,eta2,params) 
        x2_tab(i+1,j+1)   = sll_f_x2_polar_f(eta1,eta2,params) 
        jacs(i+1,j+1) = sll_f_jacobian_polar_f(eta1,eta2,params)
     end do
  end do

  ! Fill out the transformation's slopes at the borders
  do j=0,NPTS2-1
     eta1           = 0.0_f64
     eta2           = real(j,f64)*h2
     x1_eta1_min(j+1) = sll_f_deriv_x1_polar_f_eta1(eta1,eta2,params)
     x2_eta1_min(j+1) = sll_f_deriv_x2_polar_f_eta1(eta1,eta2,params)
     eta1           = 1.0_f64
     x1_eta1_max(j+1) = sll_f_deriv_x1_polar_f_eta1(eta1,eta2,params)
     x2_eta1_max(j+1) = sll_f_deriv_x2_polar_f_eta1(eta1,eta2,params)
  end do

  print *, '**********************************************************'
  print *, '              TESTING THE ANALYTIC TRANSFORMATION '
  print *, '**********************************************************'


  mesh => sll_f_new_cartesian_mesh_2d( NPTS1-1, NPTS2-1 )

  ! Need to do something about these variables being always on the stack...

  print *, sll_f_x1_polar_f(1.0_f64,1.0_f64,params)
  call t_a%initialize( &
       "map_a", &
       mesh, &
       sll_f_x1_polar_f, &
       sll_f_x2_polar_f, &
       sll_f_deriv_x1_polar_f_eta1, &
       sll_f_deriv_x1_polar_f_eta2, &
       sll_f_deriv_x2_polar_f_eta1, &
       sll_f_deriv_x2_polar_f_eta2, &
       params)
  print *, 'initialized analytic map'

  ! The following pointer is not used but wanted to test the 'new' function
  ! wrapper.
  t_a_ptr => sll_f_new_coordinate_transformation_2d_analytic( &
       "map_a", &
       mesh, &
       sll_f_x1_polar_f, &
       sll_f_x2_polar_f, &
       sll_f_deriv_x1_polar_f_eta1, &
       sll_f_deriv_x1_polar_f_eta2, &
       sll_f_deriv_x2_polar_f_eta1, &
       sll_f_deriv_x2_polar_f_eta2, &
       params )

  print *, 'jacobian_2d(t_a, 0.5, 0.5) = ', t_a%jacobian(0.5_f64,0.5_f64)

  acc  = 0.0_f64
  acc1 = 0.0_f64
  do j=0,NPTS2-1
     do i=0,NPTS1-1
        eta1    = real(i,f64)*h1
        eta2    = real(j,f64)*h2
        node_a  = t_a%x1_at_node(i+1,j+1)
        val_a   = t_a%x1(eta1,eta2)
        acc     = acc + abs(node_a-val_a)
        node_a  = t_a%x2_at_node(i+1,j+1)
        val_a   = t_a%x2(eta1,eta2)
        acc1    = acc1 + abs(node_a-val_a)
     end do
  end do
  print *, 'Average error in nodes, x1 transformation = ', acc /real(NPTS1*NPTS2,f64)
  print *, 'Average error in nodes, x2 transformation = ', acc1/real(NPTS1*NPTS2,f64)

  call t_a%write_to_file()
  !call t_a%write_to_file(sll_p_io_mtv)

  print *, '**********************************************************'
  print *, '              TESTING THE DISCRETE TRANSFORMATION         '
  print *, '**********************************************************'

  print *, 'initializing the interpolators: '

  call x1_interp%initialize( &
       NPTS1, &
       NPTS2, &
       0.0_f64, &
       1.0_f64, &      
       0.0_f64, &
       1.0_f64, &
       sll_p_hermite, &
       sll_p_periodic, &
       eta1_min_slopes=x1_eta1_min, &
       eta1_max_slopes=x1_eta1_max )

  call x2_interp%initialize( &
       NPTS1, &
       NPTS2, &
       0.0_f64, &
       1.0_f64, &
       0.0_f64, &
       1.0_f64, &
       sll_p_hermite, &
       sll_p_periodic, &
       eta1_min_slopes=x2_eta1_min, &
       eta1_max_slopes=x2_eta1_max )

  call j_interp%initialize( &
       NPTS1, &
       NPTS2, &
       0.0_f64, &
       1.0_f64, &
       0.0_f64, &
       1.0_f64, &
       sll_p_hermite, &
       sll_p_periodic, &
       const_eta1_min_slope=sll_f_deriv1_jacobian_polar_f(0.0_f64,0.0_f64,params), &
       const_eta1_max_slope=sll_f_deriv1_jacobian_polar_f(1.0_f64,0.0_f64,params) )

  print *, 'Initialized interpolators...'

  call t_d%initialize( &
       mesh, &
       "transf_d", &
       x1_interp, &
       x2_interp, &
       j_interp, &
       x1_tab, &
       x2_tab, &
       jacobians_node=jacs )

 ! print *, 'x1: '
 ! print *, map_d%x1_node(:,:)

  print *, 'Compare the values of the transformation at the nodes: '
  acc  = 0.0_f64
  acc1 = 0.0_f64

  do j=1,NPTS2
     do i=1,NPTS1
        node_a   = t_a%x1_at_node(i,j)
        node_d   = t_d%x1_at_node(i,j)
        acc      = acc + abs(node_a-node_d)
        node_a   = t_a%x2_at_node(i,j)
        node_d   = t_d%x2_at_node(i,j)
        acc1     = acc1 + abs(node_a-node_d)
     end do
  end do
  print *, 'Average error in nodes, x1 transformation = ', acc/real(NPTS1*NPTS2,f64)
  print *, 'Average error in nodes, x2 transformation = ', acc1/real(NPTS1*NPTS2,f64)

  print *, 'Compare the values of the jacobian at the nodes, resulting from ',&
       'calls to map_2d_jacobian_node() and jacobian_2D(map, eta1, eta2)'
  acc = 0.0_f64
  do j=0,NPTS2-1
     do i=0,NPTS1-1
        eta1   = real(i,f64)*h1
        eta2   = real(j,f64)*h2
!        print *, 'values: ', i, j, eta1, eta2
!        print *, 'about to call map_a%jacobian(eta1,eta2)'
        node   = t_a%jacobian(eta1,eta2)
!        node   = map_2d_jacobian_node(map_d,i+1,j+1)
!        print *, 'about to call map_d%jacobian(eta1,eta2)'
        interp = t_d%jacobian(eta1,eta2)
        delta  =  node - interp
        ! for inspecting/debugging:
!!$        print *, 'eta1 = ', eta1, 'eta2 = ', eta2
!!$        print *, '(',i+1,j+1,'): ANALYT = ', node, ', DISCR = ', interp, &
!!$             '. DIFFERENCE  = ', delta
!        print *, '(',i+1,j+1,'): NODE = ', node, ', ANALYT = ', jac_analyt, &
!             '. DIFFERENCE  = ', delta2
        acc = acc + abs(delta)
     end do
  end do

  call t_d%write_to_file()
  t_d%written = .false.
  call t_d%write_to_file(sll_p_io_mtv)

  print *, 'Average error in jacobian = ', acc/real(NPTS1*NPTS2,f64)
  call sll_o_delete(t_a)
  call sll_o_delete(t_d)

  print *, 'deleted transformations'

  ! apply some more relaxed criterion for the jacobian
  if( acc/real(NPTS1*NPTS2,f64) .lt. 3.0e-5 ) then
     print *, 'PASSED' 
  else
     print *, 'FAILED'
  end if

  deallocate(x1_eta1_min)
  deallocate(x1_eta1_max)
  deallocate(x2_eta1_min)
  deallocate(x2_eta1_max)




end program unit_test_2d
