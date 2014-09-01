program unit_test_2d
#include "sll_working_precision.h"
  use sll_constants
  use sll_logical_meshes
  use sll_module_coordinate_transformations_2d
  use sll_module_coordinate_transformations_2d_nurbs
  use sll_common_coordinate_transformations
  use sll_cubic_spline_interpolator_2d
  
#include "sll_file_io.h"
  implicit none

#define NPTS1 33
#define NPTS2 33 
  type(sll_logical_mesh_2d), pointer :: mesh
  type(sll_coordinate_transformation_2d_analytic) :: t_a    ! analytic transf
  type(sll_coordinate_transformation_2d_discrete) :: t_d    ! discrete transf
  type(sll_coordinate_transformation_2d_nurbs)    :: t_n    ! nurbs transf
  type(sll_coordinate_transformation_2d_analytic), pointer :: t_a_ptr !test
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
  sll_real64, dimension(2) :: params   ! for the polar transformation
  sll_real64 :: val_approx1,val_approx2,val_exacte1,val_exacte2,val_exacte1_bis,val_exacte2_bis
  logical    :: l_exists
  sll_real64, dimension(2,2) :: val_approx_jac
  sll_real64 :: j11,j12,j21,j22,val_jac_approx
  sll_real64 :: j11_bis,j12_bis,j21_bis,j22_bis,jac,jac_bis
  sll_real64, dimension(4) :: param1
  sll_real64, dimension(10):: param2

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
  print *, '              TESTING THE ANALYTIC TRANSFORMATION '
  print *, '**********************************************************'


  mesh => new_logical_mesh_2d( NPTS1-1, NPTS2-1 )
  ! Need to do something about these variables being always on the stack...

  print *, x1_polar_f(1.0_f64,1.0_f64,params)
  call t_a%initialize( &
       "map_a", &
       mesh, &
       x1_polar_f, &
       x2_polar_f, &
       deriv_x1_polar_f_eta1, &
       deriv_x1_polar_f_eta2, &
       deriv_x2_polar_f_eta1, &
       deriv_x2_polar_f_eta2, &
       params)
  print *, 'initialized analytic map'

  ! The following pointer is not used but wanted to test the 'new' function
  ! wrapper.
  t_a_ptr => new_coordinate_transformation_2d_analytic( &
       "map_a", &
       mesh, &
       x1_polar_f, &
       x2_polar_f, &
       deriv_x1_polar_f_eta1, &
       deriv_x1_polar_f_eta2, &
       deriv_x2_polar_f_eta1, &
       deriv_x2_polar_f_eta2, &
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
  print *, 'Average error in nodes, x1 transformation = ', acc/(NPTS1*NPTS2)
  print *, 'Average error in nodes, x2 transformation = ', acc1/(NPTS1*NPTS2)

  call t_a%write_to_file()
  !call t_a%write_to_file(SLL_IO_MTV)

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
       SLL_HERMITE, &
       SLL_PERIODIC, &
       eta1_min_slopes=x1_eta1_min, &
       eta1_max_slopes=x1_eta1_max )

  call x2_interp%initialize( &
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

  call j_interp%initialize( &
       NPTS1, &
       NPTS2, &
       0.0_f64, &
       1.0_f64, &
       0.0_f64, &
       1.0_f64, &
       SLL_HERMITE, &
       SLL_PERIODIC, &
       const_eta1_min_slope=deriv1_jacobian_polar_f(0.0_f64,0.0_f64,params), &
       const_eta1_max_slope=deriv1_jacobian_polar_f(1.0_f64,0.0_f64,params) )

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
  print *, 'Average error in nodes, x1 transformation = ', acc/(NPTS1*NPTS2)
  print *, 'Average error in nodes, x2 transformation = ', acc1/(NPTS1*NPTS2)

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
  !call t_d%write_to_file(SLL_IO_MTV)

  print *, 'Average error in jacobian = ', acc/real(NPTS1*NPTS2,f64)
  call delete(t_a)
  call delete(t_d)

  print *, 'deleted transformations'

  ! *************************************************************************
  !
  ! Test of the initialization from a file of the discrete transformation.
  !
  ! *************************************************************************

  print *, 'Test of initialization from file for a nurbs transformation:'

  inquire(file="domain_patch0.nml", exist=l_exists)

  if (l_exists) then
     call t_n%read_from_file("domain_patch0.nml")
     !t_n%mesh => new_logical_mesh_2d(64,64 )
     call t_n%write_to_file()

     param1 = (/ 0.05_f64,0.05_f64,1.0_f64,1.0_f64 /)
     param2 = (/ 0.05_f64,0.05_f64,0.0_f64,1.0_f64,-1.0_f64,1.0_f64,0.0_f64,1.0_f64,-1.0_f64,1.0_f64 /)
     do i = 1,t_n%mesh%num_cells1+1
        do j = 1,t_n%mesh%num_cells2+1

           eta1 = t_n%mesh%eta1_min + (i-1) * t_n%mesh%delta_eta1 
           eta2 = t_n%mesh%eta2_min + (j-1) * t_n%mesh%delta_eta2
           val_approx1 = t_n%x1(eta1,eta2)
           val_approx2 = t_n%x2(eta1,eta2)
           val_exacte1 = 2*sinprod_x1(eta1,eta2,param1)-1
           val_exacte2 = 2*sinprod_x2(eta1,eta2,param1)-1
           val_exacte1_bis = sinprod_gen_x1(eta1,eta2,param2)
           val_exacte2_bis = sinprod_gen_x2(eta1,eta2,param2)

          ! print*, 'diff values composante eta1 ', abs(val_approx1-val_exacte1) , val_exacte1, val_exacte1_bis,val_approx1
          ! print*, 'diff values composante eta2 ', abs(val_approx2-val_exacte2) , val_exacte2, val_exacte2_bis,val_approx2 

           val_approx_jac = t_n%jacobian_matrix(eta1,eta2)
           val_jac_approx = t_n%jacobian(eta1,eta2)
           j11 = 2*sinprod_jac11(eta1,eta2,param1)
           j12 = 2*sinprod_jac12(eta1,eta2,param1)
           j21 = 2*sinprod_jac21(eta1,eta2,param1)
           j22 = 2*sinprod_jac22(eta1,eta2,param1)
           jac = 2*2*sinprod_jac(eta1,eta2,param1)
           j11_bis = sinprod_gen_jac11(eta1,eta2,param2)
           j12_bis = sinprod_gen_jac12(eta1,eta2,param2)
           j21_bis = sinprod_gen_jac21(eta1,eta2,param2)
           j22_bis = sinprod_gen_jac22(eta1,eta2,param2)
           jac_bis = sinprod_gen_jac  (eta1,eta2,param2)
          ! print*, 'Jacobian values composante j11 ', abs(val_approx_jac(1,1)-j11),j11,j11_bis, val_approx_jac(1,1) 
          ! print*, 'Jacobian values composante j12 ', abs(val_approx_jac(1,2)-j12),j12,j12_bis, val_approx_jac(1,2) 
          ! print*, 'Jacobian values composante j21 ', abs(val_approx_jac(2,1)-j21),j21,j21_bis, val_approx_jac(2,1) 
          ! print*, 'Jacobian values composante j22 ', abs(val_approx_jac(2,2)-j22),j22,j22_bis, val_approx_jac(2,2)
          ! print*, 'Jacobian values', abs(jac- val_jac_approx),jac, jac_bis, val_jac_approx
        end do
     end do

     print*, 'label t_n', t_n%label

     call delete(t_n)
     !call write_to_file(t_d)
  else
     print *, 'nml file is missing '
  end if
  print *, 'reached end of unit test'

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
