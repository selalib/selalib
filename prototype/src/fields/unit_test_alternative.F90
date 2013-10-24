program unit_test_alternative
#include "sll_working_precision.h"
#include "sll_memory.h"
  use sll_logical_meshes
  use sll_constants
  use sll_module_scalar_field_2d_alternative
  use sll_module_coordinate_transformations_2d
  use sll_common_coordinate_transformations
  implicit none

#define SPLINE_DEG1 3
#define SPLINE_DEG2 3
#define NUM_CELLS1  64
#define NUM_CELLS2  64
#define ETA1MIN  0.0_f64
#define ETA1MAX  1.0_f64
#define ETA2MIN  0.0_f64
#define ETA2MAX  1.0_f64
  
  type(sll_logical_mesh_2d), pointer               :: mesh_2d
  class(sll_coordinate_transformation_2d_base), pointer :: T
  ! either of these type declarations can be used to work. Initialization is
  ! different.
  class(sll_scalar_field_2d_base), pointer              :: field_2d
  type(sll_scalar_field_2d_analytic_alt)                :: field_2d_a
  real(8), external :: test_function
  sll_int32 :: nc1, nc2, iplot
  ! procedure(polar_x1), pointer :: px1, px2, pjac11, pjac12, pjac21, pjac22
  ! type(init_landau_2d), target :: init_landau
  ! class(scalar_field_2d_initializer_base), pointer    :: pfinit
  ! type(cubic_spline_1d_interpolator), target  :: interp_eta1
  ! type(cubic_spline_1d_interpolator), target  :: interp_eta2
  ! class(sll_interpolator_1d_base), pointer :: interp_eta1_ptr
  ! class(sll_interpolator_1d_base), pointer :: interp_eta2_ptr
  type(arb_deg_2d_interpolator), target                 :: interp_2d_term_source
  class(sll_scalar_field_2d_base), pointer              :: rho
  sll_real64, dimension(:,:), allocatable    :: calculated
  sll_real64, dimension(:,:), allocatable    :: difference
  sll_real64, dimension(:,:), allocatable    :: tab_rho
  sll_real64, dimension(:),   allocatable    :: point1
  sll_real64, dimension(:),   allocatable    :: point2
  sll_real64 :: eta1,eta2
  sll_real64  :: h1,h2
  sll_int32 :: npts1,npts2
  sll_int32 :: i,j
  sll_int32 :: ierr
  sll_real64, dimension(1) :: params_identity

  params_identity(:) = (/ 0.0_f64 /)

  ! logical mesh
  nc1 = 32
  nc2 = 32
  
  mesh_2d => new_logical_mesh_2d( nc1, nc2, &
       0.0_f64, 2.0*sll_pi, 0.0_f64,2.0*sll_pi )
  
  print *, 'initialized mesh 2D'

  ! coordinate transformation
  T => new_coordinate_transformation_2d_analytic( &
       "analytic", &
       mesh_2d, &
       identity_x1, &
       identity_x2, &
       identity_jac11, &
       identity_jac12, &
       identity_jac21, &
       identity_jac22, &
       params_identity )
  print *, 'initialized transformation'
  call field_2d_a%initialize( &
       test_function, &
       'doubly_periodic', &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC, &
       SLL_PERIODIC )
  print *, 'initialized field 2d'

  print *, 'field value at 0,0 = ', field_2d_a%value_at_point(0.0_f64,0.0_f64)
  print *, 'field value at indices 1,1 = ', &
       field_2d_a%value_at_indices(1,1)

  call field_2d_a%write_to_file(0)

   ! the following call can also be made as:
  ! call field_2d_a%delete()
  ! we leave this as follows to test if any compilers complain about this
  ! syntax.
  call delete(field_2d_a)

!!! --------------------------------------------------------------------------
!   Test case periodic-periodic non analytic
!----------------------------------------------------------------------------

  npts1 =  NUM_CELLS1 + 1
  npts2 =  NUM_CELLS2 + 1
  h1 = (ETA1MAX - ETA1MIN)/real(NPTS1-1,f64)
  h2 = (ETA2MAX - ETA2MIN)/real(NPTS2-1,f64)
  print *, 'h1 = ', h1
  print *, 'h2 = ', h2
  
  ! Table to represent the node values of phi
  SLL_ALLOCATE(calculated(npts1,npts2),ierr)
  SLL_ALLOCATE(difference(npts1,npts2),ierr)

  allocate(point1(npts1-1))
  allocate(point2(npts2-1))
  allocate(tab_rho(npts1-1,npts2-1))
  do j=0,npts2-2
     do i=0,npts1-2
        point1(i+1)       = real(i,f64)*(ETA1MAX-ETA1MIN)/(npts1-1) + ETA1MIN 
        point2(j+1)       = real(j,f64)*(ETA2MAX-ETA2MIN)/(npts2-1) + ETA2MIN 
        tab_rho(i+1,j+1)  = cos(2.0_f64*sll_pi*point2(j+1))*cos(2.0_f64*sll_pi*point1(i+1))
     end do
  end do

  !print*, tab_rho
  call initialize_ad2d_interpolator( &
       interp_2d_term_source, &
       NUM_CELLS1+1, &
       NUM_CELLS2+1, &
       ETA1MIN, &
       ETA1MAX, &
       ETA2MIN, &
       ETA2MAX, &
       SLL_PERIODIC, &
       SLL_PERIODIC,&
       SLL_PERIODIC,&
       SLL_PERIODIC,&
       SPLINE_DEG1, &
       SPLINE_DEG2)
  
 ! terme_source_interp => interp_2d_term_source

  rho => new_scalar_field_2d_discrete_alt( &
       tab_rho, &
       "rho", &
       interp_2d_term_source, &
       T, &
       SLL_PERIODIC, &
       SLL_PERIODIC,&
       SLL_PERIODIC,&
       SLL_PERIODIC,&
       point1,&
       npts1-1,&
       point2,&
       npts2-1)

   do j=0,npts2-2
     do i=0,npts1-2
        eta1 = real(i,f64)*(ETA1MAX-ETA1MIN)/(2*(npts1-1)) + ETA1MIN 
        eta2 = real(j,f64)*(ETA2MAX-ETA2MIN)/(2*(npts2-1)) + ETA2MIN
       calculated(i+1,j+1) = rho%value_at_point(eta1,eta2)
       difference(i+1,j+1) = calculated(i+1,j+1)-cos(2.0_f64*sll_pi*eta2)*cos(2.0_f64*sll_pi*eta1)
!       print*, 'point=',eta1,eta2,'difference=', difference(i+1,j+1), calculated(i+1,j+1),cos(2.0_f64*sll_pi*eta2)*cos(2.0_f64*sll_pi*eta1)
     end do
  end do

  call rho%write_to_file(0)


 

  print *, 'PASSED'
  
end program unit_test_alternative


function test_function( eta1, eta2) result(res)
  intrinsic :: cos
  real(8) :: res
  real(8), intent(in) :: eta1
  real(8), intent(in) :: eta2
  !real(8), dimension(:), intent(in), optional :: params
  res = 2.0*cos(eta1)*cos(eta2)
end function test_function


