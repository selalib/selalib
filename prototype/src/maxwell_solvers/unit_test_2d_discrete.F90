program test_maxwell_2d_discrete
!--------------------------------------------------------------------------
!  test 2D Maxwell solver based on discontinuous galerkine on a mapped mesh
!  with discrete coordinate transformation
!--------------------------------------------------------------------------
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_constants.h"
#include "sll_maxwell_solvers_macros.h"
#include "sll_file_io.h"

use sll_logical_meshes
use sll_module_coordinate_transformations_2d
use sll_common_coordinate_transformations
use sll_cubic_spline_interpolator_2d

implicit none
  
#define NPTS1 33
#define NPTS2 33 
  type(sll_logical_mesh_2d), pointer :: mesh
  type(sll_coordinate_transformation_2d_discrete) :: t_d    ! discrete transf
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
  sll_real64 :: eta1, eta2, h1, h2
  sll_real64, dimension(2) :: params   ! for the polar transformation

#define RMIN 0.1_f64
#define RMAX 1.0_f64

  params(:) = (/RMIN, RMAX/)

  print*,'filling out discrete arrays for x1 and x2 '
  print*,'needed in the discrete case'

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
        eta1            = real(i,f64)*h1
        eta2            =  real(j,f64)*h2
        x1_tab(i+1,j+1) = x1_polar_f(eta1,eta2,params) 
        x2_tab(i+1,j+1) = x2_polar_f(eta1,eta2,params) 
        jacs(i+1,j+1)   = jacobian_polar_f(eta1,eta2,params)
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

  call t_d%write_to_file()

  call delete(t_d)

  deallocate(x1_eta1_min)
  deallocate(x1_eta1_max)
  deallocate(x2_eta1_min)
  deallocate(x2_eta1_max)

end program test_maxwell_2d_discrete
