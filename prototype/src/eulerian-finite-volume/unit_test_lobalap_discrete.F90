program test_lobalap_discrete
#include "selalib.h"
#include "sll_coordinate_transformations.h"

  use map_function_module, only: set_map_function
  use sll_lobatto_poisson
  use sll_dg_fields
  use sll_logical_meshes
  implicit none

  type(lobatto_poisson_solver)        :: solver
  type(sll_logical_mesh_2d), pointer  :: mesh
  class(sll_coordinate_transformation_2d_base), pointer :: tau
  type(dg_field), pointer :: dg_rho
  type(dg_field), pointer :: dg_ex
  type(dg_field), pointer :: dg_ey

  sll_int32, parameter :: degree = 3
  real(8), external :: f_cos, f_four

  type(cubic_spline_2d_interpolator)      :: x1_interp
  type(cubic_spline_2d_interpolator)      :: x2_interp
  type(cubic_spline_2d_interpolator)      :: j_interp
  sll_real64, dimension(:,:), allocatable :: x1_tab
  sll_real64, dimension(:,:), allocatable :: x2_tab
  sll_real64, dimension(:), allocatable   :: x1_eta1_min, x1_eta1_max
  sll_real64, dimension(:), allocatable   :: x2_eta1_min, x2_eta1_max
  sll_real64, dimension(:,:), allocatable :: jacs
  
#define NPTS1 33
#define NPTS2 33
#define R_MIN  0.1_8
#define R_MAX  1.0_8
#define THETA_MIN  0.00_8
#define THETA_MAX  1.00_8
#define N 6

  sll_int32  :: i, j
  sll_real64 :: eta1, eta2, h1, h2

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
  
  mesh => new_logical_mesh_2d( NPTS1-1, NPTS2-1 )

  do j=0,NPTS2-1
     do i=0,NPTS1-1
        eta1            = real(i,f64)*h1
        eta2            = real(j,f64)*h2
        x1_tab(i+1,j+1) = x1_polar_f(eta1,eta2,[R_MIN,R_MAX]) 
        x2_tab(i+1,j+1) = x2_polar_f(eta1,eta2,[R_MIN,R_MAX]) 
        jacs(i+1,j+1)   = jacobian_polar_f(eta1,eta2,[R_MIN,R_MAX])
     end do
  end do

  do j=0,NPTS2-1
     eta1           = 0.0_f64
     eta2           = real(j,f64)*h2
     x1_eta1_min(j+1) = deriv_x1_polar_f_eta1(eta1,eta2,[R_MIN,R_MAX])
     x2_eta1_min(j+1) = deriv_x2_polar_f_eta1(eta1,eta2,[R_MIN,R_MAX])
     eta1           = 1.0_f64
     x1_eta1_max(j+1) = deriv_x1_polar_f_eta1(eta1,eta2,[R_MIN,R_MAX])
     x2_eta1_max(j+1) = deriv_x2_polar_f_eta1(eta1,eta2,[R_MIN,R_MAX])
  end do

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
       const_eta1_min_slope=deriv1_jacobian_polar_f(0.0_f64,0.0_f64,[R_MIN,R_MAX]), &
       const_eta1_max_slope=deriv1_jacobian_polar_f(1.0_f64,0.0_f64,[R_MIN,R_MAX]) )

  print *, 'Initialized interpolators...'

  tau => new_coordinate_transformation_2d_discrete( &
       mesh, &
       "polar_discrete", &
       x1_interp, &
       x2_interp, &
       j_interp, &
       x1_tab, &
       x2_tab, &
       jacobians_node=jacs )

  call tau%write_to_file()

  dg_rho => new_dg_field( degree, tau, f_four ) 
  dg_ex => new_dg_field( degree, tau ) 
  dg_ey => new_dg_field( degree, tau ) 

  call dg_rho%write_to_file('rho')

  call initialize(solver, tau, degree )
  call solve(solver, dg_rho, dg_ex, dg_ey)
  call delete(solver)

  call dg_ex%write_to_file('ex')
  call dg_ey%write_to_file('ey')

end program test_lobalap_discrete

real(8) function f_cos( r, theta )

   real(8) :: r
   real(8) :: theta

   f_cos = (r-R_MIN)*(r-R_MAX)*cos(N*theta)*r

end function f_cos

real(8) function f_four( x, y)
   real(8) :: x, y
   f_four = -4._8

end function f_four

real(8) function f_sin( r, theta )

   real(8) :: r
   real(8) :: theta

   f_sin = (r-R_MIN)*(r-R_MAX)*sin(N*theta)*r

end function f_sin

real(8) function lap_f_cos( r, theta )

   !sage: assume(r>=1)
   !sage: assume(r<=2)
   !sage: phi = (r-R_MIN)*(r-R_MAX)*r*cos(n*theta)
   !sage: diff(r*diff(phi,r),r)/r + diff(phi,theta,theta)/(r*r)

   real(8) :: r
   real(8) :: theta

   lap_f_cos = -(r-R_MAX)*(r-R_MIN)*N*N*cos(N*theta)/r &
           + ((r-R_MAX)*(r-R_MIN)*cos(N*theta)  &
           + (r-R_MAX)*r*cos(N*theta) + (r-R_MIN)*r*cos(N*theta) &
           + 2*((r-R_MAX)*cos(N*theta) + (r-R_MIN)*cos(N*theta) &
           + r*cos(N*theta))*r)/r


end function lap_f_cos

real(8) function lap_f_sin( r, theta)

   !sage: assume(r>=1)
   !sage: assume(r<=2)
   !sage: phi = (r-R_MIN)*(r-R_MAX)*r*sin(n*theta)
   !sage: diff(r*diff(phi,r),r)/r + diff(phi,theta,theta)/(r*r)

   real(8) :: r
   real(8) :: theta
   
   lap_f_sin = -(r-R_MAX)*(r-R_MIN)*N*N*sin(N*theta)/r &
         + ((r-R_MAX)*(r-R_MIN)*sin(N*theta) &
         + (r-R_MAX)*r*sin(N*theta) + (r-R_MIN)*r*sin(N*theta) &
         + 2*((r-R_MAX)*sin(N*theta) + (r-R_MIN)*sin(N*theta)  &
         + r*sin(N*theta))*r)/r

end function lap_f_sin
