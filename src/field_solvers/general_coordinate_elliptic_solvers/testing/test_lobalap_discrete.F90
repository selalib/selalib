program test_lobalap_discrete
!YG #include "selalib.h"
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_hermite, &
    sll_p_periodic

  use sll_m_cartesian_meshes, only: &
    sll_t_cartesian_mesh_2d, &
    sll_o_delete, &
    sll_o_new

  use sll_m_common_coordinate_transformations, only: &
    sll_f_deriv1_jacobian_polar_f, &
    sll_f_deriv_x1_polar_f_eta1, &
    sll_f_deriv_x2_polar_f_eta1, &
    sll_f_jacobian_polar_f, &
    sll_f_x1_polar_f, &
    sll_f_x2_polar_f

  use sll_m_coordinate_transformation_2d_base, only: &
    sll_c_coordinate_transformation_2d_base

  use sll_m_coordinate_transformations_2d, only: &
    sll_f_new_coordinate_transformation_2d_discrete

  use sll_m_cubic_spline_interpolator_2d, only: &
    sll_t_cubic_spline_interpolator_2d

  use sll_m_dg_fields, only: &
    sll_t_dg_field_2d, &
    sll_o_new

  use sll_m_lobatto_poisson, only: &
    sll_t_lobatto_poisson_solver, &
    sll_o_create, &
    sll_o_solve, &
    sll_o_delete

  use sll_m_map_function, only: &
    sll_s_set_map_function

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type(sll_t_lobatto_poisson_solver)          :: solver
  type(sll_t_cartesian_mesh_2d), pointer  :: mesh
  class(sll_c_coordinate_transformation_2d_base), pointer :: tau
  type(sll_t_dg_field_2d), pointer :: dg_rho
  type(sll_t_dg_field_2d), pointer :: dg_ex
  type(sll_t_dg_field_2d), pointer :: dg_ey

  sll_int32, parameter :: degree = 3

  type(sll_t_cubic_spline_interpolator_2d)  :: x1_interp
  type(sll_t_cubic_spline_interpolator_2d)  :: x2_interp
  type(sll_t_cubic_spline_interpolator_2d)  :: j_interp
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
  
  mesh => sll_o_new( NPTS1-1, NPTS2-1 )

  do j=0,NPTS2-1
     do i=0,NPTS1-1
        eta1            = real(i,f64)*h1
        eta2            = real(j,f64)*h2
        x1_tab(i+1,j+1) = sll_f_x1_polar_f(eta1,eta2,[R_MIN,R_MAX]) 
        x2_tab(i+1,j+1) = sll_f_x2_polar_f(eta1,eta2,[R_MIN,R_MAX]) 
        jacs(i+1,j+1)   = sll_f_jacobian_polar_f(eta1,eta2,[R_MIN,R_MAX])
     end do
  end do

  do j=0,NPTS2-1
     eta1           = 0.0_f64
     eta2           = real(j,f64)*h2
     x1_eta1_min(j+1) = sll_f_deriv_x1_polar_f_eta1(eta1,eta2,[R_MIN,R_MAX])
     x2_eta1_min(j+1) = sll_f_deriv_x2_polar_f_eta1(eta1,eta2,[R_MIN,R_MAX])
     eta1           = 1.0_f64
     x1_eta1_max(j+1) = sll_f_deriv_x1_polar_f_eta1(eta1,eta2,[R_MIN,R_MAX])
     x2_eta1_max(j+1) = sll_f_deriv_x2_polar_f_eta1(eta1,eta2,[R_MIN,R_MAX])
  end do

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
       const_eta1_min_slope=sll_f_deriv1_jacobian_polar_f(0.0_f64,0.0_f64,[R_MIN,R_MAX]), &
       const_eta1_max_slope=sll_f_deriv1_jacobian_polar_f(1.0_f64,0.0_f64,[R_MIN,R_MAX]) )

  print *, 'Initialized interpolators...'

  tau => sll_f_new_coordinate_transformation_2d_discrete( &
       mesh, &
       "polar_discrete", &
       x1_interp, &
       x2_interp, &
       j_interp, &
       x1_tab, &
       x2_tab, &
       jacobians_node=jacs )

  call tau%write_to_file()

  dg_rho => sll_o_new( degree, tau, f_four ) 
  dg_ex  => sll_o_new( degree, tau ) 
  dg_ey  => sll_o_new( degree, tau ) 

  call dg_rho%write_to_file('rho')

  call sll_o_create(solver, tau, degree )
  call sll_o_solve(solver, dg_rho, dg_ex, dg_ey)
  call sll_o_delete(solver)

  call dg_ex%write_to_file('ex')
  call dg_ey%write_to_file('ey')

contains

  real(8) function f_cos( r, theta )

    real(8) :: r
    real(8) :: theta

    f_cos = (r-R_MIN)*(r-R_MAX)*cos(N*theta)*r

  end function f_cos

  real(8) function f_four( x, y)
    real(8) :: x, y
#ifdef DEBUG
    real(8) :: dummy
    dummy = x+y
#endif
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
            + 2.0_8*((r-R_MAX)*cos(N*theta) + (r-R_MIN)*cos(N*theta) &
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
          + 2.0_8*((r-R_MAX)*sin(N*theta) + (r-R_MIN)*sin(N*theta)  &
          + r*sin(N*theta))*r)/r

  end function lap_f_sin

end program test_lobalap_discrete
