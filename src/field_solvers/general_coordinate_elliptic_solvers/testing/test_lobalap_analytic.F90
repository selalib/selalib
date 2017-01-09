module lobalap_analytic_helper_functions

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

end module lobalap_analytic_helper_functions

program test_lobalap
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use sll_m_cartesian_meshes, only: &
    sll_f_new_cartesian_mesh_2d, &
    sll_t_cartesian_mesh_2d, &
    sll_o_delete, &
    sll_o_new

  use sll_m_common_coordinate_transformations, only: &
    sll_f_identity_jac11, &
    sll_f_identity_jac12, &
    sll_f_identity_jac21, &
    sll_f_identity_jac22, &
    sll_f_identity_x1, &
    sll_f_identity_x2

  use sll_m_coordinate_transformation_2d_base, only: &
    sll_c_coordinate_transformation_2d_base

  use sll_m_coordinate_transformations_2d, only: &
    sll_f_new_coordinate_transformation_2d_analytic

  use sll_m_dg_fields, only: &
    sll_t_dg_field_2d, &
    sll_o_new

  use sll_m_lobatto_poisson, only: &
    sll_t_lobatto_poisson_solver, &
    sll_o_create, &
    sll_o_solve, &
    sll_o_delete

  use lobalap_analytic_helper_functions

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type(sll_t_lobatto_poisson_solver)                          :: solver
  type(sll_t_cartesian_mesh_2d), pointer                  :: mesh
  class(sll_c_coordinate_transformation_2d_base), pointer :: tau
  type(sll_t_dg_field_2d), pointer                        :: dg_rho
  type(sll_t_dg_field_2d), pointer                        :: dg_ex
  type(sll_t_dg_field_2d), pointer                        :: dg_ey
  sll_int32, parameter                                  :: degree = 3
  
#define NPTS1 2
#define NPTS2 2
#define R_MIN  0.0_8
#define R_MAX  1.0_8
#define THETA_MIN  0.00_8
#define THETA_MAX  1.00_8
#define N 6

  ! logical mesh for space coordinates
  mesh => sll_f_new_cartesian_mesh_2d( NPTS1, NPTS2,    & 
       eta1_min= R_MIN, eta1_max= R_MAX,    &
       eta2_min= THETA_MIN, eta2_max= THETA_MAX)

!  ! coordinate transformation associated with space coordinates
!  tau => sll_f_new_coordinate_transformation_2d_analytic( &
!       "analytic_polar_transformation", &
!       mesh, &
!       sll_f_polar_x1, &
!       sll_f_polar_x2, &
!       sll_f_polar_jac11, &
!       sll_f_polar_jac12, &
!       sll_f_polar_jac21, &
!       sll_f_polar_jac22, &
!       (/ 0.0_f64 /) ) 

  tau => sll_f_new_coordinate_transformation_2d_analytic( &
       "analytic_identity_transformation", &
       mesh, &
       sll_f_identity_x1, &
       sll_f_identity_x2, &
       sll_f_identity_jac11, &
       sll_f_identity_jac12, &
       sll_f_identity_jac21, &
       sll_f_identity_jac22, &
       (/ 0.0_f64 /) )

!  tau => sll_f_new_coordinate_transformation_2d_analytic( &
!    "analytic_colela_transformation", &
!    mesh, &
!    sll_f_sinprod_x1, &
!    sll_f_sinprod_x2, &
!    sll_f_sinprod_jac11, &
!    sll_f_sinprod_jac12, &
!    sll_f_sinprod_jac21, &
!    sll_f_sinprod_jac22, &
!    (/ 0.5_f64,0.5_f64,4.0_f64*sll_p_pi,4.0_f64*sll_p_pi /) )

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

end program test_lobalap

