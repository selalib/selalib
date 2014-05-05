program test_lobalap
#include "selalib.h"
#include "sll_coordinate_transformations.h"

  use map_function_module, only: set_map_function
  use sll_lobatto_poisson
  use sll_dg_fields
  use sll_logical_meshes
  implicit none

  type(lobatto_poisson_solver)        :: solver
  type(sll_logical_mesh_2d), pointer  :: mesh
  class(sll_coordinate_transformation_2d_analytic), pointer :: tau
  type(dg_field), pointer :: dg_rho
  type(dg_field), pointer :: dg_ex
  type(dg_field), pointer :: dg_ey

  sll_int32, parameter :: degree = 2
  real(8), external :: f_cos, f_four
  
#define NPTS1 8
#define NPTS2 8
#define R_MIN  0.0_8
#define R_MAX  1.0_8
#define THETA_MIN  0.00_8
#define THETA_MAX  1.00_8
#define N 6

  ! logical mesh for space coordinates
  mesh => new_logical_mesh_2d( NPTS1, NPTS2,    & 
       eta1_min= R_MIN, eta1_max= R_MAX,    &
       eta2_min= THETA_MIN, eta2_max= THETA_MAX)

!  ! coordinate transformation associated with space coordinates
!  tau => new_coordinate_transformation_2d_analytic( &
!       "analytic_polar_transformation", &
!       mesh, &
!       polar_x1, &
!       polar_x2, &
!       polar_jac11, &
!       polar_jac12, &
!       polar_jac21, &
!       polar_jac22, &
!       (/ 0.0_f64 /) ) 

  tau => new_coordinate_transformation_2d_analytic( &
       "analytic_identity_transformation", &
       mesh, &
       identity_x1, &
       identity_x2, &
       identity_jac11, &
       identity_jac12, &
       identity_jac21, &
       identity_jac22, &
       (/ 0.0_f64 /) )


!  tau => new_coordinate_transformation_2d_analytic( &
!    "analytic_colela_transformation", &
!    mesh, &
!    sinprod_x1, &
!    sinprod_x2, &
!    sinprod_jac11, &
!    sinprod_jac12, &
!    sinprod_jac21, &
!    sinprod_jac22, &
!    (/ 0.5_f64,0.5_f64,4.0_f64*sll_pi,4.0_f64*sll_pi /) )

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

end program test_lobalap

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
