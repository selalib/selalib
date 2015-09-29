program test_gces_dirichlet
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_utilities.h"
#include "sll_file_io.h"

use sll_cartesian_meshes
use sll_m_coordinate_transformations_2d
use sll_common_coordinate_transformations
use sll_m_scalar_field_2d
use sll_constants
use sll_m_arbitrary_degree_spline_interpolator_2d
use sll_m_deboor_splines_2d
use sll_m_gces_dirichlet

implicit none

#define SPLINE_DEG1       3
#define SPLINE_DEG2       3
#define NUM_CELLS1        10
#define NUM_CELLS2        10
#define ETA1MIN          (-1.0_f64)
#define ETA1MAX          (+1.0_f64)
#define ETA2MIN          (-1.0_f64)
#define ETA2MAX          (+1.0_f64)

type(sll_cartesian_mesh_2d),                       pointer :: mesh_2d
class(sll_coordinate_transformation_2d_base),      pointer :: tau
type(sll_gces_dirichlet)                                   :: es
type(sll_arbitrary_degree_spline_interpolator_2d), target  :: interp_phi
class(sll_scalar_field_2d_base),                   pointer :: a11_field_mat
class(sll_scalar_field_2d_base),                   pointer :: a12_field_mat
class(sll_scalar_field_2d_base),                   pointer :: a21_field_mat
class(sll_scalar_field_2d_base),                   pointer :: a22_field_mat
class(sll_scalar_field_2d_base),                   pointer :: b1_field_vect
class(sll_scalar_field_2d_base),                   pointer :: b2_field_vect
class(sll_scalar_field_2d_base),                   pointer :: c_field
class(sll_scalar_field_2d_base),                   pointer :: rho
type(sll_scalar_field_2d_discrete),                pointer :: phi

sll_real64 :: acc
sll_real64 :: normL2
sll_real64 :: normH1

sll_real64, dimension(NUM_CELLS1+1,NUM_CELLS2+1) :: values
sll_real64, dimension(NUM_CELLS1+1,NUM_CELLS2+1) :: calculated
sll_real64, dimension(NUM_CELLS1+1,NUM_CELLS2+1) :: reference

sll_int32  :: i, j, k
sll_real64 :: h1,h2,node_val,ref
sll_real64 :: eta1(NUM_CELLS1+1)
sll_real64 :: eta2(NUM_CELLS2+1)
sll_int32  :: npts1,npts2
sll_int32  :: cell 

real(8) :: integral_solution
real(8) :: integral_exact_solution

sll_real64 :: grad1_node_val,grad2_node_val,grad1ref,grad2ref

mesh_2d => new_cartesian_mesh_2d( NUM_CELLS1, &
                                  NUM_CELLS2, &
                                  ETA1MIN,    &
                                  ETA1MAX,    &
                                  ETA2MIN,    &
                                  ETA2MAX )
acc   = 0.0_f64
npts1 =  NUM_CELLS1 + 1
npts2 =  NUM_CELLS2 + 1
h1    = (ETA1MAX-ETA1MIN)/real(NPTS1-1,f64)
h2    = (ETA2MAX-ETA2MIN)/real(NPTS2-1,f64)

do j=1,npts2
  do i=1,npts1
    eta1(i)  = (i-1)*h1 + ETA1MIN
    eta2(j)  = (j-1)*h2 + ETA2MIN
  end do
end do

normL2 = 0.0_f64
normH1 = 0.0_f64

print*, "-------------------------------------------------------------"
print*, " test case witout change of coordinates"
print*, " dirichlet-dirichlet boundary conditions"
print*, "-------------------------------------------------------------"

tau => new_coordinate_transformation_2d_analytic( &
&      "analytic",                                &
&      mesh_2d,                                   &
&      identity_x1,                               &
&      identity_x2,                               &
&      identity_jac11,                            &
&      identity_jac12,                            &
&      identity_jac21,                            &
&      identity_jac22,                            &
&      [0.0_f64]                                  )

a11_field_mat => new_scalar_field_2d_analytic( &
  one,                                         &
  "a11",                                       &
  tau,                                         &
  SLL_DIRICHLET,                               &
  SLL_DIRICHLET,                               &
  SLL_DIRICHLET,                               &
  SLL_DIRICHLET,                               &
  [0.0_f64]  ) 

a12_field_mat => new_scalar_field_2d_analytic( &
  zero,                                        &
  "a12",                                       &
  tau,                                         &
  SLL_DIRICHLET,                               &
  SLL_DIRICHLET,                               &
  SLL_DIRICHLET,                               &
  SLL_DIRICHLET,                               &
  [0.0_f64] )

a21_field_mat => new_scalar_field_2d_analytic( &
  zero,                                        &
  "a21",                                       &
  tau,                                         &
  SLL_DIRICHLET,                               &
  SLL_DIRICHLET,                               &
  SLL_DIRICHLET,                               &
  SLL_DIRICHLET,                               &
  [0.0_f64] ) 

a22_field_mat => new_scalar_field_2d_analytic( &
  one,                                         &
  "a22",                                       &
  tau,                                         &
  SLL_DIRICHLET,                               &
  SLL_DIRICHLET,                               &
  SLL_DIRICHLET,                               &
  SLL_DIRICHLET,                               &
  [0.0_f64])

b1_field_vect => new_scalar_field_2d_analytic( &
  zero,                                        &
  "b1",                                        &
  tau,                                         &
  SLL_DIRICHLET,                               &
  SLL_DIRICHLET,                               &
  SLL_DIRICHLET,                               &
  SLL_DIRICHLET,                               &
  [0.0_f64],                                   & 
  first_deriv_eta1 = zero,                     &
  first_deriv_eta2 = zero) 

b2_field_vect => new_scalar_field_2d_analytic( &
  zero,                                        &
  "b2",                                        &
  tau,                                         &
  SLL_DIRICHLET,                               &
  SLL_DIRICHLET,                               &
  SLL_DIRICHLET,                               &
  SLL_DIRICHLET,                               &
  [0.0_f64],                                   &
  first_deriv_eta1 = zero,                     &
  first_deriv_eta2 = zero)

c_field => new_scalar_field_2d_analytic(       &
  zero,                                        &
  "c_field",                                   &
  tau,                                         &
  SLL_DIRICHLET,                               &
  SLL_DIRICHLET,                               &
  SLL_DIRICHLET,                               &
  SLL_DIRICHLET,                               &
  [0.0_f64]  )

call initialize_ad2d_interpolator(             &
  interp_phi,                                  &
  NUM_CELLS1+1,                                &
  NUM_CELLS2+1,                                &
  ETA1MIN,                                     &
  ETA1MAX,                                     &
  ETA2MIN,                                     &
  ETA2MAX,                                     &
  SLL_HERMITE,                               &
  SLL_HERMITE,                               &
  SLL_HERMITE,                               &
  SLL_HERMITE,                               &
  SPLINE_DEG1,                                 &
  SPLINE_DEG2 )

phi => new_scalar_field_2d_discrete(           &
  "phi",                                       &
  interp_phi,                                  &
  tau,                                         &
  SLL_HERMITE,                               &
  SLL_HERMITE,                               &
  SLL_HERMITE,                               &
  SLL_HERMITE )

rho => new_scalar_field_2d_analytic( &
&    rhs,                            &
&    "rho",                          &     
&    tau,                            &
&    SLL_DIRICHLET,                  &
&    SLL_DIRICHLET,                  &
&    SLL_DIRICHLET,                  &
&    SLL_DIRICHLET,                  &
&    [0.0_f64]                       )

integral_solution = 0.0_f64
integral_exact_solution = 0.0_f64

call sll_create( es,                  &
&                SPLINE_DEG1,         &
&                SPLINE_DEG2,         &
&                NUM_CELLS1,          &
&                NUM_CELLS2,          &
&                ES_GAUSS_LEGENDRE,   &
&                ES_GAUSS_LEGENDRE,   &
&                ETA1MIN,             &
&                ETA1MAX,             &
&                ETA2MIN,             &
&                ETA2MAX              )
 
call factorize_mat_es( es,            &
&                      a11_field_mat, &
&                      a12_field_mat, &
&                      a21_field_mat, &
&                      a22_field_mat, &
&                      b1_field_vect, &
&                      b2_field_vect, &
&                      c_field        )

do k = 1, 10
  call sll_solve( es, rho, phi)
end do

integral_solution       = 0.0_f64
integral_exact_solution = 0.0_f64

call phi%write_to_file(0)

normL2 = 0.0_f64
normH1 = 0.0_f64
do j=1,npts2
do i=1,npts1
  node_val        = phi%value_at_point(eta1(i),eta2(j))
  calculated(i,j) = node_val
  grad1_node_val  = phi%first_deriv_eta1_value_at_point(eta1(i), eta2(j))
  grad2_node_val  = phi%first_deriv_eta2_value_at_point(eta1(i), eta2(j))
  ref             = sol(eta1(i),eta2(j), [0.0d0])
  grad1ref        = 2*sll_pi*cos(2*sll_pi*eta1(i))*sin(2*sll_pi*eta2(j))
  grad2ref        = 2*sll_pi*cos(2*sll_pi*eta2(j))*sin(2*sll_pi*eta1(i))
  reference(i,j)  = ref
  normL2          = normL2 + (node_val-ref)**2*h1*h2
  normH1          = normH1 + ((grad1_node_val-grad1ref)**2+(grad2_node_val-grad2ref)**2)*h1*h2
end do
end do

integral_solution = sum(calculated)*h1*h2
integral_exact_solution = sum(reference)*h1*h2

do j=-100,100
do i=-100,100
  write(41,*) i, j, phi%first_deriv_eta1_value_at_point(i*0.01_f64,j*0.01_f64) &
                  , 2*sll_pi*cos(2*sll_pi*i*0.01_f64)*sin(2*sll_pi*j*0.01_f64)
  write(42,*) i, j, phi%value_at_point(i*0.01_f64,j*0.01_f64) &
                  , sol(i*0.01_f64,j*0.01_f64, [0.0d0])
end do
write(41,*) 
write(42,*) 
end do

call sll_delete(es)
call rho%delete()
call c_field%delete()
call phi%delete()
call a11_field_mat%delete()
call a12_field_mat%delete()
call a21_field_mat%delete()
call b1_field_vect%delete()
call b2_field_vect%delete()
call a22_field_mat%delete()
call tau%delete()

print"('integral solution       =',g15.3)", integral_solution
print"('integral exact solution =',g15.3)", integral_exact_solution
print*, ' L2 norm :', sqrt(normL2), h1**(SPLINE_DEG1-1)
print*, ' H1 norm :', sqrt(normH1), h1**(SPLINE_DEG1-2)
if (sqrt(normL2) <= h1**(SPLINE_DEG1-1)  .and. &
    sqrt(normH1) <= h1**(SPLINE_DEG1-2)) then     
   acc = sum(abs(calculated-reference))/(npts1*npts2)
   print"('L_oo =',g15.3, 4x, 'OK' )", acc
else
  stop 'FAILED'
end if
print*, 'PASSED'

contains

function four( eta1, eta2, params ) result(res)
real(8), intent(in) :: eta1
real(8), intent(in) :: eta2
real(8), dimension(:), intent(in) :: params
real(8) :: res, pi
res = 4.0_f64
end function four

function one( eta1, eta2, params ) result(res)
real(8), intent(in) :: eta1
real(8), intent(in) :: eta2
real(8), dimension(:), intent(in) :: params
real(8) :: res
res = 1.0_f64
end function one

function zero( eta1, eta2, params ) result(res)
real(8), intent(in) :: eta1
real(8), intent(in) :: eta2
real(8), dimension(:), intent(in) :: params
real(8) :: res
res = 0.0_f64
end function zero

function rhs( eta1, eta2, params ) result(res)
real(8), intent(in) :: eta1
real(8), intent(in) :: eta2
real(8), dimension(:), intent(in) :: params
real(8) :: res
real(8) :: pi

pi = 4d0*atan(1d0)
res = -8*pi*pi*sin(2*pi*eta1)*sin(2*pi*eta2)

end function rhs

function sol( eta1, eta2, params ) result(res)
real(8), intent(in) :: eta1
real(8), intent(in) :: eta2
real(8), dimension(:), intent(in) :: params
real(8) :: res
real(8) :: pi

pi = 4d0*atan(1d0)
res = sin(2*pi*eta1)*sin(2*pi*eta2)

end function sol

end program test_gces_dirichlet
