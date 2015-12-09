!Test for new version of general coordinates elliptic solver
!made by Adnane and Michel
program test_gces_full_periodic_prototype
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_errors.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_arbitrary_degree_spline_interpolator_2d, only: &
    new_arbitrary_degree_spline_interp2d

  use sll_m_boundary_condition_descriptors, only: &
    sll_periodic

  use sll_m_cartesian_meshes, only: &
    new_cartesian_mesh_2d, &
    sll_cartesian_mesh_2d

  use sll_m_common_coordinate_transformations, only: &
    sinprod_jac11, &
    sinprod_jac12, &
    sinprod_jac21, &
    sinprod_jac22, &
    sinprod_x1, &
    sinprod_x2

  use sll_m_constants, only: &
    sll_pi

  use sll_m_coordinate_transformation_2d_base, only: &
    sll_coordinate_transformation_2d_base

  use sll_m_coordinate_transformations_2d, only: &
    new_coordinate_transformation_2d_analytic

  use sll_m_cubic_spline_interpolator_2d, only: &
    new_cubic_spline_interpolator_2d

  use sll_m_general_coordinate_elliptic_solver, only: &
    es_gauss_legendre, &
    factorize_mat_es_prototype, &
    general_coordinate_elliptic_solver, &
    initialize_general_elliptic_solver_prototype, &
    solve_general_coordinates_elliptic_eq_prototype

  use sll_m_interpolators_2d_base, only: &
    sll_c_interpolator_2d

  use sll_m_scalar_field_2d, only: &
    new_scalar_field_2d_analytic, &
    new_scalar_field_2d_discrete

  use sll_m_scalar_field_2d_base, only: &
    sll_scalar_field_2d_base

  implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

type(sll_cartesian_mesh_2d),                  pointer :: mesh_2d
class(sll_coordinate_transformation_2d_base), pointer :: tau
type(general_coordinate_elliptic_solver)              :: es
class(sll_c_interpolator_2d),              pointer :: interp_rho
class(sll_scalar_field_2d_base),              pointer :: a11_field_mat
class(sll_scalar_field_2d_base),              pointer :: a12_field_mat
class(sll_scalar_field_2d_base),              pointer :: a21_field_mat
class(sll_scalar_field_2d_base),              pointer :: a22_field_mat
class(sll_scalar_field_2d_base),              pointer :: b1_field_vect
class(sll_scalar_field_2d_base),              pointer :: b2_field_vect
class(sll_scalar_field_2d_base),              pointer :: c_field
class(sll_scalar_field_2d_base),              pointer :: rho

sll_real64 :: acc
sll_real64 :: normL2

sll_real64, dimension(:,:), allocatable :: values
sll_real64, dimension(:,:), allocatable :: calculated
sll_real64, dimension(:,:), allocatable :: reference

sll_int32  :: i, j, k
sll_int32  :: npts1,npts2
sll_real64 :: h1,h2,node_val,ref
sll_real64, dimension(:), allocatable  :: eta1
sll_real64, dimension(:), allocatable  :: eta2

sll_real64 :: integral_solution
sll_real64 :: integral_exact_solution
sll_real64 :: x1 
sll_real64 :: x2 

sll_real64, dimension(:,:), allocatable :: phi
sll_real64, dimension(:,:), allocatable :: tab_rho
sll_int32 :: ierr

sll_int32 :: num_cells1
sll_int32 :: num_cells2
sll_real64 :: alpha1
sll_real64 :: alpha2
sll_real64 :: eta1_min
sll_real64 :: eta1_max
sll_real64 :: eta2_min
sll_real64 :: eta2_max
sll_int32 :: spline_degree1
sll_int32 :: spline_degree2
character(256) :: filename
character(256) :: filename_loc
sll_int32             :: IO_stat
sll_int32, parameter  :: input_file = 99
sll_real64 :: params_mesh(4)
sll_real64 :: params_rhs(6)
sll_int32 :: mode1
sll_int32 :: mode2
character(len=256) :: rho_case
character(len=256) :: rho_interp2d_case
  namelist /params/ &
    num_cells1, &
    num_cells2, &
    spline_degree1, &
    spline_degree2, &
    eta1_min, &
    eta1_max, &
    eta2_min, &
    eta2_max, &
    mode1, &
    mode2, &
    alpha1, &
    alpha2, &
    rho_case, &
    rho_interp2d_case

!default values
num_cells1        = 128
num_cells2        = 128
spline_degree1    = 3
spline_degree2    = 3
eta1_min          = 0._f64
eta1_max          = 1._f64
eta2_min          = 0._f64
eta2_max          = 1._f64
alpha1            = 0._f64
alpha2            = 0._f64
mode1             = 1
mode2             =  1
rho_case          = "SLL_ANALYTIC" 
rho_interp2d_case = "SLL_CUBIC_SPLINES"

call get_command_argument(1, filename)

if (len_trim(filename) .ne. 0)then
  filename_loc = filename
  filename_loc = adjustl(filename_loc)
  open(unit = input_file, file=trim(filename_loc)//'.nml',IOStat=IO_stat)
    if( IO_stat /= 0 ) then
      print *, '#failed to open file ', &
      trim(filename)//'.nml'
      STOP
    end if
  read(input_file, params) 
  close(input_file)
else
  print *,'#initialization with default parameters'    
endif

print *,"#min of jacobian is ", &
  1._f64-alpha1*sll_pi/(eta1_max-eta1_min) &
  -alpha2*sll_pi/(eta2_max-eta2_min)

params_mesh(1:4) = (/ alpha1, alpha2, eta1_max-eta1_min, eta2_max-eta2_min/)
params_rhs(1:4) = params_mesh(1:4)
params_rhs(5:6) = (/real(mode1,f64),real(mode2,f64)/)
SLL_ALLOCATE(phi(num_cells1+1,num_cells2+1),ierr)
SLL_ALLOCATE(values(num_cells1+1,num_cells2+1),ierr)
SLL_ALLOCATE(calculated(num_cells1+1,num_cells2+1),ierr)
SLL_ALLOCATE(reference(num_cells1+1,num_cells2+1),ierr)
SLL_ALLOCATE(tab_rho(num_cells1+1,num_cells2+1),ierr)
SLL_ALLOCATE(eta1(num_cells1+1),ierr)
SLL_ALLOCATE(eta2(num_cells2+1),ierr)

mesh_2d => new_cartesian_mesh_2d( num_cells1, &
                                  num_cells2, &
                                  eta1_min,   &
                                  eta1_max,   &
                                  eta2_min,   &
                                  eta2_max )
acc   = 0.0_f64
npts1 =  num_cells1 + 1
npts2 =  num_cells2 + 1
h1    = (eta1_max-eta1_min)/real(npts1-1,f64)
h2    = (eta1_max-eta1_min)/real(npts2-1,f64)

do j=1,npts2
  do i=1,npts1
    eta1(i)  = (i-1)*h1 + eta1_min
    eta2(j)  = (j-1)*h2 + eta2_min
  end do
end do

normL2 = 0.0_f64

print*, "-----------------------------------------------"
print*, " test case with change of coordinates sinprod  "
print*, " dirichlet-dirichlet boundary conditions       "
print*, "-----------------------------------------------"

tau => new_coordinate_transformation_2d_analytic( &
     "analytic",                                  &
     mesh_2d,                                     &
     sinprod_x1,                                  &
     sinprod_x2,                                  &
     sinprod_jac11,                               &
     sinprod_jac12,                               &
     sinprod_jac21,                               &
     sinprod_jac22,                               &
     params_mesh )

a11_field_mat => new_scalar_field_2d_analytic( &
  one,                                         &
  "a11",                                       &
  tau,                                         &
  SLL_PERIODIC,                                &
  SLL_PERIODIC,                                &
  SLL_PERIODIC,                                &
  SLL_PERIODIC,                                &
  [0.0_f64]  ) 

a12_field_mat => new_scalar_field_2d_analytic( &
  zero,                                        &
  "a12",                                       &
  tau,                                         &
  SLL_PERIODIC,                                &
  SLL_PERIODIC,                                &
  SLL_PERIODIC,                                &
  SLL_PERIODIC,                                &
  [0.0_f64] )

a21_field_mat => new_scalar_field_2d_analytic( &
  zero,                                        &
  "a21",                                       &
  tau,                                         &
  SLL_PERIODIC,                                &
  SLL_PERIODIC,                                &
  SLL_PERIODIC,                                &
  SLL_PERIODIC,                                &
  [0.0_f64] ) 

a22_field_mat => new_scalar_field_2d_analytic( &
  one,                                         &
  "a22",                                       &
  tau,                                         &
  SLL_PERIODIC,                                &
  SLL_PERIODIC,                                &
  SLL_PERIODIC,                                &
  SLL_PERIODIC,                                &
  [0.0_f64])

b1_field_vect => new_scalar_field_2d_analytic( &
  zero,                                        &
  "b1",                                        &
  tau,                                         &
  SLL_PERIODIC,                                &
  SLL_PERIODIC,                                &
  SLL_PERIODIC,                                &
  SLL_PERIODIC,                                &
  [0.0_f64],                                   & 
  first_deriv_eta1 = zero,                     &
  first_deriv_eta2 = zero) 

b2_field_vect => new_scalar_field_2d_analytic( &
  zero,                                        &
  "b2",                                        &
  tau,                                         &
  SLL_PERIODIC,                                &
  SLL_PERIODIC,                                &
  SLL_PERIODIC,                                &
  SLL_PERIODIC,                                &
  [0.0_f64],                                   &
  first_deriv_eta1 = zero,                     &
  first_deriv_eta2 = zero)

c_field => new_scalar_field_2d_analytic(       &
  zero,                                        &
  "c_field",                                   &
  tau,                                         &
  SLL_PERIODIC,                                &
  SLL_PERIODIC,                                &
  SLL_PERIODIC,                                &
  SLL_PERIODIC,                                &
  [0.0_f64]  )

select case (rho_interp2d_case)

  case ("SLL_CUBIC_SPLINES")
    interp_rho => new_cubic_spline_interpolator_2d( &
      npts1,                                        &                                
      npts2,                                        &                                
      eta1_min,                                     &                                    
      eta1_max,                                     &                                    
      eta2_min,                                     &                                    
      eta2_max,                                     &                                     
      SLL_PERIODIC,                                 &                                
      SLL_PERIODIC)                                

  case ("SLL_ARBITRARY_DEGREE_SPLINES")
    interp_rho => new_arbitrary_degree_spline_interp2d( &
      npts1,                                            &
      npts2,                                            &
      eta1_min,                                         &
      eta1_max,                                         &
      eta2_min,                                         &
      eta2_max,                                         &
      SLL_PERIODIC,                                     &
      SLL_PERIODIC,                                     &
      SLL_PERIODIC,                                     &
      SLL_PERIODIC,                                     &
      spline_degree1,                                   &
      spline_degree2 )

  case default

    SLL_ERROR("program", "bad value of rho_interp2d_case")

end select

if (rho_case=="SLL_ANALYTIC") then

  rho => new_scalar_field_2d_analytic( &
  &    rhs,                            &
  &    "rho",                          &     
  &    tau,                            &
  &    SLL_PERIODIC,                   &
  &    SLL_PERIODIC,                   &
  &    SLL_PERIODIC,                   &
  &    SLL_PERIODIC,                   &
  &    params_rhs )

else if (rho_case=="SLL_DISCRETE") then

  rho => new_scalar_field_2d_discrete( &
       "rhototo",                      &
       interp_rho,                     &
       tau,                            &
       SLL_PERIODIC,                   &
       SLL_PERIODIC,                   &
       SLL_PERIODIC,                   &
       SLL_PERIODIC,                   &
       eta1,                           &
       num_cells1,                     &
       eta2,                           &
       num_cells2)  

  do j=1,npts2
    do i=1,npts1
      tab_rho(i,j) = rhs(eta1(i),eta2(j), params_rhs)
    end do
  end do

  call rho%set_field_data(tab_rho)
  call rho%update_interpolation_coefficients()

else
  SLL_ERROR("program","bad value of rho_case")
endif

integral_solution       = 0.0_f64
integral_exact_solution = 0.0_f64

call initialize_general_elliptic_solver_prototype( es, &
&                spline_degree1,                       &
&                spline_degree2,                       &
&                npts1-1,                              &
&                npts2-1,                              &
&                ES_GAUSS_LEGENDRE,                    &
&                ES_GAUSS_LEGENDRE,                    &
&                SLL_PERIODIC,                         &
&                SLL_PERIODIC,                         &
&                SLL_PERIODIC,                         &
&                SLL_PERIODIC,                         &
&                eta1_min,                             &
&                eta1_max,                             &
&                eta2_min,                             &
&                eta2_max                              )
 
 
call factorize_mat_es_prototype( es,  &
&                      a11_field_mat, &
&                      a12_field_mat, &
&                      a21_field_mat, &
&                      a22_field_mat, &
&                      b1_field_vect, &
&                      b2_field_vect, &
&                      c_field        )

do k = 1, 1
  call solve_general_coordinates_elliptic_eq_prototype(&
    es, &
    phi, &
    rho_field = rho)
end do

integral_solution       = 0.0_f64
integral_exact_solution = 0.0_f64

normL2 = 0.0_f64
do j=1,npts2
do i=1,npts1
  x1              = tau%x1(eta1(i),eta2(j))
  x2              = tau%x2(eta1(i),eta2(j))
  node_val        = phi(i,j)
  calculated(i,j) = node_val
  ref             = sol( eta1(i), eta2(j), params_rhs )
  reference(i,j)  = ref
  normL2          = normL2 + (node_val-ref)**2*h1*h2
end do
end do

integral_solution       = sum(calculated)*h1*h2
integral_exact_solution = sum(reference)*h1*h2

call rho%delete()
call c_field%delete()
call a11_field_mat%delete()
call a12_field_mat%delete()
call a21_field_mat%delete()
call b1_field_vect%delete()
call b2_field_vect%delete()
call a22_field_mat%delete()
call tau%delete()

print"('integral solution       =',g15.3)", integral_solution
print"('integral exact solution =',g15.3)", integral_exact_solution
print*, ' L2 norm :', sqrt(normL2), &
  sqrt(normL2)/h1**(spline_degree1+1), &
  sqrt(normL2)/h1**(spline_degree1+2)
if (sqrt(normL2) <= h1**(spline_degree1-1) ) then     
   acc = sum(abs(calculated-reference))/real(npts1*npts2,f64)
   print"('L_oo =',g15.3, 4x, 'OK' )", acc
else
  stop 'FAILED'
end if
print*, 'PASSED'

contains

function one( eta1, eta2, params ) result(res)
real(8), intent(in) :: eta1
real(8), intent(in) :: eta2
real(8), dimension(:), intent(in) :: params
real(8) :: res
#ifdef DEBUG
real(8) :: dummy
dummy = eta1+eta2+params(1)
#endif
res = 1.0_f64
end function one

function zero( eta1, eta2, params ) result(res)
real(8), intent(in) :: eta1
real(8), intent(in) :: eta2
real(8), dimension(:), intent(in) :: params
real(8) :: res
#ifdef DEBUG
real(8) :: dummy
dummy = eta1+eta2+params(1)
#endif
res = 0.0_f64
end function zero

function rhs( eta1, eta2, params ) result(res)
real(8), intent(in)               :: eta1
real(8), intent(in)               :: eta2
real(8), dimension(:), intent(in) :: params
real(8)                           :: res
real(8)                           :: pi
real(8)                           :: x1
real(8)                           :: x2
real(8)                           :: L1
real(8)                           :: L2
real(8)                           :: mode1
real(8)                           :: mode2

L1 = params(3)
L2 = params(4)
mode1 = params(5)
mode2 = params(6)

x1  = sinprod_x1(eta1,eta2, params)
x2  = sinprod_x2(eta1,eta2, params)
pi  = 4._f64*atan(1._f64)
res = (mode1*2.0_f64*pi/L1)**2
res = res+(mode2*2.0_f64*pi/L2)**2
res = -res*cos(2*pi*mode1*x1/L1)*cos(2*pi*mode2*x2/L2)

end function rhs

function sol( eta1, eta2, params ) result(res)
real(8), intent(in) :: eta1
real(8), intent(in) :: eta2
real(8), dimension(:), intent(in) :: params
real(8) :: res
real(8) :: pi
real(8) :: x1
real(8) :: x2
real(8) :: L1
real(8) :: L2
real(8) :: mode1
real(8) :: mode2

L1 = params(3)
L2 = params(4)
mode1 = params(5)
mode2 = params(6)
x1 = sinprod_x1(eta1,eta2, params)
x2 = sinprod_x2(eta1,eta2, params)
pi = 4._f64*atan(1._f64)
res = cos(2.0_f64*pi*mode1*x1/L1)*cos(2.0_f64*pi*mode2*x2/L2)
end function sol

end program test_gces_full_periodic_prototype

!  L2 norm :   7.7856894257913913E-006   1.3409651175147714        27.317917051683867     
!  L2 norm :   2.2477687837729589E-005   3.8714356127450573        78.868238672688207     
!  L2 norm :   4.1032908760074900E-006  0.70672866984844762        14.397358237585992     
