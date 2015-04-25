program test_general_elliptic_solver_multipatch

#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_file_io.h"

use sll_cartesian_meshes
use sll_cartesian_meshes_multipatch
use sll_module_coordinate_transformations_2d
use sll_common_coordinate_transformations
use sll_coordinate_transformation_multipatch_module
use sll_module_scalar_field_2d
use sll_module_scalar_field_2d_multipatch
use sll_constants
use sll_module_arbitrary_degree_spline_interpolator_2d
use sll_timer
use sll_general_coordinate_elliptic_solver_multipatch_module

implicit none

#define SPLINE_DEG1       3
#define SPLINE_DEG2       3
#define NUM_CELLS1        64
#define NUM_CELLS2        64
#define ETA1MIN           0.0_f64
#define ETA1MAX           1.0_f64
#define ETA2MIN           0.0_f64
#define ETA2MAX           1.0_f64
#define PRINT_COMPARISON  .false.

type(sll_coordinate_transformation_multipatch_2d), pointer :: T
class(sll_cartesian_mesh_2d), pointer                        :: m
class(sll_coordinate_transformation_2d_nurbs), pointer     :: transf
type(general_coordinate_elliptic_solver_mp)                :: es_mp
class(sll_scalar_field_multipatch_2d), pointer             :: a11_field_mat
class(sll_scalar_field_multipatch_2d), pointer             :: a12_field_mat
class(sll_scalar_field_multipatch_2d), pointer             :: a21_field_mat
class(sll_scalar_field_multipatch_2d), pointer             :: a22_field_mat
class(sll_scalar_field_multipatch_2d), pointer             :: b1_field_vect
class(sll_scalar_field_multipatch_2d), pointer             :: b2_field_vect
class(sll_scalar_field_multipatch_2d), pointer             :: c_field_scal
class(sll_scalar_field_multipatch_2d), pointer             :: rho_field_scal
class(sll_scalar_field_multipatch_2d), pointer             :: phi_field_scal
class(sll_scalar_field_multipatch_2d), pointer             :: phi_field_ex
class(sll_scalar_field_multipatch_2d), pointer             :: phi_field_diff

sll_int32  :: num_patches
sll_int32  :: ipatch
sll_int32  :: i
sll_int32  :: j
sll_int32  :: num_pts1
sll_int32  :: num_pts2
sll_real64 :: val_a11,val_a12,val_a21,val_a22,val_b1,val_b2,val_c,val_rho,val_phi,val_phi_exacte
sll_real64 :: eta1,eta1min
sll_real64 :: eta2,eta2min
sll_real64 :: delta1
sll_real64 :: delta2
sll_real64 :: x1
sll_real64 :: x2

real(8), external    :: func_zero
real(8), external    :: func_one
real(8), external    :: func_epsi
real(8), external    :: source_term_perdir
real(8), external    :: source_term_dirper
sll_real64, external :: sol_exacte_perdir

T => new_coordinate_transformation_multipatch_2d("circle_mp5_pts12")

a11_field_mat => new_scalar_field_multipatch_2d("a11_field_multipatch", T)
call a11_field_mat%allocate_memory()
a12_field_mat => new_scalar_field_multipatch_2d("a12_field_multipatch", T)
call a12_field_mat%allocate_memory()
a21_field_mat => new_scalar_field_multipatch_2d("a21_field_multipatch", T)
call a21_field_mat%allocate_memory()
a22_field_mat => new_scalar_field_multipatch_2d("a22_field_multipatch", T)
call a22_field_mat%allocate_memory()
b1_field_vect => new_scalar_field_multipatch_2d("b1_field_multipatch", T)
call b1_field_vect%allocate_memory()
b2_field_vect => new_scalar_field_multipatch_2d("b2_field_multipatch", T)
call b2_field_vect%allocate_memory()
c_field_scal => new_scalar_field_multipatch_2d("c_field_multipatch", T)
call c_field_scal%allocate_memory()
rho_field_scal => new_scalar_field_multipatch_2d("rho_field_multipatch", T)
call rho_field_scal%allocate_memory()
phi_field_scal => new_scalar_field_multipatch_2d("phi_field_multipatch", T)
call phi_field_scal%allocate_memory()
phi_field_ex => new_scalar_field_multipatch_2d("phi_field_ex_multipatch", T)
call phi_field_ex%allocate_memory()
phi_field_diff => new_scalar_field_multipatch_2d("phi_field_diff_multipatch", T)
call phi_field_diff%allocate_memory()

num_patches = phi_field_scal%get_number_patches()

do ipatch= 0,num_patches-1
  m        => rho_field_scal%get_cartesian_mesh(ipatch)
  transf   => rho_field_scal%get_transformation(ipatch)
  num_pts1 = m%num_cells1+1
  num_pts2 = m%num_cells2+1
  delta1   = m%delta_eta1
  delta2   = m%delta_eta2
  eta1min  = m%eta1_min
  eta2min  = m%eta2_min

  do j=1,num_pts1
    eta2 = eta2min + (j-1)*delta2
    do i=1,num_pts2
       eta1 = eta1min + (i-1)*delta1
       x1 = transf%x1(eta1,eta2)
       x2 = transf%x2(eta1,eta2)
       val_a11  = func_one(  x1, x2)
       val_a12  = func_zero( x1, x2)
       val_a21  = func_zero( x1, x2)
       val_a22  = func_one(  x1, x2)
       val_b1   = func_zero( x1, x2)
       val_b2   = func_zero( x1, x2)
       val_c    = func_zero( x1, x2)
       val_rho  = -4*9*exp(-9*(x1**2+x2**2))+(2*9)**2*(x1**2+x2**2)*exp(-9*(x1**2+x2**2))
       val_phi  = 0.0_f64
       val_phi_exacte = exp(-9*(x1**2+x2**2))
       call a11_field_mat%set_value_at_indices ( i, j, ipatch, val_a11 ) 
       call a12_field_mat%set_value_at_indices ( i, j, ipatch, val_a12 ) 
       call a21_field_mat%set_value_at_indices ( i, j, ipatch, val_a21 ) 
       call a22_field_mat%set_value_at_indices ( i, j, ipatch, val_a22 ) 
       call b1_field_vect%set_value_at_indices ( i, j, ipatch, val_b1 ) 
       call b2_field_vect%set_value_at_indices ( i, j, ipatch, val_b2 ) 
       call c_field_scal%set_value_at_indices  ( i, j, ipatch, val_c ) 
       call rho_field_scal%set_value_at_indices( i, j, ipatch, val_rho ) 
       call phi_field_scal%set_value_at_indices( i, j, ipatch, val_phi ) 
       call phi_field_ex%set_value_at_indices( i, j, ipatch, val_phi_exacte ) 
    end do
  end do
end do

print *, 'updating multipatch field interpolation coefficients...'
call a11_field_mat%update_interpolation_coefficients()
print *, 'updating multipatch field interpolation coefficients...'
call a12_field_mat%update_interpolation_coefficients()
print *, 'updating multipatch field interpolation coefficients...'
call a21_field_mat%update_interpolation_coefficients()
print *, 'updating multipatch field interpolation coefficients...'
call a22_field_mat%update_interpolation_coefficients()
print *, 'updating multipatch field interpolation coefficients...'
call b1_field_vect%update_interpolation_coefficients()
print *, 'updating multipatch field interpolation coefficients...'
call b2_field_vect%update_interpolation_coefficients()
print *, 'updating multipatch field interpolation coefficients...'
call c_field_scal%update_interpolation_coefficients()
print *, 'updating multipatch field interpolation coefficients rho...'
call rho_field_scal%update_interpolation_coefficients()
print *, 'updating multipatch field interpolation coefficients...'
call phi_field_scal%update_interpolation_coefficients()
print *, 'updating multipatch field interpolation coefficients...'
call phi_field_ex%update_interpolation_coefficients()

print*, 'Initialization solver elliptic multipacth '

call initialize_general_elliptic_solver_mp( es_mp,             &
                                            ES_GAUSS_LEGENDRE, &
                                            ES_GAUSS_LEGENDRE, &
                                            T)

print*, ' factorise matrix to solve the elleptic solver'
call factorize_mat_es_mp( es_mp,         &
                          a11_field_mat, &
                          a12_field_mat, &
                          a21_field_mat, &
                          a22_field_mat, &
                          b1_field_vect, &
                          b2_field_vect, &
                          c_field_scal)

print*, 'solve the elliptic solver'

call solve_general_coordinates_elliptic_eq_mp(&
     es_mp,&
     rho_field_scal,&
     phi_field_scal)

do ipatch= 0,num_patches-1
   m        => rho_field_scal%get_cartesian_mesh(ipatch)
   transf   => rho_field_scal%get_transformation(ipatch)
   num_pts1 = m%num_cells1+1
   num_pts2 = m%num_cells2+1
   delta1   = m%delta_eta1
   delta2   = m%delta_eta2
   eta1min  = m%eta1_min
   eta2min  = m%eta2_min
   
   do j=1,num_pts1
      eta2 = eta2min + (j-1)*delta2
      do i=1,num_pts2
      ! here it is assumed that the eta_min's are = 0. This is supposed
         ! to be the case for NURBS transformations.
         eta1 = eta1min + (i-1)*delta1
         x1 = transf%x1(eta1,eta2)
         x2 = transf%x2(eta1,eta2)
         
         val_phi  = phi_field_scal%value_at_point(eta1,eta2,ipatch)
         val_phi_exacte = exp(-9*(x1**2+x2**2))!sol_exacte_perdir(x1,x2)
         call phi_field_diff%set_value_at_indices( i, j, ipatch, abs(val_phi-val_phi_exacte) ) 
      end do
   end do
end do

print *, 'updating multipatch field interpolation coefficients...'
call phi_field_diff%update_interpolation_coefficients()

print *, 'writing to file...'
call phi_field_scal%write_to_file(0)

print *, 'writing to file...'
call phi_field_ex%write_to_file(0)

print *, 'writing to file...'
call phi_field_diff%write_to_file(0)

print*, 'delete object'
call sll_delete(T)
call sll_delete(a11_field_mat)
call sll_delete(a21_field_mat)
call sll_delete(a12_field_mat)
call sll_delete(a22_field_mat)
call sll_delete(b1_field_vect)
call sll_delete(b2_field_vect)
call sll_delete(c_field_scal)
call sll_delete(rho_field_scal)
call sll_delete(phi_field_scal)
call sll_delete(phi_field_ex)
call sll_delete(phi_field_diff)
print *, "PASSED"

end program test_general_elliptic_solver_multipatch

! External functions used as parameters in the above unit test:

function func_one( eta1, eta2) result(res)
  real(8), intent(in) :: eta1
  real(8), intent(in) :: eta2
  real(8) :: res
  res = 1.0_8
end function func_one

function func_zero( eta1, eta2) result(res)
  real(8), intent(in) :: eta1
  real(8), intent(in) :: eta2
  real(8) :: res
  res = 0.0_8
end function func_zero

function func_epsi( eta1, eta2, params ) result(res)
  real(8), intent(in) :: eta1
  real(8), intent(in) :: eta2
  real(8), dimension(:), intent(in) :: params
  real(8) :: res

  res = 0.0_8
end function func_epsi

!----------------------------------------------------------
!  Solution for a identity change of coordinates 
!   and periodic-periodic conditions
!   the matrix A is equal to identity 
!   the scalar c is equal to zero 
!-------------------------------------------------------------
function source_term_perper( eta1, eta2) result(res)
  use sll_constants
  intrinsic :: cos

  real(8), intent(in) :: eta1
  real(8), intent(in) :: eta2
  ! real(8), dimension(:), intent(in), optional :: params
  real(8) :: res

  res =  0.001*cos(2*sll_pi*eta1)
  !!-2*(2.0*sll_pi)**2*cos(2.0*sll_pi*eta1)*cos(2.0*sll_pi*eta2)! 0.001*cos(2*sll_pi*eta1)!
end function source_term_perper

real(8) function sol_exacte_perper(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  
  !real(8), dimension(:), intent(in), optional :: params
  sol_exacte_perper = -0.001/((2*sll_pi)**2)*cos(2*sll_pi*eta1)!cos(2.0*sll_pi*eta1)*cos(2.0*sll_pi*eta2)!-0.001/((2*sll_pi)**2)*cos(2*sll_pi*eta1)
end function sol_exacte_perper

real(8) function sol_exacte_perper_der1(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  
  !real(8), dimension(:), intent(in), optional :: params
  sol_exacte_perper_der1 = 0.001/(2*sll_pi)*sin(2*sll_pi*eta1) !-2.0*sll_pi*sin(2.0*sll_pi*eta1)*cos(2.0*sll_pi*eta2)
end function sol_exacte_perper_der1
real(8) function sol_exacte_perper_der2(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  
  !real(8), dimension(:), intent(in), optional :: params
  sol_exacte_perper_der2 = 0.0_f64!-2.0*sll_pi*cos(2.0*sll_pi*eta1)*sin(2.0*sll_pi*eta2)
end function sol_exacte_perper_der2

!----------------------------------------------------------
!  Solution for a identity change of coordinates 
!   and periodic-dirichlet conditions
!   and also dirichlet-dirichlet conditons
!   the matrix A is equal to identity 
!   the scalar c is equal to zero 
!-------------------------------------------------------------

real(8) function source_term_perdir(eta1,eta2) ! in the path
  use sll_constants
  intrinsic :: cos
  intrinsic :: sin 
  real(8),intent(in) :: eta1,eta2

  
  source_term_perdir = -2*(0.5*sll_pi)**2* sin(0.5*sll_pi*eta1)*sin(0.5*sll_pi*eta2)
      ! -(16.0*sll_pi**2*eta2**4 &
      ! - 16.0*sll_pi**2*eta2**2 &
      ! - 12.0*eta2**2 + 2.0)*cos(2*sll_pi*eta1)*sin(2*sll_pi*eta1)
  
end function source_term_perdir


real(8) function sol_exacte_perdir(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  intrinsic :: cos
  intrinsic :: sin
  !real(8), dimension(:), intent(in), optional :: params
  sol_exacte_perdir = sin(0.5*sll_pi*eta1)*sin(0.5*sll_pi*eta2)!eta2 ** 2 * (eta2**2-1)&
      ! * cos(2.0*sll_pi*eta1)*sin(2.0*sll_pi*eta1)
  
  !print*, 'heho'
end function sol_exacte_perdir


real(8) function sol_exacte_perdir_der1(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  intrinsic :: cos
  intrinsic :: sin
  !real(8), dimension(:), intent(in), optional :: params
  sol_exacte_perdir_der1 = 2.0*sll_pi*cos(2.0*sll_pi*eta1)*sin(2.0*sll_pi*eta2)
end function sol_exacte_perdir_der1


real(8) function sol_exacte_perdir_der2(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  intrinsic :: cos
  intrinsic :: sin
  !real(8), dimension(:), intent(in), optional :: params
  sol_exacte_perdir_der2 = 2.0*sll_pi*sin(2.0*sll_pi*eta1)*cos(2.0*sll_pi*eta2)
end function sol_exacte_perdir_der2

!  Solution for a identity change of coordinates 
!   and also dirichlet-periodicconditons
!   the matrix A is equal to identity 
!   the scalar c is equal to zero 
!-------------------------------------------------------------

real(8) function source_term_dirper(eta1,eta2,params) ! in the path
  use sll_constants
  real(8),intent(in) :: eta1,eta2
  real(8), dimension(:), intent(in), optional :: params

  source_term_dirper = -2*(2*sll_pi)**2* sin(2*sll_pi*eta1)*cos(2*sll_pi*eta2)
     ! -(16.0*sll_pi**2*eta1**4 &
     ! - 16.0*sll_pi**2*eta1**2 &
     ! - 12.0*eta1**2 + 2.0)*sin(2*sll_pi*eta2)*cos(2*sll_pi*eta2)
end function source_term_dirper


real(8) function sol_exacte_dirper(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  !real(8), dimension(:), intent(in), optional :: params
  sol_exacte_dirper = sin(2.0*sll_pi*eta1)*cos(2.0*sll_pi*eta2)
       !eta1 ** 2 * (eta1**2-1)* cos(2*sll_pi*eta2)*sin(2*sll_pi*eta2)
  
  
end function sol_exacte_dirper

real(8) function sol_exacte_dirper_der1(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  !real(8), dimension(:), intent(in), optional :: params
  sol_exacte_dirper_der1 = 2*sll_pi*cos(2.0*sll_pi*eta1)*cos(2.0*sll_pi*eta2)
       !eta1 ** 2 * (eta1**2-1)* cos(2*sll_pi*eta2)*sin(2*sll_pi*eta2)
  
  
end function sol_exacte_dirper_der1

real(8) function sol_exacte_dirper_der2(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  !real(8), dimension(:), intent(in), optional :: params
  sol_exacte_dirper_der2 = -2.0*sll_pi*sin(2.0*sll_pi*eta1)*sin(2.0*sll_pi*eta2)
       !eta1 ** 2 * (eta1**2-1)* cos(2*sll_pi*eta2)*sin(2*sll_pi*eta2)
  
  
end function sol_exacte_dirper_der2

!----------------------------------------------------------
!  Solution for a r theta change of coordinates 
!   and periodic-dirichlet conditions
!   and also dirivhlet-dirichlet conditons
!   the matrix A is equal to identity 
!   the scalar c is equal to zero 
!-------------------------------------------------------------


real(8) function rho_rtheta(eta1,eta2,params) ! in the path
  use sll_constants
  intrinsic :: cos
  intrinsic :: sin
  real(8),intent(in) :: eta1,eta2
  real(8) :: x, y
  real(8), dimension(:), intent(in), optional :: params
  
  if (present(params)) print*, params

  x = eta2*cos(2*sll_pi*eta1)
  y = eta2*sin(2*sll_pi*eta1)
  
  rho_rtheta = x*y*(-32.0*x**2 - 32.0*y**2 + 15.0)  
  
  
end function rho_rtheta


real(8) function sol_exacte_rtheta(eta1,eta2,params) ! in the path
  use sll_constants
  real(8),intent(in) :: eta1,eta2
  intrinsic :: cos
  intrinsic :: sin
  real(8), dimension(:), intent(in), optional :: params
  
  if (present(params)) print*, params
  
  sol_exacte_rtheta = ( eta2**2-1)*(eta2**2-0.5**2)*eta2**2&
       *cos(2*sll_pi*eta1)*sin(2*sll_pi*eta1)
  
end function sol_exacte_rtheta

!----------------------------------------------------------
!  Solution for a colella change of coordinates 
!   and periodic-periodic conditons
!   the matrix A is equal to identity 
!   the scalar c is equal to zero 
!-------------------------------------------------------------

real(8) function source_term_chgt_perper(eta1,eta2) ! in the path
  use sll_constants
  intrinsic :: cos
  intrinsic :: sin
  real(8):: eta1,eta2
  real(8) :: x, y
  ! real(8), dimension(:), intent(in), optional :: params
  
  x =   eta1 + 0.1_f64*sin(2*sll_pi*eta1) * sin(2*sll_pi*eta2)
  y =   eta2 + 0.1_f64*sin(2*sll_pi*eta1) * sin(2*sll_pi*eta2)
  
  source_term_chgt_perper = -8.0*sll_pi**2*cos(2*sll_pi*x)*cos(2*sll_pi*y) 
  
end function source_term_chgt_perper

real(8) function sol_exacte_chgt_perper(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  real(8) :: x,y
  intrinsic :: cos
  intrinsic :: sin
  
  x =   eta1 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  y =   eta2 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  
  sol_exacte_chgt_perper = cos(2*sll_pi*x)*cos(2*sll_pi*y)
  
  
end function sol_exacte_chgt_perper

real(8) function sol_exacte_chgt_perper_der1(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  real(8) :: x,y
  intrinsic :: cos
  intrinsic :: sin
  
  x =   eta1 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  y =   eta2 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  
  sol_exacte_chgt_perper_der1 = -2*sll_pi*sin(2*sll_pi*x)*cos(2*sll_pi*y)&
       * ( 1.0_f64 + 0.1*2*sll_pi*cos(2* sll_pi*eta1) * sin(2*sll_pi*eta2) )&
       -2*sll_pi*cos(2*sll_pi*x)*sin(2*sll_pi*y)&
       * ( 0.1*2*sll_pi*cos(2* sll_pi*eta1) * sin(2*sll_pi*eta2) )
  
end function sol_exacte_chgt_perper_der1

real(8) function sol_exacte_chgt_perper_der2(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  real(8) :: x,y
  intrinsic :: cos
  intrinsic :: sin
  
  x =   eta1 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  y =   eta2 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  
  !sol_exacte_chgt_perper_der2 = -2*sll_pi*cos(2*sll_pi*x)*sin(2*sll_pi*y)
  
  sol_exacte_chgt_perper_der2 = -2*sll_pi*sin(2*sll_pi*x)*cos(2*sll_pi*y)&
       * ( 0.1*2*sll_pi*sin(2* sll_pi*eta1) * cos(2*sll_pi*eta2) )&
       -2*sll_pi*cos(2*sll_pi*x)*sin(2*sll_pi*y)&
       * ( 1.0_f64 + 0.1*2*sll_pi*sin(2* sll_pi*eta1)*cos(2*sll_pi*eta2) )
end function sol_exacte_chgt_perper_der2


!----------------------------------------------------------
!  Solution for a colella change of coordinates 
!   and periodic-dirichlet conditons
!   and dircihlet-diichlet conditions
!   the matrix A is equal to identity 
!   the scalar c is equal to zero 
!-------------------------------------------------------------


real(8) function source_term_chgt_perdir(eta1,eta2) ! in the path
  use sll_constants
  real(8),intent(in) :: eta1,eta2
  real(8) :: x, y
  intrinsic :: cos
  intrinsic :: sin
  !real(8), dimension(:), intent(in), optional :: params
    
  x =   eta1 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  y =   eta2 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  
  source_term_chgt_perdir= -2*(2*sll_pi)**2 * sin(2*sll_pi*y)*cos(2*sll_pi*x)
  
  
end function source_term_chgt_perdir

real(8) function sol_exacte_chgt_perdir(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  real(8) :: x,y
  intrinsic :: cos
  intrinsic :: sin
 ! real(8), dimension(:), intent(in), optional :: params
  x =   eta1 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  y =   eta2 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  
  sol_exacte_chgt_perdir = cos(2*sll_pi*x)*sin(2*sll_pi*y)
  
  
end function sol_exacte_chgt_perdir

real(8) function sol_exacte_chgt_perdir_der1(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  real(8) :: x,y
  intrinsic :: cos
  intrinsic :: sin
 ! real(8), dimension(:), intent(in), optional :: params
  x =   eta1 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  y =   eta2 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  
  sol_exacte_chgt_perdir_der1 = -2*sll_pi*sin(2*sll_pi*x)*sin(2*sll_pi*y)&
       * ( 1.0_f64 + 0.1*2*sll_pi*cos(2*sll_pi*eta1) * sin(2*sll_pi*eta2) )&
       + 2*sll_pi*cos(2*sll_pi*x)*cos(2*sll_pi*y)&
       * ( 2*sll_pi*0.1*cos(2* sll_pi*eta1) * sin(2*sll_pi*eta2) ) 
  
  
end function sol_exacte_chgt_perdir_der1

real(8) function sol_exacte_chgt_perdir_der2(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  real(8) :: x,y
  intrinsic :: cos
  intrinsic :: sin
 ! real(8), dimensi@on(:), intent(in), optional :: params
  x =   eta1 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  y =   eta2 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  
  sol_exacte_chgt_perdir_der2 = -2*sll_pi*sin(2*sll_pi*x)*sin(2*sll_pi*y)&
       * ( 0.1*2*sll_pi*sin(2*sll_pi*eta1) * cos(2*sll_pi*eta2) ) &
       + 2*sll_pi*cos(2*sll_pi*x)*cos(2*sll_pi*y)&
       * ( 1.0_f64 + 2*sll_pi*0.1*sin(2* sll_pi*eta1) *cos(2*sll_pi*eta2) ) 
  
  
end function sol_exacte_chgt_perdir_der2

!----------------------------------------------------------
!  Solution for a colella change of coordinates 
!   and dirchlet-periodic conditions
!   the matrix A is equal to identity 
!   the scalar c is equal to zero 
!-------------------------------------------------------------



real(8) function source_term_chgt_dirdir(eta1,eta2) ! in the path
  use sll_constants
  real(8),intent(in) :: eta1,eta2
  real(8) :: x, y
  intrinsic :: cos
  intrinsic :: sin
  ! -------------------------------------------------
  ! In the case without change of coordinates
  ! -------------------------------------------------
  x =   eta1 + 0.1_8*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  y =   eta2 + 0.1_8*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  
  
  source_term_chgt_dirdir = &
       -2*(2.0*sll_pi)**2*sin(2*sll_pi*x)*sin(2*sll_pi*y)
  
end function source_term_chgt_dirdir




real(8) function sol_exacte_chgt_dirdir(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  real(8) :: x,y
  intrinsic :: cos
  intrinsic :: sin
  
  x =   eta1 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  y =   eta2 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  
  
  sol_exacte_chgt_dirdir = sin(2* sll_pi*y)*sin(2* sll_pi*x)
  
end function sol_exacte_chgt_dirdir


real(8) function sol_exacte_chgt_dirdir_der1(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  real(8) :: x,y
  intrinsic :: cos
  intrinsic :: sin
  
  x =   eta1 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  y =   eta2 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  
  
  sol_exacte_chgt_dirdir_der1 = 2*sll_pi*cos(2* sll_pi*x)*sin(2* sll_pi*y)&
       * ( 1.0_f64 + 0.1*2*sll_pi*cos(2*sll_pi*eta1) * sin(2*sll_pi*eta2) )&
       + 2*sll_pi*sin(2* sll_pi*x)*cos(2* sll_pi*y) &
       * ( 2*sll_pi*0.1*cos(2* sll_pi*eta1) * sin(2*sll_pi*eta2) )
end function sol_exacte_chgt_dirdir_der1


real(8) function sol_exacte_chgt_dirdir_der2(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  real(8) :: x,y
  intrinsic :: cos
  intrinsic :: sin
  
  x =   eta1 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  y =   eta2 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  
  
  sol_exacte_chgt_dirdir_der2 =  2*sll_pi*cos(2* sll_pi*x)*sin(2* sll_pi*y)&
       * ( 0.1*2*sll_pi*sin(2*sll_pi*eta1) * cos(2*sll_pi*eta2)  )&
       + 2*sll_pi*sin(2* sll_pi*x)*cos(2* sll_pi*y) &
       * ( 1.0_f64 + 2*sll_pi*0.1*sin(2* sll_pi*eta1) *cos(2*sll_pi*eta2) )
  
end function sol_exacte_chgt_dirdir_der2

!----------------------------------------------------------
!  Solution for a colella change of coordinates 
!   and dirchlet-periodic conditions
!   the matrix A is equal to identity 
!   the scalar c is equal to zero 
!-------------------------------------------------------------

real(8) function source_term_chgt_dirper(eta1,eta2) ! in the path
  use sll_constants
  real(8),intent(in) :: eta1,eta2
  real(8) :: x, y
  intrinsic :: cos
  intrinsic :: sin
  !real(8), dimension(:), intent(in), optional :: params
  ! -------------------------------------------------
  ! In the case without change of coordinates
  ! -------------------------------------------------
  x =   eta1 + 0.1_8*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  y =   eta2 + 0.1_8*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  
  
  source_term_chgt_dirper = -2*(2*sll_pi)**2*sin(2*sll_pi*x)*cos(2*sll_pi*y)
  
end function source_term_chgt_dirper

real(8) function sol_exacte_chgt_dirper(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  real(8) :: x,y
  intrinsic :: cos
  intrinsic :: sin

  x =   eta1 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  y =   eta2 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  
  sol_exacte_chgt_dirper = sin(2* sll_pi*x)*cos(2* sll_pi*y)
  
end function sol_exacte_chgt_dirper

real(8) function sol_exacte_chgt_dirper_der1(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  real(8) :: x,y
  intrinsic :: cos
  intrinsic :: sin
  
  x =   eta1 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  y =   eta2 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  
  
  sol_exacte_chgt_dirper_der1 = 2*sll_pi*cos(2* sll_pi*x)*cos(2* sll_pi*y) &
       * ( 1.0_f64 + 0.1*2*sll_pi*cos(2*sll_pi*eta1) * sin(2*sll_pi*eta2) )&
       - 2*sll_pi*sin(2* sll_pi*x)*sin(2* sll_pi*y)&
       * ( 2*sll_pi*0.1*cos(2* sll_pi*eta1) * sin(2*sll_pi*eta2) ) 
end function sol_exacte_chgt_dirper_der1

real(8) function sol_exacte_chgt_dirper_der2(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  real(8) :: x,y
  intrinsic :: cos
  intrinsic :: sin
  
  x =   eta1 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  y =   eta2 + 0.1*sin(2* sll_pi*eta1) * sin(2*sll_pi*eta2)
  
  
  sol_exacte_chgt_dirper_der2 = 2*sll_pi*cos(2* sll_pi*x)*cos(2* sll_pi*y) &
       * ( 0.1*2*sll_pi*sin(2*sll_pi*eta1) * cos(2*sll_pi*eta2)  )&
       - 2*sll_pi*sin(2* sll_pi*x)*sin(2* sll_pi*y)&
       * (1.0_f64 + 2*sll_pi*0.1*sin(2* sll_pi*eta1) *cos(2*sll_pi*eta2) ) 
  
end function sol_exacte_chgt_dirper_der2


!!!!!! test case with F(theta,phi) = (2pi theta , 2pi phi)

real(8) function adimension_chgt_x(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  print*, eta1, eta2
  adimension_chgt_x = 2*sll_pi*eta1 !+ eta2)
end function adimension_chgt_x

real(8) function adimension_chgt_y(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  print*, eta1, eta2
  adimension_chgt_y = 2*sll_pi*eta2
end function adimension_chgt_y

real(8) function jac11_adimension_chgt(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  print*, eta1, eta2
  jac11_adimension_chgt = 2*sll_pi
end function jac11_adimension_chgt

real(8) function jac12_adimension_chgt(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  print*, eta1, eta2
  jac12_adimension_chgt = 0.0!sll_pi
end function jac12_adimension_chgt

real(8) function jac21_adimension_chgt(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  print*, eta1, eta2
  jac21_adimension_chgt = 0.0!2*sll_pi!0.0
end function jac21_adimension_chgt

real(8) function jac22_adimension_chgt(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  print*, eta1, eta2
  jac22_adimension_chgt = 2*sll_pi
end function jac22_adimension_chgt

real(8) function sol_exacte_chgt_adim(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  real(8) :: x,y
  intrinsic :: cos
  intrinsic :: sin
  
  x =   2*sll_pi*eta1!+eta2)
  y =   2* sll_pi*eta2
  
  sol_exacte_chgt_adim = cos(x)*cos(y)
  
end function sol_exacte_chgt_adim

real(8) function source_term_chgt_adim(eta1,eta2)
  use sll_constants
  real(8) :: eta1,eta2
  real(8) :: x,y
  intrinsic :: cos
  intrinsic :: sin
  
  x =   2*sll_pi*eta1 !+eta2)
  y =   2* sll_pi*eta2
  
  source_term_chgt_adim = -2*cos(x)*cos(y)
  
end function source_term_chgt_adim
