program test_general_elliptic_solver_multipatch

#include "sll_memory.h"
#include "sll_working_precision.h"


use sll_m_cartesian_meshes
use sll_m_cartesian_meshes_multipatch
use sll_m_coordinate_transformations_2d
use sll_m_common_coordinate_transformations
use sll_m_coordinate_transformation_multipatch
use sll_m_scalar_field_2d
use sll_m_scalar_field_2d_multipatch
use sll_m_constants
use sll_m_arbitrary_degree_spline_interpolator_2d
use sll_m_timer
use sll_m_general_coordinate_elliptic_solver_multipatch


use m_multipatch_helper_functions, only : &
     func_zero, func_one

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
class(sll_cartesian_mesh_2d), pointer                      :: m
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
       val_rho  = -4.0_f64*9.0_f64*exp(-9*(x1**2+x2**2))+(2*9)**2*(x1**2+x2**2)*exp(-9*(x1**2+x2**2))
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

call solve_general_coordinates_elliptic_eq_mp( es_mp, rho_field_scal, phi_field_scal)

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
call phi_field_scal%write_to_file(1)

print *, 'writing to file...'
call phi_field_ex%write_to_file(1)

print *, 'writing to file...'
call phi_field_diff%write_to_file(1)

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

