program test_box_splines

#include "sll_memory.h"
#include "sll_working_precision.h"

use sll_constants
use hex_mesh
use box_splines
implicit none

type(hex_mesh_2d), pointer  :: mesh
type(linear_box_spline_2D), pointer :: spline
sll_int32    :: num_cells
sll_int32    :: i
sll_int32    :: nloops
sll_int32    :: error
! initial distribution
sll_real64   :: gauss_x2
sll_real64   :: gauss_x1
sll_real64   :: gauss_sig
sll_real64,dimension(:),allocatable :: x1
sll_real64,dimension(:),allocatable :: x2
sll_real64,dimension(:),allocatable :: f_init
sll_real64,dimension(:),allocatable :: chi1
sll_real64,dimension(:),allocatable :: chi2
! distribution at time n
sll_real64,dimension(:),allocatable :: x1_char
sll_real64,dimension(:),allocatable :: x2_char
sll_real64,dimension(:),allocatable :: f_tn
! distribution at end time
sll_real64,dimension(:),allocatable :: f_fin
sll_real64   :: diff_error
! advection
sll_real64   :: advec
sll_real64   :: dt
sll_real64   :: tmax
sll_real64   :: t



! Mesh initialization 
 num_cells = 10
 mesh => new_hex_mesh_2d(num_cells)
 call sll_display(mesh)
 print *,""

! Distribution initialization
 gauss_x1  = 0.0
 gauss_x2  = 0.0
 gauss_sig = 0.1
 SLL_ALLOCATE(f_init(mesh%num_pts_tot),error)
 SLL_ALLOCATE(f_tn(mesh%num_pts_tot),error)
 SLL_ALLOCATE(f_fin(mesh%num_pts_tot),error)
 SLL_ALLOCATE(x1(mesh%num_pts_tot),error)
 SLL_ALLOCATE(x2(mesh%num_pts_tot),error)
 SLL_ALLOCATE(x1_char(mesh%num_pts_tot),error)
 SLL_ALLOCATE(x2_char(mesh%num_pts_tot),error)
 SLL_ALLOCATE(chi1(mesh%num_pts_tot),error)
 SLL_ALLOCATE(chi2(mesh%num_pts_tot),error)

 do i=0, mesh%num_pts_tot-1
      x1(i+1) = from_global_x1(mesh, i)
      x2(i+1) = from_global_x2(mesh, i)
      f_init(i+1) = 1.5*exp(-0.5*((x1(i+1)-gauss_x1)**2 / gauss_sig**2 &
                                + (x2(i+1)-gauss_x2)**2 / gauss_sig**2))
      f_tn(i+1) = f_init(i+1)
 end do

 call write_field_hex_mesh(mesh, f_init, "init_dist.txt")


! Advection initialization
 advec = 0.!1.5_f64
 tmax  = 0.025_f64
 dt    = 0.025_f64
 t = 0._f64

! Computing characteristics
 x1_char(:) = x1(:) - advec*dt
 x2_char(:) = x2(:) - advec*dt


! Time loop
 nloops = 0
 spline => new_linear_box_spline_2d(mesh, SLL_NEUMANN, SLL_NEUMANN)
 diff_error = 0._f64

! print*,"--", from_cart_index_k1(mesh,  0.00_f64,  0.00_f64), from_cart_index_k2(mesh,  0.00_f64,  0.00_f64)
! print*,"--", from_cart_index_k1(mesh,  0.01_f64,  0.01_f64), from_cart_index_k2(mesh,  0.01_f64,  0.01_f64)
! print*,"--", from_cart_index_k1(mesh, -0.01_f64,  0.01_f64), from_cart_index_k2(mesh, -0.01_f64,  0.01_f64)
! print*,"--", from_cart_index_k1(mesh, -0.01_f64, -0.01_f64), from_cart_index_k2(mesh, -0.01_f64, -0.01_f64)
! print*,"--", from_cart_index_k1(mesh, -0.05_f64,  0.05_f64), from_cart_index_k2(mesh, -0.05_f64,  0.05_f64)
! print*,"--", from_cart_index_k1(mesh, -0.05_f64, -0.05_f64), from_cart_index_k2(mesh, -0.05_f64, -0.05_f64)
! print*,"--", from_cart_index_k1(mesh,  0.05_f64,  0.05_f64), from_cart_index_k2(mesh,  0.05_f64,  0.05_f64)
! print*,"--", from_cart_index_k1(mesh, -0.09_f64,  0.09_f64), from_cart_index_k2(mesh, -0.09_f64,  0.09_f64)
! print*,"--", from_cart_index_k1(mesh, -0.09_f64, -0.09_f64), from_cart_index_k2(mesh, -0.09_f64, -0.09_f64)
! print*,"--", from_cart_index_k1(mesh,  0.09_f64,  0.09_f64), from_cart_index_k2(mesh,  0.09_f64,  0.09_f64)
!

 print*,""
 do while (t .lt. tmax)
      nloops = nloops + 1
      call compute_linear_box_spline_2d( f_tn, spline )
      t      = t + dt
      do i=1, mesh%num_pts_tot
         ! Approximation
         f_tn(i) = hex_interpolate_value(mesh, x1_char(i), x2_char(i), spline)
         ! Analytical value 
         x1(i) = from_global_x1(mesh, i-1) - advec*dt*nloops
         x2(i) = from_global_x2(mesh, i-1) - advec*dt*nloops
         f_fin(i) = 1.5*exp(-0.5*((x1(i)-gauss_x1)**2 / gauss_sig**2 &
                                + (x2(i)-gauss_x2)**2 / gauss_sig**2))
         chi1(i) = chi_gen_val(spline,x1(i), x2(i),1)
         chi2(i) = chi_gen_val(spline,x1(i), x2(i),2)
         ! Relative error
         if (diff_error .lt. abs(f_fin(i) - f_tn(i)) ) then
            diff_error = abs(f_fin(i) - f_tn(i))
         end if

      end do
      print*,"   n = ", nloops, "    | error  = ", diff_error
 end do

 print*,""
 print*, "val at origin of f_app : ", f_tn(1)
 print*, "val at origin of f_fin : ", f_fin(1)
 print*, "val at origin of diffn : ", f_fin(1) - f_tn(1)

 call write_field_hex_mesh(mesh, f_tn, "final_dist.txt")
 call write_field_hex_mesh(mesh, chi1, "chi1.txt")
 call write_field_hex_mesh(mesh, chi2, "chi2.txt")


 call write_field_hex_mesh(mesh, f_fin, "an_dist.txt")
 print *,""
 print*," *    Final error  = ", diff_error, " *"


end program test_box_splines
