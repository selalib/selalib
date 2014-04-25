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
sll_int32    :: deg = 2
sll_int32    :: nloops
sll_int32    :: ierr
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
sll_int32    :: where_error
sll_real64   :: diff_error
sll_real64   :: norm2_error
! advection
sll_real64   :: advec
sll_real64   :: dt
sll_real64   :: tmax
sll_real64   :: t
!others
sll_real64   :: x2_basis
sll_real64   :: x1_basis



print*, ""
do num_cells =50,50,10

   ! Mesh initialization    
   mesh => new_hex_mesh_2d(num_cells, 0._f64, 0._f64, &
                                    !0.5_f64, -sqrt(3.)/2._f64, &
                                    !0.5_f64,  sqrt(3.)/2._f64, &
                                    !1.0_f64,  0._f64, &
                                    radius = 1._f64)
   call sll_display(mesh)
   print*,"spl deg : ", deg
   print *,""

   SLL_ALLOCATE(f_init(mesh%num_pts_tot),ierr)
   SLL_ALLOCATE(f_tn(mesh%num_pts_tot),ierr)
   SLL_ALLOCATE(f_fin(mesh%num_pts_tot),ierr)
   SLL_ALLOCATE(x1(mesh%num_pts_tot),ierr)
   SLL_ALLOCATE(x2(mesh%num_pts_tot),ierr)
   SLL_ALLOCATE(x1_char(mesh%num_pts_tot),ierr)
   SLL_ALLOCATE(x2_char(mesh%num_pts_tot),ierr)
   SLL_ALLOCATE(chi1(mesh%num_pts_tot),ierr)
   SLL_ALLOCATE(chi2(mesh%num_pts_tot),ierr)

   ! Distribution initialization
   gauss_x1  = 0.0
   gauss_x2  = 0.0
   gauss_sig = 0.1

   do i=0, mesh%num_pts_tot-1
      x1(i+1) = from_global_x1(mesh, i)
      x2(i+1) = from_global_x2(mesh, i)
      f_init(i+1) = x2(i+1)**2! 1._f64*exp(-0.5_f64*((x1(i+1)-gauss_x1)**2 / gauss_sig**2 &
           !+ (x2(i+1)-gauss_x2)**2 / gauss_sig**2))

      f_tn(i+1) = f_init(i+1)
   end do
   call write_field_hex_mesh(mesh, f_init, "init_dist.txt")

   ! Advection initialization
   advec = 0.0!1_f64
   tmax  = 0.025_f64
   dt    = 0.025_f64
   t = 0._f64

   ! Computing characteristics
   x1_char(:) = x1(:) - advec*dt
   x2_char(:) = x2(:) - advec*dt

   ! Time loop
   nloops = 0
   spline => new_linear_box_spline_2d(mesh, SLL_NEUMANN, SLL_NEUMANN)
   diff_error  = 0._f64
   where_error = -1
   norm2_error = 0._f64

   print*,""
   do while (t .lt. tmax)
      nloops = nloops + 1
      call compute_linear_box_spline_2d( f_tn, deg, spline )
      t      = t + dt
      do i=1, mesh%num_pts_tot
         ! Approximation
         f_tn(i) = hex_interpolate_value(mesh, x1_char(i), x2_char(i), spline, deg)
         ! Analytical value 
         x1(i) = from_global_x1(mesh, i-1) + advec*dt*nloops
         x2(i) = from_global_x2(mesh, i-1) + advec*dt*nloops
         f_fin(i) = x2(i)**2! 1._f64*exp(-0.5_f64*((x1(i)-gauss_x1)**2/gauss_sig**2 &
                    !+ (x2(i)-gauss_x2)**2 / gauss_sig**2))
         x1_basis = change_basis_x1(mesh, spline, x1(i), x2(i))
         x2_basis = change_basis_x2(mesh, spline, x1(i), x2(i))
         chi1(i) = chi_gen_val(x1_basis, x2_basis, 1)
         chi2(i) = chi_gen_val(x1_basis, x2_basis, 2)
         if (get_hex_num(i) .lt. get_hex_num(mesh%num_pts_tot) - deg*2 - 1) then
            ! Relative error
            if (diff_error .lt. abs(f_fin(i) - f_tn(i)) ) then
               diff_error = abs(f_fin(i) - f_tn(i))
               where_error = i
            end if
            ! Norm2 error :
            norm2_error = norm2_error + abs(f_fin(i) - f_tn(i))**2
         end if
      end do

      ! Norm2 error :
      norm2_error = sqrt(norm2_error)

      ! Printing error
      print*,"   n = ", nloops, "    | error  = ", diff_error, &
      & "|    at hex = ", get_hex_num(where_error)
      print*,"                        | error  = ", norm2_error

      ! !WRITING ERROR REGARDING TIME STEP
      ! if (t .eq. dt) then
      !    open (unit=12,file="error_chi2_constadv_time.txt", &
      !        action="write",status="new")
      !        write (12, "(2(g13.3,1x))") t, diff_error
      !        close(12)
      !      else 
      !        open (unit=12,file="error_chi2_constadv_time.txt",&
      !        action="write",status="old", position="append")
      !        write (12, "(2(g13.3,1x))") t, diff_error
      !        close(12)
      !      end if
   end do
 

   print*,""
   print*," equises = ", x1(2), x1(3), x1(4)
   print*," igriega = ", x2(2), x2(3), x2(4)
   print*, "val at origin of f_app : ", f_tn(1)
   print*, "val at origin of f_fin : ", f_fin(1)
   print*, "val at origin of diffn : ", f_fin(1) - f_tn(1)
   print*,""
   ! print*, "val at maxerr of f_app : ", f_tn(where_error)
   ! print*, "val at maxerr of f_fin : ", f_fin(where_error)
   ! ! print*, "val at maxerr of chi1  : ", chi1(where_error)



   call write_field_hex_mesh(mesh, f_tn, "final_dist.txt")
   call write_field_hex_mesh(mesh, chi1, "chi1.txt")
   call write_field_hex_mesh(mesh, chi2, "chi2.txt")


   call write_field_hex_mesh(mesh, f_fin, "an_dist.txt")
   print *,""
   print*," *    Final error  = ", diff_error, " *"



   !WRITING ERROR REGARDING NUMBER OF POINTS
   if (num_cells .eq. 10) then 
      !NEW FILE :
      open (unit=12,file="error_file.txt",action="write",&
           status="new")
      write (12, "(3(g13.3,1x))") num_cells, diff_error, norm2_error
      close(12)
   else
      !WRITE
      open (unit=12,file="error_file.txt",action="write",&
           status="old", position="append") 
      write (12, "(3(g13.3,1x))") num_cells, diff_error, norm2_error
      close(12)
   end if

   SLL_DEALLOCATE_ARRAY(f_init,ierr)
   SLL_DEALLOCATE_ARRAY(f_tn,ierr)
   SLL_DEALLOCATE_ARRAY(f_fin,ierr)
   SLL_DEALLOCATE_ARRAY(x1,ierr)
   SLL_DEALLOCATE_ARRAY(x2,ierr)
   SLL_DEALLOCATE_ARRAY(chi1,ierr)
   SLL_DEALLOCATE_ARRAY(chi2,ierr)
   SLL_DEALLOCATE_ARRAY(x1_char,ierr)
   SLL_DEALLOCATE_ARRAY(x2_char,ierr)
   SLL_DEALLOCATE(mesh,ierr)
   SLL_DEALLOCATE(spline,ierr)


end do
end program test_box_splines
