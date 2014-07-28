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
sll_int32    :: which_advec
sll_real64   :: advec
sll_real64   :: dt
sll_real64   :: tmax
sll_real64   :: t
! timer variables
sll_real64   :: tcpu, t_init, t_end
!others
sll_real64   :: x2_basis
sll_real64   :: x1_basis
sll_real64   :: x1_temp
character(len = 50) :: filename
character(len = 50) :: filename2
character(len = 4) :: filenum

! sll_real64              :: delta
! sll_real64              :: k1_temp
! sll_real64              :: k2_temp



print*, ""
do num_cells = 70,70,1

!   t_init = getRealTimer()

   ! Mesh initialization    
   mesh => new_hex_mesh_2d(num_cells, 0._f64, 0._f64, &
                                    !0.5_f64, -sqrt(3.)/2._f64, &
                                    !0.5_f64,  sqrt(3.)/2._f64, &
                                    !1.0_f64,  0._f64, &
                                    radius = 1._f64)
  call sll_display(mesh)
  print*,"num_pts : ", mesh%num_pts_tot
  print*,"spl deg : ", deg
  print*,""

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
   gauss_x1  = -0.25_f64
   gauss_x2  = -0.25_f64
   gauss_sig = 0.05_f64

   do i=0, mesh%num_pts_tot-1
      x1(i+1) = from_global_x1(mesh, i)
      x2(i+1) = from_global_x2(mesh, i)
      f_init(i+1) = 1._f64*exp(-0.5_f64*((x1(i+1)-gauss_x1)**2 / gauss_sig**2 &
                    + (x2(i+1)-gauss_x2)**2 / gauss_sig**2))
      if (exponent(f_init(i+1)) .lt. -17) then
         f_init(i+1) = 0._f64
      end if
      f_tn(i+1) = f_init(i+1)
   end do

!   call write_field_hex_mesh(mesh, f_init, "init_dist.txt")

   ! Advection initialization
   which_advec = 0
   advec = 0.025_f64!5_f64
   tmax  = 2.5_f64
   dt    = 0.025_f64
   t     = 0._f64

   ! Computing characteristics
   if (which_advec .eq. 0) then
      ! linear advection
      x1_char(:) = x1(:) - advec*dt
      x2_char(:) = x2(:) - advec*dt
   else
      ! Circular advection
      x1_char(1) = 0._f64
      x2_char(1) = 0._f64
      x1_char(2:) = sqrt(x1(2:)**2 + x2(2:)**2) * cos(2*sll_pi*dt + atan2(x2(2:),x1(2:)))
      x2_char(2:) = sqrt(x1(2:)**2 + x2(2:)**2) * sin(2*sll_pi*dt + atan2(x2(2:),x1(2:)))
   end if

   ! Time loop
   nloops = 0
   spline => new_linear_box_spline_2d(mesh, SLL_NEUMANN, SLL_NEUMANN)

   call cpu_time(t_init)
   print*,""
   do while (t .lt. tmax)
      !Error variables
      norm2_error = 0._f64
      diff_error  = 0._f64
      where_error = -1

      nloops = nloops + 1
      call compute_linear_box_spline_2d( f_tn, deg, spline )
      t      = t + dt
      do i=1, mesh%num_pts_tot

         ! ******************
         ! Approximation
         ! ******************
         f_tn(i) = hex_interpolate_value(mesh, x1_char(i), x2_char(i), spline, deg)
         ! if (f_tn(i) < 0.) then
         !    ! print *, "Negative value at"
         !    ! print *, "      Time  :", t
         !    ! print *, "      Loop  :", nloops
         !    ! print *, "      Point :", i
         !    ! print *, "      X1char:", x1_char(i)
         !    ! print *, "      X2char:", x2_char(i)
         !    ! STOP
         !    f_tn(i) = 0._f64
         ! end if

         ! ******************
         ! Analytical value 
         ! ******************
         ! Computing characteristics
         if (which_advec .eq. 0) then
            ! linear advection
            x1(i) = from_global_x1(mesh, i-1) - advec*dt*nloops
            x2(i) = from_global_x2(mesh, i-1) - advec*dt*nloops
         else
            ! Circular advection
            x1_temp = sqrt(x1(i)**2 + x2(i)**2) * cos(2*sll_pi*dt + atan2(x2(i),x1(i)))
            x2(i)   = sqrt(x1(i)**2 + x2(i)**2) * sin(2*sll_pi*dt + atan2(x2(i),x1(i)))
            x1(i)   = x1_temp
         end if

         f_fin(i) = exp(-0.5_f64*((x1(i)-gauss_x1)**2/gauss_sig**2 &
                    + (x2(i)-gauss_x2)**2 / gauss_sig**2))
         if (exponent(f_fin(i)) .lt. -17) then
            f_fin(i) = 0._f64
         end if

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
      
      call cpu_time(t_end)
      ! Printing error
      print*,"  nt =", nloops, "    | error_Linf = ", diff_error
      print*,"                       | at hex =", get_hex_num(where_error), where_error
      print*,"                       | error_L2   = ", norm2_error
      ! print*," Center error = ", f_fin(1)-f_tn(1)




      !WRITING ERROR REGARDING TIME STEP
      if (t .eq. dt) then
         open (unit=12,file="err_chi2_gauss_cstadv.txt", &
              action="write",status="replace")
         write (12, "(3(g13.3,1x))") t, diff_error, norm2_error
         close(12)
      else 
         open (unit=12,file="err_chi2_gauss_cstadv.txt", &
              action="write",status="old", position="append")
         write (12, "(3(g13.3,1x))") t, diff_error, norm2_error
         close(12)
      end if
          
      call int2string(nloops,filenum)
      filename2 = "./time_files/analytical/ana_dist"//trim(filenum)//".txt"
      filename  = "./time_files/numerical/num_dist"//trim(filenum)//".txt"
      print*,filename
      print*,filename2
      call write_field_hex_mesh(mesh, f_tn, filename)
      call write_field_hex_mesh(mesh, f_fin,filename2)

   end do


   call write_field_hex_mesh(mesh, f_tn, "final_dist.txt")
   call write_field_hex_mesh(mesh, chi1, "chi1.txt")
   call write_field_hex_mesh(mesh, chi2, "chi2.txt")
   call write_field_hex_mesh(mesh, f_fin, "an_dist.txt")
   print *,""
   print*," *    Final error  = ", diff_error, " *"


   ! !WRITING ERROR REGARDING NUMBER OF POINTS
   ! if (num_cells .eq. 10) then 
   !    !NEW FILE :
   !    open (unit=12,file="error_file.txt",action="write",&
   !         status="replace")
   !    write (12, "(3(g13.3,1x))") num_cells, diff_error, norm2_error
   !    close(12)
   ! else
   !    !WRITE
   !    open (unit=12,file="error_file.txt",action="write",&
   !         status="old", position="append") 
   !    write (12, "(3(g13.3,1x))") num_cells, diff_error, norm2_error
   !    close(12)
   ! end if


   !WRITING CPU TIME
   if (num_cells .eq. 10 ) then 
      !NEW FILE :
      open (unit=12,file="cpu_time.txt",action="write",&
           status="replace")
      write (12, "(3(G15.3,1x))") num_cells, mesh%num_pts_tot, t_end-t_init
      close(12)
   else
      !WRITE
      open (unit=12,file="cpu_time.txt",action="write",&
           status="old", position="append") 
      write (12, "(3(G15.3,1x))") num_cells, mesh%num_pts_tot, t_end-t_init
      close(12)
   end if



   SLL_DEALLOCATE_ARRAY(spline%coeffs,ierr)
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
