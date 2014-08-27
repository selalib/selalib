program test_box_splines

#include "sll_constants.h"
#include "sll_memory.h"
#include "sll_working_precision.h"
#include "sll_assert.h"

use sll_utilities
use sll_constants
use hex_mesh
use box_splines
implicit none

type(hex_mesh_2d),   pointer  :: mesh
type(box_spline_2D), pointer :: spline
sll_int32    :: num_cells
sll_int32    :: i
sll_int32    :: deg = 3
sll_int32    :: nloops
sll_int32    :: ierr
! initial distribution
sll_real64   :: gauss_x2
sll_real64   :: gauss_x1
sll_real64   :: gauss_sig
sll_real64   :: gauss_amp
sll_real64,dimension(:),allocatable :: x1
sll_real64,dimension(:),allocatable :: x2
sll_real64,dimension(:),allocatable :: f_init
sll_real64,dimension(:),allocatable :: chi1
sll_real64,dimension(:),allocatable :: chi2
sll_real64,dimension(:),allocatable :: chi3
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
sll_real64    :: cfl
sll_real64   :: f_min
sll_real64   :: x2_basis
sll_real64   :: x1_basis
sll_real64   :: x1_temp
sll_int32    :: k1_error
sll_int32    :: k2_error
! Output variables
sll_int32    :: WRITE_TIME_ERROR = 0
sll_int32    :: WRITE_TIME_DIST = 0
sll_int32    :: WRITE_SPLINES = 0
sll_int32    :: WRITE_CELLS_ERROR = 1
sll_int32    :: WRITE_TIME_EFF = 0
character(len = 50) :: filename
character(len = 50) :: filename2
character(len = 4)  :: filenum

! sll_real64              :: delta
! sll_real64              :: k1_temp
! sll_real64              :: k2_temp



print*, " ---- BEGIN test_box_splines.F90 -----"
print*, ""
do num_cells = 80,80,20


   ! Mesh initialization
   mesh => new_hex_mesh_2d(num_cells, 0._f64, 0._f64, &
                                    !0.5_f64, -sqrt(3.)/2._f64, &
                                    !0.5_f64,  sqrt(3.)/2._f64, &
                                    !1.0_f64,  0._f64, &
                                    radius = 8._f64)
  
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
   SLL_ALLOCATE(chi3(mesh%num_pts_tot),ierr)


   ! Distribution initialization
   gauss_x1  = 2._f64
   gauss_x2  = 2._f64
   gauss_sig = 1.0_f64/sqrt(2._f64)/2._f64
   gauss_amp = 1.0_f64

   do i=1, mesh%num_pts_tot
      x1(i) = mesh%global_to_x1(i)
      x2(i) = mesh%global_to_x2(i)
      f_init(i) = gauss_amp * &
           exp(-0.5_f64*((x1(i)-gauss_x1)**2 / gauss_sig**2 &
           + (x2(i)-gauss_x2)**2 / gauss_sig**2))
      if (exponent(f_init(i)) .lt. -17) then
         f_init(i) = 0._f64
      end if
      f_tn(i) = f_init(i)
   end do

!   call write_field_hex_mesh(mesh, f_init, "init_dist.txt")

   ! Advection initialization
   ! if : which_advec = 0 => linear advection
   ! if : which_advec = 1 => circular advection
   which_advec = 1
   advec = 0.0_f64!25_f64!5_f64
   tmax  = 1.0_f64
   dt    = 0.1_f64 * 20._f64/num_cells
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

   spline => new_box_spline_2d(mesh, SLL_DIRICHLET)

   call cpu_time(t_init)

   !print*,""
   do while (t .lt. tmax)
      !print *, " --> Time loop t =", t
      !Error variables
      norm2_error = 0._f64
      diff_error  = 0._f64
      where_error = -1

      nloops = nloops + 1

      call compute_box_spline_2d( f_tn, deg, spline )
      t      = t + dt

      do i=1, mesh%num_pts_tot

         ! ******************
         ! Approximation
         ! ******************
         f_tn(i) = hex_interpolate_value(mesh, x1_char(i), x2_char(i), spline, deg)
         ! if (f_tn(i) < 0.) then
!             print *, "Negative value at"
!             print *, "      Value :", f_tn(i)
!             print *, "      Time  :", t
!             print *, "      Loop  :", nloops
!             print *, "      Point :", i
!             print *, "      X1char:", x1_char(i)
!             print *, "      X2char:", x2_char(i)
!          !   STOP
!          !   f_tn(i) = 0._f64
         ! end if
         
         ! ******************
         ! Analytical value 
         ! ******************
         ! Computing characteristics
         ! if (which_advec .eq. 0) then
!             ! linear advection
!             x1(i) = mesh%global_to_x1(i) - advec*dt*nloops
!             x2(i) = mesh%global_to_x2(i) - advec*dt*nloops
!          else
            ! Circular advection
            x1_temp = sqrt(x1(i)**2 + x2(i)**2) * cos(2*sll_pi*dt + atan2(x2(i),x1(i)))
            x2(i)   = sqrt(x1(i)**2 + x2(i)**2) * sin(2*sll_pi*dt + atan2(x2(i),x1(i)))
            x1(i)   = x1_temp
         ! end if

         f_fin(i) = gauss_amp * &
              exp(-0.5_f64*((x1(i)-gauss_x1)**2/gauss_sig**2 &
              + (x2(i)-gauss_x2)**2 / gauss_sig**2))
         if (exponent(f_fin(i)) .lt. -17) then
            f_fin(i) = 0._f64
         end if

         x1_basis = change_basis_x1(spline, x1(i), x2(i))
         x2_basis = change_basis_x2(spline, x1(i), x2(i))
         chi1(i) = chi_gen_val(x1_basis, x2_basis, 1)
         chi2(i) = chi_gen_val(x1_basis, x2_basis, 2)
         chi3(i) = chi_gen_val(x1_basis, x2_basis, 3)
         
            ! Relative error
            if (diff_error .lt. abs(f_fin(i) - f_tn(i)) ) then
               diff_error = abs(f_fin(i) - f_tn(i))
               where_error = i
            end if
            ! Norm2 error :
            norm2_error = norm2_error + abs(f_fin(i) - f_tn(i))**2
         
         end do

!    if (WRITE_SPLINES.eq.1) then 
         call write_field_hex_mesh_xmf(mesh, chi1, "chi1")
         call write_field_hex_mesh_xmf(mesh, chi2, "chi2")
         call write_field_hex_mesh_xmf(mesh, chi3, "chi3")
!    end if

      ! Norm2 error :
      norm2_error = sqrt(norm2_error)
      
      ! Printing error
      k1_error = mesh%global_to_hex1(where_error)
      k2_error = mesh%global_to_hex2(where_error)
      print*,"  nt =", nloops, "    | error_Linf = ", diff_error
      print*,"                       | at hex =", cells_to_origin(k1_error, k2_error), where_error
      print*,"                       | error_L2   = ", norm2_error
      print*," Center error = ", f_fin(1)-f_tn(1)




      !WRITING ERROR REGARDING TIME STEP
!       if (WRITE_TIME_ERROR.eq.1) then 
!          if (t .eq. dt) then
!             open (unit=12,file="err_chi2_gauss_cstadv.txt", &
!                  action="write",status="replace")
!             write (12, "(3(g13.3,1x))") t, diff_error, norm2_error
!             close(12)
!          else 
!             open (unit=12,file="err_chi2_gauss_cstadv.txt", &
!                  action="write",status="old", position="append")
!             write (12, "(3(g13.3,1x))") t, diff_error, norm2_error
!             close(12)
!          end if
!       end if
       

!       if (WRITE_TIME_DIST.eq.1) then 
!          call int2string(nloops,filenum)
!          filename2 = "./time_files/analytical/ana_dist"//trim(filenum)!//".txt"
!          filename  = "./time_files/numerical/num_dist"//trim(filenum)!//".txt"
!          print*,filename
!          print*,filename2
!          call write_field_hex_mesh_xmf(mesh, f_tn, trim(filename))
!          call write_field_hex_mesh_xmf(mesh, f_fin, trim(filename2))
!       end if

   end do


   call cpu_time(t_end)
   print*," *    Final error  = ", diff_error, " *"

   cfl = dt * num_cells
   f_min = minval(f_tn)
   !WRITING ERROR REGARDING NUMBER OF POINTS
   if (WRITE_CELLS_ERROR.eq.1) then
      if (num_cells .eq. 20) then 
         !NEW FILE :
         open (unit=12,file="error_file.txt",action="write",&
              status="replace")
         !write (12, "(3(g13.3,1x))") num_cells, diff_error, norm2_error
         write(12,*) num_cells, mesh%num_pts_tot, dt, cfl,  norm2_error, diff_error, f_min, t_end - t_init,&
                nloops * 3._f64*real(num_cells + 1, f64)*real(num_cells, f64)/(t_end - t_init)/ 1e6_f64   
      close(12)
      else
         !WRITE
         open (unit=12,file="error_file.txt",action="write",&
              status="old", position="append") 
!         write (12, "(3(g13.3,1x))") num_cells, diff_error, norm2_error
         write(12,*) num_cells, mesh%num_pts_tot, dt, cfl,  norm2_error, diff_error, f_min, t_end - t_init,&
                nloops * 3._f64*real(num_cells + 1, f64)*real(num_cells, f64)/(t_end - t_init)/ 1e6_f64
         close(12)
      end if
   end if


!    if (WRITE_TIME_EFF.eq.1) then

!       !WRITING CPU TIME
!       if (num_cells .eq. 10 ) then 
!          !NEW FILE :
!          open (unit=12,file="cpu_time.txt",action="write",&
!               status="replace")
!          write (12, "(3(G15.3,1x))") num_cells, mesh%num_pts_tot, t_end-t_init
!          close(12)
!       else
!          !WRITE
!          open (unit=12,file="cpu_time.txt",action="write",&
!               status="old", position="append") 
!          write (12, "(3(G15.3,1x))") num_cells, mesh%num_pts_tot, t_end-t_init
!          close(12)
!       end if
!    end if



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

