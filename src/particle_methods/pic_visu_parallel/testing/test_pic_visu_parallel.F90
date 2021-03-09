program pic_visu_parallel

#include "sll_working_precision.h"

   use sll_m_collective
   use sll_m_pic_visu_parallel

   implicit none

   character(len=*), parameter       :: plot_name = "visu_pic_parallel_test_file"

   sll_int32, parameter     :: n = 5
   sll_real64, dimension(n) :: x
   sll_real64, dimension(n) :: v
   sll_real64, dimension(n) :: w
   sll_real64, parameter    :: xmin = 0.0_f64
   sll_real64, parameter    :: xmax = 1.0_f64
   sll_int32, parameter    :: nx = 8
   sll_real64, parameter    :: vmin = -1.0_f64
   sll_real64, parameter    :: vmax = 1.0_f64
   sll_int32, parameter    :: nv = 8
   sll_int32, parameter    :: iplot = 0
   sll_int32                :: root_rank = 0
   logical                  :: file_exists
   sll_int32                :: prank
   sll_int32                :: psize
   sll_int32                :: i, j

   call sll_s_boot_collective()
   prank = sll_f_get_collective_rank(sll_v_world_collective)
   psize = sll_f_get_collective_size(sll_v_world_collective)

   do i = 1, n
      x(i) = real(i - 1, f64)*(xmax - xmin)/real(nx - 1, f64) + xmin
   end do

   do j = 1, n
      v(j) = real(j - 1, f64)*(vmax - vmin)/real(nv - 1, f64) + vmin
   end do
   w = real(prank, f64)

   call sll_s_distribution_xdmf_coll(plot_name, x, v, w, &
                                     xmin, xmax, nx, vmin, vmax, nv, iplot, sll_v_world_collective, root_rank)

   if (prank == 0) then
      inquire (file=plot_name//"_0000.xmf", exist=file_exists)
      if (file_exists) print *, "PASSED"
   end if

   call sll_s_halt_collective()

end program pic_visu_parallel
