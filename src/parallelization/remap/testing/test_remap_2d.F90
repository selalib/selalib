program remap_2d_unit_test
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

   use iso_fortran_env, only: &
      output_unit

   use sll_m_collective, only: &
      sll_s_boot_collective, &
      sll_s_collective_barrier, &
      sll_o_collective_reduce, &
      sll_f_get_collective_rank, &
      sll_f_get_collective_size, &
      sll_s_halt_collective, &
      sll_v_world_collective

   use sll_m_remapper, only: &
      sll_o_apply_remap_2d, &
      sll_o_compute_local_sizes, &
      sll_o_get_layout_i_max, &
      sll_o_get_layout_i_min, &
      sll_o_get_layout_j_max, &
      sll_o_get_layout_j_min, &
      sll_o_initialize_layout_with_distributed_array, &
      sll_t_layout_2d, &
      sll_o_local_to_global, &
      sll_f_new_layout_2d, &
      sll_o_new_remap_plan, &
      sll_t_remap_plan_2d_comp64, &
      sll_t_remap_plan_2d_real64, &
      sll_o_set_layout_i_max, &
      sll_o_set_layout_i_min, &
      sll_o_set_layout_j_max, &
      sll_o_set_layout_j_min, &
      sll_o_delete, &
      sll_o_get_num_nodes, &
      sll_o_view_lims

   use sll_m_utilities, only: &
      sll_f_is_power_of_two

   use mpi, only: &
      mpi_prod

   implicit none
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   ! Test of the 2D remapper takes a 2D array whose global size Nx*Ny,
   ! distributed among NPi*NPj processors.
   sll_real64, dimension(:, :), allocatable :: local_array1
   sll_real64, dimension(:, :), allocatable :: local_array2
   sll_comp64, dimension(:, :), allocatable :: local_array1c
   sll_comp64, dimension(:, :), allocatable :: local_array2c
   sll_real64, dimension(:, :), allocatable :: arrays_diff
   sll_comp64, dimension(:, :), allocatable :: arrays_diffc
   ! Take a 2D array of dimensions ni*nj where ni, nj are the dimensions of
   ! the full array.
   integer, parameter                       :: ni = 129 !512
   integer, parameter                       :: nj = 512 !512
   ! Local sizes
   integer                                   :: loc_sz_i_init
   integer                                   :: loc_sz_j_init
   integer                                   :: loc_sz_i_final
   integer                                   :: loc_sz_j_final

   ! the process mesh
   integer                                   :: npi
   integer                                   :: npj
   sll_int32                                 :: gi, gj
   integer                                   :: ierr
   integer                                   :: myrank
   sll_int64                                 :: colsz        ! collective size
   sll_real64                                :: val
   ! Remap stuff
   type(sll_t_layout_2d), pointer                  :: layout1
   type(sll_t_layout_2d), pointer                  :: layout2
   type(sll_t_remap_plan_2d_real64), pointer       :: rmp2
   type(sll_t_remap_plan_2d_comp64), pointer       :: rmp2c

   sll_real64                                :: rand_real
   integer, parameter                        :: nbtest = 25
   integer                                   :: i_test
   integer                                   :: i, j
   sll_int32, dimension(2)                   :: global_indices, g
   sll_real32, dimension(1)               :: prod4test
   logical                                   :: test_passed
   integer                                   :: ok
   sll_int32, dimension(2)                   :: tmp_array
   test_passed = .true.

   ! Boot parallel environment
   call sll_s_boot_collective()

   colsz = int(sll_f_get_collective_size(sll_v_world_collective), i64)
   myrank = sll_f_get_collective_rank(sll_v_world_collective)

   if (myrank .eq. 0) then
      print *, ' '
      print *, '--------------- REMAP test ---------------------'
      print *, ' '
      print *, 'Running a test on ', colsz, 'processes'
      flush (output_unit)
   end if

   if (.not. sll_f_is_power_of_two(colsz)) then
      print *, 'This test needs to run in a number of processes which is ', &
         'a power of 2.'
      stop
   end if

   do, i_test = 1, nbtest
   flush (output_unit)
   if (myrank .eq. 0) then
      print *, 'Iteration ', i_test, ' of ', nbtest
   end if
   layout1 => sll_f_new_layout_2d(sll_v_world_collective)
   call factorize_in_random_2powers_2d(colsz, npi, npj)
   if (i_test == 1) then
      npi = 1
      npj = int(colsz, i32)
   end if
   if (i_test == 2) then
      npi = int(colsz, i32)
      npj = 1
   end if

   if (myrank .eq. 0) then
      print *, 'source configuration: ', npi, npj
   end if

   call sll_o_initialize_layout_with_distributed_array( &
      ni, &
      nj, &
      npi, &
      npj, &
      layout1)

   call sll_o_compute_local_sizes(layout1, loc_sz_i_init, loc_sz_j_init)

   SLL_ALLOCATE(local_array1(loc_sz_i_init, loc_sz_j_init), ierr)

   ! initialize the local data
   do j = 1, loc_sz_j_init
      do i = 1, loc_sz_i_init
         tmp_array(:) = (/i, j/)
         global_indices = sll_o_local_to_global(layout1, tmp_array)
         gi = global_indices(1)
         gj = global_indices(2)
         local_array1(i, j) = real(gi + (gj - 1)*ni, f64)
      end do
   end do

   layout2 => sll_f_new_layout_2d(sll_v_world_collective)
   call factorize_in_random_2powers_2d(colsz, npi, npj)
   if (i_test == 1) then
      npi = int(colsz, i32)
      npj = 1
   end if
   if (i_test == 2) then
      npi = 1
      npj = int(colsz, i32)
   end if

   if (myrank .eq. 0) then
      print *, 'target configuration: ', npi, npj
   end if

   call sll_o_initialize_layout_with_distributed_array( &
      ni, &
      nj, &
      npi, &
      npj, &
      layout2)

   call reorganize_randomly_2d(layout2)

   call sll_o_compute_local_sizes( &
      layout2, &
      loc_sz_i_final, &
      loc_sz_j_final)

   SLL_ALLOCATE(local_array2(loc_sz_i_final, loc_sz_j_final), ierr)
   print *, 'will create new remap plan...'
   rmp2 => sll_o_new_remap_plan(layout1, layout2, local_array1)
   print *, 'created plan, proceeding to apply plan...'
   call sll_o_apply_remap_2d(rmp2, local_array1, local_array2)
   print *, 'applied plan'
   SLL_ALLOCATE(arrays_diff(loc_sz_i_final, loc_sz_j_final), ierr)
#if 0
   if (myrank .eq. 0) then
      print *, i_test, myrank, 'Printing layout1: '
      call sll_o_view_lims(layout1)
      print *, i_test, myrank, 'Printing layout2: '
      call sll_o_view_lims(layout2)
   end if
#endif
   ! compare results with expected data
   do j = 1, loc_sz_j_final
      do i = 1, loc_sz_i_final
         tmp_array(:) = (/i, j/)
         global_indices = sll_o_local_to_global(layout2, tmp_array)
         gi = global_indices(1)
         gj = global_indices(2)
         arrays_diff(i, j) = local_array2(i, j) - (gi + (gj - 1)*ni)
         if (arrays_diff(i, j) /= 0) then
            test_passed = .false.
            print *, i_test, myrank, '"remap" unit test: FAIL'
            print *, i_test, myrank, 'local indices: ', '(', i, j, ')'
            print *, 'global indices wrt target layout'
            print *, myrank, 'in global indices: (', gi, ',', gj, ')'
            print *, myrank, 'local array1(', gi, ',', gj, ')=', local_array1(i, j)
            print *, myrank, 'local array2(', gi, ',', gj, ')=', local_array2(i, j)
            g = theoretical_global_2D_indices(int(local_array2(i, j), i32), ni)
            print *, 'Theoretical indices: (', g(1), ',', g(2), ')'
            !     if(myrank .eq. 1) then
            print *, local_array2(:, :)
            !    end if
            print *, i_test, myrank, 'Printing layout1: '
            call sll_o_view_lims(layout1)
            print *, i_test, myrank, 'Printing layout2: '
            call sll_o_view_lims(layout2)

            print *, 'program stopped by failure'
            stop
         end if
         flush (output_unit)
      end do
   end do

   ! As the difference described above is null in each point, the sum of the
   ! corresponding absolute values must be null. Each processor compute a
   ! local sum and all local sums are finally added and the result is sent to
   ! processor 0 which will check if equal 0 to validate the test. (*)
   ok = 1
   call sll_o_collective_reduce( &
      sll_v_world_collective, &
      (/real(ok)/), &
      1, &
      MPI_PROD, &
      0, &
      prod4test)

   if (myrank .eq. 0) then
      print *, ' '
      print *, '-------------------------------------------'
      print *, ' '
      flush (output_unit)
   end if
   flush (output_unit)

   call sll_s_collective_barrier(sll_v_world_collective)

   call sll_o_delete(layout1)
   call sll_o_delete(layout2)
   SLL_DEALLOCATE_ARRAY(local_array1, ierr)
   SLL_DEALLOCATE_ARRAY(local_array2, ierr)
   SLL_DEALLOCATE_ARRAY(arrays_diff, ierr)
   end do

#if 0
   if (myrank == 0) then
      if (prod4test(1) == 1.) then
         print *, 'PASSED'
      end if
   end if
#endif

   !**************************************************************************
   ! Test complex data type case
   !
   !**************************************************************************

   if (myrank .eq. 0) then
      print *, ' '
      print *, '--------------- REMAP 2D test: complex case ------------------'
      print *, ' '
      print *, 'Running a test on ', colsz, 'processes'
      flush (output_unit)
   end if

   do, i_test = 1, nbtest
   flush (output_unit)
   if (myrank .eq. 0) then
      print *, 'Iteration ', i_test, ' of ', nbtest
   end if
   layout1 => sll_f_new_layout_2d(sll_v_world_collective)
   call factorize_in_random_2powers_2d(colsz, npi, npj)

   if (myrank .eq. 0) then
      print *, 'source configuration: ', npi, npj
   end if

   call sll_o_initialize_layout_with_distributed_array( &
      ni, &
      nj, &
      npi, &
      npj, &
      layout1)

   call sll_o_compute_local_sizes(layout1, loc_sz_i_init, loc_sz_j_init)

   SLL_ALLOCATE(local_array1c(loc_sz_i_init, loc_sz_j_init), ierr)

   ! initialize the local data
   do j = 1, loc_sz_j_init
      do i = 1, loc_sz_i_init
         tmp_array(:) = (/i, j/)
         global_indices = sll_o_local_to_global(layout1, tmp_array)
         gi = global_indices(1)
         gj = global_indices(2)
         val = real(gi + (gj - 1)*ni, f64)
         local_array1c(i, j) = cmplx(val, val, f64)
      end do
   end do

   layout2 => sll_f_new_layout_2d(sll_v_world_collective)
   call factorize_in_random_2powers_2d(colsz, npi, npj)
   if (myrank .eq. 0) then
      print *, 'target configuration: ', npi, npj
   end if

   call sll_o_initialize_layout_with_distributed_array( &
      ni, &
      nj, &
      npi, &
      npj, &
      layout2)

   call reorganize_randomly_2d(layout2)

   call sll_o_compute_local_sizes( &
      layout2, &
      loc_sz_i_final, &
      loc_sz_j_final)

   SLL_ALLOCATE(local_array2c(loc_sz_i_final, loc_sz_j_final), ierr)

   rmp2c => sll_o_new_remap_plan(layout1, layout2, local_array1c)

   call sll_o_apply_remap_2d(rmp2c, local_array1c, local_array2c)

   SLL_ALLOCATE(arrays_diffc(loc_sz_i_final, loc_sz_j_final), ierr)
#if 0
   if (myrank .eq. 0) then
      print *, i_test, myrank, 'Printing layout1: '
      call sll_o_view_lims(layout1)
      print *, i_test, myrank, 'Printing layout2: '
      call sll_o_view_lims(layout2)
   end if
#endif
   ! compare results with expected data
   do j = 1, loc_sz_j_final
      do i = 1, loc_sz_i_final
         tmp_array(:) = (/i, j/)
         global_indices = sll_o_local_to_global(layout2, tmp_array)
         gi = global_indices(1)
         gj = global_indices(2)
         val = real(gi + (gj - 1)*ni, f64)
         arrays_diffc(i, j) = local_array2c(i, j) - cmplx(val, val, f64)
         if (arrays_diffc(i, j) /= (0., 0.)) then
            test_passed = .false.
            print *, i_test, myrank, '"remap" unit test: FAIL'
            print *, i_test, myrank, 'local indices: ', '(', i, j, ')'
            print *, 'global indices wrt target layout'
            print *, myrank, 'in global indices: (', gi, ',', gj, ')'
            print *, myrank, 'local array1c(', gi, ',', gj, ')=', local_array1c(i, j)
            print *, myrank, 'local array2c(', gi, ',', gj, ')=', local_array2c(i, j)
            g = theoretical_global_2D_indices(int(local_array2c(i, j), i32), ni)
            print *, 'Theoretical indices: (', g(1), ',', g(2), ')'
            !     if(myrank .eq. 1) then
            print *, local_array2c(:, :)
            !    end if
            print *, i_test, myrank, 'Printing layout1: '
            call sll_o_view_lims(layout1)
            print *, i_test, myrank, 'Printing layout2: '
            call sll_o_view_lims(layout2)

            print *, 'program stopped by failure'
            stop
         end if
         flush (output_unit)
      end do
   end do

   ! As the difference described above is null in each point, the sum of the
   ! corresponding absolute values must be null. Each processor compute a
   ! local sum and all local sums are finally added and the result is sent to
   ! processor 0 which will check if equal 0 to validate the test. (*)
   ok = 1
   call sll_o_collective_reduce( &
      sll_v_world_collective, &
      (/real(ok)/), &
      1, &
      MPI_PROD, &
      0, &
      prod4test)

   if (myrank .eq. 0) then
      print *, ' '
      print *, '-------------------------------------------'
      print *, ' '
      flush (output_unit)
   end if
   flush (output_unit)

   call sll_s_collective_barrier(sll_v_world_collective)

   call sll_o_delete(layout1)
   call sll_o_delete(layout2)

   SLL_DEALLOCATE_ARRAY(local_array1c, ierr)
   SLL_DEALLOCATE_ARRAY(local_array2c, ierr)
   SLL_DEALLOCATE_ARRAY(arrays_diffc, ierr)

end do

if (myrank == 0) then
   if (prod4test(1) == 1.) then
      print *, 'PASSED'
   end if
end if

call sll_s_halt_collective()

contains

subroutine reorganize_randomly_2d(layout)
   implicit none
   type(sll_t_layout_2d), pointer   :: layout
   integer                    :: i, colsz, proc_n, proc_p
   real                       :: rand_real
   colsz = sll_o_get_num_nodes(layout)
   do i = 0, colsz - 1
      call random_number(rand_real)
      proc_n = int(rand_real*(colsz - 1))
      call random_number(rand_real)
      proc_p = int(rand_real*(colsz - 1))
      call swap_box_2D(proc_n, proc_p, layout)
   end do
end subroutine reorganize_randomly_2d

subroutine swap_box_2D(proc_n, proc_p, layout)
   implicit none
   integer                  :: proc_n, proc_p
   type(sll_t_layout_2d), pointer :: layout
   integer                  :: i_min_n, i_max_n, j_min_n, j_max_n
   integer                  :: i_min_p, i_max_p, j_min_p, j_max_p
   ! Get proc_n contents from layout
   i_min_n = sll_o_get_layout_i_min(layout, proc_n)
   i_max_n = sll_o_get_layout_i_max(layout, proc_n)
   j_min_n = sll_o_get_layout_j_min(layout, proc_n)
   j_max_n = sll_o_get_layout_j_max(layout, proc_n)
   ! Get proc_p contents from layout
   i_min_p = sll_o_get_layout_i_min(layout, proc_p)
   i_max_p = sll_o_get_layout_i_max(layout, proc_p)
   j_min_p = sll_o_get_layout_j_min(layout, proc_p)
   j_max_p = sll_o_get_layout_j_max(layout, proc_p)
   ! Set proc_n contents in layout
   call sll_o_set_layout_i_min(layout, proc_n, i_min_p)
   call sll_o_set_layout_i_max(layout, proc_n, i_max_p)
   call sll_o_set_layout_j_min(layout, proc_n, j_min_p)
   call sll_o_set_layout_j_max(layout, proc_n, j_max_p)
   ! Set proc_p contents in layout
   call sll_o_set_layout_i_min(layout, proc_p, i_min_n)
   call sll_o_set_layout_i_max(layout, proc_p, i_max_n)
   call sll_o_set_layout_j_min(layout, proc_p, j_min_n)
   call sll_o_set_layout_j_max(layout, proc_p, j_max_n)
end subroutine swap_box_2D

! d: index of the node, in linear order.
! ni: global size of array, first dimension
function theoretical_global_2D_indices(d, ni)
   integer, dimension(1:2) :: theoretical_global_2D_indices
   integer, intent(in)     :: d, ni
   integer                 :: i, j
   integer                 :: q
   sll_real64              :: tmp
   ! We know that the linear index relates with the other
   ! indices by:
   !
   !          linear = i + ni*(j-1)
   !
   ! where i,j are the indices sought. We start by working backwards.
   tmp = real(d, f64)/real(ni, f64)
   j = ceiling(tmp)
   ! reduce the size of the number we are working with
   q = d - (j - 1)*ni
   i = q
   theoretical_global_2D_indices(1:2) = (/i, j/)
end function theoretical_global_2D_indices

subroutine factorize_in_random_2powers_2d(n, n1, n2)
   sll_int64, intent(in) :: n
   integer, intent(out)  ::n1, n2
   integer   :: expo, expo1, expo2
   if (.not. sll_f_is_power_of_two(n)) then
      print *, 'The number of processors must be a power of 2'
      stop
   end if
   expo = int(log(real(n))/log(2.))
   call random_number(rand_real)
   expo1 = int(rand_real*expo)
   expo2 = expo - expo1
   n1 = 2**expo1
   n2 = 2**expo2
end subroutine factorize_in_random_2powers_2d

end program remap_2d_unit_test
