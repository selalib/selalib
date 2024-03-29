!**************************************************************
!  Copyright INRIA
!  Authors :
!     Pierre Navaro
!
!  This code SeLaLib (for Semi-Lagrangian-Library)
!  is a parallel library for simulating the plasma turbulence
!  in a tokamak.
!
!  This software is governed by the CeCILL-B license
!  under French law and abiding by the rules of distribution
!  of free software.  You can  use, modify and redistribute
!  the software under the terms of the CeCILL-B license as
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info".
!**************************************************************

!> @ingroup file_io_parallel
!> @brief
!> parallel version of sll_m_gnuplot
!> @details
!> Data files are in ASCII format so these subroutines are slow and use
!> a lot of disk space. Consider to use it for debug purpose.
!> Here an example using the sll_t_layout_2d object, check out how to compute
!> offset values before calling sll_m_gnuplot_parallel subroutines.
!> @snippet remap/unit_test_parallel.F90 example
#define MPI_MASTER 0

module sll_m_gnuplot_parallel

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_assert.h"
#include "sll_working_precision.h"

#ifdef __INTEL_COMPILER
   use ifport
#endif

   use sll_m_ascii_io, only: &
      sll_s_ascii_file_create

   use sll_m_utilities, only: &
      sll_s_int2string, &
      sll_s_new_file_id

   use sll_m_collective, only: &
      sll_f_get_collective_rank, &
      sll_f_get_collective_size, &
      sll_v_world_collective

   implicit none

   public :: &
      sll_o_gnuplot_2d_parallel, &
      sll_s_gnuplot_curv_2d_parallel, &
      sll_s_gnuplot_rect_2d_parallel

   private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!> Create a gnuplot file for a 2d mesh (cartesian or curvilinear)
   interface sll_o_gnuplot_2d_parallel
      module procedure sll_s_gnuplot_curv_2d_parallel
      module procedure sll_s_gnuplot_rect_2d_parallel
   end interface sll_o_gnuplot_2d_parallel

contains

!> write a data file plotable by gnuplot to visualize a 2d field
   subroutine sll_s_gnuplot_curv_2d_parallel(array_x, array_y, array, &
                                             array_name, iplot, error)

      sll_real64, dimension(:, :), intent(in) :: array_x    !< x mesh coordinates
      sll_real64, dimension(:, :), intent(in) :: array_y    !< y mesh coordinates
      sll_real64, dimension(:, :), intent(in) :: array      !< data
      character(len=*), intent(in) :: array_name !< field name
      sll_int32, intent(in) :: iplot      !< plot counter
      sll_int32, intent(out):: error      !< error code

      character(len=4)           :: fin
      sll_int32, save            :: gnu_id
      sll_int32                  :: file_id
      sll_int32                  :: i, j
      character(len=4)           :: cproc
      sll_int32                  :: iproc, nproc
      logical                    :: dir_e
      logical, save              :: first_call = .true.

      iproc = sll_f_get_collective_rank(sll_v_world_collective)
      nproc = sll_f_get_collective_size(sll_v_world_collective)

      call sll_s_int2string(iproc, cproc)
      call sll_s_int2string(iplot, fin)

#ifdef __INTEL_COMPILER
      if (makedirqq(cproc)) then
         print *, ' Make directory '//cproc
      end if
#else
      inquire (file=cproc//"/"".", exist=dir_e)
      if (.not. dir_e) then
         call execute_command_line("mkdir -p "//cproc)
      end if
#endif

      SLL_ASSERT_ALWAYS(size(array_x, 1) == size(array_y, 1))
      SLL_ASSERT_ALWAYS(size(array_x, 2) == size(array_y, 2))
      SLL_ASSERT_ALWAYS(size(array, 1) == size(array_x, 1))
      SLL_ASSERT_ALWAYS(size(array, 2) == size(array_x, 2))
      SLL_ASSERT_ALWAYS(size(array, 1) == size(array_y, 1))
      SLL_ASSERT_ALWAYS(size(array, 2) == size(array_y, 2))

      call sll_s_new_file_id(file_id, error)
      call sll_s_ascii_file_create(cproc//"/"//array_name//'_'//fin//'.dat', file_id, error)

      ! we need to get rid of size() calls for anything that is not an argument
      ! check...
      do i = 1, size(array, 1)
         do j = 1, size(array, 2)
            write (file_id, *) sngl(array_x(i, j)), sngl(array_y(i, j)), sngl(array(i, j))
         end do
         write (file_id, *)
      end do
      close (file_id)

      if (iproc == MPI_MASTER) then

         if (first_call) then
            call sll_s_new_file_id(gnu_id, error)
            open (gnu_id, file=array_name//".gnu")
            rewind (gnu_id)
            first_call = .false.
         else
            open (gnu_id, file=array_name//".gnu", position="append")
         end if

         write (gnu_id, *) "set title 'Time = ", iplot, "'"
         write (gnu_id, "(a)", advance='no') &
            "splot '"//"0000/"//array_name//'_'//fin//".dat' w l"
         do iproc = 1, nproc - 1
            call sll_s_int2string(iproc, cproc)
            write (gnu_id, "(a)", advance='no') &
               ",'"//cproc//"/"//array_name//'_'//fin//".dat' w l "
         end do
         write (gnu_id, *)
         close (gnu_id)

      end if

   end subroutine sll_s_gnuplot_curv_2d_parallel

!> @brief
!> write a data file plotable by gnuplot to visualize a 2d field
!> @details
!> Data are written in parallel, every processor writes in a different directory.
!> @todo
!! There is a problem with this function. This function should be explicitly
!! told which collective
!! to use, thus permitting the writing of data by a subset of the processors
!! if needed. This assumes that all the processors will contain the writable
!! data. The collective should be passed as an argument. But this may imply that
!! different collectives will write data with the same name... further changes
!! are needed.
   subroutine sll_s_gnuplot_rect_2d_parallel(x_min, delta_x, &
                                             y_min, delta_y, &
                                             npts_x, npts_y, &
                                             array, array_name, &
                                             iplot, error)

      sll_real64, intent(in)       :: x_min      !< Box corners
      sll_real64, intent(in)       :: delta_x    !< step size
      sll_real64, intent(in)       :: y_min      !< Box corners
      sll_real64, intent(in)       :: delta_y    !< step size
      sll_int32, intent(in)       :: npts_x     !< number of points to be written, x
      sll_int32, intent(in)       :: npts_y     !< number of points to be written, y
      sll_real64, intent(in)       :: array(:, :) !< data
      character(len=*), intent(in) :: array_name !< field name
      sll_int32, intent(in)       :: iplot      !< plot counter
      sll_int32, intent(out)      :: error      !< error code

      character(len=4)             :: fin
      sll_int32, save              :: gnu_id
      sll_int32                    :: file_id
      sll_int32                    :: i, j
      sll_real64                   :: x, y
      character(len=4)             :: cproc
      sll_int32                    :: iproc, nproc
      logical                      :: dir_e
      logical, save                :: first_call = .true.

      iproc = sll_f_get_collective_rank(sll_v_world_collective)
      nproc = sll_f_get_collective_size(sll_v_world_collective)

      call sll_s_int2string(iproc, cproc)
      call sll_s_int2string(iplot, fin)

#ifdef __INTEL_COMPILER
      if (makedirqq(cproc)) then
         print *, ' Make directory '//cproc
      end if
#else
      inquire (file=cproc//"/"".", exist=dir_e)
      if (.not. dir_e) then
         call execute_command_line("mkdir -p "//cproc)
      end if
#endif

      call sll_s_new_file_id(file_id, error)
      call sll_s_ascii_file_create(cproc//"/"//array_name//'_'//fin//'.dat', file_id, error)
      do j = 1, npts_y
         y = y_min + real(j - 1, f64)*delta_y
         do i = 1, npts_x
            x = x_min + real(i - 1, f64)*delta_x
            write (file_id, *) sngl(x), sngl(y), sngl(array(i, j))
         end do
         write (file_id, *)
      end do
      close (file_id)

      if (iproc == MPI_MASTER) then

         if (first_call) then
            call sll_s_new_file_id(gnu_id, error)
            open (gnu_id, file=array_name//".gnu")
            rewind (gnu_id)
            first_call = .false.
         else
            open (gnu_id, file=array_name//".gnu", position="append")
         end if

         write (gnu_id, *) "set title 'Time = ", iplot, "'"
         write (gnu_id, *) "set output '"//array_name//'_'//fin//".png'"
         write (gnu_id, "(a)", advance='no') &
            "splot '"//"0000/"//array_name//'_'//fin//".dat' w l"
         do iproc = 1, nproc - 1
            call sll_s_int2string(iproc, cproc)
            write (gnu_id, "(a)", advance='no') &
               ",'"//cproc//"/"//array_name//'_'//fin//".dat' w l "
         end do
         write (gnu_id, *)
         close (gnu_id)

      end if

   end subroutine sll_s_gnuplot_rect_2d_parallel

! it would be nice to have a function like the one commented below,
! (unfinished) which would work with layouts. But this means that the level
! of abstraction of this module would be elevated above remap. This is very
! tricky since this module is used as a reasonably low-level module. Yet some
! solution to this would be desirable.

!!$!> write a data file plotable by gnuplot to visualize a 2d field
!!$subroutine sll_gnuplot_draw_parallel_array_2d( &
!!$     sll_t_layout_2d, &
!!$     x_min_global, &
!!$     delta_x, &
!!$     y_min_global, &
!!$     delta_y, &
!!$     array, &
!!$     array_name, &
!!$     iplot, &
!!$     error)
!!$
!!$  type(sll_t_layout_2d), pointer               :: sll_t_layout_2d
!!$  sll_real64, intent(in)                 :: x_min_global !< Box corners
!!$  sll_real64, intent(in)                 :: delta_x      !< step size
!!$  sll_real64, intent(in)                 :: y_min_global !< Box corners
!!$  sll_real64, intent(in)                 :: delta_y      !< step size
!!$  sll_real64, dimension(:,:), intent(in) :: array        !< data
!!$  character(len=*), intent(in)           :: array_name   !< field name
!!$  sll_int32, intent(in)                  :: iplot        !< plot counter
!!$  sll_int32, intent(in)                  :: error        !< error code
!!$
!!$  type(sll_collective_t), pointer :: col
!!$  sll_int32, dimension(2)         :: gi  ! for storing global indices of corner
!!$  character(len=4)                :: fin
!!$  sll_int32, save                 :: gnu_id
!!$  sll_int32                       :: file_id
!!$  sll_int32                       :: i, j
!!$  sll_real64                      :: x, y
!!$  character(len=4)                :: cproc
!!$  sll_int32                       :: comm, iproc, nproc
!!$  logical                         :: dir_e
!!$
!!$  col =>get_layout_collective(sll_t_layout_2d)
!!$  nproc = sll_f_get_collective_size(col)
!!$  iproc = sll_f_get_collective_rank(col)
!!$  call sll_s_int2string(iproc, cproc)
!!$  comm  = sll_v_world_collective%comm
!!$  call sll_s_int2string(iplot, fin)
!!$
!!$  gi(1:2) = sll_o_local_to_global(sll_t_layout_2d, (/ 1, 1 /) )
!!$
!!$  !shouldn't we allow this same instruction when iplot is zero??
!!$  if (iplot == 1) then
!!$     inquire(file=cproc//"/"".", exist=dir_e)
!!$     if (dir_e) then
!!$        write(*,*) "directory "//cproc//" exists!"
!!$     else
!!$        call execute_command_line("mkdir -p "//cproc)
!!$     end if
!!$  end if
!!$
!!$  call sll_s_new_file_id(file_id, error)
!!$  call sll_s_ascii_file_create(cproc//"/"//array_name//'_'//fin//'.dat', &
!!$       file_id, error )
!!$  do i = 1, size(array,1)
!!$     x = x_min+(i-1)*delta_x
!!$   do j = 1, size(array,2)
!!$      y = y_min + (j-1)*delta_y
!!$      write(file_id,*) x, y, sngl(array(i,j))
!!$   end do
!!$   write(file_id,*)
!!$enddo
!!$close(file_id)
!!$
!!$if (iproc == MPI_MASTER) then
!!$
!!$   if (iplot == 1) call sll_s_new_file_id(gnu_id, error)
!!$   open(gnu_id,file=array_name//".gnu", position="append")
!!$   if (iplot == 1) rewind(gnu_id)
!!$
!!$   write(gnu_id,*)"set title 'Time = ",iplot,"'"
!!$   write(gnu_id,"(a)",advance='no') &
!!$   "splot '"//"0000/"//array_name//'_'//fin//".dat' w l"
!!$   do iproc = 1, nproc - 1
!!$      call sll_s_int2string(iproc, cproc)
!!$      write(gnu_id,"(a)",advance='no')  &
!!$      ",'"//cproc//"/"//array_name//'_'//fin//".dat' w l "
!!$   end do
!!$   write(gnu_id,*)
!!$   close(gnu_id)
!!$
!!$end if
!!$
!!$end subroutine sll_s_gnuplot_rect_2d_parallel

end module sll_m_gnuplot_parallel
