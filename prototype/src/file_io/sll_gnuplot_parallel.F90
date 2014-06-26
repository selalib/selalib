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

!> @brief
!> parallel version of sll_gnuplot
#define MPI_MASTER 0

module sll_gnuplot_parallel

#include "sll_working_precision.h"
#include "sll_assert.h"
use sll_ascii_io, only: sll_ascii_file_create
use sll_utilities, only: sll_new_file_id, int2string
use mpi

implicit none

!> Create a gnuplot file for a 2d mesh (cartesian or curvilinear)
interface sll_gnuplot_2d_parallel
   module procedure sll_gnuplot_curv_2d_parallel
   module procedure sll_gnuplot_rect_2d_parallel
end interface sll_gnuplot_2d_parallel

contains  

!> write a data file plotable by gnuplot to visualize a 2d field
subroutine sll_gnuplot_curv_2d_parallel(array_x, array_y, array, &
                                   array_name, iplot, error)  

sll_real64, dimension(:,:) :: array_x    !< x mesh coordinates
sll_real64, dimension(:,:) :: array_y    !< y mesh coordinates
sll_real64, dimension(:,:) :: array      !< data
character(len=*)           :: array_name !< field name
sll_int32                  :: error      !< error code
sll_int32                  :: iplot      !< plot counter

character(len=4)           :: fin   
sll_int32, save            :: gnu_id
sll_int32                  :: file_id
sll_int32                  :: i, j
character(len=4)           :: cproc 
sll_int32                  :: comm, iproc, nproc
logical                    :: dir_e
logical, save              :: first_call = .true.

call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,error)
call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,error)
comm  = MPI_COMM_WORLD
call int2string(iproc, cproc)
call int2string(iplot, fin)

inquire(file=cproc//"/"".", exist=dir_e)
if (.not. dir_e) then
   call system("mkdir -p "//cproc)
end if

SLL_ASSERT(size(array_x,1) == size(array_y,1))
SLL_ASSERT(size(array_x,2) == size(array_y,2))
SLL_ASSERT(size(array,  1) == size(array_x,1))
SLL_ASSERT(size(array,  2) == size(array_x,2))
SLL_ASSERT(size(array,  1) == size(array_y,1))
SLL_ASSERT(size(array,  2) == size(array_y,2))

call sll_new_file_id(file_id, error)
call sll_ascii_file_create(cproc//"/"//array_name//'_'//fin//'.dat', file_id, error )

! we need to get rid of size() calls for anything that is not an argument
! check...
do i = 1, size(array,1)
   do j = 1, size(array,2)
      write(file_id,*) sngl(array_x(i,j)), sngl(array_y(i,j)), sngl(array(i,j))
   end do
   write(file_id,*)
enddo
close(file_id)

if (iproc == MPI_MASTER) then

   if (first_call) then
      call sll_new_file_id(gnu_id, error)
      open(gnu_id,file=array_name//".gnu")
      rewind(gnu_id)
      first_call = .false.
   else
      open(gnu_id,file=array_name//".gnu", position="append")
   end if

   write(gnu_id,*)"set title 'Time = ",iplot,"'"
   write(gnu_id,"(a)",advance='no') &
   "splot '"//"0000/"//array_name//'_'//fin//".dat' w l"
   do iproc = 1, nproc - 1
      call int2string(iproc, cproc)
      write(gnu_id,"(a)",advance='no')  &
      ",'"//cproc//"/"//array_name//'_'//fin//".dat' w l "
   end do
   write(gnu_id,*)
   close(gnu_id)

end if

end subroutine sll_gnuplot_curv_2d_parallel


! There is a problem with this function. This function should be explicitly 
! told which collective
! to use, thus permitting the writing of data by a subset of the processors
! if needed. This assumes that all the processors will contain the writable
! data. The collective should be passed as an argument. But this may imply that
! different collectives will write data with the same name... further changes
! are needed.

!> write a data file plotable by gnuplot to visualize a 2d field
subroutine sll_gnuplot_rect_2d_parallel(x_min, delta_x, &
                                        y_min, delta_y, &
                                        npts_x, npts_y, &
                                        array, array_name, iplot, error)  

  sll_real64, intent(in)       :: x_min      !< Box corners
  sll_real64, intent(in)       :: delta_x    !< step size
  sll_real64, intent(in)       :: y_min      !< Box corners
  sll_real64, intent(in)       :: delta_y    !< step size
  sll_int32,  intent(in)       :: npts_x     !< number of points to be written, x
  sll_int32,  intent(in)       :: npts_y     !< number of points to be written, y
  sll_real64, intent(in)       :: array(:,:) !< data
  character(len=*), intent(in) :: array_name !< field name
  sll_int32                    :: error      !< error code
  sll_int32                    :: iplot      !< plot counter
  
  character(len=4)             :: fin   
  sll_int32, save              :: gnu_id
  sll_int32                    :: file_id
  sll_int32                    :: i, j
  sll_real64                   :: x, y
  character(len=4)             :: cproc 
  sll_int32                    :: comm, iproc, nproc
  logical                      :: dir_e
  logical, save                :: first_call= .true.
  

  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,error)
  call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,error)
  comm  = MPI_COMM_WORLD
  call int2string(iproc, cproc)
  call int2string(iplot, fin)
  
  inquire(file=cproc//"/"".", exist=dir_e)
  if (.not. dir_e) then
     call system("mkdir -p "//cproc)
  end if
  
  call sll_new_file_id(file_id, error)
  call sll_ascii_file_create(cproc//"/"//array_name//'_'//fin//'.dat', file_id, error )
  do j = 1, npts_y
     y = y_min + (j-1)*delta_y  
     do i = 1, npts_x
        x = x_min+(i-1)*delta_x
        write(file_id,*) sngl(x), sngl(y), sngl(array(i,j))
     end do
     write(file_id,*)
  enddo
  close(file_id)

  if (iproc == MPI_MASTER) then
     
     if (first_call) then
        call sll_new_file_id(gnu_id, error)
        open(gnu_id,file=array_name//".gnu")
        rewind(gnu_id)
        first_call = .false.
     else
        open(gnu_id,file=array_name//".gnu", position="append")
     end if
     
     write(gnu_id,*)"set title 'Time = ",iplot,"'"
     write(gnu_id,"(a)",advance='no') &
          "splot '"//"0000/"//array_name//'_'//fin//".dat' w l"
     do iproc = 1, nproc - 1
        call int2string(iproc, cproc)
        write(gnu_id,"(a)",advance='no')  &
             ",'"//cproc//"/"//array_name//'_'//fin//".dat' w l "
     end do
     write(gnu_id,*)
     close(gnu_id)
     
  end if
  
end subroutine sll_gnuplot_rect_2d_parallel





! it would be nice to have a function like the one commented below, 
! (unfinished) which would work with layouts. But this means that the level
! of abstraction of this module would be elevated above remap. This is very
! tricky since this module is used as a reasonably low-level module. Yet some
! solution to this would be desirable.


!!$!> write a data file plotable by gnuplot to visualize a 2d field
!!$subroutine sll_gnuplot_draw_parallel_array_2d( &
!!$     layout_2d, &
!!$     x_min_global, &
!!$     delta_x, &
!!$     y_min_global, &
!!$     delta_y, &
!!$     array, &
!!$     array_name, &
!!$     iplot, &
!!$     error)  
!!$
!!$  type(layout_2D), pointer               :: layout_2d
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
!!$  col =>get_layout_collective(layout_2D)
!!$  nproc = sll_get_collective_size(col)
!!$  iproc = sll_get_collective_rank(col)
!!$  call int2string(iproc, cproc)
!!$  comm  = sll_world_collective%comm
!!$  call int2string(iplot, fin)
!!$
!!$  gi(1:2) = local_to_global(layout_2D, (/ 1, 1 /) )
!!$
!!$  !shouldn't we allow this same instruction when iplot is zero??
!!$  if (iplot == 1) then
!!$     inquire(file=cproc//"/"".", exist=dir_e)
!!$     if (dir_e) then
!!$        write(*,*) "directory "//cproc//" exists!"
!!$     else
!!$        call system("mkdir -p "//cproc)
!!$     end if
!!$  end if
!!$
!!$  call sll_new_file_id(file_id, error)
!!$  call sll_ascii_file_create(cproc//"/"//array_name//'_'//fin//'.dat', &
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
!!$   if (iplot == 1) call sll_new_file_id(gnu_id, error)
!!$   open(gnu_id,file=array_name//".gnu", position="append")
!!$   if (iplot == 1) rewind(gnu_id)
!!$
!!$   write(gnu_id,*)"set title 'Time = ",iplot,"'"
!!$   write(gnu_id,"(a)",advance='no') &
!!$   "splot '"//"0000/"//array_name//'_'//fin//".dat' w l"
!!$   do iproc = 1, nproc - 1
!!$      call int2string(iproc, cproc)
!!$      write(gnu_id,"(a)",advance='no')  &
!!$      ",'"//cproc//"/"//array_name//'_'//fin//".dat' w l "
!!$   end do
!!$   write(gnu_id,*)
!!$   close(gnu_id)
!!$
!!$end if
!!$
!!$end subroutine sll_gnuplot_rect_2d_parallel


end module sll_gnuplot_parallel
