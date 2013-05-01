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

!> @author Pierre Navaro
!> @brief
!> Implements the functions to write data file plotable by GNUplot
#define MPI_MASTER 0

module sll_gnuplot_parallel

#include "sll_working_precision.h"
#include "sll_assert.h"

use sll_ascii_io, only: sll_ascii_file_create
use sll_utilities, only: sll_new_file_id, int2string
use sll_collective
use sll_remapper

implicit none

contains  

!> write a data file plotable by gnuplot to visualize a 2d field
subroutine sll_gnuplot_rect_2d_parallel(x_min, delta_x, &
                                        y_min, delta_y, &
                                        array, array_name, iplot, error)  

sll_real64                 :: x_min      !< Box corners
sll_real64                 :: delta_x    !< step size
sll_real64                 :: y_min      !< Box corners
sll_real64                 :: delta_y    !< step size
sll_real64, dimension(:,:) :: array      !< data
character(len=*)           :: array_name !< field name
sll_int32                  :: error      !< error code
sll_int32                  :: iplot      !< plot counter

character(len=4)           :: fin   
sll_int32, save            :: gnu_id
sll_int32                  :: file_id
sll_int32                  :: i, j
sll_real64                 :: x, y
character(len=4)           :: cproc 
sll_int32                  :: comm, iproc, nproc
logical                    :: dir_e


nproc = sll_get_collective_size(sll_world_collective)
iproc = sll_get_collective_rank(sll_world_collective)
call int2string(iproc, cproc)
comm  = sll_world_collective%comm
call int2string(iplot, fin)

if (iplot == 1) then
   inquire(file=cproc//"/"".", exist=dir_e)
   if (dir_e) then
     write(*,*) "directory "//cproc//" exists!"
   else
     call system("mkdir -p "//cproc)
   end if
end if


call sll_new_file_id(file_id, error)
call sll_ascii_file_create(cproc//"/"//array_name//'_'//fin//'.dat', file_id, error )
x = x_min
do i = 1, size(array,1)
   y = y_min
   do j = 1, size(array,2)
      write(file_id,*) x, y, sngl(array(i,j))
      y = y + delta_y
   end do
   x = x + delta_x
   write(file_id,*)
enddo
close(file_id)

if (iproc == MPI_MASTER) then

   if (iplot == 1) call sll_new_file_id(gnu_id, error)
   open(gnu_id,file=array_name//".gnu", position="append")
   if (iplot == 1) rewind(gnu_id)

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

end module sll_gnuplot_parallel
