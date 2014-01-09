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
module sll_gnuplot
#include "sll_working_precision.h"
#include "sll_assert.h"

use sll_ascii_io
use sll_utilities, only: sll_new_file_id, int2string

implicit none

!> write file pltable by gnuplot to visualize 2d field
interface sll_gnuplot_field_2d
module procedure sll_gnuplot_corect_2d
module procedure sll_gnuplot_rect_2d
module procedure sll_gnuplot_curv_2d
end interface

contains  

!> write an array
subroutine sll_gnuplot_write(array,array_name,error)

   sll_real64, dimension(:), intent(in) :: array      !< data
   character(len=*), intent(in)         :: array_name !< field name
   sll_int32, intent(out)               :: error      !< error code
   sll_int32                            :: file_id    !< file unit number
   sll_int32                            :: npoints
   sll_int32                            :: ipoints    
   logical                              :: lopen
   
   npoints = size(array)

   error=0

   do 100 file_id=20,99

      inquire(unit=file_id,opened=lopen)
      if(lopen) then
         cycle
      else
         open(file_id,status='SCRATCH',err=100)
         close(file_id,status='DELETE',err=100)
         goto 200
      end if
 
   100 continue
   error=1
   200 continue
   error=0

   open(file_id,FILE=trim(array_name)//".gnu",FORM='FORMATTED',IOSTAT=error)
   rewind(file_id)

   write(file_id,"(a)")"plot '-' t '"//trim(array_name)//"' with linesp"
   do ipoints = 1, npoints
      write(file_id,"(g15.3)") array(ipoints)
   end do
   write(file_id,"(a)")"e"
   write(file_id,"(a)")"pause -1 'Hit return to quit'"
   close(file_id)


   open(file_id,FILE=trim(array_name)//".dat",FORM='FORMATTED',IOSTAT=error)
   rewind(file_id)

   !write(file_id,"(a)")"plot '-' t '"//trim(array_name)//"' with linesp"
   do ipoints = 1, npoints
      write(file_id,*) array(ipoints)
   end do
   !write(file_id,"(a)")"e"
   !write(file_id,"(a)")"pause -1 'Hit return to quit'"
   close(file_id)


end subroutine sll_gnuplot_write


subroutine sll_gnuplot_write_1d( &
  array, &
  x_array, &
  array_name, &
  iplot)

   sll_real64, dimension(:), intent(in) :: array      !< data
   sll_real64, dimension(:), intent(in) :: x_array      !< data
   character(len=*), intent(in)         :: array_name !< field name
   sll_int32,intent(in),optional :: iplot
   sll_int32               :: error      !< error code
   sll_int32                            :: file_id    !< file unit number
   sll_int32                            :: npoints
   sll_int32                            :: ipoints    
   logical                              :: lopen
   character(len=4)      :: cplot
   
   npoints = size(array)

   error=0

   do 100 file_id=20,99

      inquire(unit=file_id,opened=lopen)
      if(lopen) then
         cycle
      else
         open(file_id,status='SCRATCH',err=100)
         close(file_id,status='DELETE',err=100)
         goto 200
      end if
 
   100 continue
   error=1
   200 continue
   error=0
   if(present(iplot))then
     call int2string(iplot,cplot)
     open(file_id,FILE=trim(array_name)//cplot//".dat",FORM='FORMATTED',IOSTAT=error)
   else
     open(file_id,FILE=trim(array_name)//".dat",FORM='FORMATTED',IOSTAT=error)
   endif
   rewind(file_id)

   !write(file_id,"(a)")"plot '-' t '"//trim(array_name)//"' with linesp"
   do ipoints = 1, npoints
      write(file_id,*) x_array(ipoints), array(ipoints)
   end do
   !write(file_id,"(a)")"e"
   !write(file_id,"(a)")"pause -1 'Hit return to quit'"
   close(file_id)


end subroutine sll_gnuplot_write_1d


!> write a data file plotable by gnuplot to visualize a 2d field.
!> Axis are rectangular and spacing is constant
subroutine sll_gnuplot_corect_2d(xmin, xmax, nx, ymin, ymax, ny, &
                               array, array_name, iplot, error)  
sll_real64                   :: xmin       !< Box corners
sll_real64                   :: xmax       !< Box corners
sll_real64                   :: ymin       !< Box corners
sll_real64                   :: ymax       !< Box corners
sll_int32                    :: nx         !< x points number
sll_int32                    :: ny         !< y points number
sll_real64, dimension(nx,ny) :: array      !< data
character(len=*)             :: array_name !< field name
character(len=4)             :: fin   
sll_int32                    :: iplot      !< plot counter
sll_int32, save :: gnu_id
sll_int32 :: file_id
sll_int32 :: error                         !< error code
sll_int32 :: i, j
sll_real64 :: dx, dy, x, y

call int2string(iplot, fin)

if ( iplot == 1 ) then
   call sll_new_file_id(gnu_id, error)
end if

open(gnu_id,file=array_name//".gnu", position="append")
write(gnu_id,*)"splot '"//array_name//"_"//fin//".dat' w l"
close(gnu_id)

call sll_ascii_file_create(array_name//'_'//fin//'.dat', file_id, error )
dx = (xmax-xmin)/(nx-1)
dy = (ymax-ymin)/(ny-1)
x = xmin
do i=1,nx
   y = ymin
   do j=1,ny
      write(file_id,*) sngl(x),sngl(y),sngl(array(i,j))
      y = y + dy
   end do
   x = x + dx
   write(file_id,*)
enddo
close(file_id)

end subroutine sll_gnuplot_corect_2d

!> write a data file plotable by gnuplot to visualize a 2d field on structured
!> rectangular mesh where spacing is not constant
subroutine sll_gnuplot_rect_2d( nx, xvec, ny, yvec,&
                               array, array_name, iplot, error)  
sll_int32                    :: nx         !< x points number
sll_int32                    :: ny         !< y points number
sll_real64, dimension(nx)    :: xvec       !< x coordinates
sll_real64, dimension(ny)    :: yvec       !< y coordiantes
sll_real64, dimension(nx,ny) :: array      !< data
character(len=*)             :: array_name !< field name
character(len=4)             :: fin   
sll_int32                    :: iplot      !< plot counter
sll_int32, save :: gnu_id
sll_int32 :: file_id
sll_int32 :: error                         !< error code
sll_int32 :: i, j

call int2string(iplot, fin)

if ( iplot == 1 ) then
   call sll_new_file_id(gnu_id, error)
end if

open(gnu_id,file=array_name//".gnu", position="append")
write(gnu_id,*)"splot '"//array_name//"_"//fin//".dat' w l"
close(gnu_id)

call sll_ascii_file_create(array_name//'_'//fin//'.dat', file_id, error )
do i=1,nx
   do j=1,ny
      write(file_id,*) sngl(xvec(i)), &
                       sngl(yvec(j)), &
                       sngl(array(i,j)) 
   end do
   write(file_id,*)
enddo
close(file_id)

end subroutine sll_gnuplot_rect_2d


!> write a data file plotable by gnuplot to visualize a 2d field on structured
!> curvilinear mesh
subroutine sll_gnuplot_curv_2d( nx, ny, xcoord, ycoord,&
                               array, array_name, iplot, error)  
sll_int32                    :: nx         !< x points number
sll_int32                    :: ny         !< y points number
sll_real64, dimension(nx,ny) :: xcoord     !< x coordinates
sll_real64, dimension(nx,ny) :: ycoord     !< y coordiantes
sll_real64, dimension(nx,ny) :: array      !< data
character(len=*)             :: array_name !< field name
character(len=4)             :: fin   
sll_int32                    :: iplot      !< plot counter
sll_int32, save :: gnu_id
sll_int32 :: file_id
sll_int32 :: error                         !< error code
sll_int32 :: i, j

call int2string(iplot, fin)

if ( iplot == 1 ) then
   call sll_new_file_id(gnu_id, error)
end if

open(gnu_id,file=array_name//".gnu", position="append")
write(gnu_id,*)"splot '"//array_name//"_"//fin//".dat' w l"
close(gnu_id)

call sll_ascii_file_create(array_name//'_'//fin//'.dat', file_id, error )
do i=1,nx
   do j=1,ny
      write(file_id,*) sngl(xcoord(i,j)), &
                       sngl(ycoord(i,j)), &
                       sngl(array(i,j)) 
   end do
   write(file_id,*)
enddo
close(file_id)

end subroutine sll_gnuplot_curv_2d


end module sll_gnuplot
