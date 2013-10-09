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
module sll_plotmtv
#include "sll_working_precision.h"
#include "sll_assert.h"

use sll_ascii_io
use sll_utilities, only: sll_new_file_id, int2string

implicit none

!> Create the mtv file to plot a structured mesh (cartesian or curvilinear)
interface sll_plotmtv_write
module procedure  sll_plotmtv_curv_2d
end interface sll_plotmtv_write

contains  

!> write a data file plotable by plotmtv to visualize a 2d 
!> curvilinear mesh
subroutine sll_plotmtv_curv_2d( nx, ny, xcoord, ycoord, label, error)  
sll_int32                    :: nx     !< x points number
sll_int32                    :: ny     !< y points number
sll_real64, dimension(nx,ny) :: xcoord !< x coordinates
sll_real64, dimension(nx,ny) :: ycoord !< y coordiantes
character(len=*)             :: label  !< field name
sll_int32                    :: file_id
sll_int32                    :: error                         !< error code
sll_int32                    :: i, j, k
sll_real64                   :: x1, y1

call sll_new_file_id(file_id, error)

call sll_ascii_file_create(label//'.mtv', file_id, error )
write(file_id,"(a)")"$DATA=CURVE2D"
write(file_id,"('% xmin=',f7.3,' xmax=', f7.3)") minval(xcoord), maxval(xcoord)
write(file_id,"('% ymin=',f7.3,' ymax=', f7.3)") minval(ycoord), maxval(ycoord)
write(file_id,"(a)")"% equalscale=T"
write(file_id,"(a)")"% spline=1"
write(file_id,"(a)")"% markertype=2"
write(file_id,"(a)")"% pointID=F"
write(file_id,"(a)")"% toplabel='"//label//"' "
   
do i=1,nx
   do j=1,ny
      write(file_id,*) xcoord(i,j), ycoord(i,j)
   end do
   write(file_id,*)
end do
do j=1,ny
   do i=1,nx
      write(file_id,*) xcoord(i,j), ycoord(i,j)
   end do
   write(file_id,*)
end do

!Numeros des elements
k = 0
do i=1,nx-1
   do j=1,ny-1
      k = k+1
      x1 = 0.25*(xcoord(i,j)+xcoord(i+1,j)+xcoord(i,j+1)+xcoord(i+1,j+1))
      y1 = 0.25*(ycoord(i,j)+ycoord(i+1,j)+ycoord(i,j+1)+ycoord(i+1,j+1))
      write(file_id,"(a)"   ,  advance="no")"@text x1="
      write(file_id,"(g15.3)", advance="no") x1
      write(file_id,"(a)"   ,  advance="no")" y1="
      write(file_id,"(g15.3)", advance="no") y1
      write(file_id,"(a)"   ,  advance="no")" z1=0. lc=4 ll='"
      write(file_id,"(i4)"  ,  advance="no") k
      write(file_id,"(a)")"'"
   end do
end do

write(file_id,*)"$END"
close(file_id)

end subroutine sll_plotmtv_curv_2d


end module sll_plotmtv
