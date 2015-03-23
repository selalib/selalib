!**************************************************************
!  Copyright INRIA
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

!> @ingroup file_io
!> @brief
!> Implements the functions to write data file plotable by GNUplot
!> @details
!> This is an example of how use the sll_gnuplot module.
!> More details about this example
!> @snippet file_io/unit_test_gnuplot.F90 example
!> Here a snapshot when you execute:
!> <code>$ gnuplot -persistent plot_2.gnu </code>
!
!> @image html gnuplot.png
module sll_gnuplot

#include "sll_working_precision.h"
#include "sll_assert.h"

use sll_ascii_io
use sll_utilities, only: sll_new_file_id, int2string

implicit none

private

!> write file plotable by gnuplot to visualize 2d field
interface sll_gnuplot_1d
   module procedure  sll_gnuplot_write_1d
   module procedure  sll_gnuplot_write
   module procedure  sll_gnuplot_write_two_arrays_1d
end interface

!> Write file for gnuplot to display 2d field.
interface sll_gnuplot_2d
   module procedure sll_gnuplot_corect_2d
   module procedure sll_gnuplot_rect_2d
   module procedure sll_gnuplot_curv_2d
   module procedure sll_gnuplot_mesh_2d
   module procedure write_unstructured_field
end interface

public sll_gnuplot_1d, sll_gnuplot_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> Write an array to display with gnuplot
subroutine sll_gnuplot_write(array, array_name, iplot)

  sll_real64, dimension(:), intent(in) :: array      !< data
  character(len=*),         intent(in) :: array_name !< field name
  sll_int32                            :: iplot      !< plot counter
  sll_int32                            :: file_id    !< file unit number
  sll_int32                            :: npoints
  sll_int32                            :: ipoints    
  character(len=4)                     :: cplot
  sll_int32                            :: error
  
  npoints = size(array)

  call sll_new_file_id(file_id, error)
  call int2string(iplot,cplot)

  open(file_id,file=trim(array_name)//".gnu", &
               position="append",             &
               form='FORMATTED',iostat=error)
  if (iplot == 1) rewind(file_id)
  write(file_id,"(a)")"plot '"//trim(array_name)//cplot//".dat' with linesp"
  close(file_id)

  open(file_id,file=trim(array_name)//cplot//".dat",form='FORMATTED',iostat=error)
  do ipoints = 1, npoints
     write(file_id,*) sngl(array(ipoints))
  end do
  close(file_id)

end subroutine sll_gnuplot_write

!> Write two arrays to display with gnuplot
subroutine sll_gnuplot_write_two_arrays_1d(array_name, array1, array2, iplot)

  sll_real64, dimension(:), intent(in) :: array1     !< data
  sll_real64, dimension(:), intent(in) :: array2     !< data
  character(len=*),         intent(in) :: array_name !< field name
  sll_int32                            :: iplot      !< plot counter
  sll_int32                            :: file_id    !< file unit number
  sll_int32                            :: n
  sll_int32                            :: i
  character(len=4)                     :: cplot
  sll_int32                            :: error
  
  n = size(array1)
  SLL_ASSERT(size(array2) == n)

  call sll_new_file_id(file_id, error)
  call int2string(iplot,cplot)

  open(file_id,file=trim(array_name)//".gnu", &
               position="append",             &
               form='formatted',iostat=error)
  if (iplot == 1) rewind(file_id)
  write(file_id,"(a)")"plot '"//trim(array_name)//cplot//".dat' with linesp, &
                      & '"//trim(array_name)//cplot//".dat' u 1:3 w l"
  close(file_id)

  open(file_id,file=trim(array_name)//cplot//".dat",form='formatted',iostat=error)
  do i = 1, n
     write(file_id,*) i, sngl(array1(i)), sngl(array2(i))
  end do
  close(file_id)

end subroutine sll_gnuplot_write_two_arrays_1d


!> This subroutine write a data file to plot a 1d curve
!> @param  y_array     Y data
!> @param  x_array     X data
!> @param  array_name  field name
!> @param  iplot       Plot index 
subroutine sll_gnuplot_write_1d( y_array, x_array, array_name, iplot)

  sll_real64, dimension(:), intent(in) :: y_array    !< Y data
  sll_real64, dimension(:), intent(in) :: x_array    !< X data
  character(len=*), intent(in)         :: array_name !< field name
  sll_int32,intent(in),optional        :: iplot      !< Plot index 

  sll_int32                            :: error      !< error code
  sll_int32                            :: file_id    !< file unit number
  sll_int32                            :: fgnu_id    !< file unit number
  sll_int32                            :: npoints
  sll_int32                            :: ipoints    
  character(len=4)                     :: cplot
  
  npoints = size(x_array)

  call sll_new_file_id(file_id, error)

  if(present(iplot))then
    call int2string(iplot,cplot)
    open(file_id,file=trim(array_name)//cplot//".dat",form='FORMATTED',iostat=error)
    call sll_new_file_id(fgnu_id, error)
    open(fgnu_id,file=trim(array_name)//".gnu", &
                 position="append",             &
                 form='formatted',iostat=error)
    if (iplot == 1) rewind(file_id)
    write(fgnu_id,"(a)")"plot '"//trim(array_name)//cplot//".dat' with linesp"
    close(fgnu_id)
  else
    open(file_id,file=trim(array_name)//".dat",form='FORMATTED',iostat=error)
  endif
  rewind(file_id)
  do ipoints = 1, npoints
     write(file_id,*) sngl(x_array(ipoints)), sngl(y_array(ipoints))
  end do
  close(file_id)

end subroutine sll_gnuplot_write_1d


!> @brief
!> Write a data file plotable by gnuplot to visualize a 2d field.
!> @details
!> Axis are rectangular and spacing is constant
!> @param  xmin        Box corners
!> @param  xmax        Box corners
!> @param  ymin        Box corners
!> @param  ymax        Box corners
!> @param  nx          x points number
!> @param  ny          y points number
!> @param  array(:,:)  data
!> @param  array_name  field name
!> @param  iplot       plot counter
!> @param  error       error code
subroutine sll_gnuplot_corect_2d(xmin, xmax, nx,    &
                                 ymin, ymax, ny,    &
                                 array, array_name, &
                                 iplot, error)  

  sll_real64,       intent(in)  :: xmin
  sll_real64,       intent(in)  :: xmax
  sll_real64,       intent(in)  :: ymin
  sll_real64,       intent(in)  :: ymax
  sll_int32,        intent(in)  :: nx
  sll_int32,        intent(in)  :: ny
  sll_real64,       intent(in)  :: array(:,:)
  character(len=*), intent(in)  :: array_name
  sll_int32,        intent(in)  :: iplot
  sll_int32,        intent(out) :: error
  
  character(len=4)  :: fin   
  sll_int32, save   :: gnu_id
  sll_int32         :: file_id
  sll_int32         :: i, j
  sll_real64        :: dx, dy, x, y
  
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

!> Write a data file plotable by gnuplot to visualize a 2d field on structured
!> rectangular mesh where spacing is not constant
!> @param nx            x points number
!> @param ny            y points number
!> @param xvec(nx)      x coordinates
!> @param yvec(ny)      y coordiantes
!> @param array(nx,ny)  data
!> @param array_name    field name
!> @param iplot         plot counter
!> @param error         error code
  
subroutine sll_gnuplot_rect_2d( nx, xvec, ny, yvec,&
                               array, array_name, iplot, error)  

  sll_int32,        intent(in)  :: nx
  sll_int32,        intent(in)  :: ny
  sll_real64,       intent(in)  :: xvec(nx)
  sll_real64,       intent(in)  :: yvec(ny)
  sll_real64,       intent(in)  :: array(nx,ny)
  character(len=*), intent(in)  :: array_name
  sll_int32,        intent(in)  :: iplot
  sll_int32,        intent(out) :: error
  
  character(len=4)  :: fin   
  sll_int32, save   :: gnu_id
  sll_int32         :: file_id
  sll_int32         :: i, j
  
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

!> Write a data file plotable by gnuplot to visualize a 2d curvilinear mesh
!> @param nx          x points number
!> @param ny          y points number
!> @param xcoord      x coordinates
!> @param ycoord      y coordiantes
!> @param array_name  field name
!> @param error       error code
subroutine sll_gnuplot_mesh_2d( nx, ny, xcoord, ycoord, array_name, error)  

  sll_int32,                    intent(in)  :: nx         !< x points number
  sll_int32,                    intent(in)  :: ny         !< y points number
  sll_real64, dimension(nx,ny), intent(in)  :: xcoord     !< x coordinates
  sll_real64, dimension(nx,ny), intent(in)  :: ycoord     !< y coordiantes
  character(len=*)            , intent(in)  :: array_name !< field name
  sll_int32                   , intent(out) :: error      !< error code

  sll_int32, save :: gnu_id
  sll_int32 :: file_id
  sll_int32 :: i, j
  
  call sll_new_file_id(gnu_id, error)
  
  open(gnu_id,file=array_name//".gnu")
  write(gnu_id,*)"set view 0,0"
  write(gnu_id,*)"splot '"//array_name//"_mesh.dat' with lines"
  close(gnu_id)
  
  call sll_ascii_file_create(array_name//'_mesh.dat', file_id, error )
  do i=1,nx
     do j=1,ny
        write(file_id,*) sngl(xcoord(i,j)), &
                         sngl(ycoord(i,j)), 0.0_f32
     end do
     write(file_id,*)
  enddo
  close(file_id)

end subroutine sll_gnuplot_mesh_2d


!> @brief
!> write a data file plotable by gnuplot.
!> @details
!> We visualize a 2d field on structured curvilinear mesh
!< $param nx           x points number
!> @param ny           y points number
!> @param x(nx,ny)     x coordinates
!> @param y(nx,ny)     y coordiantes
!> @param array(nx,ny) data
!> @param array_name   field name
!> @param iplot        plot counter
!> @param error        error code
subroutine sll_gnuplot_curv_2d( nx, ny, x, y, array, array_name, iplot, error)  

  sll_int32,        intent(in)  :: nx           
  sll_int32,        intent(in)  :: ny           
  sll_real64,       intent(in)  :: x(nx,ny)
  sll_real64,       intent(in)  :: y(nx,ny)
  sll_real64,       intent(in)  :: array(nx,ny) 
  character(len=*), intent(in)  :: array_name   
  sll_int32,        intent(in)  :: iplot        
  sll_int32,        intent(out) :: error        
  
  sll_int32, save  :: gnu_id
  sll_int32        :: file_id
  sll_int32        :: i, j
  character(len=4) :: fin   
  
  call int2string(iplot, fin)
  
  if ( iplot == 1 ) then
     call sll_new_file_id(gnu_id, error)
  end if
  
  open(gnu_id,file=array_name//".gnu", position="append")
  write(gnu_id,*)"set output '"//array_name//"_"//fin//".png'"
  write(gnu_id,*)"splot '"//array_name//"_"//fin//".dat' w l"
  close(gnu_id)
  
  call sll_ascii_file_create(array_name//'_'//fin//'.dat', file_id, error )

  do i=1,nx
     do j=1,ny
        write(file_id,*) sngl(x(i,j)), &
                         sngl(y(i,j)), &
                         sngl(array(i,j)) 
     end do
     write(file_id,*)
  enddo
  close(file_id)

end subroutine sll_gnuplot_curv_2d


!> Write a field on unstructures mesh of triangles
!> @param[in] field_at_nodes field value on nodes
!> @param[in] field_name     field name use as prefix for file name
!> @param[in] coord          coordinates of nodes
!> @param[in] nodes          mesh connections
!> @param[in] plot_number    plot counter used for file name
subroutine write_unstructured_field( field_at_node,  &
                                     field_name,     &
                                     coord,          &
                                     nodes,          &
                                     plot_number     )

sll_real64, dimension(:)  , intent(in) :: field_at_node
character(len=*),           intent(in) :: field_name
sll_real64, dimension(:,:), intent(in) :: coord
sll_int32 , dimension(:,:), intent(in) :: nodes
sll_int32,                  intent(in) :: plot_number

character(len=4)                       :: cplot
sll_int32                              :: gnu_id
sll_int32                              :: num_cells
sll_int32                              :: ierr
sll_int32                              :: i
sll_real32                             :: xs1, xs2, xs3
sll_real32                             :: ys1, ys2, ys3

num_cells = size(nodes,2)
SLL_ASSERT( size(nodes,1) == 3)
call int2string(plot_number, cplot)

write(*,"(/10x, 'Output file GNUplot ',a/)") field_name//'_'//cplot//'.dat'

call sll_new_file_id(gnu_id, ierr)
open(gnu_id, file = field_name//'.gnu', position="append")
if (plot_number == 1) rewind(gnu_id)
write(gnu_id,*)"set title 'Field "//field_name//"'"
write(gnu_id,*)"splot '"//field_name//'_'//cplot//".dat' w l"
close(gnu_id)

open(gnu_id, file = field_name//'_'//cplot//'.dat')

do i = 1, num_cells

  xs1 = sngl(coord(1,nodes(1,i))); ys1 = sngl(coord(2,nodes(1,i)))
  xs2 = sngl(coord(1,nodes(2,i))); ys2 = sngl(coord(2,nodes(2,i)))
  xs3 = sngl(coord(1,nodes(3,i))); ys3 = sngl(coord(2,nodes(3,i)))

  write(gnu_id,"(3e12.3)") xs1, ys1, field_at_node(nodes(1,i))
  write(gnu_id,"(3e12.3)") xs2, ys2, field_at_node(nodes(2,i))
  write(gnu_id,"(3e12.3)") xs3, ys3, field_at_node(nodes(3,i))
  write(gnu_id,"(3e12.3)") xs1, ys1, field_at_node(nodes(1,i))
  write(gnu_id,*)
  write(gnu_id,*)

end do

close(gnu_id)

end subroutine write_unstructured_field

end module sll_gnuplot

