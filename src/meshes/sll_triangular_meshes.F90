!------------------------------------------------------------------------------
! SELALIB
!------------------------------------------------------------------------------
! MODULE: sll_triangular_meshes
!
! DESCRIPTION:
!> @file sll_triangular_meshes.F90
!>
!> @author
!> - Aurore Back
!> - Pierre Navaro
!>
!> @brief
!>  This module defines a triangular mesh.
!>
!> @details
!>  Triangular mesh 
!------------------------------------------------------------------------------
module sll_triangular_meshes
#include "sll_working_precision.h"
#include "sll_memory.h"

use sll_constants
use sll_utilities
use sll_meshes_base
use sll_tri_mesh_xmf

  implicit none

  !> @brief 2d hexagonal mesh
  type ::  sll_triangular_mesh_2d

     sll_int32  :: num_nodes  
     sll_int32  :: num_cells 
     sll_int32  :: num_edges     
     sll_real64, pointer, dimension(:,:) :: coord 
     sll_int32,  pointer, dimension(:,:) :: nodes 
     sll_real64 :: eta1_min
     sll_real64 :: eta1_max
     sll_real64 :: eta2_min
     sll_real64 :: eta2_max

!   contains

!     procedure, pass(mesh) :: eta1_node => eta1_node_hex
!     procedure, pass(mesh) :: eta2_node => eta2_node_hex
!     procedure, pass(mesh) :: eta1_cell_one_arg => eta1_cell_hex
!     procedure, pass(mesh) :: eta1_cell_two_arg => eta1_cell_triangular_two_arg
!     procedure, pass(mesh) :: eta2_cell_one_arg => eta2_cell_hex
!     procedure, pass(mesh) :: eta2_cell_two_arg => eta2_cell_triangular_two_arg
!     procedure, pass(mesh) :: display => display_triangular_mesh_2d
!     procedure, pass(mesh) :: delete => delete_triangular_mesh_2d

  end type sll_triangular_mesh_2d

  interface sll_new
     module procedure initialize_triangular_mesh_2d
  end interface sll_new

  interface sll_delete
     module procedure delete_triangular_mesh_2d
  end interface sll_delete

  interface sll_display
     module procedure display_triangular_mesh_2d
  end interface sll_display

  ! interface eta1_cell
  !    module procedure eta1_cell_one_arg, eta1_cell_two_arg
  ! end interface eta1_cell

  ! interface eta2_cell
  !    module procedure eta2_cell_one_arg, eta2_cell_two_arg
  ! end interface eta2_cell


contains

  subroutine initialize_triangular_mesh_2d( mesh,     &
                                            nc_eta1,  &
                                            eta1_min, &
                                            eta1_max, &
                                            nc_eta2,  &
                                            eta2_min, &
                                            eta2_max)


    type(sll_triangular_mesh_2d) :: mesh
    sll_int32,  intent(in) :: nc_eta1
    sll_real64, intent(in) :: eta1_min
    sll_real64, intent(in) :: eta1_max
    sll_int32,  intent(in) :: nc_eta2
    sll_real64, intent(in) :: eta2_min
    sll_real64, intent(in) :: eta2_max

 end subroutine initialize_triangular_mesh_2d

!  function eta1_node_triangular(mesh, i, j) result(res)
!    ! The coordinates (i, j) correspond to the (r1, r2) basis
!    ! This function returns the 1st coordinate on the cartesian system
!    class(sll_triangular_mesh_2d), intent(in) :: mesh
!    sll_int32, intent(in)  :: i
!    sll_int32, intent(in)  :: j
!    sll_real64 :: res
!
!    res = mesh%r1_x1*i + mesh%r2_x1*j + mesh%center_x1
!  end function eta1_node_triangular
!
!  !> @brief Computes the second coordinate of a given point
!  !> @details Computes the second coordinate on the cartesian system 
!  !> of a point which has for hexagonal coordinates (i,j)
!  !> @param i integer denoting the first hexagonal coordinate of a point
!  !> @param j integer denoting the second hexagonal coordinate of a point
!  !> returns res real containing the coordinate "eta2"
!  function eta2_node_triangular(mesh, i, j) result(res)
!    ! The coordinates (k1, k2) correspond to the (r1, r2) basis
!    ! This function the 2nd coordinate on the cartesian system
!    class(sll_triangular_mesh_2d), intent(in)     :: mesh
!    sll_int32, intent(in)  :: i
!    sll_int32, intent(in)  :: j
!    sll_real64  :: res
!
!    res = mesh%r1_x2*i + mesh%r2_x2*j + mesh%center_x2
!  end function eta2_node_triangular
!
!
!  function eta1_cell_triangular(mesh, cell_num) result(res)
!    ! The index num_ele corresponds to the index of triangle
!    ! This function returns the 1st coordinate on the cartesian system
!    ! of the center of the triangle at num_ele
!    class(sll_triangular_mesh_2d),intent(in)     :: mesh
!    sll_int32, intent(in)      :: cell_num
!    sll_real64 :: res
!
!    res = mesh%center_cartesian_coord(1, cell_num)
!  end function eta1_cell_triangular
!
!  function eta2_cell_triangular(mesh, cell_num) result(res)
!    ! The index num_ele corresponds to the index of triangle
!    ! This function returns the 2nd coordinate on the cartesian system
!    ! of the center of the triangle at num_ele
!    class(sll_triangular_mesh_2d),intent(in)     :: mesh
!    sll_int32, intent(in)      :: cell_num
!    sll_real64 :: res
!
!    res = mesh%center_cartesian_coord(2, cell_num)
!  end function eta2_cell_triangular
!
!
!  function eta1_cell_triangular_two_arg(mesh, i, j) result(res)
!    class(sll_triangular_mesh_2d),intent(in)     :: mesh
!    sll_int32, intent(in)      :: i, j
!    sll_real64 :: res
!
!    res = 0.0_f64
!    print *, "Error : eta1_cell for a triangular mesh only works with ONE parameter (num_cell)"
!    STOP
!  end function eta1_cell_triangular_two_arg
!
!  function eta2_cell_triangular_two_arg(mesh, i, j) result(res)
!    class(sll_triangular_mesh_2d),intent(in)     :: mesh
!    sll_int32, intent(in)      :: i, j
!    sll_real64 :: res
!
!    res = 0.0_f64
!    print *, "Error : eta2_cell for a triangular mesh only works with ONE parameter (num_cell)"
!    STOP
!  end function eta2_cell_triangular_two_arg
!
!
!  function cells_to_origin(k1, k2) result(val)
!    ! Takes the coordinates (k1,k2) on the (r1,r2) basis and 
!    ! returns the number of cells between that point and
!    ! the origin. If (k1, k2) = 0, val = 0
!    sll_int32, intent(in)   :: k1
!    sll_int32, intent(in)   :: k2
!    sll_int32               :: val
!
!    ! We compute the number of cells from point to center 
!    if (k1*k2 .gt. 0) then
!       val = max(abs(k1),abs(k2))
!    else
!       val = abs(k1) + abs(k2)
!    end if
!
!  end function cells_to_origin
!
!
!  function hex_to_global(mesh, k1, k2) result(val)
!    ! Takes the coordinates (k1,k2) on the (r1,r2) basis and 
!    ! returns global index of that mesh point.
!    ! By default the index of the center of the mesh is 0
!    ! Then following the r1 direction and a counter-clockwise motion
!    ! we assing an index to every point of the mesh.
!    class(sll_triangular_mesh_2d)      :: mesh
!    sll_int32, intent(in)   :: k1
!    sll_int32, intent(in)   :: k2
!    sll_int32               :: distance
!    sll_int32               :: index_tab
!    sll_int32               :: val
!    
!    distance = cells_to_origin(k1,k2)
!
!    ! Test if we are in domain
!    if (distance .le. mesh%num_cells) then
!
!       call index_triangular_to_global(mesh, k1, k2,index_tab)
!       val = mesh%global_indices(index_tab)
!
!    else
!       val = -1
!    end if
!
!  end function hex_to_global
!
!
  subroutine display_triangular_mesh_2d(mesh)
    ! Displays mesh information on the terminal
    class(sll_triangular_mesh_2d), intent(in) :: mesh

    write(*,"(/,(a))") '2D mesh : num_cells   num_nodes   num_edges '
    write(*,"(10x,3(i6,9x),4(g13.3,1x))") &
         mesh%num_cells,  &
         mesh%num_nodes,  &
         mesh%num_edges,  &
         mesh%eta1_min,   &
         mesh%eta1_max,   &
         mesh%eta2_min,   &
         mesh%eta2_max

  end subroutine display_triangular_mesh_2d
!
!
!  subroutine write_triangular_mesh_2d(mesh, name)
!    ! Writes the mesh information in a file named "name"
!    type(sll_triangular_mesh_2d), pointer :: mesh
!    character(len=*) :: name
!    sll_int32  :: i
!    sll_int32  :: num_pts_tot
!    sll_int32  :: k1, k2
!    sll_int32, parameter :: out_unit=20
!
!    open (unit=out_unit,file=name,action="write",status="replace")
!
!    num_pts_tot = mesh%num_pts_tot
!
!    ! Optional writing every mesh point and its cartesian coordinates :
!    !    write(*,"(/,(a))") 'hex mesh : num_pnt    x1     x2'
!
!    do i=1, num_pts_tot
!       k1 = mesh%global_to_hex1(i)
!       k2 = mesh%global_to_hex2(i)
!       write (out_unit, "(3(i6,1x),2(g13.3,1x))") i,                &
!            k1,                      &
!            k2,                      &
!            mesh%global_to_x1(i), &
!            mesh%global_to_x2(i)
!    end do
!
!    close(out_unit)
!  end subroutine write_triangular_mesh_2d
!
!  subroutine write_field_triangular_mesh(mesh, field, name)
!    ! Writes the points cartesian coordinates and
!    ! field(vector) values in a file named "name"
!    type(sll_triangular_mesh_2d), pointer :: mesh
!    sll_real64,dimension(:) :: field
!    character(len=*) :: name
!    sll_int32  :: i
!    sll_int32  :: num_pts_tot
!    sll_real64 :: x1, x2
!    sll_int32, parameter :: out_unit=20
!
!    open (unit=out_unit,file=name,action="write",status="replace")
!
!    num_pts_tot = mesh%num_pts_tot
!    do i=1, num_pts_tot
!       x1 = mesh%global_to_x1(i)
!       x2 = mesh%global_to_x2(i)
!       write (out_unit, "(3(g13.3,1x))") x1, &
!            x2, &
!            field(i)
!    end do
!
!    close(out_unit)
!  end subroutine write_field_triangular_mesh
!
!  subroutine write_field_triangular_mesh_xmf(mesh, field, name)
!    ! Writes the points cartesian coordinates and
!    ! field(vector) values in a file named "name"
!    type(sll_triangular_mesh_2d), pointer :: mesh
!    sll_real64,dimension(:) :: field
!    character(len=*) :: name
!    sll_int32  :: i
!    sll_int32  :: num_triangles
!    sll_int32  :: num_pts_tot
!    sll_int32  :: out_unit
!    sll_real64, allocatable :: coor(:,:)
!    sll_int32,  allocatable :: ntri(:,:)
!    sll_int32  :: error
!    sll_real64 :: x1, x2
!
!    call sll_new_file_id(out_unit, error)
!
!    num_pts_tot = mesh%num_pts_tot
!    num_triangles = mesh%num_triangles
!    SLL_ALLOCATE(coor(2,num_pts_tot),error)
!    SLL_ALLOCATE(ntri(3,num_triangles),error)
!
!    do i=1, num_pts_tot
!       coor(1,i) = mesh%global_to_x1(i)
!       coor(2,i) = mesh%global_to_x2(i)
!    end do
!
!    do i=1, num_triangles
!       x1      = mesh%center_cartesian_coord(1, i)
!       x2      = mesh%center_cartesian_coord(2, i)
!       call get_cell_vertices_index( x1, x2, mesh, ntri(1,i), ntri(2,i), ntri(3,i))
!    end do
!
!    call write_tri_mesh_xmf(name, coor, ntri, num_pts_tot, num_triangles, field, 'values')
!
!    close(out_unit)
!
!  end subroutine write_field_triangular_mesh_xmf
!
!
!  subroutine write_triangular_mesh_mtv(mesh, mtv_file)
!
!    type(sll_triangular_mesh_2d), pointer :: mesh
!    sll_real64                 :: coor(2,mesh%num_pts_tot)
!    sll_int32                  :: ntri(3,mesh%num_triangles)
!    sll_real64                 :: x1
!    sll_real64                 :: y1
!    sll_int32                  :: is1
!    sll_int32                  :: is2
!    sll_int32                  :: is3
!    character(len=*)           :: mtv_file
!    sll_int32                  :: out_unit
!    sll_int32                  :: error
!    sll_int32                  :: i
!    
!    call sll_new_file_id(out_unit, error)
!
!    open( out_unit, file=mtv_file)
!    
!    !--- Trace du maillage ---
!    
!    write(out_unit,"(a)")"$DATA=CURVE3D"
!    write(out_unit,"(a)")"%equalscale=T"
!    write(out_unit,"(a)")"%toplabel='Maillage' "
!    
!    do i = 1, mesh%num_triangles
!    
!       x1 = mesh%center_cartesian_coord(1, i)
!       y1 = mesh%center_cartesian_coord(2, i)
!    
!       call get_cell_vertices_index( x1, y1, mesh, is1, is2, is3)
!    
!       ntri(1,i) = is1
!       ntri(2,i) = is2
!       ntri(3,i) = is3
!    
!       coor(1,is1) = mesh%global_to_x1(is1)
!       coor(2,is1) = mesh%global_to_x2(is1)
!       coor(1,is2) = mesh%global_to_x1(is2)
!       coor(2,is2) = mesh%global_to_x2(is2)
!       coor(1,is3) = mesh%global_to_x1(is3)
!       coor(2,is3) = mesh%global_to_x2(is3)
!    
!       write(out_unit,"(3f10.5)")coor(:,ntri(1,i)),0.
!       write(out_unit,"(3f10.5)")coor(:,ntri(2,i)),0.
!       write(out_unit,"(3f10.5)")coor(:,ntri(3,i)),0.
!       write(out_unit,"(3f10.5)")coor(:,ntri(1,i)),0.
!       write(out_unit,*)
!    
!    end do
!    
!    !--- Numeros des noeuds et des triangles
!    
!    write(out_unit,"(a)")"$DATA=CURVE3D"
!    write(out_unit,"(a)")"%equalscale=T"
!    write(out_unit,"(a)")"%toplabel='Numeros des noeuds et des triangles' "
!    
!    do i = 1, mesh%num_triangles
!       write(out_unit,"(3f10.5)")coor(:,ntri(1,i)),0.
!       write(out_unit,"(3f10.5)")coor(:,ntri(2,i)),0.
!       write(out_unit,"(3f10.5)")coor(:,ntri(3,i)),0.
!       write(out_unit,"(3f10.5)")coor(:,ntri(1,i)),0.
!       write(out_unit,*)
!    end do
!    
!    do i = 1, mesh%num_triangles
!       x1 = (  coor(1,ntri(1,i))  &
!             + coor(1,ntri(2,i))  &
!         + coor(1,ntri(3,i))    )/3.
!       y1 = (  coor(2,ntri(1,i))  &
!             + coor(2,ntri(2,i))  &
!         + coor(2,ntri(3,i))    )/3.
!       write(out_unit,"(a)"   , advance="no")"@text x1="
!       write(out_unit,"(f8.5)", advance="no") x1
!       write(out_unit,"(a)"   , advance="no")" y1="
!       write(out_unit,"(f8.5)", advance="no") y1
!       write(out_unit,"(a)"   , advance="no")" z1=0. lc=4 ll='"
!       write(out_unit,"(i4)"  , advance="no") i
!       write(out_unit,"(a)")"'"
!    end do
!    
!    do i = 1, mesh%num_pts_tot
!       x1 = coor(1,i)
!       y1 = coor(2,i)
!       write(out_unit,"(a)"   , advance="no")"@text x1="
!       write(out_unit,"(g15.3)", advance="no") x1
!       write(out_unit,"(a)"   , advance="no")" y1="
!       write(out_unit,"(g15.3)", advance="no") y1
!       write(out_unit,"(a)"   , advance="no")" z1=0. lc=5 ll='"
!       write(out_unit,"(i4)"  , advance="no") i
!       write(out_unit,"(a)")"'"
!    end do
!    
!    !--- Numeros des noeuds 
!    
!    write(out_unit,*)"$DATA=CURVE3D"
!    write(out_unit,*)"%equalscale=T"
!    write(out_unit,*)"%toplabel='Numeros des noeuds' "
!    
!    do i = 1, mesh%num_triangles
!       write(out_unit,"(3f10.5)")coor(:,ntri(1,i)),0.
!       write(out_unit,"(3f10.5)")coor(:,ntri(2,i)),0.
!       write(out_unit,"(3f10.5)")coor(:,ntri(3,i)),0.
!       write(out_unit,"(3f10.5)")coor(:,ntri(1,i)),0.
!       write(out_unit,*)
!    end do
!    
!    do i = 1, mesh%num_pts_tot
!       x1 = coor(1,i)
!       y1 = coor(2,i)
!       write(out_unit,"(a)"   , advance="no")"@text x1="
!       write(out_unit,"(g15.3)", advance="no") x1
!       write(out_unit,"(a)"   , advance="no")" y1="
!       write(out_unit,"(g15.3)", advance="no") y1
!       write(out_unit,"(a)"   , advance="no")" z1=0. lc=5 ll='"
!       write(out_unit,"(i4)"  , advance="no") i
!       write(out_unit,"(a)")"'"
!    end do
!    
!    !--- Numeros des triangles
!    
!    write(out_unit,*)"$DATA=CURVE3D"
!    write(out_unit,*)"%equalscale=T"
!    write(out_unit,*)"%toplabel='Numeros des triangles' "
!    
!    do i = 1, mesh%num_triangles
!       write(out_unit,"(3f10.5)")coor(:,ntri(1,i)),0.
!       write(out_unit,"(3f10.5)")coor(:,ntri(2,i)),0.
!       write(out_unit,"(3f10.5)")coor(:,ntri(3,i)),0.
!       write(out_unit,"(3f10.5)")coor(:,ntri(1,i)),0.
!       write(out_unit,*)
!    end do
!    
!    do i = 1, mesh%num_triangles
!       x1 = (  coor(1,ntri(1,i))  &
!             + coor(1,ntri(2,i))  &
!         + coor(1,ntri(3,i))    )/3.
!       y1 = (  coor(2,ntri(1,i))  &
!             + coor(2,ntri(2,i))  &
!         + coor(2,ntri(3,i))    )/3.
!       write(out_unit,"(a)"   , advance="no")"@text x1="
!       write(out_unit,"(g15.3)", advance="no") x1
!       write(out_unit,"(a)"   , advance="no")" y1="
!       write(out_unit,"(g15.3)", advance="no") y1
!       write(out_unit,"(a)"   , advance="no")" z1=0. lc=4 ll='"
!       write(out_unit,"(i4)"  , advance="no") i
!       write(out_unit,"(a)")"'"
!    end do
!    
!    write(out_unit,*)"$END"
!    close(out_unit)
!   
!end subroutine write_triangular_mesh_mtv


subroutine delete_triangular_mesh_2d( mesh )
  class(sll_triangular_mesh_2d), intent(inout) :: mesh
  sll_int32 :: ierr

  print*, 'delete mesh'

end subroutine delete_triangular_mesh_2d


end module sll_triangular_meshes
