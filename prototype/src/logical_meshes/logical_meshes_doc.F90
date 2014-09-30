! This file is read by doxygen software
! Change it to match with your library
! http://www.stack.nl/~dimitri/doxygen/manual/docblocks.html#fortranblocks
! To see the documentation, remove the line containing your directory 
! in file Doxyfile.in (line 691) if it is excluded. 
! Type 'make doc' in build directory.
! To check the results, open : 
! selalib/prototype/documentation/build/html/doxygen/html/namespaces.html 
! The following lines will be read by doxygen to generate documentation:


!> @namespace sll_logical_meshes 
!> @brief 
!> Logical mesh basic types
!> @author Edwin
!> @details
!>
!> <b> Modules available </b>
!>  List fortran module available
!>  - sll_logical_meshes
!>
!> <b> How to use it </b>
!> - Link with   <code>-lsll_logical_meshes</code>
!> - Add <code> use sll_logical_meshes </code>
!>
!> <b> Examples </b>
!> -Add some fortran lines to explain how ti use the library
!> \code
!!program unit_test_logical_meshes
!!  use sll_logical_meshes
!!  implicit none
!!
!!
!!  type(sll_logical_mesh_2d), pointer :: m2d
!!  type(sll_logical_mesh_3d), pointer :: m3d
!!  type(sll_logical_mesh_4d), pointer :: m4d
!!
!!  type(sll_logical_mesh_1d), pointer :: m1d_x
!!  type(sll_logical_mesh_1d), pointer :: m1d_y
!!  type(sll_logical_mesh_2d), pointer :: m2d_xy
!!
!!  type(sll_logical_mesh_2d), pointer :: mx
!!  type(sll_logical_mesh_2d), pointer :: mv
!!  type(sll_logical_mesh_4d), pointer :: mxv
!!  
!!  !default size is [0;1]
!!  m1d_x => new_logical_mesh_1d(64)
!!  m1d_y => new_logical_mesh_1d(32,eta_min=-1.0_f64,eta_max=+1.0_f64)
!!
!!  call sll_display(m1d_x)
!!  call sll_display(m1d_y)
!!
!!  !You can use two 1d meshes to create one 2d mesh
!!  m2d_xy => m1d_x * m1d_y
!!  call sll_display(m2d_xy)
!!  call delete(m2d_xy)
!!
!!  call delete(m1d_x)
!!  call delete(m1d_y)
!!  
!!  m2d => new_logical_mesh_2d(100,100)
!!  m3d => new_logical_mesh_3d(100,100,100)
!!  m4d => new_logical_mesh_4d(32,32,32,32, eta1_min=-1.0_f64, eta1_max = 2.0_f64)
!!
!!  call sll_display(m2d)
!!  call sll_display(m3d)
!!  call sll_display(m4d)
!!
!!  call delete(m2d)
!!  call delete(m3d)
!!  call delete(m4d)
!!
!!  !You can use two 2d meshes to create one 4d mesh
!!  mx => new_logical_mesh_2d(100,100, 0.0_f64, 12.56_f64, 0.0_f64, 12.56_f64)
!!  mv => new_logical_mesh_2d(64,64,-6.0_f64,6.0_f64,-6.0_f64,6.0_f64)
!!
!!  mxv => mx * mv
!!  call sll_display(mxv)
!!  call delete(mxv)
!!
!!  call delete(mx)
!!  call delete(mv)
!> \endcode
!>
