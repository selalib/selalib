module sll_m_meshes_base
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  implicit none

  public :: &
    sll_c_mesh_1d_base, &
    sll_c_mesh_2d_base, &
    sll_c_mesh_3d_base

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> @brief 1D logical mesh
  type, abstract :: sll_c_mesh_1d_base
   contains
     procedure(get_geometry_1d), deferred, pass(mesh) :: eta1_node
     procedure(get_geometry_1d), deferred, pass(mesh) :: eta1_cell
     procedure(display_mesh_1d), deferred, pass :: display
     procedure(delete_mesh_1d),  deferred, pass :: delete
  end type sll_c_mesh_1d_base

  !> @brief 2D logical mesh
  type, abstract :: sll_c_mesh_2d_base
   contains
     procedure(get_geometry_2d), deferred, pass(mesh) :: eta1_node
     procedure(get_geometry_2d), deferred, pass(mesh) :: eta2_node
     procedure(get_geometry_2d_one_arg), deferred, pass(mesh) :: eta1_cell_one_arg
     procedure(get_geometry_2d), deferred, pass(mesh) :: eta1_cell_two_arg
     procedure(get_geometry_2d_one_arg), deferred, pass(mesh) :: eta2_cell_one_arg
     procedure(get_geometry_2d), deferred, pass(mesh) :: eta2_cell_two_arg
     procedure(display_mesh_2d), deferred, pass :: display
     procedure(delete_mesh_2d),  deferred, pass :: delete
     generic, public :: eta1_cell => eta1_cell_one_arg, &
                                     eta1_cell_two_arg
     generic, public :: eta2_cell => eta2_cell_one_arg, &
                                     eta2_cell_two_arg
  end type sll_c_mesh_2d_base

  !> @brief 3D logical mesh
  type, abstract:: sll_c_mesh_3d_base
   contains
     procedure(get_geometry_3d), deferred, pass(mesh) :: eta1_node
     procedure(get_geometry_3d), deferred, pass(mesh) :: eta2_node
     procedure(get_geometry_3d), deferred, pass(mesh) :: eta3_node
     procedure(get_geometry_3d), deferred, pass(mesh) :: eta1_cell
     procedure(get_geometry_3d), deferred, pass(mesh) :: eta2_cell
     procedure(get_geometry_3d), deferred, pass(mesh) :: eta3_cell
     procedure(display_mesh_3d), deferred, pass :: display
     procedure(delete_mesh_3d),  deferred, pass :: delete
  end type sll_c_mesh_3d_base

  !> @brief 4D logical mesh
  type, abstract :: sll_mesh_4d_base
   contains
     procedure(get_geometry_4d), deferred, pass(mesh) :: eta1_node
     procedure(get_geometry_4d), deferred, pass(mesh) :: eta2_node
     procedure(get_geometry_4d), deferred, pass(mesh) :: eta3_node
     procedure(get_geometry_4d), deferred, pass(mesh) :: eta4_node
     procedure(get_geometry_4d), deferred, pass(mesh) :: eta1_cell
     procedure(get_geometry_4d), deferred, pass(mesh) :: eta2_cell
     procedure(get_geometry_4d), deferred, pass(mesh) :: eta3_cell
     procedure(get_geometry_4d), deferred, pass(mesh) :: eta4_cell
     procedure(display_mesh_4d), deferred, pass :: display
     procedure(delete_mesh_4d),  deferred, pass  :: delete
  end type sll_mesh_4d_base



  !Abstract functions for 1d
  abstract interface
     function get_geometry_1d(mesh, i) result(res)
       use sll_m_working_precision
       import sll_c_mesh_1d_base
       class(sll_c_mesh_1d_base), intent(in) :: mesh
       sll_int32, intent(in)  :: i
       sll_real64 :: res
     end function get_geometry_1d
  end interface

  abstract interface
     subroutine delete_mesh_1d(mesh)
       import sll_c_mesh_1d_base
       class(sll_c_mesh_1d_base), intent(inout) :: mesh
     end subroutine delete_mesh_1d
  end interface

  abstract interface
     subroutine display_mesh_1d(mesh)
       import sll_c_mesh_1d_base
       class(sll_c_mesh_1d_base), intent(in) :: mesh
     end subroutine display_mesh_1d
  end interface

  !Abstract functions for 2d
  abstract interface
     function get_geometry_2d(mesh, i, j) result(res)
       use sll_m_working_precision
       import sll_c_mesh_2d_base
       class(sll_c_mesh_2d_base), intent(in) :: mesh
       sll_int32, intent(in)  :: i
       sll_int32, intent(in)  :: j
       sll_real64 :: res
     end function get_geometry_2d
  end interface

  abstract interface
     function get_geometry_2d_one_arg(mesh, cell_num) result(res)
       use sll_m_working_precision
       import sll_c_mesh_2d_base
       class(sll_c_mesh_2d_base), intent(in) :: mesh
       sll_int32, intent(in)  :: cell_num
       sll_real64 :: res
     end function get_geometry_2d_one_arg
  end interface

  abstract interface
     subroutine display_mesh_2d(mesh)
       import sll_c_mesh_2d_base
       class(sll_c_mesh_2d_base), intent(in) :: mesh
     end subroutine display_mesh_2d
  end interface
  
  abstract interface
     subroutine delete_mesh_2d(mesh)
       import sll_c_mesh_2d_base
       class(sll_c_mesh_2d_base), intent(inout) :: mesh
     end subroutine delete_mesh_2d
  end interface

  !Abstract functions for 3d
  abstract interface
     function get_geometry_3d(mesh, i1, i2, i3) result(res)
       use sll_m_working_precision
       import sll_c_mesh_3d_base
       class(sll_c_mesh_3d_base), intent(in) :: mesh
       sll_int32, intent(in)  :: i1
       sll_int32, intent(in)  :: i2
       sll_int32, intent(in)  :: i3
       sll_real64 :: res
     end function get_geometry_3d
  end interface

  abstract interface
     subroutine display_mesh_3d(mesh)
       import sll_c_mesh_3d_base
       class(sll_c_mesh_3d_base), intent(in) :: mesh
     end subroutine display_mesh_3d
  end interface

  abstract interface
     subroutine delete_mesh_3d(mesh)
       import sll_c_mesh_3d_base
       class(sll_c_mesh_3d_base), intent(inout) :: mesh
     end subroutine delete_mesh_3d
  end interface

  !Abstract functions for 4d
  abstract interface
     function get_geometry_4d(mesh, k1, k2, k3, k4) result(res)
       use sll_m_working_precision
       import sll_mesh_4d_base
       class(sll_mesh_4d_base), intent(in) :: mesh
       sll_int32, intent(in)  :: k1
       sll_int32, intent(in)  :: k2
       sll_int32, intent(in)  :: k3
       sll_int32, intent(in)  :: k4
       sll_real64 :: res
     end function get_geometry_4d
  end interface

  abstract interface
     subroutine display_mesh_4d(mesh)
       import sll_mesh_4d_base
       class(sll_mesh_4d_base), intent(in) :: mesh
     end subroutine display_mesh_4d
  end interface

  abstract interface
     subroutine delete_mesh_4d(mesh)
       import sll_mesh_4d_base
       class(sll_mesh_4d_base), intent(inout) :: mesh
     end subroutine delete_mesh_4d
  end interface

end module sll_m_meshes_base
