module sll_module_mapped_meshes_2d_cartesian
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
  use sll_cubic_splines
  use sll_module_mapped_meshes_2d_base
  use sll_module_interpolators_2d_base
  implicit none
  ! This module ought to be eliminated since a cartesian transformation is
  ! an overkill. It remains in case that it can be used for testing purposes.

! Definition of a 2D uniform cartesian mesh
! In this case
!  x1(eta1) = x1_min + eta1 * (x1_max - x1_min), eta1 \in [0,1] 
!  x2(eta2) = x2_min + eta2 * (x2_max - x2_min), eta2 \in [0,1] 
!         
!  the Jacobian matrix is 
!   ( l1   0  )
!   (  0   l2 )
! and its determinant, the Jacobian, is l1*l2, with li=(xi_max - xi_min)
! 
!  HOW-TO INITIALIZE THE CARTESIAN MESH
!
!  type(sll_mapped_mesh_2d_cartesian)        :: mesh
!  call mesh%initialize(label , x1_min, x1_max, nx1 , x2_min, x2_max, nx2)
!
!  where label is a string for identify the output ("label".xmf)
!  and nx1, nx2 the numbers of points in x1 and x2-axes.
!  The computational domain is [x1_min, x1_max] x  [x2_min, x2_max]
!
! You have access to all functions and attributes defined in 
! sll_mapped_meshes_base
 
  type, public, extends(sll_mapped_mesh_2d_base)::sll_mapped_mesh_2d_cartesian
     private
     sll_real64  :: x1_min
     sll_real64  :: x2_min
     sll_real64  :: l1   ! = x1_max - x1_min
     sll_real64  :: l2   ! = x2_max - x2_min
    ! nc_eta1, nc_eta2, delta_eta1, delta_eta2 are inherited from the base class
   contains
     procedure, pass(mesh) :: initialize => initialize_mesh_2d_cartesian
     procedure, pass(mesh) :: x1_at_node => x1_node_cartesian
     procedure, pass(mesh) :: x2_at_node => x2_node_cartesian
     procedure, pass(mesh) :: jacobian_at_node =>mesh_2d_jacobian_node_cartesian
     procedure, pass(mesh) :: x1_at_cell => x1_cell_cartesian
     procedure, pass(mesh) :: x2_at_cell => x2_cell_cartesian
     procedure, pass(mesh) :: jacobian_at_cell =>mesh_2d_jacobian_cell_cartesian
     procedure, pass(mesh) :: jacobian_matrix => jacobian_matrix_2d_cartesian
     procedure, pass(mesh) :: inverse_jacobian_matrix => &
          inverse_jacobian_matrix_2d_cartesian
     procedure, pass(mesh) :: x1 => x1_cartesian 
     procedure, pass(mesh) :: x2 => x2_cartesian
     procedure, pass(mesh) :: x1_cartesian_one_arg
     procedure, pass(mesh) :: x2_cartesian_one_arg
     procedure, pass(mesh) :: jacobian   => jacobian_2d_cartesian
  end type sll_mapped_mesh_2d_cartesian

  abstract interface
     function j_matrix_f_nopass ( eta1, eta2) result(val)
       use sll_working_precision
       sll_real64, dimension(2,2)   :: val
       sll_real64   :: eta1
       sll_real64   :: eta2
     end function j_matrix_f_nopass
  end interface

   abstract interface
      function one_arg_scalar_function( eta )
        use sll_working_precision
        sll_real64             :: one_arg_scalar_function
        sll_real64, intent(in) :: eta
      end function one_arg_scalar_function
   end interface


   abstract interface
      function two_arg_message_passing_func_analyt( map, eta1, eta2 )
        use sll_working_precision
        import     :: sll_mapped_mesh_2d_cartesian
        sll_real64                      :: two_arg_message_passing_func_analyt
        class(sll_mapped_mesh_2d_cartesian) :: map
        sll_real64, intent(in)          :: eta1
        sll_real64, intent(in)          :: eta2
      end function two_arg_message_passing_func_analyt
   end interface

  ! Here we try to represent the Jacobian matrix an actual 2D array of
  ! functions. But since fortran does not allow arrays of pointers, here
  ! we define a special type that can be used as an array element.
  type jacobian_matrix_element
     procedure(two_arg_scalar_function), pointer, nopass :: f
  end type jacobian_matrix_element

contains

  !**************************************************************************
  !
  !               Functions for the cartesian general map
  !
  !**************************************************************************

  ! initialize_mapped_mesh_2D_general() allocates all the memory needed by 
  ! the 2D map. 

  subroutine initialize_mesh_2d_cartesian( &
    mesh,           &
    label,          &
    x1_min,         &
    x1_max,         &
    npts1,          &
    x2_min,         &
    x2_max,         &
    npts2           &
   )

    class(sll_mapped_mesh_2d_cartesian), intent(inout) :: mesh
    character(len=*), intent(in)                  :: label
    sll_int32, intent(in)                         :: npts1
    sll_int32, intent(in)                         :: npts2
    sll_real64, intent(in)                        :: x1_min 
    sll_real64, intent(in)                        :: x1_max 
    sll_real64, intent(in)                        :: x2_min
    sll_real64, intent(in)                        :: x2_max  

    sll_real64 :: delta_1  ! cell spacing in eta1 
    sll_real64 :: delta_2  ! cell spacing in eta2 
    !sll_real64 :: eta_1
    !sll_real64 :: eta_2
    !sll_real64 :: jacobian_val
    !sll_int32  :: i
    !sll_int32  :: j
    !sll_int32  :: ierr

    mesh%label   = trim(label)
    mesh%nc_eta1 = npts1-1
    mesh%nc_eta2 = npts2-1
    mesh%x1_min  = x1_min
    mesh%x2_min  = x2_min    
    delta_1       = 1.0_f64 / (npts1 - 1)
    delta_2       = 1.0_f64 / (npts2 - 1)
    mesh%delta_eta1   = delta_1
    mesh%delta_eta2   = delta_2
    mesh%l1 =  x1_max - x1_min
    mesh%l2 =  x2_max - x2_min

  end subroutine 

  function jacobian_2d_cartesian( mesh, eta1, eta2 ) result(val)
    class(sll_mapped_mesh_2d_cartesian) :: mesh
    sll_real64                        :: val
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    val = mesh%l1 * mesh%l2
  end function 

  function jacobian_matrix_2d_cartesian( mesh, eta1, eta2 )
    class(sll_mapped_mesh_2d_cartesian) :: mesh
    sll_real64, dimension(1:2,1:2)     :: jacobian_matrix_2d_cartesian
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    !sll_real64             :: j11
    !sll_real64             :: j12
    !sll_real64             :: j21
    !sll_real64             :: j22
   
    jacobian_matrix_2d_cartesian(1,1) = mesh%l1
    jacobian_matrix_2d_cartesian(1,2) = 0.0
    jacobian_matrix_2d_cartesian(2,1) = 0.0
    jacobian_matrix_2d_cartesian(2,2) = mesh%l2
  end function 


  function inverse_jacobian_matrix_2d_cartesian( mesh, eta1, eta2 )
    class(sll_mapped_mesh_2d_cartesian) :: mesh
    sll_real64, dimension(1:2,1:2)     :: inverse_jacobian_matrix_2d_cartesian
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    !sll_real64             :: j11
    !sll_real64             :: j12
    !sll_real64             :: j21
    !sll_real64             :: j22
    sll_real64             :: r_jac

    r_jac = 1.0_f64/(mesh%l1 * mesh%l2)
    inverse_jacobian_matrix_2d_cartesian(2,2) = mesh%l1*r_jac
    inverse_jacobian_matrix_2d_cartesian(1,2) = 0.0
    inverse_jacobian_matrix_2d_cartesian(2,1) = 0.0
    inverse_jacobian_matrix_2d_cartesian(1,1) = mesh%l2*r_jac
  end function 



  function mesh_2d_jacobian_node_cartesian( mesh, i, j )
    class(sll_mapped_mesh_2d_cartesian)   :: mesh
    sll_real64              :: mesh_2d_jacobian_node_cartesian
    sll_int32, intent(in)   :: i
    sll_int32, intent(in)   :: j
    sll_int32 :: num_pts_1
    sll_int32 :: num_pts_2
    num_pts_1 = mesh%nc_eta1 + 1
    num_pts_2 = mesh%nc_eta2 + 1
    SLL_ASSERT( (i .ge. 1) .and. (i .le. num_pts_1) )
    SLL_ASSERT( (j .ge. 1) .and. (j .le. num_pts_2) )
    mesh_2d_jacobian_node_cartesian = mesh%l1 * mesh%l2
  end function 

  function mesh_2d_jacobian_cell_cartesian( mesh, i, j )
    class(sll_mapped_mesh_2d_cartesian)   :: mesh
    sll_real64              :: mesh_2d_jacobian_cell_cartesian
    sll_int32, intent(in)   :: i
    sll_int32, intent(in)   :: j
    sll_int32 :: num_pts_1
    sll_int32 :: num_pts_2
    num_pts_1 = mesh%nc_eta1 
    num_pts_2 = mesh%nc_eta2 
    SLL_ASSERT( (i .ge. 1) .and. (i .le. num_pts_1) )
    SLL_ASSERT( (j .ge. 1) .and. (j .le. num_pts_2) )
    mesh_2d_jacobian_cell_cartesian = mesh%l1 * mesh%l2
  end function 

  function x1_cartesian( mesh, eta1, eta2 ) result(val)
    class(sll_mapped_mesh_2d_cartesian) :: mesh
    sll_real64                         :: val
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    val = mesh%x1_min + eta1 * mesh%l1
  end function 

  function x2_cartesian( mesh, eta1, eta2 ) result(val)
    class(sll_mapped_mesh_2d_cartesian) :: mesh
    sll_real64                        :: val
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    val = mesh%x2_min + eta2 * mesh%l2
  end function

  function x1_node_cartesian( mesh, i, j ) result(val)
    class(sll_mapped_mesh_2d_cartesian) :: mesh
    sll_real64             :: val
    sll_int32, intent(in) :: i
    sll_int32, intent(in) :: j
    val = mesh%x1_min + real(i-1,f64)*mesh%delta_eta1*mesh%l1
  end function

  function x2_node_cartesian( mesh, i, j ) result(val)
    class(sll_mapped_mesh_2d_cartesian) :: mesh
    sll_real64             :: val
    sll_int32, intent(in) :: i
    sll_int32, intent(in) :: j
    val = mesh%x2_min + real(j-1,f64)*mesh%delta_eta2*mesh%l2
  end function 

 function x1_cell_cartesian( mesh, i, j ) 
    sll_real64            :: x1_cell_cartesian
    class(sll_mapped_mesh_2d_cartesian) :: mesh
    sll_int32, intent(in) :: i
    sll_int32, intent(in) :: j
    x1_cell_cartesian = mesh%x1_min + (real(i-1,f64)+0.5_f64)*mesh%delta_eta1*mesh%l1
  end function 

  function x2_cell_cartesian( mesh, i, j )
    sll_real64            :: x2_cell_cartesian
    class(sll_mapped_mesh_2d_cartesian) :: mesh
    sll_int32, intent(in) :: i
    sll_int32, intent(in) :: j
    x2_cell_cartesian = mesh%x2_min + (real(j-1,f64)+0.5_f64)*mesh%delta_eta2*mesh%l2
  end function

 function x1_cartesian_one_arg( mesh, eta1) result(val)
    class(sll_mapped_mesh_2d_cartesian) :: mesh
    sll_real64                         :: val
    sll_real64, intent(in) :: eta1
    val = mesh%x1_min + eta1 * mesh%l1
  end function 

  function x2_cartesian_one_arg( mesh, eta2) result(val)
    class(sll_mapped_mesh_2d_cartesian) :: mesh
    sll_real64                         :: val
    sll_real64, intent(in) :: eta2
    val = mesh%x2_min + eta2 * mesh%l2
  end function 

end module sll_module_mapped_meshes_2d_cartesian
