module sll_module_mapped_meshes_2d_cartesian
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
  use sll_splines
  use sll_module_mapped_meshes_2d_base
  use sll_module_interpolators_2d_base
  implicit none
 
!  HOW-TO INITIALIZE THE CARTESIAN MESH
!
!  type(sll_mapped_mesh_2d_cartesian)        :: mesh
!  call mesh%initialize(label , nx , ny )
!
!  where label is a string for identify the output ("label".xmf)
!  and nx,ny the numbers of points in x and y-axes.
!
! You have access to all functions and attributes defined in 
! sll_mapped_meshes_base
 
  type, public, extends(sll_mapped_mesh_2d_base)::sll_mapped_mesh_2d_cartesian
    private
     procedure(one_arg_scalar_function), pointer, nopass    :: x1_func  ! user
     procedure(one_arg_scalar_function), pointer, nopass    :: x2_func  ! user
     type(jacobian_matrix_element), dimension(:,:), pointer :: j_matrix
     procedure(two_arg_message_passing_func_analyt), pointer, pass :: &
          jacobian_func
     procedure(j_matrix_f_nopass), pointer, nopass :: jacobian_matrix_function
   contains
     procedure, pass(mesh) :: initialize => initialize_mesh_2d_cartesian
     procedure, pass(mesh) :: x1_at_node => x1_node_cartesian
     procedure, pass(mesh) :: x2_at_node => x2_node_cartesian
     procedure, pass(mesh) :: jacobian_at_node =>mesh_2d_jacobian_node_cartesian
     procedure, pass(mesh) :: x1_at_cell => x1_cell_cartesian
     procedure, pass(mesh) :: x2_at_cell => x2_cell_cartesian
     procedure, pass(mesh) :: jacobian_at_cell =>mesh_2d_jacobian_cell_cartesian
     procedure, pass(mesh) :: jacobian_matrix => jacobian_matrix_2d_cartesian
     procedure, pass(mesh) :: x1         => x1_cartesian
     procedure, pass(mesh) :: x2         => x2_cartesian
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
    npts1,          &
    npts2           &
   )

    class(sll_mapped_mesh_2d_cartesian), intent(inout) :: mesh
    character(len=*), intent(in)                  :: label
    sll_int32, intent(in)                         :: npts1
    sll_int32, intent(in)                         :: npts2

    procedure(one_arg_scalar_function), pointer            :: x1_func 
    procedure(one_arg_scalar_function), pointer            :: x2_func
    procedure(two_arg_scalar_function), pointer            :: j11_func
    procedure(two_arg_scalar_function), pointer            :: j12_func
    procedure(two_arg_scalar_function), pointer            :: j21_func
    procedure(two_arg_scalar_function), pointer            :: j22_func

    sll_real64 :: delta_1  ! cell spacing in eta1 
    sll_real64 :: delta_2  ! cell spacing in eta2 
    sll_real64 :: eta_1
    sll_real64 :: eta_2
    sll_real64 :: jacobian_val
    sll_int32  :: i
    sll_int32  :: j
    sll_int32  :: ierr

    mesh%label   = trim(label)
    mesh%nc_eta1 = npts1-1
    mesh%nc_eta2 = npts2-1
    delta_1       = 1.0_f64/(npts1 - 1)
    delta_2       = 1.0_f64/(npts2 - 1)
    mesh%delta_eta1   = delta_1
    mesh%delta_eta2   = delta_2

    x1_func => identity_x1
    x2_func => identity_x2
    j11_func => identity_jac11
    j12_func => identity_jac12
    j21_func => identity_jac21
    j22_func => identity_jac22

    ! Allocate the arrays for precomputed jacobians.
    SLL_ALLOCATE(mesh%jacobians_n(npts1,npts2), ierr)
    SLL_ALLOCATE(mesh%jacobians_c(npts1-1, npts2-1), ierr)

    SLL_ALLOCATE(mesh%x1_cell(npts1-1, npts2-1), ierr)
    SLL_ALLOCATE(mesh%x2_cell(npts1-1, npts2-1), ierr)

    ! Assign the transformation functions
    mesh%x1_func => x1_func
    mesh%x2_func => x2_func

    ! Fill the jacobian matrix
    SLL_ALLOCATE(mesh%j_matrix(2,2), ierr)
    mesh%j_matrix(1,1)%f => j11_func
    mesh%j_matrix(1,2)%f => j12_func
    mesh%j_matrix(2,1)%f => j21_func
    mesh%j_matrix(2,2)%f => j22_func
    mesh%jacobian_func   => jacobian_2d_cartesian
    
    ! Allocate the arrays for precomputed jacobians.
    SLL_ALLOCATE(mesh%jacobians_n(npts1,npts2), ierr)
    SLL_ALLOCATE(mesh%jacobians_c(npts1-1, npts2-1), ierr)
    
    ! Fill the values of the transformation and the jacobians at the nodes
    do j=0, npts2 - 1
       eta_2 = real(j,f64)*delta_2
       do i=0, npts1 - 1
          eta_1 = real(i,f64)*delta_1
          ! for some compiler reason, the following intermediate 
          ! variable is required, else the jacobians_n array will not
          ! be filled out properly.
          jacobian_val             = mesh%jacobian_func(eta_1,eta_2)
          mesh%jacobians_n(i+1,j+1) = jacobian_val
       end do
    end do
    
    ! Fill the values at the mid-point of the cells
    do j=0, npts2 - 2
       eta_2 = delta_2*(real(j,f64) + 0.5_f64)
       do i=0, npts1 - 2
          eta_1 = delta_1*(real(i,f64) + 0.5_f64)
          mesh%x1_cell(i+1,j+1)     = x1_func(eta_1)
          mesh%x2_cell(i+1,j+1)     = x2_func(eta_2)
          mesh%jacobians_c(i+1,j+1) = mesh%jacobian_func(eta_1,eta_2)
       end do
    end do
  end subroutine initialize_mesh_2d_cartesian

  function jacobian_2d_cartesian( mesh, eta1, eta2 ) result(val)
    sll_real64                        :: val
    class(sll_mapped_mesh_2d_cartesian) :: mesh
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    val = identity_jac(eta1,eta2)
  end function jacobian_2d_cartesian

  function jacobian_matrix_2d_cartesian( mesh, eta1, eta2 )
    sll_real64, dimension(1:2,1:2)     :: jacobian_matrix_2d_cartesian
    class(sll_mapped_mesh_2d_cartesian) :: mesh
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    sll_real64             :: j11
    sll_real64             :: j12
    sll_real64             :: j21
    sll_real64             :: j22
    j11 = identity_jac11( eta1, eta2 )
    j12 = identity_jac12( eta1, eta2 )
    j21 = identity_jac21( eta1, eta2 )
    j22 = identity_jac22( eta1, eta2 )
    ! For debugging:
    !    print *, 'jacobian_2d_analytic: '
    !    print *, j11, j12
    !    print *, j21, j22
    jacobian_matrix_2d_cartesian(1,1) = j11
    jacobian_matrix_2d_cartesian(1,2) = j12
    jacobian_matrix_2d_cartesian(2,1) = j21
    jacobian_matrix_2d_cartesian(2,2) = j22
  end function jacobian_matrix_2d_cartesian

  function mesh_2d_jacobian_node_cartesian( mesh, i, j )
    sll_real64              :: mesh_2d_jacobian_node_cartesian
    class(sll_mapped_mesh_2d_cartesian)   :: mesh
    sll_int32, intent(in)   :: i
    sll_int32, intent(in)   :: j
    sll_int32 :: num_pts_1
    sll_int32 :: num_pts_2
    num_pts_1 = mesh%nc_eta1 + 1
    num_pts_2 = mesh%nc_eta2 + 1
    SLL_ASSERT( (i .ge. 1) .and. (i .le. num_pts_1) )
    SLL_ASSERT( (j .ge. 1) .and. (j .le. num_pts_2) )
    mesh_2d_jacobian_node_cartesian = mesh%jacobians_n(i,j)
  end function mesh_2d_jacobian_node_cartesian

  function mesh_2d_jacobian_cell_cartesian( mesh, i, j )
    sll_real64              :: mesh_2d_jacobian_cell_cartesian
    class(sll_mapped_mesh_2d_cartesian)   :: mesh
    sll_int32, intent(in)   :: i
    sll_int32, intent(in)   :: j
    sll_int32 :: num_pts_1
    sll_int32 :: num_pts_2
    num_pts_1 = mesh%nc_eta1 + 1
    num_pts_2 = mesh%nc_eta2 + 1
    SLL_ASSERT( (i .ge. 1) .and. (i .le. num_pts_1) )
    SLL_ASSERT( (j .ge. 1) .and. (j .le. num_pts_2) )
    mesh_2d_jacobian_cell_cartesian = mesh%jacobians_c(i,j)
  end function mesh_2d_jacobian_cell_cartesian


  function x1_cartesian( mesh, eta1, eta2 ) result(val)
    sll_real64                         :: val
    class(sll_mapped_mesh_2d_cartesian) :: mesh
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    val = mesh%x1_func(eta1)
  end function x1_cartesian

  function x2_cartesian( mesh, eta1, eta2 ) result(val)
    sll_real64                        :: val
    class(sll_mapped_mesh_2d_cartesian) :: mesh
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    val = mesh%x2_func(eta2)
  end function x2_cartesian

  function x1_node_cartesian( mesh, i, j ) result(val)
    sll_real64             :: val
    class(sll_mapped_mesh_2d_cartesian) :: mesh
    sll_int32, intent(in) :: i
    sll_int32, intent(in) :: j
    val = real(i-1,f64)*mesh%delta_eta1
  end function x1_node_cartesian

  function x2_node_cartesian( mesh, i, j ) result(val)
    sll_real64             :: val
    class(sll_mapped_mesh_2d_cartesian) :: mesh
    sll_int32, intent(in) :: i
    sll_int32, intent(in) :: j
    val = real(j-1,f64)*mesh%delta_eta2
  end function x2_node_cartesian

 function x1_cell_cartesian( mesh, i, j ) 
    sll_real64            :: x1_cell_cartesian
    class(sll_mapped_mesh_2d_cartesian) :: mesh
    sll_int32, intent(in) :: i
    sll_int32, intent(in) :: j
    x1_cell_cartesian = (real(i-1,f64)+0.5_f64)*mesh%delta_eta1
  end function x1_cell_cartesian

  function x2_cell_cartesian( mesh, i, j )
    sll_real64            :: x2_cell_cartesian
    class(sll_mapped_mesh_2d_cartesian) :: mesh
    sll_int32, intent(in) :: i
    sll_int32, intent(in) :: j
    x2_cell_cartesian = (real(j-1,f64)+0.5_f64)*mesh%delta_eta2
  end function x2_cell_cartesian

 function x1_one_var( mesh, eta1) result(val)
    sll_real64                         :: val
    class(sll_mapped_mesh_2d_cartesian) :: mesh
    sll_real64, intent(in) :: eta1
    val = mesh%x1_func(eta1)
  end function x1_one_var



  ! identity function
  !-------------------
  ! direct mapping
  function identity_x1 ( eta1 )
    sll_real64  :: identity_x1
    sll_real64, intent(in)   :: eta1
    identity_x1 = eta1
  end function identity_x1

  function identity_x2 ( eta2 )
    sll_real64  :: identity_x2
    sll_real64, intent(in)   :: eta2
    identity_x2 = eta2
  end function identity_x2

  ! inverse mapping
  function identity_eta1 ( x1, x2 )
    sll_real64  :: identity_eta1
    sll_real64, intent(in)   :: x1
    sll_real64, intent(in)   :: x2
    identity_eta1 = x1
  end function identity_eta1

  function identity_eta2 ( x1, x2 )
    sll_real64  :: identity_eta2
    sll_real64, intent(in)   :: x1
    sll_real64, intent(in)   :: x2
    identity_eta2 = x2
  end function identity_eta2

  ! jacobian maxtrix
  function identity_jac11 ( eta1, eta2 )
    sll_real64  :: identity_jac11
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    identity_jac11 = 1.0_f64
  end function identity_jac11

    function identity_jac12 ( eta1, eta2 )
    sll_real64  :: identity_jac12
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    identity_jac12 = 0.0_f64
  end function identity_jac12

  function identity_jac21 ( eta1, eta2 )
    sll_real64  :: identity_jac21
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    identity_jac21 = 0.0_f64
  end function identity_jac21

  function identity_jac22 ( eta1, eta2 )
    sll_real64  :: identity_jac22
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    identity_jac22 = 1.0_f64
  end function identity_jac22

  ! jacobian ie determinant of jacobian matrix
  function identity_jac ( eta1, eta2 )
    sll_real64  :: identity_jac
    sll_real64, intent(in)   :: eta1
    sll_real64, intent(in)   :: eta2
    identity_jac = 1.0_f64
  end function identity_jac

end module sll_module_mapped_meshes_2d_cartesian
