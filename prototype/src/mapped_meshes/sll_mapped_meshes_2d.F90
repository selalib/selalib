module sll_module_mapped_meshes_2d
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_file_io.h"

  use sll_cubic_splines
  use sll_module_mapped_meshes_2d_base
  use sll_module_interpolators_2d_base
  implicit none
  
  ! ---------------------------------------------------------------------
  !
  !                  MESH TRANSFORMATIONS: mapping
  !
  ! ---------------------------------------------------------------------
  !
  ! The data contained in the mapping_2D type has been included in previous
  ! versions in a 'geometry' type. The name of 'mapping' seems to represent
  ! more directly what this information is.
  !
  ! A map is specified by the coordinate transformation from one coordinate
  ! system (eta1, eta2) to another (x1, x2). This transformation can be
  ! specified analytically, through the two functions:
  !
  !                     x1 = x1(eta1,eta2) 
  !                     x2 = x2(eta1,eta2)
  !
  ! Where both, eta1 and eta2 are defined in [0,1]. The same transformation 
  ! can be specified by the set of transformed points x1(i,j), x2(i,j), as
  ! two 2D arrays.
  !
  ! The transformation is also represented by the Jacobian matrix:
  !
  !                   [   partial x1(eta1,eta2)     partial x1(eta1,eta2)    ]
  !                   [ -----------------------    ------------------------- ]
  !                   [      partial eta1              partial eta2          ]
  !    J(eta1,eta2) = [                                                      ]
  !                   [   partial x2(eta1,eta2)     partial x2(eta1,eta2)    ]
  !                   [ -----------------------    ------------------------- ]
  !                   [      partial eta1              partial eta2          ]
  !
  ! Which for convenience, can have its determinant pre-evaluated at a 
  ! collection of locations. The implementation below should provide this
  ! information in the 'jacobians' array.

  type, extends(sll_mapped_mesh_2d_base)::sll_mapped_mesh_2d_analytic
     sll_real64, dimension(:,:), pointer :: x1_node   ! x1(i,j) 
     sll_real64, dimension(:,:), pointer :: x2_node   ! x2(i,j)
     type(jacobian_matrix_element), dimension(:,:), pointer :: j_matrix
     procedure(transformation_func_nopass), pointer, nopass :: x1_func  ! user
     procedure(transformation_func_nopass), pointer, nopass :: x2_func  ! user
     procedure(two_arg_message_passing_func_analyt), pointer, pass :: &
          jacobian_func
     procedure(j_matrix_f_nopass), pointer, nopass :: jacobian_matrix_function
     sll_real64, dimension(:), pointer :: params
   contains
     procedure, pass(mesh) :: initialize => initialize_mesh_2d_analytic
     ! Functions with integer arguments
     procedure, pass(mesh) :: x1_at_node => x1_node_analytic
     procedure, pass(mesh) :: x2_at_node => x2_node_analytic
     procedure, pass(mesh) :: jacobian_at_node => mesh_2d_jacobian_node_analytic
     procedure, pass(mesh) :: x1_at_cell => x1_cell_analytic
     procedure, pass(mesh) :: x2_at_cell => x2_cell_analytic
     procedure, pass(mesh) :: jacobian_at_cell => jacobian_2d_cell_analytic
     ! Functions with real arguments
     procedure, pass(mesh) :: x1         => x1_analytic
     procedure, pass(mesh) :: x2         => x2_analytic
     procedure, pass(mesh) :: jacobian   => jacobian_2d_analytic
     procedure, pass(mesh) :: jacobian_matrix => jacobian_matrix_2d_analytic
     procedure, pass(mesh) :: inverse_jacobian_matrix => &
          inverse_jacobian_matrix_2d_analytic
  end type sll_mapped_mesh_2d_analytic



  ! Discrete case

  type, extends(sll_mapped_mesh_2d_base)::sll_mapped_mesh_2d_discrete
     sll_real64, dimension(:,:), pointer :: x1_node   ! x1(i,j) 
     sll_real64, dimension(:,:), pointer :: x2_node   ! x2(i,j) 
     type(jacobian_matrix_element), dimension(:,:), pointer :: j_matrix
     class(sll_interpolator_2d_base), pointer                   :: x1_interp
     class(sll_interpolator_2d_base), pointer                   :: x2_interp
     procedure(two_arg_scalar_function), pointer, nopass    :: x1_func
     procedure(two_arg_scalar_function), pointer, nopass    :: x2_func
     !procedure(two_arg_message_passing_func_discr),pointer,pass :: jacobian_func
   contains
     procedure, pass(mesh) :: initialize => initialize_mesh_2d_discrete
     procedure, pass(mesh) :: x1_at_node => x1_node_discrete
     procedure, pass(mesh) :: x2_at_node => x2_node_discrete
     procedure, pass(mesh) :: jacobian_at_node => mesh_2d_jacobian_node_discrete
     procedure, pass(mesh) :: x1         => x1_discrete
     procedure, pass(mesh) :: x2         => x2_discrete
     procedure, pass(mesh) :: x1_at_cell => x1_cell_discrete
     procedure, pass(mesh) :: x2_at_cell => x2_cell_discrete
     procedure, pass(mesh) :: jacobian_at_cell => jacobian_2d_cell_discrete
     procedure, pass(mesh) :: jacobian   => jacobian_2d_discrete
     procedure, pass(mesh) :: jacobian_matrix => jacobian_matrix_2d_discrete
     procedure, pass(mesh) :: inverse_jacobian_matrix => &
          inverse_jacobian_matrix_2d_discrete
  end type sll_mapped_mesh_2d_discrete

  abstract interface
     function j_matrix_f_nopass ( eta1, eta2) result(val)
       use sll_working_precision
       sll_real64, dimension(2,2)   :: val
       sll_real64   :: eta1
       sll_real64   :: eta2
     end function j_matrix_f_nopass
  end interface

   abstract interface
      function two_arg_message_passing_func_discr( map, eta1, eta2 )
        use sll_working_precision
        import     :: sll_mapped_mesh_2d_discrete
        sll_real64                      :: two_arg_message_passing_func_discr
        class(sll_mapped_mesh_2d_discrete)  :: map
        sll_real64, intent(in)          :: eta1
        sll_real64, intent(in)          :: eta2
      end function two_arg_message_passing_func_discr
   end interface

   abstract interface
      function two_arg_message_passing_func_analyt( map, eta1, eta2 )
        use sll_working_precision
        import     :: sll_mapped_mesh_2d_analytic
        sll_real64                      :: two_arg_message_passing_func_analyt
        class(sll_mapped_mesh_2d_analytic) :: map
        sll_real64, intent(in)          :: eta1
        sll_real64, intent(in)          :: eta2
      end function two_arg_message_passing_func_analyt
   end interface

  ! Here we try to represent the Jacobian matrix an actual 2D array of
  ! functions. But since fortran does not allow arrays of pointers, here
  ! we define a special type that can be used as an array element.
  type jacobian_matrix_element
     procedure(transformation_func_nopass), pointer, nopass :: f
  end type jacobian_matrix_element
  
#if 0
  interface delete
     module procedure delete_mapped_mesh_2D_general
  end interface
#endif
  
contains

  !**************************************************************************
  !
  !               Functions for the analytic general map
  !
  !**************************************************************************

  ! initialize_mapped_mesh_2D_general() allocates all the memory needed by 
  ! the 2D map. 

  function new_mesh_2d_analytic( &
    label,          &
    npts1,          &
    npts2,          &
    x1_func,        &
    x2_func,        &
    j11_func,       &
    j12_func,       &
    j21_func,       &
    j22_func,       &
    params )

    type(sll_mapped_mesh_2d_analytic), pointer :: new_mesh_2d_analytic
    character(len=*), intent(in)                  :: label
    sll_int32, intent(in)                         :: npts1
    sll_int32, intent(in)                         :: npts2
    procedure(transformation_func_nopass)            :: x1_func
    procedure(transformation_func_nopass)            :: x2_func
    procedure(transformation_func_nopass)            :: j11_func
    procedure(transformation_func_nopass)            :: j12_func
    procedure(transformation_func_nopass)            :: j21_func
    procedure(transformation_func_nopass)            :: j22_func
    sll_real64, dimension(:), intent(in) :: params
    sll_int32 :: ierr

    SLL_ALLOCATE(new_mesh_2d_analytic, ierr)
    call initialize_mesh_2d_analytic( &
         new_mesh_2d_analytic, &
         label,          &
         npts1,          &
         npts2,          &
         x1_func,        &
         x2_func,        &
         j11_func,       &
         j12_func,       &
         j21_func,       &
         j22_func,       &
         params )
  end function new_mesh_2d_analytic

  subroutine initialize_mesh_2d_analytic( &
    mesh,           &
    label,          &
    npts1,          &
    npts2,          &
    x1_func,        &
    x2_func,        &
    j11_func,       &
    j12_func,       &
    j21_func,       &
    j22_func,       &
    params )

    class(sll_mapped_mesh_2d_analytic), intent(inout) :: mesh
    character(len=*), intent(in)                  :: label
    sll_int32, intent(in)                         :: npts1
    sll_int32, intent(in)                         :: npts2
    procedure(transformation_func_nopass)            :: x1_func
    procedure(transformation_func_nopass)            :: x2_func
    procedure(transformation_func_nopass)            :: j11_func
    procedure(transformation_func_nopass)            :: j12_func
    procedure(transformation_func_nopass)            :: j21_func
    procedure(transformation_func_nopass)            :: j22_func
    sll_real64, dimension(:), intent(in) :: params
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

    SLL_ALLOCATE(mesh%params(size(params)),ierr)
    mesh%params(:) = params(:)

    ! Allocate the arrays for precomputed jacobians.
    SLL_ALLOCATE(mesh%jacobians_n(npts1,npts2), ierr)
    SLL_ALLOCATE(mesh%jacobians_c(npts1-1, npts2-1), ierr)

    ! Allocation for x1 and x2 at nodes, needed regardless of the type of map
    SLL_ALLOCATE(mesh%x1_node(npts1,npts2), ierr)
    SLL_ALLOCATE(mesh%x2_node(npts1,npts2), ierr)

    ! Start filling out the fields and allocating the object's memory.
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
    mesh%jacobian_func   => jacobian_2d_analytic
    
    ! Allocate the arrays for precomputed jacobians.
    SLL_ALLOCATE(mesh%jacobians_n(npts1,npts2), ierr)
    SLL_ALLOCATE(mesh%jacobians_c(npts1-1, npts2-1), ierr)
    
    ! Fill the values of the transformation and the jacobians at the nodes
    do j=0, npts2 - 1
       eta_2 = real(j,f64)*delta_2
       do i=0, npts1 - 1
          eta_1 = real(i,f64)*delta_1
          mesh%x1_node(i+1,j+1) = x1_func(eta_1, eta_2, params)
          mesh%x2_node(i+1,j+1) = x2_func(eta_1, eta_2, params)
          ! for some compiler reason, the following intermediate 
          ! variable is required, else the jacobians_n array will not
          ! be filled out properly.
          jacobian_val          = mesh%jacobian_func(eta_1,eta_2)
          mesh%jacobians_n(i+1,j+1) = jacobian_val
       end do
    end do
    
    ! Fill the values at the mid-point of the cells
    do j=0, npts2 - 2
       eta_2 = delta_2*(real(j,f64) + 0.5_f64)
       do i=0, npts1 - 2
          eta_1 = delta_1*(real(i,f64) + 0.5_f64)
          mesh%x1_cell(i+1,j+1)     = x1_func(eta_1, eta_2, params)
          mesh%x2_cell(i+1,j+1)     = x2_func(eta_1, eta_2, params)
          mesh%jacobians_c(i+1,j+1) = mesh%jacobian_func(eta_1,eta_2)
       end do
    end do
  end subroutine initialize_mesh_2d_analytic

  function jacobian_2d_analytic( mesh, eta1, eta2 ) result(val)
    sll_real64                        :: val
    class(sll_mapped_mesh_2d_analytic) :: mesh
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    sll_real64             :: j11
    sll_real64             :: j12
    sll_real64             :: j21
    sll_real64             :: j22
    j11 = (mesh%j_matrix(1,1)%f( eta1, eta2, mesh%params ))
    j12 = (mesh%j_matrix(1,2)%f( eta1, eta2, mesh%params ))
    j21 = (mesh%j_matrix(2,1)%f( eta1, eta2, mesh%params ))
    j22 = (mesh%j_matrix(2,2)%f( eta1, eta2, mesh%params ))
    ! For debugging:
    !    print *, 'jacobian_2d_analytic: '
    !    print *, j11, j12
    !    print *, j21, j22
    val = j11*j22 - j12*j21
  end function jacobian_2d_analytic

  ! The efficiency of the following function could be improved by just
  ! passing the output array rather than returning values on the stack which
  ! need to be caught by the caller.
  function jacobian_matrix_2d_analytic( mesh, eta1, eta2 )
    sll_real64, dimension(1:2,1:2)     :: jacobian_matrix_2d_analytic
    class(sll_mapped_mesh_2d_analytic) :: mesh
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    sll_real64             :: j11
    sll_real64             :: j12
    sll_real64             :: j21
    sll_real64             :: j22
    j11 = (mesh%j_matrix(1,1)%f( eta1, eta2, mesh%params ))
    j12 = (mesh%j_matrix(1,2)%f( eta1, eta2, mesh%params ))
    j21 = (mesh%j_matrix(2,1)%f( eta1, eta2, mesh%params ))
    j22 = (mesh%j_matrix(2,2)%f( eta1, eta2, mesh%params ))
    ! For debugging:
    !    print *, 'jacobian_2d_analytic: '
    !    print *, j11, j12
    !    print *, j21, j22
    jacobian_matrix_2d_analytic(1,1) = j11
    jacobian_matrix_2d_analytic(1,2) = j12
    jacobian_matrix_2d_analytic(2,1) = j21
    jacobian_matrix_2d_analytic(2,2) = j22
  end function jacobian_matrix_2d_analytic

  function inverse_jacobian_matrix_2d_analytic( mesh, eta1, eta2 )
    sll_real64, dimension(1:2,1:2)     :: inverse_jacobian_matrix_2d_analytic
    class(sll_mapped_mesh_2d_analytic) :: mesh
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    sll_real64             :: inv_j11
    sll_real64             :: inv_j12
    sll_real64             :: inv_j21
    sll_real64             :: inv_j22
    sll_real64             :: r_jacobian ! reciprocal of the jacobian

    r_jacobian = 1.0_f64/mesh%jacobian( eta1, eta2 )
    inv_j11 =  (mesh%j_matrix(2,2)%f( eta1, eta2, mesh%params ))*r_jacobian
    inv_j12 = -(mesh%j_matrix(1,2)%f( eta1, eta2, mesh%params ))*r_jacobian
    inv_j21 = -(mesh%j_matrix(2,1)%f( eta1, eta2, mesh%params ))*r_jacobian
    inv_j22 =  (mesh%j_matrix(1,1)%f( eta1, eta2, mesh%params ))*r_jacobian
    ! For debugging:
    !    print *, 'jacobian_2d_analytic: '
    !    print *, j11, j12
    !    print *, j21, j22
    inverse_jacobian_matrix_2d_analytic(1,1) = inv_j11
    inverse_jacobian_matrix_2d_analytic(1,2) = inv_j12
    inverse_jacobian_matrix_2d_analytic(2,1) = inv_j21
    inverse_jacobian_matrix_2d_analytic(2,2) = inv_j22
  end function inverse_jacobian_matrix_2d_analytic

  function x1_analytic( mesh, eta1, eta2 ) result(val)
    sll_real64                         :: val
    class(sll_mapped_mesh_2d_analytic) :: mesh
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    val = mesh%x1_func(eta1, eta2,mesh%params)
  end function x1_analytic

  function x2_analytic( mesh, eta1, eta2 ) result(val)
    sll_real64                        :: val
    class(sll_mapped_mesh_2d_analytic) :: mesh
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    val = mesh%x2_func(eta1, eta2, mesh%params)
  end function x2_analytic

  function x1_node_analytic( mesh, i, j ) result(val)
    class(sll_mapped_mesh_2d_analytic) :: mesh
    sll_real64             :: val
    sll_int32, intent(in) :: i
    sll_int32, intent(in) :: j
    val = mesh%x1_node(i,j)
  end function x1_node_analytic

  function x2_node_analytic( mesh, i, j ) result(val)
    class(sll_mapped_mesh_2d_analytic) :: mesh
    sll_real64             :: val
    sll_int32, intent(in) :: i
    sll_int32, intent(in) :: j
    val = mesh%x2_node(i,j)
  end function x2_node_analytic

  function x1_cell_analytic( mesh, i, j ) result(var)
    class(sll_mapped_mesh_2d_analytic) :: mesh
    sll_real64            :: var
    sll_int32, intent(in) :: i
    sll_int32, intent(in) :: j
    var = mesh%x1_cell(i,j)
  end function x1_cell_analytic

  function x2_cell_analytic( mesh, i, j ) result(var)
    class(sll_mapped_mesh_2d_analytic) :: mesh
    sll_real64            :: var
    sll_int32, intent(in) :: i
    sll_int32, intent(in) :: j
    var = mesh%x2_cell(i,j)
  end function x2_cell_analytic

  function jacobian_2d_cell_analytic( mesh, i, j ) result(var)
    class(sll_mapped_mesh_2d_analytic) :: mesh
    sll_real64            :: var
    sll_int32, intent(in) :: i
    sll_int32, intent(in) :: j
    var = mesh%jacobians_c(i,j)
  end function jacobian_2d_cell_analytic

  !**************************************************************************
  !
  !               Functions for the discrete general map
  !
  !**************************************************************************

  function x1_node_discrete( mesh, i, j ) result(val)
    class(sll_mapped_mesh_2d_discrete) :: mesh
    sll_real64             :: val
    sll_int32, intent(in) :: i
    sll_int32, intent(in) :: j
    val = mesh%x1_node(i,j)
  end function x1_node_discrete

  function x2_node_discrete( mesh, i, j ) result(val)
    class(sll_mapped_mesh_2d_discrete) :: mesh
    sll_real64             :: val
    sll_int32, intent(in) :: i
    sll_int32, intent(in) :: j
    val = mesh%x2_node(i,j)
  end function x2_node_discrete

  function x1_cell_discrete( mesh, i, j ) result(var)
    class(sll_mapped_mesh_2d_discrete) :: mesh
    sll_real64                         :: var
    sll_int32, intent(in)              :: i
    sll_int32, intent(in)              :: j
    var = mesh%x1_cell(i,j)
  end function x1_cell_discrete

  function x2_cell_discrete( mesh, i, j ) result(var)
    class(sll_mapped_mesh_2d_discrete) :: mesh
    sll_real64                         :: var
    sll_int32, intent(in)              :: i
    sll_int32, intent(in)              :: j
    var = mesh%x2_cell(i,j)
  end function x2_cell_discrete

  function x1_discrete( mesh, eta1, eta2 ) result(val)
    class(sll_mapped_mesh_2d_discrete) :: mesh
    sll_real64             :: val
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    val = mesh%x1_interp%interpolate_value(eta1, eta2)
  end function x1_discrete

  function x2_discrete( mesh, eta1, eta2 ) result(val)
    class(sll_mapped_mesh_2d_discrete) :: mesh
    sll_real64             :: val
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    val = mesh%x2_interp%interpolate_value(eta1, eta2)
  end function x2_discrete

  function jacobian_2d_discrete( mesh, eta1, eta2 ) result(jac)
    class(sll_mapped_mesh_2d_discrete) :: mesh
    sll_real64             :: jac
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    sll_real64             :: j11
    sll_real64             :: j12
    sll_real64             :: j21
    sll_real64             :: j22
    j11 = mesh%x1_interp%interpolate_derivative_eta1( eta1, eta2 )
    j12 = mesh%x1_interp%interpolate_derivative_eta2( eta1, eta2 )
    j21 = mesh%x2_interp%interpolate_derivative_eta1( eta1, eta2 )
    j22 = mesh%x2_interp%interpolate_derivative_eta2( eta1, eta2 )
    ! For debugging:
    !    print *, 'jacobian_2D_discrete: '
    !    print *, j11, j12
    !    print *, j21, j22
    jac = j11*j22 - j12*j21
  end function jacobian_2d_discrete

  function jacobian_2d_cell_discrete( mesh, i, j ) result(var)
    class(sll_mapped_mesh_2d_discrete) :: mesh
    sll_real64                         :: var
    sll_int32, intent(in)              :: i
    sll_int32, intent(in)              :: j
    var = mesh%jacobians_c(i,j)
  end function jacobian_2d_cell_discrete

  function jacobian_matrix_2d_discrete( mesh, eta1, eta2 )
    class(sll_mapped_mesh_2d_discrete) :: mesh
    sll_real64, dimension(1:2,1:2)     :: jacobian_matrix_2d_discrete
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    sll_real64             :: j11
    sll_real64             :: j12
    sll_real64             :: j21
    sll_real64             :: j22
    j11 = mesh%x1_interp%interpolate_derivative_eta1( eta1, eta2 )
    j12 = mesh%x1_interp%interpolate_derivative_eta2( eta1, eta2 )
    j21 = mesh%x2_interp%interpolate_derivative_eta1( eta1, eta2 )
    j22 = mesh%x2_interp%interpolate_derivative_eta2( eta1, eta2 )
    ! For debugging:
    !    print *, 'jacobian_2D_discrete: '
    !    print *, j11, j12
    !    print *, j21, j22
    jacobian_matrix_2d_discrete(1,1) = j11
    jacobian_matrix_2d_discrete(1,2) = j12
    jacobian_matrix_2d_discrete(2,1) = j21
    jacobian_matrix_2d_discrete(2,2) = j22
  end function jacobian_matrix_2d_discrete

  function inverse_jacobian_matrix_2d_discrete( mesh, eta1, eta2 )
    class(sll_mapped_mesh_2d_discrete) :: mesh
    sll_real64, dimension(1:2,1:2)     :: inverse_jacobian_matrix_2d_discrete
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    sll_real64             :: inv_j11
    sll_real64             :: inv_j12
    sll_real64             :: inv_j21
    sll_real64             :: inv_j22
    sll_real64             :: r_jac ! reciprocal of the jacobian
    r_jac = 1.0_f64/mesh%jacobian( eta1, eta2 )
    inv_j11 = mesh%x1_interp%interpolate_derivative_eta1( eta1, eta2 )
    inv_j12 = mesh%x1_interp%interpolate_derivative_eta2( eta1, eta2 )
    inv_j21 = mesh%x2_interp%interpolate_derivative_eta1( eta1, eta2 )
    inv_j22 = mesh%x2_interp%interpolate_derivative_eta2( eta1, eta2 )
    ! For debugging:
    !    print *, 'jacobian_2D_discrete: '
    !    print *, j11, j12
    !    print *, j21, j22
    inverse_jacobian_matrix_2d_discrete(1,1) =  inv_j22*r_jac
    inverse_jacobian_matrix_2d_discrete(1,2) = -inv_j12*r_jac
    inverse_jacobian_matrix_2d_discrete(2,1) = -inv_j21*r_jac
    inverse_jacobian_matrix_2d_discrete(2,2) =  inv_j11*r_jac
  end function inverse_jacobian_matrix_2d_discrete



  subroutine initialize_mesh_2d_discrete( &
    mesh,            &
    label,            &
    npts1,          &
    npts2,          &
    x1_node,        &
    x2_node,        &
    x1_interpolator, &
    x2_interpolator, &
    jacobians_n_interpolator, &
    jacobians_node, &
    x1_cell, &
    x2_cell, &
    jacobians_cell )

    class(sll_mapped_mesh_2d_discrete)    :: mesh

    character(len=*), intent(in)         :: label
    sll_int32, intent(in)                :: npts1
    sll_int32, intent(in)                :: npts2
    sll_real64, dimension(:,:)           :: x1_node
    sll_real64, dimension(:,:)           :: x2_node
    class(sll_interpolator_2d_base), target  :: x1_interpolator
    class(sll_interpolator_2d_base), target  :: x2_interpolator
    class(sll_interpolator_2d_base), target  :: jacobians_n_interpolator
    sll_real64, dimension(:,:), optional :: jacobians_node
    sll_real64, dimension(:,:), optional :: jacobians_cell
    sll_real64, dimension(:,:), optional :: x1_cell
    sll_real64, dimension(:,:), optional :: x2_cell

    sll_real64 :: eta_1
    sll_real64 :: eta_2
    sll_real64 :: jacobian_val
    sll_int32  :: i
    sll_int32  :: j
    sll_int32  :: ierr
    logical    :: x1c
    logical    :: x2c
    logical    :: jc
    logical    :: jn

    mesh%label = trim(label)
    x1c = present(x1_cell)
    x2c = present(x2_cell)
    jc  = present(jacobians_cell)
    jn  = present(jacobians_node)

    ! Check argument consistency
    ! DISCRETE_MAPs require only some of the parameters. If the mapping is
    ! defined from the nodes of the logical (eta1, eta2) mesh to the nodes
    ! of the physical mesh (x1,x2), then the node arrays are required:
    ! jacobians_node, x1_node and x2_node.
    !
    ! If the transformation is done on the points at the center of the cells
    ! then these parameters are also required: 
    ! jacobians_cell, x1_cell, x2_cell.
    ! node and cell values are not mutually exclusive, thus all 6 parameters 
    ! can be provided in the NUMERIC case. It is up to the caller to make
    ! sure that the data set is consistent.

    ! 1. Check that the discrete representation of the transformation is
    !    consistent with the size of the 2D array.
    if( &
       (size(x1_node,1) .lt. npts1) .or. &
       (size(x1_node,2) .lt. npts2) ) then
       print *, 'ERROR, initialize_mapped_mesh_2d_discrete(), ', &
            'the size of the x1_node or x2_node arrays is ', &
            'inconsistent with the number of points declared, ', &
            'npts1 or npts2.'
       STOP
    end if
    if( jn .eqv. .true. ) then
       if( &
          (size(jacobians_node,1) .lt. npts1 - 1 ) .or. &
          (size(jacobians_node,2) .lt. npts2 - 1 ) ) then
          print *, 'ERROR, initialize_mapped_mesh_2d_discrete(): ', &
               'the size of the jacobians_node array is ', &
               'inconsistent with the number of points declared, ', &
               'npts1 or npts2.'
          STOP
       end if
    end if
    if( jc .eqv. .true. ) then
       if( &
          (size(jacobians_cell,1) .lt. npts1 - 1 ) .or. &
          (size(jacobians_cell,2) .lt. npts2 - 1 ) ) then
          print *, 'ERROR, initialize_mapped_mesh_2d_discrete(): ', &
               'the size of the jacobians_cell arrays is ', &
               'inconsistent with the number of points declared, ', &
               'npts1 or npts2.'
          STOP
       end if
    end if

    mesh%nc_eta1  = npts1-1
    mesh%nc_eta2  = npts2-1
    mesh%delta_eta1    = 1.0_f64/(npts1 - 1)
    mesh%delta_eta2    = 1.0_f64/(npts2 - 1)
    mesh%x1_interp => x1_interpolator
    mesh%x2_interp => x2_interpolator

    ! Allocate the arrays for precomputed jacobians.
    SLL_ALLOCATE(mesh%jacobians_n(npts1,npts2), ierr)
    SLL_ALLOCATE(mesh%jacobians_c(npts1-1, npts2-1), ierr)

    ! Allocation for x1 and x2 at nodes
    SLL_ALLOCATE(mesh%x1_node(npts1,npts2), ierr)
    SLL_ALLOCATE(mesh%x2_node(npts1,npts2), ierr)

    ! Allocation for x1 and x2 at cells
    SLL_ALLOCATE(mesh%x1_cell(npts1-1,npts2-1), ierr)
    SLL_ALLOCATE(mesh%x2_cell(npts1-1,npts2-1), ierr)

    ! initialize the local arrays. Note that since the map has its
    ! own copies, it owns this information locally and will destroy
    ! this information when the object is deleted. The caller is
    ! thus responsible for deallocating the arrays that were passed as
    ! arguments.
    do j=1, npts2
       do i=1, npts1
          mesh%x1_node(i,j) = x1_node(i,j)
          mesh%x2_node(i,j) = x2_node(i,j)
       end do
    end do

    ! Compute the spline coefficients
    call x1_interpolator%compute_interpolants( mesh%x1_node )
    call x2_interpolator%compute_interpolants( mesh%x2_node )


    ! The splines contain all the information to compute the
    ! jacobians everywhere; however, here we explore assigning
    ! the jacobians-at-the-nodes array with the values provided
    ! by the user if available. If there are discrepancies between
    ! the user-provided values and the predictions from the splines,
    ! then this may itself be a way to look for errors.
    !
    ! Copy the values of the jacobians at the nodes if user given:
    if( jn .eqv. .true. ) then
       do j=1, npts2
          do i=1, npts1
             mesh%jacobians_n(i,j) = jacobians_node(i,j)
          end do
       end do
    else
       ! Fill the jacobian values at the nodes calculated from the splines
       do j=0, npts2 - 1
          eta_2 = real(j,f64)*mesh%delta_eta2          
          do i=0, npts1 - 1
             eta_1 = real(i,f64)*mesh%delta_eta1
             jacobian_val = mesh%jacobian(eta_1,eta_2)
             mesh%jacobians_c(i+1,j+1) = jacobian_val
          end do
       end do
    end if

    ! copy the cell-based transformation arrays if available
    if( (x1c .and. x2c) .eqv. .true. ) then
       SLL_ALLOCATE(mesh%x1_cell(npts1-1, npts2-1), ierr)
       SLL_ALLOCATE(mesh%x2_cell(npts1-1, npts2-1), ierr)
       do j=1, npts2 - 1
          do i=1, npts1 - 1
             mesh%x1_cell(i,j) = x1_cell(i,j)
             mesh%x2_cell(i,j) = x2_cell(i,j)
          end do
       end do
    end if
    ! copy the cell-based jacobians if available
    if( jc .eqv. .true. ) then
       do j=1, npts2 - 1
          do i=1, npts1 - 1
             mesh%jacobians_c(i,j) = jacobians_cell(i,j)
          end do
       end do
    else ! if cell-based jacobians are not available, compute them.
    ! Fill the values at the mid-point of the cells
       do j=0, npts2 - 2
          eta_2 = mesh%delta_eta2*(real(j,f64) + 0.5_f64)
          do i=0, npts1 - 2
             ! it is very bad practice to invoke the mesh methods while
             ! we are not even done initializing the mesh object...
             eta_1 = mesh%delta_eta1*(real(i,f64) + 0.5_f64)
             mesh%x1_cell(i+1,j+1)     = mesh%x1(eta_1, eta_2)
             mesh%x2_cell(i+1,j+1)     = mesh%x2(eta_1, eta_2)
             mesh%jacobians_c(i+1,j+1) = mesh%jacobian(eta_1,eta_2)
          end do
       end do
    end if
  end subroutine initialize_mesh_2d_discrete

  function mesh_2d_jacobian_node_analytic( mesh, i, j )
    class(sll_mapped_mesh_2d_analytic)   :: mesh
    sll_real64              :: mesh_2d_jacobian_node_analytic
    sll_int32, intent(in)   :: i
    sll_int32, intent(in)   :: j
    sll_int32 :: num_pts_1
    sll_int32 :: num_pts_2
    num_pts_1 = mesh%nc_eta1 + 1
    num_pts_2 = mesh%nc_eta2 + 1
    SLL_ASSERT( (i .ge. 1) .and. (i .le. num_pts_1) )
    SLL_ASSERT( (j .ge. 1) .and. (j .le. num_pts_2) )
    mesh_2d_jacobian_node_analytic = mesh%jacobians_n(i,j)
  end function mesh_2d_jacobian_node_analytic

  function mesh_2d_jacobian_node_discrete( mesh, i, j )
    class(sll_mapped_mesh_2d_discrete)   :: mesh
    sll_real64              :: mesh_2d_jacobian_node_discrete
    sll_int32, intent(in)   :: i
    sll_int32, intent(in)   :: j
    sll_int32 :: num_pts_1
    sll_int32 :: num_pts_2
    num_pts_1 = mesh%nc_eta1 + 1
    num_pts_2 = mesh%nc_eta2 + 1
    SLL_ASSERT( (i .ge. 1) .and. (i .le. num_pts_1) )
    SLL_ASSERT( (j .ge. 1) .and. (j .le. num_pts_2) )
    mesh_2d_jacobian_node_discrete = mesh%jacobians_n(i,j)
  end function mesh_2d_jacobian_node_discrete

#if 0
  subroutine delete_mapped_mesh_2D_general( map )
    type(mapped_mesh_2D_general), pointer :: map
    sll_int32             :: ierr
    if( .not. associated(map) ) then
       print *, 'ERROR, delete_mapped_mesh_2D_general: passed map pointer ', &
            'is not associated.'
    end if
    if( associated(map%x1_node) ) then
       SLL_DEALLOCATE( map%x1_node, ierr )
    end if
    if( associated(map%x2_node) ) then
       SLL_DEALLOCATE( map%x2_node, ierr )
    end if
    if( associated(map%x1_cell) ) then
       SLL_DEALLOCATE( map%x1_cell, ierr )
    end if
    if( associated(map%x2_cell) ) then
       SLL_DEALLOCATE( map%x2_cell, ierr )
    end if
    if( associated(map%j_matrix) ) then
       SLL_DEALLOCATE( map%j_matrix, ierr )
    end if
    SLL_DEALLOCATE( map%jacobians_n, ierr )
    SLL_DEALLOCATE( map%jacobians_c, ierr )
    if( map%map_type .eq. DISCRETE_MAP ) then
       call delete(map%x1_spline)
       call delete(map%x2_spline)
    end if
    SLL_DEALLOCATE( map, ierr )
  end subroutine delete_mapped_mesh_2D_general


  ! Access functions for the mapping. These are an overkill and can be 
  ! changed by a macro, but for now, they are at least safer.
  function mesh_2d_x1_node( map, i, j )
    sll_real64            :: mesh_2d_x1_node
    type(mapped_mesh_2D_general), pointer :: map
    sll_int32, intent(in) :: i
    sll_int32, intent(in) :: j
    sll_int32 :: num_pts_1
    sll_int32 :: num_pts_2
    SLL_ASSERT( associated(map) )
    num_pts_1 = map%num_pts_1
    num_pts_2 = map%num_pts_2
    SLL_ASSERT( (i .ge. 1) .and. (i .le. num_pts_1) )
    SLL_ASSERT( (j .ge. 1) .and. (j .le. num_pts_2) )
    mesh_2d_x1_node = map%x1_node(i,j)
  end function mesh_2d_x1_node

  function mesh_2d_x2_node( map, i, j )
    sll_real64            :: mesh_2d_x2_node
    type(mapped_mesh_2D_general), pointer :: map
    sll_int32, intent(in) :: i
    sll_int32, intent(in) :: j
    sll_int32 :: num_pts_1
    sll_int32 :: num_pts_2
    SLL_ASSERT( associated(map) )
    num_pts_1 = map%num_pts_1
    num_pts_2 = map%num_pts_2
    SLL_ASSERT( (i .ge. 1) .and. (i .le. num_pts_1) )
    SLL_ASSERT( (j .ge. 1) .and. (j .le. num_pts_2) )
    mesh_2d_x2_node = map%x2_node(i,j)
  end function mesh_2d_x2_node

  function mesh_2d_x1_cell( map, i, j )
    sll_real64            :: mesh_2d_x1_cell
    type(mapped_mesh_2D_general), pointer :: map
    sll_int32, intent(in) :: i
    sll_int32, intent(in) :: j
    sll_int32 :: num_pts_1
    sll_int32 :: num_pts_2
    SLL_ASSERT( associated(map) )
    num_pts_1 = map%num_pts_1
    num_pts_2 = map%num_pts_2
    SLL_ASSERT( (i .ge. 1) .and. (i .le. num_pts_1 - 1) )
    SLL_ASSERT( (j .ge. 1) .and. (j .le. num_pts_2 - 1) )
    mesh_2d_x1_cell = map%x1_cell(i,j)
  end function mesh_2d_x1_cell

  function mesh_2d_x2_cell( map, i, j )
    sll_real64            :: mesh_2d_x2_cell
    type(mapped_mesh_2D_general), pointer :: map
    sll_int32, intent(in) :: i
    sll_int32, intent(in) :: j
    sll_int32 :: num_pts_1
    sll_int32 :: num_pts_2
    SLL_ASSERT( associated(map) )
    num_pts_1 = map%num_pts_1
    num_pts_2 = map%num_pts_2
    SLL_ASSERT( (i .ge. 1) .and. (i .le. num_pts_1) )
    SLL_ASSERT( (j .ge. 1) .and. (j .le. num_pts_2) )
    mesh_2d_x2_cell = map%x2_cell(i,j)
  end function mesh_2d_x2_cell

  function mesh_2d_jacobian_cell( map, i, j )
    sll_real64              :: mesh_2d_jacobian_cell
    type(mapped_mesh_2D_general), pointer   :: map
    sll_int32, intent(in)   :: i
    sll_int32, intent(in)   :: j
    sll_int32 :: num_cells_1
    sll_int32 :: num_cells_2
    SLL_ASSERT( associated(map) )
    num_cells_1 = map%num_pts_1 - 1
    num_cells_2 = map%num_pts_2 - 1
    SLL_ASSERT( (i .ge. 1) .and. (i .le. num_cells_1) )
    SLL_ASSERT( (j .ge. 1) .and. (j .le. num_cells_2) )
    mesh_2d_jacobian_cell = map%jacobians_c(i,j)
  end function mesh_2d_jacobian_cell
#endif
end module sll_module_mapped_meshes_2d
