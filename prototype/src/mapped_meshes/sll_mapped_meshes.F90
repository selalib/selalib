module sll_mapped_meshes
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
  use sll_splines
  use sll_mapped_mesh_base
  use sll_interpolators_base
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
     procedure(two_arg_scalar_function), pointer, nopass    :: x1_func  ! user
     procedure(two_arg_scalar_function), pointer, nopass    :: x2_func  ! user
     type(jacobian_matrix_element), dimension(:,:), pointer :: j_matrix
     procedure(two_arg_message_passing_func_analyt), pointer, pass :: &
          jacobian_func
     procedure(j_matrix_f_nopass), pointer, nopass          :: jacobian_matrix
   contains
     procedure, pass(mesh) :: initialize => initialize_mesh_2d_analytic
     procedure, pass(mesh) :: x1_at_node => x1_node_analytic
     procedure, pass(mesh) :: x2_at_node => x2_node_analytic
     procedure, pass(mesh) :: x1         => x1_analytic
     procedure, pass(mesh) :: x2         => x2_analytic
     procedure, pass(mesh) :: jacobian   => jacobian_2d_analytic
  end type sll_mapped_mesh_2d_analytic
  
  type, extends(sll_mapped_mesh_2d_base)::sll_mapped_mesh_2d_discrete
     procedure(two_arg_scalar_function), pointer, nopass    :: x1_func
     procedure(two_arg_scalar_function), pointer, nopass    :: x2_func
     type(jacobian_matrix_element), dimension(:,:), pointer :: j_matrix
     class(interpolator_2d_base), pointer                   :: x1_interp
     class(interpolator_2d_base), pointer                   :: x2_interp
     procedure(two_arg_message_passing_func_discr),pointer,pass :: jacobian_func

   contains
     procedure, pass(mesh) :: initialize => initialize_mesh_2d_discrete
     procedure, pass(mesh) :: x1_at_node => x1_node_discrete
     procedure, pass(mesh) :: x2_at_node => x2_node_discrete
     procedure, pass(mesh) :: x1         => x1_discrete
     procedure, pass(mesh) :: x2         => x2_discrete
     procedure, pass(mesh) :: jacobian   => jacobian_2d_discrete
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
     procedure(two_arg_scalar_function), pointer, nopass :: f
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
  !
  ! This interface is ending up very awkward because of the large amount of
  ! optional parameters that it takes, much of this in account of the 
  ! splines that it initializes, which take plenty of optional parameters
  ! themselves. This is not desirable and should be reassessed critically.
  !
  ! We should offer the possibility to pass the jacobian function directly.
  ! What routines like these really need are self-checking and consistency
  ! functions for safety.

  subroutine initialize_mesh_2d_analytic( &
    mesh,           &
    npts1,          &
    npts2,          &
    x1_func,        &
    x2_func,        &
    j11_func,       &
    j12_func,       &
    j21_func,       &
    j22_func )

    class(sll_mapped_mesh_2d_analytic), intent(inout) :: mesh
    sll_int32, intent(in)                         :: npts1
    sll_int32, intent(in)                         :: npts2
    procedure(two_arg_scalar_function)            :: x1_func
    procedure(two_arg_scalar_function)            :: x2_func
    procedure(two_arg_scalar_function)            :: j11_func
    procedure(two_arg_scalar_function)            :: j12_func
    procedure(two_arg_scalar_function)            :: j21_func
    procedure(two_arg_scalar_function)            :: j22_func

    sll_real64 :: delta_1  ! cell spacing in eta1 
    sll_real64 :: delta_2  ! cell spacing in eta2 
    sll_real64 :: eta_1
    sll_real64 :: eta_2
    sll_real64 :: jacobian_val
    sll_int32  :: i
    sll_int32  :: j
    sll_int32  :: ierr

    mesh%num_pts1 = npts1
    mesh%num_pts2 = npts2
    delta_1       = 1.0_f64/(npts1 - 1)
    delta_2       = 1.0_f64/(npts2 - 1)
    mesh%delta1   = delta_1
    mesh%delta2   = delta_2

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
    mesh%x1_func => x2_func

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
          mesh%x1_node(i+1,j+1) = x1_func(eta_1, eta_2)
          mesh%x2_node(i+1,j+1) = x2_func(eta_1, eta_2)
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
          mesh%x1_cell(i+1,j+1)     = x1_func(eta_1, eta_2)
          mesh%x2_cell(i+1,j+1)     = x2_func(eta_1, eta_2)
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
    j11 = (mesh%j_matrix(1,1)%f( eta1, eta2 ))
    j12 = (mesh%j_matrix(1,2)%f( eta1, eta2 ))
    j21 = (mesh%j_matrix(2,1)%f( eta1, eta2 ))
    j22 = (mesh%j_matrix(2,2)%f( eta1, eta2 ))
    ! For debugging:
    !    print *, 'jacobian_2d_analytic: '
    !    print *, j11, j12
    !    print *, j21, j22
    val = j11*j22 - j12*j21
  end function jacobian_2d_analytic

  function x1_analytic( mesh, eta1, eta2 ) result(val)
    sll_real64                        :: val
    class(sll_mapped_mesh_2d_analytic) :: mesh
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    val = mesh%x1_func(eta1, eta2)
  end function x1_analytic

  function x2_analytic( mesh, eta1, eta2 ) result(val)
    sll_real64                        :: val
    class(sll_mapped_mesh_2d_analytic) :: mesh
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    val = mesh%x2_func(eta1, eta2)
  end function x2_analytic

  function x1_node_analytic( mesh, i, j ) result(val)
    sll_real64             :: val
    class(sll_mapped_mesh_2d_analytic) :: mesh
    sll_int32, intent(in) :: i
    sll_int32, intent(in) :: j
    val = mesh%x1_node(i,j)
  end function x1_node_analytic

  function x2_node_analytic( mesh, i, j ) result(val)
    sll_real64             :: val
    class(sll_mapped_mesh_2d_analytic) :: mesh
    sll_int32, intent(in) :: i
    sll_int32, intent(in) :: j
    val = mesh%x2_node(i,j)
  end function x2_node_analytic

  !**************************************************************************
  !
  !               Functions for the discrete general map
  !
  !**************************************************************************

  function x1_node_discrete( mesh, i, j ) result(val)
    sll_real64             :: val
    class(sll_mapped_mesh_2d_discrete) :: mesh
    sll_int32, intent(in) :: i
    sll_int32, intent(in) :: j
    val = mesh%x1_node(i,j)
  end function x1_node_discrete

  function x2_node_discrete( mesh, i, j ) result(val)
    sll_real64             :: val
    class(sll_mapped_mesh_2d_discrete) :: mesh
    sll_int32, intent(in) :: i
    sll_int32, intent(in) :: j
    val = mesh%x2_node(i,j)
  end function x2_node_discrete

  function x1_discrete( mesh, eta1, eta2 ) result(val)
    sll_real64                        :: val
    class(sll_mapped_mesh_2d_discrete) :: mesh
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    val = mesh%x1_interp%interpolate_value(eta1, eta2)
  end function x1_discrete

  function x2_discrete( mesh, eta1, eta2 ) result(val)
    sll_real64                        :: val
    class(sll_mapped_mesh_2d_discrete) :: mesh
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    val = mesh%x2_interp%interpolate_value(eta1, eta2)
  end function x2_discrete

  function jacobian_2d_discrete( mesh, eta1, eta2 ) result(jac)
    sll_real64             :: jac
    class(sll_mapped_mesh_2d_discrete) :: mesh
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


  subroutine initialize_mesh_2d_discrete( &
    mesh,            &
    npts1,          &
    npts2,          &
    x1_node,        &
    x2_node,        &
    x1_interpolator, &
    x2_interpolator, &
    jacobians_node, &
    jacobians_n_interpolator, &
    x1_cell, &
    x2_cell, &
    jacobians_cell )

    class(sll_mapped_mesh_2d_discrete)    :: mesh
    sll_int32, intent(in)                :: npts1
    sll_int32, intent(in)                :: npts2
    sll_real64, dimension(:,:)           :: x1_node
    sll_real64, dimension(:,:)           :: x2_node
    class(interpolator_2d_base), target  :: x1_interpolator
    class(interpolator_2d_base), target  :: x2_interpolator
    class(interpolator_2d_base), target  :: jacobians_n_interpolator  
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

    mesh%num_pts1  = npts1
    mesh%num_pts2  = npts2
    mesh%delta1    = 1.0_f64/(npts1 - 1)
    mesh%delta2    = 1.0_f64/(npts2 - 1)
    mesh%x1_interp => x1_interpolator
    mesh%x1_interp => x2_interpolator

    ! Allocate the arrays for precomputed jacobians.
    SLL_ALLOCATE(mesh%jacobians_n(npts1,npts2), ierr)
    SLL_ALLOCATE(mesh%jacobians_c(npts1-1, npts2-1), ierr)

    ! Allocation for x1 and x2 at nodes
    SLL_ALLOCATE(mesh%x1_node(npts1,npts2), ierr)
    SLL_ALLOCATE(mesh%x2_node(npts1,npts2), ierr)

    ! Start filling out the fields and allocating the object's memory.
    SLL_ALLOCATE(mesh%x1_node(npts1,npts2), ierr)
    SLL_ALLOCATE(mesh%x2_node(npts1,npts2), ierr)

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
          eta_2 = real(j,f64)*mesh%delta2          
          do i=0, npts1 - 1
             eta_1 = real(i,f64)*mesh%delta1
             ! FIX THIS PART!!!
             jacobian_val = mesh%jacobian_func(eta_1,eta_2)
             mesh%jacobians_n(i+1,j+1) = jacobian_val
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
    end if
  end subroutine initialize_mesh_2d_discrete

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

  function mesh_2d_jacobian_node( map, i, j )
    sll_real64              :: mesh_2d_jacobian_node
    type(mapped_mesh_2D_general), pointer   :: map
    sll_int32, intent(in)   :: i
    sll_int32, intent(in)   :: j
    sll_int32 :: num_pts_1
    sll_int32 :: num_pts_2
    SLL_ASSERT( associated(map) )
    num_pts_1 = map%num_pts_1
    num_pts_2 = map%num_pts_2
    SLL_ASSERT( (i .ge. 1) .and. (i .le. num_pts_1) )
    SLL_ASSERT( (j .ge. 1) .and. (j .le. num_pts_2) )
    mesh_2d_jacobian_node = map%jacobians_n(i,j)
  end function mesh_2d_jacobian_node

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
end module sll_mapped_meshes
