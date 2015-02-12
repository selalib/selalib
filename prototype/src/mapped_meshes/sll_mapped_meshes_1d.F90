module sll_module_mapped_meshes_1d
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
use sll_cubic_splines
  use sll_module_interpolators_1d_base
  use sll_module_mapped_meshes_1d_base
  use sll_cubic_spline_interpolator_1d
  implicit none

  ! A 1D map is specified by the coordinate transformation from eta1 to x1:
  !
  !                x1 = x1(eta1)
  !
  ! Where eta1 is in the interval [0,1]. The same transformation can be 
  ! specified by a set of transformed points x1(i), as a 1D array.
  !
  ! The transformation is also represented by the jacobian :
  !
  !                               partial x1(eta1)  
  !          jacobian(eta1) =   --------------------
  !                               partial   eta1        
  !
  ! For convenience, we can have this pre-evaluated at a collection of 
  ! locations.

  type, extends(sll_mapped_mesh_1d_base) :: sll_mapped_mesh_1d_analytic
     sll_real64, dimension(:), pointer :: x1_node   ! x1(i) 
     sll_real64, dimension(:), pointer :: x1_cell
     sll_real64, dimension(:), pointer :: jacobians_n
     sll_real64, dimension(:), pointer :: jacobians_c
     sll_real64, dimension(:), pointer :: params
     procedure(one_arg_scalar_function_nopass), pointer, nopass :: x1_func !user
     procedure(one_arg_scalar_function_nopass), pointer, nopass :: jacobian_func
   contains
     procedure, pass(mesh) :: initialize => initialize_mesh_1d_analytic
     procedure, pass(mesh) :: x1_at_node => x1_node_analytic_1d
     procedure, pass(mesh) :: x1_at_cell => x1_cell_analytic_1d
     procedure, pass(mesh) :: jacobian_at_node => mesh_1d_jacobian_node_analytic
     procedure, pass(mesh) :: jacobian_at_cell => mesh_1d_jacobian_cell_analytic
     procedure, pass(mesh) :: x1         => x1_analytic_1d
     procedure, pass(mesh) :: jacobian   => jacobian_1d_analytic
  end type sll_mapped_mesh_1d_analytic

  type, extends(sll_mapped_mesh_1d_base)::sll_mapped_mesh_1d_discrete
     sll_real64, dimension(:), pointer                    :: x1_node   ! x1(i)
     sll_real64, dimension(:), pointer                    :: x1_cell
     sll_real64, dimension(:), pointer                    :: jacobians_n
     sll_real64, dimension(:), pointer                    :: jacobians_c
     class(sll_interpolator_1d_base), pointer             :: x1_interp
!     procedure(one_arg_scalar_function), pointer, nopass  :: x1_func
    ! Samuel : this functions seems never initialize, see jacobian
    ! procedure(one_arg_message_passing_func_discr),pointer,pass :: jacobian_func
   contains
     procedure, pass(mesh) :: initialize => initialize_mesh_1d_discrete
     procedure, pass(mesh) :: x1_at_node => x1_node_discrete_1d
     procedure, pass(mesh) :: x1_at_cell => x1_cell_discrete_1d
     procedure, pass(mesh) :: jacobian_at_node => mesh_1d_jacobian_node_discrete
     procedure, pass(mesh) :: jacobian_at_cell => mesh_1d_jacobian_cell_discrete
     procedure, pass(mesh) :: x1         => x1_discrete_1d
     procedure, pass(mesh) :: jacobian   => jacobian_1d_discrete
  end type sll_mapped_mesh_1d_discrete

  abstract interface
     function one_arg_message_passing_func_analyt( map, eta1 )
       use sll_working_precision
       import :: sll_mapped_mesh_1d_analytic
       sll_real64 :: one_arg_message_passing_func_analyt
       class(sll_mapped_mesh_1d_analytic)  :: map
       sll_real64, intent(in)              :: eta1
     end function one_arg_message_passing_func_analyt
  end interface


  abstract interface
     function one_arg_message_passing_func_discr( map, eta1 )
       use sll_working_precision
       import     :: sll_mapped_mesh_1d_discrete
       sll_real64                          :: one_arg_message_passing_func_discr
       class(sll_mapped_mesh_1d_discrete)  :: map
       sll_real64, intent(in)              :: eta1
     end function one_arg_message_passing_func_discr
  end interface

contains

 !**************************************************************************
  !
  !               Functions for the analytic general map
  !
  !**************************************************************************

  ! initialize_mapped_mesh_1d_general() allocates all the memory needed by 
  ! the 1D map. 
  !

  subroutine initialize_mesh_1d_analytic( &
    mesh,           &
    label,          &
    npts,           &
    x1_func,        &
    jacobian_func,  &
    params )

    class(sll_mapped_mesh_1d_analytic), intent(inout) :: mesh
    character(len=*), intent(in)                  :: label
    sll_int32, intent(in)                         :: npts

    procedure(one_arg_scalar_function_nopass)     :: x1_func
    procedure(one_arg_scalar_function_nopass)     :: jacobian_func
    sll_real64, dimension(:), intent(in) :: params
    sll_real64 :: delta_1  ! cell spacing in eta1 
    sll_real64 :: eta_1
    sll_real64 :: jacobian_val
    sll_int32  :: i
    sll_int32  :: ierr
    sll_int32  :: npts1

    npts1 = npts
    mesh%label      = trim(label)
    mesh%nc_eta1    = npts1-1
    delta_1         = 1.0_f64/(npts1 - 1)
    mesh%delta_eta1 = delta_1

    SLL_ALLOCATE(mesh%params(size(params)),ierr)
    mesh%params(:) = params(:)
    ! Allocate the arrays for precomputed jacobians.
    SLL_ALLOCATE(mesh%jacobians_n(npts1), ierr)
    SLL_ALLOCATE(mesh%jacobians_c(npts1-1), ierr)

    ! Allocation for x1 at nodes, needed regardless of the type of map
    SLL_ALLOCATE(mesh%x1_node(npts1), ierr)

    ! Start filling out the fields and allocating the object's memory.
    SLL_ALLOCATE(mesh%x1_cell(npts1-1), ierr)

    ! Assign the transformation functions
    mesh%x1_func       => x1_func
    mesh%jacobian_func => jacobian_func
    
    ! Allocate the arrays for precomputed jacobians.
    SLL_ALLOCATE(mesh%jacobians_n(npts1), ierr)
    SLL_ALLOCATE(mesh%jacobians_c(npts1-1), ierr)
    
    ! Fill the values of the transformation and the jacobians at the nodes
    do i=0, npts1 - 1
       eta_1 = real(i,f64)*delta_1
       mesh%x1_node(i+1) = x1_func(eta_1,params)
       ! for some compiler reason, the following intermediate 
       ! variable is required, else the jacobians_n array will not
       ! be filled out properly.
       jacobian_val          = jacobian_func(eta_1,params)
       mesh%jacobians_n(i+1) = jacobian_val
    end do
    
    ! Fill the values at the mid-point of the cells
    do i=0, npts1 - 2
       eta_1 = delta_1*(real(i,f64) + 0.5_f64)
       mesh%x1_cell(i+1)     = x1_func(eta_1,params)
       mesh%jacobians_c(i+1) = jacobian_func(eta_1,params)
    end do
  end subroutine initialize_mesh_1d_analytic

  function x1_node_analytic_1d( mesh, i ) result(val)
    class(sll_mapped_mesh_1d_analytic) :: mesh
    sll_real64             :: val
    sll_int32, intent(in) :: i
    val = mesh%x1_node(i)
  end function x1_node_analytic_1d

  function x1_cell_analytic_1d( mesh, i ) result(val)
    class(sll_mapped_mesh_1d_analytic) :: mesh
    sll_real64             :: val
    sll_int32, intent(in) :: i
    val = mesh%x1_cell(i)
  end function x1_cell_analytic_1d

  function x1_analytic_1d( mesh, eta1 ) result(val)
    class(sll_mapped_mesh_1d_analytic) :: mesh
    sll_real64                         :: val
    sll_real64, intent(in) :: eta1
    val = mesh%x1_func(eta1,mesh%params)
  end function x1_analytic_1d

  function mesh_1d_jacobian_node_analytic( mesh, i ) result(val)
    class(sll_mapped_mesh_1d_analytic)   :: mesh
    sll_real64              :: val
    sll_int32, intent(in)   :: i
    sll_int32 :: num_pts_1
    num_pts_1 = mesh%nc_eta1 + 1
    SLL_ASSERT( (i .ge. 1) .and. (i .le. num_pts_1) )
    val = mesh%jacobians_n(i)
  end function mesh_1d_jacobian_node_analytic

  function mesh_1d_jacobian_cell_analytic( mesh, i ) result(val)
    class(sll_mapped_mesh_1d_analytic)   :: mesh
    sll_real64              :: val
    sll_int32, intent(in)   :: i
    sll_int32 :: num_pts_1
    num_pts_1 = mesh%nc_eta1 
    SLL_ASSERT( (i .ge. 1) .and. (i .le. num_pts_1) )
    val = mesh%jacobians_c(i)
  end function mesh_1d_jacobian_cell_analytic

  function jacobian_1d_analytic( mesh, eta1 ) result(val)
    sll_real64                         :: val
    class(sll_mapped_mesh_1d_analytic) :: mesh
    sll_real64, intent(in) :: eta1
    val = mesh%jacobian_func(eta1,mesh%params)
  end function jacobian_1d_analytic

  ! Discrete case:

  subroutine initialize_mesh_1d_discrete( &
    mesh,            &
    label,            &
    npts1,          &
    x1_node,        &
    x1_interpolator, &
    jacobians_n_interpolator, &
    jacobians_node, &
    x1_cell, &
    jacobians_cell )

    class(sll_mapped_mesh_1d_discrete)    :: mesh
    character(len=*), intent(in)          :: label
    sll_int32, intent(in)                 :: npts1
    sll_real64, dimension(:)              :: x1_node
    class(sll_interpolator_1d_base), target  :: x1_interpolator
    class(sll_interpolator_1d_base), target  :: jacobians_n_interpolator
    sll_real64, dimension(:), optional :: jacobians_node
    sll_real64, dimension(:), optional :: jacobians_cell
    sll_real64, dimension(:), optional :: x1_cell

    sll_real64 :: eta_1
    sll_real64 :: jacobian_val
    sll_int32  :: i
    sll_int32  :: ierr
    logical    :: x1c
    logical    :: jc
    logical    :: jn

    mesh%label = trim(label)
    x1c = present(x1_cell)
    jc  = present(jacobians_cell)
    jn  = present(jacobians_node)
    ! Samuel : I've commented the next line because jacobian_func isn't use
    !mesh%jacobian_func => null() !This functions is never initialize!!

    ! Check argument consistency
    ! DISCRETE_MAPs require only some of the parameters. If the mapping is
    ! defined from the nodes of the logical (eta1) mesh to the nodes
    ! of the physical mesh (x1), then the node arrays are required:
    ! jacobians_node and x1_node.
    !
    ! If the transformation is done on the points at the center of the cells
    ! then these parameters are also required: 
    ! jacobians_cell and x1_cell.
    ! node and cell values are not mutually exclusive, thus all parameters 
    ! can be provided in the DISCRETE case. It is up to the caller to make
    ! sure that the data set is consistent.

    ! 1. Check that the discrete representation of the transformation is
    !    consistent with the size of the 2D array.
    if( (size(x1_node,1) .lt. npts1) ) then
       print *, 'ERROR, initialize_mapped_mesh_1d_discrete(), ', &
            'the size of the x1_node array is ', &
            'inconsistent with the number of points declared, npts1.'
       STOP
    end if
    if( jn .eqv. .true. ) then
       if( (size(jacobians_node,1) .lt. npts1 - 1 ) ) then
          print *, 'ERROR, initialize_mapped_mesh_1d_discrete(): ', &
               'the size of the jacobians_node array is ', &
               'inconsistent with the number of points declared, npts1.'
          STOP
       end if
    end if
    if( jc .eqv. .true. ) then
       if( (size(jacobians_cell,1) .lt. npts1 - 1 ) ) then
          print *, 'ERROR, initialize_mapped_mesh_1d_discrete(): ', &
               'the size of the jacobians_cell arrays is ', &
               'inconsistent with the number of points declared, npts1.'
          STOP
       end if
    end if

    mesh%nc_eta1    = npts1-1
    mesh%delta_eta1 = 1.0_f64/(npts1 - 1)
    mesh%x1_interp  => x1_interpolator

    ! Allocate the arrays for precomputed jacobians.
    SLL_ALLOCATE(mesh%jacobians_n(npts1), ierr)
    SLL_ALLOCATE(mesh%jacobians_c(npts1-1), ierr)

    ! Allocation for x1 at nodes
    SLL_ALLOCATE(mesh%x1_node(npts1), ierr)

    ! Start filling out the fields and allocating the object's memory.
    SLL_ALLOCATE(mesh%x1_node(npts1), ierr)

    ! initialize the local arrays. Note that since the map has its
    ! own copies, it owns this information locally and will destroy
    ! this information when the object is deleted. The caller is
    ! thus responsible for deallocating the arrays that were passed as
    ! arguments.
    do i=1, npts1
       mesh%x1_node(i) = x1_node(i)
    end do
    
    ! Compute the spline coefficients
    call x1_interpolator%compute_interpolants( mesh%x1_node )

    ! The splines contain all the information to compute the
    ! jacobians everywhere; however, here we explore assigning
    ! the jacobians-at-the-nodes array with the values provided
    ! by the user if available. If there are discrepancies between
    ! the user-provided values and the predictions from the splines,
    ! then this may itself be a way to look for errors.
    !
    ! Copy the values of the jacobians at the nodes if user given:
    if( jn .eqv. .true. ) then
       do i=1, npts1
          mesh%jacobians_n(i) = jacobians_node(i)
       end do
    else
       ! Fill the jacobian values at the nodes calculated from the splines
       do i=0, npts1 - 1
          eta_1 = real(i,f64)*mesh%delta_eta1
          ! FIX THIS PART!!!
          ! Samuel : I replace the line
          !jacobian_val = mesh%jacobian_func(eta_1)
          ! by the next
          jacobian_val = mesh%jacobian(eta_1)
          mesh%jacobians_n(i+1) = jacobian_val
       end do
    end if

    ! copy the cell-based transformation arrays if available
    if( x1c .eqv. .true. ) then
       SLL_ALLOCATE(mesh%x1_cell(npts1-1), ierr)
       do i=1, npts1 - 1
          mesh%x1_cell(i) = x1_cell(i)
       end do
    end if
    ! copy the cell-based jacobians if available
    if( jc .eqv. .true. ) then
       do i=1, npts1 - 1
          mesh%jacobians_c(i) = jacobians_cell(i)
       end do
    end if
  end subroutine initialize_mesh_1d_discrete

  function x1_node_discrete_1d( mesh, i ) result(val)
    class(sll_mapped_mesh_1d_discrete) :: mesh
    sll_real64             :: val
    sll_int32, intent(in) :: i
    val = mesh%x1_node(i)
  end function x1_node_discrete_1d

  function x1_cell_discrete_1d( mesh, i ) result(val)
    class(sll_mapped_mesh_1d_discrete) :: mesh
    sll_real64             :: val
    sll_int32, intent(in) :: i
    val = mesh%x1_cell(i)
  end function x1_cell_discrete_1d

  function mesh_1d_jacobian_node_discrete( mesh, i ) result(var)
    class(sll_mapped_mesh_1d_discrete)   :: mesh
    sll_real64              :: var
    sll_int32, intent(in)   :: i
    sll_int32 :: num_pts_1
    num_pts_1 = mesh%nc_eta1 + 1
    SLL_ASSERT( (i .ge. 1) .and. (i .le. num_pts_1) )
    var = mesh%jacobians_n(i)
  end function mesh_1d_jacobian_node_discrete

  function mesh_1d_jacobian_cell_discrete( mesh, i ) result(var)
    class(sll_mapped_mesh_1d_discrete)   :: mesh
    sll_real64              :: var
    sll_int32, intent(in)   :: i
    sll_int32 :: num_pts_1
    num_pts_1 = mesh%nc_eta1
    SLL_ASSERT( (i .ge. 1) .and. (i .le. num_pts_1) )
    var = mesh%jacobians_c(i)
  end function mesh_1d_jacobian_cell_discrete

  function jacobian_1d_discrete( mesh, eta1 ) result(jac)
    class(sll_mapped_mesh_1d_discrete) :: mesh
    sll_real64             :: jac
    sll_real64, intent(in) :: eta1
    jac = mesh%x1_interp%interpolate_derivative_eta1( eta1 )
  end function jacobian_1d_discrete

  function x1_discrete_1d( mesh, eta1 ) result(val)
    class(sll_mapped_mesh_1d_discrete) :: mesh
    sll_real64                         :: val
    sll_real64, intent(in) :: eta1
    val = mesh%x1_interp%interpolate_value(eta1)
  end function x1_discrete_1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module sll_module_mapped_meshes_1d
