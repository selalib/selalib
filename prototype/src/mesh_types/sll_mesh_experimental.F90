module sll_mesh_types_experimental
  use sll_splines
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
  implicit none

  ! Enumerators used inside the module and that are also available for 
  ! clients.
  enum, bind(C)
     enumerator :: NODE_CENTERED_MESH = 0, CELL_CENTERED_MESH = 1
  end enum

  enum, bind(C)
     enumerator :: ANALYTIC_MAP = 0, NUMERIC_MAP = 1
  end enum

  enum, bind(C)
     enumerator :: PERIODIC_MESH_BC = 0, HERMITE_MESH_BC = 1
  end enum

  ! Interface to represent the basic signature of all the mappings used
  ! in the 2D case. Transformations should be of the form 
  !                  x1 = x1(eta1, eta2)
  !                  x2 = x2(eta1, eta2)
  abstract interface
     function two_arg_scalar_function( eta1, eta2 )
       use sll_working_precision
       sll_real64 :: two_arg_scalar_function
       sll_real64, intent(in) :: eta1
       sll_real64, intent(in) :: eta2
     end function two_arg_scalar_function
  end interface

  ! Here we try to represent the Jacobian matrix an actual 2D array of
  ! functions. But since fortran does not allow arrays of pointers, here
  ! we define a special type that can be used as an array element.
  type jacobian_matrix_element
     procedure(two_arg_scalar_function), pointer, nopass :: f
  end type jacobian_matrix_element

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
  !                   [   partial x1        partial x1    ]
  !                   [ ---------------    -------------- ]
  !                   [   partial eta1      partial eta2  ]
  !    J(eta1,eta2) = [                                   ]
  !                   [   partial x2        partial x2    ]
  !                   [ ---------------    -------------- ]
  !                   [   partial eta1      partial eta2  ]
  !
  ! Which for convenience, can have its determinant pre-evaluated at a 
  ! collection of locations. The implementation below should provide this
  ! information in the 'jacobians' array.
  type mapping_2D   
     sll_int32  :: map_type        ! through functions or through data
     sll_int32  :: num_pts_1
     sll_int32  :: num_pts_2
     sll_real64, dimension(:,:), pointer :: x1   ! x1 = x1(eta1,eta2)
     sll_real64, dimension(:,:), pointer :: x2   ! x2 = x2(eta1,eta2)
     procedure(two_arg_scalar_function), pointer, nopass    :: x1_func
     procedure(two_arg_scalar_function), pointer, nopass    :: x2_func
     type(jacobian_matrix_element), dimension(:,:), pointer :: j_matrix
     sll_real64, dimension(:,:), pointer                    :: jacobians
  end type mapping_2D


  ! ---------------------------------------------------------------------
  !
  !   TRANSFORMED MESHES: i.e. a mesh + coordinate transformation
  !
  ! ---------------------------------------------------------------------
  type mesh_2D_scalar
     sll_int32  :: boundary_type_1
     sll_int32  :: boundary_type_2
     sll_int32  :: mode            ! data defined on nodes or center of cells
     sll_real64, dimension(:,:), pointer :: data 
     type(sll_spline_2D), pointer        :: u_spline_2D  ! u is for uniform
     type(mapping_2D), pointer           :: map
  end type mesh_2D_scalar

  interface delete
     module procedure delete_mesh_2D_scalar, delete_mapping_2D
  end interface

contains

  ! -------------------------------------------------------------------------
  !
  !         FUNCTIONS AND SUBROUTINES FOR THE MAPPING_2D TYPE.
  !
  ! -------------------------------------------------------------------------

  ! new_mapping_2D() only allocates the memory for the object itself. The
  ! initialization routines will allocate the memory of hte internal arrays.
  function new_mapping_2D( map_type )
    type(mapping_2D), pointer  :: new_mapping_2D
    sll_int32, intent(in)      :: map_type
! por aqui
    sll_real64, dimension(:,:), pointer, optional :: x1
    sll_real64, dimension(:,:), pointer, optional :: x2
    sll_int32                  :: ierr
    logical                    :: xij_given

    if( present(x1) .and. present(x2) ) then
       xij_given = .true.
    else
       xij_given = .false.
    end if
    ! Check some argument consistency
    if( (      present(x1) .and. (.not. present(x2))) .or. &
        (.not. present(x1) .and. (      present(x2))) ) then
       print *, 'ERROR, new_mapping_2D: both or none of the x1 and x2 ', &
            'arrays must be specified when calling this function.'
       STOP
    end if

    SLL_ALLOCATE( new_mapping_2D, ierr)
    new_mapping_2D%map_type  = map_type
    new_mapping_2D%num_pts_1 = npts1
    new_mapping_2D%num_pts_2 = npts2

    select case (map_type)
       case (ANALYTIC_MAP)
          if( xij_given ) then
             print *, 'ERROR, new_mapping_2D: x1 and x2 arrays should not ', &
                  'be specified in the case of an ANALYTIC_MAP.'
             STOP
          end if
          new_mapping_2D%x1 => null()
          new_mapping_2D%x2 => null()
       case(NUMERIC_MAP)
          if( .not. xij_given ) then
             ! If x1 and x2 were not passed, then we only allocate the
             ! memory. The user is responsible for initializing these
             ! arrays later.
             SLL_ALLOCATE(new_mapping_2D%x1(npts1,npts2), ierr)
             SLL_ALLOCATE(new_mapping_2D%x2(npts1,npts2), ierr)
             new_mapping_2D%x1_func => null()
             new_mapping_2D%x2_func => null()
          end if
    end select
    SLL_ALLOCATE(new_mapping_2D%j_matrix(2,2), ierr)
    SLL_ALLOCATE(new_mapping_2D%jacobians(npts1,npts2), ierr)
  end function new_mapping_2D

  ! Convenience function to compute jacobians of a 2D transformation.
  function jacobian_2D( map, eta1, eta2 )
    sll_real64                :: jacobian_2D
    type(mapping_2D), pointer :: map
    sll_real64                :: eta1
    sll_real64                :: eta2
    sll_real64                :: j11
    sll_real64                :: j12
    sll_real64                :: j21
    sll_real64                :: j22
    j11 = (map%j_matrix(1,1)%f( eta1, eta2 ))
    j12 = (map%j_matrix(1,2)%f( eta1, eta2 ))
    j21 = (map%j_matrix(2,1)%f( eta1, eta2 ))
    j22 = (map%j_matrix(2,2)%f( eta1, eta2 ))
    jacobian_2D = j11*j22 - j12*j21
  end function jacobian_2D

  ! initialize_mapping_2D() fills out the information that was not 
  ! initialized by the new_mapping_2D, with the exception of the
  ! x1 and x2 arrays, which may need to be filled out "by hand" if the
  ! user did not provide those arrays when calling new_mapping_2D.
  subroutine initialize_mapping_2D( &
    npts1,    &
    npts2,    &
    map,      &
    j11_func, &
    j12_func, &
    j21_func, &
    j22_func, &
    x1_array, &
    x2_array, &
    x1_func,  &
    x2_func )

    sll_int32, intent(in)      :: npts1
    sll_int32, intent(in)      :: npts2
    type(mapping_2D), pointer :: map
    procedure(two_arg_scalar_function)           :: j11_func
    procedure(two_arg_scalar_function)           :: j12_func
    procedure(two_arg_scalar_function)           :: j21_func
    procedure(two_arg_scalar_function)           :: j22_func
    procedure(two_arg_scalar_function), optional :: x1_func
    procedure(two_arg_scalar_function), optional :: x2_func
    sll_int32  :: map_type
    sll_int32  :: npts_1
    sll_int32  :: npts_2
    sll_real64 :: delta_1
    sll_real64 :: delta_2
    sll_real64 :: eta_1
    sll_real64 :: eta_2
    sll_int32  :: i
    sll_int32  :: j

    if(.not. associated(map)) then
       print *, 'ERROR, initialize_mapping_2D(): passed map pointer was not ',&
            'associated.'
    end if
    map_type = map%map_type

    ! Check arguments and their consistency with the declared type of map.
    select case (map_type)
       case(ANALYTIC_MAP)
          if( (.not. present(x1_func)) .or. (.not. present(x2_func)) ) then
             print *, 'ERROR: initialize_mapping_2D(): passed map ', &
                  'is an analytic map. Initialization requires analytic', &
                  'expressions for the x1_func and x2_func arguments.'
             STOP
          end if
       case(NUMERIC_MAP)
          if( present(x1_func) .or. present(x2_func) ) then
             print *, 'ERROR: initialize_mapping_2D(): passed map ', &
                  'is a numeric map. The x1_func and x2_func arguments ', &
                  'should not be specified.'
             STOP
          end if
    end select
    ! Fill the jacobian matrix
    map%j_matrix(1,1)%f => j11_func
    map%j_matrix(1,2)%f => j12_func
    map%j_matrix(2,1)%f => j21_func
    map%j_matrix(2,2)%f => j22_func

    ! Fill out the jacobians for the points in the uniform mesh.
    npts_1 = map%num_pts_1
    npts_2 = map%num_pts_2
    delta_1 = 1.0_f64/(npts_1 - 1)
    delta_2 = 1.0_f64/(npts_2 - 1)
    do j=1, npts_2
       do i=1, npts_1
          eta_1 = real(i,f64)*delta_1
          eta_2 = real(j,f64)*delta_2
          map%jacobians(i,j) = jacobian_2D( map, eta_1, eta_2 )
       end do
    end do
  end subroutine initialize_mapping_2D

  subroutine delete_mapping_2D( map )
    type(mapping_2D), pointer :: map
    sll_int32                 :: ierr
    if( associated(map%x1) ) then
       SLL_DEALLOCATE( map%x1, ierr )
    end if
    if( associated(map%x2) ) then
       SLL_DEALLOCATE( map%x2, ierr )
    end if
    SLL_DEALLOCATE( map%j_matrix, ierr )
    SLL_DEALLOCATE( map%jacobians, ierr )
    SLL_DEALLOCATE( map, ierr )
  end subroutine delete_mapping_2D


  ! Access functions for the mapping. These are an overkill and can be 
  ! changed by a macro, but for now, they are at least safer.
  function get_map2d_x1( map, i, j )
    sll_real64                :: get_map2d_x1
    type(mapping_2D), pointer :: map
    sll_int32, intent(in)     :: i
    sll_int32, intent(in)     :: j
    sll_int32 :: num_pts_1
    sll_int32 :: num_pts_2
    SLL_ASSERT( associated(map) )
    num_pts_1 = map%num_pts_1
    num_pts_2 = map%num_pts_2
    SLL_ASSERT( (i .ge. 1) .and. (i .le. num_pts_1) )
    SLL_ASSERT( (j .ge. 1) .and. (j .le. num_pts_2) )
    get_map2d_x1 = map%x1(i,j)
  end function get_map2d_x1

  function get_map2d_x2( map, i, j )
    sll_real64                :: get_map2d_x2
    type(mapping_2D), pointer :: map
    sll_int32, intent(in)     :: i
    sll_int32, intent(in)     :: j
    sll_int32 :: num_pts_1
    sll_int32 :: num_pts_2
    SLL_ASSERT( associated(map) )
    num_pts_1 = map%num_pts_1
    num_pts_2 = map%num_pts_2
    SLL_ASSERT( (i .ge. 1) .and. (i .le. num_pts_1) )
    SLL_ASSERT( (j .ge. 1) .and. (j .le. num_pts_2) )
    get_map2d_x2 = map%x2(i,j)
  end function get_map2d_x2

  subroutine set_map2d_x1( map, i, j, val )
    type(mapping_2D), pointer :: map
    sll_int32, intent(in)     :: i
    sll_int32, intent(in)     :: j
    sll_real64, intent(in)    :: val
    sll_int32 :: num_pts_1
    sll_int32 :: num_pts_2
    SLL_ASSERT( associated(map) )
    num_pts_1 = map%num_pts_1
    num_pts_2 = map%num_pts_2
    SLL_ASSERT( (i .ge. 1) .and. (i .le. num_pts_1) )
    SLL_ASSERT( (j .ge. 1) .and. (j .le. num_pts_2) )
    map%x1(i,j) = val
  end subroutine set_map2d_x1

  subroutine set_map2d_x2( map, i, j, val )
    type(mapping_2D), pointer :: map
    sll_int32, intent(in)     :: i
    sll_int32, intent(in)     :: j
    sll_real64, intent(in)    :: val
    sll_int32 :: num_pts_1
    sll_int32 :: num_pts_2
    SLL_ASSERT( associated(map) )
    num_pts_1 = map%num_pts_1
    num_pts_2 = map%num_pts_2
    SLL_ASSERT( (i .ge. 1) .and. (i .le. num_pts_1) )
    SLL_ASSERT( (j .ge. 1) .and. (j .le. num_pts_2) )
    map%x2(i,j) = val
  end subroutine  set_map2d_x2

  subroutine set_map_j_matrix_elem( map, i, j, func )
    type(mapping_2D), pointer          :: map
    procedure(two_arg_scalar_function) :: func
    sll_int32, intent(in)              :: i
    sll_int32, intent(in)              :: j
    SLL_ASSERT( associated(map) )
    SLL_ASSERT( (i .ge. 1) .and. (i .le. 2) )
    SLL_ASSERT( (j .ge. 1) .and. (j .le. 2) )
    map%j_matrix(i,j)%f => func
  end subroutine set_map_j_matrix_elem

  function get_map2d_jacobian( map, i, j)
    sll_real64                  :: get_map2d_jacobian
    type(mapping_2D), pointer   :: map
    sll_int32, intent(in)       :: i
    sll_int32, intent(in)       :: j
    get_map2d_jacobian = map%jacobians(i,j)
  end function get_map2d_jacobian


  ! -------------------------------------------------------------------------
  !
  !     FUNCTIONS AND SUBROUTINES FOR THE TRANSFORMED MESH_2D TYPE.
  !
  ! -------------------------------------------------------------------------

  ! The mesh_2D_scalar is basically data interpreted with the aid of a 
  ! 2D coordinate map. Additionally, we need information on:
  ! mode: whether the mesh data is to be taken to represent values on nodes
  ! or centers of cells.
  ! bc_1 and bc_2: the boundary types to use when building the interpolants,
  ! usually cubic splines.
  function new_mesh_2D_scalar( &
    mode, &
    bc_1, &
    bc_2, &
    mapping, &
    eta1_min_slope, &
    eta1_max_slope, &
    eta2_min_slope, &
    eta2_max_slope )

    type(mesh_2D_scalar), pointer                :: new_mesh_2D_scalar
    sll_int32, intent(in)                        :: mode
    sll_int32, intent(in)                        :: bc_1
    sll_int32, intent(in)                        :: bc_2
    type(mapping_2D), pointer                    :: mapping
    sll_real64, intent(in), optional             :: eta1_min_slope
    sll_real64, intent(in), optional             :: eta1_max_slope
    sll_real64, intent(in), optional             :: eta2_min_slope
    sll_real64, intent(in), optional             :: eta2_max_slope

    sll_int32                     :: num_pts_1
    sll_int32                     :: num_pts_2
    sll_int32                     :: ierr
    sll_int32                     :: spline_bc_1
    sll_int32                     :: spline_bc_2

    if( .not. associated(mapping) ) then
       print *, 'ERROR, new_mesh_2D_scalar: mapping given as argument is ', &
            'not associated. '
    end if
    SLL_ALLOCATE(new_mesh_2D_scalar, ierr)
    ! The information from the mapping is extracted to build the data array.
    num_pts_1 = mapping%num_pts_1
    num_pts_2 = mapping%num_pts_2

    ! Select what kind of boundary conditions to pass along to the
    ! splines.
    select case ( bc_1 )
       case ( PERIODIC_MESH_BC )
          spline_bc_1 = PERIODIC_SPLINE
       case ( HERMITE_MESH_BC )
          spline_bc_1 = HERMITE_SPLINE
       case default
          print *, 'ERROR, new_mesh_2D_scalar(): unrecognized boundary ', &
               'condition in direction eta_1'
          STOP
    end select

    select case ( bc_2 )
       case ( PERIODIC_MESH_BC )
          spline_bc_2 = PERIODIC_SPLINE
       case ( HERMITE_MESH_BC )
          spline_bc_2 = HERMITE_SPLINE
       case default
          print *, 'ERROR, new_mesh_2D_scalar(): unrecognized boundary ', &
               'condition in direction eta_2'
          STOP
    end select

    ! Allocate the memory for the different fields
    select case (mode)
       case (NODE_CENTERED_MESH)
          new_mesh_2D_scalar%mode = NODE_CENTERED_MESH
          SLL_ALLOCATE(new_mesh_2D_scalar%data(num_pts_1, num_pts_2), ierr)
          new_mesh_2D_scalar%u_spline_2D => new_spline_2D( num_pts_1, &
                                                           num_pts_2, &
                                                           0.0_f64, &
                                                           1.0_f64, &
                                                           0.0_f64, &
                                                           1.0_f64, &
                                                           spline_bc_1, &
                                                           spline_bc_2, &
                                                           eta1_min_slope, &
                                                           eta1_max_slope, &
                                                           eta2_min_slope, &
                                                           eta2_max_slope )
       case (CELL_CENTERED_MESH)
          print *, 'ERROR, new_mesh_2D_scalar: the CELL_CENTERED_MESH ', &
               'methods have not been implemented.'
          STOP
       case default
          print *, "ERROR: unrecognized mode in new_mesh_2D_scalar()"
          STOP
    end select
  end function new_mesh_2D_scalar

  subroutine delete_mesh_2D_scalar( mesh )
    type(mesh_2d_scalar), pointer :: mesh
    sll_int32                     :: ierr
    if( .not. associated(mesh) ) then
       print *, "ERROR, delete_mesh_2D_scalar(): mesh was not associated."
       STOP
    end if
    SLL_DEALLOCATE(mesh%data, ierr)
    call delete(mesh%u_spline_2D)
    mesh%map => null()
    SLL_DEALLOCATE(mesh, ierr)
  end subroutine delete_mesh_2D_scalar

  ! Add a copy constructor here...

  ! compute_mesh2D_interpolants() updates the interpolation information
  ! (like spline coefficients, if a cubic spline interpolation is used) that
  ! the mesh use.
  subroutine compute_mesh2D_interpolants( mesh )
    type(mesh_2D_scalar), pointer :: mesh
    if( .not. associated(mesh) ) then
       print *, 'ERROR, compute_mesh2D_interpolants: mesh pointer argument ',&
            'was not associated.'
       STOP
    end if
    call compute_spline_2D( mesh%data, mesh%u_spline_2d )
  end subroutine compute_mesh2D_interpolants


  ! The following name is ugly and even bad, we need something specific
  ! to the type of mesh and also concise... so maybe a convention is needed.
  function get_m2ds_node( mesh, i, j )
    sll_real64                    :: get_m2ds_node
    type(mesh_2D_scalar), pointer :: mesh
    sll_int32, intent(in)         :: i
    sll_int32, intent(in)         :: j
    SLL_ASSERT( associated(mesh) )
    get_m2ds_node = mesh%data(i,j)
  end function get_m2ds_node

  ! Add here get_m2ds_val(mesh, x1, x2), this requires the nonuniform splines.

  subroutine set_m2ds_node( mesh, i, j, val )
    type(mesh_2D_scalar), pointer :: mesh
    sll_int32, intent(in)         :: i
    sll_int32, intent(in)         :: j
    sll_real64, intent(in)        :: val
    SLL_ASSERT( associated(mesh) )
    mesh%data(i,j) = val
  end subroutine set_m2ds_node

  ! Add a function to compute the determinant at a point.

  ! obtain the interpolated value of the mesh, as a function of the
  ! coordinates in the uniform coordinate system.
  function get_mesh2D_value( mesh, eta1, eta2 )
    sll_real64                    :: get_mesh2D_value
    type(mesh_2D_scalar), pointer :: mesh
    sll_real64                    :: eta1
    sll_real64                    :: eta2
    get_mesh2D_value = interpolate_value_2D( eta1, eta2, mesh%u_spline_2D )
  end function get_mesh2D_value

  function mesh_2D_data( mesh )
    sll_real64, dimension(:,:), pointer :: mesh_2D_data
    type(mesh_2D_scalar), pointer       :: mesh
    mesh_2D_data => mesh%data
  end function mesh_2D_data

end module sll_mesh_types_experimental
