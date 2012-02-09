module sll_coordinate_transformation

#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
  use sll_splines
  implicit none

  ! Enumerators used inside the module and that are also available for 
  ! clients.
  enum, bind(C)
     enumerator :: NODE_CENTERED_MESH = 0, CELL_CENTERED_MESH = 1
  end enum

  enum, bind(C)
     enumerator :: ANALYTIC_MAP = 0, DISCRETE_MAP = 1
  end enum

  ! Implementation note: The following enumerator decouples the description
  ! of the boundary conditions for a coordinate transformation from the
  ! analogous description that is needed for the underlying splines. In
  ! other words: in the discrete case, to represent the x1(eta1,eta2) 
  ! and x2(eta1,eta2) transformations AS IF these were continuous 
  ! transformations, we need some underlying continuous representation, like
  ! cubic splines. The specification of the splines needs the type of 
  ! boundary conditions used, and thus, this information needs to be 
  ! passed along through the interface of the map_2D. However, it would not
  ! be good to pass the same enumerator that the spline module uses. This
  ! would expose the information about the underlying representation and
  ! would couple too strongly the dependence on the underlying splines
  ! implementation. Hence the map type needs its own enumerator which,
  ! within the module should be translated into the corresponding spline
  ! boundary condition. 

  enum, bind(C)
     enumerator :: PERIODIC_MAP_BC = 0, HERMITE_MAP_BC = 1
  end enum

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
  type map_2D   
     sll_int32  :: map_type        ! through functions or through data
     sll_int32  :: num_pts_1
     sll_int32  :: num_pts_2
     sll_real64, dimension(:,:), pointer :: x1_node   ! x1 = x1(eta1,eta2)
     sll_real64, dimension(:,:), pointer :: x2_node   ! x2 = x2(eta1,eta2)
     sll_real64, dimension(:,:), pointer :: x1_cell   ! x1 = x1(eta1,eta2)
     sll_real64, dimension(:,:), pointer :: x2_cell   ! x2 = x2(eta1,eta2)
     procedure(two_arg_scalar_function), pointer, nopass    :: x1_func
     procedure(two_arg_scalar_function), pointer, nopass    :: x2_func
     type(jacobian_matrix_element), dimension(:,:), pointer :: j_matrix
     sll_real64, dimension(:,:), pointer                    :: jacobians_n
     sll_real64, dimension(:,:), pointer                    :: jacobians_c
     type(sll_spline_2D), pointer                           :: x1_spline
     type(sll_spline_2D), pointer                           :: x2_spline
     procedure(two_arg_message_passing_func),pointer,nopass :: jacobian_func
  end type map_2D

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

  abstract interface
     function two_arg_message_passing_func( map, eta1, eta2 )
       use sll_working_precision
       import     :: map_2D
       sll_real64 :: two_arg_message_passing_func
       type(map_2D), pointer  :: map
       sll_real64, intent(in) :: eta1
       sll_real64, intent(in) :: eta2
     end function two_arg_message_passing_func
  end interface

  ! Here we try to represent the Jacobian matrix an actual 2D array of
  ! functions. But since fortran does not allow arrays of pointers, here
  ! we define a special type that can be used as an array element.
  type jacobian_matrix_element
     procedure(two_arg_scalar_function), pointer, nopass :: f
  end type jacobian_matrix_element

#if 0
  ! ---------------------------------------------------------------------
  !
  !   MAPPED MESHES: i.e. a mesh + coordinate transformation
  !
  ! ---------------------------------------------------------------------
  type mapped_mesh_2D_scalar
     sll_int32  :: boundary_type_1
     sll_int32  :: boundary_type_2
     sll_int32  :: mode            ! data defined on nodes or center of cells
     sll_real64, dimension(:,:), pointer :: data 
     type(sll_spline_2D), pointer        :: u_spline_2D  ! u is for uniform
     type(map_2D), pointer           :: map
  end type mapped_mesh_2D_scalar
#endif

  interface delete
     module procedure delete_map_2D
  end interface

contains

  ! -------------------------------------------------------------------------
  !
  !         FUNCTIONS AND SUBROUTINES FOR THE MAP_2D TYPE.
  !
  ! -------------------------------------------------------------------------

  ! new_map_2D() only allocates the memory for the object itself. The
  ! initialization routines will allocate the memory of the internal arrays.
  function new_map_2D( map_type )
    type(map_2D), pointer  :: new_map_2D
    sll_int32, intent(in)  :: map_type
    sll_int32              :: ierr
    SLL_ALLOCATE( new_map_2D, ierr)
    new_map_2D%map_type  = map_type
  end function new_map_2D

  ! Wrapper to compute the jacobian at the (eta1,eta2) point regardless
  ! of the type of mapping, analytic or discrete. This is the public
  ! interface to compute the continuous jacobian. This call could well be
  ! converted into a macro call.
  function jacobian_2D( map, eta1, eta2 )
    sll_real64            :: jacobian_2D
    type(map_2D), pointer :: map
    sll_real64            :: eta1
    sll_real64            :: eta2
    ! The following looks extremely ugly but one has to be aware that
    ! in principle, the 'map' argument to jacobian_func could be a different
    ! map than the host object. This is not entirely solved by changing
    ! the atribute of the procedure pointer to pass(map), as we would have
    ! the problem that we intend to pass a map pointer, not the map itself.
    jacobian_2D = (map%jacobian_func(map, eta1, eta2))
  end function jacobian_2D

  function jacobian_2D_analytic( map, eta1, eta2 )
    sll_real64             :: jacobian_2D_analytic
    type(map_2D), pointer  :: map
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    sll_real64             :: j11
    sll_real64             :: j12
    sll_real64             :: j21
    sll_real64             :: j22
    j11 = (map%j_matrix(1,1)%f( eta1, eta2 ))
    j12 = (map%j_matrix(1,2)%f( eta1, eta2 ))
    j21 = (map%j_matrix(2,1)%f( eta1, eta2 ))
    j22 = (map%j_matrix(2,2)%f( eta1, eta2 ))
    jacobian_2D_analytic = j11*j22 - j12*j21
  end function jacobian_2D_analytic

  function jacobian_2D_discrete( map, eta1, eta2 )
    sll_real64             :: jacobian_2D_discrete
    type(map_2D), pointer  :: map
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    sll_real64             :: j11
    sll_real64             :: j12
    sll_real64             :: j21
    sll_real64             :: j22
    j11 = interpolate_x1_derivative_2D( eta1, eta2, map%x1_spline )
    j12 = interpolate_x2_derivative_2D( eta1, eta2, map%x1_spline )
    j21 = interpolate_x1_derivative_2D( eta1, eta2, map%x2_spline )
    j22 = interpolate_x2_derivative_2D( eta1, eta2, map%x2_spline )
    jacobian_2D_discrete = j11*j22 - j12*j21
  end function jacobian_2D_discrete

  ! initialize_map_2D() allocates all the memory needed by the 2D map. 
  ! This interface is ending up very awkward because of the large amount of
  ! optional parameters that it takes, much of this in account of the 
  ! splines that it initializes, which take plenty of optional parameters
  ! themselves. This is not desirable and should be reassessed critically.
  !
  ! We should offer the possibility to pass the jacobian function directly.
  subroutine initialize_map_2D( &
    map,            &
    npts1,          &
    npts2,          &
    j11_func,       &
    j12_func,       &
    j21_func,       &
    j22_func,       &
    x1_func,        &
    x2_func,        &
    jacobians_cell, &
    jacobians_node, &
    x1_node,        &
    x2_node,        &
    x1_cell,        &
    x2_cell,        &
    eta1_bc_type_x1,   &
    eta2_bc_type_x1,   &
    eta1_min_slope_x1, &
    eta1_max_slope_x1, &
    eta2_min_slope_x1, &
    eta2_max_slope_x1, &
    eta1_bc_type_x2,   &
    eta2_bc_type_x2,   &
    eta1_min_slope_x2, &
    eta1_max_slope_x2, &
    eta2_min_slope_x2, &
    eta2_max_slope_x2 )

    type(map_2D), pointer  :: map
    sll_int32, intent(in)  :: npts1
    sll_int32, intent(in)  :: npts2
    procedure(two_arg_scalar_function), optional  :: j11_func
    procedure(two_arg_scalar_function), optional  :: j12_func
    procedure(two_arg_scalar_function), optional  :: j21_func
    procedure(two_arg_scalar_function), optional  :: j22_func
    procedure(two_arg_scalar_function), optional  :: x1_func
    procedure(two_arg_scalar_function), optional  :: x2_func
    sll_real64, dimension(:,:), optional          :: jacobians_node
    sll_real64, dimension(:,:), optional          :: x1_node
    sll_real64, dimension(:,:), optional          :: x2_node
    sll_real64, dimension(:,:), optional          :: jacobians_cell
    sll_real64, dimension(:,:), optional          :: x1_cell
    sll_real64, dimension(:,:), optional          :: x2_cell
    sll_int32, intent(in), optional               :: eta1_bc_type_x1
    sll_int32, intent(in), optional               :: eta2_bc_type_x1
    sll_real64, intent(in), optional              :: eta1_min_slope_x1
    sll_real64, intent(in), optional              :: eta1_max_slope_x1
    sll_real64, intent(in), optional              :: eta2_min_slope_x1
    sll_real64, intent(in), optional              :: eta2_max_slope_x1
    sll_int32, intent(in), optional               :: eta1_bc_type_x2
    sll_int32, intent(in), optional               :: eta2_bc_type_x2
    sll_real64, intent(in), optional              :: eta1_min_slope_x2
    sll_real64, intent(in), optional              :: eta1_max_slope_x2
    sll_real64, intent(in), optional              :: eta2_min_slope_x2
    sll_real64, intent(in), optional              :: eta2_max_slope_x2


    sll_int32  :: map_type ! enumerated constant-valued
    sll_real64 :: delta_1  ! cell spacing in eta1 
    sll_real64 :: delta_2  ! cell spacing in eta2 
    sll_real64 :: eta_1
    sll_real64 :: eta_2
    sll_int32  :: i
    sll_int32  :: j
    sll_int32  :: x1_eta1_bc ! to translate the BC enumerators to the splines
    sll_int32  :: x1_eta2_bc
    sll_int32  :: x2_eta1_bc
    sll_int32  :: x2_eta2_bc
    sll_int32  :: ierr
    logical    :: x1n
    logical    :: x2n
    logical    :: jn
    logical    :: x1c
    logical    :: x2c
    logical    :: jc
    logical    :: node_data_given
    logical    :: cell_data_given
    logical    :: x1_eta1_bc_given
    logical    :: x1_eta2_bc_given
    logical    :: x2_eta1_bc_given
    logical    :: x2_eta2_bc_given


    map_type = map%map_type
    x1n = present(x1_node)
    x2n = present(x2_node)
    jn  = present(jacobians_node)
    x1c = present(x1_cell)
    x2c = present(x2_cell)
    jc  = present(jacobians_cell)
    x1_eta1_bc_given = present(eta1_bc_type_x1)
    x1_eta2_bc_given = present(eta2_bc_type_x1)
    x2_eta1_bc_given = present(eta1_bc_type_x2)
    x2_eta2_bc_given = present(eta2_bc_type_x2)

    ! Check argument consistency
    !
    ! Maps can be ANALYTIC or NUMERIC. This determines which parameters are
    ! expected by the initialization function. Some optional parameters are
    ! expected by the ANALYTIC map and others by the NUMERIC map. These are
    ! mutually exclusive.
    !
    ! ANALYTIC_MAPs require all the optional arguments that describe the
    ! coordinate transformation. To wit: jii_func parameters, x1_func and
    ! x2_func.
    ! 
    ! DISCRETE_MAPs require only some of the parameters. If the mapping is
    ! defined from the nodes of the logical (eta1, eta2) mesh to the nodes
    ! of the physical mesh (x1,x2), then the node arrays are required:
    ! jacobians_node, x1_node and x2_node.
    ! If the transformation is done on the points at the center of the cells
    ! then these parameters are also required: 
    ! jacobians_cell, x1_cell, x2_cell.
    ! node and cell values are not mutually exclusive, thus all 6 parameters 
    ! can be provided in the NUMERIC case. It is up to the caller to make
    ! sure that the data set is consistent.

    if(.not. associated(map)) then
       print *, 'ERROR, initialize_map_2D(): passed map pointer was not ',&
            'associated.'
       STOP
    end if

    select case( map_type )
       case (ANALYTIC_MAP)
          if( (.not. present(j11_func)) .or. (.not. present(j12_func)) .or. &
              (.not. present(j21_func)) .or. (.not. present(j22_func)) ) then
             print *, 'ERROR, initialize_map_2D(): ANALYTIC_MAPs ', &
                  'require all j11_func, j12_func, j21_func and j22_func ', &
                  'parameters.'
             STOP
          end if
          if( (.not. present(x1_func)) .or. (.not. present(x2_func)) ) then
             print *, 'ERROR, initialize_map_2D(): ANALYTIC_MAPs ', &
                  'require x1_func and x2_func parameters.'
             STOP
          end if
          if( jn .or. jc .or.  x1n .or. x2n .or. x1c .or. x2c ) then
             print *, 'ERROR, initialize_map_2D(): ANALYTIC_MAPs ', &
                  'do not need any of the parameters required by the ', &
                  'DISCRETE_MAPs'
             STOP
          end if
          if( &
             x1_eta1_bc_given .or. x1_eta2_bc_given .or. &
             x2_eta1_bc_given .or. x2_eta2_bc_given ) then
             print *, 'ERROR, initialize_map_2D(): ANALYTIC_MAPs ', &
                  'do not need the specification of the boundary ', &
                  'conditions for the x1 or x2 transformations.'
             ! The following is a little drastic, after all, this info
             ! would simply not be used. 
             STOP
          end if
       case (DISCRETE_MAP)
          ! check that functions relevant to the ANALYTIC_MAP have not
          ! been passed.
          if( present(j11_func) .or. present(j12_func) .or. &
              present(j21_func) .or. present(j22_func) ) then
             print *, 'ERROR, initialize_map_2D(): DISCRETE_MAPs ', &
                  'do not need to be passed the elements of the jacobian ', &
                  'matrix. '
             STOP
          end if
          if( present(x1_func) .or. present(x2_func) ) then
             print *, 'ERROR, initialize_map_2D(): DISCRETE_MAPs ', &
                  'do not need x1_func and x2_func parameters.'
             STOP
          end if
          if( x1n .and. x2n ) then
             node_data_given = .true.
          else
             node_data_given = .false.
          end if
          if( jc .and. x1c .and. x2c ) then
             cell_data_given = .true.
          else
             cell_data_given = .false.
          end if
          ! Which combinations of numeric arguments make sense? In other 
          ! words, if the user passes all the node-based information, it 
          ! may be also convenient for some purposes to pass the cell-based 
          ! jacobians. But which are the combinations of arguments against
          ! which we should protect this function?
          !
          ! 0. The node-based information of the transformation is the
          !    absolute minimum required.
          if( (.not. x1n) .or. (.not. x2n) ) then
             print *, 'ERROR, initialize_map_2D(), DISCRETE_MAP case: ', &
                  'the node-based information (x1 and x2 arrays) is ', &
                  'the minimum information required.'
             STOP
          end if
          ! 1. If either of the (cell-based) x1 or x2 arrays is passed, 
          !    the other must also be.
          if( (x1c .and. (.not. x2c)) .or. ((.not. x1c) .and. x2c  ) ) then
             print *, 'ERROR, initialize_map_2D(), DISCRETE_MAP case: ', &
               'if either of the cell-based x1 or x2 arrays is passed, ', &
               'then the other must be passed as well.'
             STOP
          end if
          ! 2. Check that the discrete representation of the transformation is
          !    consistent with the size of the 2D array.
          if( &
             (size(x1_node,1) .lt. npts1) .or. &
             (size(x1_node,2) .lt. npts2) ) then
             print *, 'ERROR, initialize_map_2D(), DISCRETE_MAP case: ', &
                  'the size of the x1_node or x2_node arrays is ', &
                  'inconsistent with the number of points declared, ', &
                  'npts1 or npts2.'
             STOP
          end if
          if( jn .eqv. .true. ) then
             if( &
                (size(jacobians_node,1) .lt. npts1 - 1 ) .or. &
                (size(jacobians_node,2) .lt. npts2 - 1 ) ) then
                print *, 'ERROR, initialize_map_2D(), DISCRETE_MAP ', &
                     'case: the size of the jacobians_node array is ', &
                     'inconsistent with the number of points declared, ', &
                     'npts1 or npts2.'
                STOP
             end if
          else
             if( &
                (.not.x1_eta1_bc_given) .or. (.not.x1_eta2_bc_given) .or. &
                (.not.x2_eta1_bc_given) .or. (.not.x2_eta2_bc_given)  ) then
                print *, 'ERROR, initialize_map_2D(), DISCRETE_MAP ', &
                     'case: if discrete values for the jacobian are not ', &
                     'provided, then the specification of what boundary ', &
                     'conditions are wished for become necessary. i.e.: ', &
                     'pass eta1_bc_type_x1, etc.'
                STOP
             end if
          end if
          if( jc .eqv. .true. ) then
             if( &
                (size(jacobians_cell,1) .lt. npts1 - 1 ) .or. &
                (size(jacobians_cell,2) .lt. npts2 - 1 ) ) then
                print *, 'ERROR, initialize_map_2D(), DISCRETE_MAP ', &
                     'case: the size of the jacobians_cell arrays is ', &
                     'inconsistent with the number of points declared, ', &
                     'npts1 or npts2.'
                STOP
             end if
          end if
          ! 3. The discrete case requires the user to specify all the 
          !    boundary conditions for the x1 and x2 transformations.
          if( &
             (.not. x1_eta1_bc_given) .or. (.not. x1_eta2_bc_given) .or. &
             (.not. x2_eta1_bc_given) .or. (.not. x2_eta2_bc_given) ) then
             print *, 'ERROR, initialize_map_2D(), DISCRETE_MAP ', &
                  'case: it is required to pass all the boundary condition ', &
                  'specifications for the x1 and x2 transformations. '
             STOP
          end if
          ! More cases for the argument consistency should be added here...
    end select

    map%num_pts_1 = npts1
    map%num_pts_2 = npts2

    ! Allocate the arrays for precomputed jacobians.
    SLL_ALLOCATE(map%jacobians_n(npts1,npts2), ierr)
    SLL_ALLOCATE(map%jacobians_c(npts1-1, npts2-1), ierr)

    ! Allocation for x1 and x2 at nodes, needed regardless of the type of map
    SLL_ALLOCATE(map%x1_node(npts1,npts2), ierr)
    SLL_ALLOCATE(map%x2_node(npts1,npts2), ierr)

    ! Start filling out the fields and allocating the object's memory.
    !
    ! Implementation notes: It would have been very desirable to be able
    ! to implement the ANALYTIC map by analytic/functional means only, i.e.:
    ! the answer questions like x1_node(map, i, j) should be computed by
    ! functional means and not by looking up values in an array. This is
    ! doable but at the price of increasing the costs of the function calls 
    ! with additional pointer dereferencings and if() tests. It may very 
    ! well be that even after folding in the extra costs, this computation 
    ! would be more efficient than the array lookup especially in those 
    ! cases in which reading the value is not accompanied by a lot of 
    ! computations. It would have been desirable to have these different
    ! implementations and to compare the relative performance in the cases
    ! of interest.
    !
    ! Here, for simplicity, we make either type of map (DISCRETE or 
    ! ANALYTIC) use the same logic to answer questions like x1_node(map, i, j).
    select case (map_type)
       case (ANALYTIC_MAP)
  
          SLL_ALLOCATE(map%x1_cell(npts1-1, npts2-1), ierr)
          SLL_ALLOCATE(map%x2_cell(npts1-1, npts2-1), ierr)

          ! Fill the jacobian matrix
          SLL_ALLOCATE(map%j_matrix(2,2), ierr)
          map%j_matrix(1,1)%f => j11_func
          map%j_matrix(1,2)%f => j12_func
          map%j_matrix(2,1)%f => j21_func
          map%j_matrix(2,2)%f => j22_func
          map%jacobian_func =>jacobian_2D_analytic

          ! Assign the transformation functions
          map%x1_func => x1_func
          map%x1_func => x2_func

          ! Allocate the arrays for precomputed jacobians.
          SLL_ALLOCATE(map%jacobians_n(npts1,npts2), ierr)
          SLL_ALLOCATE(map%jacobians_c(npts1-1, npts2-1), ierr)

          ! Fill out the jacobians for the points in the uniform mesh.
          delta_1 = 1.0_f64/(npts1 - 1)
          delta_2 = 1.0_f64/(npts2 - 1)

          ! Fill the values at the nodes
          do j=0, npts2 - 1
             do i=0, npts1 - 1
                eta_1 = real(i,f64)*delta_1
                eta_2 = real(j,f64)*delta_2
                map%x1_node(i+1,j+1) = x1_func(eta_1, eta_2)
                map%x2_node(i+1,j+1) = x2_func(eta_1, eta_2)
                map%jacobians_n(i+1,j+1) = (map%jacobian_func(map,eta_1,eta_2))
             end do
          end do

          ! Fill the values at the mid-point of the cells
          do j=0, npts2 - 2
             do i=0, npts1 - 2
                eta_1 = delta_1*(real(i,f64) + 0.5_f64)
                eta_2 = delta_2*(real(j,f64) + 0.5_f64)
                map%x1_cell(i+1,j+1) = x1_func(eta_1, eta_2)
                map%x2_cell(i+1,j+1) = x2_func(eta_1, eta_2)
                map%jacobians_c(i+1,j+1) = (map%jacobian_func(map,eta_1,eta_2))
             end do
          end do

       case(DISCRETE_MAP)
          ! none of the direct functional capabilities are thus available
          map%x1_func  => null()
          map%x1_func  => null()
          map%j_matrix => null()

          ! allocate array memory locations
          SLL_ALLOCATE(map%x1_node(npts1,npts2), ierr)
          SLL_ALLOCATE(map%x2_node(npts1,npts2), ierr)

          ! initialize the local arrays. Note that since the map has its
          ! own copies, it owns this information locally and will destroy
          ! this information when the object is deleted. The caller is
          ! thus responsible for deallocating the arrays that were passed as
          ! arguments.
          do j=1, npts2
             do i=1, npts1
                map%x1_node(i,j) = x1_node(i,j)
                map%x2_node(i,j) = x2_node(i,j)
             end do
          end do

          ! The data from the discrete transformation is used to build
          ! the splines needed to compute the jacobian at any point.
          !
          ! First translate the enumerators from the map to the enumerators 
          ! that the splines recognize.
          select case ( eta1_bc_type_x1 )
             case (PERIODIC_MAP_BC)
                x1_eta1_bc = PERIODIC_SPLINE
             case (HERMITE_MAP_BC)
                x1_eta1_bc = HERMITE_SPLINE
             case default
                print *, 'ERROR, initialize_map_2D(): unrecognized ', &
                     'boundary type for x1, eta1 direction.'
                STOP
          end select

          select case ( eta2_bc_type_x1 )
             case (PERIODIC_MAP_BC)
                x1_eta2_bc = PERIODIC_SPLINE
             case (HERMITE_MAP_BC)
                x1_eta2_bc = HERMITE_SPLINE
             case default
                print *, 'ERROR, initialize_map_2D(): unrecognized ', &
                     'boundary type for x1, eta2 direction.'
                STOP
          end select

          select case ( eta1_bc_type_x2 )
             case (PERIODIC_MAP_BC)
                x2_eta1_bc = PERIODIC_SPLINE
             case (HERMITE_MAP_BC)
                x2_eta1_bc = HERMITE_SPLINE
             case default
                print *, 'ERROR, initialize_map_2D(): unrecognized ', &
                     'boundary type for x2, eta1 direction.'
                STOP
          end select

          select case ( eta2_bc_type_x2 )
             case (PERIODIC_MAP_BC)
                x2_eta2_bc = PERIODIC_SPLINE
             case (HERMITE_MAP_BC)
                x2_eta2_bc = HERMITE_SPLINE
             case default
                print *, 'ERROR, initialize_map_2D(): unrecognized ', &
                     'boundary type for x2, eta2 direction.'
                STOP
          end select

          ! allocate the splines for computing the jacobian terms
          map%x1_spline => new_spline_2D( &
               npts1, &
               npts2, &
               0.0_f64, & 
               1.0_f64, &
               0.0_f64, &
               1.0_f64, &
               x1_eta1_bc, &
               x1_eta2_bc, &
               eta1_min_slope_x1, &
               eta1_max_slope_x1, &
               eta2_min_slope_x1, &
               eta2_max_slope_x1 )
          
          map%x2_spline => new_spline_2D( &
               npts1, &
               npts2, &
               0.0_f64, & 
               1.0_f64, &
               0.0_f64, &
               1.0_f64, &
               x2_eta1_bc, &
               x2_eta1_bc, &
               eta1_min_slope_x2, &
               eta1_max_slope_x2, &
               eta2_min_slope_x2, &
               eta2_max_slope_x2 )

          ! Compute the spline coefficients
          call compute_spline_2D( map%x1_node, map%x1_spline )
          call compute_spline_2D( map%x2_node, map%x2_spline )

          ! Initialize the jacobian function with the right entity
          map%jacobian_func =>jacobian_2D_discrete

          ! The splines contain all the information to compute the
          ! jacobians everywhere; however, here we explore assigning
          ! the jacobians-at-the-nodes array with the values provided
          ! by the user if available. If there are discrepancies between
          ! the user-provided values and the predictions from the splines,
          ! then this may itself be a way to look for errors.
          !
          ! Copy the values of the jacobians at the nodes is user given:
          if( jn .eqv. .true. ) then
             do j=1, npts2
                do i=1, npts1
                   map%jacobians_n(i,j) = jacobians_node(i,j)
                end do
             end do
          else

         ! Fill the jacobian values at the nodes calculated from the splines
             do j=0, npts2 - 1
                do i=0, npts1 - 1
                   eta_1 = real(i,f64)*delta_1
                   eta_2 = real(j,f64)*delta_2
                   map%jacobians_n(i+1,j+1) = &
                        (map%jacobian_func(map,eta_1,eta_2)) 
                end do
             end do
          end if

          ! copy the cell-based transformation arrays if available
          if( (x1c .and. x2c) .eqv. .true. ) then
             SLL_ALLOCATE(map%x1_cell(npts1-1, npts2-1), ierr)
             SLL_ALLOCATE(map%x2_cell(npts1-1, npts2-1), ierr)
             do j=1, npts2 - 1
                do i=1, npts1 - 1
                   map%x1_cell(i,j) = x1_cell(i,j)
                   map%x2_cell(i,j) = x2_cell(i,j)
                end do
             end do
          end if
          ! copy the cell-based jacobians if available
          if( jc .eqv. .true. ) then
             do j=1, npts2 - 1
                do i=1, npts1 - 1
                   map%jacobians_c(i,j) = jacobians_cell(i,j)
                end do
             end do
          end if
       end select
  end subroutine initialize_map_2D


  subroutine delete_map_2D( map )
    type(map_2D), pointer :: map
    sll_int32             :: ierr
    if( .not. associated(map) ) then
       print *, 'ERROR, delete_map_2D: passed map pointer is not associated.'
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
  end subroutine delete_map_2D


  ! Access functions for the mapping. These are an overkill and can be 
  ! changed by a macro, but for now, they are at least safer.
  function map2d_x1_node( map, i, j )
    sll_real64            :: map2d_x1_node
    type(map_2D), pointer :: map
    sll_int32, intent(in) :: i
    sll_int32, intent(in) :: j
    sll_int32 :: num_pts_1
    sll_int32 :: num_pts_2
    SLL_ASSERT( associated(map) )
    num_pts_1 = map%num_pts_1
    num_pts_2 = map%num_pts_2
    SLL_ASSERT( (i .ge. 1) .and. (i .le. num_pts_1) )
    SLL_ASSERT( (j .ge. 1) .and. (j .le. num_pts_2) )
    map2d_x1_node = map%x1_node(i,j)
  end function map2d_x1_node

  function map2d_x2_node( map, i, j )
    sll_real64            :: map2d_x2_node
    type(map_2D), pointer :: map
    sll_int32, intent(in) :: i
    sll_int32, intent(in) :: j
    sll_int32 :: num_pts_1
    sll_int32 :: num_pts_2
    SLL_ASSERT( associated(map) )
    num_pts_1 = map%num_pts_1
    num_pts_2 = map%num_pts_2
    SLL_ASSERT( (i .ge. 1) .and. (i .le. num_pts_1 - 1) )
    SLL_ASSERT( (j .ge. 1) .and. (j .le. num_pts_2 - 1) )
    map2d_x2_node = map%x2_node(i,j)
  end function map2d_x2_node

  function map2d_x1_cell( map, i, j )
    sll_real64            :: map2d_x1_cell
    type(map_2D), pointer :: map
    sll_int32, intent(in) :: i
    sll_int32, intent(in) :: j
    sll_int32 :: num_pts_1
    sll_int32 :: num_pts_2
    SLL_ASSERT( associated(map) )
    num_pts_1 = map%num_pts_1
    num_pts_2 = map%num_pts_2
    SLL_ASSERT( (i .ge. 1) .and. (i .le. num_pts_1 - 1) )
    SLL_ASSERT( (j .ge. 1) .and. (j .le. num_pts_2 - 1) )
    map2d_x1_cell = map%x1_cell(i,j)
  end function map2d_x1_cell

  function map2d_x2_cell( map, i, j )
    sll_real64            :: map2d_x2_cell
    type(map_2D), pointer :: map
    sll_int32, intent(in) :: i
    sll_int32, intent(in) :: j
    sll_int32 :: num_pts_1
    sll_int32 :: num_pts_2
    SLL_ASSERT( associated(map) )
    num_pts_1 = map%num_pts_1
    num_pts_2 = map%num_pts_2
    SLL_ASSERT( (i .ge. 1) .and. (i .le. num_pts_1) )
    SLL_ASSERT( (j .ge. 1) .and. (j .le. num_pts_2) )
    map2d_x2_cell = map%x2_cell(i,j)
  end function map2d_x2_cell

  function map_2d_jacobian_node( map, i, j )
    sll_real64              :: map_2d_jacobian_node
    type(map_2D), pointer   :: map
    sll_int32, intent(in)   :: i
    sll_int32, intent(in)   :: j
    sll_int32 :: num_pts_1
    sll_int32 :: num_pts_2
    SLL_ASSERT( associated(map) )
    num_pts_1 = map%num_pts_1
    num_pts_2 = map%num_pts_2
    SLL_ASSERT( (i .ge. 1) .and. (i .le. num_pts_1) )
    SLL_ASSERT( (j .ge. 1) .and. (j .le. num_pts_2) )
    map_2d_jacobian_node = map%jacobians_n(i,j)
  end function map_2d_jacobian_node

  function map_2d_jacobian_cell( map, i, j )
    sll_real64              :: map_2d_jacobian_cell
    type(map_2D), pointer   :: map
    sll_int32, intent(in)   :: i
    sll_int32, intent(in)   :: j
    sll_int32 :: num_cells_1
    sll_int32 :: num_cells_2
    SLL_ASSERT( associated(map) )
    num_cells_1 = map%num_pts_1 - 1
    num_cells_2 = map%num_pts_2 - 1
    SLL_ASSERT( (i .ge. 1) .and. (i .le. num_cells_1) )
    SLL_ASSERT( (j .ge. 1) .and. (j .le. num_cells_2) )
    map_2d_jacobian_cell = map%jacobians_c(i,j)
  end function map_2d_jacobian_cell

#if 0
  ! -------------------------------------------------------------------------
  !
  !     FUNCTIONS AND SUBROUTINES FOR THE MAPPED MESH_2D TYPE.
  !
  ! -------------------------------------------------------------------------

  ! The mesh_2D_scalar is basically data interpreted with the aid of a 
  ! 2D coordinate transformation (the map_2D). Additionally, we need 
  ! information on:
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

    type(mapped_mesh_2D_scalar), pointer    :: new_mesh_2D_scalar
    sll_int32, intent(in)            :: mode
    sll_int32, intent(in)            :: bc_1
    sll_int32, intent(in)            :: bc_2
    type(map_2D), pointer            :: mapping
    sll_real64, intent(in), optional :: eta1_min_slope
    sll_real64, intent(in), optional :: eta1_max_slope
    sll_real64, intent(in), optional :: eta2_min_slope
    sll_real64, intent(in), optional :: eta2_max_slope

    sll_int32 :: num_pts_1
    sll_int32 :: num_pts_2
    sll_int32 :: ierr
    sll_int32 :: spline_bc_1
    sll_int32 :: spline_bc_2

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
    type(mapped_mesh_2d_scalar), pointer :: mesh
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
    type(mapped_mesh_2D_scalar), pointer :: mesh
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
    type(mapped_mesh_2D_scalar), pointer :: mesh
    sll_int32, intent(in)         :: i
    sll_int32, intent(in)         :: j
    SLL_ASSERT( associated(mesh) )
    get_m2ds_node = mesh%data(i,j)
  end function get_m2ds_node

  ! Add here get_m2ds_val(mesh, x1, x2), this requires the nonuniform splines.

  subroutine set_m2ds_node( mesh, i, j, val )
    type(mapped_mesh_2D_scalar), pointer :: mesh
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
    type(mapped_mesh_2D_scalar), pointer :: mesh
    sll_real64                    :: eta1
    sll_real64                    :: eta2
    get_mesh2D_value = interpolate_value_2D( eta1, eta2, mesh%u_spline_2D )
  end function get_mesh2D_value

  function mesh_2D_data( mesh )
    sll_real64, dimension(:,:), pointer :: mesh_2D_data
    type(mapped_mesh_2D_scalar), pointer       :: mesh
    mesh_2D_data => mesh%data
  end function mesh_2D_data
#endif
end module sll_coordinate_transformation
