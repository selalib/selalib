module sll_module_coordinate_transformations_2d
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_file_io.h"
  use sll_cubic_splines
  use sll_xdmf
  use sll_logical_meshes
#ifdef STDF95
  use sll_cubic_spline_interpolator_2d
  use sll_gnuplot
#else
  use sll_module_interpolators_2d_base
  use sll_coordinate_transformation_2d_base_module
  use sll_module_deboor_splines_2d


#endif
  implicit none
  
  ! ---------------------------------------------------------------------
  !
  !                       2D COORDINATE TRANSFORMATIONS
  !
  ! ---------------------------------------------------------------------
  !
  ! There are two main types of coordinate transformations: analytic and
  ! discrete. In the first case the transformation can be
  ! specified analytically, through the two functions:
  !
  !                     x1 = x1(eta1,eta2) 
  !                     x2 = x2(eta1,eta2)
  !
  ! Where both, eta1 and eta2 should be defined on the intervals that define
  ! the extent of the logical mesh (default values in logical mesh are  [0,1]. 
  ! The same transformation 
  ! can be specified by the set of transformed points x1(i,j), x2(i,j), as
  ! two 2D arrays or 1D arrays that describe the transformation on each
  ! direction.
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

!> Analytic transformation
#ifdef STDF95
  type  :: sll_coordinate_transformation_2d_analytic
     ! for pre-computing values. Not clear how advantageous this is.
!!$     sll_real64, dimension(:,:), pointer :: x1_cell
!!$     sll_real64, dimension(:,:), pointer :: x2_cell
!!$     sll_real64, dimension(:,:), pointer :: jacobians_n
!!$     sll_real64, dimension(:,:), pointer :: jacobians_c
     !ES the two following need to be in the base class
     !character(len=64) :: label
     !logical           :: written! = .false.
     type(sll_logical_mesh_2d), pointer :: mesh
#else
  type, extends(sll_coordinate_transformation_2d_base):: &
       sll_coordinate_transformation_2d_analytic
#endif
!!$     sll_real64, dimension(:,:), pointer :: x1_node   ! x1(i,j) 
!!$     sll_real64, dimension(:,:), pointer :: x2_node   ! x2(i,j)
#ifdef STDF95
#else
     !character(len=64) :: label
     !logical           :: written! = .false.
     type(jacobian_matrix_element), dimension(:,:), pointer :: j_matrix
     procedure(transformation_func_nopass), pointer, nopass :: x1_func  ! user
     procedure(transformation_func_nopass), pointer, nopass :: x2_func  ! user
     procedure(two_arg_message_passing_func_analyt), pointer, pass :: &
          jacobian_func
     procedure(j_matrix_f_nopass), pointer, nopass :: jacobian_matrix_function
     sll_real64, dimension(:), pointer :: params => null() ! transf params
   contains
     procedure, pass(transf) :: initialize => initialize_coord_transf_2d_analytic
     procedure, pass(transf) :: get_logical_mesh => get_logical_mesh_analytic
     ! Functions with integer arguments
     procedure, pass(transf) :: x1_at_node => x1_node_analytic
     procedure, pass(transf) :: x2_at_node => x2_node_analytic
     procedure, pass(transf) :: jacobian_at_node => jacobian_node_analytic
     procedure, pass(transf) :: x1_at_cell => x1_cell_analytic
     procedure, pass(transf) :: x2_at_cell => x2_cell_analytic
     procedure, pass(transf) :: jacobian_at_cell => jacobian_2d_cell_analytic
     ! Functions with real arguments
     procedure, pass(transf) :: x1         => x1_analytic
     procedure, pass(transf) :: x2         => x2_analytic
     procedure, pass(transf) :: jacobian   => jacobian_2d_analytic
     procedure, pass(transf) :: jacobian_matrix => jacobian_matrix_2d_analytic
     procedure, pass(transf) :: inverse_jacobian_matrix => &
          inverse_jacobian_matrix_2d_analytic
     procedure, pass(transf) :: write_to_file => write_to_file_2d_analytic
     procedure, pass(transf) :: read_from_file => read_from_file_2d_analytic
     procedure, pass(transf) :: delete => delete_transformation_2d_analytic
#endif
  end type sll_coordinate_transformation_2d_analytic


  ! -----------------------------------------------------------------------
  !
  !                          Discrete case
  !
  ! -----------------------------------------------------------------------

#ifdef STDF95
  type                      ::sll_coordinate_transformation_2d_discrete
#else
  type, extends(sll_coordinate_transformation_2d_base) :: &
       sll_coordinate_transformation_2d_discrete
#endif
     sll_real64, dimension(:,:), pointer :: x1_node =>null()   ! x1(i,j) 
     sll_real64, dimension(:,:), pointer :: x2_node =>null()  ! x2(i,j) 
     sll_real64, dimension(:,:), pointer :: x1_cell =>null()
     sll_real64, dimension(:,:), pointer :: x2_cell =>null()
     sll_real64, dimension(:,:), pointer :: jacobians_n =>null()
     sll_real64, dimension(:,:), pointer :: jacobians_c =>null()
#ifdef STDF95
     ! so the choice is to have only cubic splines in f95 mode... 
     ! more generality would be obtained with arbitrary degree splines...
     type(cubic_spline_2d_interpolator), pointer            :: x1_interp
     type(cubic_spline_2d_interpolator), pointer            :: x2_interp
#else
!     type(jacobian_matrix_element), dimension(:,:), pointer :: j_matrix
     class(sll_interpolator_2d_base), pointer               :: x1_interp
     class(sll_interpolator_2d_base), pointer               :: x2_interp
   contains
     procedure, pass(transf) :: initialize => &
          initialize_coord_transf_2d_discrete
     procedure, pass(transf) :: get_logical_mesh => get_logical_mesh_discrete
     procedure, pass(transf) :: x1_at_node => x1_node_discrete
     procedure, pass(transf) :: x2_at_node => x2_node_discrete
     procedure, pass(transf) :: jacobian_at_node => transf_2d_jacobian_node_discrete
     procedure, pass(transf) :: x1         => x1_discrete
     procedure, pass(transf) :: x2         => x2_discrete
     procedure, pass(transf) :: x1_at_cell => x1_cell_discrete
     procedure, pass(transf) :: x2_at_cell => x2_cell_discrete
     procedure, pass(transf) :: jacobian_at_cell => jacobian_2d_cell_discrete
     procedure, pass(transf) :: jacobian   => jacobian_2d_discrete
     procedure, pass(transf) :: jacobian_matrix => jacobian_matrix_2d_discrete
     procedure, pass(transf) :: inverse_jacobian_matrix => &
          inverse_jacobian_matrix_2d_discrete
     procedure, pass(transf) :: write_to_file => write_to_file_2d_discrete
     procedure, pass(transf) :: read_from_file => read_from_file_2d_discrete
     procedure, pass(transf) :: delete => delete_transformation_2d_discrete
#endif
  end type sll_coordinate_transformation_2d_discrete

#ifdef STDF95
  interface initialize 
     module procedure initialize_coord_transf_2d_discrete, &
                      initialize_coord_transf_2d_analytic
  end interface

  interface x1_at_node 
     module procedure x1_node_analytic, x1_node_discrete
  end interface

  interface x2_at_node 
     module procedure x2_node_analytic, x2_node_discrete
  end interface

  interface x1 
     module procedure x1_discrete
  end interface

  interface x2
     module procedure x2_discrete
  end interface

  interface jacobian
     module procedure jacobian_2d_discrete 
  end interface

  interface write_to_file
     module procedure mma_write_to_file_2d_analytic, write_to_file_2d_discrete 
  end interface
#else
  abstract interface
     function j_matrix_f_nopass ( eta1, eta2, params ) result(val)
       use sll_working_precision
       sll_real64, dimension(2,2)   :: val
       sll_real64   :: eta1
       sll_real64   :: eta2
       sll_real64, dimension(:), optional, intent(in) :: params
     end function j_matrix_f_nopass
  end interface

   abstract interface
      function two_arg_message_passing_func_discr( transf, eta1, eta2 )
        use sll_working_precision
        import     :: sll_coordinate_transformation_2d_discrete
        sll_real64                      :: two_arg_message_passing_func_discr
        class(sll_coordinate_transformation_2d_discrete)  :: transf
        sll_real64, intent(in)          :: eta1
        sll_real64, intent(in)          :: eta2
      end function two_arg_message_passing_func_discr
   end interface

   abstract interface
      function two_arg_message_passing_func_analyt( transf, eta1, eta2 )
        use sll_working_precision
        import     :: sll_coordinate_transformation_2d_analytic
        sll_real64                      :: two_arg_message_passing_func_analyt
        class(sll_coordinate_transformation_2d_analytic) :: transf
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
#endif
  

  interface delete
     module procedure &
          delete_transformation_2d_analytic, &
          delete_transformation_2d_discrete
  end interface

  
contains

  !**************************************************************************
  !
  !       Functions for the analytic coordinate transformation
  !
  !**************************************************************************

  function new_coordinate_transformation_2d_analytic( &
    label,          &
    mesh_2d,        &
    x1_func,        &
    x2_func,        &
    j11_func,       &
    j12_func,       &
    j21_func,       &
    j22_func,       &
    params )

    type(sll_coordinate_transformation_2d_analytic), pointer :: &
         new_coordinate_transformation_2d_analytic
    character(len=*), intent(in)                  :: label
    type(sll_logical_mesh_2d), pointer :: mesh_2d
#ifdef STDF95
    sll_real64            :: x1_func
    sll_real64            :: x2_func
    sll_real64            :: j11_func
    sll_real64            :: j12_func
    sll_real64            :: j21_func
    sll_real64            :: j22_func
#else
    procedure(transformation_func_nopass)            :: x1_func
    procedure(transformation_func_nopass)            :: x2_func
    procedure(transformation_func_nopass)            :: j11_func
    procedure(transformation_func_nopass)            :: j12_func
    procedure(transformation_func_nopass)            :: j21_func
    procedure(transformation_func_nopass)            :: j22_func
#endif
    sll_real64, dimension(:), intent(in) :: params
    sll_int32 :: ierr

    SLL_ALLOCATE(new_coordinate_transformation_2d_analytic, ierr)
    call initialize_coord_transf_2d_analytic( &
         new_coordinate_transformation_2d_analytic, &
         label,          &
         mesh_2d,        &
         x1_func,        &
         x2_func,        &
         j11_func,       &
         j12_func,       &
         j21_func,       &
         j22_func,       &
         params )
  end function new_coordinate_transformation_2d_analytic

  subroutine initialize_coord_transf_2d_analytic( &
    transf, &
    label,          &
    mesh_2d,        &
    x1_func,        &
    x2_func,        &
    j11_func,       &
    j12_func,       &
    j21_func,       &
    j22_func,       &
    params )

#ifdef STDF95
    type(sll_coordinate_transformation_2d_analytic), intent(inout) :: &
         transf
#else
    class(sll_coordinate_transformation_2d_analytic), intent(inout) :: &
         transf
#endif
    character(len=*), intent(in)                  :: label
#ifdef STDF95
    sll_real64            :: x1_func
    sll_real64            :: x2_func
    sll_real64            :: j11_func
    sll_real64            :: j12_func
    sll_real64            :: j21_func
    sll_real64            :: j22_func
#else
    procedure(transformation_func_nopass)            :: x1_func
    procedure(transformation_func_nopass)            :: x2_func
    procedure(transformation_func_nopass)            :: j11_func
    procedure(transformation_func_nopass)            :: j12_func
    procedure(transformation_func_nopass)            :: j21_func
    procedure(transformation_func_nopass)            :: j22_func
#endif
    type(sll_logical_mesh_2d), pointer :: mesh_2d
    sll_real64, dimension(:), intent(in), optional :: params
    sll_int32  :: npts1
    sll_int32  :: npts2
    sll_real64 :: delta_1
    sll_real64 :: delta_2
    sll_int32  :: ierr

    transf%label   = trim(label)
    transf%mesh => mesh_2d
    npts1   = mesh_2d%num_cells1 + 1
    npts2   = mesh_2d%num_cells2 + 1
    delta_1 = mesh_2d%delta_eta1
    delta_2 = mesh_2d%delta_eta2

    ! Seriously consider eliminating this to have a lighter object and 
    ! doing all calculations on-the-fly.
    ! Allocate the arrays for precomputed jacobians.
!!$    SLL_ALLOCATE(transformation%jacobians_n(npts1,npts2), ierr)
!!$    SLL_ALLOCATE(transformation%jacobians_c(npts1-1, npts2-1), ierr)

    ! Allocation for x1 and x2 at nodes, needed regardless of the type of map
!!$    SLL_ALLOCATE(transformation%x1_node(npts1,npts2), ierr)
!!$    SLL_ALLOCATE(transformation%x2_node(npts1,npts2), ierr)

    ! Start filling out the fields and allocating the object's memory.
!!$    SLL_ALLOCATE(transformation%x1_cell(npts1-1, npts2-1), ierr)
!!$    SLL_ALLOCATE(transformation%x2_cell(npts1-1, npts2-1), ierr)

#ifdef STDF95
#else
    ! Assign the transformation functions and parameters
    transf%x1_func => x1_func
    transf%x2_func => x2_func
    if( present(params) ) then
       SLL_ALLOCATE(transf%params(size(params)),ierr)
       transf%params(:) = params(:)
    end if
    ! Fill the jacobian matrix
    SLL_ALLOCATE(transf%j_matrix(2,2), ierr)
    transf%j_matrix(1,1)%f => j11_func
    transf%j_matrix(1,2)%f => j12_func
    transf%j_matrix(2,1)%f => j21_func
    transf%j_matrix(2,2)%f => j22_func
    transf%jacobian_func   => jacobian_2d_analytic
#endif
    
    ! Allocate the arrays for precomputed jacobians.
!!$    SLL_ALLOCATE(transformation%jacobians_n(npts1,   npts2), ierr)
!!$    SLL_ALLOCATE(transformation%jacobians_c(npts1-1, npts2-1), ierr)
    
    ! Fill the values of the transformation and the jacobians at the nodes
!!$    do j=0, npts2 - 1
!!$       eta_2 = real(j,f64)*delta_2
!!$       do i=0, npts1 - 1
!!$          eta_1 = real(i,f64)*delta_1
!!$          transformation%x1_node(i+1,j+1) = x1_func(eta_1, eta_2)
!!$          transformation%x2_node(i+1,j+1) = x2_func(eta_1, eta_2)
!!$          ! for some compiler reason, the following intermediate 
!!$          ! variable is required, else the jacobians_n array will not
!!$          ! be filled out properly.
!!$#ifdef STDF95
!!$          ! We can't define jacobian_2d_analytic because it use procedure pointer
!!$          ! So we compute directly the jacobian with their components
!!$          jacobian_val          = j11_func(eta_1,eta_2)*j22_func(eta_1,eta_2) &
!!$                                - j12_func(eta_1,eta_2)*j21_func(eta_1,eta_2)
!!$#else
!!$          jacobian_val          = transformation%jacobian_func(eta_1,eta_2)
!!$#endif
!!$          transformation%jacobians_n(i+1,j+1) = jacobian_val
!!$       end do
!!$    end do
!!$    
!!$    ! Fill the values at the mid-point of the cells
!!$    do j=0, npts2 - 2
!!$       eta_2 = delta_2*(real(j,f64) + 0.5_f64)
!!$       do i=0, npts1 - 2
!!$          eta_1 = delta_1*(real(i,f64) + 0.5_f64)
!!$          transformation%x1_cell(i+1,j+1)     = x1_func(eta_1, eta_2)
!!$          transformation%x2_cell(i+1,j+1)     = x2_func(eta_1, eta_2)
!!$#ifdef STDF95
!!$          transformation%jacobians_c(i+1,j+1) = &
!!$               j11_func(eta_1,eta_2)*j22_func(eta_1,eta_2) &
!!$             - j12_func(eta_1,eta_2)*j21_func(eta_1,eta_2)
!!$#else
!!$          transformation%jacobians_c(i+1,j+1) = &
!!$               transformation%jacobian_func(eta_1,eta_2)
!!$#endif
!!$       end do
!!$    end do
  end subroutine initialize_coord_transf_2d_analytic

  subroutine delete_transformation_2d_analytic( transf )
    class(sll_coordinate_transformation_2d_analytic), intent(inout) :: transf
    sll_int32 :: ierr
    if(associated(transf%j_matrix)) then
       SLL_DEALLOCATE( transf%j_matrix, ierr )
    end if
    nullify( transf%x1_func )
    nullify( transf%x2_func )
    nullify( transf%jacobian_func )
    nullify( transf%jacobian_matrix_function )
  end subroutine delete_transformation_2d_analytic

  function get_logical_mesh_analytic( transf ) result(res)
    class(sll_coordinate_transformation_2d_analytic), intent(in) :: transf
    type(sll_logical_mesh_2d), pointer :: res
    res => transf%mesh
  end function get_logical_mesh_analytic

#ifdef STDF95
#else
  function jacobian_2d_analytic( transf, eta1, eta2 ) result(val)
    sll_real64                        :: val
    class(sll_coordinate_transformation_2d_analytic) :: transf
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    sll_real64             :: j11
    sll_real64             :: j12
    sll_real64             :: j21
    sll_real64             :: j22
    j11 = (transf%j_matrix(1,1)%f( eta1, eta2, transf%params ))
    j12 = (transf%j_matrix(1,2)%f( eta1, eta2, transf%params ))
    j21 = (transf%j_matrix(2,1)%f( eta1, eta2, transf%params ))
    j22 = (transf%j_matrix(2,2)%f( eta1, eta2, transf%params ))
    ! For debugging:
    !    print *, 'jacobian_2d_analytic: '
    !    print *, j11, j12
    !    print *, j21, j22
    val = j11*j22 - j12*j21
  end function jacobian_2d_analytic

  ! The efficiency of the following function could be improved by just
  ! passing the output array rather than returning values on the stack which
  ! need to be caught by the caller.
  function jacobian_matrix_2d_analytic( transf, eta1, eta2 )
    sll_real64, dimension(1:2,1:2)     :: jacobian_matrix_2d_analytic
    class(sll_coordinate_transformation_2d_analytic) :: transf
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    sll_real64             :: j11
    sll_real64             :: j12
    sll_real64             :: j21
    sll_real64             :: j22
    j11 = (transf%j_matrix(1,1)%f( eta1, eta2, transf%params ))
    j12 = (transf%j_matrix(1,2)%f( eta1, eta2, transf%params ))
    j21 = (transf%j_matrix(2,1)%f( eta1, eta2, transf%params ))
    j22 = (transf%j_matrix(2,2)%f( eta1, eta2, transf%params ))
    ! For debugging:
    !    print *, 'jacobian_2d_analytic: '
    !    print *, j11, j12
    !    print *, j21, j22
    jacobian_matrix_2d_analytic(1,1) = j11
    jacobian_matrix_2d_analytic(1,2) = j12
    jacobian_matrix_2d_analytic(2,1) = j21
    jacobian_matrix_2d_analytic(2,2) = j22
  end function jacobian_matrix_2d_analytic

  function inverse_jacobian_matrix_2d_analytic( transf, eta1, eta2 )
    sll_real64, dimension(1:2,1:2)     :: inverse_jacobian_matrix_2d_analytic
    class(sll_coordinate_transformation_2d_analytic) :: transf
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    sll_real64             :: inv_j11
    sll_real64             :: inv_j12
    sll_real64             :: inv_j21
    sll_real64             :: inv_j22
    sll_real64             :: r_jacobian ! reciprocal of the jacobian

    r_jacobian = 1.0_f64/transf%jacobian( eta1, eta2 )
    inv_j11 =  (transf%j_matrix(2,2)%f( eta1, eta2, transf%params ))*r_jacobian
    inv_j12 = -(transf%j_matrix(1,2)%f( eta1, eta2, transf%params ))*r_jacobian
    inv_j21 = -(transf%j_matrix(2,1)%f( eta1, eta2, transf%params ))*r_jacobian
    inv_j22 =  (transf%j_matrix(1,1)%f( eta1, eta2, transf%params ))*r_jacobian
    ! For debugging:
    !    print *, 'jacobian_2d_analytic: '
    !    print *, j11, j12
    !    print *, j21, j22
    inverse_jacobian_matrix_2d_analytic(1,1) = inv_j11
    inverse_jacobian_matrix_2d_analytic(1,2) = inv_j12
    inverse_jacobian_matrix_2d_analytic(2,1) = inv_j21
    inverse_jacobian_matrix_2d_analytic(2,2) = inv_j22
  end function inverse_jacobian_matrix_2d_analytic

  function x1_analytic( transf, eta1, eta2 ) result(val)
    sll_real64                         :: val
    class(sll_coordinate_transformation_2d_analytic) :: transf
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    val = transf%x1_func(eta1, eta2, transf%params)
  end function x1_analytic

  function x2_analytic( transf, eta1, eta2 ) result(val)
    sll_real64                        :: val
    class(sll_coordinate_transformation_2d_analytic) :: transf
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    val = transf%x2_func(eta1, eta2, transf%params)
  end function x2_analytic
#endif

  function x1_node_analytic( transf, i, j ) result(val)
#ifdef STDF95
    type(sll_coordinate_transformation_2d_analytic) :: transf
#else
    class(sll_coordinate_transformation_2d_analytic) :: transf
#endif
    sll_real64             :: val
    sll_int32, intent(in) :: i
    sll_int32, intent(in) :: j
    sll_real64            :: eta1
    sll_real64            :: eta2
    sll_real64            :: eta1_min
    sll_real64            :: eta2_min
    sll_real64            :: delta_eta1
    sll_real64            :: delta_eta2
    eta1_min   = transf%mesh%eta1_min
    eta2_min   = transf%mesh%eta2_min
    delta_eta1 = transf%mesh%delta_eta1
    delta_eta2 = transf%mesh%delta_eta2
    eta1       = eta1_min + real(i-1,f64)*delta_eta1
    eta2       = eta2_min + real(j-1,f64)*delta_eta2
    val = transf%x1_func(eta1,eta2,transf%params)
  end function x1_node_analytic

  function x2_node_analytic( transf, i, j ) result(val)
#ifdef STDF95
    type(sll_coordinate_transformation_2d_analytic) :: transf
#else
    class(sll_coordinate_transformation_2d_analytic) :: transf
#endif
    sll_real64             :: val
    sll_int32, intent(in) :: i
    sll_int32, intent(in) :: j
    sll_real64            :: eta1
    sll_real64            :: eta2
    sll_real64            :: eta1_min
    sll_real64            :: eta2_min
    sll_real64            :: delta_eta1
    sll_real64            :: delta_eta2
    eta1_min   = transf%mesh%eta1_min
    eta2_min   = transf%mesh%eta2_min
    delta_eta1 = transf%mesh%delta_eta1
    delta_eta2 = transf%mesh%delta_eta2
    eta1       = eta1_min + real(i-1,f64)*delta_eta1
    eta2       = eta2_min + real(j-1,f64)*delta_eta2
    val = transf%x2_func(eta1,eta2,transf%params)
  end function x2_node_analytic

  function x1_cell_analytic( transf, i, j ) result(var)
#ifdef STDF95
    type(sll_coordinate_transformation_2d_analytic) :: transf
#else
    class(sll_coordinate_transformation_2d_analytic) :: transf
#endif
    sll_real64            :: var
    sll_int32, intent(in) :: i
    sll_int32, intent(in) :: j
    sll_real64 :: eta1
    sll_real64 :: eta2
    sll_real64 :: eta1_min
    sll_real64 :: eta2_min
    sll_real64 :: delta1
    sll_real64 :: delta2

    eta1_min = transf%mesh%eta1_min
    eta2_min = transf%mesh%eta2_min
    delta1   = transf%mesh%delta_eta1
    delta2   = transf%mesh%delta_eta2
    eta1     = eta1_min + (real(i-1,f64)+0.5_f64)*delta1 
    eta2     = eta2_min + (real(j-1,f64)+0.5_f64)*delta2
    var      = transf%x1_func( eta1, eta2, transf%params )
  end function x1_cell_analytic

  function x2_cell_analytic( transf, i, j ) result(var)
#ifdef STDF95
    type(sll_coordinate_transformation_2d_analytic) :: transf
#else
    class(sll_coordinate_transformation_2d_analytic) :: transf
#endif
    sll_real64            :: var
    sll_int32, intent(in) :: i
    sll_int32, intent(in) :: j
    sll_real64 :: eta1
    sll_real64 :: eta2
    sll_real64 :: eta1_min
    sll_real64 :: eta2_min
    sll_real64 :: delta1
    sll_real64 :: delta2

    eta1_min = transf%mesh%eta1_min
    eta2_min = transf%mesh%eta2_min
    delta1   = transf%mesh%delta_eta1
    delta2   = transf%mesh%delta_eta2
    eta1     = eta1_min + (real(i-1,f64)+0.5_f64)*delta1 
    eta2     = eta2_min + (real(j-1,f64)+0.5_f64)*delta2
    var      = transf%x2_func( eta1, eta2, transf%params )
  end function x2_cell_analytic

  function jacobian_2d_cell_analytic( transf, i, j ) result(val)
#ifdef STDF95
    type(sll_coordinate_transformation_2d_analytic) :: transf
#else
    class(sll_coordinate_transformation_2d_analytic) :: transf
#endif
    sll_real64            :: val
    sll_int32, intent(in) :: i
    sll_int32, intent(in) :: j
    sll_real64 :: eta1
    sll_real64 :: eta2
    sll_real64 :: eta1_min
    sll_real64 :: eta2_min
    sll_real64 :: delta1
    sll_real64 :: delta2
    sll_real64 :: j11
    sll_real64 :: j12
    sll_real64 :: j21
    sll_real64 :: j22

    eta1_min = transf%mesh%eta1_min
    eta2_min = transf%mesh%eta2_min
    delta1   = transf%mesh%delta_eta1
    delta2   = transf%mesh%delta_eta2
    eta1     = eta1_min + (real(i-1,f64)+0.5_f64)*delta1 
    eta2     = eta2_min + (real(j-1,f64)+0.5_f64)*delta2
    j11 = (transf%j_matrix(1,1)%f( eta1, eta2, transf%params ))
    j12 = (transf%j_matrix(1,2)%f( eta1, eta2, transf%params ))
    j21 = (transf%j_matrix(2,1)%f( eta1, eta2, transf%params ))
    j22 = (transf%j_matrix(2,2)%f( eta1, eta2, transf%params ))
    ! For debugging:
    !    print *, 'jacobian_2d_analytic: '
    !    print *, j11, j12
    !    print *, j21, j22
    val = j11*j22 - j12*j21
  end function jacobian_2d_cell_analytic

  function jacobian_node_analytic( transf, i, j )
#ifdef STDF95
    type(sll_coordinate_transformation_2d_analytic)   :: transf
#else
    class(sll_coordinate_transformation_2d_analytic)   :: transf
#endif
    sll_real64              :: jacobian_node_analytic
    sll_int32, intent(in)   :: i
    sll_int32, intent(in)   :: j
    sll_int32 :: num_pts_1
    sll_int32 :: num_pts_2
    sll_real64 :: eta1
    sll_real64 :: eta2
    sll_real64 :: eta1_min
    sll_real64 :: eta2_min
    sll_real64 :: delta1
    sll_real64 :: delta2
    sll_real64 :: j11
    sll_real64 :: j12
    sll_real64 :: j21
    sll_real64 :: j22

    num_pts_1 = transf%mesh%num_cells1 + 1
    num_pts_2 = transf%mesh%num_cells2 + 1
    SLL_ASSERT( (i .ge. 1) .and. (i .le. num_pts_1) )
    SLL_ASSERT( (j .ge. 1) .and. (j .le. num_pts_2) )

    eta1_min = transf%mesh%eta1_min
    eta2_min = transf%mesh%eta2_min
    delta1   = transf%mesh%delta_eta1
    delta2   = transf%mesh%delta_eta2
    eta1     = eta1_min + real(i-1,f64)*delta1 
    eta2     = eta2_min + real(j-1,f64)*delta2
    j11 = (transf%j_matrix(1,1)%f( eta1, eta2, transf%params ))
    j12 = (transf%j_matrix(1,2)%f( eta1, eta2, transf%params ))
    j21 = (transf%j_matrix(2,1)%f( eta1, eta2, transf%params ))
    j22 = (transf%j_matrix(2,2)%f( eta1, eta2, transf%params ))
    ! For debugging:
    !    print *, 'jacobian_2d_analytic: '
    !    print *, j11, j12
    !    print *, j21, j22
    jacobian_node_analytic = j11*j22 - j12*j21
  end function jacobian_node_analytic


  subroutine write_to_file_2d_analytic(transf,output_format)
    class(sll_coordinate_transformation_2d_analytic) :: transf
    sll_int32, optional :: output_format 
    sll_int32           :: local_format 
    sll_real64, dimension(:,:), pointer :: x1mesh
    sll_real64, dimension(:,:), pointer :: x2mesh
    sll_int32  :: i1
    sll_int32  :: i2
    sll_real64 :: eta1
    sll_real64 :: eta2
    sll_int32  :: ierr
    sll_int32  :: file_id
    sll_int32  :: nc_eta1
    sll_int32  :: nc_eta2

    nc_eta1 = transf%mesh%num_cells1
    nc_eta2 = transf%mesh%num_cells2

    if (.not. present(output_format)) then
       local_format = SLL_IO_XDMF
    else
       local_format = output_format
    end if

    if ( .not. transf%written ) then
       if (local_format == SLL_IO_XDMF) then
          SLL_ALLOCATE(x1mesh(nc_eta1+1,nc_eta2+1), ierr)
          SLL_ALLOCATE(x2mesh(nc_eta1+1,nc_eta2+1), ierr)
          eta1 = transf%mesh%eta1_min
          do i1=1, nc_eta1+1
             eta2 = transf%mesh%eta2_min
             do i2=1, nc_eta2+1
                x1mesh(i1,i2) = x1_node_analytic(transf,i1,i2)
                x2mesh(i1,i2) = x2_node_analytic(transf,i1,i2)
                eta2 = eta2 + transf%mesh%delta_eta2 
             end do
             eta1 = eta1 + transf%mesh%delta_eta1
          end do
       
          call sll_xdmf_open(trim(transf%label)//".xmf",transf%label, &
               nc_eta1+1,nc_eta2+1,file_id,ierr)
          call sll_xdmf_write_array(transf%label,x1mesh,"x1",ierr)
          call sll_xdmf_write_array(transf%label,x2mesh,"x2",ierr)
          call sll_xdmf_close(file_id,ierr)

       else if (local_format == SLL_IO_MTV) then

          SLL_ALLOCATE(x1mesh(nc_eta1+1,nc_eta2+1), ierr)
          SLL_ALLOCATE(x2mesh(nc_eta1+1,nc_eta2+1), ierr)

          do i1=1, nc_eta1+1
             do i2=1, nc_eta2+1
                x1mesh(i1,i2) = x1_node_analytic(transf,i1,i2)
                x2mesh(i1,i2) = x2_node_analytic(transf,i1,i2)
             end do
          end do
       
          call sll_plotmtv_write( nc_eta1+1,nc_eta2+1, &
                                  x1mesh, x2mesh, trim(transf%label),ierr)

       else
          print*, 'Not recognized format to write this mesh'
          stop
       end if
    else
       print*,' Warning, you have already written the mesh '
    end if
    transf%written = .true.
    if( associated(x1mesh) ) then
       SLL_DEALLOCATE(x1mesh, ierr)
    end if
    if( associated(x2mesh) ) then
       SLL_DEALLOCATE(x2mesh, ierr)
    end if
  end subroutine

  subroutine read_from_file_2d_analytic( transf, filename )
    class(sll_coordinate_transformation_2d_analytic), intent(inout) :: transf
    character(len=*), intent(in) :: filename
    print *, filename
    print *, 'read_from_file_2d_analytic: not yet implemented'
    call sll_display(transf%mesh)
    ! here we could put a case select to choose which analytic transformation
    ! we would like to use.
  end subroutine read_from_file_2d_analytic


  !**************************************************************************
  !
  !        Functions for the discrete general transformation
  !
  !**************************************************************************

  function get_logical_mesh_discrete( transf ) result(res)
    class(sll_coordinate_transformation_2d_discrete), intent(in) :: transf
    type(sll_logical_mesh_2d), pointer :: res
    res => transf%mesh
  end function get_logical_mesh_discrete


  function x1_node_discrete( transf, i, j ) result(val)
#ifdef STDF95
    type(sll_coordinate_transformation_2d_discrete) :: transf
#else
    class(sll_coordinate_transformation_2d_discrete) :: transf
#endif
    sll_real64             :: val
    sll_int32, intent(in) :: i
    sll_int32, intent(in) :: j
    val = transf%x1_node(i,j)
  end function x1_node_discrete

  function x2_node_discrete( transf, i, j ) result(val)
#ifdef STDF95
    type(sll_coordinate_transformation_2d_discrete) :: transf
#else
    class(sll_coordinate_transformation_2d_discrete) :: transf
#endif
    sll_real64             :: val
    sll_int32, intent(in) :: i
    sll_int32, intent(in) :: j
    val = transf%x2_node(i,j)
  end function x2_node_discrete

  function x1_cell_discrete( transf, i, j ) result(var)
#ifdef STDF95
    type(sll_coordinate_transformation_2d_discrete) :: transf
#else
    class(sll_coordinate_transformation_2d_discrete) :: transf
#endif
    sll_real64                         :: var
    sll_int32, intent(in)              :: i
    sll_int32, intent(in)              :: j
    var = transf%x1_cell(i,j)
  end function x1_cell_discrete

  function x2_cell_discrete( transf, i, j ) result(var)
#ifdef STDF95
    type(sll_coordinate_transformation_2d_discrete) :: transf
#else
    class(sll_coordinate_transformation_2d_discrete) :: transf
#endif
    sll_real64                         :: var
    sll_int32, intent(in)              :: i
    sll_int32, intent(in)              :: j
    var = transf%x2_cell(i,j)
  end function x2_cell_discrete

  function x1_discrete( transf, eta1, eta2 ) result(val)
#ifdef STDF95
    type(sll_coordinate_transformation_2d_discrete) :: transf
#else
    class(sll_coordinate_transformation_2d_discrete) :: transf
#endif
    sll_real64             :: val
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
#ifdef STDF95
    val = cubic_spline_interpolate_value(transf%x1_interp,eta1, eta2)
#else
    val = transf%x1_interp%interpolate_value(eta1, eta2)
#endif
  end function x1_discrete

  function x2_discrete( transf, eta1, eta2 ) result(val)
#ifdef STDF95
    type(sll_coordinate_transformation_2d_discrete) :: transf
#else
    class(sll_coordinate_transformation_2d_discrete) :: transf
#endif
    sll_real64             :: val
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
#ifdef STDF95
    val = cubic_spline_interpolate_value(transf%x2_interp, eta1, eta2)
#else
    val = transf%x2_interp%interpolate_value(eta1, eta2)
#endif
  end function x2_discrete

  function jacobian_2d_discrete( transf, eta1, eta2 ) result(jac)
#ifdef STDF95
    type(sll_coordinate_transformation_2d_discrete) :: transf
#else
    class(sll_coordinate_transformation_2d_discrete) :: transf
#endif
    sll_real64             :: jac
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    sll_real64             :: j11
    sll_real64             :: j12
    sll_real64             :: j21
    sll_real64             :: j22
#ifdef STDF95
    j11 = cubic_spline_interpolate_derivative_eta1( transf%x1_interp, eta1, eta2 )
    j12 = cubic_spline_interpolate_derivative_eta2( transf%x1_interp, eta1, eta2 )
    j21 = cubic_spline_interpolate_derivative_eta1( transf%x1_interp, eta1, eta2 )
    j22 = cubic_spline_interpolate_derivative_eta2( transf%x1_interp, eta1, eta2 )
#else
    j11 = transf%x1_interp%interpolate_derivative_eta1( eta1, eta2 )
    j12 = transf%x1_interp%interpolate_derivative_eta2( eta1, eta2 )
    j21 = transf%x2_interp%interpolate_derivative_eta1( eta1, eta2 )
    j22 = transf%x2_interp%interpolate_derivative_eta2( eta1, eta2 )
#endif
    ! For debugging:
    !    print *, 'jacobian_2D_discrete: '
    !    print *, j11, j12
    !    print *, j21, j22
    jac = j11*j22 - j12*j21
  end function jacobian_2d_discrete

  function jacobian_2d_cell_discrete( transf, i, j ) result(var)
#ifdef STDF95
    type(sll_coordinate_transformation_2d_discrete)  :: transf
#else
    class(sll_coordinate_transformation_2d_discrete) :: transf
#endif
    sll_real64                         :: var
    sll_int32, intent(in)              :: i
    sll_int32, intent(in)              :: j
    var = transf%jacobians_c(i,j)
  end function jacobian_2d_cell_discrete

  function jacobian_matrix_2d_discrete( transf, eta1, eta2 )
#ifdef STDF95
    type(sll_coordinate_transformation_2d_discrete) :: transf
#else
    class(sll_coordinate_transformation_2d_discrete) :: transf
#endif
    sll_real64, dimension(1:2,1:2)     :: jacobian_matrix_2d_discrete
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    sll_real64             :: j11
    sll_real64             :: j12
    sll_real64             :: j21
    sll_real64             :: j22
#ifdef STDF95
    j11 = cubic_spline_interpolate_derivative_eta1( transf%x1_interp, eta1, eta2 )
    j12 = cubic_spline_interpolate_derivative_eta2( transf%x1_interp, eta1, eta2 )
    j21 = cubic_spline_interpolate_derivative_eta1( transf%x2_interp, eta1, eta2 )
    j22 = cubic_spline_interpolate_derivative_eta2( transf%x2_interp, eta1, eta2 )
#else
    j11 = transf%x1_interp%interpolate_derivative_eta1( eta1, eta2 )
    j12 = transf%x1_interp%interpolate_derivative_eta2( eta1, eta2 )
    j21 = transf%x2_interp%interpolate_derivative_eta1( eta1, eta2 )
    j22 = transf%x2_interp%interpolate_derivative_eta2( eta1, eta2 )
#endif
    ! For debugging:
    !    print *, 'jacobian_2D_discrete: '
    !    print *, j11, j12
    !    print *, j21, j22
    jacobian_matrix_2d_discrete(1,1) = j11
    jacobian_matrix_2d_discrete(1,2) = j12
    jacobian_matrix_2d_discrete(2,1) = j21
    jacobian_matrix_2d_discrete(2,2) = j22
  end function jacobian_matrix_2d_discrete

  function inverse_jacobian_matrix_2d_discrete( transf, eta1, eta2 )
#ifdef STDF95
    type(sll_coordinate_transformation_2d_discrete)  :: transf
#else
    class(sll_coordinate_transformation_2d_discrete) :: transf
#endif
    sll_real64, dimension(1:2,1:2)     :: inverse_jacobian_matrix_2d_discrete
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    sll_real64             :: inv_j11
    sll_real64             :: inv_j12
    sll_real64             :: inv_j21
    sll_real64             :: inv_j22
    sll_real64             :: r_jac ! reciprocal of the jacobian
#ifdef STDF95
    r_jac = jacobian_2d_discrete( transf, eta1, eta2 )
    inv_j11 = &
         cubic_spline_interpolate_derivative_eta1( transf%x1_interp, eta1, eta2 )
    inv_j12 = &
         cubic_spline_interpolate_derivative_eta2( transf%x1_interp, eta1, eta2 )
    inv_j21 = &
         cubic_spline_interpolate_derivative_eta1( transf%x2_interp, eta1, eta2 )
    inv_j22 = &
         cubic_spline_interpolate_derivative_eta2( transf%x2_interp, eta1, eta2 )
#else
    r_jac = 1.0_f64/transf%jacobian( eta1, eta2 )
    inv_j11 = transf%x1_interp%interpolate_derivative_eta1( eta1, eta2 )
    inv_j12 = transf%x1_interp%interpolate_derivative_eta2( eta1, eta2 )
    inv_j21 = transf%x2_interp%interpolate_derivative_eta1( eta1, eta2 )
    inv_j22 = transf%x2_interp%interpolate_derivative_eta2( eta1, eta2 )
#endif
    ! For debugging:
    !    print *, 'jacobian_2D_discrete: '
    !    print *, j11, j12
    !    print *, j21, j22
    inverse_jacobian_matrix_2d_discrete(1,1) =  inv_j22*r_jac
    inverse_jacobian_matrix_2d_discrete(1,2) = -inv_j12*r_jac
    inverse_jacobian_matrix_2d_discrete(2,1) = -inv_j21*r_jac
    inverse_jacobian_matrix_2d_discrete(2,2) =  inv_j11*r_jac
  end function inverse_jacobian_matrix_2d_discrete

  function new_coordinate_transformation_2d_discrete( &
       mesh_2d,        &
       label,          &
       x1_interpolator,&
       x2_interpolator,&
       jacobians_n_interpolator, &
       x1_node,        &
       x2_node,        &
       jacobians_node, &
       x1_cell, &
       x2_cell, &
       jacobians_cell )

    ! INPUT VARIABLES
    type(sll_logical_mesh_2d), pointer    :: mesh_2d
    character(len=*)         , intent(in) :: label

#ifdef STDF95
    ! no this should not be compiled under f95, the object would not have
    ! the same functionality
    type(cubic_spline_2d_interpolator), target  :: x1_interpolator
    type(cubic_spline_2d_interpolator), target  :: x2_interpolator
    type(cubic_spline_2d_interpolator), target  :: jacobians_n_interpolator

#else
    class(sll_interpolator_2d_base), target  :: x1_interpolator
    class(sll_interpolator_2d_base), target  :: x2_interpolator
    class(sll_interpolator_2d_base), target  :: jacobians_n_interpolator
#endif 
    sll_real64, dimension(:,:), intent(in), optional :: x1_node
    sll_real64, dimension(:,:), intent(in), optional :: x2_node
    sll_real64, dimension(:,:), intent(in), optional :: jacobians_node
    sll_real64, dimension(:,:), intent(in), optional :: x1_cell
    sll_real64, dimension(:,:), intent(in), optional :: x2_cell
    sll_real64, dimension(:,:), intent(in), optional :: jacobians_cell

    ! LOCAL VARIABLES
    type(sll_coordinate_transformation_2d_discrete), pointer :: &
         new_coordinate_transformation_2d_discrete
    sll_int32 :: ierr
    
    SLL_ALLOCATE(new_coordinate_transformation_2d_discrete, ierr)
    call initialize_coord_transf_2d_discrete( &
         new_coordinate_transformation_2d_discrete, &
         mesh_2d,           &
         label,            &
         x1_interpolator, &
         x2_interpolator, &
         jacobians_n_interpolator, &
         x1_node,        &
         x2_node,        &
         jacobians_node, &
         x1_cell, &
         x2_cell, &
         jacobians_cell )
  end function new_coordinate_transformation_2d_discrete

  subroutine initialize_coord_transf_2d_discrete( &
    transf,            &
    mesh_2d,           &
    label,            &
    x1_interpolator, &
    x2_interpolator, &
    jacobians_n_interpolator, &
    x1_node,        &
    x2_node,        &
    jacobians_node, &
    x1_cell, &
    x2_cell, &
    jacobians_cell )

#ifdef STDF95
    type(sll_coordinate_transformation_2d_discrete)     :: transf
#else
    class(sll_coordinate_transformation_2d_discrete)    :: transf
#endif
    type(sll_logical_mesh_2d), pointer   :: mesh_2d
    character(len=*), intent(in)         :: label

#ifdef STDF95
    ! no this should not be compiled under f95, the object would not have
    ! the same functionality
    type(cubic_spline_2d_interpolator), target  :: x1_interpolator
    type(cubic_spline_2d_interpolator), target  :: x2_interpolator
    type(cubic_spline_2d_interpolator), target  :: jacobians_n_interpolator
#else
    class(sll_interpolator_2d_base), target  :: x1_interpolator
    class(sll_interpolator_2d_base), target  :: x2_interpolator
    class(sll_interpolator_2d_base), target :: jacobians_n_interpolator
#endif 
    sll_real64, dimension(:,:), intent(in), optional :: x1_node
    sll_real64, dimension(:,:), intent(in), optional :: x2_node
    sll_real64, dimension(:,:), intent(in), optional :: jacobians_node
    sll_real64, dimension(:,:), intent(in), optional :: jacobians_cell
    sll_real64, dimension(:,:), intent(in), optional :: x1_cell
    sll_real64, dimension(:,:), intent(in), optional :: x2_cell

    sll_real64 :: eta_1
    sll_real64 :: eta_2
    sll_real64 :: eta_1_min
    sll_real64 :: eta_2_min
    sll_real64 :: delta_eta_1
    sll_real64 :: delta_eta_2
    sll_real64 :: jacobian_val
    sll_int32  :: i
    sll_int32  :: j
    sll_int32  :: ierr
    sll_int32  :: npts1
    sll_int32  :: npts2
    logical    :: x1n
    logical    :: x2n
    logical    :: x1c
    logical    :: x2c
    logical    :: jc
    logical    :: jn

    transf%mesh => mesh_2d
    transf%label = trim(label)
    x1n = present(x1_node)
    x2n = present(x2_node)
    x1c = present(x1_cell)
    x2c = present(x2_cell)
    jc  = present(jacobians_cell)
    jn  = present(jacobians_node)
    npts1 = mesh_2d%num_cells1 + 1
    npts2 = mesh_2d%num_cells2 + 1

    ! Check argument consistency
    ! DISCRETE_MAPs require only some of the parameters. If the mapping is
    ! defined from the nodes of the logical (eta1, eta2) mesh to the nodes
    ! of the physical mesh (x1,x2), then either:
    ! - the node arrays are required: jacobians_node, x1_node and x2_node. Or
    ! - the x1_interpolator and x2_interpolator must contain already the
    !   coefficient information that would permit the calculation of the
    !   x1 and x2 points.
    !
    ! If the transformation is done on the points at the center of the cells
    ! then these parameters are also required: 
    ! jacobians_cell, x1_cell, x2_cell.
    ! node and cell values are not mutually exclusive, thus all 6 parameters 
    ! can be provided in the discrete case. It is up to the caller to make
    ! sure that the data set is consistent.

    ! 1. Check that the discrete representation of the transformation is
    !    consistent with the size of the 2D array.

    if( (      x1n  .and. (.not.x2n) ) .or. &
        ((.not.x1n) .and.       x2n) ) then
       print *,'ERROR, initialize_coordinate_transformation_2d_discrete():', &
            'for the moment, this function does not support specifying ', &
            'transformation only with one of the node arrays x1_node or ', &
            'x2_node. Either pass both, or none, but with the ', &
            'corresponding interpolators having their coefficients ', &
            'already set.'
       STOP
    end if

    if( x1n ) then
       if( size(x1_node,1) .lt. npts1 ) then
          print *,'ERROR, initialize_coordinate_transformation_2d_discrete()',&
               ' the size of the x1_node arrays is ', &
               'inconsistent with the number of points declared, ', &
               'in the logical mesh.'
          STOP 
       end if
    end if

    if( x2n ) then
       if( size(x1_node,2) .lt. npts2 ) then
          print *, 'ERROR, initialize_coordinate_transformation_2d_discrete()',&
               ' the size of the x2_node arrays is ', &
               'inconsistent with the number of points declared, ', &
               'in the logical mesh.'
          STOP
       end if
    end if

    if( jn .eqv. .true. ) then
       if( &
          (size(jacobians_node,1) .lt. npts1 - 1 ) .or. &
          (size(jacobians_node,2) .lt. npts2 - 1 ) ) then
          print *, 'ERROR, initialize_coordinate_transformation_2d_discrete()',&
               ': the size of the jacobians_node array is ', &
               'inconsistent with the number of points declared, ', &
               'npts1 or npts2.'
          STOP
       end if
    end if

    if( jc .eqv. .true. ) then
       if( &
          (size(jacobians_cell,1) .lt. npts1 - 1 ) .or. &
          (size(jacobians_cell,2) .lt. npts2 - 1 ) ) then
          print *, 'ERROR, initialize_coordinate_transformation_2d_discrete()',&
               ': the size of the jacobians_cell arrays is ', &
               'inconsistent with the number of points declared, ', &
               'npts1 or npts2.'
          STOP
       end if
    end if

    transf%x1_interp => x1_interpolator
    transf%x2_interp => x2_interpolator

    ! Allocate the arrays for precomputed jacobians.
    SLL_ALLOCATE(transf%jacobians_n(npts1,npts2), ierr)
    SLL_ALLOCATE(transf%jacobians_c(npts1-1, npts2-1), ierr)

    ! Allocation for x1 and x2 at nodes
    SLL_ALLOCATE(transf%x1_node(npts1,npts2), ierr)
    SLL_ALLOCATE(transf%x2_node(npts1,npts2), ierr)

    ! Allocation for x1 and x2 at cells
    SLL_ALLOCATE(transf%x1_cell(npts1-1,npts2-1), ierr)
    SLL_ALLOCATE(transf%x2_cell(npts1-1,npts2-1), ierr)

    ! initialize the local arrays. Note that since the map has its
    ! own copies, it owns this information locally and will destroy
    ! this information when the object is deleted. The caller is
    ! thus responsible for deallocating the arrays that were passed as
    ! arguments.

    eta_1_min   = mesh_2d%eta1_min
    eta_2_min   = mesh_2d%eta2_min
    delta_eta_1 = mesh_2d%delta_eta1
    delta_eta_2 = mesh_2d%delta_eta2

    if( x1n .and. x2n ) then
       do j=1, npts2
          do i=1, npts1
             transf%x1_node(i,j) = x1_node(i,j)
             transf%x2_node(i,j) = x2_node(i,j)
          end do
       end do
    else
       if(x1_interpolator%coefficients_are_set() .eqv. .false.) then
          print *, 'ERROR, initialize_coordinate_transformation_2d_discrete()',&
               ': the x1_node array was not passed and the corresponding ', &
               'interpolator has no initialized coefficients. Exiting...'
          STOP
       endif
       if(x2_interpolator%coefficients_are_set() .eqv. .false.) then
          print *, 'ERROR, initialize_coordinate_transformation_2d_discrete()',&
               ': the x2_node array was not passed and the corresponding ', &
               'interpolator has no initialized coefficients. Exiting...'
          STOP
       endif
       ! now initialize the arrays starting from the interpolator information
       ! and the logical mesh information.
       do j=0, npts2 - 1
          eta_2 = eta_2_min + real(j,f64)*delta_eta_2          
          do i=0, npts1 - 1
             eta_1 = eta_1_min + real(i,f64)*delta_eta_1
             transf%x1_node(i+1,j+1) = &
                  x1_interpolator%interpolate_value(eta_1,eta_2)
             transf%x2_node(i+1,j+1) = &
                  x2_interpolator%interpolate_value(eta_1,eta_2)
          end do
       end do
    end if

    ! Compute the spline coefficients
!!$#ifdef STDF95
!!$    call cubic_spline_compute_interpolants( x1_interpolator, transf%x1_node )
!!$    call cubic_spline_compute_interpolants( x2_interpolator, transf%x2_node )
!!$#else
    if( x1n .and. (x1_interpolator%coefficients_are_set() .eqv. .false.) ) then
       call x1_interpolator%compute_interpolants( transf%x1_node )
    end if
    if( x2n .and. (x2_interpolator%coefficients_are_set() .eqv. .false.) ) then
       call x2_interpolator%compute_interpolants( transf%x2_node )
    end if
!!$#endif


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
             transf%jacobians_n(i,j) = jacobians_node(i,j)
          end do
       end do
    else
       ! Fill the jacobian values at the nodes calculated from the splines
       do j=0, npts2 - 1
          eta_2 = eta_2_min + real(j,f64)*delta_eta_2          
          do i=0, npts1 - 1
             eta_1 = eta_1_min + real(i,f64)*delta_eta_1
#ifdef STDF95
             jacobian_val = jacobian_2d_discrete(transf,eta_1,eta_2)
#else
             jacobian_val = transf%jacobian(eta_1,eta_2)
#endif
             transf%jacobians_n(i+1,j+1) = jacobian_val
          end do
       end do
    end if

    ! copy the cell-based transformation arrays if available
    if( (x1c .and. x2c) .eqv. .true. ) then
       SLL_ALLOCATE(transf%x1_cell(npts1-1, npts2-1), ierr)
       SLL_ALLOCATE(transf%x2_cell(npts1-1, npts2-1), ierr)
       do j=1, npts2 - 1
          do i=1, npts1 - 1
             transf%x1_cell(i,j) = x1_cell(i,j)
             transf%x2_cell(i,j) = x2_cell(i,j)
          end do
       end do
    end if
    ! copy the cell-based jacobians if available
    if( jc .eqv. .true. ) then
       do j=1, npts2 - 1
          do i=1, npts1 - 1
             transf%jacobians_c(i,j) = jacobians_cell(i,j)
          end do
       end do
    else ! if cell-based jacobians are not available, compute them.
    ! Fill the values at the mid-point of the cells
       do j=0, npts2 - 2
          eta_2 = eta_2_min + delta_eta_2*(real(j,f64) + 0.5_f64)
          do i=0, npts1 - 2
             ! it is very bad practice to invoke the mesh methods while
             ! we are not even done initializing the mesh object...
             eta_1 = eta_1_min + delta_eta_1*(real(i,f64) + 0.5_f64)
#ifdef STDF95
             transf%x1_cell(i+1,j+1)     = x1_discrete(transf, eta_1, eta_2)
             transf%x2_cell(i+1,j+1)     = x2_discrete(transf, eta_1, eta_2)
             transf%jacobians_c(i+1,j+1) = jacobian_2d_discrete(transf, eta_1,eta_2)
#else
             transf%x1_cell(i+1,j+1)     = transf%x1(eta_1, eta_2)
             transf%x2_cell(i+1,j+1)     = transf%x2(eta_1, eta_2)
             transf%jacobians_c(i+1,j+1) = transf%jacobian(eta_1,eta_2)
#endif
          end do
       end do
    end if
  end subroutine initialize_coord_transf_2d_discrete


  function transf_2d_jacobian_node_discrete( transf, i, j )
#ifdef STDF95
    type(sll_coordinate_transformation_2d_discrete)   :: transf
#else
    class(sll_coordinate_transformation_2d_discrete)   :: transf
#endif
    sll_real64              :: transf_2d_jacobian_node_discrete
    sll_int32, intent(in)   :: i
    sll_int32, intent(in)   :: j
    sll_int32 :: num_pts_1
    sll_int32 :: num_pts_2
    num_pts_1 = transf%mesh%num_cells1 + 1
    num_pts_2 = transf%mesh%num_cells2 + 1
    SLL_ASSERT( (i .ge. 1) .and. (i .le. num_pts_1) )
    SLL_ASSERT( (j .ge. 1) .and. (j .le. num_pts_2) )
    transf_2d_jacobian_node_discrete = transf%jacobians_n(i,j)
  end function transf_2d_jacobian_node_discrete

  subroutine write_to_file_2d_discrete(transf,output_format)
    class(sll_coordinate_transformation_2d_discrete) :: transf
    sll_int32, optional :: output_format 
    sll_int32           :: local_format 
    sll_real64, dimension(:,:), pointer :: x1mesh
    sll_real64, dimension(:,:), pointer :: x2mesh
    sll_int32  :: i1
    sll_int32  :: i2
    sll_real64 :: eta1
    sll_real64 :: eta2
    sll_int32  :: ierr
    sll_int32  :: file_id
    sll_int32  :: npts_eta1
    sll_int32  :: npts_eta2
    sll_real64 :: eta1_min
    sll_real64 :: eta2_min
    sll_real64 :: delta_eta1
    sll_real64 :: delta_eta2

    npts_eta1  = transf%mesh%num_cells1 +1
    npts_eta2  = transf%mesh%num_cells2 +1
    eta1_min   = transf%mesh%eta1_min
    eta2_min   = transf%mesh%eta1_min
    delta_eta1 = transf%mesh%delta_eta1
    delta_eta2 = transf%mesh%delta_eta2

    if (.not. present(output_format)) then
       local_format = SLL_IO_XDMF
    else
       local_format = output_format
    end if

    if ( .not. transf%written ) then

       if (local_format == SLL_IO_XDMF) then
          SLL_ALLOCATE(x1mesh(npts_eta1,npts_eta2), ierr)
          SLL_ALLOCATE(x2mesh(npts_eta1,npts_eta2), ierr)
          eta1 = eta1_min
          do i1=1, npts_eta1
             eta2 = eta2_min
             do i2=1, npts_eta2
                x1mesh(i1,i2) = x1_node_discrete(transf,i1,i2)
                x2mesh(i1,i2) = x2_node_discrete(transf,i1,i2)
                eta2 = eta2 + delta_eta2 
             end do
             eta1 = eta1 + delta_eta1
          end do
       
          call sll_xdmf_open(trim(transf%label)//".xmf",transf%label, &
               npts_eta1,npts_eta2,file_id,ierr)
          call sll_xdmf_write_array(transf%label,x1mesh,"x1",ierr)
          call sll_xdmf_write_array(transf%label,x2mesh,"x2",ierr)
          call sll_xdmf_close(file_id,ierr)

       else if (local_format == SLL_IO_MTV) then

          SLL_ALLOCATE(x1mesh(npts_eta1,npts_eta2), ierr)
          SLL_ALLOCATE(x2mesh(npts_eta1,npts_eta2), ierr)

          do i1=1, npts_eta1
             do i2=1, npts_eta2
                x1mesh(i1,i2) = x1_node_discrete(transf,i1,i2)
                x2mesh(i1,i2) = x2_node_discrete(transf,i1,i2)
             end do
          end do
       
          call sll_plotmtv_write( npts_eta1,npts_eta2, &
                                  x1mesh, x2mesh, trim(transf%label),ierr)

       else
          print*, 'Not recognized format to write this mesh'
          stop
       end if
    else
       print*,' Warning, you have already written the mesh '
    end if

    transf%written = .true.
  end subroutine

  subroutine delete_transformation_2d_discrete( transf )
    class(sll_coordinate_transformation_2d_discrete), intent(inout) :: transf
!!$    sll_int32 :: ierr
!!$    SLL_DEALLOCATE( transf%j_matrix, ierr )
    transf%label = ""
    transf%written = .false.
    nullify( transf%x1_node )
    nullify( transf%x2_node )
    nullify( transf%x1_cell )
    nullify( transf%x2_cell )
    nullify( transf%jacobians_n )
    nullify( transf%jacobians_c )
    
    !call delete( transf%x1_interp)
    !call delete( transf%x2_interp)
    call delete( transf%mesh)
    ! Fix: there is a dependency problem where these pointers are not recognized
    ! during the linking step. A similar nullification of an abstract class
    ! pointer is carried out in the fields_2d_alternative type without problems.
!    transf%x1_interp => null() this gives a different message.
!!$    nullify( transf%x1_interp )
!!$    nullify( transf%x2_interp )
  end subroutine delete_transformation_2d_discrete

  ! What do we need to initialize fully a discrete coordinate transformation?
  ! - logical mesh
  ! - label
  ! - array with x1 node positions
  ! - array with x2 node positions
  ! - array with x1 at cell-center positions
  ! - array with x2 at cell-center positions
  ! - array with jacobians at nodes
  ! - array with jacobians at cell-centers
  ! - interpolator 2d for x1
  ! - interpolator 2d for x2
  ! - the file used to initialize the transformation should allow us to 
  !   initialize all this data. This routine has special rights in that it
  !   is allowed to allocate and initialize a logical mesh and the two
  !   interpolators.

  ! - Issues to decide:
  ! - Will there be a single file format? Or multiple file formats? 
  !   The transformation can be specified by two 2D arrays of points, or by
  !   the spline coefficients...
  ! - The BC information is not inside the files we are currently considering,
  !   so this should be included.
 
  subroutine read_from_file_2d_discrete( transf, filename )
    use sll_arbitrary_degree_spline_interpolator_2d_module
    class(sll_coordinate_transformation_2d_discrete), intent(inout) :: transf
    character(len=*), intent(in) :: filename
    intrinsic :: trim
    character(len=256) :: filename_local
    sll_int32 :: IO_stat
    sll_int32 :: input_file_id
    sll_int32 :: ierr
    sll_int32 :: spline_deg1
    sll_int32 :: spline_deg2
    sll_int32 :: num_pts1
    sll_int32 :: num_pts2
    sll_int32 :: is_rational
    character(len=256) :: label
    sll_real64, dimension(:), allocatable :: knots1
    sll_real64, dimension(:), allocatable :: knots2
    sll_real64, dimension(:), allocatable :: control_pts1
    sll_real64, dimension(:), allocatable :: control_pts2
    sll_real64, dimension(:), allocatable :: weights
    sll_real64, dimension(:,:), allocatable :: control_pts1_2d
    sll_real64, dimension(:,:), allocatable :: control_pts2_2d
    sll_real64, dimension(:,:), allocatable :: weights_2d
    sll_real64 :: eta1_min
    sll_real64 :: eta1_max
    sll_real64 :: eta2_min
    sll_real64 :: eta2_max
    sll_int32  :: bc_left
    sll_int32  :: bc_right
    sll_int32  :: bc_bottom
    sll_int32  :: bc_top
!    sll_real64, dimension(:,:), allocatable :: nodes1
!    sll_real64, dimension(:,:), allocatable :: nodes2
    sll_int32  :: number_cells1,number_cells2
    sll_int32 :: sz_knots1,sz_knots2
    class(arb_deg_2d_interpolator), pointer :: interp2d_1
    class(arb_deg_2d_interpolator), pointer :: interp2d_2
    class(arb_deg_2d_interpolator), pointer :: interp2d_jac
    type(sll_logical_mesh_2d), pointer      :: mesh_2d
   
    namelist /transf_label/  label
    namelist /degree/   spline_deg1, spline_deg2
    namelist /shape/    num_pts1, num_pts2 ! it is not the number of points but the number of coeff sdpline in each direction !!
    namelist /rational/ is_rational
    namelist /knots_1/   knots1
    namelist /knots_2/   knots2
    namelist /control_points/ control_pts1, control_pts2
    namelist /pt_weights/  weights
    namelist /logical_mesh_2d/ number_cells1,number_cells2
  
    if(len(filename) >= 256) then
       print *, 'ERROR, read_coefficients_from_file => ',&
            'read_from_file_discrete():',&
            'filenames longer than 256 characters are not allowed.'
       STOP
    end if
    filename_local = trim(filename)

    ! get a new identifier for the file.
    call sll_new_file_id( input_file_id, ierr )
    if( ierr .ne. 0 ) then
       print *, 'ERROR while trying to obtain an unique identifier for file ',&
            filename, '. Called from read_coeffs_ad2d().'
       stop
    end if
    open(unit=input_file_id, file=filename_local, STATUS="OLD", IOStat=IO_stat)
    if( IO_Stat .ne. 0 ) then
       print *, 'ERROR while opening file ',filename, &
            '. Called from read_coeffs_ad2d().'
       stop
    end if

    read( input_file_id, transf_label )
    read( input_file_id, degree )
    read( input_file_id, shape )
    read( input_file_id, rational )
    SLL_ALLOCATE(knots1(num_pts1+spline_deg1+1),ierr)
    SLL_ALLOCATE(knots2(num_pts2+spline_deg2+1),ierr)
    read( input_file_id, knots_1 )
    read( input_file_id, knots_2 )
    SLL_ALLOCATE(control_pts1(num_pts1*num_pts2),ierr)
    SLL_ALLOCATE(control_pts2(num_pts1*num_pts2),ierr)
    SLL_ALLOCATE(weights(num_pts1*num_pts2),ierr)
    SLL_ALLOCATE(control_pts1_2d(num_pts1,num_pts2),ierr)
    SLL_ALLOCATE(control_pts2_2d(num_pts1,num_pts2),ierr)
    SLL_ALLOCATE(weights_2d(num_pts1,num_pts2),ierr)

    read( input_file_id, control_points )
    control_pts1_2d = reshape(control_pts1,(/num_pts1,num_pts2/))
    control_pts2_2d = reshape(control_pts2,(/num_pts1,num_pts2/))
    read( input_file_id, pt_weights )
    weights_2d = reshape(weights,(/num_pts1,num_pts2/))
    read( input_file_id, logical_mesh_2d )
    close( input_file_id )


    eta1_min = knots1(1)
    eta2_min = knots2(1)
    eta1_max = knots1(num_pts1+spline_deg1+1)
    eta2_max = knots2(num_pts2+spline_deg2+1)

    ! for the moment we put the boundary condition like a dirichlet 
    ! boundary condition
    ! but we must modified this part <-- this means that this info must
    ! come within the input file: ECG

    bc_left   = SLL_DIRICHLET
    bc_right  = SLL_DIRICHLET 
    bc_bottom = SLL_DIRICHLET 
    bc_top    = SLL_DIRICHLET

!!$    ! the number of points is the knots witout the multiplicity
!!$    
    sz_knots1 = size(knots1)
!!$    SLL_ALLOCATE(nodes1(sz_knots1),ierr)
    sz_knots2 = size(knots2)
!!$    SLL_ALLOCATE(nodes2(sz_knots2),ierr)
!!$    call delete_multiplicity_in_knots(knots1,nodes1,sz_nodes1)
!!$    call delete_multiplicity_in_knots(knots2,nodes2,sz_nodes2)

    ! Initialization of the interpolator spline 2D in x
    ! ACHTUNG we have not delete it   <--- What???:ECG
    print *, 'number cells = ', number_cells1, number_cells2
    interp2d_1 => new_arbitrary_degree_spline_interp2d(&
         number_cells1 + 1,  &  
         number_cells2 + 1,  &  
         eta1_min,  &  
         eta1_max,  & 
         eta2_min,  & 
         eta2_max,  & 
         bc_left,   & 
         bc_right,  & 
         bc_bottom, & 
         bc_top,    & 
         spline_deg1, & 
         spline_deg2 )  

    call interp2d_1%set_coefficients( &
         coeffs_2d     = control_pts1_2d,&
         coeff2d_size1 = num_pts1,&
         coeff2d_size2 = num_pts2,&
         knots1        = knots1,&
         size_knots1   = sz_knots1,&
         knots2        = knots2,&
         size_knots2   = sz_knots2 )
    
!    transf%x1_interp =>interp2d_transf_1
    ! Initialization of the interpolator spline 2D in y
    ! ACHTUNG we have not delete it
    ! don't understand; write in French!! : ECG
    interp2d_2 => new_arbitrary_degree_spline_interp2d(&
         number_cells1 + 1,  & 
         number_cells2 + 2,  & 
         eta1_min,  & 
         eta1_max,  & 
         eta2_min,  & 
         eta2_max,  & 
         bc_left,   & 
         bc_right,  & 
         bc_bottom, & 
         bc_top,    & 
         spline_deg1, & 
         spline_deg2 )
    
    call interp2d_2%set_coefficients( &
         coeffs_2d     = control_pts2_2d,&
         coeff2d_size1 = num_pts1,&
         coeff2d_size2 = num_pts2,&
         knots1        = knots1,&
         size_knots1   = sz_knots1,&
         knots2        = knots2,&
         size_knots2   = sz_knots2)

    interp2d_jac => new_arbitrary_degree_spline_interp2d(&
         number_cells1 + 1,  & 
         number_cells2 + 2,  & 
         eta1_min,  & 
         eta1_max,  & 
         eta2_min,  & 
         eta2_max,  & 
         bc_left,   & 
         bc_right,  & 
         bc_bottom, & 
         bc_top,    & 
         spline_deg1, & 
         spline_deg2 )

!    transf%x2_interp =>interp2d_transf_2


    ! initialization of mesh
    mesh_2d => new_logical_mesh_2d(&
         number_cells1,&
         number_cells2,&
         eta1_min = eta1_min,&
         eta1_max = eta1_max,&
         eta2_min = eta2_min,&
         eta2_max = eta2_max)

    ! initialization of name 
   ! transf%label = trim(label)
    ! All the information from the file is now in local variables. We should
    ! now be able to initialize all the necessary objects

    ! leave the default [0,1]X[0,1] domain for the logical mesh
    !    transf%mesh => new_logical_mesh_2d(num_cells1, num_cells2)
    call transf%initialize( &
         mesh_2d, &
         label, &
         interp2d_1, &
         interp2d_2, &
         interp2d_jac )
!!$    select case (interpolator_type)
!!$       case (SLL_ARBITRARY_DEGREE_INTERPOLATOR) ! where is this specified???
!!$          transf%x1_interp => new_arbitrary_degree_spline_interpolator_2d( &
!!$)
!!$          transf%x2_interp => new_arbitrary_degree_spline_interpolator_2d( &
!!$)
!!$       case default
!!$          print *, 'ERROR, coordinate transformations 2d, ', &
!!$               'read_from_file_2d_discrete() this type of interpolator ', &
!!$               'can not be initialized from a file at present.'
!!$       end select
!    SLL_DEALLOCATE_ARRAY(nodes1,ierr)
!    SLL_DEALLOCATE_ARRAY(nodes2,ierr)
  end subroutine read_from_file_2d_discrete

  subroutine delete_multiplicity_in_knots(knots,nodes,sz_nodes)
    sll_real64, dimension(:), intent(in) :: knots
    sll_real64, dimension(:), intent(out) :: nodes
    sll_int32, intent(out) ::  sz_nodes

    sll_int32 :: i,j
    sll_int32 :: sz_knots 

    nodes(1) = knots(1)
    sz_knots = size(knots)
    j = 1
    do i = 2, sz_knots
       if (knots(i) .ne. nodes(j)) then
          j = j + 1
          nodes(j) = knots(i)
       end if
    end do
       
    sz_nodes = j
  end subroutine delete_multiplicity_in_knots

#if 0
  subroutine delete_coordinate_transformation_2D_general( transf )
    type(coordinate_transformation_2D_general), pointer :: transf
    sll_int32             :: ierr
    if( .not. associated(transf) ) then
       print *, 'ERROR, delete_coordinate_transformation_2D_general: passed map pointer ', &
            'is not associated.'
    end if
    if( associated(transf%x1_node) ) then
       SLL_DEALLOCATE( transf%x1_node, ierr )
    end if
    if( associated(transf%x2_node) ) then
       SLL_DEALLOCATE( transf%x2_node, ierr )
    end if
    if( associated(transf%x1_cell) ) then
       SLL_DEALLOCATE( transf%x1_cell, ierr )
    end if
    if( associated(transf%x2_cell) ) then
       SLL_DEALLOCATE( transf%x2_cell, ierr )
    end if
    if( associated(transf%j_matrix) ) then
       SLL_DEALLOCATE( transf%j_matrix, ierr )
    end if
    SLL_DEALLOCATE( transf%jacobians_n, ierr )
    SLL_DEALLOCATE( transf%jacobians_c, ierr )
    if( transf%transf_type .eq. DISCRETE_MAP ) then
       call delete(transf%x1_spline)
       call delete(transf%x2_spline)
    end if
    
    SLL_DEALLOCATE( transf, ierr )
  end subroutine delete_coordinate_transformation_2D_general


  ! Access functions for the mapping. These are an overkill and can be 
  ! changed by a macro, but for now, they are at least safer.
  function mesh_2d_x1_node( transf, i, j )
    sll_real64            :: mesh_2d_x1_node
    type(coordinate_transformation_2D_general), pointer :: transf
    sll_int32, intent(in) :: i
    sll_int32, intent(in) :: j
    sll_int32 :: num_pts_1
    sll_int32 :: num_pts_2
    SLL_ASSERT( associated(transf) )
    num_pts_1 = transf%num_pts_1
    num_pts_2 = transf%num_pts_2
    SLL_ASSERT( (i .ge. 1) .and. (i .le. num_pts_1) )
    SLL_ASSERT( (j .ge. 1) .and. (j .le. num_pts_2) )
    mesh_2d_x1_node = transf%x1_node(i,j)
  end function mesh_2d_x1_node

  function mesh_2d_x2_node( transf, i, j )
    sll_real64            :: mesh_2d_x2_node
    type(coordinate_transformation_2D_general), pointer :: transf
    sll_int32, intent(in) :: i
    sll_int32, intent(in) :: j
    sll_int32 :: num_pts_1
    sll_int32 :: num_pts_2
    SLL_ASSERT( associated(transf) )
    num_pts_1 = transf%num_pts_1
    num_pts_2 = transf%num_pts_2
    SLL_ASSERT( (i .ge. 1) .and. (i .le. num_pts_1) )
    SLL_ASSERT( (j .ge. 1) .and. (j .le. num_pts_2) )
    mesh_2d_x2_node = transf%x2_node(i,j)
  end function mesh_2d_x2_node

  function mesh_2d_x1_cell( transf, i, j )
    sll_real64            :: mesh_2d_x1_cell
    type(coordinate_transformation_2D_general), pointer :: transf
    sll_int32, intent(in) :: i
    sll_int32, intent(in) :: j
    sll_int32 :: num_pts_1
    sll_int32 :: num_pts_2
    SLL_ASSERT( associated(transf) )
    num_pts_1 = transf%num_pts_1
    num_pts_2 = transf%num_pts_2
    SLL_ASSERT( (i .ge. 1) .and. (i .le. num_pts_1 - 1) )
    SLL_ASSERT( (j .ge. 1) .and. (j .le. num_pts_2 - 1) )
    mesh_2d_x1_cell = transf%x1_cell(i,j)
  end function mesh_2d_x1_cell

  function mesh_2d_x2_cell( transf, i, j )
    sll_real64            :: mesh_2d_x2_cell
    type(coordinate_transformation_2D_general), pointer :: transf
    sll_int32, intent(in) :: i
    sll_int32, intent(in) :: j
    sll_int32 :: num_pts_1
    sll_int32 :: num_pts_2
    SLL_ASSERT( associated(transf) )
    num_pts_1 = transf%num_pts_1
    num_pts_2 = transf%num_pts_2
    SLL_ASSERT( (i .ge. 1) .and. (i .le. num_pts_1) )
    SLL_ASSERT( (j .ge. 1) .and. (j .le. num_pts_2) )
    mesh_2d_x2_cell = transf%x2_cell(i,j)
  end function mesh_2d_x2_cell

  function mesh_2d_jacobian_cell( transf, i, j )
    sll_real64              :: mesh_2d_jacobian_cell
    type(coordinate_transformation_2D_general), pointer   :: transf
    sll_int32, intent(in)   :: i
    sll_int32, intent(in)   :: j
    sll_int32 :: num_cells_1
    sll_int32 :: num_cells_2
    SLL_ASSERT( associated(transf) )
    num_cells_1 = transf%num_pts_1 - 1
    num_cells_2 = transf%num_pts_2 - 1
    SLL_ASSERT( (i .ge. 1) .and. (i .le. num_cells_1) )
    SLL_ASSERT( (j .ge. 1) .and. (j .le. num_cells_2) )
    mesh_2d_jacobian_cell = transf%jacobians_c(i,j)
  end function mesh_2d_jacobian_cell
#endif
end module sll_module_coordinate_transformations_2d
