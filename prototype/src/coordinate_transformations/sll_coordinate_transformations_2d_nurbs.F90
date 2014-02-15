module sll_module_coordinate_transformations_2d_nurbs
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_file_io.h"
  use sll_xdmf
  use sll_logical_meshes
  use sll_cubic_spline_interpolator_2d
  use sll_gnuplot
  use sll_module_interpolators_2d_base
  use sll_coordinate_transformation_2d_base_module
  use sll_module_deboor_splines_2d

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

  ! -----------------------------------------------------------------------
  !
  !                Discrete case : Nurbs/Splines specifications
  !
  ! -----------------------------------------------------------------------

  type, extends(sll_coordinate_transformation_2d_base) :: &
       sll_coordinate_transformation_2d_nurbs

     sll_real64, dimension(:,:), pointer :: x1_node =>null()  ! x1(i,j) 
     sll_real64, dimension(:,:), pointer :: x2_node =>null()  ! x2(i,j) 
     sll_real64, dimension(:,:), pointer :: x1_cell =>null()
     sll_real64, dimension(:,:), pointer :: x2_cell =>null()
     class(sll_interpolator_2d_base), pointer :: x1_interp =>null()
     class(sll_interpolator_2d_base), pointer :: x2_interp =>null()
     class(sll_interpolator_2d_base), pointer :: x3_interp =>null()
     sll_int32 :: is_rational
     type(sll_logical_mesh_2d), pointer  :: mesh2d_minimal =>null()
!     type(sll_logical_mesh_2d), pointer :: mesh
   contains
     procedure, pass(transf) :: get_logical_mesh => get_logical_mesh_nurbs_2d
     procedure, pass(transf) :: x1_at_node => x1_node_nurbs
     procedure, pass(transf) :: x2_at_node => x2_node_nurbs
     procedure, pass(transf) :: jacobian_at_node =>transf_2d_jacobian_node_nurbs
     procedure, pass(transf) :: x1         => x1_nurbs
     procedure, pass(transf) :: x2         => x2_nurbs
     procedure, pass(transf) :: x1_at_cell => x1_cell_nurbs
     procedure, pass(transf) :: x2_at_cell => x2_cell_nurbs
     procedure, pass(transf) :: jacobian_at_cell => jacobian_2d_cell_nurbs
     procedure, pass(transf) :: jacobian   => jacobian_2d_nurbs
     procedure, pass(transf) :: jacobian_matrix => jacobian_matrix_2d_nurbs
     procedure, pass(transf) :: inverse_jacobian_matrix => &
          inverse_jacobian_matrix_2d_nurbs
     procedure, pass(transf) :: write_to_file => write_to_file_2d_nurbs
     procedure, pass(transf) :: read_from_file => read_from_file_2d_nurbs
     procedure, pass(transf) :: delete => delete_transformation_2d_nurbs
  end type sll_coordinate_transformation_2d_nurbs

  type sll_coordinate_transformation_2d_nurbs_ptr
     type(sll_coordinate_transformation_2d_nurbs), pointer :: T
  end type sll_coordinate_transformation_2d_nurbs_ptr

  interface delete
     module procedure delete_transformation_2d_nurbs
  end interface

  
contains
  !**************************************************************************
  !
  !        Functions for the general transformation defined with nurbs
  !
  !**************************************************************************


  ! -------------------------------------------------------------------------
  ! The nurbs-based coordinate transformation is associated with 2 logical 
  ! meshes: the 'minimal' logical mesh, which is implicit in the spline 
  ! description and  that is able to faitfully represent a given geometry. 
  ! The other logical mesh, which is the usual logical mesh on which we want 
  ! to represent quantities. 
  ! This needs to be present to answer questions like T%x_node(i,j). 
  ! We set this logical mesh outside of the read_from_file routine.
  ! -------------------------------------------------------------------------
  function new_nurbs_2d_transformation_from_file( filename ) result(res)
    type(sll_coordinate_transformation_2d_nurbs), pointer :: res
    character(len=*), intent(in) :: filename
    sll_int32 :: ierr
    SLL_ALLOCATE(res,ierr)
    call read_from_file_2d_nurbs( res, filename )
  end function new_nurbs_2d_transformation_from_file

  subroutine read_from_file_2d_nurbs( transf, filename )
    use sll_arbitrary_degree_spline_interpolator_2d_module
    class(sll_coordinate_transformation_2d_nurbs), intent(inout) :: transf
    character(len=*), intent(in) :: filename
    intrinsic :: trim
    !sll_int32 :: interpolator_type
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
    sll_real64 :: eta1_min_minimal
    sll_real64 :: eta1_max_minimal
    sll_real64 :: eta2_min_minimal
    sll_real64 :: eta2_max_minimal
    sll_int32  :: bc_left
    sll_int32  :: bc_right
    sll_int32  :: bc_bottom
    sll_int32  :: bc_top
    sll_int32  :: number_cells1,number_cells2
    sll_int32 :: sz_knots1,sz_knots2
    sll_int32 :: i,j
  
    namelist /transf_label/  label
    namelist /degree/   spline_deg1, spline_deg2
    namelist /shape/    num_pts1, num_pts2 ! it is not the number of points but the number of coeff sdpline in each direction !!
    namelist /rational/ is_rational
    namelist /knots_1/   knots1
    namelist /knots_2/   knots2
    namelist /control_points/ control_pts1, control_pts2
    namelist /pt_weights/  weights
    namelist /logical_mesh_2d/ number_cells1,number_cells2
    !character(len=80) :: line_buffer

    if(len(filename) >= 256) then
       print *, 'ERROR, read_coefficients_from_file => ',&
            'read_from_file_nurbs():',&
            'filenames longer than 256 characters are not allowed.'
       STOP
    end if
    filename_local = trim(filename)

    ! get a new identifier for the file.
    call sll_new_file_id( input_file_id, ierr )
    if( ierr .ne. 0 ) then
       print *, 'ERROR while trying to obtain an unique identifier for file ',&
            filename, '. Called from read_from_file_2d_nurbs().'
       stop
    end if
    open(unit=input_file_id, file=filename_local, STATUS="OLD", IOStat=IO_stat)
    if( IO_Stat .ne. 0 ) then
       print *, 'ERROR while opening file ',filename, &
            '. Called from read_from_file_2d_nurbs().'
       stop
    end if

    ! read the label
    read( input_file_id, transf_label )
    ! read the degree of spline
    read( input_file_id, degree )
    ! read ....?
    read( input_file_id, shape )
    ! read if we use NURBS or not
    ! i.e. if we have NURBS rational == 1
    ! otherwise rational == 0
    read( input_file_id, rational )
    transf%is_rational = is_rational
    ! Allocations of knots to construct the splines
    SLL_ALLOCATE(knots1(num_pts1+spline_deg1+1),ierr)
    SLL_ALLOCATE(knots2(num_pts2+spline_deg2+1),ierr)
    ! read the knots associated to each direction 
    read( input_file_id, knots_1 )
    read( input_file_id, knots_2 )
    
    ! allocations of tables containing control points in each direction 
    ! here its table 1D
    SLL_ALLOCATE(control_pts1(num_pts1*num_pts2),ierr)
    SLL_ALLOCATE(control_pts2(num_pts1*num_pts2),ierr)
    ! allocation of table containing the weights associated to each control 
    ! points. It is used only if the transformation is a NURBS 
    ! i.e. if we have rational == 1
    ! here its table 1D
    SLL_ALLOCATE(weights(num_pts1*num_pts2),ierr)
    
    ! allocations of tables containing control points in each direction 
    ! here its table 2D
    SLL_ALLOCATE(control_pts1_2d(num_pts1,num_pts2),ierr)
    SLL_ALLOCATE(control_pts2_2d(num_pts1,num_pts2),ierr)
    ! allocation of table containing the weights associated to each control 
    ! points. It is used only if the transformation is a NURBS 
    ! i.e. if we have rational == 1
    ! here its table 2D
    SLL_ALLOCATE(weights_2d(num_pts1,num_pts2),ierr)

    ! read the control points in the file
    read( input_file_id, control_points )
    ! reshape the control points to use them in the interpolator
    control_pts1_2d = reshape(control_pts1,(/num_pts1,num_pts2/))
    control_pts2_2d = reshape(control_pts2,(/num_pts1,num_pts2/))
    ! read the weight in the file associated in each control points
    read( input_file_id, pt_weights )
    ! reshape the control points to use them in the rational interpolator
    ! only if the transformation is a NURBS
    weights_2d = reshape(weights,(/num_pts1,num_pts2/))

    ! read the minimal mesh corresponding to the transformation
    read( input_file_id, logical_mesh_2d )
    ! close the file to begin the work 
    close( input_file_id )


    !! if the transformation is a NURBS i.e is_rational == 1
    !! we have  
    !!         sum_(i,j) P^x_(i,j) weight(i,j) B_i^alpha1(eta1)B_j^alpha2(eta2) 
    !!  X =   ----------------------------------------------------------------
    !!              sum_(i,j) weight(i,j) B_i^alpha1(eta1)B_j^alpha2(eta2) 
    !!
    !! we have      
    !!
    !!         sum_(i,j) P^y_(i,j) weight(i,j) B_i^alpha1(eta1)B_j^alpha2(eta2) 
    !!    Y =  ----------------------------------------------------------------
    !!              sum_(i,j) weight(i,j) B_i^alpha1(eta1)B_j^alpha2(eta2) 
    !!
    !! that's why we multiply each control points by their weights
    !! to use it in the coefficients splines directly in the first and the 
    !! second interpolator 
    !! the third interpolator is use to define the term under the fraction
    !!
    !! if the transformation is a SPLINE i.e. is_rational == 0
    !! we have
    !!      
    !!    X =       sum_(i,j) P^x_(i,j)  B_i^alpha1(eta1)B_j^alpha2(eta2) 
    !!
    !! we have
    !!      
    !!    Y =       sum_(i,j) P^y_(i,j) B_i^alpha1(eta1)B_j^alpha2(eta2)
    !!
    !! that's why we don't multiply each control points by their weights
    !! to use it in the coefficients splines directly in the first and the 
    !! second interpolator 
    

    if (transf%is_rational ==1) then
       
       do i = 1, num_pts1 
          do j = 1, num_pts2
             control_pts1_2d(i,j) = control_pts1_2d(i,j)*weights_2d(i,j)
             control_pts2_2d(i,j) = control_pts2_2d(i,j)*weights_2d(i,j)
          end do
       end do
    end if

    ! Is this worth it? Is there the expectation that the extreme values
    ! coming from CAID are anything other than 0 and 1??
    eta1_min_minimal = knots1(1)
    eta2_min_minimal = knots2(1)
    ! Aurore: discuss, we could simply put the knots array without the 
    ! duplicates in the .nml file
    eta1_max_minimal = knots1(num_pts1+spline_deg1+1) 
    eta2_max_minimal = knots2(num_pts2+spline_deg2+1)

    ! for the moment we put the boundary condition like a dirichlet 
    ! boundary condition
    ! but we must modified this part <-- this means that this info must
    ! come within the input file: ECG

    bc_left   = SLL_DIRICHLET
    bc_right  = SLL_DIRICHLET 
    bc_bottom = SLL_DIRICHLET 
    bc_top    = SLL_DIRICHLET

!!$    ! the number of points is the knots witout the multiplicity
   
    sz_knots1 = size(knots1)
    sz_knots2 = size(knots2)

    ! Initialize the first interpolator for 
    ! the first component of our change of coordinates

    transf%x1_interp => new_arbitrary_degree_spline_interp2d(&
         number_cells1 + 1,  &  
         number_cells2 + 1,  &  
         eta1_min_minimal,  &  
         eta1_max_minimal,  & 
         eta2_min_minimal,  & 
         eta2_max_minimal,  & 
         bc_left,   & 
         bc_right,  & 
         bc_bottom, & 
         bc_top,    & 
         spline_deg1, & 
         spline_deg2 )  

    ! stock all the control points for the first interpolator 
    ! to compute the first component of our change of coordinates
    call transf%x1_interp%set_coefficients( &
         coeffs_2d     = control_pts1_2d,&
         coeff2d_size1 = num_pts1,&
         coeff2d_size2 = num_pts2,&
         knots1        = knots1,&
         size_knots1   = sz_knots1,&
         knots2        = knots2,&
         size_knots2   = sz_knots2 )
    

    ! Initialize the second interpolator for 
    ! the second component of our change of coordinates

    transf%x2_interp => new_arbitrary_degree_spline_interp2d(&
         number_cells1 + 1,  & 
         number_cells2 + 2,  & 
         eta1_min_minimal,  & 
         eta1_max_minimal,  & 
         eta2_min_minimal,  & 
         eta2_max_minimal,  & 
         bc_left,   & 
         bc_right,  & 
         bc_bottom, & 
         bc_top,    & 
         spline_deg1, & 
         spline_deg2 )

    
    ! stock all the control points for the second interpolator 
    ! to compute the second component of our change of coordinates

    call transf%x2_interp%set_coefficients( &
         coeffs_2d     = control_pts2_2d,&
         coeff2d_size1 = num_pts1,&
         coeff2d_size2 = num_pts2,&
         knots1        = knots1,&
         size_knots1   = sz_knots1,&
         knots2        = knots2,&
         size_knots2   = sz_knots2)
    
    ! Initialize the third interpolator for 
    ! the rational component of our change of coordinates
    ! That is useful if the transfromation is NURBS

    transf%x3_interp => new_arbitrary_degree_spline_interp2d(&
         number_cells1 + 1,  & 
         number_cells2 + 1,  & 
         eta1_min_minimal,  & 
         eta1_max_minimal,  & 
         eta2_min_minimal,  & 
         eta2_max_minimal,  & 
         bc_left,   & 
         bc_right,  & 
         bc_bottom, & 
         bc_top,    & 
         spline_deg1, & 
         spline_deg2 )
    
    
    ! stock all the weight for the rationnal interpolator 
    ! to compute the rationnal component of our change of coordinates
    ! in the case of NURBS
    
    call transf%x3_interp%set_coefficients( &
         coeffs_2d     = weights_2d,&
         coeff2d_size1 = num_pts1,&
         coeff2d_size2 = num_pts2,&
         knots1        = knots1,&
         size_knots1   = sz_knots1,&
         knots2        = knots2,&
         size_knots2   = sz_knots2)


    ! initialization of minimal mesh given by file.                             
    ! This is temporary, we have not decided if this object should take
    ! possession of the logical mesh or not... For now we keep the minimum
    ! information related with the number of cells to at least be able to
    ! initialize a logical mesh outside of the object.

    transf%mesh2d_minimal => new_logical_mesh_2d(&
         number_cells1,&
         number_cells2,&
         eta1_min = eta1_min_minimal,&
         eta1_max = eta1_max_minimal,&
         eta2_min = eta2_min_minimal,&
         eta2_max = eta2_max_minimal)

    transf%mesh  => null()
    transf%label =  trim(label)
  end subroutine read_from_file_2d_nurbs

  function get_logical_mesh_nurbs_2d( transf ) result(res)
    type(sll_logical_mesh_2d), pointer :: res
    class(sll_coordinate_transformation_2d_nurbs), intent(in) :: transf
    res => transf%mesh2d_minimal
  end function get_logical_mesh_nurbs_2d

  function x1_node_nurbs( transf, i, j ) result(val)
    class(sll_coordinate_transformation_2d_nurbs) :: transf
    type(sll_logical_mesh_2d), pointer :: lm
    sll_real64             :: val
    sll_int32, intent(in) :: i
    sll_int32, intent(in) :: j
    sll_real64  :: eta1
    sll_real64  :: eta2
    sll_real64  :: delta1
    sll_real64  :: delta2
    sll_real64  :: eta1_min
    sll_real64  :: eta2_min

    lm => transf%get_logical_mesh()

    eta1_min = lm%eta1_min
    eta2_min = lm%eta2_min
    delta1   = lm%delta_eta1
    delta2   = lm%delta_eta2
    
    eta1 = eta1_min + (i-1) * delta1 
    eta2 = eta2_min + (j-1) * delta2 

    SLL_ASSERT( eta1 <= 1.0_f64)
    SLL_ASSERT( eta1 >= 0.0_f64)
    SLL_ASSERT( eta2 <= 1.0_f64)
    SLL_ASSERT( eta2 >= 0.0_f64)

    if (transf%is_rational == 0) then ! IN the case of SPLINE
       val = transf%x1_interp%interpolate_value(eta1,eta2)
    else ! In the case of NURBS
       val = transf%x1_interp%interpolate_value(eta1,eta2)/&
            transf%x3_interp%interpolate_value(eta1,eta2)
    end if
       

  end function x1_node_nurbs

   function x2_node_nurbs( transf, i, j ) result(val)
    class(sll_coordinate_transformation_2d_nurbs) :: transf
    type(sll_logical_mesh_2d), pointer :: lm
    sll_real64             :: val
    sll_int32, intent(in) :: i
    sll_int32, intent(in) :: j
    sll_real64  :: eta1
    sll_real64  :: eta2
    sll_real64  :: delta1
    sll_real64  :: delta2
    sll_real64  :: eta1_min
    sll_real64  :: eta2_min

    lm => transf%get_logical_mesh()

    eta1_min = lm%eta1_min
    eta2_min = lm%eta2_min
    delta1   = lm%delta_eta1
    delta2   = lm%delta_eta2
    
    eta1 = eta1_min + (i-1) * delta1 
    eta2 = eta2_min + (j-1) * delta2 

    SLL_ASSERT( eta1 <= 1.0_f64)
    SLL_ASSERT( eta1 >= 0.0_f64)
    SLL_ASSERT( eta2 <= 1.0_f64)
    SLL_ASSERT( eta2 >= 0.0_f64)

    if (transf%is_rational == 0) then ! IN the case of SPLINE
       val = transf%x2_interp%interpolate_value(eta1,eta2)
    else ! In the case of NURBS
       val = transf%x2_interp%interpolate_value(eta1,eta2)/&
            transf%x3_interp%interpolate_value(eta1,eta2)
    end if

  end function x2_node_nurbs

  function x1_cell_nurbs( transf, i, j ) result(val)
    class(sll_coordinate_transformation_2d_nurbs) :: transf
    type(sll_logical_mesh_2d), pointer :: lm
    sll_real64             :: val
    sll_int32, intent(in) :: i
    sll_int32, intent(in) :: j
    sll_real64  :: eta1
    sll_real64  :: eta2
    sll_real64  :: delta1
    sll_real64  :: delta2
    sll_real64  :: eta1_min
    sll_real64  :: eta2_min
    
    lm => transf%get_logical_mesh()

    eta1_min = lm%eta1_min
    eta2_min = lm%eta2_min
    delta1   = lm%delta_eta1
    delta2   = lm%delta_eta2

    SLL_ASSERT( i <= lm%num_cells1)
    SLL_ASSERT( j <= lm%num_cells2)
    
    eta1 = eta1_min + (i-0.5_f64) * delta1  
    eta2 = eta2_min + (j-0.5_f64) * delta2 
    
    SLL_ASSERT( eta1 <= 1.0_f64)
    SLL_ASSERT( eta1 >= 0.0_f64)
    SLL_ASSERT( eta2 <= 1.0_f64)
    SLL_ASSERT( eta2 >= 0.0_f64)

    if (transf%is_rational == 0) then ! IN the case of SPLINE
       val = transf%x1_interp%interpolate_value(eta1,eta2)
    else ! In the case of NURBS
       val = transf%x1_interp%interpolate_value(eta1,eta2)/&
            transf%x3_interp%interpolate_value(eta1,eta2)
    end if
    
  end function x1_cell_nurbs
  
   function x2_cell_nurbs( transf, i, j ) result(val)
    class(sll_coordinate_transformation_2d_nurbs) :: transf
    type(sll_logical_mesh_2d), pointer :: lm
    sll_real64             :: val
    sll_int32, intent(in) :: i
    sll_int32, intent(in) :: j
    sll_real64  :: eta1
    sll_real64  :: eta2
    sll_real64  :: delta1
    sll_real64  :: delta2
    sll_real64  :: eta1_min
    sll_real64  :: eta2_min

    lm => transf%get_logical_mesh()
    eta1_min = lm%eta1_min
    eta2_min = lm%eta2_min
    delta1   = lm%delta_eta1
    delta2   = lm%delta_eta2
    
    SLL_ASSERT( i <= lm%num_cells1)
    SLL_ASSERT( j <= lm%num_cells2)
    
    eta1 = eta1_min + (i-0.5_f64) * delta1 
    eta2 = eta2_min + (j-0.5_f64) * delta2 
    
    SLL_ASSERT( eta1 <= 1.0_f64)
    SLL_ASSERT( eta1 >= 0.0_f64)
    SLL_ASSERT( eta2 <= 1.0_f64)
    SLL_ASSERT( eta2 >= 0.0_f64)
    
    if (transf%is_rational == 0) then ! IN the case of SPLINE
       val = transf%x2_interp%interpolate_value(eta1,eta2)
    else ! In the case of NURBS
       val = transf%x2_interp%interpolate_value(eta1,eta2)/&
            transf%x3_interp%interpolate_value(eta1,eta2)
    end if
    
  end function x2_cell_nurbs
  
  function x1_nurbs( transf, eta1, eta2 ) result(val)
    class(sll_coordinate_transformation_2d_nurbs) :: transf
    sll_real64             :: val
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    
    SLL_ASSERT( eta1 <= 1.0_f64)
    SLL_ASSERT( eta1 >= 0.0_f64)
    SLL_ASSERT( eta2 <= 1.0_f64)
    SLL_ASSERT( eta2 >= 0.0_f64)
    
    if (transf%is_rational == 0) then ! IN the case of SPLINE
       val = transf%x1_interp%interpolate_value(eta1,eta2)
    else ! In the case of NURBS
       val = transf%x1_interp%interpolate_value(eta1,eta2)/&
            transf%x3_interp%interpolate_value(eta1,eta2)
    end if
  end function x1_nurbs

  
  function x2_nurbs( transf, eta1, eta2 ) result(val)
    class(sll_coordinate_transformation_2d_nurbs) :: transf
    sll_real64             :: val
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
        
    SLL_ASSERT( eta1 <= 1.0_f64)
    SLL_ASSERT( eta1 >= 0.0_f64)
    SLL_ASSERT( eta2 <= 1.0_f64)
    SLL_ASSERT( eta2 >= 0.0_f64)
    

    if (transf%is_rational == 0) then ! IN the case of SPLINE
       val = transf%x2_interp%interpolate_value(eta1,eta2)
    else ! In the case of NURBS
       val = transf%x2_interp%interpolate_value(eta1,eta2)/&
            transf%x3_interp%interpolate_value(eta1,eta2)
    end if
  end function x2_nurbs

  
  function jacobian_2d_nurbs( transf, eta1, eta2 ) result(jac)
    class(sll_coordinate_transformation_2d_nurbs) :: transf
    sll_real64             :: jac
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    sll_real64             :: j11
    sll_real64             :: j12
    sll_real64             :: j21
    sll_real64             :: j22

    SLL_ASSERT( eta1 <= 1.0_f64)
    SLL_ASSERT( eta1 >= 0.0_f64)
    SLL_ASSERT( eta2 <= 1.0_f64)
    SLL_ASSERT( eta2 >= 0.0_f64)
    
    if (transf%is_rational == 0) then ! IN the case of SPLINE
       
       j11 = transf%x1_interp%interpolate_derivative_eta1( eta1, eta2 )
       j12 = transf%x1_interp%interpolate_derivative_eta2( eta1, eta2 )
       j21 = transf%x2_interp%interpolate_derivative_eta1( eta1, eta2 )
       j22 = transf%x2_interp%interpolate_derivative_eta2( eta1, eta2 )
       jac = j11*j22 - j12*j21
       
    else 

       j11 = transf%x1_interp%interpolate_derivative_eta1( eta1, eta2 )&
            *transf%x3_interp%interpolate_value(eta1,eta2) &
            - transf%x3_interp%interpolate_derivative_eta1( eta1, eta2 )&
            *transf%x1_interp%interpolate_value(eta1,eta2)

       j12 = transf%x1_interp%interpolate_derivative_eta2( eta1, eta2 )&
            *transf%x3_interp%interpolate_value(eta1,eta2) &
            - transf%x3_interp%interpolate_derivative_eta2( eta1, eta2 )&
            *transf%x1_interp%interpolate_value(eta1,eta2)
       
       j21 = transf%x2_interp%interpolate_derivative_eta1( eta1, eta2 )&
            *transf%x3_interp%interpolate_value(eta1,eta2) &
            - transf%x3_interp%interpolate_derivative_eta1( eta1, eta2 )&
            *transf%x2_interp%interpolate_value(eta1,eta2)

       j22 = transf%x2_interp%interpolate_derivative_eta2( eta1, eta2 )&
            *transf%x3_interp%interpolate_value(eta1,eta2) &
            - transf%x3_interp%interpolate_derivative_eta2( eta1, eta2 )&
            *transf%x2_interp%interpolate_value(eta1,eta2)
       
       
       jac = (j11*j22 - j12*j21)/&
            (transf%x3_interp%interpolate_value(eta1,eta2))**4
    end if
    
  end function jacobian_2d_nurbs
  
  function transf_2d_jacobian_node_nurbs( transf, i, j )
    class(sll_coordinate_transformation_2d_nurbs)   :: transf
    type(sll_logical_mesh_2d), pointer :: lm
    sll_real64              :: transf_2d_jacobian_node_nurbs
    sll_int32, intent(in)   :: i
    sll_int32, intent(in)   :: j
    sll_real64  :: eta1
    sll_real64  :: eta2
    sll_real64  :: delta1
    sll_real64  :: delta2
    sll_real64  :: eta1_min
    sll_real64  :: eta2_min
    sll_real64  :: j11
    sll_real64  :: j12
    sll_real64  :: j21
    sll_real64  :: j22
    
    lm => transf%get_logical_mesh()

    eta1_min = lm%eta1_min
    eta2_min = lm%eta2_min
    delta1   = lm%delta_eta1
    delta2   = lm%delta_eta2
    
    eta1 = eta1_min + (i-1) * delta1 
    eta2 = eta2_min + (j-1) * delta2 
    
    SLL_ASSERT( eta1 <= 1.0_f64)
    SLL_ASSERT( eta1 >= 0.0_f64)
    SLL_ASSERT( eta2 <= 1.0_f64)
    SLL_ASSERT( eta2 >= 0.0_f64)

    if (transf%is_rational == 0) then ! IN the case of SPLINE
       
       j11 = transf%x1_interp%interpolate_derivative_eta1( eta1, eta2 )
       j12 = transf%x1_interp%interpolate_derivative_eta2( eta1, eta2 )
       j21 = transf%x2_interp%interpolate_derivative_eta1( eta1, eta2 )
       j22 = transf%x2_interp%interpolate_derivative_eta2( eta1, eta2 )
       transf_2d_jacobian_node_nurbs = j11*j22 - j12*j21
       
    else 
       
       j11 = transf%x1_interp%interpolate_derivative_eta1( eta1, eta2 )&
            *transf%x3_interp%interpolate_value(eta1,eta2) &
            - transf%x3_interp%interpolate_derivative_eta1( eta1, eta2 )&
            *transf%x1_interp%interpolate_value(eta1,eta2)
       
       j12 = transf%x1_interp%interpolate_derivative_eta2( eta1, eta2 )&
            *transf%x3_interp%interpolate_value(eta1,eta2) &
            - transf%x3_interp%interpolate_derivative_eta2( eta1, eta2 )&
            *transf%x1_interp%interpolate_value(eta1,eta2)
       
       j21 = transf%x2_interp%interpolate_derivative_eta1( eta1, eta2 )&
            *transf%x3_interp%interpolate_value(eta1,eta2) &
            - transf%x3_interp%interpolate_derivative_eta1( eta1, eta2 )&
            *transf%x2_interp%interpolate_value(eta1,eta2)
       
       j22 = transf%x2_interp%interpolate_derivative_eta2( eta1, eta2 )&
            *transf%x3_interp%interpolate_value(eta1,eta2) &
            - transf%x3_interp%interpolate_derivative_eta2( eta1, eta2 )&
            *transf%x2_interp%interpolate_value(eta1,eta2)
       
       transf_2d_jacobian_node_nurbs = (j11*j22 - j12*j21)/&
            (transf%x3_interp%interpolate_value(eta1,eta2))**2
    end if
    
  end function transf_2d_jacobian_node_nurbs

  
  function jacobian_2d_cell_nurbs( transf, i, j ) result(var)
    class(sll_coordinate_transformation_2d_nurbs) :: transf
    type(sll_logical_mesh_2d), pointer :: lm
    sll_real64                         :: var
    sll_int32, intent(in)              :: i
    sll_int32, intent(in)              :: j
    sll_real64  :: eta1
    sll_real64  :: eta2
    sll_real64  :: delta1
    sll_real64  :: delta2
    sll_real64  :: eta1_min
    sll_real64  :: eta2_min
    sll_real64  :: j11
    sll_real64  :: j12
    sll_real64  :: j21
    sll_real64  :: j22

    lm => transf%get_logical_mesh()

    eta1_min = lm%eta1_min
    eta2_min = lm%eta2_min
    delta1   = lm%delta_eta1
    delta2   = lm%delta_eta2
    
    SLL_ASSERT( i <= lm%num_cells1)
    SLL_ASSERT( j <= lm%num_cells2)
    
    eta1 = eta1_min + (i-0.5_f64) * delta1 
    eta2 = eta2_min + (j-0.5_f64) * delta2 
    
    SLL_ASSERT( eta1 <= 1.0_f64)
    SLL_ASSERT( eta1 >= 0.0_f64)
    SLL_ASSERT( eta2 <= 1.0_f64)
    SLL_ASSERT( eta2 >= 0.0_f64)
    
    if (transf%is_rational == 0) then ! IN the case of SPLINE
       
       j11 = transf%x1_interp%interpolate_derivative_eta1( eta1, eta2 )
       j12 = transf%x1_interp%interpolate_derivative_eta2( eta1, eta2 )
       j21 = transf%x2_interp%interpolate_derivative_eta1( eta1, eta2 )
       j22 = transf%x2_interp%interpolate_derivative_eta2( eta1, eta2 )
       var = j11*j22 - j12*j21
       
    else 
       
       j11 = transf%x1_interp%interpolate_derivative_eta1( eta1, eta2 )&
            *transf%x3_interp%interpolate_value(eta1,eta2) &
            - transf%x3_interp%interpolate_derivative_eta1( eta1, eta2 )&
            *transf%x1_interp%interpolate_value(eta1,eta2)
       
       j12 = transf%x1_interp%interpolate_derivative_eta2( eta1, eta2 )&
            *transf%x3_interp%interpolate_value(eta1,eta2) &
            - transf%x3_interp%interpolate_derivative_eta2( eta1, eta2 )&
            *transf%x1_interp%interpolate_value(eta1,eta2)
       
       j21 = transf%x2_interp%interpolate_derivative_eta1( eta1, eta2 )&
            *transf%x3_interp%interpolate_value(eta1,eta2) &
            - transf%x3_interp%interpolate_derivative_eta1( eta1, eta2 )&
            *transf%x2_interp%interpolate_value(eta1,eta2)
       
       j22 = transf%x2_interp%interpolate_derivative_eta2( eta1, eta2 )&
            *transf%x3_interp%interpolate_value(eta1,eta2) &
            - transf%x3_interp%interpolate_derivative_eta2( eta1, eta2 )&
            *transf%x2_interp%interpolate_value(eta1,eta2)
       
       
       var = (j11*j22 - j12*j21)/&
            (transf%x3_interp%interpolate_value(eta1,eta2))**2
    end if
    
  end function jacobian_2d_cell_nurbs

  
  function jacobian_matrix_2d_nurbs( transf, eta1, eta2 )
    class(sll_coordinate_transformation_2d_nurbs) :: transf
    sll_real64, dimension(1:2,1:2)     :: jacobian_matrix_2d_nurbs
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    sll_real64             :: j11
    sll_real64             :: j12
    sll_real64             :: j21
    sll_real64             :: j22

    SLL_ASSERT( eta1 <= 1.0_f64)
    SLL_ASSERT( eta1 >= 0.0_f64)
    SLL_ASSERT( eta2 <= 1.0_f64)
    SLL_ASSERT( eta2 >= 0.0_f64)
    
    if (transf%is_rational == 0) then ! IN the case of SPLINE
       
       j11 = transf%x1_interp%interpolate_derivative_eta1( eta1, eta2 )
       j12 = transf%x1_interp%interpolate_derivative_eta2( eta1, eta2 )
       j21 = transf%x2_interp%interpolate_derivative_eta1( eta1, eta2 )
       j22 = transf%x2_interp%interpolate_derivative_eta2( eta1, eta2 )

       jacobian_matrix_2d_nurbs(1,1) = j11
       jacobian_matrix_2d_nurbs(1,2) = j12
       jacobian_matrix_2d_nurbs(2,1) = j21
       jacobian_matrix_2d_nurbs(2,2) = j22
       
    else 

       j11 = transf%x1_interp%interpolate_derivative_eta1( eta1, eta2 )&
            *transf%x3_interp%interpolate_value(eta1,eta2) &
            - transf%x3_interp%interpolate_derivative_eta1( eta1, eta2 )&
            *transf%x1_interp%interpolate_value(eta1,eta2)

       j12 = transf%x1_interp%interpolate_derivative_eta2( eta1, eta2 )&
            *transf%x3_interp%interpolate_value(eta1,eta2) &
            - transf%x3_interp%interpolate_derivative_eta2( eta1, eta2 )&
            *transf%x1_interp%interpolate_value(eta1,eta2)
       
       j21 = transf%x2_interp%interpolate_derivative_eta1( eta1, eta2 )&
            *transf%x3_interp%interpolate_value(eta1,eta2) &
            - transf%x3_interp%interpolate_derivative_eta1( eta1, eta2 )&
            *transf%x2_interp%interpolate_value(eta1,eta2)

       j22 = transf%x2_interp%interpolate_derivative_eta2( eta1, eta2 )&
            *transf%x3_interp%interpolate_value(eta1,eta2) &
            - transf%x3_interp%interpolate_derivative_eta2( eta1, eta2 )&
            *transf%x2_interp%interpolate_value(eta1,eta2)

       jacobian_matrix_2d_nurbs(1,1) = j11/&
            (transf%x3_interp%interpolate_value(eta1,eta2))**2
       jacobian_matrix_2d_nurbs(1,2) = j12/&
            (transf%x3_interp%interpolate_value(eta1,eta2))**2
       jacobian_matrix_2d_nurbs(2,1) = j21/&
            (transf%x3_interp%interpolate_value(eta1,eta2))**2
       jacobian_matrix_2d_nurbs(2,2) = j22/&
            (transf%x3_interp%interpolate_value(eta1,eta2))**2 

      
    end if
    
  end function jacobian_matrix_2d_nurbs

  
  function inverse_jacobian_matrix_2d_nurbs( transf, eta1, eta2 )
    class(sll_coordinate_transformation_2d_nurbs) :: transf
    sll_real64, dimension(1:2,1:2)     :: inverse_jacobian_matrix_2d_nurbs
    sll_real64, intent(in) :: eta1
    sll_real64, intent(in) :: eta2
    sll_real64             :: inv_j11
    sll_real64             :: inv_j12
    sll_real64             :: inv_j21
    sll_real64             :: inv_j22
    sll_real64             :: r_jac ! reciprocal of the jacobian
    
    SLL_ASSERT( eta1 <= 1.0_f64)
    SLL_ASSERT( eta1 >= 0.0_f64)
    SLL_ASSERT( eta2 <= 1.0_f64)
    SLL_ASSERT( eta2 >= 0.0_f64)

    r_jac = 1.0_f64/transf%jacobian( eta1, eta2 )


    if (transf%is_rational == 0) then ! IN the case of SPLINE
       
       inv_j11 = transf%x1_interp%interpolate_derivative_eta1( eta1, eta2 )
       inv_j12 = transf%x1_interp%interpolate_derivative_eta2( eta1, eta2 )
       inv_j21 = transf%x2_interp%interpolate_derivative_eta1( eta1, eta2 )
       inv_j22 = transf%x2_interp%interpolate_derivative_eta2( eta1, eta2 )
       
       inverse_jacobian_matrix_2d_nurbs(1,1) =  inv_j22*r_jac
       inverse_jacobian_matrix_2d_nurbs(1,2) = -inv_j12*r_jac
       inverse_jacobian_matrix_2d_nurbs(2,1) = -inv_j21*r_jac
       inverse_jacobian_matrix_2d_nurbs(2,2) =  inv_j11*r_jac
    
    else 
       
       inv_j11 = transf%x1_interp%interpolate_derivative_eta1( eta1, eta2 )&
            *transf%x3_interp%interpolate_value(eta1,eta2) &
            - transf%x3_interp%interpolate_derivative_eta1( eta1, eta2 )&
            *transf%x1_interp%interpolate_value(eta1,eta2)
       
       inv_j12 = transf%x1_interp%interpolate_derivative_eta2( eta1, eta2 )&
            *transf%x3_interp%interpolate_value(eta1,eta2) &
            - transf%x3_interp%interpolate_derivative_eta2( eta1, eta2 )&
            *transf%x1_interp%interpolate_value(eta1,eta2)
       
       inv_j21 = transf%x2_interp%interpolate_derivative_eta1( eta1, eta2 )&
            *transf%x3_interp%interpolate_value(eta1,eta2) &
            - transf%x3_interp%interpolate_derivative_eta1( eta1, eta2 )&
            *transf%x2_interp%interpolate_value(eta1,eta2)
       
       inv_j22 = transf%x2_interp%interpolate_derivative_eta2( eta1, eta2 )&
            *transf%x3_interp%interpolate_value(eta1,eta2) &
            - transf%x3_interp%interpolate_derivative_eta2( eta1, eta2 )&
            *transf%x2_interp%interpolate_value(eta1,eta2)
       
       
       inverse_jacobian_matrix_2d_nurbs(1,1) =  inv_j22*r_jac/&
            (transf%x3_interp%interpolate_value(eta1,eta2))**2
       inverse_jacobian_matrix_2d_nurbs(1,2) = -inv_j12*r_jac/&
            (transf%x3_interp%interpolate_value(eta1,eta2))**2
       inverse_jacobian_matrix_2d_nurbs(2,1) = -inv_j21*r_jac/&
            (transf%x3_interp%interpolate_value(eta1,eta2))**2
       inverse_jacobian_matrix_2d_nurbs(2,2) =  inv_j11*r_jac/&
            (transf%x3_interp%interpolate_value(eta1,eta2))**2
    end if

    
  end function inverse_jacobian_matrix_2d_nurbs
    

  subroutine write_to_file_2d_nurbs(transf,output_format)
    class(sll_coordinate_transformation_2d_nurbs) :: transf
    sll_int32, optional :: output_format 
    type(sll_logical_mesh_2d), pointer :: lm
    sll_int32           :: local_format 
    sll_real64, dimension(:,:), pointer :: x1mesh
    sll_real64, dimension(:,:), pointer :: x2mesh
    sll_int32  :: i1
    sll_int32  :: i2
    sll_int32  :: ierr
    sll_int32  :: file_id
    sll_int32  :: npts_eta1
    sll_int32  :: npts_eta2

    lm => transf%get_logical_mesh()

    npts_eta1  = lm%num_cells1 +1
    npts_eta2  = lm%num_cells2 +1


    if (.not. present(output_format)) then
       local_format = SLL_IO_XDMF
    else
       local_format = output_format
    end if


    if ( .not. transf%written ) then

       if (local_format == SLL_IO_XDMF) then
          SLL_ALLOCATE(x1mesh(npts_eta1,npts_eta2), ierr)
          SLL_ALLOCATE(x2mesh(npts_eta1,npts_eta2), ierr)
          do i1=1, npts_eta1
             do i2=1, npts_eta2
                x1mesh(i1,i2) = transf%x1_at_node(i1,i2)
                x2mesh(i1,i2) = transf%x2_at_node(i1,i2)
             end do
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
                x1mesh(i1,i2) = transf%x1_at_node(i1,i2)
                x2mesh(i1,i2) = transf%x2_at_node(i1,i2)
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

    SLL_DEALLOCATE(x1mesh,ierr)
    SLL_DEALLOCATE(x2mesh,ierr)
  end subroutine

  ! The coordinate transformation is reserving to itself the right to 
  ! delete the interpolators that were initially passed to it at 
  ! the moment of initialization. Technically, the transformation does 
  ! not 'own' the interpolators in the sense that it was not their 
  ! creator. However, since the interpolator is not an object that can 
  ! be shared amongst multiple interpolators, here we treat them as if
  ! indeed the coordinate transformation had created them. 
  ! This is a design choice that could be reversed in case that one 
  ! would wish to 'recycle' interpolators.

  subroutine delete_transformation_2d_nurbs( transf )
    class(sll_coordinate_transformation_2d_nurbs), intent(inout) :: transf

    transf%label = ""
    transf%written = .false.
    call transf%x1_interp%delete()
    call transf%x2_interp%delete()
    call transf%x3_interp%delete()
    nullify( transf%mesh2d_minimal)
    nullify( transf%mesh)
    
  end subroutine delete_transformation_2d_nurbs


end module sll_module_coordinate_transformations_2d_nurbs
