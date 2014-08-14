module finite_elements_solver_module
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  ! ******************************
  ! ON GOING WORK  !!!!!!!!!!!!!!!!
  ! FOR THE MOMENT I WILL BE COMENTING EVERYTHING I THINK
  ! IS UNNECESSARY OR OBSOLETE. ONCE THE DOCUMENT IS STABLE
  ! ALL THIS COMMENTS SHOULD BE ELEMINATED
  ! ******************************

  ! TODO @LM : see if all of these libraries are actually being used in stabilized version
  use sll_boundary_condition_descriptors
  use sll_module_scalar_field_2d_base
  use sll_module_scalar_field_2d_alternative
  use sll_arbitrary_degree_spline_interpolator_2d_module
  use connectivity_module
  use sll_knots
  use gauss_legendre_integration
  use gauss_lobatto_integration
  use sll_timer 
  use sll_sparse_matrix_module
  use sll_logical_meshes

  implicit none

  type :: finite_elements_solver

     ! Associated mesh (contains : eta1min, eta1max, delta1, ...)
     type(sll_logical_mesh_2d),  pointer :: mesh
     ! Total number of cells, typically for cartesian meshes nc1*nc2
     sll_int32 :: num_cells 
     sll_int32 :: num_cells_plus1

     ! Boundary conditions : THIS SHOULD BE CHANGED IN SLL AS A TYPE/LIST/... (?)
     sll_int32 :: bc_left
     sll_int32 :: bc_right
     sll_int32 :: bc_bottom
     sll_int32 :: bc_top

     ! Knots points
     sll_real64, dimension(:),   pointer :: knots1
     sll_real64, dimension(:),   pointer :: knots2
     ! Knots points for charge density
     sll_real64, dimension(:),   pointer :: knots1_source
     sll_real64, dimension(:),   pointer :: knots2_source

     ! Quadrature points/weights (gaussian, fekete, ...) of all the domain, by elements
     sll_int32 :: num_quad_pts
     sll_real64, dimension(:),   pointer :: quad_pts1
     sll_real64, dimension(:),   pointer :: quad_pts2
     sll_real64, dimension(:),   pointer :: quad_weight

     ! Degree of splines in each direction
     sll_int32 :: spline_degree1
     sll_int32 :: spline_degree2
     ! Number of local splines and total number of splines
     sll_int32 :: num_splines_loc
     sll_int32 :: num_splines_tot
     ! Global indexing of splines. The indexing of the splines
     ! in this array depends on the boundary conditions and it's
     ! independent of the mesh type.
     sll_int32, dimension(:), pointer :: global_spline_indices
     ! Local indexing of splines, independent of the mesh type
     sll_int32, dimension(:), pointer :: local_spline_indices
     sll_int32 :: local_spline_index1
     sll_int32 :: local_spline_index2
     ! Same as global_spline_indices but including 
     ! the changes resulting from the boundary conditions.
     sll_int32, dimension(:,:), pointer :: local_to_global_spline_indices
     sll_int32, dimension(:,:), pointer :: local_to_global_spline_indices_source
     ! TODO : what is this one for ? 
     sll_int32, dimension(:,:), pointer :: local_to_global_spline_indices_source_bis
     ! The following tables contain the values of splines (with degree <= deg) 
     ! in all quadrature points for the test function and source
     ! dimensions = (spacial, degree)
     sll_real64, dimension(:,:,:), pointer :: values_splines 
     sll_real64, dimension(:,:,:), pointer :: values_splines_source

     ! Table with the values of the jacobian of the transformation
     ! in all quadrature points
     sll_real64, dimension(:), pointer :: values_jacobian 

     ! Table with index of point to the left of the support
     sll_int32 , dimension(:)  , pointer :: tab_index_coeff

     ! Compressed sparse row matrix for test function (with constraint) and source
     type(sll_csr_matrix), pointer :: sll_csr_mat
     type(sll_csr_matrix), pointer :: sll_csr_mat_with_constraint
     type(sll_csr_matrix), pointer :: sll_csr_mat_source

     ! Initial data, useful variables and matrices (in global indexing)
     sll_real64, dimension(:), pointer :: source_vec
     sll_real64, dimension(:), pointer :: phi_vec
     sll_real64, dimension(:), pointer :: tmp_source_vec
     sll_real64, dimension(:), pointer :: tmp_phi_vec
     sll_real64, dimension(:), pointer :: masse
     sll_real64, dimension(:), pointer :: stiff
  
     ! Integral of the jacobian. ie sum of weighted jacobians 
     sll_real64 :: intjac

  end type finite_elements_solver

  ! For the integration mode.  
  sll_int32, parameter :: ES_GAUSS_LEGENDRE = 0, ES_GAUSS_LOBATTO = 1
  
  interface sll_delete
     module procedure delete_solver
  end interface sll_delete

  interface initialize
     module procedure initialize_finite_elements_solver
  end interface initialize
  
  
contains ! =============================================================



  function new_finite_elements_solver( &
       mesh, &
       spline_degree_eta1, &
       spline_degree_eta2, &
       quadrature_type1, &
       quadrature_type2, &
       bc_left, &
       bc_right, &
       bc_bottom, &
       bc_top ) result(solv)
    
    type(finite_elements_solver),          pointer :: solv
    type(sll_logical_mesh_2d), intent(in), pointer :: mesh
    sll_int32,  intent(in) :: spline_degree_eta1
    sll_int32,  intent(in) :: spline_degree_eta2
    sll_int32,  intent(in) :: bc_left
    sll_int32,  intent(in) :: bc_right
    sll_int32,  intent(in) :: bc_bottom
    sll_int32,  intent(in) :: bc_top
    sll_int32,  intent(in) :: quadrature_type1
    sll_int32,  intent(in) :: quadrature_type2
    sll_int32 :: ierr

    SLL_ALLOCATE(solv, ierr)
    call initialize( &
         solv, &
         mesh, &
         spline_degree_eta1, &
         spline_degree_eta2, &
         quadrature_type1, &
         quadrature_type2, &
         bc_left, &
         bc_right, &
         bc_bottom, &
         bc_top)
    
  end function new_finite_elements_solver



  subroutine initialize_finite_elements_solver( &
       solv, &
       mesh, &
       spline_degree_eta1, &
       spline_degree_eta2, &
       quadrature_type1, &
       quadrature_type2, &
       bc_left, &
       bc_right, &
       bc_bottom, &
       bc_top)
    
    type(finite_elements_solver), intent(out)         :: solv
    type(sll_logical_mesh_2d),    intent(in), pointer :: mesh
    sll_int32, intent(in) :: spline_degree_eta1
    sll_int32, intent(in) :: spline_degree_eta2
    sll_int32, intent(in) :: bc_left
    sll_int32, intent(in) :: bc_right
    sll_int32, intent(in) :: bc_bottom
    sll_int32, intent(in) :: bc_top
    sll_int32, intent(in) :: quadrature_type1
    sll_int32, intent(in) :: quadrature_type2
    sll_int32 :: knots1_size
    sll_int32 :: knots2_size
    sll_int32 :: num_splines1
    sll_int32 :: num_splines2
    sll_int32 :: vec_sz ! for source_vec and phi_vec allocations
    sll_int32 :: ierr,ierr1
    sll_int32 :: i, j
    sll_int32 :: nc1, nc2
    sll_int32 :: num_ele
    sll_int32 :: num_quad1
    sll_int32 :: num_quad2
    sll_int32 :: global_index
    sll_real64:: delta1
    sll_real64:: delta2
    ! Flag to notify is all boundary conditions are periodic
    sll_int32 :: sll_perper = 0 
    ! TODO @LM : these variables probably wont work or should be at least renamed
    sll_real64, dimension(:) :: temp_pts_wgh_1
    sll_real64, dimension(:) :: temp_pts_wgh_2

    ! Logical mesh : contains information as : eta1min, eta1max, delta1, ...
    solv%mesh => mesh

    ! Typically the total number of cells is just the product of the number
    ! of cells in each direction :
    solv%num_cells = mesh%num_cells1 * mesh%num_cells2
    solv%num_cells_plus1 = (mesh%num_cells1 + 1) * (mesh%num_cells2 + 1)

    ! Splines degrees in each direction
    solv%spline_degree1 = spline_degree_eta1
    solv%spline_degree2 = spline_degree_eta2

    ! This should be changed to verify that the passed BC's are part of the
    ! recognized list described in sll_boundary_condition_descriptors...
    solv%bc_left   = bc_left
    solv%bc_right  = bc_right
    solv%bc_bottom = bc_bottom
    solv%bc_top    = bc_top
        
    ! Allocate and fill the quadrature points/weights information --------------
    ! in both directions ------------------------------------------------- BEGIN
    num_quad1 = (spline_degree_eta1+1) 
    num_quad2 = (spline_degree_eta2+1)
    solv%num_quad_pts = num_quad1 * num_quad2 * solv%num_cells
    SLL_ALLOCATE(solv%quad_pts1  (num_quad_pts),ierr)
    SLL_ALLOCATE(solv%quad_pts2  (num_quad_pts),ierr)
    SLL_ALLOCATE(solv%quad_weight(num_quad_pts),ierr)
    ! Initialization :
    solv%quad_pts1(:)   = 0.0_f64
    solv%quad_pts2(:)   = 0.0_f64
    solv%quad_weight(:) = 0.0_f64
   
    select case(quadrature_type1)
    case (ES_GAUSS_LEGENDRE)
       temp_pts_wgh_1 = gauss_legendre_points_and_weights(num_quad1)
    case (ES_GAUSS_LOBATTO)
       temp_pts_wgh_1 = gauss_lobatto_points_and_weights(num_quad2)
    case DEFAULT
       print *, 'new_finite_elements_solver(): have not type of gauss points in the first direction'
    end select
    
    select case(quadrature_type2)
    case (ES_GAUSS_LEGENDRE)
       temp_pts_wgh_2 = gauss_legendre_points_and_weights(num_quad2)
    case (ES_GAUSS_LOBATTO)
       temp_pts_wgh_2 = gauss_lobatto_points_and_weights(num_quad2)
    case DEFAULT
       print *, 'new_finite_elements_solver(): have not type of gauss points in the second direction'
    end select

    delta1 = mesh%delta_eta1
    delta2 = mesh%delta_eta2
    ! Loop over cells/elements on both directions
    do nc2 =1,mesh%num_cells2
       do nc1 =1,mesh%num_cells1
          num_ele = nc1 + (nc2-1) * mesh%num_cells2  ! Global index of elements
          eta1 = mesh%eta1_min + delta_eta1 * nc1 ! eta1 at lower-left cell corner
          eta2 = mesh%eta2_min + delta_eta2 * nc2 ! eta2 at lower-left cell corner
          ! Loop over quadrature points of every cells
          do j=1, num_quad2
             do i=1, num_quad1

                ! We compute the global index for quadrature point at cell num_ele
                global_index = i+(j-1)*num_quad2 + (num_ele-1)*num_quad1*num_quad2

                ! Rescaling of quadrature point to [eta1,eta1+delta1]x[eta2,eta2+delta2]
                solv%quad_pts1(global_index)    = &
                     eta1 + 0.5_f64 * delta1 * (1._f64 + temp_pts_wgh1(1, i))
                solv%quad_pts2(global_index)    = &
                     eta2 + 0.5_f64 * delta2 * (1._f64 + temp_pts_wgh1(1, j))

                ! We get the quadrature weight and scale it
                solv%quad_weight(global_index) = &
                     0.5_f64 * 0.5_f64 * delta1 * delta2 * &
                     temp_pts_wgh1(2, i) * temp_pts_wgh2(2, j)

                solv%local_spline_index1 = solv%spline_degree1 + 1

                ! Boundary condition evaluation 
                if ((solv%bc_left==SLL_PERIODIC).and.(solv%bc_right==SLL_PERIODIC))then
                   ! rescale gauss point in interval [0,delta2]
                   solv%quad_pts1(global_index)  = &
                        0.5_f64 * delta2 * (solv%quad_pts1(global_index) + 1.0_f64)
                   solv%local_spline_index1 = solv%spline_degree1 + 1
                end if
                if ((solv%bc_bottom==SLL_PERIODIC).and.(solv%bc_top==SLL_PERIODIC))then
                   ! rescale gauss point in interval [0,delta2]
                   solv%quad_pts2(global_index)  = &
                        0.5_f64 * delta2 * (solv%quad_pts2(global_index) + 1.0_f64)
                   solv%local_spline_index1 = solv%spline_degree1 + 1
                end if

             
             end do
          end do
       end do
    end do
    !  ---------------------------------------------- END QUADRATURE POINTS INIT
    ! --------------------------------------------------------------------------


    ! --------------------------------------------------------------------------
    ! BEGIN ALLOCATION AND INITIALIZATION OF SPLINES KNOTS ---------------------

    if( (bc_left == SLL_PERIODIC) .and. (bc_right == SLL_PERIODIC) .and. &
         (bc_bottom == SLL_PERIODIC) .and. (bc_top == SLL_PERIODIC) ) then
       solv%num_splines_tot = solv%num_cells
       knots1_size = 2 * spline_degree_eta1 + 2
       knots2_size = 2 * spline_degree_eta2 + 2
       vec_sz      = solv%num_cells
       sll_perper  = 1 
    else if( (bc_left == SLL_PERIODIC) .and. (bc_right == SLL_PERIODIC) .and.&
         (bc_bottom == SLL_DIRICHLET) .and. (bc_top == SLL_DIRICHLET) ) then
       solv%num_splines_tot = mesh%num_cells1 * &
            (mesh%num_cells2 + spline_degree_eta2 - 2)
       knots1_size = 2*spline_degree_eta1 + 2
       knots2_size = 2*spline_degree_eta2 + mesh%num_cells2+ 1
       vec_sz      = mesh%num_cells1 * (mesh%num_cells2 + spline_degree_eta2)
    else if( (bc_left == SLL_DIRICHLET) .and. (bc_right == SLL_DIRICHLET) .and.&
         (bc_bottom == SLL_PERIODIC) .and. (bc_top == SLL_PERIODIC) ) then
       solv%num_splines_tot = (mesh%num_cells1 + spline_degree_eta1 - 2) * &
            mesh%num_cells2 
       knots1_size = 2*spline_degree_eta1 + mesh%num_cells1+1
       knots2_size = 2*spline_degree_eta2 + 2
       vec_sz      = (mesh%num_cells1 + spline_degree_eta1) * mesh%num_cells2
    else if( (bc_left == SLL_DIRICHLET) .and. (bc_right == SLL_DIRICHLET) .and.&
         (bc_bottom == SLL_DIRICHLET) .and. (bc_top == SLL_DIRICHLET) ) then
       solv%num_splines_tot = (mesh%num_cells1 + spline_degree_eta1 - 2) * &
            (mesh%num_cells2 + spline_degree_eta2 - 2)
       knots1_size = 2*spline_degree_eta1 + num_cells_eta1 + 1
       knots2_size = 2*spline_degree_eta2 + num_cells_eta2 + 1
       vec_sz      = (mesh%num_cells1 + spline_degree_eta1) * &
            (mesh%num_cells2 + spline_degree_eta2)
    end if

    SLL_ALLOCATE(solv%knots1(knots1_size),ierr1)
    SLL_ALLOCATE(solv%knots2(knots2_size),ierr)

    call initialize_knots( &
         spline_degree_eta1, &
         mesh%num_cells1, &
         eta1_min, &
         eta1_max, &
         bc_left, &
         bc_right, &
         solv%knots1 )
    
    call initialize_knots( &
         spline_degree_eta2, &
         mesh%num_cells2, &
         eta2_min, &
         eta2_max, &
         bc_bottom, &
         bc_top, &
         solv%knots2 )

    ! Knots for source term : 
    SLL_ALLOCATE(solv%knots1_source(mesh%num_cells1 + spline_degree_eta1 + 2),ierr1)
    SLL_ALLOCATE(solv%knots2_source(mesh%num_cells2 + spline_degree_eta2 + 2 ),ierr)

    ! ----------------------- END ALLOCATION AND INITIALIZATION OF SPLINES KNOTS
    ! --------------------------------------------------------------------------



    ! --------------------------------------------------------------------------
    ! ------------------- BEGIN ALLOCATION AND INITIALIZATION OF BASIS FUNCTIONS

    ! Number of local splines :
    solv%num_splines_loc = (spline_degree_eta1+1)*(spline_degree_eta2+1)
    
    ! The total number of splines in a single direction is given by
    ! num_cells + spline_degree. For 2D is just a product of both directions
    num_splines1 = mesh%num_cells1 + spline_degree_eta1
    num_splines2 = mesh%num_cells2 + spline_degree_eta2
    SLL_ALLOCATE(solv%global_spline_indices(num_splines1*num_splines2),ierr)
    solv%global_spline_indices(:) = 0
    
    ! Local indexing of the splines 
    SLL_ALLOCATE(solv%local_spline_indices(solv%num_splines_loc * solv%num_cells),ierr)
    solv%local_spline_indices(:) = 0
    
    ! Connectivity between local and global indexing systems
    SLL_ALLOCATE(solv%local_to_global_spline_indices(solv%num_splines_loc, solv%num_cells), ierr)
    solv%local_to_global_spline_indices = 0
    ! Same for source
    SLL_ALLOCATE(solv%local_to_global_spline_indices_source(solv%num_splines_loc, solv%num_cells),ierr)
    solv%local_to_global_spline_indices_source = 0
    ! Same for source (bis)
    SLL_ALLOCATE(solv%local_to_global_spline_indices_source_bis(solv%num_splines_loc, solv%num_cells),ierr)
    solv%local_to_global_spline_indices_source_bis = 0

    ! TODO @LM : spline initialization here
    call initialize_basis_functions(solv)

    ! --------------------- END ALLOCATION AND INITIALIZATION OF BASIS FUNCTIONS
    ! --------------------------------------------------------------------------

    
    SLL_ALLOCATE(solv%source_vec(vec_sz),ierr)
    SLL_ALLOCATE(solv%masse  (vec_sz),ierr)
    SLL_ALLOCATE(solv%stiff  (vec_sz),ierr)
    SLL_ALLOCATE(solv%phi_vec(num_splines_tot),ierr)
    
    if(sll_perper == 1) then
       SLL_ALLOCATE(solv%tmp_source_vec(num_splines_tot + 1),ierr)
       SLL_ALLOCATE(solv%tmp_phi_vec(num_splines_tot + 1),ierr)
    else
       SLL_ALLOCATE(solv%tmp_source_vec(num_splines_tot),ierr)
       SLL_ALLOCATE(solv%tmp_phi_vec(num_splines_tot),ierr)
    endif

    solv%source_vec(:) = 0.0_f64
    solv%phi_vec(:) = 0.0_f64
    solv%masse(:)   = 0.0_f64
    solv%stiff(:)   = 0.0_f64
    solv%intjac     = 0.0_f64
    
    call initconnectivity( &
         mesh%num_cells1, &
         mesh%num_cells2, &
         spline_degree_eta1, &
         spline_degree_eta2, &
         bc_left, &
         bc_right, &
         bc_bottom, &
         bc_top, &
         solv%local_spline_indices, &
         solv%global_spline_indices, &
         solv%local_to_global_spline_indices )
    
    solv%sll_csr_mat => new_csr_matrix( &
         num_splines_tot, &
         num_splines_tot, &
         solv%num_cells, &
         solv%local_to_global_spline_indices, &
         solv%num_splines_loc, &
         solv%local_to_global_spline_indices, &
         solv%num_splines_loc)
    
    solv%knots1_source(1:spline_degree_eta1+1) = mesh%eta1_min
    solv%knots1_source(mesh%num_cells1+2:mesh%num_cells1+1+spline_degree_eta1+1) = &
         mesh%eta1_max
    
    if (mod(spline_degree_eta1 + 1, 2) == 0) then
       do i = spline_degree_eta1 +2, mesh%num_cells1 + 1
          solv%knots1_source(i) = eta1_min +  &
               (i - (spline_degree_eta1 +1)/2 - 1) * solv%mesh%delta_eta1  
       end do
    else
       do i = spline_degree_eta1 + 2, mesh%num_cells1 + 1
          solv%knots1_source ( i ) = &
               0.5_f64*(eta1_min + (i - (spline_degree_eta1)/2 - 1) * solv%mesh%delta_eta1 + &
               eta1_min +  ( i -1 - (spline_degree_eta1)/2 -1) * solv%mesh%delta_eta1 )
       end do
    end if
    
    solv%knots2_source(1:spline_degree_eta2+1) = eta2_min
    solv%knots2_source(mesh%num_cells2+2:mesh%num_cells2+1+spline_degree_eta2+1) = eta2_max
    
    if ( mod(spline_degree_eta2 +1,2) == 0 ) then
       do i = spline_degree_eta2 +1 + 1, mesh%num_cells2 + 1
          solv%knots2_source( i ) = eta2_min + &
               ( i - (spline_degree_eta2 +1)/2-1 )*solv%mesh%delta_eta2 
       end do
    else     
       do i = spline_degree_eta2 +1 + 1, mesh%num_cells2 + 1
          solv%knots2_source ( i ) = &
               0.5*( eta2_min + ( i - (spline_degree_eta2)/2 -1)*solv%mesh%delta_eta2 + &
               eta2_min +  ( i -1 - (spline_degree_eta2)/2 -1)*solv%mesh%delta_eta2 )  
       end do
    end if
    
    
    ! ------------------------------------------------------------------------------
    ! BEGIN Allocation of basis splines related tables -----------------------------
    ! values_splines : table containning all values --------------------------------
    ! of basis splines in each direction in each quadrature point ------------------

    ! The third dimension contains the combination of the product between
    ! the derivatives of the basis functions
    SLL_ALLOCATE(solv%values_splines(solv%num_cells * solv%num_splines_loc, 4), ierr)
    solv%values_splines = 0.0_f64
    SLL_ALLOCATE(solv%values_splines_source(solv%num_cells * solv%num_splines_loc, 4), ierr)
    solv%values_splines_source = 0.0_f64
    SLL_ALLOCATE(solv%values_jacobian(solv%num_cells * solv%num_splines_loc), ierr)
    solv%values_jacobian = 0.0_f64
    SLL_ALLOCATE(solv%tab_index_coeff(solv%num_cells * solv%num_splines_loc),ierr)

    ! ------------------------------------------------ END ALLOCATION SPLINES TABLES
    ! ------------------------------------------------------------------------------
 
  end subroutine initialize_finite_elements_solver
  


  subroutine initialize_basis_functions(solv)

    type(finite_elements_solver), intent(inout) :: solv

    ! These objects will contain the value of the basis function derivatives
    sll_real64, dimension(solv%spline_degree1+1,2) :: basis_deriv_x1
    sll_real64, dimension(solv%spline_degree2+1,2) :: basis_deriv_x2
    sll_real64, dimension(solv%spline_degree1+1,2) :: basis_deriv_x1_source
    sll_real64, dimension(solv%spline_degree1+1,2) :: basis_deriv_x2_source
    ! Intermediate variables
    sll_real64, dimension(solv%spline_degree1+1,solv%spline_degree1+1) :: work1
    sll_real64, dimension(solv%spline_degree2+1,solv%spline_degree2+1) :: work2
    ! Variables containing the index of the begining of the local spline
    sll_int32  :: local_spline_index1
    sll_int32  :: local_spline_index2
    ! Counters variables
    sll_int32  :: num_ele
    sll_int32  :: local_index
    sll_int32  :: global_index
    sll_int32  :: cell_i
    sll_int32  :: cell_j
    ! Index of point to the left of the support of the spline 
    sll_int32 :: left_x,left_y
    sll_int32 :: mflag_x, mflag_y ! error flags

    ! Initialization
    work1(:,:)                 = 0.0_f64
    work2(:,:)                 = 0.0_f64
    basis_deriv_x1(:,:)        = 0.0_f64
    basis_deriv_x2(:,:)        = 0.0_f64
    basis_deriv_x1_source(:,:) = 0.0_f64
    basis_deriv_x2_source(:,:) = 0.0_f64

    cell_i = 1
    cell_j = 1

    do num_ele = 1,solv%num_cells
       do local_index = 1, solv%num_splines_loc

          ! TODO @LM add function local to global here !
          global_index = (num_ele - 1) * solv%num_splines_loc + local_index
          
          qpt1  = solv%quad_pts1(global_index)
          qpt2  = solv%quad_pts2(global_index)
          wqpt  = solv%quad_weight(global_index)
          
          ! TODO @LM : careful here !! this formula is probably wrong
          ! and should probably depend on the BC
          local_spline_index = local_index
          
          ! Computing the values of the first 2 derivatives on the point qpt2
          ! Result is saved on basis_deriv_x2
          call bsplvd( &
               solv%knots2, &
               solv%spline_degree2+1,&
               qpt2,&
               local_spline_index2,&
               work2,&
               basis_deriv_x2,&
               2)
             
          ! Computing the spline interval, returns left_y, inex of the 
          ! knot point to the left of the basis function's support
          call interv(&
               solv%knots2_source,&
               solv%mesh%num_cells2 + solv%spline_degree2+ 2,& !TODO @LM
               qpt2, &
               left_y, &
               mflag_y )
             
          !TODO @LM : Add test on mflag_y
          
          ! We do the same work for the source
          call bsplvd( &
               solv%knots2_source, &
               solv%spline_degree2+1,&
               qpt2,&
               left_y,&
               work2,&
               basis_deriv_x2_source,&
               2)
          
          ! And again, we do the same work for the other direction (x1)
          call bsplvd(&
               solv%knots1,&
               solv%spline_degree1+1,&
               qpt1,&
               local_spline_index1,&
               work1,&
               basis_deriv_x1,&
               2 )
          
          call interv ( solv%knots1_source, &
               size(solv%knots1_source), &
               qpt1, &
               left_x, &
               mflag_x )
          
          !TODO @LM : Add test on mflag_y
          
          ! And the source term in the direction x1
          call bsplvd(&
               solv%knots1_source,&
               solv%spline_degree1+1,&
               qpt1,&
               left_x,&
               work1,&
               basis_deriv_x1_source,&
               2 )
          
          
          ! If we note N1 and N2 the basis functions
          ! we stock the values of : N1N2, N1'N2, N1N2', N1'N2'
          solv%values_splines(global_index, 1) = &
               basis_deriv_x1(i,1) * basis_deriv_x2(j,1)
          solv%values_splines(global_index, 2) = &
               basis_deriv_x1(i,2) * basis_deriv_x2(j,1)
          solv%values_splines(global_index, 3) = &
               basis_deriv_x1(i,1) * basis_deriv_x2(j,2)
          solv%values_splines(global_index, 4) = &
               basis_deriv_x1(i,2) * basis_deriv_x2(j,2)
             
          ! IDEM for source distribution
          solv%values_splines_source(global_index, 1) = &
               basis_deriv_x1_source(i,1) * basis_deriv_x2_source(j,1)
          solv%values_splines_source(global_index, 2) = &
               basis_deriv_x1_source(i,2) * basis_deriv_x2_source(j,1)
          solv%values_splines_source(global_index, 3) = &
               basis_deriv_x1_source(i,1) * basis_deriv_x2_source(j,2)
          solv%values_splines_source(global_index, 4) = &
               basis_deriv_x1_source(i,2) * basis_deriv_x2_source(j,2)
          
          ! Keeping the index of coefficients
          solv%tab_index_coeff(global_index) = left_x + left_y * solv%mesh%num_cells2
          
       end do
    end do
  end subroutine initialize_basis_functions


  subroutine factorize_mat_solv(&
       solv, &
       a11_field_mat, &
       a12_field_mat,&
       a21_field_mat,&
       a22_field_mat,&
       b1_field_vect,&
       b2_field_vect,&
       c_field)

    type(finite_elements_solver),intent(inout) :: solv
    class(sll_scalar_field_2d_base), pointer   :: a11_field_mat
    class(sll_scalar_field_2d_base), pointer   :: a12_field_mat
    class(sll_scalar_field_2d_base), pointer   :: a21_field_mat
    class(sll_scalar_field_2d_base), pointer   :: a22_field_mat
    class(sll_scalar_field_2d_base), pointer   :: b1_field_vect
    class(sll_scalar_field_2d_base), pointer   :: b2_field_vect
    class(sll_scalar_field_2d_base), pointer   :: c_field
    sll_real64, dimension(:,:), allocatable :: M_c_loc
    sll_real64, dimension(:,:), allocatable :: K_a11_loc
    sll_real64, dimension(:,:), allocatable :: K_a12_loc
    sll_real64, dimension(:,:), allocatable :: K_a21_loc
    sll_real64, dimension(:,:), allocatable :: K_a22_loc
    sll_real64, dimension(:,:), allocatable :: M_b_vect_loc
    sll_real64, dimension(:,:), allocatable :: S_b1_loc
    sll_real64, dimension(:,:), allocatable :: S_b2_loc  
    sll_real64, dimension(:),   allocatable :: Masse_loc
    sll_real64, dimension(:),   allocatable :: Stiff_loc
    sll_real64, dimension(:,:,:), pointer   :: Source_loc
    sll_int32 :: num_splines_loc
    sll_int32 :: ierr
    sll_int32 :: i
    sll_int32 :: j
    sll_int32 :: cell_index
    sll_int32 :: bc_left
    sll_int32 :: bc_right
    sll_int32 :: bc_bottom
    sll_int32 :: bc_top
    sll_int32 :: sll_perper
    type(sll_time_mark)  :: t0 
    double precision :: time
    character(len=*),parameter :: as_file1='mat'
    
    !call sll_set_time_mark(t0)
    
    bc_left   = solv%bc_left
    bc_right  = solv%bc_right
    bc_bottom = solv%bc_bottom
    bc_top    = solv%bc_top

    num_splines_loc = solv%num_splines_loc
    if( (bc_left == SLL_PERIODIC) .and. (bc_right == SLL_PERIODIC) .and. &
       (bc_bottom == SLL_PERIODIC) .and. (bc_top == SLL_PERIODIC) ) then
       sll_perper = 0
    else
       sll_perper = 1  
    end if   

    SLL_ALLOCATE(Source_loc(solv%num_cells,num_splines_loc,num_splines_loc),ierr)
    SLL_ALLOCATE(M_c_loc(num_splines_loc,num_splines_loc),ierr)
    SLL_ALLOCATE(K_a11_loc(num_splines_loc,num_splines_loc),ierr)
    SLL_ALLOCATE(K_a12_loc(num_splines_loc,num_splines_loc),ierr)
    SLL_ALLOCATE(K_a21_loc(num_splines_loc,num_splines_loc),ierr)
    SLL_ALLOCATE(K_a22_loc(num_splines_loc,num_splines_loc),ierr)
    SLL_ALLOCATE(S_b1_loc(num_splines_loc,num_splines_loc),ierr)
    SLL_ALLOCATE(S_b2_loc(num_splines_loc,num_splines_loc),ierr)
    SLL_ALLOCATE(M_b_vect_loc(num_splines_loc,num_splines_loc),ierr)
    SLL_ALLOCATE(Masse_loc(num_splines_loc),ierr)
    SLL_ALLOCATE(Stiff_loc(num_splines_loc),ierr)

    Masse_loc(:) = 0.0_f64
    Stiff_loc(:) = 0.0_f64
    Source_loc(:,:,:) = 0.0_f64

    do cell_index=1,solv%num_cells

       call build_local_matrices( &
            solv, &
            cell_index,&
            a11_field_mat, &
            a12_field_mat, &
            a21_field_mat, &
            a22_field_mat, &
            b1_field_vect,&
            b2_field_vect,&
            c_field, &
            Masse_loc,&
            Stiff_loc,&
            M_c_loc, &
            K_a11_loc, &
            K_a12_loc, &
            K_a21_loc, &
            K_a22_loc, &
            M_b_vect_loc, &
            S_b1_loc,  &
            S_b2_loc,&
            Source_loc)
       
       call local_to_global_matrices( &
            solv, &
            cell_index, &
            i, &
            j, &
            Masse_loc,&
            Stiff_loc,&
            M_c_loc, &
            K_a11_loc, &
            K_a12_loc, &
            K_a21_loc, &
            K_a22_loc,&
            M_b_vect_loc, &
            S_b1_loc, &
            S_b2_loc, &
            solv%masse,&
            solv%stiff)
       
    end do
    
    if (sll_perper == 0) then
       solv%sll_csr_mat_with_constraint => new_csr_matrix_with_constraint(solv%sll_csr_mat)  
       call csr_add_one_constraint( &
            solv%sll_csr_mat%opi_ia, & 
            solv%sll_csr_mat%opi_ja, &
            solv%sll_csr_mat%opr_a, &
            solv%sll_csr_mat%num_rows, &
            solv%sll_csr_mat%num_nz, &
            solv%masse, &
            solv%sll_csr_mat_with_constraint%opi_ia, &
            solv%sll_csr_mat_with_constraint%opi_ja, &
            solv%sll_csr_mat_with_constraint%opr_a)  
       call sll_factorize_csr_matrix(solv%sll_csr_mat_with_constraint)
    else
       call sll_factorize_csr_matrix(solv%sll_csr_mat)      
    end if
    
    solv%sll_csr_mat_source => new_csr_matrix( &
         size(solv%masse,1), &
         ! Second parameter should be : TODO @LM
         ! (solv%num_cells1+1)*(solv%num_cells2+1),&
         solv%num_cells_plus1,&
         solv%num_cells, &
         solv%local_to_global_spline_indices_source_bis, &
         solv%num_splines_loc, &
         solv%local_to_global_spline_indices_source, &
         solv%num_splines_loc )

    
    call compute_Source_matrice(solv,Source_loc)
    
    SLL_DEALLOCATE_ARRAY(Source_loc,ierr)
    SLL_DEALLOCATE_ARRAY(M_c_loc,ierr)
    SLL_DEALLOCATE_ARRAY(K_a11_loc,ierr)
    SLL_DEALLOCATE_ARRAY(K_a12_loc,ierr)
    SLL_DEALLOCATE_ARRAY(K_a21_loc,ierr)
    SLL_DEALLOCATE_ARRAY(K_a22_loc,ierr)
    SLL_DEALLOCATE_ARRAY(M_b_vect_loc,ierr)
    SLL_DEALLOCATE_ARRAY(S_b1_loc,ierr)
    SLL_DEALLOCATE_ARRAY(S_b2_loc,ierr)
    SLL_DEALLOCATE_ARRAY(Stiff_loc,ierr) 
    SLL_DEALLOCATE_ARRAY(Masse_loc,ierr) 
   
    !time = sll_time_elapsed_since(t0)
    !print*, '#time for factorize_mat_solv', time

  end subroutine factorize_mat_solv
  

  
  subroutine solve_general_coordinates_elliptic_eq(&
       solv,&
       source,&
       phi)
    use sll_timer
    class(finite_elements_solver) :: solv
    class(sll_scalar_field_2d_discrete_alt), intent(inout)  :: phi
    class(sll_scalar_field_2d_base), intent(in),target  :: source
    sll_int32 :: i
    sll_int32 :: j,ierr
    sll_int32 :: cell_index
    sll_int32 :: num_splines_loc
    sll_real64 :: int_source,int_jac
    sll_real64, dimension(:), allocatable   :: M_source_loc
    sll_real64, dimension(:), allocatable   :: source_at_quad
    sll_int32 :: num_pts, ig1, ig2, ig, jg
    sll_real64 :: val_jac
    class(sll_scalar_field_2d_base),pointer  :: base_field_pointer
    class(sll_interpolator_2d_base),pointer  :: base_interpolator_pointer
    sll_real64, dimension(:,:), pointer :: coeff_source
    sll_real64, dimension(:), pointer :: source_coeff_1d
    
    num_splines_loc = solv%num_splines_loc
    SLL_ALLOCATE(M_source_loc(num_splines_loc),ierr)
    

    SLL_ALLOCATE(source_at_quad(num_quad_pts),ierr)
    source_at_quad(:) = 0.0_f64
    SLL_ALLOCATE(source_coeff_1d(solv%num_cells_plus1),ierr)
    
    M_source_loc(:)    = 0.0_f64
    solv%source_vec(:) = 0.0_f64
    source_coeff_1d(:) = 0.0_f64
   
  
    !call sll_set_time_mark(t0)
    !ES Compute source at all quad points
    ig1 = 0 
    ig2 = 0 
    int_source = 0.0_f64
    int_jac = 0.0_f64
    
    ! afin d'optimiser la construction de la matrice source
    ! on va proceder de la facon qui suit 
    base_field_pointer => source
    select type( type_field => base_field_pointer)
    class is (sll_scalar_field_2d_discrete_alt)
       base_interpolator_pointer => type_field%interp_2d
       select type( type_interpolator => base_interpolator_pointer)
       class is (arb_deg_2d_interpolator)
          coeff_source => type_interpolator%get_coefficients()
          ! TODO @LM : this whole part should be re written, quick fix though:
          ! put the spline coefficients in a 1d array
          do j=1,size(coeff_source, 2)
             do i=1,size(coeff_source, 1)
                source_coeff_1d(i+size(coeff_source,2)*(j-1)) = coeff_source(i,j)
             end do
          end do

          call sll_mult_csr_matrix_vector(&
               solv%sll_csr_mat_source,&
               source_coeff_1d,solv%source_vec)

          if( ((solv%bc_bottom==SLL_PERIODIC).and.(solv%bc_top==SLL_PERIODIC)) &
               .and.((solv%bc_left==SLL_PERIODIC).and.(solv%bc_right==SLL_PERIODIC)) )then
             
             solv%source_vec = solv%source_vec - sum(solv%source_vec)/solv%intjac*solv%masse
          end if
          
       class DEFAULT
          
       do quad_index=1,solv%num_quad_pts ! Loop over quadrature points
          qpt1   = solv%quad_pts1  (quad_index)
          qpt2   = solv%quad_pts2  (quad_index)
          weight = solv%quad_weight(quad_index)
          source_at_quad(quad_index) = source%value_at_point(qpt1,qpt2)
          val_jac = solv%values_jacobian(quad_index)
          int_source = int_source + source%value_at_point(qpt1,qpt2)*weight*val_jac 
          int_jac = int_jac + weight*val_jac
       end do

       if( ((solv%bc_bottom==SLL_PERIODIC).and.(solv%bc_top==SLL_PERIODIC)) &
            .and. ((solv%bc_left==SLL_PERIODIC).and.(solv%bc_right==SLL_PERIODIC)) )then   
          source_at_quad = source_at_quad - int_source/int_jac
       end if
          
       do cell_index=1,solv%num_cells ! Loop over elements
          call build_local_matrices_source( &
               solv, &
               cell_index, &
               source_at_quad, &
               M_source_loc)
                
          call local_to_global_matrices_source( &
               solv, &
               cell_index, &
               M_source_loc)
       end do
    end select
       
    class is (sll_scalar_field_2d_analytic_alt)
       
    do quad_index=1,solv%num_quad_pts ! Loop over quadrature points
       qpt1   = solv%quad_pts1  (quad_index)
       qpt2   = solv%quad_pts2  (quad_index)
       weight = solv%quad_weight(quad_index)
       source_at_quad(quad_index) = source%value_at_point(qpt1,qpt2)
       val_jac = solv%values_jacobian(quad_index)
       int_source = int_source + source%value_at_point(qpt1,qpt2) * weight * val_jac 
       int_jac = int_jac + weight * val_jac
    end do
    if( ((solv%bc_bottom==SLL_PERIODIC).and.(solv%bc_top==SLL_PERIODIC)) &
         .and. ((solv%bc_left==SLL_PERIODIC).and.(solv%bc_right==SLL_PERIODIC)) )then
       source_at_quad = source_at_quad - int_source/int_jac
    end if
          
    do cell_index = 1, solv%num_cells
       call build_local_matrices_source( &
            solv, &
            cell_index, &
            source_at_quad, &
            M_source_loc)
             
       call local_to_global_matrices_source( &
            solv, &
            cell_index, &
            M_source_loc)
    end do
    end select
    !time = sll_time_elapsed_since(t0)

    call solve_linear_system(solv)
    
    call  phi%interp_2d%set_coefficients( solv%phi_vec)

    SLL_DEALLOCATE_ARRAY(M_source_loc,ierr)
    SLL_DEALLOCATE_ARRAY(source_at_quad,ierr)
  end subroutine solve_general_coordinates_elliptic_eq
  
  ! This is based on the assumption that all the input fields have the same
  ! boundary conditions. TO DO: put all the boundary condition parameters in
  ! a single module called 'boundary_condition_convention' or something, which
  ! can be used library-wide, this way we could extract this information 
  ! directly from the fields without any difficulties. 

  subroutine build_local_matrices( &
       solv, &
       cell_index,&
       a11_field_mat, &
       a12_field_mat, &
       a21_field_mat, &
       a22_field_mat, &
       b1_field_vect, &
       b2_field_vect, &
       c_field, &
       Masse_loc,&
       Stiff_loc,&
       M_c_loc, &
       K_a11_loc, &
       K_a12_loc, &
       K_a21_loc, &
       K_a22_loc, &
       M_b_vect_loc, &
       S_b1_loc,  &
       S_b2_loc,&
       Source_loc)

    
    class(finite_elements_solver) :: solv
    sll_int32, intent(in)         :: cell_index
    class(sll_scalar_field_2d_base), pointer :: a11_field_mat
    class(sll_scalar_field_2d_base), pointer :: a12_field_mat
    class(sll_scalar_field_2d_base), pointer :: a21_field_mat
    class(sll_scalar_field_2d_base), pointer :: a22_field_mat
    class(sll_scalar_field_2d_base), pointer :: b1_field_vect
    class(sll_scalar_field_2d_base), pointer :: b2_field_vect
    class(sll_scalar_field_2d_base), pointer :: c_field
    sll_real64, dimension(:,:,:), intent(inout) :: Source_loc
    sll_real64, dimension(:,:),   intent(out)   :: M_c_loc
    sll_real64, dimension(:,:),   intent(out)   :: K_a11_loc
    sll_real64, dimension(:,:),   intent(out)   :: K_a12_loc
    sll_real64, dimension(:,:),   intent(out)   :: K_a21_loc
    sll_real64, dimension(:,:),   intent(out)   :: K_a22_loc
    sll_real64, dimension(:,:),   intent(out)   :: M_b_vect_loc
    sll_real64, dimension(:,:),   intent(out)   :: S_b1_loc
    sll_real64, dimension(:,:),   intent(out)   :: S_b2_loc
    sll_real64, dimension(:),     intent(out)   :: Masse_loc
    sll_real64, dimension(:),     intent(out)   :: Stiff_loc

    ! Quadrature points coordinates and associated weights :
    sll_real64 :: qpt1
    sll_real64 :: qpt2
    sll_real64 :: wqpt
    ! matrix initializers
    sll_real64 :: val_c
    sll_real64 :: val_a11
    sll_real64 :: val_a12
    sll_real64 :: val_a21
    sll_real64 :: val_a22
    sll_real64 :: val_b1=0
    sll_real64 :: val_b1_der1=0
    sll_real64 :: val_b1_der2=0
    sll_real64 :: val_b2=0
    sll_real64 :: val_b2_der1=0
    sll_real64 :: val_b2_der2=0
    ! Jacobian variables
    sll_real64, dimension(2,2) :: jac_mat
    sll_real64 :: val_jac
    ! Matrices
    sll_real64 :: B11
    sll_real64 :: B12
    sll_real64 :: B21
    sll_real64 :: B22
    sll_real64 :: MC
    sll_real64 :: C1
    sll_real64 :: C2
    ! Variables for splines index
    sll_int32 :: left_x,left_y
    sll_int32 :: mflag_x, mflag_y
    ! Loop variables
    sll_int32  :: index1
    sll_int32  :: index2    


    Masse_loc(:)      = 0.0_f64
    Stiff_loc(:)      = 0.0_f64
    M_c_loc(:,:)      = 0.0_f64
    K_a11_loc(:,:)    = 0.0_f64
    K_a12_loc(:,:)    = 0.0_f64
    K_a21_loc(:,:)    = 0.0_f64
    K_a22_loc(:,:)    = 0.0_f64
    M_b_vect_loc(:,:) = 0.0_f64
    S_b1_loc(:,:)     = 0.0_f64
    S_b2_loc(:,:)     = 0.0_f64

    do local_index = 1, solv%num_splines_loc
      
       ! TODO ... add function local_to_global
       global_index = (cell_index - 1) * solv%num_splines_loc + local_index
       
       ! Getting coordinates of quadrature points and quadrature weights
       qpt1  = solv%quad_pts1(global_index)
       qpt2  = solv%quad_pts2(global_index)
       wqpt  = solv%quad_weight(global_index)
       
       ! Initialization of the matrices
       val_c        = c_field%value_at_point(qpt1,qpt2)
       val_a11      = a11_field_mat%value_at_point(qpt1,qpt2)
       val_a12      = a12_field_mat%value_at_point(qpt1,qpt2)
       val_a21      = a21_field_mat%value_at_point(qpt1,qpt2)
       val_a22      = a22_field_mat%value_at_point(qpt1,qpt2)
       val_b1       = b1_field_vect%value_at_point(qpt1,qpt2)
       val_b1_der1  = b1_field_vect%first_deriv_eta1_value_at_point(qpt1,qpt2)
       val_b1_der2  = b1_field_vect%first_deriv_eta2_value_at_point(qpt1,qpt2)
       val_b2       = b2_field_vect%value_at_point(qpt1,qpt2)
       val_b2_der1  = b2_field_vect%first_deriv_eta1_value_at_point(qpt1,qpt2)
       val_b2_der2  = b2_field_vect%first_deriv_eta2_value_at_point(qpt1,qpt2)
       
       ! Computing jacobian matrix, jacobian and integral of the jacobian
       jac_mat(:,:) = c_field%get_jacobian_matrix(qpt1,qpt2)
       val_jac = jac_mat(1,1)*jac_mat(2,2) - jac_mat(1,2)*jac_mat(2,1)
       solv%values_jacobian(global_index) = val_jac 
       solv%intjac = solv%intjac + wqpt*val_jac
       
       ! TODO : modify this part so that the jacobian intervenes only after
       ! this part, computed in the physical.
       ! The B matrix is  by (J^(-1)) A^T (J^(-1))^T 
       B11 = jac_mat(2,2)*jac_mat(2,2)*val_a11 - &
            jac_mat(2,2)*jac_mat(1,2)*(val_a12+val_a21) + &
            jac_mat(1,2)*jac_mat(1,2)*val_a22
       
       B21 = jac_mat(1,1)*jac_mat(2,2)*val_a12 - &
            jac_mat(1,1)*jac_mat(1,2)*val_a22 - &
            jac_mat(2,1)*jac_mat(2,2)*val_a11 + &
            jac_mat(1,2)*jac_mat(2,1)*val_a21
       
       B12 = jac_mat(1,1)*jac_mat(2,2)*val_a21 - &
            jac_mat(1,1)*jac_mat(1,2)*val_a22 - &
            jac_mat(2,1)*jac_mat(2,2)*val_a11 + &
            jac_mat(1,2)*jac_mat(2,1)*val_a12
       
       B22 = jac_mat(1,1)*jac_mat(1,1)*val_a22 - &
            jac_mat(1,1)*jac_mat(2,1)*(val_a21+val_a12) + &
            jac_mat(2,1)*jac_mat(2,1)*val_a11
       
       MC =   jac_mat(2,2) * val_b1_der1 &
            - jac_mat(2,1) * val_b1_der2 &
            - jac_mat(1,2) * val_b2_der1 &
            + jac_mat(1,1) * val_b2_der2
          
       C1 =   jac_mat(2,2) * val_b1 &
            - jac_mat(1,2) * val_b2 
       
       C2 =   jac_mat(1,1) * val_b2 &
            - jac_mat(2,1) * val_b1
       
       ! loop over the splines supported in the cell that are different than
       ! zero at the point (qpt1,qpt2) (there are spline_degree+1 splines in
       ! each direction.
       ! TODO : creer une variable pour weigthed volume jac
       do index1 = 1, solv%num_splines_loc
             
          Masse_loc(index1) = &
               Masse_loc(index1) + &
               val_jac * wqpt * &
               solv%values_splines(index1, 1)
          
          Stiff_loc(index1) = &
               Stiff_loc(index1) + &
               val_jac * wqpt * &
               solv%values_splines(index1, 2)* &
               solv%values_splines(index1, 3)
          
          do index2 = 1, solv%num_splines_loc
             
             Source_loc(cell_index,index1, index2) = &
                  Source_loc(cell_index,index1, index2) + &
                  val_jac * wqpt * &
                  solv%values_splines_source(index1, 1)* &
                  solv%values_splines_source(index2, 1)
             
             M_c_loc(index1, index2) = &
                  M_c_loc(index1, index2) + &
                  val_c * val_jac * wqpt * &
                  solv%values_splines(index1, 1)* &
                  solv%values_splines(index2, 1)
             
             K_a11_loc(index1, index2) = &
                  K_a11_loc(index1, index2) + &
                  B11 * wqpt / val_jac * &
                  solv%values_splines(index1, 2)* &
                  solv%values_splines(index2, 2)
                   
             K_a22_loc(index1, index2) = &
                  K_a22_loc(index1, index2) + &
                  B22 * wqpt / val_jac *  &
                  solv%values_splines(index1, 3)* &
                  solv%values_splines(index2, 3)
             
             K_a12_loc(index1, index2) = &
                  K_a12_loc(index1, index2) + &
                  B12 * wqpt / val_jac *  &
                  solv%values_splines(index1, 2)* &
                  solv%values_splines(index2, 3)
             
             K_a21_loc(index1, index2) = &
                  K_a21_loc(index1, index2) +&
                  B21 * wqpt / val_jac *  &
                  solv%values_splines(index1, 3)* &
                  solv%values_splines(index2, 2)
             
             M_b_vect_loc(index1, index2) =      &
                  M_b_vect_loc(index1, index2) + &
                  MC * wqpt *  &
                  solv%values_splines(index1, 1)* &
                  solv%values_splines(index2, 1)
             
             S_b1_loc(index1, index2) =      &
                  S_b1_loc(index1, index2) + &
                  C1 * wqpt *  &
                  solv%values_splines(index1, 1)* &
                  solv%values_splines(index2, 2)
             
             S_b2_loc(index1, index2) = &
                  S_b2_loc(index1, index2) + &
                  C2 * wqpt *  &
                  solv%values_splines(index1, 1)* &
                  solv%values_splines(index2, 3)
          end do
       end do
    end do
  end subroutine build_local_matrices
  
  

  
  subroutine build_local_matrices_source( &
       solv, &
       cell_index, &
       source_at_quad, &
       M_source_loc)

    class(finite_elements_solver) :: solv
    sll_int32, intent(in)         :: cell_index
    sll_real64, dimension(:), intent(in)   :: source_at_quad
    sll_real64, dimension(:), intent(out)  :: M_source_loc
    ! Boundary conditions
    sll_int32 :: bc_left    
    sll_int32 :: bc_right
    sll_int32 :: bc_bottom    
    sll_int32 :: bc_top
    sll_int32  :: tmp1
    sll_int32  :: tmp2
    sll_int32  :: quad_index
    sll_int32  :: ii
    sll_int32  :: jj
    sll_int32  :: index1
    sll_real64 :: val_src
    sll_real64 :: val_jac
    sll_real64 :: spline1,spline2
    
    M_source_loc(:) = 0.0_f64

    do local_index = 1, solv%num_splines_loc
       
       global_index = (cell_index - 1) * solv%num_splines_loc + local_index
       
       ! Getting weights associated to point at global_index
       weight = solv%quad_weight(global_index)
       
       val_src = source_at_quad(global_index)
       val_jac = solv%values_jacobian(quad_index)

       ! loop over the splines supported in the cell that are different than
       ! zero at the point (qpt1,qpt2) (there are spline_degree+1 splines in
       ! each direction.
       do index1 = 1, solv%num_splines_loc

          M_source_loc(index1)= M_source_loc(index1) + &
               val_src * val_jac * weight * &
               solv%values_splines(index1, 1)
          
       end do
    end do
  end subroutine build_local_matrices_source
  


  subroutine local_to_global_matrices( &
       solv,&
       cell_index, &
       Masse_loc,&
       Stiff_loc,&
       M_c_loc, &
       K_a11_loc, &
       K_a12_loc, &
       K_a21_loc, &
       K_a22_loc,&
       M_b_vect_loc, &
       S_b1_loc, &
       S_b2_loc, &
       Masse_tot,&
       Stiff_tot)
    
    class(finite_elements_solver)  :: solv
    sll_int32,                  intent(in) :: cell_index
    sll_real64, dimension(:,:), intent(in) :: M_c_loc
    sll_real64, dimension(:,:), intent(in) :: K_a11_loc
    sll_real64, dimension(:,:), intent(in) :: K_a12_loc
    sll_real64, dimension(:,:), intent(in) :: K_a21_loc
    sll_real64, dimension(:,:), intent(in) :: K_a22_loc
    sll_real64, dimension(:,:), intent(in) :: M_b_vect_loc
    sll_real64, dimension(:,:), intent(in) :: S_b1_loc
    sll_real64, dimension(:,:), intent(in) :: S_b2_loc
    
    !  Correspond to the full Matrix of linear system 
    !  It is not necessary to keep it  
    sll_real64, dimension(:), intent(in) :: Masse_loc
    sll_real64, dimension(:), intent(in) :: Stiff_loc
    sll_real64, dimension(:), intent(inout) :: Masse_tot
    sll_real64, dimension(:), intent(inout) :: Stiff_tot
    sll_real64 :: elt_mat_global    
    sll_int32 :: index1, index2, index3, index4
    sll_int32 :: index_coef1, index_coef2
    sll_int32 :: index
    sll_int32 :: i,j,mm, nn, b, bprime,x,y
    sll_int32 :: li_A, li_Aprime
    sll_int32 :: nbsp,nbsp1
    sll_int32 :: bc_left
    sll_int32 :: bc_right
    sll_int32 :: bc_bottom
    sll_int32 :: bc_top
    
    bc_left   = solv%bc_left
    bc_right  = solv%bc_right
    bc_bottom = solv%bc_bottom
    bc_top    = solv%bc_top
    
    
!     do mm = 0,solv%spline_degree2
!        index3 = cell_j + mm
       
!        if (  (bc_bottom==SLL_PERIODIC).and.(bc_top== SLL_PERIODIC) ) then 
          
!           if ( index3 > solv%total_num_splines_eta2) then
!              index3 = index3 - solv%total_num_splines_eta2
!           end if
          
!        end if
       
!        do i = 0,solv%spline_degree1
          
!           index1 = cell_i + i
!           if ( (bc_left==SLL_PERIODIC).and.(bc_right== SLL_PERIODIC)) then 
!              if ( index1 > solv%total_num_splines_eta1) then
                
!                 index1 = index1 - solv%total_num_splines_eta1
                
!              end if
!              nbsp = solv%total_num_splines_eta1
             
!           else if ( (bc_left  == SLL_DIRICHLET).and.&
!                (bc_right == SLL_DIRICHLET) ) then
!              nbsp = solv%num_cells1 + solv%spline_degree1
!           end if
!           x          =  index1 + (index3-1)*nbsp
!           b          =  mm * ( solv%spline_degree1 + 1 ) + i + 1
!           li_A       =  solv%local_to_global_spline_indices(b, cell_index)
          
!           Masse_tot(x)    = Masse_tot(x) + Masse_loc(b)
!           Stiff_tot(x)    = Stiff_tot(x) + Stiff_loc(b)
          
!           do nn = 0,solv%spline_degree2
             
!              index4 = cell_j + nn
             
!              if ( (bc_bottom==SLL_PERIODIC).and.(bc_top== SLL_PERIODIC))then
!                 if ( index4 > solv%total_num_splines_eta2) then
                   
!                    index4 = index4 - solv%total_num_splines_eta2
!                 end if
!              end if
             
!              do j = 0,solv%spline_degree1
                
!                 index2 = cell_i + j
!                 if((bc_left==SLL_PERIODIC).and.(bc_right==SLL_PERIODIC))then
                   
!                    if ( index2 > solv%total_num_splines_eta1) then
                      
!                       index2 = index2 - solv%total_num_splines_eta1
!                    end if
!                    nbsp1 = solv%total_num_splines_eta1
                   
!                 else if ( (bc_left  == SLL_DIRICHLET) .and.&
!                      (bc_right == SLL_DIRICHLET) ) then
                   
!                    nbsp1 = solv%num_cells1 + solv%spline_degree1
!                 end if
                
!                 y         = index2 + (index4-1)*nbsp1
!                 bprime    =  nn * ( solv%spline_degree1 + 1 ) + j + 1
!                 li_Aprime = solv%local_to_global_spline_indices(bprime,cell_index)
!                 elt_mat_global = &
!                      M_c_loc(b, bprime)     - &
!                      K_a11_loc(b, bprime)   - &
!                      K_a12_loc(b, bprime)   - &
!                      K_a21_loc(b, bprime)   - &
!                      K_a22_loc(b, bprime)   - &
!                      M_b_vect_loc(b,bprime) - &
!                      S_b1_loc( b, bprime)   - &
!                      S_b2_loc( b, bprime)
!                 ! index_coef1 = solv%tab_index_coeff1(cell_i)- solv%spline_degree1 + i
! !                 index_coef2 = solv%tab_index_coeff2(cell_j)- solv%spline_degree2+ mm
!                 index = tab_index_coeff(cell_index)!index_coef1 + (index_coef2-1)*(solv%num_cells1+1)
!                 solv%local_to_global_spline_indices_source(b,cell_index)= index
                
!                 solv%local_to_global_spline_indices_source_bis(bprime,cell_index)= y
                
!                 if ( (li_A > 0) .and. (li_Aprime > 0) ) then
                   
!                    call sll_add_to_csr_matrix( &
!                         solv%sll_csr_mat, &
!                         elt_mat_global, &
!                         li_A, &
!                         li_Aprime)   


!                 end if
                
!              end do
             
!           end do
!        end do
!     end do
    
  end subroutine local_to_global_matrices
  
  
  subroutine local_to_global_matrices_source( &
       solv,&
       cell_i, &
       cell_j, &
       M_source_loc)
    
     class(finite_elements_solver)  :: solv
     sll_int32 :: cell_i
     sll_int32 :: cell_j
     sll_real64, dimension(:), intent(in)   :: M_source_loc
     sll_int32 :: i,mm, b, x!,y
     sll_int32 :: nbsp!,nbsp1
     sll_int32 :: bc_left
     sll_int32 :: bc_right
     sll_int32 :: bc_bottom
     sll_int32 :: bc_top
     sll_int32 :: index1,index3
     
     bc_left   = solv%bc_left
     bc_right  = solv%bc_right
     bc_bottom = solv%bc_bottom
     bc_top    = solv%bc_top
     
     
     
!     do mm = 0,solv%spline_degree2
!        index3 = cell_j + mm
!        if (  (bc_bottom==SLL_PERIODIC).and.(bc_top== SLL_PERIODIC) ) then    
!           if ( index3 > solv%total_num_splines_eta2) then
!              index3 = index3 - solv%total_num_splines_eta2
!           end if
!        end if
! !other option for above:      index3 = mod(index3 - 1, solv%total_num_splines_eta2) + 1
       
!        do i = 0,solv%spline_degree1
          
!           index1 = cell_i + i
!           if ( (bc_left==SLL_PERIODIC).and.(bc_right== SLL_PERIODIC)) then 
!              if ( index1 > solv%total_num_splines_eta1) then
                
!                 index1 = index1 - solv%total_num_splines_eta1
                
!              end if
!              nbsp = solv%total_num_splines_eta1
             
!           else if ( (bc_left  == SLL_DIRICHLET).and.&
!                (bc_right == SLL_DIRICHLET) ) then
!              nbsp = solv%num_cells1 + solv%spline_degree1
!           end if
!           x          =  index1 + (index3-1)*nbsp
!           b          =  mm * ( solv%spline_degree1 + 1 ) + i + 1
!           solv%source_vec(x)  =  solv%source_vec(x)  + M_source_loc(b)
          
!        end do
!     end do
    
  end subroutine local_to_global_matrices_source
  
  subroutine solve_linear_system( solv )
    ! CSR_MAT*phi = source_vec is the linear system to be solved. The solution
    ! is given in terms of the spline coefficients that represent phi.
    class(finite_elements_solver) :: solv
    !type(csr_matrix)  :: csr_masse
    integer :: elt, elt1
    integer :: i,j
     character(len=*),parameter :: as_file='source', as_file1='phi',as_file2='mat'
    sll_int32 :: bc_left
    sll_int32 :: bc_right
    sll_int32 :: bc_bottom
    sll_int32 :: bc_top
    sll_int32 :: ierr

    bc_left   = solv%bc_left
    bc_right  = solv%bc_right
    bc_bottom = solv%bc_bottom
    bc_top    = solv%bc_top
    
    solv%tmp_source_vec(:) = 0.0_f64
    solv%tmp_phi_vec(:) = 0.0_f64
    
!     if( (bc_left   == SLL_PERIODIC) .and. (bc_right == SLL_PERIODIC) .and. &
!          (bc_bottom == SLL_DIRICHLET).and. (bc_top   == SLL_DIRICHLET) ) then
       
!        do i = 1, solv%total_num_splines_eta1
!           do j = 1, solv%total_num_splines_eta2
             
!              elt  = i + solv%total_num_splines_eta1 * (  j - 1)
!              elt1 = i + ( solv%total_num_splines_eta1 ) * j
!              solv%tmp_source_vec(elt) = solv%source_vec(elt1)
!           end do
!        end do
       
!     else if ( (bc_left   == SLL_DIRICHLET).and.(bc_right==SLL_DIRICHLET) .and.&
!          (bc_bottom == SLL_DIRICHLET).and.(bc_top==SLL_DIRICHLET) ) then 
       
!        do i = 1, solv%total_num_splines_eta1
!           do j = 1, solv%total_num_splines_eta2
             
!              elt  = i + solv%total_num_splines_eta1 * (  j - 1)
!              elt1 = i + 1 + ( solv%total_num_splines_eta1 + 2 ) * j 
!              solv%tmp_source_vec( elt ) = solv%source_vec( elt1 )
!           end do
!        end do
       
!     else if((bc_left   == SLL_PERIODIC) .and. (bc_right==SLL_PERIODIC) .and.&
!          (bc_bottom == SLL_PERIODIC) .and. (bc_top  ==SLL_PERIODIC)) then
       
!        solv%tmp_source_vec(1:solv%total_num_splines_eta1*solv%total_num_splines_eta2)=&
!             solv%source_vec(1:solv%total_num_splines_eta1*solv%total_num_splines_eta2) 
       
       
!     else if( (bc_left == SLL_DIRICHLET) .and. (bc_right == SLL_DIRICHLET) .and.&
!              (bc_bottom == SLL_PERIODIC).and. (bc_top   == SLL_PERIODIC) ) then
       
!        do i = 1, solv%total_num_splines_eta1
!           do j = 1, solv%total_num_splines_eta2

!              elt1 = i + 1 + ( solv%total_num_splines_eta1 + 2 ) * (  j - 1)
!              elt  = i + solv%total_num_splines_eta1 * (  j - 1)
!              solv%tmp_source_vec( elt ) = solv%source_vec( elt1 )
!           end do
!        end do

       
!     end if

    call solve_gen_elliptic_eq(solv%sll_csr_mat,solv%tmp_source_vec,solv%tmp_phi_vec)
    solv%phi_vec(1:solv%num_splines_tot)=&
         solv%tmp_phi_vec(1:solv%num_splines_tot)

  end subroutine solve_linear_system
  
  subroutine solve_gen_elliptic_eq(csr_mat,apr_B,apr_U)
    type(sll_csr_matrix) :: csr_mat
    sll_real64, dimension(:) :: apr_U
    sll_real64, dimension(:) :: apr_B 
    sll_int32  :: ai_maxIter
    sll_real64 :: ar_eps
    
    ar_eps = 1.d-13
    ai_maxIter = 100000
  
    call sll_solve_csr_matrix(csr_mat, apr_B, apr_U)
   
  end subroutine solve_gen_elliptic_eq
  
  
  
  subroutine solve_linear_system_perper( solv)
    ! CSR_MAT*phi = source_vec is the linear system to be solved. The solution
    ! is given in terms of the spline coefficients that represent phi.
    class(finite_elements_solver) :: solv
    sll_int32 :: ierr
    
    solv%tmp_source_vec(:) = 0.0_f64
    solv%tmp_phi_vec(:) = 0.0_f64
    solv%tmp_source_vec(1:solv%num_splines_tot)=&
         solv%source_vec(1:solv%num_splines_tot)
    call solve_general_es_perper(solv%sll_csr_mat_with_constraint,solv%tmp_source_vec,solv%tmp_phi_vec) 
    solv%phi_vec(1:solv%num_splines_tot) = &
       solv%tmp_phi_vec(1:solv%num_splines_tot)

  end subroutine solve_linear_system_perper


  subroutine solve_general_solv_perper(csr_mat,apr_B,apr_U)
    type(sll_csr_matrix) :: csr_mat
    sll_real64, dimension(:) :: apr_U
    sll_real64, dimension(:) :: apr_B 
    ! We use a simple conjugate gradient on the new matrice csr_mat_with_constraint (with constraint)
    call sll_solve_csr_matrix(&
         csr_mat,&
         apr_B,&
         apr_U)
    ! We use a modified conjugate gradient on the matrice csr_mat (without constraint)     
    !call sll_solve_csr_matrix_perper(&
    !     csr_mat,&
    !     apr_B,&
    !     apr_U,&
    !     solv%masse)
  end subroutine solve_general_solv_perper

  subroutine compute_Source_matrice(solv,Source_loc)
    type(finite_elements_solver),intent(inout) :: solv
    sll_real64, dimension(:,:,:), pointer :: Source_loc
    sll_int32 :: cell_j,cell_i
    sll_int32 :: cell_index
    sll_int32 :: ideg2,ideg1
    sll_int32 :: jdeg2,jdeg1
    sll_int32 :: b, bprime
    sll_int32 :: li_A,li_Aprime
    sll_real64:: elt_mat_global
    
    do cell_j=1,solv%num_cells2
       do cell_i=1,solv%num_cells1
          
          cell_index = cell_i+solv%num_cells1*(cell_j-1)
          
          do ideg2 = 0,solv%spline_degree2
             
             do ideg1 = 0,solv%spline_degree1
                
                b          =  ideg2 * ( solv%spline_degree1 + 1 ) + ideg1 + 1
                li_A       =  solv%local_to_global_spline_indices_source_bis(b, cell_index)
                
                do jdeg2 = 0,solv%spline_degree2
                   
                   do jdeg1 = 0,solv%spline_degree1
                      
                      bprime    =  jdeg2 * ( solv%spline_degree1 + 1 ) + jdeg1 + 1
                      li_Aprime = solv%local_to_global_spline_indices_source(bprime,cell_index)
                      
                      elt_mat_global = Source_loc(cell_index,bprime,b)

                      if ( (li_A > 0) .and. (li_Aprime > 0)) then

                         call sll_add_to_csr_matrix( &
                           solv%sll_csr_mat_source, &
                           elt_mat_global, &
                           li_A, &
                           li_Aprime)

                      end if
                      
                   end do
                end do
             end do
          end do
       end do
    end do


  end subroutine compute_Source_matrice


  subroutine delete_solver( solv )
   type(finite_elements_solver) :: solv
   sll_int32 :: ierr
   ! it is not good to check some cases and not others, fix...
   if(associated(solv%knots1)) then
      SLL_DEALLOCATE(solv%knots1,ierr)
   else
       print *, 'delete_solver, WARNING: knots1 array was not allocated.'
    end if
    if(associated(solv%knots2)) then
       SLL_DEALLOCATE(solv%knots2,ierr)
    else
       print *, 'delete_solver finite elements, ', &
            'WARNING: knots2 array was not allocated.'
    end if
    SLL_DEALLOCATE(solv%global_spline_indices,ierr)
    SLL_DEALLOCATE(solv%local_spline_indices,ierr)
    SLL_DEALLOCATE(solv%local_to_global_spline_indices,ierr)
    SLL_DEALLOCATE(solv%local_to_global_spline_indices_source,ierr)
    SLL_DEALLOCATE(solv%local_to_global_spline_indices_source_bis,ierr)
    SLL_DEALLOCATE(solv%quad_pts1,ierr)
    SLL_DEALLOCATE(solv%quad_pts2,ierr)
    SLL_DEALLOCATE(solv%quad_weight,ierr)
    SLL_DEALLOCATE(solv%knots1_source,ierr)
    SLL_DEALLOCATE(solv%knots2_source,ierr)
    SLL_DEALLOCATE(solv%source_vec,ierr)
    SLL_DEALLOCATE(solv%phi_vec,ierr)
    SLL_DEALLOCATE(solv%masse,ierr)
    SLL_DEALLOCATE(solv%stiff,ierr)
    SLL_DEALLOCATE(solv%tmp_source_vec,ierr)
    SLL_DEALLOCATE(solv%tmp_phi_vec,ierr)
    call sll_delete(solv%sll_csr_mat)
    call sll_delete(solv% sll_csr_mat_with_constraint)
    call sll_delete(solv%sll_csr_mat_source)
    SLL_DEALLOCATE(solv%values_splines,ierr)
    SLL_DEALLOCATE(solv%values_splines_source,ierr)
    SLL_DEALLOCATE(solv%values_jacobian,ierr)
    SLL_DEALLOCATE(solv%tab_index_coeff,ierr)
  end subroutine delete_solver


end module finite_elements_solver_module
