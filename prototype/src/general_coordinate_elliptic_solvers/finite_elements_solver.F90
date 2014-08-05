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
     sll_real64, dimension(:),   pointer :: knots1_rho
     sll_real64, dimension(:),   pointer :: knots2_rho

     ! Quadrature points/weights (gaussian, fekete, ...) of all the domain, by elements
     sll_real64, dimension(:),   pointer :: num_quad_pts
     sll_real64, dimension(:),   pointer :: quad_pts1
     sll_real64, dimension(:),   pointer :: quad_pts2
     sll_real64, dimension(:),   pointer :: quad_weight

     ! Degree of splines in each direction
     sll_int32 :: spline_degree1
     sll_int32 :: spline_degree2
     ! Number of local splines and on each direction
     sll_int32 :: total_num_splines_loc
     sll_int32 :: total_num_splines_eta1
     sll_int32 :: total_num_splines_eta2
     ! Global indexing of splines. The indexing of the splines
     ! in this array depends on the boundary conditions and it's
     ! independent of the mesh type.
     sll_int32, dimension(:), pointer :: global_spline_indices
     ! Local indexing of splines, independent of the mesh type
     sll_int32, dimension(:), pointer :: local_spline_indices
     ! Same as global_spline_indices but including 
     ! the changes resulting from the boundary conditions.
     sll_int32, dimension(:), pointer :: local_to_global_spline_indices
     sll_int32, dimension(:), pointer :: local_to_global_spline_indices_source
     ! TODO : what is this one for ? source = rho ?
     sll_int32, dimension(:), pointer :: local_to_global_spline_indices_source_bis
     ! The following tables contain the values of splines (with degree <= deg) 
     ! in all quadrature points for the test function and source (rho)
     ! dimensions = (spacial, degree)
     sll_real64, dimension(:,:), pointer :: values_splines 
     sll_real64, dimension(:,:), pointer :: values_splines_source

     ! Table with the values of the jacobian of the transformation
     ! in all quadrature points
     sll_real64, dimension(:), pointer :: values_jacobian 

     ! Table with index of point to the left of the support
     sll_int32 , dimension(:)  , pointer :: tab_index_coeff

     ! Compressed sparse row matrix for test function (with constraint) and source (rho)
     type(sll_csr_matrix), pointer :: sll_csr_mat
     type(sll_csr_matrix), pointer :: sll_csr_mat_with_constraint
     type(sll_csr_matrix), pointer :: sll_csr_mat_source

     ! Initial data, useful variables and matrices (in global indexing)
     sll_real64, dimension(:), pointer :: rho_vec
     sll_real64, dimension(:), pointer :: phi_vec
     sll_real64, dimension(:), pointer :: tmp_rho_vec
     sll_real64, dimension(:), pointer :: tmp_phi_vec
     sll_real64, dimension(:), pointer :: masse
     sll_real64, dimension(:), pointer :: stiff
  
     ! TODO : what are this variables for ??
     sll_real64 :: epsi
     sll_real64 :: intjac

  end type finite_elements_solver

  ! For the integration mode.  
  sll_int32, parameter :: ES_GAUSS_LEGENDRE = 0, ES_GAUSS_LOBATTO = 1
  
  interface delete
     module procedure delete_solver
  end interface delete

  interface initialize
     module procedure initialize_finite_elements_solver
  end interface initialize
  
  
contains ! =============================================================


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
    sll_int32 :: num_splines
    sll_int32 :: vec_sz ! for rho_vec and phi_vec allocations
    sll_int32 :: ierr,ierr1
    sll_int32 :: solution_size
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
    sll_int32 :: size1, size2, size3, size4

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
    SLL_ALLOCATE(solv%quad_pts1   (solv%num_quad_pts),ierr)
    SLL_ALLOCATE(solv%quad_pts2   (solv%num_quad_pts),ierr)
    SLL_ALLOCATE(solv%quad_weight(solv%num_quad_pts),ierr)
    ! Initialization :
    ! TODO (@LM) : Is this necesseray ?!
    solv%quad_pts1(:)    = 0.0_f64
    solv%quad_pts2(:)    = 0.0_f64
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
             end do
          end do
       end do
    end do
    !  ---------------------------------------------- END QUADRATURE POINTS INIT
    ! --------------------------------------------------------------------------


    ! --------------------------------------------------------------------------
    ! ------------------- BEGIN ALLOCATION AND INITIALIZATION OF BASIS FUNCTIONS

    ! Number of local splines :
    solv%total_num_splines_loc = (spline_degree_eta1+1)*(spline_degree_eta2+1)
    
    ! The total number of splines in a single direction is given by
    ! num_cells + spline_degree. For 2D is just a product of both directions
    num_splines1 = mesh%num_cells1 + spline_degree_eta1
    num_splines2 = mesh%num_cells2 + spline_degree_eta2
    SLL_ALLOCATE(solv%global_spline_indices(num_splines1*num_splines2),ierr)
    solv%global_spline_indices(:) = 0
    
    ! Local indexing of the splines 
    SLL_ALLOCATE(solv%local_spline_indices(solv%total_num_splines_loc * solv%num_cells),ierr)
    solv%local_spline_indices(:) = 0
    
    ! Connectivity between local and global indexing systems
    SLL_ALLOCATE(solv%local_to_global_spline_indices(solv%total_num_splines_loc * solv%num_cells), ierr)
    solv%local_to_global_spline_indices = 0
    ! Same for source
    SLL_ALLOCATE(solv%local_to_global_spline_indices_source(solv%total_num_splines_loc * solv%num_cells),ierr)
    solv%local_to_global_spline_indices_source = 0
    ! Same for source (bis)
    SLL_ALLOCATE(solv%local_to_global_spline_indices_source_bis(solv%total_num_splines_loc * solv%num_cells),ierr)
    solv%local_to_global_spline_indices_source_bis = 0

    ! --------------------- END ALLOCATION AND INITIALIZATION OF BASIS FUNCTIONS
    ! --------------------------------------------------------------------------

    ! --------------------------------------------------------------------------
    ! BEGIN ALLOCATION AND INITIALIZATION OF SPLINES KNOTS ---------------------

    if( (bc_left == SLL_PERIODIC) .and. (bc_right == SLL_PERIODIC) .and. &
         (bc_bottom == SLL_PERIODIC) .and. (bc_top == SLL_PERIODIC) ) then
       solv%total_num_splines_eta1 = mesh%num_cells1 
       solv%total_num_splines_eta2 = mesh%num_cells2
       knots1_size = 2 * spline_degree_eta1 + 2
       knots2_size = 2 * spline_degree_eta2 + 2
       vec_sz      = solv%num_cells
       sll_perper  = 1 
    else if( (bc_left == SLL_PERIODIC) .and. (bc_right == SLL_PERIODIC) .and.&
         (bc_bottom == SLL_DIRICHLET) .and. (bc_top == SLL_DIRICHLET) ) then
       solv%total_num_splines_eta1 = mesh%num_cells1 
       solv%total_num_splines_eta2 = mesh%num_cells2 + spline_degree_eta2 - 2
       knots1_size = 2*spline_degree_eta1 + 2
       knots2_size = 2*spline_degree_eta2 + mesh%num_cells2+ 1
       vec_sz      = mesh%num_cells1 * (mesh%num_cells2 + spline_degree_eta2)
    else if( (bc_left == SLL_DIRICHLET) .and. (bc_right == SLL_DIRICHLET) .and.&
         (bc_bottom == SLL_PERIODIC) .and. (bc_top == SLL_PERIODIC) ) then
       solv%total_num_splines_eta1 = mesh%num_cells1 + spline_degree_eta1 - 2
       solv%total_num_splines_eta2 = mesh%num_cells2 
       knots1_size = 2*spline_degree_eta1 + mesh%num_cells1+1
       knots2_size = 2*spline_degree_eta2 + 2
       vec_sz      = (mesh%num_cells1 + spline_degree_eta1) * mesh%num_cells2
    else if( (bc_left == SLL_DIRICHLET) .and. (bc_right == SLL_DIRICHLET) .and.&
         (bc_bottom == SLL_DIRICHLET) .and. (bc_top == SLL_DIRICHLET) ) then
       solv%total_num_splines_eta1 = mesh%num_cells1 + spline_degree_eta1 - 2
       solv%total_num_splines_eta2 = mesh%num_cells2 + spline_degree_eta2 - 2
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
    SLL_ALLOCATE(solv%knots1_rho(mesh%num_cells1 + spline_degree_eta1 + 2),ierr1)
    SLL_ALLOCATE(solv%knots2_rho(mesh%num_cells2 + spline_degree_eta2 + 2 ),ierr)

    ! ----------------------- END ALLOCATION AND INITIALIZATION OF SPLINES KNOTS
    ! --------------------------------------------------------------------------


    solution_size = solv%total_num_splines_eta1 * solv%total_num_splines_eta2
    
    SLL_ALLOCATE(solv%rho_vec(vec_sz),ierr)
    SLL_ALLOCATE(solv%masse  (vec_sz),ierr)
    SLL_ALLOCATE(solv%stiff  (vec_sz),ierr)
    SLL_ALLOCATE(solv%phi_vec(solution_size),ierr)
    
    if(sll_perper == 1) then
       SLL_ALLOCATE(solv%tmp_rho_vec(solution_size + 1),ierr)
       SLL_ALLOCATE(solv%tmp_phi_vec(solution_size + 1),ierr)
    else
       SLL_ALLOCATE(solv%tmp_rho_vec(solution_size),ierr)
       SLL_ALLOCATE(solv%tmp_phi_vec(solution_size),ierr)
    endif

    solv%rho_vec(:) = 0.0_f64
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
         solution_size, &
         solution_size, &
         solv%num_cells, &
         solv%local_to_global_spline_indices, &
         solv%total_num_splines_loc, &
         solv%local_to_global_spline_indices, &
         solv%total_num_splines_loc)
    
    solv%knots1_rho(1:spline_degree_eta1+1) = mesh%eta1_min
    solv%knots1_rho(mesh%num_cells1+2:mesh%num_cells1+1+spline_degree_eta1+1) = &
         mesh%eta1_max
    
    if (mod(spline_degree_eta1 + 1, 2) == 0) then
       do i = spline_degree_eta1 +2, mesh%num_cells1 + 1
          solv%knots1_rho(i) = eta1_min +  &
               (i - (spline_degree_eta1 +1)/2 - 1) * solv%mesh%delta_eta1  
       end do
    else
       do i = spline_degree_eta1 + 2, mesh%num_cells1 + 1
          solv%knots1_rho ( i ) = &
               0.5_f64*(eta1_min + (i - (spline_degree_eta1)/2 - 1) * solv%mesh%delta_eta1 + &
               eta1_min +  ( i -1 - (spline_degree_eta1)/2 -1) * solv%mesh%delta_eta1 )
       end do
    end if
    
    solv%knots2_rho(1:spline_degree_eta2+1) = eta2_min
    solv%knots2_rho(mesh%num_cells2+2:mesh%num_cells2+1+spline_degree_eta2+1) = eta2_max
    
    if ( mod(spline_degree_eta2 +1,2) == 0 ) then
       do i = spline_degree_eta2 +1 + 1, mesh%num_cells2 + 1
          solv%knots2_rho( i ) = eta2_min + &
               ( i - (spline_degree_eta2 +1)/2-1 )*solv%mesh%delta_eta2 
       end do
    else     
       do i = spline_degree_eta2 +1 + 1, mesh%num_cells2 + 1
          solv%knots2_rho ( i ) = &
               0.5*( eta2_min + ( i - (spline_degree_eta2)/2 -1)*solv%mesh%delta_eta2 + &
               eta2_min +  ( i -1 - (spline_degree_eta2)/2 -1)*solv%mesh%delta_eta2 )  
       end do
    end if
    
    
    ! ------------------------------------------------------------------------------
    ! BEGIN Allocation of basis splines related tables -----------------------------
    ! values_splines : table containning all values --------------------------------
    ! of basis splines in each direction in each gauss points ----------------------

    size1 = mesh%num_cells1*(spline_degree_eta1+2)
    size2 = mesh%num_cells2*(spline_degree_eta2+2)
    size3 = spline_degree_eta1+1
    size4 = spline_degree_eta2+1
    
    SLL_ALLOCATE(solv%values_splines(size1*size2, size3*size4), ierr)
    solv%values_splines = 0.0_f64
    SLL_ALLOCATE(solv%values_splines_source(size1*size2, size3*size4), ierr)
    solv%values_splines_source = 0.0_f64
    SLL_ALLOCATE(solv%values_jacobian(size1*size2), ierr)
    solv%values_jacobian = 0.0_f64
    SLL_ALLOCATE(solv%tab_index_coeff(size1*size2),ierr)

    ! ------------------------------------------------ END ALLOCATION SPLINES TABLES
    ! ------------------------------------------------------------------------------
 
  end subroutine initialize_finite_elements_solver
  

  function new_finite_elements_solver( &
   spline_degree_eta1, &
   spline_degree_eta2, &
   num_cells_eta1, &
   num_cells_eta2, &
   quadrature_type1, &
   quadrature_type2, &
   bc_left, &
   bc_right, &
   bc_bottom, &
   bc_top, &
   eta1_min, &
   eta1_max, &
   eta2_min, &
   eta2_max ) result(es)

   type(finite_elements_solver), pointer :: es
   sll_int32, intent(in) :: spline_degree_eta1
   sll_int32, intent(in) :: spline_degree_eta2
   sll_int32, intent(in) :: num_cells_eta1
   sll_int32, intent(in) :: num_cells_eta2
   sll_int32, intent(in) :: bc_left
   sll_int32, intent(in) :: bc_right
   sll_int32, intent(in) :: bc_bottom
   sll_int32, intent(in) :: bc_top
   sll_int32, intent(in) :: quadrature_type1
   sll_int32, intent(in) :: quadrature_type2
   sll_real64, intent(in) :: eta1_min
   sll_real64, intent(in) :: eta1_max
   sll_real64, intent(in) :: eta2_min
   sll_real64, intent(in) :: eta2_max
   sll_int32 :: ierr


   !call sll_set_time_mark(t0)
   SLL_ALLOCATE(es,ierr)
   call initialize( &
        es, &
        spline_degree_eta1, &
        spline_degree_eta2, &
        num_cells_eta1, &
        num_cells_eta2, &
        quadrature_type1, &
        quadrature_type2, &
        bc_left, &
        bc_right, &
        bc_bottom, &
        bc_top, &
        eta1_min, &
        eta1_max, &
        eta2_min, &
        eta2_max )
   
    !time = sll_time_elapsed_since(t0)
    !print*, '#time for new_finite_elements_solver', time
    

  end function new_finite_elements_solver

  


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
    SLL_DEALLOCATE(solv%knots1_rho,ierr)
    SLL_DEALLOCATE(solv%knots2_rho,ierr)
    SLL_DEALLOCATE(solv%rho_vec,ierr)
    SLL_DEALLOCATE(solv%phi_vec,ierr)
    SLL_DEALLOCATE(solv%masse,ierr)
    SLL_DEALLOCATE(solv%stiff,ierr)
    SLL_DEALLOCATE(solv%tmp_rho_vec,ierr)
    SLL_DEALLOCATE(solv%tmp_phi_vec,ierr)
    call sll_delete(solv%sll_csr_mat)
    call sll_delete(solv% sll_csr_mat_with_constraint)
    call sll_delete(solv%sll_csr_mat_source)
    SLL_DEALLOCATE(solv%values_splines,ierr)
    SLL_DEALLOCATE(solv%values_splines_source,ierr)
    SLL_DEALLOCATE(solv%values_jacobian,ierr)
    SLL_DEALLOCATE(solv%tab_index_coeff,ierr)
  end subroutine delete_solver
  
  subroutine factorize_mat_solv(&
       solv, &
       a11_field_mat, &
       a12_field_mat,&
       a21_field_mat,&
       a22_field_mat,&
       b1_field_vect,&
       b2_field_vect,&
       c_field)!, &
    ! rho)
    use sll_timer
    type(finite_elements_solver),intent(inout) :: solv
    class(sll_scalar_field_2d_base), pointer :: a11_field_mat
    class(sll_scalar_field_2d_base), pointer :: a12_field_mat
    class(sll_scalar_field_2d_base), pointer :: a21_field_mat
    class(sll_scalar_field_2d_base), pointer :: a22_field_mat
    class(sll_scalar_field_2d_base), pointer :: b1_field_vect
    class(sll_scalar_field_2d_base), pointer :: b2_field_vect
    class(sll_scalar_field_2d_base), pointer     :: c_field
    sll_real64, dimension(:,:), allocatable :: M_c_loc
    sll_real64, dimension(:,:), allocatable :: K_a11_loc
    sll_real64, dimension(:,:), allocatable :: K_a12_loc
    sll_real64, dimension(:,:), allocatable :: K_a21_loc
    sll_real64, dimension(:,:), allocatable :: K_a22_loc
    sll_real64, dimension(:,:), allocatable :: M_b_vect_loc
    sll_real64, dimension(:,:), allocatable :: S_b1_loc
    sll_real64, dimension(:,:), allocatable :: S_b2_loc  
    sll_real64, dimension(:), allocatable :: Masse_loc
    sll_real64, dimension(:), allocatable :: Stiff_loc
    sll_real64, dimension(:,:,:), pointer :: Source_loc
    sll_int32 :: total_num_splines_loc
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
    total_num_splines_loc = solv%total_num_splines_loc
    if( (bc_left == SLL_PERIODIC) .and. (bc_right == SLL_PERIODIC) .and. &
       (bc_bottom == SLL_PERIODIC) .and. (bc_top == SLL_PERIODIC) ) then
       sll_perper = 0
    else
       sll_perper = 1  
    end if   
    SLL_ALLOCATE(Source_loc(solv%num_cells,total_num_splines_loc,total_num_splines_loc),ierr)
    SLL_ALLOCATE(M_c_loc(total_num_splines_loc,total_num_splines_loc),ierr)
    SLL_ALLOCATE(K_a11_loc(total_num_splines_loc,total_num_splines_loc),ierr)
    SLL_ALLOCATE(K_a12_loc(total_num_splines_loc,total_num_splines_loc),ierr)
    SLL_ALLOCATE(K_a21_loc(total_num_splines_loc,total_num_splines_loc),ierr)
    SLL_ALLOCATE(K_a22_loc(total_num_splines_loc,total_num_splines_loc),ierr)
    SLL_ALLOCATE(S_b1_loc(total_num_splines_loc,total_num_splines_loc),ierr)
    SLL_ALLOCATE(S_b2_loc(total_num_splines_loc,total_num_splines_loc),ierr)
    SLL_ALLOCATE(M_b_vect_loc(total_num_splines_loc,total_num_splines_loc),ierr)
    SLL_ALLOCATE(Masse_loc(total_num_splines_loc),ierr)
    SLL_ALLOCATE(Stiff_loc(total_num_splines_loc),ierr)

    Masse_loc(:) = 0.0_f64
    Stiff_loc(:) = 0.0_f64
    Source_loc(:,:,:) = 0.0_f64

    do cell_index=1,solv%num_cells
          call build_local_matrices( &
               solv, &
               cell_index,&
            !    i, &
!                j, &
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
         solv%total_num_splines_loc, &
         solv%local_to_global_spline_indices_source, &
         solv%total_num_splines_loc )

    
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
       rho,&
       phi)
    use sll_timer
    class(finite_elements_solver) :: solv
    class(sll_scalar_field_2d_discrete_alt), intent(inout)  :: phi
    class(sll_scalar_field_2d_base), intent(in),target  :: rho
    sll_int32 :: i
    sll_int32 :: j,ierr
    sll_int32 :: cell_index
    sll_int32 :: total_num_splines_loc
    sll_real64 :: int_rho,int_jac
    sll_real64, dimension(:), allocatable   :: M_rho_loc
    sll_real64, dimension(:), allocatable   :: rho_at_quad
    sll_int32 :: num_pts, ig1, ig2, ig, jg
    sll_real64 :: wgpt1, wgpt2, gpt1, gpt2, eta1, eta2
    sll_real64 :: val_jac
    class(sll_scalar_field_2d_base),pointer  :: base_field_pointer
    class(sll_interpolator_2d_base),pointer  :: base_interpolator_pointer
    sll_real64, dimension(:,:), pointer :: coeff_rho
    sll_real64, dimension(:), pointer :: rho_coeff_1d
    
    total_num_splines_loc = solv%total_num_splines_loc
    SLL_ALLOCATE(M_rho_loc(total_num_splines_loc),ierr)
    

    num_quad_pts = solv%num_quad_pts
    SLL_ALLOCATE(rho_at_quad(num_quad_pts),ierr)
    ! TODO (@LM) : Is this necesseray ?!
    rho_at_quad(:) = 0.0_f64
    SLL_ALLOCATE(rho_coeff_1d(solv%num_cells_plus1),ierr)
    
    M_rho_loc     = 0.0_f64
    solv%rho_vec(:) = 0.0_f64
    rho_coeff_1d  = 0.0_f64
   
  
    !call sll_set_time_mark(t0)
    !ES Compute rho at all quad points
    ig1 = 0 
    ig2 = 0 
    int_rho = 0.0_f64
    int_jac = 0.0_f64
    
    ! afin d'optimiser la construction de la matrice rho
    ! on va proceder de la facon qui suit 
    base_field_pointer => rho
    select type( type_field => base_field_pointer)
    class is (sll_scalar_field_2d_discrete_alt)
       base_interpolator_pointer => type_field%interp_2d
       select type( type_interpolator => base_interpolator_pointer)
       class is (arb_deg_2d_interpolator)
          coeff_rho => type_interpolator%get_coefficients()
          ! TODO @LM : this whole part should be re written, quick fix though:
          ! put the spline coefficients in a 1d array
          do j=1,size(coeff_rho, 2)
             do i=1,size(coeff_rho, 1)
                rho_coeff_1d(i+size(coeff_rho,2)*(j-1)) = coeff_rho(i,j)
             end do
          end do

          call sll_mult_csr_matrix_vector(&
               solv%sll_csr_mat_source,&
               rho_coeff_1d,solv%rho_vec)

          if( ((solv%bc_bottom==SLL_PERIODIC).and.(solv%bc_top==SLL_PERIODIC)) &
               .and.((solv%bc_left==SLL_PERIODIC).and.(solv%bc_right==SLL_PERIODIC)) )then
             
             solv%rho_vec = solv%rho_vec - sum(solv%rho_vec)/solv%intjac*solv%masse
          end if
          
       class DEFAULT
          
       do quad_index=1,solv%num_quad_pts ! Loop over quadrature points
          qpt1   = solv%quad_pts1  (quad_index)
          qpt2   = solv%quad_pts2  (quad_index)
          weight = solv%quad_weight(quad_index)
          rho_at_quad(quad_index) = rho%value_at_point(qpt1,qpt2)
          val_jac = solv%values_jacobian(quad_index)
          int_rho = int_rho + rho%value_at_point(qpt1,qpt2)*weight*val_jac 
          int_jac = int_jac + weight*val_jac
       end do

       if( ((solv%bc_bottom==SLL_PERIODIC).and.(solv%bc_top==SLL_PERIODIC)) &
            .and. ((solv%bc_left==SLL_PERIODIC).and.(solv%bc_right==SLL_PERIODIC)) )then   
          rho_at_quad = rho_at_quad - int_rho/int_jac
       end if
          
       do cell_index=1,solv%num_cells ! Loop over elements
          call build_local_matrices_rho( &
               solv, &
               cell_index, &
               rho_at_quad, &
               M_rho_loc)
                
          call local_to_global_matrices_rho( &
               solv, &
               cell_index, &
               M_rho_loc)
       end do
    end select
       
    class is (sll_scalar_field_2d_analytic_alt)
       
    do quad_index=1,solv%num_quad_pts ! Loop over quadrature points
       qpt1   = solv%quad_pts1  (quad_index)
       qpt2   = solv%quad_pts2  (quad_index)
       weight = solv%quad_weight(quad_index)
       rho_at_quad(quad_index) = rho%value_at_point(qpt1,qpt2)
       val_jac = solv%values_jacobian(quad_index)
       int_rho = int_rho + rho%value_at_point(qpt1,qpt2) * weight * val_jac 
       int_jac = int_jac + weight * val_jac
    end do
    if( ((solv%bc_bottom==SLL_PERIODIC).and.(solv%bc_top==SLL_PERIODIC)) &
         .and. ((solv%bc_left==SLL_PERIODIC).and.(solv%bc_right==SLL_PERIODIC)) )then
       rho_at_quad = rho_at_quad - int_rho/int_jac
    end if
          
    do cell_index = 1, solv%num_cells
       call build_local_matrices_rho( &
            solv, &
            cell_index, &
            rho_at_quad, &
            M_rho_loc)
             
       call local_to_global_matrices_rho( &
            solv, &
            cell_index, &
            M_rho_loc)
    end do
    end select
    !time = sll_time_elapsed_since(t0)

    call solve_linear_system(solv)
    
    call  phi%interp_2d%set_coefficients( solv%phi_vec)

    SLL_DEALLOCATE_ARRAY(M_rho_loc,ierr)
    SLL_DEALLOCATE_ARRAY(rho_at_quad,ierr)
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
    sll_int32, intent(in) :: cell_index

    class(sll_scalar_field_2d_base), pointer :: a11_field_mat
    class(sll_scalar_field_2d_base), pointer :: a12_field_mat
    class(sll_scalar_field_2d_base), pointer :: a21_field_mat
    class(sll_scalar_field_2d_base), pointer :: a22_field_mat
    class(sll_scalar_field_2d_base), pointer :: b1_field_vect
    class(sll_scalar_field_2d_base), pointer :: b2_field_vect
    class(sll_scalar_field_2d_base), pointer :: c_field
    sll_real64, dimension(:,:,:), intent(inout) :: Source_loc
    sll_real64, dimension(:,:), intent(out) :: M_c_loc
    sll_real64, dimension(:,:), intent(out) :: K_a11_loc
    sll_real64, dimension(:,:), intent(out) :: K_a12_loc
    sll_real64, dimension(:,:), intent(out) :: K_a21_loc
    sll_real64, dimension(:,:), intent(out) :: K_a22_loc
    sll_real64, dimension(:,:), intent(out) :: M_b_vect_loc
    sll_real64, dimension(:,:), intent(out) :: S_b1_loc
    sll_real64, dimension(:,:), intent(out) :: S_b2_loc
    sll_real64, dimension(:), intent(out) :: Masse_loc
    sll_real64, dimension(:), intent(out) :: Stiff_loc

    sll_int32 :: bc_left    
    sll_int32 :: bc_right
    sll_int32 :: bc_bottom    
    sll_int32 :: bc_top  
    sll_real64 :: delta1
    sll_real64 :: delta2
    sll_real64 :: eta1_min
    sll_real64 :: eta2_min
    sll_int32  :: tmp1
    sll_int32  :: tmp2
    sll_int32  :: num_pts_quad ! total number of quadrature points
    sll_int32  :: i_ele
    sll_int32  :: i,ii,iii
    sll_int32  :: j,jj,jjj
    ! Quadrature points coordinates and associated weights :
    sll_real64 :: qpt1
    sll_real64 :: qpt2
    sll_real64 :: wgpt
    ! Temporary quadrature point for BC treatment :
    sll_real64 :: qtmp1
    sll_real64 :: qtmp2
    sll_int32  :: local_spline_index1
    sll_int32  :: local_spline_index2
    sll_int32  :: index1
    sll_int32  :: index2
    sll_real64, dimension(solv%spline_degree1+1,solv%spline_degree1+1) :: work1
    sll_real64, dimension(solv%spline_degree2+1,solv%spline_degree2+1) :: work2
    sll_real64, dimension(solv%spline_degree1+1,2) :: dbiatx1
    sll_real64, dimension(solv%spline_degree2+1,2) :: dbiatx2
    sll_real64, dimension(solv%spline_degree1+1,2) :: dbiatx1_rho
    sll_real64, dimension(solv%spline_degree2+1,2) :: dbiatx2_rho
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
    sll_real64, dimension(2,2) :: jac_mat
    sll_real64 :: val_jac
    sll_real64 :: B11
    sll_real64 :: B12
    sll_real64 :: B21
    sll_real64 :: B22
    sll_real64 :: MC
    sll_real64 :: C1
    sll_real64 :: C2    
    sll_int32 :: left_x,left_y
    sll_int32 :: mflag_x, mflag_y
    
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
    dbiatx1(:,:)      = 0.0_f64
    dbiatx2(:,:)      = 0.0_f64
    work1(:,:)        = 0.0_f64
    work2(:,:)        = 0.0_f64
    ! The supposition is that all fields use the same logical mesh
    delta1    = solv%mesh%delta_eta1 
    delta2    = solv%mesh%delta_eta2 
    eta1_min  = solv%mesh%eta1_min  
    eta2_min  = solv%mesh%eta2_min  
    tmp1      = (solv%spline_degree1 + 1)/2
    tmp2      = (solv%spline_degree2 + 1)/2
    bc_left   = solv%bc_left
    bc_right  = solv%bc_right
    bc_bottom = solv%bc_bottom
    bc_top    = solv%bc_top
    num_pts_quad = size(solv%quad_pts1,2)
      
    do i_ele = 1,num_pts_quad
       qpt1  = solv%quad_pts1(i_ele)
       qpt2  = solv%quad_pts2(i_ele)
       wgpt2 = solv%quad_weight(i_ele)
       
       ! ----- Boundary condition evaluation -----
       if ((solv%bc_bottom==SLL_PERIODIC).and.(solv%bc_top==SLL_PERIODIC))then
          ! rescale gauss point in interval [0,delta2]
          qtmp2 = 0.5_f64*delta2*( solv%quad_pts2(i_ele) + 1.0_f64) !ATTENTION 0.5
          local_spline_index2 = solv%spline_degree2 + 1
       else if ((solv%bc_bottom == SLL_DIRICHLET).and.&
            (solv%bc_top    == SLL_DIRICHLET)) then
          qtmp2 = qpt2
          local_spline_index2 = solv%spline_degree2 + cell_j
       end if
       ! ----- end of : boundary condition evaluation -----

       call bsplvd( &
            solv%knots2, &
            solv%spline_degree2+1,&
            gtmp2,&
            local_spline_index2,&
            work2,&
            dbiatx2,&
            2)

       !! TODO @LM : this shouldnt be here... in init probably
       call interv(&
            solv%knots2_rho,&
            solv%mesh%num_cells2 + solv%spline_degree2+ 2,&
            gpt2, left_y, mflag_y )

       call bsplvd( &
            solv%knots2_rho, &
            solv%spline_degree2+1,&
            gpt2,&
            left_y,&
            work2,&
            dbiatx2_rho,&
            2)
      
       ! we stocke the values of spline to construct the source term
       solv%values_splines_eta2(cell_j + solv%mesh%num_cells2*(j-1),:) = dbiatx2(:,1)
       solv%values_splines_quad2(cell_j+solv%mesh%num_cells2*(j-1),:)=dbiatx2_rho(:,1)
       solv%tab_index_coeff2(cell_j + solv%mesh%num_cells2*(j-1)) = left_y

       do quad_index=1,solv%num_quad_pts ! Loop over quadrature points
          qpt1   = solv%quad_pts1  (quad_index)
          qpt2   = solv%quad_pts2  (quad_index)
          weight = solv%quad_weight(quad_index)
           

          call bsplvd(&
               solv%knots1,&
               solv%spline_degree1+1,&
               qpt1,&
               local_spline_index1,&
               work1,&
               dbiatx1,&
               2 )

          call interv ( solv%knots1_rho, solv%num_cells1 + solv%spline_degree1+ 2, gpt1, &
               left_x, mflag_x )
     
          !! fonction : pour evaluer les splines en chaque direction
          call bsplvd(&
               solv%knots1_rho,&
               solv%spline_degree1+1,&
               gpt1,&
               left_x,&
               work1,&
               dbiatx1_rho,&
               2 )

          solv%values_splines_eta1(cell_i + solv%num_cells1*(i-1),:) = dbiatx1(:,1)
          solv%values_splines_quad1(cell_i + solv%num_cells1*(i-1),:) = dbiatx1_rho(:,1)
          solv%tab_index_coeff1(cell_i + solv%num_cells1*(i-1)) = left_x

          val_c        = c_field%value_at_point(gpt1,gpt2)
          val_a11      = a11_field_mat%value_at_point(gpt1,gpt2)
          val_a12      = a12_field_mat%value_at_point(gpt1,gpt2)
          val_a21      = a21_field_mat%value_at_point(gpt1,gpt2)
          val_a22      = a22_field_mat%value_at_point(gpt1,gpt2)

          val_b1       = b1_field_vect%value_at_point(gpt1,gpt2)
          val_b1_der1  = b1_field_vect%first_deriv_eta1_value_at_point(gpt1,gpt2)
          val_b1_der2  = b1_field_vect%first_deriv_eta2_value_at_point(gpt1,gpt2)

          val_b2       = b2_field_vect%value_at_point(gpt1,gpt2)
          val_b2_der1  = b2_field_vect%first_deriv_eta1_value_at_point(gpt1,gpt2)
          val_b2_der2  = b2_field_vect%first_deriv_eta2_value_at_point(gpt1,gpt2)
   
          jac_mat(:,:) = c_field%get_jacobian_matrix(gpt1,gpt2)
          val_jac = jac_mat(1,1)*jac_mat(2,2) - jac_mat(1,2)*jac_mat(2,1)!abs(jac_mat(1,1)*jac_mat(2,2) - jac_mat(1,2)*jac_mat(2,1))

          solv%values_jacobian(cell_i + solv%num_cells1*(i-1),cell_j + solv%num_cells2*(j-1)) = val_jac
        
          solv%intjac = solv%intjac + wgpt2*wgpt1*val_jac
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
          ! zero at the point (gpt1,gpt2) (there are spline_degree+1 splines in
          ! each direction.
          do ii = 0,solv%spline_degree1
             do jj = 0,solv%spline_degree2
                
                index1  =  jj * ( solv%spline_degree1 + 1 ) + ii + 1
                
                
                
                Masse_loc(index1) = &
                     Masse_loc(index1) + &
                     val_jac*wgpt1*wgpt2* &
                     (dbiatx1(ii+1,1)*dbiatx2(jj+1,1))

                Stiff_loc(index1) = &
                     Stiff_loc(index1) + &
                     val_jac*wgpt1*wgpt2* &
                     (dbiatx1(ii+1,2)*dbiatx2(jj+1,1)+&
                     dbiatx1(ii+1,1)*dbiatx2(jj+1,2))
                
                
         
                
                do iii = 0,solv%spline_degree1
                   do jjj = 0,solv%spline_degree2
                      
                      index2 =  jjj*(solv%spline_degree1 + 1) + iii + 1
                
                      Source_loc(cell_index,index1, index2) = &
                           Source_loc(cell_index,index1, index2) + &
                           val_jac*wgpt1*wgpt2* &
                           dbiatx1_rho(ii+1,1)*dbiatx1(iii+1,1)*  &
                           dbiatx2_rho(jj+1,1)*dbiatx2(jjj+1,1)
                      
                      
                      
                      M_c_loc(index1, index2) = &
                           M_c_loc(index1, index2) + &
                           val_c*val_jac*wgpt1*wgpt2* &
                           dbiatx1(ii+1,1)*dbiatx1(iii+1,1)*  &
                           dbiatx2(jj+1,1)*dbiatx2(jjj+1,1)
                      
                      
                      K_a11_loc(index1, index2) = &
                           K_a11_loc(index1, index2) + &
                           B11*wgpt1*wgpt2/val_jac* &
                           dbiatx1(ii+1,2)*dbiatx1(iii+1,2)* &
                           dbiatx2(jj+1,1)*dbiatx2(jjj+1,1)
                      
                      K_a22_loc(index1, index2) = &
                           K_a22_loc(index1, index2) + &
                           B22*wgpt1*wgpt2/val_jac*  &
                           dbiatx1(ii+1,1)*dbiatx1(iii+1,1)*   &
                           dbiatx2(jj+1,2)*dbiatx2(jjj+1,2)
                      
                      K_a12_loc(index1, index2) = &
                           K_a12_loc(index1, index2) + &
                           B12*wgpt1*wgpt2/val_jac*  &
                           dbiatx1(ii+1,2)*dbiatx1(iii+1,1) *&
                           dbiatx2(jj+1,1)*dbiatx2(jjj+1,2)
                      
                      K_a21_loc(index1, index2) = &
                           K_a21_loc(index1, index2) +&
                           B21*wgpt1*wgpt2/val_jac*  &
                           dbiatx1(ii+1,1)*dbiatx1(iii+1,2)*   &
                           dbiatx2(jj+1,2)*dbiatx2(jjj+1,1)


                      M_b_vect_loc(index1, index2) =      &
                           M_b_vect_loc(index1, index2) + &
                           MC*wgpt1*wgpt2 *  &
                           dbiatx1(ii+1,1)*dbiatx1(iii+1,1)*   &
                           dbiatx2(jj+1,1)*dbiatx2(jjj+1,1)
                      
                      
                      ! A revoir 
                      S_b1_loc(index1, index2) =      &
                           S_b1_loc(index1, index2) + &
                           C1*wgpt1*wgpt2 *  &
                           dbiatx1(ii+1,1)*dbiatx1(iii+1,2)*   &
                           dbiatx2(jj+1,1)*dbiatx2(jjj+1,1)

                      ! A revoir 
                      S_b2_loc(index1, index2) = &
                           S_b2_loc(index1, index2) + &
                           C2*wgpt1*wgpt2 *  &
                           dbiatx1(ii+1,1)*dbiatx1(iii+1,1)*   &
                           dbiatx2(jj+1,1)*dbiatx2(jjj+1,2)
                   end do
                end do
             end do
          end do
       end do
    end do
    
  end subroutine build_local_matrices
  
  

  
  subroutine build_local_matrices_rho( &
       obj, &
       cell_i, &
       cell_j, &
       rho_at_quad, &
       M_rho_loc)

    class(finite_elements_solver) :: obj
    sll_int32, intent(in) :: cell_i
    sll_int32, intent(in) :: cell_j
    sll_real64, dimension(:,:), intent(in)   :: rho_at_quad

    sll_real64, dimension(:), intent(out)   :: M_rho_loc
    sll_int32 :: bc_left    
    sll_int32 :: bc_right
    sll_int32 :: bc_bottom    
    sll_int32 :: bc_top  
    sll_real64 :: delta1
    sll_real64 :: delta2
    sll_real64 :: eta1_min
    sll_real64 :: eta2_min
    sll_real64 :: eta1
    sll_real64 :: eta2
    sll_int32  :: tmp1
    sll_int32  :: tmp2
    sll_int32  :: num_pts_q1 ! number of quad points in first direction 
    sll_int32  :: num_pts_q2 ! number of quad points in second direction
    sll_int32  :: quad_index
    sll_int32  :: ii
    sll_int32  :: jj
    sll_real64 :: gpt1
    sll_real64 :: gpt2
    sll_real64 :: wgpt1
    sll_real64 :: wgpt2
    sll_int32  :: index1
    sll_real64, dimension(obj%spline_degree1+1,obj%spline_degree1+1) :: work1
    sll_real64, dimension(obj%spline_degree2+1,obj%spline_degree2+1) :: work2
    sll_real64, dimension(obj%spline_degree1+1,2) :: dbiatx1
    sll_real64, dimension(obj%spline_degree2+1,2) :: dbiatx2
    sll_real64 :: val_f
    sll_real64 :: val_jac,spline1,spline2
    
    M_rho_loc(:)  = 0.0_f64
    dbiatx1(:,:)  = 0.0_f64
    dbiatx2(:,:)  = 0.0_f64
    work1(:,:)    = 0.0_f64
    work2(:,:)    = 0.0_f64
    ! The supposition is that all fields use the same logical mesh
    delta1    = obj%delta_eta1
    delta2    = obj%delta_eta2
    eta1_min  = obj%eta1_min
    eta2_min  = obj%eta2_min
    tmp1      = (obj%spline_degree1 + 1)/2
    tmp2      = (obj%spline_degree2 + 1)/2
    bc_left   = obj%bc_left
    bc_right  = obj%bc_right
    bc_bottom = obj%bc_bottom
    bc_top    = obj%bc_top
    num_pts_q1 = size(obj%quad_pts1,2)
    num_pts_q2 = size(obj%quad_pts2,2)
    
    
    eta1  = eta1_min + (cell_i-1)*delta1
    eta2  = eta2_min + (cell_j-1)*delta2
    
   
    do quad_index=1,solv%num_quad_pts ! Loop over quadrature points
       qpt1   = solv%quad_pts1  (quad_index)
       qpt2   = solv%quad_pts2  (quad_index)
       weight = solv%quad_weight(quad_index)
       
       val_f   = rho_at_quad(quad_index)
       val_jac = obj%values_jacobian(quad_index)

       ! loop over the splines supported in the cell that are different than
       ! zero at the point (qpt1,qpt2) (there are spline_degree+1 splines in
       ! each direction.
       do ii = 0,obj%spline_degree1
          do jj = 0,obj%spline_degree2

             spline = obj%values_splines(quad_index)   
             index1  =  jj * ( obj%spline_degree1 + 1 ) + ii + 1
             M_rho_loc(index1)= M_rho_loc(index1) + &
                  val_f * val_jac * weight * spline
                
          end do
       end do
    end do

    
  end subroutine build_local_matrices_rho
  
  subroutine local_to_global_matrices( &
       solv,&
       cell_index, &
       cell_i, &
       cell_j, &
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
    sll_int32 :: cell_index
    sll_int32 :: cell_i
    sll_int32 :: cell_j
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
    sll_int32 :: index1, index2, index3, index4
    sll_int32 :: index_coef1,index_coef2,index
    sll_int32 :: i,j,mm, nn, b, bprime,x,y
    sll_int32 :: li_A, li_Aprime
    sll_real64 :: elt_mat_global
    sll_int32 :: nbsp,nbsp1
    sll_int32 :: bc_left
    sll_int32 :: bc_right
    sll_int32 :: bc_bottom
    sll_int32 :: bc_top
    
    bc_left   = solv%bc_left
    bc_right  = solv%bc_right
    bc_bottom = solv%bc_bottom
    bc_top    = solv%bc_top
    
    
    
    do mm = 0,solv%spline_degree2
       index3 = cell_j + mm
       
       if (  (bc_bottom==SLL_PERIODIC).and.(bc_top== SLL_PERIODIC) ) then 
          
          if ( index3 > solv%total_num_splines_eta2) then
             index3 = index3 - solv%total_num_splines_eta2
          end if
          
       end if
       
       do i = 0,solv%spline_degree1
          
          index1 = cell_i + i
          if ( (bc_left==SLL_PERIODIC).and.(bc_right== SLL_PERIODIC)) then 
             if ( index1 > solv%total_num_splines_eta1) then
                
                index1 = index1 - solv%total_num_splines_eta1
                
             end if
             nbsp = solv%total_num_splines_eta1
             
          else if ( (bc_left  == SLL_DIRICHLET).and.&
               (bc_right == SLL_DIRICHLET) ) then
             nbsp = solv%num_cells1 + solv%spline_degree1
          end if
          x          =  index1 + (index3-1)*nbsp
          b          =  mm * ( solv%spline_degree1 + 1 ) + i + 1
          li_A       =  solv%local_to_global_spline_indices(b, cell_index)
          
          Masse_tot(x)    = Masse_tot(x) + Masse_loc(b)
          Stiff_tot(x)    = Stiff_tot(x) + Stiff_loc(b)
          
          do nn = 0,solv%spline_degree2
             
             index4 = cell_j + nn
             
             if ( (bc_bottom==SLL_PERIODIC).and.(bc_top== SLL_PERIODIC))then
                if ( index4 > solv%total_num_splines_eta2) then
                   
                   index4 = index4 - solv%total_num_splines_eta2
                end if
             end if
             
             do j = 0,solv%spline_degree1
                
                index2 = cell_i + j
                if((bc_left==SLL_PERIODIC).and.(bc_right==SLL_PERIODIC))then
                   
                   if ( index2 > solv%total_num_splines_eta1) then
                      
                      index2 = index2 - solv%total_num_splines_eta1
                   end if
                   nbsp1 = solv%total_num_splines_eta1
                   
                else if ( (bc_left  == SLL_DIRICHLET) .and.&
                     (bc_right == SLL_DIRICHLET) ) then
                   
                   nbsp1 = solv%num_cells1 + solv%spline_degree1
                end if
                
                y         = index2 + (index4-1)*nbsp1
                bprime    =  nn * ( solv%spline_degree1 + 1 ) + j + 1
                li_Aprime = solv%local_to_global_spline_indices(bprime,cell_index)
                elt_mat_global = &
                     M_c_loc(b, bprime)     - &
                     K_a11_loc(b, bprime)   - &
                     K_a12_loc(b, bprime)   - &
                     K_a21_loc(b, bprime)   - &
                     K_a22_loc(b, bprime)   - &
                     M_b_vect_loc(b,bprime) - &
                     S_b1_loc( b, bprime)   - &
                     S_b2_loc( b, bprime)
                index_coef1 = solv%tab_index_coeff1(cell_i)- solv%spline_degree1 + i
                index_coef2 = solv%tab_index_coeff2(cell_j)- solv%spline_degree2+ mm
                index = index_coef1 + (index_coef2-1)*(solv%num_cells1+1)
                solv%local_to_global_spline_indices_source(b,cell_index)= index
                
                solv%local_to_global_spline_indices_source_bis(bprime,cell_index)= y!index_phi
                
                if ( (li_A > 0) .and. (li_Aprime > 0) ) then
                   
                   call sll_add_to_csr_matrix( &
                        solv%sll_csr_mat, &
                        elt_mat_global, &
                        li_A, &
                        li_Aprime)   


                end if
                
             end do
             
          end do
       end do
    end do
    
  end subroutine local_to_global_matrices
  
  
  subroutine local_to_global_matrices_rho( &
       solv,&
       cell_i, &
       cell_j, &
       M_rho_loc)
    
     class(finite_elements_solver)  :: solv
     sll_int32 :: cell_i
     sll_int32 :: cell_j
     sll_real64, dimension(:), intent(in)   :: M_rho_loc
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
     
     
     
    do mm = 0,solv%spline_degree2
       index3 = cell_j + mm
       if (  (bc_bottom==SLL_PERIODIC).and.(bc_top== SLL_PERIODIC) ) then    
          if ( index3 > solv%total_num_splines_eta2) then
             index3 = index3 - solv%total_num_splines_eta2
          end if
       end if
!other option for above:      index3 = mod(index3 - 1, solv%total_num_splines_eta2) + 1
       
       do i = 0,solv%spline_degree1
          
          index1 = cell_i + i
          if ( (bc_left==SLL_PERIODIC).and.(bc_right== SLL_PERIODIC)) then 
             if ( index1 > solv%total_num_splines_eta1) then
                
                index1 = index1 - solv%total_num_splines_eta1
                
             end if
             nbsp = solv%total_num_splines_eta1
             
          else if ( (bc_left  == SLL_DIRICHLET).and.&
               (bc_right == SLL_DIRICHLET) ) then
             nbsp = solv%num_cells1 + solv%spline_degree1
          end if
          x          =  index1 + (index3-1)*nbsp
          b          =  mm * ( solv%spline_degree1 + 1 ) + i + 1
          solv%rho_vec(x)  =  solv%rho_vec(x)  + M_rho_loc(b)
          
       end do
    end do
    
  end subroutine local_to_global_matrices_rho
  
  subroutine solve_linear_system( solv )
    ! CSR_MAT*phi = rho_vec is the linear system to be solved. The solution
    ! is given in terms of the spline coefficients that represent phi.
    class(finite_elements_solver) :: solv
    !type(csr_matrix)  :: csr_masse
    integer :: elt, elt1
    integer :: i,j
     character(len=*),parameter :: as_file='rho', as_file1='phi',as_file2='mat'
    sll_int32 :: bc_left
    sll_int32 :: bc_right
    sll_int32 :: bc_bottom
    sll_int32 :: bc_top
    sll_int32 :: ierr

    bc_left   = solv%bc_left
    bc_right  = solv%bc_right
    bc_bottom = solv%bc_bottom
    bc_top    = solv%bc_top
    
    solv%tmp_rho_vec(:) = 0.0_f64
    solv%tmp_phi_vec(:) = 0.0_f64
    
    if( (bc_left   == SLL_PERIODIC) .and. (bc_right == SLL_PERIODIC) .and. &
         (bc_bottom == SLL_DIRICHLET).and. (bc_top   == SLL_DIRICHLET) ) then
       
       do i = 1, solv%total_num_splines_eta1
          do j = 1, solv%total_num_splines_eta2
             
             elt  = i + solv%total_num_splines_eta1 * (  j - 1)
             elt1 = i + ( solv%total_num_splines_eta1 ) * j
             solv%tmp_rho_vec(elt) = solv%rho_vec(elt1)
          end do
       end do
       
    else if ( (bc_left   == SLL_DIRICHLET).and.(bc_right==SLL_DIRICHLET) .and.&
         (bc_bottom == SLL_DIRICHLET).and.(bc_top==SLL_DIRICHLET) ) then 
       
       do i = 1, solv%total_num_splines_eta1
          do j = 1, solv%total_num_splines_eta2
             
             elt  = i + solv%total_num_splines_eta1 * (  j - 1)
             elt1 = i + 1 + ( solv%total_num_splines_eta1 + 2 ) * j 
             solv%tmp_rho_vec( elt ) = solv%rho_vec( elt1 )
          end do
       end do
       
    else if((bc_left   == SLL_PERIODIC) .and. (bc_right==SLL_PERIODIC) .and.&
         (bc_bottom == SLL_PERIODIC) .and. (bc_top  ==SLL_PERIODIC)) then
       
       solv%tmp_rho_vec(1:solv%total_num_splines_eta1*solv%total_num_splines_eta2)=&
            solv%rho_vec(1:solv%total_num_splines_eta1*solv%total_num_splines_eta2) 
       
       
    else if( (bc_left == SLL_DIRICHLET) .and. (bc_right == SLL_DIRICHLET) .and.&
             (bc_bottom == SLL_PERIODIC).and. (bc_top   == SLL_PERIODIC) ) then
       
       do i = 1, solv%total_num_splines_eta1
          do j = 1, solv%total_num_splines_eta2

             elt1 = i + 1 + ( solv%total_num_splines_eta1 + 2 ) * (  j - 1)
             elt  = i + solv%total_num_splines_eta1 * (  j - 1)
             solv%tmp_rho_vec( elt ) = solv%rho_vec( elt1 )
          end do
       end do

       
    end if

    call solve_gen_elliptic_eq(solv%sll_csr_mat,solv%tmp_rho_vec,solv%tmp_phi_vec)
    solv%phi_vec(1:solv%total_num_splines_eta1*solv%total_num_splines_eta2)=&
         solv%tmp_phi_vec(1:solv%total_num_splines_eta1*solv%total_num_splines_eta2)

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
    ! CSR_MAT*phi = rho_vec is the linear system to be solved. The solution
    ! is given in terms of the spline coefficients that represent phi.
    class(finite_elements_solver) :: solv
    sll_int32 :: ierr
    
    solv%tmp_rho_vec(:) = 0.0_f64
    solv%tmp_phi_vec(:) = 0.0_f64
    solv%tmp_rho_vec(1:solv%total_num_splines_eta1*solv%total_num_splines_eta2)=&
         solv%rho_vec(1:solv%total_num_splines_eta1*solv%total_num_splines_eta2) 
    call solve_general_es_perper(solv%sll_csr_mat_with_constraint,solv%tmp_rho_vec,solv%tmp_phi_vec) 
    solv%phi_vec(1:solv%total_num_splines_eta1*solv%total_num_splines_eta2) = &
       solv%tmp_phi_vec(1:solv%total_num_splines_eta1*solv%total_num_splines_eta2)    

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

end module finite_elements_solver_module
