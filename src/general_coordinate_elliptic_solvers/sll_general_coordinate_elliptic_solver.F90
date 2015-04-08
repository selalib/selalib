!> @ingroup general_coordinate_elliptic_solvers
!> @brief Elliptic solver on 2d curvilinear mesh
!> @details This solver works with analytical 
!> and discrete coordinate transformations.
module sll_general_coordinate_elliptic_solver_module
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

use sll_boundary_condition_descriptors
use sll_module_scalar_field_2d_base, only: sll_scalar_field_2d_base
use sll_module_scalar_field_2d, only: sll_scalar_field_2d_analytic,  &
                                      sll_scalar_field_2d_discrete
use sll_module_interpolators_2d_base, only: sll_interpolator_2d_base
use sll_module_arbitrary_degree_spline_interpolator_2d, only:        &
  sll_arbitrary_degree_spline_interpolator_2d
use connectivity_module, only: initconnectivity
use sll_knots, only: initialize_knots
use gauss_legendre_integration
use gauss_lobatto_integration
use sll_sparse_matrix_module, only : sll_csr_matrix,                 &
                                     new_csr_matrix,                 &
                                     new_csr_matrix_with_constraint, &
                                     csr_add_one_constraint,         &
                                     sll_factorize_csr_matrix,       &
                                     sll_add_to_csr_matrix,          &
                                     sll_mult_csr_matrix_vector,     &
                                     sll_solve_csr_matrix,           &
                                     sll_delete

implicit none

private

type, public :: general_coordinate_elliptic_solver

  private
  sll_int32,  public :: num_cells1
  sll_int32,  public :: num_cells2
  sll_real64, public :: delta_eta1
  sll_real64, public :: delta_eta2
  sll_real64, public :: eta1_min
  sll_real64, public :: eta2_min   

  sll_int32 :: total_num_splines_loc
  sll_int32 :: total_num_splines1
  sll_int32 :: total_num_splines2
  sll_real64, dimension(:), pointer :: knots1
  sll_real64, dimension(:), pointer :: knots2
  sll_real64, dimension(:), pointer :: knots1_rho
  sll_real64, dimension(:), pointer :: knots2_rho
  sll_real64, dimension(:,:), pointer :: gauss_pts1
  sll_real64, dimension(:,:), pointer :: gauss_pts2
  sll_int32 :: bc_left
  sll_int32 :: bc_left_interp
  sll_int32 :: bc_right
  sll_int32 :: bc_bottom
  sll_int32 :: bc_top
  sll_int32 :: spline_degree1
  sll_int32 :: spline_degree2
  sll_real64 :: epsi
  sll_real64 :: intjac

  ! the following is otherwise known as "ID" in Aurore Back's original
  ! nomenclature. The indexing of the
  ! splines in this array depends on the boundary conditions.
  sll_int32, dimension(:), pointer, public :: global_spline_indices 
  ! the following is otherwise known as "IEN"
  sll_int32, dimension(:,:), pointer :: local_spline_indices
  ! the following is otherwise known as "LM". Same as global_spline_indices
  ! but including the changes resulting from the boundary conditions.
  ! This is:
  ! local_to_lobal_spline_indices(i,j) = 
  !   global_spline_indices(local_spline_indices(i,j))
  sll_int32, dimension(:,:), pointer :: local_to_global_spline_indices
  sll_int32, dimension(:,:), pointer :: local_to_global_spline_indices_source
  sll_int32, dimension(:,:), pointer :: local_to_global_spline_indices_source_bis

  !!! contains the values of all splines in all gauss points
  sll_real64, dimension(:,:), pointer :: values_splines1
  sll_real64, dimension(:,:), pointer :: values_splines2
  sll_real64, dimension(:,:), pointer :: values_splines_gauss1
  sll_real64, dimension(:,:), pointer :: values_splines_gauss2
  sll_real64, dimension(:,:), pointer :: values_jacobian
  sll_int32 , dimension(:)  , pointer :: tab_index_coeff1
  sll_int32 , dimension(:)  , pointer :: tab_index_coeff2
  type(sll_csr_matrix), pointer       :: sll_csr_mat
  type(sll_csr_matrix), pointer       :: sll_csr_mat_with_constraint
  type(sll_csr_matrix), pointer       :: sll_csr_mat_source
  sll_real64, dimension(:), pointer   :: rho_vec
  sll_real64, dimension(:), pointer   :: phi_vec
  sll_real64, dimension(:), pointer   :: tmp_rho_vec
  sll_real64, dimension(:), pointer   :: tmp_phi_vec
  sll_real64, dimension(:), pointer   :: masse
  sll_real64, dimension(:), pointer   :: stiff
  sll_real64, dimension(:),   pointer :: rho_coeff_1d
  sll_real64, dimension(:,:), pointer :: rho_at_gauss
  sll_real64, dimension(:), pointer   :: M_rho_loc
  logical                             :: perper
  !save variables for deboor splines (bsplvd, bsplvb)
  sll_int32                           :: ilo = 1
  sll_int32                           :: jlo = 1
end type general_coordinate_elliptic_solver

!> For the integration mode.  
sll_int32, parameter, public :: ES_GAUSS_LEGENDRE = 0
!> For the integration mode.  
sll_int32, parameter, public :: ES_GAUSS_LOBATTO = 1
  
interface sll_delete
  module procedure delete_elliptic
end interface sll_delete

interface sll_create
  module procedure initialize_general_elliptic_solver
end interface sll_create

interface sll_solve
  module procedure solve_general_coordinates_elliptic_eq
end interface sll_solve

public sll_delete,                          &
       sll_create,                          &
       sll_solve,                           &
       new_general_elliptic_solver,         &
       factorize_mat_es

contains 

! *******************************************************************

!> @brief Initialization for elliptic solver.
!> @details To have the function phi such that 
!> \f[
!>  \nabla \cdot ( A \nabla \phi ) + B \nabla \phi + C \phi = \rho
!> \f]
!>  where A is a matrix of functions , B a vectorial function,
!>  and  C and rho a scalar function.  
!>  A, B, C, rho can be discret or analytic. 
!>  \f$ \phi \f$ is given with a B-spline interpolator  
!> 
!> The parameters are
!> @param      es the type general_coordinate_elliptic_solver
!> @param[in]  spline_degree1 the degre of B-spline in the direction eta1
!> @param[in]  spline_degree2 the degre of B-spline in the direction eta2
!> @param[in]  num_cells1 the number of cells in the direction eta1
!> @param[in]  num_cells2 the number of cells in the direction eta2
!> @param[in]  quadrature_type1 the type of quadrature in the direction eta1
!> @param[in]  quadrature_type2 the type of quadrature in the direction eta2
!> @param[in]  bc_left the boundary condition at left in the direction eta1
!> @param[in]  bc_right the boundary condition at right in the direction eta2
!> @param[in]  bc_bottom the boundary condition at left in the direction eta2
!> @param[in]  bc_top the boundary condition at right in the direction eta2
!> @param[in]  eta1_min the minimun in the direction eta1
!> @param[in]  eta1_max the maximun in the direction eta1
!> @param[in]  eta2_min the minimun in the direction eta2
!> @param[in]  eta2_max the maximun in the direction eta2
!> @param[out] the type general_coordinate_elliptic_solver
subroutine initialize_general_elliptic_solver( &
       es, &
       spline_degree1, &
       spline_degree2, &
       num_cells1, &
       num_cells2, &
       quadrature_type1, &
       quadrature_type2, &
       bc_left, &
       bc_right, &
       bc_bottom, &
       bc_top, &
       eta1_min, &
       eta1_max, &
       eta2_min, &
       eta2_max)
    
 type(general_coordinate_elliptic_solver), intent(out) :: es

 sll_int32,  intent(in) :: spline_degree1
 sll_int32,  intent(in) :: spline_degree2
 sll_int32,  intent(in) :: num_cells1
 sll_int32,  intent(in) :: num_cells2
 sll_int32,  intent(in) :: bc_left
 sll_int32,  intent(in) :: bc_right
 sll_int32,  intent(in) :: bc_bottom
 sll_int32,  intent(in) :: bc_top
 sll_int32,  intent(in) :: quadrature_type1
 sll_int32,  intent(in) :: quadrature_type2
 sll_real64, intent(in) :: eta1_min
 sll_real64, intent(in) :: eta1_max
 sll_real64, intent(in) :: eta2_min
 sll_real64, intent(in) :: eta2_max

 sll_int32 :: knots1_size
 sll_int32 :: knots2_size
 sll_int32 :: num_splines1
 sll_int32 :: num_splines2
 sll_int32 :: vec_sz ! for rho_vec and phi_vec allocations
 sll_int32 :: ierr
 sll_int32 :: solution_size,i
 sll_int32 :: dim1, dim2
 sll_int32 :: num_pts_g1, num_pts_g2
 
 sll_int32 :: bc_left_knots
 
 bc_left_knots = bc_left
   
 es%total_num_splines_loc = (spline_degree1+1)*(spline_degree2+1)
 ! The total number of splines in a single direction is given by
 ! num_cells + spline_degree
 num_splines1 = num_cells1 + spline_degree1
 num_splines2 = num_cells2 + spline_degree2
 SLL_ALLOCATE(es%global_spline_indices(num_splines1*num_splines2),ierr)
 es%global_spline_indices(:) = 0
   
 dim1 = (spline_degree1+1)*(spline_degree2+1)
 dim2 = (num_cells1*num_cells2)
 SLL_ALLOCATE(es%local_spline_indices(1:dim1,1:dim2),ierr)
 es%local_spline_indices = 0
 SLL_ALLOCATE(es%local_to_global_spline_indices(1:dim1,1:dim2),ierr)
 es%local_to_global_spline_indices = 0
 SLL_ALLOCATE(es%local_to_global_spline_indices_source(1:dim1,1:dim2),ierr)
 SLL_ALLOCATE(es%local_to_global_spline_indices_source_bis(1:dim1,1:dim2),ierr)

 ! This should be changed to verify that the passed BC's are part of the
 ! recognized list described in sll_boundary_condition_descriptors...

 es%bc_left        = bc_left
 es%bc_right       = bc_right
 es%bc_bottom      = bc_bottom
 es%bc_top         = bc_top
 es%spline_degree1 = spline_degree1
 es%spline_degree2 = spline_degree2
 es%num_cells1     = num_cells1
 es%num_cells2     = num_cells2
 es%delta_eta1     = (eta1_max-eta1_min)/num_cells1
 es%delta_eta2     = (eta2_max-eta2_min)/num_cells2
 es%eta1_min       = eta1_min
 es%eta2_min       = eta2_min
 
 es%bc_left_interp = bc_left
 if(bc_left==SLL_NEUMANN)then
   es%bc_left_interp = SLL_DIRICHLET
 endif

 !quadrature_type1_tmp = ES_GAUSS_LEGENDRE
 !quadrature_type2_tmp = ES_GAUSS_LEGENDRE
 ! Allocate and fill the gauss points/weights information.
 ! First direction
 select case(quadrature_type1)
 case (ES_GAUSS_LEGENDRE)
   SLL_ALLOCATE(es%gauss_pts1(2,spline_degree1+2),ierr)
   es%gauss_pts1(:,:) = 0.0_f64
   es%gauss_pts1(:,:) = gauss_legendre_points_and_weights(spline_degree1+2)
 case (ES_GAUSS_LOBATTO)
   SLL_ALLOCATE(es%gauss_pts1(2,spline_degree1+2),ierr)
   es%gauss_pts1(:,:) = gauss_lobatto_points_and_weights(spline_degree1+2)
 case default
   print *, 'new_general_qn_solver(): &
             & have not type of gauss points in the first direction'
 end select
    
 select case(quadrature_type2)
 case (ES_GAUSS_LEGENDRE)
   SLL_ALLOCATE(es%gauss_pts2(2,spline_degree2+2),ierr)
   es%gauss_pts2(:,:) = 0.0_f64
   es%gauss_pts2(:,:) = gauss_legendre_points_and_weights(spline_degree2+2)
 case (ES_GAUSS_LOBATTO)
   SLL_ALLOCATE(es%gauss_pts2(2,spline_degree2+2),ierr)
   es%gauss_pts2(:,:) = gauss_lobatto_points_and_weights(spline_degree2+2)
   case DEFAULT
     print *, 'new_general_qn_solver(): &
                & have not type of gauss points in the second direction'
      
 end select


 es%perper  = .false. 

 if( (bc_left   == SLL_PERIODIC) .and. (bc_right == SLL_PERIODIC) .and. &
     (bc_bottom == SLL_PERIODIC) .and. (bc_top   == SLL_PERIODIC) ) then

   es%total_num_splines1 = num_cells1 
   es%total_num_splines2 = num_cells2
 
   knots1_size = 2*spline_degree1+2
   knots2_size = 2*spline_degree2+2
   vec_sz      = num_cells1*num_cells2
   es%perper   = .true. 

 else if( (bc_left   == SLL_PERIODIC)  .and. (bc_right == SLL_PERIODIC) .and.&
          (bc_bottom == SLL_DIRICHLET) .and. (bc_top   == SLL_DIRICHLET) ) then

   es%total_num_splines1 = num_cells1 
   es%total_num_splines2 = num_cells2 + spline_degree2 - 2

   knots1_size = 2*spline_degree1+2
   knots2_size = 2*spline_degree2+num_cells2+1

   vec_sz      = num_cells1*(num_cells2+spline_degree2)

 else if( (bc_left   == SLL_DIRICHLET) .and. (bc_right == SLL_DIRICHLET) .and.&
          (bc_bottom == SLL_PERIODIC)  .and. (bc_top == SLL_PERIODIC) ) then

   es%total_num_splines1 = num_cells1 + spline_degree1 - 2
   es%total_num_splines2 = num_cells2 
   knots1_size = 2*spline_degree1+num_cells1+1
   knots2_size = 2*spline_degree2+2
   vec_sz      = (num_cells1 + spline_degree1)*num_cells2

 else if( (bc_left   == SLL_NEUMANN) .and. (bc_right == SLL_DIRICHLET) .and.&
          (bc_bottom == SLL_PERIODIC)  .and. (bc_top == SLL_PERIODIC) ) then

   es%total_num_splines1 = num_cells1 + spline_degree1 - 1
   es%total_num_splines2 = num_cells2 
   knots1_size = 2*spline_degree1+num_cells1+1
   knots2_size = 2*spline_degree2+2
   vec_sz      = (num_cells1 + spline_degree1)*num_cells2
   bc_left_knots = SLL_DIRICHLET

 else if( (bc_left   == SLL_DIRICHLET) .and. (bc_right == SLL_DIRICHLET) .and.&
          (bc_bottom == SLL_DIRICHLET) .and. (bc_top   == SLL_DIRICHLET) ) then

   es%total_num_splines1 = num_cells1 + spline_degree1 - 2
   es%total_num_splines2 = num_cells2 + spline_degree2 - 2
   knots1_size = 2*spline_degree1 + num_cells1+1
   knots2_size = 2*spline_degree2 + num_cells2+1
   vec_sz      = (num_cells1 + spline_degree1)*&
                 (num_cells2 + spline_degree2)

 end if

 solution_size = es%total_num_splines1*es%total_num_splines2
 SLL_ALLOCATE(es%knots1(knots1_size),ierr)
 SLL_ALLOCATE(es%knots2(knots2_size),ierr)
 SLL_ALLOCATE(es%knots1_rho(num_cells1 + spline_degree1 + 2),ierr)
 SLL_ALLOCATE(es%knots2_rho(num_cells2 + spline_degree2 + 2),ierr)
 SLL_ALLOCATE(es%rho_vec(vec_sz),ierr)
 SLL_ALLOCATE(es%phi_vec(solution_size),ierr)
 SLL_ALLOCATE(es%masse(vec_sz),ierr)
 SLL_ALLOCATE(es%stiff(vec_sz),ierr)

 ! -------------------------------------------
 ! We must add plus 1 for the dimension of the solution 
 ! in the case periodic periodic to include the periodicity in the last point.  
 !  -----------------------------------------

 !if(es%perper) solution_size = solution_size + 1
 if(es%perper) then
   SLL_ALLOCATE(es%tmp_rho_vec(solution_size+1),ierr)
   SLL_ALLOCATE(es%tmp_phi_vec(solution_size+1),ierr)
 else
   SLL_ALLOCATE(es%tmp_rho_vec(solution_size),ierr)
   SLL_ALLOCATE(es%tmp_phi_vec(solution_size),ierr) 
 endif
 es%rho_vec = 0.0_f64
 es%phi_vec = 0.0_f64
 es%masse   = 0.0_f64
 es%stiff   = 0.0_f64
 es%intjac  = 0.0_f64

 call initialize_knots( &
   spline_degree1,  &
   num_cells1,      &
   eta1_min,            &
   eta1_max,            &
   bc_left_knots,       &
   bc_right,            &
   es%knots1 )

 call initialize_knots( &
   spline_degree2,  &
   num_cells2,      &
   eta2_min,            &
   eta2_max,            &
   bc_bottom,           &
   bc_top,              &
   es%knots2 )

 call initconnectivity(               &
   num_cells1,                    &
   num_cells2,                    &
   spline_degree1,                &
   spline_degree2,                &
   bc_left,                           &
   bc_right,                          &
   bc_bottom,                         &
   bc_top,                            &
   es%local_spline_indices,           &
   es%global_spline_indices,          &
   es%local_to_global_spline_indices )
   
 es%sll_csr_mat => new_csr_matrix(    &
   solution_size,                     &
   solution_size,                     &
   num_cells1*num_cells2,     &
   es%local_to_global_spline_indices, &
   es%total_num_splines_loc,          &
   es%local_to_global_spline_indices, &
   es%total_num_splines_loc)

 es%knots1_rho(1:spline_degree1+1) = eta1_min
 es%knots1_rho(num_cells1+2:num_cells1+1+spline_degree1+1) = eta1_max
  
 if ( mod(spline_degree1 +1,2) == 0 ) then
   do i = spline_degree1 +1 + 1, num_cells1 + 1
     es%knots1_rho(i) = eta1_min + (i-(spline_degree1 +1)/2-1 )*es%delta_eta1 
   end do
 else
   do i = spline_degree1 +1 + 1, num_cells1 + 1
     es%knots1_rho ( i ) = &
       0.5*( eta1_min + ( i - (spline_degree1)/2 -1)*es%delta_eta1 + &
       eta1_min +  ( i -1 - (spline_degree1)/2 -1)*es%delta_eta1 )
   end do
   
 end if

 es%knots2_rho(1:spline_degree2+1) = eta2_min
 es%knots2_rho(num_cells2+2:num_cells2+1+spline_degree2+1)=eta2_max
     
 if (mod(spline_degree2+1,2) == 0 ) then
   do i = spline_degree2 +1 + 1, num_cells2 + 1
     es%knots2_rho(i) = eta2_min + (i-(spline_degree2+1)/2-1)*es%delta_eta2 
   end do
 else
   do i = spline_degree2 +1 + 1, num_cells2 + 1
     es%knots2_rho ( i ) = &
       0.5*( eta2_min+(i  -(spline_degree2)/2-1)*es%delta_eta2 + &
             eta2_min+(i-1-(spline_degree2)/2-1)*es%delta_eta2 )
   end do
 end if

 !! allocation of the table containning all values of splines in each direction 
 !! in each gauss points
 SLL_ALLOCATE(es%values_splines1(num_cells1*(spline_degree1+2),spline_degree1+1),ierr)
 es%values_splines1 = 0.0_f64
 SLL_ALLOCATE(es%values_splines2(num_cells2*(spline_degree2+2),spline_degree2+1),ierr)
 es%values_splines2 = 0.0_f64
 SLL_ALLOCATE(es%values_jacobian(num_cells1*(spline_degree1+2),num_cells2*(spline_degree2+2)),ierr)
 es%values_jacobian = 0.0_f64
 SLL_ALLOCATE(es%values_splines_gauss1(num_cells1*(spline_degree1+2),spline_degree1+1),ierr)
 es%values_splines_gauss1 = 0.0_f64
 SLL_ALLOCATE(es%values_splines_gauss2(num_cells2*(spline_degree2+2),spline_degree2+1),ierr)
 es%values_splines_gauss2 = 0.0_f64
 SLL_ALLOCATE(es%tab_index_coeff1(num_cells1*(spline_degree1+2)),ierr)
 SLL_ALLOCATE(es%tab_index_coeff2(num_cells2*(spline_degree2+2)),ierr)

  num_pts_g1 = size(es%gauss_pts1,2)
  num_pts_g2 = size(es%gauss_pts2,2)
  SLL_ALLOCATE(es%rho_at_gauss(num_cells1*num_pts_g1,num_cells2*num_pts_g2),ierr)
  SLL_ALLOCATE(es%rho_coeff_1d((num_cells1+1)*(num_cells2+1)),ierr)
  SLL_ALLOCATE(es%M_rho_loc(es%total_num_splines_loc),ierr)

end subroutine initialize_general_elliptic_solver
  
!> @brief Initialization for elliptic solver.
!> @details To have the function phi such that 
!> \f[
!>  \nabla \cdot( A \nabla \phi ) + B \nabla \phi + C \phi = \rho
!> \f]
!>  where A is a matrix of functions , B a vectorial function,
!>  and  C and rho a scalar function.  
!>  A, B, C, \f$ rho \f$ can be discret or analytic. 
!>  \f$ phi \f$ is given with a B-spline interpolator  
!> 
!> The parameters are
!> @param[in] spline_degree1 the degre of B-spline in the direction eta1
!> @param[in] spline_degree2 the degre of B-spline in the direction eta2
!> @param[in] num_cells1 the number of cells in the direction eta1
!> @param[in] num_cells2 the number of cells in the direction eta2
!> @param[in] quadrature_type1 the type of quadrature in the direction eta1
!> @param[in] quadrature_type2 the type of quadrature in the direction eta2
!> @param[in] bc_left the boundary condition at left in the direction eta1
!> @param[in] bc_right the boundary condition at right in the direction eta2
!> @param[in] bc_bottom the boundary condition at left in the direction eta2
!> @param[in] bc_top the boundary condition at right in the direction eta2
!> @param[in] eta1_min the minimun in the direction eta1
!> @param[in] eta1_max the maximun in the direction eta1
!> @param[in] eta2_min the minimun in the direction eta2
!> @param[in] eta2_max the maximun in the direction eta2
!> @return the type general_coordinate_elliptic_solver such that a pointer
function new_general_elliptic_solver( &
  spline_degree1, &
  spline_degree2, &
  num_cells1, &
  num_cells2, &
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

  type(general_coordinate_elliptic_solver), pointer :: es
  sll_int32,  intent(in) :: spline_degree1
  sll_int32,  intent(in) :: spline_degree2
  sll_int32,  intent(in) :: num_cells1
  sll_int32,  intent(in) :: num_cells2
  sll_int32,  intent(in) :: bc_left
  sll_int32,  intent(in) :: bc_right
  sll_int32,  intent(in) :: bc_bottom
  sll_int32,  intent(in) :: bc_top
  sll_int32,  intent(in) :: quadrature_type1
  sll_int32,  intent(in) :: quadrature_type2
  sll_real64, intent(in) :: eta1_min
  sll_real64, intent(in) :: eta1_max
  sll_real64, intent(in) :: eta2_min
  sll_real64, intent(in) :: eta2_max
  sll_int32 :: ierr


  SLL_ALLOCATE(es,ierr)
  call sll_create( es,               &
                   spline_degree1,   &
                   spline_degree2,   &
                   num_cells1,       &
                   num_cells2,       &
                   quadrature_type1, &
                   quadrature_type2, &
                   bc_left,          &
                   bc_right,         &
                   bc_bottom,        &
                   bc_top,           &
                   eta1_min,         &
                   eta1_max,         &
                   eta2_min,         &
                   eta2_max )
   
end function new_general_elliptic_solver

!> @brief Deallocate the type general_coordinate_elliptic_solver
!> @details
!> The parameters are
!> @param[in] es the type general_coordinate_elliptic_solver
  
subroutine delete_elliptic( es )
  type(general_coordinate_elliptic_solver) :: es
  sll_int32 :: ierr
  ! it is not good to check some cases and not others, fix...
  if(associated(es%knots1)) then
    SLL_DEALLOCATE(es%knots1,ierr)
  else
    print *, 'delete es, WARNING: knots1 array was not allocated.'
  end if
  if(associated(es%knots2)) then
    SLL_DEALLOCATE(es%knots2,ierr)
  else
    print *, 'delete es general coords, ', &
             'WARNING: knots2 array was not allocated.'
  end if
  SLL_DEALLOCATE(es%gauss_pts1,ierr)
  SLL_DEALLOCATE(es%gauss_pts2,ierr)
  SLL_DEALLOCATE(es%global_spline_indices,ierr)
  SLL_DEALLOCATE(es%local_spline_indices,ierr)
  SLL_DEALLOCATE(es%local_to_global_spline_indices,ierr)
  SLL_DEALLOCATE(es%local_to_global_spline_indices_source,ierr)
  SLL_DEALLOCATE(es%local_to_global_spline_indices_source_bis,ierr)
  call sll_delete(es%sll_csr_mat)
  call sll_delete(es%sll_csr_mat_with_constraint)
  call sll_delete(es%sll_csr_mat_source)
  SLL_DEALLOCATE(es%rho_vec,ierr)
  SLL_DEALLOCATE(es%phi_vec,ierr)
  SLL_DEALLOCATE(es%tmp_rho_vec,ierr)
  SLL_DEALLOCATE(es%tmp_phi_vec,ierr)
  SLL_DEALLOCATE(es%masse,ierr)
  SLL_DEALLOCATE(es%stiff,ierr)
  SLL_DEALLOCATE(es%knots1_rho,ierr)
  SLL_DEALLOCATE(es%knots2_rho,ierr)
  SLL_DEALLOCATE(es%values_splines1,ierr)
  SLL_DEALLOCATE(es%values_splines2,ierr)
  SLL_DEALLOCATE(es%values_jacobian,ierr)
  SLL_DEALLOCATE(es%values_splines_gauss1,ierr)
  SLL_DEALLOCATE(es%values_splines_gauss2,ierr)
  SLL_DEALLOCATE(es%tab_index_coeff1,ierr)
  SLL_DEALLOCATE(es%tab_index_coeff2,ierr)

end subroutine delete_elliptic



!> @brief Assemble the matrix for elliptic solver.
!> @details To have the function phi such that 
!> \f[
!>  \nabla \cdot ( A \nabla \phi ) + B \nabla \phi + C \phi = \rho
!> \f]
!>  where A is a matrix of functions , B a vectorial function,
!>  and  C and rho a scalar function.  
!>  A, B, C, rho can be discret or analytic. 
!>  phi is given with a B-spline interpolator  
!> 
!> The parameters are
!> @param es the type general_coordinate_elliptic_solver
!> @param[in] a11_field_mat the field corresponding to the matrix coefficient A(1,1)
!> @param[in] a12_field_mat the field corresponding to the matrix coefficient A(1,2)
!> @param[in] a21_field_mat the field corresponding to the matrix coefficient A(2,1)
!> @param[in] a22_field_mat the field corresponding to the matrix coefficient A(2,2)
!> @param[in] b1_field_vect the field corresponding to the vector coefficient B(1)
!> @param[in] b2_field_vect the field corresponding to the vector coefficient B(2)
!> @param[in] c_field the field corresponding to the coefficient B(1) of the scalar C
!> @return the type general_coordinate_elliptic_solver contains the matrix 
!> to solve the equation
subroutine factorize_mat_es( es,            &
                             a11_field_mat, &
                             a12_field_mat, &
                             a21_field_mat, &
                             a22_field_mat, &
                             b1_field_vect, &
                             b2_field_vect, &
                             c_field)

  type(general_coordinate_elliptic_solver),intent(inout) :: es

  class(sll_scalar_field_2d_base), pointer :: a11_field_mat
  class(sll_scalar_field_2d_base), pointer :: a12_field_mat
  class(sll_scalar_field_2d_base), pointer :: a21_field_mat
  class(sll_scalar_field_2d_base), pointer :: a22_field_mat
  class(sll_scalar_field_2d_base), pointer :: b1_field_vect
  class(sll_scalar_field_2d_base), pointer :: b2_field_vect
  class(sll_scalar_field_2d_base), pointer :: c_field

  sll_real64, dimension(:,:),   allocatable :: M_c_loc
  sll_real64, dimension(:,:),   allocatable :: K_a11_loc
  sll_real64, dimension(:,:),   allocatable :: K_a12_loc
  sll_real64, dimension(:,:),   allocatable :: K_a21_loc
  sll_real64, dimension(:,:),   allocatable :: K_a22_loc
  sll_real64, dimension(:,:),   allocatable :: M_b_vect_loc
  sll_real64, dimension(:,:),   allocatable :: S_b1_loc
  sll_real64, dimension(:,:),   allocatable :: S_b2_loc  
  sll_real64, dimension(:),     allocatable :: Masse_loc
  sll_real64, dimension(:),     allocatable :: Stiff_loc
  sll_real64, dimension(:,:,:), pointer     :: Source_loc

  sll_int32 :: total_num_splines_loc
  sll_int32 :: ierr
  sll_int32 :: i
  sll_int32 :: j
  sll_int32 :: cell_index
  sll_int32 :: bc_left
  sll_int32 :: bc_right
  sll_int32 :: bc_bottom
  sll_int32 :: bc_top

  character(len=*),parameter :: as_file1='mat'
    
  bc_left   = es%bc_left
  bc_right  = es%bc_right
  bc_bottom = es%bc_bottom
  bc_top    = es%bc_top

  total_num_splines_loc = es%total_num_splines_loc

  if( (bc_left   == SLL_PERIODIC) .and. (bc_right == SLL_PERIODIC) .and. &
      (bc_bottom == SLL_PERIODIC) .and. (bc_top   == SLL_PERIODIC) ) then
     es%perper = .true.
  else
     es%perper = .false.  
  end if   

  SLL_CLEAR_ALLOCATE(Source_loc(1:es%num_cells1*es%num_cells2,1:total_num_splines_loc,1:total_num_splines_loc),ierr)
  SLL_CLEAR_ALLOCATE(M_c_loc(1:total_num_splines_loc,1:total_num_splines_loc),ierr)
  SLL_CLEAR_ALLOCATE(K_a11_loc(1:total_num_splines_loc,1:total_num_splines_loc),ierr)
  SLL_CLEAR_ALLOCATE(K_a12_loc(1:total_num_splines_loc,1:total_num_splines_loc),ierr)
  SLL_CLEAR_ALLOCATE(K_a21_loc(1:total_num_splines_loc,1:total_num_splines_loc),ierr)
  SLL_CLEAR_ALLOCATE(K_a22_loc(1:total_num_splines_loc,1:total_num_splines_loc),ierr)
  SLL_CLEAR_ALLOCATE(S_b1_loc(1:total_num_splines_loc,1:total_num_splines_loc),ierr)
  SLL_CLEAR_ALLOCATE(S_b2_loc(1:total_num_splines_loc,1:total_num_splines_loc),ierr)
  SLL_CLEAR_ALLOCATE(M_b_vect_loc(1:total_num_splines_loc,1:total_num_splines_loc),ierr)
  SLL_CLEAR_ALLOCATE(Masse_loc(1:total_num_splines_loc),ierr)
  SLL_CLEAR_ALLOCATE(Stiff_loc(1:total_num_splines_loc),ierr)

  do j=1,es%num_cells2
    do i=1,es%num_cells1
          
      cell_index = i+es%num_cells1*(j-1)
      
      
      call build_local_matrices( &
               es,               &
               cell_index,       &
               i,                &
               j,                &
               a11_field_mat,    &
               a12_field_mat,    &
               a21_field_mat,    &
               a22_field_mat,    &
               b1_field_vect,    &
               b2_field_vect,    &
               c_field,          &
               Masse_loc,        &
               Stiff_loc,        &
               M_c_loc,          &
               K_a11_loc,        &
               K_a12_loc,        &
               K_a21_loc,        &
               K_a22_loc,        &
               M_b_vect_loc,     &
               S_b1_loc,         &
               S_b2_loc,         &
               Source_loc)


      call local_to_global_matrices( &
               es,                   &
               cell_index,           &
               i,                    &
               j,                    &
               Masse_loc,            &
               Stiff_loc,            &
               M_c_loc,              &
               K_a11_loc,            &
               K_a12_loc,            &
               K_a21_loc,            &
               K_a22_loc,            &
               M_b_vect_loc,         &
               S_b1_loc,             &
               S_b2_loc,             &
               es%masse,             &
               es%stiff)
          
    end do
  end do
 
  if (es%perper) then

   es%sll_csr_mat_with_constraint => new_csr_matrix_with_constraint(es%sll_csr_mat)  

   call csr_add_one_constraint( es%sll_csr_mat%row_ptr,                 &  
                                es%sll_csr_mat%col_ind,                 &
                                es%sll_csr_mat%val,                     &
                                es%sll_csr_mat%num_rows,                &
                                es%sll_csr_mat%num_nz,                  &
                                es%masse,                               &
                                es%sll_csr_mat_with_constraint%row_ptr, &
                                es%sll_csr_mat_with_constraint%col_ind, &
                                es%sll_csr_mat_with_constraint%val)  

!    es%sll_csr_mat_with_constraint => es%sll_csr_mat
!
!    do i = 1, es%total_num_splines1*es%total_num_splines2    
!
!      call sll_add_to_csr_matrix( &
!           es%sll_csr_mat_with_constraint, &
!           es%masse(i), &
!           es%total_num_splines1*es%total_num_splines2+1, &
!           i)
!    end do
    print *,'#begin of sll_factorize_csr_matrix'
    call sll_factorize_csr_matrix(es%sll_csr_mat_with_constraint)
    print *,'#end of sll_factorize_csr_matrix'
  else   
    print *,'#begin of sll_factorize_csr_matrix'
    call sll_factorize_csr_matrix(es%sll_csr_mat)
    print *,'#end of sll_factorize_csr_matrix'
  end if 
 
  es%sll_csr_mat_source => new_csr_matrix( size(es%masse,1),          &
                        (es%num_cells1+1)*(es%num_cells2+1),          &
                        es%num_cells1*es%num_cells2,                  &
                        es%local_to_global_spline_indices_source_bis, &
                        es%total_num_splines_loc,                     &
                        es%local_to_global_spline_indices_source,     &
                        es%total_num_splines_loc )

  call compute_source_matrice(es,source_loc)
    
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
   
end subroutine factorize_mat_es
  
!> @brief Assemble the matrix for elliptic solver.
!> @details To have the function phi such that 
!>  \f[
!>  \nabla \cdot ( A \nabla \phi ) + B \nabla \phi + C \phi = \rho
!>  \f]
!>  where A is a matrix of functions , B a vectorial function,
!>  and  C and rho a scalar function.  
!>  A, B, C, \f$ \rho \f$  can be discret or analytic. 
!>  phi is given with a B-spline interpolator  
!> 
!> The parameters are
!> @param es the type general_coordinate_elliptic_solver
!> @param[in] rho \f$ \rho \f$ the field corresponding to the source term   
!> @param[out] phi \f$ \phi \f$ the field corresponding to the solution of the equation
!> @return phi the field solution of the equation
  
subroutine solve_general_coordinates_elliptic_eq( es, rho, phi)

  class(general_coordinate_elliptic_solver) :: es
  class(sll_scalar_field_2d_discrete), intent(inout)  :: phi
  class(sll_scalar_field_2d_base), intent(in),target  :: rho
  sll_int32 :: i
  sll_int32 :: j
  sll_int32 :: cell_index
  sll_int32 :: total_num_splines_loc
  sll_real64 :: int_rho,int_jac
  sll_int32 :: num_pts_g1, num_pts_g2, ig1, ig2, ig, jg
  sll_real64 :: wgpt1, wgpt2, gpt1, gpt2, eta1, eta2
  sll_real64 :: val_jac
  class(sll_scalar_field_2d_base),pointer  :: base_field_pointer
  class(sll_interpolator_2d_base),pointer  :: base_interpolator_pointer
  sll_real64, dimension(:,:), pointer :: coeff_rho
    
  total_num_splines_loc = es%total_num_splines_loc
    
  num_pts_g1 = size(es%gauss_pts1,2)
  num_pts_g2 = size(es%gauss_pts2,2)
  es%rho_at_gauss(:,:) = 0.0_f64   

  es%M_rho_loc     = 0.0_f64
  es%rho_vec(:) = 0.0_f64
  es%rho_coeff_1d  = 0.0_f64
   
  !ES Compute rho at all Gauss points
  ig1 = 0 
  ig2 = 0 
  int_rho = 0.0_f64
  int_jac = 0.0_f64
    
  ! afin d'optimiser la construction de la matrice rho
  ! on va proceder de la facon qui suit 
  base_field_pointer => rho
  select type( type_field => base_field_pointer)
  class is (sll_scalar_field_2d_discrete)
    base_interpolator_pointer => type_field%interp_2d

    select type( type_interpolator => base_interpolator_pointer)
    class is (sll_arbitrary_degree_spline_interpolator_2d)
      coeff_rho => type_interpolator%get_coefficients()
                
      ! put the spline coefficients in a 1d array
      do j=1,es%num_cells2+1
        do i=1,es%num_cells1+1
          es%rho_coeff_1d(i+(es%num_cells1+1)*(j-1)) = coeff_rho(i,j)
        end do
      end do

      call sll_mult_csr_matrix_vector(&
            es%sll_csr_mat_source, es%rho_coeff_1d,es%rho_vec)

      if(es%perper) then
        es%rho_vec = es%rho_vec - sum(es%rho_vec)/es%intjac*es%masse
      end if
        
   class DEFAULT
          
     do j=1,es%num_cells2
       eta2  = es%eta2_min + (j-1)*es%delta_eta2
       do i=1,es%num_cells1
         eta1  = es%eta1_min + (i-1)*es%delta_eta1
         do jg=1,num_pts_g2
           ! rescale Gauss points to be in interval [eta2 ,eta2 +delta_eta2]
           ! the bottom edge of the cell.
           gpt2  = eta2  + 0.5_f64*es%delta_eta2 * &
                    ( es%gauss_pts2(1,jg) + 1.0_f64 )
           wgpt2 = 0.5_f64*es%delta_eta2*es%gauss_pts2(2,jg) !ATTENTION 0.5
               
           ig2 = jg + (j-1)*num_pts_g2
           do ig=1,num_pts_g1
             ! rescale Gauss points to be in interval [eta1,eta1+delta1]
             gpt1  = eta1  + 0.5_f64*es%delta_eta1 * &
                     ( es%gauss_pts1(1,ig) + 1.0_f64 )
             wgpt1 = 0.5_f64*es%delta_eta1*es%gauss_pts1(2,ig)
             ig1 = ig + (i-1)*num_pts_g1
             es%rho_at_gauss(ig1,ig2)   = rho%value_at_point(gpt1,gpt2)
             val_jac = &
               es%values_jacobian(i+es%num_cells1*(ig-1),j+es%num_cells2*(jg-1))
             int_rho = int_rho + &
                       rho%value_at_point(gpt1,gpt2)*wgpt2*wgpt1*val_jac 
             int_jac = int_jac + wgpt2*wgpt1*val_jac
                   
           end do
         end do
       end do
     end do

     if (es%perper) then
       es%rho_at_gauss = es%rho_at_gauss - int_rho/int_jac
     end if
          
     do j=1,es%num_cells2
       do i=1,es%num_cells1
               
         cell_index = i+es%num_cells1*(j-1)
         call build_local_matrices_rho(es, i, j)
         call local_to_global_matrices_rho(es, i, j)
       end do
     end do
   end select
       
   class is (sll_scalar_field_2d_analytic)
    
     do j=1,es%num_cells2
       eta2  = es%eta2_min + (j-1)*es%delta_eta2
       do i=1,es%num_cells1
         eta1  = es%eta1_min + (i-1)*es%delta_eta1
         do jg=1,num_pts_g2
           ! rescale Gauss points to be in interval [eta2 ,eta2 +delta_eta2]
           ! the bottom edge of the cell.
           gpt2  = eta2  + 0.5_f64*es%delta_eta2 * &
                 ( es%gauss_pts2(1,jg) + 1.0_f64 )
           wgpt2 = 0.5_f64*es%delta_eta2*es%gauss_pts2(2,jg) !ATTENTION 0.5
               
           ig2 = jg + (j-1)*num_pts_g2
           do ig=1,num_pts_g1
             ! rescale Gauss points to be in interval [eta1,eta1+delta1]
             gpt1  = eta1  + 0.5_f64*es%delta_eta1 * &
                        ( es%gauss_pts1(1,ig) + 1.0_f64 )
             wgpt1 = 0.5_f64*es%delta_eta1*es%gauss_pts1(2,ig)
             ig1 = ig + (i-1)*num_pts_g1
             es%rho_at_gauss(ig1,ig2)   = rho%value_at_point(gpt1,gpt2)
             val_jac = &
               es%values_jacobian(i+es%num_cells1*(ig-1),j+es%num_cells2*(jg-1))
             int_rho = int_rho &
               + rho%value_at_point(gpt1,gpt2)*wgpt2*wgpt1*val_jac 
             int_jac = int_jac + wgpt2*wgpt1*val_jac
                  
           end do
         end do
       end do
     end do

    if (es%perper) then
       es%rho_at_gauss = es%rho_at_gauss - int_rho/int_jac
    end if
          
    do j=1,es%num_cells2
      do i=1,es%num_cells1
             
        cell_index = i+es%num_cells1*(j-1)
        call build_local_matrices_rho(es, i, j)
        call local_to_global_matrices_rho(es, i, j)
      end do
    end do
       
  end select

  call solve_linear_system(es)
  
  
  print *,'#solve_linear_system done'
    
  call  phi%interp_2d%set_coefficients( es%phi_vec)
!  call  phi%interp_2d%set_coefficients( &
!    es%phi_vec(1:(es%total_num_splines1-1)*es%total_num_splines2))

end subroutine solve_general_coordinates_elliptic_eq
  
subroutine build_local_matrices( &
       es,                       &
       cell_index,               &
       cell_i,                   &
       cell_j,                   &
       a11_field_mat,            &
       a12_field_mat,            &
       a21_field_mat,            &
       a22_field_mat,            &
       b1_field_vect,            &
       b2_field_vect,            &
       c_field,                  &
       Masse_loc,                &
       Stiff_loc,                &
       M_c_loc,                  &
       K_a11_loc,                &
       K_a12_loc,                &
       K_a21_loc,                &
       K_a22_loc,                &
       M_b_vect_loc,             &
       S_b1_loc,                 &
       S_b2_loc,                 &
       Source_loc)

  type(general_coordinate_elliptic_solver) :: es  
  sll_int32, intent(in) :: cell_i
  sll_int32, intent(in) :: cell_j
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
  sll_real64 :: eta1
  sll_real64 :: eta2
  sll_int32  :: tmp1
  sll_int32  :: tmp2
  sll_int32  :: num_pts_g1 ! number of gauss points in first direction 
  sll_int32  :: num_pts_g2 ! number of gauss points in second direction
  sll_int32  :: i,ii,iii
  sll_int32  :: j,jj,jjj
  sll_real64 :: gpt1
  sll_real64 :: gpt2
  sll_real64 :: wgpt1
  sll_real64 :: wgpt2
  sll_real64 :: gtmp1
  sll_real64 :: gtmp2
  sll_int32  :: local_spline_index1
  sll_int32  :: local_spline_index2
  sll_int32  :: index1
  sll_int32  :: index2
  sll_real64, dimension(es%spline_degree1+1,es%spline_degree1+1) :: work1
  sll_real64, dimension(es%spline_degree2+1,es%spline_degree2+1) :: work2
  sll_real64, dimension(es%spline_degree1+1,2) :: dbiatx1
  sll_real64, dimension(es%spline_degree2+1,2) :: dbiatx2
  sll_real64, dimension(es%spline_degree1+1,2) :: dbiatx1_rho
  sll_real64, dimension(es%spline_degree2+1,2) :: dbiatx2_rho
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
  delta1    = es%delta_eta1 !mesh2d%delta_eta1
  delta2    = es%delta_eta2 !! mesh2d%delta_eta2
  eta1_min  = es%eta1_min  
  eta2_min  = es%eta2_min  
  tmp1      = (es%spline_degree1 + 1)/2
  tmp2      = (es%spline_degree2 + 1)/2
  bc_left   = es%bc_left
  bc_right  = es%bc_right
  bc_bottom = es%bc_bottom
  bc_top    = es%bc_top
  num_pts_g1 = size(es%gauss_pts1,2)
  num_pts_g2 = size(es%gauss_pts2,2)
    
  eta1  = eta1_min + (cell_i-1)*delta1
  eta2  = eta2_min + (cell_j-1)*delta2

  do j=1,num_pts_g2
    ! rescale Gauss points to be in interval [eta2 ,eta2 +delta_eta2]
    ! the bottom edge of the cell.
    gpt2  = eta2  + 0.5_f64*delta2 * ( es%gauss_pts2(1,j) + 1.0_f64 )
    wgpt2 = 0.5_f64*delta2*es%gauss_pts2(2,j) !ATTENTION 0.5
       
    if ((es%bc_bottom==SLL_PERIODIC).and.(es%bc_top==SLL_PERIODIC))then
      ! rescale gauss point in interval [0,delta2]
      gtmp2 = 0.5_f64*delta2*( es%gauss_pts2(1,j) + 1.0_f64) !ATTENTION 0.5
      local_spline_index2 = es%spline_degree2 + 1
          
    else if ((es%bc_bottom == SLL_DIRICHLET).and.&
             (es%bc_top    == SLL_DIRICHLET)) then
      gtmp2 = gpt2
      local_spline_index2 = es%spline_degree2 + cell_j
    end if
       
    call bsplvd(                  &
           es,                   & 
           es%knots2,            &
           es%spline_degree2+1,  &
           gtmp2,                 &
           local_spline_index2,   &
           work2,                 &
           dbiatx2,               &
           2)

    call interv(                                    &
           es,                                     & 
           es%knots2_rho,                          &
           es%num_cells2 + es%spline_degree2+ 2,  &
           gpt2, left_y, mflag_y )

    call bsplvd(                                  &
         es,                                     & 
         es%knots2_rho,                          &
         es%spline_degree2+1,                    &
         gpt2,                                    &
         left_y,                                  &
         work2,                                   &
         dbiatx2_rho,                             &
         2)
      
    ! we stock values of spline to construct the source term
    es%values_splines2(cell_j + es%num_cells2*(j-1),:) = dbiatx2(:,1)
    es%values_splines_gauss2(cell_j+es%num_cells2*(j-1),:)=dbiatx2_rho(:,1)
    es%tab_index_coeff2(cell_j + es%num_cells2*(j-1)) = left_y

    do i=1,num_pts_g1
      ! rescale Gauss points to be in interval [eta1,eta1+delta1]
      gpt1  = eta1  + 0.5_f64*delta1 * ( es%gauss_pts1(1,i) + 1.0_f64 )
      wgpt1 = 0.5_f64*delta1*es%gauss_pts1(2,i)
          
      if((es%bc_left_interp==SLL_PERIODIC).and.(es%bc_right==SLL_PERIODIC)) then 
             
        gtmp1   = 0.5_f64*delta1*( es%gauss_pts1(1,i) + 1.0_f64)! ATTENTION 0.5 
        local_spline_index1 = es%spline_degree1 + 1
             
      else if ((es%bc_left_interp  == SLL_DIRICHLET).and.&
        (es%bc_right == SLL_DIRICHLET) ) then
           
        gtmp1   = gpt1
        local_spline_index1 = es%spline_degree1 + cell_i
             
      end if
    
      call bsplvd(                   &
             es,                    & 
             es%knots1,             &
             es%spline_degree1+1,   &
             gtmp1,                  &
             local_spline_index1,    &
             work1,                  &
             dbiatx1,                &
             2 )

      call interv( es,              &
                   es%knots1_rho,   &
                   es%num_cells1 + es%spline_degree1+ 2, gpt1, &
                   left_x, mflag_x )
     
      call bsplvd( es,                    &
                   es%knots1_rho,         &
                   es%spline_degree1+1,   &
                   gpt1,                   &
                   left_x,                 &
                   work1,                  &
                   dbiatx1_rho,            &
                   2 )

      es%values_splines1(cell_i + es%num_cells1*(i-1),:) = dbiatx1(:,1)
      es%values_splines_gauss1(cell_i + es%num_cells1*(i-1),:) = dbiatx1_rho(:,1)
      es%tab_index_coeff1(cell_i + es%num_cells1*(i-1)) = left_x

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
      val_jac = jac_mat(1,1)*jac_mat(2,2) - jac_mat(1,2)*jac_mat(2,1)

      es%values_jacobian(cell_i + es%num_cells1*(i-1),cell_j &
                                 + es%num_cells2*(j-1)) = val_jac
        
      es%intjac = es%intjac + wgpt2*wgpt1*val_jac
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
          
      MC =  jac_mat(2,2) * val_b1_der1 &
          - jac_mat(2,1) * val_b1_der2 &
          - jac_mat(1,2) * val_b2_der1 &
          + jac_mat(1,1) * val_b2_der2
          

      C1 =  jac_mat(2,2) * val_b1 &
          - jac_mat(1,2) * val_b2 
      C2 =  jac_mat(1,1) * val_b2 &
          - jac_mat(2,1) * val_b1
         
      ! loop over the splines supported in the cell that are different than
      ! zero at the point (gpt1,gpt2) (there are spline_degree+1 splines in
      ! each direction.
      do ii = 0,es%spline_degree1
         do jj = 0,es%spline_degree2
              
           index1  =  jj * ( es%spline_degree1 + 1 ) + ii + 1
              
           Masse_loc(index1) = &
             Masse_loc(index1) + &
             val_jac*wgpt1*wgpt2* &
             (dbiatx1(ii+1,1)*dbiatx2(jj+1,1))

           Stiff_loc(index1) = &
             Stiff_loc(index1) + &
             val_jac*wgpt1*wgpt2* &
             (dbiatx1(ii+1,2)*dbiatx2(jj+1,1)+&
             dbiatx1(ii+1,1)*dbiatx2(jj+1,2))
             
           do iii = 0,es%spline_degree1
             do jjj = 0,es%spline_degree2
                    
               index2 =  jjj*(es%spline_degree1 + 1) + iii + 1
              
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
  
subroutine build_local_matrices_rho( es, cell_i, cell_j)

  class(general_coordinate_elliptic_solver) :: es
  sll_int32, intent(in) :: cell_i
  sll_int32, intent(in) :: cell_j
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
  sll_int32  :: num_pts_g1 
  sll_int32  :: num_pts_g2
  sll_int32  :: i,ii
  sll_int32  :: j,jj
  sll_real64 :: gpt1
  sll_real64 :: gpt2
  sll_real64 :: wgpt1
  sll_real64 :: wgpt2
  sll_int32  :: index1
  sll_real64, dimension(es%spline_degree1+1,es%spline_degree1+1) :: work1
  sll_real64, dimension(es%spline_degree2+1,es%spline_degree2+1) :: work2
  sll_real64, dimension(es%spline_degree1+1,2) :: dbiatx1
  sll_real64, dimension(es%spline_degree2+1,2) :: dbiatx2
  sll_real64 :: val_f
  sll_real64 :: val_jac,spline1,spline2
    
  es%M_rho_loc(:)  = 0.0_f64
  dbiatx1(:,:)  = 0.0_f64
  dbiatx2(:,:)  = 0.0_f64
  work1(:,:)    = 0.0_f64
  work2(:,:)    = 0.0_f64
  ! The supposition is that all fields use the same logical mesh
  delta1    = es%delta_eta1
  delta2    = es%delta_eta2
  eta1_min  = es%eta1_min
  eta2_min  = es%eta2_min
  tmp1      = (es%spline_degree1 + 1)/2
  tmp2      = (es%spline_degree2 + 1)/2
  bc_left   = es%bc_left
  bc_right  = es%bc_right
  bc_bottom = es%bc_bottom
  bc_top    = es%bc_top
  num_pts_g1 = size(es%gauss_pts1,2)
  num_pts_g2 = size(es%gauss_pts2,2)
    
  eta1  = eta1_min + (cell_i-1)*delta1
  eta2  = eta2_min + (cell_j-1)*delta2
    
  do j=1,num_pts_g2
    ! rescale Gauss points to be in interval [eta2 ,eta2 +delta_eta2]
    ! the bottom edge of the cell.
    gpt2  = eta2  + 0.5_f64*delta2 * ( es%gauss_pts2(1,j) + 1.0_f64 )
    wgpt2 = 0.5_f64*delta2*es%gauss_pts2(2,j) !ATTENTION 0.5

    do i=1,num_pts_g1
      ! rescale Gauss points to be in interval [eta1,eta1+delta1]
      gpt1  = eta1  + 0.5_f64*delta1 * ( es%gauss_pts1(1,i) + 1.0_f64 )
      wgpt1 = 0.5_f64*delta1*es%gauss_pts1(2,i)
  
      val_f = es%rho_at_gauss(i+(cell_i-1)*num_pts_g1, j + (cell_j-1)*num_pts_g2)
      val_jac = &
        es%values_jacobian(cell_i+es%num_cells1*(i-1),cell_j+es%num_cells2*(j-1))

     ! loop over the splines supported in the cell that are different than
     ! zero at the point (gpt1,gpt2) (there are spline_degree+1 splines in
     ! each direction.
      do ii = 0,es%spline_degree1
        do jj = 0,es%spline_degree2
                
          spline1 = es%values_splines1(cell_i + es%num_cells1*(i-1),ii+1)
          spline2 = es%values_splines2(cell_j + es%num_cells2*(j-1),jj+1)
               
          index1  =  jj * ( es%spline_degree1 + 1 ) + ii + 1
          es%M_rho_loc(index1)= es%M_rho_loc(index1) + &
                    val_f*val_jac*wgpt1*wgpt2*spline1*spline2
                
        end do
      end do
    end do
  end do
    
end subroutine build_local_matrices_rho
  
subroutine local_to_global_matrices( &
       es,                           &
       cell_index,                   &
       cell_i,                       &
       cell_j,                       &
       Masse_loc,                    &
       Stiff_loc,                    &
       M_c_loc,                      &
       K_a11_loc,                    &
       K_a12_loc,                    &
       K_a21_loc,                    &
       K_a22_loc,                    &
       M_b_vect_loc,                 &
       S_b1_loc,                     &
       S_b2_loc,                     &
       Masse_tot,                    &
       Stiff_tot)
    
  class(general_coordinate_elliptic_solver)  :: es
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


  bc_left   = es%bc_left_interp
  bc_right  = es%bc_right
  bc_bottom = es%bc_bottom
  bc_top    = es%bc_top
    
  do mm = 0,es%spline_degree2
    index3 = cell_j + mm
      
    if ((bc_bottom==SLL_PERIODIC).and.(bc_top== SLL_PERIODIC) ) then 
          
      if ( index3 > es%total_num_splines2) then
        index3 = index3 - es%total_num_splines2
      end if
          
    end if
       
    do i = 0,es%spline_degree1
          
      index1 = cell_i + i
      if ( (bc_left==SLL_PERIODIC).and.(bc_right== SLL_PERIODIC)) then 
        if ( index1 > es%total_num_splines1) then
                
          index1 = index1 - es%total_num_splines1
          
        end if
        nbsp = es%total_num_splines1
             
      else if ( (bc_left  == SLL_DIRICHLET).and.&
                (bc_right == SLL_DIRICHLET) ) then
        nbsp = es%num_cells1 + es%spline_degree1
      end if

      x            =  index1 + (index3-1)*nbsp
      b            =  mm * ( es%spline_degree1 + 1 ) + i + 1
      li_A         =  es%local_to_global_spline_indices(b, cell_index)
         
      Masse_tot(x) = Masse_tot(x) + Masse_loc(b)
      Stiff_tot(x) = Stiff_tot(x) + Stiff_loc(b)

      index_coef1 = es%tab_index_coeff1(cell_i)- es%spline_degree1 + i
      index_coef2 = es%tab_index_coeff2(cell_j)- es%spline_degree2+ mm
      index = index_coef1 + (index_coef2-1)*(es%num_cells1+1)

      es%local_to_global_spline_indices_source(b,cell_index)= index

      do nn = 0,es%spline_degree2
             
        index4 = cell_j + nn
             
        if ( (bc_bottom==SLL_PERIODIC).and.(bc_top== SLL_PERIODIC))then
          if ( index4 > es%total_num_splines2) then
            index4 = index4 - es%total_num_splines2
          end if
        end if
             
        do j = 0,es%spline_degree1
                
          index2 = cell_i + j
          if((bc_left==SLL_PERIODIC).and.(bc_right==SLL_PERIODIC))then
                   
           if ( index2 > es%total_num_splines1) then
             index2 = index2 - es%total_num_splines1
           end if
           nbsp1 = es%total_num_splines1
                   
        else if ( (bc_left  == SLL_DIRICHLET) .and.&
                  (bc_right == SLL_DIRICHLET) ) then
          nbsp1 = es%num_cells1 + es%spline_degree1
        end if
                
        y         = index2 + (index4-1)*nbsp1
        bprime    =  nn * ( es%spline_degree1 + 1 ) + j + 1
        li_Aprime = es%local_to_global_spline_indices(bprime,cell_index)
        elt_mat_global = &
                     M_c_loc(b, bprime)     - &
                     K_a11_loc(b, bprime)   - &
                     K_a12_loc(b, bprime)   - &
                     K_a21_loc(b, bprime)   - &
                     K_a22_loc(b, bprime)   - &
                     M_b_vect_loc(b,bprime) - &
                     S_b1_loc( b, bprime)   - &
                     S_b2_loc( b, bprime)
        es%local_to_global_spline_indices_source_bis(bprime,cell_index)= y!index_phi

        if ( (li_A > 0) .and. (li_Aprime > 0) ) then
                   
          call sll_add_to_csr_matrix(   &
                        es%sll_csr_mat, &
                        elt_mat_global, &
                        li_A, &
                        li_Aprime)   

        end if
                
      end do
             
    end do
  end do
end do
    
end subroutine local_to_global_matrices

subroutine local_to_global_matrices_rho( es, cell_i, cell_j) 

  class(general_coordinate_elliptic_solver)  :: es
  sll_int32 :: cell_i
  sll_int32 :: cell_j
  sll_int32 :: i,mm, b, x!,y
  sll_int32 :: nbsp!,nbsp1
  sll_int32 :: bc_left
  sll_int32 :: bc_right
  sll_int32 :: bc_bottom
  sll_int32 :: bc_top
  sll_int32 :: index1,index3
   
  bc_left   = es%bc_left_interp
  bc_right  = es%bc_right
  bc_bottom = es%bc_bottom
  bc_top    = es%bc_top
  
  do mm = 0,es%spline_degree2
    index3 = cell_j + mm
    if (  (bc_bottom==SLL_PERIODIC).and.(bc_top== SLL_PERIODIC) ) then    
       if ( index3 > es%total_num_splines2) then
          index3 = index3 - es%total_num_splines2
       end if
    end if
     
!other option for above:      index3 = mod(index3 - 1, es%total_num_splines2) + 1

    do i = 0,es%spline_degree1
          
      index1 = cell_i + i
      if ( (bc_left==SLL_PERIODIC).and.(bc_right== SLL_PERIODIC)) then 
        if ( index1 > es%total_num_splines1) then
              
          index1 = index1 - es%total_num_splines1
             
        end if
        nbsp = es%total_num_splines1
           
     else if ( (bc_left  == SLL_DIRICHLET).and.&
       (bc_right == SLL_DIRICHLET) ) then
       nbsp = es%num_cells1 + es%spline_degree1
     end if

     x          =  index1 + (index3-1)*nbsp
     b          =  mm * ( es%spline_degree1 + 1 ) + i + 1
     es%rho_vec(x)  =  es%rho_vec(x)  + es%M_rho_loc(b)
        
   end do
end do
  
end subroutine local_to_global_matrices_rho
  
!> @details
!> CSR_MAT*phi = rho_vec is the linear system to be solved. The solution
!> is given in terms of the spline coefficients that represent phi.
subroutine solve_linear_system( es )
  class(general_coordinate_elliptic_solver) :: es
  !type(csr_matrix)  :: csr_masse
  integer :: elt, elt1
  integer :: i,j
  character(len=*),parameter :: as_file  = 'rho'
  character(len=*),parameter :: as_file1 = 'phi'
  character(len=*),parameter :: as_file2 = 'mat'
  sll_int32 :: bc_left
  sll_int32 :: bc_right
  sll_int32 :: bc_bottom
  sll_int32 :: bc_top

  bc_left   = es%bc_left
  bc_right  = es%bc_right
  bc_bottom = es%bc_bottom
  bc_top    = es%bc_top
   
  es%tmp_rho_vec(:) = 0.0_f64
  es%tmp_phi_vec(:) = 0.0_f64
    
  if( (bc_left   == SLL_PERIODIC) .and. (bc_right == SLL_PERIODIC) .and. &
      (bc_bottom == SLL_DIRICHLET).and. (bc_top   == SLL_DIRICHLET) ) then
       
    do i = 1, es%total_num_splines1
      do j = 1, es%total_num_splines2
             
        elt  = i + es%total_num_splines1 * (  j - 1)
        elt1 = i + ( es%total_num_splines1 ) * j
        es%tmp_rho_vec(elt) = es%rho_vec(elt1)

      end do
    end do
       
  else if ( (bc_left   == SLL_DIRICHLET).and.(bc_right==SLL_DIRICHLET) .and.&
            (bc_bottom == SLL_DIRICHLET).and.(bc_top  ==SLL_DIRICHLET) ) then 
       
    do i = 1, es%total_num_splines1
      do j = 1, es%total_num_splines2
            
        elt  = i + es%total_num_splines1 * (  j - 1)
        elt1 = i + 1 + ( es%total_num_splines1 + 2 ) * j 
        es%tmp_rho_vec( elt ) = es%rho_vec( elt1 )

      end do
    end do
       
  else if((bc_left   == SLL_PERIODIC) .and. (bc_right==SLL_PERIODIC) .and.&
          (bc_bottom == SLL_PERIODIC) .and. (bc_top  ==SLL_PERIODIC)) then
       
    es%tmp_rho_vec(1:es%total_num_splines1*es%total_num_splines2)=&
      es%rho_vec(1:es%total_num_splines1*es%total_num_splines2) 
       
  else if( (bc_left   == SLL_DIRICHLET) .and. (bc_right == SLL_DIRICHLET) .and.&
           (bc_bottom == SLL_PERIODIC)  .and. (bc_top   == SLL_PERIODIC) ) then
     
    do i = 1, es%total_num_splines1
      do j = 1, es%total_num_splines2

        elt1 = i + 1 + ( es%total_num_splines1 + 2 ) * (  j - 1)
        elt  = i + es%total_num_splines1 * (  j - 1)
        es%tmp_rho_vec( elt ) = es%rho_vec( elt1 )

      end do
    end do

  else if( (bc_left   == SLL_NEUMANN) .and. (bc_right == SLL_DIRICHLET) .and.&
           (bc_bottom == SLL_PERIODIC)  .and. (bc_top   == SLL_PERIODIC) ) then
     
    do i = 1, es%total_num_splines1
      do j = 1, es%total_num_splines2

        elt1 = i + 1 + ( es%total_num_splines1 + 1 ) * (  j - 1)
        elt  = i + es%total_num_splines1 * (  j - 1)
        es%tmp_rho_vec( elt ) = es%rho_vec( elt1 )

      end do
    end do

     
  end if

  if((bc_left   == SLL_PERIODIC) .and. (bc_right==SLL_PERIODIC) .and.&
          (bc_bottom == SLL_PERIODIC) .and. (bc_top  ==SLL_PERIODIC)) then
    call sll_solve_csr_matrix(es%sll_csr_mat_with_constraint, &
                            es%tmp_rho_vec, &
                            es%tmp_phi_vec)

  else
    call sll_solve_csr_matrix(es%sll_csr_mat, &
                            es%tmp_rho_vec, &
                            es%tmp_phi_vec)
  endif
  es%phi_vec(1:es%total_num_splines1*es%total_num_splines2)=&
    es%tmp_phi_vec(1:es%total_num_splines1*es%total_num_splines2)

end subroutine solve_linear_system
  
subroutine compute_Source_matrice(es,Source_loc)

  type(general_coordinate_elliptic_solver),intent(inout) :: es
  sll_real64, dimension(:,:,:), pointer :: Source_loc
  sll_int32 :: cell_j,cell_i
  sll_int32 :: cell_index
  sll_int32 :: ideg2,ideg1
  sll_int32 :: jdeg2,jdeg1
  sll_int32 :: b, bprime
  sll_int32 :: li_A,li_Aprime
  sll_real64:: elt_mat_global

  do cell_j=1,es%num_cells2
    do cell_i=1,es%num_cells1
        
      cell_index = cell_i+es%num_cells1*(cell_j-1)
        
      do ideg2 = 0,es%spline_degree2
           
        do ideg1 = 0,es%spline_degree1
              
          b    = ideg2 * ( es%spline_degree1 + 1 ) + ideg1 + 1
          li_A = es%local_to_global_spline_indices_source_bis(b, cell_index)
              
          do jdeg2 = 0,es%spline_degree2
               
             do jdeg1 = 0,es%spline_degree1
                    
               bprime    =  jdeg2 * ( es%spline_degree1 + 1 ) + jdeg1 + 1
               li_Aprime = es%local_to_global_spline_indices_source(bprime,cell_index)
                    
               elt_mat_global = Source_loc(cell_index,bprime,b)

               if ( (li_A > 0) .and. (li_Aprime > 0)) then

                  call sll_add_to_csr_matrix( &
                       es%sll_csr_mat_source, &
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

subroutine interv( es, xt, lxt, x, left, mflag )
    
    type(general_coordinate_elliptic_solver) :: es

    sll_int32,intent(in):: lxt
    sll_int32,intent(out):: left
    sll_int32,intent(out):: mflag
    sll_int32:: ihi
    sll_int32:: istep
    sll_int32:: middle
    sll_real64,intent(in) ::x
    sll_real64,dimension(:):: xt!(lxt)

    
    ihi = es%ilo + 1
    
    if ( lxt <= ihi ) then
       
       if ( xt(lxt) <= x ) then
          go to 110
       end if
       
       if ( lxt <= 1 ) then
          mflag = -1
          left = 1
          return
       end if
       
       es%ilo = lxt - 1
       ihi = lxt
       
    end if
    
    if ( xt(ihi) <= x ) then
       go to 20
    end if
    
    if ( xt(es%ilo) <= x ) then
       mflag = 0
       left = es%ilo
       return
    end if
    !
    !  Now X < XT(ILO).  Decrease ILO to capture X.
    !
    istep = 1
    
10  continue
    
    ihi = es%ilo
    es%ilo = ihi - istep
    
    if ( 1 < es%ilo ) then
       if ( xt(es%ilo) <= x ) then
          go to 50
       end if
       istep = istep * 2
       go to 10
    end if
    
    es%ilo = 1
    
    if ( x < xt(1) ) then
       mflag = -1
       left = 1
       return
    end if
    
    go to 50
    !
    !  Now XT(IHI) <= X.  Increase IHI to capture X.
    !
20  continue
    
    istep = 1
    
30  continue
    
    es%ilo = ihi
    ihi = es%ilo + istep
    
    if ( ihi < lxt ) then
       
       if ( x < xt(ihi) ) then
          go to 50
       end if
       
       istep = istep * 2
       go to 30
       
    end if
    
    if ( xt(lxt) <= x ) then
       go to 110
    end if
    !
    !  Now XT(ILO) < = X < XT(IHI).  Narrow the interval.
    !
    ihi = lxt
    
50  continue
    
    do
       
       middle = ( es%ilo + ihi ) / 2
       
       if ( middle == es%ilo ) then
          mflag = 0
          left = es%ilo
          return
       end if
       !
       !  It is assumed that MIDDLE = ILO in case IHI = ILO+1.
       !
       if ( xt(middle) <= x ) then
          es%ilo = middle
       else
          ihi = middle
       end if
       
    end do
    !
    !  Set output and return.
    !
    
    
110 continue
    
    mflag = 1
    
    if ( x == xt(lxt) ) then
       mflag = 0
    end if
    
    do left = lxt, 1, -1
       if ( xt(left) < xt(lxt) ) then
          return
       end if
    end do
    
    return

  end subroutine interv

 subroutine bsplvd ( es, t, k, x, left, a, dbiatx, nderiv )

    type(general_coordinate_elliptic_solver) :: es
    sll_int32 :: k
    sll_int32 :: left
    sll_int32 :: nderiv
    
    sll_real64, dimension(k,k):: a!(k,k)
    sll_real64,dimension(k,nderiv), intent(out) :: dbiatx!(k,nderiv)
    sll_real64:: factor
    sll_real64:: fkp1mm
    sll_int32 :: i
    sll_int32 :: ideriv
    sll_int32 :: il
    sll_int32 :: j
    sll_int32 :: jlow
    sll_int32 :: jp1mid
    sll_int32 :: ldummy
    sll_int32 :: m
    sll_int32 :: mhigh
    !  sll_real64 sum1  ! this one is not used...
    sll_real64,dimension(left+k):: t ! (left+k)
    sll_real64:: x
    
    
    mhigh = max ( min ( nderiv, k ), 1 )
    !
    !  MHIGH is usually equal to NDERIV.
    !
    call bsplvb ( es, t, k+1-mhigh, 1, x, left, dbiatx )
    
    if ( mhigh == 1 ) then
       return
    end if
    !
    !  The first column of DBIATX always contains the B-spline values
    !  for the current order.  These are stored in column K+1-current
    !  order before BSPLVB is called to put values for the next
    !  higher order on top of it.
    !
    ideriv = mhigh
    do m = 2, mhigh
       jp1mid = 1
       do j = ideriv, k
          dbiatx(j,ideriv) = dbiatx(jp1mid,1)
          jp1mid = jp1mid + 1
          
       end do
       ideriv = ideriv - 1
       
       call bsplvb ( es, t, k+1-ideriv, 2, x, left, dbiatx )
       
    end do
    !
    !  At this point, B(LEFT-K+I, K+1-J)(X) is in DBIATX(I,J) for
    !  I=J,...,K and J=1,...,MHIGH ('=' NDERIV).
    !
    !  In particular, the first column of DBIATX is already in final form.
    !
    !  To obtain corresponding derivatives of B-splines in subsequent columns,
    !  generate their B-representation by differencing, then evaluate at X.
    !
    jlow = 1
    do i = 1, k
       a(jlow:k,i) = 0.0D+00
       jlow = i
       a(i,i) = 1.0D+00
    end do
    !
    !  At this point, A(.,J) contains the B-coefficients for the J-th of the
    !  K B-splines of interest here.
    !
    do m = 2, mhigh
       
       fkp1mm = real ( k + 1 - m, kind = 8 )
       il = left
       i = k
       !
       !  For J = 1,...,K, construct B-coefficients of (M-1)st derivative of
       !  B-splines from those for preceding derivative by differencing
       !  and store again in  A(.,J).  The fact that  A(I,J) = 0 for
       !  I < J is used.
       !
       do ldummy = 1, k+1-m
          
          factor = fkp1mm / ( t(il+k+1-m) - t(il) )
          !
          !  The assumption that T(LEFT) < T(LEFT+1) makes denominator
          !  in FACTOR nonzero.
          !
          a(i,1:i) = ( a(i,1:i) - a(i-1,1:i) ) * factor
          
          il = il - 1
          i = i - 1
       
       end do
       !
       !  For I = 1,...,K, combine B-coefficients A(.,I) with B-spline values
       !  stored in DBIATX(.,M) to get value of (M-1)st derivative of
       !  I-th B-spline (of interest here) at X, and store in DBIATX(I,M).
       !
       !  Storage of this value over the value of a B-spline
       !  of order M there is safe since the remaining B-spline derivatives
       !  of the same order do not use this value due to the fact
       !  that  A(J,I) = 0  for J < I.
       !
       do i = 1, k
          
          jlow = max ( i, m )
          
          dbiatx(i,m) = dot_product ( a(jlow:k,i), dbiatx(jlow:k,m) )
          
       end do
    
    end do
    return
  end subroutine bsplvd


 subroutine bsplvb ( es, t, jhigh, index, x, left, biatx )

    type(general_coordinate_elliptic_solver) :: es
    sll_int32, parameter :: jmax = 20
    
    sll_int32:: jhigh
    
    sll_real64,dimension(jhigh):: biatx !(jhigh)
    sll_real64, save, dimension ( jmax ) :: deltal
    sll_real64, save, dimension ( jmax ) :: deltar
    sll_int32:: i
    sll_int32:: index
    sll_int32:: left
    sll_real64:: saved
    sll_real64,dimension(left+jhigh):: t!() left+jhigh
    sll_real64:: term
    sll_real64:: x
    
    if ( index == 1 ) then
       es%jlo = 1
       biatx(1) = 1.0_8
       if ( jhigh <= es%jlo ) then
          return
       end if
    end if
    
    if ( t(left+1) <= t(left) ) then
       print*,'x=',x
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'BSPLVB - Fatal error!'
       write ( *, '(a)' ) '  It is required that T(LEFT) < T(LEFT+1).'
       write ( *, '(a,i8)' ) '  But LEFT = ', left
       write ( *, '(a,g14.6)' ) '  T(LEFT) =   ', t(left)
       write ( *, '(a,g14.6)' ) '  T(LEFT+1) = ', t(left+1)
       stop
    end if
    
    do
       
       deltar(es%jlo) = t(left+es%jlo) - x
       deltal(es%jlo) = x - t(left+1-es%jlo)
       
       saved = 0.0_f64
       do i = 1, es%jlo
          term = biatx(i) / ( deltar(i) + deltal(es%jlo+1-i) )
          biatx(i) = saved + deltar(i) * term
          saved = deltal(es%jlo+1-i) * term
       end do
    
       biatx(es%jlo+1) = saved
       es%jlo = es%jlo + 1
       
       if ( jhigh <= es%jlo ) exit
    
    end do
    
  end subroutine bsplvb


end module sll_general_coordinate_elliptic_solver_module


!PN I removed this subroutine because it is not used...
!PN 
!!!> @brief Assemble the matrix for elliptic solver.
!!!> @details To have the function phi such that 
!!!> \f[
!!!>  \nabla \cdot ( A \nabla \phi ) + B \nabla \phi + C \phi = \rho
!!!> \f]
!!!>  where A is a matrix of functions , B a vectorial function,
!!!>  and  C and rho a scalar function.  
!!!>  A, B, C, rho can be discret or analytic. 
!!!>  phi is given with a B-spline interpolator  
!!!> 
!!!> The parameters are
!!!> @param es the type general_coordinate_elliptic_solver
!!!> @param[in] a11_field_mat the field corresponding to the matrix coefficient A(1,1)
!!!> @param[in] a12_field_mat the field corresponding to the matrix coefficient A(1,2)
!!!> @param[in] a21_field_mat the field corresponding to the matrix coefficient A(2,1)
!!!> @param[in] a22_field_mat the field corresponding to the matrix coefficient A(2,2)
!!!> @param[in] b1_field_vect the field corresponding to the vector coefficient B(1)
!!!> @param[in] b2_field_vect the field corresponding to the vector coefficient B(2)
!!!> @param[in] c_field the field corresponding to the coefficient B(1) of the scalar C
!!!> @return the type general_coordinate_elliptic_solver contains the matrix 
!!!> to solve the equation
!!
!!subroutine assemble_mat_es( es,            &
!!                            a11_field_mat, &
!!                            a12_field_mat, &
!!                            a21_field_mat, &
!!                            a22_field_mat, &
!!                            b1_field_vect, &
!!                            b2_field_vect, &
!!                            c_field)
!!
!!  type(general_coordinate_elliptic_solver),intent(inout) :: es
!!
!!  class(sll_scalar_field_2d_base), pointer :: a11_field_mat
!!  class(sll_scalar_field_2d_base), pointer :: a12_field_mat
!!  class(sll_scalar_field_2d_base), pointer :: a21_field_mat
!!  class(sll_scalar_field_2d_base), pointer :: a22_field_mat
!!  class(sll_scalar_field_2d_base), pointer :: b1_field_vect
!!  class(sll_scalar_field_2d_base), pointer :: b2_field_vect
!!  class(sll_scalar_field_2d_base), pointer :: c_field
!!  sll_int32 :: i
!!  sll_int32 :: j
!!  sll_int32 :: cell_index
!!  sll_real64 :: eta1
!!  sll_real64 :: eta2
!!  sll_int32 :: num_pts_g1
!!  sll_int32 :: num_pts_g2
!!  sll_int32 :: ii
!!  sll_int32 :: jj
!!  sll_real64 :: gpt1
!!  sll_real64 :: gpt2
!!  sll_real64, dimension(2,2) :: jac_mat
!!
!!  num_pts_g1 = size(es%gauss_pts1,2)
!!  num_pts_g2 = size(es%gauss_pts2,2)
!!  
!!  do j=1,es%num_cells2
!!    do i=1,es%num_cells1          
!!      cell_index = i+es%num_cells1*(j-1)
!!      eta1  = es%eta1_min + real(i-1,f64)*es%delta_eta1
!!      eta2  = es%eta2_min + real(j-1,f64)*es%delta_eta2
!!
!!      do jj=1,num_pts_g2
!!        gpt2  = eta2  + 0.5_f64*es%delta_eta2 * ( es%gauss_pts2(1,jj) + 1.0_f64 )
!!        do ii=1,num_pts_g1
!!          jac_mat(:,:) = c_field%get_jacobian_matrix(gpt1,gpt2)
!!        enddo
!!      enddo  
!!    enddo
!!  enddo
!!  
!!  print *,'#not implemented for the moment'
!!  
!!end subroutine assemble_mat_es
 
