!> @ingroup general_coordinate_elliptic_solvers
!> @brief Elliptic solver on 2d curvilinear mesh
!> @details This solver works with analytical 
!> and discrete coordinate transformations.
module sll_general_coordinate_elliptic_solver_module
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_errors.h"

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
  sll_real64, dimension(:,:,:), pointer :: v_splines1
  sll_real64, dimension(:,:,:), pointer :: v_splines2
  sll_real64, dimension(:,:,:), pointer :: d_splines1
  sll_real64, dimension(:,:,:), pointer :: d_splines2
  sll_real64, dimension(:,:,:), pointer :: r_splines1
  sll_real64, dimension(:,:,:), pointer :: r_splines2

  sll_real64, dimension(:,:,:,:), pointer :: val_jac
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
  sll_real64, dimension(:,:,:,:), pointer :: rho_at_gauss
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
 sll_int32 :: solution_size
 sll_int32 :: dim1, dim2
 sll_int32 :: num_pts_g1, num_pts_g2
 
 sll_int32 :: bc_left_knots

 sll_real64, allocatable :: work1(:,:)
 sll_real64, allocatable :: work2(:,:)
 sll_real64, allocatable :: dbs1(:,:)
 sll_real64, allocatable :: dbs2(:,:)
 sll_real64 :: xg, wxg, yg, wyg

 sll_int32  :: i, j, ii, jj, ispl1, ispl2
 sll_real64 :: eta1, eta2, gspl1, gspl2
 sll_int32  :: left_1, left_2

 
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
   es%gauss_pts1 = gauss_legendre_points_and_weights(spline_degree1+2)
 case (ES_GAUSS_LOBATTO)
   SLL_ALLOCATE(es%gauss_pts1(2,spline_degree1+2),ierr)
   es%gauss_pts1 = gauss_lobatto_points_and_weights(spline_degree1+2)
 case default
   SLL_ERROR('initialize_general_elliptic_solver','unknown type of gauss points in the direction 1')
 end select
    
 select case(quadrature_type2)
 case (ES_GAUSS_LEGENDRE)
   SLL_ALLOCATE(es%gauss_pts2(2,spline_degree2+2),ierr)
   es%gauss_pts2(:,:) = gauss_legendre_points_and_weights(spline_degree2+2)
 case (ES_GAUSS_LOBATTO)
   SLL_ALLOCATE(es%gauss_pts2(2,spline_degree2+2),ierr)
   es%gauss_pts2(:,:) = gauss_lobatto_points_and_weights(spline_degree2+2)
 case DEFAULT
   SLL_ERROR('initialize_general_elliptic_solver','unknown type of gauss points in the direction 2')
 end select

 es%gauss_pts1(1,:) = 0.5_f64*es%delta_eta1*(es%gauss_pts1(1,:)+1.0_f64)
 es%gauss_pts1(2,:) = 0.5_f64*es%delta_eta1*es%gauss_pts1(2,:)
 es%gauss_pts2(1,:) = 0.5_f64*es%delta_eta2*(es%gauss_pts2(1,:)+1.0_f64)
 es%gauss_pts2(2,:) = 0.5_f64*es%delta_eta2*es%gauss_pts2(2,:)

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
   spline_degree1,      &
   num_cells1,          &
   eta1_min,            &
   eta1_max,            &
   bc_left_knots,       &
   bc_right,            &
   es%knots1 )

 call initialize_knots( &
   spline_degree2,      &
   num_cells2,          &
   eta2_min,            &
   eta2_max,            &
   bc_bottom,           &
   bc_top,              &
   es%knots2 )

 call initconnectivity(               &
   num_cells1,                        &
   num_cells2,                        &
   spline_degree1,                    &
   spline_degree2,                    &
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
   num_cells1*num_cells2,             &
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

 ! allocation of the table containning all values of splines and its
 ! derivatives in each gauss points
 SLL_ALLOCATE(es%v_splines1(spline_degree1+1,spline_degree1+2,num_cells1),ierr)
 SLL_ALLOCATE(es%d_splines1(spline_degree1+1,spline_degree1+2,num_cells1),ierr)
 SLL_ALLOCATE(es%r_splines1(spline_degree1+1,spline_degree1+2,num_cells1),ierr)
 SLL_ALLOCATE(es%v_splines2(spline_degree2+1,spline_degree2+2,num_cells2),ierr)
 SLL_ALLOCATE(es%d_splines2(spline_degree2+1,spline_degree2+2,num_cells2),ierr)
 SLL_ALLOCATE(es%r_splines2(spline_degree2+1,spline_degree2+2,num_cells2),ierr)

 es%v_splines1 = 0.0_f64
 es%v_splines1 = 0.0_f64
 es%d_splines2 = 0.0_f64
 es%d_splines2 = 0.0_f64
 es%r_splines2 = 0.0_f64
 es%r_splines2 = 0.0_f64

 SLL_ALLOCATE(es%val_jac(spline_degree1+2,spline_degree2+2,num_cells1,num_cells2),ierr)
 es%val_jac = 0.0_f64
 SLL_ALLOCATE(es%tab_index_coeff1(num_cells1),ierr)
 SLL_ALLOCATE(es%tab_index_coeff2(num_cells2),ierr)

 num_pts_g1 = size(es%gauss_pts1,2)
 num_pts_g2 = size(es%gauss_pts2,2)
 SLL_ALLOCATE(es%rho_at_gauss(num_pts_g1,num_pts_g2,num_cells1,num_cells2),ierr)
 SLL_ALLOCATE(es%rho_coeff_1d((num_cells1+1)*(num_cells2+1)),ierr)
 SLL_ALLOCATE(es%M_rho_loc(es%total_num_splines_loc),ierr)

allocate(work1(spline_degree1+1,spline_degree1+1))
allocate(work2(spline_degree2+1,spline_degree2+1))
allocate(dbs1(spline_degree1+1,2))
allocate(dbs2(spline_degree2+1,2))

do i = 1, es%num_cells1
  eta1  = eta1_min + (i-1)*es%delta_eta1
  do ii=1,num_pts_g1
    xg  = eta1  + es%gauss_pts1(1,ii)
    wxg = es%gauss_pts1(2,ii)
    if(bc_left==SLL_PERIODIC .and. bc_right==SLL_PERIODIC) then 
      gspl1 = es%gauss_pts1(1,ii)
      ispl1 = spline_degree1+1
    else 
      gspl1 = xg
      ispl1 = spline_degree1+i
    end if
    call bsplvd( es, es%knots1, spline_degree1+1, gspl1, ispl1, work1, dbs1, 2 )
    es%v_splines1(:,ii,i) = dbs1(:,1)
    es%d_splines1(:,ii,i) = dbs1(:,2)
    call interv( es, es%knots1_rho, es%num_cells1+spline_degree1+2, xg, left_1, ierr )
    call bsplvd( es, es%knots1_rho, spline_degree1+1, xg,left_1, work1, dbs1, 1 )
    es%r_splines1(:,ii,i) = dbs1(:,1)
  end do
  es%tab_index_coeff1(i) = left_1
end do

do j = 1, es%num_cells2
  eta2  = eta2_min + (j-1)*es%delta_eta2
  do jj=1,num_pts_g2
    yg  = eta2+es%gauss_pts2(1,jj)
    wyg = es%gauss_pts2(2,jj)
    if (bc_bottom==SLL_PERIODIC .and. bc_top==SLL_PERIODIC) then
      gspl2 = es%gauss_pts2(1,jj)
      ispl2 = spline_degree2+1
    else
      gspl2 = yg
      ispl2 = spline_degree2+j
    end if
    call bsplvd( es, es%knots2, spline_degree2+1, gspl2, ispl2, work2, dbs2, 2)
    es%v_splines2(:,jj,j) = dbs2(:,1)
    es%d_splines2(:,jj,j) = dbs2(:,2)
    call interv( es, es%knots2_rho, es%num_cells2+spline_degree2+2, yg, left_2, ierr )
    call bsplvd( es, es%knots2_rho, spline_degree2+1, yg, left_2, work2, dbs2, 1)
    es%r_splines2(:,jj,j) = dbs2(:,1)
  end do
  es%tab_index_coeff2(j) = left_2
end do


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

function new_general_elliptic_solver( spline_degree1,   &
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
  SLL_DEALLOCATE(es%v_splines1,ierr)
  SLL_DEALLOCATE(es%v_splines2,ierr)
  SLL_DEALLOCATE(es%d_splines1,ierr)
  SLL_DEALLOCATE(es%d_splines2,ierr)
  SLL_DEALLOCATE(es%val_jac,ierr)
  SLL_DEALLOCATE(es%tab_index_coeff1,ierr)
  SLL_DEALLOCATE(es%tab_index_coeff2,ierr)

end subroutine delete_elliptic


!> @brief Assemble the matrix for elliptic solver.
!> @details To have the function phi such that 
!> \f[
!>  \nabla \cdot ( A \nabla \phi ) + B \nabla \phi + C \phi = \rho
!> \f]
!>  where 
!>  - A is a matrix of functions , 
!>  - B a vectorial function,
!>  - C  a scalar function.  
!>  - rho a scalar function.  
!>  - A, B, C, rho can be discret or analytic. 
!>  - \f$ \phi \f$ is given with a B-spline interpolator  
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

sll_real64, dimension(:,:),   allocatable :: M_c
sll_real64, dimension(:,:),   allocatable :: K_11
sll_real64, dimension(:,:),   allocatable :: K_12
sll_real64, dimension(:,:),   allocatable :: K_21
sll_real64, dimension(:,:),   allocatable :: K_22
sll_real64, dimension(:,:),   allocatable :: M_bv
sll_real64, dimension(:,:),   allocatable :: S_b1
sll_real64, dimension(:,:),   allocatable :: S_b2  
sll_real64, dimension(:),     allocatable :: mass
sll_real64, dimension(:),     allocatable :: stiff
sll_real64, dimension(:,:,:), pointer     :: source

sll_int32 :: ierr
sll_int32 :: i
sll_int32 :: j
sll_int32 :: icell
sll_int32 :: bc_left
sll_int32 :: bc_right
sll_int32 :: bc_bottom
sll_int32 :: bc_top

character(len=*),parameter :: as_file1='mat'

sll_int32 :: cell_i
sll_int32 :: cell_j

sll_real64 :: delta1
sll_real64 :: delta2
sll_real64 :: eta1_min
sll_real64 :: eta2_min
sll_real64 :: eta1
sll_real64 :: eta2
sll_int32  :: num_pts_g1 ! number of gauss points in first direction 
sll_int32  :: num_pts_g2 ! number of gauss points in second direction
sll_int32  :: ii,kk
sll_int32  :: jj,ll
sll_real64 :: xg, wxg
sll_real64 :: yg, wyg
sll_real64 :: gspl1
sll_real64 :: gspl2
sll_int32  :: ispl1
sll_int32  :: ispl2
sll_int32  :: index1
sll_int32  :: index2
sll_real64, dimension(es%spline_degree1+1,es%spline_degree1+1) :: work1
sll_real64, dimension(es%spline_degree2+1,es%spline_degree2+1) :: work2
sll_real64, dimension(es%spline_degree1+1,2) :: dbs1
sll_real64, dimension(es%spline_degree2+1,2) :: dbs2
sll_real64, dimension(es%spline_degree1+1,1) :: dbs1_rho
sll_real64, dimension(es%spline_degree2+1,1) :: dbs2_rho
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
sll_real64 :: jac_mat(2,2)
sll_real64 :: val_jac
sll_real64 :: B11
sll_real64 :: B12
sll_real64 :: B21
sll_real64 :: B22
sll_real64 :: MC
sll_real64 :: C1
sll_real64 :: C2    
sll_int32 :: left_x,left_y
sll_int32 :: index3, index4
sll_int32 :: index_coef1,index_coef2
sll_int32 :: mm, nn, b, bprime,x,y
sll_int32 :: a, aprime
sll_real64 :: elt_mat_global
sll_int32  :: nspl, nbsp,nbsp1
sll_int32 :: spl_deg_1, spl_deg_2, nc_1, nc_2
sll_int32 :: ideg2,ideg1
sll_int32 :: jdeg2,jdeg1

bc_left    = es%bc_left
bc_right   = es%bc_right
bc_bottom  = es%bc_bottom
bc_top     = es%bc_top
delta1     = es%delta_eta1 
delta2     = es%delta_eta2 
eta1_min   = es%eta1_min  
eta2_min   = es%eta2_min  
num_pts_g1 = size(es%gauss_pts1,2)
num_pts_g2 = size(es%gauss_pts2,2)
spl_deg_1  = es%spline_degree1
spl_deg_2  = es%spline_degree2
nc_1       = es%num_cells1
nc_2       = es%num_cells2

nspl = es%total_num_splines_loc

if( bc_left   == SLL_PERIODIC .and. bc_right == SLL_PERIODIC .and. &
    bc_bottom == SLL_PERIODIC .and. bc_top   == SLL_PERIODIC ) then
   es%perper = .true.
else
   es%perper = .false.  
end if   

SLL_CLEAR_ALLOCATE(source(1:nc_1*nc_2,1:nspl,1:nspl),ierr)
SLL_CLEAR_ALLOCATE(M_c(1:nspl,1:nspl),ierr)
SLL_CLEAR_ALLOCATE(K_11(1:nspl,1:nspl),ierr)
SLL_CLEAR_ALLOCATE(K_12(1:nspl,1:nspl),ierr)
SLL_CLEAR_ALLOCATE(K_21(1:nspl,1:nspl),ierr)
SLL_CLEAR_ALLOCATE(K_22(1:nspl,1:nspl),ierr)
SLL_CLEAR_ALLOCATE(S_b1(1:nspl,1:nspl),ierr)
SLL_CLEAR_ALLOCATE(S_b2(1:nspl,1:nspl),ierr)
SLL_CLEAR_ALLOCATE(M_bv(1:nspl,1:nspl),ierr)
SLL_CLEAR_ALLOCATE(mass(1:nspl),ierr)
SLL_CLEAR_ALLOCATE(stiff(1:nspl),ierr)

do cell_j = 1, nc_2
do cell_i = 1, nc_1
        
  eta1  = eta1_min + (cell_i-1)*delta1
  eta2  = eta2_min + (cell_j-1)*delta2
  icell = cell_i+nc_1*(cell_j-1)
    
  mass  = 0.0_f64
  stiff = 0.0_f64
  M_c   = 0.0_f64
  K_11  = 0.0_f64
  K_12  = 0.0_f64
  K_21  = 0.0_f64
  K_22  = 0.0_f64
  M_bv  = 0.0_f64
  S_b1  = 0.0_f64
  S_b2  = 0.0_f64
  dbs1  = 0.0_f64
  dbs2  = 0.0_f64
  work1 = 0.0_f64
  work2 = 0.0_f64

  do j=1,num_pts_g2
  
    yg  = eta2+es%gauss_pts2(1,j)
    wyg = es%gauss_pts2(2,j)
   !  
   ! if (es%bc_bottom == SLL_PERIODIC .and. &
   !     es%bc_top    == SLL_PERIODIC) then
   !
   !   gspl2 = es%gauss_pts2(1,j)
   !   ispl2 = spl_deg_2 + 1
   !     
   ! else if (es%bc_bottom == SLL_DIRICHLET .and.&
   !          es%bc_top    == SLL_DIRICHLET) then
   !   gspl2 = yg
   !   ispl2 = spl_deg_2 + cell_j
   ! end if
     
    !call bsplvd( es, es%knots2, spl_deg_2+1, gspl2, ispl2, work2, dbs2, 2)
    !call interv( es, es%knots2_rho, nc_2+spl_deg_2+2, yg, left_y, ierr )
    !call bsplvd( es, es%knots2_rho, spl_deg_2+1, yg, left_y, work2, dbs2_rho, 1)
    
    !es%v_splines2(:,j,cell_j) = dbs2(:,1)
  
    do i=1,num_pts_g1
    
      xg  = eta1  + es%gauss_pts1(1,i)
      wxg = es%gauss_pts1(2,i)

        
    !  if(es%bc_left_interp == SLL_PERIODIC .and. &
    !     es%bc_right       == SLL_PERIODIC) then 
    !         
    !    gspl1 = es%gauss_pts1(1,i)
    !    ispl1 = spl_deg_1+1
    !         
    !  else if (es%bc_left_interp  == SLL_DIRICHLET .and.&
    !           es%bc_right        == SLL_DIRICHLET) then
    !       
    !    gspl1 = xg
    !    ispl1 = spl_deg_1 + cell_i
    !         
    !  end if
  
      !call bsplvd( es, es%knots1, spl_deg_1+1, gspl1,ispl1, work1, dbs1, 2 )
      !call interv( es, es%knots1_rho, nc_1+spl_deg_1+2, xg, left_x, ierr )
      !call bsplvd( es, es%knots1_rho, spl_deg_1+1, xg,left_x, work1, dbs1_rho, 1 )

      !es%v_splines1(:,i,cell_i)       = dbs1(:,1)

      val_c       = c_field%value_at_point(xg,yg)
      val_a11     = a11_field_mat%value_at_point(xg,yg)
      val_a12     = a12_field_mat%value_at_point(xg,yg)
      val_a21     = a21_field_mat%value_at_point(xg,yg)
      val_a22     = a22_field_mat%value_at_point(xg,yg)

      val_b1      = b1_field_vect%value_at_point(xg,yg)
      val_b1_der1 = b1_field_vect%first_deriv_eta1_value_at_point(xg,yg)
      val_b1_der2 = b1_field_vect%first_deriv_eta2_value_at_point(xg,yg)

      val_b2      = b2_field_vect%value_at_point(xg,yg)
      val_b2_der1 = b2_field_vect%first_deriv_eta1_value_at_point(xg,yg)
      val_b2_der2 = b2_field_vect%first_deriv_eta2_value_at_point(xg,yg)
 
      jac_mat = c_field%get_jacobian_matrix(xg,yg)
      val_jac = jac_mat(1,1)*jac_mat(2,2) - jac_mat(1,2)*jac_mat(2,1)

      es%val_jac(i,j,cell_i,cell_j) = val_jac
        
      es%intjac = es%intjac + wyg*wxg*val_jac

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
          
      C1 =  jac_mat(2,2) * val_b1 - jac_mat(1,2) * val_b2 
      C2 =  jac_mat(1,1) * val_b2 - jac_mat(2,1) * val_b1
         
      ! loop over the splines supported in the cell that are different than
      ! zero at the point (xg,yg) (there are spline_degree+1 splines in
      ! each direction.
      do ii = 0,spl_deg_1
      do jj = 0,spl_deg_2
              
        mm = jj*(spl_deg_1+1)+ii+1
              
        mass(mm) = mass(mm)+val_jac*wxg*wyg* &
          (es%v_splines1(ii+1,i,cell_i)*es%v_splines2(jj+1,j,cell_j))

        stiff(mm) = stiff(mm)+val_jac*wxg*wyg* &
          (es%d_splines1(ii+1,i,cell_i)*es%v_splines2(jj+1,j,cell_j)+&
           es%v_splines1(ii+1,i,cell_i)*es%d_splines2(jj+1,j,cell_j))
             
        do kk = 0,spl_deg_1
        do ll = 0,spl_deg_2
                    
          nn =  ll*(spl_deg_1+1)+kk+1
         
          source(icell,mm,nn) = source(icell,mm,nn) + val_jac*wxg*wyg * &
            es%r_splines1(ii+1,i,cell_i)*es%v_splines1(kk+1,i,cell_i)*  &
            es%r_splines2(jj+1,j,cell_j)*es%v_splines2(ll+1,j,cell_j)
               
          M_c(mm,nn) = M_c(mm,nn) + val_c*val_jac*wxg*wyg*              &
            es%v_splines1(ii+1,i,cell_i)*es%v_splines1(kk+1,i,cell_i) * &
            es%v_splines2(jj+1,j,cell_j)*es%v_splines2(ll+1,j,cell_j)
               
          K_11(mm,nn) = K_11(mm,nn) + B11*wxg*wyg/val_jac*              &
            es%d_splines1(ii+1,i,cell_i)*es%d_splines1(kk+1,i,cell_i) * &
            es%v_splines2(jj+1,j,cell_j)*es%v_splines2(ll+1,j,cell_j)
              
          K_22(mm,nn) = K_22(mm,nn) + B22*wxg*wyg/val_jac*              &
            es%v_splines1(ii+1,i,cell_i)*es%v_splines1(kk+1,i,cell_i) * &
            es%d_splines2(jj+1,j,cell_j)*es%d_splines2(ll+1,j,cell_j)
               
          K_12(mm,nn) = K_12(mm,nn) + B12*wxg*wyg/val_jac*              &
            es%d_splines1(ii+1,i,cell_i)*es%v_splines1(kk+1,i,cell_i) * &
            es%v_splines2(jj+1,j,cell_j)*es%d_splines2(ll+1,j,cell_j)
               
          K_21(mm,nn) = K_21(mm,nn) + B21*wxg*wyg/val_jac*              &
            es%v_splines1(ii+1,i,cell_i)*es%d_splines1(kk+1,i,cell_i) * &
            es%d_splines2(jj+1,j,cell_j)*es%v_splines2(ll+1,j,cell_j)

          M_bv(mm,nn) = M_bv(mm,nn) + MC*wxg*wyg *                      &
            es%v_splines1(ii+1,i,cell_i)*es%v_splines1(kk+1,i,cell_i) * &
            es%v_splines2(jj+1,j,cell_j)*es%v_splines2(ll+1,j,cell_j)
               
          ! A revoir 
          S_b1(mm,nn) = S_b1(mm,nn) + C1*wxg*wyg *                      &
            es%v_splines1(ii+1,i,cell_i)*es%d_splines1(kk+1,i,cell_i) * &
            es%v_splines2(jj+1,j,cell_j)*es%v_splines2(ll+1,j,cell_j)

          ! A revoir 
          S_b2(mm,nn) = S_b2(mm,nn) + C2*wxg*wyg *                      &
            es%v_splines1(ii+1,i,cell_i)*es%v_splines1(kk+1,i,cell_i) * &
            es%v_splines2(jj+1,j,cell_j)*es%d_splines2(ll+1,j,cell_j)

        end do
        end do
      end do
      end do

    end do
  end do

  do j = 0, spl_deg_2

    index3 = cell_j + j
    
    if (bc_bottom==SLL_PERIODIC .and. bc_top==SLL_PERIODIC) then 
      if ( index3 > es%total_num_splines2) then
        index3 = index3 - es%total_num_splines2
      end if
    end if
     
    do i = 0,spl_deg_1
        
      index1 = cell_i + i

      if (bc_left==SLL_PERIODIC .and. bc_right==SLL_PERIODIC) then 
        if ( index1 > es%total_num_splines1) then
          index1 = index1 - es%total_num_splines1
        end if
        nbsp = es%total_num_splines1
      else if (bc_left==SLL_DIRICHLET .and. bc_right==SLL_DIRICHLET) then
        nbsp = nc_1 + spl_deg_1
      end if

      x = index1 + (index3-1)*nbsp
      b = j * ( spl_deg_1 + 1 ) + i + 1
      a = es%local_to_global_spline_indices(b, icell)
         
      es%masse(x) = es%masse(x) + mass(b)
      es%stiff(x) = es%stiff(x) + stiff(b)

      index_coef1 = es%tab_index_coeff1(cell_i) - spl_deg_1 + i
      index_coef2 = es%tab_index_coeff2(cell_j) - spl_deg_2 + j

      es%local_to_global_spline_indices_source(b,icell)= &
              index_coef1 + (index_coef2-1)*(nc_1+1)

      do nn = 0,spl_deg_2
             
        index4 = cell_j + nn
             
        if ( (bc_bottom==SLL_PERIODIC).and.(bc_top== SLL_PERIODIC))then
          if ( index4 > es%total_num_splines2) then
            index4 = index4 - es%total_num_splines2
          end if
        end if
             
        do mm = 0,spl_deg_1
                
          index2 = cell_i + mm

          if((bc_left==SLL_PERIODIC).and.(bc_right==SLL_PERIODIC))then

            if ( index2 > es%total_num_splines1) then
              index2 = index2 - es%total_num_splines1
            end if

            nbsp1 = es%total_num_splines1
                   
          else if (bc_left ==SLL_DIRICHLET .and.&
                   bc_right==SLL_DIRICHLET ) then

            nbsp1 = nc_1 + spl_deg_1

          end if
                
          y      = index2 + (index4-1)*nbsp1
          bprime =  nn * ( spl_deg_1 + 1 ) + mm + 1
          aprime = es%local_to_global_spline_indices(bprime,icell)

          elt_mat_global = M_c(b,  bprime) - &
                           K_11(b, bprime) - &
                           K_12(b, bprime) - &
                           K_21(b, bprime) - &
                           K_22(b, bprime) - &
                           M_bv(b, bprime) - &
                           S_b1(b, bprime) - &
                           S_b2(b, bprime)

          es%local_to_global_spline_indices_source_bis(bprime,icell)= y

          if ( a>0 .and. aprime>0 ) then
            call sll_add_to_csr_matrix(es%sll_csr_mat, elt_mat_global, a, aprime)   
          end if
              
        end do
      end do
    end do
  end do
end do
end do

print *,'#begin of sll_factorize_csr_matrix'
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

  call sll_factorize_csr_matrix(es%sll_csr_mat_with_constraint)
else   
  call sll_factorize_csr_matrix(es%sll_csr_mat)
end if 

print *,'#end of sll_factorize_csr_matrix'

es%sll_csr_mat_source => new_csr_matrix( size(es%masse,1),                 &
                      (nc_1+1)*(nc_2+1), nc_1*nc_2,                        &
                      es%local_to_global_spline_indices_source_bis,        &
                      nspl, es%local_to_global_spline_indices_source, nspl )

do cell_j=1,es%num_cells2
do cell_i=1,es%num_cells1
      
  icell = cell_i+es%num_cells1*(cell_j-1)
      
  do ideg2 = 0,es%spline_degree2
  do ideg1 = 0,es%spline_degree1
            
    b = ideg2 * ( es%spline_degree1 + 1 ) + ideg1 + 1
    a = es%local_to_global_spline_indices_source_bis(b, icell)
        
    do jdeg2 = 0,es%spline_degree2
    do jdeg1 = 0,es%spline_degree1
              
      bprime = jdeg2 * ( es%spline_degree1 + 1 ) + jdeg1 + 1
      aprime = es%local_to_global_spline_indices_source(bprime,icell)
           
      elt_mat_global = source(icell,bprime,b)

      if ( a > 0 .and. aprime > 0) then
        call sll_add_to_csr_matrix(es%sll_csr_mat_source,elt_mat_global,a,aprime)
      end if
              
    end do
    end do

  end do
  end do

end do
end do

SLL_DEALLOCATE_ARRAY(source,ierr)
SLL_DEALLOCATE_ARRAY(M_c,ierr)
SLL_DEALLOCATE_ARRAY(K_11,ierr)
SLL_DEALLOCATE_ARRAY(K_12,ierr)
SLL_DEALLOCATE_ARRAY(K_21,ierr)
SLL_DEALLOCATE_ARRAY(K_22,ierr)
SLL_DEALLOCATE_ARRAY(M_bv,ierr)
SLL_DEALLOCATE_ARRAY(S_b1,ierr)
SLL_DEALLOCATE_ARRAY(S_b2,ierr)
SLL_DEALLOCATE_ARRAY(stiff,ierr) 
SLL_DEALLOCATE_ARRAY(mass,ierr) 
   
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
  sll_int32 :: total_num_splines_loc
  sll_real64 :: int_rho,int_jac
  sll_int32 :: num_pts_g1, num_pts_g2, ig, jg
  sll_real64 :: wgpt1, wgpt2, gpt1, gpt2, eta1, eta2
  sll_real64 :: val_jac
  class(sll_scalar_field_2d_base),pointer  :: base_field_pointer
  class(sll_interpolator_2d_base),pointer  :: base_interpolator_pointer
  sll_real64, dimension(:,:), pointer :: coeff_rho
    
  total_num_splines_loc = es%total_num_splines_loc
    
  num_pts_g1 = size(es%gauss_pts1,2)
  num_pts_g2 = size(es%gauss_pts2,2)
  es%rho_at_gauss = 0.0_f64   

  es%M_rho_loc    = 0.0_f64
  es%rho_vec      = 0.0_f64
  es%rho_coeff_1d = 0.0_f64
   
  !ES Compute rho at all Gauss points
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
                
      do j=1,es%num_cells2+1
        do i=1,es%num_cells1+1
          es%rho_coeff_1d(i+(es%num_cells1+1)*(j-1)) = coeff_rho(i,j)
        end do
      end do

      call sll_mult_csr_matrix_vector(es%sll_csr_mat_source,es%rho_coeff_1d,es%rho_vec)

      if(es%perper) then
        es%rho_vec = es%rho_vec - sum(es%rho_vec)/es%intjac*es%masse
      end if
        
    class default
          
      do j=1,es%num_cells2
        eta2  = es%eta2_min + (j-1)*es%delta_eta2
        do i=1,es%num_cells1
          eta1  = es%eta1_min + (i-1)*es%delta_eta1
          do jg=1,num_pts_g2
            gpt2  = eta2  + es%gauss_pts2(1,jg)
            wgpt2 = es%gauss_pts2(2,jg)
            do ig=1,num_pts_g1
              gpt1  = eta1  + es%gauss_pts1(1,ig)
              wgpt1 = es%gauss_pts1(2,ig)
              es%rho_at_gauss(ig,jg,i,j)   = rho%value_at_point(gpt1,gpt2)
              val_jac = es%val_jac(ig,jg,i,j)
              int_rho = int_rho + rho%value_at_point(gpt1,gpt2)*wgpt2*wgpt1*val_jac 
              int_jac = int_jac + wgpt2*wgpt1*val_jac
            end do
          end do
        end do
      end do

      if (es%perper) es%rho_at_gauss = es%rho_at_gauss - int_rho/int_jac
          
      do j=1,es%num_cells2
      do i=1,es%num_cells1
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
          gpt2  = eta2  + es%gauss_pts2(1,jg)
          wgpt2 = es%gauss_pts2(2,jg)
          do ig=1,num_pts_g1
            gpt1  = eta1  + es%gauss_pts1(1,ig)
            wgpt1 = es%gauss_pts1(2,ig)
            es%rho_at_gauss(ig,jg,i,j) = rho%value_at_point(gpt1,gpt2)
            val_jac = es%val_jac(ig,jg,i,j)
            int_rho = int_rho + rho%value_at_point(gpt1,gpt2)*wgpt2*wgpt1*val_jac 
            int_jac = int_jac + wgpt2*wgpt1*val_jac
          end do
        end do
      end do
    end do

    if (es%perper) es%rho_at_gauss = es%rho_at_gauss - int_rho/int_jac
          
    do j=1,es%num_cells2
    do i=1,es%num_cells1
      call build_local_matrices_rho(es, i, j)
      call local_to_global_matrices_rho(es, i, j)
    end do
    end do
       
  end select

  call solve_linear_system(es)
  
  print *,'#solve_linear_system done'
    
  call phi%interp_2d%set_coefficients(es%phi_vec)

end subroutine solve_general_coordinates_elliptic_eq
  
subroutine build_local_matrices_rho( es, cell_i, cell_j)

  class(general_coordinate_elliptic_solver) :: es
  sll_int32, intent(in) :: cell_i
  sll_int32, intent(in) :: cell_j
  sll_int32 :: bc_left    
  sll_int32 :: bc_right
  sll_int32 :: bc_bottom    
  sll_int32 :: bc_top  
  sll_real64 :: eta1
  sll_real64 :: eta2
  sll_int32  :: num_pts_g1 
  sll_int32  :: num_pts_g2
  sll_int32  :: i,ii
  sll_int32  :: j,jj
  sll_real64 :: gpt1
  sll_real64 :: gpt2
  sll_real64 :: wgpt1
  sll_real64 :: wgpt2
  sll_int32  :: index1
  sll_real64, dimension(es%spline_degree1+1,2) :: dbiatx1
  sll_real64, dimension(es%spline_degree2+1,2) :: dbiatx2
  sll_real64 :: val_f
  sll_real64 :: val_jac,spline1,spline2
    
  es%M_rho_loc(:)  = 0.0_f64
  dbiatx1(:,:)  = 0.0_f64
  dbiatx2(:,:)  = 0.0_f64
  ! The supposition is that all fields use the same logical mesh
  bc_left   = es%bc_left
  bc_right  = es%bc_right
  bc_bottom = es%bc_bottom
  bc_top    = es%bc_top

  num_pts_g1 = size(es%gauss_pts1,2)
  num_pts_g2 = size(es%gauss_pts2,2)
    
  eta1  = es%eta1_min + (cell_i-1)*es%delta_eta1
  eta2  = es%eta2_min + (cell_j-1)*es%delta_eta2
    
  do j=1,num_pts_g2
    gpt2  = eta2 + es%gauss_pts2(1,j)
    wgpt2 = es%gauss_pts2(2,j)

    do i=1,num_pts_g1
      gpt1  = eta1 + es%gauss_pts1(1,i)
      wgpt1 = es%gauss_pts1(2,i)
  
      val_f = es%rho_at_gauss(i,j,cell_i,cell_j)
      val_jac = es%val_jac(i,j,cell_i,cell_j)

     ! loop over the splines supported in the cell that are different than
     ! zero at the point (gpt1,gpt2) (there are spline_degree+1 splines in
     ! each direction.
      do ii = 0,es%spline_degree1
        do jj = 0,es%spline_degree2
                
          spline1 = es%v_splines1(ii+1,i,cell_i)
          spline2 = es%v_splines2(jj+1,j,cell_j)
               
          index1  =  jj * ( es%spline_degree1 + 1 ) + ii + 1
          es%M_rho_loc(index1)= es%M_rho_loc(index1) + &
                    val_f*val_jac*wgpt1*wgpt2*spline1*spline2
                
        end do
      end do
    end do
  end do
    
end subroutine build_local_matrices_rho
  

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
    !other option for above: index3 = mod(index3-1,es%total_num_splines2)+1
    if (bc_bottom==SLL_PERIODIC .and. bc_top==SLL_PERIODIC) then    
       if ( index3 > es%total_num_splines2) then
          index3 = index3 - es%total_num_splines2
       end if
    end if

    do i = 0,es%spline_degree1
      index1 = cell_i + i
      if (bc_left==SLL_PERIODIC .and. bc_right==SLL_PERIODIC) then 
        if ( index1 > es%total_num_splines1) then
          index1 = index1 - es%total_num_splines1
        end if
        nbsp = es%total_num_splines1
      else if (bc_left==SLL_DIRICHLET .and. bc_right==SLL_DIRICHLET) then
        nbsp = es%num_cells1 + es%spline_degree1
      end if

      x             =  index1 + (index3-1)*nbsp
      b             =  mm * ( es%spline_degree1 + 1 ) + i + 1
      es%rho_vec(x) =  es%rho_vec(x)  + es%M_rho_loc(b)
        
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
    
  if( bc_left  ==SLL_PERIODIC  .and. bc_right==SLL_PERIODIC .and. &
      bc_bottom==SLL_DIRICHLET .and. bc_top  ==SLL_DIRICHLET ) then
       
    do i = 1, es%total_num_splines1
      do j = 1, es%total_num_splines2
             
        elt  = i + es%total_num_splines1 * (  j - 1)
        elt1 = i + ( es%total_num_splines1 ) * j
        es%tmp_rho_vec(elt) = es%rho_vec(elt1)

      end do
    end do
       
  else if( bc_left  ==SLL_DIRICHLET .and. bc_right==SLL_DIRICHLET .and.&
           bc_bottom==SLL_DIRICHLET .and. bc_top  ==SLL_DIRICHLET) then 
       
    do i = 1, es%total_num_splines1
      do j = 1, es%total_num_splines2
            
        elt  = i + es%total_num_splines1 * (  j - 1)
        elt1 = i + 1 + ( es%total_num_splines1 + 2 ) * j 
        es%tmp_rho_vec( elt ) = es%rho_vec( elt1 )

      end do
    end do
       
  else if(bc_left  ==SLL_PERIODIC .and. bc_right==SLL_PERIODIC .and.&
          bc_bottom==SLL_PERIODIC .and. bc_top  ==SLL_PERIODIC) then
       
    es%tmp_rho_vec(1:es%total_num_splines1*es%total_num_splines2)=&
      es%rho_vec(1:es%total_num_splines1*es%total_num_splines2) 
       
  else if(bc_left  ==SLL_DIRICHLET .and. bc_right==SLL_DIRICHLET .and.&
          bc_bottom==SLL_PERIODIC  .and. bc_top  ==SLL_PERIODIC ) then
     
    do i = 1, es%total_num_splines1
    do j = 1, es%total_num_splines2

        elt1 = i+1 + (es%total_num_splines1+2)*(j-1)
        elt  = i + es%total_num_splines1*(j-1)
        es%tmp_rho_vec(elt) = es%rho_vec(elt1)

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

  if(bc_left  ==SLL_PERIODIC .and. bc_right==SLL_PERIODIC .and.&
     bc_bottom==SLL_PERIODIC .and. bc_top  ==SLL_PERIODIC) then

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
  

!*************************************************************************
!
!! INTERV brackets a real value in an ascending vector of values.
!
!  Discussion:
!
!    The XT array is a set of increasing values.  The goal of the routine
!    is to determine the largest index I so that XT(I) <= X.
!
!    The routine is designed to be efficient in the common situation
!    that it is called repeatedly, with X taken from an increasing
!    or decreasing sequence.
!
!    This will happen when a piecewise polynomial is to be graphed.
!    The first guess for LEFT is therefore taken to be the value
!    returned at the previous call and stored in the local variable ILO.
!
!    A first check ascertains that ILO < LXT.  This is necessary
!    since the present call may have nothing to do with the previous
!    call.  Then, if
!
!      XT(ILO) <= X < XT(ILO+1),
!
!    we set LEFT = ILO and are done after just three comparisons.
!
!    Otherwise, we repeatedly double the difference ISTEP = IHI - ILO
!    while also moving ILO and IHI in the direction of X, until
!
!      XT(ILO) <= X < XT(IHI)
!
!    after which we use bisection to get, in addition, ILO + 1 = IHI.
!    The value LEFT = ILO is then returned.
!
!  Modified:
!
!    14 February 2007
!
!  Author:
!
!    Carl DeBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XT(LXT), a nondecreasing sequence of values.
!
!    Input, integer LXT, the dimension of XT.
!
!    Input, real ( kind = 8 ) X, the point whose location with
!    respect to the sequence XT is to be determined.
!
!    Output, integer LEFT, the index of the bracketing value:
!      1     if             X  <  XT(1)
!      I     if   XT(I)  <= X  < XT(I+1)
!      LXT   if  XT(LXT) <= X
!
!    Output, integer MFLAG, indicates whether X lies within the
!    range of the data.
!    -1:            X  <  XT(1)
!     0: XT(I)   <= X  < XT(I+1)
!    +1: XT(LXT) <= X
!

subroutine interv( es, xt, lxt, x, left, mflag )
    
    type(general_coordinate_elliptic_solver) :: es

    sll_int32,intent(in)  :: lxt
    sll_int32,intent(out) :: left
    sll_int32,intent(out) :: mflag
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

!*************************************************************************
!
!! BSPLVD calculates the nonvanishing B-splines and derivatives at X.
!
!  Discussion:
!
!    Values at X of all the relevant B-splines of order K:K+1-NDERIV
!    are generated via BSPLVB and stored temporarily in DBIATX.
!
!    Then the B-spline coefficients of the required derivatives
!    of the B-splines of interest are generated by differencing,
!    each from the preceding one of lower order, and combined with
!    the values of B-splines of corresponding order in DBIATX
!    to produce the desired values.
!
!  Modified:
!
!    14 February 2007
!
!  Author:
!
!    Carl DeBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663.
!
!  Parameters:
!
!Input, real ( kind = 8 ) T(LEFT+K), the knot sequence.  It is assumed that
!    T(LEFT) < T(LEFT+1).  Also, the output is correct only if
!    T(LEFT) <= X <= T(LEFT+1).
!
!    Input, integer K, the order of the B-splines to be evaluated.
!
!    Input, real ( kind = 8 ) X, the point at which these values are sought.
!
!    Input, integer LEFT, indicates the left endpoint of the interval of
!    interest.  The K B-splines whose support contains the interval
!    ( T(LEFT), T(LEFT+1) ) are to be considered.
!
!    Workspace, real ( kind = 8 ) A(K,K).
!
!    Output, real ( kind = 8 ) DBIATX(K,NDERIV).  DBIATX(I,M) contains
!    the value of the (M-1)st derivative of the (LEFT-K+I)-th B-spline
!    of order K for knot sequence T, I=M,...,K, M=1,...,NDERIV.
!
!    Input, integer NDERIV, indicates that values of B-splines and their
!    derivatives up to but not including the NDERIV-th are asked for.
!

 subroutine bsplvd ( es, t, k, x, left, a, dbiatx, nderiv )

    type(general_coordinate_elliptic_solver) :: es

    sll_int32  :: k
    sll_int32  :: left
    sll_int32  :: nderiv
    
    sll_real64 :: a(:,:)
    sll_real64,dimension(:,:), intent(out) :: dbiatx!(k,nderiv)
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


!***********************************************************************
!
!! BSPLVB evaluates B-splines at a point X with a given knot sequence.
!
!  Discusion:
!
!    BSPLVB evaluates all possibly nonzero B-splines at X of order
!
!      JOUT = MAX ( JHIGH, (J+1)*(INDEX-1) )
!
!    with knot sequence T.
!
!    The recurrence relation
!
!                     X - T(I)               T(I+J+1) - X
!    B(I,J+1)(X) = ----------- * B(I,J)(X) + --------------- * B(I+1,J)(X)
!                  T(I+J)-T(I)               T(I+J+1)-T(I+1)
!
!    is used to generate B(LEFT-J:LEFT,J+1)(X) from B(LEFT-J+1:LEFT,J)(X)
!    storing the new values in BIATX over the old.
!
!    The facts that
!
!      B(I,1)(X) = 1  if  T(I) <= X < T(I+1)
!
!    and that
!
!      B(I,J)(X) = 0  unless  T(I) <= X < T(I+J)
!
!    are used.
!
!    The particular organization of the calculations follows
!    algorithm 8 in chapter X of the text.
!
!  Modified:
!
!    14 February 2007
!
!  Author:
!
!    Carl DeBoor
!
!  Reference:
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663.
!
!  Parameters:
!
!Input, real ( kind = 8 ) T(LEFT+JOUT), the knot sequence.  T is assumed to
!    be nondecreasing, and also, T(LEFT) must be strictly less than
!    T(LEFT+1).
!
!    Input, integer JHIGH, INDEX, determine the order
!    JOUT = max ( JHIGH, (J+1)*(INDEX-1) )
!    of the B-splines whose values at X are to be returned.
!    INDEX is used to avoid recalculations when several
!    columns of the triangular array of B-spline values are
!    needed, for example, in BVALUE or in BSPLVD.
!    If INDEX = 1, the calculation starts from scratch and the entire
!    triangular array of B-spline values of orders
!    1, 2, ...,JHIGH is generated order by order, that is,
!    column by column.
!    If INDEX = 2, only the B-spline values of order J+1, J+2, ..., JOUT
!    are generated, the assumption being that BIATX, J,
!    DELTAL, DELTAR are, on entry, as they were on exit
!    at the previous call.  In particular, if JHIGH = 0,
!    then JOUT = J+1, that is, just the next column of B-spline
!    values is generated.
!    Warning: the restriction  JOUT <= JMAX (= 20) is
!    imposed arbitrarily by the dimension statement for DELTAL
!    and DELTAR, but is nowhere checked for.
!
!    Input, real ( kind = 8 ) X, the point at which the B-splines
!    are to be evaluated.
!
!    Input, integer LEFT, an integer chosen so that
!    T(LEFT) <= X <= T(LEFT+1).
!
!    Output, real ( kind = 8 ) BIATX(JOUT), with BIATX(I) containing the
!    value at X of the polynomial of order JOUT which agrees
!    with the B-spline B(LEFT-JOUT+I,JOUT,T) on the interval
!    (T(LEFT),T(LEFT+1)).
!

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
!!  sll_int32 :: icell
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
!!      icell = i+es%num_cells1*(j-1)
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
 
