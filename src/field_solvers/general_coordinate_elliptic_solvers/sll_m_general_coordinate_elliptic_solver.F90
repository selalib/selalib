!> This module is an implementation of a general elliptic solver using 
!> Finite Element Methods with B-Splines used as the basis shape functions.
!>
!> Just a reminder about bsplines
!>
!>  \f$ m+1 \f$ knots \f$ t_i \f$  in \f$ [0,1] \f$ with
!> \f[ 0 \le t_0 \le t_1 \le \ldots \le t_m \le 1 \f]
!> the n degree spline curve is
!> \f[
!> \mathbf{S}(t) = \sum_{i = 0}^{m-n-1} b_{i,n} (t) . \mathbf{P}_{i} \,,\, t \in [0, 1],
!> \f]
!> 
!> where \f$ P_i \f$ is a polynomial function \f$ (m-n) \f$ points.
!> 
!> \f$ m-n \f$ B-splines n degree functions are defined recursively :
!> 
!> \f[
!> b_{j, 0}(t) := \left\{
!> \begin{array}{ll}
!> \mathrm{if} \quad t_j \leqslant t < t_{j+1} & 1 \\\
!> \mathrm{else} \quad  & 0
!> \end{array} \right.
!> \f]
!>
!> \f[
!> b_{j,n}(t) := \frac{t-t_j}{t_{j+n}-t_j} b_{j,n-1}(t)+\frac{t_{j+n+1}-t}{t_{j+n+1}-t_{j+1}}b_{j+1,n-1}(t).
!> \f]
!>
!> @ingroup general_coordinate_elliptic_solvers
!> @brief Elliptic solver on 2d curvilinear mesh
!> @details This solver works with analytical 
!> and discrete coordinate transformations.
module sll_m_general_coordinate_elliptic_solver
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_errors.h"
#include "sll_memory.h"
#include "sll_working_precision.h"

  use iso_fortran_env, only: &
    output_unit

  use sll_m_arbitrary_degree_spline_interpolator_2d, only: &
    sll_t_arbitrary_degree_spline_interpolator_2d

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_dirichlet, &
    sll_p_periodic

  use sll_m_cubic_splines, only: &
    sll_s_compute_cubic_spline_2d, &
    sll_s_get_coeff_cubic_spline_2d, &
    sll_f_new_cubic_spline_2d, &
    sll_t_cubic_spline_2d

  use sll_m_deboor_splines_1d, only: &
    sll_s_bsplvd, &
    sll_t_deboor_type, &
    sll_s_interv

  use sll_m_gauss_legendre_integration, only: &
    sll_f_gauss_legendre_points_and_weights

  use sll_m_gauss_lobatto_integration, only: &
    sll_f_gauss_lobatto_points_and_weights

  use sll_m_interpolators_2d_base, only: &
    sll_c_interpolator_2d

  use sll_m_scalar_field_2d, only: &
    sll_t_scalar_field_2d_analytic, &
    sll_t_scalar_field_2d_discrete

  use sll_m_scalar_field_2d_base, only: &
    sll_c_scalar_field_2d_base

  use sll_m_sparse_matrix, only:          &
    sll_s_csr_add_one_constraint,         &
    sll_f_new_csr_matrix,                 &
    sll_f_new_csr_matrix_with_constraint, &
    sll_s_add_to_csr_matrix,              &
    sll_t_csr_matrix,                     &
    sll_s_free_csr_matrix,                &
    sll_s_factorize_csr_matrix,           &
    sll_s_mult_csr_matrix_vector,         &
    sll_s_solve_csr_matrix,               &
    sll_p_umfpack,                        &
    sll_p_pastix

#ifdef _OPENMP
use omp_lib
#endif
  implicit none

  public :: &
    sll_p_es_gauss_legendre, &
    sll_p_es_gauss_lobatto, &
    sll_s_factorize_mat_es, &
    sll_s_factorize_mat_es_prototype, &
    sll_t_general_coordinate_elliptic_solver, &
    sll_s_initialize_general_elliptic_solver_prototype, &
    sll_f_new_general_elliptic_solver, &
    sll_f_new_general_elliptic_solver_prototype, &
    sll_s_set_rho_coefficients_coordinates_elliptic_eq_prototype, &
    sll_o_create, &
    sll_o_delete, &
    sll_o_solve, &
    sll_s_solve_general_coordinates_elliptic_eq_prototype

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!> @brief
!> General coordinate elliptic solver derived type
!> @details
!> The indexing of the
!> splines in array global_indices depends on the boundary conditions.
!> local_indices includes the changes resulting from the boundary conditions.
!> local_to_lobal_indices(i,j) = global_indices(local_indices(i,j))
type :: sll_t_general_coordinate_elliptic_solver

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
  sll_real64, dimension(:),   pointer :: knots1
  sll_real64, dimension(:),   pointer :: knots2
  sll_real64, dimension(:),   pointer :: knots1_rho
  sll_real64, dimension(:),   pointer :: knots2_rho
  sll_real64, dimension(:,:), pointer :: gauss_pts1
  sll_real64, dimension(:,:), pointer :: gauss_pts2
  sll_int32  :: bc1_min
  sll_int32  :: bc1_max
  sll_int32  :: bc2_min
  sll_int32  :: bc2_max
  sll_int32  :: spline_degree1
  sll_int32  :: spline_degree2
  sll_real64 :: epsi
  sll_real64 :: intjac

  sll_int32, dimension(:,:), pointer :: local_to_global_indices
  sll_int32, dimension(:,:), pointer :: local_to_global_indices_source
  sll_int32, dimension(:,:), pointer :: local_to_global_indices_source_bis

  !!! contains the values of all splines in all gauss points
  sll_real64, dimension(:,:,:,:), pointer :: v_splines1
  sll_real64, dimension(:,:,:,:), pointer :: v_splines2
  sll_real64, dimension(:,:,:), pointer :: v_rho_splines1
  sll_real64, dimension(:,:,:), pointer :: v_rho_splines2
  !sll_real64, dimension(:,:,:,:), pointer :: v_splines1_check
  !sll_real64, dimension(:,:,:,:), pointer :: v_splines2_check

  type(sll_t_csr_matrix),         pointer :: csr_mat
  type(sll_t_csr_matrix),         pointer :: csr_mat_with_constraint
  type(sll_t_csr_matrix),         pointer :: csr_mat_source
  sll_real64, dimension(:),       pointer :: rho_vec
  sll_real64, dimension(:),       pointer :: tmp_rho_vec
  sll_real64, dimension(:),       pointer :: phi_vec
  sll_real64, dimension(:),       pointer :: masse
  sll_real64, dimension(:),       pointer :: stiff
  sll_real64, dimension(:),       pointer :: rho_coeff_1d
  logical                                 :: perper
  type(sll_t_deboor_type)                       :: db
  logical                                 :: precompute_rhs
  type(sll_t_cubic_spline_2d),      pointer :: cubic_spline 
  logical                                 :: use_cubic_splines
  sll_int32                               :: rho_degree1
  sll_int32                               :: rho_degree2
  logical                                 :: with_constraint
  logical                                 :: zero_mean

end type sll_t_general_coordinate_elliptic_solver

!> For the integration mode.  
sll_int32, parameter :: sll_p_es_gauss_legendre = 0
!> For the integration mode.  
sll_int32, parameter :: sll_p_es_gauss_lobatto = 1
  
interface sll_o_delete
  module procedure delete_elliptic
end interface sll_o_delete

interface sll_o_create
  module procedure initialize_general_elliptic_solver
end interface sll_o_create

interface sll_o_solve
  module procedure solve_general_coordinates_elliptic_eq
end interface sll_o_solve

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
!> @param      es the type sll_t_general_coordinate_elliptic_solver
!> @param[in]  spline_degree1 the degre of B-spline in the direction eta1
!> @param[in]  spline_degree2 the degre of B-spline in the direction eta2
!> @param[in]  num_cells1 the number of cells in the direction eta1
!> @param[in]  num_cells2 the number of cells in the direction eta2
!> @param[in]  quadrature_type1 the type of quadrature in the direction eta1
!> @param[in]  quadrature_type2 the type of quadrature in the direction eta2
!> @param[in]  bc1_min the boundary condition at left in the direction eta1
!> @param[in]  bc1_max the boundary condition at right in the direction eta2
!> @param[in]  bc2_min the boundary condition at left in the direction eta2
!> @param[in]  bc2_max the boundary condition at right in the direction eta2
!> @param[in]  eta1_min the minimun in the direction eta1
!> @param[in]  eta1_max the maximun in the direction eta1
!> @param[in]  eta2_min the minimun in the direction eta2
!> @param[in]  eta2_max the maximun in the direction eta2
!> @param[out] the type sll_t_general_coordinate_elliptic_solver
subroutine initialize_general_elliptic_solver( &
       es,                                     &
       spline_degree1,                         &
       spline_degree2,                         &
       num_cells1,                             &
       num_cells2,                             &
       quadrature_type1,                       &
       quadrature_type2,                       &
       bc1_min,                                &
       bc1_max,                                &
       bc2_min,                                &
       bc2_max,                                &
       eta1_min,                               &
       eta1_max,                               &
       eta2_min,                               &
       eta2_max,                                &
       precompute_rhs)
    
type(sll_t_general_coordinate_elliptic_solver), intent(out) :: es

sll_int32,  intent(in) :: spline_degree1
sll_int32,  intent(in) :: spline_degree2
sll_int32,  intent(in) :: num_cells1
sll_int32,  intent(in) :: num_cells2
sll_int32,  intent(in) :: bc1_min
sll_int32,  intent(in) :: bc1_max
sll_int32,  intent(in) :: bc2_min
sll_int32,  intent(in) :: bc2_max
sll_int32,  intent(in) :: quadrature_type1
sll_int32,  intent(in) :: quadrature_type2
sll_real64, intent(in) :: eta1_min
sll_real64, intent(in) :: eta1_max
sll_real64, intent(in) :: eta2_min
sll_real64, intent(in) :: eta2_max
logical,    intent(in), optional :: precompute_rhs

sll_int32 :: knots1_size
sll_int32 :: knots2_size
sll_int32 :: num_splines1
sll_int32 :: num_splines2
sll_int32 :: vec_sz ! for rho_vec and phi_vec allocations
sll_int32 :: ierr
sll_int32 :: solution_size
sll_int32 :: dim1, dim2
sll_int32 :: num_pts_g1, num_pts_g2

sll_real64, allocatable :: work1(:,:)
sll_real64, allocatable :: work2(:,:)
sll_real64, allocatable :: dbs1(:,:)
sll_real64, allocatable :: dbs2(:,:)
sll_real64 :: xg, yg

sll_int32  :: i, j, ii, jj, ispl1, ispl2
sll_real64 :: eta1, eta2, gspl1, gspl2
sll_int32  :: left
sll_int32  :: kk, ll, index_coef1, index_coef2
sll_int32  :: x, y, nbsp, nbsp1
sll_int32  :: index1, index2, index3, index4
sll_int32  :: icell, b, bprime

sll_int32, dimension(:),   allocatable :: tab_index_coeff1
sll_int32, dimension(:),   allocatable :: tab_index_coeff2

sll_int32, dimension(:,:), allocatable :: local_to_global_indices_check
sll_real64, dimension(:,:,:,:), allocatable :: v_splines1_check
sll_real64, dimension(:,:,:,:), allocatable :: v_splines2_check
sll_real64, dimension(:), allocatable :: knots1
sll_real64, dimension(:), allocatable :: knots2
sll_real64 :: err

if(present(precompute_rhs))then
  es%precompute_rhs = precompute_rhs
else
  es%precompute_rhs = .true.
endif

es%total_num_splines_loc = (spline_degree1+1)*(spline_degree2+1)
! The total number of splines in a single direction is given by
! num_cells + spline_degree
num_splines1 = num_cells1 + spline_degree1
num_splines2 = num_cells2 + spline_degree2
  
dim1 = (spline_degree1+1)*(spline_degree2+1)
dim2 = (num_cells1*num_cells2)
SLL_ALLOCATE(es%local_to_global_indices(1:dim1,1:dim2),ierr)
SLL_ALLOCATE(es%local_to_global_indices_source(1:dim1,1:dim2),ierr)
SLL_ALLOCATE(es%local_to_global_indices_source_bis(1:dim1,1:dim2),ierr)

SLL_ALLOCATE(local_to_global_indices_check(1:dim1,1:dim2),ierr)

SLL_ALLOCATE(knots1(num_cells1 + 2*spline_degree1+1),ierr)
SLL_ALLOCATE(knots2(num_cells2 + 2*spline_degree2+1),ierr)

call initconnectivity_new(   &
  num_cells1,                &
  num_cells2,                &
  spline_degree1,            &
  spline_degree2,            &
  bc1_min,                   &
  bc1_max,                   &
  bc2_min,                   &
  bc2_max,                   &
  es%local_to_global_indices )

!old way to compute initconnectivity
!just for check for the moment
!when everything is checked
!this part will be removed (begin)
call initconnectivity(          &
  num_cells1,                   &
  num_cells2,                   &
  spline_degree1,               &
  spline_degree2,               &
  bc1_min,                      &
  bc1_max,                      &
  bc2_min,                      &
  bc2_max,                      &
  local_to_global_indices_check )

ierr = maxval(abs(local_to_global_indices_check-es%local_to_global_indices))
if(ierr/=0)then
  print *,'es%local_to_global_indices=',es%local_to_global_indices
  print *,'local_to_global_indices_check=',local_to_global_indices_check
  SLL_ERROR('initialize_general_elliptic_solver','problem with local_to_global_indices')
endif

!this part will be removed (end)



! This should be changed to verify that the passed BC's are part of the
! recognized list described in sll_m_boundary_condition_descriptors...

es%bc1_min        = bc1_min
es%bc1_max        = bc1_max
es%bc2_min        = bc2_min
es%bc2_max        = bc2_max
es%spline_degree1 = spline_degree1
es%spline_degree2 = spline_degree2
es%num_cells1     = num_cells1
es%num_cells2     = num_cells2
es%delta_eta1     = (eta1_max-eta1_min)/num_cells1
es%delta_eta2     = (eta2_max-eta2_min)/num_cells2
es%eta1_min       = eta1_min
es%eta2_min       = eta2_min

! Allocate and fill the gauss points/weights information.
select case(quadrature_type1)
case (sll_p_es_gauss_legendre)
  SLL_ALLOCATE(es%gauss_pts1(2,spline_degree1+2),ierr)
  es%gauss_pts1 = sll_f_gauss_legendre_points_and_weights(spline_degree1+2)
case (sll_p_es_gauss_lobatto)
  SLL_ALLOCATE(es%gauss_pts1(2,spline_degree1+2),ierr)
  es%gauss_pts1 = sll_f_gauss_lobatto_points_and_weights(spline_degree1+2)
case default
  SLL_ERROR('initialize_general_elliptic_solver&
  &','unknown type of gauss points in the direction 1')
end select
   
select case(quadrature_type2)
case (sll_p_es_gauss_legendre)
  SLL_ALLOCATE(es%gauss_pts2(2,spline_degree2+2),ierr)
  es%gauss_pts2(:,:) = sll_f_gauss_legendre_points_and_weights(spline_degree2+2)
case (sll_p_es_gauss_lobatto)
  SLL_ALLOCATE(es%gauss_pts2(2,spline_degree2+2),ierr)
  es%gauss_pts2(:,:) = sll_f_gauss_lobatto_points_and_weights(spline_degree2+2)
case default
  SLL_ERROR('initialize_general_elliptic_solver&
  &','unknown type of gauss points in the direction 2')
end select

!PN : Gauss points positions and weights are computed in [-1:1] interval
!PN : We need to rescale them.
es%gauss_pts1(1,:) = 0.5_f64*es%delta_eta1*(es%gauss_pts1(1,:)+1.0_f64)
es%gauss_pts1(2,:) = 0.5_f64*es%delta_eta1*es%gauss_pts1(2,:)
es%gauss_pts2(1,:) = 0.5_f64*es%delta_eta2*(es%gauss_pts2(1,:)+1.0_f64)
es%gauss_pts2(2,:) = 0.5_f64*es%delta_eta2*es%gauss_pts2(2,:)

!print *,'#begin'
!print *,'#quadrature',quadrature_type1,quadrature_type2
!flush( output_unit )
!print *,es%gauss_pts1(1,:)
!print *,'#end'
!flush( output_unit )
!print *,es%delta_eta1
!print *,es%gauss_pts1(1,1:spline_degree1+2)
!print *,es%gauss_pts1(1,1:spline_degree1+2)/es%delta_eta1
!flush( output_unit )

es%perper  = .false. 

if( bc1_min == sll_p_periodic .and. bc1_max == sll_p_periodic .and. &
    bc2_min == sll_p_periodic .and. bc2_max == sll_p_periodic ) then

  es%total_num_splines1 = num_cells1 
  es%total_num_splines2 = num_cells2

  knots1_size = 2*spline_degree1+2
  knots2_size = 2*spline_degree2+2
  vec_sz      = num_cells1*num_cells2
  es%perper   = .true. 

else if( bc1_min == sll_p_periodic  .and. bc1_max == sll_p_periodic .and.&
         bc2_min == sll_p_dirichlet .and. bc2_max == sll_p_dirichlet ) then

  es%total_num_splines1 = num_cells1 
  es%total_num_splines2 = num_cells2 + spline_degree2 - 2

  knots1_size = 2*spline_degree1+2
  knots2_size = 2*spline_degree2+num_cells2+1

  vec_sz = num_cells1*(num_cells2+spline_degree2)

else if( bc1_min == sll_p_dirichlet .and. bc1_max == sll_p_dirichlet .and.&
         bc2_min == sll_p_periodic  .and. bc2_max == sll_p_periodic ) then

  es%total_num_splines1 = num_cells1 + spline_degree1 - 2
  es%total_num_splines2 = num_cells2 
  knots1_size = 2*spline_degree1+num_cells1+1
  knots2_size = 2*spline_degree2+2
  vec_sz      = (num_cells1 + spline_degree1)*num_cells2

else if( bc1_min == sll_p_dirichlet .and. bc1_max == sll_p_dirichlet .and.&
         bc2_min == sll_p_dirichlet .and. bc2_max == sll_p_dirichlet ) then

  es%total_num_splines1 = num_cells1 + spline_degree1 - 2
  es%total_num_splines2 = num_cells2 + spline_degree2 - 2
  knots1_size = 2*spline_degree1 + num_cells1+1
  knots2_size = 2*spline_degree2 + num_cells2+1
  vec_sz      = (num_cells1 + spline_degree1)*&
                (num_cells2 + spline_degree2)

end if

!print*,'#now'
!print *,es%delta_eta1
!print *,es%gauss_pts1(1,1:spline_degree1+2)
!print *,es%gauss_pts1(1,1:spline_degree1+2)/es%delta_eta1
!flush( output_unit )

SLL_ALLOCATE(es%knots1(knots1_size),ierr)
SLL_ALLOCATE(es%knots2(knots2_size),ierr)
SLL_ALLOCATE(es%knots1_rho(num_cells1 + spline_degree1 + 2),ierr)
SLL_ALLOCATE(es%knots2_rho(num_cells2 + spline_degree2 + 2),ierr)
SLL_ALLOCATE(es%rho_vec(vec_sz),ierr)
SLL_ALLOCATE(es%masse(vec_sz),ierr)
SLL_ALLOCATE(es%stiff(vec_sz),ierr)

!AB : We must add plus 1 for the dimension of the 
!AB : solution in the case periodic periodic to 
!AB : include the periodicity in the last point.  

solution_size = es%total_num_splines1*es%total_num_splines2

if(es%perper) then
  SLL_ALLOCATE(es%tmp_rho_vec(solution_size+1),ierr)
  SLL_ALLOCATE(es%phi_vec(solution_size+1),ierr)
else
  SLL_ALLOCATE(es%tmp_rho_vec(solution_size),ierr)
  SLL_ALLOCATE(es%phi_vec(solution_size),ierr) 
endif
es%rho_vec = 0.0_f64
es%masse   = 0.0_f64
es%stiff   = 0.0_f64

call initialize_knots( spline_degree1, &
&                      num_cells1,     &
&                      eta1_min,       &
&                      eta1_max,       &
&                      bc1_min,        &
&                      bc1_max,        &
&                      es%knots1 )

call initialize_knots( spline_degree2, &
&                      num_cells2,     &
&                      eta2_min,       &
&                      eta2_max,       &
&                      bc2_min,        &
&                      bc2_max,        &
&                      es%knots2 )

es%csr_mat => sll_f_new_csr_matrix( solution_size,        &
&                             solution_size,              &
&                             num_cells1*num_cells2,      &
&                             es%local_to_global_indices, &
&                             es%total_num_splines_loc,   &
&                             es%local_to_global_indices, &
&                             es%total_num_splines_loc,   &
&                             sll_p_umfpack)

es%knots1_rho(1:spline_degree1+1) = eta1_min
es%knots1_rho(num_cells1+2:num_cells1+1+spline_degree1+1) = eta1_max
 
if ( mod(spline_degree1 +1,2) == 0 ) then
  do i = spline_degree1 +1 + 1, num_cells1 + 1
    es%knots1_rho(i) = eta1_min + (i-(spline_degree1+1)/2-1)*es%delta_eta1
  end do
else
  do i = spline_degree1 +1 + 1, num_cells1 + 1
    es%knots1_rho(i) = &
      0.5*( eta1_min + (i - (spline_degree1)/2-1)*es%delta_eta1 + &
      eta1_min + (i-1 - (spline_degree1)/2 -1)*es%delta_eta1)
  end do
end if

es%knots2_rho(1:spline_degree2+1) = eta2_min
es%knots2_rho(num_cells2+2:num_cells2+1+spline_degree2+1)=eta2_max
    
if (mod(spline_degree2+1,2) == 0 ) then
  do i = spline_degree2 +1 + 1, num_cells2+1
    es%knots2_rho(i) = eta2_min + (i-(spline_degree2+1)/2-1)*es%delta_eta2 
  end do
else
  do i = spline_degree2+1+1, num_cells2+1
    es%knots2_rho(i) = &
      0.5*( eta2_min+(i  -(spline_degree2)/2-1)*es%delta_eta2 + &
            eta2_min+(i-1-(spline_degree2)/2-1)*es%delta_eta2 )
  end do
end if

! allocation of the table containing 
! all values of splines and its
! derivatives in each gauss points
SLL_ALLOCATE(es%v_splines1(3,spline_degree1+1,spline_degree1+2,num_cells1),ierr)
SLL_ALLOCATE(es%v_splines2(3,spline_degree2+1,spline_degree2+2,num_cells2),ierr)

SLL_ALLOCATE(v_splines1_check(2,spline_degree1+1,spline_degree1+2,num_cells1),ierr)
SLL_ALLOCATE(v_splines2_check(2,spline_degree2+1,spline_degree2+2,num_cells2),ierr)


es%v_splines1 = 0.0_f64
es%v_splines2 = 0.0_f64

SLL_ALLOCATE(tab_index_coeff1(num_cells1),ierr)
SLL_ALLOCATE(tab_index_coeff2(num_cells2),ierr)

num_pts_g1 = size(es%gauss_pts1,2)
num_pts_g2 = size(es%gauss_pts2,2)
SLL_ALLOCATE(es%rho_coeff_1d((num_cells1+1)*(num_cells2+1)),ierr)

allocate(work1(spline_degree1+1,spline_degree1+1))
allocate(work2(spline_degree2+1,spline_degree2+1))
allocate(dbs1(spline_degree1+1,2))
allocate(dbs2(spline_degree2+1,2))

do i = 1, es%num_cells1
  eta1  = eta1_min + (i-1)*es%delta_eta1
  do ii=1,num_pts_g1
    xg  = eta1  + es%gauss_pts1(1,ii)
    if(bc1_min==sll_p_periodic .and. bc1_max==sll_p_periodic) then 
      gspl1 = es%gauss_pts1(1,ii)
      ispl1 = spline_degree1+1
    else 
      gspl1 = xg
      ispl1 = spline_degree1+i
    end if
    call sll_s_bsplvd(es%db,es%knots1,spline_degree1+1,gspl1,ispl1,work1,dbs1,2)
    es%v_splines1(1,:,ii,i) = dbs1(:,1)
    es%v_splines1(2,:,ii,i) = dbs1(:,2)
    call sll_s_interv(es%db,es%knots1_rho,es%num_cells1+spline_degree1+2,xg,left,ierr)
    call sll_s_bsplvd(es%db,es%knots1_rho,spline_degree1+1,xg,left,work1,dbs1,1)
    es%v_splines1(3,:,ii,i) = dbs1(:,1)
  end do
  tab_index_coeff1(i) = left
end do

do j = 1, es%num_cells2
  eta2  = eta2_min + (j-1)*es%delta_eta2
  do jj=1,num_pts_g2
    yg  = eta2+es%gauss_pts2(1,jj)
    if (bc2_min==sll_p_periodic .and. bc2_max==sll_p_periodic) then
      gspl2 = es%gauss_pts2(1,jj)
      ispl2 = spline_degree2+1
    else
      gspl2 = yg
      ispl2 = spline_degree2+j
    end if
    call sll_s_bsplvd(es%db,es%knots2,spline_degree2+1,gspl2,ispl2,work2,dbs2,2)
    es%v_splines2(1,:,jj,j) = dbs2(:,1)
    es%v_splines2(2,:,jj,j) = dbs2(:,2)
    call sll_s_interv(es%db,es%knots2_rho,es%num_cells2+spline_degree2+2,yg,left,ierr)
    call sll_s_bsplvd(es%db,es%knots2_rho,spline_degree2+1,yg,left,work2,dbs2,1)
    es%v_splines2(3,:,jj,j) = dbs2(:,1)
  end do
  tab_index_coeff2(j) = left
end do

!deallocate(work1)
!deallocate(work2)
!deallocate(dbs1)
!deallocate(dbs2)

!compute v_splines in another fashion


call compute_global_knots( &
  bc1_min, &
  bc1_max, &
  num_cells1, &
  eta1_min, &
  eta1_max, &
  spline_degree1, &
  knots1)

call compute_global_knots( &
  bc2_min, &
  bc2_max, &
  num_cells2, &
  eta2_min, &
  eta2_max, &
  spline_degree2, &
  knots2)

!print *,es%delta_eta1
!print *,es%gauss_pts1(1,1:spline_degree1+2)
!print *,es%gauss_pts1(1,1:spline_degree1+2)/es%delta_eta1
!stop

call compute_non_zero_splines_and_deriv_at_cell_points( &
  num_cells1, &
  spline_degree1, &
  knots1, &
  es%gauss_pts1(1,:)/es%delta_eta1, & 
  num_pts_g1, &
  v_splines1_check )

call compute_non_zero_splines_and_deriv_at_cell_points( &
  num_cells2, &
  spline_degree2, &
  knots2, &
  es%gauss_pts2(1,:)/es%delta_eta2, & 
  num_pts_g2, &
  v_splines2_check )

  
!check that v_splines is ok
err = maxval(abs(v_splines1_check(1:2,:,:,:)-es%v_splines1(1:2,:,:,:)))  
err = max(err,maxval(abs(v_splines2_check(1:2,:,:,:)-es%v_splines2(1:2,:,:,:))))  

if(err>1.e-8)then
  !print *,'#v_splines1_check=',v_splines1_check
  !print *,'#es%v_splines1=',es%v_splines1
  do i=1,num_cells1
  do j=1,num_pts_g1
  do ii=1,spline_degree1+1    
    print *,1,i,j,ii,es%v_splines1(1,ii,j,i), &
      v_splines1_check(1,ii,j,i),&
      v_splines1_check(1,ii,j,i)-es%v_splines1(1,ii,j,i)
    print *,2,i,j,ii,es%v_splines1(2,ii,j,i), &
      v_splines1_check(2,ii,j,i),&
      v_splines1_check(2,ii,j,i)-es%v_splines1(2,ii,j,i)
  enddo
  enddo
  enddo
  do i=1,num_cells2
  do j=1,num_pts_g2
  do ii=1,spline_degree2+1
    print *,1,i,j,ii,es%v_splines2(1,ii,j,i), &
    v_splines2_check(1,ii,j,i), &
    v_splines2_check(1,ii,j,i)-es%v_splines2(1,ii,j,i)
    print *,2,i,j,ii,es%v_splines2(1,ii,j,i), &
    v_splines2_check(2,ii,j,i), &
    v_splines2_check(2,ii,j,i)-es%v_splines2(2,ii,j,i)
  enddo
  enddo
  enddo
  flush( output_unit )
  gspl1 = es%gauss_pts1(1,1)
  ispl1 = spline_degree1+1
  call sll_s_bsplvd( &
    es%db, &
    es%knots1, &
    spline_degree1+1, &
    gspl1, &
    ispl1, &
    work1, &
    dbs1,&
    2)
  print *,dbs1
  print *,'#err=',err
  flush( output_unit )
  SLL_ERROR('initialize_general_elliptic_solver&
  &','v_splines1_check and es%v_splines1 differ')
endif



do j = 1, es%num_cells2
do i = 1, es%num_cells1
  icell = i + (j-1)*es%num_cells1
  b = 0
  do jj = 0, es%spline_degree2
    index3 = j + jj
    if (bc2_min==sll_p_periodic .and. bc2_max==sll_p_periodic) then 
      if ( index3 > es%total_num_splines2) then
        index3 = index3 - es%total_num_splines2
      end if
    end if
    do ii = 0,es%spline_degree1
      index1 = i + ii
      if (bc1_min==sll_p_periodic .and. bc1_max==sll_p_periodic) then 
        if ( index1 > es%total_num_splines1) then
          index1 = index1 - es%total_num_splines1
        end if
        nbsp = es%total_num_splines1
      else !if (bc1_min==sll_p_dirichlet .and. bc1_max==sll_p_dirichlet) then
        nbsp = es%num_cells1 + es%spline_degree1
      end if
      x = index1 + (index3-1)*nbsp
      b = b+1
      index_coef1 = tab_index_coeff1(i) - es%spline_degree1 + ii
      index_coef2 = tab_index_coeff2(j) - es%spline_degree2 + jj
      es%local_to_global_indices_source(b,icell)= index_coef1 + (index_coef2-1)*(es%num_cells1+1)
      bprime = 0
      do ll = 0,es%spline_degree2
        index4 = j + ll
        if ( (bc2_min==sll_p_periodic).and.(bc2_max== sll_p_periodic))then
          if ( index4 > es%total_num_splines2) then
            index4 = index4 - es%total_num_splines2
          end if
        end if
        do kk = 0,es%spline_degree1
          index2 = i + kk
          if(bc1_min==sll_p_periodic .and. bc1_max==sll_p_periodic)then
            if ( index2 > es%total_num_splines1) then
              index2 = index2 - es%total_num_splines1
            end if
            nbsp1 = es%total_num_splines1
          else !if (bc1_min==sll_p_dirichlet .and. bc1_max==sll_p_dirichlet) then
            nbsp1 = es%num_cells1 + es%spline_degree1
          end if
          y      = index2 + (index4-1)*nbsp1
          bprime = bprime+1 
          es%local_to_global_indices_source_bis(bprime,icell)= y
        end do
      end do
    end do
  end do
end do
end do
DEALLOCATE(tab_index_coeff1)
DEALLOCATE(tab_index_coeff2)

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
!> @param[in] bc1_min the boundary condition at left in the direction eta1
!> @param[in] bc1_max the boundary condition at right in the direction eta2
!> @param[in] bc2_min the boundary condition at left in the direction eta2
!> @param[in] bc2_max the boundary condition at right in the direction eta2
!> @param[in] eta1_min the minimun in the direction eta1
!> @param[in] eta1_max the maximun in the direction eta1
!> @param[in] eta2_min the minimun in the direction eta2
!> @param[in] eta2_max the maximun in the direction eta2
!> @return the type sll_t_general_coordinate_elliptic_solver such that a pointer

function sll_f_new_general_elliptic_solver( spline_degree1,   &
                                      spline_degree2,   &
                                      num_cells1,       &
                                      num_cells2,       &
                                      quadrature_type1, &
                                      quadrature_type2, &
                                      bc1_min,          &
                                      bc1_max,          &
                                      bc2_min,          &
                                      bc2_max,          &
                                      eta1_min,         &
                                      eta1_max,         &
                                      eta2_min,         &
                                      eta2_max,         &
                                      precompute_rhs ) result(es)

type(sll_t_general_coordinate_elliptic_solver), pointer :: es

sll_int32,  intent(in) :: spline_degree1
sll_int32,  intent(in) :: spline_degree2
sll_int32,  intent(in) :: num_cells1
sll_int32,  intent(in) :: num_cells2
sll_int32,  intent(in) :: quadrature_type1
sll_int32,  intent(in) :: quadrature_type2
sll_int32,  intent(in) :: bc1_min
sll_int32,  intent(in) :: bc1_max
sll_int32,  intent(in) :: bc2_min
sll_int32,  intent(in) :: bc2_max
sll_real64, intent(in) :: eta1_min
sll_real64, intent(in) :: eta1_max
sll_real64, intent(in) :: eta2_min
sll_real64, intent(in) :: eta2_max
logical,    intent(in), optional :: precompute_rhs

sll_int32 :: ierr

SLL_ALLOCATE(es,ierr)

call sll_o_create( es,               &
&                spline_degree1,   &
&                spline_degree2,   &
&                num_cells1,       &
&                num_cells2,       &
&                quadrature_type1, &
&                quadrature_type2, &
&                bc1_min,          &
&                bc1_max,          &
&                bc2_min,          &
&                bc2_max,          &
&                eta1_min,         &
&                eta1_max,         &
&                eta2_min,         &
&                eta2_max,         &
&                precompute_rhs )
   
end function sll_f_new_general_elliptic_solver


function sll_f_new_general_elliptic_solver_prototype( spline_degree1,   &
                                      spline_degree2,   &
                                      num_cells1,       &
                                      num_cells2,       &
                                      quadrature_type1, &
                                      quadrature_type2, &
                                      bc1_min,          &
                                      bc1_max,          &
                                      bc2_min,          &
                                      bc2_max,          &
                                      eta1_min,         &
                                      eta1_max,         &
                                      eta2_min,         &
                                      eta2_max,         &
                                      precompute_rhs,   &
                                      rhs_bc1,          &
                                      rhs_bc2,          &
                                      use_cubic_splines, &
                                      rho_degree1, &
                                      rho_degree2, &
                                      with_constraint, &
                                      zero_mean ) result(es)

type(sll_t_general_coordinate_elliptic_solver), pointer :: es

sll_int32,  intent(in) :: spline_degree1
sll_int32,  intent(in) :: spline_degree2
sll_int32,  intent(in) :: num_cells1
sll_int32,  intent(in) :: num_cells2
sll_int32,  intent(in) :: quadrature_type1
sll_int32,  intent(in) :: quadrature_type2
sll_int32,  intent(in) :: bc1_min
sll_int32,  intent(in) :: bc1_max
sll_int32,  intent(in) :: bc2_min
sll_int32,  intent(in) :: bc2_max
sll_real64, intent(in) :: eta1_min
sll_real64, intent(in) :: eta1_max
sll_real64, intent(in) :: eta2_min
sll_real64, intent(in) :: eta2_max
logical,    intent(in), optional :: precompute_rhs
sll_int32,  intent(in), optional :: rhs_bc1
sll_int32,  intent(in), optional :: rhs_bc2
logical, intent(in), optional :: use_cubic_splines
sll_int32, intent(in), optional :: rho_degree1
sll_int32, intent(in), optional :: rho_degree2
logical, intent(in), optional :: with_constraint
logical, intent(in), optional :: zero_mean

sll_int32 :: ierr

SLL_ALLOCATE(es,ierr)

call sll_s_initialize_general_elliptic_solver_prototype( es, &
&                spline_degree1,                       &
&                spline_degree2,                       &
&                num_cells1,                           &
&                num_cells2,                           &
&                quadrature_type1,                     &
&                quadrature_type2,                     &
&                bc1_min,                              &
&                bc1_max,                              &
&                bc2_min,                              &
&                bc2_max,                              &
&                eta1_min,                             &
&                eta1_max,                             &
&                eta2_min,                             &
&                eta2_max,                             &
&                precompute_rhs,                       &
&                rhs_bc1,                              &
&                rhs_bc2,                              &
&                use_cubic_splines,                    &
&                rho_degree1,                          &
&                rho_degree2,                          &
&                with_constraint,                      &
&                zero_mean )
   
end function sll_f_new_general_elliptic_solver_prototype



!> @brief Deallocate the type sll_t_general_coordinate_elliptic_solver
!> @details
!> The parameters are
!> @param[in] es the type sll_t_general_coordinate_elliptic_solver
  
subroutine delete_elliptic( es )
type(sll_t_general_coordinate_elliptic_solver) :: es
sll_int32 :: ierr

SLL_DEALLOCATE(es%knots1,ierr)
SLL_DEALLOCATE(es%knots2,ierr)
SLL_DEALLOCATE(es%gauss_pts1,ierr)
SLL_DEALLOCATE(es%gauss_pts2,ierr)
SLL_DEALLOCATE(es%local_to_global_indices,ierr)
SLL_DEALLOCATE(es%local_to_global_indices_source,ierr)
SLL_DEALLOCATE(es%local_to_global_indices_source_bis,ierr)
call sll_s_free_csr_matrix(es%csr_mat)
if (associated(es%csr_mat_source))  &
  call sll_s_free_csr_matrix(es%csr_mat_source)
if (associated(es%csr_mat_with_constraint))  &
  call sll_s_free_csr_matrix(es%csr_mat_with_constraint)
SLL_DEALLOCATE(es%rho_vec,ierr)
SLL_DEALLOCATE(es%tmp_rho_vec,ierr)
SLL_DEALLOCATE(es%phi_vec,ierr)
SLL_DEALLOCATE(es%masse,ierr)
SLL_DEALLOCATE(es%stiff,ierr)
SLL_DEALLOCATE(es%knots1_rho,ierr)
SLL_DEALLOCATE(es%knots2_rho,ierr)
SLL_DEALLOCATE(es%v_splines1,ierr)
SLL_DEALLOCATE(es%v_splines2,ierr)

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
!> @param es the type sll_t_general_coordinate_elliptic_solver
!> @param[in] a11_field_mat the field corresponding to the matrix coefficient A(1,1)
!> @param[in] a12_field_mat the field corresponding to the matrix coefficient A(1,2)
!> @param[in] a21_field_mat the field corresponding to the matrix coefficient A(2,1)
!> @param[in] a22_field_mat the field corresponding to the matrix coefficient A(2,2)
!> @param[in] b1_field_vect the field corresponding to the vector coefficient B(1)
!> @param[in] b2_field_vect the field corresponding to the vector coefficient B(2)
!> @param[in] c_field the field corresponding to the coefficient B(1) of the scalar C
!> @return the type sll_t_general_coordinate_elliptic_solver contains the matrix 
!> to solve the equation
  
subroutine sll_s_factorize_mat_es( es,            &
                             a11_field_mat, &
                             a12_field_mat, &
                             a21_field_mat, &
                             a22_field_mat, &
                             b1_field_vect, &
                             b2_field_vect, &
                             c_field)

type(sll_t_general_coordinate_elliptic_solver),intent(inout) :: es

class(sll_c_scalar_field_2d_base), pointer :: a11_field_mat
class(sll_c_scalar_field_2d_base), pointer :: a12_field_mat
class(sll_c_scalar_field_2d_base), pointer :: a21_field_mat
class(sll_c_scalar_field_2d_base), pointer :: a22_field_mat
class(sll_c_scalar_field_2d_base), pointer :: b1_field_vect
class(sll_c_scalar_field_2d_base), pointer :: b2_field_vect
class(sll_c_scalar_field_2d_base), pointer :: c_field

sll_real64, dimension(:,:),   allocatable :: M_c
sll_real64, dimension(:,:),   allocatable :: K_11
sll_real64, dimension(:,:),   allocatable :: K_12
sll_real64, dimension(:,:),   allocatable :: K_21
sll_real64, dimension(:,:),   allocatable :: K_22
sll_real64, dimension(:,:),   allocatable :: M_bv
sll_real64, dimension(:,:),   allocatable :: S_b1
sll_real64, dimension(:,:),   allocatable :: S_b2  
sll_real64, dimension(:),     allocatable :: mass
sll_real64, dimension(:),     allocatable :: stif
sll_real64, dimension(:,:,:), pointer     :: source

sll_int32 :: ierr
sll_int32 :: i
sll_int32 :: j
sll_int32 :: icell
sll_int32 :: bc1_min
sll_int32 :: bc1_max
sll_int32 :: bc2_min
sll_int32 :: bc2_max

sll_real64 :: delta1
sll_real64 :: delta2
sll_real64 :: eta1_min
sll_real64 :: eta2_min
sll_real64 :: eta1
sll_real64 :: eta2
sll_int32  :: num_pts_g1
sll_int32  :: num_pts_g2
sll_int32  :: ii,kk,mm
sll_int32  :: jj,ll,nn
sll_int32  :: ig, jg
sll_real64 :: xg, wxg
sll_real64 :: yg, wyg
sll_int32  :: index1
sll_int32  :: index2
sll_real64 :: val_c
sll_real64 :: val_a11
sll_real64 :: val_a12
sll_real64 :: val_a21
sll_real64 :: val_a22
sll_real64 :: val_b1=0.0_f64
sll_real64 :: val_b1_der1=0.0_f64
sll_real64 :: val_b1_der2=0.0_f64
sll_real64 :: val_b2=0.0_f64
sll_real64 :: val_b2_der1=0.0_f64
sll_real64 :: val_b2_der2=0.0_f64
sll_real64 :: jac_mat(2,2)
sll_real64 :: val_jac
sll_real64 :: B11
sll_real64 :: B12
sll_real64 :: B21
sll_real64 :: B22
sll_real64 :: MC
sll_real64 :: C1
sll_real64 :: C2    
sll_int32  :: index3, index4
!sll_int32  :: index_coef2
sll_int32  :: b, bprime,x
sll_int32  :: a, aprime
sll_real64 :: elt_mat_global
sll_int32  :: nspl, nbsp,nbsp1
sll_int32  :: spl_deg_1, spl_deg_2, nc_1, nc_2
sll_int32  :: ideg2,ideg1
sll_int32  :: jdeg2,jdeg1
sll_real64 :: v1, v2, v3, v4, d1, d2, d3, d4, r1, r2
sll_real64 :: wxy
sll_real64 :: wxy_by_val_jac 
sll_real64 :: wxy_val_jac 
sll_real64 :: r1r2, v1v2, d1v2, v1d2
sll_real64 :: d3v4, v3v4 , v3d4
sll_real64 :: intjac

!sll_real64, allocatable :: dense_matrix(:,:)

!sll_int32 :: tid=0
!sll_int32 :: nthreads=1

bc1_min    = es%bc1_min
bc1_max    = es%bc1_max
bc2_min    = es%bc2_min
bc2_max    = es%bc2_max
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

if( bc1_min == sll_p_periodic .and. bc1_max == sll_p_periodic .and. &
    bc2_min == sll_p_periodic .and. bc2_max == sll_p_periodic ) then
   es%perper = .true.
   !SLL_WARNING("sll_general_coordinate_elliptic_solver","The full periodic version is deprecated")
else
   es%perper = .false.  
end if   

intjac = 0.0_f64

SLL_CLEAR_ALLOCATE(source(1:nspl,1:nspl,1:nc_1*nc_2),ierr)
SLL_CLEAR_ALLOCATE(M_c(1:nspl,1:nspl),ierr)
SLL_CLEAR_ALLOCATE(K_11(1:nspl,1:nspl),ierr)
SLL_CLEAR_ALLOCATE(K_12(1:nspl,1:nspl),ierr)
SLL_CLEAR_ALLOCATE(K_21(1:nspl,1:nspl),ierr)
SLL_CLEAR_ALLOCATE(K_22(1:nspl,1:nspl),ierr)
SLL_CLEAR_ALLOCATE(S_b1(1:nspl,1:nspl),ierr)
SLL_CLEAR_ALLOCATE(S_b2(1:nspl,1:nspl),ierr)
SLL_CLEAR_ALLOCATE(M_bv(1:nspl,1:nspl),ierr)
SLL_CLEAR_ALLOCATE(mass(1:nspl),ierr)
SLL_CLEAR_ALLOCATE(stif(1:nspl),ierr)

!!$OMP PARALLEL DEFAULT(NONE) &
!!$OMP SHARED( es, c_field, &
!!$OMP a11_field_mat, a12_field_mat, a21_field_mat, a22_field_mat, &
!!$OMP b1_field_vect, b2_field_vect, spl_deg_1, spl_deg_2, source, intjac ) &
!!$OMP FIRSTPRIVATE(nc_1, nc_2, delta1, delta2, eta1_min, eta2_min, &
!!$OMP bc1_min, bc1_max, bc2_min, bc2_max, num_pts_g1, num_pts_g2)  &
!!$OMP PRIVATE(nthreads, tid, &
!!$OMP i,j,eta1,eta2,mass,stif, &
!!$OMP yg,wyg,jg,ig,xg,wxg, &
!!$OMP wxy,val_c,val_a11,val_a12,val_a21,val_a22, &
!!$OMP val_b1, val_b1_der1, val_b1_der2, &
!!$OMP val_b2, val_b2_der1, val_b2_der2, &
!!$OMP icell, ii, jj, kk, ll, mm, nn,   &
!!$OMP jac_mat, val_jac, wxy_by_val_jac, wxy_val_jac, &
!!$OMP B11, B12, B21, B22, MC, C1, C2, &
!!$OMP v1, v2, v3, v4, r1, r2, d1, d2, d3, d4, &
!!$OMP v3v4, d3v4, v3d4, &
!!$OMP M_c,K_11,K_12,K_21,K_22,M_bv,S_b1,S_b2, &
!!$OMP index_coef2, index_coef1, &
!!$OMP index1, index2, index3, index4, &
!!$OMP a, b, x, y, aprime, bprime, nbsp, nbsp1, &
!!$OMP r1r2, v1v2, d1v2, v1d2, elt_mat_global )
!
!!$ tid = omp_get_thread_num()
!!$ nthreads = omp_get_num_threads()
!!$ if (tid == 0) print *, 'Number of threads = ', nthreads
!!$OMP DO SCHEDULE(STATIC,nc_2/nthreads) REDUCTION(+:intjac)
do j = 1, nc_2
do i = 1, nc_1
        
  icell = i + (j-1) * nc_1
  eta1  = eta1_min + (i-1)*delta1
  eta2  = eta2_min + (j-1)*delta2
    
  mass  = 0.0_f64
  stif  = 0.0_f64
  M_c   = 0.0_f64
  K_11  = 0.0_f64
  K_12  = 0.0_f64
  K_21  = 0.0_f64
  K_22  = 0.0_f64
  M_bv  = 0.0_f64
  S_b1  = 0.0_f64
  S_b2  = 0.0_f64

  do jg=1,num_pts_g2
  
    yg  = eta2+es%gauss_pts2(1,jg)
    wyg = es%gauss_pts2(2,jg)
  
    do ig=1,num_pts_g1
    
      xg  = eta1+es%gauss_pts1(1,ig)
      wxg = es%gauss_pts1(2,ig)

      wxy = wxg*wyg

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
      val_jac = jac_mat(1,1)*jac_mat(2,2)-jac_mat(1,2)*jac_mat(2,1)

      wxy_by_val_jac = wxy/val_jac
      wxy_val_jac    = wxy*val_jac

      val_c = val_c * wxy_val_jac
        
      intjac = intjac + wxy_val_jac

      ! The B matrix is  by (J^(-1)) A^t (J^(-1))^t 
      B11 = jac_mat(2,2)*jac_mat(2,2)*val_a11 - &
            jac_mat(2,2)*jac_mat(1,2)*(val_a12+val_a21) + &
            jac_mat(1,2)*jac_mat(1,2)*val_a22

      B11 = B11*wxy_by_val_jac
          
      B21 = jac_mat(1,1)*jac_mat(2,2)*val_a12 - &
            jac_mat(1,1)*jac_mat(1,2)*val_a22 - &
            jac_mat(2,1)*jac_mat(2,2)*val_a11 + &
            jac_mat(1,2)*jac_mat(2,1)*val_a21

      B21 = B21*wxy_by_val_jac
        
      B12 = jac_mat(1,1)*jac_mat(2,2)*val_a21 - &
            jac_mat(1,1)*jac_mat(1,2)*val_a22 - &
            jac_mat(2,1)*jac_mat(2,2)*val_a11 + &
            jac_mat(1,2)*jac_mat(2,1)*val_a12

      B12 = B12*wxy_by_val_jac

      B22 = jac_mat(1,1)*jac_mat(1,1)*val_a22 - &
            jac_mat(1,1)*jac_mat(2,1)*(val_a21+val_a12) + &
            jac_mat(2,1)*jac_mat(2,1)*val_a11
          
      B22 = B22*wxy_by_val_jac

      MC =  jac_mat(2,2)*val_b1_der1 &
          - jac_mat(2,1)*val_b1_der2 &
          - jac_mat(1,2)*val_b2_der1 &
          + jac_mat(1,1)*val_b2_der2

      MC = MC*wxy
          
      C1 = jac_mat(2,2)*val_b1-jac_mat(1,2)*val_b2 
      C2 = jac_mat(1,1)*val_b2-jac_mat(2,1)*val_b1

      C1 = C1*wxy
      C2 = C2*wxy
         
      mm = 0
      do jj = 1,spl_deg_2+1

        v2 = es%v_splines2(1,jj,jg,j)
        d2 = es%v_splines2(2,jj,jg,j)
        r2 = es%v_splines2(3,jj,jg,j)

        do ii = 1,spl_deg_1+1
              
          mm = mm+1
                
          v1 = es%v_splines1(1,ii,ig,i)
          d1 = es%v_splines1(2,ii,ig,i)
          r1 = es%v_splines1(3,ii,ig,i)

          r1r2 = r1*r2*wxy_val_jac
          v1v2 = v1*v2
          d1v2 = d1*v2
          v1d2 = v1*d2

          mass(mm) = mass(mm)+v1v2*wxy_val_jac 
          stif(mm) = stif(mm)+wxy_val_jac*(d1v2+v1d2)
               
          nn = 0
          do ll = 1,spl_deg_2+1

            v4 = es%v_splines2(1,ll,jg,j)
            d4 = es%v_splines2(2,ll,jg,j)

            do kk = 1,spl_deg_1+1
                        
              v3 = es%v_splines1(1,kk,ig,i)
              d3 = es%v_splines1(2,kk,ig,i)
              v3v4 = v4*v3
              d3v4 = d3*v4
              v3d4 = v3*d4
              nn = nn+1
             
              source(nn,mm,icell) = source(nn,mm,icell) + r1r2*v3v4
                   
              M_c (nn,mm)=M_c (nn,mm) + val_c* v1v2 * v3v4
              K_11(nn,mm)=K_11(nn,mm) + B11  * d1v2 * d3v4
              K_22(nn,mm)=K_22(nn,mm) + B22  * v1d2 * v3d4
              K_12(nn,mm)=K_12(nn,mm) + B12  * d1v2 * v3d4
              K_21(nn,mm)=K_21(nn,mm) + B21  * v1d2 * d3v4
              M_bv(nn,mm)=M_bv(nn,mm) + MC   * v1v2 * v3v4
              S_b1(nn,mm)=S_b1(nn,mm) + C1   * v1v2 * d3v4
              S_b2(nn,mm)=S_b2(nn,mm) + C2   * v1v2 * v3d4

            end do
          end do
        end do
      end do
    end do
  end do

  do jj = 0, spl_deg_2

    index3 = j + jj
    if (bc2_min==sll_p_periodic .and. bc2_max==sll_p_periodic) then 
      if ( index3 > es%total_num_splines2) then
        index3 = index3 - es%total_num_splines2
      end if
    end if
     
    do ii = 0,spl_deg_1
        
      index1 = i + ii
      if (bc1_min==sll_p_periodic .and. bc1_max==sll_p_periodic) then 
        if ( index1 > es%total_num_splines1) then
          index1 = index1 - es%total_num_splines1
        end if
        nbsp = es%total_num_splines1
      else !if (bc1_min==sll_p_dirichlet .and. bc1_max==sll_p_dirichlet) then
        nbsp = nc_1 + spl_deg_1
      end if

      x = index1 + (index3-1)*nbsp
      b = ii+1+jj*(spl_deg_1+1)
      a = es%local_to_global_indices(b, icell)
         
      es%masse(x) = es%masse(x) + mass(b)
      es%stiff(x) = es%stiff(x) + stif(b)

      do ll = 0,spl_deg_2
        index4 = j + ll
        if ( (bc2_min==sll_p_periodic).and.(bc2_max== sll_p_periodic))then
          if ( index4 > es%total_num_splines2) then
            index4 = index4 - es%total_num_splines2
          end if
        end if
        do kk = 0,spl_deg_1
          index2 = i + kk
          if(bc1_min==sll_p_periodic .and. bc1_max==sll_p_periodic)then
            if ( index2 > es%total_num_splines1) then
              index2 = index2 - es%total_num_splines1
            end if
            nbsp1 = es%total_num_splines1
          else !if (bc1_min==sll_p_dirichlet .and. bc1_max==sll_p_dirichlet) then
            nbsp1 = nc_1 + spl_deg_1
          end if
                
          bprime = kk+1+ll*(spl_deg_1+1)
          aprime = es%local_to_global_indices(bprime,icell)

          if ( a>0 .and. aprime>0 ) then
            elt_mat_global = M_c (bprime,b) - &
                             K_11(bprime,b) - &
                             K_12(bprime,b) - &
                             K_21(bprime,b) - &
                             K_22(bprime,b) - &
                             M_bv(bprime,b) - &
                             S_b1(bprime,b) - &
                             S_b2(bprime,b)

            call sll_s_add_to_csr_matrix(es%csr_mat, elt_mat_global, a, aprime)   
          end if
        end do
      end do
    end do
  end do
end do
end do

!!$OMP END PARALLEL

!#ifdef DEBUG
!allocate(dense_matrix(es%csr_mat%num_rows,es%csr_mat%num_cols))
!
!do i =1, es%csr_mat%num_rows
! do k = es%csr_mat%row_ptr(i), es%csr_mat%row_ptr(i+1)-1
!    j = es%csr_mat%col_ind(k)
!    dense_matrix(i,j) = es%csr_mat%val(k)
! end do
!end do
!
!write(*,"(3x,36i7)") (k, k = 1, size(es%tmp_rho_vec))
!do i = 1, es%csr_mat%num_rows
!  write(*, "(i3,36f7.3)')") i, dense_matrix(i,:) 
!end do
!#endif /* DEBUG */

es%intjac = intjac
print *,'#begin of sll_s_factorize_csr_matrix'

if (es%perper) then

 es%csr_mat_with_constraint => &
   sll_f_new_csr_matrix_with_constraint(es%csr_mat)  

 call sll_s_csr_add_one_constraint( es%csr_mat%row_ptr,           &  
                              es%csr_mat%col_ind,                 &
                              es%csr_mat%val,                     &
                              es%csr_mat%num_rows,                &
                              es%csr_mat%num_nz,                  &
                              es%masse,                           &
                              es%csr_mat_with_constraint%row_ptr, &
                              es%csr_mat_with_constraint%col_ind, &
                              es%csr_mat_with_constraint%val)  

  call sll_s_factorize_csr_matrix(es%csr_mat_with_constraint)

else   

  call sll_s_factorize_csr_matrix(es%csr_mat)

end if 

print *,'#end of sll_s_factorize_csr_matrix'

es%csr_mat_source =>                                     &
  sll_f_new_csr_matrix( size(es%masse,1),                &
                  (nc_1+1)*(nc_2+1),                     &
                  nc_1*nc_2,                             &
                  es%local_to_global_indices_source_bis, &
                  nspl,                                  &
                  es%local_to_global_indices_source,     &
                  nspl )

icell = 0
do j=1,es%num_cells2
do i=1,es%num_cells1
      
  icell = icell+1
  b     = 0
  do ideg2 = 0,es%spline_degree2
  do ideg1 = 0,es%spline_degree1
            
    b = b+1
    a = es%local_to_global_indices_source_bis(b, icell)
        
    bprime = 0
    do jdeg2 = 0,es%spline_degree2
    do jdeg1 = 0,es%spline_degree1
              
      bprime = bprime+1
      aprime = es%local_to_global_indices_source(bprime,icell)
           
      if ( a > 0 .and. aprime > 0) then
        elt_mat_global = source(b,bprime,icell)
        call sll_s_add_to_csr_matrix(es%csr_mat_source,elt_mat_global,a,aprime)
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
SLL_DEALLOCATE_ARRAY(stif,ierr) 
SLL_DEALLOCATE_ARRAY(mass,ierr) 
   
end subroutine sll_s_factorize_mat_es

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
!> @param[in]  es  the type sll_t_general_coordinate_elliptic_solver
!> @param[in]  rho \f$ \rho \f$ the field corresponding to the source term   
!> @param[out] phi \f$ \phi \f$ the field corresponding to the solution of the equation
!> @return     phi the field solution of the equation
  
subroutine solve_general_coordinates_elliptic_eq( es, rho, phi)

class(sll_t_general_coordinate_elliptic_solver), intent(inout)      :: es
class(sll_t_scalar_field_2d_discrete),       intent(inout)      :: phi
class(sll_c_scalar_field_2d_base),           intent(in),target  :: rho

class(sll_c_scalar_field_2d_base), pointer :: base_field_pointer
sll_real64, dimension(:,:),      pointer :: coeff_rho
sll_real64, dimension(:), allocatable :: m_rho_loc

sll_int32  :: i
sll_int32  :: j
sll_int32  :: k
sll_int32  :: num_pts_g1
sll_int32  :: num_pts_g2
sll_int32  :: x, n, b
sll_int32  :: ii, jj, kk, ll, mm, nn
sll_int32  :: bc1_min, bc1_max, bc2_min, bc2_max
sll_int32  :: index1, index3, nbsp

sll_real64 :: wgpt1, wgpt2, gpt1, gpt2, eta1, eta2
sll_real64 :: val_f, val_j, valfj, jac_mat(2,2)
sll_real64 :: spline1, spline2

sll_real64 :: int_rho
sll_real64 :: int_jac
sll_int32  :: ierr
sll_int32  :: nc_1, nc_2

!!$ sll_int32  :: tid = 0
!!$ sll_int32  :: nthreads = 1
  
num_pts_g1 = size(es%gauss_pts1,2)
num_pts_g2 = size(es%gauss_pts2,2)


nc_1            = es%num_cells1
nc_2            = es%num_cells2
es%rho_vec      = 0.0_f64
es%rho_coeff_1d = 0.0_f64
 
bc1_min = es%bc1_min
bc1_max = es%bc1_max
bc2_min = es%bc2_min
bc2_max = es%bc2_max


base_field_pointer => rho




select type( type_field => base_field_pointer)

class is (sll_t_scalar_field_2d_discrete)
  
  
  if(es%precompute_rhs)then
    coeff_rho => type_field%interp_2d%get_coefficients()
            
    do j=1,es%num_cells2+1
      do i=1,es%num_cells1+1
        es%rho_coeff_1d(i+(es%num_cells1+1)*(j-1)) = coeff_rho(i,j)
      end do
    end do

    call sll_s_mult_csr_matrix_vector(es%csr_mat_source,es%rho_coeff_1d,es%rho_vec)

    if(es%perper) then
      es%rho_vec = es%rho_vec - sum(es%rho_vec)/es%intjac*es%masse
    end if
  endif
class is (sll_t_scalar_field_2d_analytic)
   es%precompute_rhs = .false. 
end select
      

if(es%precompute_rhs .eqv. .false.)then  
  int_rho = 0.0_f64
  int_jac = 0.0_f64

  SLL_CLEAR_ALLOCATE(M_rho_loc(1:es%total_num_splines_loc),ierr)

  !!$OMP PARALLEL &
  !!$OMP FIRSTPRIVATE(nc_1, nc_2, num_pts_g1, num_pts_g2, &
  !!$OMP              bc1_min,bc1_max,bc2_min,bc2_max,    &
  !!$OMP              tid, nthreads)                      &
  !!$OMP PRIVATE(i,j,ii,jj,kk,ll,mm,nn,n,m_rho_loc,x,b,   &
  !!$OMP         index1,index3,nbsp,eta1,eta2,gpt1,gpt2,  &
  !!$OMP         wgpt1,wgpt2,spline1,spline2,val_f,val_j, &
  !!$OMP         valfj,jac_mat)
  !!$ tid = omp_get_thread_num()
  !!$ nthreads = omp_get_num_threads()
  !!$OMP MASTER
  !!$ print *, 'Number of threads = ', nthreads
  !!$OMP END MASTER
  !
  !!$OMP DO SCHEDULE(STATIC,nc_2/nthreads) REDUCTION(+:int_rho,int_jac)
  do j=1, nc_2
    do i=1, nc_1
      M_rho_loc = 0.0_f64
      eta1  = es%eta1_min + (i-1)*es%delta_eta1
      eta2  = es%eta2_min + (j-1)*es%delta_eta2
      do jj=1,num_pts_g2
        gpt2  = eta2 + es%gauss_pts2(1,jj)
        wgpt2 = es%gauss_pts2(2,jj)
        do ii=1,num_pts_g1
          gpt1  = eta1 + es%gauss_pts1(1,ii)
          wgpt1 = es%gauss_pts1(2,ii)
      
          val_f   = rho%value_at_point(gpt1,gpt2)
          jac_mat = rho%get_jacobian_matrix(gpt1,gpt2)
          val_j   = jac_mat(1,1)*jac_mat(2,2)-jac_mat(1,2)*jac_mat(2,1)
          val_j   = val_j*wgpt1*wgpt2
          valfj   = val_f*val_j
          int_rho = int_rho + valfj 
          int_jac = int_jac + val_j
      
          do ll = 1,es%spline_degree2+1
            spline2 = es%v_splines2(1,ll,jj,j)*valfj
            do kk = 1,es%spline_degree1+1
              spline1 = es%v_splines1(1,kk,ii,i)
              n = kk + (ll-1)*(es%spline_degree1+1)
              M_rho_loc(n)= M_rho_loc(n) + spline1*spline2
            end do
          end do
        end do
      end do

      do mm = 0,es%spline_degree2
        index3 = j + mm
        if (bc2_min==sll_p_periodic .and. bc2_max==sll_p_periodic) then    
          if ( index3 > es%total_num_splines2) then
            index3 = index3 - es%total_num_splines2
          end if
        end if
      
        do nn = 0,es%spline_degree1
          index1 = i + nn
          if (bc1_min==sll_p_periodic .and. bc1_max==sll_p_periodic) then 
            if (index1 > es%total_num_splines1) then
              index1=index1-es%total_num_splines1
            end if
            nbsp = es%total_num_splines1
          else if (bc1_min==sll_p_dirichlet .and. bc1_max==sll_p_dirichlet) then
            nbsp = es%num_cells1 + es%spline_degree1
          end if
      
          x             =  index1 + (index3-1)*nbsp
          b             =  nn + 1 + mm * (es%spline_degree1+1)
          es%rho_vec(x) =  es%rho_vec(x)  + M_rho_loc(b)
            
        end do
      end do
    end do
  end do

  !!$OMP END PARALLEL

  !PN fix bug found by Michel
  !if (es%perper) es%rho_vec = es%rho_vec - int_rho/int_jac
  if (es%perper)  then
    es%rho_vec = es%rho_vec - sum(es%rho_vec)/es%intjac*es%masse
  end if
     
end if

bc1_min = es%bc1_min
bc1_max = es%bc1_max
bc2_min = es%bc2_min
bc2_max = es%bc2_max
 
es%tmp_rho_vec(:) = 0.0_f64  !PN: Is it useful ?
es%phi_vec(:)     = 0.0_f64  !PN: Is it useful ?
  
if( bc1_min==sll_p_periodic  .and. bc1_max==sll_p_periodic .and. &
    bc2_min==sll_p_dirichlet .and. bc2_max==sll_p_dirichlet ) then
     
  k = 0
  do j = 1, es%total_num_splines2
    do i = 1, es%total_num_splines1
      k = k+1
      es%tmp_rho_vec(k) = es%rho_vec(i+(es%total_num_splines1)*j)
    end do
  end do
     
else if( bc1_min==sll_p_dirichlet .and. bc1_max==sll_p_dirichlet .and.&
         bc2_min==sll_p_dirichlet .and. bc2_max==sll_p_dirichlet) then 
     
  k = 0
  do j = 1, es%total_num_splines2
    do i = 1, es%total_num_splines1
      k = k+1
      es%tmp_rho_vec(k) = es%rho_vec(i+1+(es%total_num_splines1+2)*j)
    end do
  end do
     
else if(bc1_min==sll_p_periodic .and. bc1_max==sll_p_periodic .and.&
        bc2_min==sll_p_periodic .and. bc2_max==sll_p_periodic) then
     
  es%tmp_rho_vec(1:es%total_num_splines1*es%total_num_splines2)=&
    es%rho_vec(1:es%total_num_splines1*es%total_num_splines2) 
     
else if(bc1_min==sll_p_dirichlet .and. bc1_max==sll_p_dirichlet .and.&
        bc2_min==sll_p_periodic  .and. bc2_max==sll_p_periodic ) then
   
  k = 0
  do j = 1, es%total_num_splines2
    do i = 1, es%total_num_splines1
      k = k+1
      es%tmp_rho_vec(k) = es%rho_vec(i+1+(es%total_num_splines1+2)*(j-1))
    end do
  end do

end if

if(bc1_min==sll_p_periodic .and. bc1_max==sll_p_periodic .and.&
   bc2_min==sll_p_periodic .and. bc2_max==sll_p_periodic) then

  call sll_s_solve_csr_matrix(es%csr_mat_with_constraint, &
                              es%tmp_rho_vec,             &
                              es%phi_vec)

  !SLL_WARNING("sll_general_coordinate_elliptic_solver","Use sll_gces_full_periodic")
else

  call sll_s_solve_csr_matrix(es%csr_mat, es%tmp_rho_vec, es%phi_vec)

endif

!print *,'#solve_linear_system done'
  
call phi%interp_2d%set_coefficients(es%phi_vec(1:es%total_num_splines1*es%total_num_splines2))

end subroutine solve_general_coordinates_elliptic_eq
  
!PN In this subroutine we compute the matrix row corresponding to the degree of
!PN freedom in the mesh
!PN local_indices and global_indices should be allocated and initializad
!PN in this subroutine

subroutine initconnectivity( num_cells1,             &
&                            num_cells2,             &
&                            spline_degree1,         &
&                            spline_degree2,         &
&                            bc1_min,                &
&                            bc1_max,                &
&                            bc2_min,                &
&                            bc2_max,                &
&                            local_to_global_indices )
  
sll_int32, intent(in) :: num_cells1
sll_int32, intent(in) :: num_cells2
sll_int32, intent(in) :: spline_degree1
sll_int32, intent(in) :: spline_degree2
sll_int32, intent(in) :: bc1_min
sll_int32, intent(in) :: bc1_max
sll_int32, intent(in) :: bc2_min
sll_int32, intent(in) :: bc2_max

sll_int32, dimension(:,:), allocatable :: local_indices
sll_int32, dimension(:)  , allocatable :: global_indices
sll_int32, dimension(:,:), intent(out) :: local_to_global_indices

sll_int32 :: nb_spl_x
sll_int32 :: nb_spl_y
sll_int32 :: ii, jj
sll_int32 :: bloc
sll_int32 :: cell
sll_int32 :: e
sll_int32 :: b
sll_int32 :: d
sll_int32 :: i
sll_int32 :: j
sll_int32 :: a
sll_int32 :: l
sll_int32 :: ierr

SLL_ALLOCATE(local_indices( 1:(spline_degree1+1)*(spline_degree2+1),1:(num_cells1*num_cells2)),ierr)
SLL_ALLOCATE(global_indices((num_cells1+spline_degree1)*(num_cells2+spline_degree2)),ierr)

global_indices          = 0
local_indices           = 0
local_to_global_indices = 0
  
do j=1, num_cells2    
  do i=1, num_cells1  
    cell = (j-1)*num_cells1 + i
    do jj = 0 , spline_degree2
      do ii = 0, spline_degree1

        bloc = jj * (spline_degree1 + 1) + ii + 1

        b = cell + (j-1)*spline_degree1 + jj*(num_cells1+spline_degree1) + ii 

        local_indices(bloc, cell) = b

      end do
    end do
  end do
end do

nb_spl_x = num_cells1 + spline_degree1
nb_spl_y = num_cells2 + spline_degree2 
  
if( bc1_min == sll_p_periodic .and. bc1_max == sll_p_periodic .and.&
    bc2_min == sll_p_periodic .and. bc2_max == sll_p_periodic ) then

  d = 0
  do j = 1, nb_spl_y
    do i = 1, nb_spl_x
    
      a = i + nb_spl_x*(j-1)
      if (i /= nb_spl_x .and. j/= nb_spl_y) then
        if (global_indices(a) == 0) then
          d = d+1
          global_indices(a) = d
        end if
        if ( 1 <= i .and. i <= spline_degree1) then
          global_indices(a+nb_spl_x-spline_degree1) = global_indices(a)
        end if
        if ( 1 <= j .and. j <= spline_degree2) then
          l = (nb_spl_y - spline_degree2) * nb_spl_x
          global_indices(a+l) = global_indices(a)
        end if
      end if
    
    end do
  end do
  
  global_indices(nb_spl_y*nb_spl_x) = global_indices(nb_spl_x*(nb_spl_y-1) + spline_degree2)
      
else if (bc1_min == sll_p_periodic  .and. bc1_max == sll_p_periodic .and. &
         bc2_min == sll_p_dirichlet .and. bc2_max == sll_p_dirichlet) then
     
  d = 0
  do j = 1, nb_spl_y
    do i = 1, nb_spl_x
      a = i + nb_spl_x*(j-1)
      if ( j == 1 .OR. j == nb_spl_y) then
        global_indices(a) = 0
      else
        if (i /= nb_spl_x) then
          if (global_indices(a) == 0) then
            d = d + 1
            global_indices(a) = d
          end if
          if ( 1 <= i .and. i <= spline_degree1 ) then
            global_indices(a+nb_spl_x-spline_degree1) = global_indices(a)
          end if
        end if
      end if
    end do
  end do

else if(bc1_min == sll_p_dirichlet .and. bc1_max == sll_p_dirichlet .and. &
        bc2_min == sll_p_periodic  .and. bc2_max == sll_p_periodic) then

  d = 0
  do j = 1, nb_spl_y
    do i = 1, nb_spl_x
      a = i + nb_spl_x*(j-1)
      if (i == 1 .OR. i == nb_spl_x) then
        global_indices(a) = 0
      else
        if (j /= nb_spl_y) then
          if (global_indices(a) == 0) then
            d = d + 1
            global_indices(a) = d
          end if
          if ( 1<=j .and. j<=spline_degree2 ) then
            l = (nb_spl_y - spline_degree2) * nb_spl_x
            global_indices(a+l) = global_indices(a)
          end if
        end if
      end if
    end do
  end do

else if(bc1_min == sll_p_dirichlet .and. bc1_max==sll_p_dirichlet .and. &
        bc2_min == sll_p_dirichlet .and. bc2_max==sll_p_dirichlet) then

  !PN : Since points on boundary are not computed, 
  !PN : i guess Dirichlet is homogeneous only
  !PN : Some values of global_indices are not set  
  d = 0
  do j = 2, nb_spl_y-1
    do i = 2, nb_spl_x-1
      a = i + nb_spl_x*(j-1)
      d = d + 1
      global_indices(a) = d
    end do
  end do

end if
  
do e = 1, num_cells1*num_cells2
  do b = 1, (spline_degree1+1)*(spline_degree2+1)
    local_to_global_indices(b,e) = global_indices(local_indices(b,e))
  end do
end do

deallocate(local_indices)
deallocate(global_indices)

end subroutine initconnectivity



  subroutine initconnectivity_new( &
    num_cells1, &
    num_cells2, &
    spline_degree1, &
    spline_degree2, &
    bc1_min, &
    bc1_max, &
    bc2_min, &
    bc2_max, &
    local_to_global_indices )

    sll_int32, intent(in) :: num_cells1
    sll_int32, intent(in) :: num_cells2
    sll_int32, intent(in) :: spline_degree1
    sll_int32, intent(in) :: spline_degree2
    sll_int32, intent(in) :: bc1_min
    sll_int32, intent(in) :: bc1_max
    sll_int32, intent(in) :: bc2_min
    sll_int32, intent(in) :: bc2_max
    sll_int32, dimension(:,:), intent(out) :: local_to_global_indices
    sll_int32, dimension(:,:), allocatable :: local_indices
    sll_int32, dimension(:)  , allocatable :: global_indices
    
    sll_int32, dimension(:), allocatable :: bc_indices1
    sll_int32, dimension(:), allocatable :: bc_indices2
    sll_int32 :: ierr
    
    SLL_ALLOCATE(bc_indices1(num_cells1+spline_degree1),ierr)
    SLL_ALLOCATE(bc_indices2(num_cells2+spline_degree2),ierr)
    SLL_ALLOCATE(local_indices( 1:(spline_degree1+1)*(spline_degree2+1),1:(num_cells1*num_cells2)),ierr)
    SLL_ALLOCATE(global_indices((num_cells1+spline_degree1)*(num_cells2+spline_degree2)),ierr)
    
    
    call compute_local_splines_indices( &
      num_cells1, &
      num_cells2, &
      spline_degree1, &
      spline_degree2, &
      local_indices)
   
    call compute_bc_indices( &
      num_cells1, &
      spline_degree1, &
      bc1_min, &
      bc1_max, &
      bc_indices1)

    call compute_bc_indices( &
      num_cells2, &
      spline_degree2, &
      bc2_min, &
      bc2_max, &
      bc_indices2)


    call compute_global_splines_indices( &
      num_cells1+spline_degree1, &
      num_cells2+spline_degree2, &
      bc_indices1, &
      bc_indices2, &
      global_indices)
    
    call compute_local_to_global_splines_indices( &
      num_cells1, &
      num_cells2, &
      spline_degree1,&
      spline_degree2,&
      local_to_global_indices, &
      local_spline_indices = local_indices, &
      global_spline_indices = global_indices)  
    
  end subroutine initconnectivity_new


  subroutine compute_local_splines_indices( &
    num_cells1, &
    num_cells2, &
    spline_degree1, &
    spline_degree2, &
    local_spline_indices)
  
    sll_int32, intent(in) :: num_cells1
    sll_int32, intent(in) :: num_cells2
    sll_int32, intent(in) :: spline_degree1
    sll_int32, intent(in) :: spline_degree2
    sll_int32, dimension(:,:), intent(out) :: local_spline_indices
    sll_int32 :: i
    sll_int32 :: j
    sll_int32 :: iloc
    sll_int32 :: jloc
    sll_int32 :: maille
    sll_int32 :: li_B
    sll_int32 :: li_Bloc
  
    do j = 1, num_cells2    
      do i =1, num_cells1  
        maille = num_cells1*(j-1) + i
        do jloc = 0 , spline_degree2
          do iloc = 0, spline_degree1
            li_B = (j-1+jloc)*(num_cells1+spline_degree1) + i +iloc
            li_Bloc = jloc * (spline_degree1 + 1) + iloc + 1
            local_spline_indices(li_Bloc, maille) = li_B
          end do
        end do
      end do
    end do  
  
  end subroutine compute_local_splines_indices

  subroutine compute_bc_indices(num_cells,spline_degree,bc_min,bc_max,index)
    sll_int32, intent(in) :: num_cells
    sll_int32, intent(in) :: spline_degree
    sll_int32, intent(in) :: bc_min
    sll_int32, intent(in) :: bc_max
    sll_int32, dimension(:), intent(out) :: index
  
    sll_int32 :: i
    
    do i = 1, num_cells+spline_degree
      index(i) = i
    enddo
  
    !eta_min boundary correction  
    select case (bc_min)
      case (sll_p_dirichlet)
        index(1) = 0
      case default     
    end select
      
    !eta_max boundary correction
    select case (bc_max)
      case (sll_p_dirichlet)
        index(num_cells+spline_degree) = 0
      case (sll_p_periodic)
        do i=1,spline_degree
          index(num_cells+i)=index(i)
        end do  
      case default     
    end select
    
  end subroutine compute_bc_indices

subroutine compute_global_splines_indices( &
  num_spl1, &
  num_spl2, &
  bc_indices1, &
  bc_indices2, &
  global_spline_indices)
  
  sll_int32, intent(in) :: num_spl1
  sll_int32, intent(in) :: num_spl2
  sll_int32, dimension(:), intent(in) :: bc_indices1
  sll_int32, dimension(:), intent(in) :: bc_indices2
  sll_int32, dimension(:), intent(out) :: global_spline_indices
  
  sll_int32 :: i
  sll_int32 :: j
  
  sll_int32 :: ii
  sll_int32 :: jj
  
  
  sll_int32 :: li_A
  sll_int32 :: count
  
  count = 0
  
  global_spline_indices(1:num_spl1*num_spl2) = 0
  
  do j=1,num_spl2
    jj = bc_indices2(j)
    if(jj /=0) then
      do i=1,num_spl1
        ii = bc_indices1(i)
        if(ii /=0) then
          li_A = ii + num_spl1*(jj-1)
          if(global_spline_indices(li_A)==0)then
            count = count+1
            global_spline_indices(li_A) = count
          else  
            global_spline_indices(i+num_spl1*(j-1)) = global_spline_indices(li_A)  
          endif  
        endif  
      enddo
    endif  
  enddo  

end subroutine compute_global_splines_indices



subroutine compute_local_to_global_splines_indices( &
  num_cells1,&
  num_cells2,&
  spline_degree1,&
  spline_degree2,&
  local_to_global_spline_indices, &
  local_spline_indices,&
  global_spline_indices, &
  bc_indices1, &
  bc_indices2)
  
  sll_int32, intent(in) :: num_cells1
  sll_int32, intent(in) :: num_cells2
  sll_int32, intent(in) :: spline_degree1
  sll_int32, intent(in) :: spline_degree2
  sll_int32, dimension(:,:), intent(out) :: local_to_global_spline_indices
  sll_int32, dimension(:,:), intent(in), optional :: local_spline_indices
  sll_int32, dimension(:), intent(in), optional :: global_spline_indices
  sll_int32, dimension(:), intent(in), optional :: bc_indices1
  sll_int32, dimension(:), intent(in), optional :: bc_indices2
  
  sll_int32 :: i
  sll_int32 :: iloc
  
  if(present(local_spline_indices) &
    .and.present(global_spline_indices) )then
  
    do i= 1, num_cells1*num_cells2
      do iloc = 1, (spline_degree1+1)*(spline_degree2+1)
        local_to_global_spline_indices(iloc, i) = &
          global_spline_indices(local_spline_indices(iloc, i))
      enddo
    enddo
  else
    if(.not.(present(bc_indices1) &
      .and.present(bc_indices2)) )then
      
      SLL_ERROR('compute_local_to_global_splines_indices&
      &', 'bc_indices1 and bc_indices2 should be present&
      & as it is not the case for local/global_spline_indices')
      
    else  
      call compute_local_to_global_splines_indices_from_bc_indices( &
        num_cells1,&
        num_cells2,&
        spline_degree1,&
        spline_degree2,&
        bc_indices1, &
        bc_indices2, &
        local_to_global_spline_indices)
    endif  
  endif  
end subroutine compute_local_to_global_splines_indices


subroutine compute_local_to_global_splines_indices_from_bc_indices( &
  num_cells1,&
  num_cells2,&
  spline_degree1,&
  spline_degree2,&
  bc_indices1, &
  bc_indices2, &
  local_to_global_spline_indices)
  
  sll_int32, intent(in) :: num_cells1
  sll_int32, intent(in) :: num_cells2
  sll_int32, intent(in) :: spline_degree1
  sll_int32, intent(in) :: spline_degree2
  sll_int32, dimension(:), intent(in) :: bc_indices1
  sll_int32, dimension(:), intent(in) :: bc_indices2
  sll_int32, dimension(:,:), intent(out) :: local_to_global_spline_indices


    sll_int32, dimension(:,:), allocatable :: local_indices
    sll_int32, dimension(:)  , allocatable :: global_indices
    
    sll_int32 :: ierr
    
    SLL_ALLOCATE(local_indices( 1:(spline_degree1+1)*(spline_degree2+1),1:(num_cells1*num_cells2)),ierr)
    SLL_ALLOCATE(global_indices((num_cells1+spline_degree1)*(num_cells2+spline_degree2)),ierr)
    
    
    call compute_local_splines_indices( &
      num_cells1, &
      num_cells2, &
      spline_degree1, &
      spline_degree2, &
      local_indices)
   

    call compute_global_splines_indices( &
      num_cells1+spline_degree1, &
      num_cells2+spline_degree2, &
      bc_indices1, &
      bc_indices2, &
      global_indices)
    
    call compute_local_to_global_splines_indices( &
      num_cells1, &
      num_cells2, &
      spline_degree1,&
      spline_degree2,&
      local_to_global_spline_indices, &
      local_spline_indices = local_indices, &
      global_spline_indices = global_indices)  


end subroutine compute_local_to_global_splines_indices_from_bc_indices



!> @brief
!> Assembling knots array
!> @details
!> it is intentional that eta_min not used, one intends to consider
!> only the [0,1] interval...
subroutine initialize_knots( spline_degree, &
&                            num_cells,     &
&                            eta_min,       & 
&                            eta_max,       &
&                            bc_left,       &
&                            bc_right,      &
&                            knots)
    
sll_int32,  intent(in)  :: spline_degree 
sll_int32,  intent(in)  :: num_cells
sll_real64, intent(in)  :: eta_min
sll_real64, intent(in)  :: eta_max
sll_int32,  intent(in)  :: bc_left
sll_int32,  intent(in)  :: bc_right
sll_real64, intent(out) :: knots(:)
sll_int32               :: i
sll_real64              :: eta
sll_real64              :: delta

delta = (eta_max - eta_min)/num_cells

if (bc_left==sll_p_periodic .and. bc_right==sll_p_periodic) then 

  do i = -spline_degree, spline_degree+1
    knots(i+spline_degree+1) = delta*i 
  end do

else if (bc_left==sll_p_dirichlet .and. bc_right==sll_p_dirichlet) then 

  do i = 1, spline_degree + 1
    knots(i) = eta_min
  enddo
  eta = eta_min
  do i = spline_degree+2, num_cells+1+spline_degree
    eta = eta + delta
    knots(i) = eta
  enddo
  do i = num_cells+spline_degree+1, num_cells+1+2*spline_degree
    knots(i) = eta_max
  enddo

end if

end subroutine initialize_knots



!next routine is for doing the general elliptic solver
!without scalar fields
!the aim is to work only with arrays
! commented for the moment
!subroutine solve_general_coordinates_elliptic_eq_discrete_vec( es, rho, phi)
!
!  class(sll_t_general_coordinate_elliptic_solver), intent(inout) :: es
!  sll_real64, dimension(:,:), intent(in) :: rho
!  sll_int32  :: i
!  sll_int32  :: j
!  sll_real64, dimension(:,:), pointer :: coeff_rho
!
!  call interp_2d%compute_interpolants(rho)  
!  coeff_rho => phi%interp_2d%get_coefficients()
!
!  do j=1,es%num_cells2+1
!    do i=1,es%num_cells1+1
!      es%rho_coeff_1d(i+(es%num_cells1+1)*(j-1)) = rho(i,j)
!    end do
!  end do
!
!  call sll_s_mult_csr_matrix_vector(es%csr_mat_source,es%rho_coeff_1d,es%rho_vec)
!
!  if(es%perper) then
!    es%rho_vec = es%rho_vec - sum(es%rho_vec)/es%intjac*es%masse
!  end if
!
!  call solve_linear_system(es)
!
!  call phi%interp_2d%set_coefficients(es%phi_vec(1:es%total_num_splines1*es%total_num_splines2))
!  
!
!end subroutine solve_general_coordinates_elliptic_eq_discrete_vec
!

  
  
!!PN DEFINED BUT NOT USED
!!> @details
!!> CSR_MAT*phi = rho_vec is the linear system to be solved. The solution
!!> is given in terms of the spline coefficients that represent phi.
!subroutine solve_linear_system( es )
!
!class(sll_t_general_coordinate_elliptic_solver) :: es
!integer                                   :: elt
!integer                                   :: i,j,k
!character(len=*), parameter               :: as_file  = 'rho'
!character(len=*), parameter               :: as_file1 = 'phi'
!character(len=*), parameter               :: as_file2 = 'mat'
!sll_int32 :: bc1_min
!sll_int32 :: bc1_max
!sll_int32 :: bc2_min
!sll_int32 :: bc2_max
!
!end subroutine


subroutine compute_non_zero_splines_and_deriv_at_cell_points( &
  num_cells, &
  degree, &
  knots, &
  cell_points, & 
  num_cell_points, &
  v_splines )
  sll_int32, intent(in) :: num_cells
  sll_int32, intent(in) :: degree
  sll_real64, dimension(:), intent(in) :: knots
  sll_real64, dimension(:), intent(in) :: cell_points
  sll_int32, intent(in) :: num_cell_points
  sll_real64, dimension(:,:,:,:), intent(out) :: v_splines
  
  sll_int32 :: ierr
  sll_int32 :: i
  sll_int32 :: j
  sll_real64 :: eta_min_loc
  sll_real64 :: eta_max_loc
  sll_real64 :: delta_eta_loc
  sll_real64 :: eta
  sll_real64, dimension(:,:), allocatable :: dbiatx
  if(num_cells<1)then
    print *,'#num_cells=',num_cells
    flush( output_unit )
    SLL_ERROR('compute_non_zero_splines_and_deriv_at_cell_points&
    &','Problem with num_cells')
  endif
  if(degree<0)then
    print *,'#degree=',degree
    flush( output_unit )
    SLL_ERROR('compute_non_zero_splines_and_deriv_at_cell_points&
    &','Problem with degree')
  endif
  if(num_cell_points<1)then
    print *,'#num_cell_points=',num_cell_points
    flush( output_unit )
    SLL_ERROR('compute_non_zero_splines_and_deriv_at_cell_points&
    &','Problem with num_cell_points')
  endif
  if(minval(cell_points)<0.)then
    print *,'#minval(cell_points)=',minval(cell_points)
    flush( output_unit )
    SLL_ERROR('compute_non_zero_splines_and_deriv_at_cell_points&
    &','Problem with minval(cell_points)')
  endif
  if(maxval(cell_points)>=1.)then
    print *,'#maxval(cell_points)=',maxval(cell_points)
    flush( output_unit )
    SLL_ERROR('compute_non_zero_splines_and_deriv_at_cell_points&
    &','Problem with maxval(cell_points)')
  endif
  
  if(size(knots)<num_cells+2*degree+1)then
    print *,'#size(knots)=',size(knots)
    flush( output_unit )
    SLL_ERROR('compute_non_zero_splines_and_deriv_at_cell_points&
    &','Problem with knots')
  endif
  if(size(v_splines,1)<2)then
    print *,'#size(v_splines,1)=',size(v_splines,1),2
    flush( output_unit )
    SLL_ERROR('compute_non_zero_splines_and_deriv_at_cell_points&
    &','Problem with size of v_splines')
  endif
  if(size(v_splines,2)<degree+1)then
    print *,'#size(v_splines,2)=',size(v_splines,2),degree+1
    flush( output_unit )
    SLL_ERROR('compute_non_zero_splines_and_deriv_at_cell_points&
    &','Problem with size of v_splines')
  endif
  if(size(v_splines,3)<num_cell_points)then
    print *,'#size(v_splines,3)=',size(v_splines,3),num_cell_points
    flush( output_unit )
    SLL_ERROR('compute_non_zero_splines_and_deriv_at_cell_points&
    &','Problem with size of v_splines')
  endif
  if(size(v_splines,4)<num_cells)then
    print *,'#size(v_splines,4)=',size(v_splines,4),num_cells
    flush( output_unit )
    SLL_ERROR('compute_non_zero_splines_and_deriv_at_cell_points&
    &','Problem with size of v_splines')
  endif
  
  SLL_ALLOCATE(dbiatx(degree+1,2),ierr)
  
  do i=1,num_cells
    eta_min_loc = knots(i+degree)
    eta_max_loc = knots(i+1+degree)
    delta_eta_loc = eta_max_loc-eta_min_loc
    do j=1,num_cell_points
      eta = eta_min_loc+cell_points(j)*delta_eta_loc
      !print *,'#eta=',eta_min_loc,cell_points(j),delta_eta_loc
      !flush( output_unit )
      call compute_b_spline_and_deriv_at_x( &
        knots, &
        i+degree, &
        eta, &
        degree, &
        dbiatx)
      v_splines(1,1:degree+1,j,i) = dbiatx(:,1)
      v_splines(2,1:degree+1,j,i) = dbiatx(:,2)
      !print *,  v_splines(1,1:degree+1,j,i)
      !flush( output_unit )
      !stop
    enddo  
  enddo
  
  
end subroutine compute_non_zero_splines_and_deriv_at_cell_points





subroutine compute_non_zero_splines_at_cell_points( &
  num_cells, &
  degree, &
  knots, &
  cell_points, & 
  num_cell_points, &
  v_splines )
  sll_int32, intent(in) :: num_cells
  sll_int32, intent(in) :: degree
  sll_real64, dimension(:), intent(in) :: knots
  sll_real64, dimension(:), intent(in) :: cell_points
  sll_int32, intent(in) :: num_cell_points
  sll_real64, dimension(:,:,:), intent(out) :: v_splines
  
  sll_int32 :: ierr
  sll_int32 :: i
  sll_int32 :: j
  sll_real64 :: eta_min_loc
  sll_real64 :: eta_max_loc
  sll_real64 :: delta_eta_loc
  sll_real64 :: eta
  sll_real64, dimension(:), allocatable :: dbiatx
  if(num_cells<1)then
    print *,'#num_cells=',num_cells
    flush( output_unit )
    SLL_ERROR('compute_non_zero_splines_and_deriv_at_cell_points&
    &','Problem with num_cells')
  endif
  if(degree<0)then
    print *,'#degree=',degree
    flush( output_unit )
    SLL_ERROR('compute_non_zero_splines_and_deriv_at_cell_points&
    &','Problem with degree')
  endif
  if(num_cell_points<1)then
    print *,'#num_cell_points=',num_cell_points
    flush( output_unit )
    SLL_ERROR('compute_non_zero_splines_and_deriv_at_cell_points&
    &','Problem with num_cell_points')
  endif
  if(minval(cell_points)<0.)then
    print *,'#minval(cell_points)=',minval(cell_points)
    flush( output_unit )
    SLL_ERROR('compute_non_zero_splines_and_deriv_at_cell_points&
    &','Problem with minval(cell_points)')
  endif
  if(maxval(cell_points)>=1.)then
    print *,'#maxval(cell_points)=',maxval(cell_points)
    flush( output_unit )
    SLL_ERROR('compute_non_zero_splines_and_deriv_at_cell_points&
    &','Problem with maxval(cell_points)')
  endif
  
  if(size(knots)<num_cells+2*degree+1)then
    print *,'#size(knots)=',size(knots)
    flush( output_unit )
    SLL_ERROR('compute_non_zero_splines_and_deriv_at_cell_points&
    &','Problem with knots')
  endif
  if(size(v_splines,1)<degree+1)then
    print *,'#size(v_splines,1)=',size(v_splines,1),degree+1
    flush( output_unit )
    SLL_ERROR('compute_non_zero_splines_and_deriv_at_cell_points&
    &','Problem with size of v_splines')
  endif
  if(size(v_splines,2)<num_cell_points)then
    print *,'#size(v_splines,2)=',size(v_splines,2),num_cell_points
    flush( output_unit )
    SLL_ERROR('compute_non_zero_splines_and_deriv_at_cell_points&
    &','Problem with size of v_splines')
  endif
  if(size(v_splines,3)<num_cells)then
    print *,'#size(v_splines,3)=',size(v_splines,3),num_cells
    flush( output_unit )
    SLL_ERROR('compute_non_zero_splines_and_deriv_at_cell_points&
    &','Problem with size of v_splines')
  endif
  
  SLL_ALLOCATE(dbiatx(degree+1),ierr)
  
  do i=1,num_cells
    eta_min_loc = knots(i+degree)
    eta_max_loc = knots(i+1+degree)
    delta_eta_loc = eta_max_loc-eta_min_loc
    do j=1,num_cell_points
      eta = eta_min_loc+cell_points(j)*delta_eta_loc
      !print *,'#eta=',eta_min_loc,cell_points(j),delta_eta_loc
      !flush( output_unit )
      call compute_b_spline_at_x( &
        knots, &
        i+degree, &
        eta, &
        degree, &
        dbiatx)
      v_splines(1:degree+1,j,i) = dbiatx
      !print *,  v_splines(1,1:degree+1,j,i)
      !flush( output_unit )
      !stop
    enddo  
  enddo
  
  
end subroutine compute_non_zero_splines_at_cell_points



subroutine compute_global_knots( &
  bc_min, &
  bc_max, &
  num_cells, &
  eta_min, &
  eta_max, &
  spline_degree, &
  knots)
  sll_int32, intent(in) :: bc_min
  sll_int32, intent(in) :: bc_max
  sll_int32, intent(in) :: num_cells
  sll_real64, intent(in) :: eta_min
  sll_real64, intent(in) :: eta_max
  sll_int32, intent(in) :: spline_degree
  sll_real64, dimension(:), intent(out) :: knots
  sll_real64 :: delta
  sll_int32 :: i
  
  delta = (eta_max - eta_min)/real(num_cells,f64)
  
  
  !repeated left point
  do i = 1, spline_degree + 1
    knots(i) = eta_min
  enddo
  !interior points
  do i = spline_degree+2, num_cells+spline_degree
    knots(i) = eta_min+real(i-spline_degree-1,f64)*delta
  enddo
  !repeated right point
  do i = num_cells+spline_degree+1, num_cells+1+2*spline_degree
    knots(i) = eta_max
  enddo


  if((bc_min==sll_p_periodic).and.(bc_max==sll_p_periodic))then
!    do i = 1, spline_degree+num_cells+1+2*spline_degree
    do i = 1, num_cells+1+2*spline_degree
      knots(i) = eta_min+real(i-spline_degree-1,f64)*delta
    enddo
  endif

end subroutine compute_global_knots

!knots(cell+1-degree)<= ..<= knots(cell)<=x <knots(cell+1) <=..<=knots(cell+degree)
subroutine compute_b_spline_at_x( &
  knots, &
  cell, &
  x, &
  degree, &
  out)
  sll_real64, dimension(:), intent(in) :: knots
  sll_int32, intent(in) :: cell
  sll_real64, intent(in) :: x
  sll_int32, intent(in) :: degree
  sll_real64, dimension(:), intent(out) :: out
  
  sll_real64 :: tmp1
  sll_real64 :: tmp2
  sll_int32 :: ell
  sll_int32 :: k
  !sll_real64 :: res
  
  out(1) = 1._f64
  do ell=1,degree
    tmp1 = (x-knots(cell+1-ell))/(knots(cell+1)-knots(cell+1-ell))*out(1)
    out(1) = out(1) -tmp1
    do k=2,ell
      tmp2 = (x-knots(cell+k-ell))/(knots(cell+k)-knots(cell+k-ell))*out(k)
      out(k) = out(k)+tmp1-tmp2
      tmp1 = tmp2
    enddo
    out(ell+1) = tmp1
  enddo
    
  
end subroutine compute_b_spline_at_x

!!PN DEFINED BUT NOT USED
!function check_compute_b_spline_at_x( &
!  knots, &
!  cell, &
!  x, &
!  degree &
!  ) result(res)
!  sll_real64, dimension(:), intent(in) :: knots
!  sll_int32, intent(in) :: cell
!  sll_real64, intent(in) :: x
!  sll_int32, intent(in) :: degree
!  sll_real64 :: res
!  type(sll_t_deboor_type)                       :: db
!  sll_real64, dimension(:), allocatable :: splines1
!  sll_real64, dimension(:), allocatable :: splines2
!  sll_real64, dimension(:,:), allocatable :: work
!  sll_real64, dimension(:,:), allocatable :: dbs1
!  sll_int32 :: ierr
!  sll_int32 :: i
!  sll_real64 :: xx
!  
!  res = 0._f64
!  SLL_ALLOCATE(splines1(degree+1),ierr)
!  SLL_ALLOCATE(splines2(degree+1),ierr)
!  SLL_ALLOCATE(work(degree+1,degree+1),ierr)
!  SLL_ALLOCATE(dbs1(degree+1,2),ierr)
!
!!  call bsplvb_db2( &
!!    knots, &
!!    degree+1, &
!!    cell, &
!!    x, &
!!    degree+1, &
!!    splines2)
!  xx = x
!  call sll_s_bsplvd( &
!    db, &
!    knots(1:cell+degree+1), &
!    degree+1, &
!    xx, &
!    cell, &
!    work, &
!    dbs1, &
!    2)
!    !call sll_s_bsplvd(es%db,es%knots1,spline_degree1+1,gspl1,ispl1,work1,dbs1,2)
!  splines2(:) = dbs1(:,1)
!
!  call compute_b_spline_at_x( &
!    knots, &
!    cell, &
!    xx, &
!    degree, &
!    splines1)
!  res=maxval(abs(splines1-splines2))
!  
!  if(res>1.e-14)then
!    print *,'In check_compute_b_spline_at_x:'
!    print *,'#res=',res
!    print *,'#cell=',cell
!    print *,'#degree=',degree
!    print *,'#xx=',xx
!    !print *,'#knots=',knots
!    
!    do i=1,degree+1
!      print *,splines1(i),splines2(i)
!    enddo
!    flush( output_unit )
!  endif
!end function check_compute_b_spline_at_x


!knots(cell+1-degree)<= ..<= knots(cell)<=x <knots(cell+1) <=..<=knots(cell+degree)
subroutine compute_b_spline_and_deriv_at_x( &
  knots, &
  cell, &
  x, &
  degree, &
  out)
  sll_real64, dimension(:), intent(in) :: knots
  sll_int32, intent(in) :: cell
  sll_real64, intent(in) :: x
  sll_int32, intent(in) :: degree
  sll_real64, dimension(:,:), intent(out) :: out
  
  sll_real64 :: tmp1
  sll_real64 :: tmp2
  sll_int32 :: ell
  sll_int32 :: k
  
  out(1,1) = 1._f64
  do ell=1,degree
    tmp1 = (x-knots(cell+1-ell))/(knots(cell+1)-knots(cell+1-ell))*out(1,1)
    out(1,1) = out(1,1) -tmp1
    do k=2,ell
      tmp2 = (x-knots(cell+k-ell))/(knots(cell+k)-knots(cell+k-ell))*out(k,1)
      out(k,1) = out(k,1)+tmp1-tmp2
      tmp1 = tmp2
    enddo
    out(ell+1,1) = tmp1
    if(ell==degree-1)then
      !compute the derivatives
      tmp1 = real(degree,f64)/(knots(cell+1)-knots(cell+1-degree))*out(1,1)
      out(1,2) = -tmp1
      do k=2,degree
        out(k,2) = tmp1
        tmp1 = real(degree,f64)/(knots(cell+k)-knots(cell+k-degree))*out(k,1)
        out(k,2) = out(k,2)-tmp1
      enddo
      out(degree+1,2) = tmp1
    endif
  enddo
   
  
end subroutine compute_b_spline_and_deriv_at_x


subroutine sll_s_initialize_general_elliptic_solver_prototype( &
       es,                                     &
       spline_degree1,                         &
       spline_degree2,                         &
       num_cells1,                             &
       num_cells2,                             &
       quadrature_type1,                       &
       quadrature_type2,                       &
       bc1_min,                                &
       bc1_max,                                &
       bc2_min,                                &
       bc2_max,                                &
       eta1_min,                               &
       eta1_max,                               &
       eta2_min,                               &
       eta2_max,                               &
       precompute_rhs,                         &
       rhs_bc1,                                &
       rhs_bc2,                                &
       use_cubic_splines,                      &
       rho_degree1,                            &
       rho_degree2,                            &
       with_constraint,                        &
       zero_mean)
    
type(sll_t_general_coordinate_elliptic_solver), intent(out) :: es

sll_int32,  intent(in) :: spline_degree1
sll_int32,  intent(in) :: spline_degree2
sll_int32,  intent(in) :: num_cells1
sll_int32,  intent(in) :: num_cells2
sll_int32,  intent(in) :: bc1_min
sll_int32,  intent(in) :: bc1_max
sll_int32,  intent(in) :: bc2_min
sll_int32,  intent(in) :: bc2_max
sll_int32,  intent(in) :: quadrature_type1
sll_int32,  intent(in) :: quadrature_type2
sll_real64, intent(in) :: eta1_min
sll_real64, intent(in) :: eta1_max
sll_real64, intent(in) :: eta2_min
sll_real64, intent(in) :: eta2_max
logical,    intent(in), optional :: precompute_rhs
sll_int32,  intent(in), optional :: rhs_bc1
sll_int32,  intent(in), optional :: rhs_bc2
logical, intent(in), optional :: use_cubic_splines
sll_int32, intent(in), optional :: rho_degree1
sll_int32, intent(in), optional :: rho_degree2
logical, intent(in), optional :: with_constraint
logical, intent(in), optional :: zero_mean
sll_int32 :: num_splines1
sll_int32 :: num_splines2
sll_int32 :: vec_sz ! for rho_vec allocation
sll_int32 :: ierr
sll_int32 :: solution_size
sll_int32 :: dim1, dim2
sll_int32 :: num_pts_g1, num_pts_g2




sll_int32 :: vec_sz1
sll_int32 :: vec_sz2
sll_int32, dimension(:), allocatable :: bc_indices1
sll_int32, dimension(:), allocatable :: bc_indices2
sll_int32, dimension(:), allocatable :: bc_rho_indices1
sll_int32, dimension(:), allocatable :: bc_rho_indices2
sll_int32 :: i


if(present(precompute_rhs))then
  es%precompute_rhs = precompute_rhs
else
  es%precompute_rhs = .false.
endif

if(present(use_cubic_splines))then
  es%use_cubic_splines = use_cubic_splines
else
  es%use_cubic_splines = .false.
endif



if(es%use_cubic_splines)then
  if(present(rhs_bc1) .and. present(rhs_bc2))then
  es%cubic_spline => sll_f_new_cubic_spline_2d( &
    num_cells1+1, &
    num_cells2+1, &
    eta1_min, &
    eta1_max, &
    eta2_min, &
    eta2_max, &
    rhs_bc1, &
    rhs_bc2)
  else
    SLL_ERROR('sll_s_initialize_general_elliptic_solver_prototype&
    &','provide rhs_bc1 and rhs_bc2')
  endif  
endif

if(present(rho_degree1))then
  es%rho_degree1 = rho_degree1
else
  es%rho_degree1 = 3   
endif

if(present(rho_degree2))then
  es%rho_degree2 = rho_degree2
else
  es%rho_degree2 = 3   
endif

if(es%use_cubic_splines)then
  es%rho_degree1 = 3
  es%rho_degree2 = 3
endif

es%total_num_splines_loc = (spline_degree1+1)*(spline_degree2+1)
! The total number of splines in a single direction is given by
! num_cells + spline_degree
num_splines1 = num_cells1 + spline_degree1
num_splines2 = num_cells2 + spline_degree2
  
dim1 = (spline_degree1+1)*(spline_degree2+1)
dim2 = (num_cells1*num_cells2)
SLL_ALLOCATE(es%local_to_global_indices(1:dim1,1:dim2),ierr)

if(es%precompute_rhs)then
  dim1 = (spline_degree1+1)*(spline_degree2+1)
  !source_bis is for rows
  SLL_ALLOCATE(es%local_to_global_indices_source_bis(1:dim1,1:dim2),ierr)
  !source is for columns
  dim1 = (es%rho_degree1+1)*(es%rho_degree2+1)
  SLL_ALLOCATE(es%local_to_global_indices_source(1:dim1,1:dim2),ierr)
  dim1 = num_cells1 + es%rho_degree1
  dim2 = num_cells2 + es%rho_degree2
  SLL_ALLOCATE(es%rho_coeff_1d(dim1*dim2),ierr)
  SLL_ALLOCATE(bc_rho_indices1(num_cells1+es%rho_degree1),ierr)
  SLL_ALLOCATE(bc_rho_indices2(num_cells2+es%rho_degree2),ierr)
endif


SLL_ALLOCATE(es%knots1(num_cells1 + 2*spline_degree1+1),ierr)
SLL_ALLOCATE(es%knots2(num_cells2 + 2*spline_degree2+1),ierr)
SLL_ALLOCATE(es%knots1_rho(num_cells1 + 2*es%rho_degree1+1),ierr)
SLL_ALLOCATE(es%knots2_rho(num_cells2 + 2*es%rho_degree2+1),ierr)
SLL_ALLOCATE(bc_indices1(num_cells1+spline_degree1),ierr)
SLL_ALLOCATE(bc_indices2(num_cells2+spline_degree2),ierr)


call compute_bc_indices( &
  num_cells1, &
  spline_degree1, &
  bc1_min, &
  bc1_max, &
  bc_indices1 )
call compute_bc_indices( &
  num_cells2, &
  spline_degree2, &
  bc2_min, &
  bc2_max, &
  bc_indices2 )


call compute_local_to_global_splines_indices( &
  num_cells1, &
  num_cells2, &
  spline_degree1, &
  spline_degree2, &
  es%local_to_global_indices, &
  bc_indices1 = bc_indices1, &
  bc_indices2 = bc_indices2 )

if(es%precompute_rhs)then
print *,'#begin connectivity rho'
flush( output_unit )

call compute_index( &
  num_cells1, &
  spline_degree1, &
  bc1_min, &
  bc1_max, &
  bc_indices1 )
call compute_index( &
  num_cells2, &
  spline_degree2, &
  bc2_min, &
  bc2_max, &
  bc_indices2 )
call compute_local_to_global_splines_indices( &
  num_cells1, &
  num_cells2, &
  spline_degree1, &
  spline_degree2, &
  es%local_to_global_indices_source_bis, &
  bc_indices1 = bc_indices1, &
  bc_indices2 = bc_indices2 )

do i=1, num_cells1+es%rho_degree1
  bc_rho_indices1(i) = i  
enddo  
do i=1, num_cells2+es%rho_degree2
  bc_rho_indices2(i) = i  
enddo  
  
call compute_local_to_global_splines_indices( &
  num_cells1, &
  num_cells2, &
  es%rho_degree1, &
  es%rho_degree2, &
  es%local_to_global_indices_source, &
  bc_indices1 = bc_rho_indices1, &
  bc_indices2 = bc_rho_indices2 )
print *,'#end connectivity rho'
flush( output_unit )
endif



print *,'#initconnectivity_new done'
flush( output_unit )



! This should be changed to verify that the passed BC's are part of the
! recognized list described in sll_m_boundary_condition_descriptors...

es%bc1_min        = bc1_min
es%bc1_max        = bc1_max
es%bc2_min        = bc2_min
es%bc2_max        = bc2_max
es%spline_degree1 = spline_degree1
es%spline_degree2 = spline_degree2
es%num_cells1     = num_cells1
es%num_cells2     = num_cells2
es%delta_eta1     = (eta1_max-eta1_min)/num_cells1
es%delta_eta2     = (eta2_max-eta2_min)/num_cells2
es%eta1_min       = eta1_min
es%eta2_min       = eta2_min

! Allocate and fill the gauss points/weights information.
select case(quadrature_type1)
case (sll_p_es_gauss_legendre)
  SLL_ALLOCATE(es%gauss_pts1(2,spline_degree1+2),ierr)
  es%gauss_pts1 = sll_f_gauss_legendre_points_and_weights(spline_degree1+2)
case (sll_p_es_gauss_lobatto)
  SLL_ALLOCATE(es%gauss_pts1(2,spline_degree1+2),ierr)
  es%gauss_pts1 = sll_f_gauss_lobatto_points_and_weights(spline_degree1+2)
case default
  SLL_ERROR('initialize_general_elliptic_solver&
  &','unknown type of gauss points in the direction 1')
end select
   
select case(quadrature_type2)
case (sll_p_es_gauss_legendre)
  SLL_ALLOCATE(es%gauss_pts2(2,spline_degree2+2),ierr)
  es%gauss_pts2(:,:) = sll_f_gauss_legendre_points_and_weights(spline_degree2+2)
case (sll_p_es_gauss_lobatto)
  SLL_ALLOCATE(es%gauss_pts2(2,spline_degree2+2),ierr)
  es%gauss_pts2(:,:) = sll_f_gauss_lobatto_points_and_weights(spline_degree2+2)
case default
  SLL_ERROR('initialize_general_elliptic_solver&
  &','unknown type of gauss points in the direction 2')
end select

!PN : Gauss points positions and weights are computed in [-1:1] interval
!PN : We need to rescale them.
es%gauss_pts1(1,:) = 0.5_f64*es%delta_eta1*(es%gauss_pts1(1,:)+1.0_f64)
es%gauss_pts1(2,:) = 0.5_f64*es%delta_eta1*es%gauss_pts1(2,:)
es%gauss_pts2(1,:) = 0.5_f64*es%delta_eta2*(es%gauss_pts2(1,:)+1.0_f64)
es%gauss_pts2(2,:) = 0.5_f64*es%delta_eta2*es%gauss_pts2(2,:)


es%perper  = .false. 
if( bc1_min == sll_p_periodic .and. bc1_max == sll_p_periodic .and. &
    bc2_min == sll_p_periodic .and. bc2_max == sll_p_periodic ) then
  es%perper = .true.
else
  es%perper  = .false.  
endif

if(present(with_constraint))then
  es%with_constraint = with_constraint
else
  es%with_constraint = es%perper
endif
if(present(zero_mean))then
  es%zero_mean = zero_mean
else
  es%zero_mean = es%perper
endif


es%total_num_splines1=compute_total_num_splines( &
  num_cells1, &    
  spline_degree1, &
  bc1_min, &
  bc1_max)
es%total_num_splines2=compute_total_num_splines( &
  num_cells2, &    
  spline_degree2, &
  bc2_min, &
  bc2_max)
vec_sz1 = compute_vec_size( &
  num_cells1, &    
  spline_degree1, &
  bc1_min, &
  bc1_max)
vec_sz2 = compute_vec_size( &
  num_cells2, &    
  spline_degree2, &
  bc2_min, &
  bc2_max)
vec_sz = vec_sz1*vec_sz2  
!es%perper = .true. 

SLL_ALLOCATE(es%rho_vec(vec_sz),ierr)
SLL_ALLOCATE(es%masse(vec_sz),ierr)
SLL_ALLOCATE(es%stiff(vec_sz),ierr)

!AB : We must add plus 1 for the dimension of the 
!AB : solution in the case periodic periodic to 
!AB : include the periodicity in the last point.  

solution_size = es%total_num_splines1*es%total_num_splines2

if(es%with_constraint) then
  SLL_ALLOCATE(es%tmp_rho_vec(solution_size+1),ierr)
  SLL_ALLOCATE(es%phi_vec(solution_size+1),ierr)
else
  SLL_ALLOCATE(es%tmp_rho_vec(solution_size),ierr)
  SLL_ALLOCATE(es%phi_vec(solution_size),ierr) 
endif
es%rho_vec = 0.0_f64
es%masse   = 0.0_f64
es%stiff   = 0.0_f64

es%csr_mat => sll_f_new_csr_matrix( solution_size,        &
&                             solution_size,              &
&                             num_cells1*num_cells2,      &
&                             es%local_to_global_indices, &
&                             es%total_num_splines_loc,   &
&                             es%local_to_global_indices, &
&                             es%total_num_splines_loc,   &
&                             sll_p_umfpack)

if(es%precompute_rhs)then
es%csr_mat_source => sll_f_new_csr_matrix( &
  vec_sz, &
  (num_cells1+es%rho_degree1)*(num_cells2+es%rho_degree2), &
  num_cells1*num_cells2, &
  es%local_to_global_indices_source_bis, &
  (es%spline_degree1+1)*(es%spline_degree2+1), &
  es%local_to_global_indices_source, &
  (es%rho_degree1+1)*(es%rho_degree2+1))
endif
 
print *,'#sll_f_new_csr_matrix done'
flush( output_unit )

! allocation of the table containing 
! all values of splines and its
! derivatives in each gauss points
SLL_ALLOCATE(es%v_splines1(2,spline_degree1+1,spline_degree1+2,num_cells1),ierr)
SLL_ALLOCATE(es%v_splines2(2,spline_degree2+1,spline_degree2+2,num_cells2),ierr)
SLL_ALLOCATE(es%v_rho_splines1(es%rho_degree1+1,spline_degree1+2,num_cells1),ierr)
SLL_ALLOCATE(es%v_rho_splines2(es%rho_degree2+1,spline_degree2+2,num_cells2),ierr)



es%v_splines1 = 0.0_f64
es%v_splines2 = 0.0_f64


num_pts_g1 = size(es%gauss_pts1,2)
num_pts_g2 = size(es%gauss_pts2,2)




call compute_global_knots( &
  bc1_min, &
  bc1_max, &
  num_cells1, &
  eta1_min, &
  eta1_max, &
  spline_degree1, &
  es%knots1)

call compute_global_knots( &
  bc2_min, &
  bc2_max, &
  num_cells2, &
  eta2_min, &
  eta2_max, &
  spline_degree2, &
  es%knots2)

call compute_global_knots( &
  bc1_min, &
  bc1_max, &
  num_cells1, &
  eta1_min, &
  eta1_max, &
  es%rho_degree1, &
  es%knots1_rho)

call compute_global_knots( &
  bc2_min, &
  bc2_max, &
  num_cells2, &
  eta2_min, &
  eta2_max, &
  es%rho_degree2, &
  es%knots2_rho)



print *,'#compute_global_knots done'
flush( output_unit )


call compute_non_zero_splines_and_deriv_at_cell_points( &
  num_cells1, &
  spline_degree1, &
  es%knots1, &
  es%gauss_pts1(1,:)/es%delta_eta1, & 
  num_pts_g1, &
  es%v_splines1 )

call compute_non_zero_splines_and_deriv_at_cell_points( &
  num_cells2, &
  spline_degree2, &
  es%knots2, &
  es%gauss_pts2(1,:)/es%delta_eta2, & 
  num_pts_g2, &
  es%v_splines2 )
print *,'#compute_non_zero_splines_and_deriv_at_cell_points done'

call compute_non_zero_splines_at_cell_points( &
  num_cells1, &
  es%rho_degree1, &
  es%knots1_rho, &
  es%gauss_pts1(1,:)/es%delta_eta1, & 
  num_pts_g1, &
  es%v_rho_splines1 )


call compute_non_zero_splines_at_cell_points( &
  num_cells2, &
  es%rho_degree2, &
  es%knots2_rho, &
  es%gauss_pts2(1,:)/es%delta_eta2, & 
  num_pts_g2, &
  es%v_rho_splines2 )





flush( output_unit )

end subroutine sll_s_initialize_general_elliptic_solver_prototype



subroutine sll_s_factorize_mat_es_prototype( es,            &
                             a11_field_mat, &
                             a12_field_mat, &
                             a21_field_mat, &
                             a22_field_mat, &
                             b1_field_vect, &
                             b2_field_vect, &
                             c_field)

type(sll_t_general_coordinate_elliptic_solver),intent(inout) :: es

class(sll_c_scalar_field_2d_base), pointer :: a11_field_mat
class(sll_c_scalar_field_2d_base), pointer :: a12_field_mat
class(sll_c_scalar_field_2d_base), pointer :: a21_field_mat
class(sll_c_scalar_field_2d_base), pointer :: a22_field_mat
class(sll_c_scalar_field_2d_base), pointer :: b1_field_vect
class(sll_c_scalar_field_2d_base), pointer :: b2_field_vect
class(sll_c_scalar_field_2d_base), pointer :: c_field


sll_real64, dimension(:,:),   allocatable :: elt_mat_loc
sll_real64, dimension(:),     allocatable :: mass
sll_real64, dimension(:),     allocatable :: stif
sll_real64, dimension(:,:), allocatable     :: source

sll_int32 :: ierr
sll_int32 :: i
sll_int32 :: j
sll_int32 :: icell
sll_int32 :: bc1_min
sll_int32 :: bc1_max
sll_int32 :: bc2_min
sll_int32 :: bc2_max

sll_real64 :: delta1
sll_real64 :: delta2
sll_real64 :: eta1_min
sll_real64 :: eta2_min
sll_real64 :: eta1
sll_real64 :: eta2
sll_int32  :: num_pts_g1
sll_int32  :: num_pts_g2
sll_int32  :: ii,kk,mm
sll_int32  :: jj,ll,nn
sll_int32  :: ig, jg
sll_real64 :: xg, wxg
sll_real64 :: yg, wyg
sll_int32  :: index1
!sll_int32  :: index2
sll_real64 :: val_c
sll_real64 :: val_a11
sll_real64 :: val_a12
sll_real64 :: val_a21
sll_real64 :: val_a22
sll_real64 :: val_b1=0.0_f64
sll_real64 :: val_b1_der1=0.0_f64
sll_real64 :: val_b1_der2=0.0_f64
sll_real64 :: val_b2=0.0_f64
sll_real64 :: val_b2_der1=0.0_f64
sll_real64 :: val_b2_der2=0.0_f64
sll_real64 :: jac_mat(2,2)
sll_real64 :: val_jac
sll_real64 :: B11
sll_real64 :: B12
sll_real64 :: B21
sll_real64 :: B22
sll_real64 :: MC
sll_real64 :: C1
sll_real64 :: C2    
sll_int32  :: index3!, index4
!sll_int32  :: index_coef1!,index_coef2
sll_int32  :: b, bprime,x!,y
sll_int32  :: a, aprime
!sll_real64 :: elt_mat_global
sll_int32  :: nspl!, nbsp,nbsp1
sll_int32  :: spl_deg_1, spl_deg_2, nc_1, nc_2
!sll_int32  :: ideg2,ideg1
!sll_int32  :: jdeg2,jdeg1
sll_real64 :: v1, v2, v3, v4, d1, d2, d3, d4, r1, r2
sll_real64 :: wxy
sll_real64 :: wxy_by_val_jac 
sll_real64 :: wxy_val_jac 
sll_real64 :: r1r2, v1v2, d1v2, v1d2
sll_real64 :: d3v4, v3v4 , v3d4
!sll_real64 :: intjac
sll_int32, dimension(:), allocatable :: index_eta1
sll_int32, dimension(:), allocatable :: index_eta2
sll_int32 :: vec_sz1
sll_real64 :: val

!sll_real64, allocatable :: dense_matrix(:,:)

!!$ sll_int32 :: tid=0
!!$ sll_int32 :: nthreads=1

bc1_min    = es%bc1_min
bc1_max    = es%bc1_max
bc2_min    = es%bc2_min
bc2_max    = es%bc2_max
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

!intjac = 0.0_f64


SLL_ALLOCATE(index_eta1(nc_1+spl_deg_1),ierr)
SLL_ALLOCATE(index_eta2(nc_2+spl_deg_2),ierr)

call compute_index( &
  nc_1, &
  spl_deg_1, &
  bc1_min, &
  bc1_max, &
  index_eta1)
call compute_index( &
  nc_2, &
  spl_deg_2, &
  bc2_min, &
  bc2_max, &
  index_eta2)

vec_sz1 = compute_vec_size( &
  nc_1, &
  spl_deg_1, &
  bc1_min, &
  bc1_max)

print *,'begin sll_s_factorize_mat_es_prototype'
flush( output_unit )
SLL_CLEAR_ALLOCATE(elt_mat_loc(1:nspl,1:nspl),ierr)
SLL_CLEAR_ALLOCATE(mass(1:nspl),ierr)
SLL_CLEAR_ALLOCATE(stif(1:nspl),ierr)
SLL_CLEAR_ALLOCATE(source(1:nspl,(es%rho_degree1+1)*(es%rho_degree2+1)),ierr)

do j = 1, nc_2
do i = 1, nc_1
        
  icell = i + (j-1) * nc_1
  eta1  = eta1_min + (i-1)*delta1
  eta2  = eta2_min + (j-1)*delta2
    
  mass  = 0.0_f64
  stif  = 0.0_f64
  elt_mat_loc = 0._f64
  source = 0._f64


  do jg=1,num_pts_g2
  
    yg  = eta2+es%gauss_pts2(1,jg)
    wyg = es%gauss_pts2(2,jg)
  
    do ig=1,num_pts_g1
    
      xg  = eta1+es%gauss_pts1(1,ig)
      wxg = es%gauss_pts1(2,ig)

      wxy = wxg*wyg

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
      val_jac = jac_mat(1,1)*jac_mat(2,2)-jac_mat(1,2)*jac_mat(2,1)

      wxy_by_val_jac = wxy/val_jac
      wxy_val_jac    = wxy*val_jac

      val_c = val_c * wxy_val_jac
        
      !intjac = intjac + wxy_val_jac

      ! The B matrix is  by (J^(-1)) A^t (J^(-1))^t 
      B11 = jac_mat(2,2)*jac_mat(2,2)*val_a11 - &
            jac_mat(2,2)*jac_mat(1,2)*(val_a12+val_a21) + &
            jac_mat(1,2)*jac_mat(1,2)*val_a22

      B11 = B11*wxy_by_val_jac
          
      B21 = jac_mat(1,1)*jac_mat(2,2)*val_a12 - &
            jac_mat(1,1)*jac_mat(1,2)*val_a22 - &
            jac_mat(2,1)*jac_mat(2,2)*val_a11 + &
            jac_mat(1,2)*jac_mat(2,1)*val_a21

      B21 = B21*wxy_by_val_jac
        
      B12 = jac_mat(1,1)*jac_mat(2,2)*val_a21 - &
            jac_mat(1,1)*jac_mat(1,2)*val_a22 - &
            jac_mat(2,1)*jac_mat(2,2)*val_a11 + &
            jac_mat(1,2)*jac_mat(2,1)*val_a12

      B12 = B12*wxy_by_val_jac

      B22 = jac_mat(1,1)*jac_mat(1,1)*val_a22 - &
            jac_mat(1,1)*jac_mat(2,1)*(val_a21+val_a12) + &
            jac_mat(2,1)*jac_mat(2,1)*val_a11
          
      B22 = B22*wxy_by_val_jac

      MC =  jac_mat(2,2)*val_b1_der1 &
          - jac_mat(2,1)*val_b1_der2 &
          - jac_mat(1,2)*val_b2_der1 &
          + jac_mat(1,1)*val_b2_der2

      MC = MC*wxy
          
      C1 = jac_mat(2,2)*val_b1-jac_mat(1,2)*val_b2 
      C2 = jac_mat(1,1)*val_b2-jac_mat(2,1)*val_b1

      C1 = C1*wxy
      C2 = C2*wxy
      
      
      !fill local matrix contribution 
      !for the gauss point ig,jg   
      mm = 0
      do jj = 1,spl_deg_2+1
        v2 = es%v_splines2(1,jj,jg,j)
        d2 = es%v_splines2(2,jj,jg,j)
        do ii = 1,spl_deg_1+1
          mm = mm+1
          v1 = es%v_splines1(1,ii,ig,i)
          d1 = es%v_splines1(2,ii,ig,i)
          v1v2 = v1*v2
          d1v2 = d1*v2
          v1d2 = v1*d2
          mass(mm) = mass(mm)+v1v2*wxy_val_jac 
          stif(mm) = stif(mm)+wxy_val_jac*(d1v2+v1d2)
          nn = 0
          do ll = 1,spl_deg_2+1
            v4 = es%v_splines2(1,ll,jg,j)
            d4 = es%v_splines2(2,ll,jg,j)
            do kk = 1,spl_deg_1+1
              v3 = es%v_splines1(1,kk,ig,i)
              d3 = es%v_splines1(2,kk,ig,i)
              v3v4 = v4*v3
              d3v4 = d3*v4
              v3d4 = v3*d4
              nn = nn+1
              val = val_c*v1v2*v3v4                  
              val = val-B11*d1v2*d3v4                 
              val = val-B22*v1d2*v3d4                 
              val = val-B12*d1v2*v3d4
              val = val-B21*v1d2*d3v4                 
              val = val-MC*v1v2*v3v4                 
              val = val-C1*v1v2*d3v4                 
              val = val-C2*v1v2*v3d4
              elt_mat_loc(nn,mm) = elt_mat_loc(nn,mm)+val                   
            end do
          end do
        end do
      end do

      !fill source local matrix contribution 
      !for the gauss point ig,jg   
      if(es%precompute_rhs)then
      mm = 0
      do jj = 1,es%rho_degree2+1
        r2 = es%v_rho_splines2(jj,jg,j)
        do ii = 1,es%rho_degree1+1
          mm = mm+1
          r1 = es%v_rho_splines1(ii,ig,i)
          r1r2 = r1*r2*wxy_val_jac
          nn = 0
          do ll = 1,spl_deg_2+1
            v4 = es%v_splines2(1,ll,jg,j)
            do kk = 1,spl_deg_1+1
              v3 = es%v_splines1(1,kk,ig,i)
              v3v4 = v4*v3
              nn = nn+1
              source(nn,mm) = source(nn,mm) + r1r2*v3v4
            end do
          end do
        end do
      end do
      endif

    end do
  end do
  !local matrix and source matrix are computed
  !now, add contribution to csr matrix

  do jj = 0, spl_deg_2
    index3 = index_eta2(j + jj)
    do ii = 0,spl_deg_1
      index1 = index_eta1(i + ii)

      x = index1 + (index3-1)*vec_sz1
      b = ii+1+jj*(spl_deg_1+1)
      es%masse(x) = es%masse(x) + mass(b)
      es%stiff(x) = es%stiff(x) + stif(b)

      a = es%local_to_global_indices(b, icell)
      do ll = 0,spl_deg_2
        do kk = 0,spl_deg_1
          bprime = kk+1+ll*(spl_deg_1+1)
          aprime = es%local_to_global_indices(bprime,icell)
          if ( a>0 .and. aprime>0 ) then
            call sll_s_add_to_csr_matrix( &
              es%csr_mat, &
              elt_mat_loc(bprime,b), &
              a, &
              aprime)   
          end if
        end do
      end do
    end do
  end do
  
  !add contribution to csr source matrix
  if(es%precompute_rhs)then
  do jj = 0, spl_deg_2
    do ii = 0,spl_deg_1
      b = ii+1+jj*(spl_deg_1+1)
      a = es%local_to_global_indices_source_bis(b, icell)
      do ll = 0,es%rho_degree2
        do kk = 0,es%rho_degree1
          bprime = kk+1+ll*(es%rho_degree1+1)
          aprime = es%local_to_global_indices_source(bprime,icell)
          if ( a>0 .and. aprime>0 ) then
            call sll_s_add_to_csr_matrix( &
              es%csr_mat_source, &
              source(b,bprime), &
              a, &
              aprime)              
            !print *,'val=',a,aprime,source(b,bprime)     
          end if
        end do
      end do
    end do
  end do
  endif
  
  
  
end do
end do

print *,'end matrix construction sll_s_factorize_mat_es_prototype'
flush( output_unit )

!es%intjac = intjac
es%intjac = sum(es%masse)
print *,'#begin of sll_s_factorize_csr_matrix'

if (es%with_constraint) then

 es%csr_mat_with_constraint => sll_f_new_csr_matrix_with_constraint(es%csr_mat)  

 call sll_s_csr_add_one_constraint( es%csr_mat%row_ptr,                 &  
                                    es%csr_mat%col_ind,                 &
                                    es%csr_mat%val,                     &
                                    es%csr_mat%num_rows,                &
                                    es%csr_mat%num_nz,                  &
                                    es%masse,                           &
                                    es%csr_mat_with_constraint%row_ptr, &
                                    es%csr_mat_with_constraint%col_ind, &
                                    es%csr_mat_with_constraint%val)  

  call sll_s_factorize_csr_matrix(es%csr_mat_with_constraint)

else   

  call sll_s_factorize_csr_matrix(es%csr_mat)

end if 

print *,'#end of sll_s_factorize_csr_matrix'
flush( output_unit )


!SLL_DEALLOCATE_ARRAY(M_c,ierr)
!SLL_DEALLOCATE_ARRAY(K_11,ierr)
!SLL_DEALLOCATE_ARRAY(K_12,ierr)
!SLL_DEALLOCATE_ARRAY(K_21,ierr)
!SLL_DEALLOCATE_ARRAY(K_22,ierr)
!SLL_DEALLOCATE_ARRAY(M_bv,ierr)
!SLL_DEALLOCATE_ARRAY(S_b1,ierr)
!SLL_DEALLOCATE_ARRAY(S_b2,ierr)
SLL_DEALLOCATE_ARRAY(stif,ierr) 
SLL_DEALLOCATE_ARRAY(mass,ierr) 

print *,'#dealloc ok sll_s_factorize_csr_matrix'
flush( output_unit )
   
end subroutine sll_s_factorize_mat_es_prototype

subroutine sll_s_set_rho_coefficients_coordinates_elliptic_eq_prototype( &
  es, &
  rho_values, &
  rho_coeff)
class(sll_t_general_coordinate_elliptic_solver), intent(inout) :: es
sll_real64, dimension(:,:), intent(in), optional :: rho_values
sll_real64, dimension(:,:), intent(in), optional :: rho_coeff

sll_int32 :: num_spl1
sll_int32 :: num_spl2
sll_int32 :: i
sll_int32 :: j

num_spl1 = es%num_cells1+es%spline_degree1
num_spl2 = es%num_cells2+es%spline_degree2

  
if(.not. es%precompute_rhs)then
  SLL_ERROR('sll_s_set_rho_coefficients_coordinates_elliptic_eq_prototype','should not be used as es%precompute_rhs is set to false')             
endif
 
if(present(rho_values).and.present(rho_coeff))then
  SLL_ERROR('sll_s_set_rho_coefficients_coordinates_elliptic_eq_prototype','provide either rho_values or rho_coeff, not both')         
endif
  
if(present(rho_values))then
  if(es%use_cubic_splines)then
    call sll_s_compute_cubic_spline_2d( rho_values, es%cubic_spline )
    call sll_s_get_coeff_cubic_spline_2d(es%cubic_spline,es%rho_coeff_1d)
  else
  SLL_ERROR('sll_s_set_rho_coefficients_coordinates_elliptic_eq_prototype','spline coeff computed only with cubic splines for the moment if precompute_rhs is true')                 
  endif
else
  if(present(rho_coeff))then
    if(size(rho_coeff,1)<num_spl1)then
      print *,'#size(rho_coeff,1)=',size(rho_coeff,1)
      flush( output_unit )
      SLL_ERROR('sll_s_set_rho_coefficients_coordinates_elliptic_eq_prototype','bad size(rho_coeff,1)')
    endif
    if(size(rho_coeff,2)<num_spl2)then
      print *,'#size(rho_coeff,2)=',size(rho_coeff,2)
      flush( output_unit )
      SLL_ERROR('sll_s_set_rho_coefficients_coordinates_elliptic_eq_prototype','bad size(rho_coeff,2)')
    endif
    do j = 1,num_spl2
      do i = 1,num_spl1
        es%rho_coeff_1d(i+num_spl1*(j-1)) = rho_coeff(i,j)
      enddo
    enddo
  else
    SLL_ERROR('sll_s_set_rho_coefficients_coordinates_elliptic_eq_prototype','provides rho_values or rho_coeff')     
  endif
endif

end subroutine sll_s_set_rho_coefficients_coordinates_elliptic_eq_prototype

subroutine sll_s_solve_general_coordinates_elliptic_eq_prototype( &
  es, &
  phi, &
  rho_field)

class(sll_t_general_coordinate_elliptic_solver), intent(inout)      :: es
sll_real64, dimension(:,:), intent(out) :: phi
class(sll_c_scalar_field_2d_base),           intent(in),target, optional  :: rho_field
!class(sll_t_scalar_field_2d_discrete),       intent(inout)      :: phi
!
!sll_real64, dimension(:,:),      pointer :: coeff_rho
sll_real64, dimension(:), allocatable :: m_rho_loc

sll_int32  :: i
sll_int32  :: j
sll_int32  :: k
sll_int32  :: num_pts_g1
sll_int32  :: num_pts_g2
sll_int32  :: x, n, b
sll_int32  :: ii, jj, kk, ll, mm, nn
sll_int32  :: bc1_min, bc1_max, bc2_min, bc2_max
sll_int32  :: index1, index3!, nbsp

sll_real64 :: wgpt1, wgpt2, gpt1, gpt2, eta1, eta2
sll_real64 :: val_f, val_j, valfj, jac_mat(2,2)
sll_real64 :: spline1, spline2

sll_real64 :: int_rho
sll_real64 :: int_jac
sll_int32  :: ierr
sll_int32  :: nc_1, nc_2

sll_int32, dimension(:), allocatable :: index_eta1
sll_int32, dimension(:), allocatable :: index_eta2
sll_int32 :: vec_sz1
sll_int32 :: vec_sz2
sll_int32, dimension(:), allocatable :: index_sol1
sll_int32, dimension(:), allocatable :: index_sol2
sll_real64, dimension(:,:), allocatable :: coeff_phi
sll_int32, dimension(:), allocatable :: index_coeff1
sll_int32, dimension(:), allocatable :: index_coeff2



!!$ sll_int32  :: tid = 0
!!$ sll_int32  :: nthreads = 1

print *,'#solve begin'
flush( output_unit )

  
num_pts_g1 = size(es%gauss_pts1,2)
num_pts_g2 = size(es%gauss_pts2,2)


nc_1            = es%num_cells1
nc_2            = es%num_cells2
es%rho_vec      = 0.0_f64
 
bc1_min = es%bc1_min
bc1_max = es%bc1_max
bc2_min = es%bc2_min
bc2_max = es%bc2_max




SLL_ALLOCATE(index_eta1(nc_1+es%spline_degree1),ierr)
SLL_ALLOCATE(index_eta2(nc_2+es%spline_degree2),ierr)
SLL_ALLOCATE(index_sol1(es%total_num_splines1),ierr)
SLL_ALLOCATE(index_sol2(es%total_num_splines2),ierr)
SLL_ALLOCATE(index_coeff1(nc_1+es%spline_degree1),ierr)
SLL_ALLOCATE(index_coeff2(nc_2+es%spline_degree2),ierr)

call compute_index( &
  nc_1, &
  es%spline_degree1, &
  bc1_min, &
  bc1_max, &
  index_eta1)
call compute_index( &
  nc_2, &
  es%spline_degree2, &
  bc2_min, &
  bc2_max, &
  index_eta2)
vec_sz1 = compute_vec_size( &
  nc_1, &
  es%spline_degree1, &
  bc1_min, &
  bc1_max)
vec_sz2 = compute_vec_size( &
  nc_2, &
  es%spline_degree2, &
  bc2_min, &
  bc2_max)
call compute_index_sol( &
  es%total_num_splines1, &
  bc1_min, &
  bc1_max, &
  index_sol1)
call compute_index_sol( &
  es%total_num_splines2, &
  bc2_min, &
  bc2_max, &
  index_sol2)
call compute_index_coeff( &
  nc_1, &
  es%spline_degree1, &
  bc1_min, &
  bc1_max, &
  index_coeff1)
call compute_index_coeff( &
  nc_2, &
  es%spline_degree2, &
  bc2_min, &
  bc2_max, &
  index_coeff2)

int_rho = 0.0_f64
int_jac = 0.0_f64

SLL_CLEAR_ALLOCATE(M_rho_loc(1:es%total_num_splines_loc),ierr)

if(es%precompute_rhs)then
  call sll_s_mult_csr_matrix_vector( &
    es%csr_mat_source,&
    es%rho_coeff_1d, &
    es%rho_vec)    
else
  if(present(rho_field))then
do j=1, nc_2
  do i=1, nc_1
    M_rho_loc = 0.0_f64
    eta1  = es%eta1_min + real(i-1,f64)*es%delta_eta1
    eta2  = es%eta2_min + real(j-1,f64)*es%delta_eta2
    do jj=1,num_pts_g2
      gpt2  = eta2 + es%gauss_pts2(1,jj)
      wgpt2 = es%gauss_pts2(2,jj)
      do ii=1,num_pts_g1
        gpt1  = eta1 + es%gauss_pts1(1,ii)
        wgpt1 = es%gauss_pts1(2,ii)
    
        val_f   = rho_field%value_at_point(gpt1,gpt2)
        
        jac_mat = rho_field%get_jacobian_matrix(gpt1,gpt2)
        val_j   = jac_mat(1,1)*jac_mat(2,2)-jac_mat(1,2)*jac_mat(2,1)
        val_j   = val_j*wgpt1*wgpt2
        valfj   = val_f*val_j
        int_rho = int_rho + valfj 
        int_jac = int_jac + val_j
    
        do ll = 1,es%spline_degree2+1
          spline2 = es%v_splines2(1,ll,jj,j)*valfj
          do kk = 1,es%spline_degree1+1
            spline1 = es%v_splines1(1,kk,ii,i)
            n = kk + (ll-1)*(es%spline_degree1+1)
            M_rho_loc(n)= M_rho_loc(n) + spline1*spline2
          end do
        end do
      end do
    end do

    do mm = 0,es%spline_degree2
      index3 = index_eta2(j + mm)
    
      do nn = 0,es%spline_degree1
        index1 = index_eta1(i + nn)      
        x             =  index1 + (index3-1)*vec_sz1
        b             =  nn + 1 + mm * (es%spline_degree1+1)
        es%rho_vec(x) =  es%rho_vec(x)  + M_rho_loc(b)
          
      end do
    end do
  end do
end do
  else
    SLL_ERROR('sll_s_solve_general_coordinates_elliptic_eq_prototype&
    &','rho_field should be given')
  endif  
endif

!if (es%perper) es%rho_vec = es%rho_vec - int_rho/int_jac
print *,'sum(es%rho_vec)=',int_rho,sum(es%rho_vec)
if (es%zero_mean) es%rho_vec = es%rho_vec - sum(es%rho_vec)/es%intjac*es%masse
print *,'after sum(es%rho_vec)=',sum(es%rho_vec)

bc1_min = es%bc1_min
bc1_max = es%bc1_max
bc2_min = es%bc2_min
bc2_max = es%bc2_max
 
es%tmp_rho_vec(:) = 0.0_f64  !PN: Is it useful ?
es%phi_vec(:)     = 0.0_f64  !PN: Is it useful ?


  k = 0
  do j = 1, es%total_num_splines2
    jj = index_sol2(j)
    do i = 1, es%total_num_splines1
      k = k+1
      ii = index_sol1(i)
      es%tmp_rho_vec(k) = es%rho_vec(ii+vec_sz1*(jj-1))
    end do
  end do

  

if(es%with_constraint)then
  call sll_s_solve_csr_matrix(es%csr_mat_with_constraint, &
                              es%tmp_rho_vec,             &
                              es%phi_vec)
else

  call sll_s_solve_csr_matrix(es%csr_mat, es%tmp_rho_vec, es%phi_vec)

endif


!SLL_ALLOCATE(coeff_phi(nc_1+es%total_num_splines1,nc_2+es%total_num_splines2),ierr)
SLL_ALLOCATE(coeff_phi(nc_1+es%spline_degree1,nc_2+es%spline_degree2),ierr)
coeff_phi(:,:) = 0._f64

  k = 0
  do j = 1, nc_2+es%spline_degree2
    jj = index_coeff2(j)
    if(jj/=0)then
      do i = 1, nc_1+es%spline_degree1
        k = k+1
        ii = index_coeff1(i)
        if(ii/=0)then          
          coeff_phi(i,j) = &
            es%phi_vec(ii+es%total_num_splines1*(jj-1))
          !coeff_phi(i,j) = es%phi_vec(ii+vec_sz1*(jj-1))
        endif  
      end do
    endif  
  end do

  call compute_values_at_knots( &
    es%knots1, &
    es%knots2, &
    nc_1, &
    es%spline_degree1, &
    nc_2, &
    es%spline_degree2, &
    coeff_phi, &
    phi)


  

print *,'#solve done'
flush( output_unit )

end subroutine sll_s_solve_general_coordinates_elliptic_eq_prototype

function compute_total_num_splines( &
  num_cells, &
  degree, &
  bc_min, &
  bc_max) &
  result(res)
  
  sll_int32, intent(in) :: num_cells
  sll_int32, intent(in) :: degree
  sll_int32, intent(in) :: bc_min
  sll_int32, intent(in) :: bc_max
  sll_int32 :: res
  
  if(bc_min==sll_p_periodic)then
    res = num_cells
  else
    res = num_cells+degree
  endif
  
  
  if(bc_min==sll_p_dirichlet)then
    res = res-1
  endif
  if(bc_max==sll_p_dirichlet)then
    res = res-1
  endif
  
  
end function compute_total_num_splines


function compute_vec_size( &
  num_cells, &
  degree, &
  bc_min, &
  bc_max) &
  result(res)
  
  sll_int32, intent(in) :: num_cells
  sll_int32, intent(in) :: degree
  sll_int32, intent(in) :: bc_min
  sll_int32, intent(in) :: bc_max
  sll_int32 :: res
  
  if(bc_min==sll_p_periodic)then
    res = num_cells
  else
    res = num_cells+degree
  endif
  
  if(bc_max==sll_p_dirichlet)then
  end if
  
  
end function compute_vec_size

subroutine compute_index( &
  num_cells, &
  degree, &
  bc_min, &
  bc_max, &
  index)
  sll_int32, intent(in) :: num_cells
  sll_int32, intent(in) :: degree
  sll_int32, intent(in) :: bc_min
  sll_int32, intent(in) :: bc_max
  sll_int32, dimension(:), intent(out) :: index
  
  sll_int32 :: i
  
  do i=1,num_cells+degree
    index(i) = i
  enddo
  
  if(bc_min==sll_p_periodic)then
    do i=1,degree
      index(num_cells+i) = index(i)
    enddo
  endif
  if(bc_max==sll_p_dirichlet)then
  end if

end subroutine compute_index


subroutine compute_index_sol( &
  total_num_splines, &
  bc_min, &
  bc_max, &
  index)
  sll_int32, intent(in) :: total_num_splines
  sll_int32, intent(in) :: bc_min
  sll_int32, intent(in) :: bc_max
  sll_int32, dimension(:), intent(out) :: index
  
  sll_int32 :: i
  
  do i=1,total_num_splines
    index(i) = i
  enddo
  
  if(bc_min==sll_p_dirichlet)then
    do i=1,total_num_splines
      index(i) = i+1
    enddo
  endif
  if(bc_max==sll_p_dirichlet)then
  end if

end subroutine compute_index_sol


subroutine compute_index_coeff( &
  num_cells, &
  degree, &
  bc_min, &
  bc_max, &
  index)
  sll_int32, intent(in) :: num_cells
  sll_int32, intent(in) :: degree
  sll_int32, intent(in) :: bc_min
  sll_int32, intent(in) :: bc_max
  sll_int32, dimension(:), intent(out) :: index
  
  sll_int32 :: i
  
  do i=1,num_cells+degree
    index(i) = i
  enddo
  
  
  if(bc_min==sll_p_dirichlet)then
    index(1) = 0
    do i=2,num_cells+degree
      index(i) = i-1
    enddo
  endif
  if(bc_max==sll_p_dirichlet)then
    index(num_cells+degree) = 0
  endif
  
  if(bc_min==sll_p_periodic)then
    do i=1,degree
      index(i+num_cells) = index(i)
    enddo
  endif

end subroutine compute_index_coeff


subroutine compute_values_at_knots( &
  knots1, &
  knots2, &
  num_cells1, &
  degree1, &
  num_cells2, &
  degree2, &
  input, &
  output)
  sll_real64, dimension(:), intent(in) :: knots1
  sll_real64, dimension(:), intent(in) :: knots2
  sll_int32, intent(in) :: num_cells1
  sll_int32, intent(in) :: degree1
  sll_int32, intent(in) :: num_cells2
  sll_int32, intent(in) :: degree2
  sll_real64, dimension(:,:), intent(in) :: input 
  sll_real64, dimension(:,:), intent(out) :: output

  sll_real64, dimension(:,:), allocatable :: tmp   
  sll_real64, dimension(:), allocatable :: biat1  
  sll_real64, dimension(:), allocatable :: biat2  
  sll_int32 :: ierr
  sll_int32 :: i
  sll_int32 :: j
  !sll_real64 :: err
  
  
  if(num_cells1<1)then
    print *,'#num_cells1=',num_cells1
    SLL_ERROR('compute_values_at_knots','bad value of num_cells1')
  endif
  if(num_cells2<1)then
    print *,'#num_cells2=',num_cells2
    SLL_ERROR('compute_values_at_knots','bad value of num_cells2')
  endif
  if(degree1<0)then
    print *,'#degree1=',degree1
    SLL_ERROR('compute_values_at_knots','bad value of degree1')
  endif
  if(degree2<0)then
    print *,'#degree2=',degree2
    SLL_ERROR('compute_values_at_knots','bad value of degree2')
  endif
  
  if(size(knots1)<num_cells1+2*degree1+1)then
    print *,'#size(knots1)=',size(knots1)
    SLL_ERROR('compute_values_at_knots&
    &','bad value of size(knots1)')
  endif
  if(size(knots2)<num_cells2+2*degree2+1)then
    print *,'#size(knots2)=',size(knots2)
    SLL_ERROR('compute_values_at_knots&
    &','bad value of size(knots2)')
  endif

  if(size(input,1)<num_cells1+degree1)then
    print *,'#size(input,1)=',size(input,1)
    SLL_ERROR('compute_values_at_knots&
    &','bad value of size(input,1)')
  endif
  if(size(input,2)<num_cells2+degree2)then
    print *,'#size(input,2)=',size(input,2)
    SLL_ERROR('compute_values_at_knots&
    &','bad value of size(input,2)')
  endif

  if(size(output,1)<num_cells1+1)then
    print *,'#size(output,1)=',size(output,1)
    SLL_ERROR('compute_values_at_knots&
    &','bad value of size(output,1)')
  endif

  if(size(output,2)<num_cells2+1)then
    print *,'#size(output,2)=',size(output,2)
    SLL_ERROR('compute_values_at_knots&
    &','bad value of size(output,2)')
  endif

  SLL_ALLOCATE(tmp(num_cells1+1,num_cells2+degree2),ierr)
  SLL_ALLOCATE(biat1(degree1+1),ierr)
  SLL_ALLOCATE(biat2(degree2+1),ierr)
  
  do j=1,num_cells2+degree2
    do i=1,num_cells1
      call compute_b_spline_at_x( &
        knots1, &
        degree1+i, &
        knots1(degree1+i), &
        degree1, &
        biat1)
      tmp(i,j) = dot_product( &
        biat1(1:degree1+1),&
        input(i:i+degree1,j))  
    enddo
    call compute_b_spline_at_x( &
      knots1, &
      degree1+num_cells1, &
      knots1(degree1+num_cells1+1), &
      degree1, &
      biat1)
    tmp(num_cells1+1,j) = dot_product( &
      biat1(1:degree1+1),&
      input(num_cells1:num_cells1+degree1,j))  
  enddo


  do i=1,num_cells1+1
    do j=1,num_cells2
      call compute_b_spline_at_x( &
        knots2, &
        degree2+j, &
        knots2(degree2+j), &
        degree2, &
        biat2)
      output(i,j) = dot_product( &
        biat2(1:degree2+1),&
        tmp(i,j:j+degree2))  
    enddo
    call compute_b_spline_at_x( &
      knots2, &
      degree2+num_cells2, &
      knots2(degree2+num_cells2+1), &
      degree2, &
      biat2)
    output(i,num_cells2+1) = dot_product( &
      biat2(1:degree2+1),&
      tmp(i,num_cells2:num_cells2+degree2))  
  enddo
  
!  err = 0._f64
!  !check periodicity
!  do i=1,num_cells1+1
!    err = max(err,abs(output(i,1)-output(i,num_cells2+1)))
!  enddo
!  do i=1,num_cells2+1
!    err = max(err,abs(output(1,i)-output(num_cells1+1,i)))
!  enddo
  
  !print *,'#err periodicity=',errg
  
end subroutine compute_values_at_knots



!subroutine compute_local_to_global_indices_source( &
!  num_cells1, &
!  num_cells2, &
!  degree1, &
!  degree2, &
!  degree1_bis, &
!  degree2_bis, &
!  local_to_global_indices_source, &
!  local_to_global_indices_source_bis, &
!  bc_index1, & !index1
!  bc_index2, & !index1_bis 
!  bc_index3, & !index2
!  bc_index4, & !index2_bis
!  
!
!
!do j = 1, es%num_cells2
!do i = 1, es%num_cells1
!  icell = i + (j-1)*es%num_cells1
!  b = 0
!  do jj = 0, es%spline_degree2
!    index3 = j + jj
!    if (bc2_min==sll_p_periodic .and. bc2_max==sll_p_periodic) then 
!      if ( index3 > es%total_num_splines2) then
!        index3 = index3 - es%total_num_splines2
!      end if
!    end if
!    do ii = 0,es%spline_degree1
!      index1 = i + ii
!      if (bc1_min==sll_p_periodic .and. bc1_max==sll_p_periodic) then 
!        if ( index1 > es%total_num_splines1) then
!          index1 = index1 - es%total_num_splines1
!        end if
!        nbsp = es%total_num_splines1
!      else !if (bc1_min==sll_p_dirichlet .and. bc1_max==sll_p_dirichlet) then
!        nbsp = es%num_cells1 + es%spline_degree1
!      end if
!      x = index1 + (index3-1)*nbsp
!      b = b+1
!      index_coef1 = tab_index_coeff1(i) - es%spline_degree1 + ii
!      index_coef2 = tab_index_coeff2(j) - es%spline_degree2 + jj
!      es%local_to_global_indices_source(b,icell)= index_coef1 + (index_coef2-1)*(es%num_cells1+1)
!      bprime = 0
!      do ll = 0,es%spline_degree2
!        index4 = j + ll
!        if ( (bc2_min==sll_p_periodic).and.(bc2_max== sll_p_periodic))then
!          if ( index4 > es%total_num_splines2) then
!            index4 = index4 - es%total_num_splines2
!          end if
!        end if
!        do kk = 0,es%spline_degree1
!          index2 = i + kk
!          if(bc1_min==sll_p_periodic .and. bc1_max==sll_p_periodic)then
!            if ( index2 > es%total_num_splines1) then
!              index2 = index2 - es%total_num_splines1
!            end if
!            nbsp1 = es%total_num_splines1
!          else !if (bc1_min==sll_p_dirichlet .and. bc1_max==sll_p_dirichlet) then
!            nbsp1 = es%num_cells1 + es%spline_degree1
!          end if
!          y      = index2 + (index4-1)*nbsp1
!          bprime = bprime+1 
!          es%local_to_global_indices_source_bis(bprime,icell)= y
!        end do
!      end do
!    end do
!  end do
!end do
!end do



end module sll_m_general_coordinate_elliptic_solver
