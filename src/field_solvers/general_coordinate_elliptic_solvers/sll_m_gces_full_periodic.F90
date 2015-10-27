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
!> où les \f$ P_i \f$ is a polynomial function \f$ (m-n) \f$ points.
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
module sll_m_gces_full_periodic
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
#include "sll_errors.h"

use sll_m_boundary_condition_descriptors
use sll_m_scalar_field_2d_base, only: sll_scalar_field_2d_base
use sll_m_scalar_field_2d, only: sll_scalar_field_2d_analytic,  &
                                      sll_scalar_field_2d_discrete
use sll_m_interpolators_2d_base, only: sll_interpolator_2d_base
use sll_m_arbitrary_degree_spline_interpolator_2d, only:        &
  sll_arbitrary_degree_spline_interpolator_2d
use sll_m_arbitrary_degree_spline_interpolator_1d, only:        &
  interv, deboor_type, bsplvd
use sll_m_gauss_legendre_integration
use sll_m_gauss_lobatto_integration
use sll_m_sparse_matrix, only : sll_csr_matrix,                 &
                                     new_csr_matrix,                 &
                                     new_csr_matrix_with_constraint, &
                                     sll_factorize_csr_matrix,       &
                                     sll_add_to_csr_matrix,          &
                                     sll_mult_csr_matrix_vector,     &
                                     sll_solve_csr_matrix,           &
                                     sll_delete
#ifdef _OPENMP
use omp_lib
#endif

implicit none

private

!> @brief
!> General coordinate elliptic solver derived type
!> @details
!> The indexing of the
!> splines in array global_indices depends on the boundary conditions.
!> local_indices includes the changes resulting from the boundary conditions.
!> local_to_lobal_indices(i,j) = global_indices(local_indices(i,j))
type, public :: sll_gces_full_periodic

  private
  sll_int32,  public :: num_cells1
  sll_int32,  public :: num_cells2
  sll_real64, public :: delta_eta1
  sll_real64, public :: delta_eta2
  sll_real64, public :: eta1_min
  sll_real64, public :: eta2_min   

  sll_real64, dimension(:),   pointer :: knots1
  sll_real64, dimension(:),   pointer :: knots2
  sll_real64, dimension(:),   pointer :: knots1_rho
  sll_real64, dimension(:),   pointer :: knots2_rho
  sll_real64, dimension(:,:), pointer :: gauss_pts1
  sll_real64, dimension(:,:), pointer :: gauss_pts2
  sll_int32  :: spline_degree1
  sll_int32  :: spline_degree2
  sll_real64 :: intjac

  sll_int32, dimension(:,:), pointer :: local_to_global_indices
  sll_int32, dimension(:,:), pointer :: local_to_global_indices_source
  sll_int32, dimension(:,:), pointer :: local_to_global_indices_source_bis

  !!! contains the values of all splines in all gauss points
  sll_real64, dimension(:,:,:,:), pointer :: v_splines1
  sll_real64, dimension(:,:,:,:), pointer :: v_splines2

  type(sll_csr_matrix),           pointer :: csr_mat_with_constraint
  type(sll_csr_matrix),           pointer :: csr_mat_source
  sll_real64, dimension(:),       pointer :: rho_vec
  sll_real64, dimension(:),       pointer :: phi_vec
  sll_real64, dimension(:),       pointer :: masse
  sll_real64, dimension(:),       pointer :: stiff
  sll_real64, dimension(:),       pointer :: rho_coeff_1d
  type(deboor_type)                       :: db

end type sll_gces_full_periodic

sll_int32, parameter, public :: ES_GAUSS_LEGENDRE = 0
sll_int32, parameter, public :: ES_GAUSS_LOBATTO = 1
  
interface sll_delete
  module procedure delete_gces_full_periodic
end interface sll_delete

interface sll_create
  module procedure initialize_gces_full_periodic
end interface sll_create

interface sll_solve
  module procedure solve_gces_full_periodic
end interface sll_solve

public sll_delete,                     &
       sll_create,                     &
       sll_solve,                      &
       new_general_elliptic_solver,    &
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
!> @param[in]  eta1_min the minimun in the direction eta1
!> @param[in]  eta1_max the maximun in the direction eta1
!> @param[in]  eta2_min the minimun in the direction eta2
!> @param[in]  eta2_max the maximun in the direction eta2
!> @param[out] the type general_coordinate_elliptic_solver
subroutine initialize_gces_full_periodic( &
       es,                                &
       spline_degree1,                    &
       spline_degree2,                    &
       num_cells1,                        &
       num_cells2,                        &
       quadrature_type1,                  &
       quadrature_type2,                  &
       eta1_min,                          &
       eta1_max,                          &
       eta2_min,                          &
       eta2_max)
    
type(sll_gces_full_periodic), intent(out) :: es

sll_int32,  intent(in) :: spline_degree1
sll_int32,  intent(in) :: spline_degree2
sll_int32,  intent(in) :: num_cells1
sll_int32,  intent(in) :: num_cells2
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
sll_int32 :: ierr
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
sll_int32  :: index1, index2, index3, index4
sll_int32  :: icell, bprime

sll_int32, dimension(:),   allocatable :: global_indices 
sll_int32, dimension(:,:), allocatable :: local_indices
sll_int32, dimension(:),   allocatable :: tab_index_coeff1
sll_int32, dimension(:),   allocatable :: tab_index_coeff2

sll_int32 :: nb_spl_x
sll_int32 :: nb_spl_y
sll_int32 :: bloc
sll_int32 :: cell
sll_int32 :: e
sll_int32 :: b
sll_int32 :: d
sll_int32 :: a
sll_int32 :: l

dim1 = (spline_degree1+1)*(spline_degree2+1)
dim2 = (num_cells1*num_cells2)
SLL_ALLOCATE(es%local_to_global_indices(1:dim1,1:dim2),ierr)
SLL_ALLOCATE(es%local_to_global_indices_source(1:dim1,1:dim2),ierr)
SLL_ALLOCATE(es%local_to_global_indices_source_bis(1:dim1,1:dim2),ierr)

SLL_ALLOCATE(local_indices(1:dim1,1:dim2),ierr)
num_splines1 = num_cells1 + spline_degree1
num_splines2 = num_cells2 + spline_degree2
SLL_ALLOCATE(global_indices(num_splines1*num_splines2),ierr)

global_indices          = 0
local_indices           = 0
es%local_to_global_indices = 0
  
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
      
do e = 1, num_cells1*num_cells2
  do b = 1, (spline_degree1+1)*(spline_degree2+1)
    es%local_to_global_indices(b,e) = global_indices(local_indices(b,e))
  end do
end do

deallocate(local_indices)
deallocate(global_indices)

! This should be changed to verify that the passed BC's are part of the
! recognized list described in sll_m_boundary_condition_descriptors...

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
case default
  SLL_ERROR('initialize_general_elliptic_solver','unknown type of gauss points in the direction 2')
end select

!PN : Gauss points positions and weights are computed in [-1:1] interval
!PN : We need to rescale them.
es%gauss_pts1(1,:) = 0.5_f64*es%delta_eta1*(es%gauss_pts1(1,:)+1.0_f64)
es%gauss_pts1(2,:) = 0.5_f64*es%delta_eta1*es%gauss_pts1(2,:)
es%gauss_pts2(1,:) = 0.5_f64*es%delta_eta2*(es%gauss_pts2(1,:)+1.0_f64)
es%gauss_pts2(2,:) = 0.5_f64*es%delta_eta2*es%gauss_pts2(2,:)

knots1_size = 2*spline_degree1+2
knots2_size = 2*spline_degree2+2

SLL_ALLOCATE(es%knots1(knots1_size),ierr)
SLL_ALLOCATE(es%knots2(knots2_size),ierr)
SLL_ALLOCATE(es%knots1_rho(num_cells1 + spline_degree1 + 2),ierr)
SLL_ALLOCATE(es%knots2_rho(num_cells2 + spline_degree2 + 2),ierr)
SLL_ALLOCATE(es%masse(num_cells1*num_cells2),ierr)
SLL_ALLOCATE(es%stiff(num_cells1*num_cells2),ierr)

!AB : We must add plus 1 for the dimension of the 
!AB : solution in the case periodic periodic to 
!AB : include the periodicity in the last point.  

SLL_ALLOCATE(es%rho_vec(num_cells1*num_cells2+1),ierr)
SLL_ALLOCATE(es%phi_vec(num_cells1*num_cells2+1),ierr)
es%rho_vec = 0.0_f64
es%masse   = 0.0_f64
es%stiff   = 0.0_f64

do i = -spline_degree1, spline_degree1+1
  es%knots1(i+spline_degree1+1) = es%delta_eta1*i 
end do
do i = -spline_degree2, spline_degree2+1
  es%knots2(i+spline_degree2+1) = es%delta_eta2*i 
end do
  
es%knots1_rho(1:spline_degree1+1) = eta1_min
es%knots1_rho(num_cells1+2:num_cells1+1+spline_degree1+1) = eta1_max
 
if ( mod(spline_degree1+1,2) == 0 ) then
  do i = spline_degree1+1+1, num_cells1+1
    es%knots1_rho(i) = eta1_min + (i-(spline_degree1+1)/2-1)*es%delta_eta1
  end do
else
  do i = spline_degree1 +1 + 1, num_cells1 + 1
    es%knots1_rho(i) = &
      0.5*(eta1_min+(i-(spline_degree1)/2-1)*es%delta_eta1 + &
      eta1_min+(i-1-(spline_degree1)/2 -1)*es%delta_eta1)
  end do
end if

es%knots2_rho(1:spline_degree2+1) = eta2_min
es%knots2_rho(num_cells2+2:num_cells2+1+spline_degree2+1)=eta2_max
    
if (mod(spline_degree2+1,2) == 0 ) then
  do i = spline_degree2+1+1,num_cells2+1
    es%knots2_rho(i) = eta2_min+(i-(spline_degree2+1)/2-1)*es%delta_eta2 
  end do
else
  do i = spline_degree2+1+1, num_cells2+1
    es%knots2_rho(i) = &
      0.5*( eta2_min+(i  -(spline_degree2)/2-1)*es%delta_eta2 + &
            eta2_min+(i-1-(spline_degree2)/2-1)*es%delta_eta2 )
  end do
end if

! allocation of the table containning 
! all values of splines and its
! derivatives in each gauss points
SLL_ALLOCATE(es%v_splines1(3,spline_degree1+1,spline_degree1+2,num_cells1),ierr)
SLL_ALLOCATE(es%v_splines2(3,spline_degree2+1,spline_degree2+2,num_cells2),ierr)

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
    gspl1 = es%gauss_pts1(1,ii)
    ispl1 = spline_degree1+1
    call bsplvd(es%db,es%knots1,spline_degree1+1,gspl1,ispl1,work1,dbs1,2)
    es%v_splines1(1,:,ii,i) = dbs1(:,1)
    es%v_splines1(2,:,ii,i) = dbs1(:,2)
    call interv(es%db,es%knots1_rho,es%num_cells1+spline_degree1+2,xg,left,ierr)
    call bsplvd(es%db,es%knots1_rho,spline_degree1+1,xg,left,work1,dbs1,1)
    es%v_splines1(3,:,ii,i) = dbs1(:,1)
  end do
  tab_index_coeff1(i) = left
end do

do j = 1, es%num_cells2
  eta2  = eta2_min + (j-1)*es%delta_eta2
  do jj=1,num_pts_g2
    yg  = eta2+es%gauss_pts2(1,jj)
    gspl2 = es%gauss_pts2(1,jj)
    ispl2 = spline_degree2+1
    call bsplvd(es%db,es%knots2,spline_degree2+1,gspl2,ispl2,work2,dbs2,2)
    es%v_splines2(1,:,jj,j) = dbs2(:,1)
    es%v_splines2(2,:,jj,j) = dbs2(:,2)
    call interv(es%db,es%knots2_rho,es%num_cells2+spline_degree2+2,yg,left,ierr)
    call bsplvd(es%db,es%knots2_rho,spline_degree2+1,yg,left,work2,dbs2,1)
    es%v_splines2(3,:,jj,j) = dbs2(:,1)
  end do
  tab_index_coeff2(j) = left
end do

deallocate(work1)
deallocate(work2)
deallocate(dbs1)
deallocate(dbs2)

do j = 1, es%num_cells2
do i = 1, es%num_cells1
  icell = i + (j-1)*es%num_cells1
  b = 0
  do jj = 0, es%spline_degree2
    index3 = j + jj
    if ( index3 > es%num_cells2) index3 = index3 - es%num_cells2
    do ii = 0,es%spline_degree1
      index1 = i + ii
      if ( index1 > es%num_cells1) index1 = index1 - es%num_cells1
      b = b+1
      index_coef1 = tab_index_coeff1(i) - es%spline_degree1 + ii
      index_coef2 = tab_index_coeff2(j) - es%spline_degree2 + jj
      es%local_to_global_indices_source(b,icell)= index_coef1+(index_coef2-1)*(es%num_cells1+1)
      bprime = 0
      do ll = 0,es%spline_degree2
        index4 = j + ll
        if ( index4 > es%num_cells2) index4 = index4 - es%num_cells2
        do kk = 0,es%spline_degree1
          index2 = i + kk
          if ( index2 > es%num_cells1) index2 = index2 - es%num_cells1
          bprime = bprime+1 
          es%local_to_global_indices_source_bis(bprime,icell)=index2+(index4-1)*es%num_cells1
        end do
      end do
    end do
  end do
end do
end do
DEALLOCATE(tab_index_coeff1)
DEALLOCATE(tab_index_coeff2)

end subroutine initialize_gces_full_periodic
  
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
                                      eta1_min,         &
                                      eta1_max,         &
                                      eta2_min,         &
                                      eta2_max ) result(es)

type(sll_gces_full_periodic), pointer :: es

sll_int32,  intent(in) :: spline_degree1
sll_int32,  intent(in) :: spline_degree2
sll_int32,  intent(in) :: num_cells1
sll_int32,  intent(in) :: num_cells2
sll_int32,  intent(in) :: quadrature_type1
sll_int32,  intent(in) :: quadrature_type2
sll_real64, intent(in) :: eta1_min
sll_real64, intent(in) :: eta1_max
sll_real64, intent(in) :: eta2_min
sll_real64, intent(in) :: eta2_max

sll_int32 :: ierr

SLL_ALLOCATE(es,ierr)

call sll_create( es,               &
&                spline_degree1,   &
&                spline_degree2,   &
&                num_cells1,       &
&                num_cells2,       &
&                quadrature_type1, &
&                quadrature_type2, &
&                eta1_min,         &
&                eta1_max,         &
&                eta2_min,         &
&                eta2_max )
   
end function new_general_elliptic_solver

!> @brief Deallocate the type general_coordinate_elliptic_solver
!> @details
!> The parameters are
!> @param[in] es the type general_coordinate_elliptic_solver
  
subroutine delete_gces_full_periodic( es )
type(sll_gces_full_periodic) :: es
sll_int32 :: ierr

SLL_DEALLOCATE(es%knots1,ierr)
SLL_DEALLOCATE(es%knots2,ierr)
SLL_DEALLOCATE(es%gauss_pts1,ierr)
SLL_DEALLOCATE(es%gauss_pts2,ierr)
SLL_DEALLOCATE(es%local_to_global_indices,ierr)
SLL_DEALLOCATE(es%local_to_global_indices_source,ierr)
SLL_DEALLOCATE(es%local_to_global_indices_source_bis,ierr)
call sll_delete(es%csr_mat_with_constraint)
call sll_delete(es%csr_mat_source)
SLL_DEALLOCATE(es%rho_vec,ierr)
SLL_DEALLOCATE(es%phi_vec,ierr)
SLL_DEALLOCATE(es%masse,ierr)
SLL_DEALLOCATE(es%stiff,ierr)
SLL_DEALLOCATE(es%knots1_rho,ierr)
SLL_DEALLOCATE(es%knots2_rho,ierr)
SLL_DEALLOCATE(es%v_splines1,ierr)
SLL_DEALLOCATE(es%v_splines2,ierr)

end subroutine delete_gces_full_periodic


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

type(sll_gces_full_periodic),intent(inout) :: es

type(sll_csr_matrix),            pointer :: csr_mat

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
sll_real64, dimension(:),     allocatable :: stif
sll_real64, dimension(:,:,:), pointer     :: source

sll_int32 :: ierr
sll_int32 :: i
sll_int32 :: j
sll_int32 :: k
sll_int32 :: s
sll_int32 :: icell

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
sll_int32  :: b, bprime,x
sll_int32  :: a, aprime
sll_real64 :: elt_mat_global
sll_int32  :: nspl
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

!$ sll_int32 :: tid=0
!$ sll_int32 :: nthreads=1

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

nspl = (es%spline_degree1+1)*(es%spline_degree2+1)

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

csr_mat => new_csr_matrix( es%num_cells1*es%num_cells2,               &
&                          es%num_cells1*es%num_cells2,               &
&                          es%num_cells1*es%num_cells2,               &
&                          es%local_to_global_indices,                &
&                          nspl,                                      &
&                          es%local_to_global_indices,                &
&                          (es%spline_degree1+1)*(es%spline_degree2+1))

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED( es, c_field, csr_mat, &
!$OMP a11_field_mat, a12_field_mat, a21_field_mat, a22_field_mat, &
!$OMP b1_field_vect, b2_field_vect, spl_deg_1, spl_deg_2, source, intjac ) &
!$OMP FIRSTPRIVATE(nc_1, nc_2, delta1, delta2, eta1_min, eta2_min, &
!$OMP num_pts_g1, num_pts_g2)  &
!$OMP PRIVATE(nthreads, tid, &
!$OMP i,j,eta1,eta2,mass,stif, &
!$OMP yg,wyg,jg,ig,xg,wxg, &
!$OMP wxy,val_c,val_a11,val_a12,val_a21,val_a22, &
!$OMP val_b1, val_b1_der1, val_b1_der2, &
!$OMP val_b2, val_b2_der1, val_b2_der2, &
!$OMP icell, ii, jj, kk, ll, mm, nn,   &
!$OMP jac_mat, val_jac, wxy_by_val_jac, wxy_val_jac, &
!$OMP B11, B12, B21, B22, MC, C1, C2, &
!$OMP v1, v2, v3, v4, r1, r2, d1, d2, d3, d4, &
!$OMP v3v4, d3v4, v3d4, &
!$OMP M_c,K_11,K_12,K_21,K_22,M_bv,S_b1,S_b2, &
!$OMP index1, index2, index3, index4, &
!$OMP a, b, x, aprime, bprime, &
!$OMP r1r2, v1v2, d1v2, v1d2, elt_mat_global )

!$ tid = omp_get_thread_num()
!$ nthreads = omp_get_num_threads()
!$ if (tid == 0) print *, 'Number of threads = ', nthreads
!$OMP DO SCHEDULE(STATIC,nc_2/nthreads) REDUCTION(+:intjac)
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
          stif(mm) = stif(mm)+wxy_val_jac*(d1v2+v1*d2)
               
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
    if ( index3 > es%num_cells2) index3 = index3 - es%num_cells2
     
    do ii = 0,spl_deg_1
        
      index1 = i + ii
      if ( index1 > es%num_cells1) index1 = index1 - es%num_cells1
      x = index1 + (index3-1)*es%num_cells1
      b = ii+1+jj*(spl_deg_1+1)
      a = es%local_to_global_indices(b, icell)
         
      es%masse(x) = es%masse(x) + mass(b)
      es%stiff(x) = es%stiff(x) + stif(b)

      do ll = 0,spl_deg_2
        index4 = j + ll
        if ( index4 > es%num_cells2) index4 = index4 - es%num_cells2
        do kk = 0,spl_deg_1
          index2 = i + kk
          if ( index2 > es%num_cells1) index2 = index2 - es%num_cells1
          bprime = kk+1+ll*(spl_deg_1+1)
          aprime = es%local_to_global_indices(bprime,icell)

          elt_mat_global = M_c (bprime,b) - &
                           K_11(bprime,b) - &
                           K_12(bprime,b) - &
                           K_21(bprime,b) - &
                           K_22(bprime,b) - &
                           M_bv(bprime,b) - &
                           S_b1(bprime,b) - &
                           S_b2(bprime,b)

          if ( a>0 .and. aprime>0 ) then
            call sll_add_to_csr_matrix(csr_mat, elt_mat_global, a, aprime)   
          end if
        end do
      end do
    end do
  end do
end do
end do

!$OMP END PARALLEL

es%intjac = intjac
print *,'#begin of sll_factorize_csr_matrix'

es%csr_mat_with_constraint => new_csr_matrix_with_constraint(csr_mat)

es%csr_mat_with_constraint%num_rows = csr_mat%num_rows+1
es%csr_mat_with_constraint%num_nz   = csr_mat%num_nz+2*csr_mat%num_rows

s = 1
do i=1, csr_mat%num_rows
  es%csr_mat_with_constraint%row_ptr(i) = s
  do k =  csr_mat%row_ptr(i),  csr_mat%row_ptr(i+1)-1
    es%csr_mat_with_constraint%val(s) =  csr_mat%val(k)
    es%csr_mat_with_constraint%col_ind(s) =  csr_mat%col_ind(k)
    s = s+1
  enddo
  es%csr_mat_with_constraint%val(s) = es%masse(i)
  es%csr_mat_with_constraint%col_ind(s) = es%csr_mat_with_constraint%num_rows
  s = s+1
enddo
es%csr_mat_with_constraint%row_ptr(csr_mat%num_rows+1) = s
do i=1, csr_mat%num_rows
  es%csr_mat_with_constraint%val(s) = es%masse(i)
  es%csr_mat_with_constraint%col_ind(s) = i
  s = s+1      
enddo
es%csr_mat_with_constraint%row_ptr(csr_mat%num_rows+2) = s

call sll_delete(csr_mat)
 
call sll_factorize_csr_matrix(es%csr_mat_with_constraint)

print *,'#end of sll_factorize_csr_matrix'

es%csr_mat_source =>                                     &
  new_csr_matrix( size(es%masse,1),                      &
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
           
      elt_mat_global = source(b,bprime,icell)

      if ( a > 0 .and. aprime > 0) then
        call sll_add_to_csr_matrix(es%csr_mat_source,elt_mat_global,a,aprime)
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
!> @param[in]  es  the type general_coordinate_elliptic_solver
!> @param[in]  rho \f$ \rho \f$ the field corresponding to the source term   
!> @param[out] phi \f$ \phi \f$ the field corresponding to the solution of the equation
!> @return     phi the field solution of the equation
  
subroutine solve_gces_full_periodic( es, rho, phi)

class(sll_gces_full_periodic),       intent(inout)     :: es
class(sll_scalar_field_2d_discrete), intent(inout)     :: phi
class(sll_scalar_field_2d_base),     intent(in),target :: rho

class(sll_scalar_field_2d_base), pointer :: base_field_pointer
sll_real64, dimension(:,:),      pointer :: coeff_rho
sll_real64, dimension(:),    allocatable :: m_rho_loc

sll_int32  :: i
sll_int32  :: j
sll_int32  :: num_pts_g1
sll_int32  :: num_pts_g2
sll_int32  :: x, n, b
sll_int32  :: ii, jj, kk, ll, mm, nn
sll_int32  :: index1, index3

sll_real64 :: wgpt1, wgpt2, gpt1, gpt2, eta1, eta2
sll_real64 :: val_f, val_j, valfj, jac_mat(2,2)
sll_real64 :: spline1, spline2

sll_real64 :: int_rho
sll_real64 :: int_jac
sll_int32  :: ierr
sll_int32  :: nc_1, nc_2

!$ sll_int32  :: tid = 0
!$ sll_int32  :: nthreads = 1
  
num_pts_g1 = size(es%gauss_pts1,2)
num_pts_g2 = size(es%gauss_pts2,2)


nc_1            = es%num_cells1
nc_2            = es%num_cells2
es%rho_vec      = 0.0_f64
es%rho_coeff_1d = 0.0_f64
 
base_field_pointer => rho

select type( type_field => base_field_pointer)

class is (sll_scalar_field_2d_discrete)
  
  
  coeff_rho => type_field%interp_2d%get_coefficients()
            
  do j=1,es%num_cells2+1
    do i=1,es%num_cells1+1
      es%rho_coeff_1d(i+(es%num_cells1+1)*(j-1)) = coeff_rho(i,j)
    end do
  end do

  call sll_mult_csr_matrix_vector(es%csr_mat_source,  &
                                  es%rho_coeff_1d,    &
                                  es%rho_vec(1:es%num_cells1*es%num_cells2))

  es%rho_vec(1:es%num_cells1*es%num_cells2) =     &
     es%rho_vec(1:es%num_cells1*es%num_cells2)    &
     - sum(es%rho_vec(1:es%num_cells1*es%num_cells2))/es%intjac*es%masse
      
class is (sll_scalar_field_2d_analytic)
  
  int_rho = 0.0_f64
  int_jac = 0.0_f64

  SLL_CLEAR_ALLOCATE(M_rho_loc(1:(es%spline_degree1+1)*(es%spline_degree2+1)),ierr)

  !$OMP PARALLEL &
  !$OMP FIRSTPRIVATE(nc_1, nc_2, num_pts_g1, num_pts_g2, &
  !$OMP              tid, nthreads)                      &
  !$OMP PRIVATE(i,j,ii,jj,kk,ll,mm,nn,n,m_rho_loc,x,b,   &
  !$OMP         index1,index3,eta1,eta2,gpt1,gpt2,  &
  !$OMP         wgpt1,wgpt2,spline1,spline2,val_f,val_j, &
  !$OMP         valfj,jac_mat)
  !$ tid = omp_get_thread_num()
  !$ nthreads = omp_get_num_threads()
  !$OMP MASTER
  !$ print *, 'Number of threads = ', nthreads
  !$OMP END MASTER
  
  !$OMP DO SCHEDULE(STATIC,nc_2/nthreads) REDUCTION(+:int_rho,int_jac)
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
        if ( index3 > es%num_cells2) index3 = index3 - es%num_cells2
        do nn = 0,es%spline_degree1
          index1 = i + nn
          if (index1 > es%num_cells1) index1=index1-es%num_cells1
      
          x             = index1 + (index3-1)*es%num_cells1
          b             = nn + 1 + mm * (es%spline_degree1+1)
          es%rho_vec(x) = es%rho_vec(x)  + M_rho_loc(b)
            
        end do
      end do
    end do
  end do

  !$OMP END PARALLEL

  es%rho_vec(1:es%num_cells1*es%num_cells2) =  &
    es%rho_vec(1:es%num_cells1*es%num_cells2) - sum(es%rho_vec)/int_jac*es%masse
     
end select

es%phi_vec(:)     = 0.0_f64  !PN: Is it useful ?
  
call sll_solve_csr_matrix(es%csr_mat_with_constraint, &
                          es%rho_vec,                 &
                          es%phi_vec)

print *,'#solve_linear_system done'
  
call phi%interp_2d%set_coefficients(es%phi_vec(1:es%num_cells1*es%num_cells2))

end subroutine solve_gces_full_periodic
  
end module sll_m_gces_full_periodic
