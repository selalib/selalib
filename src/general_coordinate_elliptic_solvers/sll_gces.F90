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
module sll_module_gces
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
use sll_module_arbitrary_degree_spline_interpolator_1d, only:        &
  interv, deboor_type, bsplvd
use gauss_legendre_integration
use gauss_lobatto_integration
use sll_sparse_matrix_module, only : sll_csr_matrix,                 &
                                     new_csr_matrix,                 &
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
type, public :: sll_gces

  private
  sll_int32,  public :: n1
  sll_int32,  public :: n2
  sll_real64, public :: delta_eta1
  sll_real64, public :: delta_eta2
  sll_real64, public :: eta1_min
  sll_real64, public :: eta2_min   

  sll_real64, dimension(:),   pointer :: t1
  sll_real64, dimension(:),   pointer :: t2
  sll_real64, dimension(:),   pointer :: t1_rho
  sll_real64, dimension(:),   pointer :: t2_rho
  sll_real64, dimension(:,:), pointer :: gauss_pts1
  sll_real64, dimension(:,:), pointer :: gauss_pts2
  sll_int32  :: k1
  sll_int32  :: k2
  sll_real64 :: epsi
  sll_real64 :: intjac

  sll_int32, dimension(:,:), pointer :: local_to_global_indices
  sll_int32, dimension(:,:), pointer :: local_to_global_indices_source
  sll_int32, dimension(:,:), pointer :: local_to_global_indices_source_bis

  !!! contains the values of all splines in all gauss points
  sll_real64, dimension(:,:,:,:), pointer :: v_splines1
  sll_real64, dimension(:,:,:,:), pointer :: v_splines2

  sll_real64, dimension(:,:,:,:), pointer :: val_jac
  sll_int32 , dimension(:),       pointer :: tab_index_coeff1
  sll_int32 , dimension(:),       pointer :: tab_index_coeff2
  sll_real64, dimension(:),       pointer :: rho_vec
  sll_real64, dimension(:),       pointer :: tmp_rho_vec
  sll_real64, dimension(:),       pointer :: phi_vec
  sll_real64, dimension(:),       pointer :: masse
  sll_real64, dimension(:),       pointer :: stiff
  sll_real64, dimension(:),       pointer :: rho_coeff_1d
  type(deboor_type)                       :: db

  type(sll_csr_matrix),           pointer :: csr_mat
  type(sll_csr_matrix),           pointer :: csr_mat_with_constraint
  type(sll_csr_matrix),           pointer :: csr_mat_source

end type sll_gces

!> For the integration mode.  
sll_int32, parameter, public :: ES_GAUSS_LEGENDRE = 0
!> For the integration mode.  
sll_int32, parameter, public :: ES_GAUSS_LOBATTO = 1
  
interface sll_delete
  module procedure delete_gces
end interface sll_delete

interface sll_create
  module procedure initialize_gces
end interface sll_create

interface sll_solve
  module procedure solve_gces
end interface sll_solve

public sll_delete,       &
       sll_create,       &
       sll_solve,        &
       new_gces,         &
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
!> @param[in]  k1 the degre of B-spline in the direction eta1
!> @param[in]  k2 the degre of B-spline in the direction eta2
!> @param[in]  n1 the number of cells in the direction eta1
!> @param[in]  n2 the number of cells in the direction eta2
!> @param[in]  quadrature_type1 the type of quadrature in the direction eta1
!> @param[in]  quadrature_type2 the type of quadrature in the direction eta2
!> @param[in]  eta1_min the minimun in the direction eta1
!> @param[in]  eta1_max the maximun in the direction eta1
!> @param[in]  eta2_min the minimun in the direction eta2
!> @param[in]  eta2_max the maximun in the direction eta2
!> @param[out] the type general_coordinate_elliptic_solver
subroutine initialize_gces( es,                  &
                            k1,      &
                            k2,      &
                            n1,          &
                            n2,          &
                            quadrature_type1,    &
                            quadrature_type2,    &
                            eta1_min,            &
                            eta1_max,            &
                            eta2_min,            &
                            eta2_max)
    
type(sll_gces), intent(out) :: es

sll_int32,  intent(in)  :: k1
sll_int32,  intent(in)  :: k2
sll_int32,  intent(in)  :: n1
sll_int32,  intent(in)  :: n2
sll_int32,  intent(in)  :: quadrature_type1
sll_int32,  intent(in)  :: quadrature_type2
sll_real64, intent(in)  :: eta1_min
sll_real64, intent(in)  :: eta1_max
sll_real64, intent(in)  :: eta2_min
sll_real64, intent(in)  :: eta2_max

sll_int32 :: t1_size
sll_int32 :: t2_size
sll_int32 :: vec_sz ! for rho_vec and phi_vec allocations
sll_int32 :: ierr
sll_int32 :: solution_size
sll_int32 :: num_pts_g1, num_pts_g2

sll_real64, allocatable :: work1(:,:)
sll_real64, allocatable :: work2(:,:)
sll_real64, allocatable :: dbs1(:,:)
sll_real64, allocatable :: dbs2(:,:)
sll_real64 :: x1, x2

sll_int32  :: i, j, ii, jj, ispl1, ispl2
sll_real64 :: eta1, eta2, gspl1, gspl2
sll_int32  :: left_rho
sll_int32  :: cell
sll_int32  :: a, b, c, d
sll_int32, dimension(:),   pointer :: global_indices 
sll_int32, dimension(:,:), pointer :: local_indices

SLL_ALLOCATE(es%local_to_global_indices(1:(k1+1)*(k2+1),1:(n1*n2)),ierr)
SLL_ALLOCATE(es%local_to_global_indices_source(1:(k1+1)*(k2+1),1:(n1*n2)),ierr)
SLL_ALLOCATE(es%local_to_global_indices_source_bis(1:(k1+1)*(k2+1),1:(n1*n2)),ierr)

SLL_ALLOCATE(global_indices((n1+k1)*(n2+k2)),ierr)
SLL_ALLOCATE(local_indices(1:(k1+1)*(k2+1),1:(n1*n2)),ierr)
global_indices             = 0
local_indices              = 0
es%local_to_global_indices = 0
  
c = 0
do j = 1, n2    
do i = 1, n1  
  c = c+1
  b = 0
  do jj = 0, k2
  do ii = 0, k1
    b    =  b+1
    local_indices(b, c) = c + (j-1)*k1 + jj*(n1+k1) + ii 
     print"(5i4)", i, j, c, b, local_indices(b, c)
  end do
  end do
end do
end do

d = 0
do j = 2, n2+k2-1
  do i = 2, n1+k1-1
    a = i + (n1+k1)*(j-1)
    d = d + 1
    global_indices(a) = d
  end do
end do

do c = 1, n1*n2
  do b = 1, (k1+1)*(k2+1)
    es%local_to_global_indices(b,c) = global_indices(local_indices(b,c))
  end do
end do


es%k1         = k1
es%k2         = k2
es%n1         = n1
es%n2         = n2
es%delta_eta1 = (eta1_max-eta1_min)/n1
es%delta_eta2 = (eta2_max-eta2_min)/n2
es%eta1_min   = eta1_min
es%eta2_min   = eta2_min

! Allocate and fill the gauss points/weights information.
select case(quadrature_type1)
case (ES_GAUSS_LEGENDRE)
  SLL_ALLOCATE(es%gauss_pts1(2,k1+2),ierr)
  es%gauss_pts1 = gauss_legendre_points_and_weights(k1+2)
case (ES_GAUSS_LOBATTO)
  SLL_ALLOCATE(es%gauss_pts1(2,k1+2),ierr)
  es%gauss_pts1 = gauss_lobatto_points_and_weights(k1+2)
case default
  SLL_ERROR('initialize_general_elliptic_solver','unknown type of gauss points in the direction 1')
end select
   
select case(quadrature_type2)
case (ES_GAUSS_LEGENDRE)
  SLL_ALLOCATE(es%gauss_pts2(2,k2+2),ierr)
  es%gauss_pts2(:,:) = gauss_legendre_points_and_weights(k2+2)
case (ES_GAUSS_LOBATTO)
  SLL_ALLOCATE(es%gauss_pts2(2,k2+2),ierr)
  es%gauss_pts2(:,:) = gauss_lobatto_points_and_weights(k2+2)
case default
  SLL_ERROR('initialize_general_elliptic_solver','unknown type of gauss points in the direction 2')
end select

!PN : Gauss points positions and weights are computed in [-1:1] interval
!PN : We need to rescale them.
es%gauss_pts1(1,:) = 0.5_f64*es%delta_eta1*(es%gauss_pts1(1,:)+1.0_f64)
es%gauss_pts1(2,:) = 0.5_f64*es%delta_eta1*es%gauss_pts1(2,:)
es%gauss_pts2(1,:) = 0.5_f64*es%delta_eta2*(es%gauss_pts2(1,:)+1.0_f64)
es%gauss_pts2(2,:) = 0.5_f64*es%delta_eta2*es%gauss_pts2(2,:)

t1_size = 2*k1 + n1+1
t2_size = 2*k2 + n2+1
vec_sz      = (n1+k1)*(n2+k2)


SLL_ALLOCATE(es%t1(t1_size),ierr)
SLL_ALLOCATE(es%t2(t2_size),ierr)
SLL_ALLOCATE(es%t1_rho(n1+k1+2),ierr)
SLL_ALLOCATE(es%t2_rho(n2+k2+2),ierr)
SLL_ALLOCATE(es%rho_vec(vec_sz),ierr)
SLL_ALLOCATE(es%masse(vec_sz),ierr)
SLL_ALLOCATE(es%stiff(vec_sz),ierr)

!AB : We must add plus 1 for the dimension of the 
!AB : solution in the case periodic periodic to 
!AB : include the periodicity in the last point.  

solution_size = (n1+k1-2)*(n2+k2-2)

SLL_ALLOCATE(es%tmp_rho_vec(solution_size),ierr)
SLL_ALLOCATE(es%phi_vec(solution_size),ierr) 
es%rho_vec = 0.0_f64
es%masse   = 0.0_f64
es%stiff   = 0.0_f64

do i = 1, k1 + 1
  es%t1(i) = eta1_min
enddo
do i = k1+2, n1+1+k1
  es%t1(i) = es%t1(i-1) + es%delta_eta1
enddo
do i = n1+k1+1, n1+1+2*k1
  es%t1(i) = eta1_max
enddo

do j = 1, k2 + 1
  es%t2(j) = eta2_min
enddo
do j = k2+2, n2+1+k2
  es%t2(j) = es%t2(j-1) + es%delta_eta2
enddo
do j = n2+k2+1, n2+1+2*k2
  es%t2(j) = eta2_max
enddo
  
es%csr_mat => new_csr_matrix( solution_size,              &
&                             solution_size,              &
&                             n1*n2,                      &
&                             es%local_to_global_indices, &
&                             (k1+1)*(k2+1),              &
&                             es%local_to_global_indices, &
&                             (k1+1)*(k2+1))

es%t1_rho(1:k1+1) = eta1_min
es%t1_rho(n1+2:n1+1+k1+1) = eta1_max
 
if ( mod(k1+1,2) == 0 ) then
  do i = k1+1+1, n1+1
    es%t1_rho(i) = eta1_min + (i-(k1+1)/2-1)*es%delta_eta1
  end do
else
  do i = k1+1+1, n1+1
    es%t1_rho(i) = 0.5*(eta1_min+(i  -(k1)/2-1)*es%delta_eta1 + &
                            eta1_min+(i-1-(k1)/2-1)*es%delta_eta1)
  end do
end if

es%t2_rho(1:k2+1) = eta2_min
es%t2_rho(n2+2:n2+1+k2+1)=eta2_max
    
if (mod(k2+1,2) == 0 ) then
  do i = k2+1+1, n2+1
    es%t2_rho(i) = eta2_min + (i-(k2+1)/2-1)*es%delta_eta2 
  end do
else
  do i = k2+1+1, n2+1
    es%t2_rho(i) = 0.5*(eta2_min+(i  -(k2)/2-1)*es%delta_eta2 + &
                            eta2_min+(i-1-(k2)/2-1)*es%delta_eta2 )
  end do
end if

! allocation of the table containning 
! all values of splines and its
! derivatives in each gauss points
SLL_ALLOCATE(es%v_splines1(3,k1+1,k1+2,n1),ierr)
SLL_ALLOCATE(es%v_splines2(3,k2+1,k2+2,n2),ierr)

es%v_splines1 = 0.0_f64
es%v_splines2 = 0.0_f64

SLL_ALLOCATE(es%tab_index_coeff1(n1),ierr)
SLL_ALLOCATE(es%tab_index_coeff2(n2),ierr)

num_pts_g1 = size(es%gauss_pts1,2)
num_pts_g2 = size(es%gauss_pts2,2)
SLL_ALLOCATE(es%rho_coeff_1d((n1+1)*(n2+1)),ierr)

allocate(work1(k1+1,k1+1))
allocate(work2(k2+1,k2+1))
allocate(dbs1(k1+1,2))
allocate(dbs2(k2+1,2))

!PN: Compute values of basis functions on the mesh
!PN: We use two functions from Carl DeBoor
!PN:
!PN: interv computes  left = max(i: xt(i) <= x < xt(lxt)).
!PN: 
!PN:   xt  assumed to be nondecreasing
!PN:   x   the point whose location with respect to the 
!PN        sequence xt is to be determined.
!PN: 
!PN: output
!PN:
!PN:  left, mflag.....both integers, whose value is
!PN: 
!PN:  1     -1      if             x <  xt(1)
!PN:  i      0      if   xt(i)  <= x < xt(i+1)
!PN:  i      0      if   xt(i)  <  x .eq. xt(i+1) .eq. xt(lxt)
!PN:  i      1      if   xt(i)  <         xt(i+1) .eq. xt(lxt) .lt. x
!PN: 
!PN:  In particular,  mflag = 0  is the 'usual' case.  mflag .ne. 0
!PN:  indicates that  x  lies outside the CLOSED interval
!PN:  xt(1) .le. y .le. xt(lxt) . The asymmetric treatment of the
!PN:  intervals is due to the decision to make all pp functions cont-
!PN:  inuous from the right, but, by returning  mflag = 0  even if
!PN:  x = xt(lxt), there is the option of having the computed pp function
!PN:  continuous from the left at  xt(lxt) .
!PN: 
!PN:  The program is designed to be efficient in the common situation that
!PN:  it is called repeatedly, with  x  taken from an increasing or decrea-
!PN:  sing sequence. The first guess for left is therefore taken to be the val-
!PN:  ue returned at the previous call and stored in the local variable  ilo . 
!PN   A first check ascertains that  ilo .lt. lxt (this is nec-
!PN:  essary since the present call may have nothing to do with the previ-
!PN:  ous call). Then, if  xt(ilo) .le. x .lt. xt(ilo+1), we set  left =
!PN:  ilo  and are done after just three comparisons.
!PN:     Otherwise, we repeatedly double the difference  istep = ihi - ilo
!PN:  while also moving  ilo  and  ihi  in the direction of  x , until
!PN:                      xt(ilo) .le. x .lt. xt(ihi) ,
!PN:  after which we use bisection to get, in addition, ilo+1 = ihi .
!PN:  left = ilo  is then returned.
!PN: 
!PN:  splvb calculates value and deriv.s of all b-splines which do not vanish at x
!PN: 
!PN:  nput
!PN:  t     the knot array, of length left+k (at least)
!PN:  k     the order of the b-splines to be evaluated
!PN:  x     the point at which these values are sought
!PN:  left  an integer indicating the left endpoint of the interval of
!PN:        interest. the  k  b-splines whose support contains the interval
!PN:               (t(left), t(left+1))
!PN:        are to be considered.
!PN:  assumption  ---  it is assumed that
!PN:               t(left) .lt. t(left+1)
!PN:        division by zero will result otherwise (in  bsplvb ).
!PN:        also, the output is as advertised only if
!PN:               t(left) .le. x .le. t(left+1) .
!PN:  nderiv   an integer indicating that values of b-splines and their
!PN:        derivatives up to but not including the  nderiv-th  are asked
!PN:        for. ( nderiv  is replaced internally by the integer  mhigh
!PN:        in  (1,k)  closest to it.)
!PN: 
!PN:  work area 
!PN:       an array of order (k,k), to contain b-coeff.s of the derivat-
!PN:       ives of a certain order of the  k  b-splines of interest.
!PN: 
!PN:  output
!PN:  dbiatx an array of order (k,nderiv). its entry  (i,m)  contains
!PN:        value of  (m-1)st  derivative of  (left-k+i)-th  b-spline of
!PN:        order  k  for knot sequence  t , i=1,...,k, m=1,...,nderiv.
!PN: 
!PN:  values at  x  of all the relevant b-splines of order k,k-1,...,
!PN:  k+1-nderiv  are generated via  bsplvb  and stored temporarily in
!PN:  dbiatx .  then, the b-coeffs of the required derivatives of the b-
!PN:  splines of interest are generated by differencing, each from the pre-
!PN:  ceding one of lower order, and combined with the values of b-splines
!PN:  of corresponding order in  dbiatx  to produce the desired values .

do i = 1, n1
  eta1  = eta1_min + (i-1)*es%delta_eta1
  do ii=1,num_pts_g1
    x1    = eta1  + es%gauss_pts1(1,ii)
    gspl1 = x1
    ispl1 = k1+i
    call bsplvd(es%db,es%t1,k1+1,gspl1,ispl1,work1,dbs1,2)
    es%v_splines1(1,:,ii,i) = dbs1(:,1)
    es%v_splines1(2,:,ii,i) = dbs1(:,2)
    call interv(es%db,es%t1_rho,n1+k1+2,x1,left_rho,ierr)
    call bsplvd(es%db,es%t1_rho,k1+1,x1,left_rho,work1,dbs1,1)
    es%v_splines1(3,:,ii,i) = dbs1(:,1)
  end do
  es%tab_index_coeff1(i) = left_rho
end do

do j = 1, n2
  eta2  = eta2_min + (j-1)*es%delta_eta2
  do jj=1,num_pts_g2
    x2  = eta2+es%gauss_pts2(1,jj)
    gspl2 = x2
    ispl2 = k2+j
    call bsplvd(es%db,es%t2,k2+1,gspl2,ispl2,work2,dbs2,2)
    es%v_splines2(1,:,jj,j) = dbs2(:,1)
    es%v_splines2(2,:,jj,j) = dbs2(:,2)
    call interv(es%db,es%t2_rho,n2+k2+2,x2,left_rho,ierr)
    call bsplvd(es%db,es%t2_rho,k2+1,x2,left_rho,work2,dbs2,1)
    es%v_splines2(3,:,jj,j) = dbs2(:,1)
  end do
  es%tab_index_coeff2(j) = left_rho
end do

deallocate(work1)
deallocate(work2)
deallocate(dbs1)
deallocate(dbs2)

open(10, file="gces.mtv")
write(10,*)"$DATA=CURVE3D"
write(10,*)"%xmin=", -1.1, " xmax = ", 1.1
write(10,*)"%ymin=", -1.1, " ymax = ", 1.1
write(10,*)"%equalscale=T"
write(10,*)"%toplabel='Elements number ' "
   
do i=1,n1
  do j=1,n2
    write(10,*) eta1_min+(i-1)*es%delta_eta1, eta2_min+(j-1)*es%delta_eta2 , 0.
    write(10,*) eta1_min+(i  )*es%delta_eta1, eta2_min+(j-1)*es%delta_eta2 , 0.
    write(10,*) eta1_min+(i  )*es%delta_eta1, eta2_min+(j  )*es%delta_eta2 , 0.
    write(10,*) eta1_min+(i-1)*es%delta_eta1, eta2_min+(j  )*es%delta_eta2 , 0.
    write(10,*) eta1_min+(i-1)*es%delta_eta1, eta2_min+(j-1)*es%delta_eta2 , 0.
    write(10,*)
  end do
end do

!Numeros des elements
do i=1,n1
  do j=1,n2
    cell = (j-1)*n1 + i
    write(10,"(a)"   ,  advance="no")"@text x1="
    write(10,"(g15.3)", advance="no") eta1_min + (i-0.5)*es%delta_eta1
    write(10,"(a)"   ,  advance="no")" y1="
    write(10,"(g15.3)", advance="no") eta2_min + (j-0.5)*es%delta_eta2
    write(10,"(a)"   ,  advance="no")" z1=0. lc=4 ll='"
    write(10,"(i4)"  ,  advance="no") cell
    write(10,"(a)")"'"
  end do
end do

!Numeros des noeud
do i=1,n2
  do j=1,n1
    cell = (j-1)*n1 + i
    write(10,"(a)"   ,  advance="no")"@text x1="
    write(10,"(g15.3)", advance="no") eta1_min+(i-1)*es%delta_eta1
    write(10,"(a)"   ,  advance="no")" y1="
    write(10,"(g15.3)", advance="no") eta2_min+(j-1)*es%delta_eta2
    write(10,"(a)"   ,  advance="no")" z1=0. lc=5 ll='"
    write(10,"(i4)"  ,  advance="no") local_indices(:,cell)
    write(10,"(a)")"'"
  end do
end do
   
write(10,*)"$END"
close(10)

c = 0
do j=1,n2
  do i=1,n1
  end do
end do

SLL_DEALLOCATE(global_indices,ierr)
SLL_DEALLOCATE(local_indices,ierr)

end subroutine initialize_gces
  
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
!> @param[in] k1 the degre of B-spline in the direction eta1
!> @param[in] k2 the degre of B-spline in the direction eta2
!> @param[in] n1 the number of cells in the direction eta1
!> @param[in] n2 the number of cells in the direction eta2
!> @param[in] quadrature_type1 the type of quadrature in the direction eta1
!> @param[in] quadrature_type2 the type of quadrature in the direction eta2
!> @param[in] eta1_min the minimun in the direction eta1
!> @param[in] eta1_max the maximun in the direction eta1
!> @param[in] eta2_min the minimun in the direction eta2
!> @param[in] eta2_max the maximun in the direction eta2
!> @return    the type general_coordinate_elliptic_solver such that a pointer

function new_gces( k1,               &
                   k2,               &
                   n1,               &
                   n2,               &
                   quadrature_type1, &
                   quadrature_type2, &
                   eta1_min,         &
                   eta1_max,         &
                   eta2_min,         &
                   eta2_max ) result(es)

type(sll_gces), pointer :: es

sll_int32,  intent(in) :: k1
sll_int32,  intent(in) :: k2
sll_int32,  intent(in) :: n1
sll_int32,  intent(in) :: n2
sll_int32,  intent(in) :: quadrature_type1
sll_int32,  intent(in) :: quadrature_type2
sll_real64, intent(in) :: eta1_min
sll_real64, intent(in) :: eta1_max
sll_real64, intent(in) :: eta2_min
sll_real64, intent(in) :: eta2_max

sll_int32 :: ierr

SLL_ALLOCATE(es,ierr)

call sll_create( es,               &
&                k1,               &
&                k2,               &
&                n1,               &
&                n2,               &
&                quadrature_type1, &
&                quadrature_type2, &
&                eta1_min,         &
&                eta1_max,         &
&                eta2_min,         &
&                eta2_max )
   
end function new_gces

!> @brief Deallocate the type general_coordinate_elliptic_solver
!> @details
!> The parameters are
!> @param[in] es the type general_coordinate_elliptic_solver
  
subroutine delete_gces( es )
type(sll_gces) :: es
sll_int32 :: ierr

SLL_DEALLOCATE(es%t1,ierr)
SLL_DEALLOCATE(es%t2,ierr)
SLL_DEALLOCATE(es%gauss_pts1,ierr)
SLL_DEALLOCATE(es%gauss_pts2,ierr)
SLL_DEALLOCATE(es%local_to_global_indices,ierr)
SLL_DEALLOCATE(es%local_to_global_indices_source,ierr)
SLL_DEALLOCATE(es%local_to_global_indices_source_bis,ierr)
SLL_DEALLOCATE(es%rho_vec,ierr)
SLL_DEALLOCATE(es%tmp_rho_vec,ierr)
SLL_DEALLOCATE(es%phi_vec,ierr)
SLL_DEALLOCATE(es%masse,ierr)
SLL_DEALLOCATE(es%stiff,ierr)
SLL_DEALLOCATE(es%t1_rho,ierr)
SLL_DEALLOCATE(es%t2_rho,ierr)
SLL_DEALLOCATE(es%v_splines1,ierr)
SLL_DEALLOCATE(es%v_splines2,ierr)
SLL_DEALLOCATE(es%tab_index_coeff1,ierr)
SLL_DEALLOCATE(es%tab_index_coeff2,ierr)

call sll_delete(es%csr_mat)
call sll_delete(es%csr_mat_with_constraint)
call sll_delete(es%csr_mat_source)

end subroutine delete_gces

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
&                            a11_field_mat, &
&                            a12_field_mat, &
&                            a21_field_mat, &
&                            a22_field_mat, &
&                            b1_field_vect, &
&                            b2_field_vect, &
&                            c_field)

type(sll_gces),intent(inout) :: es

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

sll_int32  :: ierr
sll_int32  :: i
sll_int32  :: j
sll_int32  :: icell
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
sll_real64 :: x1, w1
sll_real64 :: x2, w2
sll_int32  :: index1
sll_int32  :: index2
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
sll_int32  :: index3, index4
sll_int32  :: index_coef1,index_coef2
sll_int32  :: b, bprime,x,y
sll_int32  :: a, aprime
sll_real64 :: elt_mat_global
sll_int32  :: k1, k2, n1, n2
sll_int32  :: ideg2,ideg1
sll_int32  :: jdeg2,jdeg1
sll_real64 :: v1, v2, v3, v4, d1, d2, d3, d4, r1, r2
sll_real64 :: w1w2
sll_real64 :: w1w2_by_val_jac 
sll_real64 :: w1w2_val_jac 
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
k1         = es%k1
k2         = es%k2
n1         = es%n1
n2         = es%n2

intjac = 0.0_f64

SLL_CLEAR_ALLOCATE(source(1:(k1+1)*(k2+1),1:(k1+1)*(k2+1),1:n1*n2),ierr)
SLL_CLEAR_ALLOCATE(M_c(1:(k1+1)*(k2+1),1:(k1+1)*(k2+1)),ierr)
SLL_CLEAR_ALLOCATE(K_11(1:(k1+1)*(k2+1),1:(k1+1)*(k2+1)),ierr)
SLL_CLEAR_ALLOCATE(K_12(1:(k1+1)*(k2+1),1:(k1+1)*(k2+1)),ierr)
SLL_CLEAR_ALLOCATE(K_21(1:(k1+1)*(k2+1),1:(k1+1)*(k2+1)),ierr)
SLL_CLEAR_ALLOCATE(K_22(1:(k1+1)*(k2+1),1:(k1+1)*(k2+1)),ierr)
SLL_CLEAR_ALLOCATE(S_b1(1:(k1+1)*(k2+1),1:(k1+1)*(k2+1)),ierr)
SLL_CLEAR_ALLOCATE(S_b2(1:(k1+1)*(k2+1),1:(k1+1)*(k2+1)),ierr)
SLL_CLEAR_ALLOCATE(M_bv(1:(k1+1)*(k2+1),1:(k1+1)*(k2+1)),ierr)
SLL_CLEAR_ALLOCATE(mass(1:(k1+1)*(k2+1)),ierr)
SLL_CLEAR_ALLOCATE(stif(1:(k1+1)*(k2+1)),ierr)

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED( es, c_field, &
!$OMP a11_field_mat, a12_field_mat, a21_field_mat, a22_field_mat, &
!$OMP b1_field_vect, b2_field_vect, k1, k2, source, intjac ) &
!$OMP FIRSTPRIVATE(n1, n2, delta1, delta2, eta1_min, eta2_min, &
!$OMP num_pts_g1, num_pts_g2)  &
!$OMP PRIVATE(nthreads, tid, &
!$OMP i,j,eta1,eta2,mass,stif, &
!$OMP x2,w2,jg,ig,x1,w1, &
!$OMP w1w2,val_c,val_a11,val_a12,val_a21,val_a22, &
!$OMP val_b1, val_b1_der1, val_b1_der2, &
!$OMP val_b2, val_b2_der1, val_b2_der2, &
!$OMP icell, ii, jj, kk, ll, mm, nn,   &
!$OMP jac_mat, val_jac, w1w2_by_val_jac, w1w2_val_jac, &
!$OMP B11, B12, B21, B22, MC, C1, C2, &
!$OMP v1, v2, v3, v4, r1, r2, d1, d2, d3, d4, &
!$OMP v3v4, d3v4, v3d4, &
!$OMP M_c,K_11,K_12,K_21,K_22,M_bv,S_b1,S_b2, &
!$OMP index_coef2, index_coef1, &
!$OMP index1, index2, index3, index4, &
!$OMP a, b, x, y, aprime, bprime,  &
!$OMP r1r2, v1v2, d1v2, v1d2, elt_mat_global )

!$ tid = omp_get_thread_num()
!$ nthreads = omp_get_num_threads()
!$ if (tid == 0) print *, 'Number of threads = ', nthreads
!$OMP DO SCHEDULE(STATIC,n2/nthreads) REDUCTION(+:intjac)
do j = 1, n2
do i = 1, n1
        
  icell = i + (j-1) * n1
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
  
    x2  = eta2+es%gauss_pts2(1,jg)
    w2 = es%gauss_pts2(2,jg)
  
    do ig=1,num_pts_g1
    
      x1  = eta1+es%gauss_pts1(1,ig)
      w1 = es%gauss_pts1(2,ig)

      w1w2 = w1*w2

      val_c       = c_field%value_at_point(x1,x2)
      val_a11     = a11_field_mat%value_at_point(x1,x2)
      val_a12     = a12_field_mat%value_at_point(x1,x2)
      val_a21     = a21_field_mat%value_at_point(x1,x2)
      val_a22     = a22_field_mat%value_at_point(x1,x2)

      val_b1      = b1_field_vect%value_at_point(x1,x2)
      val_b1_der1 = b1_field_vect%first_deriv_eta1_value_at_point(x1,x2)
      val_b1_der2 = b1_field_vect%first_deriv_eta2_value_at_point(x1,x2)

      val_b2      = b2_field_vect%value_at_point(x1,x2)
      val_b2_der1 = b2_field_vect%first_deriv_eta1_value_at_point(x1,x2)
      val_b2_der2 = b2_field_vect%first_deriv_eta2_value_at_point(x1,x2)
 
      jac_mat = c_field%get_jacobian_matrix(x1,x2)
      val_jac = jac_mat(1,1)*jac_mat(2,2)-jac_mat(1,2)*jac_mat(2,1)

      w1w2_by_val_jac = w1w2/val_jac
      w1w2_val_jac    = w1w2*val_jac

      val_c = val_c * w1w2_val_jac
        
      intjac = intjac + w1w2_val_jac

      ! The B matrix is  by (J^(-1)) A^t (J^(-1))^t 
      B11 = jac_mat(2,2)*jac_mat(2,2)*val_a11 - &
            jac_mat(2,2)*jac_mat(1,2)*(val_a12+val_a21) + &
            jac_mat(1,2)*jac_mat(1,2)*val_a22

      B11 = B11*w1w2_by_val_jac
          
      B21 = jac_mat(1,1)*jac_mat(2,2)*val_a12 - &
            jac_mat(1,1)*jac_mat(1,2)*val_a22 - &
            jac_mat(2,1)*jac_mat(2,2)*val_a11 + &
            jac_mat(1,2)*jac_mat(2,1)*val_a21

      B21 = B21*w1w2_by_val_jac
        
      B12 = jac_mat(1,1)*jac_mat(2,2)*val_a21 - &
            jac_mat(1,1)*jac_mat(1,2)*val_a22 - &
            jac_mat(2,1)*jac_mat(2,2)*val_a11 + &
            jac_mat(1,2)*jac_mat(2,1)*val_a12

      B12 = B12*w1w2_by_val_jac

      B22 = jac_mat(1,1)*jac_mat(1,1)*val_a22 - &
            jac_mat(1,1)*jac_mat(2,1)*(val_a21+val_a12) + &
            jac_mat(2,1)*jac_mat(2,1)*val_a11
          
      B22 = B22*w1w2_by_val_jac

      MC =  jac_mat(2,2)*val_b1_der1 &
          - jac_mat(2,1)*val_b1_der2 &
          - jac_mat(1,2)*val_b2_der1 &
          + jac_mat(1,1)*val_b2_der2

      MC = MC*w1w2
          
      C1 = jac_mat(2,2)*val_b1-jac_mat(1,2)*val_b2 
      C2 = jac_mat(1,1)*val_b2-jac_mat(2,1)*val_b1

      C1 = C1*w1w2
      C2 = C2*w1w2
         
      mm = 0
      do jj = 1,k2+1

        v2 = es%v_splines2(1,jj,jg,j)
        d2 = es%v_splines2(2,jj,jg,j)
        r2 = es%v_splines2(3,jj,jg,j)

        do ii = 1,k1+1
              
          mm = mm+1
                
          v1 = es%v_splines1(1,ii,ig,i)
          d1 = es%v_splines1(2,ii,ig,i)
          r1 = es%v_splines1(3,ii,ig,i)

          r1r2 = r1*r2*w1w2_val_jac
          v1v2 = v1*v2
          d1v2 = d1*v2
          v1d2 = v1*d2

          mass(mm) = mass(mm)+v1v2*w1w2_val_jac 
          stif(mm) = stif(mm)+w1w2_val_jac*(d1v2+v1*d2)
               
          nn = 0
          do ll = 1,k2+1

            v4 = es%v_splines2(1,ll,jg,j)
            d4 = es%v_splines2(2,ll,jg,j)

            do kk = 1,k1+1
                        
              v3 = es%v_splines1(1,kk,ig,i)
              d3 = es%v_splines1(2,kk,ig,i)
              v3v4 = v4*v3
              d3v4 = d3*v4
              v3d4 = v3*d4
              nn = nn+1
             
              source(nn,mm,icell) = source(nn,mm,icell) + r1r2*v3v4
                   
              M_c (nn,mm) = M_c (nn,mm) + val_c * v1v2 * v3v4
              K_11(nn,mm) = K_11(nn,mm) + B11   * d1v2 * d3v4
              K_22(nn,mm) = K_22(nn,mm) + B22   * v1d2 * v3d4
              K_12(nn,mm) = K_12(nn,mm) + B12   * d1v2 * v3d4
              K_21(nn,mm) = K_21(nn,mm) + B21   * v1d2 * d3v4
              M_bv(nn,mm) = M_bv(nn,mm) + MC    * v1v2 * v3v4
              S_b1(nn,mm) = S_b1(nn,mm) + C1    * v1v2 * d3v4
              S_b2(nn,mm) = S_b2(nn,mm) + C2    * v1v2 * v3d4

            end do
          end do
        end do
      end do
    end do
  end do

  b = 0
  do jj = 0, k2

    index3 = j + jj
    
    do ii = 0,k1
        
      index1 = i + ii

      x = index1 + (index3-1)*(n1+k1)
      b = b+1
      a = es%local_to_global_indices(b, icell)
         
      es%masse(x) = es%masse(x) + mass(b)
      es%stiff(x) = es%stiff(x) + stif(b)

      index_coef1 = es%tab_index_coeff1(i) - k1 + ii
      index_coef2 = es%tab_index_coeff2(j) - k2 + jj

      es%local_to_global_indices_source(b,icell)= &
              index_coef1 + (index_coef2-1)*(n1+1)

      bprime = 0
      do ll = 0,k2
             
        index4 = j + ll
             
        do kk = 0,k1
                
          index2 = i + kk

          y      = index2 + (index4-1)*(n1+k1)
          bprime = bprime+1 
          aprime = es%local_to_global_indices(bprime,icell)

          elt_mat_global = M_c (bprime,b) - &
                           K_11(bprime,b) - &
                           K_12(bprime,b) - &
                           K_21(bprime,b) - &
                           K_22(bprime,b) - &
                           M_bv(bprime,b) - &
                           S_b1(bprime,b) - &
                           S_b2(bprime,b)

          es%local_to_global_indices_source_bis(bprime,icell)= y

          if ( a>0 .and. aprime>0 ) then
            call sll_add_to_csr_matrix(es%csr_mat, elt_mat_global, a, aprime)   
          end if
              
        end do
      end do
    end do
  end do
end do
end do

!$OMP END PARALLEL

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
print *,'#begin of sll_factorize_csr_matrix'

call sll_factorize_csr_matrix(es%csr_mat)

print *,'#end of sll_factorize_csr_matrix'

es%csr_mat_source =>                                            &
  new_csr_matrix( size(es%masse,1),                             &
                  (n1+1)*(n2+1),                            &
                  n1*n2,                                    &
                  es%local_to_global_indices_source_bis, &
                  (k1+1)*(k2+1),                                         &
                  es%local_to_global_indices_source,     &
                  (k1+1)*(k2+1) )

icell = 0
do j=1,n2
do i=1,n1
      
  icell = icell+1
  b     = 0
  do ideg2 = 0,k2
  do ideg1 = 0,k1
            
    b = b+1
    a = es%local_to_global_indices_source_bis(b, icell)
        
    bprime = 0
    do jdeg2 = 0,k2
    do jdeg1 = 0,k1
              
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
  
subroutine solve_gces( es, rho, phi)

class(sll_gces), intent(inout)                         :: es
class(sll_scalar_field_2d_discrete), intent(inout)     :: phi
class(sll_scalar_field_2d_base),     intent(in),target :: rho

class(sll_scalar_field_2d_base), pointer :: base_field_pointer
sll_real64, dimension(:,:),      pointer :: coeff_rho
sll_real64, dimension(:), allocatable    :: m_rho_loc

sll_int32  :: i
sll_int32  :: j
sll_int32  :: k
sll_int32  :: num_pts_g1
sll_int32  :: num_pts_g2
sll_int32  :: x, n, b
sll_int32  :: ii, jj, kk, ll, mm, nn
sll_int32  :: index1, index3

sll_real64 :: w1, w2, x1, x2, eta1, eta2
sll_real64 :: val_f, val_j, valfj, jac_mat(2,2)
sll_real64 :: spline1, spline2

sll_real64 :: int_rho
sll_real64 :: int_jac
sll_int32  :: ierr
sll_int32  :: k1, k2, n1, n2

!$ sll_int32  :: tid = 0
!$ sll_int32  :: nthreads = 1
  
num_pts_g1 = size(es%gauss_pts1,2)
num_pts_g2 = size(es%gauss_pts2,2)


k1              = es%k1
k2              = es%k2
n1              = es%n1
n2              = es%n2
es%rho_vec      = 0.0_f64
es%rho_coeff_1d = 0.0_f64
 
base_field_pointer => rho

select type( type_field => base_field_pointer)

class is (sll_scalar_field_2d_discrete)

  coeff_rho => type_field%interp_2d%get_coefficients()
            
  !Loop over point mesh
  do j=1,n2+1
    do i=1,n1+1
      es%rho_coeff_1d(i+(n1+1)*(j-1)) = coeff_rho(i,j)
    end do
  end do

  call sll_mult_csr_matrix_vector(es%csr_mat_source,es%rho_coeff_1d,es%rho_vec)

class is (sll_scalar_field_2d_analytic)
  
  int_rho = 0.0_f64
  int_jac = 0.0_f64

  SLL_CLEAR_ALLOCATE(M_rho_loc(1:(k1+1)*(k2+1)),ierr)

  !$OMP PARALLEL &
  !$OMP FIRSTPRIVATE(n1, n2, num_pts_g1, num_pts_g2, &
  !$OMP              tid, nthreads)                      &
  !$OMP PRIVATE(i,j,ii,jj,kk,ll,mm,nn,n,m_rho_loc,x,b,   &
  !$OMP         index1,index3,eta1,eta2,x1,x2,  &
  !$OMP         w1,w2,spline1,spline2,val_f,val_j, &
  !$OMP         valfj, jac_mat)
  !$ tid = omp_get_thread_num()
  !$ nthreads = omp_get_num_threads()
  !$OMP MASTER
  !$ print *, 'Number of threads = ', nthreads
  !$OMP END MASTER
  
  !$OMP DO SCHEDULE(STATIC,n2/nthreads) REDUCTION(+:int_rho,int_jac)
  do j=1, n2
    do i=1, n1
      M_rho_loc = 0.0_f64
      eta1  = es%eta1_min + (i-1)*es%delta_eta1
      eta2  = es%eta2_min + (j-1)*es%delta_eta2
      do jj=1,num_pts_g2
        x2  = eta2 + es%gauss_pts2(1,jj)
        w2 = es%gauss_pts2(2,jj)
        do ii=1,num_pts_g1
          x1      = eta1 + es%gauss_pts1(1,ii)
          w1      = es%gauss_pts1(2,ii)
      
          val_f   = rho%value_at_point(x1,x2)
          jac_mat = rho%get_jacobian_matrix(x1,x2)
          val_j   = jac_mat(1,1)*jac_mat(2,2)-jac_mat(1,2)*jac_mat(2,1)
          val_j   = val_j*w1*w2
          valfj   = val_f*val_j
          int_rho = int_rho + valfj 
          int_jac = int_jac + val_j
      
          do ll = 1,k2+1
            spline2 = es%v_splines2(1,ll,jj,j)*valfj
            do kk = 1,k1+1
              spline1 = es%v_splines1(1,kk,ii,i)
              n = kk + (ll-1)*(k1+1)
              M_rho_loc(n)= M_rho_loc(n) + spline1*spline2
            end do
          end do
        end do
      end do

      do mm = 0,k2
        index3 = j + mm
        do nn = 0,k1
          index1 = i + nn
          x             =  index1 + (index3-1)*(n1+k1)
          b             =  nn + 1 + mm * (k1+1)
          es%rho_vec(x) =  es%rho_vec(x)  + M_rho_loc(b)
            
        end do
      end do
    end do
  end do

  !$OMP END PARALLEL

end select

es%tmp_rho_vec(:) = 0.0_f64  !PN: Is it useful ?
es%phi_vec(:)     = 0.0_f64  !PN: Is it useful ?
  
k = 0
do j = 1, n1+k1-2
  do i = 1, n2+k2-2
    k = k+1
    es%tmp_rho_vec(k) = es%rho_vec(i+1+(n1+k1)*j)
  end do
end do
     
call sll_solve_csr_matrix(es%csr_mat, es%tmp_rho_vec, es%phi_vec)
  
call phi%interp_2d%set_coefficients(es%phi_vec(1:(n1+k1-2)*(n2+k2-2)))

end subroutine solve_gces
  

end module sll_module_gces
