!> @ingroup poisson2d curvilinear solver
!> @brief Poisson solver on 2d curvilinear mesh
!> @details This solver works with geometry defined on gauss points

module sll_m_poisson_2d_curvilinear
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_memory.h"
#include "sll_working_precision.h"

  use sll_m_boundary_condition_descriptors, only: &
    sll_p_dirichlet, &
    sll_p_periodic

  use sll_m_gauss_legendre_integration, only: &
    sll_f_gauss_legendre_points_and_weights

  use sll_m_gauss_lobatto_integration, only: &
    sll_f_gauss_lobatto_points_and_weights

  use sll_m_sparse_matrix, only: &
    sll_f_new_csr_matrix, &
    sll_t_csr_matrix, &
    sll_s_factorize_csr_matrix

  implicit none

  public :: &
    sll_f_new_poisson_2d_curvilinear, &
    sll_t_poisson_2d_curvilinear, &
    sll_p_poisson_gauss_legendre, &
    sll_p_poisson_open_knots, &
    sll_p_poisson_periodic_knots

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

type :: sll_t_poisson_2d_curvilinear

  private
  sll_int32 :: num_cells1
  sll_int32 :: num_cells2
  sll_real64, dimension(:), pointer :: node_positions1
  sll_real64, dimension(:), pointer :: node_positions2
  sll_real64, dimension(:), pointer, public :: quadrature_points1
  sll_real64, dimension(:), pointer :: quadrature_points2
  sll_real64, dimension(:), pointer, public :: quadrature_weights1
  sll_real64, dimension(:), pointer :: quadrature_weights2
  sll_int32, public :: num_quadrature_points1
  sll_int32 :: num_quadrature_points2  
  sll_int32 :: spline_degree1
  sll_int32 :: spline_degree2
  sll_real64, dimension(:), pointer :: knots1
  sll_real64, dimension(:), pointer :: knots2
  sll_int32, dimension(:,:), pointer :: local_spline_indices
  sll_int32, dimension(:), pointer :: bc_indices1
  sll_int32, dimension(:), pointer :: bc_indices2
  sll_int32, dimension(:), pointer, public :: global_spline_indices
  sll_int32, dimension(:,:), pointer :: local_to_global_spline_indices
  sll_int32 :: solution_size1
  sll_int32 :: solution_size2
  type(sll_t_csr_matrix), pointer :: sll_csr_mat




end type sll_t_poisson_2d_curvilinear

!> For the integration mode.  
sll_int32, parameter :: sll_p_poisson_gauss_legendre = 0
!> For the integration mode.  
sll_int32, parameter :: POISSON_GAUSS_LOBATTO = 1

!> For the knots defined on the boundary  
sll_int32, parameter :: sll_p_poisson_periodic_knots = 0
!> For the knots defined on the boundary  
sll_int32, parameter :: sll_p_poisson_open_knots = 1




! *******************************************************************

contains 


!> @brief Initialization for poisson 2d curvilinear solver.
!> @details 
!> The parameters are
!> @param[in] spline_degree1 the degre of B-spline in the direction eta1
!> @param[in] spline_degree2 the degre of B-spline in the direction eta2
!> @param[in] num_cells1 the number of cells in the direction eta1
!> @param[in] num_cells2 the number of cells in the direction eta2
!> @param[in] bc_min1 the boundary condition at eta1_min
!> @param[in] bc_max1 the boundary condition at eta1_max
!> @param[in] bc_min2 the boundary condition at eta2_min
!> @param[in] bc_max2 the boundary condition at eta2_max
!> @param[in] eta1_min the minimun in the direction eta1
!> @param[in] eta1_max the maximun in the direction eta1
!> @param[in] eta2_min the minimun in the direction eta2
!> @param[in] eta2_max the maximun in the direction eta2
!> @param[in] node_positions1 the node positions in direction eta1
!> @param[in] node_positions2 the node positions in direction eta2
!> @param[in] quadrature_type1 the type of quadrature in the direction eta1
!> @param[in] quadrature_type2 the type of quadrature in the direction eta2
!> @param[in] quadrature_points1 the quadrature points 
!>            normalized to [0,1] in the direction eta1
!> @param[in] quadrature_points2 the quadrature points 
!>            normalized to [0,1] in the direction eta2
!> @param[in] quadrature_weights1 the quadrature weights 
!>            normalized to [0,1] in the direction eta1
!> @param[in] quadrature_weights2 the quadrature weights 
!>            normalized to [0,1] in the direction eta2
!> @param[in] num_quadrature_points1 the number of quadrature weights in direction eta1
!> @param[in] num_quadrature_points2 the number of quadrature weights in direction eta2
!> @param[in] bc_knots_min1 type of boundary condition at eta1_min for knots dir eta1
!> @param[in] bc_knots_max1 type of boundary condition at eta1_max for knots dir eta1
!> @param[in] bc_knots_min2 type of boundary condition at eta2_min for knots dir eta2
!> @param[in] bc_knots_max2 type of boundary condition at eta2_max for knots dir eta2
!> @param[in] knots_min1 knots at eta1_min in direction eta1
!> @param[in] knots_max1 knots at eta1_max in direction eta1
!> @param[in] knots_min2 knots at eta2_min in direction eta2
!> @param[in] knots_max2 knots at eta2_max in direction eta2
!> @return the type sll_t_poisson_2d_curvilinear as a pointer

function sll_f_new_poisson_2d_curvilinear( &
  spline_degree1, &
  spline_degree2, &
  num_cells1, &
  num_cells2, &  
  bc_min1, &
  bc_max1, &
  bc_min2, &
  bc_max2, &
  eta1_min, &
  eta1_max, &
  eta2_min, &
  eta2_max, &
  node_positions1, &
  node_positions2, &
  quadrature_type1, &
  quadrature_type2, &
  quadrature_points1, &
  quadrature_points2, &
  quadrature_weights1, &
  quadrature_weights2, &
  num_quadrature_points1, &
  num_quadrature_points2, &
  bc_knots_min1, &
  bc_knots_max1, &
  bc_knots_min2, &
  bc_knots_max2, &
  knots_min1, &
  knots_max1, &
  knots_min2, &
  knots_max2 &  
  ) result(poisson)

  type(sll_t_poisson_2d_curvilinear), pointer :: poisson
  sll_int32,  intent(in) :: spline_degree1
  sll_int32,  intent(in) :: spline_degree2
  sll_int32,  intent(in) :: num_cells1
  sll_int32,  intent(in) :: num_cells2
  sll_int32,  intent(in) :: bc_min1
  sll_int32,  intent(in) :: bc_max1
  sll_int32,  intent(in) :: bc_min2
  sll_int32,  intent(in) :: bc_max2
  sll_real64, intent(in), optional :: eta1_min
  sll_real64, intent(in), optional :: eta1_max
  sll_real64, intent(in), optional :: eta2_min
  sll_real64, intent(in), optional :: eta2_max
  sll_real64, dimension(:), intent(in), optional :: node_positions1    
  sll_real64, dimension(:), intent(in), optional :: node_positions2    
  sll_int32,  intent(in), optional :: quadrature_type1
  sll_int32,  intent(in), optional :: quadrature_type2
  sll_real64, dimension(:), intent(in), optional :: quadrature_points1
  sll_real64, dimension(:), intent(in), optional :: quadrature_points2
  sll_real64, dimension(:), intent(in), optional :: quadrature_weights1
  sll_real64, dimension(:), intent(in), optional :: quadrature_weights2
  sll_int32, intent(in), optional :: num_quadrature_points1
  sll_int32, intent(in), optional :: num_quadrature_points2
  sll_int32, intent(in), optional :: bc_knots_min1
  sll_int32, intent(in), optional :: bc_knots_max1
  sll_int32, intent(in), optional :: bc_knots_min2
  sll_int32, intent(in), optional :: bc_knots_max2
  sll_real64, dimension(:), intent(in), optional :: knots_min1
  sll_real64, dimension(:), intent(in), optional :: knots_max1
  sll_real64, dimension(:), intent(in), optional :: knots_min2
  sll_real64, dimension(:), intent(in), optional :: knots_max2
  
  sll_int32 :: ierr


  SLL_ALLOCATE(poisson,ierr)
  call initialize_poisson_2d_curvilinear( &
    poisson, &
    spline_degree1, &
    spline_degree2, &
    num_cells1, &
    num_cells2, &
    bc_min1, &
    bc_max1, &
    bc_min2, &
    bc_max2, &
    eta1_min, &
    eta1_max, &
    eta2_min, &
    eta2_max, &     
    node_positions1, &
    node_positions2, &
    quadrature_type1, &
    quadrature_type2, &
    quadrature_points1, &
    quadrature_points2, &
    quadrature_weights1, &
    quadrature_weights2, &
    num_quadrature_points1, &
    num_quadrature_points2, &
    bc_knots_min1, &
    bc_knots_max1, &
    bc_knots_min2, &
    bc_knots_max2, &
    knots_min1, &
    knots_max1, &
    knots_min2, &
    knots_max2 )
   
end function sll_f_new_poisson_2d_curvilinear



! *******************************************************************

!> @brief Initialization of poisson 2d curvilinear solver
!> @details 
!> The parameters are
!> @param[out] poisson the type general_coordinate_elliptic_solver
!> @param[in] spline_degree1 the degre of B-spline in the direction eta1
!> @param[in] spline_degree2 the degre of B-spline in the direction eta2
!> @param[in] num_cells1 the number of cells in the direction eta1
!> @param[in] num_cells2 the number of cells in the direction eta2
!> @param[in] bc_min1 the boundary condition at eta1_min
!> @param[in] bc_max1 the boundary condition at eta1_max
!> @param[in] bc_min2 the boundary condition at eta2_min
!> @param[in] bc_max2 the boundary condition at eta2_max
!> @param[in] eta1_min the minimun in the direction eta1
!> @param[in] eta1_max the maximun in the direction eta1
!> @param[in] eta2_min the minimun in the direction eta2
!> @param[in] eta2_max the maximun in the direction eta2
!> @param[in] node_positions1 the node positions in direction eta1
!> @param[in] node_positions2 the node positions in direction eta2
!> @param[in] quadrature_type1 the type of quadrature in the direction eta1
!> @param[in] quadrature_type2 the type of quadrature in the direction eta2
!> @param[in] quadrature_points1 the quadrature points 
!>            normalized to [0,1] in the direction eta1
!> @param[in] quadrature_points2 the quadrature points 
!>            normalized to [0,1] in the direction eta2
!> @param[in] quadrature_weights1 the quadrature weights 
!>            normalized to [0,1] in the direction eta1
!> @param[in] quadrature_weights2 the quadrature weights 
!>            normalized to [0,1] in the direction eta2
!> @param[in] num_quadrature_points1 the number of quadrature weights in direction eta1
!> @param[in] num_quadrature_points2 the number of quadrature weights in direction eta2
!> @param[in] bc_knots_min1 type of boundary condition at eta1_min for knots dir eta1
!> @param[in] bc_knots_max1 type of boundary condition at eta1_max for knots dir eta1
!> @param[in] bc_knots_min2 type of boundary condition at eta2_min for knots dir eta2
!> @param[in] bc_knots_max2 type of boundary condition at eta2_max for knots dir eta2
!> @param[in] knots_min1 knots at eta1_min in direction eta1
!> @param[in] knots_max1 knots at eta1_max in direction eta1
!> @param[in] knots_min2 knots at eta2_min in direction eta2
!> @param[in] knots_max2 knots at eta2_max in direction eta2
subroutine initialize_poisson_2d_curvilinear( &
  poisson, &
  spline_degree1, &
  spline_degree2, &
  num_cells1, &
  num_cells2, &  
  bc_min1, &
  bc_max1, &
  bc_min2, &
  bc_max2, &
  eta1_min, &
  eta1_max, &
  eta2_min, &
  eta2_max, &
  node_positions1, &
  node_positions2, &
  quadrature_type1, &
  quadrature_type2, &
  quadrature_points1, &
  quadrature_points2, &
  quadrature_weights1, &
  quadrature_weights2, &
  num_quadrature_points1, &
  num_quadrature_points2, &
  bc_knots_min1, &
  bc_knots_max1, &
  bc_knots_min2, &
  bc_knots_max2, &
  knots_min1, &
  knots_max1, &
  knots_min2, &
  knots_max2 )
    
  type(sll_t_poisson_2d_curvilinear), intent(out) :: poisson
  sll_int32,  intent(in) :: spline_degree1
  sll_int32,  intent(in) :: spline_degree2
  sll_int32,  intent(in) :: num_cells1
  sll_int32,  intent(in) :: num_cells2
  sll_int32,  intent(in) :: bc_min1
  sll_int32,  intent(in) :: bc_max1
  sll_int32,  intent(in) :: bc_min2
  sll_int32,  intent(in) :: bc_max2
  sll_real64, intent(in), optional :: eta1_min
  sll_real64, intent(in), optional :: eta1_max
  sll_real64, intent(in), optional :: eta2_min
  sll_real64, intent(in), optional :: eta2_max
  sll_real64, dimension(:), intent(in), optional :: node_positions1    
  sll_real64, dimension(:), intent(in), optional :: node_positions2    
  sll_int32,  intent(in), optional :: quadrature_type1
  sll_int32,  intent(in), optional :: quadrature_type2
  sll_real64, dimension(:), intent(in), optional :: quadrature_points1
  sll_real64, dimension(:), intent(in), optional :: quadrature_points2
  sll_real64, dimension(:), intent(in), optional :: quadrature_weights1
  sll_real64, dimension(:), intent(in), optional :: quadrature_weights2
  sll_int32, intent(in), optional :: num_quadrature_points1
  sll_int32, intent(in), optional :: num_quadrature_points2
  sll_int32, intent(in), optional :: bc_knots_min1
  sll_int32, intent(in), optional :: bc_knots_max1
  sll_int32, intent(in), optional :: bc_knots_min2
  sll_int32, intent(in), optional :: bc_knots_max2
  sll_real64, dimension(:), intent(in), optional :: knots_min1
  sll_real64, dimension(:), intent(in), optional :: knots_max1
  sll_real64, dimension(:), intent(in), optional :: knots_min2
  sll_real64, dimension(:), intent(in), optional :: knots_max2
  
  sll_int32 :: ierr
  
  sll_int32 :: dim1
  sll_int32 :: dim2
  sll_int32 :: num_spl1
  sll_int32 :: num_spl2
  sll_int32 :: solution_size
  sll_int32 :: total_num_splines_loc
  
  !temporary
  !sll_real64, dimension(:,:), allocatable :: rhs_at_points
  !sll_real64, dimension(:,:), allocatable :: splines_at_points1
  !sll_real64, dimension(:,:), allocatable :: splines_at_points2
  !sll_real64, dimension(:,:), allocatable :: output
  
  sll_int32 :: num_nz
  
  
  poisson%spline_degree1 = spline_degree1
  poisson%spline_degree2 = spline_degree2
  
  poisson%num_cells1 = num_cells1
  poisson%num_cells2 = num_cells2
  
  poisson%node_positions1 => new_node_positions( &
    num_cells1, &
    node_positions1, &
    eta1_min, &
    eta1_max)

  poisson%node_positions2 => new_node_positions( &
    num_cells2, &
    node_positions2, &
    eta2_min, &
    eta2_max)

  call new_quadrature( &
    spline_degree1, &
    quadrature_type1, &
    quadrature_points1, &
    quadrature_weights1, &
    num_quadrature_points1, &
    poisson%quadrature_points1, &
    poisson%quadrature_weights1, &
    poisson%num_quadrature_points1)  

  call new_quadrature( &
    spline_degree2, &
    quadrature_type2, &
    quadrature_points2, &
    quadrature_weights2, &
    num_quadrature_points2, &
    poisson%quadrature_points2, &
    poisson%quadrature_weights2, &
    poisson%num_quadrature_points2)  

  poisson%knots1 => new_knots( &
    poisson%node_positions1, &
    num_cells1, &
    spline_degree1, &
    bc_knots_min1, &
    bc_knots_max1, &
    knots_min1, &
    knots_max1 )  
    
  poisson%knots2 => new_knots( &
    poisson%node_positions2, &
    num_cells2, &
    spline_degree2, &
    bc_knots_min2, &
    bc_knots_max2, &
    knots_min2, &
    knots_max2 )  




  dim1 = (spline_degree1+1)*(spline_degree2+1)
  dim2 = num_cells1*num_cells2
  SLL_ALLOCATE(poisson%local_spline_indices(1:dim1,1:dim2),ierr)
  poisson%local_spline_indices=0
  SLL_ALLOCATE(poisson%local_to_global_spline_indices(1:dim1,1:dim2),ierr)
  poisson%local_to_global_spline_indices=0
  num_spl1=num_cells1+spline_degree1
  num_spl2=num_cells2+spline_degree2
  SLL_ALLOCATE(poisson%global_spline_indices(num_spl1*num_spl2),ierr)
  poisson%global_spline_indices=0

  call compute_local_splines_indices( &
    num_cells1, &
    num_cells2, &
    spline_degree1, &
    spline_degree2, &
    poisson%local_spline_indices)

  SLL_ALLOCATE(poisson%bc_indices1(num_spl1),ierr)
  SLL_ALLOCATE(poisson%bc_indices2(num_spl2),ierr)

  call compute_bc_indices(num_cells1,spline_degree1,bc_min1,bc_max1,poisson%bc_indices1)
  call compute_bc_indices(num_cells2,spline_degree2,bc_min2,bc_max2,poisson%bc_indices2)
  call compute_global_splines_indices( &
    num_spl1, &
    num_spl2, &
    poisson%bc_indices1, &
    poisson%bc_indices2, &
    poisson%global_spline_indices)
  call compute_local_to_global_splines_indices( &
    num_cells1, &
    num_cells2, &
    spline_degree1, &
    spline_degree2, &
    poisson%local_spline_indices, &
    poisson%global_spline_indices, &
    poisson%local_to_global_spline_indices)

  poisson%solution_size1 = compute_solution_size( &
    bc_min1, &
    bc_max1, &
    num_cells1, &
    spline_degree1)
  poisson%solution_size2 = compute_solution_size( &
    bc_min2, &
    bc_max2, &
    num_cells2, &
    spline_degree2)
  
  solution_size = poisson%solution_size1*poisson%solution_size2
  total_num_splines_loc = (spline_degree1+1)*(spline_degree2+1)

  num_nz = compute_nz( &
    num_cells1, &
    num_cells2, &
    spline_degree1, &
    spline_degree2, &
    poisson%global_spline_indices)
 
  print*,'#num_nz=',num_nz

  poisson%sll_csr_mat => sll_f_new_csr_matrix( &
    solution_size, &
    solution_size, &
    num_cells1*num_cells2, &
    poisson%local_to_global_spline_indices, &
    total_num_splines_loc, &
    poisson%local_to_global_spline_indices, &
    total_num_splines_loc)


  print *,'#begin of sll_s_factorize_csr_matrix'
  call sll_s_factorize_csr_matrix(poisson%sll_csr_mat)
  print *,'#end of sll_s_factorize_csr_matrix'

!  if (es%perper) then
!
!   es%sll_csr_mat_with_constraint => sll_f_new_csr_matrix_with_constraint(es%sll_csr_mat)  
!
!
!   call sll_s_csr_add_one_constraint( &
!    es%sll_csr_mat%row_ptr, & 
!    es%sll_csr_mat%col_ind, &
!    es%sll_csr_mat%val, &
!    es%sll_csr_mat%num_rows, &
!    es%sll_csr_mat%num_nz, &
!    es%masse, &
!    es%sll_csr_mat_with_constraint%row_ptr, &
!    es%sll_csr_mat_with_constraint%col_ind, &
!    es%sll_csr_mat_with_constraint%val)  
!
!    print *,'#begin of sll_s_factorize_csr_matrix'
!    call sll_s_factorize_csr_matrix(es%sll_csr_mat_with_constraint)
!    print *,'#end of sll_s_factorize_csr_matrix'
!  else   
!    print *,'#begin of sll_s_factorize_csr_matrix'
!    call sll_s_factorize_csr_matrix(es%sll_csr_mat)
!    print *,'#end of sll_s_factorize_csr_matrix'
!  end if 

!
!
!!  call compute_composite_quadrature_points( &
!!    eta1_min, &
!!    eta1_max, &
!!    num_cells1, &
!!    poisson%gauss_pts1(1,:), &
!!    spline_degree1+2, &
!!    -1._f64, &
!!    1._f64, &
!!    node_positions)
!!
!!
!!  call compute_dbiatx( &
!!    knots, &
!!    knot_size, &
!!    node_positions, &
!!    num_points, &  
!!    spline_degree, &
!!    bc, &
!!    dbiatx, &
!!    dbiatx_indices, &
!!    num_deriv)
!
!  SLL_ALLOCATE(rhs_at_points((spline_degree1+2)*num_cells1,(spline_degree2+2)*num_cells2),ierr)
!  SLL_ALLOCATE(splines_at_points1(spline_degree1+1,(spline_degree1+2)*num_cells1),ierr)
!  SLL_ALLOCATE(splines_at_points2(spline_degree2+1,(spline_degree2+2)*num_cells2),ierr)
!  SLL_ALLOCATE(output(num_cells1+spline_degree1+1,num_cells2+spline_degree2+1),ierr)
!
!
!  rhs_at_points = 1._f64
!  
!  splines_at_points1 = 0.5_f64
!  splines_at_points2 = 0.6_f64
!
!  call compute_source_at_points( &
!    rhs_at_points, &
!    splines_at_points1, &
!    splines_at_points2, &
!    num_cells1, &
!    num_cells2, &
!    spline_degree1+2, &
!    spline_degree2+2, &
!    spline_degree1, &
!    spline_degree2, &
!    output)
!  
!  print *,'#output=',output(1,num_cells2+spline_degree2),maxval(abs(output))

end subroutine initialize_poisson_2d_curvilinear

function new_node_positions( &
  num_cells, &
  node_positions, &  
  eta_min, &
  eta_max) result(res)
  sll_real64, dimension(:), pointer :: res
  sll_int32, intent(in) :: num_cells
  sll_real64, dimension(:), intent(in), optional :: node_positions
  sll_real64, intent(in), optional :: eta_min
  sll_real64, intent(in), optional :: eta_max
  
  sll_int32 :: ierr
  sll_int32 :: i
  sll_real64 :: eta_min_value
  sll_real64 :: eta_max_value
  sll_real64 :: delta

  SLL_ALLOCATE(res(num_cells+1),ierr)
  
  if(present(node_positions))then
    if(size(node_positions)<num_cells+1)then
      print *,'#bad size for node_positions=',size(node_positions)
      print *,'#num_cells+1=',num_cells+1
      print *,'#at line/file ',__LINE__,__FILE__
      stop
    endif
    res(1:num_cells+1) = node_positions(1:num_cells+1)
    if(present(eta_min)) then
      if(eta_min /= node_positions(1))then
        print *,'#Warning bad value for eta_min=',eta_min
        print *,'#should be equal to node_positions(1)',node_positions(1)
        print *,'#at line/file ',__LINE__,__FILE__
        print *,'#eta1_min is not used'       
      endif
    endif
    if(present(eta_max)) then
      if(eta_max /= node_positions(num_cells+1))then
        print *,'#Warning bad value for eta_max=',eta_max
        print *,'#should be equal to node_positions(num_cells+1)', &
          node_positions(num_cells+1)
        print *,'#at line/file ',__LINE__,__FILE__        
        print *,'#eta_max is not used'       
      endif
    endif    
  else
    eta_min_value = 0._f64
    if(present(eta_min))then
      eta_min_value = eta_min
    endif
    eta_max_value = 1._f64  
    if(present(eta_max))then
      eta_max_value = eta_max
    endif
    delta = (eta_max_value-eta_min_value)/real(num_cells,f64)
    do i=1,num_cells+1
      res(i) = eta_min_value+real(i-1,f64)*delta
    enddo  
  endif

end function new_node_positions
  

subroutine new_quadrature( &
  degree, &
  quadrature_type, &
  quadrature_points, &
  quadrature_weights, &
  num_quadrature_points, &
  quadrature_points_out, &
  quadrature_weights_out, &
  num_quadrature_points_out)  
  sll_int32, intent(in) :: degree
  sll_int32,  intent(in), optional :: quadrature_type
  sll_real64, dimension(:), intent(in), optional :: quadrature_points
  sll_real64, dimension(:), intent(in), optional :: quadrature_weights
  sll_int32, intent(in), optional :: num_quadrature_points
  sll_real64, dimension(:), pointer :: quadrature_points_out
  sll_real64, dimension(:), pointer :: quadrature_weights_out
  sll_int32, intent(out) :: num_quadrature_points_out
  
  sll_real64, dimension(:,:), allocatable :: pts_and_weights
  sll_int32 :: ierr
  sll_int32 :: quadrature_type_value
  
  num_quadrature_points_out = degree+2
  if(present(num_quadrature_points))then
    num_quadrature_points_out = num_quadrature_points
  endif
  
  
  SLL_ALLOCATE(quadrature_points_out(num_quadrature_points_out),ierr)
  SLL_ALLOCATE(quadrature_weights_out(num_quadrature_points_out),ierr)

  if((present(quadrature_points).or.present(quadrature_weights))) then
    
    if(.not.(present(quadrature_points).and.present(quadrature_weights)))then
      print *,'#quadrature_points and quadrature_weights'
      print *,'#should be provided'
      print *,'#at line/file',__LINE__,__FILE__
      stop
    endif
    if(size(quadrature_points)<num_quadrature_points_out) then
      print *,'#bad size for quadrature_points'
      print *,'#at line/file',__LINE__,__FILE__
    endif
    if(size(quadrature_weights)<num_quadrature_points_out) then
      print *,'#bad size for quadrature_weights'
      print *,'#at line/file',__LINE__,__FILE__
      stop
    endif
    quadrature_points_out(1:num_quadrature_points_out) &
      = quadrature_points(1:num_quadrature_points_out)
    quadrature_weights_out(1:num_quadrature_points_out) &
      = quadrature_weights(1:num_quadrature_points_out)
    if(present(quadrature_type))then
      print *,'#Warning: quadrature_type is ignored',quadrature_type
      print *,'#as quadrature_points and quadrature_weights'
      print *,'#are provided'
    endif
  else
    quadrature_type_value = sll_p_poisson_gauss_legendre
    if(present(quadrature_type))then
      quadrature_type_value = quadrature_type
    endif
    SLL_ALLOCATE(pts_and_weights(2,num_quadrature_points_out),ierr)
    select case(quadrature_type_value)
      case (sll_p_poisson_gauss_legendre)
        pts_and_weights(1:2,1:num_quadrature_points_out) = &
          sll_f_gauss_legendre_points_and_weights(num_quadrature_points_out,0._f64,1._f64)
      case (POISSON_GAUSS_LOBATTO)
        pts_and_weights(1:2,1:num_quadrature_points_out) = &
          sll_f_gauss_lobatto_points_and_weights(num_quadrature_points_out,0._f64,1._f64)
      case default
        print *, '#bad choice for quadrature_type'
        print *,'#error at line',__LINE__,'in file',__FILE__
        stop       
    end select
    quadrature_points_out(1:num_quadrature_points_out) &
      = pts_and_weights(1,1:num_quadrature_points_out)
    quadrature_weights_out(1:num_quadrature_points_out) &
      = pts_and_weights(2,1:num_quadrature_points_out)
    
  endif  

  

end subroutine new_quadrature


function new_knots( &
  node_positions, &
  num_cells, &
  spline_degree, &
  bc_min, &
  bc_max, &
  knots_min, &
  knots_max ) &  
  result(knots)
  sll_real64, dimension(:), intent(in) :: node_positions
  sll_int32, intent(in) :: num_cells
  sll_int32, intent(in) :: spline_degree
  sll_int32, intent(in), optional :: bc_min
  sll_int32, intent(in), optional :: bc_max
  sll_real64, dimension(:), intent(in), optional :: knots_min
  sll_real64, dimension(:), intent(in), optional :: knots_max
  sll_real64, dimension(:), pointer :: knots

  !sll_real64 :: delta
  sll_int32 :: i
  sll_int32 :: ierr
  sll_real64 :: L
  sll_real64 :: val
  sll_int32 :: bc_min_value
  sll_int32 :: bc_max_value
  
  if(num_cells<1)then
    print *,'#num_cells should be>0',num_cells
    print *,'#at line/file',__LINE__,__FILE__
    stop
  endif
  if(spline_degree<1)then
    print *,'#spline_degree should be>0',spline_degree
  endif
  
  if(size(node_positions)<num_cells+1)then
    print *,'#Problem for size of node_positions',size(node_positions),num_cells+1
    print *,'#at line/file',__LINE__,__FILE__
    stop
  endif
  
  do i=1,num_cells-1
    val = node_positions(i+1)-node_positions(i)
    if(val<=0)then
      print *,'#node_positions should be strictly increasing'
      stop
    endif
  enddo
  
  
  
  L = node_positions(num_cells+1)-node_positions(1)
  
  
  SLL_ALLOCATE(knots(num_cells+2*spline_degree+1),ierr)
  
  knots(spline_degree+1:num_cells+spline_degree+1) = node_positions(1:num_cells+1)
  
  if(present(knots_min))then
    if(size(knots_min)<spline_degree)then
      print *,'#Problem of size of knots_min',size(knots_min),spline_degree
      print *,'#at line/file',__LINE__,__FILE__
      stop
    endif
    knots(1:spline_degree) = knots_min(1:spline_degree)
    if(present(bc_min))then
      print *,'#Warning knots_min is present'
      print *,'#bc_min=',bc_min
      print *,'is ignored at line/file',__LINE__,__FILE__
    endif
  else
    bc_min_value = sll_p_poisson_open_knots
    if(present(bc_min))then
      bc_min_value = bc_min
    endif
    select case (bc_min_value)
      case (sll_p_poisson_open_knots)
        knots(1:spline_degree) = knots(spline_degree+1)
      case (sll_p_poisson_periodic_knots)
        knots(1:spline_degree) = knots(num_cells+1:num_cells+spline_degree)-L
      case default
        print *, '#bad choice for bc_min'
        print *,'#error at line',__LINE__,'in file',__FILE__
        stop       
    end select  
  endif
  

  if(present(knots_max))then
    if(size(knots_max)<spline_degree)then
      print *,'#Problem of size of knots_max',size(knots_max),spline_degree
      print *,'#at line/file',__LINE__,__FILE__
      stop
    endif
    knots(num_cells+spline_degree+2:num_cells+2*spline_degree+1) &
      = knots_max(1:spline_degree)
    if(present(bc_max))then
      print *,'#Warning knots_max is present'
      print *,'#bc_max=',bc_max
      print *,'is ignored at line/file',__LINE__,__FILE__
    endif
  else
    bc_max_value = sll_p_poisson_open_knots
    if(present(bc_max))then
      bc_max_value = bc_max
    endif
    select case (bc_max_value)
      case (sll_p_poisson_open_knots)
        knots(num_cells+spline_degree+2:num_cells+2*spline_degree+1) &
          = knots(num_cells+spline_degree+1)
      case (sll_p_poisson_periodic_knots)
        knots(num_cells+spline_degree+2:num_cells+2*spline_degree+1) &
          = knots(spline_degree+2:2*spline_degree+1)+L
      case default
        print *, '#bad choice for bc_max'
        print *,'#error at line',__LINE__,'in file',__FILE__
        stop       
    end select  
  endif

end function new_knots





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
  local_spline_indices,&
  global_spline_indices,&
  local_to_global_spline_indices)
  
  sll_int32, intent(in) :: num_cells1
  sll_int32, intent(in) :: num_cells2
  sll_int32, intent(in) :: spline_degree1
  sll_int32, intent(in) :: spline_degree2
  sll_int32, dimension(:,:), intent(in) :: local_spline_indices
  sll_int32, dimension(:), intent(in) :: global_spline_indices
  sll_int32, dimension(:,:), intent(out) :: local_to_global_spline_indices
  
  sll_int32 :: i
  sll_int32 :: iloc
  
  do i= 1, num_cells1*num_cells2
    do iloc = 1, (spline_degree1+1)*(spline_degree2+1)
      local_to_global_spline_indices(iloc, i) = &
        global_spline_indices(local_spline_indices(iloc, i))
    enddo
  enddo
end subroutine compute_local_to_global_splines_indices


function compute_nz( &
  num_cells1, &
  num_cells2, &
  spline_degree1, &
  spline_degree2, &
  global_spline_indices &
  ) result(nz)
  sll_int32, intent(in) :: num_cells1
  sll_int32, intent(in) :: num_cells2
  sll_int32, intent(in) :: spline_degree1
  sll_int32, intent(in) :: spline_degree2
  sll_int32, dimension(:), intent(in) :: global_spline_indices
  sll_int32 :: nz

  sll_int32 :: i
  sll_int32 :: j
  sll_int32 :: i_loc
  sll_int32 :: j_loc
  sll_int32 :: i_glob
  sll_int32 :: j_glob
  sll_int32 :: num_spl1
  sll_int32 :: num_spl2
  !sll_int32 :: index1
  !sll_int32 :: index2
  sll_int32 :: li_A
  sll_int32 :: li_B
  !sll_int32, dimension(:,:,:), allocatable :: work
  sll_int32, dimension(:), allocatable :: flag
  sll_int32 :: ierr
  sll_int32 :: count
  
  num_spl1 = num_cells1+spline_degree1
  num_spl2 = num_cells2+spline_degree2

  !SLL_CLEAR_ALLOCATE(work(-degree_spline1:degree_spline1,-degree_spline2:degree_spline2,num_spl1*num_spl2),ierr)
  SLL_ALLOCATE(flag(0:num_spl1*num_spl2),ierr)
  flag=0

  count = 0
  do j_glob=1,num_spl2
    do i_glob=1,num_spl1
      li_A = global_spline_indices(i_glob+num_spl1*(j_glob-1))
      if(li_A/=0)then
      !if((flag(li_A)==0))then
        do j_loc = -spline_degree2,spline_degree2
          j = j_glob+j_loc
          if(j>0.and.j<=num_spl2)then
            do i_loc = -spline_degree1,spline_degree1
              i = i_glob+i_loc
              if(i>0.and.i<=num_spl1)then
                li_B = global_spline_indices(i+num_spl1*(j-1))
                if(li_B/=0)then
                  flag(li_A) = flag(li_A)+1
                  count = count+1
                endif
              endif          
            enddo
          endif  
        enddo
      !endif
      endif  
    enddo
  enddo  
  
  nz = count
  
  print *,maxval(flag),minval(flag)
  print *,'#val=',num_spl1*num_spl2*(2*spline_degree1+1)*(2*spline_degree2+1)

end function  compute_nz 
  

function compute_solution_size(bc_min,bc_max,num_cells, spline_degree) result(res)
  sll_int32, intent(in) :: bc_min
  sll_int32, intent(in) :: bc_max
  sll_int32, intent(in) :: num_cells
  sll_int32, intent(in) :: spline_degree
  sll_int32 :: res   
  
  if((bc_min==sll_p_periodic).and.(bc_max==sll_p_periodic))then
    res = num_cells
  else
    res = num_cells+spline_degree
  endif
  if(bc_min==sll_p_dirichlet)then
    res = res-1
  endif
  if(bc_max==sll_p_dirichlet)then
    res = res-1
  endif
  

end function compute_solution_size


!knots(cell+1-degree)<= ..<= knots(cell)<=x <knots(cell+1) <=..<=knots(cell+degree)
!!PN DEFINED BUT NOT USED
!subroutine compute_b_spline_at_x( &
!  knots, &
!  cell, &
!  x, &
!  degree, &
!  out)
!  sll_real64, dimension(:), intent(in) :: knots
!  sll_int32, intent(in) :: cell
!  sll_real64, intent(in) :: x
!  sll_int32, intent(in) :: degree
!  sll_real64, dimension(:), intent(out) :: out
!  
!  sll_real64 :: tmp1
!  sll_real64 :: tmp2
!  sll_int32 :: ell
!  sll_int32 :: k
!  
!  out(1) = 1._f64
!  do ell=1,degree
!    tmp1 = (x-knots(cell+1-degree))/(knots(cell+1)-knots(cell+1-degree))*out(1)
!    out(1) = out(1) -tmp1
!    do k=2,ell
!      tmp2 = (x-knots(cell+k-degree))/(knots(cell+k)-knots(cell+k-degree))*out(k)
!      out(k) = out(k)+tmp1-tmp2
!      tmp1 = tmp2
!    enddo
!    out(ell+1) = tmp1
!  enddo
!  
!  
!  
!  
!end subroutine compute_b_spline_at_x
!
!PN function check_compute_b_spline_at_x( &
!PN   knots, &
!PN   cell, &
!PN   x, &
!PN   degree &
!PN   ) result(res)
!PN   sll_real64, dimension(:), intent(in) :: knots
!PN   sll_int32, intent(in) :: cell
!PN   sll_real64, intent(in) :: x
!PN   sll_int32, intent(in) :: degree
!PN   sll_real64 :: res
!PN   
!PN   sll_real64, dimension(:), allocatable :: splines1
!PN   sll_real64, dimension(:), allocatable :: splines2
!PN   sll_int32 :: ierr
!PN   
!PN   res = 0._f64
!PN   SLL_ALLOCATE(splines1(degree+1),ierr)
!PN   SLL_ALLOCATE(splines2(degree+1),ierr)
!PN 
!PN   call bsplvb(knots,degree+1,cell,x,degree+1,splines2)
!PN   call compute_b_spline_at_x( &
!PN     knots, &
!PN     cell, &
!PN     x, &
!PN     degree, &
!PN     splines1)
!PN   res=maxval(abs(splines1-splines2))
!PN 
!PN end function check_compute_b_spline_at_x
!PN 
!PN 
!PN 
!PN 
!PN 
!PN subroutine compute_dbiatx( &
!PN   knots, &
!PN   knot_size, &
!PN   node_positions, &
!PN   num_points, &  
!PN   spline_degree, &
!PN   bc, &
!PN   dbiatx, &
!PN   dbiatx_indices, &
!PN   num_deriv)
!PN   sll_real64, dimension(:), intent(in) :: knots
!PN   sll_int32, intent(in) :: knot_size
!PN   sll_real64, dimension(:), intent(in) :: node_positions
!PN   sll_int32, intent(in) :: num_points
!PN   sll_int32, intent(in) :: spline_degree
!PN   sll_int32, intent(in) :: bc
!PN   sll_real64, dimension(:,:,:), intent(out) :: dbiatx
!PN   sll_int32, dimension(:), intent(out) :: dbiatx_indices
!PN   sll_int32, intent(in), optional :: num_deriv 
!PN   
!PN   sll_int32 :: i
!PN   sll_int32 :: iflag
!PN   !sll_real64 :: x
!PN   sll_int32 :: n_deriv
!PN   sll_real64, dimension(:,:), allocatable :: work
!PN   sll_real64, dimension(:,:), allocatable :: dbiatx_loc
!PN   sll_int32 :: ierr
!PN   
!PN   SLL_ALLOCATE(work(spline_degree+1,spline_degree+1),ierr)
!PN   
!PN   
!PN   n_deriv = 2
!PN   if(present(num_deriv)) then
!PN     n_deriv = num_deriv
!PN   endif
!PN 
!PN   SLL_ALLOCATE(dbiatx_loc(spline_degree+1,n_deriv),ierr)
!PN 
!PN   do i=1,num_points
!PN     call interv( &
!PN       knots, &
!PN       knot_size,  &
!PN       node_positions(i), &
!PN       dbiatx_indices(i), &
!PN       iflag)
!PN     call bsplvd( &
!PN       knots, &
!PN       spline_degree+1, &
!PN       node_positions(i), &
!PN       dbiatx_indices(i), &
!PN       work, &
!PN       dbiatx_loc, &
!PN       n_deriv)
!PN     dbiatx(1:spline_degree+1,1:n_deriv,i) = dbiatx_loc(1:spline_degree+1,1:n_deriv)   
!PN   enddo
!PN 
!PN   return
!PN   print*,bc
!PN   
!PN end subroutine compute_dbiatx
!PN 
!PN subroutine compute_composite_quadrature_points( &
!PN   eta_min, &
!PN   eta_max, &
!PN   num_cells, &
!PN   quadrature_points, &
!PN   num_quadrature_points, &
!PN   eta_loc_min, &
!PN   eta_loc_max, &
!PN   node_positions)
!PN   sll_real64, intent(in) :: eta_min
!PN   sll_real64, intent(in) :: eta_max
!PN   sll_int32, intent(in) :: num_cells
!PN   sll_real64, dimension(:), intent(in) :: quadrature_points
!PN   sll_int32, intent(in) :: num_quadrature_points
!PN   sll_real64, intent(in) :: eta_loc_min
!PN   sll_real64, intent(in) :: eta_loc_max
!PN   sll_real64, dimension(:), intent(out) :: node_positions
!PN   
!PN   sll_int32 :: i
!PN   sll_int32 :: j
!PN   sll_real64 :: eta
!PN   sll_real64 :: delta
!PN   sll_real64, dimension(:), allocatable :: loc_points
!PN   sll_int32 :: ierr
!PN   
!PN   SLL_ALLOCATE(loc_points(num_quadrature_points),ierr)
!PN   
!PN   delta = (eta_max-eta_min)/real(num_cells,f64)
!PN   loc_points(1:num_quadrature_points) = &
!PN     (quadrature_points(1:num_quadrature_points) - eta_loc_min) &
!PN     /(eta_loc_max-eta_loc_min)
!PN   
!PN   
!PN   do i=1,num_cells
!PN     eta = eta_min+real(i-1,f64)*delta
!PN     do j=1,num_quadrature_points
!PN       node_positions(j+num_quadrature_points*(i-1)) = eta+loc_points(i)*delta
!PN     enddo
!PN   enddo
!PN   
!PN   
!PN end subroutine compute_composite_quadrature_points
!PN 
!PN subroutine compute_global_knots( &
!PN   bc_left, &
!PN   bc_right, &
!PN   num_cells, &
!PN   eta_min, &
!PN   eta_max, &
!PN   spline_degree, &
!PN   knots)
!PN   sll_int32, intent(in) :: bc_left
!PN   sll_int32, intent(in) :: bc_right
!PN   sll_int32, intent(in) :: num_cells
!PN   sll_real64, intent(in) :: eta_min
!PN   sll_real64, intent(in) :: eta_max
!PN   sll_int32, intent(in) :: spline_degree
!PN   sll_real64, dimension(:), intent(out) :: knots
!PN   sll_real64 :: delta
!PN   sll_int32 :: i
!PN   
!PN   delta = (eta_max - eta_min)/real(num_cells,f64)
!PN   
!PN   
!PN   !repeated left point
!PN   do i = 1, spline_degree + 1
!PN     knots(i) = eta_min
!PN   enddo
!PN   !interior points
!PN   do i = spline_degree+2, num_cells+spline_degree
!PN     knots(i) = eta_min+real(i-spline_degree-1,f64)*delta
!PN   enddo
!PN   !repeated right point
!PN   do i = num_cells+spline_degree+1, num_cells+1+2*spline_degree
!PN     knots(i) = eta_max
!PN   enddo
!PN 
!PN 
!PN   if((bc_left==sll_p_periodic).and.(bc_right==sll_p_periodic))then
!PN     do i = 1, spline_degree+num_cells+1+2*spline_degree
!PN       knots(i) = eta_min+real(i-spline_degree-1,f64)*delta
!PN     enddo
!PN   endif
!PN 
!PN end subroutine compute_global_knots
!PN 
!PN subroutine compute_source_at_points( &
!PN   rhs_at_points, &
!PN   splines_at_points1, &
!PN   splines_at_points2, &
!PN   num_cells1, &
!PN   num_cells2, &
!PN   num_loc_points1, &
!PN   num_loc_points2, &
!PN   spline_degree1, &
!PN   spline_degree2, &
!PN   output)
!PN   
!PN   sll_real64, dimension(:,:), intent(in) :: rhs_at_points
!PN   sll_real64, dimension(:,:), intent(in) :: splines_at_points1
!PN   sll_real64, dimension(:,:), intent(in) :: splines_at_points2
!PN   sll_int32, intent(in) :: num_cells1
!PN   sll_int32, intent(in) :: num_cells2
!PN   sll_int32, intent(in) :: num_loc_points1
!PN   sll_int32, intent(in) :: num_loc_points2
!PN   sll_int32, intent(in) :: spline_degree1
!PN   sll_int32, intent(in) :: spline_degree2
!PN   sll_real64, dimension(:,:), intent(out) :: output
!PN   
!PN   sll_int32 :: num_spl1
!PN   sll_int32 :: num_spl2
!PN   sll_int32 :: i
!PN   sll_int32 :: j
!PN   sll_int32 :: iloc
!PN   sll_int32 :: jloc
!PN   sll_int32 :: ii
!PN   sll_int32 :: jj
!PN   
!PN   num_spl1 = num_cells1+spline_degree1+1
!PN   num_spl2 = num_cells2+spline_degree2+1
!PN   
!PN   output(1:num_spl1,1:num_spl2) = 0._f64
!PN   
!PN   do j=1,num_cells2
!PN     do i=1,num_cells1
!PN       do ii = 1,spline_degree1+1
!PN         do iloc = 1,num_loc_points1
!PN           do jj = 1,spline_degree2+1
!PN             do jloc = 1,num_loc_points2
!PN             output(i+ii,j+jj) = &
!PN               output(i+ii,j+jj) &
!PN               +splines_at_points1(ii,iloc+num_loc_points1*(i-1)) &
!PN               *rhs_at_points(iloc+num_loc_points1*(i-1),jloc+num_loc_points2*(j-1)) &
!PN               *splines_at_points2(jj,jloc+num_loc_points2*(j-1))
!PN             enddo
!PN           enddo    
!PN         enddo    
!PN       enddo  
!PN     enddo
!PN   enddo  
!PN 
!PN !  do j=1,num_cells2+spline_degree2+1
!PN !    do i=1,num_cells1+spline_degree1+1
!PN !      do ii = 1,spline_degree1+1 !iold+ii=i i-ii>=1 i-ii<=num_cells1
!PN !        if((ii<=1-i).and.(ii>=i-num_cells1))
!PN !        do iloc = 1,num_loc_points1
!PN !          output(i,j) = &
!PN !            output(i,j) &
!PN !            +spline_at_points1(ii,iloc+num_cells1*(i-ii-1)) &
!PN !            *rhs(iloc+num_cells1*(i-ii-1),jloc+num_cells2*(j-jj-1)) &
!PN !            *spline_at_points2(jj,jloc+num_cells2*(j-jj-1))
!PN !        enddo    
!PN !      enddo  
!PN !    enddo
!PN !  enddo  
!PN 
!PN   
!PN end subroutine compute_source_at_points


end module sll_m_poisson_2d_curvilinear
