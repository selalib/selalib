!> @ingroup poisson2d curvilinear solver
!> @brief Poisson solver on 2d curvilinear mesh
!> @details This solver works with geometry defined on gauss points

module sll_module_poisson_2d_curvilinear
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

use gauss_legendre_integration
use gauss_lobatto_integration
use sll_boundary_condition_descriptors
use sll_sparse_matrix_module

implicit none

type, public :: poisson_2d_curvilinear

  private
  sll_int32,  public :: num_cells1
  sll_int32,  public :: num_cells2
  sll_real64, public :: eta1_min
  sll_real64, public :: eta2_min   
  sll_real64, public :: eta1_max
  sll_real64, public :: eta2_max
  sll_real64, dimension(:,:), pointer :: gauss_pts1
  sll_real64, dimension(:,:), pointer :: gauss_pts2
  sll_int32 :: spline_degree1
  sll_int32 :: spline_degree2
  sll_int32, public :: knot_size1
  sll_int32, public :: knot_size2
  sll_real64, dimension(:), pointer :: knots1
  sll_real64, dimension(:), pointer :: knots2
  sll_int32, dimension(:,:), pointer :: local_spline_indices
  sll_int32, dimension(:), pointer :: bc_indices1
  sll_int32, dimension(:), pointer :: bc_indices2
  sll_int32, dimension(:), pointer, public :: global_spline_indices
  sll_int32, dimension(:,:), pointer :: local_to_global_spline_indices
  sll_int32 :: solution_size1
  sll_int32 :: solution_size2
  type(sll_csr_matrix), pointer :: sll_csr_mat
  sll_real64, dimension(:), pointer :: global_knots1
  sll_real64, dimension(:), pointer :: global_knots2




end type poisson_2d_curvilinear

!> For the integration mode.  
sll_int32, parameter, public :: POISSON_GAUSS_LEGENDRE = 0
!> For the integration mode.  
sll_int32, parameter, public :: POISSON_GAUSS_LOBATTO = 1


public new_poisson_2d_curvilinear

private

! *******************************************************************

contains 


!> @brief Initialization for poisson 2d curvilinear solver.
!> @details 
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
!> @return the type poisson_2d_curvilinear as a pointer
function new_poisson_2d_curvilinear( &
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
  eta2_max ) result(poisson)

  type(poisson_2d_curvilinear), pointer :: poisson
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


  SLL_ALLOCATE(poisson,ierr)
  call initialize_poisson_2d_curvilinear( &
    poisson, &
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
    eta2_max )
   
end function new_poisson_2d_curvilinear



! *******************************************************************

!> @brief Initialization of poisson 2d curvilinear solver
!> @details 
!> The parameters are
!> @param      poisson the type poisson_2d_curvilinear
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
subroutine initialize_poisson_2d_curvilinear( &
  poisson, &
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
    
  type(poisson_2d_curvilinear), intent(out) :: poisson
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
  sll_int32 :: dim1
  sll_int32 :: dim2
  sll_int32 :: num_spl1
  sll_int32 :: num_spl2
  sll_int32 :: solution_size
  sll_int32 :: total_num_splines_loc
  
  
  
  
  poisson%num_cells1 = num_cells1
  poisson%num_cells2 = num_cells2
  poisson%eta1_min = eta1_min
  poisson%eta1_max = eta1_max
  poisson%eta2_min = eta2_min
  poisson%eta2_max = eta2_max
  
  poisson%spline_degree1 = spline_degree1
  poisson%spline_degree2 = spline_degree2
  

  select case(quadrature_type1)
    case (POISSON_GAUSS_LEGENDRE)
      SLL_ALLOCATE(poisson%gauss_pts1(2,spline_degree1+2),ierr)
      poisson%gauss_pts1(:,:) = 0.0_f64
      poisson%gauss_pts1(:,:) = gauss_legendre_points_and_weights(spline_degree1+2)
    case (POISSON_GAUSS_LOBATTO)
      SLL_ALLOCATE(poisson%gauss_pts1(2,spline_degree1+2),ierr)
      poisson%gauss_pts1(:,:) = gauss_lobatto_points_and_weights(spline_degree1+2)
    case default
      print *, '#initialize_poisson_2d_curvilinear: '
      print *, '#bad choice for quadrature_type1'
      print *,'#error at line',__LINE__,'in file',__FILE__
      stop       
  end select

  select case(quadrature_type2)
    case (POISSON_GAUSS_LEGENDRE)
      SLL_ALLOCATE(poisson%gauss_pts2(2,spline_degree2+2),ierr)
      poisson%gauss_pts2(:,:) = 0.0_f64
      poisson%gauss_pts2(:,:) = gauss_legendre_points_and_weights(spline_degree2+2)
    case (POISSON_GAUSS_LOBATTO)
      SLL_ALLOCATE(poisson%gauss_pts1(2,spline_degree2+2),ierr)
      poisson%gauss_pts2(:,:) = gauss_lobatto_points_and_weights(spline_degree2+2)
    case default
      print *, '#initialize_poisson_2d_curvilinear: '
      print *, '#bad choice for quadrature_type2'
      print *,'#error at line',__LINE__,'in file',__FILE__
      stop       
  end select
  
  
  poisson%knot_size1 = compute_knot_size(bc_left, bc_right, num_cells1, spline_degree1)
  poisson%knot_size2 = compute_knot_size(bc_bottom, bc_top, num_cells2, spline_degree2)
  
  SLL_ALLOCATE(poisson%knots1(poisson%knot_size1),ierr)
  SLL_ALLOCATE(poisson%knots2(poisson%knot_size2),ierr)
  
  call compute_knots(bc_left, &
    bc_right, &
    num_cells1, &
    eta1_min, &
    eta1_max, &
    spline_degree1, &
    poisson%knots1)

  call compute_knots(bc_bottom, &
    bc_top, &
    num_cells2, &
    eta2_min, &
    eta2_max, &
    spline_degree2, &
    poisson%knots2)


  dim1 = (spline_degree1+1)*(spline_degree2+1)
  dim2 = num_cells1*num_cells2
  SLL_CLEAR_ALLOCATE(poisson%local_spline_indices(1:dim1,1:dim2),ierr)
  SLL_CLEAR_ALLOCATE(poisson%local_to_global_spline_indices(1:dim1,1:dim2),ierr)
  num_spl1=num_cells1+spline_degree1
  num_spl2=num_cells2+spline_degree2
  SLL_CLEAR_ALLOCATE(poisson%global_spline_indices(num_spl1*num_spl2),ierr)

  call compute_local_splines_indices( &
    num_cells1, &
    num_cells2, &
    spline_degree1, &
    spline_degree2, &
    poisson%local_spline_indices)

  SLL_ALLOCATE(poisson%bc_indices1(num_spl1),ierr)
  SLL_ALLOCATE(poisson%bc_indices2(num_spl2),ierr)

  call compute_bc_indices(num_cells1,spline_degree1,bc_left,bc_right,poisson%bc_indices1)
  call compute_bc_indices(num_cells2,spline_degree2,bc_bottom,bc_top,poisson%bc_indices2)
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
    bc_left, &
    bc_right, &
    num_cells1, &
    spline_degree1)
  poisson%solution_size2 = compute_solution_size( &
    bc_bottom, &
    bc_top, &
    num_cells2, &
    spline_degree2)
  
  solution_size = poisson%solution_size1*poisson%solution_size2
  total_num_splines_loc = (spline_degree1+1)*(spline_degree2+1)

  poisson%sll_csr_mat => new_csr_matrix( &
    solution_size, &
    solution_size, &
    num_cells1*num_cells2, &
    poisson%local_to_global_spline_indices, &
    total_num_splines_loc, &
    poisson%local_to_global_spline_indices, &
    total_num_splines_loc)




end subroutine initialize_poisson_2d_curvilinear

function compute_knot_size(bc_left,bc_right,num_cells, spline_degree) result(res)
  sll_int32, intent(in) :: bc_left
  sll_int32, intent(in) :: bc_right
  sll_int32, intent(in) :: num_cells
  sll_int32, intent(in) :: spline_degree
  sll_int32 :: res   
  
  if((bc_left==SLL_PERIODIC).and.(bc_right==SLL_PERIODIC))then
    res = 2*spline_degree+2
  else
    res = 2*spline_degree+num_cells+1
  endif

end function compute_knot_size

subroutine compute_knots( &
  bc_left, &
  bc_right, &
  num_cells, &
  eta_min, &
  eta_max, &
  spline_degree, &
  knots)
  sll_int32, intent(in) :: bc_left
  sll_int32, intent(in) :: bc_right
  sll_int32, intent(in) :: num_cells
  sll_real64, intent(in) :: eta_min
  sll_real64, intent(in) :: eta_max
  sll_int32, intent(in) :: spline_degree
  sll_real64, dimension(:), intent(out) :: knots
  sll_real64 :: delta
  sll_int32 :: i
  
  delta = (eta_max - eta_min)/real(num_cells,f64)
  
  
  if((bc_left==SLL_PERIODIC).and.(bc_right==SLL_PERIODIC))then
    do i = -spline_degree, spline_degree+1
      knots(i+spline_degree+1) = real(i,f64)*delta 
    end do    
  else
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
  endif

end subroutine compute_knots


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

subroutine compute_bc_indices(num_cells,spline_degree,bc_left,bc_right,index)
  sll_int32, intent(in) :: num_cells
  sll_int32, intent(in) :: spline_degree
  sll_int32, intent(in) :: bc_left
  sll_int32, intent(in) :: bc_right
  sll_int32, dimension(:), intent(out) :: index
  
  sll_int32 :: i
    
  do i = 1, num_cells+spline_degree
    index(i) = i
  enddo
  
  !left boundary correction  
  select case (bc_left)
    case (SLL_DIRICHLET)
      index(1) = 0
    case default     
  end select
      
  !right boundary correction
  select case (bc_right)
    case (SLL_DIRICHLET)
      index(num_cells+spline_degree) = 0
    case (SLL_PERIODIC)
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


function compute_solution_size(bc_left,bc_right,num_cells, spline_degree) result(res)
  sll_int32, intent(in) :: bc_left
  sll_int32, intent(in) :: bc_right
  sll_int32, intent(in) :: num_cells
  sll_int32, intent(in) :: spline_degree
  sll_int32 :: res   
  
  if((bc_left==SLL_PERIODIC).and.(bc_right==SLL_PERIODIC))then
    res = num_cells
  else
    res = num_cells+spline_degree
  endif
  if(bc_left==SLL_DIRICHLET)then
    res = res-1
  endif
  if(bc_right==SLL_DIRICHLET)then
    res = res-1
  endif
  

end function compute_solution_size


subroutine compute_dbiatx( &
  knots, &
  knot_size, &
  node_positions, &
  num_points, &  
  spline_degree, &
  bc, &
  dbiatx, &
  dbiatx_indices, &
  num_deriv)
  sll_real64, dimension(:), intent(in) :: knots
  sll_int32, intent(in) :: knot_size
  sll_real64, dimension(:), intent(in) :: node_positions
  sll_int32, intent(in) :: num_points
  sll_int32, intent(in) :: spline_degree
  sll_int32, intent(in) :: bc
  sll_real64, dimension(:,:,:), intent(out) :: dbiatx
  sll_int32, dimension(:), intent(out) :: dbiatx_indices
  sll_int32, intent(in), optional :: num_deriv 
  
  sll_int32 :: i
  sll_int32 :: iflag
  sll_real64 :: x
  sll_int32 :: n_deriv
  sll_real64, dimension(:,:), allocatable :: work
  sll_real64, dimension(:,:), allocatable :: dbiatx_loc
  sll_int32 :: ierr
  
  SLL_ALLOCATE(work(spline_degree+1,spline_degree+1),ierr)
  
  
  n_deriv = 2
  if(present(num_deriv)) then
    n_deriv = num_deriv
  endif

  SLL_ALLOCATE(dbiatx_loc(spline_degree+1,n_deriv),ierr)

  do i=1,num_points
    call interv( &
      knots, &
      knot_size,  &
      node_positions(i), &
      dbiatx_indices(i), &
      iflag)
    call bsplvd( &
      knots, &
      spline_degree+1, &
      node_positions(i), &
      dbiatx_indices(i), &
      work, &
      dbiatx_loc, &
      n_deriv)
    dbiatx(1:spline_degree+1,1:n_deriv,i) = dbiatx_loc(1:spline_degree+1,1:n_deriv)   
  enddo
  
end subroutine compute_dbiatx

subroutine compute_composite_quadrature_points( &
  eta_min, &
  eta_max, &
  num_cells, &
  quadrature_points, &
  num_quadrature_points, &
  eta_loc_min, &
  eta_loc_max, &
  node_positions)
  sll_real64, intent(in) :: eta_min
  sll_real64, intent(in) :: eta_max
  sll_int32, intent(in) :: num_cells
  sll_real64, dimension(:), intent(in) :: quadrature_points
  sll_int32, intent(in) :: num_quadrature_points
  sll_real64, intent(in) :: eta_loc_min
  sll_real64, intent(in) :: eta_loc_max
  sll_real64, dimension(:), intent(out) :: node_positions
  
  sll_int32 :: i
  sll_int32 :: j
  sll_real64 :: eta
  sll_real64 :: delta
  sll_real64, dimension(:), allocatable :: loc_points
  sll_int32 :: ierr
  
  SLL_ALLOCATE(loc_points(num_quadrature_points),ierr)
  
  delta = (eta_max-eta_min)/real(num_cells,f64)
  loc_points(1:num_quadrature_points) = &
    (quadrature_points(1:num_quadrature_points) - eta_loc_min) &
    /(eta_loc_max-eta_loc_min)
  
  
  do i=1,num_cells
    eta = eta_min+real(i-1,f64)*delta
    do j=1,num_quadrature_points
      node_positions(j+num_quadrature_points*(i-1)) = eta+loc_points(i)*delta
    enddo
  enddo
  
  
end subroutine compute_composite_quadrature_points

subroutine compute_global_knots( &
  bc_left, &
  bc_right, &
  num_cells, &
  eta_min, &
  eta_max, &
  spline_degree, &
  knots)
  sll_int32, intent(in) :: bc_left
  sll_int32, intent(in) :: bc_right
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


  if((bc_left==SLL_PERIODIC).and.(bc_right==SLL_PERIODIC))then
    do i = 1, spline_degree+num_cells+1+2*spline_degree
      knots(i) = eta_min+real(i-spline_degree-1,f64)*delta
    enddo
  endif

end subroutine compute_global_knots



end module sll_module_poisson_2d_curvilinear
