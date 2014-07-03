module sll_general_coordinate_elliptic_solver_multipatch_module
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

  use sll_boundary_condition_descriptors
  use sll_module_scalar_field_2d_base
  use sll_module_scalar_field_2d_alternative
  use sll_arbitrary_degree_spline_interpolator_2d_module
  use sparsematrix_module
  use connectivity_module
  use sll_knots
  use gauss_legendre_integration
  use gauss_lobatto_integration
  use sll_timer 
  use sll_sparse_matrix_module
  use sll_sparse_matrix_mp_module
  use sll_module_scalar_field_2d_multipatch
  use sll_general_coordinate_elliptic_solver_module

  implicit none

  type :: general_coordinate_elliptic_solver_mp
     sll_int32,dimension(:),pointer :: total_num_splines_loc
     sll_int32 :: total_num_splines_eta1
     sll_int32 :: total_num_splines_eta2
     type(sll_coordinate_transformation_multipatch_2d), pointer :: T
     sll_real64, dimension(:,:), pointer :: knots1
     sll_real64, dimension(:,:), pointer :: knots2
     sll_real64, dimension(:,:), pointer :: knots1_rho
     sll_real64, dimension(:,:), pointer :: knots2_rho
     sll_real64, dimension(:,:,:), pointer :: gauss_pts1
     sll_real64, dimension(:,:,:), pointer :: gauss_pts2!
     sll_int32,dimension(:),pointer :: spline_degree1
     sll_int32,dimension(:),pointer :: spline_degree2
     sll_int32 :: solution_size_final,solution_size
     sll_int32, dimension(:),pointer :: num_elts_no_null_patchs
     sll_int32, dimension(:,:,:),pointer :: local_to_global_row, local_to_global_col
     sll_int32, dimension(:,:,:),pointer :: global_indices_row, global_indices_col
     sll_real64 :: intjac
     sll_int32, dimension(:,:,:), pointer :: local_to_global_spline_indices_source_col
     sll_int32, dimension(:,:,:), pointer :: local_to_global_spline_indices_source_row
     ! contains the values of all splines in all gauss points
     sll_real64, dimension(:,:,:), pointer :: values_jacobian
     sll_int32 , dimension(:,:)  , pointer :: tab_index_coeff1
     sll_int32 , dimension(:,:)  , pointer :: tab_index_coeff2
     type(sll_csr_matrix), pointer :: sll_csr_mat
     type(sll_csr_matrix), pointer :: sll_csr_mat_source
     sll_real64, dimension(:), pointer :: rho_vec
     sll_real64, dimension(:), pointer :: phi_vec
  end type general_coordinate_elliptic_solver_mp
 
  
  interface delete
     module procedure delete_elliptic_mp
  end interface delete
  
  interface initialize
     module procedure initialize_general_elliptic_solver_mp
  end interface initialize
  
  
contains ! *******************************************************************
  
  
  
  !> @brief Initialization for elleptic solver multipatch
  !> @details To have the function phi such that 
  !>  div( A grad phi ) + B grad phi + C phi = rho
  !>  where A is a matrix of functions , B a vectorial function,
  !>  and  C and rho a scalar function.  
  !>  A, B, C, rho can be discret or analytic. 
  !>  phi is given with a B-spline interpolator  
  !> 
  !> The parameters are
  !> @param es the type general_coordinate_elliptic_solver_mp
  !> @param[in] spline_degree_eta1 the degre of B-spline in the direction eta1
  !> @param[in] spline_degree_eta2 the degre of B-spline in the direction eta2
  !> @param[in] quadrature_type1 the type of quadrature in the direction eta1
  !> @param[in] quadrature_type2 the type of quadrature in the direction eta2
  !> @param[in] T the transformation multipatch
  !> @return the type general_coordinate_elliptic_solver_mp
  
  subroutine initialize_general_elliptic_solver_mp( &
       es_mp, &
       quadrature_type1, &
       quadrature_type2, &
       T)
    
    type(general_coordinate_elliptic_solver_mp), intent(out) :: es_mp
    type(sll_coordinate_transformation_multipatch_2d), pointer :: T
    type(sll_logical_mesh_2d), pointer                         :: lm
    sll_int32, intent(in) :: quadrature_type1
    sll_int32, intent(in) :: quadrature_type2
    sll_int32 :: knots1_size
    sll_int32 :: knots2_size
    sll_int32 :: max_degspline1
    sll_int32 :: max_degspline2
    sll_int32 :: ierr,ierr1
    sll_int32 :: i,j,k,l
    sll_int32 :: max_num_cells_eta1, max_num_cells_eta2
    sll_int32 :: max_num_elts_no_null_patchs
    sll_int32 :: max_total_num_splines_loc
    sll_int32 :: num_cells_eta1,num_cells_eta2
    sll_int32 :: num_patches

    es_mp%T => T

    num_patches = es_mp%T%number_patches
    SLL_ALLOCATE(es_mp%total_num_splines_loc(num_patches),ierr)
    SLL_ALLOCATE(es_mp%spline_degree1(num_patches),ierr)
    SLL_ALLOCATE(es_mp%spline_degree2(num_patches),ierr)

    do i = 1, es_mp%T%number_patches
       es_mp%spline_degree1(i) = es_mp%T%transfs(i)%T%spline_deg1
       es_mp%spline_degree2(i) = es_mp%T%transfs(i)%T%spline_deg2
       es_mp%total_num_splines_loc(i) = (es_mp%spline_degree1(i)+1)*(es_mp%spline_degree2(i)+1)
    end do
  
   
    ! ------------------------------------------------------------------
    ! Allocate and fill the gauss points/weights information.
    ! First direction
    max_degspline1 = MAXVAL(es_mp%spline_degree1(:))
    max_degspline2 = MAXVAL(es_mp%spline_degree2(:))
    SLL_ALLOCATE(es_mp%gauss_pts1(num_patches,2,max_degspline1+2),ierr)
    SLL_ALLOCATE(es_mp%gauss_pts2(num_patches,2,max_degspline2+2),ierr)

    select case(quadrature_type1)
    case (ES_GAUSS_LEGENDRE)
       
       
       do i = 1, es_mp%T%number_patches
          es_mp%gauss_pts1(i,:,1:es_mp%spline_degree1(i)+2) = gauss_legendre_points_and_weights(es_mp%spline_degree1(i)+2)
       end do

    case (ES_GAUSS_LOBATTO)

       do i = 1, es_mp%T%number_patches
          es_mp%gauss_pts1(i,:,1:es_mp%spline_degree1(i)+2) = gauss_lobatto_points_and_weights(es_mp%spline_degree1(i)+2)
       end do
       
    case DEFAULT
       print *, 'new_general_qn_solver(): have not type of gauss points in the first direction'
    end select
    
    select case(quadrature_type2)
    case (ES_GAUSS_LEGENDRE)
       
       do i = 1, es_mp%T%number_patches
          es_mp%gauss_pts2(i,:,1:es_mp%spline_degree2(i)+2) = gauss_legendre_points_and_weights(es_mp%spline_degree2(i)+2)
       end do
       
    case (ES_GAUSS_LOBATTO)
       
       do i = 1, es_mp%T%number_patches
          es_mp%gauss_pts2(i,:,1:es_mp%spline_degree2(i)+2) = gauss_lobatto_points_and_weights(es_mp%spline_degree2(i)+2)
       end do

    case DEFAULT
       print *, 'new_general_qn_solver(): have not type of gauss points in the second direction'
       
    end select
   
    ! ------------------------------------------------------------------


    lm => es_mp%T%get_logical_mesh(0)

    ! ------------------------------------------------------------------
    ! coeff's size for the solution phi in all the domain
    es_mp%solution_size_final = MAXVAL(es_mp%T%global_indices)
    es_mp%solution_size = size(es_mp%T%global_indices)
    SLL_ALLOCATE(es_mp%rho_vec(es_mp%solution_size_final),ierr)
    SLL_ALLOCATE(es_mp%phi_vec(es_mp%solution_size_final),ierr)
    ! ------------------------------------------------------------------

    
    ! ------------------------------------------------------------------
    ! size of knots and allocation 
    knots1_size = size(es_mp%T%transfs(1)%T%knots1)
    knots2_size = size(es_mp%T%transfs(1)%T%knots2)
    if (es_mp%T%number_patches >1) then
       do i = 2, es_mp%T%number_patches
          if (knots1_size > size(es_mp%T%transfs(i-1)%T%knots1)) then
             knots1_size = size(es_mp%T%transfs(i-1)%T%knots1)
          end if
          if (knots2_size > size(es_mp%T%transfs(i-1)%T%knots2)) then
             knots2_size = size(es_mp%T%transfs(i-1)%T%knots2)
          end if
       end do
    end if

    SLL_ALLOCATE(es_mp%knots1(num_patches,knots1_size),ierr1)
    SLL_ALLOCATE(es_mp%knots2(num_patches,knots2_size),ierr)

    ! ------------------------------------------------------------------


   es_mp%rho_vec(:) = 0.0_f64
   es_mp%phi_vec(:) = 0.0_f64
   es_mp%intjac     = 0.0_f64
   ! Knots to define be careful

   do i = 1, num_patches

      es_mp%knots1(i,:) = es_mp%T%transfs(i)%T%knots1
      es_mp%knots2(i,:) = es_mp%T%transfs(i)%T%knots2
   end do
   
   
   SLL_ALLOCATE( es_mp%num_elts_no_null_patchs(num_patches),ierr)
   
   do i=1,num_patches
      
      lm => es_mp%T%get_logical_mesh(i-1)
      num_cells_eta1 = lm%num_cells1
      num_cells_eta2 = lm%num_cells2
      es_mp%num_elts_no_null_patchs(i) = num_cells_eta1*num_cells_eta2

   end do

   max_num_elts_no_null_patchs = MAXVAL(es_mp%num_elts_no_null_patchs)
   
   max_total_num_splines_loc   = MAXVAL(es_mp%total_num_splines_loc)
   SLL_ALLOCATE(es_mp%local_to_global_row(num_patches,max_total_num_splines_loc,max_num_elts_no_null_patchs),ierr)
   SLL_ALLOCATE(es_mp%local_to_global_col(num_patches,max_total_num_splines_loc,max_num_elts_no_null_patchs),ierr)
   SLL_ALLOCATE(es_mp%local_to_global_spline_indices_source_row(num_patches,max_total_num_splines_loc,max_num_elts_no_null_patchs),ierr)
   SLL_ALLOCATE(es_mp%local_to_global_spline_indices_source_col(num_patches,max_total_num_splines_loc,max_num_elts_no_null_patchs),ierr)
   !SLL_ALLOCATE(global_indices_row(),ierr)
   !SLL_ALLOCATE(global_indices_col(),ierr)
   
   es_mp%local_to_global_spline_indices_source_row = 0
   es_mp%local_to_global_spline_indices_source_col = 0

   do i = 1, num_patches
      lm => es_mp%T%get_logical_mesh(i-1)
      do j =1,es_mp%total_num_splines_loc(i)
         do k=1,lm%num_cells1
            do l=1,lm%num_cells2
              
               es_mp%local_to_global_row(i,j,k + (l-1) *lm%num_cells1 ) =&
                    T%get_spline_local_to_global_index(i-1,j,k,l)
               es_mp%local_to_global_col(i,j,k + (l-1) *lm%num_cells1 ) = &
                    es_mp%local_to_global_row(i,j,k + (l-1) *lm%num_cells1 ) 
            end do
         end do
      enddo
   end do
         
   es_mp%intjac = 0.0_f64

   es_mp%sll_csr_mat => new_csr_matrix_mp( &
        es_mp%solution_size_final, &
        es_mp%solution_size_final, &
        es_mp%T%number_patches,&
        es_mp%num_elts_no_null_patchs, &
        es_mp%local_to_global_row, &
        es_mp%total_num_splines_loc, &
        es_mp%local_to_global_col, &
        es_mp%total_num_splines_loc ) 


!!$   do i = 1, num_patches
!!$      lm => es_mp%T%get_logical_mesh(i-1)
!!$      es_mp%knots1_rho (i, 1 : es_mp%spline_degree1(i) +1 ) = lm%eta1_min
!!$      es_mp%knots1_rho (i,lm%num_cells1+2:lm%num_cells1+2+es_mp%spline_degree1(i) ) = lm%eta1_max
!!$   
!!$   
!!$      if ( mod(es_mp%spline_degree1(i) +1,2) == 0 ) then
!!$         do j = es_mp%spline_degree1(i) +1 + 1, lm%num_cells1 + 1
!!$            es_mp%knots1_rho ( i,j ) = lm%eta1_min +  ( j - (es_mp%spline_degree1(i) +1)/2-1 )*lm%delta_eta1 
!!$         
!!$         end do
!!$      else
!!$      
!!$         do j = es_mp%spline_degree1(i) +1 + 1, lm%num_cells1 + 1
!!$            es_mp%knots1_rho ( i,j ) = &
!!$                0.5*( lm%eta1_min + ( j - (es_mp%spline_degree1(i))/2 -1)*lm%delta_eta1 + &
!!$                lm%eta1_min +  ( j -1 - (es_mp%spline_degree1(i))/2 -1)*lm%delta_eta1 )
!!$           
!!$         end do
!!$        
!!$      end if
!!$
!!$      es_mp%knots2_rho (i, 1 : es_mp%spline_degree2(i) +1 ) = lm%eta2_min
!!$      es_mp%knots2_rho (i,lm%num_cells2 + 2 : lm%num_cells2 + 2 + es_mp%spline_degree2(i) ) = lm%eta2_max
!!$     
!!$     
!!$      if ( mod(es_mp%spline_degree2(i) +1,2) == 0 ) then
!!$         do j = es_mp%spline_degree2(i) +1 + 1, lm%num_cells2 + 1
!!$            es_mp%knots2_rho ( i,j ) = lm%eta2_min +  ( j - (es_mp%spline_degree2(i) +1)/2-1 )*lm%delta_eta2 
!!$           
!!$         end do
!!$      else
!!$        
!!$        do j = es_mp%spline_degree2(i) +1 + 1, lm%num_cells2 + 1
!!$           es_mp%knots2_rho ( i,j ) = &
!!$                0.5*( lm%eta2_min + ( j - (es_mp%spline_degree2(i))/2 -1)*lm%delta_eta2 + &
!!$                lm%eta2_min +  ( j -1 - (es_mp%spline_degree2(i))/2 -1)*lm%delta_eta2 )
!!$           
!!$        end do
!!$        
!!$     end if
!!$
!!$  enddo
     
     
     
   !! allocation of the table containning all values of splines in each direction 
   !! in each gauss points
   
   lm=>es_mp%T%get_logical_mesh(0)

   max_num_cells_eta1 = lm%num_cells1
   max_num_cells_eta2 = lm%num_cells2

   SLL_ALLOCATE(es_mp%values_jacobian(num_patches,max_num_cells_eta1*(max_degspline1+2),max_num_cells_eta2*(max_degspline2+2)),ierr)
   es_mp%values_jacobian = 0.0_f64
   
   SLL_ALLOCATE(es_mp%tab_index_coeff1(num_patches,max_num_cells_eta1*(max_degspline1+2)),ierr)
   SLL_ALLOCATE(es_mp%tab_index_coeff2(num_patches,max_num_cells_eta2*(max_degspline2+2)),ierr)
   
 end subroutine initialize_general_elliptic_solver_mp
 

 
 !> @brief Initialization for elleptic solver.
 !> @details To have the function phi such that 
 !>  div( A grad phi ) + B grad phi + C phi = rho
 !>  where A is a matrix of functions , B a vectorial function,
 !>  and  C and rho a scalar function.  
 !>  A, B, C, rho can be discret or analytic. 
 !>  phi is given with a B-spline interpolator  
 !> 
 !> The parameters are
 !> @param[in] spline_degree_eta1 the degre of B-spline in the direction eta1
 !> @param[in] spline_degree_eta2 the degre of B-spline in the direction eta2
 !> @param[in] quadrature_type1 the type of quadrature in the direction eta1
 !> @param[in] quadrature_type2 the type of quadrature in the direction eta2
 !> @param[in] T the transformation multipatch
 !> @return the type general_coordinate_elliptic_solver such that a pointer
 function new_general_elliptic_solver_mp( &
      quadrature_type1, &
      quadrature_type2, &
      T) result(es_mp)
   
   type(general_coordinate_elliptic_solver_mp), pointer :: es_mp
   type(sll_coordinate_transformation_multipatch_2d), pointer :: T
   sll_int32, intent(in) :: quadrature_type1
   sll_int32, intent(in) :: quadrature_type2
   sll_int32 :: ierr
   type(sll_time_mark)  :: t0 
   double precision :: time
   
   
   call sll_set_time_mark(t0)
   SLL_ALLOCATE(es_mp,ierr)
   call initialize( &
        es_mp, &
        quadrature_type1, &
        quadrature_type2,&
        T)
   
   time = sll_time_elapsed_since(t0)
   print*, '#time for new_general_elliptic_solver', time
   
   
 end function new_general_elliptic_solver_mp
 
 
 !> @brief Deallocate the type general_coordinate_elliptic_solver
 !> 
 !> 
 !> The parameters are
 !> @param[in] es the type general_coordinate_elliptic_solver
 
 subroutine delete_elliptic_mp( es_mp )
   type(general_coordinate_elliptic_solver_mp) :: es_mp
   sll_int32 :: ierr
   ! it is not good to check some cases and not others, fix...
   if(associated(es_mp%knots1)) then
      SLL_DEALLOCATE(es_mp%knots1,ierr)
   else
      print *, 'delete es, WARNING: knots1 array was not allocated.'
   end if
   if(associated(es_mp%knots2)) then
      SLL_DEALLOCATE(es_mp%knots2,ierr)
   else
      print *, 'delete es general coords, ', &
            'WARNING: knots2 array was not allocated.'
   end if
   SLL_DEALLOCATE(es_mp%gauss_pts1,ierr)
   SLL_DEALLOCATE(es_mp%gauss_pts2,ierr)
   SLL_DEALLOCATE(es_mp%num_elts_no_null_patchs,ierr)
   SLL_DEALLOCATE(es_mp%local_to_global_row,ierr)
   SLL_DEALLOCATE(es_mp%local_to_global_col,ierr)
   !SLL_DEALLOCATE(es_mp%global_spline_indices,ierr)
   !SLL_DEALLOCATE(es_mp%local_spline_indices,ierr)
   !SLL_DEALLOCATE(es_mp%local_to_global_spline_indices,ierr)
    SLL_DEALLOCATE(es_mp%local_to_global_spline_indices_source_row,ierr)
    SLL_DEALLOCATE(es_mp%local_to_global_spline_indices_source_col,ierr)
   call sll_delete(es_mp%sll_csr_mat)
   call sll_delete(es_mp%sll_csr_mat_source)
   SLL_DEALLOCATE(es_mp%rho_vec,ierr)
   SLL_DEALLOCATE(es_mp%phi_vec,ierr)
   SLL_DEALLOCATE(es_mp%knots1_rho,ierr)
   SLL_DEALLOCATE(es_mp%knots2_rho,ierr)
   SLL_DEALLOCATE(es_mp%values_jacobian,ierr)
   SLL_DEALLOCATE(es_mp%tab_index_coeff1,ierr)
   SLL_DEALLOCATE(es_mp%tab_index_coeff2,ierr)
   SLL_DEALLOCATE(es_mp%total_num_splines_loc,ierr)
   SLL_DEALLOCATE(es_mp%spline_degree1,ierr)
   SLL_DEALLOCATE(es_mp%spline_degree2,ierr)
   !call delete(es_mp%logical_mesh)
 end subroutine delete_elliptic_mp
 
  
 !> @brief Assemble the matrix for elliptic solver.
 !> @details To have the function phi such that 
 !>  div( A grad phi ) + B grad phi + C phi = rho
 !>  where A is a matrix of functions , B a vectorial function,
 !>  and  C and rho a scalar function.  
 !>  A, B, C, rho can be discret or analytic. 
 !>  phi is given with a B-spline interpolator  
 !> 
 !> The parameters are
 !> @param es the type general_coordinate_elliptic_solver
 !> @param[in] a11_field_mat the field corresponding to the coefficient A(1,1) of the matrix A  
 !> @param[in] a12_field_mat the field corresponding to the coefficient A(1,2) of the matrix A 
 !> @param[in] a21_field_mat the field corresponding to the coefficient A(2,1) of the matrix A  
 !> @param[in] a22_field_mat the field corresponding to the coefficient A(2,2) of the matrix A 
 !> @param[in] b1_field_vect the field corresponding to the coefficient B(1) of the vector B  
 !> @param[in] b2_field_vect the field corresponding to the coefficient B(2) of the vector B
 !> @param[in] c_field the field corresponding to the coefficient B(1) of the scalar C
 !> @return the type general_coordinate_elliptic_solver contains the matrix to solve the equation
 subroutine factorize_mat_es_mp(&
      es_mp, &
      a11_field_mat, &
      a12_field_mat,&
      a21_field_mat,&
      a22_field_mat,&
      b1_field_vect,&
      b2_field_vect,&
      c_field)!, &
   ! rho)
   use sll_timer
   type(general_coordinate_elliptic_solver_mp),intent(inout) :: es_mp
   class(sll_scalar_field_multipatch_2d), pointer :: a11_field_mat
   class(sll_scalar_field_multipatch_2d), pointer :: a12_field_mat
   class(sll_scalar_field_multipatch_2d), pointer :: a21_field_mat
   class(sll_scalar_field_multipatch_2d), pointer :: a22_field_mat
   class(sll_scalar_field_multipatch_2d), pointer :: b1_field_vect
   class(sll_scalar_field_multipatch_2d), pointer :: b2_field_vect
   class(sll_scalar_field_multipatch_2d), pointer     :: c_field
   sll_real64, dimension(:,:), allocatable :: M_c_loc
   sll_real64, dimension(:,:), allocatable :: K_a11_loc
   sll_real64, dimension(:,:), allocatable :: K_a12_loc
   sll_real64, dimension(:,:), allocatable :: K_a21_loc
   sll_real64, dimension(:,:), allocatable :: K_a22_loc
   sll_real64, dimension(:,:), allocatable :: M_b_vect_loc
   sll_real64, dimension(:,:), allocatable :: S_b1_loc
   sll_real64, dimension(:,:), allocatable :: S_b2_loc 
   sll_real64, dimension(:,:,:,:), pointer :: Source_loc
   sll_int32 :: total_num_splines_loc
   sll_int32 :: ierr
   sll_int32 :: j,i
   sll_int32 :: cell_index
   sll_int32 :: patch
   type(sll_time_mark)  :: t0 
   double precision :: time
   type(sll_logical_mesh_2d), pointer:: lm
   character(len=*),parameter :: as_file1='mat'
   
   call sll_set_time_mark(t0)
   
   total_num_splines_loc = maxval(es_mp%total_num_splines_loc)
   
   SLL_ALLOCATE(M_c_loc(total_num_splines_loc,total_num_splines_loc),ierr)
   SLL_ALLOCATE(K_a11_loc(total_num_splines_loc,total_num_splines_loc),ierr)
   SLL_ALLOCATE(K_a12_loc(total_num_splines_loc,total_num_splines_loc),ierr)
   SLL_ALLOCATE(K_a21_loc(total_num_splines_loc,total_num_splines_loc),ierr)
   SLL_ALLOCATE(K_a22_loc(total_num_splines_loc,total_num_splines_loc),ierr)
   SLL_ALLOCATE(S_b1_loc(total_num_splines_loc,total_num_splines_loc),ierr)
   SLL_ALLOCATE(S_b2_loc(total_num_splines_loc,total_num_splines_loc),ierr)
   SLL_ALLOCATE(M_b_vect_loc(total_num_splines_loc,total_num_splines_loc),ierr)

   lm => es_mp%T%get_logical_mesh(0)
   SLL_ALLOCATE(Source_loc(es_mp%T%number_patches,lm%num_cells1*lm%num_cells2,total_num_splines_loc,total_num_splines_loc),ierr)
   Source_loc(:,:,:,:) = 0.0_f64

   
   
   do patch = 1,es_mp%T%number_patches

      lm => es_mp%T%get_logical_mesh(patch-1)
      
      do j=1,lm%num_cells2
         do i=1,lm%num_cells1
            
            M_c_loc   = 0.0_f64
            K_a11_loc = 0.0_f64
            K_a12_loc = 0.0_f64
            K_a21_loc = 0.0_f64 
            K_a22_loc = 0.0_f64
            S_b1_loc  = 0.0_f64
            S_b2_loc  = 0.0_f64

            cell_index = i+lm%num_cells1*(j-1)
            
            call build_local_matrices_mp( &
                 es_mp, &
                 cell_index,&
                 i, &
                 j, &
                 patch,&
                 a11_field_mat, &
                 a12_field_mat, &
                 a21_field_mat, &
                 a22_field_mat, &
                 b1_field_vect,&
                 b2_field_vect,&
                 c_field, &
                 M_c_loc, &
                 K_a11_loc, &
                 K_a12_loc, &
                 K_a21_loc, &
                 K_a22_loc, &
                 M_b_vect_loc, &
                 S_b1_loc,  &
                 S_b2_loc,&
                 Source_loc)
 
            call local_to_global_matrices_mp( &
                 es_mp, &
                 cell_index, &
                 i, &
                 j, &
                 patch,&
                 M_c_loc, &
                 K_a11_loc, &
                 K_a12_loc, &
                 K_a21_loc, &
                 K_a22_loc,&
                 M_b_vect_loc, &
                 S_b1_loc, &
                 S_b2_loc)
            
         end do
      end do
   end do

   
   call sll_factorize_csr_matrix(es_mp%sll_csr_mat)



   es_mp%sll_csr_mat_source => new_csr_matrix_mp( &
        es_mp%solution_size_final, &
        es_mp%solution_size_final, &
        es_mp%T%number_patches,&
        es_mp%num_elts_no_null_patchs, &
        es_mp%local_to_global_row, &
        es_mp%total_num_splines_loc, &
        es_mp%local_to_global_col, &
        es_mp%total_num_splines_loc ) 


    call compute_Source_matrice_mp(es_mp,Source_loc)
    print*,' ok'

    SLL_DEALLOCATE_ARRAY(Source_loc,ierr)
    SLL_DEALLOCATE_ARRAY(M_c_loc,ierr)
    SLL_DEALLOCATE_ARRAY(K_a11_loc,ierr)
    SLL_DEALLOCATE_ARRAY(K_a12_loc,ierr)
    SLL_DEALLOCATE_ARRAY(K_a21_loc,ierr)
    SLL_DEALLOCATE_ARRAY(K_a22_loc,ierr)
    SLL_DEALLOCATE_ARRAY(M_b_vect_loc,ierr)
    SLL_DEALLOCATE_ARRAY(S_b1_loc,ierr)
    SLL_DEALLOCATE_ARRAY(S_b2_loc,ierr)
   
    time = sll_time_elapsed_since(t0)
    print*, '#time for factorize_mat_es', time

  end subroutine factorize_mat_es_mp
  


  !> @brief Assemble the matrix for elliptic solver.
  !> @details To have the function phi such that 
  !>  div( A grad phi ) + B grad phi + C phi = rho
  !>  where A is a matrix of functions , B a vectorial function,
  !>  and  C and rho a scalar function.  
  !>  A, B, C, rho can be discret or analytic. 
  !>  phi is given with a B-spline interpolator  
  !> 
  !> The parameters are
  !> @param es the type general_coordinate_elliptic_solver
  !> @param[in] rho the field corresponding to the source term   
  !> @param[out] phi the field corresponding to the solution of the equation
  !> @return phi the field solution of the equation
  
  subroutine solve_general_coordinates_elliptic_eq_mp(&
       es_mp,&
       rho,&
       phi)
    use sll_timer
    class(general_coordinate_elliptic_solver_mp) :: es_mp
    class(sll_scalar_field_multipatch_2d), intent(inout)  :: phi
    class(sll_scalar_field_multipatch_2d), intent(in),target  :: rho
    type(sll_logical_mesh_2d), pointer                         :: lm
    sll_int32 :: i,k
    sll_int32 :: j,ierr
    sll_int32 :: patch,sz_coef2,sz_coef1
    type(sll_time_mark)  :: t0 
    double precision :: time
    sll_int32 :: num_pts_g1, num_pts_g2
    class(sll_scalar_field_multipatch_2d),pointer  :: base_field_pointer
    sll_real64, dimension(:,:), pointer :: coeff_rho
    sll_real64, dimension(:), pointer :: rho_coeff_1d
    sll_real64, dimension(:), pointer :: tmp_phi_vec
    sll_int32, dimension(2) :: connectivity
    sll_int32 :: num_coef,i_spline_tot,index
    

    num_pts_g1 = size(es_mp%gauss_pts1,2)
    num_pts_g2 = size(es_mp%gauss_pts2,2)

    SLL_ALLOCATE(rho_coeff_1d(es_mp%solution_size_final),ierr)
    es_mp%rho_vec(:) = 0.0_f64
    rho_coeff_1d(:)  = 0.0_f64
   
    
    print*, maxval(es_mp%T%global_indices),es_mp%solution_size_final 

    call sll_set_time_mark(t0)
    
    k = 0

    base_field_pointer => rho
    select type( type_field => base_field_pointer)
    class is (sll_scalar_field_multipatch_2d)
       
       ! put the spline coefficients in a 1d array
       do patch = 1,es_mp%T%number_patches
          coeff_rho => type_field%fields(patch)%f%interp_2d%get_coefficients()
          lm => es_mp%T%get_logical_mesh(patch-1)
        
          sz_coef2 = type_field%interps(patch)%interp%size_coeffs2
          sz_coef1 = type_field%interps(patch)%interp%size_coeffs1

          do j=1,sz_coef2
             do i=1,sz_coef1

                num_coef = i + sz_coef1*(j-1)
                i_spline_tot = num_coef + (patch- 1)*sz_coef2*sz_coef1
                index = es_mp%T%get_spline_global_index(i_spline_tot)
                if (index >0) then
                   k = k + 1
                   rho_coeff_1d(index) = coeff_rho(i,j)
                end if
             end do
          end do
       end do

       es_mp%rho_vec(1:es_mp%solution_size_final) = 0.0_f64
       call sll_mult_csr_matrix_vector(&
            es_mp%sll_csr_mat_source,&
            rho_coeff_1d(1:es_mp%solution_size_final),es_mp%rho_vec(1:es_mp%solution_size_final))
       
       
    end select

       
    time = sll_time_elapsed_since(t0)

    print*, 'time to construct the rho', time

  

    es_mp%phi_vec(1:es_mp%solution_size_final) = 0.0_f64
    call solve_linear_system_mp(es_mp,es_mp%solution_size_final)
    !!    end if



    do patch = 1,es_mp%T%number_patches
       
       lm => es_mp%T%get_logical_mesh(patch-1)
       sz_coef2 = phi%interps(patch)%interp%size_coeffs2
       sz_coef1 = phi%interps(patch)%interp%size_coeffs1
       SLL_ALLOCATE(tmp_phi_vec(sz_coef1*sz_coef2),ierr)
       do j=1,sz_coef2
          do i=1,sz_coef1
             
             num_coef = i + sz_coef1*(j-1)
             i_spline_tot = num_coef + (patch- 1)*sz_coef2*sz_coef1
             index = es_mp%T%get_spline_global_index(i_spline_tot)

             if ( index == 0) then 
                tmp_phi_vec(num_coef) = 0
             else 
                tmp_phi_vec(num_coef) = es_mp%phi_vec(index)
             end if
             
          end do
       end do
       call phi%fields(patch)%f%interp_2d%set_coefficients(tmp_phi_vec(1:sz_coef1*sz_coef2))
       SLL_DEALLOCATE(tmp_phi_vec,ierr)
    end do


    
    SLL_DEALLOCATE_ARRAY(rho_coeff_1d,ierr)
    
!    SLL_DEALLOCATE_ARRAY(M_rho_loc,ierr)
!!$    SLL_DEALLOCATE_ARRAY(rho_at_gauss,ierr)
  end subroutine solve_general_coordinates_elliptic_eq_mp
  
  ! This is based on the assumption that all the input fields have the same
  ! boundary conditions. TO DO: put all the boundary condition parameters in
  ! a single module called 'boundary_condition_convention' or something, which
  ! can be used library-wide, this way we could extract this information 
  ! directly from the fields without any difficulties. 

  subroutine build_local_matrices_mp( &
       es_mp, &
       cell_index,&
       cell_i, &
       cell_j, &
       patch,&
       a11_field_mat, &
       a12_field_mat, &
       a21_field_mat, &
       a22_field_mat, &
       b1_field_vect, &
       b2_field_vect, &
       c_field, &
       M_c_loc, &
       K_a11_loc, &
       K_a12_loc, &
       K_a21_loc, &
       K_a22_loc, &
       M_b_vect_loc, &
       S_b1_loc,  &
       S_b2_loc,&
       Source_loc)

    
    class(general_coordinate_elliptic_solver_mp) :: es_mp
    type(sll_logical_mesh_2d), pointer           :: lm
    sll_int32, intent(in) :: cell_i
    sll_int32, intent(in) :: cell_j
    sll_int32, intent(in) :: cell_index
    sll_int32,intent(in) :: patch 
    class(sll_scalar_field_multipatch_2d), pointer :: a11_field_mat
    class(sll_scalar_field_multipatch_2d), pointer :: a12_field_mat
    class(sll_scalar_field_multipatch_2d), pointer :: a21_field_mat
    class(sll_scalar_field_multipatch_2d), pointer :: a22_field_mat
    class(sll_scalar_field_multipatch_2d), pointer :: b1_field_vect
    class(sll_scalar_field_multipatch_2d), pointer :: b2_field_vect
    class(sll_scalar_field_multipatch_2d), pointer :: c_field
    sll_real64, dimension(:,:,:,:), intent(inout) :: Source_loc
    sll_real64, dimension(:,:), intent(out) :: M_c_loc
    sll_real64, dimension(:,:), intent(out) :: K_a11_loc
    sll_real64, dimension(:,:), intent(out) :: K_a12_loc
    sll_real64, dimension(:,:), intent(out) :: K_a21_loc
    sll_real64, dimension(:,:), intent(out) :: K_a22_loc
    sll_real64, dimension(:,:), intent(out) :: M_b_vect_loc
    sll_real64, dimension(:,:), intent(out) :: S_b1_loc
    sll_real64, dimension(:,:), intent(out) :: S_b2_loc

!!$    sll_int32 :: bc_left    
!!$    sll_int32 :: bc_right
!!$    sll_int32 :: bc_bottom    
!!$    sll_int32 :: bc_top  
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
    sll_real64, dimension(es_mp%spline_degree1(patch)+1,es_mp%spline_degree1(patch)+1) :: work1
    sll_real64, dimension(es_mp%spline_degree2(patch)+1,es_mp%spline_degree2(patch)+1) :: work2
    sll_real64, dimension(es_mp%spline_degree1(patch)+1,2) :: dbiatx1
    sll_real64, dimension(es_mp%spline_degree2(patch)+1,2) :: dbiatx2
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
    sll_int32 :: num_cells2,num_cells1
    
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
    lm => es_mp%T%get_logical_mesh(patch-1)
    delta1    = lm%delta_eta1 !mesh2d%delta_eta1
    delta2    = lm%delta_eta2 !! mesh2d%delta_eta2
    eta1_min  = lm%eta1_min  
    eta2_min  = lm%eta2_min
    num_cells2= lm%num_cells2
    num_cells1= lm%num_cells1
    tmp1      = (es_mp%spline_degree1(patch) + 1)/2
    tmp2      = (es_mp%spline_degree2(patch) + 1)/2

    num_pts_g1 = size(es_mp%gauss_pts1,3)
    num_pts_g2 = size(es_mp%gauss_pts2,3)
    
  
    eta1  = eta1_min + (cell_i-1)*delta1
    eta2  = eta2_min + (cell_j-1)*delta2
    

    do j=1,num_pts_g2
       ! rescale Gauss points to be in interval [eta2 ,eta2 +delta_eta2]
       ! the bottom edge of the cell.
       gpt2  = eta2  + 0.5_f64*delta2 * ( es_mp%gauss_pts2(patch,1,j) + 1.0_f64 )
       wgpt2 = 0.5_f64*delta2*es_mp%gauss_pts2(patch,2,j) !ATTENTION 0.5
       
          gtmp2 = gpt2
          local_spline_index2 = es_mp%spline_degree2(patch) + cell_j

       call interv(es_mp%knots2(patch,:),size(es_mp%T%transfs(1)%T%knots2), gpt2, left_y, mflag_y )


       call bsplvd( &
            es_mp%knots2(patch,:), &
            es_mp%spline_degree2(patch)+1,&
            gtmp2,&
            left_y,&
            work2,&
            dbiatx2,&
            2)

       ! we stocke the values of spline to construct the source term
       es_mp%tab_index_coeff2(patch,cell_j + num_cells2*(j-1)) = left_y

       
       do i=1,num_pts_g1
          ! rescale Gauss points to be in interval [eta1,eta1+delta1]
          gpt1  = eta1  + 0.5_f64*delta1 * ( es_mp%gauss_pts1(patch,1,i) + 1.0_f64 )
          wgpt1 = 0.5_f64*delta1*es_mp%gauss_pts1(patch,2,i)
          
             gtmp1   = gpt1
             local_spline_index1 = es_mp%spline_degree1(patch) + cell_i
             
      
          call interv ( es_mp%knots1(patch,:), size(es_mp%T%transfs(1)%T%knots1), gpt1,&
               left_x, mflag_x )

          call bsplvd(&
               es_mp%knots1(patch,:),&
               es_mp%spline_degree1(patch)+1,&
               gtmp1,&
               left_x,&
               work1,&
               dbiatx1,&
               2 )

          es_mp%tab_index_coeff1(patch,cell_i + num_cells1*(i-1)) = left_x

          val_c        = c_field%fields(patch)%f%value_at_point(gpt1,gpt2)
          val_a11      = a11_field_mat%fields(patch)%f%value_at_point(gpt1,gpt2)
          val_a12      = a12_field_mat%fields(patch)%f%value_at_point(gpt1,gpt2)
          val_a21      = a21_field_mat%fields(patch)%f%value_at_point(gpt1,gpt2)
          val_a22      = a22_field_mat%fields(patch)%f%value_at_point(gpt1,gpt2)

          val_b1       = b1_field_vect%fields(patch)%f%value_at_point(gpt1,gpt2)
          val_b1_der1  = b1_field_vect%fields(patch)%f%first_deriv_eta1_value_at_point(gpt1,gpt2)
          val_b1_der2  = b1_field_vect%fields(patch)%f%first_deriv_eta2_value_at_point(gpt1,gpt2)

          val_b2       = b2_field_vect%fields(patch)%f%value_at_point(gpt1,gpt2)
          val_b2_der1  = b2_field_vect%fields(patch)%f%first_deriv_eta1_value_at_point(gpt1,gpt2)
          val_b2_der2  = b2_field_vect%fields(patch)%f%first_deriv_eta2_value_at_point(gpt1,gpt2)
   
          jac_mat(:,:) = c_field%fields(patch)%f%get_jacobian_matrix(gpt1,gpt2)
          val_jac = jac_mat(1,1)*jac_mat(2,2) - jac_mat(1,2)*jac_mat(2,1)

          es_mp%values_jacobian(patch,cell_i + num_cells1*(i-1),cell_j + num_cells2*(j-1)) = val_jac
        
          es_mp%intjac = es_mp%intjac + wgpt2*wgpt1*val_jac
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
          do ii = 0,es_mp%spline_degree1(patch)
             do jj = 0,es_mp%spline_degree2(patch)
                
                index1  =  jj * ( es_mp%spline_degree1(patch) + 1 ) + ii + 1
                
                
                do iii = 0,es_mp%spline_degree1(patch)
                   do jjj = 0,es_mp%spline_degree2(patch)
                      
                      index2 =  jjj*(es_mp%spline_degree1(patch) + 1) + iii + 1
                
                      Source_loc(patch,cell_index,index1, index2) = &
                           Source_loc(patch,cell_index,index1, index2) + &
                           val_jac*wgpt1*wgpt2* &
                           dbiatx1(ii+1,1)*dbiatx1(iii+1,1)*  &
                           dbiatx2(jj+1,1)*dbiatx2(jjj+1,1)
                      
                      
                      
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
   
  end subroutine build_local_matrices_mp
  
  
  
  subroutine local_to_global_matrices_mp( &
       es_mp,&
       cell_index, &
       cell_i, &
       cell_j, &
       patch,&
       M_c_loc, &
       K_a11_loc, &
       K_a12_loc, &
       K_a21_loc, &
       K_a22_loc,&
       M_b_vect_loc, &
       S_b1_loc, &
       S_b2_loc)
    
    class(general_coordinate_elliptic_solver_mp)  :: es_mp
    type(sll_logical_mesh_2d), pointer           :: lm
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
    
  
    sll_int32 :: index1, index2, index3, index4
    sll_int32 :: index_coef1,index_coef2,index
    !sll_int32 :: index_coef1_phi,index_coef2_phi,index_phi
    sll_int32 :: i,j,mm, nn, b, bprime,x,y
    sll_int32 :: li_A, li_Aprime
    sll_real64 :: elt_mat_global
    sll_int32 :: nbsp,nbsp1
    sll_int32 :: patch
    sll_int32 :: num_cells1,num_cells2

    lm => es_mp%T%get_logical_mesh(patch-1)
    num_cells1 = lm%num_cells1
    num_cells2 = lm%num_cells2
    
    
    do mm = 0,es_mp%spline_degree2(patch)
       index3 = cell_j + mm
       

       
       do i = 0,es_mp%spline_degree1(patch)
          
          index1 = cell_i + i

             nbsp = num_cells1 + es_mp%spline_degree1(patch)

          x          =  index1 + (index3-1)*nbsp
          b          =  mm * ( es_mp%spline_degree1(patch) + 1 ) + i + 1
          li_A       =  es_mp%local_to_global_row(patch,b, cell_index)
          !li_AA= es_mp%local_to_global_spline_indices_test(b, cell_index)
         
          
          do nn = 0,es_mp%spline_degree2(patch)
             
             index4 = cell_j + nn
             

             
             do j = 0,es_mp%spline_degree1(patch)
                
                index2 = cell_i + j

                   nbsp1 = num_cells1 + es_mp%spline_degree1(patch)

                
                y         = index2 + (index4-1)*nbsp1
                bprime    =  nn * ( es_mp%spline_degree1(patch) + 1 ) + j + 1
                li_Aprime = es_mp%local_to_global_col(patch,bprime,cell_index)
                !li_AAprime = es_mp%local_to_global_spline_indices_test(bprime,cell_index)
                
                elt_mat_global = &
                     M_c_loc(b, bprime)     - &
                     K_a11_loc(b, bprime)   - &
                     K_a12_loc(b, bprime)   - &
                     K_a21_loc(b, bprime)   - &
                     K_a22_loc(b, bprime)   - &
                     M_b_vect_loc(b,bprime) - &
                     S_b1_loc( b, bprime)   - &
                     S_b2_loc( b, bprime)
                
    
                
                index_coef1 = es_mp%tab_index_coeff1(patch,cell_i)- es_mp%spline_degree1(patch) + i
                index_coef2 = es_mp%tab_index_coeff2(patch,cell_j)- es_mp%spline_degree2(patch) + mm
                index = index_coef1 + (index_coef2-1)*nbsp!(num_cells1+1)
                es_mp%local_to_global_spline_indices_source_row(patch,b,cell_index)= index
!!$                
                es_mp%local_to_global_spline_indices_source_col(patch,bprime,cell_index)= y!index_phi
                
                if ( (li_A > 0) .and. (li_Aprime > 0) ) then
                   
                   call sll_add_to_csr_matrix( &
                        es_mp%sll_csr_mat, &
                        elt_mat_global, &
                        li_A, &
                        li_Aprime)
                end if
               
                
             end do
             
          end do
       end do
    end do
    
  end subroutine local_to_global_matrices_mp
  
  
  
  subroutine solve_linear_system_mp( es_mp,k )
    ! CSR_MAT*phi = rho_vec is the linear system to be solved. The solution
    ! is given in terms of the spline coefficients that represent phi.
    class(general_coordinate_elliptic_solver_mp) :: es_mp
    !type(csr_matrix)  :: csr_masse
    integer :: k
   
    call solve_gen_elliptic_eq_mp(es_mp,es_mp%rho_vec(1:k),es_mp%phi_vec(1:k))


  end subroutine solve_linear_system_mp
  
  subroutine solve_gen_elliptic_eq_mp(es_mp,apr_B,apr_U)
    class(general_coordinate_elliptic_solver_mp) :: es_mp
    sll_real64, dimension(:) :: apr_U
    sll_real64, dimension(:) :: apr_B 
    sll_int32  :: ai_maxIter
    sll_real64 :: ar_eps
    
    ar_eps = 1.d-13
    ai_maxIter = 100000
  
    call sll_solve_csr_matrix(es_mp%sll_csr_mat, apr_B, apr_U)
   
  end subroutine solve_gen_elliptic_eq_mp
  

  subroutine compute_Source_matrice_mp(es_mp,Source_loc)
    type(general_coordinate_elliptic_solver_mp),intent(inout) :: es_mp
    type(sll_logical_mesh_2d), pointer           :: lm
    sll_real64, dimension(:,:,:,:), pointer :: Source_loc
    sll_int32 :: cell_j,cell_i
    sll_int32 :: cell_index
    sll_int32 :: ideg2,ideg1
    sll_int32 :: jdeg2,jdeg1
    sll_int32 :: b, bprime
    sll_int32 :: li_A,li_Aprime
    sll_real64:: elt_mat_global
    sll_int32 :: num_cells2,num_cells1
    sll_int32 :: patch

    do patch = 1, es_mp%T%number_patches
       lm => es_mp%T%get_logical_mesh(patch-1)

       num_cells2 = lm%num_cells2
       num_cells1 = lm%num_cells1

       
       do cell_j=1,num_cells2
          do cell_i=1,num_cells1
             
             cell_index = cell_i+num_cells1*(cell_j-1)
             
             do ideg2 = 0,es_mp%spline_degree2(patch)
                
                do ideg1 = 0,es_mp%spline_degree1(patch)
                   
                   b          =  ideg2 * ( es_mp%spline_degree1(patch) + 1 ) + ideg1 + 1
                   li_A       = es_mp%local_to_global_row(patch,b, cell_index)
                  
                   
                   do jdeg2 = 0,es_mp%spline_degree2(patch)
                      
                      do jdeg1 = 0,es_mp%spline_degree1(patch)
                         
                         bprime    =  jdeg2 * ( es_mp%spline_degree1(patch) + 1 ) + jdeg1 + 1
                         li_Aprime = es_mp%local_to_global_col(patch,bprime,cell_index)
                         
                         
                         elt_mat_global = Source_loc(patch,cell_index,bprime,b)
                  
                         
                         if ( (li_A > 0) .and. (li_Aprime > 0)) then
                            
                            call sll_add_to_csr_matrix( &
                                 es_mp%sll_csr_mat_source, &
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
    end do

end subroutine compute_Source_matrice_mp

end module sll_general_coordinate_elliptic_solver_multipatch_module


 
