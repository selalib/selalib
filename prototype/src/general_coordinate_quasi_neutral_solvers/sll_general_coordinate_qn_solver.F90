module sll_general_coordinate_qn_solver_module
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
  use GC
  use sll_boundary_condition_descriptors
  use sll_module_scalar_field_2d_base
  use sll_module_scalar_field_2d_alternative
  use sparsematrix_module
  use connectivity_module
  use sll_knots
  use gauss_legendre_integration
  use sll_timer
  !use LU
  implicit none


  type :: general_coordinate_qn_solver
     sll_int32 :: total_num_splines_loc
     sll_int32 :: total_num_splines_eta1
     sll_int32 :: total_num_splines_eta2
     sll_int32 :: num_cells1
     sll_int32 :: num_cells2
     sll_real64, dimension(:), pointer :: knots1
     sll_real64, dimension(:), pointer :: knots2
     sll_real64, dimension(:,:), pointer :: gauss_pts1
     sll_real64, dimension(:,:), pointer :: gauss_pts2
     sll_int32 :: bc_left
     sll_int32 :: bc_right
     sll_int32 :: bc_bottom
     sll_int32 :: bc_top
     sll_int32 :: spline_degree1
     sll_int32 :: spline_degree2
     sll_real64 :: epsi
     ! the following is otherwise known as "ID" in Aurore Back's original
     ! nomenclature. The indexing of the
     ! splines in this array depends on the boundary conditions.
     sll_int32, dimension(:), pointer :: global_spline_indices 
     ! the following is otherwise known as "IEN"
     sll_int32, dimension(:,:), pointer :: local_spline_indices
     ! the following is otherwise known as "LM". Same as global_spline_indices
     ! but including the changes resulting from the boundary conditions.
     ! This is:
     ! local_to_global_spline_indices(i,j) = 
     !   global_spline_indices(local_spline_indices(i,j))
     sll_int32, dimension(:,:), pointer :: local_to_global_spline_indices
     type(csr_matrix) :: csr_mat
     sll_real64, dimension(:), pointer :: rho_vec
     sll_real64, dimension(:), pointer :: phi_vec
     sll_real64, dimension(:), pointer :: tmp_rho_vec
     sll_real64, dimension(:), pointer :: masse
     sll_real64, dimension(:), pointer :: stiff
  end type general_coordinate_qn_solver
  
  ! For the integration mode.  
  integer, parameter :: QNS_GAUSS_LEGENDRE = 0, QNS_GAUSS_LOBATTO = 1
  
  interface delete
     module procedure delete_qns
  end interface delete

  interface initialize
     module procedure initialize_general_qn_solver
  end interface initialize
  
  
contains ! *******************************************************************
  
  subroutine initialize_general_qn_solver( &
       qns, &
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
       eta2_max)
    
   type(general_coordinate_qn_solver), intent(out) :: qns
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
   sll_int32 :: knots1_size
   sll_int32 :: knots2_size
   sll_int32 :: num_splines1
   sll_int32 :: num_splines2
   sll_int32 :: vec_sz ! for rho_vec and phi_vec allocations
   sll_int32 :: ierr,ierr1
   sll_int32 :: solution_size,i

   qns%total_num_splines_loc = (spline_degree_eta1+1)*(spline_degree_eta2+1)
   ! The total number of splines in a single direction is given by
   ! num_cells + spline_degree
   num_splines1 = num_cells_eta1 + spline_degree_eta1
   num_splines2 = num_cells_eta2 + spline_degree_eta2
   SLL_ALLOCATE(qns%global_spline_indices(num_splines1*num_splines2),ierr)
   qns%global_spline_indices(:) = 0
   SLL_ALLOCATE(qns%local_spline_indices((spline_degree_eta1+1)*(spline_degree_eta2+1),(num_cells_eta1*num_cells_eta2)),ierr)
   qns%local_spline_indices(:,:) = 0
   SLL_ALLOCATE(qns%local_to_global_spline_indices((spline_degree_eta1+1)*(spline_degree_eta2+1),(num_cells_eta1*num_cells_eta2)),ierr)
   qns%local_to_global_spline_indices = 0
   ! This should be changed to verify that the passed BC's are part of the
   ! recognized list described in sll_boundary_condition_descriptors...
   qns%bc_left   = bc_left
   qns%bc_right  = bc_right
   qns%bc_bottom = bc_bottom
   qns%bc_top    = bc_top
   qns%spline_degree1 = spline_degree_eta1
   qns%spline_degree2 = spline_degree_eta2
   qns%num_cells1 = num_cells_eta1
   qns%num_cells2 = num_cells_eta2

   ! Allocate and fill the gauss points/weights information.
   ! First direction
   select case(quadrature_type1)
      case (QNS_GAUSS_LEGENDRE)
         SLL_ALLOCATE(qns%gauss_pts1(2,spline_degree_eta1+2),ierr)
         qns%gauss_pts1(:,:) = gauss_points(spline_degree_eta1+2)
      case (QNS_GAUSS_LOBATTO)
         print *, 'new_general_qn_solver(): not implemented gauss_lobatto ',&
              'because the interface of that function is not good.'
   end select

   select case(quadrature_type2)
      case (QNS_GAUSS_LEGENDRE)
         SLL_ALLOCATE(qns%gauss_pts2(2,spline_degree_eta2+2),ierr)
         qns%gauss_pts2(:,:) = gauss_points(spline_degree_eta2+2)
      case (QNS_GAUSS_LOBATTO)
         print *, 'new_general_qn_solver(): not implemented gauss_lobatto ',&
              'because the interface of that function is not good.'

   end select


   if( (bc_left == SLL_PERIODIC) .and. (bc_right == SLL_PERIODIC) .and. &
       (bc_bottom == SLL_PERIODIC) .and. (bc_top == SLL_PERIODIC) ) then

      qns%total_num_splines_eta1 = num_cells_eta1 
      qns%total_num_splines_eta2 = num_cells_eta2
   !   print*, 'ACHTUNG=',2*spline_degree_eta1+2,2*spline_degree_eta2+2
      knots1_size = 2*spline_degree_eta1+2
      knots2_size = 2*spline_degree_eta2+2
      vec_sz      = num_cells_eta1*num_cells_eta2
   else if( (bc_left == SLL_PERIODIC) .and. (bc_right == SLL_PERIODIC) .and.&
       (bc_bottom == SLL_DIRICHLET) .and. (bc_top == SLL_DIRICHLET) ) then

      qns%total_num_splines_eta1 = num_cells_eta1 
      qns%total_num_splines_eta2 = num_cells_eta2 + &
                                   spline_degree_eta2 - 2
      knots1_size = 2*spline_degree_eta1+2
      knots2_size = 2*spline_degree_eta2+num_cells_eta2+1
      vec_sz      = num_cells_eta1*(num_cells_eta2+spline_degree_eta2)
   else if( (bc_left == SLL_DIRICHLET) .and. (bc_right == SLL_DIRICHLET) .and.&
            (bc_bottom == SLL_PERIODIC) .and. (bc_top == SLL_PERIODIC) ) then

      qns%total_num_splines_eta1 = num_cells_eta1 + spline_degree_eta1 - 2
      qns%total_num_splines_eta2 = num_cells_eta2 
      knots1_size = 2*spline_degree_eta1+num_cells_eta1+1
      knots2_size = 2*spline_degree_eta2+2
      vec_sz      = (num_cells_eta1 + spline_degree_eta1)*num_cells_eta2
   else if( (bc_left == SLL_DIRICHLET) .and. (bc_right == SLL_DIRICHLET) .and.&
       (bc_bottom == SLL_DIRICHLET) .and. (bc_top == SLL_DIRICHLET) ) then

      qns%total_num_splines_eta1 = num_cells_eta1 + spline_degree_eta1 - 2
      qns%total_num_splines_eta2 = num_cells_eta2 + spline_degree_eta2 - 2
      knots1_size = 2*spline_degree_eta1 + num_cells_eta1+1
      knots2_size = 2*spline_degree_eta2 + num_cells_eta2+1
      vec_sz      = (num_cells_eta1 + spline_degree_eta1)*&
                    (num_cells_eta2 + spline_degree_eta2)
   end if

  ! print*, 'ACHTUNG=',knots1_size,knots2_size
   solution_size = qns%total_num_splines_eta1*qns%total_num_splines_eta2
   SLL_ALLOCATE(qns%knots1(knots1_size),ierr1)
   SLL_ALLOCATE(qns%knots2(knots2_size),ierr)
   SLL_ALLOCATE(qns%rho_vec(vec_sz),ierr)
   SLL_ALLOCATE(qns%phi_vec(solution_size),ierr)
   SLL_ALLOCATE(qns%tmp_rho_vec(solution_size),ierr)
   SLL_ALLOCATE(qns%masse(solution_size),ierr)
   SLL_ALLOCATE(qns%stiff(solution_size),ierr)
   qns%rho_vec(:) = 0.0
   qns%phi_vec(:) = 0.0
   qns%masse(:) = 0.0
   qns%stiff(:) = 0.0
   do i = 1, knots1_size
      qns%knots1(i) = 0.0
   !   print*, i,  qns%knots1(i)
   end do

   call initialize_knots( &
        spline_degree_eta1, &
        num_cells_eta1, &
        eta1_min, &
        eta1_max, &
        bc_left, &
        bc_right, &
        qns%knots1 )

   call initialize_knots( &
        spline_degree_eta2, &
        num_cells_eta2, &
        eta2_min, &
        eta2_max, &
        bc_bottom, &
        bc_top, &
        qns%knots2 )
 
   call initconnectivity( &
        num_cells_eta1, &
        num_cells_eta2, &
        spline_degree_eta1, &
        spline_degree_eta2, &
        bc_left, &
        bc_right, &
        bc_bottom, &
        bc_top, &
        qns%local_spline_indices, &
        qns%global_spline_indices, &
        qns%local_to_global_spline_indices )

    call create_CSR( &
        qns%csr_mat, &
        solution_size, &
        solution_size, &
        num_cells_eta1*num_cells_eta2, &
        qns%local_to_global_spline_indices, &
        qns%total_num_splines_loc, &
        qns%local_to_global_spline_indices, &
        qns%total_num_splines_loc )
  end subroutine initialize_general_qn_solver

  function new_general_qn_solver( &
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
   eta2_max ) result(qns)

   type(general_coordinate_qn_solver), pointer :: qns
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
   sll_int32 :: knots1_size
   sll_int32 :: knots2_size
   sll_int32 :: num_splines1
   sll_int32 :: num_splines2
   sll_int32 :: vec_sz ! for rho_vec and phi_vec allocations
   sll_int32 :: ierr
   sll_int32 :: solution_size

   SLL_ALLOCATE(qns,ierr)
   call initialize( &
        qns, &
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
 end function new_general_qn_solver
 
 subroutine delete_qns( qns )
   type(general_coordinate_qn_solver) :: qns
   sll_int32 :: ierr
    ! it is not good to check some cases and not others, fix...
   if(associated(qns%knots1)) then
      SLL_DEALLOCATE(qns%knots1,ierr)
   else
       print *, 'delete qns, WARNING: knots1 array was not allocated.'
    end if
    if(associated(qns%knots2)) then
       SLL_DEALLOCATE(qns%knots2,ierr)
    else
       print *, 'delete qns general coords, ', &
            'WARNING: knots2 array was not allocated.'
    end if
    SLL_DEALLOCATE(qns%gauss_pts1,ierr)
    SLL_DEALLOCATE(qns%gauss_pts2,ierr)
    SLL_DEALLOCATE(qns%global_spline_indices,ierr)
    SLL_DEALLOCATE(qns%local_spline_indices,ierr)
    SLL_DEALLOCATE(qns%local_to_global_spline_indices,ierr)
    call free_csr(qns%csr_mat)
    SLL_DEALLOCATE(qns%rho_vec,ierr)
    SLL_DEALLOCATE(qns%phi_vec,ierr)
    SLL_DEALLOCATE(qns%tmp_rho_vec,ierr)
    SLL_DEALLOCATE(qns%masse,ierr)
    SLL_DEALLOCATE(qns%stiff,ierr)
  end subroutine delete_qns



  subroutine factorize_mat_qns(&
       qns, &
       a11_field_mat, &
       a12_field_mat,&
       a21_field_mat,&
       a22_field_mat,&
       c_field, &
       rho)
    
    type(general_coordinate_qn_solver),intent(inout) :: qns
    class(sll_scalar_field_2d_base_ptr), intent(in) :: a11_field_mat
    class(sll_scalar_field_2d_base_ptr), intent(in) :: a12_field_mat
    class(sll_scalar_field_2d_base_ptr), intent(in) :: a21_field_mat
    class(sll_scalar_field_2d_base_ptr), intent(in) :: a22_field_mat
    class(sll_scalar_field_2d_base), intent(in)     :: c_field
    class(sll_scalar_field_2d_base), intent(in)     :: rho
    ! local
     sll_real64, dimension(:), allocatable   :: M_rho_loc
    sll_real64, dimension(:,:), allocatable :: M_c_loc
    sll_real64, dimension(:,:), allocatable :: K_a11_loc
    sll_real64, dimension(:,:), allocatable :: K_a12_loc
    sll_real64, dimension(:,:), allocatable :: K_a21_loc
    sll_real64, dimension(:,:), allocatable :: K_a22_loc
    sll_real64, dimension(:,:), allocatable :: full_Matrix
    sll_real64, dimension(:), pointer  :: Masse_tot
    sll_real64, dimension(:), pointer  :: Stiff_tot
    !sll_real64, dimension(:,:), allocatable :: Masse_loc
    sll_real64, dimension(:), allocatable :: Masse_loc
    sll_real64, dimension(:), allocatable :: Stiff_loc
    !sll_int32, dimension(:), allocatable :: ipvt
    !sll_real64, dimension(:), allocatable :: z
    sll_int32 :: total_num_splines_eta1
    sll_int32 :: total_num_splines_eta2
    sll_int32 :: spline_degree_eta1
    sll_int32 :: spline_degree_eta2
    sll_int32 :: num_cells_eta1
    sll_int32 :: num_cells_eta2
    sll_int32 :: total_num_splines_loc
    sll_int32 :: ierr,ierr1
    sll_int32 :: i
    sll_int32 :: j
    sll_int32 :: cell_index
    type(sll_logical_mesh_2d), pointer :: mesh
    !    type(sll_time_mark) :: timer
    sll_real64 :: time,res,eta1,eta2
    character(len=*),parameter :: as_file1='mat'
    integer :: li_ios,li_ios1

    total_num_splines_loc = qns%total_num_splines_loc
    SLL_ALLOCATE(M_rho_loc(total_num_splines_loc),ierr)
    SLL_ALLOCATE(M_c_loc(total_num_splines_loc,total_num_splines_loc),ierr)
    SLL_ALLOCATE(K_a11_loc(total_num_splines_loc,total_num_splines_loc),ierr)
    SLL_ALLOCATE(K_a12_loc(total_num_splines_loc,total_num_splines_loc),ierr)
    SLL_ALLOCATE(K_a21_loc(total_num_splines_loc,total_num_splines_loc),ierr)
    SLL_ALLOCATE(K_a22_loc(total_num_splines_loc,total_num_splines_loc),ierr)
    !SLL_ALLOCATE(Masse_loc(total_num_splines_loc,total_num_splines_loc),ierr)
    SLL_ALLOCATE(Masse_loc(total_num_splines_loc),ierr)
    SLL_ALLOCATE(Stiff_loc(total_num_splines_loc),ierr)

      !   Allocation full_Matrix 
    SLL_ALLOCATE(full_Matrix(qns%total_num_splines_eta1*qns%total_num_splines_eta2,qns%total_num_splines_eta2*qns%total_num_splines_eta1),ierr1)
    SLL_ALLOCATE( Masse_tot(qns%total_num_splines_eta1*qns%total_num_splines_eta2),ierr)
    SLL_ALLOCATE( Stiff_tot(qns%total_num_splines_eta1*qns%total_num_splines_eta2),ierr)
    full_Matrix(:,:) = 0.0_f64
    Masse_tot(:) = 0.0_f64
    Stiff_tot(:) = 0.0_f64
    Masse_loc(:) = 0.0_f64
    Stiff_loc(:) = 0.0_f64

    

    full_Matrix(:,:) = 0.0_f64
    mesh => c_field%get_logical_mesh( )

    do j=1,qns%num_cells2
       do i=1,qns%num_cells1

          
          cell_index = i+qns%num_cells1*(j-1)
          
          call build_local_matrices( &
               qns, &
               i, &
               j, &
               mesh, &
               a11_field_mat, &
               a12_field_mat, &
               a21_field_mat, &
               a22_field_mat, &
               c_field, &
               rho, &
               M_rho_loc, &
               Masse_loc,&
               Stiff_loc,&
               M_c_loc, &
               K_a11_loc, &
               K_a12_loc, &
               K_a21_loc, &
               K_a22_loc )
          
          call local_to_global_matrices( &
               qns, &
               cell_index, &
               i, &
               j, &
               M_rho_loc, &
               Masse_loc,&
               Stiff_loc,&
               M_c_loc, &
               K_a11_loc, &
               K_a12_loc, &
               K_a21_loc, &
               K_a22_loc,&
               full_Matrix,&
               qns%masse,&
               qns%stiff)
          
       end do
    end do


   
    SLL_DEALLOCATE_ARRAY(M_rho_loc,ierr)
    SLL_DEALLOCATE_ARRAY(M_c_loc,ierr)
    SLL_DEALLOCATE_ARRAY(K_a11_loc,ierr)
    SLL_DEALLOCATE_ARRAY(K_a12_loc,ierr)
    SLL_DEALLOCATE_ARRAY(K_a21_loc,ierr)
    SLL_DEALLOCATE_ARRAY(K_a22_loc,ierr) 
    SLL_DEALLOCATE_ARRAY(Stiff_loc,ierr) 
    SLL_DEALLOCATE_ARRAY(full_Matrix,ierr)
  end subroutine factorize_mat_qns
  
  
  subroutine solve_qns(&
       qns,&
       rho,&
       phi)
    
    type(general_coordinate_qn_solver) :: qns
    type(sll_scalar_field_2d_discrete_alt), intent(inout)  :: phi
    class(sll_scalar_field_2d_base), intent(in)     :: rho
    type(sll_logical_mesh_2d), pointer :: mesh
    sll_int32 :: i
    sll_int32 :: j,ierr
    sll_int32 :: cell_index
    sll_int32 :: total_num_splines_loc
    sll_real64, dimension(:), allocatable   :: M_rho_loc
    
    total_num_splines_loc = qns%total_num_splines_loc
    SLL_ALLOCATE(M_rho_loc(total_num_splines_loc),ierr)

    
    M_rho_loc = 0.0
    qns%rho_vec(:) = 0.0
    mesh => phi%get_logical_mesh( )
    !    call set_time_mark(timer) ! comment this
    ! loop over domain cells build local matrices M_c_loc 
    do j=1,qns%num_cells2
       do i=1,qns%num_cells1
          ! cells are numbered in a linear fashion, convert from (i,j) indexing
          ! to the linear array index.
          
          cell_index = i+qns%num_cells1*(j-1)
          
          call build_local_matrices_rho( &
               qns, &
               i, &
               j, &
               mesh, &
               rho, &
               M_rho_loc)
          
          call local_to_global_matrices_rho( &
               qns, &
               cell_index, &
               i, &
               j, &
               M_rho_loc)
          
       end do
    end do
    
    
    if ((qns%bc_bottom==SLL_PERIODIC).and.(qns%bc_top==SLL_PERIODIC) &
         .and. (qns%bc_bottom==SLL_PERIODIC).and.(qns%bc_top==SLL_PERIODIC)) then
       
       
       call solve_linear_system_perper(qns,mesh,qns%masse)
    else 
       
       call solve_linear_system(qns,mesh)
    end if
    
    call  phi%interp_2d%set_coefficients( qns%phi_vec)
  end subroutine solve_qns
  
  subroutine solve_quasi_neutral_eq_general_coords( &
       qns, &
       a11_field_mat, &
       a12_field_mat,&
       a21_field_mat,&
       a22_field_mat,&
       c_field, &
       rho, &
       phi )
    
    type(general_coordinate_qn_solver) :: qns
    class(sll_scalar_field_2d_base_ptr), intent(in) :: a11_field_mat
    class(sll_scalar_field_2d_base_ptr), intent(in) :: a12_field_mat
    class(sll_scalar_field_2d_base_ptr), intent(in) :: a21_field_mat
    class(sll_scalar_field_2d_base_ptr), intent(in) :: a22_field_mat
    class(sll_scalar_field_2d_base), intent(in)     :: c_field
    class(sll_scalar_field_2d_base), intent(in)     :: rho
    type(sll_scalar_field_2d_discrete_alt), intent(inout)  :: phi
    type(csr_matrix)  :: csr_masse
    sll_real64 :: epsi
    sll_real64, dimension(:), allocatable   :: M_rho_loc
    sll_real64, dimension(:,:), allocatable :: M_c_loc
    sll_real64, dimension(:,:), allocatable :: K_a11_loc
    sll_real64, dimension(:,:), allocatable :: K_a12_loc
    sll_real64, dimension(:,:), allocatable :: K_a21_loc
    sll_real64, dimension(:,:), allocatable :: K_a22_loc
    sll_real64, dimension(:,:), allocatable :: full_Matrix
    sll_real64, dimension(:), pointer  :: Masse_tot
    sll_real64, dimension(:), pointer  :: Stiff_tot
    !sll_real64, dimension(:,:), allocatable :: Masse_loc
    sll_real64, dimension(:), allocatable :: Masse_loc
    sll_real64, dimension(:), allocatable :: Stiff_loc
    !sll_int32, dimension(:), allocatable :: ipvt
    !sll_real64, dimension(:), allocatable :: z
    sll_int32 :: total_num_splines_eta1
    sll_int32 :: total_num_splines_eta2
    sll_int32 :: spline_degree_eta1
    sll_int32 :: spline_degree_eta2
    sll_int32 :: num_cells_eta1
    sll_int32 :: num_cells_eta2
    sll_int32 :: total_num_splines_loc
    sll_int32 :: ierr,ierr1
    sll_int32 :: i
    sll_int32 :: j
    sll_int32 :: cell_index
    type(sll_logical_mesh_2d), pointer :: mesh
    !    type(sll_time_mark) :: timer
    sll_real64 :: time,res,eta1,eta2
    character(len=*),parameter :: as_file1='mat'
    integer :: li_ios,li_ios1
   ! sll_real64 :: cond
   ! sll_int32  :: D, code
    ! This function buils adnd solves a system:
    !
    !      A*phi_vec = rho_vec
    !
    ! Where A is a matrix in the CSR (compressed sparse row) format.
    
    ! Check arguments for consistency, errors, etc.
    
    ! First step: Build the stiffness matrix and the mass matrix, which are
    ! computed at the same time.
    ! total number of splines should come in the field...
    
    ! The quadrature degree is the number of splines that intersect a cell.
 !   call set_time_mark(timer)
    total_num_splines_loc = qns%total_num_splines_loc
    SLL_ALLOCATE(M_rho_loc(total_num_splines_loc),ierr)
    SLL_ALLOCATE(M_c_loc(total_num_splines_loc,total_num_splines_loc),ierr)
    SLL_ALLOCATE(K_a11_loc(total_num_splines_loc,total_num_splines_loc),ierr)
    SLL_ALLOCATE(K_a12_loc(total_num_splines_loc,total_num_splines_loc),ierr)
    SLL_ALLOCATE(K_a21_loc(total_num_splines_loc,total_num_splines_loc),ierr)
    SLL_ALLOCATE(K_a22_loc(total_num_splines_loc,total_num_splines_loc),ierr)
    !SLL_ALLOCATE(Masse_loc(total_num_splines_loc,total_num_splines_loc),ierr)
    SLL_ALLOCATE(Masse_loc(total_num_splines_loc),ierr)
    SLL_ALLOCATE(Stiff_loc(total_num_splines_loc),ierr)

      !   Allocation full_Matrix 
    SLL_ALLOCATE(full_Matrix(qns%total_num_splines_eta1*qns%total_num_splines_eta2,qns%total_num_splines_eta2*qns%total_num_splines_eta1),ierr1)
    SLL_ALLOCATE( Masse_tot(qns%total_num_splines_eta1*qns%total_num_splines_eta2),ierr)
    SLL_ALLOCATE( Stiff_tot(qns%total_num_splines_eta1*qns%total_num_splines_eta2),ierr)
    full_Matrix(:,:) = 0.0_f64
    Masse_tot(:) = 0.0_f64
    Stiff_tot(:) = 0.0_f64
    Masse_loc(:) = 0.0_f64
    Stiff_loc(:) = 0.0_f64

    !!!!!!!!!!!! atttention rajout !!!!!!!
    qns%rho_vec(:) = 0.0
    !print*, qns%rho_vec(:)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    full_Matrix(:,:) = 0.0_f64
    mesh => c_field%get_logical_mesh( )
    !    call set_time_mark(timer) ! comment this
    ! loop over domain cells build local matrices M_c_loc 
    do j=1,qns%num_cells2
       do i=1,qns%num_cells1
          ! cells are numbered in a linear fashion, convert from (i,j) indexing
          ! to the linear array index.
          
          cell_index = i+qns%num_cells1*(j-1)
          
          call build_local_matrices( &
               qns, &
               i, &
               j, &
               mesh, &
               a11_field_mat, &
               a12_field_mat, &
               a21_field_mat, &
               a22_field_mat, &
               c_field, &
               rho, &
               M_rho_loc, &
               Masse_loc,&
               Stiff_loc,&
               M_c_loc, &
               K_a11_loc, &
               K_a12_loc, &
               K_a21_loc, &
               K_a22_loc )
          
          call local_to_global_matrices( &
               qns, &
               cell_index, &
               i, &
               j, &
               M_rho_loc, &
               Masse_loc,&
               Stiff_loc,&
               M_c_loc, &
               K_a11_loc, &
               K_a12_loc, &
               K_a21_loc, &
               K_a22_loc,&
               full_Matrix,&
               Masse_tot,&
               Stiff_tot)
          
       end do
    end do


    if ((qns%bc_bottom==SLL_PERIODIC).and.(qns%bc_top==SLL_PERIODIC) &
         .and. (qns%bc_bottom==SLL_PERIODIC).and.(qns%bc_top==SLL_PERIODIC)) then
       
       
       call solve_linear_system_perper(qns,mesh,Masse_tot)
    else 

       call solve_linear_system(qns,mesh)
    end if


    

    call  phi%interp_2d%set_coefficients( qns%phi_vec)
    
    !SLL_DEALLOCATE_ARRAY(full_Matrix,ierr)
    ! apr_B is the source, apr_U is the solution
    SLL_DEALLOCATE_ARRAY(M_rho_loc,ierr)
    SLL_DEALLOCATE_ARRAY(M_c_loc,ierr)
    SLL_DEALLOCATE_ARRAY(K_a11_loc,ierr)
    SLL_DEALLOCATE_ARRAY(K_a12_loc,ierr)
    SLL_DEALLOCATE_ARRAY(K_a21_loc,ierr)
    SLL_DEALLOCATE_ARRAY(K_a22_loc,ierr)
    SLL_DEALLOCATE_ARRAY(Masse_loc,ierr)
    SLL_DEALLOCATE_ARRAY(Masse_tot,ierr) 
    SLL_DEALLOCATE_ARRAY(Stiff_tot,ierr) 
    SLL_DEALLOCATE_ARRAY(Stiff_loc,ierr) 
    SLL_DEALLOCATE_ARRAY(full_Matrix,ierr)
  end subroutine solve_quasi_neutral_eq_general_coords
  
  ! This is based on the assumption that all the input fields have the same
  ! boundary conditions. TO DO: put all the boundary condition parameters in
  ! a single module called 'boundary_condition_convention' or something, which
  ! can be used library-wide, this way we could extract this information 
  ! directly from the fields without any difficulties. 
  subroutine build_local_matrices( &
       obj, &
       cell_i, &
       cell_j, &
       mesh2d, &
       a11_field_mat, &
       a12_field_mat, &
       a21_field_mat, &
       a22_field_mat, &
       c_field, &
       rho, &
       M_rho_loc, &
       Masse_loc,&
       Stiff_loc,&
       M_c_loc, &
       K_a11_loc, &
       K_a12_loc, &
       K_a21_loc, &
       K_a22_loc )
    !    use sll_constants

    type(general_coordinate_qn_solver) :: obj
    sll_int32, intent(in) :: cell_i
    sll_int32, intent(in) :: cell_j
    type(sll_logical_mesh_2d), pointer :: mesh2d
    class(sll_scalar_field_2d_base_ptr), intent(in) :: a11_field_mat
    class(sll_scalar_field_2d_base_ptr), intent(in) :: a12_field_mat
    class(sll_scalar_field_2d_base_ptr), intent(in) :: a21_field_mat
    class(sll_scalar_field_2d_base_ptr), intent(in) :: a22_field_mat
    class(sll_scalar_field_2d_base), intent(in)     :: c_field
    class(sll_scalar_field_2d_base), intent(in)     :: rho
    sll_real64 :: epsi
    sll_real64, dimension(:), intent(out)   :: M_rho_loc
    sll_real64, dimension(:,:), intent(out) :: M_c_loc
    sll_real64, dimension(:,:), intent(out) :: K_a11_loc
    sll_real64, dimension(:,:), intent(out) :: K_a12_loc
    sll_real64, dimension(:,:), intent(out) :: K_a21_loc
    sll_real64, dimension(:,:), intent(out) :: K_a22_loc
    sll_real64, dimension(:), intent(out) :: Masse_loc
    sll_real64, dimension(:), intent(out) :: Stiff_loc
    !sll_real64, dimension(:), intent(out) :: Masse_loc
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
    sll_real64, dimension(obj%spline_degree1+1,obj%spline_degree1+1) :: work1
    sll_real64, dimension(obj%spline_degree2+1,obj%spline_degree2+1) :: work2
    sll_real64, dimension(obj%spline_degree1+1,2) :: dbiatx1
    sll_real64, dimension(obj%spline_degree2+1,2) :: dbiatx2
    sll_real64 :: val_f
    sll_real64 :: val_c
    sll_real64 :: val_a11
    sll_real64 :: val_a12
    sll_real64 :: val_a21
    sll_real64 :: val_a22
    sll_real64, dimension(2,2) :: jac_mat
    sll_real64 :: val_jac
    sll_real64 :: B11
    sll_real64 :: B12
    sll_real64 :: B21
    sll_real64 :: B22

    Masse_loc(:) = 0.0
    Stiff_loc(:) = 0.0
    M_rho_loc(:) = 0.0
    M_c_loc(:,:) = 0.0
    K_a11_loc(:,:) = 0.0
    K_a12_loc(:,:) = 0.0
    K_a21_loc(:,:) = 0.0
    K_a22_loc(:,:) = 0.0
    dbiatx1(:,:)  = 0.0_f64
    dbiatx2(:,:)  = 0.0_f64
    work1(:,:) = 0.0_f64
    work2(:,:) = 0.0_f64
    ! The supposition is that all fields use the same logical mesh
    delta1    = mesh2d%delta_eta1
    delta2    = mesh2d%delta_eta2
    eta1_min  = mesh2d%eta1_min
    eta2_min  = mesh2d%eta2_min
    tmp1      = (obj%spline_degree1 + 1)/2
    tmp2      = (obj%spline_degree2 + 1)/2
    bc_left   = obj%bc_left
    bc_right  = obj%bc_right
    bc_bottom = obj%bc_bottom
    bc_top    = obj%bc_top
    num_pts_g1 = size(obj%gauss_pts1,2) !obj%spline_degree1+2
    num_pts_g2 = size(obj%gauss_pts2,2)!obj%spline_degree2+2

   ! print*, delta1,delta2,eta1_min,eta2_min,mesh2d%eta1_max,mesh2d%eta2_max
   ! print*, bc_left,bc_right,bc_bottom,bc_top
  
    eta1  = eta1_min + (cell_i-1)*delta1
    eta2  = eta2_min + (cell_j-1)*delta2

  !  print*, 'point base',eta1,eta2,num_pts_g1,num_pts_g2
    do j=1,num_pts_g2
       ! rescale Gauss points to be in interval [eta2 ,eta2 +delta_eta2]
       ! the bottom edge of the cell.
       gpt2  = eta2  + 0.5_f64*delta2 * ( obj%gauss_pts2(1,j) + 1.0_f64 )
       wgpt2 = 0.5_f64*delta2*obj%gauss_pts2(2,j) !ATTENTION 0.5

       if ((obj%bc_bottom==SLL_PERIODIC).and.(obj%bc_top==SLL_PERIODIC))then
          ! rescale gauss point in interval [0,delta2]
          gtmp2 = 0.5_f64*delta2*( obj%gauss_pts2(1,j) + 1.0_f64) !ATTENTION 0.5
          local_spline_index2 = obj%spline_degree2 + 1

       else if ((obj%bc_bottom == SLL_DIRICHLET).and.&
            (obj%bc_top    == SLL_DIRICHLET)) then
          gtmp2 = gpt2
          local_spline_index2 = obj%spline_degree2 + cell_j
       end if

!
!!$
       call bsplvd( &
            obj%knots2, &
            obj%spline_degree2+1,&
            gtmp2,&
            local_spline_index2,&
            work2,&
            dbiatx2,&
            2)

       !print*, 'splin2=',dbiatx2
       !print*, 'knot2',obj%knots2
       do i=1,num_pts_g1
          ! rescale Gauss points to be in interval [eta1,eta1+delta1]
          gpt1  = eta1  + 0.5_f64*delta1 * ( obj%gauss_pts1(1,i) + 1.0_f64 )
          wgpt1 = 0.5_f64*delta1*obj%gauss_pts1(2,i)
          
          if((obj%bc_left==SLL_PERIODIC).and.(obj%bc_right==SLL_PERIODIC)) then 
             
             gtmp1   = 0.5_f64*delta1*( obj%gauss_pts1(1,i) + 1.0_f64)! ATTENTION 0.5 
             local_spline_index1 = obj%spline_degree1 + 1

          else if ((obj%bc_left  == SLL_DIRICHLET).and.&
               (obj%bc_right == SLL_DIRICHLET) ) then
             
             gtmp1   = gpt1
             local_spline_index1 = obj%spline_degree1 + cell_i
             
          end if
     !     print*,  'gauss',obj%gauss_pts1(1,i)

          
          call bsplvd(&
               obj%knots1,&
               obj%spline_degree1+1,&
               gtmp1,&
               local_spline_index1,&
               work1,&
               dbiatx1,&
               2 )

      !    print*, 'point=',gtmp1,gtmp2,local_spline_index1
          !print*, 'splin1=',dbiatx1
          !print*, 'knot1',obj%knots1
         ! print*, 'rho'
          val_f   =rho%value_at_point(gpt1,gpt2)! 0.05*cos(0.5*gpt1)
         ! print*, 'valeur rho=',val_f,2*(2.0*sll_pi)**2*cos(2.0*sll_pi*gpt1)*cos(2.0*sll_pi*gpt2)
         ! print*, 'val',gpt1,gpt2,val_f, 0.05*cos(0.5*gpt1), val_f-0.05*cos(0.5*gpt1)
       
          val_c   = c_field%value_at_point(gpt1,gpt2)
          !print*, 'val,',val_c
          val_a11 = a11_field_mat%base%value_at_point(gpt1,gpt2)
          val_a12 = a12_field_mat%base%value_at_point(gpt1,gpt2)
          val_a21 = a21_field_mat%base%value_at_point(gpt1,gpt2)
          val_a22 = a22_field_mat%base%value_at_point(gpt1,gpt2)
          !print*,'matrix values', val_a11,val_a12,val_a21,val_a22
          jac_mat(:,:) = c_field%get_jacobian_matrix(gpt1,gpt2)
          val_jac = jac_mat(1,1)*jac_mat(2,2) - jac_mat(1,2)*jac_mat(2,1)!abs(jac_mat(1,1)*jac_mat(2,2) - jac_mat(1,2)*jac_mat(2,1))
          !print*, 'determinant', val_jac
          ! The B matrix is  by (J^(-1)) A^T (J^(-1))^T 
          B11 = jac_mat(2,2)*jac_mat(2,2)*val_a11 - &
               jac_mat(2,2)*jac_mat(1,2)*(val_a12+val_a21) + &
               jac_mat(1,2)*jac_mat(1,2)*val_a22
         ! print*, 'test11', B11
!!$          B12 = jac_mat(1,1)*jac_mat(2,2)*val_a21 - &
!!$               jac_mat(1,1)*jac_mat(1,2)*val_a22 - &
!!$               jac_mat(2,1)*jac_mat(2,2)*val_a11 + &
!!$               jac_mat(1,2)*jac_mat(2,1)*val_a12

          B12 = jac_mat(1,1)*jac_mat(2,2)*val_a12 - &
               jac_mat(1,1)*jac_mat(1,2)*val_a22 - &
               jac_mat(2,1)*jac_mat(2,2)*val_a11 + &
               jac_mat(1,2)*jac_mat(2,1)*val_a21

          !print*, 'test12', B12
!!$          B21 = jac_mat(1,1)*jac_mat(2,2)*val_a12 - &
!!$               jac_mat(1,1)*jac_mat(1,2)*val_a22 - &
!!$               jac_mat(2,1)*jac_mat(2,2)*val_a11 + &
!!$               jac_mat(1,2)*jac_mat(2,1)*val_a21

          B21 = jac_mat(1,1)*jac_mat(2,2)*val_a21 - &
               jac_mat(1,1)*jac_mat(1,2)*val_a22 - &
               jac_mat(2,1)*jac_mat(2,2)*val_a11 + &
               jac_mat(1,2)*jac_mat(2,1)*val_a12
         !print*, 'test21', B21
          B22 = jac_mat(1,1)*jac_mat(1,1)*val_a22 - &
               jac_mat(1,1)*jac_mat(2,1)*(val_a21+val_a12) + &
               jac_mat(2,1)*jac_mat(2,1)*val_a11
          !print*, 'test22', B22
          ! loop over the splines supported in the cell that are different than
          ! zero at the point (gpt1,gpt2) (there are spline_degree+1 splines in
          ! each direction.
          do ii = 0,obj%spline_degree1
             do jj = 0,obj%spline_degree2
                
                index1  =  jj * ( obj%spline_degree1 + 1 ) + ii + 1
                M_rho_loc(index1)= M_rho_loc(index1) + &
                     val_f*val_jac*wgpt1*wgpt2* &
                     dbiatx1(ii+1,1)*dbiatx2(jj+1,1)

                

                Masse_loc(index1) = &
                     Masse_loc(index1) + &
                     val_jac*wgpt1*wgpt2* &
                     (dbiatx1(ii+1,1)*dbiatx2(jj+1,1))

                Stiff_loc(index1) = &
                     Stiff_loc(index1) + &
                     val_jac*wgpt1*wgpt2* &
                     (dbiatx1(ii+1,2)*dbiatx2(jj+1,1)+dbiatx1(ii+1,1)*dbiatx2(jj+1,2))
                
                
                !print*, 'ethop',M_rho_loc(index1),dbiatx1(ii+1,1),dbiatx2(jj+1,1)
                
                do iii = 0,obj%spline_degree1
                   do jjj = 0,obj%spline_degree2
                      
                      index2 =  jjj*(obj%spline_degree1 + 1) + iii + 1
                      
!!$                      Masse_loc(index1, index2) = &
!!$                           Masse_loc(index1, index2) + &
!!$                           val_jac*wgpt1*wgpt2* &
!!$                           dbiatx1(ii+1,1)*dbiatx1(iii+1,1)*  &
!!$                           dbiatx2(jj+1,1)*dbiatx2(jjj+1,1)
!!$                      
                      
                      
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
                   end do
                end do
             end do
          end do
       end do
    end do
    !print*,  '????????????????????????'
    !print*,  'rez',   Masse_loc
    ! print*, 'rho_loc=',cell_i,cell_j, M_rho_loc
  end subroutine build_local_matrices
  

  subroutine build_local_matrices_rho( &
       obj, &
       cell_i, &
       cell_j, &
       mesh2d, &
       rho, &
       M_rho_loc)
    !    use sll_constants

    type(general_coordinate_qn_solver) :: obj
    sll_int32, intent(in) :: cell_i
    sll_int32, intent(in) :: cell_j
    type(sll_logical_mesh_2d), pointer :: mesh2d
    class(sll_scalar_field_2d_base), intent(in)     :: rho
    sll_real64 :: epsi
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
    sll_real64, dimension(obj%spline_degree1+1,obj%spline_degree1+1) :: work1
    sll_real64, dimension(obj%spline_degree2+1,obj%spline_degree2+1) :: work2
    sll_real64, dimension(obj%spline_degree1+1,2) :: dbiatx1
    sll_real64, dimension(obj%spline_degree2+1,2) :: dbiatx2
    sll_real64 :: val_f
    sll_real64, dimension(2,2) :: jac_mat
    sll_real64 :: val_jac
   
   
    M_rho_loc(:) = 0.0
    dbiatx1(:,:)  = 0.0_f64
    dbiatx2(:,:)  = 0.0_f64
    work1(:,:) = 0.0_f64
    work2(:,:) = 0.0_f64
    ! The supposition is that all fields use the same logical mesh
    delta1    = mesh2d%delta_eta1
    delta2    = mesh2d%delta_eta2
    eta1_min  = mesh2d%eta1_min
    eta2_min  = mesh2d%eta2_min
    tmp1      = (obj%spline_degree1 + 1)/2
    tmp2      = (obj%spline_degree2 + 1)/2
    bc_left   = obj%bc_left
    bc_right  = obj%bc_right
    bc_bottom = obj%bc_bottom
    bc_top    = obj%bc_top
    num_pts_g1 = size(obj%gauss_pts1,2) !obj%spline_degree1+2
    num_pts_g2 = size(obj%gauss_pts2,2)!obj%spline_degree2+2

   ! print*, delta1,delta2,eta1_min,eta2_min,mesh2d%eta1_max,mesh2d%eta2_max
   ! print*, bc_left,bc_right,bc_bottom,bc_top
  
    eta1  = eta1_min + (cell_i-1)*delta1
    eta2  = eta2_min + (cell_j-1)*delta2

  !  print*, 'point base',eta1,eta2,num_pts_g1,num_pts_g2
    do j=1,num_pts_g2
       ! rescale Gauss points to be in interval [eta2 ,eta2 +delta_eta2]
       ! the bottom edge of the cell.
       gpt2  = eta2  + 0.5_f64*delta2 * ( obj%gauss_pts2(1,j) + 1.0_f64 )
       wgpt2 = 0.5_f64*delta2*obj%gauss_pts2(2,j) !ATTENTION 0.5
       
       if ((obj%bc_bottom==SLL_PERIODIC).and.(obj%bc_top==SLL_PERIODIC))then
          ! rescale gauss point in interval [0,delta2]
          gtmp2 = 0.5_f64*delta2*( obj%gauss_pts2(1,j) + 1.0_f64) !ATTENTION 0.5
          local_spline_index2 = obj%spline_degree2 + 1
          
       else if ((obj%bc_bottom == SLL_DIRICHLET).and.&
            (obj%bc_top    == SLL_DIRICHLET)) then
          gtmp2 = gpt2
          local_spline_index2 = obj%spline_degree2 + cell_j
       end if
       ! gtmp2 = eta2  + 0.5_f64*delta2 * ( obj%gauss_pts2(1,j) + 1.0_f64 ) !- 0.0001
   !    print*, 'test',gtmp2, eta2, delta2, obj%gauss_pts2(1,j)
       !    print*, 'splin2=',dbiatx2,obj%spline_degree2,local_spline_index2
       !    print*, 'knot',obj%knots2
       ! dbiatx2(:,:) = 0.0
       !
!!$
       call bsplvd( &
            obj%knots2, &
            obj%spline_degree2+1,&
            gtmp2,&
            local_spline_index2,&
            work2,&
            dbiatx2,&
            2)
       
       !print*, 'splin2=',dbiatx2
       !print*, 'knot2',obj%knots2
       do i=1,num_pts_g1
          ! rescale Gauss points to be in interval [eta1,eta1+delta1]
          gpt1  = eta1  + 0.5_f64*delta1 * ( obj%gauss_pts1(1,i) + 1.0_f64 )
          wgpt1 = 0.5_f64*delta1*obj%gauss_pts1(2,i)
          
          if((obj%bc_left==SLL_PERIODIC).and.(obj%bc_right==SLL_PERIODIC)) then 
             
             gtmp1   = 0.5_f64*delta1*( obj%gauss_pts1(1,i) + 1.0_f64)! ATTENTION 0.5 
             local_spline_index1 = obj%spline_degree1 + 1
             
          else if ((obj%bc_left  == SLL_DIRICHLET).and.&
               (obj%bc_right == SLL_DIRICHLET) ) then
             
             gtmp1   = gpt1
             local_spline_index1 = obj%spline_degree1 + cell_i
             
          end if
          !     print*,  'gauss',obj%gauss_pts1(1,i)
          
          
          call bsplvd(&
               obj%knots1,&
               obj%spline_degree1+1,&
               gtmp1,&
               local_spline_index1,&
               work1,&
               dbiatx1,&
               2 )
          val_f   =rho%value_at_point(gpt1,gpt2)! 0.05*cos(0.5*gpt1)
          ! print*, 'valeur rho=',val_f,2*(2.0*sll_pi)**2*cos(2.0*sll_pi*gpt1)*cos(2.0*sll_pi*gpt2)
          ! print*, 'val',gpt1,gpt2,val_f, 0.05*cos(0.5*gpt1), val_f-0.05*cos(0.5*gpt1)
          
          jac_mat(:,:) = rho%get_jacobian_matrix(gpt1,gpt2)
          val_jac = jac_mat(1,1)*jac_mat(2,2) - jac_mat(1,2)*jac_mat(2,1)!abs(jac_mat(1,1)*jac_mat(2,2) - jac_mat(1,2)*jac_mat(2,1))
          
          ! loop over the splines supported in the cell that are different than
          ! zero at the point (gpt1,gpt2) (there are spline_degree+1 splines in
          ! each direction.
          do ii = 0,obj%spline_degree1
             do jj = 0,obj%spline_degree2
                
                index1  =  jj * ( obj%spline_degree1 + 1 ) + ii + 1
                M_rho_loc(index1)= M_rho_loc(index1) + &
                     val_f*val_jac*wgpt1*wgpt2* &
                     dbiatx1(ii+1,1)*dbiatx2(jj+1,1)
                
             end do
          end do
       end do
    end do
    
  end subroutine build_local_matrices_rho

  subroutine local_to_global_matrices( &
       qns,&
       cell_index, &
       cell_i, &
       cell_j, &
       M_rho_loc, &
       Masse_loc,&
       Stiff_loc,&
       M_c_loc, &
       K_a11_loc, &
       K_a12_loc, &
       K_a21_loc, &
       K_a22_loc,&
       full_Matrix,&
       Masse_tot,&
       Stiff_tot)
    
    type(general_coordinate_qn_solver)  :: qns
    sll_int32 :: cell_index
    sll_int32 :: cell_i
    sll_int32 :: cell_j
    sll_real64, dimension(:), intent(in)   :: M_rho_loc
    sll_real64, dimension(:,:), intent(in) :: M_c_loc
    sll_real64, dimension(:,:), intent(in) :: K_a11_loc
    sll_real64, dimension(:,:), intent(in) :: K_a12_loc
    sll_real64, dimension(:,:), intent(in) :: K_a21_loc
    sll_real64, dimension(:,:), intent(in) :: K_a22_loc
    !  Correspond to the full Matrix of linear system 
    !  It is not necessary to keep it  
    sll_real64, dimension(:,:), intent(inout) :: full_Matrix
    sll_real64, dimension(:), intent(in) :: Masse_loc
    sll_real64, dimension(:), intent(in) :: Stiff_loc
    sll_real64, dimension(:), intent(inout) :: Masse_tot
    sll_real64, dimension(:), intent(inout) :: Stiff_tot
    sll_int32 :: index1, index2, index3, index4
    sll_int32 :: i,j,mm, nn, b, bprime,x,y
    sll_int32 :: li_A, li_Aprime
    sll_real64 :: elt_mat_global
    sll_real64 :: elt_masse
    sll_int32 :: nbsp,nbsp1
    sll_int32 :: bc_left
    sll_int32 :: bc_right
    sll_int32 :: bc_bottom
    sll_int32 :: bc_top
    sll_int32 :: ierr 
    
    bc_left   = qns%bc_left
    bc_right  = qns%bc_right
    bc_bottom = qns%bc_bottom
    bc_top    = qns%bc_top

  

    do mm = 0,qns%spline_degree2
       index3 = cell_j + mm
       
       if (  (bc_bottom==SLL_PERIODIC).and.(bc_top== SLL_PERIODIC) ) then 
          
          if ( index3 > qns%total_num_splines_eta2) then
             index3 = index3 - qns%total_num_splines_eta2
          end if
          
       end if
       
       do i = 0,qns%spline_degree1
          
          index1 = cell_i + i
          if ( (bc_left==SLL_PERIODIC).and.(bc_right== SLL_PERIODIC)) then 
             if ( index1 > qns%total_num_splines_eta1) then
                
                index1 = index1 - qns%total_num_splines_eta1
                
             end if
             nbsp = qns%total_num_splines_eta1
             
          else if ( (bc_left  == SLL_DIRICHLET).and.&
               (bc_right == SLL_DIRICHLET) ) then
             nbsp = qns%num_cells1 + qns%spline_degree1
          end if
          x          =  index1 + (index3-1)*nbsp
          b          =  mm * ( qns%spline_degree1 + 1 ) + i + 1
          li_A       =  qns%local_to_global_spline_indices(b, cell_index)
          qns%rho_vec(x)  =  qns%rho_vec(x)  + M_rho_loc(b)
          !print*, x,b
          Masse_tot(x)    = Masse_tot(x) + Masse_loc(b)
          Stiff_tot(x)    = Stiff_tot(x) + Stiff_loc(b)
          do nn = 0,qns%spline_degree2
             
             index4 = cell_j + nn
             
             if ( (bc_bottom==SLL_PERIODIC).and.(bc_top== SLL_PERIODIC))then
                if ( index4 > qns%total_num_splines_eta2) then
                   
                   index4 = index4 - qns%total_num_splines_eta2
                end if
             end if
             
             do j = 0,qns%spline_degree1
                
                index2 = cell_i + j
                if((bc_left==SLL_PERIODIC).and.(bc_right==SLL_PERIODIC))then
                   
                   if ( index2 > qns%total_num_splines_eta1) then
                      
                      index2 = index2 - qns%total_num_splines_eta1
                   end if
                   nbsp1 = qns%total_num_splines_eta1
                   
                else if ( (bc_left  == SLL_DIRICHLET) .and.&
                     (bc_right == SLL_DIRICHLET) ) then
                   
                   nbsp1 = qns%num_cells1 + qns%spline_degree1
                end if
                   
                y         = index2 + (index4-1)*nbsp1
                bprime    =  nn * ( qns%spline_degree1 + 1 ) + j + 1
                li_Aprime = qns%local_to_global_spline_indices(bprime, &
                                                               cell_index)
                elt_mat_global = &
                     M_c_loc(b, bprime) - &
                     K_a11_loc(b, bprime) - &
                     K_a12_loc(b, bprime) - &
                     K_a21_loc(b, bprime) - &
                     K_a22_loc(b, bprime)

!!$                if  ((qns%bc_bottom==SLL_PERIODIC).and.(qns%bc_top==SLL_PERIODIC) &
!!$                     .and. (qns%bc_bottom==SLL_PERIODIC).and.(qns%bc_top==SLL_PERIODIC)) then
!!$                   
!!$                   elt_mat_global= elt_mat_global - Masse_loc(bprime)
!!$                  ! print*,  Masse_loc(bprime)
!!$                end if

               ! print*, 'elt_mat_global=', M_c_loc(b, bprime), K_a11_loc(b, bprime), K_a12_loc(b, bprime),K_a21_loc(b, bprime),K_a22_loc(b, bprime)

!!$                full_Matrix(x,y) = &
!!$                     full_Matrix(x,y) + &
!!$                     M_c_loc(b, bprime) - &
!!$                     K_a11_loc(b, bprime) - &
!!$                     K_a12_loc(b, bprime) - &
!!$                     K_a21_loc(b, bprime) - &
!!$                     K_a22_loc(b, bprime)
                !print*, 'elt', full_Matrix(x,y)
                ! elt_masse = Masse_loc(b,bprime)
                
                if ( (li_A > 0) .and. (li_Aprime > 0) ) then
                   call add_MVal(qns%csr_mat,elt_mat_global,li_A,li_Aprime)
                   ! call add_MVal(csr_masse,elt_masse,li_A,li_Aprime)
                end if
                
             end do
          end do
       end do
    end do
    !print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    !print *, 'a = ', csr_masse%opr_a(1:csr_masse%opi_ia(2)-1)
    
  end subroutine local_to_global_matrices


   subroutine local_to_global_matrices_rho( &
       qns,&
       cell_index, &
       cell_i, &
       cell_j, &
       M_rho_loc)
     
     type(general_coordinate_qn_solver)  :: qns
     sll_int32 :: cell_index
     sll_int32 :: cell_i
     sll_int32 :: cell_j
     sll_real64, dimension(:), intent(in)   :: M_rho_loc
    !  Correspond to the full Matrix of linear system 
    !  It is not necessary to keep it  
     sll_int32 :: index1, index3
    sll_int32 :: i,j,mm, nn, b, bprime,x,y
    sll_int32 :: nbsp,nbsp1
    sll_int32 :: bc_left
    sll_int32 :: bc_right
    sll_int32 :: bc_bottom
    sll_int32 :: bc_top
    sll_int32 :: ierr 
    
    bc_left   = qns%bc_left
    bc_right  = qns%bc_right
    bc_bottom = qns%bc_bottom
    bc_top    = qns%bc_top

    
    
    do mm = 0,qns%spline_degree2
       index3 = cell_j + mm
       
       if (  (bc_bottom==SLL_PERIODIC).and.(bc_top== SLL_PERIODIC) ) then 
          
          if ( index3 > qns%total_num_splines_eta2) then
             index3 = index3 - qns%total_num_splines_eta2
          end if
          
       end if
       
       do i = 0,qns%spline_degree1
          
          index1 = cell_i + i
          if ( (bc_left==SLL_PERIODIC).and.(bc_right== SLL_PERIODIC)) then 
             if ( index1 > qns%total_num_splines_eta1) then
                
                index1 = index1 - qns%total_num_splines_eta1
                
             end if
             nbsp = qns%total_num_splines_eta1
             
          else if ( (bc_left  == SLL_DIRICHLET).and.&
               (bc_right == SLL_DIRICHLET) ) then
             nbsp = qns%num_cells1 + qns%spline_degree1
          end if
          x          =  index1 + (index3-1)*nbsp
          b          =  mm * ( qns%spline_degree1 + 1 ) + i + 1
          qns%rho_vec(x)  =  qns%rho_vec(x)  + M_rho_loc(b)
          
       end do
    end do
    
  end subroutine local_to_global_matrices_rho

  subroutine solve_linear_system( qns,mesh2d )
    ! CSR_MAT*phi = rho_vec is the linear system to be solved. The solution
    ! is given in terms of the spline coefficients that represent phi.
    type(general_coordinate_qn_solver) :: qns
    type(sll_logical_mesh_2d), pointer :: mesh2d
    !type(csr_matrix)  :: csr_masse
    integer :: elt, elt1
    integer :: i,j
    sll_real64, dimension(:), allocatable :: resul
    integer :: li_ios
    character(len=*),parameter :: as_file='rho', as_file1='phi',as_file2='mat'
    sll_int32 :: bc_left
    sll_int32 :: bc_right
    sll_int32 :: bc_bottom
    sll_int32 :: bc_top
    qns%tmp_rho_vec = 0.0_f64
    bc_left   = qns%bc_left
    bc_right  = qns%bc_right
    bc_bottom = qns%bc_bottom
    bc_top    = qns%bc_top

    qns%tmp_rho_vec(:) = 0.0_f64
    !print*, 'source coeff=', qns%rho_vec
    if( (bc_left   == SLL_PERIODIC) .and. (bc_right == SLL_PERIODIC) .and. &
        (bc_bottom == SLL_DIRICHLET).and. (bc_top   == SLL_DIRICHLET) ) then
       
       do i = 1, qns%total_num_splines_eta1
          do j = 1, qns%total_num_splines_eta2
             
             elt  = i + qns%total_num_splines_eta1 * (  j - 1)
             elt1 = i + ( qns%total_num_splines_eta1 ) * j
             qns%tmp_rho_vec(elt) = qns%rho_vec(elt1)
          end do
       end do
       
    else if ( (bc_left   == SLL_DIRICHLET).and.(bc_right==SLL_DIRICHLET) .and.&
              (bc_bottom == SLL_DIRICHLET).and.(bc_top==SLL_DIRICHLET) ) then 
       
       do i = 1, qns%total_num_splines_eta1
          do j = 1, qns%total_num_splines_eta2

             elt  = i + qns%total_num_splines_eta1 * (  j - 1)
             elt1 = i + 1 + ( qns%total_num_splines_eta1 + 2 ) * j 
             qns%tmp_rho_vec( elt ) = qns%rho_vec( elt1 )
          end do
       end do
       
    else if((bc_left   == SLL_PERIODIC) .and. (bc_right==SLL_PERIODIC) .and.&
         (bc_bottom == SLL_PERIODIC) .and. (bc_top  ==SLL_PERIODIC)) then
       
       qns%tmp_rho_vec(1:qns%total_num_splines_eta1*qns%total_num_splines_eta2)=&
            qns%rho_vec(1:qns%total_num_splines_eta1*qns%total_num_splines_eta2) 
       
   
       
    else if( (bc_left == SLL_DIRICHLET) .and. (bc_right == SLL_DIRICHLET) .and.&
             (bc_bottom == SLL_PERIODIC).and. (bc_top   == SLL_PERIODIC) ) then
       
       do i = 1, qns%total_num_splines_eta1
          do j = 1, qns%total_num_splines_eta2

             elt1 = i + 1 + ( qns%total_num_splines_eta1 + 2 ) * (  j - 1)
             elt  = i + qns%total_num_splines_eta1 * (  j - 1)
             qns%tmp_rho_vec( elt ) = qns%rho_vec( elt1 )
          end do
       end do


    end if
    
   

    
    !print *, 'a = ', qns%csr_mat%opr_a(1:qns%csr_mat%opi_ia(2)-1)
    call solve_general_qn(qns,qns%csr_mat,qns%tmp_rho_vec,qns%phi_vec)
  
!!$    allocate(resul(qns%total_num_splines_eta1*qns%total_num_splines_eta2))
!!$    resul = 0.0
!!$    call  Mult_CSR_Matrix_Vector(qns%csr_mat, qns%phi_vec, resul)
!!$   
!!$    print*, 'difference A phi - rho=', sum(resul-qns%tmp_rho_vec)
!!$    print*, dot_product(resul,qns%phi_vec)
!!$    deallocate(resul) 
!!$
!!$    !call print_csrMatrix_mmFormat(qns%csr_mat,as_file2)
!!$    
!!$    open(unit=1, file=as_file//'.txt', iostat=li_ios)
!!$    
!!$    do i =1,qns%total_num_splines_eta2*qns%total_num_splines_eta1
!!$       write(1, *)qns%tmp_rho_vec(i)
!!$    end do
!!$    close(1)

!!$    open(unit=1, file=as_file1//'.txt', iostat=li_ios)
!!$    
!!$    do i =1,qns%total_num_splines_eta2*qns%total_num_splines_eta1
!!$       write(1, *)qns%phi_vec(i)
!!$       !print*, qns%phi_vec(i)
!!$    end do
!!$    close(1)
!!$    do i = 1,qns%total_num_splines_eta2
!!$       print*, '---------------',qns%phi_vec(qns%total_num_splines_eta1*(i-1)+1:qns%total_num_splines_eta1*i)
!!$       print*, ' source term :', qns%tmp_rho_vec(qns%total_num_splines_eta1*(i-1)+1:qns%total_num_splines_eta1*i)
!!$    end do
    
  end subroutine solve_linear_system

  subroutine solve_general_qn(qns,csr_mat,apr_B,apr_U)
    type(general_coordinate_qn_solver) :: qns
    type(csr_matrix) :: csr_mat
    sll_real64, dimension(:) :: apr_U
    sll_real64, dimension(:) :: apr_B 
    sll_int32  :: ai_maxIter
    sll_real64 :: ar_eps
    
    ar_eps = 1.d-13
    ai_maxIter = 100000
    !print*, ai_maxIter 
    !print *, 'a = ', csr_mat % opr_a(1:csr_mat%opi_ia(2)-1)
    !print*, 'tttttttttttttt'
    !print *, 'a = ', csr_masse % opr_a(1:csr_masse%opi_ia(2)-1)
    if ( (qns%bc_left == SLL_PERIODIC) .and. (qns%bc_right == SLL_PERIODIC) .and. &
         (qns%bc_bottom == SLL_PERIODIC) .and. (qns%bc_top == SLL_PERIODIC) ) then

!!$       call Gradient_conj_projet(&
!!$             csr_mat,&
!!$            apr_B,&
!!$            apr_U,&
!!$            ai_maxIter,&
!!$            ar_eps )
!!$       call Gradient_conj_perper(&
!!$            csr_mat,&
!!$            apr_B,&
!!$            apr_U,&
!!$            ai_maxIter,&
!!$            ar_eps )
       call Gradient_conj(&
            csr_mat,&
            apr_B,&
            apr_U,&
            ai_maxIter,&
            ar_eps )

!!$       call Gradient_conj_perper_variante(&
!!$            csr_mat,&
!!$            apr_B,&
!!$            apr_U,&
!!$            ai_maxIter,&
!!$            ar_eps,&
!!$            apr_B)
    else 
       call Gradient_conj(csr_mat,&
            apr_B,&
            apr_U,&
            ai_maxIter,&
            ar_eps )
    end if
    !print*,'u', apr_U
  end subroutine solve_general_qn



  subroutine solve_linear_system_perper( qns,mesh2d,Masse_tot )
    ! CSR_MAT*phi = rho_vec is the linear system to be solved. The solution
    ! is given in terms of the spline coefficients that represent phi.
    type(general_coordinate_qn_solver) :: qns
    type(sll_logical_mesh_2d), pointer :: mesh2d
    sll_real64, dimension(:),pointer :: Masse_tot

    
    qns%tmp_rho_vec = 0.0_f64
    
    qns%tmp_rho_vec(:) = 0.0_f64
    qns%tmp_rho_vec(1:qns%total_num_splines_eta1*qns%total_num_splines_eta2)=&
         qns%rho_vec(1:qns%total_num_splines_eta1*qns%total_num_splines_eta2) 
    !print*, 'INTEGRALE RHO', sum(qns%tmp_rho_vec)
    call solve_general_qn_perper(qns,qns%csr_mat,qns%tmp_rho_vec,qns%phi_vec,Masse_tot)
    
    !print*, 'INTEGRALE DE PHI= ', dot_product(Masse_tot,qns%phi_vec)
  end subroutine solve_linear_system_perper


  subroutine solve_general_qn_perper(qns,csr_mat,apr_B,apr_U,Masse_tot)
    type(general_coordinate_qn_solver) :: qns
    type(csr_matrix) :: csr_mat
    sll_real64, dimension(:) :: apr_U
    sll_real64, dimension(:) :: apr_B 
    sll_real64, dimension(:),pointer :: Masse_tot
    sll_int32  :: ai_maxIter
    sll_real64 :: ar_eps
    
    ar_eps = 1.d-13
    ai_maxIter = 100000
    call Gradient_conj_adjusted(&
         csr_mat,&
         apr_B,&
         apr_U,&
         Masse_tot,&
         ai_maxIter,&
         ar_eps )
    
  end subroutine solve_general_qn_perper

end module sll_general_coordinate_qn_solver_module


 
