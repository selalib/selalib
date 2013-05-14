module sll_general_coordinate_qn_solver_module
#include "sll_working_precision.h"
  use sll_module_scalar_field_2d_base
  implicit none

contains

  subroutine solve_quasi_neutral_eq_general_coords( &
       a_field_mat, &
       c_field, &
       rho, &
       quadrature_type, &    ! ugly, should have more general way
        phi )

    class(sll_scalar_field_2d_base), dimension(:,:), intent(in) :: a_field_mat
    class(sll_scalar_field_2d_base), intent(in)                 :: c_field
    class(sll_scalar_field_2d_base), intent(in)                 :: rho
    sll_int32, intent(in)                                :: quadrature_type 
    class(sll_scalar_field_2d_base), intent(out)                :: phi
    sll_int32                                 :: quadrature_degree
    sll_real64, dimension(:,:), allocatable :: M_rho_loc
    sll_real64, dimension(:,:), allocatable :: M_c_loc
    sll_real64, dimension(:,:), allocatable :: K_a11_loc
    sll_real64, dimension(:,:), allocatable :: K_a12_loc
    sll_real64, dimension(:,:), allocatable :: K_a21_loc
    sll_real64, dimension(:,:), allocatable :: K_a22_loc
    sll_int32 :: ierr
    ! Check arguments for consistency, errors, etc.

    ! First step: Build the stiffness matrix and the mass matrix, which are
    ! computed at the same time.

    ! The quadrature degree is the number of splines that intersect a cell.
    quadrature_degree = c_field%interpolation_degree()+1

    SLL_ALLOCATE(M_rho_loc(quadrature_degree,quadrature_degree),ierr)
    SLL_ALLOCATE(M_c_loc(quadrature_degree,quadrature_degree),ierr)
    SLL_ALLOCATE(K_a11_loc(quadrature_degree,quadrature_degree),ierr)
    SLL_ALLOCATE(K_a12_loc(quadrature_degree,quadrature_degree),ierr)
    SLL_ALLOCATE(K_a21_loc(quadrature_degree,quadrature_degree),ierr)
    SLL_ALLOCATE(K_a22_loc(quadrature_degree,quadrature_degree),ierr)

  end subroutine solve_quasi_neutral_eq_general_coords

 ! subroutine build_global_mass_stiffness_matrices( &


end module sll_general_coordinate_qn_solver_module
