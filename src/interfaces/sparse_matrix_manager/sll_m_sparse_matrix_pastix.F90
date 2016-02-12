!**************************************************************
!  Copyright INRIA
!  Authors : 
!     CALVI project team
!  
!  This code SeLaLib (for Semi-Lagrangian-Library) 
!  is a parallel library for simulating the plasma turbulence 
!  in a tokamak.
!  
!  This software is governed by the CeCILL-B license 
!  under French law and abiding by the rules of distribution 
!  of free software.  You can  use, modify and redistribute 
!  the software under the terms of the CeCILL-B license as 
!  circulated by CEA, CNRS and INRIA at the following URL
!  "http://www.cecill.info". 
!**************************************************************

module sll_m_sparse_matrix
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"

use sll_m_pastix
use sll_m_qsort_partition, only : sll_s_qsortc

implicit none

public :: &
  sll_f_new_csr_matrix, &
  sll_f_new_csr_matrix_with_constraint, &
  sll_s_csr_add_one_constraint, &
  sll_s_delete_csr_matrix, &
  sll_s_initialize_csr_matrix, &
  sll_s_initialize_csr_matrix_with_constraint, &
  sll_s_add_to_csr_matrix, &
  sll_t_csr_matrix, &
  sll_s_factorize_csr_matrix, &
  sll_s_mult_csr_matrix_vector, &
  sll_s_solve_csr_matrix, &
  sll_s_solve_csr_matrix_perper, &
  sll_o_delete

private

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!> @brief type for CSR format
type sll_t_csr_matrix
  sll_int32                         :: num_rows !< number of rows
  sll_int32                         :: num_cols !< number of columns
  sll_int32                         :: num_nz   !< number of non zero elements
  sll_int32,  dimension(:), pointer :: row_ptr
  sll_int32,  dimension(:), pointer :: col_ind
  sll_real64, dimension(:), pointer :: val
  type(pastix_solver),      pointer :: linear_solver
end type sll_t_csr_matrix

interface sll_o_delete
   module procedure sll_s_delete_csr_matrix
end interface sll_o_delete

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine sll_s_delete_csr_matrix(csr_mat)

  type(sll_t_csr_matrix), pointer :: csr_mat

  nullify(csr_mat)
    
end subroutine sll_s_delete_csr_matrix

!> @brief allocates the memory space for a new CSR type on the heap,
!> initializes it with the given arguments and returns a pointer to the
!> object.
!> param[in] num_rows :  number of rows
!> param[in] num_cols :  number of columns
!> param[in] num_element :  number of elements
!> param[in] local_to_global_row : local_to_global_row(\ell,i) gives the global 
!> row index of the matrix, for the element i and local degree of freedom \ell
!> param[in] num_local_dof_row : number of local degrees of freedom for the rows
!> param[in] local_to_global_col : local_to_global_col(\ell,i) gives the global 
!> column index of the matrix, for the element i and local degree of freedom \ell
!> param[in] num_local_dof_col : number of local degrees of freedom for the columns
!> return a pointer to the newly allocated object.

function sll_f_new_csr_matrix( num_rows,            &
                               num_cols,            &
                               num_elts,            &
                               local_to_global_row, &
                               num_local_dof_row,   &
                               local_to_global_col, &
                               num_local_dof_col)   result(mat)

  type(sll_t_csr_matrix), pointer           :: mat
  sll_int32,                     intent(in) :: num_rows
  sll_int32,                     intent(in) :: num_cols
  sll_int32,                     intent(in) :: num_elts
  sll_int32, dimension(:,:),     intent(in) :: local_to_global_row
  sll_int32,                     intent(in) :: num_local_dof_row
  sll_int32, dimension(:,:),     intent(in) :: local_to_global_col
  sll_int32,                     intent(in) :: num_local_dof_col

  sll_int32 :: ierr

  SLL_ALLOCATE(mat, ierr)

  call sll_s_initialize_csr_matrix( mat,                 &
                                    num_rows,            &
                                    num_cols,            &
                                    num_elts,        &
                                    local_to_global_row, &
                                    num_local_dof_row,   &
                                    local_to_global_col, &
                                    num_local_dof_col)
      
end function sll_f_new_csr_matrix

!> @brief initialization of CSR matrix type
!> thanks to the global index of each local dof of each element
!> param[inout] mat : CSR matrix structure
!> param[in]    num_rows :  number of rows
!> param[in]    num_cols :  number of columns
!> param[in]    num_element :  number of elements
!> param[in]    local_to_global_row : local_to_global_row(\ell,i) gives the global 
!> row index of the matrix, for the element i and local degree of freedom \ell
!> param[in] num_local_dof_row : number of local degrees of freedom for the rows
!> param[in] local_to_global_col : local_to_global_col(\ell,i) gives the global 
!> column index of the matrix, for the element i and local degree of freedom \ell
!> param[in] num_local_dof_col : number of local degrees of freedom for the columns
subroutine sll_s_initialize_csr_matrix(  mat,                 &
                                         num_rows,            &
                                         num_cols,            &
                                         num_elts,            &
                                         local_to_global_row, &
                                         num_local_dof_row,   &
                                         local_to_global_col, &
                                         num_local_dof_col)

  type(sll_t_csr_matrix),    intent(inout) :: mat
  sll_int32,                 intent(in)    :: num_rows
  sll_int32,                 intent(in)    :: num_cols
  sll_int32,                 intent(in)    :: num_elts
  sll_int32, dimension(:,:), intent(in)    :: local_to_global_row
  sll_int32,                 intent(in)    :: num_local_dof_row
  sll_int32, dimension(:,:), intent(in)    :: local_to_global_col
  sll_int32,                 intent(in)    :: num_local_dof_col

  sll_int32          :: num_nz
  sll_int32, pointer :: lpi_columns(:,:)
  sll_int32, pointer :: lpi_occ(:)
  sll_int32          :: coef
  sll_int32          :: ierr
  
  coef = 10
  SLL_ALLOCATE(lpi_columns(num_rows, 0:coef * num_local_dof_col),ierr)
  SLL_ALLOCATE(lpi_occ(num_rows + 1),ierr)

  lpi_columns(:,:) = 0
  lpi_occ(:) = 0
  ! COUNTING NON ZERO ELEMENTS

  num_nz = sll_count_non_zero_elts( num_rows,            &
                                    num_cols,            &
                                    num_elts,            &
                                    local_to_global_row, &
                                    num_local_dof_row,   &
                                    local_to_global_col, &
                                    num_local_dof_col,   &
                                    lpi_columns,         &
                                    lpi_occ              )
  mat%num_rows = num_rows
  mat%num_cols = num_cols
  mat%num_nz   = num_nz

  call init_sparsematrix( num_elts,                 &
                          local_to_global_row,      &
                          num_local_dof_row,        &  
                          local_to_global_col,      &
                          num_local_dof_col,        &
                          lpi_columns,              &
                          lpi_occ,                  &
                          colptr,                   &
                          row,                      &
                          num_rows                  )
  stop

  mat%linear_solver%iparm(IPARM_TRANSPOSE_SOLVE) = API_YES
  SLL_ALLOCATE(mat%linear_solver, ierr)
  call initialize(mat%linear_solver,num_rows,num_nz)

  mat%linear_solver%avals = 0._f64   
  
  SLL_DEALLOCATE_ARRAY(lpi_columns,ierr)
  SLL_DEALLOCATE_ARRAY(lpi_occ,ierr)
  
end subroutine sll_s_initialize_csr_matrix

subroutine sll_s_initialize_csr_matrix_with_constraint( mat, mat_a)

  type(sll_t_csr_matrix), intent(inout) :: mat
  type(sll_t_csr_matrix), intent(in) :: mat_a

  !print*,' COUNTING NON ZERO ELEMENTS'
  mat%num_nz = mat_a%num_nz + 2*mat_a%num_rows       
  print*,'num_nz mat, num_nz mat_tot', mat_a%num_nz,mat%num_nz 
  mat%num_rows = mat_a%num_rows  +  1
  print*,'num_rows mat, num_rows mat_tot',mat_a%num_rows , mat%num_rows
  mat%num_cols = mat_a%num_cols  +  1
  print*,'num_cols mat, num_cols mat_tot',mat_a%num_cols , mat%num_cols 

  mat%linear_solver%avals = 0._f64
  
end subroutine sll_s_initialize_csr_matrix_with_constraint

function sll_f_new_csr_matrix_with_constraint(mat_a) result(mat)

  type(sll_t_csr_matrix), pointer :: mat
  type(sll_t_csr_matrix)          :: mat_a

  sll_int32                     :: ierr

  SLL_ALLOCATE(mat, ierr)

  call sll_s_initialize_csr_matrix_with_constraint( mat, mat_a)

end function sll_f_new_csr_matrix_with_constraint
 
subroutine sll_s_csr_add_one_constraint( ia_in,          &
                                   ja_in,          &
                                   a_in,           &
                                   num_rows_in,    &
                                   num_nz_in,      &
                                   constraint_vec, &
                                   ia_out,         &
                                   ja_out,         &
                                   a_out)

  sll_int32, dimension(:), intent(in)  :: ia_in  
  sll_int32, dimension(:), intent(in)  :: ja_in  
  real(8), dimension(:), intent(in)  :: a_in
  sll_int32,               intent(in)  :: num_rows_in
  sll_int32,               intent(in)  :: num_nz_in
  real(8), dimension(:), intent(in)  :: constraint_vec
  sll_int32, dimension(:), intent(out) :: ia_out
  sll_int32, dimension(:), intent(out) :: ja_out
  real(8), dimension(:), intent(out) :: a_out

  sll_int32                            :: num_rows_out
  sll_int32                            :: num_nz_out
  sll_int32                            :: i
  sll_int32                            :: s
  sll_int32                            :: k
  
  
  num_rows_out = num_rows_in+1
  num_nz_out   = num_nz_in+2*num_rows_in
  
  if(size(ia_in)<num_rows_in+1) then
    print *, '#problem of size of ia_in', size(ia_in),num_rows_in+1
    stop
  endif
  if(size(ja_in)<num_nz_in) then
    print *, '#problem of size of ja_in', size(ja_in),num_nz_in
    stop
  endif
  if(size(a_in)<num_nz_in) then
    print *, '#problem of size of a_in', size(a_in),num_nz_in
    stop
  endif
  if(size(ia_out)<num_rows_out+1) then
    print *, '#problem of size of ia_out', size(ia_out),num_rows_out+1
    stop
  endif
  if(size(ja_out)<num_nz_out) then
    print *, '#problem of size of ja_out', size(ja_out),num_nz_out
    stop
  endif
  if(size(a_out)<num_nz_out) then
    print *, '#problem of size of a_out', size(a_out),num_nz_out
    stop
  endif
  if(ia_in(num_rows_in+1).ne.num_nz_in+1)then
    print *,'#bad value of ia_in(num_rows_in+1)', ia_in(num_rows_in+1),num_nz_in+1
    stop
  endif
  
  s = 1
  do i=1,num_rows_in
    ia_out(i) = s
    do k = ia_in(i), ia_in(i+1)-1
      a_out(s) = a_in(k)
      ja_out(s) = ja_in(k)
      s = s+1
    enddo
    a_out(s) = constraint_vec(i)
    ja_out(s) = num_rows_out
    s = s+1
  enddo
  ia_out(num_rows_in+1) = s
  do i=1,num_rows_in
    a_out(s) = constraint_vec(i)
    ja_out(s) = i
    s = s+1      
  enddo
  ia_out(num_rows_in+2) = s
   
  if(ia_out(num_rows_out+1).ne.num_nz_out+1)then
    print *,'#bad value of ia_out(num_rows_out+1)',ia_out(num_rows_out+1),num_nz_out+1
    stop
  endif
  
end subroutine sll_s_csr_add_one_constraint

subroutine sll_s_factorize_csr_matrix(mat)
  type(sll_t_csr_matrix), intent(inout) :: mat
  
  call factorize(mat%linear_solver)
  
end subroutine sll_s_factorize_csr_matrix


subroutine sll_s_mult_csr_matrix_vector(mat, input, output)
  
  type(sll_t_csr_matrix) :: mat
  sll_real64, dimension(:), intent(in) :: input
  sll_real64, dimension(:), intent(out) :: output
  !local var
  sll_int32 :: i
  sll_int32 :: k_1
  sll_int32 :: k_2

  do i = 1, mat % num_rows

    k_1 = mat%linear_solver%colptr(i)
    k_2 = mat%linear_solver%colptr(i + 1) - 1
    output(i) = &
      dot_product(mat%linear_solver%avals(k_1: k_2), &
        input(mat%linear_solver%row(k_1: k_2)))
          
  end do

end subroutine sll_s_mult_csr_matrix_vector


subroutine sll_s_add_to_csr_matrix(mat, val, a, a_prime)
  
  type(sll_t_csr_matrix), intent(inout) :: mat
  sll_real64,           intent(in)    :: val
  sll_int32,            intent(in)    :: a
  sll_int32,            intent(in)    :: a_prime
  !local var
  sll_int32 :: j
  sll_int32 :: k


  ! THE CURRENT LINE IS self%row_ptr(a)
  do k = mat%linear_solver%colptr(a), mat%linear_solver%colptr(a + 1) - 1
    j = mat%linear_solver%row(k)
    if (j == a_prime) then
      mat%linear_solver%avals(k) = mat%linear_solver%avals(k) + val 
      exit
    end if
  end do

end subroutine sll_s_add_to_csr_matrix


subroutine set_values_csr_matrix(mat, val)
  implicit none
  type(sll_t_csr_matrix) :: mat
  sll_real64, dimension(:), intent(in) :: val
  !local variables
  
  if(size(val)<mat%num_nz)then
    print *,'#Problem of size of val',size(val),mat%num_nz
    print *,'#at line',__LINE__
    print *,'#in file',__FILE__
    print *,'#in subroutine set_values_csr_matrix'
    stop
  endif
  
  mat%linear_solver%avals(1:mat%num_nz) = val(1:mat%num_nz)
      
end subroutine set_values_csr_matrix

subroutine sll_s_solve_csr_matrix(mat, apr_B, apr_U)
  implicit none
  type(sll_t_csr_matrix) :: mat
  sll_real64, dimension(:) :: apr_U
  sll_real64, dimension(:) :: apr_B
  !local var
  sll_int32  :: sys
  sys = 0
  apr_U = apr_B
  call solve(mat%linear_solver,apr_U)
  
end subroutine sll_s_solve_csr_matrix


    sll_int32 function sll_count_non_zero_elts( &
      ai_nR, ai_nC, ai_nel, api_LM_1, ai_nen_1, api_LM_2, ai_nen_2, api_columns, api_occ)
        ! _1 FOR ROWS
        ! _2 FOR COLUMNS
        implicit none
        sll_int32 :: ai_nR, ai_nC
        sll_int32, dimension(:,:), intent(in) :: api_LM_1, api_LM_2
        sll_int32 :: ai_nel, ai_nen_1, ai_nen_2
        sll_int32, dimension(:,:), pointer :: api_columns
        sll_int32, dimension(:), pointer :: api_occ
        !local var
        sll_int32 :: e, b_1, A_1, b_2, A_2, i
        sll_int32 :: result
        sll_int32, dimension(2) :: lpi_size
        logical :: ll_done
        sll_int32, dimension(:,:), pointer :: lpi_columns

        ! WE FIRST COMPUTE, FOR EACH ROW, THE NUMBER OF COLUMNS THAT WILL BE USED
        do e = 1, ai_nel

            do b_1 = 1, ai_nen_1

                A_1 = api_LM_1(b_1, e)
                if (A_1 == 0) then
                    cycle
                end if

                do b_2 = 1, ai_nen_2

                    A_2 = api_LM_2(b_2, e)
                    if (A_2 == 0) then
                        cycle
                    end if

                    ll_done = .false.
                    ! WE CHECK IF IT IS THE FIRST OCCURANCE OF THE COUPLE (A_1, A_2)
                    do i = 1, api_columns(A_1, 0)

                        if (api_columns(A_1, i) /= A_2) then
                            cycle
                        end if

                        ll_done = .true.
                        exit

                    end do

                    if (.not.ll_done) then

                        api_occ(A_1) = api_occ(A_1) + 1

                        ! A_1 IS THE ROW NUM, A_2 THE COLUMN NUM
                        ! INITIALIZATION OF THE SPARSE MATRIX
                        api_columns(A_1, 0) = api_columns(A_1, 0) + 1
                        api_columns(A_1, api_columns(A_1, 0)) = A_2

                        ! resizing the array
                        lpi_size(1) = SIZE(api_columns, 1)
                        lpi_size(2) = SIZE(api_columns, 2)
                        if (lpi_size(2) < api_columns(A_1, 0)) then
                            ALLOCATE(lpi_columns(lpi_size(1), lpi_size(2)))
                            lpi_columns = api_columns

                            DEALLOCATE(api_columns)

                            ALLOCATE(api_columns(lpi_size(1), 2 * lpi_size(2)))
                            api_columns(1:lpi_size(1), 1:lpi_size(2)) = lpi_columns(1:lpi_size(1), 1:lpi_size(2))

                            DEALLOCATE(lpi_columns)
                        end if


                    end if

                end do

            end do

        end do

        ! COUNT NON ZERO ELEMENTS
        result = SUM(api_occ(1: ai_nR))

        sll_count_non_zero_elts = result
    end function sll_count_non_zero_elts

subroutine init_sparsematrix( ai_nel,      &
                              api_LM_1,    &
                              ai_nen_1,    &
                              api_LM_2,    &
                              ai_nen_2,    &
                              api_columns, &
                              api_occ,     &
                              row_ptr,     &
                              col_ind,     &
                              num_rows     )

! _1 FOR ROWS
! _2 FOR COLUMNS
sll_int32, dimension(:,:), intent(in)  :: api_LM_1
sll_int32, dimension(:,:), intent(in)  :: api_LM_2
sll_int32                              :: ai_nel
sll_int32                              :: ai_nen_1
sll_int32                              :: ai_nen_2
sll_int32, dimension(:,:), pointer     :: api_columns
sll_int32, dimension(:),   pointer     :: api_occ
sll_int32, dimension(:),   intent(out) :: row_ptr
sll_int32, dimension(:),   intent(out) :: col_ind
sll_int32, intent(in) :: num_rows

sll_int32 :: e, b_1, A_1, i, size
sll_int32 :: err, flag
sll_int32, dimension(:), pointer :: lpr_tmp

! INITIALIZING ia
row_ptr(1) = 1
do i = 1, num_rows
  row_ptr(i+1) = row_ptr(1) + sum(api_occ(1: i))
end do

! INITIALIZING ja
do e = 1, ai_nel
  do b_1 = 1, ai_nen_1

    A_1 = api_LM_1(b_1, e)
    if (A_1 == 0) cycle
    if (api_columns(A_1, 0) == 0) cycle
    size = api_columns(A_1, 0)
    allocate ( lpr_tmp(size), stat = err)
    if (err .ne. 0) flag = 10
    lpr_tmp(1: size) = api_columns(A_1, 1: size)
    call sll_s_qsortc(lpr_tmp)
    do i = 1, size
      col_ind(row_ptr(A_1) + i - 1) = int ( lpr_tmp(i))
    end do
    api_columns(A_1, 0) = 0
    deallocate ( lpr_tmp)

  end do
end do

end subroutine init_sparsematrix

subroutine sll_s_solve_csr_matrix_perper(mat, apr_B, apr_U,Masse_tot)

  type(sll_t_csr_matrix)            :: mat
  sll_real64, dimension(:)          :: apr_U
  sll_real64, dimension(:)          :: apr_B
  sll_real64, dimension(:), pointer :: Masse_tot

  sll_int32  :: sys
  sys = 0
    
end subroutine sll_s_solve_csr_matrix_perper

end module sll_m_sparse_matrix
