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
  use sll_m_umfpack
  use sll_m_qsort_partition
  use iso_fortran_env, only: output_unit


  implicit none

 public :: &
    new_csr_matrix, &
    new_csr_matrix_with_constraint, &
    csr_add_one_constraint, &
    csr_todense, &
    delete_csr_matrix, &
    initialize_csr_matrix, &
    initialize_csr_matrix_with_constraint, &
    sll_add_to_csr_matrix, &
    sll_csr_matrix, &
    sll_factorize_csr_matrix, &
    sll_mult_csr_matrix_vector, &
    sll_solve_csr_matrix, &
    sll_solve_csr_matrix_perper, &
    sll_delete

  private

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!> @brief type for CSR format
type sll_csr_matrix

  sll_int32                         :: num_rows !< number of rows
  sll_int32                         :: num_cols !< number of columns
  sll_int32                         :: num_nz   !< number of non zero elements
  sll_int32,  dimension(:), pointer :: row_ptr
  sll_int32,  dimension(:), pointer :: col_ind
  sll_real64, dimension(:), pointer :: val
  sll_int32,  dimension(:), pointer :: opi_i
  sll_int32,  dimension(:), pointer :: Ai
  sll_int32,  dimension(:), pointer :: Ap

  integer(umf_void)                 :: umf_symbolic
  integer(umf_void)                 :: umf_numeric
  sll_real64, dimension(:), pointer :: umf_control

end type sll_csr_matrix

interface sll_delete
  module procedure delete_csr_matrix
end interface sll_delete

contains

subroutine delete_csr_matrix(csr_mat)

  type(sll_csr_matrix),pointer :: csr_mat

  nullify(csr_mat)
 ! SLL_DEALLOCATE_ARRAY(csr_mat%row_ptr,ierr)
 ! SLL_DEALLOCATE_ARRAY(csr_mat%col_ind,ierr)
 ! SLL_DEALLOCATE_ARRAY(csr_mat%val,ierr)
 ! SLL_DEALLOCATE_ARRAY(csr_mat%opi_i,ierr)
 
end subroutine delete_csr_matrix

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
function new_csr_matrix( &
  num_rows,              &
  num_cols,              &
  num_elements,          &
  local_to_global_row,   &
  num_local_dof_row,     &
  local_to_global_col,   &
  num_local_dof_col)     &
  result(mat)

  type(sll_csr_matrix), pointer         :: mat

  sll_int32,                 intent(in) :: num_rows
  sll_int32,                 intent(in) :: num_cols
  sll_int32,                 intent(in) :: num_elements
  sll_int32, dimension(:,:), intent(in) :: local_to_global_row
  sll_int32,                 intent(in) :: num_local_dof_row
  sll_int32, dimension(:,:), intent(in) :: local_to_global_col
  sll_int32,                 intent(in) :: num_local_dof_col

  sll_int32 :: ierr

  SLL_ALLOCATE(mat, ierr)

  call initialize_csr_matrix( &
    mat,                      &
    num_rows,                 &
    num_cols,                 &
    num_elements,             &
    local_to_global_row,      &
    num_local_dof_row,        &
    local_to_global_col,      &
    num_local_dof_col)
    
end function new_csr_matrix

!> @brief initialization of CSR matrix type
!> thanks to the global index of each local dof of each element
!> param[inout] mat : CSR matrix structure
!> param[in] num_rows :  number of rows
!> param[in] num_cols :  number of columns
!> param[in] num_element :  number of elements
!> param[in] local_to_global_row : local_to_global_row(\ell,i) gives the global 
!> row index of the matrix, for the element i and local degree of freedom \ell
!> param[in] num_local_dof_row : number of local degrees of freedom for the rows
!> param[in] local_to_global_col : local_to_global_col(\ell,i) gives the global 
!> column index of the matrix, for the element i and local degree of freedom \ell
!> param[in] num_local_dof_col : number of local degrees of freedom for the columns
subroutine initialize_csr_matrix( &
  mat,                            &
  num_rows,                       &
  num_cols,                       &
  num_elements,                   &
  local_to_global_row,            &
  num_local_dof_row,              &
  local_to_global_col,            &
  num_local_dof_col)

  type(sll_csr_matrix), intent(inout) :: mat

  sll_int32,                 intent(in) :: num_rows
  sll_int32,                 intent(in) :: num_cols
  sll_int32,                 intent(in) :: num_elements
  sll_int32, dimension(:,:), intent(in) :: local_to_global_row
  sll_int32,                 intent(in) :: num_local_dof_row
  sll_int32, dimension(:,:), intent(in) :: local_to_global_col
  sll_int32,                 intent(in) :: num_local_dof_col

  sll_int32                             :: num_nz
  sll_int32, dimension(:,:), pointer    :: lpi_columns
  sll_int32, dimension(:),   pointer    :: lpi_occ
  sll_int32                             :: li_COEF
  sll_int32                             :: ierr
  
  li_COEF = 10
  SLL_ALLOCATE(lpi_columns(num_rows, 0:li_COEF * num_local_dof_col),ierr)
  SLL_ALLOCATE(lpi_occ(num_rows + 1),ierr)

  lpi_columns(:,:) = 0
  lpi_occ(:) = 0
  ! COUNTING NON ZERO ELEMENTS
  
  num_nz = sll_count_non_zero_elts( &
    num_rows,                       &
    num_cols,                       &
    num_elements,                   &
    local_to_global_row,            &
    num_local_dof_row,              &
    local_to_global_col,            &
    num_local_dof_col,              &
    lpi_columns,                    &
    lpi_occ)

  mat%num_rows = num_rows
  mat%num_cols = num_cols
  mat%num_nz = num_nz
  SLL_ALLOCATE(mat%row_ptr(num_rows+1),ierr)
  SLL_ALLOCATE(mat%col_ind(num_nz),ierr)
  SLL_ALLOCATE(mat%val(num_nz),ierr)
  
  call sll_init_SparseMatrix( &
    mat,                      &
    num_elements,             &
    local_to_global_row,      &
    num_local_dof_row,        &
    local_to_global_col,      &
    num_local_dof_col,        &
    lpi_columns,              &
    lpi_occ)
  
  mat%val(:) = 0.0_f64
  
  print *,'#num_nz=',mat%num_nz
  print *,'#num_rows=',mat%num_rows
  print *,'#num_cols=',mat%num_cols
  
  SLL_ALLOCATE(mat%umf_control(umfpack_control),ierr)
  SLL_ALLOCATE(mat%Ai(num_nz),ierr)
  SLL_ALLOCATE(mat%Ap(num_rows+1),ierr)

  call umf4def(mat%umf_control)  ! get the default configuration
  ! modify the csr matrix to have indices starting at 0 as is the C convention
  
  SLL_DEALLOCATE_ARRAY(lpi_columns,ierr)
  SLL_DEALLOCATE_ARRAY(lpi_occ,ierr)
 
end subroutine initialize_csr_matrix

subroutine initialize_csr_matrix_classic( mat,                 &
                                          num_rows,            &
                                          num_cols,            &
                                          num_elements,        &
                                          local_to_global_row, &
                                          num_local_dof_row,   &
                                          local_to_global_col, & 
                                          num_local_dof_col)

  type(sll_csr_matrix),      intent(inout) :: mat
  sll_int32,                 intent(in)    :: num_rows
  sll_int32,                 intent(in)    :: num_cols
  sll_int32,                 intent(in)    :: num_elements
  sll_int32, dimension(:,:), intent(in)    :: local_to_global_row
  sll_int32,                 intent(in)    :: num_local_dof_row
  sll_int32, dimension(:,:), intent(in)    :: local_to_global_col
  sll_int32,                 intent(in)    :: num_local_dof_col
  
  sll_int32                                :: num_nz
  sll_int32                                :: ierr
  sll_int32,  dimension(:,:), allocatable  :: lpi_col
  sll_int32,  dimension(:),   allocatable  :: lpi_occ
  sll_int32                                :: COEF
  sll_int32                                :: elt
  sll_int32                                :: ii
  sll_int32                                :: jj
  sll_int32                                :: row
  sll_int32                                :: col
  sll_int32                                :: i
  sll_int32                                :: sz
  logical                                  :: ll_done
  
  print *,'#initialize_csr_matrix'
  COEF = 6
  
  SLL_ALLOCATE(lpi_col(num_rows, 0:COEF*num_local_dof_col),ierr)
  SLL_ALLOCATE(lpi_occ(num_rows+1),ierr)
  
  lpi_col(:,:) = 0
  lpi_occ(:) = 0
  
  do elt = 1, num_elements  !Loop over cells
    do ii = 1, num_local_dof_row
      row = local_to_global_row(ii, elt) !Row number in matrix
      if (row /= 0) then
        do jj = 1, num_local_dof_col
          col = local_to_global_col(jj, elt) !Column number in matrix
          if (col /= 0) then
            ll_done = .false.
            ! WE CHECK IF IT IS THE FIRST OCCURANCE OF THE COUPLE (row, col)
            if(lpi_col(row, 0)>COEF*num_local_dof_col)then
              print *,'Pb of size:',lpi_col(row, 0),COEF*num_local_dof_col
              stop
            endif
            do i = 1, lpi_col(row, 0)
              if (lpi_col(row, i) == col) then
                ll_done = .true.
                exit
              end if
            end do
            if (.not.ll_done) then
              lpi_occ(row)                  = lpi_occ(row) + 1
              lpi_col(row, 0)               = lpi_col(row, 0) + 1
              lpi_col(row, lpi_col(row, 0)) = col
            end if
          end if
        end do
      end if
    end do
  end do
  print *,'Size:',maxval(lpi_col(1:num_rows, 0)),COEF*num_local_dof_col
  
  ! COUNT NON ZERO ELEMENTS
  num_nz = SUM(lpi_occ(1:num_rows))
  
  mat%num_rows = num_rows
  mat%num_cols = num_cols
  mat%num_nz   = num_nz   
  
  print *,'#num_rows=',num_rows
  print *,'#num_nz=',num_nz
  
  SLL_ALLOCATE(mat%row_ptr(num_rows + 1),ierr)
  SLL_ALLOCATE(mat%col_ind(num_nz),ierr)
  SLL_ALLOCATE(mat%val(num_nz),ierr)
  
  mat%row_ptr(1) = 1
  
  do i = 1, mat%num_rows
    mat%row_ptr(i+1) = mat%row_ptr(1) + sum(lpi_occ(1:i))
  end do
  
  do elt = 1, num_elements
    do ii = 1, num_local_dof_row
      row = local_to_global_row(ii, elt)
      if (row /= 0) then
        if (lpi_col(row,0) /= 0) then
          sz = lpi_col(row, 0)
          call QsortC(lpi_col(row,1:sz))
          do i = 1, sz
            mat%col_ind(mat%row_ptr(row)+i-1) = lpi_col(row,i)
          end do
          lpi_col(row, 0) = 0
        end if
      end if
    end do
  end do
  
  mat%val(:) = 0.0_f64
  SLL_DEALLOCATE_ARRAY(lpi_col,ierr)
  SLL_DEALLOCATE_ARRAY(lpi_occ,ierr)

end subroutine initialize_csr_matrix_classic

subroutine initialize_csr_matrix_with_constraint( mat, mat_a)

  type(sll_csr_matrix), intent(inout) :: mat
  type(sll_csr_matrix), intent(in) :: mat_a
  sll_int32 :: ierr
  
  mat%num_nz = mat_a%num_nz + 2*mat_a%num_rows       
  print*,'num_nz mat, num_nz mat_tot', mat_a%num_nz,mat%num_nz 
  mat%num_rows = mat_a%num_rows  +  1
  print*,'num_rows mat, num_rows mat_tot',mat_a%num_rows , mat%num_rows
  mat%num_cols = mat_a%num_cols  +  1
  print*,'num_cols mat, num_cols mat_tot',mat_a%num_cols , mat%num_cols 
  
  SLL_ALLOCATE(mat%row_ptr(mat%num_rows+1),ierr)
  SLL_ALLOCATE(mat%col_ind(mat%num_nz),ierr)
  SLL_ALLOCATE(mat%val(mat%num_nz),ierr)
  mat%val(:) = 0.0_f64
  
  SLL_ALLOCATE(mat%umf_control(umfpack_control),ierr)
  SLL_ALLOCATE(mat%Ai(mat%num_nz),ierr)
  SLL_ALLOCATE(mat%Ap(mat%num_rows+1),ierr)
  ! get the default configuration
  call umf4def(mat%umf_control)  

end subroutine initialize_csr_matrix_with_constraint

function new_csr_matrix_with_constraint(mat_a) result(mat)

  type(sll_csr_matrix), pointer :: mat
  type(sll_csr_matrix) :: mat_a
  sll_int32 :: ierr
  SLL_ALLOCATE(mat, ierr)
  call initialize_csr_matrix_with_constraint( mat, mat_a)

end function new_csr_matrix_with_constraint
   
subroutine csr_add_one_constraint( &
  ia_in,                           &
  ja_in,                           &
  a_in,                            &
  num_rows_in,                     &
  num_nz_in,                       &
  constraint_vec,                  &
  ia_out,                          &
  ja_out,                          &
  a_out)

  integer, dimension(:), intent(in)  :: ia_in  
  integer, dimension(:), intent(in)  :: ja_in  
  real(8), dimension(:), intent(in)  :: a_in
  integer,               intent(in)  :: num_rows_in
  integer,               intent(in)  :: num_nz_in
  real(8), dimension(:), intent(in)  :: constraint_vec
  integer, dimension(:), intent(out) :: ia_out
  integer, dimension(:), intent(out) :: ja_out
  real(8), dimension(:), intent(out) :: a_out
  integer :: num_rows_out
  integer :: num_nz_out
  integer :: i
  integer :: s
  integer :: k
  
  num_rows_out = num_rows_in+1
  num_nz_out = num_nz_in+2*num_rows_in
  
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
  
end subroutine csr_add_one_constraint

subroutine sll_factorize_csr_matrix(mat)
  type(sll_csr_matrix), intent(inout) :: mat
  sll_real64, dimension(umfpack_info) :: info
  
  mat%Ap = mat%row_ptr(:) - 1
  mat%Ai = mat%col_ind(:) - 1

  ! pre-order and symbolic analysis
  call umf4sym( &
    mat%num_rows, &
    mat%num_cols, &
    mat%Ap, &
    mat%Ai, &
    mat%val, &
    mat%umf_symbolic, &
    mat%umf_control, &
    info)

  if (info(1) .lt. 0) then
    print *, '#Error occurred in umf4sym: ', info(1)
    stop
  endif

  ! numeric factorization
  call umf4num( &
    mat%Ap, &
    mat%Ai, &
    mat%val, &
    mat%umf_symbolic, &
    mat%umf_numeric, &
    mat%umf_control, &
    info)
    
  if (info(1) .lt. 0) then
    print *, '#Error occurred in umf4num: ', info(1)
    stop
  endif
  
end subroutine sll_factorize_csr_matrix

subroutine sll_mult_csr_matrix_vector(mat, input, output)
  type(sll_csr_matrix) :: mat
  sll_real64, dimension(:), intent(in) :: input
  sll_real64, dimension(:), intent(out) :: output
  !local var
  sll_int32 :: li_i
  sll_int32 :: li_k_1
  sll_int32 :: li_k_2

  do li_i = 1, mat % num_rows

    li_k_1 = mat % row_ptr(li_i)
    li_k_2 = mat % row_ptr(li_i + 1) - 1
    output(li_i) = &
      DOT_PRODUCT(mat % val(li_k_1: li_k_2), input(mat % col_ind(li_k_1: li_k_2)))
          
  end do

end subroutine sll_mult_csr_matrix_vector

subroutine sll_add_to_csr_matrix(mat, val, ai_A, ai_Aprime)
  type(sll_csr_matrix) :: mat
  sll_real64, intent(in) :: val
  sll_int32, intent(in) :: ai_A
  sll_int32, intent(in) :: ai_Aprime
  !local var
  sll_int32 :: li_j
  sll_int32 :: li_k

  ! THE CURRENT LINE IS self%row_ptr(ai_A)
  do li_k = mat % row_ptr(ai_A), mat % row_ptr(ai_A + 1) - 1
    li_j = mat % col_ind(li_k)
    if (li_j == ai_Aprime) then
      mat % val(li_k) = mat % val(li_k) + val 
      exit
    end if
  end do

end subroutine sll_add_to_csr_matrix

subroutine sll_solve_csr_matrix(mat, apr_B, apr_U)
  type(sll_csr_matrix) :: mat
  sll_real64, dimension(:) :: apr_U
  sll_real64, dimension(:) :: apr_B
  !local var
  sll_int32  :: sys
  sll_real64, dimension(umfpack_info) :: info
  sys = 0
  call umf4sol(sys,apr_U,apr_B,mat%umf_numeric,mat%umf_control,info)
end subroutine sll_solve_csr_matrix

integer function sll_count_non_zero_elts( &
  ai_nR, ai_nC, ai_nel, api_LM_1, ai_nen_1, api_LM_2, ai_nen_2, api_columns, api_occ)
  ! _1 FOR ROWS
  ! _2 FOR COLUMNS
  integer :: ai_nR, ai_nC
  integer, dimension(:,:), intent(in) :: api_LM_1, api_LM_2
  integer :: ai_nel, ai_nen_1, ai_nen_2
  integer, dimension(:,:), pointer :: api_columns
  integer, dimension(:), pointer :: api_occ
  !local var
  integer :: li_e, li_b_1, li_A_1, li_b_2, li_A_2, li_i
  integer :: li_result
  integer, dimension(2) :: lpi_size
  logical :: ll_done
  integer, dimension(:,:), pointer :: lpi_columns

  ! WE FIRST COMPUTE, FOR EACH ROW, THE NUMBER OF COLUMNS THAT WILL BE USED
  do li_e = 1, ai_nel
      do li_b_1 = 1, ai_nen_1
          li_A_1 = api_LM_1(li_b_1, li_e)
          if (li_A_1 == 0) then
              cycle
          end if

          do li_b_2 = 1, ai_nen_2

              li_A_2 = api_LM_2(li_b_2, li_e)
              if (li_A_2 == 0) then
                  cycle
              end if

              ll_done = .false.
              ! WE CHECK IF IT IS THE FIRST OCCURANCE OF THE COUPLE (li_A_1, li_A_2)
              do li_i = 1, api_columns(li_A_1, 0)

                  if (api_columns(li_A_1, li_i) /= li_A_2) then
                      cycle
                  end if

                  ll_done = .true.
                  exit

              end do

              if (.not.ll_done) then

                  api_occ(li_A_1) = api_occ(li_A_1) + 1

                  ! li_A_1 IS THE ROW NUM, li_A_2 THE COLUMN NUM
                  ! INITIALIZATION OF THE SPARSE MATRIX
                  api_columns(li_A_1, 0) = api_columns(li_A_1, 0) + 1
          
                  ! resizing the array
                  lpi_size(1) = SIZE(api_columns, 1)
                  lpi_size(2) = SIZE(api_columns, 2)
                      
                  do while (SIZE(api_columns, 2) < api_columns(li_A_1, 0)+1)

                  !if (lpi_size(2) < api_columns(li_A_1, 0)+1) then
                      ALLOCATE(lpi_columns(lpi_size(1), 0:lpi_size(2)-1))
                      lpi_columns = api_columns

                      DEALLOCATE(api_columns)

                      ALLOCATE(api_columns(lpi_size(1), 0:2 * lpi_size(2)))
                      api_columns(1:lpi_size(1), 0:lpi_size(2)-1) = lpi_columns(1:lpi_size(1), 0:lpi_size(2)-1)

                      DEALLOCATE(lpi_columns)
                  enddo
                  !end if
                  !print *,'api_columns(li_A_1,0)=',api_columns(li_A_1,0)
                  !print *,'lpi_size(2)=',lpi_size(2)
                  flush( output_unit )
                  api_columns(li_A_1, api_columns(li_A_1, 0)) = li_A_2
              end if
          end do
      end do
  end do

  ! COUNT NON ZERO ELEMENTS
  li_result = SUM(api_occ(1: ai_nR))

  sll_count_non_zero_elts = li_result

end function sll_count_non_zero_elts

subroutine sll_init_SparseMatrix(self, ai_nel, api_LM_1, ai_nen_1, api_LM_2, &
  ai_nen_2, api_columns, api_occ)
    ! _1 FOR ROWS
    ! _2 FOR COLUMNS
    implicit none
    type(sll_csr_matrix) :: self
    integer, dimension(:,:), intent(in) :: api_LM_1
    integer, dimension(:,:), intent(in) :: api_LM_2
    integer :: ai_nel, ai_nen_1, ai_nen_2
    integer, dimension(:,:), pointer :: api_columns
    integer, dimension(:), pointer :: api_occ
    !local var
    integer :: li_e, li_b_1, li_A_1, li_i, li_size
    integer :: li_err, li_flag
    sll_int32, dimension(:), pointer :: lpr_tmp

    ! INITIALIZING ia
    self % row_ptr(1) = 1

    do li_i = 1, self % num_rows

        self % row_ptr(li_i + 1) = self % row_ptr(1) + SUM(api_occ(1: li_i))

    end do

    ! INITIALIZING ja
    do li_e = 1, ai_nel

        do li_b_1 = 1, ai_nen_1

            li_A_1 = api_LM_1(li_b_1, li_e)

            if (li_A_1 == 0) then
                cycle
            end if

            if (api_columns(li_A_1, 0) == 0) then
                cycle
            end if

            li_size = api_columns(li_A_1, 0)

            allocate ( lpr_tmp(li_size), stat = li_err)
            if (li_err .ne. 0) li_flag = 10

            lpr_tmp(1: li_size) = api_columns(li_A_1, 1: li_size)

            call QsortC(lpr_tmp)

            do li_i = 1, li_size

                self % col_ind(self % row_ptr(li_A_1) + li_i - 1) =  ( lpr_tmp(li_i))

            end do

            api_columns(li_A_1, 0) = 0
            deallocate ( lpr_tmp)

        end do

    end do

end subroutine sll_init_SparseMatrix

subroutine sll_solve_csr_matrix_perper(mat, apr_B, apr_U,Masse_tot)

  type(sll_csr_matrix)                :: mat
  sll_real64, dimension(:)            :: apr_U
  sll_real64, dimension(:)            :: apr_B
  sll_real64, dimension(:), pointer   :: Masse_tot
  sll_real64, dimension(umfpack_info) :: info
  sll_int32  :: sys

  sys = 0
  call umf4sol(sys,apr_U,apr_B,mat%umf_numeric,mat%umf_control,info)
  
end subroutine sll_solve_csr_matrix_perper

!> @brief
!> Test function to initialize a CSR matrix
!> @details
!> Fill a matrix in CSR format corresponding to a constant coefficient
!> five-point stencil on a square grid
subroutine uni2d(this,f)

type(sll_csr_matrix) :: this
sll_real64           :: f(:)
sll_real64, pointer  :: a(:)
sll_int32            :: m
sll_int32, pointer   :: ia(:),ja(:)
integer              :: k,l,i,j

real (kind(0d0)), parameter :: zero =  0.0d0
real (kind(0d0)), parameter :: cx   = -1.0d0
real (kind(0d0)), parameter :: cy   = -1.0d0
real (kind(0d0)), parameter :: cd   = +4.0d0

a  => this%val
ia => this%row_ptr
ja => this%col_ind

m = this%num_rows

k=0
l=0
ia(1)=1
do i=1,m
  do j=1,m
    k=k+1
    l=l+1
    a(l)=cd
    ja(l)=k
    f(k)=zero
    if(j < m) then
       l=l+1
       a(l)=cx
       ja(l)=k+1
      else
       f(k)=f(k)-cx
    end if
    if(i < m) then
       l=l+1
       a(l)=cy
       ja(l)=k+m
      else
       f(k)=f(k)-cy
    end if
    if(j > 1) then
       l=l+1
       a(l)=cx
       ja(l)=k-1
      else
       f(k)=f(k)-cx
    end if
    if(i >  1) then
       l=l+1
       a(l)=cy
       ja(l)=k-m
      else
       f(k)=f(k)-cy
    end if
    ia(k+1)=l+1
  end do
end do

this%num_nz = l

return
end subroutine uni2d

end module sll_m_sparse_matrix
