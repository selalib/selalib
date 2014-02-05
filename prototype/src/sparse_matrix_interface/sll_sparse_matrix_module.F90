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

!> @file sll_sparse_matrix_module.F90
!> @namespace sll_sparse_matrix
!> @brief sparse matrix


module sll_sparse_matrix_module
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"


  !> @brief type for CSR format
  type sll_csr_matrix
    sll_int32 :: num_rows !< number of rows
    sll_int32 :: num_cols !< number of columns
    sll_int32 :: num_elements !< number of non zero elements
    !integer, dimension(:), pointer :: opi_ia
    !integer, dimension(:), pointer :: opi_ja
    !real(wp), dimension(:), pointer :: opr_a
        !................
    !logical :: ol_use_mm_format
    !integer, dimension(:), pointer :: opi_i
        !................
  end type sll_csr_matrix


contains


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
    num_rows, &
    num_cols, &
    num_elements, &
    local_to_global_row, &
    num_local_dof_row, &
    local_to_global_col, &
    num_local_dof_col ) &
    result(mat)
    type(sll_csr_matrix), pointer :: mat
    sll_int32, intent(in) :: num_rows
    sll_int32, intent(in) :: num_cols
    sll_int32, intent(in) :: num_elements
    sll_int32, dimension(:,:), intent(in) :: local_to_global_row
    sll_int32, intent(in) :: num_local_dof_row
    sll_int32, dimension(:,:), intent(in) :: local_to_global_col
    sll_int32, intent(in) :: num_local_dof_col
    sll_int32 :: ierr
    SLL_ALLOCATE(mat, ierr)
    call initialize_csr_matrix( &
      mat, &
      num_rows, &
      num_cols, &
      num_elements, &
      local_to_global_row, &
      num_local_dof_row, &
      local_to_global_col, &
      num_local_dof_col )
      
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
    mat, &
    num_rows, &
    num_cols, &
    num_elements, &
    local_to_global_row, &
    num_local_dof_row, &
    local_to_global_col, &
    num_local_dof_col )
    type(sll_csr_matrix), intent(inout) :: mat
    sll_int32, intent(in) :: num_rows
    sll_int32, intent(in) :: num_cols
    sll_int32, intent(in) :: num_elements
    sll_int32, dimension(:,:), intent(in) :: local_to_global_row
    sll_int32, intent(in) :: num_local_dof_row
    sll_int32, dimension(:,:), intent(in) :: local_to_global_col
    sll_int32, intent(in) :: num_local_dof_col
    !local variables
    sll_int32 :: li_err, li_flag
    sll_int32 :: li_nnz
    sll_int32, dimension(:,:), pointer :: lpi_columns
    sll_int32, dimension(:), pointer :: lpi_occ
    sll_int32 :: li_COEF
    sll_int32 :: ierr
    
    li_COEF = 10
    SLL_ALLOCATE(lpi_columns(num_rows, 0:li_COEF * num_local_dof_col),ierr)
    SLL_ALLOCATE(lpi_occ(num_rows + 1),ierr)
    lpi_columns(:,:) = 0
    lpi_occ(:) = 0
    ! COUNTING NON ZERO ELEMENTS
    li_nnz = sll_count_non_zero_elts( &
      num_rows, &
      num_cols, &
      num_elements, &
      local_to_global_row, &
      num_local_dof_row, &
      local_to_global_col, &
      num_local_dof_col, &
      lpi_columns, &
      lpi_occ)


 
    
!    subroutine create_SparseMatrix_with_LM(this, ai_nR, ai_nC, ai_nel, api_LM_1, 
!  ai_nen_1, api_LM_2, ai_nen_2, ai_COEF)
!        implicit none
!        !> param[inout] this : CSR MATRIX STRUCTURE
!        type(csr_matrix) :: this
!        !> param[in] ai_nC : NUMBER OF COLUMNS, IT IS THE DIMENSION OF THE 1st SPACE
!        integer :: ai_nC
!        !> param[in] ai_nR : NUMBER OF ROWS, IT IS THE DIMENSION OF THE 2nd SPACE
!        integer :: ai_nR
!        !> param[in] ai_nel : NUMBER OF ELEMENTS (WITH NON ZERO MEASURE) IN THE PATCH
!        integer :: ai_nel
!        !> param[in] api_LM_1 : LM ARRAY FOR ROWS
!        integer, dimension(:,:), pointer :: api_LM_1
!        !> param[in] api_LM_1 : LM ARRAY FOR COLUMNS
!        integer, dimension(:,:), pointer :: api_LM_2
!        !> param[in] ai_nen_1 : NUMBER OF NON VANISHING FUNCTIONS PER ELEMENT, IN THE 1st SPACE
!        integer :: ai_nen_1
!        !> param[in] ai_nen_2 : NUMBER OF NON VANISHING FUNCTIONS PER ELEMENT, IN THE 2nd SPACE
!        integer :: ai_nen_2
!        !>
!        integer, optional :: ai_COEF
!        !local var
!        integer :: li_err, li_flag
!        integer :: li_nnz
!        integer, dimension(:,:), pointer :: lpi_columns
!        integer, dimension(:), pointer :: lpi_occ
!        integer :: li_COEF
!
!        li_COEF = 10
!        if (present(ai_COEF)) then
!            li_COEF = ai_COEF
!        end if
!
!        allocate(lpi_columns(ai_nR, 0:li_COEF * ai_nen_2))
!
!        allocate(lpi_occ(ai_nR + 1))
!        lpi_columns(:,:) = 0
!        lpi_occ(:) = 0
!
!        ! COUNTING NON ZERO ELEMENTS
!
!        li_nnz = count_non_zero_elts(ai_nR, ai_nC, ai_nel, api_LM_1, ai_nen_1, api_LM_2, ai_nen_2, lpi_columns, lpi_occ)
!
!
!        this % ol_use_mm_format = .FALSE.
!
!        this % oi_nR = ai_nR
!        this % oi_nC = ai_nC
!        this % oi_nel = li_nnz
!
!
!        allocate(this % opi_ia(this % oi_nR + 1))
!
!
!        allocate(this % opi_ja(this % oi_nel))
!
!
!        allocate(this % opr_a(this % oi_nel))
!
!
!        call init_SparseMatrix(this, ai_nel, api_LM_1, ai_nen_1, api_LM_2, ai_nen_2, lpi_columns, lpi_occ)
!
!        this % opr_a(:) = 0.0_wp
!
!        deallocate(lpi_columns)
!        deallocate(lpi_occ)
!    end subroutine create_SparseMatrix_with_LM





  end subroutine initialize_csr_matrix


    integer function sll_count_non_zero_elts( &
      ai_nR, ai_nC, ai_nel, api_LM_1, ai_nen_1, api_LM_2, ai_nen_2, api_columns, api_occ)
        ! _1 FOR ROWS
        ! _2 FOR COLUMNS
        implicit none
        integer :: ai_nR, ai_nC
        integer, dimension(:,:), intent(in) :: api_LM_1, api_LM_2
        integer :: ai_nel, ai_nen_1, ai_nen_2
        integer, dimension(:,:), pointer :: api_columns
        integer, dimension(:), pointer :: api_occ
        !local var
        integer :: li_e, li_b_1, li_A_1, li_b_2, li_A_2, li_i
        integer :: li_err, li_flag
        real(f64), dimension(:), pointer :: lpr_tmp
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
                        api_columns(li_A_1, api_columns(li_A_1, 0)) = li_A_2

                        ! resizing the array
                        lpi_size(1) = SIZE(api_columns, 1)
                        lpi_size(2) = SIZE(api_columns, 2)
                        if (lpi_size(2) < api_columns(li_A_1, 0)) then
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
        li_result = SUM(api_occ(1: ai_nR))

        sll_count_non_zero_elts = li_result
    end function sll_count_non_zero_elts

!    subroutine sll_init_SparseMatrix(self, ai_nel, api_LM_1, ai_nen_1, api_LM_2, &
!      ai_nen_2, api_columns, api_occ)
!        ! _1 FOR ROWS
!        ! _2 FOR COLUMNS
!        implicit none
!        type(sll_csr_matrix) :: self
!        integer, dimension(:,:), pointer :: api_LM_1, api_LM_2
!        integer :: ai_nel, ai_nen_1, ai_nen_2
!        integer, dimension(:,:), pointer :: api_columns
!        integer, dimension(:), pointer :: api_occ
!        !local var
!        integer :: li_e, li_b_1, li_A_1, li_b_2, li_A_2, li_index, li_i, li_size
!        integer :: li_err, li_flag
!        real(f64), dimension(:), pointer :: lpr_tmp
!
!        ! INITIALIZING ia
!        self % opi_ia(1) = 1
!
!        do li_i = 1, self % oi_nR
!
!            self % opi_ia(li_i + 1) = self % opi_ia(1) + SUM(api_occ(1: li_i))
!
!        end do
!
!        ! INITIALIZING ja
!        do li_e = 1, ai_nel
!
!            do li_b_1 = 1, ai_nen_1
!
!                li_A_1 = api_LM_1(li_b_1, li_e)
!
!                if (li_A_1 == 0) then
!                    cycle
!                end if
!
!                if (api_columns(li_A_1, 0) == 0) then
!                    cycle
!                end if
!
!                li_size = api_columns(li_A_1, 0)
!
!                allocate ( lpr_tmp(li_size), stat = li_err)
!                if (li_err .ne. 0) li_flag = 10
!
!                lpr_tmp(1: li_size) = real( api_columns(li_A_1, 1: li_size))
!
!                call QsortC(lpr_tmp)
!
!                do li_i = 1, li_size
!
!                    self % opi_ja(self % opi_ia(li_A_1) + li_i - 1) = int ( lpr_tmp(li_i))
!
!                end do
!
!                api_columns(li_A_1, 0) = 0
!                deallocate ( lpr_tmp)
!
!            end do
!
!        end do
!
!    end subroutine init_SparseMatrix











!    call create_CSR( &
!        es%csr_mat, &
!        solution_size, &
!        solution_size, &
!        num_cells_eta1*num_cells_eta2, &
!        es%local_to_global_spline_indices, &
!        es%total_num_splines_loc, &
!        es%local_to_global_spline_indices, &
!        es%total_num_splines_loc )




end module sll_sparse_matrix_module