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
use mod_umfpack


  !> @brief type for CSR format
  type sll_csr_matrix
    sll_int32 :: num_rows !< number of rows
    sll_int32 :: num_cols !< number of columns
    sll_int32 :: num_nz !< number of non zero elements
    sll_int32, dimension(:), pointer :: opi_ia
    sll_int32, dimension(:), pointer :: opi_ja
    sll_real64, dimension(:), pointer :: opr_a
        !................
    !logical :: ol_use_mm_format
    sll_int32, dimension(:), pointer :: opi_i
        !................
    ! work arrays for Umfpack
    sll_int32, dimension(:), pointer :: Ai, Ap
    integer(umf_void) :: umf_symbolic
    integer(umf_void) :: umf_numeric
    sll_real64, dimension(:), pointer :: umf_control
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
    sll_int32 :: num_nz
    sll_int32, dimension(:,:), pointer :: lpi_columns
    sll_int32, dimension(:), pointer :: lpi_occ
    sll_int32 :: li_COEF
    sll_int32 :: ierr
    !sll_real64, dimension(umfpack_info) :: info
    
    
    li_COEF = 10
    SLL_ALLOCATE(lpi_columns(num_rows, 0:li_COEF * num_local_dof_col),ierr)
    SLL_ALLOCATE(lpi_occ(num_rows + 1),ierr)
    lpi_columns(:,:) = 0
    lpi_occ(:) = 0
    ! COUNTING NON ZERO ELEMENTS
    num_nz = sll_count_non_zero_elts( &
      num_rows, &
      num_cols, &
      num_elements, &
      local_to_global_row, &
      num_local_dof_row, &
      local_to_global_col, &
      num_local_dof_col, &
      lpi_columns, &
      lpi_occ)

    mat%num_rows = num_rows
    mat%num_cols = num_cols
    mat%num_nz = num_nz
    SLL_ALLOCATE(mat%opi_ia(num_rows + 1),ierr)
    SLL_ALLOCATE(mat%opi_ja(num_nz),ierr)
    SLL_ALLOCATE(mat%opr_a(num_nz),ierr)
    
    call sll_init_SparseMatrix( &
      mat, &
      num_elements, &
      local_to_global_row, &
      num_local_dof_row, &
      local_to_global_col, &
      num_local_dof_col, &
      lpi_columns, &
      lpi_occ)
    
    mat%opr_a(:) = 0.0_f64
    
    !print *,umfpack_control
    !stop
    
    SLL_ALLOCATE(mat%umf_control(umfpack_control),ierr)
    SLL_ALLOCATE(mat%Ai(num_nz),ierr)
    SLL_ALLOCATE(mat%Ap(num_rows+1),ierr)

    call umf4def(mat%umf_control)  ! get the default configuration
!    es%umf_control(umfpack_prl) = real( 2 , umf_dp ) ! change verbosity
!    call umf4pcon(es%umf_control) ! update the umfpack configuration
    ! modify the csr matrix to have indices starting at 0 as is the C convention
    
    
    


    
    
    SLL_DEALLOCATE_ARRAY(lpi_columns,ierr)
    SLL_DEALLOCATE_ARRAY(lpi_occ,ierr)
    
   
  end subroutine initialize_csr_matrix
  
  subroutine sll_factorize_csr_matrix(mat)
    type(sll_csr_matrix), intent(inout) :: mat
    sll_real64, dimension(umfpack_info) :: info

    mat%Ap = mat%opi_ia(:) - 1
    mat%Ai = mat%opi_ja(:) - 1


    ! pre-order and symbolic analysis
    call umf4sym( &
      mat%num_rows, &
      mat%num_cols, &
      mat%Ap, &
      mat%Ai, &
      mat%opr_a, &
      mat%umf_symbolic, &
      mat%umf_control, &
      info)

    
    ! numeric factorization
    call umf4num( &
      mat%Ap, &
      mat%Ai, &
      mat%opr_a, &
      mat%umf_symbolic, &
      mat%umf_numeric, &
      mat%umf_control, &
      info)

    
  end subroutine sll_factorize_csr_matrix
  

  subroutine sll_mult_csr_matrix_vector(mat, input, output)
    implicit none
    type(sll_csr_matrix) :: mat
    sll_real64, dimension(:), intent(in) :: input
    sll_real64, dimension(:), intent(out) :: output
    !local var
    sll_int32 :: li_i
    sll_int32 :: li_k_1
    sll_int32 :: li_k_2

    do li_i = 1, mat % num_rows

      li_k_1 = mat % opi_ia(li_i)
      li_k_2 = mat % opi_ia(li_i + 1) - 1
      output(li_i) = &
        DOT_PRODUCT(mat % opr_a(li_k_1: li_k_2), input(mat % opi_ja(li_k_1: li_k_2)))
            
    end do

  end subroutine sll_mult_csr_matrix_vector



  subroutine sll_add_to_csr_matrix(mat, val, ai_A, ai_Aprime)
    implicit none
    type(sll_csr_matrix) :: mat
    sll_real64, intent(in) :: val
    sll_int32, intent(in) :: ai_A
    sll_int32, intent(in) :: ai_Aprime
    !local var
    sll_int32 :: li_j
    sll_int32 :: li_k


    ! THE CURRENT LINE IS self%opi_ia(ai_A)
    do li_k = mat % opi_ia(ai_A), mat % opi_ia(ai_A + 1) - 1
      li_j = mat % opi_ja(li_k)
      if (li_j == ai_Aprime) then
        mat % opr_a(li_k) = mat % opr_a(li_k) + val
        exit
      end if
    end do

  end subroutine sll_add_to_csr_matrix
 
  subroutine sll_sub_to_csr_matrix(mat, val, ai_A, ai_Aprime,Masse_tot)
    implicit none
    type(sll_csr_matrix) :: mat
    sll_real64, intent(in) :: val
    sll_int32, intent(in) :: ai_A
    sll_int32, intent(in) :: ai_Aprime
    sll_real64, dimension(:) :: Masse_tot
    !local var
    sll_int32 :: li_j
    sll_int32 :: li_k


    ! THE CURRENT LINE IS self%opi_ia(ai_A)
    do li_k = mat % opi_ia(ai_A), mat % opi_ia(ai_A + 1) - 1
      li_j = mat % opi_ja(li_k)
      if (li_j == ai_Aprime) then
        mat % opr_a(li_k) = mat % opr_a(li_k) + val - Masse_tot(li_k)
        exit
      end if
    end do

  end subroutine sll_sub_to_csr_matrix

  subroutine sll_solve_csr_matrix(mat, apr_B, apr_U)
    implicit none
    type(sll_csr_matrix) :: mat
    sll_real64, dimension(:) :: apr_U
    sll_real64, dimension(:) :: apr_B
    !local var
    sll_int32  :: sys
    sll_real64, dimension(umfpack_info) :: info
    sys = 0
    call umf4sol(sys,apr_U,apr_B,mat%umf_numeric,mat%umf_control,info)
    !print *,apr_U
    !stop
  end subroutine sll_solve_csr_matrix


!    sll_int32  :: sys
!    sll_real64, dimension(umfpack_info) :: info
!      call umf4sol(sys,apr_U,apr_B,es%umf_numeric,es%umf_control,info)
!       call Gradient_conj(&
!            csr_mat,&
!            apr_B,&
!            apr_U,&
!            ai_maxIter,&
!            ar_eps )




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
        integer :: li_e, li_b_1, li_A_1, li_b_2, li_A_2, li_index, li_i, li_size
        integer :: li_err, li_flag
        real(f64), dimension(:), pointer :: lpr_tmp

        ! INITIALIZING ia
        self % opi_ia(1) = 1

        do li_i = 1, self % num_rows

            self % opi_ia(li_i + 1) = self % opi_ia(1) + SUM(api_occ(1: li_i))

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

                lpr_tmp(1: li_size) = real( api_columns(li_A_1, 1: li_size))

                call QsortC(lpr_tmp)

                do li_i = 1, li_size

                    self % opi_ja(self % opi_ia(li_A_1) + li_i - 1) = int ( lpr_tmp(li_i))

                end do

                api_columns(li_A_1, 0) = 0
                deallocate ( lpr_tmp)

            end do

        end do

    end subroutine sll_init_SparseMatrix





recursive subroutine QsortC(A)
  real(f64), intent(in out), dimension(:) :: A
  integer :: iq

  if(size(A) > 1) then
     call Partition(A, iq)
     call QsortC(A(:iq-1))
     call QsortC(A(iq:))
  endif
end subroutine QsortC

subroutine Partition(A, marker)
  real(f64), intent(in out), dimension(:) :: A
  integer, intent(out) :: marker
  integer :: i, j
  real(f64) :: temp
  real(f64) :: x      ! pivot point
  x = A(1)
  i= 0
  j= size(A) + 1

  do
     j = j-1
     do
        if (A(j) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (A(i) >= x) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i) and A(j)
        temp = A(i)
        A(i) = A(j)
        A(j) = temp
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do

end subroutine Partition






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