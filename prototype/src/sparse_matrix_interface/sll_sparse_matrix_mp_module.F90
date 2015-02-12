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

!> @file sll_sparse_matrix_mp_module.F90
!> @namespace sll_sparse_matrix
!> @brief sparse matrix


module sll_sparse_matrix_mp_module
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_assert.h"
  
  use sll_sparse_matrix_module
  
  implicit none
   

contains


  !> @brief allocates the memory space for a new CSR type on the heap,
  !> initializes it with the given arguments and returns a pointer to the
  !> object.
  !> param[in] num_rows :  number of rows
  !> param[in] num_cols :  number of columns
  !> parma[in] num_patch:  number of patch  
  !> param[in] num_element :  num_element(p) number of elements non zero for each patch
  !> param[in] local_to_global_row : local_to_global_row(p,\ell,i) gives the global 
  !> row index of the matrix, for the element i and local degree of freedom \ell for patch p
  !> param[in] num_local_dof_row : num_local_dof_row(p) number of local degrees of freedom for the rows for each patch
  !> param[in] local_to_global_col : local_to_global_col(p,\ell,i) gives the global 
  !> column index of the matrix, for the element i and local degree of freedom \ell for patch p
  !> param[in] num_local_dof_col : num_local_dof_col(p) number of local degrees of freedom for the columns for each patch
  !> return a pointer to the newly allocated object.
  function new_csr_matrix_mp( &
    num_rows, &
    num_cols, &
    num_patch,&
    num_elements, &
    local_to_global_row, &
    num_local_dof_row, &
    local_to_global_col, &
    num_local_dof_col ) &
    result(mat)
    type(sll_csr_matrix), pointer :: mat
    sll_int32, intent(in) :: num_rows
    sll_int32, intent(in) :: num_cols
    sll_int32, intent(in) :: num_patch
    sll_int32, dimension(:), intent(in) :: num_elements
    sll_int32, dimension(:,:,:), intent(in) :: local_to_global_row
    sll_int32, dimension(:), intent(in)     :: num_local_dof_row
    sll_int32, dimension(:,:,:), intent(in) :: local_to_global_col
    sll_int32, dimension(:),intent(in)      :: num_local_dof_col
    sll_int32 :: ierr
    SLL_ALLOCATE(mat, ierr)
    call initialize_csr_matrix_mp( &
      mat, &
      num_rows, &
      num_cols, &
      num_patch,&
      num_elements, &
      local_to_global_row, &
      num_local_dof_row, &
      local_to_global_col, &
      num_local_dof_col )
      
  end function new_csr_matrix_mp

  !> @brief initialization of CSR matrix type
  !> thanks to the global index of each local dof of each element
  !> param[inout] mat : CSR matrix structure
  !> param[in] num_rows :  number of rows
  !> param[in] num_cols :  number of columns
  !> param[in] num_patch:  number of patchs
  !> param[in] num_element : num_element(p)  number of non zero elements in each patch 
  !> param[in] local_to_global_row : local_to_global_row(p,\ell,i) gives the global 
  !> row index of the matrix, for the element i and local degree of freedom \ell for patch p
  !> param[in] num_local_dof_row : num_local_dof_row(p) number of local degrees of freedom for the rows for each patch
  !> param[in] local_to_global_col : local_to_global_col(p,\ell,i) gives the global 
  !> column index of the matrix, for the element i and local degree of freedom \ell for patch p
  !> param[in] num_local_dof_col : num_local_dof_col(p) number of local degrees of freedom for the columns for each patch
  subroutine initialize_csr_matrix_mp( &
    mat, &
    num_rows, &
    num_cols, &
    num_patch,&
    num_elements, &
    local_to_global_row, &
    num_local_dof_row, &
    local_to_global_col, &
    num_local_dof_col)
    type(sll_csr_matrix), intent(inout) :: mat
    sll_int32, intent(in) :: num_rows
    sll_int32, intent(in) :: num_cols
    sll_int32, intent(in) :: num_patch
    sll_int32, dimension(:), intent(in) :: num_elements
    sll_int32, dimension(:,:,:), intent(in) :: local_to_global_row
    sll_int32, dimension(:), intent(in)     :: num_local_dof_row
    sll_int32, dimension(:,:,:), intent(in) :: local_to_global_col
    sll_int32, dimension(:), intent(in)     :: num_local_dof_col
    !local variables
    sll_int32 :: num_nz
    sll_int32, dimension(:,:), pointer :: lpi_columns
    sll_int32, dimension(:), pointer :: lpi_occ
    sll_int32 :: li_COEF
    sll_int32 :: ierr
    sll_int32 :: max_nen_C
    !print *,'#num_rows=',num_rows
    !print *,'#num_nz=',num_nz

    
    li_COEF = 20
    ! for we take the maximun of number of element basis in all patches 
    max_nen_C= MAXVAL(num_local_dof_col(:)) 
 
    SLL_ALLOCATE(lpi_columns(num_rows, 0:li_COEF * max_nen_C),ierr)
    SLL_ALLOCATE(lpi_occ(num_rows + 1),ierr)
    lpi_columns(:,:) = 0
    lpi_occ(:) = 0

    ! COUNTING NON ZERO ELEMENTS
    num_nz = sll_count_non_zero_elts_mp( &
      num_rows, &
      num_cols, &
      num_patch,&
      num_elements, &
      local_to_global_col, &
      num_local_dof_col, &
      local_to_global_row, &
      num_local_dof_row, &
      lpi_columns, &
      lpi_occ)

    mat%num_rows = num_rows
    mat%num_cols = num_cols
    mat%num_nz = num_nz

    SLL_ALLOCATE(mat%opi_ia(num_rows + 1),ierr)
    SLL_ALLOCATE(mat%opi_ja(num_nz),ierr)
    SLL_ALLOCATE(mat%opr_a(num_nz),ierr)
    
    print *,'#num_rows=',num_rows
    print *,'#num_nz=',num_nz
    
    call sll_init_SparseMatrix_mp( &
      mat, &
      num_rows,&
      num_cols,&
      num_patch,&
      num_elements, &
      local_to_global_col, &
      num_local_dof_col, &
      local_to_global_row, &
      num_local_dof_row, &
      lpi_columns, &
      lpi_occ)
    
    mat%opr_a(:) = 0.0_f64
    SLL_DEALLOCATE_ARRAY(lpi_columns,ierr)
    SLL_DEALLOCATE_ARRAY(lpi_occ,ierr)

  end subroutine initialize_csr_matrix_mp
  
  integer function sll_count_non_zero_elts_mp( &
       ai_nR,&
       ai_nC,&
       ai_npatch,&
       api_nel,  &
       LM_Columns,&
       nen_C,&
       LM_Rows,&
       nen_R,&
       api_columns,&
       api_occ)
    ! _C FOR ROWS
    ! _R FOR COLUMNS
    implicit none
    sll_int32 :: ai_nC
    sll_int32 :: ai_nR
    sll_int32 :: ai_npatch
    sll_int32, dimension(:) :: api_nel
    sll_int32, dimension(:,:,:), intent(in) :: LM_Columns, LM_Rows 
    sll_int32, dimension(:),intent(in) :: nen_C,nen_R
    sll_int32, dimension(:,:), pointer :: api_columns
    sll_int32, dimension(:), pointer :: api_occ
    !local var
    sll_int32 :: li_e
    sll_int32 :: li_nel
    sll_int32 :: li_id
    sll_int32 :: li_nen_C
    sll_int32 :: li_nen_R
    sll_int32 :: li_b_C
    sll_int32 :: li_A_C
    sll_int32 :: li_b_R
    sll_int32 :: li_A_R
    sll_int32 :: li_i
    sll_int32 :: li_result
    sll_int32 :: ierr
    logical :: ll_done
    sll_int32, dimension(2) :: lpi_size
    sll_int32, dimension(:,:), pointer :: lpi_columns
    
    do li_id = 1, ai_npatch
       
       li_nel = api_nel (li_id)
       
       ! WE FIRST COMPUTE, FOR EACH ROW, THE NUMBER OF COLUMNS THAT WILL BE USED
       do li_e = 1, li_nel
          !                print *,"li_id, li_e=",li_id, li_e
          li_nen_C =  nen_C(li_id) ! ao_con_C % opi_nen(li_id, li_e)  ????
          do li_b_C = 1, li_nen_C
             
             li_A_C = LM_Columns(li_id, li_b_C, li_e)  
             if (li_A_C == 0) then
                cycle
             end if
             
             li_nen_R = nen_R(li_id)!ao_con_R % opi_nen(li_id, li_e)
             do li_b_R = 1, li_nen_R
                
                li_A_R = LM_Rows(li_id, li_b_R, li_e)
                if (li_A_R == 0) then
                   cycle
                end if
                
                ll_done = .false.
                ! WE CHECK IF IT IS THE FIRST OCCURANCE OF THE COUPLE (li_A_C, li_A_R)
                do li_i = 1, api_columns(li_A_C, 0)
                   
                   if (api_columns(li_A_C, li_i) /= li_A_R) then
                      cycle
                   end if
                   
                   ll_done = .true.
                   exit
                   
                end do
                
                if (.not.ll_done) then
                   
                   api_occ(li_A_C) = api_occ(li_A_C) + 1
                   
                   ! li_A_C IS THE ROW NUM, li_A_R THE COLUMN NUM
                   ! INITIALIZATION OF THE SPARSE MATRIX
                   api_columns(li_A_C, 0) = api_columns(li_A_C, 0) + 1
                   api_columns(li_A_C, api_columns(li_A_C, 0)) = li_A_R
                   
                   ! resizing the array
                   lpi_size(1) = SIZE(api_columns, 1)
                   lpi_size(2) = SIZE(api_columns, 2)
                   if (lpi_size(2) < api_columns(li_A_C, 0)) then
                      SLL_ALLOCATE(lpi_columns(lpi_size(1), lpi_size(2)),ierr)
                      lpi_columns = api_columns
                      
                      SLL_DEALLOCATE(api_columns,ierr)
                      
                      SLL_ALLOCATE(api_columns(lpi_size(1), 2 * lpi_size(2)),ierr)
                      api_columns(1:lpi_size(1), 1:lpi_size(2)) = lpi_columns(1:lpi_size(1), 1:lpi_size(2))
                      
                      SLL_DEALLOCATE(lpi_columns,ierr)
                   end if
                   
                   
                end if
                
             end do
             
          end do
          
       end do
       
    END DO
    
    ! COUNT NON ZERO ELEMENTS
    li_result = SUM(api_occ(1: ai_nR))
    
    sll_count_non_zero_elts_mp = li_result
  end function sll_count_non_zero_elts_mp
  
  
  subroutine sll_init_SparseMatrix_mp(&
       self,&
       ai_nR,&
       ai_nC,&
       ai_npatch,&
       api_nel,  &
       LM_Columns,&
       nen_C,&
       LM_Rows,&
       nen_R,&
       api_columns,&
       api_occ)
    ! _C FOR ROWS
    ! _R FOR COLUMNS
    implicit none
    type(sll_csr_matrix) :: self
    sll_int32 :: ai_nC
    sll_int32 :: ai_nR
    sll_int32 :: ai_npatch
    sll_int32, dimension(:) :: api_nel
    sll_int32, dimension(:),intent(in) :: nen_C,nen_R
    sll_int32, dimension(:,:,:), intent(in) :: LM_Columns, LM_Rows 
    sll_int32, dimension(:,:), pointer :: api_columns
    sll_int32, dimension(:), pointer :: api_occ
    !local var
    sll_int32 :: li_nel
    sll_int32 :: li_id
    sll_int32 :: li_e
    sll_int32 :: li_b_C
    sll_int32 :: li_A_C
    !sll_int32 :: li_b_R
    !sll_int32 :: li_A_R
    !sll_int32 :: li_index
    sll_int32 :: li_i
    sll_int32 :: li_size
    sll_int32 :: li_nen_C
    sll_int32 :: li_err
    sll_int32 :: li_flag
    sll_real64, dimension(:), pointer :: lpr_tmp
    
    ! INITIALIZING ia
    self % opi_ia(1) = 1
    
    do li_i = 1, self%num_rows !self % oi_nR
       
       self % opi_ia(li_i + 1) = self % opi_ia(1) + SUM(api_occ(1: li_i))
       
    end do
    
    ! INITIALIZING ja
    DO li_id = 1, ai_npatch
       
       li_nel = api_nel (li_id)
       
       ! WE FIRST COMPUTE, FOR EACH ROW, THE NUMBER OF COLUMNS THAT WILL BE USED
       do li_e = 1, li_nel
          
          li_nen_C = nen_C(li_id)!ao_con_C % opi_nen(li_id, li_e)
          do li_b_C = 1, li_nen_C
             
             li_A_C = LM_Columns(li_id, li_b_C, li_e)
             !ao_con_C % opi_LM(li_id, li_b_C, li_e)
             
             if (li_A_C == 0) then
                cycle
             end if
             
             if (api_columns(li_A_C, 0) == 0) then
                cycle
             end if
             
             li_size = api_columns(li_A_C, 0)
             
             allocate ( lpr_tmp(li_size), stat = li_err)
             if (li_err .ne. 0) li_flag = 10
             
             lpr_tmp(1: li_size) = real( api_columns(li_A_C, 1: li_size))
             
             call QsortC(lpr_tmp)
             
             do li_i = 1, li_size
                
                self % opi_ja(self % opi_ia(li_A_C) + li_i - 1) = int ( lpr_tmp(li_i))
                
             end do
             
             api_columns(li_A_C, 0) = 0
             deallocate ( lpr_tmp)
             
          end do
          
       end do
       
    end do
    
  end subroutine sll_init_SparseMatrix_mp

      
      
      
end module sll_sparse_matrix_mp_module
