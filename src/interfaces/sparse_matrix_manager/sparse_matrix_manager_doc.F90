! This file is read by doxygen software
! Change it to match with your library
! http://www.stack.nl/~dimitri/doxygen/manual/docblocks.html#fortranblocks
! To see the documentation, remove the line containing your directory 
! in file Doxyfile.in (line 691) if it is excluded. 
! Type 'make doc' in build directory.
! To check the results, open : 
! selalib/doc/build/html/doxygen/html/defgroup sparse_matrix_managers.html 
! The following lines will be read by doxygen to generate documentation:


!> @defgroup sparse_matrix_manager sll_sparse_matrix_manager 
!> @brief 
!>  avoid to have a hard dependecy on any Linear library
!> @author Ahmed Ratnani
!> @details
!> See http://www.ratnani.org/jorek_doc/spm.html
!>
!> <b> Headers file available </b>
!>
!> <b> Modules available </b>
!>
!> <b> How to use it </b>
!>
!> <b> Examples </b>
!>
!> Proposition for sparse matrix interface (Pierre)
!>
!> Initialization
!>
!> sll_s_sparse_matrix_init(A, spfmt, num_rows, num_cols, row, col, val)
!> 
!> - A derived type with arrays that describes the matrix
!> - spfmt indicates which sparse format is to be used 
!>   - SLL_P_SPARSE_MATRIX_COO
!>   - SLL_P_SPARSE_MATRIX_CSR
!>   - SLL_P_SPARSE_MATRIX_CSC
!> 
!> num_rows indicates the number of rows in the sparse matrix.
!> num_cols indicates the number of columns in the sparse matrix.
!> 
!> row is an integer parallel array of rank 1. Its length and content can vary, 
!> depending on which sparse format is used:
!> 
!> SLL_P_SPARSE_MATRIX_COO - row is of the same size as arrays col and val and 
!>   contains row indices of the nonzero elements in array val.
!> SLL_P_SPARSE_MATRIX_CSR - row is of size num_rows+1 and contains pointers to the 
!>   beginning of each row in arrays col and val.
!> SLL_P_SPARSE_MATRIX_CSC - row is of size num_cols+1 and contains pointers to the 
!>   beginning of each column in arrays col and val.
!> col is an integer global array of rank 1 with the same length as array val. 
!>   Its use will vary, depending on which sparse format is used:
!> 
!> SLL_P_SPARSE_MATRIX_COO - col contains column indices of the corresponding 
!>   elements stored in array val.
!> SLL_P_SPARSE_MATRIX_CSR - col contains column indices of the corresponding 
!>  elements stored in array val.
!> SLL_P_SPARSE_MATRIX_CSC - col contains row indices of the corresponding 
!>  elements in array val.
!> val is a parallel array of rank 1, containing the nonzero elements of the 
!>  sparse matrix. The storage pattern varies, depending on which sparse 
!>  format is used:
!> 
!> SLL_P_SPARSE_MATRIX_COO - Nonzero elements can be stored in any order.
!> SLL_P_SPARSE_MATRIX_CSR - Nonzero elements should be stored row by row, from 
!>  row 1 to row num_rows.
!> SLL_P_SPARSE_MATRIX_CSC - Nonzero elements should be stored column by column, 
!>  from column 1 to column num_cols.
!> The length of val is nnz for all three formats. This parameter represents 
!>  the total number of nonzero elements in the sparse matrix. 
!> 
